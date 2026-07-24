!==========================================================================================!
! meds_core_cohort_fusefiss -- COHORT numerical-resolution control for the adaptive cohort      !
! discretization (the cohort half of the ED2 fuse_fiss_utils analogue) plus cohort RECRUITMENT.  !
! Driven by ARITHMETIC, not ecology: keep the per-patch cohort count bounded and the             !
! representation conservative WITHOUT changing what the ecosystem does.                          !
!                                                                                          !
!   Sorting        -- cohorts tallest-first within each patch (height desc, DBH desc tie). The  !
!                     order is load-bearing for the overtopping sweep and the fusion scans.     !
!   Cohort fusion  -- receptor (taller) absorbs same-PFT donors of similar HEIGHT unless the    !
!   /fission          merged cohort LAI exceeds the cap; split a cohort whose LAI exceeds it.   !
!                     The conserved invariant is TOTAL ABOVEGROUND BIOMASS (carbon):            !
!                     fused/split diameters are re-derived from conserved AGB, never averaged.  !
!   Termination    -- cull cohorts below the density/AGB floors.                                !
!   Recruitment    -- apply_recruitment: accumulate the supplied per-(PFT,patch) recruit        !
!                     density and spawn min-size cohorts (a HOST structural process).           !
! The PATCH half (sort_patches, patch fusion, disturbance) is the sibling module                !
! meds_core_patch_fusefiss, which depends on this one (patch fusion re-sorts cohorts).          !
!==========================================================================================!
module meds_core_cohort_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num, mon_per_yr
   use meds_allometry,  only : height_to_dbh
   use meds_config,     only : meds_config_t
   use meds_core_state_types, only : site_t, cohort_reorder, rebuild_csr, cohort_compact,        &
                                      cohort_ensure_capacity, copy_cohort_slot, init_cohort,       &
                                      set_cohort_size_from_carbon, assign_cohort_id
   use meds_column_state_types, only : necromass_to_litter
   implicit none
   private

   public :: sort_cohorts, apply_recruitment
   public :: new_fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts, max_cohort_count

contains

   !---------------------------------------------------------------------------------------!
   ! Sort cohorts within every patch slice: height descending, DBH descending on ties.      !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: perm(:)
      integer(ik) :: ip, i0, i1, i, j, key
      integer(ik) :: n

      n = site%cohort%n
      if (n < 1_ik) return
      allocate(perm(n))
      do i = 1_ik, n
         perm(i) = i
      end do
      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            !----- Insertion sort of the index window perm(i0:i1), descending. ------------!
            do i = i0 + 1_ik, i1
               key = perm(i)
               j   = i - 1_ik
               do while (j >= i0)
                  if (higher(cohort%height(perm(j)), cohort%dbh(perm(j)), cohort%height(key), cohort%dbh(key))) exit
                  perm(j + 1_ik) = perm(j)
                  j = j - 1_ik
               end do
               perm(j + 1_ik) = key
            end do
         end do
      end associate
      call cohort_reorder(site%cohort, perm, n)
   end subroutine sort_cohorts

   !----- True if (h1,d1) ranks at least as high as (h2,d2) under the sort order. ----------!
   pure logical function higher(h1, d1, h2, d2)
      real(wp), intent(in) :: h1, d1, h2, d2
      if (h1 > h2) then
         higher = .true.
      else if (h1 < h2) then
         higher = .false.
      else
         higher = (d1 >= d2)
      end if
   end function higher

   integer(ik) function max_cohort_count(site)
      type(site_t), intent(in) :: site
      if (site%patch%n < 1_ik) then
         max_cohort_count = 0_ik
      else
         max_cohort_count = maxval(site%patch%cohort_count(1:site%patch%n))
      end if
   end function max_cohort_count

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion with geometric tolerance relaxation.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus, maxc
      logical     :: force

      if (cfg%max_cohort == 0_ik) return
      force = (cfg%max_cohort < 0_ik)
      maxc  = abs(cfg%max_cohort)
      tol   = cfg%cohort_size_tol_min

      do ifus = 1_ik, cfg%n_cohort_fusion_iter
         call fuse_pass(site, cfg, tol, force)
         if (.not. force) then
            if (max_cohort_count(site) <= maxc) exit
         end if
         tol = tol * cfg%cohort_size_tol_mult
      end do
   end subroutine new_fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(site, cfg, tol, force)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_height, lai_comb
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(site%cohort%n))
      alive = .true.
      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (cohort%pft(donc) /= cohort%pft(recc)) cycle
                  lai_comb = cohort%nplant(recc) * cohort%leaf_area(recc) + cohort%nplant(donc) * cohort%leaf_area(donc)
                  if (.not. force .and. lai_comb >= cfg%cohort_lai_cap) cycle
                  diff_height = abs(cohort%height(donc) - cohort%height(recc))
                  !----- Same-PFT fusion (checked above), so recc and donc share hgt_max. -------!
                  if (force .or. diff_height < cohort%p_hgt_max(recc) * tol) then
                     call fuse_2_cohorts(site, recc, donc, cfg%conservation_tol)
                     alive(donc) = .false.
                     did_fuse    = .true.
                  end if
               end do
            end do
         end do
      end associate
      if (did_fuse) then
         call cohort_compact(site%cohort, alive)
         call rebuild_csr(site)
         call sort_cohorts(site)
      end if
   end subroutine fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor cohort into receptor, conserving plant number and total aboveground         !
   ! biomass (carbon). The fused per-plant AGB is the nplant-weighted mean; the diameter is   !
   ! re-derived from it (agb_to_dbh) and the rest of the geometry refreshed.                  !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(site, recc, donc, conservation_tol)
      type(site_t), intent(inout) :: site
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, agb_tot, agb_new, wr, wd, wtot
      real(wp) :: lc_tot, fc_tot, wc_tot, nc_tot                !< pre-merge pool totals [kgC/m2]

      associate (cohort => site%cohort)
         nr      = cohort%nplant(recc)
         nd      = cohort%nplant(donc)
         ntot    = nr + nd
         agb_tot = nr * cohort%agb(recc) + nd * cohort%agb(donc)     ! [kgC/m2] conserved
         !----- Leaf-area-weighted merge of the fast per-cohort INTENSIVE state (leaf_temp is a       !
         !      temperature; weight by each cohort's total leaf area BEFORE set_cohort_size below      !
         !      re-derives the survivor's leaf_area). Without this the donor's heat state is           !
         !      silently dropped while the AGB assert still passes. Internal water MASS is NOT here    !
         !      -- it is EXTENSIVE (a per-plant quantity, like the carbon pools below), so it belongs   !
         !      in the nplant-weighted block, not this leaf-area-weighted one (MEDS_ED2_RK45_DESIGN.md !
         !      sec 9 "Cohort fusion semantics change": leaf-area-weighting an extensive quantity        !
         !      would silently violate column water conservation on every fuse). ---------------------!
         wr = nr * cohort%leaf_area(recc) ; wd = nd * cohort%leaf_area(donc) ; wtot = wr + wd
         if (wtot > tiny_num) then
            cohort%leaf_temp(recc) = (wr * cohort%leaf_temp(recc) + wd * cohort%leaf_temp(donc)) / wtot
            cohort%wood_temp(recc) = (wr * cohort%wood_temp(recc) + wd * cohort%wood_temp(donc)) / wtot
            !----- Dynamic leaf traits are intensive (per leaf area): leaf-area-weight them too, so a  !
            !      fusion of a sun + shade cohort keeps the area-mean trait (set BEFORE the survivor's  !
            !      geometry is re-derived, since sla maps its leaf_carbon <-> leaf_area). -------------!
            cohort%sla(recc)     = (wr * cohort%sla(recc)     + wd * cohort%sla(donc))     / wtot
            cohort%vcmax25(recc) = (wr * cohort%vcmax25(recc) + wd * cohort%vcmax25(donc)) / wtot
            cohort%rd25(recc)    = (wr * cohort%rd25(recc)    + wd * cohort%rd25(donc))    / wtot
            cohort%llspan(recc)  = (wr * cohort%llspan(recc)  + wd * cohort%llspan(donc))  / wtot
         end if
         !----- Accumulated carbon fluxes are per-plant [kgC/plant] (extensive per ground); nplant- !
         !      weight so the site totals are conserved. gpp_accum + the three maintenance-resp       !
         !      accumulators are enumerated here BY HAND (fuse_2_cohorts is NOT the centralized        !
         !      reorder) -- add any new per-plant accumulator to this list.                            !
         cohort%gpp_accum(recc)       = (nr * cohort%gpp_accum(recc)       + nd * cohort%gpp_accum(donc))       / ntot
         cohort%leaf_resp_accum(recc) = (nr * cohort%leaf_resp_accum(recc) + nd * cohort%leaf_resp_accum(donc)) / ntot
         cohort%stem_resp_accum(recc) = (nr * cohort%stem_resp_accum(recc) + nd * cohort%stem_resp_accum(donc)) / ntot
         cohort%root_resp_accum(recc) = (nr * cohort%root_resp_accum(recc) + nd * cohort%root_resp_accum(donc)) / ntot
         !----- Internal water mass [kg/plant] is EXTENSIVE (like AGB): nplant-weight so the fused    !
         !      site-total TOTAL water (nplant*mass, summed over cohorts) is conserved across the      !
         !      fuse -- leaf-area-weighting it (the OLD psi treatment) would not conserve total water.  !
         cohort%leaf_water_mass(recc) = (nr * cohort%leaf_water_mass(recc) + nd * cohort%leaf_water_mass(donc)) / ntot
         cohort%wood_water_mass(recc) = (nr * cohort%wood_water_mass(recc) + nd * cohort%wood_water_mass(donc)) / ntot
         !----- Surface (interception film) water [kg/m2 GROUND] is the OPPOSITE convention from the  !
         !      internal water mass just above: it is ALREADY ground-area-referenced (not per-plant),  !
         !      so two cohorts' contributions to the SAME patch ground area simply ADD -- no nplant     !
         !      weighting (that would double-count the area normalization already baked into each term). !
         cohort%leaf_surf_water(recc) = cohort%leaf_surf_water(recc) + cohort%leaf_surf_water(donc)
         cohort%wood_surf_water(recc) = cohort%wood_surf_water(recc) + cohort%wood_surf_water(donc)
         !----- The survivor keeps its own moving-average growth history (ring buffer + accum  !
         !      + count + growth_avg are left untouched); the donor's is discarded with it. ---!
         cohort%nplant(recc) = ntot
         !----- Conserve the four PROGNOSTIC carbon pools (nplant-weighted) FIRST, THEN derive     !
         !      geometry from the conserved wood_carbon anchor. set_cohort_size_from_carbon takes   !
         !      the pools as INPUTS (does NOT overwrite them, unlike set_cohort_size).              !
         lc_tot = nr * cohort%leaf_carbon(recc)          + nd * cohort%leaf_carbon(donc)
         fc_tot = nr * cohort%fineroot_carbon(recc)      + nd * cohort%fineroot_carbon(donc)
         wc_tot = nr * cohort%wood_carbon(recc)          + nd * cohort%wood_carbon(donc)
         nc_tot = nr * cohort%nonstructural_carbon(recc) + nd * cohort%nonstructural_carbon(donc)
         cohort%leaf_carbon(recc)          = lc_tot / ntot
         cohort%fineroot_carbon(recc)      = fc_tot / ntot
         cohort%wood_carbon(recc)          = wc_tot / ntot
         cohort%nonstructural_carbon(recc) = nc_tot / ntot
         call set_cohort_size_from_carbon(site%cohort, recc)   ! dbh=wood_to_dbh(wood_carbon); pools kept
         !----- Carbon-pool conservation guard (per pool, density). ------------------------!
         if (abs(ntot*cohort%leaf_carbon(recc)          - lc_tot) > conservation_tol*max(lc_tot,tiny_num) .or. &
             abs(ntot*cohort%fineroot_carbon(recc)      - fc_tot) > conservation_tol*max(fc_tot,tiny_num) .or. &
             abs(ntot*cohort%wood_carbon(recc)          - wc_tot) > conservation_tol*max(wc_tot,tiny_num) .or. &
             abs(ntot*cohort%nonstructural_carbon(recc) - nc_tot) > conservation_tol*max(nc_tot,tiny_num))     &
            error stop 'fuse_2_cohorts: carbon-pool conservation violated'
         !----- AGB-density conservation guard (agb = p_aboveground_frac*wood_carbon, and          !
         !      nplant*wood_carbon was conserved). ------------------------------------------------!
         agb_new = cohort%nplant(recc) * cohort%agb(recc)
         if (abs(agb_new - agb_tot) > conservation_tol * max(agb_tot, tiny_num))             &
            error stop 'fuse_2_cohorts: AGB conservation violated'
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose LAI density exceeds the cap, iterating until none do.          !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: eps, agb_before, agb_after, wc0

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps

      do iter = 1_ik, cfg%n_cohort_fusion_iter
         n0 = site%cohort%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (site%cohort%nplant(i) * site%cohort%leaf_area(i) > cfg%cohort_lai_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         agb_before = sum(site%cohort%nplant(1:n0) * site%cohort%agb(1:n0))
         call cohort_ensure_capacity(site%cohort, n0 + nsplit)
         m = n0
         associate (cohort => site%cohort)
            do i = 1_ik, n0
               if (cohort%nplant(i) * cohort%leaf_area(i) <= cfg%cohort_lai_cap) cycle
               cohort%nplant(i) = 0.5_wp * cohort%nplant(i)
               m = m + 1_ik
               !----- Perturb the wood_carbon ANCHOR +/-eps (NOT dbh); the leaf/fineroot/            !
               !      nonstructural pools stay identical in both daughters. This conserves ALL four    !
               !      pools + AGB EXACTLY (0.5n*wc*(1+eps)+0.5n*wc*(1-eps)=n*wc) and is immune to the   !
               !      hgt_max cap (no dbh-renorm approximation).  ------------------------------------!
               wc0 = cohort%wood_carbon(i)
               call copy_cohort_slot(cohort, m, i)           ! copies halved nplant + params + pools
               cohort%wood_carbon(i) = wc0 * (1.0_wp + eps)
               cohort%wood_carbon(m) = wc0 * (1.0_wp - eps)
               call set_cohort_size_from_carbon(cohort, i)
               call set_cohort_size_from_carbon(cohort, m)
            end do
            cohort%n = m
            agb_after = sum(cohort%nplant(1:m) * cohort%agb(1:m))
         end associate
         !----- The '-eps' daughters (slots n0+1..m) are NEW cohorts -> fresh global ids; the !
         !      '+eps' half kept slot i and its parent id (the continuation).                 !
         do i = n0 + 1_ik, m
            call assign_cohort_id(site, i)
         end do
         if (abs(agb_after - agb_before) > cfg%conservation_tol * max(agb_before, tiny_num))       &
            error stop 'split_cohorts: AGB conservation violated'
         call rebuild_csr(site)
         call sort_cohorts(site)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! Cull cohorts below the AGB-density floor or the absolute density floor. The ENTIRE        !
   ! remaining carbon of a culled cohort (leaf/fine-root/wood/storage; it is being removed      !
   ! outright, unlike a turnover shed) becomes litter into its own patch's soil-carbon pools    !
   ! (B1, MEDS_SLOW_DYNAMICS_DESIGN.md Part II; OPT-IN [soil_carbon].soil_carbon_on -- default    !
   ! .false. keeps this bit-identical) -- added directly onto the named fields since this module  !
   ! cannot link biogeochemistry (necromass_to_litter is DAG-safe: plain scalars).                !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n, pf, ip
      real(wp)             :: lab_g, lab_s, str_g, str_s, lig_g, lig_s

      n = site%cohort%n
      if (n < 1_ik) return
      allocate(keep(n))
      associate (cohort => site%cohort, pft => cfg%pft, patch => site%patch)
         do i = 1_ik, n
            keep(i) = (cohort%nplant(i) * cohort%agb(i) >= cfg%min_cohort_agb) .and.               &
                      (cohort%nplant(i) >= cfg%negligible_nplant)
            if (keep(i) .or. .not. cfg%soil_carbon_on) cycle
            pf = cohort%pft(i)
            ip = cohort%owner_patch(i)
            call necromass_to_litter(cohort%nplant(i) * cohort%leaf_carbon(i),                    &
                     cohort%nplant(i) * cohort%fineroot_carbon(i),                                 &
                     cohort%nplant(i) * cohort%wood_carbon(i),                                     &
                     cohort%nplant(i) * cohort%nonstructural_carbon(i),                             &
                     pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                                  &
                     pft%aboveground_frac(pf), pft%struct_lignin_frac(pf),                          &
                     lab_g, lab_s, str_g, str_s, lig_g, lig_s)
            patch%soil_carbon(ip)%fast_grnd_carbon   = patch%soil_carbon(ip)%fast_grnd_carbon   + lab_g
            patch%soil_carbon(ip)%fast_soil_carbon   = patch%soil_carbon(ip)%fast_soil_carbon   + lab_s
            patch%soil_carbon(ip)%struct_grnd_carbon = patch%soil_carbon(ip)%struct_grnd_carbon + str_g
            patch%soil_carbon(ip)%struct_soil_carbon = patch%soil_carbon(ip)%struct_soil_carbon + str_s
            patch%soil_carbon(ip)%struct_grnd_lignin  = patch%soil_carbon(ip)%struct_grnd_lignin  + lig_g
            patch%soil_carbon(ip)%struct_soil_lignin  = patch%soil_carbon(ip)%struct_soil_lignin  + lig_s
         end do
      end associate
      if (all(keep)) return
      call cohort_compact(site%cohort, keep)
      call rebuild_csr(site)
   end subroutine terminate_cohorts

   !---------------------------------------------------------------------------------------!
   ! Apply the supplied per-(PFT,patch) recruitment rate: accumulate it into a carry-forward !
   ! pool and, when a pool reaches `min_recruit_size`, spawn ONE new cohort at the shared      !
   ! minimum cohort height (the smallest tracked size; the pool is reset, otherwise it carries  !
   ! over so rare recruiters still establish eventually). A HOST structural process -- it changes !
   ! the cohort count -- so it lives with the cohort fuse/fission housekeeping; the recruitment    !
   ! RATE comes from the vegetation-dynamics driver (meds_vegetation_dynamics).                    !
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(site, cfg, recruitment)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh

      np = site%patch%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the smallest tracked size -> the same diameter. ----------!
      recruit_dbh = height_to_dbh(cfg%pft%min_cohort_height)

      !----- Accumulate one month of the supplied per-YEAR recruit density into the carry- !
      !      forward pool (this routine runs once a month, so add rate / mon_per_yr). -------!
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            site%patch%recruit_pool(pf, ip) = site%patch%recruit_pool(pf, ip)                &
                                              + recruitment(pf, ip) / mon_per_yr
         end do
      end do

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            if (site%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      if (nspawn == 0_ik) return

      call cohort_ensure_capacity(site%cohort, site%cohort%n + nspawn)
      n_before = site%cohort%n
      m = site%cohort%n
      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, site%n_pft
               if (patch%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               call init_cohort(cohort, m, pft, pf, ip, patch%recruit_pool(pf, ip), recruit_dbh)
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      !----- Stamp each freshly spawned cohort with a persistent global id. ----------------!
      do k = n_before + 1_ik, site%cohort%n
         call assign_cohort_id(site, k)
      end do

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine apply_recruitment

end module meds_core_cohort_fusefiss
