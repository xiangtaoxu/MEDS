!==========================================================================================!
! meds_demography_cohort_fusefiss -- COHORT numerical-resolution control for the adaptive cohort      !
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
! meds_demography_patch_fusefiss, which depends on this one (patch fusion re-sorts cohorts).          !
!==========================================================================================!
module meds_demography_cohort_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num
   use meds_allometry,  only : height_to_dbh, min_cohort_carbon
   use meds_column_params, only : LEAF_TEMP_INIT
   use meds_config,     only : meds_config_t
   use meds_site_diag_types,  only : cohort_diag_fuse, CDIAG_FUSE, CSDIAG_FUSE
   use meds_site_state_types, only : site_t, cohort_reorder, rebuild_csr, cohort_compact,        &
                                      cohort_ensure_capacity, copy_cohort_slot, init_cohort,       &
                                      scale_cohort_ground_fields, cohort_tissue_heat_capacity,      &
                                      cohort_tissue_water, TISSUE_C_LEAF, TISSUE_C_SAPW,            &
                                      TISSUE_HCAP_MIN,                                              &
                                      set_cohort_size_from_carbon, assign_cohort_id,       &
                                      fuse_cohort_fast_state
   use meds_litter_partition, only : necromass_to_litter
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
         !----- The fast-loop-owned state fuses by its DECLARED per-field policy, in one place  !
         !      (meds_site_state_types%fuse_cohort_fast_state): leaf-area-weighted for the tissue   !
         !      temperatures, nplant-weighted for the per-plant amounts, summed for the ground-      !
         !      referenced film. Enumerating it here is how the film-water convention got it wrong    !
         !      once already (PR #119). The dynamic leaf traits stay here: they are slow-loop state,   !
         !      not fast, and they must be set BEFORE the survivor's geometry is re-derived because    !
         !      sla maps its leaf_carbon <-> leaf_area. --------------------------------------------!
         wr = nr * cohort%leaf_area(recc) ; wd = nd * cohort%leaf_area(donc) ; wtot = wr + wd
         call fuse_cohort_fast_state(cohort, recc, donc, nr, nd)
         if (wtot > tiny_num) then
            cohort%sla(recc)     = (wr * cohort%sla(recc)     + wd * cohort%sla(donc))     / wtot
            cohort%vcmax25(recc) = (wr * cohort%vcmax25(recc) + wd * cohort%vcmax25(donc)) / wtot
            cohort%rd25(recc)    = (wr * cohort%rd25(recc)    + wd * cohort%rd25(donc))    / wtot
            cohort%llspan(recc)  = (wr * cohort%llspan(recc)  + wd * cohort%llspan(donc))  / wtot
         end if
         !----- Fast-loop DIAGNOSTIC accumulators. Handed the SAME two weight pairs used above --   !
         !      leaf area for the intensive quantities, nplant for the extensive ones -- so a        !
         !      diagnostic and its prognostic twin can never be fused on different weights. Which     !
         !      of the two applies is declared per field in meds_site_diag_types (CDIAG_FUSE), not     !
         !      decided here.  -----------------------------------------------------------------!
         call cohort_diag_fuse(cohort%diag,  recc, donc, wr, wd, nr, nd, CDIAG_FUSE)
         call cohort_diag_fuse(cohort%sdiag, recc, donc, wr, wd, nr, nd, CSDIAG_FUSE)
         !----- The survivor keeps its own moving-average growth history (ring buffer + accum  !
         !      + count + growth_avg are left untouched); the donor's is discarded with it. ---!
         !----- Same convention for the predawn water status. Fusion requires the two cohorts to  !
         !      already match on height and LAI, so they sit in the same light and rooting        !
         !      environment and their predawn potentials agree closely; survivor-keeps is a       !
         !      smaller approximation here than the height/LAI similarity test already tolerates. !
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
      real(wp)    :: eps, agb_before, agb_after, wc0, film_before, film_after

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
         film_before = sum(site%cohort%leaf_surf_water(1:n0) + site%cohort%wood_surf_water(1:n0))
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
               !----- the interception films are per m2 GROUND, not per plant: each daughter keeps  !
               !      half, or the split creates film water (2026-09 review, item 1B #1). ----------!
               call scale_cohort_ground_fields(cohort, i, 0.5_wp)
               call scale_cohort_ground_fields(cohort, m, 0.5_wp)
               cohort%wood_carbon(i) = wc0 * (1.0_wp + eps)
               cohort%wood_carbon(m) = wc0 * (1.0_wp - eps)
               call set_cohort_size_from_carbon(cohort, i)
               call set_cohort_size_from_carbon(cohort, m)
            end do
            cohort%n = m
            agb_after  = sum(cohort%nplant(1:m) * cohort%agb(1:m))
            film_after = sum(cohort%leaf_surf_water(1:m) + cohort%wood_surf_water(1:m))
         end associate
         !----- The '-eps' daughters (slots n0+1..m) are NEW cohorts -> fresh global ids; the !
         !      '+eps' half kept slot i and its parent id (the continuation).                 !
         do i = n0 + 1_ik, m
            call assign_cohort_id(site, i)
         end do
         if (abs(agb_after - agb_before) > cfg%conservation_tol * max(agb_before, tiny_num))       &
            error stop 'split_cohorts: AGB conservation violated'
         if (abs(film_after - film_before) > 1.0e-12_wp * max(film_before, tiny_num))              &
            error stop 'split_cohorts: canopy film water not conserved'
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
   subroutine terminate_cohorts(site, cfg, water_shed)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      !----- What the cull HANDED ON, for the caller to declare to the site ledger (plan §10.2).   !
      !      Reported rather than declared here because meds_demography must not link the ledger:  !
      !      the engine applies, the driver accounts. --------------------------------------------!
      real(wp), optional, intent(out) :: water_shed  !< [kg/m2 site] tissue + film water -> the ground
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n, pf, ip
      real(wp)             :: lab_g, lab_s, str_g, str_s, lig_g, lig_s
      real(wp)             :: w_tot, w_cohort

      w_tot = 0.0_wp
      if (present(water_shed)) water_shed = 0.0_wp
      n = site%cohort%n
      if (n < 1_ik) return
      allocate(keep(n))
      associate (cohort => site%cohort, pft => cfg%pft, patch => site%patch)
         do i = 1_ik, n
            keep(i) = (cohort%nplant(i) * cohort%agb(i) >= cfg%min_cohort_agb) .and.               &
                      (cohort%nplant(i) >= cfg%negligible_nplant)
            if (keep(i)) cycle
            ip = cohort%owner_patch(i)
            !----- A culled cohort's WATER goes to the ground, down the same channel turnover      !
            !      shedding already uses -- one verified path, not a second mechanism. It was      !
            !      simply dropped by cohort_compact before (review item 1B #7). Its tissue HEAT is !
            !      a change in thermal MASS and is declared with the rest of the phase's, by the   !
            !      caller, since fusion and fission change it too and they are one mechanism. -----!
            w_cohort = cohort_tissue_water(cohort, i)                                              &
                     + cohort%leaf_surf_water(i) + cohort%wood_surf_water(i)
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + w_cohort / max(cfg%dt_slow, tiny_num)
            w_tot = w_tot + patch%area(ip) * w_cohort
            if (.not. cfg%soil_carbon_on) cycle
            pf = cohort%pft(i)
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
      if (present(water_shed)) water_shed = w_tot
      if (all(keep)) return
      call cohort_compact(site%cohort, keep)
      call rebuild_csr(site)
   end subroutine terminate_cohorts

   !---------------------------------------------------------------------------------------!
   ! Spawn from the carry-forward recruit pool: when a pool reaches `min_recruit_size`, ONE new  !
   ! cohort is born at the shared minimum cohort height (the smallest tracked size) and the pool  !
   ! is reset; otherwise it carries over, so rare recruiters still establish eventually. A HOST   !
   ! structural process -- it changes the cohort count -- so it lives with the cohort fuse/fission !
   ! housekeeping.                                                                                 !
   !                                                                                          !
   ! The pool is CREDITED DAILY by the driver (accumulate_recruit_pool) and CONSUMED monthly here. !
   ! It used to be credited here too, from whatever recruitment rate the driver had computed on    !
   ! this one day, scaled up to stand for the whole month -- a 12-point sample of a quantity the   !
   ! model computes 365 times a year, and one that made the pool's carbon content depend on which  !
   ! days happened to be month boundaries.  ------------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(site, cfg, carbon_drawn, heat_drawn)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      !----- What a recruit brought IN FROM OUTSIDE the modelled system, for the caller to declare  !
      !      (plan §10.2). Recruitment is a stand-in for everything that happens between a seed and !
      !      a 2 m sapling -- germination, and the seedling's own photosynthesis, transpiration and !
      !      energy balance -- none of which this model represents, because it tracks no cohort     !
      !      below `min_cohort_height`. The carbon and heat a recruit arrives with are therefore    !
      !      genuinely EXTERNAL: they were fixed and absorbed by a size class outside the model.    !
      !      Drawing them from the free atmosphere rather than from this patch's canopy air is      !
      !      deliberate -- crediting the CAS with seedling uptake while representing none of the    !
      !      seedling's respiration, transpiration or shading would add one term of a missing       !
      !      process and call it an improvement.                                                    !
      !                                                                                          !
      !      Only the part the model did NOT already pay for is external: `recruit_pool` carries    !
      !      reproduction carbon valued at `carbon_min` per plant, and that much IS debited from    !
      !      the parents (or, for seed rain, declared as it arrives). The remainder of the          !
      !      endowment is the draw.  --------------------------------------------------------!
      real(wp), optional,  intent(out)   :: carbon_drawn  !< [kgC/m2 site]
      real(wp), optional,  intent(out)   :: heat_drawn    !< [J/m2 site]
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh, c_min, endow, t_birth, cap_leaf, cap_wood, c_tot, e_tot

      c_tot = 0.0_wp ; e_tot = 0.0_wp
      if (present(carbon_drawn)) carbon_drawn = 0.0_wp
      if (present(heat_drawn))   heat_drawn   = 0.0_wp
      np = site%patch%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the smallest tracked size -> the same diameter. ----------!
      recruit_dbh = height_to_dbh(cfg%pft%min_cohort_height)

      !                                                                                     !
      !      The pool is CREDITED DAILY by the vegetation driver (accumulate_recruit_pool), not  !
      !      here. This routine only spawns from what has accumulated.  ------------------------!

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            if (site%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      !----- No pool reached the spawn threshold. The seed rain still ARRIVED this month, so the  !
      !      import must be reported before returning -- an early exit that skips the report        !
      !      silently under-declares on exactly the months nothing is born (two thirds of them in   !
      !      the run this was found on).  ---------------------------------------------------------!
      if (nspawn == 0_ik) then
         if (present(carbon_drawn)) carbon_drawn = c_tot
         return
      end if

      call cohort_ensure_capacity(site%cohort, site%cohort%n + nspawn)
      n_before = site%cohort%n
      m = site%cohort%n
      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, site%n_pft
               if (patch%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               !----- Born at its PATCH's canopy-air temperature, not at the global LEAF_TEMP_INIT !
               !      constant. Guarded: a canopy air that has never been stepped (a bare setup    !
               !      call, a restart before the first fast window) reads as unset, and the        !
               !      constant is the right answer there. ---------------------------------------!
               t_birth = patch%cas(ip)%can_temp
               if (t_birth < 100.0_wp) t_birth = LEAF_TEMP_INIT
               call init_cohort(cohort, m, pft, pf, ip, patch%recruit_pool(pf, ip), recruit_dbh,   &
                                birth_temp=t_birth)
               !----- The external draw: everything the new cohort holds, less the reproduction    !
               !      carbon the pool actually paid for at `carbon_min` per plant. ---------------!
               c_min = min_cohort_carbon(pft%min_cohort_height, pft%wood_density(pf))
               endow = cohort%leaf_carbon(m) + cohort%fineroot_carbon(m)                           &
                     + cohort%wood_carbon(m) + cohort%nonstructural_carbon(m)
               c_tot = c_tot + patch%area(ip) * patch%recruit_pool(pf, ip) * (endow - c_min)
               call cohort_tissue_heat_capacity(cohort, m, TISSUE_C_LEAF, TISSUE_C_SAPW,           &
                                                TISSUE_HCAP_MIN, cap_leaf, cap_wood)
               e_tot = e_tot + patch%area(ip) * (cap_leaf + cap_wood) * t_birth
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      if (present(carbon_drawn)) carbon_drawn = c_tot
      if (present(heat_drawn))   heat_drawn   = e_tot

      !----- Stamp each freshly spawned cohort with a persistent global id. ----------------!
      do k = n_before + 1_ik, site%cohort%n
         call assign_cohort_id(site, k)
      end do

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine apply_recruitment

end module meds_demography_cohort_fusefiss
