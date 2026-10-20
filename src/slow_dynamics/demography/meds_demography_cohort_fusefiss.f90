! SPDX-License-Identifier: Apache-2.0
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
!                     The four carbon pools are conserved, and with them AGB and leaf area:     !
!                     fused/split diameters are re-derived from the conserved wood carbon,      !
!                     never averaged.                                                           !
!   Termination    -- cull cohorts below the density/AGB floors.                                !
!   Recruitment    -- apply_recruitment: accumulate the supplied per-(PFT,patch) recruit        !
!                     density and spawn min-size cohorts (a HOST structural process).           !
! The PATCH half (sort_patches, patch fusion, disturbance) is the sibling module                !
! meds_demography_patch_fusefiss, which depends on this one (patch fusion re-sorts cohorts).          !
!==========================================================================================!
module meds_demography_cohort_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num
   use meds_numerics,   only : conserved, pair_mean
   use meds_allometry_lib, only : height_to_dbh, min_cohort_carbon
   use meds_config,     only : meds_config_t
   use meds_diag_vars,        only : fuse_cohort_diag, PA_MORT_C_CULL
   use meds_cohort_state_vars, only : reorder_cohort_block, compact_cohort_block,                        &
                                      expand_cohort_block, copy_cohort, init_cohort,           &
                                      cohort_tissue_heat_capacity, cohort_film_water, cohort_water_enthalpy, &
                                      cohort_tissue_water, set_cohort_size_from_carbon,                  &
                                      fuse_cohort_fast_state, set_cohort_tissue_temps,                   &
                                      hold_cohort_tissue_temps
   use meds_polygon_state_vars, only : polygon_t, rebuild_cohort_map, assign_cohort_id, patch_polygon_area
   use meds_patch_state_vars, only : patch_physical_state
   use meds_litter_partition, only : necromass_to_litter
   implicit none
   private

   public :: sort_cohorts, apply_recruitment
   public :: new_fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts, max_cohort_count

contains

   !---------------------------------------------------------------------------------------!
   ! Sort cohorts within every patch slice: height descending, DBH descending on ties.      !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik), allocatable :: perm(:)
      integer(ik) :: ip, i0, i1, i, j, key
      integer(ik) :: n

      n = poly%cohort%n
      if (n < 1_ik) return
      allocate(perm(n))
      do i = 1_ik, n
         perm(i) = i
      end do
      associate (cohort => poly%cohort, patch => poly%patch)
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
      call reorder_cohort_block(poly%cohort, perm, n)
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

   integer(ik) function max_cohort_count(poly)
      type(polygon_t), intent(in) :: poly
      if (poly%patch%n < 1_ik) then
         max_cohort_count = 0_ik
      else
         max_cohort_count = maxval(poly%patch%cohort_count(1:poly%patch%n))
      end if
   end function max_cohort_count

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion with geometric tolerance relaxation.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus, maxc
      logical     :: force

      if (cfg%max_cohort == 0_ik) return
      force = (cfg%max_cohort < 0_ik)
      maxc  = abs(cfg%max_cohort)
      tol   = cfg%cohort_size_tol_min

      do ifus = 1_ik, cfg%n_cohort_fusion_iter
         call fuse_pass(poly, cfg, tol, force)
         if (.not. force) then
            if (max_cohort_count(poly) <= maxc) exit
         end if
         tol = tol * cfg%cohort_size_tol_mult
      end do
   end subroutine new_fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(poly, cfg, tol, force)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_height, lai_comb, h_top
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(poly%cohort%n))
      alive = .true.
      associate (cohort => poly%cohort, patch => poly%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            if (i1 < i0) cycle
            !----- The tolerance is a share of the patch's TOP height (cohorts are sorted tallest  !
            !      first), so it scales with the stand: in a 3 m regrowth stand two cohorts merge   !
            !      only within about 6 cm, and one that has pulled ahead keeps its place in the     !
            !      light, while a 30 m forest merges within about 0.6 m. Tied to the PFT's maximum  !
            !      height instead, it would fold each month's recruits into a young stand's        !
            !      leaders and re-merge the halves of every split, and a dense stand could never    !
            !      form the size hierarchy it thins by. ------------------------------------------!
            h_top = cohort%height(i0)
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (cohort%pft(donc) /= cohort%pft(recc)) cycle
                  lai_comb = cohort%nplant(recc) * cohort%leaf_area(recc) + cohort%nplant(donc) * cohort%leaf_area(donc)
                  if (.not. force .and. lai_comb >= cfg%cohort_lai_cap) cycle
                  diff_height = abs(cohort%height(donc) - cohort%height(recc))
                  if (force .or. diff_height < h_top * tol) then
                     call fuse_2_cohorts(poly, recc, donc, cfg%conservation_tol)
                     alive(donc) = .false.
                     did_fuse    = .true.
                  end if
               end do
            end do
         end do
      end associate
      if (did_fuse) then
         call compact_cohort_block(poly%cohort, alive)
         call rebuild_cohort_map(poly)
         call sort_cohorts(poly)
      end if
   end subroutine fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor cohort into receptor, conserving plant number and the four carbon pools (so     !
   ! AGB and leaf area). The fused per-plant pools are nplant-weighted means; the diameter and   !
   ! the rest of the geometry are re-derived from them (set_cohort_size_from_carbon).            !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(poly, recc, donc, conservation_tol)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, agb_tot, agb_new, wr, wd, wtot
      real(wp) :: lc_tot, fc_tot, wc_tot, nc_tot                !< pre-merge pool totals [kgC/m2]
      real(wp) :: cr, cd                                        !< the two cohorts' leaf carbon [kgC/m2]
      real(wp) :: hlr, hld, hwr, hwd                            !< their leaf and wood heat capacities [J/m2/K]

      associate (cohort => poly%cohort)
         nr      = cohort%nplant(recc)
         nd      = cohort%nplant(donc)
         ntot    = nr + nd
         agb_tot = nr * cohort%agb(recc) + nd * cohort%agb(donc)     ! [kgC/m2] conserved
         !----- The fast-loop-owned state fuses by its DECLARED per-field policy, in one place  !
         !      (meds_cohort_state_vars%fuse_cohort_fast_state): nplant-weighted for every per-plant !
         !      amount, the tissue energies and films among them; the temperatures follow from the   !
         !      fused energies. Enumerating it here would let a field drift from its                 !
         !      declared convention. The dynamic leaf traits stay here: they are slow-loop state,      !
         !      not fast, and they must be set BEFORE the survivor's geometry is re-derived because    !
         !      sla maps its leaf_carbon <-> leaf_area. Each trait blends on the quantity it is per: !
         !      vcmax25 and rd25 are per leaf area, so leaf area weights them; sla is leaf area per    !
         !      leaf CARBON, so leaf carbon weights it and the fused cohort shows the leaf area its    !
         !      two parents showed; leaf_turnover_rate is the share of the leaf CARBON lost a year,    !
         !      so leaf carbon weights it and the fused cohort turns over the leaf carbon its two      !
         !      parents did (docs/science/plant_traits.md §5). -------------------------------------!
         wr = nr * cohort%leaf_area(recc) ; wd = nd * cohort%leaf_area(donc) ; wtot = wr + wd
         cr = nr * cohort%leaf_carbon(recc) ; cd = nd * cohort%leaf_carbon(donc)
         !----- The tissue heat capacities the diagnostic temperatures fuse on, taken BEFORE          !
         !      fuse_cohort_fast_state blends the survivor's tissue water into them: the same          !
         !      weights the prognostic temperatures fuse on. ----------------------------------------!
         call cohort_tissue_heat_capacity(cohort, recc, hlr, hwr)
         call cohort_tissue_heat_capacity(cohort, donc, hld, hwd)
         call fuse_cohort_fast_state(cohort, recc, donc, nr, nd)
         cohort%sla(recc)     = pair_mean(cohort%sla(recc),     cr, cohort%sla(donc),     cd)
         cohort%vcmax25(recc) = pair_mean(cohort%vcmax25(recc), wr, cohort%vcmax25(donc), wd)
         cohort%rd25(recc)    = pair_mean(cohort%rd25(recc),    wr, cohort%rd25(donc),    wd)
         cohort%leaf_turnover_rate(recc) = pair_mean(cohort%leaf_turnover_rate(recc), cr,                &
                                                     cohort%leaf_turnover_rate(donc), cd)
         !----- Fast-loop DIAGNOSTIC accumulators, handed leaf area (wr, wd), nplant (nr, nd) and the  !
         !      tissue heat capacities. Which weight applies is declared per field in                 !
         !      meds_diag_vars (CO_FUSION_KIND), not decided here. ------------------------------------------!
         call fuse_cohort_diag(cohort%diag,  recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
         call fuse_cohort_diag(cohort%sdiag, recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
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
         call set_cohort_size_from_carbon(poly%cohort, recc)   ! dbh=wood_to_dbh(wood_carbon); pools kept
         !----- The merged tissues' temperatures, from the energies fused above and the merged masses. !
         call set_cohort_tissue_temps(poly%cohort, recc)
         !----- Carbon-pool conservation guard (per pool, density). ------------------------!
         !      The pools are non-negative, so each density total is also its magnitude.        !
         if (.not. (conserved(lc_tot, ntot*cohort%leaf_carbon(recc),          lc_tot, conservation_tol) .and. &
                    conserved(fc_tot, ntot*cohort%fineroot_carbon(recc),      fc_tot, conservation_tol) .and. &
                    conserved(wc_tot, ntot*cohort%wood_carbon(recc),          wc_tot, conservation_tol) .and. &
                    conserved(nc_tot, ntot*cohort%nonstructural_carbon(recc), nc_tot, conservation_tol)))    &
            error stop 'fuse_2_cohorts: carbon-pool conservation violated'
         !----- AGB- and leaf-area-density conservation guards (agb = p_aboveground_frac*wood_carbon  !
         !      and leaf_area = sla*leaf_carbon, with nplant*wood_carbon, nplant*leaf_carbon and the  !
         !      leaf-carbon-weighted sla conserved above). -------------------------------------------!
         agb_new = cohort%nplant(recc) * cohort%agb(recc)
         if (.not. conserved(agb_tot, agb_new, agb_tot, conservation_tol)) then
            write(*,'(a,2es24.16)') ' fuse_2_cohorts: AGB before, after ', agb_tot, agb_new
            error stop 'fuse_2_cohorts: AGB conservation violated'
         end if
         if (.not. conserved(wtot, ntot*cohort%leaf_area(recc), wtot, conservation_tol)) then
            write(*,'(a,2es24.16)') ' fuse_2_cohorts: leaf area before, after ', wtot,          &
                                    ntot*cohort%leaf_area(recc)
            error stop 'fuse_2_cohorts: leaf-area conservation violated'
         end if
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose LAI density exceeds the cap, iterating until none do.          !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: eps, agb_before, agb_after, wc0

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps

      do iter = 1_ik, cfg%n_cohort_fusion_iter
         n0 = poly%cohort%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (poly%cohort%nplant(i) * poly%cohort%leaf_area(i) > cfg%cohort_lai_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         agb_before = sum(poly%cohort%nplant(1:n0) * poly%cohort%agb(1:n0))
         call expand_cohort_block(poly%cohort, n0 + nsplit)
         m = n0
         associate (cohort => poly%cohort)
            do i = 1_ik, n0
               if (cohort%nplant(i) * cohort%leaf_area(i) <= cfg%cohort_lai_cap) cycle
               cohort%nplant(i) = 0.5_wp * cohort%nplant(i)
               m = m + 1_ik
               !----- Perturb the wood_carbon ANCHOR +/-eps (NOT dbh); the leaf/fineroot/            !
               !      nonstructural pools stay identical in both daughters. This conserves ALL four    !
               !      pools + AGB EXACTLY (0.5n*wc*(1+eps)+0.5n*wc*(1-eps)=n*wc) and is immune to the   !
               !      hgt_max cap (no dbh-renorm approximation).  ------------------------------------!
               wc0 = cohort%wood_carbon(i)
               call copy_cohort(cohort, m, i)           ! copies halved nplant + params + pools
               cohort%wood_carbon(i) = wc0 * (1.0_wp + eps)
               cohort%wood_carbon(m) = wc0 * (1.0_wp - eps)
               call set_cohort_size_from_carbon(cohort, i)
               call set_cohort_size_from_carbon(cohort, m)
               !----- Both daughters keep their parent's temperatures: the wood carbon each gained  !
               !      or lost is held at them, and the two together hold the parent's energy. ----!
               call hold_cohort_tissue_temps(cohort, i)
               call hold_cohort_tissue_temps(cohort, m)
            end do
            cohort%n = m
            agb_after  = sum(cohort%nplant(1:m) * cohort%agb(1:m))
         end associate
         !----- The '-eps' daughters (slots n0+1..m) are NEW cohorts -> fresh global ids; the !
         !      '+eps' half kept slot i and its parent id (the continuation).                 !
         do i = n0 + 1_ik, m
            call assign_cohort_id(poly, i)
         end do
         if (.not. conserved(agb_before, agb_after, agb_before, cfg%conservation_tol))             &
            error stop 'split_cohorts: AGB conservation violated'
         call rebuild_cohort_map(poly)
         call sort_cohorts(poly)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! Cull cohorts below the AGB-density floor, and cohorts at or below the absolute density      !
   ! floor. The applier stops a cohort's decline exactly AT negligible_nplant (it never takes a   !
   ! density below it), so a cohort that has died out sits on the floor, not under it: the test  !
   ! must include equality, or such a cohort would be kept for ever. The ENTIRE                   !
   ! remaining carbon of a culled cohort (leaf/fine-root/wood/storage; it is being removed        !
   ! outright, unlike a turnover shed) becomes litter in its own patch's soil-carbon pools when   !
   ! [soil_carbon].soil_carbon_on (default on; docs/science/soil_carbon.md §5) -- added           !
   ! directly onto the named fields, since this module cannot link biogeochemistry                !
   ! (necromass_to_litter takes plain scalars).                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(poly, cfg, water_shed)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      !----- What the cull HANDED ON, for the caller to declare to the slow ledger                !
      !      (meds_slow_ledger). Reported rather than declared here because the demography library  !
      !      does not link the ledger: the engine applies, the driver accounts. -------------------!
      real(wp), optional, intent(out) :: water_shed  !< [kg/m2 site] tissue + film water -> the ground
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n, pf, ip
      real(wp)             :: lab_g, lab_s, str_g, str_s, lig_g, lig_s
      real(wp)             :: w_tot, w_cohort

      w_tot = 0.0_wp
      if (present(water_shed)) water_shed = 0.0_wp
      n = poly%cohort%n
      if (n < 1_ik) return
      allocate(keep(n))
      associate (cohort => poly%cohort, pft => cfg%pft, patch => poly%patch)
         do i = 1_ik, n
            keep(i) = (cohort%nplant(i) * cohort%agb(i) >= cfg%min_cohort_agb) .and.               &
                      (cohort%nplant(i) > cfg%negligible_nplant)
            if (keep(i)) cycle
            ip = cohort%owner_patch(i)
            !----- A culled cohort's WATER goes to the ground, down the same channel turnover      !
            !      shedding already uses -- one verified path, not a second mechanism. Its tissue  !
            !      HEAT is a change in thermal MASS and is declared with the rest of the phase's,  !
            !      by the caller, since fusion and fission change it too and they are one          !
            !      mechanism. ---------------------------------------------------------------------!
            w_cohort = cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i)
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + w_cohort / max(cfg%dt_slow, tiny_num)
            patch%shed_water_enth(ip) = patch%shed_water_enth(ip)                                   &
                                      + cohort%nplant(i) * cohort_water_enthalpy(cohort, i) / max(cfg%dt_slow, tiny_num)
            w_tot = w_tot + patch_polygon_area(poly, ip) * w_cohort
            !----- The CULL pathway's carbon. BEFORE the soil_carbon_on cycle below, for the          !
            !      same reason the water above is: a cohort dropping off the tracked size floor is a  !
            !      demographic event, and it still happened when the soil pools are off. The cull     !
            !      runs inside the cohort restructuring, which permutes the COHORT axis and leaves    !
            !      the patch axis alone, so patch slot `ip` is the one this cohort's diagnostics      !
            !      belong to and still will be after compact_cohort_block.                            !
            !      A rate [kgC/m2/yr] times the step's weight in seconds -- the block's contract.     !
            !      `dt_slow` is the weight one slow step accumulates; the cull fires on a             !
            !      month boundary and the output layer's temporal fold handles the rest.              !
            if (patch%diag%active)                                                                 &
               patch%diag%v(PA_MORT_C_CULL, ip) = patch%diag%v(PA_MORT_C_CULL, ip)                 &
                  + (cohort%nplant(i) * (cohort%leaf_carbon(i) + cohort%fineroot_carbon(i)         &
                                       + cohort%wood_carbon(i) + cohort%nonstructural_carbon(i))   &
                     / cfg%dt_years) * cfg%dt_slow
            if (.not. cfg%soil_carbon_on) cycle
            pf = cohort%pft(i)
            call necromass_to_litter(cohort%nplant(i) * cohort%leaf_carbon(i),                    &
                     cohort%nplant(i) * cohort%fineroot_carbon(i),                                 &
                     cohort%nplant(i) * cohort%wood_carbon(i),                                     &
                     cohort%nplant(i) * cohort%nonstructural_carbon(i),                             &
                     pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                                  &
                     cohort%pft_table%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf),       &
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
      call compact_cohort_block(poly%cohort, keep)
      call rebuild_cohort_map(poly)
   end subroutine terminate_cohorts

   !---------------------------------------------------------------------------------------!
   ! Spawn from the carry-forward recruit pool: when a pool reaches `min_recruit_size`, ONE new  !
   ! cohort is born at the shared minimum cohort height (the smallest tracked size) and the pool  !
   ! is reset; otherwise it carries over, so rare recruiters still establish eventually. A HOST   !
   ! structural process -- it changes the cohort count -- so it lives with the cohort fuse/fission !
   ! housekeeping.                                                                                 !
   !                                                                                          !
   ! The pool is CREDITED DAILY by the driver (accumulate_recruit_pool) and CONSUMED monthly here. !
   ! Crediting it here from the one day's recruitment rate, scaled up to the month, would take a   !
   ! 12-point sample of a quantity the model computes every day and make the pool's carbon         !
   ! depend on which days are month boundaries. ---------------------------------------------------!
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(poly, cfg, carbon_drawn, heat_drawn, water_drawn)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      !----- What a recruit brought IN FROM OUTSIDE the modelled system, for the caller to declare  !
      !      to the slow ledger. Recruitment is a stand-in for everything between a seed and         !
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
      !      endowment is the draw. A recruit is born holding its tissue water (init_cohort), drawn !
      !      from outside as its carbon and heat are, for the same reason. ---------------------------!
      real(wp), optional,  intent(out)   :: carbon_drawn  !< [kgC/m2 site]
      real(wp), optional,  intent(out)   :: heat_drawn    !< [J/m2 site]
      real(wp), optional,  intent(out)   :: water_drawn   !< [kg/m2 site]
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh, c_min, endow, t_birth, psi_birth, c_tot, e_tot, w_tot

      c_tot = 0.0_wp ; e_tot = 0.0_wp ; w_tot = 0.0_wp
      if (present(carbon_drawn)) carbon_drawn = 0.0_wp
      if (present(heat_drawn))   heat_drawn   = 0.0_wp
      if (present(water_drawn))  water_drawn  = 0.0_wp
      np = poly%patch%n
      if (np < 1_ik) return

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, poly%n_pft
            if (poly%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      !----- No pool reached the spawn threshold: nothing is born, so nothing is drawn (the seed    !
      !      rain is declared daily, as it arrives, by accumulate_recruit_pool). -------------------!
      if (nspawn == 0_ik) then
         if (present(carbon_drawn)) carbon_drawn = c_tot
         return
      end if

      call expand_cohort_block(poly%cohort, poly%cohort%n + nspawn)
      n_before = poly%cohort%n
      m = poly%cohort%n
      associate (cohort => poly%cohort, patch => poly%patch, pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, poly%n_pft
               if (patch%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               !----- Every PFT recruits at the smallest tracked height, on its own curve. -------!
               recruit_dbh = height_to_dbh(cohort%pft_table%allom(pf), pft%min_cohort_height)
               !----- Born in its patch's state: at the canopy-air temperature, holding the water     !
               !      its tissues hold at the top soil layer's potential (patch_physical_state). -------!
               call patch_physical_state(patch%cas(ip), patch%soil_water(ip), cfg,                     &
                                         poly%site%surface(patch%owner_site(ip)), t_birth, psi_birth)
               call init_cohort(cohort, m, pft, pf, ip, patch%recruit_pool(pf, ip), recruit_dbh,   &
                                cfg%hydraulics, t_birth, psi_birth)
               !----- The external draw: everything the new cohort holds, less the reproduction    !
               !      carbon the pool actually paid for at `carbon_min` per plant. ---------------!
               c_min = min_cohort_carbon(cohort%pft_table%allom(pf), pft%min_cohort_height, pft%wood_density(pf))
               endow = cohort%leaf_carbon(m) + cohort%fineroot_carbon(m)                           &
                     + cohort%wood_carbon(m) + cohort%nonstructural_carbon(m)
               c_tot = c_tot + patch_polygon_area(poly, ip) * patch%recruit_pool(pf, ip) * (endow - c_min)
               e_tot = e_tot + patch_polygon_area(poly, ip) * cohort%nplant(m)                       &
                             * (cohort%leaf_energy(m) + cohort%wood_energy(m))
               w_tot = w_tot + patch_polygon_area(poly, ip) * cohort_tissue_water(cohort, m)
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      if (present(carbon_drawn)) carbon_drawn = c_tot
      if (present(heat_drawn))   heat_drawn   = e_tot
      if (present(water_drawn))  water_drawn  = w_tot

      !----- Stamp each freshly spawned cohort with a persistent global id. ----------------!
      do k = n_before + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, k)
      end do

      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine apply_recruitment

end module meds_demography_cohort_fusefiss
