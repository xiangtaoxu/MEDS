! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_cohort_fusefiss -- the cohort operators (docs/science/demography.md §3). !
! sort_cohorts keeps each patch's cohorts tallest first. Fusion (new_fuse_cohorts,         !
! fuse_2_cohorts), fission (split_cohorts, for a cohort over the LAI cap) and the cull     !
! (terminate_cohorts) keep the cohort count bounded without changing what the stand does.  !
! apply_recruitment founds new cohorts from the recruit pool. What each conserves and      !
! hands to the slow ledger is demography.md §5; the patch operators are                    !
! meds_demography_patch_fusefiss.                                                          !
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
   ! Sort each patch's cohorts tallest first: height descending, diameter descending on    !
   ! ties. The overtopping LAI and the fusion scans rely on the order (demography.md §3).  !
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
   ! Cohort fusion in passes, the tolerance growing geometrically from                     !
   ! cohort_size_tol_min, until no patch holds more than max_cohort cohorts; a negative    !
   ! max_cohort fuses every same-PFT pair (demography.md §3).                              !
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
            !----- The tolerance is a share of the patch's top height (its first cohort), so it !
            !      scales with the stand: a cohort that pulls ahead in a young stand keeps its  !
            !      place in the light, and a dense stand can form the size hierarchy it thins   !
            !      by (demography.md §3). ------------------------------------------------------!
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
   ! Merge cohort donc into recc. Plant number and the four carbon pools are conserved     !
   ! (each pool a density-weighted mean per plant), and with them biomass and leaf area;   !
   ! the geometry is re-derived from the merged carbon, never averaged (demography.md §3). !
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
         !----- The fast-loop state fuses by the policy fuse_cohort_fast_state declares per field   !
         !      (each per-plant amount, the tissue energies among them, density-weighted). The leaf !
         !      traits are set before the geometry is re-derived, since sla maps leaf carbon to     !
         !      leaf area; each blends on the quantity it is per (plant_traits.md §5). -------------!
         wr = nr * cohort%leaf_area(recc) ; wd = nd * cohort%leaf_area(donc) ; wtot = wr + wd
         cr = nr * cohort%leaf_carbon(recc) ; cd = nd * cohort%leaf_carbon(donc)
         !----- The tissue heat capacities before fuse_cohort_fast_state blends the water into   !
         !      them: the diagnostic temperatures fuse on the weights the prognostic ones do. ---!
         call cohort_tissue_heat_capacity(cohort, recc, hlr, hwr)
         call cohort_tissue_heat_capacity(cohort, donc, hld, hwd)
         call fuse_cohort_fast_state(cohort, recc, donc, nr, nd)
         cohort%sla(recc)     = pair_mean(cohort%sla(recc),     cr, cohort%sla(donc),     cd)
         cohort%vcmax25(recc) = pair_mean(cohort%vcmax25(recc), wr, cohort%vcmax25(donc), wd)
         cohort%rd25(recc)    = pair_mean(cohort%rd25(recc),    wr, cohort%rd25(donc),    wd)
         cohort%leaf_turnover_rate(recc) = pair_mean(cohort%leaf_turnover_rate(recc), cr,                &
                                                     cohort%leaf_turnover_rate(donc), cd)
         !----- The diagnostic rows fuse on leaf area, density or heat capacity, as meds_diag_vars !
         !      declares per field (CO_FUSION_KIND). ----------------------------------------------!
         call fuse_cohort_diag(cohort%diag,  recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
         call fuse_cohort_diag(cohort%sdiag, recc, donc, wr, wd, nr, nd, hlr, hld, hwr, hwd)
         !----- The survivor keeps its own growth history and predawn water status: the two !
         !      cohorts already match in height and LAI, so these differ little. -----------!
         cohort%nplant(recc) = ntot
         !----- The four carbon pools, density-weighted, then the geometry from them;       !
         !      set_cohort_size_from_carbon reads the pools and leaves them as they are. ---!
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
         !----- Each pool's conservation, asserted; the pools are non-negative, so each total is !
         !      also its magnitude. -------------------------------------------------------------!
         if (.not. (conserved(lc_tot, ntot*cohort%leaf_carbon(recc),          lc_tot, conservation_tol) .and. &
                    conserved(fc_tot, ntot*cohort%fineroot_carbon(recc),      fc_tot, conservation_tol) .and. &
                    conserved(wc_tot, ntot*cohort%wood_carbon(recc),          wc_tot, conservation_tol) .and. &
                    conserved(nc_tot, ntot*cohort%nonstructural_carbon(recc), nc_tot, conservation_tol)))    &
            error stop 'fuse_2_cohorts: carbon-pool conservation violated'
         !----- Biomass (a share of the wood carbon) and leaf area (sla x leaf carbon) follow from !
         !      the pools and the leaf-carbon-weighted sla, conserved above; asserted too. --------!
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
   ! Split every cohort whose LAI exceeds cohort_lai_cap into two of half the density,     !
   ! repeating until none does (demography.md §3).                                         !
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
               !----- The wood carbon, not the diameter, moves by +/-eps and the other pools stay !
               !      equal, so every pool and the biomass are conserved exactly, even where the  !
               !      height is capped (demography.md §3). ---------------------------------------!
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
         !----- The '-eps' daughters (slots n0+1..m) are new cohorts and take new ids; the !
         !      '+eps' half keeps slot i and its parent's id. -----------------------------!
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
   ! The cull: remove the cohorts whose n x AGB is below min_cohort_agb, and those at or   !
   ! below the density floor; mortality stops a cohort exactly on the floor, so the test   !
   ! includes equality. A culled cohort's water goes to the ground, and all its carbon to  !
   ! its patch's litter when soil carbon is on (demography.md §3; soil_carbon.md §5).      !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(poly, cfg, water_shed)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      !----- What the cull handed on, for the caller to declare to the slow ledger, which this !
      !      library does not link (conservation.md §3). --------------------------------------!
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
            !----- The culled cohort's water goes to the ground down the channel turnover        !
            !      shedding uses (conservation.md §3). The heat of its tissues is a thermal-mass !
            !      change the caller declares with the phase's (conservation.md §4). ------------!
            w_cohort = cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i)
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + w_cohort / max(cfg%dt_slow, tiny_num)
            patch%shed_water_enth(ip) = patch%shed_water_enth(ip)                                   &
                                      + cohort%nplant(i) * cohort_water_enthalpy(cohort, i) / max(cfg%dt_slow, tiny_num)
            w_tot = w_tot + patch_polygon_area(poly, ip) * w_cohort
            !----- The cull's mortality carbon, a rate [kgC/m2/yr] times the step's seconds, on    !
            !      the cohort's patch, which the cohort restructuring leaves in place              !
            !      (diagnostics.md §4, "The slow rows"). It is written before the soil-carbon test !
            !      because the cull happens with the soil pools off too. --------------------------!
            if (patch%diag%active)                                                                 &
               patch%diag%v(PA_MORT_C_CULL, ip) = patch%diag%v(PA_MORT_C_CULL, ip)                 &
                  + (cohort%nplant(i) * (cohort%leaf_carbon(i) + cohort%fineroot_carbon(i)         &
                                       + cohort%wood_carbon(i) + cohort%nonstructural_carbon(i))   &
                     / cfg%dt_years) * cfg%dt_slow
            if (.not. cfg%soil_carbon%soil_carbon_on) cycle
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
   ! Recruitment: a recruit pool that reaches min_recruit_size founds one cohort at the    !
   ! smallest tracked height and is reset; a smaller pool carries over. The driver credits !
   ! the pool daily (accumulate_recruit_pool) and this consumes it monthly (demography.md  !
   ! §2-§3). It changes the cohort count, so it lives with the cohort operators.           !
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(poly, cfg, carbon_drawn, heat_drawn, water_drawn)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      !----- What the recruits drew from outside the model, for the caller to declare to the   !
      !      slow ledger: their carbon beyond the carbon_min per plant the pool paid for, and  !
      !      their tissues' heat and water. Recruitment stands in for the sizes below          !
      !      min_cohort_height, which the model does not represent, so the draw comes from the !
      !      free atmosphere, not the canopy air (demography.md §3; conservation.md §3). ------!
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

      !----- Count the pools that have reached min_recruit_size. --------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, poly%n_pft
            if (poly%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      !----- No pool reached it: nothing is born, so nothing is drawn (the seed rain is declared !
      !      daily, as it arrives, by accumulate_recruit_pool). ---------------------------------!
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
               call patch_physical_state(patch%cas(ip), patch%soil_water(ip),                          &
                                         poly%site%surface(patch%owner_site(ip)), t_birth, psi_birth)
               call init_cohort(cohort, m, pft, pf, ip, patch%recruit_pool(pf, ip), recruit_dbh,   &
                                t_birth, psi_birth)
               !----- The draw: everything the new cohort holds, less the reproduction carbon the !
               !      pool paid for at carbon_min per plant. -------------------------------------!
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

      !----- Each new cohort takes a new id. ----------------------------------------------!
      do k = n_before + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, k)
      end do

      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine apply_recruitment

end module meds_demography_cohort_fusefiss
