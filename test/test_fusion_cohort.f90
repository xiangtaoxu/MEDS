! SPDX-License-Identifier: Apache-2.0
!----- Cohort fusion/fission: conservation of plant number & AGB (carbon), DBH re-derivation. !
program test_fusion_cohort
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : pio4
   use meds_config,          only : meds_config_t
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use meds_cohort_state_vars,          only : set_cohort_size, cohort_tissue_heat_capacity, set_cohort_size_from_carbon, &
                                               cohort_bad_field, cohort_field_name, cohort_deriv_block, cohort_deriv_alloc, &
                                               init_cohort
   use meds_polygon_state_vars,         only : polygon_t
   use meds_demography_update, only : fill_cohort_deriv
   use meds_init,            only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_cohort_fusefiss, only : fuse_2_cohorts, new_fuse_cohorts, split_cohorts,        &
                                         max_cohort_count
   use meds_diagnostic_reduce, only : total_nplant, total_agb
   use meds_diag_vars,        only : alloc_cohort_diag, cohort_diag_value, CO_LEAF_TEMP, CO_WOOD_TEMP
   use meds_allometry,       only : storage_carbon_target
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(polygon_t)  :: poly
   real(wp)            :: agb_tot, n0, agb0, dbh_avg
   real(wp)            :: wr, wd, leafmass_exp, woodmass_exp, ltemp_exp, wtemp_exp
   real(wp)            :: cl1, cw1, cl2, cw2, clf, cwf, e_leaf0, e_wood0
   real(wp) :: np1, np2, xdiag(4), film0
   integer(ik) :: ndiag
   integer(ik)         :: j, pf

   call banner('cohort fusion/fission conservation')
   cfg = build_test_config()

   !=== 1. fuse_2_cohorts conserves N and AGB; DBH is re-derived, not averaged. ============!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
   call finalize_init(poly)                          ! sorted: index 1 = taller (dbh 12)
   agb_tot = poly%cohort%nplant(1)*poly%cohort%agb(1) + poly%cohort%nplant(2)*poly%cohort%agb(2)
   dbh_avg = 0.5_wp*(poly%cohort%dbh(1) + poly%cohort%dbh(2))
   !----- Seed distinct fast state + predict the merge. The tissue TEMPERATURES are intensive per   !
   !      unit HEAT CAPACITY, so they fuse on the capacity that carries them -- what that conserves  !
   !      is the tissue ENERGY, which is asserted directly below as well as through the formula.     !
   !      A leaf-area weighting is roughly right for leaves and wrong for wood (whose capacity       !
   !      follows wood carbon and the sapwood ring, not leaf area).                                  !
   !      leaf_water_mass/wood_water_mass are EXTENSIVE like AGB (nplant-weighted, so the fused      !
   !      site-TOTAL water is conserved -- leaf-area-weighting them would not conserve it). ---------!
   poly%cohort%leaf_water_mass(1) = 2.0_wp ; poly%cohort%leaf_water_mass(2) = 5.0_wp
   poly%cohort%wood_water_mass(1) = 3.0_wp ; poly%cohort%wood_water_mass(2) = 7.0_wp
   poly%cohort%leaf_temp(1) = 300.0_wp ; poly%cohort%leaf_temp(2) = 305.0_wp
   poly%cohort%wood_temp(1) = 298.0_wp ; poly%cohort%wood_temp(2) = 306.0_wp
   !----- The four per-plant flux accumulators are EXTENSIVE too. Assert them, so the declared     !
   !      policy in fuse_cohort_fast_state has a witness for every kind it names.                  !
   poly%cohort%gpp_accum(1)       = 1.0_wp ; poly%cohort%gpp_accum(2)       = 4.0_wp
   poly%cohort%leaf_resp_accum(1) = 0.2_wp ; poly%cohort%leaf_resp_accum(2) = 0.9_wp
   poly%cohort%stem_resp_accum(1) = 0.1_wp ; poly%cohort%stem_resp_accum(2) = 0.6_wp
   poly%cohort%root_resp_accum(1) = 0.3_wp ; poly%cohort%root_resp_accum(2) = 0.8_wp
   !----- And the canopy films [kg/plant]: nplant-weighted like the tissue water. ------------------!
   poly%cohort%leaf_surf_water(1) = 0.02_wp ; poly%cohort%leaf_surf_water(2) = 0.06_wp
   poly%cohort%wood_surf_water(1) = 0.01_wp ; poly%cohort%wood_surf_water(2) = 0.00_wp
   wr = poly%cohort%nplant(1) * poly%cohort%leaf_area(1)
   wd = poly%cohort%nplant(2) * poly%cohort%leaf_area(2)
   np1 = poly%cohort%nplant(1) ; np2 = poly%cohort%nplant(2)   ! captured BEFORE the fuse
   leafmass_exp = (poly%cohort%nplant(1)*2.0_wp + poly%cohort%nplant(2)*5.0_wp)                     &
                / (poly%cohort%nplant(1) + poly%cohort%nplant(2))
   woodmass_exp = (poly%cohort%nplant(1)*3.0_wp + poly%cohort%nplant(2)*7.0_wp)                     &
                / (poly%cohort%nplant(1) + poly%cohort%nplant(2))
   call cohort_tissue_heat_capacity(poly%cohort, 1_ik, cl1, cw1)
   call cohort_tissue_heat_capacity(poly%cohort, 2_ik, cl2, cw2)
   ltemp_exp = (cl1*300.0_wp + cl2*305.0_wp) / (cl1 + cl2)
   wtemp_exp = (cw1*298.0_wp + cw2*306.0_wp) / (cw1 + cw2)
   e_leaf0   = cl1*300.0_wp + cl2*305.0_wp
   e_wood0   = cw1*298.0_wp + cw2*306.0_wp
   !----- The cohort DIAGNOSTIC block, holding one sample of the same tissue temperatures, must fuse  !
   !      them on the same weights as the prognostic state, so the two agree after the merge. --------!
   call alloc_cohort_diag(poly%cohort%diag, poly%cohort%cap, .true.)
   poly%cohort%diag%n = 2_ik
   poly%cohort%diag%v(CO_LEAF_TEMP, 1:2) = [300.0_wp, 305.0_wp]
   poly%cohort%diag%v(CO_WOOD_TEMP, 1:2) = [298.0_wp, 306.0_wp]
   poly%cohort%diag%w(1:2) = 1.0_wp
   call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%conservation_tol)
   call check_close(poly%cohort%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
   call check_close(poly%cohort%nplant(1)*poly%cohort%agb(1), agb_tot, 1.0e-12_wp,                  &
                    'fused total AGB not conserved')
   call check_close(poly%cohort%basal_area(1), pio4*poly%cohort%dbh(1)**2, 1.0e-12_wp,              &
                    'basal area inconsistent with re-derived DBH')
   call check(abs(poly%cohort%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')
   call check_close(poly%cohort%leaf_water_mass(1), leafmass_exp, 1.0e-12_wp,                       &
                    'leaf_water_mass not nplant-weighted (extensive) on cohort fusion')
   call check_close(poly%cohort%wood_water_mass(1), woodmass_exp, 1.0e-12_wp,                       &
                    'wood_water_mass not nplant-weighted (extensive) on cohort fusion')
   call check_close(poly%cohort%leaf_surf_water(1), (np1*0.02_wp + np2*0.06_wp) / (np1 + np2), 1.0e-12_wp, &
                    'leaf_surf_water not nplant-weighted (extensive) on cohort fusion')
   call check_close(poly%cohort%wood_surf_water(1), np1*0.01_wp / (np1 + np2), 1.0e-12_wp,          &
                    'wood_surf_water not nplant-weighted (extensive) on cohort fusion')
   call check_close(poly%cohort%leaf_temp(1), ltemp_exp, 1.0e-9_wp,                                 &
                    'leaf_temp not heat-capacity-weighted on cohort fusion')
   call check_close(poly%cohort%wood_temp(1), wtemp_exp, 1.0e-9_wp,                                 &
                    'wood_temp not heat-capacity-weighted on cohort fusion')
   call cohort_diag_value(poly%cohort%diag, CO_LEAF_TEMP, xdiag, ndiag, 0.0_wp)
   call check_close(xdiag(1), poly%cohort%leaf_temp(1), 1.0e-9_wp,                                  &
                    'diagnostic leaf temperature fuses as the prognostic one (leaf heat capacity)')
   call cohort_diag_value(poly%cohort%diag, CO_WOOD_TEMP, xdiag, ndiag, 0.0_wp)
   call check_close(xdiag(1), poly%cohort%wood_temp(1), 1.0e-9_wp,                                  &
                    'diagnostic wood temperature fuses as the prognostic one (wood heat capacity)')
   !----- The formulae above are the MEANS to this end: the merged tissue must hold the energy the  !
   !      two cohorts held. Asserted separately so a future change of weighting is judged on what   !
   !      it conserves, not on whether it matches today's expression. The LEAF capacity is additive !
   !      through the merge (leaf_carbon and leaf_water_mass are both nplant-weighted), so its      !
   !      energy closes exactly; the WOOD's does not, because set_cohort_size_from_carbon re-derives!
   !      the sapwood ring from the merged diameter -- that residual is a thermal-mass change, of   !
   !      the same kind growth makes, and it is bounded here rather than asserted to zero. ---------!
   call cohort_tissue_heat_capacity(poly%cohort, 1_ik, clf, cwf)
   call check_close(clf*poly%cohort%leaf_temp(1), e_leaf0, 1.0e-12_wp,                              &
                    'cohort fusion must conserve LEAF tissue energy exactly')
   call check(abs(cwf*poly%cohort%wood_temp(1) - e_wood0) <= 1.0e-3_wp * abs(e_wood0),              &
              'cohort fusion must conserve WOOD tissue energy to the sapwood re-derivation')
   block
      real(wp) :: nt
      nt = np1 + np2
      call check_close(poly%cohort%gpp_accum(1),       (np1*1.0_wp + np2*4.0_wp)/nt, 1.0e-12_wp,      &
                       'gpp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%leaf_resp_accum(1), (np1*0.2_wp + np2*0.9_wp)/nt, 1.0e-12_wp,      &
                       'leaf_resp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%stem_resp_accum(1), (np1*0.1_wp + np2*0.6_wp)/nt, 1.0e-12_wp,      &
                       'stem_resp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(poly%cohort%root_resp_accum(1), (np1*0.3_wp + np2*0.8_wp)/nt, 1.0e-12_wp,      &
                       'root_resp_accum not nplant-weighted (extensive) on cohort fusion')
   end block

   !=== 2. new_fuse_cohorts reduces count to <= max_cohort, conserving N and AGB. ==========!
   call init_bare_ground(poly, cfg, 1_ik)
   do j = 1_ik, 200_ik
      !----- Low density so per-cohort LAI stays well under cohort_lai_cap (realistic). ---!
      call add_cohort(poly, cfg, 1_ik, 1_ik, 2.0e-4_wp, 2.0_wp + 0.55_wp*real(j-1_ik, wp))  ! 2..111 cm
   end do
   call finalize_init(poly)
   n0   = total_nplant(poly)
   agb0 = total_agb(poly)
   call new_fuse_cohorts(poly, cfg)
   call check(max_cohort_count(poly) <= abs(cfg%max_cohort), 'fusion did not reach max_cohort')
   call check_close(total_nplant(poly), n0,   1.0e-10_wp, 'fusion broke nplant conservation')
   call check_close(total_agb(poly),    agb0, cfg%conservation_tol, 'fusion broke AGB conservation')

   !=== 2b. The fusion tolerance is a share of the patch's top height. Two saplings 1.0 and 1.1 cm  !
   !        wide (about 0.15 m apart) stay apart in a stand of saplings, where 0.02 of the top is   !
   !        about 5 cm, and merge once the patch holds a 60 cm tree, where it is about 0.6 m. -----!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.1_wp)
   call finalize_init(poly)
   call new_fuse_cohorts(poly, cfg)
   call check(poly%cohort%n == 2_ik, 'saplings 0.15 m apart merged in a sapling stand')
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.05_wp, 1.1_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0e-3_wp, 60.0_wp)
   call finalize_init(poly)
   call new_fuse_cohorts(poly, cfg)
   call check(poly%cohort%n == 2_ik, 'saplings 0.15 m apart not merged under a 60 cm tree')

   !=== 3. split_cohorts conserves N and AGB and increases the count. =====================!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)   ! nplant*leaf_area >> cohort_lai_cap
   call finalize_init(poly)
   n0   = total_nplant(poly)
   agb0 = total_agb(poly)
   !----- The interception films are per plant: each daughter keeps its parent's, so the patch's     !
   !      film water is what it was. -------------------------------------------------------------!
   poly%cohort%leaf_surf_water(1) = 0.30_wp ; poly%cohort%wood_surf_water(1) = 0.05_wp
   film0 = poly%cohort%nplant(1) * 0.35_wp
   call split_cohorts(poly, cfg)
   call check(poly%cohort%n >= 2_ik, 'split did not create a second cohort')
   call check_close(total_nplant(poly), n0,   1.0e-12_wp, 'split broke nplant conservation')
   call check_close(total_agb(poly),    agb0, cfg%conservation_tol, 'split broke AGB conservation')
   call check_close(sum(poly%cohort%nplant(1:poly%cohort%n) * (poly%cohort%leaf_surf_water(1:poly%cohort%n)        &
                                                             + poly%cohort%wood_surf_water(1:poly%cohort%n))),   &
                    film0, 1.0e-12_wp, 'split broke canopy film water conservation')

   !=== 4. Carbon pools/traits thread correctly through the lockstep sort reorder. ==========!
   !     3 PFTs (sla 16/13/10) are inserted then height-sorted (reordered); re-deriving with   !
   !     each cohort's THREADED p_sla / p_aboveground_frac must reproduce its own PFT's         !
   !     on-allometry pools -- a mis-threaded per-cohort trait would misalign them.             !
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0e-3_wp,  8.0_wp)   ! PFT 1 (shortest)
   call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0e-3_wp, 40.0_wp)   ! PFT 2 (tallest -> sorts first)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 1.0e-3_wp, 20.0_wp)   ! PFT 3
   call finalize_init(poly)                                     ! sort_cohorts reorders EVERY field
   do j = 1_ik, poly%cohort%n
      call set_cohort_size(poly%cohort, j)                      ! re-derive with the threaded traits
      pf = poly%cohort%pft(j)
      call check_close(poly%cohort%leaf_carbon(j) * cfg%pft%sla(pf), poly%cohort%leaf_area(j),        &
                       1.0e-9_wp, 'leaf_carbon*sla(pft) /= leaf_area after reorder (trait mis-threaded?)')
      call check_close(poly%cohort%wood_carbon(j) * cfg%pft%aboveground_frac(pf), poly%cohort%agb(j), &
                       1.0e-9_wp, 'wood_carbon*aboveground_frac(pft) /= agb after reorder')
   end do

   !=== 5. CARBON-mode fusion conserves the FOUR prognostic pools: the storage buffer is         !
   !     nplant-weighted and CARRIED, not snapped back onto allometry (set_cohort_size_from_carbon !
   !     path); the empirical set_cohort_size would reset the pools from size. ====================!
   block
      real(wp) :: lc_t, fc_t, wc_t, nc_t, nr5, nd5, alloc_store
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
      call finalize_init(poly)
      !----- Push one cohort's storage OFF allometry (the starvation buffer a reset would lose). -!
      poly%cohort%nonstructural_carbon(1) = 3.0_wp * poly%cohort%nonstructural_carbon(1)
      nr5 = poly%cohort%nplant(1) ; nd5 = poly%cohort%nplant(2)
      lc_t = nr5*poly%cohort%leaf_carbon(1)          + nd5*poly%cohort%leaf_carbon(2)
      fc_t = nr5*poly%cohort%fineroot_carbon(1)      + nd5*poly%cohort%fineroot_carbon(2)
      wc_t = nr5*poly%cohort%wood_carbon(1)          + nd5*poly%cohort%wood_carbon(2)
      nc_t = nr5*poly%cohort%nonstructural_carbon(1) + nd5*poly%cohort%nonstructural_carbon(2)
      call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%conservation_tol)
      call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1),          lc_t, 1.0e-9_wp, &
                       'carbon fuse: leaf_carbon not conserved')
      call check_close(poly%cohort%nplant(1)*poly%cohort%fineroot_carbon(1),      fc_t, 1.0e-9_wp, &
                       'carbon fuse: fineroot_carbon not conserved')
      call check_close(poly%cohort%nplant(1)*poly%cohort%wood_carbon(1),          wc_t, 1.0e-9_wp, &
                       'carbon fuse: wood_carbon not conserved')
      call check_close(poly%cohort%nplant(1)*poly%cohort%nonstructural_carbon(1), nc_t, 1.0e-9_wp, &
                       'carbon fuse: nonstructural (storage) not conserved')
      !----- Prove storage was CARRIED, not reset to the allometric storage target. -------------!
      alloc_store = storage_carbon_target(cfg%pft%storage_cushion(1), poly%cohort%leaf_carbon(1),   &
                                          poly%cohort%fineroot_carbon(1))
      call check(abs(poly%cohort%nonstructural_carbon(1) - alloc_store) > 1.0e-6_wp,               &
                 'carbon fuse must CARRY off-allometry storage, not snap it to allometry')
   end block

   !=== 6. Fusion of cohorts with different SLA conserves leaf AREA: sla blends on leaf          !
   !     carbon, the quantity it converts, so sum(nplant*leaf_carbon*sla) is what the fused        !
   !     cohort shows. A leaf-area-weighted sla would give sum(LA*sla)/sum(LA) instead. The leaf    !
   !     turnover rate is the share of the leaf carbon lost a year, so it too blends on leaf        !
   !     carbon, and the fused cohort loses the leaf carbon its two parents did; with the two      !
   !     SLAs differing, a leaf-area weight would not. =============================================!
   block
      real(wp) :: lai0, loss0
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
      call finalize_init(poly)
      poly%cohort%sla(1) = 12.0_wp ; poly%cohort%sla(2) = 24.0_wp      ! shade-acclimated donor
      poly%cohort%leaf_turnover_rate(1) = 1.0_wp ; poly%cohort%leaf_turnover_rate(2) = 0.4_wp
      call set_cohort_size_from_carbon(poly%cohort, 1_ik)
      call set_cohort_size_from_carbon(poly%cohort, 2_ik)
      lai0 = poly%cohort%nplant(1)*poly%cohort%leaf_area(1) + poly%cohort%nplant(2)*poly%cohort%leaf_area(2)
      loss0 = poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1)*poly%cohort%leaf_turnover_rate(1)       &
            + poly%cohort%nplant(2)*poly%cohort%leaf_carbon(2)*poly%cohort%leaf_turnover_rate(2)
      call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%conservation_tol)
      call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_area(1), lai0, 1.0e-12_wp,              &
                       'fusion of cohorts with different SLA must conserve leaf area')
      call check_close(poly%cohort%nplant(1)*poly%cohort%leaf_carbon(1)*poly%cohort%leaf_turnover_rate(1), &
                       loss0, 1.0e-12_wp, 'fusion must conserve the leaf carbon turned over a year')
      !----- The same at any size: two LEAFLESS cohorts, whose leaf carbon has decayed to ~1e-193   !
      !      (as in a BCI run). A blend that skipped weights below 1e-30 would keep one SLA here. -!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
      call finalize_init(poly)
      poly%cohort%sla(1) = 12.0_wp ; poly%cohort%sla(2) = 24.0_wp
      poly%cohort%leaf_carbon(1) = 1.0e-193_wp ; poly%cohort%leaf_carbon(2) = 3.0e-193_wp
      call set_cohort_size_from_carbon(poly%cohort, 1_ik)
      call set_cohort_size_from_carbon(poly%cohort, 2_ik)
      lai0 = poly%cohort%nplant(1)*poly%cohort%leaf_area(1) + poly%cohort%nplant(2)*poly%cohort%leaf_area(2)
      call fuse_2_cohorts(poly, 1_ik, 2_ik, cfg%conservation_tol)
      call check(abs(poly%cohort%nplant(1)*poly%cohort%leaf_area(1) - lai0) <= 1.0e-12_wp * lai0,      &
                 'fusion of leafless cohorts with different SLA must conserve their leaf area too')
   end block

   !=== 7. Geometry has ONE writer: a cohort placed on the allometry from its dbh carries        !
   !     exactly the geometry its carbon gives, so a restart that re-derives it gets the same bits. !
   block
      real(wp) :: g0(5)
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.4_wp, 17.3_wp)
      call finalize_init(poly)
      g0 = [poly%cohort%dbh(1), poly%cohort%height(1), poly%cohort%basal_area(1),                    &
            poly%cohort%agb(1), poly%cohort%leaf_area(1)]
      call set_cohort_size_from_carbon(poly%cohort, 1_ik)
      call check(all([poly%cohort%dbh(1), poly%cohort%height(1), poly%cohort%basal_area(1),          &
                      poly%cohort%agb(1), poly%cohort%leaf_area(1)] == g0),                          &
                 'birth geometry must be the geometry derived from the birth carbon, bit for bit')

      !=== 8. The bounds guard names the first field out of bounds, and passes a sound cohort. ==!
      call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'a fresh cohort is in bounds')
      poly%cohort%leaf_surf_water(1) = -1.0e-30_wp
      call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'leaf_surf_water',  &
                 'a negative leaf film is out of bounds, however small')
      poly%cohort%nplant(1) = ieee_value(1.0_wp, ieee_quiet_nan)
      call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nplant',           &
                 'a NaN density is out of bounds, and checked first')
   end block

   !=== 9. The growth window's sum is the buffer's: once a large sample has been evicted, a       !
   !     window of zeros averages exactly zero; and a recruit in a reused slot starts empty. ======!
   block
      type(cohort_deriv_block) :: deriv
      integer(ik) :: k, pos, nwin
      real(wp)    :: sample
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.4_wp, 17.3_wp)
      call finalize_init(poly)
      nwin = int(size(poly%cohort%growth_hist, 1), ik)
      call cohort_deriv_alloc(deriv, poly%cohort%n)
      pos = 0_ik
      do k = 1_ik, 3_ik * nwin
         pos = mod(pos, nwin) + 1_ik
         sample = 0.0_wp
         if (k <= 2_ik * nwin) sample = 0.1_wp
         if (k == 1_ik)        sample = 1.0e8_wp
         call fill_cohort_deriv(poly%cohort, 1_ik, deriv, cfg%dt_years, 0.0_wp, sample,             &
                                poly%cohort%leaf_carbon(1), poly%cohort%fineroot_carbon(1),         &
                                poly%cohort%wood_carbon(1), poly%cohort%nonstructural_carbon(1),    &
                                nwin, pos)
      end do
      call check(poly%cohort%growth_avg(1) == 0.0_wp,                                                &
                 'a window of zero growth averages exactly zero after a large sample left it')
      call check(poly%cohort%growth_accum(1) == sum(poly%cohort%growth_hist(1:nwin, 1)),             &
                 'the window''s sum is the buffer''s')
      call init_cohort(poly%cohort, 1_ik, cfg%pft, 1_ik, 1_ik, 0.1_wp, 1.0_wp, cfg%hydraulics,      &
                       288.15_wp, -0.1_wp)                       ! any birth state: the window is the point
      call check(all(poly%cohort%growth_hist(1:nwin, 1) == 0.0_wp) .and.                             &
                 poly%cohort%growth_count(1) == 0_ik, 'a recruit in a reused slot starts with an empty window')
   end block

   write(*,'(a)') '   PASS'
end program test_fusion_cohort
