!----- Cohort fusion/fission: conservation of plant number & AGB (carbon), DBH re-derivation. !
program test_fusion_cohort
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : pio4
   use meds_config,          only : meds_config_t
   use meds_site_state_types,           only : site_t, set_cohort_size, cohort_tissue_heat_capacity, &
                                               TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN
   use meds_init,            only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_cohort_fusefiss, only : fuse_2_cohorts, new_fuse_cohorts, split_cohorts,        &
                                         max_cohort_count
   use meds_diagnostic_reduce, only : total_nplant, total_agb
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   real(wp)            :: agb_tot, n0, agb0, dbh_avg
   real(wp)            :: wr, wd, leafmass_exp, woodmass_exp, ltemp_exp, wtemp_exp
   real(wp)            :: cl1, cw1, cl2, cw2, clf, cwf, e_leaf0, e_wood0
   real(wp) :: np1, np2
   integer(ik)         :: j, pf

   call banner('cohort fusion/fission conservation')
   cfg = build_test_config()

   !=== 1. fuse_2_cohorts conserves N and AGB; DBH is re-derived, not averaged. ============!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
   call finalize_init(site)                          ! sorted: index 1 = taller (dbh 12)
   agb_tot = site%cohort%nplant(1)*site%cohort%agb(1) + site%cohort%nplant(2)*site%cohort%agb(2)
   dbh_avg = 0.5_wp*(site%cohort%dbh(1) + site%cohort%dbh(2))
   !----- Seed distinct fast state + predict the merge. The tissue TEMPERATURES are intensive per   !
   !      unit HEAT CAPACITY, so they fuse on the capacity that carries them -- what that conserves  !
   !      is the tissue ENERGY, which is asserted directly below as well as through the formula.     !
   !      They were leaf-area-weighted until 2026-09-09, which is roughly right for leaves and       !
   !      wrong for wood (whose capacity follows wood carbon and the sapwood ring, not leaf area).   !
   !      leaf_water_mass/wood_water_mass are EXTENSIVE like AGB (nplant-weighted, so the fused      !
   !      site-TOTAL water is conserved -- leaf-area-weighting them would not conserve it). ---------!
   site%cohort%leaf_water_mass(1) = 2.0_wp ; site%cohort%leaf_water_mass(2) = 5.0_wp
   site%cohort%wood_water_mass(1) = 3.0_wp ; site%cohort%wood_water_mass(2) = 7.0_wp
   site%cohort%leaf_temp(1) = 300.0_wp ; site%cohort%leaf_temp(2) = 305.0_wp
   site%cohort%wood_temp(1) = 298.0_wp ; site%cohort%wood_temp(2) = 306.0_wp
   !----- The four per-plant flux accumulators are EXTENSIVE too, and were the fields the fusion   !
   !      policy enumerated by hand with an "add any new one to this list" comment. Assert them,   !
   !      so the declared policy in fuse_cohort_fast_state has a witness for every kind it names.  !
   site%cohort%gpp_accum(1)       = 1.0_wp ; site%cohort%gpp_accum(2)       = 4.0_wp
   site%cohort%leaf_resp_accum(1) = 0.2_wp ; site%cohort%leaf_resp_accum(2) = 0.9_wp
   site%cohort%stem_resp_accum(1) = 0.1_wp ; site%cohort%stem_resp_accum(2) = 0.6_wp
   site%cohort%root_resp_accum(1) = 0.3_wp ; site%cohort%root_resp_accum(2) = 0.8_wp
   wr = site%cohort%nplant(1) * site%cohort%leaf_area(1)
   wd = site%cohort%nplant(2) * site%cohort%leaf_area(2)
   np1 = site%cohort%nplant(1) ; np2 = site%cohort%nplant(2)   ! captured BEFORE the fuse
   leafmass_exp = (site%cohort%nplant(1)*2.0_wp + site%cohort%nplant(2)*5.0_wp)                     &
                / (site%cohort%nplant(1) + site%cohort%nplant(2))
   woodmass_exp = (site%cohort%nplant(1)*3.0_wp + site%cohort%nplant(2)*7.0_wp)                     &
                / (site%cohort%nplant(1) + site%cohort%nplant(2))
   call cohort_tissue_heat_capacity(site%cohort, 1_ik, TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN, cl1, cw1)
   call cohort_tissue_heat_capacity(site%cohort, 2_ik, TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN, cl2, cw2)
   ltemp_exp = (cl1*300.0_wp + cl2*305.0_wp) / (cl1 + cl2)
   wtemp_exp = (cw1*298.0_wp + cw2*306.0_wp) / (cw1 + cw2)
   e_leaf0   = cl1*300.0_wp + cl2*305.0_wp
   e_wood0   = cw1*298.0_wp + cw2*306.0_wp
   call fuse_2_cohorts(site, 1_ik, 2_ik, cfg%conservation_tol)
   call check_close(site%cohort%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
   call check_close(site%cohort%nplant(1)*site%cohort%agb(1), agb_tot, 1.0e-12_wp,                  &
                    'fused total AGB not conserved')
   call check_close(site%cohort%basal_area(1), pio4*site%cohort%dbh(1)**2, 1.0e-12_wp,              &
                    'basal area inconsistent with re-derived DBH')
   call check(abs(site%cohort%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')
   call check_close(site%cohort%leaf_water_mass(1), leafmass_exp, 1.0e-12_wp,                       &
                    'leaf_water_mass not nplant-weighted (extensive) on cohort fusion')
   call check_close(site%cohort%wood_water_mass(1), woodmass_exp, 1.0e-12_wp,                       &
                    'wood_water_mass not nplant-weighted (extensive) on cohort fusion')
   call check_close(site%cohort%leaf_temp(1), ltemp_exp, 1.0e-9_wp,                                 &
                    'leaf_temp not heat-capacity-weighted on cohort fusion')
   call check_close(site%cohort%wood_temp(1), wtemp_exp, 1.0e-9_wp,                                 &
                    'wood_temp not heat-capacity-weighted on cohort fusion')
   !----- The formulae above are the MEANS to this end: the merged tissue must hold the energy the  !
   !      two cohorts held. Asserted separately so a future change of weighting is judged on what   !
   !      it conserves, not on whether it matches today's expression. The LEAF capacity is additive !
   !      through the merge (leaf_carbon and leaf_water_mass are both nplant-weighted), so its      !
   !      energy closes exactly; the WOOD's does not, because set_cohort_size_from_carbon re-derives!
   !      the sapwood ring from the merged diameter -- that residual is a thermal-mass change, of   !
   !      the same kind growth makes, and it is bounded here rather than asserted to zero. ---------!
   call cohort_tissue_heat_capacity(site%cohort, 1_ik, TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN, clf, cwf)
   call check_close(clf*site%cohort%leaf_temp(1), e_leaf0, 1.0e-12_wp,                              &
                    'cohort fusion must conserve LEAF tissue energy exactly')
   call check(abs(cwf*site%cohort%wood_temp(1) - e_wood0) <= 1.0e-3_wp * abs(e_wood0),              &
              'cohort fusion must conserve WOOD tissue energy to the sapwood re-derivation')
   block
      real(wp) :: nt
      nt = np1 + np2
      call check_close(site%cohort%gpp_accum(1),       (np1*1.0_wp + np2*4.0_wp)/nt, 1.0e-12_wp,      &
                       'gpp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(site%cohort%leaf_resp_accum(1), (np1*0.2_wp + np2*0.9_wp)/nt, 1.0e-12_wp,      &
                       'leaf_resp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(site%cohort%stem_resp_accum(1), (np1*0.1_wp + np2*0.6_wp)/nt, 1.0e-12_wp,      &
                       'stem_resp_accum not nplant-weighted (extensive) on cohort fusion')
      call check_close(site%cohort%root_resp_accum(1), (np1*0.3_wp + np2*0.8_wp)/nt, 1.0e-12_wp,      &
                       'root_resp_accum not nplant-weighted (extensive) on cohort fusion')
   end block

   !=== 2. new_fuse_cohorts reduces count to <= max_cohort, conserving N and AGB. ==========!
   call init_bare_ground(site, cfg, 1_ik)
   do j = 1_ik, 200_ik
      !----- Low density so per-cohort LAI stays well under cohort_lai_cap (realistic). ---!
      call add_cohort(site, cfg, 1_ik, 1_ik, 2.0e-4_wp, 2.0_wp + 0.55_wp*real(j-1_ik, wp))  ! 2..111 cm
   end do
   call finalize_init(site)
   n0   = total_nplant(site)
   agb0 = total_agb(site)
   call new_fuse_cohorts(site, cfg)
   call check(max_cohort_count(site) <= abs(cfg%max_cohort), 'fusion did not reach max_cohort')
   call check_close(total_nplant(site), n0,   1.0e-10_wp, 'fusion broke nplant conservation')
   call check_close(total_agb(site),    agb0, cfg%conservation_tol, 'fusion broke AGB conservation')

   !=== 3. split_cohorts conserves N and AGB and increases the count. =====================!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)   ! nplant*leaf_area >> cohort_lai_cap
   call finalize_init(site)
   n0   = total_nplant(site)
   agb0 = total_agb(site)
   !----- REVIEW 2026-09 (item 1B #1): the interception films are per m2 GROUND. A split that copies !
   !      them verbatim into both daughters doubles the patch's film water. -------------------------!
   site%cohort%leaf_surf_water(1) = 0.30_wp ; site%cohort%wood_surf_water(1) = 0.05_wp
   call split_cohorts(site, cfg)
   call check(site%cohort%n >= 2_ik, 'split did not create a second cohort')
   call check_close(total_nplant(site), n0,   1.0e-12_wp, 'split broke nplant conservation')
   call check_close(total_agb(site),    agb0, cfg%conservation_tol, 'split broke AGB conservation')
   call check_close(sum(site%cohort%leaf_surf_water(1:site%cohort%n) + site%cohort%wood_surf_water(1:site%cohort%n)), &
                    0.35_wp, 1.0e-12_wp, 'split broke canopy film water conservation (ground-referenced field)')

   !=== 4. Carbon pools/traits thread correctly through the sort reorder (PR3 lockstep). =====!
   !     3 PFTs (sla 16/13/10) are inserted then height-sorted (reordered); re-deriving with   !
   !     each cohort's THREADED p_sla / p_aboveground_frac must reproduce its own PFT's         !
   !     on-allometry pools -- a mis-threaded per-cohort trait would misalign them.             !
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 1.0e-3_wp,  8.0_wp)   ! PFT 1 (shortest)
   call add_cohort(site, cfg, 1_ik, 2_ik, 1.0e-3_wp, 40.0_wp)   ! PFT 2 (tallest -> sorts first)
   call add_cohort(site, cfg, 1_ik, 3_ik, 1.0e-3_wp, 20.0_wp)   ! PFT 3
   call finalize_init(site)                                     ! sort_cohorts reorders EVERY field
   do j = 1_ik, site%cohort%n
      call set_cohort_size(site%cohort, j)                      ! re-derive with the threaded traits
      pf = site%cohort%pft(j)
      call check_close(site%cohort%leaf_carbon(j) * cfg%pft%sla(pf), site%cohort%leaf_area(j),        &
                       1.0e-9_wp, 'leaf_carbon*sla(pft) /= leaf_area after reorder (trait mis-threaded?)')
      call check_close(site%cohort%wood_carbon(j) * cfg%pft%aboveground_frac(pf), site%cohort%agb(j), &
                       1.0e-9_wp, 'wood_carbon*aboveground_frac(pft) /= agb after reorder')
   end do

   !=== 5. CARBON-mode fusion conserves the FOUR prognostic pools (BUG2): the storage buffer is  !
   !     nplant-weighted and CARRIED, not snapped back onto allometry (set_cohort_size_from_carbon !
   !     path). Pre-fix, fuse called the empirical set_cohort_size and reset the pools. ===========!
   block
      real(wp) :: lc_t, fc_t, wc_t, nc_t, nr5, nd5, alloc_store
      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
      call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
      call finalize_init(site)
      !----- Push one cohort's storage OFF allometry (the starvation buffer the bug destroyed). --!
      site%cohort%nonstructural_carbon(1) = 3.0_wp * site%cohort%nonstructural_carbon(1)
      nr5 = site%cohort%nplant(1) ; nd5 = site%cohort%nplant(2)
      lc_t = nr5*site%cohort%leaf_carbon(1)          + nd5*site%cohort%leaf_carbon(2)
      fc_t = nr5*site%cohort%fineroot_carbon(1)      + nd5*site%cohort%fineroot_carbon(2)
      wc_t = nr5*site%cohort%wood_carbon(1)          + nd5*site%cohort%wood_carbon(2)
      nc_t = nr5*site%cohort%nonstructural_carbon(1) + nd5*site%cohort%nonstructural_carbon(2)
      call fuse_2_cohorts(site, 1_ik, 2_ik, cfg%conservation_tol)
      call check_close(site%cohort%nplant(1)*site%cohort%leaf_carbon(1),          lc_t, 1.0e-9_wp, &
                       'carbon fuse: leaf_carbon not conserved')
      call check_close(site%cohort%nplant(1)*site%cohort%fineroot_carbon(1),      fc_t, 1.0e-9_wp, &
                       'carbon fuse: fineroot_carbon not conserved')
      call check_close(site%cohort%nplant(1)*site%cohort%wood_carbon(1),          wc_t, 1.0e-9_wp, &
                       'carbon fuse: wood_carbon not conserved')
      call check_close(site%cohort%nplant(1)*site%cohort%nonstructural_carbon(1), nc_t, 1.0e-9_wp, &
                       'carbon fuse: nonstructural (storage) not conserved')
      !----- Prove storage was CARRIED, not reset to the allometric storage_cushion*leaf_carbon. -!
      alloc_store = cfg%pft%storage_cushion(1) * site%cohort%leaf_carbon(1)
      call check(abs(site%cohort%nonstructural_carbon(1) - alloc_store) > 1.0e-6_wp,               &
                 'carbon fuse must CARRY off-allometry storage, not snap it to allometry (BUG2)')
   end block

   write(*,'(a)') '   PASS'
end program test_fusion_cohort
