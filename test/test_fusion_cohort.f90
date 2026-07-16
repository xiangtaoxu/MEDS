!----- Cohort fusion/fission: conservation of plant number & AGB (carbon), DBH re-derivation. !
program test_fusion_cohort
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : pio4
   use meds_config,          only : meds_config_t
   use meds_ecosystem_state,           only : site_t, set_cohort_size
   use meds_init,            only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_fusefiss, only : fuse_2_cohorts, new_fuse_cohorts, split_cohorts,        &
                                         max_cohort_count
   use meds_demography_diagnostics, only : total_nplant, total_agb
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   real(wp)            :: agb_tot, n0, agb0, dbh_avg
   real(wp)            :: wr, wd, psi_exp, ltemp_exp
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
   !----- Seed distinct fast state + predict the leaf-area-weighted merge (weights fixed now). -!
   site%cohort%psi(1,1)     = -0.5_wp ; site%cohort%psi(1,2)     = -1.5_wp
   site%cohort%leaf_temp(1) = 300.0_wp ; site%cohort%leaf_temp(2) = 305.0_wp
   wr = site%cohort%nplant(1) * site%cohort%leaf_area(1)
   wd = site%cohort%nplant(2) * site%cohort%leaf_area(2)
   psi_exp   = (wr*(-0.5_wp)  + wd*(-1.5_wp))  / (wr + wd)
   ltemp_exp = (wr*300.0_wp   + wd*305.0_wp)   / (wr + wd)
   call fuse_2_cohorts(site, 1_ik, 2_ik, cfg%conservation_tol)
   call check_close(site%cohort%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
   call check_close(site%cohort%nplant(1)*site%cohort%agb(1), agb_tot, 1.0e-12_wp,                  &
                    'fused total AGB not conserved')
   call check_close(site%cohort%basal_area(1), pio4*site%cohort%dbh(1)**2, 1.0e-12_wp,              &
                    'basal area inconsistent with re-derived DBH')
   call check(abs(site%cohort%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')
   call check_close(site%cohort%psi(1,1),     psi_exp,   1.0e-12_wp, 'psi not leaf-area-weighted on cohort fusion')
   call check_close(site%cohort%leaf_temp(1), ltemp_exp, 1.0e-9_wp,  'leaf_temp not leaf-area-weighted on cohort fusion')

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
   call split_cohorts(site, cfg)
   call check(site%cohort%n >= 2_ik, 'split did not create a second cohort')
   call check_close(total_nplant(site), n0,   1.0e-12_wp, 'split broke nplant conservation')
   call check_close(total_agb(site),    agb0, cfg%conservation_tol, 'split broke AGB conservation')

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
