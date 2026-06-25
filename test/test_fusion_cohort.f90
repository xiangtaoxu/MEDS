!----- Cohort fusion/fission: conservation of plant number & AGB (carbon), DBH re-derivation. !
program test_fusion_cohort
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : pio4
   use meds_config,          only : meds_config_t, build_config
   use meds_demography_types,           only : site_t
   use meds_setup,           only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_structure, only : fuse_2_cohorts, new_fuse_cohorts, split_cohorts,        &
                                         max_cohort_count
   use meds_demography_diagnostics, only : total_nplant, total_agb
   use meds_test_support,    only : check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   real(wp)            :: agb_tot, n0, agb0, dbh_avg
   integer(ik)         :: j

   call banner('cohort fusion/fission conservation')
   cfg = build_config()

   !=== 1. fuse_2_cohorts conserves N and AGB; DBH is re-derived, not averaged. ============!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
   call finalize_init(site)                          ! sorted: index 1 = taller (dbh 12)
   agb_tot = site%cohort%nplant(1)*site%cohort%agb(1) + site%cohort%nplant(2)*site%cohort%agb(2)
   dbh_avg = 0.5_wp*(site%cohort%dbh(1) + site%cohort%dbh(2))
   call fuse_2_cohorts(site, 1_ik, 2_ik, cfg%conservation_tol)
   call check_close(site%cohort%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
   call check_close(site%cohort%nplant(1)*site%cohort%agb(1), agb_tot, 1.0e-12_wp,                  &
                    'fused total AGB not conserved')
   call check_close(site%cohort%basal_area(1), pio4*site%cohort%dbh(1)**2, 1.0e-12_wp,              &
                    'basal area inconsistent with re-derived DBH')
   call check(abs(site%cohort%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')

   !=== 2. new_fuse_cohorts reduces count to <= max_cohort, conserving N and AGB. ==========!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
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
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)   ! nplant*leaf_area >> cohort_lai_cap
   call finalize_init(site)
   n0   = total_nplant(site)
   agb0 = total_agb(site)
   call split_cohorts(site, cfg)
   call check(site%cohort%n >= 2_ik, 'split did not create a second cohort')
   call check_close(total_nplant(site), n0,   1.0e-12_wp, 'split broke nplant conservation')
   call check_close(total_agb(site),    agb0, cfg%conservation_tol, 'split broke AGB conservation')

   write(*,'(a)') '   PASS'
end program test_fusion_cohort
