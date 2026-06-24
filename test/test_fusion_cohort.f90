!----- Cohort fusion/fission: conservation of plant number & basal area, DBH re-derivation. !
program test_fusion_cohort
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : pio4
   use meds_config,          only : meds_config_t, build_config
   use meds_demography_types,           only : site
   use meds_setup,           only : init_bare_ground, add_cohort, finalize_init
   use meds_cohort_dynamics, only : fuse_2_cohorts, new_fuse_cohorts, split_cohorts, max_ccount
   use meds_diagnostics,     only : total_nplant, total_basal_area
   use meds_test_support,    only : check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site)     :: comm
   real(wp)            :: ba_tot, n0, ba0, dbh_avg
   integer(ik)         :: j

   call banner('cohort fusion/fission conservation')
   cfg = build_config()

   !=== 1. fuse_2_cohorts conserves N and BA; DBH is re-derived, not averaged. ============!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.3_wp, 12.0_wp)
   call finalize_init(comm)                          ! sorted: index 1 = taller (dbh 12)
   ba_tot  = comm%coh%nplant(1)*comm%coh%basarea(1) + comm%coh%nplant(2)*comm%coh%basarea(2)
   dbh_avg = 0.5_wp*(comm%coh%dbh(1) + comm%coh%dbh(2))
   call fuse_2_cohorts(comm, 1_ik, 2_ik, cfg%cons_tol)
   call check_close(comm%coh%nplant(1), 0.8_wp, 1.0e-12_wp, 'fused nplant must be summed')
   call check_close(comm%coh%nplant(1)*comm%coh%basarea(1), ba_tot, 1.0e-12_wp,            &
                    'fused total basal area not conserved')
   call check_close(comm%coh%dbh(1), sqrt(comm%coh%basarea(1)/pio4), 1.0e-12_wp,           &
                    'DBH must be re-derived from basal area')
   call check(abs(comm%coh%dbh(1) - dbh_avg) > 1.0e-6_wp, 'DBH must NOT be a plain average')

   !=== 2. new_fuse_cohorts reduces count to <= maxcohort, conserving N and BA. ===========!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   do j = 1_ik, 200_ik
      call add_cohort(comm, cfg, 1_ik, 1_ik, 0.05_wp, 2.0_wp + 0.55_wp*real(j-1_ik, wp))  ! 2..111 cm
   end do
   call finalize_init(comm)
   n0  = total_nplant(comm)
   ba0 = total_basal_area(comm)
   call new_fuse_cohorts(comm, cfg)
   call check(max_ccount(comm) <= abs(cfg%maxcohort), 'fusion did not reach maxcohort')
   call check_close(total_nplant(comm),     n0,  1.0e-10_wp, 'fusion broke nplant conservation')
   call check_close(total_basal_area(comm), ba0, cfg%cons_tol, 'fusion broke BA conservation')

   !=== 3. split_cohorts conserves N and BA and increases the count. ======================!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(comm, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)   ! nplant*basarea ~ 3927 > cap
   call finalize_init(comm)
   n0  = total_nplant(comm)
   ba0 = total_basal_area(comm)
   call split_cohorts(comm, cfg)
   call check(comm%coh%n >= 2_ik, 'split did not create a second cohort')
   call check_close(total_nplant(comm),     n0,  1.0e-12_wp, 'split broke nplant conservation')
   call check_close(total_basal_area(comm), ba0, cfg%cons_tol, 'split broke BA conservation')

   write(*,'(a)') '   PASS'
end program test_fusion_cohort
