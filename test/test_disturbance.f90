!----- Treefall patch disturbance: area conserved, age-0 gap created, tall die / short survive. !
program test_disturbance
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t, DIST_TREEFALL
   use meds_core_state_types, only : site_t
   use meds_core_patch_fusefiss, only : apply_patch_disturbance
   use meds_init,             only : init_bare_ground, add_cohort, finalize_init
   use meds_diagnostic_reduce, only : total_area, total_nplant
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   integer(ik)         :: ig, ip, i0, i1, i, n_gap_cohorts
   real(wp)            :: n_before, h_tall, h_short

   call banner('treefall patch disturbance')
   cfg = build_test_config()

   !----- One patch with a tall canopy cohort and a short understorey cohort. -------------!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall  (height >> threshold)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short (height <  threshold)
   call finalize_init(site)
   h_tall  = site%cohort%height(1)                                ! sorted tallest-first
   h_short = site%cohort%height(site%cohort%n)
   call check(h_tall  >= cfg%disturbance_survive_height, 'tall cohort should exceed threshold')
   call check(h_short <  cfg%disturbance_survive_height, 'short cohort should be below threshold')
   n_before = total_nplant(site)

   call apply_patch_disturbance(site, cfg, 1.0_wp)

   !----- A new age-0 treefall gap patch was opened; area is conserved. --------------------!
   call check(site%patch%n == 2_ik, 'disturbance should add exactly one gap patch')
   call check_close(total_area(site), 1.0_wp, 1.0e-9_wp, 'disturbance broke area conservation')

   ig = 0_ik
   do ip = 1_ik, site%patch%n
      if (site%patch%dist_type(ip) == DIST_TREEFALL) ig = ip
   end do
   call check(ig > 0_ik, 'no treefall gap patch found')
   call check_close(site%patch%age(ig), 0.0_wp, 1.0e-12_wp, 'gap patch should have age 0')
   call check_close(site%patch%shed_water_rate(ig), 0.0_wp, 1.0e-12_wp,                            &
                    'gap patch should have a fresh (zero) shed_water_rate, like age (P4)')

   !----- The gap holds ONLY survivors (height < threshold): the short cohort. -------------!
   i0 = site%patch%cohort_offset(ig) ; i1 = i0 + site%patch%cohort_count(ig) - 1_ik
   n_gap_cohorts = site%patch%cohort_count(ig)
   call check(n_gap_cohorts == 1_ik, 'gap should contain exactly the one survivor cohort')
   do i = i0, i1
      call check(site%cohort%height(i) < cfg%disturbance_survive_height,                        &
                 'tall canopy cohort must not survive in the gap')
   end do

   !----- Net plant number drops (the disturbed canopy fraction is killed). ----------------!
   call check(total_nplant(site) < n_before, 'canopy mortality should reduce plant number')

   write(*,'(a)') '   PASS'
end program test_disturbance
