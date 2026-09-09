!----- Treefall patch disturbance: area conserved, age-0 gap created, tall die / short survive. !
program test_disturbance
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t, DIST_TREEFALL
   use meds_site_state_types, only : site_t
   use meds_demography_patch_fusefiss, only : apply_patch_disturbance
   use meds_init,             only : init_bare_ground, add_cohort, finalize_init
   use meds_diagnostic_reduce, only : total_area, total_nplant
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   integer(ik)         :: ig, ip, i0, i1, i, n_gap_cohorts
   real(wp)            :: n_before, h_tall, h_short, film_before, film_after

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
   !----- REVIEW 2026-09 (item 1B #1): the survivor's film [kg/m2 of ITS patch] must dilute into the !
   !      gap like its density; the site total of the SURVIVOR's film is conserved across the split. -!
   site%cohort%leaf_surf_water(site%cohort%n) = 0.25_wp
   film_before = site%patch%area(1) * 0.25_wp

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
   film_after = 0.0_wp
   do ip = 1_ik, site%patch%n
      i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
      do i = i0, i1
         if (site%cohort%height(i) < cfg%disturbance_survive_height)                                &
            film_after = film_after + site%patch%area(ip) * site%cohort%leaf_surf_water(i)
      end do
   end do
   call check_close(film_after, film_before, 1.0e-12_wp, 'disturbance broke the survivor film-water conservation')
   i0 = site%patch%cohort_offset(ig) ; i1 = i0 + site%patch%cohort_count(ig) - 1_ik
   do i = i0, i1
      call check(site%cohort%height(i) < cfg%disturbance_survive_height,                        &
                 'tall canopy cohort must not survive in the gap')
   end do

   !----- Net plant number drops (the disturbed canopy fraction is killed). ----------------!
   call check(total_nplant(site) < n_before, 'canopy mortality should reduce plant number')

   !----- TWO donors of unequal area and unequal film: the gap's survivor copies must DILUTE their   !
   !      films by frac*area(d)/new_area exactly as their densities are (with one donor that factor  !
   !      is 1 and a verbatim copy passes by accident). Site total before = 0.6*0.25 = 0.15. --------!
   call init_bare_ground(site, cfg, 2_ik)
   site%patch%area(1) = 0.6_wp ; site%patch%area(2) = 0.4_wp
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall, dies
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short survivor, WET
   call add_cohort(site, cfg, 2_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall, dies
   call add_cohort(site, cfg, 2_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short survivor, dry
   call finalize_init(site)
   film_before = 0.0_wp
   do ip = 1_ik, site%patch%n
      i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
      do i = i0, i1
         if (site%cohort%height(i) < cfg%disturbance_survive_height .and. ip == 1_ik) then
            site%cohort%leaf_surf_water(i) = 0.25_wp
            film_before = film_before + site%patch%area(ip) * 0.25_wp
         end if
      end do
   end do
   call apply_patch_disturbance(site, cfg, 0.5_wp)
   film_after = 0.0_wp
   do ip = 1_ik, site%patch%n
      i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
      do i = i0, i1
         film_after = film_after + site%patch%area(ip) * site%cohort%leaf_surf_water(i)
      end do
   end do
   call check_close(total_area(site), 1.0_wp, 1.0e-9_wp, 'two-donor disturbance broke area conservation')
   call check_close(film_after, film_before, 1.0e-12_wp,                                            &
                    'two-donor disturbance: survivor films dilute into the gap like their densities')

   write(*,'(a)') '   PASS'
end program test_disturbance
