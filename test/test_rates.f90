!----- Rate application (data arrays): growth (+cap), mortality survivorship, recruitment, --!
!----- plus the empirical evaluator's light (overtopping-LAI) growth and temperature gate. --!
program test_rates
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : yr_day
   use meds_config,             only : meds_config_t, build_config
   use meds_demography_interface, only : site_t
   use meds_demography_dynamics, only : growth_step, mortality_step
   use meds_recruitment,        only : recruitment_month
   use meds_test_vital_rates,   only : test_vital_rates
   use meds_init,               only : init_bare_ground, add_cohort, finalize_init
   use meds_test_support,       only : check, check_close, banner
   implicit none

   type(meds_config_t)   :: cfg
   type(site_t)            :: site
   real(wp), allocatable :: g(:), m(:), rec(:,:)
   integer(ik)           :: istep, nday, n0
   real(wp)              :: dexp, full_light

   call banner('rate application (arrays) + empirical evaluator')
   cfg = build_config()

   !=== Growth: a constant per-cohort DBH rate over a year advances DBH by g*t. ===========!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 10.0_wp)      ! climax, dbh_critical=120
   call finalize_init(site)
   allocate(g(site%cohort%n)); g = 2.0_wp
   nday = 365_ik
   do istep = 1_ik, nday
      call growth_step(site%cohort%n, site%cohort%dbh, site%cohort%height, site%cohort%basal_area,       &
                       site%cohort%agb, site%cohort%leaf_area, site%cohort%p_dbh_critical,                    &
                       site%cohort%p_wood_density, g, cfg%dt_years)
   end do
   dexp = 10.0_wp + 2.0_wp * real(nday, wp) / yr_day
   call check_close(site%cohort%dbh(1), dexp, 1.0e-6_wp, 'constant growth did not advance DBH')
   deallocate(g)

   !=== Growth cap: DBH clamps at dbh_critical, never overshoots. =============================!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 119.0_wp)
   call finalize_init(site)
   allocate(g(site%cohort%n)); g = 100.0_wp
   do istep = 1_ik, 30_ik
      call growth_step(site%cohort%n, site%cohort%dbh, site%cohort%height, site%cohort%basal_area,       &
                       site%cohort%agb, site%cohort%leaf_area, site%cohort%p_dbh_critical,                    &
                       site%cohort%p_wood_density, g, cfg%dt_years)
   end do
   call check(site%cohort%dbh(1) <= 120.0_wp + 1.0e-12_wp, 'DBH overshot dbh_critical')
   call check_close(site%cohort%dbh(1), 120.0_wp, 1.0e-9_wp, 'DBH did not clamp at dbh_critical')
   deallocate(g)

   !=== Mortality: 30 per-step applications -> exp survivorship (prod of exp(-m*dt)). =======!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 1.0_wp, 10.0_wp)
   call finalize_init(site)
   allocate(m(site%cohort%n)); m = 0.1_wp
   do istep = 1_ik, 30_ik
      call mortality_step(site%cohort%n, site%cohort%nplant, m, cfg%dt_years, cfg%negligible_nplant)
   end do
   call check_close(site%cohort%nplant(1), exp(-0.1_wp * 30.0_wp / yr_day), 1.0e-9_wp,         &
                    'survivorship /= exp(-m*dt)')
   deallocate(m)

   !=== Recruitment: a positive density spawns one cohort per PFT; zero spawns none. =======!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call finalize_init(site)
   allocate(rec(site%n_pft, site%patch%n)); rec = 0.02_wp   ! >= min_recruit_size 0.01
   call recruitment_month(site, cfg, rec)
   call check(site%cohort%n == site%n_pft, 'recruitment should spawn one cohort per PFT')
   n0 = site%cohort%n
   rec = 0.0_wp
   call recruitment_month(site, cfg, rec)
   call check(site%cohort%n == n0, 'zero recruit density should spawn nothing')
   deallocate(rec)

   !=== Empirical evaluator: overtopping LAI reduces growth; cold shuts recruitment off. ===!
   call init_bare_ground(site, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 30.0_wp)   ! taller (sorted first), over_lai=0
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp,  5.0_wp)   ! shorter, sees overtopping LAI
   call finalize_init(site)
   call test_vital_rates(site, cfg, g, m, rec)
   !----- Taller cohort: full light -> growth equals gr_max*dbh. --------------------------!
   call check_close(g(1), cfg%pft%gr_max(2) * 30.0_wp, 1.0e-12_wp,                          &
                    'top cohort growth should be gr_max*dbh at full light')
   !----- Shorter cohort: overtopping LAI strictly reduces growth below its full-light val. !
   full_light = cfg%pft%gr_max(2) * 5.0_wp
   call check(g(2) < full_light, 'competition did not reduce the understorey growth')
   call check(g(2) > 0.0_wp, 'understorey growth should remain positive')
   !----- Recruitment gate: a cold site_t yields zero recruit density. ----------------------!
   call check(all(rec > 0.0_wp), 'warm site_t should permit recruitment')
   site%patch%min_month_temp(1) = 270.0_wp
   call test_vital_rates(site, cfg, g, m, rec)
   call check(all(rec == 0.0_wp), 'cold site_t should shut recruitment off')

   write(*,'(a)') '   PASS'
end program test_rates
