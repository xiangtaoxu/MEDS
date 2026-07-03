!----- Rate application (data arrays): growth (+cap), mortality survivorship, recruitment, --!
!----- plus the phenomenological evaluator (intrinsic+competition+reproduction growth,       !
!----- growth-dependent mortality, seed-rain+reproduction recruitment). ---------------------!
program test_rates
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : yr_day
   use meds_config,             only : meds_config_t
   use meds_demography_interface, only : site_t
   use meds_demography_dynamics, only : growth_step, mortality_step, apply_recruitment
   use meds_allometry,           only : b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2
   use meds_phenomenological_vital_rates, only : vital_rates
   use meds_init,               only : init_bare_ground, add_cohort, finalize_init
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t)   :: cfg
   type(site_t)            :: site
   real(wp), allocatable :: g(:), m(:), rec(:,:)
   integer(ik)           :: istep, nday, n0, pf
   real(wp)              :: dexp, gi_top, gi_short, expect_top

   call banner('rate application (arrays) + phenomenological evaluator')
   cfg = build_test_config()

   !=== Growth: a constant per-cohort DBH rate over a year advances DBH by g*t. ===========!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 10.0_wp)      ! climax, dbh_critical=120
   call finalize_init(site)
   allocate(g(site%cohort%n)); g = 2.0_wp
   nday = 365_ik
   do istep = 1_ik, nday
      call growth_step(site%cohort%n, site%cohort%dbh, site%cohort%height, site%cohort%basal_area,       &
                       site%cohort%agb, site%cohort%leaf_area, site%cohort%growth_avg,                    &
                       site%cohort%growth_accum, site%cohort%growth_count, site%cohort%growth_hist,       &
                       site%cohort%p_dbh_critical, site%cohort%p_wood_density, site%cohort%p_hgt_max,     &
                       g, cfg%dt_years, 90_ik, 1_ik, b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2)
   end do
   dexp = 10.0_wp + 2.0_wp * real(nday, wp) / yr_day
   call check_close(site%cohort%dbh(1), dexp, 1.0e-6_wp, 'constant growth did not advance DBH')
   deallocate(g)

   !=== Growth cap: DBH clamps at dbh_critical, never overshoots. =============================!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, cfg%pft%dbh_critical(3_ik) - 1.0_wp)  ! just below the cap
   call finalize_init(site)
   allocate(g(site%cohort%n)); g = 100.0_wp
   do istep = 1_ik, 30_ik
      call growth_step(site%cohort%n, site%cohort%dbh, site%cohort%height, site%cohort%basal_area,       &
                       site%cohort%agb, site%cohort%leaf_area, site%cohort%growth_avg,                    &
                       site%cohort%growth_accum, site%cohort%growth_count, site%cohort%growth_hist,       &
                       site%cohort%p_dbh_critical, site%cohort%p_wood_density, site%cohort%p_hgt_max,     &
                       g, cfg%dt_years, 90_ik, 1_ik, b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2)
   end do
   call check(site%cohort%dbh(1) <= cfg%pft%dbh_critical(3_ik) + 1.0e-12_wp, 'DBH overshot dbh_critical')
   call check_close(site%cohort%dbh(1), cfg%pft%dbh_critical(3_ik), 1.0e-9_wp, 'DBH did not clamp at dbh_critical')
   deallocate(g)

   !=== Mortality: 30 per-step applications -> exp survivorship (prod of exp(-m*dt)). =======!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 1.0_wp, 10.0_wp)
   call finalize_init(site)
   allocate(m(site%cohort%n)); m = 0.1_wp
   do istep = 1_ik, 30_ik
      call mortality_step(site%cohort%n, site%cohort%nplant, m, cfg%dt_years, cfg%negligible_nplant)
   end do
   call check_close(site%cohort%nplant(1), exp(-0.1_wp * 30.0_wp / yr_day), 1.0e-9_wp,         &
                    'survivorship /= exp(-m*dt)')
   deallocate(m)

   !=== Recruitment application: a positive density spawns one cohort per PFT; zero none. ===!
   call init_bare_ground(site, cfg, 1_ik)
   call finalize_init(site)
   allocate(rec(site%n_pft, site%patch%n)); rec = 0.2_wp    ! per-yr rate; /12 per month >= min_recruit_size 0.01
   call apply_recruitment(site, cfg, rec)
   call check(site%cohort%n == site%n_pft, 'recruitment should spawn one cohort per PFT')
   n0 = site%cohort%n
   rec = 0.0_wp
   call apply_recruitment(site, cfg, rec)
   call check(site%cohort%n == n0, 'zero recruit density should spawn nothing')
   deallocate(rec)

   !=== Phenomenological evaluator: intrinsic growth, competition suppression, growth-driven =!
   !    mortality, and a positive recruitment rate. Use PFT 1 (pioneer) so mortality is        !
   !    non-zero (mid/climax PFTs clamp to zero baseline with the default coefficients).        !
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 30.0_wp)   ! taller (sorted first), over_lai=0, mature
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp,  5.0_wp)   ! shorter, sees overtopping LAI
   call finalize_init(site)
   call vital_rates(site, cfg, g, m, rec)                    ! growth + mortality + recruitment
   pf = 1_ik
   !----- Top cohort at full light, above the reproduction height: growth = intrinsic *       !
   !      (1 - reproduction_investment_fraction). ----------------------------------------------!
   gi_top = cfg%pft%growth_dbh_max(pf) * exp(cfg%pft%growth_dbh_slope(pf)                      &
            * max(0.0_wp, log(cfg%pft%growth_dbh_cap(pf)) - log(30.0_wp)))
   expect_top = gi_top * (1.0_wp - cfg%pft%reproduction_investment_fraction(pf))
   call check_close(g(1), expect_top, 1.0e-9_wp, 'top cohort growth /= intrinsic*(1-repro) at full light')
   !----- Shorter cohort: overtopping LAI strictly reduces growth below its full-light value. !
   gi_short = cfg%pft%growth_dbh_max(pf) * exp(cfg%pft%growth_dbh_slope(pf)                    &
              * max(0.0_wp, log(cfg%pft%growth_dbh_cap(pf)) - log(5.0_wp)))
   call check(g(2) < gi_short, 'competition did not reduce the understorey growth')
   call check(g(2) > 0.0_wp,   'understorey growth should remain positive')
   !----- Growth-dependent mortality: the suppressed (low-growth) cohort dies faster. --------!
   call check(m(2) > m(1), 'low-growth cohort should have higher mortality')
   call check(m(1) > 0.0_wp, 'pioneer mortality should be positive')
   !----- Recruitment rate: every included PFT has a positive (seed-rain) recruitment rate. --!
   call check(all(rec > 0.0_wp), 'all included PFTs should have a positive recruitment rate')
   deallocate(g, m, rec)

   write(*,'(a)') '   PASS'
end program test_rates
