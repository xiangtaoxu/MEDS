!==========================================================================================!
! test_slow_diag_units -- the units contract of the SLOW per-patch diagnostic rows (#239).      !
!                                                                                          !
! The patch diagnostic block is (value, weight) and the reader returns value/weight. The weight  !
! is Sum(dt) in SECONDS -- accumulate_patch_diag adds cfg%dt_fast every fast step -- so EVERY      !
! row, fast or slow, must contribute (its rate, in the units the registry declares) x (dt in       !
! SECONDS). The fast rows do that with flux*dt_fast and always have. The slow rows, added later,   !
! contributed a bare per-step AMOUNT against a comment that described the rule correctly, and so   !
! emitted a per-SECOND rate under a per-YEAR label -- every one of them read a factor              !
! yr_sec = 3.1557e7 too small, for six shipped variables.                                          !
!                                                                                          !
! Nothing else in the suite could see it: a diagnostic is downstream of every conservation ledger, !
! and the slow ledger's own litter term reads patch%litter_in directly rather than this slot.      !
!                                                                                          !
! So this test asserts the one identity the per-year label MEANS, against an oracle outside the    !
! diagnostic path:                                                                                 !
!                                                                                          !
!     reported_rate [x/yr]  x  elapsed [yr]  ==  the amount that actually flowed [x]                !
!                                                                                          !
! -- for litter, site%patch%litter_in, the very accumulator biogeochemistry consumes; for          !
! recruitment, the carry-forward recruit pool the same step credits; for disturbance, the area     !
! fraction its own hazard implies. No physics number is asserted, so the test cannot go stale with !
! the parameters, and a scale error of any size in any of the five slots fails it.                  !
!                                                                                          !
! The block's window is ONE SLOW STEP: fast_dynamics resets it on entry, and the output layer's     !
! temporal integrator folds the per-step windows into the monthly and annual streams. So the        !
! identity is asserted per slow step, which is the level the contribution is written at.            !
!                                                                                          !
! THE FIXTURE MIRRORS THE DRIVER, and says so: check 1 asserts the block's weight after one real   !
! fast window equals cfg%dt_slow, because the identity above is only meaningful if the fast loop   !
! really did span the slow step. A fixture whose fast window is shorter than its slow step would   !
! make every later check pass for the wrong reason.                                                 !
!==========================================================================================!
program test_slow_diag_units
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : yr_sec
   use meds_config,              only : meds_config_t
   use meds_site_state_types,    only : site_t
   use meds_site_diag_types,     only : patch_diag_alloc, patch_diag_value,                       &
                                        PD_LITTER_LEAF, PD_LITTER_FINEROOT, PD_LITTER_STRUCT,     &
                                        PD_RECRUIT_NPLANT, PD_DISTURB_AREA
   use meds_init,                only : init_bare_ground, add_cohort, finalize_init
   use meds_column_params,       only : build_soil_hydr_params, build_soil_therm_params
   use meds_hydr_lib,            only : SOIL_RETENTION_VG
   use meds_fast_dynamics,       only : fast_context_t, init_fast_reservoirs, fast_dynamics
   use meds_fast_types,          only : apply_hydraulics_config
   use meds_fast_config,         only : build_leaf_photo_table, build_integrator_opts
   use meds_slow_dynamics,       only : advance_slow_dynamics
   use meds_test_support,        only : banner, build_test_config, check, check_close, check_true, &
                                        test_report
   implicit none

   integer(ik), parameter :: nsl = 10_ik
   type(meds_config_t)  :: cfg
   type(site_t)         :: site
   type(fast_context_t) :: ctx
   real(wp)             :: x(8), pool0, pool1, elapsed_yr, frac_expected
   integer(ik)          :: np

   call banner('slow per-patch diagnostic units (#239)')

   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   !----- ONE full slow step of fast sub-steps: 96 x 900 s = 86400 s = cfg%dt_slow. This is the    !
   !      property check 1 asserts; it is set here rather than inherited so the mirror is explicit.!
   cfg%dt_fast            = 900.0_wp
   cfg%n_fast_per_slow    = 96_ik
   cfg%soil_carbon_on     = .true.

   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,          &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ctx%col_config%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%col_config%soil_thermal)
   call apply_hydraulics_config(cfg%hydraulics, cfg%pft, ctx%col_config%hydraulics_table)
   call build_leaf_photo_table(cfg, ctx%col_config%leaf_photo)
   ctx%col_config%integrator = build_integrator_opts(cfg)
   ctx%air_temp = 295.0_wp ; ctx%rad_sw_top = 500.0_wp ; ctx%rad_sw_ground = 75.0_wp
   ctx%theta_init = 0.30_wp ; ctx%soil_temp_init = 295.0_wp

   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)     ! one climax cohort (turnover -> litter)
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)

   !----- Turn the patch diagnostic block on by hand. In production meds_output_registry does this !
   !      from the enabled-variable set; a test has no registry, and the block is a no-op at every !
   !      entry point until it is active. -----------------------------------------------------!
   call patch_diag_alloc(site%patch%diag, max(site%patch%cap, 1_ik), .true.)

   !=== 1. THE FIXTURE MIRRORS THE DRIVER: one fast window weighs exactly one slow step. =========!
   call fast_dynamics(site, ctx, cfg)
   call check_close(site%patch%diag%w(1), cfg%dt_slow, 1.0e-9_wp,                                  &
                    'fast window accumulates exactly dt_slow of weight (the fixture mirrors the driver)')
   elapsed_yr = site%patch%diag%w(1) / yr_sec

   !=== 2. LITTER: reported rate x elapsed == the amount biogeochemistry consumed. ===============!
   pool0 = sum(site%patch%recruit_pool(:, 1))
   call advance_slow_dynamics(site, cfg, .false., .false.)
   pool1 = sum(site%patch%recruit_pool(:, 1))

   call check(site%patch%litter_in(1)%labile_grnd > 0.0_wp, 'the step actually shed leaf litter')

   call patch_diag_value(site%patch%diag, PD_LITTER_LEAF, x, np)
   call check(np >= 1_ik, 'the patch diagnostic block reports at least one slot')
   call check_close(x(1) * elapsed_yr, site%patch%litter_in(1)%labile_grnd, 1.0e-10_wp,            &
                    'PD_LITTER_LEAF rate x elapsed == the leaf litter amount')
   write(*,'(a,es12.5,a)') '   (leaf litterfall: ', x(1), ' kgC/m2/yr)'

   call patch_diag_value(site%patch%diag, PD_LITTER_FINEROOT, x, np)
   call check_close(x(1) * elapsed_yr, site%patch%litter_in(1)%labile_soil, 1.0e-10_wp,            &
                    'PD_LITTER_FINEROOT rate x elapsed == the fine-root litter amount')

   call patch_diag_value(site%patch%diag, PD_LITTER_STRUCT, x, np)
   call check_close(x(1) * elapsed_yr,                                                             &
                    site%patch%litter_in(1)%struct_grnd + site%patch%litter_in(1)%struct_soil,     &
                    1.0e-10_wp, 'PD_LITTER_STRUCT rate x elapsed == the structural litter amount')

   !=== 3. RECRUITMENT: the same identity against the pool the SAME step credits. The pool is     !
   !    consumed only on a month boundary, and this step is not one, so its change IS the amount.  !
   call patch_diag_value(site%patch%diag, PD_RECRUIT_NPLANT, x, np)
   call check(pool1 - pool0 > 0.0_wp, 'the step actually credited the recruit pool')
   call check_close(x(1) * elapsed_yr, pool1 - pool0, 1.0e-10_wp,                                  &
                    'PD_RECRUIT_NPLANT rate x elapsed == the recruit-pool credit')

   !=== 4. DISTURBANCE: fires on a year boundary, in a different file with its own expression, so  !
   !    it gets its own check. The amount is the area fraction its hazard removes over the one-    !
   !    year PATCH_DYNAMICS_INTERVAL it is called with. -------------------------------------------!
   block
      call fast_dynamics(site, ctx, cfg)                         ! a fresh window (the block resets)
      elapsed_yr = site%patch%diag%w(1) / yr_sec
      call check_close(site%patch%diag%w(1), cfg%dt_slow, 1.0e-9_wp,                               &
                       'the second window also weighs exactly one slow step')
      call advance_slow_dynamics(site, cfg, .true., .true.)      ! new month AND new year
      frac_expected = 1.0_wp - exp(-cfg%patch_disturbance_rate * 1.0_wp)   ! PATCH_DYNAMICS_INTERVAL
      call check(frac_expected > 0.0_wp, 'the test config actually disturbs')
      call patch_diag_value(site%patch%diag, PD_DISTURB_AREA, x, np)
      call check_close(x(1) * elapsed_yr, frac_expected, 1.0e-9_wp,                                &
                       'PD_DISTURB_AREA rate x elapsed == the area fraction disturbed')
      write(*,'(a,es12.5,a,es12.5,a)') '   (disturbance: reported ', x(1),                         &
                                       ' 1/yr  vs configured hazard ', cfg%patch_disturbance_rate, ' 1/yr)'
   end block

   call test_report('test_slow_diag_units')

end program test_slow_diag_units
