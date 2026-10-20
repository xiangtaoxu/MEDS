! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_biogeochem_dynamics -- the slow soil-carbon driver seam                             !
! (docs/science/soil_carbon.md §7): a REAL fast loop (constant forcing, over the state-hub !
! reservoirs, same setup as test_fast_loop) feeds a frozen per-patch soil-carbon pool,     !
! then advance_slow_dynamics runs the daily soil_carbon_step against that SAME frozen      !
! pool.                                                                                    !
!   1. THE DOUBLE-COUNT GATE: rh_seam_gap (soil_carbon_step's own rh_today MINUS the fast  !
!      loop's independently accumulated rh_fast_accum) is ~0 for the default EULER solver, !
!      BY CONSTRUCTION -- both read the identical frozen pool and the identical            !
!      day-integrated xi_int. The lignin-tracer residual is ~0, and the freeze number      !
!      lambda is nonzero and far below 1.                                                  !
!   2. THE POOL MOVES, and stays non-negative, through the real driver seam. The exact     !
!      daily identity (after = before + litter_in - rh_today) is soil_carbon_step's own    !
!      audit, and test_soil_biogeochem checks it independently.                            !
!   3. OFF-PATH: soil_carbon_on = .false. leaves the pool untouched (bit-identical,        !
!      matching every other opt-in feature gate).                                          !
!==========================================================================================!
program test_biogeochem_dynamics
   use meds_kinds,               only : wp, ik
   use meds_config,              only : meds_config_t
   use meds_polygon_state_vars,  only : polygon_t
   use meds_init,                only : init_bare_ground, add_cohort, finalize_init
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_constants,     only : n_soil_layer_max
   use meds_soil_column_config, only : build_soil_hydr_params, build_soil_therm_params
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG
   use meds_fast_dynamics,       only : fast_context_t, init_fast_reservoirs, fast_dynamics
   use meds_fast_types,          only : apply_hydraulics_config, apply_canopy_film_config
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_slow_dynamics,       only : advance_slow_dynamics
   use meds_biogeochem_types, only : litter_input_t, soilc_seam_t
   use meds_test_support, only : banner, build_test_config, check, check_close, check_true, test_report
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks are computed on

   integer(ik)            :: n_budget_fail
   integer(ik), parameter :: NSL = 10_ik
   type(meds_config_t)  :: cfg
   type(polygon_t)      :: poly
   type(fast_context_t) :: ctx
   real(wp) :: total0, total1
   type(soilc_seam_t) :: seam   !< per-run seam worsts, accumulated by advance_slow_dynamics
   real(wp) :: we, ww

   call banner('slow soil-carbon biogeochemistry (double-count gate)')

   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg%dt_fast            = 150.0_wp
   cfg%n_fast_per_slow    = 48_ik   ! 48 x 150 s = a 2 h fast window
   cfg%soil_carbon_on     = .true.

   call geometric_layer_bottoms(2.0_wp, NSL, 3.0_wp, test_bottoms)
   call build_soil_hydr_params(NSL, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,            &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, ctx%col_config%soil)
   call build_soil_therm_params(NSL, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%col_config%soil_thermal)
   call apply_hydraulics_config(cfg%pft, ctx%col_config%hydraulics_table)
   call apply_canopy_film_config(cfg%pft, ctx%col_config%leaf_surf_water_max, ctx%col_config%wood_surf_water_max)
   call build_leaf_photo_table(cfg, ctx%col_config%leaf_photo)
   ctx%col_config%integrator = build_integrator_opts(cfg)
   ctx%air_temp = 295.0_wp ; ctx%rad_sw_top = 500.0_wp ; ctx%rad_sw_ground = 75.0_wp
   ctx%soil_temp_init = 295.0_wp

   call init_bare_ground(poly, cfg, 1_ik)

   call use_test_column(poly)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)     ! one climax cohort (turnover -> litter)
   call finalize_init(poly)
   call init_fast_reservoirs(poly, ctx)

   !----- Seed a nonzero starting pool so decomposition (and hence Rh) is actually active. -----!
   poly%patch%soil_carbon(1)%fast_grnd_carbon   = 2.0_wp
   poly%patch%soil_carbon(1)%fast_soil_carbon   = 1.0_wp
   poly%patch%soil_carbon(1)%struct_grnd_carbon = 20.0_wp
   poly%patch%soil_carbon(1)%struct_soil_carbon = 10.0_wp
   poly%patch%soil_carbon(1)%slow_carbon        = 500.0_wp

   total0 = poly%patch%soil_carbon(1)%fast_grnd_carbon + poly%patch%soil_carbon(1)%fast_soil_carbon &
          + poly%patch%soil_carbon(1)%struct_grnd_carbon + poly%patch%soil_carbon(1)%struct_soil_carbon &
          + poly%patch%soil_carbon(1)%microbial_carbon + poly%patch%soil_carbon(1)%slow_carbon      &
          + poly%patch%soil_carbon(1)%passive_carbon

   !=== 1. Run one real slow step: fast loop (accumulates xi_int/rh_fast_accum against the FROZEN !
   !    pool seeded above) THEN advance_slow_dynamics (vegetation turnover -> litter, then the       !
   !    daily soil_carbon_step consuming that litter + xi_int). ====================================!
   call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
   call check(n_budget_fail == 0_ik, 'fast loop whole-column budgets closed (n_fail == 0)')

   block
      real(wp) :: rh_fast_accum_ref
      rh_fast_accum_ref = poly%patch%xi_accum(1)%rh_fast_accum
      call check(rh_fast_accum_ref > 0.0_wp, 'fast loop accumulated positive matrix Rh over the day')

      call advance_slow_dynamics(poly, cfg, seam=seam)

      call check_true('rh_seam_gap ~ 0 (the double-count gate)',                                   &
                      seam%worst_rh_gap < 1.0e-9_wp * max(rh_fast_accum_ref, 1.0_wp), seam%worst_rh_gap)
      !----- The other two seam numbers the same call accumulates. Lambda is the freeze number:     !
      !      the fraction of a pool one slow step withdraws. It must be nonzero (the step DID        !
      !      decompose) and far below 1 (or freezing the pool across the step is unsound).  ---------!
      call check_true('lignin passive-tracer residual ~ 0', seam%worst_lignin < 1.0e-12_wp,        &
                      seam%worst_lignin)
      call check_true('lambda is nonzero -- the step actually decomposed', seam%worst_lambda > 0.0_wp, &
                      seam%worst_lambda)
      call check_true('lambda << 1 -- the frozen-pool approximation is sound at a daily step',      &
                      seam%worst_lambda < 0.1_wp, seam%worst_lambda)
   end block

   total1 = poly%patch%soil_carbon(1)%fast_grnd_carbon + poly%patch%soil_carbon(1)%fast_soil_carbon &
          + poly%patch%soil_carbon(1)%struct_grnd_carbon + poly%patch%soil_carbon(1)%struct_soil_carbon &
          + poly%patch%soil_carbon(1)%microbial_carbon + poly%patch%soil_carbon(1)%slow_carbon      &
          + poly%patch%soil_carbon(1)%passive_carbon

   !----- Direction is not asserted (this cohort's turnover litter can outweigh or be outweighed by !
   !      Rh depending on the seeded pool -- test_soil_biogeochem's INDEPENDENT mass check (Rh        !
   !      recomputed from er*xi_int*K*X0, not from the closure definition) already guarantees          !
   !      dC_pool == litter_in - rh_today to machine precision); this just checks the driver seam      !
   !      actually moved the pool and stayed physical. --------------------------------------------!
   call check(abs(total1 - total0) > 1.0e-6_wp, 'pool total changed (the daily step actually ran)')
   call check(total1 > 0.0_wp, 'pool total stayed physical (non-negative)')
   write(*,'(a,es12.5,a,es12.5,a)') '   (pool total: ', total0, ' -> ', total1, ' kgC/m2)'
   write(*,'(a,es12.5,a,es12.5,a)') '   (rh_seam_gap: ', seam%worst_rh_gap,                       &
                                    '   lambda: ', seam%worst_lambda, ')'

   !=== 2. OFF-PATH: soil_carbon_on = .false. must leave the pool exactly untouched. ==============!
   block
      type(polygon_t) :: poly2
      real(wp)     :: totalA, totalB
      cfg%soil_carbon_on = .false.
      call init_bare_ground(poly2, cfg, 1_ik)
      call use_test_column(poly2)
      call add_cohort(poly2, cfg, 1_ik, 3_ik, 0.5_wp, 20.0_wp)
      call finalize_init(poly2)
      call init_fast_reservoirs(poly2, ctx)
      poly2%patch%soil_carbon(1)%slow_carbon = 500.0_wp
      totalA = poly2%patch%soil_carbon(1)%slow_carbon
      call fast_dynamics(poly2, ctx, cfg)
      call advance_slow_dynamics(poly2, cfg)
      totalB = poly2%patch%soil_carbon(1)%fast_grnd_carbon + poly2%patch%soil_carbon(1)%fast_soil_carbon &
             + poly2%patch%soil_carbon(1)%struct_grnd_carbon + poly2%patch%soil_carbon(1)%struct_soil_carbon &
             + poly2%patch%soil_carbon(1)%microbial_carbon + poly2%patch%soil_carbon(1)%slow_carbon  &
             + poly2%patch%soil_carbon(1)%passive_carbon
      call check_close(totalB, totalA, 1.0e-12_wp, 'soil_carbon_on=.false.: pool exactly untouched')
   end block

   write(*,'(a)') '   PASS'

   call test_report('test_biogeochem_dynamics')

contains

   !----- The test's hand-built soil column is the site's, which the fast loop runs on, starting  !
   !      at 0.30 in every layer. ------------------------------------------------------------!
   subroutine use_test_column(p)
      type(polygon_t), intent(inout) :: p
      p%site%surface(1)%soil         = ctx%col_config%soil
      p%site%surface(1)%soil_thermal = ctx%col_config%soil_thermal
      p%site%surface(1)%theta_start  = 0.30_wp
   end subroutine use_test_column


end program test_biogeochem_dynamics
