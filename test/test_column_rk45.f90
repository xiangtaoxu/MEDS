! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_column_rk45 -- the INTEG_RK45 wiring of the ED2-faithful adaptive Cash-Karp integrator    !
! into column_fast_step (docs/science/numerical_scheme.md §3, `rk45`).                           !
!   A. GPP PARITY: build_column_frozen (shared with ARK) runs the SAME frozen pre-pass under     !
!      both schemes, so gpp_coh from the RK45 path is bit-identical to ARK's on step 1.          !
!   B. PHYSICAL: a 24-step dry-window march under INTEG_RK45 keeps the CAS / soil / leaf / plant    !
!      water-mass state finite, bounded, and sub-saturated.                                        !
!   C. WHOLE-COLUMN CONSERVATION LEDGER: the water AND energy budgets close over a 24 h dry          !
!      diurnal march to machine precision. RK45 carries every store in its tableau, so it is     !
!      gated at a tight tolerance directly.                                                      !
!   K. CAMPBELL SEED: on a Campbell soil the pre-pass seeds a new cohort's predawn potential     !
!      from the top layer through the Campbell curve.                                            !
!==========================================================================================!
program test_column_rk45
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o, grav_head, n_soil_layer_max, tissue_c_leaf, tissue_c_wood_dry
   use meds_config,              only : meds_config_t
   use meds_fast_opts,           only : INTEG_ARK, INTEG_RK45
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_temp_to_enthalpy, temp_to_internal_energy, cas_enthalpy_to_temp,  &
                                           sat_specific_humidity, tissue_temp_to_energy,                         &
                                           tissue_dry_heat_capacity, internal_energy_liquid
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_fast_types, only : patch_biophys_t, alloc_patch_biophys
   use meds_column_view,       only : column_cohort_init
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, soil_theta_to_head
   use meds_biophysics_opts, only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use meds_canopy_types, only : set_aero_env_atm
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_column_config, only : soil_params_t, build_soil_hydr_params, build_soil_therm_params
   use meds_diag_vars,       only : N_CO_FAST_ROWS, CO_PSI_PREDAWN
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config, apply_canopy_film_config
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_fast_step,          only : column_fast_step
   use meds_fast_types,         only : fast_workspace_t
   use meds_plant_hydr_lib,      only : tissue_water_to_psi, tissue_psi_to_water
   use meds_test_support, only : build_test_config, check_true, test_report
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks were written on

   integer(ik), parameter :: N = 1_ik, NSL = 10_ik
   real(wp),    parameter :: DT_FAST = 150.0_wp, LAT = 40.0_wp, T0 = 288.0_wp, THETA0 = 0.30_wp
   !----- Durations and times of day are written in hours and turned into steps here, so a change of  !
   !      DT_FAST keeps what each test covers: set_diurnal_forcing puts step i at (i - 1/2)*DT_FAST    !
   !      after midnight, and the sun follows the clock. -------------------------------------------!
   integer(ik), parameter :: STEPS_PER_HOUR = nint(3600.0_wp / DT_FAST, ik)
   !----- seed moisture, a VARIABLE so the saturated test can drive the column to theta_sat.   !
   real(wp) :: theta_seed = THETA0
   real(wp), parameter :: PSI_FIXTURE = -0.1_wp   !< [MPa] the fixture tissues' starting potential, a mild tension
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: col_config
   type(column_cohort_t)  :: col_cohort
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: biophys
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budget
   type(fast_workspace_t) :: ws          !< the fast step's storage, reused step to step
   type(meds_time_t)      :: sim_date
   real(wp)    :: gpp_split(N), gpp_rk45(N), gpp_coh(N), tcas, qsat, worst_super
   real(wp)    :: psi_leaf_diag, psi_leaf_probe(N)
   integer(ik) :: is, k
   logical     :: physical

   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- column setup (mirrors test_column_ark). ---------------------------------------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   !----- On-allometry fixture, a 20 cm tree at 224 stems/ha (LAI ~3): born through the canonical !
   !      path (init_cohort + set_cohort_size) and gathered by the production gather, so every    !
   !      field -- height, leaf area, the wood carbon the wood heat capacity is sized from -- is  !
   !      mutually consistent. -------------------------------------------------------------------!
   call column_cohort_init(col_cohort, cfg%pft, [1_ik], [20.0_wp], [0.0224_wp])
   call geometric_layer_bottoms(2.0_wp, NSL, 3.0_wp, test_bottoms)
   call build_soil_hydr_params(NSL, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, col_config%soil)
   call build_soil_therm_params(NSL, 3.0_wp, 0.15_wp, 2.0e6_wp, col_config%soil_thermal)
   call apply_hydraulics_config(cfg%pft, col_config%hydraulics_table)
   call apply_canopy_film_config(cfg%pft, col_config%leaf_surf_water_max, col_config%wood_surf_water_max)
   call build_leaf_photo_table(cfg, col_config%leaf_photo)
   col_config%integrator = build_integrator_opts(cfg)
   call alloc_aero_out(aero, N)
   allocate(forc%abs_sw(N), forc%abs_lw(N), forc%abs_par(N), forc%abs_sw_wood(N), forc%abs_lw_wood(N))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   print '(a)', '[test_column_rk45]'

   !=== A. GPP parity on step 1 (identical initial state; RK45 pre-pass == ARK pre-pass). ======!
   call set_noon_forcing()
   call reset_state()
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   gpp_split = gpp_coh
   call reset_state()
   cfg%time_integrator = INTEG_RK45
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   gpp_rk45 = gpp_coh
   call check_true('RK45 pre-pass gpp bit-identical to the ARK pre-pass (build_column_frozen)',                       &
           abs(gpp_rk45(1) - gpp_split(1)) < 1.0e-12_wp, abs(gpp_rk45(1) - gpp_split(1)))
   call check_true('RK45 midday gpp > 0', gpp_rk45(1) > 0.0_wp, gpp_rk45(1))

   !=== B. A dry-window march under INTEG_RK45 stays physical + bounded + sub-saturated. ========!
   call reset_state()
   cfg%time_integrator = INTEG_RK45
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   physical = .true. ; worst_super = -1.0_wp
   do is = 1_ik, 6_ik * STEPS_PER_HOUR
      call set_diurnal_forcing(is)
      call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      tcas = biophys%cas%can_temp
      qsat = sat_specific_humidity(tcas, aenv%press)
      worst_super = max(worst_super, biophys%cas%can_shv - qsat)
      physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. biophys%cas%can_shv > 0.0_wp
      physical = physical .and. biophys%leaf_temp(1) > 250.0_wp .and. biophys%leaf_temp(1) < 340.0_wp
      do k = 1_ik, NSL
         physical = physical .and. biophys%soil_water%theta(k) > 0.0_wp .and. biophys%soil_water%theta(k) < 0.6_wp
         physical = physical .and. biophys%soil_energy%soil_temp(k) > 260.0_wp .and. biophys%soil_energy%soil_temp(k) < 340.0_wp
      end do
      physical = physical .and. biophys%leaf_water_mass(1) > 0.0_wp .and. biophys%wood_water_mass(1) > 0.0_wp
      psi_leaf_diag = tissue_water_to_psi(biophys%leaf_water_mass(1), col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1))
      physical = physical .and. psi_leaf_diag < 0.5_wp .and. psi_leaf_diag > -12.0_wp
   end do
   call check_true('INTEG_RK45 dry-window march stays physical + bounded (00:00-06:00)', physical, biophys%cas%can_temp)
   call check_true('INTEG_RK45 CAS stays sub-saturated', worst_super <= 1.0e-4_wp, worst_super)

   !=== C. WHOLE-COLUMN CONSERVATION LEDGER: the water AND energy budgets close over a 24 h dry   !
   !       diurnal march at the ledger's tight tolerance, with the energy residual under 1 J. ====!
   call test_rk45_budgets()

   !=== D. PRECIP>0: a WET diurnal march runs and still closes the whole-column water AND energy   !
   !       budgets to the SAME tight tolerance. It exercises the boundary water-enthalpy advection !
   !       (infiltration and drainage enthalpy into the boundary soil layers, e_infil/e_drain)     !
   !       under a REAL, substantial rainfall rate, not just incidental background drainage. ======!
   call test_rk45_budgets_wet()

   !=== E. CANOPY-SURFACE WATER (opt-in): a diurnal march with a morning rain pulse actually gets        !
   !       intercepted; whole_water closes exactly, whole_energy stays BOUNDED (the film has no thermal  !
   !       state of its own, so its sensible heat at the tissue temperature offset is not counted; see   !
   !       docs/ROADMAP.md; mirrors test_column_dynamics.f90's RUN 6 and test_column_ark.f90's           !
   !       test_ark_canopy_water). Proves the WIRING end to end for RK45's integrated surf_water ODE     !
   !       (not the wetted-fraction algebra itself, unit-tested in test_surface_energy.f90). ============!
   call test_rk45_canopy_water()

   !=== F. LEAF/ROOT-TURNOVER SHED WATER: a constant shed_water_rate                                     !
   !       (distinct from rainfall) wets the soil and both whole_water/whole_energy still close. ========!
   call test_rk45_shed_water()
   call test_rk45_saturated()

   !=== F2. DRYDOWN to the RESIDUAL floor: the OTHER edge of the constitutive domain. The water network !
   !       solves each dry layer in effective saturation, which cannot go below zero, so nothing      !
   !       corrects theta after the step. This is the only scenario in the file that pushes a layer    !
   !       down toward that floor. =======================================================================!
   call test_rk45_drydown()

   !=== G. TINY NEWLY-RECRUITED COHORT: LAI/WAI at a just-recruited cohort's real magnitude (a 30-yr         !
   !       Ithaca run's first recruitment event) guards two RK45-only failure modes a mature-canopy          !
   !       cohort (LAI=3/WAI=0.5, tests A-F above) never exercises: (1) leaf/wood_temp overflowing           !
   !       through a near-zero veg_energy_balance denominator, which cannot happen because the forcing,      !
   !       the coupling and the heat capacity all scale with the tissue's area                               !
   !       (docs/science/vegetation_energy_dynamics.md eqs. (1), (8)), and (2) the canopy air or the soil    !
   !       leaving its physical domain within a single explicit stage for a near-bare patch, which the       !
   !       RK45 stage clamp on the canopy air and the frozen ground evaporation keep evaluable               !
   !       (meds_fast_rk45.f90). A mature canopy, whose coupling and heat capacity are large and whose       !
   !       surface<->soil coupling is far from stiff, trips neither -- this test is the ONLY one in the      !
   !       suite with a cohort that small.                                                                   !
   call test_rk45_tiny_cohort()

   !=== H. DENSE COLD CANOPY: a very dense canopy (LAI=8) under a hard cold snap makes the diagnostic    !
   !       leaf<->CAS coupling stiff enough that the FULLY-EXPLICIT RK45 surface can oscillate and rail  !
   !       the CAS/soil to the clamp bounds within one DT_FAST. The dispatcher's hybrid rescue detects   !
   !       the railed commit and REDOES the step on ARK (docs/science/numerical_scheme.md §3, `rk45`),   !
   !       keeping the state physical. Asserts that CAS/soil/leaf stay in [250,340] K with GPP finite    !
   !       throughout; the rescue count is reported, not asserted (see the routine). ====================!
   call test_rk45_dense_cold_canopy()
   call test_rk45_bedrock_and_aquifer()
   call test_rk45_prognostic_wood()

   !=== I. The dispatcher must report psi_leaf_coh on the RK45 success path too. Filled only after    !
   !       the RK45 block, whose success branch returns early, the daily-max accumulator would read   !
   !       whatever the caller's per-thread buffer held (another patch, or uninitialised memory) and  !
   !       compute beta_stomata from it. Seed the buffer with an impossible POSITIVE potential and    !
   !       assert the step overwrote it with a physical (<= 0) value without an ARK rescue having     !
   !       run. ======================================================================================!
   call test_rk45_reports_psi_leaf()

   !=== J. On moderately DRY soil the wood<->soil interface must cancel                                 !
   !       to machine precision: the soil loses exactly the uptake the wood gains. The existing         !
   !       dry-down sits near theta_res and the moist tests in wet soil; this one sits in between, and  !
   !       asks for closure at the ledger's own round-off rather than at its 1e-4 relative tolerance. ==!
   call test_rk45_dry_uptake_seam()

   !=== K. A new cohort's predawn potential is seeded from the top soil layer's (the leaf kernel's   !
   !       stomatal limb reads it until the cohort has a day of its own). On a Campbell soil the      !
   !       seed must come from the Campbell curve: the column's one parameter pair is (psi_sat, b)    !
   !       there, and reading it as van Genuchten's would give a wrong potential, or none. ===========!
   call test_campbell_predawn_seed()

   call test_report('test_column_rk45')

contains

   !----- Prognostic WOOD on RK45, via the SAME exact tissue-store relaxation the ARK uses          !
   !      (docs/science/vegetation_energy_dynamics.md eqs. (3)-(4); test_column_ark's               !
   !      test_wood_prognostic carries the rationale). Both adaptive schemes run one wood model. ---!
   subroutine test_rk45_prognostic_wood()
      integer(ik) :: istep
      real(wp)    :: dmax_lag
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      dmax_lag = 0.0_wp
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         dmax_lag = max(dmax_lag, abs(biophys%wood_temp(1) - biophys%cas%can_temp))
      end do
      call check_true('RK45 PROG-WOOD: whole_energy closes', budget%whole_energy%n_fail == 0_ik,                      &
              real(budget%whole_energy%n_fail, wp))
      call check_true('RK45 PROG-WOOD: whole_water closes', budget%whole_water%n_fail == 0_ik,                        &
              real(budget%whole_water%n_fail, wp))
      !----- RK45 is the path where a borrowed advective face would hide: a flux taken from one solve  !
      !      with state committed from the scheme's own stages. The per-layer residual is measured at  !
      !      RK45's own commit, on the b-weighted state (docs/science/numerical_scheme.md §5b). -------!
      call check_true('RK45 PROG-WOOD: per-layer face closure holds', budget%soil_face_mass%worst < 1.0e-8_wp, &
              budget%soil_face_mass%worst)
      call check_true('RK45 PROG-WOOD: the face check actually ran', budget%soil_face_mass%n_check > 0_ik, &
              real(budget%soil_face_mass%n_check, wp))
      call check_true('RK45 PROG-WOOD: wood temperature lags the CAS', dmax_lag > 1.0e-3_wp, dmax_lag)
      call check_true('RK45 PROG-WOOD: wood temperature physical',                                                    &
              biophys%wood_temp(1) > 200.0_wp .and. biophys%wood_temp(1) < 350.0_wp, biophys%wood_temp(1))
      print '(a,i0,a,i0)', '   RK45 PROG-WOOD last dt_fast: substeps = ', budget%integ_nsteps,        &
            ' , rescues = ', budget%rk45_rescue
   end subroutine test_rk45_prognostic_wood

   !----- All three bottom BCs are pure boundary conditions with no prognostic state behind them:  !
   !      bedrock seals the face; the aquifer is head-driven and TWO-WAY, and the water network    !
   !      takes it like any other face. -----------------------------------------------------------!
   subroutine test_rk45_bedrock_and_aquifer()
      integer(ik) :: istep
      real(wp)    :: theta_bot0
      !----- (a) sealed bedrock column on RK45. -----------------------------------------------!
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%soil_water_opts%bottom_bc = SOIL_BC_DRAIN ; col_config%soil_water_opts%drain_factor = 0.0_wp
      do istep = 1_ik, 12_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      call check_true('BEDROCK/RK45: whole_water closes', budget%whole_water%n_fail == 0_ik,                          &
              real(budget%whole_water%n_fail, wp))
      call check_true('BEDROCK/RK45: whole_energy closes', budget%whole_energy%n_fail == 0_ik,                        &
              real(budget%whole_energy%n_fail, wp))

      !----- (b) aquifer BC on RK45. A column started DRY must wet from below, and both ledgers     !
      !          must close with the bottom flux running UPWARD. -----------------------------------!
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%soil_water_opts%bottom_bc = SOIL_BC_AQUIFER
      biophys%soil_water%theta(1:col_config%soil%n_active) = 0.15_wp
      theta_bot0 = biophys%soil_water%theta(col_config%soil%n_active)
      do istep = 1_ik, 12_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      call check_true('AQUIFER/RK45: whole_water closes', budget%whole_water%n_fail == 0_ik,                          &
              real(budget%whole_water%n_fail, wp))
      call check_true('AQUIFER/RK45: whole_energy closes', budget%whole_energy%n_fail == 0_ik,                        &
              real(budget%whole_energy%n_fail, wp))
      call check_true('AQUIFER/RK45: dry column wets from below',                                                     &
              biophys%soil_water%theta(col_config%soil%n_active) > theta_bot0,                                            &
              biophys%soil_water%theta(col_config%soil%n_active) - theta_bot0)
      !----- K_bot/Delta with Delta = dz(n)/2 is a fast boundary term the implicit paths absorb and   !
      !      an explicit march may not, so report the substep cost rather than assume it. ------------!
      print '(a,i0,a,i0,a,i0)', '   AQUIFER/RK45 last dt_fast: substeps = ', budget%integ_nsteps,     &
            ' , rejects = ', budget%integ_nrej, ' , rescues = ', budget%rk45_rescue
      col_config%soil_water_opts%bottom_bc = SOIL_BC_DRAIN ; col_config%soil_water_opts%drain_factor = 1.0_wp
   end subroutine test_rk45_bedrock_and_aquifer

   !----- march 24 h of INTEG_RK45 over MOIST free-draining soil and assert the                     !
   !      whole-column water + energy budgets close. -----------------------------------------------!
   subroutine test_rk45_budgets()
      integer(ik) :: istep
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      call check_true('RK45: whole_water closes (n_fail==0)', budget%whole_water%n_fail == 0_ik,                      &
              real(budget%whole_water%n_fail, wp))
      call check_true('RK45: whole_energy closes (n_fail==0)', budget%whole_energy%n_fail == 0_ik,                    &
              real(budget%whole_energy%n_fail, wp))
      call check_true('RK45: ledger fired every dispatched dt_fast', budget%whole_water%n_check == 24_ik * STEPS_PER_HOUR, &
              real(budget%whole_water%n_check, wp))
      call check_true('RK45: whole-energy closes < 1 J', budget%whole_energy%worst < 1.0_wp,                          &
              budget%whole_energy%worst)
      call check_true('RK45: adaptive substeps bounded (~3 expected, cap 64)',                                        &
              budget%integ_nsteps >= 1_ik .and. budget%integ_nsteps <= 64_ik, real(budget%integ_nsteps, wp))
      print '(a,i0,a,i0)', '   RK45 last dt_fast: substeps = ', budget%integ_nsteps,                 &
            ' , rejects = ', budget%integ_nrej
      print '(a,es10.3,a,es10.3,a)', '   (worst whole-column resid: energy= ', budget%whole_energy%worst, &
            ' J/m2  water= ', budget%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets

   !----- march 24 h of INTEG_RK45 over free-draining soil WITH continuous rain                        !
   !      (rainfall>0), mirroring test_column_ark's test_ark_budgets_wet. Asserts the run completes,   !
   !      the soil wets, and BOTH whole-column budgets close at the tight tolerance.                   !
   subroutine test_rk45_budgets_wet()
      integer(ik) :: istep
      real(wp)    :: theta_dry, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      call march_dry_day(theta_dry)                       ! the same day without the rain
      call reset_state()
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-5_wp                         ! ~0.29 mm/hr continuous rain (rainfall>0)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(biophys%soil_water%theta(1:NSL))
      call check_true('RK45 wet: ran a full day of wet steps (no guard error stop)', &
            budget%whole_water%n_check == 24_ik * STEPS_PER_HOUR, &
              real(budget%whole_water%n_check, wp))
      call check_true('RK45 wet: rain wetted the soil column (wetter than the day without rain)',   &
              theta_col1 > theta_dry, theta_col1 - theta_dry)
      call check_true('RK45 wet: whole_energy closes (n_fail==0)', budget%whole_energy%n_fail == 0_ik,                &
              real(budget%whole_energy%n_fail, wp))
      call check_true('RK45 wet: whole_water closes (n_fail==0)', budget%whole_water%n_fail == 0_ik,                  &
              real(budget%whole_water%n_fail, wp))
      call check_true('RK45 wet: whole-energy closes < 1 J', budget%whole_energy%worst < 1.0_wp,                      &
              budget%whole_energy%worst)
      print '(a,es10.3,a,es10.3,a)', '   (RK45 wet worst whole-column resid: energy= ',                &
            budget%whole_energy%worst, ' J/m2  water= ', budget%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets_wet

   !----- march 24 h of INTEG_RK45 with canopy_water_on and a morning rain pulse                         !
   !      (04:45-06:00); assert some of it was intercepted and the budgets close (mirrors                !
   !      test_column_dynamics.f90's RUN 6 / test_column_ark.f90's test_ark_canopy_water). --------------!
   subroutine test_rk45_canopy_water()
      integer(ik) :: istep
      real(wp)    :: surf_water_peak
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%canopy_water_on = .true.
      surf_water_peak = 0.0_wp
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         if (istep > (19_ik * STEPS_PER_HOUR) / 4_ik .and. istep <= 6_ik * STEPS_PER_HOUR) &
            forc%rainfall = 5.0e-5_wp                    ! a morning rain pulse, 04:45-06:00
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         surf_water_peak = max(surf_water_peak, biophys%leaf_surf_water(1) + biophys%wood_surf_water(1))
      end do
      col_config%canopy_water_on = .false.   ! restore default for any test added after this
      call check_true('RK45 canopy water: whole-column water still closes', budget%whole_water%n_fail == 0_ik,        &
              real(budget%whole_water%n_fail, wp))
      call check_true('RK45 canopy water: the morning rain pulse was intercepted', surf_water_peak > 0.0_wp,            &
              surf_water_peak)
      call check_true('RK45 canopy water: whole-column energy stays BOUNDED (film sensible heat not counted)',        &
              budget%whole_energy%worst < 5.0e6_wp, budget%whole_energy%worst)
      print '(a,es10.3,a)', '   (RK45 canopy water peak film=', surf_water_peak, ' kg/m2)'
   end subroutine test_rk45_canopy_water

   !----- march 24 h of INTEG_RK45 with a constant leaf/root-turnover shed-water rate                        !
   !      (biophys%shed_water_rate -- a PATCH-level input, not atmospheric                                   !
   !      forcing, so it is frozen on biophys for the whole day rather than living on forc; distinct from    !
   !      rainfall, which stays 0 throughout), mirroring test_column_ark's test_ark_shed_water: the soil    !
   !      must wet from THIS input alone, and both whole_water AND whole_energy must still close at        !
   !      RK45's own tight tolerance -- energy closing needs no separate wiring (the shed water rides        !
   !      the SAME infiltration-enthalpy treatment, e_infil, every other infiltrating input gets). ----------!
   subroutine test_rk45_shed_water()
      integer(ik) :: istep
      real(wp)    :: theta_dry, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      call march_dry_day(theta_dry)                       ! the same day without the shed water
      call reset_state()
      biophys%shed_water_rate = 8.0e-5_wp                     ! frozen for the whole day (rainfall stays 0)
      biophys%shed_water_enth = 8.0e-5_wp * internal_energy_liquid(T0)   ! shed by tissues at t0
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(biophys%soil_water%theta(1:NSL))
      biophys%shed_water_rate = 0.0_wp   ! restore default for any test added after this
      biophys%shed_water_enth = 0.0_wp
      call check_true('RK45 shed water: leaf/root shed water alone wetted the soil column (wetter than without it)', &
              theta_col1 > theta_dry, theta_col1 - theta_dry)
      call check_true('RK45 shed water: whole-column WATER still closes with shed_water_rate active',                 &
              budget%whole_water%n_fail == 0_ik, real(budget%whole_water%n_fail, wp))
      call check_true('RK45 shed water: whole-column ENERGY still closes (no separate energy wiring needed)',         &
              budget%whole_energy%n_fail == 0_ik, real(budget%whole_energy%n_fail, wp))
      call check_true('RK45 shed water: whole-energy closes < 1 J', budget%whole_energy%worst < 1.0_wp,               &
              budget%whole_energy%worst)
   end subroutine test_rk45_shed_water

   !----- march 00:00-12:00 of INTEG_RK45 with the cohort shrunk to a                                   !
   !      just-recruited size (LAI/WAI/nplant/bleaf/bsap/broot at the magnitude of a 30-yr Ithaca       !
   !      run's first recruitment event). Without the two guards above this reproduces BOTH failures    !
   !      directly: leaf/wood_temp overflowing (T**4 in the RT forcing) and, with that guard alone, a   !
   !      SEPARATE ground_evaporation_rate crash from theta/soil_energy escaping their domain.          !
   !      Restores the shared col_cohort/biophys state on exit. ----------------------------------------!
   subroutine test_rk45_tiny_cohort()
      integer(ik) :: istep, k
      real(wp) :: lai0, wai0, leaf_area0, bleaf0, bsap0, broot0, sap_area0, nplant0
      logical  :: physical
      lai0 = col_cohort%lai(1) ; wai0 = col_cohort%wai(1) ; leaf_area0 = col_cohort%leaf_area(1)
      bleaf0 = col_cohort%bleaf(1) ; bsap0 = col_cohort%bsap(1) ; broot0 = col_cohort%broot(1)
      sap_area0 = col_cohort%sap_area(1) ; nplant0 = col_cohort%nplant(1)

      col_cohort%lai(1) = 0.00265_wp ; col_cohort%wai(1) = 0.000529_wp ; col_cohort%leaf_area(1) = 0.00265_wp
      col_cohort%bleaf(1) = 0.0203507_wp ; col_cohort%bsap(1) = 0.00212943_wp ; col_cohort%broot(1) = 0.0203507_wp
      col_cohort%sap_area(1) = 1.0e-4_wp ; col_cohort%nplant(1) = 0.01_wp

      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      physical = .true.
      do istep = 1_ik, 12_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         !----- The shared fixture hands EVERY cohort a full-canopy absorbed load; for a near-zero-LAI  !
         !      cohort that is unphysical (the real two-stream RT, apply_rt_forcing, scales absorption   !
         !      by leaf area, so a 0.0027-LAI seedling intercepts ~0.1% of the beam). Apply the same     !
         !      Beer-law scaling here, else the test asserts on a leaf absorbing 500 W/m2 through        !
         !      essentially zero area -- which genuinely MUST run hot, and which only a discontinuous    !
         !      floor that zeroes the whole balance would hide. -----------------------------------------!
         block
            real(wp) :: rad_frac
            rad_frac      = 1.0_wp - exp(-0.5_wp * col_cohort%lai(1))
            forc%abs_sw   = forc%abs_sw  * rad_frac
            forc%abs_par  = forc%abs_sw
            forc%abs_lw   = forc%abs_lw  * rad_frac
         end block
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         physical = physical .and. biophys%leaf_temp(1) > 250.0_wp .and. biophys%leaf_temp(1) < 340.0_wp
         physical = physical .and. biophys%wood_temp(1) > 250.0_wp .and. biophys%wood_temp(1) < 340.0_wp
         physical = physical .and. biophys%cas%can_temp > 250.0_wp .and. biophys%cas%can_temp < 340.0_wp
         do k = 1_ik, NSL
            physical = physical .and. biophys%soil_water%theta(k) >= 0.0_wp .and. biophys%soil_water%theta(k) <= 0.6_wp
            physical = physical .and. biophys%soil_energy%soil_temp(k) > 250.0_wp .and. biophys%soil_energy%soil_temp(k) < 340.0_wp
         end do
         physical = physical .and. gpp_coh(1) == gpp_coh(1)   ! not NaN
      end do
      call check_true('RK45 tiny (just-recruited) cohort stays physical + bounded (00:00-12:00)', physical,              &
              biophys%leaf_temp(1))
      call check_true('RK45 tiny cohort: gpp is not NaN', gpp_coh(1) == gpp_coh(1), gpp_coh(1))

      col_cohort%lai(1) = lai0 ; col_cohort%wai(1) = wai0 ; col_cohort%leaf_area(1) = leaf_area0
      col_cohort%bleaf(1) = bleaf0 ; col_cohort%bsap(1) = bsap0 ; col_cohort%broot(1) = broot0
      col_cohort%sap_area(1) = sap_area0 ; col_cohort%nplant(1) = nplant0
   end subroutine test_rk45_tiny_cohort

   !----- ROBUSTNESS STRESS: a very dense canopy (LAI=8) under a hard cold snap (dry -20 C air) at      !
   !      a macro-step of 2*DT_FAST -- the operating-point FAMILY that rails the explicit RK45 surface  !
   !      at scale. Whether the hybrid rescue fires on this exact synthetic point or the stage/y_out    !
   !      clamps alone hold it, the committed state MUST stay physical (CAS/leaf/soil in [250,340] K,   !
   !      GPP finite): a regression that lost surface stability here would rail. (A 30-yr Ithaca run    !
   !      collapses at LAI~1 in its winter 2050 without the rescue and completes healthy with it -- a   !
   !      decades-evolved point no short march reproduces, so the rescue count here is reported for     !
   !      information, not asserted.) Restores col_cohort geometry on exit. ----------------------------!
   subroutine test_rk45_dense_cold_canopy()
      integer(ik) :: istep, k, total_rescue
      real(wp) :: lai0, wai0, la0
      logical  :: physical
      lai0 = col_cohort%lai(1) ; wai0 = col_cohort%wai(1) ; la0 = col_cohort%leaf_area(1)
      ! very dense -> high coupling gain
      col_cohort%lai(1) = 8.0_wp
      col_cohort%wai(1) = 1.5_wp
      col_cohort%leaf_area(1) = 26.0_wp
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      physical = .true. ; total_rescue = 0_ik
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_coldsnap_forcing(istep)
         call column_fast_step(2.0_wp*DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, &
               gpp_coh=gpp_coh)
         total_rescue = total_rescue + budget%rk45_rescue ; budget%rk45_rescue = 0_ik
         physical = physical .and. biophys%cas%can_temp > 250.0_wp .and. biophys%cas%can_temp < 340.0_wp
         physical = physical .and. biophys%leaf_temp(1)  > 250.0_wp .and. biophys%leaf_temp(1)  < 340.0_wp
         do k = 1_ik, NSL
            physical = physical .and. biophys%soil_energy%soil_temp(k) > 250.0_wp .and. biophys%soil_energy%soil_temp(k) < 340.0_wp
         end do
         physical = physical .and. gpp_coh(1) == gpp_coh(1)   ! not NaN
      end do
      call check_true('RK45 dense cold-snap canopy stays physical (surface stability)', physical,                     &
              biophys%cas%can_temp)
      print '(a,i0,a)', '   (dense cold-snap: ', total_rescue, ' RK45->ARK rescues over 24 h)'
      col_cohort%lai(1) = lai0 ; col_cohort%wai(1) = wai0 ; col_cohort%leaf_area(1) = la0
   end subroutine test_rk45_dense_cold_canopy


   !----- Drive the column down toward theta_res: seed just above it, no rainfall, and a dry atmosphere   !
   !      so ground evaporation and transpiration both pull hard. The assertion is the DOMAIN INVARIANT  !
   !      -- theta never commits below theta_res -- plus closure of both whole-column ledgers. Both     !
   !      sinks self-limit before the floor: the root sink falls with the soil->root conductance (and   !
   !      stops where a layer is drier than the roots), and ground evaporation falls with the dry       !
   !      surface layer. This drydown regime is covered by nothing else in the file. ------------------!
   subroutine test_rk45_drydown()
      integer(ik) :: istep
      real(wp)    :: theta_min, res_min, e_worst, w_worst
      theta_seed = 0.085_wp                              ! just above theta_res = 0.078
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      theta_min = 1.0e9_wp
      res_min   = minval(col_config%soil%theta_res(1:NSL))
      do istep = 1_ik, 48_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         forc%shv_atm = 0.5_wp * forc%shv_atm            ! halve the atmospheric humidity: dry the air
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         theta_min = min(theta_min, minval(biophys%soil_water%theta(1:NSL)))
      end do
      e_worst = budget%whole_energy%worst ; w_worst = budget%whole_water%worst
      theta_seed = THETA0                                 ! restore for any test added after this
      call check_true('RK45 drydown: theta never commits below theta_res (Se >= 0 guaranteed)',                       &
              theta_min >= res_min - 1.0e-12_wp, theta_min)
      call check_true('RK45 drydown: whole-column WATER closes', w_worst < 1.0e-6_wp,                                  &
              w_worst)
      call check_true('RK45 drydown: whole-column ENERGY closes', e_worst < 1.0e-3_wp,                                 &
              e_worst)
      print '(a,es10.3,a,es10.3,a)', '   (RK45 drydown: theta_min= ', theta_min, '  theta_res= ',      &
            res_min, ')'
   end subroutine test_rk45_drydown

   !----- SATURATED sealed column: the ARK/RK45 twin of test_column_dynamics RUN 7. A bedrock bottom    !
   !      seeded just below theta_sat under heavy rain both SATURATES (the water network raises the     !
   !      layer's potential and drives the excess to the pond and the next layer) and OVERFLOWS the     !
   !      ponding store (surface runoff). None of the other RK45 runs reach either state -- they are    !
   !      free-draining at theta = 0.30. Both paths move mass and enthalpy as a pair: the pond holds    !
   !      its own enthalpy and runoff leaves at the pond's temperature, so no soil layer is debited     !
   !      for water it never received (docs/science/numerical_scheme.md §5b).                           !
   !      whole_energy closing is the assertion that both are booked with the right sign; the           !
   !      temperature bound is the assertion that the INTERIOR advective faces actually connect the     !
   !      boundary faces -- without them layer 1 accumulates the full infiltration enthalpy while a     !
   !      deeper layer sheds it, and the ledger cannot see it. -----------------------------------------!
   subroutine test_rk45_saturated()
      integer(ik) :: istep
      real(wp)    :: pond_peak, theta_peak, ss_min, ss_max, ood_peak
      theta_seed = 0.428_wp                              ! just below theta_sat = 0.43
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      pond_peak = 0.0_wp ; theta_peak = 0.0_wp
      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp
      ood_peak = 0.0_wp
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-3_wp                         ! ~29 mm/hr: far above the drainage capacity
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         pond_peak  = max(pond_peak,  biophys%soil_water%w_surface)
         theta_peak = max(theta_peak, maxval(biophys%soil_water%theta(1:NSL)))
         ss_min = min(ss_min, biophys%soil_energy%soil_temp(1)) ; ss_max = max(ss_max, biophys%soil_energy%soil_temp(1))
         ood_peak      = max(ood_peak,   budget%theta_ood_max)
      end do
      theta_seed = THETA0                                 ! restore for any test added after this
      call check_true('RK45 saturated: column reached theta_sat (saturation path is live)',                           &
              theta_peak >= 0.43_wp - 1.0e-12_wp, theta_peak)
      call check_true('RK45 saturated: ponding store overflowed (runoff path is live)',                               &
              pond_peak >= col_config%soil_water_opts%w_pond_max - 1.0e-9_wp, pond_peak)
      !----- ENERGY is asserted as a BOUND, not exact closure. Nothing clamps the committed state, but   !
      !      whether a stage clamp bites is trajectory- and compiler-dependent. A clamp that rebuilt the !
      !      COMMITTED soil_energy at a clamped temperature with no ledger term would show here (1.5e5   !
      !      J/m2 on nvfortran, never on ifx). ----------------------------------------------------------!
      call check_true('RK45 saturated: whole-column ENERGY closes (no unbooked clamp on the commit)',                 &
              budget%whole_energy%worst < 1.0e-3_wp, budget%whole_energy%worst)
      !----- WATER closes to MACHINE PRECISION: drainage, ponding and runoff come from the same water   !
      !      network solve whose theta RK45 commits, so every water flow and the state it moves agree   !
      !      (docs/science/numerical_scheme.md §5b). Assert closure, not a bound: there is no known     !
      !      gap to tolerate. --------------------------------------------------------------------------!
      call check_true('RK45 saturated: whole-column WATER closes (n_fail==0)', budget%whole_water%n_fail == 0_ik,     &
              real(budget%whole_water%n_fail, wp))
      call check_true('RK45 saturated: whole-column WATER closes to machine precision',                               &
              budget%whole_water%worst < 1.0e-6_wp, budget%whole_water%worst)
      !----- SOIL-SURFACE TEMPERATURE is the assertion that sees soil heat advected on faces that are     !
      !      not the scheme's own (a vertical misplacement invisible to any whole-column ledger). Two     !
      !      such errors of matched magnitude and OPPOSITE sign can nearly cancel in this number, and     !
      !      one alone moves it by ~50 K; with the faces carrying the flux this stage's own theta is      !
      !      using it sits near 295 K. A 29 mm/hr rain event cannot warm a soil surface past ~300 K, so   !
      !      310 K is a bound with real diagnostic power; 340 K would tolerate the defect. ---------------!
      call check_true('RK45 saturated: soil surface temp stays physical (interior faces on OWN theta)',               &
              ss_min > 250.0_wp .and. ss_max < 310.0_wp, ss_max)
      !----- Soil water is the water network's: a saturated layer holds a little more only through its  !
      !      specific storage, at a positive pressure head that pushes the excess up to the pond and     !
      !      down to the next layer within the step (docs/science/numerical_scheme.md §5b,              !
      !      "Saturation is the water network's"). So theta may commit a hair above theta_sat -- the    !
      !      specific storage times a head of centimetres. ---------------------------------------------!
      call check_true('RK45 saturated: theta above theta_sat only through specific storage (< 1e-4)',                &
              theta_peak <= 0.43_wp + 1.0e-4_wp, theta_peak)
      !----- The constitutive-domain excursion, bounded rather than merely inspected: how far a         !
      !      sub-step's starting theta lies outside [theta_res, theta_sat]. Under the water network it  !
      !      can only lie above saturation, by the specific storage, so the bound is far looser than    !
      !      what is reached (below 1e-4). -------------------------------------------------------------!
      call check_true('RK45 saturated: constitutive-domain excursion stays small', ood_peak < 2.0e-2_wp,              &
              ood_peak)
      print '(a,es10.3)', '   (RK45 saturated: worst constitutive-domain excursion in theta= ', ood_peak
   end subroutine test_rk45_saturated

   subroutine test_rk45_reports_psi_leaf()
      call set_noon_forcing()
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      psi_leaf_probe = 1.0_wp                            ! impossible: psi_leaf is <= 0 by construction
      call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws,           &
                            gpp_coh=gpp_coh, psi_leaf_coh=psi_leaf_probe)
      call check_true('RK45 psi report: the RK45 path itself ran (no ARK rescue)', budget%rk45_rescue == 0_ik,        &
              real(budget%rk45_rescue, wp))
      call check_true('RK45 psi report: psi_leaf_coh is filled on the RK45 success path',                             &
              psi_leaf_probe(1) <= 0.0_wp .and. psi_leaf_probe(1) > -50.0_wp, psi_leaf_probe(1))
   end subroutine test_rk45_reports_psi_leaf

   subroutine test_rk45_dry_uptake_seam()
      integer(ik) :: istep
      real(wp)    :: w_worst, psi_top
      theta_seed = 0.14_wp                               ! psi ~ -10 m: moderately dry
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      do istep = 1_ik, 12_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      w_worst = budget%whole_water%worst
      psi_top = soil_theta_to_head(col_config%soil%retention, biophys%soil_water%theta(1), col_config%soil%theta_sat(1), &
                                    col_config%soil%theta_res(1), col_config%soil%curve_par_a(1),          &
                                    col_config%soil%curve_par_n(1))
      theta_seed = THETA0
      call check_true('RK45 dry seam: the fixture really is moderately dry (-150 m < psi < -3 m)',                    &
              psi_top < -3.0_wp .and. psi_top > -150.0_wp, psi_top)
      call check_true('RK45 dry seam: the RK45 path itself ran (no ARK rescue)', budget%rk45_rescue == 0_ik,          &
              real(budget%rk45_rescue, wp))
      call check_true('RK45 dry seam: whole-column WATER closes to round-off (soil debit == wood credit)',            &
              w_worst < 1.0e-8_wp, w_worst)
   end subroutine test_rk45_dry_uptake_seam


   subroutine test_campbell_predawn_seed()
      real(wp), parameter :: THETA_SAT = 0.43_wp, THETA_RES = 0.0_wp, PSI_SAT = -0.26_wp, B = 5.65_wp
      type(soil_params_t) :: vg_soil
      real(wp) :: cdiag(N_CO_FAST_ROWS, N), seed
      vg_soil = col_config%soil
      call build_soil_hydr_params(NSL, SOIL_RETENTION_CAMPBELL, test_bottoms, THETA_SAT, THETA_RES,     &
                                  2.89e-6_wp, PSI_SAT, B, -3.37_wp, col_config%soil)
      call set_noon_forcing()
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      cdiag = 0.0_wp
      call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws,           &
                            gpp_coh=gpp_coh, cdiag=cdiag)
      col_config%soil = vg_soil
      seed = grav_head * soil_theta_to_head(SOIL_RETENTION_CAMPBELL, theta_seed, THETA_SAT, THETA_RES, PSI_SAT, B)
      call check_true('Campbell seed: the cohort has no predawn potential of its own yet',                            &
              col_cohort%dmax_psi_leaf(1) > 0.0_wp, col_cohort%dmax_psi_leaf(1))
      call check_true('Campbell seed: the top layer is unsaturated (a seed below zero)', seed < 0.0_wp, seed)
      call check_true('Campbell seed: the predawn seed is the Campbell curve''s top-layer potential',                 &
              abs(cdiag(CO_PSI_PREDAWN, 1) - seed) <= 1.0e-12_wp * abs(seed), cdiag(CO_PSI_PREDAWN, 1))
   end subroutine test_campbell_predawn_seed

   subroutine reset_state()
      integer(ik) :: kk
      if (allocated(biophys%leaf_temp)) deallocate(biophys%leaf_temp)
      call alloc_patch_biophys(biophys, N, T0, 0.008_wp, 400.0_wp, T0)
      biophys%leaf_water_mass(1:N) = tissue_psi_to_water(PSI_FIXTURE, col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1:N))
      biophys%wood_water_mass(1:N) = tissue_psi_to_water(PSI_FIXTURE, col_config%hydraulics_table%pft(1)%wood_curve, &
            col_cohort%bsap(1:N) + col_cohort%broot(1:N))
      !----- ...and the tissues' energy at t0, the state the column advances. ----------------------!
      biophys%leaf_energy(1:N) = tissue_temp_to_energy(tissue_dry_heat_capacity(col_cohort%nplant(1:N),                &
                                    col_cohort%bleaf(1:N), tissue_c_leaf),                                    &
                                    max(biophys%leaf_water_mass(1:N), 0.0_wp) * col_cohort%nplant(1:N), T0)
      biophys%wood_energy(1:N) = tissue_temp_to_energy(tissue_dry_heat_capacity(col_cohort%nplant(1:N),                &
                                    col_cohort%bwood(1:N), tissue_c_wood_dry),                                &
                                    max(biophys%wood_water_mass(1:N), 0.0_wp) * col_cohort%nplant(1:N), T0)
      budget = column_budget_t()
      biophys%soil_water%theta(1:NSL) = theta_seed
      do kk = 1_ik, NSL
         biophys%soil_energy%soil_energy(kk) = temp_to_internal_energy(col_config%soil_thermal%soil_dry_heat_capacity(kk), &
                                      theta_seed * rho_h2o, T0, 1.0_wp)
         biophys%soil_energy%soil_temp(kk)   = T0
      end do
   end subroutine reset_state

   subroutine set_noon_forcing()
      forc%abs_sw = 450.0_wp ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 70.0_wp ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_temp_to_enthalpy(295.0_wp, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, 295.0_wp, forc%shv_atm, forc%co2_atm)   ! MO reference = this air, not 298.15 K
   end subroutine set_noon_forcing

   subroutine set_diurnal_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * DT_FAST
      cosz  = solar_cosz(sim_date, t_sec, LAT)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
      forc%abs_sw = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_temp_to_enthalpy(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! MO reference = this air, not 298.15 K
   end subroutine set_diurnal_forcing

   !----- HARD cold snap: very cold, very dry reference air (~248-256 K) with a weak winter sun over the !
   !      warm (t0=288 K) soil the reset_state seeds -- a large canopy-through gradient that makes the    !
   !      dense-canopy leaf<->CAS coupling stiff enough to rail the explicit RK45 surface. --------------!
   subroutine set_coldsnap_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * DT_FAST
      cosz  = solar_cosz(sim_date, t_sec, LAT)
      t_air = 252.0_wp + 4.0_wp * (cosz - 0.3_wp)                 ! ~248-256 K (a −20 C cold snap)
      forc%abs_sw = 120.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 20.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_temp_to_enthalpy(t_air, 0.001_wp)   ! very dry cold air
      forc%shv_atm = 0.001_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! MO reference = this air, not 298.15 K
   end subroutine set_coldsnap_forcing


   !----- The same 24 h with no rain and no shed water, from the state the caller set up: the       !
   !      column's summed water content at the end. Free drainage empties a wet column whatever      !
   !      falls on it, so an input wets the soil when the column ends wetter than this day does. --!
   subroutine march_dry_day(theta_sum)
      real(wp), intent(out) :: theta_sum
      integer(ik) :: istep
      biophys%shed_water_rate = 0.0_wp
      do istep = 1_ik, 24_ik * STEPS_PER_HOUR
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         call column_fast_step(DT_FAST, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_sum = sum(biophys%soil_water%theta(1:NSL))
   end subroutine march_dry_day

end program test_column_rk45
