! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_column_ark -- the INTEG_ARK wiring of the coupled ESDIRK2 integrator (config name `ark`)        !
! into column_fast_step (docs/science/numerical_scheme.md §3).                                         !
!   A. GPP PARITY: build_column_frozen runs the SAME frozen pre-pass under both schemes, so gpp_coh    !
!      from the ARK path is bit-identical to RK45's on step 1 (identical initial state).               !
!   B. PHYSICAL: a 24-step dry-window march under INTEG_ARK keeps the CAS / soil / leaf / psi state    !
!      finite, bounded, and sub-saturated (the integrator + unpack are correct end to end).            !
!   C-E. the fixed-substep path, the leaf<->CAS Newton switch, psi_leaf convergence in dt_fast, and    !
!      the seven conservation budgets over a 24 h dry march.                                           !
!   F. PRECIP>0: a wet march runs, wets the soil, and closes the budgets (ENERGY to machine            !
!      precision through the boundary water-enthalpy advection, WATER at the ledger tolerance).        !
!==========================================================================================!
program test_column_ark
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o, cp_liq, n_soil_layer_max
   use meds_config,              only : meds_config_t, INTEG_ARK, INTEG_RK45
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_internal_energy, cas_temp_of_enthalpy, &
                                        sat_specific_humidity, internal_energy_liquid
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_fast_types, only : patch_biophys_t, alloc_patch_biophys
   use meds_column_view,       only : column_cohort_init
   use meds_water_retention, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use meds_canopy_types, only : set_aero_env_atm
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_column_config, only : build_soil_hydr_params, build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config, apply_canopy_film_config
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_fast_step,          only : column_fast_step
   use meds_fast_frozen,        only : build_column_frozen
   use meds_fast_types,         only : column_frozen_t, column_state_t, fast_workspace_t
   use meds_water_retention,     only : psi_from_water_content, water_content
   use meds_test_support, only : build_test_config, check_close, check_true, test_report
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks were written on

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik
   real(wp),    parameter :: dt_fast = 150.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   !----- Durations and times of day are written in hours and turned into steps here, so a change of  !
   !      dt_fast keeps what each test covers: set_diurnal_forcing puts step i at (i - 1/2)*dt_fast    !
   !      after midnight, and the sun follows the clock. -------------------------------------------!
   integer(ik), parameter :: steps_per_hour = nint(3600.0_wp / dt_fast, ik)
   !----- seed moisture, a VARIABLE so the saturated test can drive the column to theta_sat.   !
   real(wp) :: theta_seed = theta0
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
   real(wp)    :: gpp_rk45(n), gpp_ark(n), gpp_coh(n), tcas, qsat, worst_super, tcas_1, tcas_8
   real(wp)    :: psi_leaf_diag
   integer(ik) :: is, k
   logical     :: physical

   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- column setup (test_column_rk45 mirrors it). ----------------------------------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   !----- On-allometry fixture, a 20 cm tree at 224 stems/ha (LAI ~3): born through the canonical !
   !      path (init_cohort + set_cohort_size) and gathered by the production gather, so every    !
   !      field -- height, leaf area, the wood carbon the wood heat capacity is sized from -- is  !
   !      mutually consistent. -------------------------------------------------------------------!
   call column_cohort_init(col_cohort, cfg%pft, [1_ik], [20.0_wp], [0.0224_wp])
   !----- The gathered per-PFT aboveground fraction must reach the view the fast loop reads.   !
   !      The kernel takes the trait, and this pins the chain that delivers it (PFT table ->   !
   !      init_cohort -> cohort block -> the view).                                            !
   call check_close(col_cohort%aboveground_frac(1), cfg%pft%aboveground_frac(1), 1.0e-14_wp,   &
                    'aboveground_frac did not reach the column view from the PFT table')
   call geometric_layer_bottoms(2.0_wp, nsl, 3.0_wp, test_bottoms)
   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, col_config%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, col_config%soil_thermal)
   call apply_hydraulics_config(cfg%hydraulics, cfg%pft, col_config%hydraulics_table)
   call apply_canopy_film_config(cfg%pft, col_config%leaf_surf_water_max, col_config%wood_surf_water_max)
   call build_leaf_photo_table(cfg, col_config%leaf_photo)
   col_config%integrator = build_integrator_opts(cfg)
   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   print '(a)', '[test_column_ark]'

   !=== A. GPP parity on step 1 (identical initial state; ARK pre-pass == RK45 pre-pass). ======!
   call set_noon_forcing()
   call reset_state()
   cfg%time_integrator = INTEG_RK45
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   gpp_rk45 = gpp_coh
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   gpp_ark = gpp_coh
   call check_true('ARK pre-pass gpp bit-identical to RK45 (one shared build_column_frozen)',                         &
           abs(gpp_ark(1) - gpp_rk45(1)) < 1.0e-12_wp, abs(gpp_ark(1) - gpp_rk45(1)))
   call check_true('ARK midday gpp > 0', gpp_ark(1) > 0.0_wp, gpp_ark(1))

   !=== Emissivity. The energy balance's emission slope 4*eps*sigma*T^3 takes each cohort's PFT     !
   !    leaf_emissivity and wood_emissivity, the values the radiation solver emits with. Distinct   !
   !    values for leaf and wood catch a leaf-for-wood swap. ---------------------------------------!
   block
      type(column_frozen_t) :: fz
      type(column_state_t)  :: ys
      real(wp) :: eps_leaf0, eps_wood0
      eps_leaf0 = cfg%pft%leaf_emissivity(col_cohort%pft(1)) ; eps_wood0 = cfg%pft%wood_emissivity(col_cohort%pft(1))
      cfg%pft%leaf_emissivity = 0.93_wp ; cfg%pft%wood_emissivity = 0.87_wp
      call reset_state()
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                               ws%spac, n, nsl, fz, ys)
      call check_close(fz%tissue%leaf_emiss(1), 0.93_wp, 0.0_wp, 'energy balance uses the PFT leaf emissivity')
      call check_close(fz%tissue%wood_emiss(1), 0.87_wp, 0.0_wp, 'energy balance uses the PFT wood emissivity')
      cfg%pft%leaf_emissivity = eps_leaf0 ; cfg%pft%wood_emissivity = eps_wood0
      call reset_state()
   end block

   !=== B. A dry-window march under INTEG_ARK stays physical + bounded + sub-saturated. ========!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   physical = .true. ; worst_super = -1.0_wp
   do is = 1_ik, 6_ik * steps_per_hour
      call set_diurnal_forcing(is)
      call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      tcas = biophys%cas%can_temp
      qsat = sat_specific_humidity(tcas, aenv%press)
      worst_super = max(worst_super, biophys%cas%can_shv - qsat)
      physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. biophys%cas%can_shv > 0.0_wp
      physical = physical .and. biophys%leaf_temp(1) > 250.0_wp .and. biophys%leaf_temp(1) < 340.0_wp
      do k = 1_ik, nsl
         physical = physical .and. biophys%soil_water%theta(k) > 0.0_wp .and. biophys%soil_water%theta(k) < 0.6_wp
         physical = physical .and. biophys%soil_energy%soil_temp(k) > 260.0_wp .and. biophys%soil_energy%soil_temp(k) < 340.0_wp
      end do
      !----- psi is not persisted state: diagnose it from the persisted leaf_water_mass             !
      !      (docs/science/plant_hydraulics.md §2, eqs. (4)-(5)). ----------------------------------!
      psi_leaf_diag = psi_from_water_content(biophys%leaf_water_mass(1), col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1))
      physical = physical .and. psi_leaf_diag < 0.5_wp .and. psi_leaf_diag > -12.0_wp
   end do
   call check_true('INTEG_ARK dry-window march stays physical + bounded (00:00-06:00)', physical, biophys%cas%can_temp)
   call check_true('INTEG_ARK CAS stays sub-saturated', worst_super <= 1.0e-4_wp, worst_super)

   !=== C. Fixed-substep (GPU-lockstep) path also runs + stays physical. =======================!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .false. ; cfg%ark_fixed_substep = 4_ik
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call set_noon_forcing()
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   call check_true('INTEG_ARK fixed-substep path physical',                                                           &
           biophys%cas%can_temp > 270.0_wp .and. biophys%cas%can_temp < 325.0_wp, biophys%cas%can_temp)

   !=== D. ark_coupled reaches the inner solver (np<=1 baseline vs np>1 Newton must differ). ======!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true. ; cfg%ark_rtol = 1.0e-4_wp
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call set_noon_forcing()
   cfg%ark_coupled = .false.
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   tcas_1 = biophys%cas%can_temp
   call reset_state() ; call set_noon_forcing()
   cfg%ark_coupled = .true.
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
   tcas_8 = biophys%cas%can_temp
   call check_true('ARK ark_coupled reaches ark2 (baseline vs Newton differ)', abs(tcas_1 - tcas_8) > 1.0e-4_wp,      &
           abs(tcas_1-tcas_8))
   call check_true('both niter paths physical', tcas_1 > 270.0_wp .and. tcas_8 < 325.0_wp, tcas_8)

   !=== D2. psi_leaf must be CONVERGED IN dt_fast. ==============================================!
   call test_psi_dt_convergence()

   !=== E. WHOLE-COLUMN CONSERVATION LEDGER: the 7 budgets close over a 24 h dry diurnal march. ==!
   call test_ark_budgets(.true.)      ! adaptive path (accept/reject accumulation)
   call test_ark_budgets(.false.)     ! fixed-substep path (uniform accumulation)

   !=== F. PRECIP>0: a WET diurnal march runs, wets the soil, and still closes the budgets -- ENERGY   !
   !       to machine precision (the rain/runoff/drainage advection is a frozen fixed source), WATER   !
   !       at the ledger tolerance. ===================================================================!
   call test_ark_budgets_wet()

   !=== G. CANOPY-SURFACE WATER (opt-in): a diurnal march with a morning rain pulse actually gets         !
   !       intercepted; whole_water closes exactly, whole_energy stays BOUNDED (the store is valued at    !
   !       one fixed t_film_valuation reference rather than a prognostic surface-water temperature, so    !
   !       its sensible heat at the tissue temperature offset is not counted; see docs/ROADMAP.md;        !
   !       mirrors RUN 6 in test_column_dynamics.f90, same bound). Proves the WIRING (interception->      !
   !       film->CAS->ledgers), not the wetted-fraction algebra itself (unit-tested in                    !
   !       test_surface_energy.f90). =====================================================================!
   call test_ark_canopy_water()

   !=== H. LEAF/ROOT-TURNOVER SHED WATER: a constant shed_water_rate                                     !
   !       (distinct from rainfall) wets the soil and both whole_water/whole_energy still close. ========!
   call test_ark_shed_water()
   call test_ark_saturated()
   call test_ark_aquifer()
   call test_wood_prognostic(INTEG_ARK, 'ARK ')
   call test_refilling_wood_isothermal()
   call test_wood_capacity_follows_store()

   call test_report('test_column_ark')

contains

   !----- Prognostic WOOD: veg_energy_balance relaxes the wood's heat store (store_hcap_per_dt > 0)     !
   !      exactly (docs/science/vegetation_energy_dynamics.md eqs. (3)-(4)), on ARK and RK45.           !
   !      Two assertions:                                                                               !
   !        (a) the whole-column energy ledger still closes -- the store delta AND the wood net         !
   !            radiation both enter it, since surface_derivs does not fold wood into                   !
   !            coh_rnet when the frozen diagnostic inputs are zeroed;                                  !
   !        (b) wood LAGS the canopy air through the diel cycle, which is the feature.                  !
   subroutine test_wood_prognostic(integ, tag)
      integer(ik),      intent(in) :: integ
      character(len=4), intent(in) :: tag
      integer(ik) :: istep
      real(wp)    :: dmax_lag
      call reset_state()
      cfg%time_integrator = integ
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      dmax_lag = 0.0_wp
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         dmax_lag = max(dmax_lag, abs(biophys%wood_temp(1) - biophys%cas%can_temp))
      end do
      call check_true(trim(tag)//' PROG-WOOD: whole_energy closes', budget%whole_energy%n_fail == 0_ik,               &
              real(budget%whole_energy%n_fail, wp))
      call check_true(trim(tag)//' PROG-WOOD: whole_water closes', budget%whole_water%n_fail == 0_ik,                 &
              real(budget%whole_water%n_fail, wp))
      call check_true(trim(tag)//' PROG-WOOD: wood temperature lags the CAS', dmax_lag > 1.0e-3_wp, dmax_lag)
      call check_true(trim(tag)//' PROG-WOOD: wood temperature physical',                                             &
              biophys%wood_temp(1) > 200.0_wp .and. biophys%wood_temp(1) < 350.0_wp, biophys%wood_temp(1))
   end subroutine test_wood_prognostic

   !----- Wood refilling from the soil must not heat itself. Dried wood over moist soil at night,         !
   !      with soil, wood and leaf at one temperature, refills -- and the water it takes in arrives at    !
   !      its own temperature, so the wood's heat input must be zero. Counting the water's whole          !
   !      enthalpy (~4.4e5 J/kg from the liquid datum) as heat puts the wood 13 K above the canopy air    !
   !      at a storm front (BCI, 2017-04-17). The water stores take that enthalpy instead, and what       !
   !      the soil gives up must equal what the tissues and their water stores gain. --------------------!
   subroutine test_refilling_wood_isothermal()
      type(column_frozen_t) :: fz
      type(column_state_t)  :: ys
      real(wp) :: refill, gained, given
      call reset_state()
      call set_noon_forcing()
      forc%abs_sw = 0.0_wp ; forc%abs_par = 0.0_wp ; forc%abs_sw_ground = 0.0_wp   ! night
      forc%enthalpy_atm = cas_enthalpy_of_temp(t0, 0.008_wp)
      biophys%wood_water_mass(1:n) = 0.8_wp * biophys%wood_water_mass(1:n)        ! dried wood
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                               ws%spac, n, nsl, fz, ys)
      refill = sum((fz%plant%uptake_frozen(1:n) - fz%plant%sapflow_frozen(1:n)) * col_cohort%nplant(1:n))
      call check_true('dried wood over moist soil refills', refill > 0.0_wp, refill)
      call check_true('the refill brings the wood no heat (soil, wood and leaf at one temperature)',            &
                      maxval(abs(fz%tissue%q_wood_net(1:n))) < 1.0e-9_wp * refill * internal_energy_liquid(t0), &
                      maxval(abs(fz%tissue%q_wood_net(1:n))))
      gained = sum(fz%tissue%q_wood_net(1:n)) + sum(fz%tissue%qwflux_wl(1:n)) + fz%tissue%water_store_enth
      given  = sum(fz%roots%qloss_frozen(1:n))
      call check_close(gained, given, 1.0e-12_wp * abs(given),                                                &
                       'what the soil gives up = what the tissues and their water stores gain')
      call reset_state()
   end subroutine test_refilling_wood_isothermal

   !----- The wood's heat capacity is its dry wood at the dry specific heat plus the sapwood share of   !
   !      the water network's store (docs/science/vegetation_energy_dynamics.md eq. (8)), so a store     !
   !      drained by a fifth holds a fifth less of that water's heat, and the fine roots' share is not  !
   !      counted at all. ------------------------------------------------------------------------------!
   subroutine test_wood_capacity_follows_store()
      type(column_frozen_t) :: fz
      type(column_state_t)  :: ys
      real(wp) :: w_store, w_sap, cap_dry, cap_full, cap_drained
      call reset_state()
      call set_noon_forcing()
      w_store = biophys%wood_water_mass(1)
      w_sap   = w_store * col_cohort%bsap(1) / (col_cohort%bsap(1) + col_cohort%broot(1))
      cap_dry = col_cohort%bwood(1) * col_cohort%nplant(1) * 2.0_wp * col_config%veg_thermal%c_wood_dry
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                               ws%spac, n, nsl, fz, ys)
      cap_full = fz%tissue%wood_hcap_per_dt(1) * dt_fast
      call check_close(cap_full, cap_dry + col_cohort%nplant(1) * w_sap * cp_liq, 1.0e-12_wp * cap_full,       &
                       'wood heat capacity = dry wood + the sapwood share of the store''s water')
      call check_true('the fine roots hold part of the store', col_cohort%broot(1) > 0.0_wp, col_cohort%broot(1))
      biophys%wood_water_mass(1) = 0.8_wp * w_store
      call build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                               ws%spac, n, nsl, fz, ys)
      cap_drained = fz%tissue%wood_hcap_per_dt(1) * dt_fast
      call check_close(cap_full - cap_drained, 0.2_wp * col_cohort%nplant(1) * w_sap * cp_liq,                  &
                       1.0e-9_wp * cap_full, 'a store drained by a fifth holds a fifth less water heat')
      call reset_state()
   end subroutine test_wood_capacity_follows_store

   !----- The aquifer bottom BC is a head-driven, two-way boundary with no prognostic state, and the   !
   !      ARK commits the water network's theta, so it takes the aquifer like any other face. A column !
   !      started DRY must wet from below with both ledgers closed. -----------------------------------!
   subroutine test_ark_aquifer()
      integer(ik) :: istep
      real(wp)    :: theta_bot0
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%soil_water_opts%bottom_bc = SOIL_BC_AQUIFER
      biophys%soil_water%theta(1:col_config%soil%n_active) = 0.15_wp
      theta_bot0 = biophys%soil_water%theta(col_config%soil%n_active)
      do istep = 1_ik, 12_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      call check_true('AQUIFER/ARK: whole_water closes', budget%whole_water%n_fail == 0_ik,                           &
              real(budget%whole_water%n_fail, wp))
      call check_true('AQUIFER/ARK: whole_energy closes', budget%whole_energy%n_fail == 0_ik,                         &
              real(budget%whole_energy%n_fail, wp))
      call check_true('AQUIFER/ARK: dry column wets from below',                                                      &
              biophys%soil_water%theta(col_config%soil%n_active) > theta_bot0,                                            &
              biophys%soil_water%theta(col_config%soil%n_active) - theta_bot0)
      col_config%soil_water_opts%bottom_bc = SOIL_BC_DRAIN ; col_config%soil_water_opts%drain_factor = 1.0_wp
   end subroutine test_ark_aquifer

   !----- march 24 h of INTEG_ARK over MOIST free-draining soil (rainfall==0) and                   !
   !       assert the 7 conservation budgets close. ------------------------------------------------!
   subroutine test_ark_budgets(adaptive)
      logical, intent(in) :: adaptive
      integer(ik) :: istep
      character(len=6) :: tag
      tag = merge('adapt ', 'fixed ', adaptive)
      call reset_state()                               ! theta0=0.30 (moist), budget zeroed
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = adaptive
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_fixed_substep = 4_ik ; cfg%ark_coupled = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)               ! rainfall==0 always (dry); diurnal SW
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      call check_true('ARK '//trim(tag)//': all 7 budgets close (n_fail==0)',                                         &
              budget%cas_energy%n_fail == 0_ik .and. budget%cas_water%n_fail == 0_ik                                  &
              .and.               budget%cas_co2%n_fail == 0_ik .and. budget%soil_energy%n_fail == 0_ik               &
              .and.               budget%soil_water%n_fail == 0_ik .and. budget%whole_energy%n_fail == 0_ik           &
              .and.               budget%whole_water%n_fail == 0_ik, real(budget%whole_energy%n_fail, wp))
      call check_true('ARK '//trim(tag)//': ledger fired every dispatched dt_fast',                                   &
              budget%cas_energy%n_check == 24_ik * steps_per_hour, real(budget%cas_energy%n_check, wp))
      call check_true('ARK '//trim(tag)//': CAS energy machine-precision closure',                                    &
              budget%cas_energy%worst < 1.0e-3_wp, budget%cas_energy%worst)
      call check_true('ARK '//trim(tag)//': whole-energy closes < 1 J', budget%whole_energy%worst < 1.0_wp,           &
              budget%whole_energy%worst)
   end subroutine test_ark_budgets

   !----- march 24 h of INTEG_ARK over free-draining soil WITH continuous rain                         !
   !       (rainfall>0: boundary water-enthalpy advection + the pond). Assert the run completes, the   !
   !       soil wets, and the budgets close (ENERGY to machine precision, WATER at the ledger          !
   !       tolerance). --------------------------------------------------------------------------------!
   subroutine test_ark_budgets_wet()
      integer(ik) :: istep
      real(wp)    :: theta_dry, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_coupled = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      call march_dry_day(theta_dry)                       ! the same day without the rain
      call reset_state()
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-5_wp                         ! ~0.29 mm/hr continuous rain (rainfall>0)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(biophys%soil_water%theta(1:nsl))
      call check_true('ARK wet: ran a full day of wet steps (no guard error stop)', &
            budget%cas_energy%n_check == 24_ik * steps_per_hour, &
              real(budget%cas_energy%n_check, wp))
      call check_true('ARK wet: rain wetted the soil column (wetter than the day without rain)',   &
              theta_col1 > theta_dry, theta_col1 - theta_dry)
      call check_true('ARK wet: ENERGY budgets close (incl. advection)',                                              &
              budget%cas_energy%n_fail == 0_ik .and. budget%soil_energy%n_fail == 0_ik                                &
              .and.               budget%whole_energy%n_fail == 0_ik, budget%whole_energy%worst)
      call check_true('ARK wet: WATER budgets close (soil and whole-column)',                                         &
              budget%soil_water%n_fail == 0_ik .and. budget%whole_water%n_fail == 0_ik, budget%whole_water%worst)
      call check_true('ARK wet: soil energy machine-precision closure', budget%soil_energy%worst < 1.0e-3_wp,         &
              budget%soil_energy%worst)
   end subroutine test_ark_budgets_wet

   !----- march 24 h of INTEG_ARK with canopy_water_on and a morning rain pulse                         !
   !      (04:45-06:00); assert some of it was intercepted and the budgets close (mirrors               !
   !      test_column_dynamics.f90's RUN 6). -----------------------------------------------------------!
   subroutine test_ark_canopy_water()
      integer(ik) :: istep
      real(wp)    :: surf_water_peak
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_coupled = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%canopy_water_on = .true.
      surf_water_peak = 0.0_wp
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         if (istep > (19_ik * steps_per_hour) / 4_ik .and. istep <= 6_ik * steps_per_hour) &
            forc%rainfall = 5.0e-5_wp                    ! a morning rain pulse, 04:45-06:00
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         surf_water_peak = max(surf_water_peak, biophys%leaf_surf_water(1) + biophys%wood_surf_water(1))
      end do
      col_config%canopy_water_on = .false.   ! restore default for any test added after this
      call check_true('ARK canopy water: whole-column water still closes', budget%whole_water%n_fail == 0_ik,         &
              real(budget%whole_water%n_fail, wp))
      call check_true('ARK canopy water: the morning rain pulse was intercepted', surf_water_peak > 0.0_wp,             &
              surf_water_peak)
      call check_true('ARK canopy water: whole-column energy stays BOUNDED (film sensible heat not counted)',         &
              budget%whole_energy%worst < 5.0e6_wp, budget%whole_energy%worst)
      print '(a,es10.3,a)', '   (ARK canopy water peak film=', surf_water_peak, ' kg/m2)'
   end subroutine test_ark_canopy_water

   !----- march 24 h of INTEG_ARK with a constant leaf/root-turnover shed-water rate                         !
   !      (biophys%shed_water_rate -- a PATCH-level input, not atmospheric                                   !
   !      forcing, so it is frozen on biophys for the whole day rather than living on forc; distinct from    !
   !      rainfall, which stays 0 throughout): the soil must wet from THIS input alone, and both              !
   !      whole_water AND whole_energy must still close -- energy closing needs NO separate wiring of        !
   !      its own: the shed water's enthalpy rides the SAME infiltration-enthalpy treatment (e_infil)        !
   !      every other infiltrating input gets, once mixed into hforc%precip_ground                           !
   !      by build_column_frozen. ---------------------------------------------------------------------------!
   subroutine test_ark_shed_water()
      integer(ik) :: istep
      real(wp)    :: theta_dry, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_coupled = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      call march_dry_day(theta_dry)                       ! the same day without the shed water
      call reset_state()
      biophys%shed_water_rate = 8.0e-5_wp                     ! frozen for the whole day (rainfall stays 0)
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(biophys%soil_water%theta(1:nsl))
      biophys%shed_water_rate = 0.0_wp   ! restore default for any test added after this
      call check_true('ARK shed water: leaf/root shed water alone wetted the soil column (wetter than without it)', &
              theta_col1 > theta_dry, theta_col1 - theta_dry)
      call check_true('ARK shed water: whole-column WATER still closes with shed_water_rate active',                  &
              budget%whole_water%n_fail == 0_ik, real(budget%whole_water%n_fail, wp))
      call check_true('ARK shed water: whole-column ENERGY still closes (no separate energy wiring needed)',          &
              budget%whole_energy%n_fail == 0_ik, real(budget%whole_energy%n_fail, wp))
   end subroutine test_ark_shed_water


   !----- SATURATED column: the ARK twin of test_column_dynamics RUN 7. A column seeded just below       !
   !      theta_sat under heavy rain both SATURATES (the water network raises the layer's potential      !
   !      and drives the excess to the pond and the next layer) and OVERFLOWS the ponding store          !
   !      (surface runoff). None of the other ARK runs reach either state -- they are free-draining at   !
   !      theta = 0.30. Both paths move mass and enthalpy as a pair: the pond holds its own enthalpy     !
   !      and runoff leaves at the pond's temperature, so no soil layer is debited for water it never    !
   !      received (docs/science/numerical_scheme.md §5b).                                               !
   !      whole_energy closing is the assertion that both are booked with the right sign; the            !
   !      temperature bound is the assertion that the INTERIOR advective faces actually connect the      !
   !      boundary faces -- without them layer 1 accumulates the full infiltration enthalpy while a      !
   !      deeper layer sheds it, and the ledger cannot see it. ------------------------------------------!
   subroutine test_ark_saturated()
      integer(ik) :: istep
      real(wp)    :: pond_peak, theta_peak, ss_min, ss_max
      !----- FREE-DRAIN, not bedrock (test_column_rk45 covers the sealed column): the route to         !
      !      saturation here is rain far above the column's drainage capacity                          !
      !      (ksat = 2.89e-6 m/s ~ 2.9e-3 kg/m2/s) rather than a sealed bottom. -----------------------!
      theta_seed = 0.428_wp                              ! just below theta_sat = 0.43
      call reset_state()
      cfg%time_integrator = INTEG_ARK
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      pond_peak = 0.0_wp ; theta_peak = 0.0_wp
      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-3_wp                         ! ~29 mm/hr: far above the drainage capacity
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         pond_peak  = max(pond_peak,  biophys%soil_water%w_surface)
         theta_peak = max(theta_peak, maxval(biophys%soil_water%theta(1:nsl)))
         ss_min = min(ss_min, biophys%soil_energy%soil_temp(1)) ; ss_max = max(ss_max, biophys%soil_energy%soil_temp(1))
      end do
      theta_seed = theta0                                 ! restore for any test added after this
      call check_true('ARK saturated: column reached theta_sat (saturation path is live)',                            &
              theta_peak >= 0.43_wp - 1.0e-12_wp, theta_peak)
      call check_true('ARK saturated: ponding store overflowed (runoff path is live)',                                &
              pond_peak >= col_config%soil_water_opts%w_pond_max - 1.0e-9_wp, pond_peak)
      call check_true('ARK saturated: whole-column ENERGY still closes through saturation + runoff',                  &
              budget%whole_energy%n_fail == 0_ik, budget%whole_energy%worst)
      call check_true('ARK saturated: whole-column WATER still closes', budget%whole_water%n_fail == 0_ik,            &
              budget%whole_water%worst)
      call check_true('ARK saturated: soil surface temp stays physical (interior faces connected)',                   &
              ss_min > 250.0_wp .and. ss_max < 340.0_wp, ss_max)
   end subroutine test_ark_saturated

   !=======================================================================================!
   ! psi_leaf must be CONVERGED IN dt_fast, and the plant total must still conserve.
   !
   ! WHAT THIS CATCHES. The leaf is the smallest store in the column, so any mismatch between the
   ! transpiration charged to it and the sapflow credited to it lands there, and psi_leaf runs away
   ! with dt_fast. An explicit plant-water step charged the ENDPOINT-refreshed transpiration while
   ! credited the STATE-n sapflow carries a flux inconsistency growing LINEARLY in dt: psi_leaf
   ! goes from -0.334 MPa at 5 s to -1.788 MPa at 900 s on the 3 h midday probe while every OTHER
   ! state and every flux converges. The water network solves the leaf implicitly with the soil,
   ! and the energy stages hold transpiration at the network's value, so the leaf is charged and
   ! credited the same step (docs/science/numerical_scheme.md §5a).
   !
   ! No conservation ledger sees this. The whole-column water budget closes either way, because
   ! the plant store absorbing a flux mismatch IS a closed budget -- which is exactly why this
   ! needs its own test rather than relying on section E.
   !
   ! The guard: march the SAME 3 h window at a short and a long dt_fast and require the leaf
   ! potentials to agree. The mismatch makes them differ by ~1 MPa; the tolerance below is far
   ! tighter than that and far looser than the ~0.005 MPa actually achieved, so it fails loudly on
   ! a regression without pinning the exact discretisation.
   !=======================================================================================!
   subroutine test_psi_dt_convergence()
      !----- the window must be in DAYLIGHT: the defect is driven by transpiration, so a night      !
      !      window has no signal at all and the test would pass on a broken scheme. Start at solar !
      !      noon. ---------------------------------------------------------------------------------!
      real(wp), parameter :: t_start  = 43200.0_wp        ! solar noon
      real(wp), parameter :: t_window = 10800.0_wp        ! 3 h
      real(wp) :: psi_short, psi_long, w_tot_short, w_tot_long
      call run_window(75.0_wp,  t_start, t_window, psi_short, w_tot_short)
      call run_window(900.0_wp, t_start, t_window, psi_long,  w_tot_long)
      call check_true('psi convergence: psi_leaf converged in dt_fast (75 s vs 900 s over 3 h)',                      &
              abs(psi_long - psi_short) < 0.05_wp, abs(psi_long - psi_short))
      !----- the leaf/wood PARTITION is what a charge/credit mismatch moves; the plant TOTAL is set  !
      !      by uptake and transpiration alone and must stay dt-insensitive.                         !
      call check_true('psi convergence: plant TOTAL water dt-insensitive (set by uptake and transpiration)',          &
              abs(w_tot_long - w_tot_short) / max(w_tot_short, tiny(1.0_wp)) < 0.05_wp,                               &
              abs(w_tot_long - w_tot_short) / max(w_tot_short, tiny(1.0_wp)))
   end subroutine test_psi_dt_convergence

   subroutine run_window(dt, t0_sec, t_len, psi_out, w_tot_out)
      real(wp), intent(in)  :: dt, t0_sec, t_len
      real(wp), intent(out) :: psi_out, w_tot_out
      real(wp)    :: t
      integer(ik) :: istep
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      t = 0.0_wp ; istep = 0_ik
      do while (t < t_len - 1.0e-6_wp)
         istep = istep + 1_ik
         !----- set_diurnal_forcing indexes by STEP, so drive it from absolute time instead: the two !
         !      cadences must see the SAME forcing or this measures the forcing, not the scheme. ----!
         call set_forcing_at(t0_sec + t + 0.5_wp * dt)
         call column_fast_step(dt, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
         t = t + dt
      end do
      psi_out = psi_from_water_content(biophys%leaf_water_mass(1), col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1))
      w_tot_out = sum(biophys%leaf_water_mass(1:n) + biophys%wood_water_mass(1:n))
   end subroutine run_window

   subroutine set_forcing_at(t_sec)
      real(wp), intent(in) :: t_sec
      real(wp) :: cosz, t_air
      cosz  = max(0.0_wp, solar_cosz(sim_date, t_sec, lat))
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
      forc%abs_sw = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! MO reference = this air, not 298.15 K
   end subroutine set_forcing_at


   subroutine reset_state()
      integer(ik) :: kk
      if (allocated(biophys%leaf_temp)) deallocate(biophys%leaf_temp)
      call alloc_patch_biophys(biophys, n, t0, 0.008_wp, 400.0_wp, t0)
      !----- alloc_patch_biophys starts leaf_water_mass/wood_water_mass at a scratch 0, which the    !
      !      production gather overwrites from the cohort block -- give them a physical store at a  !
      !      mild tension, water_content(PSI_FIXTURE,...), or psi_from_water_content would          !
      !      diagnose an unphysical psi from an empty pool. -----------------------------------------!
      biophys%leaf_water_mass(1:n) = water_content(PSI_FIXTURE, col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1:n))
      biophys%wood_water_mass(1:n) = water_content(PSI_FIXTURE, col_config%hydraulics_table%pft(1)%wood_curve, &
            col_cohort%bsap(1:n) + col_cohort%broot(1:n))
      budget = column_budget_t()
      biophys%soil_water%theta(1:nsl) = theta_seed
      do kk = 1_ik, nsl
         biophys%soil_energy%soil_energy(kk) = temp_to_internal_energy(col_config%soil_thermal%soil_dry_heat_capacity(kk), &
                                      theta_seed * rho_h2o, t0, 1.0_wp)
         biophys%soil_energy%soil_temp(kk)   = t0
      end do
   end subroutine reset_state

   subroutine set_noon_forcing()
      forc%abs_sw = 450.0_wp ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 70.0_wp ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp                                    ! no rain
      forc%enthalpy_atm = cas_enthalpy_of_temp(295.0_wp, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
   end subroutine set_noon_forcing

   subroutine set_diurnal_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
      forc%abs_sw = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp                                    ! dry window (no rainfall)
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! MO reference = this air, not 298.15 K
   end subroutine set_diurnal_forcing


   !----- The same 24 h with no rain and no shed water, from the state the caller set up: the       !
   !      column's summed water content at the end. Free drainage empties a wet column whatever      !
   !      falls on it, so an input wets the soil when the column ends wetter than this day does. --!
   subroutine march_dry_day(theta_sum)
      real(wp), intent(out) :: theta_sum
      integer(ik) :: istep
      biophys%shed_water_rate = 0.0_wp
      do istep = 1_ik, 24_ik * steps_per_hour
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, gpp_coh=gpp_coh)
      end do
      theta_sum = sum(biophys%soil_water%theta(1:nsl))
   end subroutine march_dry_day

end program test_column_ark
