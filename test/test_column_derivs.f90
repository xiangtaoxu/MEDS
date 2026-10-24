! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_column_derivs -- unit tests for the pure fast-loop right-hand side f(y) = dy/dt               !
! (surface_derivs, march_derivs) and the steppers built on it (docs/science/numerical_scheme.md     !
! §3). Both integrators evaluate this side-effect-free RHS, so these tests pin it, starting with:    !
!   1. LEAF DIAGNOSTIC CLOSURE   -- the diagnosed leaf temperature zeroes the linearized leaf        !
!                                    energy balance it was solved from (multi-cohort).               !
!   2. ANALYTIC LEAF TEMPERATURE -- a hand-computable no-latent / no-LW-forcing case matches.        !
!   3. CAS BACKWARD-EULER CONSISTENCY -- the CAS box's committed enth1/shv1/co21 (the BE-in-atm      !
!                                    solution of the tendencies) satisfy the implicit relation        !
!                                    (y1-y0)/dt = f(y1); i.e. d_cas_* is the correct RHS.             !
!   4. CONSERVATION MARCH        -- marching the CAS twins with the RHS closes the energy/water        !
!                                    budgets to round-off and warms the CAS toward a warm atmosphere. !
!==========================================================================================!
program test_column_derivs
   use meds_test_assert, only : check, check_true, check_close, test_report
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : latent_heat_vap, stefan_boltzmann, cp_air, tiny_num, rho_h2o, n_soil_layer_max, &
                                   h2o_dry_ratio
   use meds_therm_lib,         only : cas_temp_to_enthalpy, cas_enthalpy_to_temp, sat_specific_humidity,  &
                                      sat_vapor_pressure, sat_vapor_pressure_temp_deriv, enthalpy_vapor,  &
                                      temp_to_internal_energy, internal_energy_liquid
   use meds_budget_check,   only : budget_t, budget_accumulate
   use meds_soil_heat, only : soil_heat_env_t, soil_heat_out_t
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_params,      only : soil_profile_t, soil_params_t
   use meds_soil_derive,        only : build_profile_hydraulics, build_profile_thermal
   use meds_patch_state_vars,   only : soil_energy_state_t
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG
   use meds_soil_heat,      only : soil_energy_step_implicit, soil_energy_time_deriv, soil_layer_temp
   use meds_canopy_radiation, only : radiative_temperature
   use meds_pft_params,  only : hydraulics_pft_t
   use meds_plant_hydr_lib,   only : tissue_psi_to_water, water_curve_t
   use meds_fast_integ_derivs, only : surface_derivs, march_derivs
   use meds_fast_integ_types,       only : surface_state_t, surface_tend_t, march_state_t, frozen_t,            &
                                     march_tend_t
   use meds_site_state_vars,  only : site_surface_t
   use meds_fast_integ_backward_euler,    only : backward_euler_stage, advance_plant_water
   use meds_fast_integ_esdirk2,         only : esdirk2_step
   use meds_fast_integ_stepper,        only : adaptive_march
   use meds_fast_integ_rk45,        only : rk45_step
   use meds_fast_integ_control,     only : default_error_control
   use meds_fast_integ_types,       only : error_control_t, fast_workspace_t, state_init
   use meds_fast_params,      only : CTRL_PI, INTEG_ESDIRK2
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks were written on
   !----- The site every make_column fixture runs on (the frozen record holds only what freezes     !
   !      each step; the site's soil travels separately, as in the model). ------------------------!
   type(site_surface_t)  :: site
   type(error_control_t) :: ec_default   !< the defaults: no budget stop
   type(fast_workspace_t) :: ws           !< the fast step's storage, reused by every test as in the model

   call test_leaf_closure()
   call test_leaf_analytic()
   call test_falling_leaf_film()
   call test_cas_be_consistency()
   call test_conservation_march()
   call test_soil_energy_tendency()
   call test_column_assembler()
   call test_be_euler()
   call test_be_coupled()
   call test_arrowhead()
   call test_adaptive_march()
   call test_esdirk2()
   call test_rk45_order()

   !----- 16. CONSTITUTIVE-DOMAIN SAFETY of the explicit RHS. ---------------------------------------!
   call test_rhs_domain_safety()

   call test_report('test_column_derivs')

contains

   !----- A fixed root profile for the soil-only tests, e^(-2 d) over each layer, normalized: the  !
   !      roots belong to the cohorts in the model (cohort_root_profile); a soil test needs only a  !
   !      sink that sums to one. --------------------------------------------------------------------!
   pure function exp_root_profile(soil, nsl) result(frac)
      type(soil_profile_t), intent(in) :: soil
      integer(ik),          intent(in) :: nsl
      real(wp) :: frac(nsl)
      integer(ik) :: k
      do k = 1_ik, nsl
         frac(k) = exp(2.0_wp * soil%soil_layer_z(k)) - exp(2.0_wp * soil%soil_layer_z(k+1_ik))
      end do
      frac = frac / sum(frac)
   end function exp_root_profile



   !----- A representative 3-cohort daytime surface setup (frozen pre-pass + aerodynamics). -----!
   subroutine make_frozen(frozen, n)
      type(frozen_t),  intent(out) :: frozen
      integer(ik),            intent(in)  :: n
      integer(ik) :: i
      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n), frozen%tissue%abs_lw(n), &
               frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      !----- Tissue store: a = cap/dt_fast, relaxing from t_*0. ZERO capacity here keeps every       !
      !      assertion in this file on the zero-inertia balance it was written against. ------------!
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      frozen%tissue%leaf_hcap_per_dt = 0.0_wp ; frozen%tissue%wood_hcap_per_dt = 0.0_wp
      frozen%tissue%t_leaf0 = 298.5_wp ; frozen%tissue%t_wood0 = 298.5_wp   ! the emission base, t_rad
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n), frozen%film%transfer_lw(n))
      frozen%film%transfer_lw = 0.0_wp
      frozen%film%intercept_leaf = 0.0_wp ; frozen%film%intercept_wood = 0.0_wp
      frozen%tissue%h_coeff_w = 0.0_wp ; frozen%tissue%abs_sw_wood = 0.0_wp
      frozen%tissue%abs_lw_wood = 0.0_wp ; frozen%tissue%wai = 0.0_wp
      frozen%tissue%qwflux_wl = 0.0_wp ; frozen%tissue%q_wood_net = 0.0_wp   ! advective enthalpy: no-op unless populated
      frozen%film%f_wet_c = 0.0_wp ; frozen%film%g_film_leaf = 0.0_wp
      frozen%film%g_film_w = 0.0_wp   ! canopy water: no-op unless set
      do i = 1_ik, n
         frozen%tissue%lai(i)       = 2.0_wp - 0.4_wp * real(i - 1_ik, wp)          ! 2.0, 1.6, 1.2
         frozen%tissue%abs_sw(i)    = 250.0_wp - 40.0_wp * real(i - 1_ik, wp)       ! more light at the top
         frozen%tissue%abs_lw(i)    = -30.0_wp
         frozen%tissue%h_coeff_leaf(i) = 2.0_wp * frozen%tissue%lai(i) * 0.03_wp * 1.2_wp * cp_air  ! effarea*lai*gbh*rho*cp
         frozen%tissue%g_transp_leaf(i)    = 0.004_wp * frozen%tissue%lai(i)                          ! series conductance [m/s]
      end do
      allocate(frozen%tissue%leaf_emiss(n), frozen%tissue%wood_emiss(n))
      frozen%tissue%leaf_emiss    = 0.95_wp ; frozen%tissue%wood_emiss = 0.95_wp
      frozen%cas%rho           = 1.2_wp
      frozen%cas%press         = 101325.0_wp
      associate (b => frozen%cas%box)
         b%air_mass_capacity        = frozen%cas%rho * 20.0_wp                                 ! rho * can_depth
         b%air_molar_capacity       = (frozen%cas%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 20.0_wp
         b%atm_conductance_enthalpy = frozen%cas%rho * 0.3_wp * 0.02_wp                        ! rho * ustar * heat_transfer_factor
         b%atm_conductance_vapor    = frozen%cas%rho * 0.3_wp * 0.02_wp
         b%atm_conductance_co2      = (frozen%cas%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 0.3_wp * 0.02_wp
         b%atm_enthalpy             = cas_temp_to_enthalpy(300.0_wp, 0.011_wp)
         b%atm_specific_humidity    = 0.011_wp
         b%atm_co2                  = 400.0_wp
      end associate
      frozen%cas%nee_biotic    = -5.0_wp                                          ! net CO2 uptake [umol/m2/s]
      frozen%ground%rnet_soil = 50.0_wp
      frozen%ground%ggnet     = 0.02_wp
      frozen%ground%soil_evap = 2.0e-5_wp
   end subroutine make_frozen

   !----- 1. The diagnosed leaf temperature must zero the (linearized) leaf energy balance. ------!
   subroutine test_leaf_closure()
      type(frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat, dtl, lw_slope, le_slope, le_ref, resid, worst, t_rad
      integer(ik) :: i, n
      n = 3_ik
      print '(a)', 'test_leaf_closure:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_temp_to_enthalpy(298.0_wp, 0.012_wp)
      y%cas_shv      = 0.012_wp
      y%cas_co2      = 410.0_wp
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, n, f)

      tcas   = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, frozen%cas%press)
      esat   = sat_vapor_pressure(tcas)
      dqdt   = h2o_dry_ratio * frozen%cas%press                                                                  &
               / max((frozen%cas%press - (1.0_wp - h2o_dry_ratio) * esat) ** 2, tiny_num)                       &
               * sat_vapor_pressure_temp_deriv(tcas)
      worst = 0.0_wp
      do i = 1_ik, n
         dtl      = f%leaf_temp(i) - tcas
         !----- The longwave about the emission the two-stream's abs_lw is net of, at t_rad. ---------!
         t_rad    = radiative_temperature(frozen%tissue%lai(i), frozen%tissue%wai(i), frozen%tissue%t_leaf0(i),     &
                                          frozen%tissue%t_wood0(i))
         lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(i) * stefan_boltzmann * t_rad ** 3 * frozen%tissue%lai(i)
         !----- the leaf pays the FULL vapour enthalpy at the linearization temperature             !
         !      (docs/science/vegetation_energy_dynamics.md eq. (2)). ------------------------------!
         le_slope = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * dqdt
         le_ref   = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * (qsat_c - qcas)
         !----- Rnet - sensible - latent - LW-emission, all at the diagnosed leaf temperature. ---!
         resid = frozen%tissue%abs_sw(i) + frozen%tissue%abs_lw(i) - frozen%tissue%h_coeff_leaf(i) * dtl                          &
                 - (le_ref + le_slope * dtl) - lw_slope * (f%leaf_temp(i) - t_rad)
         worst = max(worst, abs(resid))
      end do
      call check_true('diagnosed leaf T zeroes the linearized balance', worst < 1.0e-9_wp, worst)
      call check_true('leaf temperatures physical (270-330 K)',                                 &
                      minval(f%leaf_temp) > 270.0_wp .and. maxval(f%leaf_temp) < 330.0_wp,       &
                      f%leaf_temp(1))
   end subroutine test_leaf_closure

   !----- 2. No latent (g_transp_leaf = 0) and no LW forcing (abs_lw = 0): the leaf balances its       !
   !      shortwave against sensible heat and its emission about t_rad, dtl = (abs_sw +                !
   !      lw_slope*(t_rad - tcas)) / (h + lw_slope). --------------------------------------------------!
   subroutine test_leaf_analytic()
      type(frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp) :: tcas, lw_slope, expect, t_rad
      print '(a)', 'test_leaf_analytic:'
      call make_frozen(frozen, 1_ik)
      frozen%tissue%g_transp_leaf(1) = 0.0_wp ; frozen%tissue%abs_lw(1) = 0.0_wp ; frozen%tissue%abs_sw(1) = 200.0_wp
      frozen%tissue%lai(1) = 1.5_wp
      frozen%tissue%h_coeff_leaf(1) = 2.0_wp * frozen%tissue%lai(1) * 0.03_wp * 1.2_wp * cp_air
      y%cas_enthalpy = cas_temp_to_enthalpy(299.0_wp, 0.010_wp)
      y%cas_shv      = 0.010_wp
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, 1_ik, f)
      tcas     = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      t_rad    = radiative_temperature(frozen%tissue%lai(1), frozen%tissue%wai(1), frozen%tissue%t_leaf0(1),       &
                                       frozen%tissue%t_wood0(1))
      lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(1) * stefan_boltzmann * t_rad ** 3 * frozen%tissue%lai(1)
      expect   = tcas + (frozen%tissue%abs_sw(1) + lw_slope * (t_rad - tcas)) / (frozen%tissue%h_coeff_leaf(1) + lw_slope)
      call check('leaf T = tcas + (Rn + lw_slope (t_rad - tcas))/(h + lw_slope) (no latent, no LW)', f%leaf_temp(1), &
                 expect, 1.0e-10_wp)
   end subroutine test_leaf_analytic

   !----- 2b. A falling leaf's film passed to its wood (Palo Verde's dry season): the leaf, nearly    !
   !      gone, hands its whole film to the wood as the cohort's film is re-shared by capacity. The   !
   !      giver is charged nothing for the water it passes, so the leaf balances its sensible heat    !
   !      against its emission about t_rad alone; charged the rain's valuation, it sat 4e10 K off. ---!
   subroutine test_falling_leaf_film()
      type(frozen_t)        :: frozen
      type(surface_state_t) :: y
      type(surface_tend_t)  :: f
      real(wp) :: tcas, t_rad, lw_slope, expect
      print '(a)', 'test_falling_leaf_film:'
      call make_frozen(frozen, 1_ik)
      frozen%tissue%lai(1) = 1.0e-20_wp
      frozen%tissue%abs_sw(1) = 0.0_wp ; frozen%tissue%abs_lw(1) = 0.0_wp
      frozen%tissue%h_coeff_leaf(1)  = 2.0_wp * frozen%tissue%lai(1) * 0.03_wp * 1.2_wp * cp_air
      frozen%tissue%g_transp_leaf(1) = 0.004_wp * frozen%tissue%lai(1)
      frozen%tissue%wai(1) = 0.5_wp ; frozen%tissue%h_coeff_w(1) = 10.0_wp
      frozen%tissue%wood_hcap_per_dt(1) = 1.0e4_wp / 900.0_wp
      frozen%tissue%t_leaf0(1) = 303.9_wp ; frozen%tissue%t_wood0(1) = 304.5_wp
      frozen%film%f_wet_c(1) = 1.0_wp
      frozen%film%intercept_leaf(1) = -1.07e-12_wp            ! the leaf's whole film, passed on
      frozen%film%transfer_lw(1)    =  1.07e-12_wp
      frozen%film%intercept_wood(1) =  1.07e-12_wp
      frozen%film%rain_liquid_enthalpy = internal_energy_liquid(297.0_wp)
      y%cas_enthalpy = cas_temp_to_enthalpy(299.0_wp, 0.010_wp)
      y%cas_shv      = 0.010_wp
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, 1_ik, f)
      tcas     = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      t_rad    = radiative_temperature(frozen%tissue%lai(1), frozen%tissue%wai(1), frozen%tissue%t_leaf0(1),       &
                                       frozen%tissue%t_wood0(1))
      lw_slope = 4.0_wp * frozen%tissue%leaf_emiss(1) * stefan_boltzmann * t_rad ** 3 * frozen%tissue%lai(1)
      expect   = tcas + lw_slope * (t_rad - tcas) / (frozen%tissue%h_coeff_leaf(1) + lw_slope)
      call check_true('a falling leaf passing its film to its wood is charged nothing for it',                   &
                      abs(f%leaf_temp(1) - expect) < 1.0e-8_wp, f%leaf_temp(1) - expect)
   end subroutine test_falling_leaf_film

   !----- 3. The CAS box commits the BE-in-atm solution of the tendencies; verify (y1-y0)/dt = f(y1). !
   subroutine test_cas_be_consistency()
      type(frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: dt, enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: r_enth, r_shv, r_co2
      integer(ik) :: n
      n = 3_ik ; dt = 900.0_wp
      print '(a)', 'test_cas_be_consistency:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_temp_to_enthalpy(296.0_wp, 0.013_wp)
      y%cas_shv      = 0.013_wp
      y%cas_co2      = 415.0_wp
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, n, f)
      !----- Reconstruct the implicit CAS update the box kernel commits (cas_box_step_implicit). !
      associate (b => frozen%cas%box)
         enth1 = (b%air_mass_capacity * enth0 + dt * (f%src_enth + b%atm_conductance_enthalpy * b%atm_enthalpy))   &
                 / (b%air_mass_capacity + dt * b%atm_conductance_enthalpy)
         shv1  = (b%air_mass_capacity * shv0 + dt * (f%src_vap + b%atm_conductance_vapor * b%atm_specific_humidity)) &
                 / (b%air_mass_capacity + dt * b%atm_conductance_vapor)
         co21  = (b%air_molar_capacity * co20 + dt * (frozen%cas%nee_biotic + b%atm_conductance_co2 * b%atm_co2))   &
                 / (b%air_molar_capacity + dt * b%atm_conductance_co2)
         !----- BE consistency: (y1-y0)/dt must equal the tendency evaluated with the ATM term at y1  !
         !      (source frozen). This is exactly what an implicit stage solves, so it verifies        !
         !      d_cas_* is the correct RHS of the box's implicit update                               !
         !      (docs/science/canopy_air_space_biophysics.md eq. (2)).                                !
         r_enth = (enth1 - enth0) / dt - (f%src_enth + b%atm_conductance_enthalpy * (b%atm_enthalpy - enth1))      &
                  / b%air_mass_capacity
         r_shv  = (shv1  - shv0 ) / dt - (f%src_vap  + b%atm_conductance_vapor * (b%atm_specific_humidity - shv1)) &
                  / b%air_mass_capacity
         r_co2  = (co21  - co20 ) / dt - (frozen%cas%nee_biotic + b%atm_conductance_co2 * (b%atm_co2 - co21))      &
                  / b%air_molar_capacity
         call check_true('CAS enthalpy BE-consistent with d_cas_enthalpy', abs(r_enth) < 1.0e-9_wp, r_enth)
         call check_true('CAS humidity BE-consistent with d_cas_shv',      abs(r_shv)  < 1.0e-15_wp, r_shv)
         call check_true('CAS CO2 BE-consistent with d_cas_co2',           abs(r_co2)  < 1.0e-9_wp, r_co2)
         !----- At t=0 the tendency must equal (src + g*(atm - y0))/cap (sanity on the returned RHS). !
         call check('d_cas_enthalpy = (src+g_atm_heat*(atm-y0))/cas_mass_capacity',                               &
                    f%d_cas_enthalpy, (f%src_enth + b%atm_conductance_enthalpy * (b%atm_enthalpy - enth0))     &
                                      / b%air_mass_capacity, 1.0e-9_wp)
      end associate
   end subroutine test_cas_be_consistency

   !----- 4. March the CAS twins with the RHS (BE-in-atm); the closed budgets must stay tight. ---!
   subroutine test_conservation_march()
      type(frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      type(budget_t)         :: be, bw
      real(wp)    :: dt, enth0, shv0, enth1, shv1, t0, t1
      integer(ik) :: step, n
      n = 3_ik ; dt = 120.0_wp
      print '(a)', 'test_conservation_march:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_temp_to_enthalpy(294.0_wp, 0.010_wp)     ! CAS starts cool + dry
      y%cas_shv      = 0.010_wp
      y%cas_co2      = 400.0_wp
      t0 = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      do step = 1_ik, 60_ik
         call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, 297.0_wp, n, f)
         enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
         associate (b => frozen%cas%box)
            enth1 = (b%air_mass_capacity * enth0 + dt * (f%src_enth + b%atm_conductance_enthalpy * b%atm_enthalpy))   &
                    / (b%air_mass_capacity + dt * b%atm_conductance_enthalpy)
            shv1  = (b%air_mass_capacity * shv0 + dt * (f%src_vap + b%atm_conductance_vapor * b%atm_specific_humidity)) &
                    / (b%air_mass_capacity + dt * b%atm_conductance_vapor)
            !----- The same closed-budget accounting the CAS box kernel performs. -------------------!
            call budget_accumulate(be, b%air_mass_capacity * enth0, b%air_mass_capacity * enth1,                     &
                                   f%src_enth + b%atm_conductance_enthalpy * b%atm_enthalpy,                         &
                                   b%atm_conductance_enthalpy * enth1, dt, abs(b%air_mass_capacity * enth1),          &
                                   1.0e-8_wp, 1.0e-3_wp)
            call budget_accumulate(bw, b%air_mass_capacity * shv0, b%air_mass_capacity * shv1,                       &
                                   f%src_vap + b%atm_conductance_vapor * b%atm_specific_humidity,                    &
                                   b%atm_conductance_vapor * shv1, dt, max(abs(b%air_mass_capacity * shv1), 1.0e-6_wp), &
                                   1.0e-8_wp, 1.0e-10_wp)
         end associate
         y%cas_enthalpy = enth1 ; y%cas_shv = shv1
      end do
      t1 = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      call check_true('CAS energy budget closes over the march', be%n_fail == 0_ik, real(be%worst, wp))
      call check_true('CAS water  budget closes over the march', bw%n_fail == 0_ik, real(bw%worst, wp))
      call check_true('CAS warms toward the warm atmosphere + sunlit canopy', t1 > t0, t1 - t0)
   end subroutine test_conservation_march


   !----- 6. soil_energy_time_deriv == the dt->0 limit of the BE kernel soil_energy_step_implicit. ----------!
   subroutine test_soil_energy_tendency()
      type(soil_profile_t)        :: soil
      type(soil_heat_env_t)      :: forcing
      type(soil_params_t)         :: eopts
      type(soil_energy_state_t)   :: col, col0
      type(soil_heat_out_t)         :: eflux
      real(wp)    :: dedt(n_soil_layer_max), dt, worst, fd
      integer(ik) :: k, nsl
      nsl = 10_ik
      print '(a)', 'test_soil_energy_tendency:'
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,      &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, soil)
      call build_profile_thermal(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, soil)
      forcing%soil_water(1:10) = 0.30_wp ; forcing%w_flux = 0.0_wp
      forcing%g_top = 120.0_wp ; forcing%geothermal = 0.0_wp
      forcing%water_heat_out(1:10) = 3.0_wp                       ! nonzero interior sink
      do k = 1_ik, 10_ik                                          ! a temperature gradient -> nonzero faces
         col%soil_energy(k) = temp_to_internal_energy(soil%dry_heat_capacity(k),                      &
                              forcing%soil_water(k) * rho_h2o, 291.0_wp - 0.6_wp * real(k-1_ik, wp), 1.0_wp)
      end do
      col0 = col
      call soil_energy_time_deriv(col0, forcing, soil, eopts, dedt)
      dt = 1.0e-2_wp
      call soil_energy_step_implicit(col, forcing, soil, eopts, dt, eflux)   ! advances col (BE)
      worst = 0.0_wp
      do k = 1_ik, nsl
         fd    = (col%soil_energy(k) - col0%soil_energy(k)) / dt
         worst = max(worst, abs(fd - dedt(k)) / max(abs(dedt(k)), 1.0_wp))
      end do
      call check_true('soil-energy tendency = BE-limit of soil_energy_step_implicit (O(dt))', worst < 1.0e-2_wp, worst)
   end subroutine test_soil_energy_tendency


   !----- 9. march_derivs assembles a finite, physically-signed whole-column RHS + its CAS part   !
   !         matches a standalone surface_derivs call (correct wiring). -------------------------!
   subroutine test_column_assembler()
      type(march_state_t)       :: y
      type(frozen_t)      :: frozen
      type(march_tend_t)        :: f
      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_heat_env_t)     :: eforc_chk
      type(soil_energy_state_t) :: se_chk
      real(wp)    :: dedt_chk(n_soil_layer_max), worst
      integer(ik) :: k, i, n, nsl
      logical     :: all_finite
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_column_assembler:'
      call make_column(y, frozen, n, nsl)
      !----- The water network's flows, nonzero, so the wiring check below is live for every term: a  !
      !      downward drift through every face, drainage, infiltration, the uptake's enthalpy leaving   !
      !      each layer, and water pushed up into the pond from the top layer. -----------------------!
      frozen%hydrology%w_flux_frozen(1:nsl) = 2.0e-8_wp ; frozen%hydrology%w_flux_frozen(nsl) = 0.0_wp
      frozen%hydrology%drainage     = 2.0e-8_wp * rho_h2o
      frozen%hydrology%infiltration = frozen%hydrology%q_top * rho_h2o
      frozen%hydrology%exfil_enth   = 1.0e-7_wp * internal_energy_liquid(296.0_wp)
      do k = 1_ik, nsl
         frozen%roots%root_heat(k) = frozen%roots%uptake * frozen%roots%root_share(k) * internal_energy_liquid(295.0_wp)
      end do
      call march_derivs(y, frozen, site, n, nsl, f)

      all_finite = ieee_ok(f%d_cas_enthalpy) .and. ieee_ok(f%d_cas_shv) .and. ieee_ok(f%d_cas_co2)
      do k = 1_ik, nsl ; all_finite = all_finite .and. ieee_ok(f%dedt(k)) .and. ieee_ok(f%dtheta_dt(k)) ; end do
      do i = 1_ik, n
         all_finite = all_finite .and. ieee_ok(f%d_leaf_water_mass(i)) .and. ieee_ok(f%d_wood_water_mass(i))
      end do
      call check_true('march_derivs: all tendencies finite', all_finite, 0.0_wp)

      !----- the CAS part must equal a standalone surface_derivs at the state's diagnosed t_ground. !
      y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground,                  &
                          tground_of(y, frozen), n, surf_tend)
      call check('march_derivs CAS enthalpy tendency = surface_derivs', f%d_cas_enthalpy, surf_tend%d_cas_enthalpy, 1.0e-12_wp)

      !----- The CANOPY is energy-neutral. With no advected enthalpy (qwflux_wl = q_wood_net = 0 in     !
      !      this fixture) the tissues hold no store on this path, so what the canopy absorbs           !
      !      (coh_rnet) must equal what it hands the CAS: sensible + the FULL enthalpy of the vapour it !
      !      sheds. Guards against handing the CAS transp*enthalpy_vapor while the leaf pays only       !
      !      latent_heat_vap and charging the difference to the soil: that breaks this identity by the  !
      !      proxy (~30 W/m2 at 3 mm/day) and makes the soil pay the transpired water's liquid enthalpy !
      !      twice, once through the proxy and once through qloss. -------------------------------------!
      block
         real(wp) :: canopy_to_cas, tcas_chk
         tcas_chk = cas_enthalpy_to_temp(y_stage%cas_enthalpy, y_stage%cas_shv)
         canopy_to_cas = surf_tend%src_enth - surf_tend%h_ground - surf_tend%le_ground + surf_tend%cond * &
               internal_energy_liquid(tcas_chk)
         call check_true('canopy transpires at all in this fixture (test is live)',                          &
               sum(surf_tend%transp_c(1:n)) > 1.0e-7_wp, sum(surf_tend%transp_c(1:n)))
         call check('canopy is energy-neutral: coh_rnet == sensible + full vapour enthalpy to the CAS (no soil proxy)', &
                    canopy_to_cas, surf_tend%coh_rnet, 1.0e-9_wp * max(abs(surf_tend%coh_rnet), 1.0_wp))
      end block

      !----- wiring check: the assembled soil-heat tendency == a standalone soil_energy_time_deriv    !
      !      built from the SAME surface coupling (g_top) and the water network's flows: the per-layer  !
      !      enthalpy leaving with uptake, the time-mean faces, the infiltration at the pond's          !
      !      temperature and the drainage enthalpy at the bottom. Every term is nonzero in this         !
      !      fixture: a check that omits a term is only green while that term is zero. -----------------!
      se_chk%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc_chk%g_top = surf_tend%g_top
      do k = 1_ik, nsl
         eforc_chk%soil_water(k)     = y%theta(k)
         eforc_chk%water_heat_out(k) = frozen%roots%root_heat(k)
         eforc_chk%w_flux(k)         = -frozen%hydrology%w_flux_frozen(k)
      end do
      eforc_chk%water_heat_out(1)   = eforc_chk%water_heat_out(1) + frozen%hydrology%exfil_enth
      eforc_chk%water_heat_out(nsl) = eforc_chk%water_heat_out(nsl)                                  &
                                    + frozen%hydrology%drainage * internal_energy_liquid(frozen%hydrology%t_bot)
      eforc_chk%w_flux_top  = -frozen%hydrology%infiltration / rho_h2o
      eforc_chk%t_water_top = frozen%hydrology%t_infil
      call soil_energy_time_deriv(se_chk, eforc_chk, site%soil_profile, site%soil, dedt_chk)
      worst = maxval(abs(f%dedt(1:nsl) - dedt_chk(1:nsl)))
      call check_true('march_derivs wires the soil-heat tendency correctly', worst < 1.0e-9_wp, worst)
   end subroutine test_column_assembler

   !----- 11. First-order BE (be_euler_step): stable at dt = 900 s, where an explicit step blows up,  !
   !          and close to the explicit RK45 reference at small dt (both integrate the same RHS). -----!
   subroutine test_be_euler()
      type(march_state_t)  :: y, yi, yr, ytmp
      type(frozen_t) :: frozen
      real(wp)    :: tcas, dcas, dtheta, dmass
      integer(ik) :: step, n, nsl, k, i
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_be_euler:'
      call make_column(y, frozen, n, nsl)
      !----- Well-ventilated canopy (stronger CAS<->atm exchange) so the frozen-per-step source at   !
      !      the full 900 s does not out-run venting into supersaturation -- the large-dt operator-   !
      !      split coupling error the arrowhead Newton removes (test_arrowhead); here we exercise    !
      !      L-stability + physical.                                                                 !
      frozen%cas%box%atm_conductance_enthalpy = frozen%cas%box%atm_conductance_enthalpy * 4.0_wp
      frozen%cas%box%atm_conductance_vapor    = frozen%cas%box%atm_conductance_vapor * 4.0_wp
      frozen%cas%box%atm_conductance_co2      = frozen%cas%box%atm_conductance_co2 * 4.0_wp

      !----- (a) STABLE + physical at the full fast timestep dt = 900 s over a 6-hour march. --------!
      call copy_state(y, yi, n)
      physical = .true.
      do step = 1_ik, 24_ik                    ! 24 * 900 s = 6 h
         call be_euler_step(yi, frozen, n, nsl, 900.0_wp, ytmp)
         call copy_state(ytmp, yi, n)
         tcas = cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. yi%cas_shv > 0.0_wp
         do k = 1_ik, nsl
            physical = physical .and. yi%theta(k) > 0.0_wp .and. yi%theta(k) < 0.6_wp
         end do
         do i = 1_ik, n
            physical = physical .and. yi%leaf_water_mass(i) > 0.0_wp
         end do
      end do
      call check_true('first-order BE stable + physical at dt = 900 s (6 h)', physical, tcas)

      !----- (b) cross-validation: the BE step against the explicit RK45 reference at dt = 4 s over      !
      !          4 min agrees to first order. RK45 shares only the right-hand side with the BE stage, so   !
      !          agreement checks the implicit machinery against it. make_column leaves theta_rate at 0,  !
      !          so theta stays put in both; plant water moves at constant rates in both (in the RK45     !
      !          stages, and through advance_plant_water), a small O(dt) difference. ---------------------!
      call copy_state(y, yi, n)
      do step = 1_ik, 60_ik                    ! 60 * 4 s = 4 min
         call be_euler_step(yi, frozen, n, nsl, 4.0_wp, ytmp)
         call copy_state(ytmp, yi, n)
      end do
      call march_rk45(y, frozen, n, nsl, 4.0_wp, 60_ik, yr)
      dcas   = abs(cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      dtheta = maxval(abs(yi%theta(1:nsl) - yr%theta(1:nsl)))
      dmass  = max(maxval(abs(yi%leaf_water_mass(1:n) - yr%leaf_water_mass(1:n))),                &
                  maxval(abs(yi%wood_water_mass(1:n) - yr%wood_water_mass(1:n))))
      call check_true('first-order BE ~ RK45 reference: CAS temp agree < 5e-2 K at dt=4 s', dcas < 5.0e-2_wp, dcas)
      call check_true('first-order BE ~ RK45 reference: theta unchanged in both (theta_rate = 0)', dtheta < 1.0e-12_wp, &
                      dtheta)
      call check_true('first-order BE ~ RK45 reference: mass agree < 1e-4 kg/plant', dmass < 1.0e-4_wp, dmass)
   end subroutine test_be_euler

   !----- 12. the leaf<->CAS Newton coupling keeps the canopy air physical and sub-saturated under a   !
   !          harsh constant-noon forcing at 900 s, where an uncoupled pass over-humidifies.          !
   subroutine test_be_coupled()
      type(march_state_t)  :: y, yc, ytmp
      type(frozen_t) :: frozen
      real(wp)    :: tcas_coup, qsat_coup
      integer(ik) :: step, n, nsl
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_be_coupled:'
      !----- HARSH forcing: constant noon, modest venting (the case that over-humidifies at 900 s). -!
      call make_column(y, frozen, n, nsl)
      !----- The leaves transpire at the canopy air's humidity here, as the stages do whenever the    !
      !      transpiration is not held at the water network's: that leaf-vapour feedback is the        !
      !      coupling the Newton resolves. -------------------------------------------------------!
      frozen%tissue%transp_held = .false.
      call copy_state(y, yc, n)
      do step = 1_ik, 12_ik                 ! 12 * 900 s = 3 h of constant noon
         call be_euler_step(yc, frozen, n, nsl, 900.0_wp, ytmp)
         call copy_state(ytmp, yc, n)
      end do
      tcas_coup = cas_enthalpy_to_temp(yc%cas_enthalpy, yc%cas_shv)
      qsat_coup = sat_specific_humidity(tcas_coup, frozen%cas%press)
      call check_true('coupled CAS stays physical (270-325 K) at 900 s constant noon',            &
                      tcas_coup > 270.0_wp .and. tcas_coup < 325.0_wp, tcas_coup)
      call check_true('coupled CAS stays sub-saturated (shv <= qsat)', yc%cas_shv <= qsat_coup + 1.0e-4_wp, &
                      yc%cas_shv - qsat_coup)
   end subroutine test_be_coupled

   !----- 12b. the arrowhead Newton surface solve: robust near saturation (clamp + line search), and  !
   !           close to the RK45 reference at small dt (converges to the coupled surface solution). ---!
   subroutine test_arrowhead()
      type(march_state_t)  :: y, yi, yr, ytmp
      type(frozen_t) :: frozen
      real(wp)    :: tcas, qsat, dcas, dtheta, worst_super
      integer(ik) :: step, n, nsl, k
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_arrowhead:'

      !----- (a) START the CAS at 99% saturation under harsh noon: exercises the supersaturation      !
      !          clamp + line search. the Newton must keep it physical + sub-saturated at 900 s. !
      call make_column(y, frozen, n, nsl)
      tcas = 297.0_wp ; qsat = sat_specific_humidity(tcas, frozen%cas%press)
      y%cas_enthalpy = cas_temp_to_enthalpy(tcas, 0.99_wp*qsat) ; y%cas_shv = 0.99_wp*qsat
      call copy_state(y, yi, n)
      physical = .true. ; worst_super = -1.0_wp
      do step = 1_ik, 12_ik
         call be_euler_step(yi, frozen, n, nsl, 900.0_wp, ytmp)
         call copy_state(ytmp, yi, n)
         tcas = cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv)
         qsat = sat_specific_humidity(tcas, frozen%cas%press)
         worst_super = max(worst_super, yi%cas_shv - qsat)          ! must stay <= 0 (sub-saturated)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp
      end do
      call check_true('Newton stays sub-saturated from a 99%-saturated start (900 s)', worst_super <= 1.0e-4_wp, worst_super)
      call check_true('Newton stays physical (270-325 K) near saturation', physical, tcas)

      !----- (b) Newton matches the RK45 reference at small dt (converges to the coupled    !
      !          solution); theta stays put through the fixture's zero theta_rate, as in             !
      !          test_be_euler (b). ----------------------------------------------------------------!
      call make_column(y, frozen, n, nsl)
      call copy_state(y, yi, n)
      do step = 1_ik, 60_ik                     ! 60 * 4 s = 4 min
         call be_euler_step(yi, frozen, n, nsl, 4.0_wp, ytmp)
         call copy_state(ytmp, yi, n)
      end do
      call march_rk45(y, frozen, n, nsl, 4.0_wp, 60_ik, yr)
      dcas   = abs(cas_enthalpy_to_temp(yi%cas_enthalpy, yi%cas_shv) - cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv))
      dtheta = maxval(abs(yi%theta(1:nsl) - yr%theta(1:nsl)))
      call check_true('Newton ~ RK45 reference: CAS temp < 5e-2 K at dt=4 s', dcas < 5.0e-2_wp, dcas)
      call check_true('Newton ~ RK45 reference: theta unchanged in both (theta_rate = 0)', dtheta < 1.0e-12_wp, &
                      dtheta)
   end subroutine test_arrowhead

   !----- 13. the production adaptive controller (adaptive_march): tighter rtol takes more steps     !
   !          and both agree with a fine fixed reference. -------------------------------------------!
   subroutine test_adaptive_march()
      type(march_state_t)  :: y, y1, y2, yr
      type(frozen_t) :: frozen
      real(wp)    :: tc1, tc2, tcr
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, ns1, ns2, nr1, nr2
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_adaptive_march:'
      call make_column(y, frozen, n, nsl)
      frozen%cas%box%atm_conductance_enthalpy   = frozen%cas%box%atm_conductance_enthalpy   * 4.0_wp
      frozen%cas%box%atm_conductance_vapor = frozen%cas%box%atm_conductance_vapor * 4.0_wp
      frozen%cas%box%atm_conductance_co2    = frozen%cas%box%atm_conductance_co2    * 4.0_wp
      !----- start the CAS cool + dry so there is a real transient to resolve adaptively. -----------!
      y%cas_enthalpy = cas_temp_to_enthalpy(290.0_wp, 0.009_wp) ; y%cas_shv = 0.009_wp

      !----- The controller under test is the PRODUCTION one, so the test fails for the reasons a
      !      real run would hit.
      ec = default_error_control(1.0e-3_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, site, n, nsl, 1800.0_wp, ec, 50.0_wp, y1, ws, ns1, nr1)
      ec = default_error_control(1.0e-6_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, site, n, nsl, 1800.0_wp, ec, 50.0_wp, y2, ws, ns2, nr2)
      !----- fine fixed reference: ESDIRK2 at dt = 2 s, 900 steps. ----------------------------------!
      call march_esdirk2(y, frozen, n, nsl, 2.0_wp, 900_ik, yr)
      tc1 = cas_enthalpy_to_temp(y1%cas_enthalpy, y1%cas_shv)
      tc2 = cas_enthalpy_to_temp(y2%cas_enthalpy, y2%cas_shv)
      tcr = cas_enthalpy_to_temp(yr%cas_enthalpy, yr%cas_shv)
      print '(a,i0,a,i0,a)', '   adaptive-ESDIRK2 steps: rtol=1e-3 -> ', ns1, ' , rtol=1e-6 -> ', ns2, ' ; fixed dt=2s -> 900'
      call check_true('adaptive-ESDIRK2 rtol=1e-3 agrees with the fine ESDIRK2 reference (< 0.5 K)',                 &
                      abs(tc1 - tcr) < 0.5_wp, tc1 - tcr)
      call check_true('adaptive-ESDIRK2 rtol=1e-6 agrees with the fine ESDIRK2 reference (< 0.5 K)',                 &
                      abs(tc2 - tcr) < 0.5_wp, tc2 - tcr)
      call check_true('tighter rtol takes more steps (the production controller responds to tol)', ns2 > ns1, real(ns2 - ns1, wp))
      call check_true('adaptive is cheaper than the fixed fine march (steps < 900)', ns1 < 900_ik, real(ns1, wp))
      call check_true('the march made real progress (steps > 0)', ns1 > 0_ik, real(ns1, wp))
   end subroutine test_adaptive_march

   !----- soil-top temperature diagnosed from a march_state_t (for the ESDIRK2 order test). -------!
   function soil_top_temp(y, frozen) result(t)
      type(march_state_t),  intent(in) :: y
      type(frozen_t), intent(in) :: frozen
      real(wp) :: t, fl
      call soil_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                 &
                           site%soil_profile%layer_curve(1), t, fl)
   end function soil_top_temp

   !----- march the ESDIRK2 fixed-step from y for nstep steps of dt (embedded error discarded). -----!
   subroutine march_esdirk2(y0, frozen, n, nsl, dt, nstep, y_out)
      type(march_state_t),  intent(in)  :: y0
      type(frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(march_state_t),  intent(out) :: y_out
      type(march_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call esdirk2_step(y, frozen, site, ec_default, n, nsl, dt, ytmp, yerr, ws%stages)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_esdirk2


   !----- ONE first-order backward-Euler column step: the gamma = 1 DEGRADED CONFIGURATION of the    !
   !      production ESDIRK2, built from the same two kernels esdirk2_step composes -- one         !
   !      backward_euler_stage, then the closed-form mass split over the full dt on that stage's own   !
   !      transp_c (weight 1.0, because at gamma = 1 there is only one stage).                          !
   !                                                                                          !
   !      It is a probe, not an oracle: it returns no reference trajectory, and it shares              !
   !      backward_euler_stage with the scheme it degrades, so it is independent in the tableau only,  !
   !      not in the machinery. What it IS good for is exactly what the tests below use it for --      !
   !      deliberately breaking the production scheme in ONE known way (drop to first order) and       !
   !      watching what fails. A degraded-configuration                                                !
   !      probe belongs in the test that degrades it, not in the shipped source.                       !
   subroutine be_euler_step(y, frozen, n, nsl, dt, y_out)
      type(march_state_t),  intent(in)  :: y
      type(frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(march_state_t),  intent(out) :: y_out
      type(surface_tend_t) :: surf_tend
      call backward_euler_stage(y, frozen, site, ec_default, n, nsl, dt, y_out, surf_tend)
      call advance_plant_water(y, frozen, n, dt, surf_tend%transp_c(1:n), y_out)
   end subroutine be_euler_step

   !----- march the first-order BE step fixed-step (for the ESDIRK2-vs-1st-order comparison). --------!
   subroutine march_be_euler(y0, frozen, n, nsl, dt, nstep, y_out)
      type(march_state_t),  intent(in)  :: y0
      type(frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(march_state_t),  intent(out) :: y_out
      type(march_state_t) :: y, ytmp
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call be_euler_step(y, frozen, n, nsl, dt, ytmp)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_be_euler

   !----- march the Cash-Karp RK45 fixed-step (5th-order commit, embedded 4th discarded): the explicit !
   !      reference the BE checks compare against, and the order-of-accuracy study below. ----------!
   subroutine march_rk45(y0, frozen, n, nsl, dt, nstep, y_out)
      type(march_state_t),  intent(in)  :: y0
      type(frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),               intent(in)  :: dt
      type(march_state_t),  intent(out) :: y_out
      type(march_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call rk45_step(y, frozen, site, n, nsl, dt, ytmp, yerr)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_rk45

   !----- 14. ESDIRK2 (ARS(2,2,2) gamma): 2nd order on the decoupled CO2 twin; at least first order !
   !      on the coupled soil-top temperature; plant water mass stays physical at production dt;    !
   !      embedded estimate bounded; stiff-stable (docs/science/numerical_scheme.md eqs. (1)-(2)).  !
   !      ------------------------------------------------------------------------------------------!
   subroutine test_esdirk2()
      type(march_state_t)  :: y, yref, y1, y2, y4, ynew, yerr, ytmp
      type(frozen_t) :: frozen
      type(error_control_t) :: ec
      real(wp)    :: e1, e2, e4, p_lo, p_hi, tref, errnorm, tcas
      integer(ik) :: n, nsl, step, ns, nr, ns2, nr2
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_esdirk2:'
      call make_column(y, frozen, n, nsl)

      !----- (a) ORDER-2 self-convergence on a CLEANLY-integrated tableau variable: CAS CO2 (a         !
      !          decoupled affine scalar ODE -- no operator split, so the ARS(2,2,2) order is exposed).!
      !          reference = ESDIRK2 at dt = 1800/64; solutions at dt = 450, 225, 112.5. -------------!
      call march_esdirk2(y, frozen, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = yref%cas_co2
      call march_esdirk2(y, frozen, n, nsl, 450.0_wp,   4_ik,  y1) ; e1 = abs(y1%cas_co2 - tref)
      call march_esdirk2(y, frozen, n, nsl, 225.0_wp,   8_ik,  y2) ; e2 = abs(y2%cas_co2 - tref)
      call march_esdirk2(y, frozen, n, nsl, 112.5_wp,  16_ik,  y4) ; e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)', '   CAS CO2 errors: ', e1, ' / ', e2, &
            ' / ', e4, ' ; observed order p = ', p_lo, ' , ', p_hi
      call check_true('ESDIRK2 tableau is 2nd order on the decoupled CO2 twin (p >= 1.9)', &
                      p_lo >= 1.9_wp .and. p_hi >= 1.9_wp, p_hi)

      !----- (a2) On the operator-split-COUPLED soil-top temperature (soil water and plant water mass  !
      !      are split out of the ESDIRK stages; see docs/ROADMAP.md) the effective order is about 1:  !
      !      ESDIRK2's error against the fine reference falls at least first order with dt (450, 225,  !
      !      112.5 s). The first-order BE step's error is printed beside it: its constant is the       !
      !      fixture's and can nearly cancel (it does with the longwave about the tissue's emission    !
      !      temperature), so it is no bound. ---------------------------------------------------------!
      call march_esdirk2(y, frozen, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = soil_top_temp(yref, frozen)
      call march_esdirk2(y, frozen, n, nsl, 450.0_wp,  4_ik, y1) ; e1 = abs(soil_top_temp(y1, frozen) - tref)
      call march_esdirk2(y, frozen, n, nsl, 225.0_wp,  8_ik, y2) ; e2 = abs(soil_top_temp(y2, frozen) - tref)
      call march_esdirk2(y, frozen, n, nsl, 112.5_wp, 16_ik, y4) ; e4 = abs(soil_top_temp(y4, frozen) - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      call march_be_euler(y, frozen, n, nsl, 225.0_wp,  8_ik, ytmp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2,a,es10.3)', '   soil-top T errors: ', e1, ' / ', e2, ' / ', &
            e4, ' ; observed order p = ', p_lo, ' , ', p_hi, ' ; BE at 225 s ', abs(soil_top_temp(ytmp, frozen) - tref)
      call check_true('ESDIRK2 converges at least first order on the coupled soil-top temperature (p >= 0.9)', &
                      p_lo >= 0.9_wp .and. p_hi >= 0.9_wp, min(p_lo, p_hi))

      !----- (b) plant water mass stays physical (finite, positive) at production dt=900. ------------!
      call copy_state(y, y1, n) ; physical = .true.
      do step = 1_ik, 24_ik
         call esdirk2_step(y1, frozen, site, ec_default, n, nsl, 900.0_wp, ytmp, yerr, ws%stages)
         call copy_state(ytmp, y1, n)
         physical = physical .and. all(y1%leaf_water_mass(1:n) == y1%leaf_water_mass(1:n)) .and.   &
                    all(y1%wood_water_mass(1:n) == y1%wood_water_mass(1:n)) .and.                  &
                    minval(y1%leaf_water_mass(1:n)) > 0.0_wp .and. minval(y1%wood_water_mass(1:n)) > 0.0_wp
      end do
      call check_true('ESDIRK2 mass stays physical at dt=900 (split-out plant water)', physical,    &
                      minval(y1%leaf_water_mass(1:n)))

      !----- (c) embedded error estimate is bounded (not detonating). ---------------------------------!
      call esdirk2_step(y, frozen, site, ec_default, n, nsl, 900.0_wp, ynew, yerr, ws%stages)
      call state_err_norm(ynew, yerr, y, n, nsl, errnorm)
      call check_true('ESDIRK2 embedded estimate bounded (WRMS < 50) at dt=900', errnorm < 50.0_wp, errnorm)

      !----- (d) stiffness: 24 h adaptive-ESDIRK2 march stays physical + bounded. -------------------!
      ec = default_error_control(1.0e-3_wp)
      call adaptive_march(INTEG_ESDIRK2, y, frozen, site, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ws, ns, nr)
      tcas = cas_enthalpy_to_temp(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   adaptive-ESDIRK2 24 h: steps = ', ns, ' , rejects = ', nr
      call check_true('adaptive-ESDIRK2 24 h stays bounded (280 < tcas < 320 K)', tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)

      !----- (e) ERROR-CONTROL: the PI controller marches the SAME stiff window and stays               !
      !      physical -- proves the meds_fast_integ_control PI path is wired + functional. On a multi-substep !
      !      march it takes a different (typically smoother) step sequence than the I-controller.       !
      ec%controller = CTRL_PI
      call adaptive_march(INTEG_ESDIRK2, y, frozen, site, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ws, ns2, nr2)
      tcas = cas_enthalpy_to_temp(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   PI-controller  24 h: steps = ', ns2, ' , rejects = ', nr2
      call check_true('PI-controller ESDIRK2 24 h stays bounded (280 < tcas < 320 K)', tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)
   end subroutine test_esdirk2

   !----- 15. RK45 (Cash-Karp 5(4)): order-of-accuracy self-convergence.                                    !
   !          Unlike ESDIRK2's split-degraded ~1.2 order on the coupled soil-top temperature                !
   !          (test_esdirk2's part a2 -- soil water/mass are operator-split OUT of the ESDIRK stages), RK45 !
   !          carries every store in its tableau (docs/science/numerical_scheme.md §3), so the FULL         !
   !          5th-order should show up on a genuinely-coupled variable too, not just a decoupled scalar.    !
   !          Same reference-vs-halvings technique as test_esdirk2's part (a): observed order               !
   !          p = log2(e(dt)/e(dt/2)) should approach 5 as dt shrinks (some slack for pre-asymptotic         !
   !          effects: p >= 4.5). ---------------------------------------------------------------------!
   subroutine test_rk45_order()
      type(march_state_t)  :: y, yref, y1, y2, y4, y8
      type(frozen_t) :: frozen
      real(wp)    :: e1, e2, e4, e8, p_lo, p_hi, p_fine, tref, noise
      integer(ik) :: n, nsl
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rk45_order:'
      call make_column(y, frozen, n, nsl)

      !----- A 5th-order method's error ~ C*h^5 collapses to double-precision noise (~1e-13) at MUCH      !
      !      coarser h than a 2nd-order one -- h^5 falls off so fast that h must stay large enough for    !
      !      the TRUNCATION term to still dominate roundoff at the FINEST test point. Window T=2400 s     !
      !      (8 * the coarsest dt below), test points dt=300/150/75 s (8/16/32 steps; 300 s stays under   !
      !      the ~360 s explicit limit the CAS's ~130 s timescale sets -- this is a short fixed-step      !
      !      accuracy comparison, not a stability march), reference dt=2400/1024 s (>=32x finer than the  !
      !      finest test point). -------------------------------------------------------------------------!
      !----- (a) CLEAN scalar: CAS CO2 (decoupled affine ODE, mirrors test_esdirk2's part a) -- baseline  !
      !          proof the Cash-Karp tableau itself is 5th order on this RHS. ------------------------!
      call march_rk45(y, frozen, n, nsl, 2400.0_wp/1024.0_wp, 1024_ik, yref) ; tref = yref%cas_co2
      call march_rk45(y, frozen, n, nsl, 300.0_wp,   8_ik,  y1) ; e1 = abs(y1%cas_co2 - tref)
      call march_rk45(y, frozen, n, nsl, 150.0_wp,  16_ik,  y2) ; e2 = abs(y2%cas_co2 - tref)
      call march_rk45(y, frozen, n, nsl,  75.0_wp,  32_ik,  y4) ; e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)', '   CAS CO2 errors: ', e1, ' / ', e2, &
            ' / ', e4, ' ; observed order p = ', p_lo, ' , ', p_hi
      call check_true('RK45 tableau is 5th order on the decoupled CO2 twin (p >= 4.5)', &
                      p_lo >= 4.5_wp .and. p_hi >= 4.5_wp, p_hi)

      !----- (b) COUPLED: soil-top temperature (the SAME variable test_esdirk2 part a2 shows degraded to !
      !          ~1.2 order under ESDIRK2's operator split). RK45 integrates soil energy/water/mass ALL  !
      !          in its stages (no split), so this should ALSO be ~5th order. ---------------------------!
      !----- FOUR test points here, not three, and the assertion is on the two FINEST pairs. The soil-  !
      !      top temperature is the most strongly coupled variable in the column, so its coarse-h error  !
      !      carries a visible h^6 contamination and p measured across the coarsest pair is simply not   !
      !      asymptotic: depending on how strong the coupling is it reads below 5 or above it, and a     !
      !      super-convergent p > 5, which a 5th-order tableau cannot actually achieve, is contamination !
      !      too. Asserting on the finest pairs measures the tableau; asserting on the coarsest measures !
      !      whichever direction the h^6 term happens to point. -----------------------------------------!
      call march_rk45(y, frozen, n, nsl, 2400.0_wp/1024.0_wp, 1024_ik, yref) ; tref = soil_top_temp(yref, frozen)
      call march_rk45(y, frozen, n, nsl, 300.0_wp,   8_ik,  y1) ; e1 = abs(soil_top_temp(y1, frozen) - tref)
      call march_rk45(y, frozen, n, nsl, 150.0_wp,  16_ik,  y2) ; e2 = abs(soil_top_temp(y2, frozen) - tref)
      call march_rk45(y, frozen, n, nsl,  75.0_wp,  32_ik,  y4) ; e4 = abs(soil_top_temp(y4, frozen) - tref)
      call march_rk45(y, frozen, n, nsl,  37.5_wp,  64_ik,  y8) ; e8 = abs(soil_top_temp(y8, frozen) - tref)
      p_lo   = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi   = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      p_fine = log(e4/max(e8, tiny_num)) / log(2.0_wp)
      print '(a,4(es10.3,a),3(f5.2,a))', '   soil-top T errors: ', e1, ' / ', e2, ' / ', e4, ' / ',  &
            e8, ' ; observed order p = ', p_lo, ' , ', p_hi, ' , ', p_fine, ''
      !----- Soil water moves at the water network's steady rate within a step, so the soil-top        !
      !      temperature is smooth, and its error reaches the 1024-step reference's own round-off      !
      !      (~1e3 eps T) by 75 s. The order is asserted wherever the error is above that floor. ------!
      noise = 1.0e3_wp * epsilon(1.0_wp) * tref
      call check_true('RK45 stays 5th order on the COUPLED soil-top temperature (p >= 4.5 above round-off)', &
                      (p_hi >= 4.5_wp .or. e4 < noise) .and. (p_fine >= 4.5_wp .or. e8 < noise), p_fine)
   end subroutine test_rk45_order

   !----- WRMS of the embedded error estimate (mirrors adaptive_march's ESDIRK2 accept test). ------!
   subroutine state_err_norm(ynew, yerr, yref, n, nsl, errnorm)
      type(march_state_t), intent(in)  :: ynew, yerr, yref
      integer(ik),          intent(in)  :: n, nsl
      real(wp),             intent(out) :: errnorm
      real(wp), parameter :: ATE = 5.0e1_wp, ATH = 1.0e-5_wp
      real(wp)    :: s
      integer(ik) :: k, cnt
      s = (yerr%cas_enthalpy/(ATE + 1.0e-3_wp*abs(yref%cas_enthalpy)))**2 ; cnt = 1_ik
      do k = 1_ik, nsl
         s = s + (yerr%soil_energy(k)/(1.0e3_wp + 1.0e-3_wp*abs(yref%soil_energy(k))))**2
         s = s + (yerr%theta(k)/(ATH + 1.0e-3_wp*abs(yref%theta(k))))**2
         cnt = cnt + 2_ik
      end do
      errnorm = sqrt(s/real(cnt, wp))
   end subroutine state_err_norm

   !----- 16. The explicit RHS is DOMAIN-SAFE for any theta, in or out of [theta_res, theta_sat]. ------!
   ! An explicit stage's theta can sit above theta_sat: theta is never clamped (budget%theta_ood_max    !
   ! reports the excursion per sub-step). Soil psi and K evaluated there must not steer the             !
   ! trajectory, and they cannot: a stage's water tendency is the soil-plant water network's steady     !
   ! rate over the step, which does not depend on the stage's theta at all. So the water tendency at    !
   ! theta_sat + eps must equal the one at theta_sat EXACTLY -- the bit-identity asserted below, at     !
   ! eps = 0.05, far beyond any excursion the suite produces.                                           !
   !                                                                                                  !
   ! The soil-ENERGY tendency is checked separately (part c) and for a different reason: it is a pure   !
   ! flux divergence over prognostic internal energy, so theta reaches it only through the temperature  !
   ! diagnosis and the self-clamped Kersten number. It SHOULD respond to the excess water -- that       !
   ! water's heat capacity is real -- so the assertion there is finite-and-responsive, not identity.    !
   !-------------------------------------------------------------------------------------------------!
   subroutine test_rhs_domain_safety()
      type(march_state_t)  :: y, y_hi, y_lo, y_sat, y_res
      type(frozen_t) :: frozen
      type(march_tend_t)   :: f_sat, f_hi, f_res, f_lo
      real(wp), parameter   :: EPS_OOD = 0.05_wp        ! far beyond the worst excursion the suite produces
      real(wp)    :: dwater, denergy
      integer(ik) :: n, nsl, k
      logical     :: finite_all
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rhs_domain_safety:'
      call make_column(y, frozen, n, nsl)

      !----- four states: exactly at each domain edge, and far outside each. ------------------------!
      call state_init(y, n, nsl, y_sat) ; call state_init(y, n, nsl, y_hi)
      call state_init(y, n, nsl, y_res) ; call state_init(y, n, nsl, y_lo)
      do k = 1_ik, nsl
         y_sat%theta(k) = site%soil_profile%theta_sat(k)
         y_hi%theta(k)  = site%soil_profile%theta_sat(k) + EPS_OOD
         y_res%theta(k) = site%soil_profile%theta_res(k)
         y_lo%theta(k)  = max(site%soil_profile%theta_res(k) - EPS_OOD, 0.0_wp)
      end do
      call march_derivs(y_sat, frozen, site, n, nsl, f_sat)
      call march_derivs(y_hi,  frozen, site, n, nsl, f_hi)
      call march_derivs(y_res, frozen, site, n, nsl, f_res)
      call march_derivs(y_lo,  frozen, site, n, nsl, f_lo)

      !----- (a) nothing goes non-finite anywhere, on either side. ---------------------------------!
      finite_all = .true.
      do k = 1_ik, nsl
         finite_all = finite_all .and. abs(f_hi%dtheta_dt(k)) < huge(1.0_wp)                          &
                      .and. abs(f_hi%dedt(k)) < huge(1.0_wp)                                          &
                      .and. abs(f_lo%dtheta_dt(k)) < huge(1.0_wp)                                     &
                      .and. abs(f_lo%dedt(k)) < huge(1.0_wp)                                          &
                      .and. f_hi%dtheta_dt(k) == f_hi%dtheta_dt(k)                                    &
                      .and. f_lo%dedt(k) == f_lo%dedt(k)
      end do
      call check_true('RHS stays finite for theta far outside [theta_res, theta_sat]', finite_all, 0.0_wp)

      !----- (b) BIT-IDENTICAL water tendency: the network's rate does not read the stage theta. ---!
      dwater = 0.0_wp
      do k = 1_ik, nsl
         dwater = max(dwater, abs(f_hi%dtheta_dt(k) - f_sat%dtheta_dt(k)))
      end do
      call check_true('oversaturated water tendency is BIT-IDENTICAL to the saturated one', &
                      dwater == 0.0_wp, dwater)
      dwater = 0.0_wp
      do k = 1_ik, nsl
         dwater = max(dwater, abs(f_lo%dtheta_dt(k) - f_res%dtheta_dt(k)))
      end do
      call check_true('sub-residual water tendency is BIT-IDENTICAL to the residual one', &
                      dwater == 0.0_wp, dwater)

      !----- (c) the soil-ENERGY tendency is domain-safe too, but for a DIFFERENT reason worth        !
      !      recording, because it is easy to get wrong: soil_energy is prognostic INTERNAL ENERGY,   !
      !      so its tendency is a pure flux divergence -- heat capacity never divides it, and         !
      !      soil_volumetric_heat_capacity is not called on this path at all (only the implicit       !
      !      sibling uses it).                                                                        !
      !      theta reaches dedt through exactly two places, both safe:                                !
      !        * the temperature diagnosis internal_energy_to_temp(uext, theta*rho_w, ...) -- more water means    !
      !          more heat capacity means a lower T for the same internal energy, which is correct    !
      !          and is what SHOULD respond to the excess water;                                     !
      !        * soil_thermal_conductivity's Kersten number, which self-clamps at s_r = 1.            !
      !                                                                                              !
      !      make_column happens to sit on the mixed-phase MELT PLATEAU, where internal_energy_to_temp returns   !
      !      t_3ple regardless of water mass, so dedt there is bit-identical as well. That is a       !
      !      degenerate case, not the general rule, so re-seed to an all-liquid state before          !
      !      asserting the sensitivity -- otherwise this check would pass for the wrong reason. ------!
      call state_init(y, n, nsl, y_sat) ; call state_init(y, n, nsl, y_hi)
      do k = 1_ik, nsl
         y_sat%theta(k) = site%soil_profile%theta_sat(k)
         y_hi%theta(k)  = site%soil_profile%theta_sat(k) + EPS_OOD
         y_sat%soil_energy(k) = temp_to_internal_energy(site%soil_profile%dry_heat_capacity(k),                        &
                                             y_sat%theta(k)*rho_h2o, 290.0_wp, 1.0_wp)
         y_hi%soil_energy(k)  = y_sat%soil_energy(k)      ! SAME internal energy, more water
      end do
      call march_derivs(y_sat, frozen, site, n, nsl, f_sat)
      call march_derivs(y_hi,  frozen, site, n, nsl, f_hi)
      denergy = 0.0_wp
      do k = 1_ik, nsl
         denergy = max(denergy, abs(f_hi%dedt(k) - f_sat%dedt(k)))
      end do
      call check_true('oversaturated soil-energy tendency stays finite and responds through T only',   &
                      denergy > 0.0_wp .and. denergy < 1.0e6_wp, denergy)
      !----- and the water tendency is STILL bit-identical on this warm state, i.e. (b) was not an     !
      !      artefact of the plateau either. ---------------------------------------------------------!
      dwater = 0.0_wp
      do k = 1_ik, nsl
         dwater = max(dwater, abs(f_hi%dtheta_dt(k) - f_sat%dtheta_dt(k)))
      end do
      call check_true('water tendency bit-identical on an all-liquid state too (not a plateau artefact)', &
                      dwater == 0.0_wp, dwater)
   end subroutine test_rhs_domain_safety

   !----- shared full-column setup (2 cohorts, 10 soil layers, a representative daytime state). ----!
   subroutine make_column(y, frozen, n, nsl)
      type(march_state_t),  intent(out) :: y
      type(frozen_t), intent(out) :: frozen
      integer(ik),           intent(in)  :: n, nsl
      !----- LOCAL hydraulics traits, needed only to seed y%*_water_mass at a representative psi:     !
      !      march_derivs' mass ODE takes only the frozen sapflow/uptake, not the PV-curve or        !
      !      conductance parameters, so frozen_t does not carry them. -------------------------!
      type(hydraulics_pft_t)   :: hp
      integer(ik) :: i, k
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, test_bottoms)
      call build_profile_hydraulics(10_ik, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,      &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, site%soil_profile)
      call build_profile_thermal(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, site%soil_profile)
      site%soil = soil_params_t()
      hp%leaf_curve = water_curve_t(pi0 = -1.5_wp, elastic_mod = 12.0_wp, apoplast_frac = 0.30_wp, water_sat = 2.0_wp)
      hp%wood_curve = water_curve_t(pi0 = -1.0_wp, elastic_mod =  8.0_wp, apoplast_frac = 0.20_wp, water_sat = 1.0_wp)
      hp%wood_psi50 = -2.0_wp
      hp%wood_kexp = 2.0_wp ; hp%wood_kmax = 8.0_wp ; hp%vessel_curl = 1.5_wp
      frozen%hydrology%q_top = 1.0e-6_wp
      allocate(frozen%roots%root_share(nsl), frozen%plant%nplant(n))
      allocate(frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n), frozen%roots%qloss_frozen(n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n), frozen%film%transfer_lw(n))
      frozen%film%transfer_lw = 0.0_wp
      frozen%roots%root_share(1:nsl) = exp_root_profile(site%soil_profile, nsl)   ! per-layer sink placement
      frozen%plant%nplant = 0.3_wp
      !----- FROZEN sapflow/uptake (docs/science/numerical_scheme.md §2): a representative,             !
      !      state-independent pair for these RHS and reference-march tests -- not taken from the water  !
      !      network here (these tests exercise march_derivs/the tableau machinery directly, not       !
      !      build_frozen, which has its own coverage via test_column_esdirk2.f90). -------------!
      frozen%plant%sapflow_frozen = 1.0e-4_wp
      frozen%plant%uptake_frozen = 1.0e-4_wp
      !----- The stages hold each cohort's transpiration at the water network's, as in the model      !
      !      (build_frozen): here the sapflow, so the leaves neither fill nor empty over the     !
      !      hours these marches run on one frozen record. -----------------------------------------!
      allocate(frozen%tissue%transp_net(n))
      frozen%tissue%transp_net(1:n) = frozen%plant%sapflow_frozen(1:n) * frozen%plant%nplant(1:n)
      frozen%tissue%transp_held = .true.
      frozen%roots%uptake = frozen%plant%uptake_frozen(1)*sum(frozen%plant%nplant(1:n))
      frozen%roots%qloss_frozen = 0.0_wp   ! advective enthalpy: no-op unless populated (see build_frozen)
      !----- The water network's step, as build_frozen hands it on. Held still here -- no face  !
      !      flow, drainage, uptake heat or theta change -- because these marches run for hours on one   !
      !      frozen record, which the model never does: the network sets the rates afresh every step.   !
      !      test_column_assembler sets them nonzero for its wiring check. ---------------------------!
      allocate(frozen%roots%root_heat(nsl))
      frozen%roots%root_heat(1:nsl) = 0.0_wp
      frozen%hydrology%t_infil = 296.0_wp ; frozen%hydrology%t_bot = 292.0_wp
      frozen%film%intercept_leaf = 0.0_wp ; frozen%film%intercept_wood = 0.0_wp   ! canopy water: no-op unless populated
      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n), frozen%tissue%abs_lw(n), &
               frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      frozen%tissue%leaf_hcap_per_dt = 0.0_wp ; frozen%tissue%wood_hcap_per_dt = 0.0_wp
      frozen%tissue%t_leaf0 = 298.5_wp ; frozen%tissue%t_wood0 = 298.5_wp   ! the emission base, t_rad
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      allocate(frozen%tissue%leaf_emiss(n), frozen%tissue%wood_emiss(n))
      frozen%tissue%leaf_emiss = 0.95_wp ; frozen%tissue%wood_emiss = 0.95_wp
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n))
      frozen%tissue%h_coeff_w = 0.0_wp
      frozen%tissue%abs_sw_wood = 0.0_wp
      frozen%tissue%abs_lw_wood = 0.0_wp
      frozen%tissue%wai = 0.0_wp
      frozen%tissue%qwflux_wl = 0.0_wp ; frozen%tissue%q_wood_net = 0.0_wp
      frozen%film%f_wet_c = 0.0_wp ; frozen%film%g_film_leaf = 0.0_wp ; frozen%film%g_film_w = 0.0_wp
      do i = 1_ik, n
         frozen%tissue%lai(i) = 2.0_wp - 0.5_wp * real(i-1_ik, wp) ; frozen%tissue%abs_sw(i) = 250.0_wp - 50.0_wp*real(i-1_ik, wp)
         frozen%tissue%abs_lw(i) = -30.0_wp
         frozen%tissue%h_coeff_leaf(i) = 2.0_wp * frozen%tissue%lai(i) * 0.03_wp * 1.2_wp * cp_air
         frozen%tissue%g_transp_leaf(i) = 0.004_wp * frozen%tissue%lai(i)
      end do
      frozen%cas%rho = 1.2_wp ; frozen%cas%press = 101325.0_wp ; frozen%cas%nee_biotic = -5.0_wp
      associate (b => frozen%cas%box)
         b%air_mass_capacity = 1.2_wp*20.0_wp
         b%air_molar_capacity = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*20.0_wp
         b%atm_conductance_enthalpy = 1.2_wp*0.3_wp*0.02_wp ; b%atm_conductance_vapor = b%atm_conductance_enthalpy
         b%atm_conductance_co2 = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*0.3_wp*0.02_wp
         b%atm_enthalpy = cas_temp_to_enthalpy(300.0_wp, 0.011_wp) ; b%atm_specific_humidity = 0.011_wp
         b%atm_co2 = 400.0_wp
      end associate
      frozen%ground%rnet_soil = 50.0_wp
      frozen%ground%ggnet = 0.02_wp ; frozen%ground%soil_evap = 2.0e-5_wp
      y%cas_enthalpy = cas_temp_to_enthalpy(297.0_wp, 0.012_wp) ; y%cas_shv = 0.012_wp ; y%cas_co2 = 410.0_wp
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water = 0.0_wp ; y%wood_surf_water = 0.0_wp   ! canopy water: no-op unless populated
      do k = 1_ik, nsl
         y%soil_energy(k) = temp_to_internal_energy(site%soil_profile%dry_heat_capacity(k), 0.30_wp*rho_h2o,     &
                            296.0_wp - 0.4_wp*real(k-1_ik, wp), 1.0_wp)
         y%theta(k) = 0.30_wp
      end do
      !----- seed mass at a representative point (leaf psi=-1.0, wood psi=-0.5 MPa), via the forward !
      !      tissue_psi_to_water map (hp above). ----------------------------------------------------!
      do i = 1_ik, n
         y%leaf_water_mass(i) = tissue_psi_to_water(-1.0_wp, hp%leaf_curve, 0.5_wp)          ! 0.5 kgC of leaf
         y%wood_water_mass(i) = tissue_psi_to_water(-0.5_wp, hp%wood_curve, 7.0_wp)          ! 5 sapwood + 2 fine root
      end do
   end subroutine make_column

   subroutine copy_state(a, b, n)
      type(march_state_t), intent(in)  :: a
      type(march_state_t), intent(out) :: b
      integer(ik),          intent(in)  :: n
      b%cas_enthalpy = a%cas_enthalpy ; b%cas_shv = a%cas_shv ; b%cas_co2 = a%cas_co2
      b%soil_energy = a%soil_energy ; b%theta = a%theta
      if (allocated(b%leaf_water_mass)) deallocate(b%leaf_water_mass, b%wood_water_mass)
      allocate(b%leaf_water_mass(n), b%wood_water_mass(n))
      b%leaf_water_mass(1:n) = a%leaf_water_mass(1:n) ; b%wood_water_mass(1:n) = a%wood_water_mass(1:n)
      if (allocated(b%leaf_surf_water)) deallocate(b%leaf_surf_water, b%wood_surf_water)
      allocate(b%leaf_surf_water(n), b%wood_surf_water(n))
      b%leaf_surf_water(1:n) = a%leaf_surf_water(1:n) ; b%wood_surf_water(1:n) = a%wood_surf_water(1:n)
   end subroutine copy_state

   !----- helper: the soil-top temperature diagnosed from the state (mirrors      !
   !      what march_derivs does internally) so the standalone CAS comparison lines up. -----------!
   function tground_of(y, frozen) result(tg)
      type(march_state_t),   intent(in) :: y
      type(frozen_t),  intent(in) :: frozen
      real(wp) :: tg, fl
      call soil_layer_temp(y%soil_energy(1), y%theta(1), site%soil_profile%dry_heat_capacity(1),                 &
                           site%soil_profile%layer_curve(1), tg, fl)
   end function tground_of

   logical function ieee_ok(x)
      real(wp), intent(in) :: x
      ieee_ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function ieee_ok

end program test_column_derivs
