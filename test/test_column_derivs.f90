!==========================================================================================!
! test_column_derivs -- unit tests for the pure fast-loop RHS (P0 increment 1: surface subsystem). !
! The IMEX-ARK overhaul (docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md) needs a side-effect-free f(y)=dy/dt.     !
! These tests prove `surface_derivs` is a faithful extraction of the split fast step's surface      !
! path, so an ARK stepper can trust it:                                                            !
!   1. LEAF DIAGNOSTIC CLOSURE   -- the diagnosed leaf temperature zeroes the linearized leaf        !
!                                    energy balance it was solved from (multi-cohort).               !
!   2. ANALYTIC LEAF TEMPERATURE -- a hand-computable no-latent / no-LW-forcing case matches.        !
!   3. CAS BACKWARD-EULER CONSISTENCY -- the split's committed enth1/shv1/co21 (the BE-in-atm        !
!                                    solution of the tendencies) satisfy the implicit relation        !
!                                    (y1-y0)/dt = f(y1); i.e. d_cas_* is the correct RHS.             !
!   4. CONSERVATION MARCH        -- marching the CAS twins with the RHS closes the energy/water        !
!                                    budgets to round-off and warms the CAS toward a warm atmosphere. !
!==========================================================================================!
program test_column_derivs
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : latent_heat_vap, stefan, cp_air, tiny_num, rho_h2o
   use meds_therm_lib,         only : cas_enthalpy_of_temp, cas_temp_of_enthalpy,                   &
                                   sat_specific_humidity, sat_vapor_pressure,                    &
                                   sat_vapor_pressure_temp_deriv, enthalpy_vapor, temp_to_internal_energy, internal_energy_to_temp
   use meds_budget_check,   only : budget_t, budget_accumulate
   use meds_biophysics_types, only : energy_forcing_t, energy_flux_t
   use meds_column_constants, only : n_soil_layer_max
   use meds_column_reservoirs, only : soil_energy_column_t
   use meds_column_params, only : soil_params_t, soil_thermal_params_t
   use meds_hydr_lib, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : energy_opts_t, soil_opts_t
   use meds_column_params, only : build_soil_hydr_params
   use meds_column_params, only : build_soil_therm_params
   use meds_soil_energy,      only : soil_energy_step_implicit, soil_energy_time_deriv
   use meds_soil_water,       only : soil_water_time_deriv
   use meds_plant_types, only : hydro_params_t, hydro_opts_t
   use meds_hydr_lib,         only : water_content
   use meds_fast_time_derivs, only : surface_derivs, column_derivs
   use meds_therm_lib,        only : internal_energy_liquid
   use meds_fast_types,       only : surface_state_t, surface_tend_t,           &
                                   column_state_t, column_frozen_t, column_tend_t
   use meds_fast_rk4_oracle,  only : rk4_column_step, imex_euler_column_step, adaptive_imex_march
   use meds_fast_ark,         only : ark2_column_step, adaptive_ark_march
   use meds_column_state_ops, only : state_init
   use meds_fast_rk45,        only : rk45_column_step
   use meds_fast_control,     only : default_error_control
   use meds_fast_types,       only : error_control_t
   use meds_config,           only : CTRL_PI
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_leaf_closure()
   call test_leaf_analytic()
   call test_cas_be_consistency()
   call test_conservation_march()
   call test_soil_energy_tendency()
   call test_soil_water_tendency()
   call test_column_assembler()
   call test_rk4_march()
   call test_imex_euler()
   call test_imex_coupled()
   call test_arrowhead()
   call test_adaptive_march()
   call test_ark2()
   call test_rk45_order()

   !----- 15. CONSTITUTIVE-DOMAIN SAFETY of the explicit RHS (issue #78 item 2). ------------------!
   call test_rhs_domain_safety()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_derivs: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_derivs: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (', val, ')'
      else
         nfail = nfail + 1_ik ; print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

   !----- A representative 3-cohort daytime surface setup (frozen pre-pass + aerodynamics). -----!
   subroutine make_frozen(frozen, n)
      type(column_frozen_t),  intent(out) :: frozen
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
      frozen%tissue%t_leaf0 = 0.0_wp ; frozen%tissue%t_wood0 = 0.0_wp
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n))
      frozen%tissue%h_coeff_w = 0.0_wp ; frozen%tissue%abs_sw_wood = 0.0_wp
      frozen%tissue%abs_lw_wood = 0.0_wp ; frozen%tissue%wai = 0.0_wp
      frozen%tissue%qwflux_wl = 0.0_wp ; frozen%tissue%q_wood_net = 0.0_wp   ! P2 advective enthalpy: no-op unless populated
      frozen%film%f_wet_c = 0.0_wp ; frozen%film%g_film_leaf = 0.0_wp
      frozen%film%g_film_w = 0.0_wp   ! canopy water: no-op unless set
      do i = 1_ik, n
         frozen%tissue%lai(i)       = 2.0_wp - 0.4_wp * real(i - 1_ik, wp)          ! 2.0, 1.6, 1.2
         frozen%tissue%abs_sw(i)    = 250.0_wp - 40.0_wp * real(i - 1_ik, wp)       ! more light at the top
         frozen%tissue%abs_lw(i)    = -30.0_wp
         frozen%tissue%h_coeff_leaf(i) = 2.0_wp * frozen%tissue%lai(i) * 0.03_wp * 1.2_wp * cp_air  ! effarea*lai*gbh*rho*cp
         frozen%tissue%g_transp_leaf(i)    = 0.004_wp * frozen%tissue%lai(i)                          ! series conductance [m/s]
      end do
      frozen%tissue%leaf_emiss    = 0.95_wp
      frozen%cas%rho           = 1.2_wp
      frozen%cas%press         = 101325.0_wp
      frozen%cas%cas_mass_capacity          = frozen%cas%rho * 20.0_wp                                ! rho * can_depth
      frozen%cas%cas_molar_capacity          = (frozen%cas%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 20.0_wp
      frozen%cas%g_atm_heat           = frozen%cas%rho * 0.3_wp * 0.02_wp                       ! rho * ustar * temp1
      frozen%cas%g_atm_vapour           = frozen%cas%rho * 0.3_wp * 0.02_wp
      frozen%cas%g_atm_co2           = (frozen%cas%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 0.3_wp * 0.02_wp
      frozen%cas%enthalpy_atm      = cas_enthalpy_of_temp(300.0_wp, 0.011_wp)
      frozen%cas%shv_atm       = 0.011_wp
      frozen%cas%co2_atm       = 400.0_wp
      frozen%cas%nee_biotic    = -5.0_wp                                          ! net CO2 uptake [umol/m2/s]
      frozen%ground%abs_sw_ground = 60.0_wp
      frozen%ground%abs_lw_ground = -10.0_wp
      frozen%ground%ggnet         = 0.02_wp
      frozen%ground%soil_evap     = 2.0e-5_wp
   end subroutine make_frozen

   !----- 1. The diagnosed leaf temperature must zero the (linearized) leaf energy balance. ------!
   subroutine test_leaf_closure()
      type(column_frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat, dtl, lw_slope, le_slope, le_ref, resid, worst
      integer(ik) :: i, n
      n = 3_ik
      print '(a)', 'test_leaf_closure:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(298.0_wp, 0.012_wp)
      y%cas_shv      = 0.012_wp
      y%cas_co2      = 410.0_wp
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, 297.0_wp, n, f)

      tcas   = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, frozen%cas%press)
      esat   = sat_vapor_pressure(tcas)
      dqdt   = 0.622_wp * frozen%cas%press / max((frozen%cas%press - 0.378_wp * esat) ** 2, tiny_num)          &
               * sat_vapor_pressure_temp_deriv(tcas)
      worst = 0.0_wp
      do i = 1_ik, n
         dtl      = f%leaf_temp(i) - tcas
         lw_slope = 4.0_wp * frozen%tissue%leaf_emiss * stefan * tcas ** 3 * frozen%tissue%lai(i)
         !----- the leaf pays the FULL vapour enthalpy at the linearization temperature (2026-09). -!
         le_slope = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * dqdt
         le_ref   = enthalpy_vapor(tcas) * frozen%cas%rho * frozen%tissue%g_transp_leaf(i) * (qsat_c - qcas)
         !----- Rnet - sensible - latent - LW-emission, all at the diagnosed leaf temperature. ---!
         resid = frozen%tissue%abs_sw(i) + frozen%tissue%abs_lw(i) - frozen%tissue%h_coeff_leaf(i) * dtl                          &
                 - (le_ref + le_slope * dtl) - lw_slope * dtl
         worst = max(worst, abs(resid))
      end do
      call check_true('diagnosed leaf T zeroes the linearized balance', worst < 1.0e-9_wp, worst)
      call check_true('leaf temperatures physical (270-330 K)',                                 &
                      minval(f%leaf_temp) > 270.0_wp .and. maxval(f%leaf_temp) < 330.0_wp,       &
                      f%leaf_temp(1))
   end subroutine test_leaf_closure

   !----- 2. No latent (g_transp_leaf = 0) and no LW forcing (abs_lw = 0): dtl = abs_sw / (h + lw_slope). !
   subroutine test_leaf_analytic()
      type(column_frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp) :: tcas, lw_slope, expect
      print '(a)', 'test_leaf_analytic:'
      call make_frozen(frozen, 1_ik)
      frozen%tissue%g_transp_leaf(1) = 0.0_wp ; frozen%tissue%abs_lw(1) = 0.0_wp ; frozen%tissue%abs_sw(1) = 200.0_wp
      frozen%tissue%lai(1) = 1.5_wp
      frozen%tissue%h_coeff_leaf(1) = 2.0_wp * frozen%tissue%lai(1) * 0.03_wp * 1.2_wp * cp_air
      y%cas_enthalpy = cas_enthalpy_of_temp(299.0_wp, 0.010_wp)
      y%cas_shv      = 0.010_wp
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, 297.0_wp, 1_ik, f)
      tcas     = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      lw_slope = 4.0_wp * frozen%tissue%leaf_emiss * stefan * tcas ** 3 * frozen%tissue%lai(1)
      expect   = tcas + frozen%tissue%abs_sw(1) / (frozen%tissue%h_coeff_leaf(1) + lw_slope)
      call check('leaf T = tcas + Rn/(h+lw_slope) (no latent, no LW)', f%leaf_temp(1), expect, 1.0e-10_wp)
   end subroutine test_leaf_analytic

   !----- 3. The split commits the BE-in-atm solution of the tendencies; verify (y1-y0)/dt = f(y1). !
   subroutine test_cas_be_consistency()
      type(column_frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: dt, enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: r_enth, r_shv, r_co2
      integer(ik) :: n
      n = 3_ik ; dt = 900.0_wp
      print '(a)', 'test_cas_be_consistency:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(296.0_wp, 0.013_wp)
      y%cas_shv      = 0.013_wp
      y%cas_co2      = 415.0_wp
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, 297.0_wp, n, f)
      !----- Reconstruct the split's committed state (meds_fast_split.f90). ------------------------!
      associate (c => frozen%cas)
         enth1 = (c%cas_mass_capacity * enth0 + dt * (f%src_enth  + c%g_atm_heat * c%enthalpy_atm))          &
                 / (c%cas_mass_capacity + dt * c%g_atm_heat)
         shv1  = (c%cas_mass_capacity * shv0  + dt * (f%src_vap   + c%g_atm_vapour * c%shv_atm ))        &
                 / (c%cas_mass_capacity + dt * c%g_atm_vapour)
         co21  = (c%cas_molar_capacity * co20  + dt * (c%nee_biotic + c%g_atm_co2 * c%co2_atm))          &
                 / (c%cas_molar_capacity + dt * c%g_atm_co2)
      end associate
      !----- BE consistency: (y1-y0)/dt must equal the tendency evaluated with the ATM term at y1  !
      !      (source frozen). This is exactly what an IMEX/BE stage solves, so it verifies d_cas_* !
      !      is the correct RHS of the split's implicit update.                                     !
      r_enth = (enth1 - enth0) / dt - (f%src_enth + frozen%cas%g_atm_heat * (frozen%cas%enthalpy_atm - enth1))  &
               / frozen%cas%cas_mass_capacity
      r_shv  = (shv1  - shv0 ) / dt - (f%src_vap  + frozen%cas%g_atm_vapour * (frozen%cas%shv_atm  - shv1 )) &
               / frozen%cas%cas_mass_capacity
      r_co2  = (co21  - co20 ) / dt - (frozen%cas%nee_biotic + frozen%cas%g_atm_co2 * (frozen%cas%co2_atm - co21)) &
               / frozen%cas%cas_molar_capacity
      call check_true('CAS enthalpy BE-consistent with d_cas_enthalpy', abs(r_enth) < 1.0e-9_wp, r_enth)
      call check_true('CAS humidity BE-consistent with d_cas_shv',      abs(r_shv)  < 1.0e-15_wp, r_shv)
      call check_true('CAS CO2 BE-consistent with d_cas_co2',           abs(r_co2)  < 1.0e-9_wp, r_co2)
      !----- At t=0 the tendency must equal (src + g*(atm - y0))/cap (sanity on the returned RHS). !
      call check('d_cas_enthalpy = (src+g_atm_heat*(atm-y0))/cas_mass_capacity',                                      &
                 f%d_cas_enthalpy, (f%src_enth + frozen%cas%g_atm_heat * (frozen%cas%enthalpy_atm - enth0))     &
                                   / frozen%cas%cas_mass_capacity, 1.0e-9_wp)
   end subroutine test_cas_be_consistency

   !----- 4. March the CAS twins with the RHS (BE-in-atm); the closed budgets must stay tight. ---!
   subroutine test_conservation_march()
      type(column_frozen_t)  :: frozen
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      type(budget_t)         :: be, bw
      real(wp)    :: dt, enth0, shv0, enth1, shv1, t0, t1
      integer(ik) :: step, n
      n = 3_ik ; dt = 120.0_wp
      print '(a)', 'test_conservation_march:'
      call make_frozen(frozen, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(294.0_wp, 0.010_wp)     ! CAS starts cool + dry
      y%cas_shv      = 0.010_wp
      y%cas_co2      = 400.0_wp
      t0 = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      do step = 1_ik, 60_ik
         call surface_derivs(y, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, 297.0_wp, n, f)
         enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
         associate (c => frozen%cas)
            enth1 = (c%cas_mass_capacity * enth0 + dt * (f%src_enth + c%g_atm_heat * c%enthalpy_atm))        &
                    / (c%cas_mass_capacity + dt * c%g_atm_heat)
            shv1  = (c%cas_mass_capacity * shv0  + dt * (f%src_vap  + c%g_atm_vapour * c%shv_atm ))      &
                    / (c%cas_mass_capacity + dt * c%g_atm_vapour)
         end associate
         !----- Same closed-budget accounting the split uses (meds_fast_split.f90). ------------------!
         call budget_accumulate(be, frozen%cas%cas_mass_capacity * enth0, frozen%cas%cas_mass_capacity * enth1, &
                                f%src_enth + frozen%cas%g_atm_heat * frozen%cas%enthalpy_atm,                          &
                                frozen%cas%g_atm_heat * enth1, dt, abs(frozen%cas%cas_mass_capacity * enth1), 1.0e-8_wp, 1.0e-3_wp)
         call budget_accumulate(bw, frozen%cas%cas_mass_capacity * shv0, frozen%cas%cas_mass_capacity * shv1, &
                                f%src_vap + frozen%cas%g_atm_vapour * frozen%cas%shv_atm,                            &
                                frozen%cas%g_atm_vapour * shv1, dt, max(abs(frozen%cas%cas_mass_capacity * shv1), 1.0e-6_wp), &
                                     1.0e-8_wp, 1.0e-10_wp)
         y%cas_enthalpy = enth1 ; y%cas_shv = shv1
      end do
      t1 = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      call check_true('CAS energy budget closes over the march', be%n_fail == 0_ik, real(be%worst, wp))
      call check_true('CAS water  budget closes over the march', bw%n_fail == 0_ik, real(bw%worst, wp))
      call check_true('CAS warms toward the warm atmosphere + sunlit canopy', t1 > t0, t1 - t0)
   end subroutine test_conservation_march


   !----- 6. soil_energy_time_deriv == the dt->0 limit of the BE kernel soil_energy_step_implicit. ----------!
   subroutine test_soil_energy_tendency()
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: therm
      type(energy_forcing_t)      :: forcing
      type(energy_opts_t)         :: eopts
      type(soil_energy_column_t)  :: col, col0
      type(energy_flux_t)         :: eflux
      real(wp)    :: dedt(n_soil_layer_max), dt, worst, fd
      integer(ik) :: k, nsl
      nsl = 10_ik
      print '(a)', 'test_soil_energy_tendency:'
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,        &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, soil)
      call build_soil_therm_params(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, therm)
      forcing%soil_water(1:10) = 0.30_wp ; forcing%w_flux = 0.0_wp
      forcing%g_top = 120.0_wp ; forcing%geothermal = 0.0_wp
      forcing%root_heat_sink(1:10) = 3.0_wp                       ! nonzero interior sink
      do k = 1_ik, 10_ik                                          ! a temperature gradient -> nonzero faces
         col%soil_energy(k) = temp_to_internal_energy(therm%soil_dry_heat_capacity(k),                      &
                              forcing%soil_water(k) * rho_h2o, 291.0_wp - 0.6_wp * real(k-1_ik, wp), 1.0_wp)
      end do
      col0 = col
      call soil_energy_time_deriv(col0, forcing, therm, soil, eopts, dedt)
      dt = 1.0e-2_wp
      call soil_energy_step_implicit(col, forcing, therm, soil, eopts, dt, eflux)   ! advances col (BE)
      worst = 0.0_wp
      do k = 1_ik, nsl
         fd    = (col%soil_energy(k) - col0%soil_energy(k)) / dt
         worst = max(worst, abs(fd - dedt(k)) / max(abs(dedt(k)), 1.0_wp))
      end do
      call check_true('soil-energy tendency = BE-limit of soil_energy_step_implicit (O(dt))', worst < 1.0e-2_wp, worst)
   end subroutine test_soil_energy_tendency

   !----- 7. soil_water_time_deriv closes the column water balance (telescoping faces). ------------!
   subroutine test_soil_water_tendency()
      type(soil_params_t) :: soil
      type(soil_opts_t)   :: hopts
      real(wp)    :: theta(n_soil_layer_max), psi_e(n_soil_layer_max), root_uptake(n_soil_layer_max)
      real(wp)    :: dtheta(n_soil_layer_max), drain, uptk, q_top, net, colsum
      real(wp)    :: qface(n_soil_layer_max), face_sum
      integer(ik) :: k, nsl
      nsl = 10_ik
      print '(a)', 'test_soil_water_tendency:'
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,        &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, soil)
      hopts = soil_opts_t()
      do k = 1_ik, 10_ik ; theta(k) = 0.26_wp + 0.015_wp * real(k-1_ik, wp) ; end do   ! moist gradient
      psi_e = 0.0_wp ; root_uptake = 0.0_wp
      root_uptake(1:5) = 1.0e-5_wp                     ! [kg/m2/s] root sink in the top half
      q_top = 1.0e-6_wp                                ! [m/s] gentle infiltration
      call soil_water_time_deriv(theta, soil, hopts, nsl, q_top, root_uptake, dtheta, drain,   &
                                 uptk, qface)
      colsum = 0.0_wp
      do k = 1_ik, nsl ; colsum = colsum + dtheta(k) * soil%dz(k) ; end do
      net = q_top - drain / rho_h2o - uptk / rho_h2o   ! d(storage)/dt = in - drainage - uptake
      call check('soil-water column balance closes (telescoping)', colsum, net, 1.0e-12_wp)
      call check_true('root sink active and column free-drains', uptk > 0.0_wp .and. drain > 0.0_wp, drain)
      !----- The EXPORTED interior faces (#78 item 3) must be the SAME faces this dtheta_dt was built   !
      !      from, since an explicit integrator advects soil enthalpy on them: a face flux measured on     !
      !      one trajectory and a mass change on another is exactly the defect that issue is about.        !
      !      Verify per layer rather than in aggregate -- a column SUM telescopes and would stay green     !
      !      with the whole profile shifted by a constant, which is the vertical-only error mode that      !
      !      whole-column ledgers cannot see. -----------------------------------------------------------!
      !      Layer k's storage change is (face in - face out) minus its own root sink, and the sink is    !
      !      not exported per layer -- so check the layers where it is zero (root_uptake is set on 1..5    !
      !      in this fixture), which still covers both an interior pair and the bottom face. -------------!
      face_sum = 0.0_wp
      do k = 6_ik, nsl
         net = qface(k-1_ik)
         if (k <= nsl-1_ik) then
            net = net - qface(k)
         else
            net = net - drain / rho_h2o                     ! bottom face = the reported drainage
         end if
         face_sum = max(face_sum, abs(dtheta(k) * soil%dz(k) - net))
      end do
      call check_true('exported interior faces reproduce dtheta_dt layer by layer',                     &
                      face_sum < 1.0e-16_wp, face_sum)

      !----- REVIEW 2026-09 (item 2 #3): a caller whose root_uptake is ALREADY the realized,        !
      !      psi-limited sink (advance_soil_water_column's uptake_total, which the plant water ODE debits   !
      !      from wood) must be able to hand it over as-is. On soil inside the wilting ramp the         !
      !      default path limits it (uptk < sum), the apply_wilt_limit=.false. path must not. ---------!
      block
         real(wp) :: theta_dry(n_soil_layer_max), uptk_limited, uptk_asis, requested
         theta_dry = 0.10_wp                                  ! psi ~ -39 m: psi_open (-3.37) > psi > psi_wilt (-153)
         requested = sum(root_uptake(1:nsl))
         call soil_water_time_deriv(theta_dry, soil, hopts, nsl, q_top, root_uptake, dtheta, drain, &
                                    uptk_limited, qface)
         call soil_water_time_deriv(theta_dry, soil, hopts, nsl, q_top, root_uptake, dtheta, drain, &
                                    uptk_asis, qface, apply_wilt_limit=.false.)
         call check_true('dry soil: default path applies the wilting ramp (uptake < requested)',     &
                         uptk_limited > 0.0_wp .and. uptk_limited < 0.99_wp * requested, uptk_limited / requested)
         call check('dry soil: apply_wilt_limit=.false. takes the sink as-is (uptake == requested)',   &
                    uptk_asis, requested, 1.0e-12_wp * requested)
         colsum = 0.0_wp
         do k = 1_ik, nsl ; colsum = colsum + dtheta(k) * soil%dz(k) ; end do
         net = q_top - drain / rho_h2o - uptk_asis / rho_h2o
         call check('dry soil: as-is sink still telescopes into the column balance', colsum, net, 1.0e-12_wp)
      end block
   end subroutine test_soil_water_tendency

   !----- 9. column_derivs assembles a finite, physically-signed whole-column RHS + its CAS part   !
   !         matches a standalone surface_derivs call (correct wiring). -------------------------!
   subroutine test_column_assembler()
      type(column_state_t)       :: y
      type(column_frozen_t)      :: frozen
      type(column_tend_t)        :: f
      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(energy_forcing_t)     :: eforc_chk
      type(soil_energy_column_t) :: se_chk
      real(wp)    :: dedt_chk(n_soil_layer_max), worst
      real(wp)    :: uptake_chk(n_soil_layer_max), dtheta_chk(n_soil_layer_max)
      real(wp)    :: qface_chk(n_soil_layer_max)
      real(wp)    :: drain_chk, uptk_chk
      integer(ik) :: k, i, n, nsl
      logical     :: all_finite
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_column_assembler:'
      call make_column(y, frozen, n, nsl)
      call column_derivs(y, frozen, n, nsl, f)

      all_finite = ieee_ok(f%d_cas_enthalpy) .and. ieee_ok(f%d_cas_shv) .and. ieee_ok(f%d_cas_co2)
      do k = 1_ik, nsl ; all_finite = all_finite .and. ieee_ok(f%dedt(k)) .and. ieee_ok(f%dtheta_dt(k)) ; end do
      do i = 1_ik, n
         all_finite = all_finite .and. ieee_ok(f%d_leaf_water_mass(i)) .and. ieee_ok(f%d_wood_water_mass(i))
      end do
      call check_true('column_derivs: all tendencies finite', all_finite, 0.0_wp)

      !----- the CAS part must equal a standalone surface_derivs at the state's diagnosed t_ground. !
      y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,     &
                          tground_of(y, frozen), n, surf_tend)
      call check('column_derivs CAS enthalpy tendency = surface_derivs', f%d_cas_enthalpy, surf_tend%d_cas_enthalpy, 1.0e-12_wp)

      !----- REVIEW 2026-09 (item 1A #10): the CANOPY is energy-neutral. With no advected enthalpy    !
      !      (qwflux_wl = q_wood_net = 0 in this fixture) the tissues hold no store on this path, so    !
      !      what the canopy absorbs (coh_rnet) must equal what it hands the CAS: sensible + the FULL   !
      !      enthalpy of the vapour it sheds. The old code handed the CAS transp*enthalpy_vapor while   !
      !      the leaf paid only latent_heat_vap, and charged the difference to the SOIL (coh_qsoil), so  !
      !      this identity was off by exactly that proxy (~30 W/m2 at 3 mm/day) and the soil paid the   !
      !      transpired water's liquid enthalpy twice once P2 added qloss. --------------------------!
      block
         real(wp) :: canopy_to_cas, tcas_chk
         tcas_chk = cas_temp_of_enthalpy(y_stage%cas_enthalpy, y_stage%cas_shv)
         canopy_to_cas = surf_tend%src_enth - surf_tend%h_ground - surf_tend%le_ground + surf_tend%cond * &
               internal_energy_liquid(tcas_chk)
         call check_true('canopy transpires at all in this fixture (test is live)', surf_tend%coh_transp > 1.0e-7_wp, &
               surf_tend%coh_transp)
         call check('canopy is energy-neutral: coh_rnet == sensible + full vapour enthalpy to the CAS (no soil proxy)', &
                    canopy_to_cas, surf_tend%coh_rnet, 1.0e-9_wp * max(abs(surf_tend%coh_rnet), 1.0_wp))
      end block

      !----- wiring check: the assembled soil-heat tendency == a standalone soil_energy_time_deriv    !
      !      built from the SAME surface coupling (g_top, qloss * root_share) PLUS the bottom-face   !
      !      drainage enthalpy AND the interior advective faces. Both of those extra terms were once    !
      !      omitted here and the check still passed, each time for the same reason -- column_derivs    !
      !      advected them on a FROZEN quantity that happens to be 0 in this fixture (frozen%hydrology%drainage,     !
      !      then frozen%hydrology%w_flux_frozen). C2 made the bottom face ride the state-dependent f%drainage_rate  !
      !      and #78 item 3 made the interior faces ride this stage's own theta trajectory, so both are !
      !      now nonzero and the reproduction has to carry them. A check that omits a term is only      !
      !      green while that term is zero, and this one has now taught that lesson twice. ------------!
      se_chk%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc_chk%g_top = surf_tend%g_top ; eforc_chk%geothermal = frozen%hydrology%geothermal
      do k = 1_ik, nsl
         uptake_chk(k) = frozen%roots%uptake * frozen%params%soil%root_frac(k)
      end do
      call soil_water_time_deriv(y%theta, frozen%params%soil, frozen%params%hydro_opts, nsl, frozen%hydrology%q_top,        &
                                 uptake_chk, dtheta_chk, drain_chk, uptk_chk, qface_chk)
      do k = 1_ik, nsl
         eforc_chk%soil_water(k)     = y%theta(k)
         eforc_chk%root_heat_sink(k) = sum(frozen%roots%qloss_frozen(1:n)) * frozen%roots%root_share(k)
         eforc_chk%w_flux(k)         = -qface_chk(k)
      end do
      eforc_chk%root_heat_sink(nsl) = eforc_chk%root_heat_sink(nsl)                                  &
                                    + f%drainage_rate * internal_energy_liquid(frozen%hydrology%t_bot)
      call soil_energy_time_deriv(se_chk, eforc_chk, frozen%params%therm, frozen%params%soil, frozen%params%energy_opts, dedt_chk)
      worst = maxval(abs(f%dedt(1:nsl) - dedt_chk(1:nsl)))
      call check_true('column_derivs wires the soil-heat tendency correctly', worst < 1.0e-9_wp, worst)
   end subroutine test_column_assembler

   !----- 10. march the whole column with explicit RK4 over column_derivs (the reference oracle):   !
   !          it stays finite + physical, and self-converges (dt vs dt/2 agree to RK4 order). -------!
   subroutine test_rk4_march()
      type(column_state_t)  :: y, y1, y2, ytmp
      type(column_frozen_t) :: frozen
      real(wp)    :: dt, tcas, tsoil, dcas, dmass
      integer(ik) :: step, n, nsl, nstep, i, k
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rk4_march:'
      call make_column(y, frozen, n, nsl)

      !----- (a) coarse march at dt = 10 s (< the ~47 s RK4 limit of the 17 s hydraulic mode). ----!
      dt = 10.0_wp ; nstep = 180_ik            ! 30 minutes
      call copy_state(y, y1, n)
      physical = .true.
      do step = 1_ik, nstep
         call rk4_column_step(y1, frozen, n, nsl, dt, ytmp)
         call copy_state(ytmp, y1, n)
         tcas = cas_temp_of_enthalpy(y1%cas_enthalpy, y1%cas_shv)
         physical = physical .and. tcas > 260.0_wp .and. tcas < 330.0_wp .and. y1%cas_shv > 0.0_wp
         do k = 1_ik, nsl
            physical = physical .and. y1%theta(k) > 0.0_wp .and. y1%theta(k) < 0.6_wp
         end do
         do i = 1_ik, n
            physical = physical .and. y1%leaf_water_mass(i) > 0.0_wp .and. y1%wood_water_mass(i) > 0.0_wp
         end do
      end do
      call check_true('RK4 march stays finite + physical (30 min @ 10 s)', physical, dt)

      !----- (b) self-convergence: dt = 8 s vs dt = 4 s over the SAME 8 min window agree tightly. --!
      call copy_state(y, y1, n) ; call copy_state(y, y2, n)
      do step = 1_ik, 60_ik                    ! 60 * 8 s = 8 min
         call rk4_column_step(y1, frozen, n, nsl, 8.0_wp, ytmp) ; call copy_state(ytmp, y1, n)
      end do
      do step = 1_ik, 120_ik                   ! 120 * 4 s = 8 min
         call rk4_column_step(y2, frozen, n, nsl, 4.0_wp, ytmp) ; call copy_state(ytmp, y2, n)
      end do
      tcas  = cas_temp_of_enthalpy(y1%cas_enthalpy, y1%cas_shv)
      tsoil = cas_temp_of_enthalpy(y2%cas_enthalpy, y2%cas_shv)
      dcas  = abs(tcas - tsoil)
      dmass = max(maxval(abs(y1%leaf_water_mass(1:n) - y2%leaf_water_mass(1:n))),                 &
                  maxval(abs(y1%wood_water_mass(1:n) - y2%wood_water_mass(1:n))))
      call check_true('RK4 self-convergence: CAS temp dt=8 vs dt=4 agree < 1e-3 K', dcas < 1.0e-3_wp, dcas)
      call check_true('RK4 self-convergence: mass dt=8 vs dt=4 agree < 1e-6 kg/plant', dmass < 1.0e-6_wp, dmass)
   end subroutine test_rk4_march

   !----- 11. IMEX-Euler: stable at dt = 900 s (where RK4 blows up), and agrees with the RK4        !
   !          oracle in the small-dt limit (both integrate the same column_derivs RHS). ------------!
   subroutine test_imex_euler()
      type(column_state_t)  :: y, yi, yr, ytmp
      type(column_frozen_t) :: frozen
      real(wp)    :: tcas, dcas, dtheta, dmass
      integer(ik) :: step, n, nsl, k, i
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_imex_euler:'
      call make_column(y, frozen, n, nsl)
      !----- Well-ventilated canopy (stronger CAS<->atm exchange) so the frozen-per-step source at   !
      !      the full 900 s does not out-run venting into supersaturation -- the large-dt operator-   !
      !      split coupling error the P2 arrowhead removes; here we exercise L-stability + physical.  !
      frozen%cas%g_atm_heat = frozen%cas%g_atm_heat * 4.0_wp ; frozen%cas%g_atm_vapour = frozen%cas%g_atm_vapour * 4.0_wp
      frozen%cas%g_atm_co2 = frozen%cas%g_atm_co2 * 4.0_wp

      !----- (a) STABLE + physical at the full fast timestep dt = 900 s over a 6-hour march. --------!
      call copy_state(y, yi, n)
      physical = .true.
      do step = 1_ik, 24_ik                    ! 24 * 900 s = 6 h
         call imex_euler_column_step(yi, frozen, n, nsl, 900.0_wp, ytmp)
         call copy_state(ytmp, yi, n)
         tcas = cas_temp_of_enthalpy(yi%cas_enthalpy, yi%cas_shv)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. yi%cas_shv > 0.0_wp
         do k = 1_ik, nsl
            physical = physical .and. yi%theta(k) > 0.0_wp .and. yi%theta(k) < 0.6_wp
         end do
         do i = 1_ik, n
            physical = physical .and. yi%leaf_water_mass(i) > 0.0_wp
         end do
      end do
      call check_true('IMEX-Euler stable + physical at dt = 900 s (6 h)', physical, tcas)

      !----- (b) cross-validation: IMEX-Euler vs the explicit RK4 oracle at dt = 4 s over 4 min agree to  !
      !          first order (both solve the same RHS; O(dt) split-vs-coupled difference). Soil water is   !
      !          OPERATOR-SPLIT out of the ARK stepper (theta frozen across the stages -- the scratch      !
      !          advance_soil_water_column is the sole authority at the column_fast_step level), so the oracle !
      !          is run with freeze_theta=.true.: same reduced system -> theta is trivially equal and the  !
      !          CAS/soil-energy core agrees tightly. Mass stays operator-split in BOTH (in-vector oracle   !
      !          vs advance_water_mass_full) -> a small O(dt) difference (both are exact closed-form Euler  !
      !          now, so this should be tighter than psi's old O(dt) split-vs-exact-exp gap). --------------!
      call copy_state(y, yi, n) ; call copy_state(y, yr, n)
      do step = 1_ik, 60_ik                    ! 60 * 4 s = 4 min
         call imex_euler_column_step(yi, frozen, n, nsl, 4.0_wp, ytmp)              ; call copy_state(ytmp, yi, n)
         call rk4_column_step(yr, frozen, n, nsl, 4.0_wp, ytmp, freeze_theta=.true.); call copy_state(ytmp, yr, n)
      end do
      dcas   = abs(cas_temp_of_enthalpy(yi%cas_enthalpy, yi%cas_shv) - cas_temp_of_enthalpy(yr%cas_enthalpy, yr%cas_shv))
      dtheta = maxval(abs(yi%theta(1:nsl) - yr%theta(1:nsl)))
      dmass  = max(maxval(abs(yi%leaf_water_mass(1:n) - yr%leaf_water_mass(1:n))),                &
                  maxval(abs(yi%wood_water_mass(1:n) - yr%wood_water_mass(1:n))))
      call check_true('IMEX-Euler ~ RK4 oracle (core, theta frozen): CAS temp agree < 5e-2 K at dt=4 s', dcas < 5.0e-2_wp, dcas)
      call check_true('IMEX-Euler ~ RK4 oracle: theta frozen in both (dtheta == 0)', dtheta < 1.0e-12_wp, dtheta)
      call check_true('IMEX-Euler ~ RK4 oracle: mass (operator-split both) agree < 1e-4 kg/plant', dmass < 1.0e-4_wp, dmass)
   end subroutine test_imex_euler

   !----- 12. the leaf<->CAS Picard coupling (niter>1) removes the large-dt over-humidification that  !
   !          the uncoupled baseline (niter=1) suffers under a harsh constant-noon forcing at 900 s.  !
   subroutine test_imex_coupled()
      type(column_state_t)  :: y, yb, yc, ytmp
      type(column_frozen_t) :: frozen
      real(wp)    :: tcas_base, tcas_coup, qsat_base, qsat_coup
      integer(ik) :: step, n, nsl
      logical     :: base_diverged
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_imex_coupled:'
      !----- HARSH forcing: constant noon, modest venting (the case that over-humidifies at 900 s). -!
      call make_column(y, frozen, n, nsl)      ! no ventilation boost -> baseline over-humidifies

      !----- STOP the baseline the moment it has demonstrably collapsed.  It is DELIBERATELY unstable,  !
      !      and once it has left the physical band there is nothing further to learn from it -- but      !
      !      marching it on regardless takes it to -3.0e79 J/kg by step 10 and then into IEEE overflow,    !
      !      which is a landmine rather than a test: on a Release build the assertion below silently       !
      !      compared against NaN, and on the -fpe0 Debug build the overflow TRAPS inside the integrator    !
      !      before any assertion runs.  Which of those happens depends on rounding in an already-          !
      !      diverged trajectory, so ANY perturbation elsewhere in the surface kernel can flip it (one       !
      !      did: a change that is bit-identical here through step 10).  Detecting the collapse and          !
      !      stopping states the intent exactly and removes the dependence on how the collapse ends. -------!
      call copy_state(y, yb, n) ; call copy_state(y, yc, n)
      base_diverged = .false.
      do step = 1_ik, 12_ik                 ! 12 * 900 s = 3 h of constant noon
         if (.not. base_diverged) then
            call imex_euler_column_step(yb, frozen, n, nsl, 900.0_wp, ytmp)                     ! niter=1 (baseline)
            call copy_state(ytmp, yb, n)
            tcas_base = cas_temp_of_enthalpy(yb%cas_enthalpy, yb%cas_shv)
            base_diverged = .not. (tcas_base > 270.0_wp .and. tcas_base < 325.0_wp)
         end if
         call imex_euler_column_step(yc, frozen, n, nsl, 900.0_wp, ytmp, niter=12_ik)  ! coupled
         call copy_state(ytmp, yc, n)
      end do
      tcas_base = cas_temp_of_enthalpy(yb%cas_enthalpy, yb%cas_shv)
      tcas_coup = cas_temp_of_enthalpy(yc%cas_enthalpy, yc%cas_shv)
      qsat_base = sat_specific_humidity(tcas_base, frozen%cas%press)
      qsat_coup = sat_specific_humidity(tcas_coup, frozen%cas%press)
      !----- coupled stays physical (sub-saturated, physical temperature); baseline collapses. ------!
      call check_true('coupled CAS stays physical (270-325 K) at 900 s constant noon',            &
                      tcas_coup > 270.0_wp .and. tcas_coup < 325.0_wp, tcas_coup)
      call check_true('coupled CAS stays sub-saturated (shv <= qsat)', yc%cas_shv <= qsat_coup + 1.0e-4_wp, &
                      yc%cas_shv - qsat_coup)
      !----- COUPLING MATTERS: assert the baseline LEAVES the physical band, rather than comparing    !
      !      the two temperatures.  The baseline is deliberately unstable -- by step 9 its specific     !
      !      humidity is already NEGATIVE (-9.7e-2 kg/kg) and by step 10 its enthalpy is -3.0e79 J/kg,  !
      !      so tcas_base is not a number any arithmetic comparison can rest on.  Whether that garbage  !
      !      reaches IEEE overflow (-> NaN, which makes every comparison .false.) within these 12 steps  !
      !      depends on rounding in an already-diverged trajectory, so the previous form                 !
      !      `tcas_coup > tcas_base + 1` was green by accident of overflow TIMING, not by physics: an    !
      !      unrelated, bit-identical-until-step-10 change to the surface kernel flipped it to NaN.      !
      !      The band test states the intent directly (the coupled run stays physical where the          !
      !      uncoupled one collapses) and is immune to how the collapse ends. -------------------------!
      call check_true('coupling matters: the uncoupled baseline leaves the physical band',         &
                      base_diverged, tcas_base)
   end subroutine test_imex_coupled

   !----- 12b. the arrowhead Newton surface solve: robust near saturation (clamp + line search), and  !
   !           faithful to the RK4 oracle at small dt (converges to the coupled surface solution). ----!
   subroutine test_arrowhead()
      type(column_state_t)  :: y, yi, yr, ytmp
      type(column_frozen_t) :: frozen
      real(wp)    :: tcas, qsat, dcas, dtheta, worst_super
      integer(ik) :: step, n, nsl, k
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_arrowhead:'

      !----- (a) START the CAS at 99% saturation under harsh noon: exercises the supersaturation      !
      !          clamp + line search. Newton (niter>1) must keep it physical + sub-saturated at 900 s. !
      call make_column(y, frozen, n, nsl)
      tcas = 297.0_wp ; qsat = sat_specific_humidity(tcas, frozen%cas%press)
      y%cas_enthalpy = cas_enthalpy_of_temp(tcas, 0.99_wp*qsat) ; y%cas_shv = 0.99_wp*qsat
      call copy_state(y, yi, n)
      physical = .true. ; worst_super = -1.0_wp
      do step = 1_ik, 12_ik
         call imex_euler_column_step(yi, frozen, n, nsl, 900.0_wp, ytmp, niter=8_ik)
         call copy_state(ytmp, yi, n)
         tcas = cas_temp_of_enthalpy(yi%cas_enthalpy, yi%cas_shv)
         qsat = sat_specific_humidity(tcas, frozen%cas%press)
         worst_super = max(worst_super, yi%cas_shv - qsat)          ! must stay <= 0 (sub-saturated)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp
      end do
      call check_true('Newton stays sub-saturated from a 99%-saturated start (900 s)', worst_super <= 1.0e-4_wp, worst_super)
      call check_true('Newton stays physical (270-325 K) near saturation', physical, tcas)

      !----- (b) Newton (niter=8) matches the RK4 oracle at small dt (converges to the coupled soln). -!
      call make_column(y, frozen, n, nsl)
      call copy_state(y, yi, n) ; call copy_state(y, yr, n)
      do step = 1_ik, 60_ik                     ! 60 * 4 s = 4 min
         call imex_euler_column_step(yi, frozen, n, nsl, 4.0_wp, ytmp, niter=8_ik)   ; call copy_state(ytmp, yi, n)
         call rk4_column_step(yr, frozen, n, nsl, 4.0_wp, ytmp, freeze_theta=.true.) ; call copy_state(ytmp, yr, n)
      end do
      dcas   = abs(cas_temp_of_enthalpy(yi%cas_enthalpy, yi%cas_shv) - cas_temp_of_enthalpy(yr%cas_enthalpy, yr%cas_shv))
      dtheta = maxval(abs(yi%theta(1:nsl) - yr%theta(1:nsl)))
      call check_true('Newton (niter=8) ~ RK4 oracle (core, theta frozen): CAS temp < 5e-2 K at dt=4 s', dcas < 5.0e-2_wp, dcas)
      call check_true('Newton (niter=8) ~ RK4 oracle: theta frozen in both (dtheta == 0)', dtheta < 1.0e-12_wp, dtheta)
   end subroutine test_arrowhead

   !----- 13. the step-doubling adaptive controller: tighter rtol takes more steps and both agree     !
   !          with a fine fixed reference -- the P3 adaptive time-stepping contract. ---------------!
   subroutine test_adaptive_march()
      type(column_state_t)  :: y, y1, y2, yr, ytmp
      type(column_frozen_t) :: frozen
      real(wp)    :: tc1, tc2, tcr
      integer(ik) :: step, n, nsl, ns1, ns2, nr1, nr2
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_adaptive_march:'
      call make_column(y, frozen, n, nsl)
      frozen%cas%g_atm_heat   = frozen%cas%g_atm_heat   * 4.0_wp
      frozen%cas%g_atm_vapour = frozen%cas%g_atm_vapour * 4.0_wp
      frozen%cas%g_atm_co2    = frozen%cas%g_atm_co2    * 4.0_wp
      !----- start the CAS cool + dry so there is a real transient to resolve adaptively. -----------!
      y%cas_enthalpy = cas_enthalpy_of_temp(290.0_wp, 0.009_wp) ; y%cas_shv = 0.009_wp

      call adaptive_imex_march(y, frozen, n, nsl, 1800.0_wp, 1.0e-3_wp, 50.0_wp, y1, ns1, nr1)
      call adaptive_imex_march(y, frozen, n, nsl, 1800.0_wp, 1.0e-5_wp, 50.0_wp, y2, ns2, nr2)
      !----- fine fixed reference (dt = 2 s, coupled): 900 steps. ----------------------------------!
      call copy_state(y, yr, n)
      do step = 1_ik, 900_ik
         call imex_euler_column_step(yr, frozen, n, nsl, 2.0_wp, ytmp, niter=8_ik)
         call copy_state(ytmp, yr, n)
      end do
      tc1 = cas_temp_of_enthalpy(y1%cas_enthalpy, y1%cas_shv)
      tc2 = cas_temp_of_enthalpy(y2%cas_enthalpy, y2%cas_shv)
      tcr = cas_temp_of_enthalpy(yr%cas_enthalpy, yr%cas_shv)
      print '(a,i0,a,i0,a,i0)', '   adaptive steps: rtol=1e-3 -> ', ns1, ' , rtol=1e-5 -> ', ns2, ' ; fixed dt=2s -> 900'
      call check_true('adaptive rtol=1e-3 agrees with fine reference (< 0.5 K)', abs(tc1 - tcr) < 0.5_wp, tc1 - tcr)
      call check_true('adaptive rtol=1e-5 agrees with fine reference (< 0.5 K)', abs(tc2 - tcr) < 0.5_wp, tc2 - tcr)
      call check_true('tighter rtol takes more steps (controller responds to tol)', ns2 > ns1, real(ns2 - ns1, wp))
      call check_true('adaptive is cheaper than the fixed fine march (steps < 900)', ns1 < 900_ik, real(ns1, wp))
      call check_true('the march made real progress (steps > 0)', ns1 > 0_ik, real(ns1, wp))
   end subroutine test_adaptive_march

   !----- soil-top temperature diagnosed from a column_state_t (helper for the ARK2 order test). ---!
   function soil_top_temp(y, frozen) result(t)
      type(column_state_t),  intent(in) :: y
      type(column_frozen_t), intent(in) :: frozen
      real(wp) :: t, fl
      call internal_energy_to_temp(y%soil_energy(1), y%theta(1)*rho_h2o, frozen%params%therm%soil_dry_heat_capacity(1), t, fl)
   end function soil_top_temp

   !----- march the ARK2 fixed-step from y for nstep steps of dt (embedded error discarded). --------!
   subroutine march_ark2(y0, frozen, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call ark2_column_step(y, frozen, n, nsl, dt, ytmp, yerr, niter=8_ik)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_ark2

   !----- march coupled IMEX-Euler (niter=8 Newton) fixed-step (for the ARK2-vs-1st-order comparison). !
   subroutine march_imex(y0, frozen, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call imex_euler_column_step(y, frozen, n, nsl, dt, ytmp, niter=8_ik)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_imex

   !----- march the Cash-Karp RK45 fixed-step (5th-order commit, embedded 4th discarded) -- for the   !
   !      order-of-accuracy self-convergence study below. ------------------------------------------!
   subroutine march_rk45(y0, frozen, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp, yerr
      real(wp) :: w_out, e_in, e_out
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call rk45_column_step(y, frozen, n, nsl, dt, ytmp, yerr, w_out, e_in, e_out)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_rk45

   !----- 14. ARK2 (ARS(2,2,2)): 2nd-order on the soil-top temperature (differential var); mass stays !
   !          physical at production dt (FATAL-1 guard); embedded estimate bounded; stiff-stable. -----!
   subroutine test_ark2()
      type(column_state_t)  :: y, yref, y1, y2, y4, ynew, yerr, ytmp
      type(column_frozen_t) :: frozen
      type(error_control_t) :: ec
      real(wp)    :: e1, e2, e4, p_lo, p_hi, tref, errnorm, tcas
      integer(ik) :: n, nsl, step, ns, nr, ns2, nr2
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_ark2:'
      call make_column(y, frozen, n, nsl)

      !----- (a) ORDER-2 self-convergence on a CLEANLY-integrated tableau variable: CAS CO2 (a         !
      !          decoupled affine scalar ODE -- no operator split, so the ARS(2,2,2) order is exposed).!
      !          reference = ARK2 at dt = 1800/64; solutions at dt = 450, 225, 112.5. ----------------!
      call march_ark2(y, frozen, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = yref%cas_co2
      call march_ark2(y, frozen, n, nsl, 450.0_wp,   4_ik,  y1) ; e1 = abs(y1%cas_co2 - tref)
      call march_ark2(y, frozen, n, nsl, 225.0_wp,   8_ik,  y2) ; e2 = abs(y2%cas_co2 - tref)
      call march_ark2(y, frozen, n, nsl, 112.5_wp,  16_ik,  y4) ; e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)', '   CAS CO2 errors: ', e1, ' / ', e2, &
            ' / ', e4, ' ; observed order p = ', p_lo, ' , ', p_hi
      call check_true('ARK2 tableau is 2nd order on the decoupled CO2 twin (p >= 1.9)', &
                      p_lo >= 1.9_wp .and. p_hi >= 1.9_wp, p_hi)

      !----- (a2) On the operator-split-COUPLED soil-top temperature, the full 3x3 arrowhead is        !
      !           deferred (spec) so the effective order is ~1.2, but ARK2 is still MORE ACCURATE than  !
      !           the 1st-order IMEX-Euler at the same dt. -------------------------------------------!
      call march_ark2(y, frozen, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = soil_top_temp(yref, frozen)
      call march_ark2(y, frozen, n, nsl, 225.0_wp,  8_ik, y1) ; e1 = abs(soil_top_temp(y1, frozen) - tref)
      call march_imex(y, frozen, n, nsl, 225.0_wp,  8_ik, y2) ; e2 = abs(soil_top_temp(y2, frozen) - tref)
      call check_true('ARK2 beats IMEX-Euler on the coupled soil-top temperature (same dt)', e1 < e2, e2 - e1)

      !----- (b) FATAL-1 guard: mass stays physical (finite, positive) at production dt=900. ---------!
      call copy_state(y, y1, n) ; physical = .true.
      do step = 1_ik, 24_ik
         call ark2_column_step(y1, frozen, n, nsl, 900.0_wp, ytmp, yerr, niter=8_ik)
         call copy_state(ytmp, y1, n)
         physical = physical .and. all(y1%leaf_water_mass(1:n) == y1%leaf_water_mass(1:n)) .and.   &
                    all(y1%wood_water_mass(1:n) == y1%wood_water_mass(1:n)) .and.                  &
                    minval(y1%leaf_water_mass(1:n)) > 0.0_wp .and. minval(y1%wood_water_mass(1:n)) > 0.0_wp
      end do
      call check_true('ARK2 mass stays physical at dt=900 (FATAL-1 split-out guard)', physical,     &
                      minval(y1%leaf_water_mass(1:n)))

      !----- (c) embedded error estimate is bounded (not detonating -- the pre-fix failure mode). ----!
      call ark2_column_step(y, frozen, n, nsl, 900.0_wp, ynew, yerr, niter=8_ik)
      call state_err_norm(ynew, yerr, y, n, nsl, errnorm)
      call check_true('ARK2 embedded estimate bounded (WRMS < 50) at dt=900', errnorm < 50.0_wp, errnorm)

      !----- (d) stiffness: 24 h adaptive-ARK march stays physical + bounded. -----------------------!
      ec = default_error_control(1.0e-3_wp)
      call adaptive_ark_march(y, frozen, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ns, nr)
      tcas = cas_temp_of_enthalpy(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   adaptive-ARK 24 h: steps = ', ns, ' , rejects = ', nr
      call check_true('adaptive-ARK 24 h stays bounded (280 < tcas < 320 K)', tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)

      !----- (e) ERROR-CONTROL: the PI controller (goal a) marches the SAME stiff window and stays     !
      !      physical -- proves the meds_fast_control PI path is wired + functional. On a multi-substep !
      !      march it takes a different (typically smoother) step sequence than the I-controller.       !
      ec%controller = CTRL_PI
      call adaptive_ark_march(y, frozen, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ns2, nr2)
      tcas = cas_temp_of_enthalpy(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   PI-controller  24 h: steps = ', ns2, ' , rejects = ', nr2
      call check_true('PI-controller ARK 24 h stays bounded (280 < tcas < 320 K)', tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)
   end subroutine test_ark2

   !----- 15. RK45 (Cash-Karp 5(4)): order-of-accuracy self-convergence (design doc sec 8 gate 1).      !
   !          Unlike ARK2's split-degraded ~1.2 order on the coupled soil-top temperature (test_ark2's    !
   !          part a2 -- soil water/mass are operator-split OUT of the ESDIRK stages), RK45 has NO         !
   !          operator split at all (design doc sec 6/P2), so the FULL 5th-order tableau accuracy          !
   !          should show up on a genuinely-coupled variable too, not just a decoupled scalar. Same         !
   !          reference-vs-halvings technique as test_ark2's part (a): observed order                       !
   !          p = log2(e(dt)/e(dt/2)) should approach 5 as dt shrinks (some slack for pre-asymptotic         !
   !          effects: p >= 4.5). ---------------------------------------------------------------------!
   subroutine test_rk45_order()
      type(column_state_t)  :: y, yref, y1, y2, y4, y8
      type(column_frozen_t) :: frozen
      real(wp)    :: e1, e2, e4, e8, p_lo, p_hi, p_fine, tref
      integer(ik) :: n, nsl
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rk45_order:'
      call make_column(y, frozen, n, nsl)

      !----- A 5th-order method's error ~ C*h^5 collapses to double-precision noise (~1e-13) at MUCH      !
      !      coarser h than a 2nd-order one (test_ark2's own dt=200/8..32 s range) -- h^5 falls off so    !
      !      fast that h must stay large enough for the TRUNCATION term to still dominate roundoff at      !
      !      the FINEST test point. Window T=2400 s (8 * the coarsest dt below), test points dt=300/       !
      !      150/75 s (8/16/32 steps; 300 s stays under the ~360 s single-step stability estimate,          !
      !      design doc sec 6 -- this is a short fixed-step accuracy comparison, not a stability march),     !
      !      reference dt=2400/1024 s (>=32x finer than the finest test point). ---------------------------!
      !----- (a) CLEAN scalar: CAS CO2 (decoupled affine ODE, mirrors test_ark2's part a) -- baseline   !
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

      !----- (b) COUPLED: soil-top temperature (the SAME variable test_ark2 part a2 shows degraded to  !
      !          ~1.2 order under ARK's operator split). RK45 integrates soil energy/water/mass ALL      !
      !          genuinely (no split), so this should ALSO be ~5th order -- the key accuracy advantage    !
      !          the design doc's "no operator split at all" claim (sec 6) predicts. -------------------!
      !----- FOUR test points here, not three, and the assertion is on the two FINEST pairs. The soil-  !
      !      top temperature is the most strongly coupled variable in the column, so its coarse-h error  !
      !      carries a visible h^6 contamination and p measured across the coarsest pair is simply not   !
      !      asymptotic. That was true before #78 item 3 as well, just on the flattering side: with the  !
      !      interior advective enthalpy FROZEN the coupling was weaker and the coarse pair read p =     !
      !      6.05 -- super-convergence, which a 5th-order tableau cannot actually achieve, so it was     !
      !      measuring contamination too. Wiring the faces to the stage's own theta trajectory made the  !
      !      absolute error ~14x SMALLER at every dt tested while moving the asymptotic window finer.    !
      !      Asserting on the finest pairs measures the tableau; asserting on the coarsest measured      !
      !      whichever direction the h^6 term happened to point. --------------------------------------!
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
      call check_true('RK45 stays 5th order on the COUPLED soil-top temperature (p >= 4.5, no split)', &
                      p_hi >= 4.5_wp .and. p_fine >= 4.5_wp, p_fine)
   end subroutine test_rk45_order

   !----- WRMS of the embedded error estimate (mirrors adaptive_ark_march's accept test). ----------!
   subroutine state_err_norm(ynew, yerr, yref, n, nsl, errnorm)
      type(column_state_t), intent(in)  :: ynew, yerr, yref
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

   !----- shared full-column setup (2 cohorts, 10 soil layers, a representative daytime state). ----!
   !----- 15. The explicit RHS is DOMAIN-SAFE for any theta, in or out of [theta_res, theta_sat]      !
   !         (issue #78 item 2). ---------------------------------------------------------------------!
   !                                                                                                 !
   ! Item 2 as filed said RK45's stage states can sit above theta_sat, so soil_psi_from_theta and     !
   ! soil_hydr_cond_from_theta "are then evaluated at effective saturation Se > 1 -- outside the van   !
   ! Genuchten domain", and that the resulting K and psi influence the trajectory. The first half is   !
   ! true and measurable (the k1 stage reads the previous SUB-step's unclamped commit -- clamp_theta    !
   ! runs before stages 2-6 only, and C1 deliberately stopped clamping the committed state). The       !
   ! SECOND half is not: every constitutive kernel on this path clamps its own argument.               !
   !                                                                                                  !
   !   soil_psi_from_theta        se = min(max(se, SE_MIN), 1); Se >= 1 returns psi = 0 outright       !
   !   soil_hydr_cond_from_theta  se = min(max(se, SE_MIN), 1), then max(kcond, K_MIN)                 !
   !   soil_moist_cap_from_psi    psi >= 0 returns C_MIN; max(cap, C_MIN)                              !
   !   soil_thermal_cond          s_r = min(max(theta/theta_sat, SR_FLOOR), 1); Kersten clamped [0,1]  !
   !   internal_energy_to_temp               algebraic in water mass -- no domain to leave                        !
   !                                                                                                  !
   ! face_and_sink adds nothing unguarded: kface is an upstream pick among those K, gface and the      !
   ! psi-limited sink f_wilt_ramp are functions of the already-clamped psi, and cc is unused by the    !
   ! explicit RHS. So dtheta_dt depends on theta ONLY through those two clamped curves, which makes    !
   ! this a BIT-IDENTITY rather than a tolerance: the water tendency at theta_sat + eps must equal     !
   ! the water tendency at theta_sat exactly. That is the assertion below, at eps = 0.05 -- about six   !
   ! times the worst excursion measured anywhere in the suite (7.9e-3 in theta, Se = 1.022, on the      !
   ! 29 mm/h saturated fixture; a month-long forced Ithaca cell never leaves the domain at all, and     !
   ! budget%theta_ood_max now reports it per sub-step so it cannot drift unnoticed).                      !
   !                                                                                                  !
   ! The soil-ENERGY tendency is checked separately (part c) and for a different reason: it is a pure   !
   ! flux divergence over prognostic internal energy, so theta reaches it only through the temperature  !
   ! diagnosis and the self-clamped Kersten number. It SHOULD respond to the excess water -- that       !
   ! water's heat capacity is real -- so the assertion there is finite-and-responsive, not identity.    !
   !-------------------------------------------------------------------------------------------------!
   subroutine test_rhs_domain_safety()
      type(column_state_t)  :: y, y_hi, y_lo, y_sat, y_res
      type(column_frozen_t) :: frozen
      type(column_tend_t)   :: f_sat, f_hi, f_res, f_lo
      real(wp), parameter   :: EPS_OOD = 0.05_wp        ! ~6x the worst excursion measured in the suite
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
         y_sat%theta(k) = frozen%params%soil%theta_sat(k)
         y_hi%theta(k)  = frozen%params%soil%theta_sat(k) + EPS_OOD
         y_res%theta(k) = frozen%params%soil%theta_res(k)
         y_lo%theta(k)  = max(frozen%params%soil%theta_res(k) - EPS_OOD, 0.0_wp)
      end do
      call column_derivs(y_sat, frozen, n, nsl, f_sat)
      call column_derivs(y_hi,  frozen, n, nsl, f_hi)
      call column_derivs(y_res, frozen, n, nsl, f_res)
      call column_derivs(y_lo,  frozen, n, nsl, f_lo)

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

      !----- (b) BIT-IDENTICAL water tendency: the Se clamping inside the curves is COMPLETE, so an   !
      !      oversaturated cell behaves as exactly saturated (psi = 0, K = K_sat) -- which is the     !
      !      physically right answer for one, not merely a safe one. -------------------------------!
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
      !      soil_heat_cap_vol is not called on this path at all (only the implicit sibling uses it). !
      !      theta reaches dedt through exactly two places, both safe:                                !
      !        * the temperature diagnosis internal_energy_to_temp(uext, theta*rho_w, ...) -- more water means    !
      !          more heat capacity means a lower T for the same internal energy, which is correct    !
      !          and is what SHOULD respond to the excess water;                                     !
      !        * soil_thermal_cond's Kersten number, which self-clamps at s_r = 1.                    !
      !                                                                                              !
      !      make_column happens to sit on the mixed-phase MELT PLATEAU, where internal_energy_to_temp returns   !
      !      t_3ple regardless of water mass, so dedt there is bit-identical as well. That is a       !
      !      degenerate case, not the general rule, so re-seed to an all-liquid state before          !
      !      asserting the sensitivity -- otherwise this check would pass for the wrong reason. ------!
      call state_init(y, n, nsl, y_sat) ; call state_init(y, n, nsl, y_hi)
      do k = 1_ik, nsl
         y_sat%theta(k) = frozen%params%soil%theta_sat(k)
         y_hi%theta(k)  = frozen%params%soil%theta_sat(k) + EPS_OOD
         y_sat%soil_energy(k) = temp_to_internal_energy(frozen%params%therm%soil_dry_heat_capacity(k),                      &
                                             y_sat%theta(k)*rho_h2o, 290.0_wp, 1.0_wp)
         y_hi%soil_energy(k)  = y_sat%soil_energy(k)      ! SAME internal energy, more water
      end do
      call column_derivs(y_sat, frozen, n, nsl, f_sat)
      call column_derivs(y_hi,  frozen, n, nsl, f_hi)
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

   subroutine make_column(y, frozen, n, nsl)
      type(column_state_t),  intent(out) :: y
      type(column_frozen_t), intent(out) :: frozen
      integer(ik),           intent(in)  :: n, nsl
      !----- LOCAL hydraulics traits: only needed to seed y%*_water_mass at a representative psi     !
      !      (column_frozen_t no longer carries hydraulics_params/hydraulics_opts -- column_derivs' mass ODE needs      !
      !      only the frozen sapflow/uptake, not the PV-curve/conductance params any more). -----------!
      type(hydro_params_t) :: hp
      integer(ik) :: i, k
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,        &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, frozen%params%soil)
      call build_soil_therm_params(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, frozen%params%therm)
      frozen%params%hydro_opts = soil_opts_t()
      hp%leaf_pi0 = -1.5_wp ; hp%leaf_elastic_mod = 12.0_wp ; hp%leaf_apoplast_frac = 0.30_wp
      hp%leaf_water_sat = 2.0_wp ; hp%wood_pi0 = -1.0_wp ; hp%wood_elastic_mod = 8.0_wp
      hp%wood_apoplast_frac = 0.20_wp ; hp%wood_water_sat = 1.0_wp ; hp%wood_psi50 = -2.0_wp
      hp%wood_kexp = 2.0_wp ; hp%k_plant_max = 6.0e-4_wp ; hp%wood_kmax = 8.0_wp ; hp%vessel_curl = 1.5_wp
      frozen%hydrology%geothermal = 0.0_wp ; frozen%hydrology%q_top = 1.0e-6_wp
      allocate(frozen%roots%root_share(nsl), frozen%plant%nplant(n), frozen%plant%bleaf(n), frozen%plant%bsap(n), &
               frozen%plant%broot(n),                                                                            &
               frozen%plant%sap_area(n))
      allocate(frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n), frozen%roots%qloss_frozen(n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n))
      frozen%roots%root_share(1:nsl) = frozen%params%soil%root_frac(1:nsl)   ! Phase 1: per-layer sink placement
      frozen%plant%nplant = 0.3_wp ; frozen%plant%bleaf = 0.5_wp ; frozen%plant%bsap = 5.0_wp ; frozen%plant%broot = 2.0_wp
      frozen%plant%sap_area = 0.01_wp
      !----- FROZEN sapflow/uptake (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2): a representative,            !
      !      state-independent pair for these RHS/oracle-march tests -- not re-derived from a            !
      !      solve_plant_water_batch pre-pass here (these tests exercise column_derivs/the tableau         !
      !      machinery directly, not build_column_frozen's own pre-pass, which has its own coverage        !
      !      via test_column_ark.f90/test_picard_coupling.f90). -----------------------------------------!
      frozen%plant%sapflow_frozen = 1.0e-4_wp
      frozen%plant%uptake_frozen = 1.0e-4_wp
      frozen%roots%uptake = frozen%plant%uptake_frozen(1)*sum(frozen%plant%nplant(1:n))
      frozen%roots%qloss_frozen = 0.0_wp   ! P2 advective enthalpy: no-op unless populated (see build_column_frozen)
      frozen%film%intercept_leaf = 0.0_wp ; frozen%film%intercept_wood = 0.0_wp   ! P2c canopy water: no-op unless populated
      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n), frozen%tissue%abs_lw(n), &
               frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      frozen%tissue%leaf_hcap_per_dt = 0.0_wp ; frozen%tissue%wood_hcap_per_dt = 0.0_wp
      frozen%tissue%t_leaf0 = 0.0_wp ; frozen%tissue%t_wood0 = 0.0_wp
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
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
      frozen%cas%rho = 1.2_wp ; frozen%cas%press = 101325.0_wp ; frozen%cas%cas_mass_capacity = 1.2_wp*20.0_wp
      frozen%cas%cas_molar_capacity = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*20.0_wp
      frozen%cas%g_atm_heat = 1.2_wp*0.3_wp*0.02_wp ; frozen%cas%g_atm_vapour = frozen%cas%g_atm_heat
      frozen%cas%g_atm_co2 = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*0.3_wp*0.02_wp
      frozen%cas%enthalpy_atm = cas_enthalpy_of_temp(300.0_wp, 0.011_wp) ; frozen%cas%shv_atm = 0.011_wp
      frozen%cas%co2_atm = 400.0_wp ; frozen%cas%nee_biotic = -5.0_wp
      frozen%ground%abs_sw_ground = 60.0_wp ; frozen%ground%abs_lw_ground = -10.0_wp
      frozen%ground%ggnet = 0.02_wp ; frozen%ground%soil_evap = 2.0e-5_wp
      y%cas_enthalpy = cas_enthalpy_of_temp(297.0_wp, 0.012_wp) ; y%cas_shv = 0.012_wp ; y%cas_co2 = 410.0_wp
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water = 0.0_wp ; y%wood_surf_water = 0.0_wp   ! P2c canopy water: no-op unless populated
      do k = 1_ik, nsl
         y%soil_energy(k) = temp_to_internal_energy(frozen%params%therm%soil_dry_heat_capacity(k), 0.30_wp*rho_h2o,   &
                            296.0_wp - 0.4_wp*real(k-1_ik, wp), 1.0_wp)
         y%theta(k) = 0.30_wp
      end do
      !----- seed mass at the SAME representative (leaf psi=-1.0, wood psi=-0.5 MPa) point the old   !
      !      psi-based fixture used, via the forward water_content map (hp above). -------------------!
      do i = 1_ik, n
         y%leaf_water_mass(i) = water_content(-1.0_wp, hp%leaf_pi0, hp%leaf_elastic_mod,          &
              hp%leaf_apoplast_frac, hp%leaf_water_sat, frozen%plant%bleaf(i))
         y%wood_water_mass(i) = water_content(-0.5_wp, hp%wood_pi0, hp%wood_elastic_mod,          &
              hp%wood_apoplast_frac, hp%wood_water_sat, frozen%plant%bsap(i) + frozen%plant%broot(i))
      end do
   end subroutine make_column

   subroutine copy_state(a, b, n)
      type(column_state_t), intent(in)  :: a
      type(column_state_t), intent(out) :: b
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
   !      what column_derivs does internally) so the standalone CAS comparison lines up. -----------!
   function tground_of(y, frozen) result(tg)
      type(column_state_t),   intent(in) :: y
      type(column_frozen_t),  intent(in) :: frozen
      real(wp) :: tg, fl
      call internal_energy_to_temp(y%soil_energy(1), y%theta(1)*rho_h2o, frozen%params%therm%soil_dry_heat_capacity(1), tg, fl)
   end function tground_of

   logical function ieee_ok(x)
      real(wp), intent(in) :: x
      ieee_ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function ieee_ok

end program test_column_derivs
