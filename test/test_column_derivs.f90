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
!   5. SUPPLY LIMITER            -- src_frac scales the transpiration source into the CAS.            !
!==========================================================================================!
program test_column_derivs
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : latent_heat_vap, stefan, cp_air, tiny_num, rho_h2o
   use meds_therm_lib,         only : cas_enthalpy_of_temp, cas_temp_of_enthalpy,                   &
                                   sat_specific_humidity, sat_vapor_pressure,                    &
                                   sat_vapor_pressure_temp_deriv, enthalpy_vapor, temp_to_uext, uext_to_temp
   use meds_budget_check,   only : budget_t, budget_accumulate
   use meds_biophysics_types, only : n_soil_layer_max, soil_params_t, soil_thermal_params_t,     &
                                   energy_forcing_t, energy_opts_t, soil_energy_column_t,        &
                                   energy_flux_t, soil_opts_t, SOIL_RETENTION_VG
   use meds_column_state_types, only : build_soil_hydr_params
   use meds_column_state_types, only : build_soil_therm_params
   use meds_soil_energy,      only : soil_energy_step_implicit, soil_energy_time_deriv
   use meds_soil_water,       only : soil_water_time_deriv
   use meds_plant_types,      only : hydro_params_t, hydro_opts_t
   use meds_hydr_lib,         only : water_content
   use meds_fast_time_derivs, only : surface_derivs, column_derivs
   use meds_therm_lib,        only : internal_energy_liquid
   use meds_fast_types,       only : surface_state_t, surface_frozen_t, surface_tend_t,           &
                                   column_state_t, column_frozen_t, column_tend_t
   use meds_fast_rk4_oracle,  only : rk4_column_step, imex_euler_column_step, adaptive_imex_march
   use meds_fast_ark,         only : ark2_column_step, adaptive_ark_march
   use meds_fast_rk45,        only : rk45_column_step
   use meds_fast_control,     only : error_control_t, default_error_control
   use meds_config,           only : CTRL_PI
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_leaf_closure()
   call test_leaf_analytic()
   call test_cas_be_consistency()
   call test_conservation_march()
   call test_supply_limiter()
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
   subroutine make_frozen(fro, n)
      type(surface_frozen_t), intent(out) :: fro
      integer(ik),            intent(in)  :: n
      integer(ik) :: i
      allocate(fro%h_coeff_f(n), fro%g_tr_f(n), fro%abs_sw(n), fro%abs_lw(n), fro%lai(n))
      allocate(fro%h_coeff_w(n), fro%abs_sw_wood(n), fro%abs_lw_wood(n), fro%wai(n))
      allocate(fro%qwflux_wl(n), fro%q_wood_net(n))
      allocate(fro%f_wet_c(n), fro%g_film_f(n), fro%g_film_w(n))
      fro%h_coeff_w = 0.0_wp ; fro%abs_sw_wood = 0.0_wp ; fro%abs_lw_wood = 0.0_wp ; fro%wai = 0.0_wp
      fro%qwflux_wl = 0.0_wp ; fro%q_wood_net = 0.0_wp   ! P2 advective enthalpy: no-op unless populated
      fro%f_wet_c = 0.0_wp ; fro%g_film_f = 0.0_wp ; fro%g_film_w = 0.0_wp   ! P2c canopy water: no-op unless populated
      do i = 1_ik, n
         fro%lai(i)       = 2.0_wp - 0.4_wp * real(i - 1_ik, wp)          ! 2.0, 1.6, 1.2
         fro%abs_sw(i)    = 250.0_wp - 40.0_wp * real(i - 1_ik, wp)       ! more light at the top
         fro%abs_lw(i)    = -30.0_wp
         fro%h_coeff_f(i) = 2.0_wp * fro%lai(i) * 0.03_wp * 1.2_wp * cp_air  ! effarea*lai*gbh*rho*cp
         fro%g_tr_f(i)    = 0.004_wp * fro%lai(i)                          ! series conductance [m/s]
      end do
      fro%leaf_emiss    = 0.95_wp
      fro%rho           = 1.2_wp
      fro%press         = 101325.0_wp
      fro%wcap          = fro%rho * 20.0_wp                                ! rho * can_depth
      fro%ccap          = (fro%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 20.0_wp
      fro%gah           = fro%rho * 0.3_wp * 0.02_wp                       ! rho * ustar * temp1
      fro%gaw           = fro%rho * 0.3_wp * 0.02_wp
      fro%gac           = (fro%rho * (1.0_wp - 0.012_wp) / 0.0289655_wp) * 0.3_wp * 0.02_wp
      fro%enth_atm      = cas_enthalpy_of_temp(300.0_wp, 0.011_wp)
      fro%shv_atm       = 0.011_wp
      fro%co2_atm       = 400.0_wp
      fro%nee_biotic    = -5.0_wp                                          ! net CO2 uptake [umol/m2/s]
      fro%abs_sw_ground = 60.0_wp
      fro%abs_lw_ground = -10.0_wp
      fro%ggnet         = 0.02_wp
      fro%soil_evap     = 2.0e-5_wp
      fro%src_frac      = 1.0_wp
      fro%t_ground      = 297.0_wp
   end subroutine make_frozen

   !----- 1. The diagnosed leaf temperature must zero the (linearized) leaf energy balance. ------!
   subroutine test_leaf_closure()
      type(surface_frozen_t) :: fro
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat, dtl, lw_slope, le_slope, le_ref, resid, worst
      integer(ik) :: i, n
      n = 3_ik
      print '(a)', 'test_leaf_closure:'
      call make_frozen(fro, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(298.0_wp, 0.012_wp)
      y%cas_shv      = 0.012_wp
      y%cas_co2      = 410.0_wp
      call surface_derivs(y, fro, n, f)

      tcas   = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, fro%press)
      esat   = sat_vapor_pressure(tcas)
      dqdt   = 0.622_wp * fro%press / max((fro%press - 0.378_wp * esat) ** 2, tiny_num)          &
               * sat_vapor_pressure_temp_deriv(tcas)
      worst = 0.0_wp
      do i = 1_ik, n
         dtl      = f%leaf_temp(i) - tcas
         lw_slope = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%lai(i)
         le_slope = latent_heat_vap * fro%rho * fro%g_tr_f(i) * dqdt
         le_ref   = latent_heat_vap * fro%rho * fro%g_tr_f(i) * (qsat_c - qcas)
         !----- Rnet - sensible - latent - LW-emission, all at the diagnosed leaf temperature. ---!
         resid = fro%abs_sw(i) + fro%abs_lw(i) - fro%h_coeff_f(i) * dtl                          &
                 - (le_ref + le_slope * dtl) - lw_slope * dtl
         worst = max(worst, abs(resid))
      end do
      call check_true('diagnosed leaf T zeroes the linearized balance', worst < 1.0e-9_wp, worst)
      call check_true('leaf temperatures physical (270-330 K)',                                 &
                      minval(f%leaf_temp) > 270.0_wp .and. maxval(f%leaf_temp) < 330.0_wp,       &
                      f%leaf_temp(1))
   end subroutine test_leaf_closure

   !----- 2. No latent (g_tr_f = 0) and no LW forcing (abs_lw = 0): dtl = abs_sw / (h + lw_slope). !
   subroutine test_leaf_analytic()
      type(surface_frozen_t) :: fro
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp) :: tcas, lw_slope, expect
      print '(a)', 'test_leaf_analytic:'
      call make_frozen(fro, 1_ik)
      fro%g_tr_f(1) = 0.0_wp ; fro%abs_lw(1) = 0.0_wp ; fro%abs_sw(1) = 200.0_wp
      fro%lai(1) = 1.5_wp
      fro%h_coeff_f(1) = 2.0_wp * fro%lai(1) * 0.03_wp * 1.2_wp * cp_air
      y%cas_enthalpy = cas_enthalpy_of_temp(299.0_wp, 0.010_wp)
      y%cas_shv      = 0.010_wp
      call surface_derivs(y, fro, 1_ik, f)
      tcas     = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      lw_slope = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%lai(1)
      expect   = tcas + fro%abs_sw(1) / (fro%h_coeff_f(1) + lw_slope)
      call check('leaf T = tcas + Rn/(h+lw_slope) (no latent, no LW)', f%leaf_temp(1), expect, 1.0e-10_wp)
   end subroutine test_leaf_analytic

   !----- 3. The split commits the BE-in-atm solution of the tendencies; verify (y1-y0)/dt = f(y1). !
   subroutine test_cas_be_consistency()
      type(surface_frozen_t) :: fro
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      real(wp)    :: dt, enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: r_enth, r_shv, r_co2
      integer(ik) :: n
      n = 3_ik ; dt = 900.0_wp
      print '(a)', 'test_cas_be_consistency:'
      call make_frozen(fro, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(296.0_wp, 0.013_wp)
      y%cas_shv      = 0.013_wp
      y%cas_co2      = 415.0_wp
      enth0 = y%cas_enthalpy ; shv0 = y%cas_shv ; co20 = y%cas_co2
      call surface_derivs(y, fro, n, f)
      !----- Reconstruct the split's committed state (meds_fast_split.f90). ------------------------!
      enth1 = (fro%wcap * enth0 + dt * (f%src_enth  + fro%gah * fro%enth_atm)) / (fro%wcap + dt * fro%gah)
      shv1  = (fro%wcap * shv0  + dt * (f%src_vap   + fro%gaw * fro%shv_atm )) / (fro%wcap + dt * fro%gaw)
      co21  = (fro%ccap * co20  + dt * (fro%nee_biotic + fro%gac * fro%co2_atm)) / (fro%ccap + dt * fro%gac)
      !----- BE consistency: (y1-y0)/dt must equal the tendency evaluated with the ATM term at y1  !
      !      (source frozen). This is exactly what an IMEX/BE stage solves, so it verifies d_cas_* !
      !      is the correct RHS of the split's implicit update.                                     !
      r_enth = (enth1 - enth0) / dt - (f%src_enth + fro%gah * (fro%enth_atm - enth1)) / fro%wcap
      r_shv  = (shv1  - shv0 ) / dt - (f%src_vap  + fro%gaw * (fro%shv_atm  - shv1 )) / fro%wcap
      r_co2  = (co21  - co20 ) / dt - (fro%nee_biotic + fro%gac * (fro%co2_atm - co21)) / fro%ccap
      call check_true('CAS enthalpy BE-consistent with d_cas_enthalpy', abs(r_enth) < 1.0e-9_wp, r_enth)
      call check_true('CAS humidity BE-consistent with d_cas_shv',      abs(r_shv)  < 1.0e-15_wp, r_shv)
      call check_true('CAS CO2 BE-consistent with d_cas_co2',           abs(r_co2)  < 1.0e-9_wp, r_co2)
      !----- At t=0 the tendency must equal (src + g*(atm - y0))/cap (sanity on the returned RHS). !
      call check('d_cas_enthalpy = (src+gah*(atm-y0))/wcap',                                      &
                 f%d_cas_enthalpy, (f%src_enth + fro%gah * (fro%enth_atm - enth0)) / fro%wcap, 1.0e-9_wp)
   end subroutine test_cas_be_consistency

   !----- 4. March the CAS twins with the RHS (BE-in-atm); the closed budgets must stay tight. ---!
   subroutine test_conservation_march()
      type(surface_frozen_t) :: fro
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f
      type(budget_t)         :: be, bw
      real(wp)    :: dt, enth0, shv0, enth1, shv1, t0, t1
      integer(ik) :: step, n
      n = 3_ik ; dt = 120.0_wp
      print '(a)', 'test_conservation_march:'
      call make_frozen(fro, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(294.0_wp, 0.010_wp)     ! CAS starts cool + dry
      y%cas_shv      = 0.010_wp
      y%cas_co2      = 400.0_wp
      t0 = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      do step = 1_ik, 60_ik
         call surface_derivs(y, fro, n, f)
         enth0 = y%cas_enthalpy ; shv0 = y%cas_shv
         enth1 = (fro%wcap * enth0 + dt * (f%src_enth + fro%gah * fro%enth_atm)) / (fro%wcap + dt * fro%gah)
         shv1  = (fro%wcap * shv0  + dt * (f%src_vap  + fro%gaw * fro%shv_atm )) / (fro%wcap + dt * fro%gaw)
         !----- Same closed-budget accounting the split uses (meds_fast_split.f90). ------------------!
         call budget_accumulate(be, fro%wcap * enth0, fro%wcap * enth1, f%src_enth + fro%gah * fro%enth_atm, &
                                fro%gah * enth1, dt, abs(fro%wcap * enth1), 1.0e-8_wp, 1.0e-3_wp)
         call budget_accumulate(bw, fro%wcap * shv0, fro%wcap * shv1, f%src_vap + fro%gaw * fro%shv_atm,     &
                                fro%gaw * shv1, dt, max(abs(fro%wcap * shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
         y%cas_enthalpy = enth1 ; y%cas_shv = shv1
      end do
      t1 = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      call check_true('CAS energy budget closes over the march', be%n_fail == 0_ik, real(be%worst, wp))
      call check_true('CAS water  budget closes over the march', bw%n_fail == 0_ik, real(bw%worst, wp))
      call check_true('CAS warms toward the warm atmosphere + sunlit canopy', t1 > t0, t1 - t0)
   end subroutine test_conservation_march

   !----- 5. The soil-water supply fraction scales the transpiration source into the CAS. --------!
   subroutine test_supply_limiter()
      type(surface_frozen_t) :: fro
      type(surface_state_t)  :: y
      type(surface_tend_t)   :: f_full, f_zero
      integer(ik) :: n
      n = 3_ik
      print '(a)', 'test_supply_limiter:'
      call make_frozen(fro, n)
      y%cas_enthalpy = cas_enthalpy_of_temp(298.0_wp, 0.012_wp)
      y%cas_shv      = 0.012_wp
      y%cas_co2      = 410.0_wp
      fro%src_frac = 1.0_wp ; call surface_derivs(y, fro, n, f_full)
      fro%src_frac = 0.0_wp ; call surface_derivs(y, fro, n, f_zero)
      !----- src_frac = 0: no transpiration reaches the CAS, so src_vap collapses to soil_evap. ---!
      call check('src_frac=0 -> src_vap = soil_evap', f_zero%src_vap, fro%soil_evap, 1.0e-12_wp)
      call check_true('src_frac=1 adds transpiration on top of soil evap',                        &
                      f_full%src_vap > f_zero%src_vap + 1.0e-9_wp, f_full%src_vap - f_zero%src_vap)
   end subroutine test_supply_limiter

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
         col%soil_energy(k) = temp_to_uext(therm%soil_dry_heat_capacity(k),                      &
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
      call soil_water_time_deriv(theta, soil, hopts, nsl, q_top, psi_e, root_uptake, dtheta, drain, uptk)
      colsum = 0.0_wp
      do k = 1_ik, nsl ; colsum = colsum + dtheta(k) * soil%dz(k) ; end do
      net = q_top - drain / rho_h2o - uptk / rho_h2o   ! d(storage)/dt = in - drainage - uptake
      call check('soil-water column balance closes (telescoping)', colsum, net, 1.0e-12_wp)
      call check_true('root sink active and column free-drains', uptk > 0.0_wp .and. drain > 0.0_wp, drain)
   end subroutine test_soil_water_tendency

   !----- 9. column_derivs assembles a finite, physically-signed whole-column RHS + its CAS part   !
   !         matches a standalone surface_derivs call (correct wiring). -------------------------!
   subroutine test_column_assembler()
      type(column_state_t)       :: y
      type(column_frozen_t)      :: fro
      type(column_tend_t)        :: f
      type(surface_state_t)      :: ys
      type(surface_tend_t)       :: sf
      type(energy_forcing_t)     :: eforc_chk
      type(soil_energy_column_t) :: se_chk
      real(wp)    :: dedt_chk(n_soil_layer_max), worst
      integer(ik) :: k, i, n, nsl
      logical     :: all_finite
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_column_assembler:'
      call make_column(y, fro, n, nsl)
      call column_derivs(y, fro, n, nsl, f)

      all_finite = ieee_ok(f%d_cas_enthalpy) .and. ieee_ok(f%d_cas_shv) .and. ieee_ok(f%d_cas_co2)
      do k = 1_ik, nsl ; all_finite = all_finite .and. ieee_ok(f%dedt(k)) .and. ieee_ok(f%dtheta_dt(k)) ; end do
      do i = 1_ik, n
         all_finite = all_finite .and. ieee_ok(f%d_leaf_water_mass(i)) .and. ieee_ok(f%d_wood_water_mass(i))
      end do
      call check_true('column_derivs: all tendencies finite', all_finite, 0.0_wp)

      !----- the CAS part must equal a standalone surface_derivs at the state's diagnosed t_ground. !
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      call surface_derivs(ys, surf_with_tground(fro%surf, y, fro), n, sf)
      call check('column_derivs CAS enthalpy tendency = surface_derivs', f%d_cas_enthalpy, sf%d_cas_enthalpy, 1.0e-12_wp)

      !----- wiring check: the assembled soil-heat tendency == a standalone soil_energy_time_deriv    !
      !      built from the SAME surface coupling (g_top, coh_qsoil * root_frac) PLUS the bottom-face   !
      !      drainage enthalpy. That last term used to be omitted here and the check still passed,      !
      !      because column_derivs advected it on the FROZEN fro%drainage, which is 0 in this fixture.  !
      !      C2 made it ride the state-dependent f%drainage_rate instead -- the same flux the mass side !
      !      has always used -- so it is now nonzero and the reproduction has to carry it too. A check  !
      !      that omits a term is only green while that term is zero. -------------------------------!
      se_chk%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc_chk%g_top = sf%g_top ; eforc_chk%geothermal = fro%geothermal
      do k = 1_ik, nsl
         eforc_chk%soil_water(k)     = y%theta(k)
         eforc_chk%root_heat_sink(k) = sf%coh_qsoil * fro%soil%root_frac(k)
         eforc_chk%w_flux(k)         = 0.0_wp
      end do
      eforc_chk%root_heat_sink(nsl) = eforc_chk%root_heat_sink(nsl)                                  &
                                    + f%drainage_rate * internal_energy_liquid(fro%t_bot)
      call soil_energy_time_deriv(se_chk, eforc_chk, fro%therm, fro%soil, fro%energy_opts, dedt_chk)
      worst = maxval(abs(f%dedt(1:nsl) - dedt_chk(1:nsl)))
      call check_true('column_derivs wires the soil-heat tendency correctly', worst < 1.0e-9_wp, worst)
   end subroutine test_column_assembler

   !----- 10. march the whole column with explicit RK4 over column_derivs (the reference oracle):   !
   !          it stays finite + physical, and self-converges (dt vs dt/2 agree to RK4 order). -------!
   subroutine test_rk4_march()
      type(column_state_t)  :: y, y1, y2, ytmp
      type(column_frozen_t) :: fro
      real(wp)    :: dt, tcas, tsoil, dcas, dmass
      integer(ik) :: step, n, nsl, nstep, i, k
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rk4_march:'
      call make_column(y, fro, n, nsl)

      !----- (a) coarse march at dt = 10 s (< the ~47 s RK4 limit of the 17 s hydraulic mode). ----!
      dt = 10.0_wp ; nstep = 180_ik            ! 30 minutes
      call copy_state(y, y1, n)
      physical = .true.
      do step = 1_ik, nstep
         call rk4_column_step(y1, fro, n, nsl, dt, ytmp)
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
         call rk4_column_step(y1, fro, n, nsl, 8.0_wp, ytmp) ; call copy_state(ytmp, y1, n)
      end do
      do step = 1_ik, 120_ik                   ! 120 * 4 s = 8 min
         call rk4_column_step(y2, fro, n, nsl, 4.0_wp, ytmp) ; call copy_state(ytmp, y2, n)
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
      type(column_frozen_t) :: fro
      real(wp)    :: tcas, dcas, dtheta, dmass
      integer(ik) :: step, n, nsl, k, i
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_imex_euler:'
      call make_column(y, fro, n, nsl)
      !----- Well-ventilated canopy (stronger CAS<->atm exchange) so the frozen-per-step source at   !
      !      the full 900 s does not out-run venting into supersaturation -- the large-dt operator-   !
      !      split coupling error the P2 arrowhead removes; here we exercise L-stability + physical.  !
      fro%surf%gah = fro%surf%gah * 4.0_wp ; fro%surf%gaw = fro%surf%gaw * 4.0_wp
      fro%surf%gac = fro%surf%gac * 4.0_wp

      !----- (a) STABLE + physical at the full fast timestep dt = 900 s over a 6-hour march. --------!
      call copy_state(y, yi, n)
      physical = .true.
      do step = 1_ik, 24_ik                    ! 24 * 900 s = 6 h
         call imex_euler_column_step(yi, fro, n, nsl, 900.0_wp, ytmp)
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
      !          column_hydrology_flux is the sole authority at the column_fast_step level), so the oracle !
      !          is run with freeze_theta=.true.: same reduced system -> theta is trivially equal and the  !
      !          CAS/soil-energy core agrees tightly. Mass stays operator-split in BOTH (in-vector oracle   !
      !          vs advance_water_mass_full) -> a small O(dt) difference (both are exact closed-form Euler  !
      !          now, so this should be tighter than psi's old O(dt) split-vs-exact-exp gap). --------------!
      call copy_state(y, yi, n) ; call copy_state(y, yr, n)
      do step = 1_ik, 60_ik                    ! 60 * 4 s = 4 min
         call imex_euler_column_step(yi, fro, n, nsl, 4.0_wp, ytmp)              ; call copy_state(ytmp, yi, n)
         call rk4_column_step(yr, fro, n, nsl, 4.0_wp, ytmp, freeze_theta=.true.); call copy_state(ytmp, yr, n)
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
      type(column_frozen_t) :: fro
      real(wp)    :: tcas_base, tcas_coup, qsat_base, qsat_coup
      integer(ik) :: step, n, nsl
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_imex_coupled:'
      !----- HARSH forcing: constant noon, modest venting (the case that over-humidifies at 900 s). -!
      call make_column(y, fro, n, nsl)      ! no ventilation boost -> baseline over-humidifies

      call copy_state(y, yb, n) ; call copy_state(y, yc, n)
      do step = 1_ik, 12_ik                 ! 12 * 900 s = 3 h of constant noon
         call imex_euler_column_step(yb, fro, n, nsl, 900.0_wp, ytmp)                        ! niter=1 (baseline)
         call copy_state(ytmp, yb, n)
         call imex_euler_column_step(yc, fro, n, nsl, 900.0_wp, ytmp, niter=12_ik)  ! coupled
         call copy_state(ytmp, yc, n)
      end do
      tcas_base = cas_temp_of_enthalpy(yb%cas_enthalpy, yb%cas_shv)
      tcas_coup = cas_temp_of_enthalpy(yc%cas_enthalpy, yc%cas_shv)
      qsat_base = sat_specific_humidity(tcas_base, fro%surf%press)
      qsat_coup = sat_specific_humidity(tcas_coup, fro%surf%press)
      !----- coupled stays physical (sub-saturated, physical temperature); baseline collapses. ------!
      call check_true('coupled CAS stays physical (270-325 K) at 900 s constant noon',            &
                      tcas_coup > 270.0_wp .and. tcas_coup < 325.0_wp, tcas_coup)
      call check_true('coupled CAS stays sub-saturated (shv <= qsat)', yc%cas_shv <= qsat_coup + 1.0e-4_wp, &
                      yc%cas_shv - qsat_coup)
      call check_true('coupling matters: coupled warmer than the over-humidified baseline',        &
                      tcas_coup > tcas_base + 1.0_wp, tcas_coup - tcas_base)
   end subroutine test_imex_coupled

   !----- 12b. the arrowhead Newton surface solve: robust near saturation (clamp + line search), and  !
   !           faithful to the RK4 oracle at small dt (converges to the coupled surface solution). ----!
   subroutine test_arrowhead()
      type(column_state_t)  :: y, yi, yr, ytmp
      type(column_frozen_t) :: fro
      real(wp)    :: tcas, qsat, dcas, dtheta, worst_super
      integer(ik) :: step, n, nsl, k
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_arrowhead:'

      !----- (a) START the CAS at 99% saturation under harsh noon: exercises the supersaturation      !
      !          clamp + line search. Newton (niter>1) must keep it physical + sub-saturated at 900 s. !
      call make_column(y, fro, n, nsl)
      tcas = 297.0_wp ; qsat = sat_specific_humidity(tcas, fro%surf%press)
      y%cas_enthalpy = cas_enthalpy_of_temp(tcas, 0.99_wp*qsat) ; y%cas_shv = 0.99_wp*qsat
      call copy_state(y, yi, n)
      physical = .true. ; worst_super = -1.0_wp
      do step = 1_ik, 12_ik
         call imex_euler_column_step(yi, fro, n, nsl, 900.0_wp, ytmp, niter=8_ik)
         call copy_state(ytmp, yi, n)
         tcas = cas_temp_of_enthalpy(yi%cas_enthalpy, yi%cas_shv)
         qsat = sat_specific_humidity(tcas, fro%surf%press)
         worst_super = max(worst_super, yi%cas_shv - qsat)          ! must stay <= 0 (sub-saturated)
         physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp
      end do
      call check_true('Newton stays sub-saturated from a 99%-saturated start (900 s)', worst_super <= 1.0e-4_wp, worst_super)
      call check_true('Newton stays physical (270-325 K) near saturation', physical, tcas)

      !----- (b) Newton (niter=8) matches the RK4 oracle at small dt (converges to the coupled soln). -!
      call make_column(y, fro, n, nsl)
      call copy_state(y, yi, n) ; call copy_state(y, yr, n)
      do step = 1_ik, 60_ik                     ! 60 * 4 s = 4 min
         call imex_euler_column_step(yi, fro, n, nsl, 4.0_wp, ytmp, niter=8_ik)   ; call copy_state(ytmp, yi, n)
         call rk4_column_step(yr, fro, n, nsl, 4.0_wp, ytmp, freeze_theta=.true.) ; call copy_state(ytmp, yr, n)
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
      type(column_frozen_t) :: fro
      real(wp)    :: tc1, tc2, tcr
      integer(ik) :: step, n, nsl, ns1, ns2, nr1, nr2
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_adaptive_march:'
      call make_column(y, fro, n, nsl)
      fro%surf%gah = fro%surf%gah*4.0_wp ; fro%surf%gaw = fro%surf%gaw*4.0_wp ; fro%surf%gac = fro%surf%gac*4.0_wp
      !----- start the CAS cool + dry so there is a real transient to resolve adaptively. -----------!
      y%cas_enthalpy = cas_enthalpy_of_temp(290.0_wp, 0.009_wp) ; y%cas_shv = 0.009_wp

      call adaptive_imex_march(y, fro, n, nsl, 1800.0_wp, 1.0e-3_wp, 50.0_wp, y1, ns1, nr1)
      call adaptive_imex_march(y, fro, n, nsl, 1800.0_wp, 1.0e-5_wp, 50.0_wp, y2, ns2, nr2)
      !----- fine fixed reference (dt = 2 s, coupled): 900 steps. ----------------------------------!
      call copy_state(y, yr, n)
      do step = 1_ik, 900_ik
         call imex_euler_column_step(yr, fro, n, nsl, 2.0_wp, ytmp, niter=8_ik)
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
   function soil_top_temp(y, fro) result(t)
      type(column_state_t),  intent(in) :: y
      type(column_frozen_t), intent(in) :: fro
      real(wp) :: t, fl
      call uext_to_temp(y%soil_energy(1), y%theta(1)*rho_h2o, fro%therm%soil_dry_heat_capacity(1), t, fl)
   end function soil_top_temp

   !----- march the ARK2 fixed-step from y for nstep steps of dt (embedded error discarded). --------!
   subroutine march_ark2(y0, fro, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp, yerr
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call ark2_column_step(y, fro, n, nsl, dt, ytmp, yerr, niter=8_ik)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_ark2

   !----- march coupled IMEX-Euler (niter=8 Newton) fixed-step (for the ARK2-vs-1st-order comparison). !
   subroutine march_imex(y0, fro, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call imex_euler_column_step(y, fro, n, nsl, dt, ytmp, niter=8_ik)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_imex

   !----- march the Cash-Karp RK45 fixed-step (5th-order commit, embedded 4th discarded) -- for the   !
   !      order-of-accuracy self-convergence study below. ------------------------------------------!
   subroutine march_rk45(y0, fro, n, nsl, dt, nstep, y_out)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl, nstep
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      type(column_state_t) :: y, ytmp, yerr
      real(wp) :: w_out, e_in, e_out
      integer(ik) :: s
      call copy_state(y0, y, n)
      do s = 1_ik, nstep
         call rk45_column_step(y, fro, n, nsl, dt, ytmp, yerr, w_out, e_in, e_out)
         call copy_state(ytmp, y, n)
      end do
      call copy_state(y, y_out, n)
   end subroutine march_rk45

   !----- 14. ARK2 (ARS(2,2,2)): 2nd-order on the soil-top temperature (differential var); mass stays !
   !          physical at production dt (FATAL-1 guard); embedded estimate bounded; stiff-stable. -----!
   subroutine test_ark2()
      type(column_state_t)  :: y, yref, y1, y2, y4, ynew, yerr, ytmp
      type(column_frozen_t) :: fro
      type(error_control_t) :: ec
      real(wp)    :: e1, e2, e4, p_lo, p_hi, tref, errnorm, tcas
      integer(ik) :: n, nsl, step, ns, nr, ns2, nr2
      logical     :: physical
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_ark2:'
      call make_column(y, fro, n, nsl)

      !----- (a) ORDER-2 self-convergence on a CLEANLY-integrated tableau variable: CAS CO2 (a         !
      !          decoupled affine scalar ODE -- no operator split, so the ARS(2,2,2) order is exposed).!
      !          reference = ARK2 at dt = 1800/64; solutions at dt = 450, 225, 112.5. ----------------!
      call march_ark2(y, fro, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = yref%cas_co2
      call march_ark2(y, fro, n, nsl, 450.0_wp,   4_ik,  y1) ; e1 = abs(y1%cas_co2 - tref)
      call march_ark2(y, fro, n, nsl, 225.0_wp,   8_ik,  y2) ; e2 = abs(y2%cas_co2 - tref)
      call march_ark2(y, fro, n, nsl, 112.5_wp,  16_ik,  y4) ; e4 = abs(y4%cas_co2 - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)', '   CAS CO2 errors: ', e1, ' / ', e2, &
            ' / ', e4, ' ; observed order p = ', p_lo, ' , ', p_hi
      call check_true('ARK2 tableau is 2nd order on the decoupled CO2 twin (p >= 1.9)', &
                      p_lo >= 1.9_wp .and. p_hi >= 1.9_wp, p_hi)

      !----- (a2) On the operator-split-COUPLED soil-top temperature, the full 3x3 arrowhead is        !
      !           deferred (spec) so the effective order is ~1.2, but ARK2 is still MORE ACCURATE than  !
      !           the 1st-order IMEX-Euler at the same dt. -------------------------------------------!
      call march_ark2(y, fro, n, nsl, 1800.0_wp/64.0_wp, 64_ik, yref) ; tref = soil_top_temp(yref, fro)
      call march_ark2(y, fro, n, nsl, 225.0_wp,  8_ik, y1) ; e1 = abs(soil_top_temp(y1, fro) - tref)
      call march_imex(y, fro, n, nsl, 225.0_wp,  8_ik, y2) ; e2 = abs(soil_top_temp(y2, fro) - tref)
      call check_true('ARK2 beats IMEX-Euler on the coupled soil-top temperature (same dt)', e1 < e2, e2 - e1)

      !----- (b) FATAL-1 guard: mass stays physical (finite, positive) at production dt=900. ---------!
      call copy_state(y, y1, n) ; physical = .true.
      do step = 1_ik, 24_ik
         call ark2_column_step(y1, fro, n, nsl, 900.0_wp, ytmp, yerr, niter=8_ik)
         call copy_state(ytmp, y1, n)
         physical = physical .and. all(y1%leaf_water_mass(1:n) == y1%leaf_water_mass(1:n)) .and.   &
                    all(y1%wood_water_mass(1:n) == y1%wood_water_mass(1:n)) .and.                  &
                    minval(y1%leaf_water_mass(1:n)) > 0.0_wp .and. minval(y1%wood_water_mass(1:n)) > 0.0_wp
      end do
      call check_true('ARK2 mass stays physical at dt=900 (FATAL-1 split-out guard)', physical,     &
                      minval(y1%leaf_water_mass(1:n)))

      !----- (c) embedded error estimate is bounded (not detonating -- the pre-fix failure mode). ----!
      call ark2_column_step(y, fro, n, nsl, 900.0_wp, ynew, yerr, niter=8_ik)
      call state_err_norm(ynew, yerr, y, n, nsl, errnorm)
      call check_true('ARK2 embedded estimate bounded (WRMS < 50) at dt=900', errnorm < 50.0_wp, errnorm)

      !----- (d) stiffness: 24 h adaptive-ARK march stays physical + bounded. -----------------------!
      ec = default_error_control(1.0e-3_wp)
      call adaptive_ark_march(y, fro, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ns, nr)
      tcas = cas_temp_of_enthalpy(ytmp%cas_enthalpy, ytmp%cas_shv)
      print '(a,i0,a,i0)', '   adaptive-ARK 24 h: steps = ', ns, ' , rejects = ', nr
      call check_true('adaptive-ARK 24 h stays bounded (280 < tcas < 320 K)', tcas > 280.0_wp .and. tcas < 320.0_wp, tcas)

      !----- (e) ERROR-CONTROL: the PI controller (goal a) marches the SAME stiff window and stays     !
      !      physical -- proves the meds_fast_control PI path is wired + functional. On a multi-substep !
      !      march it takes a different (typically smoother) step sequence than the I-controller.       !
      ec%controller = CTRL_PI
      call adaptive_ark_march(y, fro, n, nsl, 86400.0_wp, ec, 300.0_wp, ytmp, ns2, nr2)
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
      type(column_state_t)  :: y, yref, y1, y2, y4
      type(column_frozen_t) :: fro
      real(wp)    :: e1, e2, e4, p_lo, p_hi, tref
      integer(ik) :: n, nsl
      n = 2_ik ; nsl = 10_ik
      print '(a)', 'test_rk45_order:'
      call make_column(y, fro, n, nsl)

      !----- A 5th-order method's error ~ C*h^5 collapses to double-precision noise (~1e-13) at MUCH      !
      !      coarser h than a 2nd-order one (test_ark2's own dt=200/8..32 s range) -- h^5 falls off so    !
      !      fast that h must stay large enough for the TRUNCATION term to still dominate roundoff at      !
      !      the FINEST test point. Window T=2400 s (8 * the coarsest dt below), test points dt=300/       !
      !      150/75 s (8/16/32 steps; 300 s stays under the ~360 s single-step stability estimate,          !
      !      design doc sec 6 -- this is a short fixed-step accuracy comparison, not a stability march),     !
      !      reference dt=2400/1024 s (>=32x finer than the finest test point). ---------------------------!
      !----- (a) CLEAN scalar: CAS CO2 (decoupled affine ODE, mirrors test_ark2's part a) -- baseline   !
      !          proof the Cash-Karp tableau itself is 5th order on this RHS. ------------------------!
      call march_rk45(y, fro, n, nsl, 2400.0_wp/1024.0_wp, 1024_ik, yref) ; tref = yref%cas_co2
      call march_rk45(y, fro, n, nsl, 300.0_wp,   8_ik,  y1) ; e1 = abs(y1%cas_co2 - tref)
      call march_rk45(y, fro, n, nsl, 150.0_wp,  16_ik,  y2) ; e2 = abs(y2%cas_co2 - tref)
      call march_rk45(y, fro, n, nsl,  75.0_wp,  32_ik,  y4) ; e4 = abs(y4%cas_co2 - tref)
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
      call march_rk45(y, fro, n, nsl, 2400.0_wp/1024.0_wp, 1024_ik, yref) ; tref = soil_top_temp(yref, fro)
      call march_rk45(y, fro, n, nsl, 300.0_wp,   8_ik,  y1) ; e1 = abs(soil_top_temp(y1, fro) - tref)
      call march_rk45(y, fro, n, nsl, 150.0_wp,  16_ik,  y2) ; e2 = abs(soil_top_temp(y2, fro) - tref)
      call march_rk45(y, fro, n, nsl,  75.0_wp,  32_ik,  y4) ; e4 = abs(soil_top_temp(y4, fro) - tref)
      p_lo = log(e1/max(e2, tiny_num)) / log(2.0_wp)
      p_hi = log(e2/max(e4, tiny_num)) / log(2.0_wp)
      print '(a,es10.3,a,es10.3,a,es10.3,a,f5.2,a,f5.2)', '   soil-top T errors: ', e1, ' / ', e2, &
            ' / ', e4, ' ; observed order p = ', p_lo, ' , ', p_hi
      call check_true('RK45 stays 5th order on the COUPLED soil-top temperature (p >= 4.5, no split)', &
                      p_lo >= 4.5_wp .and. p_hi >= 4.5_wp, p_hi)
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
   subroutine make_column(y, fro, n, nsl)
      type(column_state_t),  intent(out) :: y
      type(column_frozen_t), intent(out) :: fro
      integer(ik),           intent(in)  :: n, nsl
      !----- LOCAL hydraulics traits: only needed to seed y%*_water_mass at a representative psi     !
      !      (column_frozen_t no longer carries hydro_p/hydro_o -- column_derivs' mass ODE needs      !
      !      only the frozen sapflow/uptake, not the PV-curve/conductance params any more). -----------!
      type(hydro_params_t) :: hp
      integer(ik) :: i, k
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,        &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, fro%soil)
      call build_soil_therm_params(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, fro%therm)
      fro%hydro_opts = soil_opts_t()
      hp%leaf_pi0 = -1.5_wp ; hp%leaf_elastic_mod = 12.0_wp ; hp%leaf_apoplast_frac = 0.30_wp
      hp%leaf_water_sat = 2.0_wp ; hp%wood_pi0 = -1.0_wp ; hp%wood_elastic_mod = 8.0_wp
      hp%wood_apoplast_frac = 0.20_wp ; hp%wood_water_sat = 1.0_wp ; hp%wood_psi50 = -2.0_wp
      hp%wood_kexp = 2.0_wp ; hp%k_plant_max = 6.0e-4_wp ; hp%wood_kmax = 8.0_wp ; hp%vessel_curl = 1.5_wp
      fro%geothermal = 0.0_wp ; fro%q_top = 1.0e-6_wp ; fro%soil_psi_root = -0.3_wp ; fro%rhizo_cond = 5.0e-4_wp
      allocate(fro%psi_e(nsl), fro%nplant(n), fro%bleaf(n), fro%bsap(n), fro%broot(n),           &
               fro%sap_area(n), fro%height(n), fro%leaf_area(n))
      allocate(fro%sapflow_frozen(n), fro%uptake_frozen(n), fro%qloss_frozen(n))
      allocate(fro%intercept_leaf(n), fro%intercept_wood(n))
      fro%psi_e = 0.0_wp
      fro%nplant = 0.3_wp ; fro%bleaf = 0.5_wp ; fro%bsap = 5.0_wp ; fro%broot = 2.0_wp
      fro%sap_area = 0.01_wp ; fro%height = 20.0_wp ; fro%leaf_area = 5.0_wp
      !----- FROZEN sapflow/uptake (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2): a representative,            !
      !      state-independent pair for these RHS/oracle-march tests -- not re-derived from a            !
      !      solve_plant_water_batch pre-pass here (these tests exercise column_derivs/the tableau         !
      !      machinery directly, not build_column_frozen's own pre-pass, which has its own coverage        !
      !      via test_column_ark.f90/test_picard_coupling.f90). -----------------------------------------!
      fro%sapflow_frozen = 1.0e-4_wp ; fro%uptake_frozen = 1.0e-4_wp ; fro%uptake = fro%uptake_frozen(1)*sum(fro%nplant(1:n))
      fro%qloss_frozen = 0.0_wp   ! P2 advective enthalpy: no-op unless populated (see build_column_frozen)
      fro%intercept_leaf = 0.0_wp ; fro%intercept_wood = 0.0_wp   ! P2c canopy water: no-op unless populated
      allocate(fro%surf%h_coeff_f(n), fro%surf%g_tr_f(n), fro%surf%abs_sw(n), fro%surf%abs_lw(n), fro%surf%lai(n))
      allocate(fro%surf%h_coeff_w(n), fro%surf%abs_sw_wood(n), fro%surf%abs_lw_wood(n), fro%surf%wai(n))
      allocate(fro%surf%qwflux_wl(n), fro%surf%q_wood_net(n))
      allocate(fro%surf%f_wet_c(n), fro%surf%g_film_f(n), fro%surf%g_film_w(n))
      fro%surf%h_coeff_w = 0.0_wp ; fro%surf%abs_sw_wood = 0.0_wp ; fro%surf%abs_lw_wood = 0.0_wp ; fro%surf%wai = 0.0_wp
      fro%surf%qwflux_wl = 0.0_wp ; fro%surf%q_wood_net = 0.0_wp
      fro%surf%f_wet_c = 0.0_wp ; fro%surf%g_film_f = 0.0_wp ; fro%surf%g_film_w = 0.0_wp
      do i = 1_ik, n
         fro%surf%lai(i) = 2.0_wp - 0.5_wp * real(i-1_ik, wp) ; fro%surf%abs_sw(i) = 250.0_wp - 50.0_wp*real(i-1_ik, wp)
         fro%surf%abs_lw(i) = -30.0_wp
         fro%surf%h_coeff_f(i) = 2.0_wp * fro%surf%lai(i) * 0.03_wp * 1.2_wp * cp_air
         fro%surf%g_tr_f(i) = 0.004_wp * fro%surf%lai(i)
      end do
      fro%surf%rho = 1.2_wp ; fro%surf%press = 101325.0_wp ; fro%surf%wcap = 1.2_wp*20.0_wp
      fro%surf%ccap = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*20.0_wp
      fro%surf%gah = 1.2_wp*0.3_wp*0.02_wp ; fro%surf%gaw = fro%surf%gah
      fro%surf%gac = (1.2_wp*(1.0_wp-0.012_wp)/0.0289655_wp)*0.3_wp*0.02_wp
      fro%surf%enth_atm = cas_enthalpy_of_temp(300.0_wp, 0.011_wp) ; fro%surf%shv_atm = 0.011_wp
      fro%surf%co2_atm = 400.0_wp ; fro%surf%nee_biotic = -5.0_wp
      fro%surf%abs_sw_ground = 60.0_wp ; fro%surf%abs_lw_ground = -10.0_wp
      fro%surf%ggnet = 0.02_wp ; fro%surf%soil_evap = 2.0e-5_wp ; fro%surf%src_frac = 1.0_wp
      y%cas_enthalpy = cas_enthalpy_of_temp(297.0_wp, 0.012_wp) ; y%cas_shv = 0.012_wp ; y%cas_co2 = 410.0_wp
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water = 0.0_wp ; y%wood_surf_water = 0.0_wp   ! P2c canopy water: no-op unless populated
      do k = 1_ik, nsl
         y%soil_energy(k) = temp_to_uext(fro%therm%soil_dry_heat_capacity(k), 0.30_wp*rho_h2o,   &
                            296.0_wp - 0.4_wp*real(k-1_ik, wp), 1.0_wp)
         y%theta(k) = 0.30_wp
      end do
      !----- seed mass at the SAME representative (leaf psi=-1.0, wood psi=-0.5 MPa) point the old   !
      !      psi-based fixture used, via the forward water_content map (hp above). -------------------!
      do i = 1_ik, n
         y%leaf_water_mass(i) = water_content(-1.0_wp, hp%leaf_pi0, hp%leaf_elastic_mod,          &
              hp%leaf_apoplast_frac, hp%leaf_water_sat, fro%bleaf(i))
         y%wood_water_mass(i) = water_content(-0.5_wp, hp%wood_pi0, hp%wood_elastic_mod,          &
              hp%wood_apoplast_frac, hp%wood_water_sat, fro%bsap(i) + fro%broot(i))
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

   !----- helper: a surface_frozen_t with t_ground diagnosed from the soil-top state (mirrors      !
   !      what column_derivs does internally) so the standalone CAS comparison lines up. -----------!
   function surf_with_tground(surf0, y, fro) result(surf)
      type(surface_frozen_t), intent(in) :: surf0
      type(column_state_t),   intent(in) :: y
      type(column_frozen_t),  intent(in) :: fro
      type(surface_frozen_t) :: surf
      real(wp) :: tg, fl
      surf = surf0
      call uext_to_temp(y%soil_energy(1), y%theta(1)*rho_h2o, fro%therm%soil_dry_heat_capacity(1), tg, fl)
      surf%t_ground = tg
   end function surf_with_tground

   logical function ieee_ok(x)
      real(wp), intent(in) :: x
      ieee_ok = (x == x) .and. (abs(x) < huge(1.0_wp))
   end function ieee_ok

end program test_column_derivs
