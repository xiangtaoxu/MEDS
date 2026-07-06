!==========================================================================================!
! test_column_energy -- unit tests for the soil THERMAL column + shared enthalpy inverter.    !
!   1. INVERTER round-trip (liquid + ice) and CONTINUITY at u_freeze/u_melt (temp = t_3ple).   !
!   2. CLAUSIUS slope d(e_sat)/dT vs finite difference.                                         !
!   3. SOIL energy CONSERVATION to ~round-off (a G_top warming; and with advective w_flux).     !
!   4. STEADY STATE: a sealed isothermal column does not drift.                                 !
!==========================================================================================!
program test_column_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, t_3ple, cp_ice, latent_heat_fusion
   use meds_biophysics_types, only : soil_energy_column_t, energy_forcing_t, soil_thermal_params_t, &
                                     soil_params_t, energy_opts_t, energy_flux_t, SOIL_RETENTION_VG
   use meds_soil_parameters,  only : build_soil_params
   use meds_soil_thermal,     only : build_soil_thermal
   use meds_thermo,           only : temp_to_uext, uext_to_temp, sat_vapor_pressure,           &
                                     d_sat_vapor_pressure_dt
   use meds_column_energy,    only : soil_energy_flux
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_inverter()
   call test_clausius()
   call test_soil_conserve()
   call test_soil_steady()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_energy: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_energy: ', nfail, ' FAILED'
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

   subroutine test_inverter()
      real(wp), parameter :: dh = 1000.0_wp, wm = 0.3_wp
      real(wp) :: uext, t, fl, u_freeze, u_melt
      print '(a)', 'test_inverter:'
      uext = temp_to_uext(dh, wm, 290.0_wp, 1.0_wp)
      call uext_to_temp(uext, wm, dh, t, fl)
      call check('liquid T round-trip', t, 290.0_wp, 1.0e-9_wp)
      call check_true('liquid fliq = 1', abs(fl - 1.0_wp) < 1.0e-12_wp, fl)
      uext = temp_to_uext(dh, wm, 260.0_wp, 0.0_wp)
      call uext_to_temp(uext, wm, dh, t, fl)
      call check('ice T round-trip', t, 260.0_wp, 1.0e-9_wp)
      call check_true('ice fliq = 0', abs(fl) < 1.0e-12_wp, fl)
      u_freeze = (dh + wm * cp_ice) * t_3ple
      u_melt   = u_freeze + wm * latent_heat_fusion
      call uext_to_temp(u_freeze, wm, dh, t, fl) ; call check('continuity at u_freeze', t, t_3ple, 1.0e-9_wp)
      call uext_to_temp(u_melt,   wm, dh, t, fl) ; call check('continuity at u_melt',   t, t_3ple, 1.0e-9_wp)
      call uext_to_temp(0.5_wp * (u_freeze + u_melt), wm, dh, t, fl)
      call check('plateau temp = t_3ple', t, t_3ple, 1.0e-9_wp)
      call check('plateau fliq = 0.5', fl, 0.5_wp, 1.0e-9_wp)
   end subroutine test_inverter

   subroutine test_clausius()
      real(wp) :: ana, fd
      real(wp), parameter :: dt = 1.0e-3_wp
      print '(a)', 'test_clausius:'
      ana = d_sat_vapor_pressure_dt(295.0_wp)
      fd  = (sat_vapor_pressure(295.0_wp + dt) - sat_vapor_pressure(295.0_wp - dt)) / (2.0_wp * dt)
      call check('d(e_sat)/dT vs FD', ana, fd, 1.0e-3_wp * abs(fd) + 1.0e-6_wp)
   end subroutine test_clausius

   subroutine soil_setup(soil, therm, forcing)
      type(soil_params_t),         intent(out) :: soil
      type(soil_thermal_params_t), intent(out) :: therm
      type(energy_forcing_t),      intent(out) :: forcing
      call build_soil_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,      &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, soil)
      call build_soil_thermal(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, therm)
      forcing%soil_water(1:10) = 0.30_wp ; forcing%w_flux = 0.0_wp ; forcing%root_heat_sink = 0.0_wp
      forcing%g_top = 0.0_wp ; forcing%geothermal = 0.0_wp
   end subroutine soil_setup

   subroutine init_col(col, therm, forcing, t_init)
      type(soil_energy_column_t),  intent(out) :: col
      type(soil_thermal_params_t), intent(in)  :: therm
      type(energy_forcing_t),      intent(in)  :: forcing
      real(wp),                    intent(in)  :: t_init
      integer(ik) :: k
      do k = 1_ik, 10_ik
         col%soil_energy(k) = temp_to_uext(therm%soil_dry_heat_capacity(k),                    &
                              forcing%soil_water(k) * rho_h2o, t_init, 1.0_wp)
      end do
   end subroutine init_col

   subroutine test_soil_conserve()
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: therm
      type(energy_forcing_t)      :: forcing
      type(soil_energy_column_t)  :: col
      type(energy_opts_t)         :: opts
      type(energy_flux_t)         :: flux
      integer(ik) :: step
      real(wp)    :: worst, worst_adv
      print '(a)', 'test_soil_conserve:'
      call soil_setup(soil, therm, forcing) ; call init_col(col, therm, forcing, 285.0_wp)
      forcing%g_top = 100.0_wp                                    ! W/m2 surface warming
      worst = 0.0_wp
      do step = 1_ik, 20_ik
         call soil_energy_flux(col, forcing, therm, soil, opts, 1800.0_wp, flux)
         worst = max(worst, abs(flux%energy_resid))
      end do
      call check_true('soil energy residual ~ 0', worst < 1.0e-4_wp, worst)
      call check_true('top layer warmed by G_top', col%soil_temp(1) > 285.0_wp, col%soil_temp(1))
      !----- With inter-layer advective water flux (energy still conserves). --------------!
      call init_col(col, therm, forcing, 285.0_wp)
      forcing%w_flux(1:9) = -1.0e-7_wp                            ! gentle downward flow
      worst_adv = 0.0_wp
      do step = 1_ik, 20_ik
         call soil_energy_flux(col, forcing, therm, soil, opts, 1800.0_wp, flux)
         worst_adv = max(worst_adv, abs(flux%energy_resid))
      end do
      call check_true('soil energy residual ~ 0 (advective)', worst_adv < 1.0e-4_wp, worst_adv)
   end subroutine test_soil_conserve

   subroutine test_soil_steady()
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: therm
      type(energy_forcing_t)      :: forcing
      type(soil_energy_column_t)  :: col
      type(energy_opts_t)         :: opts
      type(energy_flux_t)         :: flux
      integer(ik) :: step, k
      real(wp)    :: drift
      print '(a)', 'test_soil_steady:'
      call soil_setup(soil, therm, forcing) ; call init_col(col, therm, forcing, 285.0_wp)
      ! g_top = geothermal = sink = w_flux = 0  =>  isothermal column must not drift
      do step = 1_ik, 10_ik
         call soil_energy_flux(col, forcing, therm, soil, opts, 1800.0_wp, flux)
      end do
      drift = 0.0_wp
      do k = 1_ik, 10_ik
         drift = max(drift, abs(col%soil_temp(k) - 285.0_wp))
      end do
      call check_true('sealed isothermal column steady', drift < 1.0e-6_wp, drift)
   end subroutine test_soil_steady

end program test_column_energy
