!==========================================================================================!
! test_column_energy -- unit tests for the soil THERMAL column + shared enthalpy inverter.    !
!   1. INVERTER round-trip (liquid + ice) and CONTINUITY at u_freeze/u_melt (temp = t_3ple).   !
!   2. CLAUSIUS slope d(e_sat)/dT vs finite difference.                                         !
!   3. SOIL energy CONSERVATION to ~round-off (a G_top warming; and with advective w_flux).     !
!   4. STEADY STATE: a sealed isothermal column does not drift.                                 !
!   5. ICE-AWARE conductivity: frozen saturated soil conducts more than liquid (k_ice>k_water). !
!   6. FREEZE plateau (P2a zero-curtain): cooling a wet layer pins soil_temp at t_3ple while     !
!      soil_fliq falls 1->0, absorbing exactly wmass*L_f; energy conserved through the plateau.   !
!   7. THAW plateau (P2a): warming a frozen layer pins at t_3ple while soil_fliq rises 0->1.       !
!==========================================================================================!
program test_column_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, t_3ple, cp_ice, latent_heat_fusion, k_water, k_ice
   use meds_biophysics_types, only : soil_energy_column_t, energy_forcing_t, soil_thermal_params_t, &
                                     soil_params_t, energy_opts_t, energy_flux_t, SOIL_RETENTION_VG, &
                                     ENERGY_PHASE_ON
   use meds_column_state_types, only : build_soil_hydr_params
   use meds_column_state_types, only : build_soil_therm_params
   use meds_therm_lib,           only : soil_thermal_cond
   use meds_therm_lib,           only : temp_to_uext, uext_to_temp, sat_vapor_pressure,           &
                                     sat_vapor_pressure_temp_deriv, internal_energy_liquid
   use meds_soil_energy,      only : soil_energy_step_implicit
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_inverter()
   call test_clausius()
   call test_soil_conserve()
   call test_soil_steady()
   call test_ice_conductivity()
   call test_freeze_plateau()
   call test_thaw_plateau()
   call test_mass_correction_neutrality()

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
      ana = sat_vapor_pressure_temp_deriv(295.0_wp)
      fd  = (sat_vapor_pressure(295.0_wp + dt) - sat_vapor_pressure(295.0_wp - dt)) / (2.0_wp * dt)
      call check('d(e_sat)/dT vs FD', ana, fd, 1.0e-3_wp * abs(fd) + 1.0e-6_wp)
   end subroutine test_clausius

   subroutine soil_setup(soil, therm, forcing)
      type(soil_params_t),         intent(out) :: soil
      type(soil_thermal_params_t), intent(out) :: therm
      type(energy_forcing_t),      intent(out) :: forcing
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,      &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, soil)
      call build_soil_therm_params(10_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, therm)
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
         call soil_energy_step_implicit(col, forcing, therm, soil, opts, 1800.0_wp, flux)
         worst = max(worst, abs(flux%energy_resid))
      end do
      call check_true('soil energy residual ~ 0', worst < 1.0e-4_wp, worst)
      call check_true('top layer warmed by G_top', col%soil_temp(1) > 285.0_wp, col%soil_temp(1))
      !----- With inter-layer advective water flux (energy still conserves). --------------!
      call init_col(col, therm, forcing, 285.0_wp)
      forcing%w_flux(1:9) = -1.0e-7_wp                            ! gentle downward flow
      worst_adv = 0.0_wp
      do step = 1_ik, 20_ik
         call soil_energy_step_implicit(col, forcing, therm, soil, opts, 1800.0_wp, flux)
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
         call soil_energy_step_implicit(col, forcing, therm, soil, opts, 1800.0_wp, flux)
      end do
      drift = 0.0_wp
      do k = 1_ik, 10_ik
         drift = max(drift, abs(col%soil_temp(k) - 285.0_wp))
      end do
      call check_true('sealed isothermal column steady', drift < 1.0e-6_wp, drift)
   end subroutine test_soil_steady

   !----- Ice-aware saturated conductivity: frozen soil conducts more (k_ice > k_water). -----!
   subroutine test_ice_conductivity()
      real(wp) :: k_liq, k_frz, ksat_liq, ksat_frz, ke_liq, ke_frz, sr
      real(wp), parameter :: ts = 0.5_wp, th = 0.15_wp, ksol = 3.0_wp, kdry = 0.15_wp
      print '(a)', 'test_ice_conductivity:'
      k_liq = soil_thermal_cond(0.40_wp, 1.0_wp, 0.43_wp, 3.0_wp, 0.15_wp)   ! all-liquid saturation
      k_frz = soil_thermal_cond(0.40_wp, 0.0_wp, 0.43_wp, 3.0_wp, 0.15_wp)   ! all-ice saturation
      call check_true('frozen soil conducts more than liquid', k_frz > k_liq, k_frz - k_liq)
      !----- Kersten number: Johansen log10 form when liquid, Farouki linear s_r form when frozen. --!
      sr       = th / ts                                                     ! 0.30 (intermediate sat.)
      ksat_liq = ksol ** (1.0_wp - ts) * k_water ** ts                       ! all-liquid geometric mean
      ksat_frz = ksol ** (1.0_wp - ts) * k_ice   ** ts                       ! all-ice geometric mean
      ke_liq   = (soil_thermal_cond(th, 1.0_wp, ts, ksol, kdry) - kdry) / (ksat_liq - kdry)
      ke_frz   = (soil_thermal_cond(th, 0.0_wp, ts, ksol, kdry) - kdry) / (ksat_frz - kdry)
      call check('unfrozen Kersten is the log10 form',  ke_liq, log10(sr) + 1.0_wp, 1.0e-9_wp)
      call check('frozen Kersten is the linear s_r form', ke_frz, sr,               1.0e-9_wp)
   end subroutine test_ice_conductivity

   !----- Single wet control volume, sealed bottom; the top Neumann G_top is the ONLY flux, so !
   !      energy removed per step is exactly |G_top|*dt regardless of the BE predictor.         !
   subroutine setup_1layer(soil, therm, forcing, theta, depth)
      type(soil_params_t),         intent(out) :: soil
      type(soil_thermal_params_t), intent(out) :: therm
      type(energy_forcing_t),      intent(out) :: forcing
      real(wp),                    intent(in)  :: theta, depth
      call build_soil_hydr_params(1_ik, SOIL_RETENTION_VG, depth, 3.0_wp, 0.43_wp, 0.078_wp,        &
           2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, soil)
      call build_soil_therm_params(1_ik, 3.0_wp, 0.15_wp, 2.0e6_wp, therm)
      forcing%soil_water(1) = theta ; forcing%w_flux = 0.0_wp ; forcing%root_heat_sink = 0.0_wp
      forcing%g_top = 0.0_wp ; forcing%geothermal = 0.0_wp
   end subroutine setup_1layer

   subroutine test_freeze_plateau()
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: therm
      type(energy_forcing_t)      :: forcing
      type(soil_energy_column_t)  :: col
      type(energy_opts_t)         :: opts
      type(energy_flux_t)         :: flux
      real(wp), parameter :: theta = 0.40_wp, depth = 0.05_wp
      real(wp) :: latent_expect, e_enter, e_exit, worst_pin, worst_res, fl_prev, fl, tp, e_col
      integer(ik) :: step
      logical :: got_enter, got_exit
      print '(a)', 'test_freeze_plateau:'
      call setup_1layer(soil, therm, forcing, theta, depth)
      opts%phase_change = ENERGY_PHASE_ON
      col%soil_energy(1) = temp_to_uext(therm%soil_dry_heat_capacity(1), theta * rho_h2o,       &
                                        t_3ple + 0.2_wp, 1.0_wp)               ! all-liquid, just above 0 C
      forcing%g_top = -100.0_wp                                               ! steady surface cooling
      latent_expect = theta * rho_h2o * latent_heat_fusion * soil%dz(1)       ! [J/m2] to freeze the layer
      worst_pin = 0.0_wp ; worst_res = 0.0_wp ; fl_prev = 1.0_wp
      e_enter = 0.0_wp ; e_exit = 0.0_wp ; got_enter = .false. ; got_exit = .false.
      do step = 1_ik, 100_ik
         call soil_energy_step_implicit(col, forcing, therm, soil, opts, 1800.0_wp, flux)
         worst_res = max(worst_res, abs(flux%energy_resid))
         fl = col%soil_fliq(1) ; tp = col%soil_temp(1) ; e_col = col%soil_energy(1) * soil%dz(1)
         if (fl > 1.0e-3_wp .and. fl < 1.0_wp - 1.0e-3_wp) worst_pin = max(worst_pin, abs(tp - t_3ple))
         if (.not. got_enter .and. fl < 1.0_wp - 1.0e-3_wp) then ; e_enter = e_col ; got_enter = .true. ; end if
         if (.not. got_exit  .and. fl < 1.0e-3_wp)          then ; e_exit  = e_col ; got_exit  = .true. ; end if
         fl_prev = fl
      end do
      call check_true('freeze: soil_temp pinned at t_3ple on plateau', worst_pin < 1.0e-6_wp, worst_pin)
      call check_true('freeze: energy residual ~ 0 through plateau',   worst_res < 1.0e-4_wp, worst_res)
      call check_true('freeze: layer fully frozen (fliq -> 0)',        got_exit,              col%soil_fliq(1))
      call check_true('freeze: all-ice layer cools below t_3ple',      col%soil_temp(1) < t_3ple, col%soil_temp(1))
      call check('freeze: latent heat absorbed = wmass*L_f', e_enter - e_exit, latent_expect,  &
                 0.08_wp * latent_expect)                                     ! within ~1 step granularity
   end subroutine test_freeze_plateau

   subroutine test_thaw_plateau()
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: therm
      type(energy_forcing_t)      :: forcing
      type(soil_energy_column_t)  :: col
      type(energy_opts_t)         :: opts
      type(energy_flux_t)         :: flux
      real(wp), parameter :: theta = 0.40_wp, depth = 0.05_wp
      real(wp) :: worst_pin, fl, tp
      integer(ik) :: step
      logical :: melted
      print '(a)', 'test_thaw_plateau:'
      call setup_1layer(soil, therm, forcing, theta, depth)
      opts%phase_change = ENERGY_PHASE_ON
      col%soil_energy(1) = temp_to_uext(therm%soil_dry_heat_capacity(1), theta * rho_h2o,       &
                                        t_3ple - 0.2_wp, 0.0_wp)              ! all-ice, just below 0 C
      forcing%g_top = 100.0_wp                                               ! steady surface warming
      worst_pin = 0.0_wp ; melted = .false.
      do step = 1_ik, 100_ik
         call soil_energy_step_implicit(col, forcing, therm, soil, opts, 1800.0_wp, flux)
         fl = col%soil_fliq(1) ; tp = col%soil_temp(1)
         if (fl > 1.0e-3_wp .and. fl < 1.0_wp - 1.0e-3_wp) worst_pin = max(worst_pin, abs(tp - t_3ple))
         if (fl > 1.0_wp - 1.0e-3_wp) melted = .true.
      end do
      call check_true('thaw: soil_temp pinned at t_3ple on plateau', worst_pin < 1.0e-6_wp, worst_pin)
      call check_true('thaw: layer fully melted (fliq -> 1)',        melted,                col%soil_fliq(1))
      call check_true('thaw: all-liquid layer warms above t_3ple',   col%soil_temp(1) > t_3ple, col%soil_temp(1))
   end subroutine test_thaw_plateau

   !=======================================================================================!
   !  TEMPERATURE-NEUTRALITY of an unfaced mass correction -- the identity the driver's       !
   !  clip / theta_res-floor compensation rests on (meds_fast_split.f90 sec 3d').              !
   !                                                                                          !
   !  The hydrology applies post-solve corrections to theta that no face accounts for. The      !
   !  energy column then inverts the CORRECTED water mass against the OLD internal energy, so    !
   !  an uncompensated correction lands entirely in the diagnosed temperature. It is large:      !
   !  internal_energy_liquid carries the tsupercool_liq datum, so liquid water is ~1.0 MJ/kg in   !
   !  ABSOLUTE terms and only differences are physical. Moving the matching enthalpy with the      !
   !  mass, at the layer's OWN temperature, must leave the temperature exactly where it was --     !
   !  a numerical fix-up may not heat or cool the soil.                                             !
   !                                                                                                !
   !  Both branches are asserted: uncompensated must MOVE the temperature (so the test would          !
   !  fail if the datum were ever rebased to make u_liq ~ 0 and the bug went quiet), compensated       !
   !  must not. --------------------------------------------------------------------------------------!
   subroutine test_mass_correction_neutrality()
      real(wp) :: dry_hcap, wmass0, wmass1, dw, temp0, temp_raw, temp_fix, fliq0, fliq, uext0, u_liq
      print '(a)', 'test_mass_correction_neutrality:'
      dry_hcap = 1.3e6_wp                                   ! [J/m3/K] typical dry-soil volumetric capacity
      temp0    = 291.0_wp
      wmass0   = 0.40_wp * rho_h2o                          ! [kg/m3] theta = 0.40 in the layer
      uext0    = temp_to_uext(dry_hcap, wmass0, temp0, 1.0_wp)
      call uext_to_temp(uext0, wmass0, dry_hcap, temp0, fliq0)   ! exact round-trip seed

      !----- A 1 kg/m2 clip out of a 0.1 m layer: dtheta = -0.01, i.e. -10 kg/m3. -----------------!
      dw     = -10.0_wp
      wmass1 = wmass0 + dw

      !----- (a) mass removed, energy untouched: the layer's temperature jumps. -------------------!
      call uext_to_temp(uext0, wmass1, dry_hcap, temp_raw, fliq)
      call check_true('uncompensated clip visibly moves soil temperature',                        &
                      abs(temp_raw - temp0) > 1.0_wp, temp_raw - temp0)

      !----- (b) mass AND its enthalpy removed at the layer's own temperature: T is unchanged. ----!
      !      Uses the very function the driver's compensation calls. ------------------------------!
      u_liq = internal_energy_liquid(temp0)                      ! [J/kg] ~1.0e6, NOT ~0 -- that is the point
      call uext_to_temp(uext0 + dw * u_liq, wmass1, dry_hcap, temp_fix, fliq)
      call check('compensated clip leaves soil temperature exactly unchanged', temp_fix, temp0,   &
                 1.0e-9_wp)
      call check_true('compensated clip leaves the layer all-liquid', abs(fliq - 1.0_wp) < 1.0e-12_wp, fliq)

      !----- Same for the OPPOSITE sign (the theta_res hard floor CREATES water in the layer). ----!
      dw     = 10.0_wp
      wmass1 = wmass0 + dw
      call uext_to_temp(uext0 + dw * u_liq, wmass1, dry_hcap, temp_fix, fliq)
      call check('compensated floor add is temperature-neutral too', temp_fix, temp0, 1.0e-9_wp)
   end subroutine test_mass_correction_neutrality

end program test_column_energy
