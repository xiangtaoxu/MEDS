!==========================================================================================!
! test_snow -- unit tests for the MEDS snow store (meds_snow_mass + meds_snow_energy, P0).      !
!   1. snow-cover fraction is monotone in [0,1].                                                 !
!   2. datum round-trip: temp_to_uext/uext_to_temp + internal_energy_ice are exact inverses.     !
!   3. cold-snap ACCUMULATION: swe grows, snow_temp < t_3ple, fliq = 0; mass+energy ledgers close.!
!   4. MELT event: snow_temp pins at t_3ple, fliq rises, free liquid drains as melt->infiltration, !
!      swe decreases to melt-out; both ledgers close through the plateau.                          !
!   5. ISOTHERMAL cold pack: with no radiative/turbulent forcing + t_base=t_soil, state is steady  !
!      and g_base ~ 0.                                                                             !
!   6. RAIN-ON-SNOW: warm rain on a sub-freezing pack adds mass + refreezes (energy rises).         !
!==========================================================================================!
program test_snow
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_3ple
   use meds_thermo,             only : temp_to_uext, uext_to_temp, internal_energy_ice,            &
                                       internal_energy_liquid, sat_specific_humidity
   use meds_column_state_types, only : snow_column_t
   use meds_biophysics_types,   only : snow_params_t, snow_env_t, snow_flux_t, snow_melt_t
   use meds_snow_mass,          only : snow_cover_fraction, snow_accumulate, snow_drain_meltwater
   use meds_snow_energy,        only : snow_energy_step
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_snow_cover()
   call test_datum_roundtrip()
   call test_accumulate()
   call test_melt_event()
   call test_isothermal()
   call test_rain_on_snow()

   if (nfail == 0_ik) then
      print '(a)', 'test_snow: ALL PASSED'
   else
      print '(a,i0,a)', 'test_snow: ', nfail, ' FAILED'
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

   !----- Cold-night forcing (accumulation): no SW, net LW loss, cold dry CAS, warmer soil below. -!
   subroutine cold_env(env)
      type(snow_env_t), intent(out) :: env
      env%abs_sw = 0.0_wp ; env%abs_lw = -40.0_wp
      env%can_temp = 265.0_wp ; env%can_shv = 1.0e-3_wp
      env%ggnet = 0.01_wp ; env%rho_air = 1.3_wp ; env%press = 101325.0_wp
      env%t_soil_top = 273.0_wp ; env%k_soil_top = 1.5_wp ; env%dz_soil_top = 0.05_wp
   end subroutine cold_env

   subroutine default_params(p)
      type(snow_params_t), intent(out) :: p
      p = snow_params_t()
   end subroutine default_params

   !=====================================================================================!
   subroutine test_snow_cover()
      type(snow_params_t) :: p
      real(wp) :: f0, f1, f2
      print '(a)', 'test_snow_cover:'
      call default_params(p)
      f0 = snow_cover_fraction(0.0_wp,  0.0_wp,  p)
      f1 = snow_cover_fraction(1.0_wp,  1.0_wp/p%rho_snow,  p)
      f2 = snow_cover_fraction(50.0_wp, 50.0_wp/p%rho_snow, p)
      call check('snowfac(0)=0', f0, 0.0_wp, 1.0e-12_wp)
      call check_true('snowfac monotone 0<f1<f2<=1', f1 > 0.0_wp .and. f1 < f2 .and. f2 <= 1.0_wp, f2)
      call check_true('deep pack fully covered (~1)', f2 > 0.99_wp, f2)
   end subroutine test_snow_cover

   subroutine test_datum_roundtrip()
      real(wp) :: u_ice, u_liq, t, fl
      print '(a)', 'test_datum_roundtrip:'
      ! ice at 260 K: temp_to_uext(0,swe,T,0) == swe*internal_energy_ice(T), inverts back to T
      u_ice = temp_to_uext(0.0_wp, 5.0_wp, 260.0_wp, 0.0_wp)
      call check('temp_to_uext(ice) = swe*iu_ice', u_ice, 5.0_wp * internal_energy_ice(260.0_wp), 1.0e-6_wp)
      call uext_to_temp(u_ice, 5.0_wp, 0.0_wp, t, fl)
      call check('ice round-trip temperature', t, 260.0_wp, 1.0e-9_wp)
      call check('ice round-trip fliq=0', fl, 0.0_wp, 1.0e-12_wp)
      ! melting 5 kg ice -> liquid at t_3ple costs latent_heat_fusion (u_liq - u_ice at t_3ple)
      u_ice = temp_to_uext(0.0_wp, 5.0_wp, t_3ple, 0.0_wp)
      u_liq = temp_to_uext(0.0_wp, 5.0_wp, t_3ple, 1.0_wp)
      call check('fusion cost = swe*L_f', u_liq - u_ice, 5.0_wp * 3.34e5_wp, 1.0e-3_wp)
   end subroutine test_datum_roundtrip

   !----- Build a cold pack by accumulation; assert growth, sub-freezing, and closed ledgers. -----!
   subroutine test_accumulate()
      type(snow_column_t) :: snow
      type(snow_env_t)    :: env
      type(snow_params_t) :: p
      type(snow_flux_t)   :: fx
      type(snow_melt_t)   :: melt
      real(wp) :: dt, snowf, mass_in, mass_out, e_in, e_out, worst_resid
      integer(ik) :: k
      print '(a)', 'test_accumulate:'
      call default_params(p) ; call cold_env(env)
      snow = snow_column_t()
      dt = 900.0_wp ; snowf = 1.0e-3_wp                 ! [kg/m2/s] ~ 3.6 mm SWE/hr
      mass_in = 0.0_wp ; mass_out = 0.0_wp ; e_in = 0.0_wp ; e_out = 0.0_wp ; worst_resid = 0.0_wp
      do k = 1_ik, 40_ik
         call snow_accumulate(snow, snowf, 0.0_wp, 264.0_wp, dt, p)
         mass_in = mass_in + snowf * dt
         e_in    = e_in    + snowf * dt * internal_energy_ice(min(t_3ple, 264.0_wp))
         call snow_energy_step(snow, env, p, dt, 1.0_wp, fx)
         mass_out = mass_out + fx%w_flux * dt
         e_in     = e_in     + dt * (fx%rnet - fx%h_snow - fx%le_snow - fx%g_base)
         worst_resid = max(worst_resid, abs(fx%energy_resid))
         call snow_drain_meltwater(snow, p, melt)
         mass_out = mass_out + melt%melt_mass + melt%dump_mass
         e_out    = e_out    + melt%melt_enth + melt%dump_enth
      end do
      call check_true('pack grew (swe > 30 kg/m2)', snow%swe(1) > 30.0_wp, snow%swe(1))
      call check_true('snow sub-freezing (temp < t_3ple)', snow%snow_temp(1) < t_3ple, snow%snow_temp(1))
      call check_true('fully frozen (fliq = 0)', snow%snow_fliq(1) < 1.0e-9_wp, snow%snow_fliq(1))
      call check('MASS closure: swe = in - out', snow%swe(1), mass_in - mass_out, 1.0e-9_wp)
      call check('ENERGY closure: energy = in - out', snow%snow_energy(1), e_in - e_out, 1.0e-4_wp)
      call check('per-step energy_resid ~ 0', worst_resid, 0.0_wp, 1.0e-6_wp)
   end subroutine test_accumulate

   !----- Warm a cold pack: temp pins at t_3ple, fliq rises, melt drains, swe -> melt-out. --------!
   subroutine test_melt_event()
      type(snow_column_t) :: snow
      type(snow_env_t)    :: env
      type(snow_params_t) :: p
      type(snow_flux_t)   :: fx
      type(snow_melt_t)   :: melt
      real(wp) :: dt, swe0, mass_in, mass_out, e_in, e_out, plateau_temp, tot_melt
      integer(ik) :: k
      logical :: saw_plateau
      print '(a)', 'test_melt_event:'
      call default_params(p) ; call cold_env(env)
      ! seed a 20 kg/m2 cold pack directly (all ice at 270 K)
      snow = snow_column_t()
      snow%nlayer = 1_ik ; snow%swe(1) = 20.0_wp
      snow%snow_energy(1) = temp_to_uext(0.0_wp, 20.0_wp, 270.0_wp, 0.0_wp)
      snow%snow_depth(1)  = 20.0_wp / p%rho_snow
      swe0 = snow%swe(1)
      ! strong warm radiative + warm CAS forcing to drive melt
      env%abs_sw = 400.0_wp ; env%abs_lw = 20.0_wp
      env%can_temp = 280.0_wp ; env%can_shv = 5.0e-3_wp ; env%t_soil_top = 275.0_wp
      dt = 900.0_wp
      mass_in = 0.0_wp ; mass_out = 0.0_wp ; e_in = 0.0_wp ; e_out = 0.0_wp ; tot_melt = 0.0_wp
      e_in = snow%snow_energy(1)   ! initial store folded into the ledger baseline
      saw_plateau = .false. ; plateau_temp = 0.0_wp
      do k = 1_ik, 120_ik
         call snow_energy_step(snow, env, p, dt, 1.0_wp, fx)
         mass_out = mass_out + fx%w_flux * dt
         e_in     = e_in     + dt * (fx%rnet - fx%h_snow - fx%le_snow - fx%g_base)
         if (snow%nlayer == 1_ik .and. snow%snow_fliq(1) > 0.05_wp .and. snow%snow_fliq(1) < 0.95_wp) then
            saw_plateau = .true. ; plateau_temp = snow%snow_temp(1)
         end if
         call snow_drain_meltwater(snow, p, melt)
         mass_out = mass_out + melt%melt_mass + melt%dump_mass
         tot_melt = tot_melt + melt%melt_mass
         e_out    = e_out    + melt%melt_enth + melt%dump_enth
         if (snow%nlayer == 0_ik) exit
      end do
      call check_true('melted out (nlayer = 0)', snow%nlayer == 0_ik, real(snow%nlayer, wp))
      call check_true('saw the melt plateau (0<fliq<1)', saw_plateau, merge(1.0_wp, 0.0_wp, saw_plateau))
      call check('plateau temperature = t_3ple', plateau_temp, t_3ple, 1.0e-6_wp)
      call check_true('meltwater drained to soil (> half of swe0)', tot_melt > 0.5_wp * swe0, tot_melt)
      call check('MASS closure: 0 = swe0 - out', 0.0_wp, swe0 - mass_out, 1.0e-7_wp)
      call check('ENERGY closure: 0 = in - out', 0.0_wp, e_in - e_out, 1.0e-2_wp)
   end subroutine test_melt_event

   !----- Isothermal cold pack in radiative/turbulent equilibrium with the CAS + soil: steady. ----!
   subroutine test_isothermal()
      type(snow_column_t) :: snow
      type(snow_env_t)    :: env
      type(snow_params_t) :: p
      type(snow_flux_t)   :: fx
      type(snow_melt_t)   :: melt
      real(wp) :: swe0, e0, tsnow0
      integer(ik) :: k
      print '(a)', 'test_isothermal:'
      call default_params(p)
      snow = snow_column_t()
      snow%nlayer = 1_ik ; snow%swe(1) = 15.0_wp
      snow%snow_energy(1) = temp_to_uext(0.0_wp, 15.0_wp, 268.0_wp, 0.0_wp)
      snow%snow_depth(1)  = 15.0_wp / p%rho_snow
      ! forcing that exactly matches the snow surface (no gradients): CAS at snow temp + saturated, soil at snow temp
      env%abs_sw = 0.0_wp ; env%abs_lw = 0.0_wp
      env%can_temp = 268.0_wp ; env%ggnet = 0.01_wp ; env%rho_air = 1.3_wp ; env%press = 101325.0_wp
      env%t_soil_top = 268.0_wp ; env%k_soil_top = 1.5_wp ; env%dz_soil_top = 0.05_wp
      env%can_shv = sat_specific_humidity(268.0_wp, env%press)   ! SATURATED CAS -> no vapour gradient -> no sublimation
      swe0 = snow%swe(1) ; e0 = snow%snow_energy(1) ; tsnow0 = 268.0_wp
      do k = 1_ik, 50_ik
         call snow_energy_step(snow, env, p, 900.0_wp, 1.0_wp, fx)
         call snow_drain_meltwater(snow, p, melt)
      end do
      call check_true('g_base ~ 0 (no soil gradient)', abs(fx%g_base) < 1.0e-6_wp, fx%g_base)
      call check_true('swe steady (drift < 0.5 kg/m2)', abs(snow%swe(1) - swe0) < 0.5_wp, snow%swe(1) - swe0)
      call check_true('snow_temp steady (within 1 K)', abs(snow%snow_temp(1) - tsnow0) < 1.0_wp, snow%snow_temp(1))
      call check_true('stayed frozen (no spurious melt)', snow%snow_fliq(1) < 1.0e-6_wp, snow%snow_fliq(1))
   end subroutine test_isothermal

   !----- Rain on a sub-freezing pack: adds mass + liquid enthalpy, refreezes (energy rises). -----!
   subroutine test_rain_on_snow()
      type(snow_column_t) :: snow
      type(snow_params_t) :: p
      real(wp) :: e0, swe0
      print '(a)', 'test_rain_on_snow:'
      call default_params(p)
      snow = snow_column_t()
      snow%nlayer = 1_ik ; snow%swe(1) = 10.0_wp
      snow%snow_energy(1) = temp_to_uext(0.0_wp, 10.0_wp, 265.0_wp, 0.0_wp)   ! cold pack
      snow%snow_depth(1)  = 10.0_wp / p%rho_snow
      e0 = snow%snow_energy(1) ; swe0 = snow%swe(1)
      call snow_accumulate(snow, 0.0_wp, 2.0e-3_wp, 278.0_wp, 900.0_wp, p)    ! warm rain (2 mm/900s)
      call check('rain adds mass', snow%swe(1), swe0 + 2.0e-3_wp * 900.0_wp, 1.0e-9_wp)
      call check_true('rain adds energy (warms/refreezes pack)', snow%snow_energy(1) > e0, snow%snow_energy(1) - e0)
      call check('rain enthalpy exact', snow%snow_energy(1) - e0,                                   &
                 2.0e-3_wp * 900.0_wp * internal_energy_liquid(278.0_wp), 1.0e-6_wp)
   end subroutine test_rain_on_snow

end program test_snow
