!==========================================================================================!
! test_column_hydrology -- unit tests for the soil-water column (src/biophysics).           !
!                                                                                          !
! Physical invariants the MVP must satisfy:                                                  !
!   1. CONSTITUTIVE ROUND-TRIP : theta_of_psi(psi_of_theta(theta)) == theta (vG + Campbell);   !
!                                C(psi) matches the finite-difference dtheta/dpsi.              !
!   2. MASS CONSERVATION        : total-store balance closes to ~round-off every step.          !
!   3. NO-FLUX (bedrock)         : a sealed column with no sources conserves total water.         !
!   4. FREE DRAINAGE             : a wet column drains (drainage > 0, water decreases), closes.    !
!   5. INTERCEPTION              : per-cohort bucket conserves; capacity respected, drip on top.    !
!   6. INFILTRATION CAP          : heavy rain on a low-K soil is conductivity-limited; excess ponds. !
!==========================================================================================!
program test_column_hydrology
   use meds_kinds,            only : wp, ik
   use meds_biophysics_types, only : soil_column_t, chydro_forcing_t, soil_params_t,          &
                                     soil_opts_t, chydro_flux_t, n_soil_layer_max,             &
                                     SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL,               &
                                     SOIL_BC_FREE_DRAIN, SOIL_BC_BEDROCK
   use meds_soil_parameters,  only : soil_theta_of_psi, soil_psi_of_theta, soil_moist_cap,     &
                                     build_soil_params
   use meds_column_hydrology, only : column_hydrology_flux, intercept_canopy_layer
   implicit none

   integer(ik) :: nfail
   nfail = 0_ik

   call test_constitutive()
   call test_mass_conservation()
   call test_bedrock_conserve()
   call test_free_drain()
   call test_interception()
   call test_infiltration_cap()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_hydrology: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_hydrology: ', nfail, ' FAILED'
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
         print '(a,a,a,es13.5,a,es13.5,a,es10.2)', '  FAIL : ', name, '  got ', got,          &
               ' expected ', expect, '  |diff|>', atol
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (', val, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

   !----- A 10-layer loam column (2 m) with a chosen retention curve. -------------------------!
   subroutine loam_column(retention, params, col)
      integer(ik),         intent(in)  :: retention
      type(soil_params_t), intent(out) :: params
      type(soil_column_t), intent(out) :: col
      if (retention == SOIL_RETENTION_CAMPBELL) then
         call build_soil_params(10_ik, retention, 2.0_wp, 3.0_wp, 0.44_wp, 0.0_wp,            &
              4.53e-6_wp, -0.26_wp, 5.65_wp, 2.0_wp, -3.37_wp, params)
      else
         call build_soil_params(10_ik, retention, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,          &
              2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, params)
      end if
      col%theta(1:10) = 0.30_wp
      col%w_surface = 0.0_wp ; col%w_aquifer = 0.0_wp ; col%z_wt = 0.0_wp
   end subroutine loam_column

   !=======================================================================================!
   subroutine test_constitutive()
      real(wp) :: th, ps, th_rt, c_ana, c_fd, dps
      print '(a)', 'test_constitutive:'
      !----- van Genuchten loam round-trip + capacity. -----!
      th    = 0.30_wp
      ps    = soil_psi_of_theta(SOIL_RETENTION_VG, th, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      th_rt = soil_theta_of_psi(SOIL_RETENTION_VG, ps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      call check('vG theta round-trip', th_rt, th, 1.0e-9_wp)
      dps  = 1.0e-3_wp
      c_ana = soil_moist_cap(SOIL_RETENTION_VG, ps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      c_fd  = (soil_theta_of_psi(SOIL_RETENTION_VG, ps+dps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp) &
             - soil_theta_of_psi(SOIL_RETENTION_VG, ps-dps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp))&
             / (2.0_wp * dps)
      call check('vG C = dtheta/dpsi (FD)', c_ana, c_fd, 1.0e-3_wp * abs(c_fd) + 1.0e-9_wp)
      !----- Campbell loam round-trip. -----!
      th    = 0.30_wp
      ps    = soil_psi_of_theta(SOIL_RETENTION_CAMPBELL, th, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      th_rt = soil_theta_of_psi(SOIL_RETENTION_CAMPBELL, ps, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      call check('Campbell theta round-trip', th_rt, th, 1.0e-9_wp)
      c_ana = soil_moist_cap(SOIL_RETENTION_CAMPBELL, ps, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      c_fd  = (soil_theta_of_psi(SOIL_RETENTION_CAMPBELL, ps+dps, 0.44_wp,0.0_wp,-0.26_wp,5.65_wp)&
             - soil_theta_of_psi(SOIL_RETENTION_CAMPBELL, ps-dps, 0.44_wp,0.0_wp,-0.26_wp,5.65_wp))&
             / (2.0_wp * dps)
      call check('Campbell C = dtheta/dpsi (FD)', c_ana, c_fd, 1.0e-2_wp * abs(c_fd) + 1.0e-9_wp)
   end subroutine test_constitutive

   !=======================================================================================!
   subroutine test_mass_conservation()
      type(soil_params_t)  :: params
      type(soil_column_t)  :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)    :: opts
      type(chydro_flux_t)  :: flux
      integer(ik) :: step, k
      real(wp)    :: worst
      print '(a)', 'test_mass_conservation:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      forcing%precip_ground = 5.0e-5_wp
      do k = 1_ik, 10_ik
         forcing%root_uptake(k) = 2.0e-5_wp * params%root_frac(k)
      end do
      forcing%t_ground = 298.15_wp ; forcing%q_air = 0.010_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      worst = 0.0_wp
      do step = 1_ik, 40_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
         worst = max(worst, abs(flux%mass_resid))
      end do
      call check_true('mass residual ~ 0 over 40 steps', worst < 1.0e-9_wp, worst)
   end subroutine test_mass_conservation

   !=======================================================================================!
   subroutine test_bedrock_conserve()
      type(soil_params_t)  :: params
      type(soil_column_t)  :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)    :: opts
      type(chydro_flux_t)  :: flux
      real(wp) :: w_before, w_after
      integer(ik) :: step, k
      print '(a)', 'test_bedrock_conserve:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      forcing%precip_ground = 0.0_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp        ! high humidity -> no evap
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_BEDROCK
      w_before = 0.0_wp
      do k = 1_ik, 10_ik
         w_before = w_before + col%theta(k) * params%dz(k)
      end do
      do step = 1_ik, 10_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
      end do
      w_after = 0.0_wp
      do k = 1_ik, 10_ik
         w_after = w_after + col%theta(k) * params%dz(k)
      end do
      call check('sealed column conserves water', w_after, w_before, 1.0e-10_wp)
      call check_true('bedrock drainage = 0', abs(flux%drainage) < 1.0e-30_wp, flux%drainage)
   end subroutine test_bedrock_conserve

   !=======================================================================================!
   subroutine test_free_drain()
      type(soil_params_t)  :: params
      type(soil_column_t)  :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)    :: opts
      type(chydro_flux_t)  :: flux
      real(wp) :: w_before, w_after, worst
      integer(ik) :: step, k
      print '(a)', 'test_free_drain:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      col%theta(1:10) = 0.40_wp                          ! wet, near saturation
      forcing%precip_ground = 0.0_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      w_before = 0.0_wp
      do k = 1_ik, 10_ik
         w_before = w_before + col%theta(k) * params%dz(k)
      end do
      worst = 0.0_wp
      do step = 1_ik, 20_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
         worst = max(worst, abs(flux%mass_resid))
      end do
      w_after = 0.0_wp
      do k = 1_ik, 10_ik
         w_after = w_after + col%theta(k) * params%dz(k)
      end do
      call check_true('free drainage removes water', w_after < w_before, w_before - w_after)
      call check_true('drainage flux > 0', flux%drainage > 0.0_wp, flux%drainage)
      call check_true('mass residual ~ 0 (drain)', worst < 1.0e-9_wp, worst)
   end subroutine test_free_drain

   !=======================================================================================!
   subroutine test_interception()
      real(wp) :: lw, tf, dr, sw, bal, rain, dt
      print '(a)', 'test_interception:'
      dt = 600.0_wp ; rain = 1.0e-4_wp
      !----- Dry canopy fills below capacity: no drip, exact balance. -----!
      lw = 0.0_wp
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, 0.0_wp, dt, 0.1_wp, 0.5_wp, 1.0_wp,&
                                  tf, dr, sw)
      bal = tf + lw / dt + 0.0_wp                       ! throughfall + storage-rate + evap
      call check('interception water balance', bal, rain, 1.0e-12_wp)
      call check_true('no drip below capacity', dr < 1.0e-30_wp, dr)
      call check_true('leaf_water within capacity', lw <= 0.1_wp * 2.5_wp + 1.0e-12_wp, lw)
      !----- Saturated canopy overflows: drip appears, storage capped. -----!
      lw = 0.25_wp                                       ! = dewmx * pai (full)
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, 0.0_wp, dt, 0.1_wp, 0.5_wp, 1.0_wp,&
                                  tf, dr, sw)
      call check_true('drip when full', dr > 0.0_wp, dr)
      call check_true('capacity respected', lw <= 0.25_wp + 1.0e-12_wp, lw)
      call check_true('sigma_w saturates to 1', abs(sw - 1.0_wp) < 1.0e-9_wp, sw)
   end subroutine test_interception

   !=======================================================================================!
   subroutine test_infiltration_cap()
      type(soil_params_t)  :: params
      type(soil_column_t)  :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)    :: opts
      type(chydro_flux_t)  :: flux
      print '(a)', 'test_infiltration_cap:'
      !----- WET clay (low Ksat, small suction) under a downpour: conductivity-limited so     !
      !      infiltration is capped and the excess ponds/runs off (Hortonian). A bone-dry clay !
      !      would instead have huge suction-driven capacity (Green-Ampt) -- not the cap case.  !
      call build_soil_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.38_wp, 0.068_wp,     &
           5.6e-7_wp, 0.8_wp, 1.09_wp, 2.0_wp, -3.37_wp, params)
      col%theta(1:10) = 0.36_wp
      col%w_surface = 0.0_wp ; col%w_aquifer = 0.0_wp ; col%z_wt = 0.0_wp
      forcing%precip_ground = 1.0e-2_wp                  ! 36 mm/hr downpour
      forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      call column_hydrology_flux(col, forcing, params, opts, 60.0_wp, flux)
      call check_true('infiltration < precip (capped)', flux%infiltration < forcing%precip_ground, &
                      flux%infiltration)
      call check_true('excess ponds or runs off',                                             &
                      (col%w_surface + flux%runoff_surf * 60.0_wp) > 0.0_wp,                   &
                      col%w_surface + flux%runoff_surf * 60.0_wp)
      call check_true('mass residual ~ 0 (downpour)', abs(flux%mass_resid) < 1.0e-9_wp,        &
                      abs(flux%mass_resid))
   end subroutine test_infiltration_cap

end program test_column_hydrology
