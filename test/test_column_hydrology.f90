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
                                     SOIL_BC_FREE_DRAIN, SOIL_BC_BEDROCK, SOIL_BC_AQUIFER,     &
                                     SOIL_LIN_PICARD, SOIL_SUBSTEP_FIXED, SOIL_SUBSTEP_ADAPTIVE
   use meds_hydr_lib, only : soil_theta_from_psi, soil_psi_from_theta, soil_moist_cap_from_psi
   use meds_column_state_types, only : build_soil_hydr_params
   use meds_soil_water,       only : column_hydrology_flux
   use meds_vegetation_biophysics, only : intercept_canopy_layer
   implicit none

   integer(ik) :: nfail
   nfail = 0_ik

   call test_constitutive()
   call test_mass_conservation()
   call test_bedrock_conserve()
   call test_free_drain()
   call test_interception()
   call test_infiltration_cap()
   call test_zeng_decker()
   call test_picard()
   call test_adaptive_substep()
   call test_aquifer()
   call test_dunne()
   call test_snow_free_evap()
   call test_clip_layer_decomposition()

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
         call build_soil_hydr_params(10_ik, retention, 2.0_wp, 3.0_wp, 0.44_wp, 0.0_wp,            &
              4.53e-6_wp, -0.26_wp, 5.65_wp, 2.0_wp, -3.37_wp, params)
      else
         call build_soil_hydr_params(10_ik, retention, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,          &
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
      ps    = soil_psi_from_theta(SOIL_RETENTION_VG, th, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      th_rt = soil_theta_from_psi(SOIL_RETENTION_VG, ps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      call check('vG theta round-trip', th_rt, th, 1.0e-9_wp)
      dps  = 1.0e-3_wp
      c_ana = soil_moist_cap_from_psi(SOIL_RETENTION_VG, ps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      c_fd  = (soil_theta_from_psi(SOIL_RETENTION_VG, ps+dps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp) &
             - soil_theta_from_psi(SOIL_RETENTION_VG, ps-dps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp))&
             / (2.0_wp * dps)
      call check('vG C = dtheta/dpsi (FD)', c_ana, c_fd, 1.0e-3_wp * abs(c_fd) + 1.0e-9_wp)
      !----- Campbell loam round-trip. -----!
      th    = 0.30_wp
      ps    = soil_psi_from_theta(SOIL_RETENTION_CAMPBELL, th, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      th_rt = soil_theta_from_psi(SOIL_RETENTION_CAMPBELL, ps, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      call check('Campbell theta round-trip', th_rt, th, 1.0e-9_wp)
      c_ana = soil_moist_cap_from_psi(SOIL_RETENTION_CAMPBELL, ps, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      c_fd  = (soil_theta_from_psi(SOIL_RETENTION_CAMPBELL, ps+dps, 0.44_wp,0.0_wp,-0.26_wp,5.65_wp)&
             - soil_theta_from_psi(SOIL_RETENTION_CAMPBELL, ps-dps, 0.44_wp,0.0_wp,-0.26_wp,5.65_wp))&
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
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.38_wp, 0.068_wp,     &
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

   !=======================================================================================!
   !  P2 numerics + boundary physics                                                        !
   !=======================================================================================!

   !----- Retention-integral Zeng-Decker holds a hydrostatic-equilibrium column at q->0;      !
   !      plain gravity drifts off it (design 3e). -----------------------------------------!
   subroutine test_zeng_decker()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col_zd, col_pg, cdum
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      real(wp) :: z_wt, zc, theta_e, drift_zd, drift_pg
      real(wp) :: theta_init(n_soil_layer_max)
      integer(ik) :: k, j, step, n
      print '(a)', 'test_zeng_decker:'
      call loam_column(SOIL_RETENTION_VG, params, cdum)
      n = 10_ik
      z_wt = params%soil_layer_z(n+1)                       ! water table at column bottom
      do k = 1_ik, n                                        ! layer-mean equilibrium moisture
         theta_e = 0.0_wp
         do j = 1_ik, 5_ik
            zc = params%soil_layer_z(k+1) + (real(j,wp) - 0.5_wp) / 5.0_wp * params%dz(k)
            theta_e = theta_e + soil_theta_from_psi(SOIL_RETENTION_VG, z_wt - zc,                &
               params%theta_sat(k), params%theta_res(k), params%vg_alpha(k), params%vg_n(k)) / 5.0_wp
         end do
         theta_init(k) = theta_e
      end do
      col_zd%theta(1:n) = theta_init(1:n)
      col_zd%w_surface = 0.0_wp ; col_zd%w_aquifer = 0.0_wp ; col_zd%z_wt = 0.0_wp
      col_pg = col_zd
      forcing%precip_ground = 0.0_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_BEDROCK                      ! no-flow (hydrostatic equilibrium test)
      opts%zeng_decker = .true.
      do step = 1_ik, 10_ik
         call column_hydrology_flux(col_zd, forcing, params, opts, 3600.0_wp, flux)
      end do
      opts%zeng_decker = .false.
      do step = 1_ik, 10_ik
         call column_hydrology_flux(col_pg, forcing, params, opts, 3600.0_wp, flux)
      end do
      drift_zd = 0.0_wp ; drift_pg = 0.0_wp
      do k = 1_ik, n
         drift_zd = max(drift_zd, abs(col_zd%theta(k) - theta_init(k)))
         drift_pg = max(drift_pg, abs(col_pg%theta(k) - theta_init(k)))
      end do
      call check_true('ZD holds hydrostatic equilibrium (q->0)', drift_zd < 1.0e-6_wp, drift_zd)
      call check_true('plain gravity drifts off equilibrium', drift_pg > drift_zd, drift_pg)
   end subroutine test_zeng_decker

   !----- Celia modified-Picard: mass closes and every step converges. --------------------!
   subroutine test_picard()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      integer(ik) :: step, k
      real(wp) :: worst
      logical  :: allconv
      print '(a)', 'test_picard:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      forcing%precip_ground = 5.0e-5_wp
      do k = 1_ik, 10_ik
         forcing%root_uptake(k) = 2.0e-5_wp * params%root_frac(k)
      end do
      forcing%t_ground = 298.15_wp ; forcing%q_air = 0.010_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      opts%linearize = SOIL_LIN_PICARD
      worst = 0.0_wp ; allconv = .true.
      do step = 1_ik, 40_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
         worst = max(worst, abs(flux%mass_resid))
         allconv = allconv .and. flux%converged
      end do
      call check_true('Picard mass residual ~ 0', worst < 1.0e-9_wp, worst)
      call check_true('Picard converged every step', allconv, 1.0_wp)
   end subroutine test_picard

   !----- Adaptive step-doubling matches a fixed fine-substep reference. -------------------!
   subroutine test_adaptive_substep()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col_a, col_f
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      integer(ik) :: k
      real(wp) :: diff
      print '(a)', 'test_adaptive_substep:'
      call loam_column(SOIL_RETENTION_VG, params, col_a)
      col_a%theta(1:10) = 0.20_wp                          ! dry-ish -> an infiltration front
      col_f = col_a
      forcing%precip_ground = 2.0e-4_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      opts%substep = SOIL_SUBSTEP_ADAPTIVE ; opts%h_init = 1800.0_wp
      call column_hydrology_flux(col_a, forcing, params, opts, 1800.0_wp, flux)
      opts%substep = SOIL_SUBSTEP_FIXED ; opts%h_init = 10.0_wp       ! ~180 fixed substeps
      call column_hydrology_flux(col_f, forcing, params, opts, 1800.0_wp, flux)
      diff = 0.0_wp
      do k = 1_ik, 10_ik
         diff = max(diff, abs(col_a%theta(k) - col_f%theta(k)))
      end do
      call check_true('adaptive matches fixed-fine reference', diff < 1.0e-3_wp, diff)
   end subroutine test_adaptive_substep

   !----- SIMTOP aquifer: recharge fills the store, the water table rises, mass closes. ----!
   subroutine test_aquifer()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      integer(ik) :: step
      real(wp) :: aqf0, worst
      print '(a)', 'test_aquifer:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      col%theta(1:10) = 0.38_wp                            ! wet -> recharge to aquifer
      col%z_wt = params%soil_layer_z(11)                   ! start at the column bottom
      forcing%precip_ground = 1.0e-4_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_AQUIFER
      aqf0 = col%w_aquifer ; worst = 0.0_wp
      do step = 1_ik, 30_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
         worst = max(worst, abs(flux%mass_resid))
      end do
      call check_true('aquifer mass residual ~ 0', worst < 1.0e-9_wp, worst)
      call check_true('aquifer store fills from recharge', col%w_aquifer > aqf0, col%w_aquifer)
      call check_true('water table rises above bottom', col%z_wt > params%soil_layer_z(11), col%z_wt)
   end subroutine test_aquifer

   !----- Dunne saturation-excess runoff appears under a shallow water table. --------------!
   subroutine test_dunne()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      real(wp) :: f_sat_expect
      print '(a)', 'test_dunne:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      col%z_wt = -0.1_wp                                   ! shallow water table -> large f_sat
      forcing%precip_ground = 1.0e-4_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 290.0_wp ; forcing%q_air = 0.05_wp
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      opts%bottom_bc = SOIL_BC_AQUIFER
      call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
      f_sat_expect = 0.4_wp * exp(0.5_wp * 0.5_wp * (-0.1_wp)) * 1.0e-4_wp
      call check_true('Dunne runoff > 0 (shallow WT)', flux%runoff_surf > 0.0_wp, flux%runoff_surf)
      call check('Dunne runoff ~ f_sat*precip', flux%runoff_surf, f_sat_expect, 1.0e-5_wp)
      !----- BUG7: under FREE_DRAIN (no diagnosed water table) Dunne runoff must be EXACTLY ZERO,   !
      !      even for a sub-infiltration-capacity precip. The old code pinned z_wt at the geometric  !
      !      column bottom and shed f_max*exp(0.5*f_over*z_bottom)*precip (~5% here) from an          !
      !      unsaturated column every step. Only Horton ponding-overflow may run off now. -----------!
      call loam_column(SOIL_RETENTION_VG, params, col)       ! unsaturated (theta=0.30), z_wt reset
      forcing%precip_ground = 1.0e-5_wp                      ! gentle, below infiltration capacity
      opts%bottom_bc = SOIL_BC_FREE_DRAIN
      call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
      call check('free-drain Dunne runoff is exactly zero (BUG7)', flux%runoff_surf, 0.0_wp, 1.0e-15_wp)
   end subroutine test_dunne

   !=======================================================================================!
   !  Ground evaporation is an AREA-weighted TILE flux: E = snow_free_frac * rho*dq/(r_aero  !
   !  + r_soil), so it is LINEAR in the snow-free fraction and exactly 0 under full cover.    !
   !                                                                                          !
   !  This is the invariant that distinguishes it from the form it replaced, which throttled   !
   !  partial snow cover by dividing r_aero by (1-snowfac) instead. That scales only the        !
   !  AERODYNAMIC leg of the series, so it is equivalent ONLY when r_soil = 0 -- and it is       !
   !  deliberately NOT zero here: theta_1 = 0.30 sits below dsl_theta_init*theta_sat = 0.344,     !
   !  so a dry surface layer exists and r_soil > 0. Under the old form E(0.5) would exceed         !
   !  0.5*E(1) (a smaller total resistance than area-halving implies), so the linearity check      !
   !  below fails it. g_top's energy budget depends on this: it area-weights radiation and          !
   !  sensible heat by (1-snowfac) but subtracts le_soil whole, which is only self-consistent        !
   !  once le_soil is itself area-integrated. ------------------------------------------------------!
   subroutine test_snow_free_evap()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      real(wp) :: e_full, e_half, e_none
      print '(a)', 'test_snow_free_evap:'
      forcing%precip_ground = 0.0_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 305.0_wp ; forcing%q_air = 0.001_wp        ! warm ground, dry air -> strong evap
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 50.0_wp
      opts%bottom_bc = SOIL_BC_FREE_DRAIN

      call loam_column(SOIL_RETENTION_VG, params, col)
      forcing%snow_free_frac = 1.0_wp
      call column_hydrology_flux(col, forcing, params, opts, 60.0_wp, flux)
      e_full = flux%soil_evap

      call loam_column(SOIL_RETENTION_VG, params, col)              ! same initial state
      forcing%snow_free_frac = 0.5_wp
      call column_hydrology_flux(col, forcing, params, opts, 60.0_wp, flux)
      e_half = flux%soil_evap

      call loam_column(SOIL_RETENTION_VG, params, col)
      forcing%snow_free_frac = 0.0_wp
      call column_hydrology_flux(col, forcing, params, opts, 60.0_wp, flux)
      e_none = flux%soil_evap

      call check_true('bare ground evaporates (non-trivial case)', e_full > 1.0e-8_wp, e_full)
      call check('half snow cover halves E exactly (area weighting)', e_half, 0.5_wp * e_full,   &
                 1.0e-14_wp * max(e_full, 1.0e-12_wp))
      call check('full snow cover gives exactly zero soil evaporation', e_none, 0.0_wp, 1.0e-30_wp)
   end subroutine test_snow_free_evap

   !=======================================================================================!
   !  The post-solve saturation clip moves water with NO face, so the soil ENERGY column      !
   !  cannot see it -- it consumes the corrected theta but advects enthalpy only on the         !
   !  faces. The caller compensates each layer's enthalpy at that layer's own temperature,       !
   !  which requires the clip to be resolved PER LAYER. Assert the decomposition contract:        !
   !  clip_layer sums to the scalar clip_excess, is non-negative, and actually fires here.         !
   !  (A sealed bedrock column under steady rain must saturate: at saturation psi_1 -> 0 so         !
   !  the infiltration cap tends to ksat = 2.89e-6 m/s, still above this precip, so water           !
   !  keeps entering a column that cannot drain and the solver overshoots theta_sat.) --------------!
   subroutine test_clip_layer_decomposition()
      type(soil_params_t)    :: params
      type(soil_column_t)    :: col
      type(chydro_forcing_t) :: forcing
      type(soil_opts_t)      :: opts
      type(chydro_flux_t)    :: flux
      real(wp)    :: worst_gap, worst_clip, most_negative, floor_seen
      integer(ik) :: step, k
      print '(a)', 'test_clip_layer_decomposition:'
      call loam_column(SOIL_RETENTION_VG, params, col)
      col%theta(1:10) = 0.42_wp                                     ! just below theta_sat = 0.43
      forcing%precip_ground = 5.0e-4_wp ; forcing%root_uptake = 0.0_wp
      forcing%t_ground = 285.0_wp ; forcing%q_air = 0.05_wp         ! humid -> negligible evaporation
      forcing%rho_air = 1.2_wp ; forcing%r_aero = 100.0_wp
      forcing%snow_free_frac = 1.0_wp
      opts%bottom_bc = SOIL_BC_BEDROCK                              ! sealed: the water has nowhere to go
      worst_gap = 0.0_wp ; worst_clip = 0.0_wp
      most_negative = 0.0_wp ; floor_seen = 0.0_wp
      do step = 1_ik, 200_ik
         call column_hydrology_flux(col, forcing, params, opts, 600.0_wp, flux)
         worst_gap  = max(worst_gap, abs(sum(flux%clip_layer(1:10)) - flux%clip_excess))
         worst_clip = max(worst_clip, flux%clip_excess)
         floor_seen = max(floor_seen, sum(flux%floor_layer(1:10)))
         do k = 1_ik, 10_ik
            most_negative = min(most_negative, flux%clip_layer(k), flux%floor_layer(k))
         end do
      end do
      call check_true('saturation clip actually fires in this scenario', worst_clip > 0.0_wp, worst_clip)
      call check('sum(clip_layer) == clip_excess every step', worst_gap, 0.0_wp, 1.0e-18_wp)
      call check_true('clip_layer / floor_layer are non-negative', most_negative >= 0.0_wp, most_negative)
      call check('theta_res hard floor never fires on a saturating column', floor_seen, 0.0_wp, 1.0e-30_wp)
   end subroutine test_clip_layer_decomposition

end program test_column_hydrology
