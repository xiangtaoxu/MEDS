!==========================================================================================!
! test_surface_energy -- unit tests for the leaf/wood, ground, and canopy-air-space kernels.  !
!   1. LEAF/WOOD energy CONSERVATION to ~round-off over a march.                                !
!   2. LEAF RELAXATION: with no radiation and no evaporation, a leaf relaxes to can_temp.        !
!   3. GROUND balance: with can_temp = t_ground and saturated air, G_top = Rn (H = LE = 0).       !
!   4. CAS enthalpy update CONSERVES and moves can_temp toward the warmer atmosphere.             !
!==========================================================================================!
program test_surface_energy
   use meds_kinds,            only : wp, ik
   use meds_biophysics_types, only : leaf_energy_env_t, leaf_energy_flux_t, veg_thermal_params_t
   use meds_therm_lib,           only : temp_to_uext, sat_specific_humidity, cas_enthalpy_of_temp
   use meds_vegetation_biophysics, only : veg_energy_step_implicit
   use meds_ground_biophysics, only : ground_surface_balance
   use meds_cas_biophysics,   only : canopy_air_update
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_leaf_conserve()
   call test_leaf_relax()
   call test_wood_conserve()
   call test_ground_balance()
   call test_cas()

   if (nfail == 0_ik) then
      print '(a)', 'test_surface_energy: ALL PASSED'
   else
      print '(a,i0,a)', 'test_surface_energy: ', nfail, ' FAILED'
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

   subroutine make_leaf_env(env)
      type(leaf_energy_env_t), intent(out) :: env
      env%abs_sw = 400.0_wp ; env%abs_lw = -50.0_wp
      env%can_temp = 298.0_wp ; env%can_shv = 0.012_wp
      env%gbh = 0.03_wp ; env%gbw = 0.03_wp ; env%gsw = 0.005_wp ; env%fs_open = 1.0_wp
      env%area_index = 3.0_wp ; env%leaf_water = 0.05_wp ; env%wmass = 0.30_wp
      env%dry_hcap = 1000.0_wp ; env%rho_air = 1.2_wp ; env%press = 101325.0_wp
   end subroutine make_leaf_env

   subroutine test_leaf_conserve()
      type(leaf_energy_env_t)  :: env
      type(veg_thermal_params_t) :: tp
      type(leaf_energy_flux_t) :: flux
      real(wp) :: se, worst
      integer(ik) :: step
      print '(a)', 'test_leaf_conserve:'
      call make_leaf_env(env)
      se = temp_to_uext(env%dry_hcap, env%wmass, 300.0_wp, 1.0_wp)
      worst = 0.0_wp
      do step = 1_ik, 50_ik
         call veg_energy_step_implicit(se, env, tp, 60.0_wp, .true., flux)
         worst = max(worst, abs(flux%energy_resid))
      end do
      call check_true('leaf energy residual ~ 0', worst < 1.0e-6_wp, worst)
      call check_true('leaf temperature physical (250-350 K)', flux%temp > 250.0_wp .and.      &
                      flux%temp < 350.0_wp, flux%temp)
   end subroutine test_leaf_conserve

   subroutine test_leaf_relax()
      type(leaf_energy_env_t)  :: env
      type(veg_thermal_params_t) :: tp
      type(leaf_energy_flux_t) :: flux
      real(wp) :: se
      integer(ik) :: step
      print '(a)', 'test_leaf_relax:'
      call make_leaf_env(env)
      env%abs_sw = 0.0_wp ; env%abs_lw = 0.0_wp                   ! no radiative source
      env%gbw = 0.0_wp ; env%gsw = 0.0_wp                         ! no evaporation -> pure sensible
      env%can_temp = 295.0_wp
      se = temp_to_uext(env%dry_hcap, env%wmass, 305.0_wp, 1.0_wp)  ! start hot
      do step = 1_ik, 200_ik
         call veg_energy_step_implicit(se, env, tp, 60.0_wp, .true., flux)
      end do
      call check('leaf relaxes to can_temp', flux%temp, 295.0_wp, 0.2_wp)
   end subroutine test_leaf_relax

   subroutine test_wood_conserve()
      type(leaf_energy_env_t)  :: env
      type(veg_thermal_params_t) :: tp
      type(leaf_energy_flux_t) :: flux
      real(wp) :: se, worst
      integer(ik) :: step
      print '(a)', 'test_wood_conserve:'
      call make_leaf_env(env)
      env%area_index = 1.0_wp ; env%dry_hcap = 3000.0_wp ; env%wmass = 0.6_wp    ! wood-like
      se = temp_to_uext(env%dry_hcap, env%wmass, 299.0_wp, 1.0_wp)
      worst = 0.0_wp
      do step = 1_ik, 50_ik
         call veg_energy_step_implicit(se, env, tp, 60.0_wp, .false., flux)             ! is_leaf = .false.
         worst = max(worst, abs(flux%energy_resid))
      end do
      call check_true('wood energy residual ~ 0', worst < 1.0e-6_wp, worst)
      call check_true('wood has no transpiration', abs(flux%q_transp) < 1.0e-30_wp, flux%q_transp)
   end subroutine test_wood_conserve

   subroutine test_ground_balance()
      type(leaf_energy_env_t) :: env
      real(wp) :: t_ground, g_top, h_ground, le_ground
      print '(a)', 'test_ground_balance:'
      call make_leaf_env(env)
      env%abs_sw = 300.0_wp ; env%abs_lw = -40.0_wp
      t_ground = 300.0_wp
      env%can_temp = t_ground                                     ! no sensible gradient
      env%can_shv  = sat_specific_humidity(t_ground, env%press)   ! saturated air -> no evap gradient
      call ground_surface_balance(t_ground, env, g_top, h_ground, le_ground)
      call check_true('H_ground = 0 (no gradient)', abs(h_ground) < 1.0e-9_wp, h_ground)
      call check_true('LE_ground = 0 (saturated)',  abs(le_ground) < 1.0e-6_wp, le_ground)
      call check('G_top = Rn', g_top, env%abs_sw + env%abs_lw, 1.0e-6_wp)
   end subroutine test_ground_balance

   subroutine test_cas()
      real(wp) :: enth, shv, temp, resid, enth_atm, worst
      integer(ik) :: step
      print '(a)', 'test_cas:'
      shv  = 0.012_wp
      enth = cas_enthalpy_of_temp(295.0_wp, shv)                  ! CAS starts at 295 K
      temp = 295.0_wp
      enth_atm = cas_enthalpy_of_temp(300.0_wp, shv)              ! warmer atmosphere
      worst = 0.0_wp
      do step = 1_ik, 30_ik
         call canopy_air_update(enth, shv, temp, 20.0_wp,                                      &
              50.0_wp, 20.0_wp, 0.0_wp, 0.0_wp,          &  ! cohort sensible + vapour-enthalpy fluxes
              10.0_wp, 5.0_wp, 0.0_wp, 0.0_wp,           &  ! ground fluxes; dew = 0
              0.3_wp, 1.0_wp, enth_atm, 0.0_wp, 1.2_wp, 60.0_wp, resid)   ! ustar, temp1(c3), enth_atm, w_flux_ac, ...
         worst = max(worst, abs(resid))
      end do
      call check_true('CAS energy residual ~ 0 (30 steps)', worst < 1.0e-6_wp, worst)
      call check_true('CAS warms toward atmosphere + surfaces', temp > 295.0_wp, temp)
   end subroutine test_cas

end program test_surface_energy
