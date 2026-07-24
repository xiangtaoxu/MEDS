!==========================================================================================!
! test_surface_energy -- unit tests for the leaf/wood, ground, and canopy-air-space kernels.  !
!   1. LEAF/WOOD energy CONSERVATION to ~round-off over a march.                                !
!   2. LEAF RELAXATION: with no radiation and no evaporation, a leaf relaxes to can_temp.        !
!   3. GROUND fluxes: with t_cas = t_ground and no soil evaporation, H = LE = 0.                   !
!   4. CAS enthalpy update CONSERVES and moves can_temp toward the warmer atmosphere.             !
!   5. veg_energy_diagnostic's sec 3.4 (P1) wetted-canopy extension (f_wet/le_slope_wet/            !
!      le_ref_wet/film_evap) -- no direct unit test existed anywhere for this kernel before.        !
!==========================================================================================!
program test_surface_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : cp_air, latent_heat_vap
   use meds_biophysics_types, only : leaf_energy_env_t, leaf_energy_flux_t, veg_thermal_params_t
   use meds_therm_lib,           only : temp_to_uext, sat_specific_humidity, cas_enthalpy_of_temp, &
                                        cas_temp_of_enthalpy
   use meds_vegetation_biophysics, only : veg_energy_step_implicit, veg_energy_diagnostic
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_leaf_conserve()
   call test_leaf_relax()
   call test_wood_conserve()
   call test_ground_balance()
   call test_cas()
   call test_veg_energy_diagnostic_wetted()

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
      real(wp) :: t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil, g_top, rn
      print '(a)', 'test_ground_balance:'
      ggnet = 0.02_wp ; rho = 1.2_wp ; t_ground = 300.0_wp ; rn = 300.0_wp - 40.0_wp
      !----- no sensible gradient (t_cas = t_ground) and no soil evaporation -> both vanish. -!
      t_cas = t_ground ; soil_evap = 0.0_wp
      call ground_surface_fluxes(t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil)
      call check_true('H_bare = 0 (no gradient)', abs(h_bare) < 1.0e-9_wp, h_bare)
      call check_true('LE_soil = 0 (no evap)',    abs(le_soil) < 1.0e-30_wp, le_soil)
      !----- warm ground over cooler CAS -> positive sensible; caller assembles G_top = Rn-H-LE. -!
      t_cas = 298.0_wp
      call ground_surface_fluxes(t_ground, t_cas, ggnet, rho, soil_evap, h_bare, le_soil)
      call check('H_bare = ggnet*rho*cp*dT', h_bare, ggnet*rho*cp_air*(t_ground - t_cas), 1.0e-9_wp)
      g_top = rn - h_bare - le_soil
      call check('G_top = Rn - H - LE', g_top, rn - ggnet*rho*cp_air*(t_ground - t_cas), 1.0e-6_wp)
   end subroutine test_ground_balance

   subroutine test_cas()
      !----- Migrated onto the production cas_column_step_implicit (enthalpy twin; the retired      !
      !      canopy_air_update advanced the same box math). Vapour + CO2 twins held inert (0        !
      !      sources / 0 vapour conductance; a nonzero molar capacity avoids a 0/0 in the unused    !
      !      CO2 branch). The closed-budget residual is assembled here from the production result.   !
      type(cas_column_t) :: column
      type(cas_source_t) :: source
      real(wp) :: enth, shv, temp, resid, enth_atm, worst, enth_new, shv_new, co2_new
      real(wp), parameter :: wcap = 1.2_wp * 20.0_wp     ! rho_air * can_depth  [kg/m2]
      real(wp), parameter :: gatm = 1.2_wp * 0.3_wp * 1.0_wp  ! rho_air * ustar * temp1(c3)  [kg/m2/s]
      integer(ik) :: step
      print '(a)', 'test_cas:'
      shv  = 0.012_wp
      enth = cas_enthalpy_of_temp(295.0_wp, shv)                  ! CAS starts at 295 K
      temp = 295.0_wp
      enth_atm = cas_enthalpy_of_temp(300.0_wp, shv)              ! warmer atmosphere
      column%air_mass_capacity        = wcap
      column%air_molar_capacity       = 1.0_wp                    ! inert CO2 twin (avoid 0/0)
      column%atm_conductance_enthalpy = gatm
      column%atm_conductance_vapor    = 0.0_wp
      column%atm_conductance_co2      = 0.0_wp
      column%atm_enthalpy             = enth_atm
      column%atm_specific_humidity    = 0.0_wp
      column%atm_co2                  = 0.0_wp
      source%surface_enthalpy_source  = 50.0_wp + 20.0_wp + 10.0_wp + 5.0_wp  ! cohort + ground sensible + vapour-enthalpy
      source%surface_vapor_source     = 0.0_wp
      source%biotic_co2_source        = 0.0_wp
      worst = 0.0_wp
      do step = 1_ik, 30_ik
         call cas_column_step_implicit(enth, shv, 0.0_wp, source, column, 60.0_wp, enth_new, shv_new, co2_new)
         resid = column%air_mass_capacity * (enth_new - enth)                                   &
                 - 60.0_wp * (source%surface_enthalpy_source                                    &
                              + column%atm_conductance_enthalpy * (column%atm_enthalpy - enth_new))
         worst = max(worst, abs(resid))
         enth = enth_new ; shv = shv_new
         temp = cas_temp_of_enthalpy(enth, shv)
      end do
      call check_true('CAS energy residual ~ 0 (30 steps)', worst < 1.0e-6_wp, worst)
      call check_true('CAS warms toward atmosphere + surfaces', temp > 295.0_wp, temp)
   end subroutine test_cas

   subroutine test_veg_energy_diagnostic_wetted()
      !----- Direct unit test of the sec 3.4 (P1) wet/dry extension -- veg_energy_diagnostic had NO   !
      !      direct test anywhere before this (only ever exercised via the full column integration     !
      !      tests). f_wet splits the single latent pathway into a (1-f_wet) DRY (stomatal) share and   !
      !      an f_wet WET (boundary-layer film) share; le_slope_dry/le_ref_dry model the former,          !
      !      le_slope_w/le_ref_w the latter -- two DIFFERENT conductances, so the three properties        !
      !      below are non-trivial (a bug that swapped or mixed them would still "look" plausible).       !
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_temp2, t_store2, transp2, dh2, drnet2
      real(wp), parameter :: abs_sw = 400.0_wp, abs_lw = -50.0_wp, h_coeff = 12.0_wp
      real(wp), parameter :: lw_slope = 6.0_wp, t_cas = 298.0_wp, t_emit = 298.0_wp
      real(wp), parameter :: le_slope_dry = 8.0_wp,  le_ref_dry = 15.0_wp   ! stomatal (dry) pathway
      real(wp), parameter :: le_slope_w   = 20.0_wp, le_ref_w   = 40.0_wp   ! boundary-layer (wet) pathway
      print '(a)', 'test_veg_energy_diagnostic_wetted:'

      !----- 1. ABSENT f_wet must be bit-identical to explicitly passing f_wet=0 -- the contract every   !
      !      existing caller (ARK's surface_derivs, the split's own wood branch pre-P1) relies on.  -----!
      call veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_emit, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet)
      call veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_emit, 0.0_wp, t_cas, dt_temp2, t_store2, transp2, dh2, drnet2, &
                                 0.0_wp, le_slope_w, le_ref_w, film_evap)
      call check('f_wet=0 explicit == f_wet absent (dt_temp)', dt_temp2, dt_temp, 1.0e-14_wp)
      call check('f_wet=0 explicit == f_wet absent (transp)',  transp2,  transp,  1.0e-14_wp)
      call check_true('f_wet=0 -> film_evap = 0 exactly', film_evap == 0.0_wp, film_evap)

      !----- 2. FULLY WET (f_wet=1): the dry (stomatal) pathway must vanish -- a saturated film          !
      !      transpires ~0 and evaporates its whole latent capacity through the wet pathway instead.      !
      !      SWAP-EQUIVALENCE: this must be algebraically IDENTICAL to calling the kernel with the wet     !
      !      conductance as the PRIMARY (le_slope,le_ref) argument and no wet extension at all -- not       !
      !      just approximately equal (both reduce to the same 2-unknown linear solve). -------------------!
      call veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_emit, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet,      &
                                 1.0_wp, le_slope_w, le_ref_w, film_evap)
      call veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope_w, lw_slope, le_ref_w,               &
                                 t_cas, t_emit, 0.0_wp, t_cas, dt_temp2, t_store2, transp2, dh2, drnet2)
      call check_true('f_wet=1: dry pathway (transp) vanishes', abs(transp) < 1.0e-14_wp, transp)
      call check('f_wet=1: film_evap == swap-equivalent transp', film_evap, transp2, 1.0e-12_wp)
      call check('f_wet=1: dt_temp == swap-equivalent dt_temp',  dt_temp,   dt_temp2, 1.0e-12_wp)

      !----- 3. PARTIAL wetting (f_wet=0.4): the energy balance must still close -- absorbed radiation     !
      !      net of LW emission (drnet) splits EXACTLY into sensible + BOTH latent pathways, each at the    !
      !      CONSTANT latent_heat_vap (not enthalpy_vapor(T)) -- the identity meds_fast_split.f90 relies     !
      !      on to make the CAS's temperature-dependent vapour-enthalpy credit balance against a store's     !
      !      liquid-enthalpy debit (coh_qsoil for transpiration; the analogous surface-water accounting       !
      !      for film_evap). ------------------------------------------------------------------------------!
      call veg_energy_diagnostic(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_emit, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet,      &
                                 0.4_wp, le_slope_w, le_ref_w, film_evap)
      call check('partial f_wet: energy balance closes (drnet = dh + latent)', drnet,                  &
                 dh + (transp + film_evap) * latent_heat_vap, 1.0e-9_wp)
      call check_true('partial f_wet: film_evap > 0 (wet pathway active)', film_evap > 0.0_wp, film_evap)
      call check_true('partial f_wet: transp > 0 (dry pathway still active)', transp > 0.0_wp, transp)
   end subroutine test_veg_energy_diagnostic_wetted

end program test_surface_energy
