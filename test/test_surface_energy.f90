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
   use meds_constants,        only : cp_air, latent_heat_vap, stefan, pi, cp_liq
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
   call test_wood_stiffness_spread()
   call test_veg_exponential()
   call test_ground_balance()
   call test_cas()
   call test_veg_energy_diagnostic_wetted()
   call test_veg_energy_diagnostic_floor()

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

   !=======================================================================================!
   ! PHASE 4 PREREQUISITE (MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md): MEASURE the wood thermal !
   ! timescale before deciding how to integrate a prognostic wood store.                      !
   !                                                                                          !
   ! MEDS_LEAF_WOOD_ENERGY_DESIGN.md sec 3 asserts "WOOD is NON-STIFF, tau ~ minutes-hours     !
   ! >> dt_fast" and phases the whole feature on it. MEASURED, that is WRONG -- and its own     !
   ! stated range contains the contradiction (minutes is not >> 1800 s).                        !
   !                                                                                            !
   ! It also refutes the replacement estimate offered in the parity plan (tau ~ 1e4*r, spanning  !
   ! orders of magnitude within a patch). That scaling assumes WAI is true bole surface area;    !
   ! in MEDS the driver sets wai = 0.20*lai (meds_fast_dynamics.f90), so WAI tracks LEAF area    !
   ! rather than stem geometry and the sapwood-mass-per-wood-area ratio barely moves across      !
   ! stand development. Both estimates are recorded here as refuted so neither propagates again. !
   !                                                                                            !
   ! What is actually true is simpler and STRONGER for the design decision: with                 !
   ! cap = dbio_w*(c_sapw + f*cp_liq) and |drdt| = WAI*(8*eps*sigma*T^3 + pi*gbh*rho*cp),        !
   ! tau_wood is 55-199 s across sapling -> mature structure -- only ~3.6x of spread, but EVERY  !
   ! cohort is stiff at dt_fast = 1800 s -- by 9x for the LARGEST cohort and 33x for the         !
   ! smallest. Wood is not "non-stiff with a stiff tail"; it is uniformly stiff.                 !
   ! veg_energy_step_implicit's own comment already names tau~35 s for a young-stand cohort.     !
   !=======================================================================================!
   subroutine test_wood_stiffness_spread()
      type(veg_thermal_params_t) :: tp
      real(wp), parameter :: C2B = 2.0_wp, WOOD_MOIST_FRAC = 1.0_wp
      real(wp), parameter :: DT_FAST = 1800.0_wp, TREF = 295.0_wp
      real(wp), parameter :: GBH = 0.03_wp, RHO = 1.2_wp
      !----- (wood_carbon [kgC/plant], lai [m2/m2], nplant [pl/m2]) : sapling -> pole -> mature -!
      real(wp), parameter :: wc(4)  = [0.5_wp,   20.0_wp,  200.0_wp, 1200.0_wp]
      real(wp), parameter :: lai(4) = [1.0_wp,    3.0_wp,    4.5_wp,    5.0_wp]
      real(wp), parameter :: npl(4) = [2.0_wp,    0.3_wp,   0.06_wp,   0.015_wp]
      real(wp) :: dbio_w, wai, cap, drdt, tau, tau_min, tau_max
      integer(ik) :: i
      print '(a)', 'test_wood_stiffness_spread:'
      tau_min = huge(1.0_wp) ; tau_max = 0.0_wp
      do i = 1_ik, 4_ik
         wai    = 0.20_wp * lai(i)
         dbio_w = 0.10_wp * wc(i) * npl(i) * C2B                       ! [kg dry/m2 ground]
         cap    = max(dbio_w * tp%c_sapw, tp%veg_hcap_min) + dbio_w * WOOD_MOIST_FRAC * cp_liq
         drdt   = 8.0_wp * tp%leaf_emiss * stefan * TREF**3 * wai                                &
                + pi * wai * GBH * RHO * cp_air
         tau    = cap / drdt
         print '(a,i0,a,f9.1,a,f8.1,a)', '   cohort ', i, ': cap = ', cap,                       &
               ' J/m2/K   tau_wood = ', tau, ' s'
         tau_min = min(tau_min, tau) ; tau_max = max(tau_max, tau)
      end do
      !----- The fact Phase 4 turns on: EVERY cohort is stiff, not just a small-wood tail. -------!
      call check_true('EVERY cohort is stiff at dt_fast = 1800 s (design doc claims non-stiff)',      &
           tau_max < DT_FAST, tau_max)
      call check_true('the spread is modest (~3.6x) -- WAI tracks LAI, not bole surface area',        &
           tau_max / tau_min < 10.0_wp, tau_max / tau_min)
      !----- ...which is why prognostic wood is operator-split behind the L-stable kernel on ALL  !
      !      three schemes rather than placed in RK45's explicit tableau: an explicit march would  !
      !      be step-limited by tau for EVERY cohort, not just the worst one. --------------------!
      call check_true('an explicit tableau would be step-limited for the WHOLE stand (>9x)',          &
           DT_FAST / tau_max > 9.0_wp, DT_FAST / tau_max)
   end subroutine test_wood_stiffness_spread

   !=======================================================================================!
   ! EXACT EXPONENTIAL RELAXATION in veg_energy_diagnostic.                                  !
   !                                                                                          !
   ! Under the Category-0 freeze the tissue equation is linear with tau = cap/denom, so the     !
   ! step has a closed form and the kernel uses it: endpoint weight exp(-x), flux weight        !
   ! (1-exp(-x))/x, x = denom/a_store. These four checks pin the properties that makes it        !
   ! usable -- the two limits, the conservation identity, and (the one that matters) that        !
   ! marching it is EXACT, not merely convergent.                                                !
   !=======================================================================================!
   subroutine test_veg_exponential()
      type(veg_thermal_params_t) :: tp
      real(wp) :: dtt, ts, tr, dh, drn, dt_diag_ref, ts_ref
      real(wp) :: t_cas, t0, sw, lw, h, les, lws, ler, cap, dt, denom, tau, a_st
      real(wp) :: t_exact, t_march, worst
      integer(ik) :: k, nstep
      print '(a)', 'test_veg_exponential:'
      t_cas = 295.0_wp ; sw = 300.0_wp ; lw = -40.0_wp
      h = 120.0_wp ; les = 30.0_wp ; lws = 5.0_wp ; ler = 25.0_wp
      denom = h + les + lws
      dt    = 1800.0_wp

      !----- (a) a_store -> 0 is the DIAGNOSTIC limit, exactly. This is the property that lets the   !
      !          selector be deleted: diagnostic is a limit of one formula, not a second branch. ----!
      call veg_energy_diagnostic(sw, lw, h, les, lws, ler, t_cas, t_cas, 0.0_wp, t_cas,             &
                                 dtt, ts, tr, dh, drn)
      dt_diag_ref = (sw + lw - ler) / denom
      call check('a_store = 0 gives the exact diagnostic offset', dtt, dt_diag_ref, 1.0e-12_wp)
      ts_ref = t_cas + dt_diag_ref

      !----- (b) a_store -> huge FREEZES the store at its entry temperature. -----------------------!
      t0 = 288.0_wp
      call veg_energy_diagnostic(sw, lw, h, les, lws, ler, t_cas, t_cas, 1.0e12_wp, t0,             &
                                 dtt, ts, tr, dh, drn)
      call check('a_store -> infinity holds the store at t_store0', ts, t0, 1.0e-6_wp)

      !----- (c) THE CONSERVATION IDENTITY. The endpoint and step-average weights are different      !
      !          numbers, and pairing them correctly is exactly what makes the balance close:        !
      !            a_store*(dt_end - dt_prev) + denom*dt_avg == numer.                               !
      !          dt_avg is not returned, so recover it from dh (which carries denom's sensible       !
      !          share plus g_slave; here denom_true > floor so g_slave = 0 and dh = h*dt_avg). -----!
      cap  = 4.0e4_wp                       ! a big cohort: tau = cap/denom is order dt
      a_st = cap / dt
      tau  = cap / denom
      call veg_energy_diagnostic(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t0,                  &
                                 dtt, ts, tr, dh, drn)
      call check('energy balance closes with the endpoint/average pair',                            &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-9_wp)

      !----- (d) THE POINT: marching the kernel reproduces the analytic relaxation EXACTLY, not      !
      !          approximately. A backward-Euler storage term would drift here; this does not.       !
      !          T(t) = T_eq + (T0 - T_eq)*exp(-t/tau).                                              !
      nstep = 20_ik ; t_march = t0 ; worst = 0.0_wp
      do k = 1_ik, nstep
         call veg_energy_diagnostic(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t_march,          &
                                    dtt, ts, tr, dh, drn)
         t_march = ts
         t_exact = ts_ref + (t0 - ts_ref) * exp(-real(k, wp) * dt / tau)
         worst   = max(worst, abs(t_march - t_exact))
      end do
      call check_true('marching 20 steps is EXACT against the analytic relaxation', worst < 1.0e-10_wp, worst)

      !----- (e) The small-x series branch (a very large capacity) must not lose precision where     !
      !          1 - exp(-x) would cancel. tau = 1e6 s against dt = 1800 s puts x ~ 2.7e-4. ---------!
      a_st = 1.0e6_wp * denom / dt
      call veg_energy_diagnostic(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t0,                  &
                                 dtt, ts, tr, dh, drn)
      call check('large-capacity limit stays accurate (series branch)',                             &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-6_wp)
   end subroutine test_veg_exponential

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


   subroutine test_veg_energy_diagnostic_floor()
      !----- A DIAGNOSTIC tissue is a ROUTER, not a reservoir: within a step everything that enters   !
      !      must leave, so                                                                           !
      !                                                                                          !
      !         drnet + q_extra  ==  dh + Lv*(transp + film_evap) + a_store*(t_store - t_store0)       !
      !                                                                                          !
      !      to round-off, for ANY inputs.  meds_fast_split / surface_derivs both rely on this: the     !
      !      whole-column ledger books coh_rnet on the input side and coh_h/coh_qw as they reach the    !
      !      CAS, and NOTHING re-checks the tissue in between -- there is no per-kernel budget for a     !
      !      store that holds nothing.                                                                  !
      !                                                                                          !
      !      test_veg_energy_diagnostic_wetted already asserts this identity, but at h_coeff = 12 and    !
      !      lw_slope = 6 -- a coupling ~30x the kernel's veg_coupling_floor, i.e. ONLY on the side of    !
      !      the clamp where the kernel is right.  A check that samples one side of a clamp is green      !
      !      only on that side.  This sweeps the area index ACROSS the floor, where dt_temp comes from     !
      !      the FLOORED denominator while every flux keeps its own TRUE conductance.                     !
      !                                                                                          !
      !      Sweep design, each axis for a reason:                                                        !
      !        * area index log-spaced 1e-4..10 -- straddles both crossings (leaf ~0.03, wood ~0.02);      !
      !        * every coupling term scales WITH the area index, as a real cohort's does, so a failure     !
      !          is attributable to the floor rather than to an unphysical corner;                         !
      !        * q_extra deliberately does NOT scale with area -- sapflow enthalpy scales with sap flux,   !
      !          not leaf area, and it is the term that keeps the numerator alive as the denominator       !
      !          collapses.  The wetted test never passes q_extra at all;                                  !
      !        * leaf AND wood: wood has no latent pathway, so its denominator is h_coeff + lw_slope        !
          !      alone and it crosses the floor at a LARGER area index than the leaf on the same cohort;   !
      !        * a_store 0 and veg_hcap_min/dt_fast: the PROGNOSTIC-leaf path calls this same kernel, and   !
      !          for the small cohorts that trip the floor cap is pinned at veg_hcap_min, contributing      !
      !          only ~0.011 W/m2/K -- so the storage term must not be assumed to rescue the balance.       !
      !                                                                                          !
      !      Asserted on CONSERVATION only.  Whether the floor is active, and what dt_temp it yields, are   !
      !      implementation choices a fix is free to change; energy closure is not. -------------------------!
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: ai, qx, ast, resid, les, ler, lesw, lerw, hc, lws, worst
      integer(ik) :: ia, iq, il, is
      real(wp), parameter :: T_CAS = 298.0_wp, GBH = 0.03_wp, RHO = 1.2_wp, EMISS = 0.97_wp
      logical :: is_leaf
      print '(a)', 'test_veg_energy_diagnostic_floor:'
      worst = 0.0_wp
      do iq = 0_ik, 1_ik
         qx = 0.0_wp ; if (iq == 1_ik) qx = 25.0_wp            ! sapflow enthalpy: area-INDEPENDENT
         do is = 0_ik, 1_ik
            ast = 0.0_wp ; if (is == 1_ik) ast = 20.0_wp / 1800.0_wp   ! veg_hcap_min / dt_fast
            do il = 1_ik, 2_ik
               is_leaf = (il == 1_ik)
               do ia = 1_ik, 11_ik
                  ai = 1.0e-4_wp * (10.0_wp ** (0.5_wp * real(ia - 1_ik, wp)))
                  hc  = merge(2.0_wp, pi, is_leaf) * ai * GBH * RHO * cp_air
                  lws = 4.0_wp * EMISS * stefan * T_CAS**3 * ai
                  les  = merge( 8.0_wp * ai, 0.0_wp, is_leaf)    ! stomatal (wood does not transpire)
                  ler  = merge(15.0_wp * ai, 0.0_wp, is_leaf)
                  lesw = 20.0_wp * ai                            ! film pathway: both tissues have one
                  lerw = 40.0_wp * ai
                  call veg_energy_diagnostic(400.0_wp*ai, -50.0_wp*ai, hc, les, lws, ler,           &
                                             T_CAS, T_CAS, ast, T_CAS - 0.5_wp,                     &
                                             dt_temp, t_store, transp, dh, drnet,                   &
                                             f_wet=0.4_wp, le_slope_wet=lesw, le_ref_wet=lerw,      &
                                             film_evap=film_evap, q_extra=qx)
                  resid = drnet + qx - dh - latent_heat_vap*(transp + film_evap)                    &
                          - ast * (t_store - (T_CAS - 0.5_wp))
                  worst = max(worst, abs(resid) / max(1.0_wp, abs(drnet)))
               end do
            end do
         end do
      end do
      call check_true('diagnostic tissue conserves across the coupling floor (88 cases)',           &
                      worst < 1.0e-9_wp, worst)
   end subroutine test_veg_energy_diagnostic_floor

end program test_surface_energy
