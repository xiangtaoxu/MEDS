! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_surface_energy -- unit tests for the leaf/wood, ground, and canopy-air-space kernels.  !
!   1. WOOD THERMAL TIMESCALE from the real allometry: the wood store carries memory at dt_fast.   !
!   2. EXACT RELAXATION: veg_energy_balance's closed-form tissue step, its two limits and its      !
!      energy identity (docs/science/vegetation_energy_dynamics.md eqs. (3)-(7)).                  !
!   3. GROUND fluxes: with t_cas = t_ground and no soil evaporation, H = LE = 0.                   !
!   4. CAS enthalpy update CONSERVES and moves can_temp toward the warmer atmosphere.              !
!   5. veg_energy_balance's wetted-canopy extension (f_wet/le_slope_wet/le_ref_wet/film_evap).     !
!   6. A THIN tissue: energy closes at every area index, halving the area halves every flux and    !
!      keeps the temperature, and a tissue with no area takes the temperature of what it grows     !
!      from.                                                                                       !
!==========================================================================================!
program test_surface_energy
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : cp_air, latent_heat_vap, stefan, pi, cp_liq, TISSUE_C_WOOD_DRY
   use meds_allometry,        only : dbh_to_wai, sapwood_fraction
   use meds_plant_types, only : veg_thermal_params_t
   use meds_therm_lib,           only : temp_to_internal_energy, sat_specific_humidity, cas_enthalpy_of_temp, &
                                        cas_temp_of_enthalpy
   use meds_plant_biophysics, only : veg_energy_balance
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_step_implicit
   implicit none

   call test_wood_stiffness_spread()
   call test_veg_exponential()
   call test_ground_balance()
   call test_cas()
   call test_veg_energy_diagnostic_wetted()
   call test_veg_energy_thin_tissue()

   call test_report('test_surface_energy')

contains







   !=======================================================================================!
   ! WOOD THERMAL TIMESCALE, measured from the REAL allometry.                               !
   !                                                                                          !
   ! wai comes from the ED2 WAI allometry (a STEM-SIZE power law) and bsap from the sapwood   !
   ! fraction of basal area, the masses the model sizes the store from                        !
   ! (docs/science/vegetation_energy_dynamics.md eq. (8)). A test that hardwired fixed shares !
   ! of leaf area or wood carbon would only re-measure those shares.                          !
   !                                                                                          !
   ! The allometry matters qualitatively, because tau_wood ~ cap/(coupling) ~ (wood mass) /   !
   ! (wood AREA), and tying wood area to LEAF area makes the ratio nearly size-independent.   !
   ! With a stem-size WAI the ratio grows with stem size, as a length scale should.           !
   !=======================================================================================!
   subroutine test_wood_stiffness_spread()
      real(wp), parameter :: C2B = 2.0_wp, DT_FAST = 900.0_wp, TREF = 295.0_wp
      !----- The default wood store's saturated water, [hydraulics].wood_water_sat [kg/kgC]: a      !
      !      saturated store, the most water the wood's heat capacity can count. ------------------!
      real(wp), parameter :: WOOD_WATER_SAT = 1.0_wp
      real(wp), parameter :: GBH = 0.03_wp, RHO = 1.2_wp
      real(wp), parameter :: WOOD_EMISS = 0.90_wp       !< the default PFT wood_emissivity
      !----- ED2 ed_params.f90: b1WAI = 0.0192*0.5 (note the *0.5), b2WAI = 2.0947 for tropical  !
      !      broadleaf; b1SA = 1.582, b2SA = 1.764 (iallom default). ----------------------------!
      real(wp), parameter :: WAI_B1 = 0.0096_wp, WAI_B2 = 2.0947_wp
      real(wp), parameter :: SA_B1  = 1.582_wp,  SA_B2  = 1.764_wp
      !----- (dbh [cm], wood_carbon [kgC/plant], nplant [pl/m2]) : sapling -> pole -> mature ---!
      real(wp), parameter :: dbh(4) = [  2.0_wp,   12.0_wp,  35.0_wp,  70.0_wp]
      real(wp), parameter :: wc(4)  = [  0.5_wp,   20.0_wp, 200.0_wp,1200.0_wp]
      real(wp), parameter :: npl(4) = [  2.0_wp,    0.3_wp,   0.06_wp,  0.015_wp]
      real(wp) :: wai, f_sap, dbio_w, dbio_sap, cap, drdt, tau, tau_min, tau_max, w_end
      real(wp) :: f_sap_small, f_sap_big, sap_big
      integer(ik) :: i
      print '(a)', 'test_wood_stiffness_spread:'
      tau_min = huge(1.0_wp) ; tau_max = 0.0_wp
      f_sap_small = 0.0_wp ; f_sap_big = 0.0_wp ; sap_big = 0.0_wp
      do i = 1_ik, 4_ik
         wai    = dbh_to_wai(dbh(i), npl(i), WAI_B1, WAI_B2)
         f_sap  = sapwood_fraction(dbh(i), SA_B1, SA_B2)
         !----- The two halves of the capacity take DIFFERENT masses (design decision):          !
         !        DRY tissue   -> ALL the wood. Heartwood is dead structure but still stores      !
         !                        sensible heat, and a bole-defined sapwood fraction under-counts !
         !                        branch wood, which is thermally active throughout.              !
         !        INTERNAL WATER -> the water network's wood store, here a saturated sapwood ring !
         !                        with no fine roots. Heartwood holds none, which is what makes   !
         !                        it heartwood.                                                   !
         !      ONE enthalpy pool at ONE temperature either way -- a bole/branch partition would  !
         !      mean two wood temperatures and is deliberately not taken. --------------------!
         dbio_w   = wc(i) * npl(i) * C2B                                ! [kg dry/m2] all wood
         dbio_sap = f_sap * dbio_w                                      ! [kg dry/m2] sapwood ring
         cap    = dbio_w * TISSUE_C_WOOD_DRY                                                    &
                + WOOD_WATER_SAT * dbio_sap / C2B * cp_liq
         drdt   = 8.0_wp * WOOD_EMISS * stefan * TREF**3 * wai + pi * wai * GBH * RHO * cp_air
         tau    = cap / drdt
         !----- the EXACT endpoint memory the kernel uses (eq. (4)), not backward Euler's 1/(1+x) !
         w_end  = exp(-DT_FAST / tau)
         print '(a,i0,a,f6.2,a,f6.3,a,f6.3,a,f10.1,a,f8.1,a,f7.4)', '   cohort ', i, ': dbh = ', &
               dbh(i), '  wai = ', wai, '  f_sap = ', f_sap, '  cap = ', cap,                    &
               ' J/m2/K  tau = ', tau, ' s  w_end = ', w_end
         tau_min = min(tau_min, tau) ; tau_max = max(tau_max, tau)
         if (i == 1_ik) f_sap_small = f_sap
         if (i == 4_ik) then
            f_sap_big = f_sap ; sap_big = dbio_sap
         end if
      end do
      !----- WHAT THIS MEASURES. tau_wood = 600-760 s and w_end = 0.22-0.31 at dt_fast = 900 s, with a  !
      !      saturated store, so wood carries GENUINE memory at production cadence -- the whole reason  !
      !      the store exists. Most of the capacity is the sapwood and its water: dry wood holds 1217   !
      !      J/kg/K against water's 4186, so the heartwood adds little. ED2's b1WAI (written            !
      !      0.0192*0.5) halves WAI, hence the coupling, which lengthens tau. --------------------------!
      call check_true('tau_wood is several hundred seconds', tau_min > 500.0_wp .and. tau_max < 1000.0_wp, tau_max)
      call check_true('wood carries GENUINE memory at dt_fast = 900 s (w_end > 0.2)',               &
           exp(-DT_FAST / tau_min) > 0.2_wp, exp(-DT_FAST / tau_min))
      !----- The DRY mass must be the FULL wood mass: 36 kg/m2 = 1200 kgC/plant * 0.015 pl/m2 * C2B, !
      !      which is what an independent mass audit of that cohort gives. Re-pointing the dry half  !
      !      at bsap would drop it to ~27 and fire this. --------------------------------------!
      call check('largest cohort dry mass is ALL the wood (audit: 36 kg dry/m2)',                   &
           wc(4) * npl(4) * C2B, 36.0_wp, 1.0e-9_wp)
      !----- ...but its WATER is only the sapwood ring, which is strictly less. A small stem is      !
      !      entirely sapwood (the fraction saturates at 1), a 70 cm bole is not -- that contrast is !
      !      the physical content of splitting the two masses. --------------------------------!
      call check('a 2 cm stem is ALL sapwood', f_sap_small, 1.0_wp, 1.0e-12_wp)
      call check_true('a 70 cm bole is NOT all sapwood', f_sap_big > 0.5_wp .and.                   &
           f_sap_big < 0.95_wp, f_sap_big)
      call check_true('sapwood water mass < full wood mass', sap_big < 36.0_wp, sap_big)
      !----- Modest spread: cap and WAI both grow with size, so they largely cancel. ----------!
      call check_true('spread is modest (cap and WAI grow together)', tau_max / tau_min < 2.0_wp,   &
           tau_max / tau_min)
   end subroutine test_wood_stiffness_spread

   !=======================================================================================!
   ! EXACT EXPONENTIAL RELAXATION in veg_energy_balance.                                  !
   !                                                                                          !
   ! Under the Category-0 freeze the tissue equation is linear with tau = cap/denom, so the     !
   ! step has a closed form and the kernel uses it: endpoint weight exp(-x), flux weight        !
   ! (1-exp(-x))/x, x = denom/store_hcap_per_dt (docs/science/vegetation_energy_dynamics.md     !
   ! eqs. (3)-(7)). These checks pin the properties that make it usable -- the two limits, the  !
   ! conservation identity, the small-x branch, and (the one that matters) that marching it is  !
   ! EXACT, not merely convergent.                                                              !
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

      !----- (a) store_hcap_per_dt -> 0 is the DIAGNOSTIC limit, exactly. This is why there is no     !
      !          diagnostic/prognostic selector: diagnostic is a limit of one formula, not a second   !
      !          branch. -----------------------------------------------------------------------------!
      call veg_energy_balance(sw, lw, h, les, lws, ler, t_cas, t_cas, 0.0_wp, t_cas,             &
                                 dtt, ts, tr, dh, drn)
      dt_diag_ref = (sw + lw - ler) / denom
      call check('store_hcap_per_dt = 0 gives the exact diagnostic offset', dtt, dt_diag_ref, 1.0e-12_wp)
      ts_ref = t_cas + dt_diag_ref

      !----- (b) store_hcap_per_dt -> huge FREEZES the store at its entry temperature. -----------------------!
      t0 = 288.0_wp
      call veg_energy_balance(sw, lw, h, les, lws, ler, t_cas, t_cas, 1.0e12_wp, t0,             &
                                 dtt, ts, tr, dh, drn)
      call check('store_hcap_per_dt -> infinity holds the store at t_store0', ts, t0, 1.0e-6_wp)

      !----- (c) THE CONSERVATION IDENTITY. The endpoint and step-average weights are different      !
      !          numbers, and pairing them correctly is exactly what makes the balance close:        !
      !            store_hcap_per_dt*(dt_end - dt_prev) + denom*dt_avg == numer.                     !
      !          dt_avg is not returned, so recover it from dh = h*dt_avg. ------------------------!
      cap  = 4.0e4_wp                       ! a big cohort: tau = cap/denom is order dt
      a_st = cap / dt
      tau  = cap / denom
      call veg_energy_balance(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t0,                  &
                                 dtt, ts, tr, dh, drn)
      call check('energy balance closes with the endpoint/average pair',                            &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-9_wp)

      !----- (d) THE POINT: marching the kernel reproduces the analytic relaxation EXACTLY, not      !
      !          approximately. A backward-Euler storage term would drift here; this does not.       !
      !          T(t) = T_eq + (T0 - T_eq)*exp(-t/tau).                                              !
      nstep = 20_ik ; t_march = t0 ; worst = 0.0_wp
      do k = 1_ik, nstep
         call veg_energy_balance(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t_march,          &
                                    dtt, ts, tr, dh, drn)
         t_march = ts
         t_exact = ts_ref + (t0 - ts_ref) * exp(-real(k, wp) * dt / tau)
         worst   = max(worst, abs(t_march - t_exact))
      end do
      call check_true('marching 20 steps is EXACT against the analytic relaxation', worst < 1.0e-10_wp, worst)

      !----- (e) The small-x series branch (a very large capacity) must not lose precision where     !
      !          1 - exp(-x) would cancel. tau = 1e6 s against dt = 1800 s puts x ~ 2.7e-4. ---------!
      a_st = 1.0e6_wp * denom / dt
      call veg_energy_balance(sw, lw, h, les, lws, ler, t_cas, t_cas, a_st, t0,                  &
                                 dtt, ts, tr, dh, drn)
      call check('large-capacity limit stays accurate (series branch)',                             &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-6_wp)

      !----- (f) A store with NO coupling (x = 0) warms by its forcing over the step, numer*dt/cap,  !
      !          and passes nothing on: the relaxation needs no division by the coupling. ---------!
      a_st = cap / dt
      call veg_energy_balance(sw, lw, 0.0_wp, 0.0_wp, 0.0_wp, ler, t_cas, t_cas, a_st, t0,         &
                                 dtt, ts, tr, dh, drn)
      call check('no coupling: the store warms by numer*dt/cap', ts, t0 + (sw + lw - ler) / a_st, 1.0e-12_wp)
      call check_true('no coupling: no sensible heat leaves', dh == 0.0_wp, dh)
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
      !----- Runs the production cas_column_step_implicit (enthalpy twin;                           !
      !      docs/science/canopy_air_space_biophysics.md eq. (2)). Vapour + CO2 twins held inert (0 !
      !      sources / 0 vapour conductance; a nonzero molar capacity avoids a 0/0 in the unused    !
      !      CO2 branch). The closed-budget residual is assembled here from the production result.   !
      type(cas_column_t) :: column
      type(cas_source_t) :: source
      real(wp) :: enth, shv, temp, resid, enthalpy_atm, worst, enth_new, shv_new, co2_new
      real(wp), parameter :: cas_mass_capacity = 1.2_wp * 20.0_wp     ! rho_air * can_depth  [kg/m2]
      real(wp), parameter :: gatm = 1.2_wp * 0.3_wp * 1.0_wp  ! rho_air * ustar * temp1(c3)  [kg/m2/s]
      integer(ik) :: step
      print '(a)', 'test_cas:'
      shv  = 0.012_wp
      enth = cas_enthalpy_of_temp(295.0_wp, shv)                  ! CAS starts at 295 K
      temp = 295.0_wp
      enthalpy_atm = cas_enthalpy_of_temp(300.0_wp, shv)              ! warmer atmosphere
      column%air_mass_capacity        = cas_mass_capacity
      column%air_molar_capacity       = 1.0_wp                    ! inert CO2 twin (avoid 0/0)
      column%atm_conductance_enthalpy = gatm
      column%atm_conductance_vapor    = 0.0_wp
      column%atm_conductance_co2      = 0.0_wp
      column%atm_enthalpy             = enthalpy_atm
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
      !----- Direct unit test of veg_energy_balance's wet/dry extension (the column tests exercise it   !
      !      only through the full integration). f_wet splits the single latent pathway into a          !
      !      (1-f_wet) DRY (stomatal) share and an f_wet WET (boundary-layer film) share;               !
      !      le_slope_dry/le_ref_dry model the former, le_slope_w/le_ref_w the latter -- two DIFFERENT  !
      !      conductances, so the three properties below are non-trivial (a bug that swapped or mixed   !
      !      them would still "look" plausible).                                                        !
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_temp2, t_store2, transp2, dh2, drnet2
      real(wp), parameter :: abs_sw = 400.0_wp, abs_lw = -50.0_wp, h_coeff = 12.0_wp
      real(wp), parameter :: lw_slope = 6.0_wp, t_cas = 298.0_wp, t_parent = 298.0_wp
      real(wp), parameter :: le_slope_dry = 8.0_wp,  le_ref_dry = 15.0_wp   ! stomatal (dry) pathway
      real(wp), parameter :: le_slope_w   = 20.0_wp, le_ref_w   = 40.0_wp   ! boundary-layer (wet) pathway
      print '(a)', 'test_veg_energy_diagnostic_wetted:'

      !----- 1. ABSENT f_wet must be bit-identical to explicitly passing f_wet=0 -- the contract that   !
      !      lets a caller with no film omit the wet extension. ----------------------------------------!
      call veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_parent, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet)
      call veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_parent, 0.0_wp, t_cas, dt_temp2, t_store2, transp2, dh2, drnet2, &
                                 0.0_wp, le_slope_w, le_ref_w, film_evap)
      call check('f_wet=0 explicit == f_wet absent (dt_temp)', dt_temp2, dt_temp, 1.0e-14_wp)
      call check('f_wet=0 explicit == f_wet absent (transp)',  transp2,  transp,  1.0e-14_wp)
      call check_true('f_wet=0 -> film_evap = 0 exactly', film_evap == 0.0_wp, film_evap)

      !----- 2. FULLY WET (f_wet=1): the dry (stomatal) pathway must vanish -- a saturated film          !
      !      transpires ~0 and evaporates its whole latent capacity through the wet pathway instead.      !
      !      SWAP-EQUIVALENCE: this must be algebraically IDENTICAL to calling the kernel with the wet     !
      !      conductance as the PRIMARY (le_slope,le_ref) argument and no wet extension at all -- not       !
      !      just approximately equal (both reduce to the same 2-unknown linear solve). -------------------!
      call veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_parent, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet,      &
                                 1.0_wp, le_slope_w, le_ref_w, film_evap)
      call veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope_w, lw_slope, le_ref_w,               &
                                 t_cas, t_parent, 0.0_wp, t_cas, dt_temp2, t_store2, transp2, dh2, drnet2)
      call check_true('f_wet=1: dry pathway (transp) vanishes', abs(transp) < 1.0e-14_wp, transp)
      call check('f_wet=1: film_evap == swap-equivalent transp', film_evap, transp2, 1.0e-12_wp)
      call check('f_wet=1: dt_temp == swap-equivalent dt_temp',  dt_temp,   dt_temp2, 1.0e-12_wp)

      !----- 3. PARTIAL wetting (f_wet=0.4): the energy balance must still close -- absorbed radiation     !
      !      net of LW emission (drnet) splits EXACTLY into sensible + BOTH latent pathways, each at the    !
      !      CONSTANT latent_heat_vap (not enthalpy_vapor(T)) -- the identity the surface path relies        !
      !      on to make the CAS's temperature-dependent vapour-enthalpy credit balance against a store's     !
      !      liquid-enthalpy debit (for transpiration the leaf pays the full vapour enthalpy itself; the     !
      !      surface-water accounting does the same for film_evap). -----------------------------------------!
      call veg_energy_balance(abs_sw, abs_lw, h_coeff, le_slope_dry, lw_slope, le_ref_dry,          &
                                 t_cas, t_parent, 0.0_wp, t_cas, dt_temp, t_store, transp, dh, drnet,      &
                                 0.4_wp, le_slope_w, le_ref_w, film_evap)
      call check('partial f_wet: energy balance closes (drnet = dh + latent)', drnet,                  &
                 dh + (transp + film_evap) * latent_heat_vap, 1.0e-9_wp)
      call check_true('partial f_wet: film_evap > 0 (wet pathway active)', film_evap > 0.0_wp, film_evap)
      call check_true('partial f_wet: transp > 0 (dry pathway still active)', transp > 0.0_wp, transp)
   end subroutine test_veg_energy_diagnostic_wetted


   subroutine test_veg_energy_thin_tissue()
      !----- Every term of a tissue's balance -- absorbed radiation, the sapflow's enthalpy, the        !
      !      conductances -- and its heat capacity scale with the tissue's own area                     !
      !      (docs/science/vegetation_energy_dynamics.md eqs. (1), (8)), so a thin tissue is simply a    !
      !      small one. Three properties, swept over an area index from 1e-6 to 10 for leaf and wood,    !
      !      with and without a heat store:                                                             !
      !                                                                                          !
      !        * CONSERVATION: within the step everything that enters leaves or is stored,               !
      !            drnet + q_extra == dh + Lv*(transp + film_evap) + store_hcap_per_dt*(t_store - t0),   !
      !          to round-off. surface_derivs relies on it: the whole-column ledger books coh_rnet on    !
      !          the input side and coh_h/coh_qw as they reach the CAS, and nothing re-checks the        !
      !          tissue in between.                                                                      !
      !        * SPLIT INVARIANCE: halving the area (one cohort of two at half density) keeps the        !
      !          temperature and halves every flux. A floor in ground units on the coupling or the       !
      !          capacity breaks this for a thin tissue, pulling it toward the canopy air.               !
      !        * NO AREA: a tissue with no coupling and no heat capacity takes the temperature of what   !
      !          it grows from (t_parent) and carries no flux. -------------------------------------------!
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_h, t_h, transp_h, dh_h, drnet_h, film_h
      real(wp) :: ai, ast, resid, worst, worst_split
      integer(ik) :: ia, il, is
      real(wp), parameter :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp), parameter :: CAP_PER_AREA = 1500.0_wp   !< [J/m2 tissue/K] a leaf's dry tissue and water
      logical :: is_leaf
      print '(a)', 'test_veg_energy_thin_tissue:'
      worst = 0.0_wp ; worst_split = 0.0_wp
      do is = 0_ik, 1_ik
         do il = 1_ik, 2_ik
            is_leaf = (il == 1_ik)
            do ia = 1_ik, 15_ik
               ai = 1.0e-6_wp * (10.0_wp ** (0.5_wp * real(ia - 1_ik, wp)))
               ast = 0.0_wp ; if (is == 1_ik) ast = CAP_PER_AREA * ai / 900.0_wp
               call thin_tissue(ai, ast, is_leaf, dt_temp, t_store, transp, dh, drnet, film_evap)
               resid = drnet + 25.0_wp*ai - dh - latent_heat_vap*(transp + film_evap) - ast * (t_store - T0)
               worst = max(worst, abs(resid) / (abs(drnet) + 25.0_wp*ai))
               call thin_tissue(0.5_wp*ai, 0.5_wp*ast, is_leaf, dt_h, t_h, transp_h, dh_h, drnet_h, film_h)
               worst_split = max(worst_split, abs(t_h - t_store), abs(2.0_wp*dh_h - dh) / abs(dh),          &
                                 abs(2.0_wp*drnet_h - drnet) / abs(drnet))
            end do
         end do
      end do
      call check_true('a thin tissue conserves energy at every area index (60 cases)', worst < 1.0e-12_wp, worst)
      call check_true('half the area: the same temperature, half of every flux', worst_split < 1.0e-12_wp,       &
                      worst_split)
      !----- No area: no coupling and no heat capacity. -----------------------------------------!
      call thin_tissue(0.0_wp, 0.0_wp, .true., dt_temp, t_store, transp, dh, drnet, film_evap)
      call check('a tissue with no area takes the temperature of what it grows from', t_store, T_PARENT, 0.0_wp)
      call check_true('a tissue with no area carries no flux',                                                  &
                      dh == 0.0_wp .and. drnet == 0.0_wp .and. transp == 0.0_wp .and. film_evap == 0.0_wp, dh)
   end subroutine test_veg_energy_thin_tissue

   !----- One tissue of area index a at 298 K canopy air, starting 0.5 K below it, whose parent sits at  !
   !      301 K: every coupling and forcing term proportional to a. -----------------------------------!
   subroutine thin_tissue(a, a_store, leaf, dtt, ts, tr, h_out, rn, fe)
      real(wp), intent(in)  :: a, a_store
      logical,  intent(in)  :: leaf
      real(wp), intent(out) :: dtt, ts, tr, h_out, rn, fe
      real(wp), parameter   :: GBH = 0.03_wp, RHO = 1.2_wp, EMISS = 0.97_wp
      real(wp), parameter   :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp) :: hc, lws, les, ler
      hc  = merge(2.0_wp, pi, leaf) * a * GBH * RHO * cp_air
      lws = 4.0_wp * EMISS * stefan * T_CAS**3 * a
      les = merge( 8.0_wp * a, 0.0_wp, leaf)    ! stomatal (wood does not transpire)
      ler = merge(15.0_wp * a, 0.0_wp, leaf)
      call veg_energy_balance(400.0_wp*a, -50.0_wp*a, hc, les, lws, ler, T_CAS, T_PARENT, a_store, T0,    &
                              dtt, ts, tr, h_out, rn, f_wet=0.4_wp, le_slope_wet=20.0_wp*a,              &
                              le_ref_wet=40.0_wp*a, film_evap=fe, q_extra=25.0_wp*a)
   end subroutine thin_tissue

end program test_surface_energy
