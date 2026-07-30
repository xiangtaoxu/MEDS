!==========================================================================================!
! test_column_dynamics -- integration test for the fast-loop coupling core: aerodynamics ->     !
! {LEAF GAS EXCHANGE (real GPP + stomata + Rd) + PLANT HYDRAULICS, ground balance} -> soil WATER  !
! column -> CAS three-twin update -> soil THERMAL column, with autotrophic (leaf+stem+root) +      !
! heterotrophic respiration assembling the CAS-CO2 NEE. Driven by a diurnal cycle + morning rain.   !
!   1. CONSERVATION: every fast step closes the CAS energy/water/CO2, soil thermal, soil water,     !
!      AND the WHOLE-COLUMN energy + water budgets (meds_budget_check n_fail == 0).                 !
!   2. PHYSICAL SANITY: CAS temp tracks a diurnal cycle; leaf warms; the soil-surface swing damps   !
!      with depth; soil moisture responds to rain; real photosynthesis draws the CAS CO2 down in    !
!      daylight; the leaf water potential is under tension and more negative at midday.             !
!   3. INTER-LAYER ADVECTION: the hydrology per-face Darcy flux advects liquid enthalpy (always --   !
!      it is a conservation requirement, not an option), and the whole-column budgets still close     !
!      with soil temperatures bounded, so the moisture<->energy coupling conserves.                   !
!==========================================================================================!
program test_column_dynamics
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t, INTEG_SPLIT, INTEG_ARK, INTEG_RK4
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_uext
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG,  &
                                        SOIL_BC_BEDROCK, SOIL_BC_FREE_DRAIN
   use meds_column_state_types, only : build_soil_hydr_params, PSI_INIT
   use meds_column_state_types, only : build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config
   use meds_fast_split,          only : column_fast_step
   use meds_fast_ark,            only : aero_bottom_to_top
   use meds_fast_control,        only : tol_set_t, build_tol_set, GRP_ENTH, GRP_THETA, GRP_SOIL_T
   use meds_fast_dynamics,       only : fast_context_t, build_fast_context
   use meds_hydr_lib,            only : psi_from_water_content, water_content
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik, nstep = 96_ik    ! 96 x 900 s = 24 h
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   !----- Seed moisture + rain-pulse rate, VARIABLES so RUN 7 can drive the column to saturation.  !
   !      Initialized to RUN 1-6's values, which therefore stay byte-identical. --------------------!
   real(wp) :: theta_seed = theta0, rain_pulse = 1.5e-4_wp
   real(wp) :: pond_peak, theta_peak_col
   !----- RUN 8 (snow): seed mass for the pack (0 = no pack, the default for RUNS 1-7) plus the  !
   !      diagnostics integrate_day fills while it is on. ---------------------------------------!
   real(wp) :: snow_seed = 0.0_wp, snow_swe_end, snow_temp_end
   logical  :: snow_physical, snowfall_on = .false.
   real(wp) :: snow_swe_split
   !----- RUN 9 (snowfall with the snow STORE off): the snowfall rate is a VARIABLE so the sub-      !
   !      freezing air and the frozen-precip flux can be toggled INDEPENDENTLY -- the run compares    !
   !      snowf-on against snowf-off at the SAME air temperature, so an evaporation difference        !
   !      cannot masquerade as the water this test is looking for. Default matches RUN 8's rate, so   !
   !      RUN 8 is untouched. `col_water_end` is the whole-column liquid store (soil + pond).         !
   real(wp) :: snowf_rate = 2.0e-5_wp, col_water_end
   real(wp) :: cw_on, cw_off, cw_split_gain, cw_gain, snowf_total
   integer(ik) :: isch
   character(len=9) :: schnm
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: ccfg
   type(column_cohort_t)  :: coh
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: bio
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budg
   type(meds_time_t)      :: sim_date
   real(wp) :: ct_night, ct_noon, co2_night, co2_noon, tleaf_noon, tleaf_night
   real(wp) :: ss_min, ss_max, sd_min, sd_max, th_min, th_max, gpp_noon, nee_noon
   real(wp) :: psileaf_noon, psileaf_night, psileaf_single
   real(wp) :: surf_water_peak     !< RUN 6: running max of leaf+wood interception film over the day
   integer(ik) :: nfail

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- Full model config (PFT traits for leaf gas exchange) + geometry. ----------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp

   !----- Cohort SoA (single cohort): LAI = nplant*leaf_area = 0.3*10 = 3. -----------------!
   call alloc_column_cohort(coh, n)
   coh%pft(1) = 1_ik ; coh%lai(1) = 3.0_wp ; coh%wai(1) = 0.5_wp
   coh%vcmax25(1) = cfg%pft%vcmax25(1) ; coh%rd25(1) = cfg%pft%rd25(1)   ! leaf capacities (plastic trait state)
   coh%height(1) = 16.0_wp ; coh%crown(1) = 0.9_wp
   coh%leaf_width(1) = 0.04_wp ; coh%branch_diam(1) = 0.02_wp
   coh%leaf_area(1) = 10.0_wp ; coh%nplant(1) = 0.3_wp ; coh%dbh(1) = 20.0_wp ; coh%broot(1) = 0.5_wp
   coh%bleaf(1) = 0.5_wp ; coh%bsap(1) = 5.0_wp ; coh%sap_area(1) = 0.01_wp     ! for hydraulic capacitance

   !----- Static column config: soil column + respiration parameters. ---------------------!
   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ccfg%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ccfg%soil_thermal)
   ccfg%wood%is_woody = .true. ; ccfg%wood%stem_resp_factor25 = 0.06_wp ; ccfg%wood%agf_bs = 0.7_wp
   ccfg%root%root_resp_factor25 = 0.30_wp
   ccfg%co2%rh_k_base = 0.01_wp                        ! nonzero decomposition rate so Rh > 0
   ccfg%fast_soil_carbon = 5.0_wp

   !----- Plant hydraulics: flatten cfg%hydraulics -> hydro_p + rhizo + build vuln table. ---!
   call apply_hydraulics_config(cfg%hydraulics, ccfg%hydro_p, ccfg%rhizo_cond)

   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   !=====================================================================================!
   !  RUN 1 -- the full physical-sanity suite. This is now the ONLY soil-thermal coupling      !
   !  there is: interior-face liquid-enthalpy advection used to sit behind                      !
   !  ccfg%advect_soil_heat and RUN 1 ran with it OFF, which is why RUN 2 existed to re-run     !
   !  the same day with it ON. The flag is deleted (it was a conservation requirement, not an   !
   !  accuracy knob, and ARK/RK45 never consulted it), so RUN 2's assertions are folded in      !
   !  below and its slot is retired. RUN numbering is UNCHANGED -- "RUN 7" is referenced from   !
   !  meds_fast_split, meds_fast_types, test_column_ark and test_column_rk45. -------------------!
   !=====================================================================================!
   call integrate_day()
   psileaf_single = psileaf_noon                 ! single-layer (root-frac-weighted BC) baseline for RUN 3

   !----- 1. Conservation: all seven budgets closed every step. ----------------------------!
   call ck(budg%cas_energy%n_fail    == 0_ik, 'CAS energy budget closed',   real(budg%cas_energy%n_fail, wp))
   call ck(budg%cas_water%n_fail     == 0_ik, 'CAS water budget closed',    real(budg%cas_water%n_fail, wp))
   call ck(budg%cas_co2%n_fail       == 0_ik, 'CAS CO2 budget closed',      real(budg%cas_co2%n_fail, wp))
   call ck(budg%soil_energy%n_fail   == 0_ik, 'soil thermal budget closed', real(budg%soil_energy%n_fail, wp))
   call ck(budg%soil_water%n_fail    == 0_ik, 'soil water budget closed',   real(budg%soil_water%n_fail, wp))
   call ck(budg%whole_water%n_fail   == 0_ik, 'WHOLE-COLUMN water budget closes',  real(budg%whole_water%n_fail, wp))
   call ck(budg%whole_energy%n_fail  == 0_ik, 'WHOLE-COLUMN energy budget closes', real(budg%whole_energy%n_fail, wp))

   !----- 2. Physical sanity. -------------------------------------------------------------!
   call ck(ct_noon > ct_night, 'CAS warmer near solar noon than at night', ct_noon - ct_night)
   call ck(tleaf_noon > tleaf_night, 'leaf warms under absorbed shortwave', tleaf_noon - tleaf_night)
   call ck((ss_max - ss_min) > (sd_max - sd_min), 'soil diurnal swing damped with depth', &
           (ss_max - ss_min) - (sd_max - sd_min))
   call ck(th_max - th_min > 1.0e-4_wp, 'soil moisture responds to the rain pulse', th_max - th_min)
   call ck(th_min > 0.05_wp .and. th_max < 0.43_wp, 'soil moisture stays physical', th_max)
   call ck(gpp_noon > 1.0_wp, 'daytime GPP is active (real photosynthesis)', gpp_noon)
   call ck(nee_noon < 0.0_wp, 'daytime NEE is net uptake (GPP > respiration)', nee_noon)
   call ck(co2_noon < co2_night, 'CAS CO2 lower at midday than at night', co2_night - co2_noon)
   call ck(psileaf_noon < 0.0_wp, 'leaf water potential under tension in daylight', psileaf_noon)
   call ck(psileaf_noon < psileaf_night, 'leaf more tensioned at midday than at night',  &
           psileaf_noon - psileaf_night)
   !----- folded in from the retired RUN 2: with the interior faces advecting, soil temperatures  !
   !      stay bounded at BOTH ends of the column. Forcing that advection off reached 361 K here. -!
   call ck(ss_min > 270.0_wp .and. ss_max < 330.0_wp, 'soil surface temp stays bounded', ss_max)
   call ck(sd_min > 270.0_wp .and. sd_max < 330.0_wp, 'deep soil temp stays bounded',    sd_max)

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: RUN 1 PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS temp night=', ct_night, ' K  noon=', ct_noon, ' K)'
      print '(a,f7.2,a,f7.2,a)', '   (CAS CO2  night=', co2_night, '     noon=', co2_noon, ' umol/mol)'
      print '(a,f7.2,a,f7.2,a)', '   (noon GPP=', gpp_noon, ' umol/m2/s  NEE=', nee_noon, ' umol/m2/s)'
      print '(a,f7.3,a,f7.3,a)', '   (leaf psi night=', psileaf_night, ' MPa  noon=', psileaf_noon, ' MPa)'
      print '(a,es10.3,a,es10.3,a)', '   (whole-column worst resid: energy=', budg%whole_energy%worst,      &
                                     ' J/m2  water=', budg%whole_water%worst, ' kg/m2)'
   end if

   !=====================================================================================!
   !  RUN 2 -- RETIRED, slot deliberately left empty.                                          !
   !                                                                                           !
   !  It re-ran RUN 1's day with ccfg%advect_soil_heat = .true. to show the interior-face        !
   !  advection conserves. With the flag deleted RUN 1 already runs that configuration, so this  !
   !  was the identical day integrated twice; its two unique assertions (soil surface and deep   !
   !  temperature bounded) moved into RUN 1. Kept as a numbered gap rather than renumbering       !
   !  RUNS 3-8, because "RUN 7" is cross-referenced from four other files.                        !
   !                                                                                             !
   !  What it established, worth keeping: with the interior faces LUMPED the same scenario        !
   !  reaches a soil surface of 361 K against 297 K with them connected, while the whole-column   !
   !  ledger closes eitherway -- the error is purely vertical. That is the same blind spot that   !
   !  hid issue #78 item 3; see docs/science/numerical_scheme.md sec 7. -------------------------!
   !=====================================================================================!

   !=====================================================================================!
   !  RUN 3 -- opt-in multi-layer root coupling: per-layer soil psi + rhizosphere conductance !
   !           feed the plant boundary. Conservation must still hold and the plant psi must    !
   !           respond vs the single root-frac-weighted BC (RUN 1 baseline).                   !
   !=====================================================================================!
   ccfg%multilayer_roots   = .true.
   ccfg%specific_root_area = cfg%hydraulics%specific_root_area
   call integrate_day()
   call ck(budg%whole_water%n_fail  == 0_ik, 'MULTILAYER: whole-column water still closes',  &
           real(budg%whole_water%n_fail, wp))
   call ck(budg%whole_energy%n_fail == 0_ik, 'MULTILAYER: whole-column energy still closes', &
           real(budg%whole_energy%n_fail, wp))
   call ck(psileaf_noon < 0.0_wp, 'MULTILAYER: leaf psi still under tension', psileaf_noon)
   call ck(abs(psileaf_noon - psileaf_single) > 1.0e-9_wp,                                   &
           'MULTILAYER: per-layer coupling shifts leaf psi vs single-BC', psileaf_noon - psileaf_single)
   ccfg%multilayer_roots = .false.                ! restore for the cohort-order test below

   !=====================================================================================!
   !  RUN 4 -- caller-side cohort ORDER: aero_bottom_to_top must respect the wind cascade.  !
   !=====================================================================================!
   call test_aero_order()

   !=====================================================================================!
   !  RUN 5 -- goal (a) Layer 1: the sub-solver TOLERANCE UNIFICATION plumbing.            !
   !=====================================================================================!
   call test_tolerance_unification()

   !=====================================================================================!
   !  RUN 6 -- opt-in CANOPY-SURFACE WATER (MEDS_ED2_RK45_DESIGN.md sec 3.4, P1): reruns the   !
   !           SAME diurnal-cycle-plus-morning-rain-pulse forcing (RUN 1's day, istep 12-28)     !
   !           with ccfg%canopy_water_on = .true. whole_water must close exactly (the headline     !
   !           gate, sec 8 gate 2) and the canopy must actually have intercepted some of the rain    !
   !           (surf_water_peak > 0) -- proving the wiring, not just that the feature is a silent     !
   !           no-op. The wetted-fraction transpiration-suppression ALGEBRA itself is unit-tested      !
   !           directly in test_surface_energy.f90 (no direct test existed for veg_energy_diagnostic   !
   !           before P1); this run only needs to prove the WIRING (interception->film->CAS->ledgers). !
   !                                                                                          !
   !           whole_energy is checked against a LOOSE, explicit bound rather than n_fail==0: making    !
   !           it close to machine precision needs the surface water's OWN prognostic temperature/       !
   !           heat capacity (sec 3.4's own "d(leaf_energy)/dt gains the film's storage term"), which     !
   !           is real thermal-inertia machinery this pass does not build -- this pass instead values      !
   !           the store at ONE fixed reference (rain_temp) at both endpoints, which correctly closes      !
   !           MASS and the LATENT-heat exchange (verified: film_evap*latent_heat_vap balances against      !
   !           coh_rnet exactly, same identity test_surface_energy.f90 unit-tests) but leaves a residual      !
   !           proportional to the water's SENSIBLE heat at the leaf/wood dt_temp offset uncounted -- the      !
   !           same category of deferred upwind-temperature imprecision as sec 2's qloss/qwflux_wl         !
   !           coupling (also explicitly P2-deferred energy-advection work, not this pass's gate).           !
   !=====================================================================================!
   ccfg%canopy_water_on = .true.
   call integrate_day()
   ccfg%canopy_water_on = .false.                 ! restore default for any future test added after this
   call ck(budg%whole_water%n_fail  == 0_ik, 'CANOPY WATER: whole-column water still closes',   &
           real(budg%whole_water%n_fail, wp))
   call ck(surf_water_peak > 0.0_wp, 'CANOPY WATER: the morning rain pulse was actually intercepted', &
           surf_water_peak)
   call ck(budg%whole_energy%worst < 5.0e6_wp,                                                    &
           'CANOPY WATER: whole-column energy stays BOUNDED (known deferred sensible-heat approx)', &
           budg%whole_energy%worst)
   if (nfail == 0_ik) then
      print '(a,es10.3,a)', '   (RUN 6 peak canopy film water=', surf_water_peak, ' kg/m2)'
      print '(a,es10.3,a)', '   (RUN 6 worst whole_energy resid=', budg%whole_energy%worst, ' J/m2)'
   end if

   !=====================================================================================!
   !  RUN 7 -- SATURATED column: exercises the two water-enthalpy paths RUNS 1-6 never take. !
   !                                                                                          !
   !  A sealed (bedrock) column seeded just below theta_sat under heavy rain must both          !
   !  SATURATE -- firing meds_soil_water's post-solve theta clip, which moves water with no      !
   !  face -- and OVERFLOW its ponding store, producing surface runoff. Neither happens in       !
   !  RUNS 1-6 (free-draining, theta0 = 0.30, gentle pulse), so the enthalpy bookkeeping for      !
   !  both was previously unexercised by any test:                                               !
   !    * the clip's per-layer enthalpy compensation (soil_energy debited at the layer's OWN      !
   !      temperature, paired into e_out), and                                                   !
   !    * runoff, which must carry NO enthalpy term at all -- the ponding store is a mass         !
   !      buffer with no enthalpy state, so debiting soil layer 1 for it (as the code used to)    !
   !      removed ~1 MJ per kg the layer never received.                                          !
   !  whole_energy closing to n_fail == 0 here is the assertion that both are booked with the     !
   !  right sign and magnitude: an unpaired term of either kind shows up directly in the residual.!
   !                                                                                          !
   !  The interior faces advect liquid enthalpy, unconditionally -- and THIS run is what showed     !
   !  why that must not be optional. Back when it sat behind a flag, lumping the interior faces      !
   !  took this same scenario to a soil surface of 361 K against 297 K with them connected. The      !
   !  boundary faces put liquid enthalpy into layer 1 and take it out of whichever layer the         !
   !  clip fires in; with no interior advection there is no path between them, so layer 1 gains      !
   !  the full infiltration enthalpy (~2.5e6 J/m2/step here) while a deeper layer sheds it. The      !
   !  WHOLE-column ledger still closes -- the error is purely in the vertical distribution, which    !
   !  is exactly what a column-vs-boundary sum cannot see, and is the same blind spot that hid      !
   !  issue #78 item 3. A conservation requirement, not an accuracy knob; the flag is deleted.      !
   !=====================================================================================!
   theta_seed = 0.425_wp                          ! just below theta_sat = 0.43
   rain_pulse = 2.0e-3_wp                         ! heavy: far above what a sealed column can absorb
   ccfg%hydro%bottom_bc = SOIL_BC_BEDROCK         ! sealed: the water has nowhere to drain
   call integrate_day()
   call ck(theta_peak_col >= 0.43_wp - 1.0e-12_wp,                                              &
           'SATURATED: the column reached theta_sat (clip path is live)', theta_peak_col)
   call ck(pond_peak >= ccfg%hydro%w_pond_max - 1.0e-9_wp,                                      &
           'SATURATED: ponding store filled and overflowed (runoff path is live)', pond_peak)
   call ck(budg%whole_water%n_fail  == 0_ik, 'SATURATED: whole-column water still closes',      &
           real(budg%whole_water%n_fail, wp))
   call ck(budg%whole_energy%n_fail == 0_ik, 'SATURATED: whole-column energy still closes',     &
           real(budg%whole_energy%n_fail, wp))
   call ck(ss_min > 250.0_wp .and. ss_max < 340.0_wp,                                           &
           'SATURATED: soil surface temp stays physical through clip + runoff', ss_max)
   if (nfail == 0_ik) then
      print '(a,f6.3,a,f6.3,a)', '   (RUN 7 peak theta=', theta_peak_col, '  peak pond=', pond_peak, ' kg/m2)'
      print '(a,es10.3,a,f7.2,a,f7.2,a)', '   (RUN 7 worst whole_energy resid=',                 &
            budg%whole_energy%worst, ' J/m2  soil surf ', ss_min, '-', ss_max, ' K)'
   end if

   !=====================================================================================!
   !  RUN 8 -- SNOW ON. The split path's snow COUPLING had no integration test at all:          !
   !  test_snow.f90 exercises the kernels standalone and nothing in the suite ever set           !
   !  ccfg%snow_on, so the operator-split snow stage inside column_fast_step -- accumulate ->     !
   !  surface balance -> meltwater drain -> snowfac-blended ground -> sublimation into the CAS     !
   !  -> the swe and snow_acc_enth terms in the whole-column ledgers -- was entirely uncovered.    !
   !                                                                                          !
   !  That matters beyond the usual reasons. MEDS_INTEGRATOR_PARITY.md row 2 (C4) plans to HOIST   !
   !  this stage out of column_fast_step into a routine all three integrators share, and migrating !
   !  coupling code with no regression net is migrating blind -- "snow-off bit-identical" would     !
   !  only prove the OFF path survived, which is the half that cannot break. The two budget         !
   !  assertions below are the net: if the hoist ever drops the swe delta, the pack's precip        !
   !  enthalpy, or the sublimation vapour, an unpaired term shows up directly in the residual.      !
   !=====================================================================================!
   theta_seed = theta0 ; rain_pulse = 0.0_wp      ! no rain: snowfall is the only water input
   ccfg%hydro%bottom_bc = SOIL_BC_FREE_DRAIN
   !----- seeded deep enough that the pack SURVIVES the day under every scheme. At 20 kg/m2 it sat on !
   !      a knife-edge: split ended at 1.44 kg/m2 while ARK/RK45 exhausted it exactly, so a "pack      !
   !      still present" assertion was really testing the last ~7% of the melt energy rather than the  !
   !      wiring. A surviving pack lets the schemes be compared on melt TOTAL, which is the quantity   !
   !      that actually says whether they are driving the same stage. --------------------------------!
   snow_seed    = 60.0_wp
   snowfall_on  = .true.
   call integrate_day()
   snow_swe_split = snow_swe_end
   snow_seed = 0.0_wp ; snowfall_on = .false.
   call ck(snow_physical, 'SNOW ON: pack + soil stay physical over a 24 h march', snow_temp_end)
   !----- The pack MELTS over this day rather than growing (75 W/m2 of ground shortwave with no       !
   !      longwave loss in this fixture beats 2e-5 kg/m2/s of snowfall), so assert the melt path is    !
   !      live rather than accumulation. That is the better coverage of the two: melt is what routes   !
   !      a PAIRED (mass, enthalpy) transfer into the soil top, which is exactly the seam a hoist is   !
   !      most likely to drop. It stays strictly positive, so a pack is present for the whole march    !
   !      and the snowfac-blended surface is exercised throughout. -----------------------------------!
   call ck(snow_swe_end < 60.0_wp .and. snow_swe_end > 0.0_wp,                                      &
           'SNOW ON: melt/sublimation drained the pack without exhausting it', snow_swe_end)
   call ck(budg%whole_water%n_fail  == 0_ik, 'SNOW ON: whole-column water closes with a pack',      &
           real(budg%whole_water%n_fail, wp))
   call ck(budg%whole_energy%n_fail == 0_ik, 'SNOW ON: whole-column energy closes with a pack',     &
           real(budg%whole_energy%n_fail, wp))
   if (nfail == 0_ik) then
      print '(a,f7.3,a,f7.2,a)', '   (RUN 8 snow: swe=', snow_swe_end, ' kg/m2  T_snow=', snow_temp_end, ' K)'
      print '(a,es10.3,a,es10.3,a)', '   (RUN 8 worst whole-column resid: energy=',                 &
            budg%whole_energy%worst, ' J/m2  water=', budg%whole_water%worst, ' kg/m2)'
   end if

   !----- C4: the SAME scenario under ARK and RK45. Before the hoist both imported the snow kernels   !
   !      and never called them, so this run was silently snow-free under each. These assertions are  !
   !      the SUFFICIENT half of the migration's success criterion -- snow-off bit-identity (which    !
   !      the rest of the suite covers) only proves the OFF path survived, which is the half that     !
   !      cannot break. ------------------------------------------------------------------------!
   do isch = 1_ik, 2_ik
      if (isch == 1_ik) then
         cfg%time_integrator = INTEG_ARK  ; schnm = 'SNOW ARK '
      else
         cfg%time_integrator = INTEG_RK4  ; schnm = 'SNOW RK45'
      end if
      !----- re-arm the scenario: split's block above switches snow back OFF when it finishes, so     !
      !      without this the ARK/RK45 runs are silently snow-FREE -- which looks like success (their !
      !      ledgers close trivially) and is the exact failure mode these assertions exist to catch.  !
      snow_seed = 60.0_wp ; snowfall_on = .true.
      call integrate_day()
      call ck(snow_physical, trim(schnm)//': pack + soil stay physical with a pack', snow_temp_end)
      call ck(snow_swe_end < 60.0_wp .and. snow_swe_end > 0.0_wp,                                    &
              trim(schnm)//': melt/sublimation drained the pack without exhausting it', snow_swe_end)
      call ck(budg%whole_water%n_fail  == 0_ik,                                                      &
              trim(schnm)//': whole-column water closes with a pack', real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik,                                                      &
              trim(schnm)//': whole-column energy closes with a pack', real(budg%whole_energy%n_fail, wp))
      !----- and the three must AGREE on the pack, not merely each conserve on their own: a shared     !
      !      stage fed different inputs by one caller would still close every ledger independently.    !
      !----- compare the MELT TOTAL, not the remaining mass: melt is what the stage computed, while    !
      !      the remainder is melt subtracted from a common seed and so exaggerates small differences  !
      !      once the pack is nearly gone. 10% allows for the schemes genuinely driving the surface    !
      !      balance along different CAS trajectories, while a wiring error (a stage fed the wrong     !
      !      inputs, or not called) is order-100%. ------------------------------------------------!
      call ck(abs((60.0_wp - snow_swe_end) - (60.0_wp - snow_swe_split))                             &
              < 0.10_wp * (60.0_wp - snow_swe_split),                                                 &
              trim(schnm)//': melt total agrees with split (one shared stage, same inputs)',          &
              abs(snow_swe_end - snow_swe_split))
      if (nfail == 0_ik) print '(3a,f7.3,a,es10.3)', '   (RUN 8 ', trim(schnm), ' swe=', snow_swe_end, &
            ' kg/m2  worst water resid=', budg%whole_water%worst
   end do
   cfg%time_integrator = INTEG_SPLIT

   !=====================================================================================!
   !  RUN 9 -- SNOWFALL IS CONSERVED, ON EVERY INTEGRATOR.                                      !
   !                                                                                          !
   !  Snowfall is a boundary water input like rain, and every kg of it must end up somewhere in  !
   !  the column: the pack, the soil, or the pond. `build_column_frozen` used to drop it (it       !
   !  routed only `forc%precip`, never `forc%snowf`), so ARK and RK45 lost `snowf*dt` every         !
   !  sub-freezing step while the whole-column ledger counted it as an input                        !
   !  (MEDS_INTEGRATOR_PARITY.md sec 3e, E-4).                                                       !
   !                                                                                          !
   !  RUN 8 cannot catch that. It asserts on the PACK -- melt totals and pack mass -- which stay     !
   !  right even when the routing of what LEAVES the pack is wrong. This run asserts on the whole    !
   !  column instead: soil liquid + pond + pack.                                                     !
   !                                                                                          !
   !  THE ASSERTION IS A BOUNDARY-INPUT IDENTITY, NOT A LEDGER RESIDUAL, and that is deliberate.     !
   !  A per-step `snowf*dt` sits below the whole-column closure tolerance even while it integrates   !
   !  into a large seasonal error -- a forced month on the broken code ran to completion without     !
   !  tripping the budget hard stop. So this run asks a question the ledger cannot: integrate the    !
   !  SAME day twice, snowfall on and off, holding the air temperature (and hence melt, sublimation  !
   !  and drainage) common, and require the column to gain what fell. A scheme that discards it      !
   !  gains ~0 -- an order-100% failure with no tolerance to tune. A pack is seeded in BOTH runs so   !
   !  that sublimation, a genuine boundary OUTPUT, is common to the pair and cancels in the           !
   !  difference.                                                                                     !
   !=====================================================================================!
   snow_seed = 60.0_wp ; rain_pulse = 0.0_wp ; theta_seed = theta0
   ccfg%hydro%bottom_bc = SOIL_BC_FREE_DRAIN
   snowf_total = 2.0e-5_wp * real(nstep, wp) * dt_fast     ! [kg/m2] the day's frozen-precip input
   do isch = 1_ik, 3_ik
      select case (isch)
      case (1_ik) ; cfg%time_integrator = INTEG_SPLIT ; schnm = 'SNOWF spl'
      case (2_ik) ; cfg%time_integrator = INTEG_ARK   ; schnm = 'SNOWF ark'
      case default; cfg%time_integrator = INTEG_RK4   ; schnm = 'SNOWF r45'
      end select
      !----- snowfall ON: cold air + frozen precip. -------------------------------------------!
      snowfall_on = .true. ; snowf_rate = 2.0e-5_wp
      call integrate_day()
      cw_on = col_water_end
      call ck(budg%whole_water%n_fail  == 0_ik, trim(schnm)//': whole-column water closes',        &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik, trim(schnm)//': whole-column energy closes',       &
              real(budg%whole_energy%n_fail, wp))
      !----- snowfall OFF: the SAME cold air and the SAME seeded pack, so melt, sublimation and    !
      !      drainage are common to the pair and cancel in the difference. What is left is the      !
      !      frozen-precip input alone. ---------------------------------------------------------!
      snowf_rate = 0.0_wp
      call integrate_day()
      cw_off = col_water_end
      snowf_rate = 2.0e-5_wp ; snowfall_on = .false.
      cw_gain = cw_on - cw_off
      if (isch == 1_ik) cw_split_gain = cw_gain
      !----- 10% allows for the second-order response of melt and drainage to the slightly deeper   !
      !      pack; the defect this guards against is order-100% (the water never arrives). ---------!
      call ck(abs(cw_gain - snowf_total) < 0.10_wp * snowf_total,                                   &
              trim(schnm)//': snowfall is conserved into the column', cw_gain)
      call ck(abs(cw_gain - cw_split_gain) < 0.10_wp * snowf_total,                                 &
              trim(schnm)//': snowfall gain agrees with split', cw_gain - cw_split_gain)
      if (nfail == 0_ik) print '(3a,f8.4,a,f8.4,a)', '   (RUN 9 ', trim(schnm), ' column water gain=', &
            cw_gain, ' kg/m2  expected=', snowf_total, ' kg/m2)'
   end do
   cfg%time_integrator = INTEG_SPLIT
   snow_seed = 0.0_wp

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: RUNS 3-8 PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS noon=', ct_noon, ' K  soil surf max=', ss_max, ' K)'
      print '(a,es10.3,a,es10.3,a)', '   (whole-column worst resid: energy=', budg%whole_energy%worst,       &
                                     ' J/m2  water=', budg%whole_water%worst, ' kg/m2)'
      print '(a)', 'test_column_dynamics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_dynamics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- One 24 h diurnal integration from a freshly seeded column state. Fills the host    !
   !      diagnostic variables (min/max ranges, noon/night captures) + budg. ----------------!
   subroutine integrate_day()
      real(wp)    :: t_sec, cosz, t_air
      integer(ik) :: istep, k


      !----- (Re)seed the prognostic column state. --------------------------------------!
      if (allocated(bio%leaf_temp)) deallocate(bio%leaf_temp)
      call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
      !----- alloc_patch_biophys seeds leaf_water_mass/wood_water_mass at a scratch 0 (the real     !
      !      lazy-init lives in meds_fast_dynamics.f90's site-level gather loop, which this driver-  !
      !      level test bypasses) -- seed the same water_content(PSI_INIT,...) a freshly-created     !
      !      cohort gets there, or psi_from_water_content would diagnose an unphysical psi from an   !
      !      empty pool. -------------------------------------------------------------------------!
      bio%leaf_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%leaf_pi0, ccfg%hydro_p%leaf_elastic_mod, &
           ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, coh%bleaf(1:n))
      bio%wood_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%wood_pi0, ccfg%hydro_p%wood_elastic_mod, &
           ccfg%hydro_p%wood_apoplast_frac, ccfg%hydro_p%wood_water_sat, coh%bsap(1:n) + coh%broot(1:n))
      budg = column_budget_t()
      !----- RUN 8: seed a snow pack when asked. snow_seed = 0 (RUNS 1-7) leaves nlayer = 0, which is  !
      !      exactly the no-pack state alloc_patch_biophys already produces, so those runs are          !
      !      untouched. rho_snow = 250 kg/m3 matches meds_main's own seeding. ------------------------!
      if (snow_seed > 0.0_wp) then
         bio%snow%swe(1)         = snow_seed
         bio%snow%snow_energy(1) = temp_to_uext(0.0_wp, snow_seed, 270.0_wp, 0.0_wp)
         bio%snow%snow_depth(1)  = snow_seed / 250.0_wp
         bio%snow%nlayer         = 1_ik
      end if
      snow_physical = .true.
      bio%soil_w%theta(1:nsl) = theta_seed
      do k = 1_ik, nsl
         bio%soil_e%soil_energy(k) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(k),  &
                                     theta_seed * rho_h2o, t0, 1.0_wp)
         bio%soil_e%soil_temp(k)   = t0
      end do

      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp ; sd_min = 1.0e9_wp ; sd_max = -1.0e9_wp
      th_min = 1.0e9_wp ; th_max = -1.0e9_wp
      surf_water_peak = -1.0e9_wp
      pond_peak = 0.0_wp ; theta_peak_col = 0.0_wp

      do istep = 1_ik, nstep
         t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
         cosz  = solar_cosz(sim_date, t_sec, lat)
         t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)

         forc%abs_sw   = 500.0_wp * cosz                            ! leaf-absorbed shortwave [W/m2]
         forc%abs_par  = forc%abs_sw                                 ! MVP: absorbed PAR == absorbed SW (par_per_w=2.1)
         forc%abs_lw   = 0.0_wp
         forc%abs_sw_ground = 75.0_wp * cosz
         forc%abs_lw_ground = 0.0_wp
         forc%precip   = 0.0_wp
         if (istep >= 12_ik .and. istep <= 28_ik) forc%precip = rain_pulse     ! morning rain pulse
         !----- RUN 8: steady light snowfall so the pack GROWS (tests the accumulate path and the    !
         !      pack's precip-enthalpy boundary term). 0 for every other run. ----------------------!
         forc%snowf = 0.0_wp
         if (snowfall_on) forc%snowf = snowf_rate
         if (snowfall_on) t_air = 268.0_wp + 3.0_wp * (cosz - 0.3_wp)   ! sub-freezing: the pack must survive
         forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
         forc%shv_atm      = 0.008_wp
         forc%co2_atm      = 400.0_wp

         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg)

         ss_min = min(ss_min, bio%soil_e%soil_temp(1))   ; ss_max = max(ss_max, bio%soil_e%soil_temp(1))
         sd_min = min(sd_min, bio%soil_e%soil_temp(nsl)) ; sd_max = max(sd_max, bio%soil_e%soil_temp(nsl))
         th_min = min(th_min, bio%soil_w%theta(1))       ; th_max = max(th_max, bio%soil_w%theta(1))
         surf_water_peak = max(surf_water_peak, bio%leaf_surf_water(1) + bio%wood_surf_water(1))
         pond_peak       = max(pond_peak, bio%soil_w%w_surface)
         theta_peak_col  = max(theta_peak_col, maxval(bio%soil_w%theta(1:nsl)))
         if (bio%snow%nlayer >= 1_ik) then
            snow_physical = snow_physical .and. bio%snow%swe(1) >= 0.0_wp                            &
                            .and. bio%snow%snow_temp(1) > 200.0_wp                                   &
                            .and. bio%snow%snow_temp(1) < 290.0_wp
            do k = 1_ik, nsl
               snow_physical = snow_physical .and. bio%soil_e%soil_temp(k) > 230.0_wp                &
                               .and. bio%soil_e%soil_temp(k) < 340.0_wp
            end do
         end if
         if (istep == 54_ik) then
            ct_noon = bio%cas%can_temp ; tleaf_noon = bio%leaf_temp(1) ; co2_noon = bio%cas%can_co2
            gpp_noon = budg%gpp_last ; nee_noon = budg%nee_last
            !----- psi is no longer persisted state (MEDS_ED2_RK45_DESIGN.md sec 4): diagnose it   !
            !      from the persisted leaf_water_mass. --------------------------------------------!
            psileaf_noon = psi_from_water_content(bio%leaf_water_mass(1), ccfg%hydro_p%leaf_pi0,     &
                 ccfg%hydro_p%leaf_elastic_mod, ccfg%hydro_p%leaf_apoplast_frac,                       &
                 ccfg%hydro_p%leaf_water_sat, coh%bleaf(1))
         end if
         if (istep == 2_ik) then
            ct_night = bio%cas%can_temp ; tleaf_night = bio%leaf_temp(1) ; co2_night = bio%cas%can_co2
            psileaf_night = psi_from_water_content(bio%leaf_water_mass(1), ccfg%hydro_p%leaf_pi0,    &
                 ccfg%hydro_p%leaf_elastic_mod, ccfg%hydro_p%leaf_apoplast_frac,                       &
                 ccfg%hydro_p%leaf_water_sat, coh%bleaf(1))
         end if
      end do
      snow_swe_end = bio%snow%swe(1) ; snow_temp_end = bio%snow%snow_temp(1)
      !----- whole-column LIQUID store at the end of the march (soil layers + ponding), for RUN 9's  !
      !      boundary-input accounting. Not a budget -- just the store the arriving water must land   !
      !      in if it lands anywhere at all. ----------------------------------------------------!
      col_water_end = bio%soil_w%w_surface + bio%snow%swe(1)
      do k = 1_ik, nsl
         col_water_end = col_water_end + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
   end subroutine integrate_day

   subroutine ck(cond, name, val)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: val
      if (.not. cond) then
         print '(a,a,a,es14.6)', '  FAIL ', name, ': val = ', val
         nfail = nfail + 1_ik
      end if
   end subroutine ck

   !----- Guard the caller-side cohort-order fix. aero_bottom_to_top feeds canopy_aerodynamics the !
   !      BOTTOM(1)->TOP(n) order it contracts for, starting from a height-DESCENDING column buffer  !
   !      (index 1 = tallest). So the TALL cohort (gather index 1 = canopy top) must come back with   !
   !      the top-of-canopy wind (=> higher boundary-layer gb); the short understory gets the         !
   !      attenuated wind. The old direct call (gather order straight into the kernel) inverts this,   !
   !      and reverting the fix must fail here.                                                        !
   subroutine test_aero_order()
      type(column_cohort_t) :: c2
      type(aero_out_t)      :: a2
      type(aero_env_t)      :: e2      ! defaults are fine (u_ref/zref/press/rho_air/can_* all sensible)
      type(aero_geom_t)     :: g2
      real(wp)              :: lt(2)
      call alloc_column_cohort(c2, 2_ik)
      c2%height = [18.0_wp, 6.0_wp] ; c2%lai = [2.0_wp, 2.0_wp] ; c2%crown = [0.9_wp, 0.8_wp]   ! DESCENDING
      c2%leaf_width = [0.04_wp, 0.04_wp] ; c2%branch_diam = [0.02_wp, 0.02_wp]
      lt = [300.0_wp, 300.0_wp]
      e2%u_ref = 3.0_wp
      g2%veg_height = 18.0_wp ; g2%opencan_frac = 0.0_wp ; g2%snowfac = 0.0_wp
      call alloc_aero_out(a2, 2_ik)
      call aero_bottom_to_top(ccfg%aero, e2, g2, 2_ik, c2, lt, a2)
      call ck(a2%wind(1) > a2%wind(2), 'aero order: tall cohort (gather idx1=top) gets more wind',    &
              a2%wind(1) - a2%wind(2))
      call ck(a2%leaf_gbw(1) > a2%leaf_gbw(2), 'aero order: tall cohort gets higher leaf gb',         &
              a2%leaf_gbw(1) - a2%leaf_gbw(2))
   end subroutine test_aero_order

   !----- Goal (a) Layer 1. The fast loop has FOUR adaptive error estimates (the ARK march plus the  !
   !      nested soil-water / soil-energy / plant-hydraulics sub-solvers), and before the unification  !
   !      each carried its own tolerance -- two of them unreachable from config. build_tol_set is the   !
   !      single authority; build_fast_context pushes its groups down into the sub-solver opts.         !
   !      This is a PLUMBING test on purpose: a whole-model A/B cannot see it (a tolerance change that   !
   !      does not flip a substep decision is silently byte-identical), so the wiring is asserted here   !
   !      directly. It also pins the two accuracy dials, including the atol_scale companion -- rtol_all  !
   !      alone saturates once atol dominates the atol + rtol*|y| denominator.                           !
   subroutine test_tolerance_unification()
      type(meds_config_t)   :: c
      type(fast_context_t)  :: fx
      type(tol_set_t)       :: t
      c = build_test_config()
      !----- (a) DEFAULT: each group keeps the tolerance that governs it today (byte-identical). ------!
      t = build_tol_set(c)
      call ck(t%rtol(GRP_THETA)  == c%soil%rtol,   'tol: theta group seeded from [soil].rtol',   t%rtol(GRP_THETA))
      call ck(t%atol(GRP_THETA)  == c%soil%atol,   'tol: theta group seeded from [soil].atol',   t%atol(GRP_THETA))
      call ck(t%rtol(GRP_SOIL_T) == c%energy%rtol, 'tol: soil-T group seeded from [energy].rtol', t%rtol(GRP_SOIL_T))
      call ck(t%rtol(GRP_ENTH)   == c%ark_rtol,    'tol: ARK groups seeded from ark_rtol',        t%rtol(GRP_ENTH))
      !----- (b) MASTER DIALS: rtol_all overrides every group; atol_scale multiplies every atol. ------!
      c%rtol_all   = 1.0e-7_wp
      c%atol_scale = 1.0e-2_wp
      t = build_tol_set(c)
      call ck(all(t%rtol == 1.0e-7_wp), 'tol: rtol_all overrides ALL groups', maxval(abs(t%rtol - 1.0e-7_wp)))
      call ck(t%atol(GRP_THETA) == c%soil%atol * 1.0e-2_wp, 'tol: atol_scale scales atol', t%atol(GRP_THETA))
      !----- (c) PUSH-DOWN: the dials must reach the nested sub-solvers, not just the ARK march. Note:  !
      !      the plant-hydraulics sub-solver (hydro_o) is DELIBERATELY no longer pushed from here          !
      !      (MEDS_ED2_RK45_DESIGN.md sec 4/6, P2): its retired outer group (GRP_PSI, psi-space [MPa])      !
      !      became GRP_LEAF_W/GRP_WOOD_W (mass-space [kg/plant]) when internal water mass replaced        !
      !      psi as the fast-loop prognostic state, but hydro_o's OWN internal step-doubling still          !
      !      operates in psi space (solve_plant_water's matrix exponential) -- feeding it from a mass-       !
      !      space group would be a unit mismatch, not a unification, so it now keeps its own type          !
      !      default (build_fast_context no longer touches it at all). --------------------------------!
      call build_fast_context(c, fx)
      call ck(fx%ccfg%hydro%rtol   == 1.0e-7_wp, 'tol: dial reaches the soil-WATER sub-solver',  fx%ccfg%hydro%rtol)
      call ck(fx%ccfg%energy%rtol  == 1.0e-7_wp, 'tol: dial reaches the soil-ENERGY sub-solver', fx%ccfg%energy%rtol)
   end subroutine test_tolerance_unification

end program test_column_dynamics
