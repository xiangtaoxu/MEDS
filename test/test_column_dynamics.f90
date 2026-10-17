! SPDX-License-Identifier: Apache-2.0
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
   use meds_constants,           only : latent_heat_fusion, rho_h2o, grav, pio4
   use meds_config,              only : meds_config_t, INTEG_ARK, INTEG_RK45
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_internal_energy
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_fast_types, only : patch_biophys_t, alloc_patch_biophys
   use meds_column_view,       only : column_cohort_init
   use meds_water_retention, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : SOIL_BC_BEDROCK, SOIL_BC_FREE_DRAIN
   use meds_canopy_types, only : set_aero_env_atm, set_aero_env_canopy
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_column_params, only : n_soil_layer_max
   use meds_column_params, only : PSI_INIT, build_soil_hydr_params, build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config, apply_canopy_film_config, &
                                        fast_workspace_t
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts, build_tol_set
   use meds_fast_step,          only : column_fast_step
   use meds_fast_prepass,        only : aero_bottom_to_top, cohort_leaf_par
   use meds_fast_types,          only : tol_set_t, GRP_ENTH, GRP_THETA
   use meds_fast_dynamics,       only : fast_context_t, build_fast_context
   use meds_water_retention,     only : psi_from_water_content, water_content
   use meds_test_support, only : build_test_config, check_true, test_report
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks were written on

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik
   real(wp),    parameter :: dt_fast = 150.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   !----- Durations and times of day are written in hours and turned into steps here, so a change of  !
   !      dt_fast keeps what each test covers: the march puts step i at (i - 1/2)*dt_fast after       !
   !      midnight, and the sun follows the clock. ---------------------------------------------------!
   integer(ik), parameter :: steps_per_hour = nint(3600.0_wp / dt_fast, ik)
   integer(ik), parameter :: nstep = 24_ik * steps_per_hour      ! one day
   !----- Seed moisture + rain-pulse rate, VARIABLES so RUN 7 can drive the column to saturation.  !
   !      Initialized to the values RUNS 1-6 run on. ----------------------------------------------!
   real(wp) :: theta_seed = theta0, rain_pulse = 1.5e-4_wp
   real(wp) :: pond_peak, theta_peak_col
   !----- RUN 8 (snow): seed mass for the pack (0 = no pack, the default for RUNS 1-7) plus the  !
   !      diagnostics integrate_day fills while it is on. ---------------------------------------!
   real(wp) :: snow_seed = 0.0_wp, snow_swe_end, snow_temp_end
   real(wp) :: shed_seed = 0.0_wp   !< [kg/m2/s] daily turnover shed water handed to the ground (RUN 8b)
   real(wp) :: rain_rate = 0.0_wp   !< [kg/m2/s] all-day liquid precipitation (RUN 9b), 0 = off
   real(wp) :: e_col_snow, e_col_rain, snow_mass_total, fusion_expect
   logical  :: cold_air = .false.   !< RUN 9b: the sub-freezing air of the snowfall runs, without snowfall
   logical  :: snow_physical, snowfall_on = .false.
   !----- RUN 5b: a dark, isothermal, rain-free day, and the first step's reported H with the   !
   !      conductance it ran on (rho*ustar*temp1). ------------------------------------------------!
   logical  :: iso_dark = .false.
   real(wp) :: h_first, gah_first
   real(wp) :: snow_swe_split
   !----- RUN 9 (snowfall conservation): the snowfall rate is a VARIABLE so the sub-freezing air and  !
   !      the frozen-rainfall flux can be toggled INDEPENDENTLY -- the run compares snowfall-on       !
   !      against snowfall-off at the SAME air temperature, so an evaporation difference cannot       !
   !      masquerade as the water this test is looking for. Default is RUN 8's rate. `col_water_end`  !
   !      is the whole-column water store (soil + pond + pack).                                       !
   real(wp) :: snowf_rate = 2.0e-5_wp, col_water_end
   real(wp) :: cw_on, cw_off, cw_split_gain, cw_gain, snowf_total
   integer(ik) :: isch
   character(len=9) :: schnm
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: col_config
   type(column_cohort_t)  :: col_cohort
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: biophys
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budget
   type(fast_workspace_t) :: ws          !< the fast step's storage, reused step to step
   type(meds_time_t)      :: sim_date
   real(wp) :: ct_night, ct_noon, co2_night, co2_noon, tleaf_noon, tleaf_night
   real(wp) :: ustar_noon, dtheta_noon      ! midday turbulence regime (see the check below)
   real(wp) :: ss_min, ss_max, sd_min, sd_max, th_min, th_max, gpp_noon, nee_noon
   real(wp) :: psileaf_noon, psileaf_night, psileaf_single
   real(wp) :: surf_water_peak     !< RUN 6: running max of leaf+wood interception film over the day

   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- Full model config (PFT traits for leaf gas exchange) + geometry. ----------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp

   !----- Cohort SoA (single cohort): a 20 cm tree at 224 stems/ha, which gives an LAI of ~3.    !
   !      On-allometry fixture: born through the canonical path (init_cohort + set_cohort_size)  !
   !      and gathered by the production gather, so every field -- height, leaf area, the wood   !
   !      carbon the wood heat capacity is sized from -- is mutually consistent. ----------------!
   call column_cohort_init(col_cohort, cfg%pft, [1_ik], [20.0_wp], [0.0224_wp])
   !----- The hydraulics take the sapwood area in m2; the cohort keeps it in cm2. Passed through  !
   !      unconverted it would be 10^4 times the stem's own cross-section. -----------------------!
   call check_true('the sapwood area reaches the hydraulics in m2, within the stem',             &
                   col_cohort%sap_area(1) > 0.0_wp .and. col_cohort%sap_area(1) <= pio4*0.2_wp**2, &
                   col_cohort%sap_area(1))

   !----- Static column config: soil column + respiration parameters. ---------------------!
   call geometric_layer_bottoms(2.0_wp, nsl, 3.0_wp, test_bottoms)
   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, test_bottoms, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, col_config%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, col_config%soil_thermal)
   !----- Plant hydraulics: flatten cfg%hydraulics -> hydraulics_params + rhizo + build vuln table. ---!
   call apply_hydraulics_config(cfg%hydraulics, cfg%pft, col_config%hydraulics_table)
   call apply_canopy_film_config(cfg%pft, col_config%leaf_surf_water_max, col_config%wood_surf_water_max)
   call build_leaf_photo_table(cfg, col_config%leaf_photo)
   col_config%integrator = build_integrator_opts(cfg)

   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   !----- A cohort's leaf PAR is its absorbed PAR per unit of its own leaf area: a thin cohort's   !
   !      leaves see the same light as a thick one's when each absorbs in proportion to its leaf   !
   !      area, and a leafless cohort sees none. ---------------------------------------------------!
   call check_true('leaf PAR of a thin cohort (LAI 0.05) is its own', abs(cohort_leaf_par(5.0_wp, 0.05_wp, 4.6_wp) &
                   - 460.0_wp) < 1.0e-9_wp, cohort_leaf_par(5.0_wp, 0.05_wp, 4.6_wp))
   call check_true('leaf PAR does not depend on the cohort''s thickness',                                  &
                   abs(cohort_leaf_par(0.05_wp * 100.0_wp, 0.05_wp, 4.6_wp)                              &
                       - cohort_leaf_par(0.5_wp * 100.0_wp, 0.5_wp, 4.6_wp)) < 1.0e-9_wp, 0.0_wp)
   call check_true('a leafless cohort absorbs nothing and sees no PAR', cohort_leaf_par(0.0_wp, 0.0_wp, 4.6_wp) == 0.0_wp, &
                   cohort_leaf_par(0.0_wp, 0.0_wp, 4.6_wp))

   !=====================================================================================!
   !  RUN 1 -- the full physical-sanity suite. The soil thermal column is coupled to the water  !
   !  column by interior-face liquid-enthalpy advection, unconditionally: it is a conservation  !
   !  requirement, not an accuracy knob. The RUN numbers stay fixed (RUN 2 is an empty slot)    !
   !  because test_column_ark and test_column_rk45 cite them. ----------------------------------!
   !=====================================================================================!
   call integrate_day()
   psileaf_single = psileaf_noon                 ! RUN 1's midday leaf psi (recorded, not compared)

   !----- 1. Conservation: all seven budgets closed every step. ----------------------------!
   call check_true('CAS energy budget closed', budget%cas_energy%n_fail    == 0_ik,                                   &
           real(budget%cas_energy%n_fail, wp))
   call check_true('CAS water budget closed', budget%cas_water%n_fail     == 0_ik, real(budget%cas_water%n_fail, wp))
   call check_true('CAS CO2 budget closed', budget%cas_co2%n_fail       == 0_ik, real(budget%cas_co2%n_fail, wp))
   call check_true('soil thermal budget closed', budget%soil_energy%n_fail   == 0_ik,                                 &
           real(budget%soil_energy%n_fail, wp))
   call check_true('soil water budget closed', budget%soil_water%n_fail    == 0_ik,                                   &
           real(budget%soil_water%n_fail, wp))
   call check_true('WHOLE-COLUMN water budget closes', budget%whole_water%n_fail   == 0_ik,                           &
           real(budget%whole_water%n_fail, wp))
   call check_true('WHOLE-COLUMN energy budget closes', budget%whole_energy%n_fail  == 0_ik,                          &
           real(budget%whole_energy%n_fail, wp))

   !----- 2. Physical sanity. -------------------------------------------------------------!
   call check_true('CAS warmer near solar noon than at night', ct_noon > ct_night, ct_noon - ct_night)
   !----- Guard the TURBULENCE REGIME itself, not just the state it produces. Every other              !
   !      assertion in this file is an invariant (conservation, sign, ordering) and so passes in ANY   !
   !      regime, including one where `theta_atm` sits at its 298.15 K type default instead of the     !
   !      forcing's air temperature. With a 282-294 K forcing against a fixed 298.15 K reference, MO   !
   !      manufactures a permanent STABLE layer; correctly driven, midday is UNSTABLE and `ustar`      !
   !      sits off the floor. The stratification check is the DISCRIMINATING one: it fails (-0.56 K)   !
   !      when the set_aero_env_atm call in integrate_day is removed. The ustar check does not fire at !
   !      this test's wind/LAI; it bounds the decoupled regime a fixed reference produces. ------------!
   call check_true('midday stratification is UNSTABLE (canopy warmer than ref air)', dtheta_noon > 0.0_wp,            &
           dtheta_noon)
   call check_true('midday ustar is off the MO floor (not decoupled)', ustar_noon > 1.5_wp * cfg%aero%ustmin,         &
           ustar_noon)
   call check_true('leaf warms under absorbed shortwave', tleaf_noon > tleaf_night, tleaf_noon - tleaf_night)
   call check_true('soil diurnal swing damped with depth', (ss_max - ss_min) > (sd_max - sd_min),                     &
           (ss_max - ss_min) - (sd_max - sd_min))
   call check_true('soil moisture responds to the rain pulse', th_max - th_min > 1.0e-4_wp, th_max - th_min)
   call check_true('soil moisture stays physical', th_min > 0.05_wp .and. th_max < 0.43_wp, th_max)
   call check_true('daytime GPP is active (real photosynthesis)', gpp_noon > 1.0_wp, gpp_noon)
   call check_true('daytime NEE is net uptake (GPP > respiration)', nee_noon < 0.0_wp, nee_noon)
   call check_true('CAS CO2 lower at midday than at night', co2_noon < co2_night, co2_night - co2_noon)
   call check_true('leaf water potential under tension in daylight', psileaf_noon < 0.0_wp, psileaf_noon)
   call check_true('leaf more tensioned at midday than at night', psileaf_noon < psileaf_night,                       &
           psileaf_noon - psileaf_night)
   !----- With the interior faces advecting, soil temperatures stay bounded at BOTH ends of the   !
   !      column. With the faces lumped this day reaches a 361 K soil surface (297 K connected)   !
   !      while the whole-column ledger still closes: the error is purely vertical                !
   !      (docs/science/numerical_scheme.md §5b, "Every soil-water flow comes from one solve"). --!
   call check_true('soil surface temp stays bounded', ss_min > 270.0_wp .and. ss_max < 330.0_wp, ss_max)
   call check_true('deep soil temp stays bounded', sd_min > 270.0_wp .and. sd_max < 330.0_wp, sd_max)

   call test_report('test_column_dynamics')

   !=====================================================================================!
   !  RUN 2 -- an empty slot, so the later RUN numbers stay the ones other files cite.   !
   !=====================================================================================!

   !=====================================================================================!
   !  RUN 3 -- per-layer root coupling. Roots always take up water layer by layer, so RUN 1   !
   !           already exercises the placement and there is no single-BC side to compare      !
   !           against. What is asserted is the invariant the placement must satisfy: the     !
   !           per-layer sink shares sum to 1, so the COLUMN-TOTAL uptake is independent of   !
   !           how it is distributed vertically. The vertical distribution itself is checked  !
   !           in test_plant_hydraulics (test_root_profile).                                  !
   !=====================================================================================!
   call integrate_day()
   call check_true('PER-LAYER ROOTS: whole-column water still closes', budget%whole_water%n_fail  == 0_ik,            &
           real(budget%whole_water%n_fail, wp))
   call check_true('PER-LAYER ROOTS: whole-column energy still closes', budget%whole_energy%n_fail == 0_ik,           &
           real(budget%whole_energy%n_fail, wp))
   call check_true('PER-LAYER ROOTS: leaf psi still under tension', psileaf_noon < 0.0_wp, psileaf_noon)

   !=====================================================================================!
   !  RUN 4 -- caller-side cohort ORDER: aero_bottom_to_top must respect the wind cascade.  !
   !=====================================================================================!
   call test_aero_order()

   !=====================================================================================!
   !  RUN 5 -- the TOLERANCE plumbing: one tolerance set reaches every sub-solver.       !
   !=====================================================================================!
   call test_tolerance_unification()

   !=====================================================================================!
   !  RUN 5b -- the reference level. The surface layer compares the canopy air with the air at  !
   !  the canopy-air top (zref), both referenced to that level, and the reported sensible heat  !
   !  is the flux the whole-column ledger exchanges. A dark column at the air's own temperature !
   !  has nothing driving a sensible flux, so its first step reports |H| far below              !
   !  g_atm_heat*g*zref under both integrators. Referencing only the air's potential            !
   !  temperature to the ground, +(g/cp)*zref, would report H of about -g_atm_heat*g*zref       !
   !  here (docs/science/canopy_aerodynamics.md §2, "The reference level").                     !
   !=====================================================================================!
   call test_reference_level()
   iso_dark = .true.
   do isch = 1_ik, 2_ik
      if (isch == 1_ik) then
         cfg%time_integrator = INTEG_ARK  ; schnm = 'ISO ARK  '
      else
         cfg%time_integrator = INTEG_RK45 ; schnm = 'ISO RK45 '
      end if
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      call integrate_day()
      call check_true(trim(schnm)//': a dark isothermal column reports |H| < 5% of g_atm_heat*g*zref',               &
              abs(h_first) < 0.05_wp * gah_first * grav * aenv%zref, h_first)
      print '(3a,f8.3,a,f8.3,a)', '   (RUN 5b ', trim(schnm), ' H = ', h_first, ' W/m2; g_atm_heat*g*zref = ',       &
            gah_first * grav * aenv%zref, ' W/m2)'
   end do
   iso_dark = .false.
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)

   !=====================================================================================!
   !  RUN 6 -- opt-in CANOPY-SURFACE WATER: reruns the SAME diurnal-cycle-plus-morning-rain-  !
   !           pulse forcing (RUN 1's day, rain over 03:00-07:00) with                        !
   !           col_config%canopy_water_on = .true. whole_water must close exactly and the     !
   !           canopy must actually have intercepted some of the rain (surf_water_peak > 0),  !
   !           proving the WIRING (interception->film->CAS->ledgers) rather than a silent     !
   !           no-op. The wetted-fraction transpiration-suppression algebra itself is         !
   !           unit-tested in test_surface_energy.f90.                                        !
   !                                                                                          !
   !           whole_energy is checked against a LOOSE, explicit bound rather than n_fail==0: !
   !           the film has no thermal state of its own and is valued at ONE fixed reference  !
   !           (t_film_valuation) at both endpoints. That closes MASS and the LATENT-heat     !
   !           exchange (film_evap*latent_heat_vap balances coh_rnet exactly, the identity    !
   !           test_surface_energy.f90 unit-tests) but leaves the water's SENSIBLE heat at    !
   !           the leaf/wood temperature offset uncounted. A film store with its own thermal  !
   !           state is planned (see docs/ROADMAP.md).                                        !
   !=====================================================================================!
   col_config%canopy_water_on = .true.
   call integrate_day()
   col_config%canopy_water_on = .false.                 ! restore default for any future test added after this
   call check_true('CANOPY WATER: whole-column water still closes', budget%whole_water%n_fail  == 0_ik,               &
           real(budget%whole_water%n_fail, wp))
   call check_true('CANOPY WATER: the morning rain pulse was actually intercepted', surf_water_peak > 0.0_wp,         &
           surf_water_peak)
   call check_true('CANOPY WATER: whole-column energy stays BOUNDED (film sensible heat not counted)',               &
           budget%whole_energy%worst < 5.0e6_wp, budget%whole_energy%worst)
   call test_report('test_column_dynamics')

   !=====================================================================================!
   !  RUN 7 -- SATURATED column: exercises the two water-enthalpy paths RUNS 1-6 never take. !
   !                                                                                          !
   !  A sealed (bedrock) column seeded just below theta_sat under heavy rain must both          !
   !  SATURATE -- the water network then raises the layer's potential and drives the excess   !
   !  to the pond and the next layer within the step -- and OVERFLOW its ponding store,       !
   !  producing surface runoff. Neither happens in RUNS 1-6 (free-draining, theta0 = 0.30,    !
   !  gentle pulse). Both paths move mass and enthalpy as a pair: the pond holds its own      !
   !  enthalpy and runoff leaves at the pond's temperature, so no soil layer is debited for   !
   !  water it never received (docs/science/numerical_scheme.md §5b, "The pond holds heat").  !
   !  whole_energy closing to n_fail == 0 here is the assertion that both are booked with the !
   !  right sign and magnitude: an unpaired term of either kind shows up directly in the      !
   !  residual.                                                                               !
   !                                                                                          !
   !  The interior faces advect liquid enthalpy. Without that, the boundary faces would put   !
   !  liquid enthalpy into layer 1 and take it out of a deeper layer with no path between     !
   !  them, so layer 1 would gain the full infiltration enthalpy (~2.5e6 J/m2/step here)      !
   !  while the deeper layer shed it. The WHOLE-column ledger would still close: the error is !
   !  purely in the vertical distribution, which a column-vs-boundary sum cannot see. The     !
   !  per-layer face check below is the one that can.                                         !
   !=====================================================================================!
   theta_seed = 0.425_wp                          ! just below theta_sat = 0.43
   rain_pulse = 2.0e-3_wp                         ! heavy: far above what a sealed column can absorb
   col_config%soil_water_opts%bottom_bc = SOIL_BC_BEDROCK         ! sealed: the water has nowhere to drain
   call integrate_day()
   call check_true('SATURATED: the column reached theta_sat (saturation path is live)',                               &
           theta_peak_col >= 0.43_wp - 1.0e-12_wp, theta_peak_col)
   call check_true('SATURATED: ponding store filled and overflowed (runoff path is live)',                            &
           pond_peak >= col_config%soil_water_opts%w_pond_max - 1.0e-9_wp, pond_peak)
   call check_true('SATURATED: whole-column water still closes', budget%whole_water%n_fail  == 0_ik,                  &
           real(budget%whole_water%n_fail, wp))
   call check_true('SATURATED: whole-column energy still closes', budget%whole_energy%n_fail == 0_ik,                 &
           real(budget%whole_energy%n_fail, wp))
   !----- The VERTICAL check the two above cannot make: the mass the soil-ENERGY equation was           !
   !      charged for, minus the mass the committed theta actually moved, per interior layer. Both      !
   !      budgets above stay at 0 fails through that error, even when it takes the soil surface to      !
   !      345 K, because a misplacement between LAYERS still sums correctly against the boundary.       !
   !                                                                                                    !
   !      This run is on the ARK path, which commits the water network's theta -- so its faces          !
   !      ARE this column's own and the invariant holds by construction. What is asserted               !
   !      here is that construction: the network's own per-layer face closure, at the saturated        !
   !      corner where its pond and runoff paths are live. The RK45 statement of the same invariant is  !
   !      asserted in test_column_rk45. ----------------------------------------------------------------!
   call check_true('SATURATED: per-layer face closure holds (no borrowed advective face)',                            &
           budget%soil_face_mass%worst < 1.0e-8_wp, budget%soil_face_mass%worst)
   call check_true('SATURATED: the face check actually ran', budget%soil_face_mass%n_check > 0_ik,                    &
           real(budget%soil_face_mass%n_check, wp))
   call check_true('SATURATED: soil surface temp stays physical through saturation + runoff',                         &
           ss_min > 250.0_wp .and. ss_max < 340.0_wp, ss_max)
   call test_report('test_column_dynamics')

   !=====================================================================================!
   !  RUN 8 -- SNOW ON: the snow stage's COUPLING, integrated. test_snow.f90 exercises the    !
   !  kernels standalone; this run covers the stage both integrators share                    !
   !  (advance_snow_stage: accumulate -> surface balance -> meltwater drain ->                !
   !  snowfac-blended ground -> sublimation into the CAS) and the pack's mass and energy      !
   !  terms in the whole-column ledgers.                                                      !
   !                                                                                          !
   !  A snow-free run only proves the snow-off path, which is the half that cannot break. The !
   !  two budget assertions below are the regression net: a dropped swe delta, pack rainfall  !
   !  enthalpy or sublimation vapour shows up directly in the residual as an unpaired term.   !
   !=====================================================================================!
   theta_seed = theta0 ; rain_pulse = 0.0_wp      ! no rain: snowfall is the only water input
   col_config%soil_water_opts%bottom_bc = SOIL_BC_FREE_DRAIN
   !----- seeded deep enough that the pack SURVIVES the day under every scheme. A pack that ends   !
   !      near zero makes a "pack still present" assertion test the last few percent of the melt   !
   !      energy rather than the wiring. A surviving pack lets the schemes be compared on melt     !
   !      TOTAL, which is the quantity that says whether they are driving the same stage. ---------!
   snow_seed    = 60.0_wp
   snowfall_on  = .true.
   call integrate_day()
   snow_swe_split = snow_swe_end
   snow_seed = 0.0_wp ; snowfall_on = .false.
   call check_true('SNOW ON: pack + soil stay physical over a 24 h march', snow_physical, snow_temp_end)
   !----- The pack MELTS over this day rather than growing (75 W/m2 of ground shortwave with no       !
   !      longwave loss in this fixture beats 2e-5 kg/m2/s of snowfall), so assert the melt path is    !
   !      live rather than accumulation. That is the better coverage of the two: melt is what routes a !
   !      PAIRED (mass, enthalpy) transfer into the soil top, which is exactly the seam a coupling     !
   !      error is most likely to drop. It stays strictly positive, so a pack is present for the whole !
   !      march and the snowfac-blended surface is exercised throughout. ------------------------------!
   call check_true('SNOW ON: melt/sublimation drained the pack without exhausting it',                                &
           snow_swe_end < 60.0_wp .and. snow_swe_end > 0.0_wp, snow_swe_end)
   call check_true('SNOW ON: whole-column water closes with a pack', budget%whole_water%n_fail  == 0_ik,              &
           real(budget%whole_water%n_fail, wp))
   call check_true('SNOW ON: whole-column energy closes with a pack', budget%whole_energy%n_fail == 0_ik,             &
           real(budget%whole_energy%n_fail, wp))
   call test_report('test_column_dynamics')

   !----- The SAME scenario under ARK and RK45 explicitly. A scheme that imports the snow kernels  !
   !      but never calls them runs silently snow-free, and these assertions catch that; snow-off  !
   !      bit-identity (which the rest of the suite covers) only proves the OFF path, which is the !
   !      half that cannot break. -----------------------------------------------------------------!
   do isch = 1_ik, 2_ik
      if (isch == 1_ik) then
         cfg%time_integrator = INTEG_ARK  ; schnm = 'SNOW ARK '
         col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      else
         cfg%time_integrator = INTEG_RK45  ; schnm = 'SNOW RK45'
         col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      end if
      !----- re-arm the scenario: the block above switches snow back OFF when it finishes, so         !
      !      without this the ARK/RK45 runs are silently snow-FREE -- which looks like success (their !
      !      ledgers close trivially) and is the exact failure mode these assertions exist to catch.  !
      snow_seed = 60.0_wp ; snowfall_on = .true.
      call integrate_day()
      call check_true(trim(schnm)//': pack + soil stay physical with a pack', snow_physical, snow_temp_end)
      call check_true(trim(schnm)//': melt/sublimation drained the pack without exhausting it',                       &
              snow_swe_end < 60.0_wp .and. snow_swe_end > 0.0_wp, snow_swe_end)
      call check_true(trim(schnm)//': whole-column water closes with a pack', budget%whole_water%n_fail  == 0_ik,     &
              real(budget%whole_water%n_fail, wp))
      call check_true(trim(schnm)//': whole-column energy closes with a pack', budget%whole_energy%n_fail == 0_ik,    &
              real(budget%whole_energy%n_fail, wp))
      !----- and each must AGREE with the reference run above (the default scheme) on the pack, not    !
      !      merely conserve on its own: a shared stage fed different inputs by one caller would       !
      !      still close every ledger independently.                                                   !
      !----- compare the MELT TOTAL, not the remaining mass: melt is what the stage computed, while    !
      !      the remainder is melt subtracted from a common seed and so exaggerates small differences  !
      !      once the pack is nearly gone. 10% allows for the schemes genuinely driving the surface    !
      !      balance along different CAS trajectories, while a wiring error (a stage fed the wrong     !
      !      inputs, or not called) is order-100%. ------------------------------------------------!
      call check_true(trim(schnm)//': melt total agrees with the reference run (one shared stage, same inputs)',      &
              abs((60.0_wp - snow_swe_end) - (60.0_wp - snow_swe_split))               < 0.10_wp * (60.0_wp - snow_swe_split),&
              abs(snow_swe_end - snow_swe_split))
      print '(3a,f7.3,a,es10.3)', '   (RUN 8 ', trim(schnm), ' swe=', snow_swe_end, &
            ' kg/m2  worst water resid=', budget%whole_water%worst
   end do
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg

   !=====================================================================================!
   !  RUN 8b -- SHED WATER UNDER A PACK. The daily leaf/root-turnover shed water is routed into   !
   !  the pond together with meltwater at the melt temperature whenever a pack exists, and both   !
   !  whole-column ledgers must book that ground-inflow enthalpy under a pack too. Left unbooked  !
   !  it is a per-patch CONSTANT residual of shed*dt*u_liq(t_melt) on every winter step (~1.8     !
   !  J/m2 at 150 s in the Ithaca example, +6.75 MJ/m2 over a 50-year spin-up), small enough to   !
   !  pass the tolerance, so closure is asked at round-off, on both schemes. =====================!
   do isch = 1_ik, 2_ik
      if (isch == 1_ik) then
         cfg%time_integrator = INTEG_ARK  ; schnm = 'SHED ARK '
         col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      else
         cfg%time_integrator = INTEG_RK45  ; schnm = 'SHED RK45'
         col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      end if
      snow_seed = 60.0_wp ; snowfall_on = .false. ; shed_seed = 2.0e-8_wp
      call integrate_day()
      call check_true(trim(schnm)//': whole-column water closes with shed water under a pack',                        &
              budget%whole_water%n_fail == 0_ik, real(budget%whole_water%n_fail, wp))
      call check_true(trim(schnm)//': whole-column ENERGY closes to round-off with shed water under a pack',          &
              budget%whole_energy%worst < 1.0e-3_wp, budget%whole_energy%worst)
      print '(3a,es10.3,a,f7.3,a)', '   (RUN 8b ', trim(schnm), ' worst energy resid=', budget%whole_energy%worst, &
            ' J/m2  swe=', snow_swe_end, ' kg/m2)'
   end do
   snow_seed = 0.0_wp ; shed_seed = 0.0_wp
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg

   !=====================================================================================!
   !  RUN 9 -- SNOWFALL IS CONSERVED, ON EVERY INTEGRATOR.                                      !
   !                                                                                          !
   !  Snowfall is a boundary water input like rain, and every kg of it must end up somewhere in  !
   !  the column: the pack, the soil, or the pond. A frozen pre-pass that routes only                !
   !  `forc%rainfall`, never `forc%snowfall`, loses `snowfall*dt` every sub-freezing step while      !
   !  the whole-column ledger counts it as an input.                                                 !
   !                                                                                          !
   !  RUN 8 cannot catch that. It asserts on the PACK -- melt totals and pack mass -- which stay     !
   !  right even when the routing of what LEAVES the pack is wrong. This run asserts on the whole    !
   !  column instead: soil liquid + pond + pack.                                                     !
   !                                                                                          !
   !  THE ASSERTION IS A BOUNDARY-INPUT IDENTITY, NOT A LEDGER RESIDUAL, and that is deliberate.     !
   !  A per-step `snowfall*dt` sits below the whole-column closure tolerance even while it integrates   !
   !  into a large seasonal error -- a forced month that drops it runs to completion without         !
   !  tripping the budget hard stop. So this run asks a question the ledger cannot: integrate the    !
   !  SAME day twice, snowfall on and off, holding the air temperature (and hence melt, sublimation  !
   !  and drainage) common, and require the column to gain what fell. A scheme that discards it      !
   !  gains ~0 -- an order-100% failure with no tolerance to tune. A pack is seeded in BOTH runs so   !
   !  that sublimation, a genuine boundary OUTPUT, is common to the pair and cancels in the           !
   !  difference. The first pass (the default scheme) is the reference the others are compared to.   !
   !=====================================================================================!
   snow_seed = 60.0_wp ; rain_pulse = 0.0_wp ; theta_seed = theta0
   col_config%soil_water_opts%bottom_bc = SOIL_BC_FREE_DRAIN
   snowf_total = 2.0e-5_wp * real(nstep, wp) * dt_fast     ! [kg/m2] the day's frozen-rainfall input
   do isch = 1_ik, 3_ik
      select case (isch)
      case (1_ik) ; cfg%time_integrator = INTEG_ARK ; schnm = 'SNOWF ref'
      case (2_ik) ; cfg%time_integrator = INTEG_ARK   ; schnm = 'SNOWF ark'
      case default; cfg%time_integrator = INTEG_RK45   ; schnm = 'SNOWF r45'
      end select
      !----- snowfall ON: cold air + frozen rainfall. -------------------------------------------!
      snowfall_on = .true. ; snowf_rate = 2.0e-5_wp
      call integrate_day()
      cw_on = col_water_end
      call check_true(trim(schnm)//': whole-column water closes', budget%whole_water%n_fail  == 0_ik,                 &
              real(budget%whole_water%n_fail, wp))
      call check_true(trim(schnm)//': whole-column energy closes', budget%whole_energy%n_fail == 0_ik,                &
              real(budget%whole_energy%n_fail, wp))
      !----- snowfall OFF: the SAME cold air and the SAME seeded pack, so melt, sublimation and    !
      !      drainage are common to the pair and cancel in the difference. What is left is the      !
      !      frozen-rainfall input alone. ---------------------------------------------------------!
      snowf_rate = 0.0_wp
      call integrate_day()
      cw_off = col_water_end
      snowf_rate = 2.0e-5_wp ; snowfall_on = .false.
      cw_gain = cw_on - cw_off
      if (isch == 1_ik) cw_split_gain = cw_gain
      !----- 10% allows for the second-order response of melt and drainage to the slightly deeper   !
      !      pack; the defect this guards against is order-100% (the water never arrives). ---------!
      call check_true(trim(schnm)//': snowfall is conserved into the column',                                         &
              abs(cw_gain - snowf_total) < 0.10_wp * snowf_total, cw_gain)
      call check_true(trim(schnm)//': snowfall gain agrees with the reference run',                                   &
              abs(cw_gain - cw_split_gain) < 0.10_wp * snowf_total, cw_gain - cw_split_gain)
      print '(3a,f8.4,a,f8.4,a)', '   (RUN 9 ', trim(schnm), ' column water gain=', &
            cw_gain, ' kg/m2  expected=', snowf_total, ' kg/m2)'
   end do
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   snow_seed = 0.0_wp

   !=====================================================================================!
   !  RUN 9b -- SUB-THRESHOLD SNOWFALL ONTO BARE GROUND ARRIVES AS ICE.                              !
   !  Snowfall too light to start a pack (snowfall*dt < min_new_snow_mass) is routed to the ground   !
   !  with the rain, as ice at the air temperature. Valued as LIQUID at the canopy-air temperature   !
   !  the column would receive the fusion enthalpy of that snow from nowhere -- ledger-consistent    !
   !  (the boundary term uses the same number), physically wrong. Integrate the SAME cold day        !
   !  twice, once with the mass as sub-threshold snow and once as rain; the soil + pond must end     !
   !  the snow day poorer by ~L_f per kg (melting the snow costs soil heat). Liquid valuation gives  !
   !  ~0 difference.                                                                                 !
   !=====================================================================================!
   snow_seed = 0.0_wp ; rain_pulse = 0.0_wp ; theta_seed = theta0
   snowfall_on = .true. ; snowf_rate = 5.0e-6_wp          ! 7.5e-4 kg/m2 per 150 s step < min_new_snow_mass (1e-3)
   call integrate_day()
   e_col_snow = sum(biophys%soil_e%soil_energy(1:nsl) * col_config%soil%dz(1:nsl)) + biophys%soil_w%w_surface_enth
   call check_true('RUN 9b: sub-threshold snowfall never started a pack', biophys%snow%nlayer == 0_ik,                &
           real(biophys%snow%nlayer, wp))
   call check_true('RUN 9b snow: whole-column energy closes to round-off', budget%whole_energy%worst < 1.0e-3_wp,     &
           budget%whole_energy%worst)
   snow_mass_total = snowf_rate * real(nstep, wp) * dt_fast
   snowfall_on = .false. ; cold_air = .true. ; rain_rate = snowf_rate
   call integrate_day()
   e_col_rain = sum(biophys%soil_e%soil_energy(1:nsl) * col_config%soil%dz(1:nsl)) + biophys%soil_w%w_surface_enth
   cold_air = .false. ; rain_rate = 0.0_wp ; snowf_rate = 2.0e-5_wp
   fusion_expect = latent_heat_fusion * snow_mass_total
   !----- ~0.56 L_f on this fixture: the colder, wetter surface loses less to the air over the day and    !
   !      recovers part of the melt cost; the defect this guards against is a ratio of exactly 0. -------!
   call check_true('RUN 9b: snow day leaves the soil+pond poorer than the rain day by ~L_f '//              &
                   'per kg of snow (melt paid by the ground)',                                             &
           e_col_rain - e_col_snow > 0.3_wp * fusion_expect .and. e_col_rain - e_col_snow < 1.5_wp * fusion_expect,   &
           (e_col_rain - e_col_snow) / fusion_expect)
   print '(a,es10.3,a,es10.3,a)', '   (RUN 9b soil+pond energy, rain day - snow day = ', e_col_rain - e_col_snow,      &
         ' J/m2  vs L_f*snow = ', fusion_expect, ' J/m2)'

   call test_report('test_column_dynamics')

contains

   !----- One 24 h diurnal integration from a freshly seeded column state. Fills the host    !
   !      diagnostic variables (min/max ranges, noon/night captures) + budget. ----------------!
   subroutine integrate_day()
      real(wp)    :: t_sec, cosz, t_air, h_step
      integer(ik) :: istep, k


      !----- (Re)seed the prognostic column state. --------------------------------------!
      if (allocated(biophys%leaf_temp)) deallocate(biophys%leaf_temp)
      call alloc_patch_biophys(biophys, n, t0, 0.008_wp, 400.0_wp, t0)
      !----- alloc_patch_biophys seeds leaf_water_mass/wood_water_mass at a scratch 0 (the real     !
      !      seed is reconcile_tissue_water_capacity, run by fast_dynamics, which this driver-level !
      !      test bypasses) -- seed the same water_content(PSI_INIT,...) a freshly-created          !
      !      cohort gets there, or psi_from_water_content would diagnose an unphysical psi from an   !
      !      empty pool. -------------------------------------------------------------------------!
      biophys%leaf_water_mass(1:n) = water_content(PSI_INIT, col_config%hydraulics_table%pft(1)%leaf_curve, &
            col_cohort%bleaf(1:n))
      biophys%wood_water_mass(1:n) = water_content(PSI_INIT, col_config%hydraulics_table%pft(1)%wood_curve, &
            col_cohort%bsap(1:n) + col_cohort%broot(1:n))
      budget = column_budget_t()
      biophys%shed_water_rate = shed_seed
      !----- RUN 8: seed a snow pack when asked. snow_seed = 0 (RUNS 1-7) leaves nlayer = 0, which is  !
      !      exactly the no-pack state alloc_patch_biophys already produces, so those runs are          !
      !      untouched. rho_snow = 250 kg/m3 matches the model's own initial-pack seeding. ------------!
      if (snow_seed > 0.0_wp) then
         biophys%snow%swe(1)         = snow_seed
         biophys%snow%snow_energy(1) = temp_to_internal_energy(0.0_wp, snow_seed, 270.0_wp, 0.0_wp)
         biophys%snow%snow_depth(1)  = snow_seed / 250.0_wp
         biophys%snow%nlayer         = 1_ik
      end if
      snow_physical = .true.
      biophys%soil_w%theta(1:nsl) = theta_seed
      do k = 1_ik, nsl
         biophys%soil_e%soil_energy(k) = temp_to_internal_energy(col_config%soil_thermal%soil_dry_heat_capacity(k),  &
                                     theta_seed * rho_h2o, t0, 1.0_wp)
         biophys%soil_e%soil_temp(k)   = t0
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
         forc%abs_par  = forc%abs_sw                                 ! fixture: absorbed PAR == absorbed SW (par_per_w=2.1)
         forc%abs_lw   = 0.0_wp
         forc%abs_sw_ground = 75.0_wp * cosz
         forc%abs_lw_ground = 0.0_wp
         forc%rainfall   = 0.0_wp
         if (istep >= 3_ik * steps_per_hour .and. istep <= 7_ik * steps_per_hour) &
            forc%rainfall = rain_pulse                   ! morning rain pulse, 03:00-07:00
         if (rain_rate > 0.0_wp) forc%rainfall = rain_rate                        ! RUN 9b: all-day rain
         !----- RUNS 8-9b: steady light snowfall (tests the accumulate path and the pack's          !
         !      rainfall-enthalpy boundary term). 0 for every other run. ---------------------------!
         forc%snowfall = 0.0_wp
         if (snowfall_on) forc%snowfall = snowf_rate
         if (snowfall_on .or. cold_air) t_air = 268.0_wp + 3.0_wp * (cosz - 0.3_wp)   ! sub-freezing: the pack must survive
         if (iso_dark) then     ! RUN 5b: dark, rain-free, the air at the column's own temperature
            t_air = t0 ; forc%abs_sw = 0.0_wp ; forc%abs_par = 0.0_wp ; forc%abs_sw_ground = 0.0_wp
            forc%rainfall = 0.0_wp
         end if
         forc%air_temp         = t_air                                    ! values frozen precipitation as ice at this T
         forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
         forc%shv_atm      = 0.008_wp
         forc%co2_atm      = 400.0_wp
         !----- `theta_atm` is the reference MO measures the canopy against, and it has a            !
         !      plausible 298.15 K type DEFAULT -- setting only forc%air_temp leaves the stability   !
         !      solve comparing the canopy to a fixed 298.15 K, which inverts the sign of            !
         !      stratification and can pin `ustar` on its floor. set_aero_env_atm gives MO the       !
         !      forcing's air temperature, through the SAME routine fill_aenv uses. -----------------!
         call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)

         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, &
                               h_flux=h_step)
         if (iso_dark) then
            h_first = h_step ; gah_first = aenv%rho_air * aero%ustar * aero%temp1
            exit
         end if

         ss_min = min(ss_min, biophys%soil_e%soil_temp(1))   ; ss_max = max(ss_max, biophys%soil_e%soil_temp(1))
         sd_min = min(sd_min, biophys%soil_e%soil_temp(nsl)) ; sd_max = max(sd_max, biophys%soil_e%soil_temp(nsl))
         th_min = min(th_min, biophys%soil_w%theta(1))       ; th_max = max(th_max, biophys%soil_w%theta(1))
         surf_water_peak = max(surf_water_peak, biophys%leaf_surf_water(1) + biophys%wood_surf_water(1))
         pond_peak       = max(pond_peak, biophys%soil_w%w_surface)
         theta_peak_col  = max(theta_peak_col, maxval(biophys%soil_w%theta(1:nsl)))
         if (biophys%snow%nlayer >= 1_ik) then
            snow_physical = snow_physical .and. biophys%snow%swe(1) >= 0.0_wp                            &
                            .and. biophys%snow%snow_temp(1) > 200.0_wp                                   &
                            .and. biophys%snow%snow_temp(1) < 290.0_wp
            do k = 1_ik, nsl
               snow_physical = snow_physical .and. biophys%soil_e%soil_temp(k) > 230.0_wp                &
                               .and. biophys%soil_e%soil_temp(k) < 340.0_wp
            end do
         end if
         if (istep == 13_ik * steps_per_hour + steps_per_hour / 2_ik) then     ! 13:30
            ct_noon = biophys%cas%can_temp ; tleaf_noon = biophys%leaf_temp(1) ; co2_noon = biophys%cas%can_co2
            ustar_noon = aero%ustar ; dtheta_noon = biophys%cas%can_temp - aenv%theta_atm
            gpp_noon = budget%gpp_last ; nee_noon = budget%nee_last
            !----- psi is not persisted state: diagnose it from the persisted leaf_water_mass      !
            !      (docs/science/plant_hydraulics.md §2, eqs. (4)-(5)). ---------------------------!
            psileaf_noon = psi_from_water_content(biophys%leaf_water_mass(1), col_config%hydraulics_table%pft(1)%leaf_curve, &
                  col_cohort%bleaf(1))
         end if
         if (istep == steps_per_hour / 2_ik) then                             ! 00:30
            ct_night = biophys%cas%can_temp ; tleaf_night = biophys%leaf_temp(1) ; co2_night = biophys%cas%can_co2
            psileaf_night = psi_from_water_content(biophys%leaf_water_mass(1), &
                  col_config%hydraulics_table%pft(1)%leaf_curve, col_cohort%bleaf(1))
         end if
      end do
      snow_swe_end = biophys%snow%swe(1) ; snow_temp_end = biophys%snow%snow_temp(1)
      !----- whole-column WATER store at the end of the march (soil layers + pond + pack), for        !
      !      RUN 9's boundary-input accounting. Not a budget -- just the store the arriving water     !
      !      must land in if it lands anywhere at all. -----------------------------------------------!
      col_water_end = biophys%soil_w%w_surface + biophys%snow%swe(1)
      do k = 1_ik, nsl
         col_water_end = col_water_end + biophys%soil_w%theta(k) * col_config%soil%dz(k) * rho_h2o
      end do
   end subroutine integrate_day


   !----- Guard the caller-side cohort order. aero_bottom_to_top feeds canopy_aerodynamics the        !
   !      BOTTOM(1)->TOP(n) order it contracts for, starting from a height-DESCENDING column buffer  !
   !      (index 1 = tallest). So the TALL cohort (gather index 1 = canopy top) must come back with   !
   !      the top-of-canopy wind (=> higher boundary-layer gb); the short understory gets the         !
   !      attenuated wind. Passing the gather order straight into the kernel inverts this and must    !
   !      fail here.                                                                                  !
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
      call aero_bottom_to_top(col_config%aero, e2, g2, 2_ik, c2%height, c2%lai, c2%crown, c2%leaf_width,    &
                              c2%branch_diam, lt, lt, a2)
      call check_true('aero order: tall cohort (gather idx1=top) gets more wind', a2%wind(1) > a2%wind(2),            &
              a2%wind(1) - a2%wind(2))
      call check_true('aero order: tall cohort gets higher leaf gb', a2%leaf_gbw(1) > a2%leaf_gbw(2),                 &
              a2%leaf_gbw(1) - a2%leaf_gbw(2))
      !----- The WOOD boundary layer follows the wood's own temperature (free convection is driven by  !
      !      the surface-air difference), not the leaves'. Guards against passing one temperature      !
      !      array for both, which lets a warm leaf raise the wood's conductance. ---------------------!
      block
         type(aero_out_t) :: a_ref, a_warm_leaf, a_warm_wood
         real(wp) :: tair(2), twarm(2)
         tair = e2%can_temp ; twarm = e2%can_temp + 8.0_wp
         call alloc_aero_out(a_ref, 2_ik) ; call alloc_aero_out(a_warm_leaf, 2_ik) ; call alloc_aero_out(a_warm_wood, 2_ik)
         call aero_bottom_to_top(col_config%aero, e2, g2, 2_ik, c2%height, c2%lai, c2%crown, c2%leaf_width, &
                                 c2%branch_diam, tair, tair, a_ref)
         call aero_bottom_to_top(col_config%aero, e2, g2, 2_ik, c2%height, c2%lai, c2%crown, c2%leaf_width, &
                                 c2%branch_diam, twarm, tair, a_warm_leaf)
         call aero_bottom_to_top(col_config%aero, e2, g2, 2_ik, c2%height, c2%lai, c2%crown, c2%leaf_width, &
                                 c2%branch_diam, tair, twarm, a_warm_wood)
         call check_true('wood gb ignores the leaf temperature', all(a_warm_leaf%wood_gbh == a_ref%wood_gbh), &
                 maxval(abs(a_warm_leaf%wood_gbh - a_ref%wood_gbh)))
         call check_true('warm leaves raise the leaf gb (free convection)', all(a_warm_leaf%leaf_gbh > a_ref%leaf_gbh), &
                 minval(a_warm_leaf%leaf_gbh - a_ref%leaf_gbh))
         call check_true('warm wood raises the wood gb (free convection)', all(a_warm_wood%wood_gbh > a_ref%wood_gbh), &
                 minval(a_warm_wood%wood_gbh - a_ref%wood_gbh))
      end block
   end subroutine test_aero_order

   !----- The reference level (RUN 5b). Canopy air at the temperature and humidity of the air at   !
   !      the canopy-air top is a neutral surface layer with no temperature scale, whatever zref is.  !
   !      Referencing only the air to the ground would solve it as stable, with                       !
   !      tstar = temp1*(g/cp)*zref (docs/science/canopy_aerodynamics.md §2, "The reference level").  !
   subroutine test_reference_level()
      type(column_cohort_t) :: c3
      type(aero_out_t)      :: a3
      type(aero_env_t)      :: e3
      type(aero_geom_t)     :: g3
      real(wp)              :: lt(1)
      call alloc_column_cohort(c3, 1_ik)
      c3%height = [38.0_wp] ; c3%lai = [5.0_wp] ; c3%crown = [0.9_wp]
      c3%leaf_width = [0.04_wp] ; c3%branch_diam = [0.02_wp]
      g3%veg_height = 38.0_wp ; g3%opencan_frac = 0.0_wp ; g3%snowfac = 0.0_wp
      e3%u_ref = 3.0_wp ; e3%zref = 43.0_wp
      call set_aero_env_atm(e3, 300.0_wp, 0.015_wp, 400.0_wp)
      call set_aero_env_canopy(e3, 300.0_wp, 0.015_wp, 400.0_wp, 300.0_wp)
      lt = 300.0_wp
      call alloc_aero_out(a3, 1_ik)
      call aero_bottom_to_top(col_config%aero, e3, g3, 1_ik, c3%height, c3%lai, c3%crown, c3%leaf_width,    &
                              c3%branch_diam, lt, lt, a3)
      call check_true('reference level: equal canopy and reference air has no temperature scale',               &
              a3%tstar == 0.0_wp, a3%tstar)
      call check_true('reference level: equal canopy and reference air is neutral', a3%zeta == 0.0_wp, a3%zeta)
   end subroutine test_reference_level

   !----- One tolerance set for every adaptive error estimate in the fast loop (the ARK march plus       !
   !      the nested sub-solvers). build_tol_set is the single authority; build_fast_context pushes      !
   !      its groups down into the sub-solver opts, so no solver carries a tolerance config cannot       !
   !      reach.                                                                                         !
   !      This is a PLUMBING test on purpose: a whole-model A/B cannot see it (a tolerance change that   !
   !      does not flip a substep decision is silently byte-identical), so the wiring is asserted here   !
   !      directly. It also pins the two accuracy dials, including the atol_scale companion -- rtol_all  !
   !      alone saturates once atol dominates the atol + rtol*|y| denominator.                           !
   subroutine test_tolerance_unification()
      type(meds_config_t)   :: c
      type(fast_context_t)  :: fx
      type(tol_set_t)       :: t
      c = build_test_config()
      !----- (a) DEFAULT: each group is seeded from the config key that governs it. -------------------!
      t = build_tol_set(c)
      call check_true('tol: theta group seeded from [soil].rtol', t%rtol(GRP_THETA)  == c%soil%rtol,                  &
              t%rtol(GRP_THETA))
      call check_true('tol: theta group seeded from [soil].atol', t%atol(GRP_THETA)  == c%soil%atol,                  &
              t%atol(GRP_THETA))
      call check_true('tol: ARK groups seeded from ark_rtol', t%rtol(GRP_ENTH)   == c%ark_rtol, t%rtol(GRP_ENTH))
      !----- (b) MASTER DIALS: rtol_all overrides every group; atol_scale multiplies every atol. ------!
      c%rtol_all   = 1.0e-7_wp
      c%atol_scale = 1.0e-2_wp
      t = build_tol_set(c)
      call check_true('tol: rtol_all overrides ALL groups', all(t%rtol == 1.0e-7_wp), maxval(abs(t%rtol - 1.0e-7_wp)))
      call check_true('tol: atol_scale scales atol', t%atol(GRP_THETA) == c%soil%atol * 1.0e-2_wp, t%atol(GRP_THETA))
      !----- (c) PUSH-DOWN: the dials must reach the nested sub-solver, not just the ARK march: the     !
      !      soil-plant water network takes the GRP_THETA pair for soil and plant stores alike. -------!
      call build_fast_context(c, fx)
      call check_true('tol: dial reaches the soil-WATER sub-solver',                                                  &
              fx%col_config%soil_water_opts%rtol   == 1.0e-7_wp, fx%col_config%soil_water_opts%rtol)
      !----- [energy] has NO integrator tolerance to push into: the soil-thermal solve is one           !
      !      unconditionally-stable BE step, and the accuracy of the integrated soil energy belongs to   !
      !      GRP_SE. What the dial must still reach is [energy].atol, the BUDGET-CLOSURE threshold,      !
      !      which atol_scale scales directly, outside the tolerance groups. This assertion catches a    !
      !      dropped path.                                                                               !
      call check_true('tol: atol_scale reaches the soil-energy BUDGET threshold',                                     &
              fx%col_config%energy%atol == c%energy%atol * 1.0e-2_wp, fx%col_config%energy%atol)
   end subroutine test_tolerance_unification

end program test_column_dynamics
