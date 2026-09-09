!==========================================================================================!
! meds_fast_frozen -- build the FROZEN work record one fast step marches against.            !
!                                                                                          !
! `column_frozen_t` is everything the march treats as constant across a dt_fast: the canopy-  !
! air boundary, the tissue coefficients, the film capacities, the ground boundary, the soil    !
! hydrology and root zone, the plant-water coefficients and the column parameters. PR #124     !
! decomposed it into those nine content-named pieces; this is the routine that fills them.     !
!                                                                                          !
! It lives apart from `meds_fast_ark` because it is not ARK's: RK45 builds the same record      !
! from the same pre-pass, and importing it from the ARK module made the two schemes look        !
! related when they are only siblings. "Frozen" is a statement about the record's LIFETIME      !
! within one step, made by the container and by `intent(in)` -- not a name any field carries.   !
!==========================================================================================!
module meds_fast_frozen
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air, rho_h2o, pi, tsupercool_liq, grav_head, cp_liq, t_3ple
   use meds_plant_hydraulics, only : rhizosphere_cond, solve_plant_water_batch
   use meds_site_diag_types, only : CD_PSI_WOOD, CD_PLC, CD_SAPFLOW, CD_ROOT_UPTAKE
   use meds_hydr_lib, only : soil_hydr_cond_from_theta, soil_psi_from_theta, psi_from_water_content
   use meds_config, only : meds_config_t, CTRL_L2_STRICT
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_soil_types, only : chydro_forcing_t, chydro_flux_t, snow_melt_t
   use meds_fast_types, only : patch_biophys_t
   use meds_column_state_types, only : soil_column_t
   use meds_fast_time_derivs, only : surface_derivs
   use meds_numerics, only : weighted_mean
   use meds_column_state_ops, only : clamp_soil_energy, soil_water_store, soil_energy_store, plant_water_store, &
                                     canopy_film_store, deposit_condensate, clamp_canopy_film, unpack_column_state, &
                                     diagnose_soil_temps, assemble_soil_energy_forcing, apply_process_mask
   use meds_fast_prepass, only : column_prepass
   use meds_fast_types, only : column_config_t, column_cohort_t, column_forcing_t, column_budget_t, column_state_t, &
                               column_frozen_t, surface_state_t, cas_boundary_t, surface_tend_t, stage_bflux_t, &
                               column_bflux_t, error_control_t, column_tend_t, mask_is_full
   use meds_plant_biophysics, only : sensible_heat_coeff, leaf_film_coeff, intercept_canopy_layer
   use meds_soil_water, only : advance_soil_water_column
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction
   use meds_plant_types, only : N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_therm_lib, only : internal_energy_liquid, internal_energy_ice, temp_of_liquid_enthalpy
   use meds_soil_types, only : snow_env_t, snow_flux_t
   use meds_biophysics_opts, only : snow_params_t
   use meds_column_state_types, only : snow_column_t
   use meds_ground_biophysics, only : snow_energy_step
   use meds_fast_types, only : snow_stage_t
   implicit none
   private

   public :: build_column_frozen


   !----- Per-cohort plant-hydraulics sub-step count above which the solve is judged pathological     !
   !      rather than merely stiff (issue #104). Measured band: 1.0-1.2 ordinary, ~136 collapsed. ----!


   !=========================================================================================!
   ! TISSUE HEAT STORE -- ACTIVATION SWITCH. 0 = zero-inertia tissue; 1 = the store live.          !
   ! Currently 1: the store is ON.                                                              !
   !                                                                                          !
   ! Everything the store needs is BUILT AND VERIFIED: the exact exponential relaxation, real WAI + !
   ! sapwood allometry, the dry-wood/sapwood-water capacity split, and -- the hard part -- EXACT    !
   ! conservation on both schemes via the b-weighted tissue-temperature time integrals in           !
   ! column_bflux_t. With this scale at 0 every path reproduces the pre-store answers bit for bit,  !
   ! which is the property the whole design was built around ("diagnostic is the store_hcap_per_dt -> 0 limit !
   ! of one formula, not a separate mode") and it is checked by the full suite passing at 0.        !
   !                                                                                          !
   ! WHY IT IS NOT ON YET. Turning it on flips four PHYSICS assertions on ARK that had only ever    !
   ! been exercised on the retired split path: daytime NEE goes net-release, and the leaf water     !
   ! potential comes out POSITIVE (+0.72 MPa) instead of under tension, i.e. transpiration is being !
   ! suppressed and leaf water accumulates. That is not a conservation failure -- every budget still !
   ! closes -- but it is unexplained, and the leaf store is far too small to explain it directly     !
   ! (cap_leaf = 2205 J/m2/K against h_coeff = 100 W/m2/K, so tau <~ 22 s and w_avg <~ 0.012).       !
   ! Turning the LEAF store on alone is measurably WORSE than both together, which is non-monotonic  !
   ! in the store size and therefore points at something other than the store's own inertia.         !
   !                                                                                          !
   ! Deliberately a source constant, not a TOML knob: this is an unfinished feature, not a supported !
   ! configuration choice, and it must not look like one. Flip to 1.0 to resume the investigation.   !
   !=========================================================================================!


   !----- Prognostic-wood constants (retained; the split twin that these mirrored is retired). --------!
   !      build the same store from the same biomass. ----------------------------------------------!
   real(wp), parameter :: TISSUE_STORE_SCALE  = 1.0_wp  !< tissue heat store ON (see the banner above)
   !----- Per-cohort plant-hydraulics sub-step count above which the solve is judged pathological   !
   !      rather than merely stiff (issue #104). Measured band: 1.0-1.2 ordinary, ~136 collapsed. --!
   integer(ik), parameter :: HYDRO_NSUB_THRASH  = 16_ik
   real(wp), parameter :: C2B_WOOD            = 2.0_wp  !< carbon -> biomass (carbon fraction 0.5)
   real(wp), parameter :: WOOD_MOIST_FRAC_ARK = 1.0_wp  !< [kg water/kg dry] fresh-sapwood moisture (MVP)

contains

   !----- Build the frozen ARK inputs: the shared column_prepass (meds_fast_prepass: leaf gas         !
   !      respiration / CAS caps / aero, bit-identical to the split) + this integrator's own per-cohort !
   !      geometry/radiation/wood packing + the frozen hydrology BCs, into a column_frozen_t; also      !
   !      packs the prognostic state into a column_state_t.                                             !
   subroutine build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, n, nsl, &
                                  frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      !----- intent(inout) since C4: the shared snow stage ADVANCES the pack (biophys%snow) and hands its  !
      !      melt enthalpy to biophys%soil_e%soil_energy(1) as a paired transfer, exactly as on the split   !
      !      path. Every other use of biophys here is still read-only. ------------------------------------!
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: frozen
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture

      type(snow_stage_t)     :: snow_st   !< shared pre-column snow stage (C4, issue #76)
      type(chydro_forcing_t) :: hforc ; type(chydro_flux_t) :: hflux
      type(soil_column_t)    :: soil_w_scratch
      type(surface_state_t)  :: y_stage ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, nee_biotic, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
           g_atm_co2
      integer(ik) :: i, k
      !----- Act-1 hydraulics pre-pass scratch (MEDS_ED2_RK45_DESIGN.md sec 1/3/5, P2): mirrors        !
      !      the original pre-pass order exactly (hydraulics BEFORE the soil solve, psi                   !
      !      diagnosed from state^n theta, the plant's own aggregate REQUEST becomes the soil's          !
      !      root-sink forcing, a post-hoc rescale if the soil can't honour it in full). --------------!
      real(wp) :: psi_soil_pre(nsl), psi_scratch(N_HYDRO, n), transp_pp(n)
      real(wp) :: sapflow_b(n), root_uptake_b(n), root_uptake_layer_b(nsl, n)
      real(wp) :: psi_leaf_b(n), psi_wood_b(n), plc_b(n)   !< batch outputs (unused downstream, complete SoA API)
      real(wp) :: rhizo_cond_all(nsl, n), k_theta_layer(nsl), total_uptake_b, scale, share_tot
      real(wp) :: t_up_wl, soil_temp_root, u_liq_soil, u_liq_up
      real(wp) :: sapflow_gnd(n), uptake_gnd(n)
      integer(ik) :: nsub_b(n)
      logical     :: converged_b(n)
      !----- Canopy-SURFACE water pre-pass scratch (MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c). ------------!
      real(wp) :: rain_above, combined_w, pai_i, throughfall_i, drip_i, throughfall_total
      real(wp) :: avail_leaf, avail_wood

      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n),                  &
         frozen%tissue%abs_lw(n), frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      allocate(frozen%tissue%wood_dry_hcap(n), frozen%tissue%wood_wmass(n))
      allocate(frozen%tissue%leaf_dry_hcap(n), frozen%tissue%leaf_wmass(n))
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      !----- ZERO the advective-enthalpy terms AT ALLOCATION. They are only given their real values    !
      !      further down (after the plant-hydraulics batch supplies sapflow/uptake), but the sf0       !
      !      surface_derivs evaluations ABOVE that point already read them as veg_energy_balance's    !
      !      q_extra. Without this they are read UNINITIALISED -- undefined behaviour that injected       !
      !      whatever garbage the freshly-allocated heap happened to hold straight into the leaf energy   !
      !      balance. That is the origin of the 30-yr run's hang: garbage q_extra -> |dt_temp| ~1e64..     !
      !      1e292 -> absurd transp_c -> the hydraulics returned a non-finite root_uptake -> the soil-     !
      !      water adaptive loop spun forever on err = NaN (neither t nor nsub can advance). It also       !
      !      explains why the failure moved between builds and even between runs of the SAME binary: the   !
      !      garbage depends on prior heap contents, not on the physics. 0 is exactly the documented       !
      !      default for q_extra ("absent/0 for every caller but the P2 advective-enthalpy pre-pass"). ----!
      frozen%tissue%qwflux_wl(1:n)  = 0.0_wp
      frozen%tissue%q_wood_net(1:n) = 0.0_wp
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n))
      allocate(frozen%roots%root_share(nsl), frozen%plant%nplant(n), frozen%plant%bleaf(n),                   &
         frozen%plant%bsap(n), frozen%plant%broot(n),            &
               frozen%plant%sap_area(n), frozen%plant%height(n), frozen%plant%leaf_area(n))
      allocate(frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n), frozen%roots%qloss_frozen(n))
      allocate(frozen%roots%psi_soil_pre(nsl), frozen%roots%rhizo_cond(nsl, n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n))
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water(1:n) = biophys%leaf_surf_water(1:n) ; y%wood_surf_water(1:n) = biophys%wood_surf_water(1:n)

      !----- the SHARED pre-pass (meds_fast_prepass%column_prepass): gas exchange / respiration / CAS   !
      !      aero -- writes directly into the frozen struct's h_coeff_leaf/g_transp_leaf arrays. ------------------!
      call column_prepass(cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,                       &
                          tcas, qcas, press, rho, t_ground, frozen%tissue%h_coeff_leaf, frozen%tissue%g_transp_leaf,      &
                          cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, nee_biotic, &
                          gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      !----- SHARED SNOW STAGE (C4, issue #76). AFTER column_prepass: it needs the CAS state           !
      !      (tcas/qcas), air properties (rho/press) and the aerodynamic conductance (aero%ggnet) that !
      !      the pre-pass establishes -- placing it earlier reads those undefined and yields a NaN     !
      !      pack. Split calls it after its own column_prepass for exactly this reason. Still BEFORE   !
      !      the hydrology forcing below, so meltwater reaches infiltration this step and the melt     !
      !      enthalpy is inside the soil column the state^n snapshot takes. No-op without a pack.      !
      call advance_snow_stage(biophys%snow, col_config%snow, max(-col_config%soil%z_node(1), tiny_num),    &
                              forc%abs_sw_ground, forc%abs_lw_ground, forc%snowfall, forc%rainfall, forc%air_temp,  &
                              aero%ggnet, biophys%soil_e%soil_temp(1), dt_fast, tcas, qcas, rho, press,    &
                              snow_st)
      frozen%snow = snow_st                       ! the stage's whole outcome, one record

      !----- Canopy INTERCEPTION (sec 3.4, P2c): frozen ONCE per dt_fast, frozen once per dt_fast:       !
      !      own "2c. CANOPY INTERCEPTION" sweep. ONE combined leaf+wood bucket per cohort, top-to-       !
      !      bottom over col_cohort's native height-DESCENDING gather order (the SAME direction the split path's  !
      !      own i=1..n loop already assumes is top-first), e_canopy=0 (capture/capacity only -- film       !
      !      evaporation is the SEPARATE, per-stage flux surface_derivs computes below from the frozen      !
      !      f_wet_c/g_film_leaf/g_film_w this block also sets). Converts the one-shot bucket update into an    !
      !      EQUIVALENT frozen RATE (matching sapflow_frozen/uptake_frozen's own "one frozen number, no      !
      !      per-stage re-solve" convention): integrating this rate by explicit Euler over dt_fast exactly    !
      !      reproduces the split's one-shot commit. Gated behind canopy_water_on so the untouched default    !
      !      path stays byte-identical (throughfall_total defaults to the FULL rainfall + snowfall sum, matching   !
      !      the pre-P2c hforc%precip_ground line below unconditionally).                                    !
      !                                                                                                      !
      !      PLACEMENT (E-6, MEDS_INTEGRATOR_PARITY.md [RETIRED] sec 3e): this block MUST run after                     !
      !      advance_snow_stage, because it branches on snow_st%exists -- the pack OWNS the surface and        !
      !      liquid-rain interception is mutually exclusive with it, exactly as on the split path              !
      !      (the pack owns the surface, so liquid interception is skipped). It used to sit above            !
      !      column_prepass on the (correct, but insufficient) grounds that it has no aero dependency, and     !
      !      so read snow_st%exists BEFORE its only writer ran. snow_stage_t default-initialises exists to     !
      !      .false., so the read was defined rather than undefined -- but it was ALWAYS .false., i.e. the     !
      !      guard never fired. Nothing between here and the first surface_derivs call reads f_wet_c or        !
      !      intercept_leaf/wood, so moving the block down is otherwise inert.                                !
      !                                                                                                       !
      !      SNOWFALL (E-4): throughfall_total carries forc%rainfall + forc%snowfall, as the split path always      !
      !      has. precip_phase splits the met precipitation into rain and snow WITHOUT consulting snow_on,     !
      !      so with the snow store off (the default) a sub-freezing step still delivers forc%snowfall > 0.       !
      !      Omitting it here dropped that water on the floor -- it never reached hforc%precip_ground, while   !
      !      the whole-column ledger's w_in counted it (column_fast_step_ark / _rk45), so ARK and RK45 leaked  !
      !      exactly snowfall*dt per step and diverged from split all winter (measured: split retained 39% more   !
      !      soil water over a January month; the gap vanished with snow_on = .true.). Under a pack this term  !
      !      is unused -- snow_accumulate has already taken BOTH rainfall and snowfall, so precip_ground becomes    !
      !      meltwater only (the branch below). ---------------------------------------------------------------!
      frozen%film%intercept_leaf = 0.0_wp ; frozen%film%intercept_wood = 0.0_wp
      frozen%film%f_wet_c = 0.0_wp ; frozen%film%g_film_leaf = 0.0_wp ; frozen%film%g_film_w = 0.0_wp
      throughfall_total = forc%rainfall + forc%snowfall
      if (col_config%canopy_water_on .and. .not. snow_st%exists) then
         rain_above = forc%rainfall + forc%snowfall
         do i = 1_ik, n
            pai_i      = col_cohort%lai(i) + col_cohort%wai(i)
            combined_w = biophys%leaf_surf_water(i) + biophys%wood_surf_water(i)
            call intercept_canopy_layer(combined_w, rain_above, col_cohort%lai(i), col_cohort%wai(i), 0.0_wp, dt_fast, &
                                        col_config%soil_water_opts%dewmx, col_config%soil_water_opts%intercept_k, &
                                             col_config%soil_water_opts%intercept_alpha, &
                                        throughfall_i, drip_i, frozen%film%f_wet_c(i))
            if (pai_i > tiny_num) then
               frozen%film%intercept_leaf(i) = (combined_w*col_cohort%lai(i)/pai_i - biophys%leaf_surf_water(i)) / dt_fast
               frozen%film%intercept_wood(i) = (combined_w*col_cohort%wai(i)/pai_i - biophys%wood_surf_water(i)) / dt_fast
            else
               frozen%film%intercept_leaf(i) = -biophys%leaf_surf_water(i) / dt_fast
               frozen%film%intercept_wood(i) = -biophys%wood_surf_water(i) / dt_fast
            end if
            rain_above = throughfall_i   ! cascades to the next (shorter) cohort
         end do
         throughfall_total = rain_above   ! whatever survives the shortest (last) cohort
      end if

      !----- per-cohort geometry + radiation + WOOD frozen inputs the ARK path needs (not shared with   !
      !      the split, which reads col_cohort%/forc% directly instead of packing a frozen struct). -----------!
      do i = 1_ik, n
         frozen%tissue%lai(i)    = col_cohort%lai(i)
         frozen%tissue%abs_sw(i) = forc%abs_sw(i) ; frozen%tissue%abs_lw(i) = forc%abs_lw(i)
         !----- WOOD frozen inputs: real diagnostic values, or ZERO when wood is not diagnostic (so   !
         !      surface_derivs' wood branch is a no-op; prognostic wood is operator-split in P2). -----!
         !----- PROGNOSTIC-WOOD frozen set (Phase 4). Built unconditionally -- it is simply unread when  !
         !      wood is diagnostic -- so the two wood authorities never both feed the CAS: whichever mode !
         !      is active, exactly one of {frozen%surf's diagnostic inputs, frozen%wood_*} is non-zero. -------!
         !----- The two halves of the wood heat capacity take DIFFERENT masses. DRY tissue is ALL   !
         !      the wood (col_cohort%bwood) -- heartwood is dead structure but it still stores sensible    !
         !      heat, and a sapwood fraction defined on the bole under-counts branch wood, which is  !
         !      thin enough to be thermally active throughout. INTERNAL WATER is the sapwood ring    !
         !      ONLY (col_cohort%bsap): heartwood is taken as dry, which is what makes it heartwood. bsap   !
         !      is also the HYDRAULIC capacitance -- same ring, two consumers, consistently. --------!
         frozen%tissue%wood_dry_hcap(i) = max(col_cohort%bwood(i) * col_cohort%nplant(i) * C2B_WOOD             &
                                              * col_config%veg_thermal%c_sapw,                              &
                                    col_config%veg_thermal%veg_hcap_min)
         frozen%tissue%wood_wmass(i)    = col_cohort%bsap(i)  * col_cohort%nplant(i) * C2B_WOOD * WOOD_MOIST_FRAC_ARK
         !----- LEAF capacity, same two-part construction: dry leaf tissue + the INTERNAL (symplast) !
         !      water the hydraulics actually carries, biophys%leaf_water_mass [kg/plant] -> per m2 ground.!
         !      Read from biophys, NOT from y: y%leaf_water_mass is not filled until the very end of this  !
         !      routine, so reading it here would be an uninitialised read -- the same trap the        !
         !      advective-enthalpy zeroing note above records. ---------------------------------------!
         !      The intercepted FILM is deliberately NOT here: it is a separate store with its own     !
         !      phase (Step B), because folding it into a temperature-based capacity would commit to   !
         !      a film that can never freeze. -------------------------------------------------------!
         frozen%tissue%leaf_dry_hcap(i) = max(col_cohort%bleaf(i) * col_cohort%nplant(i) * C2B_WOOD             &
                                              * col_config%veg_thermal%c_leaf,                              &
                                    col_config%veg_thermal%veg_hcap_min)
         frozen%tissue%leaf_wmass(i)    = max(biophys%leaf_water_mass(i), 0.0_wp) * col_cohort%nplant(i)
         !----- The wood's diagnostic (zero-inertia) inputs are now filled UNCONDITIONALLY. They used  !
         !      to be zeroed whenever wood was "prognostic", because the prognostic store was a wholly  !
         !      separate operator-split solve that owned the wood's radiation and sensible flux. With   !
         !      the post-commit store (apply_tissue_store) there is only ONE wood solve: the diagnostic !
         !      balance here, relaxed afterwards by its own heat capacity. Zeroing these would leave    !
         !      apply_tissue_store relaxing towards a wood temperature of tcas with no absorbed         !
         !      radiation, i.e. silently deleting the wood's energy budget. -----------------------------!
         frozen%tissue%wai(i)         = col_cohort%wai(i)
         frozen%tissue%h_coeff_w(i)   = sensible_heat_coeff(pi * col_cohort%wai(i), aero%wood_gbh(i), rho, cp_air)
         frozen%tissue%abs_sw_wood(i) = forc%abs_sw_wood(i)
         frozen%tissue%abs_lw_wood(i) = forc%abs_lw_wood(i)
         !----- TISSUE STORE, frozen for the whole dt_fast. a = cap/dt_fast; the relaxation origin is  !
         !      the START-of-step tissue temperature. Every stage evaluation therefore returns the      !
         !      same dt_fast-averaged flux and dt_fast-endpoint temperature -- the store is an algebraic !
         !      closure, not a tableau DOF. --------------------------------------------------------------!
         frozen%tissue%leaf_hcap_per_dt(i)  = TISSUE_STORE_SCALE                                                   &
                               * (frozen%tissue%leaf_dry_hcap(i) + frozen%tissue%leaf_wmass(i) * cp_liq) / dt_fast
         frozen%tissue%wood_hcap_per_dt(i)  = TISSUE_STORE_SCALE                                                   &
                               * (frozen%tissue%wood_dry_hcap(i) + frozen%tissue%wood_wmass(i) * cp_liq) / dt_fast
         frozen%tissue%t_leaf0(i) = biophys%leaf_temp(i)
         frozen%tissue%t_wood0(i) = biophys%wood_temp(i)
         frozen%plant%nplant(i)   = col_cohort%nplant(i)
         frozen%plant%bleaf(i)    = col_cohort%bleaf(i)
         frozen%plant%bsap(i) = col_cohort%bsap(i)
         frozen%plant%broot(i)    = col_cohort%broot(i)   ; frozen%plant%sap_area(i) = col_cohort%sap_area(i)
         frozen%plant%height(i)   = col_cohort%height(i)  ; frozen%plant%leaf_area(i) = col_cohort%leaf_area(i)
         !----- Canopy-SURFACE water film-evap conductances (sec 3.4, P2c): need aero%leaf_gbw/wood_gbw,   !
         !      so these run HERE (after column_prepass's aero solve above), not in the interception        !
         !      block before it. Gated behind canopy_water_on (not just "harmless via f_wet_c=0") per the    !
         !      P1 nvfortran lesson (this doc's own "P1 implementation notes": gate new ledger/diagnostic     !
         !      arithmetic behind its own flag from the start, rather than relying on it telescoping to a      !
         !      no-op). ------------------------------------------------------------------------------------!
         if (col_config%canopy_water_on) then
            frozen%film%g_film_leaf(i) = leaf_film_coeff(col_config%veg_thermal%effarea_evap, col_cohort%lai(i), aero%leaf_gbw(i))
            frozen%film%g_film_w(i) = leaf_film_coeff(col_config%veg_thermal%effarea_evap, col_cohort%wai(i), aero%wood_gbw(i))
         end if
      end do

      !----- the rest of the frozen surface inputs: CAS caps/conductances from column_prepass + atm     !
      !      state + NEE. ---------------------------------------------------------------------------!
      frozen%tissue%leaf_emiss = col_config%veg_thermal%leaf_emiss
      frozen%cas%cas_mass_capacity = cas_mass_capacity ; frozen%cas%cas_molar_capacity = cas_molar_capacity
      frozen%cas%g_atm_heat  = g_atm_heat  ; frozen%cas%g_atm_vapour  = g_atm_vapour ; frozen%cas%g_atm_co2 = g_atm_co2
      !----- Everything refresh_cas_conductances needs that is NOT the live CAS state. The         !
      !      geometry pair (displacement, roughness) is taken from the pre-pass's OWN aero output   !
      !      rather than recomputed, so a stage re-solve starts from the identical surface the      !
      !      state^n solve used, and the two agree exactly when the state has not moved. ----------!
      frozen%cas%aero_cfg     = col_config%aero
      frozen%cas%mo_u_ref     = aenv%u_ref     ; frozen%cas%mo_zref      = aenv%zref
      frozen%cas%mo_displace  = aero%displace  ; frozen%cas%mo_rough     = aero%rough
      frozen%cas%mo_theta_atm = aenv%theta_atm ; frozen%cas%mo_shv_atm   = aenv%shv_atm
      frozen%cas%mo_rho       = rho
      !----- ...and declare the inputs live, which is what licenses a per-stage re-solve. A bundle  !
      !      that has NOT been through here (a unit-test fixture, the RK4 oracle) leaves this false !
      !      and its supplied g_atm_heat/g_atm_vapour/g_atm_co2 are used verbatim. ----------------------------------------!
      frozen%cas%mo_live      = .true.
      frozen%cas%enthalpy_atm = forc%enthalpy_atm ; frozen%cas%shv_atm = forc%shv_atm ; frozen%cas%co2_atm = forc%co2_atm
      frozen%cas%nee_biotic = nee_biotic
      frozen%ground%abs_sw_ground = forc%abs_sw_ground ; frozen%ground%abs_lw_ground = forc%abs_lw_ground
      frozen%ground%ggnet = aero%ggnet ; frozen%cas%rho = rho ; frozen%cas%press = press

      !----- params + hydraulics BCs. -----------------------------------------------------------!
      frozen%params%soil = col_config%soil ; frozen%params%therm = col_config%soil_thermal
      frozen%params%energy_opts = col_config%energy
      frozen%params%hydro_opts = col_config%soil_water_opts
      frozen%cas%cas_condensation = col_config%integrator%cas_condensation      ! §8g scheme-asymmetry guard
      frozen%hydrology%geothermal = 0.0_wp

      !----- Act 1 (MEDS_ED2_RK45_DESIGN.md sec 1/3/5, P2): plant hydraulics runs BEFORE the soil     !
      !      solve, using psi diagnosed from state^n theta and the FULL transpiration demand -- no     !
      !      supply pre-throttle (the plant's own leaf/wood water MASS storage buffers any step-to-     !
      !      step soil-supply/demand mismatch instead).                                                  !
      !      Hydraulics runs BEFORE the soil solve so the soil sees the realized uptake. ------------------!
      y_stage%cas_enthalpy = biophys%cas%can_enthalpy
      y_stage%cas_shv = biophys%cas%can_shv
      y_stage%cas_co2 = biophys%cas%can_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, t_ground, n, sf0)
      !----- Canopy-SURFACE water (sec 3.4, P2c): rescale the frozen film-evap conductance -- like        !
      !      uptake_frozen's own soil-limiting rescale above -- so a WORST-CASE potential evaporation       !
      !      over the FULL dt_fast (sf0's state^n film_evap, using the FULL unscaled g_film_leaf/w just         !
      !      computed) cannot drain MORE than will be available (current store + this step's frozen           !
      !      interception). Without this, film_evap (a state^n-frozen conductance, oblivious to               !
      !      depletion mid-step) can overdraw the store, driving leaf/wood_surf_water NEGATIVE -- which         !
      !      corrupts the NEXT step's frozen interception rate (intercept_canopy_layer's own internal            !
      !      floor silently "fixes" a negative starting bucket, fabricating mass that was never really            !
      !      lost) rather than genuinely closing the ledger. Re-evaluates sf0 afterward so the SAME             !
      !      call's transp_c (feeding the hydraulics batch solve below) sees the final, scaled conductance.      !
      !      Zero-effect when canopy_water_on is off (film_evap is 0 there, scale stays 1). -------------------!
      if (col_config%canopy_water_on) then
         do i = 1_ik, n
            avail_leaf = biophys%leaf_surf_water(i) + max(0.0_wp, frozen%film%intercept_leaf(i)) * dt_fast
            if (sf0%film_evap_leaf(i) > tiny_num) frozen%film%g_film_leaf(i) = frozen%film%g_film_leaf(i)          &
                 * min(1.0_wp, avail_leaf / (sf0%film_evap_leaf(i)*dt_fast))
            avail_wood = biophys%wood_surf_water(i) + max(0.0_wp, frozen%film%intercept_wood(i)) * dt_fast
            if (sf0%film_evap_wood(i) > tiny_num) frozen%film%g_film_w(i) = frozen%film%g_film_w(i)          &
                 * min(1.0_wp, avail_wood / (sf0%film_evap_wood(i)*dt_fast))
         end do
         call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, t_ground, n, sf0)
      end if
      !----- UNITS: grav_head converts soil_psi_from_theta's METRES of head to the MPa the hydraulics   !
      !      seam expects. -----------------------------------------------------------------------------!
      psi_soil_pre(1:nsl) = grav_head * soil_psi_from_theta(col_config%soil%retention, biophys%soil_w%theta(1:nsl), &
           col_config%soil%theta_sat(1:nsl), col_config%soil%theta_res(1:nsl), col_config%soil%vg_alpha(1:nsl),          &
           col_config%soil%vg_n(1:nsl))
      !----- Per-layer rhizosphere conductance, UNCONDITIONAL since Phase 1 retired multilayer_roots.  !
      !      K(theta) is layer-only, so it is hoisted out of the cohort loop (hot-path work now). ------!
      do k = 1_ik, nsl
         k_theta_layer(k) = soil_hydr_cond_from_theta(col_config%soil%retention, biophys%soil_w%theta(k),        &
              col_config%soil%theta_sat(k), col_config%soil%theta_res(k), col_config%soil%vg_alpha(k),                   &
              col_config%soil%vg_n(k), col_config%soil%ksat(k))
      end do
      do i = 1_ik, n
         do k = 1_ik, nsl
            rhizo_cond_all(k, i) = rhizosphere_cond(rho_h2o*k_theta_layer(k)/grav_head, col_cohort%broot(i),  &
                 col_config%specific_root_area, col_config%soil%root_frac(k), col_config%soil%dz(k), col_cohort%nplant(i))
         end do
      end do
      !----- Hand the kernel's own frozen boundary inputs to the post-stage corrector, so it can       !
      !      re-solve on the SAME Category-0 coefficients this pre-pass used. -----------------------!
      frozen%roots%psi_soil_pre(1:nsl)        = psi_soil_pre(1:nsl)
      frozen%roots%rhizo_cond(1:nsl, 1:n)     = rhizo_cond_all(1:nsl, 1:n)
      frozen%params%hydraulics_params                    = col_config%hydraulics_params
      frozen%params%hydraulics_opts                    = col_config%hydraulics_opts
      psi_scratch(NODE_LEAF, 1:n) = psi_from_water_content(biophys%leaf_water_mass(1:n),                   &
           col_config%hydraulics_params%leaf_pi0, col_config%hydraulics_params%leaf_elastic_mod, &
                col_config%hydraulics_params%leaf_apoplast_frac,      &
           col_config%hydraulics_params%leaf_water_sat, col_cohort%bleaf(1:n))
      psi_scratch(NODE_WOOD, 1:n) = psi_from_water_content(biophys%wood_water_mass(1:n),                   &
           col_config%hydraulics_params%wood_pi0, col_config%hydraulics_params%wood_elastic_mod, &
                col_config%hydraulics_params%wood_apoplast_frac,      &
           col_config%hydraulics_params%wood_water_sat, col_cohort%bsap(1:n) + col_cohort%broot(1:n))
      transp_pp(1:n) = sf0%transp_c(1:n) / max(col_cohort%nplant(1:n), tiny_num)   ! [kg/plant/s] FULL demand
      call solve_plant_water_batch(n, nsl, transp_pp(1:n), col_cohort%bleaf(1:n),                             &
                                   col_cohort%bsap(1:n), col_cohort%broot(1:n), col_cohort%sap_area(1:n), &
                                         col_cohort%height(1:n),   &
                                   col_cohort%leaf_area(1:n),                                                  &
                                   psi_soil_pre(1:nsl), col_config%soil%z_node(1:nsl), rhizo_cond_all(1:nsl, 1:n), &
                                   col_config%hydraulics_params, col_config%hydraulics_opts, dt_fast, psi_scratch(:, 1:n), &
                                   sapflow_b(1:n), root_uptake_b(1:n), root_uptake_layer_b(1:nsl, 1:n),  &
                                   psi_leaf_b(1:n), psi_wood_b(1:n), plc_b(1:n), nsub_b(1:n), converged_b(1:n))
      budget%hydro_nsub    = sum(nsub_b(1:n))            ! section 5.3 work counter (same seam as split)
      budget%hydro_nonconv = count(.not. converged_b(1:n))
      !----- THRASH check (issue #104, plan E1b/E3). solve_plant_water takes 1.0-1.2 sub-steps per     !
      !      cohort in every physically ordinary state and ~136 once a wood store is pinned on its     !
      !      water floor. The threshold sits an order of magnitude above the ordinary band and an       !
      !      order below the collapsed one, so it fires on the pathology and nothing else. Under L2     !
      !      this is fatal -- a run that silently goes 13x slower AND returns floored stores is not a   !
      !      run anyone wants to keep. Under L1/L0 it is counted and reported through the existing      !
      !      work-counter output path (work_hydro_thrash_site). ------------------------------------!
      budget%hydro_thrash = merge(1_ik, 0_ik, budget%hydro_nsub > n * HYDRO_NSUB_THRASH)
      if (budget%hydro_thrash == 1_ik .and. col_config%integrator%error_control%level == CTRL_L2_STRICT) &
         error stop 'column_fast_step_ark: plant-hydraulics sub-stepping is pathological (issue #104) &
                    &-- a tissue store has almost certainly collapsed onto its water floor. Inspect &
                    &wood_water_mass; see docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md sec 5c(v).'
      !----- HR (root efflux) intentionally NOT enabled anywhere in this model -- floor the aggregate  !
      !      (the same floor the hydrology kernel applies). ---------------------------------------------!
      root_uptake_b(1:n) = max(root_uptake_b(1:n), 0.0_wp)
      root_uptake_layer_b(1:nsl, 1:n) = max(root_uptake_layer_b(1:nsl, 1:n), 0.0_wp)
      !----- Per-cohort HYDRAULIC diagnostics, captured from the SAME solve the physics commits.    !
      !      psi_wood / plc were computed on every dt_fast and discarded; plc in particular is the    !
      !      one number that says whether a cohort is losing conductance, which nothing downstream    !
      !      of here could previously see. Taken AFTER the HR floor, so the reported uptake is the    !
      !      realized (non-negative) one the soil sink actually used.  ------------------------------!
      if (present(cdiag)) then
         cdiag(CD_SAPFLOW,     1:n) = sapflow_b(1:n)
         cdiag(CD_ROOT_UPTAKE, 1:n) = root_uptake_b(1:n)
         cdiag(CD_PSI_WOOD,    1:n) = psi_wood_b(1:n)
         cdiag(CD_PLC,         1:n) = plc_b(1:n)
      end if
      total_uptake_b = sum(root_uptake_b(1:n) * col_cohort%nplant(1:n))
      !----- THIS dt_fast's per-layer sink shares (Phase 1) -- the identical construction the split path  !
      !      does inline, so the three schemes place the root mass AND heat sink in the same layers. -----!
      frozen%roots%root_share(1:nsl) = 0.0_wp
      do i = 1_ik, n
         do k = 1_ik, nsl
            frozen%roots%root_share(k) = frozen%roots%root_share(k) + root_uptake_layer_b(k, i) * col_cohort%nplant(i)
         end do
      end do
      share_tot = sum(frozen%roots%root_share(1:nsl))
      if (share_tot > tiny_num) then
         frozen%roots%root_share(1:nsl) = frozen%roots%root_share(1:nsl) / share_tot
      else
         frozen%roots%root_share(1:nsl) = col_config%soil%root_frac(1:nsl)
      end if

      !----- FROZEN hydrology BCs: the plant's OWN aggregate uptake REQUEST becomes the soil's root-   !
      !      sink forcing (not the raw transpiration demand), then a SCRATCH advance_soil_water_column      !
      !      for soil_evap / infiltration / uptake_total (the soil's TRUE, possibly fwilt-limited        !
      !      realized supply). throughfall_total (sec 3.4, P2c) is what SURVIVES the canopy -- the         !
      !      raw forc%rainfall when canopy_water_on is off (unchanged), or rainfall minus what the             !
      !      interception sweep above caught otherwise; feeding the soil the UN-reduced forc%rainfall         !
      !      here while ALSO crediting the intercepted share to the canopy surface store would create        !
      !      water from nothing (double-counted at the whole-column boundary). ------------------------!
      !----- PRECIP ROUTING under a pack (C4). snow_accumulate has                                      !
      !      ALREADY taken forc%snowfall AND forc%rainfall into the pack, so ONLY meltwater may reach the   !
      !      ground -- adding throughfall on top double-counts the rainfall at the boundary. -----------!
      if (snow_st%exists) then
         hforc%precip_ground   = snow_st%melt_rate + biophys%shed_water_rate
      else
         hforc%precip_ground   = throughfall_total + biophys%shed_water_rate
      end if
      hforc%snow_free_frac     = 1.0_wp - snow_st%snowfac
      hforc%root_uptake(1:nsl) = total_uptake_b * frozen%roots%root_share(1:nsl)
      hforc%t_ground           = t_ground ; hforc%q_air = qcas ; hforc%rho_air = rho
      !----- The hydrology kernel owns the ponding store's ENTHALPY too (#78 item 4): it needs each     !
      !      layer's temperature to value the saturation clip, and the temperature of the water         !
      !      entering the pond. Under a pack that is the MELTWATER temperature, not frozen%hydrology%t_film_valuation --     !
      !      t_film_valuation is pinned to tsupercool_liq so the ledger books no boundary input for melt. -----!
      hforc%soil_temp(1:nsl)   = biophys%soil_e%soil_temp(1:nsl)
      !----- Temperature that VALUES the ground inflow. Under a pack it is the meltwater's. On bare      !
      !      ground it is the EFFECTIVE liquid temperature of the rain + sub-threshold-snowfall mixture:  !
      !      rain arrives as liquid at the canopy-air temperature, snow as ICE at min(t_3ple, air_temp) --   !
      !      the same valuation snow_accumulate gives snowfall that does form a pack -- and the mixture   !
      !      enthalpy per kg is expressed through temp_of_liquid_enthalpy (exact inverse of              !
      !      internal_energy_liquid; below t_3ple it represents water that must still melt, which the    !
      !      pond/soil plateau then does with soil heat). Valuing the snow as liquid at tcas, as this     !
      !      used to, created the fusion enthalpy L_f per kg of sub-threshold snow at the boundary        !
      !      (ledger-consistent, physically wrong; 2026-09 review). -----------------------------------!
      hforc%t_pond_inflow = tcas
      if (snow_st%exists) then
         hforc%t_pond_inflow = snow_st%t_melt
      else if (forc%rainfall + forc%snowfall > tiny_num) then
         hforc%t_pond_inflow = temp_of_liquid_enthalpy(                                                    &
              (forc%rainfall * internal_energy_liquid(tcas)                                            &
               + forc%snowfall * internal_energy_ice(min(t_3ple, forc%air_temp))) / (forc%rainfall + forc%snowfall))
      end if
      !----- Bare-soil aerodynamic resistance, AREA-weighted by the snow-free fraction set above. This !
      !      path used to pin snow_free_frac at 1.0 because it modelled no snow at all; C4's shared     !
      !      stage removed that limitation, so the weighting is real here now. ------------------------!
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)
      soil_w_scratch = biophys%soil_w
      call advance_soil_water_column(soil_w_scratch, hforc, col_config%soil, col_config%soil_water_opts, dt_fast, hflux)
      budget%soil_nsub = hflux%nsub                 ! section 5.3 work counter (same seam on both schemes)
      frozen%ground%soil_evap = hflux%soil_evap
      frozen%hydrology%q_top          = (hflux%infiltration - hflux%soil_evap) / rho_h2o

      !----- Soil-limiting rescale (MEDS_ED2_RK45_DESIGN.md sec 3): scale down ONLY the credit         !
      !      applied to wood_water_mass (not sapflow, an internal wood->leaf transfer) so the whole-    !
      !      column ledger closes to the soil's TRUE realized supply -- scale==1 exactly in the         !
      !      common (non-limited) case. sapflow_frozen/uptake_frozen are what column_derivs' mass ODE    !
      !      reads directly; frozen%roots%uptake (== hflux%uptake_total) is ALSO the soil-water tendency's root    !
      !      sink authority in column_derivs, so the two sides of the wood<->soil interface use the       !
      !      identical number by construction. ---------------------------------------------------------!
      scale = 1.0_wp
      if (total_uptake_b > tiny_num) scale = min(1.0_wp, hflux%uptake_total / total_uptake_b)
      frozen%plant%sapflow_frozen(1:n) = sapflow_b(1:n)
      frozen%plant%uptake_frozen(1:n)  = root_uptake_b(1:n) * scale

      !----- ADVECTIVE ENTHALPY (MEDS_ED2_RK45_DESIGN.md sec 2/6, P2 -- ED2's qwflux_wl/qloss): both     !
      !      the mass FLUX (sapflow_b/uptake_frozen, already frozen above) AND the upwind reference        !
      !      TEMPERATURE are frozen at state^n here -- computing the reference temperature from the        !
      !      CURRENT (this-stage) leaf/wood/soil temperature would be circular (surface_derivs hasn't       !
      !      run yet this call, and it's what PRODUCES those temperatures). Freezing both still closes      !
      !      the energy ledger exactly (the same "one number, both sides" principle sapflow_frozen/          !
      !      uptake_frozen already rely on for water) -- it trades some fidelity in the thermal upwind       !
      !      choice for tractability, a documented, bounded approximation (see the design doc P2 notes).     !
      !      HR is disabled project-wide (uptake floored >=0), so qloss's upwind is unconditionally the       !
      !      root-frac-weighted mean soil temperature (weighted_mean is a generic weighted-sum, not      !
      !      psi-specific, so it is reused verbatim for temperature here). -------------------------------!
      soil_temp_root = weighted_mean(biophys%soil_e%soil_temp(1:nsl), col_config%soil%root_frac, nsl)
      u_liq_soil = internal_energy_liquid(soil_temp_root)
      do i = 1_ik, n
         t_up_wl   = merge(biophys%wood_temp(i), biophys%leaf_temp(i), sapflow_b(i) >= 0.0_wp)
         u_liq_up  = internal_energy_liquid(t_up_wl)
         sapflow_gnd(i) = frozen%plant%sapflow_frozen(i) * col_cohort%nplant(i)   ! [kg/m2 ground/s]
         uptake_gnd(i)  = frozen%plant%uptake_frozen(i)  * col_cohort%nplant(i)   ! [kg/m2 ground/s]
         frozen%tissue%qwflux_wl(i)  = sapflow_gnd(i) * u_liq_up
         frozen%roots%qloss_frozen(i)    = uptake_gnd(i)  * u_liq_soil
         frozen%tissue%q_wood_net(i) = frozen%roots%qloss_frozen(i) - frozen%tissue%qwflux_wl(i)
      end do

      !----- FROZEN boundary hydrology for the guard-lift: the rain/drainage/runoff water-enthalpy       !
      !      advection (state^n temps, matching the split) + the scratch's end-of-step ponding/aquifer/  !
      !      water-table (soil_w_scratch was advanced in place by advance_soil_water_column). ---------------!
      frozen%hydrology%infiltration = hflux%infiltration ; frozen%hydrology%drainage    = hflux%drainage
      frozen%hydrology%precip_ground = hforc%precip_ground
      frozen%hydrology%runoff_surf  = hflux%runoff_surf
      !----- t_film_valuation = tsupercool_liq under a pack makes internal_energy_liquid vanish, so meltwater !
      !      infiltrates its MASS at zero enthalpy -- the enthalpy already moved, paired, inside        !
      !      advance_snow_stage. Without this the melt energy is counted twice at soil layer 1. -------!
      frozen%hydrology%t_film_valuation = hforc%t_pond_inflow   ! one valuation for the boundary inflow AND the film
      if (snow_st%exists) frozen%hydrology%t_film_valuation = tsupercool_liq
      !----- what the film is valued at (see surface_derivs). --------------------------------!
      frozen%film%film_liquid_enthalpy = internal_energy_liquid(frozen%hydrology%t_film_valuation)
      !----- The infiltrating water comes OUT OF THE POND, so the soil top-face advection is        !
      !      referenced to the pond temperature the kernel just reported (#78 item 4). -----------!
      frozen%hydrology%t_infil = hflux%t_infil
      frozen%hydrology%runoff_enth = hflux%runoff_enth
      frozen%roots%uptake       = hflux%uptake_total
      !----- Interior face fluxes + the post-solve mass corrections, from the SAME scratch solve. The   !
      !      clip/floor enthalpies are valued HERE, at each layer's state^n temperature, because that   !
      !      is the temperature the correction has to be neutral against -- and it is the only place    !
      !      the per-layer soil temperature is in scope. ---------------------------------------------!
      frozen%hydrology%w_flux_frozen(1:nsl) = hflux%w_flux(1:nsl)
      do k = 1_ik, nsl
         frozen%hydrology%clip_enth(k)  = hflux%clip_layer(k)  * internal_energy_liquid(biophys%soil_e%soil_temp(k))
         frozen%hydrology%floor_enth(k) = hflux%floor_layer(k) * internal_energy_liquid(biophys%soil_e%soil_temp(k))
      end do
      frozen%hydrology%clip_mass  = sum(hflux%clip_layer(1:nsl))
      frozen%hydrology%floor_mass = sum(hflux%floor_layer(1:nsl))
      frozen%hydrology%t_bot        = biophys%soil_e%soil_temp(nsl)
      frozen%hydrology%w_surface1   = soil_w_scratch%w_surface
      frozen%hydrology%w_surface_enth1 = soil_w_scratch%w_surface_enth
      frozen%hydrology%t_pond_inflow        = hforc%t_pond_inflow
      !----- the AUTHORITATIVE committed soil moisture: soil_w_scratch was advanced IN PLACE by the robust  !
      !      advance_soil_water_column above, so its theta IS the end-of-step (relieved) soil water. -----------!
      allocate(frozen%hydrology%theta1(nsl))
      frozen%hydrology%theta1(1:nsl) = soil_w_scratch%theta(1:nsl)

      !----- pack the prognostic state: plant water MASS is now NATIVE (MEDS_ED2_RK45_DESIGN.md sec 4, !
      !      P2) -- a direct copy from the persisted state, no psi round-trip needed any more. ----------!
      y%cas_enthalpy = biophys%cas%can_enthalpy ; y%cas_shv = biophys%cas%can_shv ; y%cas_co2 = biophys%cas%can_co2
      y%soil_energy(1:nsl) = biophys%soil_e%soil_energy(1:nsl)
      y%theta(1:nsl)       = biophys%soil_w%theta(1:nsl)
      !----- pond packed onto the state vector (#93 Phase 0). Still committed from the scratch      !
      !      hydrology below, so this is carriage only -- no behaviour change. --------------------!
      y%w_surface          = biophys%soil_w%w_surface
      y%w_surface_enth     = biophys%soil_w%w_surface_enth
      y%leaf_water_mass(1:n) = biophys%leaf_water_mass(1:n)
      y%wood_water_mass(1:n) = biophys%wood_water_mass(1:n)
   end subroutine build_column_frozen

   !---------------------------------------------------------------------------------------!
   ! advance_snow_stage -- accumulate snowfall + rain-on-snow, advance the snow-surface energy     !
   ! balance at the LAGGED CAS, and drain meltwater to the PONDING store as a PAIRED (mass,        !
   ! enthalpy) transfer. Mutates the pack store only; everything else it reports through st, so the  !
   ! caller decides how the frozen results reach its own stepper. Inputs are the physical boundary   !
   ! quantities, not the driver's aggregates (2026-09 review, item 4 #7).                           !
   !                                                                                          !
   ! The meltwater's enthalpy is NOT handed to the soil here (it was, before issue #78 item 4 gave   !
   ! the pond a thermal state). It leaves the pack via snow_energy and is reported as melt_enth      !
   ! together with t_melt, the temperature that values it; the caller passes t_melt to the           !
   ! hydrology kernel as chydro_forcing_t%t_pond_inflow, and the ONE pond inflow carries both halves.     !
   ! Pack and pond are both tracked stores, so the transfer telescopes out of the whole-column       !
   ! ledger rather than needing a boundary term -- and no consumer has to rebase a soil baseline.    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_snow_stage(snow, snow_params, dz_soil_top, abs_sw_ground, abs_lw_ground,        &
                                 snowfall, rainfall, t_air, ggnet, t_soil_top, dt_fast, tcas, qcas,    &
                                 rho, press, st)
      type(snow_column_t),  intent(inout) :: snow           !< the pack store (the ONLY state mutated here)
      type(snow_params_t),  intent(in)    :: snow_params
      real(wp),             intent(in)    :: dz_soil_top    !< [m]       top soil-node depth |z_node(1)|
      real(wp),             intent(in)    :: abs_sw_ground  !< [W/m2]    shortwave reaching the ground
      real(wp),             intent(in)    :: abs_lw_ground  !< [W/m2]    net longwave at the ground
      real(wp),             intent(in)    :: snowfall       !< [kg/m2/s] frozen precipitation
      real(wp),             intent(in)    :: rainfall       !< [kg/m2/s] liquid precipitation (rain-on-snow)
      real(wp),             intent(in)    :: t_air          !< [K]       reference-level air temperature
      real(wp),             intent(in)    :: ggnet          !< [m/s]     ground <-> CAS conductance
      real(wp),             intent(in)    :: t_soil_top     !< [K]       top soil-node temperature
      real(wp),             intent(in)    :: dt_fast, tcas, qcas, rho, press
      type(snow_stage_t),   intent(out)   :: st

      type(snow_env_t)  :: senv
      type(snow_flux_t) :: sfx
      type(snow_melt_t) :: smelt
      real(wp)          :: snow_e0

      !----- default = the bare-ground boundary the snow-free column expects (snowfac = 0). --------!
      !                                                                                             !
      !      ALWAYS-ON. There is no `snow_on` switch any more: snowfall is a boundary water input    !
      !      like rain, and a model that receives it must have somewhere to put it. The flag existed  !
      !      because snow was split-only (C4 shared the stage across all three integrators), and      !
      !      while it existed the DEFAULT (.false.) silently discarded frozen precipitation on the     !
      !      ARK/RK45 paths -- precip_phase splits rain from snow without consulting it, so `off`      !
      !      never meant "no snow", it meant "snow with nowhere to go".                                !
      !                                                                                                !
      !      Always-on costs nothing on a snow-free column: snow_accumulate returns immediately unless  !
      !      a pack exists or the snowfall clears params%min_new_snow_mass, so st stays at the bare-     !
      !      ground defaults set just above, snowfac = 0, and surface_derivs' snow blend reduces         !
      !      EXACTLY to its pre-C4 form. Sub-threshold snowfall onto bare ground still reaches the       !
      !      soil as liquid via the caller's throughfall routing -- nothing is dropped either way. ------!
      st%ground_rad = abs_sw_ground + abs_lw_ground
      st%swe0       = snow%swe(1)         ; st%swe1  = snow%swe(1)
      st%enth0      = snow%snow_energy(1) ; st%enth1 = snow%snow_energy(1)
      snow_e0 = snow%snow_energy(1)
      call snow_accumulate(snow, snowfall, rainfall, t_air, dt_fast, snow_params)
      st%acc_enth = snow%snow_energy(1) - snow_e0   ! rainfall enthalpy into the pack (boundary in)
      st%exists   = snow%nlayer >= 1_ik             ! accumulate took snow+rain -> rainfall routing

      if (st%exists .and. snow%swe(1) > snow_params%tiny_snow_mass) then
         !----- SUB-COLUMN: snowfac is snow, (1-snowfac) is bare soil. The pack's boundary exchange   !
         !      is SCALED by snowfac inside snow_energy_step, so a thin/patchy pack barely exchanges  !
         !      -- continuous and stable, with no threshold cliff -- and its returned fluxes are      !
         !      already snowfac-weighted. The bare-soil share is blended by the consumer. -----------!
         st%snowfac       = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), snow_params)
         senv%abs_sw      = abs_sw_ground ; senv%abs_lw = abs_lw_ground
         senv%can_temp    = tcas ; senv%can_shv = qcas ; senv%ggnet = ggnet
         senv%rho_air     = rho ; senv%press = press
         senv%t_soil_top  = t_soil_top
         senv%dz_soil_top = dz_soil_top
         call snow_energy_step(snow, senv, snow_params, dt_fast, st%snowfac, sfx)
         call snow_drain_meltwater(snow, snow_params, smelt)
         st%h_snow  = sfx%h_snow ; st%le_snow = sfx%le_snow ; st%g_base = sfx%g_base
         st%snowfac = sfx%snowfac                         ! the clamped fraction the kernel actually used
         !----- ground radiation boundary in = snow's snowfac-weighted net + bare's (1-snowfac) share !
         st%ground_rad = sfx%rnet + (1.0_wp - st%snowfac) * (abs_sw_ground + abs_lw_ground)
         st%subl_rate  = sfx%w_flux
         st%melt_rate  = (smelt%melt_mass + smelt%dump_mass) / dt_fast
         !----- PAIRED enthalpy: snow store -> soil top (extensive J/m2 -> volumetric J/m3). The mass !
         !      half rides melt_rate into infiltration, and the caller MUST infiltrate it at zero     !
         !      enthalpy (t_film_valuation = tsupercool_liq) or this enthalpy is counted twice. ------------!
         !----- MELTWATER GOES TO THE POND, not straight into soil layer 1 (issue #78 item 4).           !
         !                                                                                              !
         !      The pack used to hand its melt enthalpy directly to soil_energy(1), and the caller then  !
         !      set t_film_valuation = tsupercool_liq so the meltwater MASS infiltrated carrying zero enthalpy  !
         !      -- "the energy already moved, paired, here". That worked while the soil was the only     !
         !      place surface water could go. Once the ponding store has a real thermal state the        !
         !      meltwater ponds FIRST and infiltrates from the pond, so the direct transfer would be     !
         !      counted twice: once here, and again in the pond->layer-1 advection at the mixed pond     !
         !      temperature. (That is the C4 double-count in a new guise.)                               !
         !                                                                                              !
         !      So report the enthalpy and the temperature that values it, and let the ONE pond inflow   !
         !      carry both. The pack still loses it (snow_energy was already debited in                  !
         !      snow_drain_meltwater), the pond gains it, and both are tracked stores -- so it           !
         !      telescopes out of the whole-column ledger instead of needing a boundary term. ----------!
         st%melt_enth = smelt%melt_enth + smelt%dump_enth
         st%t_melt    = t_3ple
         if (smelt%melt_mass + smelt%dump_mass > 0.0_wp)                                            &
            st%t_melt = temp_of_liquid_enthalpy(st%melt_enth / (smelt%melt_mass + smelt%dump_mass))
      end if
      st%swe1  = snow%swe(1)
      st%enth1 = snow%snow_energy(1)
   end subroutine advance_snow_stage

end module meds_fast_frozen
