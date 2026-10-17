! SPDX-License-Identifier: Apache-2.0
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
   use meds_plant_hydraulics, only : cohort_root_profile, crown_gravity_head
   use meds_site_diag_types, only : CD_PSI_WOOD, CD_PLC, CD_SAPFLOW, CD_ROOT_UPTAKE
   use meds_water_retention, only : rwc_from_psi
   use meds_hydr_lib, only : plc_retained, xylem_max_conductance, rhizosphere_conductance
   use meds_spac_solver, only : spac_solver_t, spac_allocate, spac_set_state, spac_advance, POND, WOOD
   use meds_config, only : meds_config_t, pft_wstress_psi_close
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_soil_types, only : chydro_forcing_t, snow_melt_t
   use meds_fast_types, only : patch_biophys_t
   use meds_fast_time_derivs, only : surface_derivs
   use meds_column_state_ops, only : clamp_soil_energy, soil_water_store, soil_energy_store, plant_water_store, &
                                     canopy_film_store, deposit_condensate, clamp_canopy_film, unpack_column_state, &
                                     diagnose_soil_temps, assemble_soil_energy_forcing, apply_process_mask
   use meds_fast_prepass, only : column_prepass
   use meds_fast_types, only : column_config_t, column_cohort_t, column_forcing_t, column_budget_t, column_state_t, &
                               column_frozen_t, surface_state_t, cas_boundary_t, surface_tend_t, stage_bflux_t, &
                               column_bflux_t, error_control_t, column_tend_t, mask_is_full
   use meds_plant_biophysics, only : sensible_heat_coeff, leaf_film_coeff, intercept_canopy_layer
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction
   use meds_therm_lib, only : internal_energy_liquid, temp_of_liquid_enthalpy, precip_enthalpy
   use meds_soil_types, only : snow_env_t, snow_flux_t
   use meds_biophysics_opts, only : snow_params_t
   use meds_column_state_types, only : snow_column_t
   use meds_ground_biophysics, only : snow_energy_step
   use meds_fast_types, only : snow_stage_t
   implicit none
   private

   public :: build_column_frozen

   !----- Prognostic-wood constants. ----------------------------------------------------------------!
   real(wp), parameter :: C2B_WOOD            = 2.0_wp  !< carbon -> biomass (carbon fraction 0.5)
   real(wp), parameter :: WOOD_MOIST_FRAC_ARK = 1.0_wp  !< [kg water/kg dry] fresh-sapwood moisture (MVP)

contains

   !----- Build the frozen ARK inputs: the shared column_prepass (meds_fast_prepass: leaf gas         !
   !      respiration / CAS caps / aero, bit-identical to the split) + this integrator's own per-cohort !
   !      geometry/radiation/wood packing + the frozen hydrology BCs, into a column_frozen_t; also      !
   !      packs the prognostic state into a column_state_t.                                             !
   subroutine build_column_frozen(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, spac, &
                                  n, nsl, frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
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
      type(spac_solver_t),     intent(inout) :: spac     !< the thread's soil-plant water solver, reused step to step
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: frozen
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) per-cohort diagnostic capture

      type(snow_stage_t)     :: snow_st   !< shared pre-column snow stage (C4, issue #76)
      type(chydro_forcing_t) :: hforc
      type(surface_state_t)  :: y_stage ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, nee_biotic, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
           g_atm_co2
      integer(ik) :: i, k
      !----- The water network's step, read into the frozen record. ------------------------------!
      real(wp) :: share_tot, face_resid, f_in, f_out, u_liq_soil_k(nsl)
      real(wp) :: root_frac(nsl, n), root_frac_patch(nsl), z_top(nsl), z_bot(nsl), root_carbon
      real(wp) :: t_up_wl, u_liq_up, u_liq_leaf, u_liq_wood, sapflow_gnd, uptake_gnd, transp_gnd
      !----- Canopy-SURFACE water pre-pass scratch (MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c). ------------!
      real(wp) :: rain_above, combined_w, leaf_cap_i, wood_cap_i, throughfall_i, drip_i, throughfall_total
      real(wp) :: avail_leaf, avail_wood

      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n),                  &
         frozen%tissue%abs_lw(n), frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      allocate(frozen%tissue%wood_dry_hcap(n), frozen%tissue%wood_wmass(n))
      allocate(frozen%tissue%leaf_dry_hcap(n), frozen%tissue%leaf_wmass(n))
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      allocate(frozen%tissue%leaf_emiss(n), frozen%tissue%wood_emiss(n))
      !----- ZERO the advective-enthalpy terms AT ALLOCATION. They are only given their real values    !
      !      further down (after the water network supplies sapflow/uptake), but the sf0                !
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
      allocate(frozen%roots%root_share(nsl), frozen%plant%nplant(n))
      allocate(frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n), frozen%roots%qloss_frozen(n))
      allocate(frozen%roots%root_heat(nsl), frozen%tissue%transp_net(n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n))
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water(1:n) = biophys%leaf_surf_water(1:n) ; y%wood_surf_water(1:n) = biophys%wood_surf_water(1:n)

      !----- Each cohort's fine roots by layer, within its own rooting depth, and the patch's profile:   !
      !      the cohorts' weighted by their fine-root carbon. A patch with no roots puts its weight in the !
      !      top layer (it has no root respiration to weight). -----------------------------------------!
      z_top(1:nsl) = -col_config%soil%soil_layer_z(1:nsl)
      z_bot(1:nsl) = -col_config%soil%soil_layer_z(2:nsl+1)
      root_frac_patch(1:nsl) = 0.0_wp
      do i = 1_ik, n
         call cohort_root_profile(col_cohort%height(i), col_config%hydraulics_table%pft(col_cohort%pft(i)),       &
                                  z_top(1:nsl), z_bot(1:nsl), root_frac(1:nsl, i))
         root_frac_patch(1:nsl) = root_frac_patch(1:nsl) + col_cohort%nplant(i)*col_cohort%broot(i)*root_frac(1:nsl, i)
      end do
      root_carbon = sum(root_frac_patch(1:nsl))
      if (root_carbon > tiny_num) then
         root_frac_patch(1:nsl) = root_frac_patch(1:nsl) / root_carbon
      else
         root_frac_patch(1:nsl) = 0.0_wp ; root_frac_patch(1) = 1.0_wp
      end if

      !----- the SHARED pre-pass (meds_fast_prepass%column_prepass): gas exchange / respiration / CAS   !
      !      aero -- writes directly into the frozen struct's h_coeff_leaf/g_transp_leaf arrays. ------------------!
      call column_prepass(cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, root_frac_patch(1:nsl), &
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
            combined_w = biophys%leaf_surf_water(i) + biophys%wood_surf_water(i)
            leaf_cap_i = col_config%leaf_surf_water_max(col_cohort%pft(i)) * col_cohort%lai(i)
            wood_cap_i = col_config%wood_surf_water_max(col_cohort%pft(i)) * col_cohort%wai(i)
            call intercept_canopy_layer(combined_w, rain_above, col_cohort%lai(i), col_cohort%wai(i), 0.0_wp, dt_fast, &
                                        col_config%leaf_surf_water_max(col_cohort%pft(i)),                  &
                                        col_config%wood_surf_water_max(col_cohort%pft(i)),                  &
                                        col_config%soil_water_opts%intercept_k, &
                                             col_config%soil_water_opts%intercept_alpha, &
                                        throughfall_i, drip_i, frozen%film%f_wet_c(i))
            !----- the combined film shared by capacity, so neither film is left above its own cap -------!
            if (leaf_cap_i + wood_cap_i > tiny_num) then
               frozen%film%intercept_leaf(i) = (combined_w*leaf_cap_i/(leaf_cap_i + wood_cap_i)                 &
                                               - biophys%leaf_surf_water(i)) / dt_fast
               frozen%film%intercept_wood(i) = (combined_w*wood_cap_i/(leaf_cap_i + wood_cap_i)                 &
                                               - biophys%wood_surf_water(i)) / dt_fast
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
         !----- TISSUE HEAT STORE, frozen for the whole dt_fast (docs/science/vegetation_energy_dynamics.md). !
         !      a = cap/dt_fast; the relaxation origin is the START-of-step tissue temperature. Every     !
         !      stage evaluation therefore returns the same dt_fast-averaged flux and dt_fast-endpoint    !
         !      temperature -- the store is an algebraic closure, not a tableau DOF. -------------------!
         frozen%tissue%leaf_hcap_per_dt(i)  = (frozen%tissue%leaf_dry_hcap(i) + frozen%tissue%leaf_wmass(i) * cp_liq) &
                                              / dt_fast
         frozen%tissue%wood_hcap_per_dt(i)  = (frozen%tissue%wood_dry_hcap(i) + frozen%tissue%wood_wmass(i) * cp_liq) &
                                              / dt_fast
         frozen%tissue%t_leaf0(i) = biophys%leaf_temp(i)
         frozen%tissue%t_wood0(i) = biophys%wood_temp(i)
         frozen%tissue%leaf_emiss(i) = cfg%pft%leaf_emissivity(col_cohort%pft(i))
         frozen%tissue%wood_emiss(i) = cfg%pft%wood_emissivity(col_cohort%pft(i))
         frozen%plant%nplant(i)   = col_cohort%nplant(i)
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
      !      lost) rather than genuinely closing the ledger. sf0 is evaluated again below, so the            !
      !      transpiration demand the water network is given sees the final, scaled conductance.               !
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
      end if
      !----- The transpiration demand the water network is given is this solve's, so the leaf here    !
      !      must be charged what the stages will charge it. The leaf pays the full vapour enthalpy of  !
      !      what it transpires and is credited the liquid enthalpy of the sapflow that replaces it      !
      !      (qwflux_wl, set from the network's flows below). Without the credit this solve charges the  !
      !      leaf ~4.4e5 J/kg too much, a cooler leaf transpires ~4% less, and the stages hold that.     !
      !      So the credit is estimated from the first solve's transpiration, as if the sapflow matched  !
      !      it, and the solve is repeated. ----------------------------------------------------------!
      do i = 1_ik, n
         t_up_wl = merge(biophys%wood_temp(i), biophys%leaf_temp(i), sf0%transp_c(i) >= 0.0_wp)
         frozen%tissue%qwflux_wl(i) = sf0%transp_c(i) * internal_energy_liquid(t_up_wl)
      end do
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, t_ground, n, sf0)
      !----- THE BOUNDARY WATER reaching the ground (C4). Under a pack, snow_accumulate has already      !
      !      taken forc%snowfall AND forc%rainfall into the pack, so only meltwater reaches the ground;     !
      !      on bare ground it is what SURVIVES the canopy, throughfall_total (sec 3.4, P2c) -- feeding     !
      !      the soil the un-intercepted rainfall while also crediting the intercepted share to the         !
      !      canopy film would create water at the whole-column boundary. ------------------------------!
      if (snow_st%exists) then
         hforc%precip_ground   = snow_st%melt_rate + biophys%shed_water_rate
      else
         hforc%precip_ground   = throughfall_total + biophys%shed_water_rate
      end if
      hforc%snow_free_frac     = 1.0_wp - snow_st%snowfac
      hforc%t_ground           = t_ground ; hforc%q_air = qcas ; hforc%rho_air = rho
      !----- Phase of the evaporating surface, for the ice saturation branch (#89). ------------------!
      hforc%ground_fliq        = biophys%soil_e%soil_fliq(1)
      !----- Temperature that VALUES the ground inflow. Under a pack it is the meltwater's. On bare      !
      !      ground it is the EFFECTIVE liquid temperature of the rain + sub-threshold-snowfall mixture:  !
      !      precip_enthalpy values both at the air temperature above the canopy (forc%air_temp, the      !
      !      forcing moved to the canopy-air top), exactly as snow_accumulate does for a pack, and the    !
      !      mixture enthalpy per kg is expressed through temp_of_liquid_enthalpy (exact inverse of       !
      !      internal_energy_liquid; below t_3ple it represents water that must still melt, which the     !
      !      pond/soil plateau then does with soil heat). Valuing the snow as liquid at tcas created the  !
      !      fusion enthalpy L_f per kg of sub-threshold snow at the boundary (2026-09 review). Rain was    !
      !      valued at the canopy-air temperature until #355: the rain falls through the air above the     !
      !      canopy, and at a cold front the canopy air is still warm (BCI 2017-04-17: rain at 297.2 K     !
      !      into 294.9 K air), which handed the canopy and soil the rain's missing cooling. ----------!
      hforc%t_pond_inflow = tcas
      if (snow_st%exists) then
         hforc%t_pond_inflow = snow_st%t_melt
      else if (forc%rainfall + forc%snowfall > tiny_num) then
         hforc%t_pond_inflow = temp_of_liquid_enthalpy(precip_enthalpy(forc%rainfall, forc%snowfall, forc%air_temp) &
                                                       / (forc%rainfall + forc%snowfall))
      end if
      !----- Bare-soil aerodynamic resistance, AREA-weighted by the snow-free fraction set above. -----!
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)

      !----- THE SOIL-PLANT WATER (meds_spac_solver; MEDS_SPAC_WATER_NETWORK_PLAN.md §4.4): the pond, the !
      !      soil layers and every cohort's wood and leaves, solved together over the step with each      !
      !      cohort's transpiration demand from the start-of-step surface solve above. Everything the     !
      !      march needs about water is read from it below; the stages hold it fixed, as they always       !
      !      held ground evaporation. --------------------------------------------------------------!
      call run_spac_solver(dt_fast, cfg, col_config, col_cohort, biophys, hforc, root_frac(1:nsl, 1:n),        &
                           sf0%transp_c(1:n), n, nsl, spac)

      !----- Work and convergence (section 5.3 work counters): the network's accepted sub-steps, its    !
      !      Newton iterations, and the sub-steps it had to take although they did not converge or meet  !
      !      the tolerance. -------------------------------------------------------------------------!
      budget%water_nsub     = spac%nsub
      budget%water_newton    = spac%nnewton
      budget%water_nfail = spac%nfail
      budget%water_reject  = spac%nreject

      !----- Ground evaporation, the soil's boundary flows and the pond, as the network moved them. ----!
      frozen%ground%soil_evap         = spac%evap / dt_fast
      frozen%hydrology%infiltration   = spac%infil / dt_fast
      frozen%hydrology%exfiltration   = spac%exfil / dt_fast
      frozen%hydrology%exfil_enth     = spac%exfil_enth / dt_fast
      frozen%hydrology%q_top          = (spac%infil - spac%exfil - spac%evap) / (rho_h2o * dt_fast)
      frozen%hydrology%drainage       = spac%drain / dt_fast
      frozen%hydrology%runoff_surf    = spac%spill / dt_fast
      budget%runoff_surf_export       = spac%spill
      budget%drainage_export          = spac%drain
      frozen%hydrology%runoff_enth    = spac%spill_enth / dt_fast
      frozen%hydrology%precip_ground  = hforc%precip_ground
      frozen%hydrology%t_pond_inflow  = hforc%t_pond_inflow
      !----- The infiltrating water left the POND, so the soil top face is valued at the temperature    !
      !      that carries exactly the enthalpy the pond gave up (#78 item 4). -------------------------!
      frozen%hydrology%t_infil = hforc%t_pond_inflow
      if (spac%infil > tiny_num) frozen%hydrology%t_infil = temp_of_liquid_enthalpy(spac%infil_enth / spac%infil)
      !----- t_film_valuation = tsupercool_liq under a pack makes internal_energy_liquid vanish, so meltwater !
      !      infiltrates its MASS at zero enthalpy -- the enthalpy already moved, paired, inside        !
      !      advance_snow_stage. Without this the melt energy is counted twice at soil layer 1. -------!
      frozen%hydrology%t_film_valuation = hforc%t_pond_inflow   ! one valuation for the boundary inflow AND the film
      if (snow_st%exists) frozen%hydrology%t_film_valuation = tsupercool_liq
      frozen%film%film_liquid_enthalpy = internal_energy_liquid(frozen%hydrology%t_film_valuation)
      frozen%hydrology%t_bot           = biophys%soil_e%soil_temp(nsl)
      frozen%hydrology%w_surface1      = spac%water_soil(POND)
      frozen%hydrology%w_surface_enth1 = spac%e_pond
      frozen%hydrology%w_flux_frozen(1:nsl)   = 0.0_wp
      frozen%hydrology%w_flux_frozen(1:nsl-1) = spac%face(1:nsl-1) / (rho_h2o * dt_fast)
      allocate(frozen%hydrology%theta1(nsl))
      do k = 1_ik, nsl
         frozen%hydrology%theta1(k) = spac%water_soil(1_ik + k) / (rho_h2o * col_config%soil%dz(k))
      end do
      frozen%hydrology%theta_rate(1:nsl) = (frozen%hydrology%theta1(1:nsl) - biophys%soil_w%theta(1:nsl)) / dt_fast

      !----- PER-LAYER FACE CLOSURE (#189): what each layer's water changed by, against the flows the   !
      !      soil energy will advect -- faces, root uptake, and at the top infiltration, exfiltration    !
      !      and evaporation. The network moves every store's water by these flows, so this reads        !
      !      round-off; it is kept as the check that the faces handed on are the ones that moved the    !
      !      water. A network sub-step taken unconverged counts as a failed check, so the end-of-run      !
      !      report says so. -----------------------------------------------------------------------!
      face_resid = 0.0_wp
      do k = 1_ik, nsl
         f_in = 0.0_wp ; f_out = 0.0_wp
         if (k > 1_ik)   f_in  = spac%face(k-1)
         if (k < nsl)    f_out = spac%face(k)
         if (k == nsl)   f_out = spac%drain
         if (k == 1_ik)  f_in  = spac%infil - spac%exfil - spac%evap
         face_resid = face_resid + abs((frozen%hydrology%theta1(k) - biophys%soil_w%theta(k))              &
                                       * col_config%soil%dz(k) * rho_h2o - (f_in - f_out - sum(spac%uptake(k, 1:n))))
      end do
      budget%soil_face_mass%resid   = face_resid
      budget%soil_face_mass%worst   = max(budget%soil_face_mass%worst, face_resid)
      budget%soil_face_mass%abs_sum = budget%soil_face_mass%abs_sum + face_resid
      budget%soil_face_mass%n_check = budget%soil_face_mass%n_check + 1_ik
      if (spac%nfail > 0_ik) budget%soil_face_mass%n_fail = budget%soil_face_mass%n_fail + 1_ik

      !----- Root uptake: in total, by layer, and by cohort; and the enthalpy it takes out of each      !
      !      layer at that layer's start-of-step temperature (plan §12.2 M4) -- the same number the    !
      !      cohort's wood is credited below. A patch whose roots took nothing places its (zero) sink    !
      !      by its root profile. --------------------------------------------------------------------!
      frozen%roots%uptake = sum(spac%uptake(1:nsl, 1:n)) / dt_fast
      do k = 1_ik, nsl
         u_liq_soil_k(k) = internal_energy_liquid(biophys%soil_e%soil_temp(k))
         frozen%roots%root_share(k) = sum(spac%uptake(k, 1:n))
         frozen%roots%root_heat(k)  = sum(spac%uptake(k, 1:n)) * u_liq_soil_k(k) / dt_fast
      end do
      share_tot = sum(frozen%roots%root_share(1:nsl))
      if (share_tot > tiny_num) then
         frozen%roots%root_share(1:nsl) = frozen%roots%root_share(1:nsl) / share_tot
      else
         frozen%roots%root_share(1:nsl) = root_frac_patch(1:nsl)
      end if
      do i = 1_ik, n
         frozen%roots%qloss_frozen(i)     = sum(spac%uptake(1:nsl, i) * u_liq_soil_k(1:nsl)) / dt_fast
         frozen%plant%sapflow_frozen(i)  = spac%sapflow(i) / (max(col_cohort%nplant(i), tiny_num) * dt_fast)
         frozen%plant%uptake_frozen(i)   = sum(spac%uptake(1:nsl, i)) / (max(col_cohort%nplant(i), tiny_num) * dt_fast)
         frozen%tissue%transp_net(i)     = spac%transp(i) / dt_fast
      end do
      frozen%tissue%transp_held = .true.

      !----- Per-cohort HYDRAULIC diagnostics, from the same network the physics commits: sapflow and   !
      !      root uptake per plant, and the wood's end-of-step potential and loss of conductance. -------!
      if (present(cdiag)) then
         do i = 1_ik, n
            cdiag(CD_SAPFLOW,     i) = frozen%plant%sapflow_frozen(i)
            cdiag(CD_ROOT_UPTAKE, i) = frozen%plant%uptake_frozen(i)
            cdiag(CD_PSI_WOOD,    i) = spac%psi_coh(WOOD, i)
            associate (hp => col_config%hydraulics_table%pft(col_cohort%pft(i)))
               cdiag(CD_PLC, i) = 1.0_wp - plc_retained(spac%psi_coh(WOOD, i), hp%wood_psi50, hp%wood_kexp)
            end associate
         end do
      end if

      !----- ADVECTIVE ENTHALPY of the plant water (MEDS_ED2_RK45_DESIGN.md sec 2/6, P2 -- ED2's        !
      !      qwflux_wl/qloss): the mass flows are the network's, and the upwind temperatures are frozen   !
      !      at state^n -- computing them from the current (stage) leaf/wood temperature would be         !
      !      circular, since surface_derivs is what produces those temperatures. Freezing both still      !
      !      closes the energy ledger exactly (one number, both sides). The root water arrives at each     !
      !      layer's own temperature (qloss_frozen, above).                                               !
      !                                                                                                  !
      !      Each tissue's own WATER is valued at the tissue's start-of-step temperature, the way the        !
      !      canopy film is valued at the liquid enthalpy its water arrived with. The tissue heat store is   !
      !      cap*T with the step's frozen capacity, which has no term for water mass coming or going; so    !
      !      the water the tissue gains or loses is booked at that value (water_store_enth, in the ledgers), !
      !      and the tissue's temperature sees only how far the arriving water's temperature is from its   !
      !      own. Counting the arriving water's whole enthalpy (~4.4e5 J/kg from the liquid datum) as heat   !
      !      warmed the wood whenever it refilled: ~10 W/m2 every night at BCI, and at the 2017-04-17 storm  !
      !      front, when wet soil refilled the dry-season wood with ~7 kg/m2 in a quarter hour, it put the  !
      !      wood 13 K above the canopy air and drove H to +600 W/m2 at zero net radiation (#355). The leaf  !
      !      still pays the full vapour enthalpy of what it transpires (surface_derivs). -------------------!
      frozen%tissue%water_store_enth = 0.0_wp
      do i = 1_ik, n
         t_up_wl    = merge(biophys%wood_temp(i), biophys%leaf_temp(i), spac%sapflow(i) >= 0.0_wp)
         u_liq_up   = internal_energy_liquid(t_up_wl)
         u_liq_leaf = internal_energy_liquid(biophys%leaf_temp(i))
         u_liq_wood = internal_energy_liquid(biophys%wood_temp(i))
         sapflow_gnd = spac%sapflow(i) / dt_fast                                  ! [kg/m2 ground/s]
         uptake_gnd  = sum(spac%uptake(1:nsl, i)) / dt_fast
         transp_gnd  = spac%transp(i) / dt_fast
         frozen%tissue%qwflux_wl(i)  = sapflow_gnd * u_liq_up - (sapflow_gnd - transp_gnd) * u_liq_leaf
         frozen%tissue%q_wood_net(i) = frozen%roots%qloss_frozen(i) - uptake_gnd * u_liq_wood                &
                                       - sapflow_gnd * (u_liq_up - u_liq_wood)
         frozen%tissue%water_store_enth = frozen%tissue%water_store_enth                                     &
                                        + (uptake_gnd - sapflow_gnd) * u_liq_wood                             &
                                        + (sapflow_gnd - transp_gnd) * u_liq_leaf
      end do

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
   ! run_spac_solver -- one patch's soil and plant water over the fast step, as one network      !
   ! (meds_spac_solver; MEDS_SPAC_WATER_NETWORK_PLAN.md §4.4). Builds the solver's inputs from the !
   ! patch's columns -- each cohort's xylem conductance from its sapwood, its rooting from its      !
   ! height and fine-root carbon, and each layer's root spacing from the patch's root-length        !
   ! density -- starts it from the stores' state^n water, and advances it. e_demand is each         !
   ! cohort's transpiration demand from the start-of-step surface solve [kg/m2 ground/s]; the        !
   ! network takes less where a leaf runs short of water (§12.2 C1).                                 !
   !---------------------------------------------------------------------------------------!
   subroutine run_spac_solver(dt_fast, cfg, col_config, col_cohort, biophys, hforc, root_frac, e_demand, &
                              n, nsl, spac)
      real(wp),               intent(in)    :: dt_fast
      type(meds_config_t),    intent(in)    :: cfg
      type(column_config_t),  intent(in)    :: col_config
      type(column_cohort_t),  intent(in)    :: col_cohort
      type(patch_biophys_t),  intent(in)    :: biophys
      type(chydro_forcing_t), intent(in)    :: hforc
      integer(ik),            intent(in)    :: n, nsl
      real(wp),               intent(in)    :: root_frac(nsl, n)   !< [-] each cohort's fine roots by layer
      real(wp),               intent(in)    :: e_demand(n)         !< [kg/m2 ground/s]
      type(spac_solver_t),    intent(inout) :: spac
      real(wp)    :: root_len(nsl, n), half_spacing(nsl)
      integer(ik) :: i, k
      call spac_allocate(spac, nsl, n)
      do i = 1_ik, n
         associate (hp => col_config%hydraulics_table%pft(col_cohort%pft(i)))
            spac%pft(i)      = col_cohort%pft(i)
            spac%nplant(i)   = col_cohort%nplant(i)
            spac%bleaf(i)    = col_cohort%bleaf(i)
            spac%bwood(i)    = col_cohort%bsap(i) + col_cohort%broot(i)
            spac%grav(i)     = crown_gravity_head(col_cohort%height(i), col_config%hydraulics_opts)
            spac%k_xylem(i)  = col_cohort%nplant(i) * xylem_max_conductance(hp%wood_kmax, col_cohort%sap_area(i), &
                                                                           col_cohort%height(i), hp%vessel_curl)
            spac%e_demand(i) = e_demand(i)
            spac%r_close(i)  = rwc_from_psi(pft_wstress_psi_close(cfg, col_cohort%pft(i)), hp%leaf_curve%pi0,      &
                                           hp%leaf_curve%elastic_mod)
            root_len(1:nsl, i) = col_cohort%broot(i) * root_frac(1:nsl, i) * hp%specific_root_length
         end associate
      end do
      !----- Half the distance between a layer's roots, (pi L_v)^-1/2 from the root-length density of  !
      !      every cohort there, and each cohort's rhizosphere conductance in the layer. ------------!
      do k = 1_ik, nsl
         half_spacing(k) = 1.0_wp / sqrt(pi * max(sum(col_cohort%nplant(1:n) * root_len(k, 1:n))                &
                                                  / col_config%soil%dz(k), tiny_num))
         do i = 1_ik, n
            spac%g_rhizo(k, i) = rhizosphere_conductance(col_cohort%nplant(i), root_len(k, i), half_spacing(k),   &
                                                         col_config%hydraulics_table%pft(col_cohort%pft(i))%fine_root_radius)
         end do
      end do
      spac%rain       = hforc%precip_ground
      spac%t_rain     = hforc%t_pond_inflow
      spac%t_soil(1:nsl) = biophys%soil_e%soil_temp(1:nsl)
      spac%w_pond_max = col_config%soil_water_opts%w_pond_max
      spac%rtol       = col_config%soil_water_opts%rtol
      spac%atol       = col_config%soil_water_opts%atol
      call spac_set_state(spac, col_config%soil, col_config%hydraulics_table, biophys%soil_w%theta(1:nsl),      &
                          biophys%soil_w%w_surface, biophys%soil_w%w_surface_enth,                              &
                          biophys%leaf_water_mass(1:n), biophys%wood_water_mass(1:n))
      call spac_advance(spac, col_config%soil, col_config%soil_water_opts, hforc, col_config%hydraulics_table,    &
                        dt_fast)
   end subroutine run_spac_solver

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
