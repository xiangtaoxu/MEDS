! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_frozen -- build the frozen work record one fast step marches against.          !
!                                                                                          !
! build_column_frozen fills column_frozen_t, everything both schemes (ARK, RK45) hold      !
! constant across a dt_fast: the canopy-air boundary, tissue coefficients, canopy film,    !
! ground boundary, snow stage, soil hydrology, root zone and plant water, from the shared  !
! pre-pass, the snow stage and the soil-plant water network's step.                        !
! docs/science/numerical_scheme.md §2; docs/science/soil_plant_water.md §4-§5.              !
!==========================================================================================!
module meds_fast_frozen
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air, rho_h2o, pi, tsupercool_liq, cp_liq, t_3ple, tissue_c_leaf,       &
                              tissue_c_wood_dry
   use meds_plant_hydraulics, only : cohort_root_profile, crown_gravity_head
   use meds_diag_vars,       only : CO_PSI_WOOD, CO_PLC, CO_SAPFLOW, CO_ROOT_UPTAKE
   use meds_plant_hydr_lib,  only : conductance_retained_fraction, xylem_max_conductance,                         &
                                    rhizosphere_conductance
   use meds_spac_solver, only : spac_solver_t, spac_allocate, spac_set_state, spac_advance, POND, WOOD
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_soil_types, only : chydro_forcing_t, snow_melt_t, snow_env_t, snow_flux_t
   use meds_fast_types, only : patch_biophys_t, column_config_t, column_cohort_t, column_forcing_t,               &
                               column_budget_t, column_state_t, column_frozen_t, surface_state_t,                 &
                               surface_tend_t, snow_stage_t, column_site_t
   use meds_fast_time_derivs, only : surface_derivs
   use meds_fast_prepass, only : column_prepass
   use meds_plant_biophysics, only : sensible_heat_coeff, leaf_film_coeff, intercept_canopy_layer, EFFAREA_EVAP
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction,                 &
                                      snow_energy_step
   use meds_therm_lib, only : internal_energy_liquid, liquid_energy_to_temp, precip_enthalpy, tissue_dry_heat_capacity
   use meds_biophysics_params, only : snow_params_t
   use meds_patch_state_vars,   only : snow_column_t
   implicit none
   private

   public :: build_column_frozen

contains

   !----- Build the frozen record both schemes march against (docs/science/numerical_scheme.md §2):     !
   !      the shared column_prepass, the snow stage, canopy interception, the per-cohort tissue         !
   !      coefficients, the soil-plant water network's step and the plant water's advective enthalpy;   !
   !      also packs the prognostic state into y. -------------------------------------------------------!
   subroutine build_column_frozen(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, spac, &
                                  n, nsl, frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(column_config_t),   intent(in)    :: col_config
      type(column_site_t),     intent(in)    :: col_site
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      !----- inout: the snow stage advances the pack (biophys%snow); every other use is read-only. -------!
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      type(spac_solver_t),     intent(inout) :: spac     !< the thread's soil-plant water solver, reused step to step
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: frozen
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CO_FAST_ROWS, ncoh) per-cohort diagnostic capture

      type(snow_stage_t)     :: snow_st   !< the shared pre-column snow stage's outcome
      type(chydro_forcing_t) :: hforc
      type(surface_state_t)  :: y_stage
      type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, nee_biotic, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
           g_atm_co2
      integer(ik) :: i, k
      !----- The water network's step, read into the frozen record. ------------------------------!
      real(wp) :: share_tot, face_resid, f_in, f_out, u_liq_soil_k(nsl)
      real(wp) :: root_frac(nsl, n), root_frac_patch(nsl), z_top(nsl), z_bot(nsl), root_carbon
      real(wp) :: t_up_wl, u_liq_up, u_liq_leaf, u_liq_wood, sapflow_gnd, uptake_gnd, transp_gnd
      !----- Canopy interception scratch. -------------------------------------------------------------!
      real(wp) :: rain_above, combined_w, leaf_cap_i, wood_cap_i, throughfall_i, drip_i, throughfall_total
      real(wp) :: avail_leaf, avail_wood
      !----- The water reaching the ground besides the shed water [kg/m2/s], its temperature [K], and   !
      !      the enthalpy all of it brings [W/m2]. ---------------------------------------------------!
      real(wp) :: water_base, t_water_base, inflow_enth, precip_ground, t_pond_inflow

      allocate(frozen%tissue%h_coeff_leaf(n), frozen%tissue%g_transp_leaf(n), frozen%tissue%abs_sw(n),                  &
         frozen%tissue%abs_lw(n), frozen%tissue%lai(n))
      allocate(frozen%tissue%h_coeff_w(n), frozen%tissue%abs_sw_wood(n), frozen%tissue%abs_lw_wood(n), frozen%tissue%wai(n))
      allocate(frozen%tissue%wood_dry_hcap(n), frozen%tissue%wood_wmass(n))
      allocate(frozen%tissue%leaf_dry_hcap(n), frozen%tissue%leaf_wmass(n))
      allocate(frozen%tissue%leaf_hcap_per_dt(n), frozen%tissue%wood_hcap_per_dt(n),                &
               frozen%tissue%t_leaf0(n), frozen%tissue%t_wood0(n))
      allocate(frozen%tissue%qwflux_wl(n), frozen%tissue%q_wood_net(n))
      allocate(frozen%tissue%leaf_water_enth(n), frozen%tissue%wood_water_enth(n))
      allocate(frozen%tissue%leaf_emiss(n), frozen%tissue%wood_emiss(n))
      !----- Zeroed now: the start-of-step surface solves below read them (as q_extra, whose default is   !
      !      0) before the water network's flows give them their values. ----------------------------------!
      frozen%tissue%qwflux_wl(1:n)  = 0.0_wp
      frozen%tissue%q_wood_net(1:n) = 0.0_wp
      allocate(frozen%film%f_wet_c(n), frozen%film%g_film_leaf(n), frozen%film%g_film_w(n))
      allocate(frozen%roots%root_share(nsl), frozen%plant%nplant(n))
      allocate(frozen%plant%sapflow_frozen(n), frozen%plant%uptake_frozen(n), frozen%roots%qloss_frozen(n))
      allocate(frozen%roots%root_heat(nsl), frozen%tissue%transp_net(n))
      allocate(frozen%film%intercept_leaf(n), frozen%film%intercept_wood(n))
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n))
      allocate(y%leaf_surf_water(n), y%wood_surf_water(n))
      y%leaf_surf_water(1:n) = biophys%leaf_surf_water(1:n)
      y%wood_surf_water(1:n) = biophys%wood_surf_water(1:n)

      !----- Each cohort's fine roots by layer, within its own rooting depth, and the patch's profile:   !
      !      the cohorts' weighted by their fine-root carbon. A patch with no roots puts its weight in the !
      !      top layer (it has no root respiration to weight). -----------------------------------------!
      z_top(1:nsl) = -col_site%soil%soil_layer_z(1:nsl)
      z_bot(1:nsl) = -col_site%soil%soil_layer_z(2:nsl+1)
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
         root_frac_patch(1:nsl) = 0.0_wp
         root_frac_patch(1) = 1.0_wp
      end if

      !----- The shared pre-pass: gas exchange, respiration, aerodynamics, CAS capacities. --------------!
      call column_prepass(col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, root_frac_patch(1:nsl), &
                          tcas, qcas, press, rho, t_ground, frozen%tissue%h_coeff_leaf, frozen%tissue%g_transp_leaf,      &
                          cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, nee_biotic, &
                          gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      !----- The shared snow stage (docs/science/snow_biophysics.md §The snow store): after the pre-pass, !
      !      whose CAS state, air properties and ground conductance it reads; before the water network, !
      !      so this step's meltwater reaches the pond this step. A no-op without a pack. -------------!
      call advance_snow_stage(biophys%snow, col_config%snow, col_config%aero%snow_rough,                   &
                              max(-col_site%soil%z_node(1), tiny_num),                                      &
                              forc%abs_sw_ground, forc%abs_lw_ground, forc%snowfall, forc%rainfall, forc%air_temp,  &
                              aero%ggnet, biophys%soil_energy%soil_temp(1), dt_fast, tcas, qcas, rho, press,    &
                              snow_st)
      frozen%snow = snow_st                       ! the stage's whole outcome, one record

      !----- Canopy interception (docs/science/soil_biophysics.md §Canopy interception): one leaf+wood     !
      !      bucket per cohort, tallest first, capture only (film evaporation is surface_derivs' per-stage !
      !      flux), held as the frozen rate whose Euler step over dt_fast is the bucket update. None under  !
      !      a pack; throughfall carries sub-threshold snowfall, which forms no pack, as liquid. ----------!
      frozen%film%intercept_leaf = 0.0_wp
      frozen%film%intercept_wood = 0.0_wp
      frozen%film%f_wet_c = 0.0_wp
      frozen%film%g_film_leaf = 0.0_wp
      frozen%film%g_film_w = 0.0_wp
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
                                        col_site%soil_water_opts%intercept_k, &
                                             col_site%soil_water_opts%intercept_alpha, &
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

      !----- Per-cohort geometry, radiation and tissue heat-capacity inputs. ---------------------------!
      do i = 1_ik, n
         frozen%tissue%lai(i)    = col_cohort%lai(i)
         frozen%tissue%abs_sw(i) = forc%abs_sw(i)
         frozen%tissue%abs_lw(i) = forc%abs_lw(i)
         !----- Wood heat capacity: the dry heat of all the wood plus the water network's whole wood store  !
         !      (docs/science/vegetation_energy_dynamics.md eq. (8), §The store is ON). -------------------!
         frozen%tissue%wood_dry_hcap(i) = tissue_dry_heat_capacity(col_cohort%nplant(i), col_cohort%bwood(i),     &
                                                                   tissue_c_wood_dry)
         frozen%tissue%wood_wmass(i)    = max(biophys%wood_water_mass(i), 0.0_wp) * col_cohort%nplant(i)
         !----- Leaf heat capacity, the same two parts, its water read from biophys (y's is filled at the  !
         !      end of this routine). The intercepted film is a separate store, not in it. ---------------!
         frozen%tissue%leaf_dry_hcap(i) = tissue_dry_heat_capacity(col_cohort%nplant(i), col_cohort%bleaf(i),     &
                                                                   tissue_c_leaf)
         frozen%tissue%leaf_wmass(i)    = max(biophys%leaf_water_mass(i), 0.0_wp) * col_cohort%nplant(i)
         !----- The wood's radiation and sensible coefficient: its whole energy budget in surface_derivs   !
         !      (docs/science/vegetation_energy_dynamics.md §Leaf and wood energy balance). -------------!
         frozen%tissue%wai(i)         = col_cohort%wai(i)
         frozen%tissue%h_coeff_w(i)   = sensible_heat_coeff(pi * col_cohort%wai(i), aero%wood_gbh(i), rho, cp_air)
         frozen%tissue%abs_sw_wood(i) = forc%abs_sw_wood(i)
         frozen%tissue%abs_lw_wood(i) = forc%abs_lw_wood(i)
         !----- Tissue heat store a = cap/dt_fast, relaxing from the start-of-step temperature             !
         !      (docs/science/vegetation_energy_dynamics.md eqs. (3)-(5)): an algebraic closure in every  !
         !      stage evaluation, not a tableau unknown. -------------------------------------------------!
         frozen%tissue%leaf_hcap_per_dt(i)  = (frozen%tissue%leaf_dry_hcap(i) + frozen%tissue%leaf_wmass(i) * cp_liq) &
                                              / dt_fast
         frozen%tissue%wood_hcap_per_dt(i)  = (frozen%tissue%wood_dry_hcap(i) + frozen%tissue%wood_wmass(i) * cp_liq) &
                                              / dt_fast
         frozen%tissue%t_leaf0(i) = biophys%leaf_temp(i)
         frozen%tissue%t_wood0(i) = biophys%wood_temp(i)
         frozen%tissue%leaf_emiss(i) = col_cohort%leaf_emissivity(i)
         frozen%tissue%wood_emiss(i) = col_cohort%wood_emissivity(i)
         frozen%plant%nplant(i)   = col_cohort%nplant(i)
         !----- Film-evaporation conductances, from the pre-pass's boundary layers; none with the canopy  !
         !      water store off. ---------------------------------------------------------------------------!
         if (col_config%canopy_water_on) then
            frozen%film%g_film_leaf(i) = leaf_film_coeff(EFFAREA_EVAP, col_cohort%lai(i), aero%leaf_gbw(i))
            frozen%film%g_film_w(i) = leaf_film_coeff(EFFAREA_EVAP, col_cohort%wai(i), aero%wood_gbw(i))
         end if
      end do

      !----- The rest of the canopy-air boundary: the pre-pass's capacities and conductances, the air above, NEE. -!
      frozen%cas%cas_mass_capacity = cas_mass_capacity
      frozen%cas%cas_molar_capacity = cas_molar_capacity
      frozen%cas%g_atm_heat  = g_atm_heat
      frozen%cas%g_atm_vapour  = g_atm_vapour
      frozen%cas%g_atm_co2 = g_atm_co2
      !----- What a per-stage surface-layer re-solve needs besides the live CAS state. Displacement and !
      !      roughness are the pre-pass's own, so a re-solve at the state^n CAS reproduces it exactly. -!
      frozen%cas%aero_cfg     = col_config%aero
      frozen%cas%mo_u_ref     = aenv%u_ref
      frozen%cas%mo_zref      = aenv%zref
      frozen%cas%mo_displace  = aero%displace
      frozen%cas%mo_rough     = aero%rough
      frozen%cas%mo_theta_atm = aenv%theta_atm
      frozen%cas%mo_shv_atm   = aenv%shv_atm
      frozen%cas%mo_rho       = rho
      !----- Live inputs license the per-stage re-solve; a bundle not built here (a test fixture) leaves  !
      !      mo_live false, and its g_atm_heat/g_atm_vapour/g_atm_co2 are used as given. ----------------!
      frozen%cas%mo_live      = .true.
      frozen%cas%enthalpy_atm = forc%enthalpy_atm
      frozen%cas%shv_atm = forc%shv_atm
      frozen%cas%co2_atm = forc%co2_atm
      frozen%cas%nee_biotic = nee_biotic
      frozen%ground%abs_sw_ground = forc%abs_sw_ground
      frozen%ground%abs_lw_ground = forc%abs_lw_ground
      frozen%ground%ggnet = aero%ggnet
      frozen%cas%rho = rho
      frozen%cas%press = press

      frozen%cas%cas_condensation = col_config%integrator%cas_condensation      ! the per-stage condensation sink

      !----- The start-of-step surface solve, at the state^n canopy air: it gives the water network      !
      !      each cohort's transpiration demand (docs/science/soil_plant_water.md §4). ------------------!
      y_stage%cas_enthalpy = biophys%cas%can_enthalpy
      y_stage%cas_shv = biophys%cas%can_shv
      y_stage%cas_co2 = biophys%cas%can_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, t_ground, n, sf0)
      !----- Scale each film-evaporation conductance so the film's potential evaporation over dt_fast  !
      !      (this solve's) cannot exceed what the film holds plus this step's interception; the solve !
      !      below then sees the scaled value (docs/science/soil_plant_water.md §4). ------------------!
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
      !----- Credit the leaf the liquid enthalpy of the sapflow that replaces its transpiration, as the  !
      !      stages will, estimated from the first solve as if the sapflow matched it, and solve again  !
      !      (docs/science/soil_plant_water.md §4). --------------------------------------------------!
      do i = 1_ik, n
         t_up_wl = merge(biophys%wood_temp(i), biophys%leaf_temp(i), sf0%transp_c(i) >= 0.0_wp)
         frozen%tissue%qwflux_wl(i) = sf0%transp_c(i) * internal_energy_liquid(t_up_wl)
      end do
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow, t_ground, n, sf0)
      !----- The water reaching the ground: the meltwater under a pack, which took the rain and snow; the  !
      !      throughfall on bare ground; to either the shed water adds itself                               !
      !      (docs/science/soil_plant_water.md §5, "The ground inflow"). --------------------------------!
      if (snow_st%exists) then
         water_base = snow_st%melt_rate
      else
         water_base = throughfall_total
      end if
      precip_ground            = water_base + biophys%shed_water_rate
      hforc%snow_free_frac     = 1.0_wp - snow_st%snowfac
      hforc%t_ground           = t_ground
      hforc%q_air = qcas
      hforc%rho_air = rho
      !----- Phase of the evaporating surface, for the ice saturation branch                             !
      !      (docs/science/snow_biophysics.md eq. (1)). -------------------------------------------------!
      hforc%ground_fliq        = biophys%soil_energy%soil_fliq(1)
      !----- The temperature that values it: the meltwater's under a pack; on bare ground that of the rain  !
      !      and sub-threshold snowfall valued at the air temperature above the canopy; then mixed with the  !
      !      shed water's enthalpy (docs/science/soil_plant_water.md §5, "The ground inflow"). ------------!
      t_water_base = tcas
      if (snow_st%exists) then
         t_water_base = snow_st%t_melt
      else if (forc%rainfall + forc%snowfall > tiny_num) then
         t_water_base = liquid_energy_to_temp(precip_enthalpy(forc%rainfall, forc%snowfall, forc%air_temp) &
                                                / (forc%rainfall + forc%snowfall))
      end if
      inflow_enth = water_base * internal_energy_liquid(t_water_base) + biophys%shed_water_enth
      t_pond_inflow = t_water_base
      if (biophys%shed_water_rate > 0.0_wp) t_pond_inflow = liquid_energy_to_temp(inflow_enth / precip_ground)
      !----- The enthalpy it brings across the column's boundary: all of it on bare ground; under a pack  !
      !      only the shed water's, since the meltwater's comes from the pack, inside the column. --------!
      frozen%hydrology%ground_inflow_enth = inflow_enth
      if (snow_st%exists) frozen%hydrology%ground_inflow_enth = biophys%shed_water_enth
      !----- Bare-soil aerodynamic resistance; snow_free_frac above weights the flux, not this. -------!
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)

      !----- The soil-plant water network (docs/science/soil_plant_water.md): the pond, the soil layers and !
      !      every cohort's wood and leaves, solved together on the transpiration demand above. The march  !
      !      reads all its water from it and holds its flows fixed. ------------------------------------!
      call run_spac_solver(dt_fast, col_config, col_site, col_cohort, biophys, hforc, precip_ground, t_pond_inflow, &
                           root_frac(1:nsl, 1:n), sf0%transp_c(1:n), n, nsl, spac)

      !----- The network's work counters (docs/science/soil_plant_water.md §6). -----------------------!
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
      !----- Infiltration leaves the pond, so the top face is valued at the temperature that carries the  !
      !      enthalpy the pond gave up (docs/science/numerical_scheme.md §5b, "The pond holds heat"). ----!
      frozen%hydrology%t_infil = t_pond_inflow
      if (spac%infil > tiny_num) frozen%hydrology%t_infil = liquid_energy_to_temp(spac%infil_enth / spac%infil)
      !----- The canopy film's valuation temperature, a reference its store, interception and evaporation  !
      !      share: the rain's on bare ground; under a pack, where nothing is intercepted, tsupercool_liq    !
      !      (docs/science/soil_plant_water.md §5, "The ground inflow"). --------------------------------!
      frozen%hydrology%t_film_valuation = t_water_base
      if (snow_st%exists) frozen%hydrology%t_film_valuation = tsupercool_liq
      frozen%film%film_liquid_enthalpy = internal_energy_liquid(frozen%hydrology%t_film_valuation)
      frozen%hydrology%t_bot           = biophys%soil_energy%soil_temp(nsl)
      frozen%hydrology%w_surface1      = spac%water_soil(POND)
      frozen%hydrology%w_surface_energy1 = spac%e_pond
      frozen%hydrology%w_flux_frozen(1:nsl)   = 0.0_wp
      frozen%hydrology%w_flux_frozen(1:nsl-1) = spac%face(1:nsl-1) / (rho_h2o * dt_fast)
      allocate(frozen%hydrology%theta1(nsl))
      do k = 1_ik, nsl
         frozen%hydrology%theta1(k) = spac%water_soil(1_ik + k) / (rho_h2o * col_site%soil%dz(k))
      end do
      frozen%hydrology%theta_rate(1:nsl) = (frozen%hydrology%theta1(1:nsl) - biophys%soil_water%theta(1:nsl)) / dt_fast

      !----- Per-layer face closure: each layer's water change against the flows the soil energy will      !
      !      advect. It reads round-off, and checks that the faces handed on are the ones that moved the     !
      !      water; a network sub-step taken unconverged counts as a failed check. -------------------------!
      face_resid = 0.0_wp
      do k = 1_ik, nsl
         f_in = 0.0_wp
         f_out = 0.0_wp
         if (k > 1_ik)   f_in  = spac%face(k-1)
         if (k < nsl)    f_out = spac%face(k)
         if (k == nsl)   f_out = spac%drain
         if (k == 1_ik)  f_in  = spac%infil - spac%exfil - spac%evap
         face_resid = face_resid + abs((frozen%hydrology%theta1(k) - biophys%soil_water%theta(k))              &
                                       * col_site%soil%dz(k) * rho_h2o - (f_in - f_out - sum(spac%uptake(k, 1:n))))
      end do
      budget%soil_face_mass%resid   = face_resid
      budget%soil_face_mass%worst   = max(budget%soil_face_mass%worst, face_resid)
      budget%soil_face_mass%abs_sum = budget%soil_face_mass%abs_sum + face_resid
      budget%soil_face_mass%n_check = budget%soil_face_mass%n_check + 1_ik
      if (spac%nfail > 0_ik) budget%soil_face_mass%n_fail = budget%soil_face_mass%n_fail + 1_ik

      !----- Root uptake in total, by layer and by cohort, and its enthalpy at each layer's start-of-step  !
      !      temperature (docs/science/soil_plant_water.md §5), the number the wood is credited below. A   !
      !      patch whose roots took nothing places its zero sink by its root profile. --------------------!
      frozen%roots%uptake = sum(spac%uptake(1:nsl, 1:n)) / dt_fast
      do k = 1_ik, nsl
         u_liq_soil_k(k) = internal_energy_liquid(biophys%soil_energy%soil_temp(k))
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

      !----- Per-cohort hydraulic diagnostics, from the same network the physics commits: sapflow and   !
      !      root uptake per plant, and the wood's end-of-step potential and loss of conductance. -------!
      if (present(cdiag)) then
         do i = 1_ik, n
            cdiag(CO_SAPFLOW,     i) = frozen%plant%sapflow_frozen(i)
            cdiag(CO_ROOT_UPTAKE, i) = frozen%plant%uptake_frozen(i)
            cdiag(CO_PSI_WOOD,    i) = spac%psi_coh(WOOD, i)
            associate (hp => col_config%hydraulics_table%pft(col_cohort%pft(i)))
               cdiag(CO_PLC, i) = 1.0_wp - conductance_retained_fraction(spac%psi_coh(WOOD, i), hp%wood_psi50, hp%wood_kexp)
            end associate
         end do
      end if

      !----- Advective enthalpy of the plant water (docs/science/soil_plant_water.md eqs. (4)-(5)): the      !
      !      network's flows, with upwind temperatures frozen at state^n, since the stage temperatures are   !
      !      computed from these terms; each tissue's water change carries its start-of-step enthalpy. -----!
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
         frozen%tissue%wood_water_enth(i) = (uptake_gnd - sapflow_gnd) * u_liq_wood
         frozen%tissue%leaf_water_enth(i) = (sapflow_gnd - transp_gnd) * u_liq_leaf
      end do

      !----- Pack the prognostic state; plant water mass is itself prognostic, so it is copied directly. -!
      y%cas_enthalpy = biophys%cas%can_enthalpy
      y%cas_shv = biophys%cas%can_shv
      y%cas_co2 = biophys%cas%can_co2
      y%soil_energy(1:nsl) = biophys%soil_energy%soil_energy(1:nsl)
      y%theta(1:nsl)       = biophys%soil_water%theta(1:nsl)
      !----- The pond rides the state: the stages pass it through; each scheme commits the network's. -!
      y%w_surface          = biophys%soil_water%w_surface
      y%w_surface_energy     = biophys%soil_water%w_surface_energy
      y%leaf_water_mass(1:n) = biophys%leaf_water_mass(1:n)
      y%wood_water_mass(1:n) = biophys%wood_water_mass(1:n)
   end subroutine build_column_frozen

   !---------------------------------------------------------------------------------------!
   ! run_spac_solver -- one patch's soil and plant water over the fast step, as one network            !
   ! (meds_spac_solver; docs/science/soil_plant_water.md): builds the solver's inputs from the         !
   ! patch's columns (xylem conductance from the sapwood, rooting from height and fine-root carbon,    !
   ! root spacing from the patch's root-length density), starts it from the stores' state^n water and !
   ! advances it. e_demand is each cohort's transpiration demand; the network takes less where a leaf  !
   ! runs short (eq. (3)).                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine run_spac_solver(dt_fast, col_config, col_site, col_cohort, biophys, hforc, precip_ground, t_pond_inflow, &
                              root_frac, e_demand, n, nsl, spac)
      real(wp),               intent(in)    :: dt_fast
      type(column_config_t),  intent(in)    :: col_config
      type(column_site_t),    intent(in)    :: col_site
      type(column_cohort_t),  intent(in)    :: col_cohort
      type(patch_biophys_t),  intent(in)    :: biophys
      type(chydro_forcing_t), intent(in)    :: hforc
      real(wp),               intent(in)    :: precip_ground   !< [kg/m2/s] water reaching the pond
      real(wp),               intent(in)    :: t_pond_inflow   !< [K]       the temperature that values it
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
            spac%grav(i)     = crown_gravity_head(col_cohort%height(i))
            spac%k_xylem(i)  = col_cohort%nplant(i) * xylem_max_conductance(hp%wood_kmax, col_cohort%sap_area(i), &
                                                                           col_cohort%height(i), hp%vessel_curl)
            spac%e_demand(i) = e_demand(i)
            spac%r_close(i)  = hp%rwc_close
            root_len(1:nsl, i) = col_cohort%broot(i) * root_frac(1:nsl, i) * hp%specific_root_length
         end associate
      end do
      !----- Half the distance between a layer's roots, (pi L_v)^-1/2 over every cohort there, and each  !
      !      cohort's rhizosphere conductance in it (docs/science/plant_hydraulics.md eq. (14)). ---------!
      do k = 1_ik, nsl
         half_spacing(k) = 1.0_wp / sqrt(pi * max(sum(col_cohort%nplant(1:n) * root_len(k, 1:n))                &
                                                  / col_site%soil%dz(k), tiny_num))
         do i = 1_ik, n
            spac%g_rhizo(k, i) = rhizosphere_conductance(col_cohort%nplant(i), root_len(k, i), half_spacing(k),   &
                                                         col_config%hydraulics_table%pft(col_cohort%pft(i))%fine_root_radius)
         end do
      end do
      spac%rain       = precip_ground
      spac%t_rain     = t_pond_inflow
      spac%t_soil(1:nsl) = biophys%soil_energy%soil_temp(1:nsl)
      spac%w_pond_max = col_site%soil_water_opts%w_pond_max
      spac%rtol       = col_site%soil_water_opts%rtol
      spac%atol       = col_site%soil_water_opts%water_atol
      call spac_set_state(spac, col_site%soil, col_config%hydraulics_table, biophys%soil_water%theta(1:nsl),        &
                          biophys%soil_water%w_surface, biophys%soil_water%w_surface_energy,                              &
                          biophys%leaf_water_mass(1:n), biophys%wood_water_mass(1:n))
      call spac_advance(spac, col_site%soil, col_site%soil_water_opts, hforc, col_config%hydraulics_table,        &
                        dt_fast)
   end subroutine run_spac_solver

   !---------------------------------------------------------------------------------------!
   ! advance_snow_stage -- accumulate snowfall and rain-on-snow, advance the snow-surface energy      !
   ! balance at the state^n canopy air, and drain meltwater to the pond as a paired (mass, enthalpy) !
   ! transfer, reported as melt_rate and t_melt (docs/science/snow_biophysics.md §The snow store).    !
   ! Mutates only the pack; everything else it reports through st.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_snow_stage(snow, snow_params, snow_rough, dz_soil_top, abs_sw_ground, abs_lw_ground, &
                                 snowfall, rainfall, t_air, ggnet, t_soil_top, dt_fast, tcas, qcas,    &
                                 rho, press, st)
      type(snow_column_t),  intent(inout) :: snow           !< the pack store (the only state mutated here)
      type(snow_params_t),  intent(in)    :: snow_params
      real(wp),             intent(in)    :: snow_rough     !< [m]       snow-surface roughness ([aerodynamics])
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

      !----- Default: the bare-ground boundary (snowfac = 0). The stage always runs; with no pack and    !
      !      snowfall under min_new_snow_mass, snow_accumulate returns at once, and these defaults make  !
      !      surface_derivs' snow blend reduce exactly to bare ground. ---------------------------------!
      st%ground_rad = abs_sw_ground + abs_lw_ground
      st%swe0       = snow%swe(1)
      st%swe1  = snow%swe(1)
      st%enth0      = snow%snow_energy(1)
      st%enth1 = snow%snow_energy(1)
      snow_e0 = snow%snow_energy(1)
      call snow_accumulate(snow, snowfall, rainfall, t_air, dt_fast, snow_params)
      st%acc_enth = snow%snow_energy(1) - snow_e0   ! rainfall enthalpy into the pack (boundary in)
      st%exists   = snow%nlayer >= 1_ik             ! accumulate took snow+rain -> rainfall routing

      if (st%exists .and. snow%swe(1) > snow_params%tiny_snow_mass) then
         !----- snowfac is snow, 1 - snowfac bare soil. snow_energy_step scales the pack's exchange by   !
         !      snowfac, so a patchy pack barely exchanges and its fluxes return snowfac-weighted. ---!
         st%snowfac       = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), snow_params, snow_rough)
         senv%abs_sw      = abs_sw_ground
         senv%abs_lw = abs_lw_ground
         senv%can_temp    = tcas
         senv%can_shv = qcas
         senv%ggnet = ggnet
         senv%rho_air     = rho
         senv%press = press
         senv%t_soil_top  = t_soil_top
         senv%dz_soil_top = dz_soil_top
         call snow_energy_step(snow, senv, snow_params, dt_fast, st%snowfac, sfx)
         call snow_drain_meltwater(snow, snow_params, smelt)
         st%h_snow  = sfx%h_snow
         st%le_snow = sfx%le_snow
         st%g_base = sfx%g_base
         st%snowfac = sfx%snowfac                         ! the clamped fraction the kernel actually used
         !----- ground radiation boundary in = snow's snowfac-weighted net + bare's (1-snowfac) share !
         st%ground_rad = sfx%rnet + (1.0_wp - st%snowfac) * (abs_sw_ground + abs_lw_ground)
         st%subl_rate  = sfx%w_flux
         st%melt_rate  = (smelt%melt_mass + smelt%dump_mass) / dt_fast
         !----- Meltwater goes to the pond with its enthalpy, valued at t_melt, and nowhere else. Pack and !
         !      pond are both tracked stores, so the transfer cancels in the whole-column ledger         !
         !      (docs/science/numerical_scheme.md §5b, "The pond holds heat"). -------------------------!
         st%melt_enth = smelt%melt_enth + smelt%dump_enth
         st%t_melt    = t_3ple
         if (smelt%melt_mass + smelt%dump_mass > 0.0_wp)                                            &
            st%t_melt = liquid_energy_to_temp(st%melt_enth / (smelt%melt_mass + smelt%dump_mass))
      end if
      st%swe1  = snow%swe(1)
      st%enth1 = snow%snow_energy(1)
   end subroutine advance_snow_stage

end module meds_fast_frozen
