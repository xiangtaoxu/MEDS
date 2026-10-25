# Output variable names, and what they were

This page is for readers of output written before the naming rule below, and for scripts that still use
the old names. It maps each current variable to the name or names it replaces.

The diagnostic and state files' variables are named `[mean_]<quantity>[_<breakdown>]_<level>`
([diagnostics](science/diagnostics.md) §3):

- `<level>` is `_co` (cohort), `_pa` (patch), `_si` (site) or `_py` (polygon);
- `<breakdown>` is `_pft`, `_size` (DBH class) or `_layer` (soil layer);
- `mean_` marks a mean over the record's period. A name without it is the value at the period's end,
  its snapshot ([diagnostics](science/diagnostics.md) §4); this page lists the means only.

The tables are grouped as `meds_main --dump-io-config` groups the variables. A row's second column holds
every name the variable replaces, including a sub-daily `*_fast` twin merged into its daily variable.
The tiers a variable can be written at, its dimensions and its units are in `meds_io_config.toml`, which
that command generates, and in [diagnostics](science/diagnostics.md) §5.

Three old variables have no successor: `gpp_rate_site` and `et_rate_site`, duplicates of `mean_gpp_py`
and `mean_et_py`, and `leaf_temp_cohort_fast`, the leaf temperature at each step's end. On the F tier
the per-cohort leaf temperature is `mean_leaf_temp_co`, the one each step's leaf solve used.

## structure

| variable | was |
|---|---|
| `nplant_co` | `nplant_cohort` |
| `mean_dbh_co` | `dbh_cohort` |
| `mean_height_co` | `height_cohort`, `height_cohort_fast` |
| `mean_basal_area_co` | `basal_area_cohort` |
| `mean_agb_co` | `agb_cohort` |
| `mean_leaf_area_co` | `leaf_area_cohort` |
| `mean_lai_co` | `lai_cohort` |
| `mean_wai_co` | `wai_cohort` |
| `mean_growth_avg_co` | `growth_avg_cohort` |
| `mean_overtopping_lai_co` | `overtopping_lai_cohort` |
| `pft_co` | `pft_cohort` |
| `owner_patch` | `owner_patch` |
| `global_cohort_id` | `global_cohort_id` |
| `mean_leaf_carbon_co` | `leaf_carbon_cohort` |
| `mean_fineroot_carbon_co` | `fineroot_carbon_cohort` |
| `mean_wood_carbon_co` | `wood_carbon_cohort` |
| `mean_storage_carbon_co` | `storage_carbon_cohort` |
| `mean_sla_co` | `sla_cohort` |
| `mean_vcmax25_co` | `vcmax25_cohort` |
| `mean_rd25_co` | `rd25_cohort` |
| `mean_leaf_turnover_rate_co` | `leaf_turnover_rate_cohort` |
| `mean_leaf_flush_tendency_co` | `leaf_flush_tendency_cohort` |
| `mean_leaf_shed_tendency_co` | `leaf_shed_tendency_cohort` |
| `mean_dmax_psi_leaf_co` | `dmax_psi_leaf_cohort` |
| `mean_dbh_growth_co` | `dbh_growth_cohort` |
| `mean_agb_growth_co` | `agb_growth_cohort` |
| `mean_mort_rate_co` | `mort_rate_cohort` |
| `area_pa` | `area_patch` |
| `age_pa` | `age_patch` |
| `dist_type_pa` | `dist_type_patch` |
| `cohort_offset` | `cohort_offset` |
| `cohort_count` | `cohort_count` |
| `global_patch_id` | `global_patch_id` |
| `patch_site` | `patch_site` |
| `mean_nplant_pa` | `nplant_patch` |
| `mean_lai_pa` | `lai_patch` |
| `mean_agb_pa` | `agb_patch` |
| `mean_basal_area_pa` | `basal_area_patch` |
| `mean_nplant_py` | `nplant_site` |
| `mean_basal_area_py` | `basal_area_site` |
| `mean_agb_py` | `agb_site` |
| `mean_lai_py` | `lai_site` |
| `mean_dbh_py` | `mean_dbh_site` |
| `mean_height_py` | `mean_height_site` |
| `mean_canopy_height_py` | `canopy_height_site` |
| `n_cohort_py` | `n_cohort_site` |
| `n_patch_py` | `n_patch_site` |
| `mean_nplant_pft_py` | `nplant_pft` |
| `mean_agb_pft_py` | `agb_pft` |
| `mean_lai_pft_py` | `lai_pft` |
| `mean_basal_area_pft_py` | `basal_area_pft` |
| `mean_dbh_pft_py` | `mean_dbh_pft` |
| `mean_nplant_size_py` | `nplant_size` |
| `mean_agb_size_py` | `agb_size` |
| `mean_lai_size_py` | `lai_size` |
| `mean_basal_area_size_py` | `basal_area_size` |
| `mean_agb_growth_size_py` | `agb_growth_size` |
| `mean_mort_rate_size_py` | `mort_rate_size` |
| `mean_agb_growth_py` | `agb_growth_site` |
| `mean_agb_mort_py` | `agb_mort_site` |
| `mean_mort_rate_py` | `mort_rate_site` |
| `mean_nplant_recruit_py` | `nplant_recruit_site` |
| `mean_disturb_area_py` | `disturb_area_site` |
| `mean_mort_carbon_background_py` | `mort_carbon_background_site` |
| `mean_mort_carbon_cull_py` | `mort_carbon_cull_site` |
| `mean_mort_carbon_disturb_py` | `mort_carbon_disturb_site` |
| `mean_mort_carbon_background_pa` | `mort_carbon_background_patch` |
| `mean_agb_growth_pft_py` | `agb_growth_pft` |
| `mean_agb_mort_pft_py` | `agb_mort_pft` |

## carbon

| variable | was |
|---|---|
| `mean_gpp_py` | `gpp_rate_site`, `gpp_site`, `gpp_rate_fast` |
| `mean_npp_py` | `npp_site`, `npp_rate_fast` |
| `mean_leaf_resp_py` | `leaf_resp_site` |
| `mean_stem_resp_py` | `stem_resp_site` |
| `mean_root_resp_py` | `root_resp_site` |
| `mean_leaf_carbon_py` | `leaf_carbon_site` |
| `mean_fineroot_carbon_py` | `fineroot_carbon_site` |
| `mean_wood_carbon_py` | `wood_carbon_site` |
| `mean_storage_carbon_py` | `storage_carbon_site` |
| `mean_bgb_py` | `bgb_site` |
| `mean_veg_carbon_py` | `veg_carbon_site` |
| `mean_gpp_pft_py` | `gpp_pft` |
| `mean_npp_pft_py` | `npp_pft` |
| `mean_gpp_pa` | `gpp_patch`, `gpp_rate_patch_fast` |
| `mean_npp_pa` | `npp_patch`, `npp_rate_patch_fast` |
| `mean_gpp_co` | `gpp_cohort`, `gpp_cohort_fast` |
| `mean_npp_co` | `npp_cohort` |
| `mean_npp_leaf_py` | `npp_leaf_site` |
| `mean_npp_fineroot_py` | `npp_fineroot_site` |
| `mean_npp_wood_py` | `npp_wood_site` |
| `mean_npp_storage_py` | `npp_storage_site` |
| `mean_npp_repro_py` | `npp_repro_site` |
| `mean_growth_resp_py` | `growth_resp_site` |
| (retired: storage maintenance was removed) | `storage_resp_site` |
| `mean_root_exudate_py` | `root_exudate_site` |
| `mean_nee_py` | `nee_site`, `nee_fast` |
| `mean_nee_pa` | `nee_patch`, `nee_patch_fast` |
| `mean_reco_py` | `reco_fast` |
| `mean_reco_pa` | `reco_patch_fast` |
| `mean_cas_co2_pa` | `cas_co2_patch_fast` |

## water

| variable | was |
|---|---|
| `mean_et_py` | `et_rate_site`, `et_site` |
| `mean_runoff_py` | `runoff_site` |
| `mean_runoff_surf_py` | `runoff_surf_site` |
| `mean_drainage_py` | `drainage_site` |
| `mean_soil_water_layer_si` | `soil_water_site`, `soil_water_site_fast` |
| `mean_soil_psi_layer_si` | `soil_psi_site` |
| `mean_soil_wetness_layer_si` | `soil_wetness_site` |
| `mean_w_surface_py` | `w_surface_site` |
| `mean_swe_py` | `swe_site` |
| `mean_snow_depth_py` | `snow_depth_site` |
| `mean_w_surface_pa` | `w_surface_patch` |
| `mean_swe_pa` | `swe_patch` |
| `mean_soil_water_layer_pa` | `soil_water_layer_patch`, `soil_water_layer_patch_fast` |
| `mean_soil_psi_layer_pa` | `soil_psi_layer_patch` |
| `mean_soil_wetness_layer_pa` | `soil_wetness_layer_patch` |
| `mean_root_uptake_py` | `root_uptake_site` |
| `mean_transp_leaf_py` | `transp_leaf_site` |

## energy

| variable | was |
|---|---|
| `mean_cas_temp_py` | `cas_temp_site`, `cas_temp_fast` |
| `mean_cas_shv_py` | `cas_shv_site` |
| `mean_cas_co2_py` | `cas_co2_site`, `cas_co2_fast` |
| `mean_cas_vpd_py` | `cas_vpd_site` |
| `mean_cas_depth_py` | `cas_depth_site` |
| `mean_soil_temp_top_py` | `soil_temp_top_site`, `soil_temp_top_fast` |
| `mean_skin_temp_py` | `skin_temp_site`, `skin_temp_fast` |
| `mean_soil_temp_layer_si` | `soil_temp_site`, `soil_temp_site_fast` |
| `mean_soil_fliq_layer_si` | `soil_fliq_site` |
| `mean_leaf_temp_py` | `leaf_temp_site` |
| `mean_wood_temp_py` | `wood_temp_site` |
| `mean_leaf_temp_co` | `leaf_temp_cohort`, `gx_leaf_temp_cohort_fast` |
| `mean_wood_temp_co` | `wood_temp_cohort` |
| `mean_cas_temp_pa` | `cas_temp_patch`, `cas_temp_patch_fast` |
| `mean_soil_temp_top_pa` | `soil_temp_top_patch`, `soil_temp_top_patch_fast` |
| `mean_skin_temp_pa` | `skin_temp_patch`, `skin_temp_patch_fast` |
| `mean_soil_temp_layer_pa` | `soil_temp_layer_patch`, `soil_temp_layer_patch_fast` |
| `mean_soil_fliq_layer_pa` | `soil_fliq_layer_patch` |
| `mean_le_py` | `le_site`, `le_flux_fast` |
| `mean_h_py` | `h_site`, `h_flux_fast` |
| `mean_rnet_py` | `rnet_site`, `rnet_fast` |
| `mean_sw_in_vis_py` | `sw_in_vis_site` |
| `mean_sw_in_nir_py` | `sw_in_nir_site` |
| `mean_sw_up_vis_py` | `sw_up_vis_site` |
| `mean_sw_up_nir_py` | `sw_up_nir_site` |
| `mean_lw_up_py` | `lw_up_site`, `lw_up_fast` |
| `mean_ustar_py` | `ustar_site`, `ustar_fast` |
| `mean_ggnet_py` | `ggnet_site` |
| `mean_rough_py` | `rough_site` |
| `mean_displace_py` | `displace_site` |
| `mean_rough_pa` | `rough_patch` |
| `mean_displace_pa` | `displace_patch` |
| `mean_le_pa` | `le_patch`, `le_patch_fast` |
| `mean_h_pa` | `h_patch`, `h_patch_fast` |
| `mean_rnet_pa` | `rnet_patch`, `rnet_patch_fast` |
| `mean_sw_up_py` | `sw_up_fast` |
| `mean_sw_up_pa` | `sw_up_patch_fast` |
| `mean_lw_up_pa` | `lw_up_patch_fast` |
| `mean_ustar_pa` | `ustar_patch_fast` |

## biogeochem

| variable | was |
|---|---|
| `mean_soilc_fast_grnd_py` | `soilc_fast_grnd_site` |
| `mean_soilc_fast_soil_py` | `soilc_fast_soil_site` |
| `mean_soilc_struct_grnd_py` | `soilc_struct_grnd_site` |
| `mean_soilc_struct_soil_py` | `soilc_struct_soil_site` |
| `mean_soilc_microbial_py` | `soilc_microbial_site` |
| `mean_soilc_slow_py` | `soilc_slow_site` |
| `mean_soilc_passive_py` | `soilc_passive_site` |
| `mean_soilc_total_py` | `soilc_total_site` |
| `mean_rh_py` | `rh_site` |
| `mean_soilc_total_pa` | `soilc_total_patch` |
| `mean_rh_pa` | `rh_patch` |
| `mean_litter_leaf_py` | `litter_leaf_site` |
| `mean_litter_fineroot_py` | `litter_fineroot_site` |
| `mean_litter_wood_py` | `litter_wood_site` |
| `mean_litter_repro_py` | `litter_repro_site` |
| `mean_litter_struct_py` | `litter_struct_site` |
| `mean_litter_leaf_pa` | `litter_leaf_patch` |

## numerics

| variable | was |
|---|---|
| `mean_work_integ_steps_py` | `work_integ_steps_site` |
| `mean_work_integ_rej_py` | `work_integ_rej_site` |
| `mean_work_water_nsub_py` | `work_water_nsub_site` |
| `mean_work_water_newton_py` | `work_water_newton_site` |
| `mean_work_water_reject_py` | `work_water_reject_site` |
| `mean_work_water_nfail_py` | `work_water_nfail_site` |
| `mean_work_rk45_rescue_py` | `work_rk45_rescue_site` |
| `mean_work_clamp_stage_py` | `work_clamp_stage_site` |
| `mean_resid_energy_py` | `resid_energy_site` |
| `mean_resid_water_py` | `resid_water_site` |

## ecophys

| variable | was |
|---|---|
| `mean_anet_co` | `anet_cohort` |
| `mean_agross_co` | `agross_cohort`, `gx_agross_cohort_fast` |
| `mean_gsw_co` | `gsw_cohort` |
| `mean_gbw_co` | `gbw_cohort` |
| `mean_ci_co` | `ci_cohort` |
| `mean_cs_co` | `cs_cohort` |
| `mean_rd_co` | `rd_cohort` |
| `mean_transp_co` | `transp_cohort` |
| `mean_beta_stomata_co` | `beta_stomata_cohort` |
| `mean_beta_nonstomata_co` | `beta_nonstomata_cohort` |
| `mean_leaf_vpd_co` | `leaf_vpd_cohort`, `gx_vpd_cohort_fast` |
| `mean_psi_leaf_co` | `psi_leaf_cohort`, `gx_psi_leaf_cohort_fast` |
| `mean_psi_wood_co` | `psi_wood_cohort` |
| `mean_plc_co` | `plc_cohort` |
| `mean_sapflow_co` | `sapflow_cohort` |
| `mean_root_uptake_co` | `root_uptake_cohort` |
| `mean_leaf_water_co` | `leaf_water_cohort` |
| `mean_wood_water_co` | `wood_water_cohort` |
| `mean_wind_co` | `wind_cohort` |
| `mean_anet_py` | `anet_site` |
| `mean_gsw_py` | `gsw_site` |
| `mean_psi_leaf_py` | `psi_leaf_site` |
| `mean_plc_py` | `plc_site` |
| `mean_beta_stomata_py` | `beta_stomata_site` |
| `mean_leaf_par_co` | `gx_par_cohort_fast` |
| `mean_ca_co` | `gx_ca_cohort_fast` |
| `mean_pressure_co` | `gx_pressure_cohort_fast` |
| `mean_psi_predawn_co` | `gx_psi_predawn_cohort_fast` |
| `mean_gb_mol_co` | `gx_gb_cohort_fast` |

## radiation

| variable | was |
|---|---|
| `mean_abs_par_co` | `abs_par_cohort`, `abs_par_cohort_fast` |
| `mean_abs_sw_co` | `abs_sw_cohort` |
| `mean_abs_lw_co` | `abs_lw_cohort` |
| `mean_sw_ground_py` | `sw_ground_site` |
| `mean_lw_ground_py` | `lw_ground_site` |

## forcing

| variable | was |
|---|---|
| `mean_sw_in_py` | `sw_in_site`, `sw_in_fast` |
| `mean_precip_py` | `precip_site`, `precip_fast` |
| `mean_air_temp_si` | `air_temp_by_site` |
| `mean_qair_si` | `qair_by_site` |
| `mean_psurf_si` | `psurf_by_site` |
| `mean_rho_air_si` | `rho_air_by_site` |
| `mean_lwdown_si` | `lwdown_by_site` |
| `mean_sw_in_si` | `sw_in_by_site` |
| `mean_par_beam_si` | `par_beam_by_site` |
| `mean_nir_beam_si` | `nir_beam_by_site` |
| `mean_growth_temp_si` | `growth_temp_by_site` |
| `mean_air_temp_py` | `air_temp_site`, `air_temp_fast` |
| `mean_qair_py` | `qair_site`, `qair_fast` |
| `mean_psurf_py` | `psurf_site`, `psurf_fast` |
| `mean_wind_py` | `wind_site`, `wind_fast` |
| `mean_lwdown_py` | `lwdown_site`, `lwdown_fast` |
| `mean_par_beam_py` | `par_beam_site`, `par_beam_fast` |
| `mean_par_diffuse_py` | `par_diffuse_site`, `par_diffuse_fast` |
| `mean_nir_beam_py` | `nir_beam_site`, `nir_beam_fast` |
| `mean_nir_diffuse_py` | `nir_diffuse_site`, `nir_diffuse_fast` |
| `mean_snowfall_py` | `snowfall_site`, `snowfall_fast` |
| `mean_atm_co2_py` | `atm_co2_site`, `atm_co2_fast` |
| `mean_cosz_py` | `cosz_site`, `cosz_fast` |
| `mean_rho_air_py` | `rho_air_site`, `rho_air_fast` |
| `mean_wind_cas_top_pa` | `wind_cas_top_patch` |
| `mean_air_temp_cas_top_pa` | `air_temp_cas_top_patch` |
| `mean_cas_depth_pa` | `cas_depth_patch` |
| `mean_rainf_py` | `rainf_fast` |

## Per-site statics and coordinates

Both file kinds hold these, except `chave_e_py`, which only the diagnostic files hold, and only with a
`chave2014` height allometry.

| variable | was |
|---|---|
| `area_si` | `site_area` |
| `aspect_si` | `site_aspect` |
| `bottom_bc_si` | `site_bottom_bc` |
| `deep_temp_si` | `site_deep_temp` |
| `drain_factor_si` | `site_drain_factor` |
| `elevation_si` | `site_elevation` |
| `n_soil_si` | `site_n_soil` |
| `slope_si` | `site_slope` |
| `soil_depth_si` | `site_soil_depth` |
| `chave_e_py` | `chave_e` |

## The state file

Every variable is the state at the file's instant, so none carries `mean_`. Unchanged: `n_cohort`,
`n_patch`, `n_site`, `n_window`, `meta_int`, `meta_real`, `global_cohort_id`, `owner_patch`,
`global_patch_id`, `owner_site`. `kt_last_day_py` has no earlier name.

| variable | was |
|---|---|
| `adapt_dt_last_pa` | `adapt_dt_last` |
| `agb_co` | `agb` |
| `age_pa` | `patch_age` |
| `area_pa` | `patch_area` |
| `basal_area_co` | `basal_area` |
| `cas_can_co2_pa` | `cas_can_co2` |
| `cas_can_enthalpy_pa` | `cas_can_enthalpy` |
| `cas_can_shv_pa` | `cas_can_shv` |
| `cas_can_temp_pa` | `cas_can_temp` |
| `cold_degree_days_co` | `cold_degree_days` |
| `dbh_co` | `dbh` |
| `dist_type_pa` | `dist_type` |
| `dmax_psi_leaf_accum_co` | `dmax_psi_leaf_accum` |
| `dmax_psi_leaf_co` | `dmax_psi_leaf` |
| `dry_psi_sum_co` | `dry_psi_sum` |
| `fineroot_carbon_co` | `fineroot_carbon` |
| `growing_degree_days_co` | `growing_degree_days` |
| `growth_accum_co` | `growth_accum` |
| `growth_avg_co` | `growth_avg` |
| `growth_count_co` | `growth_count` |
| `growth_hist_co` | `growth_hist` |
| `growth_hist_pos_py` | `growth_hist_pos` |
| `height_co` | `height` |
| `leaf_area_co` | `leaf_area` |
| `leaf_carbon_co` | `leaf_carbon` |
| `leaf_energy_co` | `leaf_energy` |
| `leaf_flush_tendency_co` | `leaf_flush_tendency` |
| `leaf_shed_tendency_co` | `leaf_shed_tendency` |
| `leaf_surf_water_co` | `leaf_surf_water` |
| `leaf_temp_co` | `leaf_temp` |
| `leaf_turnover_rate_co` | `leaf_turnover_rate` |
| `leaf_water_mass_co` | `leaf_water_mass` |
| `light_hours_avg_co` | `light_hours_mean` |
| `nonstructural_carbon_co` | `nonstructural_carbon` |
| `nplant_co` | `nplant` |
| `overtopping_lai_co` | `overtopping_lai` |
| `pft_co` | `pft` |
| `rd25_co` | `rd25` |
| `recruit_pool_pa` | `recruit_pool` |
| `shed_water_enth_pa` | `shed_water_enth` |
| `shed_water_rate_pa` | `shed_water_rate` |
| `sla_co` | `sla` |
| `slow_co2_rate_pa` | `slow_co2_rate` |
| `snow_depth_layer_pa` | `snow_depth` |
| `snow_energy_layer_pa` | `snow_energy` |
| `snow_fliq_layer_pa` | `snow_fliq` |
| `snow_nlayer_pa` | `snow_nlayer` |
| `snow_swe_layer_pa` | `snow_swe` |
| `snow_temp_layer_pa` | `snow_temp` |
| `soil_energy_layer_pa` | `soil_energy` |
| `soil_fliq_layer_pa` | `soil_fliq` |
| `soil_temp_layer_pa` | `soil_temp` |
| `soil_theta_layer_pa` | `soil_theta` |
| `soil_w_surface_enth_pa` | `soil_w_surface_enth` |
| `soil_w_surface_pa` | `soil_w_surface` |
| `soilc_fast_grnd_pa` | `soilc_fast_grnd` |
| `soilc_fast_soil_pa` | `soilc_fast_soil` |
| `soilc_lignin_grnd_pa` | `soilc_lignin_grnd` |
| `soilc_lignin_soil_pa` | `soilc_lignin_soil` |
| `soilc_microbial_pa` | `soilc_microbial` |
| `soilc_passive_pa` | `soilc_passive` |
| `soilc_slow_pa` | `soilc_slow` |
| `soilc_struct_grnd_pa` | `soilc_struct_grnd` |
| `soilc_struct_soil_pa` | `soilc_struct_soil` |
| `t_growth_avg_si` | `t_growth_avg` |
| `vcmax25_co` | `vcmax25` |
| `wet_psi_sum_co` | `wet_psi_sum` |
| `wood_carbon_co` | `wood_carbon` |
| `wood_energy_co` | `wood_energy` |
| `wood_surf_water_co` | `wood_surf_water` |
| `wood_temp_co` | `wood_temp` |
| `wood_water_mass_co` | `wood_water_mass` |
