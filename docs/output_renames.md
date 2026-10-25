# Output variable names, and what they were

The diagnostic and restart files' variables are named `[mean_]<quantity>[_<breakdown>]_<level>`
([diagnostics](science/diagnostics.md) §3). This page maps each to its earlier name: one row per
variable, with every name it replaces (a fast `*_fast` twin merged into its daily variable, or a
dropped duplicate), its default tiers and its dimensions. Three variables go without a successor:
`gpp_rate_site` and `et_rate_site`, duplicates of `mean_gpp_py` and `mean_et_py`, and
`leaf_temp_cohort_fast`, the leaf temperature at each step's end: on the F tier the per-cohort leaf
temperature is `mean_leaf_temp_co`, the one each step's leaf solve used. Each state's `mean_` variable
also has a snapshot, its name without the prefix ([diagnostics](science/diagnostics.md) §4), which this
page does not list.

## structure

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `nplant_co` | `nplant_cohort` | D M | cohort | plant/m2 |
| `mean_dbh_co` | `dbh_cohort` | D M | cohort | cm |
| `mean_height_co` | `height_cohort`, `height_cohort_fast` | F D M | cohort | m |
| `mean_basal_area_co` | `basal_area_cohort` | M | cohort | cm2/plant |
| `mean_agb_co` | `agb_cohort` | D M | cohort | kgC/plant |
| `mean_leaf_area_co` | `leaf_area_cohort` | M | cohort | m2/plant |
| `mean_lai_co` | `lai_cohort` | M | cohort | m2/m2 |
| `mean_wai_co` | `wai_cohort` | M | cohort | m2/m2 |
| `mean_growth_avg_co` | `growth_avg_cohort` | M | cohort | cm/yr |
| `mean_overtopping_lai_co` | `overtopping_lai_cohort` | M | cohort | m2/m2 |
| `pft_co` | `pft_cohort` | D M | cohort | - |
| `owner_patch` | `owner_patch` | D M | cohort | - |
| `global_cohort_id` | `global_cohort_id` | D M | cohort | - |
| `mean_leaf_carbon_co` | `leaf_carbon_cohort` | M | cohort | kgC/plant |
| `mean_fineroot_carbon_co` | `fineroot_carbon_cohort` | M | cohort | kgC/plant |
| `mean_wood_carbon_co` | `wood_carbon_cohort` | M | cohort | kgC/plant |
| `mean_storage_carbon_co` | `storage_carbon_cohort` | M | cohort | kgC/plant |
| `mean_sla_co` | `sla_cohort` | M | cohort | m2/kgC |
| `mean_vcmax25_co` | `vcmax25_cohort` | M | cohort | umol/m2/s |
| `mean_rd25_co` | `rd25_cohort` | M | cohort | umol/m2/s |
| `mean_leaf_turnover_rate_co` | `leaf_turnover_rate_cohort` | M | cohort | 1/yr |
| `mean_leaf_flush_tendency_co` | `leaf_flush_tendency_cohort` | M | cohort | - |
| `mean_leaf_shed_tendency_co` | `leaf_shed_tendency_cohort` | M | cohort | - |
| `mean_dmax_psi_leaf_co` | `dmax_psi_leaf_cohort` | D M | cohort | MPa |
| `mean_dbh_growth_co` | `dbh_growth_cohort` | M | cohort | cm/yr |
| `mean_agb_growth_co` | `agb_growth_cohort` | M | cohort | kgC/plant/yr |
| `mean_mort_rate_co` | `mort_rate_cohort` | M | cohort | 1/yr |
| `area_pa` | `area_patch` | D M | patch | - |
| `age_pa` | `age_patch` | M | patch | yr |
| `dist_type_pa` | `dist_type_patch` | M | patch | - |
| `cohort_offset` | `cohort_offset` | D M | patch | - |
| `cohort_count` | `cohort_count` | D M | patch | - |
| `global_patch_id` | `global_patch_id` | D M | patch | - |
| `patch_site` | `patch_site` | D M | patch | - |
| `mean_nplant_pa` | `nplant_patch` | D M | patch | plant/m2 |
| `mean_lai_pa` | `lai_patch` | D M | patch | m2/m2 |
| `mean_agb_pa` | `agb_patch` | D M | patch | kgC/m2 |
| `mean_basal_area_pa` | `basal_area_patch` | M | patch | m2/m2 |
| `mean_nplant_py` | `nplant_site` | D M Y | polygon | plant/m2 |
| `mean_basal_area_py` | `basal_area_site` | D M Y | polygon | m2/m2 |
| `mean_agb_py` | `agb_site` | D M Y | polygon | kgC/m2 |
| `mean_lai_py` | `lai_site` | D M Y | polygon | m2/m2 |
| `mean_dbh_py` | `mean_dbh_site` | D M Y | polygon | cm |
| `mean_height_py` | `mean_height_site` | D M Y | polygon | m |
| `mean_canopy_height_py` | `canopy_height_site` | D M Y | polygon | m |
| `n_cohort_py` | `n_cohort_site` | D M Y | polygon | - |
| `n_patch_py` | `n_patch_site` | D M Y | polygon | - |
| `mean_nplant_pft_py` | `nplant_pft` | M Y | polygon × PFT | plant/m2 |
| `mean_agb_pft_py` | `agb_pft` | M Y | polygon × PFT | kgC/m2 |
| `mean_lai_pft_py` | `lai_pft` | M Y | polygon × PFT | m2/m2 |
| `mean_basal_area_pft_py` | `basal_area_pft` | M Y | polygon × PFT | m2/m2 |
| `mean_dbh_pft_py` | `mean_dbh_pft` | M Y | polygon × PFT | cm |
| `mean_nplant_size_py` | `nplant_size` | M Y | polygon × size class | plant/m2 |
| `mean_agb_size_py` | `agb_size` | M Y | polygon × size class | kgC/m2 |
| `mean_lai_size_py` | `lai_size` | M Y | polygon × size class | m2/m2 |
| `mean_basal_area_size_py` | `basal_area_size` | M Y | polygon × size class | m2/m2 |
| `mean_agb_growth_size_py` | `agb_growth_size` | M Y | polygon × size class | kgC/m2/yr |
| `mean_mort_rate_size_py` | `mort_rate_size` | M Y | polygon × size class | 1/yr |
| `mean_agb_growth_py` | `agb_growth_site` | M Y | polygon | kgC/m2/yr |
| `mean_agb_mort_py` | `agb_mort_site` | M Y | polygon | kgC/m2/yr |
| `mean_mort_rate_py` | `mort_rate_site` | M Y | polygon | 1/yr |
| `mean_nplant_recruit_py` | `nplant_recruit_site` | M Y | polygon | plant/m2/yr |
| `mean_disturb_area_py` | `disturb_area_site` | M Y | polygon | 1/yr |
| `mean_mort_carbon_background_py` | `mort_carbon_background_site` | M Y | polygon | kgC/m2/yr |
| `mean_mort_carbon_cull_py` | `mort_carbon_cull_site` | M Y | polygon | kgC/m2/yr |
| `mean_mort_carbon_disturb_py` | `mort_carbon_disturb_site` | M Y | polygon | kgC/m2/yr |
| `mean_mort_carbon_background_pa` | `mort_carbon_background_patch` | M | patch | kgC/m2/yr |
| `mean_agb_growth_pft_py` | `agb_growth_pft` | M Y | polygon × PFT | kgC/m2/yr |
| `mean_agb_mort_pft_py` | `agb_mort_pft` | M Y | polygon × PFT | kgC/m2/yr |

## carbon

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_gpp_py` | `gpp_rate_site`, `gpp_site`, `gpp_rate_fast` | F D M Y | polygon | umol/m2/s |
| `mean_npp_py` | `npp_site`, `npp_rate_fast` | F D M Y | polygon | umol/m2/s |
| `mean_leaf_resp_py` | `leaf_resp_site` | D M Y | polygon | umol/m2/s |
| `mean_stem_resp_py` | `stem_resp_site` | D M Y | polygon | umol/m2/s |
| `mean_root_resp_py` | `root_resp_site` | D M Y | polygon | umol/m2/s |
| `mean_leaf_carbon_py` | `leaf_carbon_site` | D M Y | polygon | kgC/m2 |
| `mean_fineroot_carbon_py` | `fineroot_carbon_site` | D M Y | polygon | kgC/m2 |
| `mean_wood_carbon_py` | `wood_carbon_site` | D M Y | polygon | kgC/m2 |
| `mean_storage_carbon_py` | `storage_carbon_site` | D M Y | polygon | kgC/m2 |
| `mean_bgb_py` | `bgb_site` | D M Y | polygon | kgC/m2 |
| `mean_veg_carbon_py` | `veg_carbon_site` | D M Y | polygon | kgC/m2 |
| `mean_gpp_pft_py` | `gpp_pft` | M Y | polygon × PFT | umol/m2/s |
| `mean_npp_pft_py` | `npp_pft` | M Y | polygon × PFT | umol/m2/s |
| `mean_gpp_pa` | `gpp_patch`, `gpp_rate_patch_fast` | F D M | patch | umol/m2/s |
| `mean_npp_pa` | `npp_patch`, `npp_rate_patch_fast` | F D M | patch | umol/m2/s |
| `mean_gpp_co` | `gpp_cohort`, `gpp_cohort_fast` | F D M | cohort | umol/plant/s |
| `mean_npp_co` | `npp_cohort` | D M | cohort | umol/plant/s |
| `mean_npp_leaf_py` | `npp_leaf_site` | M Y | polygon | kgC/m2/yr |
| `mean_npp_fineroot_py` | `npp_fineroot_site` | M Y | polygon | kgC/m2/yr |
| `mean_npp_wood_py` | `npp_wood_site` | M Y | polygon | kgC/m2/yr |
| `mean_npp_storage_py` | `npp_storage_site` | M Y | polygon | kgC/m2/yr |
| `mean_npp_repro_py` | `npp_repro_site` | M Y | polygon | kgC/m2/yr |
| `mean_growth_resp_py` | `growth_resp_site` | M Y | polygon | kgC/m2/yr |
| (retired: storage maintenance was removed) | `storage_resp_site` | M Y | polygon | kgC/m2/yr |
| `mean_root_exudate_py` | `root_exudate_site` | M Y | polygon | kgC/m2/yr |
| `mean_nee_py` | `nee_site`, `nee_fast` | F D M Y | polygon | umol/m2/s |
| `mean_nee_pa` | `nee_patch`, `nee_patch_fast` | F D M | patch | umol/m2/s |
| `mean_reco_py` | `reco_fast` | F | polygon | umol/m2/s |
| `mean_reco_pa` | `reco_patch_fast` | F | patch | umol/m2/s |
| `mean_cas_co2_pa` | `cas_co2_patch_fast` | F | patch | umol/mol |

## water

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_et_py` | `et_rate_site`, `et_site` | D M Y | polygon | kg/m2/s |
| `mean_runoff_py` | `runoff_site` | D M Y | polygon | kg/m2/s |
| `mean_runoff_surf_py` | `runoff_surf_site` | D M Y | polygon | kg/m2/s |
| `mean_drainage_py` | `drainage_site` | D M Y | polygon | kg/m2/s |
| `mean_soil_water_layer_si` | `soil_water_site`, `soil_water_site_fast` | F D M Y | site × soil layer | m3/m3 |
| `mean_soil_psi_layer_si` | `soil_psi_site` | D M Y | site × soil layer | MPa |
| `mean_soil_wetness_layer_si` | `soil_wetness_site` | D M Y | site × soil layer | - |
| `mean_w_surface_py` | `w_surface_site` | D M Y | polygon | kg/m2 |
| `mean_swe_py` | `swe_site` | D M Y | polygon | kg/m2 |
| `mean_snow_depth_py` | `snow_depth_site` | D M Y | polygon | m |
| `mean_w_surface_pa` | `w_surface_patch` | M | patch | kg/m2 |
| `mean_swe_pa` | `swe_patch` | M | patch | kg/m2 |
| `mean_soil_water_layer_pa` | `soil_water_layer_patch`, `soil_water_layer_patch_fast` | F D M | patch × soil layer | m3/m3 |
| `mean_soil_psi_layer_pa` | `soil_psi_layer_patch` | D M | patch × soil layer | MPa |
| `mean_soil_wetness_layer_pa` | `soil_wetness_layer_patch` | M | patch × soil layer | - |
| `mean_root_uptake_py` | `root_uptake_site` | D M Y | polygon | kg/m2/s |
| `mean_transp_leaf_py` | `transp_leaf_site` | D M Y | polygon | mol/m2 leaf/s |

## energy

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_cas_temp_py` | `cas_temp_site`, `cas_temp_fast` | F D M Y | polygon | K |
| `mean_cas_shv_py` | `cas_shv_site` | D M Y | polygon | kg/kg |
| `mean_cas_co2_py` | `cas_co2_site`, `cas_co2_fast` | F D M Y | polygon | umol/mol |
| `mean_cas_vpd_py` | `cas_vpd_site` | D M Y | polygon | Pa |
| `mean_cas_depth_py` | `cas_depth_site` | M Y | polygon | m |
| `mean_soil_temp_top_py` | `soil_temp_top_site`, `soil_temp_top_fast` | F D M Y | polygon | K |
| `mean_skin_temp_py` | `skin_temp_site`, `skin_temp_fast` | F D M Y | polygon | K |
| `mean_soil_temp_layer_si` | `soil_temp_site`, `soil_temp_site_fast` | F D M Y | site × soil layer | K |
| `mean_soil_fliq_layer_si` | `soil_fliq_site` | D M Y | site × soil layer | - |
| `mean_leaf_temp_py` | `leaf_temp_site` | D M Y | polygon | K |
| `mean_wood_temp_py` | `wood_temp_site` | D M Y | polygon | K |
| `mean_leaf_temp_co` | `leaf_temp_cohort`, `gx_leaf_temp_cohort_fast` | F D M | cohort | K |
| `mean_wood_temp_co` | `wood_temp_cohort` | D M | cohort | K |
| `mean_cas_temp_pa` | `cas_temp_patch`, `cas_temp_patch_fast` | F D M | patch | K |
| `mean_soil_temp_top_pa` | `soil_temp_top_patch`, `soil_temp_top_patch_fast` | F M | patch | K |
| `mean_skin_temp_pa` | `skin_temp_patch`, `skin_temp_patch_fast` | F D M | patch | K |
| `mean_soil_temp_layer_pa` | `soil_temp_layer_patch`, `soil_temp_layer_patch_fast` | F D M | patch × soil layer | K |
| `mean_soil_fliq_layer_pa` | `soil_fliq_layer_patch` | M | patch × soil layer | - |
| `mean_le_py` | `le_site`, `le_flux_fast` | F D M Y | polygon | W/m2 |
| `mean_h_py` | `h_site`, `h_flux_fast` | F D M Y | polygon | W/m2 |
| `mean_rnet_py` | `rnet_site`, `rnet_fast` | F D M Y | polygon | W/m2 |
| `mean_sw_in_vis_py` | `sw_in_vis_site` | D M Y | polygon | W/m2 |
| `mean_sw_in_nir_py` | `sw_in_nir_site` | D M Y | polygon | W/m2 |
| `mean_sw_up_vis_py` | `sw_up_vis_site` | D M Y | polygon | W/m2 |
| `mean_sw_up_nir_py` | `sw_up_nir_site` | D M Y | polygon | W/m2 |
| `mean_lw_up_py` | `lw_up_site`, `lw_up_fast` | F D M Y | polygon | W/m2 |
| `mean_ustar_py` | `ustar_site`, `ustar_fast` | F D M Y | polygon | m/s |
| `mean_ggnet_py` | `ggnet_site` | M | polygon | m/s |
| `mean_rough_py` | `rough_site` | M Y | polygon | m |
| `mean_displace_py` | `displace_site` | M Y | polygon | m |
| `mean_rough_pa` | `rough_patch` | D M | patch | m |
| `mean_displace_pa` | `displace_patch` | D M | patch | m |
| `mean_le_pa` | `le_patch`, `le_patch_fast` | F D M | patch | W/m2 |
| `mean_h_pa` | `h_patch`, `h_patch_fast` | F D M | patch | W/m2 |
| `mean_rnet_pa` | `rnet_patch`, `rnet_patch_fast` | F D M | patch | W/m2 |
| `mean_sw_up_py` | `sw_up_fast` | F | polygon | W/m2 |
| `mean_sw_up_pa` | `sw_up_patch_fast` | F | patch | W/m2 |
| `mean_lw_up_pa` | `lw_up_patch_fast` | F | patch | W/m2 |
| `mean_ustar_pa` | `ustar_patch_fast` | F | patch | m/s |

## biogeochem

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_soilc_fast_grnd_py` | `soilc_fast_grnd_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_fast_soil_py` | `soilc_fast_soil_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_struct_grnd_py` | `soilc_struct_grnd_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_struct_soil_py` | `soilc_struct_soil_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_microbial_py` | `soilc_microbial_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_slow_py` | `soilc_slow_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_passive_py` | `soilc_passive_site` | D M Y | polygon | kgC/m2 |
| `mean_soilc_total_py` | `soilc_total_site` | D M Y | polygon | kgC/m2 |
| `mean_rh_py` | `rh_site` | D M Y | polygon | umol/m2/s |
| `mean_soilc_total_pa` | `soilc_total_patch` | M | patch | kgC/m2 |
| `mean_rh_pa` | `rh_patch` | M | patch | umol/m2/s |
| `mean_litter_leaf_py` | `litter_leaf_site` | M Y | polygon | kgC/m2/yr |
| `mean_litter_fineroot_py` | `litter_fineroot_site` | M Y | polygon | kgC/m2/yr |
| `mean_litter_wood_py` | `litter_wood_site` | M Y | polygon | kgC/m2/yr |
| `mean_litter_repro_py` | `litter_repro_site` | M Y | polygon | kgC/m2/yr |
| `mean_litter_struct_py` | `litter_struct_site` | M Y | polygon | kgC/m2/yr |
| `mean_litter_leaf_pa` | `litter_leaf_patch` | M | patch | kgC/m2/yr |

## numerics

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_work_integ_steps_py` | `work_integ_steps_site` | D M Y | polygon | 1 |
| `mean_work_integ_rej_py` | `work_integ_rej_site` | D M Y | polygon | 1 |
| `mean_work_water_nsub_py` | `work_water_nsub_site` | D M Y | polygon | 1 |
| `mean_work_water_newton_py` | `work_water_newton_site` | D M Y | polygon | 1 |
| `mean_work_water_reject_py` | `work_water_reject_site` | D M Y | polygon | 1 |
| `mean_work_water_nfail_py` | `work_water_nfail_site` | D M Y | polygon | 1 |
| `mean_work_rk45_rescue_py` | `work_rk45_rescue_site` | D M Y | polygon | 1 |
| `mean_work_clamp_stage_py` | `work_clamp_stage_site` | D M Y | polygon | 1 |
| `mean_resid_energy_py` | `resid_energy_site` | D M Y | polygon | W/m2 |
| `mean_resid_water_py` | `resid_water_site` | D M Y | polygon | kg/m2/s |

## ecophys

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_anet_co` | `anet_cohort` | D M | cohort | umol/m2 leaf/s |
| `mean_agross_co` | `agross_cohort`, `gx_agross_cohort_fast` | F D M | cohort | umol/m2 leaf/s |
| `mean_gsw_co` | `gsw_cohort` | D M | cohort | mol/m2 leaf/s |
| `mean_gbw_co` | `gbw_cohort` | M | cohort | m/s |
| `mean_ci_co` | `ci_cohort` | D M | cohort | umol/mol |
| `mean_cs_co` | `cs_cohort` | M | cohort | umol/mol |
| `mean_rd_co` | `rd_cohort` | M | cohort | umol/m2 leaf/s |
| `mean_transp_co` | `transp_cohort` | D M | cohort | mol/m2 leaf/s |
| `mean_beta_stomata_co` | `beta_stomata_cohort` | D M | cohort | - |
| `mean_beta_nonstomata_co` | `beta_nonstomata_cohort` | D M | cohort | - |
| `mean_leaf_vpd_co` | `leaf_vpd_cohort`, `gx_vpd_cohort_fast` | F M | cohort | Pa |
| `mean_psi_leaf_co` | `psi_leaf_cohort`, `gx_psi_leaf_cohort_fast` | F D M | cohort | MPa |
| `mean_psi_wood_co` | `psi_wood_cohort` | D M | cohort | MPa |
| `mean_plc_co` | `plc_cohort` | D M | cohort | - |
| `mean_sapflow_co` | `sapflow_cohort` | D M | cohort | kg/plant/s |
| `mean_root_uptake_co` | `root_uptake_cohort` | D M | cohort | kg/plant/s |
| `mean_leaf_water_co` | `leaf_water_cohort` | M | cohort | kg/plant |
| `mean_wood_water_co` | `wood_water_cohort` | M | cohort | kg/plant |
| `mean_wind_co` | `wind_cohort` | M | cohort | m/s |
| `mean_anet_py` | `anet_site` | D M Y | polygon | umol/m2 leaf/s |
| `mean_gsw_py` | `gsw_site` | D M Y | polygon | mol/m2 leaf/s |
| `mean_psi_leaf_py` | `psi_leaf_site` | D M Y | polygon | MPa |
| `mean_plc_py` | `plc_site` | D M Y | polygon | - |
| `mean_beta_stomata_py` | `beta_stomata_site` | D M Y | polygon | - |
| `mean_leaf_par_co` | `gx_par_cohort_fast` | F | cohort | umol/m2/s |
| `mean_ca_co` | `gx_ca_cohort_fast` | F | cohort | umol/mol |
| `mean_pressure_co` | `gx_pressure_cohort_fast` | F | cohort | Pa |
| `mean_psi_predawn_co` | `gx_psi_predawn_cohort_fast` | F | cohort | MPa |
| `mean_gb_mol_co` | `gx_gb_cohort_fast` | F | cohort | mol/m2/s |

## radiation

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_abs_par_co` | `abs_par_cohort`, `abs_par_cohort_fast` | F D M | cohort | W/m2 |
| `mean_abs_sw_co` | `abs_sw_cohort` | D M | cohort | W/m2 |
| `mean_abs_lw_co` | `abs_lw_cohort` | M | cohort | W/m2 |
| `mean_sw_ground_py` | `sw_ground_site` | D M | polygon | W/m2 |
| `mean_lw_ground_py` | `lw_ground_site` | D M | polygon | W/m2 |

## forcing

| variable | was | default tiers | dimensions | units |
|---|---|---|---|---|
| `mean_sw_in_py` | `sw_in_site`, `sw_in_fast` | F D M Y | polygon | W/m2 |
| `mean_precip_py` | `precip_site`, `precip_fast` | F D M Y | polygon | kg/m2/s |
| `mean_air_temp_si` | `air_temp_by_site` | D M Y | site | K |
| `mean_qair_si` | `qair_by_site` | D M Y | site | kg/kg |
| `mean_psurf_si` | `psurf_by_site` | D M Y | site | Pa |
| `mean_rho_air_si` | `rho_air_by_site` | D M Y | site | kg/m3 |
| `mean_lwdown_si` | `lwdown_by_site` | D M Y | site | W/m2 |
| `mean_sw_in_si` | `sw_in_by_site` | D M Y | site | W/m2 |
| `mean_par_beam_si` | `par_beam_by_site` | D M Y | site | W/m2 |
| `mean_nir_beam_si` | `nir_beam_by_site` | D M Y | site | W/m2 |
| `mean_growth_temp_si` | `growth_temp_by_site` | D M Y | site | K |
| `mean_air_temp_py` | `air_temp_site`, `air_temp_fast` | F D M Y | polygon | K |
| `mean_qair_py` | `qair_site`, `qair_fast` | F D M Y | polygon | kg/kg |
| `mean_psurf_py` | `psurf_site`, `psurf_fast` | F D M Y | polygon | Pa |
| `mean_wind_py` | `wind_site`, `wind_fast` | F D M Y | polygon | m/s |
| `mean_lwdown_py` | `lwdown_site`, `lwdown_fast` | F D M Y | polygon | W/m2 |
| `mean_par_beam_py` | `par_beam_site`, `par_beam_fast` | F D M Y | polygon | W/m2 |
| `mean_par_diffuse_py` | `par_diffuse_site`, `par_diffuse_fast` | F D M Y | polygon | W/m2 |
| `mean_nir_beam_py` | `nir_beam_site`, `nir_beam_fast` | F D M Y | polygon | W/m2 |
| `mean_nir_diffuse_py` | `nir_diffuse_site`, `nir_diffuse_fast` | F D M Y | polygon | W/m2 |
| `mean_snowfall_py` | `snowfall_site`, `snowfall_fast` | F D M Y | polygon | kg/m2/s |
| `mean_atm_co2_py` | `atm_co2_site`, `atm_co2_fast` | F D M Y | polygon | umol/mol |
| `mean_cosz_py` | `cosz_site`, `cosz_fast` | F D M Y | polygon | 1 |
| `mean_rho_air_py` | `rho_air_site`, `rho_air_fast` | F D M Y | polygon | kg/m3 |
| `mean_wind_cas_top_pa` | `wind_cas_top_patch` | D M | patch | m/s |
| `mean_air_temp_cas_top_pa` | `air_temp_cas_top_patch` | D M | patch | K |
| `mean_cas_depth_pa` | `cas_depth_patch` | D M | patch | m |
| `mean_rainf_py` | `rainf_fast` | F | polygon | kg/m2/s |

## Per-site statics and coordinates (both files)

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

## The restart file

Every variable is the state at the file's instant, so none carries `mean_`. Unchanged: `n_cohort`,
`n_patch`, `n_site`, `n_window`, `meta_int`, `meta_real`, `global_cohort_id`, `owner_patch`,
`global_patch_id`, `owner_site`.

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
