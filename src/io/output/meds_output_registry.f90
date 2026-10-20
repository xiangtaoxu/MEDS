! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_registry -- the ONE source-defined registration list of diagnostic variables, and !
! the config resolution that turns it into per-tier live-variable index lists                   !
! (docs/science/diagnostics.md §5, §8).                                                         !
!                                                                                          !
! netCDF-FREE (links types + config + demography state + the integrate SRC_* codes). Adding a      !
! variable is TWO coupled edits: one add_variable() line here + one case in extract_variable        !
! (meds_output_integrate). Resolution order (docs/science/diagnostics.md §5): source defaults        !
! (checked by the annual and fast guards) -> axis toggles -> group toggles -> per-tier enables ->    !
! freq index; the per-variable overrides (meds_io_config.toml) are applied separately, after.       !
! Design record: MEDS_IO_DESIGN.md §3, §6.                                                          !
!==========================================================================================!
module meds_output_registry
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_forcing_opts, only : MET_BACKEND_ED_ERA5LAND
   use meds_polygon_state_vars, only : polygon_t
   use meds_diag_vars,          only : N_CO_FAST_ROWS, N_CO_ROWS, N_PA_ROWS, alloc_cohort_diag, alloc_patch_diag, &
                                       N_SI_ROWS, alloc_site_diag, SI_SW_IN, SI_PAR_BEAM, SI_NIR_BEAM, SI_LWDOWN, &
                                       SI_TAIR, SI_QAIR, SI_PSURF, SI_RHO_AIR
   use meds_output_opts,    only : output_config_t, N_FREQ, N_GRP, N_AXIS,                       &
                                   N_DBH_CLASS_DEFAULT, DBH_EDGES_DEFAULT,                       &
                                   FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE,   &
                                   GRP_STRUCTURE, GRP_CARBON, GRP_WATER, GRP_ENERGY,              &
                                   GRP_RADIATION, GRP_ECOPHYS, GRP_BIOGEOCHEM, GRP_NUMERICS,      &
                                   GRP_FORCING
   use meds_output_types,   only : var_desc_t, output_registry_t, output_files_t, output_buffers_t,              &
                                   MAX_OUTPUT_VARS, MAX_DBH_CLASS,                                &
                                   AGG_MEAN, AGG_LAST, AGG_TMEAN, AGG_SUM, DIM_SCALAR, DIM_COHORT,&
                                   DIM_PATCH, DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, &
                                   XTYPE_DOUBLE, XTYPE_INT, MISSING_VALUE, MISSING_INT, ragged_dim
   use meds_diagnostic_reduce, only : W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB
   use meds_constants,         only : cm2_to_m2, n_soil_layer_max
   use meds_output_integrate, only : alloc_integ_buffer, src_class, SRCK_FAST,                                         &
        FLD_C_NPLANT, FLD_C_DBH, FLD_C_HEIGHT, FLD_C_BASAL_AREA, FLD_C_AGB, FLD_C_LEAF_AREA,     &
        FLD_C_GROWTH_AVG, FLD_C_PFT, FLD_C_OWNER_PATCH, FLD_C_GLOBAL_ID, FLD_C_LAI, FLD_C_WAI,   &
        FLD_C_LEAF_CARBON, FLD_C_FINEROOT_CARBON, FLD_C_WOOD_CARBON, FLD_C_STORAGE_CARBON,       &
        FLD_C_BGB, FLD_C_VEG_CARBON, FLD_C_SLA, FLD_C_VCMAX25, FLD_C_RD25, FLD_C_LEAF_TURNOVER,  &
        FLD_C_OVERTOP_LAI, FLD_C_GPP_ACCUM, FLD_C_NPP_ACCUM, FLD_C_LEAF_RESP, FLD_C_STEM_RESP,   &
        FLD_C_ROOT_RESP, FLD_C_DMAX_PSI_LEAF, FLD_C_PHENO_FLUSH, FLD_C_PHENO_SHED,               &
        FLD_P_AREA, FLD_P_AGE, FLD_P_DIST_TYPE, FLD_P_COHORT_OFFSET, FLD_P_COHORT_COUNT,         &
        FLD_P_GLOBAL_ID, FLD_P_SITE, FLD_P_CAS_DEPTH, FLD_P_SWE, FLD_P_SNOW_DEPTH,               &
        FLD_P_SOILC_FAST_GRND, FLD_P_SOILC_FAST_SOIL, FLD_P_SOILC_STRUCT_GRND,                   &
        FLD_P_SOILC_STRUCT_SOIL, FLD_P_SOILC_MICROBIAL, FLD_P_SOILC_SLOW,                        &
        FLD_P_SOILC_PASSIVE, FLD_P_SOILC_TOTAL, FLD_P_RH,                                        &
        FLD_L_SOIL_TEMP, FLD_L_SOIL_WATER, FLD_L_SOIL_PSI, FLD_L_SOIL_WETNESS, FLD_L_SOIL_FLIQ,  &
        SRC_S_ET, SRC_S_N_COHORT, SRC_S_N_PATCH, SRC_S_CANOPY_HEIGHT,                            &
        SRC_S_RUNOFF, SRC_S_RUNOFF_SURF, SRC_S_DRAINAGE,                                         &
        SRC_S_WORK_STEPS, SRC_S_WORK_REJ, SRC_S_WORK_WATER_NSUB, SRC_S_WORK_WATER_NEWTON,           &
        SRC_S_WORK_WATER_NFAIL, SRC_S_WORK_WATER_REJECT, SRC_S_WORK_RK45_RESCUE,                     &
        SRC_S_WORK_CLAMP_STAGE,                                                                  &
        SRC_F_SOIL_TEMP, SRC_F_SOIL_WATER, SRC_F_PY0, SRC_F_PA0,                                &
        SRC_F_COH_LEAF_TEMP, SRC_F_COH_GPP, SRC_F_COH_HEIGHT, SRC_F_CO0, FLD_C_DIAG0, FLD_P_DIAG0, &
        FLD_PY_DIAG0, FLD_SI_DIAG0, FLD_SI_GROWTH_TEMP
   use meds_diag_vars,       only : CO_ANET, CO_AGROSS, CO_GSW, CO_GBW, CO_CI, CO_CS, CO_RD, CO_TRANSP, CO_BETA_STOM, &
                                    CO_BETA_NONSTOM, CO_LEAF_VPD, CO_PSI_LEAF, CO_PSI_WOOD, CO_PLC, CO_SAPFLOW, CO_ROOT_UPTAKE, &
                                    CO_ABS_PAR, CO_ABS_SW, CO_ABS_LW, CO_WIND, CO_LEAF_WATER, CO_WOOD_WATER, CO_LEAF_TEMP, &
                                    CO_WOOD_TEMP, CO_LEAF_PAR, CO_GB_MOL, CO_PSI_PREDAWN, CO_CA, CO_PRESSURE, PA_LE, PA_H, &
                                    PA_RNET, PA_SW_GROUND, PA_LW_GROUND, PA_SW_IN_VIS, PA_SW_IN_NIR, PA_SW_UP_VIS, PA_SW_UP_NIR, &
                                    PA_LW_UP, PA_USTAR, PA_GGNET, PA_ROUGH, PA_DISPLACE, PA_GPP, PA_NEE, PA_TRANSP, &
                                    PA_SOIL_TEMP_TOP, PA_RESID_ENERGY, PA_CAS_VPD, PA_W_SURFACE, PA_CAS_TEMP, PA_CAS_SHV, &
                                    PA_CAS_CO2, PA_WIND_CAS_TOP, PA_TAIR_CAS_TOP, PA_Z_CAS_TOP, PA_RESID_WATER, PA_NPP, PA_RECO, &
                                    PA_SW_UP, PA_SKIN_TEMP, PY_SW_IN, PY_PRECIP, PY_TAIR, PY_QAIR, PY_PSURF, PY_WIND, PY_LWDOWN, &
                                    PY_PAR_BEAM, PY_PAR_DIFFUSE, PY_NIR_BEAM, PY_NIR_DIFFUSE, PY_SNOWFALL, PY_CO2, PY_COSZ, &
                                    PY_RHO_AIR, PY_RAINF, N_PY_ROWS, CO_DDBH_DT, CO_DAGB_DT, CO_MORT_RATE, CO_NPP_LEAF, &
                                    CO_NPP_FINEROOT, CO_NPP_WOOD, CO_NPP_STORAGE, CO_NPP_REPRO, CO_GROWTH_RESP, CO_STORAGE_RESP, &
                                    CO_EXUDATE, PA_LITTER_LEAF, PA_LITTER_FINEROOT, PA_MORT_C_BACKGROUND, PA_MORT_C_CULL, &
                                    PA_MORT_C_DISTURB, PA_LITTER_STRUCT, PA_LITTER_WOOD, PA_LITTER_REPRO, PA_RECRUIT_NPLANT, &
                                    PA_DISTURB_AREA
   implicit none
   private

   public :: build_output_registry, build_freq_index, find_var_index, parse_stream_mask
   public :: apply_group_toggles, apply_freq_enables, apply_variable_override, apply_axis_toggles
   public :: freq_bit, dim_axis_index, OVR_TRUE, OVR_FALSE, OVR_MASK
   public :: manager_alloc, manager_setup, manager_finalize, manager_alloc_buffers
   public :: manager_set_polygons
   public :: activate_site_diag, dump_io_config

   !----- Named default stream masks (readable `ior` combinations). DAY_MON_YR deliberately     !
   !      EXCLUDES the fast bit -- the fast tier is always opt-in.                                !
   integer(ik), parameter :: DAY_MON    = ior(FREQ_DAILY, FREQ_MONTHLY)
   integer(ik), parameter :: MON        = FREQ_MONTHLY
   integer(ik), parameter :: MON_YR     = ior(FREQ_MONTHLY, FREQ_ANNUAL)
   integer(ik), parameter :: DAY_MON_YR = ior(ior(FREQ_DAILY, FREQ_MONTHLY), FREQ_ANNUAL)
   integer(ik), parameter :: FAST_ONLY  = FREQ_FAST

   !----- Per-variable override kinds (meds_io_config.toml value grammar,                    !
   !      docs/science/diagnostics.md §5). -------------------------------------------------!
   integer(ik), parameter :: OVR_TRUE  = 1_ik  !< bool true  -> restore registry-default streams
   integer(ik), parameter :: OVR_FALSE = 2_ik  !< bool false -> disable everywhere
   integer(ik), parameter :: OVR_MASK  = 3_ik  !< stream-string -> replace the stream mask

contains

   !----- FREQ_* bit for a 1-based tier index (1=FAST,2=DAILY,3=MONTHLY,4=ANNUAL). ------------!
   pure integer(ik) function freq_bit(tier) result(b)
      integer(ik), intent(in) :: tier
      b = ishft(1_ik, tier - 1_ik)
   end function freq_bit

   !=======================================================================================!
   !  Build the registry: the source list, then apply the main-config group/tier resolution. !
   !=======================================================================================!
   subroutine build_output_registry(reg, cfg)
      type(output_registry_t), intent(out) :: reg
      type(meds_config_t),     intent(in)  :: cfg

      allocate(reg%var(MAX_OUTPUT_VARS)) ; reg%nvar = 0_ik
      allocate(reg%idx_freq(MAX_OUTPUT_VARS, N_FREQ)) ; reg%idx_freq = 0_ik ; reg%nidx = 0_ik

      call register_structure_cohort(reg)
      call register_structure_patch(reg)
      call register_structure_site(reg)
      call register_structure_pft(reg)
      call register_structure_size(reg)
      call register_carbon(reg)
      call register_water(reg)
      call register_energy(reg)
      call register_biogeochem(reg)
      call register_numerics(reg)
      call register_ecophys(reg)
      call register_allocation(reg)
      call register_patch_fluxes(reg)
      call register_fast(reg)
      call register_forcing(reg)

      call enforce_annual_guard(reg)              ! cohort/patch var MUST NOT declare FREQ_ANNUAL (docs/science/diagnostics.md §3)
      call enforce_fast_guard(reg)                ! FAST tier <-> FAST-staged source, both ways
      call apply_axis_toggles(reg, cfg%output)    ! whole-axis suppression (docs/science/diagnostics.md §5, step 2)
      call apply_group_toggles(reg, cfg%output)   ! high-level group toggles (docs/science/diagnostics.md §5, step 3)
      call apply_freq_enables(reg, cfg%output)    ! per-tier enable (docs/science/diagnostics.md §5, step 4)
      call build_freq_index(reg)                  ! precompute idx_freq (per-variable overrides applied later)
   end subroutine build_output_registry

   !=======================================================================================!
   !  THE REGISTRATION LIST, grouped one subroutine per (group, axis) so a reviewer can diff   !
   !  one domain at a time. Adding a variable is still TWO coupled edits: an add_variable()      !
   !  line here + the matching FIELD case in meds_output_integrate. But because the reduction     !
   !  is DATA (dim + weight + mean + scale), a field that already has a case emits its patch /     !
   !  site / PFT / size-class twins for one more line EACH, with no new extraction code.           !
   !                                                                                          !
   !  add_variable(reg, name, long_name, units, dim, agg, group, streams, source                  !
   !               [, xtype][, weight][, mean][, scale])                                           !
   !=======================================================================================!

   !----- Per-cohort structure: the demographic atom's own state, written raw. --------------!
   subroutine register_structure_cohort(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'nplant_cohort', 'plant number density', 'plant/m2',                &
                        DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_C_NPLANT)
      call add_variable(reg, 'dbh_cohort', 'diameter at breast height', 'cm',                    &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_DBH)
      call add_variable(reg, 'height_cohort', 'height', 'm',                                     &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_HEIGHT)
      call add_variable(reg, 'basal_area_cohort', 'basal area per plant', 'cm2/plant',           &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_BASAL_AREA)
      call add_variable(reg, 'agb_cohort', 'aboveground biomass per plant', 'kgC/plant',         &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_AGB)
      call add_variable(reg, 'leaf_area_cohort', 'leaf area per plant', 'm2/plant',              &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_LEAF_AREA)
      call add_variable(reg, 'lai_cohort', 'cohort leaf area index', 'm2/m2',                    &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_LAI)
      call add_variable(reg, 'wai_cohort', 'cohort wood (stem and branch) area index', 'm2/m2',  &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_WAI)
      call add_variable(reg, 'growth_avg_cohort', 'moving-average growth', 'cm/yr',              &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_GROWTH_AVG)
      call add_variable(reg, 'overtopping_lai_cohort', 'cumulative LAI of taller cohorts', 'm2/m2', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_OVERTOP_LAI)
      call add_variable(reg, 'pft_cohort', 'plant functional type index', '-',                   &
                        DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_C_PFT, xt=XTYPE_INT)
      call add_variable(reg, 'owner_patch', 'owning patch (1-based)', '-',                       &
                        DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_C_OWNER_PATCH, xt=XTYPE_INT)
      call add_variable(reg, 'global_cohort_id', 'persistent cohort id', '-',                    &
                        DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_C_GLOBAL_ID, xt=XTYPE_INT)
      !----- Per-plant carbon pools (the four PARTEH-style pools the size anchor rides on). ---!
      call add_variable(reg, 'leaf_carbon_cohort', 'leaf carbon per plant', 'kgC/plant',         &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_LEAF_CARBON)
      call add_variable(reg, 'fineroot_carbon_cohort', 'fine-root carbon per plant', 'kgC/plant', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_FINEROOT_CARBON)
      call add_variable(reg, 'wood_carbon_cohort', 'wood carbon per plant', 'kgC/plant',         &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_WOOD_CARBON)
      call add_variable(reg, 'storage_carbon_cohort', 'non-structural carbon per plant', 'kgC/plant', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_STORAGE_CARBON)
      !----- Dynamic (light-acclimating) leaf traits. -----------------------------------------!
      call add_variable(reg, 'sla_cohort', 'specific leaf area', 'm2/kgC',                       &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_SLA)
      call add_variable(reg, 'vcmax25_cohort', 'max carboxylation rate at 25 degC', 'umol/m2/s', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_VCMAX25)
      call add_variable(reg, 'rd25_cohort', 'leaf dark respiration at 25 degC', 'umol/m2/s',     &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_RD25)
      call add_variable(reg, 'leaf_turnover_rate_cohort', 'baseline leaf turnover rate', '1/yr', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_LEAF_TURNOVER)
      call add_variable(reg, 'leaf_flush_tendency_cohort', 'smoothed leaf flush tendency', '-',          &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_PHENO_FLUSH)
      call add_variable(reg, 'leaf_shed_tendency_cohort', 'smoothed leaf senescence tendency', '-',      &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_PHENO_SHED)
      call add_variable(reg, 'dmax_psi_leaf_cohort', 'predawn (daily-max) leaf water potential', 'MPa', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_DMAX_PSI_LEAF)
      !----- Slow-loop RATES off the per-step diagnostic block (poly%cohort%sdiag). -----------!
      call add_variable(reg, 'dbh_growth_cohort', 'diameter growth rate', 'cm/yr',               &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_DIAG0 + CO_DDBH_DT)
      call add_variable(reg, 'agb_growth_cohort', 'per-plant AGB growth rate', 'kgC/plant/yr',   &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_DIAG0 + CO_DAGB_DT)
      call add_variable(reg, 'mort_rate_cohort', 'mortality hazard rate', '1/yr',                &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_DIAG0 + CO_MORT_RATE)
   end subroutine register_structure_cohort

   !----- Per-patch structure: geometry + identity, plus cohort fields reduced to the patch. -!
   subroutine register_structure_patch(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'area_patch', 'patch area fraction', '-',                           &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_AREA)
      call add_variable(reg, 'age_patch', 'time since disturbance', 'yr',                        &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, MON, FLD_P_AGE)
      call add_variable(reg, 'dist_type_patch', 'disturbance type', '-',                         &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, MON, FLD_P_DIST_TYPE, xt=XTYPE_INT)
      call add_variable(reg, 'cohort_offset', 'first cohort of patch (CSR)', '-',                &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_COHORT_OFFSET, xt=XTYPE_INT)
      call add_variable(reg, 'cohort_count', 'cohorts in patch (CSR)', '-',                      &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_COHORT_COUNT, xt=XTYPE_INT)
      call add_variable(reg, 'global_patch_id', 'persistent patch id', '-',                      &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_GLOBAL_ID, xt=XTYPE_INT)
      call add_variable(reg, 'patch_site', 'the patch''s site (1-based)', '-',                   &
                        DIM_PATCH, AGG_LAST, GRP_STRUCTURE, DAY_MON, FLD_P_SITE, xt=XTYPE_INT)
      !----- Cohort fields reduced onto the patch axis. EXTENSIVE (weighted SUM over nplant),   !
      !      so these are per m2 of THAT PATCH's ground -- the gap-vs-closed-canopy contrast the  !
      !      site scalars average away.  -------------------------------------------------------!
      call add_variable(reg, 'nplant_patch', 'patch plant number density', 'plant/m2',           &
                        DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_NPLANT, w=W_NONE)
      call add_variable(reg, 'lai_patch', 'patch leaf area index', 'm2/m2',                      &
                        DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_LEAF_AREA, w=W_NPLANT)
      call add_variable(reg, 'agb_patch', 'patch aboveground biomass', 'kgC/m2',                 &
                        DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_AGB, w=W_NPLANT)
      call add_variable(reg, 'basal_area_patch', 'patch basal area', 'm2/m2',                    &
                        DIM_PATCH, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_BASAL_AREA,               &
                        w=W_NPLANT, sc=cm2_to_m2)
   end subroutine register_structure_patch

   !----- Site structure scalars. -----------------------------------------------------------!
   subroutine register_structure_site(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'nplant_site', 'site total plant number', 'plant/m2',               &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_NPLANT, w=W_NONE)
      call add_variable(reg, 'basal_area_site', 'site total basal area', 'm2/m2',                &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_BASAL_AREA,       &
                        w=W_NPLANT, sc=cm2_to_m2)
      call add_variable(reg, 'agb_site', 'site aboveground biomass', 'kgC/m2',                   &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_AGB, w=W_NPLANT)
      call add_variable(reg, 'lai_site', 'site leaf area index', 'm2/m2',                        &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_LEAF_AREA, w=W_NPLANT)
      !----- INTENSIVE structural means. The weight is the physical statement: dbh is basal-    !
      !      area-weighted (the forestry convention, dominated by the trees that hold the stand), !
      !      NOT stem-weighted, which a regenerating understory would swamp.  -------------------!
      call add_variable(reg, 'mean_dbh_site', 'basal-area-weighted mean DBH', 'cm',              &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_DBH,             &
                        w=W_BASAL_AREA, mn=.true.)
      call add_variable(reg, 'mean_height_site', 'basal-area-weighted mean height', 'm',         &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, FLD_C_HEIGHT,          &
                        w=W_BASAL_AREA, mn=.true.)
      call add_variable(reg, 'canopy_height_site', 'height of the tallest cohort', 'm',          &
                        DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR, SRC_S_CANOPY_HEIGHT)
      call add_variable(reg, 'n_cohort_site', 'live cohort count', '-',                          &
                        DIM_SCALAR, AGG_LAST, GRP_STRUCTURE, DAY_MON_YR, SRC_S_N_COHORT, xt=XTYPE_INT)
      call add_variable(reg, 'n_patch_site', 'live patch count', '-',                            &
                        DIM_SCALAR, AGG_LAST, GRP_STRUCTURE, DAY_MON_YR, SRC_S_N_PATCH, xt=XTYPE_INT)
   end subroutine register_structure_site

   !----- PFT axis. The composition view: is any PFT actually competing, or has one won? ----!
   subroutine register_structure_pft(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'nplant_pft', 'plant number density by PFT', 'plant/m2',            &
                        DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_NPLANT, w=W_NONE)
      call add_variable(reg, 'agb_pft', 'aboveground biomass by PFT', 'kgC/m2',                  &
                        DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_AGB, w=W_NPLANT)
      call add_variable(reg, 'lai_pft', 'leaf area index by PFT', 'm2/m2',                       &
                        DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_LEAF_AREA, w=W_NPLANT)
      call add_variable(reg, 'basal_area_pft', 'basal area by PFT', 'm2/m2',                     &
                        DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_BASAL_AREA,             &
                        w=W_NPLANT, sc=cm2_to_m2)
      call add_variable(reg, 'mean_dbh_pft', 'basal-area-weighted mean DBH by PFT', 'cm',        &
                        DIM_PFT, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_DBH,                    &
                        w=W_BASAL_AREA, mn=.true.)
   end subroutine register_structure_pft

   !----- DBH size-class axis (ED2-style, a size class of PLANTS). The demographic core's      !
   !      actual output: a stem-density distribution directly comparable to forest inventory.   !
   subroutine register_structure_size(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'nplant_size', 'plant number density by DBH class', 'plant/m2',     &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_NPLANT, w=W_NONE)
      call add_variable(reg, 'agb_size', 'aboveground biomass by DBH class', 'kgC/m2',           &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_AGB, w=W_NPLANT)
      call add_variable(reg, 'lai_size', 'leaf area index by DBH class', 'm2/m2',                &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_LEAF_AREA, w=W_NPLANT)
      call add_variable(reg, 'basal_area_size', 'basal area by DBH class', 'm2/m2',              &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_BASAL_AREA,            &
                        w=W_NPLANT, sc=cm2_to_m2)
      call add_variable(reg, 'agb_growth_size', 'AGB growth rate by DBH class', 'kgC/m2/yr',     &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'mort_rate_size', 'nplant-weighted mortality rate by DBH class', '1/yr', &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE,             &
                        w=W_NPLANT, mn=.true.)
   end subroutine register_structure_size

   !=======================================================================================!
   !  CARBON. gpp_accum / *_resp_accum are per-slow-step accumulators the fast loop fills, so  !
   !  AGG_SUM over the period gives a period carbon TOTAL (all 0 when the fast loop is off).    !
   !=======================================================================================!
   subroutine register_carbon(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'gpp_site', 'site gross primary productivity (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_CARBON, DAY_MON_YR, FLD_C_GPP_ACCUM, w=W_NPLANT)
      call add_variable(reg, 'npp_site', 'site net primary productivity (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_CARBON, DAY_MON_YR, FLD_C_NPP_ACCUM, w=W_NPLANT)
      call add_variable(reg, 'leaf_resp_site', 'site leaf dark respiration (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_CARBON, DAY_MON_YR, FLD_C_LEAF_RESP, w=W_NPLANT)
      call add_variable(reg, 'stem_resp_site', 'site stem maintenance respiration (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_CARBON, DAY_MON_YR, FLD_C_STEM_RESP, w=W_NPLANT)
      call add_variable(reg, 'root_resp_site', 'site fine-root maintenance respiration (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_CARBON, DAY_MON_YR, FLD_C_ROOT_RESP, w=W_NPLANT)
      !----- Vegetation carbon STOCKS (AGG_MEAN -- a stock, not a flux). ----------------------!
      call add_variable(reg, 'leaf_carbon_site', 'site leaf carbon', 'kgC/m2',                   &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_LEAF_CARBON, w=W_NPLANT)
      call add_variable(reg, 'fineroot_carbon_site', 'site fine-root carbon', 'kgC/m2',          &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_FINEROOT_CARBON, w=W_NPLANT)
      call add_variable(reg, 'wood_carbon_site', 'site wood carbon', 'kgC/m2',                   &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_WOOD_CARBON, w=W_NPLANT)
      call add_variable(reg, 'storage_carbon_site', 'site non-structural carbon', 'kgC/m2',      &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_STORAGE_CARBON, w=W_NPLANT)
      call add_variable(reg, 'bgb_site', 'site belowground biomass', 'kgC/m2',                   &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_BGB, w=W_NPLANT)
      call add_variable(reg, 'veg_carbon_site', 'site total live vegetation carbon', 'kgC/m2',   &
                        DIM_SCALAR, AGG_MEAN, GRP_CARBON, DAY_MON_YR, FLD_C_VEG_CARBON, w=W_NPLANT)
      !----- PFT-resolved carbon fluxes. ------------------------------------------------------!
      call add_variable(reg, 'gpp_pft', 'gross primary productivity by PFT (period total)', 'kgC/m2', &
                        DIM_PFT, AGG_SUM, GRP_CARBON, MON_YR, FLD_C_GPP_ACCUM, w=W_NPLANT)
      call add_variable(reg, 'npp_pft', 'net primary productivity by PFT (period total)', 'kgC/m2', &
                        DIM_PFT, AGG_SUM, GRP_CARBON, MON_YR, FLD_C_NPP_ACCUM, w=W_NPLANT)
      !----- Per-patch carbon fluxes. ---------------------------------------------------------!
      call add_variable(reg, 'gpp_patch', 'patch gross primary productivity (period total)', 'kgC/m2', &
                        DIM_PATCH, AGG_SUM, GRP_CARBON, DAY_MON, FLD_C_GPP_ACCUM, w=W_NPLANT)
      call add_variable(reg, 'npp_patch', 'patch net primary productivity (period total)', 'kgC/m2', &
                        DIM_PATCH, AGG_SUM, GRP_CARBON, DAY_MON, FLD_C_NPP_ACCUM, w=W_NPLANT)
      !----- Per-cohort carbon fluxes (per plant). --------------------------------------------!
      call add_variable(reg, 'gpp_cohort', 'per-plant GPP (period total)', 'kgC/plant',          &
                        DIM_COHORT, AGG_SUM, GRP_CARBON, DAY_MON, FLD_C_GPP_ACCUM)
      call add_variable(reg, 'npp_cohort', 'per-plant NPP (period total)', 'kgC/plant',          &
                        DIM_COHORT, AGG_SUM, GRP_CARBON, DAY_MON, FLD_C_NPP_ACCUM)
   end subroutine register_carbon

   !=======================================================================================!
   !  WATER.                                                                                 !
   !=======================================================================================!
   subroutine register_water(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'et_site', 'site evapotranspiration (period total)', 'kg/m2',       &
                        DIM_SCALAR, AGG_SUM, GRP_WATER, DAY_MON_YR, SRC_S_ET)
      !----- The water leaving the soil column: surface runoff (the pond's spill above its capacity) !
      !      and drainage (the bottom flux: free drainage, 0 on bedrock, net and signed with an       !
      !      aquifer), and their sum. ------------------------------------------------------------------!
      call add_variable(reg, 'runoff_site', 'site runoff: surface runoff + drainage (period total)', 'kg/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_WATER, DAY_MON_YR, SRC_S_RUNOFF)
      call add_variable(reg, 'runoff_surf_site', 'site surface runoff (period total)', 'kg/m2',  &
                        DIM_SCALAR, AGG_SUM, GRP_WATER, DAY_MON_YR, SRC_S_RUNOFF_SURF)
      call add_variable(reg, 'drainage_site', 'site drainage at the soil column''s bottom (period total)', 'kg/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_WATER, DAY_MON_YR, SRC_S_DRAINAGE)
      !----- Each site's soil column, its patches area-weighted (DIM_SOIL). --------------------!
      call add_variable(reg, 'soil_water_site', 'each site''s area-weighted volumetric soil moisture', 'm3/m3', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_WATER)
      call add_variable(reg, 'soil_psi_site', 'each site''s area-weighted soil matric potential', 'MPa', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_PSI)
      call add_variable(reg, 'soil_wetness_site', 'each site''s area-weighted relative saturation', '-', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_WETNESS)
      !----- Surface water stores. ------------------------------------------------------------!
      call add_variable(reg, 'w_surface_site', 'ponded surface water', 'kg/m2',                  &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_W_SURFACE)
      call add_variable(reg, 'swe_site', 'snow water equivalent', 'kg/m2',                       &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_SWE)
      call add_variable(reg, 'snow_depth_site', 'snow depth', 'm',                               &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_SNOW_DEPTH)
      call add_variable(reg, 'w_surface_patch', 'patch ponded surface water', 'kg/m2',           &
                        DIM_PATCH, AGG_TMEAN, GRP_WATER, MON, FLD_P_DIAG0 + PA_W_SURFACE)
      call add_variable(reg, 'swe_patch', 'patch snow water equivalent', 'kg/m2',                &
                        DIM_PATCH, AGG_TMEAN, GRP_WATER, MON, FLD_P_SWE)
      !----- 2-D (soil layer x patch) columns. Gated by axes_soil_patch (default off): the      !
      !      gap-vs-closed-canopy drydown contrast the area-weighted site column hides. --------!
      call add_variable(reg, 'soil_water_layer_patch', 'volumetric soil moisture by layer and patch', &
                        'm3/m3', DIM_SOIL_PATCH, AGG_TMEAN, GRP_WATER, DAY_MON, FLD_L_SOIL_WATER)
      call add_variable(reg, 'soil_psi_layer_patch', 'soil matric potential by layer and patch', &
                        'MPa', DIM_SOIL_PATCH, AGG_TMEAN, GRP_WATER, DAY_MON, FLD_L_SOIL_PSI)
      call add_variable(reg, 'soil_wetness_layer_patch', 'relative saturation by layer and patch', &
                        '-', DIM_SOIL_PATCH, AGG_TMEAN, GRP_WATER, MON, FLD_L_SOIL_WETNESS)
   end subroutine register_water

   !=======================================================================================!
   !  ENERGY.                                                                                !
   !=======================================================================================!
   subroutine register_energy(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'cas_temp_site', 'site canopy-air-space temperature', 'K',          &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_TEMP)
      call add_variable(reg, 'cas_shv_site', 'site canopy-air specific humidity', 'kg/kg',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_SHV)
      call add_variable(reg, 'cas_co2_site', 'site canopy-air CO2 mixing ratio', 'umol/mol',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_CO2)
      call add_variable(reg, 'cas_vpd_site', 'site canopy-air vapour-pressure deficit', 'Pa',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_CAS_VPD)
      call add_variable(reg, 'cas_depth_site', 'site canopy-air-space depth', 'm',               &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_CAS_DEPTH)
      call add_variable(reg, 'soil_temp_top_site', 'site soil-top (layer 1) temperature', 'K',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SOIL_TEMP_TOP)
      !----- The skin temperature (docs/science/diagnostics.md §4): what an infrared thermometer   !
      !      or a satellite land-surface temperature sees, not the soil. Over a forest it is mostly  !
      !      the canopy. ------------------------------------------------------------------------------!
      call add_variable(reg, 'skin_temp_site', 'skin temperature: black-body temperature of the longwave up', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SKIN_TEMP)
      call add_variable(reg, 'soil_temp_site', 'each site''s area-weighted soil temperature', 'K', &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_L_SOIL_TEMP)
      call add_variable(reg, 'soil_fliq_site', 'each site''s area-weighted soil liquid fraction', '-', &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_L_SOIL_FLIQ)

      !----- Canopy temperatures, LEAF-AREA-weighted (the intensive rule: a bare sapling must    !
      !      not pull the canopy mean as hard as a closed overstory).  --------------------------!
      call add_variable(reg, 'leaf_temp_site', 'leaf-area-weighted canopy leaf temperature', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_C_DIAG0 + CO_LEAF_TEMP, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'wood_temp_site', 'leaf-area-weighted wood temperature', 'K',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_C_DIAG0 + CO_WOOD_TEMP, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'leaf_temp_cohort', 'per-cohort leaf temperature', 'K',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_C_DIAG0 + CO_LEAF_TEMP)
      call add_variable(reg, 'wood_temp_cohort', 'per-cohort wood temperature', 'K',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_C_DIAG0 + CO_WOOD_TEMP)
      call add_variable(reg, 'cas_temp_patch', 'patch canopy-air-space temperature', 'K',        &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_CAS_TEMP)
      call add_variable(reg, 'soil_temp_top_patch', 'patch soil-top temperature', 'K',           &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, MON, FLD_P_DIAG0 + PA_SOIL_TEMP_TOP)
      call add_variable(reg, 'skin_temp_patch', 'patch skin temperature (longwave up)', 'K',     &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_SKIN_TEMP)
      call add_variable(reg, 'soil_temp_layer_patch', 'soil temperature by layer and patch', 'K', &
                        DIM_SOIL_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_L_SOIL_TEMP)
      call add_variable(reg, 'soil_fliq_layer_patch', 'soil liquid fraction by layer and patch', '-', &
                        DIM_SOIL_PATCH, AGG_TMEAN, GRP_ENERGY, MON, FLD_L_SOIL_FLIQ)
   end subroutine register_energy

   !=======================================================================================!
   !  BIOGEOCHEMISTRY (the CENTURY matrix pools + Rh). All 0 when soil_carbon_on is off.      !
   !=======================================================================================!
   subroutine register_biogeochem(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'soilc_fast_grnd_site', 'site fast/metabolic litter carbon, above-ground', &
                        'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_FAST_GRND)
      call add_variable(reg, 'soilc_fast_soil_site', 'site fast/metabolic litter carbon, below-ground', &
                        'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_FAST_SOIL)
      call add_variable(reg, 'soilc_struct_grnd_site', 'site structural litter + CWD carbon, above-ground', &
                        'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_STRUCT_GRND)
      call add_variable(reg, 'soilc_struct_soil_site', 'site structural litter + CWD carbon, below-ground', &
                        'kgC/m2', DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_STRUCT_SOIL)
      call add_variable(reg, 'soilc_microbial_site', 'site microbial SOM carbon', 'kgC/m2',      &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_MICROBIAL)
      call add_variable(reg, 'soilc_slow_site', 'site slow/humified SOM carbon', 'kgC/m2',       &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_SLOW)
      call add_variable(reg, 'soilc_passive_site', 'site passive SOM carbon', 'kgC/m2',          &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_PASSIVE)
      call add_variable(reg, 'soilc_total_site', 'site total soil organic carbon', 'kgC/m2',     &
                        DIM_SCALAR, AGG_MEAN, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_SOILC_TOTAL)
      call add_variable(reg, 'rh_site', 'site heterotrophic respiration (period total)', 'kgC/m2', &
                        DIM_SCALAR, AGG_SUM, GRP_BIOGEOCHEM, DAY_MON_YR, FLD_P_RH)
      call add_variable(reg, 'soilc_total_patch', 'patch total soil organic carbon', 'kgC/m2',   &
                        DIM_PATCH, AGG_MEAN, GRP_BIOGEOCHEM, MON, FLD_P_SOILC_TOTAL)
      call add_variable(reg, 'rh_patch', 'patch heterotrophic respiration (period total)', 'kgC/m2', &
                        DIM_PATCH, AGG_SUM, GRP_BIOGEOCHEM, MON, FLD_P_RH)
   end subroutine register_biogeochem

   !=======================================================================================!
   !  NUMERICS: integrator WORK + health. AGG_SUM throughout -- these are counts of work done  !
   !  over the period, so they add; a period MEAN would hide the single day that went wrong,   !
   !  which is exactly the signal being sought.                                                !
   !=======================================================================================!
   subroutine register_numerics(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'work_integ_steps_site', 'accepted integrator sub-steps (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_STEPS)
      call add_variable(reg, 'work_integ_rej_site', 'rejected integrator steps (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_REJ)
      call add_variable(reg, 'work_water_nsub_site', 'soil-plant water network sub-steps (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_WATER_NSUB)
      call add_variable(reg, 'work_water_newton_site', 'soil-plant water network Newton iterations (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_WATER_NEWTON)
      call add_variable(reg, 'work_water_reject_site',                                          &
                        'soil-plant water network sub-steps rejected and retried shorter (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_WATER_REJECT)
      call add_variable(reg, 'work_water_nfail_site',                                           &
                        'soil-plant water network sub-steps taken although unconverged (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_WATER_NFAIL)
      call add_variable(reg, 'work_rk45_rescue_site', 'RK45 steps rescued to the implicit path (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_RK45_RESCUE)
      call add_variable(reg, 'work_clamp_stage_site', 'stage-input clamp activations (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_CLAMP_STAGE)
   end subroutine register_numerics


   !=======================================================================================!
   !  GRP_ECOPHYS -- the per-cohort leaf gas-exchange + hydraulics set. Every row reads the       !
   !  fast-loop diagnostic block, so every row is a quantity the model computes ~48 times a day    !
   !  and keeps nowhere else.                                                                       !
   !                                                                                          !
   !  OFF by default: this is the highest-volume group (cohort axis x ~20 variables), and it is     !
   !  the one a production run most often wants off. `ecophys = true` is the single switch.         !
   !                                                                                          !
   !  The canopy-level twins (leaf-area-weighted means over the cohorts) are registered alongside    !
   !  the raw slabs, because a leaf-area-weighted canopy gsw and a canopy-mean psi_leaf are the      !
   !  standard flux-tower comparison quantities and a reader should not have to re-derive them.      !
   !=======================================================================================!
   subroutine register_ecophys(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'anet_cohort', 'net leaf assimilation', 'umol/m2 leaf/s',           &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_ANET)
      call add_variable(reg, 'agross_cohort', 'gross leaf assimilation', 'umol/m2 leaf/s',       &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_AGROSS)
      call add_variable(reg, 'gsw_cohort', 'stomatal conductance to water vapour', 'mol/m2 leaf/s', &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_GSW)
      call add_variable(reg, 'gbw_cohort', 'leaf boundary-layer conductance', 'm/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_GBW)
      call add_variable(reg, 'ci_cohort', 'intercellular CO2 mole fraction', 'umol/mol',         &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_CI)
      call add_variable(reg, 'cs_cohort', 'leaf-surface CO2 mole fraction', 'umol/mol',          &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_CS)
      call add_variable(reg, 'rd_cohort', 'leaf dark respiration', 'umol/m2 leaf/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_RD)
      call add_variable(reg, 'transp_cohort', 'leaf transpiration (the water the leaves lose)', 'mol/m2 leaf/s', &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_TRANSP)
      call add_variable(reg, 'beta_stomata_cohort', 'stomatal water-stress factor', '-',         &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_BETA_STOM)
      call add_variable(reg, 'beta_nonstomata_cohort', 'non-stomatal (capacity) water-stress factor', '-', &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_BETA_NONSTOM)
      call add_variable(reg, 'leaf_vpd_cohort', 'leaf-to-air vapour-pressure deficit', 'Pa',     &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_LEAF_VPD)
      call add_variable(reg, 'psi_leaf_cohort', 'leaf water potential', 'MPa',                   &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_PSI_LEAF)
      call add_variable(reg, 'psi_wood_cohort', 'wood water potential', 'MPa',                   &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_PSI_WOOD)
      call add_variable(reg, 'plc_cohort', 'plant loss of hydraulic conductance', '-',           &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_PLC)
      call add_variable(reg, 'sapflow_cohort', 'wood-to-leaf sapflow', 'kg/plant/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_SAPFLOW)
      call add_variable(reg, 'root_uptake_cohort', 'soil-to-root water uptake', 'kg/plant/s',    &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CO_ROOT_UPTAKE)
      call add_variable(reg, 'leaf_water_cohort', 'internal leaf water', 'kg/plant',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_LEAF_WATER)
      call add_variable(reg, 'wood_water_cohort', 'internal wood water', 'kg/plant',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_WOOD_WATER)
      call add_variable(reg, 'abs_par_cohort', 'absorbed photosynthetically active radiation', 'W/m2', &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_C_DIAG0 + CO_ABS_PAR)
      call add_variable(reg, 'abs_sw_cohort', 'absorbed shortwave radiation', 'W/m2',            &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_C_DIAG0 + CO_ABS_SW)
      call add_variable(reg, 'abs_lw_cohort', 'net longwave radiation', 'W/m2',                  &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, MON, FLD_C_DIAG0 + CO_ABS_LW)
      call add_variable(reg, 'wind_cohort', 'in-canopy wind speed', 'm/s',                       &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CO_WIND)
      !----- CANOPY-LEVEL twins: leaf-area-weighted means, the standard flux-tower comparands.  --!
      call add_variable(reg, 'anet_site', 'leaf-area-weighted canopy net assimilation', 'umol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CO_ANET,   &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'gsw_site', 'leaf-area-weighted canopy stomatal conductance', 'mol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CO_GSW,    &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'psi_leaf_site', 'leaf-area-weighted canopy leaf water potential', 'MPa', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CO_PSI_LEAF, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'plc_site', 'leaf-area-weighted canopy loss of conductance', '-',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CO_PLC,    &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'beta_stomata_site', 'leaf-area-weighted stomatal water-stress factor', '-', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CO_BETA_STOM, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'root_uptake_site', 'site root water uptake', 'kg/m2/s',            &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CO_ROOT_UPTAKE, &
                        w=W_NPLANT)
      call add_variable(reg, 'transp_leaf_site', 'leaf-area-weighted canopy transpiration', 'mol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CO_TRANSP,   &
                        w=W_LEAF_AREA, mn=.true.)
   end subroutine register_ecophys

   !=======================================================================================!
   !  SLOW-loop carbon allocation: where compute_carbon_allocation sends each cohort's carbon   !
   !  on every step, captured so the carbon budget closes from an output file -- npp_site       !
   !  reconciles against the sum of its five destinations plus growth respiration.              !
   !=======================================================================================!
   subroutine register_allocation(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'npp_leaf_site', 'NPP allocated to leaf', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_LEAF, w=W_NPLANT)
      call add_variable(reg, 'npp_fineroot_site', 'NPP allocated to fine root', 'kgC/m2/yr',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_FINEROOT, w=W_NPLANT)
      call add_variable(reg, 'npp_wood_site', 'NPP allocated to wood', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_WOOD, w=W_NPLANT)
      call add_variable(reg, 'npp_storage_site', 'NPP to non-structural storage', 'kgC/m2/yr',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_STORAGE, w=W_NPLANT)
      call add_variable(reg, 'npp_repro_site', 'NPP allocated to reproduction', 'kgC/m2/yr',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_NPP_REPRO, w=W_NPLANT)
      call add_variable(reg, 'growth_resp_site', 'growth respiration', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_GROWTH_RESP, w=W_NPLANT)
      call add_variable(reg, 'storage_resp_site', 'storage maintenance respiration', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_STORAGE_RESP, w=W_NPLANT)
      call add_variable(reg, 'root_exudate_site', 'root exudate: carbon the growth sink could not use', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_DIAG0 + CO_EXUDATE, w=W_NPLANT)
      !----- Site-level demographic RATES: the three terms that MAKE the AGB trajectory, so a      !
      !      reader can check d(agb)/dt against growth - mortality directly from the file.  -------!
      call add_variable(reg, 'agb_growth_site', 'site AGB growth rate', 'kgC/m2/yr',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'agb_mort_site', 'site AGB mortality loss rate', 'kgC/m2/yr',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, w=W_AGB)
      call add_variable(reg, 'mort_rate_site', 'nplant-weighted mortality rate', '1/yr',         &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, &
                        w=W_NPLANT, mn=.true.)
      !----- Litterfall: the vegetation -> soil carbon seam, from the SAME accumulator the        !
      !      biogeochemistry consumes, by the tissue the carbon came from. Leaf + fine root          !
      !      + wood + fruit/seed + root_exudate_site is the carbon turnover and background mortality !
      !      send to the soil pools; cohorts culled or killed by disturbance add to the pools        !
      !      directly (mort_carbon_cull_site, mort_carbon_disturb_site). With soil carbon off they   !
      !      leave the model. litter_struct_site is a destination (the structural pools), so it     !
      !      overlaps the others. --------------------------------------------------------------------!
      call add_variable(reg, 'litter_leaf_site', 'leaf litterfall carbon (after resorption)', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_LEAF)
      call add_variable(reg, 'litter_fineroot_site', 'fine-root litter carbon (after resorption)', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_FINEROOT)
      call add_variable(reg, 'litter_wood_site', 'wood and storage litter carbon of dead plants', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_WOOD)
      call add_variable(reg, 'litter_repro_site', 'fruit and seed litter carbon', 'kgC/m2/yr',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_REPRO)
      call add_variable(reg, 'litter_struct_site', 'structural litter + CWD carbon', 'kgC/m2/yr', &
                        DIM_SCALAR, AGG_TMEAN, GRP_BIOGEOCHEM, MON_YR, FLD_P_DIAG0 + PA_LITTER_STRUCT)
      call add_variable(reg, 'nplant_recruit_site', 'recruitment density flux', 'plant/m2/yr',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_RECRUIT_NPLANT)
   call add_variable(reg, 'disturb_area_site', 'patch area fraction disturbed', '1/yr',           &
                     DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_P_DIAG0 + PA_DISTURB_AREA)
      !----- Mortality carbon SPLIT BY PATHWAY. The three are mutually exclusive and jointly       !
      !      exhaustive, so they sum to the whole-individual mortality carbon -- which is what makes !
      !      them worth emitting separately: a stand thinning continuously and a stand being knocked !
      !      over lose the same carbon and mean completely different things. Every pool (leaf, fine  !
      !      root, wood, storage) is included, because a whole individual dies with all of it.       !
      !      Emitted whatever [soil_carbon].soil_carbon_on says: this is demography.  ---------------!
      call add_variable(reg, 'mort_carbon_background_site', 'mortality carbon, background hazard',  &
                        'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR,                  &
                        FLD_P_DIAG0 + PA_MORT_C_BACKGROUND)
      call add_variable(reg, 'mort_carbon_cull_site', 'mortality carbon, cohorts culled below the tracking floor', &
                        'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR,                  &
                        FLD_P_DIAG0 + PA_MORT_C_CULL)
      call add_variable(reg, 'mort_carbon_disturb_site', 'mortality carbon, canopy killed by disturbance', &
                        'kgC/m2/yr', DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR,                  &
                        FLD_P_DIAG0 + PA_MORT_C_DISTURB)
      call add_variable(reg, 'mort_carbon_background_patch', 'patch mortality carbon, background hazard', &
                        'kgC/m2/yr', DIM_PATCH, AGG_TMEAN, GRP_STRUCTURE, MON,                      &
                        FLD_P_DIAG0 + PA_MORT_C_BACKGROUND)
      call add_variable(reg, 'litter_leaf_patch', 'patch leaf litterfall carbon', 'kgC/m2/yr',   &
                        DIM_PATCH, AGG_TMEAN, GRP_BIOGEOCHEM, MON, FLD_P_DIAG0 + PA_LITTER_LEAF)
      call add_variable(reg, 'agb_growth_pft', 'AGB growth rate by PFT', 'kgC/m2/yr',            &
                        DIM_PFT, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'agb_mort_pft', 'AGB mortality loss rate by PFT', 'kgC/m2/yr',      &
                        DIM_PFT, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_DIAG0 + CO_MORT_RATE, w=W_AGB)
   end subroutine register_allocation

   !=======================================================================================!
   !  The per-PATCH fast-loop fluxes and states, and their area-weighted site twins. These are  !
   !  the sub-daily energy / water / carbon exchange terms -- the ones a run is judged on --      !
   !  available at EVERY tier, not only in the opt-in FAST stream.                                !
   !=======================================================================================!
   subroutine register_patch_fluxes(reg)
      type(output_registry_t), intent(inout) :: reg
      !--- energy ---!
      call add_variable(reg, 'le_site', 'latent heat flux', 'W/m2',                              &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_LE)
      call add_variable(reg, 'h_site', 'sensible heat flux', 'W/m2',                             &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_H)
      call add_variable(reg, 'rnet_site', 'net all-wave radiation', 'W/m2',                      &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_RNET)
      !----- TOP-OF-CANOPY radiative fluxes. The two-stream solves these every step; with them a    !
      !      run can be compared against a radiometer or a satellite product without re-deriving    !
      !      them offline.                                                                          !
      !                                                                                          !
      !      FLUXES, not a time-averaged albedo, and deliberately so: a period-mean albedo is the    !
      !      mean of a RATIO, which is not the ratio of the means, and at night the shortwave ratio  !
      !      is 0/0. The albedo a reader wants -- and what a satellite product IS -- is                !
      !      sum(up)/sum(down) over the period, which these make computable exactly                   !
      !      (docs/science/canopy_radiation_transfer.md eq. (15)).                                     !
      call add_variable(reg, 'sw_in_vis_site', 'incident VIS shortwave at canopy top', 'W/m2',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_IN_VIS)
      call add_variable(reg, 'sw_in_nir_site', 'incident NIR shortwave at canopy top', 'W/m2',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_IN_NIR)
      call add_variable(reg, 'sw_up_vis_site', 'upwelling VIS shortwave at canopy top', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_UP_VIS)
      call add_variable(reg, 'sw_up_nir_site', 'upwelling NIR shortwave at canopy top', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_SW_UP_NIR)
      call add_variable(reg, 'lw_up_site', 'upwelling longwave at canopy top (emission included)', &
                        'W/m2', DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_LW_UP)
      call add_variable(reg, 'sw_in_site', 'incident shortwave at canopy top', 'W/m2',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_SW_IN)
      call add_variable(reg, 'sw_ground_site', 'shortwave absorbed at the ground', 'W/m2',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_P_DIAG0 + PA_SW_GROUND)
      call add_variable(reg, 'lw_ground_site', 'net longwave at the ground', 'W/m2',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_P_DIAG0 + PA_LW_GROUND)
      call add_variable(reg, 'ustar_site', 'friction velocity', 'm/s',                           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PA_USTAR)
      call add_variable(reg, 'ggnet_site', 'ground conductance', 'm/s',                          &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON, FLD_P_DIAG0 + PA_GGNET)
      call add_variable(reg, 'rough_site', 'aerodynamic roughness length', 'm',                  &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_DIAG0 + PA_ROUGH)
      call add_variable(reg, 'displace_site', 'zero-plane displacement height', 'm',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_DIAG0 + PA_DISPLACE)
      !----- Per patch too: with its canopy-air depth these are what move the forcing to the patch's  !
      !      canopy-air top (docs/science/forcing.md §8), so each patch's forcing can be rebuilt. ----!
      call add_variable(reg, 'rough_patch', 'patch aerodynamic roughness length', 'm',           &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_ROUGH)
      call add_variable(reg, 'displace_patch', 'patch zero-plane displacement height', 'm',      &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_DISPLACE)
      !--- carbon ---!
      call add_variable(reg, 'gpp_rate_site', 'gross primary productivity (mean rate)', 'umol/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, DAY_MON_YR, FLD_P_DIAG0 + PA_GPP)
      call add_variable(reg, 'nee_site', 'net ecosystem exchange (mean rate)', 'umol/m2/s',      &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, DAY_MON_YR, FLD_P_DIAG0 + PA_NEE)
      !--- water ---!
      call add_variable(reg, 'et_rate_site', 'evapotranspiration (mean rate)', 'kg/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PA_TRANSP)
      call add_variable(reg, 'precip_site', 'total precipitation (mean rate)', 'kg/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_PRECIP)
      !--- budget health (GRP_NUMERICS): the numbers that say whether anything above is real. ---!
      call add_variable(reg, 'resid_energy_site', 'mean signed whole-column energy-budget residual (+ = appearing)', 'W/m2', &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PA_RESID_ENERGY)
      call add_variable(reg, 'resid_water_site', 'mean signed whole-column water-budget residual (+ = appearing)', 'kg/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PA_RESID_WATER)
      !--- per-patch twins of the four that vary most between a gap and a closed canopy. ---!
      call add_variable(reg, 'le_patch', 'patch latent heat flux', 'W/m2',                       &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_LE)
      call add_variable(reg, 'h_patch', 'patch sensible heat flux', 'W/m2',                      &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_H)
      call add_variable(reg, 'rnet_patch', 'patch net all-wave radiation', 'W/m2',               &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PA_RNET)
      call add_variable(reg, 'nee_patch', 'patch net ecosystem exchange', 'umol/m2/s',           &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, DAY_MON, FLD_P_DIAG0 + PA_NEE)
   end subroutine register_patch_fluxes

   !=======================================================================================!
   !  FAST-tier (sub-daily) diagnostics. These read the STAGED per-sub-step samples, not live   !
   !  site state, so they carry the FAST bit ONLY -- putting one on a coarse tier would read a   !
   !  staging buffer that is not populated there.                                                 !
   !=======================================================================================!
   subroutine register_fast(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'cas_temp_fast', 'site canopy-air-space temperature (sub-daily)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_CAS_TEMP)
      call add_variable(reg, 'soil_temp_top_fast', 'site soil-top temperature (sub-daily)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SOIL_TEMP_TOP)
      call add_variable(reg, 'skin_temp_fast', 'site skin temperature (longwave up, sub-daily)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SKIN_TEMP)
      call add_variable(reg, 'gpp_rate_fast', 'site instantaneous GPP rate', 'umol/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_GPP)
      call add_variable(reg, 'le_flux_fast', 'site latent-heat (ET) flux', 'W/m2',               &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_PA0 + PA_LE)
      call add_variable(reg, 'h_flux_fast', 'site sensible-heat flux', 'W/m2',                   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_H)
      call add_variable(reg, 'rnet_fast', 'site net all-wave radiation', 'W/m2',                 &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_RNET)
      call add_variable(reg, 'sw_in_fast', 'incident shortwave at canopy top (forcing)', 'W/m2', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_SW_IN)
      call add_variable(reg, 'sw_up_fast', 'shortwave leaving the canopy top (VIS + NIR)', 'W/m2', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SW_UP)
      call add_variable(reg, 'lw_up_fast', 'longwave leaving the canopy top (emission included)', 'W/m2', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_LW_UP)
      call add_variable(reg, 'ustar_fast', 'friction velocity', 'm/s',                           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_USTAR)
      call add_variable(reg, 'air_temp_fast', 'reference-level forcing air temperature', 'K',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_TAIR)
      call add_variable(reg, 'nee_fast', 'net ecosystem exchange (+ to atmosphere)', 'umol/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_NEE)
      call add_variable(reg, 'npp_rate_fast', 'net primary productivity (GPP - maintenance Ra)',  &
                        'umol/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_NPP)
      call add_variable(reg, 'reco_fast', 'ecosystem respiration (autotrophic + heterotrophic)',  &
                        'umol/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_RECO)
      call add_variable(reg, 'cas_co2_fast', 'canopy-air CO2 mixing ratio', 'umol/mol',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_CAS_CO2)
      call add_variable(reg, 'atm_co2_fast', 'free-atmosphere CO2 (forcing)', 'umol/mol',         &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_CO2)
      call add_variable(reg, 'soil_temp_site_fast', 'each site''s area-weighted soil temperature (sub-daily)', 'K', &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_SOIL_TEMP)
      call add_variable(reg, 'soil_water_site_fast', 'each site''s area-weighted soil moisture (sub-daily)', 'm3/m3', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_SOIL_WATER)
      call add_variable(reg, 'leaf_temp_cohort_fast', 'per-cohort leaf temperature (sub-daily)', 'K', &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_COH_LEAF_TEMP)
      call add_variable(reg, 'gpp_cohort_fast', 'per-cohort GPP rate (sub-daily)', 'umol/plant/s', &
                        DIM_COHORT, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_COH_GPP)
      call add_variable(reg, 'height_cohort_fast', 'per-cohort height (tallest-cohort selection)', 'm', &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_COH_HEIGHT)
      !----- The leaf solve's inputs and output per cohort, sub-daily, as the solve saw them (rows of  !
      !      the per-cohort diagnostic block, SRC_F_CO0 + CO_*): enough to repeat the solve outside  !
      !      the model, which the fast calibration's canopy of leaf kernels does (calibrate_fast). ---!
      call add_variable(reg, 'gx_par_cohort_fast', 'leaf solve: incident-equivalent PAR per leaf area', &
                        'umol/m2/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_LEAF_PAR)
      call add_variable(reg, 'gx_leaf_temp_cohort_fast', 'leaf solve: leaf temperature', 'K',         &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_LEAF_TEMP)
      call add_variable(reg, 'gx_vpd_cohort_fast', 'leaf solve: leaf-to-air vapour pressure deficit', &
                        'Pa', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_LEAF_VPD)
      call add_variable(reg, 'gx_ca_cohort_fast', 'leaf solve: canopy-air CO2', 'umol/mol',           &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_CA)
      call add_variable(reg, 'gx_pressure_cohort_fast', 'leaf solve: air pressure', 'Pa',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_PRESSURE)
      call add_variable(reg, 'gx_psi_leaf_cohort_fast', 'leaf solve: leaf water potential', 'MPa',    &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_PSI_LEAF)
      call add_variable(reg, 'gx_psi_predawn_cohort_fast', 'leaf solve: predawn leaf water potential', &
                        'MPa', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_PSI_PREDAWN)
      call add_variable(reg, 'gx_gb_cohort_fast', 'leaf solve: boundary-layer conductance',           &
                        'mol/m2/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_GB_MOL)
      call add_variable(reg, 'gx_agross_cohort_fast', 'leaf solve: gross assimilation per leaf area', &
                        'umol/m2/s', DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, FAST_ONLY, SRC_F_CO0 + CO_AGROSS)
      call add_variable(reg, 'abs_par_cohort_fast', 'absorbed PAR (incident-equivalent, per ground area)', &
                        'W/m2', DIM_COHORT, AGG_TMEAN, GRP_RADIATION, FAST_ONLY, SRC_F_CO0 + CO_ABS_PAR)
      !----- Per-patch twins (docs/science/diagnostics.md §4, The FAST tier by patch): the same    !
      !      samples before the area-weighted sum over patches. A site mean over a closed canopy    !
      !      and a gap describes neither -- at midday the two can differ in sign -- so a patchy      !
      !      stand's sub-daily behaviour needs the patch axis. Each is its coarse patch variable's   !
      !      name with _fast (a quantity only this tier has keeps its site stem), and each follows   !
      !      [output].axes_patch or axes_soil_patch as the coarse ones do. --------------------------!
      call add_variable(reg, 'cas_temp_patch_fast', 'patch canopy-air-space temperature (sub-daily)', 'K', &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_CAS_TEMP)
      call add_variable(reg, 'soil_temp_top_patch_fast', 'patch soil-top temperature (sub-daily)', 'K', &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SOIL_TEMP_TOP)
      call add_variable(reg, 'skin_temp_patch_fast', 'patch skin temperature (longwave up, sub-daily)', 'K', &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SKIN_TEMP)
      call add_variable(reg, 'gpp_rate_patch_fast', 'patch instantaneous GPP rate', 'umol/m2/s',  &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_GPP)
      call add_variable(reg, 'le_patch_fast', 'patch latent-heat (ET) flux', 'W/m2',              &
                        DIM_PATCH, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_PA0 + PA_LE)
      call add_variable(reg, 'h_patch_fast', 'patch sensible-heat flux', 'W/m2',                  &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_H)
      call add_variable(reg, 'rnet_patch_fast', 'patch net all-wave radiation', 'W/m2',           &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_RNET)
      call add_variable(reg, 'sw_up_patch_fast', 'patch shortwave leaving the canopy top (VIS + NIR)', 'W/m2', &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_SW_UP)
      call add_variable(reg, 'lw_up_patch_fast', 'patch longwave leaving the canopy top (emission included)', 'W/m2', &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_LW_UP)
      call add_variable(reg, 'ustar_patch_fast', 'patch friction velocity', 'm/s',                &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_PA0 + PA_USTAR)
      call add_variable(reg, 'nee_patch_fast', 'patch net ecosystem exchange (+ to atmosphere)', 'umol/m2/s', &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_NEE)
      call add_variable(reg, 'npp_rate_patch_fast', 'patch net primary productivity (GPP - maintenance Ra)', 'umol/m2/s', &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_NPP)
      call add_variable(reg, 'reco_patch_fast', 'patch ecosystem respiration (autotrophic + heterotrophic)', 'umol/m2/s', &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_RECO)
      call add_variable(reg, 'cas_co2_patch_fast', 'patch canopy-air CO2 mixing ratio', 'umol/mol', &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_PA0 + PA_CAS_CO2)
      call add_variable(reg, 'soil_temp_layer_patch_fast', 'soil temperature by layer and patch (sub-daily)', 'K', &
                        DIM_SOIL_PATCH, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_SOIL_TEMP)
      call add_variable(reg, 'soil_water_layer_patch_fast', 'volumetric soil moisture by layer and patch (sub-daily)', 'm3/m3', &
                        DIM_SOIL_PATCH, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_SOIL_WATER)
   end subroutine register_fast

   !=======================================================================================!
   !  The FORCING the run used (docs/science/diagnostics.md §5): the atmospheric boundary after the !
   !  reader's shortwave partition, rain/snow split and optional height and lapse corrections --   !
   !  write-only provenance. The coarse tiers are dt-weighted means over the fast sub-steps, kept   !
   !  once per polygon (poly%diag, PY_*) with sw_in_site and precip_site, which carry the total     !
   !  shortwave and precipitation; the FAST tier is each record's mean of the sub-step samples,     !
   !  which is what checks the diurnal reconstruction against a tower. cosz is 0 while the sun is   !
   !  down, so a daily mean counts the night as zero.                                               !
   !=======================================================================================!
   subroutine register_forcing(reg)
      type(output_registry_t), intent(inout) :: reg
      !----- Each site's forcing: the polygon's after the site's terrain (docs/science/forcing.md  !
      !      §8), and the growth temperature its leaves acclimate to. -----------------------------!
      call add_variable(reg, 'air_temp_by_site', 'each site''s air temperature at the forcing''s height', 'K', &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_TAIR)
      call add_variable(reg, 'qair_by_site', 'each site''s specific humidity', 'kg/kg',            &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_QAIR)
      call add_variable(reg, 'psurf_by_site', 'each site''s surface pressure', 'Pa',               &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_PSURF)
      call add_variable(reg, 'rho_air_by_site', 'each site''s air density', 'kg/m3',               &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_RHO_AIR)
      call add_variable(reg, 'lwdown_by_site', 'each site''s downward longwave', 'W/m2',           &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_LWDOWN)
      call add_variable(reg, 'sw_in_by_site', 'each site''s incident shortwave, its beam on its slope', 'W/m2', &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_SW_IN)
      call add_variable(reg, 'par_beam_by_site', 'each site''s direct-beam PAR on its slope', 'W/m2', &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_PAR_BEAM)
      call add_variable(reg, 'nir_beam_by_site', 'each site''s direct-beam NIR on its slope', 'W/m2', &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_DIAG0 + SI_NIR_BEAM)
      call add_variable(reg, 'growth_temp_by_site', 'each site''s running-mean growth temperature', 'K', &
                        DIM_SITE, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_SI_GROWTH_TEMP)
      call add_variable(reg, 'air_temp_site', 'air temperature at the forcing''s own height (forcing)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_TAIR)
      call add_variable(reg, 'qair_site', 'specific humidity at the forcing''s own height (forcing)', 'kg/kg', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_QAIR)
      call add_variable(reg, 'psurf_site', 'surface pressure (forcing)', 'Pa',                    &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_PSURF)
      call add_variable(reg, 'wind_site', 'wind speed at the forcing''s own height (forcing)', 'm/s',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_WIND)
      call add_variable(reg, 'lwdown_site', 'downward longwave at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_LWDOWN)
      call add_variable(reg, 'par_beam_site', 'direct-beam PAR at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_PAR_BEAM)
      call add_variable(reg, 'par_diffuse_site', 'diffuse PAR at canopy top (forcing)', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_PAR_DIFFUSE)
      call add_variable(reg, 'nir_beam_site', 'direct-beam NIR at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_NIR_BEAM)
      call add_variable(reg, 'nir_diffuse_site', 'diffuse NIR at canopy top (forcing)', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_NIR_DIFFUSE)
      call add_variable(reg, 'snowfall_site', 'frozen precipitation (mean rate, split at each patch''s canopy-air top)', &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_SNOWFALL)
      call add_variable(reg, 'atm_co2_site', 'free-atmosphere CO2 (forcing)', 'umol/mol',         &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_CO2)
      call add_variable(reg, 'cosz_site', 'cosine of the solar zenith angle (0 with the sun down)', '1', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_COSZ)
      call add_variable(reg, 'rho_air_site', 'air density at the forcing''s own height', 'kg/m3',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, DAY_MON_YR, FLD_PY_DIAG0 + PY_RHO_AIR)
      !----- The patch's own forcing: the sample moved to its canopy-air top                         !
      !      (docs/science/forcing.md §8). Only the wind and the air temperature differ between       !
      !      patches; the rest is the polygon's, above. -----------------------------------------------!
      call add_variable(reg, 'wind_cas_top_patch', 'wind speed at the patch canopy-air top', 'm/s', &
                        DIM_PATCH, AGG_TMEAN, GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_WIND_CAS_TOP)
      call add_variable(reg, 'air_temp_cas_top_patch', 'air temperature at the patch canopy-air top', 'K', &
                        DIM_PATCH, AGG_TMEAN, GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_TAIR_CAS_TOP)
      call add_variable(reg, 'cas_depth_patch', 'patch canopy-air depth (the top its forcing is moved to)', &
                        'm', DIM_PATCH, AGG_TMEAN, GRP_FORCING, DAY_MON, FLD_P_DIAG0 + PA_Z_CAS_TOP)
      !----- The FAST tier: each sub-step's own sample (air_temp_fast, sw_in_fast and atm_co2_fast  !
      !      are registered with the fast fluxes, in this group). --------------------------------------!
      call add_variable(reg, 'qair_fast', 'specific humidity at the forcing''s own height (forcing)', 'kg/kg', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_QAIR)
      call add_variable(reg, 'psurf_fast', 'surface pressure (forcing)', 'Pa',                    &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_PSURF)
      call add_variable(reg, 'wind_fast', 'wind speed at the forcing''s own height (forcing)', 'm/s',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_WIND)
      call add_variable(reg, 'lwdown_fast', 'downward longwave at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_LWDOWN)
      call add_variable(reg, 'par_beam_fast', 'direct-beam PAR at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_PAR_BEAM)
      call add_variable(reg, 'par_diffuse_fast', 'diffuse PAR at canopy top (forcing)', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_PAR_DIFFUSE)
      call add_variable(reg, 'nir_beam_fast', 'direct-beam NIR at canopy top (forcing)', 'W/m2',  &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_NIR_BEAM)
      call add_variable(reg, 'nir_diffuse_fast', 'diffuse NIR at canopy top (forcing)', 'W/m2',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_NIR_DIFFUSE)
      call add_variable(reg, 'precip_fast', 'total precipitation (mean rate, forcing)', 'kg/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_PRECIP)
      call add_variable(reg, 'rainf_fast', 'liquid precipitation (mean rate, split at each patch''s canopy-air top)', &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_RAINF)
      call add_variable(reg, 'snowfall_fast', 'frozen precipitation (mean rate, split at each patch''s canopy-air top)', &
                        'kg/m2/s', DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_SNOWFALL)
      call add_variable(reg, 'cosz_fast', 'cosine of the solar zenith angle (0 with the sun down)', '1', &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_COSZ)
      call add_variable(reg, 'rho_air_fast', 'air density at the forcing''s own height', 'kg/m3',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_FORCING, FAST_ONLY, SRC_F_PY0 + PY_RHO_AIR)
   end subroutine register_forcing

   !----- Append one descriptor (the single "add a variable" registry edit).                   !
   !                                                                                          !
   !      w / mn / sc are the AGGREGATION contract (see var_desc_t): the per-cohort weight kind, !
   !      whether the reduction is a weighted MEAN (intensive) or a weighted SUM (extensive), and  !
   !      any unit conversion. Defaults (W_NONE, sum, 1.0) are right for a raw slab, which is why  !
   !      the cohort/patch rows above mostly omit them.                                            !
   subroutine add_variable(reg, nm, ln, un, dm, ag, gr, st, src, xt, w, mn, sc)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: nm, ln, un
      integer(ik),             intent(in)    :: dm, ag, gr, st, src
      integer(ik), optional,   intent(in)    :: xt, w
      logical,     optional,   intent(in)    :: mn
      real(wp),    optional,   intent(in)    :: sc
      integer(ik) :: k
      if (reg%nvar >= MAX_OUTPUT_VARS) error stop 'meds_output_registry: MAX_OUTPUT_VARS exceeded'
      if (find_var_index(reg, nm) > 0_ik)                                                        &
         error stop 'meds_output_registry: duplicate variable name ('//trim(nm)//')'
      reg%nvar = reg%nvar + 1_ik ; k = reg%nvar
      reg%var(k)%name            = nm
      reg%var(k)%long_name       = ln
      reg%var(k)%units           = un
      reg%var(k)%dim             = dm
      reg%var(k)%agg             = ag
      reg%var(k)%group           = gr
      reg%var(k)%streams         = st
      reg%var(k)%streams_default = st
      reg%var(k)%enabled         = .true.
      reg%var(k)%source_id       = src
      reg%var(k)%xtype           = XTYPE_DOUBLE ; if (present(xt)) reg%var(k)%xtype = xt
      reg%var(k)%weight          = W_NONE       ; if (present(w))  reg%var(k)%weight = w
      reg%var(k)%mean            = .false.      ; if (present(mn)) reg%var(k)%mean   = mn
      reg%var(k)%scale           = 1.0_wp       ; if (present(sc)) reg%var(k)%scale  = sc
   end subroutine add_variable

   !----- Axis index (1..N_AXIS) of a DIM_* code, for the [output].axes_* toggles. DIM_SCALAR   !
   !      has no axis toggle (a site scalar is always available) and returns 0.                  !
   pure integer(ik) function dim_axis_index(dm) result(a)
      integer(ik), intent(in) :: dm
      select case (dm)
      case (DIM_COHORT)     ; a = 1_ik
      case (DIM_PATCH)      ; a = 2_ik
      case (DIM_PFT)        ; a = 3_ik
      case (DIM_SIZE)       ; a = 4_ik
      case (DIM_SOIL_PATCH) ; a = 5_ik
      case default          ; a = 0_ik
      end select
   end function dim_axis_index

   !----- Step 2: suppress a whole AXIS without naming variables. With ~66 cohort-dimensioned   !
   !      variables (docs/science/diagnostics.md §7), "site level only, nothing per-cohort" is the !
   !      most common request; it is also the single biggest lever on output volume.              !
   subroutine apply_axis_toggles(reg, out_cfg)
      type(output_registry_t), intent(inout) :: reg
      type(output_config_t),   intent(in)    :: out_cfg
      integer(ik) :: k, a
      do k = 1_ik, reg%nvar
         a = dim_axis_index(reg%var(k)%dim)
         if (a >= 1_ik .and. a <= N_AXIS) then
            if (.not. out_cfg%axis_on(a)) reg%var(k)%streams = FREQ_NONE
         end if
      end do
   end subroutine apply_axis_toggles

   !----- Guard: a cohort/patch-dimensioned variable may NOT be on the annual stream          !
   !      (docs/science/diagnostics.md §3). ------------------------------------------------!
   subroutine enforce_annual_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         if (ragged_dim(reg%var(k)%dim) .and. iand(reg%var(k)%streams_default, FREQ_ANNUAL) /= 0_ik) &
            error stop 'meds_output_registry: cohort/patch variable on the annual stream ('//trim(reg%var(k)%name)//')'
      end do
   end subroutine enforce_annual_guard

   !----- Guard: the FAST tier and the FAST-staged sources are a CLOSED PAIR, both directions.  !
   !                                                                                          !
   !      Sub-daily resolution exists only inside the fast loop's sub-step, so a FAST record is   !
   !      replayed from the STAGED samples (output_integrate_fast), never from live site state --  !
   !      which by replay time is the end-of-slow-step snapshot. Two things follow, and neither    !
   !      announces itself if violated:                                                            !
   !                                                                                          !
   !        * a variable with a NON-fast source on the FAST tier resolves against a staging buffer  !
   !          that was never filled for it, and writes _FillValue for the whole run;                 !
   !        * a FAST-staged source on a COARSE tier resolves against live site state, which does not !
   !          carry it, and does the same.                                                            !
   !                                                                                          !
   !      Both are reachable from meds_io_config.toml -- `cas_temp_site = "F"` is exactly the kind of !
   !      thing a user writes, and unguarded it would silently produce a file full of fill. Rejecting !
   !      the pairing at start-up, naming the offender and the variable they probably meant, is the   !
   !      whole point: a diagnostic that quietly reports nothing is worse than one that refuses to     !
   !      start. The meds_io_config.toml overrides are applied after this, so apply_variable_override !
   !      runs the same check on each stream mask it sets.                                             !
   subroutine enforce_fast_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         call check_fast_tier(reg%var(k)%name, reg%var(k)%source_id, reg%var(k)%streams)
      end do
   end subroutine enforce_fast_guard

   !----- One variable's streams against its source: the FAST tier <-> a FAST-staged source. ---!
   subroutine check_fast_tier(name, source_id, streams)
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: source_id, streams
      logical :: is_fast_src
      is_fast_src = (src_class(source_id) == SRCK_FAST)
      if (iand(streams, FREQ_FAST) /= 0_ik .and. .not. is_fast_src)                              &
         error stop 'meds_output_registry: "'//trim(name)//'" is not a FAST-staged '//            &
                    'variable and cannot be put on the F stream (look for a *_fast twin)'
      if (is_fast_src .and. iand(streams, not(FREQ_FAST)) /= 0_ik)                                &
         error stop 'meds_output_registry: "'//trim(name)//'" is FAST-staged and can '//           &
                    'only be on the F stream (the coarse tiers read live state, not the staging)'
   end subroutine check_fast_tier

   !----- Step 3: disable every variable in a group whose main-config toggle is .false. -------!
   subroutine apply_group_toggles(reg, out_cfg)
      type(output_registry_t), intent(inout) :: reg
      type(output_config_t),   intent(in)    :: out_cfg
      integer(ik) :: k, g
      do k = 1_ik, reg%nvar
         g = reg%var(k)%group
         if (g >= 1_ik .and. g <= N_GRP) then
            if (.not. out_cfg%grp_on(g)) reg%var(k)%streams = FREQ_NONE
         end if
      end do
   end subroutine apply_group_toggles

   !----- Step 4: a disabled tier clears its FREQ_* bit from every variable's mask. -----------!
   subroutine apply_freq_enables(reg, out_cfg)
      type(output_registry_t), intent(inout) :: reg
      type(output_config_t),   intent(in)    :: out_cfg
      integer(ik) :: k, t, bit
      do t = 1_ik, N_FREQ
         if (.not. out_cfg%freq_on(t)) then
            bit = freq_bit(t)
            do k = 1_ik, reg%nvar
               reg%var(k)%streams = iand(reg%var(k)%streams, not(bit))
            end do
         end if
      end do
   end subroutine apply_freq_enables

   !----- Step 5 primitive: apply one per-variable override (from meds_io_config.toml).        !
   !      found=.false. => the name matched no registry variable (the unknown-key trap caller). !
   subroutine apply_variable_override(reg, name, kind, mask, found)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: name
      integer(ik),             intent(in)    :: kind   !< OVR_TRUE | OVR_FALSE | OVR_MASK
      integer(ik),             intent(in)    :: mask   !< used for OVR_MASK
      logical,                 intent(out)   :: found
      integer(ik) :: k, m
      k = find_var_index(reg, name)
      found = (k > 0_ik)
      if (.not. found) return
      select case (kind)
      case (OVR_TRUE)  ; reg%var(k)%streams = reg%var(k)%streams_default
      case (OVR_FALSE) ; reg%var(k)%streams = FREQ_NONE
      case (OVR_MASK)
         m = mask
         if (ragged_dim(reg%var(k)%dim) .and. iand(m, FREQ_ANNUAL) /= 0_ik)                      &
            error stop 'meds_output_registry: annual (Y) stream forbidden on cohort/patch variable ('//trim(name)//')'
         call check_fast_tier(reg%var(k)%name, reg%var(k)%source_id, m)
         reg%var(k)%streams = m
      end select
   end subroutine apply_variable_override

   !----- 1-based registry index of a variable by name; 0 if absent. -------------------------!
   pure integer(ik) function find_var_index(reg, name) result(idx)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik) :: k
      idx = 0_ik
      do k = 1_ik, reg%nvar
         if (trim(reg%var(k)%name) == trim(name)) then ; idx = k ; return ; end if
      end do
   end function find_var_index

   !----- Parse a "F D M Y" stream-string into an ior(FREQ_*) mask (the value grammar of      !
   !      docs/science/diagnostics.md §5).                                                   !
   !      status = 0 ok; status = 1 unrecognized token (the offender returned in bad_token).    !
   subroutine parse_stream_mask(s, mask, status, bad_token)
      character(len=*), intent(in)  :: s
      integer(ik),      intent(out) :: mask, status
      character(len=*), intent(out) :: bad_token
      integer(ik) :: i
      character(len=1) :: c
      mask = FREQ_NONE ; status = 0_ik ; bad_token = ''
      do i = 1_ik, int(len_trim(s), ik)
         c = s(i:i)
         select case (c)
         case ('F','f') ; mask = ior(mask, FREQ_FAST)
         case ('D','d') ; mask = ior(mask, FREQ_DAILY)
         case ('M','m') ; mask = ior(mask, FREQ_MONTHLY)
         case ('Y','y') ; mask = ior(mask, FREQ_ANNUAL)
         case (' ')     ! separator
         case default   ; status = 1_ik ; bad_token = c ; return
         end select
      end do
   end subroutine parse_stream_mask

   !=======================================================================================!
   !  Build a file set's REGISTRY + config (netCDF-free; NO buffers yet). Split from           !
   !  allocation so a caller can apply meds_io_config.toml per-variable overrides to files%reg    !
   !  (which change which (var,tier) buffers are needed) BEFORE manager_finalize and             !
   !  manager_alloc_buffers.                                                                      !
   !=======================================================================================!
   subroutine manager_setup(files, cfg)
      type(output_files_t),   intent(out) :: files
      type(meds_config_t),    intent(in)  :: cfg
      integer(ik) :: ne
      files%enabled = cfg%output%enabled
      files%fast_loop_on = cfg%fast_biophysics_on
      call build_output_registry(files%reg, cfg)
      files%cohort_max = cfg%output%cohort_max
      files%patch_max  = cfg%output%patch_max
      files%dir        = cfg%output%dir
      files%prefix     = cfg%output%prefix
      files%file_chunk = cfg%output%file_chunk
      files%sync_every = cfg%output%sync_every
      files%fast_interval_steps = cfg%output%fast_interval_steps
      if (cfg%forcing%forcing_on .and. cfg%forcing%format == MET_BACKEND_ED_ERA5LAND)                &
         files%forcing_qair = 'computed by MEDS from the ED_ERA5land 2 m dewpoint Td and surface '//         &
                            'pressure P: q = 0.622 e / (P - 0.378 e), e = e_sat(Td), the Bolton (1980) '// &
                            'liquid-water saturation vapour pressure of meds_therm_lib'

      !----- Run-dependent parameters the DERIVED diagnostics need (diag_setup_t). n_pft is the   !
      !      RUN-TIME PFT count, so the netCDF `pft` dimension is run-dependent -- which is why      !
      !      the serializer also writes a `pft` coordinate variable, so a file stays self-describing !
      !      when compared across runs with different PFT tables.                                    !
      files%diag%n_pft  = cfg%pft%n
      !----- DBH size classes: TOML edges if given, else the built-in inventory set. -----------!
      if (cfg%output%n_dbh_class > 0_ik) then
         files%diag%n_dbh_class = min(cfg%output%n_dbh_class, MAX_DBH_CLASS)
         ne = files%diag%n_dbh_class + 1_ik
         files%diag%dbh_edges(1:ne) = cfg%output%dbh_class_edges(1:ne)
      else
         files%diag%n_dbh_class = N_DBH_CLASS_DEFAULT
         files%diag%dbh_edges(1:N_DBH_CLASS_DEFAULT + 1_ik) = DBH_EDGES_DEFAULT
      end if
      call check_dbh_edges(files%diag%dbh_edges, files%diag%n_dbh_class)
      !----- attach_output sets soil_ready on a polygon's buffers when its fast loop runs; until !
      !      then the psi/wetness diagnostics emit _FillValue. ----------------------------------!
      files%diag%soil_ready = .false.
   end subroutine manager_setup

   !----- Largest slab length any LIVE variable can produce (see manager_finalize). ------------!
   pure integer(ik) function live_max_slab(files) result(cap)
      type(output_files_t),   intent(in) :: files
      integer(ik) :: k
      cap = 1_ik
      do k = 1_ik, files%reg%nvar
         if (.not. files%reg%var(k)%enabled)      cycle
         if (files%reg%var(k)%streams == FREQ_NONE) cycle
         cap = max(cap, dim_capacity(files, files%reg%var(k)%dim))
      end do
   end function live_max_slab

   !----- Slab capacity of one axis. ---------------------------------------------------------!
   pure integer(ik) function dim_capacity(files, dm) result(cap)
      type(output_files_t),   intent(in) :: files
      integer(ik),            intent(in) :: dm
      select case (dm)
      case (DIM_COHORT)     ; cap = files%cohort_max
      case (DIM_PATCH)      ; cap = files%patch_max
      case (DIM_SOIL)       ; cap = files%site_max * n_soil_layer_max
      case (DIM_PFT)        ; cap = max(files%diag%n_pft, 1_ik)
      case (DIM_SIZE)       ; cap = max(files%diag%n_dbh_class, 1_ik)
      case (DIM_SOIL_PATCH) ; cap = files%patch_max * n_soil_layer_max
      case (DIM_SITE)
         cap = files%site_max
      case default          ; cap = 0_ik
      end select
   end function dim_capacity

   !----- The size-class edges must be strictly ascending, or dbh_class_index silently mis-bins  !
   !      (and sum_class == site would quietly stop holding). Fail at start-up instead.          !
   subroutine check_dbh_edges(edges, n_class)
      real(wp),    intent(in) :: edges(:)
      integer(ik), intent(in) :: n_class
      integer(ik) :: j
      do j = 1_ik, n_class
         if (.not. (edges(j+1_ik) > edges(j)))                                                   &
            error stop 'meds_output_registry: [output].dbh_class_edges must be strictly ascending'
      end do
   end subroutine check_dbh_edges

   !----- Finish the file set once the registry is FINALIZED (after manager_setup and any          !
   !      per-variable overrides): the per-tier stream handles and the slab size. -------------------!
   subroutine manager_finalize(files)
      type(output_files_t),  intent(inout) :: files
      integer(ik) :: t
      do t = 1_ik, N_FREQ
         files%stream(t)%freq  = freq_bit(t)
         allocate(files%stream(t)%vid(files%reg%nvar)) ; files%stream(t)%vid = -1_ik
      end do
      !----- SLAB SIZING. Size the shared pending-record slab to the largest axis that is        !
      !      ACTUALLY LIVE, not to the largest axis that exists. With ~66 cohort-dimensioned        !
      !      variables the pending records are the dominant memory term, and a run with            !
      !      axes_cohort = false (or the 2-D soil axis off, the default) should not pay for the     !
      !      axis it switched off. It is computed HERE, after every override: an [output].io_config  !
      !      can switch on a slab variable that manager_setup left off, and a slab sized before     !
      !      that is too short for it. ---------------------------------------------------------!
      files%max_slab = live_max_slab(files)
      !----- Whether the cohort and patch caps bind: only a live variable on that axis fills a slab of  !
      !      that size. ------------------------------------------------------------------------------!
      files%cohort_axis = live_axis(files, [DIM_COHORT])
      files%patch_axis  = live_axis(files, [DIM_PATCH, DIM_SOIL_PATCH])
   end subroutine manager_finalize

   !----- Does any live variable have one of these axes? ---------------------------------------!
   pure logical function live_axis(files, dims) result(yes)
      type(output_files_t), intent(in) :: files
      integer(ik),          intent(in) :: dims(:)
      integer(ik) :: k
      yes = .false.
      do k = 1_ik, files%reg%nvar
         if (.not. files%reg%var(k)%enabled)        cycle
         if (files%reg%var(k)%streams == FREQ_NONE) cycle
         if (any(dims == files%reg%var(k)%dim)) then ; yes = .true. ; return ; end if
      end do
   end function live_axis

   !----- Allocate one polygon's buffers for a finalized file set: the integrator buffers of       !
   !      every live (variable, tier) pair and the per-tier scratch records. ------------------------!
   subroutine manager_alloc_buffers(files, bufs)
      type(output_files_t),  intent(in)  :: files
      type(output_buffers_t), intent(out) :: bufs
      integer(ik) :: t, k, nv, cap, bit
      nv = files%reg%nvar
      allocate(bufs%buf(nv, N_FREQ))
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, nv
            if (files%reg%var(k)%enabled .and. iand(files%reg%var(k)%streams, bit) /= 0_ik) then
               cap = dim_capacity(files, files%reg%var(k)%dim)
               call alloc_integ_buffer(bufs%buf(k,t), files%reg%var(k), bit, cap)
               bufs%buf(k,t)%var_id = k
            end if
         end do
      end do
      do t = 1_ik, N_FREQ
         allocate(bufs%pending(t)%sval(nv), bufs%pending(t)%svalid(nv), bufs%pending(t)%nslab(nv))
         allocate(bufs%pending(t)%slab(files%max_slab, nv), bufs%pending(t)%slabvalid(files%max_slab, nv))
         bufs%pending(t)%used = .false.
      end do
      bufs%fast_on = files%enabled .and. files%reg%nidx(1) > 0_ik
      bufs%diag = files%diag
   end subroutine manager_alloc_buffers

   !----- Convenience: registry, stream handles and one polygon's buffers in one call (no       !
   !      overrides). ---------------------------------------------------------------------------!
   subroutine manager_alloc(files, bufs, cfg)
      type(output_files_t),  intent(out) :: files
      type(output_buffers_t), intent(out) :: bufs
      type(meds_config_t),   intent(in)  :: cfg
      call manager_setup(files, cfg)
      call manager_set_polygons(files, [0_ik], [1_ik], [0.0_wp], [0.0_wp], [0_ik], [0_ik], .false.)
      call manager_finalize(files)
      call manager_alloc_buffers(files, bufs)
   end subroutine manager_alloc

   !----- The file set's polygon and site axes: one row per polygon, its id, location, grid indices !
   !      and site count; the site axis is the largest polygon's. The rest of a row -- its sites'   !
   !      areas, terrain and soil columns -- is filled when the polygon attaches to the set         !
   !      (attach_output); a polygon's missing sites keep the fill value. --------------------------!
   subroutine manager_set_polygons(files, polygon_id, n_site, lat, lon, row, col, with_chave_e)
      type(output_files_t), intent(inout) :: files
      integer(ik),          intent(in)    :: polygon_id(:), n_site(:), row(:), col(:)
      real(wp),             intent(in)    :: lat(:), lon(:)
      logical,              intent(in)    :: with_chave_e   !< the Chave height allometry's site index
      integer(ik) :: n, ns
      n  = size(polygon_id, kind=ik)
      ns = max(maxval(n_site), 1_ik)
      files%n_polygon      = n
      files%site_max       = ns
      files%polygon_id     = polygon_id
      files%polygon_n_site = n_site
      files%polygon_row    = row
      files%polygon_col    = col
      files%polygon_lat    = lat
      files%polygon_lon    = lon
      allocate(files%site_area(ns, n), files%site_elevation(ns, n), files%site_slope(ns, n),        &
               files%site_aspect(ns, n), files%site_soil_depth(ns, n), files%site_drain_factor(ns, n), &
               files%site_deep_temp(ns, n), source=MISSING_VALUE)
      allocate(files%site_n_soil(ns, n), source=0_ik)
      allocate(files%site_bottom_bc(ns, n), source=MISSING_INT)
      allocate(files%site_soil_z(n_soil_layer_max, ns, n), source=0.0_wp)
      allocate(files%site_soil_face(n_soil_layer_max + 1_ik, ns, n), source=0.0_wp)
      if (with_chave_e) allocate(files%polygon_chave_e(n), source=0.0_wp)
   end subroutine manager_set_polygons

   !----- Last: precompute, per tier, the list of live (enabled + in-tier) variable indices. -!
   subroutine build_freq_index(reg)
      type(output_registry_t), intent(inout) :: reg
      integer(ik) :: t, k, bit
      reg%idx_freq = 0_ik ; reg%nidx = 0_ik
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, reg%nvar
            if (reg%var(k)%enabled .and. iand(reg%var(k)%streams, bit) /= 0_ik) then
               reg%nidx(t) = reg%nidx(t) + 1_ik
               reg%idx_freq(reg%nidx(t), t) = k
            end if
         end do
      end do
   end subroutine build_freq_index

   !=======================================================================================!
   !  Turn the per-cohort / per-patch fast diagnostic blocks ON iff the finalized registry has     !
   !  at least one LIVE variable that reads them. This is what makes the capture free for a run     !
   !  that does not report it: `active = .false.` leaves both blocks unallocated and every entry     !
   !  point in meds_diag_vars a no-op, and the fast loop never even asks the leaf kernel for         !
   !  the extra flux fields.                                                                         !
   !=======================================================================================!
   subroutine activate_site_diag(files, poly)
      type(output_files_t),   intent(in)    :: files
      type(polygon_t),        intent(inout) :: poly
      logical     :: need_c, need_p, need_s, need_y, need_si
      integer(ik) :: t, j, k, src
      need_c = .false. ; need_p = .false. ; need_s = .false. ; need_y = .false.
      need_si = .false.
      if (files%enabled) then
         do t = 1_ik, N_FREQ
            do j = 1_ik, files%reg%nidx(t)
               k   = files%reg%idx_freq(j, t)
               src = files%reg%var(k)%source_id
               if (src > FLD_C_DIAG0 .and. src <= FLD_C_DIAG0 + N_CO_FAST_ROWS) need_c = .true.
               if (src > SRC_F_CO0 .and. src <= SRC_F_CO0 + N_CO_FAST_ROWS) need_c = .true.
               if (src > FLD_P_DIAG0 .and. src <= FLD_P_DIAG0 + N_PA_ROWS) need_p = .true.
               if (src > FLD_C_DIAG0 + N_CO_FAST_ROWS .and. src <= FLD_C_DIAG0 + N_CO_ROWS) need_s = .true.
               if (src > FLD_PY_DIAG0 .and. src <= FLD_PY_DIAG0 + N_PY_ROWS) need_y = .true.
               if (src > FLD_SI_DIAG0 .and. src <= FLD_SI_DIAG0 + N_SI_ROWS) need_si = .true.
            end do
         end do
      end if
      !----- Only the fast loop fills the per-cohort and polygon blocks, and the fast rows of the     !
      !      patch block: in a slow-only run they stay off, and their variables read as missing. -----!
      need_c = need_c .and. files%fast_loop_on ; need_y = need_y .and. files%fast_loop_on
      need_si = (need_si .and. files%fast_loop_on) .or. poly%site%diag%active
      poly%patch%diag%fast_rows = files%fast_loop_on
      !----- Additive: a block already on stays on. -----------------------------------------------!
      need_c = need_c .or. poly%cohort%diag%active  ; need_s = need_s .or. poly%cohort%sdiag%active
      need_p = need_p .or. poly%patch%diag%active   ; need_y = need_y .or. poly%diag%active
      call alloc_cohort_diag(poly%cohort%diag,  max(poly%cohort%cap, 1_ik), need_c)
      call alloc_cohort_diag(poly%cohort%sdiag, max(poly%cohort%cap, 1_ik), need_s, slow_rows=.true.)
      call alloc_patch_diag (poly%patch%diag,   max(poly%patch%cap,  1_ik), need_p)
      poly%diag%active = need_y                    ! the polygon block: fixed size, nothing to allocate
      call alloc_site_diag(poly%site%diag, poly%site%n, need_si)
   end subroutine activate_site_diag

   !=======================================================================================!
   !  Write a fully-commented meds_io_config.toml listing EVERY registered variable with its    !
   !  default stream mask, group and axis, ready to uncomment.                                  !
   !                                                                                          !
   !  This is how a user learns which names the per-variable overrides accept. It doubles as     !
   !  the regression guard that the registry and its documentation cannot drift apart -- the    !
   !  listing is GENERATED from the registry, so it is right by construction.                    !
   !=======================================================================================!
   subroutine dump_io_config(reg, path)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: path
      integer     :: u
      integer(ik) :: k
      character(len=52) :: lhs      !< left-hand side, blank-padded so the comments line up
      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(a)') '# MEDS per-variable output overrides -- GENERATED by `meds_main --dump-io-config`.'
      write(u,'(a)') '#'
      write(u,'(a)') '# Point the main config at this file with [output].io_config = "meds_io_config.toml".'
      write(u,'(a)') '# Entries here are the FINEST granularity and OVERRIDE the axis / group / per-tier'
      write(u,'(a)') '# switches in the main config.'
      write(u,'(a)') '#'
      write(u,'(a)') '#   name = true      -> force ON at the variable''s registry-default streams'
      write(u,'(a)') '#   name = false     -> force OFF everywhere'
      write(u,'(a)') '#   name = "F D M Y" -> replace the stream mask (F fast, D daily, M monthly, Y annual;'
      write(u,'(a)') '#                       Y is forbidden on a cohort- or patch-dimensioned variable)'
      write(u,'(a)') '#'
      write(u,'(a)') '# A name that matches no registered variable is a hard error, not a silent no-op.'
      write(u,'(a)') '# Every line below is commented out, so this file as generated changes nothing.'
      write(u,'(a)') ''
      write(u,'(a)') '[variables]'
      do k = 1_ik, reg%nvar
         lhs = '# '//trim(reg%var(k)%name)//' = "'//trim(mask_to_string(reg%var(k)%streams_default))//'"'
         write(u,'(a)') lhs//'  # ['//trim(group_name(reg%var(k)%group))//', '//                  &
                        trim(dim_name(reg%var(k)%dim))//'] '//trim(reg%var(k)%long_name)//        &
                        ' ('//trim(reg%var(k)%units)//')'
      end do
      close(u)
      write(*,'(2a)') ' wrote ', trim(path)
   end subroutine dump_io_config

   pure function mask_to_string(mask) result(s)
      integer(ik), intent(in) :: mask
      character(len=8) :: s
      integer :: p
      s = '' ; p = 0
      if (iand(mask, FREQ_FAST)    /= 0_ik) then ; p = p+1 ; s(p:p) = 'F' ; p = p+1 ; s(p:p) = ' ' ; end if
      if (iand(mask, FREQ_DAILY)   /= 0_ik) then ; p = p+1 ; s(p:p) = 'D' ; p = p+1 ; s(p:p) = ' ' ; end if
      if (iand(mask, FREQ_MONTHLY) /= 0_ik) then ; p = p+1 ; s(p:p) = 'M' ; p = p+1 ; s(p:p) = ' ' ; end if
      if (iand(mask, FREQ_ANNUAL)  /= 0_ik) then ; p = p+1 ; s(p:p) = 'Y' ; p = p+1 ; s(p:p) = ' ' ; end if
   end function mask_to_string

   pure function group_name(g) result(s)
      integer(ik), intent(in) :: g
      character(len=10) :: s
      select case (g)
      case (GRP_STRUCTURE)  ; s = 'structure'
      case (GRP_CARBON)     ; s = 'carbon'
      case (GRP_WATER)      ; s = 'water'
      case (GRP_ENERGY)     ; s = 'energy'
      case (GRP_NUMERICS)   ; s = 'numerics'
      case (GRP_RADIATION)  ; s = 'radiation'
      case (GRP_ECOPHYS)    ; s = 'ecophys'
      case (GRP_BIOGEOCHEM) ; s = 'biogeochem'
      case (GRP_FORCING)    ; s = 'forcing'
      case default          ; s = '?'
      end select
   end function group_name

   pure function dim_name(d) result(s)
      integer(ik), intent(in) :: d
      character(len=10) :: s
      select case (d)
      case (DIM_SCALAR)     ; s = 'site'
      case (DIM_COHORT)     ; s = 'cohort'
      case (DIM_PATCH)      ; s = 'patch'
      case (DIM_SOIL)       ; s = 'soil'
      case (DIM_PFT)        ; s = 'pft'
      case (DIM_SIZE)       ; s = 'dbh_class'
      case (DIM_SOIL_PATCH) ; s = 'soil,patch'
      case (DIM_SITE)
         s = 'by site'
      case default          ; s = '?'
      end select
   end function dim_name

end module meds_output_registry
