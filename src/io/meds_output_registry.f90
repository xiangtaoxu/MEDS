!==========================================================================================!
! meds_output_registry -- the ONE source-defined registration list of diagnostic variables, and !
! the config resolution that turns it into per-tier live-variable index lists (§3.2, §3.4, §6.1). !
!                                                                                          !
! netCDF-FREE (links types + config + demography state + the integrate SRC_* codes). Adding a      !
! variable is TWO coupled edits: one add_variable() line here + one case in extract_variable        !
! (meds_output_integrate). Resolution order (§6.1): source defaults -> group toggles -> per-tier      !
! enables -> per-variable overrides (meds_io_config.toml, applied separately) -> annual guard +        !
! freq index. Design: docs/dev_plans/MEDS_IO_DESIGN.md §3, §6.                                               !
!==========================================================================================!
module meds_output_registry
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_column_state_types, only : n_soil_layer_max, soil_params_t, curve_a, curve_n
   use meds_core_state_types,   only : site_t
   use meds_core_diag_types,    only : N_CDIAG, N_PDIAG, N_CSDIAG, cohort_diag_alloc,           &
                                       patch_diag_alloc
   use meds_output_config,  only : output_config_t, N_FREQ, N_GRP, N_AXIS,                       &
                                   N_DBH_CLASS_DEFAULT, DBH_EDGES_DEFAULT,                       &
                                   FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE,   &
                                   GRP_STRUCTURE, GRP_CARBON, GRP_WATER, GRP_ENERGY,              &
                                   GRP_RADIATION, GRP_ECOPHYS, GRP_BIOGEOCHEM, GRP_NUMERICS
   use meds_output_types,   only : var_desc_t, output_registry_t, output_manager_t,              &
                                   MAX_OUTPUT_VARS, MAX_DBH_CLASS,                                &
                                   AGG_MEAN, AGG_LAST, AGG_TMEAN, AGG_SUM, DIM_SCALAR, DIM_COHORT,&
                                   DIM_PATCH, DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH,        &
                                   XTYPE_DOUBLE, XTYPE_INT
   use meds_diagnostic_reduce, only : W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB,        &
                                      cm2_to_m2
   use meds_output_integrate, only : alloc_integ_buffer,                                         &
        FLD_C_NPLANT, FLD_C_DBH, FLD_C_HEIGHT, FLD_C_BASAL_AREA, FLD_C_AGB, FLD_C_LEAF_AREA,     &
        FLD_C_GROWTH_AVG, FLD_C_PFT, FLD_C_OWNER_PATCH, FLD_C_GLOBAL_ID, FLD_C_LAI,              &
        FLD_C_LEAF_CARBON, FLD_C_FINEROOT_CARBON, FLD_C_WOOD_CARBON, FLD_C_STORAGE_CARBON,       &
        FLD_C_BGB, FLD_C_VEG_CARBON, FLD_C_SLA, FLD_C_VCMAX25, FLD_C_RD25, FLD_C_LLSPAN,         &
        FLD_C_OVERTOP_LAI, FLD_C_GPP_ACCUM, FLD_C_NPP_ACCUM, FLD_C_LEAF_RESP, FLD_C_STEM_RESP,   &
        FLD_C_ROOT_RESP, FLD_C_DMAX_PSI_LEAF, FLD_C_PHENO_FLUSH, FLD_C_PHENO_SHED,               &
        FLD_C_LEAF_TEMP, FLD_C_WOOD_TEMP, FLD_C_SDIAG0,                                          &
        FLD_P_AREA, FLD_P_AGE, FLD_P_DIST_TYPE, FLD_P_COHORT_OFFSET, FLD_P_COHORT_COUNT,         &
        FLD_P_GLOBAL_ID, FLD_P_CAS_TEMP, FLD_P_CAS_SHV, FLD_P_CAS_CO2, FLD_P_CAS_VPD,            &
        FLD_P_CAS_DEPTH, FLD_P_SOIL_TEMP_TOP, FLD_P_SWE, FLD_P_SNOW_DEPTH, FLD_P_W_SURFACE,      &
        FLD_P_SOILC_FAST_GRND, FLD_P_SOILC_FAST_SOIL, FLD_P_SOILC_STRUCT_GRND,                   &
        FLD_P_SOILC_STRUCT_SOIL, FLD_P_SOILC_MICROBIAL, FLD_P_SOILC_SLOW,                        &
        FLD_P_SOILC_PASSIVE, FLD_P_SOILC_TOTAL, FLD_P_RH,                                        &
        FLD_L_SOIL_TEMP, FLD_L_SOIL_WATER, FLD_L_SOIL_PSI, FLD_L_SOIL_WETNESS, FLD_L_SOIL_FLIQ,  &
        SRC_S_ET, SRC_S_N_COHORT, SRC_S_N_PATCH, SRC_S_CANOPY_HEIGHT,                            &
        SRC_S_WORK_STEPS, SRC_S_WORK_REJ, SRC_S_WORK_SOIL_NSUB, SRC_S_WORK_HYDRO_NSUB,           &
        SRC_S_WORK_NONCONV, SRC_S_WORK_HYDRO_THRASH, SRC_S_WORK_RK45_RESCUE,                     &
        SRC_S_WORK_CLAMP_STAGE, SRC_S_WORK_CLAMP_COMMIT, SRC_S_WORK_CLAMP_MASS,                  &
        SRC_S_WORK_CLAMP_ENERGY,                                                                 &
        SRC_F_CAS_TEMP, SRC_F_SOIL_TEMP_TOP, SRC_F_GPP_RATE, SRC_F_LE, SRC_F_H, SRC_F_RNET,      &
        SRC_F_SW_IN, SRC_F_USTAR, SRC_F_AIR_TEMP, SRC_F_SOIL_TEMP, SRC_F_SOIL_WATER,             &
        SRC_F_COH_LEAF_TEMP, SRC_F_COH_GPP, SRC_F_COH_HEIGHT, FLD_C_DIAG0, FLD_P_DIAG0
   use meds_core_diag_types, only : CD_ANET, CD_AGROSS, CD_GSW, CD_GBW, CD_CI, CD_CS, CD_RD,     &
                                    CD_TRANSP, CD_BETA_STOM, CD_BETA_NONSTOM, CD_LEAF_VPD,       &
                                    CD_PSI_LEAF, CD_PSI_WOOD, CD_PLC, CD_SAPFLOW,                &
                                    CD_ROOT_UPTAKE, CD_ABS_PAR, CD_ABS_SW, CD_ABS_LW, CD_WIND,   &
                                    CD_LEAF_WATER, CD_WOOD_WATER,                                &
                                    PD_LE, PD_H, PD_RNET, PD_SW_IN, PD_SW_GROUND, PD_LW_GROUND,  &
                                    PD_USTAR, PD_GGNET, PD_ROUGH, PD_DISPLACE, PD_GPP, PD_NEE,   &
                                    PD_TRANSP, PD_PRECIP, PD_GROUND_TEMP, PD_RESID_ENERGY,       &
                                    PD_RESID_WATER,                                              &
                                    CS_DDBH_DT, CS_DAGB_DT, CS_MORT_RATE, CS_NPP_LEAF,            &
                                    CS_NPP_FINEROOT, CS_NPP_WOOD, CS_NPP_STORAGE, CS_NPP_REPRO,   &
                                    CS_GROWTH_RESP
   implicit none
   private

   public :: build_output_registry, build_freq_index, find_var_index, parse_stream_mask
   public :: apply_group_toggles, apply_freq_enables, apply_variable_override, apply_axis_toggles
   public :: freq_bit, dim_axis_index, OVR_TRUE, OVR_FALSE, OVR_MASK
   public :: manager_alloc, manager_setup, manager_alloc_buffers, manager_set_soil_params
   public :: activate_site_diag

   !----- Named default stream masks (readable `ior` combinations). DAY_MON_YR deliberately     !
   !      EXCLUDES the fast bit -- the fast tier is always opt-in.                                !
   integer(ik), parameter :: DAY_MON    = ior(FREQ_DAILY, FREQ_MONTHLY)
   integer(ik), parameter :: MON        = FREQ_MONTHLY
   integer(ik), parameter :: MON_YR     = ior(FREQ_MONTHLY, FREQ_ANNUAL)
   integer(ik), parameter :: DAY_MON_YR = ior(ior(FREQ_DAILY, FREQ_MONTHLY), FREQ_ANNUAL)
   integer(ik), parameter :: FAST_ONLY  = FREQ_FAST

   !----- Per-variable override kinds (meds_io_config.toml value grammar, §6.1). --------------!
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

      call enforce_annual_guard(reg)              ! cohort/patch var MUST NOT declare FREQ_ANNUAL (§3.1)
      call apply_axis_toggles(reg, cfg%output)    ! whole-axis suppression (§6.1 step 2a)
      call apply_group_toggles(reg, cfg%output)   ! high-level flux-group toggles (§6.1 step 2)
      call apply_freq_enables(reg, cfg%output)    ! per-tier enable (§6.1 step 3)
      call build_freq_index(reg)                  ! precompute idx_freq (per-variable overrides applied later)
   end subroutine build_output_registry

   !=======================================================================================!
   !  THE REGISTRATION LIST, grouped one subroutine per (group, axis) so the source order     !
   !  matches the documentation order (MEDS_IO_V01_PLAN.md section 4) and a reviewer can diff   !
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
      call add_variable(reg, 'llspan_cohort', 'leaf lifespan', 'yr',                             &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_LLSPAN)
      call add_variable(reg, 'pheno_flush_cohort', 'phenology flush governor', '-',              &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_PHENO_FLUSH)
      call add_variable(reg, 'pheno_shed_cohort', 'phenology active-shed governor', '-',         &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_PHENO_SHED)
      call add_variable(reg, 'dmax_psi_leaf_cohort', 'predawn (daily-max) leaf water potential', 'MPa', &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, FLD_C_DMAX_PSI_LEAF)
      !----- Slow-loop RATES off the per-step deriv bundle (P0: no new state needed). ---------!
      call add_variable(reg, 'dbh_growth_cohort', 'diameter growth rate', 'cm/yr',               &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_SDIAG0 + CS_DDBH_DT)
      call add_variable(reg, 'agb_growth_cohort', 'per-plant AGB growth rate', 'kgC/plant/yr',   &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_SDIAG0 + CS_DAGB_DT)
      call add_variable(reg, 'mort_rate_cohort', 'mortality hazard rate', '1/yr',                &
                        DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, MON, FLD_C_SDIAG0 + CS_MORT_RATE)
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
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'mort_rate_size', 'nplant-weighted mortality rate by DBH class', '1/yr', &
                        DIM_SIZE, AGG_MEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_MORT_RATE,             &
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
      !----- Area-weighted SITE soil columns (DIM_SOIL). --------------------------------------!
      call add_variable(reg, 'soil_water_site', 'area-weighted volumetric soil moisture', 'm3/m3', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_WATER)
      call add_variable(reg, 'soil_psi_site', 'area-weighted soil matric potential', 'MPa',      &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_PSI)
      call add_variable(reg, 'soil_wetness_site', 'area-weighted relative saturation', '-',      &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_L_SOIL_WETNESS)
      !----- Surface water stores. ------------------------------------------------------------!
      call add_variable(reg, 'w_surface_site', 'ponded surface water', 'kg/m2',                  &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_W_SURFACE)
      call add_variable(reg, 'swe_site', 'snow water equivalent', 'kg/m2',                       &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_SWE)
      call add_variable(reg, 'snow_depth_site', 'snow depth', 'm',                               &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_SNOW_DEPTH)
      call add_variable(reg, 'w_surface_patch', 'patch ponded surface water', 'kg/m2',           &
                        DIM_PATCH, AGG_TMEAN, GRP_WATER, MON, FLD_P_W_SURFACE)
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
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_CAS_TEMP)
      call add_variable(reg, 'cas_shv_site', 'site canopy-air specific humidity', 'kg/kg',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_CAS_SHV)
      call add_variable(reg, 'cas_co2_site', 'site canopy-air CO2 mixing ratio', 'umol/mol',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_CAS_CO2)
      call add_variable(reg, 'cas_vpd_site', 'site canopy-air vapour-pressure deficit', 'Pa',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_CAS_VPD)
      call add_variable(reg, 'cas_depth_site', 'site canopy-air-space depth', 'm',               &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_CAS_DEPTH)
      call add_variable(reg, 'soil_temp_top_site', 'site soil-top (layer 1) temperature', 'K',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_SOIL_TEMP_TOP)
      call add_variable(reg, 'soil_temp_site', 'area-weighted soil temperature', 'K',            &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_L_SOIL_TEMP)
      call add_variable(reg, 'soil_fliq_site', 'area-weighted soil liquid fraction', '-',        &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_L_SOIL_FLIQ)
      !----- Canopy temperatures, LEAF-AREA-weighted (the intensive rule: a bare sapling must    !
      !      not pull the canopy mean as hard as a closed overstory).  --------------------------!
      call add_variable(reg, 'leaf_temp_site', 'leaf-area-weighted canopy leaf temperature', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_C_LEAF_TEMP,         &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'wood_temp_site', 'leaf-area-weighted wood temperature', 'K',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_C_WOOD_TEMP,         &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'leaf_temp_cohort', 'per-cohort leaf temperature', 'K',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_C_LEAF_TEMP)
      call add_variable(reg, 'wood_temp_cohort', 'per-cohort wood temperature', 'K',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_C_WOOD_TEMP)
      call add_variable(reg, 'cas_temp_patch', 'patch canopy-air-space temperature', 'K',        &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_CAS_TEMP)
      call add_variable(reg, 'soil_temp_top_patch', 'patch soil-top temperature', 'K',           &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, MON, FLD_P_SOIL_TEMP_TOP)
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
      call add_variable(reg, 'work_soil_nsub_site', 'soil-water solver sub-steps (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_SOIL_NSUB)
      call add_variable(reg, 'work_hydro_nsub_site', 'plant-hydraulics sub-steps (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_HYDRO_NSUB)
      call add_variable(reg, 'work_nonconv_site', 'solver non-convergence events (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_NONCONV)
      call add_variable(reg, 'work_hydro_thrash_site',                                          &
                        'dt_fast steps with pathological hydraulics sub-stepping (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_HYDRO_THRASH)
      call add_variable(reg, 'work_rk45_rescue_site', 'RK45 steps rescued to the implicit path (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_RK45_RESCUE)
      call add_variable(reg, 'work_clamp_stage_site', 'stage-input clamp activations (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_CLAMP_STAGE)
      call add_variable(reg, 'work_clamp_commit_site', 'committed-state clamp activations (period total)', &
                        '--', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_CLAMP_COMMIT)
      call add_variable(reg, 'work_clamp_mass_site', 'water moved by committed-state clamps (period total)', &
                        'kg/m2', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_CLAMP_MASS)
      call add_variable(reg, 'work_clamp_energy_site', 'energy moved by committed-state clamps (period total)', &
                        'J/m2', DIM_SCALAR, AGG_SUM, GRP_NUMERICS, DAY_MON_YR, SRC_S_WORK_CLAMP_ENERGY)
   end subroutine register_numerics


   !=======================================================================================!
   !  GRP_ECOPHYS -- the per-cohort leaf gas-exchange + hydraulics set (MEDS_IO_V01_PLAN.md      !
   !  section 4.7). Every row reads the fast-loop diagnostic block, so every row is a quantity     !
   !  the model already computed ~48 times a day and previously threw away.                        !
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
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_ANET)
      call add_variable(reg, 'agross_cohort', 'gross leaf assimilation', 'umol/m2 leaf/s',       &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_AGROSS)
      call add_variable(reg, 'gsw_cohort', 'stomatal conductance to water vapour', 'mol/m2 leaf/s', &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_GSW)
      call add_variable(reg, 'gbw_cohort', 'leaf boundary-layer conductance', 'm/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_GBW)
      call add_variable(reg, 'ci_cohort', 'intercellular CO2 mole fraction', 'umol/mol',         &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_CI)
      call add_variable(reg, 'cs_cohort', 'leaf-surface CO2 mole fraction', 'umol/mol',          &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_CS)
      call add_variable(reg, 'rd_cohort', 'leaf dark respiration', 'umol/m2 leaf/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_RD)
      call add_variable(reg, 'transp_cohort', 'leaf transpiration', 'mol/m2 leaf/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_TRANSP)
      call add_variable(reg, 'beta_stomata_cohort', 'stomatal water-stress factor', '-',         &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_BETA_STOM)
      call add_variable(reg, 'beta_nonstomata_cohort', 'non-stomatal (capacity) water-stress factor', '-', &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_BETA_NONSTOM)
      call add_variable(reg, 'leaf_vpd_cohort', 'leaf-to-air vapour-pressure deficit', 'Pa',     &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_LEAF_VPD)
      call add_variable(reg, 'psi_leaf_cohort', 'leaf water potential', 'MPa',                   &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_PSI_LEAF)
      call add_variable(reg, 'psi_wood_cohort', 'wood water potential', 'MPa',                   &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_PSI_WOOD)
      call add_variable(reg, 'plc_cohort', 'plant loss of hydraulic conductance', '-',           &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_PLC)
      call add_variable(reg, 'sapflow_cohort', 'wood-to-leaf sapflow', 'kg/plant/s',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_SAPFLOW)
      call add_variable(reg, 'root_uptake_cohort', 'soil-to-root water uptake', 'kg/plant/s',    &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, DAY_MON, FLD_C_DIAG0 + CD_ROOT_UPTAKE)
      call add_variable(reg, 'leaf_water_cohort', 'internal leaf water', 'kg/plant',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_LEAF_WATER)
      call add_variable(reg, 'wood_water_cohort', 'internal wood water', 'kg/plant',             &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_WOOD_WATER)
      call add_variable(reg, 'abs_par_cohort', 'absorbed photosynthetically active radiation', 'W/m2', &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_C_DIAG0 + CD_ABS_PAR)
      call add_variable(reg, 'abs_sw_cohort', 'absorbed shortwave radiation', 'W/m2',            &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_C_DIAG0 + CD_ABS_SW)
      call add_variable(reg, 'abs_lw_cohort', 'net longwave radiation', 'W/m2',                  &
                        DIM_COHORT, AGG_TMEAN, GRP_RADIATION, MON, FLD_C_DIAG0 + CD_ABS_LW)
      call add_variable(reg, 'wind_cohort', 'in-canopy wind speed', 'm/s',                       &
                        DIM_COHORT, AGG_TMEAN, GRP_ECOPHYS, MON, FLD_C_DIAG0 + CD_WIND)
      !----- CANOPY-LEVEL twins: leaf-area-weighted means, the standard flux-tower comparands.  --!
      call add_variable(reg, 'anet_site', 'leaf-area-weighted canopy net assimilation', 'umol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CD_ANET,   &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'gsw_site', 'leaf-area-weighted canopy stomatal conductance', 'mol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CD_GSW,    &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'psi_leaf_site', 'leaf-area-weighted canopy leaf water potential', 'MPa', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CD_PSI_LEAF, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'plc_site', 'leaf-area-weighted canopy loss of conductance', '-',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CD_PLC,    &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'beta_stomata_site', 'leaf-area-weighted stomatal water-stress factor', '-', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ECOPHYS, DAY_MON_YR, FLD_C_DIAG0 + CD_BETA_STOM, &
                        w=W_LEAF_AREA, mn=.true.)
      call add_variable(reg, 'root_uptake_site', 'site root water uptake', 'kg/m2/s',            &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CD_ROOT_UPTAKE, &
                        w=W_NPLANT)
      call add_variable(reg, 'transp_leaf_site', 'leaf-area-weighted canopy transpiration', 'mol/m2 leaf/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_C_DIAG0 + CD_TRANSP,   &
                        w=W_LEAF_AREA, mn=.true.)
   end subroutine register_ecophys

   !=======================================================================================!
   !  SLOW-loop carbon allocation. These were pure locals inside compute_carbon_allocation --  !
   !  computed for every cohort on every step and discarded -- so the carbon budget could not   !
   !  be closed from an output file at all. With them, npp_site should reconcile against the     !
   !  sum of its five destinations plus growth respiration.                                      !
   !=======================================================================================!
   subroutine register_allocation(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'npp_leaf_site', 'NPP allocated to leaf', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_NPP_LEAF, w=W_NPLANT)
      call add_variable(reg, 'npp_fineroot_site', 'NPP allocated to fine root', 'kgC/m2/yr',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_NPP_FINEROOT, w=W_NPLANT)
      call add_variable(reg, 'npp_wood_site', 'NPP allocated to wood', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_NPP_WOOD, w=W_NPLANT)
      call add_variable(reg, 'npp_storage_site', 'NPP to non-structural storage', 'kgC/m2/yr',   &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_NPP_STORAGE, w=W_NPLANT)
      call add_variable(reg, 'npp_repro_site', 'NPP allocated to reproduction', 'kgC/m2/yr',     &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_NPP_REPRO, w=W_NPLANT)
      call add_variable(reg, 'growth_resp_site', 'growth respiration', 'kgC/m2/yr',              &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, MON_YR, FLD_C_SDIAG0 + CS_GROWTH_RESP, w=W_NPLANT)
      !----- Site-level demographic RATES: the three terms that MAKE the AGB trajectory, so a      !
      !      reader can check d(agb)/dt against growth - mortality directly from the file.  -------!
      call add_variable(reg, 'agb_growth_site', 'site AGB growth rate', 'kgC/m2/yr',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'agb_mort_site', 'site AGB mortality loss rate', 'kgC/m2/yr',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_MORT_RATE, w=W_AGB)
      call add_variable(reg, 'mort_rate_site', 'nplant-weighted mortality rate', '1/yr',         &
                        DIM_SCALAR, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_MORT_RATE, &
                        w=W_NPLANT, mn=.true.)
      call add_variable(reg, 'agb_growth_pft', 'AGB growth rate by PFT', 'kgC/m2/yr',            &
                        DIM_PFT, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_DAGB_DT, w=W_NPLANT)
      call add_variable(reg, 'agb_mort_pft', 'AGB mortality loss rate by PFT', 'kgC/m2/yr',      &
                        DIM_PFT, AGG_TMEAN, GRP_STRUCTURE, MON_YR, FLD_C_SDIAG0 + CS_MORT_RATE, w=W_AGB)
   end subroutine register_allocation

   !=======================================================================================!
   !  The per-PATCH fast-loop fluxes and states, and their area-weighted site twins. These are  !
   !  the sub-daily energy / water / carbon exchange terms -- the ones a run is judged on -- now  !
   !  available at EVERY tier rather than only in the opt-in FAST stream.                         !
   !=======================================================================================!
   subroutine register_patch_fluxes(reg)
      type(output_registry_t), intent(inout) :: reg
      !--- energy ---!
      call add_variable(reg, 'le_site', 'latent heat flux', 'W/m2',                              &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_LE)
      call add_variable(reg, 'h_site', 'sensible heat flux', 'W/m2',                             &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_H)
      call add_variable(reg, 'rnet_site', 'net all-wave radiation', 'W/m2',                      &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_RNET)
      call add_variable(reg, 'sw_in_site', 'incident shortwave at canopy top', 'W/m2',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_SW_IN)
      call add_variable(reg, 'sw_ground_site', 'shortwave absorbed at the ground', 'W/m2',       &
                        DIM_SCALAR, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_P_DIAG0 + PD_SW_GROUND)
      call add_variable(reg, 'lw_ground_site', 'net longwave at the ground', 'W/m2',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_RADIATION, DAY_MON, FLD_P_DIAG0 + PD_LW_GROUND)
      call add_variable(reg, 'ustar_site', 'friction velocity', 'm/s',                           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_USTAR)
      call add_variable(reg, 'ggnet_site', 'ground conductance', 'm/s',                          &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON, FLD_P_DIAG0 + PD_GGNET)
      call add_variable(reg, 'rough_site', 'aerodynamic roughness length', 'm',                  &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_DIAG0 + PD_ROUGH)
      call add_variable(reg, 'displace_site', 'zero-plane displacement height', 'm',             &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, MON_YR, FLD_P_DIAG0 + PD_DISPLACE)
      call add_variable(reg, 'ground_temp_site', 'ground (skin) temperature', 'K',               &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, DAY_MON_YR, FLD_P_DIAG0 + PD_GROUND_TEMP)
      !--- carbon ---!
      call add_variable(reg, 'gpp_rate_site', 'gross primary productivity (mean rate)', 'umol/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, DAY_MON_YR, FLD_P_DIAG0 + PD_GPP)
      call add_variable(reg, 'nee_site', 'net ecosystem exchange (mean rate)', 'umol/m2/s',      &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, DAY_MON_YR, FLD_P_DIAG0 + PD_NEE)
      !--- water ---!
      call add_variable(reg, 'et_rate_site', 'evapotranspiration (mean rate)', 'kg/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PD_TRANSP)
      call add_variable(reg, 'precip_site', 'total precipitation (mean rate)', 'kg/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, DAY_MON_YR, FLD_P_DIAG0 + PD_PRECIP)
      !--- budget health (GRP_NUMERICS): the numbers that say whether anything above is real. ---!
      call add_variable(reg, 'resid_energy_site', 'worst whole-column energy-budget residual', 'W/m2', &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PD_RESID_ENERGY)
      call add_variable(reg, 'resid_water_site', 'worst whole-column water-budget residual', 'kg/m2/s', &
                        DIM_SCALAR, AGG_TMEAN, GRP_NUMERICS, DAY_MON_YR, FLD_P_DIAG0 + PD_RESID_WATER)
      !--- per-patch twins of the four that vary most between a gap and a closed canopy. ---!
      call add_variable(reg, 'le_patch', 'patch latent heat flux', 'W/m2',                       &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PD_LE)
      call add_variable(reg, 'h_patch', 'patch sensible heat flux', 'W/m2',                      &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PD_H)
      call add_variable(reg, 'rnet_patch', 'patch net all-wave radiation', 'W/m2',               &
                        DIM_PATCH, AGG_TMEAN, GRP_ENERGY, DAY_MON, FLD_P_DIAG0 + PD_RNET)
      call add_variable(reg, 'nee_patch', 'patch net ecosystem exchange', 'umol/m2/s',           &
                        DIM_PATCH, AGG_TMEAN, GRP_CARBON, DAY_MON, FLD_P_DIAG0 + PD_NEE)
   end subroutine register_patch_fluxes

   !=======================================================================================!
   !  FAST-tier (sub-daily) diagnostics. These read the STAGED per-sub-step samples, not live   !
   !  site state, so they carry the FAST bit ONLY -- putting one on a coarse tier would read a   !
   !  staging buffer that is not populated there. (P2 removes this restriction by moving the      !
   !  sub-step capture into site state.)                                                          !
   !=======================================================================================!
   subroutine register_fast(reg)
      type(output_registry_t), intent(inout) :: reg
      call add_variable(reg, 'cas_temp_fast', 'site canopy-air-space temperature (sub-daily)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_CAS_TEMP)
      call add_variable(reg, 'soil_temp_top_fast', 'site soil-top temperature (sub-daily)', 'K', &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_SOIL_TEMP_TOP)
      call add_variable(reg, 'gpp_rate_fast', 'site instantaneous GPP rate', 'umol/m2/s',        &
                        DIM_SCALAR, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_GPP_RATE)
      call add_variable(reg, 'le_flux_fast', 'site latent-heat (ET) flux', 'W/m2',               &
                        DIM_SCALAR, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_LE)
      call add_variable(reg, 'h_flux_fast', 'site sensible-heat flux', 'W/m2',                   &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_H)
      call add_variable(reg, 'rnet_fast', 'site net all-wave radiation', 'W/m2',                 &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_RNET)
      call add_variable(reg, 'sw_in_fast', 'incident shortwave at canopy top', 'W/m2',           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_SW_IN)
      call add_variable(reg, 'ustar_fast', 'friction velocity', 'm/s',                           &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_USTAR)
      call add_variable(reg, 'air_temp_fast', 'reference-level forcing air temperature', 'K',    &
                        DIM_SCALAR, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_AIR_TEMP)
      call add_variable(reg, 'soil_temp_site_fast', 'area-weighted soil temperature (sub-daily)', 'K', &
                        DIM_SOIL, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_SOIL_TEMP)
      call add_variable(reg, 'soil_water_site_fast', 'area-weighted soil moisture (sub-daily)', 'm3/m3', &
                        DIM_SOIL, AGG_TMEAN, GRP_WATER, FAST_ONLY, SRC_F_SOIL_WATER)
      call add_variable(reg, 'leaf_temp_cohort_fast', 'per-cohort leaf temperature (sub-daily)', 'K', &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_COH_LEAF_TEMP)
      call add_variable(reg, 'gpp_cohort_fast', 'per-cohort GPP rate (sub-daily)', 'umol/plant/s', &
                        DIM_COHORT, AGG_TMEAN, GRP_CARBON, FAST_ONLY, SRC_F_COH_GPP)
      call add_variable(reg, 'height_cohort_fast', 'per-cohort height (tallest-cohort selection)', 'm', &
                        DIM_COHORT, AGG_TMEAN, GRP_ENERGY, FAST_ONLY, SRC_F_COH_HEIGHT)
   end subroutine register_fast

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

   !----- Step 2a: suppress a whole AXIS without naming variables. With ~55 cohort-dimensioned  !
   !      variables, "site level only, nothing per-cohort" is the most common request and was      !
   !      previously inexpressible; it is also the single biggest lever on output volume.          !
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

   !----- Guard: a cohort/patch-dimensioned variable may NOT be on the annual stream (§3.1). --!
   subroutine enforce_annual_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         if ((reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) .and.               &
             iand(reg%var(k)%streams_default, FREQ_ANNUAL) /= 0_ik)                              &
            error stop 'meds_output_registry: cohort/patch variable on the annual stream ('//trim(reg%var(k)%name)//')'
      end do
   end subroutine enforce_annual_guard

   !----- Step 2: disable every variable in a group whose main-config toggle is .false. -------!
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

   !----- Step 3: a disabled tier clears its FREQ_* bit from every variable's mask. -----------!
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

   !----- Step 4 primitive: apply one per-variable override (from meds_io_config.toml, §6.1).  !
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
         if ((reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) .and.               &
             iand(m, FREQ_ANNUAL) /= 0_ik)                                                       &
            error stop 'meds_output_registry: annual (Y) stream forbidden on cohort/patch variable ('//trim(name)//')'
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

   !----- Parse a "F D M Y" stream-string into an ior(FREQ_*) mask (§6.1 value grammar). ------!
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
   !  Build the manager's REGISTRY + config (netCDF-free; NO buffers yet). Split from buffer    !
   !  allocation so a caller can apply meds_io_config.toml per-variable overrides to mgr%reg      !
   !  (which change which (var,tier) buffers are needed) BEFORE manager_alloc_buffers (§6.1).     !
   !=======================================================================================!
   subroutine manager_setup(mgr, cfg)
      type(output_manager_t), intent(out) :: mgr
      type(meds_config_t),    intent(in)  :: cfg
      integer(ik) :: ne
      mgr%enabled = cfg%output%enabled
      call build_output_registry(mgr%reg, cfg)
      mgr%cohort_max = cfg%output%cohort_max
      mgr%patch_max  = cfg%output%patch_max
      mgr%dir        = cfg%output%dir
      mgr%prefix     = cfg%output%prefix
      mgr%file_chunk = cfg%output%file_chunk
      mgr%sync_every = cfg%output%sync_every
      mgr%fast_interval_steps = cfg%output%fast_interval_steps

      !----- Run-dependent parameters the DERIVED diagnostics need (diag_params_t). n_pft is the  !
      !      RUN-TIME PFT count, so the netCDF `pft` dimension is run-dependent -- which is why      !
      !      the serializer also writes a `pft` coordinate variable, so a file stays self-describing !
      !      when compared across runs with different PFT tables.                                    !
      mgr%diag%n_pft  = cfg%pft%n
      mgr%diag%n_soil = n_soil_layer_max
      !----- DBH size classes: TOML edges if given, else the built-in inventory set. -----------!
      if (cfg%output%n_dbh_class > 0_ik) then
         mgr%diag%n_dbh_class = min(cfg%output%n_dbh_class, MAX_DBH_CLASS)
         ne = mgr%diag%n_dbh_class + 1_ik
         mgr%diag%dbh_edges(1:ne) = cfg%output%dbh_edges(1:ne)
      else
         mgr%diag%n_dbh_class = N_DBH_CLASS_DEFAULT
         mgr%diag%dbh_edges(1:N_DBH_CLASS_DEFAULT + 1_ik) = DBH_EDGES_DEFAULT
      end if
      call check_dbh_edges(mgr%diag%dbh_edges, mgr%diag%n_dbh_class)
      !----- The soil retention parameters are NOT set here: meds_main copies them from the FAST  !
      !      CONTEXT the physics actually ran (manager_set_soil_params), so a reported psi and the  !
      !      psi the roots saw are the same curve by construction. Until that call, soil_ready is   !
      !      .false. and the psi/wetness diagnostics emit _FillValue rather than a plausible        !
      !      wrong number from an assumed texture.                                                  !
      mgr%diag%soil_ready = .false.

      !----- SLAB SIZING. Size the shared pending-record slab to the largest axis that is        !
      !      ACTUALLY LIVE, not to the largest axis that exists. With ~55 cohort-dimensioned        !
      !      variables the pending records are the dominant memory term, and a run with            !
      !      axes_cohort = false (or the 2-D soil axis off, the default) should not pay for the     !
      !      axis it switched off. Computed AFTER the registry is finalized, for exactly that       !
      !      reason.                                                                                !
      mgr%max_slab = live_max_slab(mgr)
   end subroutine manager_setup

   !----- Largest slab length any LIVE variable can produce (see the sizing note above). ------!
   pure integer(ik) function live_max_slab(mgr) result(cap)
      type(output_manager_t), intent(in) :: mgr
      integer(ik) :: k
      cap = 1_ik
      do k = 1_ik, mgr%reg%nvar
         if (.not. mgr%reg%var(k)%enabled)      cycle
         if (mgr%reg%var(k)%streams == FREQ_NONE) cycle
         cap = max(cap, dim_capacity(mgr, mgr%reg%var(k)%dim))
      end do
   end function live_max_slab

   !----- Slab capacity of one axis. ---------------------------------------------------------!
   pure integer(ik) function dim_capacity(mgr, dm) result(cap)
      type(output_manager_t), intent(in) :: mgr
      integer(ik),            intent(in) :: dm
      select case (dm)
      case (DIM_COHORT)     ; cap = mgr%cohort_max
      case (DIM_PATCH)      ; cap = mgr%patch_max
      case (DIM_SOIL)       ; cap = n_soil_layer_max
      case (DIM_PFT)        ; cap = max(mgr%diag%n_pft, 1_ik)
      case (DIM_SIZE)       ; cap = max(mgr%diag%n_dbh_class, 1_ik)
      case (DIM_SOIL_PATCH) ; cap = mgr%patch_max * n_soil_layer_max
      case default          ; cap = 0_ik
      end select
   end function dim_capacity

   !----- Install the soil retention parameters the psi/wetness diagnostics need. Called by     !
   !      meds_main with the SAME soil_params_t the fast loop integrates on.                     !
   subroutine manager_set_soil_params(mgr, params)
      type(output_manager_t), intent(inout) :: mgr
      type(soil_params_t),    intent(in)    :: params
      integer(ik) :: k
      mgr%diag%retention = params%retention
      do k = 1_ik, n_soil_layer_max
         mgr%diag%theta_sat(k) = params%theta_sat(k)
         mgr%diag%theta_res(k) = params%theta_res(k)
         !----- The generic (a, n) pair means (alpha, n) for van Genuchten and (psi_sat, b) for  !
         !      Campbell. Resolved through the SAME accessors the Richards solver uses, so the    !
         !      diagnostic psi cannot come from a different curve than the one integrated.  ------!
         mgr%diag%par_a(k)     = curve_a(params, k)
         mgr%diag%par_n(k)     = curve_n(params, k)
      end do
      mgr%diag%n_soil     = params%n_active
      mgr%diag%soil_ready = .true.
   end subroutine manager_set_soil_params

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

   !----- Allocate the integrator buffers + pending records + stream handles from the (now        !
   !      FINALIZED) registry. Call after manager_setup [+ overrides].                             !
   subroutine manager_alloc_buffers(mgr)
      type(output_manager_t), intent(inout) :: mgr
      integer(ik) :: t, k, nv, cap, bit
      nv = mgr%reg%nvar
      allocate(mgr%buf(nv, N_FREQ))
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, nv
            if (mgr%reg%var(k)%enabled .and. iand(mgr%reg%var(k)%streams, bit) /= 0_ik) then
               cap = dim_capacity(mgr, mgr%reg%var(k)%dim)
               call alloc_integ_buffer(mgr%buf(k,t), mgr%reg%var(k), bit, cap)
               mgr%buf(k,t)%var_id = k
            end if
         end do
      end do
      do t = 1_ik, N_FREQ
         allocate(mgr%pending(t)%sval(nv), mgr%pending(t)%svalid(nv), mgr%pending(t)%nslab(nv))
         allocate(mgr%pending(t)%slab(mgr%max_slab, nv), mgr%pending(t)%slabvalid(mgr%max_slab, nv))
         mgr%pending(t)%used = .false.
         mgr%stream(t)%freq  = freq_bit(t)
         allocate(mgr%stream(t)%vid(nv)) ; mgr%stream(t)%vid = -1_ik
      end do
   end subroutine manager_alloc_buffers

   !----- Convenience: registry + buffers in one call (no per-variable overrides). ------------!
   subroutine manager_alloc(mgr, cfg)
      type(output_manager_t), intent(out) :: mgr
      type(meds_config_t),    intent(in)  :: cfg
      call manager_setup(mgr, cfg)
      call manager_alloc_buffers(mgr)
   end subroutine manager_alloc

   !----- Step 5: precompute, per tier, the list of live (enabled + in-tier) variable indices.-!
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
   !  point in meds_core_diag_types a no-op, and the fast loop never even asks the leaf kernel for   !
   !  the extra flux fields.                                                                         !
   !=======================================================================================!
   subroutine activate_site_diag(mgr, site)
      type(output_manager_t), intent(in)    :: mgr
      type(site_t),           intent(inout) :: site
      logical     :: need_c, need_p, need_s
      integer(ik) :: t, j, k, src
      need_c = .false. ; need_p = .false. ; need_s = .false.
      if (mgr%enabled) then
         do t = 1_ik, N_FREQ
            do j = 1_ik, mgr%reg%nidx(t)
               k   = mgr%reg%idx_freq(j, t)
               src = mgr%reg%var(k)%source_id
               if (src > FLD_C_DIAG0 .and. src <= FLD_C_DIAG0 + N_CDIAG) need_c = .true.
               if (src > FLD_P_DIAG0 .and. src <= FLD_P_DIAG0 + N_PDIAG) need_p = .true.
               if (src > FLD_C_SDIAG0 .and. src <= FLD_C_SDIAG0 + N_CSDIAG) need_s = .true.
            end do
         end do
      end if
      call cohort_diag_alloc(site%cohort%diag,  max(site%cohort%cap, 1_ik), need_c)
      call cohort_diag_alloc(site%cohort%sdiag, max(site%cohort%cap, 1_ik), need_s, nfield=N_CSDIAG)
      call patch_diag_alloc (site%patch%diag,   max(site%patch%cap,  1_ik), need_p)
   end subroutine activate_site_diag

end module meds_output_registry
