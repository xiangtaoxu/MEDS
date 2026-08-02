!==========================================================================================!
! meds_pft_params -- plant functional type (PFT) trait table.                              !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Size allometry is    !
! pan-tropical and PFT-independent (see meds_allometry); the per-PFT allometric inputs are    !
! `dbh_critical` (the maximum diameter) and `wood_density` (rho, which enters AGB).            !
!                                                                                          !
! The PFTs carry the parameters of the PHENOMENOLOGICAL vital rates                            !
! (meds_demography_rates, assembled by meds_vegetation_dynamics):                                                        !
!   * GROWTH: an intrinsic capped log-linear function of dbh (growth_dbh_slope/cap/max),       !
!     suppressed multiplicatively by neighbourhood competition (growth_lai_slope on overtopping!
!     LAI) and by reproductive allocation (reproduction_investment_fraction above maturity).   !
!   * MORTALITY: the Camac et al. (2018, PNAS) additive hazard, simplified to                  !
!     rate = mort_gamma + mort_alpha*exp(-mort_beta*growth_avg), with the three parameters     !
!     DERIVED from wood density (low rho => higher baseline & low-growth hazard).              !
!   * RECRUITMENT: a baseline external `seed_rain_recruits` plus a reproduction flux computed  !
!     in the rate module from the carbon diverted to reproduction.                             !
!                                                                                          !
! Two shared height thresholds: `min_cohort_height` (the smallest tracked cohort; recruits are !
! born here) and `min_reproduction_height` (the height a cohort must exceed to reproduce). The !
! default table seeds three strategies (pioneer/mid/climax) differing only in wood density and !
! maximum diameter.                                                                          !
!==========================================================================================!
module meds_pft_params
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pft_table_t, alloc_pft_table, derive_pft_rates, derive_leaf_params
   public :: PATH_C3, PATH_C4

   !----- Photosynthetic-pathway flags (per-PFT trait values). ----------------------------!
   integer(ik), parameter :: PATH_C3 = 1_ik   !< C3 (Farquhar-von Caemmerer-Berry)
   integer(ik), parameter :: PATH_C4 = 2_ik   !< C4 (Collatz et al. 1992)

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- Size limits + the wood-density axis (allometry itself is global). ------------!
      real(wp), allocatable :: dbh_critical(:)   !< [cm]    maximum diameter (growth clamp)
      real(wp), allocatable :: hgt_max(:)        !< [m]     asymptotic max height (per-PFT allometry cap)
      real(wp), allocatable :: wood_density(:)   !< [g/cm3] rho: AGB + mortality anchor
      !----- Reproductive allocation (carbon path). The empirical growth-curve params           !
      !       (growth_dbh_slope/cap/max, growth_lai_slope) were REMOVED: the phenomenological      !
      !       growth law is an experiment "example" and now lives in the Python example.           !
      real(wp), allocatable :: reproduction_investment_fraction(:) !< [--] growth fraction diverted to reproduction
      real(wp), allocatable :: repro_carbon_efficiency(:)          !< [--] reproduction carbon -> establishable recruits
      !----- Wood-density -> mortality-hazard DERIVATION coefficients (Camac et al. 2018 PNAS,  !
      !       a power law centred on mort_rho_ref: param = param_0 * (rho/rho_ref)^exp, always   !
      !       positive). Shared scalars for now (PFT-specific later). --------------------------!
      real(wp) :: mort_rho_ref            !< [g/cm3] wood-density reference for the power law
      real(wp) :: mort_gamma_0, mort_gamma_exp   !< baseline hazard scale + exponent
      real(wp) :: mort_alpha_0, mort_alpha_exp   !< low-growth hazard scale + exponent
      real(wp) :: mort_beta_0,  mort_beta_exp    !< growth-sensitivity scale + exponent
      !----- Wood-density-derived mortality-hazard parameters (DERIVED from the above). -------!
      real(wp), allocatable :: mort_gamma(:)     !< [1/yr]  growth-independent baseline hazard
      real(wp), allocatable :: mort_alpha(:)     !< [1/yr]  low-growth hazard magnitude
      real(wp), allocatable :: mort_beta(:)      !< [yr/cm] growth sensitivity of the hazard
      !----- Recruitment. -----------------------------------------------------------------!
      real(wp),    allocatable :: seed_rain_recruits(:) !< [plant/m2/yr] baseline external seed rain
      integer(ik), allocatable :: include_pft(:)        !< 1 = PFT may recruit, 0 = excluded
      !----- Shared height thresholds. ----------------------------------------------------!
      real(wp) :: min_cohort_height          !< [m] smallest tracked cohort; recruits born here
      real(wp) :: min_reproduction_height    !< [m] height a cohort must exceed to reproduce
      !----- Leaf photosynthesis / stomatal traits (per PFT; see meds_plant_interface). -------!
      integer(ik), allocatable :: photosynthetic_pathway(:) !< PATH_C3 | PATH_C4
      real(wp), allocatable :: vcmax25(:)            !< [umol/m2/s] max carboxylation at 25 degC
      real(wp), allocatable :: jmax_vcmax_ratio(:)   !< [--]    Jmax25 / Vcmax25
      real(wp), allocatable :: tpu_vcmax_ratio(:)    !< [--]    TPU25 / Vcmax25 (C3 product limit)
      real(wp), allocatable :: rd_vcmax_ratio(:)     !< [--]    Rd25 / Vcmax25 (leaf respiration)
      real(wp), allocatable :: kp25(:)               !< [mol/m2/s] C4 PEPcase initial slope at 25 degC
      real(wp), allocatable :: stomatal_g0(:)        !< [mol/m2/s] cuticular/residual conductance
      real(wp), allocatable :: stomatal_g1(:)        !< Leuning [--] / Medlyn [kPa^0.5] slope
      real(wp), allocatable :: stomatal_d0(:)        !< [Pa]    Leuning humidity sensitivity
      real(wp), allocatable :: quantum_yield_c4(:)   !< [mol CO2/mol photon] C4 light-limited slope
      real(wp), allocatable :: theta_j(:)            !< [--]    non-rectangular hyperbola curvature (C3 J)
      real(wp), allocatable :: theta_cj_c4(:)        !< [--]    C4 co-limitation curvature 1
      real(wp), allocatable :: theta_ic_c4(:)        !< [--]    C4 co-limitation curvature 2
      real(wp), allocatable :: katul_lambda25(:)     !< [umol CO2/mol H2O] Katul marginal water-use efficiency
      real(wp), allocatable :: wstress_psi_open(:)   !< [MPa]   leaf potential at which beta = 1 (<= 0)
      real(wp), allocatable :: wstress_psi_close(:)  !< [MPa]   leaf potential at which beta = 0 (< psi_open)
      real(wp), allocatable :: wstress_lambda_exp(:) !< [--]    Katul lambda water-stress exponent
      real(wp), allocatable :: wstress_sref_stomata(:) !< [1/MPa] Sabot stomatal-stress sensitivity (beta_stomata = exp(sref*psi), psi = predawn leaf psi)
      !----- Leaf photosynthesis DERIVED per-PFT (derive_leaf_params). ------------------------!
      real(wp), allocatable :: jmax25(:)             !< [umol/m2/s] DERIVED = jmax_vcmax_ratio * vcmax25
      real(wp), allocatable :: tpu25(:)              !< [umol/m2/s] DERIVED = tpu_vcmax_ratio  * vcmax25
      real(wp), allocatable :: rd25(:)               !< [umol/m2/s] DERIVED = rd_vcmax_ratio   * vcmax25
      !----- Carbon-dynamics per-PFT traits (meds_plant_carbon_dynamics). All state is CARBON;   !
      !       every biomass<->carbon conversion is folded into these traits here at init. Consumed !
      !       once the carbon-allocation engine is wired (a later PR); the size-target allometry    !
      !       size2leaf_carbon / size2wood_carbon uses sla / aboveground_frac.                      !
      real(wp),    allocatable :: sla(:)                    !< [m2/kgC] specific leaf area (un-folds SLA from lai_b1)
      !----- WOOD AREA INDEX allometry (ED2 b1WAI/b2WAI): wai = nplant*b1wai*dbh**b2wai. This        !
      !      REPLACES the wai = 0.20*lai placeholder, which tied wood area to LEAF area and made     !
      !      the wood thermal timescale ~6-20x too short. ------------------------------------------!
      real(wp),    allocatable :: wai_b1(:)                 !< [--] WAI intercept  (ED2 b1WAI)
      real(wp),    allocatable :: wai_b2(:)                 !< [--] WAI exponent   (ED2 b2WAI)
      !----- SAPWOOD-AREA allometry (ED2 b1SA/b2SA, dbh2sf): sapwood area = b1sa*dbh**b2sa, and the  !
      !      sapwood FRACTION of basal area sets bsap. REPLACES bsap = 0.10*wood_carbon. ------------!
      real(wp),    allocatable :: sapwood_area_b1(:)        !< [cm2/cm^b2] sapwood-area intercept (ED2 b1SA)
      real(wp),    allocatable :: sapwood_area_b2(:)        !< [--]        sapwood-area exponent  (ED2 b2SA)
      real(wp),    allocatable :: root_to_leaf_ratio(:)     !< [--]     fine-root:leaf target ratio (ED2 q)
      real(wp),    allocatable :: huber_value(:)            !< [m2 sap/m2 leaf] sapwood-area:leaf-area (sapwood + hydraulics)
      real(wp),    allocatable :: aboveground_frac(:)       !< [--]     aboveground fraction of woody carbon (ED2 agf_bs)
      real(wp),    allocatable :: storage_cushion(:)        !< [--]     storage target as a multiple of the leaf target
      real(wp),    allocatable :: growth_resp_factor(:)     !< [--]     construction cost (fraction of metabolic NPP)
      real(wp),    allocatable :: leaf_lifespan_toc(:)      !< [yr]     top-of-canopy leaf lifespan; baseline leaf
                                                            !<          turnover = 1/llspan (was leaf_turnover_rate)
      real(wp),    allocatable :: fineroot_turnover_rate(:) !< [1/yr]   baseline fine-root turnover
      real(wp),    allocatable :: wood_carbon_density(:)    !< [kgC/m3] wood carbon density (Huber sapwood carbon)
      integer(ik), allocatable :: evergreen(:)              !< 1 = evergreen (cold-suppress turnover), 0 = deciduous
      !----- Necromass -> litter-destination split (meds_column_state_types%necromass_to_litter,   !
      !       MEDS_SLOW_DYNAMICS_DESIGN.md Part II B1). Consumed wherever plant carbon dies (leaf/    !
      !       fine-root turnover shed, continuous background mortality, cull-termination, treefall    !
      !       disturbance) to route it into the per-patch soil_carbon pools. aboveground_frac (above) !
      !       doubles as the necromass above/below-ground split for BOTH the leaf/storage and wood      !
      !       streams (fine roots are always below-ground). ------------------------------------------!
      real(wp),    allocatable :: f_labile_leaf(:)          !< [--] labile (vs structural) fraction of leaf/storage/fine-root necromass
      real(wp),    allocatable :: f_labile_stem(:)          !< [--] labile (vs structural) fraction of wood/CWD necromass
      real(wp),    allocatable :: struct_lignin_frac(:)     !< [--] lignin fraction of the structural-litter stream
      !----- Light trait-PLASTICITY slopes (meds_plant_trait_dynamics): the per-cohort leaf traits    !
      !       sla / vcmax25 / rd25 / llspan acclimate to cumulative LAI above as trait_toc*exp(k*LAI). !
      !       DERIVED from the base traits in derive_pft_rates (ED2 trait_plasticity_scheme=2, Lloyd    !
      !       et al. 2010); consumed only when [trait_dynamics].trait_plasticity_on. --------------!
      real(wp),    allocatable :: kplastic_sla(:)           !< [1/(m2/m2)] SLA light-response slope (>0)
      real(wp),    allocatable :: kplastic_vm0(:)           !< [1/(m2/m2)] Vcmax light-response slope (<0)
      real(wp),    allocatable :: kplastic_rd(:)            !< [1/(m2/m2)] Rd light-response slope (defaults to vm0)
      real(wp),    allocatable :: kplastic_llspan(:)        !< [1/(m2/m2)] leaf-lifespan light-response slope
      !----- Leaf-phenology cue params (per PFT; consumed by the slow-loop phenology advance, which  !
      !       flattens these into a meds_plant pheno_params_t). Phenology is UNCONDITIONAL now         !
      !       (docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part I): these are read only if a PFT file   !
      !       supplies a [phenology] override; otherwise left at the literature defaults installed by    !
      !       alloc_pft_table (a vanilla evergreen -- masks=CUE_NONE, realistic ~15-day flush). The       !
      !       two masks (flush/shed) select which cues drive each side; only TEMP(1)+PHOTO(8) are         !
      !       permitted until their WATER(2)/HYDRO(4)/LIGHT(16) drivers are threaded (validate_config).   !
      integer(ik), allocatable :: pheno_flush_cue_mask(:)    !< OR of CUE_* driving the flush side (min)
      integer(ik), allocatable :: pheno_shed_cue_mask(:)     !< OR of CUE_* driving the shed side (max)
      real(wp),    allocatable :: pheno_cue_sharpness(:)      !< [--]    logistic slope (large => ED2-sharp)
      real(wp),    allocatable :: pheno_k_flush_max(:)        !< [1/day] max relative flush rate (~full in 1/k days)
      real(wp),    allocatable :: pheno_k_shed_max(:)         !< [1/day] max relative active-shed rate
      real(wp),    allocatable :: pheno_tau_flush(:)          !< [day]   flush governor low-pass timescale
      real(wp),    allocatable :: pheno_tau_shed(:)           !< [day]   shed governor low-pass timescale
      real(wp),    allocatable :: pheno_gdd_base_temp(:)      !< [K]     GDD accumulation base
      real(wp),    allocatable :: pheno_chill_base_temp(:)    !< [K]     chilling-day base
      real(wp),    allocatable :: pheno_phen_a(:)             !< [K day] GDD threshold intercept (Botta 2000)
      real(wp),    allocatable :: pheno_phen_b(:)             !< [K day] GDD threshold amplitude
      real(wp),    allocatable :: pheno_phen_c(:)             !< [1/day] chilling exponent
      real(wp),    allocatable :: pheno_cold_drop_daylength(:)  !< [h] autumn short-day drop trigger (White 1997)
      real(wp),    allocatable :: pheno_cold_drop_soiltemp1(:)  !< [K] cool-soil drop (with short days)
      real(wp),    allocatable :: pheno_cold_drop_soiltemp2(:)  !< [K] very-cold-soil drop (unconditional)
      real(wp),    allocatable :: pheno_water_width(:)        !< water logistic transition width (CUE_WATER; P3)
      real(wp),    allocatable :: pheno_photo_crit(:)         !< [h]     critical daylength (CUE_PHOTO)
      real(wp),    allocatable :: pheno_photo_slope(:)        !< [1/h]   daylength logistic slope (CUE_PHOTO)
      real(wp),    allocatable :: pheno_light_on_threshold(:) !< [W/m2]  light-shed onset (CUE_LIGHT; P3)
      real(wp),    allocatable :: pheno_light_width(:)        !< [W/m2]  light-shed transition width (CUE_LIGHT; P3)
      real(wp),    allocatable :: pheno_light_window(:)       !< [day]   radiation running-mean window (CUE_LIGHT; P3)
      !----- Per-cue transition WIDTHS (were module constants in meds_phenology). --------------!
      real(wp),    allocatable :: pheno_gdd_width(:)          !< [K day] GDD flush transition width
      real(wp),    allocatable :: pheno_daylen_width(:)       !< [h]     autumn daylength transition width
      real(wp),    allocatable :: pheno_soiltemp_width(:)     !< [K]     autumn soil-temperature transition width
      !----- Baseline-turnover (degenerate phenology) controls: evergreen cold-suppression of the  !
      !       leaf/fine-root turnover shed rate, and the dormant-canopy snap-to-bare leaf fraction. !
      !       (Were module constants in meds_plant_carbon_dynamics: evg_ref_temp, evg_slope, ELONGF_MIN.)!
      real(wp),    allocatable :: pheno_evg_ref_temp(:)       !< [K]   evergreen cold-suppression reference (~5 degC)
      real(wp),    allocatable :: pheno_evg_slope(:)          !< [1/K] evergreen cold-suppression sharpness
      real(wp),    allocatable :: pheno_bare_snap_frac(:)     !< [--]  leaf fraction below which a dormant canopy snaps to bare
   end type pft_table_t

contains

   !---------------------------------------------------------------------------------------!
   ! Allocate every trait array to n PFTs.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_pft_table(pft, n)
      type(pft_table_t), intent(inout) :: pft
      integer(ik),       intent(in)    :: n
      pft%n = n
      allocate(pft%dbh_critical(n), pft%hgt_max(n), pft%wood_density(n))
      allocate(pft%reproduction_investment_fraction(n), pft%repro_carbon_efficiency(n))
      allocate(pft%mort_gamma(n), pft%mort_alpha(n), pft%mort_beta(n))
      allocate(pft%seed_rain_recruits(n), pft%include_pft(n))
      allocate(pft%photosynthetic_pathway(n), pft%vcmax25(n), pft%jmax_vcmax_ratio(n),       &
               pft%tpu_vcmax_ratio(n), pft%rd_vcmax_ratio(n), pft%kp25(n))
      allocate(pft%stomatal_g0(n), pft%stomatal_g1(n), pft%stomatal_d0(n),                   &
               pft%quantum_yield_c4(n), pft%theta_j(n), pft%theta_cj_c4(n), pft%theta_ic_c4(n))
      allocate(pft%katul_lambda25(n), pft%wstress_psi_open(n), pft%wstress_psi_close(n),     &
               pft%wstress_lambda_exp(n), pft%wstress_sref_stomata(n))
      allocate(pft%jmax25(n), pft%tpu25(n), pft%rd25(n))
      allocate(pft%wai_b1(n), pft%wai_b2(n), pft%sapwood_area_b1(n), pft%sapwood_area_b2(n))
      allocate(pft%sla(n), pft%root_to_leaf_ratio(n), pft%huber_value(n),                    &
               pft%aboveground_frac(n), pft%storage_cushion(n), pft%growth_resp_factor(n),   &
               pft%leaf_lifespan_toc(n), pft%fineroot_turnover_rate(n),                      &
               pft%wood_carbon_density(n), pft%evergreen(n))
      allocate(pft%f_labile_leaf(n), pft%f_labile_stem(n), pft%struct_lignin_frac(n))
      allocate(pft%kplastic_sla(n), pft%kplastic_vm0(n), pft%kplastic_rd(n), pft%kplastic_llspan(n))
      pft%kplastic_sla = 0.0_wp ; pft%kplastic_vm0 = 0.0_wp     ! derived in derive_pft_rates;
      pft%kplastic_rd  = 0.0_wp ; pft%kplastic_llspan = 0.0_wp  ! 0 => static (plasticity off)
      !----- Leaf-phenology cue params: allocate + install the meds_plant pheno_params_t literature   !
      !       defaults (both masks = 0 => permissive flush / no active shed = a vanilla evergreen with  !
      !       a realistic ~15-day flush -- the standard default now that phenology is UNCONDITIONAL,     !
      !       docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part I). The config loader overwrites the       !
      !       active subset per-PFT only if a [phenology] override is present in the PFT file (the        !
      !       P3-only WATER/HYDRO/LIGHT fields keep these defaults either way). ----------------------!
      allocate(pft%pheno_flush_cue_mask(n), pft%pheno_shed_cue_mask(n), pft%pheno_cue_sharpness(n),  &
               pft%pheno_k_flush_max(n), pft%pheno_k_shed_max(n), pft%pheno_tau_flush(n),            &
               pft%pheno_tau_shed(n), pft%pheno_gdd_base_temp(n), pft%pheno_chill_base_temp(n),      &
               pft%pheno_phen_a(n), pft%pheno_phen_b(n), pft%pheno_phen_c(n),                        &
               pft%pheno_cold_drop_daylength(n), pft%pheno_cold_drop_soiltemp1(n),                   &
               pft%pheno_cold_drop_soiltemp2(n), pft%pheno_water_width(n),                           &
               pft%pheno_photo_crit(n), pft%pheno_photo_slope(n),                                    &
               pft%pheno_light_on_threshold(n), pft%pheno_light_width(n), pft%pheno_light_window(n),  &
               pft%pheno_gdd_width(n), pft%pheno_daylen_width(n), pft%pheno_soiltemp_width(n),         &
               pft%pheno_evg_ref_temp(n), pft%pheno_evg_slope(n), pft%pheno_bare_snap_frac(n))
      pft%pheno_flush_cue_mask      = 0_ik         ! CUE_NONE (permissive flush)
      pft%pheno_shed_cue_mask       = 0_ik         ! CUE_NONE (no active shed)
      pft%pheno_cue_sharpness       = 2.0_wp
      pft%pheno_k_flush_max         = 0.06667_wp   ! ~ full in 15 days
      pft%pheno_k_shed_max          = 0.05_wp      ! ~ bare in 20 days
      pft%pheno_tau_flush           = 5.0_wp
      pft%pheno_tau_shed            = 5.0_wp
      pft%pheno_gdd_base_temp       = 278.15_wp
      pft%pheno_chill_base_temp     = 278.15_wp
      pft%pheno_phen_a              = -68.0_wp
      pft%pheno_phen_b              = 638.0_wp
      pft%pheno_phen_c              = -0.01_wp
      pft%pheno_cold_drop_daylength = 10.9_wp
      pft%pheno_cold_drop_soiltemp1 = 284.3_wp
      pft%pheno_cold_drop_soiltemp2 = 275.15_wp
      pft%pheno_water_width         = 0.1_wp
      pft%pheno_photo_crit          = 11.0_wp
      pft%pheno_photo_slope         = 2.0_wp
      pft%pheno_light_on_threshold  = 200.0_wp
      pft%pheno_light_width         = 50.0_wp
      pft%pheno_light_window        = 10.0_wp
      pft%pheno_gdd_width           = 50.0_wp      ! [K day] (was meds_phenology module const)
      pft%pheno_daylen_width        = 1.0_wp       ! [h]
      pft%pheno_soiltemp_width      = 2.0_wp       ! [K]
      pft%pheno_evg_ref_temp        = 278.15_wp    ! [K]   5 degC (was carbon-dynamics evg_ref_temp)
      pft%pheno_evg_slope           = 0.4_wp       ! [1/K] (was evg_slope)
      pft%pheno_bare_snap_frac      = 0.02_wp      ! [--]  (was ELONGF_MIN)
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Derive the Camac (2018) mortality-hazard parameters from wood density as a power law      !
   ! centred on mort_rho_ref (param = param_0 * (rho/rho_ref)^exp). Always positive. Call after !
   ! wood_density and the Camac coefficients are set (e.g. from the PFT config file).           !
   !---------------------------------------------------------------------------------------!
   subroutine derive_pft_rates(pft)
      type(pft_table_t), intent(inout) :: pft
      real(wp) :: lnexp(pft%n)
      real(wp), parameter :: lnexp_min = -30.0_wp, lnexp_max = 30.0_wp
      pft%mort_gamma = pft%mort_gamma_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_gamma_exp
      pft%mort_alpha = pft%mort_alpha_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_alpha_exp
      pft%mort_beta  = pft%mort_beta_0  * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_beta_exp
      !----- Light trait-plasticity slopes, ED2 trait_plasticity_scheme=2 (Lloyd et al. 2010).      !
      !       Vcmax/Rd decrease in shade (slope < 0); SLA and leaf lifespan increase. Consumed only  !
      !       when trait_plasticity_on; overridable from the [pft] config. --------------------------!
      lnexp            = max(lnexp_min, min(lnexp_max, -2.788_wp + 0.01439_wp * pft%vcmax25))
      pft%kplastic_vm0 = -exp(lnexp)                                   ! Vcmax down in the understorey
      pft%kplastic_rd  = pft%kplastic_vm0                              ! ED2 default: Rd tracks Vcmax
      pft%kplastic_sla = -pft%kplastic_vm0 * (1.18_wp / 1.10_wp)       ! SLA up (eplastic_sla/eplastic_vm0)
      !----- Leaf lifespan lengthens in shade for short-lived PFTs; the ED2 fit crosses 0 near     !
      !       ~2.6 yr. Cap at 0 so long-lived PFTs keep ONE lifespan through the canopy (never       !
      !       shorten in shade) rather than following the fit negative. ----------------------------!
      pft%kplastic_llspan = max(0.0_wp, 0.2126_wp - 0.062_wp * log(12.0_wp * pft%leaf_lifespan_toc))
   end subroutine derive_pft_rates

   !---------------------------------------------------------------------------------------!
   ! Derive the per-PFT leaf-photosynthesis capacities (Jmax25, TPU25, Rd25) as fixed ratios !
   ! of Vcmax25. Call after vcmax25 and the ratio arrays are set (from the PFT config file).   !
   !---------------------------------------------------------------------------------------!
   subroutine derive_leaf_params(pft)
      type(pft_table_t), intent(inout) :: pft
      pft%jmax25 = pft%jmax_vcmax_ratio * pft%vcmax25
      pft%tpu25  = pft%tpu_vcmax_ratio  * pft%vcmax25
      pft%rd25   = pft%rd_vcmax_ratio   * pft%vcmax25
   end subroutine derive_leaf_params

end module meds_pft_params
