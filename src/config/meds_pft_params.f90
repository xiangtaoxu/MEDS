! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_pft_params -- plant functional type (PFT) trait table.                              !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Each PFT's size    !
! allometry is one record, allom(pft) (meds_allometry_lib): its height curve and cap,      !
! biomass, leaf area, wood area, sapwood and the share of its wood above ground. Wood      !
! density (rho) enters its biomass and is a trait of its own.                              !
!                                                                                          !
! The PFTs carry the parameters of the vital rates (meds_demography_rates, called from     !
! meds_vegetation_dynamics):                                                               !
!                                                                                          !
! GROWTH is carbon-driven, by the wood NPP the cohort allocates                            !
! (meds_plant_carbon_allocation), within the sink limits max_*_growth_*. MORTALITY is the  !
! Camac et al. (2018, PNAS) additive hazard, rate = mort_gamma +                           !
! mort_alpha*exp(-mort_beta*growth_avg), with the three parameters DERIVED from wood       !
! density (low rho => higher baseline and low-growth hazard). RECRUITMENT is a baseline    !
! external `seed_rain_recruits` plus the recruits bought by the carbon diverted to         !
! reproduction (recruit_carbon_efficiency).                                                !
!                                                                                          !
! Two height thresholds: `min_cohort_height`, shared, the smallest tracked cohort, where    !
! recruits are born; and each PFT's `min_repro_height`, the height it reproduces from.      !
!==========================================================================================!
module meds_pft_params
   use meds_kinds,         only : wp, ik
   use meds_allometry_lib, only : allometry_t
   implicit none
   private

   public :: pft_table_t, alloc_pft_table, derive_pft_rates, derive_leaf_params
   public :: HYD_UNSET
   public :: PATH_C3, PATH_C4

   !----- Photosynthetic-pathway flags (per-PFT trait values). ----------------------------!
   integer(ik), parameter :: PATH_C3 = 1_ik   !< C3 (Farquhar-von Caemmerer-Berry)
   integer(ik), parameter :: PATH_C4 = 2_ik   !< C4 (Collatz et al. 1992)

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- The size allometry and the wood-density axis. -------------------------------------!
      type(allometry_t), allocatable :: allom(:) !< each PFT's allometry ([allometry])
      real(wp), allocatable :: wood_density(:)   !< [g/cm3] rho: AGB + mortality anchor
      !----- Reproductive allocation and recruitment (carbon path). An empirical growth law is     !
      !      not a trait here; it lives in the Python example. ------------------------------------!
      real(wp), allocatable :: reproduction_investment_fraction(:) !< [--] growth fraction diverted to reproduction
      real(wp), allocatable :: min_repro_height(:)                 !< [m]  the height a cohort reproduces from
      real(wp), allocatable :: recruit_carbon_efficiency(:)        !< [--] share of reproduction carbon that becomes
                                                                   !<      recruits: seed set, germination and seedling
                                                                   !<      survival up to min_cohort_height together
      real(wp), allocatable :: recruit_shade_decay(:)              !< [1/(m2/m2)] that share falls as exp(-k * patch LAI)
      !----- A near-bare-ground start's seedlings (ED2's init_density): plants per m2 of each patch, at  !
      !      min_cohort_height on the polygon's height curve. 0 = no cohort of the PFT (the empty start). !
      real(wp), allocatable :: nbg_init_density(:)                 !< [plant/m2]
      !----- Wood-density -> mortality-hazard DERIVATION coefficients (Camac et al. 2018         !
      !      PNAS, a power law centred on mort_rho_ref: param = param_0 * (rho/rho_ref)^exp,     !
      !      always positive). Shared scalars, one set for every PFT. ---------------------------!
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
      !----- Leaf photosynthesis / stomatal traits (per PFT; docs/science/leaf_gas_exchange.md). -!
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
      !----- CO-LIMITATION curvatures, one pair per pathway. These describe how sharply the        !
      !      leaf transitions between two LIMITING PROCESSES and are physically distinct from      !
      !      theta_j, the curvature of the electron-transport hyperbola. A J curvature of 0.85     !
      !      on both costs 27 % of gross assimilation at saturating light and Ci 280 (31-34 % of   !
      !      net assimilation at PAR 1500-2000), because A_c and A_j sit close together at ambient !
      !      CO2 (docs/science/leaf_gas_exchange.md §2, "The co-limitation curvatures are not      !
      !      theta_J"). ---------------------------------------------------------------------------!
      real(wp), allocatable :: theta_cj_c3(:)        !< [--]    C3 Ac/Aj co-limitation curvature
      real(wp), allocatable :: theta_ip_c3(:)        !< [--]    C3 (Ac,Aj)/Ap co-limitation curvature
      real(wp), allocatable :: theta_cj_c4(:)        !< [--]    C4 co-limitation curvature 1
      real(wp), allocatable :: theta_ic_c4(:)        !< [--]    C4 co-limitation curvature 2
      real(wp), allocatable :: katul_lambda25(:)     !< [umol CO2/mol H2O] Katul marginal water-use efficiency
      !< [MPa] current leaf potential where the fast limb beta_nonstomata starts to fall (open) and where it
      !<       reaches 0 (close), docs/science/leaf_gas_exchange.md eq. (9). OPTIONAL: HYD_UNSET means the
      !<       PFT's turgor-loss point and twice it (pft_wstress_psi_*)
      real(wp), allocatable :: wstress_psi_open(:)
      real(wp), allocatable :: wstress_psi_close(:)
      real(wp), allocatable :: wstress_lambda_exp(:) !< [--]    Katul lambda water-stress exponent
      !< [1/MPa] Sabot stomatal-stress sensitivity: beta_stomata = exp(sref*(psi - psi_onset)) below the
      !<         onset, psi = predawn leaf psi (docs/science/leaf_gas_exchange.md eq. (10))
      real(wp), allocatable :: wstress_sref_stomata(:)
      !< [MPa] predawn leaf psi where the stomatal stress begins (beta_stomata = 1 above it). OPTIONAL:
      !<       HYD_UNSET means half the PFT's leaf turgor-loss point (pft_stomata_psi_onset)
      real(wp), allocatable :: stomata_psi_onset(:)
      !----- The canopy films' capacities (interception and dew): the water a leaf or a woody surface holds
      !      before it drips. OPTIONAL, 0.1 by default (CLM5's dewmx, per unit leaf or stem area).
      real(wp), allocatable :: leaf_surf_water_max(:)   !< [kg/m2 leaf]
      real(wp), allocatable :: wood_surf_water_max(:)   !< [kg/m2 wood]
      !----- Leaf photosynthesis DERIVED per-PFT (derive_leaf_params). ------------------------!
      real(wp), allocatable :: jmax25(:)             !< [umol/m2/s] DERIVED = jmax_vcmax_ratio * vcmax25
      real(wp), allocatable :: tpu25(:)              !< [umol/m2/s] DERIVED = tpu_vcmax_ratio  * vcmax25
      real(wp), allocatable :: rd25(:)               !< [umol/m2/s] DERIVED = rd_vcmax_ratio   * vcmax25
      !----- Carbon-allocation per-PFT traits (meds_plant_carbon_allocation;                       !
      !      docs/science/plant_carbon_allocation.md). All state is CARBON; every                  !
      !      biomass<->carbon conversion is folded into these traits at init. -------------------!
      real(wp),    allocatable :: sla(:)                    !< [m2/kgC] specific leaf area (un-folds SLA from lai_b1)
      !----- Canopy-element geometry the aerodynamic boundary layers read                          !
      !      (docs/science/canopy_aerodynamics.md §4). Per PFT, because leaf size is one of        !
      !      the strongest interspecific axes there is. -------------------------------------------!
      real(wp),    allocatable :: leaf_width(:)             !< [m]  characteristic leaf width (leaf boundary layer)
      real(wp),    allocatable :: branch_diameter(:)        !< [m]  characteristic branch diameter (wood boundary layer)
      real(wp),    allocatable :: root_to_leaf_ratio(:)     !< [--]     fine-root:leaf target ratio (ED2 q)
      real(wp),    allocatable :: huber_value(:)            !< [m2 sap/m2 leaf] sapwood-area:leaf-area (sapwood + hydraulics)
      real(wp),    allocatable :: storage_cushion(:)        !< [--]     storage target as a multiple of the leaf + fine-root target
      !----- STORAGE MAINTENANCE. Fractional turnover of the non-structural pool [1/yr], charged as
      !      respiration and taken from the pool -- ED2's growth_balive.f90 form, with no temperature
      !      dependence, as in ED2 (docs/science/plant_carbon_allocation.md eq. (2)). DEFAULT 0.0; ED2's
      !      own values range from 0 to 0.6243 by PFT, so zero is a legitimate ED2 setting, not an absence
      !      of physics. Whether to default it on is open (see docs/ROADMAP.md).
      !----- RESORPTION. Fraction of all leaf and fine-root loss, senescence and background turnover
      !      alike, returned to the non-structural pool instead of entering litter
      !      (docs/science/plant_carbon_allocation.md eq. (3)). DEFAULT 0.3.
      !----- PER-PFT HYDRAULIC TRAITS, so PFTs differ hydraulically on more than wood density. All
      !      OPTIONAL: an absent key takes the [hydraulics] scalar (wood_kmax, which has none, follows
      !      the PFT's wood density), so a config can make just ONE trait per-PFT without restating
      !      the others.
      real(wp),    allocatable :: hyd_leaf_pi0(:), hyd_leaf_elastic_mod(:), hyd_leaf_apoplast_frac(:)
      real(wp),    allocatable :: hyd_leaf_water_sat(:)
      real(wp),    allocatable :: hyd_wood_pi0(:), hyd_wood_elastic_mod(:), hyd_wood_apoplast_frac(:)
      real(wp),    allocatable :: hyd_wood_water_sat(:)
      real(wp),    allocatable :: hyd_wood_psi50(:), hyd_wood_kexp(:)
      real(wp),    allocatable :: hyd_wood_kmax(:), hyd_vessel_curl(:)
      real(wp),    allocatable :: hyd_root_beta(:), hyd_root_depth_b1(:), hyd_root_depth_b2(:)
      real(wp),    allocatable :: hyd_specific_root_length(:), hyd_fine_root_radius(:)
      real(wp),    allocatable :: retained_carbon_fraction(:)  !< [-] share of leaf and fine-root loss resorbed
      real(wp),    allocatable :: storage_turnover_rate(:)  !< [1/yr] non-structural pool turnover
      !----- SINK limitation (meds_sink_limitation): wood carbon grows by at most this fraction of  !
      !      itself a year, and the stem diameter by at most rate x dbh^exponent a year; the tighter !
      !      one applies, and the carbon left over is exuded by the roots. A rate of 0 turns that    !
      !      limit off (default); the exponent defaults to 0, a size-independent diameter cap.        !
      real(wp),    allocatable :: max_relative_growth_rate(:) !< [1/yr] maximum relative growth rate of wood
      real(wp),    allocatable :: max_absolute_growth_rate(:) !< [cm/yr] maximum diameter growth of a 1 cm stem
      real(wp),    allocatable :: max_absolute_growth_exponent(:) !< [--] its scaling with dbh [cm]
      real(wp),    allocatable :: growth_resp_factor(:)     !< [--]     construction cost (fraction of metabolic NPP)
      !----- MAINTENANCE respiration of the non-leaf tissues                                      !
      !      (docs/science/plant_respiration.md §2-§3). Per-PFT, as in ED2, so PFTs can           !
      !      respire their stems and roots differently. (Leaf dark respiration is per-PFT         !
      !      through rd_vcmax_ratio.) ------------------------------------------------------------!
      logical,     allocatable :: is_woody(:)               !< .false. (e.g. grass) => stem respiration is 0
      real(wp),    allocatable :: stem_resp_factor25(:)     !< [umol CO2/m2 stem/s @25C] stem rate per surface
      real(wp),    allocatable :: stem_resp_sapwood25(:)    !< [umol CO2/m3 sapwood/s @25C] stem rate per sapwood volume
      real(wp),    allocatable :: root_resp_factor25(:)     !< [umol CO2/kgC fine root/s @25C] baseline root rate
      !----- CANOPY OPTICS (the two-stream's per-PFT table;                                       !
      !      docs/science/canopy_radiation_transfer.md §2). Per PFT, so two PFTs can differ       !
      !      in how they intercept or scatter light -- what the two-stream solver exists to       !
      !      resolve. Shortwave is given as reflectance + transmittance per band; LONGWAVE is     !
      !      given as EMISSIVITY, because that is the number anyone has for a leaf, and the       !
      !      band's reflectance is 1 - emissivity with zero transmittance (a leaf is opaque       !
      !      at thermal wavelengths -- physics, not a knob). -------------------------------------!
      real(wp),    allocatable :: leaf_reflect_vis(:), leaf_transmit_vis(:)   !< [--] leaf VIS
      real(wp),    allocatable :: leaf_reflect_nir(:), leaf_transmit_nir(:)   !< [--] leaf NIR
      real(wp),    allocatable :: leaf_emissivity(:)                          !< [--] leaf thermal emissivity
      real(wp),    allocatable :: wood_reflect_vis(:), wood_transmit_vis(:)   !< [--] wood/bark VIS
      real(wp),    allocatable :: wood_reflect_nir(:), wood_transmit_nir(:)   !< [--] wood/bark NIR
      real(wp),    allocatable :: wood_emissivity(:)                          !< [--] wood thermal emissivity
      real(wp),    allocatable :: leaf_clumping(:), wood_clumping(:)          !< [--] clumping factor in (0,1]
      real(wp),    allocatable :: leaf_angle_mean(:)                          !< [deg] mean leaf inclination
      real(wp),    allocatable :: leaf_angle_std(:)                           !< [deg] its standard deviation
      real(wp),    allocatable :: leaf_turnover_rate_toc(:) !< [1/yr]   top-of-canopy baseline leaf turnover
      real(wp),    allocatable :: fineroot_turnover_rate(:) !< [1/yr]   baseline fine-root turnover
      real(wp),    allocatable :: wood_carbon_density(:)    !< [kgC/m3] wood carbon density (Huber sapwood carbon)
      !----- Necromass -> litter-destination split (meds_litter_partition%necromass_to_litter;      !
      !      docs/science/soil_carbon.md §5). Consumed wherever plant carbon dies (leaf/fine-       !
      !      root turnover shed, continuous background mortality, cull-termination, treefall        !
      !      disturbance) to route it into the per-patch soil_carbon pools. The allometry's         !
      !      aboveground_frac splits the wood and storage streams above/below ground; leaves and    !
      !      reproduction are all above, fine roots all below. -------------------------------------!
      !< [--] labile (vs structural) fraction of leaf/storage/fine-root necromass
      real(wp),    allocatable :: f_labile_leaf(:)
      real(wp),    allocatable :: f_labile_stem(:)          !< [--] labile (vs structural) fraction of wood/CWD necromass
      real(wp),    allocatable :: struct_lignin_frac(:)     !< [--] lignin fraction of the structural-litter stream
      !----- Light trait-PLASTICITY slopes (meds_plant_trait_dynamics): the per-cohort leaf        !
      !      traits sla / vcmax25 / rd25 / leaf_turnover_rate acclimate to cumulative LAI above as !
      !      trait_toc*exp(k*LAI) (docs/science/plant_traits.md eq. (1)). DERIVED from the         !
      !      base traits in derive_pft_rates (ED2 trait_plasticity_scheme=2, Lloyd et al.          !
      !      2010; docs/science/plant_traits.md eqs. (3)-(4)); consumed only when                  !
      !      [trait_dynamics].trait_plasticity_on. ------------------------------------------------!
      real(wp),    allocatable :: kplastic_sla(:)           !< [1/(m2/m2)] SLA light-response slope (>0)
      real(wp),    allocatable :: kplastic_vm0(:)           !< [1/(m2/m2)] Vcmax light-response slope (<0)
      real(wp),    allocatable :: kplastic_rd(:)            !< [1/(m2/m2)] Rd light-response slope (defaults to vm0)
      real(wp),    allocatable :: kplastic_leaf_turnover(:) !< [1/(m2/m2)] leaf-turnover light-response slope (<=0)
      !----- Leaf-phenology cue params (per PFT; flattened into a pheno_params_t by the slow-loop  !
      !       phenology advance). Read only from a PFT file's [phenology] section; without one a PFT  !
      !       keeps the defaults alloc_pft_table installs: both masks CUE_NONE, i.e. always flushing  !
      !       and never senescing, with leaf loss from background turnover alone. Every switch is    !
      !       sigma(s (x - x*)): a centre and a signed sharpness (docs/science/plant_phenology.md). !
      !       The water cue's turgor-loss point is derived (pft_leaf_psi_tlp), not a key here. ----!
      integer(ik), allocatable :: pheno_flush_cue_mask(:)        !< cues of the flush side (TEMP 1, LIGHT 2, WATER 4)
      integer(ik), allocatable :: pheno_shed_cue_mask(:)         !< cues of the shed side
      real(wp),    allocatable :: pheno_flush_cue_timescale(:)   !< [day]   smoothing of the flush tendency
      real(wp),    allocatable :: pheno_shed_cue_timescale(:)    !< [day]   smoothing of the shed tendency
      real(wp),    allocatable :: pheno_flush_rate_max(:)        !< [1/day] leaf growth at full tendency
      real(wp),    allocatable :: pheno_shed_rate_max(:)         !< [1/day] senescence at full tendency
      real(wp),    allocatable :: pheno_flush_base_temp(:)       !< [K]     warmth base
      real(wp),    allocatable :: pheno_flush_degree_days(:)     !< [K day] warmth requirement
      real(wp),    allocatable :: pheno_flush_temp_sharpness(:)  !< [1/(K day)]
      real(wp),    allocatable :: pheno_shed_base_temp(:)        !< [K]     cold base
      real(wp),    allocatable :: pheno_shed_degree_days(:)      !< [K day] cold requirement after midsummer
      real(wp),    allocatable :: pheno_shed_temp_sharpness(:)   !< [1/(K day)]
      real(wp),    allocatable :: pheno_par_min(:)               !< [umol/m2/s] PAR at the cohort's top that counts as light
      real(wp),    allocatable :: pheno_flush_light_hours(:)     !< [h/day] hours of light that permit flushing
      real(wp),    allocatable :: pheno_flush_light_sharpness(:) !< [1/h]
      real(wp),    allocatable :: pheno_shed_light_hours(:)      !< [h/day] hours of light that trigger shedding
      real(wp),    allocatable :: pheno_shed_light_sharpness(:)  !< [1/h]; < 0: short or dim days trigger, > 0 bright
      real(wp),    allocatable :: pheno_light_window(:)          !< [day]   running mean of the hours of light
      real(wp),    allocatable :: pheno_flush_water_sum(:)       !< [MPa day] wet sum that permits flushing
      real(wp),    allocatable :: pheno_flush_water_sharpness(:) !< [1/(MPa day)]
      real(wp),    allocatable :: pheno_shed_water_sum(:)        !< [MPa day] dry sum that triggers shedding
      real(wp),    allocatable :: pheno_shed_water_sharpness(:)  !< [1/(MPa day)]
      real(wp),    allocatable :: pheno_min_leaf_cover(:)        !< [--] leaf cover senescence stops at
      real(wp),    allocatable :: pheno_bare_leaf_cover(:)       !< [--] a dormant canopy below this goes bare
   end type pft_table_t

   !----- 'take the [hydraulics] scalar' sentinel. Large and negative so it cannot be mistaken for
   !      any physical hydraulic trait -- psi50 and pi0 are the only negative ones and neither
   !      reaches -1e30.
   real(wp), parameter :: HYD_UNSET = -1.0e30_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Allocate every trait array to n PFTs.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_pft_table(pft, n)
      type(pft_table_t), intent(inout) :: pft
      integer(ik),       intent(in)    :: n
      pft%n = n
      allocate(pft%allom(n), pft%wood_density(n))
      allocate(pft%reproduction_investment_fraction(n), pft%min_repro_height(n), pft%recruit_carbon_efficiency(n))
      allocate(pft%recruit_shade_decay(n))
      pft%recruit_shade_decay = 0.0_wp      ! optional key; 0 = the same share under any canopy
      allocate(pft%nbg_init_density(n))
      pft%nbg_init_density = 0.0_wp         ! optional key; 0 = bare ground starts empty
      allocate(pft%stem_resp_sapwood25(n))
      pft%stem_resp_sapwood25 = 0.0_wp      ! optional key; 0 = stem respiration by surface alone
      allocate(pft%mort_gamma(n), pft%mort_alpha(n), pft%mort_beta(n))
      allocate(pft%seed_rain_recruits(n), pft%include_pft(n))
      allocate(pft%photosynthetic_pathway(n), pft%vcmax25(n), pft%jmax_vcmax_ratio(n),       &
               pft%tpu_vcmax_ratio(n), pft%rd_vcmax_ratio(n), pft%kp25(n))
      allocate(pft%stomatal_g0(n), pft%stomatal_g1(n), pft%stomatal_d0(n),                   &
               pft%quantum_yield_c4(n), pft%theta_j(n), pft%theta_cj_c4(n), pft%theta_ic_c4(n),  &
               pft%theta_cj_c3(n), pft%theta_ip_c3(n))
      allocate(pft%katul_lambda25(n), pft%wstress_psi_open(n), pft%wstress_psi_close(n),     &
               pft%wstress_lambda_exp(n), pft%wstress_sref_stomata(n), pft%stomata_psi_onset(n))
      allocate(pft%leaf_surf_water_max(n), pft%wood_surf_water_max(n))
      allocate(pft%jmax25(n), pft%tpu25(n), pft%rd25(n))
      allocate(pft%leaf_width(n), pft%branch_diameter(n))
      allocate(pft%sla(n), pft%root_to_leaf_ratio(n), pft%huber_value(n),                    &
               pft%storage_cushion(n), pft%growth_resp_factor(n),                            &
               pft%storage_turnover_rate(n), pft%retained_carbon_fraction(n),               &
               pft%hyd_leaf_pi0(n), pft%hyd_leaf_elastic_mod(n), pft%hyd_leaf_apoplast_frac(n), &
               pft%hyd_leaf_water_sat(n), pft%hyd_wood_pi0(n), pft%hyd_wood_elastic_mod(n),     &
               pft%hyd_wood_apoplast_frac(n), pft%hyd_wood_water_sat(n),                        &
               pft%hyd_wood_psi50(n), pft%hyd_wood_kexp(n),                                     &
               pft%hyd_wood_kmax(n), pft%hyd_vessel_curl(n),                                    &
               pft%is_woody(n), pft%stem_resp_factor25(n), pft%root_resp_factor25(n),        &
               pft%leaf_reflect_vis(n), pft%leaf_transmit_vis(n), pft%leaf_reflect_nir(n),   &
               pft%leaf_transmit_nir(n), pft%leaf_emissivity(n),                             &
               pft%wood_reflect_vis(n), pft%wood_transmit_vis(n), pft%wood_reflect_nir(n),   &
               pft%wood_transmit_nir(n), pft%wood_emissivity(n),                             &
               pft%leaf_clumping(n), pft%wood_clumping(n),                                   &
               pft%leaf_angle_mean(n), pft%leaf_angle_std(n),                                &
               pft%leaf_turnover_rate_toc(n), pft%fineroot_turnover_rate(n),                  &
               pft%wood_carbon_density(n))
      pft%storage_turnover_rate = 0.0_wp   ! optional key; 0 = no storage maintenance
      allocate(pft%max_relative_growth_rate(n), pft%max_absolute_growth_rate(n),                     &
               pft%max_absolute_growth_exponent(n))
      pft%max_relative_growth_rate     = 0.0_wp   ! optional keys; a rate of 0 = that sink limit off
      pft%max_absolute_growth_rate     = 0.0_wp
      pft%max_absolute_growth_exponent = 0.0_wp
      pft%retained_carbon_fraction = 0.3_wp ! optional key; 30% of leaf and root loss resorbed
      !----- HYD_UNSET marks "no per-PFT value given"; the table builder then takes the            !
      !      [hydraulics] scalar. A sentinel rather than pre-filling from [hydraulics] here,       !
      !      because alloc_pft_table runs before the [hydraulics] block is read. ------------------!
      pft%hyd_leaf_pi0 = HYD_UNSET ; pft%hyd_leaf_elastic_mod = HYD_UNSET
      pft%stomata_psi_onset = HYD_UNSET   ! the same sentinel: half the PFT's turgor-loss point
      pft%wstress_psi_open  = HYD_UNSET   ! ... and the turgor-loss point / twice it
      pft%wstress_psi_close = HYD_UNSET
      pft%leaf_surf_water_max = 0.1_wp ; pft%wood_surf_water_max = 0.1_wp   ! optional keys: CLM5's dewmx
      pft%hyd_leaf_apoplast_frac = HYD_UNSET ; pft%hyd_leaf_water_sat = HYD_UNSET
      pft%hyd_wood_pi0 = HYD_UNSET ; pft%hyd_wood_elastic_mod = HYD_UNSET
      pft%hyd_wood_apoplast_frac = HYD_UNSET ; pft%hyd_wood_water_sat = HYD_UNSET
      pft%hyd_wood_psi50 = HYD_UNSET ; pft%hyd_wood_kexp = HYD_UNSET
      pft%hyd_wood_kmax = HYD_UNSET ; pft%hyd_vessel_curl = HYD_UNSET
      allocate(pft%hyd_root_beta(n), pft%hyd_root_depth_b1(n), pft%hyd_root_depth_b2(n),               &
               pft%hyd_specific_root_length(n), pft%hyd_fine_root_radius(n))
      pft%hyd_root_beta = HYD_UNSET ; pft%hyd_root_depth_b1 = HYD_UNSET ; pft%hyd_root_depth_b2 = HYD_UNSET
      pft%hyd_specific_root_length = HYD_UNSET ; pft%hyd_fine_root_radius = HYD_UNSET
      allocate(pft%f_labile_leaf(n), pft%f_labile_stem(n), pft%struct_lignin_frac(n))
      allocate(pft%kplastic_sla(n), pft%kplastic_vm0(n), pft%kplastic_rd(n), pft%kplastic_leaf_turnover(n))
      pft%kplastic_sla = 0.0_wp ; pft%kplastic_vm0 = 0.0_wp            ! derived in derive_pft_rates;
      pft%kplastic_rd  = 0.0_wp ; pft%kplastic_leaf_turnover = 0.0_wp  ! 0 => static (plasticity off)
      !----- Leaf-phenology cue params: the defaults a PFT without a [phenology] section keeps --   !
      !       both masks CUE_NONE, so only the flush rate, the cue timescales and the background    !
      !       turnover act. They equal the pheno_params_t defaults (meds_phenology_types). ----------!
      allocate(pft%pheno_flush_cue_mask(n), pft%pheno_shed_cue_mask(n),                              &
               pft%pheno_flush_cue_timescale(n), pft%pheno_shed_cue_timescale(n),                    &
               pft%pheno_flush_rate_max(n), pft%pheno_shed_rate_max(n),                              &
               pft%pheno_flush_base_temp(n), pft%pheno_flush_degree_days(n),                         &
               pft%pheno_flush_temp_sharpness(n), pft%pheno_shed_base_temp(n),                       &
               pft%pheno_shed_degree_days(n), pft%pheno_shed_temp_sharpness(n),                      &
               pft%pheno_par_min(n), pft%pheno_flush_light_hours(n), pft%pheno_flush_light_sharpness(n), &
               pft%pheno_shed_light_hours(n), pft%pheno_shed_light_sharpness(n), pft%pheno_light_window(n), &
               pft%pheno_flush_water_sum(n), pft%pheno_flush_water_sharpness(n),                     &
               pft%pheno_shed_water_sum(n), pft%pheno_shed_water_sharpness(n),                       &
               pft%pheno_min_leaf_cover(n), pft%pheno_bare_leaf_cover(n))
      pft%pheno_flush_cue_mask        = 0_ik         ! CUE_NONE: always flushing
      pft%pheno_shed_cue_mask         = 0_ik         ! CUE_NONE: no senescence
      pft%pheno_flush_cue_timescale   = 5.0_wp
      pft%pheno_shed_cue_timescale    = 5.0_wp
      pft%pheno_flush_rate_max        = 0.06667_wp   ! a bare canopy fills in ~15 days
      pft%pheno_shed_rate_max         = 0.05_wp
      pft%pheno_flush_base_temp       = 278.15_wp
      pft%pheno_flush_degree_days     = 100.0_wp
      pft%pheno_flush_temp_sharpness  = 0.04_wp
      pft%pheno_shed_base_temp        = 290.15_wp
      pft%pheno_shed_degree_days      = 50.0_wp
      pft%pheno_shed_temp_sharpness   = 0.1_wp
      pft%pheno_par_min               = 5.0_wp       ! a low threshold: the hours of light ~ the day length
      pft%pheno_flush_light_hours     = 12.0_wp
      pft%pheno_flush_light_sharpness = 1.0_wp
      pft%pheno_shed_light_hours      = 11.0_wp
      pft%pheno_shed_light_sharpness  = -1.0_wp
      pft%pheno_light_window          = 10.0_wp
      pft%pheno_flush_water_sum       = 10.0_wp
      pft%pheno_flush_water_sharpness = 0.5_wp
      pft%pheno_shed_water_sum        = 10.0_wp
      pft%pheno_shed_water_sharpness  = 0.5_wp
      pft%pheno_min_leaf_cover        = 0.0_wp
      pft%pheno_bare_leaf_cover       = 0.02_wp
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Derive the Camac (2018) mortality-hazard parameters from wood density as a power law      !
   ! centred on mort_rho_ref (param = param_0 * (rho/rho_ref)^exp). Always positive. Call after !
   ! wood_density and the Camac coefficients are set (e.g. from the PFT config file).           !
   !---------------------------------------------------------------------------------------!
   subroutine derive_pft_rates(pft)
      type(pft_table_t), intent(inout) :: pft
      real(wp) :: lnexp(pft%n)
      real(wp), parameter :: lnexp_min_pft = -30.0_wp, lnexp_max = 30.0_wp
      pft%mort_gamma = pft%mort_gamma_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_gamma_exp
      pft%mort_alpha = pft%mort_alpha_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_alpha_exp
      pft%mort_beta  = pft%mort_beta_0  * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_beta_exp
      !----- Light trait-plasticity slopes, ED2 trait_plasticity_scheme=2 (Lloyd et al. 2010).      !
      !       Vcmax/Rd and leaf turnover decrease in shade (slope < 0); SLA increases. Consumed only !
      !       when trait_plasticity_on; with override_derived the PFT file's [derived] table can    !
      !       pin each one. ------------------------------------------------------------------------!
      lnexp            = max(lnexp_min_pft, min(lnexp_max, -2.788_wp + 0.01439_wp * pft%vcmax25))
      pft%kplastic_vm0 = -exp(lnexp)                                   ! Vcmax down in the understorey
      pft%kplastic_rd  = pft%kplastic_vm0                              ! ED2 default: Rd tracks Vcmax
      pft%kplastic_sla = -pft%kplastic_vm0 * (1.18_wp / 1.10_wp)       ! SLA up (eplastic_sla/eplastic_vm0)
      !----- Leaf turnover slows in shade for fast-turnover PFTs: the slope is minus ED2's fitted     !
      !       lifespan slope, which falls with ln(12 / rate) (the lifespan in months) and crosses 0 near !
      !       a rate of 0.39 /yr. Cap at 0 so slow-turnover PFTs keep ONE rate through the canopy        !
      !       (never turn over faster in shade). docs/science/plant_traits.md eq. (4). ----------------!
      pft%kplastic_leaf_turnover = min(0.0_wp, 0.062_wp * log(12.0_wp / pft%leaf_turnover_rate_toc) - 0.2126_wp)
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
