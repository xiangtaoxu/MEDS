! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_pft_params -- the plant functional type (PFT) trait table, filled from the PFT file !
! (docs/configuration.md, "The PFT file").                                                 !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Each PFT's size    !
! allometry (allom, meds_allometry_lib), hydraulics (hydro, meds_plant_hydr_lib) and       !
! phenology (pheno) are records. Every trait is the PFT file's: a trait that follows a     !
! fitted relation (the sapwood conductivity or the mortality hazards from wood density,    !
! the plasticity slopes from Vcmax) is computed outside the model, by scripts/pft_traits.  !
!==========================================================================================!
module meds_pft_params
   use meds_kinds,          only : wp, ik
   use meds_allometry_lib,  only : allometry_t
   use meds_plant_hydr_lib, only : plant_hydraulics_t
   implicit none
   private

   public :: pft_table_t, phenology_traits_t
   public :: HYD_UNSET
   public :: PATH_C3, PATH_C4

   !----- Photosynthetic-pathway flags (per-PFT trait values). ----------------------------!
   integer(ik), parameter :: PATH_C3 = 1_ik   !< C3 (Farquhar-von Caemmerer-Berry)
   integer(ik), parameter :: PATH_C4 = 2_ik   !< C4 (Collatz et al. 1992)

   !----- A PFT's [phenology] block (docs/science/plant_phenology.md). The defaults equal  !
   !      pheno_params_t's (meds_phenology_types), into which the slow-loop driver copies  !
   !      them with the turgor-loss point. ------------------------------------------------!
   type :: phenology_traits_t
      integer(ik) :: flush_cue_mask        = 0_ik        !< cues of the flush side (TEMP 1, LIGHT 2, WATER 4)
      integer(ik) :: shed_cue_mask         = 0_ik        !< cues of the shed side
      real(wp)    :: flush_cue_timescale   = 5.0_wp      !< [day] smoothing of the flush tendency
      real(wp)    :: shed_cue_timescale    = 5.0_wp      !< [day] smoothing of the shed tendency
      real(wp)    :: flush_rate_max        = 0.06667_wp  !< [1/day] leaf growth at full tendency
      real(wp)    :: shed_rate_max         = 0.05_wp     !< [1/day] senescence at full tendency
      real(wp)    :: min_leaf_cover        = 0.0_wp      !< [-] leaf cover senescence stops at
      real(wp)    :: bare_leaf_cover       = 0.02_wp     !< [-] a dormant canopy below this goes bare
      real(wp)    :: flush_base_temp       = 278.15_wp   !< [K] warmth base
      real(wp)    :: flush_degree_days     = 100.0_wp    !< [K day] warmth requirement
      real(wp)    :: flush_temp_sharpness  = 0.04_wp     !< [1/(K day)]
      real(wp)    :: shed_base_temp        = 290.15_wp   !< [K] cold base
      real(wp)    :: shed_degree_days      = 50.0_wp     !< [K day] cold requirement after midsummer
      real(wp)    :: shed_temp_sharpness   = 0.1_wp      !< [1/(K day)]
      real(wp)    :: par_min               = 5.0_wp      !< [umol/m2/s] PAR that counts as an hour of light
      real(wp)    :: flush_light_hours     = 12.0_wp     !< [h/day] hours of light that permit flushing
      real(wp)    :: flush_light_sharpness = 1.0_wp      !< [1/h]
      real(wp)    :: shed_light_hours      = 11.0_wp     !< [h/day] hours of light that trigger shedding
      real(wp)    :: shed_light_sharpness  = -1.0_wp     !< [1/h]; < 0: short or dim days trigger
      real(wp)    :: light_window          = 10.0_wp     !< [day] running mean of the hours of light
      real(wp)    :: flush_water_sum       = 10.0_wp     !< [MPa day] wet sum that permits flushing
      real(wp)    :: flush_water_sharpness = 0.5_wp      !< [1/(MPa day)]
      real(wp)    :: shed_water_sum        = 10.0_wp     !< [MPa day] dry sum that triggers shedding
      real(wp)    :: shed_water_sharpness  = 0.5_wp      !< [1/(MPa day)]
   end type phenology_traits_t

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- The size allometry and the wood-density axis. --------------------------------!
      type(allometry_t),        allocatable :: allom(:) !< each PFT's allometry ([allometry])
      type(plant_hydraulics_t), allocatable :: hydro(:) !< each PFT's hydraulic traits ([hydraulics])
      real(wp), allocatable :: wood_density(:)   !< [g/cm3] rho: the wood's dry mass per volume
      !----- Reproduction and recruitment. ------------------------------------------------!
      real(wp), allocatable :: reproduction_investment_fraction(:) !< [--] growth fraction diverted to reproduction
      real(wp), allocatable :: min_repro_height(:)                 !< [m]  the height a cohort reproduces from
      real(wp), allocatable :: recruit_carbon_efficiency(:)        !< [--] share of reproduction carbon that becomes recruits
      real(wp), allocatable :: recruit_shade_decay(:)              !< [1/(m2/m2)] that share falls as exp(-k * patch LAI)
      real(wp), allocatable :: nbg_init_density(:)                 !< [plant/m2] seedlings of a bare-ground start, per patch
      !----- Mortality hazard gamma + alpha*exp(-beta*growth) (Camac et al. 2018). --------!
      real(wp), allocatable :: mort_gamma(:)     !< [1/yr]  growth-independent baseline hazard
      real(wp), allocatable :: mort_alpha(:)     !< [1/yr]  low-growth hazard magnitude
      real(wp), allocatable :: mort_beta(:)      !< [yr/cm] growth sensitivity of the hazard
      !----- Recruitment. -----------------------------------------------------------------!
      real(wp),    allocatable :: seed_rain_recruits(:) !< [plant/m2/yr] baseline external seed rain
      integer(ik), allocatable :: include_pft(:)        !< 1 = PFT may recruit, 0 = excluded
      !----- Shared height thresholds. ----------------------------------------------------!
      real(wp) :: min_cohort_height          !< [m] smallest tracked cohort; recruits born here
      !----- Leaf photosynthesis and stomata (docs/science/leaf_gas_exchange.md). ---------!
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
      !----- Co-limitation, not theta_j (docs/science/leaf_gas_exchange.md §2). -----------!
      real(wp), allocatable :: theta_cj_c3(:)        !< [--]    C3 Ac/Aj co-limitation curvature
      real(wp), allocatable :: theta_ip_c3(:)        !< [--]    C3 (Ac,Aj)/Ap co-limitation curvature
      real(wp), allocatable :: theta_cj_c4(:)        !< [--]    C4 co-limitation curvature 1
      real(wp), allocatable :: theta_ic_c4(:)        !< [--]    C4 co-limitation curvature 2
      real(wp), allocatable :: katul_lambda25(:)     !< [umol CO2/mol H2O] Katul marginal water-use efficiency
      !----- Water stress (docs/science/leaf_gas_exchange.md §4, eqs. (9)-(10)). ----------!
      real(wp), allocatable :: wstress_psi_open(:)   !< [MPa] beta_nonstomata starts to fall here; unset: the turgor-loss point
      real(wp), allocatable :: wstress_psi_close(:)  !< [MPa] ... and reaches 0 here; unset: twice the turgor-loss point
      real(wp), allocatable :: wstress_lambda_exp(:) !< [--]    Katul lambda water-stress exponent
      real(wp), allocatable :: wstress_sref_stomata(:) !< [1/MPa] Sabot slope, beta_stomata = exp(sref*(psi_pd - onset))
      real(wp), allocatable :: stomata_psi_onset(:)  !< [MPa] predawn onset of beta_stomata; unset: half the turgor-loss point
      !----- Canopy-film capacities (rain and dew), 0.1 by default (CLM5's dewmx). --------!
      real(wp), allocatable :: leaf_surf_water_max(:)   !< [kg/m2 leaf]
      real(wp), allocatable :: wood_surf_water_max(:)   !< [kg/m2 wood]
      !----- Leaf capacities derived from vcmax25 (derive_leaf_params). -------------------!
      real(wp), allocatable :: jmax25(:)             !< [umol/m2/s] DERIVED = jmax_vcmax_ratio * vcmax25
      real(wp), allocatable :: tpu25(:)              !< [umol/m2/s] DERIVED = tpu_vcmax_ratio  * vcmax25
      real(wp), allocatable :: rd25(:)               !< [umol/m2/s] DERIVED = rd_vcmax_ratio   * vcmax25
      !----- Carbon allocation, in carbon (docs/science/plant_carbon_allocation.md). ------!
      real(wp),    allocatable :: sla(:)                    !< [m2/kgC] specific leaf area
      !----- Leaf and branch sizes (docs/science/canopy_aerodynamics.md §4). --------------!
      real(wp),    allocatable :: leaf_width(:)             !< [m]  characteristic leaf width (leaf boundary layer)
      real(wp),    allocatable :: branch_diameter(:)        !< [m]  characteristic branch diameter (wood boundary layer)
      !----- Allocation targets (docs/science/plant_carbon_allocation.md §2). -------------!
      real(wp),    allocatable :: root_to_leaf_ratio(:)     !< [--]     fine-root:leaf target ratio
      real(wp),    allocatable :: storage_cushion(:)        !< [--]     storage target as a multiple of the leaf + fine-root target
      !----- Resorption, storage upkeep (docs/science/plant_carbon_allocation.md §1). -----!
      real(wp),    allocatable :: retained_carbon_fraction(:)  !< [-] share of leaf and fine-root loss resorbed
      real(wp),    allocatable :: storage_turnover_rate(:)  !< [1/yr] non-structural pool turnover, charged as respiration
      !----- Sink limits, 0 = off (docs/science/plant_carbon_allocation.md eq. (4)). ------!
      real(wp),    allocatable :: max_relative_growth_rate(:) !< [1/yr] maximum relative growth rate of wood
      real(wp),    allocatable :: max_absolute_growth_rate(:) !< [cm/yr] maximum diameter growth of a 1 cm stem
      real(wp),    allocatable :: max_absolute_growth_exponent(:) !< [--] its scaling with dbh [cm]
      real(wp),    allocatable :: growth_resp_factor(:)     !< [--]     construction cost (fraction of metabolic NPP)
      !----- Stem and root respiration (docs/science/plant_respiration.md §2-§3). ---------!
      logical,     allocatable :: is_woody(:)               !< .false. (e.g. grass) => stem respiration is 0
      real(wp),    allocatable :: stem_resp_factor25(:)     !< [umol CO2/m2 stem/s @25C] stem rate per surface
      real(wp),    allocatable :: stem_resp_sapwood25(:)    !< [umol CO2/m3 sapwood/s @25C] stem rate per sapwood volume
      real(wp),    allocatable :: root_resp_factor25(:)     !< [umol CO2/kgC fine root/s @25C] baseline root rate
      !----- Optics, LW as emissivity (docs/science/canopy_radiation_transfer.md §2). -----!
      real(wp),    allocatable :: leaf_reflect_vis(:), leaf_transmit_vis(:)   !< [--] leaf VIS
      real(wp),    allocatable :: leaf_reflect_nir(:), leaf_transmit_nir(:)   !< [--] leaf NIR
      real(wp),    allocatable :: leaf_emissivity(:)                          !< [--] leaf thermal emissivity
      real(wp),    allocatable :: wood_reflect_vis(:), wood_transmit_vis(:)   !< [--] wood/bark VIS
      real(wp),    allocatable :: wood_reflect_nir(:), wood_transmit_nir(:)   !< [--] wood/bark NIR
      real(wp),    allocatable :: wood_emissivity(:)                          !< [--] wood thermal emissivity
      real(wp),    allocatable :: leaf_clumping(:), wood_clumping(:)          !< [--] clumping factor in (0,1]
      real(wp),    allocatable :: leaf_angle_mean(:)                          !< [deg] mean leaf inclination
      real(wp),    allocatable :: leaf_angle_std(:)                           !< [deg] its standard deviation
      !----- Tissue turnover. -------------------------------------------------------------!
      real(wp),    allocatable :: leaf_turnover_rate_toc(:) !< [1/yr]   top-of-canopy baseline leaf turnover
      real(wp),    allocatable :: fineroot_turnover_rate(:) !< [1/yr]   baseline fine-root turnover
      !----- Necromass split into litter pools (docs/science/soil_carbon.md §5). ----------!
      real(wp),    allocatable :: f_labile_leaf(:)          !< [--] labile fraction of leaf, storage and fine-root necromass
      real(wp),    allocatable :: f_labile_stem(:)          !< [--] labile (vs structural) fraction of wood/CWD necromass
      real(wp),    allocatable :: struct_lignin_frac(:)     !< [--] lignin fraction of the structural-litter stream
      !----- Light-plasticity slopes (docs/science/plant_traits.md eq. (1)). --------------!
      real(wp),    allocatable :: kplastic_sla(:)           !< [1/(m2/m2)] SLA light-response slope (>0)
      real(wp),    allocatable :: kplastic_vcmax(:)         !< [1/(m2/m2)] Vcmax light-response slope (<0)
      real(wp),    allocatable :: kplastic_rd(:)            !< [1/(m2/m2)] Rd light-response slope
      real(wp),    allocatable :: kplastic_leaf_turnover(:) !< [1/(m2/m2)] leaf-turnover light-response slope (<=0)
      type(phenology_traits_t), allocatable :: pheno(:) !< each PFT's phenology ([phenology])
   end type pft_table_t

   !----- Sentinel of an optional stomatal stress point the PFT file leaves out: large and !
   !      negative, so no leaf potential is mistaken for it. ------------------------------!
   real(wp), parameter :: HYD_UNSET = -1.0e30_wp

end module meds_pft_params
