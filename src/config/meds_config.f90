! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config -- the run configuration, threaded read-only through the engine              !
! (docs/configuration.md).                                                                 !
!                                                                                          !
! One field per key, named as the key, in the order of meds_config_main.toml: a key of a   !
! block with its own record is cfg%<record>%<key> (cfg%init%soil_temp), and a key of the   !
! other blocks is cfg%<key>. The two exceptions are the [output] group, axis               !
! and tier switches, which are arrays indexed by their codes, and the [hydraulics] tissue  !
! curves (leaf_pi0 is hydro(pft)%leaf_curve%pi0). The reference config describes each key; !
! a field without an initializer is a required key or a derived value (marked DERIVED).    !
!                                                                                          !
! The routines derive what the parameters imply: the run's derived values, the PFT table's !
! arrays and leaf capacities, and the soil column (its layers and the bundles the soil     !
! kernels read), built once per column: each site and each map cell.                       !
!==========================================================================================!
module meds_config
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : yr_day, yr_sec, n_soil_layer_max
   use meds_time,               only : meds_time_t
   use meds_pft_params,         only : pft_table_t, HYD_UNSET
   use meds_plant_hydr_lib,     only : tissue_psi_tlp
   use meds_soil_hydr_lib,      only : soil_head_to_theta, soil_matric_potential
   use meds_fast_opts,          only : process_mask_t, INTEG_ARK, CTRL_I, CTRL_L1_ADAPTIVE
   use meds_soil_column_config, only : soil_column_config_t, soil_texture_layers_t, soil_params_t,   &
                                       soil_thermal_params_t
   use meds_biophysics_opts,    only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biogeochem_opts,    only : decomp_opts_t
   use meds_forcing_opts,       only : forcing_config_t
   use meds_output_opts,        only : output_config_t
   use meds_region_opts,        only : region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   use meds_maps_opts,          only : maps_opts_t, soil_maps_on
   use meds_sites_opts,         only : site_spec_t
   implicit none
   private

   public :: meds_config_t, init_opts_t, state_opts_t
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: derive_config, derive_parameters, growth_window_steps, site_config
   public :: pft_leaf_psi_tlp, pft_stomata_psi_onset, pft_wstress_psi_open, pft_wstress_psi_close
   public :: forcing_config_t, output_config_t
   public :: decomp_opts_t
   public :: region_opts_t, RUN_MODE_SITE, RUN_MODE_REGION
   public :: maps_opts_t, soil_maps_on
   public :: alloc_pft_table, derive_leaf_params
   public :: soil_column_layers, soil_column_texture, soil_column_params, soil_column_start_theta
   public :: top_soil_psi, SOIL_LAYER_TEMPLATE, soil_layer_grid
   public :: build_soil_hydr_params, build_soil_therm_params

   !----- [init].init_mode. -------------------------------------------------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< a cohort census (census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< a state file (restart_file)

   !----- The default layer bottoms [m, down], cut at each column's depth: SoilGrids 2.0's !
   !      depths to 1 m, then thickening to 27 m (docs/science/soil_biophysics.md, "The    !
   !      soil layers"). ------------------------------------------------------------------!
   real(wp), parameter :: SOIL_LAYER_TEMPLATE(15) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp,     &
                          1.50_wp, 2.00_wp, 3.00_wp, 4.50_wp, 6.50_wp, 9.00_wp, 12.0_wp, 16.0_wp, 21.0_wp, 27.0_wp]
   !----- A remainder below the cut thinner than this share of the list layer it falls in  !
   !      joins the layer above, so no column ends in a sliver. ---------------------------!
   real(wp), parameter :: LAYER_MERGE_SHARE = 0.5_wp
   real(wp), parameter :: GRID_TOL = 1.0e-6_wp     !< [m] a depth this close to a list bottom ends on it

   !----- [init]: how the run starts. -----------------------------------------------------------!
   type :: init_opts_t
      integer(ik)        :: init_mode                     !< INIT_BARE | INIT_CENSUS | INIT_RESTART
      character(len=256) :: restart_file = ''             !< the state file init_mode = 2 reads
      character(len=256) :: census_file  = ''             !< the cohort CSV init_mode = 1 reads
      character(len=256) :: pft_config                    !< the PFT file
      real(wp)           :: soil_temp    = 288.0_wp       !< [K] every soil layer's start
      real(wp)           :: soil_theta   = 0.30_wp        !< [m3/m3] every soil layer's start
      logical            :: reacclimate_traits = .false.  !< a restart's plastic traits from this PFT file
   end type init_opts_t

   !----- [state]: the restart stream. ------------------------------------------------------------!
   type :: state_opts_t
      character(len=256) :: output_dir                    !< directory of the state files and the PFT CSV
      character(len=256) :: output_prefix                 !< their file-name stem
      logical            :: write_state                   !< write state checkpoints
      integer(ik)        :: interval_years                !< [yr] between checkpoints
   end type state_opts_t

   type :: meds_config_t
      !----- [run]. ---------------------------------------------------------------------------!
      real(wp)          :: dt_slow                        !< [s] slow step
      type(meds_time_t) :: start_time, end_time           !< the run's span
      integer(ik)       :: mode           = RUN_MODE_SITE !< RUN_MODE_SITE | RUN_MODE_REGION
      integer(ik)       :: n_threads      = 1_ik          !< host threads
      logical           :: slow_on        = .true.        !< .false.: the slow tier is skipped
      logical           :: slow_ledger_on = .true.        !< the slow-step conservation ledger
      real(wp)          :: dt_years                       !< [yr] DERIVED dt_slow / yr_sec

      !----- [fast]. --------------------------------------------------------------------------!
      logical     :: fast_biophysics_on                   !< the sub-daily loop
      real(wp)    :: dt_fast                              !< [s] sub-daily step
      integer(ik) :: time_integrator     = INTEG_ARK      !< INTEG_ARK | INTEG_RK45
      logical     :: ark_adaptive        = .true.         !< adaptive sub-steps; .false.: ark_fixed_substep
      integer(ik) :: ark_fixed_substep   = 4_ik           !< sub-steps per dt_fast when not adaptive
      real(wp)    :: ark_dt_init         = 0.0_wp         !< [s] first adaptive sub-step (0: warm start)
      real(wp)    :: rtol                = 1.0e-3_wp      !< [-] relative tolerance of every group
      real(wp)    :: atol_scale          = 1.0_wp         !< [-] multiplier on every group's atol
      integer(ik) :: error_level         = CTRL_L1_ADAPTIVE !< CTRL_L0_FIXED | CTRL_L1_ADAPTIVE | CTRL_L2_STRICT
      integer(ik) :: step_controller     = CTRL_I         !< CTRL_I | CTRL_PI
      logical     :: ark_coupled         = .true.         !< the coupled leaf <-> canopy-air Newton solve
      real(wp)    :: forcing_sample_frac = 0.5_wp         !< [-] where in the sub-step the forcing is sampled
      logical     :: cas_condensation    = .true.         !< condense canopy-air supersaturation
      logical     :: canopy_water_on     = .false.        !< the canopy interception films
      real(wp)    :: snow_init_swe       = 0.0_wp         !< [kg/m2] snow seeded at the run start
      real(wp)    :: snow_init_temp      = 270.0_wp       !< [K] its temperature
      logical     :: fast_probe          = .false.        !< the sub-daily probe CSV
      character(len=256) :: fast_probe_file = 'fast_probe.csv'
      type(process_mask_t) :: mask                        !< the stores that evolve (all on: the full column)
      integer(ik) :: n_fast_per_slow                      !< DERIVED nint(dt_slow / dt_fast), at least 1

      !----- [soil]: the ground; its solver settings come from [fast]. ------------------------!
      type(soil_column_config_t) :: soil_column           !< its layers and texture
      type(soil_opts_t)          :: soil                  !< its water: surface, bottom, the network's atol
      type(energy_opts_t)        :: energy                !< its heat: the thermal bottom, the budget check
      !----- [snow], [aerodynamics]. -----------------------------------------------------------!
      type(snow_params_t) :: snow
      type(aero_cfg_t)    :: aero

      !----- [demography]. ---------------------------------------------------------------------!
      logical     :: demography_on                        !< .false.: the stand's structure is frozen
      logical     :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse
      integer(ik) :: max_cohort, n_cohort_fusion_iter
      real(wp)    :: cohort_size_tol_min, cohort_size_tol_max
      real(wp)    :: cohort_lai_cap                       !< [m2/m2]
      real(wp)    :: min_cohort_agb                       !< [kgC/m2]
      real(wp)    :: negligible_nplant                    !< [plant/m2]
      real(wp)    :: split_eps
      logical     :: enable_cohort_fission
      integer(ik) :: n_height_layers, max_patch, n_patch_fusion_iter
      real(wp)    :: patch_light_tol, patch_light_tol_max
      real(wp)    :: light_ext_coeff_fusefiss = 0.5_wp    !< [-]
      real(wp)    :: patch_light_maxdev_factor
      real(wp)    :: min_patch_area, conservation_tol
      real(wp)    :: growth_memory_days                   !< [day]
      real(wp)    :: cohort_size_tol_mult                 !< DERIVED geometric step of the cohort tolerance
      real(wp)    :: patch_light_tol_mult                 !< DERIVED geometric step of the patch tolerance
      real(wp), allocatable :: height_edges(:)            !< [m] DERIVED interior edges of the height layers
      !----- [disturbance], [recruitment]. ------------------------------------------------------!
      real(wp) :: patch_disturbance_rate                  !< [1/yr]
      real(wp) :: disturbance_survive_height              !< [m]
      real(wp) :: min_recruit_size                        !< [plant/m2]

      !----- [init], [state], [output]. ---------------------------------------------------------!
      type(init_opts_t)     :: init
      type(state_opts_t)    :: state
      type(output_config_t) :: output

      !----- [carbon], [trait_dynamics]. --------------------------------------------------------!
      real(wp) :: gpp_ref = 0.0_wp                        !< [kgC/m2 leaf/yr] stub GPP, fast loop off
      logical  :: trait_plasticity_on       = .false.     !< light acclimation of the leaf traits
      logical  :: fineroot_turnover_plastic = .false.     !< fine roots turn over with their leaves

      !----- [leaf_physiology]. ------------------------------------------------------------------!
      integer(ik) :: stomatal_model                       !< SM_LEUNING | SM_MEDLYN | SM_KATUL
      integer(ik) :: temp_response_form                   !< TEMP_RESPONSE_ARRHENIUS | TEMP_RESPONSE_PEAKED
      integer(ik) :: colimitation                         !< COLIM_MIN | COLIM_QUADRATIC
      logical     :: use_boundary_layer                   !< couple through the leaf boundary layer
      real(wp)    :: kc25, ko25, gstar25                  !< [Pa]
      real(wp)    :: ea_kc, ea_ko, ea_gstar               !< [J/mol]
      real(wp)    :: ea_vcmax, ea_jmax, ea_rd             !< [J/mol]
      real(wp)    :: hd_vcmax, hd_jmax, hd_rd             !< [J/mol]
      real(wp)    :: ds_vcmax, ds_jmax, ds_rd             !< [J/mol/K]
      real(wp)    :: o2_mol_frac                          !< [mol/mol]
      real(wp)    :: leaf_absorptance                     !< [-]
      real(wp)    :: phi_psii                             !< [-]
      real(wp)    :: medlyn_vpd_min = 0.05_wp             !< [kPa]
      logical     :: thermal_acclimation = .false.        !< Kattge & Knorr (2007), peaked form only
      real(wp)    :: acclim_ds_vcmax_a = 668.39_wp, acclim_ds_vcmax_b = 1.07_wp   !< [J/mol/K], [J/mol/K2]
      real(wp)    :: acclim_ds_jmax_a  = 659.70_wp, acclim_ds_jmax_b  = 0.75_wp   !< [J/mol/K], [J/mol/K2]
      real(wp)    :: acclim_jv_a       = 2.59_wp,   acclim_jv_b       = 0.035_wp  !< [-], [1/K]
      real(wp)    :: acclim_window_days = 30.0_wp      !< [day]

      !----- [soil_carbon], [site]/[forcing], [region], [maps], [[sites]]. ----------------------!
      type(decomp_opts_t)    :: soil_carbon
      type(forcing_config_t) :: forcing
      type(region_opts_t)    :: region
      type(maps_opts_t)      :: maps
      type(site_spec_t), allocatable :: sites(:)          !< not allocated: one site, the polygon

      !----- The PFT file. ------------------------------------------------------------------------!
      type(pft_table_t) :: pft
   end type meds_config_t

contains

   !---------------------------------------------------------------------------------------!
   ! The derived run values: the step in years, the fast sub-steps per slow step, the       !
   ! geometric steps of the fusion tolerances, and the height-layer edges, evenly spaced     !
   ! from 0 to the tallest PFT's hgt_max.                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine derive_config(cfg)
      type(meds_config_t), intent(inout) :: cfg
      integer(ik) :: i
      cfg%dt_years = cfg%dt_slow / yr_sec
      if (cfg%dt_fast > 0.0_wp) then
         cfg%n_fast_per_slow = max(1_ik, nint(cfg%dt_slow / cfg%dt_fast, ik))
      else
         cfg%n_fast_per_slow = 1_ik
      end if
      if (cfg%n_cohort_fusion_iter > 1_ik) then
         cfg%cohort_size_tol_mult = (cfg%cohort_size_tol_max / cfg%cohort_size_tol_min)             &
                               ** (1.0_wp / real(cfg%n_cohort_fusion_iter - 1_ik, wp))
      else
         cfg%cohort_size_tol_mult = 1.0_wp
      end if
      if (cfg%n_patch_fusion_iter > 1_ik .and. cfg%patch_light_tol > 0.0_wp) then
         cfg%patch_light_tol_mult = (cfg%patch_light_tol_max / cfg%patch_light_tol)                  &
                               ** (1.0_wp / real(cfg%n_patch_fusion_iter - 1_ik, wp))
      else
         cfg%patch_light_tol_mult = 1.0_wp
      end if
      if (allocated(cfg%height_edges)) deallocate(cfg%height_edges)
      allocate(cfg%height_edges(cfg%n_height_layers - 1_ik))
      do i = 1_ik, cfg%n_height_layers - 1_ik
         cfg%height_edges(i) = real(i, wp) * maxval(cfg%pft%allom(1:cfg%pft%n)%hgt_max)          &
                               / real(cfg%n_height_layers, wp)
      end do
   end subroutine derive_config

   !----- Every derived value, in one place for the loader and the test builder: the run's     !
   !      (derive_config) and the leaf capacities from their ratios to Vcmax25. Callers validate  !
   !      after this returns. -----------------------------------------------------------------!
   subroutine derive_parameters(cfg)
      type(meds_config_t), intent(inout) :: cfg
      call derive_config(cfg)
      call derive_leaf_params(cfg%pft)
   end subroutine derive_parameters

   !----- The steps the growth-memory window spans (at least 1): the length of each cohort's   !
   !      moving-average ring. --------------------------------------------------------------!
   pure integer(ik) function growth_window_steps(cfg) result(nw)
      type(meds_config_t), intent(in) :: cfg
      nw = max(1_ik, nint(cfg%growth_memory_days / (cfg%dt_years * yr_day), ik))
   end function growth_window_steps

   !----- A PFT's leaf turgor-loss point [MPa], from its own leaf pressure-volume curve, so the      !
   !      stomatal closure and the phenology's drought counter see the curve the solver stores leaf  !
   !      water on. -----------------------------------------------------------------------------------!
   pure real(wp) function pft_leaf_psi_tlp(cfg, ipft) result(psi_tlp)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      associate (leaf_curve => cfg%pft%hydro(ipft)%leaf_curve)
         psi_tlp = tissue_psi_tlp(leaf_curve%pi0, leaf_curve%elastic_mod)
      end associate
   end function pft_leaf_psi_tlp

   !----- Where a PFT's stomatal water stress begins [MPa of predawn leaf potential]: the [pft]       !
   !      stomata_psi_onset where given, else half the PFT's turgor-loss point, so the gravity head  !
   !      of a tall tree in wet soil does not read as drought (docs/science/leaf_gas_exchange.md     !
   !      eq. (10)). ----------------------------------------------------------------------------------!
   pure real(wp) function pft_stomata_psi_onset(cfg, ipft) result(psi_onset)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_onset = 0.5_wp * pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%stomata_psi_onset)) then
         if (cfg%pft%stomata_psi_onset(ipft) > HYD_UNSET) psi_onset = cfg%pft%stomata_psi_onset(ipft)
      end if
   end function pft_stomata_psi_onset

   !----- The leaf potential [MPa] at which the fast water-stress limb on Vcmax/Jmax/TPU starts and   !
   !      completes: the [pft] wstress_psi_open / wstress_psi_close where given, else the PFT's      !
   !      turgor-loss point and twice it (docs/science/leaf_gas_exchange.md eq. (9)). ----------------!
   pure real(wp) function pft_wstress_psi_open(cfg, ipft) result(psi_open)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_open = pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%wstress_psi_open)) then
         if (cfg%pft%wstress_psi_open(ipft) > HYD_UNSET) psi_open = cfg%pft%wstress_psi_open(ipft)
      end if
   end function pft_wstress_psi_open

   pure real(wp) function pft_wstress_psi_close(cfg, ipft) result(psi_close)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      psi_close = 2.0_wp * pft_leaf_psi_tlp(cfg, ipft)
      if (allocated(cfg%pft%wstress_psi_close)) then
         if (cfg%pft%wstress_psi_close(ipft) > HYD_UNSET) psi_close = cfg%pft%wstress_psi_close(ipft)
      end if
   end function pft_wstress_psi_close

   !----- Site `is`'s configuration: the polygon's, with the site's ground when its table gives  !
   !      one. make_site_surface builds the site's surface from it. -----------------------------!
   function site_config(cfg, is) result(scfg)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: is
      type(meds_config_t)             :: scfg
      scfg = cfg
      if (.not. cfg%sites(is)%own_ground) return
      scfg%soil_column            = cfg%sites(is)%soil_column
      scfg%soil%water_bottom_bc   = cfg%sites(is)%water_bottom_bc
      scfg%soil%drain_factor      = cfg%sites(is)%drain_factor
      scfg%energy%deep_temp       = cfg%sites(is)%deep_temp
      scfg%energy%deep_depth      = cfg%sites(is)%deep_depth
   end function site_config

   !=======================================================================================!
   !  The PFT table.                                                                       !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Allocate every trait array to n PFTs.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_pft_table(pft, n)
      type(pft_table_t), intent(inout) :: pft
      integer(ik),       intent(in)    :: n
      pft%n = n
      allocate(pft%allom(n), pft%hydro(n), pft%wood_density(n))
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
      allocate(pft%sla(n), pft%root_to_leaf_ratio(n),                                        &
               pft%storage_cushion(n), pft%growth_resp_factor(n),                            &
               pft%storage_turnover_rate(n), pft%retained_carbon_fraction(n),               &
               pft%is_woody(n), pft%stem_resp_factor25(n), pft%root_resp_factor25(n),        &
               pft%leaf_reflect_vis(n), pft%leaf_transmit_vis(n), pft%leaf_reflect_nir(n),   &
               pft%leaf_transmit_nir(n), pft%leaf_emissivity(n),                             &
               pft%wood_reflect_vis(n), pft%wood_transmit_vis(n), pft%wood_reflect_nir(n),   &
               pft%wood_transmit_nir(n), pft%wood_emissivity(n),                             &
               pft%leaf_clumping(n), pft%wood_clumping(n),                                   &
               pft%leaf_angle_mean(n), pft%leaf_angle_std(n),                                &
               pft%leaf_turnover_rate_toc(n), pft%fineroot_turnover_rate(n))
      pft%storage_turnover_rate = 0.0_wp   ! optional key; 0 = no storage maintenance
      allocate(pft%max_relative_growth_rate(n), pft%max_absolute_growth_rate(n),                     &
               pft%max_absolute_growth_exponent(n))
      pft%max_relative_growth_rate     = 0.0_wp   ! optional keys; a rate of 0 = that sink limit off
      pft%max_absolute_growth_rate     = 0.0_wp
      pft%max_absolute_growth_exponent = 0.0_wp
      pft%retained_carbon_fraction = 0.3_wp ! optional key; 30% of leaf and root loss resorbed
      !----- HYD_UNSET marks a stomatal stress point the PFT file leaves out; the PFT's    !
      !      turgor-loss point then sets it (pft_stomata_psi_onset,                        !
      !      pft_wstress_psi_open/close). -------------------------------------------------!
      pft%stomata_psi_onset = HYD_UNSET   ! half the PFT's turgor-loss point
      pft%wstress_psi_open  = HYD_UNSET   ! the turgor-loss point
      pft%wstress_psi_close = HYD_UNSET   ! twice it
      pft%leaf_surf_water_max = 0.1_wp
      pft%wood_surf_water_max = 0.1_wp   ! optional keys: CLM5's dewmx
      allocate(pft%f_labile_leaf(n), pft%f_labile_stem(n), pft%struct_lignin_frac(n))
      allocate(pft%kplastic_sla(n), pft%kplastic_vcmax(n), pft%kplastic_rd(n), pft%kplastic_leaf_turnover(n))
      allocate(pft%pheno(n))                ! a PFT without a [phenology] block keeps the defaults
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Derive the per-PFT leaf-photosynthesis capacities (Jmax25, TPU25, Rd25) as fixed      !
   ! ratios of Vcmax25. Call after vcmax25 and the ratio arrays are set from the PFT file. !
   !---------------------------------------------------------------------------------------!
   subroutine derive_leaf_params(pft)
      type(pft_table_t), intent(inout) :: pft
      pft%jmax25 = pft%jmax_vcmax_ratio * pft%vcmax25
      pft%tpu25  = pft%tpu_vcmax_ratio  * pft%vcmax25
      pft%rd25   = pft%rd_vcmax_ratio   * pft%vcmax25
   end subroutine derive_leaf_params

   !=======================================================================================!
   !  The soil column: its layers, its texture per layer, and the bundles the soil         !
   !  kernels read, built once per column (each site, each map cell).                      !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! The column's layers: the bottoms [m, positive down] of `layer_bottoms` when the       !
   ! config sets it, else SOIL_LAYER_TEMPLATE, cut at `depth` (soil_layer_grid).           !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_layers(sc, bottoms, n_layer)
      type(soil_column_config_t), intent(in)  :: sc
      real(wp),                   intent(out) :: bottoms(n_soil_layer_max)
      integer(ik),                intent(out) :: n_layer
      if (sc%n_layer_bottoms > 0_ik) then
         call soil_layer_grid(sc%depth, sc%layer_bottoms, sc%n_layer_bottoms, bottoms, n_layer)
      else
         call soil_layer_grid(sc%depth, SOIL_LAYER_TEMPLATE, size(SOIL_LAYER_TEMPLATE, kind=ik),   &
                              bottoms, n_layer)
      end if
   end subroutine soil_column_layers

   !----- The texture per layer: the maps' per-layer values, or the one texture. ----------!
   pure subroutine soil_column_texture(sc, tex)
      type(soil_column_config_t),  intent(in)  :: sc
      type(soil_texture_layers_t), intent(out) :: tex
      if (sc%n_layer_values > 0_ik) then
         tex%theta_sat = sc%layer_theta_sat
         tex%theta_res = sc%layer_theta_res
         tex%ksat = sc%layer_ksat
         tex%par_a     = sc%layer_curve_par_a
         tex%par_n   = sc%layer_curve_par_n
         tex%k_solid   = sc%layer_solid_conductivity
         tex%k_dry = sc%layer_dry_conductivity
         tex%dry_cvol  = sc%layer_dry_heat_capacity
      else
         tex%theta_sat = sc%theta_sat
         tex%theta_res = sc%theta_res
         tex%ksat = sc%ksat
         tex%par_a     = sc%curve_par_a
         tex%par_n   = sc%curve_par_n
         tex%k_solid   = sc%solid_conductivity
         tex%k_dry = sc%dry_conductivity
         tex%dry_cvol  = sc%dry_heat_capacity
      end if
   end subroutine soil_column_texture

   !---------------------------------------------------------------------------------------!
   ! The soil column [soil] describes: its layers (soil_column_layers),                    !
   ! each layer's retention curve and Ks, and, when asked, each layer's thermal texture.   !
   ! The fast loop, the slow ledger and validate_config all build the column here.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_params(sc, soil, therm)
      type(soil_column_config_t),            intent(in)  :: sc
      type(soil_params_t),                   intent(out) :: soil
      type(soil_thermal_params_t), optional, intent(out) :: therm
      type(soil_texture_layers_t) :: tex
      real(wp)    :: bottoms(n_soil_layer_max)
      integer(ik) :: n_layer
      call soil_column_layers(sc, bottoms, n_layer)
      call soil_column_texture(sc, tex)
      call build_soil_hydr_params_layers(n_layer, sc%retention, bottoms, tex%theta_sat, tex%theta_res,  &
                                         tex%ksat, tex%par_a, tex%par_n, sc%psi_fc, soil)
      if (present(therm)) call build_soil_therm_params_layers(n_layer, tex%k_solid, tex%k_dry,        &
                                                              tex%dry_cvol, therm)
   end subroutine soil_column_params

   !---------------------------------------------------------------------------------------!
   ! The water each layer starts the run with [m3/m3]: [init].soil_theta, or each layer's  !
   ! field capacity for a column with a texture per layer (the soil maps).                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_column_start_theta(sc, soil, init_soil_theta, theta)
      type(soil_column_config_t), intent(in)  :: sc
      type(soil_params_t),        intent(in)  :: soil
      real(wp),                   intent(in)  :: init_soil_theta   !< [m3/m3] [init].soil_theta
      real(wp),                   intent(out) :: theta(n_soil_layer_max)
      if (sc%n_layer_values > 0_ik) then
         theta = soil%theta_fc
      else
         theta = init_soil_theta
      end if
   end subroutine soil_column_start_theta

   !---------------------------------------------------------------------------------------!
   ! The water potential [MPa] of the column's top layer at water content theta_top: the   !
   ! predawn potential of a plant that has relaxed to the soil overnight. It seeds the     !
   ! stomatal predawn potential and a new cohort's tissue water.                           !
   !---------------------------------------------------------------------------------------!
   pure real(wp) function top_soil_psi(soil, theta_top) result(psi)
      type(soil_params_t), intent(in) :: soil
      real(wp),            intent(in) :: theta_top   !< [m3/m3]
      psi = soil_matric_potential(soil%retention, theta_top, soil%theta_sat(1), soil%theta_res(1),   &
                                  soil%curve_par_a(1), soil%curve_par_n(1))
   end function top_soil_psi

   !=======================================================================================!
   !  Per-column soil parameter assemblers: pure, state-free constructors that fill a      !
   !  parameter struct from texture and geometry; not part of the fast loop.               !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! A column's layer bottoms [m, positive down]: `list` (SOIL_LAYER_TEMPLATE or a         !
   ! config list) cut at `depth`. The last layer ends at the depth; a remainder thinner    !
   ! than LAYER_MERGE_SHARE of the list layer it falls in joins the layer above. A list    !
   ! that ends above the depth gets one more layer down to it (the config refuses one).    !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_layer_grid(depth, list, n_list, bottoms, n_layer)
      real(wp),    intent(in)  :: depth, list(:)
      integer(ik), intent(in)  :: n_list
      real(wp),    intent(out) :: bottoms(n_soil_layer_max)
      integer(ik), intent(out) :: n_layer
      integer(ik) :: k
      real(wp)    :: top
      bottoms = 0.0_wp
      n_layer = 0_ik
      top     = 0.0_wp
      do k = 1_ik, n_list
         if (list(k) >= depth - GRID_TOL) exit                 ! list layer k holds the cut
         n_layer = n_layer + 1_ik
         bottoms(n_layer) = list(k)
         top = list(k)
      end do
      if (k <= n_list .and. n_layer > 0_ik) then
         if (depth - top < LAYER_MERGE_SHARE * (list(k) - top)) n_layer = n_layer - 1_ik
      end if
      n_layer = n_layer + 1_ik
      bottoms(n_layer) = depth
   end subroutine soil_layer_grid

   !---------------------------------------------------------------------------------------!
   ! Assemble a per-column soil_params_t: the negative-z geometry from the layer bottoms   !
   ! (soil_layer_grid), each layer's own texture, and the derived theta_fc (from each      !
   ! layer's retention curve). The arrays hold one value per active layer: from a soil     !
   ! map's layers, or one value broadcast (build_soil_hydr_params). The roots are the      !
   ! cohorts' own (meds_plant_hydraulics::cohort_root_profile), not the column's.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_hydr_params_layers(n_active, retention, layer_bottoms, theta_sat,       &
                                                 theta_res, ksat, par_a, par_n, psi_fc_m, params)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: layer_bottoms(:)    !< [m] bottoms of the active layers, positive down
      real(wp),            intent(in)  :: theta_sat(:), theta_res(:), ksat(:), par_a(:), par_n(:)
      real(wp),            intent(in)  :: psi_fc_m
      type(soil_params_t), intent(out) :: params
      integer(ik) :: k

      params%n_active  = n_active
      params%retention = retention

      !----- Interface elevations soil_layer_z(k) <= 0, k = 1 at the surface. -------------!
      params%soil_layer_z(1)              = 0.0_wp
      params%soil_layer_z(2:n_active + 1) = -layer_bottoms(1:n_active)

      !----- Thicknesses, node elevations, internode spacings (dz, dz_node > 0). ----------!
      do k = 1_ik, n_active
         params%dz(k)     = params%soil_layer_z(k) - params%soil_layer_z(k+1)
         params%z_node(k) = 0.5_wp * (params%soil_layer_z(k) + params%soil_layer_z(k+1))
      end do
      do k = 1_ik, n_active - 1_ik
         params%dz_node(k) = params%z_node(k) - params%z_node(k+1)
      end do
      params%dz_node(n_active) = params%dz(n_active)      ! unused at the bottom face; kept > 0

      !----- Each layer's texture. --------------------------------------------------------!
      params%theta_sat(1:n_active) = theta_sat(1:n_active)
      params%theta_res(1:n_active) = theta_res(1:n_active)
      params%ksat(1:n_active)      = ksat(1:n_active)
      params%curve_par_a(1:n_active) = par_a(1:n_active)
      params%curve_par_n(1:n_active) = par_n(1:n_active)

      !----- Derived threshold (from each layer's retention curve). -----------------------!
      params%theta_fc(1:n_active) = soil_head_to_theta(retention, psi_fc_m, theta_sat(1:n_active),    &
                                        theta_res(1:n_active), par_a(1:n_active), par_n(1:n_active))
   end subroutine build_soil_hydr_params_layers

   !----- One texture for every layer: the per-layer build with each value broadcast. -----!
   pure subroutine build_soil_hydr_params(n_active, retention, layer_bottoms, theta_sat, theta_res,  &
                                ksat, par_a, par_n, psi_fc_m, params)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: layer_bottoms(:)    !< [m] bottoms of the active layers, positive down
      real(wp),            intent(in)  :: theta_sat, theta_res, ksat, par_a, par_n, psi_fc_m
      type(soil_params_t), intent(out) :: params
      real(wp) :: ts(n_active), tr(n_active), ks(n_active), pa(n_active), pn(n_active)
      ts = theta_sat
      tr = theta_res
      ks = ksat
      pa = par_a
      pn = par_n
      call build_soil_hydr_params_layers(n_active, retention, layer_bottoms, ts, tr, ks, pa, pn,        &
                                         psi_fc_m, params)
   end subroutine build_soil_hydr_params

   !----- The column's thermal texture, one value per active layer. -----------------------!
   pure subroutine build_soil_therm_params_layers(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid(:), k_dry(:), dry_cvol(:)
      type(soil_thermal_params_t), intent(out) :: therm
      therm%nzg_active = n_active
      therm%soil_solid_conductivity(1:n_active) = k_solid(1:n_active)
      therm%soil_dry_conductivity(1:n_active)   = k_dry(1:n_active)
      therm%soil_dry_heat_capacity(1:n_active)  = dry_cvol(1:n_active)
   end subroutine build_soil_therm_params_layers

   !----- One thermal texture for every layer. --------------------------------------------!
   pure subroutine build_soil_therm_params(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid, k_dry, dry_cvol
      type(soil_thermal_params_t), intent(out) :: therm
      real(wp) :: ks(n_active), kd(n_active), cv(n_active)
      ks = k_solid
      kd = k_dry
      cv = dry_cvol
      call build_soil_therm_params_layers(n_active, ks, kd, cv, therm)
   end subroutine build_soil_therm_params

end module meds_config
