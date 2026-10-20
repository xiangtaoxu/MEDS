! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config -- the run configuration, threaded read-only through the engine             !
! (docs/configuration.md).                                                                 !
!                                                                                          !
! One field per key, named as the key, in the order of meds_config_main.toml: a key of a   !
! block with its own record is cfg%<record>%<key> (cfg%init%soil_temp, cfg%forcing%format); !
! a key of the other blocks is cfg%<key>. The two exceptions are the [output] group, axis  !
! and tier switches, which are arrays indexed by their codes, and the [hydraulics] tissue  !
! curves (leaf_pi0 is hydro(pft)%leaf_curve%pi0). The reference config describes each key; !
! a field without an initializer is a required key or a derived value (marked DERIVED).    !
!==========================================================================================!
module meds_config
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : yr_day, yr_sec
   use meds_time,               only : meds_time_t
   use meds_pft_params,         only : pft_table_t, derive_leaf_params, HYD_UNSET
   use meds_plant_hydr_lib,     only : tissue_psi_tlp
   use meds_fast_opts,          only : process_mask_t, INTEG_ARK, CTRL_I, CTRL_L1_ADAPTIVE
   use meds_soil_column_config, only : soil_column_config_t
   use meds_biophysics_opts,    only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t
   use meds_biogeochem_opts,    only : decomp_opts_t
   use meds_forcing_config,     only : forcing_config_t
   use meds_output_config,      only : output_config_t
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

   !----- [init].init_mode. -------------------------------------------------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< a cohort census (census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< a state file (restart_file)

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

end module meds_config
