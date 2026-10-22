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
! The routines derive what the parameters imply: the run's derived values, and through    !
! meds_pft_derive the PFT table's leaf capacities and the tables the fast loop reads;      !
! apply_map_inputs writes a map cell's inputs into its polygon's configuration. A site's   !
! soil profile is meds_soil_derive's.                                                      !
!==========================================================================================!
module meds_config
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : yr_day, yr_sec, n_soil_layer_max
   use meds_time,               only : meds_time_t
   use meds_pft_params,         only : pft_table_t, HYD_UNSET, leaf_photo_pft_t
   use meds_pft_derive,         only : derive_leaf_params, derive_pft_tables
   use meds_soil_hydr_lib,      only : soil_texture_t, soil_texture_to_hydraulics, SOIL_RETENTION_VG
   use meds_therm_lib,          only : soil_texture_to_thermal
   use meds_fast_params,        only : INTEG_ARK, CTRL_I
   use meds_soil_params,        only : soil_params_t, soil_profile_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER, &
                                       ANNUAL_DAMPING_DEPTH
   use meds_soil_derive,        only : soil_layers, build_soil_profile, map_layer_texture
   use meds_biophysics_params,  only : snow_params_t, aero_params_t
   use meds_biogeochem_params,  only : soil_carbon_params_t
   use meds_forcing_params,     only : forcing_params_t
   use meds_output_params,      only : output_params_t
   use meds_region_params,      only : region_params_t, RUN_MODE_SITE, RUN_MODE_REGION
   use meds_maps_params,        only : maps_params_t, soil_maps_on, map_inputs_t, MAP_BC_LEAKY,        &
                                       MAP_BC_SEALED, MAP_BC_AQUIFER
   use meds_sites_params,       only : site_spec_t
   implicit none
   private

   public :: meds_config_t, init_params_t
   public :: INIT_BARE, INIT_CENSUS, INIT_RESTART
   public :: derive_config, derive_parameters, growth_window_steps, site_config
   public :: forcing_params_t, output_params_t
   public :: soil_carbon_params_t
   public :: region_params_t, RUN_MODE_SITE, RUN_MODE_REGION
   public :: maps_params_t, soil_maps_on
   public :: alloc_pft_table, apply_map_inputs

   !----- [init].init_mode. -------------------------------------------------------------------!
   integer(ik), parameter :: INIT_BARE    = 0_ik    !< near-bare ground
   integer(ik), parameter :: INIT_CENSUS  = 1_ik    !< a cohort census (census_file)
   integer(ik), parameter :: INIT_RESTART = 2_ik    !< a state file (restart_file)

   !----- [init]: how the run starts. -----------------------------------------------------------!
   type :: init_params_t
      integer(ik)        :: init_mode                     !< INIT_BARE | INIT_CENSUS | INIT_RESTART
      character(len=256) :: restart_file = ''             !< the state file init_mode = 2 reads
      character(len=256) :: census_file  = ''             !< the cohort CSV init_mode = 1 reads
      character(len=256) :: pft_config                    !< the PFT file
      real(wp)           :: soil_temp    = 288.0_wp       !< [K] every soil layer's start
      real(wp)           :: soil_theta   = 0.30_wp        !< [m3/m3] every soil layer's start
      logical            :: reacclimate_traits = .false.  !< a restart's plastic traits from this PFT file
   end type init_params_t

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
      real(wp)    :: rtol                = 1.0e-3_wp      !< [-] relative tolerance of every group
      real(wp)    :: atol_scale          = 1.0_wp         !< [-] multiplier on every group's atol
      integer(ik) :: step_controller     = CTRL_I         !< CTRL_I | CTRL_PI
      real(wp)    :: water_atol          = 1.0e-4_wp      !< [-] the water network's absolute tolerance
      real(wp)    :: soil_energy_budget_atol = 1.0e-2_wp  !< [K] the soil heat's budget stop; not a step tolerance
      logical     :: debug_error         = .false.        !< a budget breach stops the run; else it is reported
      real(wp)    :: forcing_sample_frac = 0.5_wp         !< [-] where in the sub-step the forcing is sampled
      logical     :: canopy_water_on     = .false.        !< the canopy interception films
      real(wp)    :: intercept_alpha     = 1.0_wp         !< [-] Beer interception efficiency
      real(wp)    :: intercept_k         = 0.5_wp         !< [1/PAI] Beer interception extinction
      real(wp)    :: snow_init_swe       = 0.0_wp         !< [kg/m2] snow seeded at the run start
      real(wp)    :: snow_init_temp      = 270.0_wp       !< [K] its temperature
      integer(ik) :: n_fast_per_slow                      !< DERIVED nint(dt_slow / dt_fast), at least 1

      !----- [soil]: the ground: its column, surface and bottom boundaries. ---------------------!
      type(soil_params_t) :: soil
      !----- [snow], [aerodynamics]. -----------------------------------------------------------!
      type(snow_params_t) :: snow
      type(aero_params_t) :: aero

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

      !----- [init], [output] (the diagnostics and the state files). -----------------------------!
      type(init_params_t)   :: init
      type(output_params_t) :: output

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
      type(soil_carbon_params_t) :: soil_carbon
      type(forcing_params_t)     :: forcing
      type(region_params_t)      :: region
      type(maps_params_t)        :: maps
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
   !      (derive_config), the leaf capacities from their ratios to Vcmax25, and the PFT tables   !
   !      the fast loop reads. A caller that changes a PFT trait or a [leaf_physiology] key      !
   !      afterwards calls it again. Callers validate after this returns. ----------------------!
   subroutine derive_parameters(cfg)
      type(meds_config_t), intent(inout) :: cfg
      call derive_config(cfg)
      call derive_leaf_params(cfg%pft)
      call derive_pft_tables(cfg%pft, leaf_photo_shared(cfg), cfg%stomatal_model, cfg%temp_response_form,   &
                             cfg%colimitation, cfg%use_boundary_layer)
   end subroutine derive_parameters

   !----- [leaf_physiology]'s constants: the part of every PFT's leaf parameters the PFTs share. --!
   pure function leaf_photo_shared(cfg) result(p)
      type(meds_config_t), intent(in) :: cfg
      type(leaf_photo_pft_t)          :: p
      p%kc25 = cfg%kc25
      p%ko25 = cfg%ko25
      p%gstar25 = cfg%gstar25
      p%ea_kc = cfg%ea_kc
      p%ea_ko = cfg%ea_ko
      p%ea_gstar = cfg%ea_gstar
      p%ea_vcmax = cfg%ea_vcmax
      p%ea_jmax = cfg%ea_jmax
      p%ea_rd = cfg%ea_rd
      p%hd_vcmax = cfg%hd_vcmax
      p%hd_jmax = cfg%hd_jmax
      p%hd_rd = cfg%hd_rd
      p%ds_vcmax = cfg%ds_vcmax
      p%ds_jmax = cfg%ds_jmax
      p%ds_rd = cfg%ds_rd
      p%o2_mol_frac = cfg%o2_mol_frac
      p%absorptance = cfg%leaf_absorptance
      p%phi_psii = cfg%phi_psii
      p%medlyn_vpd_min = cfg%medlyn_vpd_min
   end function leaf_photo_shared

   !----- The steps the growth-memory window spans (at least 1): the length of each cohort's   !
   !      moving-average ring. --------------------------------------------------------------!
   pure integer(ik) function growth_window_steps(cfg) result(nw)
      type(meds_config_t), intent(in) :: cfg
      nw = max(1_ik, nint(cfg%growth_memory_days / (cfg%dt_years * yr_day), ik))
   end function growth_window_steps

   !----- Site `is`'s configuration: the polygon's, with the site's ground when its table gives  !
   !      one. make_site_surface builds the site's surface from it. -----------------------------!
   function site_config(cfg, is) result(scfg)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: is
      type(meds_config_t)             :: scfg
      scfg = cfg
      if (cfg%sites(is)%own_ground) scfg%soil = cfg%sites(is)%soil
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
   ! apply_map_inputs -- cell k's inputs into the configuration its polygon runs           !
   ! (docs/science/soil_biophysics.md §The column from the soil maps): the soil column on  !
   ! the default layers to the map's depth, each layer's curve, Ks and thermal properties  !
   ! from the map layers it overlaps; the bottom water boundary from bc_class; the thermal !
   ! anchor at deep_temp, ANNUAL_DAMPING_DEPTH / sqrt(2) below the bottom node (eq. (13)); !
   ! and the Chave allometry's E.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine apply_map_inputs(inp, k, cfg)
      type(map_inputs_t), intent(in)    :: inp
      integer(ik),            intent(in)    :: k
      type(meds_config_t),    intent(inout) :: cfg
      type(soil_texture_t) :: tex
      type(soil_profile_t) :: profile
      real(wp)    :: bottoms(n_soil_layer_max), top
      integer(ik) :: n, l

      if (inp%soil) then
         associate (sc => cfg%soil)
            sc%depth = inp%soil_depth(k)
            sc%n_layer_bottoms = 0_ik
            sc%retention = SOIL_RETENTION_VG
            call soil_layers(sc, bottoms, n)
            top = 0.0_wp
            do l = 1_ik, n
               tex = map_layer_texture(inp, k, top, bottoms(l))
               call soil_texture_to_hydraulics(tex, sc%layer_theta_sat(l), sc%layer_theta_res(l),      &
                                                 sc%layer_curve_par_a(l), sc%layer_curve_par_n(l),       &
                                                 sc%layer_ksat(l))
               call soil_texture_to_thermal(tex%sand, tex%clay, sc%layer_theta_sat(l),                     &
                                              sc%layer_solid_conductivity(l), sc%layer_dry_conductivity(l),  &
                                              sc%layer_dry_heat_capacity(l))
               top = bottoms(l)
            end do
            sc%n_layer_values = n
         end associate
         select case (nint(inp%bc_class(k), ik))
         case (MAP_BC_AQUIFER)
            cfg%soil%water_bottom_bc = SOIL_BC_AQUIFER
         case (MAP_BC_LEAKY)
            cfg%soil%water_bottom_bc = SOIL_BC_DRAIN
            cfg%soil%drain_factor = inp%f_leak(k)
         case (MAP_BC_SEALED)
            cfg%soil%water_bottom_bc = SOIL_BC_DRAIN
            cfg%soil%drain_factor = 0.0_wp
         case default
            cfg%soil%water_bottom_bc = SOIL_BC_DRAIN
            cfg%soil%drain_factor = 1.0_wp
         end select
         call build_soil_profile(cfg%soil, profile)
         cfg%soil%deep_temp  = inp%deep_temp(k)
         cfg%soil%deep_depth = abs(profile%z_node(profile%n_active)) + ANNUAL_DAMPING_DEPTH / sqrt(2.0_wp)
      end if
      if (inp%chave) cfg%pft%allom(:)%chave_e = inp%chave_e(k)
   end subroutine apply_map_inputs

end module meds_config
