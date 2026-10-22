! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_config -- the fast loop's setup from meds_config_t. build_fast_context builds  !
! the run's fast_context_t: the column configuration (the soil column, the per-PFT        !
! hydraulics, film and leaf-photosynthesis tables, the integrator options) and the canopy !
! optics table. site_column_configs gives each site its copy, with its own soil and its   !
! leaf table, which acclimate_leaf_photo_table refreshes each slow step. Nothing below     !
! this module reads meds_config_t: the kernels and the numerics take these records.       !
!==========================================================================================!
module meds_fast_config
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t, soil_column_params, pft_stomata_psi_onset, pft_wstress_psi_open,  &
                                pft_wstress_psi_close
   use meds_pft_params,  only : pft_table_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_args_t, leaf_photo_table_t, hydro_args_table_t
   use meds_plant_hydr_lib, only : build_hydro_table, tissue_psi_to_rwc
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_therm_lib,     only : kattge_knorr_entropy, kattge_knorr_jv_ratio
   use meds_fast_types,  only : column_config_t, column_site_t, tol_set_t, error_control_t, integrator_opts_t, GRP_THETA
   use meds_canopy_types,     only : rad_pft_optics_t, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use meds_canopy_radiation, only : derive_rad_optics
   use meds_optics_lib,       only : moments_to_beta_params
   implicit none
   private

   public :: fast_context_t, build_fast_context, site_column_configs
   public :: build_leaf_photo_table
   public :: acclimate_leaf_photo_table
   public :: build_rad_optics
   public :: build_integrator_opts
   public :: apply_hydraulics_config, apply_canopy_film_config
   public :: leaf_photo_params_for_pft, leaf_gas_exchange, build_tol_set   ! used only by the tests

   !----- What the fast loop needs beyond the polygon and cfg (build_fast_context): the column      !
   !      configuration, and the constant reference climate a run without a forcing source uses. ---!
   type :: fast_context_t
      !----- The column configuration: the run part every site shares, and the site part [soil] and  !
      !      the PFT table set, from which each site's own is made (site_column_configs). ------------!
      type(column_config_t) :: col_config
      type(column_site_t)   :: col_site
      !----- Each site's leaf table, acclimated to its growth temperature; until then col_site's. ---!
      type(leaf_photo_table_t), allocatable :: site_leaf_photo(:)
      real(wp) :: u_ref   = 2.0_wp                  !< [m/s] reference wind (at the canopy-air top)
      real(wp) :: press   = 101325.0_wp             !< [Pa]
      real(wp) :: rho_air = 1.2_wp                  !< [kg/m3]
      real(wp) :: air_temp = 288.0_wp               !< [K]        reference-level air temperature
      real(wp) :: shv_atm  = 0.008_wp               !< [kg/kg]    reference-level specific humidity
      real(wp) :: co2_atm  = 400.0_wp               !< [umol/mol] free-atmosphere CO2
      real(wp) :: rad_sw_top    = 400.0_wp          !< [W/m2] shortwave into the canopy (leaves)
      real(wp) :: rad_sw_ground = 60.0_wp           !< [W/m2] shortwave reaching the ground
      real(wp) :: precip          = 0.0_wp            !< [kg/m2/s] total precipitation (split like any record)
      !----- Patch-loop threads: [run].n_threads in a one-polygon run, else 1 (polygon_prepare). --!
      integer(ik) :: patch_threads = 1_ik
      real(wp) :: veg_height_bare = 1.0_wp          !< [m] canopy height for a cohort-free patch
      !----- Canopy-RT optics and ground surface, built once from the [pft] table and [soil]. -----!
      type(rad_pft_optics_t) :: rad_opt             !< per-PFT, per-band leaf/wood scattering + LIDF
      real(wp) :: soil_albedo(3) = [0.15_wp, 0.30_wp, 0.0_wp]  !< [VIS,NIR,LW] ground albedo (LW=0)
      real(wp) :: soil_emiss     = 0.95_wp          !< [-] ground longwave emissivity
   end type fast_context_t

contains

   !=======================================================================================!
   !  Build the fast-loop context from the run configuration: the soil column, the per-PFT !
   !  hydraulics, leaf and optics tables, and the solver options. The reference climate    !
   !  keeps the fast_context_t defaults.                                                   !
   !=======================================================================================!
   subroutine build_fast_context(cfg, ctx)
      type(meds_config_t),  intent(in)  :: cfg
      type(fast_context_t), intent(out) :: ctx
      !----- The soil column from [soil]: one texture for the column or one per layer (soil maps). -!
      call soil_column_params(cfg%soil_column, ctx%col_site%soil, ctx%col_site%soil_thermal)
      !----- Per-PFT plant hydraulics from each PFT's [hydraulics] record. ------------------------!
      call apply_hydraulics_config(cfg, ctx%col_config%hydraulics_table)
      call apply_canopy_film_config(cfg%pft, ctx%col_config%leaf_surf_water_max, ctx%col_config%wood_surf_water_max)
      call build_leaf_photo_table(cfg, ctx%col_site%leaf_photo)    ! per-PFT leaf parameters, once per run
      !----- The canopy interception film (+ film evaporation, dew) on/off, [fast] block. --------!
      ctx%col_config%canopy_water_on    = cfg%canopy_water_on
      !----- The [soil], [snow], [aerodynamics] and [fast] options, copied as they are. -----------!
      ctx%col_site%soil_water_opts  = cfg%soil        ! [soil], [fast] -> soil water and the water network
      ctx%soil_albedo(1) = cfg%soil%ground_albedo_vis   ! [soil] ground optics -> the canopy radiation solver
      ctx%soil_albedo(2) = cfg%soil%ground_albedo_nir
      ctx%soil_emiss     = cfg%soil%ground_emissivity
      ctx%col_site%energy = cfg%energy      ! [soil], [fast] -> soil-thermal solver opts
      ctx%col_config%snow   = cfg%snow        ! [snow]         -> snow physical parameter table
      ctx%col_config%aero   = cfg%aero        ! [aerodynamics] -> canopy-aerodynamics constants
      ctx%col_config%soil_carbon = cfg%soil_carbon   ! [soil_carbon] -> heterotrophic respiration
      ctx%col_config%ea_rd = cfg%ea_rd             ! [leaf_physiology] Rd response -> stem and root respiration
      ctx%col_config%hd_rd = cfg%hd_rd
      ctx%col_config%ds_rd = cfg%ds_rd

      !----- One tolerance set for the whole fast loop (build_tol_set); the water network takes the   !
      !      GRP_THETA pair (docs/science/soil_plant_water.md §3). --------------------------------!
      ctx%col_config%integrator = build_integrator_opts(cfg)
      associate (tols => ctx%col_config%integrator%error_control%tols)
         ctx%col_site%soil_water_opts%rtol   = tols%rtol(GRP_THETA)
         ctx%col_site%soil_water_opts%water_atol   = tols%atol(GRP_THETA)
         !----- A budget-closure threshold, not a tolerance group: atol_scale scales it directly. -!
         ctx%col_site%energy%soil_energy_budget_atol  = cfg%energy%soil_energy_budget_atol * cfg%atol_scale
      end associate

      ctx%col_config%mask   = cfg%mask        ! [fast].mask    -> the stores that evolve (all on: the column)

      !----- The canopy-RT optics table, built once and read-only downstream. --------------------!
      call build_rad_optics(cfg, ctx%rad_opt)
   end subroutine build_fast_context

   !----- Each site's part of the column configuration: the run's, with the site's soil column,       !
   !      bottom and acclimated leaf table. Made every slow step, so each has the current table. -----!
   subroutine site_column_configs(ctx, poly, site_col)
      type(fast_context_t),  intent(in)  :: ctx
      type(polygon_t),       intent(in)  :: poly
      type(column_site_t),   intent(out) :: site_col(:)
      integer(ik) :: is
      do is = 1_ik, poly%site%n
         site_col(is) = ctx%col_site
         associate (col => site_col(is), surf => poly%site%surface(is))
            col%soil                         = surf%soil
            col%soil_thermal                 = surf%soil_thermal
            col%soil_water_opts%water_bottom_bc    = surf%water_bottom_bc
            col%soil_water_opts%drain_factor = surf%drain_factor
            col%energy%deep_temp             = surf%deep_temp
            col%energy%deep_depth            = surf%deep_depth
            if (allocated(ctx%site_leaf_photo)) col%leaf_photo = ctx%site_leaf_photo(is)
         end associate
      end do
   end subroutine site_column_configs

   !----- PFT ipft's leaf-photosynthesis parameters, with the table's Vcmax25/Jmax25/TPU25/Rd25: the  !
   !      one place the configuration is flattened into the kernel's record. leaf_gas_exchange       !
   !      applies a cohort's plastic overrides on top of it. --------------------------------------!
   pure subroutine leaf_photo_params_for_pft(cfg, ipft, p)
      type(meds_config_t),       intent(in)  :: cfg
      integer(ik),               intent(in)  :: ipft
      type(leaf_photo_args_t),   intent(out) :: p
      associate (t => cfg%pft)
         p%pathway        = t%photosynthetic_pathway(ipft)
         p%vcmax25        = t%vcmax25(ipft)
         p%jmax25         = t%jmax25(ipft)
         p%tpu25          = t%tpu25(ipft)
         p%rd25           = t%rd25(ipft)
         p%kp25           = t%kp25(ipft)
         p%g0             = t%stomatal_g0(ipft)
         p%g1             = t%stomatal_g1(ipft)
         p%d0             = t%stomatal_d0(ipft)
         p%quantum_yield  = t%quantum_yield_c4(ipft)
         p%theta_j        = t%theta_j(ipft)
         p%theta_cj_c3    = t%theta_cj_c3(ipft)
         p%theta_ip_c3    = t%theta_ip_c3(ipft)
         p%theta_cj_c4    = t%theta_cj_c4(ipft)
         p%theta_ic_c4    = t%theta_ic_c4(ipft)
         p%lambda25       = t%katul_lambda25(ipft)
         p%psi_open       = pft_wstress_psi_open(cfg, ipft)    ! turgor loss unless the PFT file sets it
         p%psi_close      = pft_wstress_psi_close(cfg, ipft)   ! twice it unless the PFT file sets it
         p%lambda_psi_exp = t%wstress_lambda_exp(ipft)
         p%sref_stomata   = t%wstress_sref_stomata(ipft)
         p%psi_onset      = pft_stomata_psi_onset(cfg, ipft)
      end associate
      !----- Copy the shared biochemistry constants. --------------------------------------!
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
   end subroutine leaf_photo_params_for_pft

   !----- Every PFT's leaf parameters plus the run-level solver selectors, for the column config. --!
   subroutine build_leaf_photo_table(cfg, table)
      type(meds_config_t),      intent(in)  :: cfg
      type(leaf_photo_table_t), intent(out) :: table
      integer(ik) :: ipft
      table%n_pft = cfg%pft%n
      allocate(table%pft(table%n_pft), table%jmax_vcmax_ratio(table%n_pft), table%tpu_vcmax_ratio(table%n_pft))
      do ipft = 1_ik, table%n_pft
         call leaf_photo_params_for_pft(cfg, ipft, table%pft(ipft))
         table%jmax_vcmax_ratio(ipft) = cfg%pft%jmax_vcmax_ratio(ipft)
         table%tpu_vcmax_ratio(ipft)  = cfg%pft%tpu_vcmax_ratio(ipft)
      end do
      table%stomatal_model     = cfg%stomatal_model
      table%temp_response_form = cfg%temp_response_form
      table%colimitation       = cfg%colimitation
      table%use_boundary_layer = cfg%use_boundary_layer
   end subroutine build_leaf_photo_table

   !=======================================================================================!
   !  Refresh a site's leaf table, the peaked-form entropy terms and the Jmax:Vcmax ratio, !
   !  from its running-mean growth temperature (Kattge & Knorr 2007;                       !
   !  docs/science/leaf_gas_exchange.md §2, Thermal acclimation, eq. (4)). The reference   !
   !  coefficients stay in cfg, so the refresh is idempotent.                              !
   !=======================================================================================!
   subroutine acclimate_leaf_photo_table(cfg, t_growth, table)
      type(meds_config_t),      intent(in)    :: cfg
      real(wp),                 intent(in)    :: t_growth     !< [K] running-mean growth temperature
      type(leaf_photo_table_t), intent(inout) :: table
      integer(ik) :: ipft
      if (.not. cfg%thermal_acclimation) return
      if (t_growth <= 0.0_wp) return                          ! not yet seeded: leave the fixed values
      do ipft = 1_ik, table%n_pft
         table%pft(ipft)%ds_vcmax = kattge_knorr_entropy(cfg%acclim_ds_vcmax_a, cfg%acclim_ds_vcmax_b, t_growth)
         table%pft(ipft)%ds_jmax  = kattge_knorr_entropy(cfg%acclim_ds_jmax_a,  cfg%acclim_ds_jmax_b,  t_growth)
         table%jmax_vcmax_ratio(ipft) = kattge_knorr_jv_ratio(cfg%acclim_jv_a, cfg%acclim_jv_b, t_growth)
      end do
   end subroutine acclimate_leaf_photo_table

   !----- The canopy two-stream's per-PFT optics table from the [pft] traits                         !
   !      (docs/science/canopy_radiation_transfer.md §1-§2); the canopy C API builds the same. ------!
   subroutine build_rad_optics(cfg, optics)
      type(meds_config_t),    intent(in)  :: cfg
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT
      integer(ik) :: np, ipf
      real(wp), allocatable :: rl(:,:), tl(:,:), rw(:,:), tw(:,:), cl(:), cw(:), bp(:), bq(:)
      logical  :: hb(NB), he(NB)
      real(wp) :: bpp, bqq
      np = cfg%pft%n
      allocate(rl(NB,np), tl(NB,np), rw(NB,np), tw(NB,np), cl(np), cw(np), bp(np), bq(np))
      !----- Shortwave comes as reflectance and transmittance; longwave as emissivity, whose band   !
      !      reflectance is 1 - emissivity with zero transmittance (a leaf is opaque there). ---!
      associate (t => cfg%pft)
         rl(RAD_VIS,1:np) = t%leaf_reflect_vis(1:np)
         tl(RAD_VIS,1:np) = t%leaf_transmit_vis(1:np)
         rl(RAD_NIR,1:np) = t%leaf_reflect_nir(1:np)
         tl(RAD_NIR,1:np) = t%leaf_transmit_nir(1:np)
         rl(RAD_LW ,1:np) = 1.0_wp - t%leaf_emissivity(1:np)
         tl(RAD_LW,1:np) = 0.0_wp
         rw(RAD_VIS,1:np) = t%wood_reflect_vis(1:np)
         tw(RAD_VIS,1:np) = t%wood_transmit_vis(1:np)
         rw(RAD_NIR,1:np) = t%wood_reflect_nir(1:np)
         tw(RAD_NIR,1:np) = t%wood_transmit_nir(1:np)
         rw(RAD_LW ,1:np) = 1.0_wp - t%wood_emissivity(1:np)
         tw(RAD_LW,1:np) = 0.0_wp
         cl(1:np) = t%leaf_clumping(1:np)
         cw(1:np) = t%wood_clumping(1:np)
         !----- The Beta leaf-angle shape parameters, from the angle's measured mean and std. ---!
         do ipf = 1_ik, np
            call moments_to_beta_params(t%leaf_angle_mean(ipf), t%leaf_angle_std(ipf), bpp, bqq)
            bp(ipf) = bpp
            bq(ipf) = bqq
         end do
      end associate
      hb = [.true.,  .true.,  .false.]                         ! VIS/NIR have a beam; LW does not
      he = [.false., .false., .true. ]                         ! only LW emits
      call derive_rad_optics(NB, np, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, optics)
   end subroutine build_rad_optics

   !----- The per-PFT hydraulics table: each PFT's [hydraulics] traits, the Kirchhoff lookup built   !
   !      from its own wood_kexp (docs/science/plant_hydraulics.md eq. (8)), which makes a          !
   !      non-integer vulnerability exponent affordable on the hot path, and the leaf water content   !
   !      at which its stomata close, the water network's supply limit (soil_plant_water.md eq. (3)). !
   subroutine apply_hydraulics_config(cfg, table)
      type(meds_config_t),        intent(in)  :: cfg
      type(hydro_args_table_t),   intent(out) :: table
      integer(ik) :: i
      table%n_pft = cfg%pft%n
      allocate(table%pft(max(cfg%pft%n, 1_ik)))
      do i = 1_ik, cfg%pft%n
         table%pft(i)%plant_hydraulics_t = cfg%pft%hydro(i)
         call build_hydro_table(table%pft(i)%vuln_table, table%pft(i)%wood_kexp)
         table%pft(i)%rwc_close = tissue_psi_to_rwc(pft_wstress_psi_close(cfg, i), table%pft(i)%leaf_curve%pi0,   &
                                                    table%pft(i)%leaf_curve%elastic_mod)
      end do
   end subroutine apply_hydraulics_config

   !----- The canopy films' per-PFT capacities, from the PFT table (both default to 0.1). -------!
   subroutine apply_canopy_film_config(pft, leaf_max, wood_max)
      type(pft_table_t),     intent(in)  :: pft
      real(wp), allocatable, intent(out) :: leaf_max(:), wood_max(:)
      leaf_max = pft%leaf_surf_water_max(1:pft%n)
      wood_max = pft%wood_surf_water_max(1:pft%n)
   end subroutine apply_canopy_film_config

   !----- One leaf's gas exchange for PFT ipft, with a cohort's optional plastic overrides. -------!
   subroutine leaf_gas_exchange(env, cfg, ipft, flux, vcmax25, rd25)
      type(leaf_env_t),    intent(in)  :: env
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ipft
      type(leaf_flux_t),   intent(out) :: flux
      real(wp), optional,  intent(in)  :: vcmax25   !< per-cohort (plastic) Vcmax25 override; jmax/tpu scale with it
      real(wp), optional,  intent(in)  :: rd25      !< per-cohort (plastic) Rd25 override
      type(leaf_photo_args_t)          :: p
      call leaf_photo_params_for_pft(cfg, ipft, p)
      if (present(vcmax25)) then
         p%vcmax25 = vcmax25
         p%jmax25  = cfg%pft%jmax_vcmax_ratio(ipft) * vcmax25
         p%tpu25   = cfg%pft%tpu_vcmax_ratio(ipft)  * vcmax25
      end if
      if (present(rd25)) p%rd25 = rd25
      call solve_leaf_gas_exchange(env, p, cfg%stomatal_model, cfg%temp_response_form,         &
                                   cfg%colimitation, cfg%use_boundary_layer, flux)
   end subroutine leaf_gas_exchange

   !---------------------------------------------------------------------------------------!
   ! The single tolerance set of the fast loop, the per-group tolerances of its error norm !
   ! (docs/science/numerical_scheme.md eq. (3)): [fast].rtol for every group, each group's !
   ! default atol, and [fast].water_atol for GRP_THETA. ARK advances plant water outside   !
   ! its stages, so GRP_LEAF_W and GRP_WOOD_W act only under RK45.                         !
   !---------------------------------------------------------------------------------------!
   pure function build_tol_set(cfg) result(tols)
      type(meds_config_t), intent(in) :: cfg
      type(tol_set_t)                 :: tols
      tols%rtol             = cfg%rtol
      tols%atol(GRP_THETA)  = cfg%soil%water_atol
      !----- Every atol scales with atol_scale: rtol alone stops tightening once atol dominates the   !
      !      WRMS denominator atol + rtol*|y| (docs/science/numerical_scheme.md §6). ----------------!
      tols%atol = tols%atol * cfg%atol_scale
   end function build_tol_set

   !----- The error-control bundle: the tolerance set, the step controller and the strictness. -------!
   pure function build_error_control(cfg) result(ec)
      type(meds_config_t), intent(in) :: cfg
      type(error_control_t)           :: ec
      ec%tols       = build_tol_set(cfg)
      ec%controller = cfg%step_controller
      ec%level      = cfg%error_level
   end function build_error_control

   !----- The integrator's whole configuration, once per run (carried on column_config_t%integrator). -!
   pure function build_integrator_opts(cfg) result(opts)
      type(meds_config_t), intent(in) :: cfg
      type(integrator_opts_t)         :: opts
      opts%scheme           = cfg%time_integrator
      opts%adaptive         = cfg%ark_adaptive
      opts%dt_init          = cfg%ark_dt_init
      opts%coupled_newton   = cfg%ark_coupled
      opts%fixed_substeps   = cfg%ark_fixed_substep
      opts%cas_condensation = cfg%cas_condensation
      opts%error_control    = build_error_control(cfg)
   end function build_integrator_opts

end module meds_fast_config
