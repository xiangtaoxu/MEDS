! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_config -- the fast loop's setup from meds_config_t. build_fast_context builds  !
! the run's fast_context_t: the column configuration (the soil column, the PFT table's    !
! hydraulics and film tables, the integrator options) and the PFT table's canopy optics.  !
! site_column_configs gives each site its copy, with its own soil and its own leaf table, !
! acclimated to its growth temperature. Nothing below this module reads meds_config_t:     !
! the kernels and the numerics take these records.                                         !
!==========================================================================================!
module meds_fast_config
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t
   use meds_soil_derive, only : build_soil_profile
   use meds_pft_params,  only : pft_table_t, rad_pft_optics_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_fast_types,  only : column_config_t, column_site_t, tol_set_t, error_control_t, integrator_opts_t, GRP_THETA
   implicit none
   private

   public :: fast_context_t, build_fast_context, site_column_configs
   public :: build_integrator_opts
   public :: apply_canopy_film_config
   public :: build_tol_set   ! used only by the tests

   !----- What the fast loop needs beyond the polygon and cfg (build_fast_context): the column      !
   !      configuration, the canopy optics and the ground surface. ------------------------------!
   type :: fast_context_t
      !----- The column configuration: the run part every site shares, and the site part [soil] and  !
      !      the PFT table set, from which each site's own is made (site_column_configs). ------------!
      type(column_config_t) :: col_config
      type(column_site_t)   :: col_site
      !----- [kg/m3] The air density the slow step values the canopy air at when its depth changes. --!
      real(wp) :: rho_air = 1.2_wp
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
      call build_soil_profile(cfg%soil, ctx%col_site%soil_profile)
      !----- The PFT table's hydraulics and films (derive_pft_tables); each site's leaf table is the  !
      !      site's (site_column_configs). ------------------------------------------------------------!
      ctx%col_config%hydraulics_table = cfg%pft%hydraulics
      call apply_canopy_film_config(cfg%pft, ctx%col_config%leaf_surf_water_max, ctx%col_config%wood_surf_water_max)
      !----- The canopy interception film (+ film evaporation, dew) on/off, [fast] block. --------!
      ctx%col_config%canopy_water_on    = cfg%canopy_water_on
      !----- The [soil], [snow], [aerodynamics] and [fast] options, copied as they are. -----------!
      ctx%col_site%soil = cfg%soil        ! [soil] -> soil water and heat (each site's own: site_column_configs)
      ctx%col_config%intercept_alpha = cfg%intercept_alpha   ! [fast] -> the canopy interception
      ctx%col_config%intercept_k     = cfg%intercept_k
      ctx%soil_albedo(1) = cfg%soil%ground_albedo_vis   ! [soil] ground optics -> the canopy radiation solver
      ctx%soil_albedo(2) = cfg%soil%ground_albedo_nir
      ctx%soil_emiss     = cfg%soil%ground_emissivity
      ctx%col_site%debug_error = cfg%debug_error         ! [fast] -> the budget stops
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
         !----- A budget-closure threshold, not a tolerance group: atol_scale scales it directly. -!
         ctx%col_site%soil_energy_budget_atol  = cfg%soil_energy_budget_atol * cfg%atol_scale
      end associate


      !----- The PFT table's canopy-RT optics, read-only downstream. -----------------------------!
      ctx%rad_opt = cfg%pft%optics
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
            col%soil_profile = surf%soil_profile
            col%soil         = surf%soil
            col%leaf_photo   = poly%site%leaf_photo(is)
         end associate
      end do
   end subroutine site_column_configs

   !----- The canopy films' per-PFT capacities, from the PFT table (both default to 0.1). -------!
   subroutine apply_canopy_film_config(pft, leaf_max, wood_max)
      type(pft_table_t),     intent(in)  :: pft
      real(wp), allocatable, intent(out) :: leaf_max(:), wood_max(:)
      leaf_max = pft%leaf_surf_water_max(1:pft%n)
      wood_max = pft%wood_surf_water_max(1:pft%n)
   end subroutine apply_canopy_film_config

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
      tols%atol(GRP_THETA)  = cfg%water_atol
      !----- Every atol scales with atol_scale: rtol alone stops tightening once atol dominates the   !
      !      WRMS denominator atol + rtol*|y| (docs/science/numerical_scheme.md §6). ----------------!
      tols%atol = tols%atol * cfg%atol_scale
   end function build_tol_set

   !----- The error-control bundle: the tolerance set and the step controller. ----------------------!
   pure function build_error_control(cfg) result(ec)
      type(meds_config_t), intent(in) :: cfg
      type(error_control_t)           :: ec
      ec%tols       = build_tol_set(cfg)
      ec%controller = cfg%step_controller
   end function build_error_control

   !----- The integrator's whole configuration, once per run (carried on column_config_t%integrator). -!
   pure function build_integrator_opts(cfg) result(opts)
      type(meds_config_t), intent(in) :: cfg
      type(integrator_opts_t)         :: opts
      opts%scheme           = cfg%time_integrator
      opts%adaptive         = cfg%ark_adaptive
      opts%fixed_substeps   = cfg%ark_fixed_substep
      opts%error_control    = build_error_control(cfg)
   end function build_integrator_opts

end module meds_fast_config
