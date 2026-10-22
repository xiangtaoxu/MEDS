! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_config -- builds the fast loop's option records from meds_config_t, once per   !
! run from build_fast_context: the per-PFT leaf-photosynthesis table                       !
! (build_leaf_photo_table; acclimate_leaf_photo_table refreshes it each slow step), the    !
! canopy optics table (build_rad_optics) and the integrator options with their tolerance   !
! set (build_integrator_opts). The kernels take these records; column_fast_step and the    !
! numerics below it also take cfg, for the PFT leaf and wood emissivities, the stomatal    !
! closing potential and the [soil_carbon] settings.                                        !
!==========================================================================================!
module meds_fast_config
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t, pft_stomata_psi_onset, pft_wstress_psi_open,           &
                                pft_wstress_psi_close
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_args_t, leaf_photo_table_t
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_therm_lib,     only : kattge_knorr_entropy, kattge_knorr_jv_ratio
   use meds_fast_types,  only : tol_set_t, error_control_t, integrator_opts_t, GRP_THETA
   use meds_canopy_types,     only : rad_pft_optics_t, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use meds_canopy_radiation, only : derive_rad_optics
   use meds_optics_lib,       only : moments_to_beta_params
   implicit none
   private

   public :: build_leaf_photo_table
   public :: acclimate_leaf_photo_table
   public :: build_rad_optics
   public :: build_integrator_opts
   public :: leaf_photo_params_for_pft, leaf_gas_exchange, build_tol_set   ! used only by the tests

contains

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
