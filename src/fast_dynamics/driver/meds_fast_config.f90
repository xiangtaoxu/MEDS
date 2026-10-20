! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_config -- the ONE seam from meds_config_t to the fast loop's options leaves.   !
!                                                                                          !
! Every procedure here reads the whole run configuration and produces a self-contained     !
! record the fast loop carries: the per-PFT leaf-photosynthesis table, the canopy optics   !
! table, the tolerance set, the error controller, the integrator knobs. The builders run   !
! ONCE per run, from build_fast_context; acclimate_leaf_photo_table refreshes the leaf     !
! table each slow step, and leaf_gas_exchange is the single-leaf wrapper the tests call.   !
!                                                                                          !
! They are driver code by definition (MEDS_CODE_STRUCTURE_DESIGN.md §6, rule 5): a kernel  !
! takes an OPTIONS LEAF, never a config handle, so a procedure that takes meds_config_t    !
! lives here, not in a kernel, a facade (rule 6) or a numerics module. NOTHING below this  !
! module reads meds_config_t: the kernels take leaves, and the leaves are built here.      !
!==========================================================================================!
module meds_fast_config
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t, pft_stomata_psi_onset, pft_wstress_psi_open,           &
                                pft_wstress_psi_close
   use meds_fast_opts,   only : CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_therm_lib,     only : kattge_knorr_entropy, kattge_knorr_jv_ratio
   use meds_fast_types,  only : tol_set_t, error_control_t, integrator_opts_t,                     &
                                GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W,        &
                                GRP_THETA, N_TOL_GROUP
   use meds_fast_control, only : default_tol_set, default_error_control
   use meds_canopy_types,     only : rad_pft_optics_t, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use meds_canopy_radiation, only : derive_rad_optics
   use meds_optics_lib,       only : moments_to_beta_params
   implicit none
   private

   public :: leaf_photo_params_for_pft, build_leaf_photo_table, leaf_gas_exchange
   public :: acclimate_leaf_photo_table
   public :: build_rad_optics
   public :: build_tol_set, build_error_control, build_integrator_opts

contains

   !---------------------------------------------------------------------------------------!
   !----- leaf_photo_params_for_pft -- PFT ipft's leaf-photosynthesis parameters, with the TABLE's  !
   !      Vcmax25/Jmax25/TPU25/Rd25. The ONE place the configuration is flattened into the kernel's !
   !      parameter record: leaf_gas_exchange applies a cohort's plastic overrides on top of it, and !
   !      build_leaf_photo_table calls it once per PFT. ---------------------------------------------!
   pure subroutine leaf_photo_params_for_pft(cfg, ipft, p)
      type(meds_config_t),       intent(in)  :: cfg
      integer(ik),               intent(in)  :: ipft
      type(leaf_photo_params_t), intent(out) :: p
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
      p%kc25 = cfg%kc25 ; p%ko25 = cfg%ko25 ; p%gstar25 = cfg%gstar25
      p%ea_kc = cfg%ea_kc ; p%ea_ko = cfg%ea_ko ; p%ea_gstar = cfg%ea_gstar
      p%ea_vcmax = cfg%ea_vcmax ; p%ea_jmax = cfg%ea_jmax ; p%ea_rd = cfg%ea_rd
      p%hd_vcmax = cfg%hd_vcmax ; p%hd_jmax = cfg%hd_jmax ; p%hd_rd = cfg%hd_rd
      p%ds_vcmax = cfg%ds_vcmax ; p%ds_jmax = cfg%ds_jmax ; p%ds_rd = cfg%ds_rd
      p%o2_mol_frac = cfg%o2_mol_frac ; p%absorptance = cfg%leaf_absorptance ; p%phi_psii = cfg%phi_psii
      p%medlyn_vpd_min = cfg%medlyn_vpd_min
   end subroutine leaf_photo_params_for_pft

   !----- build_leaf_photo_table -- every PFT's parameters plus the run-level solver selectors,     !
   !      once per run. The fast loop's column configuration carries the result. -----------------!
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
   !  acclimate_leaf_photo_table -- refresh the peaked-form entropy terms and the          !
   !  Jmax:Vcmax ratio from a running-mean GROWTH temperature (Kattge & Knorr 2007;        !
   !  docs/science/leaf_gas_exchange.md eq. (4)).                                          !
   !                                                                                       !
   !  WHY THE TABLE AND NOT THE KERNEL. dS and the capacity ratio are already table        !
   !  fields, so acclimation is a parameter refresh, not a new quantity the solver has to  !
   !  be told about -- no kernel signature changes, no per-cohort argument, and the law    !
   !  stays in one place. The REFERENCE values live in cfg and are never overwritten, so   !
   !  this is idempotent: calling it twice with the same temperature gives the same table, !
   !  and it can be called every slow step without drifting.                               !
   !                                                                                       !
   !  SITE-LEVEL, not per cohort. Kattge & Knorr fit dS against the mean AIR temperature   !
   !  of the preceding weeks, so a per-cohort leaf temperature would apply the relation    !
   !  well outside its calibration. Per-cohort acclimation is a separate question (see     !
   !  docs/ROADMAP.md).                                                                    !
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

   !----- build_rad_optics -- the canopy two-stream's per-PFT optics table                          !
   !      (docs/science/canopy_radiation_transfer.md §1-§2), from the [pft] trait table, once per   !
   !      run. Per-PFT, so two PFTs can differ in how they intercept light. The fast loop builds it !
   !      in build_fast_context; the canopy C API (meds.canopy) builds the same. -------------------!
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
      !----- PER-PFT, from the [pft] table. Shortwave arrives as reflectance and transmittance;  !
      !      LONGWAVE arrives as emissivity, and the band's reflectance is 1 - emissivity with   !
      !      zero transmittance, because a leaf is opaque at thermal wavelengths. That is        !
      !      physics, so it is derived here rather than offered as two more knobs a user could   !
      !      set inconsistently. ----------------------------------------------------------------!
      associate (t => cfg%pft)
         rl(RAD_VIS,1:np) = t%leaf_reflect_vis(1:np) ; tl(RAD_VIS,1:np) = t%leaf_transmit_vis(1:np)
         rl(RAD_NIR,1:np) = t%leaf_reflect_nir(1:np) ; tl(RAD_NIR,1:np) = t%leaf_transmit_nir(1:np)
         rl(RAD_LW ,1:np) = 1.0_wp - t%leaf_emissivity(1:np) ; tl(RAD_LW,1:np) = 0.0_wp
         rw(RAD_VIS,1:np) = t%wood_reflect_vis(1:np) ; tw(RAD_VIS,1:np) = t%wood_transmit_vis(1:np)
         rw(RAD_NIR,1:np) = t%wood_reflect_nir(1:np) ; tw(RAD_NIR,1:np) = t%wood_transmit_nir(1:np)
         rw(RAD_LW ,1:np) = 1.0_wp - t%wood_emissivity(1:np) ; tw(RAD_LW,1:np) = 0.0_wp
         cl(1:np) = t%leaf_clumping(1:np) ; cw(1:np) = t%wood_clumping(1:np)
         !----- The Beta leaf-angle shape is a TRANSFORM of (mean, std), so it is derived per  !
         !      PFT rather than configured: the two shape parameters are not quantities anyone  !
         !      measures, and an inconsistent pair has no leaf-angle distribution behind it. ---!
         do ipf = 1_ik, np
            call moments_to_beta_params(t%leaf_angle_mean(ipf), t%leaf_angle_std(ipf), bpp, bqq)
            bp(ipf) = bpp ; bq(ipf) = bqq
         end do
      end associate
      hb = [.true.,  .true.,  .false.]                         ! VIS/NIR have a beam; LW does not
      he = [.false., .false., .true. ]                         ! only LW emits
      call derive_rad_optics(NB, np, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, optics)
   end subroutine build_rad_optics

   subroutine leaf_gas_exchange(env, cfg, ipft, flux, vcmax25, rd25)
      type(leaf_env_t),    intent(in)  :: env
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ipft
      type(leaf_flux_t),   intent(out) :: flux
      real(wp), optional,  intent(in)  :: vcmax25   !< per-cohort (plastic) Vcmax25 override; jmax/tpu scale with it
      real(wp), optional,  intent(in)  :: rd25      !< per-cohort (plastic) Rd25 override
      type(leaf_photo_params_t)        :: p
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
   ! build_tol_set -- THE single tolerance source for the whole fast loop: the per-group   !
   ! tolerances of the grouped error norm (docs/science/numerical_scheme.md eq. (3)).      !
   !   * every group's rtol <- [fast].rtol, one relative accuracy for the whole loop;      !
   !   * ARK/RK45-integrated groups (enthalpy, shv, CO2, soil energy, leaf_w, wood_w) keep !
   !     each group's default atol;                                                        !
   !   * GRP_THETA's atol <- [fast].water_atol, the soil-plant water network's control;    !
   !   * GRP_LEAF_W/GRP_WOOD_W: only RK45 integrates plant water mass in its stages; ARK   !
   !     advances it outside the ESDIRK tableau, so these terms are structurally zero      !
   !     there. Seeded regardless, so every group exists uniformly.                        !
   !                                                                                       !
   !---------------------------------------------------------------------------------------!
   pure function build_tol_set(cfg) result(tols)
      type(meds_config_t), intent(in) :: cfg
      type(tol_set_t)                 :: tols
      !----- One relative tolerance for every group; each group keeps its default atol, and the   !
      !      water network's is [fast].water_atol. ----------------------------------------------------!
      tols%rtol             = cfg%rtol
      tols%atol(GRP_THETA)  = cfg%soil%water_atol
      !----- ...and the ABSOLUTE companion. The WRMS denominator is atol + rtol*|y|, so rtol         !
      !      alone SATURATES once atol dominates (docs/science/numerical_scheme.md §6): tightening    !
      !      rtol from 1e-3 to 1e-6 raises the soil-water error estimate only ~4x (4e-4 -> 1e-4       !
      !      denominator) -- never enough to force a substep -- while scaling BOTH does. atol is      !
      !      dimensional and differs per group, so it scales rather than broadcasts. The default 1.0  !
      !      is an exact IEEE identity, so it changes nothing. ---------------------------------------!
      tols%atol = tols%atol * cfg%atol_scale
   end function build_tol_set


   !----- The full error-control bundle from config: unified tolerances + controller + strictness. ----!
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
