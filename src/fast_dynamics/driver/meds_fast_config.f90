!==========================================================================================!
! meds_fast_config -- the ONE seam from meds_config_t to the fast loop's options leaves.     !
!                                                                                          !
! Every procedure here reads the whole run configuration and produces a self-contained record !
! the fast loop carries: the per-PFT leaf-photosynthesis table, the tolerance set, the error   !
! controller, the integrator knobs. They run ONCE per run, from build_fast_context.             !
!                                                                                          !
! They are here rather than where they used to live for one reason (structure-plan rule 5 and   !
! rule 6): a kernel takes an OPTIONS LEAF, never a config handle, so a procedure that takes      !
! meds_config_t is by definition driver code. Two of these sat in meds_fast_config, which     !
! made a facade that was supposed to be pure re-export also carry logic; three sat in              !
! meds_fast_control, a numerics module, which now keeps only the norm and the step controller.     !
!                                                                                          !
! Consequence worth stating: NOTHING below this module reads meds_config_t. The kernels take       !
! leaves, and the leaves are built here.                                                            !
!==========================================================================================!
module meds_fast_config
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t, CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT,    &
                                CTRL_I, CTRL_PI
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_hydr_lib,    only : pv_psi_tlp
   use meds_fast_types,  only : tol_set_t, error_control_t, integrator_opts_t,                     &
                                GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W,        &
                                GRP_THETA, GRP_SOIL_T, N_TOL_GROUP
   use meds_fast_control, only : default_tol_set, default_error_control
   implicit none
   private

   public :: leaf_photo_params_for_pft, build_leaf_photo_table, leaf_gas_exchange
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
         p%theta_cj       = t%theta_cj_c4(ipft)
         p%theta_ic       = t%theta_ic_c4(ipft)
         p%lambda25       = t%katul_lambda25(ipft)
         p%psi_open       = t%wstress_psi_open(ipft)
         p%psi_close      = t%wstress_psi_close(ipft)
         p%lambda_psi_exp = t%wstress_lambda_exp(ipft)
         p%sref_stomata   = t%wstress_sref_stomata(ipft)
         p%psi_tlp        = pv_psi_tlp(cfg%hydraulics%leaf_pi0, cfg%hydraulics%leaf_elastic_mod)
      end associate
      p%wstress_nonstomatal = cfg%leaf_wstress_nonstomatal
      p%stress_arrestor     = cfg%leaf_stress_arrestor
      !----- Copy the shared biochemistry constants. --------------------------------------!
      p%kc25 = cfg%kc25 ; p%ko25 = cfg%ko25 ; p%gstar25 = cfg%gstar25
      p%ea_kc = cfg%ea_kc ; p%ea_ko = cfg%ea_ko ; p%ea_gstar = cfg%ea_gstar
      p%ea_vcmax = cfg%ea_vcmax ; p%ea_jmax = cfg%ea_jmax ; p%ea_rd = cfg%ea_rd
      p%hd_vcmax = cfg%hd_vcmax ; p%hd_jmax = cfg%hd_jmax ; p%hd_rd = cfg%hd_rd
      p%ds_vcmax = cfg%ds_vcmax ; p%ds_jmax = cfg%ds_jmax ; p%ds_rd = cfg%ds_rd
      p%o2_mol_frac = cfg%o2_mol_frac ; p%absorptance = cfg%leaf_absorptance ; p%phi_psii = cfg%phi_psii
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
      table%use_boundary_layer = cfg%leaf_use_boundary_layer
   end subroutine build_leaf_photo_table

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
                                   cfg%colimitation, cfg%leaf_use_boundary_layer, flux)
   end subroutine leaf_gas_exchange


   !---------------------------------------------------------------------------------------!
   ! build_tol_set -- THE single tolerance source for the whole fast loop (§8c Layer 1). Each group   !
   ! is SEEDED from the setting that governs it today, so the result is byte-identical to the         !
   ! pre-unification behaviour:                                                                        !
   !   * ARK/RK45-integrated groups (enthalpy/shv/CO2/soil-energy/leaf_w/wood_w) <- [fast].ark_rtol +   !
   !     historical atols;                                                                              !
   !   * GRP_THETA   <- the [soil]   sub-solver's own (rtol, atol)  -- soil-water Richards step-doubling; !
   !   * GRP_SOIL_T  <- the [energy] sub-solver's own (rtol, atol)  -- soil-energy substepping;           !
   !   * GRP_LEAF_W/GRP_WOOD_W (MEDS_ED2_RK45_DESIGN.md sec 6, P2, replaces the retired GRP_PSI): only     !
   !     RK45 actually folds these into its embedded-error WRMS (mass is operator-split out of ARK's      !
   !     ESDIRK tableau, like psi was, via with_mass=.false.) -- seeded here regardless so the group        !
   !     exists uniformly. --------------------------------------------------------------------------------!
   !                                                                                          !
   ! ONE MASTER DIAL: when cfg%rtol_all > 0 it OVERRIDES every group's rtol, so a single number sets the  !
   ! relative accuracy of the entire hierarchy (the "target accuracy" axis goals (b)/(c) need). Left at   !
   ! its 0 default, each group keeps its own per-sub-solver value => byte-identical, even for a config    !
   ! that already customised [soil]/[energy] tolerances.                                                  !
   !---------------------------------------------------------------------------------------!
   pure function build_tol_set(cfg) result(tols)
      type(meds_config_t), intent(in) :: cfg
      type(tol_set_t)                 :: tols
      !----- ARK/RK45-integrated groups: the single ark_rtol, historical atols (atol defaults kept). --!
      tols%rtol(GRP_ENTH)   = cfg%ark_rtol
      tols%rtol(GRP_SHV)    = cfg%ark_rtol
      tols%rtol(GRP_CO2)    = cfg%ark_rtol
      tols%rtol(GRP_SE)     = cfg%ark_rtol
      tols%rtol(GRP_LEAF_W) = cfg%ark_rtol
      tols%rtol(GRP_WOOD_W) = cfg%ark_rtol
      !----- Sub-solver groups: seed from the opts that drive them today. ----------------------------!
      tols%rtol(GRP_THETA)  = cfg%soil%rtol   ; tols%atol(GRP_THETA)  = cfg%soil%atol
      tols%rtol(GRP_SOIL_T) = cfg%energy%rtol ; tols%atol(GRP_SOIL_T) = cfg%energy%atol
      !----- The one master accuracy dial (0 => unset => keep the per-group values above). ------------!
      if (cfg%rtol_all > 0.0_wp) tols%rtol = cfg%rtol_all
      !----- ...and its ABSOLUTE companion. The WRMS denominator is atol + rtol*|y|, so rtol_all alone   !
      !      SATURATES once atol dominates: measured on the split path, rtol 1e-3 -> 1e-6 raises the      !
      !      soil-water error estimate only ~4x (4e-4 -> 1e-4 denominator) -- never enough to force a     !
      !      substep, while scaling BOTH does. atol is dimensional and differs per group, so it scales    !
      !      rather than broadcasts. The default 1.0 is an exact IEEE identity => byte-identical. --------!
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
