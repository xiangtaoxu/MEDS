!==========================================================================================!
! meds_plant_interface -- the public façade of the plant-ecophysiology library (the plant-  !
! level analogue of meds_core_interface): it RE-EXPORTS the compute kernels + their public   !
! types so production callers can `use` this one module. The only genuine WRAPPER is           !
! leaf_gas_exchange, which flattens cfg%pft into a self-contained leaf_photo_params_t before    !
! driving the coupled Ci solver; every other symbol (solve_plant_water, phenology_kernel,        !
! pheno_drives_to_rates, plant_carbon_allocation, ...) is re-exported verbatim.                   !
!                                                                                          !
! ORCHESTRATION is NOT here -- it lives in the drivers, which sequence these kernels + the        !
! cohort/patch state: the SLOW loop in meds_vegetation_dynamics (carbon growth, phenology,         !
! turnover) and the FAST loop in meds_fast_split/meds_fast_ark (the coupled leaf<->CAS<->soil<->    !
! hydraulics fixed point). The former coarse "get_plant_flux_{fast,slow}" seams were removed: the   !
! coupling is a whole-column concern, not a per-plant call.                                         !
!                                                                                          !
! This is a CONVENIENCE façade, not a sealed wall (DAG hygiene is enforced by the library link      !
! graph, not here). Two legitimate access patterns coexist: BLACK-BOX callers `use` this module     !
! for the common types + solve-style kernels; WHITE-BOX callers -- the fast-loop numerical           !
! integrators (meds_fast_ark, meds_fast_time_derivs) -- `use` the kernel modules directly            !
! (meds_plant_hydraulics, meds_hydr_lib) because they need the RHS/tendency + constitutive        !
! curves at each stage, which a per-call solve seam cannot expose. meds_plant_vital_rates is         !
! likewise imported directly by the slow driver (its single consumer).                              !
!==========================================================================================!
module meds_plant_interface
   use meds_kinds,       only : wp, ik
   use meds_config,      only : meds_config_t
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_params_t,                   &
                                hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,        &
                                N_HYDRO, NODE_LEAF, NODE_STEM, NODE_WOOD, NODE_ROOT,            &
                                HYDRO_NODES_2, HYDRO_NODES_3,                                   &
                                HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE,                             &
                                HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT,                          &
                                HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED,                    &
                                pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,        &
                                CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO, CUE_LIGHT, &
                                wood_params_t, root_params_t
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_plant_hydraulics,  only : solve_plant_water, solve_plant_water_batch
   use meds_phenology,         only : phenology_kernel, pheno_drives_to_rates, turnover_shed_rates
   use meds_plant_respiration, only : stem_maintenance_respiration,                            &
                                      fine_root_maintenance_respiration
   use meds_plant_carbon_allocation, only : plant_carbon_allocation, growth_respiration
   use meds_hydr_lib, only : pv_psi_tlp
   implicit none
   private

   !----- Re-exported public types + constants, so callers need only this module. ----------!
   public :: leaf_env_t, leaf_flux_t
   public :: hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t
   public :: N_HYDRO, NODE_LEAF, NODE_STEM, NODE_WOOD, NODE_ROOT
   public :: HYDRO_NODES_2, HYDRO_NODES_3, HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE
   public :: HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT, HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED
   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO, CUE_LIGHT
   public :: wood_params_t, root_params_t
   !----- The seams. leaf_gas_exchange is a genuine wrapper (it flattens cfg%pft into a self-  !
   !      contained leaf_photo_params_t); everything else is a plain RE-EXPORT of a kernel --   !
   !      the orchestration lives in the drivers (slow: meds_vegetation_dynamics; fast:          !
   !      meds_fast_split/meds_fast_ark), so a per-plant "flux seam" wrapper would only add        !
   !      indirection.                                                                          !
   public :: leaf_gas_exchange, leaf_gas_exchange_batch
   public :: solve_plant_water, solve_plant_water_batch, phenology_kernel, pheno_drives_to_rates, turnover_shed_rates
   public :: stem_maintenance_respiration, fine_root_maintenance_respiration
   public :: plant_carbon_allocation, growth_respiration

contains

   !---------------------------------------------------------------------------------------!
   ! Leaf gas exchange: flatten the per-PFT traits + shared biochemistry constants into a    !
   ! self-contained leaf_photo_params_t, read the model selectors from cfg, and drive the     !
   ! coupled A-gs-Ci solver for PFT `ipft` under the environment `env`.                       !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_gas_exchange(env, cfg, ipft, flux, vcmax25, rd25)
      type(leaf_env_t),    intent(in)  :: env
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ipft
      type(leaf_flux_t),   intent(out) :: flux
      real(wp), optional,  intent(in)  :: vcmax25   !< per-cohort (plastic) Vcmax25 override; jmax/tpu scale with it
      real(wp), optional,  intent(in)  :: rd25      !< per-cohort (plastic) Rd25 override
      type(leaf_photo_params_t)        :: p

      !----- Flatten per-PFT traits. Plastic Vcmax (if supplied) carries Jmax/TPU with it via the  !
      !      fixed ratios; plastic Rd is its own value. Absent => the static PFT table (identical). !
      associate (t => cfg%pft)
         p%pathway        = t%photosynthetic_pathway(ipft)
         if (present(vcmax25)) then
            p%vcmax25 = vcmax25
            p%jmax25  = t%jmax_vcmax_ratio(ipft) * vcmax25
            p%tpu25   = t%tpu_vcmax_ratio(ipft)  * vcmax25
         else
            p%vcmax25 = t%vcmax25(ipft)
            p%jmax25  = t%jmax25(ipft)
            p%tpu25   = t%tpu25(ipft)
         end if
         if (present(rd25)) then
            p%rd25 = rd25
         else
            p%rd25 = t%rd25(ipft)
         end if
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
         !----- Turgor-loss point from the SAME PV curve the hydraulics solver uses, so the two      !
         !      cannot describe different leaves. ---------------------------------------------!
         p%psi_tlp        = pv_psi_tlp(cfg%hydraulics%leaf_pi0, cfg%hydraulics%leaf_elastic_mod)
      end associate

      p%wstress_nonstomatal = cfg%leaf_wstress_nonstomatal

      !----- Copy the shared biochemistry constants. --------------------------------------!
      p%kc25 = cfg%kc25 ; p%ko25 = cfg%ko25 ; p%gstar25 = cfg%gstar25
      p%ea_kc = cfg%ea_kc ; p%ea_ko = cfg%ea_ko ; p%ea_gstar = cfg%ea_gstar
      p%ea_vcmax = cfg%ea_vcmax ; p%ea_jmax = cfg%ea_jmax ; p%ea_rd = cfg%ea_rd
      p%hd_vcmax = cfg%hd_vcmax ; p%hd_jmax = cfg%hd_jmax ; p%hd_rd = cfg%hd_rd
      p%ds_vcmax = cfg%ds_vcmax ; p%ds_jmax = cfg%ds_jmax ; p%ds_rd = cfg%ds_rd
      p%o2_mol_frac = cfg%o2_mol_frac ; p%absorptance = cfg%leaf_absorptance ; p%phi_psii = cfg%phi_psii

      call solve_leaf_gas_exchange(env, p, cfg%stomatal_model, cfg%temp_response_form,         &
                                   cfg%colimitation, cfg%leaf_use_boundary_layer, flux)
   end subroutine leaf_gas_exchange

   !---------------------------------------------------------------------------------------!
   ! leaf_gas_exchange_batch -- BARE-ARRAY entry point over n leaves (MEDS_NUMERICS_SCOPING.md    !
   ! "bare-array process kernels": one call solves a whole array of leaf environments, so the      !
   ! per-cohort fast-loop driver loops need no longer thread the leaf_env_t/leaf_flux_t derived     !
   ! types, and a Python/ctypes wrapper (meds_plant_capi) can vectorise over numpy arrays instead   !
   ! of calling one leaf at a time). The per-leaf PHYSICS is UNCHANGED: this loops `do i=1,n`        !
   ! calling the SAME leaf_gas_exchange above, so it is bit-identical to an inline caller loop.      !
   !                                                                                          !
   ! CONVENTION (shared by every *_batch kernel): genuinely PER-ELEMENT quantities are bare arrays   !
   ! of length n; PATCH/RUN-UNIFORM quantities (here ca, pressure, and the whole cfg trait table)    !
   ! are passed as scalars/one config object (broadcast to every element); outputs are bare arrays.  !
   ! psi_soil is OPTIONAL (absent => 0, the leaf_env_t default = the drought-stomata limb inert,     !
   ! matching the current fast-loop wiring). Only the outputs the fast loop consumes (A_gross, gs,   !
   ! rd) are surfaced; add more leaf_flux_t fields here when a caller needs them.                    !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_gas_exchange_batch(n, par, leaf_temp, vpd, ca, pressure, psi_leaf, gb,      &
                                      cfg, pft, vcmax25, rd25, a_gross, gs, rd, psi_soil)
      integer(ik),         intent(in)  :: n
      real(wp),            intent(in)  :: par(n), leaf_temp(n), vpd(n), psi_leaf(n), gb(n)  !< per-leaf env
      real(wp),            intent(in)  :: ca, pressure                                      !< patch-uniform (broadcast)
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: pft(n)
      real(wp),            intent(in)  :: vcmax25(n), rd25(n)                               !< per-leaf plastic capacities
      real(wp),            intent(out) :: a_gross(n), gs(n), rd(n)
      real(wp), optional,  intent(in)  :: psi_soil(n)                                       !< absent => 0 (well-watered)
      type(leaf_env_t)  :: env
      type(leaf_flux_t) :: flux
      integer(ik) :: i
      do i = 1_ik, n
         env%par = par(i) ; env%leaf_temp = leaf_temp(i) ; env%vpd = vpd(i)
         env%ca = ca ; env%pressure = pressure ; env%psi_leaf = psi_leaf(i) ; env%gb = gb(i)
         if (present(psi_soil)) then ; env%psi_soil = psi_soil(i) ; else ; env%psi_soil = 0.0_wp ; end if
         call leaf_gas_exchange(env, cfg, pft(i), flux, vcmax25=vcmax25(i), rd25=rd25(i))
         a_gross(i) = flux%A_gross ; gs(i) = flux%gs ; rd(i) = flux%rd
      end do
   end subroutine leaf_gas_exchange_batch

end module meds_plant_interface
