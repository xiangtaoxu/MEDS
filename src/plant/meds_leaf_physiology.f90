!==========================================================================================!
! meds_leaf_physiology -- THE sealed public interface of the leaf-physiology library (the  !
! leaf-level analogue of meds_demography_interface).                                        !
!                                                                                          !
! One public routine, leaf_gas_exchange, takes a leaf environment, the run configuration and !
! a PFT index, and returns the leaf-level fluxes. It flattens the per-PFT traits (cfg%pft) and !
! the shared biochemistry constants (cfg) into a self-contained leaf_photo_params_t, reads the !
! model selectors (stomatal model, temperature-response form, co-limitation, boundary layer)   !
! from cfg, and drives the coupled solver. Everything else in src/plant/leaf/ is internal:      !
! production callers use ONLY this module. The leaf types are re-exported for convenience.      !
!==========================================================================================!
module meds_leaf_physiology
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_params_t
   use meds_leaf_solver,only : solve_leaf_gas_exchange
   implicit none
   private

   public :: leaf_env_t, leaf_flux_t                 ! re-exported so callers need only this module
   public :: leaf_gas_exchange

contains

   !---------------------------------------------------------------------------------------!
   ! Compute leaf-level net/gross assimilation, stomatal conductance, intercellular CO2 and !
   ! transpiration for PFT `ipft` under the environment `env`, using the models selected in   !
   ! the run configuration.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_gas_exchange(env, cfg, ipft, flux)
      type(leaf_env_t),    intent(in)  :: env
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ipft
      type(leaf_flux_t),   intent(out) :: flux
      type(leaf_photo_params_t)        :: p

      !----- Flatten per-PFT traits. ------------------------------------------------------!
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
      end associate

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

end module meds_leaf_physiology
