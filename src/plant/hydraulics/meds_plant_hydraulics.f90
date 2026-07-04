!==========================================================================================!
! meds_plant_hydraulics -- THE sealed public interface of the plant-hydraulics library (the  !
! hydraulic analogue of meds_leaf_physiology / meds_canopy_radiation).                        !
!                                                                                          !
! One public routine, plant_water_flux, advances a cohort's node water potentials over one     !
! fast step given the boundary conditions (transpiration, soil potential, rhizosphere            !
! conductance, per-plant geometry/biomass), the per-PFT hydraulic traits, and the run options.    !
! It is the stateless compute half of the *Mem/compute split: the prognostic psi lives in the       !
! cohort SoA and is passed by argument. Everything else in src/plant/hydraulics/ is internal;         !
! production callers use ONLY this module. The public types are re-exported for convenience.            !
!                                                                                          !
! (The config-driven seam plant_water_flux(env, cfg, ipft, ...) that flattens cfg%pft into            !
! hydro_params_t is added once meds_config gains a [hydraulics] block; for now the flat params/opts    !
! are passed explicitly, exactly as the leaf SOLVER takes flat params + selectors.)                     !
!==========================================================================================!
module meds_plant_hydraulics
   use meds_kinds,        only : wp
   use meds_hydro_types,  only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,      &
                                 N_HYDRO, NODE_LEAF, NODE_STEM, NODE_WOOD, NODE_ROOT,          &
                                 HYDRO_NODES_2, HYDRO_NODES_3,                                 &
                                 HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE,                           &
                                 HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT,                        &
                                 HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED
   use meds_hydro_solver, only : solve_plant_water
   implicit none
   private

   !----- Re-exported so callers need only this module. ------------------------------------!
   public :: hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t
   public :: N_HYDRO, NODE_LEAF, NODE_STEM, NODE_WOOD, NODE_ROOT
   public :: HYDRO_NODES_2, HYDRO_NODES_3, HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE
   public :: HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT, HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED
   public :: plant_water_flux

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the node water potentials psi(1:N_HYDRO) over dt for one cohort/individual.      !
   !---------------------------------------------------------------------------------------!
   subroutine plant_water_flux(env, params, opts, dt, psi, flux)
      type(hydro_env_t),    intent(in)    :: env
      type(hydro_params_t), intent(in)    :: params
      type(hydro_opts_t),   intent(in)    :: opts
      real(wp),             intent(in)    :: dt
      real(wp),             intent(inout) :: psi(N_HYDRO)
      type(hydro_flux_t),   intent(out)   :: flux
      call solve_plant_water(env, params, opts, dt, psi, flux)
   end subroutine plant_water_flux

end module meds_plant_hydraulics
