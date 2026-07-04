!==========================================================================================!
! meds_plant_phenology -- THE sealed public interface of the leaf-phenology library (the    !
! phenological analogue of meds_leaf_physiology / meds_plant_hydraulics).                      !
!                                                                                          !
! One public routine, update_phenology, advances a cohort's cue accumulators over one daily     !
! step and returns the phenological STATUS -- the direction of leaf-display change             !
! {ON,OFF,DORMANT} and the governing cue. It is the stateless-in-leaf-mass compute half of the   !
! *Mem/compute split: the prognostic cue memory lives in the cohort SoA (pheno_state_t) and is    !
! passed by argument; the displayed leaf fraction is NOT here -- it emerges downstream (a future   !
! leaf-dynamics module) from this status plus leaf growth/shed rates. Everything else in            !
! src/plant/phenology/ is internal; production callers use ONLY this module. The public types are   !
! re-exported for convenience.                                                                       !
!                                                                                          !
! (The config-driven seam update_phenology(env, cfg, ipft, ...) that flattens cfg%pft into           !
! pheno_params_t is added once meds_config gains a [phenology] block; for now the flat params are     !
! passed explicitly, exactly as the leaf/hydraulics kernels take flat params.)                         !
!==========================================================================================!
module meds_plant_phenology
   use meds_kinds,        only : wp
   use meds_pheno_types,  only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,      &
                                 CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO,          &
                                 PHEN_ON, PHEN_DORMANT, PHEN_OFF
   use meds_pheno_engine, only : phenology_kernel, daylength
   implicit none
   private

   !----- Re-exported so callers need only this module. ------------------------------------!
   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO
   public :: PHEN_ON, PHEN_DORMANT, PHEN_OFF
   public :: update_phenology, daylength

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one cohort's phenology over one daily step dt: update the cue accumulators in    !
   ! `state` and return the directional status in `out`.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine update_phenology(env, params, dt, state, out)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      type(pheno_out_t),    intent(out)   :: out
      call phenology_kernel(env, params, dt, state, out)
   end subroutine update_phenology

end module meds_plant_phenology
