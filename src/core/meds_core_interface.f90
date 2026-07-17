!==========================================================================================!
! meds_core_interface -- the public API of the CORE ecosystem-structure ENGINE: it re-exports    !
! the state type `site_t` and the apply-PRIMITIVES (grow / die / recruit / disturb + fuse/fission !
! + the competition sweep) so a single `use meds_core_interface` gives a caller the whole engine  !
! surface.                                                                                        !
!                                                                                          !
! The engine is PURE MECHANISM: every routine here APPLIES a supplied change to `site_t` (or  !
! restructures it) -- NONE computes a rate and NONE calls a plant kernel, so `core _|_ plant`     !
! holds. The ORCHESTRATION (compute the carbon rates via the plant vital-rate kernels, then       !
! sequence these primitives + the cadence) lives in the driver, meds_vegetation_dynamics.         !
! The former `update_demography` seam was dissolved into that driver.                             !
!==========================================================================================!
module meds_core_interface
   use meds_core_state_types,      only : site_t, carbon_flux_block, cohort_deriv_block,         &
                                          cohort_deriv_alloc
   use meds_core_state_update,     only : update_cohort_states, update_overtopping_lai
   use meds_core_cohort_fusefiss,  only : new_fuse_cohorts, terminate_cohorts, split_cohorts,   &
                                          sort_cohorts, apply_recruitment
   use meds_core_patch_fusefiss,   only : new_fuse_patches, terminate_patches, sort_patches,    &
                                          apply_patch_disturbance
   implicit none
   private

   !----- The core engine's public surface (state type + apply-primitives). ----------------!
   public :: site_t, carbon_flux_block, cohort_deriv_block, cohort_deriv_alloc
   public :: update_cohort_states, apply_recruitment, apply_patch_disturbance
   public :: new_fuse_cohorts, terminate_cohorts, split_cohorts, new_fuse_patches,             &
             terminate_patches, sort_cohorts, sort_patches
   public :: update_overtopping_lai

end module meds_core_interface
