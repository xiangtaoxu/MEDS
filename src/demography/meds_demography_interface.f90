!==========================================================================================!
! meds_demography_interface -- the public API of the demography ENGINE: it re-exports the state !
! type `site_t` and the apply-PRIMITIVES (grow / die / recruit / disturb + fuse/fission + the    !
! competition sweep) so a single `use meds_demography_interface` gives a caller the whole engine  !
! surface.                                                                                        !
!                                                                                          !
! The engine is now PURE MECHANISM: every routine here APPLIES a supplied change to `site_t` (or  !
! restructures it) -- NONE computes a rate and NONE calls a plant kernel, so `demography _|_      !
! plant` holds. The ORCHESTRATION (compute the carbon rates via the plant vital-rate kernels,     !
! then sequence these primitives + the cadence) lives in the driver, meds_vegetation_dynamics.    !
! The former `update_demography` seam was dissolved into that driver.                             !
!==========================================================================================!
module meds_demography_interface
   use meds_ecosystem_state,      only : site_t, carbon_flux_block
   use meds_demography_dynamics,  only : growth_step, mortality_step, apply_growth,             &
                                         apply_patch_disturbance, apply_recruitment
   use meds_demography_fusefiss,  only : new_fuse_cohorts, terminate_cohorts, split_cohorts,    &
                                         new_fuse_patches, terminate_patches,                    &
                                         sort_cohorts, sort_patches
   use meds_competition,          only : compute_overtopping_lai
   implicit none
   private

   !----- The demography engine's public surface (state type + apply-primitives). ----------!
   public :: site_t, carbon_flux_block
   public :: growth_step, mortality_step, apply_growth, apply_recruitment, apply_patch_disturbance
   public :: new_fuse_cohorts, terminate_cohorts, split_cohorts, new_fuse_patches,             &
             terminate_patches, sort_cohorts, sort_patches
   public :: compute_overtopping_lai

end module meds_demography_interface
