!==========================================================================================!
! meds_stepper -- the master stepper (the ED2 ed_model analogue; seed of an all-process        !
! top-level loop). Lives in src/driver/, the home of the top-level utilities that wire the      !
! process modules together.                                                                    !
!                                                                                          !
! It owns only the CADENCE: each step it is told whether a month/year rolled over, and it       !
! drives the process modules on the appropriate timescale. Today the only process is the        !
! slow-loop VEGETATION DYNAMICS (meds_vegetation_dynamics), which assembles the demographic     !
! rates and applies them through the demography seam. A future master step would, in addition,  !
! call the fast-loop biophysics processes here.                                                 !
!==========================================================================================!
module meds_stepper
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t
   use meds_vegetation_dynamics,  only : vegetation_dynamics
   implicit none
   private

   public :: advance_one_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one step. The caller's calendar supplies the cadence flags; this routine passes  !
   ! them to the process drivers. (Future: fast-loop biophysics before/after the slow loop.)  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_one_step(site, cfg, is_new_month, is_new_year)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year

      !----- Slow loop: vegetation dynamics (rate assembly + demographic application). -------!
      call vegetation_dynamics(site, cfg, is_new_month, is_new_year)
   end subroutine advance_one_step

end module meds_stepper
