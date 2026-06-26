!==========================================================================================!
! meds_recruitment -- the recruitment-RATE provider (physiology layer).                     !
!                                                                                          !
! Determines the potential per-(PFT, patch) recruitment density [plant m-2 month-1]. This is  !
! a RATE: like the other physiology evaluators it READS the site and returns plain data, which !
! the master stepper then hands to the demography engine through `update_demography`; the      !
! engine accumulates it and spawns new cohorts (meds_demography_dynamics::apply_recruitment).  !
!                                                                                          !
! The current rule is a constant per-PFT density, gated only by include_pft (no environmental  !
! control). This module is the seam where future, physiology-driven recruitment (seed          !
! production, carbon allocation, light/seedling dynamics, ...) will live.                      !
!==========================================================================================!
module meds_recruitment
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_demography_types, only : site_t
   implicit none
   private

   public :: recruitment_rate

contains

   subroutine recruitment_rate(site, cfg, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: ip, pf, np

      np = site%patch%n
      allocate(recruitment(cfg%pft%n, max(np, 1_ik)))
      recruitment = 0.0_wp
      associate (pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, pft%n
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%recruit_dens(pf)
            end do
         end do
      end associate
   end subroutine recruitment_rate

end module meds_recruitment
