! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_params -- the selector codes of the [fast] block, as a config leaf: the time   !
! integrator and its step controller (docs/science/numerical_scheme.md §3, eq. (4)), and   !
! RK45's tolerance. The fast-loop numerics switch on them and import them from here.       !
!==========================================================================================!
module meds_fast_params
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: INTEG_ESDIRK2, INTEG_RK45
   public :: CTRL_I, CTRL_PI
   public :: RK45_RTOL

   !----- [fast].time_integrator (docs/science/numerical_scheme.md §3). -------------------!
   integer(ik), parameter :: INTEG_ESDIRK2 = 2_ik   !< "esdirk2": the coupled implicit column (default)
   integer(ik), parameter :: INTEG_RK45    = 3_ik   !< "rk45": adaptive Cash-Karp, the explicit alternative

   !----- RK45's [fast].rtol when the key is absent, and the loosest it accepts. Its step is set   !
   !      by stability, not by rtol, so a looser rtol saves almost no work while its fluxes drift  !
   !      far from the converged answer (docs/science/numerical_scheme.md §6). -------------------!
   real(wp), parameter :: RK45_RTOL = 1.0e-4_wp

   !----- [fast].step_controller. ---------------------------------------------------------!
   integer(ik), parameter :: CTRL_I  = 1_ik   !< integral (default)
   integer(ik), parameter :: CTRL_PI = 2_ik   !< Gustafsson's proportional-integral: damps step-size hunting

end module meds_fast_params
