! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_opts -- the selector codes of the [fast] block, as a low-level config leaf:    !
! the time integrator and its error control (docs/science/numerical_scheme.md §3, eqs.    !
! (3)-(4)). The fast-loop numerics switch on them and import them from here.               !
!==========================================================================================!
module meds_fast_opts
   use meds_kinds, only : ik
   implicit none
   private

   public :: INTEG_ARK, INTEG_RK45
   public :: CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI

   !----- [fast].time_integrator: two schemes. "ark" is a 2-solve ESDIRK2, not an IMEX method   !
   !      (docs/science/numerical_scheme.md §3). ----------------------------------------------!
   integer(ik), parameter :: INTEG_ARK  = 2_ik   !< ESDIRK2 coupled implicit column (default)
   integer(ik), parameter :: INTEG_RK45 = 3_ik   !< adaptive Cash-Karp RK45, the accuracy baseline

   !----- [fast].error_level: fixed sub-steps, adaptive (default), or adaptive with a hard stop   !
   !      when the floor step misses the tolerance. [fast].step_controller: integral (default) or !
   !      Gustafsson's proportional-integral, which damps step-size hunting. ----------------------!
   integer(ik), parameter :: CTRL_L0_FIXED    = 0_ik
   integer(ik), parameter :: CTRL_L1_ADAPTIVE = 1_ik
   integer(ik), parameter :: CTRL_L2_STRICT   = 2_ik
   integer(ik), parameter :: CTRL_I  = 1_ik
   integer(ik), parameter :: CTRL_PI = 2_ik

end module meds_fast_opts
