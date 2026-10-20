! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_opts -- the selector codes and the process mask of the [fast] block, as a      !
! config leaf: the time integrator and its error control (docs/science/numerical_scheme.md !
! §3, eqs. (3)-(4)). The fast-loop numerics switch on them and import them from here.      !
!==========================================================================================!
module meds_fast_opts
   use meds_kinds, only : ik
   implicit none
   private

   public :: INTEG_ARK, INTEG_RK45
   public :: CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI
   public :: process_mask_t

   !----- [fast].time_integrator (docs/science/numerical_scheme.md §3). -------------------!
   integer(ik), parameter :: INTEG_ARK  = 2_ik   !< ESDIRK2 coupled implicit column (default)
   integer(ik), parameter :: INTEG_RK45 = 3_ik   !< adaptive Cash-Karp RK45, the accuracy baseline

   !----- [fast].error_level. -------------------------------------------------------------!
   integer(ik), parameter :: CTRL_L0_FIXED    = 0_ik   !< fixed sub-steps
   integer(ik), parameter :: CTRL_L1_ADAPTIVE = 1_ik   !< adaptive (default)
   integer(ik), parameter :: CTRL_L2_STRICT   = 2_ik   !< adaptive; a miss at the floor step stops the run
   !----- [fast].step_controller. ---------------------------------------------------------!
   integer(ik), parameter :: CTRL_I  = 1_ik   !< integral (default)
   integer(ik), parameter :: CTRL_PI = 2_ik   !< Gustafsson's proportional-integral: damps step-size hunting

   !----- [fast].mask: a masked-off store is held at its start-of-step state while its     !
   !      fluxes still act on the others, so a reduced column costs as much as the full    !
   !      one and does not conserve (mask_is_full gates the budget halts). ----------------!
   type :: process_mask_t
      logical :: veg_energy = .true.   !< leaf + wood energy stores (prognostic modes)
      logical :: cas_energy = .true.   !< canopy-air-space enthalpy (temperature)
      logical :: cas_vapour = .true.   !< canopy-air-space specific humidity
      logical :: cas_co2    = .true.   !< canopy-air-space CO2 twin
      logical :: soil_heat  = .true.   !< soil thermal column
      logical :: soil_water = .true.   !< soil water column (the water network's soil stores)
      logical :: hydraulics = .true.   !< plant hydraulics (psi)
   end type process_mask_t

end module meds_fast_opts
