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
   public :: process_mask_t

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

   !----- [fast].mask: the uniform PROCESS MASK (MEDS_NUMERICS_SCOPING.md §5.1). One switch set that every scheme     !
   !      honors, so the SAME driver can run a REDUCED column ODE: turning off energy leaves a water +   !
   !      CO2 column, turning off all but hydraulics isolates the stiffest mode. It is the process-      !
   !      complexity axis of a numerics sweep, and it is scientifically useful on its own (attribute     !
   !      column behaviour to a process).                                                                !
   !                                                                                                     !
   !      SEMANTICS: a masked-OFF process is FROZEN -- its store does not evolve, so the ODE loses that  !
   !      dimension while the process still supplies its couplings to the others as a CONSTANT. The      !
   !      kernel is still invoked and its store restored afterwards (rather than the call being skipped) !
   !      so no downstream consumer is ever handed an unset flux; the cost of a reduced system is        !
   !      therefore NOT lower, which matters when reading harness WORK metrics.                          !
   !                                                                                                     !
   !      CONSERVATION: a reduced column is deliberately NOT closed -- freezing a store while its fluxes !
   !      still act on its neighbours breaks the ledger by construction. mask_is_full() reports whether  !
   !      the budget halts are meaningful, and the driver suppresses them when they are not.             !
   !      All-true (the default) is the full column.                                                     !
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
