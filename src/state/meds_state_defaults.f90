! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_state_defaults -- the values a state field holds before it has a real one, in one   !
! place: the "no value yet" markers the model tests for, a new cohort's phenology at       !
! birth, and the canopy-air placeholders a patch holds until the loops first set them. A   !
! default of 0 (or of 1 for a liquid fraction) is written where it is used, not named.     !
!==========================================================================================!
module meds_state_defaults
   use meds_kinds, only : wp
   implicit none
   private

   !----- Cohort: "no value yet" markers, each outside its field's real range. --------------!
   real(wp), parameter, public :: GROWTH_AVG_UNSET          = -1.0_wp     !< [cm/yr] before the first growth sample (a real one is >= 0)
   real(wp), parameter, public :: DMAX_PSI_LEAF_UNSET       =  1.0_wp     !< [MPa] no predawn history (a real one is <= 0); the fast loop seeds it
   real(wp), parameter, public :: DMAX_PSI_LEAF_ACCUM_RESET = -1.0e30_wp  !< [MPa] start of the daily max: below any potential, so max() replaces it
   real(wp), parameter, public :: LIGHT_HOURS_MEAN_UNSET    = -1.0_wp     !< [h/day] no light memory; the cohort's first day sets it

   !----- Cohort: a newborn's phenology -- flushing, not shedding. ---------------------------!
   real(wp), parameter, public :: PHENO_FLUSH_INIT = 1.0_wp
   real(wp), parameter, public :: PHENO_SHED_INIT  = 0.0_wp

   !----- Patch: the canopy air until the fast loop seeds it and the slow loop sizes it. -----!
   real(wp), parameter, public :: CAN_CO2_UNSEEDED = 400.0_wp   !< [umol/mol]
   real(wp), parameter, public :: CAN_DEPTH_UNSET  =  20.0_wp   !< [m]

end module meds_state_defaults
