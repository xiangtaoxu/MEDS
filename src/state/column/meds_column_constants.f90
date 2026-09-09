!==========================================================================================!
! meds_column_constants -- the compile-time DIMENSIONS of a soil/snow column and the         !
! fresh-cohort INITIAL values. Split out of the former meds_column_state_types so that the   !
! reservoirs and the parameter bundles can each depend on the dimensions independently.      !
!                                                                                          !
! Fixed sizes are what keep the columns allocatable-free and GPU-eligible; the init values   !
! are the state a cohort is BORN with, which is neither a reservoir nor a parameter.          !
!==========================================================================================!
module meds_column_constants
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: n_soil_layer_max, n_snow_layer_max, N_HYDRO_NODE, LEAF_TEMP_INIT, PSI_INIT

   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time soil-column-depth ceiling
   integer(ik), parameter :: n_snow_layer_max = 1_ik        !< MVP single bulk layer; P1 raises to ED2 nzs~8

   !----- Per-cohort fast state carried on cohort_block (rides the cohort lockstep). N_HYDRO_NODE !
   !      MUST equal meds_plant_types%N_HYDRO (the psi node count); a fresh cohort starts at a     !
   !      mild tension + a neutral leaf temperature (both relax within one fast step).            !
   integer(ik), parameter :: N_HYDRO_NODE   = 3_ik          !< == N_HYDRO (leaf/wood/root psi nodes)
   real(wp),    parameter :: LEAF_TEMP_INIT = 288.15_wp     !< [K]   fresh-cohort leaf temperature
   real(wp),    parameter :: PSI_INIT       = -0.1_wp       !< [MPa] fresh-cohort node water potential

end module meds_column_constants
