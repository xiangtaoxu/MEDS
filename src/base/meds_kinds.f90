!==========================================================================================!
! meds_kinds -- single source of numeric kinds for the whole demographic engine.           !
!                                                                                          !
! Retarget the working precision by editing ONE line. `wp` is used for all per-plant and   !
! demographic state (and for area accumulation, which needs the headroom of real64); `sp`  !
! is available for reduced-precision GPU experiments; `ik` is the standard integer kind.    !
!==========================================================================================!
module meds_kinds
   use, intrinsic :: iso_fortran_env, only : real64, real32, int32, int64
   implicit none
   public

   integer, parameter :: wp = real64   !< working precision (CPU kernels + area accumulation)
   integer, parameter :: sp = real32   !< optional reduced precision for GPU throughput
   integer, parameter :: ik = int32    !< index / count integer kind
   integer, parameter :: lk = int64    !< wide integer (step counters)
end module meds_kinds
