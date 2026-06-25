!==========================================================================================!
! meds_demography_kernels -- the dominant per-cohort daily kernels.                        !
!                                                                                          !
! The demographic RATES arrive precomputed from outside the engine (see update_demography), !
! so these are PURE ARRAY kernels: branch-light arithmetic over plain 1-D arrays, with NO    !
! derived types and NO polymorphism. Each is wrapped in an explicit OpenMP `target` region    !
! with `map` clauses, so the data is copied to/from the device per call and the host keeps    !
! ALL state in normal memory. (This deliberately avoids `-stdpar=gpu`'s global CUDA-managed   !
! allocator, whose deep-copy/finalize of the allocatable-component `site` type double-frees.) !
!                                                                                          !
! One directive serves three backends via the build flag: nvfortran `-mp=gpu` -> GPU,        !
! `-mp` -> CPU threads, no flag -> the `!$omp` lines are comments and the loop runs serially. !
!==========================================================================================!
module meds_demography_kernels
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pio4, tiny_num, lnexp_min, lnexp_max
   implicit none
   private

   public :: growth_step, mortality_step

contains

   !---------------------------------------------------------------------------------------!
   ! Daily growth: advance DBH by the supplied per-cohort growth rate and re-derive the     !
   ! cached geometry. Arrays are indexed in the current cohort order (1:n).                  !
   !---------------------------------------------------------------------------------------!
   subroutine growth_step(n, dbh, height, basal_area, dbh_crit, height_min, b1_height,      &
                          b2_height, growth, dt_yr)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: dbh(:), height(:), basal_area(:)
      real(wp),    intent(in)    :: dbh_crit(:), height_min(:), b1_height(:), b2_height(:)
      real(wp),    intent(in)    :: growth(:)
      real(wp),    intent(in)    :: dt_yr
      integer(ik) :: i

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: dbh_crit, height_min, b1_height, b2_height, growth)              &
      !$omp&        map(tofrom: dbh, height, basal_area)
      do i = 1_ik, n
         dbh(i)        = min(dbh(i) + growth(i) * dt_yr, dbh_crit(i))
         height(i)     = height_min(i) + b1_height(i) * (1.0_wp - exp(-b2_height(i) * dbh(i)))
         basal_area(i) = pio4 * dbh(i) * dbh(i)
      end do
   end subroutine growth_step

   !---------------------------------------------------------------------------------------!
   ! Per-step mortality: apply the supplied per-cohort total mortality rate directly to       !
   ! nplant for one step of length dt_yr, floored (downward only) so density never rises and  !
   ! never drops below `negligible`.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_step(n, nplant, mortality, dt_yr, negligible)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: nplant(:)
      real(wp),    intent(in)    :: mortality(:)
      real(wp),    intent(in)    :: dt_yr, negligible
      integer(ik) :: i
      real(wp)    :: dln, floor_dln

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: mortality) map(tofrom: nplant) private(dln, floor_dln)
      do i = 1_ik, n
         dln         = -mortality(i) * dt_yr
         floor_dln   = min(0.0_wp, log(max(negligible, tiny_num) / max(nplant(i), tiny_num)))
         dln         = max(dln, floor_dln)
         nplant(i)   = nplant(i) * exp(min(max(dln, lnexp_min), lnexp_max))
      end do
   end subroutine mortality_step

end module meds_demography_kernels
