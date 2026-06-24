!==========================================================================================!
! meds_demography_kernels -- the dominant per-cohort daily kernels.                        !
!                                                                                          !
! The demographic RATES arrive precomputed from outside the engine (see update_demography), !
! so these are now PURE ARRAY kernels: branch-light arithmetic over plain 1-D arrays with no !
! polymorphism and no host-side rate evaluation. They are written as `do concurrent` so the  !
! compiler maps them to CPU SIMD (nvfortran -stdpar=multicore / ifx) or the GPU (nvfortran   !
! -stdpar=gpu) with no fallback to the host.                                                 !
!==========================================================================================!
module meds_demography_kernels
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : pio4, tiny_num, lnexp_min, lnexp_max
   use meds_demography_types,  only : site
   implicit none
   private

   public :: growth_step, mortality_accumulate, mortality_apply_month

contains

   !---------------------------------------------------------------------------------------!
   ! Daily growth: advance DBH by the supplied per-cohort growth rate and re-derive the     !
   ! cached geometry. `growth` is indexed in the current cohort order (1:site%coh%n).        !
   !---------------------------------------------------------------------------------------!
   subroutine growth_step(asite, growth, dt_years)
      type(site), intent(inout) :: asite
      real(wp),   intent(in)    :: growth(:)        !< [cm/yr] per cohort
      real(wp),   intent(in)    :: dt_years
      integer(ik) :: i, n
      n = asite%coh%n
      associate (c => asite%coh)
         do concurrent (i = 1_ik:n)
            c%dbh(i)     = min(c%dbh(i) + growth(i) * dt_years, c%p_dbh_crit(i))
            c%hite(i)    = c%p_hgt_min(i) + c%p_b1ht(i) * (1.0_wp - exp(-c%p_b2ht(i) * c%dbh(i)))
            c%basarea(i) = pio4 * c%dbh(i) * c%dbh(i)
         end do
      end associate
   end subroutine growth_step

   !---------------------------------------------------------------------------------------!
   ! Daily mortality accumulation: integrate the supplied per-cohort total mortality rate    !
   ! into the monthly log-survival accumulator. nplant is NOT changed daily (applied once a   !
   ! month, mirroring ED2 bookkeeping).                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_accumulate(asite, mortality, dt_years)
      type(site), intent(inout) :: asite
      real(wp),   intent(in)    :: mortality(:)     !< [1/yr] per cohort, total
      real(wp),   intent(in)    :: dt_years
      integer(ik) :: i, n
      n = asite%coh%n
      associate (c => asite%coh)
         do concurrent (i = 1_ik:n)
            c%monthly_dlnndt(i) = c%monthly_dlnndt(i) - mortality(i) * dt_years
         end do
      end associate
   end subroutine mortality_accumulate

   !---------------------------------------------------------------------------------------!
   ! Monthly mortality application: apply accumulated survivorship, floored (downward only)  !
   ! so density never rises and never drops below `negligible`. Resets the accumulator.       !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_apply_month(asite, negligible)
      type(site), intent(inout) :: asite
      real(wp),   intent(in)    :: negligible
      integer(ik) :: i, n
      real(wp)    :: dln, floor_dln
      n = asite%coh%n
      associate (c => asite%coh)
         do concurrent (i = 1_ik:n)
            floor_dln = min(0.0_wp, log(max(negligible, tiny_num) / max(c%nplant(i), tiny_num)))
            dln       = max(c%monthly_dlnndt(i), floor_dln)
            c%nplant(i)         = c%nplant(i) * exp(min(max(dln, lnexp_min), lnexp_max))
            c%monthly_dlnndt(i) = 0.0_wp
         end do
      end associate
   end subroutine mortality_apply_month

end module meds_demography_kernels
