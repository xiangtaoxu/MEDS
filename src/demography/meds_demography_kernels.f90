!==========================================================================================!
! meds_demography_kernels -- the dominant per-cohort daily kernels.                                   !
!                                                                                          !
! HARD host/device split (see meds_demography_interface):                                        !
!   * Rate FILL loops call the polymorphic provider -> they run on the HOST and write the   !
!     plain arrays ddbh_dt(:) and a transient total-mortality array.                       !
!   * APPLICATION kernels are `do concurrent` over ALL cohorts of the site, touching only   !
!     plain arrays with branch-light arithmetic -> the compiler maps them to CPU SIMD       !
!     (nvfortran -stdpar=multicore / ifx) or the GPU (nvfortran -stdpar=gpu). No            !
!     polymorphism, no derived-type access, ever reaches the device.                       !
!                                                                                          !
! The competition signal driving growth is the overtopping basal area per ground area      !
! [m2 m-2], computed top-down within each patch (cohorts are kept sorted tallest-first).    !
!==========================================================================================!
module meds_demography_kernels
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : pio4, tiny_num, lnexp_min, lnexp_max
   use meds_demography_types,          only : community
   use meds_demography_interface, only : rate_provider_t, n_mort_cause
   implicit none
   private

   public :: compute_competition, growth_step, mortality_accumulate, mortality_apply_month

   real(wp), parameter :: cm2_to_m2 = 1.0e-4_wp   !< cm2 -> m2 (basarea is cm2/plant)

contains

   !---------------------------------------------------------------------------------------!
   ! Overtopping basal area [m2 m-2] for every cohort, computed per patch from the          !
   ! height-sorted slice (tallest first => index increases downward). HOST; the inner       !
   ! cumulation is sequential within a patch but patches are independent.                   !
   !---------------------------------------------------------------------------------------!
   subroutine compute_competition(comm)
      type(community), intent(inout) :: comm
      integer(ik) :: ip, i0, i1, i
      real(wp)    :: cum
      associate (coh => comm%coh, pat => comm%pat)
         do ip = 1_ik, pat%n
            i0 = pat%coff(ip)
            i1 = i0 + pat%ccount(ip) - 1_ik
            cum = 0.0_wp
            do i = i0, i1
               coh%comp(i) = cum                                   ! BA strictly above this cohort
               cum = cum + coh%nplant(i) * coh%basarea(i) * cm2_to_m2
            end do
         end do
      end associate
   end subroutine compute_competition

   !---------------------------------------------------------------------------------------!
   ! Daily growth: fill rates (host, polymorphic) then advance DBH and re-derive geometry   !
   ! (device kernel, arithmetic only).                                                      !
   !---------------------------------------------------------------------------------------!
   subroutine growth_step(comm, provider, dt_years)
      type(community),        intent(inout) :: comm
      class(rate_provider_t), intent(in)    :: provider
      real(wp),               intent(in)    :: dt_years
      integer(ik) :: i, n

      n = comm%coh%n
      associate (c => comm%coh)
         !----- HOST rate fill (polymorphic dispatch). ------------------------------------!
         do i = 1_ik, n
            c%ddbh_dt(i) = provider%growth(c%pft(i), c%dbh(i), c%hite(i), c%nplant(i), c%comp(i))
         end do
         !----- DEVICE application kernel (offloadable). ----------------------------------!
         do concurrent (i = 1_ik:n)
            c%dbh(i)     = min(c%dbh(i) + c%ddbh_dt(i) * dt_years, c%p_dbh_crit(i))
            c%hite(i)    = c%p_hgt_min(i) + c%p_b1ht(i) * (1.0_wp - exp(-c%p_b2ht(i) * c%dbh(i)))
            c%basarea(i) = pio4 * c%dbh(i) * c%dbh(i)
         end do
      end associate
   end subroutine growth_step

   !---------------------------------------------------------------------------------------!
   ! Daily mortality accumulation: fill total per-capita rate (host, polymorphic) then       !
   ! accumulate log-survival into monthly_dlnndt (device kernel). nplant is NOT changed      !
   ! daily -- it is applied once per month (mirrors ED2 dbalive_dt bookkeeping).             !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_accumulate(comm, provider, dt_years)
      type(community),        intent(inout) :: comm
      class(rate_provider_t), intent(in)    :: provider
      real(wp),               intent(in)    :: dt_years
      real(wp), allocatable :: total_mort(:)
      real(wp)              :: mc(n_mort_cause)
      integer(ik)           :: i, ip, n

      n = comm%coh%n
      allocate(total_mort(n))
      associate (c => comm%coh, p => comm%pat)
         !----- HOST rate fill. -----------------------------------------------------------!
         do i = 1_ik, n
            ip = c%owner_patch(i)
            mc = provider%mortality(c%pft(i), c%dbh(i), c%hite(i), c%nplant(i), c%ddbh_dt(i), &
                                    p%avg_daily_temp(ip), p%age(ip))
            total_mort(i) = sum(mc)
         end do
         !----- DEVICE accumulation kernel. -----------------------------------------------!
         do concurrent (i = 1_ik:n)
            c%monthly_dlnndt(i) = c%monthly_dlnndt(i) - total_mort(i) * dt_years
         end do
      end associate
   end subroutine mortality_accumulate

   !---------------------------------------------------------------------------------------!
   ! Monthly mortality application: apply accumulated survivorship, floored so density does  !
   ! not drop below `negligible` (cohorts that hit the floor are culled by termination).     !
   ! Device kernel (arithmetic only); resets the accumulator.                                !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_apply_month(comm, negligible)
      type(community), intent(inout) :: comm
      real(wp),        intent(in)    :: negligible
      integer(ik) :: i, n
      real(wp)    :: dln, floor_dln

      n = comm%coh%n
      associate (c => comm%coh)
         !----- exp() is inlined (clamped) because non-intrinsic calls cannot offload. ----!
         do concurrent (i = 1_ik:n)
            !----- min(0,..) keeps the floor a DOWNWARD clamp: it must never raise nplant. -!
            floor_dln = min(0.0_wp, log(max(negligible, tiny_num) / max(c%nplant(i), tiny_num)))
            dln       = max(c%monthly_dlnndt(i), floor_dln)
            c%nplant(i)         = c%nplant(i) * exp(min(max(dln, lnexp_min), lnexp_max))
            c%monthly_dlnndt(i) = 0.0_wp
         end do
      end associate
   end subroutine mortality_apply_month

end module meds_demography_kernels
