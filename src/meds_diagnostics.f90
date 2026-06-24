!==========================================================================================!
! meds_diagnostics -- pure, side-effect-free reductions over the community state.          !
!                                                                                          !
! Site-level totals weight per-patch quantities by patch area. These are used by the demo   !
! and the test suite to check conservation, detect NaNs, and report the evolving stand      !
! structure.                                                                               !
!==========================================================================================!
module meds_diagnostics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_types,     only : community
   implicit none
   private

   public :: total_nplant, total_basal_area, total_area, mean_dbh, count_cohorts, has_nan
   public :: print_summary

   real(wp), parameter :: cm2_to_m2 = 1.0e-4_wp

contains

   !----- Site plant number [plant m-2 ground]. -------------------------------------------!
   pure real(wp) function total_nplant(comm) result(tot)
      type(community), intent(in) :: comm
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, comm%pat%n
         i0 = comm%pat%coff(ip) ; i1 = i0 + comm%pat%ccount(ip) - 1_ik
         if (i1 >= i0) tot = tot + comm%pat%area(ip) * sum(comm%coh%nplant(i0:i1))
      end do
   end function total_nplant

   !----- Site basal area [m2 m-2]. -------------------------------------------------------!
   pure real(wp) function total_basal_area(comm) result(tot)
      type(community), intent(in) :: comm
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, comm%pat%n
         i0 = comm%pat%coff(ip) ; i1 = i0 + comm%pat%ccount(ip) - 1_ik
         if (i1 >= i0) tot = tot + comm%pat%area(ip)                                        &
                            * sum(comm%coh%nplant(i0:i1) * comm%coh%basarea(i0:i1)) * cm2_to_m2
      end do
   end function total_basal_area

   pure real(wp) function total_area(comm) result(tot)
      type(community), intent(in) :: comm
      tot = sum(comm%pat%area(1:comm%pat%n))
   end function total_area

   !----- Basal-area-weighted mean diameter [cm] (a simple structural index). -------------!
   pure real(wp) function mean_dbh(comm) result(dm)
      type(community), intent(in) :: comm
      real(wp)    :: wsum, dsum, w
      integer(ik) :: ip, i0, i1, i
      wsum = 0.0_wp ; dsum = 0.0_wp
      do ip = 1_ik, comm%pat%n
         i0 = comm%pat%coff(ip) ; i1 = i0 + comm%pat%ccount(ip) - 1_ik
         do i = i0, i1
            w    = comm%pat%area(ip) * comm%coh%nplant(i) * comm%coh%basarea(i)
            wsum = wsum + w
            dsum = dsum + w * comm%coh%dbh(i)
         end do
      end do
      dm = dsum / max(wsum, tiny_num)
   end function mean_dbh

   pure integer(ik) function count_cohorts(comm) result(nc)
      type(community), intent(in) :: comm
      nc = comm%coh%n
   end function count_cohorts

   !----- True if any diameter or density is non-finite. ----------------------------------!
   pure logical function has_nan(comm) result(bad)
      type(community), intent(in) :: comm
      integer(ik) :: i
      bad = .false.
      do i = 1_ik, comm%coh%n
         if (.not. (comm%coh%dbh(i) == comm%coh%dbh(i))) bad = .true.
         if (.not. (comm%coh%nplant(i) == comm%coh%nplant(i))) bad = .true.
      end do
   end function has_nan

   !----- Human-readable one-line summary. ------------------------------------------------!
   subroutine print_summary(comm, label)
      type(community),  intent(in) :: comm
      character(len=*), intent(in) :: label
      write(*,'(a,t18,a,i6,a,i4,a,f9.4,a,f10.5,a,f7.2)')                                   &
         trim(label), 'cohorts=', comm%coh%n, '  patches=', comm%pat%n,                    &
         '  N=', total_nplant(comm), '  BA=', total_basal_area(comm),                       &
         '  Dmean=', mean_dbh(comm)
   end subroutine print_summary

end module meds_diagnostics
