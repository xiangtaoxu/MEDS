!==========================================================================================!
! meds_test_support -- shared assert helpers for the CTest programs.                       !
!==========================================================================================!
module meds_test_support
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: check, check_close, banner

contains

   subroutine check(cond, msg)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         write(*,'(2a)') 'FAIL: ', msg
         error stop 1
      end if
   end subroutine check

   subroutine check_close(a, b, rtol, msg)
      real(wp),         intent(in) :: a, b, rtol
      character(len=*), intent(in) :: msg
      real(wp) :: tol
      tol = rtol * max(abs(b), 1.0e-30_wp) + 1.0e-12_wp
      if (abs(a - b) > tol) then
         write(*,'(2a)') 'FAIL: ', msg
         write(*,'(a,es16.8,a,es16.8,a,es10.2)') '   got=', a, ' expected=', b, ' rtol=', rtol
         error stop 1
      end if
   end subroutine check_close

   subroutine banner(name)
      character(len=*), intent(in) :: name
      write(*,'(2a)') '[test] ', name
   end subroutine banner

end module meds_test_support
