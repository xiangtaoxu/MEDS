!==========================================================================================!
! meds_test_support -- shared assert helpers and a configurable analytic rate provider for   !
! the CTest programs. The analytic provider lets a test impose KNOWN growth/mortality/        !
! recruitment so results can be checked against closed-form expectations, and doubles as the  !
! second provider in the provider-swap test.                                                 !
!==========================================================================================!
module meds_test_support
   use meds_kinds,          only : wp, ik
   use meds_rate_interface, only : rate_provider_t, n_mort_cause
   implicit none
   private

   public :: const_provider_t, check, check_close, banner

   !----- Provider with constant growth g, constant mortality m, recruit r above rmin. -----!
   type, extends(rate_provider_t) :: const_provider_t
      real(wp) :: g    = 0.0_wp
      real(wp) :: m    = 0.0_wp
      real(wp) :: r    = 0.0_wp
      real(wp) :: rmin = -1.0e30_wp     !< recruit only if min_month_temp >= rmin
   contains
      procedure :: growth    => cp_growth
      procedure :: mortality => cp_mort
      procedure :: recruit   => cp_recruit
   end type const_provider_t

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

   pure function cp_growth(self, pft, dbh, hite, nplant, comp) result(ddbh)
      class(const_provider_t), intent(in) :: self
      integer(ik),             intent(in) :: pft
      real(wp),                intent(in) :: dbh, hite, nplant, comp
      real(wp)                            :: ddbh
      ddbh = self%g
   end function cp_growth

   pure function cp_mort(self, pft, dbh, hite, nplant, growth_pyr, avg_daily_temp,         &
                         patch_age) result(mrate)
      class(const_provider_t), intent(in) :: self
      integer(ik),             intent(in) :: pft
      real(wp),                intent(in) :: dbh, hite, nplant, growth_pyr, avg_daily_temp, &
                                             patch_age
      real(wp)                            :: mrate(n_mort_cause)
      mrate    = 0.0_wp
      mrate(1) = self%m
   end function cp_mort

   pure function cp_recruit(self, pft, dist_type, min_month_temp) result(rdens)
      class(const_provider_t), intent(in) :: self
      integer(ik),             intent(in) :: pft, dist_type
      real(wp),                intent(in) :: min_month_temp
      real(wp)                            :: rdens
      if (min_month_temp >= self%rmin) then
         rdens = self%r
      else
         rdens = 0.0_wp
      end if
   end function cp_recruit

end module meds_test_support
