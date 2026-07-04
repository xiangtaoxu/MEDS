!==========================================================================================!
! meds_constants -- compile-time constants and a few pure scalar helpers.                  !
!                                                                                          !
! No mutable state lives here. Calendar constants follow ED2 (proleptic Gregorian mean     !
! year). Tolerances default to the values used by the fusion/fission conservation checks.   !
!==========================================================================================!
module meds_constants
   use meds_kinds, only : wp
   implicit none
   public

   !----- Geometry. -----------------------------------------------------------------------!
   real(wp), parameter :: pi   = 3.14159265358979323846_wp
   real(wp), parameter :: pio4 = 0.25_wp * pi                 !< basal area = pio4 * dbh^2

   !----- Calendar (mean Gregorian year). -------------------------------------------------!
   real(wp), parameter :: yr_day      = 365.2425_wp           !< days per year
   real(wp), parameter :: day_sec     = 86400.0_wp            !< seconds per day
   real(wp), parameter :: yr_sec      = yr_day * day_sec      !< seconds per year
   real(wp), parameter :: mon_per_yr  = 12.0_wp
   real(wp), parameter :: day_per_mon = yr_day / mon_per_yr   !< ~30.4 days per month

   !----- Physical constants (leaf biophysics: photosynthesis / stomatal conductance). ----!
   real(wp), parameter :: r_gas       = 8.314462618_wp       !< [J/mol/K] universal gas constant
   real(wp), parameter :: t_kelvin    = 273.15_wp            !< [K] 0 degC in Kelvin
   real(wp), parameter :: p_std       = 101325.0_wp          !< [Pa] standard atmospheric pressure
   real(wp), parameter :: t_ref_photo = 298.15_wp            !< [K] 25 degC reference for photosynthesis
   real(wp), parameter :: stefan      = 5.670374419e-8_wp    !< [W/m2/K4] Stefan-Boltzmann (canopy RT longwave)
   real(wp), parameter :: halfpi      = 0.5_wp * pi          !< [rad] pi/2 (leaf-inclination domain upper bound)
   real(wp), parameter :: grav_head   = 9.804e-3_wp          !< [MPa/m] hydrostatic head (rho_w*g); plant hydraulics

   !----- Numerical safety. ---------------------------------------------------------------!
   real(wp), parameter :: tiny_num   = 1.0e-30_wp
   real(wp), parameter :: almost_one = 1.0_wp - 1.0e-7_wp
   real(wp), parameter :: lnexp_min  = -38.0_wp               !< clamp for exp() arguments
   real(wp), parameter :: lnexp_max  =  38.0_wp

   !----- Default conservation tolerances. ------------------------------------------------!
   real(wp), parameter :: size_tol = 0.01_wp                  !< 1% basal-area conservation
   real(wp), parameter :: area_tol = 1.0e-5_wp                !< patch-area sum tolerance

contains

   !---------------------------------------------------------------------------------------!
   ! Numerically safe exponential: clamps the argument to avoid under/overflow.            !
   !---------------------------------------------------------------------------------------!
   elemental pure function safe_exp(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = exp(min(max(x, lnexp_min), lnexp_max))
   end function safe_exp

end module meds_constants
