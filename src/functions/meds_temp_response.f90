!==========================================================================================!
! meds_temp_response -- temperature scaling of the photosynthetic biochemistry.       !
!                                                                                          !
! Two forms, both relative to a 25 degC reference (t_ref_photo):                            !
!   * Arrhenius:  k(T) = k25 * exp[ Ea/(R*Tref) * (1 - Tref/T) ]                            !
!   * Peaked (Medlyn et al. 2002 / FATES): the Arrhenius rise multiplied by a high-          !
!     temperature deactivation envelope, normalised to 1 at the reference:                   !
!         k(T) = k_arr(T) * fH(Tref)/fH(T),   fH(T) = 1 + exp[ (dS*T - Hd)/(R*T) ]           !
!     which produces a thermal optimum and a roll-off above it.                              !
!                                                                                          !
! The Michaelis constants Kc, Ko and the CO2 compensation point Gamma* always use the plain  !
! Arrhenius form (Bernacchi et al. 2001 activation energies); Vcmax, Jmax and Rd use the      !
! caller-selected form. All routines are pure elemental and use the clamped safe_exp.         !
!==========================================================================================!
module meds_temp_response
   use meds_kinds,     only : wp, ik
   use meds_constants, only : r_gas, t_ref_photo, safe_exp
   implicit none
   private

   public :: arrhenius_scale, peaked_arrhenius_scale, temp_response
   public :: TRESP_ARRHENIUS, TRESP_PEAKED

   !----- Temperature-response form selectors (owned here; re-exported by meds_config). ----!
   integer(ik), parameter :: TRESP_ARRHENIUS = 1_ik  !< plain Arrhenius
   integer(ik), parameter :: TRESP_PEAKED    = 2_ik  !< Arrhenius with high-temperature deactivation

contains

   !---------------------------------------------------------------------------------------!
   ! Plain Arrhenius scaling of a rate constant from its 25 degC value to leaf temperature. !
   !---------------------------------------------------------------------------------------!
   elemental pure function arrhenius_scale(k25, ea, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol] activation energy
      real(wp), intent(in) :: t_leaf    !< [K] leaf temperature
      real(wp)             :: k
      k = k25 * safe_exp(ea / (r_gas * t_ref_photo) * (1.0_wp - t_ref_photo / t_leaf))
   end function arrhenius_scale

   !---------------------------------------------------------------------------------------!
   ! Peaked Arrhenius: the Arrhenius rise with a high-temperature deactivation envelope,    !
   ! normalised to the reference temperature so k(Tref) = k25.                              !
   !---------------------------------------------------------------------------------------!
   elemental pure function peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol]   activation energy
      real(wp), intent(in) :: hd        !< [J/mol]   deactivation energy
      real(wp), intent(in) :: ds        !< [J/mol/K] entropy term
      real(wp), intent(in) :: t_leaf    !< [K] leaf temperature
      real(wp)             :: k, fh_ref, fh_leaf
      fh_ref  = 1.0_wp + safe_exp((ds * t_ref_photo - hd) / (r_gas * t_ref_photo))
      fh_leaf = 1.0_wp + safe_exp((ds * t_leaf      - hd) / (r_gas * t_leaf))
      k = arrhenius_scale(k25, ea, t_leaf) * fh_ref / fh_leaf
   end function peaked_arrhenius_scale

   !---------------------------------------------------------------------------------------!
   ! Dispatch on the configured form (TRESP_ARRHENIUS uses Ea only; TRESP_PEAKED adds the   !
   ! Hd/dS deactivation). The deactivation arguments are ignored for the Arrhenius form.    !
   !---------------------------------------------------------------------------------------!
   elemental pure function temp_response(form, k25, ea, hd, ds, t_leaf) result(k)
      integer(ik), intent(in) :: form   !< TRESP_ARRHENIUS | TRESP_PEAKED
      real(wp),    intent(in) :: k25, ea, hd, ds, t_leaf
      real(wp)                :: k
      select case (form)
      case (TRESP_PEAKED) ; k = peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf)
      case default        ; k = arrhenius_scale(k25, ea, t_leaf)
      end select
   end function temp_response

end module meds_temp_response
