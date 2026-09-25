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
   use meds_constants, only : r_gas, t_ref_photo, safe_exp, t_kelvin
   implicit none
   private

   public :: arrhenius_scale, peaked_arrhenius_scale, temp_response
   public :: kattge_knorr_entropy, kattge_knorr_jv_ratio
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
   ! THERMAL ACCLIMATION, Kattge & Knorr (2007) -- issue #176.                                 !
   !                                                                                          !
   ! The peaked form's entropy term dS sets where the response PEAKS: raising dS moves the     !
   ! optimum up. Kattge & Knorr fit dS as a falling linear function of the GROWTH temperature   !
   ! (the mean air temperature of the preceding weeks), so a warm-grown plant runs a higher     !
   ! optimum than a cold-grown one of the same PFT:                                             !
   !                                                                                          !
   !     dS(T_growth) = a - b * T_growth[degC]                                                  !
   !                                                                                          !
   ! with a = 668.39, b = 1.07 for Vcmax and a = 659.70, b = 0.75 for Jmax. The SAME fit also   !
   ! acclimates the capacity RATIO, Jmax25/Vcmax25 = 2.59 - 0.035*T_growth -- a warm-grown      !
   ! plant invests relatively less in electron transport. Applying the dS shift without the     !
   ! ratio would acclimate the shape of the response while leaving its two branches in a fixed  !
   ! proportion, which is not what the study measured.                                          !
   !                                                                                          !
   ! T_growth is the GROWTH temperature, not the leaf temperature: the fit is calibrated on the !
   ! mean AIR temperature of the preceding month, so feeding it an instantaneous leaf value     !
   ! would use the relation far outside what it was fitted to.                                  !
   !---------------------------------------------------------------------------------------!
   elemental pure function kattge_knorr_entropy(a_coef, b_coef, t_growth) result(ds)
      real(wp), intent(in) :: a_coef, b_coef   !< [J/mol/K], [J/mol/K2] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: ds
      ds = a_coef - b_coef * (t_growth - t_kelvin)
   end function kattge_knorr_entropy

   elemental pure function kattge_knorr_jv_ratio(a_coef, b_coef, t_growth) result(jv)
      real(wp), intent(in) :: a_coef, b_coef   !< [-], [1/K] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: jv
      jv = max(0.1_wp, a_coef - b_coef * (t_growth - t_kelvin))
   end function kattge_knorr_jv_ratio

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
