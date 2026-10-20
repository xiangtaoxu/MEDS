! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_temp_response -- how a rate scales from its 25 degC value (t_ref_photo) to the      !
! tissue temperature: the plain Arrhenius form, and the peaked form, which adds a high-    !
! temperature deactivation (Medlyn et al. 2002; docs/science/plant_respiration.md          !
! eq. (5)). Leaf gas exchange and stem and root respiration all scale through these        !
! routines, so the tissues share one thermal shape (docs/science/leaf_gas_exchange.md §2   !
! says which rate takes which form). Also the Kattge & Knorr (2007) growth-temperature     !
! acclimation of the peaked form (docs/science/leaf_gas_exchange.md eq. (4)).              !
!==========================================================================================!
module meds_temp_response
   use meds_kinds,     only : wp, ik
   use meds_constants, only : r_gas, t_ref_photo, safe_exp, t_kelvin
   implicit none
   private

   public :: arrhenius_scale, peaked_arrhenius_scale, temp_response
   public :: kattge_knorr_entropy, kattge_knorr_jv_ratio
   public :: TRESP_ARRHENIUS, TRESP_PEAKED

   !----- The forms [leaf_physiology].temp_response_form selects. -------------------------!
   integer(ik), parameter :: TRESP_ARRHENIUS = 1_ik  !< plain Arrhenius
   integer(ik), parameter :: TRESP_PEAKED    = 2_ik  !< Arrhenius with high-temperature deactivation

contains

   !----- Plain Arrhenius: k25 exp[Ea/(R Tref) (1 - Tref/T)], Tref = t_ref_photo. ---------!
   elemental pure function arrhenius_scale(k25, ea, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol] activation energy
      real(wp), intent(in) :: t_leaf    !< [K] tissue temperature
      real(wp)             :: k
      k = k25 * safe_exp(ea / (r_gas * t_ref_photo) * (1.0_wp - t_ref_photo / t_leaf))
   end function arrhenius_scale

   !----- Peaked Arrhenius: the Arrhenius rise times a high-temperature deactivation,      !
   !      normalised so k(t_ref_photo) = k25. ---------------------------------------------!
   elemental pure function peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf) result(k)
      real(wp), intent(in) :: k25       !< value at the reference temperature
      real(wp), intent(in) :: ea        !< [J/mol]   activation energy
      real(wp), intent(in) :: hd        !< [J/mol]   deactivation energy
      real(wp), intent(in) :: ds        !< [J/mol/K] entropy term
      real(wp), intent(in) :: t_leaf    !< [K] tissue temperature
      real(wp)             :: k, fh_ref, fh_leaf
      fh_ref  = 1.0_wp + safe_exp((ds * t_ref_photo - hd) / (r_gas * t_ref_photo))
      fh_leaf = 1.0_wp + safe_exp((ds * t_leaf      - hd) / (r_gas * t_leaf))
      k = arrhenius_scale(k25, ea, t_leaf) * fh_ref / fh_leaf
   end function peaked_arrhenius_scale

   !----- Kattge & Knorr (2007) acclimation of the peaked form's entropy term,             !
   !      dS = a - b T_growth[degC] (docs/science/leaf_gas_exchange.md eq. (4)), with      !
   !      T_growth a running mean of air temperature, the quantity the fit was made on. ---!
   elemental pure function kattge_knorr_entropy(a_coef, b_coef, t_growth) result(ds)
      real(wp), intent(in) :: a_coef, b_coef   !< [J/mol/K], [J/mol/K2] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: ds
      ds = a_coef - b_coef * (t_growth - t_kelvin)
   end function kattge_knorr_entropy

   !----- Kattge & Knorr (2007) acclimation of Jmax25/Vcmax25 = a - b T_growth[degC]. -----!
   elemental pure function kattge_knorr_jv_ratio(a_coef, b_coef, t_growth) result(jv)
      real(wp), intent(in) :: a_coef, b_coef   !< [-], [1/K] intercept and slope
      real(wp), intent(in) :: t_growth         !< [K] growth temperature
      real(wp)             :: jv
      jv = max(0.1_wp, a_coef - b_coef * (t_growth - t_kelvin))   ! clamp-ok: positive far outside the 11-35 C fit
   end function kattge_knorr_jv_ratio

   !----- The configured form; the Arrhenius form ignores hd and ds. ----------------------!
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
