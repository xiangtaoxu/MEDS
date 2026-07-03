!==========================================================================================!
! meds_leaf_stomata -- stomatal-conductance models.                                        !
!                                                                                          !
! Two semi-empirical gs(A) laws return stomatal conductance to water vapour [mol H2O/m2/s]   !
! from net assimilation and the leaf environment; the third (Katul) is an OPTIMIZATION model  !
! that fixes intercellular CO2 directly, so it has no gs(A) law -- the solver finds the         !
! optimal Ci and back-computes gs. This module supplies the two empirical laws plus the         !
! water-stress-scaled marginal water-use efficiency lambda that the Katul residual needs.       !
!                                                                                          !
!   Leuning (1995):  gs = g0 + g1 A / [ (Cs - Gamma*) (1 + VPD/D0) ]                         !
!   Medlyn (2011):   gs = g0 + 1.6 (1 + g1/sqrt(VPD_kPa)) A / Cs                             !
!                                                                                          !
! Net assimilation A and the CO2 mole fractions Cs, Gamma* share units [umol/mol]; VPD is in   !
! Pa (converted to kPa internally for Medlyn). When A <= 0 (respiring leaf) both laws collapse  !
! to the cuticular conductance g0.                                                             !
!==========================================================================================!
module meds_leaf_stomata
   use meds_kinds,     only : wp
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: stomata_gs_leuning, stomata_gs_medlyn, katul_lambda

   real(wp), parameter :: vpd_floor_pa = 50.0_wp     !< [Pa] VPD floor (avoid 1/sqrt(0) in Medlyn)
   real(wp), parameter :: beta_floor   = 1.0e-4_wp   !< [--] water-stress floor (bound lambda as beta->0)

contains

   !---------------------------------------------------------------------------------------!
   ! Leuning (1995) semi-empirical stomatal conductance.                                   !
   !---------------------------------------------------------------------------------------!
   pure function stomata_gs_leuning(a_net, cs, gstar, vpd, g0, g1, d0) result(gs)
      real(wp), intent(in) :: a_net   !< [umol/m2/s] net assimilation
      real(wp), intent(in) :: cs      !< [umol/mol]  leaf-surface CO2
      real(wp), intent(in) :: gstar   !< [umol/mol]  CO2 compensation point
      real(wp), intent(in) :: vpd     !< [Pa]        leaf-to-air VPD
      real(wp), intent(in) :: g0, g1  !< [mol/m2/s], [--]
      real(wp), intent(in) :: d0      !< [Pa]        humidity sensitivity
      real(wp)             :: gs, denom
      if (a_net <= 0.0_wp) then
         gs = g0 ; return
      end if
      denom = max(cs - gstar, tiny_num) * (1.0_wp + vpd / d0)
      gs = g0 + g1 * a_net / denom
   end function stomata_gs_leuning

   !---------------------------------------------------------------------------------------!
   ! Medlyn et al. (2011) unified stomatal optimization (USO) conductance.                 !
   !---------------------------------------------------------------------------------------!
   pure function stomata_gs_medlyn(a_net, cs, vpd, g0, g1) result(gs)
      real(wp), intent(in) :: a_net   !< [umol/m2/s] net assimilation
      real(wp), intent(in) :: cs      !< [umol/mol]  leaf-surface CO2
      real(wp), intent(in) :: vpd     !< [Pa]        leaf-to-air VPD
      real(wp), intent(in) :: g0, g1  !< [mol/m2/s], [kPa^0.5]
      real(wp)             :: gs, vpd_kpa
      if (a_net <= 0.0_wp) then
         gs = g0 ; return
      end if
      vpd_kpa = max(vpd, vpd_floor_pa) * 1.0e-3_wp
      gs = g0 + 1.6_wp * (1.0_wp + g1 / sqrt(vpd_kpa)) * a_net / max(cs, tiny_num)
   end function stomata_gs_medlyn

   !---------------------------------------------------------------------------------------!
   ! Katul marginal water-use efficiency lambda [umol CO2/mol H2O], scaled by the water-     !
   ! stress factor beta in (0,1]: drier leaves carry a larger marginal water cost, so lambda  !
   ! rises (lambda = lambda25 * beta^(-exp)) and the optimal stomata close. exp = 0 disables.  !
   !---------------------------------------------------------------------------------------!
   pure function katul_lambda(lambda25, beta, lambda_psi_exp) result(lambda)
      real(wp), intent(in) :: lambda25        !< [umol CO2/mol H2O] well-watered marginal WUE
      real(wp), intent(in) :: beta            !< [--] water-stress factor (0 = closed, 1 = open)
      real(wp), intent(in) :: lambda_psi_exp  !< [--] water-stress exponent
      real(wp)             :: lambda
      lambda = lambda25 * max(beta, beta_floor) ** (-lambda_psi_exp)
   end function katul_lambda

end module meds_leaf_stomata
