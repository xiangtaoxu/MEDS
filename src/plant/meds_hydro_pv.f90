!==========================================================================================!
! meds_hydro_pv -- pressure-volume (p-v) kernels: water content and capacitance as functions  !
! of tissue water potential, for the nonlinear symplast model.                               !
!                                                                                          !
! The curve is the classical Tyree & Hammel (1972) turgor + osmotic model, in the form         !
! analysed by Bartlett, Scoffoni & Sack (2012, Ecol. Lett. 15:393), extended with a constant     !
! apoplastic reservoir (Christoffersen et al. 2016, GMD 9:4227). Symplastic relative water        !
! content R in (0,1], R=1 at full turgor. The two fundamental traits are the osmotic potential     !
! at full turgor pi0 [MPa,<0] and the bulk modulus of elasticity eps [MPa,>0]; the turgor-loss     !
! point and RWC at turgor loss follow (Bartlett eqns 1-2):                                          !
!     psi_tlp = pi0*eps/(pi0+eps),     rwc_tlp = (pi0+eps)/eps                                        !
!                                                                                          !
!   turgid  (R >= rwc_tlp):  psi(R) = eps*(R - rwc_tlp) + pi0/R        (turgor + osmotic)             !
!   flaccid (R <  rwc_tlp):  psi(R) = pi0/R                            (turgor lost)                   !
!                                                                                          !
! The turgid branch inverts in closed form (a quadratic in R), so psi->R needs no Newton iteration.  !
! Capacitance C = dW/dpsi = W_sat_sym * dR/dpsi; the apoplast is constant so it drops out of C.         !
! All kernels are elemental / pure and reference only meds_kinds.                                       !
!==========================================================================================!
module meds_hydro_pv
   use meds_kinds, only : wp
   implicit none
   private

   public :: pv_psi_tlp, pv_rwc_tlp, rwc_from_psi, psi_from_rwc
   public :: water_content, capacitance, pv_water_cap_from_traits

   real(wp), parameter :: rwc_floor = 1.0e-4_wp   !< keep R strictly positive in the flaccid tail

contains

   !----- Turgor loss point [MPa] (Bartlett eqn 1). ---------------------------------------!
   elemental real(wp) function pv_psi_tlp(pi0, eps) result(psi_tlp)
      real(wp), intent(in) :: pi0, eps
      psi_tlp = pi0*eps / (pi0 + eps)
   end function pv_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett eqn 2). --------------!
   elemental real(wp) function pv_rwc_tlp(pi0, eps) result(rwc_tlp)
      real(wp), intent(in) :: pi0, eps
      rwc_tlp = (pi0 + eps) / eps
   end function pv_rwc_tlp

   !----- Water potential [MPa] from symplastic RWC. --------------------------------------!
   elemental real(wp) function psi_from_rwc(rwc, pi0, eps) result(psi)
      real(wp), intent(in) :: rwc, pi0, eps
      real(wp) :: rwc_tlp, r
      rwc_tlp = (pi0 + eps) / eps
      r       = max(rwc, rwc_floor)
      if (r >= rwc_tlp) then
         psi = eps*(r - rwc_tlp) + pi0/r      ! turgor + osmotic
      else
         psi = pi0/r                          ! turgor lost
      end if
   end function psi_from_rwc

   !----- Symplastic RWC from water potential (closed-form turgid inverse). ----------------!
   elemental real(wp) function rwc_from_psi(psi, pi0, eps) result(rwc)
      real(wp), intent(in) :: psi, pi0, eps
      real(wp) :: psi_tlp, b, disc
      psi_tlp = pi0*eps / (pi0 + eps)
      if (psi >= psi_tlp) then                ! turgid (psi less negative than the TLP)
         !----- eps*R^2 - (psi+eps+pi0)*R + pi0 = 0; take the + root in [rwc_tlp, 1]. ------!
         b    = psi + eps + pi0
         disc = max(b*b - 4.0_wp*eps*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*eps)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, rwc_floor)
   end function rwc_from_psi

   !----- Total tissue water [kg, per plant] at potential psi. -----------------------------!
   !      W = W_sat_sym*R + W_apoplast, with W_sat = water_sat*biomass, W_sat_sym =           !
   !      (1-af)*W_sat and W_apoplast = af*W_sat (a constant reservoir).                       !
   elemental real(wp) function water_content(psi, pi0, eps, af, water_sat, biomass) result(w)
      real(wp), intent(in) :: psi, pi0, eps, af, water_sat, biomass
      real(wp) :: w_sat, r
      w_sat = water_sat*biomass
      r     = rwc_from_psi(psi, pi0, eps)
      w     = (1.0_wp - af)*w_sat*r + af*w_sat
   end function water_content

   !----- Capacitance C = dW/dpsi [kg/MPa, per plant] at potential psi. ---------------------!
   elemental real(wp) function capacitance(psi, pi0, eps, af, water_sat, biomass) result(c)
      real(wp), intent(in) :: psi, pi0, eps, af, water_sat, biomass
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      psi_tlp   = pi0*eps / (pi0 + eps)
      w_sat_sym = (1.0_wp - af)*water_sat*biomass
      r         = rwc_from_psi(psi, pi0, eps)
      if (psi >= psi_tlp) then
         cr = 1.0_wp / (eps - pi0/(r*r))      ! dR/dpsi, turgid  (-> 1/(eps+|pi0|) at R=1)
      else
         cr = -(r*r)/pi0                       ! dR/dpsi, flaccid (= |pi0|/psi^2)
      end if
      c = w_sat_sym*cr
   end function capacitance

   !----- Legacy-linear mapping: the constant capacitance per biomass [kg/kgC/MPa] that a       !
   !      pure-elastic (X16) reservoir would need to match this curve's full-turgor capacitance. !
   elemental real(wp) function pv_water_cap_from_traits(pi0, eps, af, water_sat) result(water_cap)
      real(wp), intent(in) :: pi0, eps, af, water_sat
      water_cap = (1.0_wp - af)*water_sat / (eps + abs(pi0))
   end function pv_water_cap_from_traits

end module meds_hydro_pv
