!==========================================================================================!
! meds_hydro_conductance -- xylem vulnerability and the Kirchhoff-integrated conductance.    !
!                                                                                          !
! The pointwise loss of conductance is the standard sigmoid                                  !
!     plc_retained(psi) = 1 / (1 + (psi/psi50)^kexp)      (fraction of conductance retained)   !
! but it is NEVER sampled at a single potential to form an edge conductance. Instead every       !
! edge uses the Kirchhoff (matric flux) potential                                                !
!     Phi(psi) = integral_0^psi plc_retained(s) ds       [MPa]                                     !
! so the exact steady flux across a storage-free edge is q = k_cond * (Phi(up) - Phi(down)), and    !
! the equivalent secant conductance is                                                               !
!     K_eff = k_cond * (Phi(up) - Phi(down)) / (up - down) = k_cond * <plc>                            !
! i.e. the potential-AVERAGED retained fraction over the drop. This is more physically exact for a     !
! storage-free edge than a pointwise value, and -- because the flux is monotone in psi with a bounded  !
! Jacobian d(q)/d(psi) = +/- k_cond*plc -- it is far better conditioned near cavitation. Phi is closed   !
! form for kexp in {1,2} (log / arctan) and a fixed Gauss-Legendre quadrature otherwise; its inverse     !
! (needed by the vertical-profile diagnostic) is a monotone bisection.                                    !
!==========================================================================================!
module meds_hydro_conductance
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi
   implicit none
   private

   public :: plc_retained, dplc_dpsi, flux_potential, phi_inverse, kirchhoff_edge
   public :: rhizosphere_cond

   real(wp), parameter :: dpsi_eps = 1.0e-6_wp   !< |up-down| below which K_eff -> pointwise limit
   integer(ik), parameter :: NG = 7_ik           !< Gauss-Legendre points for general-exponent Phi

   !----- 7-point Gauss-Legendre nodes/weights on [-1,1]. ----------------------------------!
   real(wp), parameter :: gl_x(NG) = [ -0.9491079123427585_wp, -0.7415311855993945_wp,        &
                                       -0.4058451513773972_wp,  0.0000000000000000_wp,        &
                                        0.4058451513773972_wp,  0.7415311855993945_wp,        &
                                        0.9491079123427585_wp ]
   real(wp), parameter :: gl_w(NG) = [  0.1294849661688697_wp,  0.2797053914892766_wp,        &
                                        0.3818300505051189_wp,  0.4179591836734694_wp,        &
                                        0.3818300505051189_wp,  0.2797053914892766_wp,        &
                                        0.1294849661688697_wp ]

contains

   !----- Fraction of conductance retained (1 - PLC). psi<0, psi50<0 => r>0; clamp for psi>0. -!
   elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      f = 1.0_wp / (1.0_wp + r**kexp)
   end function plc_retained

   !----- d(plc_retained)/d(psi); finite at psi=0 for kexp>1 (avoids the -(a/psi)f(1-f) NaN). --!
   elemental real(wp) function dplc_dpsi(psi, psi50, kexp) result(df)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r  = max(psi/psi50, 0.0_wp)
      df = -kexp * r**(kexp - 1.0_wp) / (psi50 * (1.0_wp + r**kexp)**2)
   end function dplc_dpsi

   !----- Kirchhoff (matric flux) potential Phi(psi) = integral_0^psi plc ds [MPa]. ----------!
   !      Closed form for kexp in {1,2}; fixed Gauss-Legendre on [0,r] otherwise. Phi(0)=0,     !
   !      Phi(psi<0)<0, strictly increasing in psi.                                             !
   pure real(wp) function flux_potential(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r, s, u, acc
      integer(ik) :: g
      r = max(psi/psi50, 0.0_wp)                 ! r >= 0
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = psi50 * atan(r)
      else
         !----- integral_0^r du/(1+u^kexp) by 7-pt Gauss-Legendre (map [-1,1] -> [0,r]). ----!
         acc = 0.0_wp
         do g = 1_ik, NG
            u   = 0.5_wp*r*(gl_x(g) + 1.0_wp)
            acc = acc + gl_w(g) / (1.0_wp + u**kexp)
         end do
         s   = 0.5_wp*r*acc
         phi = psi50 * s
      end if
   end function flux_potential

   !----- Inverse of Phi: find psi in [psi_lo, 0] with flux_potential(psi)=phi_target. --------!
   !      Monotone => robust bisection. Used by the vertical-profile diagnostic (3-node).       !
   pure real(wp) function phi_inverse(phi_target, psi50, kexp, psi_lo) result(psi)
      real(wp), intent(in) :: phi_target, psi50, kexp, psi_lo
      real(wp) :: lo, hi, mid, fmid
      integer(ik) :: it
      lo = psi_lo ; hi = 0.0_wp
      do it = 1_ik, 60_ik
         mid  = 0.5_wp*(lo + hi)
         fmid = flux_potential(mid, psi50, kexp)
         if (fmid < phi_target) then ; lo = mid ; else ; hi = mid ; end if
      end do
      psi = 0.5_wp*(lo + hi)
   end function phi_inverse

   !----- Kirchhoff edge conductance K_eff = k_cond * <plc> [kg/s/MPa]. k_cond is the maximum    !
   !      (plc=1) whole-plant/segment conductance already scaled to per-plant [kg/s/MPa]. The     !
   !      |dpsi|->0 limit is the pointwise value (L'Hopital of DeltaPhi/Deltapsi).                 !
   pure real(wp) function kirchhoff_edge(psi_up, psi_down, k_cond, psi50, kexp) result(keff)
      real(wp), intent(in) :: psi_up, psi_down, k_cond, psi50, kexp
      real(wp) :: dpsi
      dpsi = psi_up - psi_down
      if (abs(dpsi) > dpsi_eps) then
         keff = k_cond * ( flux_potential(psi_up,   psi50, kexp)                              &
                         - flux_potential(psi_down, psi50, kexp) ) / dpsi
      else
         keff = k_cond * plc_retained(0.5_wp*(psi_up + psi_down), psi50, kexp)
      end if
   end function kirchhoff_edge

   !----- Optional TEST helper: Katul-2003 rhizosphere (soil->root) conductance [kg/s/MPa] per   !
   !      plant, from soil conductance, fine-root biomass and specific root area. Production      !
   !      fills rhizo_cond from a future soil module; do not use on the hot path.                 !
   pure real(wp) function rhizosphere_cond(soil_cond, broot, sra, root_frac, dz, nplant) result(gw)
      real(wp), intent(in) :: soil_cond   !< [kg/m2/s] soil hydraulic conductance
      real(wp), intent(in) :: broot       !< [kgC] fine-root biomass (per plant)
      real(wp), intent(in) :: sra         !< [m2/kgC] specific root area
      real(wp), intent(in) :: root_frac   !< [-] fraction of roots in this layer
      real(wp), intent(in) :: dz          !< [m] layer thickness
      real(wp), intent(in) :: nplant      !< [pl/m2] plant density
      real(wp) :: rai
      rai = broot * sra * root_frac * nplant                     ! root area index [m2/m2]
      gw  = soil_cond * sqrt(max(rai, 0.0_wp)) / (pi * dz) / max(nplant, tiny(1.0_wp))
   end function rhizosphere_cond

end module meds_hydro_conductance
