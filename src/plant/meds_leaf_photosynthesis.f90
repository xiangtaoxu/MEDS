!==========================================================================================!
! meds_leaf_photosynthesis -- the CO2 demand functions: gross assimilation as a function of !
! intercellular CO2, for C3 (Farquhar-von Caemmerer-Berry) and C4 (Collatz et al. 1992).   !
!                                                                                          !
! All concentrations (Ci, Gamma*, Kc, Ko, O2) are MOLE FRACTIONS [umol/mol]; capacities      !
! (Vcmax, J, TPU) and the returned rates are [umol CO2/m2/s]. The functions are pure and       !
! stateless -- they take already-temperature-scaled inputs (the temperature response and       !
! unit conversion happen in the solver), so the root-finder can call them cheaply.             !
!                                                                                          !
!   C3:  Ac = Vcmax (Ci-G*)/(Ci + Kc(1+O2/Ko))   [Rubisco]                                  !
!        Aj = J (Ci-G*)/(4 Ci + 8 G*)            [RuBP / light, via electron_transport_j]    !
!        Ap = 3 TPU                              [triose-phosphate use]                       !
!   C4:  Ac = Vcmax,  Aj = light-limited slope,  Ap = kp Ci (PEPcase),  G* ~ 0               !
!                                                                                          !
! Gross assimilation combines the limitation rates either by a sharp minimum (COLIM_MIN) or   !
! by smoothing quadratics (COLIM_QUADRATIC, C1-continuous -> a smooth solver residual).        !
!==========================================================================================!
module meds_leaf_photosynthesis
   use meds_kinds,     only : wp, ik
   use meds_config,    only : COLIM_MIN, COLIM_QUADRATIC
   implicit none
   private

   public :: assim_demand_c3, assim_demand_c4, electron_transport_j

contains

   !---------------------------------------------------------------------------------------!
   ! Actual electron transport rate J from absorbed PAR via the non-rectangular hyperbola   !
   ! theta J^2 - (I2 + Jmax) J + I2 Jmax = 0 (smaller root). I2 = 0.5 phi_psii absorptance PAR.!
   !---------------------------------------------------------------------------------------!
   elemental pure function electron_transport_j(par, absorptance, phi_psii, jmax, theta) result(j)
      real(wp), intent(in) :: par         !< [umol photon/m2/s] incident PAR
      real(wp), intent(in) :: absorptance !< [--] leaf PAR absorptance
      real(wp), intent(in) :: phi_psii    !< [--] PSII quantum yield (electrons/photon)
      real(wp), intent(in) :: jmax        !< [umol/m2/s] electron-transport capacity (T-scaled)
      real(wp), intent(in) :: theta       !< [--] curvature (0 < theta < 1)
      real(wp)             :: j, i2
      i2 = 0.5_wp * phi_psii * absorptance * par
      j  = smaller_root(theta, i2, jmax)
   end function electron_transport_j

   !---------------------------------------------------------------------------------------!
   ! C3 demand: gross assimilation a_gross and the three raw limitation rates (ac/aj/ap).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine assim_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colim, theta,        &
                                   a_gross, ac, aj, ap)
      real(wp),    intent(in)  :: ci, vcmax, j, tpu, gstar, kc, ko, o2, theta
      integer(ik), intent(in)  :: colim
      real(wp),    intent(out) :: a_gross, ac, aj, ap
      ac = vcmax * (ci - gstar) / (ci + kc * (1.0_wp + o2 / ko))
      aj = j     * (ci - gstar) / (4.0_wp * ci + 8.0_wp * gstar)
      ap = 3.0_wp * tpu
      a_gross = combine_limits(ac, aj, ap, colim, theta)
   end subroutine assim_demand_c3

   !---------------------------------------------------------------------------------------!
   ! C4 demand (Collatz 1992): ac = Vcmax, aj = light-limited slope (passed in), ap = PEPcase !
   ! CO2 limitation kp_eff*Ci. Gamma* ~ 0 (CO2-concentrating mechanism suppresses photoresp). !
   !---------------------------------------------------------------------------------------!
   pure subroutine assim_demand_c4(ci, vcmax, aj_light, kp_eff, colim, theta_cj, theta_ic,    &
                                   a_gross, ac, aj, ap)
      real(wp),    intent(in)  :: ci, vcmax, aj_light, kp_eff, theta_cj, theta_ic
      integer(ik), intent(in)  :: colim
      real(wp),    intent(out) :: a_gross, ac, aj, ap
      real(wp) :: ai
      ac = vcmax
      aj = aj_light
      ap = kp_eff * ci
      if (colim == COLIM_QUADRATIC) then
         ai      = smaller_root(theta_cj, ac, aj)        ! co-limit Rubisco & light
         a_gross = smaller_root(theta_ic, ai, ap)        ! co-limit with PEPcase CO2
      else
         a_gross = min(ac, aj, ap)
      end if
   end subroutine assim_demand_c4

   !---------------------------------------------------------------------------------------!
   ! Combine the three C3 limitation rates: sharp min, or two nested smoothing quadratics.  !
   !---------------------------------------------------------------------------------------!
   pure function combine_limits(ac, aj, ap, colim, theta) result(a)
      real(wp),    intent(in) :: ac, aj, ap, theta
      integer(ik), intent(in) :: colim
      real(wp)                :: a, ai
      if (colim == COLIM_QUADRATIC) then
         ai = smaller_root(theta, ac, aj)                ! co-limit Rubisco & RuBP
         a  = smaller_root(theta, ai, ap)                ! co-limit with product (TPU)
      else
         a = min(ac, aj, ap)
      end if
   end function combine_limits

   !---------------------------------------------------------------------------------------!
   ! Smaller root of the co-limitation quadratic theta x^2 - (a+b) x + a b = 0. The smaller  !
   ! root is the smooth analogue of min(a,b): it approaches min(a,b) as theta -> 1 and the     !
   ! harmonic-like co-limited mean for smaller theta. Discriminant is guarded >= 0.            !
   !---------------------------------------------------------------------------------------!
   elemental pure function smaller_root(theta, a, b) result(x)
      real(wp), intent(in) :: theta, a, b
      real(wp)             :: x, s, disc
      s    = a + b
      disc = max(s * s - 4.0_wp * theta * a * b, 0.0_wp)
      x    = (s - sqrt(disc)) / (2.0_wp * theta)
   end function smaller_root

end module meds_leaf_photosynthesis
