! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_water_retention -- how much water a tissue or a soil layer holds at a given water        !
! potential, and the inverse. Two families of stateless, elemental curves:                     !
!                                                                                          !
!   * tissue (leaf, wood): the pressure-volume curve (Bartlett / Tyree-Hammel). A tissue's four  !
!     curve traits travel together in one record, water_curve_t, so every caller passes the      !
!     record and the tissue's biomass rather than the traits one by one.                        !
!   * soil: the retention curve (van Genuchten-Mualem by default, Campbell as an option) with     !
!     the conductivity and capacity that follow from it.                                        !
!                                                                                          !
! How water moves between the stores -- the xylem vulnerability curve, the Kirchhoff flux and   !
! the root profile -- is in meds_hydr_lib.                                                       !
!==========================================================================================!
module meds_water_retention
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   !----- Tissue pressure-volume family. ---------------------------------------------------!
   public :: water_curve_t
   public :: pv_psi_tlp, pv_rwc_tlp, rwc_from_psi, psi_from_rwc
   public :: water_content, capacitance, psi_from_water_content, clamp_water_to_capacity
   !----- Soil retention family: the soil-water analogue of the tissue curve -- theta(psi) and    !
   !      psi(theta), K(theta) and the capacity C(psi).                                           !
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   public :: soil_theta_from_psi, soil_psi_from_theta, soil_hydr_cond_from_theta, soil_moist_cap_from_psi

   !----- One tissue's pressure-volume curve: the four traits every tissue-water routine reads. ---!
   type :: water_curve_t
      real(wp) :: pi0           = 0.0_wp   !< [MPa] osmotic potential at full turgor (< 0)
      real(wp) :: elastic_mod   = 0.0_wp   !< [MPa] bulk elastic modulus
      real(wp) :: apoplast_frac = 0.0_wp   !< [-]   share of the saturated water outside the living cells
      real(wp) :: water_sat     = 0.0_wp   !< [kg H2O / kgC] water per unit biomass at saturation
   end type water_curve_t

   !----- Retention-curve family selector codes (the soil_params_t%retention field). -------!
   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   !----- Numerical floors (a property of the curve regularization, not of uptake). --------!
   real(wp), parameter :: K_MIN  = 1.16e-13_wp    !< [m/s] ED2 hydcond_min (no zero-conductance lock)
   real(wp), parameter :: C_MIN  = 1.0e-9_wp      !< [1/m] specific-capacity floor (vG C -> 0 at saturation)
   real(wp), parameter :: SE_MIN = 1.0e-9_wp      !< effective-saturation floor (residual/air-dry)
   real(wp), parameter :: rwc_floor = 1.0e-4_wp   !< keep R strictly positive in the flaccid tail

contains

   !=======================================================================================!
   !     Pressure-volume (Bartlett / Tyree-Hammel) constitutive curves.                     !
   !=======================================================================================!

   !----- Turgor loss point [MPa] (Bartlett eqn 1). ---------------------------------------!
   elemental real(wp) function pv_psi_tlp(pi0, elastic_mod) result(psi_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
   end function pv_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett eqn 2). --------------!
   elemental real(wp) function pv_rwc_tlp(pi0, elastic_mod) result(rwc_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
   end function pv_rwc_tlp

   !----- Water potential [MPa] from symplastic RWC. --------------------------------------!
   elemental real(wp) function psi_from_rwc(rwc, pi0, elastic_mod) result(psi)
      real(wp), intent(in) :: rwc, pi0, elastic_mod
      real(wp) :: rwc_tlp, r
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
      r       = max(rwc, rwc_floor)
      if (r >= rwc_tlp) then
         psi = elastic_mod*(r - rwc_tlp) + pi0/r      ! turgor + osmotic
      else
         psi = pi0/r                          ! turgor lost
      end if
   end function psi_from_rwc

   !----- Symplastic RWC from water potential (closed-form turgid inverse). ----------------!
   elemental real(wp) function rwc_from_psi(psi, pi0, elastic_mod) result(rwc)
      real(wp), intent(in) :: psi, pi0, elastic_mod
      real(wp) :: psi_tlp, b, disc
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
      if (psi >= psi_tlp) then                ! turgid (psi less negative than the TLP)
         !----- elastic_mod*R^2 - (psi+elastic_mod+pi0)*R + pi0 = 0; take the + root in [rwc_tlp, 1]. ------!
         b    = psi + elastic_mod + pi0
         disc = max(b*b - 4.0_wp*elastic_mod*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*elastic_mod)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, rwc_floor)
   end function rwc_from_psi

   !----- Total tissue water [kg, per plant] at potential psi. -----------------------------!
   !      W = W_sat_sym*R + W_apoplast, with W_sat = water_sat*biomass, W_sat_sym =           !
   !      (1-apoplast_frac)*W_sat and W_apoplast = apoplast_frac*W_sat (a constant reservoir). !
   elemental real(wp) function water_content(psi, curve, biomass) result(w)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, r
      w_sat = curve%water_sat*biomass
      r     = rwc_from_psi(psi, curve%pi0, curve%elastic_mod)
      w     = (1.0_wp - curve%apoplast_frac)*w_sat*r + curve%apoplast_frac*w_sat
   end function water_content

   !----- Capacitance C = dW/dpsi [kg/MPa, per plant] at potential psi. ---------------------!
   elemental real(wp) function capacitance(psi, curve, biomass) result(c)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      associate (pi0 => curve%pi0, elastic_mod => curve%elastic_mod)
         psi_tlp   = pi0*elastic_mod / (pi0 + elastic_mod)
         w_sat_sym = (1.0_wp - curve%apoplast_frac)*curve%water_sat*biomass
         r         = rwc_from_psi(psi, pi0, elastic_mod)
         if (psi >= psi_tlp) then
            cr = 1.0_wp / (elastic_mod - pi0/(r*r))      ! dR/dpsi, turgid  (-> 1/(elastic_mod+|pi0|) at R=1)
         else
            cr = -(r*r)/pi0                       ! dR/dpsi, flaccid (= |pi0|/psi^2)
         end if
      end associate
      c = w_sat_sym*cr
   end function capacitance

   !----- Exact inverse of water_content: water potential [MPa] from total tissue water [kg,      !
   !      per plant] (MEDS_ED2_RK45_DESIGN.md sec 4 -- diagnoses psi from the prognostic mass      !
   !      state). Recovers the symplastic RWC by removing the constant apoplastic reservoir, then   !
   !      composes the existing exact psi_from_rwc inverse -- every constituent already exists,     !
   !      this only inverts water_content's own forward map (w = (1-apo)*w_sat*rwc + apo*w_sat).    !
   elemental real(wp) function psi_from_water_content(w, curve, biomass) result(psi)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, rwc
      w_sat = curve%water_sat*biomass
      rwc   = (w - curve%apoplast_frac*w_sat) / max((1.0_wp - curve%apoplast_frac)*w_sat, tiny_num)
      psi   = psi_from_rwc(rwc, curve%pi0, curve%elastic_mod)
   end function psi_from_water_content

   !----- Cap tissue water mass at the saturation ceiling W_sat = water_sat*biomass (rwc=1): the     !
   !      slow/fast SEAM guard (MEDS_ED2_RK45_DESIGN.md P3). Mass, not psi, is the seam-continuous     !
   !      quantity across a biomass update -- a persisted leaf/wood_water_mass carries forward         !
   !      UNCHANGED into new biomass, so a capacity GROWTH simply reads as a lower rwc/psi next touch  !
   !      (the physically-correct signal that draws more water from the soil, not a defect to patch    !
   !      over). Only a capacity SHRINK needs a guard: a discontinuous biomass drop (e.g. the           !
   !      phenology dormant-canopy leaf snap-to-bare in meds_vegetation_dynamics%update_biomass_        !
   !      turnover) can leave more mass than the new, smaller ceiling admits -- a tissue state that     !
   !      is not reachable (rwc>1 has no valid inverse on the turgid PV branch). The caller bookkeeps  !
   !      the released excess (w - result) rather than silently retaining a supersaturated pool.        !
   elemental real(wp) function clamp_water_to_capacity(w, curve, biomass) result(w_capped)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      w_capped = min(w, curve%water_sat*biomass)
   end function clamp_water_to_capacity

   !=======================================================================================!
   !  SOIL retention curves (van Genuchten-Mualem default / Campbell-Clapp-Hornberger option). !
   !  Closed-form theta(psi)/psi(theta) inverses + K(theta) + C(psi), all elemental, branch-    !
   !  light and FPE-safe under -fpe0 / -Ktrap=fp (Se clamped, K floored, C guarded). Coupling to !
   !  the plant hydraulics kernel is through the POTENTIAL psi (curve-independent), so the curve  !
   !  is a free run-time choice. par_a/par_n carry {alpha,n} for vG or {psi_sat,b} for Campbell.  !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! theta(psi): volumetric water content from matric potential [m, <= 0].                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_theta_from_psi(retention, psi, theta_sat, theta_res, par_a, par_n)  &
                      result(theta)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: psi, theta_sat, theta_res, par_a, par_n
      real(wp)                :: theta, se, m, ah
      if (psi >= 0.0_wp) then
         theta = theta_sat
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         !----- par_a = psi_sat (< 0), par_n = b. --------------------------------------!
         se = min(1.0_wp, (psi / par_a) ** (-1.0_wp / par_n))
      else
         !----- van Genuchten: par_a = alpha [1/m], par_n = n. -------------------------!
         m  = 1.0_wp - 1.0_wp / par_n
         ah = par_a * abs(psi)
         se = (1.0_wp + ah ** par_n) ** (-m)
      end if
      se    = min(max(se, SE_MIN), 1.0_wp)
      theta = theta_res + (theta_sat - theta_res) * se
   end function soil_theta_from_psi

   !---------------------------------------------------------------------------------------!
   ! psi(theta): matric potential [m, <= 0] from water content -- closed-form inverse.      !
   !---------------------------------------------------------------------------------------!
   elemental function soil_psi_from_theta(retention, theta, theta_sat, theta_res, par_a, par_n) &
                      result(psi)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      real(wp)                :: psi, se, m
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, SE_MIN), 1.0_wp)
      if (se >= 1.0_wp) then
         psi = 0.0_wp
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         psi = par_a * se ** (-par_n)                       ! par_a = psi_sat, par_n = b
      else
         m   = 1.0_wp - 1.0_wp / par_n
         psi = -(1.0_wp / par_a) * (se ** (-1.0_wp / m) - 1.0_wp) ** (1.0_wp / par_n)
      end if
   end function soil_psi_from_theta

   !---------------------------------------------------------------------------------------!
   ! K(theta): unsaturated hydraulic conductivity [m/s], floored at K_MIN.                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_hydr_cond_from_theta(retention, theta, theta_sat, theta_res, par_a, par_n,    &
                      ksat) result(kcond)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n, ksat
      real(wp)                :: kcond, se, m, tmp
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, SE_MIN), 1.0_wp)
      if (retention == SOIL_RETENTION_CAMPBELL) then
         kcond = ksat * se ** (2.0_wp * par_n + 3.0_wp)     ! par_n = b
      else
         m     = 1.0_wp - 1.0_wp / par_n
         tmp   = 1.0_wp - (1.0_wp - se ** (1.0_wp / m)) ** m
         kcond = ksat * sqrt(se) * tmp * tmp                ! l = 0.5 Mualem pore-connectivity
      end if
      kcond = max(kcond, K_MIN)
   end function soil_hydr_cond_from_theta

   !---------------------------------------------------------------------------------------!
   ! C(psi) = dtheta/dpsi: specific moisture capacity [1/m], >= 0, floored at C_MIN.        !
   !---------------------------------------------------------------------------------------!
   elemental function soil_moist_cap_from_psi(retention, psi, theta_sat, theta_res, par_a, par_n)      &
                      result(cap)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: psi, theta_sat, theta_res, par_a, par_n
      real(wp)                :: cap, m, ah, dtr
      dtr = theta_sat - theta_res
      if (psi >= 0.0_wp) then
         cap = C_MIN
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         !----- C = -(theta_sat-theta_res)/(b*psi_sat) * (psi/psi_sat)^(-1/b - 1). ------!
         cap = -dtr / (par_n * par_a) * (psi / par_a) ** (-1.0_wp / par_n - 1.0_wp)
      else
         m   = 1.0_wp - 1.0_wp / par_n
         ah  = par_a * abs(psi)
         cap = par_a * m * par_n * dtr * ah ** (par_n - 1.0_wp)                            &
               * (1.0_wp + ah ** par_n) ** (-m - 1.0_wp)
      end if
      cap = max(cap, C_MIN)
   end function soil_moist_cap_from_psi

end module meds_water_retention
