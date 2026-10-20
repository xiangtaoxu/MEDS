! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_water_retention -- how much water a store holds at a given water potential, and the !
! inverse: a plant tissue's pressure-volume curve (docs/science/plant_hydraulics.md §2), a !
! soil layer's retention curve and conductivity (docs/science/soil_biophysics.md eq. (2))  !
! and the pond's depression store.                                                         !
!                                                                                          !
! Each comes in two forms. The read-offs (rwc_from_psi, soil_psi_from_theta, ...) take a   !
! state and return a potential or a water content, with floors that keep a diagnostic      !
! finite. The water relations (*_water_relations) give the soil-plant water solver         !
! (meds_spac_solver) a store's potential, water and conductivity with their derivatives in !
! the variable it solves for, with no floors (docs/science/soil_plant_water.md §2). How    !
! water moves between the stores is in meds_hydr_lib, and between soil layers in           !
! meds_soil_water.                                                                         !
!==========================================================================================!
module meds_water_retention
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav_head, rho_h2o, tiny_num
   implicit none
   private

   !----- Tissue pressure-volume curve. ---------------------------------------------------!
   public :: water_curve_t
   public :: pv_psi_tlp, rwc_from_psi, psi_from_rwc
   public :: water_content, capacitance, psi_from_water_content, clamp_water_to_capacity, shared_potential
   !----- Soil retention curve: theta(psi) and psi(theta) read-offs. ----------------------!
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, VG_PORE_CONNECTIVITY
   public :: soil_theta_from_psi, soil_psi_from_theta
   !----- The water relations the soil-plant water solver uses, in the variable it solves  !
   !      for: a soil layer's saturation (dry) or potential (wet), a tissue's relative     !
   !      water content, the pond's head. -------------------------------------------------!
   public :: SOIL_SPECIFIC_STORAGE, POND_DEPRESSION_DEPTH
   public :: soil_layer_curve_t, build_soil_layer_curve, soil_water_relations, soil_unknown_from_theta
   public :: tissue_water_relations, tissue_unknown_from_rwc, tissue_rwc_from_water
   public :: pond_water_relations, pond_head_from_water

   !----- One tissue's pressure-volume curve: the four traits every tissue-water routine   !
   !      reads, passed together with the tissue's biomass. -------------------------------!
   type :: water_curve_t
      real(wp) :: pi0           = 0.0_wp   !< [MPa] osmotic potential at full turgor (< 0)
      real(wp) :: elastic_mod   = 0.0_wp   !< [MPa] bulk elastic modulus
      real(wp) :: apoplast_frac = 0.0_wp   !< [-]   share of the saturated water outside the living cells
      real(wp) :: water_sat     = 0.0_wp   !< [kg H2O / kgC] water per unit biomass at saturation
   end type water_curve_t

   !----- Retention-curve selector codes (soil_params_t%retention). -----------------------!
   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   !----- Mualem's pore-connectivity exponent L in the van Genuchten conductivity          !
   !      (docs/science/soil_biophysics.md eq. (2)), set here only: -1 after Schaap & Leij !
   !      (2000, SSSA J. 64:843), which keeps a soil of n near 1.4 conducting down to      !
   !      about -2 MPa. -------------------------------------------------------------------!
   real(wp), parameter :: VG_PORE_CONNECTIVITY = -1.0_wp

   !----- Read-off floors that keep a diagnostic finite; the solver's curves have none. ---!
   real(wp), parameter :: SE_MIN = 1.0e-9_wp      !< [-] effective-saturation floor
   real(wp), parameter :: rwc_floor = 1.0e-4_wp   !< [-] symplastic relative water content floor
   !----- Cap on shared_potential's halvings. The loop stops when its bracket is as narrow !
   !      as the arithmetic resolves, about 70 halvings for a root near -1 MPa; only a     !
   !      root at 0, where that width has no floor, reaches the cap. ----------------------!
   integer(ik), parameter :: BISECTION_STEPS_MAX = 200_ik

   !----- [1/m] Water a saturated soil layer gains per metre of pressure head, from the    !
   !      compressibility of the soil matrix and of water: a typical value for clays and   !
   !      loams (Freeze and Cherry 1979, table 2.5). It acts only above saturation         !
   !      (docs/science/soil_plant_water.md §2). ------------------------------------------!
   real(wp), parameter :: SOIL_SPECIFIC_STORAGE = 1.0e-4_wp
   !----- [m] Depth of the surface's small hollows, which the pond fills before a head     !
   !      builds, so its storage stays smooth as it fills or empties                       !
   !      (docs/science/soil_plant_water.md §2). ------------------------------------------!
   real(wp), parameter :: POND_DEPRESSION_DEPTH = 1.0e-3_wp
   !----- The largest logarithm the dry-soil head is evaluated at: e^700 m of head is far  !
   !      beyond any physical soil and still finite in double precision. ------------------!
   real(wp), parameter :: LOG_HEAD_MAX = 700.0_wp

   !----- One soil layer's curve, set up once per layer: the retention parameters and the  !
   !      point where the solver's unknown changes from saturation to potential. ----------!
   type :: soil_layer_curve_t
      integer(ik) :: retention = 1_ik   !< SOIL_RETENTION_VG or SOIL_RETENTION_CAMPBELL
      real(wp) :: theta_sat = 0.0_wp    !< [m3/m3]
      real(wp) :: theta_res = 0.0_wp    !< [m3/m3]
      real(wp) :: par_a     = 0.0_wp    !< vG alpha [1/m], or Campbell's air-entry potential [m, < 0]
      real(wp) :: par_n     = 0.0_wp    !< vG n [-], or Campbell's b [-]
      real(wp) :: ksat      = 0.0_wp    !< [m/s] saturated conductivity
      real(wp) :: dz        = 0.0_wp    !< [m] layer thickness
      real(wp) :: se_switch  = 1.0_wp   !< [-] saturation where the unknown changes branch
      real(wp) :: psi_switch = 0.0_wp   !< [m] the potential there
      real(wp) :: slope_switch = 1.0_wp !< [1/m] dSe/dpsi there, the unknown's scale on the wet side
   end type soil_layer_curve_t

contains

   !=======================================================================================!
   !     Pressure-volume (Bartlett / Tyree-Hammel) curves                                  !
   !     (docs/science/plant_hydraulics.md §2, eqs. (3)-(6)).                              !
   !=======================================================================================!

   !----- Turgor loss point [MPa] (Bartlett et al. 2012 eq. 1;                             !
   !      docs/science/plant_hydraulics.md eq. (3)). --------------------------------------!
   elemental real(wp) function pv_psi_tlp(pi0, elastic_mod) result(psi_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
   end function pv_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett et al. 2012 eq. 2;    !
   !      docs/science/plant_hydraulics.md eq. (3)). --------------------------------------!
   elemental real(wp) function pv_rwc_tlp(pi0, elastic_mod) result(rwc_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
   end function pv_rwc_tlp

   !----- Water potential [MPa] at symplastic relative water content rwc, floored at       !
   !      rwc_floor (docs/science/plant_hydraulics.md eq. (4)). ---------------------------!
   elemental real(wp) function psi_from_rwc(rwc, pi0, elastic_mod) result(psi)
      real(wp), intent(in) :: rwc, pi0, elastic_mod
      real(wp) :: rwc_tlp, r
      rwc_tlp = pv_rwc_tlp(pi0, elastic_mod)
      r       = max(rwc, rwc_floor)
      if (r >= rwc_tlp) then
         psi = elastic_mod*(r - rwc_tlp) + pi0/r      ! turgor + osmotic
      else
         psi = pi0/r                          ! turgor lost
      end if
   end function psi_from_rwc

   !----- Symplastic relative water content at water potential psi, floored at rwc_floor   !
   !      (docs/science/plant_hydraulics.md eq. (4)). -------------------------------------!
   elemental real(wp) function rwc_from_psi(psi, pi0, elastic_mod) result(rwc)
      real(wp), intent(in) :: psi, pi0, elastic_mod
      real(wp) :: psi_tlp, b, disc
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
      if (psi >= psi_tlp) then                ! turgid
         !----- The larger root of elastic_mod R^2 - b R + pi0 = 0. -----------------------!
         b    = psi + elastic_mod + pi0
         disc = max(b*b - 4.0_wp*elastic_mod*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*elastic_mod)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, rwc_floor)
   end function rwc_from_psi

   !----- Total tissue water [kg per plant] at potential psi: the symplast's saturated     !
   !      water times R, plus the apoplast's, a constant reservoir                         !
   !      (docs/science/plant_hydraulics.md eq. (5)). -------------------------------------!
   elemental real(wp) function water_content(psi, curve, biomass) result(w)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, r
      w_sat = curve%water_sat*biomass
      r     = rwc_from_psi(psi, curve%pi0, curve%elastic_mod)
      w     = (1.0_wp - curve%apoplast_frac)*w_sat*r + curve%apoplast_frac*w_sat
   end function water_content

   !----- Capacitance dW/dpsi [kg/MPa per plant] at potential psi                          !
   !      (docs/science/plant_hydraulics.md eqs. (5)-(6)). --------------------------------!
   elemental real(wp) function capacitance(psi, curve, biomass) result(c)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      associate (pi0 => curve%pi0, elastic_mod => curve%elastic_mod)
         psi_tlp   = pi0*elastic_mod / (pi0 + elastic_mod)
         w_sat_sym = (1.0_wp - curve%apoplast_frac)*curve%water_sat*biomass
         r         = rwc_from_psi(psi, pi0, elastic_mod)
         if (psi >= psi_tlp) then
            cr = 1.0_wp / (elastic_mod - pi0/(r*r))      ! dR/dpsi, turgid
         else
            cr = -(r*r)/pi0                       ! dR/dpsi, flaccid
         end if
      end associate
      c = w_sat_sym*cr
   end function capacitance

   !----- Water potential [MPa] of a tissue holding w [kg per plant], the exact inverse of !
   !      water_content (docs/science/plant_hydraulics.md eq. (5)). -----------------------!
   elemental real(wp) function psi_from_water_content(w, curve, biomass) result(psi)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, rwc
      w_sat = curve%water_sat*biomass
      rwc   = (w - curve%apoplast_frac*w_sat) / max((1.0_wp - curve%apoplast_frac)*w_sat, tiny_num)
      psi   = psi_from_rwc(rwc, curve%pi0, curve%elastic_mod)
   end function psi_from_water_content

   !----- The potential [MPa] at which two tissues holding w_total [kg] between them are   !
   !      at rest, a at psi and b at psi - head_b (docs/science/plant_hydraulics.md        !
   !      eq. (2) with no flow), by bisection between the curves' dry floor and a's        !
   !      saturation; a total outside returns the nearer end. The caller moves water by    !
   !      the curves at the result, so what it moves is exact whatever the residual. ------!
   elemental real(wp) function shared_potential(w_total, curve_a, biomass_a, curve_b, biomass_b, head_b)  &
                                                result(psi)
      real(wp),            intent(in) :: w_total               !< [kg] the water the two hold together
      type(water_curve_t), intent(in) :: curve_a, curve_b
      real(wp),            intent(in) :: biomass_a, biomass_b  !< [kgC]
      real(wp),            intent(in) :: head_b                !< [MPa] gravity head of b above a: at rest psi_b = psi_a - head_b
      real(wp)    :: lo, hi
      integer(ik) :: k
      lo = min(curve_a%pi0, curve_b%pi0 + head_b) / rwc_floor
      hi = 0.0_wp
      if (w_total <= held(lo)) then
         psi = lo
         return
      end if
      if (w_total >= held(hi)) then
         psi = hi
         return
      end if
      do k = 1_ik, BISECTION_STEPS_MAX
         psi = 0.5_wp * (lo + hi)
         if (held(psi) < w_total) then
            lo = psi
         else
            hi = psi
         end if
         if (hi - lo <= epsilon(psi) * max(abs(lo), abs(hi))) exit
      end do
      psi = 0.5_wp * (lo + hi)
   contains
      pure real(wp) function held(p)
         real(wp), intent(in) :: p
         held = water_content(p, curve_a, biomass_a) + water_content(p - head_b, curve_b, biomass_b)
      end function held
   end function shared_potential

   !----- Tissue water [kg per plant] capped at saturation, water_sat * biomass, for a     !
   !      store that lost biomass without its water; the caller sheds the excess           !
   !      (w - result) to the ground (docs/science/soil_plant_water.md §2). ---------------!
   elemental real(wp) function clamp_water_to_capacity(w, curve, biomass) result(w_capped)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      w_capped = min(w, curve%water_sat*biomass)
   end function clamp_water_to_capacity

   !=======================================================================================!
   !     Soil retention curves: van Genuchten-Mualem (default) or Campbell                 !
   !     (docs/science/soil_biophysics.md eq. (2)). par_a and par_n are alpha [1/m] and n  !
   !     for van Genuchten, psi_sat [m, < 0] and b for Campbell; heads are in metres of    !
   !     water, below zero when unsaturated. Every soil-water routine evaluates the curve  !
   !     through the four kernels below, so a diagnostic potential is the solver's own.    !
   !=======================================================================================!

   !----- Volumetric water content [m3/m3] at matric head psi [m]: theta_sat at or above   !
   !      zero, and the saturation floored at SE_MIN. -------------------------------------!
   elemental function soil_theta_from_psi(retention, psi, theta_sat, theta_res, par_a, par_n)  &
                      result(theta)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: psi, theta_sat, theta_res, par_a, par_n
      real(wp)                :: theta, se, dse_dhead
      if (psi >= 0.0_wp) then
         theta = theta_sat
         return
      end if
      call saturation_of_head(retention, par_a, par_n, psi, se, dse_dhead)
      se    = min(max(se, SE_MIN), 1.0_wp)
      theta = theta_res + (theta_sat - theta_res) * se
   end function soil_theta_from_psi

   !----- Matric head [m, <= 0] of a layer holding theta [m3/m3], its saturation floored   !
   !      at SE_MIN; zero at saturation. --------------------------------------------------!
   elemental function soil_psi_from_theta(retention, theta, theta_sat, theta_res, par_a, par_n) &
                      result(psi)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      real(wp)                :: psi, se, dse_dhead
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, SE_MIN), 1.0_wp)
      if (se >= 1.0_wp) then
         psi = 0.0_wp
         return
      end if
      call head_of_saturation(retention, par_a, par_n, se, psi, dse_dhead)
   end function soil_psi_from_theta

   !=======================================================================================!
   !     Water relations for the soil-plant water solver: soil layers.                     !
   !=======================================================================================!

   !----- Set up one layer's curve and its branch point, where the curve holds the most    !
   !      water per unit potential: (alpha |psi|)^n = m, Se = (1+m)^-m, for van Genuchten; !
   !      the air-entry point for Campbell (docs/science/soil_plant_water.md §2). ---------!
   pure subroutine build_soil_layer_curve(retention, theta_sat, theta_res, par_a, par_n, ksat, dz, c)
      integer(ik),              intent(in)  :: retention
      real(wp),                 intent(in)  :: theta_sat, theta_res, par_a, par_n, ksat, dz
      type(soil_layer_curve_t), intent(out) :: c
      real(wp) :: m, x
      c%retention = retention ; c%theta_sat = theta_sat ; c%theta_res = theta_res
      c%par_a = par_a ; c%par_n = par_n ; c%ksat = ksat ; c%dz = dz
      if (retention == SOIL_RETENTION_CAMPBELL) then
         c%se_switch    = 1.0_wp
         c%psi_switch   = par_a
         c%slope_switch = -1.0_wp / (par_n * par_a)
      else
         m = 1.0_wp - 1.0_wp / par_n
         x = m ** (1.0_wp / par_n)
         c%se_switch    = (1.0_wp + m) ** (-m)
         c%psi_switch   = -x / par_a
         c%slope_switch = par_a * m * par_n * c%se_switch * (m / (1.0_wp + m)) / x
      end if
   end subroutine build_soil_layer_curve

   !----- A soil layer's potential [MPa], water [kg/m2] and conductivity [m/s], with their !
   !      derivatives in the solver's unknown tau (docs/science/soil_plant_water.md §2).   !
   !      Dry of the branch point tau is the effective saturation Se; wet of it tau runs   !
   !      linearly with the potential, at the slope that joins the two branches smoothly.  !
   !      Above saturation the layer stores water through SOIL_SPECIFIC_STORAGE instead of !
   !      being clipped. ------------------------------------------------------------------!
   elemental subroutine soil_water_relations(c, tau, psi, dpsi, w, dw, kcond, dkcond)
      type(soil_layer_curve_t), intent(in)  :: c
      real(wp),                 intent(in)  :: tau
      real(wp),                 intent(out) :: psi, dpsi    !< [MPa], [MPa per unit tau]
      real(wp),                 intent(out) :: w, dw        !< [kg/m2], [kg/m2 per unit tau]
      real(wp),                 intent(out) :: kcond, dkcond !< [m/s], [m/s per unit tau]
      real(wp) :: head, dhead, se, dse_dhead, dk_dse, dk_dhead, dtheta
      dtheta = c%theta_sat - c%theta_res
      if (tau <= c%se_switch) then
         !----- Dry branch: tau is the saturation itself. ---------------------------------!
         se = tau
         call head_of_saturation(c%retention, c%par_a, c%par_n, se, head, dse_dhead)
         call conductivity_of_saturation(c%retention, c%par_n, c%ksat, se, kcond, dk_dse)
         psi    = grav_head * head
         dpsi   = grav_head / max(dse_dhead, tiny_num)   ! clamp-ok: a slope underflowed in extreme dry soil
         w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
         dw     = rho_h2o * c%dz * dtheta
         dkcond = dk_dse
      else
         !----- Wet branch: tau measures the potential. -----------------------------------!
         head  = c%psi_switch + (tau - c%se_switch) / c%slope_switch
         dhead = 1.0_wp / c%slope_switch
         psi   = grav_head * head
         dpsi  = grav_head * dhead
         if (head < 0.0_wp) then
            call saturation_of_head(c%retention, c%par_a, c%par_n, head, se, dse_dhead)
            call conductivity_of_head(c%retention, c%par_a, c%par_n, c%ksat, head, kcond, dk_dhead)
            w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
            dw     = rho_h2o * c%dz * dtheta * dse_dhead * dhead
            dkcond = dk_dhead * dhead
         else
            w      = rho_h2o * c%dz * (c%theta_sat + SOIL_SPECIFIC_STORAGE * head)
            dw     = rho_h2o * c%dz * SOIL_SPECIFIC_STORAGE * dhead
            kcond  = c%ksat
            dkcond = 0.0_wp
         end if
      end if
   end subroutine soil_water_relations

   !----- The unknown tau of a layer holding theta [m3/m3], the inverse of                 !
   !      soil_water_relations's water. A layer at or below residual starts just inside    !
   !      the curve, at SE_MIN; a Campbell layer at saturation sits at the branch point, !
   !      its air-entry potential. --------------------------------------------------------!
   elemental function soil_unknown_from_theta(c, theta) result(tau)
      type(soil_layer_curve_t), intent(in) :: c
      real(wp),                 intent(in) :: theta
      real(wp) :: tau
      real(wp) :: se, head, dse_dhead
      se = (theta - c%theta_res) / max(c%theta_sat - c%theta_res, tiny_num)
      if (se <= c%se_switch) then
         tau = max(se, SE_MIN)     ! clamp-ok: a state at residual starts just inside the curve
         return
      end if
      if (theta > c%theta_sat) then
         head = (theta - c%theta_sat) / SOIL_SPECIFIC_STORAGE
      else if (se >= 1.0_wp) then
         head = 0.0_wp
      else
         call head_of_saturation(c%retention, c%par_a, c%par_n, se, head, dse_dhead)
      end if
      tau = c%se_switch + c%slope_switch * (head - c%psi_switch)
   end function soil_unknown_from_theta

   !----- The curve's kernels. ------------------------------------------------------------!

   !----- Matric head [m] at saturation se in (0, 1), and dSe/dhead [1/m]. The van         !
   !      Genuchten head is evaluated through its logarithm, so a very dry layer gives a   !
   !      large but finite number, and its slope Se (1 - Se^(1/m)) / x is written without  !
   !      a cancellation. -----------------------------------------------------------------!
   elemental subroutine head_of_saturation(retention, par_a, par_n, se, head, dse_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, se
      real(wp),    intent(out) :: head, dse_dhead
      real(wp) :: m, u, logx, x
      if (retention == SOIL_RETENTION_CAMPBELL) then
         head      = par_a * se ** (-par_n)
         dse_dhead = -se ** (1.0_wp + par_n) / (par_n * par_a)
      else
         m    = 1.0_wp - 1.0_wp / par_n
         u    = 1.0_wp - se ** (1.0_wp / m)
         logx = (-log(se) / m + log(u)) / par_n
         x    = exp(min(logx, LOG_HEAD_MAX))   ! clamp-ok: keeps an extreme dry head finite
         head = -x / par_a
         dse_dhead = par_a * m * par_n * se * u / x
      end if
   end subroutine head_of_saturation

   !----- Saturation and dSe/dhead [1/m] at a head below zero, written in x = alpha |head| !
   !      so that 1 - Se^(1/m) = x^n / (1 + x^n) stays exact near saturation. A Campbell   !
   !      layer between its air-entry point and zero is saturated. ------------------------!
   elemental subroutine saturation_of_head(retention, par_a, par_n, head, se, dse_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, head
      real(wp),    intent(out) :: se, dse_dhead
      real(wp) :: m, x, xn, u
      if (retention == SOIL_RETENTION_CAMPBELL) then
         if (head < par_a) then
            se        = (head / par_a) ** (-1.0_wp / par_n)
            dse_dhead = -se ** (1.0_wp + par_n) / (par_n * par_a)
         else
            se = 1.0_wp ; dse_dhead = 0.0_wp
         end if
         return
      end if
      m  = 1.0_wp - 1.0_wp / par_n
      x  = par_a * abs(head)
      xn = x ** par_n
      se = (1.0_wp + xn) ** (-m)
      u  = xn / (1.0_wp + xn)                      ! 1 - Se^(1/m)
      dse_dhead = par_a * m * par_n * se * u / max(x, tiny_num)
   end subroutine saturation_of_head

   !----- Conductivity [m/s] and dK/dSe at saturation se: van Genuchten-Mualem with        !
   !      L = VG_PORE_CONNECTIVITY (docs/science/soil_biophysics.md eq. (2)), or           !
   !      Campbell's Ks Se^(2b+3). No floor: K goes to zero with Se. The van Genuchten     !
   !      slope is unbounded at Se = 1; the solver calls this only dry of the branch       !
   !      point. --------------------------------------------------------------------------!
   elemental subroutine conductivity_of_saturation(retention, par_n, ksat, se, kcond, dk_dse)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_n, ksat, se
      real(wp),    intent(out) :: kcond, dk_dse
      real(wp) :: m, u, t
      if (retention == SOIL_RETENTION_CAMPBELL) then
         kcond  = ksat * se ** (2.0_wp * par_n + 3.0_wp)
         dk_dse = (2.0_wp * par_n + 3.0_wp) * kcond / se
      else
         m  = 1.0_wp - 1.0_wp / par_n
         u  = 1.0_wp - se ** (1.0_wp / m)
         t  = 1.0_wp - u ** m
         kcond  = ksat * se ** VG_PORE_CONNECTIVITY * t * t
         dk_dse = ksat * (VG_PORE_CONNECTIVITY * se ** (VG_PORE_CONNECTIVITY - 1.0_wp) * t * t                  &
                  + 2.0_wp * se ** VG_PORE_CONNECTIVITY * t * u ** (m - 1.0_wp) * se ** (1.0_wp / m - 1.0_wp))
      end if
   end subroutine conductivity_of_saturation

   !----- Conductivity [m/s] and dK/dhead [1/s] at a head below zero, in x = alpha |head|  !
   !      like saturation_of_head. The slope is taken in the head directly: near           !
   !      saturation dK/dSe grows without bound while dSe/dhead goes to zero. -------------!
   elemental subroutine conductivity_of_head(retention, par_a, par_n, ksat, head, kcond, dk_dhead)
      integer(ik), intent(in)  :: retention
      real(wp),    intent(in)  :: par_a, par_n, ksat, head
      real(wp),    intent(out) :: kcond, dk_dhead
      real(wp) :: m, x, xn, u, t, se, dse_dhead, dk_dse
      if (retention == SOIL_RETENTION_CAMPBELL) then
         call saturation_of_head(retention, par_a, par_n, head, se, dse_dhead)
         if (head < par_a) then
            call conductivity_of_saturation(retention, par_n, ksat, se, kcond, dk_dse)
            dk_dhead = dk_dse * dse_dhead
         else
            kcond = ksat ; dk_dhead = 0.0_wp
         end if
         return
      end if
      m  = 1.0_wp - 1.0_wp / par_n
      x  = par_a * abs(head)
      xn = x ** par_n
      se = (1.0_wp + xn) ** (-m)
      u  = xn / (1.0_wp + xn)                      ! 1 - Se^(1/m)
      dse_dhead = par_a * m * par_n * se * u / max(x, tiny_num)
      t     = 1.0_wp - u ** m
      kcond = ksat * se ** VG_PORE_CONNECTIVITY * t * t
      !----- dK/dhead = dK/dSe dSe/dhead, the Mualem term rewritten in u and x. -----------!
      dk_dhead = ksat * (VG_PORE_CONNECTIVITY * t * t * se ** (VG_PORE_CONNECTIVITY - 1.0_wp) * dse_dhead          &
                 + 2.0_wp * par_a * m * par_n * t * se ** (VG_PORE_CONNECTIVITY + 1.0_wp / m) * u ** m             &
                   / max(x, tiny_num))
   end subroutine conductivity_of_head

   !=======================================================================================!
   !     Plant tissues.                                                                    !
   !=======================================================================================!

   !----- A tissue at symplastic relative water content r: potential [MPa] and water       !
   !      [kg/m2 of ground] for nplant plants of `biomass` each, with their derivatives in !
   !      r. The pressure-volume curve without psi_from_rwc's floor: the solver keeps r    !
   !      above zero (docs/science/plant_hydraulics.md eqs. (4)-(5)). ---------------------!
   elemental subroutine tissue_water_relations(curve, biomass, nplant, r, psi, dpsi, w, dw)
      type(water_curve_t), intent(in)  :: curve
      real(wp),            intent(in)  :: biomass   !< [kgC per plant]
      real(wp),            intent(in)  :: nplant    !< [plants/m2]
      real(wp),            intent(in)  :: r         !< [-] symplastic relative water content
      real(wp),            intent(out) :: psi, dpsi, w, dw
      real(wp) :: r_tlp, w_sat
      associate (pi0 => curve%pi0, emod => curve%elastic_mod, apo => curve%apoplast_frac)
         r_tlp = pv_rwc_tlp(pi0, emod)
         if (r >= r_tlp) then
            psi  = emod * (r - r_tlp) + pi0 / r
            dpsi = emod - pi0 / (r * r)
         else
            psi  = pi0 / r
            dpsi = -pi0 / (r * r)
         end if
         w_sat = nplant * curve%water_sat * biomass
         w     = w_sat * ((1.0_wp - apo) * r + apo)
         dw    = w_sat * (1.0_wp - apo)
      end associate
   end subroutine tissue_water_relations

   !----- Symplastic relative water content of a tissue holding w [kg per plant], the      !
   !      inverse of tissue_water_relations's water. Undefined with no biomass: the caller !
   !      then starts the tissue from a potential instead (rwc_from_psi). -----------------!
   elemental function tissue_rwc_from_water(curve, biomass, w) result(r)
      type(water_curve_t), intent(in) :: curve
      real(wp),            intent(in) :: biomass, w
      real(wp) :: r, w_sat
      w_sat = curve%water_sat * biomass
      r = (w / w_sat - curve%apoplast_frac) / (1.0_wp - curve%apoplast_frac)
   end function tissue_rwc_from_water

   !----- A tissue's unknown is its relative water content, kept strictly positive: zero   !
   !      is an infinitely dry tissue. ----------------------------------------------------!
   elemental function tissue_unknown_from_rwc(r) result(tau)
      real(wp), intent(in) :: r
      real(wp) :: tau
      real(wp), parameter :: R_START = 1.0e-6_wp
      tau = max(r, R_START)   ! clamp-ok: an emptied store starts just inside the curve
   end function tissue_unknown_from_rwc

   !=======================================================================================!
   !     The pond.                                                                         !
   !=======================================================================================!

   !----- The pond at head h [m]: potential [MPa], water [kg/m2] and their derivatives     !
   !      (docs/science/soil_plant_water.md §2). At or below zero head it holds nothing    !
   !      and its potential follows the soil it drains into. ------------------------------!
   elemental subroutine pond_water_relations(h, psi, dpsi, w, dw)
      real(wp), intent(in)  :: h
      real(wp), intent(out) :: psi, dpsi, w, dw
      psi  = grav_head * h
      dpsi = grav_head
      if (h <= 0.0_wp) then
         w = 0.0_wp ; dw = 0.0_wp
      else if (h < POND_DEPRESSION_DEPTH) then
         w  = rho_h2o * h * h / (2.0_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o * h / POND_DEPRESSION_DEPTH
      else
         w  = rho_h2o * (h - 0.5_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o
      end if
   end subroutine pond_water_relations

   !----- The head [m] of a pond holding w [kg/m2]. An empty pond has no head to name; it  !
   !      returns zero, and the solver finds the head from the soil below. ----------------!
   elemental function pond_head_from_water(w) result(h)
      real(wp), intent(in) :: w
      real(wp) :: h, depth
      depth = w / rho_h2o
      if (depth <= 0.0_wp) then
         h = 0.0_wp
      else if (depth < 0.5_wp * POND_DEPRESSION_DEPTH) then
         h = sqrt(2.0_wp * POND_DEPRESSION_DEPTH * depth)
      else
         h = depth + 0.5_wp * POND_DEPRESSION_DEPTH
      end if
   end function pond_head_from_water

end module meds_water_retention
