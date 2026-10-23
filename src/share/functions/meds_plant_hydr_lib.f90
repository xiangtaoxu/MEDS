! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_hydr_lib -- the water relations of a plant: how much water a tissue holds at  !
! a given potential, its pressure-volume curve (docs/science/plant_hydraulics.md §2), and  !
! how water moves into and through the plant: the vulnerability curve and its Kirchhoff    !
! potential, the xylem's conductance and sapflow, and root water uptake                    !
! (plant_hydraulics.md §3-§4, docs/science/soil_plant_water.md §1).                        !
!                                                                                          !
! The read-offs (tissue_psi_to_rwc, tissue_water_to_psi, ...) carry floors that keep a     !
! diagnostic finite. The solver's relations (tissue_water_relations and the flows) have    !
! none and come with their derivatives in the solver's unknowns, so the water network     !
! builds its Jacobian from the numbers it moves water with (soil_plant_water.md §2). The   !
! soil's side is meds_soil_hydr_lib.                                                       !
!==========================================================================================!
module meds_plant_hydr_lib
   use meds_kinds,     only : wp, ik
   use meds_numerics,  only : smooth_max0
   use meds_constants, only : grav_head, rho_h2o, pi, tiny_num
   implicit none
   private

   !----- A PFT's hydraulic traits, from its PFT file's [hydraulics] block. ---------------!
   public :: plant_hydraulics_t
   !----- Tissue pressure-volume curve: read-offs and the solver's relations. -------------!
   public :: water_curve_t
   public :: tissue_psi_tlp, tissue_psi_to_rwc, tissue_rwc_to_psi
   public :: tissue_psi_to_water, tissue_capacitance, tissue_water_to_psi, tissue_capped_water
   public :: tissue_water_relations, tissue_rwc_to_unknown, tissue_water_to_rwc
   !----- Vulnerability curve and Kirchhoff potential. ------------------------------------!
   public :: conductance_retained_fraction, kirchhoff_potential, hydro_table_t, build_hydro_table
   !----- Xylem and roots: the conductances and the flows they carry. ---------------------!
   public :: xylem_max_conductance, xylem_sapflow
   public :: root_fraction_profile, rhizosphere_conductance, root_water_uptake, root_water_uptake_profile

   !----- One tissue's pressure-volume curve: the four traits every tissue-water routine   !
   !      reads, passed together with the tissue's biomass. -------------------------------!
   type :: water_curve_t
      real(wp) :: pi0           = 0.0_wp   !< [MPa] osmotic potential at full turgor (< 0)
      real(wp) :: elastic_mod   = 0.0_wp   !< [MPa] bulk elastic modulus
      real(wp) :: apoplast_frac = 0.0_wp   !< [-]   share of the saturated water outside the living cells
      real(wp) :: water_sat     = 0.0_wp   !< [kg H2O / kgC] water per unit biomass at saturation
   end type water_curve_t

   !----- One PFT's hydraulic traits, from its PFT file's [hydraulics] block               !
   !      (docs/science/plant_hydraulics.md, Parameters). The traits have no defaults: a   !
   !      PFT file gives every one, for every PFT. ----------------------------------------!
   type :: plant_hydraulics_t
      type(water_curve_t) :: leaf_curve, wood_curve      !< the tissues' pressure-volume curves
      real(wp) :: wood_psi50           = 0.0_wp  !< [MPa, < 0] xylem potential at 50 % loss of conductance
      real(wp) :: wood_kexp            = 0.0_wp  !< [-] vulnerability-curve shape
      real(wp) :: wood_kmax            = 0.0_wp  !< [kg/m/s/MPa] sapwood specific conductivity
      real(wp) :: vessel_curl          = 0.0_wp  !< [-] xylem path length per unit height
      real(wp) :: root_beta            = 0.0_wp  !< [-] the untruncated profile's share below the rooting depth (eq. (13))
      real(wp) :: root_depth_b1        = 0.0_wp  !< [m] rooting depth = b1 * height[m]^b2
      real(wp) :: root_depth_b2        = 0.0_wp  !< [-]
      real(wp) :: specific_root_length = 0.0_wp  !< [m/kgC] fine-root length per unit fine-root carbon
      real(wp) :: fine_root_radius     = 0.0_wp  !< [m] fine-root radius
   end type plant_hydraulics_t

   !----- [-] Floor on the symplastic relative water content in the read-offs; the         !
   !      solver's relations have none. ---------------------------------------------------!
   real(wp), parameter :: RWC_FLOOR = 1.0e-4_wp

   !----- The Kirchhoff lookup table for a general vulnerability exponent: G(r) on a       !
   !      uniform grid (docs/science/plant_hydraulics.md §3). A fixed-size record, so a    !
   !      PFT's parameters copy as plain data. --------------------------------------------!
   integer,  parameter :: HYDRO_TABLE_NTAB = 512        !< intervals over [0, HYDRO_TABLE_RMAX]
   real(wp), parameter :: HYDRO_TABLE_RMAX = 8.0_wp     !< [-] table domain in r = psi/psi50
   type :: hydro_table_t
      real(wp)    :: kexp  = 0.0_wp                       !< [-] vulnerability exponent the table holds; its tail uses it too
      real(wp)    :: r_max = HYDRO_TABLE_RMAX             !< [-] table domain [0, r_max]
      real(wp)    :: inv_h = 0.0_wp                       !< [-] 1 / grid spacing
      integer(ik) :: n     = int(HYDRO_TABLE_NTAB, ik)    !< number of intervals
      real(wp)    :: g(0:HYDRO_TABLE_NTAB) = 0.0_wp       !< g(j) = integral_0^{j*h} du/(1+u^kexp)
   end type hydro_table_t

   !----- Floors that keep a cohort's leaves joined to its wood however dry or small it    !
   !      is, so its equations stay solvable; neither moves water that matters. -----------!
   real(wp), parameter :: XYLEM_RETAINED_MIN = 1.0e-12_wp   !< [-] share of the intact conductance a fully embolised xylem keeps
   real(wp), parameter :: XYLEM_COND_MIN = 1.0e-15_wp       !< [kg/s/MPa per plant] intact conductance of a plant without sapwood
   !----- [MPa] Below this difference of potentials the Kirchhoff mean conductance is      !
   !      taken at the midpoint, where the difference quotient would lose its digits. -----!
   real(wp), parameter :: KIRCHHOFF_EPS = 1.0e-6_wp

   !----- [-] Floor of half-spacing / root radius in the single-root (Gardner)             !
   !      conductance: closer, the roots all but touch and the form fails. ----------------!
   real(wp), parameter :: MIN_SPACING_RATIO = 2.0_wp
   !----- [MPa] Width over which smooth_max0 rounds root uptake's corner. At night the wood  !
   !      settles towards its wettest layers, which then sit in the corner, where uptake's  !
   !      slope is small; a narrower corner makes Newton see those layers as almost         !
   !      unresponsive, overshoot the wood's potential and backtrack. This width is still   !
   !      10-100 times below the 0.1-1 MPa a root works across, and changes uptake only     !
   !      where the driving force is below it (docs/science/soil_plant_water.md eq. (1)). --!
   real(wp), parameter :: UPTAKE_RAMP_WIDTH = 1.0e-2_wp

contains

   !=======================================================================================!
   !     Tissue pressure-volume curve (Bartlett / Tyree-Hammel;                            !
   !     docs/science/plant_hydraulics.md §2, eqs. (3)-(6)): the read-offs diagnostics and !
   !     a store's set-up use, floored at RWC_FLOOR so they stay finite.                   !
   !=======================================================================================!

   !----- Turgor loss point [MPa] (Bartlett et al. 2012 eq. 1;                             !
   !      docs/science/plant_hydraulics.md eq. (3)). --------------------------------------!
   elemental real(wp) function tissue_psi_tlp(pi0, elastic_mod) result(psi_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
   end function tissue_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett et al. 2012 eq. 2;    !
   !      docs/science/plant_hydraulics.md eq. (3)). --------------------------------------!
   elemental real(wp) function tissue_rwc_tlp(pi0, elastic_mod) result(rwc_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
   end function tissue_rwc_tlp

   !----- Water potential [MPa] at symplastic relative water content rwc, floored at       !
   !      RWC_FLOOR (docs/science/plant_hydraulics.md eq. (4)). ---------------------------!
   elemental real(wp) function tissue_rwc_to_psi(rwc, pi0, elastic_mod) result(psi)
      real(wp), intent(in) :: rwc, pi0, elastic_mod
      real(wp) :: rwc_tlp, r
      rwc_tlp = tissue_rwc_tlp(pi0, elastic_mod)
      r       = max(rwc, RWC_FLOOR)
      if (r >= rwc_tlp) then
         psi = elastic_mod*(r - rwc_tlp) + pi0/r   ! turgor + osmotic
      else
         psi = pi0/r                               ! turgor lost
      end if
   end function tissue_rwc_to_psi

   !----- Symplastic relative water content at water potential psi, floored at RWC_FLOOR   !
   !      (docs/science/plant_hydraulics.md eq. (4)). -------------------------------------!
   elemental real(wp) function tissue_psi_to_rwc(psi, pi0, elastic_mod) result(rwc)
      real(wp), intent(in) :: psi, pi0, elastic_mod
      real(wp) :: psi_tlp, b, disc
      psi_tlp = tissue_psi_tlp(pi0, elastic_mod)
      if (psi >= psi_tlp) then                ! turgid
         !----- The larger root of elastic_mod R^2 - b R + pi0 = 0. -----------------------!
         b    = psi + elastic_mod + pi0
         disc = max(b*b - 4.0_wp*elastic_mod*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*elastic_mod)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, RWC_FLOOR)
   end function tissue_psi_to_rwc

   !----- Total tissue water [kg per plant] at potential psi: the symplast's saturated     !
   !      water times R, plus the apoplast's, a constant reservoir                         !
   !      (docs/science/plant_hydraulics.md eq. (5)). -------------------------------------!
   elemental real(wp) function tissue_psi_to_water(psi, curve, biomass) result(w)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, r
      w_sat = curve%water_sat*biomass
      r     = tissue_psi_to_rwc(psi, curve%pi0, curve%elastic_mod)
      w     = (1.0_wp - curve%apoplast_frac)*w_sat*r + curve%apoplast_frac*w_sat
   end function tissue_psi_to_water

   !----- Capacitance dW/dpsi [kg/MPa per plant] at potential psi                          !
   !      (docs/science/plant_hydraulics.md eqs. (5)-(6)). --------------------------------!
   elemental real(wp) function tissue_capacitance(psi, curve, biomass) result(c)
      real(wp),            intent(in) :: psi, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      associate (pi0 => curve%pi0, elastic_mod => curve%elastic_mod)
         psi_tlp   = tissue_psi_tlp(pi0, elastic_mod)
         w_sat_sym = (1.0_wp - curve%apoplast_frac)*curve%water_sat*biomass
         r         = tissue_psi_to_rwc(psi, pi0, elastic_mod)
         if (psi >= psi_tlp) then
            cr = 1.0_wp / (elastic_mod - pi0/(r*r))   ! dR/dpsi, turgid
         else
            cr = -(r*r)/pi0                           ! dR/dpsi, flaccid
         end if
      end associate
      c = w_sat_sym*cr
   end function tissue_capacitance

   !----- Water potential [MPa] of a tissue holding w [kg per plant], the exact inverse of !
   !      tissue_psi_to_water (docs/science/plant_hydraulics.md eq. (5)). -----------------!
   elemental real(wp) function tissue_water_to_psi(w, curve, biomass) result(psi)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      real(wp) :: w_sat, rwc
      w_sat = curve%water_sat*biomass
      rwc   = (w - curve%apoplast_frac*w_sat) / max((1.0_wp - curve%apoplast_frac)*w_sat, tiny_num)
      psi   = tissue_rwc_to_psi(rwc, curve%pi0, curve%elastic_mod)
   end function tissue_water_to_psi

   !----- Tissue water [kg per plant] capped at saturation, water_sat * biomass, for a     !
   !      store that lost biomass without its water; the caller sheds the excess (w -      !
   !      result) to the ground (docs/science/soil_plant_water.md §2). --------------------!
   elemental real(wp) function tissue_capped_water(w, curve, biomass) result(w_capped)
      real(wp),            intent(in) :: w, biomass
      type(water_curve_t), intent(in) :: curve
      w_capped = min(w, curve%water_sat*biomass)
   end function tissue_capped_water

   !=======================================================================================!
   !     Tissue water relations for the soil-plant water solver, without the read-offs'    !
   !     floors.                                                                           !
   !=======================================================================================!

   !----- A tissue at symplastic relative water content r: potential [MPa] and water       !
   !      [kg/m2 of ground] for nplant plants of `biomass` each, with their derivatives in !
   !      r. The pressure-volume curve without tissue_rwc_to_psi's floor: the solver keeps !
   !      r above zero (docs/science/plant_hydraulics.md eqs. (4)-(5)). -------------------!
   elemental subroutine tissue_water_relations(curve, biomass, nplant, r, psi, dpsi, w, dw)
      type(water_curve_t), intent(in)  :: curve
      real(wp),            intent(in)  :: biomass   !< [kgC per plant]
      real(wp),            intent(in)  :: nplant    !< [plants/m2]
      real(wp),            intent(in)  :: r         !< [-] symplastic relative water content
      real(wp),            intent(out) :: psi, dpsi, w, dw
      real(wp) :: r_tlp, w_sat
      associate (pi0 => curve%pi0, emod => curve%elastic_mod, apo => curve%apoplast_frac)
         r_tlp = tissue_rwc_tlp(pi0, emod)
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
   !      then starts the tissue from a potential instead (tissue_psi_to_rwc). ------------!
   elemental function tissue_water_to_rwc(curve, biomass, w) result(r)
      type(water_curve_t), intent(in) :: curve
      real(wp),            intent(in) :: biomass, w
      real(wp) :: r, w_sat
      w_sat = curve%water_sat * biomass
      r = (w / w_sat - curve%apoplast_frac) / (1.0_wp - curve%apoplast_frac)
   end function tissue_water_to_rwc

   !----- A tissue's unknown is its relative water content, kept strictly positive: zero   !
   !      is an infinitely dry tissue. ----------------------------------------------------!
   elemental function tissue_rwc_to_unknown(r) result(tau)
      real(wp), intent(in) :: r
      real(wp) :: tau
      real(wp), parameter :: R_START = 1.0e-6_wp
      tau = max(r, R_START)   ! clamp-ok: an emptied store starts just inside the curve
   end function tissue_rwc_to_unknown

   !=======================================================================================!
   !     Vulnerability curve and Kirchhoff potential (docs/science/plant_hydraulics.md §3, !
   !     eqs. (7)-(8)): closed forms for vulnerability exponents 1 and 2, otherwise the    !
   !     PFT's lookup table, filled by 7-point Gauss-Legendre panels.                      !
   !=======================================================================================!

   !----- Retained fraction of xylem conductance, 1 - PLC = 1 / (1 + r^kexp), r =          !
   !      psi/psi50, and one at or above zero potential (docs/science/plant_hydraulics.md  !
   !      eq. (7)). r = 0 is kept off the general power: under floating-point trapping a   !
   !      compiler may take 0**kexp through log(0). Exponent 2, tested as in               !
   !      has_closed_form, squares r by a product: the general power is a library call,    !
   !      and this fraction is the sapflow's slope at every evaluation of the water        !
   !      network. ------------------------------------------------------------------------!
   elemental real(wp) function conductance_retained_fraction(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (r <= 0.0_wp) then
         f = 1.0_wp
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         f = 1.0_wp / (1.0_wp + r * r)
      else
         f = 1.0_wp / (1.0_wp + r**kexp)
      end if
   end function conductance_retained_fraction

   !----- Kirchhoff (matric flux) potential [MPa], Phi(psi) = integral of the retained     !
   !      fraction from 0 to psi (docs/science/plant_hydraulics.md eq. (8)): closed form   !
   !      for kexp 1 and 2, the PFT's table otherwise. Above zero potential nothing is     !
   !      embolised, so Phi = psi there; a flat Phi would read a supersaturated tissue's   !
   !      xylem as blocked. ---------------------------------------------------------------!
   pure real(wp) function kirchhoff_potential(psi, psi50, kexp, table) result(phi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      if (psi > 0.0_wp) then
         phi = psi
      else if (has_closed_form(kexp)) then
         phi = kirchhoff_potential_closed(psi, psi50, kexp)
      else
         phi = kirchhoff_potential_table(psi, psi50, table)
      end if
   end function kirchhoff_potential

   !----- Phi(psi_a) - Phi(psi_b) [MPa], the Kirchhoff potential's difference across a     !
   !      path. For exponent 2 with both ends at or below zero it is one arctangent, psi50 !
   !      atan((r_a - r_b) / (1 + r_a r_b)), not two: half the cost, and it keeps the      !
   !      digits two close potentials share. ----------------------------------------------!
   pure real(wp) function kirchhoff_difference(psi_a, psi_b, psi50, kexp, table) result(dphi)
      real(wp),            intent(in) :: psi_a, psi_b, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      real(wp) :: r_a, r_b
      if (psi_a <= 0.0_wp .and. psi_b <= 0.0_wp .and. abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         r_a  = psi_a / psi50
         r_b  = psi_b / psi50
         dphi = psi50 * atan((r_a - r_b) / (1.0_wp + r_a * r_b))
      else
         dphi = kirchhoff_potential(psi_a, psi50, kexp, table) - kirchhoff_potential(psi_b, psi50, kexp, table)
      end if
   end function kirchhoff_difference

   !----- dPhi/dpsi of exactly the potential kirchhoff_potential computes (the retained    !
   !      fraction, the slope of the table's linear piece, or the tail's), so the solver's !
   !      Jacobian is the derivative of the flow it moves. --------------------------------!
   pure real(wp) function kirchhoff_slope(psi, psi50, kexp, table) result(dphi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      real(wp) :: r
      integer  :: idx
      if (psi > 0.0_wp) then
         dphi = 1.0_wp
      else if (has_closed_form(kexp)) then
         dphi = conductance_retained_fraction(psi, psi50, kexp)
      else
         r = psi / psi50
         if (r >= table%r_max) then
            dphi = r ** (-table%kexp)
         else
            idx  = min(int(r * table%inv_h), table%n - 1)
            dphi = (table%g(idx + 1) - table%g(idx)) * table%inv_h
         end if
      end if
   end function kirchhoff_slope

   !----- d2Phi/dpsi2, the slope of kirchhoff_slope: d/dpsi of 1/(1 + r^kexp) for the      !
   !      closed forms; zero on the table's linear pieces and above zero potential; the    !
   !      tail's slope beyond the table. r = 0 is guarded as in                            !
   !      conductance_retained_fraction. --------------------------------------------------!
   pure real(wp) function kirchhoff_curvature(psi, psi50, kexp, table) result(d2phi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      real(wp) :: r
      r = psi / psi50
      if (r <= 0.0_wp) then
         d2phi = 0.0_wp
      else if (has_closed_form(kexp)) then
         d2phi = -kexp * r ** (kexp - 1.0_wp) / (psi50 * (1.0_wp + r ** kexp) ** 2)
      else if (r >= table%r_max) then
         d2phi = -table%kexp * r ** (-table%kexp - 1.0_wp) / psi50
      else
         d2phi = 0.0_wp
      end if
   end function kirchhoff_curvature

   !----- Whether the vulnerability exponent is 1 or 2, where Phi has a closed form. ------!
   pure logical function has_closed_form(kexp)
      real(wp), intent(in) :: kexp
      has_closed_form = abs(kexp - 1.0_wp) < 1.0e-9_wp .or. abs(kexp - 2.0_wp) < 1.0e-9_wp
   end function has_closed_form

   !----- The closed forms: psi50 ln(1 + r) for kexp 1, psi50 atan(r) for kexp 2. ---------!
   pure real(wp) function kirchhoff_potential_closed(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else
         phi = psi50 * atan(r)
      end if
   end function kirchhoff_potential_closed

   !----- Kirchhoff potential [MPa] read from the table by linear interpolation; psi50 is  !
   !      a run-time scale, so only kexp sets the table. Beyond the table the analytic     !
   !      tail is added to its last node: extrapolating the last interval could turn a     !
   !      dried leaf's potential positive. ------------------------------------------------!
   pure real(wp) function kirchhoff_potential_table(psi, psi50, tab) result(phi)
      real(wp),            intent(in) :: psi, psi50
      type(hydro_table_t), intent(in) :: tab
      real(wp) :: r, x, t
      integer  :: idx
      r   = max(psi/psi50, 0.0_wp)
      if (r >= tab%r_max) then
         phi = psi50 * ( tab%g(HYDRO_TABLE_NTAB) + kirchhoff_tail(tab%r_max, r, tab%kexp) )
         return
      end if
      x   = r * tab%inv_h
      idx = min(int(x), HYDRO_TABLE_NTAB - 1)
      t   = x - real(idx, wp)
      phi = psi50 * ( tab%g(idx)*(1.0_wp - t) + tab%g(idx+1)*t )
   end function kirchhoff_potential_table

   !----- Fill the table, one 7-point panel per interval added to the node below, so it    !
   !      rises monotonically. A subroutine, not an array-valued function                  !
   !      (docs/building.md, Compilers). --------------------------------------------------!
   pure subroutine build_hydro_table(tab, kexp)
      type(hydro_table_t), intent(out) :: tab
      real(wp),            intent(in)  :: kexp
      real(wp) :: h
      integer  :: j
      tab%kexp  = kexp
      tab%r_max = HYDRO_TABLE_RMAX
      tab%n     = int(HYDRO_TABLE_NTAB, ik)
      h         = HYDRO_TABLE_RMAX / real(HYDRO_TABLE_NTAB, wp)
      tab%inv_h = 1.0_wp / h
      tab%g(0) = 0.0_wp
      do j = 1, HYDRO_TABLE_NTAB
         tab%g(j) = tab%g(j-1) + kirchhoff_panel(real(j-1, wp)*h, real(j, wp)*h, kexp)
      end do
   end subroutine build_hydro_table

   !----- integral_a^b du/(1+u^kexp) on one panel (0 <= a <= b) by 7-point Gauss-Legendre, !
   !      exact for polynomials up to degree 13. The sum is written out rather than handed !
   !      a function (docs/building.md, no procedure arguments); u = 0 is guarded as in    !
   !      conductance_retained_fraction. --------------------------------------------------!
   pure real(wp) function kirchhoff_panel(a, b, kexp) result(g)
      real(wp), intent(in) :: a, b, kexp
      real(wp), parameter :: NODE(7)   = [ -0.9491079123427585_wp, -0.7415311855993945_wp,     &
                                           -0.4058451513773972_wp,  0.0000000000000000_wp,     &
                                            0.4058451513773972_wp,  0.7415311855993945_wp,     &
                                            0.9491079123427585_wp ]
      real(wp), parameter :: WEIGHT(7) = [  0.1294849661688697_wp,  0.2797053914892766_wp,     &
                                            0.3818300505051189_wp,  0.4179591836734694_wp,     &
                                            0.3818300505051189_wp,  0.2797053914892766_wp,     &
                                            0.1294849661688697_wp ]
      real(wp)    :: mid, half, u, acc
      integer(ik) :: k
      mid  = 0.5_wp * (a + b)                    ! the nodes on [-1,1], mapped onto [a,b]
      half = 0.5_wp * (b - a)
      acc  = 0.0_wp
      do k = 1_ik, 7_ik
         u = mid + half * NODE(k)
         if (u <= 0.0_wp) then
            acc = acc + WEIGHT(k)
         else
            acc = acc + WEIGHT(k) / (1.0_wp + u**kexp)
         end if
      end do
      g = half * acc
   end function kirchhoff_panel

   !----- integral_a^b du/(1+u^kexp) beyond the table (a >= HYDRO_TABLE_RMAX), where the   !
   !      integrand is u^-kexp to leading order: positive and, for kexp > 1, bounded as b  !
   !      grows, so the potential keeps rising toward its asymptote. ----------------------!
   pure real(wp) function kirchhoff_tail(a, b, kexp) result(g)
      real(wp), intent(in) :: a, b, kexp
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         g = log(b / a)
      else
         g = (b**(1.0_wp - kexp) - a**(1.0_wp - kexp)) / (1.0_wp - kexp)
      end if
   end function kirchhoff_tail

   !=======================================================================================!
   !     Xylem (docs/science/plant_hydraulics.md eqs. (9)-(10)).                           !
   !=======================================================================================!

   !----- One plant's intact-xylem conductance [kg/s/MPa]: the sapwood's specific          !
   !      conductivity over its cross-section and its path, height x vessel_curl, floored  !
   !      for a plant with no sapwood yet (docs/science/plant_hydraulics.md eq. (10)). ----!
   elemental function xylem_max_conductance(wood_kmax, sap_area, height, vessel_curl) result(k)
      real(wp), intent(in) :: wood_kmax    !< [kg/m/s/MPa] sapwood specific conductivity
      real(wp), intent(in) :: sap_area     !< [m2] sapwood cross-section
      real(wp), intent(in) :: height       !< [m]
      real(wp), intent(in) :: vessel_curl  !< [-] path length per unit height
      real(wp) :: k
      k = max(wood_kmax * sap_area / max(height * vessel_curl, tiny_num), XYLEM_COND_MIN)
   end function xylem_max_conductance

   !----- Sapflow [kg/m2/s] from a cohort's wood (psi_w) up to its leaves (psi_l), `grav`  !
   !      [MPa] of head above it, and its derivatives. k_cond [kg/s/MPa per m2] is the     !
   !      intact conductance of all the cohort's plants; the vulnerability curve reduces   !
   !      it through the Kirchhoff mean over the two ends, taken at the midpoint when they !
   !      nearly meet (docs/science/plant_hydraulics.md eqs. (7)-(9)). --------------------!
   pure subroutine xylem_sapflow(k_cond, psi_w, psi_l, grav, psi50, kexp, table, s, ds_dpsi_w, ds_dpsi_l)
      real(wp),            intent(in)  :: k_cond, psi_w, psi_l, grav, psi50, kexp
      type(hydro_table_t), intent(in)  :: table
      real(wp),            intent(out) :: s, ds_dpsi_w, ds_dpsi_l
      real(wp) :: dpsi, keff, dk_dw, dk_dl, mid, drive
      dpsi = psi_w - psi_l
      if (abs(dpsi) > KIRCHHOFF_EPS) then
         keff  = k_cond * kirchhoff_difference(psi_w, psi_l, psi50, kexp, table) / dpsi
         dk_dw = (k_cond * kirchhoff_slope(psi_w, psi50, kexp, table) - keff) / dpsi
         dk_dl = (keff - k_cond * kirchhoff_slope(psi_l, psi50, kexp, table)) / dpsi
      else
         mid   = 0.5_wp * (psi_w + psi_l)
         keff  = k_cond * kirchhoff_slope(mid, psi50, kexp, table)
         dk_dw = 0.5_wp * k_cond * kirchhoff_curvature(mid, psi50, kexp, table)
         dk_dl = dk_dw
      end if
      keff = keff + k_cond * XYLEM_RETAINED_MIN
      drive = dpsi - grav
      s         = keff * drive
      ds_dpsi_w = dk_dw * drive + keff
      ds_dpsi_l = dk_dl * drive - keff
   end subroutine xylem_sapflow

   !=======================================================================================!
   !     Roots (docs/science/plant_hydraulics.md eqs. (13)-(15)).                          !
   !=======================================================================================!

   !----- The share of a cohort's roots in the soil layer between depths z_top and z_bot,  !
   !      beta^(z_top/D) - beta^(z_bot/D) with D = root_depth and the depths clamped to    !
   !      [0, D], so a layer below the rooting depth gets none. Over [0, D] the shares sum !
   !      to 1 - beta, which the caller divides by (docs/science/plant_hydraulics.md eq.   !
   !      (13)). --------------------------------------------------------------------------!
   elemental real(wp) function root_fraction_profile(root_beta, root_depth, z_top, z_bot) result(frac)
      real(wp), intent(in) :: root_beta   !< [-]  root-profile decay (0,1); smaller => shallower
      real(wp), intent(in) :: root_depth  !< [m]  maximum rooting depth (> 0)
      real(wp), intent(in) :: z_top       !< [m]  depth of the layer top    (>= 0)
      real(wp), intent(in) :: z_bot       !< [m]  depth of the layer bottom (> z_top)
      real(wp) :: inv_d, a, b
      inv_d = 1.0_wp / max(root_depth, tiny_num)
      a = min(max(z_top, 0.0_wp), root_depth) * inv_d
      b = min(max(z_bot, 0.0_wp), root_depth) * inv_d
      frac = root_beta**a - root_beta**b
   end function root_fraction_profile

   !----- A cohort's rhizosphere conductance in one soil layer per unit soil conductivity  !
   !      [kg/s/MPa per m2 of ground, per m/s]: steady radial flow to a single root        !
   !      (Gardner 1960), 2 pi L / ln(r_half / r_root) for each plant, made a potential    !
   !      gradient by rho/grav_head. It depends on the roots alone, so it is set once a    !
   !      step (docs/science/plant_hydraulics.md eq. (14)). -------------------------------!
   elemental function rhizosphere_conductance(nplant, root_length, half_spacing, root_radius)      &
                      result(g_rhizo)
      real(wp), intent(in) :: nplant        !< [plants/m2]
      real(wp), intent(in) :: root_length   !< [m per plant] fine-root length in the layer
      real(wp), intent(in) :: half_spacing  !< [m] half the distance between the layer's roots
      real(wp), intent(in) :: root_radius   !< [m] absorbing-root radius
      real(wp) :: g_rhizo
      g_rhizo = nplant * 2.0_wp * pi * root_length * rho_h2o                                      &
               / (grav_head * log(max(half_spacing / root_radius, MIN_SPACING_RATIO)))
   end function rhizosphere_conductance

   !----- A cohort's root water uptake [kg/m2/s] from one soil layer: g_rhizo K_soil times !
   !      the driving force psi_soil + grav_head z_soil - psi_wood (z_soil < 0), through   !
   !      the ramp, so a root never loses water to drier soil                              !
   !      (docs/science/plant_hydraulics.md eq. (15)). Outside its rounded corner the ramp !
   !      is the force itself or zero, written out here: this runs for every cohort and    !
   !      layer at every evaluation of the water network, and smooth_max0 is a call into   !
   !      another module. -----------------------------------------------------------------!
   elemental subroutine root_water_uptake(g_rhizo, k_soil, psi_soil, z_soil, psi_wood, u,             &
                                          du_dpsi_soil, du_dpsi_wood, du_dk)
      real(wp), intent(in)  :: g_rhizo, k_soil, psi_soil, z_soil, psi_wood
      real(wp), intent(out) :: u, du_dpsi_soil, du_dpsi_wood, du_dk
      real(wp) :: force, drive, dramp
      force = psi_soil + grav_head * z_soil - psi_wood
      if (force >= UPTAKE_RAMP_WIDTH) then
         drive = force
         dramp = 1.0_wp
      else if (force <= 0.0_wp) then
         drive = 0.0_wp
         dramp = 0.0_wp
      else
         call smooth_max0(force, UPTAKE_RAMP_WIDTH, drive, dramp)
      end if
      u            = g_rhizo * k_soil * drive
      du_dpsi_soil = g_rhizo * k_soil * dramp
      du_dpsi_wood = -du_dpsi_soil
      du_dk        = g_rhizo * drive
   end subroutine root_water_uptake

   !----- root_water_uptake for every layer of a cohort's profile, nothing from a layer    !
   !      without roots (g_rhizo <= 0). The wood's derivative is minus the soil's, so it   !
   !      is not returned. ----------------------------------------------------------------!
   pure subroutine root_water_uptake_profile(nsl, g_rhizo, k_soil, psi_soil, z_soil, psi_wood, u,      &
                                             du_dpsi_soil, du_dk)
      integer(ik), intent(in)  :: nsl
      real(wp),    intent(in)  :: g_rhizo(nsl), k_soil(nsl), psi_soil(nsl), z_soil(nsl), psi_wood
      real(wp),    intent(out) :: u(nsl), du_dpsi_soil(nsl), du_dk(nsl)
      real(wp)    :: du_dpsi_wood
      integer(ik) :: k
      do k = 1_ik, nsl
         if (g_rhizo(k) <= 0.0_wp) then
            u(k) = 0.0_wp ; du_dpsi_soil(k) = 0.0_wp ; du_dk(k) = 0.0_wp
            cycle
         end if
         call root_water_uptake(g_rhizo(k), k_soil(k), psi_soil(k), z_soil(k), psi_wood, u(k),          &
                                du_dpsi_soil(k), du_dpsi_wood, du_dk(k))
      end do
   end subroutine root_water_uptake_profile

end module meds_plant_hydr_lib
