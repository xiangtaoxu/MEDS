! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_hydraulics -- how water moves through a cohort's plants, and the              !
! plant-water facts the fast loop needs about each cohort (docs/science/                   !
! plant_hydraulics.md §2-§4): the crown gravity head, the leaf node's water potential,     !
! the rooting depth and root profile (eqs. (12)-(13)), the sapflow through the             !
! Kirchhoff mean of the vulnerability curve (eqs. (7)-(9)), root water uptake              !
! (eq. (15)), and the leaves' supply limit on transpiration (docs/science/                 !
! soil_plant_water.md eq. (3)). Each flow comes with its derivatives in the water          !
! network's unknowns. The tissue and vulnerability curves and the conductances' laws       !
! are meds_plant_hydr_lib's; each PFT's Kirchhoff lookup is built at load                  !
! (meds_pft_derive); the integrator solves the network.                                    !
!==========================================================================================!
module meds_plant_hydraulics
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : grav_head, tiny_num
   use meds_pft_params,     only : hydraulics_pft_t, kirchhoff_table_t, KIRCHHOFF_TABLE_NTAB
   use meds_plant_hydr_lib, only : root_fraction_profile, tissue_water_to_psi, conductance_retained_fraction
   use meds_numerics,       only : smooth_max0
   implicit none
   private

   public :: leaf_node_psi, crown_gravity_head
   public :: cohort_root_profile
   public :: xylem_sapflow, root_water_uptake, root_water_uptake_profile
   public :: supply_limited_transpiration
   public :: cohort_root_depth, kirchhoff_potential   ! used only by the tests

   !----- [-] Width, in relative water content, over which smooth_max0 rounds the leaves'  !
   !      supply limit (docs/science/soil_plant_water.md eq. (3)). ------------------------!
   real(wp), parameter :: TRANSP_CLOSE_WIDTH = 2.0e-2_wp

   !----- [-] Share of the intact conductance a fully embolised xylem keeps, so a cohort's     !
   !      leaves stay joined to its wood however dry it is; it moves no water that matters. --!
   real(wp), parameter :: XYLEM_RETAINED_MIN = 1.0e-12_wp
   !----- [MPa] Below this difference of potentials the Kirchhoff mean conductance is      !
   !      taken at the midpoint, where the difference quotient would lose its digits. -----!
   real(wp), parameter :: KIRCHHOFF_EPS = 1.0e-6_wp

   !----- [MPa] Width over which smooth_max0 rounds root uptake's corner. At night the wood  !
   !      settles towards its wettest layers, which then sit in the corner, where uptake's  !
   !      slope is small; a narrower corner makes Newton see those layers as almost         !
   !      unresponsive, overshoot the wood's potential and backtrack. This width is still   !
   !      10-100 times below the 0.1-1 MPa a root works across, and changes uptake only     !
   !      where the driving force is below it (docs/science/soil_plant_water.md eq. (1)). --!
   real(wp), parameter :: UPTAKE_RAMP_WIDTH = 1.0e-2_wp

contains

   !----- The water-potential drop from the wood node to the crown [MPa]: the hydrostatic head of  !
   !      the cohort's height. -------------------------------------------------------------------!
   pure real(wp) function crown_gravity_head(height) result(grav)
      real(wp), intent(in) :: height   !< [m] cohort height
      grav = grav_head*height
   end function crown_gravity_head

   !----- The leaf node's water potential [MPa], from the leaf water and its pressure-volume curve   !
   !      (docs/science/plant_hydraulics.md eqs. (4)-(5)). A leafless cohort has no curve to read: its !
   !      leaf node sits at the wood potential less the crown gravity head (§2). ---------------------!
   pure real(wp) function leaf_node_psi(leaf_water, bleaf, wood_water, bwood, height, p) result(psi_l)
      real(wp),               intent(in) :: leaf_water  !< [kg H2O/plant] leaf water
      real(wp),               intent(in) :: bleaf       !< [kgC/plant] leaf carbon
      real(wp),               intent(in) :: wood_water  !< [kg H2O/plant] wood (sapwood + fine-root) water
      real(wp),               intent(in) :: bwood       !< [kgC/plant] sapwood + fine-root carbon
      real(wp),               intent(in) :: height      !< [m] cohort height
      type(hydraulics_pft_t), intent(in) :: p
      if (p%leaf_curve%water_sat*bleaf > tiny_num) then
         psi_l = tissue_water_to_psi(leaf_water, p%leaf_curve, bleaf)
      else
         psi_l = tissue_water_to_psi(wood_water, p%wood_curve, bwood) - crown_gravity_head(height)
      end if
   end function leaf_node_psi

   !----- A cohort's rooting depth [m], docs/science/plant_hydraulics.md eq. (12): root_depth_b1 *     !
   !      height^root_depth_b2, capped at the soil column. --------------------------------------------!
   pure real(wp) function cohort_root_depth(height, hp, soil_depth) result(depth)
      real(wp),               intent(in) :: height      !< [m] cohort height
      type(hydraulics_pft_t), intent(in) :: hp          !< the cohort's PFT
      real(wp),               intent(in) :: soil_depth  !< [m] depth of the soil column's bottom
      depth = min(hp%root_depth_b1 * max(height, 0.0_wp)**hp%root_depth_b2, soil_depth)
   end function cohort_root_depth

   !----- A cohort's fine roots in each soil layer [-], summing to one: the exponential profile      !
   !      within its rooting depth (docs/science/plant_hydraulics.md eq. (13)). ---------------------!
   pure subroutine cohort_root_profile(height, hp, z_top, z_bot, frac)
      real(wp),               intent(in)  :: height      !< [m] cohort height
      type(hydraulics_pft_t), intent(in)  :: hp          !< the cohort's PFT
      real(wp),               intent(in)  :: z_top(:)    !< [m] depth of each layer's top (>= 0)
      real(wp),               intent(in)  :: z_bot(:)    !< [m] depth of each layer's bottom (> z_top)
      real(wp),               intent(out) :: frac(:)     !< [-] the cohort's fine roots in each layer
      real(wp) :: depth
      depth = cohort_root_depth(height, hp, z_bot(size(z_bot)))
      frac = root_fraction_profile(hp%root_beta, depth, z_top, z_bot) / (1.0_wp - hp%root_beta)
   end subroutine cohort_root_profile

   !=======================================================================================!
   !     The Kirchhoff mean of the vulnerability curve (docs/science/plant_hydraulics.md   !
   !     §3, eqs. (7)-(8)): closed forms for vulnerability exponents 1 and 2, otherwise    !
   !     the PFT's lookup table (meds_pft_derive).                                         !
   !=======================================================================================!

   !----- Kirchhoff (matric flux) potential [MPa], Phi(psi) = integral of the retained     !
   !      fraction from 0 to psi (docs/science/plant_hydraulics.md eq. (8)): closed form   !
   !      for kexp 1 and 2, the PFT's table otherwise. Above zero potential nothing is     !
   !      embolised, so Phi = psi there; a flat Phi would read a supersaturated tissue's   !
   !      xylem as blocked. ---------------------------------------------------------------!
   pure real(wp) function kirchhoff_potential(psi, psi50, kexp, table) result(phi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(kirchhoff_table_t), intent(in) :: table
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
      type(kirchhoff_table_t), intent(in) :: table
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
      type(kirchhoff_table_t), intent(in) :: table
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
      type(kirchhoff_table_t), intent(in) :: table
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
      type(kirchhoff_table_t), intent(in) :: tab
      real(wp) :: r, x, t
      integer  :: idx
      r   = max(psi/psi50, 0.0_wp)
      if (r >= tab%r_max) then
         phi = psi50 * ( tab%g(KIRCHHOFF_TABLE_NTAB) + kirchhoff_tail(tab%r_max, r, tab%kexp) )
         return
      end if
      x   = r * tab%inv_h
      idx = min(int(x), KIRCHHOFF_TABLE_NTAB - 1)
      t   = x - real(idx, wp)
      phi = psi50 * ( tab%g(idx)*(1.0_wp - t) + tab%g(idx+1)*t )
   end function kirchhoff_potential_table

   !----- integral_a^b du/(1+u^kexp) beyond the table (a >= KIRCHHOFF_TABLE_RMAX), where the   !
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
   !     The xylem's and the roots' flows (docs/science/plant_hydraulics.md eqs. (9),      !
   !     (15)), with their derivatives in the water network's unknowns.                   !
   !=======================================================================================!

   !----- Sapflow [kg/m2/s] from a cohort's wood (psi_w) up to its leaves (psi_l), `grav`  !
   !      [MPa] of head above it, and its derivatives. k_cond [kg/s/MPa per m2] is the     !
   !      intact conductance of all the cohort's plants; the vulnerability curve reduces   !
   !      it through the Kirchhoff mean over the two ends, taken at the midpoint when they !
   !      nearly meet (docs/science/plant_hydraulics.md eqs. (7)-(9)). --------------------!
   pure subroutine xylem_sapflow(k_cond, psi_w, psi_l, grav, psi50, kexp, table, s, ds_dpsi_w, ds_dpsi_l)
      real(wp),            intent(in)  :: k_cond, psi_w, psi_l, grav, psi50, kexp
      type(kirchhoff_table_t), intent(in)  :: table
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

   !----- A drying leaf's transpiration [kg/m2/s]: the stomata's demand while the leaf's     !
   !      symplastic relative water content r is above r_close (the closing potential's), the !
   !      demand times r/r_close below it, so the implicit step always has a solution and the !
   !      leaf never empties (docs/science/soil_plant_water.md eq. (3)). --------------------!
   elemental subroutine supply_limited_transpiration(e_demand, r, r_close, e, de_dr)
      real(wp), intent(in)  :: e_demand, r, r_close
      real(wp), intent(out) :: e, de_dr
      real(wp) :: excess, dexcess, share
      call smooth_max0(r - r_close, TRANSP_CLOSE_WIDTH, excess, dexcess)
      share = (r - excess) / r_close
      e     = e_demand * share
      de_dr = e_demand * (1.0_wp - dexcess) / r_close
   end subroutine supply_limited_transpiration

end module meds_plant_hydraulics
