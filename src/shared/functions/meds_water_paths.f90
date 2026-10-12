! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_water_paths -- how fast water moves along each kind of path in the soil-plant water     !
! network (MEDS_SPAC_WATER_NETWORK_PLAN.md §4.2): between two soil layers or the pond and the  !
! top layer (Darcy), from a soil layer into a cohort's roots (the rhizosphere), and up a        !
! cohort's xylem from its wood to its leaves; and how much a leaf can transpire when it is      !
! short of water.                                                                            !
!                                                                                          !
! Every routine returns the flow [kg per m2 of ground per s] and its derivatives with respect to  !
! the potentials at its two ends [MPa] and, where it carries one, the soil conductivity, so the  !
! solver builds its Jacobian from the same numbers it moves water with. Every switch is a smooth  !
! ramp (smooth_max0), so Newton never meets a corner it cannot see.                               !
!==========================================================================================!
module meds_water_paths
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav_head, rho_h2o, pi, tiny_num
   use meds_hydr_lib,  only : plc_retained, flux_potential, flux_potential_lin, hydro_table_t
   implicit none
   private

   public :: ROOT_FLOW_WIDTH, TRANSP_CLOSE_WIDTH
   public :: smooth_max0, darcy_flow, root_conductance_geometry, root_uptake_flow, root_uptake_profile
   public :: xylem_max_conductance, xylem_flow, supply_limited_transpiration

   !----- [MPa] Width of the ramp that keeps water from flowing back out of the roots: below a     !
   !      zero driving force the root takes nothing (no hydraulic redistribution, plan D2); above  !
   !      this width it takes the full Gardner flow. Small against the 0.1-1 MPa a root works      !
   !      across. -------------------------------------------------------------------------------!
   real(wp), parameter :: ROOT_FLOW_WIDTH = 1.0e-3_wp
   !----- [-] Width, in relative water content, of the corner where a leaf's transpiration starts  !
   !      to follow its remaining water (supply_limited_transpiration). -----------------------!
   real(wp), parameter :: TRANSP_CLOSE_WIDTH = 2.0e-2_wp
   !----- Xylem conductance kept when the vulnerability curve has taken the rest, as a share of   !
   !      the intact xylem's: it keeps a fully embolised cohort's leaf joined to its wood, so its   !
   !      equations stay solvable, and moves no water that matters. --------------------------!
   real(wp), parameter :: XYLEM_RETAINED_MIN = 1.0e-12_wp
   !----- [kg/s/MPa per plant] The least maximum conductance a plant is given, for a cohort with no  !
   !      sapwood yet: the same guard the plant solver used. -------------------------------------!
   real(wp), parameter :: XYLEM_COND_MIN = 1.0e-15_wp
   !----- Root spacing below which the single-root (Gardner) form no longer holds: closer, the    !
   !      roots all but touch. -------------------------------------------------------------!
   real(wp), parameter :: MIN_SPACING_RATIO = 2.0_wp
   !----- [MPa] Below this difference of potentials the Kirchhoff mean is taken at the midpoint,   !
   !      where the difference quotient would lose its digits. ---------------------------------!
   real(wp), parameter :: KIRCHHOFF_EPS = 1.0e-6_wp

contains

   !----- max(x, 0) with its corner rounded over [0, width]: exactly zero below zero and exactly x  !
   !      above the width, joined by the cubic 2x^2/w - x^3/w^2, which matches both in value and     !
   !      slope and rises monotonically. Outside the corner the flow is the unrounded one. ---------!
   elemental subroutine smooth_max0(x, width, f, df)
      real(wp), intent(in)  :: x, width
      real(wp), intent(out) :: f, df
      real(wp) :: t
      if (x <= 0.0_wp) then
         f = 0.0_wp ; df = 0.0_wp
      else if (x < width) then
         t  = x / width
         f  = x * t * (2.0_wp - t)
         df = t * (4.0_wp - 3.0_wp * t)
      else
         f = x ; df = 1.0_wp
      end if
   end subroutine smooth_max0

   !----- Darcy flow [kg/m2/s] from node a to node b, positive from a to b. The face conductivity is !
   !      the thickness-weighted geometric mean K_a^(1-w_b) K_b^w_b (ED2's rule for soil layers;     !
   !      the pond counts as a saturated layer as thick as the top one, w_b = 1/2). The driving      !
   !      force is the difference of potentials plus the drop in elevation z_a - z_b [m], over the   !
   !      node distance [m]. Potentials in MPa, conductivities in m/s. ----------------------------!
   elemental subroutine darcy_flow(k_a, k_b, w_b, psi_a, psi_b, z_drop, dist, q,                    &
                                   dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b)
      real(wp), intent(in)  :: k_a, k_b, w_b, psi_a, psi_b, z_drop, dist
      real(wp), intent(out) :: q, dq_dpsi_a, dq_dpsi_b, dq_dk_a, dq_dk_b
      real(wp) :: k_face, drive
      k_face = k_a ** (1.0_wp - w_b) * k_b ** w_b
      drive  = ((psi_a - psi_b) / grav_head + z_drop) / dist
      q      = rho_h2o * k_face * drive
      dq_dpsi_a =  rho_h2o * k_face / (grav_head * dist)
      dq_dpsi_b = -dq_dpsi_a
      dq_dk_a   = (1.0_wp - w_b) * q / max(k_a, tiny_num)
      dq_dk_b   = w_b * q / max(k_b, tiny_num)
   end subroutine darcy_flow

   !----- The geometric part of a layer's soil-to-root conductance, per m2 of ground [kg/s/MPa per  !
   !      m/s of soil conductivity]: steady radial flow to a single root (Gardner 1960),            !
   !      2 pi L / ln(r_half / r_root) for each of nplant plants with root length L in the layer.   !
   !      Times the soil's conductivity in m/s (made a potential gradient by rho/grav_head), it is  !
   !      meds_plant_hydraulics' rhizosphere_cond for the whole cohort. -------------------------!
   elemental function root_conductance_geometry(nplant, root_length, half_spacing, root_radius)       &
                      result(g_geom)
      real(wp), intent(in) :: nplant        !< [plants/m2]
      real(wp), intent(in) :: root_length   !< [m per plant] fine-root length in the layer
      real(wp), intent(in) :: half_spacing  !< [m] half the distance between the layer's roots
      real(wp), intent(in) :: root_radius   !< [m] absorbing-root radius
      real(wp) :: g_geom
      g_geom = nplant * 2.0_wp * pi * root_length * rho_h2o                                       &
               / (grav_head * log(max(half_spacing / root_radius, MIN_SPACING_RATIO)))
   end function root_conductance_geometry

   !----- Uptake [kg/m2/s] from a soil layer into a cohort's wood: conductance g_geom K times the     !
   !      driving force psi_soil + grav_head z_soil - psi_wood (z_soil the layer's depth, < 0), passed !
   !      through the ramp so a root never loses water to drier soil. ----------------------------!
   elemental subroutine root_uptake_flow(g_geom, k_soil, psi_soil, z_soil, psi_wood, u,               &
                                         du_dpsi_soil, du_dpsi_wood, du_dk)
      real(wp), intent(in)  :: g_geom, k_soil, psi_soil, z_soil, psi_wood
      real(wp), intent(out) :: u, du_dpsi_soil, du_dpsi_wood, du_dk
      real(wp) :: drive, dramp
      call smooth_max0(psi_soil + grav_head * z_soil - psi_wood, ROOT_FLOW_WIDTH, drive, dramp)
      u            = g_geom * k_soil * drive
      du_dpsi_soil = g_geom * k_soil * dramp
      du_dpsi_wood = -du_dpsi_soil
      du_dk        = g_geom * drive
   end subroutine root_uptake_flow

   !----- A cohort's uptake from every layer of its root profile in one call: root_uptake_flow for   !
   !      each layer it roots in, and nothing from a layer without roots (g_geom <= 0). The wood's   !
   !      derivative is minus the soil's, so it is not returned. ----------------------------------!
   pure subroutine root_uptake_profile(nsl, g_geom, k_soil, psi_soil, z_soil, psi_wood, u,           &
                                       du_dpsi_soil, du_dk)
      integer(ik), intent(in)  :: nsl
      real(wp),    intent(in)  :: g_geom(nsl), k_soil(nsl), psi_soil(nsl), z_soil(nsl), psi_wood
      real(wp),    intent(out) :: u(nsl), du_dpsi_soil(nsl), du_dk(nsl)
      real(wp)    :: du_dpsi_wood
      integer(ik) :: k
      do k = 1_ik, nsl
         if (g_geom(k) <= 0.0_wp) then
            u(k) = 0.0_wp ; du_dpsi_soil(k) = 0.0_wp ; du_dk(k) = 0.0_wp
            cycle
         end if
         call root_uptake_flow(g_geom(k), k_soil(k), psi_soil(k), z_soil(k), psi_wood, u(k),               &
                               du_dpsi_soil(k), du_dpsi_wood, du_dk(k))
      end do
   end subroutine root_uptake_profile

   !----- One plant's intact-xylem conductance [kg/s/MPa]: the sapwood's specific conductivity over  !
   !      its cross-section and its path, height x vessel_curl. A plant with no sapwood yet keeps the  !
   !      least conductance, so its leaf stays joined to its wood. -------------------------------!
   elemental function xylem_max_conductance(wood_kmax, sap_area, height, vessel_curl) result(k)
      real(wp), intent(in) :: wood_kmax    !< [kg/m/s/MPa] sapwood specific conductivity
      real(wp), intent(in) :: sap_area     !< [m2] sapwood cross-section
      real(wp), intent(in) :: height       !< [m]
      real(wp), intent(in) :: vessel_curl  !< [-] path length per unit height
      real(wp) :: k
      k = max(wood_kmax * sap_area / max(height * vessel_curl, tiny_num), XYLEM_COND_MIN)
   end function xylem_max_conductance

   !----- Sapflow [kg/m2/s] from a cohort's wood (potential psi_w) to its leaves (psi_l), and its   !
   !      derivatives. k_cond [kg/s/MPa per m2] is the intact conductance of all the cohort's plants; !
   !      the vulnerability curve enters through the Kirchhoff mean over the two ends,                  !
   !      k (Phi(psi_w) - Phi(psi_l)) / (psi_w - psi_l), and the leaves sit `grav` [MPa] of           !
   !      hydrostatic head above the wood. Phi is closed-form for vulnerability exponents 1 and 2 and  !
   !      read from the PFT's table otherwise, as in the plant solver; the derivative of Phi is the     !
   !      retained fraction itself. Near equal potentials the mean is taken at the midpoint, value and  !
   !      slopes, where the quotient would lose its digits. ---------------------------------------!
   pure subroutine xylem_flow(k_cond, psi_w, psi_l, grav, psi50, kexp, table, s, ds_dpsi_w, ds_dpsi_l)
      real(wp),            intent(in)  :: k_cond, psi_w, psi_l, grav, psi50, kexp
      type(hydro_table_t), intent(in)  :: table
      real(wp),            intent(out) :: s, ds_dpsi_w, ds_dpsi_l
      real(wp) :: dpsi, keff, dk_dw, dk_dl, mid, drive
      dpsi = psi_w - psi_l
      if (abs(dpsi) > KIRCHHOFF_EPS) then
         keff  = k_cond * (kirchhoff_potential(psi_w, psi50, kexp, table)                          &
                          - kirchhoff_potential(psi_l, psi50, kexp, table)) / dpsi
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
   end subroutine xylem_flow

   !----- The Kirchhoff potential Phi(psi) = integral of the retained fraction from 0 to psi [MPa]:  !
   !      closed-form for vulnerability exponents 1 and 2, the PFT's lookup table otherwise -- the   !
   !      choice the plant solver's edge_cond makes. Above zero potential nothing is embolised, so  !
   !      Phi rises one for one (meds_hydr_lib holds it flat there, which reads a supersaturated    !
   !      tissue's xylem as partly blocked). -------------------------------------------------!
   pure real(wp) function kirchhoff_potential(psi, psi50, kexp, table) result(phi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      if (psi > 0.0_wp) then
         phi = psi
      else if (abs(kexp - 1.0_wp) < 1.0e-9_wp .or. abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = flux_potential(psi, psi50, kexp)
      else
         phi = flux_potential_lin(psi, psi50, table)
      end if
   end function kirchhoff_potential

   !----- dPhi/dpsi, the slope of exactly the potential kirchhoff_potential computes: the retained  !
   !      fraction for the closed forms, the slope of the table's linear piece (or of its analytic   !
   !      tail beyond the table) otherwise, so the solver's Jacobian is the derivative of the flow  !
   !      it moves. --------------------------------------------------------------------------!
   pure real(wp) function kirchhoff_slope(psi, psi50, kexp, table) result(dphi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      real(wp) :: r
      integer  :: idx
      if (psi > 0.0_wp) then
         dphi = 1.0_wp
      else if (abs(kexp - 1.0_wp) < 1.0e-9_wp .or. abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         dphi = plc_retained(psi, psi50, kexp)
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

   !----- d2Phi/dpsi2, the slope of kirchhoff_slope: for the closed forms the slope of the retained  !
   !      fraction 1/(1 + r^kexp), r = psi/psi50; zero on the table's linear pieces and above zero    !
   !      potential; the tail's own slope beyond the table. r = 0 is kept off the general power for   !
   !      the reason plc_retained keeps it off (0**x through log(0) on strict compilers). ----------!
   pure real(wp) function kirchhoff_curvature(psi, psi50, kexp, table) result(d2phi)
      real(wp),            intent(in) :: psi, psi50, kexp
      type(hydro_table_t), intent(in) :: table
      real(wp) :: r
      r = psi / psi50
      if (r <= 0.0_wp) then
         d2phi = 0.0_wp
      else if (abs(kexp - 1.0_wp) < 1.0e-9_wp .or. abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         d2phi = -kexp * r ** (kexp - 1.0_wp) / (psi50 * (1.0_wp + r ** kexp) ** 2)
      else if (r >= table%r_max) then
         d2phi = -table%kexp * r ** (-table%kexp - 1.0_wp) / psi50
      else
         d2phi = 0.0_wp
      end if
   end function kirchhoff_curvature

   !----- Transpiration [kg/m2/s] of a leaf whose stomata ask for e_demand, given its symplastic   !
   !      relative water content r (plan §12.2 C1). Above r_close -- the water content at which the   !
   !      fast stomatal limb has shut the stomata -- the demand stands. Below it, a leaf with closed  !
   !      stomata loses its residual water in proportion to the water it has left, r/r_close, so an   !
   !      emptying leaf stops transpiring before it empties and every step has a solution. The       !
   !      corner is rounded over [r_close, r_close + TRANSP_CLOSE_WIDTH]. --------------------------!
   elemental subroutine supply_limited_transpiration(e_demand, r, r_close, e, de_dr)
      real(wp), intent(in)  :: e_demand, r, r_close
      real(wp), intent(out) :: e, de_dr
      real(wp) :: excess, dexcess, share
      call smooth_max0(r - r_close, TRANSP_CLOSE_WIDTH, excess, dexcess)
      share = (r - excess) / r_close
      e     = e_demand * share
      de_dr = e_demand * (1.0_wp - dexcess) / r_close
   end subroutine supply_limited_transpiration

end module meds_water_paths
