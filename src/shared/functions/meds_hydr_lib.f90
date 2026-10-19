! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_hydr_lib -- how water moves into and through a plant: the conductances, the flows   !
! they carry and the leaves' supply limit on transpiration                                 !
! (docs/science/plant_hydraulics.md §3-§4, docs/science/soil_plant_water.md §1 and §4),    !
! and what they are built from: the vulnerability curve and its Kirchhoff potential, the   !
! sapwood conductivity of a wood density and the root profile. Every flow comes with its   !
! derivatives in the potentials at its two ends, so the solver (meds_spac_solver) builds   !
! its Jacobian from the numbers it moves water with. How much water a store holds is in    !
! meds_water_retention; flow between soil layers in meds_soil_water.                       !
!==========================================================================================!
module meds_hydr_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav_head, rho_h2o, pi, tiny_num
   implicit none
   private

   !----- Vulnerability curve, Kirchhoff potential and its lookup table, sapwood           !
   !      conductivity. -------------------------------------------------------------------!
   public :: plc_retained, flux_potential, wood_kmax_from_density
   public :: hydro_table_t, build_hydro_table, flux_potential_lin
   !----- The conductances and flows: rhizosphere, xylem, transpiration supply. -----------!
   public :: smooth_max0
   public :: rhizosphere_conductance, root_water_uptake, root_water_uptake_profile
   public :: xylem_max_conductance, xylem_sapflow, supply_limited_transpiration
   !----- The root profile. ---------------------------------------------------------------!
   public :: root_fraction_profile

   !----- Widths over which smooth_max0 rounds a one-way flow's corner                     !
   !      (docs/science/soil_plant_water.md eqs. (1) and (3)). ----------------------------!
   real(wp), parameter :: UPTAKE_RAMP_WIDTH = 1.0e-3_wp    !< [MPa] root uptake; small against the 0.1-1 MPa a root works across
   real(wp), parameter :: TRANSP_CLOSE_WIDTH = 2.0e-2_wp   !< [-] the leaves' supply limit, in relative water content
   !----- Floors that keep a cohort's leaves joined to its wood however dry or small it    !
   !      is, so its equations stay solvable; neither moves water that matters. -----------!
   real(wp), parameter :: XYLEM_RETAINED_MIN = 1.0e-12_wp   !< [-] share of the intact conductance a fully embolised xylem keeps
   real(wp), parameter :: XYLEM_COND_MIN = 1.0e-15_wp       !< [kg/s/MPa per plant] intact conductance of a plant without sapwood
   !----- [-] Floor of half-spacing / root radius in the single-root (Gardner)             !
   !      conductance: closer, the roots all but touch and the form fails. ----------------!
   real(wp), parameter :: MIN_SPACING_RATIO = 2.0_wp
   !----- [MPa] Below this difference of potentials the Kirchhoff mean conductance is      !
   !      taken at the midpoint, where the difference quotient would lose its digits. -----!
   real(wp), parameter :: KIRCHHOFF_EPS = 1.0e-6_wp

   !----- Sapwood specific conductivity against wood density, ln Ks = a + b rho: a         !
   !      lognormal fit to 72 neotropical dry-forest species (Xu et al. 2016, New Phytol   !
   !      212:80; docs/science/plant_hydraulics.md eq. (11)). -----------------------------!
   real(wp), parameter :: kmax_ln_a     =  2.348_wp   !< [ln(kg/m/s/MPa)] intercept
   real(wp), parameter :: kmax_ln_b     = -2.455_wp   !< [ln(kg/m/s/MPa) per g/cm3] slope
   real(wp), parameter :: kmax_ln_var   =  0.6186_wp  !< [-] residual variance of ln Ks
   real(wp), parameter :: kmax_rho_min  =  0.35_wp    !< [g/cm3] lowest wood density fitted; a lower one takes it
   real(wp), parameter :: kmax_rho_max  =  0.95_wp    !< [g/cm3] highest wood density fitted; a higher one takes it

   !----- The Kirchhoff lookup table for a general vulnerability exponent: G(r) on a       !
   !      uniform grid (docs/science/plant_hydraulics.md §3). A fixed-size record, so a    !
   !      PFT's parameters copy as plain data. --------------------------------------------!
   integer,  parameter :: HYDRO_TABLE_NTAB = 512        !< intervals over [0, HYDRO_TABLE_RMAX]
   real(wp), parameter :: HYDRO_TABLE_RMAX = 8.0_wp     !< [-] table domain in r = psi/psi50
   real(wp), parameter :: KIRCH_PANEL      = 0.25_wp    !< [-] widest quadrature panel in r
   type :: hydro_table_t
      real(wp)    :: kexp  = 0.0_wp                       !< [-] vulnerability exponent the table holds; its tail uses it too
      real(wp)    :: r_max = HYDRO_TABLE_RMAX             !< [-] table domain [0, r_max]
      real(wp)    :: inv_h = 0.0_wp                       !< [-] 1 / grid spacing
      integer(ik) :: n     = int(HYDRO_TABLE_NTAB, ik)    !< number of intervals
      real(wp)    :: g(0:HYDRO_TABLE_NTAB) = 0.0_wp       !< g(j) = integral_0^{j*h} du/(1+u^kexp)
   end type hydro_table_t

contains

   !----- Sapwood specific conductivity [kg/m/s/MPa] of a wood density rho [g/cm3]: the    !
   !      mean of the lognormal fit, exp(a + b rho + s2/2)                                 !
   !      (docs/science/plant_hydraulics.md eq. (11)). ------------------------------------!
   elemental real(wp) function wood_kmax_from_density(rho) result(kmax)
      real(wp), intent(in) :: rho
      kmax = exp(kmax_ln_a + kmax_ln_b*min(max(rho, kmax_rho_min), kmax_rho_max) + 0.5_wp*kmax_ln_var)
   end function wood_kmax_from_density

   !----- Retained fraction of xylem conductance, 1 - PLC = 1 / (1 + r^kexp),              !
   !      r = psi/psi50, and one at or above zero potential                                !
   !      (docs/science/plant_hydraulics.md eq. (7)). r = 0 is kept off the general power: !
   !      under floating-point trapping a compiler may take 0**kexp through log(0). -------!
   elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (r <= 0.0_wp) then
         f = 1.0_wp
      else
         f = 1.0_wp / (1.0_wp + r**kexp)
      end if
   end function plc_retained

   !----- Kirchhoff (matric flux) potential [MPa], Phi(psi) = integral of the retained     !
   !      fraction from 0 to psi: closed form for kexp 1 and 2, panel quadrature           !
   !      otherwise; zero at and above zero potential (docs/science/plant_hydraulics.md    !
   !      eq. (8)). -----------------------------------------------------------------------!
   pure real(wp) function flux_potential(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = psi50 * atan(r)
      else
         phi = psi50 * kirchhoff_integral(r, kexp)
      end if
   end function flux_potential

   !----- integral_0^r du/(1+u^kexp): 7-point panels no wider than KIRCH_PANEL up to the   !
   !      table's edge, and the analytic tail beyond it. One rule over a long range would  !
   !      miss the integrand's fall near u = 1 and could give a dried leaf's xylem a       !
   !      negative conductance (docs/science/plant_hydraulics.md §3). ---------------------!
   pure real(wp) function kirchhoff_integral(r, kexp) result(g)
      real(wp), intent(in) :: r, kexp
      real(wp)    :: rc, w
      integer(ik) :: k, npan
      rc   = min(max(r, 0.0_wp), HYDRO_TABLE_RMAX)
      npan = max(1_ik, int(ceiling(rc / KIRCH_PANEL), ik))
      w    = rc / real(npan, wp)
      g    = 0.0_wp
      do k = 1_ik, npan
         g = g + kirchhoff_panel(real(k - 1_ik, wp)*w, real(k, wp)*w, kexp)
      end do
      if (r > HYDRO_TABLE_RMAX) g = g + kirchhoff_tail(HYDRO_TABLE_RMAX, r, kexp)
   end function kirchhoff_integral

   !----- integral_a^b du/(1+u^kexp) on one panel (0 <= a <= b) by 7-point Gauss-Legendre, !
   !      exact for polynomials up to degree 13. The sum is written out rather than handed !
   !      a function (docs/building.md, no procedure arguments); u = 0 is guarded as in    !
   !      plc_retained. -------------------------------------------------------------------!
   pure real(wp) function kirchhoff_panel(a, b, kexp) result(g)
      real(wp), intent(in) :: a, b, kexp
      real(wp), parameter :: node(7)   = [ -0.9491079123427585_wp, -0.7415311855993945_wp,     &
                                           -0.4058451513773972_wp,  0.0000000000000000_wp,     &
                                            0.4058451513773972_wp,  0.7415311855993945_wp,     &
                                            0.9491079123427585_wp ]
      real(wp), parameter :: weight(7) = [  0.1294849661688697_wp,  0.2797053914892766_wp,     &
                                            0.3818300505051189_wp,  0.4179591836734694_wp,     &
                                            0.3818300505051189_wp,  0.2797053914892766_wp,     &
                                            0.1294849661688697_wp ]
      real(wp)    :: mid, half, u, acc
      integer(ik) :: k
      mid  = 0.5_wp * (a + b)                    ! the nodes on [-1,1], mapped onto [a,b]
      half = 0.5_wp * (b - a)
      acc  = 0.0_wp
      do k = 1_ik, 7_ik
         u = mid + half * node(k)
         if (u <= 0.0_wp) then
            acc = acc + weight(k)
         else
            acc = acc + weight(k) / (1.0_wp + u**kexp)
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
   !     Fixed-grid lookup table for the general-exponent Kirchhoff integral.              !
   !=======================================================================================!

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

   !----- Kirchhoff potential [MPa] read from the table by linear interpolation; psi50 is  !
   !      a run-time scale, so only kexp sets the table. Beyond the table the analytic     !
   !      tail is added to its last node: extrapolating the last interval could turn a     !
   !      dried leaf's potential positive. ------------------------------------------------!
   pure real(wp) function flux_potential_lin(psi, psi50, tab) result(phi)
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
   end function flux_potential_lin

   !----- The share of a cohort's roots in the soil layer between depths z_top and z_bot,  !
   !      beta^(z_top/D) - beta^(z_bot/D) with D = root_depth and the depths clamped to    !
   !      [0, D], so a layer below the rooting depth gets none. Over [0, D] the shares sum !
   !      to 1 - beta, which the caller divides by (docs/science/plant_hydraulics.md       !
   !      eq. (13)). ----------------------------------------------------------------------!
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

   !=======================================================================================!
   !     The conductances and flows.                                                       !
   !=======================================================================================!

   !----- max(x, 0) with its corner rounded over [0, width] by the cubic 2x^2/w - x^3/w^2, !
   !      which matches both sides in value and slope, so Newton never meets a corner.     !
   !      Every one-way water flow uses it: root uptake, the leaves' supply limit and      !
   !      ground evaporation (docs/science/soil_plant_water.md eq. (1)). ------------------!
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
   !      (docs/science/plant_hydraulics.md eq. (15)). ------------------------------------!
   elemental subroutine root_water_uptake(g_rhizo, k_soil, psi_soil, z_soil, psi_wood, u,             &
                                          du_dpsi_soil, du_dpsi_wood, du_dk)
      real(wp), intent(in)  :: g_rhizo, k_soil, psi_soil, z_soil, psi_wood
      real(wp), intent(out) :: u, du_dpsi_soil, du_dpsi_wood, du_dk
      real(wp) :: drive, dramp
      call smooth_max0(psi_soil + grav_head * z_soil - psi_wood, UPTAKE_RAMP_WIDTH, drive, dramp)
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
   end subroutine xylem_sapflow

   !----- Transpiration [kg/m2/s] of a leaf whose stomata ask for e_demand, at symplastic  !
   !      relative water content r: the demand above r_close, the water content at the     !
   !      closing potential, and below it the demand scaled by r/r_close, so an emptying   !
   !      leaf stops transpiring before it empties. The corner is rounded over             !
   !      TRANSP_CLOSE_WIDTH (docs/science/soil_plant_water.md eq. (3)). ------------------!
   elemental subroutine supply_limited_transpiration(e_demand, r, r_close, e, de_dr)
      real(wp), intent(in)  :: e_demand, r, r_close
      real(wp), intent(out) :: e, de_dr
      real(wp) :: excess, dexcess, share
      call smooth_max0(r - r_close, TRANSP_CLOSE_WIDTH, excess, dexcess)
      share = (r - excess) / r_close
      e     = e_demand * share
      de_dr = e_demand * (1.0_wp - dexcess) / r_close
   end subroutine supply_limited_transpiration

   !=======================================================================================!
   !     The Kirchhoff mean conductance's potential and its slopes.                        !
   !=======================================================================================!

   !----- The Kirchhoff potential [MPa] xylem_sapflow uses: closed form for kexp 1 and 2,  !
   !      the PFT's table otherwise. Above zero potential nothing is embolised, so         !
   !      Phi = psi there, where flux_potential's flat Phi would read a supersaturated     !
   !      tissue's xylem as blocked. ------------------------------------------------------!
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

   !----- d2Phi/dpsi2, the slope of kirchhoff_slope: d/dpsi of 1/(1 + r^kexp) for the      !
   !      closed forms; zero on the table's linear pieces and above zero potential; the    !
   !      tail's slope beyond the table. r = 0 is guarded as in plc_retained. -------------!
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

end module meds_hydr_lib
