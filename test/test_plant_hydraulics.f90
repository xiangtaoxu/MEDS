! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_hydraulics -- unit tests for the stateless plant-hydraulics library: the curves and  !
! laws the soil-plant water network is built from (its solve is test_water_network).             !
!                                                                                          !
!   1. PV CURVE       : turgor-loss point matches Bartlett eqn 1; psi<->RWC round-trips;        !
!                       capacitance C matches a finite-difference dW/dpsi.                        !
!   2. KIRCHHOFF EDGE : the integrated conductance uses <plc> in (plc(leaf), plc(wood)); the      !
!                       |dpsi|->0 limit reduces to the pointwise value; table and quadrature.      !
!   3. SAPWOOD        : the intact xylem conductance is the sapwood's; wood_kmax from density.   !
!   4. ROOTS          : ED2's root profile and rooting depth; Gardner's root conductance.        !
!   5. SEAM           : growth at conserved water lowers psi; the capacity clamp on a shrink.    !
!   6. LEAFLESS       : a leafless cohort's leaf node sits at the wood less the crown head.      !
!==========================================================================================!
program test_plant_hydraulics
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : grav_head, pi
   use meds_hydr_lib,      only : plc_retained, flux_potential, kirchhoff_edge, hydro_table_t, build_hydro_table, &
                                  flux_potential_lin, kirchhoff_edge_tab, wood_kmax_from_density
   use meds_water_retention, only : pv_psi_tlp, rwc_from_psi, psi_from_rwc, water_content, capacitance, &
                                    psi_from_water_content, clamp_water_to_capacity, water_curve_t
   use meds_plant_types, only : hydro_params_t, hydro_opts_t
   use meds_plant_hydraulics, only : leaf_node_psi
   use meds_plant_hydraulics, only : root_fraction_profile, cohort_root_depth, cohort_root_profile
   use meds_water_paths, only : xylem_max_conductance, root_conductance_geometry
   implicit none

   call test_pv_curve()
   call test_kirchhoff_edge()
   call test_kirchhoff_quadrature()
   call test_hydro_table()
   call test_sapwood_conductance()
   call test_root_profile()
   call test_growth_seam()
   call test_cohort_roots()
   call test_seam_capacity_clamp()
   call test_leafless_leaf_psi()
   call test_kirchhoff_far_tail()

   call test_report('test_plant_hydraulics')

contains

   !----- Assertion. -----------------------------------------------------------------------!

   !----- A pioneer-like default trait set + cohort environment. ---------------------------!
   subroutine defaults(p, o)
      type(hydro_params_t), intent(out) :: p
      type(hydro_opts_t),   intent(out) :: o
      p%leaf_curve = water_curve_t(pi0 = -1.5_wp, elastic_mod = 12.0_wp, apoplast_frac = 0.30_wp, water_sat = 2.0_wp)
      p%wood_curve = water_curve_t(pi0 = -1.0_wp, elastic_mod =  8.0_wp, apoplast_frac = 0.20_wp, water_sat = 1.0_wp)
      p%wood_psi50 = -2.0_wp ; p%wood_kexp = 2.0_wp
      p%wood_kmax = 9.0_wp ; p%vessel_curl = 1.5_wp      ! 3.0e-3 kg/s/MPa over the 20 m path
      o = hydro_opts_t()   ! gravity on
   end subroutine defaults

   !=======================================================================================!
   subroutine test_pv_curve()
      real(wp), parameter :: pi0 = -1.5_wp, eps = 12.0_wp, af = 0.3_wp, ws = 2.0_wp, bm = 0.5_wp
      type(water_curve_t), parameter :: curve = water_curve_t(pi0, eps, af, ws)
      real(wp) :: psi_tlp, r, psi, cnum, cana, dp, w1, w2
      print '(a)', '-- PV curve --'
      !----- Turgor loss point (Bartlett eqn 1). --------------------------------------!
      psi_tlp = pv_psi_tlp(pi0, eps)
      call check('psi_tlp = pi0*eps/(pi0+eps)', psi_tlp, pi0*eps/(pi0+eps), 1.0e-12_wp)
      !----- psi<->RWC round trip in the turgid region. -------------------------------!
      r   = 0.9_wp
      psi = psi_from_rwc(r, pi0, eps)
      call check('rwc round-trip (turgid)', rwc_from_psi(psi, pi0, eps), r, 1.0e-9_wp)
      !----- ... and in the flaccid region. -------------------------------------------!
      r   = 0.5_wp*((pi0+eps)/eps)
      psi = psi_from_rwc(r, pi0, eps)
      call check('rwc round-trip (flaccid)', rwc_from_psi(psi, pi0, eps), r, 1.0e-9_wp)
      !----- Full-turgor water = saturated water. -------------------------------------!
      call check('W(psi=0) = water_sat*biomass', water_content(0.0_wp, curve, bm),  &
                 ws*bm, 1.0e-9_wp)
      !----- Capacitance = finite-difference dW/dpsi. ---------------------------------!
      psi  = -0.8_wp ; dp = 1.0e-5_wp
      w1   = water_content(psi - dp, curve, bm)
      w2   = water_content(psi + dp, curve, bm)
      cnum = (w2 - w1)/(2.0_wp*dp)
      cana = capacitance(psi, curve, bm)
      call check('C(psi) ~ dW/dpsi (turgid)', cana, cnum, 1.0e-6_wp)
   end subroutine test_pv_curve

   !=======================================================================================!
   subroutine test_kirchhoff_edge()
      real(wp), parameter :: psi50 = -2.0_wp, kexp = 2.0_wp, kcond = 3.0e-3_wp
      real(wp) :: psi_l, psi_w, avg_plc, keff, keff0, pmid
      print '(a)', '-- Kirchhoff edge --'
      psi_l = -1.2_wp ; psi_w = -0.4_wp
      keff  = kirchhoff_edge(psi_w, psi_l, kcond, psi50, kexp)
      avg_plc = keff/kcond
      !----- <plc> lies between the two endpoint retained fractions. -------------------!
      call check_true('<plc> in (plc(leaf), plc(wood))',                                       &
           avg_plc > plc_retained(psi_l,psi50,kexp) .and. avg_plc < plc_retained(psi_w,psi50,kexp))
      !----- <plc> = (Phi(up)-Phi(down))/(up-down). -----------------------------------!
      call check('<plc> = DeltaPhi/Deltapsi', avg_plc,                                          &
           (flux_potential(psi_w,psi50,kexp)-flux_potential(psi_l,psi50,kexp))/(psi_w-psi_l),  &
           1.0e-12_wp)
      !----- |dpsi|->0 limit is the pointwise value at the mean. -----------------------!
      pmid  = -0.9_wp
      keff0 = kirchhoff_edge(pmid, pmid - 1.0e-9_wp, kcond, psi50, kexp)
      call check('Kirchhoff dpsi->0 limit', keff0/kcond, plc_retained(pmid,psi50,kexp), 1.0e-6_wp)
   end subroutine test_kirchhoff_edge

   !=======================================================================================!
   ! flux_potential's general-exponent branch (kexp not 1 or 2, a 7-point quadrature) matches a   !
   ! fine composite-Simpson reference of the same integral, and at kexp = 2 the same quadrature    !
   ! reproduces the closed form atan(r) it would otherwise replace.                                !
   !=======================================================================================!
   subroutine test_kirchhoff_quadrature()
      real(wp), parameter :: psi50 = -2.0_wp
      real(wp) :: psi, r, ref
      print '(a)', '-- Kirchhoff potential quadrature --'
      psi = -1.3_wp ; r = psi / psi50
      ref = psi50 * simpson_reference(r, 3.0_wp)
      call check('quadrature matches the Simpson reference (kexp=3)', flux_potential(psi, psi50, 3.0_wp), &
                 ref, 1.0e-6_wp)
      psi = -0.7_wp ; r = psi / psi50
      call check('quadrature branch at kexp~2 reproduces atan (closed form)',                         &
                 flux_potential(psi, psi50, 2.0_wp + 2.0e-9_wp), psi50 * atan(r), 1.0e-6_wp)
   end subroutine test_kirchhoff_quadrature

   !----- integral_0^r du/(1+u^kexp) by composite Simpson on 2000 intervals: the reference. ---!
   pure real(wp) function simpson_reference(r, kexp) result(g)
      real(wp), intent(in) :: r, kexp
      integer(ik), parameter :: nint = 2000_ik
      real(wp)    :: h, u
      integer(ik) :: k
      h = r / real(nint, wp)
      g = 1.0_wp + 1.0_wp / (1.0_wp + r**kexp)
      do k = 1_ik, nint - 1_ik
         u = real(k, wp) * h
         g = g + merge(4.0_wp, 2.0_wp, mod(k, 2_ik) == 1_ik) / (1.0_wp + u**kexp)
      end do
      g = g * h / 3.0_wp
   end function simpson_reference

   !=======================================================================================!
   ! The precomputed lookup table (build_hydro_table + flux_potential_lin) reproduces the      !
   ! exact 7-pt-quadrature Kirchhoff potential for a general (non-{1,2}) exponent, and the       !
   ! table-backed edge conductance matches the exact edge. kexp is stamped for the staleness      !
   ! guard, and psi50 is a runtime scale (one table serves any psi50).                             !
   !=======================================================================================!
   subroutine test_hydro_table()
      real(wp), parameter :: kexp = 3.0_wp, psi50 = -2.0_wp, kcond = 3.0e-3_wp
      type(hydro_table_t) :: tab
      real(wp) :: psi, err, emax, psi_l, psi_w
      integer(ik) :: i
      print '(a)', '-- Kirchhoff lookup table --'
      call build_hydro_table(tab, kexp)
      call check('table stamps kexp (staleness guard)', tab%kexp, kexp, 0.0_wp)
      !----- flux_potential_lin(psi) matches the exact quadrature across the operating range. --!
      emax = 0.0_wp
      do i = 0_ik, 80_ik
         psi = psi50 * (0.05_wp + real(i,wp)*0.05_wp)          ! r in [0.05, 4.05]
         err = abs(flux_potential_lin(psi, psi50, tab) - flux_potential(psi, psi50, kexp))
         emax = max(emax, err)
      end do
      call check_true('linear table within 1e-4 MPa of quadrature (kexp=3)', emax < 1.0e-4_wp)
      !----- Table-backed edge conductance matches the exact Kirchhoff edge. -------------------!
      psi_l = -1.7_wp ; psi_w = -0.6_wp
      call check('kirchhoff_edge_tab ~ kirchhoff_edge (kexp=3)',                                  &
           kirchhoff_edge_tab(psi_w, psi_l, kcond, psi50, tab),                                   &
           kirchhoff_edge(psi_w, psi_l, kcond, psi50, kexp), 1.0e-6_wp)
      !----- psi50 is a runtime scale: the SAME table serves a different psi50. ----------------!
      call check('table serves a different psi50 (r-normalized)',                                &
           flux_potential_lin(-3.0_wp, -3.0_wp, tab), flux_potential(-3.0_wp, -3.0_wp, kexp),    &
           1.0e-4_wp)
   end subroutine test_hydro_table

   !=======================================================================================!
   !----- The maximum conductance is the sapwood's, wood_kmax * sap_area / (height * vessel_curl):  !
   !      doubling the conductivity, doubling the sapwood area and halving the path factor give the !
   !      same step, and a faster one than the base. wood_kmax follows wood density (Xu et al.    !
   !      2016): 3.27 at 0.6 g/cm3, falling with density, and held at the fit's 0.35-0.95. -------!
   subroutine test_sapwood_conductance()
      real(wp) :: k0
      print '(a)', '-- Sapwood conductance --'
      k0 = xylem_max_conductance(9.0_wp, 0.01_wp, 20.0_wp, 1.5_wp)
      call check('wood_kmax * sapwood area / (height * vessel_curl)', k0, 9.0_wp * 0.01_wp / 30.0_wp, 1.0e-15_wp)
      call check('2 x sapwood area = 2 x wood_kmax', xylem_max_conductance(9.0_wp, 0.02_wp, 20.0_wp, 1.5_wp), &
                 xylem_max_conductance(18.0_wp, 0.01_wp, 20.0_wp, 1.5_wp), 1.0e-15_wp)
      call check('vessel_curl / 2 = 2 x wood_kmax', xylem_max_conductance(9.0_wp, 0.01_wp, 20.0_wp, 0.75_wp), &
                 2.0_wp * k0, 1.0e-15_wp)
      call check_true('a plant with no sapwood keeps a positive conductance',                         &
                      xylem_max_conductance(9.0_wp, 0.0_wp, 20.0_wp, 1.5_wp) > 0.0_wp)
      call check('wood_kmax at 0.6 g/cm3 (Xu et al. 2016)', wood_kmax_from_density(0.6_wp),       &
                 exp(2.348_wp - 2.455_wp*0.6_wp + 0.5_wp*0.6186_wp), 1.0e-12_wp)
      call check_true('denser wood conducts less',                                               &
                      wood_kmax_from_density(0.8_wp) < wood_kmax_from_density(0.4_wp))
      call check('below the fitted range: the 0.35 value', wood_kmax_from_density(0.2_wp),        &
                 wood_kmax_from_density(0.35_wp), 0.0_wp)
      call check('above the fitted range: the 0.95 value', wood_kmax_from_density(1.1_wp),        &
                 wood_kmax_from_density(0.95_wp), 0.0_wp)
   end subroutine test_sapwood_conductance

   !=======================================================================================!
   ! ED2's root-fraction profile, and how a drying layer stops feeding the roots: its Gardner        !
   ! conductance falls with its own K(theta) (test_water_network: no backflow, the network's       !
   ! uptake by layer).                                                                              !
   !=======================================================================================!
   subroutine test_root_profile()
      real(wp) :: f1, f2, f3
      print '(a)', '-- Root profile --'
      f1 = root_fraction_profile(0.96_wp, 2.0_wp, 0.0_wp, 0.5_wp)   ! equal-thickness (0.5 m) layers
      f2 = root_fraction_profile(0.96_wp, 2.0_wp, 0.5_wp, 1.0_wp)
      f3 = root_fraction_profile(0.96_wp, 2.0_wp, 1.0_wp, 1.5_wp)
      call check_true('shallower (equal-thickness) layer holds more roots', f1 > f2 .and. f2 > f3)
      call check('full column [0,D] telescopes to 1-beta',                                       &
           root_fraction_profile(0.96_wp,2.0_wp,0.0_wp,1.0_wp) +                                  &
           root_fraction_profile(0.96_wp,2.0_wp,1.0_wp,2.0_wp), 1.0_wp-0.96_wp, 1.0e-12_wp)
      call check_true('below rooting depth => 0 roots',                                          &
           root_fraction_profile(0.96_wp,2.0_wp,2.5_wp,3.0_wp) == 0.0_wp)
   end subroutine test_root_profile

   !=======================================================================================!
   ! Biomass fast/slow SEAM (MEDS_ED2_RK45_DESIGN.md P3): mass, not psi, is conserved across a     !
   ! slow-loop biomass update, so a day of wood growth at FIXED water mass reads as a drier wood,  !
   ! which then draws more water from the soil (the network's uptake follows its potential). ----!
   !=======================================================================================!
   subroutine test_growth_seam()
      type(hydro_params_t) :: p ; type(hydro_opts_t) :: o
      real(wp) :: bwood0, bwood1, psi_wood0, psi_wood1, w0
      print '(a)', '-- Biomass fast/slow seam (P3): mass-conserving growth --'
      call defaults(p, o)
      bwood0    = 7.0_wp
      psi_wood0 = psi_from_rwc(0.90_wp, p%wood_curve%pi0, p%wood_curve%elastic_mod)
      w0 = water_content(psi_wood0, p%wood_curve, bwood0)
      bwood1    = bwood0 * 1.05_wp
      psi_wood1 = psi_from_water_content(w0, p%wood_curve, bwood1)
      call check_true('growth at conserved mass lowers psi_wood (drier)', psi_wood1 < psi_wood0)
   end subroutine test_growth_seam

   !=======================================================================================!
   ! Cohort roots: the rooting depth follows ED2's height allometry (5 m for a 35 m tree) and   !
   ! stops at the soil column; a cohort's profile sums to one within its own depth and holds     !
   ! ED2's (1 - beta^(d/D)) / (1 - beta) above depth d; and the single-root conductance is        !
   ! Gardner's value, so a layer split in two keeps its total.                                     !
   !=======================================================================================!
   subroutine test_cohort_roots()
      type(hydro_params_t) :: hp                       ! the defaults: ED2 IALLOM 1, beta 0.1
      real(wp)    :: edges(13), frac(12), depth, lv, half, g, g_half
      integer(ik) :: k
      print '(a)', '-- Cohort roots --'
      call check('a 35 m tree roots to 5 m (ED2 IALLOM 1)', cohort_root_depth(35.0_wp, hp, 10.0_wp), 5.0_wp, 1.0e-3_wp)
      call check('the rooting depth stops at the soil column', cohort_root_depth(35.0_wp, hp, 2.0_wp), 2.0_wp, 0.0_wp)
      call check_true('a 3 m sapling roots shallower than 2 m', cohort_root_depth(3.0_wp, hp, 10.0_wp) < 2.0_wp)
      edges = [(0.5_wp * k, k = 0, 12)]                ! 0 to 6 m in 0.5 m layers
      call cohort_root_profile(35.0_wp, hp, edges(1:12), edges(2:13), frac)
      depth = cohort_root_depth(35.0_wp, hp, 6.0_wp)
      call check('the profile sums to one', sum(frac), 1.0_wp, 1.0e-12_wp)
      call check_true('no roots below the rooting depth', all(frac(11:12) == 0.0_wp))
      call check('ED2''s share above 1 m', sum(frac(1:2)),                                         &
                 (1.0_wp - hp%root_beta**(1.0_wp / depth)) / (1.0_wp - hp%root_beta), 1.0e-12_wp)
      !----- 0.01 plants/m2 with 1000 m of root each in a 0.4 m layer; the same density in two 0.2 m !
      !      halves gives each plant 500 m per half, and the same half-spacing. ------------------------!
      lv   = 0.01_wp * 1000.0_wp / 0.4_wp
      half = 1.0_wp / sqrt(pi * lv)
      !----- root_conductance_geometry is per m2 of ground and per m/s of soil conductivity: for one  !
      !      plant on a m2 with soil conductivity K [m/s] it is Gardner's 2 pi (rho K / grav_head) L /  !
      !      ln(r_half / r_root). --------------------------------------------------------------!
      g      = root_conductance_geometry(1.0_wp, 1000.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp
      g_half = root_conductance_geometry(1.0_wp,  500.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp
      call check('Gardner: 2 pi (rho K / grav_head) L / ln(r_half / r_root)', g,                    &
                 2.0_wp * pi * (1000.0_wp * 1.0e-9_wp / grav_head) * 1000.0_wp / log(half / 1.5e-4_wp), 1.0e-15_wp)
      call check('a layer split in two keeps its conductance', 2.0_wp * g_half, g, 1.0e-15_wp)
      call check('twice the plants at the same spacing, twice the conductance',                    &
                 root_conductance_geometry(2.0_wp, 1000.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp, 2.0_wp * g, 1.0e-15_wp)
   end subroutine test_cohort_roots

   !=======================================================================================!
   ! Saturation-CEILING clamp (P3, the guard mass-conservation needs): a persisted mass can only   !
   ! exceed the current capacity after biomass SHRINKS (e.g. the phenology dormant-canopy leaf       !
   ! snap-to-bare), never after it grows. clamp_water_to_capacity must be an exact no-op on growth,   !
   ! must cap exactly at the new ceiling on a shrink, and must reduce a fully-bare tissue to exactly  !
   ! 0 (no NaN/negative) -- see meds_water_retention.f90.                                              !
   !=======================================================================================!
   subroutine test_seam_capacity_clamp()
      real(wp), parameter :: pi0 = -1.0_wp, eps = 8.0_wp, af = 0.20_wp, ws = 1.0_wp
      type(water_curve_t), parameter :: curve = water_curve_t(pi0, eps, af, ws)
      real(wp) :: biomass0, biomass_shrunk, biomass_grown, w0, ceiling_shrunk, w_capped
      print '(a)', '-- Biomass fast/slow seam (P3): saturation-ceiling clamp --'
      biomass0 = 7.0_wp
      w0 = water_content(psi_from_rwc(0.97_wp, pi0, eps), curve, biomass0)   ! near-saturated

      !----- (a) GROWING biomass: the clamp must be a strict no-op. ---------------------------!
      biomass_grown = biomass0 * 1.05_wp
      call check('clamp is a no-op when biomass grows', clamp_water_to_capacity(w0, curve, biomass_grown), &
                 w0, 1.0e-12_wp)

      !----- (b) SHRINKING biomass (a snap-to-bare-like 10x drop): the carried-forward mass would   !
      !          now be wildly supersaturated (w0 above the new ceiling) without the clamp. --------!
      biomass_shrunk = biomass0 * 0.10_wp
      ceiling_shrunk = ws * biomass_shrunk
      call check_true('unclamped carried-forward mass would exceed the new ceiling', w0 > ceiling_shrunk)
      w_capped = clamp_water_to_capacity(w0, curve, biomass_shrunk)
      call check('clamp caps mass exactly at the new ceiling', w_capped, ceiling_shrunk, 1.0e-12_wp)
      call check_true('released excess is non-negative', (w0 - w_capped) >= 0.0_wp)
      !----- The clamped mass reads back at exactly rwc=1 (sane, bounded), not the pathological     !
      !      super-turgid extrapolation the unclamped mass would imply. ----------------------------!
      call check('clamped mass diagnoses psi at exactly rwc=1',                                      &
                 psi_from_water_content(w_capped, curve, biomass_shrunk),                  &
                 psi_from_rwc(1.0_wp, pi0, eps), 1.0e-9_wp)

      !----- (c) Full snap-to-bare (biomass -> 0): clamps to exactly 0, not NaN/negative. ----------!
      w_capped = clamp_water_to_capacity(w0, curve, 0.0_wp)
      call check('fully-bare tissue clamps to exactly 0', w_capped, 0.0_wp, 0.0_wp)
   end subroutine test_seam_capacity_clamp

   !=======================================================================================!
   ! A leafless cohort (leaf carbon 0, an empty leaf store) has no leaf curve to read: its leaf    !
   ! node takes the wood potential less the gravity head to the crown, not the floored reading    !
   ! (about pi0 * 1e4 MPa). A cohort with leaves still reads its own leaf water.                    !
   !=======================================================================================!
   subroutine test_leafless_leaf_psi()
      type(hydro_params_t) :: p ; type(hydro_opts_t) :: o
      real(wp) :: bwood, w_wood, w_leaf, psi_w, psi_l
      real(wp), parameter :: bleaf = 0.5_wp, height = 20.0_wp
      print '(a)', '-- Leafless cohort: leaf node at the wood potential less gravity --'
      call defaults(p, o)
      bwood  = 7.0_wp
      w_wood = water_content(-1.2_wp, p%wood_curve, bwood)
      psi_w  = psi_from_water_content(w_wood, p%wood_curve, bwood)
      psi_l  = leaf_node_psi(0.0_wp, 0.0_wp, w_wood, bwood, height, p, o)
      call check('leafless: wood psi less gravity', psi_l, psi_w - grav_head*height, 1.0e-12_wp)
      o%gravity_on = .false.
      psi_l  = leaf_node_psi(0.0_wp, 0.0_wp, w_wood, bwood, height, p, o)
      call check('leafless, gravity off: wood psi', psi_l, psi_w, 1.0e-12_wp)
      o%gravity_on = .true.
      w_leaf = water_content(-0.8_wp, p%leaf_curve, bleaf)
      psi_l  = leaf_node_psi(w_leaf, bleaf, w_wood, bwood, height, p, o)
      call check('with leaves: the leaf curve', psi_l, psi_from_water_content(w_leaf, p%leaf_curve, bleaf), 0.0_wp)
      call check('with leaves: -0.8 MPa recovered', psi_l, -0.8_wp, 1.0e-9_wp)
   end subroutine test_leafless_leaf_psi

   !=======================================================================================!
   ! Far beyond the table (a leaf dried to the apoplastic floor reads about pi0 * 1e4 MPa), the     !
   ! Kirchhoff potential approaches its asymptote psi50 * (pi/a)/sin(pi/a), and the wood-to-leaf    !
   ! edge keeps a positive conductance, so the leaf can refill. The table and the quadrature agree. !
   !=======================================================================================!
   subroutine test_kirchhoff_far_tail()
      real(wp), parameter :: kexp = 4.0_wp, psi50 = -3.209_wp, kcond = 5.0e-3_wp
      type(hydro_table_t) :: tab
      real(wp) :: g_inf, psi_far
      print '(a)', '-- Kirchhoff potential far beyond the table --'
      call build_hydro_table(tab, kexp)
      g_inf   = (pi/kexp) / sin(pi/kexp)
      psi_far = -25040.0_wp
      call check('quadrature reaches the asymptote', flux_potential(psi_far, psi50, kexp), psi50*g_inf, 1.0e-4_wp)
      call check('table reaches the asymptote', flux_potential_lin(psi_far, psi50, tab), psi50*g_inf, 1.0e-4_wp)
      call check('table at its edge = quadrature', flux_potential_lin(8.0_wp*psi50, psi50, tab),         &
                 flux_potential(8.0_wp*psi50, psi50, kexp), 1.0e-9_wp)
      call check_true('desiccated-leaf edge conductance is positive',                                    &
                      kirchhoff_edge_tab(-1.33_wp, psi_far, kcond, psi50, tab) > 0.0_wp)
      call check('desiccated-leaf flux = k (Phi_wood - Phi_inf)',                                        &
                 kirchhoff_edge_tab(-1.33_wp, psi_far, kcond, psi50, tab) * (-1.33_wp - psi_far),        &
                 kcond * (flux_potential(-1.33_wp, psi50, kexp) - psi50*g_inf), 1.0e-6_wp)
   end subroutine test_kirchhoff_far_tail

end program test_plant_hydraulics
