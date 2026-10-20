! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_hydraulics -- unit tests for the stateless plant-hydraulics library: the curves and  !
! laws the soil-plant water network is built from (its solve is test_spac_solver).               !
!                                                                                          !
!   1. PV CURVE       : turgor-loss point matches Bartlett eqn 1; psi<->RWC round-trips;        !
!                       capacitance C matches a finite-difference dW/dpsi.                        !
!   2. KIRCHHOFF EDGE : the integrated conductance uses <plc> in (plc(leaf), plc(wood)); the      !
!                       |dpsi|->0 limit reduces to the pointwise value; the table and a Simpson    !
!                       reference.                                                                 !
!   3. SAPWOOD        : the intact xylem conductance is the sapwood's.                             !
!   4. ROOTS          : ED2's root profile and rooting depth; Gardner's root conductance.        !
!   5. SEAM           : growth at conserved water lowers psi; the capacity clamp on a shrink.    !
!   6. LEAFLESS       : a leafless cohort's leaf node sits at the wood less the crown head.      !
! The laws are docs/science/plant_hydraulics.md §2-§4 (eqs. (3)-(15)).                     !
!==========================================================================================!
program test_plant_hydraulics
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : grav_head, pi
   use meds_plant_hydr_lib, only : conductance_retained_fraction, kirchhoff_potential, hydro_table_t,            &
                                   build_hydro_table, root_fraction_profile,                                     &
                                   xylem_max_conductance, rhizosphere_conductance, xylem_sapflow,                &
                                   tissue_psi_tlp, tissue_psi_to_rwc, tissue_rwc_to_psi, tissue_psi_to_water,    &
                                   tissue_capacitance, tissue_water_to_psi, tissue_capped_water, water_curve_t
   use meds_plant_types, only : hydro_params_t, hydro_opts_t
   use meds_plant_hydraulics, only : leaf_node_psi
   use meds_plant_hydraulics, only : cohort_root_depth, cohort_root_profile
   implicit none

   call test_pv_curve()
   call test_xylem_kirchhoff_mean()
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

   !----- A pioneer-like default trait set + cohort environment. ---------------------------!
   subroutine defaults(p, o)
      type(hydro_params_t), intent(out) :: p
      type(hydro_opts_t),   intent(out) :: o
      p%leaf_curve = water_curve_t(pi0 = -1.5_wp, elastic_mod = 12.0_wp, apoplast_frac = 0.30_wp, water_sat = 2.0_wp)
      p%wood_curve = water_curve_t(pi0 = -1.0_wp, elastic_mod =  8.0_wp, apoplast_frac = 0.20_wp, water_sat = 1.0_wp)
      p%wood_psi50 = -2.0_wp ; p%wood_kexp = 2.0_wp
      p%wood_kmax = 9.0_wp ; p%vessel_curl = 1.5_wp      ! 3.0e-3 kg/s/MPa over the 20 m path
      p%root_beta = 0.1_wp
      p%root_depth_b1 = 1.1140580_wp
      p%root_depth_b2 = 0.4223014_wp
      p%specific_root_length = 2.0e4_wp
      p%fine_root_radius = 2.5e-4_wp
      o = hydro_opts_t()   ! gravity on
   end subroutine defaults

   !=======================================================================================!
   subroutine test_pv_curve()
      real(wp), parameter :: PI0 = -1.5_wp, EPS = 12.0_wp, AF = 0.3_wp, WS = 2.0_wp, BM = 0.5_wp
      type(water_curve_t), parameter :: CURVE = water_curve_t(PI0, EPS, AF, WS)
      real(wp) :: psi_tlp, r, psi, cnum, cana, dp, w1, w2
      print '(a)', '-- PV curve --'
      !----- Turgor loss point (Bartlett eqn 1). --------------------------------------!
      psi_tlp = tissue_psi_tlp(PI0, EPS)
      call check('psi_tlp = pi0*eps/(pi0+eps)', psi_tlp, PI0*EPS/(PI0+EPS), 1.0e-12_wp)
      !----- psi<->RWC round trip in the turgid region. -------------------------------!
      r   = 0.9_wp
      psi = tissue_rwc_to_psi(r, PI0, EPS)
      call check('rwc round-trip (turgid)', tissue_psi_to_rwc(psi, PI0, EPS), r, 1.0e-9_wp)
      !----- ... and in the flaccid region. -------------------------------------------!
      r   = 0.5_wp*((PI0+EPS)/EPS)
      psi = tissue_rwc_to_psi(r, PI0, EPS)
      call check('rwc round-trip (flaccid)', tissue_psi_to_rwc(psi, PI0, EPS), r, 1.0e-9_wp)
      !----- Full-turgor water = saturated water. -------------------------------------!
      call check('W(psi=0) = water_sat*biomass', tissue_psi_to_water(0.0_wp, CURVE, BM),  &
                 WS*BM, 1.0e-9_wp)
      !----- Capacitance = finite-difference dW/dpsi. ---------------------------------!
      psi  = -0.8_wp ; dp = 1.0e-5_wp
      w1   = tissue_psi_to_water(psi - dp, CURVE, BM)
      w2   = tissue_psi_to_water(psi + dp, CURVE, BM)
      cnum = (w2 - w1)/(2.0_wp*dp)
      cana = tissue_capacitance(psi, CURVE, BM)
      call check('C(psi) ~ dW/dpsi (turgid)', cana, cnum, 1.0e-6_wp)
   end subroutine test_pv_curve

   !=======================================================================================!
   subroutine test_xylem_kirchhoff_mean()
      real(wp), parameter :: PSI50 = -2.0_wp, KEXP = 2.0_wp, KCOND = 3.0e-3_wp
      type(hydro_table_t) :: tab
      real(wp) :: psi_l, psi_w, avg_plc, s, s0, dum1, dum2, pmid
      print '(a)', '-- Kirchhoff mean conductance of the xylem --'
      call build_hydro_table(tab, KEXP)
      psi_l = -1.2_wp ; psi_w = -0.4_wp
      !----- With no gravity head the sapflow over the drop is the Kirchhoff mean conductance. ---!
      call xylem_sapflow(KCOND, psi_w, psi_l, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      avg_plc = s / ((psi_w - psi_l) * KCOND)
      !----- <plc> lies between the two endpoint retained fractions. -------------------!
      call check_true('<plc> in (plc(leaf), plc(wood))',                                       &
           avg_plc > conductance_retained_fraction(psi_l,PSI50,KEXP)                                &
           .and. avg_plc < conductance_retained_fraction(psi_w,PSI50,KEXP))
      !----- <plc> = (Phi(up)-Phi(down))/(up-down). -----------------------------------!
      call check('<plc> = DeltaPhi/Deltapsi', avg_plc,                                          &
           (kirchhoff_potential(psi_w,PSI50,KEXP,tab)-kirchhoff_potential(psi_l,PSI50,KEXP,tab))/(psi_w-psi_l),  &
           1.0e-11_wp)
      !----- |dpsi|->0 limit is the pointwise value at the mean. -----------------------!
      pmid  = -0.9_wp
      call xylem_sapflow(KCOND, pmid, pmid - 1.0e-9_wp, 0.0_wp, PSI50, KEXP, tab, s0, dum1, dum2)
      call check('Kirchhoff dpsi->0 limit', s0 / (1.0e-9_wp * KCOND), conductance_retained_fraction(pmid,PSI50,KEXP), 1.0e-6_wp)
   end subroutine test_xylem_kirchhoff_mean

   !----- integral_0^r du/(1+u^kexp) by composite Simpson on 2000 intervals: the reference. ---!
   pure real(wp) function simpson_reference(r, kexp) result(g)
      real(wp), intent(in) :: r, kexp
      integer(ik), parameter :: NINT = 2000_ik
      real(wp)    :: h, u
      integer(ik) :: k
      h = r / real(NINT, wp)
      g = 1.0_wp + 1.0_wp / (1.0_wp + r**kexp)
      do k = 1_ik, NINT - 1_ik
         u = real(k, wp) * h
         g = g + merge(4.0_wp, 2.0_wp, mod(k, 2_ik) == 1_ik) / (1.0_wp + u**kexp)
      end do
      g = g * h / 3.0_wp
   end function simpson_reference

   !----- The Kirchhoff potential [MPa] from the Simpson reference, for r = psi/psi50 <= 8. -!
   pure real(wp) function kirchhoff_reference(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      phi = psi50 * simpson_reference(psi / psi50, kexp)
   end function kirchhoff_reference

   !=======================================================================================!
   ! The precomputed lookup table (build_hydro_table, read through kirchhoff_potential)            !
   ! reproduces a Simpson reference of the Kirchhoff potential for a general (non-{1,2})           !
   ! exponent, and the table-backed xylem conductance matches the reference's. kexp is stamped for !
   ! the staleness guard, and psi50 is a runtime scale (one table serves any psi50).               !
   !=======================================================================================!
   subroutine test_hydro_table()
      real(wp), parameter :: KEXP = 3.0_wp, PSI50 = -2.0_wp, KCOND = 3.0e-3_wp
      type(hydro_table_t) :: tab
      real(wp) :: psi, err, emax, psi_l, psi_w, s, dum1, dum2
      integer(ik) :: i
      print '(a)', '-- Kirchhoff lookup table --'
      call build_hydro_table(tab, KEXP)
      call check('table stamps kexp (staleness guard)', tab%kexp, KEXP, 0.0_wp)
      !----- The table matches the reference across the operating range. -----------------------!
      emax = 0.0_wp
      do i = 0_ik, 80_ik
         psi = PSI50 * (0.05_wp + real(i,wp)*0.05_wp)          ! r in [0.05, 4.05]
         err = abs(kirchhoff_potential(psi, PSI50, KEXP, tab) - kirchhoff_reference(psi, PSI50, KEXP))
         emax = max(emax, err)
      end do
      call check_true('linear table within 1e-4 MPa of the reference (kexp=3)', emax < 1.0e-4_wp)
      !----- The table-backed xylem conductance matches the reference's Kirchhoff mean. --------!
      psi_l = -1.7_wp ; psi_w = -0.6_wp
      call xylem_sapflow(KCOND, psi_w, psi_l, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      call check('table-backed xylem conductance ~ reference (kexp=3)', s / (psi_w - psi_l),      &
           KCOND * (kirchhoff_reference(psi_w, PSI50, KEXP) - kirchhoff_reference(psi_l, PSI50, KEXP)) &
           / (psi_w - psi_l), 1.0e-6_wp)
      !----- psi50 is a runtime scale: the SAME table serves a different psi50. ----------------!
      call check('table serves a different psi50 (r-normalized)',                                &
           kirchhoff_potential(-3.0_wp, -3.0_wp, KEXP, tab), kirchhoff_reference(-3.0_wp, -3.0_wp, KEXP), &
           1.0e-4_wp)
   end subroutine test_hydro_table

   !=======================================================================================!
   !----- The maximum conductance is the sapwood's, wood_kmax * sap_area / (height * vessel_curl):  !
   !      doubling the conductivity, doubling the sapwood area and halving the path factor give the !
   !      same step, and a faster one than the base. ------------------------------------------!
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
   end subroutine test_sapwood_conductance

   !=======================================================================================!
   ! ED2's root-fraction profile, and how a drying layer stops feeding the roots: its Gardner        !
   ! conductance falls with its own K(theta) (test_spac_solver: no backflow, the solver's          !
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
   ! Biomass fast/slow SEAM: mass, not psi, is conserved across a slow-loop biomass update, so a !
   ! day of wood growth at FIXED water mass reads as a drier wood, which then draws more water   !
   ! from the soil (the network's uptake follows its potential). --------------------------------!
   !=======================================================================================!
   subroutine test_growth_seam()
      type(hydro_params_t) :: p ; type(hydro_opts_t) :: o
      real(wp) :: bwood0, bwood1, psi_wood0, psi_wood1, w0
      print '(a)', '-- Biomass fast/slow seam: mass-conserving growth --'
      call defaults(p, o)
      bwood0    = 7.0_wp
      psi_wood0 = tissue_rwc_to_psi(0.90_wp, p%wood_curve%pi0, p%wood_curve%elastic_mod)
      w0 = tissue_psi_to_water(psi_wood0, p%wood_curve, bwood0)
      bwood1    = bwood0 * 1.05_wp
      psi_wood1 = tissue_water_to_psi(w0, p%wood_curve, bwood1)
      call check_true('growth at conserved mass lowers psi_wood (drier)', psi_wood1 < psi_wood0)
   end subroutine test_growth_seam

   !=======================================================================================!
   ! Cohort roots: the rooting depth follows ED2's height allometry (5 m for a 35 m tree) and   !
   ! stops at the soil column; a cohort's profile sums to one within its own depth and holds     !
   ! ED2's (1 - beta^(d/D)) / (1 - beta) above depth d; and the single-root conductance is        !
   ! Gardner's value, so a layer split in two keeps its total.                                     !
   !=======================================================================================!
   subroutine test_cohort_roots()
      type(hydro_params_t) :: hp
      type(hydro_opts_t)   :: o
      real(wp)    :: edges(13), frac(12), depth, lv, half, g, g_half
      integer(ik) :: k
      print '(a)', '-- Cohort roots --'
      call defaults(hp, o)                             ! ED2 IALLOM 1, beta 0.1
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
      !----- rhizosphere_conductance is per m2 of ground and per m/s of soil conductivity: for one  !
      !      plant on a m2 with soil conductivity K [m/s] it is Gardner's 2 pi (rho K / grav_head) L /  !
      !      ln(r_half / r_root). --------------------------------------------------------------!
      g      = rhizosphere_conductance(1.0_wp, 1000.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp
      g_half = rhizosphere_conductance(1.0_wp,  500.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp
      call check('Gardner: 2 pi (rho K / grav_head) L / ln(r_half / r_root)', g,                    &
                 2.0_wp * pi * (1000.0_wp * 1.0e-9_wp / grav_head) * 1000.0_wp / log(half / 1.5e-4_wp), 1.0e-15_wp)
      call check('a layer split in two keeps its conductance', 2.0_wp * g_half, g, 1.0e-15_wp)
      call check('twice the plants at the same spacing, twice the conductance',                    &
                 rhizosphere_conductance(2.0_wp, 1000.0_wp, half, 1.5e-4_wp) * 1.0e-9_wp, 2.0_wp * g, 1.0e-15_wp)
   end subroutine test_cohort_roots

   !=======================================================================================!
   ! Saturation-CEILING clamp (the guard mass conservation needs): a persisted mass can only       !
   ! exceed the current capacity after biomass SHRINKS (e.g. the phenology dormant-canopy leaf       !
   ! snap-to-bare), never after it grows. tissue_capped_water must be an exact no-op on growth,       !
   ! must cap exactly at the new ceiling on a shrink, and must reduce a fully-bare tissue to exactly  !
   ! 0 (no NaN/negative) -- see meds_plant_hydr_lib.f90.                                               !
   !=======================================================================================!
   subroutine test_seam_capacity_clamp()
      real(wp), parameter :: PI0 = -1.0_wp, EPS = 8.0_wp, AF = 0.20_wp, WS = 1.0_wp
      type(water_curve_t), parameter :: CURVE = water_curve_t(PI0, EPS, AF, WS)
      real(wp) :: biomass0, biomass_shrunk, biomass_grown, w0, ceiling_shrunk, w_capped
      print '(a)', '-- Biomass fast/slow seam: saturation-ceiling clamp --'
      biomass0 = 7.0_wp
      w0 = tissue_psi_to_water(tissue_rwc_to_psi(0.97_wp, PI0, EPS), CURVE, biomass0)   ! near-saturated

      !----- (a) GROWING biomass: the clamp must be a strict no-op. ---------------------------!
      biomass_grown = biomass0 * 1.05_wp
      call check('clamp is a no-op when biomass grows', tissue_capped_water(w0, CURVE, biomass_grown), &
                 w0, 1.0e-12_wp)

      !----- (b) SHRINKING biomass (a snap-to-bare-like 10x drop): the carried-forward mass would   !
      !          now be wildly supersaturated (w0 above the new ceiling) without the clamp. --------!
      biomass_shrunk = biomass0 * 0.10_wp
      ceiling_shrunk = WS * biomass_shrunk
      call check_true('unclamped carried-forward mass would exceed the new ceiling', w0 > ceiling_shrunk)
      w_capped = tissue_capped_water(w0, CURVE, biomass_shrunk)
      call check('clamp caps mass exactly at the new ceiling', w_capped, ceiling_shrunk, 1.0e-12_wp)
      call check_true('released excess is non-negative', (w0 - w_capped) >= 0.0_wp)
      !----- The clamped mass reads back at exactly rwc=1 (sane, bounded), not the pathological     !
      !      super-turgid extrapolation the unclamped mass would imply. ----------------------------!
      call check('clamped mass diagnoses psi at exactly rwc=1',                                      &
                 tissue_water_to_psi(w_capped, CURVE, biomass_shrunk),                  &
                 tissue_rwc_to_psi(1.0_wp, PI0, EPS), 1.0e-9_wp)

      !----- (c) Full snap-to-bare (biomass -> 0): clamps to exactly 0, not NaN/negative. ----------!
      w_capped = tissue_capped_water(w0, CURVE, 0.0_wp)
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
      real(wp), parameter :: BLEAF = 0.5_wp, HEIGHT = 20.0_wp
      print '(a)', '-- Leafless cohort: leaf node at the wood potential less gravity --'
      call defaults(p, o)
      bwood  = 7.0_wp
      w_wood = tissue_psi_to_water(-1.2_wp, p%wood_curve, bwood)
      psi_w  = tissue_water_to_psi(w_wood, p%wood_curve, bwood)
      psi_l  = leaf_node_psi(0.0_wp, 0.0_wp, w_wood, bwood, HEIGHT, p, o)
      call check('leafless: wood psi less gravity', psi_l, psi_w - grav_head*HEIGHT, 1.0e-12_wp)
      o%gravity_on = .false.
      psi_l  = leaf_node_psi(0.0_wp, 0.0_wp, w_wood, bwood, HEIGHT, p, o)
      call check('leafless, gravity off: wood psi', psi_l, psi_w, 1.0e-12_wp)
      o%gravity_on = .true.
      w_leaf = tissue_psi_to_water(-0.8_wp, p%leaf_curve, BLEAF)
      psi_l  = leaf_node_psi(w_leaf, BLEAF, w_wood, bwood, HEIGHT, p, o)
      call check('with leaves: the leaf curve', psi_l, tissue_water_to_psi(w_leaf, p%leaf_curve, BLEAF), 0.0_wp)
      call check('with leaves: -0.8 MPa recovered', psi_l, -0.8_wp, 1.0e-9_wp)
   end subroutine test_leafless_leaf_psi

   !=======================================================================================!
   ! Far beyond the table (a leaf dried to the apoplastic floor reads about pi0 * 1e4 MPa), the     !
   ! Kirchhoff potential approaches its asymptote psi50 * (pi/a)/sin(pi/a), and the wood-to-leaf    !
   ! edge keeps a positive conductance, so the leaf can refill. At its edge the table equals the    !
   ! Simpson reference.                                                                            !
   !=======================================================================================!
   subroutine test_kirchhoff_far_tail()
      real(wp), parameter :: KEXP = 4.0_wp, PSI50 = -3.209_wp, KCOND = 5.0e-3_wp
      type(hydro_table_t) :: tab
      real(wp) :: g_inf, psi_far, s, dum1, dum2
      print '(a)', '-- Kirchhoff potential far beyond the table --'
      call build_hydro_table(tab, KEXP)
      g_inf   = (pi/KEXP) / sin(pi/KEXP)
      psi_far = -25040.0_wp
      call check('table reaches the asymptote', kirchhoff_potential(psi_far, PSI50, KEXP, tab), PSI50*g_inf,   &
                 1.0e-4_wp)
      call check('table at its edge = reference', kirchhoff_potential(8.0_wp*PSI50, PSI50, KEXP, tab),      &
                 kirchhoff_reference(8.0_wp*PSI50, PSI50, KEXP), 1.0e-9_wp)
      call xylem_sapflow(KCOND, -1.33_wp, psi_far, 0.0_wp, PSI50, KEXP, tab, s, dum1, dum2)
      call check_true('desiccated-leaf sapflow is positive', s > 0.0_wp)
      call check('desiccated-leaf sapflow = k (Phi_wood - Phi_inf)', s,                                   &
                 KCOND * (kirchhoff_reference(-1.33_wp, PSI50, KEXP) - PSI50*g_inf), 1.0e-6_wp)
   end subroutine test_kirchhoff_far_tail

end program test_plant_hydraulics
