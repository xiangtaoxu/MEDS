! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_canopy_aerodynamics -- the canopy's turbulent exchange (meds_canopy_aerodynamics;   !
! docs/science/canopy_aerodynamics.md §2-§4, eqs. (3), (6)-(10)).                          !
!   1. NEUTRAL (mo_surface_layer): with the air and the canopy air alike, u* and the       !
!      scalar factor are the log law's, kappa u / ln(zldis/z0m) and kappa / ln(zldis/z0m), !
!      and the same a hair either side of neutral.                                         !
!   2. STABILITY (mo_surface_layer): unstable gives a larger u* than neutral and stable a  !
!      smaller one; a stable night over a 30 m forest keeps zeta within zeta_max_stable    !
!      and its coupling rises with the wind; under a weak wind and hot canopy air the      !
!      sweeps settle on the root a 200-sweep solve finds.                                  !
!   3. CANOPY-TOP WIND (canopy_aerodynamics, eq. (7)): neutral, it is the log law's share  !
!      of the wind at zref, ln((h - d)/z0) / ln((zref - d)/z0).                            !
!   4. BOUNDARY LAYER (canopy_aerodynamics, eqs. (9)-(10)): a leaf's conductance is above  !
!      its floor and rises with the wind; an isothermal 5 cm leaf at 20 C has Monteith &   !
!      Unsworth's laminar flat-plate value with the air's tabulated viscosity and          !
!      diffusivity.                                                                        !
!   5. TWO COHORTS (canopy_aerodynamics, eq. (8)): u* and the ground conductance are       !
!      positive, a closed canopy couples the ground less than bare soil, and the wind      !
!      falls top to bottom by each closed crown's extinction.                              !
!   6. OPEN CANOPY (canopy_aerodynamics, eq. (6)): the ground conductance is the bare      !
!      value at W = 1, the dense-canopy value at W = 0, and rises monotonically between;   !
!      W = exp(-(LAI + WAI)); a bare patch couples as CLM5 couples bare soil.              !
!   7. TIED COHORTS (canopy_aerodynamics, eq. (8)): a crown split into two cohorts of half !
!      the leaf area meets the wind the whole crown does, and so does the cohort below.    !
!==========================================================================================!
program test_canopy_aerodynamics
   use meds_test_assert,         only : check, check_true, test_report
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : von_karman
   use meds_biophysics_params,   only : aero_params_t
   use meds_canopy_aerodynamics, only : aero_env_t, aero_args_t, aero_out_t, ensure_aero_out_capacity,          &
                                        canopy_aerodynamics, mo_surface_layer, open_canopy_fraction
   implicit none

   call test_neutral()
   call test_stability_order()
   call test_canopy_top_wind()
   call test_boundary_layer()
   call test_two_cohorts()
   call test_open_canopy()
   call test_tied_layer()

   call test_report('test_canopy_aerodynamics')

contains

   !----- 1. Neutral log-law recovery. -------------------------------------------------------!
   subroutine test_neutral()
      type(aero_params_t) :: cfg
      real(wp) :: ustar, heat_transfer_factor, zeta, rib, obu, uref, z0m, zldis, ln, us_side
      uref  = 5.0_wp
      z0m   = 0.1_wp
      zldis = 10.0_wp                                       ! displace = 0 => zldis = zref
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp, 0.010_wp, ustar, heat_transfer_factor, zeta, rib, obu)
      ln = log(zldis / z0m)
      call check('neutral ustar', ustar, von_karman * uref / ln, 1.0e-5_wp)
      call check('neutral heat_transfer_factor', heat_transfer_factor, von_karman / ln, 1.0e-6_wp)
      call check('neutral zeta ~ 0', zeta, 0.0_wp, 1.0e-6_wp)
      !----- A hair either side of neutral gives the neutral answer: the convective velocity     !
      !      fades with the buoyancy flux, so the wind has no jump there. --------------------------!
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp + 1.0e-9_wp, 0.010_wp, us_side, heat_transfer_factor, zeta, rib, obu)
      call check('ustar a hair unstable of neutral', us_side, ustar, 1.0e-6_wp)
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp - 1.0e-9_wp, 0.010_wp, us_side, heat_transfer_factor, zeta, rib, obu)
      call check('ustar a hair stable of neutral', us_side, ustar, 1.0e-6_wp)
   end subroutine test_neutral

   !----- 2. Stability ordering: unstable > neutral > stable ustar. --------------------------!
   subroutine test_stability_order()
      type(aero_params_t) :: cfg
      real(wp) :: us_u, us_n, us_s, t1, z, r, o, g_1
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 298.15_wp, 0.010_wp,          &
                            298.15_wp, 0.010_wp, us_n, t1, z, r, o)                 ! neutral
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 296.15_wp, 0.010_wp,          &
                            300.15_wp, 0.010_wp, us_u, t1, z, r, o)                 ! CAS warmer -> unstable
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 300.15_wp, 0.010_wp,          &
                            296.15_wp, 0.010_wp, us_s, t1, z, r, o)                 ! air warmer -> stable
      call check_true('unstable ustar > neutral', us_u > us_n, us_u - us_n)
      call check_true('stable ustar < neutral',   us_s < us_n, us_n - us_s)
      !----- Canopy air 3 K colder than the air at zref over a 30 m forest (d 18.9 m, z0 3.9 m). ---!
      call mo_surface_layer(cfg, 1.0_wp, 40.0_wp, 18.9_wp, 3.9_wp, 290.15_wp, 0.008_wp,         &
                            287.15_wp, 0.008_wp, us_s, t1, z, r, o)
      g_1 = us_s * t1
      call check_true('stable sweeps within zeta_max_stable', z <= cfg%zeta_max_stable, z)
      call mo_surface_layer(cfg, 2.0_wp, 40.0_wp, 18.9_wp, 3.9_wp, 290.15_wp, 0.008_wp,         &
                            287.15_wp, 0.008_wp, us_s, t1, z, r, o)
      call check_true('stable coupling rises with wind', us_s * t1 > 1.5_wp * g_1, us_s * t1 - g_1)
      !----- Free convection: 0.5 m/s over a 30 m forest, the canopy air 5 K warmer, where four fixed  !
      !      sweeps missed the coupling by 19 %. --------------------------------------------------------!
      call mo_surface_layer(cfg, 0.5_wp, 35.0_wp, 18.9_wp, 3.9_wp, 298.15_wp, 0.010_wp,         &
                            303.15_wp, 0.010_wp, us_s, t1, z, r, o)
      g_1 = us_s * t1
      cfg%n_iter_mo = 200_ik
      call mo_surface_layer(cfg, 0.5_wp, 35.0_wp, 18.9_wp, 3.9_wp, 298.15_wp, 0.010_wp,         &
                            303.15_wp, 0.010_wp, us_s, t1, z, r, o)
      call check('free convection: the sweeps settle on the root', g_1, us_s * t1, 1.0e-3_wp * us_s * t1)
   end subroutine test_stability_order

   !----- A neutral patch at 20 C: the air and the canopy air alike, and the ground at their      !
   !      temperature, under a wind u_ref at zref = 30 m. ----------------------------------------!
   function neutral_env(u_ref) result(env)
      real(wp), intent(in) :: u_ref   !< [m/s]
      type(aero_env_t)     :: env
      env%u_ref     = u_ref
      env%zref      = 30.0_wp
      env%theta_atm = 293.15_wp
      env%can_temp  = 293.15_wp
      env%t_ground  = 293.15_wp
      env%shv_atm   = 0.010_wp
      env%can_shv   = 0.010_wp
   end function neutral_env

   !----- One cohort of 5 cm leaves at the top of a 20 m canopy, at the canopy air's temperature,   !
   !      under the neutral wind u_ref. ---------------------------------------------------------!
   subroutine one_leaf_cohort(u_ref, out)
      real(wp),         intent(in)    :: u_ref   !< [m/s]
      type(aero_out_t), intent(inout) :: out
      type(aero_params_t) :: cfg
      type(aero_args_t)   :: geom
      integer(ik), parameter :: ORDER(1) = [1_ik]
      geom%veg_height   = 20.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call ensure_aero_out_capacity(out, 1_ik)
      call canopy_aerodynamics(cfg, neutral_env(u_ref), geom, 1_ik, ORDER, [20.0_wp], [0.5_wp], [293.15_wp],  &
                               [293.15_wp], [0.05_wp], [0.02_wp], out)
   end subroutine one_leaf_cohort

   !----- 3. The canopy-top wind (eq. (7)): neutral, u* and the profile are the log law's, so the  !
   !      wind at the canopy top is the wind at zref scaled by ln((h - d)/z0) / ln((zref - d)/z0). !
   subroutine test_canopy_top_wind()
      type(aero_params_t) :: cfg
      type(aero_out_t)    :: out
      real(wp), parameter :: U_REF = 3.0_wp, H = 20.0_wp, ZREF = 30.0_wp
      real(wp) :: z0, d, uh_log
      call one_leaf_cohort(U_REF, out)
      z0 = cfg%z0m_ratio * H
      d  = cfg%d_ratio * H
      uh_log = U_REF * log((H - d) / z0) / log((ZREF - d) / z0)
      call check('neutral canopy-top wind is the log law''s', out%uh, uh_log, 1.0e-12_wp * uh_log)
      call check_true('the canopy-top wind is positive and below the wind at zref',                 &
                      out%uh > 0.0_wp .and. out%uh < U_REF, out%uh)
   end subroutine test_canopy_top_wind

   !----- 4. The leaf boundary layer (eqs. (9)-(10)), through the kernel. An isothermal leaf has no   !
   !      free convection, so at the low Reynolds number of a 5 cm leaf in ~1.5 m/s its conductance is !
   !      the laminar flat plate's, Nu = 0.60 Re^0.5 (Monteith & Unsworth 2013, Principles of        !
   !      Environmental Physics, ch. 10), g = Nu kappa_T / l, with the air's viscosity 1.51e-5 and   !
   !      thermal diffusivity 2.15e-5 m2/s at 20 C (their Appendix A.3). Within 1 %: the linear     !
   !      fits of the two air properties sit within 0.4 % of the table. ---------------------------!
   subroutine test_boundary_layer()
      type(aero_params_t) :: cfg
      type(aero_out_t)    :: lo, hi, mid
      real(wp), parameter :: NU_20 = 1.51e-5_wp, KAPPA_20 = 2.15e-5_wp, WIDTH = 0.05_wp   ! [m2/s], [m2/s], [m]
      real(wp) :: g_plate
      call one_leaf_cohort(1.0_wp, lo)
      call one_leaf_cohort(6.0_wp, hi)
      call check_true('gbh positive above floor', lo%leaf_gbh(1) > cfg%gbhmos_min, lo%leaf_gbh(1))
      call check_true('gbh increases with wind', hi%leaf_gbh(1) > lo%leaf_gbh(1), hi%leaf_gbh(1) - lo%leaf_gbh(1))
      call one_leaf_cohort(3.0_wp, mid)
      g_plate = 0.60_wp * sqrt(mid%wind(1) * WIDTH / NU_20) * KAPPA_20 / WIDTH
      call check('an isothermal 5 cm leaf: the laminar flat plate at 20 C', mid%leaf_gbh(1), g_plate,         &
                 1.0e-2_wp * g_plate)
   end subroutine test_boundary_layer

   !----- 5. Two cohorts under the default forcing. ------------------------------------------!
   subroutine test_two_cohorts()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp), parameter :: HEIGHT(N) = [8.0_wp, 18.0_wp], LAI(N) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL(N) = [299.0_wp, 300.0_wp], TW(N) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW(N) = [0.04_wp, 0.04_wp], BD(N) = [0.02_wp, 0.02_wp]
      geom%veg_height   = 18.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp

      call ensure_aero_out_capacity(out, N)
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)

      call check_true('ustar > 0', out%ustar > 0.0_wp, out%ustar)
      call check_true('ggnet > 0', out%ggnet > 0.0_wp, out%ggnet)
      call check_true('ggnet <= ggbare (closed canopy)', out%ggnet <= out%ggbare + 1.0e-12_wp, out%ggbare - out%ggnet)
      call check_true('wind decreases top->bottom', out%wind(2) > out%wind(1), out%wind(2) - out%wind(1))
      !----- Each closed crown takes the wind down by exp(-LAI/4) to its mid-depth and exp(-LAI/2)  !
      !      to its base (docs/science/canopy_aerodynamics.md eq. (8)). ---------------------------!
      call check('top cohort: mid-crown wind = uh exp(-LAI/4)', out%wind(2),                        &
                 max(cfg%ugbmin, out%uh * exp(-0.25_wp * 2.0_wp)), 1.0e-12_wp)
      call check('bottom cohort: under the top crown, then to its own mid-depth', out%wind(1),       &
                 max(cfg%ugbmin, out%uh * exp(-0.50_wp * 2.0_wp) * exp(-0.25_wp * 1.5_wp)), 1.0e-12_wp)
   end subroutine test_two_cohorts

   !----- 6. The open-canopy blend of the ground conductance (eq. (6)). ----------------------!
   subroutine test_open_canopy()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik, NW = 11_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp), parameter :: HEIGHT(N) = [8.0_wp, 18.0_wp], LAI(N) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL(N) = [299.0_wp, 300.0_wp], TW(N) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW(N) = [0.04_wp, 0.04_wp], BD(N) = [0.02_wp, 0.02_wp]
      real(wp) :: gg(NW), r_meds, r_clm, nu
      integer(ik) :: k
      geom%veg_height = 18.0_wp
      geom%snowfac    = 0.0_wp
      env%t_ground = env%can_temp + 2.0_wp           ! a warm ground: the dense term's stability is 0
      call ensure_aero_out_capacity(out, N)

      !----- The open fraction: 1 on bare ground, exp(-PAI) under plants. --------------------------!
      call check('a patch with no leaf or wood area is fully open', open_canopy_fraction(0.0_wp, 0.0_wp), &
                 1.0_wp, 0.0_wp)
      call check('W = exp(-(LAI + WAI))', open_canopy_fraction(5.0_wp, 1.0_wp), exp(-6.0_wp), 1.0e-15_wp)

      !----- The two limits, exactly, and a monotone rise between them. ----------------------------!
      do k = 1_ik, NW
         geom%opencan_frac = real(k - 1_ik, wp) / real(NW - 1_ik, wp)
         call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
         gg(k) = out%ggnet
         if (k == 1_ik) call check('W = 0: the dense-canopy conductance', out%ggnet, out%ggveg, 0.0_wp)
         if (k == NW)   call check('W = 1: the bare-ground conductance', out%ggnet, out%ggbare, 0.0_wp)
      end do
      call check_true('the bare-ground conductance exceeds the dense one', out%ggbare > out%ggveg,  &
                      out%ggbare - out%ggveg)
      call check_true('the ground conductance rises monotonically with W', all(gg(2:) > gg(:NW - 1)), &
                      minval(gg(2:) - gg(:NW - 1)))
      !----- At LAI 6 (W = 0.0025) the blend is the linear weighting of the two terms, and the       !
      !      ground stays far closer to the dense value than to the bare one. ------------------------!
      geom%opencan_frac = open_canopy_fraction(6.0_wp, 0.0_wp)
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
      call check('LAI 6: W ggbare + (1 - W) ggveg', out%ggnet,                                       &
                 geom%opencan_frac * out%ggbare + (1.0_wp - geom%opencan_frac) * out%ggveg, 1.0e-15_wp)
      call check_true('LAI 6: the ground couples nearly as under the dense form alone',                &
                      out%ggnet - out%ggveg < 0.01_wp * (out%ggbare - out%ggveg), out%ggnet / out%ggveg)
      !----- A bare patch at the soil's roughness, neutral: 1/(u* temp1) + 1/ggbare is CLM5's        !
      !      (ln(zldis/z0m) + 0.13 (z0m u*/nu)^0.45) / (kappa u*), the resistance from the soil to zref. !
      geom%veg_height   = cfg%soil_rough / cfg%z0m_ratio
      geom%opencan_frac = 1.0_wp
      env%can_temp  = 298.15_wp
      env%theta_atm = 298.15_wp
      env%t_ground  = 298.15_wp
      env%can_shv   = env%shv_atm
      call canopy_aerodynamics(cfg, env, geom, 0_ik, ORDER, HEIGHT, LAI, TL, TW, LW, BD, out)
      r_meds = 1.0_wp / (out%ustar * out%heat_transfer_factor) + 1.0_wp / out%ggnet
      nu     = cfg%kin_visc0 * (1.0_wp + cfg%dkin_visc * (298.15_wp - cfg%t_ref_air))
      r_clm  = (log((env%zref - out%displace) / out%rough)                                        &
                + 0.13_wp * (out%rough * out%ustar / nu) ** 0.45_wp) / (von_karman * out%ustar)
      call check('a bare patch is coupled as CLM5 couples bare soil', r_meds, r_clm, 1.0e-10_wp * r_clm)
   end subroutine test_open_canopy

   !----- 7. Tied cohorts share one wind layer (docs/science/canopy_aerodynamics.md eq. (8)), as they  !
   !      share one radiation layer: test_two_cohorts' stand with its top crown split into two cohorts !
   !      of half the leaf area each must give every cohort the wind the unsplit stand gives it. -----!
   subroutine test_tied_layer()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: one, split
      real(wp), parameter :: H2(2) = [8.0_wp, 18.0_wp], LAI2(2) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: TL2(2) = [299.0_wp, 300.0_wp], TW2(2) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: LW2(2) = [0.04_wp, 0.04_wp], BD2(2) = [0.02_wp, 0.02_wp]
      real(wp), parameter :: H3(3) = [8.0_wp, 18.0_wp, 18.0_wp], LAI3(3) = [1.5_wp, 1.0_wp, 1.0_wp]
      real(wp), parameter :: TL3(3) = [299.0_wp, 300.0_wp, 300.0_wp], TW3(3) = [298.5_wp, 299.0_wp, 299.0_wp]
      real(wp), parameter :: LW3(3) = [0.04_wp, 0.04_wp, 0.04_wp], BD3(3) = [0.02_wp, 0.02_wp, 0.02_wp]
      integer(ik), parameter :: ORD2(2) = [1_ik, 2_ik], ORD3(3) = [1_ik, 2_ik, 3_ik]   ! bottom to top
      geom%veg_height   = 18.0_wp
      geom%opencan_frac = 0.0_wp
      geom%snowfac      = 0.0_wp
      call ensure_aero_out_capacity(one, 2_ik)
      call ensure_aero_out_capacity(split, 3_ik)
      call canopy_aerodynamics(cfg, env, geom, 2_ik, ORD2, H2, LAI2, TL2, TW2, LW2, BD2, one)
      call canopy_aerodynamics(cfg, env, geom, 3_ik, ORD3, H3, LAI3, TL3, TW3, LW3, BD3, split)
      call check('tied: both halves of the top crown meet its wind', split%wind(2), one%wind(2), 1.0e-14_wp)
      call check('tied: the two halves meet the same wind', split%wind(3), split%wind(2), 0.0_wp)
      call check('tied: the cohort below meets the wind it met under the whole crown', split%wind(1),     &
                 one%wind(1), 1.0e-14_wp)
      call check('tied: the same leaf boundary layer', split%leaf_gbh(3), one%leaf_gbh(2), 1.0e-14_wp)
   end subroutine test_tied_layer

end program test_canopy_aerodynamics
