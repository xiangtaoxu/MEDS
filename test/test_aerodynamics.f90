! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_aerodynamics -- unit tests for the stateless canopy-aerodynamics kernel.                !
!   1. NEUTRAL log-law: mo_surface_layer recovers ustar = vonk*u/ln(zldis/z0m) and             !
!      heat_transfer_factor = vonk/ln(zldis/z0m) when the air and canopy air have identical    !
!      (theta, q), and the same                                                                  !
!      ustar a hair either side of neutral.                                                       !
!   2. STABILITY ordering: unstable (warm surface) gives a LARGER ustar than neutral, stable      !
!      (warm air aloft) a SMALLER one, at the same wind. A stable night over a 30 m forest       !
!      keeps zeta within zeta_max_stable, and its coupling rises with the wind.                   !
!   3. reduced_wind is positive and monotonically INCREASING with height (neutral); in the very   !
!      stable range the profile through the surface layer's u* passes through the wind at zref.  !
!   4. boundary_layer_conductance is positive, above the floor, and increases with wind       !
!      (forced convection); an isothermal 5 cm leaf at 20 C and 1 m/s has the laminar          !
!      flat-plate value 0.6 Re^0.5 alpha/l, Re = u l/nu, with the air's 20 C nu and alpha.       !
!   5. MASTER kernel on a 2-cohort patch: positive ustar/ggnet; wind DECREASES top->bottom, by     !
!      the closed-crown extinction; gbw = gbh_2_gbw*gbh; can_depth floored; and ggnet == ggbare   !
!      in the open-canopy limit.                                                                 !
!   6. OPEN CANOPY: the ground conductance is the bare-ground value at an open fraction W = 1    !
!      (a patch with no leaf or wood area), the dense-canopy value at W = 0, and rises            !
!      monotonically between; W = exp(-(LAI + WAI)).                                             !
!   7. TIED cohorts (one height) are one wind layer: a crown split into two cohorts of half the    !
!      leaf area meets the wind the whole crown does, and so does the cohort below.               !
! The equations are docs/science/canopy_aerodynamics.md §2-§4 (eqs. (3), (6)-(10)).        !
!==========================================================================================!
program test_aerodynamics
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : von_karman
   use meds_canopy_aerodynamics, only : aero_env_t, aero_args_t, aero_out_t, alloc_aero_out
   use meds_biophysics_params,   only : aero_params_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics, mo_surface_layer,            &
                                        reduced_wind, boundary_layer_conductance, open_canopy_fraction
   implicit none

   call test_neutral()
   call test_stability_order()
   call test_wind_profile()
   call test_boundary_layer()
   call test_master()
   call test_open_canopy()
   call test_tied_layer()

   call test_report('test_aerodynamics')

contains



   !----- 1. Neutral log-law recovery. -------------------------------------------------------!
   subroutine test_neutral()
      type(aero_params_t) :: cfg
      real(wp) :: ustar, heat_transfer_factor, zeta, rib, obu, uref, z0m, zldis, ln, us_side
      uref = 5.0_wp ; z0m = 0.1_wp ; zldis = 10.0_wp    ! displace = 0 => zldis = zref
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp, 0.010_wp, ustar, heat_transfer_factor, zeta, rib, obu)
      ln = log(zldis / z0m)
      call check('neutral ustar', ustar, von_karman * uref / ln, 1.0e-5_wp)
      call check('neutral heat_transfer_factor', heat_transfer_factor, von_karman / ln,        1.0e-6_wp)
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
   end subroutine test_stability_order

   !----- 3. Wind profile monotonic in height (neutral). -------------------------------------!
   subroutine test_wind_profile()
      type(aero_params_t) :: cfg
      real(wp) :: u_lo, u_hi, ustar, t1, zeta, r, o
      u_lo = reduced_wind(cfg, 0.4_wp, 0.0_wp,  3.0_wp, 0.0_wp, 0.1_wp, 10.0_wp)
      u_hi = reduced_wind(cfg, 0.4_wp, 0.0_wp, 12.0_wp, 0.0_wp, 0.1_wp, 10.0_wp)
      call check_true('wind positive', u_lo > 0.0_wp, u_lo)
      call check_true('wind increases with height', u_hi > u_lo, u_hi - u_lo)
      !----- A 4 K inversion under 3 m/s over a 30 m forest, the stable cap raised to 5: zeta 3.5, !
      !      u* above its floor. ------------------------------------------------------------------!
      cfg%zeta_max_stable = 5.0_wp
      call mo_surface_layer(cfg, 3.0_wp, 40.0_wp, 18.9_wp, 3.9_wp, 290.15_wp, 0.008_wp,         &
                            286.15_wp, 0.008_wp, ustar, t1, zeta, r, o)
      call check_true('very stable case', zeta > 1.0_wp .and. ustar > cfg%ustmin, zeta)
      call check('very stable profile through the wind at zref',                                &
                 reduced_wind(cfg, ustar, zeta, 40.0_wp, 18.9_wp, 3.9_wp, 21.1_wp), 3.0_wp, 1.0e-9_wp)
   end subroutine test_wind_profile

   !----- 4. Boundary-layer conductance positive + increases with wind. ----------------------!
   subroutine test_boundary_layer()
      type(aero_params_t) :: cfg
      real(wp) :: g_lo, g_hi, g_20, nu_20, alpha_20
      g_lo = boundary_layer_conductance(cfg, 0.5_wp, 300.0_wp, 298.0_wp, 0.05_wp, .true.)
      g_hi = boundary_layer_conductance(cfg, 5.0_wp, 300.0_wp, 298.0_wp, 0.05_wp, .true.)
      call check_true('gbh positive above floor', g_lo > cfg%gbhmos_min, g_lo)
      call check_true('gbh increases with wind',  g_hi > g_lo, g_hi - g_lo)
      nu_20    = 1.33e-5_wp * 1.14_wp                    ! [m2/s] air at 20 C
      alpha_20 = 1.89e-5_wp * 1.14_wp
      g_20 = boundary_layer_conductance(cfg, 1.0_wp, 293.15_wp, 293.15_wp, 0.05_wp, .true.)
      call check('gbh laminar flat plate at 20 C', g_20,                                         &
                 0.6_wp * sqrt(1.0_wp * 0.05_wp / nu_20) * alpha_20 / 0.05_wp, 1.0e-12_wp)
   end subroutine test_boundary_layer

   !----- 5. Master kernel on a 2-cohort patch. ----------------------------------------------!
   subroutine test_master()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp) :: height(N), lai(N), tl(N), tw(N), lw(N), bd(N)
      real(wp) :: ggbare_open
      !----- Bottom(1) -> Top(2) ordering (RT convention). ----------------------------------!
      height = [8.0_wp, 18.0_wp] ; lai = [1.5_wp, 2.0_wp]
      tl = [299.0_wp, 300.0_wp] ; tw = [298.5_wp, 299.0_wp]
      lw = [0.04_wp, 0.04_wp]   ; bd = [0.02_wp, 0.02_wp]
      geom%veg_height = 18.0_wp ; geom%opencan_frac = 0.0_wp ; geom%snowfac = 0.0_wp

      call alloc_aero_out(out, N)
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, height, lai, tl, tw, lw, bd, out)

      call check_true('master ustar > 0', out%ustar > 0.0_wp, out%ustar)
      call check_true('master ggnet > 0', out%ggnet > 0.0_wp, out%ggnet)
      call check_true('master ggnet <= ggbare (closed canopy)', out%ggnet <= out%ggbare + 1.0e-12_wp, &
                      out%ggbare - out%ggnet)
      call check_true('wind decreases top->bottom', out%wind(2) > out%wind(1), out%wind(2) - out%wind(1))
      !----- Each closed crown takes the wind down by exp(-LAI/4) to its mid-depth and exp(-LAI/2)  !
      !      to its base (docs/science/canopy_aerodynamics.md eq. (8)). ---------------------------!
      call check('top cohort: mid-crown wind = uh exp(-LAI/4)', out%wind(2),                        &
                 max(cfg%ugbmin, out%uh * exp(-0.25_wp * 2.0_wp)), 1.0e-12_wp)
      call check('bottom cohort: under the top crown, then to its own mid-depth', out%wind(1),       &
                 max(cfg%ugbmin, out%uh * exp(-0.50_wp * 2.0_wp) * exp(-0.25_wp * 1.5_wp)), 1.0e-12_wp)
      call check('leaf gbw = ratio*gbh', out%leaf_gbw(1), cfg%gbh_2_gbw * out%leaf_gbh(1), 1.0e-12_wp)
      call check_true('r_aero = 1/ggnet finite', out%ggnet > 0.0_wp .and. 1.0_wp/out%ggnet < 1.0e6_wp, &
                      1.0_wp / out%ggnet)

      !----- Open-canopy limit: ggnet collapses to ggbare. ----------------------------------!
      geom%opencan_frac = 1.0_wp
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, height, lai, tl, tw, lw, bd, out)
      ggbare_open = out%ggbare
      call check('ggnet == ggbare when open', out%ggnet, ggbare_open, 1.0e-12_wp)
   end subroutine test_master

   !----- 6. The open-canopy blend of the ground conductance (eq. (6)). ----------------------!
   subroutine test_open_canopy()
      type(aero_params_t) :: cfg
      type(aero_env_t)    :: env
      type(aero_args_t)   :: geom
      type(aero_out_t)    :: out
      integer(ik), parameter :: N = 2_ik, NW = 11_ik
      integer(ik), parameter :: ORDER(N) = [1_ik, 2_ik]   ! the arrays are bottom to top already
      real(wp) :: height(N), lai(N), tl(N), tw(N), lw(N), bd(N), gg(NW)
      integer(ik) :: k
      height = [8.0_wp, 18.0_wp] ; lai = [1.5_wp, 2.0_wp]
      tl = [299.0_wp, 300.0_wp] ; tw = [298.5_wp, 299.0_wp]
      lw = [0.04_wp, 0.04_wp]   ; bd = [0.02_wp, 0.02_wp]
      geom%veg_height = 18.0_wp ; geom%snowfac = 0.0_wp
      env%t_ground = env%can_temp + 2.0_wp           ! a warm ground: the dense term's stability is 0
      call alloc_aero_out(out, N)

      !----- The open fraction: 1 on bare ground, exp(-PAI) under plants. --------------------------!
      call check('a patch with no leaf or wood area is fully open', open_canopy_fraction(0.0_wp, 0.0_wp), &
                 1.0_wp, 0.0_wp)
      call check('W = exp(-(LAI + WAI))', open_canopy_fraction(5.0_wp, 1.0_wp), exp(-6.0_wp), 1.0e-15_wp)

      !----- The two limits, exactly, and a monotone rise between them. ----------------------------!
      do k = 1_ik, NW
         geom%opencan_frac = real(k - 1_ik, wp) / real(NW - 1_ik, wp)
         call canopy_aerodynamics(cfg, env, geom, N, ORDER, height, lai, tl, tw, lw, bd, out)
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
      call canopy_aerodynamics(cfg, env, geom, N, ORDER, height, lai, tl, tw, lw, bd, out)
      call check('LAI 6: W ggbare + (1 - W) ggveg', out%ggnet,                                       &
                 geom%opencan_frac * out%ggbare + (1.0_wp - geom%opencan_frac) * out%ggveg, 1.0e-15_wp)
      call check_true('LAI 6: the ground couples nearly as under the dense form alone',                &
                      out%ggnet - out%ggveg < 0.01_wp * (out%ggbare - out%ggveg), out%ggnet / out%ggveg)
   end subroutine test_open_canopy

   !----- 7. Tied cohorts share one wind layer (docs/science/canopy_aerodynamics.md eq. (8)), as they  !
   !      share one radiation layer: test_master's stand with its top crown split into two cohorts of  !
   !      half the leaf area each must give every cohort the wind the unsplit stand gives it. ------!
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
      geom%veg_height = 18.0_wp ; geom%opencan_frac = 0.0_wp ; geom%snowfac = 0.0_wp
      call alloc_aero_out(one, 2_ik) ; call alloc_aero_out(split, 3_ik)
      call canopy_aerodynamics(cfg, env, geom, 2_ik, ORD2, H2, LAI2, TL2, TW2, LW2, BD2, one)
      call canopy_aerodynamics(cfg, env, geom, 3_ik, ORD3, H3, LAI3, TL3, TW3, LW3, BD3, split)
      call check('tied: both halves of the top crown meet its wind', split%wind(2), one%wind(2), 1.0e-14_wp)
      call check('tied: the two halves meet the same wind', split%wind(3), split%wind(2), 0.0_wp)
      call check('tied: the cohort below meets the wind it met under the whole crown', split%wind(1),     &
                 one%wind(1), 1.0e-14_wp)
      call check('tied: the same leaf boundary layer', split%leaf_gbh(3), one%leaf_gbh(2), 1.0e-14_wp)
   end subroutine test_tied_layer

end program test_aerodynamics
