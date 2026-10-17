! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_aerodynamics -- unit tests for the stateless canopy-aerodynamics kernel.                !
!   1. NEUTRAL log-law: mo_surface_layer recovers ustar = vonk*u/ln(zldis/z0m) and temp1 =       !
!      vonk/ln(zldis/z0m) when the air and canopy air have identical (theta, q).                  !
!   2. STABILITY ordering: unstable (warm surface) gives a LARGER ustar than neutral, stable      !
!      (warm air aloft) a SMALLER one, at the same wind.                                          !
!   3. reduced_wind is positive and monotonically INCREASING with height (neutral).               !
!   4. boundary_gbh_mos is positive, above the floor, and increases with wind (forced convection).!
!   5. MASTER kernel on a 2-cohort patch: positive ustar/ggnet; wind DECREASES top->bottom, by     !
!      the closed-crown extinction; gbw = gbh_2_gbw*gbh; can_depth floored; and ggnet == ggbare   !
!      in the open-canopy limit.                                                                 !
!   6. TIED cohorts (one height) are one wind layer: a crown split into two cohorts of half the    !
!      leaf area meets the wind the whole crown does, and so does the cohort below.               !
! The equations are docs/science/canopy_aerodynamics.md §2 and §4 (eqs. (3), (7)-(10)).    !
!==========================================================================================!
program test_aerodynamics
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : vonkarman
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_biophysics_opts, only : aero_cfg_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics, mo_surface_layer,            &
                                        reduced_wind, boundary_gbh_mos
   implicit none

   call test_neutral()
   call test_stability_order()
   call test_wind_profile()
   call test_boundary_layer()
   call test_master()
   call test_tied_layer()

   call test_report('test_aerodynamics')

contains



   !----- 1. Neutral log-law recovery. -------------------------------------------------------!
   subroutine test_neutral()
      type(aero_cfg_t) :: cfg
      real(wp) :: ustar, temp1, zeta, rib, obu, uref, z0m, zldis, ln
      uref = 5.0_wp ; z0m = 0.1_wp ; zldis = 10.0_wp    ! displace = 0 => zldis = zref
      call mo_surface_layer(cfg, uref, zldis, 0.0_wp, z0m, 298.15_wp, 0.010_wp,                 &
                            298.15_wp, 0.010_wp, ustar, temp1, zeta, rib, obu)
      ln = log(zldis / z0m)
      call check('neutral ustar', ustar, vonkarman * uref / ln, 1.0e-5_wp)
      call check('neutral temp1', temp1, vonkarman / ln,        1.0e-6_wp)
      call check('neutral zeta ~ 0', zeta, 0.0_wp, 1.0e-6_wp)
   end subroutine test_neutral

   !----- 2. Stability ordering: unstable > neutral > stable ustar. --------------------------!
   subroutine test_stability_order()
      type(aero_cfg_t) :: cfg
      real(wp) :: us_u, us_n, us_s, t1, z, r, o
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 298.15_wp, 0.010_wp,          &
                            298.15_wp, 0.010_wp, us_n, t1, z, r, o)                 ! neutral
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 296.15_wp, 0.010_wp,          &
                            300.15_wp, 0.010_wp, us_u, t1, z, r, o)                 ! CAS warmer -> unstable
      call mo_surface_layer(cfg, 3.0_wp, 10.0_wp, 0.0_wp, 0.1_wp, 300.15_wp, 0.010_wp,          &
                            296.15_wp, 0.010_wp, us_s, t1, z, r, o)                 ! air warmer -> stable
      call check_true('unstable ustar > neutral', us_u > us_n, us_u - us_n)
      call check_true('stable ustar < neutral',   us_s < us_n, us_n - us_s)
   end subroutine test_stability_order

   !----- 3. Wind profile monotonic in height (neutral). -------------------------------------!
   subroutine test_wind_profile()
      type(aero_cfg_t) :: cfg
      real(wp) :: u_lo, u_hi
      u_lo = reduced_wind(cfg, 0.4_wp, 0.0_wp,  3.0_wp, 0.0_wp, 0.1_wp, 10.0_wp)
      u_hi = reduced_wind(cfg, 0.4_wp, 0.0_wp, 12.0_wp, 0.0_wp, 0.1_wp, 10.0_wp)
      call check_true('wind positive', u_lo > 0.0_wp, u_lo)
      call check_true('wind increases with height', u_hi > u_lo, u_hi - u_lo)
   end subroutine test_wind_profile

   !----- 4. Boundary-layer conductance positive + increases with wind. ----------------------!
   subroutine test_boundary_layer()
      type(aero_cfg_t) :: cfg
      real(wp) :: g_lo, g_hi
      g_lo = boundary_gbh_mos(cfg, 0.5_wp, 300.0_wp, 298.0_wp, 0.05_wp, .true.)
      g_hi = boundary_gbh_mos(cfg, 5.0_wp, 300.0_wp, 298.0_wp, 0.05_wp, .true.)
      call check_true('gbh positive above floor', g_lo > cfg%gbhmos_min, g_lo)
      call check_true('gbh increases with wind',  g_hi > g_lo, g_hi - g_lo)
   end subroutine test_boundary_layer

   !----- 5. Master kernel on a 2-cohort patch. ----------------------------------------------!
   subroutine test_master()
      type(aero_cfg_t)  :: cfg
      type(aero_env_t)  :: env
      type(aero_geom_t) :: geom
      type(aero_out_t)  :: out
      integer(ik), parameter :: n = 2_ik
      real(wp) :: height(n), lai(n), tl(n), tw(n), lw(n), bd(n)
      real(wp) :: ggbare_open
      !----- Bottom(1) -> Top(2) ordering (RT convention). ----------------------------------!
      height = [8.0_wp, 18.0_wp] ; lai = [1.5_wp, 2.0_wp]
      tl = [299.0_wp, 300.0_wp] ; tw = [298.5_wp, 299.0_wp]
      lw = [0.04_wp, 0.04_wp]   ; bd = [0.02_wp, 0.02_wp]
      geom%veg_height = 18.0_wp ; geom%opencan_frac = 0.0_wp ; geom%snowfac = 0.0_wp

      call alloc_aero_out(out, n)
      call canopy_aerodynamics(cfg, env, geom, n, height, lai, tl, tw, lw, bd, out)

      call check_true('master ustar > 0', out%ustar > 0.0_wp, out%ustar)
      call check_true('master ggnet > 0', out%ggnet > 0.0_wp, out%ggnet)
      call check_true('master ggnet <= ggbare (closed canopy)', out%ggnet <= out%ggbare + 1.0e-12_wp, &
                      out%ggbare - out%ggnet)
      !----- CAS depth = tallest cohort + freeboard, floored. The canopy air space extends ABOVE  !
      !      the crowns; it is not the canopy volume. -------------------------------------------!
      call check('master can_depth = h_top + freeboard', out%can_depth,                          &
                 max(cfg%min_canopy_depth, 18.0_wp + cfg%canopy_freeboard), 1.0e-12_wp)
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
      call canopy_aerodynamics(cfg, env, geom, n, height, lai, tl, tw, lw, bd, out)
      ggbare_open = out%ggbare
      call check('ggnet == ggbare when open', out%ggnet, ggbare_open, 1.0e-12_wp)
   end subroutine test_master

   !----- 6. Tied cohorts share one wind layer (docs/science/canopy_aerodynamics.md eq. (8)), as they  !
   !      share one radiation layer: test_master's stand with its top crown split into two cohorts of  !
   !      half the leaf area each must give every cohort the wind the unsplit stand gives it. ------!
   subroutine test_tied_layer()
      type(aero_cfg_t)  :: cfg
      type(aero_env_t)  :: env
      type(aero_geom_t) :: geom
      type(aero_out_t)  :: one, split
      real(wp), parameter :: h2(2) = [8.0_wp, 18.0_wp], lai2(2) = [1.5_wp, 2.0_wp]
      real(wp), parameter :: tl2(2) = [299.0_wp, 300.0_wp], tw2(2) = [298.5_wp, 299.0_wp]
      real(wp), parameter :: lw2(2) = [0.04_wp, 0.04_wp], bd2(2) = [0.02_wp, 0.02_wp]
      real(wp), parameter :: h3(3) = [8.0_wp, 18.0_wp, 18.0_wp], lai3(3) = [1.5_wp, 1.0_wp, 1.0_wp]
      real(wp), parameter :: tl3(3) = [299.0_wp, 300.0_wp, 300.0_wp], tw3(3) = [298.5_wp, 299.0_wp, 299.0_wp]
      real(wp), parameter :: lw3(3) = [0.04_wp, 0.04_wp, 0.04_wp], bd3(3) = [0.02_wp, 0.02_wp, 0.02_wp]
      geom%veg_height = 18.0_wp ; geom%opencan_frac = 0.0_wp ; geom%snowfac = 0.0_wp
      call alloc_aero_out(one, 2_ik) ; call alloc_aero_out(split, 3_ik)
      call canopy_aerodynamics(cfg, env, geom, 2_ik, h2, lai2, tl2, tw2, lw2, bd2, one)
      call canopy_aerodynamics(cfg, env, geom, 3_ik, h3, lai3, tl3, tw3, lw3, bd3, split)
      call check('tied: both halves of the top crown meet its wind', split%wind(2), one%wind(2), 1.0e-14_wp)
      call check('tied: the two halves meet the same wind', split%wind(3), split%wind(2), 0.0_wp)
      call check('tied: the cohort below meets the wind it met under the whole crown', split%wind(1),     &
                 one%wind(1), 1.0e-14_wp)
      call check('tied: the same leaf boundary layer', split%leaf_gbh(3), one%leaf_gbh(2), 1.0e-14_wp)
   end subroutine test_tied_layer

end program test_aerodynamics
