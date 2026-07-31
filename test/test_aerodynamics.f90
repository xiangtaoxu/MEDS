!==========================================================================================!
! test_aerodynamics -- unit tests for the stateless canopy-aerodynamics kernel.                !
!   1. NEUTRAL log-law: mo_surface_layer recovers ustar = vonk*u/ln(zldis/z0m) and temp1 =       !
!      vonk/ln(zldis/z0m) when the air and canopy air have identical (theta, q).                  !
!   2. STABILITY ordering: unstable (warm surface) gives a LARGER ustar than neutral, stable      !
!      (warm air aloft) a SMALLER one, at the same wind.                                          !
!   3. reduced_wind is positive and monotonically INCREASING with height (neutral).               !
!   4. boundary_gbh_mos is positive, above the floor, and increases with wind (forced convection).!
!   5. MASTER kernel on a 2-cohort patch: positive ustar/ggnet; wind DECREASES top->bottom;        !
!      gbw = gbh_2_gbw*gbh; can_depth floored; and ggnet == ggbare in the open-canopy limit.       !
!==========================================================================================!
program test_aerodynamics
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : vonkarman
   use meds_biophysics_types,    only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_canopy_aerodynamics, only : canopy_aerodynamics, mo_surface_layer,            &
                                        reduced_wind, boundary_gbh_mos
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_neutral()
   call test_stability_order()
   call test_wind_profile()
   call test_boundary_layer()
   call test_master()

   if (nfail == 0_ik) then
      print '(a)', 'test_aerodynamics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_aerodynamics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) > atol) then
         print '(a,a,a,es15.7,a,es15.7)', '  FAIL ', name, ': got ', got, ' expected ', expect
         nfail = nfail + 1_ik
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (.not. cond) then
         print '(a,a,a,es15.7)', '  FAIL ', name, ': condition false, val = ', val
         nfail = nfail + 1_ik
      end if
   end subroutine check_true

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
      real(wp) :: height(n), lai(n), crown(n), tl(n), tw(n), lw(n), bd(n)
      real(wp) :: ggbare_open
      !----- Bottom(1) -> Top(2) ordering (RT convention). ----------------------------------!
      height = [8.0_wp, 18.0_wp] ; lai = [1.5_wp, 2.0_wp] ; crown = [0.8_wp, 0.9_wp]
      tl = [299.0_wp, 300.0_wp] ; tw = [298.5_wp, 299.0_wp]
      lw = [0.04_wp, 0.04_wp]   ; bd = [0.02_wp, 0.02_wp]
      geom%veg_height = 18.0_wp ; geom%opencan_frac = 0.0_wp ; geom%snowfac = 0.0_wp

      call alloc_aero_out(out, n)
      call canopy_aerodynamics(cfg, env, geom, n, height, lai, crown, tl, tw, lw, bd, out)

      call check_true('master ustar > 0', out%ustar > 0.0_wp, out%ustar)
      call check_true('master ggnet > 0', out%ggnet > 0.0_wp, out%ggnet)
      call check_true('master ggnet <= ggbare (closed canopy)', out%ggnet <= out%ggbare + 1.0e-12_wp, &
                      out%ggbare - out%ggnet)
      !----- CAS depth = tallest cohort + freeboard, floored. The canopy air space extends ABOVE  !
      !      the crowns; it is not the canopy volume. -------------------------------------------!
      call check('master can_depth = h_top + freeboard', out%can_depth,                          &
                 max(cfg%min_canopy_depth, 18.0_wp + cfg%canopy_freeboard), 1.0e-12_wp)
      call check_true('wind decreases top->bottom', out%wind(2) > out%wind(1), out%wind(2) - out%wind(1))
      call check('leaf gbw = ratio*gbh', out%leaf_gbw(1), cfg%gbh_2_gbw * out%leaf_gbh(1), 1.0e-12_wp)
      call check_true('r_aero = 1/ggnet finite', out%ggnet > 0.0_wp .and. 1.0_wp/out%ggnet < 1.0e6_wp, &
                      1.0_wp / out%ggnet)

      !----- Open-canopy limit: ggnet collapses to ggbare. ----------------------------------!
      geom%opencan_frac = 1.0_wp
      call canopy_aerodynamics(cfg, env, geom, n, height, lai, crown, tl, tw, lw, bd, out)
      ggbare_open = out%ggbare
      call check('ggnet == ggbare when open', out%ggnet, ggbare_open, 1.0e-12_wp)
   end subroutine test_master

end program test_aerodynamics
