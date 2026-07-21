!==========================================================================================!
! meds_optics_lib -- shared optical-property kernels: the leaf-inclination distribution        !
! (LIDF) and the two-stream single-scatter coefficients. Pure, semi-empirical, scalar/array-in  !
! constitutive functions (SCOPE / 4SAIL; van der Tol et al. 2009, Verhoef 1984, Goel & Strebel  !
! 1984) -- the optical analogue of meds_hydr_lib / meds_therm_lib -- so they live in             !
! shared/functions and can be evaluated at config-load (per-PFT optics-table build) and on the   !
! GPU hot path alike. The per-PFT/per-cohort ASSEMBLY that fills rad_pft_optics_t and the surface !
! optics live in src/biophysics/meds_canopy_radiation; that module `use`s this one.               !
!==========================================================================================!
module meds_optics_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi, halfpi, tiny_num
   implicit none
   private

   public :: N_LEAF_CLASS, leaf_class_angle
   public :: beta_lidf, beta_params_from_mean, leaf_bf, gfun_direct
   public :: scatter_pair

   !----- SCOPE 13-class inclination grid: 10-deg bins to 80 deg, then 2-deg bins to 89 deg. --!
   integer(ik), parameter :: N_LEAF_CLASS = 13_ik
   !----- Class midpoints [deg]: [5:10:75, 81:2:89]. -----------------------------------------!
   real(wp), parameter :: leaf_class_deg(N_LEAF_CLASS) =                                      &
      [  5.0_wp, 15.0_wp, 25.0_wp, 35.0_wp, 45.0_wp, 55.0_wp, 65.0_wp, 75.0_wp,               &
        81.0_wp, 83.0_wp, 85.0_wp, 87.0_wp, 89.0_wp ]
   !----- Class UPPER boundaries [deg] (lower of class 1 is 0); last is 90. -------------------!
   real(wp), parameter :: leaf_bnd_deg(N_LEAF_CLASS) =                                        &
      [ 10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp, 70.0_wp, 80.0_wp,               &
        82.0_wp, 84.0_wp, 86.0_wp, 88.0_wp, 90.0_wp ]

contains

   !=======================================================================================!
   !  LEAF ANGLE -- leaf-inclination distribution (LIDF) and its integrals (SCOPE / 4SAIL). !
   !                                                                                       !
   !  The LIDF is a two-parameter BETA distribution on the normalized inclination          !
   !  t = theta/(pi/2) in [0,1] (Goel & Strebel 1984) -- a single generic family that       !
   !  subsumes the Campbell ellipsoidal / Verhoef archetypes (planophile / erectophile /     !
   !  plagiophile / extremophile / uniform = Beta(1,1)). Discretized onto the SCOPE 13-class  !
   !  inclination grid; class weights are differences of the regularized incomplete beta.      !
   !                                                                                          !
   !  Two quantities feed the two-stream optics below:                                         !
   !    * bf     = <cos^2(theta_leaf)>  -- 2nd moment; MU-INDEPENDENT, drives the back/forward   !
   !               scattering split. bf=1 horizontal, 0 vertical, 1/3 spherical.                 !
   !    * G(mu)  = <chi_s(theta_leaf, mu)> -- the exact Ross G-function (mean leaf-area            !
   !               projection toward the sun) via the SCOPE `volscatt` projection; MU-DEPENDENT.   !
   !               The direct-beam extinction coefficient is k = G(mu)/mu.                          !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Class midpoint inclination in RADIANS (public helper for tests / optics).             !
   !---------------------------------------------------------------------------------------!
   pure function leaf_class_angle(k) result(theta)
      integer(ik), intent(in) :: k
      real(wp)                :: theta
      theta = leaf_class_deg(k) * pi / 180.0_wp
   end function leaf_class_angle

   !---------------------------------------------------------------------------------------!
   ! Beta-distribution LIDF: area fraction in each of the 13 inclination classes, from the  !
   ! two shape parameters (p, q). Each weight is the integral of the (unnormalized) Beta pdf   !
   ! t^(p-1)*(1-t)^(q-1) over the class, on t = theta/(pi/2); the class weights are renormalized !
   ! to sum to 1 (so the Beta normalizing constant cancels). A midpoint rule is used -- it is    !
   ! compiler-robust (no iterative special function) and never evaluates the possibly-singular    !
   ! endpoints t = 0, 1.                                                                          !
   !---------------------------------------------------------------------------------------!
   pure function beta_lidf(p, q) result(lidf)
      real(wp), intent(in) :: p, q
      real(wp)             :: lidf(N_LEAF_CLASS)
      integer(ik), parameter :: nsub = 64_ik           ! midpoint panels per class
      real(wp)    :: tlo, thi, dt, t, acc, tot
      integer(ik) :: k, j
      tlo = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         thi = leaf_bnd_deg(k) / 90.0_wp
         dt  = (thi - tlo) / real(nsub, wp)
         acc = 0.0_wp
         do j = 1_ik, nsub
            t   = tlo + (real(j, wp) - 0.5_wp) * dt     ! panel midpoint
            acc = acc + beta_pdf_kernel(t, p, q) * dt
         end do
         lidf(k) = acc
         tlo     = thi
      end do
      tot = sum(lidf)
      if (tot > tiny_num) lidf = lidf / tot             ! renormalize
   end function beta_lidf

   !----- Unnormalized Beta density t^(p-1)*(1-t)^(q-1), evaluated via exp/log (any real p,q, --!
   !      t clamped off 0 and 1 so integrable endpoint singularities stay finite). -------------!
   pure function beta_pdf_kernel(t, p, q) result(f)
      real(wp), intent(in) :: t, p, q
      real(wp)             :: f, tt
      tt = min(max(t, 1.0e-9_wp), 1.0_wp - 1.0e-9_wp)
      f  = exp((p - 1.0_wp) * log(tt) + (q - 1.0_wp) * log(1.0_wp - tt))
   end function beta_pdf_kernel

   !---------------------------------------------------------------------------------------!
   ! Map an interpretable (mean leaf angle, standard deviation) [deg] to Beta (p, q) on the  !
   ! normalized inclination t = theta/90. Goel & Strebel (1984) moment matching.            !
   !---------------------------------------------------------------------------------------!
   pure subroutine beta_params_from_mean(mean_deg, std_deg, p, q)
      real(wp), intent(in)  :: mean_deg, std_deg
      real(wp), intent(out) :: p, q
      real(wp) :: mt, vt, kappa
      real(wp), parameter :: VT_MIN = 1.0e-6_wp   ! variance floor (std_deg floor ~0.09 deg)
      mt = min(max(mean_deg / 90.0_wp, 1.0e-3_wp), 1.0_wp - 1.0e-3_wp)   ! mean of t in (0,1)
      vt = (std_deg / 90.0_wp) ** 2                                       ! variance of t
      vt = min(vt, mt * (1.0_wp - mt) * (1.0_wp - 1.0e-6_wp))             ! keep < mt(1-mt)
      vt = max(vt, VT_MIN)                                                ! keep kappa finite at std_deg=0
      kappa = mt * (1.0_wp - mt) / vt - 1.0_wp                            ! concentration
      p = mt * kappa
      q = (1.0_wp - mt) * kappa
   end subroutine beta_params_from_mean

   !---------------------------------------------------------------------------------------!
   ! bf = <cos^2(theta_leaf)> over the LIDF (2nd moment of leaf inclination). MU-independent.!
   !---------------------------------------------------------------------------------------!
   pure function leaf_bf(lidf) result(bf)
      real(wp), intent(in) :: lidf(:)
      real(wp)             :: bf
      integer(ik)          :: k
      bf = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         bf = bf + lidf(k) * cos(leaf_class_angle(k)) ** 2
      end do
   end function leaf_bf

   !---------------------------------------------------------------------------------------!
   ! Ross G-function G(mu) = <chi_s(theta_leaf, theta_sun)> over the LIDF, the mean leaf-area !
   ! projection toward the sun. `cosz` is the cosine of the solar zenith (mu, floored > 0).   !
   ! The beam extinction coefficient is k = G(mu)/mu. Uses the SCOPE `volscatt` projection.    !
   !---------------------------------------------------------------------------------------!
   pure function gfun_direct(lidf, cosz) result(gee)
      real(wp), intent(in) :: lidf(:)
      real(wp), intent(in) :: cosz
      real(wp)             :: gee
      real(wp)             :: tts, sinz, cs, ss, as, bts, chi_s
      integer(ik)          :: k
      tts  = acos(min(max(cosz, -1.0_wp), 1.0_wp))       ! solar zenith [rad]
      sinz = sqrt(max(0.0_wp, 1.0_wp - cosz * cosz))
      gee  = 0.0_wp
      do k = 1_ik, N_LEAF_CLASS
         cs = cos(leaf_class_angle(k)) * cosz
         ss = sin(leaf_class_angle(k)) * sinz
         as = max(ss, cs)
         if (as < tiny_num) then
            chi_s = 0.0_wp
         else
            bts   = acos(min(max(-cs / as, -1.0_wp), 1.0_wp))
            chi_s = (2.0_wp / pi) * ((bts - halfpi) * cs + sin(bts) * ss)
         end if
         gee = gee + lidf(k) * chi_s
      end do
   end function gfun_direct

   !----- omega = rho+tau; g = bf*(rho-tau)/omega (0 if omega ~ 0, e.g. a perfect absorber). --!
   pure subroutine scatter_pair(rho, tau, bf, omega, g)
      real(wp), intent(in)  :: rho, tau, bf
      real(wp), intent(out) :: omega, g
      omega = rho + tau
      if (omega > tiny_num) then
         g = bf * (rho - tau) / omega
      else
         g = 0.0_wp
      end if
   end subroutine scatter_pair

end module meds_optics_lib
