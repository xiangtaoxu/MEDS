! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_optics_lib -- the optical kernels of the canopy radiation: the leaf-inclination     !
! distribution (LIDF), a Beta distribution on the SCOPE 13-class grid; its moments bf and  !
! G(mu); and the single-scatter pair of a leaf or wood surface (Goel & Strebel 1984;       !
! Verhoef 1984; van der Tol et al. 2009). The equations are                                !
! docs/science/canopy_radiation_transfer.md §1-§2, eqs. (1)-(4); meds_canopy_radiation     !
! builds the per-PFT and per-cohort optics from these kernels.                             !
!==========================================================================================!
module meds_optics_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi, pio2, tiny_num
   implicit none
   private

   public :: n_leaf_class
   public :: beta_lidf, moments_to_beta_params, leaf_bf, gfun_direct
   public :: scatter_pair

   !----- SCOPE's 13 inclination classes: 10 degrees wide to 80, then 2 degrees to 90. ----!
   integer(ik), parameter :: n_leaf_class = 13_ik
   !----- Class midpoints [deg]. ----------------------------------------------------------!
   real(wp), parameter :: LEAF_CLASS_DEG(n_leaf_class) =                                      &
      [  5.0_wp, 15.0_wp, 25.0_wp, 35.0_wp, 45.0_wp, 55.0_wp, 65.0_wp, 75.0_wp,               &
        81.0_wp, 83.0_wp, 85.0_wp, 87.0_wp, 89.0_wp ]
   !----- Class upper bounds [deg]; class 1 starts at 0. ----------------------------------!
   real(wp), parameter :: LEAF_BND_DEG(n_leaf_class) =                                        &
      [ 10.0_wp, 20.0_wp, 30.0_wp, 40.0_wp, 50.0_wp, 60.0_wp, 70.0_wp, 80.0_wp,               &
        82.0_wp, 84.0_wp, 86.0_wp, 88.0_wp, 90.0_wp ]

contains

   !----- The midpoint inclination of class k [rad]. --------------------------------------!
   pure function leaf_class_angle(k) result(theta)
      integer(ik), intent(in) :: k
      real(wp)                :: theta
      theta = LEAF_CLASS_DEG(k) * pi / 180.0_wp
   end function leaf_class_angle

   !----- The Beta LIDF: the leaf-area share of each inclination class, for shape          !
   !      parameters (p, q) (docs/science/canopy_radiation_transfer.md eq. (1)), by a      !
   !      midpoint rule: no special function, and the singular endpoints are never hit. ---!
   pure function beta_lidf(p, q) result(lidf)
      real(wp), intent(in) :: p, q
      real(wp)             :: lidf(n_leaf_class)
      integer(ik), parameter :: NSUB = 64_ik           ! midpoint panels per class
      real(wp)    :: tlo, thi, dt, t, acc, tot
      integer(ik) :: k, j
      tlo = 0.0_wp
      do k = 1_ik, n_leaf_class
         thi = LEAF_BND_DEG(k) / 90.0_wp
         dt  = (thi - tlo) / real(NSUB, wp)
         acc = 0.0_wp
         do j = 1_ik, NSUB
            t   = tlo + (real(j, wp) - 0.5_wp) * dt
            acc = acc + beta_pdf_kernel(t, p, q) * dt
         end do
         lidf(k) = acc
         tlo     = thi
      end do
      tot = sum(lidf)
      if (tot > tiny_num) lidf = lidf / tot             ! the Beta normalizing constant cancels
   end function beta_lidf

   !----- The unnormalized Beta density t^(p-1) (1-t)^(q-1), with t kept off 0 and 1 so an !
   !      integrable endpoint singularity stays finite. -----------------------------------!
   pure function beta_pdf_kernel(t, p, q) result(f)
      real(wp), intent(in) :: t, p, q
      real(wp)             :: f, tt
      tt = min(max(t, 1.0e-9_wp), 1.0_wp - 1.0e-9_wp)   ! clamp-ok: the beta density's domain (0, 1)
      f  = exp((p - 1.0_wp) * log(tt) + (q - 1.0_wp) * log(1.0_wp - tt))
   end function beta_pdf_kernel

   !----- Beta (p, q) from the mean and standard deviation of the leaf angle [deg], by     !
   !      moment matching on t = theta/90 (Goel & Strebel 1984). --------------------------!
   pure subroutine moments_to_beta_params(mean_deg, std_deg, p, q)
      real(wp), intent(in)  :: mean_deg, std_deg
      real(wp), intent(out) :: p, q
      real(wp) :: mt, vt, kappa
      real(wp), parameter :: VT_MIN = 1.0e-6_wp   ! variance floor (std_deg floor ~0.09 deg)
      mt = min(max(mean_deg / 90.0_wp, 1.0e-3_wp), 1.0_wp - 1.0e-3_wp)   ! clamp-ok: the beta mean inside (0, 1)
      vt = (std_deg / 90.0_wp) ** 2
      vt = min(vt, mt * (1.0_wp - mt) * (1.0_wp - 1.0e-6_wp))             ! a Beta needs vt < mt (1 - mt)
      vt = max(vt, VT_MIN)                                                ! keep kappa finite at std_deg = 0
      kappa = mt * (1.0_wp - mt) / vt - 1.0_wp                            ! concentration, p + q
      p = mt * kappa
      q = (1.0_wp - mt) * kappa
   end subroutine moments_to_beta_params

   !----- bf = <cos^2 theta_leaf> over the LIDF (docs/science/canopy_radiation_transfer.md !
   !      eq. (2)). -----------------------------------------------------------------------!
   pure function leaf_bf(lidf) result(bf)
      real(wp), intent(in) :: lidf(:)
      real(wp)             :: bf
      integer(ik)          :: k
      bf = 0.0_wp
      do k = 1_ik, n_leaf_class
         bf = bf + lidf(k) * cos(leaf_class_angle(k)) ** 2
      end do
   end function leaf_bf

   !----- The Ross G-function G(mu), the mean projection of leaf area toward the sun, by   !
   !      SCOPE's volscatt projection (docs/science/canopy_radiation_transfer.md eq. (3)). !
   !      cosz is mu, the cosine of the solar zenith angle. -------------------------------!
   pure function gfun_direct(lidf, cosz) result(gee)
      real(wp), intent(in) :: lidf(:)
      real(wp), intent(in) :: cosz
      real(wp)             :: gee
      real(wp)             :: sinz, cs, ss, as, bts, chi_s
      integer(ik)          :: k
      sinz = sqrt(max(0.0_wp, 1.0_wp - cosz * cosz))
      gee  = 0.0_wp
      do k = 1_ik, n_leaf_class
         cs = cos(leaf_class_angle(k)) * cosz
         ss = sin(leaf_class_angle(k)) * sinz
         as = max(ss, cs)
         if (as < tiny_num) then
            chi_s = 0.0_wp
         else
            bts   = acos(min(max(-cs / as, -1.0_wp), 1.0_wp))
            chi_s = (2.0_wp / pi) * ((bts - pio2) * cs + sin(bts) * ss)
         end if
         gee = gee + lidf(k) * chi_s
      end do
   end function gfun_direct

   !----- Single-scatter albedo omega = rho + tau and asymmetry g = bf (rho - tau)/omega,  !
   !      with g = 0 for a perfect absorber                                                !
   !      (docs/science/canopy_radiation_transfer.md eq. (4)). ----------------------------!
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
