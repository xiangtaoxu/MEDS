!==========================================================================================!
! meds_optics -- optical properties throughout the ecosystem column: the leaf-angle           !
! distribution, the canopy (leaf + wood) two-stream scattering coefficients, and the ground /  !
! surface reflectance & emission that closes the two-stream at the bottom. Consolidates the     !
! former meds_leaf_angle + meds_canopy_optics + meds_surface_optics.                            !
!                                                                                          !
! Three optical layers of the column live here, top to bottom:                                !
!   * LEAF ANGLE  -- the two-parameter Beta leaf-inclination distribution (SCOPE / 4SAIL; van    !
!                    der Tol et al. 2009, Verhoef 1984) and its integrals bf = <cos^2 theta>      !
!                    (mu-independent) and the exact Ross G(mu) (mu-dependent, per timestep).       !
!   * CANOPY      -- the mu-independent per-PFT scattering table (omega, g) derived from leaf /    !
!                    wood spectra + the LIDF, and the per-cohort leaf/wood blend with the solar     !
!                    geometry (diffuse backscatter beta, beam upscatter beta0, extinction k).        !
!   * SURFACE     -- ground (soil / litter / water / snow) reflectance and thermal emission, the     !
!                    lower boundary of the two-stream (bare-soil placeholder for now).                !
!==========================================================================================!
module meds_optics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : pi, halfpi, tiny_num, stefan
   use meds_biophysics_types, only : rad_pft_optics_t, alloc_rad_pft_optics
   implicit none
   private

   !----- Leaf-angle distribution (SCOPE / 4SAIL). -------------------------------------------!
   public :: N_LEAF_CLASS, leaf_class_angle
   public :: beta_lidf, beta_params_from_mean, leaf_bf, gfun_direct
   !----- Canopy two-stream scattering coefficients. ----------------------------------------!
   public :: derive_rad_optics, blend_cohort_optics
   !----- Ground / surface optics. ----------------------------------------------------------!
   public :: surface_state_t, ground_optics

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

   !---------------------------------------------------------------------------------------!
   ! Surface state. Only the bare-soil fields are used now; the rest are reserved.          !
   !---------------------------------------------------------------------------------------!
   type :: surface_state_t
      integer(ik)           :: n_band = 0_ik
      real(wp), allocatable :: soil_albedo(:)    !< (band) shortwave soil albedo; unused for emission bands
      real(wp)              :: soil_emiss = 0.96_wp   !< thermal emissivity of the ground
      real(wp)              :: soil_temp  = 298.0_wp  !< [K] ground (skin) temperature
      !----- Reserved for the full implementation (not consulted yet). -------------------!
      real(wp)              :: soil_moisture = 0.0_wp !< [m3/m3] top-layer volumetric water
      real(wp)              :: snow_frac     = 0.0_wp !< [--] snow cover fraction
      real(wp)              :: water_frac    = 0.0_wp !< [--] standing-water fraction
   end type surface_state_t

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

   !=======================================================================================!
   !  CANOPY OPTICS -- two-stream scattering coefficients from leaf/wood spectra + the LIDF, !
   !  following SCOPE / 4SAIL (Verhoef 1984).                                                !
   !                                                                                         !
   !  MU-INDEPENDENT (once, at setup) -- derive_rad_optics fills the per-PFT table:            !
   !    omega = rho + tau                              single-scatter albedo                    !
   !    g     = bf*(rho - tau)/omega                   leaf-angle scattering asymmetry           !
   !  from which the diffuse backscatter is  beta  = 0.5*(1 + g)  = sigb/omega.                  !
   !                                                                                          !
   !  MU-DEPENDENT (per timestep, per cohort) -- blend_cohort_optics blends leaf & wood by       !
   !  clumping-corrected area and adds the solar geometry:                                        !
   !    k     = G(mu)/mu                               direct-beam extinction (exact Ross G)       !
   !    beta0 = 0.5*(1 + g/k)                          direct-beam upscatter (normalised by k*omega)!
   !  plus the leaf/wood absorption split weight leaf_frac = (1-omega_leaf)*elai /                  !
   !    [(1-omega_leaf)*elai + (1-omega_wood)*ewai].                                                !
   !                                                                                          !
   !  For the thermal band the caller passes rho = 1 - emissivity, tau = 0 (ED2/SCOPE           !
   !  convention), so omega = 1 - emissivity and g = bf; nothing special is needed here.          !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Fill the per-PFT optics table. Spectral inputs are dimensioned (band, pft): for shortwave !
   ! bands pass leaf/wood reflectance & transmittance; for the thermal band pass reflect = 1 -  !
   ! emissivity and trans = 0. `beta_p`,`beta_q` are the per-PFT Beta-LIDF shape parameters.     !
   !---------------------------------------------------------------------------------------!
   subroutine derive_rad_optics(n_band, n_pft, reflect_leaf, trans_leaf, reflect_wood,       &
                                trans_wood, clumping_leaf, clumping_wood, beta_p, beta_q,     &
                                has_beam, has_emission, opt)
      integer(ik), intent(in) :: n_band, n_pft
      real(wp),    intent(in) :: reflect_leaf(n_band, n_pft), trans_leaf(n_band, n_pft)
      real(wp),    intent(in) :: reflect_wood(n_band, n_pft), trans_wood(n_band, n_pft)
      real(wp),    intent(in) :: clumping_leaf(n_pft), clumping_wood(n_pft)
      real(wp),    intent(in) :: beta_p(n_pft), beta_q(n_pft)
      logical,     intent(in) :: has_beam(n_band), has_emission(n_band)
      type(rad_pft_optics_t), intent(out) :: opt
      integer(ik) :: b, p

      call alloc_rad_pft_optics(opt, n_band, n_pft, N_LEAF_CLASS)
      opt%clumping_leaf = clumping_leaf
      opt%clumping_wood = clumping_wood
      opt%has_beam      = has_beam
      opt%has_emission  = has_emission

      do p = 1_ik, n_pft
         opt%lidf(:,p) = beta_lidf(beta_p(p), beta_q(p))
         opt%bf(p)     = leaf_bf(opt%lidf(:,p))
         do b = 1_ik, n_band
            call scatter_pair(reflect_leaf(b,p), trans_leaf(b,p), opt%bf(p),                  &
                              opt%omega_leaf(b,p), opt%g_leaf(b,p))
            call scatter_pair(reflect_wood(b,p), trans_wood(b,p), opt%bf(p),                  &
                              opt%omega_wood(b,p), opt%g_wood(b,p))
         end do
      end do
   end subroutine derive_rad_optics

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

   !---------------------------------------------------------------------------------------!
   ! Blend leaf & wood optics for every cohort of a patch, for ONE band, and add the solar     !
   ! geometry. Returns per-cohort omega, beta (diffuse backscatter), beta0 (beam upscatter), k  !
   ! (beam extinction) and leaf_frac (leaf share of absorption). `elai`,`ewai` are the clumping- !
   ! corrected leaf and wood area indices (elai = clumping_leaf*lai, ewai = clumping_wood*wai).  !
   !---------------------------------------------------------------------------------------!
   subroutine blend_cohort_optics(opt, band, cosz, ncoh, pft, elai, ewai,                    &
                                  omega, beta, beta0, kdir, leaf_frac)
      type(rad_pft_optics_t), intent(in)  :: opt
      integer(ik),            intent(in)  :: band, ncoh
      real(wp),               intent(in)  :: cosz
      integer(ik),            intent(in)  :: pft(ncoh)
      real(wp),               intent(in)  :: elai(ncoh), ewai(ncoh)
      real(wp),               intent(out) :: omega(ncoh), beta(ncoh), beta0(ncoh)
      real(wp),               intent(out) :: kdir(ncoh), leaf_frac(ncoh)
      integer(ik) :: i, ip
      real(wp)    :: etai, wl, ww, gc, gee, aleaf, awood

      do i = 1_ik, ncoh
         ip   = pft(i)
         etai = elai(i) + ewai(i)
         if (etai > tiny_num) then
            wl = elai(i) / etai ; ww = ewai(i) / etai
         else
            wl = 1.0_wp ; ww = 0.0_wp
         end if
         !----- Area-weighted single-scatter albedo and leaf-angle asymmetry. --------------!
         omega(i) = wl * opt%omega_leaf(band,ip) + ww * opt%omega_wood(band,ip)
         gc       = wl * opt%g_leaf(band,ip)     + ww * opt%g_wood(band,ip)
         beta(i)  = 0.5_wp * (1.0_wp + gc)
         !----- Direct-beam extinction and upscatter (only meaningful for beam bands). ------!
         if (opt%has_beam(band)) then
            gee      = gfun_direct(opt%lidf(:,ip), cosz)
            kdir(i)  = gee / max(cosz, tiny_num)
            beta0(i) = 0.5_wp * (1.0_wp + gc / max(kdir(i), tiny_num))
         else
            kdir(i)  = 0.0_wp
            beta0(i) = 0.0_wp
         end if
         !----- Leaf share of absorbed radiation (absorptivity 1-omega, clumping-weighted). -!
         aleaf = (1.0_wp - opt%omega_leaf(band,ip)) * elai(i)
         awood = (1.0_wp - opt%omega_wood(band,ip)) * ewai(i)
         if (aleaf + awood > tiny_num) then
            leaf_frac(i) = aleaf / (aleaf + awood)
         else
            leaf_frac(i) = 1.0_wp
         end if
      end do
   end subroutine blend_cohort_optics

   !=======================================================================================!
   !  SURFACE OPTICS -- ground (soil / litter / water / snow) reflectance and thermal        !
   !  emission, the lower boundary condition of the canopy two-stream.                       !
   !                                                                                         !
   !  PLACEHOLDER: this first cut models BARE SOIL only -- a configured per-band soil albedo   !
   !  for shortwave bands, and reflectance = 1 - emissivity with a blackbody source for the     !
   !  thermal band. The `surface_state_t` reserves the fields a full implementation will use     !
   !  (soil moisture, litter, standing water, snow), but they are not consulted yet. When soil    !
   !  state exists this is the one place to grow -- and the place to implement ED2's soil-          !
   !  moisture / colour / snow albedo CORRECTLY (ED2's radiate_driver.f90:639 has an                 !
   !  `albedo_damp_nir = albedo_damp_nir` self-assignment bug in its bedrock branch; do not          !
   !  reproduce it).                                                                                  !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Fill the per-band ground reflectance and thermal emission. For a shortwave band the      !
   ! reflectance is the soil albedo and the emission is zero; for a thermal (emission) band    !
   ! the reflectance is 1 - emissivity and the emission is emissivity * sigma * T_soil^4.       !
   !---------------------------------------------------------------------------------------!
   subroutine ground_optics(surf, n_band, has_emission, grnd_refl, grnd_emiss)
      type(surface_state_t), intent(in)  :: surf
      integer(ik),           intent(in)  :: n_band
      logical,               intent(in)  :: has_emission(n_band)
      real(wp),              intent(out) :: grnd_refl(n_band)
      real(wp),              intent(out) :: grnd_emiss(n_band)
      integer(ik) :: b
      do b = 1_ik, n_band
         if (has_emission(b)) then
            grnd_refl(b)  = 1.0_wp - surf%soil_emiss
            grnd_emiss(b) = surf%soil_emiss * stefan * surf%soil_temp ** 4
         else
            grnd_refl(b)  = surf%soil_albedo(b)
            grnd_emiss(b) = 0.0_wp
         end if
      end do
   end subroutine ground_optics

end module meds_optics
