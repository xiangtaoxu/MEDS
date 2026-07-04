!==========================================================================================!
! meds_canopy_optics -- derive the two-stream optical coefficients from leaf/wood spectral    !
! properties and the leaf-angle distribution, following SCOPE / 4SAIL (Verhoef 1984).         !
!                                                                                          !
! MU-INDEPENDENT (once, at setup) -- derive_rad_optics fills the per-PFT table:               !
!   omega = rho + tau                              single-scatter albedo                       !
!   g     = bf*(rho - tau)/omega                   leaf-angle scattering asymmetry              !
! from which the diffuse backscatter is  beta  = 0.5*(1 + g)  = sigb/omega.                     !
!                                                                                          !
! MU-DEPENDENT (per timestep, per cohort) -- blend_cohort_optics blends leaf & wood by          !
! clumping-corrected area and adds the solar geometry:                                          !
!   k     = G(mu)/mu                               direct-beam extinction (exact Ross G)         !
!   beta0 = 0.5*(1 + g/k)                          direct-beam upscatter (normalised by k*omega)  !
! plus the leaf/wood absorption split weight leaf_frac = (1-omega_leaf)*elai /                    !
!   [(1-omega_leaf)*elai + (1-omega_wood)*ewai].                                                  !
!                                                                                          !
! For the thermal band the caller passes rho = 1 - emissivity, tau = 0 (ED2/SCOPE convention),   !
! so omega = 1 - emissivity and g = bf; nothing special is needed here.                           !
!==========================================================================================!
module meds_canopy_optics
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num
   use meds_leaf_angle, only : N_LEAF_CLASS, beta_lidf, leaf_bf, gfun_direct
   use meds_rad_types,  only : rad_pft_optics_t, alloc_rad_pft_optics
   implicit none
   private

   public :: derive_rad_optics, blend_cohort_optics

contains

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

end module meds_canopy_optics
