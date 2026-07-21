!==========================================================================================!
! meds_canopy_radiation -- the canopy radiative-transfer module: the per-PFT / per-cohort       !
! optics ASSEMBLY, the unified single-band two-stream SOLVER, and the sealed public SEAM         !
! canopy_radiation (the RT analogue of meds_plant_interface / meds_core_interface). The pure     !
! optical-property kernels (leaf-angle LIDF, single-scatter coefficients) live in the shared     !
! meds_optics_lib; this module `use`s them and adds the RT-specific assembly + solve + boundary. !
!                                                                                          !
!   * derive_rad_optics    -- fill the MU-INDEPENDENT per-PFT optics table (once, at setup).      !
!   * blend_cohort_optics  -- per-cohort leaf/wood optics blend + solar geometry (per timestep).   !
!   * ground_optics        -- ground reflectance + thermal emission (the two-stream lower BC).      !
!   * solve_band / layer_rt -- the O(N) adding-method single-band two-stream solve.                  !
!   * canopy_radiation     -- the seam: loop bands, blend, solve, return per-cohort absorbed rad.    !
!                                                                                          !
! State-free and device-eligible: plain arrays + value types, never a site_t. Cohorts are ordered  !
! BOTTOM (index 1) to TOP (index ncoh). Everything is absolute W m-2 -- no normalize-to-one.         !
!==========================================================================================!
module meds_canopy_radiation
   use meds_kinds,        only : wp, ik
   use meds_constants,    only : stefan, tiny_num, lnexp_max
   use meds_biophysics_types, only : rad_pft_optics_t, rad_forcing_t, rad_flux_t, alloc_rad_flux, &
                                     alloc_rad_pft_optics, surface_state_t
   use meds_optics_lib,   only : N_LEAF_CLASS, beta_lidf, leaf_bf, gfun_direct, scatter_pair
   implicit none
   private

   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t, surface_state_t   ! re-export for callers
   public :: derive_rad_optics, blend_cohort_optics, ground_optics
   public :: solve_band, canopy_radiation

contains

   !=======================================================================================!
   !  CANOPY OPTICS ASSEMBLY -- two-stream scattering coefficients from leaf/wood spectra +  !
   !  the LIDF, following SCOPE / 4SAIL (Verhoef 1984).                                       !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Fill the per-PFT optics table. Spectral inputs are dimensioned (band, pft): for shortwave !
   ! bands pass leaf/wood reflectance & transmittance; for the thermal band pass reflect = 1 -  !
   ! emissivity and trans = 0. `beta_p`,`beta_q` are the per-PFT Beta-LIDF shape parameters.     !
   !---------------------------------------------------------------------------------------!
   subroutine derive_rad_optics(n_band, n_pft, reflect_leaf, trans_leaf, reflect_wood,       &
                                trans_wood, clumping_leaf, clumping_wood, beta_p, beta_q,     &
                                has_beam, has_emission, optics)
      integer(ik), intent(in) :: n_band, n_pft
      real(wp),    intent(in) :: reflect_leaf(n_band, n_pft), trans_leaf(n_band, n_pft)
      real(wp),    intent(in) :: reflect_wood(n_band, n_pft), trans_wood(n_band, n_pft)
      real(wp),    intent(in) :: clumping_leaf(n_pft), clumping_wood(n_pft)
      real(wp),    intent(in) :: beta_p(n_pft), beta_q(n_pft)
      logical,     intent(in) :: has_beam(n_band), has_emission(n_band)
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik) :: b, p

      call alloc_rad_pft_optics(optics, n_band, n_pft, N_LEAF_CLASS)
      optics%clumping_leaf = clumping_leaf
      optics%clumping_wood = clumping_wood
      optics%has_beam      = has_beam
      optics%has_emission  = has_emission

      do p = 1_ik, n_pft
         optics%lidf(:,p) = beta_lidf(beta_p(p), beta_q(p))
         optics%bf(p)     = leaf_bf(optics%lidf(:,p))
         do b = 1_ik, n_band
            call scatter_pair(reflect_leaf(b,p), trans_leaf(b,p), optics%bf(p),               &
                              optics%omega_leaf(b,p), optics%g_leaf(b,p))
            call scatter_pair(reflect_wood(b,p), trans_wood(b,p), optics%bf(p),               &
                              optics%omega_wood(b,p), optics%g_wood(b,p))
         end do
      end do
   end subroutine derive_rad_optics

   !---------------------------------------------------------------------------------------!
   ! Blend leaf & wood optics for every cohort of a patch, for ONE band, and add the solar     !
   ! geometry. Returns per-cohort omega, beta (diffuse backscatter), beta0 (beam upscatter), k  !
   ! (beam extinction) and leaf_frac (leaf share of absorption). `elai`,`ewai` are the clumping- !
   ! corrected leaf and wood area indices (elai = clumping_leaf*lai, ewai = clumping_wood*wai).  !
   !---------------------------------------------------------------------------------------!
   subroutine blend_cohort_optics(optics, band, ncoh, pft, elai, ewai, kdir_in,              &
                                  omega, beta, beta0, kdir, leaf_frac)
      type(rad_pft_optics_t), intent(in)  :: optics
      integer(ik),            intent(in)  :: band, ncoh
      integer(ik),            intent(in)  :: pft(ncoh)
      real(wp),               intent(in)  :: elai(ncoh), ewai(ncoh)
      real(wp),               intent(in)  :: kdir_in(ncoh)   ! band-independent beam extinction, precomputed
      real(wp),               intent(out) :: omega(ncoh), beta(ncoh), beta0(ncoh)
      real(wp),               intent(out) :: kdir(ncoh), leaf_frac(ncoh)
      integer(ik) :: i, ip
      real(wp)    :: etai, wl, ww, gc, aleaf, awood

      do i = 1_ik, ncoh
         ip   = pft(i)
         etai = elai(i) + ewai(i)
         if (etai > tiny_num) then
            wl = elai(i) / etai
            ww = ewai(i) / etai
         else
            wl = 1.0_wp
            ww = 0.0_wp
         end if
         !----- Area-weighted single-scatter albedo and leaf-angle asymmetry. --------------!
         omega(i) = wl * optics%omega_leaf(band,ip) + ww * optics%omega_wood(band,ip)
         gc       = wl * optics%g_leaf(band,ip)     + ww * optics%g_wood(band,ip)
         beta(i)  = 0.5_wp * (1.0_wp + gc)
         !----- Direct-beam upscatter (only meaningful for beam bands); kdir precomputed once. ----!
         if (optics%has_beam(band)) then
            kdir(i)  = kdir_in(i)                              ! band-independent Ross-G extinction
            beta0(i) = 0.5_wp * (1.0_wp + gc / max(kdir(i), tiny_num))
         else
            kdir(i)  = 0.0_wp
            beta0(i) = 0.0_wp
         end if
         !----- Leaf share of absorbed radiation (absorptivity 1-omega, clumping-weighted). -!
         aleaf = (1.0_wp - optics%omega_leaf(band,ip)) * elai(i)
         awood = (1.0_wp - optics%omega_wood(band,ip)) * ewai(i)
         if (aleaf + awood > tiny_num) then
            leaf_frac(i) = aleaf / (aleaf + awood)
         else
            leaf_frac(i) = 1.0_wp
         end if
      end do
   end subroutine blend_cohort_optics

   !=======================================================================================!
   !  SURFACE OPTICS -- ground (soil / litter / water / snow) reflectance and thermal        !
   !  emission, the lower boundary condition of the canopy two-stream. Bare-soil placeholder. !
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

   !=======================================================================================!
   !  TWO-STREAM SOLVER -- the unified single-band adding-method solve over the cohort stack. !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Solve ONE band. Inputs (per cohort, bottom->top): etai (effective total area index),     !
   ! omega/beta/beta0/kdir (blended optics), emission (blackbody sigma*T^4, 0 for shortwave).   !
   ! Scalars: incident beam & diffuse at canopy top, ground reflectance & emission, band flags.  !
   ! Outputs: absorbed(ncoh) [W/m2 per cohort, beam+diffuse], dn_ground/up_ground [W/m2] and     !
   ! albedo (upward fraction leaving the top).                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine solve_band(ncoh, etai, omega, beta, beta0, kdir, emission, incid_beam,         &
                         incid_diff, grnd_refl, grnd_emiss, has_beam, has_emission,           &
                         absorbed, dn_ground, up_ground, albedo)
      integer(ik), intent(in)  :: ncoh
      real(wp),    intent(in)  :: etai(ncoh), omega(ncoh), beta(ncoh), beta0(ncoh), kdir(ncoh)
      real(wp),    intent(in)  :: emission(ncoh)
      real(wp),    intent(in)  :: incid_beam, incid_diff, grnd_refl, grnd_emiss
      logical,     intent(in)  :: has_beam, has_emission
      real(wp),    intent(out) :: absorbed(ncoh)
      real(wp),    intent(out) :: dn_ground, up_ground, albedo

      real(wp), dimension(ncoh)   :: rdd, tdd, sup, sdn          ! layer R/T + up/down sources
      real(wp), dimension(ncoh+1) :: down0                       ! direct-beam profile (interfaces)
      real(wp), dimension(0:ncoh) :: rbel, sbel                  ! reflectance & up-source of stack below
      real(wp), dimension(ncoh)   :: din, dbot, uup, utop        ! diffuse fluxes per layer
      real(wp) :: intercept, scat, emit_face, denom, incid_tot
      integer(ik) :: i

      !----- Direct-beam profile: Beer's law from the top down (cai = 1). -----------------!
      if (has_beam) then
         down0(ncoh+1) = incid_beam
         do i = ncoh, 1_ik, -1_ik
            down0(i) = down0(i+1) * exp(-min(kdir(i) * etai(i), lnexp_max))
         end do
      else
         down0(1:ncoh+1) = 0.0_wp                    ! no beam: zero the whole profile, incl. top face
      end if

      !----- Per-layer diffuse reflectance/transmittance and source fluxes. ----------------!
      do i = 1_ik, ncoh
         call layer_rt(omega(i), beta(i), etai(i), rdd(i), tdd(i))
         intercept = down0(i+1) - down0(i)                        ! beam intercepted by layer i (>=0)
         scat      = omega(i) * intercept                         ! beam single-scattered into diffuse
         emit_face = 0.0_wp
         if (has_emission) emit_face = (1.0_wp - rdd(i) - tdd(i)) * emission(i)
         sup(i) = beta0(i)          * scat + emit_face            ! up-source (top of layer)
         sdn(i) = (1.0_wp - beta0(i)) * scat + emit_face          ! down-source (bottom of layer)
      end do

      !----- Bottom-up sweep: reflectance & upward source of {ground + layers 1..i}. The       !
      !      ground's upward source is its thermal emission PLUS the reflection of the direct    !
      !      beam that reaches it (down0(1)); the diffuse reflection is handled by rbel(0).       !
      rbel(0) = grnd_refl
      sbel(0) = grnd_emiss + grnd_refl * down0(1)
      do i = 1_ik, ncoh
         denom   = 1.0_wp - rbel(i-1) * rdd(i)
         rbel(i) = rdd(i) + tdd(i) * tdd(i) * rbel(i-1) / denom
         sbel(i) = sup(i) + tdd(i) * (sbel(i-1) + rbel(i-1) * sdn(i)) / denom
      end do

      !----- Top-down sweep: recover diffuse fluxes at every interface. --------------------!
      din(ncoh) = incid_diff
      do i = ncoh, 1_ik, -1_ik
         denom   = 1.0_wp - rdd(i) * rbel(i-1)
         dbot(i) = (tdd(i) * din(i) + sdn(i) + rdd(i) * sbel(i-1)) / denom
         uup(i)  = sbel(i-1) + rbel(i-1) * dbot(i)                ! up flux at bottom of layer i
         utop(i) = rdd(i) * din(i) + tdd(i) * uup(i) + sup(i)     ! up flux at top of layer i
         if (i > 1_ik) din(i-1) = dbot(i)
      end do

      !----- Absorbed per layer = beam divergence + diffuse divergence (in - out). ---------!
      do i = 1_ik, ncoh
         absorbed(i) = (down0(i+1) - down0(i)) + (din(i) + uup(i)) - (utop(i) + dbot(i))
      end do

      !----- Ground and top-of-canopy diagnostics. ----------------------------------------!
      dn_ground = dbot(1) + down0(1)                              ! diffuse + beam reaching the ground
      up_ground = uup(1)                                          ! upwelling leaving the ground
      albedo    = 0.0_wp
      incid_tot = incid_beam + incid_diff
      if (incid_tot > tiny_num) albedo = utop(ncoh) / incid_tot
   end subroutine solve_band

   !---------------------------------------------------------------------------------------!
   ! Analytic two-stream diffuse reflectance (rdd) and transmittance (tdd) of a homogeneous    !
   ! layer of effective area index `etai`, single-scatter albedo `omega`, backscatter `beta`.   !
   ! Uses the classical Stokes formulas with rinf = (att-m)/sigb, e = exp(-m*etai).             !
   !---------------------------------------------------------------------------------------!
   pure subroutine layer_rt(omega, beta, etai, rdd, tdd)
      real(wp), intent(in)  :: omega, beta, etai
      real(wp), intent(out) :: rdd, tdd
      real(wp) :: w, sigb, sigf, att, m, rinf, e, e2, denom
      w = min(omega, 1.0_wp - 1.0e-6_wp)          ! guard the conservative-scattering singularity
      sigb = w * beta                              ! backward diffuse scattering
      sigf = w * (1.0_wp - beta)                   ! forward diffuse scattering
      att  = 1.0_wp - sigf
      if (sigb < tiny_num) then                    ! no backscatter -> no diffuse reflectance
         rdd = 0.0_wp
         tdd = exp(-min(att * etai, lnexp_max))
         return
      end if
      m    = sqrt(max(att * att - sigb * sigb, 0.0_wp))
      rinf = (att - m) / sigb
      e    = exp(-min(m * etai, lnexp_max))
      e2   = e * e
      denom = 1.0_wp - rinf * rinf * e2
      tdd  = (1.0_wp - rinf * rinf) * e / denom
      rdd  = rinf * (1.0_wp - e2) / denom
   end subroutine layer_rt

   !=======================================================================================!
   !  THE SEAM -- solve every band for one patch.                                            !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Solve every band for one patch. Cohorts are ordered BOTTOM (1) to TOP (ncoh). `lai`,`wai`  !
   ! are the true (unclumped) leaf & wood area indices; clumping is applied here from the optics. !
   ! `canopy_temp` [K] drives thermal emission (emission bands only). `flux` is allocated here.   !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_radiation(optics, forcing, ncoh, pft, lai, wai, canopy_temp, flux)
      type(rad_pft_optics_t), intent(in)  :: optics
      type(rad_forcing_t),    intent(in)  :: forcing
      integer(ik),            intent(in)  :: ncoh
      integer(ik),            intent(in)  :: pft(ncoh)
      real(wp),               intent(in)  :: lai(ncoh), wai(ncoh), canopy_temp(ncoh)
      type(rad_flux_t),       intent(out) :: flux

      real(wp), dimension(max(ncoh,1)) :: elai, ewai, etai, omega, beta, beta0, kdir, leaf_frac
      real(wp), dimension(max(ncoh,1)) :: emission, absorbed, kdir_coh
      real(wp)    :: dn_ground, up_ground, albedo, incid_tot
      integer(ik) :: b, i, ip

      call alloc_rad_flux(flux, forcing%n_band, ncoh)

      !----- No resolvable cohorts: all radiation reaches the ground. ---------------------!
      if (ncoh < 1_ik) then
         do b = 1_ik, forcing%n_band
            incid_tot = forcing%incid_beam(b) + forcing%incid_diff(b)
            flux%dn_ground(b) = incid_tot
            flux%up_ground(b) = forcing%grnd_emiss(b) + forcing%grnd_refl(b) * incid_tot
            if (incid_tot > tiny_num) flux%albedo(b) = flux%up_ground(b) / incid_tot
         end do
         return
      end if

      !----- Clumping-corrected leaf & wood area indices (per cohort); total etai. ---------!
      do i = 1_ik, ncoh
         ip      = pft(i)
         elai(i) = optics%clumping_leaf(ip) * lai(i)
         ewai(i) = optics%clumping_wood(ip) * wai(i)
         etai(i) = elai(i) + ewai(i)
      end do

      !----- Ross-G beam extinction is band-independent: compute once per cohort, reuse each band. -!
      kdir_coh(1:ncoh) = 0.0_wp
      if (any(optics%has_beam(1:forcing%n_band))) then
         do i = 1_ik, ncoh
            kdir_coh(i) = gfun_direct(optics%lidf(:, pft(i)), forcing%cosz) / max(forcing%cosz, tiny_num)
         end do
      end if

      !----- Loop the configured bands. ---------------------------------------------------!
      do b = 1_ik, forcing%n_band
         call blend_cohort_optics(optics, b, ncoh, pft, elai(1:ncoh), ewai(1:ncoh), kdir_coh(1:ncoh), &
                                  omega(1:ncoh), beta(1:ncoh), beta0(1:ncoh), kdir(1:ncoh),         &
                                  leaf_frac(1:ncoh))

         if (optics%has_emission(b)) then
            emission(1:ncoh) = stefan * canopy_temp(1:ncoh) ** 4     ! blackbody source sigma*T^4
         else
            emission(1:ncoh) = 0.0_wp
         end if

         call solve_band(ncoh, etai(1:ncoh), omega(1:ncoh),                                    &
                         beta(1:ncoh), beta0(1:ncoh), kdir(1:ncoh), emission(1:ncoh),           &
                         forcing%incid_beam(b), forcing%incid_diff(b), forcing%grnd_refl(b),    &
                         forcing%grnd_emiss(b), optics%has_beam(b), optics%has_emission(b),      &
                         absorbed(1:ncoh), dn_ground, up_ground, albedo)

         !----- Split absorbed radiation between leaves and wood (same weighting the solver used). -!
         do i = 1_ik, ncoh
            flux%abs_leaf(b,i) = absorbed(i) * leaf_frac(i)
            flux%abs_wood(b,i) = absorbed(i) * (1.0_wp - leaf_frac(i))
         end do
         flux%albedo(b)    = albedo
         flux%dn_ground(b) = dn_ground
         flux%up_ground(b) = up_ground
      end do
   end subroutine canopy_radiation

end module meds_canopy_radiation
