!==========================================================================================!
! meds_canopy_radiation -- THE sealed public seam of the canopy radiative-transfer library    !
! (the RT analogue of meds_plant_interface / meds_core_interface).                       !
!                                                                                          !
! One public routine, canopy_radiation, takes the precomputed per-PFT optics, the per-band     !
! incident forcing + ground boundary, and a patch's cohorts (PFT, leaf & wood area, temperature),!
! and returns the per-cohort absorbed radiation -- split into leaves and wood, per band -- plus    !
! the patch albedo and below-canopy fluxes. It loops the configured bands, blending leaf/wood     !
! optics per cohort and dispatching each band to the unified two-stream solver. The absorbed PAR   !
! (leaf, RAD_VIS) is the field the leaf-physiology module will consume.                             !
!                                                                                          !
! State-free and device-eligible: it takes plain arrays + value types, never a site_t. A thin      !
! future orchestration layer will walk site_t patches and call this per patch (the analogue of      !
! how the stepper drives the demography interface).                                                  !
!==========================================================================================!
module meds_canopy_radiation
   use meds_kinds,        only : wp, ik
   use meds_constants,    only : stefan, tiny_num
   use meds_biophysics_types, only : rad_pft_optics_t, rad_forcing_t, rad_flux_t, alloc_rad_flux
   use meds_optics,           only : blend_cohort_optics, gfun_direct
   use meds_twostream_band,   only : solve_band
   implicit none
   private

   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t     ! re-export for callers
   public :: canopy_radiation

contains

   !---------------------------------------------------------------------------------------!
   ! Solve every band for one patch. Cohorts are ordered BOTTOM (1) to TOP (ncoh). `lai`,`wai`  !
   ! are the true (unclumped) leaf & wood area indices; clumping is applied here from the optics. !
   ! `canopy_temp` [K] drives thermal emission (emission bands only). `flux` is allocated here.   !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_radiation(opt, forcing, ncoh, pft, lai, wai, canopy_temp, flux)
      type(rad_pft_optics_t), intent(in)  :: opt
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
         elai(i) = opt%clumping_leaf(ip) * lai(i)
         ewai(i) = opt%clumping_wood(ip) * wai(i)
         etai(i) = elai(i) + ewai(i)
      end do

      !----- Ross-G beam extinction is band-independent: compute once per cohort, reuse each band. -!
      kdir_coh(1:ncoh) = 0.0_wp
      if (any(opt%has_beam(1:forcing%n_band))) then
         do i = 1_ik, ncoh
            kdir_coh(i) = gfun_direct(opt%lidf(:, pft(i)), forcing%cosz) / max(forcing%cosz, tiny_num)
         end do
      end if

      !----- Loop the configured bands. ---------------------------------------------------!
      do b = 1_ik, forcing%n_band
         call blend_cohort_optics(opt, b, ncoh, pft, elai(1:ncoh), ewai(1:ncoh), kdir_coh(1:ncoh), &
                                  omega(1:ncoh), beta(1:ncoh), beta0(1:ncoh), kdir(1:ncoh),         &
                                  leaf_frac(1:ncoh))

         if (opt%has_emission(b)) then
            emission(1:ncoh) = stefan * canopy_temp(1:ncoh) ** 4     ! blackbody source sigma*T^4
         else
            emission(1:ncoh) = 0.0_wp
         end if

         call solve_band(ncoh, etai(1:ncoh), omega(1:ncoh),                                    &
                         beta(1:ncoh), beta0(1:ncoh), kdir(1:ncoh), emission(1:ncoh),           &
                         forcing%incid_beam(b), forcing%incid_diff(b), forcing%grnd_refl(b),    &
                         forcing%grnd_emiss(b), opt%has_beam(b), opt%has_emission(b),            &
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
