!==========================================================================================!
! meds_rad_types -- data structures of the canopy radiative-transfer interface.             !
!                                                                                          !
! Pure DATA (no methods, no hidden state), the RT analogue of meds_leaf_types. Three types:  !
!   * rad_pft_optics_t -- the PRECOMPUTED, MU-INDEPENDENT per-PFT optics table (single-scatter !
!                         albedo omega and the leaf-angle scattering asymmetry g per band and   !
!                         tissue, plus the leaf-angle distribution lidf and its 2nd moment bf).   !
!                         Filled ONCE by meds_canopy_optics%derive_rad_optics.                    !
!   * rad_forcing_t    -- the per-band incident fluxes (beam & diffuse), the solar cosine, and    !
!                         the ground reflectance/emission (from meds_surface_optics). Absolute      !
!                         W/m2 throughout -- no ED2-style normalize-to-one.                          !
!   * rad_flux_t       -- the returned per-cohort absorbed radiation (leaf & wood, per band) plus    !
!                         the patch-level albedo and below-canopy transmission diagnostics.           !
!                                                                                          !
! Bands are configurable; the default set is {VIS, NIR, LW}. Every band carries a thermal-      !
! emission term, identically zero for VIS/NIR (has_emission = .false.).                          !
!==========================================================================================!
module meds_rad_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux

   !----- Default three-band layout (indices into the band dimension). --------------------!
   integer(ik), parameter :: RAD_VIS = 1_ik   !< visible / PAR (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_NIR = 2_ik   !< near-infrared (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_LW  = 3_ik   !< thermal / longwave (diffuse + emission, no beam)
   integer(ik), parameter :: N_RAD_BAND_DEFAULT = 3_ik

   !---------------------------------------------------------------------------------------!
   ! Precomputed per-PFT optics. omega = rho + tau (single-scatter albedo); g = bf*(rho-tau)/  !
   ! omega is the SCOPE leaf-angle asymmetry, so the diffuse backscatter is beta = 0.5*(1+g)   !
   ! and the direct-beam upscatter is beta0 = 0.5*(1 + g/k), k = G(mu)/mu. Leaf and wood are   !
   ! stored separately and blended per cohort by clumping-corrected area. `lidf` (per PFT) is    !
   ! kept so the mu-dependent G(mu) can be evaluated each timestep.                              !
   !---------------------------------------------------------------------------------------!
   type :: rad_pft_optics_t
      integer(ik) :: n_pft  = 0_ik
      integer(ik) :: n_band = 0_ik
      real(wp), allocatable :: omega_leaf(:,:)     !< (band,pft) single-scatter albedo, leaves
      real(wp), allocatable :: omega_wood(:,:)     !< (band,pft) single-scatter albedo, wood
      real(wp), allocatable :: g_leaf(:,:)         !< (band,pft) bf*(rho-tau)/omega, leaves
      real(wp), allocatable :: g_wood(:,:)         !< (band,pft) bf*(rho-tau)/omega, wood
      real(wp), allocatable :: clumping_leaf(:)    !< (pft) leaf clumping factor (0,1]
      real(wp), allocatable :: clumping_wood(:)    !< (pft) wood clumping factor (0,1]
      real(wp), allocatable :: lidf(:,:)           !< (class,pft) leaf-angle distribution weights
      real(wp), allocatable :: bf(:)               !< (pft) <cos^2(theta_leaf)>
      logical,  allocatable :: has_beam(:)         !< (band) band has a collimated (solar) beam
      logical,  allocatable :: has_emission(:)     !< (band) band has thermal emission (LW)
   end type rad_pft_optics_t

   !---------------------------------------------------------------------------------------!
   ! Per-band incident forcing + ground boundary (absolute W/m2). For emission bands the      !
   ! incident diffuse is the atmospheric downwelling (rlong) and grnd_emiss is the surface      !
   ! thermal source; for shortwave bands grnd_emiss = 0 and grnd_refl is the albedo.            !
   !---------------------------------------------------------------------------------------!
   type :: rad_forcing_t
      integer(ik) :: n_band = 0_ik
      real(wp)    :: cosz   = 1.0_wp               !< cosine of solar zenith / incidence (floored > 0)
      real(wp), allocatable :: incid_beam(:)       !< (band) [W/m2] direct-beam incident at canopy top
      real(wp), allocatable :: incid_diff(:)       !< (band) [W/m2] diffuse incident at canopy top
      real(wp), allocatable :: grnd_refl(:)        !< (band) ground reflectance (albedo, or 1-emiss)
      real(wp), allocatable :: grnd_emiss(:)       !< (band) [W/m2] ground thermal emission (0 for SW)
   end type rad_forcing_t

   !---------------------------------------------------------------------------------------!
   ! Returned fluxes. Per cohort, per band: radiation absorbed by leaves and by wood [W/m2 of   !
   ! ground]. Patch level: albedo (upward/incident) and the below-canopy downwelling fluxes.     !
   !---------------------------------------------------------------------------------------!
   type :: rad_flux_t
      integer(ik) :: n_band = 0_ik, n_coh = 0_ik
      real(wp), allocatable :: abs_leaf(:,:)       !< (band,coh) [W/m2] absorbed by leaves
      real(wp), allocatable :: abs_wood(:,:)       !< (band,coh) [W/m2] absorbed by wood
      real(wp), allocatable :: albedo(:)           !< (band) canopy+ground albedo (SW) / upward frac
      real(wp), allocatable :: dn_ground(:)        !< (band) [W/m2] downwelling below canopy (to ground)
      real(wp), allocatable :: up_ground(:)        !< (band) [W/m2] upwelling from ground into canopy
   end type rad_flux_t

contains

   subroutine alloc_rad_pft_optics(opt, n_band, n_pft, n_class)
      type(rad_pft_optics_t), intent(out) :: opt
      integer(ik),            intent(in)  :: n_band, n_pft, n_class
      opt%n_band = n_band ; opt%n_pft = n_pft
      allocate(opt%omega_leaf(n_band, n_pft), opt%omega_wood(n_band, n_pft))
      allocate(opt%g_leaf(n_band, n_pft),     opt%g_wood(n_band, n_pft))
      allocate(opt%clumping_leaf(n_pft),      opt%clumping_wood(n_pft))
      allocate(opt%lidf(n_class, n_pft),      opt%bf(n_pft))
      allocate(opt%has_beam(n_band),          opt%has_emission(n_band))
      opt%omega_leaf = 0.0_wp ; opt%omega_wood = 0.0_wp
      opt%g_leaf = 0.0_wp ; opt%g_wood = 0.0_wp
      opt%clumping_leaf = 1.0_wp ; opt%clumping_wood = 1.0_wp
      opt%lidf = 0.0_wp ; opt%bf = 0.0_wp
      opt%has_beam = .false. ; opt%has_emission = .false.
   end subroutine alloc_rad_pft_optics

   subroutine alloc_rad_forcing(f, n_band)
      type(rad_forcing_t), intent(out) :: f
      integer(ik),         intent(in)  :: n_band
      f%n_band = n_band
      allocate(f%incid_beam(n_band), f%incid_diff(n_band), f%grnd_refl(n_band), f%grnd_emiss(n_band))
      f%incid_beam = 0.0_wp ; f%incid_diff = 0.0_wp ; f%grnd_refl = 0.0_wp ; f%grnd_emiss = 0.0_wp
   end subroutine alloc_rad_forcing

   subroutine alloc_rad_flux(flux, n_band, n_coh)
      type(rad_flux_t), intent(out) :: flux
      integer(ik),      intent(in)  :: n_band, n_coh
      flux%n_band = n_band ; flux%n_coh = n_coh
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp ; flux%abs_wood = 0.0_wp
      flux%albedo = 0.0_wp ; flux%dn_ground = 0.0_wp ; flux%up_ground = 0.0_wp
   end subroutine alloc_rad_flux

end module meds_rad_types
