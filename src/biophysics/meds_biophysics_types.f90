!==========================================================================================!
! meds_biophysics_types -- shared derived types of the biophysics domain (fast, sub-daily     !
! physical processes over the ecosystem column). Currently the canopy radiative-transfer        !
! interface; the intended home for future energy-balance and hydrology types as those land.     !
!                                                                                          !
! Pure DATA (no methods, no hidden state), the RT analogue of meds_plant_types. Three types:  !
!   * rad_pft_optics_t -- the PRECOMPUTED, MU-INDEPENDENT per-PFT optics table (single-scatter !
!                         albedo omega and the leaf-angle scattering asymmetry g per band and   !
!                         tissue, plus the leaf-angle distribution lidf and its 2nd moment bf).   !
!                         Filled ONCE by meds_optics%derive_rad_optics.                            !
!   * rad_forcing_t    -- the per-band incident fluxes (beam & diffuse), the solar cosine, and    !
!                         the ground reflectance/emission (from meds_optics%ground_optics).          !
!                         Absolute W/m2 throughout -- no ED2-style normalize-to-one.                  !
!   * rad_flux_t       -- the returned per-cohort absorbed radiation (leaf & wood, per band) plus    !
!                         the patch-level albedo and below-canopy transmission diagnostics.           !
!                                                                                          !
! Bands are configurable; the default set is {VIS, NIR, LW}. Every band carries a thermal-      !
! emission term, identically zero for VIS/NIR (has_emission = .false.).                          !
!==========================================================================================!
module meds_biophysics_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux

   !----- Soil-column hydrology (see meds_column_hydrology / meds_soil_parameters). ----!
   public :: n_soil_layer_max
   public :: SOIL_SOLVER_BE
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   public :: SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER, SOIL_BC_BEDROCK, SOIL_BC_SLOPE
   public :: SOIL_LIN_FROZEN, SOIL_LIN_PICARD
   public :: SOIL_SUBSTEP_ADAPTIVE, SOIL_SUBSTEP_FIXED
   public :: soil_column_t, chydro_forcing_t, soil_params_t, soil_opts_t, chydro_flux_t

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

   !=======================================================================================!
   !  Soil-column hydrology types + selector codes (meds_column_hydrology, design §4).  !
   !  Fixed-size (n_soil_layer_max) so the kernel stays allocatable-free and GPU-eligible.   !
   !  ED2 negative-z convention: elevation z <= 0 below ground; dz, dz_node are positive       !
   !  magnitudes. The SOIL_* codes live here for the standalone P0/P1 cut and migrate to        !
   !  meds_config when the TOML string->enum mapping lands at P3.                                !
   !=======================================================================================!
   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time column-depth ceiling

   integer(ik), parameter :: SOIL_SOLVER_BE = 1_ik         !< linearly-implicit backward-Euler (only solver)

   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   integer(ik), parameter :: SOIL_BC_FREE_DRAIN = 1_ik     !< unit-gradient q = K(theta_N)  (MVP default)
   integer(ik), parameter :: SOIL_BC_AQUIFER    = 2_ik     !< SIMTOP aquifer / water table (P2)
   integer(ik), parameter :: SOIL_BC_BEDROCK    = 3_ik     !< no-flow  q = 0
   integer(ik), parameter :: SOIL_BC_SLOPE      = 4_ik     !< slope-reduced drainage (P2)

   integer(ik), parameter :: SOIL_LIN_FROZEN = 1_ik        !< frozen-coefficient single solve (MVP)
   integer(ik), parameter :: SOIL_LIN_PICARD = 2_ik        !< Celia modified-Picard (P2)

   integer(ik), parameter :: SOIL_SUBSTEP_ADAPTIVE = 1_ik  !< adaptive step-doubling
   integer(ik), parameter :: SOIL_SUBSTEP_FIXED    = 2_ik  !< fixed count (GPU warp-uniform)

   !----- Prognostic per-patch soil column (the one mutable value the kernel updates). -----!
   type :: soil_column_t
      real(wp) :: theta(n_soil_layer_max) = 0.0_wp   !< [m3/m3] volumetric soil moisture (PROGNOSTIC)
      real(wp) :: w_surface = 0.0_wp                 !< [kg/m2] ponded surface water
      real(wp) :: w_aquifer = 0.0_wp                 !< [kg/m2] lumped aquifer store (SOIL_BC_AQUIFER, P2)
      real(wp) :: z_wt      = 0.0_wp                 !< [m] water-table elevation (<= 0; P2)
   end type soil_column_t

   !----- Soil-column boundary conditions (read-only). ------------------------------------!
   type :: chydro_forcing_t
      real(wp) :: precip_ground = 0.0_wp                  !< [kg/m2/s] ground-reaching liquid (post interception)
      real(wp) :: root_uptake(n_soil_layer_max) = 0.0_wp  !< [kg/m2/s] per-layer transpiration DEMAND (x nplant)
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temp (FORCED = T_air until soil energy)
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance (soil evap series)
   end type chydro_forcing_t

   !----- Per-column geometry + texture (assembled once per site; ED2 negative-z). ---------!
   type :: soil_params_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0, ED2 slz)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      real(wp) :: vg_alpha(n_soil_layer_max)  = 0.0_wp       !< [1/m] van Genuchten inverse air-entry
      real(wp) :: vg_n(n_soil_layer_max)      = 0.0_wp       !< [-] van Genuchten pore-size index (> 1)
      real(wp) :: psi_sat(n_soil_layer_max)   = 0.0_wp       !< [m] Campbell air-entry potential (option)
      real(wp) :: b_camp(n_soil_layer_max)    = 0.0_wp       !< [-] Campbell exponent (option)
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
      real(wp) :: theta_wp(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] wilting point (DERIVED)
      real(wp) :: root_frac(n_soil_layer_max) = 0.0_wp       !< [-] normalized root fraction (sum = 1)
   end type soil_params_t

   !----- Pre-extracted solver selectors + tolerances (NOT the whole config). --------------!
   type :: soil_opts_t
      integer(ik) :: solver    = SOIL_SOLVER_BE
      integer(ik) :: bottom_bc = SOIL_BC_FREE_DRAIN
      integer(ik) :: linearize = SOIL_LIN_FROZEN
      integer(ik) :: substep   = SOIL_SUBSTEP_ADAPTIVE
      logical     :: zeng_decker = .false.               !< P2 only; MVP plain-gravity flux
      real(wp)    :: rtol   = 1.0e-3_wp
      real(wp)    :: atol   = 1.0e-4_wp
      real(wp)    :: h_init = 900.0_wp                   !< [s] initial substep
      integer(ik) :: max_substep = 200_ik
      integer(ik) :: max_picard  = 5_ik
      real(wp)    :: w_pond_max = 5.0_wp                 !< [kg/m2] ponding capacity before surface runoff
      real(wp)    :: dewmx      = 0.1_wp                 !< [kg/m2 per PAI] canopy storage capacity
      real(wp)    :: intercept_alpha = 1.0_wp            !< Beer interception efficiency
      real(wp)    :: intercept_k     = 0.5_wp            !< [1/PAI] Beer extinction
      real(wp)    :: dsl_dmax       = 0.015_wp           !< [m] DSL max thickness (soil evap)
      real(wp)    :: dsl_theta_init = 0.8_wp             !< DSL initiation theta_1/theta_sat
      real(wp)    :: psi_wilt = -152.96_wp               !< [m] wilting head (~ -1.5 MPa) for the sink f_wilt
      real(wp)    :: psi_open = -3.37_wp                 !< [m] onset head (~ -0.033 MPa) for the f_wilt ramp
      logical     :: debug_error = .false.
   end type soil_opts_t

   !----- Soil-column outputs + diagnostics. ----------------------------------------------!
   type :: chydro_flux_t
      real(wp) :: infiltration   = 0.0_wp                !< [kg/m2/s] top-face infiltration
      real(wp) :: drainage       = 0.0_wp                !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf    = 0.0_wp                !< [kg/m2/s] surface runoff
      real(wp) :: soil_evap      = 0.0_wp                !< [kg/m2/s] ground evaporation
      real(wp) :: uptake_total   = 0.0_wp                !< [kg/m2/s] realized root uptake (after theta_wp cap)
      real(wp) :: uptake_deficit = 0.0_wp                !< [kg/m2/s] capped (unmet) sink
      real(wp) :: clip_excess    = 0.0_wp                !< [kg/m2/s] theta-clip water routed to ponding
      real(wp) :: psi_soil(n_soil_layer_max) = 0.0_wp    !< [MPa] per-layer matric potential (EXPORTED to hydraulics)
      real(wp) :: mass_resid     = 0.0_wp                !< [kg/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik                         !< sub-steps taken
      logical  :: converged = .true.                     !< .false. on any cap-hit
   end type chydro_flux_t

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

end module meds_biophysics_types
