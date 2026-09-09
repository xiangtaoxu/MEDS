!==========================================================================================!
! meds_canopy_types -- the argument records of the canopy kernels: radiative transfer and    !
! aerodynamics. Pure DATA, no methods, no hidden state.                                      !
!                                                                                          !
! RT: the precomputed per-PFT optics table, the per-band incident forcing, the returned      !
! per-cohort absorbed radiation, and the ground optical state that closes the two-stream's    !
! lower boundary. Aerodynamics: the free-atmosphere + canopy-air forcing, the canopy geometry, !
! and the returned conductances and in-canopy wind profile.                                    !
!                                                                                          !
! Split out of the former meds_canopy_types/meds_soil_types, which held four domains' worth of records in  !
! one 581-line module -- the last instance of the grab-bag decision #8 exists to prevent. The   !
! run-config bundles these kernels switch on (aero_cfg_t, the selector codes) live one layer     !
! down in meds_biophysics_opts; the prognostic stores live in state/column.                      !
!==========================================================================================!
module meds_canopy_types
   use meds_kinds,     only : wp, ik
   use meds_constants, only : grav, cp_air
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t, ground_optics_state_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux
   public :: aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out, ensure_aero_out_capacity
   public :: set_aero_env_atm, set_aero_env_canopy

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
      real(wp), allocatable :: abs_leaf(:,:)       !< (band,col_cohort) [W/m2] absorbed by leaves
      real(wp), allocatable :: abs_wood(:,:)       !< (band,col_cohort) [W/m2] absorbed by wood
      real(wp), allocatable :: albedo(:)           !< (band) canopy+ground albedo (SW) / upward frac
      real(wp), allocatable :: dn_ground(:)        !< (band) [W/m2] downwelling below canopy (to ground)
      real(wp), allocatable :: up_ground(:)        !< (band) [W/m2] upwelling from ground into canopy
   end type rad_flux_t

   !---------------------------------------------------------------------------------------!
   ! Ground / surface optical state -- the two-stream lower boundary (bare-soil placeholder;   !
   ! only the soil fields are consulted now, the rest are reserved for the full surface model). !
   !---------------------------------------------------------------------------------------!
   type :: ground_optics_state_t
      integer(ik)           :: n_band = 0_ik
      real(wp), allocatable :: soil_albedo(:)    !< (band) shortwave soil albedo; unused for emission bands
      real(wp)              :: soil_emiss = 0.96_wp   !< thermal emissivity of the ground
      real(wp)              :: soil_temp  = 298.0_wp  !< [K] ground (skin) temperature
   end type ground_optics_state_t

   !----- Per-patch forcing + canopy-air-space state (read-only). ---------------------------!
   type :: aero_env_t
      real(wp) :: u_ref     = 2.0_wp                !< [m/s]      wind at reference height
      real(wp) :: zref      = 30.0_wp               !< [m]        reference (measurement) height
      real(wp) :: theta_atm = 298.15_wp             !< [K]        potential temp at zref
      real(wp) :: shv_atm   = 0.010_wp              !< [kg/kg]    specific humidity at zref
      real(wp) :: co2_atm   = 400.0_wp              !< [umol/mol] free-atmosphere CO2
      real(wp) :: press     = 101325.0_wp           !< [Pa]
      real(wp) :: rho_air   = 1.2_wp                !< [kg/m3]    (diagnostic passthrough)
      real(wp) :: can_theta = 298.15_wp             !< [K]        CAS potential temp
      real(wp) :: can_temp  = 298.15_wp             !< [K]        CAS actual temp (buoyancy Grashof)
      real(wp) :: can_shv   = 0.010_wp              !< [kg/kg]    CAS specific humidity
      real(wp) :: can_co2   = 400.0_wp              !< [umol/mol] CAS CO2
      real(wp) :: t_ground  = 298.15_wp             !< [K]        ground skin temp (ground-conductance stability)
   end type aero_env_t

   !----- Per-patch canopy geometry. -------------------------------------------------------!
   type :: aero_geom_t
      real(wp) :: veg_height   = 20.0_wp            !< [m]  canopy top height
      real(wp) :: opencan_frac = 0.0_wp             !< [-]  open-sky fraction (0 = closed canopy)
      real(wp) :: snowfac      = 0.0_wp             !< [-]  snow burial fraction of the canopy
   end type aero_geom_t

   !----- Outputs: per-patch scalars + per-cohort arrays (caller-owned; written in place). ---!
   type :: aero_out_t
      integer(ik) :: n_coh = 0_ik
      real(wp) :: ustar = 0.0_wp, tstar = 0.0_wp, qstar = 0.0_wp, cstar = 0.0_wp   !< [m/s],[K],[kg/kg],[umol/mol]
      real(wp) :: temp1 = 0.0_wp, temp2 = 0.0_wp    !< scalar profile factors (g_atm_heat=rho*ustar*temp1, g_atm_vapour=..temp2)
      real(wp) :: zeta = 0.0_wp, rib = 0.0_wp, obu = 0.0_wp   !< stability diagnostics
      real(wp) :: ggbare = 0.0_wp, ggveg = 0.0_wp, ggnet = 0.0_wp   !< [m/s] ground conductances (r_aero = 1/ggnet)
      real(wp) :: rough = 0.0_wp, displace = 0.0_wp, can_depth = 0.0_wp   !< [m]
      real(wp) :: uh = 0.0_wp                        !< [m/s] canopy-top wind (diagnostic)
      real(wp), allocatable :: wind(:)               !< [m/s] per-cohort in-canopy wind
      real(wp), allocatable :: leaf_gbh(:), leaf_gbw(:)   !< [m/s] leaf boundary-layer heat/vapour conductance
      real(wp), allocatable :: wood_gbh(:), wood_gbw(:)   !< [m/s] wood boundary-layer heat/vapour conductance
   end type aero_out_t

contains

   subroutine alloc_rad_pft_optics(optics, n_band, n_pft, n_class)
      type(rad_pft_optics_t), intent(out) :: optics
      integer(ik),            intent(in)  :: n_band, n_pft, n_class
      optics%n_band = n_band
      optics%n_pft  = n_pft
      allocate(optics%omega_leaf(n_band, n_pft), optics%omega_wood(n_band, n_pft))
      allocate(optics%g_leaf(n_band, n_pft),     optics%g_wood(n_band, n_pft))
      allocate(optics%clumping_leaf(n_pft),      optics%clumping_wood(n_pft))
      allocate(optics%lidf(n_class, n_pft),      optics%bf(n_pft))
      allocate(optics%has_beam(n_band),          optics%has_emission(n_band))
      optics%omega_leaf = 0.0_wp
      optics%omega_wood = 0.0_wp
      optics%g_leaf = 0.0_wp
      optics%g_wood = 0.0_wp
      optics%clumping_leaf = 1.0_wp
      optics%clumping_wood = 1.0_wp
      optics%lidf = 0.0_wp
      optics%bf = 0.0_wp
      optics%has_beam = .false.
      optics%has_emission = .false.
   end subroutine alloc_rad_pft_optics

   subroutine alloc_rad_forcing(f, n_band)
      type(rad_forcing_t), intent(out) :: f
      integer(ik),         intent(in)  :: n_band
      f%n_band = n_band
      allocate(f%incid_beam(n_band), f%incid_diff(n_band), f%grnd_refl(n_band), f%grnd_emiss(n_band))
      f%incid_beam = 0.0_wp
      f%incid_diff = 0.0_wp
      f%grnd_refl  = 0.0_wp
      f%grnd_emiss = 0.0_wp
   end subroutine alloc_rad_forcing

   subroutine alloc_rad_flux(flux, n_band, n_coh)
      type(rad_flux_t), intent(out) :: flux
      integer(ik),      intent(in)  :: n_band, n_coh
      flux%n_band = n_band
      flux%n_coh  = n_coh
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp
      flux%abs_wood = 0.0_wp
      flux%albedo = 0.0_wp
      flux%dn_ground = 0.0_wp
      flux%up_ground = 0.0_wp
   end subroutine alloc_rad_flux

   !----- Allocate the per-cohort output arrays of an aero_out_t (the kernel writes in place). !
   subroutine alloc_aero_out(out, n_coh)
      type(aero_out_t), intent(out) :: out
      integer(ik),      intent(in)  :: n_coh
      out%n_coh = n_coh
      allocate(out%wind(n_coh), out%leaf_gbh(n_coh), out%leaf_gbw(n_coh),                      &
               out%wood_gbh(n_coh), out%wood_gbw(n_coh))
      out%wind = 0.0_wp
      out%leaf_gbh = 0.0_wp
      out%leaf_gbw = 0.0_wp
      out%wood_gbh = 0.0_wp
      out%wood_gbw = 0.0_wp
   end subroutine alloc_aero_out

   !----- Grow-only capacity check for aero_out_t (mirrors ensure_column_cohort_capacity,       !
   !      MEDS_NUMERICS_SCOPING.md BB1 phase 1): aero is a per-call scratch (canopy_aerodynamics  !
   !      overwrites every active element each call), so reusing over-sized capacity is safe.    !
   subroutine ensure_aero_out_capacity(out, n_coh)
      type(aero_out_t), intent(inout) :: out
      integer(ik),       intent(in)    :: n_coh
      if (.not. allocated(out%wind)) then
         call alloc_aero_out(out, n_coh)
      else if (size(out%wind) < n_coh) then
         call alloc_aero_out(out, n_coh)
      else
         out%n_coh = n_coh
      end if
   end subroutine ensure_aero_out_capacity

   !---------------------------------------------------------------------------------------!
   ! Fill `aero_env_t`'s ATMOSPHERIC reference state from the met driver's ACTUAL air        !
   ! temperature. THE one place the potential-temperature conversion lives.                  !
   !                                                                                          !
   ! WHY THIS EXISTS (issue #97). `theta_atm` is the reference the Monin-Obukhov solve measures !
   ! the canopy against, and `aero_env_t` gives it a plausible 298.15 K DEFAULT. Setting only    !
   ! `forc%air_temp` therefore leaves MO comparing the canopy to a fixed 298.15 K -- which is not an  !
   ! error, just wrong: with a 282-294 K forcing it manufactures a permanent STABLE layer and     !
   ! FLOORS `ustar` at 0.1, a ~44x suppression of turbulent exchange with the wrong sign of       !
   ! stratification. Every column test and every hand-built probe had exactly this omission, and  !
   ! the artifact is convincing -- it looks like a cadence-dependent decoupling bifurcation in    !
   ! the integrator, not like a missing assignment.                                               !
   !                                                                                          !
   ! The defect was that fixtures assembled `aenv` by a DIFFERENT route than production, so the   !
   ! fix is a shared routine rather than a corrected copy: `meds_fast_dynamics::fill_aenv` and    !
   ! every test now go through this. `zref` must already be set.                                  !
   !                                                                                          !
   ! Approximation (inherited from fill_aenv): the shallow-layer dry-adiabatic form theta =        !
   ! T + (g/cp)*z, ignoring displacement height. A proper met driver would use zref - displace.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_aero_env_atm(aenv, air_temp, shv_atm, co2_atm)
      type(aero_env_t), intent(inout) :: aenv
      real(wp),         intent(in)    :: air_temp   !< [K] ACTUAL air temperature at zref (not potential)
      real(wp),         intent(in)    :: shv_atm    !< [kg/kg]    specific humidity at zref
      real(wp),         intent(in)    :: co2_atm    !< [umol/mol] free-atmosphere CO2
      aenv%theta_atm = air_temp + (grav / cp_air) * aenv%zref
      aenv%shv_atm   = shv_atm
      aenv%co2_atm   = co2_atm
   end subroutine set_aero_env_atm

   !----- ...and the CANOPY/ground half, refreshed from the prognostic state each sub-step.   !
   !      Split from the atmospheric half because they have different cadences and different     !
   !      sources: the atmosphere comes from the met driver, this comes from the model's own      !
   !      evolving state. Both must be refreshed -- a stale CAS temperature mis-prices the same   !
   !      stability solve that a stale `theta_atm` does. -----------------------------------------!
   pure subroutine set_aero_env_canopy(aenv, can_temp, can_shv, can_co2, t_ground)
      type(aero_env_t), intent(inout) :: aenv
      real(wp),         intent(in)    :: can_temp, can_shv, can_co2, t_ground
      aenv%can_theta = can_temp ; aenv%can_temp = can_temp
      aenv%can_shv   = can_shv  ; aenv%can_co2  = can_co2
      aenv%t_ground  = t_ground
   end subroutine set_aero_env_canopy

end module meds_canopy_types
