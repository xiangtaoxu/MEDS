! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_canopy_types -- the argument records of the canopy kernels, pure data. Radiative    !
! transfer: the per-PFT optics table, the per-band forcing, the returned fluxes and the    !
! ground's optical state (docs/science/canopy_radiation_transfer.md). Aerodynamics: the    !
! atmosphere and canopy-air state, the canopy geometry, and the returned conductances and  !
! wind (docs/science/canopy_aerodynamics.md). set_aero_env_atm and set_aero_env_canopy     !
! are the one way to fill aero_env_t; the aerodynamics options (aero_params_t) are in      !
! meds_biophysics_params.                                                                  !
!==========================================================================================!
module meds_canopy_types
   use meds_kinds,     only : wp, ik
   implicit none
   private

   public :: RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   public :: rad_pft_optics_t, rad_forcing_t, rad_flux_t, ground_optics_state_t
   public :: alloc_rad_pft_optics, alloc_rad_forcing, alloc_rad_flux
   public :: aero_env_t, aero_geom_t, aero_out_t, ensure_aero_out_capacity
   public :: set_aero_env_atm, set_aero_env_canopy
   public :: alloc_aero_out   ! used only by the tests

   !----- Default three-band layout (indices into the band dimension). --------------------!
   integer(ik), parameter :: RAD_VIS = 1_ik   !< visible / PAR (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_NIR = 2_ik   !< near-infrared (beam + diffuse, no emission)
   integer(ik), parameter :: RAD_LW  = 3_ik   !< thermal / longwave (diffuse + emission, no beam)
   integer(ik), parameter :: N_RAD_BAND_DEFAULT = 3_ik

   !---------------------------------------------------------------------------------------!
   ! Per-PFT optics (docs/science/canopy_radiation_transfer.md eqs. (4)-(6)), leaf and     !
   ! wood apart, blended per cohort by clumped area; lidf is kept to evaluate G(mu).       !
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

   !----- Per-band incident forcing and ground boundary [W/m2]. In the thermal band the diffuse is   !
   !      the sky's longwave and grnd_emiss the ground's emission; in shortwave grnd_emiss = 0. -----!
   type :: rad_forcing_t
      integer(ik) :: n_band = 0_ik
      real(wp)    :: cosz   = 1.0_wp               !< cosine of solar zenith / incidence (floored > 0)
      real(wp), allocatable :: incid_beam(:)       !< (band) [W/m2] direct-beam incident at canopy top
      real(wp), allocatable :: incid_diff(:)       !< (band) [W/m2] diffuse incident at canopy top
      real(wp), allocatable :: grnd_refl(:)        !< (band) ground reflectance (albedo, or 1-emiss)
      real(wp), allocatable :: grnd_emiss(:)       !< (band) [W/m2] ground thermal emission (0 for SW)
   end type rad_forcing_t

   !----- Returned fluxes: per cohort and band, absorbed by leaves and wood [W/m2 of ground]; per  !
   !      band, the albedo and the fluxes below the canopy. ----------------------------------!
   type :: rad_flux_t
      integer(ik) :: n_band = 0_ik, n_coh = 0_ik
      real(wp), allocatable :: abs_leaf(:,:)       !< (band,col_cohort) [W/m2] absorbed by leaves
      real(wp), allocatable :: abs_wood(:,:)       !< (band,col_cohort) [W/m2] absorbed by wood
      real(wp), allocatable :: incid_top(:,:)      !< (band,col_cohort) [W/m2] beam + down diffuse at its layer's top
      real(wp), allocatable :: albedo(:)           !< (band) canopy+ground albedo (SW) / upward frac
      real(wp), allocatable :: dn_ground(:)        !< (band) [W/m2] downwelling below canopy (to ground)
      real(wp), allocatable :: up_ground(:)        !< (band) [W/m2] upwelling from ground into canopy
   end type rad_flux_t

   !----- The ground's optical state, the two-stream's lower boundary                              !
   !      (docs/science/canopy_radiation_transfer.md §4): bare soil, blended toward the snow. ------!
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
      real(wp) :: theta_atm = 298.15_wp             !< [K]        potential temp, referenced to zref (= the air temperature there)
      real(wp) :: shv_atm   = 0.010_wp              !< [kg/kg]    specific humidity at zref
      real(wp) :: co2_atm   = 400.0_wp              !< [umol/mol] free-atmosphere CO2
      real(wp) :: press     = 101325.0_wp           !< [Pa]
      real(wp) :: rho_air   = 1.2_wp                !< [kg/m3]    (diagnostic passthrough)
      real(wp) :: can_theta = 298.15_wp             !< [K]        CAS potential temp, referenced to zref (= can_temp)
      real(wp) :: can_temp  = 298.15_wp             !< [K]        CAS actual temp (buoyancy Grashof)
      real(wp) :: can_shv   = 0.010_wp              !< [kg/kg]    CAS specific humidity
      real(wp) :: can_co2   = 400.0_wp              !< [umol/mol] CAS CO2
      real(wp) :: t_ground  = 298.15_wp             !< [K]        ground skin temp (ground-conductance stability)
   end type aero_env_t

   !----- Per-patch canopy geometry. -------------------------------------------------------!
   type :: aero_geom_t
      real(wp) :: veg_height   = 20.0_wp            !< [m]  canopy top height
      real(wp) :: opencan_frac = 0.0_wp             !< [-]  open fraction W = exp(-(LAI+WAI)): 1 bare, 0 dense
      real(wp) :: snowfac      = 0.0_wp             !< [-]  snow burial fraction of the canopy
   end type aero_geom_t

   !----- Outputs: per-patch scalars + per-cohort arrays (caller-owned; written in place). ---!
   type :: aero_out_t
      integer(ik) :: n_coh = 0_ik
      real(wp) :: ustar = 0.0_wp, tstar = 0.0_wp   !< [m/s],[K]
      real(wp) :: temp1 = 0.0_wp, temp2 = 0.0_wp    !< scalar profile factors (g_atm_heat=rho*ustar*temp1, g_atm_vapour=..temp2)
      real(wp) :: zeta = 0.0_wp, rib = 0.0_wp, obu = 0.0_wp   !< stability diagnostics
      real(wp) :: ggbare = 0.0_wp, ggveg = 0.0_wp, ggnet = 0.0_wp   !< [m/s] ground conductances (r_aero = 1/ggnet)
      real(wp) :: rough = 0.0_wp, displace = 0.0_wp   !< [m]
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
      allocate(flux%abs_leaf(n_band, n_coh), flux%abs_wood(n_band, n_coh), flux%incid_top(n_band, n_coh))
      allocate(flux%albedo(n_band), flux%dn_ground(n_band), flux%up_ground(n_band))
      flux%abs_leaf = 0.0_wp
      flux%abs_wood = 0.0_wp
      flux%incid_top = 0.0_wp
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

   !----- Grow-only capacity of aero_out_t: canopy_aerodynamics overwrites every active element. --!
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
   ! Fill aero_env_t's atmospheric state from the forcing at zref, the canopy-air top,     !
   ! where potential temperature equals the actual air temperature                         !
   ! (docs/science/canopy_aerodynamics.md §2, The reference level). The driver and the     !
   ! column tests all fill it here: a theta_atm left at its 298.15 K default gives the     !
   ! Monin-Obukhov solve a false stratification.                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_aero_env_atm(aenv, air_temp, shv_atm, co2_atm)
      type(aero_env_t), intent(inout) :: aenv
      real(wp),         intent(in)    :: air_temp   !< [K] actual air temperature at zref
      real(wp),         intent(in)    :: shv_atm    !< [kg/kg]    specific humidity at zref
      real(wp),         intent(in)    :: co2_atm    !< [umol/mol] free-atmosphere CO2
      aenv%theta_atm = air_temp
      aenv%shv_atm   = shv_atm
      aenv%co2_atm   = co2_atm
   end subroutine set_aero_env_atm

   !----- The canopy-air and ground half, from the model's own state each sub-step; a stale        !
   !      canopy-air temperature mis-prices the stability solve as a stale theta_atm does. -------!
   pure subroutine set_aero_env_canopy(aenv, can_temp, can_shv, can_co2, t_ground)
      type(aero_env_t), intent(inout) :: aenv
      real(wp),         intent(in)    :: can_temp, can_shv, can_co2, t_ground
      aenv%can_theta = can_temp
      aenv%can_temp = can_temp
      aenv%can_shv   = can_shv
      aenv%can_co2  = can_co2
      aenv%t_ground  = t_ground
   end subroutine set_aero_env_canopy

end module meds_canopy_types
