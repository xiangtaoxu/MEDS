!==========================================================================================!
! meds_forcing_config -- the [forcing]/[site] configuration for meteorological forcing, plus  !
! the selector codes shared by the config, the reader, and the disaggregation kernels.         !
!                                                                                          !
! Placed in src/shared (NOT src/forcing) so meds_config -- the DAG ROOT -- can carry            !
! forcing_config_t as a plain-scalar component with NO backward `shared -> forcing` edge (the    !
! same rule the energy design used for soil_thermal_params_t). Pure DATA + parameters; links      !
! meds_kinds only. Defaults are the Ithaca NY / ERA5-Land reference site (design                   !
! MEDS_FORCING_DESIGN.md sections 3.3, 6.6).                                                        !
!==========================================================================================!
module meds_forcing_config
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: forcing_config_t
   public :: MET_BACKEND_CONST, MET_BACKEND_NETCDF
   public :: METAVG_INSTANT, METAVG_END, METAVG_BEGIN, METAVG_CENTER
   public :: SWPART_PASSTHROUGH, SWPART_WEISS_NORMAN, SWPART_SIB, SWPART_CLEARIDX
   public :: LW_FILE, LW_SYNTHESIZE
   public :: CLAMP_ERROR, CLAMP_HOLD
   public :: INTERP_LINEAR, INTERP_STEP, INTERP_COSZ
   public :: GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST

   !----- Reader backend: a MEDS forcing NetCDF, or a no-file reference-climate box. --------!
   integer(ik), parameter :: MET_BACKEND_CONST  = 0_ik   !< no file: met_forcing_t defaults (reference climate)
   integer(ik), parameter :: MET_BACKEND_NETCDF = 1_ik   !< the MEDS multi-grid forcing NetCDF (the file format)

   !----- Timestamp semantics of a forcing record (avg_convention). ------------------------!
   integer(ik), parameter :: METAVG_INSTANT = 0_ik   !< value is instantaneous AT the stamp
   integer(ik), parameter :: METAVG_END     = 1_ik   !< flux mean over the interval ENDING at the stamp (ERA5-Land)
   integer(ik), parameter :: METAVG_BEGIN   = 2_ik   !< mean over the interval BEGINNING at the stamp
   integer(ik), parameter :: METAVG_CENTER  = 3_ik   !< mean over the interval CENTERED on the stamp

   !----- Shortwave partition of total SWdown into the four (beam/diffuse)x(PAR/NIR) streams. !
   integer(ik), parameter :: SWPART_PASSTHROUGH  = 0_ik   !< file already carries the four streams
   integer(ik), parameter :: SWPART_CLEARIDX     = 1_ik   !< clearness-index (Erbs) split -- the ERA5-Land P0 default
   integer(ik), parameter :: SWPART_WEISS_NORMAN = 2_ik   !< Weiss-Norman band-specific (P1)
   integer(ik), parameter :: SWPART_SIB          = 3_ik   !< SiB (Sellers 1986) (reserved)

   !----- Downwelling longwave source. -----------------------------------------------------!
   integer(ik), parameter :: LW_FILE       = 0_ik   !< read from file (ERA5-Land has strd)
   integer(ik), parameter :: LW_SYNTHESIZE = 1_ik   !< Brutsaert/Idso clear-sky synthesis (source lacking LW)

   !----- Model-start-before-first-record policy. ------------------------------------------!
   integer(ik), parameter :: CLAMP_ERROR = 0_ik   !< hard stop
   integer(ik), parameter :: CLAMP_HOLD  = 1_ik   !< clamp to record #1 (w_next = 0)

   !----- Per-variable temporal-interpolation policy (used by interpolate_forcing). ---------!
   integer(ik), parameter :: INTERP_LINEAR = 0_ik   !< linear in the window (state vars)
   integer(ik), parameter :: INTERP_STEP   = 1_ik   !< step-constant (precip; hold prev)
   integer(ik), parameter :: INTERP_COSZ   = 2_ik   !< cosz-weighted (shortwave; handled by disaggregate_shortwave)

   !----- How a polygon binds to the file's `grid` dimension (multi-polygon P2 subset). -----!
   integer(ik), parameter :: GRIDMATCH_EXPLICIT = 0_ik  !< use grid_index verbatim (default; current behaviour)
   integer(ik), parameter :: GRIDMATCH_NEAREST  = 1_ik  !< pick the grid cell nearest [site] lat/lon (great-circle)

   !==========================================================================================!
   !  The [forcing]/[site] block. Plain scalars (no allocatables), so meds_config carries it     !
   !  trivially. NOTE: there is NO gap_policy -- MEDS never gap-fills; a missing required value    !
   !  is a hard error (design comment 2 / §5.5).                                                   !
   !==========================================================================================!
   type :: forcing_config_t
      logical            :: forcing_on   = .false.               !< master gate (indep. of fast_biophysics_on)
      integer(ik)        :: backend      = MET_BACKEND_NETCDF    !< "netcdf" | "const"
      character(len=256) :: path         = ''                    !< forcing NetCDF path
      integer(ik)        :: grid_index   = 1_ik                  !< which (time,grid) location this polygon reads
      real(wp)           :: dt_forcing   = 3600.0_wp             !< [s] native interval (hourly ERA5-Land)
      integer(ik)        :: avg_convention = METAVG_END          !< flux vars mean over the hour ENDING at the stamp
      integer(ik)        :: sw_partition = SWPART_CLEARIDX        !< ERA5-Land total SW -> partition required
      integer(ik)        :: lwdown_source = LW_FILE              !< file (ERA5-Land strd) | synthesize
      real(wp)           :: co2_const    = 420.0_wp              !< [umol/mol] ERA5-Land has no CO2 (single authority)
      real(wp)           :: rad_sw_ground_const = 60.0_wp        !< [W/m2] CONST-backend ground SW
      logical            :: recycle      = .true.                !< cycle the record when the run outruns the file
      integer(ik)        :: start_clamp  = CLAMP_ERROR           !< model start < base_time: error | hold record #1
      integer(ik)        :: grid_match   = GRIDMATCH_EXPLICIT    !< explicit grid_index | nearest [site] lat/lon (§4.1)
      !----- site geolocation ([site]) -- solar geometry + lapse. ----------------------------!
      real(wp)           :: latitude_deg  = 42.44_wp             !< [deg +N] Ithaca NY
      real(wp)           :: longitude_deg = -76.50_wp            !< [deg +E] Ithaca NY (solar time, §5.1)
      real(wp)           :: utc_offset_h  = 0.0_wp               !< [h] ERA5-Land is UTC -> offset 0
      logical            :: apply_solar_longitude = .true.       !< UTC file -> longitude gives local solar time
      real(wp)           :: elevation_m   = 320.0_wp             !< [m] site elevation (Ithaca ~320 m; elevation lapse)
      real(wp)           :: reference_height = 40.0_wp           !< [m] forcing reference height (> every PFT hgt_max)
      real(wp)           :: wind_meas_height = 10.0_wp           !< [m] ERA5-Land wind is 10 m (§5.2/§10)
      !----- wind-height + elevation-lapse corrections (§5.2/§10-Q2; both OFF by default). ----!
      logical            :: apply_wind_profile    = .false.      !< lift wind from wind_meas_height to reference_height
      logical            :: apply_elevation_lapse = .false.      !< lapse T/P from the grid-cell elevation to the site
      real(wp)           :: wind_roughness_z0 = 0.1_wp           !< [m] roughness length for the neutral-log wind profile
      real(wp)           :: lapse_rate_tair   = 0.0065_wp        !< [K/m] environmental lapse (positive = cooling upward)
      real(wp)           :: grid_elevation_m  = 320.0_wp         !< [m] elevation of the forcing source grid cell
   end type forcing_config_t

end module meds_forcing_config
