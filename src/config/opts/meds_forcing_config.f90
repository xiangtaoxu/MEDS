! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_forcing_config -- the [forcing]/[site] configuration for meteorological forcing     !
! (docs/science/forcing.md), plus the selector codes shared by the config, the reader, and !
! the disaggregation kernels.                                                              !
!                                                                                          !
! Placed in src/config (NOT src/io/input/forcing) so meds_config can carry                 !
! forcing_config_t as a plain-scalar component with no edge from the config layer back to  !
! the forcing layer. Pure data and parameters; uses meds_kinds and meds_time only. The     !
! defaults are an Ithaca NY / ERA5-Land reference site.                                    !
!==========================================================================================!
module meds_forcing_config
   use meds_kinds, only : wp, ik
   use meds_time,  only : meds_time_t
   implicit none
   private

   public :: forcing_config_t
   public :: MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT, MET_BACKEND_ED_ERA5LAND
   public :: MET_PATH_LEN, ARCHIVE_DT_SEC, MAX_RECYCLE_YEARS
   public :: METAVG_END, METAVG_BEGIN, METAVG_CENTER
   public :: SWPART_PASSTHROUGH, SWPART_WEISS_NORMAN, SWPART_CLEARIDX
   public :: LW_FILE, LW_SYNTHESIZE, LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO
   public :: CLAMP_ERROR, CLAMP_HOLD
   public :: INTERP_LINEAR, INTERP_STEP
   public :: GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   public :: CO2_SOURCE_CONST, CO2_SOURCE_FILE
   public :: HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND, WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL

   !----- Reader backend ([forcing].format): the single MEDS forcing file, the global ED_ERA5land  !
   !      archive, or a no-file reference-climate box. ------------------------------------------!
   !----- The ED_ERA5land archive's record spacing: hourly, end-stamped. -----------------------!
   real(wp), parameter :: ARCHIVE_DT_SEC = 3600.0_wp        !< [s]
   !----- Upper bound on the declared forcing recycle window, in whole calendar years: the search  !
   !      bound of the whole-year test (no met record spans centuries). -------------------------------!
   integer(ik), parameter :: MAX_RECYCLE_YEARS = 200_ik
   integer(ik), parameter :: MET_BACKEND_CONST       = 0_ik  !< no file: met_forcing_t defaults (reference climate)
   integer(ik), parameter :: MET_BACKEND_ED_DEFAULT  = 1_ik  !< the multi-grid forcing NetCDF (format = "ED_default")
   integer(ik), parameter :: MET_BACKEND_ED_ERA5LAND = 2_ik  !< the monthly ED_ERA5land archive (format = "ED_ERA5land")

   integer, parameter :: MET_PATH_LEN = 1024                !< length of every forcing path field

   !----- Timestamp semantics of a forcing record (avg_convention; docs/science/forcing.md §5). -!
   integer(ik), parameter :: METAVG_END     = 1_ik   !< flux mean over the interval ENDING at the stamp (ERA5-Land)
   integer(ik), parameter :: METAVG_BEGIN   = 2_ik   !< mean over the interval BEGINNING at the stamp
   integer(ik), parameter :: METAVG_CENTER  = 3_ik   !< mean over the interval CENTERED on the stamp

   !----- Shortwave partition of total SWdown into the four (beam/diffuse)x(PAR/NIR) streams       !
   !      (docs/science/forcing.md §6). -----------------------------------------------------------!
   integer(ik), parameter :: SWPART_PASSTHROUGH  = 0_ik   !< file already carries the four streams
   integer(ik), parameter :: SWPART_CLEARIDX     = 1_ik   !< clearness-index (Erbs) split, eq. (7)
   integer(ik), parameter :: SWPART_WEISS_NORMAN = 2_ik   !< Weiss-Norman (1985) band-specific split, eqs. (8)-(10)

   !----- Downwelling longwave source (docs/science/forcing.md §11). -------------------------!
   integer(ik), parameter :: LW_FILE       = 0_ik   !< read from file (ERA5-Land has strd)
   integer(ik), parameter :: LW_SYNTHESIZE = 1_ik   !< Brutsaert/Idso clear-sky synthesis (source lacking LW)
   !----- Clear-sky emissivity forms for that synthesis. They live HERE, beside the source         !
   !      selector they qualify, rather than with the kernel that evaluates them: meds_forcing     !
   !      links meds_config, so putting them in the kernel would point the config layer at the     !
   !      forcing layer and close a cycle. --------------------------------------------------------!
   integer(ik), parameter :: LW_CLEAR_BRUTSAERT = 0_ik   !< Brutsaert (1975), vapour-pressure power law
   integer(ik), parameter :: LW_CLEAR_IDSO      = 1_ik   !< Idso & Jackson (1969), temperature only

   !----- Model-start-before-first-record policy. ------------------------------------------!
   integer(ik), parameter :: CLAMP_ERROR = 0_ik   !< hard stop
   integer(ik), parameter :: CLAMP_HOLD  = 1_ik   !< clamp to the first record (w_next = 0)

   !----- Per-variable temporal-interpolation policy (used by interpolate_forcing). ---------!
   integer(ik), parameter :: INTERP_LINEAR = 0_ik   !< linear in the window (state vars)
   integer(ik), parameter :: INTERP_STEP   = 1_ik   !< step-constant (rainfall; hold prev)

   !----- How a polygon binds to the file's `grid` dimension (docs/science/forcing.md §1). ---!
   integer(ik), parameter :: GRIDMATCH_EXPLICIT = 0_ik  !< use grid_index verbatim (default)
   integer(ik), parameter :: GRIDMATCH_NEAREST  = 1_ik  !< pick the grid cell nearest [site] lat/lon (great-circle)

   !----- Where the free-atmosphere CO2 comes from ([forcing].co2_source;                         !
   !      docs/science/forcing.md §12). One source for every backend, looked up by MODEL time     !
   !      -- never the met file, so it does not repeat when the met is recycled. -----------------!
   integer(ik), parameter :: CO2_SOURCE_CONST = 0_ik   !< co2_const, held for the whole run
   integer(ik), parameter :: CO2_SOURCE_FILE  = 1_ik   !< a MEDS CO2 file (format: docs/science/forcing.md §12)

   !----- The forcing's own vertical frame ([forcing].height_above, .wind_exposure). The forcing   !
   !      is moved from its heights to the top of each patch's canopy air space (meds_lapse_rate;  !
   !      docs/science/forcing.md §8 Step 2).                                                      !
   integer(ik), parameter :: HEIGHT_ABOVE_ZERO_PLANE = 0_ik  !< heights above the patch's displacement height
                                                            !< (a reanalysis: its model has no d; CLM, JULES)
   integer(ik), parameter :: HEIGHT_ABOVE_GROUND     = 1_ik  !< heights above the ground (a flux tower)
   integer(ik), parameter :: WIND_EXPOSURE_OPEN_TERRAIN = 0_ik  !< the wind is an open-terrain diagnostic
                                                               !< (ERA5: from a blending height with z0 = 0.03 m)
   integer(ik), parameter :: WIND_EXPOSURE_LOCAL        = 1_ik  !< the wind was measured over this canopy

   !==========================================================================================!
   !  The [forcing]/[site] block. Plain scalars (no allocatables), so meds_config carries it  !
   !  trivially. There is NO gap_policy: MEDS never gap-fills, and a missing required value   !
   !  is a hard error (docs/science/forcing.md §10).                                          !
   !==========================================================================================!
   type :: forcing_config_t
      logical            :: forcing_on   = .false.               !< master gate (indep. of fast_biophysics_on)
      integer(ik)        :: format       = MET_BACKEND_ED_DEFAULT !< format: "ED_default" | "ED_ERA5land" | "const"
      character(len=MET_PATH_LEN) :: path = ''                   !< forcing NetCDF path (format = "ED_default")
      !----- The ED_ERA5land archive (format = "ED_ERA5land"; docs/science/forcing.md §1). The      !
      !      file template and static file are derived from data_path when left empty, which is     !
      !      what the archive's own layout needs; override them only for an archive laid out        !
      !      differently. Template tokens: {data_path}, {var}, {yyyy}, {mm}. -----------------------!
      character(len=MET_PATH_LEN) :: data_path     = ''          !< archive folder
      character(len=MET_PATH_LEN) :: file_template = ''          !< '' -> {data_path}/ED_ERA5land_{var}_{yyyy}{mm}.nc
      character(len=MET_PATH_LEN) :: static_file   = ''          !< '' -> {data_path}/ED_ERA5land_static.nc
      real(wp)           :: max_distance_km = 15.0_wp            !< [km] a site on a no-data cell takes the nearest
                                                                 !<      valid cell within this, else it is an error
      integer(ik)        :: grid_index   = 1_ik                  !< which (time,grid) location this polygon reads
      real(wp)           :: timestep     = 3600.0_wp             !< [s] native interval (hourly ERA5-Land)
      integer(ik)        :: avg_convention = METAVG_END          !< flux vars mean over the hour ENDING at the stamp
      integer(ik)        :: sw_partition = SWPART_CLEARIDX        !< ERA5-Land total SW -> partition required
      integer(ik)        :: lwdown_source = LW_FILE              !< file (ERA5-Land strd) | synthesize
      !----- Longwave SYNTHESIS parameters (docs/science/forcing.md eq. (17)), read only when       !
      !      lwdown_source = synthesize. `lw_clear_form` picks the clear-sky emissivity:            !
      !      Brutsaert (1975) uses the screen-level vapour pressure, Idso & Jackson (1969)          !
      !      temperature alone -- the fallback when a source carries humidity you do not trust.     !
      !      `lw_cloud_a` is the coefficient in the cloud term (1 + a(1-kt)); 0 gives a pure        !
      !      clear-sky sky, which UNDERESTIMATES under cloud. --------------------------------------!
      integer(ik)        :: lw_clear_form = 0_ik                  !< LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp)           :: lw_cloud_a    = 0.22_wp               !< [-] cloud-correction coefficient
      integer(ik)        :: co2_source   = CO2_SOURCE_CONST      !< "const" | "file"
      real(wp)           :: co2_const    = 420.0_wp              !< [umol/mol] co2_source = "const"
      character(len=MET_PATH_LEN) :: co2_file = ''               !< the MEDS CO2 file (co2_source = "file")
      !----- Recycling is OPT-IN (default off). It cannot be defaulted on: it REQUIRES a declared   !
      !      recycle_start/recycle_end below, and a default-constructed config has no meaningful     !
      !      window to offer. The TOML reader requires the key explicitly in any case.               !
      logical            :: recycle      = .false.               !< cycle the record when the run outruns the file
      !----- The RECYCLE WINDOW, DECLARED (never inferred; docs/science/forcing.md §9).             !
      !      recycle_start/recycle_end are required whenever recycle = .true. and are validated     !
      !      to span an exact whole number of calendar years [config check] and to land on the      !
      !      file's record stamps [met_open]; a mismatch is a HARD ERROR. A window that is not      !
      !      a whole number of years would slide the sub-daily phase and the season on every        !
      !      wrap while the daily mean stays right. The window is HALF-OPEN [recycle_start,         !
      !      recycle_end), recycle_end being the same instant one cycle later, and may start at     !
      !      any record stamp. A region's window starts at 00:00 or 01:00 on the 1st, because a     !
      !      region loads a month's forcing once (validate_config). --------------------------------!
      type(meds_time_t)  :: recycle_start                        !< first instant of the cycle (inclusive)
      type(meds_time_t)  :: recycle_end                          !< first instant AFTER the cycle (exclusive)
      integer(ik)        :: start_clamp  = CLAMP_ERROR           !< model start < base_time: error | hold the first record
      integer(ik)        :: grid_match   = GRIDMATCH_EXPLICIT    !< explicit grid_index | nearest [site] lat/lon
      !----- Site geolocation ([site]) -- solar geometry and lapse. Every forcing clock is UTC,     !
      !      so the longitude alone gives local solar time (docs/science/forcing.md eq. (2)). ------!
      real(wp)           :: latitude      = 42.44_wp             !< [deg +N] Ithaca NY
      real(wp)           :: longitude     = -76.50_wp            !< [deg +E] Ithaca NY (solar time)
      real(wp)           :: elevation     = 320.0_wp             !< [m] site elevation (Ithaca ~320 m; elevation lapse)
      !----- The forcing's own heights ([forcing]; defaults are ERA5-Land's). Every sample is moved  !
      !      from them to the top of each patch's canopy air space, per patch (meds_lapse_rate).     !
      real(wp)           :: tq_height            = 2.0_wp       !< [m] height of air temperature and humidity
      real(wp)           :: wind_height          = 10.0_wp      !< [m] height of the wind
      integer(ik)        :: height_above         = HEIGHT_ABOVE_ZERO_PLANE   !< what those heights are above
      integer(ik)        :: wind_exposure        = WIND_EXPOSURE_OPEN_TERRAIN !< open-terrain diagnostic | local
      real(wp)           :: wind_exposure_z0     = 0.03_wp      !< [m] roughness the open-terrain wind was made for
      real(wp)           :: wind_blending_height = 40.0_wp      !< [m] height it was brought down from
      !----- The terrain lapse, from the forcing cell's elevation to the site's ([site];            !
      !      docs/science/forcing.md eqs. (13)-(14)). The two keys below are read only when the     !
      !      lapse is on. --------------------------------------------------------------------------!
      logical            :: apply_elevation_lapse = .false.      !< lapse T, P, q, LW from the grid-cell elevation to the site
      real(wp)           :: lapse_rate_tair(12) = 0.0065_wp      !< [K/m] environmental lapse, January .. December
                                                                 !<       (positive = cooling upward)
      real(wp)           :: grid_elevation    = 320.0_wp         !< [m] elevation of the forcing source grid cell
   end type forcing_config_t

end module meds_forcing_config
