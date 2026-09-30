! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_forcing_config -- the [forcing]/[site] configuration for meteorological forcing, plus  !
! the selector codes shared by the config, the reader, and the disaggregation kernels.         !
!                                                                                          !
! Placed in src/config (NOT src/forcing) so meds_config -- the DAG ROOT -- can carry            !
! forcing_config_t as a plain-scalar component with NO backward `config -> forcing` edge (the    !
! same rule the energy design used for soil_thermal_params_t). Pure DATA + parameters; uses       !
! meds_kinds and meds_time only. Defaults are the Ithaca NY / ERA5-Land reference site (design     !
! MEDS_FORCING_DESIGN.md sections 3.3, 6.6).                                                        !
!==========================================================================================!
module meds_forcing_config
   use meds_kinds, only : wp, ik
   use meds_time,  only : meds_time_t
   implicit none
   private

   public :: forcing_config_t
   public :: MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT, MET_BACKEND_ED_ERA5LAND
   public :: MET_PATH_LEN
   public :: METAVG_INSTANT, METAVG_END, METAVG_BEGIN, METAVG_CENTER
   public :: SWPART_PASSTHROUGH, SWPART_WEISS_NORMAN, SWPART_CLEARIDX
   public :: LW_FILE, LW_SYNTHESIZE, LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO
   public :: CLAMP_ERROR, CLAMP_HOLD
   public :: INTERP_LINEAR, INTERP_STEP, INTERP_COSZ
   public :: GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   public :: CO2_SOURCE_CONST, CO2_SOURCE_FILE
   public :: HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND, WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL

   !----- Reader backend ([forcing].format): the single MEDS forcing file, the global ED_ERA5land  !
   !      archive, or a no-file reference-climate box. ------------------------------------------!
   integer(ik), parameter :: MET_BACKEND_CONST       = 0_ik  !< no file: met_forcing_t defaults (reference climate)
   integer(ik), parameter :: MET_BACKEND_ED_DEFAULT  = 1_ik  !< the multi-grid forcing NetCDF (format = "ED_default")
   integer(ik), parameter :: MET_BACKEND_ED_ERA5LAND = 2_ik  !< the monthly ED_ERA5land archive (format = "ED_ERA5land")

   integer, parameter :: MET_PATH_LEN = 1024                !< length of every forcing path field (§15.2)

   !----- Timestamp semantics of a forcing record (avg_convention). ------------------------!
   integer(ik), parameter :: METAVG_INSTANT = 0_ik   !< value is instantaneous AT the stamp
   integer(ik), parameter :: METAVG_END     = 1_ik   !< flux mean over the interval ENDING at the stamp (ERA5-Land)
   integer(ik), parameter :: METAVG_BEGIN   = 2_ik   !< mean over the interval BEGINNING at the stamp
   integer(ik), parameter :: METAVG_CENTER  = 3_ik   !< mean over the interval CENTERED on the stamp

   !----- Shortwave partition of total SWdown into the four (beam/diffuse)x(PAR/NIR) streams. !
   integer(ik), parameter :: SWPART_PASSTHROUGH  = 0_ik   !< file already carries the four streams
   integer(ik), parameter :: SWPART_CLEARIDX     = 1_ik   !< clearness-index (Erbs) split -- the ERA5-Land P0 default
   integer(ik), parameter :: SWPART_WEISS_NORMAN = 2_ik   !< Weiss-Norman band-specific (P1)

   !----- Downwelling longwave source. -----------------------------------------------------!
   integer(ik), parameter :: LW_FILE       = 0_ik   !< read from file (ERA5-Land has strd)
   integer(ik), parameter :: LW_SYNTHESIZE = 1_ik   !< Brutsaert/Idso clear-sky synthesis (source lacking LW)
   !----- Clear-sky emissivity forms for that synthesis (#182). They live HERE, beside the source  !
   !      selector they qualify, rather than with the kernel that evaluates them: meds_forcing      !
   !      links meds_config, so putting them in the kernel would point the config layer at the      !
   !      forcing layer and close a cycle.  ----------------------------------------------------------!
   integer(ik), parameter :: LW_CLEAR_BRUTSAERT = 0_ik   !< Brutsaert (1975), vapour-pressure power law
   integer(ik), parameter :: LW_CLEAR_IDSO      = 1_ik   !< Idso & Jackson (1969), temperature only

   !----- Model-start-before-first-record policy. ------------------------------------------!
   integer(ik), parameter :: CLAMP_ERROR = 0_ik   !< hard stop
   integer(ik), parameter :: CLAMP_HOLD  = 1_ik   !< clamp to record #1 (w_next = 0)

   !----- Per-variable temporal-interpolation policy (used by interpolate_forcing). ---------!
   integer(ik), parameter :: INTERP_LINEAR = 0_ik   !< linear in the window (state vars)
   integer(ik), parameter :: INTERP_STEP   = 1_ik   !< step-constant (rainfall; hold prev)
   integer(ik), parameter :: INTERP_COSZ   = 2_ik   !< cosz-weighted (shortwave; handled by disaggregate_shortwave)

   !----- How a polygon binds to the file's `grid` dimension (multi-polygon P2 subset). -----!
   integer(ik), parameter :: GRIDMATCH_EXPLICIT = 0_ik  !< use grid_index verbatim (default; current behaviour)
   integer(ik), parameter :: GRIDMATCH_NEAREST  = 1_ik  !< pick the grid cell nearest [site] lat/lon (great-circle)

   !----- Where the free-atmosphere CO2 comes from ([forcing].co2_source, #184). One source for   !
   !      every backend, looked up by MODEL time -- never the met file, so it does not repeat when  !
   !      the met is recycled. ------------------------------------------------------------------!
   integer(ik), parameter :: CO2_SOURCE_CONST = 0_ik   !< co2_const, held for the whole run
   integer(ik), parameter :: CO2_SOURCE_FILE  = 1_ik   !< a MEDS CO2 file (format: docs/science/forcing.md, "CO2")

   !----- The forcing's own vertical frame ([forcing].height_above, .wind_exposure). The forcing  !
   !      is moved from its heights to the top of each patch's canopy air space (meds_lapse_rate). !
   integer(ik), parameter :: HEIGHT_ABOVE_ZERO_PLANE = 0_ik  !< heights above the patch's displacement height
                                                            !< (a reanalysis: its model has no d; CLM, JULES)
   integer(ik), parameter :: HEIGHT_ABOVE_GROUND     = 1_ik  !< heights above the ground (a flux tower)
   integer(ik), parameter :: WIND_EXPOSURE_OPEN_TERRAIN = 0_ik  !< the wind is an open-terrain diagnostic
                                                               !< (ERA5: from a blending height with z0 = 0.03 m)
   integer(ik), parameter :: WIND_EXPOSURE_LOCAL        = 1_ik  !< the wind was measured over this canopy

   !==========================================================================================!
   !  The [forcing]/[site] block. Plain scalars (no allocatables), so meds_config carries it     !
   !  trivially. NOTE: there is NO gap_policy -- MEDS never gap-fills; a missing required value    !
   !  is a hard error (design comment 2 / §5.5).                                                   !
   !==========================================================================================!
   type :: forcing_config_t
      logical            :: forcing_on   = .false.               !< master gate (indep. of fast_biophysics_on)
      integer(ik)        :: backend      = MET_BACKEND_ED_DEFAULT !< format: "ED_default" | "ED_ERA5land" | "const"
      character(len=MET_PATH_LEN) :: path = ''                   !< forcing NetCDF path (format = "ED_default")
      !----- The ED_ERA5land archive (format = "ED_ERA5land", MEDS_FORCING_DESIGN.md §15.2). The file !
      !      template and static file are derived from data_path when left empty, which is what the  !
      !      archive's own layout needs; override them only for an archive laid out differently.    !
      !      Template tokens: {data_path}, {var}, {yyyy}, {mm}.                                        !
      character(len=MET_PATH_LEN) :: data_path     = ''          !< archive folder
      character(len=MET_PATH_LEN) :: file_template = ''          !< '' -> {data_path}/ED_ERA5land_{var}_{yyyy}{mm}.nc
      character(len=MET_PATH_LEN) :: static_file   = ''          !< '' -> {data_path}/ED_ERA5land_static.nc
      real(wp)           :: max_distance_km = 15.0_wp            !< [km] a site on a no-data cell takes the nearest
                                                                 !<      valid cell within this, else it is an error
      integer(ik)        :: grid_index   = 1_ik                  !< which (time,grid) location this polygon reads
      real(wp)           :: dt_forcing   = 3600.0_wp             !< [s] native interval (hourly ERA5-Land)
      integer(ik)        :: avg_convention = METAVG_END          !< flux vars mean over the hour ENDING at the stamp
      integer(ik)        :: sw_partition = SWPART_CLEARIDX        !< ERA5-Land total SW -> partition required
      integer(ik)        :: lwdown_source = LW_FILE              !< file (ERA5-Land strd) | synthesize
      !----- Longwave SYNTHESIS parameters (#182), read only when lwdown_source = synthesize.        !
      !      `lw_clear_form` picks the clear-sky emissivity: Brutsaert (1975) uses the screen-level   !
      !      vapour pressure, Idso & Jackson (1969) temperature alone -- the fallback when a source   !
      !      carries humidity you do not trust. `lw_cloud_a` is the coefficient in the cloud term     !
      !      (1 + a(1-kt)); 0 gives a pure clear-sky sky, which UNDERESTIMATES under cloud.           !
      integer(ik)        :: lw_clear_form = 0_ik                  !< LW_CLEAR_BRUTSAERT | LW_CLEAR_IDSO
      real(wp)           :: lw_cloud_a    = 0.22_wp               !< [-] cloud-correction coefficient
      integer(ik)        :: co2_source   = CO2_SOURCE_CONST      !< "const" | "file"
      real(wp)           :: co2_const    = 420.0_wp              !< [umol/mol] co2_source = "const"
      character(len=MET_PATH_LEN) :: co2_file = ''               !< the MEDS CO2 file (co2_source = "file")
      !----- Recycling is OPT-IN (default off). It cannot be defaulted on: it REQUIRES a declared   !
      !      recycle_start/recycle_end below, and a default-constructed config has no meaningful     !
      !      window to offer. The TOML reader requires the key explicitly in any case.               !
      logical            :: recycle      = .false.               !< cycle the record when the run outruns the file
      !----- The RECYCLE WINDOW, DECLARED (never inferred). MEDS does not sniff the file to guess   !
      !      where a cycle starts or how long it is: recycle_start/recycle_end are required whenever  !
      !      recycle=.true., and are validated to (a) span an exact whole number of calendar years     !
      !      [config check] and (b) match the file's actual record timestamps [met_open]. A mismatch    !
      !      -- e.g. a config declaring 00:00:00 against an ERA5-Land file whose records are stamped     !
      !      01:00:00 -- is a HARD ERROR, not a silent fallback. Inferring the window from the file      !
      !      instead drives a 30-yr run with phase-scrambled sub-daily shortwave: the ERA5-Land span is    !
      !      366 d 22 h, so every wrap on it shifts hour-of-day, while the daily MEAN stays right and the  !
      !      slow demography still looks sane. See MEDS_FORCING_DESIGN.md.                                 !
      !      The window is HALF-OPEN [recycle_start, recycle_end): recycle_end is the exclusive upper      !
      !      bound, i.e. the same instant one cycle later, so N years of records are covered exactly once. !
      !      The anchor may sit at any record stamp (mid-year windows are fine) -- only the whole-year     !
      !      SPAN is required, not a Jan-1 start. A region's window starts at 00:00 or 01:00 on the 1st,   !
      !      because a region loads a month's forcing once (validate_config).                             !
      type(meds_time_t)  :: recycle_start                        !< first instant of the cycle (inclusive)
      type(meds_time_t)  :: recycle_end                          !< first instant AFTER the cycle (exclusive)
      integer(ik)        :: start_clamp  = CLAMP_ERROR           !< model start < base_time: error | hold record #1
      integer(ik)        :: grid_match   = GRIDMATCH_EXPLICIT    !< explicit grid_index | nearest [site] lat/lon (§4.1)
      !----- site geolocation ([site]) -- solar geometry + lapse. Every forcing clock is UTC       !
      !      (MEDS_FLUX_TOWER_FORCING_PLAN.md D1), so the longitude alone gives local solar time. ----!
      real(wp)           :: latitude_deg  = 42.44_wp             !< [deg +N] Ithaca NY
      real(wp)           :: longitude_deg = -76.50_wp            !< [deg +E] Ithaca NY (solar time, §5.1)
      real(wp)           :: elevation_m   = 320.0_wp             !< [m] site elevation (Ithaca ~320 m; elevation lapse)
      !----- The forcing's own heights ([forcing]; defaults are ERA5-Land's). Every sample is moved  !
      !      from them to the top of each patch's canopy air space, per patch (meds_lapse_rate).     !
      real(wp)           :: tq_height            = 2.0_wp       !< [m] height of air temperature and humidity
      real(wp)           :: wind_height          = 10.0_wp      !< [m] height of the wind
      integer(ik)        :: height_above         = HEIGHT_ABOVE_ZERO_PLANE   !< what those heights are above
      integer(ik)        :: wind_exposure        = WIND_EXPOSURE_OPEN_TERRAIN !< open-terrain diagnostic | local
      real(wp)           :: wind_exposure_z0     = 0.03_wp      !< [m] roughness the open-terrain wind was made for
      real(wp)           :: wind_blending_height = 40.0_wp      !< [m] height it was brought down from
      !----- The terrain lapse, from the forcing cell's elevation to the site's ([site]). The two  !
      !      keys below are read only when the lapse is on. --------------------------------------!
      logical            :: apply_elevation_lapse = .false.      !< lapse T, P, q, LW from the grid-cell elevation to the site
      real(wp)           :: lapse_rate_tair(12) = 0.0065_wp      !< [K/m] environmental lapse, January .. December
                                                                 !<       (positive = cooling upward)
      real(wp)           :: grid_elevation_m  = 320.0_wp         !< [m] elevation of the forcing source grid cell
   end type forcing_config_t

end module meds_forcing_config
