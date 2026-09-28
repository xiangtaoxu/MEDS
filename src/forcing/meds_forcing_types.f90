! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_forcing_types -- the RUNTIME meteorological-forcing types (design MEDS_FORCING_DESIGN.md   !
! sections 3.1-3.2; MEDS_POLYGON_RUNTIME_PLAN.md §10.3): the instantaneous per-site record         !
! met_forcing_t (a read-only boundary condition, the analogue of rad_forcing_t / chydro_forcing_t), !
! one raw file record met_record_t, the reader's shared state met_source_t (one per run) and each   !
! polygon's met_cursor_t, the archive's cell list met_cells_t and month buffer met_month_t, and the  !
! prescribed CO2 series co2_series_t.                                                               !
!                                                                                          !
! Lives in src/forcing (libmeds_forcing); links meds_shared only (meds_kinds, meds_time for the     !
! meds_time_t timestamp, meds_forcing_config for the [forcing]/[site] config). No physics library    !
! ever names these types; the reader (meds_met_driver) owns the only mutable forcing state.           !
!==========================================================================================!
module meds_forcing_types
   use meds_kinds,          only : wp, sp, ik
   use meds_time,           only : meds_time_t
   use meds_forcing_config, only : forcing_config_t, MET_BACKEND_CONST, MET_PATH_LEN
   implicit none
   private

   public :: met_forcing_t, met_record_t, met_source_t, met_cursor_t, met_cells_t, met_month_t
   public :: co2_series_t

   !==========================================================================================!
   !  met_forcing_t -- the instantaneous per-SITE atmospheric state the fast loop consumes.       !
   !  Per-site scalars + four SW-band reals (no cohort-sized allocatables) -> trivially copyable,   !
   !  GPU-mappable. All defaults are a valid "reference climate" box (the four SW streams SUM to      !
   !  400 W/m2 = today's fast_context_t%rad_sw_top), so met_forcing_t() is usable as the CONST         !
   !  backend value. cosz and rho_air are DERIVED each substep (not read).                             !
   !==========================================================================================!
   type :: met_forcing_t
      real(wp) :: tair_k       = 288.0_wp     !< [K]        air temperature at reference height
      real(wp) :: qair         = 0.008_wp     !< [kg/kg]    specific humidity
      real(wp) :: psurf_pa     = 101325.0_wp  !< [Pa]       surface pressure
      real(wp) :: rainf        = 0.0_wp       !< [kg/m2/s]  liquid precipitation rate (post phase-split)
      real(wp) :: snowfall        = 0.0_wp       !< [kg/m2/s]  frozen rainfall
      real(wp) :: wind         = 2.0_wp       !< [m/s]      wind speed at reference height
      !----- The wind VECTOR (§3.1, §5.3), carried beside the speed for a direction-aware consumer.  !
      !      Filled only when the source supplies components (has_wind_vector); aerodynamics reads   !
      !      `wind`, whose energy-form interpolation makes it at least the vector's length.          !
      real(wp) :: wind_u       = 2.0_wp       !< [m/s]      eastward component at reference height
      real(wp) :: wind_v       = 0.0_wp       !< [m/s]      northward component at reference height
      logical  :: has_wind_vector = .false.   !< .true. when the source supplied components, not just speed
      real(wp) :: lwdown       = 380.0_wp     !< [W/m2]     downwelling longwave (positive down)
      real(wp) :: par_beam     = 180.0_wp     !< [W/m2]     direct-beam PAR at canopy top
      real(wp) :: par_diffuse  = 40.0_wp      !< [W/m2]     diffuse PAR
      real(wp) :: nir_beam     = 150.0_wp     !< [W/m2]     direct-beam NIR
      real(wp) :: nir_diffuse  = 30.0_wp      !< [W/m2]     diffuse NIR   (Sigma = 400 W/m2)
      real(wp) :: co2          = 420.0_wp     !< [umol/mol] free-atmosphere CO2 ([forcing].co2_source)
      real(wp) :: cosz         = 0.0_wp       !< [-]        cos(solar zenith); DERIVED each substep
      real(wp) :: rho_air      = 1.2_wp       !< [kg/m3]    DERIVED from air_temp/psurf/qair
   contains
      procedure :: swdown         => met_swdown          !< total downward SW = sum of the four streams
      procedure :: rshort_diffuse => met_rshort_diffuse  !< diffuse SW = par_diffuse + nir_diffuse
   end type met_forcing_t

   !----- One raw timestamped record as read (pre-interpolation). The four SW streams are      !
   !      already split at INGEST (partition_shortwave) from the file's total SWdown. No CO2: it is  !
   !      not meteorology, and comes from co2_series_t or [forcing].co2_const on MODEL time (#184). !
   type :: met_record_t
      type(meds_time_t) :: when
      real(wp) :: tair_k = 288.0_wp, qair = 0.008_wp, psurf_pa = 101325.0_wp
      real(wp) :: rainf = 0.0_wp, wind = 2.0_wp, lwdown = 380.0_wp
      real(wp) :: wind_u = 2.0_wp, wind_v = 0.0_wp          !< wind vector (only when the source has it)
      real(wp) :: par_beam = 180.0_wp, par_diffuse = 40.0_wp, nir_beam = 150.0_wp, nir_diffuse = 30.0_wp
   end type met_record_t

   !==========================================================================================!
   !  met_cells_t -- the cells of a regular-grid archive (ED_ERA5land) that a run reads: one for  !
   !  a site, the valid cells of a box for the polygon runtime (§15.5). Grid indices are 0-based,  !
   !  as the netCDF C API counts. `by_chunk` lists the cells grouped by the archive's 16 x 16      !
   !  spatial chunk, `chunk_first` indexes into it (length nchunk+1), so a month load reads each    !
   !  touched chunk column once (§15.3), and only the box its cells occupy within that chunk.       !
   !==========================================================================================!
   type :: met_cells_t
      integer(ik) :: ncell = 0_ik
      integer(ik) :: nlat = 0_ik, nlon = 0_ik               !< archive grid size
      integer(ik), allocatable :: row(:), col(:)            !< 0-based (lat, lon) index of each cell
      real(wp),    allocatable :: lat(:), lon(:)            !< [deg] cell centre
      real(wp),    allocatable :: elevation(:)              !< [m] static-file orography
      real(wp),    allocatable :: land_fraction(:)          !< [-] static-file land fraction (box selection only)
      integer(ik) :: nchunk = 0_ik
      integer(ik), allocatable :: chunk_row(:), chunk_col(:)!< 0-based first row/col of the box read in each chunk
      integer(ik), allocatable :: chunk_nrow(:), chunk_ncol(:) !< that box's extent (1 x 1 for a site)
      integer(ik), allocatable :: chunk_first(:)            !< by_chunk(chunk_first(k):chunk_first(k+1)-1)
      integer(ik), allocatable :: by_chunk(:)               !< cell numbers grouped by chunk
   end type met_cells_t

   !----- One archive month for every cell of a domain: values(hour, cell, variable), float32 as  !
   !      the archive stores it (§15.3 memory table). Hours run 01:00 on the 1st .. 00:00 on the  !
   !      1st of the next month, the archive's end-stamped convention. ----------------------------!
   type :: met_month_t
      integer(ik) :: year = 0_ik, month = 0_ik, nt = 0_ik
      real(sp), allocatable :: values(:,:,:)
   end type met_month_t

   !==========================================================================================!
   !  co2_series_t -- a MEDS CO2 file, as read (co2_series_read). Each row is the mean over one    !
   !  period; the value sits at the period's middle, and the lookup interpolates linearly between  !
   !  middles on MODEL time, holding the end values over the first and last half-periods. Times   !
   !  are seconds after the first period's start, the series' own origin.                           !
   !==========================================================================================!
   type :: co2_series_t
      integer(ik)       :: n = 0_ik                         !< rows (periods)
      type(meds_time_t) :: first_start                      !< start of the first period: the time origin
      real(wp)          :: span_end_sec = 0.0_wp            !< [s] end of the last period, after first_start
      real(wp), allocatable :: mid_sec(:)                   !< [s] each period's middle, after first_start
      real(wp), allocatable :: co2(:)                       !< [umol/mol] each period's mean
   end type co2_series_t

   !==========================================================================================!
   !  met_source_t -- the SHARED reader state: one per run, read by every polygon (MEDS_POLYGON_     !
   !  RUNTIME_PLAN.md §10.3). Holds the config, the time axis and recycle window, the cells, and    !
   !  what is in memory: the archive's loaded month plus its carried record, or a MEDS forcing       !
   !  file's series. Only met_open, met_prefetch and met_close change it, all outside the compute    !
   !  phase, so stepping only reads it -- which is what lets polygons share it across threads.      !
   !==========================================================================================!
   type :: met_source_t
      type(forcing_config_t)  :: fcfg                       !< [forcing]/[site] config
      integer(ik) :: backend    = MET_BACKEND_CONST         !< NETCDF | ERA5LAND | CONST
      integer(ik) :: ncid       = -1_ik                     !< NetCDF handle (via meds_netcdf_c), -1 if closed
      integer(ik) :: grid_index = 1_ik                      !< MEDS forcing file: the location (1..ngrid) read
      integer(ik) :: ngrid      = 1_ik                      !< total locations in the file (cells for the archive)
      integer(ik) :: nrec       = 0_ik                      !< total time records on the axis
      integer(ik) :: rec_first  = 1_ik                      !< the record a cursor's first bracket starts at
      real(wp)    :: dt_forcing  = 3600.0_wp                !< [s] native forcing interval
      type(meds_time_t)  :: base_time                       !< time-coordinate anchor ("seconds since base_time")
      real(wp), allocatable :: time_sec(:)                  !< [s] cached time coordinate (seconds since base_time)
      !----- CALENDAR recycling over the window DECLARED in [forcing] (recycle_start/recycle_end)  !
      !      and validated against this file at open (validate_recycle_window). Nothing here is     !
      !      sniffed from the file: the anchor is the config's recycle_start, the cycle length is    !
      !      the config's whole-year span, and a config that disagrees with the file's actual record !
      !      timestamps is a hard error rather than a silent fallback. -------------------------------!
      type(meds_time_t) :: cycle_anchor                     !< = fcfg%recycle_start; the cycle's first instant
      integer(ik) :: n_cycle_years   = 0_ik                 !< whole calendar years the declared window spans
      integer(ik) :: irec_cycle_first = 0_ik                !< record index of cycle_anchor (exact match, 1-based)
      integer(ik) :: irec_cycle_last  = 0_ik                !< record index of the LAST record strictly inside the window
      logical     :: has_wind_vector = .false.              !< the source supplies u10/v10 (§15.4)
      !----- ED_ERA5land archive backend (§15). The months the run needs form one continuous hourly   !
      !      axis (time_sec, seconds since 1970-01-01, the archive's own epoch); month k holds records !
      !      month_rec0(k)+1 .. month_rec0(k)+hours, and one month of every cell sits in `buffer`.    !
      character(len=MET_PATH_LEN) :: file_template = ''     !< resolved file template
      character(len=MET_PATH_LEN) :: static_file   = ''     !< resolved static file
      type(met_cells_t) :: cells                            !< the cells read: one for a site, many for a region
      integer(ik), allocatable :: month_year(:), month_month(:), month_rec0(:)
      type(met_month_t) :: buffer                           !< the loaded month (year = 0 before the first)
      !----- Everything a step reads is in memory before the step (MEDS_POLYGON_RUNTIME_PLAN.md §4, R1): !
      !      the archive's loaded month plus the one record read before it (00:00 on the 1st, from the   !
      !      previous month's file, or the window's last record at the recycle wrap), loaded by         !
      !      met_prefetch, and the recycle window's first day, kept from open; a MEDS forcing file's    !
      !      records for the run, read at open. --------------------------------------------------------!
      real(sp),    allocatable :: carry(:,:)                !< (cell, variable) the record before the loaded month
      integer(ik) :: carry_rec = 0_ik                       !< its record index on the axis (0 = none)
      !----- The window's first record through the next midnight (recycling only), from axis month 1. !
      !      The daily step that crosses the seam reads the window's last record and then these, and    !
      !      the two lie in different months, so the one-month buffer cannot hold both. ---------------!
      real(sp),    allocatable :: head(:,:,:)               !< (record, cell, variable) from irec_cycle_first on
      integer(ik) :: n_head    = 0_ik                       !< records in head (0 = none)
      integer(ik) :: n_loads   = 0_ik                       !< archive month loads so far (tests: none inside a step)
      real(wp),    allocatable :: series(:,:)               !< (record, field) MEDS forcing file at grid_index
      character(len=24), allocatable :: series_name(:)      !< the fields present in that file
      !----- The prescribed CO2 (co2_source = "file"), read at open. One global series for every   !
      !      polygon, looked up on model time, so it is shared read-only like the rest. ------------!
      type(co2_series_t) :: co2
   end type met_source_t

   !==========================================================================================!
   !  met_cursor_t -- ONE polygon's reader state: its cell and location, the two records bracketing !
   !  the model time, and the cursor. Plain scalars and records, so a region holds one per polygon   !
   !  cheaply. The location drives solar geometry and the elevation lapse: a site run's comes from   !
   !  [site], a region polygon's from its cell (the ED2 polygon keeps its lat/lon the same way).     !
   !==========================================================================================!
   type :: met_cursor_t
      integer(ik) :: cell = 1_ik                            !< this polygon's cell in the source (1 for a site)
      real(wp)    :: latitude_deg     = 42.44_wp            !< [deg +N] solar geometry
      real(wp)    :: longitude_deg    = -76.50_wp           !< [deg +E] local solar time
      real(wp)    :: utc_offset_h     = 0.0_wp              !< [h] of the forcing clock
      real(wp)    :: elevation_m      = 320.0_wp            !< [m] site elevation (elevation lapse target)
      real(wp)    :: grid_elevation_m = 320.0_wp            !< [m] forcing-cell elevation (lapse origin)
      integer(ik) :: irec_prev  = 0_ik                      !< cursor: index of the "previous" bracketing record
      type(met_record_t) :: rec_prev, rec_next              !< the two records bracketing the model time
      logical     :: at_wrap_seam    = .false.              !< current bracket is the cycle-boundary seam (last -> first)
      !----- LONGWAVE SYNTHESIS (#182): the last DAYTIME clearness index, held through the night.   !
      !      kt is undefined after dark -- there is no shortwave to divide -- and treating that as   !
      !      kt = 0 would apply the maximum cloud correction to every night. Holding dusk's value    !
      !      is the better guess, and the one scalar of state it needs lives here rather than in a   !
      !      module variable. Seeded to 1 (clear) so a run that starts at night starts clear-sky.    !
      real(wp)    :: kt_last_day    = 1.0_wp                !< [-] last daytime clearness index
   end type met_cursor_t

contains

   pure real(wp) function met_swdown(met) result(sw)
      class(met_forcing_t), intent(in) :: met
      sw = met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse
   end function met_swdown

   pure real(wp) function met_rshort_diffuse(met) result(diff)
      class(met_forcing_t), intent(in) :: met
      diff = met%par_diffuse + met%nir_diffuse
   end function met_rshort_diffuse

end module meds_forcing_types
