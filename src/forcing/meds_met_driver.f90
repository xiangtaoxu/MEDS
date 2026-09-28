! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_driver -- the meteorological-forcing READER (design MEDS_FORCING_DESIGN.md sections !
! 4 and 15). Opens the multi-grid forcing NetCDF (ED_default), the global ED_ERA5land archive, or the !
! no-file CONST backend into a met_source_t shared by every polygon of a run; each polygon has a   !
! met_cursor_t holding its cell, its location and the two records bracketing the model time,      !
! slides that window as the model marches, and produces an instantaneous met_forcing_t via the     !
! pure disaggregation kernels (MEDS_POLYGON_RUNTIME_PLAN.md §10.3). Owns file I/O (via meds_netcdf_c !
! and meds_era5land_reader) -- the ED2 cgrid%metinput analogue, threaded by the driver, never a      !
! global. Stepping (met_advance, met_instant) only reads the source.                                  !
!                                                                                          !
! The archive's months form one continuous hourly axis over the months a run needs (the recycle   !
! window, or the run period), so bracketing, recycling and the seam are the same code for both      !
! file backends; only where a record's values come from differs (read_record).                      !
!                                                                                          !
! MEDS NEVER gap-fills (design §5.5): a missing/NaN required value is a HARD ERROR. netCDF is a       !
! hard dependency, so the reader is always the real thing (no stub).                                   !
!                                                                                          !
! CO2 is not meteorology (#184). It never comes from the met file: [forcing].co2_source gives one      !
! value (co2_const) or a MEDS CO2 file (meds_co2_series), read at open, looked up on MODEL time for     !
! every backend alike -- so it does not repeat when the met is recycled.                               !
!==========================================================================================!
module meds_met_driver
   use iso_c_binding,       only : c_int, c_size_t, c_double
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : tiny_num
   use meds_therm_lib,         only : air_density
   use meds_time,           only : meds_time_t, time_from_string, time_advance_seconds,        &
                                   seconds_between, seconds_into_day, time_lt, time_le,        &
                                   is_leap_year, days_in_year, time_to_string,                 &
                                   time_advance_years, whole_years_between
   use meds_forcing_config, only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT,  &
                                   MET_BACKEND_ED_ERA5LAND, MET_PATH_LEN,                       &
                                   METAVG_END, METAVG_BEGIN, SWPART_PASSTHROUGH,                &
                                   CLAMP_ERROR, INTERP_LINEAR,                                  &
                                   GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST, LW_SYNTHESIZE,        &
                                   CO2_SOURCE_FILE, HEIGHT_ABOVE_GROUND
   use meds_forcing_types,  only : met_forcing_t, met_record_t, met_source_t, met_cursor_t, met_cells_t, &
                                   HUMIDITY_QAIR, HUMIDITY_RHAIR, HUMIDITY_TDEW
   use meds_config,         only : MAX_RECYCLE_YEARS   ! the config's bound: one definition
   use meds_lapse_rate,     only : lapse_air_temperature, lapse_pressure, monthly_lapse_rate,    &
                                   lapse_specific_humidity, lapse_longwave
   use meds_co2_series,     only : co2_series_read, co2_series_at, co2_series_covers,          &
                                   co2_series_end, co2_series_free
   use meds_forcing_kernels, only : interpolate_forcing, interpolate_wind_energy,              &
                                   met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,   &
                                   partition_shortwave, precip_phase, nearest_grid_index,       &
                                   great_circle_distance,                                      &
                                   clearness_index, synthesize_lwdown, dewpoint_to_specific_humidity, &
                                   rh_to_specific_humidity
   use meds_era5land_reader, only : era5land_path, era5land_default_template, era5land_default_static, &
                                   era5land_select_site, era5land_load_month, era5land_month_hours, &
                                   ERA_NVAR, ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10,        &
                                   ERA_RAINF, ERA_SWDOWN, ERA_LWDOWN, ERA_VAR_NAME, ERA_EPOCH,       &
                                   ERA_OK, ERA_ERR_OPEN, ERA_ERR_NO_CELL
   use meds_netcdf_c,       only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f,                  &
                                   nc_get_att_text_f, nc_get_att_double_f, nc_get_vara_double,  &
                                   nc_close, nc_check,                                          &
                                   NC_NOERR, NC_NOWRITE, NC_GLOBAL
   implicit none
   private

   public :: met_open, met_cursor_init, met_advance, met_instant, met_close, met_prefetch
   public :: MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD,                &
             MET_ERR_WINDOW_NOT_COVERED, MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,         &
             MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED,  &
             MET_ERR_CO2_IN_MET_FILE, MET_ERR_NOT_UTC, MET_ERR_HUMIDITY

   real(wp), parameter :: U_MIN     = 0.1_wp     !< [m/s] wind floor (M-O similarity stability)
   integer(ik), parameter :: N_COSZ_SUB = 10_ik  !< sub-samples per forcing interval for <cosz>_win

   !----- met_open status codes (see the `stat` argument). ---------------------------------!
   integer(ik), parameter :: MET_OK                          = 0_ik
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_WHOLE_YEARS  = 1_ik   !< recycle_end - recycle_start /= N years
   integer(ik), parameter :: MET_ERR_START_NOT_A_RECORD      = 2_ik   !< recycle_start is not a record stamp
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_COVERED      = 3_ik   !< file stops short of recycle_end
   integer(ik), parameter :: MET_ERR_DT_MISMATCH             = 4_ik   !< dt_forcing /= the file's record spacing
   integer(ik), parameter :: MET_ERR_AXIS_NOT_UNIFORM        = 5_ik   !< the file's time axis is ragged
   integer(ik), parameter :: MET_ERR_ATTR_MISMATCH           = 6_ik   !< a file global attribute contradicts the config
   integer(ik), parameter :: MET_ERR_ARCHIVE                 = 7_ik   !< archive: static file, site cell or a month missing
   integer(ik), parameter :: MET_ERR_CO2_FILE                = 8_ik   !< the CO2 file is unreadable or breaks format 1
   integer(ik), parameter :: MET_ERR_CO2_NOT_COVERED         = 9_ik   !< the CO2 file does not cover the run
   integer(ik), parameter :: MET_ERR_CO2_IN_MET_FILE         = 10_ik  !< the met file carries CO2air
   integer(ik), parameter :: MET_ERR_NOT_UTC                 = 11_ik  !< the file does not say time_zone = "UTC"
   integer(ik), parameter :: MET_ERR_HUMIDITY                = 12_ik  !< no, two, or out-of-range humidity variables

   !----- Upper bound on the declared recycle window, in whole calendar years (search bound only). !
   !----- Tolerance for "this record stamp IS that instant" [s]. The time axis is float seconds,   !
   !      so an exact == would be brittle; sub-second slack is far below any real forcing dt. -------!
   real(wp), parameter :: REC_MATCH_TOL = 0.5_wp
   !----- A file's declared measurement height agrees with the config within this [m]. -----------!
   real(wp), parameter :: HEIGHT_MATCH_TOL = 0.01_wp
   !----- RHair is a fraction; a value above this is a percentage written into a fractional field. !
   real(wp), parameter :: RH_FRACTION_MAX = 1.5_wp

contains

   !=======================================================================================!
   !  OPEN: CONST -> reference climate; ED_DEFAULT -> read the grid/time dims, the time axis, the !
   !  base-time anchor from the `time:units` attribute, and load records #1-2 at grid_index.     !
   !=======================================================================================!
   !  `stat` (optional) reports a rejected recycle window instead of halting, so the validation    !
   !  is exercisable from a test (CLAUDE.md: "errors via error stop / status codes ... so failures  !
   !  are catchable in tests"). Absent -> a rejection is a hard error, which is what a production    !
   !  run wants: a forcing file that does not match its declared window must never be guessed at.   !
   !  `run_start`/`run_end` bound the months a non-recycling ED_ERA5land run reads; a recycling    !
   !  run reads the declared window instead, and the other backends ignore them.                  !
   subroutine met_open(src, fcfg, stat, run_start, run_end, cells)
      type(met_source_t),     intent(inout)  :: src
      type(forcing_config_t), intent(in)    :: fcfg
      integer(ik), optional,  intent(out)   :: stat
      type(meds_time_t), optional, intent(in) :: run_start, run_end
      type(met_cells_t), optional, intent(in) :: cells   !< a region's cells (archive only); absent: the site's
      integer(c_int)    :: st, ncid
      integer(c_size_t) :: dlen
      character(len=256):: units
      logical           :: ok
      integer(ik)       :: vstat

      if (present(stat)) stat = MET_OK

      src%fcfg       = fcfg
      src%backend    = fcfg%backend
      src%grid_index = fcfg%grid_index
      src%dt_forcing = fcfg%dt_forcing
      src%has_wind_vector = .false.

      !----- CO2 first, before the backends branch, so every backend takes it the same way. ------!
      call open_co2(src, run_start, run_end, vstat)
      if (vstat /= MET_OK) then
         call met_close(src)
         if (present(stat)) then ; stat = vstat ; return ; end if
         error stop 'met_open: the CO2 file cannot drive this run (see the message above)'
      end if

      if (fcfg%backend == MET_BACKEND_CONST) then
         src%ngrid = 1_ik ; src%nrec = 0_ik
         return
      end if

      if (fcfg%backend == MET_BACKEND_ED_ERA5LAND) then
         call open_archive(src, run_start, run_end, vstat, cells)
         if (vstat == MET_OK) call validate_recycle_window(src, vstat)
         if (vstat /= MET_OK) then
            call met_close(src)
            if (present(stat)) then ; stat = vstat ; return ; end if
            error stop 'met_open: the ED_ERA5land archive cannot drive this run (see the message above)'
         end if
         call load_axis_month(src, 1_ik)
         call keep_window_head(src)
         src%rec_first = 1_ik
         return
      end if

      !----- NetCDF backend. -----------------------------------------------------------------!
      st = nc_open_f(trim(fcfg%path), NC_NOWRITE, ncid)
      call nc_check(st, 'met_open: nc_open '//trim(fcfg%path))
      src%ncid = int(ncid, ik)

      st = nc_inq_dimlen_f(ncid, 'grid', dlen) ; call nc_check(st, 'met_open: grid dim')
      src%ngrid = int(dlen, ik)
      !----- bind this polygon to a grid slice: explicit index, or nearest [site] lat/lon (§4.1). !
      if (fcfg%grid_match == GRIDMATCH_NEAREST) then
         call resolve_grid_index(src)                        ! sets src%grid_index in 1..ngrid
      else                                                   ! GRIDMATCH_EXPLICIT
         if (fcfg%grid_index < 1_ik .or. fcfg%grid_index > src%ngrid) then
            write(*,'(a,i0,a,i0)') 'met_open: grid_index ', fcfg%grid_index, ' out of range 1..', src%ngrid
            error stop 'met_open: grid_index out of range'
         end if
      end if

      st = nc_inq_dimlen_f(ncid, 'time', dlen) ; call nc_check(st, 'met_open: time dim')
      src%nrec = int(dlen, ik)
      if (src%nrec < 1_ik) error stop 'met_open: forcing file has no time records'

      !----- base time from "seconds since <base>" units attribute. --------------------------!
      block
         integer(c_int) :: vid
         st = nc_inq_varid_f(ncid, 'time', vid) ; call nc_check(st, 'met_open: time varid')
         st = nc_get_att_text_f(ncid, vid, 'units', units) ; call nc_check(st, 'met_open: time:units')
         call parse_time_units(units, src%base_time, ok)
         if (.not. ok) error stop 'met_open: could not parse time:units "seconds since <base>"'
      end block

      !----- cache the whole time axis (seconds since base_time). ----------------------------!
      call read_time_axis(src)

      !----- The wind VECTOR (§7.1, §15.4): a file carrying u10 and v10 supplies it, and the speed  !
      !      is then derived from the components; a file with Wind alone gives the speed only.     !
      block
         integer(c_int) :: vu, vv
         src%has_wind_vector = nc_inq_varid_f(ncid, 'u10', vu) == NC_NOERR .and.                  &
                               nc_inq_varid_f(ncid, 'v10', vv) == NC_NOERR
      end block

      !----- The humidity the file carries: exactly one of RHair, Tdew and Qair. ------------------!
      call detect_humidity(src, ncid, vstat)
      if (vstat /= MET_OK) then
         call met_close(src)
         if (present(stat)) then ; stat = vstat ; return ; end if
         error stop 'met_open: the forcing file''s humidity is ambiguous or absent (see the message above)'
      end if

      !----- V4 (#185): the file's own record spacing against [forcing].dt_forcing, and the two   !
      !      global attributes the prep script writes against the config that claims to describe    !
      !      the same file. A file that disagreed with its config would otherwise be silently       !
      !      mis-timed or mis-partitioned. ------------------------------------------------------------!
      call validate_file_against_config(src, ncid, vstat)
      if (vstat /= MET_OK) then
         call met_close(src)
         if (present(stat)) then ; stat = vstat ; return ; end if
         error stop 'met_open: the forcing file contradicts [forcing] (see the message above)'
      end if

      !----- V2/V3: check the DECLARED recycle window against this file's actual record stamps.  !
      !      Cannot live in the config sanity check (V1, whole-year span) -- it needs the time      !
      !      axis, which is only known here. -------------------------------------------------------!
      call validate_recycle_window(src, vstat)
      if (vstat /= MET_OK) then
         !----- Release the file before bailing out. A `stat` return is a normal (if unhappy) exit  !
         !      the caller may retry from, and an un-closed netCDF handle keeps an HDF5 lock on the  !
         !      path -- a later create/open of the same file then fails with a bare "permission      !
         !      denied" that points nowhere near the real cause. ------------------------------------!
         call met_close(src)
         if (present(stat)) then ; stat = vstat ; return ; end if
         error stop 'met_open: forcing file does not match the declared [forcing] recycle window'
      end if

      !----- The records the run can read, into memory, so no step reads the file (R1): the recycle !
      !      window, or the run period (the whole file when neither is known, as in unit tests). Only  !
      !      those: a file written in 1 x 1 chunks makes a whole-file read cost HDF5 metadata for      !
      !      every chunk. ----------------------------------------------------------------------------!
      block
         integer(ik) :: r0, r1
         r0 = 1_ik ; r1 = src%nrec
         if (src%fcfg%recycle .and. src%n_cycle_years >= 1_ik) then
            r0 = max(1_ik, src%irec_cycle_first - 1_ik) ; r1 = min(src%nrec, src%irec_cycle_last + 1_ik)
         else if (present(run_start) .and. present(run_end)) then
            r0 = max(1_ik, first_record_after(src, seconds_between(src%base_time, run_start)) - 2_ik)
            r1 = min(src%nrec, first_record_after(src, seconds_between(src%base_time, run_end)) + 1_ik)
         end if
         if (r1 <= r0) r1 = min(src%nrec, r0 + 1_ik)
         call read_series(src, ncid, r0, r1)
         st = nc_close(ncid) ; src%ncid = -1_ik
         src%rec_first = r0                          ! a cursor's first bracket: the first records in range
      end block

      !----- RHair is a fraction. Written as a percentage it would clip every record to saturation, !
      !      a silent 100 % sky, so a value no fraction can reach stops the run here. ---------------!
      if (src%humidity == HUMIDITY_RHAIR) then
         block
            integer(ik) :: j
            j = series_field(src, 'RHair')
            if (maxval(src%series(:, j)) > RH_FRACTION_MAX) then
               write(*,'(a,f0.2,a)') ' met_open: RHair reaches ', maxval(src%series(:, j)),                &
                                     ', but it is a fraction (units "1"), not a percentage.'
               call met_close(src)
               if (present(stat)) then ; stat = MET_ERR_HUMIDITY ; return ; end if
               error stop 'met_open: RHair is in percent (see the message above)'
            end if
         end block
      end if
   end subroutine met_open

   !----- CO2 (#184): with co2_source = "file", read the MEDS CO2 file and check it covers the  !
   !      run. MEDS does not extrapolate CO2, so a run the file does not cover stops here rather  !
   !      than decades in. The run span is absent only in unit tests. ------------------------------!
   subroutine open_co2(src, run_start, run_end, stat)
      type(met_source_t),          intent(inout) :: src
      type(meds_time_t), optional, intent(in)    :: run_start, run_end
      integer(ik),                 intent(out)   :: stat
      character(len=512) :: msg
      logical :: ok
      stat = MET_OK
      if (src%fcfg%co2_source /= CO2_SOURCE_FILE) return
      call co2_series_read(src%fcfg%co2_file, src%co2, ok, msg)
      if (.not. ok) then
         write(*,'(2a)') ' met_open: [forcing].co2_file = ', trim(src%fcfg%co2_file)
         write(*,'(2a)') '   ', trim(msg)
         stat = MET_ERR_CO2_FILE ; return
      end if
      if (present(run_start) .and. present(run_end)) then
         if (.not. co2_series_covers(src%co2, run_start, run_end)) then
            write(*,'(2a)') ' met_open: the CO2 file does not cover the run: ', trim(src%fcfg%co2_file)
            write(*,'(4a)') '   run  ', time_to_string(run_start), ' .. ', time_to_string(run_end)
            write(*,'(4a)') '   file ', time_to_string(src%co2%first_start), ' .. ',               &
                            time_to_string(co2_series_end(src%co2))
            write(*,'(a)')  '   MEDS does not extrapolate CO2: extend the file, or shorten the run.'
            stat = MET_ERR_CO2_NOT_COVERED ; return
         end if
      end if
   end subroutine open_co2

   !=======================================================================================!
   !  CURSOR: bind one polygon to the open source -- its cell and location -- and load its first  !
   !  bracket. A site run binds cell 1 at the [site] location; a region binds each polygon at its  !
   !  cell centre. The grid elevation (the lapse origin) is the cell's orography for the archive,  !
   !  the configured grid_elevation for a MEDS forcing file.                                       !
   !=======================================================================================!
   subroutine met_cursor_init(src, cur, cell, latitude_deg, longitude_deg, elevation_m)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(out) :: cur
      integer(ik),        intent(in)  :: cell
      real(wp),           intent(in)  :: latitude_deg, longitude_deg, elevation_m
      cur%cell = cell
      cur%latitude_deg = latitude_deg ; cur%longitude_deg = longitude_deg
      cur%elevation_m  = elevation_m
      cur%grid_elevation_m = src%fcfg%grid_elevation_m
      if (src%backend == MET_BACKEND_ED_ERA5LAND) cur%grid_elevation_m = src%cells%elevation(cell)
      if (src%backend == MET_BACKEND_CONST) then
         cur%rec_prev = met_record_t() ; cur%rec_next = met_record_t()
         return
      end if
      call load_bracket(src, cur, src%rec_first)
   end subroutine met_cursor_init

   !=======================================================================================!
   !  PREFETCH (R1): load what the step starting at `step_start` reads, before its compute      !
   !  phase, so met_advance never touches a file. A daily step from midnight reads one archive  !
   !  month plus the record before it: 00:00 on the 1st, which lives in the previous month's    !
   !  file, or the window's last record at the recycle wrap. Moving into the next month, that   !
   !  record comes from the outgoing buffer, so a run loads each month once. A step the recycle  !
   !  seam falls inside also reads the window's first day, which met_open keeps. validate_config !
   !  restricts format = "era5land" to daily steps from midnight, the shape this assumes.       !
   !=======================================================================================!
   subroutine met_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(meds_time_t),  intent(in)    :: step_start
      real(wp)    :: s0
      integer(ik) :: r, k, kp, p
      logical     :: wrap
      if (src%backend /= MET_BACKEND_ED_ERA5LAND) return
      s0 = file_lookup_sec(src, step_start)
      wrap = .false.
      if (src%n_cycle_years >= 1_ik .and. src%irec_cycle_last > src%irec_cycle_first) then
         wrap = s0 >= src%time_sec(src%irec_cycle_last) - REC_MATCH_TOL
      end if
      if (wrap) then                                    ! the window's first month, after its last record
         k = axis_month_of(src, src%irec_cycle_first) ; p = src%irec_cycle_last
      else
         r = first_record_after(src, s0)                ! clamped to the axis
         k = axis_month_of(src, r) ; p = src%month_rec0(k)
      end if
      if (p > 0_ik .and. src%carry_rec /= p) then       ! the record before the month
         if (.not. allocated(src%carry)) allocate(src%carry(src%cells%ncell, ERA_NVAR))
         kp = axis_month_of(src, p)
         if (.not. month_loaded(src, kp) .and. in_window_head(src, p)) then
            src%carry = src%head(p - src%irec_cycle_first + 1_ik, :, :)   ! a midnight anchor's record
         else
            if (.not. month_loaded(src, kp)) call load_axis_month(src, kp)
            src%carry = src%buffer%values(p - src%month_rec0(kp), :, :)
         end if
         src%carry_rec = p
      end if
      if (.not. month_loaded(src, k)) call load_axis_month(src, k)
   end subroutine met_prefetch

   !=======================================================================================!
   !  Validate the DECLARED recycle window ([forcing].recycle_start/recycle_end) against the    !
   !  file that was just opened. MEDS does NOT infer the window: it is told where the cycle       !
   !  starts and how long it is, and this routine's only job is to confirm the file agrees.       !
   !                                                                                              !
   !  Inferring the window is what this refuses to do. Wrapping a file on its own span shifts        !
   !  hour-of-day on EVERY wrap unless the span is whole years: the real ERA5-Land record (first     !
   !  stamp 01:00 under the end-of-interval convention) spans 366 d 22 h, and a 29-yr run wrapped    !
   !  on it reads late May at a ~10 h offset while its daily-mean shortwave stays correct, so the    !
   !  slow demography looks healthy and nothing surfaces the problem.                                !
   !=======================================================================================!
   subroutine validate_recycle_window(src, stat)
      type(met_source_t), intent(inout)  :: src
      integer(ik),        intent(out)   :: stat
      real(wp)    :: start_sec, end_sec
      integer(ik) :: i
      logical     :: covered

      stat = MET_OK
      src%cycle_anchor = src%fcfg%recycle_start
      src%n_cycle_years = 0_ik ; src%irec_cycle_first = 0_ik ; src%irec_cycle_last = 0_ik
      if (.not. src%fcfg%recycle) return

      !----- Whole-year span. Re-checked here (not only in validate_config) so a driver built     !
      !      directly from a forcing_config_t -- as the unit tests do -- cannot bypass it. --------!
      src%n_cycle_years = whole_years_between(src%fcfg%recycle_start, src%fcfg%recycle_end,     &
                                              MAX_RECYCLE_YEARS)
      if (src%n_cycle_years < 1_ik) then
         write(*,'(4a)') ' met_open: recycle window ', time_to_string(src%fcfg%recycle_start),   &
                         ' .. ', time_to_string(src%fcfg%recycle_end)
         write(*,'(a)')  '   is not an exact whole number of calendar years.'
         stat = MET_ERR_WINDOW_NOT_WHOLE_YEARS ; return
      end if

      start_sec = seconds_between(src%base_time, src%fcfg%recycle_start)
      end_sec   = seconds_between(src%base_time, src%fcfg%recycle_end)

      !----- V2: recycle_start must land EXACTLY on a record stamp. This is the check that stops  !
      !      MEDS guessing the start of the day: a config saying 00:00:00 against an ERA5-Land      !
      !      file stamped 01:00:00 is reported, with both values, rather than quietly re-derived.   !
      do i = 1_ik, src%nrec
         if (abs(src%time_sec(i) - start_sec) <= REC_MATCH_TOL) then
            src%irec_cycle_first = i ; exit
         end if
      end do
      if (src%irec_cycle_first == 0_ik) then
         write(*,'(2a)') ' met_open: forcing.recycle_start = ', time_to_string(src%fcfg%recycle_start)
         write(*,'(a)')  '   does not match any record in the forcing file. The file record stamps run'
         write(*,'(4a)') '   from ', time_to_string(time_advance_seconds(src%base_time, src%time_sec(1))), &
                         ' to ',    time_to_string(time_advance_seconds(src%base_time, src%time_sec(src%nrec)))
         write(*,'(a)')  '   Set recycle_start to an actual record stamp -- MEDS does not guess it. Note'
         write(*,'(a)')  '   avg_convention="end" files (ERA5-Land) stamp the END of each interval, so the'
         write(*,'(a)')  '   first record of a calendar year is 01:00:00 for hourly data, not 00:00:00.'
         stat = MET_ERR_START_NOT_A_RECORD ; return
      end if

      !----- V3: the file must cover the whole window, i.e. hold a record in the LAST interval    !
      !      below recycle_end. Without that the cycle has a hole at its far edge and the seam      !
      !      bracket below would interpolate across it. --------------------------------------------!
      do i = src%nrec, src%irec_cycle_first, -1_ik
         if (src%time_sec(i) < end_sec - REC_MATCH_TOL) then ; src%irec_cycle_last = i ; exit ; end if
      end do
      !----- NESTED, not a single .or.: Fortran does not guarantee short-circuit evaluation, so a  !
      !      combined test would index time_sec(irec_cycle_last) even when the loop above left it 0. !
      covered = src%irec_cycle_last > src%irec_cycle_first
      if (covered) covered = src%time_sec(src%irec_cycle_last) >= end_sec - src%dt_forcing - REC_MATCH_TOL
      if (.not. covered) then
         write(*,'(4a)') ' met_open: the forcing file does not cover the declared recycle window ', &
                         time_to_string(src%fcfg%recycle_start), ' .. ',                         &
                         time_to_string(src%fcfg%recycle_end)
         write(*,'(2a)') '   last record in the file: ',                                         &
                         time_to_string(time_advance_seconds(src%base_time, src%time_sec(src%nrec)))
         stat = MET_ERR_WINDOW_NOT_COVERED ; return
      end if

      write(*,'(a,i0,4a)') ' force : recycling ', src%n_cycle_years, ' calendar year(s), ',      &
            time_to_string(src%fcfg%recycle_start), ' .. ', time_to_string(src%fcfg%recycle_end)
   end subroutine validate_recycle_window

   !=======================================================================================!
   !  ADVANCE: slide the window so rec_prev%when <= now < rec_next%when (at this grid_index).    !
   !  Handles a start before the first record (hold or error), the recycle seam, and a run past  !
   !  the file's last record (clamped to the last interval).                                     !
   !=======================================================================================!
   subroutine met_advance(src, cur, now)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      type(meds_time_t),  intent(in)    :: now
      real(wp)    :: now_sec
      integer(ik) :: i

      if (src%backend == MET_BACKEND_CONST) return

      !----- the SAME effective (recycle-mapped) seconds on the FILE axis met_instant uses. ---!
      now_sec = file_lookup_sec(src, now)

      !----- start before the first record (never reached once recycle-wrapped into span). ---!
      !      Hold: load records 1-2 and let met_instant clamp w_next=0 for now < t1 -- do NOT     !
      !      overwrite rec_next: a stale copy there suppresses the reload once now reaches t1.    !
      if (now_sec < src%time_sec(1)) then
         if (src%fcfg%start_clamp == CLAMP_ERROR) then
            error stop 'met_advance: model start precedes the first forcing record (start_clamp=error)'
         end if
         if (cur%irec_prev /= 1_ik) call load_bracket(src, cur, 1_ik)
         return
      end if

      !----- CYCLE-BOUNDARY seam: the final interval of the declared window wraps back to its own  !
      !      first record, so the cycle is exactly periodic. The seam is the WINDOW's edge, not the  !
      !      file's -- a window may be a sub-range of a longer file, and trailing records past        !
      !      recycle_end belong to the next file year, not to this cycle. -----------------------------!
      !      NESTED, not one .and.: Fortran does not guarantee short-circuit evaluation, so a
      !      combined test would index time_sec(irec_cycle_last) on a non-recycling driver, where
      !      it is still 0.
      if (src%n_cycle_years >= 1_ik .and. src%irec_cycle_last > src%irec_cycle_first) then
         if (now_sec >= src%time_sec(src%irec_cycle_last)) then
            if (.not. cur%at_wrap_seam) call load_wrap_bracket(src, cur)
            return
         end if
      end if

      !----- non-recycle run past the end: clamp to the last interval (recycle wraps above). --!
      if (.not. src%fcfg%recycle .and. src%nrec >= 2_ik .and. now_sec >= src%time_sec(src%nrec)) then
         if (cur%irec_prev /= src%nrec - 1_ik) call load_bracket(src, cur, src%nrec - 1_ik)
         return
      end if

      !----- incremental cursor from the current bracket (cheap for a forward march). --------!
      i = max(1_ik, cur%irec_prev)
      do while (i < src%nrec - 1_ik)
         if (src%time_sec(i + 1_ik) <= now_sec) then ; i = i + 1_ik ; else ; exit ; end if
      end do
      do while (i > 1_ik)
         if (src%time_sec(i) > now_sec) then ; i = i - 1_ik ; else ; exit ; end if
      end do
      if (i /= cur%irec_prev) call load_bracket(src, cur, i)
      !----- #182: refresh the DAYTIME clearness memory the longwave synthesis falls back on after   !
      !      dark. Done HERE, not in met_instant, for two reasons: met_instant is intent(in) and is   !
      !      called per sub-step, and this block sits outside the patch loop, so the one scalar of    !
      !      state cannot become a data race once the patch axis is threaded.  ------------------------!
      if (src%fcfg%lwdown_source == LW_SYNTHESIZE) call remember_clearness(src, cur, now)
   end subroutine met_advance

   !----- Update cur%kt_last_day when the sun is up. The clearness index is formed from the        !
   !      bracketing record's TOTAL shortwave against the instantaneous TOA irradiance -- a mean    !
   !      over an interval against an instant at its sample point. That mismatch is deliberate and  !
   !      second-order here: the cloud term it feeds is a one-coefficient empirical correction on a !
   !      FALLBACK synthesis, and the alternative is to rebuild the interval-mean-cosz machinery    !
   !      for a factor that multiplies a 0.22 coefficient.  -----------------------------------------!
   subroutine remember_clearness(src, cur, now)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      type(meds_time_t),  intent(in)    :: now
      real(wp) :: cosz_now, sw_total, kt
      type(met_record_t) :: r
      r = interval_mean_record(src, cur)
      cosz_now = met_solar_cosz(now, seconds_into_day(now), cur%latitude_deg, cur%longitude_deg)
      sw_total = r%par_beam + r%par_diffuse + r%nir_beam + r%nir_diffuse
      kt = clearness_index(sw_total, cosz_now)
      if (kt >= 0.0_wp) cur%kt_last_day = kt         ! negative = night: keep dusk's value
   end subroutine remember_clearness

   !----- The record carrying the means over the interval that contains the model instant: the  !
   !      later of the bracket on an avg_convention = "end" file, the earlier on a "begin" one.     !
   !      Every consumer of an interval mean -- rain, shortwave, the clearness the longwave         !
   !      synthesis remembers -- takes it from here, so none of them reads a neighbouring interval.  !
   pure function interval_mean_record(src, cur) result(r)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(met_record_t) :: r
      select case (src%fcfg%avg_convention)
      case (METAVG_BEGIN) ; r = cur%rec_prev
      case default        ; r = cur%rec_next            ! METAVG_END (ERA5-Land) + fallback
      end select
   end function interval_mean_record

   !=======================================================================================!
   !  INSTANT: interpolate/disaggregate the loaded window to the model instant `now`.           !
   !  State vars linear; wind energy-form; rainfall step-constant then phase-split; shortwave via    !
   !  the reciprocal-mean-cosz reconstruction of the interval-mean streams; cosz + rho_air derived. !
   !=======================================================================================!
   function met_instant(src, cur, now) result(met)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      type(meds_time_t),  intent(in) :: now
      type(met_forcing_t) :: met
      type(met_record_t)  :: mean_rec
      type(meds_time_t)   :: mws
      real(wp) :: now_sec, tprev, tnext, w_next, cosz_now, factor, win_start_sec, precip_total
      associate (f => src%fcfg, p => cur%rec_prev, n => cur%rec_next)

      !----- solar zenith at `now` (UTC -> apparent solar seconds inside met_solar_cosz). ------!
      cosz_now = met_solar_cosz(now, seconds_into_day(now), cur%latitude_deg, cur%longitude_deg)
      met%cosz = cosz_now

      !----- CO2 on MODEL time, one way for every backend (#184). ------------------------------!
      if (f%co2_source == CO2_SOURCE_FILE) then
         met%co2 = co2_series_at(src%co2, now)
      else
         met%co2 = f%co2_const
      end if

      !----- The reference climate (CONST) is the forcing type's defaults, held flat; a file backend  !
      !      interpolates its bracket to `now`.                                                    !
      if (src%backend /= MET_BACKEND_CONST) then
         !----- interpolation weight within the loaded window (SAME effective seconds as advance, !
         !      so a recycle-mapped now interpolates within the recycled interval, not clamps to 1). !
         now_sec = file_lookup_sec(src, now)
         if (cur%at_wrap_seam) then                    ! cycle-boundary: window's last rec -> its first
            tprev = src%time_sec(src%irec_cycle_last) ; tnext = tprev + src%dt_forcing
         else
            tprev = src%time_sec(cur%irec_prev)
            tnext = src%time_sec(min(cur%irec_prev + 1_ik, src%nrec))
         end if
         if (tnext > tprev) then
            w_next = min(1.0_wp, max(0.0_wp, (now_sec - tprev) / (tnext - tprev)))
         else
            w_next = 0.0_wp                                          ! degenerate (start-clamp hold)
         end if

         !----- state variables: linear; wind: energy-conserving. ------------------------------!
         met%tair_k   = interpolate_forcing(INTERP_LINEAR, p%tair_k,   n%tair_k,   w_next)
         met%qair     = interpolate_forcing(INTERP_LINEAR, p%qair,     n%qair,     w_next)
         met%psurf_pa = interpolate_forcing(INTERP_LINEAR, p%psurf_pa, n%psurf_pa, w_next)
         met%lwdown   = interpolate_forcing(INTERP_LINEAR, p%lwdown,   n%lwdown,   w_next)
         met%wind     = interpolate_wind_energy(p%wind, n%wind, w_next, U_MIN)
         !----- The vector interpolates linearly, which keeps its direction (§5.3); never floored. ---!
         if (src%has_wind_vector) then
            met%wind_u = interpolate_forcing(INTERP_LINEAR, p%wind_u, n%wind_u, w_next)
            met%wind_v = interpolate_forcing(INTERP_LINEAR, p%wind_v, n%wind_v, w_next)
            met%has_wind_vector = .true.
         end if

         !----- The fluxes are interval means: the interval CONTAINING now is [prev, next], whose     !
         !      mean is rec_next's on an avg_convention = "end" file and rec_prev's on a "begin" one.    !
         !      Rain and shortwave both come from that record, so neither lags the other.              !
         mean_rec = interval_mean_record(src, cur)

         !----- rainfall: the interval's total rate, held across it (never smeared), then split by  !
         !      phase on the interpolated temperature. ------------------------------------------------!
         precip_total = mean_rec%rainf
         call precip_phase(precip_total, met%tair_k, met%rainf, met%snowfall)

         !----- shortwave: that interval's mean streams, disaggregated by cosz(now)/<cosz>_win. ----!
         !----- reconstruction factor anchored on the MODEL window start (mws), so <cosz>_win aligns  !
         !      with cosz_now on the model calendar (identity = rec_prev%when when not recycling; under  !
         !      calendar recycling it follows the model sun, so the interval-mean identity still holds).  !
         mws           = time_advance_seconds(now, tprev - now_sec)   ! model instant at the window start
         win_start_sec = seconds_into_day(mws)
         factor = cosz_reconstruct_factor(mws, win_start_sec,                                       &
                                          (tnext - tprev) / real(N_COSZ_SUB, wp), tnext - tprev,   &
                                          cur%latitude_deg, cur%longitude_deg)
         met%par_beam    = disaggregate_sw(mean_rec%par_beam,    cosz_now, factor)
         met%par_diffuse = disaggregate_sw(mean_rec%par_diffuse, cosz_now, factor)
         met%nir_beam    = disaggregate_sw(mean_rec%nir_beam,    cosz_now, factor)
         met%nir_diffuse = disaggregate_sw(mean_rec%nir_diffuse, cosz_now, factor)
      end if

      !----- LONGWAVE SYNTHESIS (#182), for a source that carries no LWdown. Placed AFTER the      !
      !      shortwave block so the instantaneous streams are available: by day the cloud term uses  !
      !      this instant's clearness index, and after dark it falls back to the last daytime value  !
      !      met_advance remembered. Overwrites whatever the file read, which is the point -- the    !
      !      selector says the file's longwave is not to be trusted or is not there.  ----------------!
      if (f%lwdown_source == LW_SYNTHESIZE) then
         block
            real(wp) :: sw_now, kt_now
            sw_now = met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse
            kt_now = clearness_index(sw_now, cosz_now)
            if (kt_now < 0.0_wp) kt_now = cur%kt_last_day
            met%lwdown = synthesize_lwdown(f%lw_clear_form, met%tair_k, met%qair, met%psurf_pa,   &
                                           kt_now, f%lw_cloud_a)
         end block
      end if

      met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)
      end associate
   end function met_instant

   !=======================================================================================!
   !  CLOSE.                                                                                     !
   !=======================================================================================!
   subroutine met_close(src)
      type(met_source_t), intent(inout)  :: src
      integer(c_int) :: st
      if (src%backend == MET_BACKEND_ED_DEFAULT .and. src%ncid >= 0_ik) then
         st = nc_close(int(src%ncid, c_int)) ; call nc_check(st, 'met_close: nc_close')
         src%ncid = -1_ik
      end if
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      if (allocated(src%month_year))  deallocate(src%month_year, src%month_month, src%month_rec0)
      if (allocated(src%buffer%values)) deallocate(src%buffer%values)
      if (allocated(src%carry))  deallocate(src%carry)
      if (allocated(src%head))   deallocate(src%head)
      if (allocated(src%series)) deallocate(src%series, src%series_name)
      call co2_series_free(src%co2)
      src%buffer%year = 0_ik ; src%carry_rec = 0_ik ; src%n_head = 0_ik ; src%n_loads = 0_ik
   end subroutine met_close

   !----- The effective seconds-since-base on the FILE time axis, used by BOTH bracket selection  !
   !      (met_advance) and the interpolation weight (met_instant) so they stay consistent. Under   !
   !      calendar recycling the model date maps into the declared window (recycle_model_to_file),   !
   !      keeping month/day/hour so day-of-year is exact across leap boundaries; otherwise it is the !
   !      model date itself.                                                                          !
   pure function file_lookup_sec(src, now) result(s)
      type(met_source_t), intent(in)  :: src
      type(meds_time_t),  intent(in) :: now
      real(wp) :: s
      if (src%fcfg%recycle .and. src%n_cycle_years >= 1_ik) then
         s = seconds_between(src%base_time, recycle_model_to_file(src, now))
      else
         s = seconds_between(src%base_time, now)
      end if
   end function file_lookup_sec

   !----- Map a model instant into the declared recycle window (ED2 year_use substitution), keeping !
   !      month/day/hour/minute/second EXACT -- that exactness is the whole point: it is what makes   !
   !      the sub-daily phase and the day-of-year survive an arbitrary number of wraps.               !
   !                                                                                                  !
   !      ANCHOR-RELATIVE, so the window may begin anywhere in the calendar rather than only on Jan-1: !
   !                                                                                                   !
   !         off = 1 when the model's (month,day,hh:mm:ss) falls BEFORE the window's own anchor         !
   !               month/day/time-of-day -- that instant belongs to the NEXT file year of the cycle.    !
   !         yf  = anchor_year + off + modulo(model_year - anchor_year - off, n_cycle_years)            !
   !                                                                                                   !
   !      For a Jan-1 00:00:00 anchor `off` is identically 0 and this reduces, term for term, to the    !
   !      plain year substitution. Without `off`, a mid-year window (or an end-of-interval file whose   !
   !      first stamp is 01:00) maps instants OUTSIDE the window it is supposed to cycle over.         !
   !                                                                                                   !
   !      LEAP DAY: Feb-29 -> Feb-28 when the target file year is non-leap (ED2 read_ol_file repeats    !
   !      Feb 28). A non-leap model year never asks for Feb-29, so the substitution is one-directional. !
   pure function recycle_model_to_file(src, now) result(now_file)
      type(met_source_t), intent(in)  :: src
      type(meds_time_t),  intent(in) :: now
      type(meds_time_t) :: now_file
      integer(ik) :: yf, off
      off = merge(0_ik, 1_ik, mdhms_at_or_after(now, src%cycle_anchor))
      yf  = src%cycle_anchor%year + off                                                          &
          + modulo(now%year - src%cycle_anchor%year - off, src%n_cycle_years)
      now_file = now
      now_file%year = yf
      if (now%month == 2_ik .and. now%day == 29_ik .and. .not. is_leap_year(yf)) now_file%day = 28_ik
   end function recycle_model_to_file

   !----- Compare two instants on (month, day, hour, minute, second) ONLY -- i.e. their position    !
   !      within the calendar year, ignoring which year they fall in. -------------------------------!
   pure logical function mdhms_at_or_after(a, b) result(yes)
      type(meds_time_t), intent(in) :: a, b
      integer(ik) :: ka, kb
      ka = ((a%month*32_ik + a%day)*24_ik + a%hour)*3600_ik + a%minute*60_ik + a%second
      kb = ((b%month*32_ik + b%day)*24_ik + b%hour)*3600_ik + b%minute*60_ik + b%second
      yes = ka >= kb
   end function mdhms_at_or_after


   !=======================================================================================!
   !  Helpers (NetCDF backend).                                                                 !
   !=======================================================================================!
   !----- Read the whole time coordinate (seconds since base_time) into src%time_sec. ---------!
   subroutine read_time_axis(src)
      type(met_source_t), intent(inout)  :: src
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start1(1), count1(1)
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      allocate(src%time_sec(src%nrec))
      st = nc_inq_varid_f(int(src%ncid, c_int), 'time', vid) ; call nc_check(st, 'read_time_axis: time varid')
      start1(1) = 0_c_size_t ; count1(1) = int(src%nrec, c_size_t)
      st = nc_get_vara_double(int(src%ncid, c_int), vid, start1, count1, src%time_sec)
      call nc_check(st, 'read_time_axis: time values')
   end subroutine read_time_axis

   !----- Load rec_prev = record(irec), rec_next = record(irec+1) (clamped at EOF). -----------!
   subroutine load_bracket(src, cur, irec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      integer(ik),        intent(in)    :: irec
      call read_record(src, cur, irec, cur%rec_prev)
      call read_record(src, cur, min(irec + 1_ik, src%nrec), cur%rec_next)
      cur%irec_prev    = irec
      cur%at_wrap_seam = .false.
   end subroutine load_bracket

   !----- Cycle-boundary seam bracket: rec_prev = record nrec, rec_next = record 1, so the last  !
   !      dt interval interpolates across the wrap instead of clamping to the final record.        !
   subroutine load_wrap_bracket(src, cur)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      call read_record(src, cur, src%irec_cycle_last,  cur%rec_prev)
      call read_record(src, cur, src%irec_cycle_first, cur%rec_next)
      cur%irec_prev    = src%irec_cycle_last
      cur%at_wrap_seam = .true.
   end subroutine load_wrap_bracket

   !=======================================================================================!
   !  ED_ERA5land OPEN (§15): resolve the paths, bind the site to its cell (whose orography       !
   !  replaces a hand-copied grid_elevation), and lay the months the run needs end to end as one  !
   !  hourly axis. Every month file is checked for existence here, so a gap in the archive stops  !
   !  the run before it starts rather than decades into it.                                       !
   !=======================================================================================!
   subroutine open_archive(src, run_start, run_end, stat, cells)
      type(met_source_t),          intent(inout)  :: src
      type(meds_time_t), optional, intent(in)    :: run_start, run_end
      integer(ik),                 intent(out)   :: stat
      type(met_cells_t), optional, intent(in)    :: cells
      character(len=MET_PATH_LEN) :: path
      real(wp)       :: distance_km, first
      integer(ik)    :: est, y0, m0, y1, m1, nmonth, k, v, h, y, m, nmiss
      integer(c_int) :: st, ncid
      logical        :: exists

      stat = MET_OK
      associate (f => src%fcfg)
      src%file_template = f%file_template
      if (len_trim(src%file_template) == 0) src%file_template = era5land_default_template()
      path = f%static_file
      if (len_trim(path) == 0) path = era5land_default_static()
      src%static_file = era5land_path(path, f%data_path, '', 0_ik, 0_ik)

      if (f%sw_partition == SWPART_PASSTHROUGH) then
         write(*,'(a)') ' met_open: the ED_ERA5land archive carries TOTAL shortwave, but'
         write(*,'(a)') '   [forcing].sw_partition = "passthrough" expects the four component streams.'
         stat = MET_ERR_ATTR_MISMATCH ; return
      end if

      !----- The cells: a region's, selected by its caller, or the site's own (§15.5). ------------!
      est = ERA_OK
      if (present(cells)) then
         src%cells = cells
      else
         call era5land_select_site(src%static_file, f%latitude_deg, f%longitude_deg, f%max_distance_km, &
                                   src%cells, distance_km, est)
      end if
      if (est /= ERA_OK) then
         select case (est)
         case (ERA_ERR_NO_CELL)
            write(*,'(a,f0.1,a,f0.4,a,f0.4,a)') ' met_open: no valid ERA5-Land cell within ',          &
                  f%max_distance_km, ' km of the site (', f%latitude_deg, ', ', f%longitude_deg, ')'
         case (ERA_ERR_OPEN)
            write(*,'(2a)') ' met_open: cannot open the ED_ERA5land static file ', trim(src%static_file)
         case default
            write(*,'(3a)') ' met_open: the static file ', trim(src%static_file),                      &
                            ' does not hold a regular lat/lon grid with valid and elevation'
         end select
         stat = MET_ERR_ARCHIVE ; return
      end if
      src%grid_index = 1_ik ; src%ngrid = src%cells%ncell
      src%has_wind_vector = .true.
      src%humidity = HUMIDITY_TDEW                           ! the archive stores the 2 m dewpoint
      if (present(cells)) then
         write(*,'(a,i0,a)') ' force : ED_ERA5land region, ', src%cells%ncell, ' cells'
      else
         write(*,'(a,f0.3,a,f0.3,a,f7.2,a,f0.1,a)') ' force : ED_ERA5land cell (', src%cells%lat(1), ', ',  &
               src%cells%lon(1), '), ', distance_km, ' km from the site; orography ', src%cells%elevation(1), ' m'
      end if

      !----- The months: the declared recycle window, or the run period. ----------------------!
      if (f%recycle) then
         call archive_month_of(f%recycle_start, .false., y0, m0)
         call archive_month_of(time_advance_seconds(f%recycle_end, -src%dt_forcing), .false., y1, m1)
      else
         if (.not. (present(run_start) .and. present(run_end)))                                     &
            error stop 'met_open: format = "era5land" without recycling needs the run period'
         call archive_month_of(run_start, .false., y0, m0)
         call archive_month_of(run_end,   .true.,  y1, m1)
      end if
      nmonth = (y1 - y0) * 12_ik + (m1 - m0) + 1_ik
      if (nmonth < 1_ik) error stop 'met_open: the ED_ERA5land month range is empty'
      allocate(src%month_year(nmonth), src%month_month(nmonth), src%month_rec0(nmonth))
      y = y0 ; m = m0 ; src%nrec = 0_ik
      do k = 1_ik, nmonth
         src%month_year(k) = y ; src%month_month(k) = m ; src%month_rec0(k) = src%nrec
         src%nrec = src%nrec + era5land_month_hours(y, m)
         m = m + 1_ik
         if (m > 12_ik) then ; m = 1_ik ; y = y + 1_ik ; end if
      end do
      src%base_time = ERA_EPOCH
      if (allocated(src%time_sec)) deallocate(src%time_sec)
      allocate(src%time_sec(src%nrec))
      do k = 1_ik, nmonth
         first = seconds_between(ERA_EPOCH, meds_time_t(src%month_year(k), src%month_month(k), 1_ik, 1_ik))
         do h = 1_ik, era5land_month_hours(src%month_year(k), src%month_month(k))
            src%time_sec(src%month_rec0(k) + h) = first + 3600.0_wp * real(h - 1_ik, wp)
         end do
      end do

      nmiss = 0_ik
      do k = 1_ik, nmonth
         do v = 1_ik, ERA_NVAR
            path = era5land_path(src%file_template, f%data_path, ERA_VAR_NAME(v),                  &
                                 src%month_year(k), src%month_month(k))
            inquire(file=trim(path), exist=exists)
            if (.not. exists) then
               nmiss = nmiss + 1_ik
               if (nmiss == 1_ik) write(*,'(2a)') ' met_open: ED_ERA5land file missing: ', trim(path)
            end if
         end do
      end do
      if (nmiss > 0_ik) then
         write(*,'(a,i0,a,i0,a)') '   ', nmiss, ' of ', nmonth * ERA_NVAR, ' files the run needs are missing.'
         stat = MET_ERR_ARCHIVE ; return
      end if

      !----- The first month file describes the archive; check it against the config (#185). ---!
      path = era5land_path(src%file_template, f%data_path, 'Tair', y0, m0)
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) then
         write(*,'(2a)') ' met_open: cannot open ', trim(path)
         stat = MET_ERR_ARCHIVE ; return
      end if
      call validate_file_against_config(src, ncid, stat)
      st = nc_close(ncid)
      src%buffer%year = 0_ik ; src%carry_rec = 0_ik ; src%n_loads = 0_ik
      end associate
   end subroutine open_archive

   !----- The archive month whose file holds the record at t, taken to the hour below (or above,   !
   !      round_up) on the hourly grid. Files are end-stamped, so 00:00 on the 1st belongs to the   !
   !      month before. -----------------------------------------------------------------------------!
   subroutine archive_month_of(t, round_up, year, month)
      type(meds_time_t), intent(in)  :: t
      logical,           intent(in)  :: round_up
      integer(ik),       intent(out) :: year, month
      type(meds_time_t) :: before
      real(wp) :: hours
      hours = seconds_between(ERA_EPOCH, t) / 3600.0_wp
      if (round_up) then ; hours = real(ceiling(hours), wp) ; else ; hours = real(floor(hours), wp) ; end if
      before = time_advance_seconds(ERA_EPOCH, hours * 3600.0_wp - 1.0_wp)
      year = before%year ; month = before%month
   end subroutine archive_month_of

   !----- The axis month holding record irec (binary search: an axis can span decades). ----------!
   pure integer(ik) function axis_month_of(src, irec) result(k)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      integer(ik) :: lo, hi, mid
      lo = 1_ik ; hi = size(src%month_rec0, kind=ik)          ! month_rec0(lo) < irec always holds
      do while (lo < hi)
         mid = (lo + hi + 1_ik) / 2_ik
         if (src%month_rec0(mid) < irec) then ; lo = mid ; else ; hi = mid - 1_ik ; end if
      end do
      k = lo
   end function axis_month_of

   !----- The first record strictly after file time s (the last record if none). ------------------!
   pure integer(ik) function first_record_after(src, s) result(r)
      type(met_source_t), intent(in)  :: src
      real(wp),           intent(in) :: s
      integer(ik) :: lo, hi, mid
      lo = 1_ik ; hi = src%nrec
      do while (lo < hi)
         mid = (lo + hi) / 2_ik
         if (src%time_sec(mid) > s + REC_MATCH_TOL) then ; hi = mid ; else ; lo = mid + 1_ik ; end if
      end do
      r = lo
   end function first_record_after

   pure logical function month_loaded(src, k) result(yes)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: k
      yes = src%buffer%year == src%month_year(k) .and. src%buffer%month == src%month_month(k)
   end function month_loaded

   pure logical function in_window_head(src, irec) result(yes)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      yes = irec >= src%irec_cycle_first .and. irec < src%irec_cycle_first + src%n_head
   end function in_window_head

   !----- Keep the recycle window's first day, from axis month 1 just loaded: its first record    !
   !      through the next midnight, all in that month's file. The seam falls inside a daily step  !
   !      unless the window starts at 01:00, and that step reads the window's last record, in the  !
   !      axis's last month, and then these. A midnight anchor needs its one record.               !
   subroutine keep_window_head(src)
      type(met_source_t), intent(inout) :: src
      real(wp)    :: time_of_day
      integer(ik) :: h1
      src%n_head = 0_ik
      if (allocated(src%head)) deallocate(src%head)
      if (src%n_cycle_years < 1_ik) return
      time_of_day = seconds_into_day(src%cycle_anchor)
      src%n_head = 1_ik
      if (time_of_day > REC_MATCH_TOL)                                                            &
         src%n_head = 1_ik + nint((86400.0_wp - time_of_day) / src%dt_forcing, ik)
      h1 = src%irec_cycle_first - src%month_rec0(1)
      allocate(src%head(src%n_head, src%cells%ncell, ERA_NVAR))
      src%head = src%buffer%values(h1:h1 + src%n_head - 1_ik, :, :)
   end subroutine keep_window_head

   !----- Read axis month k into the buffer: the ONLY place the archive is read after open. -----!
   subroutine load_axis_month(src, k)
      type(met_source_t), intent(inout)  :: src
      integer(ik),        intent(in)    :: k
      character(len=MET_PATH_LEN + 128) :: message
      integer(ik) :: est
      call era5land_load_month(src%file_template, src%fcfg%data_path, src%cells, src%month_year(k), &
                               src%month_month(k), src%buffer, est, message)
      if (est /= ERA_OK) then
         write(*,'(2a)') ' met_driver: ED_ERA5land month failed its checks: ', trim(message)
         error stop 'met_driver: an ED_ERA5land month cannot be read (MEDS does not gap-fill)'
      end if
      src%n_loads = src%n_loads + 1_ik
   end subroutine load_axis_month

   !----- Where record irec's values are: its hour h > 0 in the loaded month, h = 0 for the carried !
   !      record, h < 0 for the -h-th record of the window's first day. Anything else was not       !
   !      prefetched, which is a programming error, not a data gap. -------------------------------!
   subroutine locate_record(src, irec, h)
      type(met_source_t), intent(in)   :: src
      integer(ik),        intent(in)  :: irec
      integer(ik),        intent(out) :: h
      integer(ik) :: k
      k = axis_month_of(src, irec)
      if (month_loaded(src, k)) then
         h = irec - src%month_rec0(k)
      else if (irec == src%carry_rec) then
         h = 0_ik
      else if (in_window_head(src, irec)) then
         h = src%irec_cycle_first - 1_ik - irec
      else
         write(*,'(a,i0,2a)') ' met_driver: forcing record ', irec, ' at ',                            &
               time_to_string(time_advance_seconds(src%base_time, src%time_sec(irec)))
         error stop 'met_driver: a forcing record was not prefetched before the step (internal error)'
      end if
   end subroutine locate_record

   pure real(wp) function archive_value(src, cur, h, var) result(val)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in) :: h, var
      if (h == 0_ik) then
         val = real(src%carry(cur%cell, var), wp)
      else if (h < 0_ik) then
         val = real(src%head(-h, cur%cell, var), wp)
      else
         val = real(src%buffer%values(h, cur%cell, var), wp)
      end if
   end function archive_value

   !----- Resolve src%grid_index by nearest great-circle distance from [site] lat/lon to the file's !
   !      latitude(grid)/longitude(grid) coordinate vectors (multi-polygon P2 subset, §4.1).        !
   subroutine resolve_grid_index(src)
      type(met_source_t), intent(inout)  :: src
      integer(c_int)    :: st, vlat, vlon
      integer(c_size_t) :: start1(1), count1(1)
      real(wp), allocatable :: lat_grid(:), lon_grid(:)
      allocate(lat_grid(src%ngrid), lon_grid(src%ngrid))
      start1(1) = 0_c_size_t ; count1(1) = int(src%ngrid, c_size_t)
      st = nc_inq_varid_f(int(src%ncid, c_int), 'latitude', vlat)  ; call nc_check(st, 'resolve_grid_index: latitude varid')
      st = nc_get_vara_double(int(src%ncid, c_int), vlat, start1, count1, lat_grid)
      call nc_check(st, 'resolve_grid_index: latitude values')
      st = nc_inq_varid_f(int(src%ncid, c_int), 'longitude', vlon) ; call nc_check(st, 'resolve_grid_index: longitude varid')
      st = nc_get_vara_double(int(src%ncid, c_int), vlon, start1, count1, lon_grid)
      call nc_check(st, 'resolve_grid_index: longitude values')
      src%grid_index = nearest_grid_index(src%fcfg%longitude_deg, src%fcfg%latitude_deg, lon_grid, lat_grid)
      write(*,'(a,i0,a,f8.3,a,f8.3,a)') 'met_open: grid_match=nearest resolved grid_index=', src%grid_index, &
            ' (site lon=', src%fcfg%longitude_deg, ' lat=', src%fcfg%latitude_deg, ')'
      deallocate(lat_grid, lon_grid)
   end subroutine resolve_grid_index

   !----- Read one record at (time=irec, this polygon's cell); partition SW at ingest. --------!
   subroutine read_record(src, cur, irec, rec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(out)   :: rec
      real(wp)    :: sw_total, cosz_mid, mid_sec, humidity_value
      integer(ik) :: h
      rec%when   = time_advance_seconds(src%base_time, src%time_sec(irec))
      if (src%backend == MET_BACKEND_ED_ERA5LAND) then
         !----- The archive stores dewpoint and the wind components (§15.4): humidity comes from the !
         !      model's own saturation curve, the speed from the vector at each stamp. -------------!
         call locate_record(src, irec, h)
         rec%tair_k   = archive_value(src, cur, h, ERA_TAIR)
         rec%psurf_pa = archive_value(src, cur, h, ERA_PSURF)
         humidity_value = archive_value(src, cur, h, ERA_TDEW)
         rec%wind_u   = archive_value(src, cur, h, ERA_U10)
         rec%wind_v   = archive_value(src, cur, h, ERA_V10)
         rec%wind     = sqrt(rec%wind_u**2 + rec%wind_v**2)
         rec%rainf    = archive_value(src, cur, h, ERA_RAINF)            ! total rainfall rate [kg/m2/s]
         rec%lwdown   = archive_value(src, cur, h, ERA_LWDOWN)
      else
         rec%tair_k   = read_scalar(src, 'Tair',  irec)
         humidity_value = read_scalar(src, humidity_name(src%humidity), irec)
         rec%psurf_pa = read_scalar(src, 'PSurf', irec)
         if (src%has_wind_vector) then
            rec%wind_u = read_scalar(src, 'u10', irec)
            rec%wind_v = read_scalar(src, 'v10', irec)
            rec%wind   = sqrt(rec%wind_u**2 + rec%wind_v**2)
         else
            rec%wind   = read_scalar(src, 'Wind',  irec)
         end if
         rec%rainf    = read_scalar(src, 'Rainf', irec)              ! total rainfall rate [kg/m2/s]
         !----- LWdown is OPTIONAL when we are synthesizing it (#182): a source without longwave is    !
         !      exactly the case lwdown_source = "synthesize" exists for, so demanding the variable      !
         !      would defeat the feature. Read it when present either way -- it costs nothing and keeps  !
         !      the record complete for diagnostics.  ---------------------------------------------------!
         if (src%fcfg%lwdown_source == LW_SYNTHESIZE) then
            rec%lwdown = read_scalar_default(src, 'LWdown', irec, 0.0_wp)
         else
            rec%lwdown = read_scalar(src, 'LWdown', irec)
         end if
      end if
      call assert_finite(rec%tair_k, 'Tair', irec, src%grid_index)
      call assert_finite(humidity_value, humidity_name(src%humidity), irec, src%grid_index)
      call assert_finite(rec%psurf_pa, 'PSurf', irec, src%grid_index)
      !----- Specific humidity at the forcing's own temperature and pressure. A dewpoint becomes q  !
      !      exactly as the archive reader always made it; relative humidity goes through the same   !
      !      liquid-water curve (docs/science/forcing.md §7). --------------------------------------!
      select case (src%humidity)
      case (HUMIDITY_RHAIR) ; rec%qair = rh_to_specific_humidity(humidity_value, rec%tair_k, rec%psurf_pa)
      case (HUMIDITY_TDEW)  ; rec%qair = dewpoint_to_specific_humidity(humidity_value, rec%psurf_pa)
      case default          ; rec%qair = humidity_value
      end select
      call assert_finite(rec%wind, 'Wind', irec, src%grid_index)
      call assert_finite(rec%rainf, 'Rainf', irec, src%grid_index)
      if (src%fcfg%lwdown_source /= LW_SYNTHESIZE)                                                  &
         call assert_finite(rec%lwdown, 'LWdown', irec, src%grid_index)

      !----- TERRAIN lapse at ingest (docs/science/forcing.md §8), from the cell's elevation to the  !
      !      site's: T by the month's lapse rate, P hydrostatically with the same linear T(z), q at     !
      !      constant relative humidity, file longwave by the clear-sky eps*T^4 ratio. The move to each !
      !      patch's canopy-air top is the fast loop's (per patch), not the reader's.                  !
      if (src%fcfg%apply_elevation_lapse) then
         block
            real(wp) :: dz, gamma, t_grid, p_grid, q_grid
            dz     = cur%elevation_m - cur%grid_elevation_m       ! + when site is higher
            gamma  = monthly_lapse_rate(src%fcfg%lapse_rate_tair, rec%when%month)
            t_grid = rec%tair_k ; p_grid = rec%psurf_pa ; q_grid = rec%qair   ! capture BEFORE overwrite
            rec%psurf_pa = lapse_pressure(p_grid, t_grid, dz, gamma)
            rec%tair_k   = lapse_air_temperature(t_grid, dz, gamma)
            !----- Relative humidity is held across the lapse: a measured one is used directly, and a  !
            !      dewpoint or specific humidity gives its own at the cell. ----------------------------!
            if (src%humidity == HUMIDITY_RHAIR) then
               rec%qair  = rh_to_specific_humidity(humidity_value, rec%tair_k, rec%psurf_pa)
            else
               rec%qair  = lapse_specific_humidity(q_grid, t_grid, p_grid, rec%tair_k, rec%psurf_pa)
            end if
            !----- A synthesized longwave is built later from the lapsed T and q, so only a file's   !
            !      longwave is scaled here. ------------------------------------------------------------!
            if (src%fcfg%lwdown_source /= LW_SYNTHESIZE)                                           &
               rec%lwdown = lapse_longwave(rec%lwdown, src%fcfg%lw_clear_form, t_grid, q_grid, p_grid, &
                                           rec%tair_k, rec%qair, rec%psurf_pa)
         end block
      end if

      if (src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
         rec%par_beam    = read_scalar(src, 'SWdown_par_beam',    irec)
         rec%par_diffuse = read_scalar(src, 'SWdown_par_diffuse', irec)
         rec%nir_beam    = read_scalar(src, 'SWdown_nir_beam',    irec)
         rec%nir_diffuse = read_scalar(src, 'SWdown_nir_diffuse', irec)
         call assert_finite(rec%par_beam,    'SWdown_par_beam',    irec, src%grid_index)   ! required source
         call assert_finite(rec%par_diffuse, 'SWdown_par_diffuse', irec, src%grid_index)   ! fields -> no gap-fill
         call assert_finite(rec%nir_beam,    'SWdown_nir_beam',    irec, src%grid_index)
         call assert_finite(rec%nir_diffuse, 'SWdown_nir_diffuse', irec, src%grid_index)
      else
         if (src%backend == MET_BACKEND_ED_ERA5LAND) then
            sw_total = archive_value(src, cur, h, ERA_SWDOWN)
         else
            sw_total = read_scalar(src, 'SWdown', irec)
         end if
         call assert_finite(sw_total, 'SWdown', irec, src%grid_index)
         !----- interval-mean cosz for the partition (avg_convention=end -> midpoint = when - dt/2). !
         mid_sec  = seconds_into_day(rec%when)
         if (src%fcfg%avg_convention == METAVG_END)   mid_sec = mid_sec - 0.5_wp * src%dt_forcing
         if (src%fcfg%avg_convention == METAVG_BEGIN)  mid_sec = mid_sec + 0.5_wp * src%dt_forcing
         cosz_mid = met_solar_cosz(rec%when, mid_sec, cur%latitude_deg, cur%longitude_deg)
         call partition_shortwave(sw_total, cosz_mid, rec%psurf_pa, src%fcfg%sw_partition,       &
                                  rec%par_beam, rec%par_diffuse, rec%nir_beam, rec%nir_diffuse)
      end if
   end subroutine read_record

   !----- A MEDS forcing file's fields at the cell, records r0..r1, read at open. ---------------!
   subroutine read_series(src, ncid, r0, r1)
      type(met_source_t), intent(inout)  :: src
      integer(c_int),     intent(in)    :: ncid
      integer(ik),        intent(in)    :: r0, r1
      character(len=24), parameter :: FIELDS(15) = [character(len=24) ::                            &
         'Tair', 'Qair', 'RHair', 'Tdew', 'PSurf', 'Wind', 'u10', 'v10', 'Rainf', 'LWdown', 'SWdown', &
         'SWdown_par_beam', 'SWdown_par_diffuse', 'SWdown_nir_beam', 'SWdown_nir_diffuse']
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start2(2), count2(2)
      logical     :: present_(size(FIELDS))
      integer(ik) :: j, n
      do j = 1, size(FIELDS)
         present_(j) = nc_inq_varid_f(ncid, trim(FIELDS(j)), vid) == NC_NOERR
      end do
      n = int(count(present_), ik)
      if (allocated(src%series)) deallocate(src%series, src%series_name)
      allocate(src%series(r0:r1, n), src%series_name(n))              ! indexed by record number
      src%series_name = pack(FIELDS, present_)
      start2 = [int(r0 - 1_ik, c_size_t), int(src%grid_index - 1_ik, c_size_t)]    ! [time, grid]
      count2 = [int(r1 - r0 + 1_ik, c_size_t), 1_c_size_t]
      do j = 1_ik, n
         st = nc_inq_varid_f(ncid, trim(src%series_name(j)), vid)
         st = nc_get_vara_double(ncid, vid, start2, count2, src%series(:, j))
         call nc_check(st, 'read_series: get '//trim(src%series_name(j)))
      end do
   end subroutine read_series

   !----- The value of field `name` at record irec (0 when the file does not carry it). ---------!
   pure integer(ik) function series_field(src, name) result(j)
      type(met_source_t), intent(in)  :: src
      character(len=*),   intent(in) :: name
      do j = 1_ik, size(src%series_name, kind=ik)
         if (trim(src%series_name(j)) == trim(name)) return
      end do
      j = 0_ik
   end function series_field

   !----- The file variable that carries each humidity form. -----------------------------------!
   pure function humidity_name(form) result(name)
      integer(ik), intent(in) :: form
      character(len=5) :: name
      select case (form)
      case (HUMIDITY_RHAIR) ; name = 'RHair'
      case (HUMIDITY_TDEW)  ; name = 'Tdew'
      case default          ; name = 'Qair'
      end select
   end function humidity_name

   function read_scalar(src, name, irec) result(val)
      type(met_source_t), intent(in)  :: src
      character(len=*),   intent(in) :: name
      integer(ik),        intent(in) :: irec
      real(wp)    :: val
      integer(ik) :: j
      j = series_field(src, name)
      if (j == 0_ik) then
         write(*,'(2a)') ' met_driver: the forcing file has no variable ', trim(name)
         error stop 'met_driver: a required forcing variable is missing'
      end if
      call check_in_series(src, irec)
      val = src%series(irec, j)
   end function read_scalar

   !----- A record outside the range read at open was not prefetched: a programming error. ------!
   subroutine check_in_series(src, irec)
      type(met_source_t), intent(in)  :: src
      integer(ik),        intent(in) :: irec
      if (irec < lbound(src%series, 1) .or. irec > ubound(src%series, 1)) then
         write(*,'(a,i0,a,i0,a,i0)') ' met_driver: forcing record ', irec, ' outside the range read, ',  &
               lbound(src%series, 1), '..', ubound(src%series, 1)
         error stop 'met_driver: a forcing record was not read before the step (internal error)'
      end if
   end subroutine check_in_series

   !----- `name` if the file carries it, else the supplied default (e.g. LWdown under synthesis). -!
   function read_scalar_default(src, name, irec, default) result(val)
      type(met_source_t), intent(in)  :: src
      character(len=*),   intent(in) :: name
      integer(ik),        intent(in) :: irec
      real(wp),           intent(in) :: default
      real(wp)    :: val
      integer(ik) :: j
      j = series_field(src, name)
      if (j == 0_ik) then
         val = default
      else
         call check_in_series(src, irec)
         val = src%series(irec, j)
      end if
   end function read_scalar_default

   !----- MEDS never gap-fills: a NaN in a required field halts the run (design §5.5). ---------!
   subroutine assert_finite(x, name, irec, grid)
      real(wp),         intent(in) :: x
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: irec, grid
      if (x /= x) then                                   ! NaN test (no ieee dependency)
         write(*,'(4a,i0,a,i0)') 'met_driver: missing/NaN value in required field ', trim(name), &
               ' at time record ', ' ', irec, ', grid ', grid
         error stop 'met_driver: forcing has a data gap (MEDS does not gap-fill; fix the source upstream)'
      end if
   end subroutine assert_finite

   !----- Parse "seconds since YYYY-MM-DD HH:MM:SS" -> base_time. ------------------------------!
   subroutine parse_time_units(units, base_time, ok)
      character(len=*),  intent(in)  :: units
      type(meds_time_t), intent(out) :: base_time
      logical,           intent(out) :: ok
      integer :: idx
      idx = index(units, 'since')
      if (idx == 0) then ; ok = .false. ; return ; end if
      call time_from_string(adjustl(units(idx + 5:)), base_time, ok)
   end subroutine parse_time_units


   !---------------------------------------------------------------------------------------!
   ! validate_file_against_config -- the forcing file describes itself; the config also describes  !
   ! it. When they disagree, the run is wrong in a way no later check can see, so stop here (#185). !
   !                                                                                          !
   ! THE SPACING CHECK IS THE IMPORTANT ONE. dt_forcing is taken verbatim from the config and used  !
   ! to place interval midpoints, disaggregate shortwave and bracket the recycle seam; nothing ever  !
   ! compared it against the file. A config saying 3600 s against a half-hourly file mis-times every !
   ! one of those. Checking the ACTUAL record spacing is strictly stronger than checking the         !
   ! `timestep_seconds` attribute, which is why that attribute stays provenance and this does the     !
   ! work.                                                                                            !
   !                                                                                          !
   ! The two text attributes are checked when PRESENT and skipped when absent, so a file written     !
   ! before the prep script emitted them still loads. `elevation(grid)` stays unread on purpose:     !
   ! site elevation is a [site] property, and the variable is provenance about the source grid.      !
   !---------------------------------------------------------------------------------------!
   subroutine validate_file_against_config(src, ncid, vstat)
      type(met_source_t), intent(in)   :: src
      integer(c_int),     intent(in)  :: ncid
      integer(ik),        intent(out) :: vstat
      character(len=64) :: attr
      real(wp)          :: dt_file
      integer(c_int)    :: st
      integer(ik)       :: i

      vstat = MET_OK

      !----- (a) record spacing. Uniform by construction for every source MEDS reads, so the first  !
      !      interval is the file's cadence; a ragged axis is a different (unsupported) thing and    !
      !      shows up as a mismatch on whichever interval differs. ------------------------------------!
      if (src%nrec >= 2_ik) then
         dt_file = src%time_sec(2) - src%time_sec(1)
         if (abs(dt_file - src%dt_forcing) > REC_MATCH_TOL) then
            write(*,'(a,f12.3,a)') ' met_open: [forcing].dt_forcing = ', src%dt_forcing, ' s'
            write(*,'(a,f12.3,a)') '   but the file record spacing is ', dt_file, ' s.'
            write(*,'(a)')         '   dt_forcing places interval midpoints, disaggregates shortwave and'
            write(*,'(a)')         '   brackets the recycle seam. A wrong value mis-times all three.'
            vstat = MET_ERR_DT_MISMATCH ; return
         end if
         do i = 3_ik, src%nrec
            if (abs((src%time_sec(i) - src%time_sec(i-1)) - dt_file) > REC_MATCH_TOL) then
               write(*,'(a,i0,a)') ' met_open: forcing record spacing changes at record ', i, '.'
               write(*,'(a)')      '   MEDS assumes a uniform time axis.'
               vstat = MET_ERR_AXIS_NOT_UNIFORM ; return
            end if
         end do
      end if

      !----- (b) avg_convention: the file says which end of the interval its flux means belong to.  !
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'avg_convention', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) /= trim(avg_convention_name(src%fcfg%avg_convention))) then
            write(*,'(4a)') ' met_open: file avg_convention = "', trim(attr),                     &
                            '", config = "', trim(avg_convention_name(src%fcfg%avg_convention))//'"'
            vstat = MET_ERR_ATTR_MISMATCH ; return
         end if
      end if

      !----- (c) sw_input_kind: "total" needs a partition, "components" must not be partitioned.   !
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) == 'total' .and. src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
            write(*,'(a)') ' met_open: the file carries TOTAL shortwave (sw_input_kind = "total")'
            write(*,'(a)') '   but [forcing].sw_partition = "passthrough", which expects the four'
            write(*,'(a)') '   component streams. The run would read fields the file does not have.'
            vstat = MET_ERR_ATTR_MISMATCH ; return
         end if
         if (trim(attr) == 'components' .and. src%fcfg%sw_partition /= SWPART_PASSTHROUGH) then
            write(*,'(a)') ' met_open: the file already carries the four shortwave components'
            write(*,'(a)') '   (sw_input_kind = "components") but [forcing].sw_partition asks for a'
            write(*,'(a)') '   split. Partitioning an already-split stream double-counts the beam.'
            vstat = MET_ERR_ATTR_MISMATCH ; return
         end if
      end if

      !----- (d) The clock. Every forcing file is in UTC (MEDS_FLUX_TOWER_FORCING_PLAN.md D1) and   !
      !      says so. A local-time file read as UTC keeps its daily totals and moves its sun by the   !
      !      offset, which nothing downstream notices, so the attribute is required. The archive is   !
      !      UTC by construction and is not asked. -------------------------------------------------!
      if (src%backend == MET_BACKEND_ED_DEFAULT) then
         st = nc_get_att_text_f(ncid, NC_GLOBAL, 'time_zone', attr)
         if (st /= NC_NOERR) attr = '(absent)'
         if (trim(attr) /= 'UTC') then
            write(*,'(3a)') ' met_open: the forcing file''s time_zone attribute is "', trim(attr), '".'
            write(*,'(a)')  '   MEDS runs in UTC: build the file on a UTC clock and set time_zone = "UTC".'
            vstat = MET_ERR_NOT_UTC ; return
         end if
      end if

      !----- (e) The heights the file was measured at, when it states them, against the heights the  !
      !      config moves the forcing from (docs/science/forcing.md §8). A disagreement would move   !
      !      every sample from the wrong height. -------------------------------------------------!
      call check_height_attribute(ncid, 'tq_height_m',        src%fcfg%tq_height,   vstat)
      if (vstat /= MET_OK) return
      call check_height_attribute(ncid, 'wind_height_m',      src%fcfg%wind_height, vstat)
      if (vstat /= MET_OK) return
      call check_height_attribute(ncid, 'wind_meas_height_m', src%fcfg%wind_height, vstat)
      if (vstat /= MET_OK) return
      st = nc_get_att_text_f(ncid, NC_GLOBAL, 'height_above', attr)
      if (st == NC_NOERR .and. len_trim(attr) > 0) then
         if (trim(attr) /= trim(height_above_name(src%fcfg%height_above))) then
            write(*,'(4a)') ' met_open: file height_above = "', trim(attr), '", config = "',        &
                            trim(height_above_name(src%fcfg%height_above))//'"'
            vstat = MET_ERR_ATTR_MISMATCH ; return
         end if
      end if

      !----- (f) CO2air (#184): CO2 comes from [forcing].co2_source, never the met file, so a file  !
      !      carrying it would be read by nothing. Rejected rather than ignored. -------------------!
      block
         integer(c_int) :: vid
         if (nc_inq_varid_f(ncid, 'CO2air', vid) == NC_NOERR) then
            write(*,'(a)') ' met_open: the forcing file carries CO2air, but CO2 comes from'
            write(*,'(a)') '   [forcing].co2_source ("const" or "file"), never from the met file.'
            write(*,'(a)') '   Remove it (ncks -x -v CO2air in.nc out.nc) or rebuild the file with'
            write(*,'(a)') '   make_forcing_file.py, which no longer writes it.'
            vstat = MET_ERR_CO2_IN_MET_FILE ; return
         end if
      end block
   end subroutine validate_file_against_config

   !----- A file's stated measurement height (a global attribute in m), when present, against the   !
   !      configured one.                                                                             !
   subroutine check_height_attribute(ncid, name, configured, vstat)
      integer(c_int),   intent(in)  :: ncid
      character(len=*), intent(in)  :: name
      real(wp),         intent(in)  :: configured
      integer(ik),      intent(out) :: vstat
      real(c_double) :: stated
      vstat = MET_OK
      if (nc_get_att_double_f(ncid, NC_GLOBAL, name, stated) /= NC_NOERR) return
      if (abs(real(stated, wp) - configured) > HEIGHT_MATCH_TOL) then
         write(*,'(3a,f0.3,a,f0.3,a)') ' met_open: the forcing file states ', name, ' = ', stated,     &
                                       ' m, but [forcing] says ', configured, ' m.'
         vstat = MET_ERR_ATTR_MISMATCH
      end if
   end subroutine check_height_attribute

   !----- The humidity form of a MEDS forcing file: exactly one of its three variables (D2). A file  !
   !      with none cannot drive the canopy air; a file with two leaves the reader to pick one, and   !
   !      the two may disagree through the provider's own saturation curve. ------------------------!
   subroutine detect_humidity(src, ncid, stat)
      type(met_source_t), intent(inout) :: src
      integer(c_int),     intent(in)    :: ncid
      integer(ik),        intent(out)   :: stat
      integer(ik), parameter :: FORMS(3) = [HUMIDITY_QAIR, HUMIDITY_RHAIR, HUMIDITY_TDEW]
      integer(c_int) :: vid
      integer(ik)    :: i, n
      stat = MET_OK ; n = 0_ik
      do i = 1_ik, size(FORMS, kind=ik)
         if (nc_inq_varid_f(ncid, trim(humidity_name(FORMS(i))), vid) == NC_NOERR) then
            n = n + 1_ik ; src%humidity = FORMS(i)
         end if
      end do
      if (n /= 1_ik) then
         write(*,'(a,i0,a)') ' met_open: the forcing file carries ', n,                                &
                             ' of the humidity variables RHair, Tdew and Qair; it must carry exactly one.'
         write(*,'(a)')      '   Store the humidity the source measured (RHair for a tower, Tdew for'
         write(*,'(a)')      '   ERA5-Land); MEDS converts it with its own saturation curve.'
         stat = MET_ERR_HUMIDITY
      end if
   end subroutine detect_humidity

   pure function height_above_name(code) result(nm)
      integer(ik), intent(in) :: code
      character(len=10) :: nm
      nm = 'zero_plane'
      if (code == HEIGHT_ABOVE_GROUND) nm = 'ground'
   end function height_above_name

   !----- The config code's own spelling, so the mismatch message quotes both sides in one vocabulary. !
   pure function avg_convention_name(code) result(nm)
      integer(ik), intent(in) :: code
      character(len=8) :: nm
      select case (code)
      case (METAVG_BEGIN) ; nm = 'begin'
      case default        ; nm = 'end'
      end select
   end function avg_convention_name

end module meds_met_driver
