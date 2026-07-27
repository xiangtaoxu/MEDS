!==========================================================================================!
! meds_met_driver -- the meteorological-forcing READER (design MEDS_FORCING_DESIGN.md section  !
! 4). Opens the MEDS multi-grid forcing NetCDF (or the no-file CONST backend), holds the two      !
! records bracketing the model time at this polygon's grid_index (met_driver_t, the ONLY mutable   !
! forcing state), slides that window as the model marches, and produces an instantaneous            !
! met_forcing_t via the pure disaggregation kernels. Owns file I/O (via meds_netcdf_c) -- the ED2    !
! cgrid%metinput analogue, threaded by the driver, never a global.                                    !
!                                                                                          !
! MEDS NEVER gap-fills (design §5.5): a missing/NaN required value is a HARD ERROR. netCDF is a       !
! hard dependency, so the reader is always the real thing (no stub).                                   !
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
   use meds_forcing_config, only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_NETCDF,     &
                                   METAVG_END, METAVG_BEGIN, SWPART_PASSTHROUGH,                &
                                   CLAMP_ERROR, INTERP_LINEAR, INTERP_STEP,                     &
                                   GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   use meds_forcing_types,  only : met_forcing_t, met_record_t, met_driver_t
   use meds_forcing_kernels, only : interpolate_forcing, interpolate_wind_energy,              &
                                   met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,   &
                                   partition_shortwave, precip_phase, nearest_grid_index,       &
                                   great_circle_distance, wind_log_profile,                    &
                                   lapse_air_temperature, lapse_pressure
   use meds_netcdf_c,       only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f,                  &
                                   nc_get_att_text_f, nc_get_vara_double, nc_close, nc_check,   &
                                   NC_NOERR, NC_NOWRITE
   implicit none
   private

   public :: met_open, met_advance, met_instant, met_close
   public :: MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD,                &
             MET_ERR_WINDOW_NOT_COVERED

   real(wp), parameter :: U_MIN     = 0.1_wp     !< [m/s] wind floor (M-O similarity stability)
   integer(ik), parameter :: N_COSZ_SUB = 10_ik  !< sub-samples per forcing interval for <cosz>_win

   !----- met_open status codes (see the `stat` argument). ---------------------------------!
   integer(ik), parameter :: MET_OK                          = 0_ik
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_WHOLE_YEARS  = 1_ik   !< recycle_end - recycle_start /= N years
   integer(ik), parameter :: MET_ERR_START_NOT_A_RECORD      = 2_ik   !< recycle_start is not a record stamp
   integer(ik), parameter :: MET_ERR_WINDOW_NOT_COVERED      = 3_ik   !< file stops short of recycle_end

   !----- Upper bound on the declared recycle window, in whole calendar years (search bound only). !
   integer(ik), parameter :: MAX_RECYCLE_YEARS = 200_ik
   !----- Tolerance for "this record stamp IS that instant" [s]. The time axis is float seconds,   !
   !      so an exact == would be brittle; sub-second slack is far below any real forcing dt. -------!
   real(wp), parameter :: REC_MATCH_TOL = 0.5_wp

contains

   !=======================================================================================!
   !  OPEN: CONST -> reference climate; NETCDF -> read the grid/time dims, the time axis, the    !
   !  base-time anchor from the `time:units` attribute, and load records #1-2 at grid_index.     !
   !=======================================================================================!
   !  `stat` (optional) reports a rejected recycle window instead of halting, so the validation    !
   !  is exercisable from a test (CLAUDE.md: "errors via error stop / status codes ... so failures  !
   !  are catchable in tests"). Absent -> a rejection is a hard error, which is what a production    !
   !  run wants: a forcing file that does not match its declared window must never be guessed at.   !
   subroutine met_open(drv, fcfg, stat)
      type(met_driver_t),     intent(inout) :: drv
      type(forcing_config_t), intent(in)    :: fcfg
      integer(ik), optional,  intent(out)   :: stat
      integer(c_int)    :: st, ncid
      integer(c_size_t) :: dlen
      character(len=256):: units
      logical           :: ok
      integer(ik)       :: vstat

      if (present(stat)) stat = MET_OK

      drv%fcfg       = fcfg
      drv%backend    = fcfg%backend
      drv%grid_index = fcfg%grid_index
      drv%dt_forcing = fcfg%dt_forcing
      drv%irec_prev  = 0_ik

      if (fcfg%backend == MET_BACKEND_CONST) then
         drv%ngrid = 1_ik ; drv%nrec = 0_ik
         drv%rec_prev = met_record_t() ; drv%rec_next = met_record_t()
         return
      end if

      !----- NetCDF backend. -----------------------------------------------------------------!
      st = nc_open_f(trim(fcfg%path), NC_NOWRITE, ncid)
      call nc_check(st, 'met_open: nc_open '//trim(fcfg%path))
      drv%ncid = int(ncid, ik)

      st = nc_inq_dimlen_f(ncid, 'grid', dlen) ; call nc_check(st, 'met_open: grid dim')
      drv%ngrid = int(dlen, ik)
      !----- bind this polygon to a grid slice: explicit index, or nearest [site] lat/lon (§4.1). !
      if (fcfg%grid_match == GRIDMATCH_NEAREST) then
         call resolve_grid_index(drv)                        ! sets drv%grid_index in 1..ngrid
      else                                                   ! GRIDMATCH_EXPLICIT
         if (fcfg%grid_index < 1_ik .or. fcfg%grid_index > drv%ngrid) then
            write(*,'(a,i0,a,i0)') 'met_open: grid_index ', fcfg%grid_index, ' out of range 1..', drv%ngrid
            error stop 'met_open: grid_index out of range'
         end if
      end if

      st = nc_inq_dimlen_f(ncid, 'time', dlen) ; call nc_check(st, 'met_open: time dim')
      drv%nrec = int(dlen, ik)
      if (drv%nrec < 1_ik) error stop 'met_open: forcing file has no time records'

      !----- base time from "seconds since <base>" units attribute. --------------------------!
      block
         integer(c_int) :: vid
         st = nc_inq_varid_f(ncid, 'time', vid) ; call nc_check(st, 'met_open: time varid')
         st = nc_get_att_text_f(ncid, vid, 'units', units) ; call nc_check(st, 'met_open: time:units')
         call parse_time_units(units, drv%base_time, ok)
         if (.not. ok) error stop 'met_open: could not parse time:units "seconds since <base>"'
      end block

      !----- cache the whole time axis (seconds since base_time). ----------------------------!
      call read_time_axis(drv)

      !----- V2/V3: check the DECLARED recycle window against this file's actual record stamps.  !
      !      Cannot live in the config sanity check (V1, whole-year span) -- it needs the time      !
      !      axis, which is only known here. -------------------------------------------------------!
      call validate_recycle_window(drv, vstat)
      if (vstat /= MET_OK) then
         !----- Release the file before bailing out. A `stat` return is a normal (if unhappy) exit  !
         !      the caller may retry from, and an un-closed netCDF handle keeps an HDF5 lock on the  !
         !      path -- a later create/open of the same file then fails with a bare "permission      !
         !      denied" that points nowhere near the real cause. ------------------------------------!
         call met_close(drv)
         if (present(stat)) then ; stat = vstat ; return ; end if
         error stop 'met_open: forcing file does not match the declared [forcing] recycle window'
      end if

      call load_bracket(drv, 1_ik)      ! records #1-2
   end subroutine met_open

   !=======================================================================================!
   !  Validate the DECLARED recycle window ([forcing].recycle_start/recycle_end) against the    !
   !  file that was just opened. MEDS does NOT infer the window: it is told where the cycle       !
   !  starts and how long it is, and this routine's only job is to confirm the file agrees.       !
   !                                                                                              !
   !  The predecessor of this routine did the opposite -- it classified the file and, when the      !
   !  file did not look Jan-1-aligned, silently dropped to an absolute-seconds span-wrap. For the   !
   !  real ERA5-Land record (first stamp 01:00 under the end-of-interval convention) that fallback   !
   !  wrapped on a 366 d 22 h span, shifting hour-of-day on EVERY wrap: a 29-yr run ended up reading  !
   !  late May at a ~10 h offset while its daily-mean shortwave stayed correct, so the slow            !
   !  demography looked healthy and nothing surfaced the problem for 30 simulated years.                !
   !=======================================================================================!
   subroutine validate_recycle_window(drv, stat)
      type(met_driver_t), intent(inout) :: drv
      integer(ik),        intent(out)   :: stat
      real(wp)    :: start_sec, end_sec
      integer(ik) :: i
      logical     :: covered

      stat = MET_OK
      drv%cycle_anchor = drv%fcfg%recycle_start
      drv%n_cycle_years = 0_ik ; drv%irec_cycle_first = 0_ik ; drv%irec_cycle_last = 0_ik
      if (.not. drv%fcfg%recycle) return

      !----- Whole-year span. Re-checked here (not only in validate_config) so a driver built     !
      !      directly from a forcing_config_t -- as the unit tests do -- cannot bypass it. --------!
      drv%n_cycle_years = whole_years_between(drv%fcfg%recycle_start, drv%fcfg%recycle_end,     &
                                              MAX_RECYCLE_YEARS)
      if (drv%n_cycle_years < 1_ik) then
         write(*,'(4a)') ' met_open: recycle window ', time_to_string(drv%fcfg%recycle_start),   &
                         ' .. ', time_to_string(drv%fcfg%recycle_end)
         write(*,'(a)')  '   is not an exact whole number of calendar years.'
         stat = MET_ERR_WINDOW_NOT_WHOLE_YEARS ; return
      end if

      start_sec = seconds_between(drv%base_time, drv%fcfg%recycle_start)
      end_sec   = seconds_between(drv%base_time, drv%fcfg%recycle_end)

      !----- V2: recycle_start must land EXACTLY on a record stamp. This is the check that stops  !
      !      MEDS guessing the start of the day: a config saying 00:00:00 against an ERA5-Land      !
      !      file stamped 01:00:00 is reported, with both values, rather than quietly re-derived.   !
      do i = 1_ik, drv%nrec
         if (abs(drv%time_sec(i) - start_sec) <= REC_MATCH_TOL) then
            drv%irec_cycle_first = i ; exit
         end if
      end do
      if (drv%irec_cycle_first == 0_ik) then
         write(*,'(2a)') ' met_open: forcing.recycle_start = ', time_to_string(drv%fcfg%recycle_start)
         write(*,'(a)')  '   does not match any record in the forcing file. The file record stamps run'
         write(*,'(4a)') '   from ', time_to_string(time_advance_seconds(drv%base_time, drv%time_sec(1))), &
                         ' to ',    time_to_string(time_advance_seconds(drv%base_time, drv%time_sec(drv%nrec)))
         write(*,'(a)')  '   Set recycle_start to an actual record stamp -- MEDS does not guess it. Note'
         write(*,'(a)')  '   avg_convention="end" files (ERA5-Land) stamp the END of each interval, so the'
         write(*,'(a)')  '   first record of a calendar year is 01:00:00 for hourly data, not 00:00:00.'
         stat = MET_ERR_START_NOT_A_RECORD ; return
      end if

      !----- V3: the file must cover the whole window, i.e. hold a record in the LAST interval    !
      !      below recycle_end. Without that the cycle has a hole at its far edge and the seam      !
      !      bracket below would interpolate across it. --------------------------------------------!
      do i = drv%nrec, drv%irec_cycle_first, -1_ik
         if (drv%time_sec(i) < end_sec - REC_MATCH_TOL) then ; drv%irec_cycle_last = i ; exit ; end if
      end do
      !----- NESTED, not a single .or.: Fortran does not guarantee short-circuit evaluation, so a  !
      !      combined test would index time_sec(irec_cycle_last) even when the loop above left it 0. !
      covered = drv%irec_cycle_last > drv%irec_cycle_first
      if (covered) covered = drv%time_sec(drv%irec_cycle_last) >= end_sec - drv%dt_forcing - REC_MATCH_TOL
      if (.not. covered) then
         write(*,'(4a)') ' met_open: the forcing file does not cover the declared recycle window ', &
                         time_to_string(drv%fcfg%recycle_start), ' .. ',                         &
                         time_to_string(drv%fcfg%recycle_end)
         write(*,'(2a)') '   last record in the file: ',                                         &
                         time_to_string(time_advance_seconds(drv%base_time, drv%time_sec(drv%nrec)))
         stat = MET_ERR_WINDOW_NOT_COVERED ; return
      end if

      write(*,'(a,i0,4a)') ' force : recycling ', drv%n_cycle_years, ' calendar year(s), ',      &
            time_to_string(drv%fcfg%recycle_start), ' .. ', time_to_string(drv%fcfg%recycle_end)
   end subroutine validate_recycle_window

   !=======================================================================================!
   !  ADVANCE: slide the window so rec_prev%when <= now < rec_next%when (at this grid_index).    !
   !  Handles start-before-base_time (clamp/error) and EOF (recycle by whole file spans).        !
   !=======================================================================================!
   subroutine met_advance(drv, now)
      type(met_driver_t), intent(inout) :: drv
      type(meds_time_t),  intent(in)    :: now
      real(wp)    :: now_sec
      integer(ik) :: i

      if (drv%backend == MET_BACKEND_CONST) return

      !----- the SAME effective (recycle-mapped) seconds on the FILE axis met_instant uses. ---!
      now_sec = file_lookup_sec(drv, now)

      !----- start before the first record (never reached once recycle-wrapped into span). ---!
      !      Hold: load records 1-2 and let met_instant clamp w_next=0 for now < t1 -- do NOT     !
      !      overwrite rec_next (that stale-copy bug suppressed the reload once now reached t1).  !
      if (now_sec < drv%time_sec(1)) then
         if (drv%fcfg%start_clamp == CLAMP_ERROR) then
            error stop 'met_advance: model start precedes the first forcing record (start_clamp=error)'
         end if
         if (drv%irec_prev /= 1_ik) call load_bracket(drv, 1_ik)
         return
      end if

      !----- CYCLE-BOUNDARY seam: the final interval of the declared window wraps back to its own  !
      !      first record, so the cycle is exactly periodic. The seam is the WINDOW's edge, not the  !
      !      file's -- a window may be a sub-range of a longer file, and trailing records past        !
      !      recycle_end belong to the next file year, not to this cycle. -----------------------------!
      !      NESTED, not one .and.: Fortran does not guarantee short-circuit evaluation, so a
      !      combined test would index time_sec(irec_cycle_last) on a non-recycling driver, where
      !      it is still 0.
      if (drv%n_cycle_years >= 1_ik .and. drv%irec_cycle_last > drv%irec_cycle_first) then
         if (now_sec >= drv%time_sec(drv%irec_cycle_last)) then
            if (.not. drv%at_wrap_seam) call load_wrap_bracket(drv)
            return
         end if
      end if

      !----- non-recycle run past the end: clamp to the last interval (recycle wraps above). --!
      if (.not. drv%fcfg%recycle .and. drv%nrec >= 2_ik .and. now_sec >= drv%time_sec(drv%nrec)) then
         if (drv%irec_prev /= drv%nrec - 1_ik) call load_bracket(drv, drv%nrec - 1_ik)
         return
      end if

      !----- incremental cursor from the current bracket (cheap for a forward march). --------!
      i = max(1_ik, drv%irec_prev)
      do while (i < drv%nrec - 1_ik)
         if (drv%time_sec(i + 1_ik) <= now_sec) then ; i = i + 1_ik ; else ; exit ; end if
      end do
      do while (i > 1_ik)
         if (drv%time_sec(i) > now_sec) then ; i = i - 1_ik ; else ; exit ; end if
      end do
      if (i /= drv%irec_prev) call load_bracket(drv, i)
   end subroutine met_advance

   !=======================================================================================!
   !  INSTANT: interpolate/disaggregate the loaded window to the model instant `now`.           !
   !  State vars linear; wind energy-form; precip step-constant then phase-split; shortwave via    !
   !  the reciprocal-mean-cosz reconstruction of the interval-mean streams; cosz + rho_air derived. !
   !=======================================================================================!
   function met_instant(drv, now) result(met)
      type(met_driver_t), intent(in) :: drv
      type(meds_time_t),  intent(in) :: now
      type(met_forcing_t) :: met
      type(met_record_t)  :: mean_rec
      type(meds_time_t)   :: mws
      real(wp) :: now_sec, tprev, tnext, w_next, cosz_now, factor, win_start_sec, precip_total
      associate (f => drv%fcfg, p => drv%rec_prev, n => drv%rec_next)

      !----- solar zenith at `now` (file clock -> apparent solar seconds inside met_solar_cosz). !
      cosz_now = met_solar_cosz(now, seconds_into_day(now), f%latitude_deg, f%longitude_deg,   &
                                f%utc_offset_h, f%apply_solar_longitude)
      met%cosz = cosz_now

      if (drv%backend == MET_BACKEND_CONST) then                  ! reference climate held flat
         met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)
         return
      end if

      !----- interpolation weight within the loaded window (SAME effective seconds as advance, !
      !      so a recycle-mapped now interpolates within the recycled interval, not clamps to 1). !
      now_sec = file_lookup_sec(drv, now)
      if (drv%at_wrap_seam) then                    ! cycle-boundary: window's last rec -> its first
         tprev = drv%time_sec(drv%irec_cycle_last) ; tnext = tprev + drv%dt_forcing
      else
         tprev = drv%time_sec(drv%irec_prev)
         tnext = drv%time_sec(min(drv%irec_prev + 1_ik, drv%nrec))
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
      met%co2      = interpolate_forcing(INTERP_LINEAR, p%co2,      n%co2,      w_next)
      met%wind     = interpolate_wind_energy(p%wind, n%wind, w_next, U_MIN)

      !----- precip: step-constant total (never smeared), then phase-split. -----------------!
      precip_total = interpolate_forcing(INTERP_STEP, p%rainf, n%rainf, w_next)
      call precip_phase(precip_total, met%tair_k, met%rainf, met%snowf)

      !----- shortwave: the interval-mean streams of the interval CONTAINING now, disaggregated  !
      !      by cosz(now)/<cosz>_win. avg_convention=end -> the interval [prev,next] mean is        !
      !      rec_next's; begin -> rec_prev's.                                                        !
      select case (f%avg_convention)
      case (METAVG_BEGIN) ; mean_rec = p
      case default        ; mean_rec = n            ! METAVG_END (ERA5-Land) + fallback
      end select
      !----- reconstruction factor anchored on the MODEL window start (mws), so <cosz>_win aligns  !
      !      with cosz_now on the model calendar (identity = rec_prev%when when not recycling; under  !
      !      calendar recycling it follows the model sun, so the interval-mean identity still holds).  !
      mws           = time_advance_seconds(now, tprev - now_sec)   ! model instant at the window start
      win_start_sec = seconds_into_day(mws)
      factor = cosz_reconstruct_factor(mws, win_start_sec,                                       &
                                       (tnext - tprev) / real(N_COSZ_SUB, wp), tnext - tprev,   &
                                       f%latitude_deg, f%longitude_deg, f%utc_offset_h,         &
                                       f%apply_solar_longitude)
      met%par_beam    = disaggregate_sw(mean_rec%par_beam,    cosz_now, factor)
      met%par_diffuse = disaggregate_sw(mean_rec%par_diffuse, cosz_now, factor)
      met%nir_beam    = disaggregate_sw(mean_rec%nir_beam,    cosz_now, factor)
      met%nir_diffuse = disaggregate_sw(mean_rec%nir_diffuse, cosz_now, factor)

      met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)
      end associate
   end function met_instant

   !=======================================================================================!
   !  CLOSE.                                                                                     !
   !=======================================================================================!
   subroutine met_close(drv)
      type(met_driver_t), intent(inout) :: drv
      integer(c_int) :: st
      if (drv%backend == MET_BACKEND_NETCDF .and. drv%ncid >= 0_ik) then
         st = nc_close(int(drv%ncid, c_int)) ; call nc_check(st, 'met_close: nc_close')
         drv%ncid = -1_ik
      end if
      if (allocated(drv%time_sec)) deallocate(drv%time_sec)
   end subroutine met_close

   !----- The effective seconds-since-base on the FILE time axis, used by BOTH bracket selection  !
   !      (met_advance) and the interpolation weight (met_instant) so they stay consistent. Three   !
   !      regimes: (a) CALENDAR recycle (whole-year Jan-1 file) maps the model date to its file      !
   !      calendar year, preserving month/day/hour so day-of-year is exact across leap boundaries;   !
   !      (b) LEGACY absolute-seconds span-wrap (non-calendar recyclable file, e.g. an idealized      !
   !      diurnal repeat) once past EOF; (c) identity while the model date is within the file range.  !
   pure function file_lookup_sec(drv, now) result(s)
      type(met_driver_t), intent(in) :: drv
      type(meds_time_t),  intent(in) :: now
      real(wp) :: s
      if (drv%fcfg%recycle .and. drv%n_cycle_years >= 1_ik) then
         s = seconds_between(drv%base_time, recycle_model_to_file(drv, now))
      else
         s = seconds_between(drv%base_time, now)
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
   !      plain year substitution -- so a window that the old Jan-1-only classifier would have accepted !
   !      maps bit-for-bit as before. Without `off`, a mid-year window (or an end-of-interval file      !
   !      whose first stamp is 01:00) maps instants OUTSIDE the window it is supposed to cycle over.    !
   !                                                                                                   !
   !      LEAP DAY: Feb-29 -> Feb-28 when the target file year is non-leap (ED2 read_ol_file repeats    !
   !      Feb 28). A non-leap model year never asks for Feb-29, so the substitution is one-directional. !
   pure function recycle_model_to_file(drv, now) result(now_file)
      type(met_driver_t), intent(in) :: drv
      type(meds_time_t),  intent(in) :: now
      type(meds_time_t) :: now_file
      integer(ik) :: yf, off
      off = merge(0_ik, 1_ik, mdhms_at_or_after(now, drv%cycle_anchor))
      yf  = drv%cycle_anchor%year + off                                                          &
          + modulo(now%year - drv%cycle_anchor%year - off, drv%n_cycle_years)
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
   !----- Read the whole time coordinate (seconds since base_time) into drv%time_sec. ---------!
   subroutine read_time_axis(drv)
      type(met_driver_t), intent(inout) :: drv
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start1(1), count1(1)
      if (allocated(drv%time_sec)) deallocate(drv%time_sec)
      allocate(drv%time_sec(drv%nrec))
      st = nc_inq_varid_f(int(drv%ncid, c_int), 'time', vid) ; call nc_check(st, 'read_time_axis: time varid')
      start1(1) = 0_c_size_t ; count1(1) = int(drv%nrec, c_size_t)
      st = nc_get_vara_double(int(drv%ncid, c_int), vid, start1, count1, drv%time_sec)
      call nc_check(st, 'read_time_axis: time values')
   end subroutine read_time_axis

   !----- Load rec_prev = record(irec), rec_next = record(irec+1) (clamped at EOF). -----------!
   subroutine load_bracket(drv, irec)
      type(met_driver_t), intent(inout) :: drv
      integer(ik),        intent(in)    :: irec
      call read_record(drv, irec, drv%rec_prev)
      call read_record(drv, min(irec + 1_ik, drv%nrec), drv%rec_next)
      drv%irec_prev    = irec
      drv%at_wrap_seam = .false.
   end subroutine load_bracket

   !----- Cycle-boundary seam bracket: rec_prev = record nrec, rec_next = record 1, so the last  !
   !      dt interval interpolates across the wrap instead of clamping to the final record.        !
   subroutine load_wrap_bracket(drv)
      type(met_driver_t), intent(inout) :: drv
      call read_record(drv, drv%irec_cycle_last,  drv%rec_prev)
      call read_record(drv, drv%irec_cycle_first, drv%rec_next)
      drv%irec_prev    = drv%irec_cycle_last
      drv%at_wrap_seam = .true.
   end subroutine load_wrap_bracket

   !----- Resolve drv%grid_index by nearest great-circle distance from [site] lat/lon to the file's !
   !      latitude(grid)/longitude(grid) coordinate vectors (multi-polygon P2 subset, §4.1).        !
   subroutine resolve_grid_index(drv)
      type(met_driver_t), intent(inout) :: drv
      integer(c_int)    :: st, vlat, vlon
      integer(c_size_t) :: start1(1), count1(1)
      real(wp), allocatable :: lat_grid(:), lon_grid(:)
      allocate(lat_grid(drv%ngrid), lon_grid(drv%ngrid))
      start1(1) = 0_c_size_t ; count1(1) = int(drv%ngrid, c_size_t)
      st = nc_inq_varid_f(int(drv%ncid, c_int), 'latitude', vlat)  ; call nc_check(st, 'resolve_grid_index: latitude varid')
      st = nc_get_vara_double(int(drv%ncid, c_int), vlat, start1, count1, lat_grid)
      call nc_check(st, 'resolve_grid_index: latitude values')
      st = nc_inq_varid_f(int(drv%ncid, c_int), 'longitude', vlon) ; call nc_check(st, 'resolve_grid_index: longitude varid')
      st = nc_get_vara_double(int(drv%ncid, c_int), vlon, start1, count1, lon_grid)
      call nc_check(st, 'resolve_grid_index: longitude values')
      drv%grid_index = nearest_grid_index(drv%fcfg%longitude_deg, drv%fcfg%latitude_deg, lon_grid, lat_grid)
      write(*,'(a,i0,a,f8.3,a,f8.3,a)') 'met_open: grid_match=nearest resolved grid_index=', drv%grid_index, &
            ' (site lon=', drv%fcfg%longitude_deg, ' lat=', drv%fcfg%latitude_deg, ')'
      deallocate(lat_grid, lon_grid)
   end subroutine resolve_grid_index

   !----- Read one record at (time=irec, grid=grid_index); partition SW at ingest. ------------!
   subroutine read_record(drv, irec, rec)
      type(met_driver_t), intent(in)  :: drv
      integer(ik),        intent(in)  :: irec
      type(met_record_t), intent(out) :: rec
      real(wp) :: sw_total, cosz_mid, mid_sec
      rec%when   = time_advance_seconds(drv%base_time, drv%time_sec(irec))
      rec%tair_k   = read_scalar(drv, 'Tair',  irec)
      rec%qair     = read_scalar(drv, 'Qair',  irec)
      rec%psurf_pa = read_scalar(drv, 'PSurf', irec)
      rec%wind     = read_scalar(drv, 'Wind',  irec)
      rec%rainf    = read_scalar(drv, 'Rainf', irec)                 ! total precip rate [kg/m2/s]
      rec%lwdown   = read_scalar(drv, 'LWdown', irec)
      rec%co2      = read_scalar_default(drv, 'CO2air', irec, drv%fcfg%co2_const)
      call assert_finite(rec%tair_k, 'Tair', irec, drv%grid_index)
      call assert_finite(rec%qair, 'Qair', irec, drv%grid_index)
      call assert_finite(rec%psurf_pa, 'PSurf', irec, drv%grid_index)
      call assert_finite(rec%wind, 'Wind', irec, drv%grid_index)
      call assert_finite(rec%rainf, 'Rainf', irec, drv%grid_index)
      call assert_finite(rec%lwdown, 'LWdown', irec, drv%grid_index)

      !----- wind-height + elevation-lapse corrections at ingest (opt-in; both default OFF, so    !
      !      CONST and un-flagged runs are untouched). Wind is a per-record height rescale (commutes !
      !      with the downstream energy-form interpolation). Lapse moves T then P (hydrostatic,       !
      !      consistent) from the grid-cell elevation to the site; qair is held (rho re-derived).     !
      if (drv%fcfg%apply_wind_profile) then
         rec%wind = wind_log_profile(rec%wind, drv%fcfg%wind_meas_height,                        &
                                     drv%fcfg%reference_height, drv%fcfg%wind_roughness_z0)
      end if
      if (drv%fcfg%apply_elevation_lapse) then
         block
            real(wp) :: dz, t_grid, p_grid
            dz     = drv%fcfg%elevation_m - drv%fcfg%grid_elevation_m       ! + when site is higher
            t_grid = rec%tair_k ; p_grid = rec%psurf_pa                     ! capture BEFORE overwrite
            rec%psurf_pa = lapse_pressure(p_grid, t_grid, dz, drv%fcfg%lapse_rate_tair)
            rec%tair_k   = lapse_air_temperature(t_grid, dz, drv%fcfg%lapse_rate_tair)
         end block
      end if

      if (drv%fcfg%sw_partition == SWPART_PASSTHROUGH) then
         rec%par_beam    = read_scalar(drv, 'SWdown_par_beam',    irec)
         rec%par_diffuse = read_scalar(drv, 'SWdown_par_diffuse', irec)
         rec%nir_beam    = read_scalar(drv, 'SWdown_nir_beam',    irec)
         rec%nir_diffuse = read_scalar(drv, 'SWdown_nir_diffuse', irec)
         call assert_finite(rec%par_beam,    'SWdown_par_beam',    irec, drv%grid_index)   ! required source
         call assert_finite(rec%par_diffuse, 'SWdown_par_diffuse', irec, drv%grid_index)   ! fields -> no gap-fill
         call assert_finite(rec%nir_beam,    'SWdown_nir_beam',    irec, drv%grid_index)
         call assert_finite(rec%nir_diffuse, 'SWdown_nir_diffuse', irec, drv%grid_index)
      else
         sw_total = read_scalar(drv, 'SWdown', irec)
         call assert_finite(sw_total, 'SWdown', irec, drv%grid_index)
         !----- interval-mean cosz for the partition (avg_convention=end -> midpoint = when - dt/2). !
         mid_sec  = seconds_into_day(rec%when)
         if (drv%fcfg%avg_convention == METAVG_END)   mid_sec = mid_sec - 0.5_wp * drv%dt_forcing
         if (drv%fcfg%avg_convention == METAVG_BEGIN)  mid_sec = mid_sec + 0.5_wp * drv%dt_forcing
         cosz_mid = met_solar_cosz(rec%when, mid_sec, drv%fcfg%latitude_deg,                    &
                                   drv%fcfg%longitude_deg, drv%fcfg%utc_offset_h,               &
                                   drv%fcfg%apply_solar_longitude)
         call partition_shortwave(sw_total, cosz_mid, rec%psurf_pa, drv%fcfg%sw_partition,       &
                                  rec%par_beam, rec%par_diffuse, rec%nir_beam, rec%nir_diffuse)
      end if
   end subroutine read_record

   !----- Read a single value of `name` at (time=irec-1, grid=grid_index-1) (0-based C index). !
   function read_scalar(drv, name, irec) result(val)
      type(met_driver_t), intent(in) :: drv
      character(len=*),   intent(in) :: name
      integer(ik),        intent(in) :: irec
      real(wp) :: val
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start2(2), count2(2)
      real(c_double)    :: buf(1)
      st = nc_inq_varid_f(int(drv%ncid, c_int), name, vid)
      call nc_check(st, 'read_scalar: varid '//trim(name))
      start2 = [int(irec - 1_ik, c_size_t), int(drv%grid_index - 1_ik, c_size_t)]   ! [time, grid]
      count2 = [1_c_size_t, 1_c_size_t]
      st = nc_get_vara_double(int(drv%ncid, c_int), vid, start2, count2, buf)
      call nc_check(st, 'read_scalar: get '//trim(name))
      val = buf(1)
   end function read_scalar

   !----- Read `name` if present, else return the supplied default (e.g. CO2air absent). -------!
   function read_scalar_default(drv, name, irec, default) result(val)
      type(met_driver_t), intent(in) :: drv
      character(len=*),   intent(in) :: name
      integer(ik),        intent(in) :: irec
      real(wp),           intent(in) :: default
      real(wp) :: val
      integer(c_int) :: st, vid
      st = nc_inq_varid_f(int(drv%ncid, c_int), name, vid)
      if (st /= NC_NOERR) then
         val = default
      else
         val = read_scalar(drv, name, irec)
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

end module meds_met_driver
