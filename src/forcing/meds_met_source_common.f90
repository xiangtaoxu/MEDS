! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_source_common -- what the two forcing FILE sources share (MEDS_FORCING_DESIGN.md   !
! sections 4 and 15): the reader's status codes, the file time axis with its declared recycle  !
! window, and the check of a file against the [forcing] config that claims to describe it.     !
!                                                                                          !
! A MEDS forcing file (meds_met_file_source) and the ED_ERA5land archive                        !
! (meds_met_archive_source) both lay their records on one axis of seconds since a base time,   !
! so bracketing, recycling and the seam are one code for both (meds_met_driver).              !
!==========================================================================================!
module meds_met_source_common
   use iso_c_binding, only : c_int, c_double
   use meds_kinds, only : wp, ik
   use meds_time, only : meds_time_t, time_advance_seconds, seconds_between, is_leap_year,        &
                         time_to_string, whole_years_between
   use meds_forcing_config, only : MET_BACKEND_ED_DEFAULT, METAVG_BEGIN, SWPART_PASSTHROUGH,      &
                                   HEIGHT_ABOVE_GROUND
   use meds_forcing_types, only : met_source_t, HUMIDITY_RHAIR, HUMIDITY_TDEW, FLD_QAIR,          &
                                  FLD_RHAIR, FLD_TDEW
   use meds_config, only : MAX_RECYCLE_YEARS   ! the config's bound: one definition
   use meds_netcdf_c, only : nc_inq_varid_f, nc_get_att_text_f, nc_get_att_double_f, NC_NOERR,    &
                             NC_GLOBAL
   implicit none
   private

   public :: MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD,                  &
             MET_ERR_WINDOW_NOT_COVERED, MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,           &
             MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED,   &
             MET_ERR_CO2_IN_MET_FILE, MET_ERR_NOT_UTC, MET_ERR_HUMIDITY
   public :: REC_MATCH_TOL, validate_recycle_window, file_lookup_sec, first_record_after,         &
             validate_file_against_config, humidity_field

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

   !----- Tolerance for "this record stamp IS that instant" [s]. The time axis is float seconds,   !
   !      so an exact == would be brittle; sub-second slack is far below any real forcing dt. -------!
   real(wp), parameter :: REC_MATCH_TOL = 0.5_wp
   !----- A file's declared measurement height agrees with the config within this [m]. -----------!
   real(wp), parameter :: HEIGHT_MATCH_TOL = 0.01_wp

contains

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

   !----- The file variable that carries each humidity form. -----------------------------------!
   pure integer(ik) function humidity_field(form) result(fld)
      integer(ik), intent(in) :: form
      select case (form)
      case (HUMIDITY_RHAIR) ; fld = FLD_RHAIR
      case (HUMIDITY_TDEW)  ; fld = FLD_TDEW
      case default          ; fld = FLD_QAIR
      end select
   end function humidity_field

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

end module meds_met_source_common
