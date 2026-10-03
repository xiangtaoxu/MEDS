! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_driver -- the meteorological-forcing READER (design MEDS_FORCING_DESIGN.md sections !
! 4 and 15). Opens the multi-grid forcing NetCDF (ED_default), the global ED_ERA5land archive, or the !
! no-file CONST backend into a met_source_t shared by every polygon of a run; each polygon has a   !
! met_cursor_t holding its cell, its location and the two records bracketing the model time,      !
! slides that window as the model marches, and produces an instantaneous met_forcing_t via the     !
! pure disaggregation kernels (MEDS_POLYGON_RUNTIME_PLAN.md §10.3) -- the ED2 cgrid%metinput       !
! analogue, threaded by the driver, never a global. Stepping (met_advance, met_instant) only reads  !
! the source.                                                                                       !
!                                                                                          !
! The file I/O is the two sources': meds_met_file_source (a MEDS forcing file) and                !
! meds_met_archive_source (the archive), each an open and a record's values as stored; what they   !
! share is meds_met_source_common. This module keeps what is the same for every backend: opening,  !
! the cursor, the stepping, and the one ingest (read_record) that checks, converts, lapses and    !
! partitions a record.                                                                              !
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
   use iso_c_binding, only : c_int
   use meds_kinds, only : wp, ik
   use meds_therm_lib, only : air_density
   use meds_time, only : meds_time_t, time_advance_seconds, seconds_into_day, time_to_string
   use meds_forcing_config, only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT,   &
                                   MET_BACKEND_ED_ERA5LAND, METAVG_END, METAVG_BEGIN,             &
                                   SWPART_PASSTHROUGH, CLAMP_ERROR, INTERP_LINEAR,                &
                                   LW_SYNTHESIZE, CO2_SOURCE_FILE
   use meds_forcing_types, only : met_forcing_t, met_record_t, met_source_t, met_cursor_t,        &
                                  met_cells_t, HUMIDITY_RHAIR, HUMIDITY_TDEW, MEDS_FIELD
   use meds_lapse_rate, only : lapse_air_temperature, lapse_pressure, monthly_lapse_rate,         &
                               lapse_specific_humidity, lapse_longwave
   use meds_co2_series, only : co2_series_read, co2_series_at, co2_series_covers,                 &
                               co2_series_end, co2_series_free
   use meds_forcing_kernels, only : interpolate_forcing, interpolate_wind_energy,                 &
                                    met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,     &
                                    partition_shortwave, precip_phase, clearness_index,           &
                                    synthesize_lwdown, dewpoint_to_specific_humidity,             &
                                    rh_to_specific_humidity
   use meds_netcdf_c, only : nc_close, nc_check
   use meds_met_source_common, only : MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS,                     &
                                      MET_ERR_START_NOT_A_RECORD, MET_ERR_WINDOW_NOT_COVERED,     &
                                      MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,              &
                                      MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE,   &
                                      MET_ERR_CO2_NOT_COVERED, MET_ERR_CO2_IN_MET_FILE,           &
                                      MET_ERR_NOT_UTC, MET_ERR_HUMIDITY, file_lookup_sec,         &
                                      humidity_field
   use meds_met_file_source, only : file_source_open, file_source_record
   use meds_met_archive_source, only : archive_source_open, archive_prefetch,                     &
                                       archive_source_record
   implicit none
   private

   public :: met_open, met_cursor_init, met_advance, met_instant, met_close, met_prefetch
   public :: MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS, MET_ERR_START_NOT_A_RECORD,                  &
             MET_ERR_WINDOW_NOT_COVERED, MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,           &
             MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE, MET_ERR_CO2_NOT_COVERED,   &
             MET_ERR_CO2_IN_MET_FILE, MET_ERR_NOT_UTC, MET_ERR_HUMIDITY

   integer(ik), parameter :: N_COSZ_SUB = 10_ik  !< sub-samples per forcing interval for <cosz>_win

contains

   !=======================================================================================!
   !  OPEN: CONST -> the reference climate; ED_DEFAULT -> a MEDS forcing file                 !
   !  (file_source_open); ED_ERA5land -> the archive (archive_source_open). CO2 first, the   !
   !  same way for every backend.                                                            !
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
      integer(ik)       :: vstat
      character(len=160):: why

      if (present(stat)) stat = MET_OK

      src%fcfg       = fcfg
      src%backend    = fcfg%backend
      src%grid_index = fcfg%grid_index
      src%dt_forcing = fcfg%dt_forcing
      src%has_wind_vector = .false.

      !----- CO2 first, before the backends branch, so every backend takes it the same way. ------!
      call open_co2(src, run_start, run_end, vstat)
      why = 'met_open: the CO2 file cannot drive this run (see the message above)'
      if (vstat == MET_OK) then
         select case (fcfg%backend)
         case (MET_BACKEND_CONST)                                  ! the reference climate: no file
            src%ngrid = 1_ik ; src%nrec = 0_ik
         case (MET_BACKEND_ED_ERA5LAND)
            call archive_source_open(src, vstat, why, run_start, run_end, cells)
         case default                                              ! MET_BACKEND_ED_DEFAULT
            call file_source_open(src, vstat, why, run_start, run_end)
         end select
      end if
      if (vstat /= MET_OK) then
         !----- Release the file before bailing out. A `stat` return is a normal (if unhappy) exit  !
         !      the caller may retry from, and an un-closed netCDF handle keeps an HDF5 lock on the  !
         !      path -- a later create/open of the same file then fails with a bare "permission      !
         !      denied" that points nowhere near the real cause. ------------------------------------!
         call met_close(src)
         if (present(stat)) then ; stat = vstat ; return ; end if
         !----- The reason is printed, then a fixed stop code: ifx 2026 garbles a stop code that is  !
         !      not a constant (stray bytes, a truncated text, exit status 128). -----------------!
         write(*,'(1x,a)') trim(why)
         error stop 'met_open: the forcing cannot drive this run (see the messages above)'
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
   !  phase, so met_advance never touches a file. A MEDS forcing file was read at open; the     !
   !  archive loads a month at a time (archive_prefetch).                                       !
   !=======================================================================================!
   subroutine met_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(meds_time_t),  intent(in)    :: step_start
      if (src%backend == MET_BACKEND_ED_ERA5LAND) call archive_prefetch(src, step_start)
   end subroutine met_prefetch

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
   pure function met_instant(src, cur, now) result(met)
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
         met%wind     = interpolate_wind_energy(p%wind, n%wind, w_next)   ! floored at aerodynamics.ubmin where used
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
      if (allocated(src%series)) deallocate(src%series)
      src%series_col = 0_ik
      call co2_series_free(src%co2)
      src%buffer%year = 0_ik ; src%carry_rec = 0_ik ; src%n_head = 0_ik ; src%n_loads = 0_ik
   end subroutine met_close

   !----- Load rec_prev = record(irec), rec_next = record(irec+1) (clamped at EOF). When the      !
   !      bracket slides by one record, as it does once a forcing interval, the old upper record is !
   !      the new lower one and only the upper is read. ------------------------------------------!
   subroutine load_bracket(src, cur, irec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      integer(ik),        intent(in)    :: irec
      if (cur%irec_prev >= 1_ik .and. irec == cur%irec_prev + 1_ik .and. .not. cur%at_wrap_seam) then
         cur%rec_prev = cur%rec_next
      else
         call read_record(src, cur, irec, cur%rec_prev)
      end if
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

   !----- Read one record at (time=irec, this polygon's cell); partition SW at ingest. --------!
   subroutine read_record(src, cur, irec, rec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(out)   :: rec
      real(wp)    :: sw_total, cosz_mid, mid_sec, humidity_value
      rec%when   = time_advance_seconds(src%base_time, src%time_sec(irec))
      !----- The values as the source stores them; from here on, one ingest for both. -----------!
      if (src%backend == MET_BACKEND_ED_ERA5LAND) then
         call archive_source_record(src, cur, irec, rec, humidity_value, sw_total)
      else
         call file_source_record(src, irec, rec, humidity_value, sw_total)
      end if
      call assert_finite(rec%tair_k, 'Tair', irec, src%grid_index)
      call assert_finite(humidity_value, MEDS_FIELD(humidity_field(src%humidity)), irec, src%grid_index)
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
         call assert_finite(rec%par_beam,    'SWdown_par_beam',    irec, src%grid_index)   ! required source
         call assert_finite(rec%par_diffuse, 'SWdown_par_diffuse', irec, src%grid_index)   ! fields -> no gap-fill
         call assert_finite(rec%nir_beam,    'SWdown_nir_beam',    irec, src%grid_index)
         call assert_finite(rec%nir_diffuse, 'SWdown_nir_diffuse', irec, src%grid_index)
      else
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

end module meds_met_driver
