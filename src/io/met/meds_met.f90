! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met -- the meteorological-forcing reader (docs/science/forcing.md §1-§3, §10-§12).   !
! met_open opens a MEDS forcing file, the ED_ERA5land archive or the no-file CONST backend into the !
! met_source_t a run shares; met_cursor_init binds a polygon's met_cursor_t to it; met_prefetch     !
! loads what a step reads; met_advance and met_instant step and sample a cursor with no file        !
! access. read_record is the one ingest for both file sources (meds_met_default,                !
! meds_met_era5land): it checks, converts, lapses and partitions a record, and a missing      !
! value stops the run. CO2 comes from [forcing].co2_source on model time, never from the met file.  !
!==========================================================================================!
module meds_met
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use iso_c_binding, only : c_int
   use meds_kinds, only : wp, ik
   use meds_therm_lib, only : air_density
   use meds_time, only : meds_time_t, time_advance_seconds, seconds_into_day, time_to_string
   use meds_forcing_params, only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT,   &
                                   MET_BACKEND_ED_ERA5LAND, METAVG_END, METAVG_BEGIN, METAVG_CENTER, &
                                   SWPART_PASSTHROUGH, CLAMP_ERROR, INTERP_LINEAR,                &
                                   LW_SYNTHESIZE, CO2_SOURCE_FILE
   use meds_met_types, only : met_forcing_t, met_record_t, met_source_t, met_cursor_t,        &
                                  met_cells_t, HUMIDITY_RHAIR, HUMIDITY_TDEW, MEDS_FIELD
   use meds_lapse_rate, only : lapse_air_temperature, lapse_pressure, monthly_lapse_rate,         &
                               lapse_specific_humidity, lapse_longwave
   use meds_met_co2, only : co2_series_read, co2_series_at, co2_series_covers,                 &
                               co2_series_end, co2_series_free
   use meds_forcing_kernels, only : interpolate_forcing, interpolate_wind_energy,                 &
                                    met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,     &
                                    partition_shortwave, clearness_index,                         &
                                    synthesize_lwdown, dewpoint_to_specific_humidity,             &
                                    rh_to_specific_humidity
   use meds_netcdf_c, only : nc_close, nc_check
   use meds_met_share, only : MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS,                     &
                                      MET_ERR_START_NOT_A_RECORD, MET_ERR_WINDOW_NOT_COVERED,     &
                                      MET_ERR_DT_MISMATCH, MET_ERR_AXIS_NOT_UNIFORM,              &
                                      MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE, MET_ERR_CO2_FILE,   &
                                      MET_ERR_CO2_NOT_COVERED, MET_ERR_CO2_IN_MET_FILE,           &
                                      MET_ERR_NOT_UTC, MET_ERR_HUMIDITY, file_lookup_sec,         &
                                      humidity_field
   use meds_met_default, only : default_file_open, default_file_record
   use meds_met_era5land, only : era5land_open, era5land_prefetch,                     &
                                       era5land_record
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
   !  OPEN: the CO2, then the backend: CONST (the reference climate), a MEDS forcing file     !
   !  (default_file_open) or the archive (era5land_open). With `stat` present a          !
   !  rejection returns its code, so the tests can exercise it; absent, it stops the run.     !
   !  run_start and run_end bound what a non-recycling run reads and the CO2 must cover.      !
   !=======================================================================================!
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
      src%backend    = fcfg%format
      src%grid_index = fcfg%grid_index
      src%dt_forcing = fcfg%timestep
      src%has_wind_vector = .false.

      !----- CO2 first, so every backend takes it the same way. ---------------------------------!
      call open_co2(src, run_start, run_end, vstat)
      why = 'met_open: the CO2 file cannot drive this run (see the message above)'
      if (vstat == MET_OK) then
         select case (fcfg%format)
         case (MET_BACKEND_CONST)                                  ! the reference climate: no file
            src%ngrid = 1_ik
            src%nrec = 0_ik
         case (MET_BACKEND_ED_ERA5LAND)
            call era5land_open(src, vstat, why, run_start, run_end, cells)
         case default                                              ! MET_BACKEND_ED_DEFAULT
            call default_file_open(src, vstat, why, run_start, run_end)
         end select
      end if
      if (vstat /= MET_OK) then
         !----- Close the file first: an open netCDF handle keeps an HDF5 lock on the path, and a  !
         !      later open of the same file fails with a bare "permission denied". ---------------!
         call met_close(src)
         if (present(stat)) then
            stat = vstat
            return
         end if
         !----- The reason, then a fixed stop code: ifx garbles a stop code that is not a constant. -!
         write(*,'(1x,a)') trim(why)
         error stop 'met_open: the forcing cannot drive this run (see the messages above)'
      end if
   end subroutine met_open

   !----- CO2 (docs/science/forcing.md §12): with co2_source = "file", read the MEDS CO2 file   !
   !      and check it covers the run. MEDS does not extrapolate CO2, so a run the file does not  !
   !      cover stops here rather than decades in. The run span is absent only in unit tests. ---!
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
         stat = MET_ERR_CO2_FILE
         return
      end if
      if (present(run_start) .and. present(run_end)) then
         if (.not. co2_series_covers(src%co2, run_start, run_end)) then
            write(*,'(2a)') ' met_open: the CO2 file does not cover the run: ', trim(src%fcfg%co2_file)
            write(*,'(4a)') '   run  ', time_to_string(run_start), ' .. ', time_to_string(run_end)
            write(*,'(4a)') '   file ', time_to_string(src%co2%first_start), ' .. ',               &
                            time_to_string(co2_series_end(src%co2))
            write(*,'(a)')  '   MEDS does not extrapolate CO2: extend the file, or shorten the run.'
            stat = MET_ERR_CO2_NOT_COVERED
            return
         end if
      end if
   end subroutine open_co2

   !=======================================================================================!
   !  CURSOR: bind one polygon to the open source -- its cell and location -- and load its first  !
   !  bracket. The grid elevation (the lapse origin) is the cell's orography for the archive and   !
   !  [site].grid_elevation for a MEDS forcing file.                                               !
   !=======================================================================================!
   subroutine met_cursor_init(src, cur, cell, latitude_deg, longitude_deg, elevation_m)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(out) :: cur
      integer(ik),        intent(in)  :: cell
      real(wp),           intent(in)  :: latitude_deg, longitude_deg, elevation_m
      cur%cell = cell
      cur%latitude_deg = latitude_deg
      cur%longitude_deg = longitude_deg
      cur%elevation_m  = elevation_m
      cur%grid_elevation_m = src%fcfg%grid_elevation
      if (src%backend == MET_BACKEND_ED_ERA5LAND) cur%grid_elevation_m = src%cells%elevation(cell)
      if (src%backend == MET_BACKEND_CONST) then
         cur%rec_prev = met_record_t()
         cur%rec_next = met_record_t()
         return
      end if
      call load_bracket(src, cur, src%rec_first)
   end subroutine met_cursor_init

   !=======================================================================================!
   !  PREFETCH: load what the step starting at `step_start` reads, before its compute phase.   !
   !  A MEDS forcing file was read at open; the archive loads a month at a time.               !
   !=======================================================================================!
   subroutine met_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(meds_time_t),  intent(in)    :: step_start
      if (src%backend == MET_BACKEND_ED_ERA5LAND) call era5land_prefetch(src, step_start)
   end subroutine met_prefetch

   !=======================================================================================!
   !  ADVANCE: slide this cursor's bracket so rec_prev%when <= now < rec_next%when: a start     !
   !  before the first record holds or stops (start_clamp), the recycle seam wraps, and a run   !
   !  past the last record clamps to the last interval (docs/science/forcing.md §2).            !
   !=======================================================================================!
   subroutine met_advance(src, cur, now)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      type(meds_time_t),  intent(in)    :: now
      real(wp)    :: now_sec
      integer(ik) :: i

      if (src%backend == MET_BACKEND_CONST) return

      !----- The same recycle-mapped seconds on the file axis that met_instant uses. ----------!
      now_sec = file_lookup_sec(src, now)

      !----- A start before the first record. Hold loads records 1-2 and met_instant clamps     !
      !      w_next to 0; rec_next is left alone, or a stale copy would stop the reload at t1. ---!
      if (now_sec < src%time_sec(1)) then
         if (src%fcfg%start_clamp == CLAMP_ERROR) then
            error stop 'met_advance: model start precedes the first forcing record (start_clamp=error)'
         end if
         if (cur%irec_prev /= 1_ik) call load_bracket(src, cur, 1_ik)
         return
      end if

      !----- The recycle seam: the window's last interval wraps to its first record              !
      !      (docs/science/forcing.md §9). The seam is the window's edge, not the file's. --------!
      !      Nested, not one .and.: Fortran need not short-circuit, and irec_cycle_last is 0 unless recycling.
      if (src%n_cycle_years >= 1_ik .and. src%irec_cycle_last > src%irec_cycle_first) then
         if (now_sec >= src%time_sec(src%irec_cycle_last)) then
            if (.not. cur%at_wrap_seam) call load_wrap_bracket(src, cur)
            return
         end if
      end if

      !----- A non-recycling run past the end: clamp to the last interval. ---------------------!
      if (.not. src%fcfg%recycle .and. src%nrec >= 2_ik .and. now_sec >= src%time_sec(src%nrec)) then
         if (cur%irec_prev /= src%nrec - 1_ik) call load_bracket(src, cur, src%nrec - 1_ik)
         return
      end if

      !----- Search from the current bracket (cheap for a forward march). ---------------------!
      i = max(1_ik, cur%irec_prev)
      do while (i < src%nrec - 1_ik)
         if (src%time_sec(i + 1_ik) <= now_sec) then
            i = i + 1_ik
         else
            exit
         end if
      end do
      do while (i > 1_ik)
         if (src%time_sec(i) > now_sec) then
            i = i - 1_ik
         else
            exit
         end if
      end do
      if (i /= cur%irec_prev) call load_bracket(src, cur, i)
      !----- The daytime clearness the longwave synthesis holds after dark: updated here, outside !
      !      the patch loop, because met_instant is intent(in) (docs/science/forcing.md §11). ----!
      if (src%fcfg%lwdown_source == LW_SYNTHESIZE) call remember_clearness(src, cur, now)
   end subroutine met_advance

   !----- Update cur%kt_last_day when the sun is up, from the interval's mean shortwave against   !
   !      the sun at `now` (an approximation, docs/science/forcing.md §11). -----------------------!
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
   !      later of the bracket on an avg_convention = "end" file, the earlier on a "begin" one     !
   !      (docs/science/forcing.md §3). Every interval mean is taken from here. -----------------!
   pure function interval_mean_record(src, cur) result(r)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(met_record_t) :: r
      select case (src%fcfg%avg_convention)
      case (METAVG_BEGIN)
         r = cur%rec_prev
      case default
         r = cur%rec_next            ! METAVG_END (ERA5-Land) and any other
      end select
   end function interval_mean_record

   !=======================================================================================!
   !  INSTANT: the bracket at the model instant `now` (docs/science/forcing.md §3-§5): states     !
   !  linear, wind in energy form, precipitation held, shortwave reconstructed from the interval   !
   !  means, cosz and rho_air derived. A pure function, so it is safe to sample ahead.             !
   !=======================================================================================!
   pure function met_instant(src, cur, now) result(met)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      type(meds_time_t),  intent(in) :: now
      type(met_forcing_t) :: met
      type(met_record_t)  :: mean_rec
      type(meds_time_t)   :: mws
      real(wp) :: now_sec, tprev, tnext, w_next, cosz_now, factor, win_start_sec
      associate (f => src%fcfg, p => cur%rec_prev, n => cur%rec_next)

      !----- The solar zenith at `now` (met_solar_cosz turns UTC into apparent solar time). ---!
      cosz_now = met_solar_cosz(now, seconds_into_day(now), cur%latitude_deg, cur%longitude_deg)
      met%cosz = cosz_now

      !----- CO2 on model time, one way for every backend (docs/science/forcing.md eq. (18)). ----!
      if (f%co2_source == CO2_SOURCE_FILE) then
         met%co2 = co2_series_at(src%co2, now)
      else
         met%co2 = f%co2_const
      end if

      !----- CONST keeps the type's defaults; a file backend interpolates its bracket. ----------!
      if (src%backend /= MET_BACKEND_CONST) then
         !----- The weight within the bracket, on the same recycle-mapped seconds as met_advance. -!
         now_sec = file_lookup_sec(src, now)
         if (cur%at_wrap_seam) then                    ! cycle-boundary: window's last rec -> its first
            tprev = src%time_sec(src%irec_cycle_last)
            tnext = tprev + src%dt_forcing
         else
            tprev = src%time_sec(cur%irec_prev)
            tnext = src%time_sec(min(cur%irec_prev + 1_ik, src%nrec))
         end if
         if (tnext > tprev) then
            w_next = min(1.0_wp, max(0.0_wp, (now_sec - tprev) / (tnext - tprev)))
         else
            w_next = 0.0_wp                                          ! degenerate (start-clamp hold)
         end if

         !----- States linear; wind in energy form (docs/science/forcing.md eq. (1)). ----------!
         met%tair_k   = interpolate_forcing(INTERP_LINEAR, p%tair_k,   n%tair_k,   w_next)
         met%qair     = interpolate_forcing(INTERP_LINEAR, p%qair,     n%qair,     w_next)
         met%psurf_pa = interpolate_forcing(INTERP_LINEAR, p%psurf_pa, n%psurf_pa, w_next)
         met%lwdown   = interpolate_forcing(INTERP_LINEAR, p%lwdown,   n%lwdown,   w_next)
         met%wind     = interpolate_wind_energy(p%wind, n%wind, w_next)   ! floored at aerodynamics.ubmin where used
         !----- The vector interpolates linearly, which keeps its direction; never floored. --------!
         if (src%has_wind_vector) then
            met%wind_u = interpolate_forcing(INTERP_LINEAR, p%wind_u, n%wind_u, w_next)
            met%wind_v = interpolate_forcing(INTERP_LINEAR, p%wind_v, n%wind_v, w_next)
            met%has_wind_vector = .true.
         end if

         !----- Rain and shortwave both come from the record whose interval contains `now`. -----!
         mean_rec = interval_mean_record(src, cur)

         !----- Precipitation: the interval's total rate, held across it; each patch splits it. -!
         met%precip = mean_rec%rainf

         !----- Shortwave: each mean stream times cosz(now)/<cosz>_win (docs/science/forcing.md      !
         !      eq. (5)), the window anchored on the model calendar (mws), so both see the same sun. !
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

      !----- The longwave synthesis (docs/science/forcing.md §11), after the shortwave: this      !
      !      instant's clearness index by day, the remembered one after dark. ----------------------!
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
         st = nc_close(int(src%ncid, c_int))
         call nc_check(st, 'met_close: nc_close')
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
      src%buffer%year = 0_ik
      src%carry_rec = 0_ik
      src%n_head = 0_ik
      src%n_loads = 0_ik
   end subroutine met_close

   !----- Load rec_prev = record(irec), rec_next = record(irec+1) (clamped at the end). A slide   !
   !      by one record reuses the old upper record as the new lower one. ------------------------!
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

   !----- The seam bracket: rec_prev = the window's last record, rec_next = its first, so the   !
   !      last interval interpolates across the wrap instead of clamping. -----------------------!
   subroutine load_wrap_bracket(src, cur)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(inout) :: cur
      call read_record(src, cur, src%irec_cycle_last,  cur%rec_prev)
      call read_record(src, cur, src%irec_cycle_first, cur%rec_next)
      cur%irec_prev    = src%irec_cycle_last
      cur%at_wrap_seam = .true.
   end subroutine load_wrap_bracket

   !----- The one ingest: record irec at this polygon's cell, checked, converted, lapsed and    !
   !      partitioned (docs/science/forcing.md §6-§8). ------------------------------------------!
   subroutine read_record(src, cur, irec, rec)
      type(met_source_t), intent(in)  :: src
      type(met_cursor_t), intent(in) :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(out)   :: rec
      real(wp)    :: sw_total, humidity_value
      rec%when   = time_advance_seconds(src%base_time, src%time_sec(irec))
      !----- The values as the source stores them; from here on, one ingest for both. -----------!
      if (src%backend == MET_BACKEND_ED_ERA5LAND) then
         call era5land_record(src, cur, irec, rec, humidity_value, sw_total)
      else
         call default_file_record(src, irec, rec, humidity_value, sw_total)
      end if
      call assert_finite(rec%tair_k, 'Tair', irec, src%grid_index)
      call assert_finite(humidity_value, MEDS_FIELD(humidity_field(src%humidity)), irec, src%grid_index)
      call assert_finite(rec%psurf_pa, 'PSurf', irec, src%grid_index)
      !----- Specific humidity at the forcing's own temperature and pressure: a dewpoint and a     !
      !      relative humidity both go through the model's liquid-water saturation curve           !
      !      (docs/science/forcing.md eq. (11)). ---------------------------------------------------!
      select case (src%humidity)
      case (HUMIDITY_RHAIR)
         rec%qair = rh_to_specific_humidity(humidity_value, rec%tair_k, rec%psurf_pa)
      case (HUMIDITY_TDEW)
         rec%qair = dewpoint_to_specific_humidity(humidity_value, rec%psurf_pa)
      case default
         rec%qair = humidity_value
      end select
      call assert_finite(rec%wind, 'Wind', irec, src%grid_index)
      call assert_finite(rec%rainf, 'Rainf', irec, src%grid_index)
      if (src%fcfg%lwdown_source /= LW_SYNTHESIZE)                                                  &
         call assert_finite(rec%lwdown, 'LWdown', irec, src%grid_index)

      !----- The terrain lapse, from the cell's elevation to the polygon's                        !
      !      (docs/science/forcing.md §8 Step 1, eqs. (13)-(14)). ----------------------------------!
      if (src%fcfg%apply_elevation_lapse) then
         block
            real(wp) :: dz, gamma, t_grid, p_grid, q_grid
            dz     = cur%elevation_m - cur%grid_elevation_m       ! + when the polygon is higher
            gamma  = monthly_lapse_rate(src%fcfg%lapse_rate_tair, rec%when%month)
            t_grid = rec%tair_k
            p_grid = rec%psurf_pa
            q_grid = rec%qair   ! before it is overwritten
            rec%psurf_pa = lapse_pressure(p_grid, t_grid, dz, gamma)
            rec%tair_k   = lapse_air_temperature(t_grid, dz, gamma)
            !----- Relative humidity is held: a measured one is used as is, any other is formed at the cell. !
            if (src%humidity == HUMIDITY_RHAIR) then
               rec%qair  = rh_to_specific_humidity(humidity_value, rec%tair_k, rec%psurf_pa)
            else
               rec%qair  = lapse_specific_humidity(q_grid, t_grid, p_grid, rec%tair_k, rec%psurf_pa)
            end if
            !----- Only a file's longwave: a synthesized one is built later from the lapsed T and q. -!
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
         call partition_shortwave(sw_total, record_mean_cosz(src, cur, rec%when), rec%psurf_pa,    &
                                  src%fcfg%sw_partition,                                          &
                                  rec%par_beam, rec%par_diffuse, rec%nir_beam, rec%nir_diffuse)
      end if
   end subroutine read_record

   !----- The sun the shortwave partition uses for one record: the mean of max(cosz, 0) over the  !
   !      interval it averages, or the sun at the stamp for an instantaneous record             !
   !      (docs/science/forcing.md §6). --------------------------------------------------------!
   real(wp) function record_mean_cosz(src, cur, when) result(cosz)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(meds_time_t),  intent(in) :: when
      real(wp) :: stamp_sec, start_sec, factor
      stamp_sec = seconds_into_day(when)
      select case (src%fcfg%avg_convention)
      case (METAVG_END)
         start_sec = stamp_sec - src%dt_forcing
      case (METAVG_BEGIN)
         start_sec = stamp_sec
      case (METAVG_CENTER)
         start_sec = stamp_sec - 0.5_wp * src%dt_forcing
      case default
         cosz = met_solar_cosz(when, stamp_sec, cur%latitude_deg, cur%longitude_deg)
         return
      end select
      factor = cosz_reconstruct_factor(when, start_sec, src%dt_forcing / real(N_COSZ_SUB, wp),       &
                                       src%dt_forcing, cur%latitude_deg, cur%longitude_deg)
      cosz = 0.0_wp
      if (factor > 0.0_wp) cosz = 1.0_wp / factor              ! factor = 1 / <cosz>_interval
   end function record_mean_cosz

   !----- MEDS never gap-fills: a NaN in a required field halts the run (docs/science/forcing.md !
   !      §10). ------------------------------------------------------------------------------!
   subroutine assert_finite(x, name, irec, grid)
      real(wp),         intent(in) :: x
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: irec, grid
      if (ieee_is_nan(x)) then
         write(*,'(4a,i0,a,i0)') 'met_driver: missing/NaN value in required field ', trim(name), &
               ' at time record ', ' ', irec, ', grid ', grid
         error stop 'met_driver: forcing has a data gap (MEDS does not gap-fill; fix the source upstream)'
      end if
   end subroutine assert_finite

end module meds_met
