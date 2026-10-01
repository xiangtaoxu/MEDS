! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_archive_source -- the ED_ERA5land ARCHIVE source (MEDS_FORCING_DESIGN.md section    !
! 15): monthly files of hourly ERA5-Land fields over a grid of cells, laid end to end as one    !
! hourly axis over the months a run needs. It opens the archive (archive_source_open), loads   !
! a month before the step that reads it (archive_prefetch), and hands the shared ingest one    !
! record's values as stored (archive_source_record). File I/O is meds_era5land_reader's.      !
!==========================================================================================!
module meds_met_archive_source
   use iso_c_binding, only : c_int
   use meds_kinds, only : wp, ik
   use meds_time, only : meds_time_t, time_advance_seconds, seconds_between, seconds_into_day,    &
                         time_to_string
   use meds_forcing_config, only : MET_PATH_LEN, ARCHIVE_DT_SEC, SWPART_PASSTHROUGH
   use meds_forcing_types, only : met_record_t, met_source_t, met_cursor_t, met_cells_t,          &
                                  HUMIDITY_TDEW
   use meds_era5land_reader, only : era5land_path, era5land_default_template,                     &
                                    era5land_default_static, era5land_select_site,                &
                                    era5land_load_month, era5land_month_hours, ERA_NVAR,          &
                                    ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10, ERA_RAINF,   &
                                    ERA_SWDOWN, ERA_LWDOWN, ERA_VAR_NAME, ERA_EPOCH, ERA_OK,      &
                                    ERA_ERR_OPEN, ERA_ERR_NO_CELL
   use meds_netcdf_c, only : nc_open_f, nc_close, NC_NOERR, NC_NOWRITE
   use meds_met_source_common, only : MET_OK, MET_ERR_ATTR_MISMATCH, MET_ERR_ARCHIVE,             &
                                      REC_MATCH_TOL, validate_recycle_window, file_lookup_sec,    &
                                      first_record_after, validate_file_against_config
   implicit none
   private

   public :: archive_source_open, archive_prefetch, archive_source_record

contains

   !=======================================================================================!
   !  OPEN the ED_ERA5land archive (§15): its cells and month axis (open_archive), the        !
   !  declared recycle window checked against that axis, and the axis's first month loaded,  !
   !  with the window's first day kept for the recycle seam. A rejection returns its status   !
   !  and the reason; the caller closes the source.                                           !
   !=======================================================================================!
   subroutine archive_source_open(src, stat, why, run_start, run_end, cells)
      type(met_source_t),          intent(inout) :: src
      integer(ik),                 intent(out)   :: stat
      character(len=*),            intent(out)   :: why          !< what stops the run when stat /= MET_OK
      type(meds_time_t), optional, intent(in)    :: run_start, run_end
      type(met_cells_t), optional, intent(in)    :: cells        !< a region's cells; absent: the site's
      why = ''
      call open_archive(src, run_start, run_end, stat, cells)
      if (stat == MET_OK) call validate_recycle_window(src, stat)
      if (stat /= MET_OK) then
         why = 'met_open: the ED_ERA5land archive cannot drive this run (see the message above)' ; return
      end if
      call load_axis_month(src, 1_ik)
      call keep_window_head(src)
      src%rec_first = 1_ik
   end subroutine archive_source_open

   !=======================================================================================!
   !  PREFETCH (R1) from the archive: load what the step starting at `step_start` reads. A    !
   !  daily step from midnight reads one archive                                             !
   !  month plus the record before it: 00:00 on the 1st, which lives in the previous month's    !
   !  file, or the window's last record at the recycle wrap. Moving into the next month, that   !
   !  record comes from the outgoing buffer, so a run loads each month once. A step the recycle  !
   !  seam falls inside also reads the window's first day, which met_open keeps. validate_config !
   !  restricts format = "era5land" to daily steps from midnight, the shape this assumes.       !
   !=======================================================================================!
   subroutine archive_prefetch(src, step_start)
      type(met_source_t), intent(inout)  :: src
      type(meds_time_t),  intent(in)    :: step_start
      real(wp)    :: s0
      integer(ik) :: r, k, kp, p
      logical     :: wrap
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
   end subroutine archive_prefetch

   !----- One record's values at the cursor's cell, as the archive stores them (§15.4): the     !
   !      dewpoint, whose humidity comes from the model's own saturation curve in the shared    !
   !      ingest, the wind components, whose speed is the vector's at each stamp, and the total   !
   !      shortwave. The archive never carries the four shortwave streams (archive_source_open    !
   !      refuses "passthrough"). -----------------------------------------------------------------!
   subroutine archive_source_record(src, cur, irec, rec, humidity_value, sw_total)
      type(met_source_t), intent(in)    :: src
      type(met_cursor_t), intent(in)    :: cur
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(inout) :: rec
      real(wp),           intent(out)   :: humidity_value, sw_total
      integer(ik) :: h
      call locate_record(src, irec, h)
      rec%tair_k   = archive_value(src, cur, h, ERA_TAIR)
      rec%psurf_pa = archive_value(src, cur, h, ERA_PSURF)
      humidity_value = archive_value(src, cur, h, ERA_TDEW)
      rec%wind_u   = archive_value(src, cur, h, ERA_U10)
      rec%wind_v   = archive_value(src, cur, h, ERA_V10)
      rec%wind     = sqrt(rec%wind_u**2 + rec%wind_v**2)
      rec%rainf    = archive_value(src, cur, h, ERA_RAINF)            ! total rainfall rate [kg/m2/s]
      rec%lwdown   = archive_value(src, cur, h, ERA_LWDOWN)
      sw_total     = archive_value(src, cur, h, ERA_SWDOWN)
   end subroutine archive_source_record

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
            src%time_sec(src%month_rec0(k) + h) = first + ARCHIVE_DT_SEC * real(h - 1_ik, wp)
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
      hours = seconds_between(ERA_EPOCH, t) / ARCHIVE_DT_SEC
      if (round_up) then ; hours = real(ceiling(hours), wp) ; else ; hours = real(floor(hours), wp) ; end if
      before = time_advance_seconds(ERA_EPOCH, hours * ARCHIVE_DT_SEC - 1.0_wp)
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

end module meds_met_archive_source
