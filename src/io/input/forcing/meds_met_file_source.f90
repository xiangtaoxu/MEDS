! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_met_file_source -- the MEDS forcing FILE source (format = "ED_default",                 !
! MEDS_FORCING_DESIGN.md section 4): a multi-grid netCDF file whose records the run needs are  !
! read into memory at open, so no step reads the file. It opens the file (file_source_open)   !
! and hands the shared ingest one record's values as stored (file_source_record).             !
!==========================================================================================!
module meds_met_file_source
   use iso_c_binding, only : c_int, c_size_t
   use meds_kinds, only : wp, ik
   use meds_time, only : meds_time_t, time_units_base, seconds_between
   use meds_forcing_config, only : SWPART_PASSTHROUGH, GRIDMATCH_NEAREST, LW_SYNTHESIZE
   use meds_forcing_types, only : met_record_t, met_source_t, HUMIDITY_QAIR, HUMIDITY_RHAIR,      &
                                  HUMIDITY_TDEW, N_MEDS_FIELD, MEDS_FIELD, FLD_TAIR, FLD_RHAIR,   &
                                  FLD_PSURF, FLD_WIND, FLD_U10, FLD_V10, FLD_RAINF, FLD_LWDOWN,   &
                                  FLD_SWDOWN, FLD_PAR_BEAM, FLD_PAR_DIFFUSE, FLD_NIR_BEAM,        &
                                  FLD_NIR_DIFFUSE
   use meds_forcing_kernels, only : nearest_grid_index
   use meds_netcdf_c, only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f, nc_get_att_text_f,       &
                             nc_get_vara_double, nc_close, nc_check, nc_inq_var_chunking,         &
                             NC_NOERR, NC_NOWRITE, NC_CHUNKED
   use meds_met_source_common, only : MET_OK, MET_ERR_HUMIDITY, validate_recycle_window,          &
                                      first_record_after, validate_file_against_config,           &
                                      humidity_field
   implicit none
   private

   public :: file_source_open, file_source_record

   !----- RHair is a fraction; a value above this is a percentage written into a fractional field. !
   real(wp), parameter :: RH_FRACTION_MAX = 1.5_wp

contains

   !=======================================================================================!
   !  OPEN a MEDS forcing file (ED_default): the grid and time dims, the time axis and its     !
   !  base time from the `time:units` attribute, the checks of the file against the config and !
   !  the declared recycle window, and the records the run needs, read into memory. A         !
   !  rejection returns its status and the reason; the caller closes the source.              !
   !=======================================================================================!
   subroutine file_source_open(src, stat, why, run_start, run_end)
      type(met_source_t),          intent(inout) :: src
      integer(ik),                 intent(out)   :: stat
      character(len=*),            intent(out)   :: why          !< what stops the run when stat /= MET_OK
      type(meds_time_t), optional, intent(in)    :: run_start, run_end
      integer(c_int)    :: st, ncid
      integer(c_size_t) :: dlen
      character(len=256):: units
      logical           :: ok

      stat = MET_OK ; why = ''
      st = nc_open_f(trim(src%fcfg%path), NC_NOWRITE, ncid)
      call nc_check(st, 'met_open: nc_open '//trim(src%fcfg%path))
      src%ncid = int(ncid, ik)

      st = nc_inq_dimlen_f(ncid, 'grid', dlen) ; call nc_check(st, 'met_open: grid dim')
      src%ngrid = int(dlen, ik)
      !----- bind this polygon to a grid slice: explicit index, or nearest [site] lat/lon (§4.1). !
      if (src%fcfg%grid_match == GRIDMATCH_NEAREST) then
         call resolve_grid_index(src)                        ! sets src%grid_index in 1..ngrid
      else                                                   ! GRIDMATCH_EXPLICIT
         if (src%fcfg%grid_index < 1_ik .or. src%fcfg%grid_index > src%ngrid) then
            write(*,'(a,i0,a,i0)') 'met_open: grid_index ', src%fcfg%grid_index, ' out of range 1..', src%ngrid
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
         call time_units_base(units, src%base_time, ok)
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
      call detect_humidity(src, ncid, stat)
      if (stat /= MET_OK) then
         why = 'met_open: the forcing file''s humidity is ambiguous or absent (see the message above)' ; return
      end if

      !----- V4 (#185): the file's own record spacing against [forcing].dt_forcing, and the two   !
      !      global attributes the prep script writes against the config that claims to describe    !
      !      the same file. A file that disagreed with its config would otherwise be silently       !
      !      mis-timed or mis-partitioned. ------------------------------------------------------------!
      call validate_file_against_config(src, ncid, stat)
      if (stat /= MET_OK) then
         why = 'met_open: the forcing file contradicts [forcing] (see the message above)' ; return
      end if

      !----- V2/V3: check the DECLARED recycle window against this file's actual record stamps.  !
      !      Cannot live in the config sanity check (V1, whole-year span) -- it needs the time      !
      !      axis, which is only known here. -------------------------------------------------------!
      call validate_recycle_window(src, stat)
      if (stat /= MET_OK) then
         why = 'met_open: forcing file does not match the declared [forcing] recycle window' ; return
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
            j = src%series_col(FLD_RHAIR)
            if (maxval(src%series(:, j)) > RH_FRACTION_MAX) then
               write(*,'(a,f0.2,a)') ' met_open: RHair reaches ', maxval(src%series(:, j)),                &
                                     ', but it is a fraction (units "1"), not a percentage.'
               stat = MET_ERR_HUMIDITY ; why = 'met_open: RHair is in percent (see the message above)'
               return
            end if
         end block
      end if
   end subroutine file_source_open

   !----- One record's values at the cell, as the file stores them: the humidity in its own form,  !
   !      and the shortwave total, or its four streams under sw_partition = "passthrough". The    !
   !      shared ingest (meds_met_driver's read_record) converts, checks and lapses them. -------!
   subroutine file_source_record(src, irec, rec, humidity_value, sw_total)
      type(met_source_t), intent(in)    :: src
      integer(ik),        intent(in)    :: irec
      type(met_record_t), intent(inout) :: rec
      real(wp),           intent(out)   :: humidity_value, sw_total
      rec%tair_k   = series_value(src, FLD_TAIR,  irec)
      humidity_value = series_value(src, humidity_field(src%humidity), irec)
      rec%psurf_pa = series_value(src, FLD_PSURF, irec)
      if (src%has_wind_vector) then
         rec%wind_u = series_value(src, FLD_U10, irec)
         rec%wind_v = series_value(src, FLD_V10, irec)
         rec%wind   = sqrt(rec%wind_u**2 + rec%wind_v**2)
      else
         rec%wind   = series_value(src, FLD_WIND,  irec)
      end if
      rec%rainf    = series_value(src, FLD_RAINF, irec)              ! total rainfall rate [kg/m2/s]
      !----- LWdown is OPTIONAL when we are synthesizing it (#182): a source without longwave is    !
      !      exactly the case lwdown_source = "synthesize" exists for, so demanding the variable      !
      !      would defeat the feature. Read it when present either way -- it costs nothing and keeps  !
      !      the record complete for diagnostics.  ---------------------------------------------------!
      if (src%fcfg%lwdown_source == LW_SYNTHESIZE) then
         rec%lwdown = series_value(src, FLD_LWDOWN, irec, default=0.0_wp)
      else
         rec%lwdown = series_value(src, FLD_LWDOWN, irec)
      end if
      sw_total = 0.0_wp
      if (src%fcfg%sw_partition == SWPART_PASSTHROUGH) then
         rec%par_beam    = series_value(src, FLD_PAR_BEAM,    irec)
         rec%par_diffuse = series_value(src, FLD_PAR_DIFFUSE, irec)
         rec%nir_beam    = series_value(src, FLD_NIR_BEAM,    irec)
         rec%nir_diffuse = series_value(src, FLD_NIR_DIFFUSE, irec)
      else
         sw_total = series_value(src, FLD_SWDOWN, irec)
      end if
   end subroutine file_source_record

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
         if (nc_inq_varid_f(ncid, trim(MEDS_FIELD(humidity_field(FORMS(i)))), vid) == NC_NOERR) then
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

   !----- A MEDS forcing file's fields at the cell, records r0..r1, read at open. Each field's    !
   !      column in the series is found here, once (series_col, 0 for a field the file lacks). -----!
   subroutine read_series(src, ncid, r0, r1)
      type(met_source_t), intent(inout)  :: src
      integer(c_int),     intent(in)    :: ncid
      integer(ik),        intent(in)    :: r0, r1
      integer(c_int)    :: st, vid(N_MEDS_FIELD)
      integer(c_size_t) :: start2(2), count2(2)
      logical     :: present_(N_MEDS_FIELD)
      integer(ik) :: j, n
      do j = 1_ik, N_MEDS_FIELD
         present_(j) = nc_inq_varid_f(ncid, trim(MEDS_FIELD(j)), vid(j)) == NC_NOERR
      end do
      if (allocated(src%series)) deallocate(src%series)
      allocate(src%series(r0:r1, count(present_)))                   ! indexed by record number
      start2 = [int(r0 - 1_ik, c_size_t), int(src%grid_index - 1_ik, c_size_t)]    ! [time, grid]
      count2 = [int(r1 - r0 + 1_ik, c_size_t), 1_c_size_t]
      src%series_col = 0_ik ; n = 0_ik
      do j = 1_ik, N_MEDS_FIELD
         if (.not. present_(j)) cycle
         n = n + 1_ik ; src%series_col(j) = n
         if (n == 1_ik) call note_one_record_chunks(ncid, vid(j), r1 - r0 + 1_ik, src%fcfg%path)
         st = nc_get_vara_double(ncid, vid(j), start2, count2, src%series(:, n))
         call nc_check(st, 'read_series: get '//trim(MEDS_FIELD(j)))
      end do
   end subroutine read_series

   !----- A file stored one time record per chunk is read chunk by chunk at open: about 3 s for  !
   !      BCI's 90,528 records, most of a short run's start-up. Say so once, with the fix. The     !
   !      values are the same either way. ------------------------------------------------------!
   subroutine note_one_record_chunks(ncid, varid, nrec, path)
      integer(c_int),   intent(in) :: ncid, varid
      integer(ik),      intent(in) :: nrec
      character(len=*), intent(in) :: path
      integer(c_int)    :: storage
      integer(c_size_t) :: chunks(8)
      if (nrec < 1000_ik) return
      if (nc_inq_var_chunking(ncid, varid, storage, chunks) /= NC_NOERR) return
      if (storage /= NC_CHUNKED .or. chunks(1) /= 1_c_size_t) return      ! chunks(1) is the time axis
      write(*,'(3a)') ' note: ', trim(path), ' is stored one time record per chunk, so reading it'
      write(*,'(a)')  '       at start-up is slow (about 3 s per 90,000 records). Rewriting it once with'
      write(*,'(a)')  '       "nccopy -c time/8760,grid/1 <file> <new file>" keeps every value and reads fast.'
   end subroutine note_one_record_chunks

   !----- Field fld at record irec. A field the file does not carry stops the run, unless the     !
   !      caller gives a default (LWdown when it is synthesized). -------------------------------!
   function series_value(src, fld, irec, default) result(val)
      type(met_source_t), intent(in) :: src
      integer(ik),        intent(in) :: fld, irec
      real(wp), optional, intent(in) :: default
      real(wp)    :: val
      integer(ik) :: j
      j = src%series_col(fld)
      if (j == 0_ik) then
         if (present(default)) then ; val = default ; return ; end if
         write(*,'(2a)') ' met_driver: the forcing file has no variable ', trim(MEDS_FIELD(fld))
         error stop 'met_driver: a required forcing variable is missing'
      end if
      call check_in_series(src, irec)
      val = src%series(irec, j)
   end function series_value

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

end module meds_met_file_source
