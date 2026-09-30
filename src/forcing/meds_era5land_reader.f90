! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_era5land_reader -- file access for the global per-variable monthly ED_ERA5land archive !
! (MEDS_FORCING_DESIGN.md sections 14-15): path templates, the static file, site and box cell  !
! selection (§15.5), and a month of every domain cell read one chunk column at a time into a   !
! float32 buffer (§15.3). The unit conversions live in the reader (meds_met_driver), which      !
! turns every nonzero status here into a hard error; returning a status instead of stopping is  !
! what lets the tests exercise the rejections. MEDS never gap-fills: a NaN in a selected cell   !
! is an error, never a fill.                                                                     !
!==========================================================================================!
module meds_era5land_reader
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use iso_c_binding,        only : c_int, c_size_t, c_double, c_float
   use meds_kinds,           only : wp, sp, ik
   use meds_time,            only : meds_time_t, days_in_month, seconds_between, time_from_string
   use meds_forcing_config,  only : MET_PATH_LEN
   use meds_forcing_types,   only : met_cells_t, met_month_t
   use meds_forcing_kernels, only : great_circle_distance
   use meds_netcdf_c,        only : nc_open_f, nc_inq_varid_f, nc_inq_dimlen_f, nc_get_att_text_f, &
                                    nc_get_att_double_f, &
                                    nc_get_vara_double, nc_get_vara_float, nc_close, NC_NOERR, NC_NOWRITE
   implicit none
   private

   public :: era5land_path, era5land_default_template, era5land_default_static
   public :: era5land_select_site, era5land_select_box, era5land_load_month, era5land_month_hours
   public :: ERA_NVAR, ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10, ERA_RAINF, ERA_SWDOWN, ERA_LWDOWN
   public :: ERA_VAR_NAME, ERA_VAR_UNITS, ERA_EPOCH
   public :: ERA_OK, ERA_ERR_OPEN, ERA_ERR_GRID, ERA_ERR_TIME, ERA_ERR_UNITS, ERA_ERR_NAN, ERA_ERR_NO_CELL
   public :: ERA_ERR_FILL

   !----- The archive's eight variables, their file names and the units the builder writes. -----!
   integer(ik), parameter :: ERA_NVAR = 8_ik
   integer(ik), parameter :: ERA_TAIR = 1_ik, ERA_TDEW = 2_ik, ERA_PSURF = 3_ik, ERA_U10 = 4_ik,     &
                             ERA_V10 = 5_ik, ERA_RAINF = 6_ik, ERA_SWDOWN = 7_ik, ERA_LWDOWN = 8_ik
   character(len=6),  parameter :: ERA_VAR_NAME(ERA_NVAR) =                                        &
      ['Tair  ', 'Tdew  ', 'PSurf ', 'u10   ', 'v10   ', 'Rainf ', 'SWdown', 'LWdown']
   character(len=10), parameter :: ERA_VAR_UNITS(ERA_NVAR) =                                       &
      ['K         ', 'K         ', 'Pa        ', 'm s-1     ', 'm s-1     ', 'kg m-2 s-1',           &
       'W m-2     ', 'W m-2     ']
   !----- The archive's time coordinate is seconds since this instant (§14.2). -------------------!
   type(meds_time_t), parameter :: ERA_EPOCH = meds_time_t(1970_ik, 1_ik, 1_ik, 0_ik, 0_ik, 0_ik)

   integer(ik), parameter :: CHUNK = 16_ik       !< the archive's spatial chunk edge (§14.2)
   real(wp),    parameter :: GRID_TOL = 1.0e-4_wp !< [deg] slack for "these coordinates are equal"
   real(wp),    parameter :: STAMP_TOL = 0.5_wp  !< [s] slack for "this record is that stamp"

   !----- Status codes. --------------------------------------------------------------------------!
   integer(ik), parameter :: ERA_OK          = 0_ik
   integer(ik), parameter :: ERA_ERR_OPEN    = 1_ik   !< a file is missing or unreadable
   integer(ik), parameter :: ERA_ERR_GRID    = 2_ik   !< a file's grid differs from the static file's
   integer(ik), parameter :: ERA_ERR_TIME    = 3_ik   !< a month file's stamps are not that month's hours
   integer(ik), parameter :: ERA_ERR_UNITS   = 4_ik   !< a variable's units attribute is not the expected one
   integer(ik), parameter :: ERA_ERR_NAN     = 5_ik   !< a selected cell has a missing value
   integer(ik), parameter :: ERA_ERR_NO_CELL = 6_ik   !< no valid cell within the distance limit / in the box
   integer(ik), parameter :: ERA_ERR_FILL    = 7_ik   !< a variable's _FillValue is a number, not NaN

contains

   !=======================================================================================!
   !  Paths. Tokens {data_path}, {var}, {yyyy} and {mm}; the defaults are the archive's own    !
   !  flat layout (§14.1).                                                                     !
   !=======================================================================================!
   function era5land_path(template, data_path, var, year, month) result(path)
      character(len=*), intent(in) :: template, data_path, var
      integer(ik),      intent(in) :: year, month
      character(len=MET_PATH_LEN)  :: path
      character(len=4) :: yyyy
      character(len=2) :: mm
      character(len=:), allocatable :: s
      write(yyyy, '(i4.4)') year
      write(mm, '(i2.2)') month
      s = replace_all(trim(template), '{data_path}', trim(data_path))
      s = replace_all(s, '{var}', trim(var))
      s = replace_all(s, '{yyyy}', yyyy)
      s = replace_all(s, '{mm}', mm)
      path = s
   end function era5land_path

   pure function era5land_default_template() result(t)
      character(len=MET_PATH_LEN) :: t
      t = '{data_path}/ED_ERA5land_{var}_{yyyy}{mm}.nc'
   end function era5land_default_template

   pure function era5land_default_static() result(t)
      character(len=MET_PATH_LEN) :: t
      t = '{data_path}/ED_ERA5land_static.nc'
   end function era5land_default_static

   !----- Every occurrence of `token` in `s` replaced by `value`. ----------------------------------!
   pure function replace_all(s, token, value) result(out)
      character(len=*), intent(in)  :: s, token, value
      character(len=:), allocatable :: out
      integer :: k, pos
      out = ''
      pos = 1
      do
         k = index(s(pos:), token)
         if (k == 0) exit
         out = out // s(pos:pos + k - 2) // value
         pos = pos + k - 1 + len(token)
      end do
      out = out // s(pos:)
   end function replace_all

   !----- Hours in an archive month file: 01:00 on the 1st .. 00:00 on the 1st of the next month. --!
   elemental integer(ik) function era5land_month_hours(year, month) result(nt)
      integer(ik), intent(in) :: year, month
      nt = 24_ik * days_in_month(year, month)
   end function era5land_month_hours

   !=======================================================================================!
   !  SITE (§15.5): the cell containing the site by regular-grid arithmetic; if it has no data, !
   !  the nearest valid cell by great-circle distance, lowest (row-major) index on a tie, up to  !
   !  max_km. Rows farther than max_km in latitude alone cannot qualify, so only the rows within !
   !  that band are read, across all longitudes (which handles the antimeridian for free).       !
   !=======================================================================================!
   subroutine era5land_select_site(static_path, site_lat, site_lon, max_km, dom, distance_km, stat)
      character(len=*),   intent(in)  :: static_path
      real(wp),           intent(in)  :: site_lat, site_lon, max_km
      type(met_cells_t), intent(out) :: dom
      real(wp),           intent(out) :: distance_km
      integer(ik),        intent(out) :: stat
      integer(c_int) :: ncid, st
      real(wp), allocatable :: lat(:), lon(:), valid(:,:), elev(:,:)
      real(wp)    :: row_km, d, dlat
      integer(ik) :: j0, rj, r0, r1, j, i, best_j, best_i

      distance_km = huge(1.0_wp)
      call open_static(static_path, ncid, lat, lon, stat)
      if (stat /= ERA_OK) return
      dlat = lat(2) - lat(1)
      j0 = min(size(lat, kind=ik) - 1_ik, max(0_ik, nint((site_lat - lat(1)) / dlat, ik)))
      row_km = great_circle_distance(0.0_wp, lat(1), 0.0_wp, lat(2)) / 1000.0_wp
      rj = ceiling(max_km / row_km, ik) + 1_ik
      r0 = max(0_ik, j0 - rj) ; r1 = min(size(lat, kind=ik) - 1_ik, j0 + rj)
      call read_2d(ncid, 'valid', r0, r1 - r0 + 1_ik, 0_ik, size(lon, kind=ik), valid, stat)
      if (stat /= ERA_OK) then ; st = nc_close(ncid) ; return ; end if

      best_j = -1_ik ; best_i = -1_ik
      do j = r0, r1
         do i = 0_ik, size(lon, kind=ik) - 1_ik
            if (valid(i + 1_ik, j - r0 + 1_ik) < 0.5_wp) cycle
            d = great_circle_distance(site_lon, site_lat, lon(i + 1_ik), lat(j + 1_ik)) / 1000.0_wp
            if (d < distance_km) then
               distance_km = d ; best_j = j ; best_i = i
            end if
         end do
      end do
      if (best_j < 0_ik .or. distance_km > max_km) then
         st = nc_close(ncid) ; stat = ERA_ERR_NO_CELL ; return
      end if
      call read_2d(ncid, 'elevation', best_j, 1_ik, best_i, 1_ik, elev, stat)
      st = nc_close(ncid)
      if (stat /= ERA_OK) return

      dom%nlat = size(lat, kind=ik) ; dom%nlon = size(lon, kind=ik) ; dom%ncell = 1_ik
      dom%row = [best_j] ; dom%col = [best_i]
      dom%lat = [lat(best_j + 1_ik)] ; dom%lon = [lon(best_i + 1_ik)] ; dom%elevation = [elev(1, 1)]
      call group_by_chunk(dom)
   end subroutine era5land_select_site

   !=======================================================================================!
   !  BOX (§15.5): the valid cells whose centres lie in [south, north] x [west, east], in       !
   !  row-major order north to south and west to east. A box with west > east crosses the       !
   !  antimeridian: its columns run from west to the grid's east edge, then on from the west    !
   !  edge to east. nwse = [north, west, south, east] in degrees.                                !
   !=======================================================================================!
   subroutine era5land_select_box(static_path, nwse, dom, stat, land_fraction_min)
      character(len=*),   intent(in)  :: static_path
      real(wp),           intent(in)  :: nwse(4)
      type(met_cells_t), intent(out) :: dom
      integer(ik),        intent(out) :: stat
      !----- Keep only cells whose static land fraction is at least this (a region's selection      !
      !      rule, MEDS_POLYGON_RUNTIME_PLAN.md §5); absent, every valid cell is kept. ------------!
      real(wp), optional, intent(in)  :: land_fraction_min
      integer(c_int) :: ncid, st
      real(wp), allocatable :: lat(:), lon(:), valid(:,:), elev(:,:), lfrac(:,:)
      real(wp)    :: lf_min
      integer(ik), allocatable :: cols(:), col_index(:)
      real(wp)    :: west, east
      integer(ik) :: nlat, nlon, r0, r1, j, i, k, n

      call open_static(static_path, ncid, lat, lon, stat)
      if (stat /= ERA_OK) return
      nlat = size(lat, kind=ik) ; nlon = size(lon, kind=ik)
      r0 = -1_ik ; r1 = -1_ik
      do j = 0_ik, nlat - 1_ik
         if (lat(j + 1_ik) <= nwse(1) + GRID_TOL .and. lat(j + 1_ik) >= nwse(3) - GRID_TOL) then
            if (r0 < 0_ik) r0 = j
            r1 = j
         end if
      end do
      west = wrap_longitude(nwse(2), lon(1)) ; east = wrap_longitude(nwse(4), lon(1))
      !----- The 0-based column indices are a NAMED array: nvfortran 25.11 miscompiles an implied-do  !
      !      inside a pack inside an array constructor -- wrong indices at -O2, a segfault at -O0. ---!
      col_index = [(i, i = 0_ik, nlon - 1_ik)]
      if (west <= east) then
         cols = pack(col_index, lon >= west - GRID_TOL .and. lon <= east + GRID_TOL)
      else
         cols = [pack(col_index, lon >= west - GRID_TOL), pack(col_index, lon <= east + GRID_TOL)]
      end if
      if (r0 < 0_ik .or. size(cols) == 0) then
         st = nc_close(ncid) ; stat = ERA_ERR_NO_CELL ; return
      end if
      call read_2d(ncid, 'valid', r0, r1 - r0 + 1_ik, 0_ik, nlon, valid, stat)
      if (stat == ERA_OK) call read_2d(ncid, 'elevation', r0, r1 - r0 + 1_ik, 0_ik, nlon, elev, stat)
      if (stat == ERA_OK) call read_2d(ncid, 'land_fraction', r0, r1 - r0 + 1_ik, 0_ik, nlon, lfrac, stat)
      st = nc_close(ncid)
      if (stat /= ERA_OK) return

      lf_min = -1.0_wp ; if (present(land_fraction_min)) lf_min = land_fraction_min
      where (lfrac < lf_min) valid = 0.0_wp
      n = int(count(valid(cols + 1_ik, :) > 0.5_wp), ik)
      if (n == 0_ik) then ; stat = ERA_ERR_NO_CELL ; return ; end if
      dom%nlat = nlat ; dom%nlon = nlon ; dom%ncell = n
      allocate(dom%row(n), dom%col(n), dom%lat(n), dom%lon(n), dom%elevation(n), dom%land_fraction(n))
      k = 0_ik
      do j = r0, r1
         do i = 1_ik, size(cols, kind=ik)
            if (valid(cols(i) + 1_ik, j - r0 + 1_ik) < 0.5_wp) cycle
            k = k + 1_ik
            dom%row(k) = j ; dom%col(k) = cols(i)
            dom%lat(k) = lat(j + 1_ik) ; dom%lon(k) = lon(cols(i) + 1_ik)
            dom%elevation(k) = elev(cols(i) + 1_ik, j - r0 + 1_ik)
            dom%land_fraction(k) = lfrac(cols(i) + 1_ik, j - r0 + 1_ik)
         end do
      end do
      call group_by_chunk(dom)
   end subroutine era5land_select_box

   !----- A longitude moved into [lon0, lon0 + 360), the grid's own range. ------------------------!
   pure real(wp) function wrap_longitude(x, lon0) result(y)
      real(wp), intent(in) :: x, lon0
      y = lon0 + modulo(x - lon0, 360.0_wp)
   end function wrap_longitude

   !----- Group the domain's cells by 16 x 16 chunk with a stable counting sort, so each touched  !
   !      chunk column is read once per variable-month and the cells keep their domain order. Each !
   !      chunk's read box is the rows and columns its cells occupy: HDF5 decompresses the whole    !
   !      chunk either way, but a site then copies one cell instead of 256. ------------------------!
   subroutine group_by_chunk(dom)
      type(met_cells_t), intent(inout) :: dom
      integer(ik), allocatable :: id(:), cnt(:), next(:)
      integer(ik) :: ncx, ncy, c, k, m
      ncx = (dom%nlon + CHUNK - 1_ik) / CHUNK ; ncy = (dom%nlat + CHUNK - 1_ik) / CHUNK
      allocate(id(dom%ncell), cnt(0:ncx*ncy))
      id = (dom%row / CHUNK) * ncx + dom%col / CHUNK
      cnt = 0_ik
      do c = 1_ik, dom%ncell
         cnt(id(c)) = cnt(id(c)) + 1_ik
      end do
      dom%nchunk = int(count(cnt(0:ncx*ncy - 1_ik) > 0_ik), ik)
      allocate(dom%chunk_row(dom%nchunk), dom%chunk_col(dom%nchunk), dom%chunk_first(dom%nchunk + 1_ik))
      allocate(dom%chunk_nrow(dom%nchunk), dom%chunk_ncol(dom%nchunk))
      allocate(dom%by_chunk(dom%ncell), next(0:ncx*ncy - 1_ik))
      k = 0_ik ; m = 1_ik
      do c = 0_ik, ncx*ncy - 1_ik
         if (cnt(c) == 0_ik) cycle
         k = k + 1_ik
         dom%chunk_first(k) = m ; next(c) = m
         m = m + cnt(c)
      end do
      dom%chunk_first(dom%nchunk + 1_ik) = m
      do c = 1_ik, dom%ncell
         dom%by_chunk(next(id(c))) = c
         next(id(c)) = next(id(c)) + 1_ik
      end do
      do k = 1_ik, dom%nchunk
         associate (cells => dom%by_chunk(dom%chunk_first(k):dom%chunk_first(k + 1_ik) - 1_ik))
            dom%chunk_row(k)  = minval(dom%row(cells)) ; dom%chunk_col(k) = minval(dom%col(cells))
            dom%chunk_nrow(k) = maxval(dom%row(cells)) - dom%chunk_row(k) + 1_ik
            dom%chunk_ncol(k) = maxval(dom%col(cells)) - dom%chunk_col(k) + 1_ik
         end associate
      end do
   end subroutine group_by_chunk

   !=======================================================================================!
   !  MONTH LOAD (§15.3): every variable of one archive month for every domain cell. Each file's  !
   !  grid, stamps and units are checked against what the reader assumes; each touched chunk     !
   !  column is one hyperslab read of (month hours, <=16, <=16). `message` names the offending   !
   !  file, cell or hour on a nonzero status.                                                     !
   !=======================================================================================!
   subroutine era5land_load_month(template, data_path, dom, year, month, buf, stat, message)
      character(len=*),   intent(in)    :: template, data_path
      type(met_cells_t), intent(in)    :: dom
      integer(ik),        intent(in)    :: year, month
      type(met_month_t),  intent(inout) :: buf
      integer(ik),        intent(out)   :: stat
      character(len=*),   intent(out)   :: message
      character(len=MET_PATH_LEN) :: path
      character(len=64)  :: units
      integer(c_int)     :: ncid, vid, st
      integer(c_size_t)  :: start3(3), count3(3)
      real(c_float), allocatable :: block(:,:,:)
      real(c_double)     :: fill
      integer(ik) :: nt, v, k, m, c, nr, nc
      integer(ik) :: bad(3)

      stat = ERA_OK ; message = ''
      nt = era5land_month_hours(year, month)
      buf%year = 0_ik
      if (allocated(buf%values)) then
         if (any(shape(buf%values) /= [nt, dom%ncell, ERA_NVAR])) deallocate(buf%values)
      end if
      if (.not. allocated(buf%values)) allocate(buf%values(nt, dom%ncell, ERA_NVAR))

      do v = 1_ik, ERA_NVAR
         path = era5land_path(template, data_path, ERA_VAR_NAME(v), year, month)
         st = nc_open_f(trim(path), NC_NOWRITE, ncid)
         if (st /= NC_NOERR) then
            stat = ERA_ERR_OPEN ; message = trim(path) ; return
         end if
         call check_month_file(ncid, dom, year, month, nt, stat)
         if (stat /= ERA_OK) then
            st = nc_close(ncid) ; message = trim(path) ; return
         end if
         st = nc_inq_varid_f(ncid, trim(ERA_VAR_NAME(v)), vid)
         if (st == NC_NOERR) st = nc_get_att_text_f(ncid, vid, 'units', units)
         if (st /= NC_NOERR .or. trim(units) /= trim(ERA_VAR_UNITS(v))) then
            st = nc_close(ncid) ; stat = ERA_ERR_UNITS
            message = trim(path)//': units "'//trim(units)//'", expected "'//trim(ERA_VAR_UNITS(v))//'"'
            return
         end if
         !----- The reader recognises a missing value as NaN and nothing else (below). A file whose   !
         !      _FillValue is a number would slip that number into the forcing as data, so it is     !
         !      refused. The archive stores NaN, and a file without the attribute has no fill.  -----!
         st = nc_get_att_double_f(ncid, vid, '_FillValue', fill)
         if (st == NC_NOERR) then
            if (.not. ieee_is_nan(fill)) then
               st = nc_close(ncid) ; stat = ERA_ERR_FILL
               write(message, '(2a,es12.4,a)') trim(path), ': _FillValue is ', fill, ', not NaN'
               return
            end if
         end if
         do k = 1_ik, dom%nchunk
            nr = dom%chunk_nrow(k) ; nc = dom%chunk_ncol(k)
            allocate(block(nc, nr, nt))                      ! C order [time][lat][lon]
            start3 = [0_c_size_t, int(dom%chunk_row(k), c_size_t), int(dom%chunk_col(k), c_size_t)]
            count3 = [int(nt, c_size_t), int(nr, c_size_t), int(nc, c_size_t)]
            st = nc_get_vara_float(ncid, vid, start3, count3, block)   ! the archive stores float32
            if (st /= NC_NOERR) then
               deallocate(block) ; st = nc_close(ncid) ; stat = ERA_ERR_OPEN
               message = trim(path)//': read failed' ; return
            end if
            do m = dom%chunk_first(k), dom%chunk_first(k + 1_ik) - 1_ik
               c = dom%by_chunk(m)
               buf%values(:, c, v) = block(dom%col(c) - dom%chunk_col(k) + 1_ik,                     &
                                           dom%row(c) - dom%chunk_row(k) + 1_ik, :)
            end do
            deallocate(block)
         end do
         st = nc_close(ncid)
      end do

      !----- MEDS never gap-fills (§5.5): the static mask promises data in every selected cell.    !
      !      ieee_is_nan, not x /= x: a vectorized self-comparison may signal "invalid" on a quiet   !
      !      NaN, which a -fpe0 build traps before the message below can name the cell. The first NaN  !
      !      is found through an integer mask: nvfortran 25.11 has no findloc for LOGICAL arrays.    !
      if (any(ieee_is_nan(buf%values))) then
         bad = int(findloc(merge(1, 0, ieee_is_nan(buf%values)), 1), ik)
         write(message, '(a,a,i4.4,a,i2.2,a,i0,a,f9.3,a,f9.3,a)') trim(ERA_VAR_NAME(bad(3))), ' ',         &
               year, '-', month, ' hour ', bad(1), ' at cell (', dom%lat(bad(2)), ', ', dom%lon(bad(2)), ')'
         stat = ERA_ERR_NAN ; return
      end if
      buf%year = year ; buf%month = month ; buf%nt = nt
   end subroutine era5land_load_month

   !----- A month file must have the static grid and exactly the month's end-stamped hours. -------!
   subroutine check_month_file(ncid, dom, year, month, nt, stat)
      integer(c_int),     intent(in)  :: ncid
      type(met_cells_t), intent(in)  :: dom
      integer(ik),        intent(in)  :: year, month, nt
      integer(ik),        intent(out) :: stat
      integer(c_int)     :: st, vid
      integer(c_size_t)  :: n, start1(1), count1(1)
      character(len=256) :: units
      type(meds_time_t)  :: base
      real(c_double), allocatable :: tsec(:)
      real(wp)           :: first
      logical            :: ok
      integer(ik)        :: k

      stat = ERA_ERR_GRID
      st = nc_inq_dimlen_f(ncid, 'lat', n) ; if (st /= NC_NOERR .or. int(n, ik) /= dom%nlat) return
      st = nc_inq_dimlen_f(ncid, 'lon', n) ; if (st /= NC_NOERR .or. int(n, ik) /= dom%nlon) return
      stat = ERA_ERR_TIME
      st = nc_inq_dimlen_f(ncid, 'time', n) ; if (st /= NC_NOERR .or. int(n, ik) /= nt) return
      st = nc_inq_varid_f(ncid, 'time', vid) ; if (st /= NC_NOERR) return
      st = nc_get_att_text_f(ncid, vid, 'units', units) ; if (st /= NC_NOERR) return
      k = int(index(units, 'since'), ik) ; if (k == 0_ik) return
      call time_from_string(adjustl(units(k + 5_ik:)), base, ok) ; if (.not. ok) return
      allocate(tsec(nt))
      start1 = 0_c_size_t ; count1 = int(nt, c_size_t)
      st = nc_get_vara_double(ncid, vid, start1, count1, tsec) ; if (st /= NC_NOERR) return
      first = seconds_between(base, meds_time_t(year, month, 1_ik, 1_ik, 0_ik, 0_ik))
      do k = 1_ik, nt
         if (abs(tsec(k) - (first + 3600.0_wp * real(k - 1_ik, wp))) > STAMP_TOL) return
      end do
      stat = ERA_OK
   end subroutine check_month_file

   !----- Open the static file and read its (regular, 1-D) lat and lon axes. -----------------------!
   subroutine open_static(path, ncid, lat, lon, stat)
      character(len=*),      intent(in)  :: path
      integer(c_int),        intent(out) :: ncid
      real(wp), allocatable, intent(out) :: lat(:), lon(:)
      integer(ik),           intent(out) :: stat
      integer(c_int) :: st
      stat = ERA_ERR_OPEN
      st = nc_open_f(trim(path), NC_NOWRITE, ncid) ; if (st /= NC_NOERR) return
      stat = ERA_ERR_GRID
      call read_axis(ncid, 'lat', lat, st)
      if (st == NC_NOERR) call read_axis(ncid, 'lon', lon, st)
      if (st == NC_NOERR) then
         if (size(lat) >= 2 .and. size(lon) >= 2) then            ! a regular grid (§14.2)
            if (all(abs((lat(2:) - lat(:size(lat) - 1)) - (lat(2) - lat(1))) <= GRID_TOL) .and.        &
                all(abs((lon(2:) - lon(:size(lon) - 1)) - (lon(2) - lon(1))) <= GRID_TOL)) stat = ERA_OK
         end if
      end if
      if (stat /= ERA_OK) st = nc_close(ncid)
   end subroutine open_static

   subroutine read_axis(ncid, name, x, st)
      integer(c_int),        intent(in)  :: ncid
      character(len=*),      intent(in)  :: name
      real(wp), allocatable, intent(out) :: x(:)
      integer(c_int),        intent(out) :: st
      integer(c_int)    :: vid
      integer(c_size_t) :: n, start1(1), count1(1)
      st = nc_inq_dimlen_f(ncid, name, n) ; if (st /= NC_NOERR) return
      st = nc_inq_varid_f(ncid, name, vid) ; if (st /= NC_NOERR) return
      allocate(x(n))
      start1 = 0_c_size_t ; count1 = n
      st = nc_get_vara_double(ncid, vid, start1, count1, x)
   end subroutine read_axis

   !----- rows r0 .. r0+nr-1 and columns c0 .. c0+nc-1 (0-based) of a (lat, lon) variable, as x(lon, lat). !
   subroutine read_2d(ncid, name, r0, nr, c0, nc, x, stat)
      integer(c_int),        intent(in)  :: ncid
      character(len=*),      intent(in)  :: name
      integer(ik),           intent(in)  :: r0, nr, c0, nc
      real(wp), allocatable, intent(out) :: x(:,:)
      integer(ik),           intent(out) :: stat
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start2(2), count2(2)
      stat = ERA_ERR_GRID
      st = nc_inq_varid_f(ncid, name, vid) ; if (st /= NC_NOERR) return
      allocate(x(nc, nr))
      start2 = [int(r0, c_size_t), int(c0, c_size_t)] ; count2 = [int(nr, c_size_t), int(nc, c_size_t)]
      st = nc_get_vara_double(ncid, vid, start2, count2, x) ; if (st /= NC_NOERR) return
      stat = ERA_OK
   end subroutine read_2d

end module meds_era5land_reader
