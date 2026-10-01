! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_test_era5land_archive -- a tiny synthetic ED_ERA5land archive for the tests, written       !
! through the netCDF C API: a 4 x 18 grid (10 deg x 20 deg) spanning the antimeridian and two     !
! 16-cell chunk columns wide, every month of 2021 and January 2022 (so a recycle window may start  !
! after January 1st), two no-data cells, a mostly-water column, and values that are analytic       !
! functions of cell and hour, so every check has a known answer. Shared by test_met_era5land (the  !
! reader) and test_region (a region against site runs).                                            !
!==========================================================================================!
module meds_test_era5land_archive
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use meds_kinds,           only : wp, sp, ik
   use meds_time,            only : meds_time_t, seconds_between
   use meds_era5land_reader, only : era5land_path, era5land_default_template, era5land_month_hours, &
                                    ERA_NVAR, ERA_VAR_NAME, ERA_VAR_UNITS
   use meds_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   implicit none
   private

   public :: NLAT, NLON, LAT0, DLAT, LON0, DLON, T0
   public :: is_valid, field, air_temperature, hour_index, write_archive, write_static, write_month

   integer(ik), parameter :: NLAT = 4_ik, NLON = 18_ik
   real(wp),    parameter :: LAT0 = 45.0_wp, DLAT = -10.0_wp, LON0 = -180.0_wp, DLON = 20.0_wp
   real(wp),    parameter :: TWO_PI = 6.283185307179586_wp
   !----- The archive's first stamp; g below counts hours from it (0-based). ---------------------!
   type(meds_time_t), parameter :: T0 = meds_time_t(2021_ik, 1_ik, 1_ik, 1_ik, 0_ik, 0_ik)

contains

   !----- No data at (row 0, col 5) and (row 2, col 0); everything else is valid. ----------------!
   pure logical function is_valid(r, c)
      integer(ik), intent(in) :: r, c
      is_valid = .not. ((r == 0_ik .and. c == 5_ik) .or. (r == 2_ik .and. c == 0_ik))
   end function is_valid

   !----- The value of variable v at grid (row r, col c) and hour g, as the archive stores it. ---!
   pure real(wp) function field(v, r, c, g) result(x)
      integer(ik), intent(in) :: v, r, c, g
      real(wp) :: hr, day
      hr = real(g, wp) ; day = real(modulo(g, 24_ik), wp)
      select case (v)
      case (1) ; x = air_temperature(r, c, hr)                                                   ! Tair
      case (2) ; x = air_temperature(r, c, hr) - 4.0_wp - 0.5_wp * sin(TWO_PI * hr / 24.0_wp)    ! Tdew
      case (3) ; x = 90000.0_wp + 100.0_wp * real(r, wp) + 10.0_wp * real(c, wp) + day         ! PSurf
      case (4) ; x = 2.0_wp + 3.0_wp * sin(TWO_PI * hr / 24.0_wp) + 0.1_wp * real(c, wp)       ! u10
      case (5) ; x = -1.0_wp + 2.0_wp * cos(TWO_PI * hr / 17.0_wp) - 0.1_wp * real(r, wp)      ! v10
      case (6) ; x = merge(1.0e-4_wp, 0.0_wp, modulo(g, 7_ik) == 3_ik)                          ! Rainf
      case (7) ; x = max(0.0_wp, 600.0_wp * sin(TWO_PI * (day - 6.0_wp) / 24.0_wp))              ! SWdown
      case default ; x = 300.0_wp + 0.001_wp * hr + real(c, wp)                                  ! LWdown
      end select
      x = real(real(x, sp), wp)                              ! the archive holds float32
   end function field

   pure real(wp) function air_temperature(r, c, hr) result(x)
      integer(ik), intent(in) :: r, c
      real(wp),    intent(in) :: hr
      x = 260.0_wp + 0.002_wp * hr + real(r, wp) + 0.1_wp * real(c, wp)
   end function air_temperature

   !----- Hours from T0 to t. -----------------------------------------------------------------------!
   integer(ik) function hour_index(t) result(g)
      type(meds_time_t), intent(in) :: t
      g = nint(seconds_between(T0, t) / 3600.0_wp, ik)
   end function hour_index

   !=======================================================================================!
   !  The synthetic archive.                                                                     !
   !=======================================================================================!
   subroutine write_archive(dir)
      character(len=*), intent(in) :: dir
      integer(ik) :: m, v
      call execute_command_line('mkdir -p '//dir)
      call write_static(dir//'/ED_ERA5land_static.nc', .false.)
      call write_static(dir//'/static_allvalid.nc', .true.)
      do m = 1_ik, 12_ik
         do v = 1_ik, ERA_NVAR
            call write_month(dir, v, 2021_ik, m)
         end do
      end do
      do v = 1_ik, ERA_NVAR
         call write_month(dir, v, 2022_ik, 1_ik)
      end do
   end subroutine write_archive

   subroutine define_grid(ncid, dlat_id, dlon_id)
      integer(c_int), intent(in)  :: ncid
      integer(c_int), intent(out) :: dlat_id, dlon_id
      integer(c_int) :: st, vla, vlo, dims1(1)
      st = nc_def_dim_f(ncid, 'lat', int(NLAT, c_size_t), dlat_id) ; call nc_check(st, 'def lat dim')
      st = nc_def_dim_f(ncid, 'lon', int(NLON, c_size_t), dlon_id) ; call nc_check(st, 'def lon dim')
      dims1 = dlat_id
      st = nc_def_var_f(ncid, 'lat', NC_DOUBLE, 1, dims1, vla) ; call nc_check(st, 'def lat')
      dims1 = dlon_id
      st = nc_def_var_f(ncid, 'lon', NC_DOUBLE, 1, dims1, vlo) ; call nc_check(st, 'def lon')
   end subroutine define_grid

   subroutine put_grid(ncid)
      integer(c_int), intent(in) :: ncid
      integer(c_int)    :: st, vid
      integer(c_size_t) :: start1(1), count1(1)
      real(c_double)    :: lat(NLAT), lon(NLON)
      integer(ik) :: k
      lat = [(LAT0 + DLAT * real(k, wp), k = 0_ik, NLAT - 1_ik)]
      lon = [(LON0 + DLON * real(k, wp), k = 0_ik, NLON - 1_ik)]
      start1 = 0_c_size_t
      st = nc_inq_varid_f(ncid, 'lat', vid) ; count1 = int(NLAT, c_size_t)
      st = nc_put_vara_double(ncid, vid, start1, count1, lat) ; call nc_check(st, 'put lat')
      st = nc_inq_varid_f(ncid, 'lon', vid) ; count1 = int(NLON, c_size_t)
      st = nc_put_vara_double(ncid, vid, start1, count1, lon) ; call nc_check(st, 'put lon')
   end subroutine put_grid

   subroutine write_static(path, all_valid)
      character(len=*), intent(in) :: path
      logical,          intent(in) :: all_valid
      integer(c_int)    :: st, ncid, dla, dlo, vv, ve, vf, dims2(2)
      integer(c_size_t) :: start2(2), count2(2)
      real(c_double)    :: valid(NLON, NLAT), elev(NLON, NLAT), lfrac(NLON, NLAT)
      integer(ik) :: r, c
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create static')
      call define_grid(ncid, dla, dlo)
      dims2 = [dla, dlo]
      st = nc_def_var_f(ncid, 'valid', NC_BYTE, 2, dims2, vv) ; call nc_check(st, 'def valid')
      st = nc_def_var_f(ncid, 'elevation', NC_FLOAT, 2, dims2, ve) ; call nc_check(st, 'def elevation')
      st = nc_def_var_f(ncid, 'land_fraction', NC_FLOAT, 2, dims2, vf) ; call nc_check(st, 'def land_fraction')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef static')
      call put_grid(ncid)
      do r = 0_ik, NLAT - 1_ik
         do c = 0_ik, NLON - 1_ik
            valid(c + 1_ik, r + 1_ik) = merge(1.0_c_double, 0.0_c_double, all_valid .or. is_valid(r, c))
            elev(c + 1_ik, r + 1_ik)  = 100.0_c_double * real(r, c_double) + real(c, c_double)
            lfrac(c + 1_ik, r + 1_ik) = merge(0.3_c_double, 1.0_c_double, c == 1_ik)   ! column 1: mostly water
         end do
      end do
      start2 = 0_c_size_t ; count2 = [int(NLAT, c_size_t), int(NLON, c_size_t)]
      st = nc_put_vara_double(ncid, vv, start2, count2, valid) ; call nc_check(st, 'put valid')
      st = nc_put_vara_double(ncid, ve, start2, count2, elev) ; call nc_check(st, 'put elevation')
      st = nc_put_vara_double(ncid, vf, start2, count2, lfrac) ; call nc_check(st, 'put land_fraction')
      st = nc_close(ncid) ; call nc_check(st, 'close static')
   end subroutine write_static

   subroutine write_month(dir, v, y, m, fill)
      character(len=*), intent(in) :: dir
      integer(ik),      intent(in) :: v, y, m
      real(c_double), optional, intent(in) :: fill   !< a numeric _FillValue to declare (a bad file)
      character(len=:), allocatable :: path
      integer(c_int)    :: st, ncid, dt, dla, dlo, vt, vx, dims1(1), dims3(3)
      integer(c_size_t) :: start1(1), count1(1), start3(3), count3(3)
      real(c_double), allocatable :: tsec(:), dat(:,:,:)
      real(c_double) :: nan
      integer(ik) :: nt, g0, h, r, c
      nan = ieee_value(0.0_c_double, ieee_quiet_nan)
      nt = era5land_month_hours(y, m)
      g0 = hour_index(meds_time_t(y, m, 1_ik, 1_ik))
      path = trim(era5land_path(era5land_default_template(), dir, ERA_VAR_NAME(v), y, m))
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create '//path)
      st = nc_def_dim_f(ncid, 'time', int(nt, c_size_t), dt) ; call nc_check(st, 'def time dim')
      call define_grid(ncid, dla, dlo)
      dims1 = dt
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt) ; call nc_check(st, 'def time')
      st = nc_put_att_text_f(ncid, vt, 'units', 33_c_size_t, 'seconds since 1970-01-01 00:00:00')
      call nc_check(st, 'time units')
      dims3 = [dt, dla, dlo]
      st = nc_def_var_f(ncid, trim(ERA_VAR_NAME(v)), NC_FLOAT, 3, dims3, vx) ; call nc_check(st, 'def var')
      st = nc_put_att_text_f(ncid, vx, 'units', int(len_trim(ERA_VAR_UNITS(v)), c_size_t), trim(ERA_VAR_UNITS(v)))
      call nc_check(st, 'var units')
      if (present(fill)) then
         st = nc_put_att_double_f(ncid, vx, '_FillValue', NC_FLOAT, fill) ; call nc_check(st, 'var fill')
      end if
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', 3_c_size_t, 'end') ; call nc_check(st, 'avg_convention')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef month')
      call put_grid(ncid)
      allocate(tsec(nt), dat(NLON, NLAT, nt))
      do h = 1_ik, nt
         tsec(h) = seconds_between(meds_time_t(1970_ik, 1_ik, 1_ik), T0) + 3600.0_wp * real(g0 + h - 1_ik, wp)
         do r = 0_ik, NLAT - 1_ik
            do c = 0_ik, NLON - 1_ik
               dat(c + 1_ik, r + 1_ik, h) = merge(field(v, r, c, g0 + h - 1_ik), nan, is_valid(r, c))
            end do
         end do
      end do
      start1 = 0_c_size_t ; count1 = int(nt, c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec) ; call nc_check(st, 'put time')
      start3 = 0_c_size_t ; count3 = [int(nt, c_size_t), int(NLAT, c_size_t), int(NLON, c_size_t)]
      st = nc_put_vara_double(ncid, vx, start3, count3, dat)
      call nc_check(st, 'put '//path)
      st = nc_close(ncid) ; call nc_check(st, 'close '//path)
   end subroutine write_month

end module meds_test_era5land_archive
