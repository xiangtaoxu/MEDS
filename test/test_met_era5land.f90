! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_met_era5land -- the ED_ERA5land archive backend of the met reader (MEDS_FORCING_DESIGN  !
! .md §15.6). A tiny synthetic archive is written through the netCDF C API: a 4 x 18 grid      !
! (10 deg x 20 deg) spanning the antimeridian and two 16-cell chunk columns wide, every month   !
! of 2021, two no-data cells, and values that are analytic functions of cell and hour, so every  !
! check below has a known answer. It is read back through met_open / met_advance / met_instant   !
! and through the reader's own selection and month-load routines.                                !
!==========================================================================================!
program test_met_era5land
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use meds_test_assert,     only : check, check_true, test_report
   use meds_kinds,           only : wp, sp, ik
   use meds_time,            only : meds_time_t, seconds_between
   use meds_forcing_config,  only : forcing_config_t, MET_BACKEND_ERA5LAND, METAVG_END,          &
                                    SWPART_CLEARIDX, SWPART_PASSTHROUGH, CLAMP_ERROR
   use meds_forcing_types,   only : met_driver_t, met_forcing_t, met_domain_t, met_month_t
   use meds_forcing_kernels, only : dewpoint_to_specific_humidity, wind_log_profile
   use meds_met_driver,      only : met_open, met_advance, met_instant, met_close,             &
                                    MET_OK, MET_ERR_ARCHIVE, MET_ERR_ATTR_MISMATCH
   use meds_era5land_reader, only : era5land_path, era5land_default_template, era5land_select_site, &
                                    era5land_select_box, era5land_load_month, era5land_month_hours, &
                                    ERA_NVAR, ERA_VAR_NAME, ERA_VAR_UNITS, ERA_OK, ERA_ERR_NAN,     &
                                    ERA_ERR_NO_CELL, ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10
   use meds_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   implicit none

   character(len=*), parameter :: DIR = 'era5land_test_tmp'
   integer(ik), parameter :: NLAT = 4_ik, NLON = 18_ik
   real(wp),    parameter :: LAT0 = 45.0_wp, DLAT = -10.0_wp, LON0 = -180.0_wp, DLON = 20.0_wp
   real(wp),    parameter :: TWO_PI = 6.283185307179586_wp
   !----- The archive's first stamp; g below counts hours from it (0-based). ---------------------!
   type(meds_time_t), parameter :: T0 = meds_time_t(2021_ik, 1_ik, 1_ik, 1_ik, 0_ik, 0_ik)

   call execute_command_line('mkdir -p '//DIR)
   call write_archive()

   call test_paths()
   call test_site_selection()
   call test_box_selection()
   call test_nan_rejected()
   call test_driver_months()
   call test_driver_recycle()
   call test_driver_rejections()

   call test_report('test_met_era5land')

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
   !  1. Templates.                                                                              !
   !=======================================================================================!
   subroutine test_paths()
      print '(a)', '-- era5land 1: path templates --'
      call check_true('default template',                                                        &
         trim(era5land_path(era5land_default_template(), '/data/x', 'Tair', 2021_ik, 3_ik))       &
         == '/data/x/ED_ERA5land_Tair_202103.nc')
      call check_true('custom template with a repeated token',                                   &
         trim(era5land_path('{data_path}/{yyyy}/{var}_{yyyy}{mm}.nc', '/d', 'u10', 2021_ik, 12_ik)) &
         == '/d/2021/u10_202112.nc')
   end subroutine test_paths

   !=======================================================================================!
   !  2. Site -> cell: its own cell when valid, else the nearest valid one within the limit.     !
   !=======================================================================================!
   subroutine test_site_selection()
      type(met_domain_t) :: dom
      real(wp)    :: d
      integer(ik) :: st
      print '(a)', '-- era5land 2: site to cell --'
      call era5land_select_site(DIR//'/ED_ERA5land_static.nc', 35.0_wp, -120.0_wp, 50.0_wp, dom, d, st)
      call check_true('a site on a valid cell takes that cell', st == ERA_OK .and. dom%ncell == 1_ik &
                      .and. dom%row(1) == 1_ik .and. dom%col(1) == 3_ik)
      call check('its distance is 0', d, 0.0_wp, 1.0e-6_wp)
      call check('its elevation comes from the static file', dom%elevation(1), 103.0_wp, 1.0e-6_wp)
      call check_true('a site reads a 1 x 1 box', dom%nchunk == 1_ik .and. dom%chunk_nrow(1) == 1_ik  &
                      .and. dom%chunk_ncol(1) == 1_ik .and. dom%chunk_row(1) == 1_ik .and. dom%chunk_col(1) == 3_ik)

      !----- (0,5) has no data. Its nearest valid neighbour is (1,5), 10 deg of latitude south   !
      !      (1112 km), not (0,4) or (0,6), 20 deg of longitude away at 45 N (about 1570 km). ------!
      call era5land_select_site(DIR//'/ED_ERA5land_static.nc', 45.0_wp, -80.0_wp, 2000.0_wp, dom, d, st)
      call check_true('a no-data site takes the nearest valid cell', st == ERA_OK                &
                      .and. dom%row(1) == 1_ik .and. dom%col(1) == 5_ik)
      call check('and reports its distance', d, 1111.95_wp, 0.1_wp)
      call era5land_select_site(DIR//'/ED_ERA5land_static.nc', 45.0_wp, -80.0_wp, 1000.0_wp, dom, d, st)
      call check_true('no valid cell within the limit is an error', st == ERA_ERR_NO_CELL)
   end subroutine test_site_selection

   !=======================================================================================!
   !  3. Box: valid cells in row-major order, across the antimeridian; and a month load of      !
   !     them, which reads one chunk column per touched chunk.                                   !
   !=======================================================================================!
   subroutine test_box_selection()
      type(met_domain_t) :: dom
      type(met_month_t)  :: buf
      character(len=512) :: msg
      real(wp)    :: across(4), inside(4)
      integer(ik) :: st, g
      print '(a)', '-- era5land 3: box selection and month load --'
      across = [40.0_wp, 150.0_wp, 20.0_wp, -150.0_wp] ; inside = [50.0_wp, -130.0_wp, 10.0_wp, -90.0_wp]
      !----- [N, W, S, E] = [40, 150, 20, -150]: rows 35 and 25 N; columns 160 E, then -180 and  !
      !      -160 across 180. (2,0) has no data, so 5 cells.                                        !
      call era5land_select_box(DIR//'/ED_ERA5land_static.nc', across, dom, st)
      call check_true('box across 180: 5 valid cells', st == ERA_OK .and. dom%ncell == 5_ik)
      call check_true('row-major, west to east across 180',                                      &
                      all(dom%row == [1_ik, 1_ik, 1_ik, 2_ik, 2_ik]) .and.                        &
                      all(dom%col == [17_ik, 0_ik, 1_ik, 17_ik, 1_ik]))
      call check_true('two chunk columns touched', dom%nchunk == 2_ik)
      !----- Each chunk reads only its cells' box: rows 1-2 x cols 0-1 west of 180, rows 1-2 x col 17 east. !
      call check_true('each chunk reads only the box of its cells',                                &
                      all(dom%chunk_row == [1_ik, 1_ik]) .and. all(dom%chunk_col == [0_ik, 17_ik]) .and. &
                      all(dom%chunk_nrow == [2_ik, 2_ik]) .and. all(dom%chunk_ncol == [2_ik, 1_ik]))
      call era5land_select_box(DIR//'/ED_ERA5land_static.nc', inside, dom, st)
      call check_true('box not across 180: 4 rows x 2 columns', st == ERA_OK .and. dom%ncell == 8_ik)

      call era5land_select_box(DIR//'/ED_ERA5land_static.nc', across, dom, st)
      call era5land_load_month(era5land_default_template(), DIR, dom, 2021_ik, 2_ik, buf, st, msg)
      call check_true('February loads', st == ERA_OK .and. buf%nt == 672_ik)
      g = hour_index(meds_time_t(2021_ik, 2_ik, 1_ik, 1_ik))
      call check('first hour, cell (2,17), Tair', real(buf%values(1, 4, ERA_TAIR), wp),          &
                 field(ERA_TAIR, 2_ik, 17_ik, g), 1.0e-9_wp)
      call check('last hour, cell (1,0), u10', real(buf%values(672, 2, ERA_U10), wp),            &
                 field(ERA_U10, 1_ik, 0_ik, g + 671_ik), 1.0e-9_wp)
      call check('mid-month, cell (1,1), PSurf', real(buf%values(300, 3, ERA_PSURF), wp),         &
                 field(ERA_PSURF, 1_ik, 1_ik, g + 299_ik), 1.0e-9_wp)
   end subroutine test_box_selection

   !=======================================================================================!
   !  4. MEDS never gap-fills: a selected cell with no data is an error, never a fill. The      !
   !     second static file claims data at (0,5), where the month files hold NaN.                !
   !=======================================================================================!
   subroutine test_nan_rejected()
      type(met_domain_t) :: dom
      type(met_month_t)  :: buf
      character(len=512) :: msg
      real(wp)    :: d
      integer(ik) :: st
      print '(a)', '-- era5land 4: NaN in a selected cell --'
      call era5land_select_site(DIR//'/static_allvalid.nc', 45.0_wp, -80.0_wp, 50.0_wp, dom, d, st)
      call era5land_load_month(era5land_default_template(), DIR, dom, 2021_ik, 1_ik, buf, st, msg)
      call check_true('a NaN in a selected cell is rejected', st == ERA_ERR_NAN)
      call check_true('and the message names the variable', index(msg, 'Tair') > 0)
   end subroutine test_nan_rejected

   !----- A reader config on the synthetic archive at the valid cell (1,3), 35 N 120 W. ----------!
   function archive_config() result(fc)
      type(forcing_config_t) :: fc
      fc%backend = MET_BACKEND_ERA5LAND ; fc%data_path = DIR ; fc%max_distance_km = 50.0_wp
      fc%dt_forcing = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude_deg = 35.0_wp ; fc%longitude_deg = -120.0_wp ; fc%utc_offset_h = 0.0_wp
      fc%apply_solar_longitude = .true. ; fc%recycle = .false. ; fc%start_clamp = CLAMP_ERROR
   end function archive_config

   !=======================================================================================!
   !  5. The driver over a month boundary: conversions, the seam bracket, the wind vector, and  !
   !     the static elevation.                                                                   !
   !=======================================================================================!
   subroutine test_driver_months()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      real(wp)    :: u1, u2, v1, v2, factor
      print '(a)', '-- era5land 5: driver across a month boundary --'
      fc = archive_config()
      call met_open(drv, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call check_true('opens (MET_OK)', st == MET_OK)
      call check_true('reads January and February: 1416 hourly records', drv%nrec == 1416_ik)
      call check('grid_elevation from the static file', drv%fcfg%grid_elevation_m, 103.0_wp, 1.0e-6_wp)
      call check_true('the archive supplies the wind vector', drv%has_wind_vector)

      !----- At a stamp every weight is 0, so the record comes through untouched. ---------------!
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      call met_advance(drv, t) ; met = met_instant(drv, t)
      call check('Tair at a stamp', met%tair_k, field(ERA_TAIR, 1_ik, 3_ik, g), 1.0e-9_wp)
      call check('Qair = the kernel of the stored Tdew and PSurf', met%qair,                      &
                 dewpoint_to_specific_humidity(field(ERA_TDEW, 1_ik, 3_ik, g),                    &
                                               field(ERA_PSURF, 1_ik, 3_ik, g)), 1.0e-15_wp)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; v1 = field(ERA_V10, 1_ik, 3_ik, g)
      call check('wind_u = u10', met%wind_u, u1, 1.0e-12_wp)
      call check('wind_v = v10', met%wind_v, v1, 1.0e-12_wp)
      call check('speed from the components', met%wind, sqrt(u1**2 + v1**2), 1.0e-12_wp)
      call check_true('the instant carries the vector', met%has_wind_vector)

      !----- Inside January's file: 31 Jan 23:30 sits between 23:00 and 00:00 (both January). --!
      t = meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik, 30_ik) ; g = hour_index(meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik))
      call met_advance(drv, t) ; met = met_instant(drv, t)
      call check('31 Jan 23:30 interpolates inside the January file', met%tair_k,                 &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)

      !----- The seam: 1 Feb 00:00 is January's last record, 01:00 February's first. -----------!
      t = meds_time_t(2021_ik, 2_ik, 1_ik, 0_ik, 30_ik) ; g = hour_index(meds_time_t(2021_ik, 2_ik, 1_ik))
      call met_advance(drv, t) ; met = met_instant(drv, t)
      call check('1 Feb 00:30 brackets January''s last and February''s first record', met%tair_k, &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; u2 = field(ERA_U10, 1_ik, 3_ik, g + 1_ik)
      v1 = field(ERA_V10, 1_ik, 3_ik, g) ; v2 = field(ERA_V10, 1_ik, 3_ik, g + 1_ik)
      call check('the vector interpolates linearly (u)', met%wind_u, 0.5_wp * (u1 + u2), 1.0e-12_wp)
      call check('the vector interpolates linearly (v)', met%wind_v, 0.5_wp * (v1 + v2), 1.0e-12_wp)
      call check_true('the energy-form speed is at least the vector''s length',                    &
                      met%wind >= sqrt(met%wind_u**2 + met%wind_v**2) - 1.0e-12_wp)
      call met_close(drv)

      !----- The height correction scales both components by the speed's factor. ---------------!
      fc%apply_wind_profile = .true. ; fc%wind_meas_height = 10.0_wp ; fc%reference_height = 40.0_wp
      fc%wind_roughness_z0 = 0.1_wp
      factor = log(400.0_wp) / log(100.0_wp)
      call met_open(drv, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      call met_advance(drv, t) ; met = met_instant(drv, t)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; v1 = field(ERA_V10, 1_ik, 3_ik, g)
      call check('height correction on u', met%wind_u, u1 * factor, 1.0e-12_wp)
      call check('height correction on v', met%wind_v, v1 * factor, 1.0e-12_wp)
      call check('direction preserved', atan2(met%wind_v, met%wind_u), atan2(v1, u1), 1.0e-12_wp)
      call met_close(drv)
   end subroutine test_driver_months

   !=======================================================================================!
   !  6. Recycling across the months of a declared window, including the year-end seam, which   !
   !     brackets December's last record (2022-01-01 00:00) and January's first.                 !
   !=======================================================================================!
   subroutine test_driver_recycle()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      print '(a)', '-- era5land 6: recycling across months --'
      fc = archive_config()
      fc%recycle = .true.
      fc%recycle_start = T0 ; fc%recycle_end = meds_time_t(2022_ik, 1_ik, 1_ik, 1_ik)
      call met_open(drv, fc, stat=st)
      call check_true('a whole-year window opens without a run period', st == MET_OK)
      call check_true('the window is 8760 records', drv%nrec == 8760_ik)
      t = meds_time_t(2035_ik, 6_ik, 15_ik, 12_ik) ; g = hour_index(meds_time_t(2021_ik, 6_ik, 15_ik, 12_ik))
      call met_advance(drv, t) ; met = met_instant(drv, t)
      call check('2035-06-15 12:00 reads 2021-06-15 12:00', met%tair_k, field(ERA_TAIR, 1_ik, 3_ik, g), 1.0e-9_wp)
      t = meds_time_t(2036_ik, 1_ik, 1_ik, 0_ik, 30_ik)
      call met_advance(drv, t) ; met = met_instant(drv, t)
      call check('the year-end seam brackets December''s last and January''s first record',       &
                 met%tair_k, 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, 8759_ik) + field(ERA_TAIR, 1_ik, 3_ik, 0_ik)), &
                 1.0e-9_wp)
      call met_close(drv)
   end subroutine test_driver_recycle

   !=======================================================================================!
   !  7. Rejected at open: a month the archive lacks, total shortwave read as components, and  !
   !     a site with no valid cell in reach.                                                     !
   !=======================================================================================!
   subroutine test_driver_rejections()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      integer(ik) :: st
      print '(a)', '-- era5land 7: rejections at open --'
      fc = archive_config()
      call met_open(drv, fc, stat=st, run_start=meds_time_t(2021_ik, 11_ik, 1_ik),               &
                    run_end=meds_time_t(2022_ik, 2_ik, 1_ik))
      call check_true('a run past the archive''s last month is rejected', st == MET_ERR_ARCHIVE)
      fc%sw_partition = SWPART_PASSTHROUGH
      call met_open(drv, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('passthrough against total shortwave is rejected', st == MET_ERR_ATTR_MISMATCH)
      fc = archive_config()
      fc%latitude_deg = 45.0_wp ; fc%longitude_deg = -80.0_wp
      call met_open(drv, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('a site with no valid cell within max_distance_km is rejected', st == MET_ERR_ARCHIVE)
   end subroutine test_driver_rejections

   !=======================================================================================!
   !  The synthetic archive.                                                                     !
   !=======================================================================================!
   subroutine write_archive()
      integer(ik) :: m, v
      call write_static(DIR//'/ED_ERA5land_static.nc', .false.)
      call write_static(DIR//'/static_allvalid.nc', .true.)
      do m = 1_ik, 12_ik
         do v = 1_ik, ERA_NVAR
            call write_month(v, m)
         end do
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
      integer(c_int)    :: st, ncid, dla, dlo, vv, ve, dims2(2)
      integer(c_size_t) :: start2(2), count2(2)
      real(c_double)    :: valid(NLON, NLAT), elev(NLON, NLAT)
      integer(ik) :: r, c
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'create static')
      call define_grid(ncid, dla, dlo)
      dims2 = [dla, dlo]
      st = nc_def_var_f(ncid, 'valid', NC_BYTE, 2, dims2, vv) ; call nc_check(st, 'def valid')
      st = nc_def_var_f(ncid, 'elevation', NC_FLOAT, 2, dims2, ve) ; call nc_check(st, 'def elevation')
      st = nc_enddef(ncid) ; call nc_check(st, 'enddef static')
      call put_grid(ncid)
      do r = 0_ik, NLAT - 1_ik
         do c = 0_ik, NLON - 1_ik
            valid(c + 1_ik, r + 1_ik) = merge(1.0_c_double, 0.0_c_double, all_valid .or. is_valid(r, c))
            elev(c + 1_ik, r + 1_ik)  = 100.0_c_double * real(r, c_double) + real(c, c_double)
         end do
      end do
      start2 = 0_c_size_t ; count2 = [int(NLAT, c_size_t), int(NLON, c_size_t)]
      st = nc_put_vara_double(ncid, vv, start2, count2, valid) ; call nc_check(st, 'put valid')
      st = nc_put_vara_double(ncid, ve, start2, count2, elev) ; call nc_check(st, 'put elevation')
      st = nc_close(ncid) ; call nc_check(st, 'close static')
   end subroutine write_static

   subroutine write_month(v, m)
      integer(ik), intent(in) :: v, m
      character(len=:), allocatable :: path
      integer(c_int)    :: st, ncid, dt, dla, dlo, vt, vx, dims1(1), dims3(3)
      integer(c_size_t) :: start1(1), count1(1), start3(3), count3(3)
      real(c_double), allocatable :: tsec(:), dat(:,:,:)
      real(c_double) :: nan
      integer(ik) :: nt, g0, h, r, c
      nan = ieee_value(0.0_c_double, ieee_quiet_nan)
      nt = era5land_month_hours(2021_ik, m)
      g0 = hour_index(meds_time_t(2021_ik, m, 1_ik, 1_ik))
      path = trim(era5land_path(era5land_default_template(), DIR, ERA_VAR_NAME(v), 2021_ik, m))
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

end program test_met_era5land
