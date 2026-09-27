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
   use meds_test_assert,     only : check, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_time,            only : meds_time_t
   use meds_forcing_config,  only : forcing_config_t, MET_BACKEND_ERA5LAND, METAVG_END,          &
                                    SWPART_CLEARIDX, SWPART_PASSTHROUGH, CLAMP_ERROR
   use meds_forcing_types,   only : met_source_t, met_cursor_t, met_forcing_t, met_cells_t, met_month_t
   use meds_forcing_kernels, only : dewpoint_to_specific_humidity
   use meds_lapse_rate,      only : wind_log_profile
   use meds_met_driver,      only : met_open, met_cursor_init, met_advance, met_instant, met_close, &
                                   met_prefetch, &
                                    MET_OK, MET_ERR_ARCHIVE, MET_ERR_ATTR_MISMATCH
   use meds_era5land_reader, only : era5land_path, era5land_default_template, era5land_select_site, &
                                    era5land_select_box, era5land_load_month, ERA_OK, ERA_ERR_NAN,  &
                                    ERA_ERR_NO_CELL, ERA_TAIR, ERA_TDEW, ERA_PSURF, ERA_U10, ERA_V10
   use meds_test_era5land_archive, only : T0, field, hour_index, write_archive
   implicit none

   character(len=*), parameter :: DIR = 'era5land_test_tmp'

   call write_archive(DIR)

   call test_paths()
   call test_site_selection()
   call test_box_selection()
   call test_nan_rejected()
   call test_driver_months()
   call test_driver_recycle()
   call test_driver_rejections()
   call test_region_cursors()

   call test_report('test_met_era5land')

contains

   !----- What the driver does around a step: prefetch the day's forcing (the step's I/O), then   !
   !      advance and sample inside it with no file access. Returns the instantaneous forcing.  ----!
   function step_sample(src, cur, t) result(met)
      type(met_source_t), intent(inout) :: src
      type(met_cursor_t), intent(inout) :: cur
      type(meds_time_t),  intent(in)    :: t
      type(met_forcing_t) :: met
      integer(ik) :: loads
      call met_prefetch(src, meds_time_t(t%year, t%month, t%day))
      loads = src%n_loads
      call met_advance(src, cur, t) ; met = met_instant(src, cur, t)
      if (src%n_loads /= loads) error stop 'test_met_era5land: the reader loaded a month inside a step'
   end function step_sample

   !----- The site's cursor into an opened source: cell 1, at the config's location. ---------------!
   subroutine site_cursor(src, cur, fc)
      type(met_source_t),     intent(in)  :: src
      type(met_cursor_t),     intent(out) :: cur
      type(forcing_config_t), intent(in)  :: fc
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%utc_offset_h,  &
                           fc%elevation_m)
   end subroutine site_cursor

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
      type(met_cells_t) :: dom
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
      type(met_cells_t) :: dom
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
      !----- A region's selection rule: column 1 is 30% land, so land_fraction_min = 0.5 drops     !
      !      its two cells and keeps the order of the rest. -----------------------------------------!
      call era5land_select_box(DIR//'/ED_ERA5land_static.nc', across, dom, st, land_fraction_min=0.5_wp)
      call check_true('land_fraction_min drops the mostly-water cells',                            &
                      st == ERA_OK .and. dom%ncell == 3_ik .and.                                   &
                      all(dom%row == [1_ik, 1_ik, 2_ik]) .and. all(dom%col == [17_ik, 0_ik, 17_ik]))
      call check_true('the box keeps each cell''s land fraction', all(dom%land_fraction == 1.0_wp))

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
      type(met_cells_t) :: dom
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
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      real(wp)    :: u1, u2, v1, v2, factor
      print '(a)', '-- era5land 5: driver across a month boundary --'
      fc = archive_config()
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call check_true('opens (MET_OK)', st == MET_OK)
      call site_cursor(src, cur, fc)
      call check_true('reads January and February: 1416 hourly records', src%nrec == 1416_ik)
      call check('grid_elevation from the static file', cur%grid_elevation_m, 103.0_wp, 1.0e-6_wp)
      call check_true('the archive supplies the wind vector', src%has_wind_vector)

      !----- At a stamp every weight is 0, so the record comes through untouched. ---------------!
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      met = step_sample(src, cur, t)
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
      met = step_sample(src, cur, t)
      call check('31 Jan 23:30 interpolates inside the January file', met%tair_k,                 &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)

      !----- The seam: 1 Feb 00:00 is January's last record, 01:00 February's first. -----------!
      t = meds_time_t(2021_ik, 2_ik, 1_ik, 0_ik, 30_ik) ; g = hour_index(meds_time_t(2021_ik, 2_ik, 1_ik))
      met = step_sample(src, cur, t)
      call check('1 Feb 00:30 brackets January''s last and February''s first record', met%tair_k, &
                 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, g) + field(ERA_TAIR, 1_ik, 3_ik, g + 1_ik)), 1.0e-9_wp)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; u2 = field(ERA_U10, 1_ik, 3_ik, g + 1_ik)
      v1 = field(ERA_V10, 1_ik, 3_ik, g) ; v2 = field(ERA_V10, 1_ik, 3_ik, g + 1_ik)
      call check('the vector interpolates linearly (u)', met%wind_u, 0.5_wp * (u1 + u2), 1.0e-12_wp)
      call check('the vector interpolates linearly (v)', met%wind_v, 0.5_wp * (v1 + v2), 1.0e-12_wp)
      call check_true('the energy-form speed is at least the vector''s length',                    &
                      met%wind >= sqrt(met%wind_u**2 + met%wind_v**2) - 1.0e-12_wp)
      !----- R1: January loaded at open, February at its first day's prefetch; the seam's 00:00      !
      !      record came from January's buffer, not from a second read. -----------------------------!
      call check_true('each month read once: January, then February', src%n_loads == 2_ik,          &
                      real(src%n_loads, wp))
      call met_close(src)

      !----- The height correction scales both components by the speed's factor. ---------------!
      fc%apply_wind_profile = .true. ; fc%wind_meas_height = 10.0_wp ; fc%reference_height = 40.0_wp
      fc%wind_roughness_z0 = 0.1_wp
      factor = log(400.0_wp) / log(100.0_wp)
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),               &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call site_cursor(src, cur, fc)
      t = meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik) ; g = hour_index(t)
      met = step_sample(src, cur, t)
      u1 = field(ERA_U10, 1_ik, 3_ik, g) ; v1 = field(ERA_V10, 1_ik, 3_ik, g)
      call check('height correction on u', met%wind_u, u1 * factor, 1.0e-12_wp)
      call check('height correction on v', met%wind_v, v1 * factor, 1.0e-12_wp)
      call check('direction preserved', atan2(met%wind_v, met%wind_u), atan2(v1, u1), 1.0e-12_wp)
      call met_close(src)
   end subroutine test_driver_months

   !=======================================================================================!
   !  6. Recycling across the months of a declared window, including the year-end seam, which   !
   !     brackets December's last record (2022-01-01 00:00) and January's first.                 !
   !=======================================================================================!
   subroutine test_driver_recycle()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t
      integer(ik) :: st, g
      print '(a)', '-- era5land 6: recycling across months --'
      fc = archive_config()
      fc%recycle = .true.
      fc%recycle_start = T0 ; fc%recycle_end = meds_time_t(2022_ik, 1_ik, 1_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('a whole-year window opens without a run period', st == MET_OK)
      call site_cursor(src, cur, fc)
      call check_true('the window is 8760 records', src%nrec == 8760_ik)
      t = meds_time_t(2035_ik, 6_ik, 15_ik, 12_ik) ; g = hour_index(meds_time_t(2021_ik, 6_ik, 15_ik, 12_ik))
      met = step_sample(src, cur, t)
      call check('2035-06-15 12:00 reads 2021-06-15 12:00', met%tair_k, field(ERA_TAIR, 1_ik, 3_ik, g), 1.0e-9_wp)
      t = meds_time_t(2036_ik, 1_ik, 1_ik, 0_ik, 30_ik)
      met = step_sample(src, cur, t)
      call check('the year-end seam brackets December''s last and January''s first record',       &
                 met%tair_k, 0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, 8759_ik) + field(ERA_TAIR, 1_ik, 3_ik, 0_ik)), &
                 1.0e-9_wp)
      !----- Walk the model from 2036-12-30 across the wrap day by day: the wrap step needs December's  !
      !      last record and January, and the reader reads January once more, not December again. ---!
      block
         integer(ik) :: loads
         t = meds_time_t(2036_ik, 12_ik, 30_ik, 12_ik) ; met = step_sample(src, cur, t)
         t = meds_time_t(2036_ik, 12_ik, 31_ik, 23_ik, 30_ik) ; met = step_sample(src, cur, t)
         loads = src%n_loads
         t = meds_time_t(2037_ik, 1_ik, 1_ik, 0_ik, 30_ik) ; met = step_sample(src, cur, t)
         call check('across the wrap: December''s last and January''s first record', met%tair_k,   &
                    0.5_wp * (field(ERA_TAIR, 1_ik, 3_ik, 8759_ik) + field(ERA_TAIR, 1_ik, 3_ik, 0_ik)), 1.0e-9_wp)
         call check_true('the wrap reads January once, December''s record comes from the buffer',  &
                         src%n_loads == loads + 1_ik, real(src%n_loads - loads, wp))
      end block
      call met_close(src)
   end subroutine test_driver_recycle

   !=======================================================================================!
   !  7. Rejected at open: a month the archive lacks, total shortwave read as components, and  !
   !     a site with no valid cell in reach.                                                     !
   !=======================================================================================!
   subroutine test_driver_rejections()
      type(met_source_t)     :: src
      type(forcing_config_t) :: fc
      integer(ik) :: st
      print '(a)', '-- era5land 7: rejections at open --'
      fc = archive_config()
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 11_ik, 1_ik),               &
                    run_end=meds_time_t(2022_ik, 2_ik, 1_ik))
      call check_true('a run past the archive''s last month is rejected', st == MET_ERR_ARCHIVE)
      fc%sw_partition = SWPART_PASSTHROUGH
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('passthrough against total shortwave is rejected', st == MET_ERR_ATTR_MISMATCH)
      fc = archive_config()
      fc%latitude_deg = 45.0_wp ; fc%longitude_deg = -80.0_wp
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik, 3_ik, 1_ik),                &
                    run_end=meds_time_t(2021_ik, 4_ik, 1_ik))
      call check_true('a site with no valid cell within max_distance_km is rejected', st == MET_ERR_ARCHIVE)
   end subroutine test_driver_rejections

   !=======================================================================================!
   !  8. One source, many cursors (MEDS_POLYGON_RUNTIME_PLAN.md R2): a region opens the archive  !
   !     once for its cells, and each polygon walks it with its own cursor. Each cursor must read !
   !     what a site run at its cell reads, bit for bit, and one prefetch must serve them all.    !
   !=======================================================================================!
   subroutine test_region_cursors()
      integer(ik), parameter :: PICK(2) = [2_ik, 5_ik], PR(2) = [1_ik, 2_ik], PC(2) = [0_ik, 1_ik]
      type(met_source_t)     :: region, site(2)
      type(met_cursor_t)     :: cur(2), cur_site(2)
      type(met_cells_t)      :: cells
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met, met_site
      type(meds_time_t)      :: t, times(3)
      integer(ik) :: st, k, j, loads
      logical     :: same
      print '(a)', '-- era5land 8: one region source, one cursor per polygon --'
      !----- The box across 180 from test 3: (1,17), (1,0), (1,1), (2,17), (2,1). Polygons at the   !
      !      second and fifth cells, (1,0) and (2,1); a site run at each cell for comparison. ------!
      call era5land_select_box(DIR//'/ED_ERA5land_static.nc', [40.0_wp, 150.0_wp, 20.0_wp, -150.0_wp], &
                               cells, st)
      fc = archive_config()
      call met_open(region, fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),            &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik), cells=cells)
      call check_true('a region opens on its cells (MET_OK)', st == MET_OK .and. region%cells%ncell == 5_ik)
      do k = 1_ik, 2_ik
         call met_cursor_init(region, cur(k), PICK(k), cells%lat(PICK(k)), cells%lon(PICK(k)),    &
                              fc%utc_offset_h, fc%elevation_m)
         fc = archive_config()
         fc%latitude_deg = cells%lat(PICK(k)) ; fc%longitude_deg = cells%lon(PICK(k))
         call met_open(site(k), fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),        &
                       run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
         call site_cursor(site(k), cur_site(k), fc)
      end do
      call check('the first cursor takes its cell''s orography', cur(1)%grid_elevation_m, 100.0_wp, 1.0e-6_wp)
      call check('the second cursor takes its cell''s orography', cur(2)%grid_elevation_m, 201.0_wp, 1.0e-6_wp)

      !----- A stamp, a mid-month interpolation and the February seam. ----------------------------!
      times = [meds_time_t(2021_ik, 1_ik, 20_ik, 6_ik), meds_time_t(2021_ik, 1_ik, 31_ik, 23_ik, 30_ik), &
               meds_time_t(2021_ik, 2_ik, 1_ik, 0_ik, 30_ik)]
      same = .true.
      do j = 1_ik, size(times)
         t = times(j)
         call met_prefetch(region, meds_time_t(t%year, t%month, t%day))
         loads = region%n_loads
         do k = 1_ik, 2_ik
            call met_advance(region, cur(k), t) ; met = met_instant(region, cur(k), t)
            met_site = step_sample(site(k), cur_site(k), t)
            same = same .and. met%tair_k == met_site%tair_k .and. met%qair == met_site%qair       &
                 .and. met%psurf_pa == met_site%psurf_pa .and. met%rainf == met_site%rainf      &
                 .and. met%snowfall == met_site%snowfall .and. met%wind == met_site%wind        &
                 .and. met%wind_u == met_site%wind_u .and. met%wind_v == met_site%wind_v        &
                 .and. met%lwdown == met_site%lwdown .and. met%par_beam == met_site%par_beam    &
                 .and. met%par_diffuse == met_site%par_diffuse                                  &
                 .and. met%nir_beam == met_site%nir_beam                                        &
                 .and. met%nir_diffuse == met_site%nir_diffuse .and. met%cosz == met_site%cosz
            if (j == 1_ik) call check('at a stamp the cursor reads its own cell''s Tair', met%tair_k, &
                                      field(ERA_TAIR, PR(k), PC(k), hour_index(t)), 1.0e-9_wp)
         end do
         call check_true('one prefetch serves every cursor (no read inside the step)',           &
                         region%n_loads == loads)
      end do
      call check_true('each polygon reads what a site run at its cell reads, bit for bit', same)
      call check_true('the region read each month once', region%n_loads == 2_ik)
      call met_close(region) ; call met_close(site(1))

      !----- A source closed and opened again starts clean (the C API reuses a run's source): on the  !
      !      second cell, a fresh cursor's first bracket at 1 Feb 00:30 needs January's last record,  !
      !      which would be the first cell's if met_close left the carried record behind.  -------------!
      fc = archive_config()
      fc%latitude_deg = cells%lat(PICK(2)) ; fc%longitude_deg = cells%lon(PICK(2))
      call met_open(site(1), fc, stat=st, run_start=meds_time_t(2021_ik, 1_ik, 15_ik),             &
                    run_end=meds_time_t(2021_ik, 3_ik, 1_ik))
      call site_cursor(site(1), cur_site(1), fc)
      t = times(3)
      met = step_sample(site(1), cur_site(1), t)
      call check('a reopened source reads its new cell''s carried record', met%tair_k,               &
                 0.5_wp * (field(ERA_TAIR, PR(2), PC(2), hour_index(meds_time_t(2021_ik, 2_ik, 1_ik)))  &
                 + field(ERA_TAIR, PR(2), PC(2), hour_index(meds_time_t(2021_ik, 2_ik, 1_ik)) + 1_ik)), 1.0e-9_wp)
      call met_close(site(1)) ; call met_close(site(2))
   end subroutine test_region_cursors

end program test_met_era5land
