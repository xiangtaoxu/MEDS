! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_region -- a region run is its site runs (MEDS_POLYGON_RUNTIME_PLAN.md §10.3, R2).          !
!                                                                                          !
! On the synthetic ED_ERA5land archive, a region over a box that crosses 180 degrees, holds a      !
! no-data cell and a mostly-water column, runs its three polygons for twelve days across a month   !
! boundary. Then each polygon is run again as a site at its cell (centre, orography, UTC). The     !
! region's files must hold, in every variable's `polygon = p` slice, exactly the site run's        !
! series, bit for bit; the detail polygon's own files must equal its site run's files; and the     !
! polygon axis must carry the cells' ids and coordinates. The region is run again on four threads  !
! and must give the same files. Last, meds_main (the second argument)                               !
! must refuse each config that breaks a region-mode rule (§9), with that rule's message.            !
!                                                                                          !
! The configs are DERIVED from test/region/meds_config_site.toml (the test runs from that          !
! directory): an override block is prepended -- meds_toml returns the FIRST match of a              !
! key -- and the keys the archive or region mode reject are dropped. The monthly tier is chunked by  !
! month: a site file caps a tier with cohort variables at a month, and the region's (which has none)  !
! would otherwise keep the template's year, so the same records would sit in differently named       !
! files. Everything the test writes goes under the work directory given as the first argument.       !
!==========================================================================================!
program test_region
   use iso_c_binding,              only : c_int, c_size_t, c_double
   use meds_kinds,                 only : wp, ik
   use meds_test_assert,           only : check_true, test_report
   use meds_test_era5land_archive, only : write_archive
   use meds_netcdf_c
   use meds_output_types,          only : output_registry_t, DIM_SCALAR, DIM_COHORT, DIM_PATCH,   &
                                          DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, XTYPE_INT,  &
                                          MISSING_INT
   use meds_driver,                only : meds_run_t, driver_open, driver_step, driver_finalize,  &
                                          driver_free, driver_done, DRIVER_OK, DRIVER_ERR_SOILC
   use meds_region,                only : meds_region_t, region_open, region_step_month,          &
                                          region_done, region_finalize, region_free
   implicit none

   character(len=*), parameter :: TEMPLATE = 'meds_config_site.toml'
   character(len=1), parameter :: LETTER(4) = ['F', 'D', 'M', 'Y']
   !----- The box [N, W, S, E] = [40, 150, 20, -150] holds rows 1-2 (35 N, 25 N) and columns 17, 0, 1  !
   !      (160 E, 180, 160 W). (2,0) has no data and column 1 is 30% land, so land_fraction_min = 0.5 !
   !      leaves (1,17), (1,0), (2,17): ids row * 18 + col, cells at the orography 100 row + col.  -----!
   integer(ik), parameter :: NP = 3_ik
   integer(ik), parameter :: IDS(NP) = [35_ik, 18_ik, 53_ik]
   real(wp),    parameter :: LAT(NP) = [35.0_wp, 35.0_wp, 25.0_wp], LON(NP) = [160.0_wp, -180.0_wp, 160.0_wp]
   real(wp),    parameter :: ELEV(NP) = [117.0_wp, 100.0_wp, 217.0_wp]
   integer(ik), parameter :: DETAIL = 18_ik

   character(len=512)  :: work, meds_main
   type(meds_region_t) :: reg
   type(meds_run_t)    :: run
   type(output_registry_t) :: reg_region, reg_site
   character(len=64), allocatable :: files(:)
   integer(ik) :: p, st, f, nf, nbad
   logical     :: ok
   character(len=24) :: idstr

   call get_command_argument(1, work)
   call get_command_argument(2, meds_main)
   if (len_trim(work) == 0 .or. len_trim(meds_main) == 0)                                           &
      error stop 'test_region: pass the work directory and the meds_main executable'
   call execute_command_line('rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'/out"')
   call write_archive(trim(work)//'/archive')

   !----- The region. ----------------------------------------------------------------------------!
   call derive(trim(work)//'/region.toml', region_block(), region=.true.)
   call region_open(trim(work)//'/region.toml', reg, ok, verbose=.false.)
   call check_true('the region opens', ok)
   if (.not. ok) error stop 'test_region: the region did not open'
   call check_true('three polygons: a no-data cell and a mostly-water column are not simulated',     &
                   size(reg%poly) == NP)
   call check_true('polygon ids are the cells'' row-major indices, in row-major order',              &
                   all(reg%poly(:)%id == IDS))
   do while (.not. region_done(reg))
      call region_step_month(reg, st)
      if (st /= DRIVER_OK) error stop 'test_region: a region month failed'
   end do
   call region_finalize(reg, st)
   call check_true('every polygon conserves area', st == DRIVER_OK)
   reg_region = reg%out_files%reg
   call region_free(reg)

   !----- The same cells as sites. ---------------------------------------------------------------!
   do p = 1_ik, NP
      write(idstr,'(i0)') IDS(p)
      call derive(trim(work)//'/site'//trim(idstr)//'.toml', site_block(p), region=.false.)
      call driver_open(trim(work)//'/site'//trim(idstr)//'.toml', run, ok, verbose=.false.)
      if (.not. ok) error stop 'test_region: a site run did not open'
      do while (.not. driver_done(run))
         call driver_step(run, st)
         if (st /= DRIVER_OK) error stop 'test_region: a site step failed'
      end do
      call driver_finalize(run, st)
      if (p == 1_ik) reg_site = run%out_files%reg
      call driver_free(run)
   end do

   !----- The polygon axis. ----------------------------------------------------------------------!
   call list_files(trim(work)//'/out', files, nf)
   call check_polygon_axis(trim(work)//'/out/region-D-202101.nc')

   !----- Every region file against the site files: each polygon's slice, every variable. --------!
   nbad = 0_ik
   do f = 1_ik, nf
      if (index(files(f), 'region-') /= 1 .or. index(files(f), 'region-p') == 1) cycle
      do p = 1_ik, NP
         write(idstr,'(i0)') IDS(p)
         nbad = nbad + compare_region_file(trim(work)//'/out/'//trim(files(f)), p,                &
                                           trim(work)//'/out/site'//trim(idstr)//files(f)(7:))
      end do
   end do
   call check_true('each region variable''s polygon slice equals its site run, bit for bit', nbad == 0_ik)
   call check_true('region files: two months each of the daily and monthly tiers, no fast tier',       &
                   count_prefix('region-D-') == 2_ik .and. count_prefix('region-M-') == 2_ik .and.     &
                   count_prefix('region-F-') == 0_ik)

   !----- The detail polygon's own files against its site run's files. -----------------------------!
   nbad = 0_ik
   write(idstr,'(i0)') DETAIL
   do f = 1_ik, nf
      if (index(files(f), 'site'//trim(idstr)//'-') /= 1) cycle
      nbad = nbad + compare_site_file(trim(work)//'/out/region-p'//trim(idstr)//                   &
                                      files(f)(len_trim('site'//idstr) + 1:),                        &
                                      trim(work)//'/out/'//trim(files(f)))
   end do
   call check_true('the detail polygon writes its site run''s files, fast tier included',             &
                   count_prefix('region-p'//trim(idstr)//'-') == count_prefix('site'//trim(idstr)//'-') &
                   .and. count_prefix('region-p'//trim(idstr)//'-F-') == 12_ik)
   call check_true('the detail polygon''s files equal its site run''s, every variable', nbad == 0_ik)

   !----- The same region on four threads (#183 R3): the polygons run side by side, and every file   !
   !      must still equal the site runs bit for bit -- region files slice by slice, the detail      !
   !      polygon's file by file. In a build without OpenMP this is the serial run again. -------------!
   call derive(trim(work)//'/regt.toml', '[run]'//nl()//'n_threads = 4'//nl()//'[output]'//nl()//     &
               'prefix = "regt"'//nl()//region_block(), region=.true.)
   call region_open(trim(work)//'/regt.toml', reg, ok, verbose=.false.)
   if (.not. ok) error stop 'test_region: the threaded region did not open'
   do while (.not. region_done(reg))
      call region_step_month(reg, st)
      if (st /= DRIVER_OK) error stop 'test_region: a threaded region month failed'
   end do
   call region_finalize(reg, st)
   call region_free(reg)
   nbad = 0_ik
   do f = 1_ik, nf
      if (index(files(f), 'region-') /= 1) cycle
      if (index(files(f), 'region-p') == 1) then
         write(idstr,'(i0)') DETAIL
         nbad = nbad + compare_site_file(trim(work)//'/out/regt'//files(f)(7:),                      &
                                         trim(work)//'/out/site'//trim(idstr)//                       &
                                         files(f)(len_trim('region-p'//idstr) + 1:))
      else
         do p = 1_ik, NP
            write(idstr,'(i0)') IDS(p)
            nbad = nbad + compare_region_file(trim(work)//'/out/regt'//files(f)(7:), p,             &
                                              trim(work)//'/out/site'//trim(idstr)//files(f)(7:))
         end do
      end if
   end do
   call check_true('on four threads every region and detail file still equals the site runs', nbad == 0_ik)

   !----- A reused run owes no restructuring from the run before it. A run ending on the 1st leaves  !
   !      its boundary's restructuring pending; driver_open must clear that, or the next stand is     !
   !      restructured at its first step. -------------------------------------------------------------!
   call derive(trim(work)//'/reuse.toml', '[run]'//nl()//'end_time = "2021-02-01"'//nl()//            &
               '[output]'//nl()//'enabled = false'//nl()//site_block(1_ik), region=.false.)
   call driver_open(trim(work)//'/reuse.toml', run, ok, verbose=.false.)
   if (.not. ok) error stop 'test_region: the reuse run did not open'
   do while (.not. driver_done(run))
      call driver_step(run, st)
   end do
   call check_true('a run ending on the 1st leaves its boundary''s restructuring owed',              &
                   run%poly%restructure_pending)
   call driver_finalize(run, st) ; call driver_free(run)
   call driver_open(trim(work)//'/reuse.toml', run, ok, verbose=.false.)
   call check_true('a reused run opens owing no restructuring', ok .and. .not. run%poly%restructure_pending)
   !----- ...and its canopy air starts at the bare-ground floor, not the type's default depth (#306). -!
   call check_true('a bare-ground stand opens at the canopy-air floor',                               &
                   all(run%poly%site%patch%cas(1:run%poly%site%patch%n)%can_depth                     &
                       == run%cfg%aero%min_canopy_depth))
   call driver_finalize(run, st) ; call driver_free(run)

   !----- A polygon that fails in mid-month. Polygon 2's soil carbon is made impossible after the    !
   !      first month, so its first February step fails. The others finish the month: polygons 1     !
   !      and 3 close February's five days and polygon 2 its first. The month's I/O phase writes      !
   !      what closed, with the fill value where polygon 2 closed nothing, and the region moves on. ---!
   call derive(trim(work)//'/fail.toml', '[output]'//nl()//'prefix = "fail"'//nl()//                  &
               '[region]'//nl()//'detail_polygons = [35]'//nl()//region_block(), region=.true.)
   call region_open(trim(work)//'/fail.toml', reg, ok, verbose=.false.)
   if (.not. ok) error stop 'test_region: the failing region did not open'
   call region_step_month(reg, st)
   call check_true('the month before the failure steps', st == DRIVER_OK)
   reg%poly(2)%site%patch%soil_carbon(1)%slow_carbon = -1.0_wp
   call region_step_month(reg, st)
   call check_true('the month a polygon fails in returns its status', st == DRIVER_ERR_SOILC)
   call check_true('the failed month is consumed: the region moves on', region_done(reg))
   call region_finalize(reg, st)
   call region_free(reg)
   call check_failed_month(trim(work)//'/out/fail-D-202102.nc')

   !----- Region-mode rules (MEDS_POLYGON_RUNTIME_PLAN.md §9): each bad config is refused with its  !
   !      own message. The bad keys come first, so they win over the good region block. -------------!
   call refused('the [site] location is refused in a region', region_block(), .false., .false.,     &
                'do not apply to [run].mode = "region"')
   call refused('forcing.max_distance_km is refused in a region',                                  &
                '[forcing]'//nl()//'max_distance_km = 50.0'//nl()//region_block(), .true., .false.,   &
                'do not apply to [run].mode = "region"')
   call refused('a region needs the archive', '[forcing]'//nl()//'format = "ED_default"'//nl()//region_block(), &
                .true., .true., 'needs forcing.format = "ED_ERA5land"')
   call refused('a region starts from bare ground', '[init]'//nl()//'init_mode = 2'//nl()//region_block(), &
                .true., .false., 'starts from bare ground')
   call refused('a region writes no checkpoints', '[state]'//nl()//'write_state = true'//nl()//region_block(), &
                .true., .false., 'writes no checkpoints')
   call refused('a fast-tier window may not straddle two slow steps',                               &
                '[output]'//nl()//'fast_interval_steps = 5'//nl()//site_block(1_ik), .false., .false., &
                'output.fast_interval_steps must divide')
   call refused('output.strict_caps is refused: the caps always stop the run',                     &
                '[output]'//nl()//'strict_caps = true'//nl()//site_block(1_ik), .false., .false.,     &
                'is retired: a run always stops at the step')
   call refused('a [region] block needs region mode', '[run]'//nl()//'mode = "site"'//nl()//          &
                site_block(1_ik)//'[region]'//nl()//'land_fraction_min = 0.5'//nl(), .false., .false.,  &
                'a [region] block needs [run].mode = "region"')
   call refused('the box must have S < N', '[region]'//nl()//'box_nwse = [20.0, 150.0, 40.0, -150.0]'//nl()// &
                region_block(), .true., .false., 'region.box_nwse needs')
   call refused('a detail polygon must be one of the region''s',                                    &
                '[region]'//nl()//'detail_polygons = [0]'//nl()//region_block(), .true., .false.,       &
                'is not a polygon of this region')
   call refused('the mode is "site" or "region"', '[run]'//nl()//'mode = "globe"'//nl()//region_block(), &
                .true., .false., '[run].mode must be "site" or "region"')
   call refused('the [io] block is refused, naming [state]', '[io]'//nl()//'write_state = false'//nl()// &
                site_block(1_ik), .false., .false., 'the [io] block is now [state]')
   !----- The root profile is a plant trait: the old [soil_column] key is refused, naming [hydraulics]. -!
   call refused('soil_column.root_beta is refused, naming [hydraulics]',                            &
                '[soil_column]'//nl()//'root_beta = 2.0'//nl()//site_block(1_ik), .false., .false.,  &
                'is retired: the root profile is a plant trait: set [hydraulics].root_beta')
   call refused('hydraulics.conductance is whole_plant or segment',                                 &
                '[hydraulics]'//nl()//'conductance = "xylem"'//nl()//site_block(1_ik), .false., .false., &
                'hydraulics.conductance must be')
   call refused('leaf_physiology.low_water_potential_control is linear_decline',                    &
                '[leaf_physiology]'//nl()//'low_water_potential_control = "clamp"'//nl()//          &
                site_block(1_ik), .false., .false., 'low_water_potential_control must be')
   !----- The forcing-format rules (MEDS_FLUX_TOWER_FORCING_PLAN.md §5): the old format names stop  !
   !      naming their replacements, every clock is UTC, and the terrain-lapse keys need the lapse. --!
   call refused('format = "netcdf" is refused, naming ED_default',                                  &
                '[forcing]'//nl()//'format = "netcdf"'//nl()//site_block(1_ik), .false., .false.,     &
                'is now "ED_default"')
   call refused('format = "era5land" is refused, naming ED_ERA5land',                               &
                '[forcing]'//nl()//'format = "era5land"'//nl()//site_block(1_ik), .false., .false.,   &
                'is now "ED_ERA5land"')
   call refused('site.utc_offset is refused: every clock is UTC',                                   &
                '[site]'//nl()//'utc_offset = -5.0'//nl()//site_block(1_ik), .false., .false.,         &
                'is retired: every forcing file is in UTC')
   call refused('site.apply_solar_longitude is refused',                                            &
                '[site]'//nl()//'apply_solar_longitude = true'//nl()//site_block(1_ik), .false., .false., &
                'is retired: every forcing file is in UTC')
   call refused('the lapse keys are refused with the lapse off',                                    &
                '[site]'//nl()//'apply_elevation_lapse = false'//nl()//site_block(1_ik), .false., .false., &
                'apply only with site.apply_elevation_lapse = true')
   call refused('a region''s recycle window starts at 00:00 or 01:00 on the 1st',                   &
                '[forcing]'//nl()//'recycle = true'//nl()//'recycle_start = "2021-01-01 06:00:00"'//nl()// &
                'recycle_end = "2022-01-01 06:00:00"'//nl()//region_block(), .true., .false.,          &
                'needs forcing.recycle_start at 00:00 or')

   call test_report('test_region')

contains

   !----- The override blocks. -------------------------------------------------------------------!
   function common_block() result(b)
      character(len=:), allocatable :: b
      b = '[run]'//nl()//'start_time = "2021-01-25"'//nl()//'end_time = "2021-02-06"'//nl()//        &
          '[fast]'//nl()//'fast_probe = false'//nl()//                                              &
          '[init]'//nl()//'init_mode = 0'//nl()//                                                   &
          '[state]'//nl()//'write_state = false'//nl()//                                            &
          '[forcing]'//nl()//'format = "ED_ERA5land"'//nl()//'data_path = "'//trim(work)//'/archive"'//nl()// &
          'recycle = false'//nl()//                                                                 &
          '[output]'//nl()//'dir = "'//trim(work)//'/out"'//nl()//                                  &
          '[output.monthly]'//nl()//'enabled = true'//nl()//'file_chunk = "month"'//nl()
   end function common_block

   function region_block() result(b)
      character(len=:), allocatable :: b
      b = '[run]'//nl()//'mode = "region"'//nl()//'[output]'//nl()//'prefix = "region"'//nl()//     &
          '[region]'//nl()//'box_nwse = [40.0, 150.0, 20.0, -150.0]'//nl()//                        &
          'land_fraction_min = 0.5'//nl()//'detail_polygons = [18]'//nl()//common_block()
   end function region_block

   function site_block(p) result(b)
      integer(ik), intent(in) :: p
      character(len=:), allocatable :: b
      character(len=24)  :: id
      write(id,'(i0)') IDS(p)
      b = '[site]'//nl()//'latitude = '//fmt(LAT(p))//nl()//'longitude = '//fmt(LON(p))//nl()//     &
          'elevation = '//fmt(ELEV(p))//nl()//                                                     &
          '[forcing]'//nl()//'max_distance_km = 50.0'//nl()//                                      &
          '[output]'//nl()//'prefix = "site'//trim(id)//'"'//nl()//common_block()
   end function site_block

   function fmt(x) result(s)
      real(wp), intent(in) :: x
      character(len=:), allocatable :: s
      character(len=32) :: buf
      write(buf,'(f0.4)') x
      s = trim(adjustl(buf))
      if (s(1:1) == '.') s = '0'//s
      if (s(1:2) == '-.') s = '-0'//s(2:)
   end function fmt

   pure function nl() result(c)
      character(len=1) :: c
      c = new_line('a')
   end function nl

   !----- meds_main must fail on this config and say why. --------------------------------------------!
   subroutine refused(what, overrides, region, keep_file_keys, message)
      character(len=*), intent(in) :: what, overrides, message
      logical,          intent(in) :: region, keep_file_keys
      character(len=:), allocatable :: cfg, log
      character(len=1024) :: line
      integer :: es, u, ios
      logical :: said
      cfg = trim(work)//'/refused.toml' ; log = trim(work)//'/refused.log'
      call derive(cfg, overrides, region, keep_file_keys)
      call execute_command_line('"'//trim(meds_main)//'" "'//cfg//'" > "'//log//'" 2>&1', exitstat=es)
      said = .false.
      open(newunit=u, file=log, status='old', action='read')
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, message) > 0) said = .true.
      end do
      close(u)
      call check_true(what, es /= 0 .and. said)
   end subroutine refused

   !----- Write the override block, then the template without the keys that do not apply. ---------!
   subroutine derive(path, overrides, region, keep_file_keys)
      character(len=*), intent(in) :: path, overrides
      logical,          intent(in) :: region
      logical, optional, intent(in) :: keep_file_keys
      character(len=1024) :: line, sec, key
      integer :: uin, uout, ios, i
      open(newunit=uout, file=path, status='replace', action='write')
      write(uout,'(a)') overrides
      open(newunit=uin, file=TEMPLATE, status='old', action='read')
      sec = ''
      do
         read(uin, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) == '[') then
            sec = line(2:index(line, ']') - 1)
         else
            i = index(line, '=')
            if (i > 1) then
               key = adjustl(line(1:i - 1))
               if (drop(trim(sec), trim(key), region, keep_file_keys)) cycle
            end if
         end if
         write(uout,'(a)') trim(line)
      end do
      close(uin) ; close(uout)
   end subroutine derive

   logical function drop(sec, key, region, keep_file_keys)
      character(len=*), intent(in) :: sec, key
      logical,          intent(in) :: region
      logical, optional, intent(in) :: keep_file_keys
      drop = (sec == 'forcing' .and. (key == 'path' .or. key == 'grid_index' .or. key == 'grid_match')) &
             .or. (sec == 'site' .and. key == 'grid_elevation')
      if (present(keep_file_keys)) then
         if (keep_file_keys) drop = .false.
      end if
      if (region) drop = drop .or. (sec == 'site' .and. (key == 'latitude' .or. key == 'longitude' .or. &
                                                          key == 'elevation'))
   end function drop

   !----- The files in a directory, sorted. -------------------------------------------------------!
   subroutine list_files(dir, names, n)
      character(len=*),               intent(in)  :: dir
      character(len=64), allocatable, intent(out) :: names(:)
      integer(ik),                    intent(out) :: n
      character(len=64) :: buf(512)
      integer :: u, ios
      call execute_command_line('ls "'//dir//'" > "'//trim(work)//'/files.txt"')
      open(newunit=u, file=trim(work)//'/files.txt', status='old', action='read')
      n = 0_ik
      do
         read(u, '(a)', iostat=ios) buf(n + 1_ik)
         if (ios /= 0) exit
         n = n + 1_ik
      end do
      close(u)
      names = buf(1:n)
   end subroutine list_files

   integer(ik) function count_prefix(prefix) result(c)
      character(len=*), intent(in) :: prefix
      integer(ik) :: i
      c = 0_ik
      do i = 1_ik, nf
         if (index(files(i), prefix) == 1) c = c + 1_ik
      end do
   end function count_prefix

   !----- The tier of a file, from its name: <prefix>-<letter>-<stamp>.nc. ------------------------!
   integer(ik) function tier_of(name) result(t)
      character(len=*), intent(in) :: name
      integer :: i
      i = index(name, '-', back=.true.)
      if (i < 2) then ; t = 0_ik ; return ; end if
      do t = 1_ik, 4_ik
         if (name(i - 1:i - 1) == LETTER(t) .and. name(i - 2:i - 2) == '-') return
      end do
      t = 0_ik
      do t = 1_ik, 4_ik                              ! a stampless run file: <prefix>-<letter>.nc
         if (name(i + 1:i + 1) == LETTER(t) .and. name(i + 2:i + 4) == '.nc') return
      end do
      t = 0_ik
   end function tier_of

   integer(c_int) function open_nc(path) result(ncid)
      character(len=*), intent(in) :: path
      call nc_check(nc_open_f(path, NC_NOWRITE, ncid), 'open '//path)
   end function open_nc

   integer(ik) function dim_len(ncid, name) result(n)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name
      integer(c_size_t) :: len
      n = 0_ik
      if (nc_inq_dimlen_f(ncid, name, len) == NC_NOERR) n = int(len, ik)
   end function dim_len

   !----- Read a whole variable of the given C-order shape into a flat array. -----------------------!
   subroutine read_var(ncid, name, shape_c, x, found)
      integer(c_int),        intent(in)  :: ncid
      character(len=*),      intent(in)  :: name
      integer(ik),           intent(in)  :: shape_c(:)
      real(c_double), allocatable, intent(out) :: x(:)
      logical,               intent(out) :: found
      integer(c_int)    :: vid
      integer(c_size_t) :: start(size(shape_c)), count(size(shape_c))
      found = nc_inq_varid_f(ncid, name, vid) == NC_NOERR
      allocate(x(product(shape_c)))
      if (.not. found .or. size(x) == 0) return
      start = 0_c_size_t ; count = int(shape_c, c_size_t)
      call nc_check(nc_get_vara_double(ncid, vid, start, count, x), 'read '//name)
   end subroutine read_var

   !----- A variable's shape (C order, time first) from its registry axis and the file's dims. -----!
   function var_shape(ncid, dm, with_polygon) result(sh)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: dm
      logical,        intent(in) :: with_polygon
      integer(ik), allocatable :: sh(:)
      integer(ik) :: nt
      nt = dim_len(ncid, 'time')
      select case (dm)
      case (DIM_SCALAR)     ; sh = [integer(ik) ::]
      case (DIM_COHORT)     ; sh = [dim_len(ncid, 'cohort')]
      case (DIM_PATCH)      ; sh = [dim_len(ncid, 'patch')]
      case (DIM_SOIL)       ; sh = [dim_len(ncid, 'soil')]
      case (DIM_PFT)        ; sh = [dim_len(ncid, 'pft')]
      case (DIM_SIZE)       ; sh = [dim_len(ncid, 'dbh_class')]
      case (DIM_SOIL_PATCH) ; sh = [dim_len(ncid, 'patch'), dim_len(ncid, 'soil')]
      end select
      if (with_polygon) sh = [dim_len(ncid, 'polygon'), sh]
      sh = [nt, sh]
   end function var_shape

   !----- Region file vs one polygon's site file: the polygon's slice of every live variable. ------!
   integer(ik) function compare_region_file(rpath, p, spath) result(nbad)
      character(len=*), intent(in) :: rpath, spath
      integer(ik),      intent(in) :: p
      integer(c_int) :: rn, sn
      integer(ik)    :: t, j, k, nt, na, np_, it
      integer(ik), allocatable :: rs(:), ss(:)
      real(c_double), allocatable :: xr(:), xs(:)
      logical :: fr, fs
      nbad = 0_ik
      t = tier_of(rpath)
      rn = open_nc(rpath) ; sn = open_nc(spath)
      do j = 1_ik, reg_region%nidx(t)
         k = reg_region%idx_freq(j, t)
         rs = var_shape(rn, reg_region%var(k)%dim, .true.)
         ss = var_shape(sn, reg_region%var(k)%dim, .false.)
         call read_var(rn, trim(reg_region%var(k)%name), rs, xr, fr)
         call read_var(sn, trim(reg_region%var(k)%name), ss, xs, fs)
         if (.not. (fr .and. fs) .or. rs(1) /= ss(1)) then
            write(*,'(5a)') '  missing or mismatched: ', trim(reg_region%var(k)%name), ' in ', rpath, ''
            nbad = nbad + 1_ik ; cycle
         end if
         nt = rs(1) ; np_ = rs(2) ; na = product(rs(3:))      ! C order: axis fastest, then polygon, time
         do it = 0_ik, nt - 1_ik
            if (any(xr(it*np_*na + (p - 1_ik)*na + 1_ik : it*np_*na + p*na) /= xs(it*na + 1_ik : (it + 1_ik)*na))) then
               write(*,'(3a,i0,2a)') '  DIFFERS: ', trim(reg_region%var(k)%name), ' polygon ', IDS(p), ' in ', rpath
               nbad = nbad + 1_ik ; exit
            end if
         end do
      end do
      call nc_check(nc_close(rn), 'close') ; call nc_check(nc_close(sn), 'close')
   end function compare_region_file

   !----- Two single-site files of the same registry: every live variable of the tier. ------------!
   integer(ik) function compare_site_file(apath, bpath) result(nbad)
      character(len=*), intent(in) :: apath, bpath
      integer(c_int) :: an, bn
      integer(ik)    :: t, j, k
      integer(ik), allocatable :: as(:), bs(:)
      real(c_double), allocatable :: xa(:), xb(:)
      logical :: fa, fb, ex
      nbad = 0_ik
      inquire(file=apath, exist=ex)
      if (.not. ex) then
         write(*,'(2a)') '  missing detail file ', apath ; nbad = 1_ik ; return
      end if
      t = tier_of(bpath)
      an = open_nc(apath) ; bn = open_nc(bpath)
      do j = 1_ik, reg_site%nidx(t)
         k = reg_site%idx_freq(j, t)
         as = var_shape(an, reg_site%var(k)%dim, .false.)
         bs = var_shape(bn, reg_site%var(k)%dim, .false.)
         call read_var(an, trim(reg_site%var(k)%name), as, xa, fa)
         call read_var(bn, trim(reg_site%var(k)%name), bs, xb, fb)
         if (.not. (fa .and. fb) .or. any(as /= bs)) then
            write(*,'(4a)') '  missing or reshaped: ', trim(reg_site%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik ; cycle
         end if
         if (any(xa /= xb)) then
            write(*,'(4a)') '  DIFFERS: ', trim(reg_site%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik
         end if
      end do
      call nc_check(nc_close(an), 'close') ; call nc_check(nc_close(bn), 'close')
   end function compare_site_file

   !----- The daily file of the month polygon 2 failed in: five records, polygons 1 and 3's days,   !
   !      polygon 2's first, and the fill value in polygon 2's other days, in every site variable. ---!
   subroutine check_failed_month(path)
      character(len=*), intent(in) :: path
      integer(c_int) :: ncid
      real(c_double), allocatable :: x(:)
      real(c_double) :: fill
      integer(ik) :: nt, j, k, nbad
      logical :: found
      ncid = open_nc(path)
      nt = dim_len(ncid, 'time')
      call check_true('a failed month writes every day the other polygons closed', nt == 5_ik)
      call read_var(ncid, 'cas_temp_site', [nt, NP], x, found)
      call check_true('polygons 1 and 3 hold their five days, polygon 2 its first',                  &
                      found .and. all(x(1:5*NP:NP) < 1.0e30_c_double) .and. all(x(3:5*NP:NP) < 1.0e30_c_double) &
                      .and. x(2) < 1.0e30_c_double)
      nbad = 0_ik
      do j = 1_ik, reg_region%nidx(2)
         k = reg_region%idx_freq(j, 2)
         if (reg_region%var(k)%dim /= DIM_SCALAR) cycle
         call read_var(ncid, trim(reg_region%var(k)%name), [nt, NP], x, found)
         if (.not. found) cycle
         if (reg_region%var(k)%xtype == XTYPE_INT) then
            fill = real(MISSING_INT, c_double)
            if (any(x(NP + 2:nt*NP:NP) /= fill)) nbad = nbad + 1_ik
         else if (any(x(NP + 2:nt*NP:NP) < 1.0e30_c_double)) then
            nbad = nbad + 1_ik
         end if
      end do
      call check_true('the fill value where polygon 2 had failed', nbad == 0_ik)
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_failed_month

   !----- The polygon axis: ids, centres and grid indices of the selected cells. -------------------!
   subroutine check_polygon_axis(path)
      character(len=*), intent(in) :: path
      integer(c_int) :: ncid
      real(c_double), allocatable :: x(:)
      logical :: found
      ncid = open_nc(path)
      call check_true('the region file has a polygon dimension of 3', dim_len(ncid, 'polygon') == NP)
      call read_var(ncid, 'polygon_id', [NP], x, found)
      call check_true('polygon_id holds the cells'' ids', found .and. all(nint(x, ik) == IDS))
      call read_var(ncid, 'lat', [NP], x, found)
      call check_true('lat holds the cell centres', found .and. all(abs(x - LAT) < 1.0e-9_wp))
      call read_var(ncid, 'lon', [NP], x, found)
      call check_true('lon holds the cell centres', found .and. all(abs(x - LON) < 1.0e-9_wp))
      call read_var(ncid, 'row', [NP], x, found)
      call check_true('row and col map back onto the grid', found .and. all(nint(x, ik) == IDS / 18_ik))
      call read_var(ncid, 'col', [NP], x, found)
      call check_true('col is the column index', found .and. all(nint(x, ik) == mod(IDS, 18_ik)))
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_polygon_axis

end program test_region
