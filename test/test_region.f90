! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_region -- a region run is its site runs (docs/configuration.md §Regional runs).     !
!                                                                                          !
! On the synthetic ED_ERA5land archive, a region over a box that crosses 180 degrees,      !
! holds a no-data cell and a mostly-water column, runs its three polygons for twelve days  !
! across a month boundary. Then each polygon is run again as a site at its cell (centre,   !
! orography, UTC), one slow step at a time as the C API steps a run, where the region      !
! steps a month at a time. The region's files must hold, in every variable's `polygon = p` !
! slice, exactly the site run's series, bit for bit, on every tier, cohort and patch       !
! variables included, with the fill value past the polygon's own count; and the polygon    !
! axis must carry the cells' ids and coordinates. The region is run again on four threads  !
! and must give the same files. A region run in two halves, the second restarted from the  !
! first's state file, must write the second half's files exactly as the whole run does,    !
! and a region with a cell masked out must resume from that file, its polygons running on  !
! as they do in the second half. A region freed without finalizing must open again in the  !
! same variable. Every run takes its soil from the archive's [maps] files                  !
! (docs/configuration.md §Input maps), so the three polygons have their own soil depth,    !
! layers and bottom boundaries (a leaky base drains at its cell's f_leak), and each site   !
! run must find the same soil at its cell; the region files carry each polygon's own       !
! layers and map inputs, its drainage factor among them. Last, meds_main (the second       !
! argument) must refuse each config that breaks a region-mode rule, with that rule's       !
! message.                                                                                 !
!                                                                                          !
! The configs are DERIVED from test/region/meds_config_site.toml (the test runs from that  !
! directory): an override block is prepended -- meds_toml returns the FIRST match of a key !
! -- and the keys the archive or region mode reject are dropped. Everything the test       !
! writes goes under the work directory given as the first argument.                        !
!==========================================================================================!
program test_region
   use iso_c_binding,              only : c_int, c_size_t, c_double
   use meds_kinds,                 only : wp, ik
   use meds_test_assert,           only : check_true, test_report
   use meds_test_era5land_archive, only : write_archive, write_map, write_surface_maps, deep_temp_at, &
                                          f_leak_at, NLAT, NLON
   use meds_soil_params,           only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use meds_constants,             only : n_soil_layer_max
   use meds_netcdf_c
   use meds_output_types,          only : output_registry_t, DIM_SCALAR, DIM_COHORT, DIM_PATCH,   &
                                          DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, DIM_SITE, XTYPE_INT, &
                                          MISSING_INT
   use meds_region,                only : meds_region_t, region_open, region_step, region_step_month, &
                                          region_done, region_finalize, region_free, DRIVER_OK,   &
                                          DRIVER_ERR_BOUNDS
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
   integer(ik), parameter :: ROW(NP) = [1_ik, 1_ik, 2_ik], COL(NP) = [17_ik, 0_ik, 17_ik]
   !----- Their soil (meds_test_era5land_archive): 1.95, 1.10 and 2.25 m deep, so 7, 5 and 7 layers on  !
   !      the default ones; sealed (drainage at f = 0), leaky (drainage at the cell's f_leak) and      !
   !      aquifer bottoms. ---------------------------------------------------------------------------!
   integer(ik), parameter :: N_LAYER(NP) = [7_ik, 5_ik, 7_ik]
   real(wp),    parameter :: COLUMN_DEPTH(NP) = [1.95_wp, 1.10_wp, 2.25_wp]
   integer(ik), parameter :: BOTTOM(NP)  = [SOIL_BC_DRAIN, SOIL_BC_DRAIN, SOIL_BC_AQUIFER]

   character(len=512)  :: work, meds_main
   type(meds_region_t) :: reg, site_run
   type(output_registry_t) :: reg_region, reg_site
   character(len=64), allocatable :: files(:), files2(:)
   integer(ik) :: p, st, f, nf, nbad, nf2, f2
   logical     :: ok, ex
   character(len=24) :: idstr
   real(wp)          :: mask(NLON, NLAT, 1)

   call get_command_argument(1, work)
   call get_command_argument(2, meds_main)
   if (len_trim(work) == 0 .or. len_trim(meds_main) == 0)                                           &
      error stop 'test_region: pass the work directory and the meds_main executable'
   call execute_command_line('rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'/out"')
   call write_archive(trim(work)//'/archive')
   call write_surface_maps(trim(work)//'/archive')
   call write_surface_maps(trim(work)//'/archive/gap', gap=[1_ik, 0_ik])     ! no sand at polygon 18

   !----- The region. ----------------------------------------------------------------------------!
   call derive(trim(work)//'/region.toml', region_block(), region=.true.)
   call region_open(trim(work)//'/region.toml', reg, ok, verbose=.false.)
   call check_true('the region opens', ok)
   if (.not. ok) error stop 'test_region: the region did not open'
   call check_true('three polygons: a no-data cell and a mostly-water column are not simulated',     &
                   size(reg%poly) == NP)
   call check_true('polygon ids are the cells'' row-major indices, in row-major order',              &
                   all(reg%poly(:)%id == IDS))
   do p = 1_ik, NP
      associate (s => reg%poly(p)%site%surface(1))
         call check_true('each polygon''s site takes its soil column, bottom boundary and deep '//         &
                         'temperature from its own cell', s%soil%n_active == N_LAYER(p)                     &
                         .and. s%water_bottom_bc == BOTTOM(p) .and. s%soil_column%n_layer_values == N_LAYER(p) &
                         .and. abs(s%deep_temp - deep_temp_at(ROW(p), COL(p))) < 1.0e-9_wp)
      end associate
   end do
   call check_true('the sealed polygon drains at f = 0 and the leaky one at its cell''s f_leak',       &
                   reg%poly(1)%site%surface(1)%drain_factor == 0.0_wp .and.                           &
                   reg%poly(2)%site%surface(1)%drain_factor == f_leak_at(ROW(2), COL(2)))
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
      call region_open(trim(work)//'/site'//trim(idstr)//'.toml', site_run, ok, verbose=.false.)
      if (.not. ok) error stop 'test_region: a site run did not open'
      call check_true('a site run is a region of one polygon', size(site_run%poly) == 1)
      do while (.not. region_done(site_run))
         call region_step(site_run, st)
         if (st /= DRIVER_OK) error stop 'test_region: a site step failed'
      end do
      call region_finalize(site_run, st)
      if (p == 1_ik) reg_site = site_run%out_files%reg
      call region_free(site_run)
   end do

   !----- Two [[sites]] at the first cell (docs/science/forcing.md §8): the first flat at the cell, !
   !      the second 500 m above it on a 20-degree slope facing south. The archive has no sun, so   !
   !      the slope is test_site_terrain's; here the lapse. ----------------------------------------!
   call derive(trim(work)//'/terrain.toml', '[output]'//nl()//'prefix = "terrain"'//nl()//            &
               '[[sites]]'//nl()//'area_fraction = 0.5'//nl()//'[[sites]]'//nl()//                      &
               'area_fraction = 0.5'//nl()//'elevation = '//fmt(ELEV(1) + 500.0_wp)//nl()//             &
               'slope = 20.0'//nl()//'aspect = 180.0'//nl()//site_block(1_ik), region=.false.)
   call run_region(trim(work)//'/terrain.toml')
   call check_terrain(trim(work)//'/out/terrain-D-202101.nc')

   !----- The polygon axis. ----------------------------------------------------------------------!
   call list_files(trim(work)//'/out', files, nf)
   call check_polygon_axis(trim(work)//'/out/region-D-202101.nc')
   call check_polygon_inputs('-D-202101.nc')

   !----- Every region file against the site files: each polygon's slice, every variable. --------!
   nbad = 0_ik
   do f = 1_ik, nf
      if (index(files(f), 'region-') /= 1) cycle
      do p = 1_ik, NP
         write(idstr,'(i0)') IDS(p)
         nbad = nbad + compare_region_file(trim(work)//'/out/'//trim(files(f)), p,                &
                                           trim(work)//'/out/site'//trim(idstr)//files(f)(7:))
      end do
   end do
   call check_true('each region variable''s polygon slice equals its site run, bit for bit, the '//      &
                   'padding past its count the fill value', nbad == 0_ik)
   call check_true('region files: every tier, the sub-daily one included, file for file as a site run',  &
                   count_prefix('region-D-') == 2_ik .and. count_prefix('region-M-') == 2_ik .and.     &
                   count_prefix('region-F-') == count_prefix('site18-F-') .and. count_prefix('region-F-') > 0_ik)

   !----- The same region on four threads: the polygons run side by side, and every file            !
   !      must still equal the site runs bit for bit, slice by slice. In a build without OpenMP this  !
   !      is the serial run again. --------------------------------------------------------------------!
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
      do p = 1_ik, NP
         write(idstr,'(i0)') IDS(p)
         nbad = nbad + compare_region_file(trim(work)//'/out/regt'//files(f)(7:), p,                &
                                           trim(work)//'/out/site'//trim(idstr)//files(f)(7:))
      end do
   end do
   call check_true('on four threads every region file still equals the site runs', nbad == 0_ik)

   !----- A region in two halves, on the archive's 2021 recycled: 2022-12-25 to 2023-01-01            !
   !      writing its state file, then on from it to 2023-01-06. The seam crosses the recycle wrap and  !
   !      the year boundary. The second half writes January's files exactly as the whole run does. ---!
   call derive(trim(work)//'/whole.toml', recycled('whole', '"2023-01-06"', .false.), region=.true.)
   call run_region(trim(work)//'/whole.toml')
   call derive(trim(work)//'/half1.toml', '[output]'//nl()//'state_write = true'//nl()//'state_file_dir = "'// &
               trim(work)//'/state"'//nl()//'state_file_prefix = "half"'//nl()//'state_interval = "1 year"'//nl()// &
               recycled('half1', '"2023-01-01"', .false.), region=.true.)
   call execute_command_line('mkdir -p "'//trim(work)//'/state"')
   call run_region(trim(work)//'/half1.toml')
   call check_state_file(trim(work)//'/state/half-S-20230101000000.nc')
   call derive(trim(work)//'/half2.toml', recycled('half2', '"2023-01-06"', .true.), region=.true.)
   call run_region(trim(work)//'/half2.toml')
   call list_files(trim(work)//'/out', files2, nf2)
   nbad = 0_ik ; f2 = 0_ik
   do f = 1_ik, nf2
      if (index(files2(f), 'whole-') /= 1 .or. index(files2(f), '202301') == 0) cycle
      f2 = f2 + 1_ik
      nbad = nbad + compare_whole_file(trim(work)//'/out/half2'//files2(f)(6:), trim(work)//'/out/'//trim(files2(f)))
   end do
   call check_true('a region restarted from its state file across the recycle wrap writes the whole run''s files', &
                   nbad == 0_ik .and. f2 >= 2_ik)

   !----- A reused run owes no restructuring from the run before it. A run ending on the 1st leaves  !
   !      its boundary's restructuring pending; region_open must clear that, or the next stand is     !
   !      restructured at its first step. It writes a one-polygon state file, which a region refuses  !
   !      below. -------------------------------------------------------------------------------------!
   call derive(trim(work)//'/reuse.toml', '[run]'//nl()//'end_time = "2021-02-01"'//nl()//            &
               '[output]'//nl()//'enabled = false'//nl()//'state_write = true'//nl()//                  &
               'state_file_dir = "'//trim(work)//'/state"'//nl()//'state_file_prefix = "site"'//nl()//  &
               site_block(1_ik), region=.false.)
   call region_open(trim(work)//'/reuse.toml', site_run, ok, verbose=.false.)
   if (.not. ok) error stop 'test_region: the reuse run did not open'
   do while (.not. region_done(site_run))
      call region_step(site_run, st)
   end do
   call check_true('a run ending on the 1st leaves its boundary''s restructuring owed',              &
                   site_run%poly(1)%restructure_pending)
   call region_finalize(site_run, st) ; call region_free(site_run)
   call region_open(trim(work)//'/reuse.toml', site_run, ok, verbose=.false.)
   call check_true('a reused run opens owing no restructuring', ok .and. .not. site_run%poly(1)%restructure_pending)
   !----- ...and its canopy air starts at the bare-ground floor, not the type's default depth. ------!
   call check_true('a bare-ground stand opens at the canopy-air floor',                               &
                   all(site_run%poly(1)%patch%cas(1:site_run%poly(1)%patch%n)%can_depth               &
                       == site_run%cfg%aero%min_canopy_depth))
   call region_finalize(site_run, st) ; call region_free(site_run)

   !----- A state file that cannot start the run stops it, naming the file: one that is not there,  !
   !      and one without a row for each of the run's polygons (a region's in a site run, whose      !
   !      polygon is not among the region's, and a site's in a region). -----------------------------!
   call refused('a region restart from a missing state file stops, saying it cannot be opened',     &
                '[init]'//nl()//'init_mode = 2'//nl()//'restart_file = "'//trim(work)//'/none.nc"'//nl()// &
                region_block(), .true., .false., trim(work)//'/none.nc cannot be opened')
   call refused('a site run refuses a region state file', '[init]'//nl()//'init_mode = 2'//nl()//      &
                'restart_file = "'//trim(work)//'/state/half-S-20230101000000.nc"'//nl()//site_block(1_ik), &
                .false., .false., 'does not hold polygon 0')
   call refused('a region refuses a site state file', '[init]'//nl()//'init_mode = 2'//nl()//          &
                'restart_file = "'//trim(work)//'/state/site-S-20210201000000.nc"'//nl()//region_block(), &
                .true., .false., 'does not hold polygon 35')

   !----- A polygon that fails in mid-month. Polygon 2's soil carbon is made impossible after the    !
   !      first month, so its first February step fails. The others finish the month: polygons 1     !
   !      and 3 close February's five days and polygon 2 its first. The month's I/O phase writes      !
   !      what closed, with the fill value where polygon 2 closed nothing, and the region moves on.   !
   !      The end of the run then writes nothing more: no partial February. ---------------------------!
   call derive(trim(work)//'/fail.toml', '[output]'//nl()//'prefix = "fail"'//nl()//region_block(), &
               region=.true.)
   call region_open(trim(work)//'/fail.toml', reg, ok, verbose=.false.)
   if (.not. ok) error stop 'test_region: the failing region did not open'
   call region_step_month(reg, st)
   call check_true('the month before the failure steps', st == DRIVER_OK)
   reg%poly(2)%patch%soil_carbon(1)%slow_carbon = -1.0_wp
   call region_step_month(reg, st)
   call check_true('the month a polygon fails in returns its status', st == DRIVER_ERR_BOUNDS)
   call check_true('the failed month is consumed: the region moves on', region_done(reg))
   call region_finalize(reg, st)
   call region_free(reg)
   call check_failed_month(trim(work)//'/out/fail-D-202102.nc')
   inquire(file=trim(work)//'/out/fail-M-202102.nc', exist=ex)
   call check_true('after a failure the end of the run writes nothing more', .not. ex)

   !----- The cell mask: a 0/1 map that drops polygon 35 leaves 18 and 53, in their order. ----------!
   mask = 1.0_wp
   mask(18, 2, 1) = 0.0_wp                                   ! (row 1, col 17): polygon 35
   call write_map(trim(work)//'/archive/keep.nc', 'keep', mask, -1.0_wp)
   call derive(trim(work)//'/masked.toml', '[output]'//nl()//'prefix = "masked"'//nl()//'[region]'//nl()// &
               'cell_mask = "'//trim(work)//'/archive/keep.nc"'//nl()//'cell_mask_var = "keep"'//nl()//     &
               region_block(), region=.true.)
   call region_open(trim(work)//'/masked.toml', reg, ok, verbose=.false.)
   call check_true('a region with a cell mask opens', ok)
   if (ok) call check_true('the mask drops polygon 35 and keeps 18 and 53, in order',                    &
                           size(reg%poly) == 2 .and. all(reg%poly(:)%id == [18_ik, 53_ik]))
   if (ok) call region_free(reg)
   !----- Freed without finalizing, the region closed its forcing source and files, so the same   !
   !      variable opens again and steps. ---------------------------------------------------------!
   call region_open(trim(work)//'/masked.toml', reg, ok, verbose=.false.)
   st = -1_ik
   if (ok) call region_step_month(reg, st)
   call check_true('a region freed without finalizing opens again in the same variable, and steps',   &
                   ok .and. st == DRIVER_OK)
   call region_free(reg)
   call write_map(trim(work)//'/archive/offgrid.nc', 'keep', mask, -1.0_wp, lat_shift=0.05_wp)

   !----- A region resumes without a cell: the first half's state file holds polygons 35, 18 and 53,  !
   !      the region masked as above 18 and 53. Each starts from its own row of the file and runs on  !
   !      exactly as it does in the second half, which resumed all three. ----------------------------!
   call derive(trim(work)//'/half3.toml', '[region]'//nl()//'cell_mask = "'//trim(work)//'/archive/keep.nc"'// &
               nl()//'cell_mask_var = "keep"'//nl()//recycled('half3', '"2023-01-06"', .true.), region=.true.)
   call region_open(trim(work)//'/half3.toml', reg, ok, verbose=.false.)
   call check_true('a region resumes from a state file that holds its polygons and more', ok)
   if (.not. ok) error stop 'test_region: the masked region did not resume'
   call check_true('the resumed region holds polygons 18 and 53',                                     &
                   size(reg%poly) == 2 .and. all(reg%poly(:)%id == [18_ik, 53_ik]))
   call check_true('each resumed polygon''s stand is its row of the state file',                     &
                   rows_resumed(trim(work)//'/state/half-S-20230101000000.nc', reg))
   do while (.not. region_done(reg))
      call region_step_month(reg, st)
      if (st /= DRIVER_OK) error stop 'test_region: the resumed masked region failed'
   end do
   call region_finalize(reg, st)
   call region_free(reg)
   nbad = 0_ik ; f2 = 0_ik
   do f = 1_ik, nf2
      if (index(files2(f), 'half2-') /= 1) cycle
      f2 = f2 + 1_ik
      do p = 2_ik, NP                                    ! half3's polygon p - 1 is half2's polygon p
         nbad = nbad + compare_polygons(trim(work)//'/out/half3'//files2(f)(6:), p - 1_ik,           &
                                        trim(work)//'/out/'//trim(files2(f)), p)
      end do
   end do
   call check_true('the resumed polygons run on as they do in the region that resumed all three',     &
                   nbad == 0_ik .and. f2 >= 2_ik)

   !----- Region-mode rules (docs/configuration.md §Regional runs, Rules): each bad config is      !
   !      refused with its own message. The bad keys come first, so they win over the good region  !
   !      block. -------------------------------------------------------------------------------------!
   call refused('the [site] location is refused in a region', region_block(), .false., .false.,     &
                'does not apply to [run].mode = "region"')
   call refused('forcing.max_distance_km is refused in a region',                                  &
                '[forcing]'//nl()//'max_distance_km = 50.0'//nl()//region_block(), .true., .false.,   &
                'does not apply to [run].mode = "region"')
   call refused('a region needs the archive', '[forcing]'//nl()//'format = "default"'//nl()//region_block(),    &
                .true., .true., 'needs forcing.format = "ERA5land"')
   call refused('a region does not start from a census', '[init]'//nl()//'init_mode = 1'//nl()//           &
                'census_file = "census.csv"'//nl()//region_block(),                                           &
                .true., .false., 'not a census')
   call refused('a fast-tier window may not straddle two slow steps',                               &
                '[output]'//nl()//'fast_interval_steps = 5'//nl()//site_block(1_ik), .false., .false., &
                'output.fast_interval_steps must divide')
   !----- An io_config override is held to the tiers' sources as the defaults are: a variable    !
   !      with no sub-daily value on the F stream, one that exists only sub-daily on a coarse tier,  !
   !      and `true` on a variable off by default, which has no default tiers to restore. -----------!
   call write_text(trim(work)//'/io_f_coarse.toml', '[variables]'//nl()//'mean_dbh_co = "F D"')
   call refused('an io_config override may not put a coarse variable on the F stream',            &
                '[output]'//nl()//'io_config = "'//trim(work)//'/io_f_coarse.toml"'//nl()//          &
                site_block(1_ik), .false., .false., 'has no sub-daily value')
   call write_text(trim(work)//'/io_fast_daily.toml', '[variables]'//nl()//'mean_reco_pa = "D"')
   call refused('an io_config override may not put a FAST-staged variable on a coarse tier',      &
                '[output]'//nl()//'io_config = "'//trim(work)//'/io_fast_daily.toml"'//nl()//        &
                site_block(1_ik), .false., .false., 'exists only sub-daily')
   call write_text(trim(work)//'/io_true_off.toml', '[variables]'//nl()//'cas_temp_py = true')
   call refused('an io_config true may not switch on a variable that is off by default',          &
                '[output]'//nl()//'io_config = "'//trim(work)//'/io_true_off.toml"'//nl()//          &
                site_block(1_ik), .false., .false., 'is off by default')
   !----- The selection's own rules (docs/science/diagnostics.md §5): a line that matches nothing, !
   !      two lines as specific as each other that disagree, a tier none of a line's variables can be !
   !      on, and the axes' bookkeeping named; and the retired switches are unknown keys. -----------!
   call refuse_io('an io_config pattern that matches no variable is refused', '"*_xyz" = false',            &
                  'matches no variable')
   call refuse_io('two io_config lines as specific as each other that disagree are refused',              &
                  '"mean_a*" = "D"'//nl()//'"*agb_py" = "M"', 'make one more specific')
   call refuse_io('an io_config line asking for a tier none of its variables can be on is refused',       &
                  '"*_co" = "Y"', 'none of its variables can be on Y')
   call refuse_io('the axes'' bookkeeping is not selected', 'owner_patch = "D"', 'is written with its axis')
   call refused('a retired [output] switch is an unknown key', '[output]'//nl()//'axes_cohort = false'//nl()// &
                site_block(1_ik), .false., .false., 'is not a key MEDS reads')
   call refused('a cell mask must be on the forcing grid', '[region]'//nl()//'cell_mask = "'//            &
                trim(work)//'/archive/offgrid.nc"'//nl()//'cell_mask_var = "keep"'//nl()//region_block(),  &
                .true., .false., 'is not on the forcing grid')
   call refused('a [region] block needs region mode', '[run]'//nl()//'mode = "site"'//nl()//          &
                site_block(1_ik)//'[region]'//nl()//'land_fraction_min = 0.5'//nl(), .false., .false.,  &
                'needs [run].mode = "region"')
   call refused('the box must have S < N', '[region]'//nl()//'box_nwse = [20.0, 150.0, 40.0, -150.0]'//nl()// &
                region_block(), .true., .false., 'region.box_nwse needs')
   call refused('the mode is "site" or "region"', '[run]'//nl()//'mode = "globe"'//nl()//region_block(), &
                .true., .false., 'run.mode = "globe" is not one of site | region')
   !----- The terrain-lapse keys need the lapse (docs/configuration.md §The forcing file). ------!
   call refused('the lapse keys are refused with the lapse off',                                    &
                '[site]'//nl()//'apply_elevation_lapse = false'//nl()//site_block(1_ik), .false., .false., &
                'applies only with site.apply_elevation_lapse = true')
   !----- The [maps] files (docs/configuration.md §Input maps): a region needs the soil maps, a   !
   !      key they provide is refused, a run cell they leave empty is named, and Chave's E serves !
   !      only the Chave heights. ----------------------------------------------------------------!
   call refused('a region needs the soil maps', '[maps]'//nl()//'soil_depth_map = ""'//nl()//        &
                'soil_texture_map = ""'//nl()//'soil_chemistry_map = ""'//nl()//'[soil]'//nl()//         &
                'deep_temp = 284.75'//nl()//region_block(), .true., .false., 'needs the soil maps')
   call refused('a soil key the maps provide is refused', '[soil]'//nl()//'theta_sat = 0.40'//nl()//    &
                site_block(1_ik), .false., .false., 'comes from the [maps] soil maps')
   call refused('a deep temperature beside the maps is refused', '[soil]'//nl()//'deep_temp = 290.0'//nl()// &
                region_block(), .true., .false., 'comes from the [maps] soil maps')
   call refused('a run cell the soil maps leave empty is named', '[maps]'//nl()//'soil_texture_map = "'// &
                trim(work)//'/archive/gap/soil_texture.nc"'//nl()//region_block(), .true., .false.,         &
                'sand has no value at the cell (35.00, -180.00)')
   call refused('a chave_e map needs the Chave height allometry', '[maps]'//nl()//'chave_e_map = "'//    &
                trim(work)//'/archive/chave_e.nc"'//nl()//site_block(1_ik), .false., .false.,               &
                'maps.chave_e_map needs a PFT with allometry.height_allometry = "chave2014"')
   !----- [[sites]] (docs/configuration.md §Several sites): a region takes one site per cell, and  !
   !      with the soil maps a site takes the cell's soil. ----------------------------------------!
   call refused('a region refuses [[sites]]', '[[sites]]'//nl()//'area_fraction = 1.0'//nl()//         &
                region_block(), .true., .false., 'a region refuses [[sites]]')
   call refused('a [[sites]] soil key beside the maps is refused', '[[sites]]'//nl()//               &
                'area_fraction = 1.0'//nl()//'depth = 1.0'//nl()//site_block(1_ik), .false., .false.,   &
                'every site takes its cell''s soil')
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
          '[forcing]'//nl()//'format = "ERA5land"'//nl()//'data_path = "'//trim(work)//'/archive"'//nl()//    &
          'recycle = false'//nl()//                                                                 &
          '[output]'//nl()//'dir = "'//trim(work)//'/out"'//nl()//'state_write = false'//nl()//       &
          '[output.monthly]'//nl()//'enabled = true'//nl()//'file_chunk = "month"'//nl()//            &
          '[maps]'//nl()//'soil_depth_map = "'//trim(work)//'/archive/soil_depth.nc"'//nl()//      &
          'soil_texture_map = "'//trim(work)//'/archive/soil_texture.nc"'//nl()//                     &
          'soil_chemistry_map = "'//trim(work)//'/archive/soil_chemistry.nc"'//nl()
   end function common_block

   function region_block() result(b)
      character(len=:), allocatable :: b
      b = '[run]'//nl()//'mode = "region"'//nl()//'[output]'//nl()//'prefix = "region"'//nl()//     &
          '[region]'//nl()//'box_nwse = [40.0, 150.0, 20.0, -150.0]'//nl()//                        &
          'land_fraction_min = 0.5'//nl()//common_block()
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

   !----- A small text file in the work directory. --------------------------------------------------!
   subroutine write_text(path, text)
      character(len=*), intent(in) :: path, text
      integer :: u
      open(newunit=u, file=path, status='replace', action='write')
      write(u, '(a)') text
      close(u)
   end subroutine write_text

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

   !----- An io_config of these [variables] lines, refused with the message given. -----------------!
   subroutine refuse_io(what, lines, message)
      character(len=*), intent(in) :: what, lines, message
      call write_text(trim(work)//'/io_refused.toml', '[variables]'//nl()//lines)
      call refused(what, '[output]'//nl()//'io_config = "'//trim(work)//'/io_refused.toml"'//nl()//     &
                   site_block(1_ik), .false., .false., message)
   end subroutine refuse_io

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
      drop = drop .or. (sec == 'soil' .and. key == 'deep_temp')       ! the soil maps give it
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
      case (DIM_SOIL)       ; sh = [dim_len(ncid, 'site'), dim_len(ncid, 'soil')]
      case (DIM_PFT)        ; sh = [dim_len(ncid, 'pft')]
      case (DIM_SIZE)       ; sh = [dim_len(ncid, 'dbh_class')]
      case (DIM_SOIL_PATCH) ; sh = [dim_len(ncid, 'patch'), dim_len(ncid, 'soil')]
      case (DIM_SITE)
         sh = [dim_len(ncid, 'site')]
      end select
      if (with_polygon) sh = [dim_len(ncid, 'polygon'), sh]
      sh = [nt, sh]
   end function var_shape

   !----- Region file vs one polygon's site file: the polygon's slice of every live variable. Its    !
   !      first entries are the site's (a cohort or patch axis as long as the site's own count); the  !
   !      rest of a padded axis is the fill value. The region file has every variable the site has;  !
   !      the site's file has a polygon axis of one. --------------------------------------------------!
   integer(ik) function compare_region_file(rpath, p, spath) result(nbad)
      character(len=*), intent(in) :: rpath, spath
      integer(ik),      intent(in) :: p
      integer(c_int) :: rn, sn
      integer(ik)    :: t, j, k, nt, na, ns, np_, it, base
      integer(ik), allocatable :: rs(:), ss(:)
      real(c_double), allocatable :: xr(:), xs(:)
      logical :: fr, fs, pad_ok
      nbad = 0_ik
      t = tier_of(rpath)
      rn = open_nc(rpath) ; sn = open_nc(spath)
      do j = 1_ik, reg_site%nidx(t)
         k = reg_site%idx_freq(j, t)
         rs = var_shape(rn, reg_site%var(k)%dim, .true.)
         ss = var_shape(sn, reg_site%var(k)%dim, .true.)
         call read_var(rn, trim(reg_site%var(k)%name), rs, xr, fr)
         call read_var(sn, trim(reg_site%var(k)%name), ss, xs, fs)
         nt = rs(1) ; np_ = rs(2) ; na = product(rs(3:)) ; ns = product(ss(3:))   ! C order: axis fastest
         if (.not. (fr .and. fs) .or. rs(1) /= ss(1) .or. ns > na) then
            write(*,'(5a)') '  missing or mismatched: ', trim(reg_site%var(k)%name), ' in ', rpath, ''
            nbad = nbad + 1_ik ; cycle
         end if
         do it = 0_ik, nt - 1_ik
            base = it*np_*na + (p - 1_ik)*na
            if (reg_site%var(k)%xtype == XTYPE_INT) then
               pad_ok = all(xr(base + ns + 1_ik : base + na) == real(MISSING_INT, c_double))
            else
               pad_ok = all(xr(base + ns + 1_ik : base + na) > 1.0e30_c_double)
            end if
            if (any(xr(base + 1_ik : base + ns) /= xs(it*ns + 1_ik : (it + 1_ik)*ns)) .or. .not. pad_ok) then
               write(*,'(3a,i0,2a)') '  DIFFERS: ', trim(reg_site%var(k)%name), ' polygon ', IDS(p), ' in ', rpath
               nbad = nbad + 1_ik ; exit
            end if
         end do
      end do
      call nc_check(nc_close(rn), 'close') ; call nc_check(nc_close(sn), 'close')
   end function compare_region_file

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
      call read_var(ncid, 'mean_cas_temp_py', [nt, NP], x, found)
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

   !----- A region on the archive's 2021, recycled, from 2022-12-25 to `last`, its files named       !
   !      `prefix`; `restart` starts it from the first half's state file at 2023-01-01. -----------------!
   function recycled(prefix, last, restart) result(b)
      character(len=*), intent(in) :: prefix, last
      logical,          intent(in) :: restart
      character(len=:), allocatable :: b
      b = '[run]'//nl()//'start_time = "2022-12-25"'//nl()//'end_time = '//last//nl()//                   &
          '[forcing]'//nl()//'recycle = true'//nl()//'recycle_start = "2021-01-01 01:00:00"'//nl()//       &
          'recycle_end = "2022-01-01 01:00:00"'//nl()//'[output]'//nl()//'prefix = "'//prefix//'"'//nl()
      if (restart) b = b//'[init]'//nl()//'init_mode = 2'//nl()//'restart_file = "'//trim(work)//          &
                       '/state/half-S-20230101000000.nc"'//nl()
      b = b//region_block()
   end function recycled

   !----- Open, step and close a region from its config. ------------------------------------------!
   subroutine run_region(path)
      character(len=*), intent(in) :: path
      type(meds_region_t) :: r
      integer(ik) :: s
      logical     :: opened
      call region_open(path, r, opened, verbose=.false.)
      if (.not. opened) error stop 'test_region: a region did not open'
      do while (.not. region_done(r))
         call region_step_month(r, s)
         if (s /= DRIVER_OK) error stop 'test_region: a region month failed'
      end do
      call region_finalize(r, s)
      call region_free(r)
   end subroutine run_region

   !----- Two region files of the same registry: every live variable of the tier, every polygon. ----!
   integer(ik) function compare_whole_file(apath, bpath) result(nbad)
      character(len=*), intent(in) :: apath, bpath
      integer(c_int) :: an, bn
      integer(ik)    :: t, j, k
      integer(ik), allocatable :: as(:), bs(:)
      real(c_double), allocatable :: xa(:), xb(:)
      logical :: fa, fb, ex
      nbad = 0_ik
      inquire(file=apath, exist=ex)
      if (.not. ex) then
         write(*,'(2a)') '  missing restarted file ', apath ; nbad = 1_ik ; return
      end if
      t = tier_of(bpath)
      an = open_nc(apath) ; bn = open_nc(bpath)
      do j = 1_ik, reg_region%nidx(t)
         k = reg_region%idx_freq(j, t)
         as = var_shape(an, reg_region%var(k)%dim, .true.)
         bs = var_shape(bn, reg_region%var(k)%dim, .true.)
         call read_var(an, trim(reg_region%var(k)%name), as, xa, fa)
         call read_var(bn, trim(reg_region%var(k)%name), bs, xb, fb)
         if (.not. (fa .and. fb) .or. any(as /= bs)) then
            write(*,'(4a)') '  missing or reshaped: ', trim(reg_region%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik ; cycle
         end if
         if (any(xa /= xb)) then
            write(*,'(4a)') '  DIFFERS after the restart: ', trim(reg_region%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik
         end if
      end do
      call nc_check(nc_close(an), 'close') ; call nc_check(nc_close(bn), 'close')
   end function compare_whole_file

   !----- Polygon pa of one region file against polygon pb of another of the same tier: every live   !
   !      variable, entry for entry where both files have the entry, and the fill value where only one !
   !      does (each file's cohort, patch and soil axes are its largest polygon's). --------------------!
   integer(ik) function compare_polygons(apath, pa, bpath, pb) result(nbad)
      character(len=*), intent(in) :: apath, bpath
      integer(ik),      intent(in) :: pa, pb
      integer(c_int) :: an, bn
      integer(ik)    :: t, j, k, it, i1, i2, ea(2), eb(2)
      integer(ik), allocatable :: as(:), bs(:)
      real(c_double), allocatable :: xa(:), xb(:)
      real(c_double) :: va, vb
      logical :: fa, fb, ex, in_a, in_b, int_var, same
      nbad = 0_ik
      inquire(file=apath, exist=ex)
      if (.not. ex) then
         write(*,'(2a)') '  missing resumed file ', apath ; nbad = 1_ik ; return
      end if
      t = tier_of(bpath)
      an = open_nc(apath) ; bn = open_nc(bpath)
      do j = 1_ik, reg_region%nidx(t)
         k = reg_region%idx_freq(j, t)
         as = var_shape(an, reg_region%var(k)%dim, .true.)
         bs = var_shape(bn, reg_region%var(k)%dim, .true.)
         call read_var(an, trim(reg_region%var(k)%name), as, xa, fa)
         call read_var(bn, trim(reg_region%var(k)%name), bs, xb, fb)
         if (.not. (fa .and. fb) .or. as(1) /= bs(1)) then
            write(*,'(4a)') '  missing or mismatched: ', trim(reg_region%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik ; cycle
         end if
         ea = 1_ik ; eb = 1_ik                           ! the axes after (time, polygon), at most two
         ea(1:size(as) - 2) = as(3:) ; eb(1:size(bs) - 2) = bs(3:)
         int_var = reg_region%var(k)%xtype == XTYPE_INT
         same = .true.
         do it = 0_ik, as(1) - 1_ik
            do i1 = 0_ik, max(ea(1), eb(1)) - 1_ik
               do i2 = 0_ik, max(ea(2), eb(2)) - 1_ik
                  in_a = i1 < ea(1) .and. i2 < ea(2) ; in_b = i1 < eb(1) .and. i2 < eb(2)
                  if (in_a) va = xa(((it*as(2) + pa - 1_ik)*ea(1) + i1)*ea(2) + i2 + 1_ik)
                  if (in_b) vb = xb(((it*bs(2) + pb - 1_ik)*eb(1) + i1)*eb(2) + i2 + 1_ik)
                  if (in_a .and. in_b) then
                     same = same .and. va == vb
                  else if (in_a) then
                     same = same .and. is_fill(va, int_var)
                  else
                     same = same .and. is_fill(vb, int_var)
                  end if
               end do
            end do
         end do
         if (.not. same) then
            write(*,'(4a)') '  DIFFERS after the resume: ', trim(reg_region%var(k)%name), ' in ', apath
            nbad = nbad + 1_ik
         end if
      end do
      call nc_check(nc_close(an), 'close') ; call nc_check(nc_close(bn), 'close')
   end function compare_polygons

   logical function is_fill(x, int_var)
      real(c_double), intent(in) :: x
      logical,        intent(in) :: int_var
      if (int_var) then
         is_fill = x == real(MISSING_INT, c_double)
      else
         is_fill = x > 1.0e30_c_double
      end if
   end function is_fill

   !----- Each polygon of a resumed region against its row of the state file it resumed from, the   !
   !      file's polygons being IDS: its counts, each cohort's plant density, diameter and carbon     !
   !      pools, and each patch's area and soil water, matched by their persistent ids. --------------!
   logical function rows_resumed(path, r) result(same)
      character(len=*),    intent(in) :: path
      type(meds_region_t), intent(in) :: r
      integer(c_int) :: ncid
      real(c_double), allocatable :: meta(:), gid(:), nplant(:), dbh(:), leafc(:), woodc(:)
      real(c_double), allocatable :: pgid(:), area(:), theta(:)
      integer(ik), allocatable :: row_ids(:)
      integer(ik) :: q, row, nc, np_, ncoh, npat, nsoil, i, j, base
      logical :: found(9)
      ncid = open_nc(path)
      nc = dim_len(ncid, 'cohort') ; np_ = dim_len(ncid, 'patch') ; nsoil = dim_len(ncid, 'soil_layer')
      call read_var(ncid, 'meta_int',         [NP, 11_ik],       meta,   found(1))
      call read_var(ncid, 'global_cohort_id', [NP, nc],          gid,    found(2))
      call read_var(ncid, 'nplant_co',           [NP, nc],          nplant, found(3))
      call read_var(ncid, 'dbh_co',              [NP, nc],          dbh,    found(4))
      call read_var(ncid, 'leaf_carbon_co',      [NP, nc],          leafc,  found(5))
      call read_var(ncid, 'wood_carbon_co',      [NP, nc],          woodc,  found(6))
      call read_var(ncid, 'global_patch_id',  [NP, np_],         pgid,   found(7))
      call read_var(ncid, 'area_pa',       [NP, np_],         area,   found(8))
      call read_var(ncid, 'soil_theta_layer_pa',       [NP, np_, nsoil],  theta,  found(9))
      call nc_check(nc_close(ncid), 'close')
      same = all(found)
      if (.not. same) return
      do q = 1_ik, size(r%poly, kind=ik)
         row = findloc(IDS, r%poly(q)%id, dim=1) - 1_ik
         ncoh = nint(meta(row*11_ik + 1_ik), ik) ; npat = nint(meta(row*11_ik + 2_ik), ik)
         associate (c => r%poly(q)%cohort, pt => r%poly(q)%patch)
            same = same .and. c%n == ncoh .and. pt%n == npat
            if (.not. same) return
            row_ids = nint(gid(row*nc + 1_ik : row*nc + ncoh), ik)
            do i = 1_ik, ncoh
               j = findloc(row_ids, c%global_id(i), dim=1)
               if (j == 0_ik) then ; same = .false. ; return ; end if
               base = row*nc + j
               same = same .and. c%nplant(i) == nplant(base) .and. c%dbh(i) == dbh(base)            &
                           .and. c%leaf_carbon(i) == leafc(base) .and. c%wood_carbon(i) == woodc(base)
            end do
            row_ids = nint(pgid(row*np_ + 1_ik : row*np_ + npat), ik)
            do i = 1_ik, npat
               j = findloc(row_ids, pt%global_id(i), dim=1)
               if (j == 0_ik) then ; same = .false. ; return ; end if
               base = (row*np_ + j - 1_ik)*nsoil
               same = same .and. pt%area(i) == area(row*np_ + j)                                    &
                           .and. all(pt%soil_water(i)%theta(1:nsoil) == theta(base + 1_ik : base + nsoil))
            end do
         end associate
      end do
   end function rows_resumed

   !----- The region's state file: a polygon axis holding the three polygons, each with its own     !
   !      counts in meta_int, on cohort and patch axes as long as the largest polygon's. -------------!
   subroutine check_state_file(path)
      character(len=*), intent(in) :: path
      integer(c_int) :: ncid
      real(c_double), allocatable :: pid(:), meta(:)
      logical :: found
      integer(ik) :: p, ncoh_max, npat_max
      ncid = open_nc(path)
      call check_true('the region state file has a polygon axis of 3', dim_len(ncid, 'polygon') == NP)
      call read_var(ncid, 'polygon_id', [NP], pid, found)
      call check_true('its polygon_id holds the polygons'' ids', found .and. all(nint(pid, ik) == IDS))
      call read_var(ncid, 'meta_int', [NP, 11_ik], meta, found)
      ncoh_max = 0_ik ; npat_max = 0_ik
      do p = 1_ik, NP
         ncoh_max = max(ncoh_max, nint(meta((p - 1_ik) * 11_ik + 1_ik), ik))
         npat_max = max(npat_max, nint(meta((p - 1_ik) * 11_ik + 2_ik), ik))
      end do
      call check_true('each polygon''s counts are its row of meta_int; the axes are the largest polygon''s', &
                      found .and. dim_len(ncid, 'cohort') == max(ncoh_max, 1_ik) .and.                      &
                      dim_len(ncid, 'patch') == max(npat_max, 1_ik))
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_state_file

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

   !----- The two sites' forcing: the flat one at the cell has the polygon's air exactly, the one !
   !      500 m up is 500 m x 0.0065 K/m colder, at a lower pressure, under less longwave. -------!
   subroutine check_terrain(path)
      character(len=*), intent(in) :: path
      integer(c_int) :: ncid
      integer(ik)    :: nt
      real(c_double), allocatable :: t(:), tp(:), p(:), lw(:)
      logical :: ft, ftp, fp, flw
      ncid = open_nc(path)
      nt = dim_len(ncid, 'time')
      call check_true('a run of two sites has a site axis of two', dim_len(ncid, 'site') == 2_ik)
      call read_var(ncid, 'mean_air_temp_si', [nt, 1_ik, 2_ik], t, ft)
      call read_var(ncid, 'mean_air_temp_py', [nt, 1_ik], tp, ftp)
      call read_var(ncid, 'mean_psurf_si', [nt, 1_ik, 2_ik], p, fp)
      call read_var(ncid, 'mean_lwdown_si', [nt, 1_ik, 2_ik], lw, flw)
      call check_true('the flat site at the cell has the polygon''s air temperature exactly',        &
                      ft .and. ftp .and. all(t(1::2) == tp))
      call check_true('the site 500 m up is 3.25 K colder',                                          &
                      ft .and. all(abs(t(2::2) - t(1::2) + 3.25_c_double) < 1.0e-6_c_double))
      call check_true('at a lower pressure', fp .and. all(p(2::2) < p(1::2)))
      call check_true('under less longwave', flw .and. all(lw(2::2) < lw(1::2)))
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_terrain

   !----- Each polygon's own layers and map inputs, on its one site. Its row of soil_z is its site     !
   !      run's soil_z down to its own bottom layer, and the fill value below. Its soil_z_bnds run     !
   !      from the surface, layer to layer, down to its soil depth, with each node midway. Polygons 1  !
   !      and 3 both have seven layers, ending at 1.95 and 2.25 m, so their bottom nodes differ. ------!
   subroutine check_polygon_inputs(suffix)
      character(len=*), intent(in) :: suffix
      integer(c_int) :: ncid, sn, vid
      real(c_double), allocatable :: z(:), b(:), zs(:), depth(:), bc(:), tdeep(:), drain(:)
      real(c_double) :: top(n_soil_layer_max), bot(n_soil_layer_max)
      logical     :: fz, fb, fd, fc, ft, fs, fr, same_z, edges
      integer(ik) :: p, ns, n, k, base
      character(len=24) :: id
      ncid = open_nc(trim(work)//'/out/region'//suffix)
      ns = dim_len(ncid, 'soil')
      call check_true('each polygon has one site', dim_len(ncid, 'site') == 1_ik)
      call read_var(ncid, 'soil_z', [NP, 1_ik, ns], z, fz)
      call read_var(ncid, 'soil_z_bnds', [NP, 1_ik, ns, 2_ik], b, fb)
      call read_var(ncid, 'soil_depth_si', [NP, 1_ik], depth, fd)
      call read_var(ncid, 'bottom_bc_si', [NP, 1_ik], bc, fc)
      call read_var(ncid, 'deep_temp_si', [NP, 1_ik], tdeep, ft)
      call read_var(ncid, 'drain_factor_si', [NP, 1_ik], drain, fr)
      same_z = fz .and. ns == n_soil_layer_max
      edges  = fb .and. fd .and. ns == n_soil_layer_max
      do p = 1_ik, NP
         if (.not. (same_z .and. edges)) exit
         n = N_LAYER(p) ; base = (p - 1_ik) * ns
         write(id,'(i0)') IDS(p)
         sn = open_nc(trim(work)//'/out/site'//trim(id)//suffix)
         call read_var(sn, 'soil_z', [1_ik, 1_ik, ns], zs, fs)
         call nc_check(nc_close(sn), 'close')
         same_z = fs .and. all(z(base + 1_ik : base + n) == zs(1:n)) .and.                           &
                  all(z(base + n + 1_ik : base + ns) > 1.0e30_c_double)      ! the loop stops at a miss
         do k = 1_ik, ns
            top(k) = b(2_ik * (base + k) - 1_ik) ; bot(k) = b(2_ik * (base + k))
         end do
         edges = top(1) == 0.0_c_double .and. all(top(2:n) == bot(1:n - 1_ik)) .and.                 &
                 all(z(base + 1_ik : base + n) == 0.5_c_double * (top(1:n) + bot(1:n))) .and.         &
                 bot(n) == -depth(p) .and. abs(depth(p) - COLUMN_DEPTH(p)) < 1.0e-6_wp .and.          &
                 all(top(n + 1_ik : ns) > 1.0e30_c_double) .and. all(bot(n + 1_ik : ns) > 1.0e30_c_double)
      end do
      call check_true('each polygon''s row of soil_z is its site run''s, down to its own bottom layer, '// &
                      'and the fill value below', same_z)
      call check_true('the two seven-layer polygons keep their own bottom nodes',                     &
                      fz .and. z(7) /= z(2_ik * ns + 7_ik))
      call check_true('soil_z_bnds run from the surface, layer to layer, to each polygon''s soil_depth', edges)
      call check_true('bottom_bc and deep_temp are each polygon''s own',                             &
                      fc .and. ft .and. all(nint(bc, ik) == BOTTOM) .and.                             &
                      all([(abs(tdeep(p) - deep_temp_at(ROW(p), COL(p))) < 1.0e-9_wp, p = 1_ik, NP)]))
      call check_true('drain_factor is each polygon''s own: 0 sealed, its f_leak leaky, the fill value '// &
                      'on an aquifer', fr .and. drain(1) == 0.0_c_double .and.                        &
                      drain(2) == real(f_leak_at(ROW(2), COL(2)), c_double) .and. drain(3) > 1.0e30_c_double)
      call check_true('no chave_e without the Chave height allometry',                               &
                      nc_inq_varid_f(ncid, 'chave_e_py', vid) /= NC_NOERR)
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_polygon_inputs

end program test_region
