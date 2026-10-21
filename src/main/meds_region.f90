! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_region -- a run: one polygon or many, in one process (docs/science/order_of_processes.md  !
! §1; the configuration is docs/configuration.md). Every run goes through here: `meds_main`      !
! steps it a month at a time, and the C-API shim `meds_c_api_run` (Python's `meds.model`) one    !
! slow step at a time.                                                                            !
!                                                                                          !
!   type(meds_region_t) :: reg                                                                   !
!   call region_open('meds_config_main.toml', reg, ok)                                            !
!   do while (.not. region_done(reg)) ; call region_step_month(reg, status) ; end do               !
!   call region_finalize(reg, status) ; call region_free(reg)                                      !
!                                                                                          !
! [run].mode chooses the polygons: "site" (the default) is one polygon at the [site] location,     !
! on either forcing format; "region" is every selected cell of a box of the ED_ERA5land archive.  !
! Each polygon is stepped by polygon_step. The forcing source loads a month for every cell once;   !
! then each polygon runs through it, reading its own cell; then the I/O phase writes every         !
! polygon's queued records, one hyperslab per variable per record. Nothing inside the polygon      !
! loop touches a file, so the polygons of a month run side by side on [run].n_threads threads,     !
! each with a single-threaded patch loop; one polygon gives its patch loop the threads instead.   !
!                                                                                          !
! Output and state files always carry a `polygon` axis, one row per polygon, the cohort and       !
! patch axes padded to the largest polygon's. A NaN state or a state out of its bounds is a        !
! STATUS RETURN, not an `error stop`: this module is compiled into the shared library Python       !
! loads, and `error stop` there takes the interpreter down with it. `meds_main` turns the status  !
! into an `error stop`.                                                                             !
!==========================================================================================!
module meds_region
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : day_sec, yr_day
   use meds_numerics,               only : ascending_order
   use meds_config,                 only : meds_config_t, INIT_RESTART,                                &
                                           RUN_MODE_REGION, soil_maps_on, site_config, apply_map_inputs
   use meds_config_check,           only : validate_config
   use meds_config_io,              only : load_meds_config, write_pft_params_csv
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days, time_to_string, &
                                           years_between
   use meds_polygon_state_vars,     only : polygon_t, free_polygon
   use meds_init_restart,           only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                           restart_file_slots
   use meds_init,                   only : init_polygon
   use meds_met_types,              only : met_source_t, met_grid_t, met_cells_t
   use meds_forcing_params,         only : MET_BACKEND_ED_ERA5LAND
   use meds_met,             only : met_open, met_close, met_prefetch, MET_OK
   use meds_maps,                   only : apply_cell_mask, map_status_message, MAP_OK, read_map_inputs
   use meds_maps_params,            only : map_inputs_t
   use meds_met_grid,               only : select_box_cells, select_site_cell, GRID_OK, GRID_ERR_NO_CELL
   use meds_met_era5land,           only : era5land_grid
   use meds_output_reduce,      only : print_summary, total_area, total_agb, total_lai
   use meds_polygon,                only : polygon_run_t, polygon_prepare, polygon_step, polygon_report,  &
                                           open_output_files, attach_output, ensure_output_dir,   &
                                           DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, &
                                           DRIVER_ERR_BOUNDS
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_manager,         only : output_serialize, output_close
   use meds_toml,                   only : toml_write_record
   implicit none
   private

   public :: meds_region_t, region_open, region_step, region_step_month, region_done, region_finalize
   public :: region_free
   !----- The status codes live with the step (meds_polygon); re-exported for the run's callers. ---!
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_BOUNDS

   !----- How far the patch areas of a polygon may sum from one at the end of a run before the run  !
   !      reports DRIVER_ERR_AREA. ------------------------------------------------------------------!
   real(wp), parameter :: AREA_TOL = 1.0e-5_wp

   !----- Everything a run needs between steps. No state hides in module scope, so two runs can be  !
   !      open at once (which is what the C-API's handle registry does). --------------------------!
   type :: meds_region_t
      type(meds_config_t)   :: cfg
      type(met_source_t)    :: met_src          !< one reader for every polygon (opened only if forcing_on)
      type(output_files_t)  :: out_files           !< the run's files (built only if output.enabled)
      !----- Each polygon as two records side by side: its state (the stand) and what the run needs  !
      !      to step and check it (meds_polygon). -------------------------------------------------!
      type(polygon_t),     allocatable :: poly(:)
      type(polygon_run_t), allocatable :: run(:)
      !----- Each polygon's share of the files, one per polygon. A contiguous array of its own, so  !
      !      the I/O phase passes it whole: a section run(:)%bufs of a type with allocatable         !
      !      components makes gfortran copy it through a temporary whose copy-out dangles them. ----!
      type(output_buffers_t), allocatable :: out_bufs(:)
      type(meds_time_t)     :: now
      integer(ik)           :: step_days = 1_ik
      integer(ik)           :: istep = 0_ik       !< slow steps taken
      integer(ik)           :: iyear = 0_ik       !< years stepped, for [state].interval_years
      !----- The time spent loading the forcing before each month (met_prefetch), over the run and   !
      !      this year, and the archive month files read before this year: the year summary and the  !
      !      end of the run report them, so a run shows what its forcing reads cost. -----------------!
      real(wp)              :: read_seconds = 0.0_wp, year_read_seconds = 0.0_wp
      integer(ik)           :: year_loads0 = 0_ik
      !----- Last month's wall time per polygon: the next month starts the costliest first, so a   !
      !      slow polygon does not leave the other threads idle at the end of the month. ------------!
      real(wp),    allocatable :: cost(:)
      logical               :: verbose = .true.
   end type meds_region_t

contains

   pure logical function region_done(reg)
      type(meds_region_t), intent(in) :: reg
      region_done = .not. time_lt(reg%now, reg%cfg%end_time)
   end function region_done

   !---------------------------------------------------------------------------------------!
   ! region_open -- config -> the polygons (one at the [site], or the box's cells) -> the shared  !
   ! forcing source -> each polygon's initial community and preparation -> the output files:      !
   ! everything before the first step.                                                            !
   !---------------------------------------------------------------------------------------!
   subroutine region_open(path, reg, ok, verbose)
      character(len=*),    intent(in)    :: path
      type(meds_region_t), intent(inout) :: reg
      logical,             intent(out)   :: ok
      logical, optional,   intent(in)    :: verbose
      type(met_cells_t)          :: cells
      type(met_grid_t)           :: grid
      integer(ik)                :: st, p, n
      integer(ik), allocatable   :: ids(:), slot(:)
      type(map_inputs_t)         :: inputs
      type(meds_config_t)        :: pcfg
      type(meds_config_t), allocatable :: site_cfg(:)     !< [[sites]]: each site's configuration ...
      real(wp),            allocatable :: site_area(:)    !< ... and share of the polygon
      character(len=:), allocatable :: msg
      logical                    :: fast_found, region, started

      ok = .false.
      if (present(verbose)) reg%verbose = verbose
      !----- The C-API reuses a freed registry slot, so the counters start again here. ------------!
      reg%istep = 0_ik
      reg%iyear = 0_ik
      reg%read_seconds = 0.0_wp
      reg%year_read_seconds = 0.0_wp
      reg%year_loads0 = 0_ik

      call load_meds_config(trim(path), reg%cfg)
      if (reg%verbose) write(*,'(2a)') ' config: ', trim(path)
      associate (cfg => reg%cfg)
      region = cfg%mode == RUN_MODE_REGION

      if (region) then
         !----- A region steps each polygon on an OpenMP worker thread, whose stack is OMP_STACKSIZE   !
         !      (4 MB by default with ifx), not the process's. A grown stand's automatic arrays        !
         !      overflow 4 MB after a few decades of regrowth and the run dies with SIGSEGV (the        !
         !      stack-ceiling note in CMakeLists.txt). Say so at open when neither variable sets it. ---!
         if (cfg%n_threads > 1_ik .and. .not. (env_set('OMP_STACKSIZE') .or. env_set('KMP_STACKSIZE'))) &
            write(*,'(2a)') ' WARNING: region: OMP_STACKSIZE is not set, and a grown stand overflows a ', &
                            'worker thread''s default stack (4 MB with ifx); set it, e.g. OMP_STACKSIZE=512M'
         !----- The cells: the box's valid cells with enough land, in row-major order, less those the  !
         !      cell mask drops; each one's soil and site index from the [maps] files. ----------------!
         grid = era5land_grid(cfg%forcing)
         if (.not. select_cells(cfg, grid, cells, reg%verbose)) return
         call read_map_inputs(cfg%maps, grid, cells, inputs, msg)
         if (len(msg) > 0) then
            write(*,'(2a)') ' region: ', msg
            return
         end if
         n = cells%ncell
         allocate(ids(n))
         ids = cells%row * cells%nlon + cells%col        ! each cell's row-major index on the global grid
      else
         !----- One polygon at the [site]: its soil column and height allometry's site index from   !
         !      the [maps] files at the forcing cell it reads, before anything is built on them. -------!
         if (soil_maps_on(cfg%maps) .or. len_trim(cfg%maps%chave_e_map) > 0) then
            call site_map_inputs(reg%cfg, msg)
            if (len(msg) > 0) then
               write(*,'(2a)') ' maps  : ', msg
               return
            end if
            if (reg%verbose) write(*,'(a)') ' maps  : soil and site index from the [maps] files'
         end if
         n = 1_ik
         allocate(ids(1), source=0_ik)
      end if

      !----- A restart (init_mode = 2): each polygon's row of the state file, by its id. A file may  !
      !      hold more polygons than the run: a region whose cell failed resumes with that cell      !
      !      masked out. A one-polygon run also reads a state file without a polygon axis, as every   !
      !      site run wrote one before runs became regions. -----------------------------------------!
      allocate(slot(n), source=0_ik)
      if (cfg%init%init_mode == INIT_RESTART) then
         call restart_file_slots(trim(cfg%init%restart_file), ids, slot, msg)
         if (len(msg) > 0) then
            write(*,'(4a)') ' init  : restart (mode 2): the state file ', trim(cfg%init%restart_file), ' ', msg
            return
         end if
      end if

      !----- The forcing source: a region's cells, or the cell or grid nearest the [site]. ----------!
      if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) then
         if (region) then
            call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,              &
                          run_end=cfg%end_time, cells=cells)
         else
            call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,              &
                          run_end=cfg%end_time)
         end if
         if (st /= MET_OK) return
         if (reg%verbose) then
            if (cfg%forcing%format == MET_BACKEND_ED_ERA5LAND) then
               write(*,'(3a)') ' force : met forcing ON (ED_ERA5land archive ', trim(cfg%forcing%data_path), ')'
            else
               write(*,'(3a)') ' force : met forcing ON (', trim(cfg%forcing%path), ')'
            end if
         end if
      end if
      if (.not. region) call site_cell(reg%met_src, cfg, cells)

      !----- Each polygon at its location, from its initial community, with its own configuration: !
      !      the run's, with a region cell's own soil and site index. --------------------------------!
      allocate(reg%poly(n), reg%run(n), reg%out_bufs(n))
      reg%now = cfg%start_time
      do p = 1_ik, n
         associate (poly => reg%poly(p), run => reg%run(p))
            run%cell = p
            pcfg = cfg
            if (region) then
               write(run%label,'(a,i0,a,f0.2,a,f0.2,a)') 'polygon ', ids(p), ' (', cells%lat(p), ', ',   &
                                                       cells%lon(p), ')'
               call apply_map_inputs(inputs, p, pcfg)
               call validate_config(pcfg)
            end if
            fast_found = .false.
            !----- Its sites: [[sites]]'s, each with its configuration and share, or one site. ----------!
            if (allocated(site_cfg)) deallocate(site_cfg, site_area)
            if (allocated(pcfg%sites)) call polygon_sites(pcfg, site_cfg, site_area)
            call init_polygon(poly, pcfg, polygon_name(run), reg%verbose .and. n == 1_ik, slot(p), reg%now,     &
                              fast_found, started, site_cfg=site_cfg, site_area=site_area)
            if (.not. started) return
            !----- The polygon's id and place, on the stand the initial community built. -------------!
            poly%id            = ids(p)
            poly%latitude_deg  = cells%lat(p)
            poly%longitude_deg = cells%lon(p)
            poly%elevation_m   = cells%elevation(p)
            call set_site_terrain(poly, pcfg)

            !----- The fast context and reservoirs, the forcing cursor, snow, soil-carbon spin-up. One  !
            !      polygon gives its patch loop the run's threads; a region gives them to its polygons. -!
            call polygon_prepare(cfg, reg%met_src, poly, run, keep_fast_state=fast_found,              &
                                 keep_soil_carbon=cfg%init%init_mode == INIT_RESTART,                     &
                                 patch_threads=merge(cfg%n_threads, 1_ik, n == 1_ik),               &
                                 verbose=reg%verbose .and. n == 1_ik)
         end associate
      end do

      reg%step_days = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
      allocate(reg%cost(n), source=0.0_wp)

      if (reg%verbose) then
         write(*,'(a)') '==================== MEDS run ===================='
         if (region) write(*,'(a,4(1x,f0.3),a,i0,a,f0.2,a)') ' region: box_nwse =', cfg%region%box_nwse, ', ', &
               n, ' polygons (land fraction >= ', cfg%region%land_fraction_min, ')'
         write(*,'(5a)') ' run   : ', time_to_string(cfg%start_time), ' -> ', time_to_string(cfg%end_time)
         write(*,'(a,f0.2,a,i0)') ' span  : ', years_between(reg%now, cfg%end_time), ' yr  steps/yr~',    &
               max(1_ik, nint(yr_day / real(reg%step_days, wp), ik))
      end if

      !----- The STATE stream and the run's PFT provenance table. ------------------------------------!
      if (cfg%state%write_state) then
         call ensure_output_dir(trim(cfg%state%output_dir))
         call write_pft_params_csv(cfg, trim(cfg%state%output_dir)//'/'//trim(cfg%state%output_prefix)// &
                                   '_pft_parameters.csv')
      end if

      !----- DIAGNOSTIC output ([output].enabled): every configured variable and tier, behind the     !
      !      polygon axis. Each polygon writes its own soil layers, cohorts and patches and leaves the  !
      !      rest of each axis empty; the files carry each polygon's own layer depths and inputs. -----!
      if (cfg%output%enabled) then
         call ensure_output_dir(trim(cfg%output%dir))
         call open_output_files(cfg, reg%out_files, reg%verbose, cells, ids, reg%poly(:)%site%n)
         do p = 1_ik, n
            call attach_output(cfg, reg%out_files, reg%out_bufs(p), reg%poly(p), reg%run(p))
         end do
         if (reg%verbose) write(*,'(a)') ' output: diagnostic aggregation ON ([output])'
      end if

      !----- The PARAMETER RECORD: every setting the loader read, from which file, whether it was set  !
      !      there or defaulted, and the value used -- beside the output and beside the state,         !
      !      whichever the run writes. A caller that set a key finds it here or learns that the key was !
      !      never read. ---------------------------------------------------------------------------------!
      if (cfg%output%enabled)                                                                      &
         call toml_write_record(trim(cfg%output%dir)//'/'//trim(cfg%output%prefix)//'_parameters.csv')
      if (cfg%state%write_state)                                                                   &
         call toml_write_record(trim(cfg%state%output_dir)//'/'//trim(cfg%state%output_prefix)//        &
                                '_parameters.csv')

      if (reg%verbose .and. n == 1_ik) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         call print_summary(reg%poly(1), 'start')
      end if
      end associate
      ok = .true.
   end subroutine region_open

   !----- A polygon's sites from its configuration's [[sites]]: each one's configuration (the        !
   !      polygon's, with the site's own ground when it has one) and its share of the polygon, the   !
   !      shares made to sum to exactly 1. ----------------------------------------------------------!
   subroutine polygon_sites(pcfg, site_cfg, site_area)
      type(meds_config_t),              intent(in)  :: pcfg
      type(meds_config_t), allocatable, intent(out) :: site_cfg(:)
      real(wp),            allocatable, intent(out) :: site_area(:)
      integer(ik) :: is
      allocate(site_cfg(size(pcfg%sites)))
      do is = 1_ik, size(pcfg%sites, kind=ik)
         site_cfg(is) = site_config(pcfg, is)
      end do
      site_area = pcfg%sites%area_fraction / sum(pcfg%sites%area_fraction)
   end subroutine polygon_sites

   !----- Each site's terrain: the polygon's elevation, flat, unless its [[sites]] table gives its own. !
   subroutine set_site_terrain(poly, pcfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: pcfg
      integer(ik) :: is
      poly%site%surface(1:poly%site%n)%elevation_m = poly%elevation_m
      if (.not. allocated(pcfg%sites)) return
      do is = 1_ik, poly%site%n
         associate (spec => pcfg%sites(is), surf => poly%site%surface(is))
            if (spec%has_elevation) surf%elevation_m = spec%elevation
            surf%slope_deg  = spec%slope
            surf%aspect_deg = spec%aspect
         end associate
      end do
   end subroutine set_site_terrain

   !----- A site run's [maps] inputs, at the forcing cell its [site] location reads (the archive's  !
   !      own choice, select_site_cell), into cfg, which is then validated again. msg is '' on    !
   !      success, else why the inputs could not be had. ---------------------------------------!
   subroutine site_map_inputs(cfg, msg)
      type(meds_config_t),           intent(inout) :: cfg
      character(len=:), allocatable, intent(out)   :: msg
      type(met_grid_t)   :: grid
      type(met_cells_t)  :: cell
      type(map_inputs_t) :: inp
      real(wp)    :: distance_km
      integer(ik) :: st
      grid = era5land_grid(cfg%forcing)
      call select_site_cell(grid, cfg%forcing%latitude, cfg%forcing%longitude, cfg%forcing%max_distance_km, &
                            cell, distance_km, st)
      if (st /= GRID_OK) then
         if (st == GRID_ERR_NO_CELL) then
            msg = 'no valid forcing cell near the site for the [maps] maps'
         else
            msg = 'cannot read the forcing''s static file '//trim(grid%static_file)//' for the [maps] maps'
         end if
         return
      end if
      call read_map_inputs(cfg%maps, grid, cell, inp, msg)
      if (len(msg) > 0) return
      call apply_map_inputs(inp, 1_ik, cfg)
      call validate_config(cfg)
   end subroutine site_map_inputs

   !----- A region's cells: the box's valid cells with enough land, in row-major order, less those   !
   !      the cell mask drops. .false. (after saying why) when none is left or a file cannot be read. !
   logical function select_cells(cfg, grid, cells, verbose) result(ok)
      type(meds_config_t), intent(in)  :: cfg
      type(met_grid_t),    intent(in)  :: grid
      type(met_cells_t),   intent(out) :: cells
      logical,             intent(in)  :: verbose
      integer(ik) :: st, n, j
      ok = .false.
      call select_box_cells(grid, cfg%region%box_nwse, cells, st, cfg%region%land_fraction_min)
      if (st == GRID_ERR_NO_CELL) then
         write(*,'(a,4(f0.3,1x),a,f0.2)') ' region: no valid cell in box_nwse = ', cfg%region%box_nwse, &
               'with land fraction >= ', cfg%region%land_fraction_min
         return
      else if (st /= GRID_OK) then
         write(*,'(2a)') ' region: cannot read valid, elevation and land_fraction from ', trim(grid%static_file)
         return
      end if
      if (len_trim(cfg%region%cell_mask) > 0) then
         n = cells%ncell
         call apply_cell_mask(cells, trim(cfg%region%cell_mask), trim(cfg%region%cell_mask_var),       &
                              cfg%region%cell_mask_min, grid, j, st)
         if (st /= MAP_OK) then
            write(*,'(4a)') ' region: the cell mask ', trim(cfg%region%cell_mask), ' ', map_status_message(st)
            return
         end if
         if (verbose) write(*,'(a,i0,a,i0,3a)') ' region: the cell mask keeps ', cells%ncell, ' of ', n, &
               ' cells (', trim(cfg%region%cell_mask_var), ')'
         if (cells%ncell == 0_ik) then
            write(*,'(2a)') ' region: the cell mask keeps no cell of the box: ', trim(cfg%region%cell_mask)
            return
         end if
      end if
      ok = .true.
   end function select_cells

   !----- A one-polygon run's cell: the [site] location, which the polygon runs at, and the grid      !
   !      indices of the archive cell the forcing reads (0 on the MEDS forcing file, whose grids are    !
   !      not on a global grid). ---------------------------------------------------------------------!
   subroutine site_cell(met_src, cfg, cells)
      type(met_source_t),  intent(in)  :: met_src
      type(meds_config_t), intent(in)  :: cfg
      type(met_cells_t),   intent(out) :: cells
      cells%ncell = 1_ik
      allocate(cells%row(1), cells%col(1), source=0_ik)
      allocate(cells%lat(1), source=cfg%forcing%latitude)
      allocate(cells%lon(1), source=cfg%forcing%longitude)
      allocate(cells%elevation(1), source=cfg%forcing%elevation)
      if (met_src%cells%ncell == 1_ik) then
         cells%nlat = met_src%cells%nlat
         cells%nlon = met_src%cells%nlon
         cells%row  = met_src%cells%row
         cells%col  = met_src%cells%col
      end if
   end subroutine site_cell

   !----- How messages name a polygon: its label, or "the site" for a one-polygon run. -------------!
   function polygon_name(run) result(name)
      type(polygon_run_t), intent(in) :: run
      character(len=:), allocatable :: name
      if (len_trim(run%label) > 0) then
         name = trim(run%label)
      else
         name = 'the site'
      end if
   end function polygon_name

   !----- Is the environment variable set? ---------------------------------------------------------!
   logical function env_set(name) result(yes)
      character(len=*), intent(in) :: name
      integer :: n, st
      call get_environment_variable(name, length=n, status=st)
      yes = st == 0 .and. n > 0
   end function env_set

   !---------------------------------------------------------------------------------------!
   ! region_step -- ONE slow step of every polygon: load the step's forcing (a no-op unless a new  !
   ! month starts), step the polygons, and after a month's last step run the I/O phase. The C API  !
   ! steps a run this way; region_step_month takes a month of these steps at once.                !
   !---------------------------------------------------------------------------------------!
   subroutine region_step(reg, status)
      type(meds_region_t), intent(inout) :: reg
      integer(ik),         intent(out)   :: status
      type(meds_time_t) :: step_start(1), step_end(1)
      logical           :: new_month(1), new_year(1)
      status = DRIVER_OK
      if (region_done(reg)) then
         status = DRIVER_FINISHED
         return
      end if
      step_start(1) = reg%now
      step_end(1)   = time_advance_days(reg%now, reg%step_days)
      new_year(1)   = step_end(1)%year /= step_start(1)%year
      new_month(1)  = new_year(1) .or. step_end(1)%month /= step_start(1)%month
      call step_polygons(reg, step_start, step_end, new_month, new_year, status)
   end subroutine region_step

   !---------------------------------------------------------------------------------------!
   ! region_step_month -- one month of every polygon: load the month's forcing for every cell, !
   ! advance each polygon through it, then write the queued output. The month runs from the   !
   ! current date to the next month boundary, or to end_time.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine region_step_month(reg, status)
      type(meds_region_t), intent(inout) :: reg
      integer(ik),         intent(out)   :: status
      type(meds_time_t), allocatable :: step_start(:), step_end(:)
      logical,           allocatable :: new_month(:), new_year(:)
      type(meds_time_t) :: clk
      integer(ik)       :: nstep, s

      status = DRIVER_OK
      if (region_done(reg)) then
         status = DRIVER_FINISHED
         return
      end if
      !----- The month's steps, listed once. A step's forcing lies in the month its start date      !
      !      belongs to (daily steps from midnight, and validate_config puts the recycle seam on a   !
      !      month boundary), so one prefetch at the month's start loads every step's forcing. ------!
      nstep = 0_ik
      clk = reg%now
      do
         clk = time_advance_days(clk, reg%step_days)
         nstep = nstep + 1_ik
         if (clk%year /= reg%now%year .or. clk%month /= reg%now%month) exit
         if (.not. time_lt(clk, reg%cfg%end_time)) exit
      end do
      allocate(step_start(nstep), step_end(nstep), new_month(nstep), new_year(nstep))
      clk = reg%now
      do s = 1_ik, nstep
         step_start(s) = clk
         clk = time_advance_days(clk, reg%step_days)
         step_end(s) = clk
         new_year(s)  = step_end(s)%year /= step_start(s)%year
         new_month(s) = new_year(s) .or. step_end(s)%month /= step_start(s)%month
      end do
      call step_polygons(reg, step_start, step_end, new_month, new_year, status)
   end subroutine region_step_month

   !---------------------------------------------------------------------------------------!
   ! step_polygons -- advance every polygon through the listed steps, then close what they end. !
   !                                                                                          !
   ! The forcing the steps read is loaded first. The polygons run side by side on               !
   ! [run].n_threads threads, costliest first; they share only what they read (the config, the  !
   ! forcing, the output layout), and each writes its own state and output buffers, so the      !
   ! results do not depend on the thread count. One polygon runs on the calling thread and gives !
   ! its patch loop the threads instead. A polygon that fails is reported and stops; the others   !
   ! finish, the output is written, and the run moves on either way. The status returned is the  !
   ! first failure's.                                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine step_polygons(reg, step_start, step_end, new_month, new_year, status)
      type(meds_region_t), intent(inout) :: reg
      type(meds_time_t),   intent(in)    :: step_start(:), step_end(:)
      logical,             intent(in)    :: new_month(:), new_year(:)
      integer(ik),         intent(out)   :: status
      logical,     allocatable :: was_running(:)
      integer(ik), allocatable :: order(:), failed_step(:)
      integer(ik) :: nstep, s, p, k, st, npoly
      integer(ik) :: tick0, tick1, tick_rate

      status = DRIVER_OK
      associate (cfg => reg%cfg)
      npoly = size(reg%poly, kind=ik)
      nstep = size(step_start, kind=ik)
      if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) then
         call system_clock(tick0, tick_rate)
         call met_prefetch(reg%met_src, step_start(1))
         call system_clock(tick1)
         reg%year_read_seconds = reg%year_read_seconds + real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)
      end if

      !----- The compute phase: no file is touched here; each polygon writes only its own state  !
      !      and buffers. ---------------------------------------------------------------------------!
      allocate(order(npoly), failed_step(npoly), was_running(npoly))
      call ascending_order(-reg%cost, npoly, order)
      was_running = reg%run(:)%status == DRIVER_OK
      failed_step = 0_ik
      !$omp parallel do schedule(dynamic, 1) num_threads(cfg%n_threads) default(shared) if(npoly > 1)   &
      !$omp    private(k, p, s, st, tick0, tick1, tick_rate)
      do k = 1_ik, npoly
         p = order(k)
         if (.not. was_running(p)) cycle
         call system_clock(tick0, tick_rate)
         do s = 1_ik, nstep
            call polygon_step(cfg, reg%met_src, reg%out_files, reg%out_bufs(p), reg%poly(p),         &
                              reg%run(p),                                                              &
                              step_start(s), step_end(s), reg%step_days, new_month(s), new_year(s), st)
            if (st /= DRIVER_OK) then
               failed_step(p) = s
               exit
            end if
         end do
         call system_clock(tick1)
         reg%cost(p) = real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)
      end do
      !$omp end parallel do

      do p = 1_ik, npoly
         if (failed_step(p) == 0_ik) cycle
         if (npoly > 1_ik) write(*,'(4a)') ' region: ', trim(reg%run(p)%label), ' failed on ',       &
                                           time_to_string(step_end(failed_step(p)))
         if (status == DRIVER_OK) status = reg%run(p)%status
      end do
      reg%now   = step_end(nstep)
      reg%istep = reg%istep + nstep

      !----- A year's end: the summary line(s), and the forcing-read accounting. -------------------!
      if (new_year(nstep)) then
         reg%iyear = reg%iyear + 1_ik
         if (reg%verbose) call year_summary(reg)
         reg%read_seconds = reg%read_seconds + reg%year_read_seconds
         reg%year_read_seconds = 0.0_wp
         reg%year_loads0 = reg%met_src%n_loads
      end if

      !----- The I/O phase: the records the steps closed, after a month's last step, a failure or the  !
      !      run's last step, and at a year's end every interval_years the state file. -----------------!
      if (cfg%output%enabled .and. (new_month(nstep) .or. status /= DRIVER_OK .or. region_done(reg)))  &
         call output_serialize(reg%out_files, reg%out_bufs)
      if (new_year(nstep) .and. cfg%state%write_state .and. mod(reg%iyear, cfg%state%interval_years) == 0_ik) &
         call region_write_state(reg)
      end associate
   end subroutine step_polygons

   !----- The state file: every polygon's row, at the current date (meds_init_restart). The cohort, patch and   !
   !      site axes are the largest polygon's. The polygons share the calendar, so they owe the same  !
   !      boundary restructuring. --------------------------------------------------------------------!
   subroutine region_write_state(reg)
      type(meds_region_t), intent(in) :: reg
      type(restart_file_t)       :: sf
      integer(ik), allocatable :: ids(:)
      real(wp),    allocatable :: lat(:), lon(:)
      integer(ik)              :: p, n_cohort, n_patch, n_site
      n_cohort = 0_ik
      n_patch = 0_ik
      n_site = 0_ik
      allocate(ids(size(reg%poly)), lat(size(reg%poly)), lon(size(reg%poly)))
      do p = 1_ik, size(reg%poly, kind=ik)
         ids(p)   = reg%poly(p)%id
         lat(p)   = reg%poly(p)%latitude_deg
         lon(p)   = reg%poly(p)%longitude_deg
         n_cohort = max(n_cohort, reg%poly(p)%cohort%n)
         n_patch  = max(n_patch, reg%poly(p)%patch%n)
         n_site   = max(n_site, reg%poly(p)%site%n)
      end do
      call restart_file_open(sf, trim(reg%cfg%state%output_dir), trim(reg%cfg%state%output_prefix), reg%now, &
                           reg%cfg%pft%n, int(size(reg%poly(1)%cohort%growth_hist, 1), ik), n_cohort,       &
                           n_patch, n_site, polygon_id=ids, polygon_lat=lat, polygon_lon=lon,              &
                           restructure_pending=reg%poly(1)%restructure_pending,                            &
                           restructure_new_year=reg%poly(1)%restructure_new_year)
      do p = 1_ik, size(reg%poly, kind=ik)
         call restart_file_write(sf, reg%poly(p), reg%now, p)
      end do
      call restart_file_close(sf)
   end subroutine region_write_state

   !----- A year's summary: one polygon's stand every 5 years, or a line for the polygons' mean and  !
   !      range of AGB and LAI every year, with what the forcing reads cost. ------------------------!
   subroutine year_summary(reg)
      type(meds_region_t), intent(in) :: reg
      character(len=19) :: when
      real(wp)    :: agb(size(reg%poly)), lai(size(reg%poly))
      integer(ik) :: p
      when = time_to_string(reg%now)
      if (size(reg%poly) == 1) then
         if (mod(reg%iyear, 5_ik) == 0_ik .or. reg%iyear == 1_ik) call print_summary(reg%poly(1), when(1:10))
         return
      end if
      do p = 1_ik, size(reg%poly, kind=ik)
         agb(p) = total_agb(reg%poly(p))
         lai(p) = total_lai(reg%poly(p))
      end do
      write(*,'(a,a10,a,3(1x,f0.3),a,3(1x,f0.3))') ' ', when(1:10), '  AGB kgC/m2 min/mean/max:',     &
            minval(agb), sum(agb) / size(agb), maxval(agb), '   LAI:', minval(lai),                  &
            sum(lai) / size(lai), maxval(lai)
      write(*,'(a,f0.2,a,i0,a)') '             forcing: ', reg%year_read_seconds, ' s loading, ',        &
            reg%met_src%n_loads - reg%year_loads0, ' archive month(s) read from file'
   end subroutine year_summary

   !---------------------------------------------------------------------------------------!
   ! region_finalize -- the terminal state file, the end-of-run report (one polygon's stand and   !
   ! conservation reports, or a table of the polygons), close every file and the forcing source.  !
   ! DRIVER_ERR_AREA if any polygon's patch areas do not sum to 1. Once a polygon has failed, the   !
   ! run writes nothing more: no state file and no partial periods, so a closed chunk file is      !
   ! never made again from scratch.                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine region_finalize(reg, status)
      type(meds_region_t),   intent(inout) :: reg
      integer(ik), optional, intent(out)   :: status
      integer(ik) :: p, st, n_bad
      real(wp)    :: a1
      logical     :: failed

      st = DRIVER_OK
      n_bad = 0_ik
      failed = any(reg%run(:)%status /= DRIVER_OK)
      !----- The true terminal state, so a restart resumes exactly here. -----------------------------!
      if (reg%cfg%state%write_state .and. .not. failed) call region_write_state(reg)
      if (size(reg%poly) == 1) then
         associate (poly => reg%poly(1), run => reg%run(1))
            if (reg%verbose) call print_summary(poly, 'final')
            a1 = total_area(poly)
            if (reg%verbose) then
               write(*,'(a)') '-----------------------------------------------------------------------------'
               write(*,'(a,f12.9,a,f12.9)') ' site area start=', run%area_start, '  end=', a1
            end if
            if (abs(a1 - 1.0_wp) > AREA_TOL) st = DRIVER_ERR_AREA
            if (reg%verbose) call polygon_report(reg%cfg, run)
         end associate
      else
         if (reg%verbose) then
            write(*,'(a)') '-----------------------------------------------------------------------------'
            write(*,'(a)') '  polygon      lat      lon  AGB kgC/m2    LAI  cohorts  area error  energy/water fails'
         end if
         do p = 1_ik, size(reg%poly, kind=ik)
            associate (poly => reg%poly(p), run => reg%run(p))
               a1 = total_area(poly)
               if (abs(a1 - 1.0_wp) > AREA_TOL) then
                  st = DRIVER_ERR_AREA
                  n_bad = n_bad + 1_ik
               end if
               if (reg%verbose)                                                                    &
                  write(*,'(i9,2f9.2,f12.4,f7.3,i9,es12.2,2(1x,i0))') poly%id,                      &
                        poly%latitude_deg, poly%longitude_deg, total_agb(poly),                     &
                        total_lai(poly), poly%cohort%n, a1 - 1.0_wp,                          &
                        run%energy_budget%n_fail, run%water_budget%n_fail
            end associate
         end do
         if (n_bad > 0_ik) write(*,'(a,i0,a)') ' ERROR: ', n_bad, ' polygon(s) did not conserve area'
         if (reg%verbose) write(*,'(a,f0.2,a,i0,a)') ' forcing: ',                                  &
               reg%read_seconds + reg%year_read_seconds, ' s loading in all, ', reg%met_src%n_loads,    &
               ' archive month(s) read from file'
      end if

      if (reg%cfg%output%enabled) call output_close(reg%out_files, reg%out_bufs, flush_partial=.not. failed)
      if (reg%cfg%fast_biophysics_on .and. reg%cfg%forcing%forcing_on) call met_close(reg%met_src)
      if (present(status)) status = st
   end subroutine region_finalize

   !---------------------------------------------------------------------------------------!
   ! region_free -- release a run, finalized or not: close its output files (a period still    !
   ! open is not written; region_finalize writes it) and its forcing source, free its polygons, !
   ! and leave the variable as a run never opened, which region_open takes again. Separate from  !
   ! finalize, so a caller can still read the final state (the getters the C-API exposes) after   !
   ! the files are closed.                                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine region_free(reg)
      type(meds_region_t), intent(inout) :: reg
      type(meds_region_t) :: never_opened     ! default-initialised, not built by a constructor (nvfortran)
      integer(ik) :: p
      if (allocated(reg%out_bufs)) call output_close(reg%out_files, reg%out_bufs, flush_partial=.false.)
      call met_close(reg%met_src)
      if (allocated(reg%poly)) then
         do p = 1_ik, size(reg%poly, kind=ik)
            call free_polygon(reg%poly(p))
         end do
      end if
      reg = never_opened
   end subroutine region_free

end module meds_region
