! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_region -- a region run: every selected forcing cell of a contiguous box as its own        !
! polygon, in one process (docs/science/order_of_processes.md §1; the configuration is           !
! docs/configuration.md §Regional runs).                                                         !
!                                                                                          !
!   type(meds_region_t) :: reg                                                                   !
!   call region_open('meds_config_main.toml', reg, ok)                                            !
!   do while (.not. region_done(reg)) ; call region_step_month(reg, status) ; end do               !
!   call region_finalize(reg, status) ; call region_free(reg)                                      !
!                                                                                          !
! The unit of work is a calendar month. The forcing source loads the month for every cell once;   !
! then each polygon runs through the whole month with polygon_step, the same step a site run       !
! takes, reading its own cell of the shared month buffer; then the I/O phase writes every          !
! polygon's queued records into the region files, one hyperslab per variable per record. Nothing  !
! inside the polygon loop touches a file, so the polygons of a month run side by side on           !
! [run].n_threads threads, each with a single-threaded patch loop.                                 !
!                                                                                          !
! Each polygon sits at its cell's centre, at the cell's orography, in UTC, and starts from bare    !
! ground or from its row of a region state file. The region's files hold every variable and tier,  !
! each polygon's slice behind a `polygon` axis, cohort and patch axes padded to the largest        !
! polygon's. With [state].write_state the region writes one state file, every polygon in it        !
! (meds_io), every interval_years and at the end. Once a polygon has failed, the end of the run    !
! writes nothing more.                                                                             !
!==========================================================================================!
module meds_region
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : day_sec
   use meds_numerics,               only : ascending_order
   use meds_config,                 only : meds_config_t, validate_config, INIT_RESTART
   use meds_config_io,              only : load_meds_config
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days, time_to_string, &
                                           years_between
   use meds_site_state_types,       only : site_free
   use meds_io,                     only : state_file_t, state_file_open, state_file_put, state_file_close, &
                                           state_file_polygons, io_read_state, state_file_refusal
   use meds_vegetation_dynamics,    only : reacclimate_plant_traits
   use meds_init,                   only : init_bare_ground, seed_bare_ground
   use meds_forcing_types,          only : met_source_t, met_cells_t
   use meds_forcing_config,         only : MET_PATH_LEN
   use meds_met_driver,             only : met_open, met_close, met_prefetch, MET_OK
   use meds_surface_maps,           only : apply_cell_mask, map_status_message, MAP_OK,           &
                                           map_inputs_t, read_map_inputs, apply_map_inputs
   use meds_era5land_reader,        only : era5land_select_box, era5land_static_file, ERA_OK, ERA_ERR_NO_CELL
   use meds_diagnostic_reduce,      only : total_area, total_agb, total_lai, count_cohorts
   use meds_polygon,                only : meds_polygon_t, polygon_prepare, polygon_step,        &
                                           open_output_files, attach_output, ensure_output_dir,   &
                                           DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_AREA, N_PATCH_INIT
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_manager,         only : output_serialize, output_close
   implicit none
   private

   public :: meds_region_t, region_open, region_step_month, region_done, region_finalize, region_free

   type :: meds_region_t
      type(meds_config_t)   :: cfg
      type(met_source_t)    :: met_src          !< one reader for every cell of the region
      type(output_files_t)  :: out_files           !< the region files (built only if output.enabled)
      type(meds_polygon_t), allocatable :: poly(:)
      !----- Each polygon's share of the region files, one per polygon. A contiguous array of its own,  !
      !      so the I/O phase passes it whole: a section poly(:)%bufs of a type with allocatable        !
      !      components makes gfortran copy it through a temporary whose copy-out dangles them. -------!
      type(output_buffers_t), allocatable :: out_bufs(:)
      type(meds_time_t)     :: now
      integer(ik)           :: step_days = 1_ik
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
   ! region_open -- config -> the box's cells -> the shared forcing source -> one polygon per    !
   ! cell -> the region files.                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine region_open(path, reg, ok, verbose)
      character(len=*),    intent(in)    :: path
      type(meds_region_t), intent(inout) :: reg
      logical,             intent(out)   :: ok
      logical, optional,   intent(in)    :: verbose
      type(met_cells_t)          :: cells
      character(len=MET_PATH_LEN) :: static
      integer(ik)                :: st, p, n, j
      integer(ik), allocatable   :: ids(:)
      type(map_inputs_t)     :: inputs
      type(meds_config_t)        :: pcfg
      character(len=:), allocatable :: msg
      type(meds_time_t)          :: restart_time
      integer(ik), allocatable   :: file_ids(:), slot(:)
      logical                    :: found, fast_found

      ok = .false.
      if (present(verbose)) reg%verbose = verbose

      call load_meds_config(trim(path), reg%cfg)
      if (reg%verbose) write(*,'(2a)') ' config: ', trim(path)
      associate (cfg => reg%cfg)

      !----- A region steps each polygon on an OpenMP worker thread, whose stack is OMP_STACKSIZE (4 MB  !
      !      by default with ifx), not the process's. A grown stand's automatic arrays overflow 4 MB after !
      !      a few decades of regrowth and the run dies with SIGSEGV (the stack-ceiling note in            !
      !      CMakeLists.txt). Say so at open when neither variable sets it. ------------------------------!
      if (cfg%n_threads > 1_ik .and. .not. (env_set('OMP_STACKSIZE') .or. env_set('KMP_STACKSIZE')))  &
         write(*,'(2a)') ' WARNING: region: OMP_STACKSIZE is not set, and a grown stand overflows a worker ', &
                         'thread''s default stack (4 MB with ifx); set it, e.g. OMP_STACKSIZE=512M'

      !----- The cells: the box's valid cells with enough land, in row-major order. --------------!
      static = era5land_static_file(cfg%forcing)
      call era5land_select_box(static, cfg%region%box_nwse, cells, st, cfg%region%land_fraction_min)
      if (st == ERA_ERR_NO_CELL) then
         write(*,'(a,4(f0.3,1x),a,f0.2)') ' region: no valid cell in box_nwse = ', cfg%region%box_nwse, &
               'with land fraction >= ', cfg%region%land_fraction_min
         return
      else if (st /= ERA_OK) then
         write(*,'(2a)') ' region: cannot read valid, elevation and land_fraction from ', trim(static)
         return
      end if
      !----- The cell mask: of those cells, the ones the mask keeps. ------------------------------!
      if (len_trim(cfg%region%cell_mask) > 0) then
         n = cells%ncell
         call apply_cell_mask(cells, trim(cfg%region%cell_mask), trim(cfg%region%cell_mask_var),       &
                              cfg%region%cell_mask_min, static, j, st)
         if (st /= MAP_OK) then
            write(*,'(4a)') ' region: the cell mask ', trim(cfg%region%cell_mask), ' ', map_status_message(st)
            return
         end if
         if (reg%verbose) write(*,'(a,i0,a,i0,3a)') ' region: the cell mask keeps ', cells%ncell, ' of ', n, &
               ' cells (', trim(cfg%region%cell_mask_var), ')'
         if (cells%ncell == 0_ik) then
            write(*,'(2a)') ' region: the cell mask keeps no cell of the box: ', trim(cfg%region%cell_mask)
            return
         end if
      end if

      !----- The polygon ids, each cell's row-major index on the global grid. -----------------------!
      n = cells%ncell
      allocate(ids(n))
      ids = cells%row * cells%nlon + cells%col

      !----- Each cell's soil and site index, from the [maps] files, read once for every cell. ----!
      call read_map_inputs(cfg%maps, static, cells, inputs, msg)
      if (len(msg) > 0) then
         write(*,'(2a)') ' region: ', msg
         return
      end if

      !----- A restart (init_mode = 2): each polygon's row of the region state file, by its id. The  !
      !      file must hold exactly this region's polygons. -----------------------------------------!
      if (cfg%init_mode == INIT_RESTART) then
         msg = state_file_refusal(trim(cfg%init_restart_file), region=.true.)
         if (len(msg) > 0) then
            write(*,'(4a)') ' region: the state file ', trim(cfg%init_restart_file), ' ', msg
            return
         end if
         call state_file_polygons(trim(cfg%init_restart_file), file_ids, found)
         if (.not. found) then
            write(*,'(3a)') ' region: the state file ', trim(cfg%init_restart_file), ' names no polygons'
            return
         end if
         allocate(slot(n))
         do p = 1_ik, n
            slot(p) = int(findloc(file_ids, ids(p), dim=1), ik)
         end do
         if (size(file_ids) /= n .or. any(slot == 0_ik)) then
            write(*,'(3a)') ' region: the polygons of ', trim(cfg%init_restart_file), ' are not this region''s'
            return
         end if
      end if

      call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,                  &
                    run_end=cfg%end_time, cells=cells)
      if (st /= MET_OK) return
      if (reg%verbose) write(*,'(3a)') ' force : met forcing ON (ED_ERA5land archive ', trim(cfg%forcing%data_path), ')'

      !----- One polygon per cell, each at its cell centre and orography, in UTC, from bare ground  !
      !      or its row of the restart file, with its own configuration: the run's, with the cell's  !
      !      soil and site index. ---------------------------------------------------------------------!
      allocate(reg%poly(n), reg%out_bufs(n))
      reg%now = cfg%start_time
      do p = 1_ik, n
         associate (poly => reg%poly(p))
            poly%cell = p
            poly%id   = ids(p)
            write(poly%label,'(a,i0,a,f0.2,a,f0.2,a)') 'polygon ', poly%id, ' (', cells%lat(p), ', ', &
                                                     cells%lon(p), ')'
            pcfg = cfg
            call apply_map_inputs(inputs, p, pcfg)
            call validate_config(pcfg)
            fast_found = .false.
            if (cfg%init_mode == INIT_RESTART) then
               call io_read_state(poly%site, pcfg, trim(cfg%init_restart_file), restart_time, found,   &
                                  fast_found=fast_found, restructure_pending=poly%restructure_pending, &
                                  restructure_new_year=poly%restructure_new_year, polygon=slot(p))
               if (.not. found) then
                  write(*,'(5a)') ' region: cannot read ', trim(poly%label), '''s row of the state file ', &
                                  trim(cfg%init_restart_file)
                  return
               end if
               reg%now = restart_time
               if (pcfg%init_reacclimate_traits) call reacclimate_plant_traits(poly%site, pcfg)
            else
               call init_bare_ground(poly%site, pcfg, N_PATCH_INIT)
               call seed_bare_ground(poly%site, pcfg)
            end if
            call polygon_prepare(pcfg, reg%met_src, poly, cells%lat(p), cells%lon(p),              &
                                 cells%elevation(p), keep_fast_state=fast_found,                  &
                                 keep_soil_carbon=cfg%init_mode == INIT_RESTART, patch_threads=1_ik, &
                                 verbose=.false.)
         end associate
      end do

      reg%step_days = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
      if (cfg%state_write_state) call ensure_output_dir(trim(cfg%state_output_dir))
      allocate(reg%cost(n)) ; reg%cost = 0.0_wp

      !----- Output. The region files: every configured variable and tier, with the polygon axis.   !
      !      Each polygon writes its own soil layers, cohorts and patches and leaves the rest of each  !
      !      axis empty; the files carry each polygon's own layer depths and map inputs. --------------!
      if (cfg%output%enabled) then
         call ensure_output_dir(trim(cfg%output%dir))
         call open_output_files(cfg, reg%out_files, reg%verbose, cells=cells, polygon_id=ids)
         do p = 1_ik, n
            call attach_output(reg%out_files, reg%out_bufs(p), reg%poly(p))
         end do
      end if

      if (reg%verbose) then
         write(*,'(a)') '==================== MEDS region run ===================='
         write(*,'(a,4(1x,f0.3),a,i0,a,f0.2,a)') ' region: box_nwse =', cfg%region%box_nwse, ', ',    &
               n, ' polygons (land fraction >= ', cfg%region%land_fraction_min, ')'
         write(*,'(5a)') ' run   : ', time_to_string(cfg%start_time), ' -> ', time_to_string(cfg%end_time)
         write(*,'(a,f0.2,a)') ' span  : ', years_between(reg%now, cfg%end_time), ' yr'
      end if
      end associate
      ok = .true.
   end subroutine region_open

   !----- Is the environment variable set? ---------------------------------------------------------!
   logical function env_set(name) result(yes)
      character(len=*), intent(in) :: name
      integer :: n, st
      call get_environment_variable(name, length=n, status=st)
      yes = st == 0 .and. n > 0
   end function env_set

   !---------------------------------------------------------------------------------------!
   ! region_step_month -- one month of the whole region: load the month's forcing for every cell, !
   ! advance each polygon through it, then write the queued output. The month runs from the       !
   ! current date to the next month boundary, or to end_time.                                      !
   !                                                                                          !
   ! The polygons run side by side on [run].n_threads threads. They share only what they read      !
   ! (the config, the forcing month, the output layout), and each writes its own state and output  !
   ! buffers, so the region's results do not depend on the thread count. A polygon that fails is    !
   ! reported and stops; the others finish the month, the month's output is written, and the       !
   ! region moves on to the next month either way. The status returned is the first failure's.     !
   !---------------------------------------------------------------------------------------!
   subroutine region_step_month(reg, status)
      type(meds_region_t), intent(inout) :: reg
      integer(ik),         intent(out)   :: status
      type(meds_time_t), allocatable :: step_start(:), step_end(:)
      logical,           allocatable :: new_month(:), new_year(:), was_running(:)
      integer(ik),       allocatable :: order(:), failed_step(:)
      type(meds_time_t) :: clk
      integer(ik)       :: nstep, s, p, k, st, npoly
      integer(ik)       :: tick0, tick1, tick_rate

      status = DRIVER_OK
      if (region_done(reg)) then ; status = DRIVER_FINISHED ; return ; end if
      associate (cfg => reg%cfg)
      npoly = size(reg%poly, kind=ik)

      !----- The month's steps, listed once. A step's forcing lies in the month its start date      !
      !      belongs to (daily steps from midnight, and validate_config puts the recycle seam on a   !
      !      month boundary), so one prefetch at the month's start loads every step's forcing. ------!
      nstep = 0_ik ; clk = reg%now
      do
         clk = time_advance_days(clk, reg%step_days) ; nstep = nstep + 1_ik
         if (clk%year /= reg%now%year .or. clk%month /= reg%now%month) exit
         if (.not. time_lt(clk, cfg%end_time)) exit
      end do
      allocate(step_start(nstep), step_end(nstep), new_month(nstep), new_year(nstep))
      clk = reg%now
      do s = 1_ik, nstep
         step_start(s) = clk ; clk = time_advance_days(clk, reg%step_days) ; step_end(s) = clk
         new_year(s)  = step_end(s)%year /= step_start(s)%year
         new_month(s) = new_year(s) .or. step_end(s)%month /= step_start(s)%month
      end do
      call system_clock(tick0, tick_rate)
      call met_prefetch(reg%met_src, reg%now)
      call system_clock(tick1)
      reg%year_read_seconds = reg%year_read_seconds + real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)

      !----- The compute phase: every polygon through the month, costliest first. No file is      !
      !      touched here; each polygon writes only its own state and buffers. ---------------------!
      allocate(order(npoly), failed_step(npoly), was_running(npoly))
      call ascending_order(-reg%cost, npoly, order)
      was_running = reg%poly(:)%status == DRIVER_OK
      failed_step = 0_ik
      !$omp parallel do schedule(dynamic, 1) num_threads(cfg%n_threads) default(shared)            &
      !$omp    private(k, p, s, st, tick0, tick1, tick_rate)
      do k = 1_ik, npoly
         p = order(k)
         if (.not. was_running(p)) cycle
         call system_clock(tick0, tick_rate)
         do s = 1_ik, nstep
            call polygon_step(reg%met_src, reg%out_files, reg%out_bufs(p), reg%poly(p),              &
                              step_start(s), step_end(s), reg%step_days, new_month(s), new_year(s), st)
            if (st /= DRIVER_OK) then
               failed_step(p) = s ; exit
            end if
         end do
         call system_clock(tick1)
         reg%cost(p) = real(tick1 - tick0, wp) / real(max(tick_rate, 1_ik), wp)
      end do
      !$omp end parallel do

      do p = 1_ik, npoly
         if (failed_step(p) == 0_ik) cycle
         write(*,'(4a)') ' region: ', trim(reg%poly(p)%label), ' failed on ',                      &
                         time_to_string(step_end(failed_step(p)))
         if (status == DRIVER_OK) status = reg%poly(p)%status
      end do
      reg%now = step_end(nstep)

      if (new_year(nstep)) then
         if (reg%verbose) call region_summary(reg, time_to_string(reg%now))
         reg%read_seconds = reg%read_seconds + reg%year_read_seconds
         reg%year_read_seconds = 0.0_wp ; reg%year_loads0 = reg%met_src%n_loads
      end if

      !----- The I/O phase: the month's records, and at a year's end every interval_years the       !
      !      region's state file. ------------------------------------------------------------------!
      call io_phase(reg)
      if (new_year(nstep)) then
         reg%iyear = reg%iyear + 1_ik
         if (cfg%state_write_state .and. mod(reg%iyear, cfg%state_interval_years_cfg) == 0_ik)          &
            call region_write_state(reg)
      end if
      end associate
   end subroutine region_step_month

   !----- The region's state file: every polygon's row, at the current date (meds_io). The          !
   !      cohort and patch axes are the largest polygon's. The polygons share the calendar, so they   !
   !      owe the same boundary restructuring. ---------------------------------------------------------!
   subroutine region_write_state(reg)
      type(meds_region_t), intent(in) :: reg
      type(state_file_t)       :: sf
      integer(ik), allocatable :: ids(:)
      real(wp),    allocatable :: lat(:), lon(:)
      integer(ik)              :: p, n_cohort, n_patch
      n_cohort = 0_ik ; n_patch = 0_ik
      allocate(ids(size(reg%poly)), lat(size(reg%poly)), lon(size(reg%poly)))
      do p = 1_ik, size(reg%poly, kind=ik)
         ids(p)   = reg%poly(p)%id
         lat(p)   = reg%met_src%cells%lat(reg%poly(p)%cell)
         lon(p)   = reg%met_src%cells%lon(reg%poly(p)%cell)
         n_cohort = max(n_cohort, reg%poly(p)%site%cohort%n)
         n_patch  = max(n_patch, reg%poly(p)%site%patch%n)
      end do
      call state_file_open(sf, trim(reg%cfg%state_output_dir), trim(reg%cfg%state_output_prefix), reg%now, &
                           reg%cfg%pft%n, int(size(reg%poly(1)%site%cohort%growth_hist, 1), ik), n_cohort,  &
                           n_patch, reg%poly(1)%restructure_pending, reg%poly(1)%restructure_new_year,     &
                           polygon_id=ids, polygon_lat=lat, polygon_lon=lon)
      do p = 1_ik, size(reg%poly, kind=ik)
         call state_file_put(sf, reg%poly(p)%site, reg%now, p)
      end do
      call state_file_close(sf)
   end subroutine region_write_state

   !----- Write every polygon's queued records into the region files. -----------------------------!
   subroutine io_phase(reg)
      type(meds_region_t), intent(inout) :: reg
      if (.not. reg%cfg%output%enabled) return
      call output_serialize(reg%out_files, reg%out_bufs)
   end subroutine io_phase

   !----- One line for the region at a year boundary: the polygons' mean and range of AGB and LAI. -!
   subroutine region_summary(reg, when)
      type(meds_region_t), intent(in) :: reg
      character(len=*),    intent(in) :: when
      real(wp) :: agb(size(reg%poly)), lai(size(reg%poly))
      integer(ik) :: p
      do p = 1_ik, size(reg%poly, kind=ik)
         agb(p) = total_agb(reg%poly(p)%site) ; lai(p) = total_lai(reg%poly(p)%site)
      end do
      write(*,'(a,a10,a,3(1x,f0.3),a,3(1x,f0.3))') ' ', when(1:10), '  AGB kgC/m2 min/mean/max:',   &
            minval(agb), sum(agb) / size(agb), maxval(agb), '   LAI:', minval(lai),                  &
            sum(lai) / size(lai), maxval(lai)
      write(*,'(a,f0.2,a,i0,a)') '             forcing: ', reg%year_read_seconds, ' s loading, ',        &
            reg%met_src%n_loads - reg%year_loads0, ' archive month(s) read from file'
   end subroutine region_summary

   !---------------------------------------------------------------------------------------!
   ! region_finalize -- the terminal state file, the per-polygon status table (state,             !
   ! conservation), close every file and the forcing source. DRIVER_ERR_AREA if any polygon's      !
   ! patch areas do not sum to 1. Once a polygon has failed, it writes nothing more: no state      !
   ! file and no partial periods, so a closed chunk file is never made again from scratch.         !
   !---------------------------------------------------------------------------------------!
   subroutine region_finalize(reg, status)
      type(meds_region_t),   intent(inout) :: reg
      integer(ik), optional, intent(out)   :: status
      integer(ik) :: p, st, n_bad
      real(wp)    :: a1
      logical     :: failed

      st = DRIVER_OK ; n_bad = 0_ik
      failed = any(reg%poly(:)%status /= DRIVER_OK)
      !----- The true terminal state, so a restart resumes exactly here. -----------------------------!
      if (reg%cfg%state_write_state .and. .not. failed) call region_write_state(reg)
      if (reg%verbose) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         write(*,'(a)') '  polygon      lat      lon  AGB kgC/m2    LAI  cohorts  area error  energy/water fails'
      end if
      do p = 1_ik, size(reg%poly, kind=ik)
         associate (poly => reg%poly(p))
            a1 = total_area(poly%site)
            if (abs(a1 - 1.0_wp) > 1.0e-5_wp) then
               st = DRIVER_ERR_AREA ; n_bad = n_bad + 1_ik
            end if
            if (reg%verbose)                                                                       &
               write(*,'(i9,2f9.2,f12.4,f7.3,i9,es12.2,2(1x,i0))') poly%id,                         &
                     reg%met_src%cells%lat(poly%cell), reg%met_src%cells%lon(poly%cell),           &
                     total_agb(poly%site),                                                         &
                     total_lai(poly%site), count_cohorts(poly%site), a1 - 1.0_wp,                  &
                     poly%energy_budget%n_fail, poly%water_budget%n_fail
         end associate
      end do
      if (n_bad > 0_ik) write(*,'(a,i0,a)') ' ERROR: ', n_bad, ' polygon(s) did not conserve area'
      if (reg%verbose) write(*,'(a,f0.2,a,i0,a)') ' forcing: ',                                     &
            reg%read_seconds + reg%year_read_seconds, ' s loading in all, ', reg%met_src%n_loads,       &
            ' archive month(s) read from file'

      if (reg%cfg%output%enabled) call output_close(reg%out_files, reg%out_bufs, flush_partial=.not. failed)
      call met_close(reg%met_src)
      if (present(status)) status = st
   end subroutine region_finalize

   subroutine region_free(reg)
      type(meds_region_t), intent(inout) :: reg
      integer(ik) :: p
      if (.not. allocated(reg%poly)) return
      do p = 1_ik, size(reg%poly, kind=ik)
         call site_free(reg%poly(p)%site)
      end do
      deallocate(reg%poly)
      if (allocated(reg%out_bufs)) deallocate(reg%out_bufs)
      if (allocated(reg%cost))     deallocate(reg%cost)
   end subroutine region_free

end module meds_region
