! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_region -- a region run: every selected forcing cell of a contiguous box as its own        !
! polygon, in one process (MEDS_POLYGON_RUNTIME_PLAN.md §4-6, R2).                               !
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
! [run].n_threads threads, each with a single-threaded patch loop (#183 R3).                       !
!                                                                                          !
! Each polygon sits at its cell's centre, at the cell's orography, in UTC, and starts from bare    !
! ground (OR8). Region files hold the fixed-shape variables; [region].detail_polygons also write   !
! ordinary single-site files, named <prefix>-p<polygon id>, with every variable and tier.          !
!==========================================================================================!
module meds_region
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : day_sec
   use meds_numerics,               only : ascending_order
   use meds_config,                 only : meds_config_t
   use meds_config_io,              only : load_meds_config
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days, time_to_string, &
                                           years_between
   use meds_site_state_types,       only : site_free
   use meds_init,                   only : init_bare_ground
   use meds_forcing_types,          only : met_source_t, met_cells_t
   use meds_forcing_config,         only : MET_PATH_LEN
   use meds_met_driver,             only : met_open, met_close, met_prefetch, MET_OK
   use meds_era5land_reader,        only : era5land_select_box, era5land_path, era5land_default_static, &
                                           ERA_OK, ERA_ERR_NO_CELL
   use meds_diagnostic_reduce,      only : total_area, total_agb, total_lai, count_cohorts
   use meds_polygon,                only : meds_polygon_t, polygon_prepare, polygon_step,        &
                                           DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_AREA, N_PATCH_INIT
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_registry,        only : manager_setup, manager_finalize, manager_alloc_buffers,  &
                                           manager_set_soil_params, manager_restrict_region,      &
                                           activate_site_diag
   use meds_output_manager,         only : output_serialize_pending, output_manager_close,        &
                                           output_serialize_region, output_region_close
   use meds_driver,                 only : apply_io_overrides, ensure_output_dir
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
   ! region_open -- config -> the box's cells -> the shared forcing source -> one bare-ground     !
   ! polygon per cell -> the region files and the detail polygons' own files.                    !
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
      character(len=24)          :: idstr

      ok = .false.
      if (present(verbose)) reg%verbose = verbose

      call load_meds_config(trim(path), reg%cfg)
      if (reg%verbose) write(*,'(2a)') ' config: ', trim(path)
      associate (cfg => reg%cfg)

      !----- The cells: the box's valid cells with enough land, in row-major order (§5). ---------!
      static = cfg%forcing%static_file
      if (len_trim(static) == 0) static = era5land_default_static()
      static = era5land_path(static, cfg%forcing%data_path, '', 0_ik, 0_ik)
      call era5land_select_box(static, cfg%region%box_nwse, cells, st, cfg%region%land_fraction_min)
      if (st == ERA_ERR_NO_CELL) then
         write(*,'(a,4(f0.3,1x),a,f0.2)') ' region: no valid cell in box_nwse = ', cfg%region%box_nwse, &
               'with land fraction >= ', cfg%region%land_fraction_min
         return
      else if (st /= ERA_OK) then
         write(*,'(2a)') ' region: cannot read valid, elevation and land_fraction from ', trim(static)
         return
      end if

      !----- The polygon ids, each cell's row-major index on the global grid. Every detail polygon   !
      !      must be one of them, checked before anything is opened or built. -----------------------!
      n = cells%ncell
      allocate(ids(n))
      ids = cells%row * cells%nlon + cells%col
      do j = 1_ik, cfg%region%n_detail
         if (.not. any(ids == cfg%region%detail_polygons(j))) then
            write(*,'(a,i0,a)') ' region: detail polygon ', cfg%region%detail_polygons(j),            &
                                ' is not a polygon of this region'
            return
         end if
      end do

      call met_open(reg%met_src, cfg%forcing, stat=st, run_start=cfg%start_time,                  &
                    run_end=cfg%end_time, cells=cells)
      if (st /= MET_OK) return
      if (reg%verbose) write(*,'(3a)') ' force : met forcing ON (ED_ERA5land archive ', trim(cfg%forcing%data_path), ')'

      !----- One polygon per cell, each at its cell centre and orography, in UTC, from bare ground. !
      allocate(reg%poly(n), reg%out_bufs(n))
      do p = 1_ik, n
         associate (poly => reg%poly(p))
            poly%cell = p
            poly%id   = ids(p)
            write(poly%label,'(a,i0,a,f0.2,a,f0.2,a)') 'polygon ', poly%id, ' (', cells%lat(p), ', ', &
                                                     cells%lon(p), ')'
            call init_bare_ground(poly%site, cfg, N_PATCH_INIT)
            call polygon_prepare(cfg, reg%met_src, poly, cells%lat(p), cells%lon(p),               &
                                 cells%elevation(p), keep_fast_state=.false.,                     &
                                 keep_soil_carbon=.false., patch_threads=1_ik, verbose=.false.)
         end associate
      end do

      reg%now = cfg%start_time
      reg%step_days = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
      allocate(reg%cost(n)) ; reg%cost = 0.0_wp

      !----- Output. The region files: the configured variables less the ragged and fast ones, with  !
      !      the polygon axis. A detail polygon: its own full single-site file set as well. ----------!
      if (cfg%output%enabled) then
         call ensure_output_dir(trim(cfg%output%dir))
         call manager_setup(reg%out_files, cfg)
         call manager_set_soil_params(reg%out_files, reg%poly(1)%fast_ctx%col_config%soil)
         if (len_trim(cfg%output%io_config) > 0)                                                   &
            call apply_io_overrides(reg%out_files, trim(cfg%output%io_config), reg%verbose)
         call manager_restrict_region(reg%out_files)
         reg%out_files%n_polygon = n
         reg%out_files%polygon_id  = reg%poly(:)%id
         reg%out_files%polygon_row = cells%row ; reg%out_files%polygon_col = cells%col
         reg%out_files%polygon_lat = cells%lat ; reg%out_files%polygon_lon = cells%lon
         call manager_finalize(reg%out_files)
         do p = 1_ik, n
            associate (poly => reg%poly(p))
               call manager_alloc_buffers(reg%out_files, reg%out_bufs(p))
               if (any(cfg%region%detail_polygons(1:cfg%region%n_detail) == poly%id)) then
                  allocate(poly%detail_files, poly%detail_bufs)
                  call manager_setup(poly%detail_files, cfg)
                  write(idstr,'(i0)') poly%id
                  poly%detail_files%prefix = trim(cfg%output%prefix)//'-p'//trim(idstr)
                  call manager_set_soil_params(poly%detail_files, poly%fast_ctx%col_config%soil)
                  if (len_trim(cfg%output%io_config) > 0)                                          &
                     call apply_io_overrides(poly%detail_files, trim(cfg%output%io_config), .false.)
                  call manager_finalize(poly%detail_files)
                  call manager_alloc_buffers(poly%detail_files, poly%detail_bufs)
                  call activate_site_diag(poly%detail_files, poly%site)
               else
                  call activate_site_diag(reg%out_files, poly%site)
               end if
            end associate
         end do
      end if

      if (reg%verbose) then
         write(*,'(a)') '==================== MEDS region run ===================='
         write(*,'(a,4(1x,f0.3),a,i0,a,f0.2,a)') ' region: box_nwse =', cfg%region%box_nwse, ', ',    &
               n, ' polygons (land fraction >= ', cfg%region%land_fraction_min, ')'
         if (cfg%region%n_detail > 0) write(*,'(a,i0,a)') ' region: ', cfg%region%n_detail,        &
               ' detail polygon(s) with single-site files'
         write(*,'(5a)') ' run   : ', time_to_string(cfg%start_time), ' -> ', time_to_string(cfg%end_time)
         write(*,'(a,f0.2,a)') ' span  : ', years_between(reg%now, cfg%end_time), ' yr'
      end if
      end associate
      ok = .true.
   end subroutine region_open

   !---------------------------------------------------------------------------------------!
   ! region_step_month -- one month of the whole region: load the month's forcing for every cell, !
   ! advance each polygon through it, then write the queued output (§4). The month runs from the  !
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
      call met_prefetch(reg%met_src, reg%now)

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
            call polygon_step(cfg, reg%met_src, reg%out_files, reg%out_bufs(p), reg%poly(p),         &
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

      if (new_year(nstep) .and. reg%verbose) call region_summary(reg, time_to_string(reg%now))

      !----- The I/O phase. -----------------------------------------------------------------------!
      call io_phase(reg)
      end associate
   end subroutine region_step_month

   !----- Write every polygon's queued records: the region files, then each detail polygon's own. --!
   subroutine io_phase(reg)
      type(meds_region_t), intent(inout) :: reg
      integer(ik) :: p
      if (.not. reg%cfg%output%enabled) return
      call output_serialize_region(reg%out_files, reg%out_bufs)
      do p = 1_ik, size(reg%poly, kind=ik)
         if (allocated(reg%poly(p)%detail_bufs))                                                  &
            call output_serialize_pending(reg%poly(p)%detail_files, reg%poly(p)%detail_bufs)
      end do
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
   end subroutine region_summary

   !---------------------------------------------------------------------------------------!
   ! region_finalize -- the per-polygon status table (state, conservation), close every file and   !
   ! the forcing source. DRIVER_ERR_AREA if any polygon's patch areas no longer sum to 1.          !
   !---------------------------------------------------------------------------------------!
   subroutine region_finalize(reg, status)
      type(meds_region_t),   intent(inout) :: reg
      integer(ik), optional, intent(out)   :: status
      integer(ik) :: p, st, n_bad
      real(wp)    :: a1

      st = DRIVER_OK ; n_bad = 0_ik
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

      if (reg%cfg%output%enabled) then
         call output_region_close(reg%out_files, reg%out_bufs)
         do p = 1_ik, size(reg%poly, kind=ik)
            if (allocated(reg%poly(p)%detail_bufs))                                               &
               call output_manager_close(reg%poly(p)%detail_files, reg%poly(p)%detail_bufs, .true.)
         end do
      end if
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
