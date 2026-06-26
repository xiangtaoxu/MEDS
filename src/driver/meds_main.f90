!==========================================================================================!
! meds_main -- THE entry point of MEDS. It reads the run configuration, builds the initial      !
! community, runs the demographic spin-up, saves the netCDF state output, and exits. (Merged    !
! from the former app/meds_demo + app/meds_io_demo into one driver.)                            !
!                                                                                          !
! Per step the engine does recruitment -> growth (light competition via overtopping LAI) ->     !
! shade-driven mortality -> cohort & patch fusion/fission -> treefall disturbance, all through   !
! meds_stepper and the meds_demography_interface seam; an explicit calendar lets daily/weekly/   !
! monthly modes share the one engine call. Run, initial-condition and output parameters all     !
! come from a meds_config.toml file ([run], [init], [demography], ..., [io]).                    !
!                                                                                          !
! Output goes to <[io].output_dir>/<[io].output_prefix>.nc, every [io].output_interval_years.    !
! The netCDF writer is active only when MEDS is built with -DMEDS_ENABLE_IO=ON; otherwise the    !
! I/O layer is a no-op stub and the run still proceeds (printing its summary, writing nothing).  !
!                                                                                          !
! Usage:  meds_main [config.toml]   (defaults to ./meds_config.toml; built-in defaults if absent)!
!==========================================================================================!
program meds_main
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : yr_day
   use meds_config,                 only : meds_config_t, TS_MONTHLY, TS_WEEKLY
   use meds_config_io,              only : load_meds_config
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_free
   use meds_init,                   only : init_bare_ground, init_from_census
   use meds_stepper,                only : advance_one_step
   use meds_demography_diagnostics, only : print_summary, total_area, has_nan
   use meds_io,                     only : meds_io_t, io_create, io_write_snapshot, io_close
   implicit none

   integer(ik), parameter :: N_PATCH_INIT = 6_ik    ! bare-ground patches when no census is given

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   type(meds_io_t)     :: io
   integer(ik)         :: n_years, steps_per_year, nday, step_days
   integer(ik)         :: istep, nsteps, iyear, yday, month, prev_month
   logical             :: is_new_month, is_new_year, have_cfg, census_ok
   real(wp)            :: a0, a1
   character(len=256)  :: path
   character(len=512)  :: ncfile

   !----- 1. Read the run configuration. -----------------------------------------------!
   path = 'meds_config.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)
   call load_meds_config(trim(path), cfg, n_years, have_cfg)
   if (have_cfg) then
      write(*,'(2a)') ' config: ', trim(path)
   else
      write(*,'(3a)') ' config: ', trim(path), ' not found -- using built-in defaults'
   end if

   !----- Calendar: number of steps per year for the chosen time step. -----------------!
   nday = nint(yr_day, ik)
   select case (cfg%ts_mode)
   case (TS_MONTHLY) ; step_days = 0_ik ; steps_per_year = 12_ik
   case (TS_WEEKLY)  ; step_days = 7_ik ; steps_per_year = nday / step_days
   case default      ; step_days = 1_ik ; steps_per_year = nday
   end select
   nsteps = n_years * steps_per_year

   !----- 2. Build the initial community: a cohort census if configured, else bare ground.!
   census_ok = .false.
   if (len_trim(cfg%init_census_file) > 0) then
      call init_from_census(site, cfg, trim(cfg%init_census_file), 298.15_wp, 295.15_wp, census_ok)
      if (census_ok) then
         write(*,'(2a)') ' init  : census ', trim(cfg%init_census_file)
      else
         write(*,'(3a)') ' init  : census ', trim(cfg%init_census_file),                      &
                         ' not usable -- falling back to bare ground'
      end if
   end if
   if (.not. census_ok) call init_bare_ground(site, cfg, N_PATCH_INIT,                        &
                                              avg_temp = 298.15_wp, min_temp = 295.15_wp)
   a0 = total_area(site)

   write(*,'(a)') '==================== MEDS demographic spin-up ===================='
   write(*,'(a,i0,a,i0)') ' years=', n_years, '  steps/yr=', steps_per_year

   !----- 3. Open the netCDF output (path from [io]; a no-op without -DMEDS_ENABLE_IO). -!
   if (cfg%io_write_output) then
      ncfile = trim(cfg%io_output_dir)//'/'//trim(cfg%io_output_prefix)//'.nc'
      call ensure_output_dir(trim(cfg%io_output_dir))
      call io_create(io, trim(ncfile), cfg)
      call io_write_snapshot(io, site, 0.0_wp)         ! initial state
   end if

   write(*,'(a)') '-----------------------------------------------------------------------------'
   call print_summary(site, 'year 0')

   !----- 4. Run the simulation. -------------------------------------------------------!
   prev_month = 0_ik ; yday = 0_ik ; iyear = 0_ik
   do istep = 1_ik, nsteps
      if (cfg%ts_mode == TS_MONTHLY) then
         month        = mod(istep - 1_ik, 12_ik) + 1_ik
         is_new_month = .true.
         is_new_year  = (month == 1_ik)
      else
         yday        = yday + step_days
         is_new_year = .false.
         if (yday > nday) then
            yday        = yday - nday
            is_new_year = .true.
         end if
         month        = min(12_ik, (yday - 1_ik) * 12_ik / nday + 1_ik)
         is_new_month = (month /= prev_month) .or. is_new_year .or. (istep == 1_ik)
         if (istep == 1_ik) is_new_year = .true.
         prev_month   = month
      end if

      call advance_one_step(site, cfg, is_new_month, is_new_year)

      if (is_new_year) then
         iyear = iyear + 1_ik
         if (mod(iyear, 5_ik) == 0_ik .or. iyear == 1_ik) call print_summary(site, 'year '//itoa(iyear))
         if (has_nan(site)) error stop 'meds_main: NaN detected in state'
         if (cfg%io_write_output .and. mod(iyear, cfg%io_output_interval_years) == 0_ik)       &
            call io_write_snapshot(io, site, real(iyear, wp))
      end if
   end do

   !----- 5. Finalize: summary, conservation check, close output, free state. -----------!
   call print_summary(site, 'final')
   a1 = total_area(site)
   write(*,'(a)') '-----------------------------------------------------------------------------'
   write(*,'(a,f12.9,a,f12.9)') ' site area start=', a0, '  end=', a1
   if (abs(a1 - 1.0_wp) > 1.0e-5_wp) error stop 'meds_main: site area not conserved'
   if (cfg%io_write_output) call io_close(io)
   write(*,'(a)') ' OK: simulation completed, area conserved, no NaNs.'

   call site_free(site)

contains

   function itoa(i) result(s)
      integer(ik), intent(in)       :: i
      character(len=:), allocatable :: s
      character(len=12)             :: buf
      write(buf,'(i0)') i
      s = trim(buf)
   end function itoa

   !----- Create the output directory if it does not exist (driver-level convenience).    !
   !       Filesystem access stays in the driver; the engine/library never touches it.    !
   subroutine ensure_output_dir(dir)
      character(len=*), intent(in) :: dir
      integer :: stat
      if (len_trim(dir) == 0 .or. trim(dir) == '.') return
      call execute_command_line('mkdir -p "'//trim(dir)//'"', wait=.true., exitstat=stat)
      if (stat /= 0) write(*,'(3a)') ' warning: could not create output dir "', trim(dir), '"'
   end subroutine ensure_output_dir

end program meds_main
