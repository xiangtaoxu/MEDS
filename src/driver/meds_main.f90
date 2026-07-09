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
! Two output streams (both under [io]; netCDF is always compiled):                               !
!   * DIAGNOSTIC timeseries -> <dir>/<prefix>-D-output.nc, every [io].output_interval_years.      !
!   * STATE checkpoints (restart) -> <dir>/<prefix>-S-<YYYYMMDDHHMMSS>.nc, every                  !
!     [io].state_interval_years (and the final year) when [io].write_state is true.               !
! Initialization is selected by [init].init_mode: 0 near-bare ground, 1 cohort census             !
! ([init].census_file), 2 restart from a state file ([init].restart_file). The file for the        !
! non-selected mode is read into the config but ignored; unusable input falls back to bare ground. !
!                                                                                          !
! Usage:  meds_main [main.toml]   (defaults to ./meds_config_main.toml; hard error if a config !
!         file or any required parameter is missing -- there are NO built-in defaults).         !
!==========================================================================================!
program meds_main
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : day_sec, yr_day
   use meds_config,                 only : meds_config_t, INIT_CENSUS, INIT_RESTART
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days,            &
                                           time_to_string, years_between
   use meds_config_io,              only : load_meds_config, write_pft_params_csv
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_free
   use meds_init,                   only : init_bare_ground, init_from_census
   use meds_stepper,                only : advance_one_step
   use meds_fast_loop,              only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_forcing_types,          only : met_driver_t
   use meds_met_driver,             only : met_open, met_close
   use meds_demography_diagnostics, only : print_summary, total_area, has_nan
   use meds_io,                     only : meds_io_t, io_create, io_write_snapshot, io_close,   &
                                           io_write_state, io_read_state
   implicit none

   integer(ik), parameter :: N_PATCH_INIT = 6_ik    ! bare-ground patches when no census/restart

   type(meds_config_t)  :: cfg
   type(site_t)         :: site
   type(meds_io_t)      :: io
   type(fast_context_t) :: fast_ctx        ! sub-daily biophysics context (built only if fast_biophysics_on)
   type(met_driver_t)   :: met_drv         ! live met-forcing reader (opened only if forcing_on)
   type(meds_time_t)   :: now, prev, restart_time
   integer(ik)         :: steps_per_year, istep, iyear, step_days
   logical             :: is_new_month, is_new_year, init_ok
   real(wp)            :: a0, a1
   character(len=256)  :: path
   character(len=512)  :: ncfile
   character(len=19)   :: datestr

   !----- 1. Read the run configuration. -----------------------------------------------!
   path = 'meds_config_main.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)
   call load_meds_config(trim(path), cfg)         ! hard error if a file or required key is missing
   write(*,'(2a)') ' config: ', trim(path)

   !----- 2. Build the initial community per [init].init_mode (0 bare | 1 census | 2 restart);!
   !         the file for the non-selected mode is ignored. Unusable input -> bare ground. A   !
   !         successful restart also recovers the calendar date to continue from.  ------------!
   init_ok = .false. ; now = cfg%start_time
   select case (cfg%init_mode)
   case (INIT_RESTART)
      call io_read_state(site, cfg, trim(cfg%init_restart_file), restart_time, init_ok)  ! stub -> .false.
      if (init_ok) then
         now = restart_time
         write(*,'(4a)') ' init  : restart (mode 2) ', trim(cfg%init_restart_file),           &
                         ' @ ', time_to_string(now)
      else
         write(*,'(3a)') ' init  : restart (mode 2) ', trim(cfg%init_restart_file),           &
                         ' not usable -- falling back to bare ground'
      end if
   case (INIT_CENSUS)
      call init_from_census(site, cfg, trim(cfg%init_census_file), init_ok)
      if (init_ok) then
         write(*,'(2a)') ' init  : census (mode 1) ', trim(cfg%init_census_file)
      else
         write(*,'(3a)') ' init  : census (mode 1) ', trim(cfg%init_census_file),             &
                         ' not usable -- falling back to bare ground'
      end if
   end select
   if (.not. init_ok) then
      call init_bare_ground(site, cfg, N_PATCH_INIT)
      write(*,'(a)') ' init  : bare ground (mode 0)'
   end if
   a0 = total_area(site)

   !----- 2b. Fast biophysics context (opt-in): build the static column config + seed the per-  !
   !          patch CAS/soil reservoirs ONCE. Skipped entirely when fast_biophysics_on=.false.   !
   !          (the default), so the demographic-only run is untouched.  ------------------------!
   if (cfg%fast_biophysics_on) then
      call build_fast_context(cfg, fast_ctx)
      !----- Live met forcing (opt-in): open the reader ONCE and thread the reference height into  !
      !      the aerodynamic reference. Off -> the fast loop runs the constant-forcing MVP.  -------!
      if (cfg%forcing%forcing_on) then
         call met_open(met_drv, cfg%forcing)
         fast_ctx%zref = cfg%forcing%reference_height
         write(*,'(3a)') ' force : met forcing ON (', trim(cfg%forcing%path), ')'
      end if
      call init_fast_reservoirs(site, fast_ctx)
      write(*,'(a)') ' fast  : sub-daily biophysics ON'
   end if

   step_days      = max(1_ik, nint(cfg%dt_slow / day_sec, ik))   ! calendar advance per slow step
   steps_per_year = max(1_ik, nint(yr_day / real(step_days, wp), ik))   ! header line only

   write(*,'(a)') '==================== MEDS demographic spin-up ===================='
   write(*,'(5a)')        ' run   : ', time_to_string(cfg%start_time), ' -> ',                &
                          time_to_string(cfg%end_time)
   write(*,'(a,f0.2,a,i0)') ' span  : ', years_between(now, cfg%end_time),                    &
                          ' yr  steps/yr~', steps_per_year

   !----- 3. Open the DIAGNOSTIC output (netCDF). ---------------------------------------!
   if (cfg%io_write_output .or. cfg%io_write_state) then
      call ensure_output_dir(trim(cfg%io_output_dir))
      !----- Dump the PFT parameter table next to the output (provenance; netCDF-free). ------!
      call write_pft_params_csv(cfg, trim(cfg%io_output_dir)//'/'//trim(cfg%io_output_prefix)// &
                                '_pft_parameters.csv')
   end if
   if (cfg%io_write_output) then
      ncfile = trim(cfg%io_output_dir)//'/'//trim(cfg%io_output_prefix)//'-D-output.nc'
      call io_create(io, trim(ncfile), cfg)
      call io_write_snapshot(io, site, now)            ! initial state at the start date
   end if

   write(*,'(a)') '-----------------------------------------------------------------------------'
   call print_summary(site, 'start')

   !----- 4. Run the simulation on the real calendar until the end date. Each step advances  !
   !         the date (leap years exact); a year/month roll-over sets the structural cadence. !
   istep = 0_ik ; iyear = 0_ik
   do while (time_lt(now, cfg%end_time))
      prev = now
      now  = time_advance_days(prev, step_days)      ! advance the calendar by the slow step (dt_slow)
      istep = istep + 1_ik

      is_new_year  = now%year  /= prev%year
      is_new_month = is_new_year .or. (now%month /= prev%month)

      if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) then
         call advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx,                   &
                               met_drv=met_drv, step_start=prev)   ! forcing spans [prev, now]
      else
         call advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx)
      end if

      if (is_new_year) then
         iyear = iyear + 1_ik
         if (mod(iyear, 5_ik) == 0_ik .or. iyear == 1_ik) then
            datestr = time_to_string(now)
            call print_summary(site, datestr(1:10))               ! date only (fits the label)
         end if
         if (has_nan(site)) error stop 'meds_main: NaN detected in state'
         if (cfg%io_write_output .and. mod(iyear, cfg%io_output_interval_years) == 0_ik)       &
            call io_write_snapshot(io, site, now)
         if (cfg%io_write_state .and. mod(iyear, cfg%io_state_interval_years) == 0_ik)         &
            call io_write_state(site, cfg, trim(cfg%io_output_dir), trim(cfg%io_output_prefix), now)
      end if
   end do

   !----- Always checkpoint the true terminal state so a restart resumes exactly here. ------!
   if (cfg%io_write_state)                                                                    &
      call io_write_state(site, cfg, trim(cfg%io_output_dir), trim(cfg%io_output_prefix), now)

   !----- 5. Finalize: summary, conservation check, close output, free state. -----------!
   call print_summary(site, 'final')
   a1 = total_area(site)
   write(*,'(a)') '-----------------------------------------------------------------------------'
   write(*,'(a,f12.9,a,f12.9)') ' site area start=', a0, '  end=', a1
   if (abs(a1 - 1.0_wp) > 1.0e-5_wp) error stop 'meds_main: site area not conserved'
   if (cfg%io_write_output) call io_close(io)
   if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) call met_close(met_drv)
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
