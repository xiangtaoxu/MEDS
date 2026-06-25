!==========================================================================================!
! meds_io_demo -- demographic spin-up that writes the full cohort/patch/site state to a       !
! netCDF file (one snapshot every cfg%io_output_interval_years). Run parameters come from a    !
! meds_config.toml file. Built only when configured with -DMEDS_ENABLE_IO=ON (links netCDF).   !
!                                                                                          !
! Usage:  meds_io_demo [config.toml] [output.nc]                                             !
!==========================================================================================!
program meds_io_demo
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : yr_day
   use meds_config,                 only : meds_config_t, TS_MONTHLY, TS_WEEKLY
   use meds_config_io,              only : load_meds_config
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_free
   use meds_init,                   only : init_bare_ground, init_from_census
   use meds_stepper,                only : advance_one_step
   use meds_demography_diagnostics, only : print_summary
   use meds_io,                     only : meds_io_t, io_create, io_write_snapshot, io_close
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   type(meds_io_t)     :: io
   integer(ik)         :: n_years, nday, step_days, steps_per_year, istep, nsteps
   integer(ik)         :: yday, iyear, month, prev_month
   logical             :: is_new_month, is_new_year, have_cfg, census_ok
   character(len=256)  :: path, ncpath

   path   = 'meds_config.toml'
   ncpath = 'meds_state.nc'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)
   if (command_argument_count() >= 2_ik) call get_command_argument(2, ncpath)

   call load_meds_config(trim(path), cfg, n_years, have_cfg)
   if (have_cfg) write(*,'(2a)') ' config: ', trim(path)

   nday = nint(yr_day, ik)
   select case (cfg%ts_mode)
   case (TS_MONTHLY) ; step_days = 0_ik ; steps_per_year = 12_ik
   case (TS_WEEKLY)  ; step_days = 7_ik ; steps_per_year = nday / step_days
   case default      ; step_days = 1_ik ; steps_per_year = nday
   end select
   nsteps = n_years * steps_per_year

   census_ok = .false.
   if (len_trim(cfg%init_census_file) > 0) then
      call init_from_census(site, cfg, trim(cfg%init_census_file), 298.15_wp, 295.15_wp, census_ok)
      if (census_ok) write(*,'(2a)') ' init  : census ', trim(cfg%init_census_file)
   end if
   if (.not. census_ok) call init_bare_ground(site, cfg, 6_ik, avg_temp = 298.15_wp, min_temp = 295.15_wp)
   call io_create(io, trim(ncpath), cfg)
   write(*,'(a,i0,3a)') ' MEDS IO spin-up: ', n_years, ' yr -> ', trim(ncpath),               &
                        ' (snapshot every output_interval_years)'

   yday = 0_ik ; iyear = 0_ik ; prev_month = 0_ik
   call io_write_snapshot(io, site, 0.0_wp)        ! initial (bare-ground) state

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
         if (mod(iyear, cfg%io_output_interval_years) == 0_ik)                                &
            call io_write_snapshot(io, site, real(iyear, wp))
      end if
   end do

   call print_summary(site, 'final')
   call io_close(io)
   call site_free(site)
   write(*,'(2a)') ' wrote ', trim(ncpath)
end program meds_io_demo
