!==========================================================================================!
! meds_io_demo -- demographic spin-up that writes the full cohort/patch/site state to a       !
! netCDF file (one snapshot every cfg%io_output_interval_years). Built only when the project  !
! is configured with -DMEDS_ENABLE_IO=ON (links the netCDF C library).                        !
!                                                                                          !
! Usage:  meds_io_demo [n_years] [output.nc]                                                 !
!==========================================================================================!
program meds_io_demo
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : yr_day
   use meds_config,                 only : meds_config_t, build_config, validate_config, TS_DAILY
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_free
   use meds_setup,                  only : init_bare_ground
   use meds_stepper,                only : advance_one_step
   use meds_demography_diagnostics, only : print_summary
   use meds_io,                     only : meds_io_t, io_create, io_write_snapshot, io_close
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   type(meds_io_t)     :: io
   integer(ik)         :: n_years, nday, istep, nsteps, yday, iyear, month, prev_month
   logical             :: is_new_month, is_new_year
   character(len=256)  :: arg, path

   n_years = 60_ik
   path    = 'meds_state.nc'
   if (command_argument_count() >= 1_ik) then
      call get_command_argument(1, arg) ; read(arg, *) n_years
   end if
   if (command_argument_count() >= 2_ik) call get_command_argument(2, path)

   cfg = build_config(ts_mode = TS_DAILY)
   call validate_config(cfg)
   call init_bare_ground(site, cfg, 6_ik, avg_temp = 298.15_wp, min_temp = 295.15_wp)
   call io_create(io, trim(path), cfg)

   write(*,'(a,i0,3a)') ' MEDS IO spin-up: ', n_years, ' yr -> ', trim(path), ' (annual snapshots)'

   nday = nint(yr_day, ik) ; yday = 0_ik ; iyear = 0_ik ; prev_month = 0_ik
   nsteps = n_years * nday
   call io_write_snapshot(io, site, 0.0_wp)        ! initial (bare-ground) state

   do istep = 1_ik, nsteps
      yday        = yday + 1_ik
      is_new_year = .false.
      if (yday > nday) then
         yday        = yday - nday
         is_new_year = .true.
      end if
      month        = min(12_ik, (yday - 1_ik) * 12_ik / nday + 1_ik)
      is_new_month = (month /= prev_month) .or. is_new_year .or. (istep == 1_ik)
      if (istep == 1_ik) is_new_year = .true.
      prev_month   = month

      call advance_one_step(site, cfg, is_new_month, is_new_year)

      if (is_new_year) then
         iyear = iyear + 1_ik
         if (mod(iyear, cfg%io_output_interval_years) == 0_ik) then
            call io_write_snapshot(io, site, real(iyear, wp))
         end if
      end if
   end do

   call print_summary(site, 'final')
   call io_close(io)
   call site_free(site)
   write(*,'(2a)') ' wrote ', trim(path)
end program meds_io_demo
