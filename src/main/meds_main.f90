!==========================================================================================!
! meds_main -- THE executable entry point of MEDS: parse the command line, drive the model,   !
! translate a driver status into an exit code.                                                 !
!                                                                                          !
! The model itself lives in `meds_driver` (open / step / finalize). This program used to BE     !
! that logic -- 400 lines of it, with no seam -- which meant the only way to run MEDS was to    !
! exec this binary. It is now one of two callers; the other is the C-API shim `meds_capi_run`,  !
! which is what lets examples/example_biophysics drive the full coupled model from Python.      !
! Keeping the program this thin is the point: anything added here is, by construction,          !
! unreachable from Python.                                                                      !
!                                                                                          !
! Usage:  meds_main [main.toml]           (defaults to ./meds_config_main.toml; hard error if a !
!                                          config file or any required parameter is missing --  !
!                                          there are NO built-in defaults)                       !
!         meds_main --dump-io-config [main.toml]                                                 !
!                                          write a fully-commented meds_io_config.toml listing   !
!                                          every registered output variable, then exit.          !
!==========================================================================================!
program meds_main
   use meds_kinds,            only : ik
   use meds_config,           only : meds_config_t
   use meds_config_io,        only : load_meds_config
   use meds_output_types,     only : output_registry_t
   use meds_output_registry,  only : build_output_registry, dump_io_config
   use meds_driver,           only : meds_run_t, driver_open, driver_step, driver_finalize,    &
                                     driver_free, driver_done, DRIVER_ERR_NAN, DRIVER_ERR_AREA,   &
                                     DRIVER_ERR_SOILC
   implicit none

   type(meds_run_t)   :: run
   integer(ik)        :: status
   logical            :: ok
   character(len=256) :: path

   path = 'meds_config_main.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)

   !----- `--dump-io-config [main.toml]` is the discoverability half of per-variable output      !
   !      control: the override mechanism already worked, but nothing told a user which names    !
   !      exist. It builds the registry and exits without running anything, so it stays here in  !
   !      the CLI rather than in the driver.  ----------------------------------------------------!
   if (trim(path) == '--dump-io-config') then
      block
         type(meds_config_t)     :: cfg
         type(output_registry_t) :: reg_dump
         path = 'meds_config_main.toml'
         if (command_argument_count() >= 2_ik) call get_command_argument(2, path)
         call load_meds_config(trim(path), cfg)
         write(*,'(2a)') ' config: ', trim(path)
         call build_output_registry(reg_dump, cfg)
         call dump_io_config(reg_dump, 'meds_io_config.toml')
      end block
      stop
   end if

   call driver_open(trim(path), run, ok)
   if (.not. ok) error stop 'meds_main: could not open the run'

   do while (.not. driver_done(run))
      call driver_step(run, status)
      !----- The driver RETURNS what this program used to `error stop` on, so that a library      !
      !      caller (Python) survives it. In the executable it is still fatal.  -------------------!
      if (status == DRIVER_ERR_NAN)   error stop 'meds_main: NaN detected in state'
      if (status == DRIVER_ERR_SOILC) error stop 'meds_main: impossible soil-carbon pool'
   end do

   call driver_finalize(run, status)
   if (status == DRIVER_ERR_AREA) error stop 'meds_main: site area not conserved'
   write(*,'(a)') ' OK: simulation completed, area conserved, no NaNs.'

   call driver_free(run)

end program meds_main
