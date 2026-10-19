! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_main -- THE executable entry point of MEDS: parse the command line, drive the run,       !
! translate a run status into an exit code.                                                    !
!                                                                                          !
! The model itself lives in `meds_region` (open / step / finalize), and this program is one of !
! its two callers; the other is the C-API shim `meds_c_api_run`, which lets Python drive the   !
! full coupled model (`meds.model`). Keeping the program this thin is the point: anything      !
! added here is, by construction, unreachable from Python.                                     !
!                                                                                          !
! Usage:  meds_main [main.toml]           (defaults to ./meds_config_main.toml; hard error if a !
!                                          config file or any required parameter is missing --  !
!                                          see docs/configuration.md)                            !
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
   use meds_region,           only : meds_region_t, region_open, region_step_month, region_done,  &
                                     region_finalize, region_free, DRIVER_ERR_NAN, DRIVER_ERR_AREA, &
                                     DRIVER_ERR_BOUNDS
   implicit none

   type(meds_region_t) :: reg
   integer(ik)         :: status
   logical             :: ok
   character(len=256)  :: path

   path = 'meds_config_main.toml'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)

   !----- `--dump-io-config [main.toml]` is the discoverability half of per-variable output      !
   !      control: it lists every name an override can switch. It builds the registry and exits  !
   !      without running anything, so it stays here in the CLI rather than in the run. ---------!
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

   call region_open(trim(path), reg, ok)
   if (.not. ok) error stop 'meds_main: could not open the run'
   do while (.not. region_done(reg))
      call region_step_month(reg, status)
      !----- The run RETURNS a failed step as a status, so that a library caller (Python) survives  !
      !      it. In the executable it is fatal.  ---------------------------------------------------!
      if (status == DRIVER_ERR_NAN)    error stop 'meds_main: NaN in the state or the energy/water ledger'
      if (status == DRIVER_ERR_BOUNDS) error stop 'meds_main: a cohort or soil-carbon state out of bounds'
   end do
   call region_finalize(reg, status)
   if (status == DRIVER_ERR_AREA) error stop 'meds_main: a polygon''s area was not conserved'
   write(*,'(a)') ' OK: simulation completed, area conserved, no NaNs.'
   call region_free(reg)

end program meds_main
