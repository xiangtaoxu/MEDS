! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_driver -- the coupled model as an OPEN / STEP / FINALIZE object, so something other      !
! than the `meds_main` program can drive it.                                                    !
!                                                                                          !
! This is `meds_main`'s body, lifted verbatim. That program was 400 lines of driver logic with   !
! no seam in it: configuration, initial community, fast context, met reader, output manager,     !
! the calendar loop and the closing conservation reports were all statements in one PROGRAM, so  !
! the ONLY way to run MEDS was to exec the binary. `meds_main` is now a thin shell over this      !
! module and the C-API shim `meds_c_api_run` is a second caller -- which is what lets             !
! `examples/example_biophysics` drive the full coupled model from Python.                          !
!                                                                                          !
!   type(meds_run_t) :: run                                                                      !
!   call driver_open('meds_config_main.toml', run, ok)                                            !
!   do while (.not. driver_done(run)) ; call driver_step(run, status) ; end do                     !
!   call driver_finalize(run) ; call driver_free(run)                                              !
!                                                                                          !
! ONE deliberate behaviour change: the run loop's `error stop` on a NaN state is now a STATUS      !
! RETURN. In a program those are the same thing, but this module is compiled into a shared         !
! library that a Python interpreter dlopens, and `error stop` there takes the interpreter down     !
! with it -- no traceback, no chance to inspect the state that went bad. `meds_main` re-raises      !
! it as the same `error stop`, so the executable's behaviour is unchanged.                          !
!==========================================================================================!
module meds_driver
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : day_sec, yr_day
   use meds_config,                 only : meds_config_t, INIT_CENSUS, INIT_RESTART
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days,            &
                                           time_to_string, years_between
   use meds_config_io,              only : load_meds_config, write_pft_params_csv
   use meds_site_state_types,       only : site_t, site_free
   use meds_demography_update,      only : update_overtopping_lai
   use meds_init,                   only : init_bare_ground, init_from_census
   use meds_vegetation_dynamics,    only : advance_plant_traits
   use meds_forcing_types,          only : met_source_t
   use meds_met_driver,             only : met_open, met_close, met_prefetch
   use meds_forcing_config,         only : MET_BACKEND_ERA5LAND
   use meds_diagnostic_reduce,      only : print_summary, total_area
   use meds_polygon,                only : meds_polygon_t, polygon_prepare, polygon_step,        &
                                           polygon_report, DRIVER_OK, DRIVER_FINISHED,           &
                                           DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_SOILC,      &
                                           N_PATCH_INIT
   use meds_io,                     only : state_write_state, io_read_state
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_registry,        only : manager_setup, manager_finalize, manager_alloc_buffers, &
                                           manager_set_soil_params, activate_site_diag,         &
                                           apply_variable_override, parse_stream_mask,          &
                                           build_freq_index, OVR_TRUE, OVR_FALSE, OVR_MASK
   use meds_output_manager,         only : output_serialize_pending, output_manager_close
   use meds_toml,                   only : toml_table_t, toml_parse_file
   implicit none
   private

   public :: meds_run_t, driver_open, driver_step, driver_finalize, driver_free, driver_done
   public :: apply_io_overrides, ensure_output_dir
   !----- The status codes live with the step (meds_polygon); re-exported for the driver's callers. !
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_SOILC


   !----- Everything the calendar loop needs between steps. These were meds_main's locals; making  !
   !      them components is the whole extraction -- no state hides in module scope, so two runs    !
   !      can be open at once (which is exactly what the C-API's handle registry does). A site run  !
   !      is one polygon (meds_polygon): the site, its fast context, forcing cursor, output buffers  !
   !      and ledgers live there; the run holds what a region would share.                          !
   type :: meds_run_t
      type(meds_config_t)    :: cfg
      type(meds_polygon_t)   :: poly            !< the site
      type(output_files_t)   :: out_files          !< the run's output files (built only if output.enabled)
      type(met_source_t)     :: met_src         !< opened only if forcing_on
      type(meds_time_t)      :: now, prev
      integer(ik)            :: istep = 0_ik, iyear = 0_ik
      integer(ik)            :: step_days = 1_ik, steps_per_year = 365_ik
      logical                :: is_open = .false.
      logical                :: verbose = .true.  !< the progress lines meds_main prints
   end type meds_run_t

contains

   !---------------------------------------------------------------------------------------!
   ! driver_done -- has the calendar reached end_time? The loop predicate, exposed so a caller   !
   ! (Python included) can write its own `while` rather than being handed a callback.            !
   !---------------------------------------------------------------------------------------!
   pure logical function driver_done(run)
      type(meds_run_t), intent(in) :: run
      driver_done = .not. time_lt(run%now, run%cfg%end_time)
   end function driver_done

   !---------------------------------------------------------------------------------------!
   ! driver_open -- config -> initial community -> fast context -> met reader -> soil-carbon      !
   ! spin-up -> output manager. Everything meds_main did before its `do while`.                   !
   !---------------------------------------------------------------------------------------!
   subroutine driver_open(path, run, ok, verbose)
      character(len=*),  intent(in)    :: path
      type(meds_run_t),  intent(inout) :: run
      logical,           intent(out)   :: ok
      logical, optional, intent(in)    :: verbose
      type(meds_time_t) :: restart_time
      logical           :: init_ok, fast_state_found

      ok = .false.
      if (present(verbose)) run%verbose = verbose

      !----- RESET the run's counters before anything writes to one (polygon_prepare resets the    !
      !      site's ledgers and budgets).                                                           !
      !                                                                                          !
      !      `meds_run_t` has default initialisers, which is enough for a fresh variable -- and   !
      !      NOT enough here. The C-API keeps its runs in a module-`save` registry and hands out   !
      !      freed slots again, so a second run in the same process opens on top of the first      !
      !      one's counters, ledgers and budgets. examples/example_biophysics does exactly that:    !
      !      run_example.py opens the spin-up, closes it, then opens the July stage in the SAME     !
      !      process. Caught when both stages reported the identical worst seam gap to four         !
      !      significant figures -- stage 2 was still carrying stage 1's maximum, and would have    !
      !      carried its step counters and both conservation ledgers too.  ------------------------!
      run%istep = 0_ik
      run%iyear = 0_ik

      !----- 1. Read the run configuration. --------------------------------------------------!
      call load_meds_config(trim(path), run%cfg)   ! hard error if a file or required key is missing
      if (run%verbose) write(*,'(2a)') ' config: ', trim(path)

      !----- 2. Initial community per [init].init_mode (0 bare | 1 census | 2 restart); the file  !
      !         for the non-selected mode is ignored. Unusable input -> bare ground. A successful  !
      !         restart also recovers the calendar date to continue from.  ------------------------!
      init_ok = .false. ; fast_state_found = .false. ; run%now = run%cfg%start_time
      select case (run%cfg%init_mode)
      case (INIT_RESTART)
         call io_read_state(run%poly%site, run%cfg, trim(run%cfg%init_restart_file), restart_time,   &
                            init_ok, fast_found=fast_state_found,                                &
                            restructure_pending=run%poly%restructure_pending,                    &
                            restructure_new_year=run%poly%restructure_new_year)
         if (init_ok) then
            run%now = restart_time
            if (run%verbose) then
               write(*,'(4a)') ' init  : restart (mode 2) ', trim(run%cfg%init_restart_file),   &
                               ' @ ', time_to_string(run%now)
               if (fast_state_found) then
                  write(*,'(a)') '         (CAS/soil/snow fast reservoirs restored -- exact restart)'
               else
                  write(*,'(a)') '         (state file predates fast-reservoir persistence -- re-seeding)'
               end if
               if (run%poly%restructure_pending)                                                  &
                  write(*,'(a)') '         (the boundary''s restructuring is pending: it runs before the first step)'
            end if
         else if (run%verbose) then
            write(*,'(3a)') ' init  : restart (mode 2) ', trim(run%cfg%init_restart_file),      &
                            ' not usable -- falling back to bare ground'
         end if
      case (INIT_CENSUS)
         call init_from_census(run%poly%site, run%cfg, trim(run%cfg%init_census_file), init_ok)
         if (run%verbose) then
            if (init_ok) then
               write(*,'(2a)') ' init  : census (mode 1) ', trim(run%cfg%init_census_file)
            else
               write(*,'(3a)') ' init  : census (mode 1) ', trim(run%cfg%init_census_file),     &
                               ' not usable -- falling back to bare ground'
            end if
         end if
      end select
      if (.not. init_ok) then
         call init_bare_ground(run%poly%site, run%cfg, N_PATCH_INIT)
         if (run%verbose) write(*,'(a)') ' init  : bare ground (mode 0)'
      end if

      !----- 2a. Census restart with plasticity ON: census cohorts sit in an established stand but !
      !          carry NO trait history, so acclimate their leaf traits to the current light        !
      !          environment INSTANTANEOUSLY (after the competition sweep). Bare ground legitimately !
      !          starts at top-of-canopy; a state restart already read the plastic traits from file. !
      if (run%cfg%trait_plasticity_on .and. run%cfg%init_mode == INIT_CENSUS .and. init_ok) then
         call update_overtopping_lai(run%poly%site)
         call advance_plant_traits(run%poly%site, run%cfg, run%cfg%dt_years, instantaneous=.true.)
      end if

      !----- 2b. The forcing source (opt-in), then the site as a polygon: fast context and          !
      !          reservoirs, forcing cursor at the [site] location, snow, soil-carbon spin-up.  ----!
      if (run%cfg%fast_biophysics_on .and. run%cfg%forcing%forcing_on) then
         call met_open(run%met_src, run%cfg%forcing, run_start=run%cfg%start_time,                &
                       run_end=run%cfg%end_time)
         if (run%verbose) then
            if (run%cfg%forcing%backend == MET_BACKEND_ERA5LAND) then
               write(*,'(3a)') ' force : met forcing ON (ED_ERA5land archive ', trim(run%cfg%forcing%data_path), ')'
            else
               write(*,'(3a)') ' force : met forcing ON (', trim(run%cfg%forcing%path), ')'
            end if
         end if
      end if
      call polygon_prepare(run%cfg, run%met_src, run%poly, run%cfg%forcing%latitude_deg,          &
                           run%cfg%forcing%longitude_deg, run%cfg%forcing%utc_offset_h,           &
                           run%cfg%forcing%elevation_m,                                           &
                           keep_fast_state=run%cfg%init_mode == INIT_RESTART .and. init_ok .and.  &
                                           fast_state_found,                                     &
                           keep_soil_carbon=run%cfg%init_mode == INIT_RESTART .and. init_ok,      &
                           verbose=run%verbose)

      run%step_days      = max(1_ik, nint(run%cfg%dt_slow / day_sec, ik))
      run%steps_per_year = max(1_ik, nint(yr_day / real(run%step_days, wp), ik))

      if (run%verbose) then
         write(*,'(a)') '==================== MEDS demographic spin-up ===================='
         write(*,'(5a)') ' run   : ', time_to_string(run%cfg%start_time), ' -> ',                &
                         time_to_string(run%cfg%end_time)
         write(*,'(a,f0.2,a,i0)') ' span  : ', years_between(run%now, run%cfg%end_time),         &
                         ' yr  steps/yr~', run%steps_per_year
      end if

      !----- 3. The STATE (restart) stream + the run's PFT provenance table. -------------------!
      if (run%cfg%state_write_state) then
         call ensure_output_dir(trim(run%cfg%state_output_dir))
         call write_pft_params_csv(run%cfg, trim(run%cfg%state_output_dir)//'/'//                   &
                                   trim(run%cfg%state_output_prefix)//'_pft_parameters.csv')
      end if

      !----- 3b. DIAGNOSTIC output ([output].enabled): the netCDF-free manager (registry +        !
      !          integrator buffers). The per-step tick stages closed periods; the step drains them.!
      if (run%cfg%output%enabled) then
         call ensure_output_dir(trim(run%cfg%output%dir))
         call manager_setup(run%out_files, run%cfg)
         !----- Give the DERIVED soil diagnostics the SAME retention curve the fast loop           !
         !      integrates on, rather than a second derivation from the TOML.  --------------------!
         if (run%cfg%fast_biophysics_on) call manager_set_soil_params(run%out_files, run%poly%fast_ctx%col_config%soil)
         if (len_trim(run%cfg%output%io_config) > 0)                                             &
            call apply_io_overrides(run%out_files, trim(run%cfg%output%io_config), run%verbose)
         call manager_finalize(run%out_files)
         call manager_alloc_buffers(run%out_files, run%poly%out_bufs)
         call activate_site_diag(run%out_files, run%poly%site)
         if (run%verbose) write(*,'(a)') ' output: diagnostic aggregation ON ([output])'
      end if

      if (run%verbose) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         call print_summary(run%poly%site, 'start')
      end if

      run%is_open = .true.
      ok = .true.
   end subroutine driver_open

   !---------------------------------------------------------------------------------------!
   ! driver_step -- ONE slow step: advance the calendar, run the coupled stepper (which sub-steps  !
   ! the fast loop inside it), drain the FAST diagnostic tier, tick the slower tiers, and handle    !
   ! the year roll-over (summary, NaN guard, state checkpoint).                                     !
   !---------------------------------------------------------------------------------------!
   subroutine driver_step(run, status)
      type(meds_run_t), intent(inout) :: run
      integer(ik),      intent(out)   :: status
      logical          :: is_new_month, is_new_year
      character(len=19):: datestr

      status = DRIVER_OK
      if (driver_done(run)) then ; status = DRIVER_FINISHED ; return ; end if

      run%prev  = run%now
      run%now   = time_advance_days(run%prev, run%step_days)
      run%istep = run%istep + 1_ik

      is_new_year  = run%now%year  /= run%prev%year
      is_new_month = is_new_year .or. (run%now%month /= run%prev%month)

      !----- I/O before the compute phase: the step's forcing is loaded here, so nothing below reads  !
      !      a file (MEDS_POLYGON_RUNTIME_PLAN.md §4, R1). A no-op unless a new archive month starts.  !
      if (run%cfg%fast_biophysics_on .and. run%cfg%forcing%forcing_on) call met_prefetch(run%met_src, run%prev)

      call polygon_step(run%cfg, run%met_src, run%out_files, run%poly, run%prev, run%now,           &
                        run%step_days, is_new_month, is_new_year, status)

      if (is_new_year) then
         run%iyear = run%iyear + 1_ik
         if (run%verbose .and. (mod(run%iyear, 5_ik) == 0_ik .or. run%iyear == 1_ik)) then
            datestr = time_to_string(run%now)
            call print_summary(run%poly%site, datestr(1:10))
         end if
      end if
      !----- A failed step (NaN, impossible soil carbon) still writes the output up to it. ---------!
      if (status /= DRIVER_OK) then
         if (run%cfg%output%enabled) call output_serialize_pending(run%out_files, run%poly%out_bufs)
         return
      end if

      if (is_new_month) call driver_io_phase(run, is_new_year)
   end subroutine driver_step

   !---------------------------------------------------------------------------------------!
   ! driver_io_phase -- the file work of a closed month, done outside the compute phase so a step !
   ! makes no netCDF call (MEDS_POLYGON_RUNTIME_PLAN.md §4, R1): write the queued output records, !
   ! then the yearly checkpoint. The next month's forcing is loaded by met_prefetch before the     !
   ! next step.                                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine driver_io_phase(run, is_new_year)
      type(meds_run_t), intent(inout) :: run
      logical,          intent(in)    :: is_new_year
      if (run%cfg%output%enabled) call output_serialize_pending(run%out_files, run%poly%out_bufs)
      if (is_new_year) then
         if (run%cfg%state_write_state .and. mod(run%iyear, run%cfg%state_interval_years_cfg) == 0_ik) &
            call state_write_state(run%poly%site, run%cfg, trim(run%cfg%state_output_dir),                  &
                                trim(run%cfg%state_output_prefix), run%now,                         &
                                run%poly%restructure_pending, run%poly%restructure_new_year)
      end if
   end subroutine driver_io_phase

   !---------------------------------------------------------------------------------------!
   ! driver_finalize -- terminal checkpoint, summary, the two whole-column budget reports, the    !
   ! slow ledger, and close the output streams + met reader. Safe to call once.                   !
   !---------------------------------------------------------------------------------------!
   subroutine driver_finalize(run, status)
      type(meds_run_t), intent(inout) :: run
      integer(ik), optional, intent(out) :: status
      real(wp) :: a1
      integer(ik) :: st

      st = DRIVER_OK
      !----- Always checkpoint the true terminal state so a restart resumes exactly here. --------!
      if (run%cfg%state_write_state)                                                                &
         call state_write_state(run%poly%site, run%cfg, trim(run%cfg%state_output_dir),                     &
                             trim(run%cfg%state_output_prefix), run%now,                            &
                             run%poly%restructure_pending, run%poly%restructure_new_year)

      if (run%verbose) call print_summary(run%poly%site, 'final')
      a1 = total_area(run%poly%site)
      if (run%verbose) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         write(*,'(a,f12.9,a,f12.9)') ' site area start=', run%poly%area_start, '  end=', a1
      end if
      if (abs(a1 - 1.0_wp) > 1.0e-5_wp) st = DRIVER_ERR_AREA

      if (run%verbose) call polygon_report(run%cfg, run%poly)

      if (run%cfg%output%enabled) call output_manager_close(run%out_files, run%poly%out_bufs, .true.)
      if (run%cfg%fast_biophysics_on .and. run%cfg%forcing%forcing_on) call met_close(run%met_src)
      run%is_open = .false.
      if (present(status)) status = st
   end subroutine driver_finalize

   !----- Release the site's allocatable components. Separate from finalize so a caller can still  !
   !      read the final state (the getters the C-API exposes) after the streams are closed. ------!
   subroutine driver_free(run)
      type(meds_run_t), intent(inout) :: run
      call site_free(run%poly%site)
   end subroutine driver_free

   !----- Apply the optional meds_io_config.toml per-variable override table to the manager's     !
   !      registry (§6.1 value grammar + unknown-key trap). Each `variables.<name> = <value>`      !
   !      entry: a bool force-enables / disables everywhere; a quoted "F D M Y" string replaces     !
   !      the stream mask. A name matching no registry variable is a hard error.                    !
   subroutine apply_io_overrides(files, tomlpath, verbose)
      type(output_files_t),   intent(inout) :: files
      character(len=*),       intent(in)    :: tomlpath
      logical,                intent(in)    :: verbose
      type(toml_table_t) :: tt
      logical            :: ok, found
      integer(ik)        :: i, mask, status
      character(len=256) :: raw, sval
      character(len=64)  :: name
      character(len=8)   :: bad
      call toml_parse_file(tomlpath, tt, ok)
      if (.not. ok) error stop 'meds_driver: cannot read [output].io_config = '//trim(tomlpath)
      do i = 1_ik, tt%n
         if (len_trim(tt%key(i)) <= 10) cycle
         if (tt%key(i)(1:10) /= 'variables.') cycle
         name = trim(tt%key(i)(11:))
         raw  = adjustl(tt%val(i))
         select case (trim(raw))
         case ('true', '.true.', 'True', 'TRUE')
            call apply_variable_override(files%reg, trim(name), OVR_TRUE, 0_ik, found)
         case ('false', '.false.', 'False', 'FALSE')
            call apply_variable_override(files%reg, trim(name), OVR_FALSE, 0_ik, found)
         case default
            if (raw(1:1) == '"') then                    ! quoted stream string "F D M Y"
               sval = raw(2:index(raw(2:), '"'))
               call parse_stream_mask(trim(sval), mask, status, bad)
               if (status /= 0_ik)                                                               &
                  error stop 'meds_driver: io_config unknown stream token "'//trim(bad)//        &
                             '" for variable '//trim(name)
               call apply_variable_override(files%reg, trim(name), OVR_MASK, mask, found)
            else
               error stop 'meds_driver: io_config bad value for '//trim(name)//                  &
                          ' (expected true|false or a quoted "F D M Y" string)'
            end if
         end select
         if (.not. found)                                                                        &
            error stop 'meds_driver: io_config variable "'//trim(name)//                         &
                       '" matches no registry variable (typo?)'
      end do
      call build_freq_index(files%reg)
      if (verbose) write(*,'(3a)') ' output: applied per-variable overrides from ', trim(tomlpath), ''
   end subroutine apply_io_overrides

   !----- Create the output directory if it does not exist (driver-level convenience).           !
   !      Filesystem access stays in the driver; the engine/library never touches it.            !
   subroutine ensure_output_dir(dir)
      character(len=*), intent(in) :: dir
      integer :: stat
      if (len_trim(dir) == 0 .or. trim(dir) == '.') return
      call execute_command_line('mkdir -p "'//trim(dir)//'"', wait=.true., exitstat=stat)
      if (stat /= 0) write(*,'(3a)') ' warning: could not create output dir "', trim(dir), '"'
   end subroutine ensure_output_dir

end module meds_driver
