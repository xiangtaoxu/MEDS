!==========================================================================================!
! meds_driver -- the coupled model as an OPEN / STEP / FINALIZE object, so something other      !
! than the `meds_main` program can drive it.                                                    !
!                                                                                          !
! This is `meds_main`'s body, lifted verbatim. That program was 400 lines of driver logic with   !
! no seam in it: configuration, initial community, fast context, met reader, output manager,     !
! the calendar loop and the closing conservation reports were all statements in one PROGRAM, so  !
! the ONLY way to run MEDS was to exec the binary. `meds_main` is now a thin shell over this      !
! module and the C-API shim `meds_capi_run` is a second caller -- which is what lets              !
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
   use meds_therm_lib,              only : temp_to_internal_energy, internal_energy_to_temp
   use meds_constants,              only : day_sec, yr_day
   use meds_config,                 only : meds_config_t, INIT_CENSUS, INIT_RESTART
   use meds_time,                   only : meds_time_t, time_lt, time_advance_days,            &
                                           time_to_string, years_between
   use meds_config_io,              only : load_meds_config, write_pft_params_csv
   use meds_site_state_types,       only : site_t, site_free
   use meds_demography_update,      only : update_overtopping_lai
   use meds_init,                   only : init_bare_ground, init_from_census
   use meds_stepper,                only : advance_one_step
   use meds_vegetation_dynamics,    only : advance_plant_traits
   use meds_fast_dynamics,          only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_biogeochem_types,       only : litter_input_t, n_soil_pool
   use meds_column_state_types,     only : soil_carbon_t
   use meds_soil_biogeochem,        only : assemble_transfer_matrix, solve_soil_carbon_steady_state, &
                                           build_litter_input
   use meds_forcing_types,          only : met_driver_t
   use meds_met_driver,             only : met_open, met_close
   use meds_diagnostic_reduce,      only : print_summary, total_area, has_nan
   use meds_budget_check,           only : budget_t, budget_report
   use meds_slow_ledger,            only : slow_ledger_t, slow_ledger_report
   use meds_soil_biogeochem,        only : soil_carbon_bad_pool, soil_carbon_pool_name
   use meds_io,                     only : io_write_state, io_read_state
   use meds_output_types,           only : output_manager_t
   use meds_output_registry,        only : manager_setup, manager_alloc_buffers,                &
                                           manager_set_soil_params, activate_site_diag,         &
                                           apply_variable_override, parse_stream_mask,          &
                                           build_freq_index, OVR_TRUE, OVR_FALSE, OVR_MASK
   use meds_output_integrate,       only : output_integrate, output_integrate_fast, close_tier
   use meds_output_manager,         only : output_serialize_pending, output_manager_close
   use meds_toml,                   only : toml_table_t, toml_parse_file
   implicit none
   private

   public :: meds_run_t, driver_open, driver_step, driver_finalize, driver_free, driver_done
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_SOILC

   !----- driver_step status codes. OK/DONE are normal; the two ERR codes are the conditions the   !
   !      program used to `error stop` on, returned instead so a library caller survives them.     !
   integer(ik), parameter :: DRIVER_OK      = 0_ik   !< a slow step was taken
   integer(ik), parameter :: DRIVER_FINISHED= 1_ik   !< the calendar already reached end_time; nothing done
   integer(ik), parameter :: DRIVER_ERR_NAN = 2_ik   !< NaN in the state at a year roll-over
   integer(ik), parameter :: DRIVER_ERR_AREA= 3_ik   !< patch areas no longer sum to 1 (finalize only)
   integer(ik), parameter :: DRIVER_ERR_SOILC = 4_ik !< a CENTURY pool is physically impossible

   integer(ik), parameter :: N_PATCH_INIT = 6_ik     !< bare-ground patches when no census/restart

   !----- Everything the calendar loop needs between steps. These were meds_main's locals; making  !
   !      them components is the whole extraction -- no state hides in module scope, so two runs    !
   !      can be open at once (which is exactly what the C-API's handle registry does).             !
   type :: meds_run_t
      type(meds_config_t)    :: cfg
      type(site_t)           :: site
      type(output_manager_t) :: mgr             !< built only if output.enabled
      type(fast_context_t)   :: fast_ctx        !< built only if fast_biophysics_on
      type(met_driver_t)     :: met_drv         !< opened only if forcing_on
      type(budget_t)         :: energy_budget, water_budget   !< whole-column ledgers over the run
      type(slow_ledger_t)    :: slow_ledger                   !< site store across each SLOW step
      type(meds_time_t)      :: now, prev
      integer(ik)            :: istep = 0_ik, iyear = 0_ik, fast_step_total = 0_ik
      integer(ik)            :: step_days = 1_ik, steps_per_year = 365_ik
      real(wp)               :: area_start = 0.0_wp
      !----- Worst soil-carbon SEAM gap over the run [kgC/m2]: |daily pool debit - the fast loop's   !
      !      own accumulated Rh|. Both ends read the same frozen pool and the same per-pool xi        !
      !      integral, so this is ~0 BY CONSTRUCTION; a nonzero value means the double-counting       !
      !      contract broke (a stale frozen copy, a mid-day pool write, a lost sub-step). Reported    !
      !      with the whole-column budgets rather than asserted fatally, matching them.  -------------!
      real(wp)               :: worst_rh_seam = 0.0_wp
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

      !----- 1. Read the run configuration. --------------------------------------------------!
      call load_meds_config(trim(path), run%cfg)   ! hard error if a file or required key is missing
      if (run%verbose) write(*,'(2a)') ' config: ', trim(path)

      !----- 2. Initial community per [init].init_mode (0 bare | 1 census | 2 restart); the file  !
      !         for the non-selected mode is ignored. Unusable input -> bare ground. A successful  !
      !         restart also recovers the calendar date to continue from.  ------------------------!
      init_ok = .false. ; fast_state_found = .false. ; run%now = run%cfg%start_time
      select case (run%cfg%init_mode)
      case (INIT_RESTART)
         call io_read_state(run%site, run%cfg, trim(run%cfg%init_restart_file), restart_time,   &
                            init_ok, fast_found=fast_state_found)
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
            end if
         else if (run%verbose) then
            write(*,'(3a)') ' init  : restart (mode 2) ', trim(run%cfg%init_restart_file),      &
                            ' not usable -- falling back to bare ground'
         end if
      case (INIT_CENSUS)
         call init_from_census(run%site, run%cfg, trim(run%cfg%init_census_file), init_ok)
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
         call init_bare_ground(run%site, run%cfg, N_PATCH_INIT)
         if (run%verbose) write(*,'(a)') ' init  : bare ground (mode 0)'
      end if
      run%area_start = total_area(run%site)

      !----- 2a. Census restart with plasticity ON: census cohorts sit in an established stand but !
      !          carry NO trait history, so acclimate their leaf traits to the current light        !
      !          environment INSTANTANEOUSLY (after the competition sweep). Bare ground legitimately !
      !          starts at top-of-canopy; a state restart already read the plastic traits from file. !
      if (run%cfg%trait_plasticity_on .and. run%cfg%init_mode == INIT_CENSUS .and. init_ok) then
         call update_overtopping_lai(run%site)
         call advance_plant_traits(run%site, run%cfg, run%cfg%dt_years, instantaneous=.true.)
      end if

      !----- 2b. Fast biophysics context (opt-in): build the static column config + seed the      !
      !          per-patch CAS/soil reservoirs ONCE.  ---------------------------------------------!
      if (run%cfg%fast_biophysics_on) then
         call build_fast_context(run%cfg, run%fast_ctx)
         if (run%cfg%forcing%forcing_on) then
            call met_open(run%met_drv, run%cfg%forcing)
            run%fast_ctx%zref = run%cfg%forcing%reference_height
            if (run%verbose) write(*,'(3a)') ' force : met forcing ON (', trim(run%cfg%forcing%path), ')'
         end if
         !----- Skip the generic re-seed when a restart already restored the true evolved CAS/soil/ !
         !      snow state (P5, MEDS_ED2_RK45_DESIGN.md): overwriting it here would silently discard !
         !      exactly what io_read_state just read back, reintroducing the restart discontinuity.  !
         if (.not. (run%cfg%init_mode == INIT_RESTART .and. init_ok .and. fast_state_found))     &
            call init_fast_reservoirs(run%site, run%fast_ctx)
         if (run%cfg%snow_init_swe > 0.0_wp) call seed_snow(run)
         if (run%verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if

      !----- 2c. Slow soil-carbon spin-up (opt-in): a successful STATE restart already carries the !
      !          real persisted pools (skip); otherwise the pools start at the allocation-time zero !
      !          unless spinup_steady requests the SASU steady-state solve.  ----------------------!
      if (run%cfg%soil_carbon_on .and. run%cfg%soil_carbon_spinup_steady .and.                   &
          .not. (run%cfg%init_mode == INIT_RESTART .and. init_ok)) call soil_carbon_steady(run)

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
      if (run%cfg%io_write_state) then
         call ensure_output_dir(trim(run%cfg%io_output_dir))
         call write_pft_params_csv(run%cfg, trim(run%cfg%io_output_dir)//'/'//                   &
                                   trim(run%cfg%io_output_prefix)//'_pft_parameters.csv')
      end if

      !----- 3b. DIAGNOSTIC output ([output].enabled): the netCDF-free manager (registry +        !
      !          integrator buffers). The per-step tick stages closed periods; the step drains them.!
      if (run%cfg%output%enabled) then
         call ensure_output_dir(trim(run%cfg%output%dir))
         call manager_setup(run%mgr, run%cfg)
         !----- Give the DERIVED soil diagnostics the SAME retention curve the fast loop           !
         !      integrates on, rather than a second derivation from the TOML.  --------------------!
         if (run%cfg%fast_biophysics_on) call manager_set_soil_params(run%mgr, run%fast_ctx%col_config%soil)
         if (len_trim(run%cfg%output%io_config) > 0)                                             &
            call apply_io_overrides(run%mgr, trim(run%cfg%output%io_config), run%verbose)
         call manager_alloc_buffers(run%mgr)
         call activate_site_diag(run%mgr, run%site)
         if (run%verbose) write(*,'(a)') ' output: diagnostic aggregation ON ([output])'
      end if

      if (run%verbose) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         call print_summary(run%site, 'start')
      end if

      run%slow_ledger%active = run%cfg%slow_on .and. run%cfg%slow_ledger_on
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
      logical          :: is_new_month, is_new_year, is_new_day
      integer(ik)      :: isub
      real(wp)         :: seam_gap
      character(len=19):: datestr

      status = DRIVER_OK
      if (driver_done(run)) then ; status = DRIVER_FINISHED ; return ; end if

      run%prev  = run%now
      run%now   = time_advance_days(run%prev, run%step_days)
      run%istep = run%istep + 1_ik

      is_new_year  = run%now%year  /= run%prev%year
      is_new_month = is_new_year .or. (run%now%month /= run%prev%month)

      !----- step_start is passed UNCONDITIONALLY (leaf phenology needs day-of-year every step);    !
      !      met_drv/mgr stay gated on forcing_on.  ------------------------------------------------!
      if (run%cfg%fast_biophysics_on .and. run%cfg%forcing%forcing_on) then
         call advance_one_step(run%site, run%cfg, is_new_month, is_new_year, run%fast_ctx,       &
                               met_drv=run%met_drv, step_start=run%prev, mgr=run%mgr,            &
                               run_energy_budget=run%energy_budget,                              &
                               run_water_budget=run%water_budget,                                &
                               slow_ledger=run%slow_ledger, worst_rh_seam_gap=seam_gap)
      else
         call advance_one_step(run%site, run%cfg, is_new_month, is_new_year, run%fast_ctx,       &
                               step_start=run%prev, run_energy_budget=run%energy_budget,         &
                               run_water_budget=run%water_budget, slow_ledger=run%slow_ledger,   &
                               worst_rh_seam_gap=seam_gap)
      end if
      run%worst_rh_seam = max(run%worst_rh_seam, seam_gap)

      !----- FAST (sub-daily) tier: replay the sub-step samples the fast loop staged in mgr%fast(:),!
      !      closing + draining the tier every fast_interval_steps sub-steps.  ---------------------!
      if (run%cfg%output%enabled .and. run%mgr%fast_ready) then
         do isub = 1_ik, run%mgr%n_fast_sub
            call output_integrate_fast(run%mgr, isub, run%cfg%dt_fast)
            run%fast_step_total = run%fast_step_total + 1_ik
            if (mod(run%fast_step_total, max(run%mgr%fast_interval_steps, 1_ik)) == 0_ik) then
               call close_tier(run%mgr, 1_ik)
               call output_serialize_pending(run%mgr)
            end if
         end do
         run%mgr%fast_ready = .false.
      end if

      !----- Diagnostic tick: fold this step's (post-dynamics) state into the active tiers and     !
      !      stage any closed period, then flush.  -------------------------------------------------!
      if (run%cfg%output%enabled) then
         is_new_day = is_new_month .or. (run%now%day /= run%prev%day)
         call output_integrate(run%mgr, run%site, run%now, run%cfg%dt_slow, is_new_day,          &
                               is_new_month, is_new_year)
         call output_serialize_pending(run%mgr)
      end if

      if (is_new_year) then
         run%iyear = run%iyear + 1_ik
         if (run%verbose .and. (mod(run%iyear, 5_ik) == 0_ik .or. run%iyear == 1_ik)) then
            datestr = time_to_string(run%now)
            call print_summary(run%site, datestr(1:10))
         end if
         !----- A NaN is a STATUS here, not an `error stop`: see the module header. ---------------!
         if (has_nan(run%site)) then ; status = DRIVER_ERR_NAN ; return ; end if
      end if

      !----- SOIL-CARBON PLAUSIBILITY, checked every step and NOT gated on the ledger. -----------!
      !                                                                                          !
      !      The ledger checks the same predicate at every phase boundary and can say WHICH        !
      !      operator broke it, which is far more useful -- but the ledger is a diagnostic with an  !
      !      off switch ([run].slow_ledger_on), and a safety assertion that a diagnostic flag can  !
      !      disable is a trap. So the guard itself lives here, unconditional, and the ledger's     !
      !      copy is the attribution.                                                              !
      !      Once per slow step over ~12 patches x 7 pools; the cost does not register.            !
      if (run%cfg%soil_carbon_on) then
         block
            integer(ik) :: ipp, kk
            do ipp = 1_ik, run%site%patch%n
               kk = soil_carbon_bad_pool(run%site%patch%soil_carbon(ipp))
               if (kk /= 0_ik) then
                  write(*,'(a)') ' ERROR: a CENTURY soil-carbon pool is physically impossible.'
                  write(*,'(3a,i0,a,i0)') '        pool ', trim(soil_carbon_pool_name(kk)),        &
                        ', patch ', ipp, ' of ', run%site%patch%n
                  write(*,'(2a)')  '        date ', time_to_string(run%now)
                  write(*,'(a)') '        A carbon pool is a mass: it cannot be negative, and the'
                  write(*,'(a)') '        ceiling is ~500x the richest real soil. Conservation can'
                  write(*,'(a)') '        hold perfectly while this is true -- see the slow ledger.'
                  status = DRIVER_ERR_SOILC
                  return
               end if
            end do
         end block
      end if

      if (is_new_year) then
         if (run%cfg%io_write_state .and. mod(run%iyear, run%cfg%io_state_interval_years) == 0_ik) &
            call io_write_state(run%site, run%cfg, trim(run%cfg%io_output_dir),                  &
                                trim(run%cfg%io_output_prefix), run%now)
      end if
   end subroutine driver_step

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
      if (run%cfg%io_write_state)                                                                &
         call io_write_state(run%site, run%cfg, trim(run%cfg%io_output_dir),                     &
                             trim(run%cfg%io_output_prefix), run%now)

      if (run%verbose) call print_summary(run%site, 'final')
      a1 = total_area(run%site)
      if (run%verbose) then
         write(*,'(a)') '-----------------------------------------------------------------------------'
         write(*,'(a,f12.9,a,f12.9)') ' site area start=', run%area_start, '  end=', a1
      end if
      if (abs(a1 - 1.0_wp) > 1.0e-5_wp) st = DRIVER_ERR_AREA

      !----- Whole-column conservation over the RUN: the signed cumulative residual is the number a !
      !      per-step tolerance cannot see, so it is reported whether or not any single step        !
      !      breached.  -----------------------------------------------------------------------------!
      if (run%cfg%fast_biophysics_on .and. run%verbose) then
         call budget_report(run%energy_budget, 'whole_energy', 'J/m2',  'W/m2')
         call budget_report(run%water_budget,  'whole_water',  'kg/m2', 'kg/m2/s')
         if (run%energy_budget%n_fail + run%water_budget%n_fail > 0_ik)                          &
            write(*,'(a,i0,a)') ' WARNING: ', run%energy_budget%n_fail + run%water_budget%n_fail, &
               ' whole-column budget checks breached tolerance (see [energy].debug_error to make this fatal)'
      end if
      !----- The soil-carbon seam, in the same place and spirit as the two budgets above: a number  !
      !      that should be machine-zero, reported whether or not it ever breached. Until now it was  !
      !      computed inside the slow driver and thrown away, because nothing asked for it.  ---------!
      if (run%cfg%soil_carbon_on .and. run%verbose)                                              &
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_rh]  worst |pool debit - fast Rh| = ',        &
               run%worst_rh_seam, ' kgC/m2  (0 by construction)'

      !----- The SLOW tier's ledger, over the window the two above cannot see (plan §10.2). -------!
      if (run%verbose) call slow_ledger_report(run%slow_ledger)

      if (run%cfg%output%enabled) call output_manager_close(run%mgr, .true.)
      if (run%cfg%fast_biophysics_on .and. run%cfg%forcing%forcing_on) call met_close(run%met_drv)
      run%is_open = .false.
      if (present(status)) status = st
   end subroutine driver_finalize

   !----- Release the site's allocatable components. Separate from finalize so a caller can still  !
   !      read the final state (the getters the C-API exposes) after the streams are closed. ------!
   subroutine driver_free(run)
      type(meds_run_t), intent(inout) :: run
      call site_free(run%site)
   end subroutine driver_free

   !----- Seed an initial snow pack (spin-up / test); [fast].snow_init_swe. ---------------------!
   subroutine seed_snow(run)
      type(meds_run_t), intent(inout) :: run
      integer(ik) :: ipp
      do ipp = 1_ik, run%site%patch%n
         run%site%patch%snow(ipp)%swe(1)         = run%cfg%snow_init_swe
         run%site%patch%snow(ipp)%snow_energy(1) =                                               &
            temp_to_internal_energy(0.0_wp, run%cfg%snow_init_swe, run%cfg%snow_init_temp, 0.0_wp)
         run%site%patch%snow(ipp)%snow_depth(1)  = run%cfg%snow_init_swe / 250.0_wp
         run%site%patch%snow(ipp)%nlayer         = 1_ik
         call internal_energy_to_temp(run%site%patch%snow(ipp)%snow_energy(1),                   &
                                      run%cfg%snow_init_swe, 0.0_wp,                             &
                                      run%site%patch%snow(ipp)%snow_temp(1),                     &
                                      run%site%patch%snow(ipp)%snow_fliq(1))
      end do
      if (run%verbose) write(*,'(a,f6.1,a)') ' snow  : seeded initial pack SWE = ',              &
                                             run%cfg%snow_init_swe, ' kg/m2'
   end subroutine seed_snow

   !----- SASU steady-state soil-carbon spin-up from the configured climatological xi + litter. --!
   subroutine soil_carbon_steady(run)
      type(meds_run_t), intent(inout) :: run
      type(soil_carbon_t)  :: pools0, pools_ss
      type(litter_input_t) :: lit
      real(wp)    :: a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp)    :: xi_bar(n_soil_pool), u_bar(n_soil_pool), lignin_bar(2)
      integer(ik) :: ipp
      pools0 = soil_carbon_t()                          ! zero-carbon reference (f_lignin=0)
      call assemble_transfer_matrix(pools0, run%cfg%soil_carbon, a_mat, k_diag, er)
      xi_bar = run%cfg%soil_carbon_spinup_xi
      lit%labile_grnd = run%cfg%soil_carbon_spinup_labile_grnd
      lit%labile_soil = run%cfg%soil_carbon_spinup_labile_soil
      lit%struct_grnd = run%cfg%soil_carbon_spinup_struct_grnd
      lit%struct_soil = run%cfg%soil_carbon_spinup_struct_soil
      call build_litter_input(lit, u_bar, lignin_bar)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, run%cfg%soil_carbon, pools_ss)
      do ipp = 1_ik, run%site%patch%n
         run%site%patch%soil_carbon(ipp) = pools_ss
      end do
      if (run%verbose) write(*,'(a,f8.3,a)') ' soilc : steady-state spin-up, total = ',          &
         pools_ss%fast_grnd_carbon + pools_ss%fast_soil_carbon + pools_ss%struct_grnd_carbon +   &
         pools_ss%struct_soil_carbon + pools_ss%microbial_carbon + pools_ss%slow_carbon +        &
         pools_ss%passive_carbon, ' kgC/m2'
   end subroutine soil_carbon_steady

   !----- Apply the optional meds_io_config.toml per-variable override table to the manager's     !
   !      registry (§6.1 value grammar + unknown-key trap). Each `variables.<name> = <value>`      !
   !      entry: a bool force-enables / disables everywhere; a quoted "F D M Y" string replaces     !
   !      the stream mask. A name matching no registry variable is a hard error.                    !
   subroutine apply_io_overrides(mgr, tomlpath, verbose)
      type(output_manager_t), intent(inout) :: mgr
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
            call apply_variable_override(mgr%reg, trim(name), OVR_TRUE, 0_ik, found)
         case ('false', '.false.', 'False', 'FALSE')
            call apply_variable_override(mgr%reg, trim(name), OVR_FALSE, 0_ik, found)
         case default
            if (raw(1:1) == '"') then                    ! quoted stream string "F D M Y"
               sval = raw(2:index(raw(2:), '"'))
               call parse_stream_mask(trim(sval), mask, status, bad)
               if (status /= 0_ik)                                                               &
                  error stop 'meds_driver: io_config unknown stream token "'//trim(bad)//        &
                             '" for variable '//trim(name)
               call apply_variable_override(mgr%reg, trim(name), OVR_MASK, mask, found)
            else
               error stop 'meds_driver: io_config bad value for '//trim(name)//                  &
                          ' (expected true|false or a quoted "F D M Y" string)'
            end if
         end select
         if (.not. found)                                                                        &
            error stop 'meds_driver: io_config variable "'//trim(name)//                         &
                       '" matches no registry variable (typo?)'
      end do
      call build_freq_index(mgr%reg)
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
