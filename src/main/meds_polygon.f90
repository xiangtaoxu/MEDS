! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_polygon -- one polygon of a run: everything that evolves at one forcing cell, and the one  !
! step that advances it (docs/science/order_of_processes.md §2).                                  !
!                                                                                          !
! A polygon is two records side by side in the run: its state (`polygon_t`, the stand and its     !
! sites' surfaces) and its run machinery (`polygon_run_t`: its fast context, forcing cursor,      !
! conservation ledgers and status). The machinery is not part of the state because its parts are  !
! defined in the fast and slow libraries, which use the state library. Every routine here takes   !
! both, and the run's one configuration: what differs between polygons, a cell's soil, is on its  !
! sites. What every polygon of a run shares -- the config, the forcing source, the output files   !
! -- is passed in read-only. Its output buffers are passed in too: the run holds them, one per    !
! polygon in a contiguous array, so the I/O phase takes that array whole rather than a component  !
! section. A run (meds_region) is one polygon or many, each stepped by polygon_step, so a polygon !
! of a region computes exactly what a one-polygon run at its cell computes.                       !
!                                                                                          !
! polygon_step does no file work: the caller loads the step's forcing before it (met_prefetch)     !
! and writes the queued output records after a month closes.                                       !
!==========================================================================================!
module meds_polygon
   use meds_kinds,                  only : wp, ik
   use meds_state_defaults,         only : is_unset_positive
   use meds_therm_lib,              only : temp_to_internal_energy, internal_energy_to_temp
   use meds_config,                 only : meds_config_t
   use meds_time,                   only : meds_time_t, time_to_string
   use meds_cohort_state_vars,      only : cohort_bad_field, cohort_field_name
   use meds_polygon_state_vars,     only : polygon_t, reset_step_diagnostics
   use meds_stepper,                only : advance_one_step, advance_boundary
   use meds_slow_dynamics,          only : refresh_canopy_depth
   use meds_fast_dynamics,          only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_fast_config,            only : acclimate_leaf_photo_table
   use meds_biogeochem_types,       only : litter_input_t, n_soil_pool, soilc_seam_t
   use meds_patch_state_vars,       only : soil_carbon_t
   use meds_soil_biogeochem,        only : assemble_transfer_matrix, solve_soil_carbon_steady_state, &
                                           build_litter_input, soil_carbon_bad_pool,             &
                                           soil_carbon_pool_name
   use meds_forcing_types,          only : met_source_t, met_cursor_t
   use meds_met_driver,             only : met_cursor_init
   use meds_diagnostic_reduce,      only : total_area, has_nan
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use meds_demography_cohort_fusefiss, only : max_cohort_count
   use meds_budget_check,           only : budget_t, budget_report
   use meds_slow_ledger,            only : slow_ledger_t, slow_ledger_report
   use meds_output_types,           only : output_files_t, output_buffers_t, MISSING_VALUE
   use meds_output_integrate,       only : output_integrate, output_integrate_fast, close_tier
   use meds_output_registry,        only : manager_setup, manager_finalize, manager_alloc_buffers, &
                                           manager_set_polygons,                                  &
                                           activate_site_diag, apply_variable_override,          &
                                           parse_stream_mask, build_freq_index, OVR_TRUE,         &
                                           OVR_FALSE, OVR_MASK
   use meds_constants,              only : n_soil_layer_max
   use meds_allometry_lib,          only : HEIGHT_CHAVE
   use meds_biophysics_opts,        only : SOIL_BC_AQUIFER
   use meds_forcing_types,          only : met_cells_t
   use meds_toml,                   only : toml_table_t, toml_parse_file
   implicit none
   private

   public :: polygon_run_t, polygon_prepare, polygon_step, polygon_report
   public :: open_output_files, attach_output, ensure_output_dir
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_BOUNDS
   public :: N_PATCH_INIT

   !----- Step and run status codes. OK/DONE are normal; the ERR codes are failures, returned as a  !
   !      status rather than an `error stop` so a library caller survives them.                   !
   integer(ik), parameter :: DRIVER_OK      = 0_ik   !< a slow step was taken
   integer(ik), parameter :: DRIVER_FINISHED= 1_ik   !< the calendar already reached end_time; nothing done
   integer(ik), parameter :: DRIVER_ERR_NAN = 2_ik   !< NaN in the state or the energy/water ledger
   integer(ik), parameter :: DRIVER_ERR_AREA= 3_ik   !< patch areas do not sum to 1 (finalize only)
   integer(ik), parameter :: DRIVER_ERR_BOUNDS = 4_ik !< a cohort or soil-carbon state is outside its bounds

   integer(ik), parameter :: N_PATCH_INIT = 6_ik     !< bare-ground patches when no census/restart

   !----- What the run needs to step and check one polygon; its state is the polygon_t beside it. --!
   type :: polygon_run_t
      integer(ik)            :: cell = 1_ik       !< the cell in the forcing source's cell list
      character(len=40)      :: label = ''        !< how messages name it ('' for a site run)
      type(fast_context_t)   :: fast_ctx          !< built only if fast_biophysics_on
      type(met_cursor_t)     :: met_cur           !< set only if forcing_on
      type(budget_t)         :: energy_budget, water_budget   !< whole-column ledgers over the run
      !----- The per-layer face-closure residual over the run. Kept beside the two above              !
      !      because it answers the question they cannot: not "did the column conserve" but "did the  !
      !      heat move with the water". ---------------------------------------------------------------!
      type(budget_t)         :: face_budget
      type(slow_ledger_t)    :: slow_ledger                   !< site store across each SLOW step
      !----- Worst soil-carbon SEAM gap over the run [kgC/m2]: |daily pool debit - the fast loop's   !
      !      own accumulated Rh|. Both ends read the same frozen pool and the same per-pool xi       !
      !      integral (docs/science/soil_carbon.md §7), so this is ~0 BY CONSTRUCTION; a nonzero     !
      !      value means the double-counting contract broke (a stale frozen copy, a mid-day pool     !
      !      write, a lost sub-step). Reported with the whole-column budgets rather than asserted    !
      !      fatally, matching them. ----------------------------------------------------------------!
      type(soilc_seam_t)     :: seam                        !< per-run worsts (rh gap, lignin, lambda)
      character(len=19)      :: worst_rh_seam_when = ''     !< the date the rh gap peaked
      integer(ik)            :: worst_rh_seam_npatch = 0_ik !< and the patch count then
      !----- The most patches, and the most cohorts in one patch, that any restructuring left, and    !
      !      when. Above max_patch or max_cohort, the fusion tolerance ceilings kept dissimilar       !
      !      patches or cohorts apart; the end-of-run report says so.                                  !
      integer(ik)            :: most_patches = 0_ik, most_cohorts = 0_ik
      character(len=19)      :: most_patches_when = '', most_cohorts_when = ''
      integer(ik)            :: fast_step_total = 0_ik      !< fast sub-steps replayed into the output
      real(wp)               :: area_start = 0.0_wp
      integer(ik)            :: status = DRIVER_OK          !< the last step's status
   end type polygon_run_t

contains

   !---------------------------------------------------------------------------------------!
   ! polygon_prepare -- everything after the initial community is in place: the polygon's own copy  !
   ! of the run's config, the fast context and                                                     !
   ! the fast reservoirs, the forcing cursor at the polygon's cell and location, an initial snow   !
   ! pack, and the soil-carbon spin-up. `keep_fast_state` / `keep_soil_carbon` are set when a       !
   ! restart already restored those states, which re-seeding would silently discard.                !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_prepare(cfg, met_src, poly, run, keep_fast_state, keep_soil_carbon,          &
                              patch_threads, verbose)
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t),   intent(in)    :: met_src
      type(polygon_t),      intent(inout) :: poly
      type(polygon_run_t),  intent(inout) :: run
      logical,              intent(in)    :: keep_fast_state, keep_soil_carbon, verbose
      !----- Threads for the polygon's patch loop: a site run gives it [run].n_threads, a region, whose  !
      !      threads go to the polygons, gives each polygon one. ----------------------------------!
      integer(ik),          intent(in)    :: patch_threads

      run%energy_budget = budget_t() ; run%water_budget = budget_t() ; run%face_budget = budget_t()
      run%slow_ledger   = slow_ledger_t()
      run%seam          = soilc_seam_t()
      run%worst_rh_seam_when = '' ; run%worst_rh_seam_npatch = 0_ik
      run%most_patches = 0_ik ; run%most_cohorts = 0_ik
      run%most_patches_when = '' ; run%most_cohorts_when = ''
      run%fast_step_total = 0_ik ; run%status = DRIVER_OK
      run%area_start = total_area(poly)

      if (cfg%fast_biophysics_on) then
         call build_fast_context(cfg, run%fast_ctx)
         run%fast_ctx%patch_threads = patch_threads
         if (cfg%forcing%forcing_on) then
            call met_cursor_init(met_src, run%met_cur, run%cell, poly%latitude_deg,                 &
                                 poly%longitude_deg, poly%elevation_m)
         end if
         !----- Skip the generic re-seed when a restart already restored the true evolved            !
         !      CAS/soil/snow state (docs/science/order_of_processes.md §1): overwriting it here     !
         !      would silently discard exactly what io_read_state just read back, and the run would  !
         !      jump at the restart.                                                                 !
         if (.not. keep_fast_state) call init_fast_reservoirs(poly, run%fast_ctx)
         if (cfg%snow_init_swe > 0.0_wp) call seed_snow(poly, cfg, verbose)
         if (verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if
      !----- The canopy-air depth follows the stand from the first step, whatever the stand came     !
      !      from (bare ground, census, restart); the slow step keeps it current after that. No      !
      !      ledger: a plain geometry update that keeps the canopy air's intensive state. -----------!
      call refresh_canopy_depth(poly, cfg)

      !----- Slow soil-carbon spin-up (opt-in): a successful STATE restart already carries the real !
      !      persisted pools; otherwise the pools start at the allocation-time zero unless          !
      !      spinup_steady requests the SASU steady-state solve.  ---------------------------------!
      if (cfg%soil_carbon_on .and. cfg%soil_carbon_spinup_steady .and. .not. keep_soil_carbon)    &
         call soil_carbon_steady(poly, cfg, verbose)

      run%slow_ledger%active = cfg%slow_on .and. cfg%slow_ledger_on
   end subroutine polygon_prepare

   !---------------------------------------------------------------------------------------!
   ! polygon_step -- ONE slow step of one polygon, from `prev` to `now`: the restructuring the last !
   ! step's boundary owes (if any), the growth-temperature mean, the coupled stepper (which         !
   ! sub-steps the fast loop inside it), the fast output tier's replay, the slower tiers' tick, and !
   ! the NaN and bounds guards. The step's forcing must already be loaded into `met_src`;           !
   ! closed output records are queued in the polygon's buffers, `out_bufs`.                         !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_step(cfg, met_src, out_files, out_bufs, poly, run, prev, now, step_days,      &
                           is_new_month, is_new_year, status)
      type(meds_config_t),    intent(in)    :: cfg
      type(met_source_t),     intent(in)    :: met_src
      type(output_files_t),   intent(in)    :: out_files
      type(output_buffers_t), intent(inout) :: out_bufs   !< this polygon's share of out_files
      type(polygon_t),        intent(inout) :: poly
      type(polygon_run_t),    intent(inout) :: run
      type(meds_time_t),     intent(in)    :: prev, now
      integer(ik),           intent(in)    :: step_days
      logical,               intent(in)    :: is_new_month, is_new_year
      integer(ik),           intent(out)   :: status
      logical  :: is_new_day
      real(wp) :: seam_prev
      integer(ik) :: n_cohort_max

      status = DRIVER_OK
      seam_prev = run%seam%worst_rh_gap       ! so the date below records the step the max MOVED on

      !----- The boundary the last step ended on (restructure_pending): the stand's monthly and       !
      !      yearly restructuring, before this step's fast loop sees the stand. ----------------------!
      if (poly%restructure_pending) then
         call advance_boundary(poly, cfg, .true., poly%restructure_new_year,                       &
                               slow_ledger=run%slow_ledger)
         poly%restructure_pending = .false. ; poly%restructure_new_year = .false.
         if (poly%patch%n > run%most_patches) then
            run%most_patches = poly%patch%n ; run%most_patches_when = time_to_string(prev)
         end if
         n_cohort_max = max_cohort_count(poly)
         if (n_cohort_max > run%most_cohorts) then
            run%most_cohorts = n_cohort_max ; run%most_cohorts_when = time_to_string(prev)
         end if
      end if

      !----- THERMAL ACCLIMATION (docs/science/leaf_gas_exchange.md §2 Thermal acclimation, eq.     !
      !      (4)). Advance the growth-temperature running mean from the daily mean the PREVIOUS     !
      !      step's fast loop accumulated, then refresh the leaf table. Taking last step's mean is  !
      !      not a lag to apologise for: the quantity being tracked IS the preceding weeks'         !
      !      temperature, and one slow step is a thirtieth of the window. A no-op when the flag is  !
      !      off. -----------------------------------------------------------------------------------!
      if (cfg%leaf_thermal_acclimation) call advance_growth_temperature(cfg, step_days, poly, run)

      !----- The fast loop stages sub-daily samples into the polygon's share of the run's files. --!
      call stepper(out_bufs)
      !----- Keep WHERE the worst gap happened, not just how big. A seam that is zero except on the !
      !      days a patch operator fires is telling you something quite different from one that      !
      !      drifts every day, and the date is what distinguishes them. The accumulation itself is   !
      !      done inside the biogeochem driver, which sees every patch; this only records the date   !
      !      on the step where the running maximum moved.  ------------------------------------------!
      if (run%seam%worst_rh_gap > seam_prev) then
         run%worst_rh_seam_when   = time_to_string(now)
         run%worst_rh_seam_npatch = poly%patch%n
      end if

      !----- Output: replay the staged FAST samples, then fold this step -- its fluxes from `prev`  !
      !      to `now` and its post-dynamics state at `now` -- into the slower tiers' windows that   !
      !      hold `prev`, closing each period `now` has left and queueing it for the I/O phase. ----!
      is_new_day = is_new_month .or. (now%day /= prev%day)
      if (out_files%enabled) call tick_output(out_files, out_bufs)

      !----- The step's diagnostics are read: zero them for the next window. If the step ended on a  !
      !      month boundary, the stand's restructuring is owed: it runs at the start of the next     !
      !      step, after the output has read this one, so the ending period's records describe one   !
      !      stand, and what the restructuring does -- its events, and the stand it leaves --        !
      !      belongs to the period that begins (docs/science/diagnostics.md §4). --------------------!
      call reset_step_diagnostics(poly)
      poly%restructure_pending  = is_new_month
      poly%restructure_new_year = is_new_year

      !----- NAN GUARD, every step: the stand, the tissue temperatures and the canopy air         !
      !      (has_nan), and the run's energy and water ledgers, whose running residual turns NaN    !
      !      the step any flux or store in them does. A NaN is a STATUS here, not an `error stop`:   !
      !      a library caller survives it.  ---------------------------------------------------------!
      if (has_nan(poly) .or. ledger_not_finite()) then
         write(*,'(a)') ' ERROR: a NaN reached the state or the energy/water ledger.'
         if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
         write(*,'(2a)') '        in the step ending ', time_to_string(now)
         status = DRIVER_ERR_NAN ; run%status = status ; return
      end if

      !----- BOUNDS GUARD, checked every step and NOT gated on the ledger: each cohort's state     !
      !      (cohort_bad_field: density, carbon pools, SLA, water) and each CENTURY soil-carbon     !
      !      pool. An operator keeps these in bounds by construction, so a violation is a bug;     !
      !      stopping on the step it happens, naming the field, is what makes it findable.         !
      !                                                                                          !
      !      The ledger checks the soil-carbon predicate at every phase boundary and can say WHICH  !
      !      operator broke it, which is far more useful -- but the ledger is a diagnostic with an  !
      !      off switch ([run].slow_ledger_on), and a safety assertion that a diagnostic flag can  !
      !      disable is a trap. So the guard itself lives here, unconditional, and the ledger's     !
      !      copy is the attribution.                                                              !
      !      Once per slow step over a few hundred cohorts x 10 fields; the cost does not register. !
      block
         integer(ik) :: ii, kk
         do ii = 1_ik, poly%cohort%n
            kk = cohort_bad_field(poly%cohort, ii)
            if (kk /= 0_ik) then
               write(*,'(a)') ' ERROR: a cohort''s state is outside its physical bounds.'
               if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
               write(*,'(3a,i0,a,i0,a,i0)') '        field ', trim(cohort_field_name(kk)),          &
                     ', cohort id ', poly%cohort%global_id(ii), ' (pft ',                         &
                     poly%cohort%pft(ii), ') in patch ', poly%cohort%owner_patch(ii)
               write(*,'(2a)')  '        date ', time_to_string(now)
               write(*,'(a)') '        Density, carbon and water are amounts: never negative, never'
               write(*,'(a)') '        infinite. Every operator keeps them so by construction.'
               status = DRIVER_ERR_BOUNDS ; run%status = status
               return
            end if
         end do
         if (cfg%soil_carbon_on) then
            do ii = 1_ik, poly%patch%n
               kk = soil_carbon_bad_pool(poly%patch%soil_carbon(ii))
               if (kk /= 0_ik) then
                  write(*,'(a)') ' ERROR: a CENTURY soil-carbon pool is physically impossible.'
                  if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
                  write(*,'(3a,i0,a,i0)') '        pool ', trim(soil_carbon_pool_name(kk)),        &
                        ', patch ', ii, ' of ', poly%patch%n
                  write(*,'(2a)')  '        date ', time_to_string(now)
                  write(*,'(a)') '        A carbon pool is a mass: it cannot be negative, and the'
                  write(*,'(a)') '        ceiling is ~500x the richest real soil. Conservation can'
                  write(*,'(a)') '        hold perfectly while this is true -- see the slow ledger.'
                  status = DRIVER_ERR_BOUNDS ; run%status = status
                  return
               end if
            end do
         end if
      end block
      run%status = status

   contains

      !----- The run's whole-column ledgers: one NaN flux or store poisons the running residual  !
      !      for good, so a non-finite residual means a NaN entered the energy or water budget. ---!
      logical function ledger_not_finite()
         ledger_not_finite = cfg%fast_biophysics_on .and.                                              &
                             .not. (ieee_is_finite(run%energy_budget%resid_sum) .and.             &
                                    ieee_is_finite(run%water_budget%resid_sum))
      end function ledger_not_finite

      !----- Everything is handed on; the fast loop uses the forcing only when forcing_on, and     !
      !      stages the fast output tier only then. -------------------------------------------------!
      subroutine stepper(fast_bufs)
         type(output_buffers_t), intent(inout) :: fast_bufs
         call advance_one_step(poly, cfg, run%fast_ctx, met_src=met_src, met_cur=run%met_cur,             &
                               step_start=prev, out_bufs=fast_bufs,                                 &
                               run_energy_budget=run%energy_budget,                                 &
                               run_water_budget=run%water_budget,                                   &
                               run_face_budget=run%face_budget,                                     &
                               slow_ledger=run%slow_ledger, seam=run%seam,                          &
                               latitude_deg=poly%latitude_deg)
      end subroutine stepper

      !----- FAST tier: replay the sub-step samples the fast loop staged in bufs%fast(:), closing   !
      !      the tier every fast_interval_steps sub-steps. Then the diagnostic tick of the slower   !
      !      tiers. Closed records are queued; the I/O phase writes them.  -------------------------!
      subroutine tick_output(files, bufs)
         type(output_files_t),  intent(in)    :: files
         type(output_buffers_t), intent(inout) :: bufs
         integer(ik) :: isub
         if (bufs%fast_ready) then
            do isub = 1_ik, bufs%n_fast_sub
               call output_integrate_fast(files, bufs, isub, cfg%dt_fast)
               run%fast_step_total = run%fast_step_total + 1_ik
               if (mod(run%fast_step_total, max(files%fast_interval_steps, 1_ik)) == 0_ik)          &
                  call close_tier(files, bufs, 1_ik)
            end do
            bufs%fast_ready = .false.
         end if
         call output_integrate(files, bufs, poly, prev, cfg%dt_slow, is_new_day, is_new_month,           &
                               is_new_year)
      end subroutine tick_output

   end subroutine polygon_step

   !---------------------------------------------------------------------------------------!
   ! polygon_report -- the end-of-run conservation reports of one polygon: the two whole-column    !
   ! budgets and the face closure, the soil-carbon seam and freeze, and the slow ledger.           !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_report(cfg, run)
      type(meds_config_t), intent(in) :: cfg
      type(polygon_run_t), intent(in) :: run
      !----- Whole-column conservation over the RUN: the signed cumulative residual is the number a !
      !      per-step tolerance cannot see, so it is reported whether or not any single step        !
      !      breached.  -----------------------------------------------------------------------------!
      if (cfg%fast_biophysics_on) then
         call budget_report(run%energy_budget, 'whole_energy', 'J/m2',  'W/m2')
         call budget_report(run%water_budget,  'whole_water',  'kg/m2', 'kg/m2/s')
         !----- The VERTICAL check, reported in the same place and spirit as the soil-carbon          !
         !      seam below: a number that should be machine-zero, printed whether or not it breached. !
         !      A nonzero value means soil enthalpy was advected on a mass flux the committed water    !
         !      never moved -- a misplacement between LAYERS, which both ledgers above sum away.       !
         write(*,'(a,es12.3,a,es12.3,a,i0)') ' faces[soil_layer_mass]  worst = ',                  &
               run%face_budget%worst, ' kg/m2   mean |resid| = ',                                   &
               run%face_budget%abs_sum / max(real(run%face_budget%n_check, wp), 1.0_wp),            &
               ' kg/m2   checks = ', run%face_budget%n_check
         if (run%face_budget%n_fail > 0_ik)                                                        &
            write(*,'(a,i0,a)') ' WARNING: ', run%face_budget%n_fail,                                &
               ' soil-water solves did not converge (see [fast].debug_error to make this fatal)'
         if (run%energy_budget%n_fail + run%water_budget%n_fail > 0_ik)                          &
            write(*,'(a,i0,a)') ' WARNING: ', run%energy_budget%n_fail + run%water_budget%n_fail,   &
               ' whole-column budget checks breached tolerance (see [fast].debug_error to make this fatal)'
      end if
      !----- The soil-carbon seam, in the same place and spirit as the two budgets above: a number  !
      !      that should be machine-zero, reported whether or not it ever breached.  ----------------!
      if (cfg%soil_carbon_on) then
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_rh]  worst |pool debit - fast Rh| = ',        &
               run%seam%worst_rh_gap, ' kgC/m2'
         if (run%seam%worst_rh_gap > 0.0_wp)                                                      &
            write(*,'(3a,i0,a)') '                       worst at ', run%worst_rh_seam_when,      &
                  ' with ', run%worst_rh_seam_npatch, ' patches'
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_lignin]  worst passive-tracer resid = ',      &
               run%seam%worst_lignin, ' kgC/m2/day'
         !----- LAMBDA: the fraction of a pool one slow step withdraws -- the number that says      !
         !      whether freezing the pool across the step is sound. Not a residual: it is expected   !
         !      to be small-but-nonzero (3.6e-3/day measured on the mature Ithaca stand), and it is a  !
         !      WARNING as it approaches 1, long before a pool could go negative.  ------------------!
         write(*,'(a,es12.3,2a)') ' freeze[soil_carbon]   worst lambda = dt.K.xi / pool  = ',     &
               run%seam%worst_lambda, ' [-]  pool ',                                              &
               trim(soil_carbon_pool_name(run%seam%lambda_pool))
         if (run%seam%worst_lambda > 0.1_wp)                                                      &
            write(*,'(a)') '   WARNING: a slow step withdraws >10% of a soil-carbon pool --'      &
                        // ' the frozen-pool approximation is degrading; shorten dt_slow.'
      end if
      !----- Fusion targets the tolerance ceilings overrode. Silent when every restructuring      !
      !      reached them, which is the usual case.  -------------------------------------------------!
      if (cfg%do_patch_fissfuse .and. cfg%max_patch > 0_ik .and. run%most_patches > cfg%max_patch)     &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: patch fusion left ', run%most_patches, ' patches on ',         &
               trim(run%most_patches_when), ' (max_patch = ', cfg%max_patch,                            &
               '): patch_light_tol_max kept dissimilar patches apart'
      if (cfg%do_cohort_fissfuse .and. cfg%max_cohort > 0_ik .and. run%most_cohorts > cfg%max_cohort)  &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: cohort fusion left ', run%most_cohorts, ' cohorts in one patch on ',  &
               trim(run%most_cohorts_when), ' (max_cohort = ', cfg%max_cohort,                          &
               '): cohort_size_tol_max kept dissimilar cohorts apart'
      !----- The SLOW tier's ledger, over the window the two above cannot see                      !
      !      (docs/science/order_of_processes.md §4). ---------------------------------------------!
      call slow_ledger_report(run%slow_ledger)
   end subroutine polygon_report

   !---------------------------------------------------------------------------------------!
   ! Advance each site's growth-temperature running mean and push it into the site's leaf table    !
   ! (docs/science/leaf_gas_exchange.md §2 Thermal acclimation). Exponential mean with weight      !
   ! dt/window; UNSEEDED (negative) adopts the first day's mean outright, so the acclimated        !
   ! optimum starts from the site's own climate rather than relaxing toward it from an arbitrary   !
   ! origin over the first month.                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_growth_temperature(cfg, step_days, poly, run)
      type(meds_config_t),  intent(in)    :: cfg
      integer(ik),          intent(in)    :: step_days
      type(polygon_t),      intent(inout) :: poly
      type(polygon_run_t),  intent(inout) :: run
      real(wp)    :: t_day, w
      integer(ik) :: is
      if (.not. allocated(run%fast_ctx%site_leaf_photo))                                             &
         allocate(run%fast_ctx%site_leaf_photo(poly%site%n), source=run%fast_ctx%col_config%leaf_photo)
      do is = 1_ik, poly%site%n
         if (poly%site%pheno_tair_n(is) < 1_ik) cycle      ! no fast sub-steps ran: nothing to average
         t_day = poly%site%pheno_tair_sum(is) / real(poly%site%pheno_tair_n(is), wp)
         associate (t_growth => poly%site%t_growth_avg(is))
            if (is_unset_positive(t_growth)) then
               t_growth = t_day
            else
               w = min(1.0_wp, real(step_days, wp) / max(cfg%acclim_window_days, real(step_days, wp)))
               t_growth = t_growth + w * (t_day - t_growth)
            end if
            call acclimate_leaf_photo_table(cfg, t_growth, run%fast_ctx%site_leaf_photo(is))
         end associate
      end do
   end subroutine advance_growth_temperature

   !----- Seed an initial snow pack (spin-up / test); [fast].snow_init_swe. ---------------------!
   subroutine seed_snow(poly, cfg, verbose)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      logical,              intent(in)    :: verbose
      integer(ik) :: ipp
      do ipp = 1_ik, poly%patch%n
         poly%patch%snow(ipp)%swe(1)         = cfg%snow_init_swe
         poly%patch%snow(ipp)%snow_energy(1) =                                                   &
            temp_to_internal_energy(0.0_wp, cfg%snow_init_swe, cfg%snow_init_temp, 0.0_wp)
         poly%patch%snow(ipp)%snow_depth(1)  = cfg%snow_init_swe / 250.0_wp
         poly%patch%snow(ipp)%nlayer         = 1_ik
         call internal_energy_to_temp(poly%patch%snow(ipp)%snow_energy(1),                       &
                                      cfg%snow_init_swe, 0.0_wp,                                 &
                                      poly%patch%snow(ipp)%snow_temp(1),                         &
                                      poly%patch%snow(ipp)%snow_fliq(1))
      end do
      if (verbose) write(*,'(a,f6.1,a)') ' snow  : seeded initial pack SWE = ', cfg%snow_init_swe, ' kg/m2'
   end subroutine seed_snow

   !----- SASU steady-state soil-carbon spin-up from the configured climatological xi + litter    !
   !      (docs/science/soil_carbon.md §8, eq. (7)). ---------------------------------------------!
   subroutine soil_carbon_steady(poly, cfg, verbose)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      logical,              intent(in)    :: verbose
      type(soil_carbon_t)  :: pools0, pools_ss
      type(litter_input_t) :: lit
      real(wp)    :: a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp)    :: xi_bar(n_soil_pool), u_bar(n_soil_pool), lignin_bar(2)
      integer(ik) :: ipp
      pools0 = soil_carbon_t()                          ! zero-carbon reference (f_lignin=0)
      call assemble_transfer_matrix(pools0, cfg%soil_carbon, a_mat, k_diag, er)
      xi_bar = cfg%soil_carbon_spinup_xi
      lit%labile_grnd = cfg%soil_carbon_spinup_labile_grnd
      lit%labile_soil = cfg%soil_carbon_spinup_labile_soil
      lit%struct_grnd = cfg%soil_carbon_spinup_struct_grnd
      lit%struct_soil = cfg%soil_carbon_spinup_struct_soil
      call build_litter_input(lit, u_bar, lignin_bar)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, cfg%soil_carbon, pools_ss)
      do ipp = 1_ik, poly%patch%n
         poly%patch%soil_carbon(ipp) = pools_ss
      end do
      if (verbose) write(*,'(a,f8.3,a)') ' soilc : steady-state spin-up, total = ',               &
         pools_ss%fast_grnd_carbon + pools_ss%fast_soil_carbon + pools_ss%struct_grnd_carbon +   &
         pools_ss%struct_soil_carbon + pools_ss%microbial_carbon + pools_ss%slow_carbon +        &
         pools_ss%passive_carbon, ' kgC/m2'
   end subroutine soil_carbon_steady

   !---------------------------------------------------------------------------------------!
   ! open_output_files -- lay out one output file set: the registry from the config, the   !
   ! per-variable overrides of [output].io_config, and the polygon axis with each          !
   ! polygon's cell and site count. Each polygon that writes into the set then calls       !
   ! attach_output, which gives the files that polygon's own sites.                        !
   !---------------------------------------------------------------------------------------!
   subroutine open_output_files(cfg, files, verbose, cells, polygon_id, n_site)
      type(meds_config_t),  intent(in)    :: cfg
      type(output_files_t), intent(inout) :: files
      logical,              intent(in)    :: verbose
      type(met_cells_t),    intent(in)    :: cells          !< the run's cells, one per polygon ...
      integer(ik),          intent(in)    :: polygon_id(:)  !< ... their polygon ids ...
      integer(ik),          intent(in)    :: n_site(:)      !< ... and their site counts
      call manager_setup(files, cfg)
      if (len_trim(cfg%output%io_config) > 0)                                                     &
         call apply_io_overrides(files, trim(cfg%output%io_config), verbose)
      call manager_set_polygons(files, polygon_id, n_site, cells%lat, cells%lon, cells%row, cells%col, &
                                any(cfg%pft%allom(1:cfg%pft%n)%height_form == HEIGHT_CHAVE))
      call manager_finalize(files)
   end subroutine open_output_files

   !----- A polygon's share of an output file set: its buffers, which read the soil diagnostics   !
   !      when its fast loop runs, the diagnostic blocks it must accumulate for the set's         !
   !      variables, and its row of the per-polygon variables: each site's area, terrain, layers  !
   !      and soil column, and the height allometry's site index its cell's map gave it. ---------!
   subroutine attach_output(cfg, files, bufs, poly, run)
      type(meds_config_t),    intent(in)    :: cfg
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs
      type(polygon_t),        intent(inout) :: poly
      type(polygon_run_t),    intent(in)    :: run
      integer(ik) :: p, is
      call manager_alloc_buffers(files, bufs)
      bufs%diag%soil_ready = cfg%fast_biophysics_on
      call activate_site_diag(files, poly)
      p = run%cell
      do is = 1_ik, poly%site%n
         associate (surf => poly%site%surface(is))
            files%site_area(is, p)         = poly%site%area(is)
            files%site_elevation(is, p)    = surf%elevation_m
            files%site_slope(is, p)        = surf%slope_deg
            files%site_aspect(is, p)       = surf%aspect_deg
            files%site_n_soil(is, p)       = surf%soil%n_active
            files%site_soil_z(:, is, p)    = surf%soil%z_node
            files%site_soil_face(:, is, p) = surf%soil%soil_layer_z
            files%site_soil_depth(is, p)   = surf%soil_column%depth
            files%site_bottom_bc(is, p)    = surf%water_bottom_bc
            files%site_drain_factor(is, p) = merge(MISSING_VALUE, surf%drain_factor,                &
                                                   surf%water_bottom_bc == SOIL_BC_AQUIFER)
            files%site_deep_temp(is, p)    = surf%deep_temp
         end associate
      end do
      !----- E is the polygon's, the same in every PFT's record. -----------------------------!
      if (allocated(files%polygon_chave_e)) files%polygon_chave_e(p) = poly%cohort%pft_table%allom(1)%chave_e
   end subroutine attach_output

   !----- Apply the optional meds_io_config.toml per-variable override table to the manager's   !
   !      registry (value grammar and unknown-name trap: docs/science/diagnostics.md §5). Each  !
   !      `variables.<name> = <value>` entry: a bool force-enables / disables everywhere; a     !
   !      quoted "F D M Y" string replaces the stream mask. A name matching no registry         !
   !      variable is a hard error.                                                             !
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
      if (.not. ok) error stop 'meds_polygon: cannot read [output].io_config = '//trim(tomlpath)
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
                  error stop 'meds_polygon: io_config unknown stream token "'//trim(bad)//        &
                             '" for variable '//trim(name)
               call apply_variable_override(files%reg, trim(name), OVR_MASK, mask, found)
            else
               error stop 'meds_polygon: io_config bad value for '//trim(name)//                  &
                          ' (expected true|false or a quoted "F D M Y" string)'
            end if
         end select
         if (.not. found)                                                                        &
            error stop 'meds_polygon: io_config variable "'//trim(name)//                         &
                       '" matches no registry variable (typo?)'
      end do
      call build_freq_index(files%reg)
      if (verbose) write(*,'(3a)') ' output: applied per-variable overrides from ', trim(tomlpath), ''
   end subroutine apply_io_overrides

   !----- Create the output directory if it does not exist (a run-level convenience). The engine   !
   !      and the kernel libraries never touch the filesystem.                                      !
   subroutine ensure_output_dir(dir)
      character(len=*), intent(in) :: dir
      integer :: stat
      if (len_trim(dir) == 0 .or. trim(dir) == '.') return
      call execute_command_line('mkdir -p "'//trim(dir)//'"', wait=.true., exitstat=stat)
      if (stat /= 0) write(*,'(3a)') ' warning: could not create output dir "', trim(dir), '"'
   end subroutine ensure_output_dir

end module meds_polygon
