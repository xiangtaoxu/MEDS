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
   use meds_config,                 only : meds_config_t, INIT_RESTART
   use meds_time,                   only : meds_time_t, time_to_string
   use meds_cohort_state_vars,      only : cohort_bad_field, cohort_field_finite, cohort_field_name
   use meds_patch_state_vars,       only : canopy_air_bad_field, canopy_air_field_name, soil_carbon_state_t, &
                                           litter_input_t
   use meds_polygon_state_vars,     only : polygon_t, reset_step_diagnostics
   use meds_stepper,                only : advance_one_step, advance_boundary
   use meds_slow_update,            only : refresh_canopy_depth, max_cohort_count
   use meds_slow_state,             only : acclimate_site_leaf_tables
   use meds_slow_ledger,            only : slow_ledger_t, slow_ledger_report
   use meds_biogeochem_params,      only : N_SOIL_POOL
   use meds_soil_carbon,            only : solve_soil_carbon_steady_state, build_litter_input,      &
                                           soil_carbon_bad_pool, soil_carbon_pool_name
   use meds_respiration,            only : decomposition_transfer_matrix
   use meds_met_types,          only : met_source_t, met_cursor_t
   use meds_met,             only : met_cursor_init
   use meds_output_reduce,      only : total_area
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use meds_budget_check,           only : budget_t, budget_report
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_params,          only : TIER_FAST
   use meds_output_integrate,       only : output_integrate, output_integrate_fast, close_tier
   implicit none
   private

   public :: polygon_run_t, polygon_prepare, polygon_step, polygon_report
   public :: ensure_output_dir
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_BOUNDS

   !----- Step and run status codes. OK/DONE are normal; the ERR codes are failures, returned as a  !
   !      status rather than an `error stop` so a library caller survives them.                   !
   integer(ik), parameter :: DRIVER_OK      = 0_ik   !< a slow step was taken
   integer(ik), parameter :: DRIVER_FINISHED= 1_ik   !< the calendar already reached end_time; nothing done
   integer(ik), parameter :: DRIVER_ERR_NAN = 2_ik   !< NaN in the state or the energy/water ledger
   integer(ik), parameter :: DRIVER_ERR_AREA= 3_ik   !< patch areas do not sum to 1 (finalize only)
   integer(ik), parameter :: DRIVER_ERR_BOUNDS = 4_ik !< a cohort or soil-carbon state is outside its bounds


   !----- What the run needs to step and check one polygon; its state is the polygon_t beside it. --!
   type :: polygon_run_t
      integer(ik)            :: cell = 1_ik       !< the cell in the forcing source's cell list
      character(len=40)      :: label = ''        !< how messages name it ('' for a site run)
      !----- Threads for the polygon's patch loop: [run].n_threads in a site run, 1 in a region. ----!
      integer(ik)            :: patch_threads = 1_ik
      !----- A bare-ground or census start's first window seeds the canopy air from the forcing; a  !
      !      restart's state file holds it. ----------------------------------------------------------!
      logical                :: seed_canopy_air = .false.
      type(met_cursor_t)     :: met_cur           !< set only if forcing_on
      type(budget_t)         :: energy_budget, water_budget   !< whole-column ledgers over the run
      !----- The per-layer face-closure residual over the run. Kept beside the two above              !
      !      because it answers the question they cannot: not "did the column conserve" but "did the  !
      !      heat move with the water". ---------------------------------------------------------------!
      type(budget_t)         :: face_budget
      !----- The SLOW step's ledger: the site store across each slow step, and the soil-carbon     !
      !      seam's run-worst checks (docs/science/soil_carbon.md §7). -------------------------------!
      type(slow_ledger_t)    :: slow_ledger
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
   ! polygon_prepare -- everything after the initial state is in place (init_polygon, which also   !
   ! seeds the fast loop's stores): the patch threads, the forcing cursor at the polygon's cell and  !
   ! location, the canopy-air depth and the soil-carbon spin-up. `keep_soil_carbon` is set when a   !
   ! restart already restored the pools, which the spin-up would silently discard.                 !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_prepare(cfg, met_src, poly, run, keep_soil_carbon,                           &
                              patch_threads, verbose)
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t),   intent(in)    :: met_src
      type(polygon_t),      intent(inout) :: poly
      type(polygon_run_t),  intent(inout) :: run
      logical,              intent(in)    :: keep_soil_carbon, verbose
      !----- Threads for the polygon's patch loop: a site run gives it [run].n_threads, a region, whose  !
      !      threads go to the polygons, gives each polygon one. ----------------------------------!
      integer(ik),          intent(in)    :: patch_threads
      type(slow_ledger_t) :: no_ledger   !< inactive: the first canopy-air depth books nothing

      run%energy_budget = budget_t() ; run%water_budget = budget_t() ; run%face_budget = budget_t()
      run%slow_ledger   = slow_ledger_t()
      run%most_patches = 0_ik ; run%most_cohorts = 0_ik
      run%most_patches_when = '' ; run%most_cohorts_when = ''
      run%fast_step_total = 0_ik ; run%status = DRIVER_OK
      run%area_start = total_area(poly)
      run%patch_threads = patch_threads
      run%seed_canopy_air = cfg%init%init_mode /= INIT_RESTART

      if (cfg%fast_biophysics_on) then
         if (cfg%forcing%forcing_on) then
            call met_cursor_init(met_src, run%met_cur, run%cell, poly%latitude_deg,                 &
                                 poly%longitude_deg, poly%elevation_m)
         end if
         if (verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if
      !----- The canopy-air depth follows the stand from the first step, whatever the stand came     !
      !      from (bare ground, census, restart); the slow step keeps it current after that. No      !
      !      ledger: a plain geometry update that keeps the canopy air's intensive state. -----------!
      call refresh_canopy_depth(poly, cfg, no_ledger)

      !----- Slow soil-carbon spin-up (opt-in): a successful STATE restart already carries the real !
      !      persisted pools; otherwise the pools start at the allocation-time zero unless          !
      !      spinup_steady requests the SASU steady-state solve.  ---------------------------------!
      if (cfg%soil_carbon%soil_carbon_on .and. cfg%soil_carbon%spinup_steady .and. .not. keep_soil_carbon)    &
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
      integer(ik) :: n_cohort_max

      status = DRIVER_OK

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
      !      (4)). Advance each site's growth-temperature running mean from the daily mean the       !
      !      PREVIOUS step's fast loop accumulated, then acclimate its leaf table. Taking last step's !
      !      mean is not a lag to apologise for: the quantity being tracked IS the preceding weeks'  !
      !      temperature, and one slow step is a thirtieth of the window. -------------------------!
      if (cfg%leaf_physiology%thermal_acclimation) call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, step_days)

      !----- The fast loop stages sub-daily samples into the polygon's share of the run's files. --!
      call stepper(out_bufs)

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

      !----- STATE GUARD, every step: the run's energy and water ledgers, whose running residual   !
      !      turns NaN the step any flux or store in them does; each cohort's state                  !
      !      (cohort_bad_field: density, carbon pools, SLA, water, tissue energies); each patch's    !
      !      canopy air (canopy_air_bad_field); and each CENTURY soil-carbon pool. A non-finite      !
      !      value is DRIVER_ERR_NAN and an out-of-bounds one DRIVER_ERR_BOUNDS: a STATUS, not an    !
      !      `error stop`, so a library caller survives it. An operator keeps these in bounds by     !
      !      construction, so a violation is a bug; stopping on the step it happens, naming the     !
      !      field, is what makes it findable.                                                       !
      !                                                                                          !
      !      The ledger checks the soil-carbon predicate at every phase boundary and can say WHICH  !
      !      operator broke it, which is far more useful -- but the ledger is a diagnostic with an  !
      !      off switch ([run].slow_ledger_on), and a safety assertion that a diagnostic flag can  !
      !      disable is a trap. So the guard itself lives here, unconditional, and the ledger's     !
      !      copy is the attribution.                                                              !
      !      Once per slow step over a few hundred cohorts x 10 fields; the cost does not register. !
      if (ledger_not_finite()) then
         write(*,'(a)') ' ERROR: a NaN reached the energy/water ledger.'
         if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
         write(*,'(2a)') '        in the step ending ', time_to_string(now)
         status = DRIVER_ERR_NAN
         run%status = status
         return
      end if
      block
         integer(ik) :: ii, kk
         do ii = 1_ik, poly%cohort%n
            kk = cohort_bad_field(poly%cohort, ii)
            if (kk /= 0_ik) then
               if (cohort_field_finite(poly%cohort, ii, kk)) then
                  write(*,'(a)') ' ERROR: a cohort''s state is outside its physical bounds.'
                  status = DRIVER_ERR_BOUNDS
               else
                  write(*,'(a)') ' ERROR: a cohort''s state is not finite.'
                  status = DRIVER_ERR_NAN
               end if
               if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
               write(*,'(3a,i0,a,i0,a,i0)') '        field ', trim(cohort_field_name(kk)),          &
                     ', cohort id ', poly%cohort%global_id(ii), ' (pft ',                         &
                     poly%cohort%pft(ii), ') in patch ', poly%cohort%owner_patch(ii)
               write(*,'(2a)')  '        date ', time_to_string(now)
               if (status == DRIVER_ERR_BOUNDS) then
                  write(*,'(a)') '        Density, carbon and water are amounts: never negative, never'
                  write(*,'(a)') '        infinite. Every operator keeps them so by construction.'
               end if
               run%status = status
               return
            end if
         end do
         if (allocated(poly%patch%cas)) then
            do ii = 1_ik, poly%patch%n
               kk = canopy_air_bad_field(poly%patch%cas(ii))
               if (kk /= 0_ik) then
                  write(*,'(a)') ' ERROR: a patch''s canopy air is not finite.'
                  if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
                  write(*,'(3a,i0,a,i0)') '        field ', trim(canopy_air_field_name(kk)), ', patch ', ii,  &
                        ' of ', poly%patch%n
                  write(*,'(2a)')  '        date ', time_to_string(now)
                  status = DRIVER_ERR_NAN
                  run%status = status
                  return
               end if
            end do
         end if
         if (cfg%soil_carbon%soil_carbon_on) then
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
      !      stages its sub-steps when the output reads them (bufs%fast_on). ------------------------!
      subroutine stepper(fast_bufs)
         type(output_buffers_t), intent(inout) :: fast_bufs
         call advance_one_step(poly, cfg, met_src=met_src, met_cur=run%met_cur,                           &
                               step_start=prev, out_bufs=fast_bufs,                                 &
                               run_energy_budget=run%energy_budget,                                 &
                               run_water_budget=run%water_budget,                                   &
                               run_face_budget=run%face_budget,                                     &
                               slow_ledger=run%slow_ledger, patch_threads=run%patch_threads,       &
                               seed_canopy_air=run%seed_canopy_air)
      end subroutine stepper

      !----- FAST tier: replay the sub-step samples the fast loop staged in bufs%fast(:), closing   !
      !      the tier every fast_interval_steps sub-steps. Then the diagnostic tick of the slower   !
      !      tiers, whose snapshots of fast-loop values read the last staged sub-step. Closed       !
      !      records are queued; the I/O phase writes them.  ---------------------------------------!
      subroutine tick_output(files, bufs)
         type(output_files_t),  intent(in)    :: files
         type(output_buffers_t), intent(inout) :: bufs
         integer(ik) :: isub
         if (bufs%fast_ready .and. files%reg%nidx(TIER_FAST) > 0_ik) then
            do isub = 1_ik, bufs%n_fast_sub
               call output_integrate_fast(files, bufs, isub, cfg%dt_fast)
               run%fast_step_total = run%fast_step_total + 1_ik
               if (mod(run%fast_step_total, max(files%fast_interval_steps, 1_ik)) == 0_ik)          &
                  call close_tier(files, bufs, TIER_FAST)
            end do
         end if
         call output_integrate(files, bufs, poly, prev, cfg%dt_slow, is_new_day, is_new_month,           &
                               is_new_year)
         bufs%fast_ready = .false.
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
      !----- Fusion targets the tolerance ceilings overrode. Silent when every restructuring      !
      !      reached them, which is the usual case.  -------------------------------------------------!
      if (cfg%demography%patch_fusion_on .and. run%most_patches > cfg%demography%max_patch)               &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: patch fusion left ', run%most_patches, ' patches on ',         &
               trim(run%most_patches_when), ' (max_patch = ', cfg%demography%max_patch,                            &
               '): patch_light_tol_max kept dissimilar patches apart'
      if (cfg%demography%cohort_fusion_on .and. run%most_cohorts > cfg%demography%max_cohort)             &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: cohort fusion left ', run%most_cohorts, ' cohorts in one patch on ',  &
               trim(run%most_cohorts_when), ' (max_cohort = ', cfg%demography%max_cohort,                          &
               '): cohort_size_tol_max kept dissimilar cohorts apart'
      !----- The soil-carbon seam, then the SLOW tier's ledger, over the window the two above      !
      !      cannot see (docs/science/order_of_processes.md §4). ---------------------------------!
      call slow_ledger_report(run%slow_ledger)
   end subroutine polygon_report

   !----- SASU steady-state soil-carbon spin-up from the configured climatological xi + litter    !
   !      (docs/science/soil_carbon.md §8, eq. (7)). ---------------------------------------------!
   subroutine soil_carbon_steady(poly, cfg, verbose)
      type(polygon_t),      intent(inout) :: poly
      type(meds_config_t),  intent(in)    :: cfg
      logical,              intent(in)    :: verbose
      type(soil_carbon_state_t) :: pools0, pools_ss
      type(litter_input_t) :: lit
      real(wp)    :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp)    :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL), lignin_bar(2)
      integer(ik) :: ipp
      pools0 = soil_carbon_state_t()                    ! zero-carbon reference (f_lignin=0)
      call decomposition_transfer_matrix(pools0, cfg%soil_carbon, a_mat, k_diag, er)
      xi_bar = cfg%soil_carbon%spinup_xi
      lit%labile_grnd = cfg%soil_carbon%spinup_labile_grnd
      lit%labile_soil = cfg%soil_carbon%spinup_labile_soil
      lit%struct_grnd = cfg%soil_carbon%spinup_struct_grnd
      lit%struct_soil = cfg%soil_carbon%spinup_struct_soil
      call build_litter_input(lit, u_bar, lignin_bar)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, pools_ss)
      do ipp = 1_ik, poly%patch%n
         poly%patch%soil_carbon(ipp) = pools_ss
      end do
      if (verbose) write(*,'(a,f8.3,a)') ' soilc : steady-state spin-up, total = ',               &
         pools_ss%fast_grnd_carbon + pools_ss%fast_soil_carbon + pools_ss%struct_grnd_carbon +   &
         pools_ss%struct_soil_carbon + pools_ss%microbial_carbon + pools_ss%slow_carbon +        &
         pools_ss%passive_carbon, ' kgC/m2'
   end subroutine soil_carbon_steady

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
