! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_polygon -- one polygon of a run: everything that evolves at one forcing cell, and the one  !
! step that advances it (MEDS_POLYGON_RUNTIME_PLAN.md §4, R2).                                     !
!                                                                                          !
! A polygon owns its site, fast context, forcing cursor, conservation ledgers and status. What     !
! every polygon of a run shares -- the config, the forcing source, the output files -- is passed   !
! in read-only. Its output buffers are passed in too: the run holds them, one per polygon in a     !
! contiguous array, so the I/O phase takes that array whole rather than a component section.      !
! A site run (meds_driver) is one polygon; a region (meds_region) is many, stepped by the same     !
! polygon_step, so a polygon of a region computes exactly what a site run at its cell computes.    !
!                                                                                          !
! polygon_step does no file work: the caller loads the step's forcing before it (met_prefetch)     !
! and writes the queued output records after a month closes.                                       !
!==========================================================================================!
module meds_polygon
   use meds_kinds,                  only : wp, ik
   use meds_therm_lib,              only : temp_to_internal_energy, internal_energy_to_temp
   use meds_config,                 only : meds_config_t
   use meds_time,                   only : meds_time_t, time_to_string
   use meds_site_state_types,       only : site_t, reset_step_diagnostics
   use meds_stepper,                only : advance_one_step, advance_boundary
   use meds_slow_dynamics,          only : refresh_canopy_depth
   use meds_fast_dynamics,          only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_fast_config,            only : acclimate_leaf_photo_table
   use meds_biogeochem_types,       only : litter_input_t, n_soil_pool, soilc_seam_t
   use meds_column_state_types,     only : soil_carbon_t
   use meds_soil_biogeochem,        only : assemble_transfer_matrix, solve_soil_carbon_steady_state, &
                                           build_litter_input, soil_carbon_bad_pool,             &
                                           soil_carbon_pool_name
   use meds_forcing_types,          only : met_source_t, met_cursor_t
   use meds_met_driver,             only : met_cursor_init
   use meds_diagnostic_reduce,      only : total_area, has_nan
   use meds_budget_check,           only : budget_t, budget_report
   use meds_slow_ledger,            only : slow_ledger_t, slow_ledger_report
   use meds_output_types,           only : output_files_t, output_buffers_t
   use meds_output_integrate,       only : output_integrate, output_integrate_fast, close_tier
   implicit none
   private

   public :: meds_polygon_t, polygon_prepare, polygon_step, polygon_report
   public :: DRIVER_OK, DRIVER_FINISHED, DRIVER_ERR_NAN, DRIVER_ERR_AREA, DRIVER_ERR_SOILC
   public :: N_PATCH_INIT

   !----- Step and run status codes. OK/DONE are normal; the ERR codes are failures, returned as a  !
   !      status rather than an `error stop` so a library caller survives them.                   !
   integer(ik), parameter :: DRIVER_OK      = 0_ik   !< a slow step was taken
   integer(ik), parameter :: DRIVER_FINISHED= 1_ik   !< the calendar already reached end_time; nothing done
   integer(ik), parameter :: DRIVER_ERR_NAN = 2_ik   !< NaN in the state at a year roll-over
   integer(ik), parameter :: DRIVER_ERR_AREA= 3_ik   !< patch areas no longer sum to 1 (finalize only)
   integer(ik), parameter :: DRIVER_ERR_SOILC = 4_ik !< a CENTURY pool is physically impossible

   integer(ik), parameter :: N_PATCH_INIT = 6_ik     !< bare-ground patches when no census/restart

   type :: meds_polygon_t
      integer(ik)            :: id = 0_ik         !< global id: the cell's row-major index (0 for a site run)
      integer(ik)            :: cell = 1_ik       !< the cell in the forcing source's cell list
      character(len=40)      :: label = ''        !< how messages name it ('' for a site run)
      type(site_t)           :: site
      type(fast_context_t)   :: fast_ctx          !< built only if fast_biophysics_on
      type(met_cursor_t)     :: met_cur           !< set only if forcing_on
      !----- A region's detail polygon also writes a full single-site file set of its own. ---------!
      type(output_files_t),  allocatable :: detail_files
      type(output_buffers_t), allocatable :: detail_bufs
      type(budget_t)         :: energy_budget, water_budget   !< whole-column ledgers over the run
      !----- The per-layer face-closure residual over the run (#189). Kept beside the two above       !
      !      because it answers the question they cannot: not "did the column conserve" but "did the  !
      !      heat move with the water". ---------------------------------------------------------------!
      type(budget_t)         :: face_budget
      type(slow_ledger_t)    :: slow_ledger                   !< site store across each SLOW step
      !----- Worst soil-carbon SEAM gap over the run [kgC/m2]: |daily pool debit - the fast loop's   !
      !      own accumulated Rh|. Both ends read the same frozen pool and the same per-pool xi        !
      !      integral, so this is ~0 BY CONSTRUCTION; a nonzero value means the double-counting       !
      !      contract broke (a stale frozen copy, a mid-day pool write, a lost sub-step). Reported    !
      !      with the whole-column budgets rather than asserted fatally, matching them.  -------------!
      type(soilc_seam_t)     :: seam                        !< per-run worsts (rh gap, lignin, lambda)
      character(len=19)      :: worst_rh_seam_when = ''     !< the date the rh gap peaked
      integer(ik)            :: worst_rh_seam_npatch = 0_ik !< and the patch count then
      integer(ik)            :: fast_step_total = 0_ik      !< fast sub-steps replayed into the output
      real(wp)               :: area_start = 0.0_wp
      integer(ik)            :: status = DRIVER_OK          !< the last step's status
      !----- The calendar restructuring the last step's end owes: set when a step ends on a month    !
      !      (and, with restructure_new_year, a year) boundary, run at the start of the next step,   !
      !      before its fast loop. Between the two sit the output tick, the I/O phase and any        !
      !      checkpoint, so a checkpoint holds the stand before the restructuring together with      !
      !      this flag, and a run resumed from it restructures first, exactly as the continuous run  !
      !      does. ----------------------------------------------------------------------------------!
      logical                :: restructure_pending  = .false.
      logical                :: restructure_new_year = .false.
   end type meds_polygon_t

contains

   !---------------------------------------------------------------------------------------!
   ! polygon_prepare -- everything after the initial community is in place: the fast context and   !
   ! the fast reservoirs, the forcing cursor at the polygon's cell and location, an initial snow   !
   ! pack, and the soil-carbon spin-up. `keep_fast_state` / `keep_soil_carbon` are set when a       !
   ! restart already restored those states, which re-seeding would silently discard.                !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_prepare(cfg, met_src, poly, latitude_deg, longitude_deg, utc_offset_h,     &
                              elevation_m, keep_fast_state, keep_soil_carbon, verbose)
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t),   intent(in)    :: met_src
      type(meds_polygon_t), intent(inout) :: poly
      real(wp),             intent(in)    :: latitude_deg, longitude_deg, utc_offset_h, elevation_m
      logical,              intent(in)    :: keep_fast_state, keep_soil_carbon, verbose

      poly%energy_budget = budget_t() ; poly%water_budget = budget_t() ; poly%face_budget = budget_t()
      poly%slow_ledger   = slow_ledger_t()
      poly%seam          = soilc_seam_t()
      poly%worst_rh_seam_when = '' ; poly%worst_rh_seam_npatch = 0_ik
      poly%fast_step_total = 0_ik ; poly%status = DRIVER_OK
      poly%area_start = total_area(poly%site)

      if (cfg%fast_biophysics_on) then
         call build_fast_context(cfg, poly%fast_ctx)
         if (cfg%forcing%forcing_on) then
            call met_cursor_init(met_src, poly%met_cur, poly%cell, latitude_deg, longitude_deg,     &
                                 utc_offset_h, elevation_m)
         end if
         !----- Skip the generic re-seed when a restart already restored the true evolved CAS/soil/ !
         !      snow state (P5, MEDS_ED2_RK45_DESIGN.md): overwriting it here would silently discard !
         !      exactly what io_read_state just read back, reintroducing the restart discontinuity.  !
         if (.not. keep_fast_state) call init_fast_reservoirs(poly%site, poly%fast_ctx)
         if (cfg%snow_init_swe > 0.0_wp) call seed_snow(cfg, poly, verbose)
         if (verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if
      !----- The canopy-air depth follows the stand from the first step, whatever the stand came    !
      !      from (bare ground, census, restart); the slow step keeps it current after that (#306).  !
      !      No ledger: a plain geometry update that keeps the canopy air's intensive state.  -------!
      call refresh_canopy_depth(poly%site, cfg)

      !----- Slow soil-carbon spin-up (opt-in): a successful STATE restart already carries the real !
      !      persisted pools; otherwise the pools start at the allocation-time zero unless          !
      !      spinup_steady requests the SASU steady-state solve.  ---------------------------------!
      if (cfg%soil_carbon_on .and. cfg%soil_carbon_spinup_steady .and. .not. keep_soil_carbon)    &
         call soil_carbon_steady(cfg, poly, verbose)

      poly%slow_ledger%active = cfg%slow_on .and. cfg%slow_ledger_on
   end subroutine polygon_prepare

   !---------------------------------------------------------------------------------------!
   ! polygon_step -- ONE slow step of one polygon, from `prev` to `now`: the restructuring the last !
   ! step's boundary owes (if any), the growth-temperature mean, the coupled stepper (which         !
   ! sub-steps the fast loop inside it), the fast output tier's replay, the slower tiers' tick, and !
   ! the NaN and soil-carbon guards. The step's forcing must already be loaded into `met_src`;      !
   ! closed output records are queued in the polygon's buffers, `out_bufs`.                         !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_step(cfg, met_src, out_files, out_bufs, poly, prev, now, step_days,           &
                           is_new_month, is_new_year, status)
      type(meds_config_t),    intent(in)    :: cfg
      type(met_source_t),     intent(in)    :: met_src
      type(output_files_t),   intent(in)    :: out_files
      type(output_buffers_t), intent(inout) :: out_bufs   !< this polygon's share of out_files
      type(meds_polygon_t),   intent(inout) :: poly
      type(meds_time_t),     intent(in)    :: prev, now
      integer(ik),           intent(in)    :: step_days
      logical,               intent(in)    :: is_new_month, is_new_year
      integer(ik),           intent(out)   :: status
      logical  :: is_new_day
      real(wp) :: seam_prev

      status = DRIVER_OK
      seam_prev = poly%seam%worst_rh_gap      ! so the date below records the step the max MOVED on

      !----- The boundary the last step ended on (restructure_pending): the stand's monthly and       !
      !      yearly restructuring, before this step's fast loop sees the stand. ----------------------!
      if (poly%restructure_pending) then
         call advance_boundary(poly%site, cfg, .true., poly%restructure_new_year,                  &
                               slow_ledger=poly%slow_ledger)
         poly%restructure_pending = .false. ; poly%restructure_new_year = .false.
      end if

      !----- THERMAL ACCLIMATION (#176). Advance the growth-temperature running mean from the      !
      !      daily mean the PREVIOUS step's fast loop accumulated, then refresh the leaf table.     !
      !      Taking last step's mean is not a lag to apologise for: the quantity being tracked IS   !
      !      the preceding weeks' temperature, and one slow step is a thirtieth of the window.      !
      !      A no-op when the flag is off, so the shipped path is untouched.  ------------------------!
      if (cfg%leaf_thermal_acclimation) call advance_growth_temperature(cfg, step_days, poly)

      !----- The fast loop stages sub-daily samples into the buffers that write the FAST tier: a    !
      !      detail polygon's own files, else the polygon's share of the run's files.  ----------------!
      if (allocated(poly%detail_bufs)) then
         call stepper(poly%detail_bufs)
      else
         call stepper(out_bufs)
      end if
      !----- Keep WHERE the worst gap happened, not just how big. A seam that is zero except on the !
      !      days a patch operator fires is telling you something quite different from one that      !
      !      drifts every day, and the date is what distinguishes them. The accumulation itself is   !
      !      done inside the biogeochem driver, which sees every patch; this only records the date   !
      !      on the step where the running maximum moved.  ------------------------------------------!
      if (poly%seam%worst_rh_gap > seam_prev) then
         poly%worst_rh_seam_when   = time_to_string(now)
         poly%worst_rh_seam_npatch = poly%site%patch%n
      end if

      !----- Output: replay the staged FAST samples, then fold this step -- its fluxes from `prev`  !
      !      to `now` and its post-dynamics state at `now` -- into the slower tiers' windows that   !
      !      hold `prev`, closing each period `now` has left and queueing it for the I/O phase. ----!
      is_new_day = is_new_month .or. (now%day /= prev%day)
      if (out_files%enabled) call tick_output(out_files, out_bufs)
      if (allocated(poly%detail_bufs)) call tick_output(poly%detail_files, poly%detail_bufs)

      !----- The step's diagnostics are read: zero them for the next window. If the step ended on a  !
      !      month boundary, the stand's restructuring is owed: it runs at the start of the next     !
      !      step, after the output has read this one, so the ending period's records describe one   !
      !      stand, and what the restructuring does -- its events, and the stand it leaves --        !
      !      belongs to the period that begins (docs/science/diagnostics.md §4). --------------------!
      call reset_step_diagnostics(poly%site)
      poly%restructure_pending  = is_new_month
      poly%restructure_new_year = is_new_year

      !----- A NaN is a STATUS here, not an `error stop`: a library caller survives it. -----------!
      if (is_new_year) then
         if (has_nan(poly%site)) then
            status = DRIVER_ERR_NAN ; poly%status = status ; return
         end if
      end if

      !----- SOIL-CARBON PLAUSIBILITY, checked every step and NOT gated on the ledger. -----------!
      !                                                                                          !
      !      The ledger checks the same predicate at every phase boundary and can say WHICH        !
      !      operator broke it, which is far more useful -- but the ledger is a diagnostic with an  !
      !      off switch ([run].slow_ledger_on), and a safety assertion that a diagnostic flag can  !
      !      disable is a trap. So the guard itself lives here, unconditional, and the ledger's     !
      !      copy is the attribution.                                                              !
      !      Once per slow step over ~12 patches x 7 pools; the cost does not register.            !
      if (cfg%soil_carbon_on) then
         block
            integer(ik) :: ipp, kk
            do ipp = 1_ik, poly%site%patch%n
               kk = soil_carbon_bad_pool(poly%site%patch%soil_carbon(ipp))
               if (kk /= 0_ik) then
                  write(*,'(a)') ' ERROR: a CENTURY soil-carbon pool is physically impossible.'
                  if (len_trim(poly%label) > 0) write(*,'(2a)') '        ', trim(poly%label)
                  write(*,'(3a,i0,a,i0)') '        pool ', trim(soil_carbon_pool_name(kk)),        &
                        ', patch ', ipp, ' of ', poly%site%patch%n
                  write(*,'(2a)')  '        date ', time_to_string(now)
                  write(*,'(a)') '        A carbon pool is a mass: it cannot be negative, and the'
                  write(*,'(a)') '        ceiling is ~500x the richest real soil. Conservation can'
                  write(*,'(a)') '        hold perfectly while this is true -- see the slow ledger.'
                  status = DRIVER_ERR_SOILC ; poly%status = status
                  return
               end if
            end do
         end block
      end if
      poly%status = status

   contains

      !----- step_start is passed UNCONDITIONALLY (leaf phenology needs day-of-year every step);  !
      !      the met source/cursor, the output buffers and the polygon's latitude go with forcing_on. !
      subroutine stepper(fast_bufs)
         type(output_buffers_t), intent(inout) :: fast_bufs
         if (cfg%fast_biophysics_on .and. cfg%forcing%forcing_on) then
            call advance_one_step(poly%site, cfg, poly%fast_ctx,                                   &
                                  met_src=met_src, met_cur=poly%met_cur, step_start=prev,       &
                                  out_bufs=fast_bufs,                                           &
                                  run_energy_budget=poly%energy_budget,                         &
                                  run_water_budget=poly%water_budget,                           &
                                  run_face_budget=poly%face_budget,                             &
                                  slow_ledger=poly%slow_ledger, seam=poly%seam,                 &
                                  latitude_deg=poly%met_cur%latitude_deg)
         else
            call advance_one_step(poly%site, cfg, poly%fast_ctx,                                   &
                                  step_start=prev, run_energy_budget=poly%energy_budget,        &
                                  run_water_budget=poly%water_budget,                           &
                                  run_face_budget=poly%face_budget,                             &
                                  slow_ledger=poly%slow_ledger, seam=poly%seam)
         end if
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
               poly%fast_step_total = poly%fast_step_total + 1_ik
               if (mod(poly%fast_step_total, max(files%fast_interval_steps, 1_ik)) == 0_ik)         &
                  call close_tier(files, bufs, 1_ik)
            end do
            bufs%fast_ready = .false.
         end if
         call output_integrate(files, bufs, poly%site, prev, cfg%dt_slow, is_new_day, is_new_month, &
                               is_new_year)
      end subroutine tick_output

   end subroutine polygon_step

   !---------------------------------------------------------------------------------------!
   ! polygon_report -- the end-of-run conservation reports of one polygon: the two whole-column    !
   ! budgets and the face closure, the soil-carbon seam and freeze, and the slow ledger.           !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_report(cfg, poly)
      type(meds_config_t),  intent(in) :: cfg
      type(meds_polygon_t), intent(in) :: poly
      !----- Whole-column conservation over the RUN: the signed cumulative residual is the number a !
      !      per-step tolerance cannot see, so it is reported whether or not any single step        !
      !      breached.  -----------------------------------------------------------------------------!
      if (cfg%fast_biophysics_on) then
         call budget_report(poly%energy_budget, 'whole_energy', 'J/m2',  'W/m2')
         call budget_report(poly%water_budget,  'whole_water',  'kg/m2', 'kg/m2/s')
         !----- The VERTICAL check (#189), reported in the same place and spirit as the soil-carbon   !
         !      seam below: a number that should be machine-zero, printed whether or not it breached. !
         !      A nonzero value means soil enthalpy was advected on a mass flux the committed water    !
         !      never moved -- a misplacement between LAYERS, which both ledgers above sum away.       !
         write(*,'(a,es12.3,a,es12.3,a,i0)') ' faces[soil_layer_mass]  worst = ',                  &
               poly%face_budget%worst, ' kg/m2   mean |resid| = ',                                  &
               poly%face_budget%abs_sum / max(real(poly%face_budget%n_check, wp), 1.0_wp),          &
               ' kg/m2   checks = ', poly%face_budget%n_check
         if (poly%energy_budget%n_fail + poly%water_budget%n_fail > 0_ik)                        &
            write(*,'(a,i0,a)') ' WARNING: ', poly%energy_budget%n_fail + poly%water_budget%n_fail, &
               ' whole-column budget checks breached tolerance (see [energy].debug_error to make this fatal)'
      end if
      !----- The soil-carbon seam, in the same place and spirit as the two budgets above: a number  !
      !      that should be machine-zero, reported whether or not it ever breached.  ----------------!
      if (cfg%soil_carbon_on) then
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_rh]  worst |pool debit - fast Rh| = ',        &
               poly%seam%worst_rh_gap, ' kgC/m2'
         if (poly%seam%worst_rh_gap > 0.0_wp)                                                     &
            write(*,'(3a,i0,a)') '                       worst at ', poly%worst_rh_seam_when,     &
                  ' with ', poly%worst_rh_seam_npatch, ' patches'
         write(*,'(a,es12.3,a)') ' seam[soil_carbon_lignin]  worst passive-tracer resid = ',      &
               poly%seam%worst_lignin, ' kgC/m2/day'
         !----- LAMBDA: the fraction of a pool one slow step withdraws -- the number that says      !
         !      whether freezing the pool across the step is sound. Not a residual: it is expected   !
         !      to be small-but-nonzero (3.6e-3/day measured on the mature Ithaca stand), and it is a  !
         !      WARNING as it approaches 1, long before a pool could go negative.  ------------------!
         write(*,'(a,es12.3,2a)') ' freeze[soil_carbon]   worst lambda = dt.K.xi / pool  = ',     &
               poly%seam%worst_lambda, ' [-]  pool ',                                             &
               trim(soil_carbon_pool_name(poly%seam%lambda_pool))
         if (poly%seam%worst_lambda > 0.1_wp)                                                     &
            write(*,'(a)') '   WARNING: a slow step withdraws >10% of a soil-carbon pool --'      &
                        // ' the frozen-pool approximation is degrading; shorten dt_slow.'
      end if
      !----- The SLOW tier's ledger, over the window the two above cannot see (plan §10.2). -------!
      call slow_ledger_report(poly%slow_ledger)
   end subroutine polygon_report

   !---------------------------------------------------------------------------------------!
   ! Advance the polygon's growth-temperature running mean and push it into the leaf table (#176). !
   ! Exponential mean with weight dt/window; UNSEEDED (negative) adopts the first day's mean       !
   ! outright, so the acclimated optimum starts from the site's own climate rather than relaxing   !
   ! toward it from an arbitrary origin over the first month.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine advance_growth_temperature(cfg, step_days, poly)
      type(meds_config_t),  intent(in)    :: cfg
      integer(ik),          intent(in)    :: step_days
      type(meds_polygon_t), intent(inout) :: poly
      real(wp) :: t_day, w
      if (poly%site%pheno_tair_n < 1_ik) return          ! no fast sub-steps ran: nothing to average
      t_day = poly%site%pheno_tair_sum / real(poly%site%pheno_tair_n, wp)
      if (poly%site%t_growth_avg <= 0.0_wp) then
         poly%site%t_growth_avg = t_day
      else
         w = min(1.0_wp, real(step_days, wp) / max(cfg%acclim_window_days, real(step_days, wp)))
         poly%site%t_growth_avg = poly%site%t_growth_avg + w * (t_day - poly%site%t_growth_avg)
      end if
      call acclimate_leaf_photo_table(cfg, poly%site%t_growth_avg, poly%fast_ctx%col_config%leaf_photo)
   end subroutine advance_growth_temperature

   !----- Seed an initial snow pack (spin-up / test); [fast].snow_init_swe. ---------------------!
   subroutine seed_snow(cfg, poly, verbose)
      type(meds_config_t),  intent(in)    :: cfg
      type(meds_polygon_t), intent(inout) :: poly
      logical,              intent(in)    :: verbose
      integer(ik) :: ipp
      do ipp = 1_ik, poly%site%patch%n
         poly%site%patch%snow(ipp)%swe(1)         = cfg%snow_init_swe
         poly%site%patch%snow(ipp)%snow_energy(1) =                                              &
            temp_to_internal_energy(0.0_wp, cfg%snow_init_swe, cfg%snow_init_temp, 0.0_wp)
         poly%site%patch%snow(ipp)%snow_depth(1)  = cfg%snow_init_swe / 250.0_wp
         poly%site%patch%snow(ipp)%nlayer         = 1_ik
         call internal_energy_to_temp(poly%site%patch%snow(ipp)%snow_energy(1),                  &
                                      cfg%snow_init_swe, 0.0_wp,                                 &
                                      poly%site%patch%snow(ipp)%snow_temp(1),                    &
                                      poly%site%patch%snow(ipp)%snow_fliq(1))
      end do
      if (verbose) write(*,'(a,f6.1,a)') ' snow  : seeded initial pack SWE = ', cfg%snow_init_swe, ' kg/m2'
   end subroutine seed_snow

   !----- SASU steady-state soil-carbon spin-up from the configured climatological xi + litter. --!
   subroutine soil_carbon_steady(cfg, poly, verbose)
      type(meds_config_t),  intent(in)    :: cfg
      type(meds_polygon_t), intent(inout) :: poly
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
      do ipp = 1_ik, poly%site%patch%n
         poly%site%patch%soil_carbon(ipp) = pools_ss
      end do
      if (verbose) write(*,'(a,f8.3,a)') ' soilc : steady-state spin-up, total = ',               &
         pools_ss%fast_grnd_carbon + pools_ss%fast_soil_carbon + pools_ss%struct_grnd_carbon +   &
         pools_ss%struct_soil_carbon + pools_ss%microbial_carbon + pools_ss%slow_carbon +        &
         pools_ss%passive_carbon, ' kgC/m2'
   end subroutine soil_carbon_steady

end module meds_polygon
