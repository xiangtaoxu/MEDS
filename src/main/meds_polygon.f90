! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_polygon -- one polygon of a run and the step that advances it                       !
! (docs/science/order_of_processes.md §2).                                                 !
!                                                                                          !
! A polygon is two records side by side in the run: its state (polygon_t) and what the run !
! keeps to step and check it (polygon_run_t: its forcing cursor, its conservation ledgers, !
! its status). The second is not part of the state because its types come from the fast    !
! and slow libraries, which use the state library. A run is one polygon or many, each      !
! stepped by polygon_step, so a polygon of a region computes what a one-polygon run at its !
! cell computes.                                                                           !
!                                                                                          !
! polygon_step does no file work: the caller loads the step's forcing before it            !
! (met_prefetch) and writes the queued output records after a month closes.                !
!==========================================================================================!
module meds_polygon
   use, intrinsic :: ieee_arithmetic,   only : ieee_is_finite
   use meds_kinds,                      only : wp, ik
   use meds_config,                     only : meds_config_t, INIT_RESTART
   use meds_time,                       only : meds_time_t, calendar_step_t, time_to_string
   use meds_budget_check,               only : budget_t, budget_report
   use meds_cohort_state_vars,          only : cohort_bad_field, cohort_field_finite, cohort_field_name
   use meds_patch_state_vars,           only : canopy_air_bad_field, canopy_air_field_name, soil_carbon_state_t, &
                                               litter_input_t
   use meds_polygon_state_vars,         only : polygon_t, reset_step_diagnostics
   use meds_stepper,                    only : advance_one_step, advance_boundary
   use meds_slow_update,                only : refresh_canopy_depth, max_cohort_count
   use meds_slow_state,                 only : acclimate_site_leaf_tables
   use meds_slow_ledger,                only : slow_ledger_t, slow_ledger_report
   use meds_biogeochem_params,          only : N_SOIL_POOL
   use meds_soil_carbon,                only : solve_soil_carbon_steady_state, build_litter_input,      &
                                               soil_carbon_bad_pool, soil_carbon_pool_name
   use meds_respiration,                only : decomposition_transfer_matrix
   use meds_met_types,                  only : met_source_t, met_cursor_t
   use meds_met,                        only : met_cursor_init
   use meds_output_types,               only : output_files_t, output_buffers_t
   use meds_output_integrate,           only : output_tick
   use meds_output_reduce,              only : total_area
   implicit none
   private

   public :: polygon_run_t, polygon_prepare, polygon_step, polygon_report
   public :: RUN_OK, RUN_FINISHED, RUN_ERR_NAN, RUN_ERR_AREA, RUN_ERR_BOUNDS

   !----- A step's or a run's status. A failure is a status, not an `error stop`, so a library     !
   !      caller (Python) survives it. ------------------------------------------------------------!
   integer(ik), parameter :: RUN_OK         = 0_ik   !< a slow step was taken
   integer(ik), parameter :: RUN_FINISHED   = 1_ik   !< the calendar had reached end_time; nothing done
   integer(ik), parameter :: RUN_ERR_NAN    = 2_ik   !< NaN in the state or the energy/water ledger
   integer(ik), parameter :: RUN_ERR_AREA   = 3_ik   !< patch areas do not sum to 1 (finalize only)
   integer(ik), parameter :: RUN_ERR_BOUNDS = 4_ik   !< a cohort or soil-carbon state outside its bounds

   !----- What the run keeps to step and check one polygon; its state is the polygon_t beside      !
   !      it. -------------------------------------------------------------------------------------!
   type :: polygon_run_t
      character(len=40)   :: label = ''                   !< how messages name it ('' for a site run)
      integer(ik)         :: patch_threads = 1_ik         !< the patch loop's threads ([run].n_threads, or 1)
      logical             :: seed_canopy_air = .false.    !< a start's first window seeds the canopy air
      type(met_cursor_t)  :: met_cur                      !< set only if forcing_on
      type(budget_t)      :: energy_budget, water_budget  !< the whole-column ledgers over the run
      type(budget_t)      :: face_budget                  !< each layer face's water: did the heat move with it
      type(slow_ledger_t) :: slow_ledger                  !< the site's stores across each slow step, the seam
      !----- The most patches, and the most cohorts in one patch, a restructuring left, and when:  !
      !      above max_patch or max_cohort, the fusion tolerance ceilings kept dissimilar ones     !
      !      apart. -------------------------------------------------------------------------------!
      integer(ik)         :: most_patches = 0_ik, most_cohorts = 0_ik
      character(len=19)   :: most_patches_when = '', most_cohorts_when = ''
      real(wp)            :: area_start = 0.0_wp          !< the patch areas' sum at the start
      integer(ik)         :: status = RUN_OK           !< the last step's status
   end type polygon_run_t

contains

   !---------------------------------------------------------------------------------------!
   ! polygon_prepare -- what follows the initial state (init_polygon), on a freshly        !
   ! allocated run record: the forcing cursor at the polygon's cell and location, the      !
   ! canopy-air depth for the stand, and the soil-carbon spin-up, unless a restart         !
   ! restored the pools (keep_soil_carbon).                                                !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_prepare(cfg, met_src, cell, poly, run, keep_soil_carbon, patch_threads, verbose)
      type(meds_config_t), intent(in)    :: cfg
      type(met_source_t),  intent(in)    :: met_src
      integer(ik),         intent(in)    :: cell            !< the polygon's cell in the forcing source's list
      type(polygon_t),     intent(inout) :: poly
      type(polygon_run_t), intent(inout) :: run
      logical,             intent(in)    :: keep_soil_carbon, verbose
      integer(ik),         intent(in)    :: patch_threads   !< [run].n_threads for one polygon, 1 for a region's
      type(slow_ledger_t) :: no_ledger                      !< inactive: the first canopy-air depth books nothing
      run%area_start      = total_area(poly)
      run%patch_threads   = patch_threads
      run%seed_canopy_air = cfg%init%init_mode /= INIT_RESTART
      if (cfg%fast_biophysics_on) then
         if (cfg%forcing%forcing_on) then
            call met_cursor_init(met_src, run%met_cur, cell, poly%latitude_deg, poly%longitude_deg,      &
                                 poly%elevation_m)
         end if
         if (verbose) write(*,'(a)') ' fast  : sub-daily biophysics ON'
      end if
      !----- The canopy-air depth follows the stand from the first step, whatever the stand came   !
      !      from; the slow step keeps it current after that. A plain geometry update that keeps   !
      !      the canopy air's intensive state, so no ledger books it. -----------------------------!
      call refresh_canopy_depth(poly, cfg, no_ledger)
      if (cfg%soil_carbon%soil_carbon_on .and. cfg%soil_carbon%spinup_steady .and. .not. keep_soil_carbon)    &
         call soil_carbon_steady(poly, cfg, verbose)
      run%slow_ledger%active = cfg%slow_on .and. cfg%slow_ledger_on
   end subroutine polygon_prepare

   !---------------------------------------------------------------------------------------!
   ! polygon_step -- one slow step of one polygon, `step`: the restructuring the last      !
   ! step's boundary owes, the thermal acclimation, the coupled step, the output tick and  !
   ! the state guard. The step's forcing is already loaded into `met_src`; closed output   !
   ! records are queued in the polygon's buffers, `out_bufs`.                              !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_step(cfg, met_src, out_files, out_bufs, poly, run, step, status)
      type(meds_config_t),    intent(in)    :: cfg
      type(met_source_t),     intent(in)    :: met_src
      type(output_files_t),   intent(in)    :: out_files
      type(output_buffers_t), intent(inout) :: out_bufs   !< this polygon's share of out_files
      type(polygon_t),        intent(inout) :: poly
      type(polygon_run_t),    intent(inout) :: run
      type(calendar_step_t),  intent(in)    :: step
      integer(ik),            intent(out)   :: status
      integer(ik) :: n_cohort_max

      !----- The boundary the last step ended on: the stand's monthly and yearly restructuring,    !
      !      before this step's fast loop sees the stand, and the most patches and cohorts it      !
      !      leaves. ------------------------------------------------------------------------------!
      if (poly%restructure_pending) then
         call advance_boundary(poly, cfg, .true., poly%restructure_new_year,                       &
                               slow_ledger=run%slow_ledger)
         poly%restructure_pending  = .false.
         poly%restructure_new_year = .false.
         if (poly%patch%n > run%most_patches) then
            run%most_patches      = poly%patch%n
            run%most_patches_when = time_to_string(step%start)
         end if
         n_cohort_max = max_cohort_count(poly)
         if (n_cohort_max > run%most_cohorts) then
            run%most_cohorts      = n_cohort_max
            run%most_cohorts_when = time_to_string(step%start)
         end if
      end if

      !----- Thermal acclimation (docs/science/leaf_gas_exchange.md §2): each site's               !
      !      growth-temperature mean takes the daily mean the last step's fast loop accumulated,   !
      !      and its leaf table follows. The mean tracks the preceding weeks, so one step's delay  !
      !      is a thirtieth of its window. --------------------------------------------------------!
      if (cfg%leaf_physiology%thermal_acclimation) call acclimate_site_leaf_tables(poly, cfg%leaf_physiology, step%days)

      !----- The coupled step: the fast loop, which stages its sub-steps in `out_bufs` when the    !
      !      output reads them, then the slow dynamics. -------------------------------------------!
      call advance_one_step(poly, cfg, met_src=met_src, met_cur=run%met_cur, step_start=step%start,      &
                            out_bufs=out_bufs, run_energy_budget=run%energy_budget,                    &
                            run_water_budget=run%water_budget, run_face_budget=run%face_budget,         &
                            slow_ledger=run%slow_ledger, patch_threads=run%patch_threads,              &
                            seed_canopy_air=run%seed_canopy_air)

      !----- Output: replay the staged FAST samples, then fold this step -- its fluxes over the    !
      !      step and its state at its end -- into the slower tiers' windows that hold its start,  !
      !      closing each period its end has left and queueing it for the I/O phase. --------------!
      if (out_files%enabled) call output_tick(out_files, out_bufs, poly, step, cfg%dt_slow, cfg%dt_fast)

      !----- The output has read the step: zero its diagnostics. A step that ends on a month       !
      !      boundary leaves the restructuring owed to the next step, so the ending period's       !
      !      records describe one stand and what the restructuring does belongs to the period that !
      !      begins (docs/science/diagnostics.md §4). ---------------------------------------------!
      call reset_step_diagnostics(poly)
      poly%restructure_pending  = step%new_month
      poly%restructure_new_year = step%new_year

      status = check_polygon_state(cfg, poly, run, step%end)
      run%status = status

   end subroutine polygon_step

   !---------------------------------------------------------------------------------------!
   ! check_polygon_state -- the guard every step ends with: the run's energy and water     !
   ! ledgers, each cohort's state, each patch's canopy air and, with soil carbon on, each  !
   ! soil-carbon pool (docs/science/order_of_processes.md §2). The first value found that  !
   ! is not finite is RUN_ERR_NAN, the first outside its bounds RUN_ERR_BOUNDS, reported   !
   ! with what and where. An operator keeps these in bounds by construction, so either is  !
   ! a bug, and stopping on the step it happens is what makes it findable.                 !
   !---------------------------------------------------------------------------------------!
   integer(ik) function check_polygon_state(cfg, poly, run, now) result(status)
      type(meds_config_t), intent(in) :: cfg
      type(polygon_t),     intent(in) :: poly
      type(polygon_run_t), intent(in) :: run
      type(meds_time_t),   intent(in) :: now
      character(len=160) :: where
      integer(ik)        :: i, k
      status = RUN_OK
      !----- One NaN flux or store turns a ledger's running residual NaN for good. ----------------!
      if (cfg%fast_biophysics_on) then
         if (.not. (ieee_is_finite(run%energy_budget%resid_sum) .and.                              &
                    ieee_is_finite(run%water_budget%resid_sum))) then
            call report_failure(run, 'a NaN reached the energy/water ledger.', '',                   &
                                'in the step ending '//time_to_string(now))
            status = RUN_ERR_NAN
            return
         end if
      end if
      do i = 1_ik, poly%cohort%n
         k = cohort_bad_field(poly%cohort, i)
         if (k == 0_ik) cycle
         write(where,'(3a,i0,a,i0,a,i0)') 'field ', trim(cohort_field_name(k)), ', cohort id ',     &
               poly%cohort%global_id(i), ' (pft ', poly%cohort%pft(i), ') in patch ', poly%cohort%owner_patch(i)
         if (cohort_field_finite(poly%cohort, i, k)) then
            call report_failure(run, 'a cohort''s state is outside its physical bounds.', where,       &
                                'date '//time_to_string(now),                                          &
                                [character(len=64) :: 'Density, carbon and water are amounts: never negative, never', &
                                                      'infinite. Every operator keeps them so by construction.'])
            status = RUN_ERR_BOUNDS
         else
            call report_failure(run, 'a cohort''s state is not finite.', where, 'date '//time_to_string(now))
            status = RUN_ERR_NAN
         end if
         return
      end do
      if (allocated(poly%patch%cas)) then
         do i = 1_ik, poly%patch%n
            k = canopy_air_bad_field(poly%patch%cas(i))
            if (k == 0_ik) cycle
            write(where,'(3a,i0,a,i0)') 'field ', trim(canopy_air_field_name(k)), ', patch ', i, ' of ',  &
                  poly%patch%n
            call report_failure(run, 'a patch''s canopy air is not finite.', where, 'date '//time_to_string(now))
            status = RUN_ERR_NAN
            return
         end do
      end if
      if (cfg%soil_carbon%soil_carbon_on) then
         do i = 1_ik, poly%patch%n
            k = soil_carbon_bad_pool(poly%patch%soil_carbon(i))
            if (k == 0_ik) cycle
            write(where,'(3a,i0,a,i0)') 'pool ', trim(soil_carbon_pool_name(k)), ', patch ', i, ' of ',   &
                  poly%patch%n
            call report_failure(run, 'a CENTURY soil-carbon pool is physically impossible.', where,       &
                                'date '//time_to_string(now),                                          &
                                [character(len=64) :: 'A carbon pool is a mass: it cannot be negative, and the',      &
                                                      'ceiling is ~500x the richest real soil. Conservation can',     &
                                                      'hold perfectly while this is true -- see the slow ledger.'])
            status = RUN_ERR_BOUNDS
            return
         end do
      end if
   end function check_polygon_state

   !----- The guard's report: what failed, which polygon (a region's), where in it, when, and a    !
   !      note. -----------------------------------------------------------------------------------!
   subroutine report_failure(run, what, where, when, note)
      type(polygon_run_t), intent(in) :: run
      character(len=*),    intent(in) :: what, where, when
      character(len=*),    intent(in), optional :: note(:)
      integer :: i
      write(*,'(2a)') ' ERROR: ', what
      if (len_trim(run%label) > 0) write(*,'(2a)') '        ', trim(run%label)
      if (len_trim(where) > 0) write(*,'(2a)') '        ', trim(where)
      write(*,'(2a)') '        ', when
      if (present(note)) then
         do i = 1, size(note)
            write(*,'(2a)') '        ', trim(note(i))
         end do
      end if
   end subroutine report_failure

   !---------------------------------------------------------------------------------------!
   ! polygon_report -- one polygon's end-of-run conservation reports: the two whole-column !
   ! budgets and the face closure, the fusion ceilings, and the slow ledger with the       !
   ! soil-carbon seam. Each is printed whether or not it ever breached: the run's signed   !
   ! cumulative residual is the number a per-step tolerance cannot see.                    !
   !---------------------------------------------------------------------------------------!
   subroutine polygon_report(cfg, run)
      type(meds_config_t), intent(in) :: cfg
      type(polygon_run_t), intent(in) :: run
      if (cfg%fast_biophysics_on) then
         call budget_report(run%energy_budget, 'whole_energy', 'J/m2',  'W/m2')
         call budget_report(run%water_budget,  'whole_water',  'kg/m2', 'kg/m2/s')
         !----- The face closure should be machine zero: a nonzero value means soil enthalpy moved !
         !      on a mass flux the committed water never moved, a misplacement between layers that !
         !      both whole-column budgets sum away. -----------------------------------------------!
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
      if (cfg%do_patch_fissfuse .and. cfg%max_patch > 0_ik .and. run%most_patches > cfg%max_patch)     &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: patch fusion left ', run%most_patches, ' patches on ',         &
               trim(run%most_patches_when), ' (max_patch = ', cfg%max_patch,                            &
               '): patch_light_tol_max kept dissimilar patches apart'
      if (cfg%do_cohort_fissfuse .and. cfg%max_cohort > 0_ik .and. run%most_cohorts > cfg%max_cohort)  &
         write(*,'(a,i0,3a,i0,a)') ' NOTE: cohort fusion left ', run%most_cohorts, ' cohorts in one patch on ',  &
               trim(run%most_cohorts_when), ' (max_cohort = ', cfg%max_cohort,                          &
               '): cohort_size_tol_max kept dissimilar cohorts apart'
      call slow_ledger_report(run%slow_ledger)
   end subroutine polygon_report

   !----- The SASU steady-state soil-carbon spin-up from the configured climatological xi and      !
   !      litter (docs/science/soil_carbon.md §8, eq. (7)), into every patch. ---------------------!
   subroutine soil_carbon_steady(poly, cfg, verbose)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: verbose
      type(soil_carbon_state_t) :: pools0, pools_ss
      type(litter_input_t) :: lit
      real(wp)    :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp)    :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL), lignin_bar(2)
      integer(ik) :: ipp
      pools0 = soil_carbon_state_t()                    ! zero-carbon reference (f_lignin = 0)
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

end module meds_polygon
