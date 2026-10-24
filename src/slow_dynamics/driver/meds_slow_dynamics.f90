! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_slow_dynamics -- the order of the slow loop (docs/science/order_of_processes.md     !
! §4-§5), and where the slow ledger marks its phases (conservation.md §3). It loops over   !
! nothing: meds_slow_state reads the stand and calls the processes, meds_slow_update       !
! writes the stand, and meds_slow_diag writes the slow rows.                               !
!   * advance_slow_dynamics     -- one slow step: vegetation_dynamics, patch ageing, the   !
!                                  canopy-air depth, then the soil-carbon step.            !
!   * vegetation_dynamics       -- phenology, plasticity, allocation, the vital rates, and !
!                                  the commit of growth and mortality.                     !
!   * advance_boundary_dynamics -- a calendar boundary: restructure_stand (the monthly     !
!                                  cohort operators, the yearly disturbance and patch      !
!                                  operators), then the canopy-air depth.                  !
!   * settle_tissue_water       -- before the fast loop, the tissue water above capacity.  !
!==========================================================================================!
module meds_slow_dynamics
   use meds_kinds,              only : wp, ik
   use meds_config,             only : meds_config_t
   use meds_time,               only : meds_time_t, day_of_year
   use meds_polygon_state_vars, only : polygon_t
   use meds_slow_ledger,        only : slow_ledger_t, soilc_seam_t, slow_ledger_open, slow_ledger_mark,  &
                                       slow_ledger_rebase, slow_tissue_energy, slow_ledger_declare,      &
                                       SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC, SLOW_PHASE_TISSUE_WATER,     &
                                       slow_fast_carbon_handover, SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW,  &
                                       SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT, SLOW_PHASE_DISTURB,        &
                                       SLOW_PHASE_PATCH, slow_co2_handoff, litter_carbon_total
   use meds_slow_state,         only : slow_step_data_t, advance_leaf_phenology, advance_plant_traits,   &
                                       compute_carbon_allocation, shed_turnover_water,                   &
                                       compute_vital_rates, fill_new_leaf_water,                         &
                                       accumulate_mortality_litter, shed_mortality_water,                &
                                       shed_excess_tissue_water, advance_soil_carbon
   use meds_slow_update,        only : accumulate_recruit_pool, commit_vegetation_step, sort_cohorts,    &
                                       update_overtopping_lai, age_patches, refresh_canopy_depth,        &
                                       apply_recruitment, fuse_cohorts, terminate_cohorts,               &
                                       split_cohorts, apply_patch_disturbance, sort_patches,             &
                                       fuse_patches, terminate_patches
   use meds_slow_diag,          only : open_slow_cohort_rows, write_slow_patch_rows
   implicit none
   private

   public :: advance_slow_dynamics, advance_boundary_dynamics, settle_tissue_water
   public :: vegetation_dynamics, restructure_stand

   !----- Patch disturbance runs once a year, so it integrates its yearly rate over this interval. !
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !< [yr]

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step (order_of_processes.md §4): vegetation dynamics, patch ageing,  !
   ! the canopy-air depth, then the daily soil-carbon step. Without `step_start` (a bare   !
   ! test call) phenology is not advanced.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine advance_slow_dynamics(poly, cfg, step_start, seam, ledger, rho_air, latitude_deg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(meds_time_t),   intent(in), optional :: step_start   !< the slow step's start (its day of year)
      type(soilc_seam_t),  intent(inout), optional :: seam   !< run-worst fast/slow seam checks (soil_carbon.md §7)
      !----- The slow ledger (optional: a test may step without one) and its air density;  !
      !      the stepper passes the constant rho_air_ref (conservation.md §2). ------------!
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp),            intent(in),    optional :: rho_air
      real(wp),            intent(in),    optional :: latitude_deg   !< the polygon's (default [site])

      if (present(ledger)) call slow_ledger_open(ledger, poly, cfg, rho_air)

      if (present(step_start)) then
         call vegetation_dynamics(poly, cfg, day_of_year(step_start), ledger, latitude_deg)
      else
         call vegetation_dynamics(poly, cfg, ledger=ledger, latitude_deg=latitude_deg)
      end if

      call age_patches(poly%patch, cfg%dt_years)

      !----- Resize each canopy-air volume to the stand growth and mortality left, so the  !
      !      fast loop sees one depth all day. --------------------------------------------!
      call refresh_canopy_depth(poly, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)

      !----- `seam` accumulates over the run; with soil carbon off it stays as it was. ----!
      if (cfg%soil_carbon%soil_carbon_on) call advance_soil_carbon(poly, cfg, seam, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_SOILC)
   end subroutine advance_slow_dynamics

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation one slow step: phenology, trait plasticity, carbon allocation, the  !
   ! vital rates, and the growth and mortality commit. The calendar's restructuring runs at the !
   ! boundary, in restructure_stand (docs/science/order_of_processes.md §4).                   !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(poly, cfg, doy, ledger, latitude_deg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in), optional :: doy   !< day-of-year at the step start (drives phenology)
      type(slow_ledger_t), intent(inout), optional :: ledger    !< site conservation ledger (meds_slow_ledger)
      real(wp),            intent(in),    optional :: latitude_deg !< the polygon's latitude [deg] (default [site])
      type(slow_step_data_t) :: step
      real(wp)               :: mort_water, seed_rain_carbon
      real(wp)               :: tissue_heat0, tissue_heat1, handover, plastic_heat0

      !----- 0. Leaf phenology: advance the leaf tendencies one daily step, before the allocation     !
      !         reads them (docs/science/plant_phenology.md). A caller with no calendar (the C-API's   !
      !         carbon mode) omits doy, and the tendencies are not advanced. -----------------------!
      if (present(doy)) call advance_leaf_phenology(poly, cfg, doy, latitude_deg)

      !----- 0b. Light trait plasticity (opt-in): acclimate the leaf traits toward their shaded      !
      !          targets from last step's overtopping LAI, before the allocation reads them. The     !
      !          leaf carbon an SLA rise resorbs to storage (no heat capacity) changes the tissue     !
      !          energy at unchanged temperature: declared (docs/science/conservation.md §4). -------!
      if (cfg%trait_plasticity_on) then
         if (present(ledger)) plastic_heat0 = slow_tissue_energy(poly)
         call advance_plant_traits(poly, cfg, cfg%dt_years)
         if (present(ledger)) call slow_ledger_declare(ledger, energy_in = slow_tissue_energy(poly) - plastic_heat0)
      end if

      !----- 1. Carbon allocation, which adds the turnover litter to the per-patch accumulator. ----!
      call open_slow_cohort_rows(poly)
      call compute_carbon_allocation(poly, cfg, step)

      !----- 1b. Shed the turnover's tissue water before the commit, which changes the pools the    !
      !          loss is measured against (docs/science/soil_plant_water.md §2). ------------------!
      call shed_turnover_water(poly, cfg, step, ledger)
      !----- Every movement of the allocation phase is declared, so a residual here is a movement   !
      !      nobody declared (docs/science/conservation.md §3). -------------------------------------!
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_ALLOCATE)

      !----- 2. The vital rates, before the commit: mortality reads the growth average of the steps !
      !         before this one (docs/science/demography.md §1). ------------------------------------!
      call compute_vital_rates(poly, cfg, step)

      !----- 2b. Credit the recruit pool with this step's recruitment, and declare the seed rain    !
      !          as carbon in: no parent at the site paid for it (docs/science/demography.md §2). --!
      call accumulate_recruit_pool(poly, cfg, step%recruitment, cfg%dt_years, seed_rain_carbon)
      if (present(ledger)) call slow_ledger_declare(ledger, carbon_in = seed_rain_carbon)

      !----- 3. Growth and mortality. The commit and the new leaves' water leave the temperatures  !
      !         alone, so the energy change is the tissues' thermal mass: one exchange, declared   !
      !         below, not a leak (docs/science/conservation.md §4). The fast loop's net carbon is  !
      !         valued at the pre-commit density it accumulated on; after the commit it would miss  !
      !         the dying plants' share (conservation.md §3). -----------------------------------!
      tissue_heat0 = slow_tissue_energy(poly)
      handover     = slow_fast_carbon_handover(poly, cfg)
      call commit_vegetation_step(poly, cfg, step)
      call fill_new_leaf_water(poly, cfg, step)
      tissue_heat1 = slow_tissue_energy(poly)

      !----- Mortality's litter and water, on the density the commit removed and the pools it left  !
      !      (docs/science/demography.md eq. (9)); recorded whether or not soil carbon is on. -----!
      call accumulate_mortality_litter(poly, cfg, step)
      call shed_mortality_water(poly, cfg, step, mort_water)

      !----- Growth changed heights: re-sort tallest first (docs/science/demography.md §3). -------!
      call sort_cohorts(poly)

      !----- The growth phase declares the fast loop's net carbon (in), the litter and the CO2 owed  !
      !      to the canopy air (out), the dying plants' water and the thermal mass, so it closes to  !
      !      round-off (docs/science/conservation.md §3). ------------------------------------------!
      if (present(ledger)) then
         call slow_ledger_declare(ledger, carbon_in  = handover,                                    &
                                          carbon_out = litter_carbon_total(poly, poly%patch%litter_in)  &
                                                     + slow_co2_handoff(poly, cfg),                &
                                          water_out  = mort_water,                                 &
                                          energy_in  = tissue_heat1 - tissue_heat0)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_GROW)
      end if

      !----- The slow patch rows, before the restructuring permutes the patch axis. -------------!
      call write_slow_patch_rows(poly, cfg, step%recruitment)

      !----- 4. Refresh the overtopping LAI on the sorted stand (docs/science/demography.md §6). --!
      call update_overtopping_lai(poly)
   end subroutine vegetation_dynamics

   !---------------------------------------------------------------------------------------!
   ! advance_boundary_dynamics -- the calendar boundary (order_of_processes.md §5):        !
   ! restructure the stand, then resize the canopy-air volumes; a no-op off a month        !
   ! boundary. Its marks continue the last step's ledger window.                           !
   !---------------------------------------------------------------------------------------!
   subroutine advance_boundary_dynamics(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      type(slow_ledger_t), intent(inout), optional :: ledger
      if (.not. is_new_month) return
      if (present(ledger)) then
         if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
      end if
      call restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      call refresh_canopy_depth(poly, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)
   end subroutine advance_boundary_dynamics

   !---------------------------------------------------------------------------------------!
   ! restructure_stand -- the calendar's restructuring: monthly recruitment, cohort fusion,   !
   ! cull, fission and sort; yearly patch disturbance, then patch fusion and termination. It  !
   ! runs at the boundary, after the output has read the period's last step, so each record   !
   ! describes one stand and the boundary's events belong to the period it opens              !
   ! (docs/science/order_of_processes.md §5).                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year   !< which boundary this is
      type(slow_ledger_t), intent(inout), optional :: ledger    !< site conservation ledger (meds_slow_ledger)
      real(wp) :: cull_water, th0, th1, rec_carbon, rec_heat, rec_water, dist_water
      logical  :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      !----- The calendar cadence, and the demography and fusion/fission switches. ----------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- Cohort restructuring (monthly): recruitment, fusion, cull, fission and sort. ----!
      if (do_cohort_fissfuse) then
         call apply_recruitment(poly, cfg, rec_carbon, rec_heat, rec_water)
         if (present(ledger)) call slow_ledger_declare(ledger, carbon_in = rec_carbon,             &
                                                               energy_in = rec_heat,               &
                                                               water_in  = rec_water)
         !----- The operators below change the tissues' thermal mass (fusion and fission re-derive  !
         !      the sapwood ring, nonlinear in dbh; the cull removes a cohort): one exchange,        !
         !      declared (docs/science/conservation.md §4). ----------------------------------------!
         if (present(ledger)) th0 = slow_tissue_energy(poly)
         !----- Recruitment is marked on its own: it alone draws carbon, heat and water from outside !
         !      the model, checkable only if no other operator shares its phase                     !
         !      (docs/science/conservation.md §3). -------------------------------------------------!
         if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_RECRUIT)
         call fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg, cull_water)
         call split_cohorts(poly, cfg)
         call sort_cohorts(poly)
         if (present(ledger)) then
            th1 = slow_tissue_energy(poly)
            call slow_ledger_declare(ledger, water_out = cull_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_COHORT)
         end if
      end if

      !----- Patch disturbance, then patch fusion and termination (yearly, each switchable). -!
      if (do_patch_disturbance) then
         if (present(ledger)) th0 = slow_tissue_energy(poly)
         call apply_patch_disturbance(poly, cfg, PATCH_DYNAMICS_INTERVAL, dist_water)
         if (present(ledger)) then
            th1 = slow_tissue_energy(poly)
            call slow_ledger_declare(ledger, water_out = dist_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_DISTURB)
         end if
      end if
      if (do_patch_fissfuse) then
         if (present(ledger)) th0 = slow_tissue_energy(poly)
         call sort_patches(poly)
         call fuse_patches(poly, cfg)
         call terminate_patches(poly, cfg)
         call fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg, cull_water)
         call sort_cohorts(poly)
         if (present(ledger)) then
            th1 = slow_tissue_energy(poly)
            call slow_ledger_declare(ledger, water_out = cull_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_PATCH)
         end if
      end if

      !----- Refresh the overtopping LAI for the restructured stand. ---------------------------!
      call update_overtopping_lai(poly)
   end subroutine restructure_stand

   !---------------------------------------------------------------------------------------!
   ! settle_tissue_water -- before the fast loop, shed the tissue water above what the     !
   ! biomass holds and declare it: the water as a handoff to the ground, the heat as a     !
   ! thermal-mass change (conservation.md §3-§4), in a phase of its own.                   !
   !---------------------------------------------------------------------------------------!
   subroutine settle_tissue_water(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp) :: water_shed, heat0
      heat0 = 0.0_wp
      if (present(ledger)) then
         if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
         heat0 = slow_tissue_energy(poly)
      end if
      call shed_excess_tissue_water(poly, cfg, water_shed)
      if (present(ledger)) then
         call slow_ledger_declare(ledger, water_out = water_shed, energy_in = slow_tissue_energy(poly) - heat0)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_TISSUE_WATER)
      end if
   end subroutine settle_tissue_water

end module meds_slow_dynamics
