! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_vegetation_dynamics -- the slow loop's vegetation driver (library meds_slow). Each slow  !
! step (vegetation_dynamics) it advances leaf phenology and trait plasticity, allocates the     !
! day's carbon, takes the vital rates from the laws of meds_demography_rates and has the        !
! demography engine apply them; at a calendar boundary restructure_stand runs the restructuring. !
! It is the one place plant and demography meet: the engine only applies the arrays it is given !
! (src/README.md rule 4). docs/science/order_of_processes.md §4-§5, demography.md §1-§2.         !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : day_sec, tiny_num, yr_day, grav_head
   use meds_config,               only : meds_config_t, growth_window_steps
   use meds_pft_derive,           only : pft_leaf_psi_tlp
   use meds_plant_hydr_lib,       only : water_curve_t, tissue_psi_to_water, tissue_water_to_psi, tissue_capped_water
   use meds_therm_lib,            only : internal_energy_liquid
   use meds_allometry_lib,        only : size_to_leaf_carbon, carbon_to_structure, min_cohort_carbon,  &
                                         storage_carbon_target
   use meds_cohort_state_vars,    only : cohort_deriv_alloc, set_cohort_size_from_carbon, cohort_wood_store_carbon, &
                                         hold_cohort_tissue_temps, cohort_water_enthalpy
   use meds_polygon_state_vars, only : polygon_t, patch_polygon_area
   use meds_slow_ledger,          only : slow_ledger_t, slow_ledger_mark, slow_ledger_declare,     &
                                         slow_fast_carbon_handover, SLOW_PHASE_ALLOCATE,           &
                                         SLOW_PHASE_GROW, SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT,   &
                                         SLOW_PHASE_DISTURB, SLOW_PHASE_PATCH, KGC_PER_UMOL_C,   &
                                         slow_tissue_energy
   use meds_demography_update, only : update_cohort_states, fill_cohort_deriv, advance_pool,           &
                                      update_overtopping_lai
   use meds_demography_cohort_fusefiss, only : apply_recruitment, new_fuse_cohorts, terminate_cohorts, split_cohorts, sort_cohorts
   use meds_demography_patch_fusefiss, only : apply_patch_disturbance, new_fuse_patches, terminate_patches, sort_patches
   use meds_demography_rates,    only : npp_to_growth, camac_mortality, decline_to_floor_mortality,  &
                                        npp_to_recruitment
   use meds_plant_trait_dynamics, only : light_plastic_traits, update_plastic_trait
   use meds_phenology_types, only : pheno_env_t, pheno_args_t, pheno_state_t, pheno_out_t
   use meds_phenology, only : phenology_kernel, leaf_turnover_step, turnover_gate
   use meds_plant_carbon_allocation, only : plant_carbon_allocation
   use meds_sink_limitation,      only : growth_sink_limitation
   use meds_litter_partition, only : add_necromass
   use meds_biogeochem_types, only : litter_input_t
   use meds_diag_vars,            only : CO_DDBH_DT, CO_DAGB_DT, CO_MORT_RATE, CO_NPP_LEAF, CO_NPP_FINEROOT, CO_NPP_WOOD, &
                                         CO_NPP_STORAGE, CO_NPP_REPRO, CO_GROWTH_RESP, CO_STORAGE_RESP, CO_EXUDATE, &
                                         expand_cohort_diag, PA_LITTER_LEAF, PA_LITTER_FINEROOT, PA_LITTER_STRUCT, PA_LITTER_WOOD, &
                                         PA_LITTER_REPRO, PA_MORT_C_BACKGROUND, PA_RECRUIT_NPLANT, reset_cohort_diag
   implicit none
   private

   public :: vegetation_dynamics, restructure_stand, advance_leaf_phenology, advance_plant_traits
   public :: reacclimate_plant_traits
   public :: shed_turnover_water, shed_excess_tissue_water, fill_new_leaf_water, accumulate_recruit_pool
   public :: update_cohort_derivatives
   public :: pool_change_block

   !----- Each cohort's net change of its four carbon pools over one slow step [kgC/plant/step]:    !
   !      growth less shed for leaf and fine root, growth for wood, and allocation less storage    !
   !      upkeep plus resorption for storage. Not the reported NPP (the gross CO_NPP_* rows): a     !
   !      canopy that sheds more than it grows has a negative leaf change, the sign               !
   !      shed_turnover_water reads. Public only for the tests. ----------------------------------!
   type :: pool_change_block
      real(wp), allocatable :: leaf(:), fineroot(:), wood(:), nonstructural(:)
   end type pool_change_block

   !----- Patch disturbance runs once a year, so it integrates its yearly rate over this interval. !
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !< [yr]

contains

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
      real(wp), allocatable    :: mortality(:), recruitment(:,:), npp_repro(:)
      logical,  allocatable    :: exhausted(:)
      real(wp)                 :: mort_water
      real(wp)                 :: tissue_heat0, tissue_heat1, handover, plastic_heat0
      real(wp), allocatable    :: nplant_before(:), leaf_carbon_before(:)
      real(wp)                 :: seed_rain_carbon
      integer(ik)              :: ip
      type(pool_change_block)  :: dpool
      integer(ik)              :: n_window

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

      !----- 1. Carbon allocation, which adds the turnover litter to the per-patch accumulator. Zero  !
      !         that first: it is patch state, so it arrives holding the last step's litter. ------!
      poly%patch%litter_in(1:poly%patch%n) = litter_input_t()
      !----- Size and zero the slow cohort block (one sample per step) before anything writes to it. !
      if (poly%cohort%sdiag%active) then
         call expand_cohort_diag(poly%cohort%sdiag, max(poly%cohort%n, 1_ik))
         call reset_cohort_diag(poly%cohort%sdiag)
         poly%cohort%sdiag%n = poly%cohort%n
      end if
      call compute_carbon_allocation(poly, cfg, cfg%dt_years, dpool, npp_repro, poly%patch%litter_in, exhausted)

      !----- 1b. Shed the turnover's tissue water before the commit, which changes the pools the    !
      !          loss is measured against (docs/science/soil_plant_water.md §2). ------------------!
      call shed_turnover_water(poly, cfg, dpool, ledger)
      !----- Every movement of the allocation phase is declared, so a residual here is a movement   !
      !      nobody declared (docs/science/conservation.md §3). -------------------------------------!
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_ALLOCATE)

      !----- 2. The vital rates, before the commit: mortality reads the growth average of the steps !
      !         before this one (docs/science/demography.md §1). ------------------------------------!
      call compute_vital_rates(poly, cfg, dpool%wood, npp_repro, exhausted, cfg%dt_years, mortality, recruitment)

      !----- 2b. Credit the recruit pool with this step's recruitment, and declare the seed rain    !
      !          as carbon in: no parent at the site paid for it (docs/science/demography.md §2). --!
      call accumulate_recruit_pool(poly, cfg, recruitment, cfg%dt_years, seed_rain_carbon)
      if (present(ledger)) call slow_ledger_declare(ledger, carbon_in = seed_rain_carbon)

      !----- 3. Growth and mortality. First, the growth history's slot for this step's sample. ---!
      n_window = growth_window_steps(cfg)
      poly%growth_hist_pos = mod(poly%growth_hist_pos, n_window) + 1_ik

      !----- The driver computes each cohort's tendencies and advances its growth history; the      !
      !      engine's applier commits them (docs/science/demography.md §2). ----------------------!
      call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
      call update_cohort_derivatives(poly, dpool, mortality, cfg%dt_years, n_window, poly%growth_hist_pos)
      !----- The commit and the new leaves' water leave the temperatures alone, so the energy change  !
      !      is the tissues' thermal mass: one exchange, declared below, not a leak                  !
      !      (docs/science/conservation.md §4). ---------------------------------------------------!
      tissue_heat0 = slow_tissue_energy(poly)
      !----- The fast loop's net carbon, valued at the pre-commit density it accumulated on; after  !
      !      the commit it would miss the dying plants' share (docs/science/conservation.md §3). --!
      handover = slow_fast_carbon_handover(poly, cfg)
      if (allocated(nplant_before)) deallocate(nplant_before)
      allocate(nplant_before(max(poly%cohort%n, 1_ik)))
      nplant_before(1:poly%cohort%n) = poly%cohort%nplant(1:poly%cohort%n)
      leaf_carbon_before = poly%cohort%leaf_carbon(1:poly%cohort%n)
      call update_cohort_states(poly%cohort, poly%deriv, cfg%dt_years, cfg%negligible_nplant)
      call fill_new_leaf_water(poly, cfg, leaf_carbon_before)
      tissue_heat1 = slow_tissue_energy(poly)

      !----- Mortality's litter and water, on the density the commit removed and the pools it left  !
      !      (docs/science/demography.md eq. (9)); recorded whether or not soil carbon is on. -----!
      call accumulate_mortality_litter(poly, cfg, nplant_before, poly%patch%litter_in)
      call shed_mortality_water(poly, cfg, nplant_before, mort_water)

      !----- Growth changed heights: re-sort tallest first (docs/science/demography.md §3). -------!
      call sort_cohorts(poly)

      !----- The growth phase declares the fast loop's net carbon (in), the litter and the CO2 owed  !
      !      to the canopy air (out), the dying plants' water and the thermal mass, so it closes to  !
      !      round-off (docs/science/conservation.md §3). ------------------------------------------!
      if (present(ledger)) then
         call slow_ledger_declare(ledger, carbon_in  = handover,                                    &
                                          carbon_out = litter_carbon_total(poly, poly%patch%litter_in)              &
                                                     + slow_co2_handoff(poly, cfg),                &
                                          water_out  = mort_water,                                 &
                                          energy_in  = tissue_heat1 - tissue_heat0)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_GROW)
      end if

      !----- Patch ageing runs after this routine, in meds_slow_dynamics. ------------------------!

      !----- Slow patch rows: litterfall and recruitment, written before the restructuring permutes  !
      !      the patch axis. Each adds its per-year rate times dt_slow in seconds, a litter amount   !
      !      made a rate first; the litter rows and the soil's litter input are one number           !
      !      (docs/science/diagnostics.md §4, "The slow rows"). -----------------------------------!
      if (poly%patch%diag%active) then
         !----- Without the fast loop, the step sets the block's patch count itself. ------------!
         if (.not. cfg%fast_biophysics_on) poly%patch%diag%n = poly%patch%n
         do ip = 1_ik, poly%patch%n
            !----- Litter by the tissue it came from, not by the pool it reaches. ---------------!
            poly%patch%diag%v(PA_LITTER_LEAF,     ip) = poly%patch%diag%v(PA_LITTER_LEAF,     ip) &
                                                      + (poly%patch%litter_in(ip)%leaf_c                &
                                                         / cfg%dt_years) * cfg%dt_slow
            poly%patch%diag%v(PA_LITTER_FINEROOT, ip) = poly%patch%diag%v(PA_LITTER_FINEROOT, ip) &
                                                      + (poly%patch%litter_in(ip)%fineroot_c            &
                                                         / cfg%dt_years) * cfg%dt_slow
            poly%patch%diag%v(PA_LITTER_WOOD,     ip) = poly%patch%diag%v(PA_LITTER_WOOD,     ip) &
                                                      + ((poly%patch%litter_in(ip)%wood_c               &
                                                        + poly%patch%litter_in(ip)%storage_c)           &
                                                         / cfg%dt_years) * cfg%dt_slow
            poly%patch%diag%v(PA_LITTER_REPRO,    ip) = poly%patch%diag%v(PA_LITTER_REPRO,    ip) &
                                                      + (poly%patch%litter_in(ip)%repro_c               &
                                                         / cfg%dt_years) * cfg%dt_slow
            poly%patch%diag%v(PA_LITTER_STRUCT,   ip) = poly%patch%diag%v(PA_LITTER_STRUCT,   ip) &
                                                      + ((poly%patch%litter_in(ip)%struct_grnd          &
                                                        + poly%patch%litter_in(ip)%struct_soil)         &
                                                         / cfg%dt_years) * cfg%dt_slow
            poly%patch%diag%v(PA_RECRUIT_NPLANT,  ip) = poly%patch%diag%v(PA_RECRUIT_NPLANT,  ip) &
                                                      + sum(recruitment(:, ip)) * cfg%dt_slow
            !----- Without the fast loop, the step adds the block's weight itself, or every row reads 0. !
            if (.not. cfg%fast_biophysics_on) poly%patch%diag%w(ip) = poly%patch%diag%w(ip) + cfg%dt_slow
         end do
      end if

      !----- 4. Refresh the overtopping LAI on the sorted stand (docs/science/demography.md §6). --!
      call update_overtopping_lai(poly)
   end subroutine vegetation_dynamics

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
         call new_fuse_cohorts(poly, cfg)
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
         call new_fuse_patches(poly, cfg)
         call terminate_patches(poly, cfg)
         call new_fuse_cohorts(poly, cfg)
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

   !----- The CO2 this step owes the canopy air [kgC/m2 site]: growth respiration and storage upkeep, !
   !       less the starvation deficit. A declared handoff, not a store: it reaches the canopy air   !
   !       in the next fast window (docs/science/conservation.md §3). -------------------------------!
   pure function slow_co2_handoff(poly, cfg) result(c)
      type(polygon_t),     intent(in) :: poly
      type(meds_config_t), intent(in) :: cfg
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, poly%patch%n
         c = c + patch_polygon_area(poly, ip) * poly%patch%slow_co2_rate(ip) * KGC_PER_UMOL_C * cfg%dt_slow
      end do
   end function slow_co2_handoff

   !----- Site-total carbon in this step's litter accumulator [kgC/m2 site], declared leaving at    !
   !       SLOW_PHASE_GROW and arriving at SLOW_PHASE_SOILC; with soil carbon off nothing arrives   !
   !       (docs/science/conservation.md §3). ------------------------------------------------------!
   pure function litter_carbon_total(poly, lit) result(c)
      type(polygon_t),      intent(in) :: poly
      type(litter_input_t), intent(in) :: lit(:)
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, poly%patch%n
         c = c + patch_polygon_area(poly, ip) * (lit(ip)%labile_grnd + lit(ip)%labile_soil          &
                                      + lit(ip)%struct_grnd + lit(ip)%struct_soil)
      end do
   end function litter_carbon_total

   !---------------------------------------------------------------------------------------!
   ! Each cohort's tendencies for the applier, which commits them: the four pools' rates of     !
   ! change, the log-space mortality rate, and the realized diameter rate, entered in the     !
   ! growth history that mortality averages. Nothing is committed here; the diameter comes from  !
   ! the allometry's inverse at the tentative pools (docs/science/demography.md §1-§2).          !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_derivatives(poly, dpool, mortality, dt_yr, n_window, hist_pos)
      type(polygon_t),          intent(inout) :: poly       !< fills poly%deriv and the growth history
      type(pool_change_block),  intent(in)    :: dpool
      real(wp),                 intent(in)    :: mortality(:)
      real(wp),                 intent(in)    :: dt_yr
      integer(ik),              intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: lc_new, fc_new, wc_new, nc_new
      real(wp)    :: dbh_new, height_new, ba_new, agb_new, la_new, dbh_rate

      associate (cohort => poly%cohort)
         do i = 1_ik, cohort%n
            !----- Tentative pools by the applier's own rule (advance_pool): an overdraw stays       !
            !      negative, for the bounds guard to name (docs/science/demography.md eq. (8)). -----!
            lc_new = advance_pool(cohort%leaf_carbon(i),          dpool%leaf(i))
            fc_new = advance_pool(cohort%fineroot_carbon(i),      dpool%fineroot(i))
            wc_new = advance_pool(cohort%wood_carbon(i),          dpool%wood(i))
            nc_new = advance_pool(cohort%nonstructural_carbon(i), dpool%nonstructural(i))
            !----- The tentative geometry, for the reported growth rates; update_cohort_states     !
            !      derives the same geometry from the committed pools. -----------------------------!
            call carbon_to_structure(cohort%pft_table%allom(cohort%pft(i)), wc_new, lc_new,          &
                                     cohort%pft_table%wood_density(cohort%pft(i)), cohort%sla(i),   &
                                     dbh_new, height_new, ba_new, agb_new, la_new)
            !----- The pools' tendencies, and the realized diameter rate into the growth history. -!
            dbh_rate = (dbh_new - cohort%dbh(i)) / dt_yr
            call fill_cohort_deriv(cohort, i, poly%deriv, dt_yr, mortality(i), dbh_rate,            &
                                   lc_new, fc_new, wc_new, nc_new, n_window, hist_pos)
            !----- Copy the reported rates into the slow cohort block, which follows the cohort     !
            !      reordering; poly%deriv does not (docs/science/diagnostics.md §8, "The one trap"). !
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CO_DDBH_DT,   i) = dbh_rate
               cohort%sdiag%v(CO_DAGB_DT,   i) = (agb_new - cohort%agb(i)) / dt_yr
               cohort%sdiag%v(CO_MORT_RATE, i) = mortality(i)
               cohort%sdiag%w(i)               = 1.0_wp
            end if
         end do
      end associate
   end subroutine update_cohort_derivatives

   !---------------------------------------------------------------------------------------!
   ! The vital rates from the step's carbon: each cohort's mortality [1/yr] from its wood growth, !
   ! and each PFT's recruitment per patch [plant/m2/yr] from the reproduction carbon plus the     !
   ! seed rain, by the laws of meds_demography_rates (docs/science/demography.md §1). The sweep   !
   ! over the stand stays here, so the laws keep no state.                                        !
   !---------------------------------------------------------------------------------------!
   subroutine compute_vital_rates(poly, cfg, npp_wood, npp_repro, exhausted, dt_yr, mortality, recruitment)
      type(polygon_t),       intent(in)  :: poly
      type(meds_config_t),   intent(in)  :: cfg
      real(wp),              intent(in)  :: npp_wood(:), npp_repro(:), dt_yr
      logical,               intent(in)  :: exhausted(:)  !< no leaves and no storage (compute_carbon_allocation)
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr] per cohort
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, j, pf, ip
      real(wp)    :: dbh_rate
      real(wp), allocatable :: carbon_min(:), efficiency(:)

      n    = poly%cohort%n
      np   = poly%patch%n
      npft = cfg%pft%n
      allocate(mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp
      efficiency  = cohort_recruit_efficiency(poly, cfg)

      associate (cohort => poly%cohort, pft => cfg%pft)
         allocate(carbon_min(npft))
         do pf = 1_ik, npft
            carbon_min(pf) = min_cohort_carbon(cohort%pft_table%allom(pf), pft%min_cohort_height,      &
                                               pft%wood_density(pf))
         end do
         !----- The baseline seed rain, from outside the site. ----------------------------!
         do ip = 1_ik, np
            do pf = 1_ik, npft
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%seed_rain_recruits(pf)
            end do
         end do
         do j = 1_ik, n
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !----- Prospective diameter rate and the Camac hazard (docs/science/demography.md eqs. (1), (3)). !
            dbh_rate = npp_to_growth(cohort%pft_table%allom(pf), cohort%wood_carbon(j), npp_wood(j),  &
                                     cohort%dbh(j), dt_yr, pft%wood_density(pf))
            mortality(j) = camac_mortality(cohort%growth_avg(j), dbh_rate,                          &
                                           cohort%growth_count(j) > 0_ik,                           &
                                           pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))
            !----- A cohort with no leaves and no storage cannot flush again: its density goes to the  !
            !      floor this step (docs/science/demography.md eq. (4); plant_carbon_allocation.md §2). !
            if (exhausted(j)) mortality(j) = max(mortality(j),                                      &
                    decline_to_floor_mortality(cohort%nplant(j), cfg%negligible_nplant, dt_yr))
            !----- Reproduction carbon to recruits, a yearly rate the recruit pool takes daily; a  !
            !      disabled PFT recruits nothing (docs/science/demography.md eq. (6)). ------------!
            if (pft%include_pft(pf) == 1_ik)                                                       &
            recruitment(pf, ip) = recruitment(pf, ip)                                              &
                 + npp_to_recruitment(cohort%nplant(j), npp_repro(j), dt_yr,                       &
                                      efficiency(j), carbon_min(pf))
         end do
      end associate
   end subroutine compute_vital_rates

   !---------------------------------------------------------------------------------------!
   ! Each cohort's recruit carbon efficiency this step: the share of its reproduction carbon   !
   ! that becomes recruits, reduced by the shade of its patch's LAI (docs/science/demography.md !
   ! eq. (5)). The recruits and the seed lost to litter both take it from here, so the two     !
   ! always add up to the reproduction carbon.                                                 !
   !---------------------------------------------------------------------------------------!
   function cohort_recruit_efficiency(poly, cfg) result(efficiency)
      type(polygon_t),     intent(in) :: poly
      type(meds_config_t), intent(in) :: cfg
      real(wp), allocatable :: efficiency(:)
      real(wp), allocatable :: patch_lai(:)
      integer(ik) :: j, pf, ip
      allocate(efficiency(poly%cohort%n), patch_lai(max(poly%patch%n, 1_ik)))
      patch_lai = 0.0_wp
      associate (cohort => poly%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            ip = cohort%owner_patch(j)
            patch_lai(ip) = patch_lai(ip) + cohort%nplant(j) * cohort%leaf_area(j)
         end do
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            efficiency(j) = pft%recruit_carbon_efficiency(pf)                                      &
                          * exp(-pft%recruit_shade_decay(pf) * patch_lai(cohort%owner_patch(j)))
         end do
      end associate
   end function cohort_recruit_efficiency

   !---------------------------------------------------------------------------------------!
   ! Each cohort's carbon for the step: gather the day's GPP and maintenance respiration, apply  !
   ! turnover first (this step's litter) so the growth demands are formed against the post-shed !
   ! pools, allocate with plant_carbon_allocation, and form the net per-pool change. Turnover   !
   ! first lets a canopy that is not senescing replace its loss within the step and hold its    !
   ! target. The litter and the CO2 owed go to the patch                                         !
   ! (docs/science/plant_carbon_allocation.md §1-§2, §4-§5).                                     !
   !---------------------------------------------------------------------------------------!
   subroutine compute_carbon_allocation(poly, cfg, dt_yr, dpool, npp_repro, lit, exhausted)
      type(polygon_t),         intent(inout) :: poly   !< no pool moves here; inout for the diagnostics and patch rates
      type(meds_config_t),     intent(in)  :: cfg
      real(wp),                intent(in)  :: dt_yr
      type(pool_change_block), intent(out) :: dpool
      real(wp), allocatable,   intent(out) :: npp_repro(:)
      type(litter_input_t),    intent(inout) :: lit(:)   !< per-patch litter accumulator
      logical,  allocatable,   intent(out) :: exhausted(:)  !< per cohort: no leaves and no storage
      integer(ik) :: n, j, pf, ip
      real(wp)    :: leaf_target, gross_gpp, resp_maint, dt_day, r2l, leaf_turn, root_turn
      real(wp)    :: leaf_demand, fineroot_demand, store_demand, repro_frac
      real(wp)    :: leaf_shed_c, fineroot_shed_c
      real(wp)    :: g_leaf, g_fineroot, g_wood, npp_store, g_repro, growth_resp, deficit
      real(wp)    :: storage_maint, store_post, leaf_retained_c, fineroot_retained_c
      real(wp)    :: seed_lost, exudate, window_yr, grown_in_window
      real(wp), allocatable :: co2_owed(:), efficiency(:)
      logical     :: starving
      integer(ik) :: n_window, next_pos

      n = poly%cohort%n
      dt_day = cfg%dt_slow / day_sec            ! days in the slow step (phenology rates are [1/day])
      allocate(dpool%leaf(n), dpool%fineroot(n), dpool%wood(n), dpool%nonstructural(n), npp_repro(n), exhausted(n))
      allocate(co2_owed(poly%patch%n)) ; co2_owed = 0.0_wp
      efficiency = cohort_recruit_efficiency(poly, cfg)
      !----- The growth-history slot this step's sample will fill (vegetation_dynamics advances    !
      !      growth_hist_pos the same way): the sink limit's window leaves out its old sample. ----!
      n_window = growth_window_steps(cfg)
      next_pos = mod(poly%growth_hist_pos, n_window) + 1_ik
      associate (cohort => poly%cohort, pft => cfg%pft)
         do j = 1_ik, n
            pf  = cohort%pft(j)
            r2l = pft%root_to_leaf_ratio(pf)
            !----- The leaf target at the cohort's own SLA conserves leaf area (docs/science/plant_traits.md §4). !
            leaf_target  = size_to_leaf_carbon(cohort%pft_table%allom(pf), cohort%dbh(j), cohort%height(j), &
                                            cohort%sla(j))
            leaf_turn    = cohort%leaf_turnover_rate(j)   ! baseline leaf turnover [1/yr]
            !----- Storage upkeep, charged first: netted into the storage tendency, not written to the  !
            !      pool, and the allocator sees the reserve after it. The fraction is clamped to [0,1],  !
            !      so the charge never exceeds the pool or turns negative                                !
            !      (docs/science/plant_carbon_allocation.md eq. (2), §1 "Storage maintenance"). --------!
            storage_maint = cohort%nonstructural_carbon(j)                                              &
                          * min(1.0_wp, max(0.0_wp, pft%storage_turnover_rate(pf) * dt_yr))
            store_post    = cohort%nonstructural_carbon(j) - storage_maint
            !----- Storage upkeep is owed to the canopy air with growth respiration                 !
            !      (docs/science/plant_respiration.md §6). ----------------------------------------!
            co2_owed(cohort%owner_patch(j)) = co2_owed(cohort%owner_patch(j))                           &
                                            + cohort%nplant(j) * storage_maint

            store_demand = max(0.0_wp, storage_carbon_target(pft%storage_cushion(pf), leaf_target,      &
                                                             r2l * leaf_target) - store_post)
            !----- Reproduction starts at min_repro_height (docs/science/plant_carbon_allocation.md §2). !
            repro_frac   = merge(pft%reproduction_investment_fraction(pf), 0.0_wp,                  &
                                 cohort%height(j) >= pft%min_repro_height(pf))
            !----- Fine-root turnover, slowed with the leaves' when fineroot_turnover_plastic      !
            !      (docs/science/plant_traits.md §4); then the turnover and the growth demands. ----!
            root_turn = pft%fineroot_turnover_rate(pf)
            if (cfg%fineroot_turnover_plastic)                                                     &
               root_turn = root_turn * cohort%leaf_turnover_rate(j) / pft%leaf_turnover_rate_toc(pf)
            call cohort_carbon_demand(cfg%fast_biophysics_on,                                      &
                     cohort%gpp_accum(j), cohort%leaf_resp_accum(j), cohort%stem_resp_accum(j),    &
                     cohort%root_resp_accum(j), cfg%gpp_ref, cohort%leaf_area(j), dt_yr,           &
                     cohort%leaf_flush_tendency(j), cohort%leaf_shed_tendency(j),                  &
                     pft%pheno(pf)%flush_rate_max, pft%pheno(pf)%shed_rate_max, leaf_turn,         &
                     root_turn,                                                                    &
                     pft%pheno(pf)%min_leaf_cover, pft%pheno(pf)%bare_leaf_cover,                  &
                     cohort%leaf_carbon(j), cohort%fineroot_carbon(j), leaf_target, dt_day, r2l,   &
                     gross_gpp, resp_maint, leaf_shed_c, fineroot_shed_c, leaf_demand,             &
                     fineroot_demand, pft%retained_carbon_fraction(pf), leaf_retained_c,           &
                     fineroot_retained_c)
            call growth_window_so_far(cohort%growth_count(j), cohort%growth_accum(j),              &
                                      cohort%growth_hist(next_pos, j), n_window, dt_yr,             &
                                      window_yr, grown_in_window)
            !----- Allocate the day's carbon: growth respiration on realized growth, wood capped by !
            !      the sink limit (docs/science/plant_carbon_allocation.md §2-§3). ---------------!
            call plant_carbon_allocation(gpp=gross_gpp, resp_maint=resp_maint,                     &
                     growth_resp_frac=pft%growth_resp_factor(pf), storage=store_post,                    &
                     leaf_demand=leaf_demand, fineroot_demand=fineroot_demand,                     &
                     storage_demand=store_demand, repro_frac=repro_frac,                           &
                     growth_leaf=g_leaf, growth_fineroot=g_fineroot, growth_wood=g_wood,          &
                     npp_store=npp_store, growth_repro=g_repro, growth_resp=growth_resp,          &
                     deficit=deficit, starving=starving,                                           &
                     wood_growth_max=growth_sink_limitation(                                       &
                        allom=cohort%pft_table%allom(pf),                                            &
                        wood_carbon=cohort%wood_carbon(j), dbh=cohort%dbh(j),                       &
                        wood_density=pft%wood_density(pf),                                          &
                        max_relative_growth_rate=pft%max_relative_growth_rate(pf),                  &
                        max_absolute_growth_rate=pft%max_absolute_growth_rate(pf),                  &
                        max_absolute_growth_exponent=pft%max_absolute_growth_exponent(pf),          &
                        dt_yr=dt_yr, window_yr=window_yr, grown_in_window=grown_in_window),         &
                     exudate=exudate)
            !----- The net per-pool change: growth less the shed for leaf and fine root. ---------!
            dpool%leaf(j)          = g_leaf     - leaf_shed_c
            dpool%fineroot(j)      = g_fineroot - fineroot_shed_c
            dpool%wood(j)          = g_wood
            !----- Storage: allocation less upkeep, plus the resorbed carbon, which dpool%leaf and    !
            !      dpool%fineroot remove in full (docs/science/plant_carbon_allocation.md eq. (3)). --!
            dpool%nonstructural(j) = npp_store - storage_maint + leaf_retained_c + fineroot_retained_c
            npp_repro(j)         = g_repro
            !----- No leaves and no storage: the maintenance left unpaid once storage is spent exceeds !
            !      all the leaf carbon, so no flush can be paid for                                   !
            !      (docs/science/plant_carbon_allocation.md §2, "A cohort with no leaves and no storage"). !
            exhausted(j)         = deficit > cohort%leaf_carbon(j)
            !----- Slow cohort rows: the allocation split, growth respiration, storage upkeep and     !
            !      exudate, as rates [kgC/plant/yr], so the carbon budget closes from the output. ----!
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CO_NPP_LEAF,     j) = g_leaf     / dt_yr
               cohort%sdiag%v(CO_NPP_FINEROOT, j) = g_fineroot / dt_yr
               cohort%sdiag%v(CO_NPP_WOOD,     j) = g_wood     / dt_yr
               cohort%sdiag%v(CO_NPP_STORAGE,  j) = npp_store  / dt_yr
               cohort%sdiag%v(CO_NPP_REPRO,    j) = g_repro    / dt_yr
               cohort%sdiag%v(CO_GROWTH_RESP,  j) = growth_resp / dt_yr
               cohort%sdiag%v(CO_STORAGE_RESP, j) = storage_maint / dt_yr
               cohort%sdiag%v(CO_EXUDATE,      j) = exudate / dt_yr
               cohort%sdiag%w(j)                  = 1.0_wp
            end if
            !----- Growth respiration, less the starvation deficit the fast loop over-reported, is   !
            !      owed to the canopy air (docs/science/plant_respiration.md §6). -------------------!
            ip = cohort%owner_patch(j)
            co2_owed(ip) = co2_owed(ip) + cohort%nplant(j) * (growth_resp - deficit)

            !----- The seed that does not establish is litter that falls with the canopy           !
            !      (docs/science/demography.md §1, "Recruitment"). ------------------------------!
            seed_lost = g_repro * cohort%nplant(j)                                              &
                      * (1.0_wp - min(max(efficiency(j), 0.0_wp), 1.0_wp))
            !----- The unresorbed share of the turnover and the seed lost go to the patch's litter !
            !      (docs/science/plant_carbon_allocation.md eq. (3)). ----------------------------!
            call add_necromass(lit(ip), (leaf_shed_c - leaf_retained_c) * cohort%nplant(j),       &
                     (fineroot_shed_c - fineroot_retained_c) * cohort%nplant(j), 0.0_wp, 0.0_wp,   &
                     seed_lost, pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                       &
                     cohort%pft_table%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf))
            !----- Root exudate, the carbon the sink limit left, is below-ground labile litter.  !
            lit(ip)%labile_soil = lit(ip)%labile_soil + exudate * cohort%nplant(j)
            lit(ip)%exudate_c   = lit(ip)%exudate_c   + exudate * cohort%nplant(j)
         end do
      end associate

      !----- The CO2 owed as a molar rate [umol/m2/s] the next fast window emits. -------------!
      do ip = 1_ik, poly%patch%n
         poly%patch%slow_co2_rate(ip) = co2_owed(ip) / max(cfg%dt_slow, tiny_num) / KGC_PER_UMOL_C
      end do
   end subroutine compute_carbon_allocation

   !---------------------------------------------------------------------------------------!
   ! A cohort's growth window as the sink limit sees it this step: the diameter already grown   !
   ! in it [cm] and its length [yr] with this step in it. Once the history is full, this step's !
   ! sample replaces the oldest; a new cohort's window is as long as its life                   !
   ! (docs/science/plant_carbon_allocation.md eq. (4); demography.md eq. (2)).                  !
   !---------------------------------------------------------------------------------------!
   pure subroutine growth_window_so_far(count, accum, oldest, n_window, dt_yr, window_yr, grown)
      integer(ik), intent(in)  :: count       !< samples in the ring buffer
      real(wp),    intent(in)  :: accum       !< [cm/yr] their sum
      real(wp),    intent(in)  :: oldest      !< [cm/yr] the sample in the slot this step overwrites
      integer(ik), intent(in)  :: n_window    !< the buffer's length [steps]
      real(wp),    intent(in)  :: dt_yr       !< [yr] step length
      real(wp),    intent(out) :: window_yr   !< [yr] the window this step closes
      real(wp),    intent(out) :: grown       !< [cm] diameter grown in it before this step
      if (count < n_window) then
         window_yr = real(count + 1_ik, wp) * dt_yr
         grown     = accum * dt_yr
      else
         window_yr = real(n_window, wp) * dt_yr
         grown     = (accum - oldest) * dt_yr
      end if
   end subroutine growth_window_so_far

   !---------------------------------------------------------------------------------------!
   ! The turnover's tissue water. A net leaf or fine-root loss sheds water in proportion to the !
   ! carbon lost, so the tissue left keeps its water per unit carbon; fine roots shed their      !
   ! share of the wood store. A net gain moves no water: the store keeps its mass. It runs       !
   ! before the commit, against the pre-step pools. The water and its enthalpy go to the ground  !
   ! as patch rates (shed_water_rate, shed_water_enth) the next fast window adds                 !
   ! (docs/science/soil_plant_water.md §2; conservation.md §3). ---------------------------------!
   subroutine shed_turnover_water(poly, cfg, dpool, ledger)
      type(polygon_t),         intent(inout) :: poly
      type(meds_config_t),     intent(in)    :: cfg
      type(pool_change_block), intent(in)    :: dpool
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp), allocatable :: patch_total(:), patch_enth(:)
      real(wp) :: shed_frac, shed_amount
      integer(ik) :: j, ip

      allocate(patch_total(poly%patch%n), patch_enth(poly%patch%n))
      patch_total = 0.0_wp
      patch_enth  = 0.0_wp

      associate (cohort => poly%cohort)
         do j = 1_ik, cohort%n
            ip = cohort%owner_patch(j)
            if (dpool%leaf(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -dpool%leaf(j) / max(cohort%leaf_carbon(j), tiny_num))
               shed_amount = cohort%leaf_water_mass(j) * shed_frac
               cohort%leaf_water_mass(j) = cohort%leaf_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
               patch_enth(ip)  = patch_enth(ip)  + cohort%nplant(j) * shed_amount                        &
                                                 * internal_energy_liquid(cohort%leaf_temp(j))
            end if
            if (dpool%fineroot(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -dpool%fineroot(j) / max(cohort%fineroot_carbon(j), tiny_num))
               !----- The fine roots' own share of the store, in proportion to their carbon. ------!
               shed_amount = max(cohort%wood_water_mass(j), 0.0_wp) * cohort%fineroot_carbon(j)   &
                           / max(cohort_wood_store_carbon(cohort%sapwood_carbon(j), cohort%fineroot_carbon(j)), &
                                 tiny_num) * shed_frac
               cohort%wood_water_mass(j) = cohort%wood_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
               patch_enth(ip)  = patch_enth(ip)  + cohort%nplant(j) * shed_amount                        &
                                                 * internal_energy_liquid(cohort%wood_temp(j))
            end if
            !----- The water leaves at its tissue's temperature: the temperature holds, the energy falls. !
            call hold_cohort_tissue_temps(cohort, j)
         end do
      end associate

      do ip = 1_ik, poly%patch%n
         poly%patch%shed_water_rate(ip) = patch_total(ip) / max(cfg%dt_slow, tiny_num)
         poly%patch%shed_water_enth(ip) = patch_enth(ip)  / max(cfg%dt_slow, tiny_num)
      end do

      !----- The shed water and its enthalpy: a declared handoff (docs/science/conservation.md §3). !
      if (present(ledger)) then
         do ip = 1_ik, poly%patch%n
            call slow_ledger_declare(ledger, water_out  = patch_polygon_area(poly, ip) * patch_total(ip),     &
                                             energy_out = patch_polygon_area(poly, ip) * patch_enth(ip))
         end do
      end if
   end subroutine shed_turnover_water

   !---------------------------------------------------------------------------------------!
   ! New leaves take their water from the wood: the leaf tissue the commit added is filled,  !
   ! with at most the water the wood holds, to the wood's potential less the crown's gravity !
   ! head. The plant's water is unchanged, and a canopy re-flushing from bare does not start  !
   ! at the dry end of its pressure-volume curve. The wood ends a little drier than the       !
   ! leaves at rest; the next fast step evens them out (docs/science/soil_plant_water.md §2). !
   !---------------------------------------------------------------------------------------!
   subroutine fill_new_leaf_water(poly, cfg, leaf_carbon_before)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: leaf_carbon_before(:)   !< [kgC/plant] before the commit
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: grown, head, psi_wood, moved
      integer(ik) :: i

      associate (cohort => poly%cohort)
         do i = 1_ik, cohort%n
            grown = cohort%leaf_carbon(i) - leaf_carbon_before(i)
            if (grown <= 0.0_wp) cycle
            leaf_curve = cfg%pft%hydro(cohort%pft(i))%leaf_curve
            wood_curve = cfg%pft%hydro(cohort%pft(i))%wood_curve
            head     = grav_head * cohort%height(i)          ! the crown above the wood datum [MPa]
            psi_wood = tissue_water_to_psi(cohort%wood_water_mass(i), wood_curve,                           &
                          cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
            moved    = min(tissue_psi_to_water(psi_wood - head, leaf_curve, grown), cohort%wood_water_mass(i))
            cohort%leaf_water_mass(i) = cohort%leaf_water_mass(i) + moved
            cohort%wood_water_mass(i) = cohort%wood_water_mass(i) - moved
            call hold_cohort_tissue_temps(cohort, i)
         end do
      end associate
   end subroutine fill_new_leaf_water

   !---------------------------------------------------------------------------------------!
   ! Tissue water above the capacity a cohort's biomass now allows goes to the ground, down the !
   ! turnover shed's channel (patch%shed_water_rate). An operator that shrinks a store without  !
   ! touching its water (an SLA rise, a fusion or split) can leave it above saturation; the     !
   ! caller declares the water shed (docs/science/soil_plant_water.md §2).                      !
   !---------------------------------------------------------------------------------------!
   subroutine shed_excess_tissue_water(poly, cfg, water_shed)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(out)   :: water_shed     !< [kg/m2 site] the water shed
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: leaf_kept, wood_kept, excess, excess_enth
      integer(ik) :: i, ip

      water_shed = 0.0_wp
      associate (cohort => poly%cohort, patch => poly%patch)
         do i = 1_ik, cohort%n
            leaf_curve = cfg%pft%hydro(cohort%pft(i))%leaf_curve
            wood_curve = cfg%pft%hydro(cohort%pft(i))%wood_curve
            leaf_kept = tissue_capped_water(cohort%leaf_water_mass(i), leaf_curve, cohort%leaf_carbon(i))
            wood_kept = tissue_capped_water(cohort%wood_water_mass(i), wood_curve,                &
                           cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
            excess = (cohort%leaf_water_mass(i) - leaf_kept) + (cohort%wood_water_mass(i) - wood_kept)
            if (excess <= 0.0_wp) cycle
            excess_enth = (cohort%leaf_water_mass(i) - leaf_kept) * internal_energy_liquid(cohort%leaf_temp(i)) &
                        + (cohort%wood_water_mass(i) - wood_kept) * internal_energy_liquid(cohort%wood_temp(i))
            cohort%leaf_water_mass(i) = leaf_kept
            cohort%wood_water_mass(i) = wood_kept
            call hold_cohort_tissue_temps(cohort, i)          ! the excess leaves at the tissue's temperature
            ip = cohort%owner_patch(i)
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip)                                   &
                                      + cohort%nplant(i) * excess / max(cfg%dt_slow, tiny_num)
            patch%shed_water_enth(ip) = patch%shed_water_enth(ip)                                   &
                                      + cohort%nplant(i) * excess_enth / max(cfg%dt_slow, tiny_num)
            water_shed = water_shed + patch_polygon_area(poly, ip) * cohort%nplant(i) * excess
         end do
      end associate
   end subroutine shed_excess_tissue_water

   !---------------------------------------------------------------------------------------!
   ! The water of the plants that die this step: their tissue water and canopy films go to the !
   ! ground down the turnover shed's channel. Runs after the commit, with the dying number      !
   ! n0 - n1 taken as accumulate_mortality_litter takes it (docs/science/demography.md §2).      !
   !---------------------------------------------------------------------------------------!
   subroutine shed_mortality_water(poly, cfg, nplant_before, water_shed)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: nplant_before(:) !< [plant/m2] density before the commit
      real(wp),            intent(out)   :: water_shed       !< [kg/m2 site]
      integer(ik) :: j, ip
      real(wp)    :: died_nplant, w_j

      water_shed = 0.0_wp
      associate (cohort => poly%cohort, patch => poly%patch)
         do j = 1_ik, cohort%n
            died_nplant = nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            ip  = cohort%owner_patch(j)
            !----- The density drop times the per-plant water (cohort_tissue_water would carry the !
            !      post-commit density and scale the loss by n1/n0). -----------------------------!
            w_j = died_nplant * (max(cohort%leaf_water_mass(j), 0.0_wp)                            &
                               + max(cohort%wood_water_mass(j), 0.0_wp)                            &
                               + cohort%leaf_surf_water(j) + cohort%wood_surf_water(j))
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + w_j / max(cfg%dt_slow, tiny_num)
            patch%shed_water_enth(ip) = patch%shed_water_enth(ip)                                   &
                                      + died_nplant * cohort_water_enthalpy(cohort, j) / max(cfg%dt_slow, tiny_num)
            water_shed = water_shed + patch_polygon_area(poly, ip) * w_j
         end do
      end associate
   end subroutine shed_mortality_water

   !---------------------------------------------------------------------------------------!
   ! Credit one step's recruitment into each (PFT, patch) recruit pool, and report the seed     !
   ! rain, the part no parent at the site paid for. The parents' share is debited and credited  !
   ! in the same phase, so it is internal. include_pft gates the seed rain as it gates the rate  !
   ! (docs/science/demography.md §2; conservation.md §3).                                        !
   !---------------------------------------------------------------------------------------!
   subroutine accumulate_recruit_pool(poly, cfg, recruitment, dt_yr, seed_rain_carbon)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      real(wp),            intent(in)    :: dt_yr
      real(wp),            intent(out)   :: seed_rain_carbon  !< [kgC/m2 site] the external part
      integer(ik) :: ip, pf
      real(wp)    :: c_min

      seed_rain_carbon = 0.0_wp
      associate (patch => poly%patch, pft => cfg%pft)
         do ip = 1_ik, patch%n
            do pf = 1_ik, poly%n_pft
               patch%recruit_pool(pf, ip) = patch%recruit_pool(pf, ip) + recruitment(pf, ip) * dt_yr
               if (pft%include_pft(pf) /= 1_ik) cycle
               c_min = min_cohort_carbon(poly%cohort%pft_table%allom(pf), pft%min_cohort_height,          &
                                         pft%wood_density(pf))
               seed_rain_carbon = seed_rain_carbon                                                 &
                                + patch_polygon_area(poly, ip) * pft%seed_rain_recruits(pf) * dt_yr * c_min
            end do
         end do
      end associate
   end subroutine accumulate_recruit_pool

   !---------------------------------------------------------------------------------------!
   ! Mortality litter: every pool of the plants that died this step, (n0 - n1).p1, on the pools  !
   ! and the density drop of the commit it follows, so litter is valued on the new pools and a   !
   ! cohort the floor stopped gets none (docs/science/demography.md eq. (9)). It goes to the same !
   ! per-patch accumulator as the turnover litter.                                                !
   !---------------------------------------------------------------------------------------!
   subroutine accumulate_mortality_litter(poly, cfg, nplant_before, lit)
      type(polygon_t),       intent(inout) :: poly      !< inout for the mortality-carbon diagnostic
      type(meds_config_t),   intent(in)    :: cfg
      real(wp),              intent(in)    :: nplant_before(:) !< [plant/m2] density before the commit
      type(litter_input_t),  intent(inout) :: lit(:)
      integer(ik) :: j, pf, ip
      real(wp)    :: died_nplant

      associate (cohort => poly%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            died_nplant = nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !----- The background mortality's carbon, from the same numbers as the litter below, as  !
            !      a rate [kgC/m2/yr] times the step's seconds (docs/science/diagnostics.md §4). -----!
            if (poly%patch%diag%active)                                                            &
               poly%patch%diag%v(PA_MORT_C_BACKGROUND, ip) =                                       &
                  poly%patch%diag%v(PA_MORT_C_BACKGROUND, ip)                                      &
                  + (died_nplant * (cohort%leaf_carbon(j) + cohort%fineroot_carbon(j)              &
                                  + cohort%wood_carbon(j) + cohort%nonstructural_carbon(j))        &
                     / cfg%dt_years) * cfg%dt_slow
            call add_necromass(lit(ip), died_nplant * cohort%leaf_carbon(j),                      &
                     died_nplant * cohort%fineroot_carbon(j), died_nplant * cohort%wood_carbon(j),  &
                     died_nplant * cohort%nonstructural_carbon(j), 0.0_wp,                         &
                     pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                                 &
                     cohort%pft_table%allom(pf)%aboveground_frac, pft%struct_lignin_frac(pf))
         end do
      end associate
   end subroutine accumulate_mortality_litter

   !---------------------------------------------------------------------------------------!
   ! One cohort's carbon demand: gather the step's GPP and maintenance respiration (the fast    !
   ! loop's accumulators, or the gpp_ref stub), apply turnover first (leaf_turnover_step for    !
   ! the leaves; the fine roots at their own rate under the leaves' gate, not dropping with      !
   ! them), and form the flush-capped growth demands from the post-shed pools                    !
   ! (docs/science/plant_phenology.md eq. (9) and §5).                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine cohort_carbon_demand(fast_biophysics_on,                                       &
            gpp_accum, leaf_resp_accum, stem_resp_accum, root_resp_accum, gpp_ref, leaf_area, dt_yr, &
            leaf_flush_tendency, leaf_shed_tendency, flush_rate_max, shed_rate_max, leaf_turn,     &
            fineroot_turnover_rate, min_leaf_cover, bare_leaf_cover,                               &
            leaf_carbon, fineroot_carbon, leaf_target, dt_day, r2l,                                &
            gross_gpp, resp_maint, leaf_shed_c, fineroot_shed_c, leaf_demand, fineroot_demand,     &
            retained_frac, leaf_retained_c, fineroot_retained_c)
      logical,  intent(in)  :: fast_biophysics_on
      real(wp), intent(in)  :: gpp_accum, leaf_resp_accum, stem_resp_accum, root_resp_accum
      real(wp), intent(in)  :: gpp_ref, leaf_area, dt_yr
      real(wp), intent(in)  :: leaf_flush_tendency, leaf_shed_tendency  !< [-]
      real(wp), intent(in)  :: flush_rate_max, shed_rate_max            !< [1/day]
      real(wp), intent(in)  :: leaf_turn, fineroot_turnover_rate        !< [1/yr]
      real(wp), intent(in)  :: min_leaf_cover, bare_leaf_cover          !< [-]
      real(wp), intent(in)  :: leaf_carbon, fineroot_carbon, leaf_target, dt_day, r2l
      real(wp), intent(out) :: gross_gpp, resp_maint, leaf_shed_c, fineroot_shed_c
      real(wp), intent(out) :: leaf_demand, fineroot_demand
      real(wp), intent(in)  :: retained_frac        !< [-] share of the loss resorbed; the shed stays the full loss
      real(wp), intent(out) :: leaf_retained_c      !< [kgC/plant] leaf carbon resorbed this step
      real(wp), intent(out) :: fineroot_retained_c  !< [kgC/plant] fine-root carbon resorbed this step
      real(wp) :: senescence, background, flush_cap, leaf_post, root_post, root_pool, retained

      if (fast_biophysics_on) then
         gross_gpp  = gpp_accum
         resp_maint = leaf_resp_accum + stem_resp_accum + root_resp_accum
      else
         gross_gpp  = gpp_ref * leaf_area * dt_yr
         resp_maint = 0.0_wp
      end if
      !----- Turnover first: the leaf and fine-root losses, then the post-shed pools. ---------!
      call leaf_turnover_step(leaf_carbon, leaf_target, leaf_flush_tendency, leaf_shed_tendency,    &
                              flush_rate_max, shed_rate_max, leaf_turn, min_leaf_cover,              &
                              bare_leaf_cover, dt_day, senescence, background, flush_cap)
      leaf_shed_c     = senescence + background
      root_pool       = max(fineroot_carbon, 0.0_wp)
      fineroot_shed_c = min(max(fineroot_turnover_rate, 0.0_wp) / yr_day                            &
                            * turnover_gate(leaf_flush_tendency, leaf_shed_tendency) * root_pool * dt_day, &
                            root_pool)
      leaf_post = leaf_carbon     - leaf_shed_c
      root_post = fineroot_carbon - fineroot_shed_c
      retained            = min(max(retained_frac, 0.0_wp), 1.0_wp)
      leaf_retained_c     = retained * leaf_shed_c
      fineroot_retained_c = retained * fineroot_shed_c
      !----- Flush-capped growth demands toward the target, from the post-shed pools. --------!
      leaf_demand     = min(max(0.0_wp, leaf_target       - leaf_post), flush_cap)
      fineroot_demand = min(max(0.0_wp, r2l * leaf_target - root_post), flush_cap * r2l)
   end subroutine cohort_carbon_demand

   !---------------------------------------------------------------------------------------!
   ! Advance every cohort's leaf phenology one slow step (docs/science/plant_phenology.md): build  !
   ! the daily drivers from the fast loop's accumulators (each site's air temperature, each       !
   ! cohort's PAR hours) and the cohort's predawn leaf psi, run the kernel, and write back the     !
   ! tendencies and cue memory. No carbon moves here. A step with no fast sub-steps is skipped,   !
   ! so the memory never advances on a 0/0 mean.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine advance_leaf_phenology(poly, cfg, doy, latitude_deg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: doy
      real(wp),            intent(in), optional :: latitude_deg   !< the polygon's [deg] (default: the [site] latitude)
      type(pheno_env_t)    :: env
      type(pheno_args_t)   :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: i, pf, is
      real(wp)    :: dt_days, temp_day(poly%site%n), lat
      logical     :: north

      if (all(poly%site%pheno_tair_n(1:poly%site%n) < 1_ik)) return   ! no fast sub-steps -> no drivers
      dt_days  = cfg%dt_slow / day_sec
      !----- Each site's daily-mean air temperature: its own forcing after its terrain. -----------!
      do is = 1_ik, poly%site%n
         temp_day(is) = poly%site%pheno_tair_sum(is) / real(max(poly%site%pheno_tair_n(is), 1_ik), wp)
      end do
      lat = cfg%forcing%latitude
      if (present(latitude_deg)) lat = latitude_deg
      north    = lat >= 0.0_wp

      do i = 1_ik, poly%cohort%n
         pf = poly%cohort%pft(i)
         call flatten_pheno_params(cfg, pf, params)

         !----- The daily drivers. dmax_psi_leaf is the last completed predawn value; a cohort born !
         !      today reads 0, well watered, for one day, as in the leaf kernel. ---------------!
         env%temp_day         = temp_day(poly%patch%owner_site(poly%cohort%owner_patch(i)))
         !----- Hours the PAR at the cohort's top exceeded par_min over the step, per day. -----!
         env%par_hours        = poly%cohort%light_hours_accum(i) / dt_days
         env%predawn_leaf_psi = min(0.0_wp, poly%cohort%dmax_psi_leaf(i))
         env%doy              = doy
         env%hemis_north      = north

         !----- Pack the cohort's tendencies and cue memory, advance, unpack. ------------------!
         state%leaf_flush_tendency = poly%cohort%leaf_flush_tendency(i)
         state%leaf_shed_tendency  = poly%cohort%leaf_shed_tendency(i)
         state%growing_degree_days = poly%cohort%growing_degree_days(i)
         state%cold_degree_days    = poly%cohort%cold_degree_days(i)
         state%wet_psi_sum         = poly%cohort%wet_psi_sum(i)
         state%dry_psi_sum         = poly%cohort%dry_psi_sum(i)
         state%light_hours_mean    = poly%cohort%light_hours_mean(i)
         call phenology_kernel(env, params, dt_days, state, out)
         poly%cohort%leaf_flush_tendency(i) = state%leaf_flush_tendency
         poly%cohort%leaf_shed_tendency(i)  = state%leaf_shed_tendency
         poly%cohort%growing_degree_days(i) = state%growing_degree_days
         poly%cohort%cold_degree_days(i)    = state%cold_degree_days
         poly%cohort%wet_psi_sum(i)         = state%wet_psi_sum
         poly%cohort%dry_psi_sum(i)         = state%dry_psi_sum
         poly%cohort%light_hours_mean(i)    = state%light_hours_mean
      end do
   end subroutine advance_leaf_phenology

   !---------------------------------------------------------------------------------------!
   ! Light trait plasticity over one slow step: each cohort's sla, vcmax25, rd25 and leaf        !
   ! turnover rate relax toward their targets at the LAI above, as fast as leaves are replaced,  !
   ! and leaf carbon above the new leaf target is resorbed to storage. Called only with          !
   ! trait_plasticity_on (docs/science/plant_traits.md eqs. (1)-(2), §4).                        !
   !---------------------------------------------------------------------------------------!
   subroutine advance_plant_traits(poly, cfg, dt_yr, instantaneous)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      logical, optional,   intent(in)    :: instantaneous   !< census restart: jump straight to the target
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, turnover_t, turnover_now, leaf_target, excess
      logical     :: instant
      instant = .false.
      if (present(instantaneous)) instant = instantaneous
      associate (cohort => poly%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            !----- Light-acclimated targets from this cohort's cumulative LAI above. -----------!
            call light_plastic_traits(cohort%overtopping_lai(j),                                 &
                     pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_turnover_rate_toc(pf),  &
                     pft%kplastic_sla(pf), pft%kplastic_vcmax(pf), pft%kplastic_rd(pf),             &
                     pft%kplastic_leaf_turnover(pf), sla_t, vcmax_t, rd_t, turnover_t)
            !----- Each trait toward its target, or straight to it on a census start. ----------!
            turnover_now = cohort%leaf_turnover_rate(j)
            cohort%sla(j)     = update_plastic_trait(cohort%sla(j),     sla_t,   turnover_now, dt_yr, instant)
            cohort%vcmax25(j) = update_plastic_trait(cohort%vcmax25(j), vcmax_t, turnover_now, dt_yr, instant)
            cohort%rd25(j)    = update_plastic_trait(cohort%rd25(j),    rd_t,    turnover_now, dt_yr, instant)
            cohort%leaf_turnover_rate(j) = update_plastic_trait(cohort%leaf_turnover_rate(j), turnover_t,      &
                                                                turnover_now, dt_yr, instant)
            !----- A higher SLA lowers the leaf target: resorb the excess leaf carbon to storage, so !
            !      leaf area stays at allometry; a lower SLA leaves a deficit growth fills. --------!
            leaf_target = size_to_leaf_carbon(cohort%pft_table%allom(cohort%pft(j)), cohort%dbh(j), cohort%height(j), &
                                           cohort%sla(j))
            if (cohort%leaf_carbon(j) > leaf_target) then
               excess = cohort%leaf_carbon(j) - leaf_target
               cohort%leaf_carbon(j)          = leaf_target
               cohort%nonstructural_carbon(j) = cohort%nonstructural_carbon(j) + excess
            end if
            call set_cohort_size_from_carbon(cohort, j)   ! leaf_area = sla*leaf_carbon, either way
            call hold_cohort_tissue_temps(cohort, j)      ! the resorbed leaf carbon left at the leaf's temperature
         end do
      end associate
   end subroutine advance_plant_traits

   !---------------------------------------------------------------------------------------!
   ! reacclimate_plant_traits -- a restart's plastic leaf traits from this run's PFT file      !
   ! ([init].reacclimate_traits): the census start's instant acclimation with plasticity on,   !
   ! the top-of-canopy values without. Leaf area is kept: leaf carbon scales by the SLA's       !
   ! change, storage taking up the difference, and a lower SLA takes at most the storage there  !
   ! is, so no carbon is made (docs/configuration.md, "How a run starts"; plant_traits.md §5).  !
   !---------------------------------------------------------------------------------------!
   subroutine reacclimate_plant_traits(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, turnover_t, leaf_carbon_new, sla_old
      associate (cohort => poly%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            if (cfg%trait_plasticity_on) then
               call light_plastic_traits(cohort%overtopping_lai(j),                              &
                        pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_turnover_rate_toc(pf), &
                        pft%kplastic_sla(pf), pft%kplastic_vcmax(pf), pft%kplastic_rd(pf),          &
                        pft%kplastic_leaf_turnover(pf), sla_t, vcmax_t, rd_t, turnover_t)
            else
               sla_t = pft%sla(pf) ; vcmax_t = pft%vcmax25(pf) ; rd_t = pft%rd25(pf)
               turnover_t = pft%leaf_turnover_rate_toc(pf)
            end if
            sla_old           = cohort%sla(j)
            cohort%sla(j)     = sla_t
            cohort%vcmax25(j) = vcmax_t
            cohort%rd25(j)    = rd_t
            cohort%leaf_turnover_rate(j) = turnover_t
            if (cohort%sla(j) == sla_old) cycle
            leaf_carbon_new = cohort%leaf_carbon(j) * sla_old / max(cohort%sla(j), tiny_num)
            leaf_carbon_new = min(leaf_carbon_new, cohort%leaf_carbon(j) + cohort%nonstructural_carbon(j))
            cohort%nonstructural_carbon(j) = cohort%nonstructural_carbon(j) + cohort%leaf_carbon(j)        &
                                             - leaf_carbon_new
            cohort%leaf_carbon(j) = leaf_carbon_new
            call set_cohort_size_from_carbon(cohort, j)
            call hold_cohort_tissue_temps(cohort, j)
         end do
      end associate
   end subroutine reacclimate_plant_traits

   !---------------------------------------------------------------------------------------!
   ! Copy one PFT's phenology traits (cfg%pft%pheno) into the kernel's self-contained parameter  !
   ! set, as meds_pft_derive does for the leaf traits.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine flatten_pheno_params(cfg, ipft, p)
      type(meds_config_t),  intent(in)  :: cfg
      integer(ik),          intent(in)  :: ipft
      type(pheno_args_t),   intent(out) :: p
      associate (t => cfg%pft)
         p%flush_cue_mask        = t%pheno(ipft)%flush_cue_mask
         p%shed_cue_mask         = t%pheno(ipft)%shed_cue_mask
         p%flush_cue_timescale   = t%pheno(ipft)%flush_cue_timescale
         p%shed_cue_timescale    = t%pheno(ipft)%shed_cue_timescale
         p%flush_rate_max        = t%pheno(ipft)%flush_rate_max
         p%shed_rate_max         = t%pheno(ipft)%shed_rate_max
         p%flush_base_temp       = t%pheno(ipft)%flush_base_temp
         p%flush_degree_days     = t%pheno(ipft)%flush_degree_days
         p%flush_temp_sharpness  = t%pheno(ipft)%flush_temp_sharpness
         p%shed_base_temp        = t%pheno(ipft)%shed_base_temp
         p%shed_degree_days      = t%pheno(ipft)%shed_degree_days
         p%shed_temp_sharpness   = t%pheno(ipft)%shed_temp_sharpness
         p%par_min               = t%pheno(ipft)%par_min
         p%flush_light_hours     = t%pheno(ipft)%flush_light_hours
         p%flush_light_sharpness = t%pheno(ipft)%flush_light_sharpness
         p%shed_light_hours      = t%pheno(ipft)%shed_light_hours
         p%shed_light_sharpness  = t%pheno(ipft)%shed_light_sharpness
         p%light_window          = t%pheno(ipft)%light_window
         p%flush_water_sum       = t%pheno(ipft)%flush_water_sum
         p%flush_water_sharpness = t%pheno(ipft)%flush_water_sharpness
         p%shed_water_sum        = t%pheno(ipft)%shed_water_sum
         p%shed_water_sharpness  = t%pheno(ipft)%shed_water_sharpness
      end associate
      !----- The turgor-loss point is derived from the leaf's pressure-volume curve, the threshold  !
      !      the leaf's water-stress controls use too (docs/science/plant_phenology.md §2;          !
      !      leaf_gas_exchange.md §4.1-§4.2). --------------------------------------------------!
      p%leaf_psi_tlp = pft_leaf_psi_tlp(cfg%pft, ipft)
   end subroutine flatten_pheno_params

end module meds_vegetation_dynamics
