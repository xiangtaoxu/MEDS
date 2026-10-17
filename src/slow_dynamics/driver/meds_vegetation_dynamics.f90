! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_vegetation_dynamics -- MEDS's slow-loop VEGETATION-DYNAMICS driver (the analogue of      !
! ED2's veg_dynamics_driver), in src/slow_dynamics/driver/ (library meds_slow). It is THE         !
! ORCHESTRATOR: it (1) assembles the per-cohort carbon NPP from the plant seam, (2) turns that     !
! into the demographic RATES by calling the per-individual vital-rate LAWS                         !
! (meds_demography_rates: carbon growth rate, Camac mortality, recruitment), and (3) sequences     !
! the demography engine's APPLY-PRIMITIVES (grow / die / recruit / fuse-fiss), folding the         !
! calendar cadence into the structural triggers (docs/science/order_of_processes.md §4-§5).        !
!                                                                                          !
! This is THE single place plant and demography meet: the driver calls the plant kernels and       !
! the rate laws; the demography engine only APPLIES supplied arrays (src/README.md rule 4). The    !
! empirical growth and recruitment laws live in the Python example, so the Fortran model is the    !
! carbon path.                                                                                     !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : day_sec, tiny_num, cp_liq, yr_day, grav_head
   use meds_config,               only : meds_config_t, growth_window_steps, pft_leaf_psi_tlp, pft_tissue_curves
   use meds_water_retention,      only : water_curve_t, water_content, clamp_water_to_capacity, shared_potential
   use meds_allometry,            only : size2leaf_carbon, carbon_to_structure, min_cohort_carbon,      &
                                         storage_carbon_target, local_hmax
   use meds_site_state_types,      only : carbon_flux_block, cohort_deriv_alloc, GROWTH_AVG_UNSET, &
                                         set_cohort_size_from_carbon
   use meds_site_state_types, only : site_t, cohort_tissue_heat_capacity, cohort_wood_store_carbon, &
                                     TISSUE_C_LEAF, TISSUE_C_SAPW,                            &
                                     TISSUE_HCAP_MIN
   use meds_slow_ledger,          only : slow_ledger_t, slow_ledger_mark, slow_ledger_declare,     &
                                         slow_fast_carbon_handover, SLOW_PHASE_ALLOCATE,           &
                                         SLOW_PHASE_GROW, SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT,   &
                                         SLOW_PHASE_DISTURB, SLOW_PHASE_PATCH, KGC_PER_UMOL_C,   &
                                         slow_tissue_heat
   use meds_demography_update, only : update_cohort_states, fill_cohort_deriv, advance_pool,           &
                                      update_overtopping_lai
   use meds_demography_cohort_fusefiss, only : apply_recruitment, new_fuse_cohorts, terminate_cohorts, split_cohorts, sort_cohorts
   use meds_demography_patch_fusefiss, only : apply_patch_disturbance, new_fuse_patches, terminate_patches, sort_patches
   use meds_demography_rates,    only : npp_to_growth, camac_mortality, decline_to_floor_mortality,  &
                                        npp_to_recruitment
   use meds_plant_trait_dynamics, only : light_plastic_traits, update_plastic_trait
   use meds_phenology_types, only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   use meds_phenology, only : phenology_kernel, leaf_turnover_step, turnover_gate
   use meds_plant_carbon_allocation, only : plant_carbon_allocation
   use meds_sink_limitation,      only : growth_sink_limitation
   use meds_litter_partition, only : add_necromass
   use meds_biogeochem_types, only : litter_input_t
   use meds_site_diag_types,      only : CS_DDBH_DT, CS_DAGB_DT, CS_MORT_RATE, CS_NPP_LEAF,      &
                                        CS_NPP_FINEROOT, CS_NPP_WOOD, CS_NPP_STORAGE,           &
                                        CS_NPP_REPRO, CS_GROWTH_RESP, CS_STORAGE_RESP, CS_EXUDATE, &
                                        cohort_diag_grow,                                      &
                                        PD_LITTER_LEAF, PD_LITTER_FINEROOT, PD_LITTER_STRUCT,  &
                                        PD_LITTER_WOOD, PD_LITTER_REPRO,                       &
                                        PD_MORT_C_BACKGROUND,                                &
                                        PD_RECRUIT_NPLANT,                                     &
                                        cohort_diag_reset
   implicit none
   private

   public :: vegetation_dynamics, restructure_stand, advance_leaf_phenology, advance_plant_traits
   public :: reacclimate_plant_traits
   public :: shed_turnover_water, shed_excess_tissue_water, fill_new_leaf_water, accumulate_recruit_pool
   public :: update_cohort_derivatives

   !----- Guard: baseline leaf turnover = 1/llspan; floor llspan so a degenerate value cannot divide-by-0.
   real(wp), parameter :: tiny_llspan = 1.0e-6_wp         !< [yr]
   !----- Patch structural dynamics (fusion/fission/disturbance) are an ANNUAL process, so       !
   !       disturbance integrates its yearly rate over this interval, not the per-step dt.        !
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !< [yr]

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the carbon NPP, compute the carbon !
   ! vital rates via the plant kernels, and apply growth and mortality. The calendar's         !
   ! restructuring is not part of the step: restructure_stand runs it at the boundary.         !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, doy, ledger, latitude_deg)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in), optional :: doy   !< day-of-year at the step start (drives phenology)
      !----- The per-patch litter accumulator is site%patch%litter_in, NOT a local array, so that it
      !      follows its patches when restructure_stand reorders, blends or creates them. See the
      !      field's comment in meds_site_state_types.
      type(slow_ledger_t), intent(inout), optional :: ledger    !< site conservation ledger (meds_slow_ledger)
      real(wp),            intent(in),    optional :: latitude_deg !< the polygon's (default [site])
      real(wp), allocatable    :: mortality(:), recruitment(:,:), npp_repro(:)
      logical,  allocatable    :: exhausted(:)
      real(wp)                 :: mort_water
      real(wp)                 :: tissue_heat0, tissue_heat1, handover, plastic_heat0
      real(wp), allocatable    :: nplant_before(:), leaf_carbon_before(:)
      real(wp)                 :: seed_rain_carbon
      integer(ik)              :: ip
      type(carbon_flux_block)  :: npp
      integer(ik)              :: n_window

      !----- 0. Leaf phenology (UNCONDITIONAL): advance the per-cohort leaf tendencies ONE daily      !
      !         step, BEFORE compute_carbon_allocation reads them (ED2 too runs its phenology         !
      !         driver inside the vegetation-dynamics slow loop). Needs the step-start day-of-year;   !
      !         a caller with no calendar context (e.g. the Python carbon-mode C-API path) simply     !
      !         omits doy, so the tendencies stay at their birth values (flushing, not senescing).    !
      !         advance_leaf_phenology ALSO no-ops on its own when no fast sub-step has yet           !
      !         supplied a daily-mean temperature (site%pheno_tair_n < 1). ---------------------------!
      if (present(doy)) call advance_leaf_phenology(site, cfg, doy, latitude_deg)

      !----- 0b. Light trait plasticity (opt-in): acclimate the per-cohort leaf traits toward     !
      !          their shaded targets from last step's overtopping LAI, BEFORE compute_carbon_allocation reads  !
      !          cohort%sla / cohort%llspan. OFF => traits stay at top-of-canopy. -------------------!
      !          A rising SLA resorbs leaf carbon into storage, which has no heat capacity, so the  !
      !          tissue's heat content changes at unchanged temperature: declared to the ledger,  !
      !          as the cohort and patch phases declare theirs. ----------------------------------!
      if (cfg%trait_plasticity_on) then
         if (present(ledger)) plastic_heat0 = slow_tissue_heat(site)
         call advance_plant_traits(site, cfg, cfg%dt_years)
         if (present(ledger)) call slow_ledger_declare(ledger, energy_in = slow_tissue_heat(site) - plastic_heat0)
      end if

      !----- 1. Carbon NPP from the plant seam. Also accumulates this step's leaf/fine-root         !
      !         TURNOVER litter (leaf_shed_c/fineroot_shed_c) into the per-patch litter             !
      !         accumulator. ZERO it for this step first: it is patch state, so it arrives          !
      !         carrying the last step's values. ---------------------------------------------------!
      site%patch%litter_in(1:site%patch%n) = litter_input_t()
      !----- SLOW diagnostic block: size + zero it for this step BEFORE anything writes into it.   !
      !      SET semantics (one sample per slow step), so `w` is a presence flag rather than a dt.   !
      if (site%cohort%sdiag%active) then
         call cohort_diag_grow(site%cohort%sdiag, max(site%cohort%n, 1_ik))
         call cohort_diag_reset(site%cohort%sdiag)
         site%cohort%sdiag%n = site%cohort%n
      end if
      call compute_carbon_allocation(site, cfg, cfg%dt_years, npp, npp_repro, site%patch%litter_in, exhausted)

      !----- 1b. Leaf/fine-root turnover WATER shedding: must run BEFORE update_cohort_states       !
      !          commits this step's new leaf_carbon/fineroot_carbon below, since it needs the PRE-   !
      !          step pools as its "how much was lost" denominator (see shed_turnover_water). --------!
      call shed_turnover_water(site, cfg, npp, ledger)
      !----- LEDGER: the allocation phase COMPUTES; the only state it moves is tissue water into    !
      !      patch%shed_water_rate, a handoff to the fast tier that shed_turnover_water declares,   !
      !      so this phase should close on every currency. It is marked anyway -- a phase that      !
      !      must be quiet is the one worth watching, because a residual here means something       !
      !      moved that nobody declared.                                                            !
      if (present(ledger)) call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_ALLOCATE)

      !----- 2. Carbon vital RATES via the rate laws, PRE-apply: mortality reads the growth_avg    !
      !         of the steps before this one (this step's sample enters the ring buffer below).    !
      call compute_vital_rates(site, cfg, npp%wood, npp_repro, exhausted, cfg%dt_years, mortality, recruitment)

      !----- 2b0. Credit the recruit pool with THIS STEP's recruitment, every step: the exact      !
      !           integral of a daily quantity, where a monthly credit from one day's rate scaled  !
      !           to the month would be a 12-point sample that depends on which days are month     !
      !           boundaries.                                                                      !
      !                                                                                            !
      !           Crediting daily also makes the pool an EXACT carbon quantity, with no new state. !
      !           `recruitment` is n*(npp_repro/dt_yr)*efficiency/carbon_min, so rate*dt_yr is     !
      !           n*npp_repro*efficiency/carbon_min -- the establishing share of the reproduction  !
      !           carbon the parents were debited this very step, in plant units. The pool valued  !
      !           at carbon_min is therefore that carbon exactly, and the growth phase's debit and !
      !           the pool's credit are the same number in the same phase. ------------------------!
      call accumulate_recruit_pool(site, cfg, recruitment, cfg%dt_years, seed_rain_carbon)
      if (present(ledger)) call slow_ledger_declare(ledger, carbon_in = seed_rain_carbon)

      !----- 2b. Mortality litter and mortality water are computed AFTER the commit below, on the         !
      !          density it actually removed (see accumulate_mortality_litter). The litter goes           !
      !          into patch%litter_in, which the daily soil_carbon_step (meds_biogeochem_dynamics)        !
      !          consumes as the matrix ODE's litter input, so adding it straight to the soil             !
      !          pools here would count it twice. --------------------------------------------------------!
      !----- 3. Sequence the apply-primitives. -----------------------------------------------------------!
      n_window = growth_window_steps(cfg)
      site%growth_hist_pos = mod(site%growth_hist_pos, n_window) + 1_ik

      !----- Growth + mortality via the TENDENCY seam: compute the per-cohort time-derivatives     !
      !      (carbon flip -> geometry/pool rates + the realized-dbh-rate moving-average ring       !
      !      buffer; mortality -> log-space rate) HERE in the driver, into the site-carried        !
      !      scratch bundle, then let the demography engine's pure applier advance the state.      !
      !      The ring buffer is refreshed inside update_cohort_derivatives, so mortality's         !
      !      growth_avg dependency holds.                                                          !
      call cohort_deriv_alloc(site%deriv, site%cohort%n)
      call update_cohort_derivatives(site, npp, mortality, cfg%dt_years, n_window, site%growth_hist_pos)
      !----- TISSUE THERMAL MASS across the commit. A cohort's heat capacity is a function of its  !
      !      biomass and density, so growing or dying changes `cap * T` with no flux at all: the   !
      !      model makes thermal mass out of carbon, and unmakes it, with no counterparty in       !
      !      either direction.                                                                     !
      !                                                                                            !
      !      The whole change is declared as ONE exchange because it is one mechanism, and because !
      !      splitting it into a growth part and a mortality part needs a cross-term convention    !
      !      the physics does not supply. `update_cohort_states`, the water new leaves take from    !
      !      the wood (fill_new_leaf_water) and the mortality water shed below ALL leave the        !
      !      temperatures alone, so this is purely the thermal-mass term.                          !
      !                                                                                            !
      !      WHY IT IS AN EXCHANGE AND NOT A LEAK: new tissue is assembled from CO2 and water at    !
      !      the plant's own temperature, and it arrives carrying the sensible heat of that mass.  !
      !      The model tracks no thermal content for CO2 (the canopy air's capacity is dry air     !
      !      alone), so that heat genuinely crosses the boundary of the modelled thermal system.   !
      !      A fuller treatment would give CO2 a heat capacity in the CAS and carry enthalpy on    !
      !      root water uptake; both are far larger than this, and the approximation is named here !
      !      rather than buried.  ------------------------------------------------------------!
      tissue_heat0 = slow_tissue_heat(site)
      !----- The fast tier's handover is valued at the density that FIXED it: gpp_accum is per     !
      !      plant and was accumulated over the day at the pre-commit density, and the store's     !
      !      growth term is n0.(p1-p0) by the decomposition in accumulate_mortality_litter. Taken  !
      !      after the commit it would be n1, short by n_died x net -- the same offset, in the      !
      !      other direction.  --------------------------------------------------------------!
      handover = slow_fast_carbon_handover(site, cfg)
      if (allocated(nplant_before)) deallocate(nplant_before)
      allocate(nplant_before(max(site%cohort%n, 1_ik)))
      nplant_before(1:site%cohort%n) = site%cohort%nplant(1:site%cohort%n)
      leaf_carbon_before = site%cohort%leaf_carbon(1:site%cohort%n)
      call update_cohort_states(site%cohort, site%deriv, cfg%dt_years, cfg%negligible_nplant)
      call fill_new_leaf_water(site, cfg, leaf_carbon_before)
      tissue_heat1 = slow_tissue_heat(site)

      !----- Mortality's litter and water, on the density the applier actually removed and the     !
      !      pools it left behind (the decomposition is in accumulate_mortality_litter's header).  !
      !      UNCONDITIONAL: the routine also records how much carbon died, a demographic quantity  !
      !      whether or not the soil pools are on, and the litter is recorded either way.          !
      call accumulate_mortality_litter(site, cfg, nplant_before, site%patch%litter_in)
      call shed_mortality_water(site, cfg, nplant_before, mort_water)

      !----- Re-sort every step: growth changed heights, so re-establish the tallest-first order  !
      !      the overtopping-LAI sweep + the patch-light profiles depend on (fuse/fiss also sorts, !
      !      but only on the monthly/annual cadence). -------------------------------------------!
      call sort_cohorts(site)

      !----- LEDGER: the growth phase is where every carbon movement is COMMITTED, and it needs two !
      !      declared terms because both cross a boundary this ledger's stores cannot see.          !
      !                                                                                             !
      !        carbon_in  -- the fast loop already took GPP out of the canopy air and paid the      !
      !                      maintenance respiration back into it; the net sat in gpp_accum /       !
      !                      *_resp_accum, which are accumulators, not stores. Without this term    !
      !                      the day's photosynthate would read as carbon appearing from nowhere.   !
      !        carbon_out -- necromass routed into the CENTURY pools. It leaves the live pools here !
      !                      and arrives in the soil store at SLOW_PHASE_SOILC, so declaring it at  !
      !                      both ends TESTS the litter seam rather than hiding it. Plus the CO2    !
      !                      owed to the canopy air (slow_co2_handoff), which reaches it in the     !
      !                      next fast window.                                                      !
      !                                                                                             !
      !      With these declared the phase closes to round-off: a residual here is a movement that  !
      !      nobody declared.                                                                       !
      if (present(ledger)) then
         call slow_ledger_declare(ledger, carbon_in  = handover,                                    &
                                          carbon_out = litter_carbon_total(site, site%patch%litter_in)              &
                                                     + slow_co2_handoff(site, cfg),                &
                                          water_out  = mort_water,                                 &
                                          energy_in  = tissue_heat1 - tissue_heat0)
         call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_GROW)
      end if

      !----- Patch ageing, the slow per-patch state vegetation and biogeochemistry share, runs in     !
      !      meds_slow_dynamics after this routine. --------------------------------------------------!

      !----- SLOW patch diagnostics: litterfall and recruitment. Recorded HERE, before the       !
      !      calendar's restructuring (restructure_stand) permutes the patch axis.                 !
      !                                                                                            !
      !      UNITS. The block's weight `w` is Sum(dt) in SECONDS -- accumulate_patch_diag            !
      !      adds cfg%dt_fast every fast step -- and the reader returns v/w. So EVERY row, fast or   !
      !      slow, must contribute (its rate, in the units the registry declares) x (dt in SECONDS). !
      !      The fast rows do that with flux*dt_fast. These slow rows must do it with dt_slow, and   !
      !      their declared rates are per YEAR, so the per-year rate is formed first and weighted    !
      !      second: patch%litter_in holds this step's litter AMOUNT [kgC/m2], whose per-year rate   !
      !      is amount/dt_years; `recruitment` is already [plant/m2/yr] and needs only the weight.   !
      !                                                                                            !
      !      patch%litter_in is the SAME accumulator biogeochemistry consumes, so litter_*_site and  !
      !      the soil-carbon input are one number by construction, not two derivations that agree. -!
      if (site%patch%diag%active) then
         !----- Without a fast loop the step also sets the block's patch count, which the fast loop  !
         !      otherwise sets, so a reader sees these rows. -------------------------------------------!
         if (.not. cfg%fast_biophysics_on) site%patch%diag%n = site%patch%n
         do ip = 1_ik, site%patch%n
            !----- Leaf and fine-root litter by the tissue they came from: what falls, not the        !
            !      share of it that reaches one destination pool. ----------------------------------!
            site%patch%diag%v(PD_LITTER_LEAF,     ip) = site%patch%diag%v(PD_LITTER_LEAF,     ip) &
                                                      + (site%patch%litter_in(ip)%leaf_c                &
                                                         / cfg%dt_years) * cfg%dt_slow
            site%patch%diag%v(PD_LITTER_FINEROOT, ip) = site%patch%diag%v(PD_LITTER_FINEROOT, ip) &
                                                      + (site%patch%litter_in(ip)%fineroot_c            &
                                                         / cfg%dt_years) * cfg%dt_slow
            site%patch%diag%v(PD_LITTER_WOOD,     ip) = site%patch%diag%v(PD_LITTER_WOOD,     ip) &
                                                      + ((site%patch%litter_in(ip)%wood_c               &
                                                        + site%patch%litter_in(ip)%storage_c)           &
                                                         / cfg%dt_years) * cfg%dt_slow
            site%patch%diag%v(PD_LITTER_REPRO,    ip) = site%patch%diag%v(PD_LITTER_REPRO,    ip) &
                                                      + (site%patch%litter_in(ip)%repro_c               &
                                                         / cfg%dt_years) * cfg%dt_slow
            site%patch%diag%v(PD_LITTER_STRUCT,   ip) = site%patch%diag%v(PD_LITTER_STRUCT,   ip) &
                                                      + ((site%patch%litter_in(ip)%struct_grnd          &
                                                        + site%patch%litter_in(ip)%struct_soil)         &
                                                         / cfg%dt_years) * cfg%dt_slow
            site%patch%diag%v(PD_RECRUIT_NPLANT,  ip) = site%patch%diag%v(PD_RECRUIT_NPLANT,  ip) &
                                                      + sum(recruitment(:, ip)) * cfg%dt_slow
            !----- The block's weight: the fast loop adds dt_fast each sub-step; without one, the step  !
            !      adds its own dt_slow, or a slow-only run reads every patch row back as 0. ----------!
            if (.not. cfg%fast_biophysics_on) site%patch%diag%w(ip) = site%patch%diag%w(ip) + cfg%dt_slow
         end do
      end if

      !----- 4. Refresh the overtopping-LAI competition diagnostic (the per-step sort above left   !
      !         the stand sorted). ----------------------------------------------------------------!
      call update_overtopping_lai(site)
   end subroutine vegetation_dynamics

   !---------------------------------------------------------------------------------------!
   ! restructure_stand -- the calendar's restructuring of the cohorts and patches: monthly    !
   ! recruitment, cohort fusion, culling, splitting and sort; yearly patch disturbance and    !
   ! patch fusion. It runs AT the boundary, between the last step of the period that ends and !
   ! the first step of the one that begins, after the output has read that last step: the     !
   ! period's records then describe one stand throughout, and what the restructuring does --  !
   ! its events and the stand it leaves -- belongs to the new period (docs/science/           !
   ! diagnostics.md §4). `is_new_month` / `is_new_year` say which boundary this is.           !
   !---------------------------------------------------------------------------------------!
   subroutine restructure_stand(site, cfg, is_new_month, is_new_year, ledger)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      type(slow_ledger_t), intent(inout), optional :: ledger    !< site conservation ledger (meds_slow_ledger)
      real(wp) :: cull_water, th0, th1, rec_carbon, rec_heat, rec_water, dist_water
      logical  :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      !----- The calendar cadence + demography on/off + fiss/fuse switches. -----------------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- Cohort restructuring (monthly): recruit + fuse/split + sort. -------------------!
      if (do_cohort_fissfuse) then
         call apply_recruitment(site, cfg, rec_carbon, rec_heat, rec_water)
         if (present(ledger)) call slow_ledger_declare(ledger, carbon_in = rec_carbon,             &
                                                               energy_in = rec_heat,               &
                                                               water_in  = rec_water)
         !----- Every structural operator below changes WHICH biomass exists at what temperature, !
         !      so each changes the tissue's thermal mass -- the same exchange the growth commit   !
         !      declares, for the same reason (see there). Fusion and fission re-derive the        !
         !      sapwood ring from the merged or perturbed diameter, and `sapwood_fraction` is      !
         !      NONLINEAR in dbh, so the merged capacity is not the sum of the two even though     !
         !      every carbon pool is. The cull removes its cohort's capacity outright.             !
         !                                                                                          !
         !      Cohort fusion mixes the temperatures on heat capacity, so what is declared here is !
         !      the thermal-mass term alone, not a wrong average. ---------------------------------!
         if (present(ledger)) th0 = slow_tissue_heat(site)
         !----- Marked ON ITS OWN: recruitment is the one structural operator that CREATES matter   !
         !      rather than rearranging it, and the carbon and heat it draws from outside the model !
         !      (apply_recruitment's carbon_drawn, heat_drawn) are only checkable if nothing else   !
         !      shares the phase. Fusion and fission, which must be exact, are marked together      !
         !      below.                                                                              !
         if (present(ledger)) call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_RECRUIT)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg, cull_water)
         call split_cohorts(site, cfg)
         call sort_cohorts(site)
         if (present(ledger)) then
            th1 = slow_tissue_heat(site)
            call slow_ledger_declare(ledger, water_out = cull_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_COHORT)
         end if
      end if

      !----- Patch disturbance (annual) then patch restructuring (annual, independent). ----!
      if (do_patch_disturbance) then
         if (present(ledger)) th0 = slow_tissue_heat(site)
         call apply_patch_disturbance(site, cfg, PATCH_DYNAMICS_INTERVAL, dist_water)
         if (present(ledger)) then
            th1 = slow_tissue_heat(site)
            call slow_ledger_declare(ledger, water_out = dist_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_DISTURB)
         end if
      end if
      if (do_patch_fissfuse) then
         if (present(ledger)) th0 = slow_tissue_heat(site)
         call sort_patches(site)
         call new_fuse_patches(site, cfg)
         call terminate_patches(site, cfg)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg, cull_water)
         call sort_cohorts(site)
         if (present(ledger)) then
            th1 = slow_tissue_heat(site)
            call slow_ledger_declare(ledger, water_out = cull_water, energy_in = th1 - th0)
            call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_PATCH)
         end if
      end if

      !----- Refresh the overtopping-LAI competition diagnostic for the restructured stand. -----!
      call update_overtopping_lai(site)
   end subroutine restructure_stand

   !----- The CO2 this step handed to the fast tier [kgC/m2 site]: growth respiration and           !
   !       storage maintenance, less the starvation over-report -- the mirror of                    !
   !       slow_fast_carbon_handover. It leaves the plant inside this ledger's window and reaches   !
   !       the canopy air inside the NEXT fast window, so like the shed water it is a declared      !
   !       handoff rather than a store. ------------------------------------------------------------!
   pure function slow_co2_handoff(site, cfg) result(c)
      type(site_t),        intent(in) :: site
      type(meds_config_t), intent(in) :: cfg
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, site%patch%n
         c = c + site%patch%area(ip) * site%patch%slow_co2_rate(ip) * KGC_PER_UMOL_C * cfg%dt_slow
      end do
   end function slow_co2_handoff

   !----- Site-total carbon in this step's litter accumulator [kgC/m2 site]. The ledger declares    !
   !       this leaving the live pools at SLOW_PHASE_GROW and arriving in the CENTURY pools at      !
   !       SLOW_PHASE_SOILC. With soil_carbon_on off it is still accumulated and declared leaving   !
   !       at SLOW_PHASE_GROW, and nothing arrives at SLOW_PHASE_SOILC: no soil pool takes it, so   !
   !       the necromass leaves the model as a declared outflow. -----------------------------------!
   pure function litter_carbon_total(site, lit) result(c)
      type(site_t),         intent(in) :: site
      type(litter_input_t), intent(in) :: lit(:)
      real(wp)    :: c
      integer(ik) :: ip
      c = 0.0_wp
      do ip = 1_ik, site%patch%n
         c = c + site%patch%area(ip) * (lit(ip)%labile_grnd + lit(ip)%labile_soil                   &
                                      + lit(ip)%struct_grnd + lit(ip)%struct_soil)
      end do
   end function litter_carbon_total

   !---------------------------------------------------------------------------------------!
   ! CARBON-mode TENDENCY computer, kept apart from the applier so the ENGINE only applies. Per     !
   ! cohort it adds this step's NPP to the four pools (tentatively, in locals -- state is NOT       !
   ! committed here), backs out each pool's time-derivative d_X_dt = (X_new - X_old)/dt for the     !
   ! applier, and runs the allometric FLIP wood_carbon->dbh (carbon_to_structure) only for          !
   ! the reported rates: the applier derives the committed geometry from the committed carbon       !
   ! itself. It ALSO advances the moving-average ring buffer with the realized dbh-increment RATE   !
   ! (persistent state) so Camac mortality sees carbon growth, and records the log-space            !
   ! mortality rate dln_nplant_dt = -mortality. Host-only (the flip calls the branchy wood_to_dbh). !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_derivatives(site, npp, mortality, dt_yr, n_window, hist_pos)
      type(site_t),             intent(inout) :: site       ! inout: fills site%deriv + refreshes the ring buffer
      type(carbon_flux_block),  intent(in)    :: npp
      real(wp),                 intent(in)    :: mortality(:)
      real(wp),                 intent(in)    :: dt_yr
      integer(ik),              intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: lc_new, fc_new, wc_new, nc_new
      real(wp)    :: dbh_new, height_new, ba_new, agb_new, la_new, dbh_rate

      associate (cohort => site%cohort)
         do i = 1_ik, cohort%n
            !----- Tentative new pools, in locals, by the applier's own rule (advance_pool): a pool !
            !      the step empties lands on 0 to round-off, and an overdraw stays negative, so the  !
            !      applier carries it on and the bounds guard stops the run naming the pool. A floor !
            !      here would create the carbon the overdraw took before the applier could see it. -!
            lc_new = advance_pool(cohort%leaf_carbon(i),          npp%leaf(i))
            fc_new = advance_pool(cohort%fineroot_carbon(i),      npp%fineroot(i))
            wc_new = advance_pool(cohort%wood_carbon(i),          npp%wood(i))
            nc_new = advance_pool(cohort%nonstructural_carbon(i), npp%nonstructural(i))
            !----- FLIP the geometry from the tentative pools (dbh from wood_carbon, leaf_area    !
            !      from leaf_carbon) into locals, for the growth rate the ring buffer and the     !
            !      diagnostics report; update_cohort_states commits the same flip. ---------------!
            call carbon_to_structure(cohort%allom, wc_new, lc_new, cohort%p_wood_density(i),        &
                                     cohort%p_hgt_max(i), cohort%p_aboveground_frac(i), cohort%sla(i), &
                                     dbh_new, height_new, ba_new, agb_new, la_new)
            !----- Back out the pools' tendencies + advance the ring buffer with the REALIZED dbh !
            !      rate (fill_cohort_deriv, which the C-API's empirical path shares). -------------!
            dbh_rate = (dbh_new - cohort%dbh(i)) / dt_yr
            call fill_cohort_deriv(cohort, i, site%deriv, dt_yr, mortality(i), dbh_rate,            &
                                   lc_new, fc_new, wc_new, nc_new, n_window, hist_pos)
            !----- Mirror the two REPORTED tendencies into the lockstep-riding slow diagnostic       !
            !      block. site%deriv itself is transient and deliberately NOT reordered, so reading   !
            !      it at the output tick -- which runs AFTER the monthly fiss/fuse permutes the       !
            !      cohort axis -- would pair tendency i with a different plant i.  ------------------!
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CS_DDBH_DT,   i) = dbh_rate
               cohort%sdiag%v(CS_DAGB_DT,   i) = (agb_new - cohort%agb(i)) / dt_yr
               cohort%sdiag%v(CS_MORT_RATE, i) = mortality(i)
               cohort%sdiag%w(i)               = 1.0_wp
            end if
         end do
      end associate
   end subroutine update_cohort_derivatives

   !---------------------------------------------------------------------------------------!
   ! CARBON vital-RATE assembler: turn the per-cohort carbon fluxes (npp_wood for the dbh rate,   !
   ! npp_repro for recruits) into the mortality [1/yr] + recruitment [plant/m2/yr] arrays, by     !
   ! calling the per-individual rate laws (meds_demography_rates). The site-level sweep (patch    !
   ! reduction + baseline seed rain + the GROWTH_AVG_UNSET test) stays here in the DRIVER, so     !
   ! the laws remain state-free.                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine compute_vital_rates(site, cfg, npp_wood, npp_repro, exhausted, dt_yr, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp),              intent(in)  :: npp_wood(:), npp_repro(:), dt_yr
      logical,               intent(in)  :: exhausted(:)  !< no leaves and no storage (compute_carbon_allocation)
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr] per cohort
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, j, pf, ip
      real(wp)    :: dbh_rate
      real(wp), allocatable :: carbon_min(:), efficiency(:)

      n    = site%cohort%n
      np   = site%patch%n
      npft = cfg%pft%n
      allocate(mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp
      efficiency  = cohort_recruit_efficiency(site, cfg)

      associate (cohort => site%cohort, pft => cfg%pft)
         allocate(carbon_min(npft))
         do pf = 1_ik, npft
            carbon_min(pf) = min_cohort_carbon(cohort%allom, pft%min_cohort_height, pft%wood_density(pf))
         end do
         !----- Baseline external seed rain (independent of structure). --------------------!
         do ip = 1_ik, np
            do pf = 1_ik, npft
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%seed_rain_recruits(pf)
            end do
         end do
         do j = 1_ik, n
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !----- Prospective carbon dbh-increment rate (the mortality instantaneous term). --!
            dbh_rate = npp_to_growth(cohort%allom, cohort%wood_carbon(j), npp_wood(j), cohort%dbh(j),  &
                                     dt_yr, cohort%p_wood_density(j), cohort%p_hgt_max(j),        &
                                     cohort%p_aboveground_frac(j))
            mortality(j) = camac_mortality(cohort%growth_avg(j), dbh_rate,                          &
                                           cohort%growth_avg(j) /= GROWTH_AVG_UNSET,                &
                                           pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))
            !----- A cohort with no leaves and no storage cannot flush again, so it dies: its density  !
            !      goes to the floor this step, down the path every death takes (its carbon to        !
            !      litter, its water to the ground), and the next monthly cull removes it. This is a   !
            !      special case, not a starvation hazard: a cohort that still has leaves or storage    !
            !      dies at the rate its growth sets (docs/science/plant_carbon_allocation.md §2). ----!
            if (exhausted(j)) mortality(j) = max(mortality(j),                                      &
                    decline_to_floor_mortality(cohort%nplant(j), cfg%negligible_nplant, dt_yr))
            !----- Reproduction carbon -> recruits (an annual rate; accumulate_recruit_pool       !
            !      credits rate x dt each step), gated by include_pft (a disabled PFT must not    !
            !      re-establish). ----------------------------------------------------------------!
            if (pft%include_pft(pf) == 1_ik)                                                       &
            recruitment(pf, ip) = recruitment(pf, ip)                                              &
                 + npp_to_recruitment(cohort%nplant(j), npp_repro(j), dt_yr,                       &
                                      efficiency(j), carbon_min(pf))
         end do
      end associate
   end subroutine compute_vital_rates

   !---------------------------------------------------------------------------------------!
   ! Each cohort's recruit carbon efficiency this step: the share of its reproduction carbon   !
   ! that becomes recruits in its own patch. It is its PFT's recruit_carbon_efficiency, reduced !
   ! by exp(-recruit_shade_decay * LAI) of that patch: under a deep canopy fewer seeds germinate !
   ! and fewer seedlings live to min_cohort_height. The recruits (compute_vital_rates) and the   !
   ! seed lost to litter (compute_carbon_allocation) both take it from here, so the two always  !
   ! add up to the reproduction carbon.                                                        !
   !---------------------------------------------------------------------------------------!
   function cohort_recruit_efficiency(site, cfg) result(efficiency)
      type(site_t),        intent(in) :: site
      type(meds_config_t), intent(in) :: cfg
      real(wp), allocatable :: efficiency(:)
      real(wp), allocatable :: patch_lai(:)
      integer(ik) :: j, pf, ip
      allocate(efficiency(site%cohort%n), patch_lai(max(site%patch%n, 1_ik)))
      patch_lai = 0.0_wp
      associate (cohort => site%cohort, pft => cfg%pft)
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
   ! CARBON NPP assembler + slow-carbon ORCHESTRATOR. Per cohort it (1) gathers the daily GPP /  !
   ! maintenance respiration (fast-loop accumulators, or the gpp_ref stub); (2) applies TURNOVER  !
   ! first (cohort_carbon_demand -> the leaf / fine-root loss AMOUNTS = this step's litter), so    !
   ! the growth demands are formed against the POST-SHED pool; (3) calls the ELEMENTAL,            !
   ! growth-only plant_carbon_allocation kernel; and (4) forms the NET per-pool change (growth -   !
   ! shed). Turnover-first makes the within-step leaf replacement exact (a canopy that is not      !
   ! senescing holds target). docs/science/plant_carbon_allocation.md §1-§2, §5.                   !
   !---------------------------------------------------------------------------------------!
   subroutine compute_carbon_allocation(site, cfg, dt_yr, npp, npp_repro, lit, exhausted)
      !----- intent(inout) ONLY to record the NPP allocation split into the slow diagnostic block; !
      !      no prognostic state is touched here (the appliers still own that).  ------------------!
      type(site_t),            intent(inout) :: site
      type(meds_config_t),     intent(in)  :: cfg
      real(wp),                intent(in)  :: dt_yr
      type(carbon_flux_block), intent(out) :: npp
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

      n = site%cohort%n
      dt_day = cfg%dt_slow / day_sec            ! days in the slow step (phenology rates are [1/day])
      allocate(npp%leaf(n), npp%fineroot(n), npp%wood(n), npp%nonstructural(n), npp_repro(n), exhausted(n))
      allocate(co2_owed(site%patch%n)) ; co2_owed = 0.0_wp
      efficiency = cohort_recruit_efficiency(site, cfg)
      !----- The growth ring-buffer slot this step's realized growth will be written to (step 3 of  !
      !      vegetation_dynamics advances growth_hist_pos the same way): the sink limit's window is  !
      !      the buffer without that slot's old sample, plus this step. ----------------------------!
      n_window = growth_window_steps(cfg)
      next_pos = mod(site%growth_hist_pos, n_window) + 1_ik
      associate (cohort => site%cohort, pft => cfg%pft)
         do j = 1_ik, n
            pf  = cohort%pft(j)
            r2l = pft%root_to_leaf_ratio(pf)
            !----- Leaf-AREA-conserving target: uses the cohort's (plastic) SLA, so a shaded canopy  !
            !      reaches the same allometric leaf area with less leaf carbon. --------------------!
            leaf_target  = size2leaf_carbon(cohort%allom, cohort%dbh(j), cohort%height(j), cohort%sla(j))
            leaf_turn    = 1.0_wp / max(cohort%llspan(j), tiny_llspan)   ! baseline leaf turnover [1/yr]
            !----- STORAGE MAINTENANCE (docs/science/plant_carbon_allocation.md eq. (2)), charged        !
            !      BEFORE the demand and the allocation, which is ED2's ordering (apply_maintenance      !
            !      runs ahead of growth). The form is ED2's growth_balive.f90: a fractional turnover     !
            !      of the pool per year, with no temperature dependence (the page says why).             !
            !                                                                                          !
            !      The fraction is clamped to [0,1] rather than the AMOUNT to the pool, so a large       !
            !      rate x dt cannot make the charge exceed what is there -- and cannot go negative       !
            !      and manufacture storage.                                                              !
            !                                                                                          !
            !      IT IS A TENDENCY, NOT A WRITE. The charge is netted into npp%nonstructural and        !
            !      applied with every other pool change in the GROW phase; this routine does not touch   !
            !      cohort%nonstructural_carbon. That is the standing rule -- the driver computes, the    !
            !      engine applies -- and here it is also what keeps the ledger closed, because the       !
            !      phase that DECLARES the CO2 efflux (slow_co2_handoff, GROW) has to be the phase       !
            !      where the pool drops.                                                                 !
            !                                                                                          !
            !      `store_post` is what the ALLOCATOR sees, so maintenance is still paid before growth   !
            !      is funded from the reserve (ED2's apply_maintenance-then-allocate ordering) even      !
            !      though the pool itself moves one phase later.                                         !
            storage_maint = cohort%nonstructural_carbon(j)                                              &
                          * min(1.0_wp, max(0.0_wp, pft%storage_turnover_rate(pf) * dt_yr))
            store_post    = cohort%nonstructural_carbon(j) - storage_maint
            !----- Rides co2_owed, the SAME per-patch channel growth respiration uses to reach the      !
            !      fast loop's nee_biotic, and is declared to the ledger by slow_co2_handoff with it.    !
            co2_owed(cohort%owner_patch(j)) = co2_owed(cohort%owner_patch(j))                           &
                                            + cohort%nplant(j) * storage_maint

            store_demand = max(0.0_wp, storage_carbon_target(pft%storage_cushion(pf), leaf_target,      &
                                                             r2l * leaf_target) - store_post)
            !----- Reproduction starts at a share of the local maximum height, so where the local   !
            !      allometry makes trees short they mature short (Thomas 1996's relative size at     !
            !      onset of maturity). --------------------------------------------------------------!
            repro_frac   = merge(pft%reproduction_investment_fraction(pf), 0.0_wp,                  &
                                 cohort%height(j) >= pft%reproduction_height_fraction                &
                                                     * local_hmax(cohort%allom, cohort%p_hgt_max(j)))
            !----- Gather this step's GPP/respiration, apply turnover-first, and form the flush-    !
            !      capped growth demands (the single per-cohort "how much carbon, how much litter"  !
            !      computation, isolated from the surrounding orchestration). ----------------------!
            !----- Fine roots turn over at the PFT's rate, or, with fineroot_lifespan_plastic, as  !
            !      much slower as the cohort's leaves live longer than at the top of the canopy. ---!
            root_turn = pft%fineroot_turnover_rate(pf)
            if (cfg%fineroot_lifespan_plastic .and. cohort%llspan(j) > 0.0_wp)                     &
               root_turn = root_turn * pft%leaf_lifespan_toc(pf) / cohort%llspan(j)
            call cohort_carbon_demand(cfg%fast_biophysics_on,                                      &
                     cohort%gpp_accum(j), cohort%leaf_resp_accum(j), cohort%stem_resp_accum(j),    &
                     cohort%root_resp_accum(j), cfg%gpp_ref, cohort%leaf_area(j), dt_yr,           &
                     cohort%leaf_flush_tendency(j), cohort%leaf_shed_tendency(j),                  &
                     pft%pheno_flush_rate_max(pf), pft%pheno_shed_rate_max(pf), leaf_turn,         &
                     root_turn,                                                                    &
                     pft%pheno_min_leaf_cover(pf), pft%pheno_bare_leaf_cover(pf),                  &
                     cohort%leaf_carbon(j), cohort%fineroot_carbon(j), leaf_target, dt_day, r2l,   &
                     gross_gpp, resp_maint, leaf_shed_c, fineroot_shed_c, leaf_demand,             &
                     fineroot_demand, pft%retained_carbon_fraction(pf), leaf_retained_c,           &
                     fineroot_retained_c)
            call growth_window_so_far(cohort%growth_count(j), cohort%growth_accum(j),              &
                                      cohort%growth_hist(next_pos, j), n_window, dt_yr,             &
                                      window_yr, grown_in_window)
            !----- Allocate the daily carbon to GROWTH (growth respiration charged on realized     !
            !      growth INSIDE the kernel; storage funds leaf/root growth even when net < 0). ---!
            call plant_carbon_allocation(gpp=gross_gpp, resp_maint=resp_maint,                     &
                     growth_resp_frac=pft%growth_resp_factor(pf), storage=store_post,                    &
                     leaf_demand=leaf_demand, fineroot_demand=fineroot_demand,                     &
                     storage_demand=store_demand, repro_frac=repro_frac,                           &
                     growth_leaf=g_leaf, growth_fineroot=g_fineroot, growth_wood=g_wood,          &
                     npp_store=npp_store, growth_repro=g_repro, growth_resp=growth_resp,          &
                     deficit=deficit, starving=starving,                                           &
                     wood_growth_max=growth_sink_limitation(                                       &
                        allom=cohort%allom,                                                          &
                        wood_carbon=cohort%wood_carbon(j), dbh=cohort%dbh(j),                       &
                        wood_density=cohort%p_wood_density(j), hgt_max=cohort%p_hgt_max(j),         &
                        aboveground_frac=cohort%p_aboveground_frac(j),                              &
                        max_relative_growth_rate=pft%max_relative_growth_rate(pf),                  &
                        max_absolute_growth_rate=pft%max_absolute_growth_rate(pf),                  &
                        max_absolute_growth_exponent=pft%max_absolute_growth_exponent(pf),          &
                        dt_yr=dt_yr, window_yr=window_yr, grown_in_window=grown_in_window),         &
                     exudate=exudate)
            !----- NET per-pool change = growth - shed (leaf/root); litter = shed. ---------------!
            npp%leaf(j)          = g_leaf     - leaf_shed_c
            npp%fineroot(j)      = g_fineroot - fineroot_shed_c
            npp%wood(j)          = g_wood
            !----- The storage maintenance leaves the pool HERE, as part of this step's tendency,     !
            !      and the resorbed leaf and fine-root carbon ARRIVES here for the same reason.        !
            !      npp%leaf and npp%fineroot above already carry the FULL shed removals, so this       !
            !      credit moves carbon between pools rather than creating any                          !
            !      (docs/science/plant_carbon_allocation.md eq. (3)).                                  !
            npp%nonstructural(j) = npp_store - storage_maint + leaf_retained_c + fineroot_retained_c
            npp_repro(j)         = g_repro
            !----- NO LEAVES AND NO STORAGE: the maintenance this step's photosynthesis and the whole  !
            !      storage pool left unpaid (`deficit`, > 0 only once storage is spent) is more than   !
            !      all the carbon in the cohort's leaves. Burning every leaf could not pay one step's  !
            !      upkeep, and nothing is left to pay for a flush, so the cohort cannot recover. The  !
            !      scale is the cohort's own maintenance respiration, so the test needs no threshold; !
            !      a cohort with a canopy holds far more leaf carbon than one step's upkeep, and one   !
            !      with storage pays it (docs/science/plant_carbon_allocation.md §2). ----------------!
            exhausted(j)         = deficit > cohort%leaf_carbon(j)
            !----- SLOW-loop diagnostics: the NPP allocation split, growth respiration, storage       !
            !      maintenance and exudate, so the carbon budget closes from the output file. As      !
            !      RATES [kgC/plant/yr], so they are comparable across slow-step lengths. Written     !
            !      into a block that rides the cohort lockstep, NOT into a transient scratch (see     !
            !      meds_site_diag_types). ------------------------------------------------------------!
            if (cohort%sdiag%active) then
               cohort%sdiag%v(CS_NPP_LEAF,     j) = g_leaf     / dt_yr
               cohort%sdiag%v(CS_NPP_FINEROOT, j) = g_fineroot / dt_yr
               cohort%sdiag%v(CS_NPP_WOOD,     j) = g_wood     / dt_yr
               cohort%sdiag%v(CS_NPP_STORAGE,  j) = npp_store  / dt_yr
               cohort%sdiag%v(CS_NPP_REPRO,    j) = g_repro    / dt_yr
               cohort%sdiag%v(CS_GROWTH_RESP,  j) = growth_resp / dt_yr
               cohort%sdiag%v(CS_STORAGE_RESP, j) = storage_maint / dt_yr
               cohort%sdiag%v(CS_EXUDATE,      j) = exudate / dt_yr
               cohort%sdiag%w(j)                  = 1.0_wp
            end if
            !----- EVERY ALLOCATOR OUTPUT HAS A DESTINATION (docs/science/plant_respiration.md §6).    !
            !                                                                                          !
            !      GROWTH RESPIRATION is autotrophic CO2, roughly 23% of the carbon entering growth.   !
            !      It goes to the canopy air, as a frozen daily rate the fast loop adds to             !
            !      nee_biotic -- the carbon twin of the shed-water handoff.                            !
            !                                                                                          !
            !      THE STARVATION DEFICIT rides the same channel with the opposite sign. The fast loop !
            !      has ALREADY exhaled the full maintenance respiration; when storage could not fund   !
            !      it, `deficit` is the part no pool paid for. Real maintenance respiration is         !
            !      substrate-limited, so the honest reading is that the fast loop OVER-REPORTED, and   !
            !      an operator-split model corrects an over-report on the next step rather than        !
            !      rewriting the last one. Netting it here keeps one channel and one sign convention.  !
            ip = cohort%owner_patch(j)
            co2_owed(ip) = co2_owed(ip) + cohort%nplant(j) * (growth_resp - deficit)

            !----- leaf_shed_c + fineroot_shed_c = this step's TURNOVER litter -> the per-patch      !
            !      litter accumulator, whether or not the soil-carbon pools are on: with them         !
            !      off the necromass leaves the model, declared at SLOW_PHASE_GROW, and the litter     !
            !      outputs still report it. The SEED LOSS joins them: reproduction carbon is debited  !
            !      in full from the parent, and only the recruit carbon efficiency of it ever         !
            !      establishes -- the remaining ~99.9% is dead seed and dead seedling, which is        !
            !      NECROMASS, not nothing. It falls with the canopy (above ground, f_labile_leaf).   !
            seed_lost = g_repro * cohort%nplant(j)                                              &
                      * (1.0_wp - min(max(efficiency(j), 0.0_wp), 1.0_wp))
            !----- Only the NON-resorbed share of the leaf and fine-root shed becomes litter;      !
            !      the rest went to storage above, so each pool's full removal is matched          !
            !      exactly by (storage credit + litter input). ------------------------------------!
            call add_necromass(lit(ip), (leaf_shed_c - leaf_retained_c) * cohort%nplant(j),       &
                     (fineroot_shed_c - fineroot_retained_c) * cohort%nplant(j), 0.0_wp, 0.0_wp,   &
                     seed_lost, pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                       &
                     pft%aboveground_frac(pf), pft%struct_lignin_frac(pf))
            !----- ROOT EXUDATE: the carbon the growth sink could not use enters the soil as     !
            !      labile, below-ground carbon, where the microbes respire it. -------------------!
            lit(ip)%labile_soil = lit(ip)%labile_soil + exudate * cohort%nplant(j)
            lit(ip)%exudate_c   = lit(ip)%exudate_c   + exudate * cohort%nplant(j)
         end do
      end associate

      !----- Commit the per-patch owing as a RATE the fast loop will emit over the coming day.     !
      !      KGC_PER_UMOL_C converts kgC -> umol C, since nee_biotic is a molar CO2 flux. --------!
      do ip = 1_ik, site%patch%n
         site%patch%slow_co2_rate(ip) = co2_owed(ip) / max(cfg%dt_slow, tiny_num) / KGC_PER_UMOL_C
      end do
   end subroutine compute_carbon_allocation

   !---------------------------------------------------------------------------------------!
   ! A cohort's growth window as the sink limit sees it this step: the diameter it has already  !
   ! grown in the window [cm] and the window's length [yr] with this step in it. The ring buffer  !
   ! holds one realized dbh rate [cm/yr] per step (count of them, summing to accum); once it is   !
   ! full, this step's sample replaces the oldest, which therefore leaves the window now. A new   !
   ! cohort's window is as long as its life, so it cannot draw on growth it never had time for.   !
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
   ! Leaf/fine-root turnover WATER shedding. A NET biomass LOSS this slow step                         !
   ! (npp%leaf(j)/npp%fineroot(j) < 0 -- shedding out-ran growth) sheds tissue water in the SAME       !
   ! proportion as the carbon loss, so the REMAINING tissue's rwc is unchanged: mass, not psi, is      !
   ! what a biomass change must not silently perturb. A NET GAIN step touches nothing here -- the      !
   ! fast loop's mass-conserving seam (meds_fast_dynamics) handles that direction (capacity            !
   ! grows, mass stays put, rwc reads lower next time). Root sheds into the LUMPED wood node's         !
   ! water (bwood = bsap+broot; only broot sheds under the current turnover model, so                  !
   ! fineroot_carbon's own fractional loss is used as the proxy for the lumped pool's shed             !
   ! fraction -- an approximation, since only the fine-root half of the lumped node actually           !
   ! sheds). MUST run before update_cohort_states commits this step's new                              !
   ! leaf_carbon/fineroot_carbon (the caller sequences it right after compute_carbon_allocation,       !
   ! well before that commit), since the shed FRACTION is computed against the PRE-step pool.          !
   !                                                                                          !
   ! The shed water is credited (nplant-weighted, summed per patch) to site%patch%shed_water_rate      !
   ! -- a frozen daily RATE the fast loop adds to its ground-water input exactly like throughfall      !
   ! (its own distinct variable, not merged into rainfall), so what leaves                             !
   ! leaf_water_mass/wood_water_mass here reappears in the soil in the next fast window instead        !
   ! of vanishing. Its ENERGY needs NO separate wiring: mixed into the ground-water input              !
   ! (precip_ground), it is valued like the precipitation it joins -- "similarly as                    !
   ! precipitation," not a more precise mechanism than precipitation itself receives. -----------------!
   subroutine shed_turnover_water(site, cfg, npp, ledger)
      type(site_t),            intent(inout) :: site
      type(meds_config_t),     intent(in)    :: cfg
      type(carbon_flux_block), intent(in)    :: npp
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp), allocatable :: patch_total(:)
      real(wp) :: shed_frac, shed_amount, shed_enth
      integer(ik) :: j, ip

      allocate(patch_total(site%patch%n))
      patch_total = 0.0_wp
      shed_enth   = 0.0_wp

      associate (cohort => site%cohort)
         do j = 1_ik, cohort%n
            ip = cohort%owner_patch(j)
            if (npp%leaf(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -npp%leaf(j) / max(cohort%leaf_carbon(j), tiny_num))
               shed_amount = cohort%leaf_water_mass(j) * shed_frac
               cohort%leaf_water_mass(j) = cohort%leaf_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
               !----- LEDGER: the leaf's heat capacity carries its internal water, so shedding that  !
               !      water takes sensible heat with it, AT THE TISSUE'S TEMPERATURE. The fast loop  !
               !      puts the mass back into the ground valued as the ground inflow it joins: a     !
               !      re-valuation. Declaring the departure here is what makes that difference       !
               !      measurable rather than absorbed. The WOOD shed contributes nothing: the wood's !
               !      THERMAL water is the sapwood ring at a fixed moisture fraction, not            !
               !      wood_water_mass, so draining the hydraulic store moves no heat.                !
               shed_enth = shed_enth + site%patch%area(ip) * cohort%nplant(j) * shed_amount          &
                                     * cp_liq * cohort%leaf_temp(j)
            end if
            if (npp%fineroot(j) < 0.0_wp) then
               shed_frac   = min(1.0_wp, -npp%fineroot(j) / max(cohort%fineroot_carbon(j), tiny_num))
               shed_amount = cohort%wood_water_mass(j) * shed_frac
               cohort%wood_water_mass(j) = cohort%wood_water_mass(j) - shed_amount
               patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
            end if
         end do
      end associate

      do ip = 1_ik, site%patch%n
         site%patch%shed_water_rate(ip) = patch_total(ip) / max(cfg%dt_slow, tiny_num)
      end do

      !----- LEDGER: the shed is a HANDOFF to the fast tier, not a store this ledger carries (see   !
      !      slow_site_store's header on why treating it as one double-counts at every open).       !
      if (present(ledger)) then
         do ip = 1_ik, site%patch%n
            call slow_ledger_declare(ledger, water_out = site%patch%area(ip) * patch_total(ip))
         end do
         call slow_ledger_declare(ledger, energy_out = shed_enth)
      end if
   end subroutine shed_turnover_water

   !---------------------------------------------------------------------------------------!
   ! New leaves take their water from the wood. When a cohort's leaf carbon grows, the new     !
   ! leaf tissue is filled with the water it holds at rest with the wood store, moved from    !
   ! that store, so the plant's water is unchanged: the new tissue and the wood end at the    !
   ! potential they share (shared_potential), the leaf the crown's gravity head below the     !
   ! wood, as the network holds them with no flow (docs/science/plant_hydraulics.md eqs. (2), !
   ! (5)). New tissue that started empty would read at the dry end of its pressure-volume     !
   ! curve: a canopy re-flushing from bare would start the next fast step near -1.5e4 MPa,    !
   ! its stomata shut for that step and the canopy's mean leaf potential dragged with it. The !
   ! water moved is bounded by what the wood holds. Runs right after the commit, on the leaf  !
   ! carbon it added (docs/science/soil_plant_water.md §2).                                   !
   !---------------------------------------------------------------------------------------!
   subroutine fill_new_leaf_water(site, cfg, leaf_carbon_before)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: leaf_carbon_before(:)   !< [kgC/plant] before the commit
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: grown, head, psi_wood, moved
      integer(ik) :: i

      associate (cohort => site%cohort)
         do i = 1_ik, cohort%n
            grown = cohort%leaf_carbon(i) - leaf_carbon_before(i)
            if (grown <= 0.0_wp) cycle
            call pft_tissue_curves(cfg%hydraulics, cfg%pft, cohort%pft(i), leaf_curve, wood_curve)
            head     = grav_head * cohort%height(i)          ! the crown above the wood datum [MPa]
            psi_wood = shared_potential(cohort%wood_water_mass(i), wood_curve,                       &
                          cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)), &
                          leaf_curve, grown, head)
            moved    = min(water_content(psi_wood - head, leaf_curve, grown), cohort%wood_water_mass(i))
            cohort%leaf_water_mass(i) = cohort%leaf_water_mass(i) + moved
            cohort%wood_water_mass(i) = cohort%wood_water_mass(i) - moved
         end do
      end associate
   end subroutine fill_new_leaf_water

   !---------------------------------------------------------------------------------------!
   ! Tissue water ABOVE the capacity a cohort's biomass now allows goes to the ground. A store   !
   ! keeps its water mass across a change of biomass (docs/science/soil_plant_water.md §2), and   !
   ! most changes leave it within its capacity: turnover takes its tissue's water with it, and    !
   ! new leaves are filled from the wood. A few shrink a store's capacity without touching its   !
   ! water -- a rising SLA resorbing leaf carbon to storage, a census start's acclimation, a      !
   ! fusion or a split re-deriving the sapwood ring -- and water above saturation is not a state  !
   ! a tissue can hold. The excess goes down the channel the turnover shed uses                   !
   ! (patch%shed_water_rate), so it is moved, not deleted; the caller declares it to the ledger.  !
   !---------------------------------------------------------------------------------------!
   subroutine shed_excess_tissue_water(site, cfg, water_shed)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(out)   :: water_shed     !< [kg/m2 site] the water shed
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp)    :: leaf_kept, wood_kept, excess
      integer(ik) :: i, ip

      water_shed = 0.0_wp
      associate (cohort => site%cohort, patch => site%patch)
         do i = 1_ik, cohort%n
            call pft_tissue_curves(cfg%hydraulics, cfg%pft, cohort%pft(i), leaf_curve, wood_curve)
            leaf_kept = clamp_water_to_capacity(cohort%leaf_water_mass(i), leaf_curve, cohort%leaf_carbon(i))
            wood_kept = clamp_water_to_capacity(cohort%wood_water_mass(i), wood_curve,                &
                           cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
            excess = (cohort%leaf_water_mass(i) - leaf_kept) + (cohort%wood_water_mass(i) - wood_kept)
            if (excess <= 0.0_wp) cycle
            cohort%leaf_water_mass(i) = leaf_kept
            cohort%wood_water_mass(i) = wood_kept
            ip = cohort%owner_patch(i)
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip)                                   &
                                      + cohort%nplant(i) * excess / max(cfg%dt_slow, tiny_num)
            water_shed = water_shed + patch%area(ip) * cohort%nplant(i) * excess
         end do
      end associate
   end subroutine shed_excess_tissue_water

   !---------------------------------------------------------------------------------------!
   ! Continuous background-mortality WATER. The individuals that die this step take their          !
   ! tissue water with them; without this routine it would vanish when `update_cohort_states`      !
   ! lowers nplant -- a silent sink on every cohort on every step. (Their tissue HEAT is a change  !
   ! in thermal MASS, declared with the growth commit's, since the two are one mechanism.)         !
   !                                                                                          !
   ! The water goes to the ground down the SAME channel turnover shedding already uses, so there  !
   ! is one verified path for "tissue water that stopped being tissue water" rather than a second  !
   ! mechanism to keep in step with the first.                                                      !
   !                                                                                          !
   ! Only the TISSUE water moves. The interception film is per m2 of ground and `update_cohort_     !
   ! states` does not scale it when nplant falls, so no film leaves the store here -- routing it     !
   ! would invent a flux the state never made.                                                       !
   !                                                                                          !
   ! Runs AFTER update_cohort_states commits the new nplant, with the density from BEFORE the       !
   ! commit (nplant_before): the dying number is n0 - n1, exactly as accumulate_mortality_litter    !
   ! measures its carbon.                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine shed_mortality_water(site, cfg, nplant_before, water_shed)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: nplant_before(:) !< [plant/m2] density BEFORE the commit
      real(wp),            intent(out)   :: water_shed       !< [kg/m2 site]
      integer(ik) :: j, ip
      real(wp)    :: died_nplant, w_j

      water_shed = 0.0_wp
      associate (cohort => site%cohort, patch => site%patch)
         do j = 1_ik, cohort%n
            died_nplant = nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            ip  = cohort%owner_patch(j)
            !----- The water that LEFT the store is the density drop times the PER-PLANT water,    !
            !      not a fraction of what the survivors still hold. cohort_tissue_water carries    !
            !      the post-commit nplant, so using it here would scale the loss by n1/n0. --------!
            w_j = died_nplant * (max(cohort%leaf_water_mass(j), 0.0_wp)                            &
                               + max(cohort%wood_water_mass(j), 0.0_wp))
            patch%shed_water_rate(ip) = patch%shed_water_rate(ip) + w_j / max(cfg%dt_slow, tiny_num)
            water_shed = water_shed + patch%area(ip) * w_j
         end do
      end associate
   end subroutine shed_mortality_water

   !---------------------------------------------------------------------------------------!
   ! Credit one step's recruitment into the per-(PFT, patch) carry-forward pool, and report the  !
   ! part of it that came from OUTSIDE the site.                                                  !
   !                                                                                          !
   ! `recruitment` already carries both contributions: the stand's own reproduction NPP, which    !
   ! the parents were debited for in this same step, and the baseline SEED RAIN, which no parent  !
   ! here paid for. Only the second is an external import, and only it is reported -- the first   !
   ! is an internal transfer from the plants to the pool, and both ends of it happen in the same  !
   ! phase, so the ledger sees it close rather than seeing carbon vanish in one phase and appear  !
   ! in another.                                                                                  !
   !                                                                                          !
   ! `include_pft` gates the seed rain here exactly as compute_vital_rates gates it into the rate,!
   ! so a disabled PFT contributes neither recruits nor an import.                                !
   !---------------------------------------------------------------------------------------!
   subroutine accumulate_recruit_pool(site, cfg, recruitment, dt_yr, seed_rain_carbon)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      real(wp),            intent(in)    :: dt_yr
      real(wp),            intent(out)   :: seed_rain_carbon  !< [kgC/m2 site] the external part
      integer(ik) :: ip, pf
      real(wp)    :: c_min

      seed_rain_carbon = 0.0_wp
      associate (patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, patch%n
            do pf = 1_ik, site%n_pft
               patch%recruit_pool(pf, ip) = patch%recruit_pool(pf, ip) + recruitment(pf, ip) * dt_yr
               if (pft%include_pft(pf) /= 1_ik) cycle
               c_min = min_cohort_carbon(site%cohort%allom, pft%min_cohort_height, pft%wood_density(pf))
               seed_rain_carbon = seed_rain_carbon                                                 &
                                + patch%area(ip) * pft%seed_rain_recruits(pf) * dt_yr * c_min
            end do
         end do
      end associate
   end subroutine accumulate_recruit_pool

   !---------------------------------------------------------------------------------------!
   ! Continuous background-mortality LITTER: the carbon carried by the individuals that died        !
   ! this step. Runs AFTER update_cohort_states, and both of its inputs are the reason why.         !
   !                                                                                          !
   ! The live store is nplant*pool, and its change decomposes exactly as                            !
   !                                                                                          !
   !     n1.p1 - n0.p0  =  n0.(p1 - p0)  +  (n1 - n0).p1                                            !
   !                                                                                          !
   ! -- growth at the OLD density, mortality at the NEW pools. So the carbon the applier actually   !
   ! removed for mortality is (n0-n1).p1; valuing the litter on the PRE-growth pools p0 would       !
   ! leave n_died x (this step's per-plant growth) unaccounted every step of every cohort, small    !
   ! per step and one-signed: the shape that only shows up as a season's drift.                     !
   !                                                                                          !
   ! The density drop is taken as (n0 - n1) rather than recomputed from the hazard, so when the     !
   ! `negligible_nplant` floor stops a cohort dying, litter is not credited for individuals that    !
   ! are still standing.                                                                            !
   !                                                                                          !
   ! Whole-individual death carries EVERY pool (leaf/fine-root/wood/storage), unlike a turnover     !
   ! shed. Scatters into the SAME per-patch `lit` accumulator the turnover litter uses. ------------!
   subroutine accumulate_mortality_litter(site, cfg, nplant_before, lit)
      type(site_t),          intent(inout) :: site      !< inout for the mortality-carbon diagnostic
      type(meds_config_t),   intent(in)    :: cfg
      real(wp),              intent(in)    :: nplant_before(:) !< [plant/m2] density BEFORE the commit
      type(litter_input_t),  intent(inout) :: lit(:)
      integer(ik) :: j, pf, ip
      real(wp)    :: died_nplant

      associate (cohort => site%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            died_nplant = nplant_before(j) - cohort%nplant(j)
            if (died_nplant <= 0.0_wp) cycle
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !----- The BACKGROUND pathway's carbon, on the same density and the same pools           !
            !      the litter below is valued on -- one measurement, two consumers, so the pathway    !
            !      split and the litter can never disagree: how much biomass died is demography, not  !
            !      biogeochemistry. As a rate [kgC/m2/yr] times the step's weight in seconds, which   !
            !      is the block's contract.                                                           !
            if (site%patch%diag%active)                                                            &
               site%patch%diag%v(PD_MORT_C_BACKGROUND, ip) =                                       &
                  site%patch%diag%v(PD_MORT_C_BACKGROUND, ip)                                      &
                  + (died_nplant * (cohort%leaf_carbon(j) + cohort%fineroot_carbon(j)              &
                                  + cohort%wood_carbon(j) + cohort%nonstructural_carbon(j))        &
                     / cfg%dt_years) * cfg%dt_slow
            call add_necromass(lit(ip), died_nplant * cohort%leaf_carbon(j),                      &
                     died_nplant * cohort%fineroot_carbon(j), died_nplant * cohort%wood_carbon(j),  &
                     died_nplant * cohort%nonstructural_carbon(j), 0.0_wp,                         &
                     pft%f_labile_leaf(pf), pft%f_labile_stem(pf),                                 &
                     pft%aboveground_frac(pf), pft%struct_lignin_frac(pf))
         end do
      end associate
   end subroutine accumulate_mortality_litter

   !---------------------------------------------------------------------------------------!
   ! Per-cohort carbon-demand assembler: gathers this step's GPP and maintenance respiration   !
   ! (fast-loop accumulators, or the gpp_ref stub), applies TURNOVER FIRST -- the leaf and      !
   ! fine-root losses, this step's litter -- and forms the flush-capped growth demands against   !
   ! the POST-SHED pools. The leaf loss and the flush cap come from the phenology kernel's       !
   ! leaf_turnover_step, the same rule the Python mirror applies: senescence at shed_rate_max *  !
   ! leaf_shed_tendency down to min_leaf_cover of the full canopy, plus background turnover      !
   ! while the canopy flushes and drops no leaves. Fine roots do not drop with the leaves: they   !
   ! keep the plant's water uptake (root length sets the root-soil conductance), and only their  !
   ! background turnover follows the leaves' gate, so they age and regrow with the canopy.       !
   ! Every input is a plain scalar (the caller's per-kernel-call convention), so no cohort/PFT   !
   ! SoA type needs to be named here (docs/science/plant_phenology.md eq. (9) and §5). ---------!
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
      !----- The share of this step's leaf and fine-root LOSS (senescence and background         !
      !      turnover alike) resorbed back into the non-structural pool. Reported separately from  !
      !      leaf_shed_c and fineroot_shed_c, which stay the FULL removals from the two pools --    !
      !      crediting storage while removing only the litter share would CREATE carbon. ---------!
      real(wp), intent(in)  :: retained_frac        !< [-] retained_carbon_fraction for this PFT
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
      !----- TURNOVER FIRST: the leaf and fine-root losses (litter), then the POST-SHED pools. --!
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
      !----- Flush-capped GROWTH demands toward target, from the post-shed pools. -------------!
      leaf_demand     = min(max(0.0_wp, leaf_target       - leaf_post), flush_cap)
      fineroot_demand = min(max(0.0_wp, r2l * leaf_target - root_post), flush_cap * r2l)
   end subroutine cohort_carbon_demand

   !---------------------------------------------------------------------------------------!
   ! Advance the leaf phenology of every cohort over one slow step (ED2's phenology_driv          !
   ! analogue; docs/science/plant_phenology.md). It flattens the per-PFT cue params, builds the   !
   ! daily env from what the fast loop accumulated (the site's air temperature, each cohort's     !
   ! PAR at its top) + the cohort's predawn leaf psi + latitude + day-of-year, advances the       !
   ! kernel, and writes the tendencies and cue memory back to the cohort. It touches NO           !
   ! leaf/storage carbon; compute_carbon_allocation turns the stored tendencies into leaf         !
   ! growth and loss. A step with no fast sub-steps is skipped so the memory is never advanced    !
   ! on a bogus 0/0 mean.                                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine advance_leaf_phenology(site, cfg, doy, latitude_deg)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: doy
      !----- The polygon's latitude (its hemisphere). A region passes each polygon's own;           !
      !      absent, it is the run's [site] latitude. -----------------------------------------------!
      real(wp),            intent(in), optional :: latitude_deg
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: i, pf
      real(wp)    :: dt_days, temp_day, lat
      logical     :: north

      if (site%pheno_tair_n < 1_ik) return             ! no fast sub-steps this slow step -> no drivers
      dt_days  = cfg%dt_slow / day_sec
      temp_day = site%pheno_tair_sum / real(site%pheno_tair_n, wp)
      lat = cfg%forcing%latitude_deg
      if (present(latitude_deg)) lat = latitude_deg
      north    = lat >= 0.0_wp

      do i = 1_ik, site%cohort%n
         pf = site%cohort%pft(i)
         call flatten_pheno_params(cfg, pf, params)

         !----- Daily environment. `dmax_psi_leaf` is the cohort's own published predawn value      !
         !      -- already a completed daily maximum, double-buffered against its accumulator.      !
         !      An UNSET cohort (born today) reports 0, which reads as well-watered for one day;     !
         !      that is the same seeding convention the leaf kernel uses.                            !
         env%temp_day         = temp_day
         !----- light_hours_accum is the time [h] the PAR at the cohort's top exceeded its PFT's   !
         !      par_min over the step (counted in the fast loop); per day for the cue. -----------!
         env%par_hours        = site%cohort%light_hours_accum(i) / dt_days
         env%predawn_leaf_psi = min(0.0_wp, site%cohort%dmax_psi_leaf(i))
         env%doy              = doy
         env%hemis_north      = north

         !----- Pack the cohort's tendencies + cue memory, advance, unpack. ---------------------!
         state%leaf_flush_tendency = site%cohort%leaf_flush_tendency(i)
         state%leaf_shed_tendency  = site%cohort%leaf_shed_tendency(i)
         state%growing_degree_days = site%cohort%growing_degree_days(i)
         state%cold_degree_days    = site%cohort%cold_degree_days(i)
         state%wet_psi_sum         = site%cohort%wet_psi_sum(i)
         state%dry_psi_sum         = site%cohort%dry_psi_sum(i)
         state%light_hours_mean    = site%cohort%light_hours_mean(i)
         call phenology_kernel(env, params, dt_days, state, out)
         site%cohort%leaf_flush_tendency(i) = state%leaf_flush_tendency
         site%cohort%leaf_shed_tendency(i)  = state%leaf_shed_tendency
         site%cohort%growing_degree_days(i) = state%growing_degree_days
         site%cohort%cold_degree_days(i)    = state%cold_degree_days
         site%cohort%wet_psi_sum(i)         = state%wet_psi_sum
         site%cohort%dry_psi_sum(i)         = state%dry_psi_sum
         site%cohort%light_hours_mean(i)    = state%light_hours_mean
      end do
   end subroutine advance_leaf_phenology

   !---------------------------------------------------------------------------------------!
   ! Light trait plasticity for every cohort over one slow step: acclimate the per-cohort leaf   !
   ! traits (sla / vcmax25 / rd25 / llspan) toward their shaded targets from the cumulative LAI   !
   ! above (overtopping_lai), relaxing at the leaf-replacement rate 1/llspan (turnover-limited).  !
   ! Targets come from the PFT top-of-canopy values + the derived plasticity slopes; the live     !
   ! cohort traits are the state being mutated. Called only when trait_plasticity_on.             !
   ! docs/science/plant_traits.md eqs. (1)-(2) and §4.                                             !
   !---------------------------------------------------------------------------------------!
   subroutine advance_plant_traits(site, cfg, dt_yr, instantaneous)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      logical, optional,   intent(in)    :: instantaneous   !< census restart: jump straight to the target
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, llspan_t, ll_now, leaf_target, excess
      logical     :: instant
      instant = .false.
      if (present(instantaneous)) instant = instantaneous
      associate (cohort => site%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            !----- Light-acclimated targets from this cohort's cumulative LAI above. -----------!
            call light_plastic_traits(cohort%overtopping_lai(j),                                 &
                     pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_lifespan_toc(pf),       &
                     pft%kplastic_sla(pf), pft%kplastic_vm0(pf), pft%kplastic_rd(pf),             &
                     pft%kplastic_llspan(pf), sla_t, vcmax_t, rd_t, llspan_t)
            !----- Update each live trait toward its target (gradual at the leaf-replacement rate,  !
            !      or an instantaneous jump on a census restart). ---------------------------------!
            ll_now = cohort%llspan(j)
            cohort%sla(j)     = update_plastic_trait(cohort%sla(j),     sla_t,    ll_now, dt_yr, instant)
            cohort%vcmax25(j) = update_plastic_trait(cohort%vcmax25(j), vcmax_t,  ll_now, dt_yr, instant)
            cohort%rd25(j)    = update_plastic_trait(cohort%rd25(j),    rd_t,     ll_now, dt_yr, instant)
            cohort%llspan(j)  = update_plastic_trait(cohort%llspan(j),  llspan_t, ll_now, dt_yr, instant)
            !----- SLA rose => leaf_carbon now exceeds the (lower) allometric leaf target: resorb   !
            !      the excess to storage so displayed leaf AREA stays at allometry (leaf-area-       !
            !      conserving). An UNshaded step lowers sla => a deficit, filled by normal growth.  !
            leaf_target = size2leaf_carbon(cohort%allom, cohort%dbh(j), cohort%height(j), cohort%sla(j))
            if (cohort%leaf_carbon(j) > leaf_target) then
               excess = cohort%leaf_carbon(j) - leaf_target
               cohort%leaf_carbon(j)          = leaf_target
               cohort%nonstructural_carbon(j) = cohort%nonstructural_carbon(j) + excess
            end if
            call set_cohort_size_from_carbon(cohort, j)   ! leaf_area = sla*leaf_carbon, either way
         end do
      end associate
   end subroutine advance_plant_traits

   !---------------------------------------------------------------------------------------!
   ! reacclimate_plant_traits -- a restart's plastic leaf traits from THIS run's PFT file      !
   ! ([init].reacclimate_traits). A state file carries the traits the writing run's PFT file   !
   ! gave, so without this a restart with a changed vcmax25 would run on the old one. With      !
   ! plasticity on, each cohort takes the census start's instant acclimation to the LAI above   !
   ! it as the file holds it; with plasticity off, the PFT's top-of-canopy values, as a census   !
   ! start gives. The stand's STRUCTURE is kept: leaf carbon scales by the SLA's change (old /   !
   ! new), storage taking or giving the difference, so leaf area stays as read to round-off and  !
   ! a trial that changes a trait changes the leaf physiology and not the LAI -- and a restart   !
   ! whose PFT file gives the traits it already has changes nothing. A lower SLA asks storage    !
   ! for more leaf carbon than it may hold: the leaves then take only what storage has, and the  !
   ! leaf area falls short by the rest, so no carbon is made.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine reacclimate_plant_traits(site, cfg)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: j, pf
      real(wp)    :: sla_t, vcmax_t, rd_t, llspan_t, leaf_carbon_new, sla_old
      associate (cohort => site%cohort, pft => cfg%pft)
         do j = 1_ik, cohort%n
            pf = cohort%pft(j)
            if (cfg%trait_plasticity_on) then
               call light_plastic_traits(cohort%overtopping_lai(j),                              &
                        pft%sla(pf), pft%vcmax25(pf), pft%rd25(pf), pft%leaf_lifespan_toc(pf),    &
                        pft%kplastic_sla(pf), pft%kplastic_vm0(pf), pft%kplastic_rd(pf),          &
                        pft%kplastic_llspan(pf), sla_t, vcmax_t, rd_t, llspan_t)
            else
               sla_t = pft%sla(pf) ; vcmax_t = pft%vcmax25(pf) ; rd_t = pft%rd25(pf)
               llspan_t = pft%leaf_lifespan_toc(pf)
            end if
            sla_old           = cohort%sla(j)
            cohort%sla(j)     = sla_t
            cohort%vcmax25(j) = vcmax_t
            cohort%rd25(j)    = rd_t
            cohort%llspan(j)  = llspan_t
            if (cohort%sla(j) == sla_old) cycle
            leaf_carbon_new = cohort%leaf_carbon(j) * sla_old / max(cohort%sla(j), tiny_num)
            leaf_carbon_new = min(leaf_carbon_new, cohort%leaf_carbon(j) + cohort%nonstructural_carbon(j))
            cohort%nonstructural_carbon(j) = cohort%nonstructural_carbon(j) + cohort%leaf_carbon(j)        &
                                             - leaf_carbon_new
            cohort%leaf_carbon(j) = leaf_carbon_new
            call set_cohort_size_from_carbon(cohort, j)
         end do
      end associate
   end subroutine reacclimate_plant_traits

   !---------------------------------------------------------------------------------------!
   ! Flatten the per-PFT phenology traits (cfg%pft%pheno_*) into the self-contained kernel param  !
   ! set. Mirrors meds_fast_config's leaf-trait flattening.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine flatten_pheno_params(cfg, ipft, p)
      type(meds_config_t),  intent(in)  :: cfg
      integer(ik),          intent(in)  :: ipft
      type(pheno_params_t), intent(out) :: p
      associate (t => cfg%pft)
         p%flush_cue_mask        = t%pheno_flush_cue_mask(ipft)
         p%shed_cue_mask         = t%pheno_shed_cue_mask(ipft)
         p%flush_cue_timescale   = t%pheno_flush_cue_timescale(ipft)
         p%shed_cue_timescale    = t%pheno_shed_cue_timescale(ipft)
         p%flush_rate_max        = t%pheno_flush_rate_max(ipft)
         p%shed_rate_max         = t%pheno_shed_rate_max(ipft)
         p%flush_base_temp       = t%pheno_flush_base_temp(ipft)
         p%flush_degree_days     = t%pheno_flush_degree_days(ipft)
         p%flush_temp_sharpness  = t%pheno_flush_temp_sharpness(ipft)
         p%shed_base_temp        = t%pheno_shed_base_temp(ipft)
         p%shed_degree_days      = t%pheno_shed_degree_days(ipft)
         p%shed_temp_sharpness   = t%pheno_shed_temp_sharpness(ipft)
         p%par_min               = t%pheno_par_min(ipft)
         p%flush_light_hours     = t%pheno_flush_light_hours(ipft)
         p%flush_light_sharpness = t%pheno_flush_light_sharpness(ipft)
         p%shed_light_hours      = t%pheno_shed_light_hours(ipft)
         p%shed_light_sharpness  = t%pheno_shed_light_sharpness(ipft)
         p%light_window          = t%pheno_light_window(ipft)
         p%flush_water_sum       = t%pheno_flush_water_sum(ipft)
         p%flush_water_sharpness = t%pheno_flush_water_sharpness(ipft)
         p%shed_water_sum        = t%pheno_shed_water_sum(ipft)
         p%shed_water_sharpness  = t%pheno_shed_water_sharpness(ipft)
      end associate
      !----- The turgor-loss point is DERIVED, not a phenology key: the water cue compares the       !
      !      predawn leaf psi against the same psi_tlp the leaf gas-exchange kernel builds from the  !
      !      pressure-volume curve. One authority, so the cue cannot drift from the leaf's water-    !
      !      stress controls built on the same threshold: the fast capacity ramp and the slow        !
      !      stomatal factor (docs/science/leaf_gas_exchange.md §4.1 eq. (9), §4.2 eq. (10)). -------!
      p%leaf_psi_tlp = pft_leaf_psi_tlp(cfg, ipft)
   end subroutine flatten_pheno_params

end module meds_vegetation_dynamics
