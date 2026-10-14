# The order of processes

This page is the map: which process runs when, and what each one sees. The other science pages
give each process's equations; this one gives their sequence. Three clocks nest:

- a **fast sub-step** `dt_fast` (900 s in the shipped configs) advances the surface column: canopy air, leaves and
  wood, soil and snow;
- a **slow step** `dt_slow` (one day) runs `n_fast_per_slow` fast sub-steps, then the vegetation and
  soil-carbon dynamics, then the output;
- the **calendar boundaries** between two slow steps restructure the stand: monthly for cohorts,
  yearly for patches. At each month boundary, the **I/O phase** writes the files.

Everything below is per polygon. A site run is one polygon; a region run is many, each stepped the
same way ([§1](#1-a-run)).

---

## 1. A run

**Setup** (`driver_open`):

1. Read the configuration.
2. Build the initial community from `[init].init_mode`: bare ground, a census, or a restart. A restart
   also restores the calendar date and any restructuring the checkpoint still owes ([§5](#5-the-calendar-boundary)).
3. Open the forcing source. Prepare the polygon: the fast context, the forcing cursor, and the fast
   reservoirs (unless the restart restored them); the initial snow; the canopy-air depth, from the
   stand; the steady-state soil carbon (unless restored); and the slow ledger.
4. Open the state (restart) stream and the output manager.

**The loop.** Until `end_time`, each `driver_step`:

1. Prefetches the forcing the step needs; the archive backend loads a month at a time.
2. Runs one slow step of the polygon ([§2](#2-one-slow-step)).
3. At a month boundary, runs the I/O phase ([§6](#6-output-and-io)).

**Finish** (`driver_finalize`): the terminal checkpoint, the summary, the conservation reports (the
whole-column energy and water budgets, the face closure, the soil-carbon seam, the slow ledger), and
the output's final partial windows.

**A region run** steps a month at a time (`region_step_month`):

1. Prefetch the month's forcing once.
2. The compute phase: take each polygon through all of the month's steps. It uses the same
   `polygon_step` as a site run, and touches no file.
3. The I/O phase: the region files, then each detail polygon's own files.

Each polygon's pending restructuring ([§5](#5-the-calendar-boundary)) carries from one month to the
next.

---

## 2. One slow step

`polygon_step` advances one polygon from `prev` to `now`, in this order:

1. **The boundary the last step ended on**, if it owes one (`restructure_pending`): the stand's
   monthly and yearly restructuring ([§5](#5-the-calendar-boundary)).
2. **Growth-temperature acclimation** (opt-in, `leaf_thermal_acclimation`): advance the running mean
   from the previous step's daily mean air temperature, and refresh the leaf photosynthesis table.
3. **The fast loop** ([§3](#3-the-fast-loop)), when `fast_biophysics_on`.
4. **The slow dynamics** ([§4](#4-the-slow-dynamics)), when `slow_on`.
5. **Output.**
   - Replay the fast loop's staged sub-step samples into the fast tier.
   - Fold the step into the daily, monthly and annual windows that hold `prev`, then close each period
     `now` has left. A record holds the steps that start in its period ([diagnostics](diagnostics.md) §4).
6. **Zero the step diagnostic blocks**, which the output has now read. If `now` is a month boundary,
   set `restructure_pending`, and at a year boundary `restructure_new_year` as well.
7. **Guards**:
   - a NaN check at each year boundary;
   - a soil-carbon plausibility check every step. It is unconditional, while the slow ledger's copy of
     the check can be switched off.

---

## 3. The fast loop

`fast_dynamics` runs once per slow step.

**Before the sub-steps, once per slow step:**

1. **Reset the fast-to-slow accumulators** the slow dynamics will read:
   - gross GPP and maintenance respiration per cohort;
   - the soil-carbon environment and heterotrophic-respiration integrals;
   - the daily air-temperature sum used by phenology.

   Also roll over the predawn leaf water potential.
2. **Size the diagnostic blocks** to today's stand. They are not zeroed here ([§6](#6-output-and-io)).
3. **Reconcile each cohort's stored tissue water** with the capacity today's biomass allows.
4. **Sample the forcing** for every sub-step, at `prev + (k − 1 + forcing_sample_frac)·dt_fast` for
   `k = 1 … n_fast_per_slow` ([forcing](forcing.md)). The CO₂ is looked up at the same instants,
   on model time ([forcing](forcing.md) §12). Accumulate the polygon's forcing diagnostics.

**Then, for each patch:**

The patches are independent here, so this loop is threaded ([numerical_scheme](numerical_scheme.md)
§6a). Each patch first gathers its cohorts into the column buffer with their canopy geometry, and takes
its canopy-air top, roughness and displacement for the day. It then freezes, for the day, its
soil-carbon pools and the tissue water the slow step shed.

**Each fast sub-step `k`, in order:**

1. **Forcing.** The sub-step's meteorological sample, which already carries the terrain lapse. Its air
   temperature joins the phenology sum. The patch takes a copy moved to its own canopy-air top: the wind
   along its log profile, the temperature conserving potential temperature ([forcing](forcing.md) §8).
   The aerodynamics and the canopy air's exchange with the atmosphere use that copy.
2. **Radiation.** With forcing on, the canopy radiation transfer splits absorbed shortwave and
   longwave among the cohorts and the ground ([canopy_radiation_transfer](canopy_radiation_transfer.md)).
   It runs in the driver (`apply_rt_forcing`), before the column step. With forcing off, the shortwave
   is split by LAI share.
3. **The aerodynamic environment** above the canopy ([canopy_aerodynamics](canopy_aerodynamics.md)).
4. **The coupled column step** (`column_fast_step`, [column_biophysics](column_biophysics.md)).
   - Its pre-pass computes leaf gas exchange (GPP, stomatal conductance, leaf respiration), stem and
     root respiration, the canopy-air capacities, the aerodynamics and the plant-hydraulics solve.
   - These are frozen for the sub-step, except the canopy-air-to-atmosphere conductances, which are
     re-solved at every stage.
   - Then one integrator (`ark` or `rk45`, [numerical_scheme](numerical_scheme.md)) advances the
     canopy air, leaf and wood energy and water, the soil column and the snow together.
5. **Accumulate.**
   - The patch's evapotranspiration, the phenology cues, the integrator's work counters, and the
     soil-carbon environment scalar and heterotrophic respiration.
   - The fast-tier samples.
   - The step diagnostic blocks, each weighted by `dt_fast`.

After the sub-steps, the run's whole-column energy and water budgets fold in the step's residuals.

---

## 4. The slow dynamics

`advance_slow_dynamics` opens the slow ledger's window for the step, then runs
`vegetation_dynamics`, which advances the cohorts one day:

1. **Leaf phenology**: one daily step of the flush and senescence tendencies, from the step's start
   day of year, the day's mean air temperature, and each cohort's daily PAR at its top and predawn
   leaf water potential ([plant_phenology](plant_phenology.md)).
2. **Leaf-trait plasticity** (opt-in, `trait_plasticity_on`): acclimate the traits toward their shaded
   targets ([plant_traits](plant_traits.md)).
3. **Carbon allocation**:
   - the day's net carbon from the fast loop's GPP and maintenance respiration, through the
     allocation ladder ([plant_carbon_allocation](plant_carbon_allocation.md),
     [plant_respiration](plant_respiration.md));
   - leaf and fine-root turnover into the patch's litter accumulator;
   - the turnover's tissue water shed to the ground.
4. **Vital rates**: growth, mortality and recruitment. Credit the patch's recruit pool with the day's
   recruitment.
5. **Growth and mortality commit**: the cohorts' tendencies are applied. The dead plants' carbon goes
   to litter and their water to the ground.
6. **Re-sort** the cohorts by height. Record the day's litterfall and recruitment diagnostics, and
   the overtopping LAI.

Then, in `advance_slow_dynamics`:

7. **Patch ageing** by one day.
8. **The canopy-air depth refresh**: each patch's canopy-air volume is resized to its tallest cohort,
   with the entrained or detrained air booked as an exchange.
9. **Soil biogeochemistry** (`soil_carbon_on`): the daily CENTURY step, which consumes the litter
   accumulator and the fast loop's environmental integral ([soil_carbon](soil_carbon.md)).

The slow ledger marks each phase, so a conservation residual is attributed to the operator that made
it.

---

## 5. The calendar boundary

A slow step whose end `now` lies in a new month leaves the stand's restructuring **owed**. The next
step performs it first, before its fast loop (`advance_boundary` → `restructure_stand`). Between the
two steps sit the output tick, the I/O phase and any checkpoint. The restructuring runs only with
`[run].slow_on` (the master freeze of the slow tier), and each operator has its own switch besides.

**At every month boundary** (`demography_on` and `do_cohort_fissfuse`), in order:
1. recruitment from the recruit pools;
2. cohort fusion;
3. culling of cohorts below the tracking floor;
4. splitting;
5. a re-sort.

**At a year boundary**, after the month's operators:
1. patch disturbance (`do_patch_disturbance`);
2. patch fusion and termination (`do_patch_fissfuse`), followed by the cohort fusion, culling and
   re-sort they call for.

**Then**, for the restructured stand: the overtopping LAI, and the canopy-air depth.

Three consequences:

- **Every output record describes one stand.** The restructuring happens between two records, never
  inside one.
- **A boundary's events belong to the period it opens.** The disturbance area and mortality of 1
  January are in January's record, and the biomass the disturbance removes first shows in the 1
  January daily state ([diagnostics](diagnostics.md) §4).
- **A checkpoint on a boundary holds the stand before the restructuring.** Its global attribute
  `restructure_pending` (`none`, `month` or `year`) records what is owed, and a run resumed from it
  performs that first, as the continuous run does. A state file without the attribute owes nothing.

---

## 6. Output and I/O

- **The step diagnostic blocks** (per cohort, per patch, per polygon) run from one output tick to the
  next:
  - the fast sub-steps add to them;
  - so do the slow step's rows (litter, recruitment, mortality);
  - so do the events of a boundary between two steps (disturbance, culls).

  They are zeroed right after the tick. The slow per-cohort block is set afresh each step instead.
- **The fast tier** is staged inside the fast loop, one sample per sub-step, and replayed at the tick.
- **The daily, monthly and annual tiers** fold once per slow step at the tick, and stage each closed
  period in a queue. No netCDF call is made inside a step.
- **The I/O phase** runs at each month boundary, after the step and before the next one:
  - it writes the queued records;
  - in a site run, every `[state].interval_years`, it writes the checkpoint at the year boundary
    (a region writes none until its restarts exist).

  The terminal checkpoint is written at the end of the run.

---

## 7. Code map

| Stage | Routine | Module |
|---|---|---|
| run setup, loop, finish | `driver_open`, `driver_step`, `driver_io_phase`, `driver_finalize` | `meds_driver` |
| region month | `region_step_month`, `io_phase` | `meds_region` |
| one slow step | `polygon_step` | `meds_polygon` |
| fast loop + slow dynamics | `advance_one_step` | `meds_stepper` |
| the fast loop | `fast_dynamics` | `meds_fast_dynamics` |
| vertical corrections: the terrain lapse at ingest, the move to each patch's canopy-air top | `read_record`; `met_to_cas_top` | `meds_met_driver`, `meds_lapse_rate` |
| one sub-step of one patch | `column_fast_step` | `meds_fast_step` |
| the slow dynamics | `advance_slow_dynamics` | `meds_slow_dynamics` |
| the cohorts' day | `vegetation_dynamics` | `meds_vegetation_dynamics` |
| soil carbon | `advance_biogeochem_dynamics` | `meds_biogeochem_dynamics` |
| the calendar boundary | `advance_boundary` → `advance_boundary_dynamics` → `restructure_stand` | `meds_stepper`, `meds_slow_dynamics`, `meds_vegetation_dynamics` |
| the output tick | `output_integrate`, `output_integrate_fast`, `close_tier` | `meds_output_integrate` |
| zeroing the step diagnostics | `reset_step_diagnostics` | `meds_site_state_types` |
| writing records | `output_serialize` (a site's file set and a region's) | `meds_output_manager` |
| checkpoints | `state_write_state`, `io_read_state` | `meds_io` |
