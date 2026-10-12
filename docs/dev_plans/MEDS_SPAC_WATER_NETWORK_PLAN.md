# Soil and plant water as one network: process functions and one implicit solver

Written 2026-10-11 on branch `feat/spac-water-solver` (from beta 0cbbb31). Design only: nothing here is
implemented. The owner chose (2026-10-11) to solve soil and plant water together, and asked that the
work consolidate the code rather than add fixes for special cases, with functions designed around
processes and generic numerical schemes. Decisions D1–D8 are in §9.

## 0. Summary

Soil water and plant water are solved today by two different solvers, each treating the other as a
fixed boundary for the whole 900 s fast step. Every recent water failure comes from the seam between
them, and most of the code around it exists to patch that seam. This plan replaces both solvers, and
the patches, with:

1. **Process functions**, one per physical process, each returning its value and its derivatives:
   water storage in a soil layer or a plant tissue, flow along a path (soil to soil, soil to root, wood
   to leaf, pond to soil), and the boundary fluxes (transpiration, soil evaporation, rain, runoff,
   drainage).
2. **One generic implicit solver** for a network of water stores joined by flow paths. It knows nodes
   and edges, not soil or plants: implicit Euler, Newton with a line search, a structured linear solve,
   adaptive sub-steps.
3. **One call per patch per fast step**, which hands the energy integrator the water fluxes it needs.

The plant solver, the soil solver, the uptake rescale, the residual floors and caps, the ARK plant
re-solve and the RK45 in-stage water all go (§5).

## 1. Why

### 1.1 The failures

- **BCI regrowth NaN (2026-10-11).** In the step ending 1994-04-20 (ERA5 2014), a drought-depleted
  stand refilling its wood asked for 4.3 kg/m² in 15 minutes from a 2.6 cm topsoil layer that rain had
  just wetted. The soil's cap let the roots take every drop above θr, about 3.1 kg/m². The next rain
  step's linearised Richards solve then met a layer at ψ ≈ −6e8 m: no convergence in 2000 sub-steps,
  face flows of 4e17, a θr floor that created 4e17 kg/m² of water, and the energy solve overflowed.
  The reproducer replays it in one simulated day (`meds_runs/ex06_regrow/meds_config_bci33d.toml`).
- **#352 (BCI 2017-04-17 storm front).** Same precondition: a top layer near residual meeting heavy
  rain. #354 fixed how the solver evaluates flows after its solve, not why the layer was that dry.
- **#104 and example 06's April collapse.** When the soil cannot supply the plant's request, a single
  patch-wide factor scales down every cohort's uptake but not its sapflow, so the wood drains toward its
  floor. In example 06 (Palo Verde), the nearly bare Ds and D cohorts' leaf potential reaches −25 to −100
  MPa every April. This predates #386–#388; the c2 and a_g0 design runs show it too.
- **The transpiration mismatch.** The plant re-solve in the ARK stages uses stage transpiration but
  the uptake priced in the pre-pass, so the plant ends about 1.3 % short of water at 900 s
  (`meds_fast_be_stage.f90` ~456–470).

### 1.2 The common cause

- **The plant solve** (`solve_plant_water`, 2 nodes per cohort) sees each soil layer as a reservoir of
  fixed ψ and fixed rhizosphere conductance for the whole step (`meds_fast_frozen.f90` 375–423).
- **The soil solve** (`advance_soil_water_column`) sees root uptake as a fixed demand, priced at the
  start of the step, capped at θr (`meds_soil_water.f90` 150–165). It is a single linearised implicit
  pass with C and K at θⁿ (`linearize = "frozen"`) and no lower bound on ψ.

Neither side sees the other's capacity, and the wet-soil exchange between them has a time scale of
about 0.4 s (`MEDS_PRODUCTION_INTEGRATOR_PLAN.md` P1), so no explicit hand-off can be consistent at
900 s. The caps, floors, rescale, re-solve and mass floors are each a patch on that seam.

## 2. Principles

1. **One process, one function.** Each physical process is one pure function that returns its flux or
   storage and their derivatives. No solver keeps its own copy of a process.
2. **The solver is generic.** It sees nodes with storage curves, edges with flow functions, and
   boundary terms. Nothing in it is specific to soil or plants, so new nodes (a root node, a stem
   segment, a litter layer) or edges (a root conductance in series) need no solver change.
3. **Delete the patches, do not port them.** A guard that exists only because two solvers disagree
   goes away with the disagreement. A guard that encodes physics (the air-dry limit) becomes part of
   the constitutive curve, not a post-hoc clamp.
4. **Conservation by construction.** Water mass closes to round-off because the residual is a mass
   balance and every flux is used once, on both sides. There is no residual floor, clip or rescale to
   book.
5. **One time scale for the water.** One implicit solve per patch per fast step, with its own error
   control; the energy integrator receives time-averaged water fluxes, as it does today.

## 3. What exists today, and its fate

| Component | Where | Fate |
|---|---|---|
| Plant 2-node solver (frozen-coefficient matrix exponential, step doubling) | `meds_plant_hydraulics.f90` `solve_plant_water`, `freeze_coeffs`, `advance_exact_linear`, `exact_substep` | Replaced by the network solver |
| Batch wrapper, root-boundary aggregation, per-layer split by end-of-step supply | `solve_plant_water_batch`, `effective_root_boundary` | Replaced: per-layer uptake is an edge flux |
| Declared but unused 3-node topology, `HYDRO_SOLVER_BE` | `meds_plant_types.f90` 169–182 | Deleted; a root node becomes a network option |
| Rhizosphere conductance, root profile, PV curves, Kirchhoff/PLC | `rhizosphere_cond`, `meds_hydr_lib.f90`, `meds_water_retention.f90` | Kept, moved into the process library with derivatives |
| Soil Richards solve (linearised BE / Celia Picard, step doubling) | `meds_soil_water.f90` `soil_water_advance`, `soil_water_step_implicit`, `face_and_sink` | Replaced by the network solver |
| Explicit Richards RHS for RK45 | `soil_water_time_deriv` | Deleted |
| Root-sink cap at θr; evaporation cap at θr | `meds_soil_water.f90` 150–165, 683–684 | Deleted (uptake and evaporation are implicit) |
| Saturation clip to the pond; θr floor, give-back, hard floor | `meds_soil_water.f90` 206–244 | Deleted (pond is a node; curves are bounded) |
| Infiltration capacity formula | `meds_soil_water.f90` 139–148 | Replaced by a pond–soil Darcy edge |
| Ground evaporation (Philip α + DSL) | `ground_evaporation` | Kept as a boundary-flux process, with its ψ derivative |
| Bottom boundary (free drain, aquifer, bedrock); pond overflow | `bottom_flux`, `pond_overflow` | Kept as boundary processes |
| Pre-pass plant pricing, `root_share`, uptake rescale | `meds_fast_frozen.f90` 375–534 | Replaced by one network call |
| ARK plant re-solve (transpiration corrector) and its mass floors | `meds_fast_be_stage.f90` `advance_water_mass_full` (127 lines) | Deleted |
| ARK stages moving θ at a constant rate | `meds_fast_be_stage.f90` 173 | Kept; the rate comes from the network |
| RK45 θ and plant mass in the tableau, stage clamps, commit-side clip/floor and pond rebuild | `meds_fast_rk45.f90` (≈150 lines), `meds_fast_time_derivs.f90` 399–403 | Deleted; RK45 uses the network like ARK (D3) |
| Water-enthalpy hand-off: `qloss_frozen`, `qwflux_wl`, `q_wood_net`, `water_store_enth`, clip/floor enthalpies | `meds_fast_frozen.f90` 562–621 | Replaced by one generic upwind advection over network edges |
| Work counters `hydro_nsub`, `hydro_thrash` (and its L2 stop), `soil_nsub`, `n_fail` | `meds_fast_frozen.f90` 424–437, 519, 617 | Replaced by the network solver's own counters |
| Config keys `[soil]` `linearize`, `max_picard`, `substep`, `h_init`, `max_substep`; plant solver options | `meds_biophysics_opts.f90`, `meds_plant_types.f90` | Retired with an error message; `[soil]` `rtol`/`atol` become the network's tolerances |

Measured cost of what is replaced (VTune, BCI regrowth stand with 600 cohorts in 25 patches, six
months, 9.1 s on 12 threads): plant hydraulics 11.5 % of CPU, soil water 4.2 %.

## 4. Design

### 4.1 The network

For one patch:

- **Nodes.** One per soil layer, one pond, and for each cohort one wood node and one leaf node (a root
  node later, D7). Each node holds water; its unknown is the water potential ψ [MPa or m of head,
  one convention throughout]; its mass follows from its storage curve. Each node has an elevation.
- **Edges.**
  - between adjacent soil layers (Darcy);
  - pond to the top soil layer (infiltration, Darcy from a saturated node);
  - each soil layer to each cohort's wood (rhizosphere conductance, optionally in series with a root
    radial conductance, #356);
  - each cohort's wood to its leaf (xylem with vulnerability).
- **Boundary terms.** Transpiration out of each leaf node (fixed for the step, from the energy
  pre-pass), soil evaporation out of the top layer (a function of its ψ), rain into the pond, pond
  overflow above its capacity, drainage at the bottom.

The graph has a fixed shape: the soil layers and pond are **shared nodes**, and each cohort is a
**private group** of nodes that connects only to itself and to shared nodes. The solver uses that
shape generically (§4.3), without knowing what the nodes are.

### 4.2 Process functions

Each function is pure and elemental where it can be, takes the states at its ends, and returns the
value and its derivatives with respect to those states. The derivatives are analytic; each is unit-tested
against a finite difference.

| Process | Function returns | Today's source |
|---|---|---|
| Soil water storage | θ(ψ), C = dθ/dψ; van Genuchten or Campbell, bounded below at an air-dry ψ (D5), with a small specific storage above saturation so the curve has no corner | `meds_water_retention` |
| Plant tissue storage | W(ψ), C(ψ) from the pressure–volume curve (symplast + apoplast) | `meds_water_retention` (PV) |
| Pond storage | W = ρ·h, linear | new, trivial |
| Darcy flow between soil layers | q, ∂q/∂ψ_a, ∂q/∂ψ_b; face K the thickness-weighted geometric mean, plus gravity | `face_and_sink` |
| Infiltration | the same Darcy edge, from the pond (saturated) to layer 1 | replaces the capacity formula |
| Soil to root | g(θ_k)·(ψ_k + ρg z_k − ψ_wood), with g from Gardner; a root radial conductance in series when configured (D7); a smooth floor at zero while hydraulic redistribution is off (D2) | `rhizosphere_cond` |
| Xylem, wood to leaf | Kirchhoff form with the vulnerability curve, plus the gravity head | `meds_hydr_lib` |
| Soil evaporation | E(ψ₁) with Philip's α and the dry-surface-layer resistance, and dE/dψ₁ | `ground_evaporation` |
| Transpiration | the step's demand per leaf node (fixed, D4) | pre-pass `transp_c` |
| Rain, overflow, drainage | boundary fluxes with their ψ derivatives | `pond_overflow`, `bottom_flux` |

A leafless cohort needs no special case: its leaf node has almost no storage and no transpiration, and
its ψ follows the wood through the xylem edge.

### 4.3 The solver

Generic over any network of nodes, edges and boundary terms:

- **Time stepping.** Implicit Euler, which is L-stable for the 0.4 s exchange mode, on sub-steps inside
  the fast step. Each sub-step has an error estimate from the difference between an explicit predictor
  and the implicit solution (no step doubling, so one solve per sub-step). Sub-step sizes are kept across
  calls for a warm start. A sub-step is rejected and shortened if Newton fails or the error is too large.
- **Nonlinear solve.** Newton on the mixed-form residual: for each node,
  (W(ψ) − Wⁿ)/h − Σ(edge flows in) + Σ(boundary losses) = 0. Mass enters as W(ψ), not C·Δψ, so
  conservation is exact at convergence. An Armijo line search on the residual norm and a bound on the
  ψ change per iteration keep Newton safe in very dry nodes, where C is small.
- **Linear solve.** The Jacobian is block-structured by §4.1's shape. The solver eliminates each
  private group's block (2×2 for a cohort today) onto the shared nodes it touches: a Schur complement
  that adds a low-rank update per group. The reduced system over the shared nodes (10–16 soil layers plus
  the pond) is solved densely by LU, and the groups are back-substituted. Cost per iteration grows as
  groups × shared² + shared³, a few thousand flops for a patch. The routine is written for any
  "groups attached to a shared core" graph, not for soil and cohorts.
- **Results.** End-of-step ψ and mass for every node; each edge's water moved over the step; each
  boundary term's total; Newton, sub-step and rejection counts.

### 4.4 Coupling to the rest of the fast step

1. **Pre-pass** (as today): radiation, aerodynamics, photosynthesis and stomata (reading ψ_leaf at the
   start of the step), interception, the starting `surface_derivs`, which gives each cohort's
   transpiration demand.
2. **Water network** over the fast step, with that demand.
3. **Energy integrator** (ARK, or RK45): canopy air, soil energy and tissue temperatures, with the water
   fluxes held at their time averages over the step, exactly as the frozen fluxes are held today.
   Soil θ moves at the network's rate within the stages.
4. **Water enthalpy, generically.** Every edge carries its water at the temperature of the node it
   leaves (upwind); each node's enthalpy source is the sum over its edges. That one rule replaces the
   per-pathway formulas for root uptake, sapflow, tissue refill, clip, floor, infiltration and drainage
   heat.
5. **Commit.** Water states from the network, energy states from the integrator; the ledgers read the
   network's edge and boundary totals.

## 5. What is deleted

- The plant solver and its batch wrapper, root-boundary aggregation, per-layer supply split, the unused
  3-node and BE declarations, the thrash counter and its L2 stop.
- The soil solver's step-doubling wrapper, linearised and Picard steps, `face_and_sink`, the explicit
  Richards RHS, the root-sink and evaporation caps at θr, the saturation clip, the θr floor, give-back
  and hard floor, and the infiltration-capacity formula.
- The pre-pass plant pricing, `root_share`, the uptake rescale.
- The ARK plant re-solve and its mass floors; RK45's in-stage water, stage clamps and commit-side
  corrections.
- The per-pathway water-enthalpy formulas.
- The retired config keys.

Rough size: about 1,300 lines leave and about 900 arrive, most of them moved process functions, so the
water code shrinks and the number of places that move water falls from seven to one.

## 6. Phases

Each phase is a pull request into beta that leaves the test suite green.

### P0: baselines and the test set (no model change)

- Freeze reference results for the cases in §8, on beta 0cbbb31.
- Unit tests that the new code must pass are written first, against analytic or existing results.

### P1: the process library

- Move each process in §4.2 into one module per kind (storage, edges, boundaries), with analytic
  derivatives; the current solvers call the moved functions, so model results do not change.
- Add the air-dry bound to the soil curves (D5) and the specific storage above saturation; these change
  results only where today's code is outside its valid range.
- Exit: unit tests for values and derivatives; ctest green; example 04 unchanged to round-off.

### P2: the generic network solver

- A standalone module: network description, residual and Jacobian assembly from the process library,
  the grouped Schur linear solve, damped Newton, adaptive implicit Euler.
- Tests on synthetic networks:
  - a chain of linear reservoirs (analytic solution);
  - a single plant against today's matrix-exponential solver in its linear range;
  - a soil column against today's Richards solver in wet conditions;
  - mass conservation to round-off;
  - the BCI 1994-04-20 state, dumped from the model: the solve converges and no node passes its bounds.
- Exit: unit tests green; the solver used nowhere in the model yet.

### P3: one call per patch

- Build the patch network each fast step; replace the plant solve, the soil solve, the rescale, the ARK
  re-solve and RK45's in-stage water with one network call; the generic enthalpy advection; ledgers from
  the network's totals.
- Delete everything in §5 in the same pull request, and retire the config keys.
- Exit: ctest green (tests of deleted code rewritten against the network); §8's reproducers pass with
  no NaN and no water created; water and energy ledgers close; example 04's fluxes against the tower
  reported before and after; the profile in §7 rerun.

### P4: documents and examples

- Rewrite `docs/science/plant_hydraulics.md`, the soil-water part of `soil_biophysics.md` and
  `numerical_scheme.md` around the network; CHANGELOG with before-and-after numbers.
- Rerun examples 04–06 and update their READMEs; recalibrate example 04 if its fit moves (owner's
  call).
- Move this plan to `archive/` with a tombstone; update `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` (N2d, P1)
  and `MEDS_EFFICIENCY_SWEEP_PLAN.md` Phase 3 (#104), which this supersedes.

### Later, each a small change once P3 is in

- A root radial conductance in series (#356): a different soil-to-root edge function.
- Hydraulic redistribution (#181): drop the smooth floor on that edge.
- A root node per cohort (#180): a 3-node private group; no solver change.
- Water and energy together: transpiration as a function of ψ_leaf inside the Newton (D4).

## 7. Efficiency

| | Share of CPU now | Estimate after P3 |
|---|---|---|
| Plant hydraulics (about 3 solves per cohort per step: pre-pass and ARK re-solves) | 11.5 % | |
| Soil water (3 linear solves per step, typically) | 4.2 % | |
| Water network (1–2 sub-steps × 2–4 Newton iterations; each about one of today's plant coefficient evaluations per cohort, plus a 10–16-node dense solve) | | 9–24 % |

- **Ordinary weather:** total run time changes by −7 % to +8 %, if the ARK plant re-solve is dropped
  (keeping a second network pass for transpiration consistency would add about 10 %).
- **Storm fronts and drought refill:** much cheaper than today's worst cases (hundreds of soil sub-steps,
  137 plant sub-steps per cohort in #104), and no failures.
- **Parallelism:** patches stay independent, so the OpenMP loop over patches is unchanged.
- **Separately:** 19–21 % of CPU in that profile is threads waiting on uneven patches, a larger prize
  than this plan, for the efficiency sweep.

## 8. Tests and reproducers

| Case | What it checks | Where |
|---|---|---|
| BCI regrowth, 1994-04-20 | the crash: dry thin layer, refill, first storm | `meds_runs/ex06_regrow/meds_config_bci33d.toml` (restart, one day) |
| BCI 2017-04-17 storm front | #352 | the one-day restart from trial 7fe525 (see #352) |
| Example 06, April 2009–2010 | #104: leaf potential of nearly bare cohorts stays physical | `examples/example06_dry_forest` (draft branch) |
| Example 04, one year | fluxes against the BCI tower; ledgers | `examples/example04_column_biophysics` |
| Palo Verde regrowth, 20 years | biomass trajectory; cost | `meds_runs/ex06_regrow/meds_config_lr33.toml` |
| Six months of the BCI stand | profile against §7 | `meds_runs/ex06_regrow/meds_config_bci33prof.toml`, `run_prof.sbatch` |

## 9. Decisions for the owner

- **D1 Unknown per node.** ψ, with mass from the storage curve (the mixed form). *Recommended:* conserves
  mass exactly and handles saturation through the curve.
- **D2 Hydraulic redistribution.** Keep it off at first, through a smooth floor at zero on the
  soil-to-root edge. *Recommended;* turning it on later means removing the floor.
- **D3 RK45.** Use the same network and drop its in-stage water. *Recommended,* so both integrators
  treat water identically; RK45 keeps the energy side.
- **D4 Transpiration.** Fixed per fast step from the pre-pass, as now (phase 1). Coupling it to ψ_leaf
  inside the Newton is the later water–energy step.
- **D5 Air-dry bound on the soil curves.** ψ ≥ −1e5 m (CLM's minimum soil potential), or −1e4 m (the
  air-dry head the evaporation code already uses).
- **D6 The pond as a node.** *Recommended:* it replaces the infiltration-capacity formula, the clip and
  the separate pond bookkeeping with one Darcy edge and a storage.
- **D7 Root radial conductance (#356).** Leave it out of P3 (behaviour as today, rhizosphere only), or
  add it as an optional series term in the same edge function, off by default.
- **D8 Tolerances.** The network uses `[soil]` `rtol`/`atol` (already the ARK θ tolerances) for soil
  nodes, and new plant-node tolerances in MPa; or one set for all nodes.

## 10. Related

- Issues: #352 (storm front), #104 (wood store drained), #356 (root conductance in series), #375 (dry
  season refill), #377 (absorbing root length), #180 (root nodes), #181 (hydraulic redistribution).
- Plans: `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` (N2d, the reorder so soil water sees the realised uptake,
  and P1, the in-stage root sink that was reverted); `MEDS_EFFICIENCY_SWEEP_PLAN.md` Phase 3 (#104);
  `MEDS_FROZEN_SEAM_CONTRACT.md` ("borrowing one solve's flux while committing another's state").
