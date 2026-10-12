# Soil and plant water as one network: process functions and one implicit solver

Written 2026-10-11 on branch `feat/spac-water-solver` (from beta 0cbbb31). Parked (owner, 2026-10-11)
until example 06 is finished; later that day the owner asked for the plan to be implemented on this
branch, kept up with beta. An independent review (§11) found the direction sound but the design
incomplete; §12 records how each of its items is resolved and the working choices for D1–D8. The
owner chose (2026-10-11) to solve soil and plant water together, and asked that the work consolidate
the code rather than add fixes for special cases, with functions designed around processes and
generic numerical schemes. Decisions D1–D8 are in §9.

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

## 11. Independent review (2026-10-11) and the changes it calls for

A separate reviewer read this plan against the code. Verdict: sound only with changes. The diagnosis
and the direction (one conservative implicit network, backward Euler with Newton on the mixed form)
stand; the items below must be resolved in a revision before P2.

**Critical**
- **C1 Transpiration.** A transpiration fixed for the step, drawn from a finite leaf store, has no
  solution once E·h exceeds the leaf's storage plus its maximum sapflow (stomata close only to g0, and
  the Kirchhoff supply is bounded): the nearly bare cohorts of example 06 hit this, and
  `MEDS_FROZEN_SEAM_CONTRACT.md` already rules a prescribed flux against a prognostic store inadmissible.
  Deleting the ARK re-solve (which cut the ψ_leaf error at 900 s from 1.454 to 0.0046 MPa) also needs a
  replacement that keeps the canopy air's vapour credit equal to the leaf's debit. Change: supply-limited
  transpiration inside the network, E = E_demand·s(ψ_leaf), smooth, with its derivative; choose the
  consistency mechanism by measurement before P3 (the network's E frozen in the stages, N2d's
  energy-first order, or a second warm-started pass); the dt-convergence tests of
  `numerical_scheme.md` §5a become P3 exit criteria. D4 is reopened.

**Major**
- **M1 Error estimate.** A predictor–corrector difference misbehaves on stiff modes (the leaf after each
  change of E, the 0.4 s wood–soil mode after rain). Use warm-started step doubling, a filtered estimate,
  or TR-BDF2, with a mass-weighted norm.
- **M2 Curves for Newton.** A hard air-dry bound gives C = 0 and, at −1e4 m, leaves Philip's α at 0.49, so
  evaporation never stops: use an oven-dry extension (C > 0, θ → 0 near −1e5 m; Webb 2000). For van
  Genuchten with n < 2, dK/dψ is unbounded at saturation: add the air-entry modification (Ippisch et al.
  2006). Replace every `max()` floor with a smooth form, and bound Newton steps in Se/RWC or log|ψ|, or
  let each store choose its unknown (Se dry, ψ wet). Split P1 into a bit-identical move (P1a) and a
  measured curve change (P1b).
- **M3 Pond.** An empty pond and the overflow threshold are switching points: specify a smooth
  depression-storage curve, the spill, and the face-K rule (it sets Horton runoff and moves example 04).
- **M4 Enthalpy.** One upwind rule is right for soil (total internal energy) but not for the tissue
  stores (cap·T with a frozen capacity, valued at the start-of-step temperature): applying it there
  recreates #355. State the rule per store type with its valuation time, route condensate to the pond,
  and add a #355 refill regression test.
- **M5 Cost.** Line-search and error-estimate evaluations and every cohort×layer edge per iteration push
  the realistic range to −4 % to +20 %, with a 2× tail if M1 or M2 misbehave. Use arrays of edge
  parameters per edge kind (no polymorphic edges) and a merge gate (≤ +15 % on the BCI profile).
- **M6 Size.** P3 touches about 4,800 lines of tests, the RK4 oracle, the ledgers, the IO counters and
  the config template, and the warm-start step is new persisted patch state (restart, fusion, lockstep
  reorder). Split into P3a (ARK on the network, old code still present) and P3b (RK45, deletions, keys).

**Minor:** the diagnosis should also list the per-layer split by end-of-step supply and the ~24 cohorts
pricing one frozen layer; the inventory misses the film rescale, the HR floors and the stage-3 θ clamp;
"seven places to one" overclaims (film, snow, condensate and the slow-loop reconcile remain); commit
masses through the flux form; a C¹ quadratic ramp, not softplus, for the no-HR floor; leafless cohorts
at PLC ≈ 1 make the 2×2 block singular (write the Schur term as g·c/(g+c)); Kirchhoff derivatives need
care as Δψ → 0.

**Tests to add:** restart bit-reproducibility, thread-count invariance, gfortran and nvfortran builds,
Debug `-fpe0` runs of the reproducers, the §5a dt-convergence tests, snow and frozen soil, arid cells,
example 04's calibration objective, a Jacobian-versus-finite-difference check, an infeasible-demand case,
extreme states (nplant 1e-8, leafless, extreme ψ), and bounds on Newton and sub-step statistics. Measure
the #104 premise first (is the April collapse the uptake rescale or a fixed E on a tiny leaf store?), and
prototype the solver offline on the dumped 1994-04-20 state in P0.

**The reviewer's decisions:** D1 mixed form with a well-scaled per-node unknown; D2 off, with a C¹
quadratic ramp; D3 agree, and move the RK4 oracle too; D4 disagree as written (C1); D5 an oven-dry
extension near −1e5 m; D6 agree once M3 is specified; D7 leave #356 out of P3; D8 one mass-based set
(Δθ, ΔRWC).

**Modularity (owner's question, 2026-10-11), to be written into §4 in the revision.** The solver is
generic, so a plant-only run (soil potentials prescribed) and a soil-only run (no plants, or a prescribed
uptake) are two configurations of one assembly, not separate code. Three generic features make that
real: a prescribed-potential node type; boundary terms that are either fixed fluxes or functions of a
node's ψ; and network assembly driven by the active processes, exposed through the C-API so the Python
front end can run a single plant or a single soil column. P2's single-plant and soil-column tests become
those configurations.

## 12. Implementation decisions (2026-10-11), answering §11

The owner asked for the plan to be implemented; these are the working choices, each one reversible
and the measured ones marked. They follow the reviewer's recommendations except where a reason is
given.

### 12.1 D1–D8

- **D1 The unknown.** Each kind of store picks a well-scaled unknown and the solver sees only that
  (Brenner and Cancès 2017):
  - a **soil layer** uses its effective saturation Se on the dry side and its potential on the wet
    side, joined smoothly (C¹) at the point of the curve's largest capacity: Se = (1 + m)^−m for van
    Genuchten, the air-entry point for Campbell. Dry, Se is the variable that moves smoothly; wet, ψ is;
  - a **plant tissue** uses its symplastic relative water content R everywhere: the
    pressure–volume curve has a finite dψ/dR on every branch, above saturation included, and R stays
    defined for a leafless leaf, whose storage is simply zero;
  - the **pond** uses its head.

  Mass comes from the storage curve (the mixed form), so conservation is exact at convergence.
- **D2 No hydraulic redistribution.** The soil-to-root flow is g·r(Δψ), with r a C¹ quadratic ramp
  of width 1e-3 MPa: zero below zero, linear above the width.
- **D3 RK45** uses the network as ARK does (P3b); the RK4 oracle moves with it.
- **D4 Transpiration:** see C1 below.
- **D5 No air-dry bound.** Van Genuchten and Campbell keep their θr asymptote: ψ → −∞ as θ → θr, with
  a capacity above zero at every finite ψ.
  - D1's unknown keeps Se > 0: a Newton step may remove at most 90 % of the saturation a dry layer
    has left.
  - The network's curves carry no conductivity or capacity floor; the floors stay in the diagnostic
    read-offs.
  - Ground evaporation stops by itself, since Philip's α → 0 as ψ → −∞. That is the case the
    reviewer's oven-dry extension was for: it arose only because of the hard bound.
  - The Darcy flow out of a drying layer, √(K_k K_k+1)·Δψ, vanishes as Se → 0 for every van
    Genuchten n > 1 (as Se^1.25), so no flow can empty a layer.
  - The conductivity floor that let a −6e8 m potential drive 1e-3 m s⁻¹ through a face in the
    1994-04-20 crash is gone.
- **D6 The pond is a node**, specified under M3.
- **D7 #356** stays out of P3.
- **D8 One tolerance set**, `[soil]` `rtol`/`atol`, measured in θ for soil layers and the pond (pond
  water over ρ·dz₁) and in R for plant tissues.

### 12.2 The review items

- **C1 Transpiration.** Each cohort's transpiration is E = E_demand·s(R_leaf), solved inside the
  network.
  - E_demand is today's demand: the pre-pass surface solve at the start of the step, with both stomatal
    water-stress controls.
  - s = 1 while the leaf holds more water than R_close, the symplastic water content at the PFT's
    fast-limb closing potential (`pft_wstress_psi_close`, twice the turgor-loss point by default). Below
    R_close, s = R/R_close, with a C¹ corner.
  - Above full stomatal closure the existing controls act unchanged. Below it, a leaf with its stomata
    shut loses its residual water in proportion to the symplast water it has left. So E → 0 as the leaf
    empties while the sapflow stays positive, and every step has a solution.

  The energy stages then hold each cohort's transpiration at the network's time mean, exactly as they
  already hold ground evaporation (`frozen%ground%soil_evap`): the leaf energy balance takes
  `le_slope = 0`, `le_ref = h_evap·E`. The canopy air's vapour credit is the leaf's debit by
  construction, and the ARK re-solve has nothing left to correct.

  What this gives up: within one step, transpiration no longer responds to the canopy air's changing
  humidity; the next step's pre-pass restores that response. To measure in P3a against
  `numerical_scheme.md` §5a. The fallback, if §5a degrades, is a second warm-started network pass.
- **M1 Error estimate:** ½·J⁻²·(W_{n+1} − W_n − h·F(W_n)).
  - J = I − h·∂F/∂W is the implicit-Euler matrix, applied through the last Newton factorisation, so
    the estimate costs two back-substitutions.
  - On a slow mode this is the implicit-Euler local error, h²W''/2. On a stiff mode just after a jump
    it is y₀/(λh), the same as step doubling gives, at one solve per sub-step instead of three.
  - The norm is mass-weighted (D8).
  - P2 tests the estimate against step doubling on a stiff chain.
- **M2 Curves for Newton.** D1 and D5 above, and every `max`/`min` in the network has a smooth form
  (the root ramp, the transpiration corner, the pond curve).
  - Newton steps are bounded per node (fraction-to-boundary on Se and R; |Δτ| ≤ 0.2 per iteration on
    the dry branch) with an Armijo line search on the mass-scaled residual.
  - The Kirchhoff conductance uses its midpoint expansion as Δψ → 0, for the value and both
    derivatives.
  - Above saturation a soil layer stores water through a specific storage, S_s = 1e-4 m⁻¹. A layer
    pushed past saturation therefore raises its potential and sends water up to the pond and down to
    the next layer, which is what replaces the saturation clip.
  - **Measured, P2:** whether Ippisch's air-entry modification is needed. For n < 2, dK/dψ → ∞ at
    saturation; it goes in only if the saturated-column test shows Newton needs it, since it moves
    wet-soil drainage.
  - The P1a/P1b split: P1 moves the processes value-identical (P1a). The curve changes exist only
    inside the network, so they are measured with it in P2–P3 rather than in a separate P1b.
- **M3 The pond.**
  - **Storage:** W = ρ·S(h), with S = 0 for h ≤ 0, h²/(2δ) for 0 < h < δ, and h − δ/2 above, where
    δ = 1 mm. An empty pond holds no water and has no capacity; its head follows the soil below
    through the edge.
  - **Edge to layer 1:** Darcy, with face conductivity √(K_s1·K₁) over the distance to node 1. That
    is today's infiltration-capacity formula, now implicit and two-way: a supersaturated top layer
    exfiltrates.
  - **Rain** enters the pond.
  - **Spill:** after each accepted sub-step, water above `w_pond_max` leaves as runoff at the pond's
    temperature. This is a store at its capacity, as today, not a correction for a solver.
- **M4 Enthalpy, by store.**
  - **Soil layers and pond** hold total internal energy. Water crossing a soil face, the pond face or
    the bottom carries the liquid enthalpy of the node it leaves; the soil-energy solve does this,
    given the network's time-mean flows.
  - **Root uptake** leaves each layer at that layer's own temperature, layer by layer from the edge
    flows. This replaces the root share times one root-weighted temperature.
  - **Tissues** keep #355's rule: heat is cap·T with the step's frozen capacity; the change in water
    mass is booked at the start-of-step tissue temperature (`water_store_enth`); arriving water heats
    the tissue only by its temperature difference.
  - **Condensate** stays in soil layer 1, unchanged; a top layer pushed past saturation exfiltrates to
    the pond the next step.
  - A #355 refill test joins the suite.
- **M5 Cost.**
  - Parameters are stored as arrays per edge kind, with no polymorphic edges.
  - A cohort's Schur update covers only the layers it roots in.
  - Merge gate: no more than +15 % on the BCI profile.
- **M6 Size.**
  - **P3a** puts ARK on the network, with the old solvers still compiled but no longer called from
    it. **P3b** moves RK45 and the oracle, deletes the old code and retires the keys.
  - **No warm start.** Every fast step begins with one sub-step of `dt_fast`, so there is no new
    persisted patch state and restart, fusion and the lockstep reorder are untouched. Revisit only if
    the profile shows repeated rejections.
- **Minor.**
  - A leafless cohort's group stays invertible: the xylem edge keeps a positive conductance floor (a
    domain guard, like today's `k_floor`) and the leaf storage term is never negative.
  - Masses are committed through the flux form.
  - "Seven places to one" becomes "seven to three": the film, snow and condensate keep their own
    bookkeeping, as does the slow-loop reconcile.
- **Tests.** The §11 list, except the gfortran and nvfortran builds where this cluster lacks a
  compiler: nvfortran is not installed here (`CLAUDE.local.md`).
- **Modularity.** The solver takes prescribed-potential nodes and fixed-flux boundary terms; P2's
  plant-only and soil-only tests are those configurations. The C-API exposure goes to `ROADMAP.md`.

### 12.3 Where the code goes

| Module | Holds |
|---|---|
| `src/shared/util/meds_grouped_solve.f90` | the generic linear solve: a dense LU for the shared core and the Schur update from groups attached to it |
| `src/shared/functions/meds_water_storage.f90` | storage curves in each store's unknown (D1), with derivatives |
| `src/shared/functions/meds_water_paths.f90` | flow along a path (Darcy, rhizosphere, xylem) and the smooth ramp, with derivatives |
| `src/fast_dynamics/soil/meds_soil_water.f90` | the soil's boundary processes (evaporation with its ψ derivative, bottom drainage), which need the forcing types |
| `src/fast_dynamics/water/meds_water_network.f90` | one patch's network (nodes, edges by kind, boundary terms), its assembly, Newton and sub-steps; a new folder in `meds_fast_kernels` |
