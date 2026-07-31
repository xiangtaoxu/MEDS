# MEDS production numerical scheme — the plan

> **SUPERSEDED IN PART — 2026-07-31.** Three things in this document no longer hold:
> 1. **The `split` integrator is RETIRED** and `[fast].integration_scheme` is deleted. `ark` is the
>    default; `rk45` is the accuracy baseline and its stiff rescue now redoes the step on `ark`.
>    Anything here that treats split as the default, the reference, or a comparison anchor is history.
> 2. **"IMEX-ARK" is a misnomer.** The biotic CO₂ source is folded implicit, so `f_E == 0` and the
>    scheme is a 2-solve **ESDIRK2** (γ = 1 − 1/√2). The config string stays `"ark"`.
> 3. **`dt_fast` is STABILITY-limited, not accuracy-limited.** Above ~150–225 s the frozen surface
>    coupling drives a sustained period-2 canopy-air oscillation (~8 K at 900 s) that **no
>    conservation ledger detects**. Default is now 150 s. Every measurement in this document taken at
>    900 s or 1800 s is inside that regime.
>
> Current state: `docs/science/numerical_scheme.md` §2–3 and
> `docs/dev_plans/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–14.

**Status:** design, 2026-07-29. Supersedes the `select_integrator(cfg)` sketch in
`MEDS_NUMERICS_SCOPING.md` goal (c), on the evidence of Phase D (§3d) and Phase E (§3e) of
`MEDS_INTEGRATOR_PARITY.md`.

**Question this answers:** what does MEDS run by default, what may a user change, and what does the
model do on its own when a step turns out to be hard?

---

## 1. What the evidence actually constrains

Four measured facts do most of the work. All from `b4_stand_summer` / `b3_stand_winter`, one month,
`--parity`, hourly FAST scoring.

**F1 — `dt_fast` is the freeze cadence, and it is the dominant error term at production settings.**
At `dt_fast = 1800 s` ARK's *time-stepping* error is 0.039 K while its *total* error against a
resolved reference is 1.10 K. Nearly all of it is the Category-0 coefficient freeze (photosynthesis,
conductances, radiation, hydraulics, and — on the ARK — the entire soil-water solve). Tightening the
integrator tolerance at fixed `dt_fast` cannot touch it.

**F2 — refining `dt_fast` works, and how well depends on the scheme.** 8× refinement buys 17× for
ARK and RK45 (1.10 → 0.062 K) and 2.6× for `split`, which then stalls. *(The earlier claim that
nothing converges was a reference artefact; see §3e E-0.)*

**F3 — `split` and ARK/RK45 converge to different answers**, ~0.45 K apart in canopy-air temperature
and ~0.15 K in soil-surface temperature. They are two models, not two accuracies of one model. The
cause is unattributed; the previous difference of this kind (missing sapflow advected enthalpy) was
real physics.

**F4 — ARK is on the efficiency frontier once accuracy is a constraint.** `ark @ 900 s` (7286
sub-steps, 0.42 K) beats `split @ 225 s` (11904 sub-steps, 0.46 K) on both axes at once. Below
`dt_fast ≈ 900 s`, `split` cannot reach ARK's accuracy at any cost.

**The design consequence of F3 is the one that shapes everything below.** A run that switches
integrator *family* mid-simulation switches model mid-simulation, and puts a ~0.45 K step in its own
state trajectory at the switch. That is not adaptivity; it is a discontinuity that no error
controller sees, because the controller measures error *within* a family. So:

> **Adapt the step size within a run. Choose the scheme at configuration time.**

This is also the retrospective verdict on the existing RK45→`split` hybrid rescue: it is a
mid-run family switch, and it is precisely why `work_rk45_rescue_site` had to be invented and why
"read that counter before reading any RK45 result" is a standing instruction. It should be replaced,
not extended (§4, P2).

---

## 2. Target architecture

Three layers, each with one job. The point of the separation is that each layer has a *different*
correctness criterion, so conflating them (as `select_integrator(cfg)` would have) makes the failure
modes inseparable.

```
  layer 3   COUPLING CADENCE      dt_fast          set by the science target
            (what is frozen, how often)             -> physics fidelity (F1)

  layer 2   INTEGRATOR FAMILY     time_integrator  set once per run, never switched
            (which model you solve)                 -> ark | split | rk45   (F3)

  layer 1   STEP CONTROL          adaptive march   fully automatic, within-family only
            (how well you solve it)                 -> rtol/atol + I/PI controller (F2, F4)
```

### Layer 3 — `dt_fast`, the coupling cadence

Chosen by the user from the science target, documented rather than auto-selected. It is *not* an
accuracy knob that can be left alone: it sets how often photosynthesis, radiation and hydraulics see
a refreshed column.

| target | `dt_fast` | expected CAS-T error (ARK) |
|---|---|---|
| multi-century spin-up, demography-dominated | 1800 s | ~1.1 K |
| standard production, coupled carbon–water | 900 s | ~0.4 K |
| flux-tower comparison, diel partitioning | 450 s | ~0.14 K |
| reference / verification | ≤ 225 s | ~0.06 K |

### Layer 2 — the integrator family

**`ark` becomes the production default**, once the blockers in §3 clear. Rationale, in order of
weight:

1. It converges (F2), so `dt_fast` is a working accuracy lever.
2. It is L-stable on the stiff block, so a hard step *slows down* rather than *falls over*. It has
   never needed a rescue path, at any `dt_fast`, in any Phase B/D/E cell.
3. It is on the efficiency frontier wherever accuracy matters (F4).
4. Its error control acts on the states that actually couple (canopy air + soil heat), which is where
   the stiffness is.

`split` stays available and supported as the **cheap tier**, but its documentation changes: it is a
*different model*, ~0.45 K away, not a coarser solution of the same one — until F3 is attributed.
`rk45` stays as the **ED2-comparison and reference tier**.

### Layer 1 — step control

Already built (`meds_fast_control`): per-group WRMS tolerances, I and PI controllers, warm start, a
sub-step floor of `dt_fast/64`, and L0/L1/L2 strictness. What changes is that it becomes the *only*
thing allowed to react to a hard step, and that its norm is made scheme-consistent (§3, P0-c).

---

## 3. The plan

### P0 — correctness blockers. **All five LANDED 2026-07-29.**

**P0-a. Route snowfall on the ARK/RK45 path when `snow_on = false`. ✅ DONE.** (§3e E-4.)
`build_column_frozen` now carries `forc%precip + forc%snowf`, as `meds_fast_split.f90:384` always
has. ARK/RK45 month-total Δθ went +0.009602/+0.009611 → **+0.013325/+0.013328** against split's
unchanged +0.013318 — from 39% apart to 0.08%. Covered by `test_column_dynamics` **RUN 9**, which
asserts a *boundary-input identity* (same day run twice, snowfall on and off, same sub-freezing air;
the column's liquid store must gain what fell) rather than a ledger residual — see P0-e for why.
Mutation-checked: reverting the one line gives a gain of exactly 0.000000 kg/m² on both schemes.

*Settled since* (§3f, F-1): **`[fast].snow_on` is gone — the snow store is always active.** The
switch was over a structural fact, and its default was the dangerous value: `precip_phase` splits
snowfall out of the met precipitation without consulting it, so `off` never meant "no snow", it meant
"snow with nowhere to go". Always-on is free on a snow-free column (`snow_accumulate` returns unless
a pack exists or the snowfall clears `min_new_snow_mass`), and the winter month now reproduces the
old `snow_on = true` control exactly, all three schemes within 3.2%.

**P0-b. Close the `snow_st%exists` ordering hazard. ✅ DONE.** (§3e E-6.) The interception sweep moved
below `advance_snow_stage`. Nothing between the old and new positions reads `f_wet_c` or
`intercept_leaf/wood`, so the move is otherwise inert.

**P0-c. Make the adaptive error norm scheme-consistent. ✅ DONE.** (§3e E-5, §3f F-2.)

- **θ absent from the norm — FIXED.** `state_wrms_grouped` gains `with_theta`, default `.false.`
  (byte-identical for the ARK, whose θ stage-difference is structurally zero); the RK45 march passes
  `.true.` via the existing `GRP_THETA`. Measured on `b4_stand_summer`: split and ARK byte-identical
  (the control), RK45 accuracy unchanged (0.0609 → 0.0609 K at `dt = 225`) and cost slightly lower
  (9381 → 8987 sub-steps at 450). **Honestly: nearly free rather than beneficial in this cell**, where
  the temperature terms dominate the norm. Its value is in the wetting/saturation regimes this cell
  does not exercise, where θ moves fast and previously moved uncontrolled.
- **The `with_mass` / `with_theta` switches are now GONE entirely** (§3f, F-2). One norm over the
  whole column state, no per-caller opt-outs. The ARK's θ and mass terms are structurally zero and
  therefore dilute its norm — accepted deliberately, because a switch reading "do not measure this
  state" is indistinguishable at the call site from "this state cannot move", and the two *were*
  confused. Measured cost: ARK accuracy ≤0.5%, sub-steps **−9 to −15%**. The clean way to retire the
  dilution is P1's soil-water-into-the-tableau work, not a flag.

**P0-d. Guard the reference-cell timing in the harness. ✅ DONE.** (`numerics_sweep.py`: the
reference `dt_fast` must divide both `dt_slow` and 3600 s.) Still to add: the corollary check in
`parity_fidelity.py` — **refuse to score** when a cell's FAST record count differs from the expected
hours, or when its `(day, hour)` keys are not unique. The failure this prevents was silent for a
whole phase of work and is generic to any tiered-output comparison.

**P0-e. Arm the budget hard-stop. ✅ DONE, and it immediately found something.** (§3e E-7.)
`[energy].debug_error` had **no TOML reader** — it existed only as a `.false.` default settable from
Fortran, so `budget_check_stop` (wired by QW2 specifically to stop being dead code) was dead code in
every configured run, and every "verified by a forced run with `debug_error`" claim in the project
record is vacuous *with respect to halting*. Now wired in `load_energy_opts`, default `.false.`.

With it armed, a forced January Ithaca month **halted on all three schemes** at
`whole_energy resid = -3.51e3 J/m², tol = 1.45e3` — pre-existing and shared, so it never contaminated
scheme comparisons. ~0.4% of a winter day's net radiation: exactly the size a per-step tolerance sees
and a seasonal diagnostic does not.

**RESOLVED (2026-07-29) — this P1 blocker is CLEAR.** The cause was the `veg_coupling_floor` clamp
destroying energy in the diagnostic leaf/wood balance (fixed in PR #81). Reverting that one fix on
current `main` still reproduces the halt (`resid = -1.51e3, tol = 1.42e3`), which attributes it; the
original −3.51e3 was that leak plus contributions since closed by #82–#85, so it is the dominant cause
rather than provably the only one. Because the clamp lives in `surface_derivs` — shared by all three
schemes — the failure was identical on all three, and that symmetry is what identifies it as one shared
kernel defect rather than a deep seam.

A January month now closes with ~5x margin on every scheme:

| scheme | worst abs(resid) (J/m²) | tolerance (J/m²) | margin |
|---|---|---|---|
| `split` | 269 | 1394 | 5.2x |
| `ark` | 285 | 1395 | 4.9x |
| `rk45` | 287 | 1395 | 4.9x |

Verified the halt is genuinely **armed** in that configuration by tightening its tolerance until it
fires, rather than inferring closure from a clean exit — the same trap that made `debug_error`'s missing
TOML reader invisible for so long.

**Consequence for the default:** winter closure no longer blocks `ark`. What remains before the default
moves is **P1-a** — attributing the residual `split` <-> ARK/RK45 family gap (~0.45 K), which is a
scientific question rather than a bookkeeping one.

### P1 — make the recommendation defensible

**P1-a. Attribute F3 — the residual `split` ↔ ARK/RK45 gap.** This is the scientific question the
production choice hangs on: if `split` is the one that is wrong, the default *must* move; if ARK is,
the whole convergence argument points at a different limit than the physics does. Method, reusing the
one that found the sapflow term: single-`dt_fast` step from an identical state, instrumented per-term
diff of the leaf and ground energy balances, then refinement at fixed total time to separate the
step-size part from the floor. Known-unequal candidates already visible from code reading:
`root_sink_share` (split-only, but `multilayer_roots` is pinned off under `--parity`, so not this),
the condensation quadrature (exponential relaxation vs per-stage rate — a *numerics* difference that
would shrink, and it does not), and the ARK's frozen soil-water solve versus split's re-solve inside
the Picard iterate. That last one is the leading candidate and is testable directly with the
`in_tableau` mask.

**P1-a′. The three remaining Class-1 physics rows are separately planned.** Per-layer root-sink
placement, prognostic leaf/wood energy and the non-free-drain bottom boundary all differ across
schemes, but **all three are inert at default settings** and are pinned by `--parity`, so none of them
confounds P1-a. Their closure plan is
[MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md](MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md) (design, 2026-07-30).
It is *feature coverage* work, not attribution work — do not expect it to move F3.

**P1-b. Generalise the evidence base.** Everything above is one site, one month per cell, one forcing
year, no confidence intervals. Before the default changes: at least a dry-down window (where the
transpiration↔uptake seam and the hydraulics freeze bite), a cold season with snow on, and a wet
saturated window (the only regime where the clamps and the pond composition matter). One month each,
all three schemes, `--parity`, both compilers.

**P1-c. Cut ARK's rejection storm.** 1253 rejections against 5492 accepted sub-steps at
`dt_fast = 1800 s` — 19% of the work discarded — falling to 3 rejections at 225 s. Warm start already
removed the cold-start component. What remains is a coarse-`dt_fast` phenomenon and the cheapest
remaining performance win. Candidates, in order of expected value: (i) the PI controller, which is
built and off by default and which exists precisely to damp step oscillation; (ii) clamping the
controller's growth factor `fmax` below 5 for the first sub-step of a `dt_fast`; (iii) checking
whether the `BETA = 2.414` stage-3 extrapolation is what is being rejected — its clamp counter
(`work_clamp_stage_site`) already distinguishes this and is being written.

### P2 — retire the cross-family fallback

Replace the RK45 → `split` hybrid rescue with a **within-family** degradation:

1. On `stiff_bail` or `rk45_state_railed`, re-run the `dt_fast` on the **ARK**, not on `split`. ARK is
   L-stable, so it resolves exactly the regime RK45 rails in (cold, dense canopy), *and* it is in the
   same semi-discretisation family — so the fallback does not move the model.
2. Keep `work_rk45_rescue_site` as the provenance counter. Its meaning improves: it becomes "this
   step was solved by a different *stepper*", not "this step was solved by a different *model*".
3. Once P0-c and P1-c land, re-measure whether the rescue ever fires at all. `rk45_rescue` was 0 in
   every Phase B/D/E cell except `b3_stand_winter` (0.497 patch-steps per month). If it goes to zero
   with a corrected error norm, the whole path can be deleted rather than re-pointed.

### P3 — the selection rule, once there is something to select

Only build this after P0 and P1-a. Deliberately *not* `select_integrator(cfg)` returning a different
family per condition — see §1. What it should be instead:

- **`select_dt_fast(target)`** — a documented mapping from science target to coupling cadence (the
  §2 layer-3 table), emitted in the run header so every output file records which tier it was run at.
- **A configuration-time validator**, not a runtime switch: given `time_integrator`, `dt_fast` and the
  active process set, refuse or warn on combinations the evidence says are unsound. Concretely:
  `split` with `dt_fast < 900 s` (paying for accuracy it structurally cannot deliver — F3/F4);
  `rk45` with a non-free-drain bottom boundary (already an error stop); any adaptive scheme
  with `[output].numerics = false` in a run intended for a numerics study.
- **A cost ceiling inside the march** — if a `dt_fast` needs more than *N* sub-steps, report it
  through the existing work counters and (under L2) halt, rather than silently spending. The
  machinery exists (`RK45_WORK_CAP`, the `dt_fast/64` floor); what is missing is that it is currently
  wired to a family switch rather than to a diagnostic.

---

## 4. Acceptance criteria

The default moves to `ark` when **all** of these hold:

1. ~~P0-a, P0-b, P0-c landed; `ctest` green on ifx Release, ifx Debug and nvfortran multicore.~~
   **✅ MET** — all of P0 landed, plus the two option removals in §3f. 37/37 on ifx Release, ifx Debug
   and nvfortran multicore.
2. ~~A cold-season month with `snowf > 0` closes the whole-column water ledger at `n_fail == 0` on
   all three schemes, and their soil-water trajectories agree within a few percent.~~
   **✅ MET for water** — `test_column_dynamics` RUN 9, plus the forced-month agreement to
   0.08%. **NOT met for energy**: with the hard-stop now armed (P0-e), the same forced winter month
   fails whole-column *energy* on all three schemes. Pre-existing and shared, but it must be
   attributed before a winter-validated default.
3. The P1-b scenario set (dry-down, cold, saturated, plus the existing four) shows ARK converging
   under `dt_fast` refinement in every cell, with `rk45_rescue = 0` and `clamp_mass = 0` on every
   reference used.
4. F3 is attributed, or — if it turns out to be genuine physics on the `split` side — `split` is
   fixed and the two families agree after refinement.
5. A 30-yr Ithaca run on `ark` completes, conserves, and its AGB trajectory is explicable against the
   `split` run given (1)–(4). The existing `runs/ithaca_ark30/` results predate both the recycle
   phase fix (PR #69) and P0-a, so this has to be re-derived from scratch, not compared to.

Until (1)–(5), the default stays `split` and the recommendation in
[docs/science/numerical_scheme.md](../science/numerical_scheme.md) §6 is what users should follow.

---

## 5. What is deliberately not in this plan

- **A higher-order tableau.** F1 says the freeze dominates at production `dt_fast`; an ARK4(3) over a
  first-order splitting buys nothing. Revisit only if F3 is attributed and `dt_fast` refinement
  becomes the binding constraint.
- **Folding plant hydraulics into the tableau.** Measured 4×–170× *worse* at `dt_fast = 900 s`, and
  it is a category error besides: ψ is advanced by an exact matrix exponential and the coupled
  subsystem is ψ-independent within a step. ED2 makes the same choice for the same reason.
- **Folding soil water into the ARK tableau.** Real, but it buys correctness rather than order, costs
  1.5–2 weeks, and P1-a should establish first whether it is the source of F3 — in which case the
  scoping changes.
- **GPU offload of the fast loop.** Orthogonal. The bare-array kernel conversion that enables it is
  tracked in `MEDS_NUMERICS_SCOPING.md` §11.
