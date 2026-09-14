# MEDS v0.2.0 release plan

**Status: planned, 2026-09-13. No code written. All decisions taken — see §10.** This document
groups every open issue into phased pull requests for the v0.2.0 release, and records the decisions
taken to get there.

**What v0.2.0 is.** v0.1.0 shipped a model that runs: the coupled fast/slow loop, a conservation
ledger that closes, 45 tests on two compilers, and a Python driver. It shipped explicitly
pre-benchmark. v0.2.0 is the release where **the numbers become defensible and the configuration
surface becomes real** — the silent-wrongness class is closed, every deliberate number-mover lands
in one rebaseline, and the strategies and forcing modes the model claims to support can actually be
selected.

**What v0.2.0 is not.** It is not the benchmark release. `psi_leaf` still does not converge at the
production cadence (#162), no site has been scored against observations, and the ED2 comparison
(#114) is a description of behaviour, not a validation.

Baseline at planning time: `main` at `2ee6d5f`, 45/45 green on ifx.

---

## 1. Scope

43 of the 66 open issues have work in v0.2.0; 23 are deferred to v0.3+. The deferral line is
**major enhancements** — new subsystems rather than completions of existing ones. §9 lists them.

**42 of the 43 close. #1 does not** — it carries a stale-reference fix in Phase 0 but is held open
deliberately as a standing design question; see §2.1.

Six working phases plus release. Phases are ordered by *risk and dependency*, not by issue number:

| Phase | Theme | Issues | Number-moving? |
|---|---|---|---|
| [0](#2-phase-0--clear-the-board) | Clear the board | 13 | No |
| [1](#3-phase-1--silent-wrongness) | Silent wrongness | 6 | Bit-identical at default, or new signal only |
| [2](#4-phase-2--structure-and-performance) | Structure and performance | 3 | Byte-identical, verified |
| [3](#5-phase-3--the-rebaseline-window) | The rebaseline window | 6 | **Yes — one golden re-cut** |
| [4](#6-phase-4--configurability-unlocks) | Configurability unlocks | 7 | Only on newly selectable paths |
| [5](#7-phase-5--diagnostics-and-evaluation) | Diagnostics and evaluation | 6 | No |
| [6](#8-phase-6--release) | Release | 2 | No |

A phase is one or more pull requests. Where a phase splits, the split is stated.

---

## 2. Phase 0 — clear the board

**SHIPPED 2026-09-13** (PRs #203, #204, #205, #206, #208). Twelve of thirteen closed; #1 stays
open by design, and the check it was waiting on came back live and is filed as #207.

No behaviour change. Decisions, stale records, dead code. This phase exists so the later phases are
not reviewed against a board full of noise.

**#191 goes first, before anything else in this phase.** It touches 17 test files, and every later
phase adds tests.

| # | Work | Size |
|---|---|---|
| #191 | Consolidate the per-test `check` / `check_close` into `meds_test_support`. **Measured: 13 files define their own `check(`, 4 their own `check_close(`, against 28 that already use the shared module** — the issue says 14, which understates it. | S |
| #194 | Correct the GPU claims in `CMakeLists.txt` comments and the project docs to match the measurement: the offload build runs 1.4x slower than CPU (49.8 s against 36.2 s). | S |
| #199 | Delete the `--parity` preset from `scripts/numerics_sweep.py` and its use in `parity_scenarios.py`. Confirmed dead: it pins `fast.leaf_energy_model` and `fast.wood_energy_model`, neither of which exists. | S |
| #168 | **Premise is stale — see §12.** `bsap` is no longer a placeholder. Re-measure the wood time constant against the real allometric `bsap`, delete the stale comment, correct the ROADMAP line, close. | S |
| #166 | **Premise is stale — see §12.** `veg_energy_step_implicit` was deleted in PR #120. Delete its three surviving references: the module header at `meds_plant_biophysics.f90:8`, the test comment at `test_column_ark.f90:200`, and the row in `docs/science/column_biophysics.md:130` — which also names `veg_energy_diagnostic`, equally gone. The real public kernel is `veg_energy_balance`. | S |
| #200 | `meds_capi_demography` calls `meds_allometry` instead of importing `b1Ht`/`b2Ht`/`agb_c1`/`agb_c2`/`lai_b1`/`lai_b2` and inlining. Same for the two hard-coded coefficients in `examples/example_demography/empirical_laws.py`. | S |
| #1 | **Repoint the stale file references only. DO NOT CLOSE — see §2.1.** `meds_test_vital_rates.f90` no longer exists and the co-dominant sweep is now `update_overtopping_lai` in `meds_demography_update.f90`. | S |
| #192 | **Decided (§10.1):** the year-rollover `seam[soil_carbon_rh]` residual (8.37e-4 kgC/m2) is attributed to the year-boundary environmental-scalar reset, documented, and no longer reported as an unexplained residual. | S |
| #193 | **Decided (§10.1):** `test/` stays flat. Record why — the suite is 47 files with unique names, and a mirror costs a rename per source move for no lookup gain. | S |
| #198 | **Decided (§10.1): delete the IMEX-Euler tier.** It is not an oracle — it returns no reference trajectory. `imex_euler_column_step` is a two-line wrapper around ARK's own `column_be_stage` plus `advance_water_mass_full`, i.e. the production kernel in a degraded configuration. Port its three tests as below, then remove `imex_euler_column_step` and `adaptive_imex_march`. `meds_fast_rk4_oracle` then holds only `rk4_column_step` and its name is accurate again. | M |
| #201 | Write the frozen-seam contract note: the Lambda criterion, the four-seam classification, the arbitration rule. Marked optional in the structure plan, but two real defect classes came out of not having it. | M |
| #7 | The nvfortran array-valued-function trap is already a standing rule in `CLAUDE.md`. Close with a pointer. | S |
| #6 | **Decided (§10.1):** the RT bugs are **already reported upstream** — no further work there. Fold the MEDS-to-ED2 mapping into `docs/ed2_comparison.md` and close. | S |

### 2.2 #198 — porting the three IMEX-Euler tests

`imex_euler_column_step` is `column_be_stage` + `advance_water_mass_full`, and `column_be_stage`
already takes the `niter` argument (`1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton`).
Nothing in the tier is independent of ARK, so none of its tests needs it:

**Shipped 2026-09-13.** The plan said three tests; there were **five** consumers, and none was
dropped — all five kept their coverage:

| consumer in `test_column_derivs.f90` | outcome |
|---|---|
| `test_be_euler` — L-stability at 900 s, and agreement with the RK4 oracle at dt = 4 s | kept, on a **test-local** `be_euler_step` |
| `test_be_coupled` — `niter=1` collapses / `niter=12` survives | kept; the uncoupled baseline still leaves the band at 415 K |
| `test_arrowhead` — Newton from a 99 %-saturated start, and vs the RK4 oracle | kept |
| `test_adaptive_march` — the controller responds to `rtol` | **ported to `adaptive_ark_march`**, the PRODUCTION controller. Better than the plan's "drop": it had been testing a controller that existed only to serve the tier. Measured 8 steps at rtol 1e-3, 25 at 1e-6. |
| `test_ark2` (a2) — "ARK2 beats the first-order step at the same dt" | kept; it needs a first-order comparand, which `march_be_euler` supplies |

The tier left `src/` and the two-line wrapper moved into the test that degrades it. `be_euler_step`
is one `column_be_stage` plus `advance_water_mass_full` — exactly what the deleted routine was.

**The unplanned win:** deleting the tier let `meds_fast_rk4_oracle` drop `use meds_fast_be_stage`
entirely. The RK4 oracle's independence from the implicit machinery is now **structural and visible
in its import list**, where before it was incidental — the module imported `column_be_stage` for the
tier's benefit while the oracle itself never touched it.

**The remaining accuracy oracle is `rk4_column_step`,** which is genuinely independent of the
implicit machinery — its body is `column_derivs` plus `state_axpy`/`state_accum`, with no
`column_be_stage`, no `newton_surface_solve`, no `advance_water_mass_full`. It also carries
`freeze_theta`, which makes it solve the **same reduced system ARK does** (soil water split out) —
the reason production `rk45` cannot substitute for it.

**Open follow-up, not scheduled:** RK4 is a weak accuracy oracle — 4th order, stable only below
~2.785*tau_fast with tau_fast ~ 17 s. `adaptive_rk45_march` already takes `frozen` as an argument, so a
tightened call inside one frozen record would be strictly better (5th order, error-controlled). It
needs a `freeze_theta` equivalent first. **Note the trap:** running `rk45` at a small `dt_fast`
is NOT this — it refines the freeze along with the step (§2 of `docs/science/numerical_scheme.md`),
so it answers a different question and is not a reference for the frozen system ARK solves at 900 s.
An oracle must hold the same frozen record and sub-step inside it.

---

### 2.1 #1 stays open, by design

**#1 is not a board-clearing item and v0.2.0 does not close it.** It is held open deliberately, as
a standing reminder to revisit the *philosophy and assumptions of cohort splitting under high and
low LAI* — not as a bug report about a shading sweep. Phase 0 fixes only its stale file references,
because a standing question that points at a deleted file is a standing question nobody can follow.

**The scope has grown since it was filed, and this is worth recording.** #1 was written when light
competition was a single Beer-Lambert sweep in the test rate evaluator, and it notes that a
mechanistic radiation layer "would supersede it". That layer now exists — and it does not supersede
it, it **inherits** the problem in a harder place:

- `sort_cohorts` orders cohorts **height descending, DBH descending on ties**, and `higher()`
  returns true on an exact tie, so the insertion sort is **stable**: exactly-tied cohorts keep their
  creation order, which is PFT-index order.
- `update_overtopping_lai` then treats equal-height cohorts as **co-dominant** — they share the
  overtopping LAI and do not shade each other. That is #1's fix, and it holds.
- **`canopy_radiation` does not.** Its own header says cohorts are "ordered BOTTOM (1) to TOP
  (ncoh)", and it builds one discrete RT layer per cohort with its own effective area index. Tied
  cohorts therefore become **stacked layers**, and the lower-indexed PFT is systematically placed
  above the others.

So the co-dominant fix protects the light-plastic *trait* path, while the path that actually drives
**photosynthesis** reintroduces the order- and PFT-dependent bias #1 was filed about. Under high LAI
the stacking penalty on the lower layer is largest, which is precisely the regime the standing
question is about.

**What is confirmed and what is not.** Confirmed by reading: the stable tie-break, and that
`canopy_radiation` stacks cohorts as ordered layers. **Not confirmed:** that recruits actually tie
*exactly* in both height and DBH in a running model. #1's body asserts they are born at a shared
`min_reproduction_height` with the same DBH, which would make the tie exact and the bias real, but
that path was not re-read. **That one check decides whether this is a live defect or a latent one**,
and it should be run before anything is built on it.

**Recommendation:** file the two-stream ordering bias as its own issue once that check is done, and
schedule it in [Phase 1](#3-phase-1--silent-wrongness) if it is live — it is a textbook member of
that phase's class, a systematic PFT-dependent bias with no signal. #1 itself stays open past
v0.2.0 as the design question it was written to be.

---

## 3. Phase 1 — silent wrongness

**SHIPPED 2026-09-13 except #104** (PRs #209, #210, #211, #212, #213). Six of seven closed.

**#104 is DEFERRED to v0.3.0** — it is an open physics decision the plan should not have
pre-empted. See §3.1 and §9.

Every item here produces a wrong or unmeasured number **with no signal**. Each is bit-identical on
the default configuration, or adds signal only.

| # | Work | Size |
|---|---|---|
| #117 | Scale the CO2 compensation point by O2: `gstar_ppm = ... * (p%o2_mol_frac / o2_ref_gstar)` with `o2_ref_gstar = 0.209` a named constant. **Exactly 1.0 at the shipped default, so the default path is bit-identical.** Note the O2 reference in `docs/science/leaf_gas_exchange.md` and the `gstar25` comment. | S |
| #148 | Export the tissue-water floor's clamped mass from `advance_water_mass_full` into the existing `budget%clamp_mass` channel, reduced to `site%work_clamp_mass`. Both callers (`ark2_column_step`, the RK4 oracle) already thread clamp counters. **Do not** build the per-cohort identity check the issue warns against — it fires on correct corrector behaviour (measured 4.5e-2 kg/plant over a July). | M |
| #104 | **DEFERRED to v0.3.0** (decision, 2026-09-13). E1b was already shipped with #105. E1 is an open physics decision, not a numerics fix — see §3.1. The detector ships and reports the condition, so it is visible while unfixed. | — |
| #160 | Warn from `validate_config` when `time_integrator = "rk45"` runs at the production `dt_fast`. The PR #91 transpiration corrector lives in `advance_water_mass_full`, which RK45 does not call, so RK45 carries a psi_leaf error the default path does not. Porting the corrector to RK45 is the larger alternative and is **not** in v0.2.0. | S |
| #170 | Write `PD_DISTURB_AREA` (slot 29) from the disturbance step and add the registry row. Confirmed: the slot is declared in `meds_site_diag_types.f90` and appears in no `use` list in `meds_output_registry.f90`. | S |
| #185 | Five sub-items, each decided independently: `sw_input_kind`, `timestep_seconds`, `avg_convention`, `elevation(grid)`, and the `SWPART_SIB` / `METAVG_INSTANT` / `METAVG_CENTER` codes that parse but do not route. Default disposition: **validate against the config and stop on mismatch**; where that is not meaningful, stop writing the attribute. | M |
| #153 | **Decided (§10):** delete `heterotrophic_respiration_damm`, `heterotrophic_respiration_flux`, the `hr_model` selector and their tests from mainline. **Preserve them on a branch** (`archive/damm-hr`) pushed before the deletion PR, so the DAMM implementation is recoverable for future work without carrying maintenance cost in mainline. The deletion PR must name that branch. | M |

---

### 3.1 #104 — why E1 was not built

**E1b is already done.** `work_hydro_nsub_site`, `work_nonconv_site` and `work_hydro_thrash_site`
are all registry variables; they shipped with the #105 detector. The issue's "surface the condition"
half needs nothing.

**E1 as specified is wrong when the soil is wet, and that is the open question.** The plan wrote it
down as "on the floor, uptake and transpiration are zero and psi stays put, in one closed-form
sub-step". Reading the solver, that is only half right:

- A plant whose store is on the floor genuinely cannot **transpire** — there is no water to lose.
- But it can absolutely take up water, and today it does. The floored store diagnoses
  `psi_wood ≈ −1.5×10⁴ MPa` against a soil at perhaps −2 MPa, so the root gradient is enormous, the
  uptake term is huge, and the cohort **recovers** — violently, over ~136 adaptive sub-steps, which
  is precisely the 13× wall-clock signature the issue measured.

So E1 as written does not just make the collapsed case cheap; it makes a collapsed cohort
**permanently dead**, because arresting uptake removes the only path back. That is a real ecological
consequence — a cohort that transiently desiccates during a dry spell would never re-hydrate after
rain — and it is not what the issue's cost argument was about.

`docs/ROADMAP.md` §4 already flags this as an **open question** with three candidate treatments:
clamp the potential, arrest the solve, or kill the cohort. The plan picked "arrest" without
engaging with the consequence. The three differ in what they claim:

1. **Clamp the potential** to a physical floor (the flaccid tail, not −10⁴ MPa). Keeps recovery,
   removes the enormous artificial gradient, and therefore removes most of the cost. Changes numbers
   in the collapsed regime only.
2. **Arrest the solve** (E1 as written). Cheapest, and wrong on a wet soil.
3. **Kill the cohort.** Defensible ecology — a fully desiccated plant is dead — and it removes the
   state rather than modelling it, but it is a mortality pathway, not a numerics fix.

**Recommendation was (1)** — the −10⁴ MPa potential is an artefact of `rwc_floor = 1e-4`, not a
pressure any tissue reaches, and clamping it is the only option that leaves the physics recoverable.

**Decided 2026-09-13: deferred to v0.3.0.** It is a physics choice with demographic consequences,
not a numerics cleanup, and v0.2.0 does not need to make it. The condition is not invisible in the
meantime: `work_hydro_thrash_site` reports it, and a nonzero value is both a cost flag and a
correctness flag — the states that cause it have pinned stores whose answer is not to be trusted.

## 4. Phase 2 — structure and performance

**COMPLETE 2026-09-13** (PRs #215–#219). Three closed with measurements (#163, #172, #161); four
deferred to v0.3.0 (#195, #188, #164, #190) — every one of them because measuring its premise
changed the answer.

Byte-identical, verified at 1/2/4/8 threads on both back ends.

**Why before the physics phases.** Phases 3-5 each need long Ithaca verification runs, and #195
alone is about 24% of fast-loop self time. Landing the structural shape first also means the
physics work in Phases 4-5 is written against the final types rather than migrated afterwards.
**The tradeoff** is that an invasive fast-loop refactor lands before the physics changes rather
than after; if Phase 2 destabilises, Phases 3-5 slip behind it.

| # | Work | Size |
|---|---|---|
| #188 | **DEFERRED to v0.3.0 with #195** — it was that item's prerequisite and has little standing value alone. Measured: `column_params_t` is 7 208 bytes of fixed-size arrays memcpy'd per step, and `frozen%params` is written *only* by a 5-line wholesale copy from `col_config`, so the two cannot diverge. Removing it means threading `col_config` into `column_be_stage`, `column_derivs` and the state ops — hot-path surgery in three modules — or holding a pointer, against the style rule. | — |
| #195 | **DEFERRED to v0.3.0** (decision, 2026-09-13). The "24 % of fast-loop self time" is the one premise in this plan that could not be confirmed independently: there is no `perf` or `gprof` on the dev box and `ltrace` is orders of magnitude too slow on a real run, so it can only be judged by implement-and-measure. The fix wants a persistent `column_frozen_t` in the per-thread pool (it is a local rebuilt every call), which is a real refactor to spend on an unverified number. **Baseline is measured and recorded for whoever picks it up: 9.457 s**, min-of-5, 1-year Ithaca ARK, `OMP_NUM_THREADS=1`. | — |
| #190 | **DEFERRED to v0.3.0, PAIRED WITH #146** — four of its five benefits have already landed piecemeal, and the fifth is #146's hazard class. See §4.2. | — |
| #164 | **DEFERRED to v0.3.0** — its two supports both gave way on inspection. See §4.1. | — |
| #166 | **Moved to Phase 0 — premise is stale (§12).** `veg_energy_step_implicit` was already deleted in PR #120. |
| #172 | **CLOSED — the premise does not hold.** The FAST tier is *already* on the general machinery: it shares the registry, the buffers, `close_tier` and the serializer. Only extraction differs, and it must — sub-daily values are sampled *during* the fast loop, and by fold time `site` holds the end-of-slow-step snapshot. `SRC_S_*` (4000–4999) and `SRC_F_*` (5000–5999) are disjoint, so neither switchboard can resolve the other's ids. Deleting the staging would delete sub-daily sampling. Recorded at the site. | S |
| #161 | **CLOSED — measured, no work to save.** The rescue fires **zero times** over a full simulated year at Ithaca on `rk45` at 900 s (`work_rk45_rescue_site = 0` against 70 577 sub-steps on the same run, so the counter is live). Building it would also mean one `dt_fast`'s ledger summing two schemes' flux contributions, and it is only well-defined for the `stiff_bail` trigger — `rk45_state_railed` tests the *final* state, so there is no identified sub-step to resume from. Recorded at the site. | S |
| #163 | **DONE — deleted.** Measured first, as the item demanded, and the measurement settled it: `soil_energy_step_implicit` hard-codes `flux%nsub = 1` and BE is unconditionally stable, so substepping could only buy accuracy — which the outer march already owns through `GRP_SE`. The dead surface was wider than filed: `rtol` fed only `GRP_SOIL_T`, **a tolerance group with no member in `state_wrms_grouped`**. | S |

---

### 4.2 #190 — deferred, and paired with #146

`column_cohort_t` (`meds_fast_types.f90`, 22 arrays) is a **read-only per-patch view**, rebuilt each
step by `copy_column_cohort` from one patch's CSR section of `site%cohort`. It is not state; the
per-thread scratch exists only because the gather writes into it.

**Four of the five benefits the design claimed have already landed**, in the PRs that came before:

| design benefit | status |
|---|---|
| Hand-built test views are allometrically inconsistent and never set `bwood`, so the wood heat capacity runs on uninitialized memory | **Done.** `column_cohort_init` builds through `init_cohort` → `set_cohort_size` and then runs the *production* gather; all three column tests use it. The comment at `test_column_ark.f90:62` records the bug and its fix. |
| The three hard-coded constants become PFT parameters | **Done.** `p_leaf_width`, `p_branch_diameter`, `p_crown_area_frac`. |
| Derived geometry becomes cohort-block fields | **Mostly done.** `sapwood_carbon`, `sapwood_area`, `wood_area` are cohort-block fields. Only `lai`/`wai` are formed in the gather, as `nplant × area` — a deliberate unit conversion, commented "never stored". |
| Seed and clamp leave the fast gather | **Done.** `reconcile_tissue_water_capacity` runs once per slow step. |
| One per-field policy table | **Partly — see below.** |

**On the policy.** The issue says a field's fusion or scaling kind "is decided at each call site".
That is no longer true: it is centralised in `fuse_cohort_fast_state` and
`scale_cohort_ground_fields` in `state/site`, with the kinds as explicit headings (INTENSIVE by
heat capacity / EXTENSIVE by nplant / GROUND, add-never-weight). The first calls itself "the one
place the per-field policy is declared". PR #119's film-water fix produced them.

**What a declared table would still add is exactly one thing: completeness you cannot forget.**
Adding a per-cohort fast field and omitting it from `fuse_cohort_fast_state` compiles clean and is
silently wrong at the next fusion. A table the blend *iterates* cannot omit a field that is in the
type; two hand-written routines can.

**That is #146's hazard class**, and #146 was deferred with the measurement "there are currently
zero live defects of this class" — because Fortran has no reflection and only a packed, iterated
layout makes omission a compile-time error. Spending a 38-reference hot-path refactor to buy that
property for ~14 cohort fields, when the same property for `column_state_t`'s 1 207 references was
judged not worth buying, is the wrong order to do them in.

**Paired with #146** (decision, 2026-09-13): do them together or not at all. Together, one packed,
policy-carrying layout serves both the fast state vector and the cohort slice, and the table pays
for itself across both. Separately, each is a large refactor buying a fraction of one property.

### 4.1 #164 — why the bare-array conversion was deferred

The item rests on two supports and neither holds up.

**1. The stated precedent does not transfer.** `solve_plant_water_batch` and
`leaf_gas_exchange_batch` are bare-array because they have a **cohort axis** to batch over — that
is what makes a kernel a vectorisation or offload candidate. The four remaining kernels have no
such axis: `cas_column_step_implicit` is one canopy-air box (already three scalars in, three out),
`soil_energy_step_implicit` one soil column, `soil_carbon_step` one patch's pools, the snow kernels
one pack. Flattening their derived types would produce ~25-argument signatures with nothing to
batch, against a tree rule that says readability beats succinctness when performance is equal.

**2. The justification was measured as not paying.** Bare arrays buy OpenMP `target` eligibility,
and #194 established that offload runs **1.4× slower** than the CPU. The *other* property the
convention protects — standalone-buildability, rule 1 — is **already satisfied**: no kernel folder
contains a `site_t` use (the two greps that hit are both comments).

So the conversion would cost readability across four kernels to buy a property that is already held
by other means, plus one whose payoff is negative. **Deferred rather than closed**: if #196
(cohort-axis threading and vectorisation) ever lands, the vectorisation argument returns — but it
returns for kernels with a cohort axis, which these are not.

## 5. Phase 3 — the rebaseline window

**Every deliberate number-mover in v0.2.0 lands here, in one window, with one golden re-cut and one
before/after CHANGELOG entry per item.** Splitting them across the release would mean re-cutting
goldens repeatedly and losing the attribution of each change.

Run order within the phase: #145 first (it changes the soil temperature every other item is
evaluated against), then #118, then the three small ones, then #47.

| # | Work | Size |
|---|---|---|
| #145 | **Decided (§10): Dirichlet temperature anchor.** Add a thermal-BC selector following the pattern `free_drain \| bedrock \| aquifer` already sets for water, with a Dirichlet branch anchoring the base at mean annual surface temperature. Confirmed: `ENERGY_BC_GEOTHERMAL` is the only bottom-BC code and `frozen%hydrology%geothermal` is hard-wired to zero. Measured today: fitted e-folding depth 5.01 m at a 2 m column and 3.98 m at 3 m, against a physical 2.0-2.5 m; the 3 m base layer still swings 13.1 K. Because soil temperature drives the CENTURY scalar and decomposition is exponential in it, this biases Rh seasonally **and**, by Jensen, in the annual mean. Do **not** rely on the 2 m-vs-3 m warm bias measurement — both runs restarted from a 2 m state file. | L |
| #118 | **Decided (§10): give C3 its own co-limitation curvatures.** Add `theta_cj_c3` (~0.98) and `theta_ip_c3` (~0.95), mirroring the `theta_cj_c4` / `theta_ic_c4` that C4 already has, and stop passing `p%theta_j` into `combine_limits`. Measured: at 0.85 the smoothing costs **29% of assimilation** against `min(Ac, Aj, Ap)`, nearly independent of Vcmax (30/31/31/32% at Vcmax25 = 60/90/120/150) — a systematic offset, not a regime effect. Expect **~20% higher C3 GPP at ambient CO2**. **This PR must also re-examine the PFT Vcmax presets**: if they were implicitly tuned against the current smoothing, the compensation is currently hidden inside a curvature named for a different process, and they need re-setting here. | L |
| #152 | Unify `solar_cosz` (Cooper 1969, `23.45*sin(2pi(284+doy)/365)`) and `daylength` (White 1997, `-23.44*cos(2pi(doy+9)/365)`) on one `solar_declination(doy)`. | S |
| #167 | Add the 1.25*h free-convection slope to the tissue energy linearization. Understates the conductance response at low wind today. | M |
| #89 | Saturation vapour pressure. MEDS uses a bare Bolton form (`611.2*exp(17.67*tc/(tc+243.5))`) with **no ice branch**, which matters for sublimation and the frozen canopy. Evaluate the ED2 proposal (EDmodel/ED2#442) and adopt or document. | M |
| #47 | **Decided (§10.1): keep the Sabot two-limb scheme and document the divergence** from ED2's Manzoni-style `stoma_beta` deliberately — the capacity-limb shape and targets, the stomatal-limb scope and driver, and the `stoma_beta = -sref_stomata * lambda_psi_exp` parameter split, all in `docs/science/leaf_gas_exchange.md`. Fix the related defect the issue names — the driver leaves `psi_soil` at 0 even though `soil_psi_root` is computed, just after the leaf call. Fixing the ordering is a number-mover and belongs in this window. | M |

---

## 6. Phase 4 — configurability unlocks

Things the model claims to support and cannot select. Two pull requests.

### 6.1 PR A — the daily-reduction machinery, then its two consumers

**#150 and #176 need the same machinery**: a fast-loop daily reduction feeding a running-mean
cohort structure-of-arrays column, with the lockstep reorder, every creation site, and the fusion
blend. Build it once. This is why #176 is in v0.2.0 at all — on its own it would be a phase of its
own; behind #150's machinery it is a consumer.

| # | Work | Size |
|---|---|---|
| #150 | Phenology P3. Four cohort columns (`pheno_water_avg`, `pheno_low_psi_days`, `pheno_high_psi_days`, `pheno_light_avg`) — today they are locals, so they are re-zeroed daily. Thread the real drivers: a soil-water running mean, a daily **maximum** leaf water potential, a running-mean radiation, and the shallow soil-layer temperature in place of the air-temperature proxy (`meds_vegetation_dynamics` hard-codes `avail_water = 0`, `dmax_leaf_psi = 0`, `rad = 0`). Lift the `validate_config` rejection per cue as its driver lands. **Acceptance:** the tropical drought-deciduous and light-driven leaf-exchanging strategies run from configuration alone. **Decided (§10.1): flush keeps the present fixed high `k_flush_max` for the light-driven strategy** — it is what `test_plant_phenology` pins today (`'leaf-exch: flush = k_flush_max (permissive)'`), so matching flush to shed is a deliberate behaviour change, not a wiring detail. Revisit only if the emergent leaf-area cycle is implausible once the cues run. **Ships with a stated validation caveat — see §6.1.1.** | L |

#### 6.1.1 The phenology validation caveat

**Decided (§10.1): wire the strategies, and say plainly that the phenology is not validated.**

Be precise about what is and is not covered, because "untested" would be wrong:

- **The kernel is unit-tested.** `test_plant_phenology.f90` already exercises **all four**
  strategies — evergreen, temperate-deciduous, drought-deciduous and light-exchanging — by feeding
  the cue values in directly, plus the rate mapping and the [0,1] drive bounds.
- **The driver is unit-tested for two of the four.** `test_phenology_driver.f90` covers the
  temperature-driven pair end to end, plus birth state and the no-temperature no-op.
- **What #150 adds is untested by construction:** the cue *drivers* themselves — the soil-water
  running mean, the daily-maximum leaf water potential, the running-mean radiation, and the
  shallow-layer soil temperature. Nothing has checked that the values the fast loop computes are
  the values the kernel was written against.
- **Nothing at all has validated the emergent behaviour.** No MEDS run's leaf-area cycle has been
  scored against a phenology observation, at any site, under any strategy — including the two
  temperature strategies that have shipped since v0.1.0.

**What v0.2.0 therefore claims:** the four strategies are *selectable and self-consistent*, not
*correct*. The two new drivers each need a unit test against a hand-computed cue value, and the
release notes must carry the caveat alongside #162. Do not describe the drought-deciduous or
light-driven strategies as working in the CHANGELOG or the ED2 comparison — describe them as
runnable.
| #176 | Thermal acclimation: a running-mean tissue temperature shifting the peaked-Arrhenius reference, on #150's machinery. The `t_acclim` seams were shaped in and then dropped from the environment records, so they need re-adding. | M |
| #151 | Phenology P4: `retained_carbon_fraction` on leaf shed, with the full-removal closure. Optional companion `root_phen_factor`. Today every gram of shed leaf carbon goes to litter. | M |
| #177 | Charge maintenance respiration on the storage pool. Today a large reserve is thermodynamically free to hold. | S |
| #178 | Layered fine-root maintenance respiration, replacing the single root-fraction-weighted mean soil temperature. Pairs with #145 — both are about the model using the per-layer soil temperature it already resolves. | M |

### 6.2 PR B — forcing and traits

| # | Work | Size |
|---|---|---|
| #182 | Implement the Brutsaert/Idso clear-sky LWdown synthesis and lift the `validate_config` rejection. Until this lands, a forcing source without a longwave field cannot drive MEDS at all. | M |
| #179 | Per-PFT hydraulic traits — conductance, vulnerability shape, pressure-volume curve. Today wood density is the only axis on which PFTs differ hydraulically. The constitutive curves are already in a shared library, which is what makes this cheap. | M |

---

## 7. Phase 5 — diagnostics and evaluation

What makes the release legible to someone who is not its author.

| # | Work | Size |
|---|---|---|
| #175 | An evaluation notebook against the Ithaca test bed and a PFT / size-class plotter in `post_proc/`. **The highest-value item in this phase** — the PFT and size-class output axes currently have no reference consumer. Depends on Phase 3 being complete, or it evaluates numbers that are about to change. | M |
| #169 | Mortality carbon by pathway: background, cull, disturbance. The slow-loop ledger now values mortality per phase, which is the seam this needs. | M |
| #171 | Per-band albedo and up-welling shortwave and longwave. The two-stream solves them every step and the output discards them, so a run cannot be compared against a radiometer or a satellite product. Needs a surface radiative-flux record. | M |
| #189 | Per-layer face budget imbalance on the committed path, plus per-cohort tissue residuals, plus RK45 ledgers asserted after the rail decision. This is the check that catches the **vertical-only** defect class that whole-column ledgers are blind to. | M |
| #174 | A companion output slot per variable for `AGG_MEANSQ`. The aggregation exists and nothing reads it. | M |
| #173 | **Decided (§10.1):** rename `[io]` to `[state]` — the block now carries only `output_dir`, `output_prefix`, `write_state`, `state_interval_years`. **With a deprecation path**: accept `[io]`, warn, document the rename. A 0.x minor is the cheapest moment it will ever have. | M |

---

## 8. Phase 6 — release

| # | Work |
|---|---|
| #114 | Refresh `docs/ed2_comparison.md` to v0.2.0. It currently pins MEDS at tag `v0.1.0`. Everything Phase 3 moved has to be restated. **Do not re-introduce the retracted non-stomatal-limb GPP-vs-dt_fast numbers.** |
| #162 | Not fixed. Record `psi_leaf` non-convergence at `dt_fast = 900 s` as a **stated known limitation** in the release notes: the error is inherited from the canopy air and amplified about 4x, and the residual relocates to `psi_wood` through the frozen uptake seam. Every other state and flux converges. |
| #150 | Carry the **phenology validation caveat** (§6.1.1) into the release notes beside #162: the four strategies are selectable and self-consistent, not validated. No leaf-area cycle has been scored against an observation at any site. |
| — | `CHANGELOG.md` release section; prune the shipped items out of `docs/ROADMAP.md`; version bump; tag. |
| — | Give each shipped dev plan a tombstone and move it to `archive/` in the PR that closes its last item, per the `dev_plans/README.md` rule. `MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md` and `MEDS_GPU_EVALUATION.md` both reach zero open items in this release. |

---

## 9. Deferred to v0.3+

23 issues. The line is **new subsystems, not completions**.

| # | Title | Why deferred |
|---|---|---|
| #190 | Delete `column_cohort_t` for typed slices with a policy table | **Targeted at v0.3.0, paired with #146.** Four of five benefits already landed; the fifth is #146's silent-omission class. See §4.2. |
| #164 | Bare-array forms for four more kernels | **Targeted at v0.3.0.** Both supports gave way: the `_batch` precedent is about a cohort axis these kernels lack, and the device-eligibility payoff was measured negative (#194). Standalone-buildability is already held. See §4.1. |
| #195 | Allocator traffic in `build_column_frozen` | **Targeted at v0.3.0.** The 24 % figure is unverified on this box; baseline 9.457 s is recorded in §4 for an implement-and-measure attempt. |
| #188 | `column_params_t` copied into the frozen record every step | **Targeted at v0.3.0 with #195**, whose prerequisite it was. 7 208 bytes/step, and the copy cannot diverge from its source. |
| #104 | Plant hydraulics burns 13× wall clock on a collapsed store | **Targeted at v0.3.0** by decision, 2026-09-13. E1b already shipped (#105). E1 is a physics choice — clamp the artefact potential, arrest the solve, or kill the cohort — with demographic consequences, and arresting (the plan's original pick) makes a collapsed cohort permanently dead. See §3.1. |
| #157 | Fire | New subsystem. Named by the user as explicitly deferred. |
| #155 | Vertically resolved soil-carbon pools | Named by the user as explicitly deferred. |
| #156 | Coarse woody debris pool | Pairs with #155. |
| #154 | The nitrogen twin | New subsystem: `f_decomp`, mineralization/immobilization, NPP limitation, restart serialization. |
| #184 | A transient or observed CO2 stream | **Targeted at v0.3.0** by decision, 2026-09-13. v0.2.0 keeps a constant. This costs little: `[forcing].co2_const` is already a config key (default 420.0, documented in `meds_config_main.toml`), so a run can already choose a different **fixed** CO2 — only the time-varying stream waits. |
| #183 | The multi-polygon runtime | Large and orthogonal — a state hierarchy, an array of readers, a polygon loop and MPI. |
| #186 #187 | Snow P1 multi-layer, P2 canopy interception | New subsystem each. |
| #180 | Per-layer root nodes (hydraulics Phase B) | Follows #179's per-PFT traits, not v0.2.0. |
| #181 | Hydraulic redistribution | Blocked on deciding how the ledger treats water moving between layers through the plant. |
| #146 | Packed `column_state_t` | 1,207 field references; **zero live defects of this class**. Cost already measured. **Paired with #190** — one packed, policy-carrying layout should serve the fast state vector and the cohort slice together. |
| #196 | Cohort-axis threading and vectorisation | Large; follows Phase 2's structural work. |
| #197 | Single-precision experiment | An experiment, not a fix. |
| #158 | Adaptive freeze cadence | Efficiency, and the last item in the numerics plan. |
| #159 | Soil water in the ARK tableau | Measured cost +14-26%. |
| #165 | Canopy film thermal state and phase change | New state with its own energy. |
| #74 | Condensate onto leaf/wood surface water | **Depends on #165** — and only bites once `canopy_water_on` is enabled, which is not the default. Soil layer 1 is the correct fallback today. Carries a coverage gap worth closing with it: no unit test supersaturates the CAS. |
| #96 | Kelvin `e_i` for leaf transpiration | Built, measured, removed, and documented. Its real value is foliar water uptake, which is its own piece of work. |

---

## 10. Decisions

### 10.1 Taken, 2026-09-13

1. **Scope:** all six phases, 49 issues.
2. **#118:** give C3 its own co-limitation curvatures, and re-examine the Vcmax presets in the same
   pull request.
3. **#153:** delete the unreachable DAMM and Q10 kernels and their tests from mainline; preserve
   them on a branch so the implementation is recoverable.
4. **#145:** Dirichlet temperature anchor with a thermal-BC selector, not a deeper column.
5. **#47:** keep the Sabot two-limb scheme and document the divergence from ED2 deliberately. Fix
   the `psi_soil` ordering defect regardless — it is a number-mover and belongs in Phase 3.
6. **#150:** wire the strategies, and ship a stated caveat that the phenology is not validated
   (§6.1.1). Flush keeps the present fixed high `k_flush_max` for the light-driven strategy.
7. **#173:** rename `[io]` to `[state]`, with a deprecation path.
8. **#192:** attribute the year-rollover seam residual to the environmental-scalar reset and
   document it.
9. **#193:** `test/` stays flat; record why.
10. **#198:** **delete** the IMEX-Euler tier. It returns no reference trajectory, so it is not an
    oracle; it is a wrapper around ARK's own `column_be_stage`. Port its tests (§2.2) and remove it.
11. **#184:** keep a constant atmospheric CO2 in v0.2.0; the transient/observed stream is
    **targeted at v0.3.0**. `[forcing].co2_const` already lets a run pick a different fixed value.
12. **#6:** the RT bugs are already reported to ED2. Fold the mapping into `ed2_comparison.md`
    and close; no upstream work.

**No open decisions remain.** Every question this plan raised has been answered; anything that
surfaces during implementation goes to a new issue rather than back into this section.

---

## 11. Verification discipline

Per phase, in addition to the standing rules in `CLAUDE.md`:

- **Phases 0, 1, 2:** 45/45 (rising) on ifx **and** nvfortran. Phase 2 additionally requires
  byte-identical output at 1/2/4/8 threads, which is the standard PR #109 established.
- **Phase 1:** #117 must be demonstrated bit-identical at `o2_mol_frac = 0.209`. #104 and #148 are
  measurements before they are fixes — report `hydro_nsub` and `clamp_mass` on a run that
  desiccates a cohort, so the "unobserved" claims become numbers.
- **Phase 3:** one golden re-cut at the end of the phase, and a before/after number in
  `CHANGELOG.md` for **each** item, not for the phase.
- **Phase 4:** acceptance for #150 is that the two currently unreachable strategies run from
  configuration alone. Each of the four new cue drivers additionally needs a unit test against a
  hand-computed value — the kernel is already covered for all four strategies, the drivers are not
  (§6.1.1). Selectable is not validated, and the release notes say so.
- **A green suite proves nothing** if the changed code is outside every build path or the
  assertions are regime-blind. Phase 3's items in particular need a forced run, not a unit test.
- **A closed budget proves bookkeeping, not plausibility.** #148 exists precisely because the
  whole-column water ledger closed to 4e-12 kg/m2 while the floor was creating water.

---

## 12. Premises re-measured, 2026-09-13

Filed issues were checked against the source before being scheduled. All held except one.

**#168 is stale — the issue and the ROADMAP line are both wrong.** `bsap` is no longer a
placeholder. `set_cohort_wood_geometry` derives it from ED2's real `b1SA`/`b2SA` sapwood-area
allometry via `sapwood_fraction(dbh, ...)`, and `meds_pft_params.f90` states it: "the sapwood
FRACTION of basal area sets bsap. REPLACES bsap = 0.10*wood_carbon." It landed in PR #125. What
survives is one stale comment in `meds_vegetation_dynamics.f90` saying `bsap` "is still an MVP
placeholder". The remaining question — whether the wood thermal time constant is now realistic — is
a measurement against the real value, not a fix.

**#166 is stale.** `veg_energy_step_implicit` was already deleted in PR #120 ("dead code out").
Three references survive it — a module header comment, a test comment, and a row in
`docs/science/column_biophysics.md` that also names `veg_energy_diagnostic`, which does not
exist either. The real public kernel is `veg_energy_balance`. #166 drops from a Phase 2 refactor
to a Phase 0 documentation fix.

**#191 is understated:** 13 files define their own `check(` and 4 their own `check_close(`, not 14
in total.

Confirmed as filed: #1 (behaviour intact, file references stale), #104, #117, #118, #145, #148,
#153, #160, #163, #170, #172, #174, #185, #188, #190, #195, #199, #200.
