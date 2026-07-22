# MEDS `src/driver/` reorganization — design & plan

Status: **DESIGN-ONLY (not executed).** Scope: reorganize the driver layer for legibility, retire dead
code, de-duplicate the split/ARK physics, and settle the integration-scheme roadmap. Behaviour-
preserving; the numerical schemes and the golden anchor are held invariant. Slow soil-carbon
biogeochemistry wiring is **deferred** to a later feature PR (see §8).

Grounded in a full read of the 7 driver files + call graph + adversarial evaluation (workflow
`meds-driver-refactor-plan`, 2026-07-21). All `file:line` refs verified against the working tree at
that time.

---

## 0. TL;DR

- **Two production integrators**: operator-split + Picard, and IMEX-ARK2 — cross-validated at
  *diagnostic* leaf. A **self-converged stiff reference** (tolerance-halved implicit) is the ground
  truth; **RK4 stays as a mild-window cross-check only**. ARK4 and the prognostic-leaf arrowhead are
  **deferred** behind their real prerequisites (lift the operator split; build the 3×3 leaf↔CAS solve).
- **Rename toward a level-legible, scheme-symmetric layout.** `*_dynamics` = a *tier orchestrator* that
  `advance_one_step` fans out to; `meds_fast_split` / `meds_fast_ark` = the two *scheme backends*;
  `meds_fast_time_derivs` / `meds_fast_types` = shared *machinery*. This kills the `stepper`/`loop`
  name overload (the word "stepper" survives only in `meds_stepper`, the calendar-tick conductor).
- **Extract the shared driver types to a leaf module FIRST** — it is the prerequisite that breaks the
  `dynamics ↔ ark ↔ derivs` module cycle a naive file split would create.
- **Five phases, lowest-risk first**, each behind green ifx Debug **and** nvfortran multicore + the
  `tc_split(54)=292.450065` anchor to printed digits + a fast-on carbon smoke run.

---

## 1. How MEDS runs today — the tier hierarchy

MEDS is a **two-timescale model**: a FAST (sub-daily) biophysics column and SLOW (daily/monthly/annual)
vegetation + (future) biogeochemistry. The driver is the orchestration layer (`libmeds_aux`); all edges
point one-way into `src/core` (the demographic mechanism layer) through `meds_core_interface` — verified
no `core → driver` back-edge.

```
program meds                                    [meds_main.f90]         calendar march over dt_slow
  └─ DO WHILE time_lt(now, end_time)
     ├─ advance_one_step(...)                   [meds_stepper.f90]      LEVEL 1 — one calendar TICK
     │   ├─ IF fast_biophysics_on:
     │   │     run_fast_biophysics(...)         [meds_fast_loop.f90]    LEVEL 2 — the SUB-DAILY loop
     │   │       └─ column_fast_step(...)       [meds_column_dynamics]  LEVEL 3 — one dt_fast advance (DISPATCH)
     │   │            ├─ (split) operator-split + Picard          [inline in meds_column_dynamics]
     │   │            └─ (ark)   column_fast_step_ark             [meds_column_dynamics]
     │   │                 └─ adaptive_ark_march / ark2_column_step  [meds_ark_stepper.f90]  LEVEL 4 — the ODE integrator
     │   │                      └─ surface_derivs / column_derivs     [meds_column_derivs.f90]  the RHS
     │   └─ vegetation_dynamics(...)            [meds_vegetation_dynamics.f90]   the SLOW veg tier
     ├─ FAST-tier output drain (netCDF serialize lives in main — the DAG wall)
     └─ annual: io_write_snapshot / io_write_state + has_nan guard
```

**The four levels, and why the current names confuse.** "stepper"/"loop"/"step" appear at three
different scales:

| File (today) | Level | "Steps" what | Called by | Calls |
|---|---|---|---|---|
| `meds_stepper` | 1 (model tick) | the whole model one calendar step (fans to fast+slow) | `meds_main` | `run_fast_biophysics`, `vegetation_dynamics` |
| `meds_fast_loop` | 2 (fast orchestration) | drives the sub-daily loop (patches × N sub-steps) | `advance_one_step` | `column_fast_step` |
| `meds_ark_stepper` | 4 (ODE integrator) | advances the column ODE one dt_fast (ARK2/RK4) | `column_fast_step_ark` | `meds_column_derivs` |

Two structural facts make it worse than similar names: (a) `meds_fast_loop` does **not** call
`meds_ark_stepper` directly — the intermediary `column_fast_step` (in `meds_column_dynamics`) dispatches;
(b) `meds_ark_stepper` is **scheme-specific** (only the ARK backend), while the split scheme has *no*
parallel `_stepper` module — its stepping is inline in `meds_column_dynamics`. So the names imply a
symmetry that doesn't exist. §5 fixes both.

**Two load-bearing facts that dominate the refactor risk:**

- **The fast→slow carbon handoff is shared mutable SoA state.** `run_fast_biophysics` resets/fills
  `site%cohort%gpp_accum` etc. (`meds_fast_loop.f90:215-224, 369-372`); `carbon_growth` reads them
  (`meds_vegetation_dynamics.f90:276-279`). Correctness depends on FAST-before-SLOW order
  (`meds_stepper.f90:47-67`) + CSR layout stability across monthly fusion/fission. **The default config
  is fast-off, so this contract is untested by ordinary `ctest`** — a reorder can pass green and break
  carbon mode. A fast-on carbon smoke run must be in per-phase verification.
- **Slow soil-carbon biogeochemistry is NOT wired.** The CENTURY kernel library
  (`src/biogeochemistry/meds_soil_biogeochem.f90`: `soil_carbon_step`, …) has **zero driver callers**;
  there is no per-patch `soil_carbon_t` on `site_t`. The only biogeochem symbol on the live path is
  `heterotrophic_respiration_flux` over a hard-coded `fast_soil_carbon=5.0` pool feeding the *fast*
  CAS-CO2 twin — not slow soil carbon. So `meds_biogeochem_dynamics` is a **new feature**, not a
  reorg of existing code (§8).

---

## 2. Integration schemes — decision & rationale

### 2.1 What exists

| Scheme | Status | Where |
|---|---|---|
| SPLIT + Picard, full prognostic | **REAL, complete, golden-anchored** | `column_fast_step` default (`meds_column_dynamics.f90:309`) |
| IMEX-ARK2, diagnostic leaf | **REAL, production-reachable** (`[fast].time_integrator="ark"`) | `ark2_column_step` + `adaptive_ark_march` (`meds_ark_stepper.f90:453,662`) |
| IMEX-ARK2, prognostic leaf/wood | **does NOT exist** — `error stop` (`meds_column_dynamics.f90:291-294`) | needs the arrowhead |
| ARK4 | **does NOT exist** — no Butcher tableau / stage-vector / variable-`a_ii` solve | — |
| Fixed-step RK4 | exists but **test-only, diagnostic-leaf-only, blows up at production dt** | `rk4_column_step` (`meds_ark_stepper.f90:359`) |
| IMEX-Euler (1st order) | test-only, **superseded** by ARK2 | `imex_euler_column_step` + `adaptive_imex_march` |

### 2.2 What ARK2 actually is (IMEX — with a caveat)

**IMEX = Implicit-Explicit.** An additive Runge-Kutta splits the RHS `dy/dt = f_I(y) + f_E(y)` into a
**stiff part `f_I`** (advanced implicitly, for stability) and a **non-stiff part `f_E`** (advanced
explicitly, cheaply); each gets its own Butcher tableau, sharing stage times.

**The scheme here is ARS(2,2,2)** (Ascher-Ruuth-Spiteri 1997; *2 implicit stages, 2 explicit stages,
order 2*), L-stable, `γ = 1 − 1/√2 ≈ 0.293`, `β = 1 + √2` (`meds_ark_stepper.f90:462-463`). Per step
(`ark2_column_step`): stage 2 solves `Y2 = y_n + γ·dt·f_I(Y2)` (`column_be_stage`); stage 3 extrapolates a
`β`-weighted base (`state_extrap`, with `clamp_cas`/`clamp_theta` overshoot guards) then solves
`Y3 = base + γ·dt·f_I(Y3)`; combine to 2nd order + an **embedded 1st-order estimate** (`state_err_diff`)
that drives the adaptive controller `adaptive_ark_march`. The stiff implicit block is the **CAS twins**
(backward-Euler in the atmospheric exchange; the coupled 2×2 leaf↔CAS Newton `newton_surface_solve` — the
"arrowhead", leaf folded in *diagnostically*) + **soil heat**.

**Is it IMEX? Nominally yes, but today the explicit half is degenerate.** `f_E ≡ 0`: CO2 is folded into
the implicit numerator, and the two naturally-non-stiff subsystems — **soil water and plant hydraulics —
are operator-split OUT entirely** (soil water committed once from a scratch Richards solve; hydraulics
advanced once via exact matrix-exponential over full dt; both excluded from the embedded error,
`meds_ark_stepper.f90:105-106, 447-451`). So in practice ARK2 = an **L-stable ESDIRK-2 on the stiff
CAS+soil-heat block wrapped in a 1st-order operator split**. This is why **the whole-column order is
capped at 1** regardless of tableau order, and why the whole-water budget closes only to split tolerance,
not machine (`meds_column_dynamics.f90:894-908`).

### 2.3 Critical evaluation of the original three-leg plan

- **Leg (a) SPLIT+Picard prognostic — keep.** Real, complete, golden-anchored. Prognostic leaf requires
  `integration_scheme="picard"` (explicit leaf↔CAS split oscillates).
- **Leg (b) "ARK2 or ARK4, diagnostic or prognostic veg" — wrong axis.** (1) ARK4 buys nothing measurable
  until the operator split is lifted (§2.2) — the real axis is *"is the split lifted?"*, not
  *"ARK2 vs ARK4"*. (2) Prognostic-leaf-under-ARK is `error stop`'d; it needs the leaf↔CAS **arrowhead**
  (a bordered 3-variable Newton over (H, q, leaf_temp) with `a_store = cap_leaf/dt`). And the leaf's
  thermal inertia is physically **negligible** at dt_fast~900 s (~0.015 K), so the arrowhead is a
  *quantification/validation* deliverable, not a production necessity.
- **Leg (c) fixed-step explicit RK4 at fine dt as ground truth — numerically unsound.** The loop is
  genuinely stiff (ratio ~6.5e5); RK4 is stable only below ~2.785·τ_fast and *only because it integrates
  the reduced diagnostic-leaf system* (`meds_ark_stepper.f90:6-11`). It is **stability-limited, not
  accuracy-limited** ("very fine dt" makes it stable, not a neutral oracle), and it **structurally cannot
  cover the prognostic half** (the sub-second leaf pole returns, and `column_state_t` carries no leaf
  state to integrate). RK4's correct, already-implemented role is a **mild-window cross-check** (shares no
  code with the BE path → agreement rules out a shared-bug false pass). The plan promoted a cross-check to
  a baseline.

### 2.4 The adopted scheme set

1. **Two production integrators: SPLIT+Picard and ARK2**, cross-validated at **diagnostic leaf** now
   (the achievable, rigorous milestone).
2. **Ground truth = a self-converged STIFF reference**: drive the existing implicit machinery
   (`adaptive_ark_march` / `adaptive_imex_march`) with **tolerance-halving** until the solution stops
   moving (Richardson self-convergence). Optionally add a code-disjoint stiff oracle (Radau IIA / SDIRK)
   for full independence — extra build, not required for the minimal correct reference.
3. **RK4 kept exactly as-is** — a mild-dt cross-check on the reduced system. Not scaled up.
4. **ARK4 deferred** behind lifting the operator split (soil water into the coupled solve or a
   higher-order/Strang composition; psi via a higher-order operator-split composition since the exact-exp
   form cannot enter a Butcher tableau).
5. **Leaf arrowhead deferred**, built only to *quantify* prognostic leaf given its negligible physical
   effect — a validation task, not a blocker.

This refactor **keeps both production schemes bit-behaviour-invariant**; it does not build new numerics.
The deferred items (2 stiff-reference, 4 ARK4, 5 arrowhead) are roadmap, not part of the reorg.

---

## 3. Dead code & duplication inventory

### 3.1 Delete now (low-risk, no behaviour change)
- `itoa` — `meds_main.f90:256-262`, zero call sites anywhere.
- Duplicated comment block — `meds_column_dynamics.f90:325-329` (verbatim paste, twice).
- Unused import `SCHEME_PICARD_COUPLED` — `meds_fast_loop.f90:18` (branch lives in `column_fast_step`).

### 3.2 Superseded — retire as a *separate decision* (not folded into the reorg)
- IMEX-Euler tier `imex_euler_column_step` + `adaptive_imex_march` (`meds_ark_stepper.f90:242, :55`) —
  test-only, superseded by ARK2 / `adaptive_ark_march`.

### 3.3 Test-only but LOAD-BEARING — keep, annotate, do NOT delete
- The RK4 oracle chain: `rk4_column_step` + `column_derivs` + `state_axpy`/`state_accum` + `column_tend_t`.
  It is the *independent* cross-validation of the production ARK (no shared code with the BE path).
  Relocate to a clearly test-support module (`meds_fast_rk4_oracle`) so it is not mistaken for production
  or optimized away.

### 3.4 Inert scaffolding — delete or clearly mark reserved
- `soil_water_coupling` selector + `SOILH2O_LAGGED/COUPLED` (`meds_column_dynamics.f90:116, 126-127`) —
  **no dispatch exists**; both "modes" re-solve from state^n. Advertises behaviour it doesn't have.
- `relax` argument threaded through the entire production ARK stack (`adaptive_ark_march` → `ark2` →
  `column_be_stage`) — documented "ignored"; `newton_surface_solve` doesn't accept it. A live-looking
  knob with zero effect. Remove from the API.

### 3.5 Duplication — extract (higher effort, golden-anchor-gated → Phase 5)
- **Headline dup:** `build_column_frozen` (`:968-1021`) re-implements the split pre-pass (`:366-401`)
  **verbatim** — `rho_mol`/`e_air`/`leaf_gas_exchange`/`h_coeff_f`/`g_tr_f`/respiration/NEE. Self-documented
  as load-bearing (exists so split/ARK GPP stay bit-identical). Extract ONE `pure column_prepass(...)`.
- Smaller inline dups fold into that: root-soil aggregation (`:349` vs `:961`), CAS caps (`:411` vs
  `:1028`), `soil_psi_root` (`:600` vs `:1057`).
- **Physics that belongs in kernels:** leaf/wood linearization coefficients (`h_coeff_f`/`lw_slope`/
  `le_slope`) → `meds_vegetation_biophysics`; inline inter-store enthalpy transport (snow-melt poke
  `:463-464`, boundary advection `:670-674`, the `rain_temp=tsupercool_liq` double-count guard) →
  `soil_energy_step_implicit`.
- Low-priority: pre-allocate `surface_derivs` tend buffers (heap-allocs every call in the hottest loop,
  `meds_column_derivs.f90:221`); dedup `growth_hist_pos` ring-buffer advance
  (`meds_vegetation_dynamics.f90:98` vs `meds_demography_capi.f90:131`) and `flatten_pheno_params`.

**Phase 5 execution note (post-hoc):** the headline dup, the three smaller inline dups, and the
leaf/wood linearization coefficients were all extracted as planned — `column_prepass` (`meds_fast_ark`,
called from both `column_fast_step` and `build_column_frozen`), `root_weighted_psi`
(`meds_fast_time_derivs`), and `sensible_heat_coeff`/`leaf_transp_coeff`/`lw_emission_slope`/
`le_conductance_flux` (`meds_vegetation_biophysics`) — all bit-identical (bare code motion of a
self-contained formula). The **inter-store enthalpy transport → `soil_energy_step_implicit`** move was
investigated and **deliberately deferred, not done**: the split applies it as a pre-solve STATE POKE
(`col%soil_energy(1) +=` before `soil_energy_step_implicit` even inverts `t_n`), while
`column_be_stage`'s ARK path already folds the identical physical quantity into `eforc%root_heat_sink`
as an implicit-BE SOURCE TERM (`q_src`, evaluated inside the solve). These are two different numerical
treatments of the same boundary term, not one duplicated formula — unifying them would be an algorithm
change (a real, if probably small, perturbation to the split's committed values), which contradicts this
phase's bit-identical mandate. Moving *only* the split's poke-before-solve verbatim into the kernel
(as an optional forcing field, unused by the ARK call site) was considered and rejected as added
`energy_forcing_t`/kernel-API surface for zero behavioral or duplication benefit (it's already a single
occurrence, not two). Left as a flagged follow-up: unifying split/ARK's water-enthalpy-advection
treatment is a genuine numerics design question, not a mechanical refactor.

---

## 4. Target structure

### 4.1 Naming convention

- **`meds_*_dynamics.f90` = a tier orchestrator** that `advance_one_step` fans out to. The `_dynamics`
  suffix (established by `vegetation_dynamics`, an orchestrator/policy driver — not a numerical stepper)
  now applies uniformly across the fast and slow tiers. The fast/slow naming asymmetry (fast by
  *timescale*, slow by *domain*) is deliberate: the fast processes are one tightly-coupled stiff column;
  the slow processes are separable domains.
- **`meds_fast_split` / `meds_fast_ark` = scheme backends** the fast orchestrator dispatches to (peers).
- **`meds_fast_time_derivs` / `meds_fast_types` = shared machinery** both backends use.
- Result: the word **"stepper" survives only in `meds_stepper`** (the calendar-tick conductor); no two
  files at different levels share a name; the two schemes are named as peers.

### 4.2 Target module list

| Module | was | Contents |
|---|---|---|
| `meds_main.f90` | (same) | program entry, calendar loop, I/O orchestration — anchor-safe to touch |
| `meds_stepper.f90` | (same) | `advance_one_step` cadence + its 3 load-bearing guards (fast_ctx-present, step_start-present, FAST-before-SLOW) |
| **`meds_fast_dynamics.f90`** | `meds_fast_loop.f90` | fast-tier orchestrator: the sub-daily loop, reservoir gather/writeback, forcing overlay, fast→slow handoff. Entry sub **`fast_dynamics`** (renamed from `run_fast_biophysics`, peer of `vegetation_dynamics`) |
| **`meds_fast_types.f90`** (new leaf) | — | the 4 `column_*_t` working buffers + the ARK POD types (`column_state_t`, `column_frozen_t`, `surface_*_t`, `column_tend_t`, `stage_bflux_t`, `column_bflux_t`). Composes — does not copy — the biophysics/state types (§4.3) |
| **`meds_fast_time_derivs.f90`** | `meds_column_derivs.f90` | `surface_derivs` / `column_derivs` (the RHS evaluators), minus the POD types (moved to `meds_fast_types`) |
| **`meds_fast_split.f90`** | (core of `meds_column_dynamics.f90`) | `column_fast_step` (operator-split + Picard) + the shared pre-pass |
| **`meds_fast_ark.f90`** | (ARK bits of `meds_column_dynamics` + `meds_ark_stepper`) | `column_fast_step_ark` + `build_column_frozen` + `ark2_column_step` / `adaptive_ark_march` / `column_be_stage` / `newton_surface_solve` / `jac_surface` / `advance_hydraulics_full` / `bflux_*` / `state_*` |
| **`meds_fast_rk4_oracle.f90`** (test-support) | (RK4 + IMEX-Euler bits of `meds_ark_stepper`) | `rk4_column_step` + the superseded IMEX-Euler tier — relocated so it is unmistakably test-only |
| `meds_vegetation_dynamics.f90` | (same) | slow veg-tier orchestrator (`carbon_growth` split internally) |
| ~~`meds_biogeochem_dynamics.f90`~~ | — | **deferred** (§8) |

The final `advance_one_step` story: **conductor** (`meds_stepper`) → peer **tier orchestrators**
(`meds_fast_dynamics` / `meds_vegetation_dynamics` / future `meds_biogeochem_dynamics`) → the fast
one **dispatches** to scheme backends (`meds_fast_split` / `meds_fast_ark`) → which call the shared **RHS**
(`meds_fast_time_derivs`) over the shared **types** (`meds_fast_types`).

### 4.3 `meds_fast_types` scope + overlap with the two existing type modules

`meds_fast_types` is a **driver-scope** types module. It does NOT belong in `shared/state` (persistent,
core-owned state) or `biophysics` (kernel I/O contracts). Verified relationships:

- **vs `biophysics/meds_biophysics_types.f90` → composition, not duplication.** `column_config_t` *nests*
  the biophysics param/opts types wholesale (`aero_cfg_t`, `veg_thermal_params_t`, `soil_params_t`,
  `soil_thermal_params_t`, `energy_opts_t`, `soil_opts_t`, `snow_params_t`,
  `meds_column_dynamics.f90:88-118`) plus plant types + driver scalars. It's the **aggregation seam** that
  `build_fast_context` fills (the `[soil]/[energy]/[snow]/[aerodynamics]` TOML wiring lands here).
  `column_forcing_t` / `surface_frozen_t` have a *conceptual* shape-overlap with the kernel forcing/env
  types but are driver-assembled bundles — "similar shape, different scope."
- **vs `shared/state/meds_column_state_types.f90` → one real *representational* overlap.** The persistent
  per-store structs (`cas_state_t`, `soil_column_t`, `soil_energy_column_t`, `snow_column_t`) are the
  source of truth on `site_t`. The working buffers mostly *reference* that state (`column_cohort_t` gathers
  the cohort SoA; the fast step adopts the reservoirs into a `patch_biophys_t`). BUT the ARK
  `column_state_t = CAS(3) + soil_energy(:) + theta(:) + psi(:,:)` is a **flat re-packing** of the same
  prognostic quantities the per-store structs hold — deliberate (the ARK needs a contiguous vector for
  `state_axpy`/`state_wrms`/tableau linear combinations; the per-store structs serve the persistent hub +
  the split path's per-store kernels). **Not a bug; a representation choice.** Document it in the header;
  do NOT try to unify (that's a deeper change, out of scope).

---

## 5. Invariants & what is allowed to change

**Correction to earlier framing (verified):** the "golden anchor" is a **1e-3 tolerance** check
(`test/test_picard_coupling.f90:82`), driven **directly through `column_fast_step`** (`:175`), NOT a
bit-identity check and NOT routed through `meds_main`/`meds_stepper`/`meds_fast_loop`. Therefore
reorganizing the entry/orchestration files is **provably anchor-safe**; scrutiny belongs on the split-path
body + the four `column_*_t` types. Discipline: hold the anchor to *printed digits*; treat 1e-3 as the
hard gate; the tighter guard is the per-kernel/whole-column budget closures in `test_column_dynamics` /
`test_column_ark`.

**MUST hold invariant:**
- `column_fast_step` signature + split-path body + its ~10 golden-ordering items (ARK-dispatch-before-
  mutation, `te=tcas`, `a_leaf=0`, snow-before-snapshot, `src_frac` after hydrology/before CAS,
  soil-thermal-inside-loop, …).
- The four `column_*_t` type layouts/semantics (consumed by field by tests + fast orchestration).
- `surface_derivs` — the ONE evaluator shared by both integrators; any edit moves both.
- The fast→slow SoA accumulator produce/reset/read contract + FAST-before-SLOW order.
- The nvfortran `h_flux` returned-arg contract (`meds_fast_loop.f90:341-343`) — must be surfaced as a
  returned arg, not recomputed from post-call state (nvfortran miscompiles the post-call read to 0).
- All `bind(c, name=...)` symbols in `src/capi/meds_demography_capi.f90` (Python-facing ABI).

**Allowed to change freely:** `meds_main` boot/loop structure, comment fixes, dead `itoa`; the location of
`meds_stepper`'s cadence logic *iff* the FAST-before-SLOW order + the two present-guards are preserved
verbatim; file/module names (mechanical; a partial rename fails loud via the shared `.mod` dir); stale
header comments (`meds_fast_loop.f90:8-9` claim leaf_temp/psi reseeded — code persists them).

**Structural landmines (verified):**
- `CMakeLists.txt:181` globs `src/driver/*.f90` with `file(GLOB)` — **NOT `GLOB_RECURSE`**. Any move into
  a subfolder silently drops files from `meds_aux`. Flip to `GLOB_RECURSE` (or list subfolders) before any
  subfolder move.
- The four `column_*_t` types cross the split/ARK boundary; splitting `meds_column_dynamics` without first
  extracting them to a leaf module **cycles** `dynamics ↔ ark_stepper ↔ derivs`.
- Module-rename ripple reaches `src/capi` (`use meds_vegetation_dynamics`) and every test — lockstep the
  `use` updates; keep the C ABI names frozen.

---

## 6. Phased plan (each phase behind green ifx Debug + nvfortran multicore + anchor)

**Phase 0 — Baseline capture.** Clean-build both compilers; run full `ctest`; record: (a) the printed
`tc_split(54)` to full digits, (b) the whole-energy/whole-water worst residuals from
`test_column_dynamics`/`test_column_ark`, (c) the *set* of nvfortran failures (the pre-existing
`plant_hydraulics` flake — so a new failure isn't hidden), (d) a **fast-on + carbon-mode** short run
(construct one — not in the default suite) + its summary/AGB/LAI, (e) a `meds_c`/`meds_plant_c` link +
Python import smoke.

**Phase 1 — Dead code + comments only.** `itoa`; the duplicated comment (`:327-329`); the unused
`SCHEME_PICARD_COUPLED` import; stale header comments. Zero semantics. Verify: identical `ctest`, identical
anchor digits, both compilers link.

**Phase 2 — Flat-dir module/file renames** (no subfolders yet): `meds_fast_loop → meds_fast_dynamics`
(+ entry sub `run_fast_biophysics → fast_dynamics`), `meds_column_derivs → meds_fast_time_derivs`.
Lockstep-update every `use` across `src/`, `test/`, `src/capi/`. Do NOT rename `bind(c)` names. Verify:
clean-build-dir both compilers (partial rename fails loud), `ctest` set identical, anchor identical,
`meds_c` links + Python smoke passes.

**Phase 3 — Extract `meds_fast_types` (leaf module).** Move the four `column_*_t` + the ARK POD types out
of `meds_column_dynamics`/`meds_column_derivs`. This is prerequisite plumbing for any file split. If moving
into subfolders, **first flip `CMakeLists.txt:181` to `GLOB_RECURSE`** and diff the object list. Verify:
anchor + all parity tests (`picard_coupling`, `column_dynamics`, `column_ark`, `column_derivs`,
`fast_loop`) identical.

**Phase 4 — Structural non-physics moves.** Split `meds_column_dynamics` → `meds_fast_split` (the split
stepper) + `meds_fast_ark` (the ARK dispatch + `build_column_frozen`); fold `meds_ark_stepper`'s ARK guts
into `meds_fast_ark` and relocate the RK4/IMEX-Euler oracle into `meds_fast_rk4_oracle`; split
`carbon_growth` internally; retire the inert `soil_water_coupling` selector + the vestigial `relax` arg.
Verify: anchor + budgets + the **fast-on carbon smoke** (guards the SoA handoff) + Python smoke.

**Phase 5 (last, highest scrutiny, alone) — the pre-pass dedup + any `column_fast_step` touch.** Extract
the one `pure column_prepass(...)` from the verbatim split/frozen copy; fold the smaller inline dups; move
leaf/wood linearization coeffs to `meds_vegetation_biophysics` and the inter-store enthalpy transport into
`soil_energy_step_implicit`. This is a *behaviour* change (can perturb the anchor + split/ARK GPP parity).
Gate on: anchor to printed digits, `test_column_ark` gpp/budget parity, both compilers, fast-on carbon run
bit-stable. Do not batch with anything else.

**Per-phase verification (every phase):** run BOTH ifx Debug and nvfortran multicore (the `h_flux`
post-call-recompute trap and array-temp regressions are ifx-invisible); anchor to printed digits; a fast-on
carbon smoke (the default fast-off suite can't catch a tick reorder / dropped optional-arg guard); a
`meds_c` Python import smoke; a clean build dir (stray root `*.mod` = out-of-tree pollution masking a
missed `use`).

---

## 7. IMEX-Euler retirement — separate decision

The IMEX-Euler tier (`imex_euler_column_step` + `adaptive_imex_march`) is test-only and superseded by
ARK2. Retiring it is a **separate decision from this reorg** (do not fold into a phase) — but if retired,
it moves out of `meds_fast_rk4_oracle` cleanly. The RK4 oracle chain stays (independent cross-check, §3.3).

---

## 8. Deferred / out of scope

- **`meds_biogeochem_dynamics` (slow soil carbon).** A NEW feature, not a reorg: needs a per-patch
  `soil_carbon_t` on `site_t`, the demography→litter→Rh seam (`carbon_growth` computes litter at
  `meds_vegetation_dynamics.f90:320-326` but never routes it), and a slow-tier call in `advance_one_step`.
  Its home will follow the `*_dynamics` convention when built.
- **ARK4** (needs the operator split lifted — §2.3/2.4).
- **Prognostic-leaf arrowhead** (needs the 3×3 leaf↔CAS solve; validation-only given negligible physical
  effect).
- **A self-converged stiff reference + Radau/SDIRK oracle** (the ground-truth work — separate from the
  reorg).

---

## 9. Open items for implementation kickoff

- Confirm scope of the **first PR**: recommended = Phases 1–3 (dead code + renames + types-extraction),
  since Phases 4–5 concentrate the golden-anchor risk. Phases 4 and 5 land as their own PRs.
- Decide the **entry-sub name** finalization: `fast_dynamics` (adopted) — confirm no downstream doc/CLAUDE.md
  references to `run_fast_biophysics` are missed in the lockstep rename.
- Decide whether to **retire the IMEX-Euler tier** in the same window as the reorg or separately (§7).
