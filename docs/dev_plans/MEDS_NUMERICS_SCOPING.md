# MEDS Fast-Loop Numerics — Scoping & Framework Design

**Status:** scoping (2026-07-22) + **P0 implementation started** (2026-07-22, this pass — see §8a).
Supersedes the "remaining work" framing of `MEDS_IMEX_ARK_DESIGN.md` (the ARK integrator it planned
is now *implemented and merged* — see §1.3). Companion to `MEDS_IMEX_ARK_DESIGN.md`,
`MEDS_COLUMN_DYNAMICS_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`, `MEDS_INTEGRATOR_TEST.md`.

---

## 0. Purpose — the three goals

This round of numerical development is **not** about crowning one "best" integrator. The stiffness
structure of the fast loop is dt-dependent and process-dependent, so the right scheme changes with the
configuration. The deliverable is the **infrastructure to decide**, made of three layers:

- **(a) Faithful error control.** A uniform, per-scheme machinery for *strict* error and *conservation*
  control — embedded local-error estimates, consistent WRMS tolerance norms, and machine-precision
  budget closure that can be **enforced** (a hard failure), not merely tracked. "Faithful" = the
  reported answer is provably within a stated accuracy/conservation tolerance of the reference.

- **(b) A scheme-vs-process-vs-`dt` benchmark harness.** The capacity to run the **same** column under
  different integrators, at different `dt_fast`, and with **different process complexity** (energy on/off,
  hydraulics on/off, soil water on/off, CO₂ on/off, snow on/off), and to emit structured
  accuracy / work / conservation / wall-clock data for each combination.

- **(c) A data-driven integrator-selection rule.** A `select_integrator(cfg)` decision keyed on
  `dt_fast`, the active process mask, the accuracy/strictness target, and the hardware
  (CPU-serial / CPU-multicore / GPU), *derived from the benchmark data of (b)* rather than by intuition.

Everything below serves these three. The scheme families (§3), the error-control API (§4), the harness
(§5), and the selection rule (§6) are the four concrete work products; §7–§8 fold in the supporting
robustness / GPU / accuracy pathways and phase them.

---

## 1. Where we are today

### 1.1 The default fast step

Default is **`INTEG_SPLIT` + `SCHEME_SPLIT_SEQUENTIAL`** (`meds_config.f90`; TOML `[fast]`), i.e. a
**single operator-split Gauss–Seidel sweep** per `dt_fast` (default 900 s) in
`meds_fast_split.f90:column_fast_step`. Per sub-step:

0. **Frozen pre-pass** (`column_prepass`, in `meds_fast_ark.f90`) — Monin–Obukhov aerodynamics +
   boundary layers, leaf gas exchange (GPP / gs / Rd), stem+root+heterotrophic respiration, CAS
   capacities/conductances. Held constant for the whole sub-step (ED2 DTLSM freeze — faithful, see §2.1).
1. **Snow** advanced operator-split outside the iterate at the lagged CAS.
2. **Snapshot `state^n`**, then the sequential body: leaf-T diagnostic → soil water + `src_frac` supply
   limiter → plant hydraulics → ground fluxes → **CAS 3-twin implicit-in-atm box** → **soil thermal
   BE-Thomas**.

Each store is advanced by an **L-stable backward-Euler** sub-solve; the cross-store coupling is a single
Gauss–Seidel pass, so the step is **first-order in `dt_fast`**. Whole-column energy + water budgets close
to machine precision (`meds_fast_split.f90` §7b). `dt_fast` is fixed — the only adaptivity is buried
inside the Richards and plant-hydraulics sub-solvers (step-doubling).

### 1.2 The Picard coupled-surface path

`SCHEME_PICARD_COUPLED` (opt-in) wraps an outer fixed point over
`{leaf → soil water/src_frac → CAS twins → soil thermal}`, `picard_max_iter=20`, `tol_temp=1e-3 K`,
`tol_shv=1e-6 kg/kg`, `relax=0.5`. It buys **coupling-consistency** (removes the Gauss–Seidel lag → the
fully-coupled backward-Euler solution) but stays **first-order** in `dt`. It is **required** for a
prognostic leaf (the explicit split with a leaf storage term 2·`dt`-oscillates, ~1.7 K midday spikes;
Picard damps it to ~0.2 K).

### 1.3 IMEX-ARK — implemented and merged (correcting the record)

The IMEX-ARK integrator planned in `MEDS_IMEX_ARK_DESIGN.md` is **built, tested, and on `main`** (opt-in
`[fast].time_integrator="ark"`), reachable by dispatch from `column_fast_step`. What exists in
`meds_fast_ark.f90`:

- **`ark2_column_step`** — a real **ARS(2,2,2)** L-stable, stiffly-accurate ESDIRK2 tableau
  (`γ=1−1/√2`, `β=1+√2`) with a **free embedded first-order error estimate**.
- **`newton_surface_solve` + `jac_surface`** — a **2×2 Newton "arrowhead"** on (CAS enthalpy, CAS shv)
  with leaf-T and ground static-condensed, **finite-difference** Jacobian, damped-diagonal + Armijo +
  eval-cap guards.
- **`adaptive_ark_march`** — an **embedded-error WRMS adaptive controller** (2 solves/step) plus a
  GPU-warp-uniform fixed-substep path.
- **`column_be_stage` / `advance_hydraulics_full`** — the shared BE stage (one hydraulics path) with a
  b-weighted boundary-flux ledger that closes **6–7 conservation budgets to machine precision** on the
  ARK path.
- Cross-validated against the **RK4 oracle** (`meds_fast_rk4_oracle.f90`, test-only) in
  `test_column_derivs`; exercised end-to-end in `test_column_ark`; validated on a **30-yr Ithaca run**
  that matches the SPLIT (LAI 1.33 vs 1.36, AGB 0.865 vs 0.87).

**What remains MVP / not-as-designed:**

1. **Hydraulics and soil-water are operator-split *out* of the tableau** — ψ frozen through the stages
   (advanced once by exact matrix-exp), θ committed from a scratch Richards solve. This is *why the
   coupled-state order collapses to ≈ p1*: order-2 is realized only on cleanly-integrated twins
   (CO₂ tests p ≥ 1.9). Confirmed by a reverted 3×3-arrowhead experiment: the ~1.2 order cap is the
   **whole nonlinear surface Newton at the `β`-extrapolated base**, not the CAS↔soil-top split.
2. **Prognostic leaf (P4) and prognostic wood (P2) are hard `error_stop`** under ARK.
3. The arrowhead is a **2×2 FD-Jacobian**, not the designed dense-Schur border; the shipped tableau is
   **ARS(2,2,2)**, not the design's named Kennedy–Carpenter ARK2(1)3L[2]SA.
4. **Free-drain bottom BC only** (ponding / aquifer / Zeng–Decker are non-smooth → non-ESDIRK-integrable).

**Production `dt_fast` evaluation already ran** (`MEDS_INTEGRATOR_TEST.md`, 24-run matrix): **`dt_fast=900s`**
for tight targets; SPLIT is actually *most accurate* on energy/water at coarse `dt` and ~30 % cheaper than
ARK; **ARK-fixed(1 substep) is a bad config** (≈ IMEX-Euler on a stiff system — LE/H blow up at coarse
`dt`) → use ARK-**adaptive**. This is exactly the kind of data goal (c) will systematize.

### 1.4 Slow loop

Daily 7-pool CENTURY matrix ODE (`meds_soil_biogeochem.f90`): production step is **forward Euler** on the
fast loop's accumulated `xi_int·K`; **exact augmented matrix-exponential** (`step_expm`) + SASU exist but
are **used only for spin-up**, not the daily step. Out of scope for the fast-loop families but relevant to
goal (a) (positivity) — see QW3 in §7.

### 1.5 Current-state table

| Subsystem | Integrates | Method | Order / stability | Default? |
|---|---|---|---|---|
| Split column `meds_fast_split.f90` | CAS 3 twins, soil E/θ, leaf/wood T, ψ, snow | operator-split single G–S sweep of per-store BE | **1st-order** coupling; sub-solves L-stable | **YES** |
| Picard `SCHEME_PICARD_COUPLED` | same | outer fixed point (relax 0.5) | 1st-order, coupling-consistent | opt-in |
| IMEX-ARK `meds_fast_ark.f90` | CAS+soil-heat in tableau; ψ/θ split out | ARS(2,2,2) ESDIRK2 + 2×2 Newton + embedded adaptive | L-stable; **coupled step ≈ p1** | opt-in |
| RK4 oracle `meds_fast_rk4_oracle.f90` | whole column via `column_derivs` | explicit classical RK4 | 4th-order, **not L-stable** (stiff → tiny `dt`) | test-only |
| Soil heat `meds_soil_energy.f90` | per-layer internal energy | BE-Thomas, enthalpy inverter | 1st-order, L-stable, **`nsub=1`** | in split |
| Soil water `meds_soil_water.f90` | per-layer θ (Richards) | BE-Thomas + step-doubling | 1st-order, L-stable | in split |
| Plant hydraulics `meds_plant_hydraulics.f90` | per-cohort ψ(leaf,wood) | **exact 2×2 matrix-exp** + step-doubling | unconditionally stable | always-on |

---

## 2. The RHS partition (the foundation both scheme families share)

The single most important structural asset is already in the tree: **`column_derivs`**
(`meds_fast_time_derivs.f90`) — a **pure, side-effect-free whole-column RHS** `f(y)`. The RK4 oracle
marches it explicitly, the ARK marches it IMEX, and the operator-split is the degenerate single-stage
case. **All four scheme families are tableaux over this one RHS.** Do not grow four parallel steppers;
build one harness parameterized by a scheme descriptor.

### 2.1 Three categories — frozen / implicit / explicit(+exponential)

Verified against the ED2 reference (`ED2/ED/src/dynamics/rk4_driver.F90:259`): ED2 computes
`canopy_photosynthesis` **once per DTLSM, before** `integrate_patch_rk4`, and the RK stages only *read*
the frozen `gsw_open/closed`, `A_open/closed` (→GPP), `leaf_resp` (Rd), `fs_open`. Photosynthesis is
**never re-called inside the substep loop**. Only the transpiration *flux* (frozen conductance × evolving
leaf→CAS vapour gradient) and the aerodynamic conductances (`update_diagnostic_vars → canopy_turbulence8`)
are re-evaluated per stage. **So freezing gs/GPP/Rd per `dt_fast` is ED2-faithful and deliberate** (the Ci
root-find is expensive; leaf biochemistry is quasi-stationary over the substep). MEDS's `column_prepass`
matches this exactly and **must keep it** — the partition operates only on what is *integrated*.

| Category | Components | Treatment | Rationale |
|---|---|---|---|
| **0 — Frozen** | gs, GPP, Rd, radiation absorption (SW/PAR/LW) | held constant once per `dt_fast` (in neither tableau) | ED2-faithful; expensive nonlinearity, quasi-stationary |
| **Implicit (stiff)** | leaf energy (diagnostic=algebraic / prognostic=stiff ODE, §2.2), CAS enthalpy+vapour, ground/soil-top, soil-heat column, soil water (Richards) | DIRK / ESDIRK stage (Newton arrowhead) | short timescales, tight bidirectional coupling, diffusion stiffness |
| **Exact-exponential** | plant hydraulics ψ (leaf–wood RC) | `exp(M·h)` stage-map (Lawson/ETD), **not** a DIRK component | linear per frozen coeffs → exact & cheap; fastest mode ~17 s |
| **Explicit (non-stiff)** | prognostic wood energy (§2.2), CAS-CO₂, the transpiration/aero flux couplings | ERK stage (the additive explicit tableau) | no self-stiffness; this is the "adaptive-explicit" half |

Two facts make this partition honest:

- It is **state-based and `dt`-dependent.** At `dt_fast=900 s` almost everything is stiff → the split is
  implicit-heavy (correct). At small `dt` (e.g. 60 s) wood / CO₂ / deep soil migrate explicit. The
  partition should therefore be a **per-component tag**, ideally auto-placed by a stiffness diagnostic
  (`τ_local > k·dt` ⇒ explicit). This tagging is itself a benchmark axis (§5).

- **Hydraulics *and* soil-water must be *in* the additive scheme** for order > 1 on the coupled state.
  Today they are split out (§1.3.1), which is exactly why ARK's order collapses. The fix — ψ as an
  exponential stage-map inside the composition, θ as an implicit additive component — is the single most
  important thing for the IMEX family to actually beat SPLIT (§3.4).

### 2.2 Vegetation energy — the diagnostic↔prognostic axis (new RHS `vegetation_energy_time_derivs`)

Today `column_derivs` integrates CAS + soil + hydraulics; **leaf/wood temperature is diagnostic**
(quasi-steady) in `surface_derivs`. Prognostic leaf/wood stores exist only on the SPLIT/Picard path (BE
`veg_energy_step_implicit`, `meds_vegetation_biophysics.f90`) and are **`error_stop`-gated under ARK** — so
we currently *cannot* compare diagnostic vs prognostic vegetation energy across schemes. **ED2 integrates
leaf and wood energy prognostically** (verified, `ED2/ED/src/dynamics/rk4_derivs.f90`: `dinitp%leaf_energy`,
`dinitp%wood_energy`, `dinitp%veg_energy = leaf+wood` are real RK4 tendencies from absorbed SW + net LW +
sensible + latent, coupled by an internal leaf↔wood heat flux `qwflux_wl` at :2122–2123; temperature is
diagnosed from the prognostic internal energy; slaved to 0 only when a cohort is too small to resolve).

**So this development adds `vegetation_energy_time_derivs`** — the side-effect-free RHS for
`d(leaf_energy)/dt` and `d(wood_energy)/dt` — making leaf/wood internal energy **optional prognostic
members of `column_state_t`**, integrated by *every* scheme (this removes the ARK `error_stop` gate — the
deferred "prognostic leaf under ARK arrowhead", P4 of the leaf-wood-energy work — `[RESUME]`). It is the
vegetation-energy analogue of "bring hydraulics/water into the tableau": expose the tendency so the store
becomes a first-class integrated state, and the benchmark can compare the modes:

- **Diagnostic leaf** = the algebraic (infinitely-stiff / steady-state) limit, static-condensed into the
  surface block (today's default, cheapest). Leaf thermal inertia is physically ~0.015 K, so we expect this
  to be an excellent approximation — the point is to **measure** it, not assume it.
- **Prognostic leaf** = a genuinely **stiff** ODE (heat capacity ~seconds) → **implicit**, folded into the
  surface Newton arrowhead as a real state (extends the 2×2). This is exactly why the explicit split
  2·`dt`-oscillates and needs Picard; under the arrowhead it is stable.
- **Prognostic wood** = **non-stiff** (τ tens of min–hours) → **explicit** ERK component. Wood inertia is
  *not* negligible, so the diagnostic-vs-prognostic comparison matters more for wood than for leaf.

**Per-cohort routing by heat capacity (the veg-energy instance of stiffness auto-placement, §2.3).** A
cohort's leaf/wood thermal timescale is `τ = C_total / g_total`, and `C_total` scales with `bleaf·nplant`
(or `bsap·nplant`). For a **tiny-hcap** cohort `τ → 0`, so the prognostic ODE relaxes to steady state
*within* the step — i.e. the diagnostic solution **is** its stiff limit, and it is the source of the
worst stiffness in the whole system (the ~6.5e5 ratio). So route per cohort: **large `C_total` →
prognostic** (real thermal memory, integrated); **`C_total` below a `τ < ε·dt_fast` threshold →
diagnostic** (algebraic steady state). This is strictly better than ED2's binary
resolvable/slaved-to-CAS flag: MEDS's diagnostic branch still solves a **real, exactly energy-conserving**
steady-state balance (§ below), not a slave. Benefit: removes the stiffest modes from the integrated
set (better-conditioned arrowhead Newton, bigger explicit steps, fewer states) *and* avoids the tiny-hcap
blow-ups ED2's threshold exists to dodge. The threshold is safe because at the crossover `C_total·ΔT` is
tiny, so the prognostic↔diagnostic switch is near-continuous in energy; use a small hysteresis/blend band
so the switch is smooth for the error controller. The threshold scales with `dt_fast` — a coarser step
makes more cohorts diagnostic-eligible — which is exactly why it is an *auto-placement* criterion (§2.3).

**Conserving the diagnostic↔prognostic transition.** Switching a cohort's category changes the *store
set* — a promoted cohort suddenly carries a tracked leaf/wood energy store `E = C·T`. This is **not** energy
creation: a diagnostic leaf physically *has* internal energy `C·T`; the diagnostic model just doesn't track
it (it assumes `C→0`, so the store telescopes out, below). Promotion starts tracking it. It is the **same
class of event** as snow appearing/disappearing (a store toggles) and cohort birth/death, and MEDS already
has the machinery — the snow model's **paired mass+enthalpy transfers** when its store toggles
(`meds_fast_split.f90`). Reuse that pattern:
- **Promotion (diagnostic→prognostic):** seed `E ← C·T_diag` from the cohort's current diagnostic
  (steady-state) temperature ⇒ leaf temperature is **continuous** (no spurious flux) and the seed is the
  **physically-correct** energy. Log `+C·T_diag` as an explicit transition transfer *into* the accounted
  system (the thermal energy the `C→0` diagnostic model was neglecting) so the ledger closes with a known,
  bounded term, not a silent jump.
- **Demotion (prognostic→diagnostic):** log `−C·T_prog` out at continuous temperature and stop tracking.
  Over a promote→demote cycle the transfers ≈ cancel (`C·(T_diag−T_prog) ≈ 0`) ⇒ no drift.

The transfer is **bounded and small by construction**: routing at `τ = C/g < ε·dt_fast` means `C < ε·dt·g`,
so `C·T` at the switch is small; a **hysteresis band** (promote at `τ<ε₁·dt`, demote at `τ>ε₂·dt`,
`ε₂>ε₁`) makes switches rare. L2-strict conservation (§4) treats the transition transfer as an *accounted*
term. **Fully-rigorous alternative (deferred):** tie the store's thermal-mass appearance to the
biomass-growth seam — source `ΔC·T` when leaf tissue is created on the slow loop — so the store builds up
through flux-balanced growth and there is *no* jump at the threshold at all; needs the carbon→biomass seam
to transport tissue enthalpy.

*Why the diagnostic branch conserves energy exactly:* `veg_energy_diagnostic`
(`meds_vegetation_biophysics.f90`) solves the linearized **steady-state** balance (`a_store=0`), so by
construction the net radiation the leaf absorbs equals the sensible + latent it hands the CAS:
`drnet ≡ dh + L·transp` (algebraic identity in the returned fluxes). The leaf **stores nothing**
(`leaf_store0 = leaf_store1 = 0`), so it telescopes out of the whole-column ledger — it is an
instantaneous energy **router**, not a reservoir, and `budg%whole_energy` closes to machine precision.
Conservation is exact *independent of accuracy*: the diagnosed leaf temperature is a linearization about
`t_emit`/`t_cas`, so it can sit slightly off the true nonlinear steady state, but the ledger always
closes because it uses those same linearized fluxes consistently.

**Leaf↔wood sapflow enthalpy coupling — add it in this PR (resolves open decision #5).** ED2 couples the
prognostic leaf and wood energy stores through the **sapflow** heat flux `qwflux_wl` = `wflux_wl ·
u_liq(T_upwind)` (Appendix B; upwind on the wood→leaf water flux). MEDS **already computes the water
flux** (`flux%sapflow`, `meds_plant_hydraulics.f90:208`, on `hydro_flux_t`) but **does not advect its
enthalpy** into any tissue store — because leaf/wood energy is diagnostic today (no store to feed). With
prognostic stores, add the term to the veg-energy RHS: `+ sapflow·u_liq(T_up)` to `d(leaf_energy)/dt`,
`− sapflow·u_liq(T_up)` to `d(wood_energy)/dt`, `T_up = wood_temp` if `sapflow ≥ 0` else `leaf_temp`. It
is cheap (the flux exists) and telescopes in the ledger (internal transfer). This is the physically-correct
leaf↔wood thermal link and the reason a prognostic-leaf/prognostic-wood pair should be co-integrated.

The leaf/wood energy stores join the fast-state SoA that **BB1 hoists** (§7), so
`vegetation_energy_time_derivs` and BB1 are **co-designed against one state layout** (the fast prognostic
set = CAS twins, soil E/θ, per-cohort ψ, per-cohort leaf/wood energy, snow).

### 2.3 Stiffness auto-placement (the general form of §2.2's hcap routing)

The partition (§2.1) is `dt`-dependent, so *fix it by measurement, not by hand*. Tag each component
implicit / explicit / diagnostic by its local timescale vs the step: `τ_i < ε·dt_fast` ⇒ the reservoir
equilibrates within the step ⇒ treat it **diagnostic** (algebraic) or **implicit**; `τ_i ≫ dt_fast` ⇒
**explicit**-eligible. The veg-energy heat-capacity routing (§2.2) is the first concrete instance
(`τ = C/g`). A one-off **Jacobian-spectrum probe** (per-block eigenvalue/timescale scan over the
benchmark states) places the borderline components (wood, CAS-CO₂, deep soil) empirically before P1
freezes the partition; the same `τ_i`-vs-`dt` test then drives the per-cohort veg-energy routing at
runtime. This is open decision #2, and it is indeed the same idea as the hcap routing — generalized.

---

## 3. Scheme families to ship

Four families, two per approach, over the one RHS. Each has a "done right" condition without which it does
**not** deliver its advertised order — these are the substance of the build.

### 3.1 SPLIT-Picard (have it)

Sequential BE sweep iterated to a coupled fixed point. Order 1, coupling-consistent, warp-friendly (no
global Newton). Already implemented; keep as the robust baseline and the GPU candidate.

### 3.2 SPLIT-Strang + TR-BDF2 (new)

Symmetric ½·full·½ composition of the stiff sub-operators. **Gotcha: Strang over backward-Euler pieces is
still globally first-order** — the O(`dt²`) composition error is masked by the O(`dt`) error of BE
sub-solves. To realize 2nd order, the stiff sub-solvers must be **≥2nd order**, and for stiff surface
energy that means **TR-BDF2** (L-stable), **not** Crank–Nicolson (CN rings on stiff modes). Hydraulics is
already exact-exp (fine). So "SPLIT-Strang done right" = Strang composition + TR-BDF2 sub-solvers +
exact-exp hydraulics. Modular, warp-friendly, 2nd order, but retains O(`dt²`) splitting error on the
coupling.

### 3.3 IMEX-ARK2 (have it; upgrade)

Current ARS(2,2,2). Optionally swap to **Kennedy–Carpenter ARK3(2)4L[2]SA** (the design's named target)
for a better-matched embedded pair. Order 2, embedded error, one Newton/stage (the arrowhead), no
splitting error on the additive partition — **once §3.4 lands**.

### 3.4 IMEX-ARK "RK45" (new) — the high-order member

Your "adaptive RK45 for the non-stiff part" realized correctly. **Gotcha: this is a matched additive pair,
not a nested Dormand–Prince.** The explicit and implicit tableaux must share abscissae `c_i` and satisfy
the coupling order conditions; the implicit half must carry the stability. The right object is
**Kennedy–Carpenter ARK4(3)6L[2]SA**: an L-stable stiffly-accurate ESDIRK implicit half (one Newton/stage,
reuses the arrowhead) + a matched ERK explicit half + a 3rd-order embedded estimate driving
`adaptive_ark_march`. Running an *independent* adaptive RK45 on the explicit block would trigger
stability-limited step *rejection* (not accuracy control) the moment any stiffness leaks into it — the
classic IMEX blow-up. Keep the explicit set strictly non-stiff (§2.1).

**Prerequisite for 3.3 and 3.4 to pay off:** bring hydraulics (exponential stage-map) and soil-water
(implicit additive component) *into* the tableau, so the embedded error and the order reflect the coupled
state, not just CO₂.

### 3.5 Comparison (done-right column)

| Family | Order | Splitting error | Coupled solve/step | Warp-uniformity | Status |
|---|---|---|---|---|---|
| SPLIT-Picard | 1 (coupling-consistent) | none (iterated) | fixed-point (no Newton) | good | have |
| SPLIT-Strang + TR-BDF2 | 2 | O(`dt²`) on partition | none | good | new |
| IMEX-ARK2 (ARS/ARK3(2)) | 2 | none | 1 Newton/stage | poor | have → upgrade |
| IMEX-ARK4(3) | 4 (3rd embedded) | none | 1 Newton/stage | poor | new |

These four genuinely bracket the space — modular-warp-friendly-lower-order (SPLIT) vs
coupled-higher-order-adaptive (ARK) — which is precisely what goals (b)/(c) must measure and choose
between.

---

## 4. Goal (a): faithful error-control infrastructure

A uniform machinery shared by all four families, exposed through config as **strictness levels**:

- **L0 — fixed-step, no control** (GPU / production floor): fixed `dt_fast`, fixed inner-solver counts
  (warp-uniform). Cheapest, predictable cost.
- **L1 — embedded-error adaptive** (accuracy): per-scheme embedded local-error estimate + WRMS norm +
  step controller. This is the default "accurate" mode.
- **L2 — strict / faithful** (validation / debug): L1 **plus enforced conservation** — the whole-column
  energy/water/carbon residuals become **hard asserts**, and non-convergence is a hard failure.

Concrete pieces (most primitives already exist — this is assembly + wiring):

1. **One embedded-error path per scheme.** ARK2/ARK4(3) carry native embedded estimates; SPLIT-Strang gets
   step-doubling; the RK4 oracle stays the small-`dt` ground truth. Reuse `state_wrms`
   (`meds_fast_ark.f90`) and `adaptive_step_update` (`meds_numerics.f90`); consider a PI controller for
   smoother step sequences.

2. **A tolerance API with per-state-group `rtol`/`atol`.** Temperature (K), moisture (m³/m³), potentials
   (MPa), CO₂ (µmol/mol) live on different scales; one global tol is wrong. Define named groups and thread
   them consistently through every scheme's WRMS norm. Unifies the currently-scattered
   `ark_rtol` / soil `rtol,atol` / energy `rtol,atol` knobs.

3. **Conservation as a first-class, enforceable check.** `budget_check_stop` / `budget_assert`
   (`meds_budget_check.f90`) are **defined but never called anywhere in `src/`**, and `track_resid` is
   private to `meds_fast_split.f90`. Promote the whole-column residuals to a shared, callable check that
   L2 turns into a hard failure; add a **zero-pivot guard** to `thomas_solve`. (QW2 in §7 — this is the
   safety net every later refactor rides on.)

4. **The RK4 oracle is the definition of "correct."** Extend its coverage (it already marches the full
   `column_derivs`; add process-masked and reduced-system variants, §5) so every production scheme has a
   ground-truth reference at small `dt`, plus a Richardson-extrapolated reference where the oracle is too
   expensive.

**Deliverable:** an integrator that can be asked for a *target accuracy* or a *target conservation
tolerance* and either meets it (adapting) or fails loudly (L2) — never silently force-accepts. Today the
soil-water / hydraulics solvers silently force-accept at their substep cap; L2 forbids that.

---

## 5. Goal (b): the scheme × process × `dt` benchmark harness

A repeatable, CTest-able harness (not a one-off script) that sweeps the cross-product and emits structured
data.

### 5.1 Process-complexity toggles (a uniform process mask)

The harness must reduce `column_derivs` to sub-systems by turning processes on/off:
`{leaf/wood energy, soil heat, soil water, plant hydraulics, CAS-CO₂, snow}`. Some toggles exist in config
today (energy phase, multilayer roots, snow, advect-soil-heat); this generalizes them into **one process
mask** that both the RHS and every scheme honor, so the *same* harness runs a reduced ODE. Turning off
energy dynamics (your example) leaves a water+CO₂ column; turning off everything but hydraulics isolates
the stiffest mode. This is also **scientifically** useful (attribute behavior to a process) — not just a
numerics fixture.

### 5.2 Sweep axes

- **Scheme:** {SPLIT-Picard, SPLIT-Strang, ARK2, ARK4(3)} × {fixed, adaptive}.
- **`dt_fast`:** e.g. {1800, 900, 600, 300, 120, 60} s.
- **Process mask:** the 2ᵏ (or a curated subset) of §5.1.
- **Forcing regime:** the `MEDS_INTEGRATOR_TEST` cold-Jan / warm-Jun windows + a wet/saturating window
  (the ARK's historical stress case) + a high-VPD dry-down.
- **Hardware:** CPU-serial, CPU-multicore, GPU (once the offload track, §7/BB, lands).

### 5.3 Metrics emitted per cell

- **Accuracy:** vs the RK4/Richardson reference — hourly CAS-T RMSE, month-GPP bias, LE/H RMSE, ψ error.
- **Work:** f-evals, Newton iterations, step rejections, accepted substeps, expm evals.
- **Conservation:** the whole-column energy/water/carbon residuals (already computed).
- **Cost:** wall-clock per simulated month, CPU and GPU.

Emit as a tidy table (CSV + netCDF) keyed by (scheme, `dt`, mask, regime, hardware) — this table **is** the
input to goal (c). Reuse `runs/ithaca_ark30/integ/` (`gen_configs.py` / `run_matrix.sh` / `analyze.py`)
as the skeleton; generalize it to the process mask and the four families.

---

## 6. Goal (c): the config-driven integrator-selection rule

- **Inputs:** `dt_fast`, active process mask, target (accuracy tol / strictness level), hardware.
- **Output:** `{scheme family, tableau/order, adaptive|fixed, substep budget, inner-solver counts}`.
- **Method:** derive the rule from the §5 Pareto fronts. Start with a **hand-authored decision table**
  read off the benchmark data; evolve toward a small tabulated / fitted rule as data accumulates.
  Encode as a pure `select_integrator(cfg)` in the fast driver + a documented rule table that the
  benchmark **regenerates** (so the rule is reproducible from data, not folklore).

Illustrative shape of the expected rule (to be *confirmed by data*, not assumed):

| Situation | Likely choice |
|---|---|
| GPU, energy-on, large `dt` | SPLIT-Strang fixed-count **or** warp-uniform ARK2-fixed (warp-uniform wins) |
| CPU accuracy run, full physics | ARK-adaptive (ARK2 at production `dt`; ARK4(3) for tight tol) |
| Process-reduced / non-stiff mask | explicit-heavy: adaptive RK45 half dominates, cheap |
| Strict validation (L2) | ARK-adaptive + enforced conservation, small `dt` |
| Coarse `dt`, energy/water-dominated | SPLIT-adaptive (empirically most accurate + ~30 % cheaper — §1.3) |

The existing 24-run result (SPLIT competitive and cheaper at 900 s; ARK-fixed-1 pathological) is the first
row of this table — goal (c) is that observation, generalized across the process mask and hardware.

---

## 7. Supporting pathways (folded from the broader scoping)

**Deferred-by-design**
- **`beta_stomata` stays inert (`env%psi_soil = 0`).** The stomatal drought limb must read **predawn / prior-day
  maximum** ψ (a *daily* quantity), **not** instantaneous or fast-loop ψ — so it does **not** wire into the
  fast loop; when it lands it enters as a daily frozen parameter (Category 0). Tracked in an existing
  GitHub issue; out of scope here.

**Quick wins (enable goal a; mostly bit-identical)**
- **QW2 — wire the conservation halts + `thomas_solve` zero-pivot guard.** The L2 safety net (§4.3).
- **QW3 — daily biogeochem on the exact EXPM path** (already implemented for spin-up): simultaneously more
  accurate, unconditionally stable, and **positivity-preserving** vs Euler. Breaks bit-identical.
- **QW4 — one shared, warp-bounded `expm` in `meds_numerics`** (promote the biogeochem scaling-and-squaring;
  keep the analytic 2×2 sinhc as a specialization). Removes a duplicated, divergence-prone primitive.
- **QW5 — transpose ψ to coalesced per-node SoA** (`psi_leaf/psi_wood/psi_root` vs `psi(3,ncoh)`); update
  the single lockstep reorder. Folds into BB1.

**Medium bets**
- **MB2 — soil-energy substepping.** `flux%nsub=1` is hardcoded and the `energy_opts_t` adaptive knobs are
  never read; implement them (a fixed warp-uniform path mirroring `soil_water`). Prerequisite for
  trustworthy phase change and for uniform error control across stores.
- **MB3 — Strang leaf↔CAS seam.** Part of family §3.2.
- **MB4 — lift the ARK 2×2 Newton into the split path** as the leaf↔CAS coupler (guaranteed quadratic
  convergence vs `relax=0.5` Picard). Shared-arrowhead reuse.

**State-layout + GPU track — BB1 and BB2 are COMMITTED to this development** (BB3 follows). BB1 is not
merely a GPU enabler: it is a bit-identical structural refactor that **removes the per-patch
gather-into-scratch** and gives real **code clarity** (one global state, direct CSR indexing), and it is
the single state layout that `vegetation_energy_time_derivs` (§2.2) and QW5 (ψ-SoA) both build on — so it
is co-designed with them, not bolted on later.
- **BB1 — hoist the fast prognostic state into global device-resident flat SoA + one `!$omp target data`
  day-region** (already a named CLAUDE.md follow-up). The fast prognostic set = CAS twins, soil E/θ,
  per-cohort ψ, per-cohort leaf/wood energy, snow — sized once to (n_patch)/(n_cohort_global), keyed by the
  existing CSR (`cohort_offset`/`cohort_count`/`owner_patch`); replaces per-patch alloc+gather with direct
  global indexing. **Land bit-identical and still serial first**, guarded by the whole-column budgets and
  the `tc_split(54)=292.450065` golden anchor. Root of the GPU DAG; likely double-digit CPU win + code
  clarity on its own.
- **BB2 — offload `column_fast_step` as two segmented kernels** (per-cohort over all cohorts — leaf gas
  exchange, hydraulics, leaf/wood energy; per-patch — CAS box + soil Thomas), coupled by a deterministic
  segmented reduction over the CSR — mirroring the proven `meds_core_state_update` offload (OpenMP target +
  `-gpu=mem:separate`, **never** `-stdpar=gpu`; build the nvfortran multicore back end to catch the
  issue-#7 array-temp trap). Do it once the scheme set is stable so the offload target is fixed.
- **BB3 (follows) — warp-synchronous team-collective adaptivity + stiffness binning** — one shared step
  `h=min` over a warp's lanes driven by team-max WRMS; reconciles adaptivity with warp-uniformity. Plus a
  `gpu_warp_uniform` profile that defaults the already-coded fixed-count solvers (`SOIL_SUBSTEP_FIXED`,
  `HYDRO_SUBSTEP_FIXED`, `ark_fixed_substep`, `picard_fixed_iter`).

---

## 8. Phasing

Dependency-ordered. `[RESUME]` marks work already designed/deferred in existing docs (pick up, don't
redesign).

- **P0 — Foundation for (a) + (b) + the state layout.** RHS partition cleanup + the uniform **process
  mask** (§5.1); the **tolerance/error-control API** with per-group `rtol/atol` (§4.1–4.2); **wire the
  conservation halts + `thomas_solve` guard** (QW2); **extend the RK4 oracle** to process-masked / reduced
  systems (§4.4); **BB1** — hoist the fast prognostic state into global device-resident flat SoA,
  **bit-identical and serial**, the home for QW5 (ψ-SoA) and the new leaf/wood energy state and the
  co-design point for `vegetation_energy_time_derivs`. Land QW4 here too.
- **P1 — The comparands.** Bring **hydraulics + soil-water into the additive scheme** (§2.2, §3.4); add
  **`vegetation_energy_time_derivs`** so prognostic leaf/wood integrate under every scheme (removes the ARK
  `error_stop`; `[RESUME]` the leaf-wood-energy branch); **ARK2 proper** (optionally ARK3(2));
  **SPLIT-Strang + TR-BDF2** (§3.2). MB2 (soil-energy substepping), MB4 (Newton leaf↔CAS). `[RESUME]` the
  ARK follow-ons in `MEDS_IMEX_ARK_DESIGN.md`.
- **P2 — The harness (goal b).** Assemble the sweep + metrics + structured emission (§5), generalizing
  `runs/ithaca_ark30/integ/`. Produce the first full scheme × mask (incl. diagnostic↔prognostic veg
  energy) × `dt` dataset.
- **P3 — High-order + strict mode + offload.** ARK4(3)6L[2]SA (§3.4); the L2 strict-conservation mode
  (§4); QW3 (EXPM daily biogeochem, positivity); **BB2** — offload `column_fast_step` as two segmented
  kernels (scheme set now stable → fixed offload target).
- **P4 — The selection rule (goal c) + BB3.** Fit `select_integrator(cfg)` from the P2 data (incl. the
  hardware axis); document the regenerable rule table (§6); **BB3** team-collective warp-synchronous
  adaptivity + the `gpu_warp_uniform` fixed-count profile.

### 8a. Implementation status (2026-07-22, this pass)

A first slice of P0 was implemented and verified (ifx Debug `-check all` + nvfortran multicore, full
36-test suite green on both, plus a 30-year multi-patch/multi-cohort `meds_main` run as an additional
stress test — "area conserved, no NaNs"). Uncommitted on `main`, no PR yet.

**Done:**
- **QW4 — shared bounded `expm`.** `matrix_exp`/`matrix_exp_fixed`/`matmul_sq` promoted from
  `meds_soil_biogeochem.f90` to `meds_numerics.f90` (dimension-free, issue-#7-safe); `matrix_exp`
  delegates to the new fixed-squaring-count `matrix_exp_fixed` (the GPU/warp-uniform sibling described
  in §7/QW4). `soil_carbon_step`'s `step_expm` now calls the shared primitive — bit-identical (verified:
  identical scaling-and-squaring arithmetic, just relocated).
- **QW2 — enforceable conservation + zero-pivot guard.** `budget_check_stop` wired after every
  `budget_accumulate`/`track_resid` call in both `meds_fast_split.f90` (7 budgets) and
  `meds_fast_ark.f90` (7 budgets), gated on the existing `ccfg%energy%debug_error` flag (default
  `.false.` ⇒ **no behavior change**; this is the L0→L2 enforcement switch from §4.3). `thomas_solve`
  gained a sign-preserving zero-pivot floor (mirrors `meds_soil_biogeochem`'s `gaussian_solve` style);
  a no-op for the diagonally-dominant matrices it solves today.
  **Finding (not fixed, flagged for follow-up):** stress-testing the new hard-stop with
  `debug_error=.true.` forced across the whole suite surfaced that `soil_water (split)`'s own residual
  (`hflux%mass_resid`) already breaches its stated tolerance (`rtol=1e-6, atol=1e-4`) in `test_fast_loop`
  (resid 2.35e-3) and `test_biogeochem_dynamics` (resid 2.7e-1) — a **pre-existing** silent gap (the
  tolerance literals are unchanged from before this pass; `n_fail` was already counting it, just never
  enforced). This is exactly cross-cutting weakness #4 from §3, now caught for real. Root cause not
  investigated (candidate: the Richards solver force-accepting at its substep cap under these forcing
  conditions) — left as a follow-up, since fixing it is a numerics investigation outside QW2's scope.
- **BB1 phase 1 — allocation-churn removal (NOT the full device-resident-SoA vision).** Added
  grow-only "ensure capacity" allocators (`ensure_column_cohort_capacity`, `ensure_patch_biophys_capacity`,
  `ensure_aero_out_capacity`, plus a grow-only fix to the private `alloc_forcing`) and pre-size `coh`/
  `bio`/`aero`/`forc` + the per-cohort output accumulators to the **site-wide max cohort count once per
  `fast_dynamics` call**, instead of once **per patch**. Cuts O(n_patch) heap allocations/call to O(1).
  Verified safe: no code anywhere reads `size(coh%pft)`/`size(bio%leaf_temp)`/etc. (only the active count
  `coh%n`/`ncoh` is ever used as a loop bound), and the four now-oversized-capacity accumulator arrays
  (`gpp_coh` etc., plain arrays with no own active-count field) are explicitly sliced `(1:ncoh)` at their
  two call sites to stay assumed-shape-correct regardless of backing capacity.
  **What this is NOT:** the persistent prognostic reservoirs (CAS/soil/snow per patch, leaf_temp/
  wood_temp/psi per cohort) were **already** site-wide flat SoA before this pass — not per-patch scratch
  — so BB1's "hoist persistent state" motivation was already largely satisfied by the existing
  architecture. What remained un-hoisted was the *transient* per-patch working set
  (`column_cohort_t`/`patch_biophys_t`/`column_forcing_t`/`aero_out_t`) — phase 1 stopped reallocating it
  but did not restructure it, and (as phase 2 below found) that restructuring is not what actually
  unlocks GPU offload anyway.

- **BB1 phase 2 — bare-array proof of concept (a re-scoped target, see the note below).** Re-examining
  phase 1's own "not `!$omp target`-mappable" caveat found that **pointer-aliasing into the global SoA
  would not have fixed it**: `column_fast_step` still takes derived types with allocatable/pointer
  components, and OpenMP `target` regions can't cleanly map those — the codebase's own working GPU
  kernel (`meds_core_state_update`'s `update_cohort_states`) is offload-eligible specifically *because*
  it takes bare arrays (CLAUDE.md: "no `site_t`, no derived types, so the `map` clauses are clean").
  The real prerequisite is a **bare-array calling convention** for the fast loop's per-cohort kernels —
  which is really BB2's foundational step pulled forward, not a data-layout tweak. Implemented a
  proof-of-concept on the cleanest, most independent candidate: **plant hydraulics**. Added
  **`solve_plant_water_batch`** (`meds_plant_hydraulics.f90`) — a bare-array wrapper that loops
  `do i=1,n` calling the *existing, unmodified* `solve_plant_water` once per cohort (zero changes to
  the validated per-cohort physics); every argument is a bare array or scalar broadcast (mirroring
  `update_cohort_states`'s pattern), re-exported through `meds_plant_interface`. Rewired
  `meds_fast_split.f90`'s hydraulics call site: precompute the per-(layer,cohort) rhizosphere
  conductance into an array (only exercised when `multilayer_roots`), compute the per-plant
  transpiration array, one call to the batch routine for the whole patch, then the same postprocess
  accumulation loop as before (identical `i`-outer/`k`-inner order, so the running `root_sink_share`
  sum is bit-identical). No `!$omp target` pragma added yet — `solve_plant_water` is not `pure` (it can
  `error stop` on a stale vulnerability-table guard), which blocks literal device offload; that follow-up
  (making the guard a caller-side precondition instead) is left for the dedicated BB2 pass, consistent
  with "land it bit-identical and still serial first."
  **Verification (this was the highest-risk change of the session — new 2D-array indexing with no
  existing test combining `multilayer_roots=.true.` with >1 cohort):** full 36-test suite green on ifx
  Debug + nvfortran multicore; then, since no unit test exercised the multilayer + multi-cohort
  combination together, ran a **before/after `git stash` comparison** — a 30-year multi-patch,
  multi-cohort, `multilayer_roots=true` `meds_main` run against the pre-refactor code, diagnostic netCDF
  output compared with `cmp` — **byte-for-byte identical**. Repeated for the default (non-multilayer)
  path against the phase-1 checkpoint — also byte-for-byte identical.
  **What remains for a real BB2:** extending the bare-array treatment to the rest of `column_fast_step`
  (leaf gas exchange, CAS box, soil heat/water — far more tangled with cross-store coupling than
  hydraulics' clean per-cohort independence), making the per-cohort kernels `pure` where they currently
  aren't, and only then adding actual `!$omp target` pragmas + validating on the `build-gpu` back end.

**Deferred (not started this pass), with rationale:**
- **Process mask + RK4-oracle extension (§5.1).** Substantial new type + threading through
  `column_derivs`/`surface_derivs` and the oracle; several open design choices (exactly which
  granularity, whether to also thread it into the production `column_fast_step`/`column_fast_step_ark`)
  are better resolved with the user before committing an implementation.
- **Tolerance/PI-controller module (`meds_fast_control.f90`, §4/§9.3).** ~~Deferred~~ **DONE in a
  follow-up pass — see §8b (goal a).**
- **`column_prepass` relocation (§10.2).** On inspection this is less "purely mechanical" than scoped:
  `column_prepass` pulls in `meds_config_t`, leaf gas-exchange, canopy aerodynamics, and heterotrophic
  respiration — moving it would add several new dependency edges to whatever module receives it (the
  plan's own §10.2 suggests a *new* small module, `meds_fast_frozen.f90`, rather than folding it into
  `meds_fast_time_derivs.f90`, to avoid heavying up the pure-RHS module). Deferred to keep this pass's
  diff focused on behavior-preserving numerics work rather than a pure reorganization.
- **Jacobian-spectrum stiffness probe (§2.3/§9.2).** Not started.

### 8b. Goal (a) — error-control infrastructure (2026-07-23, IMPLEMENTED)

The tolerance/PI/strictness machinery of §4 + §9.3, delivered as one module and wired into the ARK
adaptive march. Verified 36/36 on ifx Debug + nvfortran multicore; **byte-identical to the committed
baseline with the ARK path ON** (`time_integrator="ark"`, `fast_biophysics_on=true`, 2yr `meds_main`,
netCDF `cmp`) under the default controller — so it is a pure capability add, no behavior change off by
default.

- **New `src/driver/meds_fast_control.f90`:** `tol_set_t` (per-**group** rtol/atol over the 5 physical
  field classes — CAS enthalpy / shv / CO₂ / soil energy / ψ), `error_control_t` (strictness level +
  controller + step-clamp knobs + PI gains + the tol set), `state_wrms_grouped` (the grouped WRMS —
  generalizes the old single-`rtol` `state_wrms`), `step_control_factor` (dispatch: **I-controller** =
  the legacy `adaptive_step_update`, or **PI** = Gustafsson `safety·err^-α·err_prev^β`), and
  `default_tol_set`/`default_error_control` (broadcast one rtol → byte-identical to the old norm).
- **Selectors in `meds_config`** (`CTRL_I`/`CTRL_PI`, `CTRL_L0/L1/L2`) + config fields
  `[fast].step_controller` (`i`|`pi`, default `i`) and `[fast].error_level` (`fixed`|`adaptive`|`strict`,
  default `adaptive`), read defaulted in `meds_config_io` (existing TOMLs unchanged) + range-validated.
- **Wiring:** `adaptive_ark_march` now takes an `error_control_t` (was a scalar `rtol`), uses
  `state_wrms_grouped` + `step_control_factor` with `err_prev` tracking for PI, and — under
  **L2 strict** — `error stop`s on a floor-forced accept that still breaches tolerance (vs L1's graceful
  force-accept, the "faithful/validation" mode). The old `state_wrms`/`ATOL_*` were removed from
  `meds_fast_ark` and the RK4 oracle repointed to `state_wrms_grouped`.
- **PI proven functional:** on `test_column_derivs`'s stiff 24 h ARK march, the I-controller takes
  **19** substeps and PI takes **28** (both bounded, tcas ≈ 307.6 K) — different step sequences, so the
  PI path is genuinely engaged (the smooth constant-forcing whole-model MVP is single-substep per
  `dt_fast`, so there PI ≡ I, as expected).

**Still open in goal (a):** the L2 strict mode ties into QW2's already-wired `budget_check_stop` (a full
"enforce conservation everywhere" sweep is the remaining L2 piece). The scattered sub-solver tolerances
were closed by §8c below.

### 8c. Goal (a) Layer 1 — sub-solver tolerance unification (2026-07-23, IMPLEMENTED)

The fast loop carries **four** independent adaptive error estimates — the ARK march plus the nested
soil-water, soil-energy and plant-hydraulics step-doubling solvers — and each owned its own tolerance.
Two were unreachable from config (`hydraulics_config_t` has no rtol/atol at all, so the hydraulics
solver ran on a *type default*), so "the run's target accuracy" was not a well-defined quantity.

- `tol_set_t` gains **`GRP_THETA`** (soil moisture) and **`GRP_SOIL_T`** (soil temperature), so all 7
  groups now cover every adaptive estimate.
- **`build_tol_set(cfg)`** is the single authority, seeding each group from whatever governs it *today*
  (ARK groups ← `[fast].ark_rtol`; `GRP_THETA` ← `[soil]`; `GRP_SOIL_T` ← `[energy]`; `GRP_PSI` shared by
  the ARK WRMS and the hydraulics solver, whose defaults already agreed) ⇒ the default is the identity.
- **`build_fast_context` pushes the resolved groups down** into `ccfg%hydro` / `%energy` / `%hydro_o` —
  which is what finally makes the hydraulics tolerance settable at all.
- **Two master dials, and both are needed.** `[fast].rtol_all` overrides every group's rtol;
  `[fast].atol_scale` multiplies every group's atol. Because the WRMS denominator is `atol + rtol·|y|`,
  `rtol_all` **alone saturates** once atol dominates: measured on the split path, 1e-3 → 1e-6 raises the
  soil-water error estimate only ~4× (4e-4 → 1e-4 denominator), never enough to flip a substep decision.
  Scaling both does bite. Defaults (`0` / `1.0`) are exact identities.

**Verification (and a corrected instrument).** Byte-identical on both the split and ARK paths — but the
check had to be re-based first: the legacy `[io]` **`-D-output.nc` stream is demographic-only**
(nplant/dbh/AGB/…), carries no fast-loop state, and on a near-bare fixture is very nearly blind to the
integrator. It reported SPLIT ≡ ARK, which is impossible for two different algorithms. Fast-loop
verification now uses the **`[output]` daily aggregation stream** (`cas_temp_site`, `soil_temp_top_site`,
`et_site`, `gpp_site`), which does separate the schemes. `test_column_dynamics` **RUN 5** asserts the
plumbing directly — a whole-model A/B cannot see a tolerance change that does not flip a substep
decision — and was mutation-checked by disabling the hydraulics push-down (it fails).

### 8d. §5.1 process mask + §5 harness (2026-07-23, IMPLEMENTED)

- **`process_mask_t`** on `column_config_t` (so it reaches both schemes with no signature churn): seven
  switches — `veg_energy`, `cas_energy`, `cas_vapour`, `cas_co2`, `soil_heat`, `soil_water`,
  `hydraulics`, config-settable as `[fast.mask]`.
- **Semantics: frozen, not skipped.** A masked-off store is restored to state^n after its kernel runs, so
  it contributes no ODE dimension while still supplying its couplings as a constant. The freeze reuses the
  Picard `soil_w_n`/`soil_e_n`/`psi_n` snapshots already being taken (no new buffers, no alloc churn); the
  ARK applies the same freeze at its single state-commit point. **Because the kernel still runs, a reduced
  column is not cheaper** — which matters when reading harness *work* metrics.
- **Conservation:** a reduced column cannot close by construction, so `mask_is_full()` gates the budget
  **hard stops** (the soft `n_fail` counters still tally). All-on ⇒ exactly the old `debug_error` gate.
- **`scripts/numerics_sweep.py`** sweeps scheme × `dt_fast` × mask × controller × `rtol_all`, runs each
  cell, and emits one tidy `sweep.csv`. Two design points: configs are generated by **parsing and
  re-serializing** the base TOML (appending an override block to a config that already defines that table
  silently re-parents every key that follows — a duplicate `[soil]` is *not* a merge, and the symptom is
  indistinguishable from "the setting had no effect"); and there is **one reference per mask**, since a
  masked cell integrates a different physical system.

**Verified:** all-on byte-identical on both paths; on a census stand (LAI ≈ 7, real GPP) six switches
demonstrably reduce the system, and `veg_energy` is correctly a no-op under the default *diagnostic*
leaf/wood (no store exists) while biting under prognostic wood. The harness caught the mask being honored
only by the split stepper. ifx 36/36 + nvfortran multicore 36/36 throughout.

**First dataset** (census stand, Jan, split/ark × {1800, 900} s × {full, no_energy}): split-full converges
**~1st order** in `dt` (5.2e-5 → 2.8e-5 K RMSE on CAS T), ARK-full **~1.6 order** (1.5e-4 → 5.0e-5 K) and
is **less accurate than the split at equal `dt`**. This is consistent with §3.5's note that the ARK's
formal order collapses toward p1 because hydraulics + soil water are operator-split **out** of the
tableau — i.e. the evidence says the next P1 step is **removing that splitting-order barrier, not adding
a higher-order tableau** (ARK4(3) over a p1-limited splitting would buy nothing).

**Honest gap:** §5.3's *work* counters (f-evals, Newton iterations, rejections, accepted substeps) live in
`column_budget_t` and are not serialized to netCDF, so wall-clock is the only cost axis the harness
reports. Emitting them is the natural next harness increment.

### 8e. Instrument repair + the corrected measurement (2026-07-23)

A multi-agent investigation into "should hydraulics + soil water go into the tableau?" returned **no**
for hydraulics (a category error: ψ is already a state, advanced by an exact exponential map, and
`surface_derivs` never reads it — the coupled subsystem is exactly ψ-independent over a step) and
"not for order" for soil water. It also found that the §8d numbers were produced by a **broken
instrument**. Three repairs landed first:

1. **Work counters (§5.3) — the missing cost axis.** `column_budget_t` now carries `integ_nsteps` /
   `integ_nrej` / `soil_nsub` / `hydro_nsub` / `hydro_nonconv`; **both** steppers fill them at the same
   seams; the driver area-weights them onto `site_t`; an opt-in `GRP_NUMERICS` group (`[output].numerics`)
   emits five `work_*_site` variables. Rejected steps are counted, not just accepted ones.
2. **The ARK error norm was diluted ~1.4×.** `state_err_diff` zeroes `err%psi`, so `y_lo%psi == y_new%psi`
   exactly, yet `state_wrms_grouped` summed 2n structurally-zero terms **and counted them**. The march
   now passes `with_psi=.false.`; the RK4 oracle keeps ψ (it compares two fully-evolved states, where ψ
   genuinely differs). Pre-existing — the byte-identical Layer-1 refactor preserved it faithfully, which
   is exactly why it went unnoticed.
3. **Process-mask asymmetry.** On the ARK path `mask%soil_water` restored only θ, then
   `w_surface`/`w_aquifer`/`z_wt` were committed unconditionally on the next lines, so the two schemes
   ran *different* reduced systems. Fixed.

**Corrected measurement** (forced ERA5, 10 d, census stand, reference = split @ 54 s, min-of-5 wall with
the 0.065 s fixed overhead subtracted). This **supersedes §8d's ranking**: with the norm repaired the ARK
is markedly better than previously reported.

| scheme | dt | RMSE CAS-T | RMSE soil-T | net s | steps | rej |
|---|---|---|---|---|---|---|
| ark | 450 | 0.0513 | 0.0812 | 0.603 | 2808 | 21.0% |
| split | 225 | 0.0612 | 0.0179 | 0.603 | 3840 | 0% |

At **equal measured cost** the ARK wins CAS-T by 1.19× and *loses* soil-T by 4.5×. Neither scheme shows a
clean order: successive-halving gives split ≈ 0.57/1.41/1.32 and ARK ≈ 1.55/1.19/1.53.

**The headline finding — a structural rejection storm.** Rejections per `dt_fast` are 1.65 / 1.08 / 0.39 /
0.014 at dt = 1800/900/450/225: at production `dt_fast` the ARK throws away **about one full step attempt
per call**. Cause is not controller gain — the PI controller only moves 31% → 27% and costs *more* steps.
It is `dt0 = dt_fast` (meds_fast_ark.f90): the march **cold-starts at the full step every call**,
discarding the step size it just converged to, and re-learns it ~960 times a run. Seeding a smaller
initial step via the existing `ark_dt_init` collapses the rejection rate at identical net wall time:

| `ark_dt_init` | rej % | net s |
|---|---|---|
| cold (= dt_fast) | 31.1 | 0.452 |
| 450 s | 17.0 | 0.502 |
| 225 s | 8.4 | 0.452 |
| 112 s | 3.6 | 0.452 |

**Next step: persist the converged step size across `dt_fast` calls (a warm start)** — self-tuning, unlike
the hand-set constant, and it recovers ~30% of the ARK's wasted work before any scheme redesign.

**Known confound, not yet closed:** accuracy is scored against a *split-family* reference, which plausibly
flatters the split (especially on soil-T, where it "wins" 4.5×). Wiring the RK4 oracle as a selectable
sweep reference (`--ref-scheme rk4`) is the fix and should precede any scheme verdict.

### 8f. The dominant error is the COEFFICIENT FREEZE, not the integrator (2026-07-23)

Two independent measurements now say the same thing, and together they reframe the whole effort.

**(i) The tableau has no headroom left at production `dt_fast`.** Holding `dt_fast` fixed and refining
ONLY the ESDIRK step (`ark_fixed_substep` 1 → 32) does not reduce error — constant forcing, dt = 1800 s:
1.557e-4 → plateau 2.34e-4 K; forced ERA5, dt = 900 s: flat at 1.01e-1 across nsub = 1..16 — while
halving `dt_fast` does (0.128 → 0.054). Refining the *inner* march 1000× (3915 → 45971 substeps, 12×
work) moves soil-T only 0.0753 → 0.0662 K.

**(ii) The freeze is INCONSISTENT, and that inconsistency is worth ~2×.** The fast loop samples met at
the sub-interval midpoint while `column_prepass` freezes every Category-0 coefficient on the state at
tⁿ. Midpoint is the better quadrature of the forcing *alone*, but the coefficients depend on forcing AND
state, so this is a first-order mismatch. Sweeping the new `[fast].forcing_sample_frac` (forced ERA5,
10 d, census stand, dt = 900 s, RMSE CAS-T, scored against BOTH a `frac=0.5` and a `frac=0.0` fine-dt
reference so the result is not an artifact of the reference's own sampling):

| frac | vs ref(0.5) | vs ref(0.0) |
|---|---|---|
| 0.00 | 0.2049 K | 0.2015 K |
| 0.25 | 0.2469 K | 0.2593 K |
| 0.50 | **0.4041 K** | **0.4216 K** |
| 1.00 | 0.9949 K | 1.0145 K |

Monotonic on CAS-T, minimum where forcing and state agree: the production default costs ~1.97× on CAS
temperature at this site/season.

**But this is a TRADE, not a free win, and the harness nearly hid it.** Scoring the other channels on the
same runs, CAS-T and soil-T have *opposite* optima:

| frac | RMSE CAS-T | RMSE soil-T | bias ET | bias GPP |
|---|---|---|---|---|
| 0.00 | **0.2049** | 0.1707 | 0.0172 | −0.0058 |
| 0.50 | 0.4041 | 0.1392 | 0.0192 | −0.0075 |
| 0.75 | 0.5948 | **0.1228** | 0.0205 | −0.0084 |

Setting `frac = 0` buys 1.97× on CAS-T and *loses* 23% on soil-T. **So `frac = 0` must NOT be adopted as
a default.** The theory says why: for a coefficient `C(y, F)` frozen over the step with forcing sampled at
offset α, the local error is `dt²·[C_F·Ḟ·(½ − α) + C_y·ẏ/2]`. A measured optimum at α ≈ 0 means
`C_y·ẏ ≈ −C_F·Ḟ` *at this site, season and stand* — a cancellation, not a structural optimum. Two further
signs it is not structural: at dt = 1800 s the consistent-at-t case is WORSE (0.7381 vs 0.6008), and soil-T
moves the opposite way from CAS-T across the whole sweep.

The structural fix is a **consistent midpoint freeze** (forcing *and* state at t + dt/2), which annihilates
the `C_F` term for any site and season. Its predictor must be **L-stable, not explicit**: the CAS enthalpy
relaxation time is median ~128 s over the high-conductance daytime half, so at dt_fast = 900 s
`x = dt/τ ≈ 3.5`, where a forward-Euler predictor has error `|1 − x − e^{−x}| = 2.53` — worse than using
yⁿ unchanged (0.92). A half-step of the existing closed-form `cas_column_step_implicit` gives 0.19.

**Caveat on the mechanism.** `forcing_sample_frac` moves the met sample for the coefficient pre-pass AND
for the boundary sources (`apply_rt_forcing` cosz/SW, atm T/q/CO₂) simultaneously. The 2× is a measured
property of the knob; attributing it specifically to the *coefficient* freeze is an inference until the
two roles are separated.

### 8g. The two schemes were not integrating the same model

`surface_derivs` applies a smooth CAS supersaturation (condensation) sink, and it is reached **only from
the ARK stages** — `meds_fast_split` never calls it. Every split-vs-ARK number produced before this was
therefore comparing two *models*, not two integrators. `[fast].cas_condensation` (default `.true.` =
unchanged) makes it controllable. At dt = 450 s the ARK-only sink accounts for ~22% of the soil-T gap
(0.0753 → 0.0591 K). Worse than an asymmetry: `bf%whole_wat_out` includes `sf%cond`, so **the condensed
dew/fog water leaves the column entirely** instead of being deposited on the soil / interception store —
a water sink on the ARK path, and the likelier reason disabling it moves layer-1 θ by 59% while barely
moving ET. So there are two defects, not one: the sink is missing from the split, and its destination is
wrong on the ARK. Related pre-existing bug found in the same pass: **the ARK path has no `ccfg%snow_on`
guard** — the dispatch returns before the split's snow block, so a snow-enabled ARK run silently drops
the snow store.

**Corrected soil-T picture.** At **equal `dt`** the gap is only ~1.65×, and it is mostly θ contaminating
the *temperature read-off* (`soil_temp_top = uext_to_temp(soil_energy, θ·ρ, C_dry)`). Removing that
contamination leaves **0.95–1.42×** — i.e. on soil *energy* the ARK is at parity or marginally better at
equal dt. (∂T/∂θ ≈ −310 K per unit θ, using the code's `tsupercool_liq` = 56.79 K from
`meds_constants.f90`; an earlier −128 K figure used `t_3ple − L_f/c_liq`, dropping the `cp_ice` term, and
every number derived from it was wrong by 2.4×.) The κ(θⁿ) hypothesis is **refuted** as the dominant cause.

At **equal cost** the gap is ≥4.2× — but it cannot currently be quantified: split@225's soil-T score
(0.0179 K) is indistinguishable from the reference-disagreement floor (0.0181 K). Decomposing 4.2× into
"1.65× × 2.54×" is a tautology (A/C = (A/B)(B/C) for any B) whose second factor rests on that
unresolvable number, and it should not be quoted as if the equal-cost gap were explained away. A second ARK-only defect remains open: the **transp↔uptake gap** (the scratch hydrology
removes the state-ⁿ demand while the stages re-evaluate it, and the difference is never returned to the
soil) — a probe closing it cut θ error 33% and soil-T 14%.

**Reference floor, and the cheap fix.** RMSE between the split@54 s and ark@54 s references is 0.0181 K
on soil-T and 0.0137 K on CAS-T — the same size as split@225's own score, so **any ranking at dt ≤ 225 s
was reference-limited**. Two distinct causes, and only one needs Fortran:

* *Reference truncation* — fixed by the EXISTING `--ref-refine` flag, not by new code. Measured
  self-convergence of the split family against a split@6 s anchor: split@54 carries 0.0168 K (CAS-T) /
  0.0065 K (soil-T); split@18 carries 0.0042 / 0.0013 K, for 3.63 s → 6.24 s. The harness default is now
  **12** (was 4). This retires most of the case for a bespoke RK4 reference.
* *Model asymmetry* — NOT fixed by refinement. The split@54-vs-ark@54 disagreement includes a layer-1 θ
  difference of 1.4e-4 that does **not** shrink with dt, because the ARK path carries the condensation
  sink (§8g) and the transp↔uptake gap that the split does not.

**The RK4 oracle cannot arbitrate this.** `column_derivs` calls `surface_derivs`
(`meds_fast_time_derivs.f90`), so the oracle inherits the same `TAU_COND` sink and sits inside the ARK's
model family; it also lacks the split's snow terms. Unify the model first, refine the reference second,
and only then consider a bespoke oracle reference.

**Warm start (landed).** The march cold-started at the full `dt_fast` every call and paid ~1 rejected
attempt per call (1.65/1.08/0.39/0.014 rejections per call at dt = 1800/900/450/225 s). Carrying the last
*accepted* step across calls — not the grown proposal, which merely re-imports the over-estimate — cut
rejections to 18.7/13.5/5.8% and net wall-clock 8–14%, accuracy unchanged.

**Statistical resolution — applies to every number in §8e–§8g.** All RMSEs are over n = 10 daily site
values from ONE 10-day July window at ONE site with one stand. No intervals were computed. Ratios below
~1.5× (the 1.21×/1.29×/1.42×/1.65× soil-T figures, and the 1.19× CAS-T "win") should not be treated as
resolved until they are repeated over a month and a second season, or given a bootstrap interval over the
daily values.

---

## 9. Decisions (resolved 2026-07-22)

1. **RESOLVED — ship all four families this round** (SPLIT{Picard, Strang+TR-BDF2} + IMEX-ARK{2, 4(3)}).
   P1 grows accordingly, but the shared RHS + stage core (§10) makes each additional family a tableau
   descriptor, not a new stepper — so the marginal cost per family is small.
2. **RESOLVED — do the stiffness auto-placement** (§2.3); it is the general form of the per-cohort
   heat-capacity veg-energy routing (§2.2). Land the Jacobian-spectrum probe in P0 to place the borderline
   components (wood, CAS-CO₂, deep soil) empirically; the same `τ < ε·dt` test drives the runtime
   veg-energy routing.
3. **Tolerance groups + controller (recommendation — confirm).** Two coupled choices in the adaptive
   controller:
   - **Tolerance groups.** The state vector mixes incomparable units — CAS T (~290 K), soil θ (0–0.5
     m³/m³), ψ (−2–0 MPa), CO₂ (~400 µmol/mol), internal energy (J). A single global tolerance is
     meaningless. Instead the local-error norm is a **WRMS**:
     `err = sqrt(mean_i( (Δy_i / (atol_g + rtol_g·|y_i|))² ))`, where each state `i` belongs to a physical
     **group** `g` with its own `(rtol_g, atol_g)` — a temperature group, a moisture group, a potential
     group, a CO₂ group, an energy group. This is standard CVODE/LSODA practice; it is what makes "target
     accuracy" well-defined across the coupled column. **Recommend:** define the 4–5 groups now and thread
     them through every scheme's norm (unifies today's scattered `ark_rtol` / soil `rtol,atol` / energy
     `rtol,atol`).
   - **Controller: PI vs the current I / step-doubling.** How the error becomes the next step. The current
     `adaptive_step_update` is an **I (integral) controller** `h·safety·(tol/err)^(1/(p+1))`, and the
     step-doubling paths pay 3× for the estimate. **Recommend:** use the schemes' **embedded** estimates
     (free — ARK2/ARK4(3)/TR-BDF2 all carry one) and upgrade to a **PI controller**
     `h·(tol/err_n)^{k_I}·(err_{n-1}/err_n)^{k_P}` (Gustafsson), which uses the previous error to damp the
     step-size "hunting" that an I-controller shows on stiff/rough problems → fewer rejections, smoother
     sequences. Keep step-doubling only as an oracle cross-check.
4. **RESOLVED — BB1 in P0, BB2 in P3** (timing accepted; no earlier partial offload). **Amended by
   implementation (2026-07-22):** BB1 and BB2 are less separable than this staging implied — a
   `!$omp target data` device-residency region is only meaningful once there is at least one
   `!$omp target` kernel inside it, and the thing that actually makes a kernel offload-eligible is a
   **bare-array calling convention**, not the SoA-hoist framing BB1 was scoped around (see §8a's
   "BB1 phase 2" note). BB1 phase 1 (allocation-churn removal) still stands as a clean, independently
   valuable P0 deliverable; the bare-array proof of concept (plant hydraulics) that phase 2 became is,
   honestly, BB2's foundational step pulled forward — future phasing should track it as such rather
   than as a distinct "BB1 phase 2."
5. **RESOLVED — add the sapflow enthalpy coupling in this PR** (§2.2). It is exactly the ED2 `qwflux_wl`
   sapflow term (Appendix B), MEDS already computes `flux%sapflow`, and only the enthalpy-advection term is
   missing — cheap, and the physically-correct leaf↔wood link.

---

## 10. Targeted file structure after the development

**Bias: consolidate and simplify first; add modules only where a genuinely new concept needs a home.**
Bit-identical is **not** a goal — much of the current fast-loop code is initial domain-specific drafting
(two hand-rolled full-column steppers with duplicated flux/budget logic), and this development is the
right time to collapse it. The correctness gates shift accordingly (see "Correctness strategy" below).

### 10.1 The consolidation thesis

Today the fast loop carries **two independent hand-written column steppers** — the operator-split body in
`meds_fast_split.f90` (~590 lines) and the ARK stepper in `meds_fast_ark.f90` (~1030 lines) — plus a
per-patch **gather-into-scratch** driver (`meds_fast_dynamics.f90`, ~660 lines) and per-patch scratch
types. They duplicate the CAS box, the soil solves, the hydraulics call, the budget ledger, and the frozen
pre-pass. The end state is **one** additive RHS, **one** set of stage builders, **one** integrator
dispatch hosting the four families as tableaux, **one** global state (no scratch), and **one**
error-control/conservation module. The four schemes become *data* (tableau descriptors), not code.

### 10.2 Target module map

| Module (src/driver unless noted) | Action | Contents / change |
|---|---|---|
| `meds_fast_time_derivs.f90` | **grow → the one RHS** | keep `surface_derivs`/`column_derivs`; add the veg-energy tendency wiring + the **process mask**; this is the single RHS every scheme calls |
| `meds_vegetation_biophysics.f90` (biophysics) | **extend** | add `vegetation_energy_time_derivs` (leaf+wood `dE/dt`) beside the existing diagnostic/BE kernels; add the **sapflow enthalpy** term (§2.2); no new file |
| `meds_fast_stage.f90` | **NEW (extract)** | the shared implicit stage builders pulled out of `meds_fast_ark`: `column_be_stage`, `newton_surface_solve`/`jac_surface` (arrowhead), the exponential hydraulics stage-map, the new **TR-BDF2** stage. Used by **both** split and ARK |
| `meds_fast_integrator.f90` | **NEW (consolidate)** | the one dispatch + the four **tableau descriptors** {split-picard, split-strang, ark2, ark4(3)} + the adaptive march. Absorbs the ARK-specific march wrapper and replaces the hand-rolled split body |
| `meds_fast_split.f90` | **shrink / absorb** | the sequential BE sweep becomes the single-stage (+Picard) tableau in `meds_fast_integrator`; the file is retired or reduced to the split/Picard descriptor |
| `meds_fast_ark.f90` | **shrink** | stage builders → `meds_fast_stage`; tableau/march → `meds_fast_integrator`; `column_prepass` moves out (below). Little unique logic remains |
| `column_prepass` (currently in `meds_fast_ark`) | **move** | it is the shared **Category-0 frozen pre-pass** (gs/GPP/Rd/radiation/aero) — belongs beside the RHS (`meds_fast_time_derivs` or a small `meds_fast_frozen`), not inside the ARK module |
| `meds_fast_control.f90` | **NEW (consolidate)** or extend `meds_numerics` | the error-control API: per-group WRMS tolerances, embedded-error interface, **PI controller**, the **enforced conservation** check (finally *calls* `budget_check_stop`). The soil-water and hydraulics substep loops call this instead of each rolling its own |
| `meds_fast_dynamics.f90` | **simplify (BB1)** | delete the per-patch `alloc_column_cohort`/`alloc_patch_biophys` + gather/scatter; iterate the **global fast-state SoA** by CSR index |
| `meds_fast_types.f90` | **extend** | the global fast-state SoA + the process-mask type; per-cohort leaf/wood energy + coalesced ψ (QW5) added to the state (one place, the lockstep discipline) |
| `meds_fast_rk4_oracle.f90` | **keep (test), slim** | `rk4_column_step` stays the ground-truth oracle; `imex_euler` reuses `meds_fast_stage`/`meds_fast_control` |
| `meds_numerics.f90` (shared) | **extend** | one **bounded `expm`** (QW4) — remove the biogeochem + hydraulics duplicates (keep the analytic 2×2 as a specialization); shared PI controller helper |
| `meds_soil_water.f90`, `meds_soil_energy.f90`, `meds_plant_hydraulics.f90` (biophysics) | **de-duplicate** | their private substep controllers call `meds_fast_control`; soil-energy substepping (MB2) implemented via the shared path; hydraulics fixed-count option |
| per-patch scratch types (`column_cohort_t`, gather targets) | **remove/thin** | replaced by SoA views once BB1 lands |
| `runs/ithaca_ark30/integ/` + a `test_fast_integrator_matrix` | **NEW (harness)** | generalize the existing scripts to scheme × process-mask × `dt`; emit the goal-(c) dataset; CTest target |

**Net:** two steppers + a gather driver + scattered controllers collapse to RHS + stage + integrator +
control + global-state — fewer lines, one authority per concept, four schemes as descriptors.

### 10.3 Correctness strategy (bit-identical relaxed)

Reimplementing the split as a tableau over the shared RHS will differ from today's hand-rolled body at
round-off — **acceptable**. The `tc_split(54)=292.450065` golden anchor is kept only as a cheap
regression *tripwire*, not a hard contract. The real gates become: **(1)** whole-column energy/water/carbon
**conservation closure** to machine precision (already in place; now *enforced*, §4-L2); **(2)**
cross-validation against the **RK4 oracle** at small `dt`; **(3)** benchmark **accuracy vs the reference**
(§5). Migration is safe because the trusted production path is SPLIT: build the new split-tableau, verify
it against the *old* split within a stated tolerance (plus gates 1–3), then retire the hand-rolled body.

---

## 11. Bare-array process-kernel convention (BB2 foundation + Python-wrap enabler)

**Goal (user-directed 2026-07-22):** every process kernel called in a per-element loop should have a
**bare-array batch entry point**, so (a) the loop becomes `!$omp target`-eligible (the codebase's only
working GPU kernel, `update_cohort_states`, is offloadable *because* it takes bare arrays — CLAUDE.md),
and (b) a Python/ctypes wrapper (`meds_plant_capi` pattern) can vectorise over numpy arrays instead of
calling one element at a time. This is the concrete continuation of §8a's "bare-array proof of concept"
— it is what a real BB2 is built on.

### 11.1 The convention — two tiers: `elemental pure` first, `*_batch` only where forced

There are **two** ways to give a kernel a bare-array interface; prefer the first.

**Tier 1 — `elemental pure` (no wrapper).** If a kernel's per-element inputs *and* outputs are all
**scalars**, and it can be **`pure`** (no external I/O, no reachable `error stop`, all callees pure), then
flatten its per-element `env`/`flux` derived types into **bare scalar arguments** (keep run-uniform POD
config like `wood_params_t` as a scalar dummy — it broadcasts) and mark it `elemental pure`. Fortran's
elemental broadcast then makes **one procedure serve both scalar callers (ARK RHS, tests) and array
callers (the fast loop, Python)** — a scalar call does one cohort, an array call does a whole patch. **No
`*_batch` wrapper, no derived-type env, no duplication.** This is strictly better than a wrapper for GPU
(bare SoA → coalesced) and Python (numpy-friendly), and it's what to do for every eligible kernel.

**Tier 2 — explicit `*_batch` (only when Tier 1 is impossible).** Two situations force a wrapper:
- **Per-element *vector/array* state** (`solve_plant_water`'s `psi(N_HYDRO)`, soil columns) — `elemental`
  forbids array args, so an explicit `do i=1,n` batch is unavoidable.
- **Heavy per-element config gather** — `leaf_gas_exchange` flattens the big allocatable-component
  `meds_config_t`/PFT table per PFT; passing `cfg` as an elemental broadcast dummy is dubious, and
  pure-ifying the whole 3-module leaf-solver chain is a large separate refactor. So leaf stays a batch.

For a Tier-2 batch: per-element quantities are bare arrays of length `n`; run/patch-uniform ones are
scalars/one POD (broadcast); outputs are bare arrays; the wrapper lives in the process module (its
bare-array public API), keeps the single-element kernel unchanged (so it's bit-identical to an inline
loop), and is named `<kernel>_batch`.

**Purity / `!$omp target` note.** Neither tier adds a `target` pragma yet. `elemental pure` is the
compiler's "embarrassingly parallel" signal (it can SIMD-vectorise / offload the broadcast); the Tier-2
batch is the map-able loop. Actual offload still needs the kernel `pure` (`solve_plant_water`'s three
`error stop`s → caller preconditions) + the pragma + `build-gpu` validation.

### 11.2 Verification protocol (corrects an earlier gap)

**The fast loop is gated by `[fast].fast_biophysics_on` (default `false`).** The 30-year `meds_main`
runs used for BB1/BB2 verification had it **off**, so they exercised only the slow/demographic path —
they did *not* test `column_prepass`/`column_fast_step` (the unit tests `test_fast_loop` /
`test_column_ark` / `test_column_dynamics` / `test_picard_coupling` did). Every fast-loop bare-array
conversion must therefore be verified with a **fast-loop-ON** run: `fast_biophysics_on = true`
(constant-forcing MVP is enough, `forcing_on = false`, no forcing file needed), a `git stash`
before/after `meds_main` run, diagnostic netCDF `cmp`'d for **byte-for-byte identity** vs `origin/main`
— in addition to the full 36-test suite on ifx Debug **and** nvfortran multicore (issue-#7).

### 11.3 Inventory + status

| Kernel(s) | Tier | Status | Notes |
|---|---|---|---|
| `stem_maintenance_respiration` | **1 (elemental)** | **✅ `elemental pure`; wrapper + `wood_env_t`/`wood_flux_t` DELETED** | env flattened to bare scalars; `wood_params_t` kept as scalar POD; test rewritten; byte-identical (fast-loop-ON) |
| `fine_root_maintenance_respiration` | **1 (elemental)** | **✅ `elemental pure`; wrapper + `root_env_t`/`root_flux_t` DELETED** | same |
| `veg_energy_diagnostic` | **1 (elemental)** | **✅ marked `elemental`** (was already `pure`, all-scalar) | enables array calls; the tangled *caller* loop (Picard, LAI-slave, wood branch, 5 running sums) is a separate restructure |
| `ground_surface_fluxes` | **1 (elemental)** | **✅ marked `elemental`** | all-scalar |
| `peaked_arrhenius_scale`, `meds_therm_lib` (17), `meds_hydr_lib` (16), `meds_allometry` (10) | 1 | **✅ already `elemental pure`** | the constitutive/thermo layer was done |
| `solve_plant_water` | **2 (batch — vector state)** | **✅ `solve_plant_water_batch`** | `psi(N_HYDRO)` per cohort ⇒ elemental impossible; §8a; byte-identical incl. multilayer-roots |
| `leaf_gas_exchange` | **2 (batch — heavy config gather)** | **✅ `leaf_gas_exchange_batch`** | flattens the allocatable-component `cfg`/PFT table per PFT; full elemental = pure-ify the 3-module solver chain + resolve cfg-elemental + C-API ripple (a separate effort) |
| `canopy_radiation` | — | **✅ already bare-array** | takes `pft_bt/lai_bt/wai_bt/tcan_bt(:)` |
| `veg_energy_step_implicit` | 1-eligible | TODO | scalar store + `inout` — `elemental` (impure allowed) candidate |
| `cas_column_step_implicit` | 1-eligible (per-patch) | TODO | scalar twins; per-patch (few elements), low GPU payoff |
| `soil_energy_step_implicit`, `column_hydrology_flux`, `soil_carbon_step` | 2 (vector) | TODO | per-layer/pool arrays; per-patch; bare-array signature helps Python |
| snow kernels | 2 | TODO | POD `snow_column_t` |
| `aero_bottom_to_top` / `canopy_aerodynamics` | — | **NOT batchable** | cohort-coupled wind cascade — not per-element independent |

**Scope statement:** the **per-cohort physiology set is complete**, now in its *final* form —
respiration is **`elemental pure`** (one procedure, no wrapper, `env`/`flux` types deleted); hydraulics
and leaf are **`_batch`** because the language/config force it (vector `psi`, heavy `cfg` gather);
`veg_energy_diagnostic` + `ground_surface_fluxes` are marked `elemental` for future array callers. The
remainder (soil/CAS/snow/biogeochem per-patch, the tangled veg-energy caller loop, the coupled
aerodynamics) is the real staged work-list — not a uniform sweep.

---

## Appendix A — ED2 fidelity note (gs/GPP/Rd frozen)

`ED2/ED/src/dynamics/rk4_driver.F90`: per DTLSM per patch — `plant_hydro_driver` (uses previous-step
`fs_open`) → **`canopy_photosynthesis`** (:259, sets `gsw_open/closed`, `A_open/closed`, `leaf_resp`,
`fs_open`) → respiration → `copy_rk4patch_init` → **`integrate_patch_rk4`**. Photosynthesis is **not**
re-called inside the substep loop (no call in `rk4_integ_utils.f90` / `rk4_derivs.f90` / `rk4_misc.f90`).
The RK stage RHS (`rk4_derivs.f90:leaftw_derivs`, 6 RKF stages) only *reads* the frozen conductance and
biochemistry; the transpiration flux = frozen `leaf_gsw` in series with `leaf_gbw` × evolving leaf→CAS
shv gradient, and NEE uses frozen `gpp − leaf_resp`. Aerodynamic conductances (`leaf_gbw`, u*) *are*
refreshed per stage via `update_diagnostic_vars → canopy_turbulence8`. **Conclusion:** MEDS's frozen
`column_prepass` (gs/GPP/Rd once per `dt_fast`) is faithful; the IMEX partition is applied only to the
integrated stores, and gs/GPP/Rd remain Category-0 constants in both scheme families.

## Appendix B — ED2 fidelity note (prognostic leaf/wood energy)

`ED2/ED/src/dynamics/rk4_derivs.f90`: leaf and wood **internal energy are prognostic RK4 states**. Their
tendencies are assembled from the absorbed shortwave (`rshort_l` / `rshort_w`) plus net longwave, sensible,
and latent terms (`dinitp%leaf_energy(ico) = initp%rshort_l(ico) + …` at :1660; `dinitp%wood_energy` at
:1942), the two are coupled by an internal leaf↔wood conductive heat flux `qwflux_wl`
(`dinitp%leaf_energy += qwflux_wl`, `dinitp%wood_energy -= qwflux_wl` at :2122–2123), and the combined
vegetation energy is `dinitp%veg_energy = leaf_energy + wood_energy` (:2143). Leaf/wood **temperature** is
diagnosed from the prognostic internal energy + water mass via the thermodynamic inverter. When a cohort
has too little leaf/wood to resolve, the tendency is set to 0 (:1737 leaf, :2013 wood) — i.e. it is
slaved, not integrated. **Conclusion:** prognostic leaf/wood energy is ED2-faithful; MEDS's
`vegetation_energy_time_derivs` (§2.2) brings this into `column_derivs` so every scheme can integrate it
and the benchmark can quantify diagnostic vs prognostic.
