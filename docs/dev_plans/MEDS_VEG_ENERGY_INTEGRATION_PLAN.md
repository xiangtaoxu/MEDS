# Vegetation energy: the exact-exponential tissue store

**Status:** 2026-07-31. **P0 and P1a LANDED**; P1b/P1c (activation) parked on
`wip/veg-store-activation`, NOT green. Supersedes `MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §3 and §5-P4,
and the "bordered arrowhead" scoping in `MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md` Phase 5.

**Question this answers:** how should MEDS integrate leaf and wood temperature, for thousands of
cohorts, with real heat capacity, and what can serve as an accuracy reference?

---

## 1. The result: the tissue ODE has a closed form, so do not discretise it

Under the Category-0 coefficient freeze the tissue equation is **linear with a known timescale**:

```math
\mathrm{cap}\,\frac{dT}{dt} = \mathrm{numer} - \mathrm{denom}\,(T - T_{cas}),
\qquad \tau = \frac{\mathrm{cap}}{\mathrm{denom}}
```

So the step is not something to approximate. `veg_energy_diagnostic` now uses the exact solution,
with **two different weights**:

```math
x = \frac{\Delta t}{\tau} = \frac{\mathrm{denom}}{a_{\rm store}},\qquad
w_{\rm end} = e^{-x},\qquad
w_{\rm avg} = \frac{1-e^{-x}}{x}
```

`w_end` weights the old state at the **endpoint** — the committed state. `w_avg` weights it in the
**step average** — every reported flux. Using one for both is what breaks conservation. Pairing them
makes the balance close identically:

```math
a_{\rm store}(\Delta T_{\rm end} - \Delta T_{\rm prev}) + \mathrm{denom}\cdot\Delta T_{\rm avg} = \mathrm{numer}
```

which reduces to `a_store(1 − w_end) = denom·w_avg`, true by construction.

Since `a_store = cap/dt`, the ratio `denom/a_store` is `dt/τ` — **`dt` cancels**, so the kernel
signature is unchanged.

**Why this beats backward Euler, which is what the store used before.** BE is the
`w_end = w_avg = 1/(1+x)` approximation. At production `dt_fast` a leaf has `x ≈ 144`: the exact
endpoint weight is `~5e-63` against BE's `0.0069`, and an SDIRK2 tableau gives `−0.031` — a
**sign-alternating** artificial memory for a mode whose true memory is nil. The prognostic tissue's
only measurable effect under BE was the discretisation's own artefact.

**Two limits, both exact, both reached continuously.** `a_store → 0` gives `x → ∞`, both weights
→ 0, and the result is the pure diagnostic balance. So **"diagnostic" is the kernel's zero-inertia
limit, not a separate code path** — which is what lets the `leaf_energy_model` / `wood_energy_model`
selectors be deleted. `a_store → ∞` gives `x → 0`, both weights → 1, and the tissue holds its
temperature.

`denom` now **excludes** `a_store`: it is the coupling conductance, a property of the surroundings,
not of the tissue's own inertia. `veg_coupling_floor` and its `g_slave` routing apply to the coupling
alone and are unchanged.

Fortran has no `expm1`, so the small-`x` branch uses a series (relative error `~x⁴/120 < 1e-18` at
the `1e-4` threshold). Without it `1 − exp(−x)` cancels catastrophically as `x → 0` — exactly the
large-capacity limit a big cohort lives in.

---

## 2. Why NOT the bordered arrowhead

The earlier plan was a per-cohort Schur elimination into the ARK's 2×2 CAS Newton. The algebra is
sound and stays O(n), but it is the wrong tool, and it had three blockers that were not identified
when it was scoped:

1. **`t_store0` cannot live on `surface_frozen_t`.** `ark2_column_step` builds stage 3 as
   `base3 = y + (1−γ)dt·f_I(Y₂)` and then does a BE step of size `γ·dt` *from `base3`*. A leaf DOF's
   storage reference in stage 3 must be `base3%leaf_temp`, not `T^n`. Seeding it once per `dt_fast`
   — which `wip/unified-veg-energy` does — makes both stages relax from `T^n`, which is not a
   discretisation of the leaf ODE at all.
2. **The `β = 2.414` stage-3 extrapolation destroys the embedded estimator for a stiff DOF.** For a
   mode already collapsed to equilibrium in stage 2, `state_err_diff` returns `≈ 3(T⁰ − T_eq)`;
   against a `1e-2 K` atol a 1 K diurnal offset normalises to ~330 and the controller rejects to
   `dt_floor` permanently.
3. **`state_wrms_grouped` is a MEAN over `3 + 2·nsl + 2n` components.** Adding per-cohort tissue
   terms makes step-size control a function of the **cohort partition** — splitting a stand into 2×
   the cohorts at half `nplant`, physically identical, halves the CAS terms' weight. At n = 1000 the
   three CAS scalars are 3/4003 of the error norm. (This already bites `GRP_LEAF_W`/`GRP_WOOD_W`.)

The exponential form needs none of it: no new tableau state, no WRMS group, no Newton, no arrowhead.

---

## 3. Scalability: what actually blocked thousands of cohorts

Neither blocker was the integrator.

**The O(n²) sort.** `aero_bottom_to_top` built its bottom-to-top ordering with a selection sort — the
only superlinear term in the fast loop, run once per `dt_fast` on **all three schemes**. Measured
over a 6-day run: 0.97 s at n = 2000, 3.9 s at n = 4000, **16.5 s at n = 8000**. It is also
redundant: `sort_cohorts` leaves the block height-descending and the gather preserves that, so
bottom-to-top is the reverse. Detected in O(n), with the sort kept as a fallback for unit tests that
build cohorts unsorted. The two branches agree **bit-identically including ties** — the sort's `<=`
keeps the last index achieving the running minimum, which in a descending array is always the largest
remaining index.

**The stack ceiling.** ~740 B of automatic arrays per cohort, live across four nested frames
(`column_fast_step` → `build_column_frozen` → `column_prepass` → `aero_bottom_to_top`): a hard
SIGSEGV near 10–12k cohorts, and it scales *down* with `OMP_STACKSIZE` once the patch loop is
threaded.

**The fix is a runtime one** — `ulimit -s unlimited` (or `ulimit -s 65536` for ~90k cohorts), plus
`OMP_STACKSIZE` once threaded. A `-DMEDS_HEAP_ARRAYS=ON` option moves large automatic arrays to the
heap instead, but it is **OFF by default**: ifx's `-heap-arrays` segfaults `test_output_registry` at
the exit of `apply_variable_override` in a Debug build. That is a compiler codegen interaction, not a
MEDS bug — the routine assigns only plain integer components and is clean without the flag. Prefer
raising the stack; if you must use the flag, re-run the Debug suite.

**Still open:** the dispatcher pays split's ~310 B/cohort of scratch on every ARK/RK45 run
(`column_fast_step` declares its full automatic set before dispatching); `surface_derivs` and
`column_derivs` allocate 5 size-n arrays *per call*; and `fs = fro%surf` deep-copies **14 size-n
arrays per stage** to overwrite one scalar. These are refactors, not flags.

---

## 4. The accuracy reference

**`rk45_tight` cannot anchor tissue temperature, and the evidence that said it could is circular.**
With `a_store = 0` and every coefficient frozen, `T_leaf − T_cas = numer/denom` is *constant within a
`dt_fast` by construction*, so leaf-T error ≡ CAS-T error plus a fixed offset, identically. The sweep
could not have produced any other answer. It is algebra, not a measurement.

With the exponential form the tissue is **exact given the freeze**, so it cannot be the reference's
weak point at all, and the anchor question reduces to the CAS and soil columns — where `rk45_tight`
remains reasonable, subject to its own rescue-path caveat (`work_rk45_rescue_site` must be read
first).

---

## 5. Physics decisions taken (2026-07-31)

1. **Heat capacity = dry biomass + INTERNAL (xylem/symplast) water.** The canopy **surface film is a
   separate store carrying its own phase**, so intercepted snow can stay frozen and a freezing film
   can pin the tissue at the melt plateau. Folding film mass into a *temperature*-based capacity
   would commit to a film that can never freeze — which is the larger error, since fusion is
   ~2×10⁵ J/m² ≈ 80 K of the film's sensible capacity, a plateau **longer than `dt_fast`**.
2. **Wood is lumped (internally isothermal)**, accepting the Biot-number error, per the user's
   decision. Bi ≈ 7.45 for a bole skin against the `Bi ≪ 0.1` that lumped capacitance needs, so real
   wood is a diffusion problem with no single τ. This is a *sizing and structure* limitation, not an
   integrator one: the store is right, its capacity is too small, and replacing the placeholders
   later moves τ without touching the method.
3. **The tissue emits LW at its own temperature**, always. It previously emitted at `t_cas` except
   under Picard, which was defensible only while the tissue had no capacity of its own.

---

## 6. Known-wrong things found and NOT fixed

**`τ_wood = 55–199 s` is circular.** `wai = 0.20·lai` and `bsap = 0.10·wood_carbon` are hardcoded
placeholders (`meds_fast_dynamics.f90:344-346`) and `test_wood_stiffness_spread` hardwires the same
ratios, so the test re-measures the placeholder. The implied wood is a **0.6–2.3 mm shell**.
Physically sized, τ ≈ **1300 s** (branch wood) to **3690 s** (bole) — so wood memory at
`dt_fast = 1800 s` is `w ≈ 0.42`, not 0.10, and the "prognostic ≈ diagnostic for wood" conclusion is
an artefact of the placeholder. `wood_temp` drives `stem_maintenance_respiration`, so the resulting
22–60 min phase lag is not cosmetic.

**Canopy phase change is absent entirely.** Snowfall is intercepted as liquid with no fusion debit
and later leaves via `enthalpy_vapor` (evaporation, never sublimation). This is what decision (1)
above exists to enable.

**`veg_energy_step_implicit` carries the LW emission response in its Jacobian (`−8εσT³·AI`) but not
in its residual**, so its fixed point lacks the radiative restoring term `veg_energy_diagnostic` has;
the two also differ on the slope (`8εσT³` vs `4εσT³`) and *neither* matches a cylinder's `π·WAI`. It
also renormalises its reported fluxes by an **unbounded** `scale = turb_avg/turb_star` guarded only
by `tiny_num`. Unifying on `veg_energy_diagnostic` retires this kernel and the whole class with it.

**Sunlit/shaded separation does not exist** (zero hits for `sunlit|shaded|fsun`), so canopy GPP
carries the standard big-leaf concavity bias and one leaf temperature cannot span the sunlit/shaded
spread. Out of scope here, but it dominates the physics error budget this work sits inside.

**The free-convection derivative is missing from the linearization.** `boundary_gbh_mos` includes
Grashof free convection, so `H ∝ ΔT^1.25` and the true slope is `1.25·h`, but `h_coeff` is frozen and
enters `denom` as `1.00·h` — the solved `ΔT_leaf` is overstated ~20% in calm conditions.

---

## 7. What remains

`wip/veg-store-activation` turns `a_store` on for the split leaf. **It is not green** — 5 tests fail.
Two are the expected mid-conversion comparison (`SNOW ARK/RK45 agrees with split`), which cannot pass
until ARK/RK45 convert in the same change; the other three are **not diagnosed**. Candidates recorded
on that branch: the `LAI_SLAVE_MIN` branch sets `leaf_temp = tcas` and `cycle`s, skipping the store
accumulation (an unaccounted `cap·(tcas − t_emit)` once cap > 0); `bio%leaf_temp` initialises to
`LEAF_TEMP_INIT = 288.15 K` while the CAS may be ~300 K, so the first sub-step carries a large
legitimate transient; and `transp`/`coh_qsoil` now use the step-average rather than the endpoint,
changing the demand handed to plant hydraulics.

Then: the separate film store with phase, deletion of both `*_energy_model` selectors, and the
honest wood sizing.


---

## 8. The tissue store and the canopy air must be solved TOGETHER (2026-07-31)

**Status: ARK and RK45 carry the store correctly and close every budget to ~1e-6 J. The `split`
path does not, and no fixed-point scheme tried so far fixes it.**

### What the store does to the canopy-air coupling

The step-average weight `w_avg = (1-e^-x)/x` puts part of the tissue's flux on its own
start-of-step temperature. When `tau << dt_fast` that start-of-step temperature is very nearly the
**previous canopy-air temperature**, so the flux feeds a lagged `t_cas` back into the current air
balance. The loop gain is `denom*w_avg*dt/(air heat capacity)`:

| tissue | tau | `w_avg` | gain |
|---|---|---|---|
| leaf | ~20 s | 0.011 | ~0.03, harmless |
| wood | ~310 s | 0.172 | **~0.51 alone**, unstable once leaf + ground are added |

Note where the danger is: `w_end = 0.003` for wood, so its **endpoint temperature is essentially
diagnostic** — the lag lives in the FLUX, not the temperature. "The wood barely lags" does not imply
"the wood barely matters", and that is why the earlier conclusion that prognostic wood ~ diagnostic
wood was wrong: it looked only at the temperature.

The leaf looked exempt only because its heat capacity is small. Correctly-sized wood carries ~36 kg
dry/m2 against a leaf's ~1.

### Three things were tried on `split`. All measured, all insufficient.

1. **Single explicit pass** (what `SCHEME_SPLIT_SEQUENTIAL` does). Diverges: the canopy air
   alternated 291 K <-> 311 K every step with a growing envelope (`dh` = -175, +167, -195, +186,
   -214, +224 W/m2), midday CAS 304.7 -> 314.9 K, and the soil DRIED under a prescribed water input.
   This is the `theta`/NEE/melt-agreement test damage — not a re-baseline.
2. **Force the Picard iterate.** Not sufficient either: on a full-site fixture the ITERATION itself
   oscillated (+-0.4 K leaf/air alternation, growing) to floating overflow. The pre-existing note in
   `column_fast_step` had required Picard for a prognostic *leaf* and explicitly exempted wood --
   "wood has no transpiration feedback and is stable on the pure-split path". That exemption is the
   thing that broke.
3. **Schur-precondition the canopy-air solve** with `cas_cond = (h_coeff+g_slave)*w_avg`, blended
   into the existing atm conductance slot. This DOES restore contraction (the iterate converged
   cleanly where it previously overflowed) and is provably transparent at the fixed point — the added
   source and sink cancel identically, so it changes the path, not the answer. But the conservation
   budgets still fail, and **not for want of iterations** (200 made no difference): the budget terms
   are assembled from the under-relaxed iterate `tcas` while `enth1` comes from the solve, and the
   `picard_relax = 0.5` seed keeps a persistent mismatch between them.

### What this means

The tissue store with realistic wood mass is a genuinely COUPLED stiff system: `cap_wood/cap_cas ~
0.5` with `tau_wood << dt_fast`. It needs a coupled implicit solve — a Newton on (tissue, CAS)
together — not a fixed-point iteration with a fixed relaxation. **ARK already has that**
(`newton_surface_solve`), which is exactly why ARK and RK45 close every budget with the store on
while `split` does not.

Recommended paths, in order:
- Accept that `split` cannot carry the store, and gate it (ARK is the production integrator anyway).
- Or extend `newton_surface_solve`'s coupling to the split path.
- The Schur preconditioner is worth keeping in either case — it is the correct implicit form of the
  tissue's contribution and it is a no-op at the fixed point.

Do NOT re-baseline the `split` goldens: they are reporting a real instability, not new physics.


---

## 9. The non-monotonicity resolved: a PRE-EXISTING dt_fast-cadence oscillation (2026-07-31)

**The tissue store does not cause the four failing physics assertions. It shifts the PHASE of an
oscillation that is already there, and the assertions sample a single step of it.** That is the whole
explanation for why leaf-store-alone (10 failures) looked worse than both stores together (4).

### What the measurement shows

`test_column_dynamics` with the store OFF — i.e. code that passes 36/36 — canopy-air temperature over
consecutive 900 s steps through mid-morning:

```
286.5 287.1 288.2 289.6 291.3 293.2 295.2 297.0   <- smooth warm-up
294.0 297.1 294.4 298.1 290.9 296.1 300.6 291.0   <- period-2 onset
297.2 299.5 292.0 298.9 292.8 300.2 292.6 300.4   <- ~8 K, sustained
```

The forcing is a smooth `cosz` ramp (`t_air = 288 + 6*(cosz - 0.3)`, range 286.2-291.6 K), so nothing
in the boundary condition alternates. This is numerical. It also explains the second oddity seen at
the same time: sapflow runs 3-4x the transpiration demand (3.3e-4 vs 7.5e-5 kg/plant/s) and leaf
water sloshes 0.90 <-> 1.05 kg/plant in antiphase — the plant is chasing an oscillating leaf.

**This fixture had NEVER been run on ARK before.** It set `time_integrator = INTEG_SPLIT` explicitly;
retiring split is what pointed it at ARK and exposed this.

### It is the OUTER FREEZE, not the inner march

`ark_adaptive` is on by default, so the embedded-error controller is already sub-stepping inside each
`dt_fast` — and the oscillation survives it, because the Category-0 coefficients are re-frozen once
per `dt_fast` and the oscillation lives at exactly that cadence. Quartering the freeze interval
(`dt_fast` 900 -> 225 s) cuts the amplitude from ~8 K to ~4 K:

```
dt_fast = 900 s   292.6  300.4  292.9  300.9        swing ~8 K
dt_fast = 225 s   295.9  297.5  296.7  298.2  294.1 swing ~4 K
```

So no amount of inner adaptivity fixes it; only a shorter freeze does. The stiff leaf<->CAS coupling
is being evaluated against coefficients held fixed across the whole step.

### Why this matters more than the store

1. **The budgets do not catch it.** Every conservation ledger closes to ~1e-6 J while the canopy air
   swings 8 K step to step. Conservation is not stability, and this is the sharpest example yet.
2. **It is the same family as the split failure in §8** — a lagged canopy-air temperature inside the
   canopy-air balance — which is why the store makes it worse rather than causing it.
3. **The store cannot be validated until this is fixed.** Its four failing assertions are
   phase-samples, so any judgement about the store's physics made on top of this is meaningless.

### Next

Fix the freeze, not the store. Candidates, cheapest first: re-freeze the leaf<->CAS coefficients per
inner sub-step rather than per `dt_fast`; or bring the coefficients that carry the feedback (the
latent slope `le_slope`, which depends on the live canopy-air state) inside `newton_surface_solve`.
Only then re-run the `TISSUE_STORE_SCALE = 1` assertions.
