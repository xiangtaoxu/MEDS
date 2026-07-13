# MEDS fast-loop integrator evaluation & production `dt_fast` recommendation

**Status:** plan (design-only). **Test bed:** Ithaca NY, 1-year ERA5-Land hourly forcing
(`data/forcing/ithaca_forcing.nc`, recycled). **Goal:** recommend a production `dt_fast` and integrator, and
settle whether MEDS should add ED2-style met/radiation interpolation.

This plan is the concrete simulation program behind the numerical-scheme comparison. It is grounded in two prior
analyses (verified against the MEDS + ED2 source): the *accuracy-testing / Strang* protocol and the *`dt_fast` /
forcing-interpolation* evaluation.

---

## 1. Background: what actually limits fast-step accuracy

The fast step is a **Lie–Trotter operator split** wrapping the coupled core, so it is **globally first-order in
`dt_fast`** even though the ARK2 core is 2nd-order. Its O(`dt`) error has four sources:

| term | source | lever |
|---|---|---|
| `e_split` | Lie–Trotter operator ordering (water → CAS/soil-energy core → hydraulics) | Strang / P3 Picard |
| **`e_frozen`** | **zeroth-order HOLD of the derived pre-pass** (two-stream RT, gₛ/GPP, aero conductances, respiration, scratch Richards BCs) — built once per `dt_fast` in `build_column_frozen` and held | **`dt_fast`** (or an inner pre-pass mini-cycle) |
| `e_theta` | soil θ committed out-of-band from the scratch Richards, read as lagged θⁿ into soil-energy props | P3 Picard |
| `e_transp` | hydraulics advanced with endpoint transpiration over the full `dt` | P3 Picard |

**`e_frozen` is the dominant, `dt_fast`-controlled term** and it lives where the forcing changes *within* a step —
sunrise/sunset SW ramps (→ GPP) and storm onset (→ ET/soil). This is the quantity the whole evaluation exists to bound.

**Forcing-interpolation finding (verified, corrects a common premise):** ED2's `dt_radinterp` is **not** a
finer-than-`DTLSM` met refresh — it defaults to `DTLSM` and is only the sub-sampling resolution of the daytime
`⟨sec z⟩` SW-downscaling quadrature (the functional analogue of MEDS's `cosz_reconstruct_factor`/`N_COSZ_SUB`, which
MEDS already has). ED2 recomputes radiation at `radfrq ≥ DTLSM` (held across steps) and marches its RK4 against met
snapshotted once per `DTLSM` at the **left endpoint**. MEDS samples met at the **sub-interval midpoint** each
`dt_fast` (linear interp + cosz-subsampled SW) — **strictly higher-order than ED2**. So there is no un-implemented
"ED2-style interpolation" to adopt; the only genuine `e_frozen` lever is **`dt_fast`**.

---

## 2. Integrator schemes under test

| id | scheme | config | role |
|---|---|---|---|
| **SPLIT** | legacy operator-split + Picard, fixed `dt_fast` | `time_integrator="split"` | production default; may need a finer rung for equal accuracy/robustness |
| **ARK-fixed** | IMEX-ARK ARS(2,2,2), fixed substeps | `time_integrator="ark"`, `ark_adaptive=false`, `ark_fixed_substep ∈ {1,2,4}` | clean ARK-core convergence at equal work; the *water/ET* reference (refines the θ commit cadence) |
| **ARK-adaptive** | IMEX-ARK + embedded-error controller | `ark_adaptive=true`, `ark_rtol ∈ {1e-2,1e-3,1e-4}` | `dt_fast` = outer forcing-refresh/pre-pass-hold cadence; **valid only for energy/carbon metrics** — θ is operator-split (committed once per `dt_fast`, excluded from `state_wrms`), so it is **not** a fine march for water/ET |
| RK4 oracle | explicit RK4 over `column_derivs` (frozen forcing) | `rk4_column_step` in tests | isolates `e_split` only (see §4); cannot run at production `dt_fast` (stable only < ~17 s) |

---

## 3. The simulation matrix

**Held fixed across all runs:** Ithaca site + PFT community + initial cohorts, `fast_context_t`/`column_config_t`
MVP params, the real hourly ERA5-Land file, RT join on, `dt_slow = 1 day` (86400 s — the production demographic step;
all ladder values integer-divide it, so `n_fast_per_slow` auto-derives).

**`dt_fast` ladder (halving):** `1800, 900, 600, 450, 225, 112.5 s` (+ `56.25 s` on high-gradient days for
truth-confirmation). `1800` = default; `600–900` = ED2 `DTLSM` range; the fine end anchors truth.

**Run types per cell `(dt_fast × integrator)`:**
- **Full Ithaca year** — daily + annual diagnostics, wall-clock, year-long conservation/robustness.
- **High-gradient days** (fast-cadence output) — (i) clear-sky near-equinox (worst sunrise/sunset `e_frozen`),
  (ii) heavy-storm day (precip onset + SW dropout + LW swing), (iii) frontal passage (fast T/wind). Host the
  sub-daily axis + the truth spot-check here.

---

## 4. Reference ("truth") construction

There is **no fully-coupled reference in the code** — truth is built by convergence, at two levels:

- **Full-step / ecological truth = LIVE-forcing Richardson.** Run the identical Ithaca config down the halving
  `dt_fast` ladder **with forcing ON** and take the finest as truth. *Held-forcing Richardson is invalid here* —
  `e_frozen` is identically zero under a constant forcing, so it would report artificially clean 1st-order behavior
  and recommend a `dt_fast` far too coarse for the real diurnal cycle. **Confirm the finest run is self-converged**
  (e.g. 112.5 vs 56.25 s agreeing well below tolerance on the high-gradient days) before adopting it as truth.
  *Water caveat:* θ comes from `column_hydrology_flux`, which already substeps adaptively, so soil moisture is
  near-converged at coarse `dt_fast`; only the θ **commit cadence** refines — the water axis may legitimately admit
  a coarser knee (don't misread as under-resolution).
- **`e_split` isolation (column-kernel).** At a fixed `dt` and *fixed frozen coefficients*, compare the ARK2
  split-march against the RK4 oracle (`rk4_column_step` over `column_derivs`): the difference is `e_split` alone.
  **Two traps:** use `ark2_column_step` (2nd-order, matching a converged oracle), *not* `column_be_stage +
  advance_hydraulics_full` (that's 1st-order IMEX-Euler and folds in the BE core-discretization error); and hold
  θ = θⁿ on both sides (don't commit `fro%theta1`). *There is no "sub-cycled-pre-pass at fixed 1800 s march" config*
  — one `dt_fast` controls both the pre-pass hold and the march step, so refreshing the pre-pass K× just is the fast
  loop at `1800/K`.

**Cross-scheme (SPLIT vs ARK at equal `dt_fast`) is a *consistency* check only** — both share the frozen pre-pass and
commit θ out-of-band, so their agreement cancels shared bias and bounds only the Picard-vs-ESDIRK core difference,
**not** distance to truth. The standing <1% 30-yr agreement is a heuristic prior, never the arbiter.

---

## 5. Metrics & thresholds (three timescales)

- **Sub-daily / diurnal** (where `e_frozen` lives; high-gradient days, fast cadence): CAS air-temp diurnal RMSE
  **< 0.25 K**; canopy GPP diurnal RMSE **< ~5% of midday peak** + sunrise **GPP-onset timing error < one `dt_fast`**;
  ET/latent-heat RMSE **< ~5 W/m²**; top-soil skin-temp RMSE **< 0.5 K**. *CAS-temp and dawn GPP-onset are the sharpest
  discriminators and blow first.*
- **Daily-integrated** (the fast→slow handoff — **DECISIVE tier**): daily GPP & NPP **bias < 1%**, daily RMSE < ~2%;
  daily ET bias < 1%. (From the `gpp_accum`/`*_resp_accum` accumulators → `total_gpp`/`total_npp`, which reflect `dt_fast`.)
- **Annual / ecological**: annual GPP/NPP, LAI (seasonal + max), AGB each **< 1%**. Most forgiving — necessary but
  **not** sufficient (day-to-day errors + demographic buffering cancel).
- **Conservation (hard gate):** worst `whole_energy`/`whole_water` residual over the year stays at the finest run's
  closure order (SPLIT: machine energy; ARK: machine energy, split-order water). Any `dt_fast` that degrades closure
  is disqualified.
- **Stability/robustness (hard gate):** zero budget failures, Picard/Newton convergence every step, bounded temps,
  no NaN over 8760 h; for ARK-adaptive, log `nrej/nsteps`.
- **Metric:** the shipped per-reservoir `state_wrms` (scales tied to `ark_rtol`) for column-kernel work — **add θ back
  to the norm for a Richardson study; do NOT re-add ψ (already included)** — plus human-readable per-store deltas.

---

## 6. Expected comparisons

- **Column-kernel order** (fixed forcing): ARK2 shows p≈2 only on cleanly-integrated tableau variables (the decoupled
  affine CO₂ twin — `test_ark2` asserts p≥1.9); p≈1.2 on split-coupled variables; IMEX-Euler p≈1.
- **Full-step order** (live forcing, `dt_fast` ladder): **both SPLIT and ARK plateau at p≈1**, capped by the split +
  `e_frozen` hold, *not* the ARK2 core. Seeing ARK ≈2 on the CO₂ twin but ≈1 on the full step is the expected,
  self-consistent picture.
- **Daily/annual tiers:** almost certainly clear <1% at `dt_fast`=1800 s for both integrators (consistent with the
  30-yr <1% ARK-vs-SPLIT match). The **binding constraint, if any, is the sub-daily tier** on transition/storm days.
- **SPLIT vs ARK:** ARK is L-stable → `dt_fast` chosen purely by accuracy; SPLIT may need one rung finer for the same
  accuracy/robustness (Picard convergence at large `dt_fast`). ARK is ~15–20% slower per step at equal `dt_fast`
  (implicit FD-Jacobian Newton march ≈ 33% of the step; shared pre-pass ≈ 67%).

---

## 7. Decision criteria

1. **Conservation + stability gates** must pass (hard prerequisites).
2. **Production `dt_fast` = the LARGEST (cheapest) `dt_fast`** whose **decisive-metric** error (daily GPP/NPP bias
   < 1%) is within tolerance of the live-forcing truth, **subject to** the sub-daily diurnal thresholds also holding
   on the high-gradient days.
3. Report separately for SPLIT vs ARK; report the error-vs-wall-clock "knee" for each metric.

**Strang decision:** adopt a Strang split (necessarily *with* midpoint forcing) only if a probe shows `e_split/e_tot
> 0.10` **and** a midpoint+Strang run moves a headline diagnostic > 1% **and** the ~1.6–2× cost is acceptable —
expected to fail #1. Cheaper alternatives if the split error matters: **midpoint forcing** (fixes the dominant
`e_frozen` at ~1.67× — better spend than Strang), the already-designed **P3 outer Picard** (kills `e_split`+`e_theta`+
`e_transp` on the single pre-pass, cheapest), or **shorter `dt_fast`** (halves all four terms).

---

## 8. Forcing / ED2-interpolation decision

**Keep MEDS's current per-sub-step midpoint linear interpolation + cosz-subsampled SW disaggregation** — already
best-in-class and strictly higher-order than ED2's left-endpoint hold; **do not build an ED2 `dt_radinterp` knob**
(MEDS already has its analogue). The only `e_frozen` lever is `dt_fast`. A pure RT-only refinement (the literal
`dt_radinterp` idea) **won't cut the dominant GPP error** — GPP comes only from the leaf-Cᵢ solve consuming the coarse
absorbed PAR, so refining RT alone helps *energy* (leaf/CAS temp) but not GPP at dawn/dusk. The one defensible
cheap-decouple, *only if the test shows `e_frozen` dominates*: an inner **RT+aero+leaf-Cᵢ** mini-cycle leaving only
Richards frozen (soil moisture barely moves in <1800 s).

---

## 9. Prerequisite (Step 0): sub-daily instrumentation — **implement first**

The evaluation's sub-daily axis has **no output path today**: `output_integrate` ticks once per slow step; the FAST
tier is registered but deferred; and **latent-heat/ET and CAS-temperature are not registered variables at any tier**.
Before any high-gradient-day run:
- Register **CAS temperature**, **latent heat / ET**, and **top-soil temperature** as site-level output variables with
  registry extractors (daily/monthly/annual), and
- Implement a **fast-cadence (sub-daily) output** path for these + GPP on the high-gradient days (either the FAST
  aggregation tier or a bespoke per-sub-step probe), and
- Thread `run_fast_biophysics`'s `worst_energy`/`worst_water`/`n_budget_fail` out-args (computed but discarded on the
  production path) up for year-long conservation/convergence logging.

Daily GPP/NPP and annual GPP/NPP/LAI/AGB already have valid `dt_fast`-reflecting output paths.

---

## 10. The single most informative first run

One clear-sky near-equinox Ithaca day (max dSW/dt = worst-case `e_frozen`), `INTEG_ARK`, per-sub-step (60 s) output of
GPP / ET / CAS-temp, the **isolation quartet**:
- **A** = `dt_fast` ≈ 90 s validated fine reference (45 s spot-check),
- **B** = `dt_fast` = 1800 s (current default),
- **D** = cheap-refined: 1800 s coupled march, but met+two-stream RT+aero refreshed at `dt_rad` = 300 s, leaf-Cᵢ and
  Richards held at 1800 s.

`(B−A)` measures worst-case `e_frozen` at the default and whether dawn/dusk cancellation keeps 1800 s inside budget —
the whole `dt_fast` decision. `(D−A)` vs `(B−A)` settles the interpolation question: if **D collapses toward B**, the
sunrise GPP error is leaf-Cᵢ-gated → an ED2-style RT-only layer is provably useless (only `dt_fast` helps); if **D
recovers toward A**, a targeted inner mini-cycle is worth building.

---

## 11. Provisional production recommendation (pending the runs)

- **Integrator: `INTEG_ARK`** — L-stable, no stability floor, so `dt_fast` is accuracy-driven, not Picard-limited.
- **`dt_fast`: 1800 s default**, fall back to **900 s** only if transition/storm-day daily GPP or ET bias exceeds
  ~1–2% vs the fine reference. This straddles ED2's 600–900 s but for the opposite reason: MEDS already interpolates
  forcing to any `dt_fast` at the midpoint, so it should run **coarser** than ED2 for equivalent forcing fidelity.
- **Forcing handling: unchanged** (midpoint + cosz). Add an inner pre-pass mini-cycle only as a contingency if the
  isolation test proves `e_frozen` dominates *and* config D recovers it.
