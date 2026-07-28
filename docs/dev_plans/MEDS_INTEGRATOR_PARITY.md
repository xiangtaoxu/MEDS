# MEDS integrator parity — split vs ARK vs RK45

**Status:** Phase A (instrumentation) implemented. Phases B–D open.

## Why this document exists

MEDS has three fast-loop time integrators, selected by `[fast].time_integrator`:

| selector | module | character |
|---|---|---|
| `split` (default) | `src/driver/meds_fast_split.f90` | operator-split + Picard, implicit CAS, BE-Thomas soil |
| `ark` | `src/driver/meds_fast_ark.f90` | IMEX-ARK ARS(2,2,2); implicit CAS+soil heat, mass operator-split out |
| `rk45` | `src/driver/meds_fast_rk45.f90` | ED2-faithful adaptive Cash-Karp 5(4), fully explicit over `column_derivs` |

They are routinely compared against each other. **They do not implement the same model**, and
several of the differences are on by default, so an uncontrolled comparison mixes physics with
numerics. That confound has already produced at least one wrong conclusion in this project's
history (an early split-vs-ARK ranking, invalidated once the ARK-only condensation sink was found).

This file is the single home for the difference inventory: what differs, whether it is *physics*,
*bookkeeping* or *numerics*, and what the plan is for each. Update it when a row changes state.

---

## 1. The inventory

Verified by code reading against `main @ aee5d68` (2026-07-28).

### Class 1 — model-family differences (the schemes integrate different physics)

| # | difference | split | ARK | RK45 | class | plan |
|---|---|---|---|---|---|---|
| 1 | CAS supersaturation sink (`TAU_COND`, `meds_fast_time_derivs.f90:175`) | **absent** | present | present | physics | **C3** — give split the same sink; route condensate to a store instead of out of the column |
| 2 | Snow / temporary surface water | present | **imports kernels, never calls them** | **same** | physics | **C4** — hoist to a shared pre-column stage |
| 3 | Per-layer root sink placement (`root_sink_share`, under `multilayer_roots`) | present | `root_frac` only | `root_frac` only | physics | C6 — deferred; `multilayer_roots` defaults off |
| 4 | Prognostic leaf/wood energy | present | hard `error stop` | hard `error stop` | physics | C6 — deferred; the error stop makes it non-silent |
| 5 | Non-free-drain bottom BC / Zeng–Decker | present | hard `error stop` (`meds_fast_ark.f90:838`) | hard `error stop` | physics | **C5 DONE** |
| 6 | transp↔uptake seam | re-solved inside the Picard iterate (realized, supply-limited) | frozen from state^n | frozen from state^n | numerics | measured in Phase B — **not** the cause of the structural gap (`no_water` mask leaves it unchanged) |
| 12 | **ground surface energy seam under a canopy** — constant +1.10 K soil-surface offset, day and night, surviving 12× refinement (§3b B-2) | via `ground_surface_fluxes` | via `surface_derivs` | via `surface_derivs` | physics | **NEW — Phase C, highest priority**: the only difference found that dt refinement cannot remove |

Row 1 is the one that most affects existing results: `cas_condensation` defaults `.true.`, so the
30-yr split run and the 29-yr RK45 run in `runs/ithaca_ark30/` were **not the same model**. Worse,
both ARK and RK45 book `sf%cond` into `whole_wat_out` (`meds_fast_ark.f90:224`,
`meds_fast_rk45.f90:217`) — condensed dew/fog leaves the tracked column entirely rather than
reaching soil or canopy. Split neither condenses nor loses it. The two paths are wrong in
different directions.

### Class 2 — unbookkept corrections (state edited outside the conservation ledger)

| # | difference | plan |
|---|---|---|
| 7 | RK45 clamped the **committed** state, not just throwaway stage inputs. `clamp_theta` moved θ with no mass debit; `clamp_soil_energy` then re-derived T at the new water mass. | **C1 DONE** — commit clamp removed; stage clamps kept |
| 8 | RK45 commits its **own** integrated θ but takes pond / drainage / runoff from the **frozen scratch** solve (`meds_fast_rk45.f90:444-448`, and `fro%drainage`/`fro%runoff_surf` in the ledger). | **C2** — derive them from RK45's own trajectory, b-weighted like `stage_bnd` already does for the CAS boundary terms |

On #7, the in-code argument is that a clamp which bites shows up as a large 5th-vs-4th discrepancy
and the controller rejects the step. That argument **does not cover the accept condition**, which
is `err <= 1.0 .or. dt <= dt_floor`: at the sub-step floor a clamped state is committed
unconditionally and the controller has no veto left.

### Class 3 — instrumentation gaps (closed in Phase A)

| # | gap | status |
|---|---|---|
| 9 | `budg%rk45_rescue` incremented but never aggregated or written — no way to tell whether an "RK45 run" was actually RK45 | **closed** — `work_rk45_rescue_site` |
| 10 | Clamp activations invisible, though they are the trajectory- and compiler-dependent term | **closed** — 4 new counters |
| 11 | `scripts/numerics_sweep.py` had no `rk45` scheme at all | **closed** — `rk45`, `rk45_tight`, `--parity` |

---

## 2. RK45 is a hybrid, not an explicit scheme

Worth stating separately because it is easy to miss and it changes what an "RK45 run" means. When
the explicit march rails or burns its work budget (`RK45_WORK_CAP`), `column_fast_step` rolls the
state back and **redoes that `dt_fast` on the split path** (`meds_fast_split.f90:235-278`). A run
whose cold months all rescue is substantially a split run wearing an RK45 label.

`work_rk45_rescue_site` is the provenance check. **Read it before reading any RK45 result.**

---

## 3. Phase A — what was added

All in `GRP_NUMERICS`, so they ride the existing opt-in `[output].numerics = true`.

| variable | units | meaning |
|---|---|---|
| `work_rk45_rescue_site` | — | `dt_fast` steps RK45 abandoned and redid on split |
| `work_clamp_stage_site` | — | stage-input clamp activations (throwaway states) |
| `work_clamp_commit_site` | — | committed-state clamp activations |
| `work_clamp_mass_site` | kg/m² | water moved by committed-state clamps |
| `work_clamp_energy_site` | J/m² | energy moved by committed-state clamps |

**Read the magnitudes, not the commit count.** This was the first thing the instrument found, and
it was not what was expected: the committed-state clamp fires on essentially *every* sub-step even
in a well-behaved column, because the Richards solve routinely lands a whisker outside
`[theta_res, theta_sat]`. Measured on `test_column_rk45`, 96 sub-steps:

| | activations | unbookkept mass | unbookkept energy |
|---|---|---|---|
| wet, unsaturated (ifx) | 1246 | **0.0** kg/m² | 3.1e-5 J/m² |
| wet, unsaturated (nvfortran) | 1209 | **0.0** kg/m² | 2.9e-5 J/m² |
| saturated (ifx) | 2842 | 142.7 kg/m² | 7.0e-5 J/m² |
| saturated (nvfortran) | 2907 | 142.6 kg/m² | **1.513e5 J/m²** |

So the count separates healthy from broken by a factor of ~2; the magnitudes separate them by up
to nine orders. Two further readings fall straight out of the table:

- The **mass** violation is compiler-independent (142.6 vs 142.7) — `clamp_theta` is a structural
  defect, not a floating-point accident.
- The **energy** violation is entirely compiler-split (7.0e-5 vs 1.5e5) — this is the previously
  known "ifx never triggers it where nvfortran does", now attributed to `clamp_soil_energy` by
  direct measurement rather than inferred by disabling the call. It also reproduces the 1.5e5 J/m²
  figure previously obtained that way, from the inside.
- In the benign case the unbookkept mass is *exactly* zero on both compilers: every activation
  there is a round-off-scale soil-energy touch, and θ never leaves its domain at all.

Both magnitudes are **gross**: a running sum of `|correction|`, not a net residual. Opposite-signed
corrections deliberately do not cancel, so these are an upper bound on the resulting budget gap.

### Harness

`scripts/numerics_sweep.py` gains `rk45` and `rk45_tight` schemes, the `HEALTH_VARS` columns, and
a `--parity` flag that pins every Class-1 row above to the common subset:

```bash
scripts/numerics_sweep.py --base cfg.toml --out runs/parity1 --parity \
    --schemes split ark rk45 --dt 1800 900 --ref-scheme rk45_tight
```

`rk45_tight` exists for a specific methodological reason. Refining `dt_fast` refines the
**Category-0 coefficient freeze** as well as the time step, so a `dt_fast` sweep measures both
error sources at once — and the freeze is known to dominate at production `dt_fast`. Holding
`dt_fast` fixed and driving the integrator tolerance down isolates the stepper. The two references
answer two different questions and should not be substituted for one another.

---

## 3b. Phase B results — fidelity with demography frozen

Setup: `scripts/parity_scenarios.py` emits four base configs (2×2 season × stand structure), all
with `[run].slow_on = false`, `dt_fast = 1800 s`, hourly FAST output, and `--parity`. Scored by
`scripts/parity_fidelity.py` on the **hourly** stream — the sub-daily temperatures live only in the
FAST tier, and a one-month daily series averages away the diel structure the schemes differ on.
744 paired hourly records per cell. All 16 ifx cells and 16 nvfortran cells ran clean.

### B-1. The schemes fall into two families, and the split is the semi-discretisation

Pairwise RMSE at production `dt_fast`, ifx (same freeze cadence for all three, so this is
integrator-and-coupling disagreement with the dominant error source held common):

| scenario | var | split~ark | split~rk45 | **ark~rk45** |
|---|---|---|---|---|
| bare winter | CAS T | 0.101 | 0.102 | **0.014** |
| bare summer | CAS T | 0.196 | 0.210 | **0.022** |
| stand winter | CAS T | 0.451 | 0.359 | **0.235** |
| stand summer | CAS T | 0.746 | 0.684 | **0.415** |

ARK and RK45 agree with each other 7–9× more closely than either agrees with split, despite being
an implicit ESDIRK and a fully explicit Cash-Karp — about as different as two tableaux get. They
share `build_column_frozen` and `column_derivs`; split builds its own frozen set and re-solves
hydraulics inside the Picard iterate. **The family boundary is the semi-discretisation, not the
tableau** — a direct extension of the earlier finding that the Category-0 freeze dominates the
stepper, now with the freeze *cadence* held common and the schemes still splitting by *what* they
freeze.

Magnitudes rise with both coupling strength and season, as predicted: stand > bare, summer >
winter. Worst cell (stand summer) is CAS T 0.75 K, soil-top T 1.2 K, LE 23.6 W/m², GPP 0.59
µmol/m²/s.

### B-2. A structural soil-surface gap that refinement does NOT remove

Refining both families 12× (dt 1800 → 150) and comparing the two refined references directly —
`split@150` vs `rk45_tight@150` — separates time-stepping error (collapses) from structural
difference (survives). Ratio = refined family gap ÷ production-dt `split~ark` gap:

| scenario | CAS T | **soil-top T** | LE | H | GPP |
|---|---|---|---|---|---|
| bare winter | 0.18 | 0.54 | 0.20 | 0.17 | — |
| bare summer | 0.11 | **0.06** | 0.12 | 0.11 | — |
| stand winter | 0.35 | **1.01** | 0.19 | 0.12 | 0.38 |
| stand summer | 0.42 | **1.00** | 0.10 | 0.57 | 0.30 |

Everything converges except **soil-surface temperature in the vegetated cells**, where the two
families disagree at dt=150 by *exactly as much* as the schemes do at dt=1800. That is not an
accuracy difference — it is two different problems being solved, and no dt refinement will
reconcile them.

Characterising it (stand summer, refined, split-family minus ARK-family):

```
 bare  summer   mean +0.0006 K   sd 0.029 K      hourly means: +0.02 .. -0.03
 stand summer   mean +1.099  K   sd 0.431 K      hourly means: +1.06 .. +1.16  (FLAT)
```

**A constant offset, day and night, present only under a canopy.**

⚠ **Correction to a first reading of this.** The flat diel profile was initially taken to rule out
a shortwave route ("an SW term would vanish at night"). That inference is wrong: the soil column
*integrates*. Per-day means over the month run 0.15 → 0.34 → … → 1.33 K and then plateau, while the
CAS difference — which has no memory — is stationary throughout. So the flat hourly profile is
thermal inertia responding to a net flux bias, and a daytime-only mechanism produces exactly this
signature. Diel flatness constrains nothing here; the *accumulation* is the real evidence, and what
it says is that there is a persistent **net energy flux bias** into the soil column.

The energy budget locates that bias (stand summer, refined refs, monthly means):

| | split-family | ARK-family | diff |
|---|---|---|---|
| Rnet | 126.22 | 125.96 | +0.27 W/m² |
| **H (sensible)** | 42.83 | 47.79 | **−4.96 W/m²** |
| LE | 42.59 | 43.24 | −0.65 W/m² |
| soil-top T | 292.80 | 291.71 | **+1.10 K** |

Bare summer: every one of these is ~0. With Rnet essentially equal, split routes ≈5.9 W/m² **away
from sensible heat and into the ground** — which is the soil warming, self-consistently. So the gap
is in the **canopy sensible-heat partition**, not in radiation.

**Three hypotheses tested and eliminated:**

1. *transp↔uptake seam (row 6) via soil moisture* — the `no_water` mask leaves the gap unchanged
   (1.14 vs 1.18 K) and `no_hydro` makes it *worse* (3.07 K). Not soil-moisture-mediated.
2. *Leaf longwave emission base* — split's leaf `veg_energy_diagnostic` passes `te` as `t_emit`
   where `surface_derivs` passes `tcas`. But `meds_fast_split.f90:492` is
   `te = tcas ; if (picard) te = t_emit(i)`, and these runs are non-Picard, so `te == tcas`.
   Identical. (It *would* differ under `integration_scheme = "picard"` — worth a separate look.)
3. *`src_frac` asymmetry* — `surface_derivs` scales `coh_qw`/`coh_qsoil`/`coh_transp` by
   `fro%src_frac` and split has no such factor at all, but `build_column_frozen` leaves `src_frac`
   at its 1.0 default (`meds_fast_ark.f90:1372`), so it is a no-op on both. Identical.

**A note on the attribution tool:** the `no_veg` mask cannot test this. `mask%veg_energy` only
freezes the *prognostic* leaf/wood stores (`meds_fast_split.f90:599,602`), and `--parity` runs
diagnostic leaf/wood — so the mask has no store to freeze and the run comes back byte-identical.
That is a vacuous result, not an exoneration. Testing the diagnostic leaf↔CAS balance needs an
instrumented single-step comparison of the per-cohort `dh` terms, not a mask.

**Status: characterised, not localised.** The next step is that instrumented single-step diff —
dump `coh_h`, `dh` per cohort, `h_coeff_f`, `g_tr_f`, `abs_sw`/`abs_lw` from both paths at an
identical state and difference them term by term.

### B-3. RK45 silently degrades to split in the cold dense-canopy cell

`work_rk45_rescue_site` = **0.497** in stand-winter (area-weighted month total ≈ 4 patch-steps),
0.000 in the other three, and *identical on both compilers*. So the rescue is deterministic, and
the regime is exactly the one the P6 design predicted. The stand-winter RK45 row is a partial
hybrid, not a clean RK45 comparand — which is precisely what this counter was added to reveal.

### B-4. The Phase-A compiler split is a saturation phenomenon, absent here

Same scheme, same config, ifx vs nvfortran, soil-top T RMSE: 5e-14 … 7e-13 K for every scheme in
three scenarios, worst case 2.3e-4 K (RK45, stand summer). `work_clamp_mass_site` is 0.0 in all
four scenarios — none of them saturates the column. So the 9-order compiler divergence measured in
Phase A is specific to the saturated regime and does **not** contaminate Phase B. Cost/robustness
(ifx, month totals): split 1488 steps / 0 rejections; ARK 2929–5498 steps / 50–1261 rejections;
RK45 1714–4138 steps / 1.7–587 rejections.

### B-5. Honest limits

One site, one month per cell, one forcing year, no confidence intervals — these are single
realisations, and only differences of more than ~2× should be read as real. GPP is identically
zero in both bare cells, so those rows carry no carbon signal by construction (as designed). The
50-yr spin-up ran on `split`, so the state restart starts every scheme from split's own converged
fast state; it is a shared starting point, so the comparison is fair, but split starts at home and
the B-2 sign convention should be read with that in mind.

## 3c. Phase C — landed so far

**C5 — RK45 bottom-BC guard.** `column_fast_step_rk45` now refuses an aquifer / bedrock /
Zeng-Decker bottom BC, as `column_fast_step_ark` always has. RK45 takes its ponding, aquifer and
water-table stores from the Act-1 scratch solve and integrates only θ, so those configurations ran
silently and wrong.

**C1 — the committed-state clamp is gone.** `rk45_column_step` no longer clamps `y_out`; the STAGE
clamps stay, because they bound throwaway inputs so `column_derivs` stays evaluable (notably
`ground_evaporation`'s fractional `pow()` needs a non-negative base) and a discarded stage state
breaks no books. An out-of-range commit is now handled the way every other bad step is — the
embedded error grows, the controller rejects and shrinks; and a floor-forced railed commit is
caught by `rk45_state_railed` at the dispatch, where the hybrid rescue redoes the `dt_fast` on the
split path. That path already existed; the clamp was pre-empting it with a silent correction.

Measured on the saturated column, 96 sub-steps, before → after:

| | ifx | nvfortran |
|---|---|---|
| unbookkept mass | 142.7 → **0.0** kg/m² | 142.6 → **0.0** kg/m² |
| unbookkept energy | 7.0e-5 → **0.0** J/m² | **1.513e5 → 0.0** J/m² |
| whole-column energy residual | 8.8e-7 → **6.0e-7** J/m² | **1.5e5 → 5.5e-7** J/m² |
| soil-surface peak | 329.4 → **297.0** K | — |

The compiler split on energy closure is gone with it (ifx 5.96e-7 vs nvfortran 5.52e-7), and the
saturated energy assertion tightened from `< 1e6` to `< 1e-3` J/m².

**Honest consequence, asserted in the test so it cannot drift:** θ is now committed slightly above
θ_sat (0.4536 vs 0.43) on the sealed saturated column. That is the same error the clamp was hiding,
now visible and on the books rather than silently paid for in fabricated mass. It is **C2's** to
remove — RK45 takes ponding/drainage/runoff from the frozen scratch solve while integrating its own
θ, so water that should leave as runoff has nowhere to go. The whole-column water residual improved
4.73 → 4.35 kg/m²; the rest is C2.

## 4. Phase plan

- **Phase A — instrument.** Done. Rows 9–11 closed.
- **Phase B — fidelity with demography frozen.** Done; see §3b. Headline: the schemes split into
  two families by semi-discretisation, and there is one soil-surface difference that survives 12×
  refinement, appears only under a canopy, and is not soil-moisture-mediated (**new row 12**,
  below).
- **Phase B — diagnose the fast loop with demography frozen** (`[run].slow_on = false`).
  Four scenarios, 2×2 over season × stand structure:
  1. near-bare start, winter, one month
  2. near-bare start, summer, one month
  3. `example_biophysics` 50-yr spin-up state, winter, one month
  4. same state, summer, one month

  Each run on ifx *and* nvfortran (the clamp path is compiler-split, and whether that split closes
  after C1 is itself a result). Cells 1–2 are scored on soil temperature / freeze front / θ, not
  GPP — a near-bare stand has essentially none. Cells 3–4 need an identical burn-in across schemes
  so the initialization transient is not mistaken for the answer.

  Target is **biophysics fidelity**, not explaining the 29-yr AGB divergence.
- **Phase C — parity fixes.** **C1 and C5 DONE** (see §3c). C2, C3, C4 open; C6 documented and
  deferred; row 12 characterised but not yet localised.
- **Phase D — the numerical comparison.** Accuracy × cost × conservation × robustness, both
  compilers, against both references described in §3.

## 5. Note on existing results

The 29/30-yr runs in `runs/ithaca_ark30/out/` predate PR #69 (`05c0a8d`), which fixed a forcing
recycle bug that phase-scrambled sub-daily shortwave in multi-decade recycled runs. Daily-mean
shortwave was correct throughout, so slow diagnostics looked healthy while the diel cycle was
wrong. Any conclusion drawn from those outputs about sub-daily behaviour — including the ~35% AGB
divergence — needs re-deriving on post-fix runs.
