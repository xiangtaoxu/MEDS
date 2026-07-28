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
| 5 | Non-free-drain bottom BC / Zeng–Decker | present | hard `error stop` (`meds_fast_ark.f90:838`) | **no guard — runs silently** | physics | **C5** — give RK45 the same guard |
| 6 | transp↔uptake seam | re-solved inside the Picard iterate (realized, supply-limited) | frozen from state^n | frozen from state^n | numerics | measure in Phase B; knob if it matters |

Row 1 is the one that most affects existing results: `cas_condensation` defaults `.true.`, so the
30-yr split run and the 29-yr RK45 run in `runs/ithaca_ark30/` were **not the same model**. Worse,
both ARK and RK45 book `sf%cond` into `whole_wat_out` (`meds_fast_ark.f90:224`,
`meds_fast_rk45.f90:217`) — condensed dew/fog leaves the tracked column entirely rather than
reaching soil or canopy. Split neither condenses nor loses it. The two paths are wrong in
different directions.

### Class 2 — unbookkept corrections (state edited outside the conservation ledger)

| # | difference | plan |
|---|---|---|
| 7 | RK45 clamps the **committed** state (`meds_fast_rk45.f90:199`), not just throwaway stage inputs. `clamp_theta` moves θ with no mass debit; `clamp_soil_energy` then re-derives T at the new water mass. | **C1 — remove the commit clamp** (decision taken 2026-07-28) |
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

## 4. Phase plan

- **Phase A — instrument.** Done. Rows 9–11 closed.
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
- **Phase C — parity fixes.** C1, C2, C3, C4, C5 above; C6 documented and deferred.
- **Phase D — the numerical comparison.** Accuracy × cost × conservation × robustness, both
  compilers, against both references described in §3.

## 5. Note on existing results

The 29/30-yr runs in `runs/ithaca_ark30/out/` predate PR #69 (`05c0a8d`), which fixed a forcing
recycle bug that phase-scrambled sub-daily shortwave in multi-decade recycled runs. Daily-mean
shortwave was correct throughout, so slow diagnostics looked healthy while the diel cycle was
wrong. Any conclusion drawn from those outputs about sub-daily behaviour — including the ~35% AGB
divergence — needs re-deriving on post-fix runs.
