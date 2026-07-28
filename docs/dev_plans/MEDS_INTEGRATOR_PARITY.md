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
| 12 | **sapflow advected enthalpy (`qwflux_wl`/`q_wood_net`) in the leaf+wood energy balance** — ED2's `qwflux_wl`; drives a +1.10 K soil-surface and ~5 W/m² sensible-heat offset that survives 12× refinement (§3b B-2) | **now present** via `q_extra` | present via `q_extra` | present via `q_extra` | physics | **FIXED** — ~80% of the soil bias / ~75% of the H bias removed; residual overshoot from the transp-for-sapflow freeze |

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

### LOCALISED — the split path omits the sapflow advected enthalpy

The instrumented single-step diff (one `dt_fast` from an identical state, `test_column_rk45`'s
fixture, midday) shows the divergence is in the **leaf energy balance**, not the ground:

```
                   split            rk45          diff
  leaf_temp    2.9081e+02      2.9827e+02     -7.458 K      <-- one step, dt = 900 s
  cas_temp     2.9811e+02      2.9633e+02     +1.778 K
  cas_shv      9.8998e-03      1.1148e-02     -1.248e-3
```

GPP is bit-identical on that same step (the shared frozen pre-pass agrees), so everything upstream
of the leaf balance matches. Refining at fixed total time separates the two error sources:

| dt_fast | Δleaf_temp | Δcas_temp |
|---|---|---|
| 900 s | −7.458 | +1.778 |
| 300 s | −2.760 | +0.196 |
| 100 s | −1.122 | +0.290 |
| 20 s | −0.602 | +0.305 |

Most of it is step-size error and collapses; a **non-zero floor of ~0.5 K in leaf temperature**
remains, matching the 0.32 K monthly-mean leaf difference measured independently in the Phase B
stand-summer cell. That floor is row 12.

**Mechanism.** `surface_derivs` passes `q_extra = fro%qwflux_wl(i)` to the leaf
`veg_energy_diagnostic` and `q_extra = fro%q_wood_net(i)` to the wood one
(`meds_fast_time_derivs.f90:121,144`). These are the **sapflow / xylem advected enthalpy** terms —
ED2's `qwflux_wl`/`qloss` — built in `build_column_frozen` (`meds_fast_ark.f90:1495-1502`) and
therefore shared by ARK and RK45. **The split path neither computes nor consumes them**: `q_extra`
appears nowhere in `meds_fast_split.f90`, whose only mentions of `qwflux_wl` are comments noting
the coupling is missing (`:877`, `:1055`).

`q_extra` enters the `dt_temp` balance exactly like absorbed shortwave, so it shifts the leaf's
equilibrium temperature and hence `dh` — the sensible-heat term — which is precisely the ~5 W/m²
partition shift measured above. Magnitude checks out: `sapflow · u_liq` with sapflow ~ transpiration
and `internal_energy_liquid ~ 1.0 MJ/kg` gives tens of W/m², against a leaf `h_coeff` of order
50–100 W/m²/K ⇒ a few tenths of a K. Canopy-only (no plant, no sapflow), and a different *equation*
rather than a different discretisation, so refinement cannot remove it — every observed property.

**Which path is right: ARK/RK45.** ED2 carries this coupling (`rk4_derivs.f90:2122`), MEDS lacked it
project-wide, and `MEDS_ED2_RK45_DESIGN.md` P2 added it — on the ARK/RK45 side only. Split is the
one missing a real term.

### FIXED — split now carries the sapflow advected enthalpy

Implemented as the telescoping triple, so conservation depends only on all three terms sharing one
frozen flux — not on *which* flux:

| store | term | where |
|---|---|---|
| leaf | `+qwflux_wl` | `q_extra` on the leaf `veg_energy_diagnostic` |
| wood | `+(qloss − qwflux_wl)` | `q_extra` on the wood `veg_energy_diagnostic` |
| soil | `−qloss` | `root_heat_sink = (coh_qsoil + qloss_total) · profile` |

**Which frozen flux, and the approximation taken.** ARK freezes `solve_plant_water`'s exact-solve
sapflow in its Act-1 pre-pass. Split solves hydraulics *after* the leaf balance, so that number
does not exist at the call site. Rather than carry a lagged sapflow as new per-cohort lockstep
state, the leaf balance is probed once with no advective term to get a state^n transpiration and
the triple is frozen on that. Sapflow and transpiration differ by the leaf water storage change,
second-order over a `dt_fast`. The probe is closed-form algebra (`veg_energy_diagnostic` is
`elemental pure`), not a second hydraulics solve.

**Result** — the structural family gap at the refined reference (dt=150), before → after:

| | b2 bare summer | b3 stand winter | b4 stand summer |
|---|---|---|---|
| soil-top T, RMSE | 0.0286 → **0.0286** | 0.623 → **0.442** | 1.180 → **0.221** |
| soil-top T, mean bias | +0.0006 → **+0.0006** | +0.278 → **+0.162** | **+1.099 → −0.210** |
| H, mean bias | +0.081 → **+0.081** | −0.223 → +0.319 | **−4.956 → +1.308** |
| CAS T, RMSE | 0.0208 → **0.0208** | 0.160 → **0.123** | 0.314 → *0.610* |

- **The bare-summer column is byte-identical.** With no transpiration `qloss_total ≡ 0` and the
  code reduces to the previous form exactly — the control the design intended.
- **~80% of the soil-surface bias and ~75% of the sensible-heat bias are removed** in the worst
  cell, and both now sit slightly past zero.

**The planned refinement was tested and REJECTED on measurement.** The residual sign flip was
attributed to the transpiration-for-sapflow freeze over-correcting, with the fix being to freeze on
the exact-solve sapflow (requiring split's hydraulics moved ahead of the leaf balance, or a
persisted per-cohort sapflow). Instrumenting the two fluxes over a stand-summer month says
otherwise:

```
  mean sapflow 1.7886e-5   mean transp 1.7826e-5   -> 0.33% apart
  median ratio sapflow/transp = 1.0032
  per-step ratio ranges widely (-0.69 .. 8.87) at low flux -- dawn/dusk store charging --
  but those are moments when the flux itself is ~0, so they carry no energy
```

A frozen flux that is right to 0.33% in the mean cannot account for over-correcting a 1.3 K
adjustment by 19%. So the restructure would buy essentially nothing, and it is **not** being done.
The more likely reading of the residual is simply that it is *the next difference down*, of
opposite sign, previously masked by the much larger row-12 term — not an artefact of the freeze.

**Known deviation from ARK, deliberately left:** ARK builds `qloss` from `uptake_frozen`
(soil→wood) while this implementation uses the same probe flux for both `qloss` and `qwflux_wl`.
The two differ by the wood water storage change, which the measurement above bounds at well under
1%. Fixing it needs the same unavailable ordering, for the same negligible return.
- **CAS-T RMSE in b4 got worse** (0.314 → 0.610) while its soil bias collapsed. Reported rather
  than buried: the leaf equilibrium moved, so the leaf→CAS sensible flux moved with it. Whether
  that is the over-correction or a second effect is not established.

Golden anchors rebased (approved): CAS 292.543227 → **292.660995** (+0.118 K), soil-surface
292.884295 → **291.718995** (−1.165 K). The soil move matches the independently measured +1.10 K
family gap in sign and magnitude, which is the cross-check that the right term was added. Anchors
are now named constants with the realized values printed beside them, so the next rebase is a
copy-paste rather than a bisect. ifx Release / ifx Debug / nvfortran multicore all 37/37.

**Superseded scoping note** (kept for the reasoning): Split already has the ingredients: it calls
`solve_plant_water_batch` itself and so has `flux%sapflow`. But `q_extra` is only half of a pair —
the kernel's own doc-comment is explicit that it is "an internal soil↔leaf transfer already debited
from the soil's own `root_heat_sink`", and it is deliberately excluded from `drnet` so the
whole-column identity closes. So the change is: compute `qwflux_wl`/`q_wood_net` on the split path
from its own sapflow, pass them as `q_extra` to both `veg_energy_diagnostic` calls, and verify the
matching soil-side debit (split's `coh_qsoil`) pairs with them exactly. This touches the validated,
golden-anchored split energy balance and **will move `test_picard_coupling`'s anchors**, so it wants
its own verification pass: whole-column energy closure first, anchors re-pinned second.

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


### B-6. Post-fix re-measurement (after C1 + C5 + the row-12 fix)

The Phase B matrix re-run on the same four scenarios. Soil-top T, RMSE vs the refined reference —
the quantity that was structurally stuck:

| scenario | split BEFORE | split AFTER | ark | rk45 |
|---|---|---|---|---|
| b3 stand winter | 0.716 | **0.542** | 0.232 | 0.181 |
| b4 stand summer | 1.397 | **0.360** | 0.357 | 0.275 |

**Split has joined the family on soil temperature**: in the worst cell it went from 3.9× worse than
ARK to indistinguishable from it. Pairwise, `split~ark` soil-top T fell 1.185 → 0.628 and
`split~rk45` 1.202 → 0.518 in b4.

Two honest caveats:

- **CAS temperature in b4 got worse**, both pairwise (`split~ark` 0.746 → 1.046) and against the
  reference (1.229 → 1.355). The leaf equilibrium moved, so the leaf→CAS sensible flux moved with
  it. Net the fix is clearly positive — soil improved ~4×, CAS degraded ~1.1× — and it restores a
  term ED2 and MEDS's own ARK/RK45 paths already carry, but the CAS direction is unexplained and
  is the first thing to look at when refining the frozen sapflow.
- The **bare cells' `split~ark` pairs are unchanged** (byte-identical, as designed), but every pair
  *involving rk45* moved — that is C1 changing RK45's trajectory, not the row-12 fix.

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
