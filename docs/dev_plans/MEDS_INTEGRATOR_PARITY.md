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
| 1 | CAS supersaturation sink (`TAU_COND`) | **now present** (exact relaxation of the state^n excess) | present (per-stage rate) | present (per-stage rate) | physics | **C3 DONE** — same model, different quadrature; dropped from the `--parity` preset |
| 1b | condensate was DELETED — `sf%cond` booked into `whole_wat_out` on every path, so dew/fog left the tracked column | **deposited into soil layer 1** | same | same | physics | **DONE** — paired mass+enthalpy transfer, identical destination on all three paths |
| 2 | Snow / temporary surface water | **shared stage** (`meds_fast_snow`) | **same stage** | **same stage** | physics | **C4 DONE** — all three run identical snow and close both ledgers with a pack; dropped from `--parity` |
| 2b | Canopy surface water (interception film / film-evap / dew, `canopy_water_on`) | present | present (`advance_surf_water_full`) | present (native `column_state_t`) | — | **already unified** — verified, no work needed |
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
| 8 | RK45 commits its **own** integrated θ but took pond / drainage / runoff from the **frozen scratch** solve. | **C2 DONE** — drainage b-weighted from RK45's own stages; residual clip routed to the pond; pond rebuilt from RK45's own trajectory. θ overshoot gone; **both regimes now close to machine precision** |

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


### C3 — condensation unified (and canopy surface water confirmed already unified)

**Canopy surface water needed no work.** It is already on all three paths: split operator-split,
ARK via `advance_surf_water_full`, RK45 natively in `column_state_t`. Verified, recorded as row 2b.

**The CAS supersaturation sink is now on split too**, closing row 1. It is folded into
`src_vap`/`src_enth` *before* the shared CAS box, exactly as `surface_derivs` does — that placement
is what keeps every downstream ledger consistent for free. A first attempt applied it *after* the
box instead, which satisfied the state update while silently breaking the per-kernel `cas_water` /
`cas_energy` budgets; the tests caught it.

**Different quadrature, deliberately.** `surface_derivs` uses a rate,
`(wcap/TAU_COND)·max(0, q−qsat)`, re-evaluated per stage and resolved by the adaptive march. Split
takes one step of `dt_fast`, and `dt_fast/TAU_COND = 1800/300 = 6` — that rate applied over the
whole step would remove ~6× the available excess and drive the air far below saturation. Split
therefore uses the exact relaxation of the state^n excess, `(1 − exp(−dt/τ))`, which can never
remove more than the excess present and agrees with the rate form as `dt → 0`. Same model, two
quadratures: a numerics difference the harness measures, not a physics one refinement cannot close.

**Verified to actually engage** (the unit fixtures never supersaturate, so they cannot show this —
a coverage gap worth closing). Forced b3 stand-winter, split, condensation on vs off:

```
  cas_temp_site       mean +0.041 K    max |Δ| 1.537 K
  soil_temp_top_site  mean +0.018 K    max |Δ| 0.576 K
  le_flux_fast        mean -0.317 W/m2 max |Δ| 6.715 W/m2
```

**Success criterion met.** A *default-config* comparison (condensation on, as production runs have
it) now matches the parity comparison: split~ark CAS 0.434 vs 0.447, soil-top 0.473 vs 0.471. The
sink is no longer a model-family difference, and `fast.cas_condensation` has been dropped from the
`PARITY` preset.

**Row 1b, newly separated:** all three paths still book `sf%cond` into `whole_wat_out`, so
condensate leaves the tracked column rather than landing on the canopy film or the ground. That was
previously entangled with row 1; it is now a *shared* defect, identical on every path, so it no
longer contaminates scheme comparisons. Fixing it needs a paired mass+enthalpy transfer and a
destination decision (canopy film when `canopy_water_on`, else ground/pond) — its own change.


### Row 1b — condensate is no longer deleted

Every path booked `sf%cond` (and its enthalpy) into `whole_wat_out`/`whole_enth_out`, so dew and fog
left the tracked column entirely. All three now **deposit it into soil layer 1** as a paired
mass + enthalpy transfer valued at the CAS temperature it condensed at, and the corresponding
boundary terms are gone — the whole-column ledger closes with no boundary term for it at all.

ARK needed a plumbing change first: its `sf%cond` was summed into `whole_wat_out` alongside the
atmospheric vapour flux, where it could not be told apart. `stage_bflux_t`/`column_bflux_t` gained a
`whole_cond` slot so the b-weighted accumulation carries it separately. RK45's `bw_cond` was already
separable; split's `cond_rate` is a local.

**Why soil layer 1 rather than the canopy film.** Layer 1 exists on every path in every
configuration and carries a full mass+energy state, so the transfer is exact. The canopy film store
is opt-in (`canopy_water_on` defaults false, so it is absent from production runs) *and* its
enthalpy is referenced at a fixed `rain_temp` rather than tracked, so depositing there would
introduce an enthalpy mismatch at the seam. Routing dew onto foliage when `canopy_water_on` is the
more physical choice and is the natural refinement — it needs the film's enthalpy reference resolved
first, and only bites once that store is enabled. **Tracked as GitHub issue #74**, together with the
coverage gap that no unit fixture supersaturates.

The deposit happens after the step's soil solve, so the water joins the column one step before it
diffuses or drains. That is a one-step lag on a dew-sized flux; injecting it mid-solve would reopen
the time-level inconsistency PR #71 closed.

**Verified.** The unit fixtures never supersaturate, so they cannot exercise this. Forced b3
stand-winter, all three schemes, condensation on vs off, with `energy.debug_error = .true.` so any
budget breach *halts*: all six runs completed with zero budget messages. Column soil water with
condensation on, minus off — positive means the water now lands instead of vanishing:

```
  split   mean +8.86e-4    final +8.41e-4
  ark     mean +1.17e-3    final +1.05e-3
  rk45    mean +1.29e-3    final +1.26e-3
```

Same sign and same order on all three, which is the point.


### C2 (partial) — drainage now follows RK45's own θ

`soil_water_time_deriv` has always returned a **state-dependent** `drainage_rate` computed from the
stage's own θ, and the mass side used it — but the whole-column ledger booked the *frozen*
`fro%drainage`, and the soil-energy bottom face advected enthalpy on the frozen value too. So the
two halves of the same face disagreed by exactly the amount RK45's trajectory departed from the
Act-1 scratch solve. Both now use the b-weighted `drainage_rate`, with the same b-vector as the
state commit. `column_derivs` evaluates the soil-**water** tendency before the soil-**energy** one so
the enthalpy can ride the flux the mass side just produced.

| | before | after |
|---|---|---|
| unsaturated whole-water residual | 1.012e-6 | **2.433e-13** kg/m² |
| θ overshoot above θ_sat (saturated) | 0.4536 | **0.4382** (θ_sat = 0.43) |
| saturated whole-water residual | 4.35 | *5.86* kg/m² |

**The unsaturated case now closes to machine precision** — drainage was the last inconsistent term
there, and that is the regime every production run lives in (the 29-yr runs never saturate).

**The saturated case got ~35% worse, and that is expected.** Before, drainage *and* runoff both came
from the frozen scratch solve — mutually consistent with each other, even though neither matched the
committed θ. Now drainage matches θ and runoff does not, so the imbalance moved rather than shrank.
That is the cost of fixing half a pair, and it is worth paying at 1e-6 → 1e-13 in the regime that
matters. Closing it needs the ponding/runoff half — **GitHub issue #75**.

**Two things this surfaced, both filed:**

- `test_column_rk45`'s saturated case calls itself a *"sealed bedrock column"* and imports
  `SOIL_BC_BEDROCK`, but **never sets `bottom_bc`** — the default is `SOIL_BC_FREE_DRAIN`, so it is a
  free-draining column driven to saturation by 29 mm/hr precip. The unused import is the fingerprint
  of an intent that was never wired.
- `test_column_derivs`' soil-heat wiring check omitted the bottom-face drainage enthalpy from its
  reproduction and passed anyway, because the fixture's frozen `fro%drainage` was zero. A check that
  omits a term is only green while that term is zero. It now carries the term.





### C2 second pass — the residual saturation clip

Two findings narrowed this a lot more than the original scoping suggested:

- **Dunne runoff is unreachable under RK45.** `f_sat` is nonzero only for `SOIL_BC_AQUIFER`, which
  `column_fast_step_rk45` hard-errors on since C5. RK45's runoff is therefore *purely* pond overflow
  — no `f_sat` logic needs reproducing inside the stages.
- **Frozen infiltration is correct, not a compromise.** `column_hydrology_flux` computes `infl` from
  state^n *before* its own implicit solve, so the split path freezes it identically. There was never
  a discrepancy there.

That leaves one genuinely state-dependent term: the **post-solve saturation clip**. RK45 commits its
own θ, which can sit above θ_sat by however far its trajectory departed from the scratch's — exactly
the overshoot C1 exposed once the unbookkept commit clamp was removed. It now clips its own excess
per layer, sheds that layer's enthalpy at *its own* temperature, and routes the water to the ponding
store, overflowing to runoff. This is a clamp **with** bookkeeping, which is the distinction C1 drew:
the water is real and has somewhere to go, so moving it is a transfer rather than a silent
correction.

| | value |
|---|---|
| θ overshoot above θ_sat | 0.4382 → **0.43000** (exactly θ_sat) |
| saturated whole-water residual | 5.86 → **3.33** kg/m² (was 4.35 pre-C2) |
| unsaturated whole-water residual | **2.7e-13** kg/m² |

**The hypothesis was tested and confirmed.** The pond was composed as `fro%w_surface1 + clip_rk`,
but `fro%w_surface1` is the *scratch* solve's end-of-step pond and already contains the scratch's own
clip — mass RK45's θ never shed, since the explicit tendency carries no clip term. That water was
counted twice. Rebuilding the pond from RK45's own trajectory,
`w_surface0 + (precip_ground − infl)·dt + clip_rk`, with runoff derived from *that* pond rather than
`fro%runoff_surf`, closes it:

| | saturated whole-water residual |
|---|---|
| pre-C2 (all three terms frozen) | 4.35 kg/m² |
| drainage made state-consistent, runoff not | 5.86 (imbalance *moved*) |
| residual saturation clip added | 3.33 |
| **pond rebuilt from RK45's own trajectory** | **2.96e-13** |

Both the saturated and unsaturated regimes now close to machine precision, and the test asserts
`n_fail == 0` rather than a tolerated bound.

**Caveat carried forward:** the clip itself is a symptom repair, not the physical mechanism — the
design doc's Neumann→Dirichlet ponded-surface switch is the real fix and is still deferred. Its
assumptions are recorded in **issue #78**.

### C2 third pass — the interior advective faces (#78 item 3) DONE

C2 fixed the *bottom* face and the ponding store. Two more borrowed-flux terms survived on RK45, and
they turned out to be a matched cancelling pair — which is why the fixture looked healthy:

| term | what RK45 did | what its own θ does |
|---|---|---|
| interior faces `eforc%w_flux(k)` | advected on the scratch solve's frozen `hflux%w_flux` (2.2 kg/m² per step) | moves ~5.0 kg/m² per step |
| stage clip enthalpy `fro%clip_enth(k)` | paid ~2.6×10⁶ J/m² per step for the scratch's clip mass | sheds only ~0.3 kg/m² |

Both are correct for the **ARK**, which operator-splits soil water out of its stages and commits the
scratch solve's θ verbatim — the scratch's fluxes genuinely are its own. Both are wrong for RK45,
which integrates θ itself. And the two errors were ~2.6×10⁶ J/m² per step *each*, with opposite sign:

| saturated soil-surface peak | |
|---|---|
| both defects present (they cancel) | 285.3 K |
| borrowed clip enthalpy removed alone | **345.0 K** |
| interior faces on RK45's own θ, clip enthalpy removed | **295.1 K** |

`soil_water_time_deriv` now exports the interior face flux it built `dtheta_dt` from, and
`column_derivs` advects soil enthalpy on that. `fro%clip_enth`/`floor_enth`/`w_flux_frozen` remain on
`column_frozen_t` and remain correct — they are now **ARK-only**.

Two things are worth carrying forward from this one:

* **A whole-column ledger is blind to this defect class.** Enthalpy placed in the wrong *layer* still
  sums correctly against the boundary; the error is purely vertical. Both budgets closed to machine
  precision throughout, on every version above including the 345 K one. What exposed it was a
  temperature that stopped being plausible, and what localized it was printing the terms.
* **The change is a measured no-op on every Ithaca cell** (RK45 soil-temperature peak 276.49 →
  276.50 K on the winter stand) because those columns never saturate: rain there is far below `K_sat`,
  so the scratch's faces and RK45's own faces nearly coincide. It only bites at 29 mm/h. `split` and
  `ark` output is **byte-identical** before and after, as designed.

A saturation-**relief tendency** (`−max(0, θ−θ_sat)/τ` inside the stages, so θ never oversaturates and
the post-commit clip becomes vestigial) was built, measured, and **rejected**: with the faces fixed it
changed nothing observable — 295.10 vs 295.05 K, whole suite green, zero stage clamps either way. Its
only real effect was a ~27× smaller in-stage excursion past θ_sat, which is issue #78 **item 2**'s
concern (out-of-domain constitutive evaluation), not item 3's. Adding a timescale parameter and a
tendency term for no measurable gain is how the borrowed `clip_enth` got there in the first place.

RK45 also gained the **θ_res floor** on its committed state, since removing the `fro%floor_enth`
compensation left that edge of the constitutive domain guarded by nothing. It creates water, so it is
booked as an explicit boundary input rather than hidden in a residual. Honest limitation: a drydown
test bottoms out at θ = 0.0807 against θ_res = 0.078 — both sinks self-limit first (ψ-limited root
uptake, capped ground evaporation) — so the floor never fires and **its ledger booking is
unexercised**; a mutation that unbooks it passes. The test asserts the *invariant*, which is the right
assertion whichever mechanism enforces it.

### C4 DONE — all three integrators run the same snow

`advance_snow_stage` (`src/driver/meds_fast_snow.f90`) is now called by all three paths: split
directly, ARK and RK45 through `build_column_frozen`. `surface_derivs` blends the snow surface, and
both ledgers carry the pack's mass and energy stores.

**RUN 8, all three schemes, 24 h with a 60 kg/m² pack under snowfall:**

| | final swe | melt | worst water resid |
|---|---|---|---|
| split | 56.058 | 3.942 | 2.68e-13 |
| ARK | 55.955 | 4.045 | 2.68e-13 |
| RK45 | 55.968 | 4.032 | 2.37e-13 |

Both whole-column ledgers close at `n_fail == 0` on every scheme, melt totals agree within 3%, and
the numbers are identical on ifx Release, ifx Debug and nvfortran multicore. `fast.snow_on` is
dropped from the `--parity` preset.

**Four bugs the activation had to clear, each caught by an assertion rather than by inspection:**

1. **Call placement.** The stage must run *after* `column_prepass` — it needs `tcas`/`qcas`,
   `rho`/`press` and `aero%ggnet`, undefined earlier in `build_column_frozen`. Placing it before
   gave a NaN pack.
2. **`bio` had to become `intent(inout)`** in `build_column_frozen` — the stage advances the pack and
   hands melt enthalpy to the soil.
3. **Precip routing is a branch, not an addition.** `snow_accumulate` has already taken `snowf` *and*
   `precip` into the pack, so only meltwater may reach the ground; adding throughfall on top
   double-counts.
4. **`forc%snowf` was missing from `w_in`** on both ARK and RK45 (split has always carried it).
   Snowfall is a boundary water input landing in the pack, so without it the ledger saw mass appear
   with no source and leaked exactly `snowf·dt` per step. **This was the water-closure bug.**
5. **The soil baseline needed rebasing by the melt enthalpy.** Split snapshots `e_soil0` *before* the
   snow stage; ARK/RK45 build `y` *after* it, so the melt energy sat in the soil baseline *and* in
   `snow_enth0` — counted twice. `snow_stage_t%melt_enth` reports it and both paths subtract it.
   **This was the energy-closure bug**, and it fired only on melting steps (46 of 96), which is what
   pointed at it.

A methodological note worth keeping: at one point the ARK/RK45 ledgers "closed" because the test had
switched snow back off before those runs — a silently snow-free run passes every budget assertion
trivially. The melt-total-agrees-with-split assertion is what makes that visible, and it is the
reason to compare schemes against each other rather than only against zero.

### C4 progress — shared stage landed, activation does not conserve yet

**Landed and verified** (snow-off bit-identical, snow-ON RUN 8 closing to machine precision, three
compilers 37/37):

- `src/driver/meds_fast_snow.f90` — the shared pre-column stage, `advance_snow_stage` + the
  `snow_stage_t` frozen-output bundle, extracted verbatim from `column_fast_step`. Split now calls
  it and consumes the bundle; its own snow locals are gone.
- `surface_frozen_t` carries the six snow fields (`snowfac`, `h_snow`, `le_snow`, `g_base_snow`,
  `subl_rate`, `ground_rad`) plus the five ledger terms, all defaulting to zero.
- `surface_derivs` blends the snow surface — `h_ground = h_snow + (1-snowfac)·h_bare`, etc. Written
  so zeros reduce it *exactly* to the pre-C4 expression, which makes snow-off bit-identity a
  structural property rather than something each scheme must re-verify.

**NOT landed: the ARK/RK45 activation.** Wiring `advance_snow_stage` into `build_column_frozen` was
attempted and **backed out** — with a pack present the whole-column ledgers failed on essentially
every sub-step, and the pack diverged ~89% from split's. It is not committed, because a wired-but-
non-conserving snow path is worse than an honestly absent one: `snow_on` would silently produce
wrong ARK/RK45 results instead of no snow.

**What the attempt established, so the next pass starts from evidence rather than a blank page:**

1. **Call placement matters and is non-obvious.** The stage must run *after* `column_prepass` — it
   needs `tcas`/`qcas` (CAS state), `rho`/`press`, and `aero%ggnet`, none of which exist earlier in
   `build_column_frozen`. Placing it before produced a NaN pack. Split calls it after its own
   `column_prepass` for the same reason.
2. **`build_column_frozen` must take `bio` as `intent(inout)`** — the stage advances `bio%snow` and
   hands melt enthalpy to `bio%soil_e%soil_energy(1)`.
3. **Precip routing is a branch, not an addition.** Split does
   `if (snow_exists) precip_ground = melt_rate` *else* `throughfall_total` — the pack has already
   taken `snowf` *and* `precip`, so adding melt on top of throughfall double-counts. Fixing this
   alone did **not** close the ledgers, so at least one further term is missing.
4. The remaining imbalance is **not** the precip double-count and has not been isolated. The RUN 8
   assertions extended to ARK/RK45 (physicality, melt path, both ledgers, and pack-mass agreement
   with split to 1%) are the right acceptance test and were written; they are not committed since
   they currently fail.

### C4 — prerequisite landed, migration not started

Sizing C4 properly turned up a blocker worth more than the migration itself: **nothing in the suite
ever set `ccfg%snow_on`.** `test_snow.f90` exercises the snow *kernels* standalone, but the snow
*coupling* inside `column_fast_step` — accumulate → surface balance → meltwater drain →
snowfac-blended ground → sublimation into the CAS → the `swe`/`snow_acc_enth` terms in the
whole-column ledgers — had no integration test at all.

That makes the migration unsafe in a specific way: the obvious success criterion, "snow-off stays
bit-identical", only proves the OFF path survived, and the OFF path is the half that cannot break.

`test_column_dynamics` RUN 8 is that net. A seeded 20 kg/m² pack under sub-freezing air and light
snowfall, marched 24 h on the split path, asserting the pack and soil stay physical, the
melt/sublimation path is live, and **both whole-column ledgers close at `n_fail == 0`**. First
measurement of split's snow coupling closure:

```
  worst whole-column resid:  energy 5.49e-7 J/m2   water 2.25e-13 kg/m2
  pack: 20 -> 1.44 kg/m2, T_snow pinned at 273.16 K (the melting plateau)
```

It conserves to machine precision. That was believed but never asserted.

**The migration itself is [issue #76](https://github.com/xiangtaoxu/MEDS/issues/76).** Scope, as measured rather than estimated: the snow coupling
reaches seven places in split (pre-column advance, `rain_temp` routing, precip routing,
`hforc%snow_free_frac`, the ground blend, `src_vap` sublimation, three ledger terms), and enabling
ARK/RK45 additionally needs new `surface_frozen_t` fields, a blend inside `surface_derivs`, and snow
terms in two more ledgers. `surface_derivs` has no snow terms whatsoever today, so this is *adding
snow coupling to the ARK/RK45 surface block*, not relocating code.


## 3d. Phase D — the numerical comparison

Run on `b4_stand_summer` (the strongest leaf↔CAS coupling of the four Phase B cells), all cells
under `--parity`, scored on the **hourly** FAST stream. Both references were checked for provenance
before anything was read from them — `work_rk45_rescue_site = 0`, no clamp activations, and E2's
reference resolved 99200 sub-steps with **zero rejections**. That check is the whole reason the
Phase A counters exist.

### D-1. At fixed `dt_fast`, the stepper choice matters a great deal

Every cell at `dt_fast = 1800 s` — identical Category-0 freeze cadence — against `rk45_tight` at the
**same** `dt_fast`. This isolates time-stepping error with the dominant error source held common:

| scheme | CAS T | soil-top T | LE | GPP | wall (net) | steps | rejections |
|---|---|---|---|---|---|---|---|
| **ark** | **0.039** | **0.282** | **0.94** | **0.085** | 3.27 s | 5492 | 1253 |
| rk45 | 0.442 | 0.272 | 8.14 | 0.298 | 3.87 s | 5033 | 913 |
| split | 0.831 | 0.533 | 22.9 | 0.489 | **2.51 s** | 1488 | 0 |
| picard | 0.899 | 0.466 | 17.4 | 0.663 | 5.17 s | 1488 | 0 |
| ark_fixed | 1.998 | 0.751 | 155.7 | 1.089 | 2.91 s | 1488 | 0 |

**ARK dominates the Pareto front**: 21× more accurate than split on CAS temperature for 1.3× the
wall time. `ark_fixed` — one ESDIRK step per `dt_fast`, adaptivity off — is by far the worst on every
metric, which is the cleanest evidence available that the adaptive march is doing real work rather
than just spending time.

**Unexplained, and worth a look:** RK45 is 11× *less* accurate than ARK here despite the reference
being `rk45_tight` — its own integrator at tight tolerance. Two concrete leads: RK45's error norm
includes mass (`with_mass = .true.`, where ARK excludes it), which dilutes the temperature error
control; and its controller exponent uses `p_order = 4` against ARK's embedded `p = 1`, so it resizes
steps far less aggressively for the same estimate.

### D-2. Refining `dt_fast` does NOT converge the solution

> ⛔ **D-2 AND D-3 ARE WITHDRAWN — the reference cell was mis-timed. See §3e.** The tables in this
> subsection and the whole of D-3 are an artefact of `numerics_sweep.py` choosing a reference
> `dt_fast` (27 s) that does not divide 3600 s, which made the reference's "hourly" FAST records
> 3591 s long and drifting. With an hour-aligned reference **every scheme converges cleanly.** The
> text below is kept only so the correction is legible. D-1 is unaffected (its reference runs at
> `dt_fast = 1800 s`, which divides 3600 exactly), and so is Phase B (refined references at 150 s).

Same cells swept over `dt_fast`, against `rk45_tight` at 27 s:

| scheme | 1800 | 900 | 450 | 225 |
|---|---|---|---|---|
| split, CAS T | 0.919 | 0.707 | 0.774 | 0.817 |
| ark, CAS T | 0.975 | 0.675 | 0.771 | 0.829 |
| rk45, CAS T | 0.809 | 0.673 | 0.774 | 0.830 |

The error falls from 1800→900 and then **rises again**. That is not convergence.

**Control (D-3): split against its OWN refined self** (`split` @ 27 s), which removes any
family/model difference at the reference:

| dt | 1800 | 900 | 450 | 225 |
|---|---|---|---|---|
| CAS T | 0.836 | 0.677 | 0.806 | 0.874 |

Identical shape and magnitude. **So the non-convergence is intrinsic, not an artefact of comparing
across scheme families.** Refining `dt_fast` by 8× does not reduce the hourly-scored CAS-temperature
error below ~0.7 K, and past 900 s it makes it worse.

Cause not established — tracked as **GitHub issue #79**, with a three-step test plan. The
leading candidate is the forcing quadrature: `forcing_sample_frac`
samples met at the sub-interval midpoint, so different `dt_fast` values sample a different set of
points from the same hourly file, and the shortwave reconstruction is interval-mean-conserving
rather than pointwise. Averaging that to hourly output does not commute with the nonlinear canopy
response. This is a hypothesis, not a diagnosis.

### D-4. What this means for scheme selection

The naive lever — "take smaller time steps" — **is not available here**, so the practical choice is
the integrator at production `dt_fast`. On that axis ARK is the clear recommendation: it buys an
order of magnitude in accuracy for ~30% wall time, and it is the only scheme whose stepper error
(0.039 K) sits well below the ~0.7 K floor rather than at it.

Split's error at `dt_fast = 1800` is essentially *all* time-stepping error (stepper 0.831 vs total
~0.84), while ARK's is essentially all floor (stepper 0.039 vs total ~0.97). That is the sharpest
statement of the trade: **split is limited by its stepper, ARK is limited by everything else.**

### D-5. Honest limits

One site, one month, one stand, one forcing year; single realisations with no confidence intervals;
`--repeats 3` min-of-N on wall clock only. Scored as hourly TMEAN, which suppresses sub-hourly
structure. D-1 and D-2 use **different references** and answer different questions — their numbers
must not be added or compared directly. The ~0.7 K floor is a property of this configuration and has
not been traced to a mechanism.

## 3e. Phase E — the D-2 reference was mis-timed; everything converges (2026-07-29)

### E-0. The instrument fault

`numerics_sweep.py` picks the reference cell's `dt_fast` as `min(dt) // ref_refine`, then snaps it
down until it divides `dt_slow`. With `--dt 1800 900 450 225 --ref-refine 8` that yields **27 s**:
86400 / 27 = 3200, so it passed the one check there was.

It fails a second check nobody had written down. The FAST output tier is closed every
`fast_interval_steps` **sub-steps**, and the harness sets that to `round(3600 / dt_fast)` so every
cell lands on a common hourly axis. 3600 / 27 = 133.33 → **133**, and 133 × 27 = **3591 s**. The
reference's records are therefore 9 s short of an hour and drift: over a 744-hour month the drift
accumulates to ~1.9 h, the file holds **746 records for 744 hours**, the `(day, hour)` key that
`parity_fidelity` aligns on collides twice, and the record's minute stamp wanders 0 → 59 and wraps.

Direct evidence, `e3_selfconv_split`, first/last FAST record stamps as `(day, hour, minute)`:

```
  ref_split_27_full    nrec=746  unique(day,hour)=744  dupes=2
     first: (1,0,0) (1,1,0) (1,1,59) (1,2,59) ...      <-- 3591 s records, already drifting
     last : (31,20,8) (31,21,8) (31,22,8) (31,23,8)
  split_225_full_i     nrec=744  unique=744  dupes=0   minute stamp pinned at 1
  split_1800_full_i    nrec=744  unique=744  dupes=0   minute stamp pinned at 15
```

Every cell was compared against a reference sampled over a progressively mis-centred window. Since
the mis-centring belongs to the *reference*, it contaminates every cell **identically** — which is
exactly a residual that does not shrink with `dt_fast`. The diel composite of the residual showed
the fingerprint before the cause was known: GPP residual at hour 09 was **−2.60 µmol m⁻² s⁻¹ at all
four `dt_fast` values**, to three digits, while hours 01–08 were exactly 0.

**Fixed** in `numerics_sweep.py`: the reference `dt_fast` must now divide **both** `dt_slow` and
3600 s. The same `--ref-refine 8` now yields 25 s.

### E-1. Corrected D-2 — all three schemes converge

`b4_stand_summer`, `--parity`, hourly FAST, reference `rk45_tight @ 25 s` (107940 sub-steps, **0
rejections**, `rk45_rescue = 0`, `clamp_mass = 0`). CAS-T RMSE [K]:

| `dt_fast` | split | ark | rk45 |
|---|---|---|---|
| 1800 | 1.195 | 1.096 | 0.902 |
| 900 | 0.692 | 0.416 | 0.398 |
| 450 | 0.522 | **0.139** | **0.136** |
| 225 | 0.457 | **0.062** | **0.061** |

soil-top T RMSE [K]: split 0.328 / 0.222 / 0.173 / 0.151; ark 0.356 / 0.120 / 0.039 / 0.016;
rk45 0.267 / 0.103 / 0.038 / 0.019.

- **ARK and RK45 converge at ~2nd order** (successive ratios 2.6 / 3.0 / 2.2) and reach 0.06 K.
- **`split` stalls at ~0.45 K.** Against its own refined self (`split @ 25 s`) it converges perfectly
  well — **1.003 / 0.396 / 0.174 / 0.079**, ratios 2.5 / 2.3 / 2.2 — so this is not a convergence
  failure. Split converges to a **different limit**, ≈ `sqrt(0.457² − 0.079²)` ≈ **0.45 K** away in
  CAS-T and ~0.15 K in soil-top T.
- **RK45 is now indistinguishable from ARK** on accuracy (0.061 vs 0.062) and cheaper (11904 vs
  13357 sub-steps at `dt = 225`; 5033 vs 5492 at 1800). D-1's "RK45 is 11× less accurate than ARK"
  was measured against `rk45_tight` at the SAME `dt_fast` and answers a different question (stepper
  error at a fixed freeze); it is not contradicted, but it is also not the whole picture.

### E-2. Efficiency — ARK dominates once accuracy is the constraint

Sub-step counts (exact, machine-independent; wall clock at these runtimes is quantisation noise):

| cell | sub-steps | rejections | CAS-T |
|---|---|---|---|
| split @ 1800 | 1488 | 0 | 1.195 |
| split @ 225 | 11904 | 0 | 0.457 |
| **ark @ 900** | **7286** | 920 | **0.416** |
| ark @ 1800 | 5492 | 1253 | 1.096 |
| rk45 @ 225 | 11904 | 0 | 0.061 |

`ark @ 900` beats `split @ 225` on **both** axes simultaneously. `split` is ~3.7× cheaper than ARK
per `dt_fast` at 1800 s for comparable accuracy — but it cannot buy accuracy at any price, and ARK
can. ARK's rejection count collapses with `dt_fast` (1253 → 920 → 171 → 3), so the rejection storm is
a coarse-`dt_fast` phenomenon and the obvious cheap optimisation target.

### E-3. Issue #79's forcing-quadrature hypothesis is REFUTED

`forcing_sample_frac` pinned to 0.0 and 1.0, each swept against its **own** refined self at 25 s so
the sample point is held common between cell and reference. CAS-T RMSE:

| frac | 1800 | 900 | 450 | 225 |
|---|---|---|---|---|
| 0.0 | 0.753 | 0.248 | 0.102 | 0.045 |
| 0.5 (default) | 1.003 | 0.396 | 0.174 | 0.079 |
| 1.0 | 1.371 | 0.611 | 0.280 | 0.130 |

Clean monotone convergence at **every** sample point, ratios 2.2–3.0 throughout. The sample point is
a real accuracy knob — a ~3× spread at coarse `dt_fast` — but it produces **no floor**. Constant
forcing (`forcing_on = false`) likewise converges cleanly on the daily flux metric (relative ET bias
3.85e-4 → 2.11e-4 → 1.02e-4 → 4.74e-5, first order); note the FAST output tier is gated on
`do_forcing`, so a constant-forcing run cannot be scored on the sub-daily stream at all — a harness
limitation worth lifting.

**⇒ Issue #79 should be closed as an instrument fault, and the "no smaller-step lever" conclusion in
D-4 is withdrawn.** The lever exists for ARK and RK45 and buys 17× over an 8× refinement.

### E-4. NEW: snowfall is discarded by ARK/RK45 when `snow_on = false`

Found by code reading, confirmed by measurement. `precip_phase` splits the met precipitation into
`rainf`/`snowf` unconditionally — the split does not consult `snow_on` — and `fill_forcing` carries
both onto `forc`. Then:

- `meds_fast_split.f90:384` — `throughfall_total = forc%precip + forc%snowf`, which becomes
  `hforc%precip_ground`. The frozen share reaches the soil as liquid. Ledger `w_in` matches.
- `meds_fast_ark.f90:1344` — `throughfall_total = forc%precip`. **`forc%snowf` is never added to
  `hforc%precip_ground`** (`:1512`), yet `w_in` counts `(precip + snowf + shed)·dt`
  (`meds_fast_ark.f90:1106`, `meds_fast_rk45.f90:655` — added by C4 for the pack case).

So with `snow_on = .false.` (the **default**) and sub-freezing air, ARK and RK45 lose `snowf·dt` of
water per step and their whole-column ledger reports it as a leak.

Measured, `b3_stand_winter`, `--parity`, `dt_fast = 1800`, month-total change in layer-mean θ:

| | `snow_on = false` | `snow_on = true` |
|---|---|---|
| split | **+0.013318** | +0.012435 |
| ark | **+0.009602** | +0.012838 |
| rk45 | **+0.009611** | +0.012750 |

Split retains **39% more** soil water than ARK/RK45 with snow off; enabling snow closes the gap to
within 3%. That is the mechanism, isolated.

**Consequence for existing results:** the 30-yr `split` run and the 29-yr `rk45` run in
`runs/ithaca_ark30/` were driven by Ithaca meteorology with `snow_on` off. Every winter of the RK45
run was missing its snowfall water. This is a candidate contributor to the ~35% AGB divergence,
alongside the recycle phase bug (PR #69).

**FIXED.** `build_column_frozen` now carries `forc%precip + forc%snowf`, as the split path always
has. Re-measured on the same month:

| | before | **after** | control (`snow_on = true`) |
|---|---|---|---|
| split | +0.013318 | +0.013318 *(byte-identical — path untouched)* | +0.012435 |
| ark | +0.009602 | **+0.013325** | +0.012838 |
| rk45 | +0.009611 | **+0.013328** | +0.012750 |

The three now agree to **0.08%** where they were 39% apart.

**Regression net: `test_column_dynamics` RUN 9.** RUN 8 could not have caught this — it runs with the
pack ON, where `snow_accumulate` consumes the snowfall before the routing branch is reached, so the
uncovered combination was the *default* one. RUN 9 runs the same day twice, snowfall on and off, at
the **same sub-freezing air temperature** (a new `snowf_rate` variable decouples the frozen-precip
flux from the cold-air switch, so an evaporation difference cannot masquerade as the water under
test), and requires the column's liquid store to gain what fell:

```
   (RUN 9 SNOWF spl column water gain=  1.7263 kg/m2  expected=  1.7280 kg/m2)
   (RUN 9 SNOWF ark column water gain=  1.6930 kg/m2  expected=  1.7280 kg/m2)
   (RUN 9 SNOWF r45 column water gain=  1.7364 kg/m2  expected=  1.7280 kg/m2)
```

Mutation-checked: with the one-line fix reverted, ARK and RK45 gain **exactly 0.000000 kg/m²** and
their whole-column water ledger fails on all 96 sub-steps. The assertion is on the *boundary input*
rather than on a residual deliberately — see E-7 for why a ledger assertion alone would not have been
enough.

### E-5 FIXED — RK45's soil moisture is now error-controlled

`state_wrms_grouped` gains an optional `with_theta`, default `.false.` (byte-identical for the ARK,
which has no theta term to add — its stage difference is structurally zero); the RK45 march passes
`.true.` and the norm uses `GRP_THETA`, already seeded from `[soil]` by `build_tol_set`.

Measured on `b4_stand_summer`, CAS-T RMSE (sub-steps in parentheses), before → after:

| scheme | 1800 | 900 | 450 | 225 |
|---|---|---|---|---|
| split | 1.1954 (1488) → **identical** | 0.6923 (2976) | 0.5220 (5952) | 0.4571 (11904) |
| ark | 1.0957 (5492) → **identical** | 0.4155 (7286) | 0.1389 (10796) | 0.0620 (13357) |
| rk45 | 0.9024 (5033) → 0.9103 (**4948**) | 0.3978 (6249) → 0.3930 (**6206**) | 0.1363 (9381) → 0.1355 (**8987**) | 0.0609 (11904) → 0.0609 (11904) |

split and ARK are byte-identical — the intended control, since neither path is touched. RK45's
accuracy is unchanged and its cost falls slightly. **Read that honestly: in this unsaturated summer
cell the temperature terms dominate the norm, so the fix is nearly free rather than beneficial. Its
value is in the regimes this cell does not exercise** — wetting fronts, saturation, and the ponding
transition, where θ moves fast and previously moved with no controller watching it. It does **not**
close D-1's ARK/RK45 gap, so that lead now belongs entirely to the `with_mass` asymmetry.

### E-6 FIXED — interception now runs after the snow stage

The interception sweep moved below `advance_snow_stage` in `build_column_frozen`, so `snow_st%exists`
is read after its writer. Nothing between the old and new positions reads `f_wet_c` or
`intercept_leaf/wood` (the first consumer is the `surface_derivs` pre-pass call further down), so the
move is otherwise inert; `canopy_water_on` defaults false, so no covered path changes.

### E-7 NEW — `[energy].debug_error` had no TOML reader, so the halt was never armed

Found while trying to verify E-4's fix the way the C3/row-1b work documented: *"run with
`energy.debug_error = .true.` so budget breaches HALT"*. **That key was never read.** `debug_error`
appears nowhere in `meds_config_io.f90`; it exists only as a `.false.` default on `energy_opts_t`,
settable from Fortran alone. So `budget_check_stop` — wired after all seven budgets by QW2 precisely
to stop it being dead code — was dead code in every *configured* run, and the runs that "completed
with zero budget messages" completed because nothing was checking.

**Wired** (`load_energy_opts`, default `.false.` ⇒ production runs unchanged). Two things fell out
immediately:

1. The earlier statement in this session's notes that "the pre-fix leak did not trip the halting
   ledger" was **wrong** — the flag had simply never been on. With it live, the pre-fix ARK does trip,
   and so does the unit fixture (96/96 sub-steps).
2. **A pre-existing, shared whole-energy closure failure on the forced winter configuration.**
   `b3_stand_winter`, `dt_fast = 1800`, all three schemes halt at
   `whole_energy (split) resid = -3.50845E+03 J/m², tol = 1.44734E+03`. Confirmed pre-existing:
   the identical residual appears with every driver change of this pass stashed, so it is not
   introduced here. It is present on `split` too, i.e. shared rather than scheme-specific. **Not
   investigated — own issue.** Note it is ~0.4% of a winter day's net radiation and does not
   destabilise the run, which is why it survived: it is exactly the size a per-step tolerance sees
   but a seasonal diagnostic does not.

**Methodological consequence worth carrying:** every "verified by a forced run with `debug_error`"
claim in this document's history is vacuous with respect to *halting* (the accompanying
positive-quantity measurements stand on their own). E-4's own regression test is deliberately built
on a **boundary-input identity** rather than on a ledger assertion for this reason.

### E-5. NEW: RK45's soil moisture is outside its own error norm

`state_wrms_grouped` (`meds_fast_control.f90:193`) contains **no θ term at all**. Its header explains
why — "theta is operator-split out of the ESDIRK stages, so its stage-difference is identically zero"
— and that is true of the ARK. It is **false of RK45**, which is the one scheme that integrates θ
inside its tableau: `state_sub` fills `y_err%theta` with a genuine, non-zero difference, and the norm
then discards it. RK45's soil-moisture error is therefore uncontrolled.

The complementary asymmetry: RK45 passes `with_mass = .true.` while the ARK march passes `.false.`,
so RK45's norm carries `2·ncoh` plant-water-mass terms the ARK's does not. With ~6 cohorts and 10
soil layers that is 12 of 25 terms — half the norm is a state the ARK does not control at all. Both
choices are defensible in isolation; together they mean the two schemes' `rtol` do not mean the same
thing, which matters for any tolerance-matched comparison.

### E-6. NEW (latent): `snow_st%exists` is read before it is set

`meds_fast_ark.f90:1345` gates the canopy-interception sweep on `.not. snow_st%exists`, but
`advance_snow_stage` — the only writer — does not run until `:1378`. `snow_stage_t` has default
component initialisation (`exists = .false.`), so the read is defined rather than undefined, but it
is always `.false.`: **ARK/RK45 will intercept rain in the canopy under a snow pack where `split`
skips the sweep** (`meds_fast_split.f90:394`). Only bites with `canopy_water_on = .true.` (default
off) and a pack present, so it is latent — but it is the same class of ordering hazard C4 documented
and should be closed while it is cheap.

## 3f. Phase F — two options removed (2026-07-29)

Both were switches over a *structural* fact, and in both cases the default silently did the wrong
thing. Removing them is a net simplification and, in one case, a measurable speed-up.

### F-1. `[fast].snow_on` is gone — the snow store is always active

Snowfall is a boundary water input like rain, and `precip_phase` splits it out of the met
precipitation **without consulting any switch**. So `snow_on = .false.` never meant "this run has no
snow"; it meant "this run receives snow and has nowhere to put it" — which is exactly how E-4
happened, and why its default value was the dangerous one.

Always-on is free on a snow-free column: `snow_accumulate` returns immediately unless a pack exists
or the snowfall clears `params%min_new_snow_mass`, so `snow_stage_t` keeps its bare-ground defaults,
`snowfac = 0`, and `surface_derivs`' snow blend reduces exactly to its pre-C4 form. Sub-threshold
snowfall onto bare ground still reaches the soil as liquid through the caller's throughfall routing,
so nothing is dropped in that branch either.

Removed from: `meds_config_t`, `column_config_t`, `load_meds_config`, the `fast_context_t` copy, the
RT albedo guard, `meds_main`'s pack seed, `advance_snow_stage`'s early return, and the test suite.
**Verified:** RUNS 1–7 of `test_column_dynamics` are unchanged to every printed digit (CAS noon
292.66 K, GPP 17.90 µmol m⁻² s⁻¹, whole-column residuals 4.8e-7 J/m² and 2.4e-13 kg/m²), and the
summer forced sweep's `split` column is byte-identical. The winter month now reproduces the old
`snow_on = true` control exactly:

| | split | ark | rk45 | spread |
|---|---|---|---|---|
| `snow_on = false`, pre-E-4 | +0.013318 | +0.009602 | +0.009611 | 27.9% |
| `snow_on = false`, post-E-4 | +0.013318 | +0.013325 | +0.013328 | 0.1% |
| **always on (final)** | **+0.012435** | **+0.012837** | **+0.012738** | **3.2%** |

The residual 3.2% is genuine physics — the three schemes drive the surface balance along different
CAS trajectories, so they melt and sublimate slightly differently — and it is the same ≤3% RUN 8
already asserts on the melt total.

### F-2. `with_mass` / `with_theta` are gone — one norm over the whole column state

`state_wrms_grouped` had two optional switches selecting which states entered the error norm. They
encoded a real structural fact (on the ARK path both theta and plant water mass are operator-split
out of the tableau, so their embedded differences are *exactly* zero and summing them only inflates
`cnt`), but they were the wrong mechanism for it: at a call site, "do not measure this state" and
"this state cannot move" are indistinguishable, and the two **were** confused — RK45 integrates theta
inside its tableau and the norm discarded it anyway (E-5).

Now every field of `column_state_t` contributes, unconditionally. The signature loses two dummies and
the body loses two branches.

**The dilution this accepts is real but not measurable here.** ARK, `b4_stand_summer`, CAS-T RMSE
(sub-steps in parentheses), before → after:

| dt | before | after |
|---|---|---|
| 1800 | 1.0957 (5492) | 1.0956 (**4678**) |
| 900 | 0.4155 (7286) | 0.4168 (**6467**) |
| 450 | 0.1389 (10796) | 0.1396 (**10425**) |
| 225 | 0.0620 (13357) | 0.0621 (**12192**) |

Accuracy moves by ≤0.5%; cost falls **9–15%**. The reason is §3d/E-1: at production `dt_fast` the
error is dominated by the Category-0 freeze, not the stepper, so loosening the stepper slightly does
not move the total. `split` and `rk45` are byte-identical (split has no adaptive march; RK45 already
included both terms after E-5).

**Worth stating so it is not mistaken for a regression in control:** soil moisture is
error-controlled on *every* path regardless of this norm. `split` and the ARK take theta wholly from
`column_hydrology_flux`, whose own adaptive step-doubling is driven by the **same** `GRP_THETA`
tolerances `build_tol_set` seeds here. What the norm adds is control for the one scheme that took
theta out of that solver and into its own stages. The clean way to retire the remaining dilution is
to bring soil water into the ARK tableau — the P1 item, unchanged.

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
- **Phase C — parity fixes.** **C1, C3, C4, C5, row 12, row 1b DONE; C2 PARTIAL** (see §3c).
  Open: C2's ponding/runoff half (#75), the condensate destination (#74); C6 documented and deferred.
- **Phase D — the numerical comparison.** DONE; see §3d. Headline: at fixed `dt_fast` the stepper
  choice is worth ~21× in accuracy (ARK over split) for 1.3× cost. ~~while refining `dt_fast` does
  not converge at all~~ — **D-2/D-3 withdrawn, see Phase E.**
- **Phase E — the reference audit.** DONE 2026-07-29; see §3e. Headline: the D-2 non-convergence was
  an instrument fault (a reference `dt_fast` that does not divide one hour). With it fixed, ARK and
  RK45 converge at ~2nd order to 0.06 K while `split` converges to a limit ~0.45 K away — so both
  levers exist, and the remaining question is the residual family gap, not the step size. Also
  found: snowfall discarded by ARK/RK45 under `snow_on = false` (E-4, live bug); RK45's soil moisture
  outside its own error norm (E-5); a latent ordering hazard in `build_column_frozen` (E-6).

## 5. Note on existing results

The 29/30-yr runs in `runs/ithaca_ark30/out/` predate PR #69 (`05c0a8d`), which fixed a forcing
recycle bug that phase-scrambled sub-daily shortwave in multi-decade recycled runs. Daily-mean
shortwave was correct throughout, so slow diagnostics looked healthy while the diel cycle was
wrong. Any conclusion drawn from those outputs about sub-daily behaviour — including the ~35% AGB
divergence — needs re-deriving on post-fix runs.
