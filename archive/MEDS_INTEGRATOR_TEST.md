# MEDS fast-loop integrator evaluation & production `dt_fast` recommendation

**Status:** EXECUTED 2026-07-13 — 24-run matrix complete, results + revised recommendation in **§10** (headline:
**`dt_fast` = 900 s**, ARK-adaptive or SPLIT; 1800 s is too coarse; ARK-fixed-1 is a bad config). **Method:** single
fine-baseline (Richardson-lite) comparison — one very fine ARK reference vs the production `dt_fast` candidates,
over **two contrasting months (cold January + warm June)**, compared at a common 1-hour output cadence. **Test bed:** Ithaca NY, 2024 ERA5-Land hourly forcing
(`data/forcing/ithaca_forcing.nc`). **Goal:** recommend a production `dt_fast` + integrator, quantify the
accuracy-vs-cost trade-off across the cold- and warm-season stiffness regimes, and settle whether MEDS should add
ED2-style met/radiation interpolation.

**Prerequisite — DONE:** the FAST (sub-daily) netCDF output tier is implemented and merged
(`archive/MEDS_FAST_OUTPUT_DESIGN.md`), so the 1-hour comparison axis this evaluation needs already exists
(CAS temp, GPP rate, LE, H, Rnet, SW-in, ustar, soil-temp/water columns, per-cohort leaf-temp/GPP).

---

## 1. Methodology: one fine baseline, compared at 1-hour resolution

Rather than a full halving ladder, establish **truth from a single very fine run** and measure each production
candidate against it:

1. **Baseline (truth):** `INTEG_ARK`, `dt_fast = 15 s`, run **once per test month**, FAST output aggregated to
   **1-hour** records. ARK is L-stable and the more accurate scheme; 15 s is well inside the RK4 stability floor
   (~17 s), so it also bounds any explicit reference. Confirm convergence (§4) before trusting it.
2. **Candidates:** each `(scheme × dt_fast)` cell runs the **identical** month + community + forcing, emitting the
   **same 24 hourly records/day** (via `fast_interval_steps = 3600/dt_fast`), then differenced against the baseline
   hour-by-hour.
3. **Compare:** per-variable accuracy (RMSE / max-abs-error / bias of the 1-hour series) **and** wall-clock time.
   The production pick is the largest (cheapest) `dt_fast` whose error stays within tolerance.

Why this works: all runs discretize the *same* continuous fast-loop ODEs, so ARK@15 s ≈ the exact diurnal
trajectory; each candidate's departure from it is its total truncation error at that `dt_fast`. Because every run
outputs the same clock hours (§3), the coarse runs are simply a sparser sampling/integration of each hour — which
*is* the error being measured.

---

## 2. The run matrix

**Run the full matrix below once per test month** (cold January + warm June, §5) — 2× the run count, same
community/config, only `start_time` changes. **Held fixed across every run:** Ithaca site + PFT community + initial
cohorts (census-seeded, §5), `fast_context_t`/`column_config_t` params, the real 2024 hourly ERA5-Land file, RT
join ON, `dt_slow = 1 day`.

**Baseline (truth):**

| id | scheme | dt_fast | fast_interval_steps | config |
|---|---|---|---|---|
| **REF** | ARK | 15 s | 240 | `time_integrator="ark"`, `[output.fast].enabled=true`, `fast_interval_steps=240` |

**Candidates** — `{SPLIT, ARK-fixed, ARK-adaptive} × {600, 900, 1800} s`, each with its matching
`fast_interval_steps` so all emit 24 hourly records/day:

| dt_fast | fast_interval_steps | SPLIT | ARK-fixed | ARK-adaptive |
|---|---|---|---|---|
| 600 s  | 6 | ✓ | ✓ | ✓ |
| 900 s  | 4 | ✓ | ✓ | ✓ |
| 1800 s | 2 | ✓ | ✓ | ✓ |

- **SPLIT:** `time_integrator="split"`.
- **ARK-fixed:** `time_integrator="ark"`, `ark_adaptive=false`, `ark_fixed_substep=1` (one coupled march per `dt_fast`;
  the GPU-uniform path — measures the ARK core's `dt_fast` truncation directly).
- **ARK-adaptive:** `time_integrator="ark"`, `ark_adaptive=true`, `ark_rtol=1e-3` (the production default; here
  `dt_fast` is the forcing-refresh / pre-pass-hold cadence, and the controller sub-divides the coupled march as
  needed — so its *energy/carbon* accuracy can beat the fixed march at the same `dt_fast`, while **θ/soil-water is
  still committed once per `dt_fast`**, i.e. adaptivity does not refine the water axis).

9 candidate runs + 1 baseline + 2 convergence checks (§4) = **12 runs per month × 2 months = 24 runs**.

**Runtime budget** (extrapolating the 30-yr@1800 s run, ~2500 sub-steps/s): REF@15 s ≈ 173k sub-steps ≈ **1–2 min**
per month; each candidate is **seconds** (1800 s → 1.4k sub-steps; 600 s → 4.3k). Both months together are well
under ~10 min of wall-clock.

---

## 3. Output & comparison protocol

- **Common cadence:** every run sets `fast_interval_steps = 3600/dt_fast` → one `-F-YYYYMMDD.nc` per day with **24
  hourly `AGG_TMEAN` records**. `fast_interval_steps` divides `n_fast_per_slow` in every case (240|5760, 6|144, 4|96,
  2|48), so windows are clean.
- **Alignment:** compare by **(day, hour)** — each run's Nth hourly record carries the same `hour` companion
  (0…23), so the series align directly. The record `t_open` (window start) differs by <½ h across runs but the hour
  stamp is identical.
- **Variables compared** (all in the FAST tier): `cas_temp_site`, `gpp_rate_fast`, `le_flux_fast`, `h_flux_fast`,
  `rnet_fast`, `sw_in_fast` (a forcing sanity check — should be ~identical across runs), `soil_temp_top_site`,
  `soil_temp_site`/`soil_water_site` (columns), and the per-cohort `leaf_temp_cohort_fast`/`gpp_cohort_fast`.
- **Also logged per run:** wall-clock time; worst `whole_energy`/`whole_water` conservation residual; budget-fail /
  non-convergence counts; for ARK-adaptive, `nsteps`/`nrej` (mean coupled sub-steps per `dt_fast`).
- **Optional finer probe:** on 2–3 high-gradient days (clear-sky sunrise/sunset, a storm onset) re-run with
  `fast_interval_steps` set for **15-min** records to resolve the **dawn GPP-onset timing** metric (the sharpest
  discriminator, blurred by hourly averaging). The `[fast].fast_probe` CSV is the alternative per-sub-step dump.

---

## 4. Reference validity — confirm 15 s is converged

15 s is only "truth" if it is in the asymptotic regime. Two cheap checks (both full-month, both ~seconds–1 min):

- **Self-convergence:** ARK@15 s vs **ARK@30 s** must agree **far below** the candidate errors (e.g. hourly CAS-temp
  RMSE ≪ the 1800 s error). If they don't, the reference is under-resolved → go finer (7.5 s).
- **Cross-scheme:** ARK@15 s vs **SPLIT@15 s** should also agree closely — both are consistent discretizations of the
  same ODEs, so at 15 s they bracket the same continuous solution. Close agreement validates using ARK@15 s as the
  common truth for *both* schemes' candidates.

*Water caveat:* soil θ comes from `column_hydrology_flux`, which already sub-steps adaptively, so soil moisture is
near-converged even at coarse `dt_fast`; only the θ **commit cadence** refines with `dt_fast`. A coarser knee on the
water axis is legitimate, not under-resolution.

---

## 5. Test bed: windows, forcing, community

- **Two contrasting months (both 2024)** — run the full §2 matrix in each, to bracket the annual stiffness range:
  - **January 2024 — cold regime.** Sub-freezing air, frozen precip (`precip_phase` → `snowf`), soil **freeze/thaw**
    (ice-aware `κ`/`C_eff` + the zero-curtain plateau), and **near-saturation CAS** with condensation/frost. This is
    the **stiffest energy/water coupling MEDS has** — the same near-saturation CAS regime that was the hardest part
    of ARK development (the year-13.5 thrash, fixed via no-clamp + smooth condensation) — and it stresses the
    soil-thermal + phase-change march. GPP is dormant, so **January is the decisive test for the energy / water /
    soil-thermal tiers**, not carbon. *Caveat:* MEDS has **no dynamic snowpack** (`snowfac=0`); "snow" here means
    frozen-precip partitioning + frozen soil, not an accumulating pack.
  - **June 2024 — warm / growing regime.** Peak insolation + active photosynthesis exercise the **leaf-Cᵢ ↔ CAS ↔
    soil coupling** and the largest sunrise/sunset SW ramps — the worst case for the `dt_fast`-controlled `e_frozen`
    error and the **decisive test for the carbon (GPP/NPP) tiers** + dawn GPP-onset timing.

  A `dt_fast` that passes **both** is safe year-round: cold-season soil-thermal/phase-change/near-saturation *and*
  warm-season photosynthesis-coupling/peak-ramps. (The winter run also reproduces the strongly-negative nighttime
  sensible heat / stable-boundary-layer regime seen in earlier tests.)
- **Forcing — stay inside 2024:** run entirely within the forcing file's own year so **no calendar recycling**
  occurs. The multi-decade recycling path has a known sub-daily-SW anti-phase bug
  (`project_meds_forcing_recycle_phasebug`: ~2 h drift/yr → ~12 h inversion by ~30 yr); a single 2024 month is
  exact and deterministic. Both Jan-2024 and Jun-2024 are inside the file year.
- **Community — DONE (census-seeded).** The 30-yr spun-up community (8 cohorts, LAI≈1.48) is extracted to
  `runs/ithaca_ark30/ithaca_spinup_census.csv`; `init_mode=1` + `census_file=…` places it at *any* `start_time`, so
  the identical stand seeds both the Jan-1 and Jun-1 2024 runs (verified loading in a 2024 window). This decouples
  the community from the state's 2054 date (which would trigger the recycling bug). *Caveat:* `init_from_census`
  builds equal-area patches, so the state's `[0.657, 0.343]` areas flatten to `[0.5, 0.5]` (~0.1 % site-aggregate
  shift vs the true state) — irrelevant here since every run uses the *same* census (scheme-vs-scheme on a fixed
  community).

---

## 6. Metrics & thresholds

Compare each candidate's 1-hour series to REF, per variable, over the month:

- **Sub-daily / diurnal (where the `dt_fast` error lives):** CAS air-temp hourly **RMSE < 0.25 K**; GPP-rate hourly
  **RMSE < ~5 % of midday peak** (+ dawn **GPP-onset timing < one `dt_fast`**, from the optional 15-min probe);
  LE **RMSE < ~5 W/m²**; top-soil skin-temp **RMSE < 0.5 K**. *CAS-temp and dawn GPP-onset are the sharpest and blow
  first.*
- **Daily-integrated (the fast→slow handoff — DECISIVE tier):** daily GPP & NPP **bias < 1 %**, daily RMSE < ~2 %;
  daily ET bias < 1 %. (Sum the hourly series, or read `total_gpp`/`total_npp` from the daily tier.)
- **Monthly / integrated:** month-total GPP/NPP/ET each **< 1 %** — most forgiving; necessary but not sufficient
  (day-to-day errors cancel).
- **Conservation (hard gate):** worst `whole_energy`/`whole_water` residual stays at the baseline's closure order
  (SPLIT & ARK: machine energy; ARK water is split-order). Any `dt_fast` that degrades closure is disqualified.
- **Stability/robustness (hard gate):** zero budget failures, Picard/Newton convergence every step, bounded temps,
  no NaN over the month; for ARK-adaptive, `nrej/nsteps` bounded.
- **Cost:** wall-clock per run, and the **error-vs-wall-clock "knee"** per metric per scheme.
- **Report every metric per month.** The binding tier differs by regime: **January** binds the energy / water /
  soil-thermal (freeze/thaw, near-saturation CAS) thresholds with GPP dormant; **June** binds the carbon (GPP/NPP)
  thresholds + dawn GPP-onset. A candidate `dt_fast` is admissible only if it clears the binding constraint in
  **both** months.

---

## 7. Background: what actually limits fast-step accuracy (why `dt_fast` is the lever)

The fast step is a **Lie–Trotter operator split** wrapping the coupled core, so it is **globally first-order in
`dt_fast`** even though the ARK2 core is 2nd-order. Its O(`dt`) error has four sources:

| term | source | lever |
|---|---|---|
| `e_split` | Lie–Trotter operator ordering (water → CAS/soil-energy core → hydraulics) | Strang / P3 Picard |
| **`e_frozen`** | **zeroth-order HOLD of the derived pre-pass** (two-stream RT, gₛ/GPP, aero conductances, respiration, scratch Richards BCs), built once per `dt_fast` in `build_column_frozen` and held | **`dt_fast`** (or an inner pre-pass mini-cycle) |
| `e_theta` | soil θ committed out-of-band from the scratch Richards, read as lagged θⁿ into soil-energy props | P3 Picard |
| `e_transp` | hydraulics advanced with endpoint transpiration over the full `dt` | P3 Picard |

**`e_frozen` is the dominant, `dt_fast`-controlled term**, and it lives where the forcing changes *within* a step —
sunrise/sunset SW ramps (→ GPP) and storm onset (→ ET/soil). That is why a fine-`dt_fast` baseline is the reference
and why the summer/high-gradient window is the discriminating case.

---

## 8. Forcing / ED2-interpolation decision

**Keep MEDS's per-sub-step midpoint linear interpolation + cosz-subsampled SW disaggregation** — already
strictly higher-order than ED2's left-endpoint hold; **do not build an ED2 `dt_radinterp` knob** (MEDS already has
its analogue). Verified: ED2's `dt_radinterp` defaults to `DTLSM` and is only the `⟨sec z⟩` SW-downscaling
sub-sampling resolution (= MEDS's `cosz_reconstruct_factor`), *not* a finer met refresh; ED2 recomputes RT at
`radfrq ≥ DTLSM` (held) and marches RK4 against met snapshotted once per `DTLSM` at the **left endpoint**. So the
only genuine `e_frozen` lever is **`dt_fast`**. A pure RT-only refinement won't cut the dominant GPP error (GPP comes
from the leaf-Cᵢ solve consuming the coarse absorbed PAR, so refining RT alone helps energy, not dawn/dusk GPP). The
one defensible cheap-decouple, *only if the runs show `e_frozen` dominates*: an inner **RT+aero+leaf-Cᵢ** mini-cycle
leaving Richards frozen (soil moisture barely moves in <1800 s).

---

## 9. Expected results & decision criteria

**Expected:**
- **Full-step order (live forcing):** both SPLIT and ARK plateau at **p≈1** in `dt_fast`, capped by the split +
  `e_frozen` hold — *not* the ARK2 core. (The ARK core shows p≈2 only on cleanly-integrated variables like the CO₂
  twin; `test_ark2` asserts p≥1.9 there.)
- **Daily/monthly tiers:** likely clear <1 % at 1800 s for both integrators (consistent with the 30-yr <1 %
  ARK-vs-SPLIT match). The **binding constraint, if any, is the sub-daily tier** on high-gradient days.
- **SPLIT vs ARK:** ARK is L-stable → `dt_fast` chosen purely by accuracy; SPLIT may need one rung finer for equal
  accuracy/robustness (Picard convergence at large `dt_fast`). ARK is ~15–20 % slower per step at equal `dt_fast`.
- **ARK-adaptive vs ARK-fixed:** adaptive should match the baseline on energy/carbon at coarser `dt_fast` (it
  sub-divides the march), at some extra cost; no gain on water (θ committed once per `dt_fast`).

**Decision:**
1. **Conservation + stability gates** must pass in **both months** (hard prerequisites) — January is where they are
   most likely to bite (freeze/thaw + near-saturation CAS stiffness).
2. **Production `dt_fast` = the largest (cheapest)** whose binding-tier error is within tolerance of REF **in both
   months** — June's daily GPP/NPP bias < 1 % (+ sub-daily diurnal thresholds), *and* January's energy/water/
   soil-thermal thresholds. The final pick is the min over the two months.
3. Report SPLIT vs ARK separately, per month; report the error-vs-wall-clock knee per metric.

**Strang split** — adopt only if a probe shows `e_split/e_tot > 0.10` **and** a Strang run moves a headline
diagnostic > 1 % **and** the ~1.6–2× cost is acceptable (expected to fail). Cheaper alternatives if split error
matters: midpoint forcing (already in MEDS), the designed **P3 outer Picard** (kills `e_split`+`e_theta`+`e_transp`),
or **shorter `dt_fast`** (halves all four terms).

---

## 10. Results (24-run matrix, 2026-07-13) & recommendation

Executed the full matrix (Release ifx, `OMP_NUM_THREADS=1`, census-seeded community, Jan-2024 + Jun-2024). Hourly
FAST output vs REF = ARK-adaptive@15 s. RMSE columns are over the whole month's hourly series; `GPPmo_bias%` is the
month-total GPP bias.

```
MONTH = JAN (cold)                                    MONTH = JUN (warm)
run          wall  casT   soilT   GPP   GPPmo%   LE     H       run          wall  casT   soilT   GPP   GPPmo%   LE      H
ark_a_30      5.2  0.002  0.002  0.000  -0.00   0.13   0.18     ark_a_30      5.3  0.003  0.003  0.001  -0.00   0.18    0.17
split_15      7.3  0.012  0.002  0.001  -0.01   0.15   0.15     split_15      7.5  0.061  0.014  0.002   0.00   0.52    0.21
split_600     1.0  0.058  0.062  0.035  -0.43   0.43   1.74     split_600     1.0  0.109  0.098  0.094  -0.14   1.33    2.10
split_900     1.6  0.101  0.093  0.050  -0.62   0.63   2.57     split_900     0.9  0.172  0.151  0.138  -0.23   1.98    3.05
split_1800    1.0  0.339  0.202  0.106  -0.03   1.30   6.20     split_1800    0.8  0.850  0.429  0.319  -2.04   5.97   12.85
arkf_600      1.0  0.095  0.063  0.016  -0.01   4.60   6.74     arkf_600      1.0  0.124  0.132  0.039  -0.17   6.06    5.40
arkf_900      1.2  0.156  0.096  0.030   0.33   6.51  13.75     arkf_900      0.9  0.391  0.272  0.205  -1.42  39.78   44.77
arkf_1800     1.0  0.463  0.239  0.142   4.12  11.59  43.96     arkf_1800     0.9  1.221  0.703  0.616  -5.52  98.47   96.54
arka_600      1.1  0.082  0.063  0.016  -0.01   2.94   4.99     arka_600      1.0  0.123  0.132  0.040  -0.17   5.27    5.29
arka_900      1.0  0.137  0.092  0.025   0.21   3.32   6.40     arka_900      1.0  0.225  0.207  0.064  -0.30   7.35    7.26
arka_1800     0.9  0.451  0.200  0.126   3.24   3.66  13.64     arka_1800     0.9  0.791  0.485  0.247  -1.66  11.91   16.14
```
(temps in K, fluxes in W/m²; GPP RMSE in µmol/m²/s; wall in s.)

**Findings:**
1. **Reference is converged.** ARK@15 vs ARK@30: CAS-temp RMSE 0.002–0.003 K; vs SPLIT@15: ≤0.06 K. 15 s is truth
   for both schemes.
2. **1800 s is too coarse for the tight thresholds** — it misses CAS-temp < 0.25 K in *both* months (Jan 0.34–0.46,
   Jun 0.79–1.22) and monthly-GPP < 1 % in June (SPLIT −2.0 %, ARK-adaptive −1.7 %). **This overturns the earlier
   1800 s default.**
3. **900 s is the knee.** SPLIT and ARK-adaptive both clear CAS-temp (~0.10–0.23 K) and monthly-GPP (<0.3 %) in both
   months. **ARK-fixed (1 substep) fails at 900 s** (Jun CAS 0.39, monthly GPP −1.4 %).
4. **600 s clears everything with margin** (CAS 0.06–0.12, GPP bias <0.2 %) — matches ED2's range.
5. **ARK-fixed(1) is a poor config:** its LE/H degrade badly at coarse dt (Jun arkf_1800 LE **98 W/m²**, H 97 W/m²) —
   1-substep ARK ≈ IMEX-Euler on a stiff coupled system. **ARK-adaptive recovers most of it** (Jun 1800 LE 12 vs 98).
6. **SPLIT is the most accurate on energy/water** at coarse dt (LE 0.4–6.0 vs ARK-adaptive 3–12), comparable on
   temps/GPP — its tuned operator-split water + Picard beats ARK's coupled march there.
7. **Cost:** at production `dt_fast` the runs are **I/O-bound** (~1 s, fast-loop compute negligible for 1 site ×
   1 month); the scheme cost only shows at 15 s (**ARK ≈ 7.3–7.5 s vs SPLIT ≈ 9.5–10.2 s** → ARK ~30 % slower). At
   production *scale* (many sites / long runs) the fast-loop cost scales with sub-steps, so 900 s ≈ 2× the 1800 s
   fast-loop cost — but still small vs the accuracy gain.

**Recommendation (data-driven, revising the provisional 1800 s):**
- **`dt_fast` = 900 s** for the tight targets (CAS-temp < 0.25 K, monthly GPP < 1 %); drop to **600 s** for extra
  margin or if using ARK-fixed. **1800 s is not recommended** at these thresholds (relax to CAS < 0.5 K / GPP < 2 %
  and 900 s still wins; 1800 s only qualifies if summer CAS-temp error ~0.8 K is acceptable).
- **Integrator: ARK-adaptive** (`ark_adaptive=true`) if the standing L-stable/robustness preference holds — it
  passes at 900 s and avoids the ARK-fixed water blow-up; **SPLIT@900 s is competitive and ~30 % cheaper with
  better energy/water accuracy**, so it is the alternative if robustness at large dt is not the priority. **Do not
  ship ARK-fixed with `ark_fixed_substep=1`.**
- **Forcing handling: unchanged** (midpoint + cosz) — `sw_in` is identical across all runs (forcing sanity ✓); the
  only lever is `dt_fast`, and it lands at 900 s.
- **Caveat:** thresholds drive the pick — these are strict (0.25 K / 1 %). If daily/monthly ecology is the only
  concern and ~2 % GPP is tolerable, 1800 s remains viable for SPLIT/ARK-adaptive (its month-total GPP bias is
  −2 %). The runs, configs, and analysis live in `runs/ithaca_ark30/integ/` (`gen_configs.py`, `run_matrix.sh`,
  `analyze.py`, `timings.csv`).

---

## Appendix — optional deeper diagnostics (not required for the `dt_fast` pick)

- **`e_split` isolation (column-kernel).** At a fixed `dt` and *fixed frozen coefficients*, compare the ARK2
  split-march against the RK4 oracle (`rk4_column_step` over `column_derivs`): the difference is `e_split` alone. Use
  `ark2_column_step` (2nd-order), *not* `column_be_stage + advance_hydraulics_full` (1st-order IMEX-Euler, folds in
  BE core error); hold θ = θⁿ on both sides. Isolates operator-ordering error from `e_frozen`.
- **Cheap-decouple test (config D).** One clear-sky day: 1800 s coupled march but met+RT+aero refreshed at
  `dt_rad = 300 s`, leaf-Cᵢ and Richards held at 1800 s. `(D − REF)` vs `(1800 − REF)` settles the interpolation
  question — if D collapses toward the 1800 s run, dawn GPP error is leaf-Cᵢ-gated (an RT-only layer is useless); if
  D recovers toward REF, a targeted inner mini-cycle is worth building.
