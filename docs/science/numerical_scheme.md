# Numerical schemes of the fast loop

How MEDS actually advances the sub-daily biophysics in time, what the three available integrators
assume, where they differ, and which one to use. Written for ecological modellers rather than
numerical analysts: the goal is that you can read a MEDS run's `[fast]` block and know what it
committed you to.

Companion pages: [column_biophysics](column_biophysics.md) describes *what* is being integrated
(the stores and their fluxes). This page is about *how*. The engineering record — measurements,
per-difference status, open work — lives in
`docs/dev_plans/MEDS_INTEGRATOR_PARITY.md` and `docs/dev_plans/MEDS_NUMERICS_SCOPING.md`.

---

## 1. The problem: one column, wildly different clocks

Each `dt_fast` MEDS advances one patch's whole column: canopy air (heat, moisture, CO₂), leaf and
wood temperature, leaf and wood internal water, canopy surface water, a snow pack, a soil thermal
column and a soil water column. These are coupled — the leaf's temperature sets its transpiration,
which sets canopy humidity, which sets the leaf's temperature — so they cannot be advanced
independently without penalty.

The difficulty is not the coupling itself but that the pieces have **response times spanning five
orders of magnitude**:

| store | how fast it responds | why |
|---|---|---|
| leaf temperature | seconds | almost no heat capacity, strong radiative + convective coupling |
| plant water potential | ~20 s | small tissue capacitance, high xylem conductance |
| canopy air | ~2 min (daytime) | thin air layer, vented by turbulence |
| snow surface | minutes | thin, low heat capacity |
| soil surface layer | hours | real heat capacity, conduction-limited |
| deep soil, soil moisture | days to weeks | large stores, slow transport |

A system with that spread is called **stiff**. The practical consequence is blunt: an ordinary
explicit time step (predict the rate, multiply by `dt`, add) is *unstable* unless `dt` is smaller than
the fastest response time. Advancing a whole column at the fastest mode's pace — tens of seconds —
for thirty simulated years is not affordable. Every land model therefore does something cleverer, and
what it does is the single largest determinant of its sub-daily behaviour.

---

## 2. The one decision that matters most: what gets frozen

Before any integrator runs, MEDS makes a **semi-discretisation** choice: which quantities are
recomputed continuously as the column evolves, and which are computed once at the start of the step
and held fixed. This is inherited from ED2 and it is deliberate.

MEDS sorts the column into three categories.

**Category 0 — frozen once per `dt_fast`** (`column_prepass`, `build_column_frozen`):

- leaf photosynthesis, stomatal conductance and dark respiration (GPP, $g_s$, $R_d$)
- canopy radiative transfer (absorbed shortwave and net longwave per cohort)
- stem and fine-root maintenance respiration
- aerodynamic conductances
- plant hydraulics — solved *exactly* (a matrix exponential over the step) and its **time-averaged**
  sapflow and root uptake handed to the rest of the column as constants
- infiltration, and (for the ARK) the whole soil-water solve

This is ED2-faithful. ED2 calls `canopy_photosynthesis` once before its RK integrator and the RK
stages only *read* the frozen conductances; `plant_hydro.f90` freezes hydraulic fluxes per `DTLSM`
for exactly the same reason. Freezing hydraulics removes the ~20 s mode from the coupled system,
which is what makes any of the schemes below affordable.

**Category 1 — stiff, must be treated implicitly**: canopy-air enthalpy / humidity / CO₂, the
diagnostic leaf and wood energy balance, the ground skin, the soil thermal column, the soil water
column.

**Category 2 — non-stiff**: leaf and wood internal water mass, canopy surface (film) water,
canopy-air CO₂, the transpiration and aerodynamic fluxes themselves.

**The freeze is not free.** It introduces an error that no time integrator can remove, because it
changes the equations being solved, not their discretisation. At production `dt_fast = 1800 s` this
freeze error is comparable to — often larger than — the time-stepping error. Two things follow, and
they are the most useful facts on this page:

1. Refining `dt_fast` refines **both** the freeze and the time step, so a `dt_fast` sweep measures
   them together.
2. Making the integrator more accurate at a *fixed* `dt_fast` eventually stops helping, because the
   freeze error is the floor.

---

## 3. The three integrators

Selected by `[fast].time_integrator`. All three share the same physical kernels, the same frozen
pre-pass, and the same conservation ledgers.

### `split` (default) — operator splitting with backward Euler

Advance the stores **one after another** in a fixed order, each with a stable implicit (backward
Euler) solve, each seeing its neighbours at whatever time level they happen to be at. One pass per
`dt_fast`, no sub-stepping, no error estimate.

*Analogy:* a relay race. Each runner covers the whole distance in one go and hands over; nobody runs
alongside anybody.

- **Strength:** cheap and unconditionally stable. Exactly one step per `dt_fast` regardless of how
  stiff the column gets — it never thrashes, never rejects, never bails.
- **Weakness:** the splitting itself is only first-order accurate, and there is **no way to buy
  accuracy other than shrinking `dt_fast`**, which also shrinks the slow loop's coupling interval and
  costs linearly.
- **Optional refinement:** `[fast].integration_scheme = "picard"` iterates the
  {leaf → soil → canopy air} block to a fixed point instead of a single pass. Required for a
  prognostic leaf energy store (a single explicit pass with leaf heat capacity oscillates at 2·dt).

### `ark` — implicit–explicit additive Runge–Kutta (ARS(2,2,2))

Advance the coupled stiff block with a **second-order, L-stable, two-stage implicit** method, with an
**embedded error estimate** that drives adaptive sub-stepping inside each `dt_fast`.

*Analogy:* a group ride with a pace-setter. Everyone moves together, and the pace is adjusted from
how far the group drifted apart on the last stretch.

- Canopy air + soil heat are inside the tableau. Soil water, plant water mass and canopy surface
  water are **operator-split out** and advanced once over the full `dt_fast`.
- Adaptive: the step size is chosen from the difference between the second-order and first-order
  solutions, measured by **one weighted norm over the whole column state** with a per-variable-group
  tolerance (`[fast].ark_rtol`, `rtol_all`, `atol_scale`). The ARK's soil-moisture and plant-water
  terms in that norm are structurally zero (both ride operator-split maps outside the tableau), so
  they dilute it slightly — accepted deliberately, since one norm with no per-scheme opt-outs is
  what stops an *integrated* state going unmeasured, as soil moisture once did on the RK45 path.
- Warm-started: the last accepted sub-step size carries over to the next `dt_fast`, which cut the
  rejection rate from ~31% to under 20%.

### `rk45` — adaptive explicit Cash–Karp 5(4), ED2-style

A fully explicit fifth-order Runge–Kutta with an embedded fourth-order estimate, marched adaptively
over the *whole* column state — this is the closest analogue of ED2's own `rk4` integrator.

*Analogy:* the same group ride, but nobody is allowed to look ahead; if the road gets steep the only
option is to take smaller steps.

- Everything — including soil moisture — is inside the tableau. It is the only scheme that integrates
  its own soil water rather than taking it from the pre-pass solve.
- **It is a hybrid, not a pure explicit scheme.** When the explicit march cannot resolve a step (a
  cold night under a closed canopy, where the leaf↔canopy-air coupling is stiff and the state rails
  against its physical bounds), the driver **discards the RK45 step and redoes that `dt_fast` on the
  `split` path**. The counter `work_rk45_rescue_site` reports how often. *Read it before reading any
  RK45 result* — a run whose winter months all rescue is substantially a `split` run wearing an RK45
  label.

---

## 4. What differs between them, and what kind of difference it is

Three categories, because they demand different responses. A *physics* difference means the schemes
solve different equations and no refinement will reconcile them; a *numerics* difference shrinks as
you refine; a *bookkeeping* difference is a conservation defect.

| difference | split | ark | rk45 | kind |
|---|---|---|---|---|
| canopy-air condensation sink | exact exponential relaxation, once per step | per-stage rate | per-stage rate | numerics (same model, two quadratures) |
| snow | shared stage, always active | same | same | unified |
| canopy surface water | present | present | present | unified |
| sapflow advected enthalpy | present | present | present | unified |
| soil water | implicit Richards, re-solved per Picard pass | frozen pre-pass solve, outside the tableau | integrated in the tableau | numerics + assumptions |
| adaptive error norm | n/a | one norm over the whole column state | same | unified |
| per-layer root-sink placement (`multilayer_roots`) | tracks realized uptake | static root profile | static root profile | physics (defaults off) |
| prognostic leaf/wood energy | supported | hard error stop | hard error stop | physics (non-silent) |
| non-free-drain bottom boundary | supported | hard error stop | hard error stop | physics (non-silent) |
| stability clamps | none needed | extrapolation base only | every stage input | numerics |

**The family boundary is the semi-discretisation, not the tableau.** ARK and RK45 agree with each
other 7–9× more closely than either agrees with `split`, despite being an implicit method and an
explicit one — about as different as two integrators get. What they share is the frozen set and the
right-hand side built from it. `split` builds its own frozen set and re-solves hydraulics inside its
iterate. That is where the difference lives.

---

## 5. What the measurements say

All numbers below: one temperate stand (Ithaca, NY), July, one month, hourly output, demography
frozen, and every known model-family difference pinned to the common subset (the harness's
`--parity` preset). Single realisations — read factors, not digits.

### 5a. Everything converges under `dt_fast` refinement

Each scheme scored against a highly-resolved reference (`rk45_tight`, `dt_fast = 25 s`), RMSE of
hourly canopy-air temperature [K]:

| `dt_fast` | split | ark | rk45 |
|---|---|---|---|
| 1800 s | 1.20 | 1.10 | 0.90 |
| 900 s | 0.69 | 0.42 | 0.40 |
| 450 s | 0.52 | **0.14** | **0.14** |
| 225 s | 0.46 | **0.062** | **0.061** |

ARK and RK45 converge cleanly at roughly second order — each halving of `dt_fast` cuts the error by
2.2–3×. **`split` stalls at ~0.45 K.** Against its *own* refined self it converges perfectly well
(1.00 → 0.40 → 0.17 → 0.079), so this is not a failure to converge: split converges to a *different
answer*, about 0.45 K away in canopy-air temperature and 0.15 K in soil-surface temperature. That
residual is a model-family difference, and it is the top open item (§7).

> ⚠ **This supersedes an earlier, wrong result.** A previous version of this sweep reported that
> *nothing* converged and that the error rose again below 900 s. That was an artefact of the
> reference cell, not the model: the diagnostic writer closes a "hourly" record every
> *N* sub-steps with `N = round(3600/dt_fast)`, and the old reference ran at `dt_fast = 27 s`, which
> does not divide 3600. Its records were 3591 s long and drifted ~9 s per hour — about 1.9 h over the
> month — so every cell was compared against a progressively mis-timed reference. Using a reference
> whose step divides both the slow step and one hour removes the effect entirely. The lesson is
> generic and worth keeping: **when a convergence study refuses to converge, audit the reference
> before the model.**

### 5b. At a fixed `dt_fast`, the integrator choice is worth an order of magnitude

Holding `dt_fast = 1800 s` — so all cells share the identical Category-0 freeze — and scoring only
the time-stepping error (reference = `rk45_tight` at the *same* `dt_fast`):

| scheme | CAS T RMSE [K] | wall (net) | sub-steps | rejections |
|---|---|---|---|---|
| **ark** | **0.039** | 3.27 s | 5492 | 1253 |
| rk45 | 0.442 | 3.87 s | 5033 | 913 |
| split | 0.831 | **2.51 s** | 1488 | 0 |
| picard | 0.899 | 5.17 s | 1488 | 0 |
| ark, adaptivity OFF | 1.998 | 2.91 s | 1488 | 0 |

Two readings:

- **Adaptivity is doing real work**, not just spending time: one fixed ESDIRK step per `dt_fast` is
  the worst cell on every metric.
- **`split`'s error at production `dt_fast` is essentially all time-stepping error** (0.83 out of a
  ~0.84 total), while ARK's is essentially all freeze error (0.04 out of ~0.97). Split is limited by
  its stepper; ARK is limited by everything else.

### 5c. Cost

Cost is best read from the sub-step counters (`[output].numerics`), which are exact and
machine-independent, rather than wall clock. At `dt_fast = 1800 s` `split` is ~3.7× cheaper per
`dt_fast` than ARK. But cost per *unit accuracy* runs the other way: `ark` at 900 s (7286 sub-steps,
0.42 K) beats `split` at 225 s (11904 sub-steps, 0.46 K) on both axes at once.

ARK's overhead at coarse `dt_fast` is dominated by **rejected steps** — 1253 rejections at
`dt_fast = 1800 s`, falling to 3 at 225 s. Rejections are wasted work and are the obvious place to
look for cheap speed-ups.

### 5d. Conservation

All three close the whole-column water and energy ledgers to machine precision in the unsaturated
regime that every production run lives in (residuals ~1e-13 kg m⁻² and ~1e-7 J m⁻²), with a pack
present, and with condensation active. The clamp counters (`work_clamp_*`) are zero for `split` by
construction and non-zero but unbookkept-mass-free for the adaptive schemes.

---

## 6. Choosing, in practice

| you want | use |
|---|---|
| long spin-ups, multi-century runs, cost-dominated work | `split` at `dt_fast = 1800 s` |
| sub-daily fidelity — flux-tower comparison, diel cycles, energy partitioning | `ark`, `dt_fast` 900 s or finer |
| a reference solution to check anything else against | `rk45` with `rtol_all = 1e-9`, `atol_scale = 1e-3` |
| ED2 comparison at the algorithmic level | `rk45` (but check `work_rk45_rescue_site` first) |

Configuration essentials:

```toml
[fast]
time_integrator     = "ark"      # "split" (default) | "ark" | "rk45"
integration_scheme  = "split"    # "picard" iterates the leaf<->CAS block on the split path
dt_fast             = "1800s"
ark_rtol            = 1.0e-3     # adaptive schemes only
rtol_all            = 0.0        # master relative-accuracy dial; 0 = per-group defaults
atol_scale          = 1.0        # scales every absolute tolerance -- needed, rtol alone saturates
step_controller     = "i"        # "pi" damps step oscillation on stiff transients
forcing_sample_frac = 0.5        # where in the sub-step the meteorology is sampled

[output]
numerics = true                  # the work + health counters; turn this on for any numerics study

[energy]
debug_error = false              # true = HALT on a non-closing budget. Use it when validating a
                                 # change to any flux seam; leave it off for production.
```

**Two traps worth naming.**

- `rtol_all` alone saturates. Accuracy is judged against `atol + rtol·|y|`; once `atol` dominates,
  tightening `rtol` changes nothing. Scale both.
- Any comparison across schemes must pin the model-family differences in §4 first, or you are
  comparing physics and numerics at once. This has produced at least one wrong conclusion in this
  project's history.

---

## 7. Known limitations and open questions

1. **`split` converges to a different answer** than ARK/RK45 — ~0.45 K in canopy-air temperature,
   ~0.15 K in soil-surface temperature, ~7 W m⁻² in sensible heat, after refinement has removed all
   time-stepping error (§5a). The previous difference of this kind was the missing sapflow advected
   enthalpy, found and fixed; this is the next one down and is not yet attributed.
2. **The ARK's error norm carries structurally-zero terms.** Soil moisture and plant water mass are
   operator-split out of its tableau, so their contributions to the norm are exactly zero and only
   dilute it — the ARK therefore runs slightly looser than its stated `ark_rtol`. Measured cost on a
   summer stand: 0–0.5% in accuracy against 9–15% *fewer* sub-steps. The clean fix is to bring soil
   water into the ARK tableau, which is a substantial change tracked separately.
3. **A shared whole-energy closure failure on forced winter runs.** With the budget hard-stop armed,
   a January Ithaca month halts on all three schemes at a whole-column energy residual of
   −3.5×10³ J m⁻² against a 1.4×10³ J m⁻² tolerance — about 0.4% of a winter day's net radiation. It
   is pre-existing and identical across schemes, so it does not contaminate scheme comparisons, but
   it is a real unclosed seam and is not yet attributed. (It went unnoticed because
   `[energy].debug_error` had no TOML reader until now, so the hard-stop was never armed in a
   configured run.)
4. **RK45's post-commit saturation clip** is a symptom repair. The physically correct treatment is
   the deferred Neumann→Dirichlet ponded-surface switch: once the surface saturates the top boundary
   should become a prescribed head rather than a prescribed flux, so the solve never oversaturates
   and nothing needs clipping. Assumptions on the record in issue #78.
5. **The ponding store holds mass but no heat.** Water crossing into it sheds its enthalpy, which
   leaves the column ledger even though the water is physically still sitting on the surface. Books
   close because both sides agree; the approximation is bounded by the maximum pond depth.
6. **Condensate is deposited into soil layer 1**, not onto foliage, because the canopy film store is
   optional and carries its enthalpy at a fixed reference temperature. Routing dew onto leaves is the
   natural refinement (issue #74) and needs that reference fixed first.
7. **Where the forcing is sampled within the sub-step** (`forcing_sample_frac`) is a real accuracy
   knob and a *trade*, not a free win: sampling at the sub-step start improves canopy-air temperature
   and degrades soil temperature, with opposite optima. It does not cause a convergence floor —
   sweeps at sample points 0, 0.5 and 1 all converge cleanly — but it shifts the error by up to ~3×
   at coarse `dt_fast`. The structurally correct fix is a consistent midpoint freeze, whose predictor
   must be L-stable rather than explicit.
8. **Evidence base.** Everything in §5 is one site, one month per cell, one forcing year, no
   confidence intervals. Differences smaller than ~2× should not be read as real.

---

## 8. Where the code is

| piece | file |
|---|---|
| dispatch + operator-split scheme | `src/driver/meds_fast_split.f90` (`column_fast_step`) |
| IMEX-ARK scheme | `src/driver/meds_fast_ark.f90` |
| adaptive Cash–Karp scheme | `src/driver/meds_fast_rk45.f90` |
| the shared right-hand side | `src/driver/meds_fast_time_derivs.f90` (`column_derivs`) |
| the frozen pre-pass | `meds_fast_ark.f90` (`column_prepass`, `build_column_frozen`) |
| tolerances, error norm, step controller | `src/driver/meds_fast_control.f90` |
| shared snow stage | `src/driver/meds_fast_snow.f90` |
| benchmark harness | `scripts/numerics_sweep.py`, `scripts/parity_fidelity.py` |
