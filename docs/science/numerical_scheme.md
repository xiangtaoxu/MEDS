# Numerical schemes of the fast loop

How MEDS actually advances the sub-daily biophysics in time, what the two available integrators
assume, where they differ, and which one to use. Written for ecological modellers rather than
numerical analysts: the goal is that you can read a MEDS run's `[fast]` block and know what it
committed you to.

Companion pages: [column_biophysics](column_biophysics.md) describes *what* is being integrated
(the stores and their fluxes). This page is about *how*. The engineering record — measurements,
per-difference status, open work — lives in
`docs/dev_plans/docs/dev_plans/archive/MEDS_INTEGRATOR_PARITY.md [RETIRED]` and `docs/dev_plans/MEDS_NUMERICS_SCOPING.md`.

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
changes the equations being solved, not their discretisation. Three things follow, and they are the
most useful facts on this page:

1. Refining `dt_fast` refines **both** the freeze and the time step, so a `dt_fast` sweep measures
   them together.
2. Making the integrator more accurate at a *fixed* `dt_fast` eventually stops helping, because the
   freeze error is the floor.
3. **Past a threshold the freeze is not merely inaccurate — it is unstable.** The canopy air is a
   very low-capacity node (`wcap·cp ≈ 2.4×10⁴ J m⁻² K⁻¹`) driven by surface fluxes of hundreds of
   W m⁻²: 300 W m⁻² over a 900 s step is an **11 K excursion**. Holding its coupling coefficients
   fixed across a step that long feeds a lagged canopy-air temperature back into its own balance, and
   the result is a **sustained period-2 oscillation**:

   | `dt_fast` | canopy-air peak-to-peak, consecutive steps |
   |---|---|
   | 900 s | ~8 K |
   | 225 s | ~4 K |
   | 150 s | ~1 K (diurnal trend only) |
   | 100 s | smooth |

   **No conservation ledger detects this** — every budget closes to ~10⁻⁶ J throughout. That is the
   single most important caveat on this page: *conservation is not stability.* The default `dt_fast`
   is therefore **150 s**, and `meds_config` prints a warning above 300 s.

   It is the **aggregate** of the frozen couplings, not any one of them. Pinning the aerodynamic
   conductances leaves ~6 K; pinning the stomatal conductance leaves ~9 K; pinning the entire surface
   coupling still leaves ~6.7 K of the ~9.5 K baseline. Refreshing only part of the pre-pass per
   sub-step therefore does **not** buy back a longer `dt_fast` — that was measured and rejected.
   ED2 sidesteps the question by recomputing canopy turbulence at *every* RK stage
   (`update_diagnostic_vars` → `canopy_turbulence8`).

   The threshold is stand-dependent: it scales with the flux-to-capacity ratio, so a short or
   regenerating patch (smaller canopy-air volume) is **stricter** than 150 s, not looser.

---

## 3. The two integrators

Selected by `[fast].time_integrator`. Both share the same physical kernels, the same frozen pre-pass,
and the same conservation ledgers.

> **A third scheme, `split`, was retired (2026-07-31.)** It was an operator-split backward-Euler
> stepper and the historical default. Three reasons: it converged to a **different limit** than
> ARK/RK45 (~0.45 K on canopy-air temperature), so every "ARK agrees with split" test was anchored on
> something that is not the answer; it could not carry the coupled tissue heat store (with correctly
> sized wood, `cap_wood/cap_cas ≈ 0.5` while `τ_wood ≪ dt_fast`, making tissue and canopy air a
> coupled stiff pair that needs a joint implicit solve); and it made every physics addition a
> three-way wiring exercise. The `[fast].integration_scheme` selector (its Gauss–Seidel-vs-Picard
> sweep) went with it. A config still asking for `"split"` is a hard error, not a silent downgrade.
>
> For the record on provenance: ED2's own default is `INTEGRATION_SCHEME = 1` (RK4). ED2 *does* have
> an operator split — scheme 3, "hybrid", implicit canopy via `bdf2_solver` — but as an alternative,
> not the default. MEDS's `rk45` covers the ED2 default lineage.

### `ark` (default) — 2-solve ESDIRK2

Advance the coupled stiff block with a **second-order, L-stable, two-stage implicit** method
(γ = 1 − 1/√2, the ARS(2,2,2) value), with an **embedded error estimate** that drives adaptive
sub-stepping inside each `dt_fast`.

> **Naming, stated once because it propagates.** Despite the historical `ark` label this is **not an
> IMEX method.** The biotic CO₂ source is folded implicit, so the explicit tableau is empty
> (`f_E == 0`) and the scheme reduces to a clean 2-solve ESDIRK2. The config string stays `"ark"` for
> compatibility; the *scheme* is ESDIRK2. Do not rename it "IMEX-ARK" — that asserts a property the
> implementation does not have.

*Analogy:* a group ride with a pace-setter. Everyone moves together, and the pace is adjusted from
how far the group drifted apart on the last stretch.

- Canopy air + soil heat are inside the tableau. Soil water, plant water mass and canopy surface
  water are **operator-split out** and advanced once over the full `dt_fast`.
- The leaf↔canopy-air coupling is solved implicitly *within* every stage by a direct 2×2 Newton
  (`newton_surface_solve`, `[fast].ark_niter`, default 8) with a numerical Jacobian — which is what
  lets it carry the tissue heat store that the retired split path could not.
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
over the *whole* column state — this is the closest analogue of ED2's own `rk4` integrator. Kept as
the **accuracy baseline**, deliberately not optimised.

*Analogy:* the same group ride, but nobody is allowed to look ahead; if the road gets steep the only
option is to take smaller steps.

- Everything — including soil moisture — is inside the tableau. It is the only scheme that integrates
  its own soil water rather than taking it from the pre-pass solve.
- **It is a hybrid, not a pure explicit scheme.** When the explicit march cannot resolve a step (a
  cold night under a closed canopy, where the leaf↔canopy-air coupling is stiff and the state rails
  against its physical bounds), the driver **discards the RK45 step and redoes that `dt_fast` on
  `ark`**. The counter `work_rk45_rescue_site` reports how often. *Read it before reading any RK45
  result* — a run whose winter months all rescue is substantially an ARK run wearing an RK45 label.
  (The rescue target used to be `split`; ARK is the better one on its own merits, since what is being
  rescued from *is* a coupled stiff canopy-air problem, which is exactly what its Newton solves.)

---

## 4. What differs between them, and what kind of difference it is

Three categories, because they demand different responses. A *physics* difference means the schemes
solve different equations and no refinement will reconcile them; a *numerics* difference shrinks as
you refine; a *bookkeeping* difference is a conservation defect.

| difference | ark | rk45 | kind |
|---|---|---|---|
| snow | shared stage, always active | same | unified |
| canopy surface water | present | present | unified |
| sapflow advected enthalpy | present | present | unified |
| tissue (leaf + wood) heat store | exact exponential relaxation | same | unified |
| per-layer root-sink placement | tracks realized uptake | same | unified — `multilayer_roots` deleted, the per-layer path is unconditional |
| bottom boundary (free-drain / bedrock / aquifer) | supported | same | unified — head-driven two-way aquifer, no storage bucket |
| bottom thermal boundary | adiabatic (`geothermal ≡ 0`) | same | unified — see [soil_biophysics](soil_biophysics.md) |
| canopy-air depth | tallest cohort + freeboard, slow-loop owned | same | unified |
| soil water | frozen pre-pass solve, outside the tableau | integrated in the tableau | numerics + assumptions |
| canopy-air condensation sink | per-stage rate | per-stage rate | unified |
| adaptive error norm | one norm over the whole column state | same | unified |
| stability clamps | extrapolation base only | every stage input | numerics |

**There are no *physics* rows left.** The two schemes solve the same equations; what remains is one
assumptions row (soil water inside vs outside the tableau) and one numerics row (where clamps apply).
That is a deliberate outcome — the parity work that closed the last rows is recorded in
`docs/dev_plans/MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md`.

Two of those rows closed by being **rebuilt rather than propagated**, which is worth knowing because
the fix was not "make the other scheme do what this one does":

- **the aquifer BC** used to apply the deep-water-table limit `q = K(θ_n)` unconditionally — identical
  to free drainage — behind a lumped storage bucket whose water table could rise into the soil column
  without saturating it. It is now a head-driven two-way boundary against a saturated zone at the
  column base, with the bucket, its baseflow, the Dunne runoff parameterization and Zeng–Decker all
  deleted.
- **the canopy-air depth** was a hardcoded 20 m that nothing ever assigned; the aerodynamics computed
  the right value and discarded it. It is now per-patch state owned by the slow loop.

**The family boundary is the semi-discretisation, not the tableau.** ARK and RK45 agree with each
other far more closely than either agreed with the retired `split` path, despite being an implicit
method and an explicit one — about as different as two integrators get. What they share is the frozen
set and the right-hand side built from it. That is where the difference lives, and it is why §2's
freeze is the first thing on this page.

---

## 5. What the measurements say

The cross-scheme accuracy tables that used to fill this section have been **retired with the `split`
scheme they were anchored on** — see `docs/dev_plans/archive/MEDS_INTEGRATOR_PARITY.md [RETIRED]` for
the record and its tombstone. They were invalid twice over: they scored `ark` and `rk45` against
`split`, which converged to a different, never-attributed limit, and they were taken at
`dt_fast` 900–1800 s, inside the oscillating regime of §2, so they are phase-samples rather than
convergence measurements. **Re-measurement at `dt_fast = 150 s` on two schemes is outstanding.**

What is measured and current:

### 5a. `dt_fast` is a stability boundary

Canopy-air temperature, peak-to-peak over *consecutive* steps, high-LAI sunlit stand:

| `dt_fast` | p2p | verdict |
|---|---|---|
| 900 s | ~8 K | sustained period-2 oscillation |
| 450 s | ~8 K | oscillating |
| 300 s | ~5 K | oscillating |
| 225 s | ~4 K | oscillating |
| 150 s | ~1 K | diurnal trend only |
| 100 s | — | smooth |

Every conservation ledger closed to ~10⁻⁶ J at every row. See §2 for the mechanism and for the
elimination experiments that show it is the *aggregate* of the frozen couplings.

### 5b. Conservation

Both schemes close the whole-column water and energy ledgers to machine precision in the unsaturated
regime every production run lives in (residuals ~10⁻¹³ kg m⁻² and ~10⁻⁷ J m⁻²), with a snow pack
present and with condensation active. This holds with the tissue heat store active: its energy is
carried as a b-weighted time integral so the store's gain is exactly the complement of the flux the
canopy air received.

**Read 5a and 5b together.** They are the same run. A ledger that closes to 10⁻⁶ J tells you the
bookkeeping is right; it tells you nothing about whether the trajectory is physical.

---

## 6. Choosing, in practice

| you want | use |
|---|---|
| anything, including long spin-ups | `ark` at `dt_fast = 150 s` |
| sub-daily fidelity — flux-tower comparison, diel cycles, energy partitioning | `ark`, `dt_fast` 150 s or finer |
| a reference solution to check anything else against | `rk45` with `rtol_all = 1e-9`, `atol_scale = 1e-3` |
| ED2 comparison at the algorithmic level | `rk45` (but check `work_rk45_rescue_site` first) |

There is no longer a "cheap but rough" option: `dt_fast` is bounded from above by **stability**, not
by an accuracy preference (§2), so a long spin-up pays the same step a diel study does. Buying speed
by lengthening `dt_fast` buys an oscillation that no budget will report.

Configuration essentials:

```toml
[fast]
time_integrator     = "ark"      # "ark" (default, ESDIRK2) | "rk45" (accuracy baseline)
dt_fast             = "150s"     # STABILITY-limited (sec 2). Warned above 300 s.
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

**Three traps worth naming.**

- `rtol_all` alone saturates. Accuracy is judged against `atol + rtol·|y|`; once `atol` dominates,
  tightening `rtol` changes nothing. Scale both.
- Any comparison across schemes must pin the model-family differences in §4 first, or you are
  comparing physics and numerics at once. This has produced at least one wrong conclusion in this
  project's history.
- **A closed budget is not a healthy run.** Every ledger closed to ~10⁻⁶ J while canopy-air
  temperature oscillated 8 K step to step for a year of simulated time. If you change `dt_fast`, look
  at a consecutive-step trace of `cas_temp_site`, not just the residuals.

---

## 7. Known limitations and open questions

1. **The freeze-cadence oscillation has no identified single cause.** §2 records that it is the
   aggregate of the frozen surface couplings and that no partial refresh fixes it, but the specific
   loop that sustains it is *not* known. Eliminated by measurement: aerodynamic conductances,
   stomatal conductance, the whole surface coupling together, the plant-hydraulics store
   (`mask%hydraulics`), the soil thermal column (`mask%soil_heat`), and the soil water column
   (`mask%soil_water`) — every one still oscillates at `dt_fast = 900 s`. Also ruled out: the
   leaf↔canopy-air coupling is not explicit (the Newton runs by default). Untested suspects, in
   order: the leaf gas-exchange pre-pass (`gsw` is solved once per `dt_fast` at state *n* and feeds
   `g_tr_f`, so pinning `g_tr_f` removed its *variation* while the pre-pass still ran lagged),
   `t_emit`, and `f_wet`. Until this is closed, `dt_fast ≤ 150 s` is a constraint, not a tuning knob.
2. **The two schemes have never been cross-checked at a stable `dt_fast`.** Everything in the
   archived comparison was measured at 900–1800 s, inside §2's oscillating regime, and against a
   third scheme that has since been retired for converging to a different limit. A clean `ark`-vs-
   `rk45` convergence study at 150 s and finer is outstanding, and until it exists there is no
   current quantitative statement of how far apart they are.
3. **The ARK's error norm carries structurally-zero terms.** Soil moisture and plant water mass are
   operator-split out of its tableau, so their contributions to the norm are exactly zero and only
   dilute it — the ARK therefore runs slightly looser than its stated `ark_rtol`. Measured cost on a
   summer stand: 0–0.5% in accuracy against 9–15% *fewer* sub-steps. The clean fix is to bring soil
   water into the ARK tableau, which is a substantial change tracked separately.
4. **The shared whole-energy closure failure on forced winter runs is RESOLVED and attributed.** With
   the budget hard-stop armed, a January Ithaca month used to halt on every scheme at a
   whole-column energy residual of −3.5×10³ J m⁻² against a 1.4×10³ J m⁻² tolerance — about 0.4% of a
   winter day's net radiation, and shared rather than scheme-specific. The cause was the
   `veg_coupling_floor` clamp destroying energy in the diagnostic leaf/wood balance; reverting that one
   fix alone still reproduces the halt (−1.5×10³ against 1.4×10³), which is what attributes it. Because
   the clamp lives in `surface_derivs`, which both schemes share, the failure was identical on all
   both — the very symmetry that made it look like a deep seam is what identifies it as a single shared
   kernel defect. A January month now closes with roughly 5× margin on every scheme: worst residual
   285 (ARK), 287 (RK45) J m⁻² against a ~1.39×10³ J m⁻² tolerance. The guard is verified live
   in that configuration by tightening the tolerance and confirming it halts. (It went unnoticed
   originally because `[energy].debug_error` had no TOML reader, so the hard-stop was never armed in a
   configured run.)
5. **RK45's post-commit saturation clip** is a symptom repair. The physically correct treatment is
   the deferred Neumann→Dirichlet ponded-surface switch: once the surface saturates the top boundary
   should become a prescribed head rather than a prescribed flux, so the solve never oversaturates
   and nothing needs clipping. Assumptions on the record in issue #78. The *worst* consequence of the
   clip is now gone (issue #78 item 3, below), but the clip itself remains the mechanism, because an
   explicit method has no post-solve hook the way the implicit sibling does. A related worry — that the
   soil retention curves get *evaluated outside their domain* while a stage sits above saturation —
   turned out not to be real: every constitutive kernel clamps its own effective saturation, so an
   oversaturated cell is treated as exactly saturated (ψ = 0, K = K_sat), which is the physically right
   answer for one. The excursion itself is real but small (worst measured 0.008 in θ, i.e. Se = 1.022,
   on a 29 mm h⁻¹ fixture; identically zero in every forced Ithaca cell) and is now reported per
   sub-step as `theta_ood_max` so it cannot grow unnoticed.
6. **Borrowing one solve's numbers while committing another's trajectory** is the single most
   productive defect class this review found, and worth stating as a rule rather than as three
   separate bugs. Each fast-loop scheme freezes an "Act-1" scratch hydrology solve and then advances
   the column; wherever a scheme took a *flux* from the scratch while committing *state* from its own
   stages, the two disagreed by exactly however far the two trajectories had diverged — invisible while
   they nearly coincide, dominant once the column saturates. Three instances have now been found and
   fixed, all on RK45: its drainage (fixed by taking the b-weighted stage value), its interior
   advective enthalpy faces, and the saturation-clip enthalpy it paid for mass its own θ never shed.
   The last two were of matched magnitude and *opposite* sign (~2.6×10⁶ J m⁻² per step each on a
   saturating column), so each concealed the other: removing either alone drove the soil surface from
   285 K to 345 K, and only removing both together gave a physical 295 K. The general lesson is that a
   **whole-column** conservation ledger cannot detect this class at all, because the error is purely
   *vertical* — enthalpy placed in the wrong layer still sums correctly against the boundary. Only a
   per-layer face check, or a temperature that stops being plausible, exposes it.
7. **The ponding store now holds heat** (issue #78 item 4, resolved). It used to be a mass buffer with
   no thermal state, so water crossing into it shed its enthalpy and that energy left the column ledger
   while the water sat on the surface still holding it. The store now carries a prognostic enthalpy and
   every seam moves mass and enthalpy as a pair: rain in at its own temperature, infiltration out at the
   *mixed* pond temperature, the saturation clip in at each layer's temperature, and runoff out at the
   final pond temperature. Snowmelt was rerouted here too, which deleted a correction the ARK and RK45
   paths both used to need — a net simplification, not just an addition.
8. **Condensate is deposited into soil layer 1**, not onto foliage, because the canopy film store is
   optional and carries its enthalpy at a fixed reference temperature. Routing dew onto leaves is the
   natural refinement (issue #74) and needs that reference fixed first.
9. **Where the forcing is sampled within the sub-step** (`forcing_sample_frac`) is a real accuracy
   knob and a *trade*, not a free win: sampling at the sub-step start improves canopy-air temperature
   and degrades soil temperature, with opposite optima. It does not cause a convergence floor —
   sweeps at sample points 0, 0.5 and 1 all converge cleanly — but it shifts the error by up to ~3×
   at coarse `dt_fast`. The structurally correct fix is a consistent midpoint freeze, whose predictor
   must be L-stable rather than explicit.
10. **Evidence base.** Everything in §5 is one site, one month per cell, one forcing year, no
   confidence intervals. Differences smaller than ~2× should not be read as real. Two of the defects
   above were found only on a *sealed unit fixture* driven at 29 mm h⁻¹ — far above this site's rain
   rates — and are measurably no-ops in every Ithaca cell, because those columns never saturate. A
   production-forced test suite is not a substitute for a fixture that deliberately visits the corner.

---

## 8. Where the code is

| piece | file |
|---|---|
| dispatch + the RK45→ARK stiff rescue | `src/driver/meds_fast_step.f90` (`column_fast_step`) |
| ESDIRK2 scheme (config name `ark`) | `src/driver/meds_fast_ark.f90` |
| adaptive Cash–Karp scheme | `src/driver/meds_fast_rk45.f90` |
| the shared right-hand side | `src/driver/meds_fast_time_derivs.f90` (`column_derivs`) |
| the frozen pre-pass | `meds_fast_ark.f90` (`column_prepass`, `build_column_frozen`) |
| tolerances, error norm, step controller | `src/driver/meds_fast_control.f90` |
| shared snow stage | `src/driver/meds_fast_snow.f90` |
| benchmark harness | `scripts/numerics_sweep.py`, `scripts/parity_fidelity.py` |
