# Numerical schemes of the fast loop

How MEDS actually advances the sub-daily biophysics in time, what the two available integrators
assume, where they differ, and which one to use. Written for ecological modellers rather than
numerical analysts: the goal is that you can read a MEDS run's `[fast]` block and know what it
committed you to.

Companion pages: [column_biophysics](column_biophysics.md) describes *what* is being integrated
(the stores and their fluxes). This page is about *how*. Open numerics work is listed in
`docs/ROADMAP.md` (§4).

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
- soil and plant water — the pond, the soil layers and every cohort's wood and leaves, solved
  together as one network by implicit Euler with its own sub-steps
  ([soil_plant_water](soil_plant_water.md)); its **time-mean** flows (face fluxes, infiltration,
  drainage, ground evaporation, root uptake, sapflow, transpiration) handed to the rest of the column
  as constants, on both schemes

This is ED2-faithful. ED2 calls `canopy_photosynthesis` once before its RK integrator and the RK
stages only *read* the frozen conductances; `plant_hydro.f90` freezes hydraulic fluxes per `DTLSM`
for exactly the same reason. Freezing hydraulics removes the ~20 s mode from the coupled system,
which is what makes any of the schemes below affordable.

**Category 1 — stiff, must be treated implicitly**: canopy-air enthalpy / humidity / CO₂, the
leaf and wood energy balance, the ground surface (snow-free, the top soil layer), the soil thermal
column, and the soil–plant water network (pond, soil layers, leaf and wood water).

**Category 2 — non-stiff**: canopy surface (film) water, the transpiration and aerodynamic fluxes
themselves.

**The freeze is not free.** It introduces an error that no time integrator can remove, because it
changes the equations being solved, not their discretisation. Three things follow, and they are the
most useful facts on this page:

1. Refining `dt_fast` refines **both** the freeze and the time step, so a `dt_fast` sweep measures
   them together.
2. Making the integrator more accurate at a *fixed* `dt_fast` eventually stops helping, because the
   freeze error is the floor.
3. **One coefficient is deliberately NOT frozen, because freezing it is unstable.** The
   CAS↔atmosphere conductances $`g_{ah}`$, $`g_{aw}`$, $`g_{ac}`$ are re-solved at **every integrator
   stage**, against that stage's own canopy-air state. Everything else in Category 0 stays frozen.

   The canopy air is a very low-capacity node (`cas_mass_capacity·cp ≈ 2.4×10⁴ J m⁻² K⁻¹`) driven by surface
   fluxes of hundreds of W m⁻², and its ventilation depends on its own temperature through
   atmospheric stability: warmer canopy air ⇒ more unstable surface layer ⇒ larger `ustar` and
   transfer factor ⇒ harder venting. Measured, that feedback has a gain of
   `d ln g_atm_heat/dT ≈ 2.2 K⁻¹`. Evaluated a whole `dt_fast` behind the state it responds to, it turns
   into a **sustained period-2 oscillation** in canopy-air temperature — up to ~8 K step-to-step at
   `dt_fast = 900 s` — while **every conservation budget closes to ~10⁻⁶ J throughout**. That is the
   single most important caveat on this page: *conservation is not stability.*

   Re-solving per stage removes it. The one-step amplification factor of the canopy air (see §5a′)
   falls from −23.2 to −0.14 at 900 s, and step-to-step swing from 7.7 K to 0.10 K:

   | `dt_fast` | frozen conductances | re-solved per stage |
   |---|---|---|
   | 900 s | ~7.7 K | **0.10 K** |
   | 450 s | ~6.6 K | 0.03 K |
   | 150 s | ~2.0 K | 0.012 K |
   | 100 s | ~0.7 K | 0.02 K |

   This is what ED2 does (`update_diagnostic_vars` → `canopy_turbulence8` at every RK stage), and it
   is cheap: `canopy_aerodynamics` is ~2% of the frozen pre-pass against ~89% for leaf gas exchange at
   30 cohorts. It costs a few percent per sub-step and pays for itself in fewer sub-steps. There is no
   switch — the frozen alternative is unstable at production step sizes on most stand heights.

   **The other frozen coefficients were checked and do not matter here.** Holding the stomatal
   conductance, the leaf sensible coefficient, the longwave emission base or the wetted fraction at
   their unperturbed values each changes the amplification factor by ≤ 0.7%, against −23.2 → +0.80 for
   the conductances. See [canopy_aerodynamics](canopy_aerodynamics.md) §2 for the mechanism.

---

## 3. The two integrators

Selected by `[fast].time_integrator`. Both share the same physical kernels, the same frozen pre-pass,
and the same conservation ledgers.

> **Why there is no operator-split scheme.** An operator-split backward-Euler stepper, advancing the
> stores one after another, is the obvious third option, and MEDS does not offer one. It converges to
> a **different limit** than ARK and RK45 (~0.45 K on canopy-air temperature), so a test anchored on
> it is anchored on something that is not the answer. It cannot carry the coupled tissue heat store:
> with correctly sized wood, `cap_wood/cap_cas ≈ 0.5` while `τ_wood ≪ dt_fast`, so tissue and canopy
> air are a coupled stiff pair that needs a joint implicit solve. And a third scheme makes every
> physics addition a three-way wiring exercise. A config asking for `"split"` is a hard error, not a
> silent downgrade.
>
> On provenance: ED2's own default is `INTEGRATION_SCHEME = 1` (RK4). ED2 *does* have an operator
> split — scheme 3, "hybrid", implicit canopy via `bdf2_solver` — but as an alternative, not the
> default. MEDS's `rk45` covers the ED2 default lineage.

### `ark` (default) — 2-solve ESDIRK2

Advance the coupled stiff block with a **second-order, L-stable, two-stage implicit** method
(γ = 1 − 1/√2, the ARS(2,2,2) value), with an **embedded error estimate** that drives adaptive
sub-stepping inside each `dt_fast`.

> **Naming, stated once because it propagates.** Despite the `ark` label this is **not an
> IMEX method.** The biotic CO₂ source is folded implicit, so the explicit tableau is empty
> (`f_E == 0`) and the scheme reduces to a clean 2-solve ESDIRK2. The config string stays `"ark"` for
> compatibility; the *scheme* is ESDIRK2. Do not rename it "IMEX-ARK" — that asserts a property the
> implementation does not have.

*Analogy:* a group ride with a pace-setter. Everyone moves together, and the pace is adjusted from
how far the group drifted apart on the last stretch.

- Canopy air + soil heat are inside the tableau. Soil water, plant water mass and canopy surface
  water are **operator-split out**: soil water moves at the water network's steady rate, and plant
  and surface water are advanced once over the full `dt_fast`.
- The leaf↔canopy-air coupling is solved implicitly *within* every stage by a direct 2×2 Newton
  (`newton_surface_solve`, `[fast].ark_coupled`, default true) with a numerical Jacobian — which is what
  lets it carry the tissue heat store.
- Adaptive: the step size is chosen from the difference between the second-order and first-order
  solutions, measured by **one weighted norm over the whole column state** with a per-variable-group
  tolerance (`[fast].ark_rtol`, `rtol_all`, `atol_scale`). The ARK's soil-moisture and plant-water
  terms in that norm are structurally zero (both ride operator-split maps outside the tableau), so
  they dilute it slightly — accepted deliberately: one norm with no per-scheme opt-outs means no
  scheme can leave a state it *integrates* unmeasured.
- Warm-started: the last accepted sub-step size carries over to the next `dt_fast`, which keeps the
  rejection rate under 20% (about 31% when every `dt_fast` starts cold). Each patch keeps its own
  (`adapt_dt_last`), so a patch's first sub-step depends neither on the patch run before it nor,
  with the patch axis threaded, on the thread schedule.

In symbols, with $h$ the sub-step, $f$ the right-hand side of the stiff block and $`\gamma = 1-1/\sqrt{2}`$,
a step is two backward-Euler solves, the second from the base $`y_n + (Y_2-y_n)(1-\gamma)/\gamma`$:

```math
Y_2 = y_n + \gamma h\,f(Y_2), \qquad
Y_3 = y_n + (1-\gamma)\,h\,f(Y_2) + \gamma h\,f(Y_3), \qquad
y_{n+1} = Y_3 \qquad(1)
```

The step is its last stage (the method is stiffly accurate), and the difference between the two
stages' increments is the free first-order error estimate:

```math
\varepsilon_{n+1} = (Y_3 - \text{base}) - (Y_2 - y_n) = \gamma h\,\big[f(Y_3) - f(Y_2)\big] \qquad(2)
```

Both schemes judge a sub-step by the same weighted root-mean-square of their error estimate
(`meds_fast_control`), each of the $N$ counted state entries scaled by its group's tolerances:

```math
\lVert\varepsilon\rVert
= \sqrt{\frac{1}{N}\sum_{i=1}^{N}\left(\frac{\varepsilon_i}{a_{g(i)} + r_{g(i)}\,|y_i|}\right)^{2}} \qquad(3)
```

$`a_g`$ is the group's absolute tolerance, times `atol_scale`; $`r_g`$ is `ark_rtol` (soil moisture
takes `[fast].water_rtol` and `water_atol`), or `rtol_all` for every group when it is set. The canopy films and the
pond carry no group. A sub-step with $`\lVert\varepsilon\rVert\le 1`$ is accepted, and the next one is
scaled by the I controller or, with `step_controller = "pi"`, the PI controller:

```math
\frac{h_{new}}{h} = \mathrm{clamp}\!\left(s\,\lVert\varepsilon\rVert^{-1/(p+1)},\ f_{min},\ f_{max}\right)
\quad\text{or}\quad
\frac{h_{new}}{h} = \mathrm{clamp}\!\left(s\,\lVert\varepsilon\rVert^{-\alpha}\,
\lVert\varepsilon_{prev}\rVert^{\beta},\ f_{min},\ f_{max}\right)
\qquad(4)
```

with $p$ the order of the embedded estimate (1 for `ark`, 4 for `rk45`), $s = 0.9$, $`f_{min} = 0.2`$,
$`f_{max} = 5`$, $\alpha = 0.35$ and $\beta = 0.20$.

### `rk45` — adaptive explicit Cash–Karp 5(4), ED2-style

A fully explicit fifth-order Runge–Kutta with an embedded fourth-order estimate, marched adaptively
over the *whole* column state — this is the closest analogue of ED2's own `rk4` integrator. Kept as
the **accuracy baseline**, deliberately not optimised.

*Analogy:* the same group ride, but nobody is allowed to look ahead; if the road gets steep the only
option is to take smaller steps.

- Everything — including soil moisture and plant water mass — is inside the tableau. Soil moisture
  moves there at the water network's steady rate, which any consistent tableau integrates exactly, so
  RK45's committed θ is the network's, as ARK's is.
- **It is a hybrid, not a pure explicit scheme.** When the explicit march cannot resolve a step (a
  cold night under a closed canopy, where the leaf↔canopy-air coupling is stiff and the state rails
  against its physical bounds), the driver **discards the RK45 step and redoes that `dt_fast` on
  `ark`**. The counter `work_rk45_rescue_site` reports how often. *Read it before reading any RK45
  result* — a run whose winter months all rescue is substantially an ARK run wearing an RK45 label.
  ARK is the rescue target because what is being rescued from *is* a coupled stiff canopy-air
  problem, which is exactly what its Newton solves.

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
| per-layer root-sink placement | tracks realized uptake | same | unified — the per-layer path is unconditional |
| bottom boundary (drainage, free to sealed / aquifer) | supported | same | unified — head-driven two-way aquifer, no storage bucket |
| bottom thermal boundary | `geothermal` (prescribed flux, held at 0) or `dirichlet` (deep temperature anchor) | same | unified — one `bottom_heat_face` seam, see [soil_biophysics](soil_biophysics.md) |
| canopy-air depth | tallest cohort + freeboard, slow-loop owned | same | unified |
| soil water | the water network's, moved at its steady rate outside the tableau | the network's, at the same rate inside the tableau | unified — both commit the network's θ |
| canopy-air condensation sink | per-stage rate | per-stage rate | unified |
| adaptive error norm | one norm over the whole column state | same | unified |
| stability clamps | extrapolation base only (canopy air) | every stage input (canopy air) | numerics |

**There are no *physics* rows.** The two schemes solve the same equations; what remains is one
numerics row (where clamps apply).

**The family boundary is the semi-discretisation, not the tableau.** ARK and RK45 agree closely (§5)
despite being an implicit method and an explicit one — about as different as two integrators get.
What they share is the frozen set and the right-hand side built from it. That is where the difference
lives, and it is why §2's freeze is the first thing on this page.

---

## 5. What the measurements say

**Cross-scheme accuracy is measured at stable cadences.** `ark` and `rk45` agree to **1.3 mK** on
canopy-air temperature at `dt_fast = 150 s` with tightened tolerances — an implicit and an explicit
method landing on the same answer. With the conductances frozen the same comparison reads 0.49 K and
does not shrink when the tolerance is tightened, because it measures the oscillation rather than
the schemes. That is worth stating as a rule: **a cross-scheme difference measured in an unstable
regime measures the instability, not the schemes**, and a measurement at `dt_fast` 900–1800 s with
the conductances frozen is a phase-sample of the oscillation, not a convergence measurement.

Two further results from the same sweep. **The freeze dominates the residual error**: tightening the
integrator tolerance at `dt_fast = 150 s` moves `ark` from 0.0174 K to 0.0180 K, so ~97% of what is
left is the Category-0 freeze and the integrator tolerance is not the lever. And **`ark` converges at
roughly second order** in `dt_fast`; with the conductances frozen there is no order to measure at
all, only a discontinuity at the stability boundary.

### 5a. `dt_fast` is an ACCURACY parameter

With the conductances re-solved per stage the canopy air is stable at every `dt_fast` measured, and
the sensitivity lies in the **carbon and water budget** — because leaf gas exchange, and the leaf
water potential that sets stomatal conductance, are frozen per step.

The temperatures, on a high-LAI sunlit stand over 24 h, scored against a 12.5 s reference:

| `dt_fast` | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| CAS-T RMSE | 0.017 K | 0.037 K | 0.055 K | 0.085 K | 0.16 K |
| soil-T RMSE | 0.06 K | 0.14 K | 0.22 K | 0.37 K | 0.70 K |

*The quadrature is sound.* GPP is accumulated as `Σ gpp(state n)·dt`, a left rectangle in the
state, so an averaging artefact is the obvious suspect. It is not: on the same stand, with leaf
capacity held independent of ψ, that accumulation reproduces the 12.5 s answer to **0.05% at
`dt_fast = 900 s`** (ET to 1.2%).

*Leaf water potential is the sensitive input, because the fast water-stress control amplifies it.*
The non-stomatal limb `beta = (ψ − ψ_close)/(ψ_open − ψ_close)` is a linear ramp on Vcmax/Jmax/TPU in
the step-start leaf potential, from the turgor-loss point to twice it (`leaf_gas_exchange.md` §4.1):
its slope is about 0.6 MPa⁻¹ on the default `[hydraulics]` curve, so a 1 MPa error in ψ becomes a
~60% error in capacity. An explicit plant-water update with frozen sapflow and uptake makes errors of
that size, growing with the step: on the same stand it runs the daytime-mean leaf potential at
−1.19 MPa at 900 s against −0.23 MPa at 12.5 s (a per-step |Δψ| of up to 0.84 MPa), and GPP 33% and
ET 24% low. The limb is an **amplifier, not the source**, and its clamp at `beta = 0` is never
reached, so this is straight linear amplification rather than a saturation artefact.

**That is why leaf and wood water are solved implicitly, with the soil's, in one network**
([soil_plant_water](soil_plant_water.md)), and why the energy stages hold each cohort's transpiration
at the network's for the step: within a step transpiration does not follow the canopy air's humidity,
and the canopy air's vapour credit is the leaves' water debit. Measured on BCI's dry season
(example 04, 2013-01-15 to 04-15, `dt_fast` 900 against 150 s): the window-mean latent heat differs by
−0.2 W m⁻² and the mean diurnal cycle by at most 1.9 W m⁻² in any hour (refreshing transpiration at
every stage instead gives +0.3 and 2.9 W m⁻²; at 150 s the two agree to 0.1 W m⁻² in the mean and 0.6
in any hour). The 3-hour midday probe's `psi_leaf` agrees to 0.005 MPa between 75 and 900 s
(`test_column_ark`). Over a wet and a dry season at BCI (example 04 on the tower forcing, 2012-08 to
2013-05, 900 against 150 s) GPP differs by +0.08 %, LE by −0.37 W m⁻² and H by +0.25 W m⁻², and the
daily leaf and wood potentials by at most 0.004 MPa.

**A drought season carries a bias at 900 s.** Measured on Palo Verde's inventory configuration
(example 05), restarted from a 150 s state on 2008-12-01 and run at 150, 300, 450 and 900 s:

| month | GPP at 300 / 450 / 900 s vs 150 s | ET at 900 s | `psi_wood` at 900 s |
|---|---|---|---|
| Jan | +0.1 / +0.3 / +0.6 % | −0.2 % | +0.024 MPa |
| Feb | +1.0 / +2.0 / +5.1 % | +5.4 % | +0.076 MPa |
| Mar | +1.0 / +1.9 / +5.4 % | +3.9 % | +0.032 MPa |
| Apr | +0.9 / +2.3 / +6.6 % | +8.1 % | +0.051 MPa |

At 900 s February–April GPP is 5.5 % and ET 6.0 % high, the December–May season 1.7 % and the year
0.8 %, and the drought-deciduous trees shed their leaves 1–2 days late. The bias is linear in
`dt_fast`, about 1 % at 300 s. It is the per-step freeze, not the water solver (tightening the
solver's tolerances five-fold changes nothing): the gas exchange and the transpiration demand are set
from the step-start leaf potential, and the fast water-stress limb amplifies the potential's change
within the step. **For a run whose purpose is the timing of dry-season leaf habits, use
`dt_fast = 300 s`**; 900 s stays the default, since a wet site and the annual budgets barely move.

### 5a′. The stability multiplier, and how to measure it

The quantity that decides stability is the one-step amplification factor of the canopy air,
`Phi' = ∂H_{n+1}/∂H_n`, measured by perturbing canopy-air enthalpy *before* the pre-pass runs and advancing
one `dt_fast`. `|Phi'| > 1` means the freeze-cadence map is unstable, and a negative `Phi'` makes a
perturbation alternate in sign from step to step. With the conductances frozen, on five stand
heights at `dt_fast = 150 s`:

| stand | 2 m regen | 5 m | 10 m | 18 m | 30 m |
|---|---|---|---|---|---|
| `Phi'` frozen | −0.90 | −3.84 | **−8.22** | −4.84 | −4.88 |

With the conductances re-solved per stage, as both schemes run, a 30 m stand of LAI 4 at noon:

| `u_ref` | 50 s | 100 s | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|---|---|
| 2 m/s | +0.39 | +0.08 | −0.10 | −0.06 | −0.03 | +0.03 | +0.11 |
| 10 m/s | +0.26 | +0.04 | −0.04 | +0.02 | +0.03 | +0.03 | +0.02 |

`|Phi'|` is at most 0.42 in every other case tried, except a calm night at 1 m/s, where it is +0.91 to
+0.43: a slow relaxation that does not oscillate. The re-solved canopy air is stable at every
`dt_fast`, in light wind as in strong. A base state for this measurement must itself be spun up with
the conductances re-solved; one spun up frozen carries the frozen map's state into the multiplier.

Two things worth knowing from the frozen row. **The worst case is a mid-height canopy, not a short
one** — a regenerating patch has less canopy-air capacity but also much less leaf area driving it,
and `Phi'` tracks the ratio. And **amplitude is a bad proxy for stability**: the 2 m stand has the
mildest multiplier at 900 s and the largest excursion (20 K), because the nonlinearity saturates.
Rank by `Phi'`, not by peak-to-peak.

### 5b. Conservation

Both schemes close the whole-column water and energy ledgers to machine precision in the unsaturated
regime every production run lives in (residuals ~10⁻¹³ kg m⁻² and ~10⁻⁷ J m⁻²), with a snow pack
present and with condensation active. This holds with the tissue heat store active: the heat each
tissue's energy gains is taken from a b-weighted time integral of its temperature, so it is exactly
the complement of the flux the canopy air received.

**Read 5a and 5b together.** They are the same run. A ledger that closes to 10⁻⁶ J tells you the
bookkeeping is right; it tells you nothing about whether the trajectory is physical.

The rest of this section is the rules the ledgers rest on, and what each guards against.

**Every soil-water flow comes from one solve, whose θ both schemes commit.** A scheme that took a
*flux* from one solve while committing *state* from its own stages would carry an error of exactly
however far the two trajectories diverge — invisible while they nearly coincide, dominant once the
column saturates. A **whole-column** ledger cannot detect that error at all, because it is purely
*vertical*: enthalpy placed in the wrong layer still sums correctly against the boundary. Only a
per-layer face check, or a temperature that stops being plausible, exposes it. Taking the face
fluxes, infiltration, drainage and uptake from the water network, and committing the network's θ,
makes that class impossible for water.

**Saturation is the water network's, and nothing is clipped.** The network treats the pond as a store
with a head and each layer above saturation as storing water through a specific storage of
10⁻⁴ m⁻¹, so a saturated column raises its potential and drives the excess to the pond and the next
layer within the step, on both schemes. θ exceeds θ_sat only by that storage: on the sealed
29 mm h⁻¹ fixture the worst excursion `theta_ood_max` reports is 5.5×10⁻⁶. A clip would remove water
with no face to carry it, and so no enthalpy to go with it.

**The pond holds heat.** It carries a prognostic enthalpy, and every seam moves mass and enthalpy as a
pair: rain in at its own temperature, snowmelt in at the pack's, infiltration out at the *mixed* pond
temperature, exfiltration from the top layer in at that layer's temperature, and runoff out at the
pond's. A pond with no thermal state would shed the enthalpy of every kilogram crossing into it,
and that energy would leave the column ledger while the water sat on the surface still holding it.
The pond's energy (`w_surface_energy`, J m⁻²) is stored like the snow's, and its temperature and
liquid fraction are read off it the same way, so a pond freezes and thaws with no code of its own.

**A thin tissue closes its own balance.** A leaf or wood store's relaxation closes for any coupling,
however small, since it is written on the start-of-step imbalance with no division by the coupling;
and a tissue with neither coupling nor heat capacity (no area) takes the temperature of what it grows
from and carries no flux ([vegetation_energy_dynamics](vegetation_energy_dynamics.md) eqs. (5), (7)).
No tissue needs a floor that routes energy elsewhere. With the budget hard-stop armed
(`[fast].debug_error`), a forced January month at Ithaca (a deciduous stand in winter) closes on
both schemes: none of the 8928 whole-column checks fails, and the worst residual is 1.2×10⁻⁶ J m⁻²
(ARK and RK45). Tightening the tolerance makes the run halt, which shows the guard is live.

---

## 6. Choosing, in practice

`dt_fast` is not bounded by stability (§2). It is bounded by **how much flux error you can
accept** (§5a), so a spin-up and a flux-tower study genuinely want different values.

| you want | use |
|---|---|
| general production, including spin-ups and coupled carbon–water | `ark` at `dt_fast = 900 s` (the default) |
| sub-daily fidelity — flux-tower comparison, diel cycles, energy partitioning | `ark`, `dt_fast` ≤ 150 s |
| anything keyed to **leaf water potential** (hydraulic stress, ψ-driven mortality) | `ark` at 900 s for the leaf and wood potentials (§5a); check the first day after a restart at a shorter step |
| the timing of **dry-season leaf habits** (when drought-deciduous trees shed and flush) | `ark`, `dt_fast = 300 s`: 900 s runs a drought season's GPP and ET ~6 % high (§5a) |
| a reference solution to check anything else against | `rk45` with `rtol_all = 1e-9`, `atol_scale = 1e-3`, `dt_fast` ≤ 50 s |
| ED2 comparison at the algorithmic level | `rk45` (but check `work_rk45_rescue_site` first) |

**What the `dt_fast = 900 s` row costs.** Little: canopy-air temperature within 0.16 K of a resolved
reference, and the latent heat within 0.2 W m⁻² in the mean and 1.9 W m⁻² in any hour of the diurnal
cycle on BCI's dry season. Soil-surface temperature is the loosest at 0.70 K. Under drought the fast
water-stress limb reads the step-start leaf potential, which runs Palo Verde's February–April GPP and
ET 5.5–6 % high at 900 s, 1.7 % over the season and 0.8 % over the year (§5a).

Configuration essentials:

```toml
[fast]
time_integrator     = "ark"      # "ark" (default, ESDIRK2) | "rk45" (accuracy baseline)
dt_fast             = "900s"     # ACCURACY-limited (sec 5a). Warned above 900 s, the measured
                                 # range; canopy-air temperature is not the constraint. 300 s for
                                 # the timing of dry-season leaf habits.
ark_rtol            = 1.0e-3     # adaptive schemes only
rtol_all            = 0.0        # master relative-accuracy dial; 0 = per-group defaults
atol_scale          = 1.0        # scales every absolute tolerance -- needed, rtol alone saturates
step_controller     = "i"        # "pi" damps step oscillation on stiff transients
forcing_sample_frac = 0.5        # where in the sub-step the meteorology is sampled
debug_error         = false      # true = HALT on a non-closing budget. Use it when validating a
                                 # change to any flux seam; leave it off for production.

[output]
numerics = true                  # the work + health counters; turn this on for any numerics study
```

**Three traps worth naming.**

- `rtol_all` alone saturates. Accuracy is judged against `atol + rtol·|y|`; once `atol` dominates,
  tightening `rtol` changes nothing. Scale both. (And per §5, tolerance is rarely the binding lever —
  the freeze is.)
- Any comparison across schemes must be made at a *stable* cadence, or it measures the instability
  rather than the schemes.
- **A closed budget is not a healthy run.** With the conductances frozen (§2) every ledger closes to
  ~10⁻⁶ J while canopy-air temperature oscillates 8 K step to step. If you change `dt_fast`, look at a
  consecutive-step trace of `cas_temp_site` and at daily GPP, not just the residuals.

---

## 6a. Multi-core: threading the patch axis

The fast loop's columns are **independent within a `dt_fast`** — patches couple only through the slow
loop — so patches are the parallel axis, and threading them costs no accuracy at all. The default
build compiles the threading in (`-DMEDS_OPENMP=OFF` leaves it out), and a run asks for threads:

```toml
[run]
n_threads = 4                    # default 1. Ignored by a serial build.
```

**The output is byte-identical at every thread count**, and that is a design requirement rather than a
happy accident. Every site-level accumulator (ET, the daily-mean air temperature, the ten integrator
work counters, the sub-daily output staging) is written into a per-`(sub-step, patch)` staging array
and folded back **in patch order** after the loop. An OpenMP `reduction(+:)` would sum thread partials
in *arrival* order, so the last bits of every site diagnostic would drift with the thread count. Folding
in patch order reproduces the serial fold exactly, so the answer does not move with the thread count.
Verified over a 12-patch / 115-cohort forced month: 35 netCDF files across all four output tiers,
identical at 1, 2 and 4 threads; and over the five-year BCI census example at 1, 4, 8 and 16 threads.
An OpenMP build is not bit-identical to a serial build, because compiling the loop as a parallel region
changes its rounding (`docs/building.md`, "Parallel builds").

One knob is deliberately incompatible: `[fast].fast_probe` writes a single, order-significant CSV from
a saved unit, so `fast_probe` together with `n_threads > 1` is a hard configuration error rather than a
silently interleaved file.

### What it buys, measured

A 50-year spin-up at Ithaca NY on ERA5-Land forcing (bare-ground start, `ark`, 4 threads on 4 physical cores):

| `dt_fast` | wall | CPU | final stand |
|---|---|---|---|
| 150 s | **2321.8 s** (38 min 42 s) | 7622.9 s (3.28 cores avg) | 113 cohorts / 12 patches |
| 900 s | **545.1 s** (9 min 05 s) | 1510.6 s (2.77 cores avg) | 115 cohorts / 12 patches |

**The 900 s default is worth 4.26× here, not the 6× the step-size ratio suggests** — the ARK march
takes *two* sub-steps at 900 s where it takes one at 150 s, so 6× fewer steps each cost ~1.4× more, and
the slow tier does not scale at all. That matches the 4.0–4.6× measured independently per simulated
day.

Threading itself gives ~2.0× on 4 cores for this run. Read that against the machine, not the core
count: four *independent* single-threaded copies of the same workload take 55.6 s each where one alone
takes 42.3 s, so this hardware's four-core aggregate throughput is only **3.03×** its single-core
throughput (turbo, shared L3, memory). Against what the hardware can actually deliver the loop is ~67%
efficient. **Always measure that denominator before concluding a parallel speedup is poor.**

Speedup is also bounded by the patch count, and a spin-up from bare ground is the worst case for it:
this run held **2 patches** for its first ~15 simulated years, where four threads can do nothing. CPU
utilisation climbed 188% → 272% → 328% as the stand built from 2 patches to 12. A production run
starting from an established stand has the full parallelism from step one.

### The step-size difference this comparison exposes

The two runs above differ **only** in `dt_fast`, so everything below is the semi-discretisation's
50-year footprint. Every patch-area-weighted site aggregate agrees to **≤ 0.5%**: stem density 0.40%,
basal area 0.31%, AGB 0.44%, LAI 0.04%, soil moisture 0.39%, leaf temperature 0.02%, soil temperature
0.04%. (AGB is 0.44% *lower* at 900 s.)

**But the demography took a different path**, and this is the part to carry away. The runs end with 113
vs 115 cohorts, and at year 35 they held 82 vs 67 before reconverging; patch ages differ by ~8% and
cohort identities differ entirely. `dt_fast` perturbs GPP, which perturbs growth, which changes *which*
cohorts fuse, split or are culled. That is a **discrete** difference, not a truncation error that
shrinks with the step — so two runs at different `dt_fast` cannot be compared cohort-by-cohort, only
through site aggregates, and the agreement of those aggregates in one realization is reassurance rather
than a convergence proof.

---

## 7. Known limitations and open questions

1. **The per-step freeze biases a drought season, and a restart's first day is not scored.** The
   leaf and wood potentials converge in `dt_fast` at a wet site, but at 900 s Palo Verde's driest
   months run GPP and ET 5.5–6 % high, linearly in `dt_fast` (§5a). The candidate fix is a
   predictor–corrector: re-evaluate the gas exchange at the mid-step leaf potential and re-solve the
   water network once, about one extra pre-pass per step. The start-up transient on the first day
   after a restart has not been re-measured on the network (see `docs/ROADMAP.md`).

2. **Light wind is no worse.** The feedback behind the oscillation is stability-dependent and its
   gain varies with wind, so light-to-moderate wind, the common band, is the one that matters. On BCI's
   tower forcing (example 04, 2012-08 to 2013-05, fast output every step) calm wind below 2 m/s is
   13 % of the records, and at 900 s the step-to-step alternation of canopy-air temperature is no
   larger than the 150 s run's at the same instants in any wind band: the lag-1 correlation of its
   second difference is −0.06, where a period-2 oscillation gives −1. §5a′ gives the multiplier
   across wind speeds and `dt_fast`. It is scored at one site.

3. **The ARK's error norm carries structurally-zero terms.** Soil moisture and plant water mass are
   operator-split out of its tableau, so their contributions to the norm are exactly zero and only
   dilute it — the ARK therefore runs slightly looser than its stated `ark_rtol`. Measured cost on a
   summer stand: 0–0.5% in accuracy against 9–15% *fewer* sub-steps. Bringing soil water into the ARK
   tableau would remove the dilution, at a cost (see `docs/ROADMAP.md`).
4. **Condensate is deposited into soil layer 1**, not onto foliage. It is small, at most 0.5 % of
   evapotranspiration at BCI, so routing it onto the leaves would change little, and would need the
   canopy film, which is optional and carries its enthalpy at a fixed reference temperature, to carry
   its own temperature first.
5. **Where the forcing is sampled within the sub-step** (`forcing_sample_frac`) is a real accuracy
   knob and a *trade*, not a free win: sampling at the sub-step start improves canopy-air temperature
   and degrades soil temperature, with opposite optima. It does not cause a convergence floor —
   sweeps at sample points 0, 0.5 and 1 all converge cleanly — but it shifts the error by up to ~3×
   at coarse `dt_fast`. A consistent midpoint freeze would remove the trade; its predictor would have
   to be L-stable rather than explicit.
6. **Evidence base.** Most of §5 rests on one site and one forcing year per measurement, with no
   confidence intervals. Differences smaller than ~2× should not be read as real. The saturating-column
   checks of §5b run on a *sealed unit fixture* driven at 29 mm h⁻¹, far above Ithaca's rain rates,
   because a forced Ithaca column never saturates: a production-forced test suite is not a substitute
   for a fixture that deliberately visits the corner.

---

## 8. Where the code is

| piece | file |
|---|---|
| dispatch + the RK45→ARK stiff rescue | `src/fast_dynamics/numerics/meds_fast_step.f90` (`column_fast_step`) |
| ESDIRK2 scheme (config name `ark`), eqs. (1)–(2) | `src/fast_dynamics/numerics/meds_fast_ark.f90` (`ark2_column_step`, `adaptive_ark_march`) |
| adaptive Cash–Karp scheme | `src/fast_dynamics/numerics/meds_fast_rk45.f90` |
| the shared right-hand side | `src/fast_dynamics/numerics/meds_fast_time_derivs.f90` (`column_derivs`) |
| the frozen pre-pass | `src/fast_dynamics/numerics/meds_fast_prepass.f90` (`column_prepass`), `meds_fast_frozen.f90` (`build_column_frozen`) |
| the implicit stage and its leaf↔canopy-air Newton | `src/fast_dynamics/numerics/meds_fast_be_stage.f90` (`column_be_stage`, `newton_surface_solve`) |
| the soil–plant water network | `src/fast_dynamics/numerics/meds_spac_solver.f90` ([soil_plant_water](soil_plant_water.md)) |
| error norm and step controller, eqs. (3)–(4) | `src/fast_dynamics/numerics/meds_fast_control.f90` (`state_wrms_grouped`, `step_control_factor`); the tolerance set from config: `src/fast_dynamics/driver/meds_fast_config.f90` (`build_tol_set`) |
| shared snow stage | `meds_fast_frozen.f90` (`advance_snow_stage`) |
| patch loop, per-thread scratch pool, order-preserving reductions (§6a) | `src/fast_dynamics/driver/meds_fast_dynamics.f90` (`fast_dynamics`) |
