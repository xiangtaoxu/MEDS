# Vegetation energy dynamics

The per-cohort **leaf and wood** tissue energy balance: how each cohort's canopy surfaces exchange
radiation, sensible heat, and water vapour with the canopy air space, and how their temperatures are
relaxed exactly over the step toward their quasi-steady balance, with their own heat capacity.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven each `dt_fast`. Absorbed shortwave and net longwave enter from the radiation kernel
([canopy_radiation_transfer](canopy_radiation_transfer.md)); the boundary-layer conductances come from
the aerodynamics kernel ([canopy_aerodynamics](canopy_aerodynamics.md)); the transpiration demand comes
from leaf gas exchange ([leaf_gas_exchange](leaf_gas_exchange.md)) and is supply-limited by plant
hydraulics ([plant_hydraulics](plant_hydraulics.md)).

---

## Leaf and wood energy balance

`veg_energy_balance` (`meds_plant_biophysics`) is the one kernel for a cohort's leaf **or** wood; wood is
the case without transpiration. Its conductances and absorbed radiation are frozen for the step
(Category 0, [numerical_scheme](numerical_scheme.md) §2), the saturation terms follow the canopy-air
state, and every flux is linearized about the canopy-air temperature $`T_{cas}`$. The balance is then
linear in the tissue's offset $`\Delta T = T - T_{cas}`$: a forcing $`\mathrm{numer}`$ that does not
depend on $`\Delta T`$, and a coupling conductance $`D_0`$ that multiplies it,

```math
\mathrm{numer} = A_{sw} + A_{lw} + Q_{adv} - \mathrm{LE}_{ref}, \qquad
D_0 = h + \mathrm{LE}_{slope} + 4\,\varepsilon\,\sigma\,T_{cas}^3\,\mathrm{AI}, \qquad
\mathrm{denom} = \max(D_0,\ D_{min}) \qquad(1)
```

$`A_{sw}`$ and $`A_{lw}`$ are the absorbed shortwave and net longwave from the two-stream; $`Q_{adv}`$ is
the enthalpy the sapflow carries in (`q_extra`); $\mathrm{AI}$ is LAI for leaves and WAI for wood, and
$\varepsilon$ the PFT's leaf or wood emissivity. The sensible coefficient is
$`h=\rho\,c_p\,g_{bh}\cdot2\,\mathrm{LAI}`$ for leaves (flat plates, both sides, `effarea_heat`) and
$`h=\rho\,c_p\,g_{bh}\cdot\pi\,\mathrm{WAI}`$ for wood (cylinders). Each latent pathway with conductance
$g$ linearizes the saturation deficit at $`T_{cas}`$:

```math
\mathrm{LE}_{ref} = h_v\,\rho\,g\,\big(q_{sat}(T_{cas}) - q_{cas}\big), \qquad
\mathrm{LE}_{slope} = h_v\,\rho\,g\,\frac{dq_{sat}}{dT}\bigg|_{T_{cas}} \qquad(2)
```

with $`h_v`$ the specific enthalpy of the vapour at $`T_{cas}`$, so the tissue pays what the canopy air
receives. The terms in eq. (1) sum the pathways, each weighted by its share of the surface: a leaf's
dry share $`1-f_{wet}`$ transpires through boundary layer and stomata in series; its wetted share
$`f_{wet}`$ evaporates the intercepted film through the boundary layer alone, paying $`h_v`$ less the
enthalpy the film's liquid arrived with. Wood has only the film pathway. Once the water network has
fixed the step's transpiration, the dry share carries it as a fixed latent term
($`\mathrm{LE}_{slope}=0`$), so the transpiration the leaf reports is exactly the water it lost.

**The coupling floor.** $`D_{min}`$ = 1 W m⁻² K⁻¹ (`veg_coupling_floor`) bounds $`\Delta T`$ for a
cohort whose area index, and so every conductance in $`D_0`$, is tiny. The deficit
$`g_{slave}=\mathrm{denom}-D_0`$ is not dropped: it is an extra conductance to the canopy air that the
sensible flux carries, so a tissue too weakly coupled to resolve passes what it absorbs straight to the
CAS — ED2's `leaf_resolvable`/`wood_resolvable`, made continuous — and the balance still closes. The
floor acts on the denominator rather than as a branch on the result, so the right-hand side stays
continuous in the state (a jump there keeps an embedded-error controller from ever accepting a step),
and it is inactive for any established canopy.

The tissue relaxes toward the zero-inertia offset $`\Delta T_{diag}=\mathrm{numer}/\mathrm{denom}`$ with
its own heat capacity (next section), and reports its fluxes over the step average, eqs. (5)–(6).

**The emission base.** The two-stream emits each cohort at its effective radiative temperature
$`T_{rad}=\big(f_l T_{leaf}^4+(1-f_l)\,T_{wood}^4\big)^{1/4}`$, $`f_l=\mathrm{LAI}/(\mathrm{LAI}+\mathrm{WAI})`$,
from the start-of-step tissue temperatures ([canopy_radiation_transfer](canopy_radiation_transfer.md)
§3), so the absorbed longwave $`A_{lw}`$ already contains the tissue's own emission at $`T_{rad}`$.
`veg_energy_balance` linearizes the change in emission about $`T_{cas}`$ (every caller passes
`t_emit` $`=T_{cas}`$), so the loss it books is $`A_{lw}-4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T-T_{cas})`$
where the consistent one is $`A_{lw}-4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T-T_{rad})`$: a tissue warmer
than the canopy air loses an extra $`4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T_{rad}-T_{cas})`$ and a
cooler one gains it. Basing the emission on the lagged tissue temperature (`t_emit` $`=T_{rad}`$)
removes the residual but adds a lagged positive feedback (a warmer leaf raises the next step's $`T_{rad}`$), which drove a single-pass
integrator to NaN when it was tried (`docs/dev_plans/archive/MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §5) and has
not been tried under `ark`.

## One formulation: the tissue relaxes exactly

There is no diagnostic/prognostic *mode*. The tissue always carries a heat capacity and
`veg_energy_balance` always relaxes it — **exactly**, because under the Category-0 coefficient
freeze the equation is linear with a known timescale:

```math
\mathrm{cap}\,\frac{dT}{dt}=\mathrm{numer}-\mathrm{denom}\,(T-T_{cas}),\qquad \tau=\frac{\mathrm{cap}}{\mathrm{denom}}
\qquad(3)
```

The step therefore has a closed form, and the kernel uses it rather than discretising it. Two
**different** weights fall out, and both are needed:

```math
x=\frac{\Delta t}{\tau}=\frac{\mathrm{denom}}{a_{\rm store}},\qquad
w_{\rm end}=e^{-x},\qquad
w_{\rm avg}=\frac{1-e^{-x}}{x} \qquad(4)
```

with $`a_{\rm store}=\mathrm{cap}/\Delta t`$. `w_end` weights the old state at the **endpoint** — the
committed temperature. `w_avg` weights it in the **step average** — every flux the canopy air actually
receives:

```math
\Delta T_{\rm end}=(1-w_{\rm end})\,\Delta T_{\rm diag}+w_{\rm end}\,\Delta T_{\rm prev},\qquad
\Delta T_{\rm avg}=(1-w_{\rm avg})\,\Delta T_{\rm diag}+w_{\rm avg}\,\Delta T_{\rm prev} \qquad(5)
```

where $`\Delta T_{\rm prev}`$ is the start-of-step offset. The reported sensible heat, evaporation and net
radiation all use $`\Delta T_{\rm avg}`$:

```math
H = (h + g_{slave})\,\Delta T_{\rm avg}, \qquad
E = \frac{\mathrm{LE}_{ref}+\mathrm{LE}_{slope}\,\Delta T_{\rm avg}}{h_v}, \qquad
R_{net} = A_{sw} + A_{lw} - 4\,\varepsilon\,\sigma\,T_{cas}^3\,\mathrm{AI}\,\Delta T_{\rm avg} \qquad(6)
```

($E$ per latent pathway; $`Q_{adv}`$ stays out of $`R_{net}`$ because it is an internal soil-to-plant
transfer, not a radiative input at the column boundary). Using one weight for both is what breaks
conservation; pairing them closes the balance identically:

```math
a_{\rm store}\,(\Delta T_{\rm end}-\Delta T_{\rm prev})+\mathrm{denom}\cdot\Delta T_{\rm avg}=\mathrm{numer}
\qquad(7)
```

**"Diagnostic" is this kernel's zero-inertia limit, not a second code path.** As
$`a_{\rm store}\to 0`$ (no heat capacity) $`x\to\infty`$, both weights vanish, and the result is the
pure steady-state balance. As $`a_{\rm store}\to\infty`$ both weights → 1 and the tissue holds its
temperature. Nothing switches; the physics chooses, per cohort, every step, continuously.

### The store is ON, and how its capacity is sized

The store is always on (`build_column_frozen` in `meds_fast_frozen`). The capacity splits its two
halves across **different** masses, which is a physical distinction and not a detail:

- **dry tissue → ALL the wood.** Heartwood is dead structure but still stores sensible heat, and a
  sapwood fraction defined on the bole under-counts branch wood, which is thin enough to be thermally
  active throughout.
- **internal water → the SAPWOOD RING only**, at a fixed moisture. Heartwood is taken as dry, which
  is what makes it heartwood. This is not the hydraulic store: the water network's wood store is the
  sapwood **plus the fine roots**, and carries its own prognostic water
  ([plant_hydraulics](plant_hydraulics.md) §1), which the heat capacity does not read (see
  `docs/ROADMAP.md`).

Per unit ground area, with $B$ the dry biomass (carbon × 2 × $`n_{plant}`$):

```math
\mathrm{cap}_{wood} = \max\!\big(c_{wood}\,B_{wood},\ C_{min}\big) + c_{liq}\,f_{moist}\,B_{sap}, \qquad
\mathrm{cap}_{leaf} = \max\!\big(c_{leaf}\,B_{leaf},\ C_{min}\big) + c_{liq}\,W_{leaf} \qquad(8)
```

with $`c_{wood}`$ = 2700 and $`c_{leaf}`$ = 3200 J kg⁻¹ K⁻¹ (`c_sapw`, `c_leaf`),
$`C_{min}`$ = 20 J m⁻² K⁻¹ (`veg_hcap_min`), a fresh-sapwood moisture $`f_{moist}`$ = 1 kg water per kg
dry wood, and $`W_{leaf}`$ the leaf's internal water from the hydraulics. The intercepted film is a
separate store and is not in $`\mathrm{cap}_{leaf}`$. The wood masses come from the ED2 allometry
(`dbh_to_wai`, `sapwood_fraction`; PFT traits `wai_b1/b2`, `sapwood_area_b1/b2`). Sized this way,
`τ_wood ≈ 1200–1600 s` and `w_end ≈ 0.23–0.32` at `dt_fast = 1800 s` — wood carries **genuine memory**.

**`τ_wood` is set by the mass, not by the method** — the sapwood ring alone would give 74–102 s. And
because `w_end` is small for wood while `w_avg` is not, **the store's effect lives in the FLUXES, not
in the temperature**: "the wood barely lags" does not imply "the wood barely matters", and comparing
tissue temperatures alone makes the store look negligible.

Measured effect, from a 2×2 over (store, longwave) on a 2 h night window: the store leaves the night
canopy air **1.57 K warmer** (leaf + wood carry ~1.5×10⁴ J m⁻² K⁻¹ against the canopy air's ~3.0×10⁴,
so shedding a ~3 K excess should lift the air by ~0.5 × 3 K — and it does). It also absorbs most of
the net-longwave signal, which is why a test that thresholds nighttime cooling is pinned against the
store-on control.

**Why not backward Euler.** BE is the $`w_{\rm end}=w_{\rm avg}=1/(1+x)`$ approximation. A leaf at
production `dt_fast` has $`x\approx144`$: the exact endpoint weight is $`\sim 5\times10^{-63}`$ against
BE's $`0.0069`$, and an SDIRK2 tableau gives $`-0.031`$ — a *sign-alternating* artificial memory for a
mode whose true memory is nil. Under BE the prognostic tissue's only measurable effect would be the
discretisation's own artefact.

Note `denom` **excludes** $`a_{\rm store}`$: it is the coupling conductance — a property of the
tissue's surroundings, not of its own inertia. Since $`a_{\rm store}=\mathrm{cap}/\Delta t`$, the
ratio $`\mathrm{denom}/a_{\rm store}`$ is $`\Delta t/\tau`$ and the step size cancels, so the kernel
needs no `dt` argument and stays a pure tendency-safe routine.

### What the timescales actually are

$`\tau_{\rm leaf}\approx 12`$–20 s, and — a useful result — it is **LAI-independent**: `cap` and
`denom` both scale with LAI, so every leaf cohort shares one τ and there is no per-cohort structure to
resolve. A leaf is genuinely at equilibrium at any `dt_fast` the model runs at.

Wood is different. The allometric sizing above gives $`\tau_{\rm wood}\approx 1200`$–1600 s; a physical
estimate (branch wood plus the bole's diurnal damping depth) is ≈1300 s for branches and ≈3690 s for
boles — so wood memory at `dt_fast = 1800 s` is real, not negligible. Wood is also treated as
**internally isothermal**, a deliberate simplification: Bi ≈ 7.45 for a bole skin against the Bi ≪ 0.1
lumped capacitance needs, so real wood is a diffusion problem with no single τ. That is a limitation of
the one-node wood store, not of the integrator.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` [K] (heat store, exact exponential relaxation) | — |

## Code map

| Concept | Routine |
|---|---|
| leaf/wood balance + exact exponential relaxation, eqs. (1), (3)–(7) | `meds_plant_biophysics`: `veg_energy_balance` |
| coefficients $h$, $g$, emission slope | `sensible_heat_coeff`, `leaf_transp_coeff`, `leaf_film_coeff`, `lw_emission_slope` |
| latent terms, eq. (2); the per-cohort call | `meds_fast_time_derivs`: `surface_derivs` |
| tissue heat capacity, eq. (8), and $`a_{\rm store}`$ | `meds_fast_frozen`: `build_column_frozen` |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |

## References
- ED2 `../ED2/ED/src/dynamics/rk4_derivs.f90` — the leaf/wood energy budget in the reference model.
- Design records (archived, `docs/dev_plans/archive/`): `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` (the
  exact-exponential tissue store), `MEDS_ENERGY_BALANCE_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.
