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

`veg_energy_step_implicit` (`meds_vegetation_biophysics`) advances a cohort's leaf **or** wood prognostic
internal energy over one **L-stable linearized backward-Euler** step. It diagnoses $`T^n`$, forms the
net flux $`R^n = A_{sw}+A_{lw} - H - Q_w - Q_{transp}`$ (absorbed SW + net LW − sensible −
film-evaporation − transpiration), and the linearization slope $`\partial R/\partial T = \frac{dR}{dT}\le0`$
(radiative + sensible + latent-saturation terms, all $\le0$). The implicit temperature estimate stays
within one Newton step of $`T^n`$, so the energy update is **bounded** even when the store is stiff:

```math
T^* = T^n + \frac{R^n}{C/\Delta t - dR/dT}, \qquad
\Delta E = C\,(T^*-T^n) \qquad(1)
```

Using the bounded $`C(T^*-T^n)`$ rather than the endpoint flux $`\Delta t\,R^*`$ is what makes it
L-stable: applying $`\Delta t\,R^*`$ over the whole step over-removes energy from a stiff store and can
drive $T$ negative (a young-stand wood cohort has time constant $`\tau=C/|dR/dT|\sim35`$ s vs
$`\Delta t=1800`$ s). Leaf uses a $2\times$LAI flat-plate sensible area; wood uses $`\pi\times`$WAI
cylinders and no transpiration. The flux report is step-averaged and split by the shares at $`T^*`$, and
`energy_resid` closes to zero.

## One formulation: the tissue relaxes exactly

There is no diagnostic/prognostic *mode*. The tissue always carries a heat capacity and
`veg_energy_diagnostic` always relaxes it — **exactly**, because under the Category-0 coefficient
freeze the equation is linear with a known timescale:

```math
\mathrm{cap}\,\frac{dT}{dt}=\mathrm{numer}-\mathrm{denom}\,(T-T_{cas}),\qquad \tau=\frac{\mathrm{cap}}{\mathrm{denom}}
```

The step therefore has a closed form, and the kernel uses it rather than discretising it. Two
**different** weights fall out, and both are needed:

```math
x=\frac{\Delta t}{\tau}=\frac{\mathrm{denom}}{a_{\rm store}},\qquad
w_{\rm end}=e^{-x},\qquad
w_{\rm avg}=\frac{1-e^{-x}}{x}
```

`w_end` weights the old state at the **endpoint** — the committed temperature. `w_avg` weights it in
the **step average** — every flux the canopy air actually receives. Using one for both is what breaks
conservation; pairing them closes the balance identically:

```math
a_{\rm store}\,(\Delta T_{\rm end}-\Delta T_{\rm prev})+\mathrm{denom}\cdot\Delta T_{\rm avg}=\mathrm{numer}
```

**"Diagnostic" is now this kernel's zero-inertia limit, not a second code path.** As
$`a_{\rm store}\to 0`$ (no heat capacity) $`x\to\infty`$, both weights vanish, and the result is the
pure steady-state balance. As $`a_{\rm store}\to\infty`$ both weights → 1 and the tissue holds its
temperature. Nothing switches; the physics chooses, per cohort, every step, continuously.

### The store is ON, and how its capacity is sized

`TISSUE_STORE_SCALE = 1` (`meds_fast_ark`). The capacity splits its two halves across **different**
masses, which is a physical distinction and not a detail:

- **dry tissue → ALL the wood.** Heartwood is dead structure but still stores sensible heat, and a
  sapwood fraction defined on the bole under-counts branch wood, which is thin enough to be thermally
  active throughout.
- **internal water → the SAPWOOD RING only.** Heartwood is taken as dry, which is what makes it
  heartwood. The same ring is the hydraulic capacitance — one quantity, two consumers.

Both come from the ED2 allometry (`dbh_to_wai`, `sapwood_fraction`; PFT traits `wai_b1/b2`,
`sapwood_area_b1/b2`). Sized this way, `τ_wood ≈ 1200–1600 s` and `w_end ≈ 0.23–0.32` at
`dt_fast = 1800 s` — wood carries **genuine memory**. `τ_leaf ≈ 12–20 s` and is LAI-independent (cap
and denom both scale with LAI), so the leaf is near its zero-inertia limit at any production step.

**A warning that cost real time to learn.** Every earlier answer for `τ_wood` was an artefact of the
*mass*, never of the method — placeholders gave 55–199 s, the sapwood ring alone gave 74–102 s. And
because `w_end` is small for wood while `w_avg` is not, **the store's effect lives in the FLUXES, not
in the temperature**: "the wood barely lags" does not imply "the wood barely matters". An earlier
conclusion that prognostic wood ≈ diagnostic wood came from looking only at the temperature.

Measured effect, from a 2×2 over (store, longwave) on a 2 h night window: the store leaves the night
canopy air **1.57 K warmer** (leaf + wood carry ~1.5×10⁴ J m⁻² K⁻¹ against the canopy air's ~3.0×10⁴,
so shedding a ~3 K excess should lift the air by ~0.5 × 3 K — and it does). It also absorbs most of
the net-longwave signal, which is why any test thresholding on nighttime cooling has to be re-pinned
against the store-on control rather than nudged.

**`veg_energy_step_implicit` is no longer called by any driver** — it survives only as a unit-tested
kernel. It stepped a backward-Euler store and carried the longwave response in its Jacobian but not
its residual; the exact exponential above replaced it.

**Why not backward Euler**, which is what the store used before. BE is the
$`w_{\rm end}=w_{\rm avg}=1/(1+x)`$ approximation. A leaf at production `dt_fast` has $`x\approx144`$:
the exact endpoint weight is $`\sim 5\times10^{-63}`$ against BE's $`0.0069`$, and an SDIRK2 tableau
gives $`-0.031`$ — a *sign-alternating* artificial memory for a mode whose true memory is nil. Under
BE the prognostic tissue's only measurable effect was the discretisation's own artefact.

Note `denom` **excludes** $`a_{\rm store}`$: it is the coupling conductance — a property of the
tissue's surroundings, not of its own inertia. Since $`a_{\rm store}=\mathrm{cap}/\Delta t`$, the
ratio $`\mathrm{denom}/a_{\rm store}`$ is $`\Delta t/\tau`$ and the step size cancels, so the kernel
needs no `dt` argument and stays a pure tendency-safe routine.

### What the timescales actually are

$`\tau_{\rm leaf}\approx 12.5`$ s, and — a useful result — it is **LAI-independent**: `cap` and
`denom` both scale with LAI, so every leaf cohort shares one τ and there is no per-cohort structure to
resolve. A leaf is genuinely at equilibrium at any `dt_fast` the model runs at.

Wood is different, and the number in the code is not yet the physical one. The model's `wai` and
`bsap` are hardcoded placeholder ratios, which make the thermal wood a 0.6–2.3 mm shell and give
$`\tau_{\rm wood}\approx 55`$–199 s. Sized physically (branch wood plus the bole's diurnal damping
depth) it is ≈1300 s for branches and ≈3690 s for boles — so wood memory at `dt_fast = 1800 s` is
real, not negligible. Wood is also treated as **internally isothermal**, a deliberate simplification:
Bi ≈ 7.45 for a bole skin against the Bi ≪ 0.1 lumped capacitance needs, so real wood is a diffusion
problem with no single τ. Both are sizing/structure limitations, not integrator ones.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` (internal energy) | temp / fliq read-offs |

## Code map

| Concept | Routine |
|---|---|
| leaf/wood diagnostic solve (shared) | `meds_vegetation_biophysics`: `veg_energy_diagnostic` |
| leaf/wood exact exponential relaxation | `meds_vegetation_biophysics`: `veg_energy_diagnostic` (the `a_store` path) |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |

## References
- ED2 `../ED2/ED/src/dynamics/rk4_derivs.f90` — the leaf/wood energy budget in the reference model.
- Design doc: `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`,
  `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.
