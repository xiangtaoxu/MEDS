# Vegetation energy dynamics

The per-cohort **leaf and wood** tissue energy balance: how each cohort's canopy surfaces exchange
radiation, sensible heat, and water vapour with the canopy air space, and how their temperatures are
either diagnosed (quasi-steady) or advanced as prognostic internal-energy stores.

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

## Diagnostic vs prognostic modes

Two **modes** are config-selectable (`docs/dev_plans/MEDS_LEAF_WOOD_ENERGY_DESIGN.md`):

- **Diagnostic** (default): the leaf/wood takes a linearized steady-state balance
  ($`T_l=T_{cas}+\Delta T_l`$, the leaf inertia is a negligible ~0.015 K). This is the single shared
  solve `veg_energy_diagnostic`, which both integrators call — the split sweep and the ARK surface path.
  Wood is the $`\text{le\_slope}=\text{le\_ref}=0`$ (no-transpiration) case of the same kernel.
- **Prognostic**: advances the tissue internal-energy store via (1). A prognostic leaf **requires**
  `integration_scheme="picard"` — the explicit leaf↔CAS split oscillates (the leaf is stiff, its
  sensible+latent feedback on the CAS has a fixed-point slope near $-1$, giving ~1.7 K midday spikes; the
  Picard iterate damps them to ~0.2 K). Prognostic leaf under the ARK integrator is deferred (it needs
  the implicit leaf↔CAS arrowhead stage solve).

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` (internal energy) | temp / fliq read-offs |

## Code map

| Concept | Routine |
|---|---|
| leaf/wood diagnostic solve (shared) | `meds_vegetation_biophysics`: `veg_energy_diagnostic` |
| leaf/wood prognostic store (BE) | `meds_vegetation_biophysics`: `veg_energy_step_implicit` |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |

## References
- ED2 `../ED2/ED/src/dynamics/rk4_derivs.f90` — the leaf/wood energy budget in the reference model.
- Design doc: `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`,
  `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.
