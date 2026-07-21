# Soil biophysics

The soil column is two coupled prognostic sub-columns on one ED2 negative-$z$ grid: a **water** column
(mixed-form Richards) and a **thermal** column (internal-energy heat diffusion). They share the Thomas
sweep and the grid geometry, and are coupled by construction — the thermal step reads the just-updated
moisture, and the moving water carries liquid enthalpy across the soil boundaries.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (**prognostic internal energy**, not temperature; **closed-budget
discipline**) and how the stores are woven each `dt_fast`. Both principles are load-bearing here: the
thermal column is prognostic in volumetric internal energy so freeze/thaw is a read-off of the
thermodynamic inverter, and both sub-columns close a machine-precision budget every step.

---

## The soil water column

`column_hydrology_flux` (`meds_soil_water`) advances one patch's prognostic soil moisture
$`\theta_k`$ (+ ponded surface water and a lumped aquifer) over `dt_fast` with an **implicit
backward-Euler Thomas** solve of the mixed-form Richards equation on the ED2 negative-$z$ grid. Per
layer,

```math
C(\psi_k)\,\frac{\partial\psi_k}{\partial t}
= \frac{\partial}{\partial z}\!\left[K(\theta)\Big(\frac{\partial\psi}{\partial z}+1\Big)\right] - S_k
\qquad(1)
```

with matric potential $`\psi(\theta)`$, hydraulic conductivity $`K(\theta)`$, and specific moisture
capacity $`C=d\theta/d\psi`$ from the **van Genuchten–Mualem** (default) or **Campbell / Clapp-Hornberger**
retention curves (`meds_hydr_lib`). Interface conductivity is **upstream-weighted** by the
total-head gradient; the linearization is either a single **frozen-coefficient** solve or a **Celia
(1990) modified-Picard** iterate (`opts%linearize`), and the step is sub-cycled by adaptive
step-doubling (`soil_water_advance`) — BE is L-stable, so sub-stepping only buys accuracy.

The **root sink** is the transpiration demand distributed over layers and gated by a smooth wilting
ramp $`f_{wilt}(\psi)\in[0,1]`$ between `psi_wilt` and `psi_open`; its $\psi$-derivative enters the
implicit matrix. The **top boundary** is a conductivity-limited infiltration flux with ponding
overflow; the **bottom** is free-drainage (unit-gradient, MVP default), no-flow bedrock, or a SIMTOP
aquifer with a diagnosed water table. **Ground evaporation** combines a Philip pore-space relative
humidity with a Swenson-Lawrence dry-surface-layer (DSL) resistance in series with $`r_{aero}=1/g_{g,net}`$
from the aerodynamics kernel — and it is the **single authority** for the ground latent flux (it drives
both the CAS vapour twin and the ground energy balance's LE, so no double-count). Dunne saturation-excess
runoff runs only under the genuinely-diagnosed aquifer water table.

Every step closes a machine-precision water budget:

```math
\text{mass\_resid} = \Delta W_{stores}
 - \Delta t\,(P - E_{soil} - \text{drainage} - \text{uptake} - \text{runoff}) \approx 0 \qquad(2)
```

and exports the per-layer matric potential `psi_soil` [MPa] that closes the **plant-hydraulics** soil
boundary condition (see [plant_hydraulics](plant_hydraulics.md)), plus the time-mean per-face Darcy
flux for optional advective soil heat.

---

## The soil thermal column

`soil_energy_step_implicit` (`meds_soil_energy`) advances the prognostic volumetric internal energy
$`E_k`$ [J m⁻³] by an implicit BE heat-diffusion solve that **reuses the same Thomas sweep and negative-$z$
geometry** as the hydrology. At each state it inverts $`E_k\to(T_k,\text{fliq}_k)`$, forms the ice-aware
thermal conductivity $`\kappa(\theta,\text{fliq})`$ and effective volumetric heat capacity
$`C_{eff}(\theta,\text{fliq})`$, solves for $`T^{n+1}`$ with a top Neumann flux $`G_{top}`$ and a bottom
geothermal flux, then commits a **conservative** energy update from the $`T^{n+1}`$ conductive faces plus
optional upwind water-enthalpy advection:

```math
E_k^{n+1} = E_k^n + \frac{\Delta t}{\Delta z_k}\big[(hf_k-hf_{k-1}) + (qwf_k-qwf_{k-1})\big]
            + \Delta t\,q_{src,k} \qquad(3)
```

Temperature and liquid fraction are re-diagnosed from the committed energy — so **freeze/thaw is the
internal-energy plateau**, captured for free once `phase_change = ENERGY_PHASE_ON` (ice-aware
$`\kappa`$/$`C_{eff}`$; the zero-curtain is tested, cooling a wet layer pins `soil_temp` at the triple
point while `soil_fliq` absorbs the fusion enthalpy). The closed residual is `energy_resid`
$`=\Delta E - \Delta t\,(G_{top}-\text{bottom}-\sum\text{root\_heat\_sink})\approx0`$. A sibling
`soil_energy_time_deriv` exposes the same flux divergence as an explicit RHS (faces at $`T^n`$) for the
ARK integrator.

## Water–thermal coupling

The water and thermal columns are coupled by construction: the thermal step reads the just-updated soil
moisture (for $`\kappa,C_{eff}`$), and the transpiration/infiltration/drainage water carries liquid
enthalpy across the soil boundaries via `root_heat_sink` and the boundary-face advection.

## Canopy interception

Per-cohort **canopy interception** is a separate top→bottom cascade (`intercept_canopy_layer`, in
`meds_vegetation_biophysics`): a capacity-limited bucket with a Beer interception fraction, each
cohort's throughfall feeding the next, so the net throughfall reaching the soil top is what the
infiltration boundary sees.

## Soil heterotrophic respiration

The heterotrophic soil CO₂ flux $`R_h`$ — Q10 or the DAMM dual-Arrhenius/Michaelis-Menten scheme
(`meds_soil_biogeochem`: `heterotrophic_respiration_flux` / `heterotrophic_respiration_damm`) — is a
soil-carbon decomposition process (the slow biogeochemistry lives in `biogeochemistry/`), but its flux
enters the fast loop as a source term of the CAS CO₂ twin
(see [canopy_air_space_biophysics](canopy_air_space_biophysics.md)). The driver is its single authority
and passes the resulting flux into the CO₂ source.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Soil water | `soil_column_t` | `theta(:)`, `w_surface`, `w_aquifer` | `z_wt`, `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |

## Code map

| Concept | Routine |
|---|---|
| soil water (implicit Richards) | `meds_soil_water`: `column_hydrology_flux`, `soil_water_step_implicit`, `soil_water_advance` |
| explicit water tendency (ARK) | `meds_soil_water`: `soil_water_time_deriv` |
| retention curves (vG / Campbell) | `meds_hydr_lib` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| explicit thermal tendency (ARK) | `meds_soil_energy`: `soil_energy_time_deriv` |
| canopy interception | `meds_vegetation_biophysics`: `intercept_canopy_layer` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_flux`, `heterotrophic_respiration_damm` |
| prognostic soil types | `meds_column_state_types`: `soil_column_t`, `soil_energy_column_t` |

## References
- ED2 `../ED2/ED/src/dynamics/`: `lsm_hyd.f90`, `soil_respiration.f90` — the soil stores in the reference model.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — mixed-form Richards / modified-Picard.
- van Genuchten (1980); Clapp & Hornberger (1978) — soil retention curves.
- Swenson & Lawrence (2014), *JGR-Atmos.* 119:10299 — dry-surface-layer soil evaporation.
- Design docs: `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`.
