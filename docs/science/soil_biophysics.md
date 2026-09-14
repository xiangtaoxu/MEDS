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

`advance_soil_water_column` (`meds_soil_water`) advances one patch's prognostic soil moisture
$`\theta_k`$ (+ ponded surface water) over `dt_fast` with an **implicit
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

The **root sink** is the transpiration demand distributed **per layer** — tracking the realized uptake
the plant-hydraulics solve returns, unconditionally (the `multilayer_roots` switch is deleted; there is
no root-fraction-weighted fallback) — and gated by a smooth wilting ramp $`f_{wilt}(\psi)\in[0,1]`$
between `psi_wilt` and `psi_open`, whose $\psi$-derivative enters the implicit matrix. **Ground
evaporation** combines a Philip pore-space relative humidity with a Swenson-Lawrence
dry-surface-layer (DSL) resistance in series with $`r_{aero}=1/g_{g,net}`$ from the aerodynamics kernel
— and it is the **single authority** for the ground latent flux (it drives both the CAS vapour twin and
the ground energy balance's LE, so no double-count).

### Boundary conditions, and what the aquifer option actually means

The **top boundary** is a conductivity-limited infiltration flux with ponding overflow.

The **bottom** (`[soil].bottom_bc`) is one of three, and the third was rebuilt rather than tuned:

| option | flux at the base | use for |
|---|---|---|
| `free_drain` (default) | $`q = K(\theta_n)`$, unit gradient, always **downward** | a deep, unseen water table |
| `bedrock` | $`q = 0`$ | an impermeable base |
| `aquifer` | $`q = K_{bot}\,(\psi_n/\Delta + 1)`$, $`\Delta = \Delta z_n/2`$, **two-way** | a shallow water table at the column base |

**`aquifer` is a boundary condition, not a store.** The water table is *defined* to sit at the column
base, so the modelled soil column **is** the unsaturated zone above it. There is no aquifer bucket, no
baseflow, and no prognostic water-table depth — the lumped store, its baseflow, the diagnosed $`z_{wt}`$,
the Dunne saturation-excess runoff that keyed off it, and the Zeng–Decker equilibrium correction have
all been **deleted**. $`K_{bot}`$ is upstream-weighted: $`K(\theta_n)`$ when the flux is downward,
$`K_{sat}`$ when it reverses, because the upstream cell is then the saturated zone (using
$`K(\theta_n)`$ there under-predicts capillary rise).

Two things follow that are easy to get wrong:

- The flux **reverses upward** once $`|\psi_n| > \Delta`$, i.e. once the deepest layer is dry enough to
  pull against the saturated zone. So `aquifer` is the **wet-site** boundary — right for riparian,
  floodplain and wetland columns, and wrong for an upland one, where it will supply water a real
  hillslope would not.
- The old `aquifer` was *identical to free drainage* in the flux it applied: it used the deep-water-table
  limit $`q = K(\theta_n)`$ unconditionally, behind a bucket whose water table could rise into the soil
  column without ever saturating it. Any result predating the rebuild that turned on `aquifer` was
  running free drainage with extra bookkeeping.

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
internal-energy plateau**, captured for free and unconditionally (ice-aware
$`\kappa`$/$`C_{eff}`$; the zero-curtain is tested, cooling a wet layer pins `soil_temp` at the triple
point while `soil_fliq` absorbs the fusion enthalpy). The closed residual is `energy_resid`
$`=\Delta E - \Delta t\,(G_{top}-\text{bottom}-\sum\text{root\_heat\_sink})\approx0`$. A sibling
`soil_energy_time_deriv` exposes the same flux divergence as an explicit RHS (faces at $`T^n`$) for the
ARK integrator.

### The bottom thermal boundary

Two boundary conditions are selectable, through `[energy].bottom_bc`.

**`geothermal` (Neumann, the default).** A prescribed flux at the base, held at exactly zero. The term
is fully plumbed — `forcing%geothermal` reaches the kernel, is differenced into $`hf_n`$ and is debited
to the ledger as `flux%bottom_heat` — but nothing assigns it a non-zero value, so the base is an
**adiabatic wall**. Real continental geothermal flux is ~0.05–0.09 W m⁻², two to three orders below the
diurnal $`G_{top}`$ signal, so neglecting *that* is defensible; the wall itself is the consequential
half. A zero-flux base **reflects** the downward thermal wave instead of transmitting it, and with a
2.0 m column against an annual damping depth near 2 m the annual cycle has barely attenuated by the
time it arrives.

**`dirichlet` (a deep temperature anchor).** The bottom node conducts to a plane held at
`deep_temp`, a distance $`\ell = \texttt{deep\_depth} - |z_{node,n}|`$ below it, through the bottom
layer's own conductivity:

```math
hf_n = -\,g_{deep}\,(T_n - T_{deep}), \qquad g_{deep} = \frac{\kappa_n}{\ell} \qquad(4)
```

The anchor enters the BE **matrix** (it adds $`g_{deep}`$ to the bottom diagonal and
$`g_{deep}T_{deep}`$ to the residual), not the right-hand side alone, so the solve stays
unconditionally stable. The same `bottom_heat_face` expression serves the implicit step at $`T^{n+1}`$
and the explicit `soil_energy_time_deriv` at $`T^n`$, so the two paths cannot drift apart.

#### Why the anchor sits *below* the column, and how deep

Pinning the base **face** to a fixed temperature is not the answer: an adiabatic base reflects the
annual wave with coefficient $`+1`$, and a pinned face reflects it just as hard with coefficient
$`-1`$. Neither is transparent. What a semi-infinite continuation presents to the column at its base
is an impedance $`\kappa(1+i)/d`$, where $`d=\sqrt{2\alpha/\omega}`$ is the annual damping depth. A
resistive link of length $`\ell`$ presents $`\kappa/\ell`$, and the reflection coefficient is smallest
when the two magnitudes match:

```math
\frac{\kappa}{\ell} = \left|\frac{\kappa(1+i)}{d}\right| = \frac{\sqrt{2}\,\kappa}{d}
\quad\Longrightarrow\quad \ell = \frac{d}{\sqrt{2}} \qquad(5)
```

So the anchor **depth** is a physical choice, and the default is derived rather than fitted. For the
default column (2 m, 10 layers, `grid_growth = 3`) in a mid-latitude mineral soil at $`\theta = 0.3`$:
$`d = 1.973`$ m and $`|z_{node,10}| = 1.727`$ m, giving `deep_depth` $`= 1.727 + 1.973/\sqrt{2} = 3.12`$ m.

Measured against the analytic profile (`test_soil_annual_damping`, a homogeneous column driven by a
harmonic $`G_{top}`$ to periodic steady state — the constant-$`\kappa`$, constant-$`C`$ case where
$`\exp(-z/d)`$ holds exactly):

| bottom BC | amplitude at $`-1.73`$ m, relative to the surface | vs analytic | RMS error over the profile |
|---|---|---|---|
| analytic (semi-infinite) | 0.421 | — | — |
| `geothermal` (adiabatic) | 0.764 | **+82 %** | 0.145 |
| `dirichlet`, `deep_depth = 3.12` | 0.412 | **−2 %** | 0.015 |

A sweep of `deep_depth` puts the measured optimum at 3.0–3.1 m — the derivation above is right to
within one sweep step — and the minimum is flat enough that anything from 2.8 to 4.0 m is still three
to seven times better than the adiabatic wall. The harness is validated independently: a 12 m column
with the adiabatic base reproduces $`\exp(-z/d)`$ to 0.1 % over the top three damping depths, because
there the boundary is far enough away to be irrelevant.

**What this does not fix.** A purely resistive termination cannot reflect less than 0.41 in amplitude,
whatever $`\ell`$ is — matching a complex impedance with a real one leaves the phase wrong by 45°.
Closing that last gap needs heat *capacity* below the column, i.e. real layers: passive deep thermal
layers under the hydrologically active column (`docs/ROADMAP.md`, #145 follow-up). The anchor buys
roughly a factor of ten in this metric, not exactness.

**Choosing `deep_temp` is the user's job, and the loader insists on it.** It is the mean annual soil
temperature below the damping depth — close to the mean annual air temperature of the driving forcing,
and a site property like latitude. An error in it is a steady flux $`g_{deep}\,\Delta T`$ into the
column base, so it biases deep-soil temperature in the annual mean; there is no defensible default and
a silent one would reintroduce, in a new place, exactly the bias this boundary condition exists to
remove.

## Water–thermal coupling

The water and thermal columns are coupled by construction: the thermal step reads the just-updated soil
moisture (for $`\kappa,C_{eff}`$), and the transpiration/infiltration/drainage water carries liquid
enthalpy across the soil boundaries via `root_heat_sink` and the boundary-face advection.

## Canopy interception

Per-cohort **canopy interception** is a separate top→bottom cascade (`intercept_canopy_layer`, in
`meds_plant_biophysics`): a capacity-limited bucket with a Beer interception fraction, each
cohort's throughfall feeding the next, so the net throughfall reaching the soil top is what the
infiltration boundary sees.

## Soil heterotrophic respiration

The heterotrophic soil CO₂ flux $`R_h`$ comes from **one authority**: the CENTURY transfer matrix
(`meds_soil_biogeochem`: `heterotrophic_respiration_matrix`, reached through
`patch_heterotrophic_respiration`). It is a soil-carbon decomposition process (the slow
biogeochemistry lives in `biogeochemistry/`), but its flux enters the fast loop as a source term of
the CAS CO₂ twin
(see [canopy_air_space_biophysics](canopy_air_space_biophysics.md)). The driver is its single authority
and passes the resulting flux into the CO₂ source.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Soil water | `soil_column_t` | `theta(:)`, `w_surface` | `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |

## Code map

| Concept | Routine |
|---|---|
| soil water (implicit Richards) | `meds_soil_water`: `advance_soil_water_column`, `soil_water_step_implicit`, `soil_water_advance` |
| explicit water tendency (ARK) | `meds_soil_water`: `soil_water_time_deriv` |
| retention curves (vG / Campbell) | `meds_hydr_lib` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| explicit thermal tendency (ARK) | `meds_soil_energy`: `soil_energy_time_deriv` |
| canopy interception | `meds_plant_biophysics`: `intercept_canopy_layer` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_flux`, `heterotrophic_respiration_damm` |
| prognostic soil types | `meds_column_state_types`: `soil_column_t`, `soil_energy_column_t` |

## References
- ED2 `../ED2/ED/src/dynamics/`: `lsm_hyd.f90`, `soil_respiration.f90` — the soil stores in the reference model.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — mixed-form Richards / modified-Picard.
- van Genuchten (1980); Clapp & Hornberger (1978) — soil retention curves.
- Swenson & Lawrence (2014), *JGR-Atmos.* 119:10299 — dry-surface-layer soil evaporation.
- Design docs: `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`.
