# Column biophysics

The integrative fast-loop (sub-daily) surface model: how MEDS advances the coupled internal-energy,
water, and CO₂ budgets of the whole soil–vegetation–atmosphere column each `dt_fast`. This page ties
together the per-store kernels (energy, hydrology, snow, CO₂) with the aerodynamics
([canopy_aerodynamics](canopy_aerodynamics.md)) and radiation
([canopy_radiation_transfer](canopy_radiation_transfer.md)) pages, and describes the two integrators
(`meds_column_dynamics`) that weave the stores together. Two design principles run through everything:

1. **Prognostic internal energy, not temperature.** Every thermal store carries specific enthalpy or
   internal energy as its prognostic variable and diagnoses temperature (and liquid fraction) by
   inverting the thermodynamics (`meds_therm_lib%uext_to_temp`). Freeze/thaw is then a *read-off* of the
   inverter — cooling a wet layer pins its temperature at the triple point while the internal energy
   keeps falling and the liquid fraction absorbs $`w\,L_f`$. Phase change comes for free, with zero
   solver change.
2. **Closed-budget discipline.** Every store returns a residual that is the change in storage minus the
   time-integral of its boundary fluxes, algebraically zero by construction and asserted to round-off in
   Debug builds. A whole-column ledger (`budg%whole_energy`/`whole_water`) catches cross-seam leaks the
   per-kernel budgets would miss.

The **canopy air space (CAS)** is the shared coupling reservoir. It carries three implicit prognostic
twins — enthalpy, specific humidity, CO₂ mixing ratio — each advanced from its surface sources and an
implicit atmospheric-exchange term that shares the $`u_*\cdot\mathrm{temp1/temp2}`$ conductance from the
aerodynamics kernel.

---

## 1. The three CAS twins

The CAS is a well-mixed box of air of depth $D_{can}$ (from canopy height) and mass per ground area
$`W_{cap}=\rho\,D_{can}`$. Its enthalpy and humidity twins are advanced by `canopy_air_update`
(`meds_column_energy`); the CO₂ twin by `canopy_air_co2_update` (`meds_column_co2`). Each twin obeys

```math
\mathcal{C}\,\frac{dX}{dt} = F_{surf} + g_{atm}\,(X_{atm}-X) \qquad(1)
```

with capacity $\mathcal{C}$, summed surface source $F_{surf}$, and the atm↔CAS conductance $g_{atm}$.
The atmospheric term is taken **implicit** for L-stability (it is the stiff coupling); the surface
sources are explicit. For the enthalpy twin, with $`F_{sens}`$ the summed cohort+ground sensible+latent
flux and $`g_{atm}=\rho\,u_*\,\mathrm{temp1}`$:

```math
H^{n+1} = \frac{W_{cap}H^n + \Delta t\,(F_{sens} + g_{atm}\,H_{atm})}{W_{cap} + \Delta t\,g_{atm}},
\qquad
\text{resid} = W_{cap}(H^{n+1}-H^n) - \Delta t\big(F_{sens}+g_{atm}(H_{atm}-H^{n+1})\big) = 0 \qquad(2)
```

The **vapour** twin is identical with $`g_{aw}=\rho u_*\,\mathrm{temp2}`$ and a non-negativity clamp.
The **CO₂** twin differs only in units: CO₂ is a molar mixing ratio, so the capacity is the dry-air
**molar** column $`C_{cap}=\rho\,(1-q)/M_{d}\cdot D_{can}`$ [mol m⁻²] and the conductance is
$`g_{ac}=\rho_{mol}\,u_*\,\mathrm{temp2}`$. Its biotic source is
$`F_{bio}=R_a+R_h-\mathrm{GPP}`$ (net ecosystem exchange), assembled by `column_co2_step` from the
aggregated cohort GPP/leaf-respiration ($`\times`$ LAI), stem/root maintenance respiration
($`\times n_{plant}`$), and the heterotrophic soil flux (§5). CAS temperature is re-diagnosed from
$(H^{n+1},q^{n+1})$ each step. Sharing $`u_*`$ **and** the profile factor across all three twins keeps
them on one turbulence basis — the fix for the nocturnal-CO₂ over-coupling bug (§canopy_aerodynamics).

---

## 2. The soil water column

`column_hydrology_flux` (`meds_column_hydrology`) advances one patch's prognostic soil moisture
$`\theta_k`$ (+ ponded surface water and a lumped aquifer) over `dt_fast` with an **implicit
backward-Euler Thomas** solve of the mixed-form Richards equation on the ED2 negative-$z$ grid. Per
layer,

```math
C(\psi_k)\,\frac{\partial\psi_k}{\partial t}
= \frac{\partial}{\partial z}\!\left[K(\theta)\Big(\frac{\partial\psi}{\partial z}+1\Big)\right] - S_k
\qquad(3)
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
 - \Delta t\,(P - E_{soil} - \text{drainage} - \text{uptake} - \text{runoff}) \approx 0 \qquad(4)
```

and exports the per-layer matric potential `psi_soil` [MPa] that closes the **plant-hydraulics** soil
boundary condition (see [plant_hydraulics](plant_hydraulics.md)), plus the time-mean per-face Darcy
flux for optional advective soil heat. Per-cohort **canopy interception** is a separate top→bottom
cascade (`intercept_canopy_layer`): a capacity-limited bucket with a Beer interception fraction, each
cohort's throughfall feeding the next.

---

## 3. The soil thermal column

`soil_energy_flux` (`meds_column_energy`) advances the prognostic volumetric internal energy
$`E_k`$ [J m⁻³] by an implicit BE heat-diffusion solve that **reuses the same Thomas sweep and negative-$z$
geometry** as the hydrology. At each state it inverts $`E_k\to(T_k,\text{fliq}_k)`$, forms the ice-aware
thermal conductivity $`\kappa(\theta,\text{fliq})`$ and effective volumetric heat capacity
$`C_{eff}(\theta,\text{fliq})`$, solves for $`T^{n+1}`$ with a top Neumann flux $`G_{top}`$ and a bottom
geothermal flux, then commits a **conservative** energy update from the $`T^{n+1}`$ conductive faces plus
optional upwind water-enthalpy advection:

```math
E_k^{n+1} = E_k^n + \frac{\Delta t}{\Delta z_k}\big[(hf_k-hf_{k-1}) + (qwf_k-qwf_{k-1})\big]
            + \Delta t\,q_{src,k} \qquad(5)
```

Temperature and liquid fraction are re-diagnosed from the committed energy — so **freeze/thaw is the
internal-energy plateau**, captured for free once `phase_change = ENERGY_PHASE_ON` (ice-aware
$`\kappa`$/$`C_{eff}`$; the zero-curtain is tested, cooling a wet layer pins `soil_temp` at the triple
point while `soil_fliq` absorbs the fusion enthalpy). The closed residual is `energy_resid`
$`=\Delta E - \Delta t\,(G_{top}-\text{bottom}-\sum\text{root\_heat\_sink})\approx0`$. A sibling
`soil_energy_tendency` exposes the same flux divergence as an explicit RHS (faces at $`T^n`$) for the
ARK integrator.

The water and thermal columns are coupled by construction: the thermal step reads the just-updated soil
moisture (for $`\kappa,C_{eff}`$), and the transpiration/infiltration/drainage water carries liquid
enthalpy across the soil boundaries via `root_heat_sink` and the boundary-face advection.

---

## 4. Leaf and wood energy balance

`veg_energy_balance` (`meds_column_energy`) advances a cohort's leaf **or** wood prognostic internal
energy over one **L-stable linearized backward-Euler** step. It diagnoses $`T^n`$, forms the net flux
$`R^n = A_{sw}+A_{lw} - H - Q_w - Q_{transp}`$ (absorbed SW + net LW − sensible − film-evaporation −
transpiration), and the linearization slope $`\partial R/\partial T = \frac{dR}{dT}\le0`$ (radiative +
sensible + latent-saturation terms, all $\le0$). The implicit temperature estimate stays within one
Newton step of $`T^n`$, so the energy update is **bounded** even when the store is stiff:

```math
T^* = T^n + \frac{R^n}{C/\Delta t - dR/dT}, \qquad
\Delta E = C\,(T^*-T^n) \qquad(6)
```

Using the bounded $`C(T^*-T^n)`$ rather than the endpoint flux $`\Delta t\,R^*`$ is what makes it
L-stable: applying $`\Delta t\,R^*`$ over the whole step over-removes energy from a stiff store and can
drive $T$ negative (a young-stand wood cohort has time constant $`\tau=C/|dR/dT|\sim35`$ s vs
$`\Delta t=1800`$ s). Leaf uses a $2\times$LAI flat-plate sensible area; wood uses $`\pi\times`$WAI
cylinders and no transpiration. The flux report is step-averaged and split by the shares at $`T^*`$, and
`energy_resid` closes to zero.

Two **modes** (config-selectable, `docs/dev_plans/MEDS_LEAF_WOOD_ENERGY_DESIGN.md`): the default
**diagnostic** leaf/wood takes a linearized steady-state balance ($`T_l=T_{cas}+\Delta T_l`$, the leaf
inertia is a negligible ~0.015 K), while the **prognostic** mode advances the store via (6). A
prognostic leaf **requires** `integration_scheme="picard"` — the explicit leaf↔CAS split oscillates
(the leaf is stiff, its sensible+latent feedback on the CAS has a fixed-point slope near $-1$, giving
~1.7 K midday spikes; the Picard iterate damps them to ~0.2 K).

---

## 5. The snow / temporary-surface-water store

An opt-in (`[fast].snow_on`) mass+energy reservoir stacked between the CAS and the top soil layer
(`meds_snow_energy` + `meds_snow_mass`, `docs/dev_plans/MEDS_SNOW_DESIGN.md` P0 — a single bulk layer).
The prognostic state is water-equivalent mass `swe` [kg m⁻²] and **extensive** internal energy [J m⁻²];
temperature and liquid fraction are read-offs of `uext_to_temp` (`dry_hcap=0`), so **melt/refreeze is the
internal-energy plateau**, exactly as for soil. The fast-loop driver orchestrates
**accumulate → energy step → drain**:

- **Accumulate** (`snow_accumulate`): snowfall lands as ice at $`\min(T_{3ple},T_{air})`$, rain-on-snow
  as liquid at $`T_{air}`$ (refreezing later via the inverter). A layer is created only above
  `min_new_snow_mass`; sub-threshold snow is folded into the soil store by the caller.
- **Energy step** (`snow_energy_step`): the same bounded, L-stable, plateau-aware linearized BE step as
  the veg store — net SW (snow albedo) + net LW − sensible − sublimation/evaporation − base conduction,
  with the emission response made consistent with the linearization slope (the wood-store lesson). The
  **snow-base → soil-top conductance** is the series resistance of the half snow layer and the top soil
  node; $`k_{snow}\ll k_{soil}`$ throttles it as the pack deepens, the physical decoupling that caps the
  winter soil surface, and it becomes the soil's top BC.
- **Drain** (`snow_drain_meltwater`): free liquid above the holding capacity percolates out as a paired
  (mass, enthalpy) hand-off to soil infiltration; full melt-out dumps the residual and reverts to bare
  ground.

The **Niu-Yang (2007) snow-cover fraction** $`\mathrm{snowfac}=\tanh(\text{depth}/\text{scale})`$
(`snow_cover_fraction`) ramps the ground optics albedo, the soil-top-BC blend, and the aerodynamic
roughness — and area-weights ("sub-column") all boundary exchange so a thin patchy pack barely exchanges
(stable and continuous through $`\mathrm{snowfac}\to0`$). Snow-off is bit-identical to no store.

---

## 6. Weaving the stores: operator split vs IMEX-ARK

`column_fast_step` (`meds_column_dynamics`) advances one patch one `dt_fast`. It offers two integrators
(`cfg%time_integrator`).

### 6.1 Operator split (default) with leaf↔CAS Picard coupling

Per sub-step the driver runs a **pre-pass** (leaf gas exchange → GPP/$g_s$/$R_d$; stem+root respiration;
CAS capacities; aerodynamics), then a sequence: aerodynamics → {diagnostic leaf balance, ground balance}
→ soil-water column (infiltration / DSL evaporation / root uptake / drainage → `src_frac`, the
supply-limited fraction of transpiration demand) → CAS three-twin update (implicit atm, eq. 2) → soil
thermal column (implicit BE, reading the just-updated moisture). Plant hydraulics then advances each
cohort's node potentials from the realized transpiration and the root-weighted `psi_soil`, feeding
*next* step's leaf gas exchange (the soil→plant→stomata drought feedback, lagged one `dt_fast`).

Under `integration_scheme="picard"` (`SCHEME_PICARD_COUPLED`) the sequence
{leaf energy → soil-water/`src_frac` → CAS twins → soil thermal} is wrapped in an **outer Picard fixed
point**. State$`^n`$ is snapshotted once; each pass re-solves the *same* backward-Euler steps from the
snapshot, under-relaxing the CAS seed (`picard_relax`, ~0.5, to tame the slope-$-1$ oscillation) while
the committed BE solutions stay exact. `niter=1` reproduces the pure operator split **bit-identically**
(`SCHEME_SPLIT_SEQUENTIAL`), so the coupled path is a strict generalization. Convergence tests the
inter-iterate CAS+leaf+ground temperature and CAS humidity; a non-converged run clamps to the last
iterate (never a partial state).

### 6.2 IMEX-ARK (opt-in)

`column_fast_step_ark` drives the coupled column with an additive Runge-Kutta stepper
(`meds_ark_stepper`, `docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md`), which needs a side-effect-free
$`f(y)=dy/dt`$ for the whole column — supplied by `meds_column_derivs` (`column_derivs`). Each
reservoir's tendency is the explicit RHS whose BE advance reproduces its split kernel as
$`\Delta t\to0`$ (soil heat/water, via `soil_energy_tendency`/`soil_water_tendency`) or exactly (the CAS
implicit-in-atm twins; the plant-hydraulics $2\times2$ matrix-exponential). The surface hydrology BCs
(`q_top`, `psi_e`, `soil_psi_root`) and the leaf-gas-exchange pre-pass are the **frozen, explicit** part
of the additive split, held constant across the macro-step; soil water is currently operator-split out
(the robust ponding/runoff scratch solve commits `theta1`, a lagged split). An ARK **conservation
ledger** accumulates the $b$-weighted boundary-flux amounts so the same budgets the split closes also
close on the ARK path.

Whichever integrator runs, **every store closes its residual** and the whole-column ledger closes to
machine precision — the invariant that makes the coupled surface model trustworthy.

---

## Prognostic state (per patch, `patch_biophys_t` / `meds_column_state_types`)

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` | `can_temp` |
| Soil water | `soil_column_t` | `theta(:)`, `w_surface`, `w_aquifer` | `z_wt`, `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |
| Snow / surface water | `snow_column_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` (internal energy) | temp / fliq read-offs |
| Plant hydraulics | `patch_biophys_t%psi` | node water potentials (leaf/wood/root) | — |

Area-weighted `blend_*` mixers conserve these intensively when patches fuse or a disturbance gap is
carved (temp/fliq are re-diagnosed, never blended).

## References
- ED2 `../ED2/ED/src/dynamics/`: `rk4_derivs.f90`, `canopy_struct_dynamics.f90`, `soil_respiration.f90`,
  `lsm_hyd.f90` — the fast-loop stores in the reference model.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — mixed-form Richards / modified-Picard.
- van Genuchten (1980); Clapp & Hornberger (1978) — soil retention curves.
- Swenson & Lawrence (2014), *JGR-Atmos.* 119:10299 — dry-surface-layer soil evaporation.
- Niu & Yang (2007), *JGR* 112:D21101 — snow-cover fraction.
- Kennedy & Carpenter (2003) — additive Runge-Kutta (IMEX) methods.
- Design docs: `MEDS_COLUMN_DYNAMICS_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`,
  `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`, `MEDS_IMEX_ARK_DESIGN.md`,
  `MEDS_SNOW_DESIGN.md`, `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| CAS enthalpy + vapour twins | `meds_column_energy`: `canopy_air_update` |
| CAS CO₂ twin (molar) | `meds_column_co2`: `canopy_air_co2_update`, `column_co2_step` |
| soil water (implicit Richards) | `meds_column_hydrology`: `column_hydrology_flux`, `soil_be_single_step`, `soil_water_advance` |
| canopy interception | `meds_column_hydrology`: `intercept_canopy_layer` |
| soil thermal (implicit BE heat) | `meds_column_energy`: `soil_energy_flux`, `soil_heat_be_step` |
| leaf/wood energy | `meds_column_energy`: `veg_energy_balance`, `veg_surface_fluxes` |
| ground skin balance | `meds_column_energy`: `ground_surface_balance` |
| soil heterotrophic Rh | `meds_column_co2`: `heterotrophic_respiration_flux`, `heterotrophic_respiration_damm` |
| snow energy / base conductance | `meds_snow_energy`: `snow_energy_step`, `snow_base_conductance` |
| snow mass / cover / melt | `meds_snow_mass`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| operator-split + Picard step | `meds_column_dynamics`: `column_fast_step` |
| IMEX-ARK step | `meds_column_dynamics`: `column_fast_step_ark`; `meds_ark_stepper` |
| whole-column tendency RHS | `meds_column_derivs`: `column_derivs`, `surface_derivs` |
| explicit tendency siblings | `soil_energy_tendency`, `soil_water_tendency`, `plant_water_tendency` |
| prognostic column types | `meds_column_state_types`: `cas_state_t`, `soil_column_t`, `soil_energy_column_t`, `snow_column_t` |
