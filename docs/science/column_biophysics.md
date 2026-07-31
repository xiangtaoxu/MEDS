# Column biophysics

The integrative fast-loop (sub-daily) surface model: how MEDS advances the coupled internal-energy,
water, and CO₂ budgets of the whole soil–vegetation–atmosphere column each `dt_fast`. This page ties
together the per-store kernels — each documented on its own page — and describes the two integrators
(`meds_fast_step` dispatch, `meds_fast_ark` ESDIRK2, `meds_fast_rk45` Cash-Karp) that weave the stores together:

- [canopy_air_space_biophysics](canopy_air_space_biophysics.md) — the CAS coupling reservoir and its
  three prognostic twins (enthalpy, humidity, CO₂).
- [soil_biophysics](soil_biophysics.md) — the coupled soil water (Richards) and soil thermal
  (internal-energy heat) columns, canopy interception, ground evaporation, heterotrophic respiration.
- [vegetation_energy_dynamics](vegetation_energy_dynamics.md) — the leaf/wood tissue energy balance
  (diagnostic and prognostic modes).
- [snow_biophysics](snow_biophysics.md) — the opt-in snow / temporary-surface-water store.
- [canopy_aerodynamics](canopy_aerodynamics.md) — the Monin-Obukhov surface layer and boundary-layer
  conductances that set the shared $`u_*`$ atm↔CAS transfer.
- [canopy_radiation_transfer](canopy_radiation_transfer.md) — the two-stream absorbed SW / net LW that
  drive the surface energy stores.
- [plant_hydraulics](plant_hydraulics.md) — the soil→plant→stomata water pathway coupled through
  `psi_soil` and the realized transpiration.
- [numerical_scheme](numerical_scheme.md) — *how* the coupled column is advanced in time: what is
  frozen each `dt_fast`, the three selectable integrators and what differs between them, and the
  measured accuracy/cost trade. Read it before choosing `[fast].time_integrator` or comparing runs
  made with different ones.

---

## Two cross-cutting principles

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
aerodynamics kernel. See [canopy_air_space_biophysics](canopy_air_space_biophysics.md).

---

## Weaving the stores: how the column is advanced

`column_fast_step` (`meds_fast_step`) advances one patch one `dt_fast` and dispatches to one of **two**
integrators (`cfg%time_integrator`). A third, the operator-split + Picard stepper, was **retired**
(2026-07-31) — see [numerical_scheme](numerical_scheme.md) §3 for why, and note that the
`integration_scheme` selector went with it.

Whichever runs, the driver first executes a **pre-pass** (leaf gas exchange → GPP/$g_s$/$R_d$;
stem+root respiration; CAS capacities; aerodynamics; canopy RT; the plant-hydraulics matrix-exponential
solve, whose time-averaged sapflow and root uptake are handed on as constants). Those quantities are
**frozen for the whole `dt_fast`** — the Category-0 semi-discretisation MEDS inherits from ED2.

**That freeze bounds `dt_fast` from above by STABILITY, not accuracy.** The canopy air is a
low-capacity node (`wcap·cp ≈ 2.4×10⁴ J m⁻² K⁻¹`) under fluxes of hundreds of W m⁻², so holding its
coupling coefficients fixed across a long step feeds a lagged canopy-air temperature back into its own
balance and produces a sustained period-2 oscillation (~8 K peak-to-peak at 900 s, smooth at 100 s).
**No conservation ledger detects it.** Default `dt_fast` is 150 s; see numerical_scheme §2.

### `ark` (default) — 2-solve ESDIRK2

`column_fast_step_ark` (`meds_fast_ark`) drives the coupled column with a two-stage L-stable implicit
stepper (γ = 1 − 1/√2), which needs a side-effect-free $`f(y)=dy/dt`$ for the whole column — supplied by
`meds_fast_time_derivs` (`column_derivs`). Each reservoir's tendency is the explicit RHS whose BE
advance reproduces its kernel as $`\Delta t\to0`$ (soil heat/water) or exactly (the CAS
implicit-in-atm twins; the plant-hydraulics $2\times2$ matrix-exponential). Within every stage the
leaf↔CAS block is solved implicitly by a direct 2×2 Newton (`newton_surface_solve`). Soil water, plant
water mass and canopy surface water are operator-split out of the tableau and advanced once over the
full step. An ARK **conservation ledger** accumulates the $b$-weighted boundary-flux amounts.

Despite the historical name this is **not IMEX**: the biotic CO₂ source is folded implicit, so the
explicit tableau is empty (`f_E == 0`).

### `rk45` — adaptive Cash–Karp 5(4), the accuracy baseline

Fully explicit over the whole column state, including soil moisture. When a step cannot be resolved
(stiff cold-canopy nights) the driver discards it and redoes that `dt_fast` on `ark`; the counter
`work_rk45_rescue_site` reports how often, and it must be read before any RK45 result.

Whichever integrator runs, **every store closes its residual** and the whole-column ledger closes to
machine precision — the invariant that makes the coupled surface model trustworthy.

---

## Prognostic state (per patch, `patch_biophys_t` / `meds_column_state_types`)

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2`, `can_depth` | `can_temp` |
| Soil water | `soil_column_t` | `theta(:)`, `w_surface` | `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |
| Snow / surface water | `snow_column_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` (heat store, exact exponential relaxation) | — |
| Plant hydraulics | `patch_biophys_t%psi` | node water potentials (leaf/wood/root) | — |

Area-weighted `blend_*` mixers conserve these intensively when patches fuse or a disturbance gap is
carved (temp/fliq are re-diagnosed, never blended).

## References
- ED2 `../ED2/ED/src/dynamics/`: `rk4_derivs.f90`, `canopy_struct_dynamics.f90`, `soil_respiration.f90`,
  `lsm_hyd.f90` — the fast-loop stores in the reference model.
- Kennedy & Carpenter (2003) — additive Runge-Kutta (IMEX) methods.
- Per-store references (retention curves, DSL evaporation, snow-cover fraction) are on the domain pages
  linked above.
- Design docs: `MEDS_COLUMN_DYNAMICS_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`,
  `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`, `MEDS_IMEX_ARK_DESIGN.md`,
  `MEDS_SNOW_DESIGN.md`, `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.

## Code map (whole column)

| Concept | Routine |
|---|---|
| CAS enthalpy + vapour + CO₂ twins | `meds_cas_biophysics`: `cas_column_step_implicit`, `cas_column_time_deriv` |
| soil water (implicit Richards) | `meds_soil_water`: `column_hydrology_flux`, `soil_water_step_implicit`, `soil_water_advance` |
| canopy interception | `meds_vegetation_biophysics`: `intercept_canopy_layer` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| leaf/wood energy | `meds_vegetation_biophysics`: `veg_energy_diagnostic` (shared diagnostic solve), `veg_energy_step_implicit` (prognostic store) |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_flux`, `heterotrophic_respiration_damm` |
| snow energy / base conductance | `meds_ground_biophysics`: `snow_energy_step`, `snow_base_conductance` |
| snow mass / cover / melt | `meds_ground_biophysics`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| dispatch + RK45→ARK rescue | `meds_fast_step`: `column_fast_step` |
| ESDIRK2 step (config `ark`) | `meds_fast_ark`: `column_fast_step_ark`, `ark2_column_step`, `adaptive_ark_march` |
| whole-column tendency RHS | `meds_fast_time_derivs`: `column_derivs`, `surface_derivs` |
| explicit tendency siblings | `soil_energy_time_deriv`, `soil_water_time_deriv`, `plant_water_tendency` |
| prognostic column types | `meds_column_state_types`: `cas_state_t`, `soil_column_t`, `soil_energy_column_t`, `snow_column_t` |
