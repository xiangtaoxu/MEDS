# Column biophysics

The integrative fast-loop (sub-daily) surface model: how MEDS advances the coupled internal-energy,
water, and CO₂ budgets of the whole soil–vegetation–atmosphere column each `dt_fast`. This page ties
together the per-store kernels — each documented on its own page — and describes the two integrators
(`meds_fast_integ_esdirk2` ESDIRK2, `meds_fast_integ_rk45` Cash-Karp, run by `meds_fast_integ_stepper`)
that weave the stores together:

- [canopy_air_space_biophysics](canopy_air_space_biophysics.md) — the CAS coupling reservoir and its
  three prognostic twins (enthalpy, humidity, CO₂).
- [soil_biophysics](soil_biophysics.md) — the coupled soil water (Richards) and soil thermal
  (internal-energy heat) columns, canopy interception, ground evaporation, heterotrophic respiration.
- [soil_plant_water](soil_plant_water.md) — the pond, soil layers and every cohort's wood and leaves,
  solved as one water network each fast step.
- [vegetation_energy_dynamics](vegetation_energy_dynamics.md) — the leaf/wood tissue energy balance
  (diagnostic and prognostic modes).
- [snow_biophysics](snow_biophysics.md) — the always-active snow / temporary-surface-water store.
- [canopy_aerodynamics](canopy_aerodynamics.md) — the Monin-Obukhov surface layer and boundary-layer
  conductances that set the shared $`u_*`$ atm↔CAS transfer.
- [canopy_radiation_transfer](canopy_radiation_transfer.md) — the two-stream absorbed SW / net LW that
  drive the surface energy stores.
- [plant_hydraulics](plant_hydraulics.md) — the soil→plant→stomata water pathway coupled through
  `psi_soil` and the realized transpiration.
- [numerical_scheme](numerical_scheme.md) — *how* the coupled column is advanced in time: what is
  frozen each `dt_fast`, the two integrators and what differs between them, and the
  measured accuracy/cost trade. Read it before choosing `[fast].time_integrator` or comparing runs
  made with different ones.

---

## Two cross-cutting principles

1. **Prognostic internal energy, not temperature.** Every thermal store carries specific enthalpy or
   internal energy as its prognostic variable and diagnoses temperature (and liquid fraction) by
   inverting the thermodynamics (`meds_therm_lib%internal_energy_to_temp`). Freeze/thaw is then a *read-off* of the
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

`fast_step` (`meds_fast_integ_stepper`) advances one patch one `dt_fast` with one of **two**
integrators (`[fast].time_integrator`); [numerical_scheme](numerical_scheme.md) §3 compares them and
says why an operator-split stepper is not among them.

Whichever runs, the step first builds the **frozen record** (`build_frozen`: aerodynamics; leaf gas
exchange → GPP/$g_s$/$R_d$; stem, root and heterotrophic respiration; CAS capacities; the snow stage
and canopy interception) and steps the soil-plant water network (`water_step`), whose time-mean flows
are handed on as constants (`record_water_flows`) and whose θ is committed. After the march,
`update_dynamic_state` commits the step and `close_ledgers` closes the ledgers, the same for both
schemes. The canopy radiation transfer runs before it, in the fast-loop driver
([order_of_processes](order_of_processes.md) §3). Those quantities are **frozen for the whole
`dt_fast`** — the Category-0 semi-discretisation MEDS inherits from ED2.

**One coefficient is deliberately excluded from the freeze.** The CAS↔atmosphere conductances are
re-solved at every integrator stage, because the canopy air is a low-capacity node
(`cas_mass_capacity·cp ≈ 2.4×10⁴ J m⁻² K⁻¹`) whose own ventilation depends on its own temperature through
atmospheric stability. Freezing that one feedback across a step would feed a lagged canopy-air
temperature back into its own balance and drive a sustained period-2 oscillation (~8 K peak-to-peak
at 900 s) that **no conservation ledger detects**. See numerical_scheme §2 and canopy_aerodynamics §2.

**What bounds `dt_fast` is accuracy, not stability.** Leaf gas exchange and the leaf water potential
that sets stomatal conductance are frozen per step, so `dt_fast` sets how far the carbon and water
fluxes lag the state they respond to. The default 900 s suits general production, spin-ups included;
use ≤ 150 s for sub-daily fidelity (flux-tower comparison, diel cycles, energy partitioning).
numerical_scheme §5a gives the measured errors and §6 the choice.

### `esdirk2` (default) — 2-solve ESDIRK2

`esdirk2_step` (`meds_fast_integ_esdirk2`) advances the coupled column with a two-stage L-stable implicit
stepper (γ = 1 − 1/√2, the ARS(2,2,2) value). Each stage is a backward-Euler solve
(`column_be_stage`, `meds_fast_integ_backward_euler`) that calls the stores' own implicit kernels —
`cas_box_step_implicit` for the CAS twins, `soil_energy_step_implicit` for the soil heat — over
the side-effect-free surface tendencies of `meds_fast_integ_derivs` (`surface_derivs`). Within every
stage the leaf↔CAS block is solved implicitly by a direct 2×2 Newton (`newton_surface_solve`). Soil
water and plant water move at the water network's rates, and canopy surface water is operator-split
out of the tableau; each is advanced once over the full step. The **conservation ledgers**
accumulate the $b$-weighted boundary-flux amounts.

The ARS(2,2,2) tableau is usually presented as an IMEX (additive Runge–Kutta) pair, but this scheme is
**not IMEX**: the biotic CO₂ source is folded implicit, so the explicit tableau is empty (`f_E == 0`).

### `rk45` — adaptive Cash–Karp 5(4), the explicit alternative

A full alternative a run can choose, and the explicit reference the tests compare `esdirk2` against.
Fully explicit over the canopy air and the column's energy, through the whole-column tendency
`march_derivs` (`meds_fast_integ_derivs`); soil and plant water move at the water network's rates,
as on `esdirk2`. When a step cannot be resolved (stiff cold-canopy nights) the stepper discards the
march and marches that `dt_fast` again on `esdirk2`, from the same frozen record; the counter `mean_work_rk45_rescue_py` reports how often, and it must be
read before any RK45 result.

Whichever integrator runs, **every store closes its residual** and the whole-column ledger closes to
machine precision — the invariant that makes the coupled surface model trustworthy.

---

## Prognostic state (per patch, `fast_dynamic_state_t` / `meds_patch_state_vars`)

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` (`can_depth` is set by the slow step) | `can_temp` |
| Soil water | `soil_water_column_t` | `theta(:)`, `w_surface`, `w_surface_energy` [J m⁻²] | `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |
| Snow / surface water | `snow_column_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |
| Leaf / wood tissue | `fast_dynamic_state_t` | `leaf_energy`, `wood_energy` [J m⁻²; on the cohort block, J plant⁻¹] (exact exponential relaxation) | `leaf_temp`, `wood_temp` |
| Plant water | `fast_dynamic_state_t` | `leaf_water_mass`, `wood_water_mass` [kg plant⁻¹] | leaf and wood ψ |
| Canopy surface water | `fast_dynamic_state_t` | `leaf_surf_water`, `wood_surf_water` [kg m⁻²; on the cohort block, kg plant⁻¹] | — |

`fast_dynamic_state_t` (`meds_fast_integ_types`) is the fast loop's working copy of one patch: the fast
driver (`meds_fast_state`: `gather_patch`) gathers it from the persisted patch and cohort state and
writes it back after the slow step's sub-steps (`scatter_patch`). Each integrator packs it into a
flat `march_state_t`, which it can combine and measure as one vector. What the sub-steps only read
is gathered beside it, once a day, into `fast_static_state_t`: each cohort's structure and traits,
the canopy geometry and its canopy-air top, the site's acclimated leaf table, and three inputs from
the slow tier: the soil-carbon pools that heterotrophic respiration draws on
([soil_carbon](soil_carbon.md) §7), the water and enthalpy shed by leaf and root turnover, and the
slow tier's CO₂ owing from growth respiration ([plant_respiration](plant_respiration.md) §6).
`fast_step(dt_fast, cfg, site, static, forcing, dynamic, out, work)` reads the run's configuration,
the patch's site (its soil column) and the step's forcing, advances `dynamic`, and fills `out`: the
aerodynamics, each cohort's carbon fluxes and leaf potential, H and LE, and the step report
(`step_report_t`: the ledgers, the exports and the work counters).

### Mixing two patches

When two patches fuse, or a disturbance carves a gap from its donors, each store of the new patch is
mixed from the two by a `blend_*` function (`meds_patch_state_vars`). The weights are the area
fractions $`w_i = a_i/(a_1+a_2)`$:

- **Per ground area: area-weighted.** Soil water and soil energy per layer, the pond's mass and
  energy, the snow's mass and energy, the soil-carbon pools, and the day's accumulators (`decomp_xi`,
  `litter_in`) are all amounts per m² of ground, so $`w_1 v_1 + w_2 v_2`$ keeps the site total.
- **The canopy air: weighted by air mass.** Its enthalpy, humidity and CO₂ are per kg of air, and a
  patch holds $`\rho\,a\,d`$ of air over depth $d$. The merged depth is area-weighted (the merged volume
  over the merged ground), so conserving the content gives

  ```math
  v = \frac{a_1 d_1 v_1 + a_2 d_2 v_2}{a_1 d_1 + a_2 d_2} \qquad(1)
  ```

  where $\rho$, uniform over the site, cancels. Weighting on area alone would drop the covariance
  term $`w_1 w_2 (d_1-d_2)(v_1-v_2)`$: fusing a tall patch with a gap would shift the canopy air's
  energy, humidity and CO₂ by the same relative error, of order 10⁴ J per fusion.
- **The cohorts: only their density.** A cohort's amounts are per plant (its carbon, tissue water and
  energy, and canopy films), so the cohorts of the merged patch keep them, and their densities take
  the area weights.
- **Diagnosed fields** (temperatures, liquid fractions) are mixed on the same weights only so that a
  read before the caller re-diagnoses them from the blended energy is not far off.

## References
- ED2 `../ED2/ED/src/dynamics/`: `rk4_derivs.f90`, `canopy_struct_dynamics.f90`, `soil_respiration.f90`,
  `lsm_hyd.f90` — the fast-loop stores in the reference model.
- Kennedy & Carpenter (2003) — additive Runge-Kutta (IMEX) methods.
- Ascher, Ruuth & Spiteri (1997), *Appl. Numer. Math.* 25:151 — the ARS(2,2,2) tableau of `esdirk2`.
- Per-store references (retention curves, DSL evaporation, snow-cover fraction) are on the domain pages
  linked above.
- Design docs (archived, `docs/dev_plans/archive/`): `MEDS_COLUMN_DYNAMICS_DESIGN.md`,
  `MEDS_ENERGY_BALANCE_DESIGN.md`, `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`,
  `MEDS_IMEX_ARK_DESIGN.md`, `MEDS_SNOW_DESIGN.md`, `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`,
  `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.

## Code map (whole column)

| Concept | Routine |
|---|---|
| CAS enthalpy + vapour + CO₂ twins | `meds_cas_biophysics`: `cas_box_step_implicit`, `cas_box_time_deriv` |
| soil and plant water (one network) | `meds_fast_integ_backward_euler`: `water_step`; `meds_fast_integ_derivs`: `water_derivs`, `set_water_state`; `meds_fast_integ_frozen`: `record_water_flows` |
| canopy interception | `meds_plant_biophysics`: `intercept_canopy_layer` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| leaf/wood energy | `meds_plant_biophysics`: `veg_energy_balance` — one kernel for leaf and wood, diagnostic (`store_hcap_per_dt = 0`) or prognostic (`> 0`) |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_matrix`; per patch and sub-step, `meds_fast_integ_frozen`: `patch_heterotrophic_respiration` |
| snow energy / base conductance | `meds_ground_biophysics`: `snow_energy_step`, `snow_base_conductance` |
| snow mass / cover / melt | `meds_ground_biophysics`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| the step, the march + RK45→ESDIRK2 rescue | `meds_fast_integ_stepper`: `fast_step`, `adaptive_march` |
| the frozen record | `meds_fast_integ_frozen`: `build_frozen`, `record_water_flows` |
| the update + ledgers | `meds_fast_integ_update`: `update_dynamic_state`, `close_ledgers` |
| ESDIRK2 step (config `esdirk2`) | `meds_fast_integ_esdirk2`: `esdirk2_step` |
| one ESDIRK2 stage | `meds_fast_integ_backward_euler`: `column_be_stage`, `newton_surface_solve` |
| Cash–Karp step (config `rk45`) | `meds_fast_integ_rk45`: `rk45_step`, `rk45_state_railed` |
| whole-column tendency RHS | `meds_fast_integ_derivs`: `march_derivs` (`rk45`), `surface_derivs` (both) |
| explicit tendency siblings | `meds_soil_energy`: `soil_energy_time_deriv`; the soil and plant water tendencies are the water network's rates, set in `march_derivs` |
| prognostic column types | `meds_patch_state_vars`: `cas_state_t`, `snow_column_t`, `soil_water_column_t`, `soil_energy_column_t` |
