# Snow / temporary-surface-water biophysics

An opt-in (`[fast].snow_on`) mass+energy reservoir stacked between the canopy air space and the top soil
layer. When present it becomes the CAS-facing surface (its energy balance replaces the bare-ground skin)
and the soil column's top thermal boundary, blended by a cover fraction so the transition through
snow-on/snow-off is continuous.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven each `dt_fast`.

---

## The snow store

The `snow_*` kernels live in `meds_ground_biophysics` (`docs/dev_plans/MEDS_SNOW_DESIGN.md` P0 — a single
bulk layer). The prognostic state is water-equivalent mass `swe` [kg m⁻²] and **extensive** internal
energy [J m⁻²]; temperature and liquid fraction are read-offs of `uext_to_temp` (`dry_hcap=0`), so
**melt/refreeze is the internal-energy plateau**, exactly as for soil. The fast-loop driver orchestrates
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

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Snow / surface water | `snow_column_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |

## Code map

| Concept | Routine |
|---|---|
| snow energy / base conductance | `meds_ground_biophysics`: `snow_energy_step`, `snow_base_conductance` |
| snow mass / cover / melt | `meds_ground_biophysics`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| prognostic snow type | `meds_column_state_types`: `snow_column_t` |

## References
- Niu & Yang (2007), *JGR* 112:D21101 — snow-cover fraction.
- Design doc: `MEDS_SNOW_DESIGN.md`.
