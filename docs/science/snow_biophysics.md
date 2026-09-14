# Snow / temporary-surface-water biophysics

An always-active mass+energy reservoir stacked between the canopy air space and the top soil
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
energy [J m⁻²]; temperature and liquid fraction are read-offs of `internal_energy_to_temp` (`dry_hcap=0`), so
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

### Sublimation saturates over ice

The vapour flux is $`E=g_{net}\rho\,[q_{sat}(T_s,f_{liq})-q_{CAS}]`$, and $`q_{sat}`$ is taken on the
**ice** curve when the pack is frozen:

```math
e_{sat}(T_c,f_{liq}) = f_{liq}\cdot 611.2\,e^{\frac{17.67\,T_c}{T_c+243.5}}
                     + (1-f_{liq})\cdot 611.2\,e^{\frac{21.87\,T_c}{T_c+265.5}} \qquad(1)
```

Both branches carry the same 611.2 Pa constant, so they cross **exactly** at $`T_c=0`$ and the blend is
continuous in temperature *and* in $`f_{liq}`$ — a pack that freezes or melts slides between the curves
rather than stepping, which matters because a step here lands in the right-hand side an adaptive
controller integrates.

$`f_{liq}`$ is the pack's own prognostic liquid fraction, not a temperature threshold. The energy side
was already ice-aware — removing `enthalpy_vapor` from an ice-referenced layer debits sublimation
(vaporization + fusion) automatically — so before this the model treated the surface as ice for
*enthalpy* and as liquid for *vapour pressure*. Over ice $`e_{sat}`$ is **10 % lower at −10 °C, 22 % at
−20 °C and 34 % at −30 °C**, so the liquid curve overstated the driving gradient by those factors. The
ice form is within 0.1 % of Murphy & Koop (2005) at −10 °C and 0.9 % at −30 °C.

The same blend is used for ground evaporation from frozen soil, weighted by the top layer's
`soil_fliq`. It is deliberately **not** used for dewpoint conversion or diagnostic VPD: dewpoint is
*defined* over liquid, so an ice branch there would mis-convert the forcing.

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
