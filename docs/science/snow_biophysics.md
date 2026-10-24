# Snow biophysics

An always-active mass+energy reservoir stacked between the canopy air space and the top soil
layer. When present it becomes the CAS-facing surface (its energy balance replaces the bare-ground skin)
and the soil column's top thermal boundary, blended by a cover fraction so the transition through
snow-on/snow-off is continuous.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven each `dt_fast`.

---

## The snow store

The `snow_*` kernels live in `meds_surface`; the pack is a single bulk layer (see
`docs/ROADMAP.md` for multi-layer snow). The prognostic state is water-equivalent mass `swe`
[kg m⁻²] and **extensive** internal energy [J m⁻²]; temperature and liquid fraction are read-offs of
`internal_energy_to_temp` (`dry_hcap=0`), so **melt/refreeze is the internal-energy plateau**, exactly
as for soil. Each `dt_fast` runs **accumulate → drain** before the march, the **exchange** inside it,
and the **commit** after it:

- **Accumulate** (`snow_accumulate`): snowfall lands as ice at $`\min(T_{3ple},T_{air})`$, rain-on-snow
  as liquid at $`T_{air}`$ (refreezing later via the inverter). A layer is created only from
  `min_new_snow_mass` of new snow in one `dt_fast`; sub-threshold snow joins the rain as liquid (the
  canopy films, then the pond).
- **Drain** (`snow_drain_meltwater`): the free liquid the pack then holds above its holding capacity, the
  last step's melt and this step's rain, percolates out as a paired (mass, enthalpy) hand-off to the
  pond, which it reaches at the melt temperature and infiltrates from in this step's water network;
  full melt-out dumps the residual and reverts to bare ground. Meltwater therefore leaves the pack
  the step after it melts, 15 minutes at the default `dt_fast`.
- **Exchange** (`snow_surface_balance`, called through `ground_surface_exchange` at every evaluation of the march):
  a linearized backward-Euler balance of the pack's enthalpy over the step, against the canopy air and
  the top soil node the march holds at that evaluation, as each tissue's is against them — the
  shortwave it absorbs at its own albedo + the longwave it absorbs less its own emission − sensible −
  sublimation/evaporation − base conduction, each linearized about
  the pack's temperature $`T^n`$, so the net flux $R(T)$ is linear and falls as $T$ rises. The pack's
  enthalpy $`E(T)`$ rises with $T$ along the ice line, across the melt plateau and along the liquid
  line, so the balance $`E(T) = E^n + \Delta t\,f_{snow}\,R(T)`$ has one root, found in closed form on
  the plateau when its enthalpy range holds it, else on the ice or the liquid line, and the fluxes are
  taken at that root: a thin pack under an open patch's bare-ground conductance relaxes toward its
  balance and cannot overshoot it, and a melting pack under cold air freezes through the plateau and
  stops where the balance does. Sublimation takes its latent heat from the balance, and the mass that
  leaves takes the pack's specific energy with it, so the rest keeps $T$; a pack that would sublimate
  whole does, the vapour taking all the energy left. Because the canopy air's temperature and humidity
  are the march's own at each evaluation, the pack's sensible, latent and vapour exchange is as
  implicit in them as the canopy air's other sources: a deposition from warm saturated air cannot take
  more vapour than the canopy air holds. The **snow-base → soil-top conductance** is the series
  resistance of the half snow layer and the top soil node; $`k_{snow}\ll k_{soil}`$ throttles it as the
  pack deepens, the physical decoupling that caps the winter soil surface, and it becomes the soil's
  top BC against the march's own top-node temperature. The top node's $`k_{soil}`$ is the soil heat
  solver's own conductivity of layer 1 at the step's start (its water, ice fraction and texture), so a
  frozen or a dry top layer conducts as the soil column says it does.
- **Commit** (`snow_take_exchange`): the pack takes the march's time integral of the exchange, the
  energy and vapour the canopy air and the soil received from it, so the pack, the canopy air and the
  soil close together in the whole-column ledgers; its temperature reads off the new energy.

### The ground's radiation

The two-stream sees one ground, the bare soil and the pack blended by the cover fraction $`f`$: its
albedo and emissivity are the cover-weighted means, and it emits at the temperature at which that
emissivity gives the two tiles' emission,

```math
\bar\varepsilon\,\sigma T_{eff}^4 = f\,\varepsilon_{snow}\sigma T_{snow}^4 + (1-f)\,\varepsilon_{soil}\sigma T_{soil}^4 .
```

What reaches the ground in each band, $`S_b`$ and $`L^\downarrow`$, is then shared at each tile's own
optics (`ground_tile_radiation`). Per unit of its area the bare soil takes
$`\sum_b S_b(1-\alpha_{soil,b}) + \varepsilon_{soil}(L^\downarrow-\sigma T_{soil}^4)`$, at the top
soil node's temperature at the step's start, and the snow $`\sum_b S_b(1-\alpha_{snow,b})`$ of
shortwave, with $`\varepsilon_{snow}(L^\downarrow-\sigma T^4)`$ of longwave formed in its balance at its
own temperature. Weighted by the cover, the two give the blended ground's absorption, so the canopy and
the ground agree on the ground's radiation. Half-covered ground under 300 W m⁻² of shortwave splits it by
albedo, not by area: the snow's half absorbs about 38 W m⁻² of ground area and the soil's 117, where an
even split would give each 77.

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
is ice-aware — removing `enthalpy_vapor` from an ice-referenced layer debits sublimation
(vaporization + fusion) automatically — so the vapour side must be too: the liquid curve alone would
treat the surface as ice for *enthalpy* and as liquid for *vapour pressure*. Over ice $`e_{sat}`$ is
**10 % lower at −10 °C, 22 % at −20 °C and 34 % at −30 °C**, so the liquid curve would overstate the
driving gradient by those factors. The ice form is within 0.1 % of Murphy & Koop (2005) at −10 °C and
0.9 % at −30 °C.

Ground evaporation from a frozen soil does not use the blend: its liquid sits on the soil freezing
curve, whose potential puts the liquid's vapour pressure at the ice's, so Philip's humidity on that
potential over liquid saturation already carries the ice ([soil_biophysics](soil_biophysics.md),
"Frozen soil"). The blend is deliberately **not** used for dewpoint conversion or diagnostic VPD:
dewpoint is *defined* over liquid, so an ice branch there would mis-convert the forcing.

The **Niu-Yang (2007) snow-cover fraction** $`\mathrm{snowfac}=\tanh(\text{depth}/\text{scale})`$
(`snow_cover_fraction`), with scale = `ny07_a` × `[aerodynamics].soil_rough` ×
($`\rho_{snow}`$/100 kg m⁻³)^`ny07_m` (Niu & Yang's ground roughness, the micro-relief a thin pack has to
fill: 6.25 cm at the defaults, so 1 cm of snow covers 16 % of the ground and 10 cm 92 %), ramps the ground optics albedo and the soil-top-BC blend, and area-weights ("sub-column") all boundary exchange so a thin patchy pack barely exchanges
(stable and continuous through $`\mathrm{snowfac}\to0`$). Snow-off is bit-identical to no store.
Ground evaporation is such a tile flux: only the snow-free share evaporates, at the bare soil's own
resistances, $`E=(1-\mathrm{snowfac})\,\rho\,(q_g-q_{air})/(r_{aero}+r_{soil})`$. Inflating $`r_{aero}`$ by
$`1/(1-\mathrm{snowfac})`$ instead is not the same whenever the dry-surface-layer resistance $`r_{soil}`$
is positive: it divides only the aerodynamic leg, so it over-predicts bare-soil evaporation at
intermediate cover, most where a dry surface makes $`r_{soil}`$ dominant. The two agree only at no cover
and full cover. The canopy's roughness sees the pack through its own burial fraction, the pack's depth
over the canopy height, not this cover fraction ([canopy_aerodynamics](canopy_aerodynamics.md) §1).

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Snow / surface water | `snow_state_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |

## Code map

| Concept | Routine |
|---|---|
| snow energy / base conductance | `meds_surface`: `snow_surface_balance`, `snow_take_exchange`, `snow_base_conductance` (`snow_energy_step`, the two in one step, for the tests) |
| the stage before the march, the exchange in it, the commit | `meds_surface`: `prepare_snow_stage`, `ground_surface_exchange` (from `meds_fast_integ_derivs`' `surface_derivs`); `meds_fast_integ_update`: `update_dynamic_state` |
| snow mass / cover / melt | `meds_surface`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| prognostic snow type | `meds_patch_state_vars`: `snow_state_t` |

## References
- Murphy & Koop (2005), *Q. J. R. Meteorol. Soc.* 131:1539 — vapour pressure over ice.
- Niu & Yang (2007), *JGR* 112:D21101 — snow-cover fraction.
- Design record: `docs/dev_plans/archive/MEDS_SNOW_DESIGN.md`.
