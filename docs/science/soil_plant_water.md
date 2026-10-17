# Soil and plant water: one network

Each `dt_fast`, one patch's water is solved as a single network: the ponded surface water, every soil
layer, and each cohort's wood and leaves are **stores**, joined by **paths** along which water flows
down its potential gradient. The stores' curves (how much water a store holds at a given potential)
and the paths' flow laws are the process library: [soil_biophysics](soil_biophysics.md) for the soil
retention curves and the soil's boundary processes, [plant_hydraulics](plant_hydraulics.md) for the
tissue pressure–volume curves, xylem vulnerability and root geometry. This page is how they are
joined and solved. The design and its measurements are in
`docs/dev_plans/MEDS_SPAC_WATER_NETWORK_PLAN.md`.

Potentials are in **MPa** and water in **kg m⁻² of ground**. A soil layer's potential is its matric
potential; its total head adds $`\rho g z`$ ($`\rho g`$ = `grav_head`, $z<0$ below the surface).

## 1. The network

| Store | Holds | Joined to |
|---|---|---|
| pond (core node 1) | ponded water, $`W = \rho\,S(h)`$ | rain; soil layer 1 |
| soil layer $k$ (core node $k+1$) | $`\rho\,\Delta z_k\,\theta_k`$ | layers $k\pm1$; every cohort rooting in it; evaporation (layer 1); drainage (bottom layer) |
| a cohort's wood | sapwood + fine-root water of all its plants | its rooted layers; its leaves |
| a cohort's leaves | leaf water of all its plants | its wood; transpiration |

The soil layers and the pond are the **core**. Each cohort is a **group** of two stores that touch
only each other and the core, and the core only through the wood. That structure is what the linear
solve exploits (§3).

The paths:

- **Soil faces** — Darcy flow between adjacent layers, $`q = \rho\,K_{face}\,[(\psi_k-\psi_{k+1})/(\rho g)
  + \Delta z_{node}]/\Delta z_{node}`$, with ED2's thickness-weighted geometric-mean face conductivity
  $`K_k^{1-w}K_{k+1}^{w}`$.
- **Pond to top layer** — the same law, the pond counting as a saturated layer as thick as the top
  one ($w=1/2$). It is two-way: a top layer pushed past saturation exfiltrates to the pond.
- **Root uptake** — from layer $k$ into a cohort's wood, $`u_k = g_k\,K_k\,r\!\left(\psi_k + \rho g z_k -
  \psi_W\right)`$, with $`g_k`$ the cohort's Gardner single-root geometry in the layer
  ([plant_hydraulics](plant_hydraulics.md) §4) and $r$ the ramp below. **No hydraulic
  redistribution:** a layer drier than the wood gives nothing and takes nothing.
- **Xylem** — sapflow from wood to leaves, $`s = K_{eff}\,(\psi_W - \psi_L - \rho g H)`$, with the
  Kirchhoff mean conductance $`K_{eff} = k_{cond}\,[\Phi(\psi_W)-\Phi(\psi_L)]/(\psi_W-\psi_L)`$ (its
  midpoint expansion as the two potentials meet), plus a $10^{-12}$ share of $`k_{cond}`$ so a fully
  embolised cohort's leaves stay joined to its wood.
- **Boundaries** — rain into the pond; ground evaporation from layer 1 (Philip's pore humidity with
  the Swenson–Lawrence dry-surface-layer resistance); drainage out of the bottom layer
  (`[soil].bottom_bc`: drainage at a factor of free drainage, 1 to 0 for bedrock, or a head-driven
  two-way aquifer); and each cohort's transpiration out of its leaves.

Every flow is written with its derivatives in each store's unknown, so Newton's Jacobian is exact.
Every switch in a flow is smoothed (C¹), so Newton never meets a corner:

```math
r(x) = \begin{cases} 0 & x \le 0\\ x\,t\,(2-t),\ t = x/w & 0 < x < w\\ x & x \ge w \end{cases}
\qquad w = 10^{-3}\ \text{MPa (uptake)} \qquad(1)
```

## 2. The stores and their unknowns

Each store is solved in the unknown in which its curve is best behaved:

- A **soil layer** uses its effective saturation $S_e$ on the dry side and its potential on the wet
  side, joined (C¹) at the curve's point of largest capacity: $`S_e = (1+m)^{-m}`$ for van Genuchten, the
  air-entry point for Campbell. Above saturation the layer stores water through a specific storage
  $`S_s = 10^{-4}`$ m⁻¹, so a layer pushed past saturation raises its potential and sends the excess up to
  the pond and down to the next layer: there is no saturation clip. Dry, $S_e > 0$ always, and the
  potential falls without bound towards $`\theta_{res}`$, so no flow can empty a layer and ground
  evaporation stops by itself.
- A **tissue** uses its symplastic relative water content $R$, on which the pressure–volume curve is
  linear in water and has a finite slope $d\psi/dR$ everywhere. A leafless cohort's leaves hold no
  water; their $R$ stays defined.
- **Between slow steps** a tissue keeps its water **mass**, not its potential: a change in biomass leaves
  the water as it is, and the potential is read again on the new capacity (grown biomass reads as a
  slightly lower $R$). A store above what today's biomass can hold is clamped to its saturated content,
  and an empty one is seeded at its starting potential (`reconcile_tissue_water_capacity`; neither
  amount is booked yet, see `docs/ROADMAP.md`). Tissue lost to turnover takes its water with it, in
  proportion to the carbon lost (`shed_turnover_water`), and a plant that dies takes all of its own
  (`shed_mortality_water`); that water reaches the ground as `shed_water_rate`, valued like the
  precipitation it joins.
- The **pond** uses its head $h$, with a 1 mm depression store: $S(h) = 0$ for $h\le0$, $h^2/2\delta$
  for $0<h<\delta$, $h-\delta/2$ above. An empty pond holds nothing and its head follows the soil
  below.

## 3. The solve

Over one `dt_fast`, **implicit Euler on sub-steps**. A sub-step of length $h$ solves the water
balance of every store at its end,

```math
\frac{W_i(\tau) - W_i^{start}}{h} - F_i(\tau) = 0 \qquad(2)
```

with $F_i$ the net inflow, by **Newton** in the unknowns $\tau$:

- The Newton matrix $`J = \partial W/\partial\tau\,/h - \partial F/\partial\tau`$ is solved by eliminating
  each cohort onto the core (`meds_arrowhead_solve`, a block arrowhead system): each cohort's 2×2 block is inverted in closed
  form, its rank-one Schur update folded into the core, and the core (pond + layers) LU-factored. The
  cost grows as cohorts × layers, not as (cohorts + layers)³.
- A Newton step takes no soil layer or tissue more than 90 % of the way to empty, and moves no dry
  layer's $S_e$ or tissue's $R$ by more than 0.2; an Armijo line search on the sum of squared, scaled
  imbalances then backtracks it if needed.
- Newton has converged when every store's imbalance over the sub-step is below 1 % of the step's
  absolute water tolerance.

**Error and sub-steps.** A sub-step's error is implicit Euler's local error, half the change in the
inflows across it, filtered twice through $J^{-1}$ (two back-substitutions with the last
factorisation). The filter leaves slow modes alone and divides a stiff, already-settled exchange's
error by its stiffness squared — what step doubling would measure, at a third of the cost. Each
store's error is weighed against $`\text{atol}\cdot\text{scale} + \text{rtol}\cdot|W|`$: `[soil]`
`atol` and `rtol`, read as water content for soil layers and the pond and as relative water content
for tissues. A sub-step that fails Newton is cut to a quarter; one over tolerance is cut by the usual
order-one controller, which also grows the next one. Every step starts with one sub-step of
`dt_fast`; most need no more (BCI: 1.07 sub-steps and 3.5 Newton iterations per step on average).

**Conservation.** Each store's water is advanced by the flows themselves, $`W \mathrel{+}= h\,F`$, so
every kilogram that leaves one store arrives in another or crosses a boundary, whatever Newton's
residual. Two stores hang off the network by one path and can hold far less than Newton's slack —
the pond, and a near-leafless cohort's leaves. They take their water from their own curve instead,
and their single path carries exactly the flow that closes their balance; the store across that path
(the top layer, the wood) absorbs the slack, which is small against what it holds. After each sub-step
water above `[soil].w_pond_max` spills as runoff at the pond's temperature.

## 4. Transpiration and the energy step

Each cohort's transpiration is its stomata's demand, limited by what its leaves hold:

```math
E = E_{demand}\,\frac{\min(R_L,\ R_{close})}{R_{close}} \qquad(3)
```

with the corner rounded over 0.02 in $R$. $`R_{close}`$ is the leaf's symplastic water content at the
PFT's fast-limb closing potential (`pft_wstress_psi_close`, twice the turgor-loss point by default).
Above it the stomatal controls act unchanged; below it a leaf with its stomata shut loses its last
water in proportion to what it has left, so every step has a solution.

$`E_{demand}`$ comes from the start-of-step surface solve, the one that also sets the leaf energy
balance's coefficients. The leaf there pays the full vapour enthalpy of what it transpires and is
credited the liquid enthalpy of the sapflow that replaces it; the credit is estimated from a first
solve, as if the sapflow matched the transpiration, and the solve is repeated. (Without the credit,
the leaf in that solve is charged ~4.4×10⁵ J kg⁻¹ too much, runs cooler and transpires about 4 % less.)

The energy stages then **hold each cohort's transpiration at the network's** for the whole step, as
they hold ground evaporation: the leaf energy balance takes `le_slope = 0` and
$`\text{le\_ref} = h_{evap}E/(1-f_{wet})`$. The canopy air's vapour credit is then the leaf's water
debit by construction. Within a step transpiration does not follow the canopy air's changing
humidity; the next step's surface solve restores the response. Measured on the BCI dry season
(2013-01-15 to 04-15): the hourly latent heat at `dt_fast` = 900 s is within 1.9 W m⁻² of 150 s,
against 2.9 W m⁻² when transpiration is refreshed at every stage instead, and at 150 s the two
agree to 0.6 W m⁻² in any hour.

## 5. Enthalpy

The network moves water; the soil-energy solve and the tissue stores carry its enthalpy:

- Water crossing a soil face, the pond face or the base carries the liquid enthalpy of the node it
  leaves.
- Root uptake leaves each layer at that layer's own temperature (`root_heat`, per layer).
- A tissue's heat is its frozen capacity times its temperature; the change in its water is booked at
  its start-of-step temperature (`water_store_enth`), and arriving water heats it only by its
  temperature difference.
- The pond carries its own enthalpy through each sub-step: rain joins it first, what infiltrates
  leaves at the mixture's temperature, water the top layer pushes up arrives at that layer's
  temperature, and the spill leaves at the pond's.

**The plant water's heat.** With sapflow $s$, total root uptake $U=\sum_k U_k$ and transpiration $E$
(all per ground area per second), $u(T)$ the liquid internal energy (`internal_energy_liquid`),
$`T_{up}=T_W`$ when $s\ge 0$ and $T_L$ otherwise, and $`Q_{loss}=\sum_k U_k\,u(T_k)`$ the enthalpy the
roots draw from the layers, the leaf and the wood receive

```math
q_{wl} = s\,u(T_{up}) - (s-E)\,u(T_L), \qquad
q_{w} = Q_{loss} - U\,u(T_W) - s\,\big[u(T_{up}) - u(T_W)\big] \qquad(4)
```

and the change in their water is booked, in the ledgers, at the start-of-step temperatures:

```math
H_{store} = (U-s)\,u(T_W) + (s-E)\,u(T_L) \qquad(5)
```

($`q_{wl}`$, $`q_{w}`$ and $`H_{store}`$ are `qwflux_wl`, `q_wood_net` and `water_store_enth` in
`build_column_frozen`; the upwind temperatures are frozen at the start of the step, which keeps the
energy ledger exact.) The heat store, capacity × temperature, has no term for water mass coming or
going, so counting the arriving water's whole enthalpy (about 4.4×10⁵ J kg⁻¹ from the liquid datum)
as heat would warm the wood whenever it refills: about 10 W m⁻² every night at BCI, and at a storm
front, where wet soil refills dry-season wood with about 7 kg m⁻² in a quarter hour, 13 K above the
canopy air and +600 W m⁻² of sensible heat at zero net radiation.

**The ground inflow.** Under a pack the water reaching the pond is meltwater, at the melt temperature.
On bare ground it is the rain plus the snowfall too light to start a pack, both valued at the air
temperature above the canopy (the forcing moved to the canopy-air top, [forcing](forcing.md) §8), as
`snow_accumulate` values a pack's snowfall (`precip_enthalpy`). The mixture's enthalpy per kilogram
becomes a temperature through `temp_of_liquid_enthalpy`, the exact inverse of
`internal_energy_liquid`; below $`T_{3ple}`$ that temperature stands for water that must still melt,
which the pond and the soil melt with soil heat. Valuing the snow as liquid would create $L_f$ per
kilogram at the boundary, and valuing the rain at the canopy air's temperature would hand the canopy
and the soil heat the rain did not bring: at a cold front the canopy air is still warm (297.2 K under
294.9 K air at BCI on 2017-04-17).

## 6. Work counters

| Output | Counts |
|---|---|
| `work_water_nsub_site` | accepted sub-steps |
| `work_water_newton_site` | Newton iterations |
| `work_water_reject_site` | sub-steps rejected and retried shorter |
| `work_water_nfail_site` | sub-steps taken at the shortest length although unconverged or over tolerance |

## Code map

| Concept | Routine |
|---|---|
| the solver: assembly, Newton, sub-steps, commit | `meds_spac_solver` (`numerics/`): `spac_allocate`, `spac_set_state`, `spac_advance` |
| the linear solve (block arrowhead, Schur complement) | `meds_arrowhead_solve`: `arrowhead_allocate`, `arrowhead_factor`, `arrowhead_solve` |
| the stores' water relations in their unknowns | `meds_water_retention`: `soil_water_relations`, `tissue_water_relations`, `pond_water_relations` and their inverses; `VG_PORE_CONNECTIVITY` |
| plant flows | `meds_hydr_lib`: `rhizosphere_conductance`, `root_water_uptake`, `root_water_uptake_profile`, `xylem_max_conductance`, `xylem_sapflow`, `supply_limited_transpiration`, `smooth_max0` |
| soil flows | `meds_soil_water`: `soil_water_flux`, `ground_evaporation_rate`, `bottom_flux` |
| one patch's step | `meds_fast_frozen`: `run_spac_solver` (called from `build_column_frozen`) |
| plant water over the step | `meds_fast_be_stage`: `advance_plant_water` (ARK; constant rates, RK45 integrates the same rates in its stages) |

## References

- Brenner & Cancès (2017), *Adv. Water Resour.* — variable switching for Richards' equation.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — the mixed form.
- Gardner (1960), *Soil Sci.* 89:63 — single-root rhizosphere conductance.
- Hairer & Wanner (1996), *Solving ODEs II* — implicit Euler, Newton and error estimates for stiff systems.
- Xu, Medvigy, Powers, Becknell & Guan (2016), *New Phytol.* 212:80 — plant hydraulics in ED2.
