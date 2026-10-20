# Vegetation energy dynamics

The per-cohort **leaf and wood** tissue energy balance: how each cohort's canopy surfaces exchange
radiation, sensible heat, and water vapour with the canopy air space, and how their temperatures are
relaxed exactly over the step toward their quasi-steady balance, with their own heat capacity.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven each `dt_fast`. Absorbed shortwave and net longwave enter from the radiation kernel
([canopy_radiation_transfer](canopy_radiation_transfer.md)); the boundary-layer conductances come from
the aerodynamics kernel ([canopy_aerodynamics](canopy_aerodynamics.md)); the transpiration demand comes
from leaf gas exchange ([leaf_gas_exchange](leaf_gas_exchange.md)) and is supply-limited by plant
hydraulics ([plant_hydraulics](plant_hydraulics.md)).

---

## Leaf and wood energy balance

`veg_energy_balance` (`meds_plant_biophysics`) is the one kernel for a cohort's leaf **or** wood; wood is
the case without transpiration. Its conductances and absorbed radiation are frozen for the step
(Category 0, [numerical_scheme](numerical_scheme.md) §2), the saturation terms follow the canopy-air
state, and every flux is linearized about the canopy-air temperature $`T_{cas}`$. The balance is then
linear in the tissue's offset $`\Delta T = T - T_{cas}`$: a forcing $`\mathrm{numer}`$ that does not
depend on $`\Delta T`$, and a coupling conductance $`\mathrm{denom}`$ that multiplies it,

```math
\mathrm{numer} = A_{sw} + A_{lw} + Q_{adv} - \mathrm{LE}_{ref}, \qquad
\mathrm{denom} = h + \mathrm{LE}_{slope} + 4\,\varepsilon\,\sigma\,T_{cas}^3\,\mathrm{AI} \qquad(1)
```

$`A_{sw}`$ and $`A_{lw}`$ are the absorbed shortwave and net longwave from the two-stream; $`Q_{adv}`$ is
the enthalpy the sapflow carries in (`q_extra`); $\mathrm{AI}$ is LAI for leaves and WAI for wood, and
$\varepsilon$ the PFT's leaf or wood emissivity. The sensible coefficient is
$`h=\rho\,c_p\,g_{bh}\cdot2\,\mathrm{LAI}`$ for leaves (flat plates, both sides, `effarea_heat`) and
$`h=\rho\,c_p\,g_{bh}\cdot\pi\,\mathrm{WAI}`$ for wood (cylinders). Each latent pathway with conductance
$g$ linearizes the saturation deficit at $`T_{cas}`$:

```math
\mathrm{LE}_{ref} = h_v\,\rho\,g\,\big(q_{sat}(T_{cas}) - q_{cas}\big), \qquad
\mathrm{LE}_{slope} = h_v\,\rho\,g\,\frac{dq_{sat}}{dT}\bigg|_{T_{cas}} \qquad(2)
```

with $`h_v`$ the specific enthalpy of the vapour at $`T_{cas}`$, so the tissue pays what the canopy air
receives. The terms in eq. (1) sum the pathways, each weighted by its share of the surface: a leaf's
dry share $`1-f_{wet}`$ transpires through boundary layer and stomata in series; its wetted share
$`f_{wet}`$ evaporates the intercepted film through the boundary layer alone, paying $`h_v`$ less the
enthalpy the film's liquid arrived with. Wood has only the film pathway. Once the water network has
fixed the step's transpiration, the dry share carries it as a fixed latent term
($`\mathrm{LE}_{slope}=0`$), so the transpiration the leaf reports is exactly the water it lost.

**A thin tissue is a small tissue.** Every term of eq. (1), and the heat capacity of eq. (8), is a
per-plant quantity times the plant density: the absorbed radiation, the sapflow's enthalpy, the
conductances and the emission slope. A cohort split into two at half the density therefore reaches the
same temperatures and each half carries half of every flux, however thin the cohort is. For a leaf
every term also scales with its leaf area, so a small leaf sits at the balance of a large one. Nothing
floors the coupling: the heat store bounds the step (eq. (5) below). ED2 instead slaves a tissue below
a resolvable size to the canopy air (`leaf_resolvable`/`wood_resolvable`); a floor of that kind, in
ground-area units, pulls a thin cohort toward the canopy-air temperature, and two half-density cohorts
then differ from one.

**A tissue with no area** — a leafless deciduous cohort's leaves, or wood with no wood area — has no
coupling and no heat capacity, so eq. (1) gives it no temperature. It takes the temperature of what
it grows from: a leafless leaf takes its wood's (a bud sits on its twig), and wood with no area takes
the canopy air's. Nothing weights that value, since every flux the tissue carries and its heat
capacity are zero, and a new leaf reaches its own balance within a step. In the kernel, "no area" is a
coupling and a heat capacity both below its divide guard (`tiny_num`), which is where a leafless
cohort's leaf carbon, decaying toward 10⁻¹⁹⁰, puts its leaves.

The tissue relaxes toward the zero-inertia offset $`\Delta T_{diag}=\mathrm{numer}/\mathrm{denom}`$ with
its own heat capacity (next section), and reports its fluxes over the step average, eqs. (5)–(6).

**The emission base.** The two-stream emits each cohort at its effective radiative temperature
$`T_{rad}=\big(f_l T_{leaf}^4+(1-f_l)\,T_{wood}^4\big)^{1/4}`$, $`f_l=\mathrm{LAI}/(\mathrm{LAI}+\mathrm{WAI})`$,
from the start-of-step tissue temperatures ([canopy_radiation_transfer](canopy_radiation_transfer.md)
§3), so the absorbed longwave $`A_{lw}`$ already contains the tissue's own emission at $`T_{rad}`$.
`veg_energy_balance` linearizes the change in emission about $`T_{cas}`$, so the loss it books is
$`A_{lw}-4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T-T_{cas})`$
where the consistent one is $`A_{lw}-4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T-T_{rad})`$: a tissue warmer
than the canopy air loses an extra $`4\varepsilon\sigma T_{cas}^3\,\mathrm{AI}\,(T_{rad}-T_{cas})`$ and a
cooler one gains it. Linearizing about the lagged tissue temperature $`T_{rad}`$ instead
removes the residual but adds a lagged positive feedback (a warmer leaf raises the next step's $`T_{rad}`$), which drove a single-pass
integrator to NaN when it was tried (`docs/dev_plans/archive/MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §5) and has
not been tried under `ark`.

## One formulation: the tissue relaxes exactly

There is no diagnostic/prognostic *mode*. The tissue always carries a heat capacity and
`veg_energy_balance` always relaxes it — **exactly**, because under the Category-0 coefficient
freeze the equation is linear with a known timescale:

```math
\mathrm{cap}\,\frac{dT}{dt}=\mathrm{numer}-\mathrm{denom}\,(T-T_{cas}),\qquad \tau=\frac{\mathrm{cap}}{\mathrm{denom}}
\qquad(3)
```

The step therefore has a closed form, and the kernel uses it rather than discretising it. Two
**different** weights fall out, and both are needed:

```math
x=\frac{\Delta t}{\tau}=\frac{\mathrm{denom}}{a_{\rm store}},\qquad
w_{\rm end}=e^{-x},\qquad
w_{\rm avg}=\frac{1-e^{-x}}{x} \qquad(4)
```

with $`a_{\rm store}=\mathrm{cap}/\Delta t`$. `w_end` weights the old state at the **endpoint** — the
committed temperature. `w_avg` weights it in the **step average** — every flux the canopy air actually
receives:

```math
\Delta T_{\rm end}=(1-w_{\rm end})\,\Delta T_{\rm diag}+w_{\rm end}\,\Delta T_{\rm prev},\qquad
\Delta T_{\rm avg}=(1-w_{\rm avg})\,\Delta T_{\rm diag}+w_{\rm avg}\,\Delta T_{\rm prev} \qquad(5)
```

where $`\Delta T_{\rm prev}`$ is the start-of-step offset. The kernel evaluates eq. (5) on the start-of-step
imbalance $`r=\mathrm{numer}-\mathrm{denom}\,\Delta T_{\rm prev}`$, as
$`\Delta T_{\rm end}=\Delta T_{\rm prev}+w_{\rm avg}\,r/a_{\rm store}`$ and
$`\Delta T_{\rm avg}=\Delta T_{\rm prev}+\frac{1-w_{\rm avg}}{x}\,r/a_{\rm store}`$, which is the same
relaxation written with no division by denom: it holds as denom → 0, where a store with no coupling
warms by $`r\,\Delta t/\mathrm{cap}`$, and it shows that the heat store bounds the step whatever the
coupling. The reported sensible heat, evaporation and net radiation all use $`\Delta T_{\rm avg}`$:

```math
H = h\,\Delta T_{\rm avg}, \qquad
E = \frac{\mathrm{LE}_{ref}+\mathrm{LE}_{slope}\,\Delta T_{\rm avg}}{h_v}, \qquad
R_{net} = A_{sw} + A_{lw} - 4\,\varepsilon\,\sigma\,T_{cas}^3\,\mathrm{AI}\,\Delta T_{\rm avg} \qquad(6)
```

($E$ per latent pathway; $`Q_{adv}`$ stays out of $`R_{net}`$ because it is an internal soil-to-plant
transfer, not a radiative input at the column boundary). Using one weight for both is what breaks
conservation; pairing them closes the balance identically:

```math
a_{\rm store}\,(\Delta T_{\rm end}-\Delta T_{\rm prev})+\mathrm{denom}\cdot\Delta T_{\rm avg}=\mathrm{numer}
\qquad(7)
```

**"Diagnostic" is this kernel's zero-inertia limit, not a second code path.** As
$`a_{\rm store}\to 0`$ (no heat capacity) $`x\to\infty`$, both weights vanish, and the result is the
pure steady-state balance. As $`a_{\rm store}\to\infty`$ both weights → 1 and the tissue holds its
temperature. Nothing switches; the physics chooses, per cohort, every step, continuously.

### The store is ON, and how its capacity is sized

The store is always on (`build_column_frozen` in `meds_fast_frozen`). The capacity splits its two
halves across **different** masses, which is a physical distinction and not a detail:

- **dry tissue → ALL the wood**, at the dry specific heat. Heartwood is dead structure but still
  stores sensible heat, and a sapwood fraction defined on the bole under-counts branch wood, which is
  thin enough to be thermally active throughout.
- **water → the water network's whole wood store.** The network's wood store is the sapwood **plus
  the fine roots**, with one prognostic water mass $`W_{wood}`$ ([plant_hydraulics](plant_hydraulics.md)
  §1). The network moves all of that water at the wood's temperature
  ([soil_plant_water](soil_plant_water.md) eqs. (4)–(5)), so the heat capacity counts all of it, the
  fine roots' share included; heartwood holds none, which is what makes it heartwood. A store that
  drains through a dry season holds less heat.

Per unit ground area, with $B$ the dry biomass (carbon × 2 × $`n_{plant}`$):

```math
\mathrm{cap}_{wood} = c_{wood}\,B_{wood} + c_{liq}\,n_{plant}\,W_{wood}, \qquad
\mathrm{cap}_{leaf} = c_{leaf}\,B_{leaf} + c_{liq}\,W_{leaf} \qquad(8)
```

with $`c_{wood}`$ = 1217 J kg⁻¹ K⁻¹, ED2's specific heat of **dry** wood at 15 °C
(`tissue_c_wood_dry`; 103.1 + 3.867 T), $`c_{leaf}`$ = 3200 J kg⁻¹ K⁻¹ for dry leaf (`tissue_c_leaf`;
both in `meds_constants`), $`W_{wood}`$ the wood store's water per plant and $`W_{leaf}`$ the leaf's
internal water, both from the water network. Fresh wood's specific heat (about 2700 J kg⁻¹ K⁻¹ at
100 % moisture) is not used: the water is counted on its own, so a fresh-wood value would count it
twice. The fast loop and the state layer's `cohort_tissue_heat_capacity` build the capacity the same
way, so the fusion weights and the slow ledger use the same capacity. The intercepted film is a
separate store and is not in $`\mathrm{cap}_{leaf}`$.

The wood's area and its sapwood come from stem-size allometries (`meds_allometry_lib`; each PFT's
`[allometry]` `wai_b1/b2`, `sapwood_area_b1/b2`), with $D$ the diameter [cm]. The wood area index is
$`\mathrm{WAI} = n_{plant}\,b_1\,D^{b_2}`$ (`dbh_to_wai`), and the sapwood is the share
$`f_{sap} = \min\big(1,\ b_{1SA}\,D^{b_{2SA}}/(\tfrac{\pi}{4} D^2)\big)`$ of the basal area
(`sapwood_fraction`), capped at one because a small stem is sapwood throughout. WAI follows the stem
rather than the leaf area on purpose: a wood area proportional to LAI would make
$`\tau_{wood}`$ (wood mass over wood area) nearly independent of size and 6–20 times too short, and
would mis-scale the wood's boundary layer, longwave emission and sensible-heat coefficient, which all
take WAI (eq. (1)).

Sized this way, with a saturated store (`wood_water_sat` = 1 kg per kg of carbon),
`τ_wood ≈ 600–760 s` and `w_end ≈ 0.22–0.31` at `dt_fast = 900 s` — wood carries **genuine memory**.

**Most of the wood's capacity is the sapwood and its water**: dry wood holds 1217 J kg⁻¹ K⁻¹ against
water's 4186, so the heartwood adds little. And because `w_end` is small for wood while `w_avg` is
not, **the store's effect lives in the FLUXES, not in the temperature**: "the wood barely lags" does
not imply "the wood barely matters", and comparing tissue temperatures alone makes the store look
negligible.

Measured effect, from a 2×2 over (store, longwave) on a 2 h night window (`test_fast_loop`): the store
leaves the night canopy air **0.85 K warmer**, since leaf and wood carry a heat capacity comparable to
the canopy air's and pass their stored heat to it as they cool. It also absorbs a third of the
net-longwave signal (2.44 K without the store, 1.60 K with it), which is why a test that thresholds
nighttime cooling is pinned against the store-on control.

**Why not backward Euler.** BE is the $`w_{\rm end}=w_{\rm avg}=1/(1+x)`$ approximation. A leaf at
production `dt_fast` has $`x\approx144`$: the exact endpoint weight is $`\sim 5\times10^{-63}`$ against
BE's $`0.0069`$, and an SDIRK2 tableau gives $`-0.031`$ — a *sign-alternating* artificial memory for a
mode whose true memory is nil. Under BE the prognostic tissue's only measurable effect would be the
discretisation's own artefact.

Note `denom` **excludes** $`a_{\rm store}`$: it is the coupling conductance — a property of the
tissue's surroundings, not of its own inertia. Since $`a_{\rm store}=\mathrm{cap}/\Delta t`$, the
ratio $`\mathrm{denom}/a_{\rm store}`$ is $`\Delta t/\tau`$ and the step size cancels, so the kernel
needs no `dt` argument and stays a pure tendency-safe routine.

### What the timescales actually are

$`\tau_{\rm leaf}\approx 12`$–20 s, and — a useful result — it is **LAI-independent**: `cap` and
`denom` both scale with LAI, so every leaf cohort shares one τ and there is no per-cohort structure to
resolve. A leaf is genuinely at equilibrium at any `dt_fast` the model runs at.

Wood is different. The allometric sizing above gives $`\tau_{\rm wood}\approx 600`$–760 s with a
saturated store, less as it drains; a physical estimate (branch wood plus the bole's diurnal damping
depth) is ≈1300 s for branches and ≈3690 s for boles — so wood memory at `dt_fast = 900 s` is real,
not negligible. Wood is also treated as
**internally isothermal**, a deliberate simplification: Bi ≈ 7.45 for a bole skin against the Bi ≪ 0.1
lumped capacitance needs, so real wood is a diffusion problem with no single τ. That is a limitation of
the one-node wood store, not of the integrator.

## The tissue's energy is its state

Each cohort carries the internal energy of its leaves and of its wood, `leaf_energy` and `wood_energy`
[J per plant], and their temperatures are derived from them. The energy is counted in the soil's
datum, `temp_to_internal_energy` with the water liquid, which is the datum of every other store and of
every water flow:

```math
E = D\,T + W\,c_{liq}\,(T - T_{sc}) \qquad(9)
```

```math
T = \frac{E + W\,c_{liq}\,T_{sc}}{D + W\,c_{liq}} \qquad(10)
```

with $D$ the dry heat capacity and $W$ the water of eq. (8), per plant, and $`T_{sc}`$
`tsupercool_liq`; $`dE/dT`$ at fixed water is eq. (8)'s capacity (`tissue_energy`, `tissue_temp`). The
water stays liquid at any temperature: a live tissue is taken not to freeze, so a tissue below 0 °C is
supercooled liquid. Water that comes or goes at the tissue's temperature carries
$`u_{liq}(T)=c_{liq}(T-T_{sc})`$ with it, which leaves the temperature where it was.

**In the fast step.** The kernel relaxes each tissue in temperature, at the heat capacity held over the
step (eqs. (3)–(7)). The commit (`commit_tissue_energy`) then adds two things to the energy: the heat
the tissue relaxed through, $`\mathrm{cap}_0(\bar T-T_0)`$, with $`\bar T`$ the march's time-average; and
the enthalpy its own water gained ([soil_plant_water](soil_plant_water.md) eq. (5)). It takes the
temperature from that energy, eq. (10), with the water the tissue holds at the end of the step. The
water gained is thereby mixed in at the start-of-step temperature it was valued at,
$`T_{end}=(\mathrm{cap}_0\,\bar T+\Delta W\,c_{liq}\,T_0)/\mathrm{cap}_1`$. The next step starts from
that energy, and the whole-column ledger counts the tissue store as the sum of the energies.

**In the slow step.** Matter that a slow operator adds or removes does so at its tissue's temperature.
That covers tissue grown or shed, and water taken into new leaves, shed to the ground, or held above a
store's capacity. The operator changes the mass and holds the temperature (`hold_cohort_tissue_temps`),
so the energy changes by the matter's heat, which the slow ledger declares.

- **Fusion** sums the energies, nplant-weighted per plant like the tissue water. The merged tissue's
  temperature follows from the merged masses, so a fusion conserves the tissue energy exactly.
- **A split** keeps the parent's temperatures in both daughters.
- **A leafless crown** holds no energy and takes its wood's temperature.

As the fast loop reads each cohort, it checks that the stored temperatures are the ones the energies
give.

**Why energy, and in this datum.** Valuing the tissue in the same datum as the water that crosses it
makes the tissue store a quantity the model can add up across steps, fusions and splits. Valued
instead as capacity × temperature, it sits $`W c_{liq} T_{sc}`$ (about 2.4×10⁵ J per kg of water) away
from the water's datum. The store a step ends with is then not the one the next step begins with: it
is re-valued by $`\Delta W c_{liq}(T_1-T_0+T_{sc})`$ each step, about 1.7×10⁶ J m⁻² in the quarter
hour in which wet soil refills dry-season wood at a BCI storm front.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Leaf / wood tissue | `cohort_block` (per plant); `patch_biophys_t` (per m² of ground) | `leaf_energy`, `wood_energy` [J plant⁻¹], eq. (9) | `leaf_temp`, `wood_temp` [K], eq. (10) |

## Code map

| Concept | Routine |
|---|---|
| leaf/wood balance + exact exponential relaxation, eqs. (1), (3)–(7) | `meds_plant_biophysics`: `veg_energy_balance` |
| coefficients $h$, $g$, emission slope | `sensible_heat_coeff`, `leaf_transp_coeff`, `leaf_film_coeff`, `lw_emission_slope` |
| latent terms, eq. (2); the per-cohort call | `meds_fast_time_derivs`: `surface_derivs` |
| tissue heat capacity, eq. (8), and $`a_{\rm store}`$ | `meds_fast_frozen`: `build_column_frozen` |
| tissue energy and temperature, eqs. (9)–(10) | `meds_therm_lib`: `tissue_energy`, `tissue_temp` |
| the fast step's energy commit | `meds_column_state_ops`: `commit_tissue_energy` |
| a cohort's temperatures from its energies; holding them through a change of mass | `meds_cohort_state_vars`: `cohort_tissue_temps`, `set_cohort_tissue_temps`, `hold_cohort_tissue_temps` |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |

## References
- ED2 `../ED2/ED/src/dynamics/rk4_derivs.f90` — the leaf/wood energy budget in the reference model.
- Design records (archived, `docs/dev_plans/archive/`): `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` (the
  exact-exponential tissue store), `MEDS_ENERGY_BALANCE_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.
