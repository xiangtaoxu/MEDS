# Soil biophysics

The soil column is two coupled prognostic sub-columns on one ED2 negative-$z$ grid: a **water** column
(mixed-form Richards) and a **thermal** column (internal-energy heat diffusion). They are coupled by
construction — the thermal step reads the just-updated moisture, and the moving water carries liquid
enthalpy across the soil boundaries. The water column is solved together with the pond and the
plants' water as one network ([soil_plant_water](soil_plant_water.md)); this page is the soil's part.

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (**prognostic internal energy**, not temperature; **closed-budget
discipline**) and how the stores are woven each `dt_fast`. Both principles are load-bearing here: the
thermal column is prognostic in volumetric internal energy so freeze/thaw is a read-off of the
thermodynamic inverter, and both sub-columns close a machine-precision budget every step.

---

## The soil water column

One patch's prognostic soil moisture $`\theta_k`$ (+ ponded surface water) is advanced over `dt_fast`
by the soil-plant water network's step (`water_step`): implicit Euler on adaptive sub-steps, with
Newton on every store at once ([soil_plant_water](soil_plant_water.md)). For the soil layers this is the
mixed-form Richards equation on the ED2 negative-$z$ grid. Per layer,

```math
C(\psi_k)\,\frac{\partial\psi_k}{\partial t}
= \frac{\partial}{\partial z}\!\left[K(\theta)\Big(\frac{\partial\psi}{\partial z}+1\Big)\right] - S_k
\qquad(1)
```

with matric potential $`\psi(\theta)`$, hydraulic conductivity $`K(\theta)`$, and specific moisture
capacity $`C=d\theta/d\psi`$ from the **van Genuchten–Mualem** (default) or **Campbell / Clapp-Hornberger**
retention curves (`meds_soil_hydr_lib`). In effective saturation
$`S_e=(\theta-\theta_{res})/(\theta_{sat}-\theta_{res})`$ and matric head $h<0$ [m], van Genuchten's
curve is $`S_e=[1+(\alpha|h|)^n]^{-m}`$ with $`m=1-1/n`$, and Campbell's is $`h=\psi_{sat}\,S_e^{-b}`$,
saturated between its air-entry head $`\psi_{sat}<0`$ and zero, with $`K=K_{sat}\,S_e^{2b+3}`$. Both
keep their parameters in one pair, `[soil].curve_par_a` and `curve_par_n` (α [m⁻¹] and $n$, or
$`\psi_{sat}`$ [m] and $b$), so every reader of a column (the water network, the predawn seed, the
diagnostics) evaluates the curve the column was built with. The van Genuchten conductivity is Mualem's,

```math
K = K_{sat}\,S_e^{L}\,\big[1-(1-S_e^{1/m})^m\big]^2, \qquad L = -1 \qquad(2)
```

with pore connectivity $L=-1$ (Schaap & Leij 2000, `VG_PORE_CONNECTIVITY`), as ED2 uses with this
curve: with Mualem's 0.5 a soil of $n$ near 1.4 conducts less than ED2's floor of
1.16×10⁻¹³ m s⁻¹ from about −0.7 MPa, far wetter than measured (UNSODA), and −1 keeps it above that
to about −2 MPa. MEDS has no film-flow term
for drier soil (see `docs/ROADMAP.md`). Interface conductivity is ED2's rule (`rk4_derivs`): $\ln K$
interpolated linearly between the two nodes to the face, a thickness-weighted geometric mean,

```math
K_{k+1/2} = K_k^{\,1-w}\,K_{k+1}^{\,w}, \qquad w = \frac{\Delta z_k}{\Delta z_k+\Delta z_{k+1}} \qquad(3)
```

since the face lies $`\Delta z_k/2`$ below node $k$ and $`\Delta z_{k+1}/2`$ above node $k+1$; equal layers
give the plain geometric mean. Across a wet layer over a dry one the face conductivity is small, so a
wetting front, or capillary rise into a dry profile, advances more slowly than it would under an
upstream pick.

Each layer is solved in its effective saturation on the dry side and its potential on the wet side.
Above saturation it stores water through a specific storage of 10⁻⁴ m⁻¹, so a layer pushed past
saturation raises its potential and drives the excess up to the pond and down to the next layer: there
is no saturation clip and no residual floor. As a layer dries its potential falls without bound
towards $`\theta_{res}`$ and the conductivity of its faces with it, so no flow can empty it.

The **root sink** is the network's own: each cohort's uptake from each layer it roots in, solved with
the layer's water. There is no soil-side wilting point: the soil→root conductance falls with
$`K(\theta)`$ and a layer drier than the wood gives nothing, so a wilting threshold here would cut the
same water twice. **Ground evaporation** combines a Philip pore-space relative humidity with a
Swenson-Lawrence dry-surface-layer (DSL) resistance in series with $`r_{aero}=1/g_{g,net}`$ from the
aerodynamics kernel, at the end-of-sub-step state (`ground_evaporation_rate`); as the top layer dries,
Philip's humidity goes to zero and evaporation stops by itself. The DSL resistance grows as the top layer
dries below `[soil].dsl_theta_init` of saturation, its thickness and tortuosity referenced to the
air-dry water content at a matric head of −10⁴ m (CLM5's): referenced to the bulk $`\theta_1`$ instead,
the tortuosity would fall faster than the layer thickens and a wetter soil would evaporate less. The
tortuosity is Millington–Quirk's, $`\phi_{air}^{10/3}/\phi^2`$, rather than CLM5's form in the
Clapp–Hornberger $b$; with the air-filled porosity $`\phi_{air}`$ fixed the two differ by a constant
factor, so the moisture response has the same shape. It is the **single authority** for the
ground latent flux (it drives both the CAS vapour twin and the ground energy balance's LE, so no
double-count).

### Boundary conditions, and what the aquifer option actually means

The **top boundary** is the pond, a store of its own: rain enters it, and it exchanges water with the
top layer by the Darcy flux

```math
q_{pond} = \sqrt{K(\theta_1)\,K_{sat}}\;\Big(1 + \frac{h-\psi_1}{-z_1}\Big) \qquad(4)
```

from the pond, treated as a saturated layer as thick as the top one, to the top node at depth $`z_1`$:
eq. (3) with $w = 1/2$, $h$ the pond's head and $`\psi_1`$ here in metres. It is solved implicitly and is
two-way: a top layer pushed past saturation exfiltrates to the pond. The pond holds up to
`w_pond_max` (its first millimetre is a depression store that fills before any head builds); above
that, water runs off at the pond's temperature. The top layer's own $`K(\theta_1)`$ alone would seal a dried surface: near
residual water content it is about $10^{-15}$ of $`K_{sat}`$ for a van Genuchten loam, so after a dry
season the pond overflows and the soil never re-wets. With the geometric mean the capacity at that
state is about 0.1 mm/h, and it rises by an order of magnitude with the first millimetre taken up,
so the 5 mm pond buffers a storm and the top layer re-wets within hours.

The **bottom** either drains or sits on an aquifer. Drainage is the unit-gradient flux scaled by a
drainage factor $f \in [0, 1]$, $`q = f\,K(\theta_n)`$, always **downward**; free drainage, a leaky
base and a sealed one are its values of $f$ (`drain_factor`), not separate cases:

| option | flux at the base | use for |
|---|---|---|
| `free_drain` (default) | drainage, $f = 1$: $`q = K(\theta_n)`$ | a deep, unseen water table |
| leaky (soil maps only) | drainage, $`f = f_{leak}`$, the map's (0.05 in ED-SoilMaps v2) | soil on weathered rock that drains slowly |
| `bedrock` | drainage, $f = 0$: $`q = 0`$ | an impermeable base |
| `aquifer` | $`q = K_{bot}\,(\psi_n/\Delta + 1)`$, $`\Delta = \Delta z_n/2`$, **two-way** | a shallow water table at the column base |

`[soil].water_bottom_bc` names `free_drain`, `bedrock` or `aquifer`; a leaky base comes only from the soil
maps (§ The column from the soil maps), which give each cell its own $`f_{leak}`$.

**`aquifer` is a boundary condition, not a store.** The water table is *defined* to sit at the column
base, so the modelled soil column **is** the unsaturated zone above it. There is no aquifer bucket, no
baseflow, no prognostic or diagnosed water-table depth, and so no Dunne saturation-excess runoff and no
Zeng–Decker equilibrium correction to key off one. $`K_{bot}`$ is upstream-weighted: $`K(\theta_n)`$
when the flux is downward, $`K_{sat}`$ when it reverses, because the upstream cell is then the saturated zone (using
$`K(\theta_n)`$ there under-predicts capillary rise).

Two things follow that are easy to get wrong:

- The flux **reverses upward** once $`|\psi_n| > \Delta`$, i.e. once the deepest layer is dry enough to
  pull against the saturated zone. So `aquifer` is the **wet-site** boundary — right for riparian,
  floodplain and wetland columns, and wrong for an upland one, where it will supply water a real
  hillslope would not.
- The head term is what makes `aquifer` differ from `free_drain`: applied unconditionally, the
  deep-water-table limit $`q = K(\theta_n)`$ is free drainage, whatever store is put behind it. A
  water-table bucket that could rise into the column without saturating it would add bookkeeping, not
  a water table.

Every step closes a machine-precision water budget:

```math
\text{mass\_resid} = \Delta W_{stores}
 - \Delta t\,(P - E_{soil} - \text{drainage} - \text{uptake} - \text{runoff}) \approx 0 \qquad(5)
```

and exports the time-mean per-face Darcy flux and each layer's uptake for the soil-energy solve's
water enthalpy (§ Water–thermal coupling).

### The soil layers

Water and heat share one grid. Every column cuts one list of layer bottoms at its own depth
(`soil_layer_grid`): SoilGrids 2.0's depths to 1 m (5, 15, 30, 60 and 100 cm), 150 and 200 cm, then 3,
4.5, 6.5, 9, 12, 16, 21 and 27 m, each layer at most twice the one above. The last layer ends at the
column depth, and a remainder thinner than half of the list layer it falls in joins the layer above,
so no column ends in a sliver.

- **The top layer is 5 cm in every column.** That is thin enough for the diurnal heat wave, whose
  damping depth is about 0.1 m. The drying surface needs no thinner layer: the dry-surface-layer
  resistance stands in for it.
- **A shallow soil gets a few layers rather than thin ones:** 0.3 m is 5, 10 and 15 cm.
- **Per-layer soil maps line up with the layers,** and every column of a region shares one depth axis.

`[soil].layer_bottoms` replaces the list for a site with a measured profile, and the column
then ends at its last bottom unless `depth` cuts it.

### Hydraulic parameters from soil texture

`soil_texture_to_hydraulics` (`meds_soil_hydr_lib`) gives a layer's van Genuchten curve and saturated
conductivity from what soil maps report: sand $Sa$, silt $Si$ and clay $Cl$ (% of the fine earth,
< 2 mm), organic carbon $OC$ (%), cation exchange capacity CEC (cmol(c) kg⁻¹), pH in water, the fine
earth's bulk density $\rho_b$ (g cm⁻³) and the volume share of coarse fragments $f_v$ — the
quantities SoilGrids 2.0 maps per layer. The curve is Hodnett & Tomasella's (2002, Table 8), fitted to
tropical soils and kept by ED2 as `soil_hydro_scheme` 2; each polynomial is in percent:

```math
\begin{aligned}
\ln\alpha &= \big(-2.294 - 3.526\,Si + 2.440\,OC - 0.076\,\mathrm{CEC} - 11.331\,\mathrm{pH} + 0.019\,Si^2\big)/100\\
\ln n &= \big(62.986 - 0.833\,Cl - 0.529\,OC + 0.593\,\mathrm{pH} + 0.007\,Cl^2 - 0.014\,Sa\,Si\big)/100\\
\theta_{sat} &= \big(81.799 + 0.099\,Cl - 31.42\,\rho_b + 0.018\,\mathrm{CEC} + 0.451\,\mathrm{pH} - 0.0005\,Sa\,Cl\big)/100\\
\theta_{res} &= \max\Big[0,\ \big(22.733 - 0.164\,Sa + 0.235\,\mathrm{CEC} - 0.831\,\mathrm{pH} + 0.0018\,Cl^2 + 0.0026\,Sa\,Cl\big)/100\Big]
\end{aligned}
\qquad(6)
```

with $\alpha$ per kPa, turned into per metre of water head. The saturated conductivity of the fine
earth follows Ottoni et al. (2019) from the water that drains between saturation and −33 kPa:

```math
K_{sat} = 19.31\ \mathrm{m\,d^{-1}}\,\big(\theta_{sat} - \theta(-33\ \mathrm{kPa})\big)^{1.948} \qquad(7)
```

Coarse fragments hold no water, so the layer's curve is the fine earth's scaled by $1-f_v$; they slow
the flow by Saxton & Rawls (2006, eqs. 19 and 22), written for their share of the mass $f_m$, with
$\rho_p = 2.65$ g cm⁻³ the density of the grains and of the fragments:

```math
\theta_{sat,res} \leftarrow (1-f_v)\,\theta_{sat,res}, \qquad
K_{sat} \leftarrow K_{sat}\,\frac{1-f_m}{1-f_m\,(1-\tfrac{3}{2}\rho_b/\rho_p)}, \qquad
f_m = \frac{f_v}{(\rho_b/\rho_p)(1-f_v)+f_v} \qquad(8)
```

**Why this function.** Against 12,297 measured retention curves (GSHP; Gupta et al. 2022) it gives the
smallest error in plant-available water ($\theta$ from −33 to −1500 kPa) of the functions tried —
Cosby (1984), Saxton & Rawls (2006), Tomasella & Hodnett (1998) and Rosetta3: a median of 0.026
m³ m⁻³ from measured texture and 0.035 from SoilGrids at 0.1°, and 0.023 in tropical clays, where
Rosetta3 overstates it by 0.07–0.09 and Cosby by 0.04. Its $K_{sat}$ is among the least biased against
SoilKsatDB (Gupta et al. 2021) in the tropics, and at BCI its curve follows the flat measured dry
branch (Kupers et al. 2019). It is closed-form, so it costs nothing at run time.

**Where it is weak.** In boreal, organic-rich and acid soils CEC and low pH raise $`\theta_{res}`$, so
the available water comes out about 0.04 m³ m⁻³ low (Hyytiälä's topsoil: $\theta$ at −1500 kPa 0.39
against 0.07 measured). Like every texture-based function, $K_{sat}$ scatters by about a factor of ten
and runs low where soil structure adds large pores.

### The column from the soil maps

With the `[maps]` soil maps (ED-SoilMaps v2; `docs/configuration.md` "Input maps"),
`apply_map_inputs` (`meds_config`) builds each cell's column from what `meds_maps` read. A region's polygons and a site
on the ED_ERA5land grid build it the same way.

- **The layers.** The column ends at the map's `soil_depth`, on the layers above.
- **Each layer's soil.** A layer takes the texture and chemistry of the map layers it overlaps (SoilGrids'
  six, to 2 m), each weighted by the overlap. Below 2 m it takes the deepest map layer's.
- **Hydraulics and heat, per layer.** Each layer gets its own curve and $K_{sat}$ by eqs. (6)–(8). Its
  thermal properties come from the same texture by CLM5's mineral-soil formulas (Lawrence & Slater 2008), with
  $\theta_{sat}$ the layer's porosity and silt left out as in CLM5:

  ```math
  k_{solid} = \frac{8.80\,Sa + 2.92\,Cl}{Sa + Cl}, \qquad
  k_{dry} = \frac{0.135\,\rho_d + 64.7}{2700 - 0.947\,\rho_d}, \qquad
  c_{dry} = \frac{2.128\,Sa + 2.385\,Cl}{Sa + Cl}\,10^6\,(1-\theta_{sat}) \qquad(9)
  ```

  Here $`\rho_d = 2700\,(1-\theta_{sat})`$ kg m⁻³ is the dry bulk density; $k$ is in W m⁻¹ K⁻¹ and
  $c_{dry}$, the solids' heat capacity per m³ of layer, in J m⁻³ K⁻¹ (MEDS counts the water itself).
- **The start.** Each layer starts at its own field capacity, $\theta(\psi_{fc})$.
- **The bottom water boundary** is the map's `bc_class`:
  - free drain, drainage at $f = 1$;
  - leaky, drainage at the map's `f_leak` (a leaky cell whose `f_leak` is outside (0, 1] stops the
    run);
  - sealed, drainage at $f = 0$ (bedrock);
  - aquifer.
- **The thermal anchor.** The Dirichlet bottom conducts to the map's `deep_temp`, the cell's mean air
  temperature plus 1 K, at $d/\sqrt2$ below the bottom node, $d = 1.973$ m (§ The bottom thermal
  boundary).

A cell whose maps hold no value stops the run, naming the cell. The products are made complete in
their own pipelines; MEDS does not fill them.

A map-built column is layered by design rather than one curve fitted to the thickness-weighted
column: its top layers hold the surface's organic, high-CEC soil and its deep layers the subsoil's.

---

## The soil thermal column

`soil_energy_step_implicit` (`meds_soil_energy`) advances the prognostic volumetric internal energy
$`E_k`$ [J m⁻³] by an implicit BE heat-diffusion solve, a Thomas sweep (`meds_linear_algebra`: `thomas_solve`)
on **the same negative-$z$ geometry** as the water column. At each state it inverts
$`E_k\to(T_k,\text{fliq}_k)`$, forms the ice-aware thermal conductivity $`\kappa(\theta,\text{fliq})`$
and effective volumetric heat capacity $`C_{eff}(\theta,\text{fliq})`$, solves for $`T^{n+1}`$ with a
top Neumann flux $`G_{top}`$ and the bottom boundary of § The bottom thermal boundary, then commits a
**conservative** energy update from the $`T^{n+1}`$ conductive faces plus optional upwind
water-enthalpy advection:

```math
E_k^{n+1} = E_k^n + \frac{\Delta t}{\Delta z_k}\big[(hf_k-hf_{k-1}) + (qwf_k-qwf_{k-1})\big]
            + \Delta t\,q_{src,k} \qquad(10)
```

Temperature and liquid fraction are re-diagnosed from the committed energy — so **freeze/thaw is the
internal-energy plateau**, captured for free and unconditionally (ice-aware
$`\kappa`$/$`C_{eff}`$; the zero-curtain is tested, cooling a wet layer pins `soil_temp` at the triple
point while `soil_fliq` absorbs the fusion enthalpy). The closed residual is `energy_resid`
$`=\Delta E - \Delta t\,(G_{top}-\text{bottom}-\sum\text{root\_heat\_sink})\approx0`$, with the
infiltrating water's enthalpy, the top face's $`qwf_0`$, among the inputs. A sibling
`soil_energy_time_deriv` exposes the same flux divergence as an explicit RHS (faces at $`T^n`$) for the
RK45 integrator; ESDIRK2's implicit stages take the BE step.

The conductivity is Johansen's, made ice-aware (`soil_thermal_conductivity`, `meds_therm_lib`). With the
relative saturation $`S_r=\theta/\theta_{sat}`$, floored at 0.05 where the logarithm leaves its range of
validity, the Kersten number blends Johansen's unfrozen form with the linear frozen one by the liquid
fraction $`f_{liq}`$, and the saturated conductivity is the geometric mean of the solids, the liquid and
the ice:

```math
K_e = f_{liq}\,(\log_{10}S_r + 1) + (1-f_{liq})\,S_r,\quad 0\le K_e\le 1, \qquad
\kappa_{sat} = k_{solid}^{\,1-\theta_{sat}}\,k_w^{\,\theta_{sat}f_{liq}}\,k_i^{\,\theta_{sat}(1-f_{liq})}, \qquad
\kappa = K_e\,\kappa_{sat} + (1-K_e)\,k_{dry} \qquad(11)
```

$`k_{solid}`$ and $`k_{dry}`$ come from the texture (eq. (9)); $`k_w`$ and $`k_i`$ are the conductivities of
liquid water and ice.

### The bottom thermal boundary

Two boundary conditions are selectable, through `[soil].thermal_bottom_bc`.

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
hf_n = -\,g_{deep}\,(T_n - T_{deep}), \qquad g_{deep} = \frac{\kappa_n}{\ell} \qquad(12)
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
\quad\Longrightarrow\quad \ell = \frac{d}{\sqrt{2}} \qquad(13)
```

So the anchor **depth** is a physical choice, and the default is derived rather than fitted. For the
default column (2 m on the default layers) in a mid-latitude mineral soil at $`\theta = 0.3`$:
$`d = 1.973`$ m and the bottom node at $`|z_{node,7}| = 1.75`$ m, giving
`deep_depth` $`= 1.75 + 1.973/\sqrt{2} = 3.145`$ m.

Measured against the analytic profile (`test_soil_annual_damping`, a homogeneous 2 m column of ten
layers thickening with depth, bottom node at 1.727 m and so its anchor at 3.12 m, driven by a harmonic
$`G_{top}`$ to periodic steady state — the constant-$`\kappa`$, constant-$`C`$ case where
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
layers under the hydrologically active column (see `docs/ROADMAP.md`). The anchor buys
roughly a factor of ten in this metric, not exactness.

**Choosing `deep_temp` is the user's job, and the loader insists on it.** It is the mean annual soil
temperature below the damping depth — close to the mean annual air temperature of the driving forcing,
and a site property like latitude. An error in it is a steady flux $`g_{deep}\,\Delta T`$ into the
column base, so it biases deep-soil temperature in the annual mean; there is no defensible default and
a silent one would reintroduce, in a new place, exactly the bias this boundary condition exists to
remove.

## Water–thermal coupling

The water and thermal columns are coupled by construction: the thermal step reads the just-updated soil
moisture (for $`\kappa,C_{eff}`$), and the water the network moves carries liquid enthalpy: across each
face at the temperature of the node it leaves, out of each layer by root uptake at that layer's own
temperature (`root_heat`, per layer), and between the pond and the top layer (`water_heat_out`).
Infiltration crosses the top face at the pond's temperature. Drainage has no face term: its enthalpy
leaves with the bottom layer's sink, `root_heat_sink`, beside the root uptake. Soil evaporation is not a
face either: its water leaves as vapour, and its enthalpy, liquid plus latent, has already left in the
ground's latent flux inside $`G_{top}`$.

Every face, the top one included, is applied by one upwind rule at the time level of the conduction it
sits beside ($`T^{n+1}`$ in the step, $`T^n`$ in the tendency). In the supercooled-liquid datum each face
term is large, about 1 MJ per kg of water moved, so the heat the moving water delivers is a small
difference of large numbers, and a face valued at another temperature would swamp it.

## Canopy interception

Per-cohort **canopy interception** is a separate top→bottom cascade (`intercept_canopy_layer`, in
`meds_plant_biophysics`): a capacity-limited bucket with a Beer interception fraction, each
cohort's throughfall feeding the next, so the net throughfall reaching the soil top is what the
infiltration boundary sees. A film left above its capacity, when a cohort sheds leaf area under a full
film, drips the excess to the layer below, so the water stays in the column. A cohort holds its films
per plant, as it holds its tissue water, and the
fast loop works with them per m² of ground. So a split, a patch fusion or a disturbance changes a
cohort's density and leaves its films alone, and the plants that die take their films to the ground
with their tissue water.

In the fast loop the cascade runs once per `dt_fast`, at the step's start, as capture only; film
evaporation is the energy stages' own flux. Each cohort holds one combined leaf-and-wood bucket, filled
tallest first and shared between leaf and wood by capacity, and the update reaches both integrators as
a constant rate whose Euler step over `dt_fast` reproduces it. Under a snowpack nothing is intercepted:
the pack takes the rain and the snow, and snowfall below the new-snow threshold forms no pack and
passes the canopy as liquid throughfall. The films are clamped at their boundary after the march, not
inside it. A film evaporates at a conductance frozen at the step's start and only roughly rescaled by
what it holds, so it can overdraw below zero; left negative, it would create water in the next step's
interception, so it is floored at zero and the shortfall is booked as a negative outflow, while water
above capacity leaves as outflow.

## Soil heterotrophic respiration

The heterotrophic soil CO₂ flux $`R_h`$ comes from **one authority**: the CENTURY transfer matrix
(`meds_soil_biogeochem`: `heterotrophic_respiration_matrix`, reached through
`patch_heterotrophic_respiration`). It is a soil-carbon decomposition process (the slow
biogeochemistry lives in `src/slow_dynamics/soil/`; [soil_carbon](soil_carbon.md)), but its flux enters the fast loop as a source term of
the CAS CO₂ twin
(see [canopy_air_space_biophysics](canopy_air_space_biophysics.md)). The driver is its single authority
and passes the resulting flux into the CO₂ source.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Soil water | `soil_water_column_t` | `theta(:)`, `w_surface` | `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |

## Code map

| Concept | Routine |
|---|---|
| soil water (the network) | `meds_fast_integ_backward_euler`: `water_step` ([soil_plant_water](soil_plant_water.md)); a layer's water relations: `meds_soil_hydr_lib`: `soil_water_relations` |
| soil water flows | `meds_soil_water`: `soil_water_flux`, `ground_evaporation_rate`, `bottom_flux` |
| soil layers and profile | `meds_soil_derive`: `SOIL_LAYER_TEMPLATE`, `soil_layer_grid`, `soil_layers`, `build_soil_profile`; the types, `soil_params_t` ([soil]) and `soil_profile_t` (what the solvers read), in `meds_soil_params` |
| retention curves (vG / Campbell) | `meds_soil_hydr_lib`, with the matric potential `soil_matric_potential` [MPa] |
| hydraulic parameters from texture, eqs. (6)–(8) | `meds_soil_hydr_lib`: `soil_texture_to_hydraulics` |
| thermal properties from texture, eq. (9) | `meds_therm_lib`: `soil_texture_to_thermal` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| explicit thermal tendency (RK45) | `meds_soil_energy`: `soil_energy_time_deriv` |
| bottom thermal boundary, eqs. (12)–(13) | `meds_soil_energy`: `bottom_heat_face`, `deep_conductance` |
| the column from the soil maps | `meds_maps`: `read_map_inputs`; `meds_config`: `apply_map_inputs` |
| canopy interception | `meds_plant_biophysics`: `intercept_canopy_layer` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_matrix`; `meds_fast_integ_frozen`: `patch_heterotrophic_respiration` |
| prognostic soil types | `meds_patch_state_vars`: `soil_water_column_t`, `soil_energy_column_t` |

## References
- ED2 `../ED2/ED/src/dynamics/`: `lsm_hyd.f90`, `soil_respiration.f90` — the soil stores in the reference model.
- Celia, Bouloutas & Zarba (1990), *Water Resour. Res.* 26:1483 — mixed-form Richards.
- Brenner & Cancès (2017), *Adv. Water Resour.* — variable switching for Richards' equation.
- Mualem (1976), *Water Resour. Res.* 12:513; Schaap & Leij (2000), *Soil Sci. Soc. Am. J.* 64:843 — the conductivity and its pore connectivity.
- van Genuchten (1980); Clapp & Hornberger (1978) — soil retention curves.
- Hodnett & Tomasella (2002), *Geoderma* 108:155 — van Genuchten parameters from texture for tropical soils.
- Ottoni et al. (2019), *J. Hydrol.* 575:1345 — saturated conductivity from the drainable porosity.
- Saxton & Rawls (2006), *Soil Sci. Soc. Am. J.* 70:1569 — coarse fragments.
- Gupta et al. (2021), *Earth Syst. Sci. Data* 13:1593 (SoilKsatDB) and (2022), *Sci. Data* 9:444 (GSHP) — the measurements the function was chosen on.
- Swenson & Lawrence (2014), *JGR-Atmos.* 119:10299 — dry-surface-layer soil evaporation.
- Lawrence & Slater (2008), *Clim. Dyn.* 30:145 — the mineral and dry-soil thermal properties CLM uses, eq. (9).
- Design records: `docs/dev_plans/archive/MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `docs/dev_plans/archive/MEDS_ENERGY_BALANCE_DESIGN.md`.
