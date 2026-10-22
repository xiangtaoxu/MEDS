# Leaf gas exchange

How MEDS computes leaf-level photosynthesis, stomatal conductance, and their coupling. The code
lives in [`src/fast_dynamics/plant/meds_leaf_gas_exchange.f90`](../../src/fast_dynamics/plant/meds_leaf_gas_exchange.f90); the
public entry points are `solve_leaf_gas_exchange` and `leaf_gas_exchange_batch` in
`meds_leaf_gas_exchange`; the config-flattening convenience wrapper `leaf_gas_exchange` lives in
`meds_fast_config`.

Symbols: $A$ = net CO₂ assimilation, $A_g$ = gross assimilation, $R_d$ = leaf (dark) respiration
($A = A_g - R_d$); $C_a, C_s, C_i$ = CO₂ mole fraction (ppm) in ambient air / at the leaf surface /
in the intercellular space; $g_s$ = stomatal conductance **to water vapour**; $g_b$ = leaf
boundary-layer conductance to water vapour; $D$ = leaf-to-air vapour deficit; $\psi$ = water
potential. Two fixed diffusion constants appear throughout — the H₂O:CO₂ diffusivity ratios
$1.6$ through the stomata (`gsw_2_gsc`) and $1.4$ through the boundary layer (`gbw_2_gbc`).

---

## 1. The coupled A–Cᵢ–gₛ system

Assimilation is fixed by the intersection of two curves in $C_i$:

- **Biochemical demand** — the Farquhar–von-Caemmerer–Berry (C3) or Collatz (C4) rate $A(C_i)$
  (§2), which *rises* with $C_i$.
- **Diffusive supply** — CO₂ drawn down through the stomata:

```math
A = \frac{g_s}{1.6}\,(C_s - C_i), \qquad C_s = C_a - 1.4\,\frac{A}{g_b} \qquad(1)
```

  (the $C_s$ term applies only when the boundary layer is resolved, `use_boundary_layer`;
  otherwise $C_s = C_a$). This *falls* with $C_i$.

Because the stomatal models set $g_s$ as a function of $A$, $C_s$, and $D$, the whole system
reduces to **one nonlinear equation in $C_i$**, solved by bracketing and bisection (§5).

---

## 2. Photosynthetic demand $A(C_i)$

**C3 (Farquhar et al. 1980).** Three potential rates, co-limited:

```math
A_c = V_{cmax}\frac{C_i-\Gamma^*}{C_i + K_c\left(1 + O/K_o\right)},\quad
A_j = J\frac{C_i-\Gamma^*}{4C_i + 8\Gamma^*},\quad
A_p = 3\,\mathrm{TPU} \qquad(2)
```

where $A_c$ is Rubisco-limited, $A_j$ RuBP/light-limited, $A_p$ triose-phosphate-limited;
$`\Gamma^*`$ is the CO₂ compensation point without respiration and $K_c, K_o$ the Rubisco
Michaelis constants. The electron-transport rate $J$ follows the non-rectangular hyperbola
$`\theta J^2 - (I_2 + J_{max})J + I_2 J_{max}=0`$ (smaller root), with
$`I_2 = \tfrac12\,\phi_{\text{PSII}}\,\alpha_{\text{leaf}}\,\mathrm{PAR}`$. The shipped
$`\phi_{\text{PSII}}`$ = 0.74 makes $J$'s initial slope 0.37 electrons per absorbed photon, which at four
electrons per CO₂ ($A_j$'s $4c_i + 8\Gamma^*$) is the measured maximum CO₂ yield of C3 leaves, 0.093 per
absorbed photon (Long et al. 1993). Fv/Fm (0.85) in its place would fix CO₂ at the O₂-evolution yield,
0.106 (Björkman & Demmig 1987).

**$`\Gamma^*`$ is proportional to the oxygen partial pressure**, because it is set by Rubisco's
specificity for CO₂ over O₂: $`\Gamma^* = \tfrac{1}{2} O / S_{c/o}`$. The O₂ partial pressure is
the O₂ mole fraction times the air's pressure, as the CO₂ one is, so $`\Gamma^*`$ as a mole fraction
depends on the O₂ mole fraction and **not on pressure**. MEDS computes it as a mole fraction from the
O₂ and the pressure its reference value is stated at,

```math
\Gamma^*(T, O) = \frac{\Gamma^*_{25}}{P_\mathrm{ref}}\,f_{\mathrm{Arr}}(T)\;\frac{O}{O_\mathrm{ref}},
\qquad O_\mathrm{ref} = 0.209,\quad P_\mathrm{ref} = 100\ \mathrm{kPa} \qquad(3)
```

The shipped `gstar25` = 4.275 Pa is Bernacchi et al. (2001)'s 42.75 µmol mol⁻¹ at 21 % O₂, written
as a partial pressure at 100 kPa, so $`\Gamma^*`$ at 25 °C is 42.75 µmol mol⁻¹ at any elevation.
$`K_c`$ and $`K_o`$, the Michaelis constants of dissolved CO₂ and O₂, are partial pressures and are
divided by the local pressure; $`O/K_o`$ is then the ratio of the two partial pressures
(`rubisco_kinetics` in `meds_leaf_gas_exchange`, `meds.plant.leaf.rubisco_kinetics` in Python).
Dividing `gstar25` by the local pressure instead would raise $`\Gamma^*`$ by 11 % at 1000 m and by
43 % at 3000 m, cutting gross photosynthesis in dim light, where photorespiration sets the quantum
yield, by about 4 % and 14 % (25 °C, 400 µmol mol⁻¹, PAR 100).

The O₂ factor makes `o2_mol_frac` reach **both** places oxygen enters the C3 demand — the Michaelis
term $`K_c(1 + O/K_o)`$ and the compensation point — and is exactly 1 in any default run. It matters
for palaeo-atmosphere work (Carboniferous/Permian O₂ is reconstructed near 30–35 %) and for the low-O₂
sensitivity experiment that is the standard way to isolate photorespiration: measured on a
midday tropical leaf, halving O₂ raises $`A_\mathrm{net}`$ by 22 % and 35 % O₂ cuts it by 24 %.

**C4 (Collatz et al. 1992):** $A_c = V_{cmax}$, $A_j$ = a light-limited slope, $A_p = k_p C_i$
(PEP-case CO₂ limitation), with $`\Gamma^* \approx 0`$ (the CO₂-concentrating mechanism suppresses
photorespiration).

The three rates are combined into $A_g$ either as a sharp $\min(A_c,A_j,A_p)$ (`COLIM_MIN`) or as
two nested smoothing quadratics (`COLIM_QUADRATIC`, the same `quadratic_smaller_root` used for $J$).
Then $A = A_g - R_d$.

**Temperature.** $`K_c`$, $`K_o`$ and $`\Gamma^*`$ scale from 25 °C by the plain Arrhenius form, with
Bernacchi et al. (2001)'s activation energies ($`f_{\mathrm{Arr}}`$ in eq. (3)). $`V_{cmax}`$,
$`J_{max}`$, TPU, $`R_d`$ and the C4 $`k_p`$ take the form `[leaf_physiology].temp_response_form`
selects, plain or peaked (`meds_therm_lib`; `plant_respiration.md` eq. (5)). TPU and $`k_p`$ have no
temperature parameters of their own: both take $`V_{cmax}`$'s, and with them its high-temperature
deactivation under the peaked form.

### The co-limitation curvatures are not $`\theta_J`$

Each pathway has **its own pair** of co-limitation curvatures, and neither is `theta_j`:

| | first transition | second transition |
|---|---|---|
| C3 | `theta_cj_c3` — $A_c$ vs $A_j$ | `theta_ip_c3` — $\min(A_c,A_j)$ vs $A_p$ |
| C4 | `theta_cj_c4` — $A_c$ vs $A_j$ | `theta_ic_c4` — $\min(A_c,A_j)$ vs $A_p$ |

A co-limitation curvature says **how sharply the leaf switches between two limiting processes**;
`theta_j` is the curvature of the electron-transport hyperbola, describing light saturation of $J$.
They share units and a functional form, and nothing else. Implementations generally put the
co-limitation curvature much closer to 1 (CLM and FATES use ~0.98 for the C3 $A_c/A_j$ transition)
than the 0.7–0.9 that fits the $J$ hyperbola.

Using a $J$ curvature for the co-limitation would cost assimilation exactly where the model spends
most of its time, because $A_c$ and $A_j$ sit close together at ambient CO₂. Measured with
`meds.plant.leaf` at 25 °C and 101.325 kPa, $C_i = 280$ µmol mol⁻¹, with PFT 1's ratios
($`J_{max} = 1.8\,V_{cmax}`$, $`\mathrm{TPU} = 0.167\,V_{cmax}`$, $`R_d = 0.015\,V_{cmax}`$),
against the sharp $\min(A_c,A_j,A_p)$:

| `theta_cj_c3` / `theta_ip_c3` | gross $A$, saturating light ($J = J_{max}$) | net $A$, PAR 1500–2000 µmol m⁻² s⁻¹ |
|---|---|---|
| 0.85 / 0.85 (a $J$ curvature on both transitions) | **−27.4 %** | **−31 to −34 %** |
| 0.85 / 0.95 (on $A_c$/$A_j$ only) | −23.7 % | −27 to −31 % |
| 0.98 / 0.95 (shipped) | −9.8 % | −12 to −15 % |

At saturating light the shortfall is the same at every $V_{cmax,25}$ from 60 to 150, because all
three rates scale with $V_{cmax}$ (at $V_{cmax,25}$ = 60 the min is 14.51 µmol m⁻² s⁻¹ and the three
pairs give 10.53, 11.06 and 13.09). The finite-light column spans PAR 1500 and 2000 and
$V_{cmax,25}$ = 60, 90, 120, 150, with $J$ from the hyperbola (`theta_j` 0.70, absorptance 0.85,
$`\phi_{\text{PSII}}`$ 0.74); the net rate loses more than the gross because $R_d$ does not shrink
with it. Either way the loss is a systematic offset rather than a regime effect, so **no measured
$V_{cmax}$ would reproduce a measured assimilation rate**. The shipped $V_{cmax,25}$ presets sit
mid-range for their PFT descriptions and are not tuned to offset a smoothing loss.

Kernels: `rubisco_kinetics`, `assimilation_demand_c3`, `assimilation_demand_c4`, `electron_transport_j`,
`combine_limits`.

---

### Thermal acclimation (opt-in)

With `[leaf_physiology].thermal_acclimation = true`, the peaked form's entropy term and the capacity
ratio follow the growth temperature (Kattge & Knorr 2007):

```math
\Delta S_{V} = 668.39 - 1.07\,T_g, \qquad
\Delta S_{J} = 659.70 - 0.75\,T_g, \qquad
\frac{J_{max,25}}{V_{cmax,25}} = 2.59 - 0.035\,T_g \qquad(4)
```

with $T_g$ in °C. $`\Delta S`$ sets where the response **peaks**, so a warm-grown stand runs a higher
optimum: measured across a 20 K range of growth temperature, the $V_{cmax}$ optimum moves from
**28.2 °C** (grown at 10 °C) to **35.0 °C** (grown at 30 °C).

Three things about this are deliberate:

- **The ratio acclimates too.** Shifting $`\Delta S`$ alone would acclimate the *shape* of each
  response while pinning the two branches in a fixed proportion, which is not what the study
  measured — a warm-grown plant invests relatively less in electron transport.
- **$T_g$ is the growth temperature, not the leaf temperature.** The fit is calibrated against the
  mean **air** temperature of the preceding weeks, so MEDS tracks an exponential running mean of the
  daily-mean air temperature (window `acclim_window_days`, default 30 d) and feeds that. An
  instantaneous leaf temperature would use the relation far outside what it was fitted to.
- **It is one per site, not per cohort**, for the same reason: each site's running mean of its own air
  temperature, after its terrain ([forcing](forcing.md) §8, Each site of a polygon), and each site's
  leaf table. Per-cohort acclimation is a different question — it would need a relation calibrated on
  leaf rather than air temperature (see `docs/ROADMAP.md`).

The running mean is prognostic and slow, so it is written to the state file: rebuilding a month of
memory after every restart would make the acclimated optimum jump. It is applied by refreshing the
per-PFT leaf table once per slow step from the **config** reference coefficients, which are never
overwritten — so the refresh is idempotent. `thermal_acclimation` requires the peaked form; with the
Arrhenius form there is no $`\Delta S`$ to shift and the loader refuses rather than letting the flag
be silently inert.

## 3. Stomatal-conductance models

Selected by `stomatal_model` (`SM_LEUNING` | `SM_MEDLYN` | `SM_KATUL`). The first two give an
**explicit** $g_s$ law; Katul defines $g_s$ **implicitly** through an optimization.

### 3.1 Leuning (1995)

```math
g_s = g_0 + \frac{g_1\,A}{(C_s - \Gamma^*)\,(1 + D/D_0)} \qquad(5)
```

A semi-empirical law: conductance tracks assimilation, discounted by CO₂ drawdown
$`(C_s-\Gamma^*)`$ and by a hyperbolic humidity response $(1 + D/D_0)$. Parameters: `stomatal_g0`
($g_0$), `stomatal_g1` ($g_1$), `stomatal_d0` ($D_0$). Code: `stomata_gs_leuning`.

### 3.2 Medlyn et al. (2011) — unified stomatal optimization (USO)

```math
g_s = g_0 + 1.6\left(1 + \frac{g_1}{\sqrt{D}}\right)\frac{A}{C_s}, \qquad D \text{ in kPa} \qquad(6)
```

with $D$ floored at `medlyn_vpd_min` (`[leaf_physiology]`), since $g_1/\sqrt{D}$ is undefined at 0.
Derived from the same $\max(A-\lambda E)$ optimization as Katul but solved analytically under a
linearized $A(C_i)$, giving a closed form. Here $g_1 \propto \lambda^{-1/2}$ — a key relation used
by the water-stress limbs (§4). Parameters: `stomatal_g0`, `stomatal_g1`. Code: `stomata_gs_medlyn`.

### 3.3 Katul et al. (2010) — marginal-WUE optimality

Stomata are assumed to **maximize carbon gain net of a water cost**:

```math
\max_{g_s}\bigl(A - \lambda E\bigr) \qquad(7)
```

where $E = g_s D$ is transpiration and $\lambda$ is the **marginal water-use efficiency** — the
"price" of water (mol C per mol H₂O). The optimum satisfies $\partial A/\partial E = \lambda$.
Writing both fluxes in $C_i$-space (with $A = \tfrac{g_s}{1.6}(C_s - C_i)$, so
$`g_s = 1.6\,A/(C_s-C_i)`$ and $`E = 1.6\,D\,A/(C_s-C_i)`$) and setting
$`dA/dC_i = \lambda\,dE/dC_i`$ rearranges to the residual root-found on $C_i$:

```math
A'\,(C_s - C_i)^2 \;=\; 1.6\,D\,\lambda\,\bigl(A'(C_s - C_i) + A\bigr), \qquad A' \equiv \frac{dA}{dC_i} \qquad(8)
```

$A'$ is taken by central difference (the co-limited $A(C_i)$ has no clean analytic slope). Unlike
the explicit models, $g_s$ is recovered *after* solving, as $`g_s = 1.6\,A/(C_s-C_i)`$; if that
would fall below the cuticular floor $g_0$, the solver re-solves g0-pinned (§5). Parameters:
`katul_lambda25` ($\lambda_{25}$), and the water-stress terms below. Code: `residual_optimality`,
`katul_lambda`.

---

## 4. Water stress — a fast and a slow control

MEDS follows the two-factor decomposition of Sabot et al. (2022) / Zhou et al. (2013): drought
acts on **capacity** (biochemistry) and on the **stomatal** aperture independently. The two act on
different time scales, as in Xu et al. (2016)'s ED2 scheme: capacity on the leaf's own current water
potential, the aperture on the previous day's predawn potential.

### 4.1 The fast control: capacity, on the current leaf potential — $\beta_{ns}$

A linear ramp in the **current leaf** water potential (the step-start $`\psi_{leaf}`$ the plant-water
solver holds), applied to $V_{cmax}$, $J_{max}$ and TPU for **all** stomatal models (capacity
downregulation is a biochemistry effect, scheme-independent). $R_d$ is not scaled.

```math
\beta_{ns} = \mathrm{clamp}\!\left(\frac{\psi_{leaf} - \psi_{close}}{\psi_{open} - \psi_{close}},\,0,\,1\right),
\qquad \{V_{cmax}, J_{max}, \mathrm{TPU}\} \mathrel{*}= \beta_{ns} \qquad(9)
```

$`\psi_{open}`$ is the PFT's leaf turgor-loss point $`\psi_{tlp}=\pi_0\varepsilon/(\pi_0+\varepsilon)`$
and $`\psi_{close}`$ twice it, from the same pressure-volume curve the hydraulics solver stores leaf
water on (each PFT's own `[hydraulics]` `leaf_pi0` and `leaf_elastic_mod`). A PFT
file may set `wstress_psi_open` and `wstress_psi_close` instead. A leaf keeps its capacity while it
holds turgor and loses it as it wilts.

Every stomatal model ties $`g_s`$ to $`A`$, so this is what closes the stomata *within the day*: as the
leaf wilts, $`A`$ and with it $`g_s`$ fall, continuously, to the cuticular $`g_0`$. At
$`\psi_{leaf}\le\psi_{close}`$ the leaf by day has $`A_g=0`$, $`A_n=-R_d`$ and $`g_s=g_0`$; it then loses
only the cuticular water $`g_0 D`$. That makes $`g_0`$ the leaf's minimum conductance, a physical trait,
rather than a fitted intercept: a large $`g_0`$ keeps transpiration running from a wilted leaf.

### 4.2 The slow control: stomatal aperture, on the predawn potential — $\beta_s$ (Sabot 2022)

An exponential in the **predawn** water potential, downregulating the *aperture*:

```math
\beta_s = \min\!\bigl(1,\ \exp\bigl(s_{ref}\,(\psi - \psi_{onset})\bigr)\bigr) \qquad(10)
```

- **Leuning / Medlyn:** scale the slope, $`g_{1,\text{eff}} = g_1\,\beta_s`$.
- **Katul:** raise the marginal price, $`\lambda = \lambda_{25}\,\beta_s^{-e}`$ (`katul_lambda`).
  Because Medlyn gives $g_1 \propto \lambda^{-1/2}$, the exponent $e = 2$ makes the Katul and
  Leuning/Medlyn responses *identical* for the same $\beta_s$.

Parameters: `wstress_sref_stomata` ($s_{ref}$, ~2 MPa⁻¹), `wstress_lambda_exp` ($e$), and
`stomata_psi_onset` ($`\psi_{onset}`$, default half the PFT's turgor-loss point, so the gravity head a
tall tree has in wet soil does not read as drought). Above the onset $`\beta_s = 1`$: Sabot et al.
(2022, eq. 5) apply no stress while the soil is at field capacity, and $`\psi_{onset}`$ stands in for
their field-capacity threshold because a tall tree's gravity head is there even in wet soil.

> **Note.** The stomatal driver is the *supply-side* potential $\psi$ (via `leaf_env_t%psi`,
> default 0 = well-watered); the capacity driver is the current $\psi_{leaf}$, the *demand-side*
> tension. Sabot writes $\psi$ as soil/predawn potential; MEDS supplies the cohort's own **predawn
> (previous-day daily-max) leaf** potential — see below. This mirrors ED2's Katul
> `stoma_beta` (with $`s_{ref}\,e \equiv`$ ED2's `stoma_beta`, and $\psi \equiv$ ED2's
> `dmax_leaf_psi`); the divergences are in §4.3.

> **What drives $`\beta_s`$.** The **previous day's daily-maximum leaf water potential** — the
> model's predawn potential, which is what the field measurement behind $`s_{ref}`$ actually is, and
> the plant's own integration over its rooted profile (so it already carries rooting depth,
> per-PFT vulnerability, and the gravity head that no single soil layer's $`\psi`$ does). Unset
> cohorts are seeded from the surface-layer soil potential — the one place a genuine soil $`\psi`$
> enters this path.
>
> Carried on the cohort block as **`dmax_psi_leaf`** (the published value the kernel reads) plus
> **`dmax_psi_leaf_accum`** (the running accumulator). They are a *double buffer*, not a redundant
> pair — within a day one is read-only and the other write-only, and a single field cannot be both:
> reset it at day start and there is nothing left to read; never reset it and `max()` ratchets
> monotonically to the least-negative $`\psi`$ the run has ever seen, silently disabling the closure
> forever. The supply potential and soil $`\psi`$ coincide only in wet soil: the wood↔soil relaxation
> time is ~9 s at $`\theta`$ 0.25 but ~4.8 **days** at $`\theta`$ 0.10, so a droughted cohort never
> equilibrates overnight — which is why the kernel field is named plainly `psi`.
>
> Merging the two would also change the physics: the kernel would read *today's* partial maximum,
> which early in the morning is barely different from the instantaneous $`\psi_{leaf}`$ — the tight
> sub-step feedback the slow control exists to avoid. A genuine one-field design is a different
> algorithm: sample $`\psi_{leaf}`$ once a day at true predawn. It saves one value per cohort but
> adds a clock dependency to the fast loop and a sensitivity to which sub-step is caught; the daily
> maximum needs neither.

### 4.3 Where this diverges from ED2

ED2's `farq_katul.f90` (Xu et al. 2016's scheme, `H2O_PLANT_LIM = 4`) has the same two time scales:
a turgor-loss capacity term on the leaf potential at the start of each step, and a Manzoni-style
`stoma_beta` on the day's running maximum potential. The forms differ:

| | MEDS | ED2 (`farq_katul.f90`) |
|---|---|---|
| **capacity shape** | linear ramp in $\psi_{leaf}$ from $\psi_{tlp}$ to $2\psi_{tlp}$ | $`1/\bigl(1+0.1(\psi_{leaf}/\psi_{tlp})^{6}\bigr)`$, a 6th-power turgor-loss curve |
| **capacity targets** | $V_{cmax}$, $J_{max}$ **and** TPU | $J_{max}$ and quantum yield $\alpha$ only ($V_{cmax}$ and TPU are commented out) |
| **stomatal scope** | Leuning/Medlyn $g_1$ **and** Katul $\lambda$ | Katul $\lambda$ only |
| **stomatal driver** | predawn (previous-day daily-max) leaf $\psi$ | `dmax_leaf_psi`, the current day's running maximum |
| **parameter split** | $s_{ref}$ (shared with the $g_1$ limb) $\times$ $e$ (Katul-only) | one `stoma_beta` |

**Capacity shape and targets.** Both anchor the capacity term on $\psi_{tlp}$. ED2 applies it to
$J_{max}$ and $\alpha$ while leaving $V_{cmax}$ untouched, which makes it act only on the
light-limited branch. MEDS downregulates all three potential rates together, which is what the
underlying observations describe, and keeps the effect independent of which limitation a leaf happens
to be under. MEDS also applies it for **every** stomatal model, not just Katul, because capacity
downregulation is biochemistry and should not depend on the aperture scheme selected. The ramp is
linear, so its response to an error in $\psi_{leaf}$ is linear too; ED2's 6th power would amplify
the same error harder.

**Stomatal scope.** Restricting $\beta_s$ to Katul, as ED2 does, would leave the two explicit
schemes with no stomatal drought response at all. MEDS applies it to both, and the exponent
$e = 2$ makes them exactly equivalent because Medlyn gives $`g_1 \propto \lambda^{-1/2}`$ (§4.2).

**Parameter split.** MEDS's $`s_{ref}\,e`$ is precisely ED2's `stoma_beta` — the same exponential in
the same potential — so a value calibrated for one transfers to the other as
$`\texttt{stoma\_beta} = -s_{ref}\,e`$. For Katul alone only the product matters, so the split is
degenerate there; it earns its keep because $s_{ref}$ *alone* drives the Leuning/Medlyn limb, which
ED2 has no equivalent of.

### 4.4 No separate shut-off

$`\beta_s`$ scales $`g_1`$ only, so on its own it leaves a drought-stressed leaf the conductance
$`g_0 + 1.6\,A/C_s`$ (Medlyn's "1 +" term) for as long as it still photosynthesises. The fast control
removes the $`A`$-driven part as the leaf wilts, so no further factor on $`g_s`$ is needed: what is
left is $`g_0`$, the cuticular loss, which a real leaf has too. At a VPD of 3.5 kPa a $`g_0`$ of
0.28 mmol m⁻² s⁻¹ (example 04's calibration) loses 0.010 mmol m⁻² s⁻¹; the 0.01 mol m⁻² s⁻¹ of the
shipped PFT file loses 0.35 mmol m⁻² s⁻¹, which is why $`g_0`$ should be set to a measured minimum
conductance.

**The substomatal air is taken as saturated.** It is in equilibrium with leaf water at
$`\psi_{leaf}`$, so strictly its humidity is the Kelvin value
$`\mathrm{RH}=\exp\!\big(\psi/(\rho_w R_v T)\big)`$, not 1. Transpiration on
$`e_i=\mathrm{RH}\cdot e_{sat}`$ would shrink the gradient continuously as the leaf dries and
**reverse** it once $`e_i<e_a`$ — i.e. foliar water uptake — arresting transpiration with no
threshold parameter at all, at $`\psi=\rho_w R_v T\ln(\mathrm{RH}_{air})`$. MEDS applies this exact
relation to **soil** evaporation (`ground_evaporation_rate`'s `alpha_soil`) but not to leaves, so the
model assumes sub-saturated vapour over soil at −10 MPa but saturated vapour over leaf water at
−100 MPa.

The leaf RH cannot simply scale the latent flux. `le_ref` in `surface_derivs` is the latent flux
*linearised about* $`T_{cas}`$, with the leaf–air temperature difference carried by
`le_slope`·$`\Delta T_l`$, so scaling both terms by RH moves the base point and the slope
inconsistently: on a nearly turgid leaf it gives *condensation* at midday, and flooring RH does not
help, because the fault is the linearisation rather than the magnitude. Doing it properly means
rebuilding the flux *and* its temperature derivative from $`\mathrm{RH}\cdot q_{sat}(T_{leaf})`$
inside the leaf energy balance. Its likely value is as a route to **foliar water uptake**, not as a
drought response (see `docs/ROADMAP.md`).

> Note the Kelvin correction is deliberately kept out of the **stomatal** model: Medlyn/Leuning
> $`g_1`$ were calibrated against $`e_{sat}`$-based VPD, so feeding them a corrected VPD would
> silently re-tune $`g_1`$.

---

## 5. The coupled solver

`solve_leaf_gas_exchange` brackets $`C_i \in (\Gamma^*, C_a]`$ and bisects the residual to a
tolerance `CI_TOL_PPM` (the module's private `bisect_ci`):

- **Night / closed branch** — if the best-case net rate $A(C_a) \le 0$, stomata sit at $g_0$ and
  the solve returns immediately.
- **Two residual providers**, chosen by `ci_residual` — `residual_explicit_gs` (Leuning/Medlyn
  diffusion identity, eq. (1), and the g0-pinned fallback with $g_s := g_0$) and
  `residual_optimality` (Katul, eq. (8)). Both are `pure`.
- **Cuticular fallback** (`GS_CUTICULAR`) — if the model yields no sign change (e.g. Katul under
  strong stress wanting $g_s < g_0$), the solver retries g0-pinned, which always brackets when
  $A(C_a) > 0$.

Then $g_s$, $C_s$, $C_i$, and $E$ are assembled consistently.

---

## Parameters (config names)

| Symbol | Config key (`[pft]`) | Meaning |
|---|---|---|
| $g_0$ | `stomatal_g0` | cuticular / residual conductance |
| $g_1$ | `stomatal_g1` | stomatal slope (Leuning `--` / Medlyn kPa$^{0.5}$) |
| $D_0$ | `stomatal_d0` | Leuning humidity sensitivity |
| $\lambda_{25}$ | `katul_lambda25` | Katul marginal WUE (well-watered) |
| $`\psi_{open},\psi_{close}`$ | `wstress_psi_open`, `wstress_psi_close` | fast-limb ramp bounds (optional; default $`\psi_{tlp}`$, $`2\psi_{tlp}`$) |
| $`\psi_{onset}`$ | `stomata_psi_onset` | where the stomatal limb starts (optional; default $`\psi_{tlp}/2`$) |
| $s_{ref}$ | `wstress_sref_stomata` | stomatal-limb sensitivity (Sabot) |
| $e$ | `wstress_lambda_exp` | Katul $\lambda$ exponent on $\beta_s$ |

---

## References

- Farquhar, von Caemmerer & Berry (1980). *Planta* 149:78–90 — FvCB C3 photosynthesis.
- Collatz, Ribas-Carbo & Berry (1992). *Aust. J. Plant Physiol.* 19:519–538 — C4.
- Leuning (1995). *Plant Cell Environ.* 18:339–355 — semi-empirical $g_s$.
- Medlyn et al. (2011). *Glob. Change Biol.* 17:2134–2144 — USO $g_s$.
- Katul, Manzoni et al. (2010). *Ann. Bot.* 105:431–442 — optimization theory ($\lambda$).
- Manzoni et al. (2011). *Funct. Ecol.* 25:456–467 — $\lambda$ under water stress.
- Vico, Manzoni et al. (2013). *Agric. For. Meteorol.* 182–183:191–199.
- Sabot et al. (2022). *JAMES* 14(4):e2021MS002761 — two-limb (β\_stomata / β\_nonstomata).
- Zhou et al. (2013). *Agric. For. Meteorol.* 182–183:204–214 — combined stomatal + non-stomatal.
- Xu et al. (2016). *New Phytol.* 212:80–95 — the ED2 plant-hydraulics / Katul implementation.

## Code map

| Concept | Routine (`meds_leaf_gas_exchange`) |
|---|---|
| Rubisco kinetics $K_c$, $K_o$, $`\Gamma^*`$ in µmol mol⁻¹, eq. (3) | `rubisco_kinetics` |
| C3 / C4 demand | `assimilation_demand_c3`, `assimilation_demand_c4` |
| electron transport $J$ | `electron_transport_j` |
| Leuning / Medlyn / Katul-$\lambda$ | `stomata_gs_leuning`, `stomata_gs_medlyn`, `katul_lambda` |
| coupled solve | `solve_leaf_gas_exchange` |
| residuals | `ci_residual` → `residual_explicit_gs`, `residual_optimality` |
| root finder | `bisect_ci` (private) |
