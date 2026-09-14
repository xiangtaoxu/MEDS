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
A = \frac{g_s}{1.6}\,(C_s - C_i), \qquad C_s = C_a - 1.4\,\frac{A}{g_b}
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
A_p = 3\,\mathrm{TPU}
```

where $A_c$ is Rubisco-limited, $A_j$ RuBP/light-limited, $A_p$ triose-phosphate-limited;
$`\Gamma^*`$ is the CO₂ compensation point without respiration and $K_c, K_o$ the Rubisco
Michaelis constants. The electron-transport rate $J$ follows the non-rectangular hyperbola
$`\theta J^2 - (I_2 + J_{max})J + I_2 J_{max}=0`$ (smaller root), with
$`I_2 = \tfrac12\,\phi_{\text{PSII}}\,\alpha_{\text{leaf}}\,\mathrm{PAR}`$.

**$`\Gamma^*`$ is proportional to the oxygen partial pressure**, because it is set by Rubisco's
specificity for CO₂ over O₂: $`\Gamma^* = \tfrac{1}{2} O / S_{c/o}`$. MEDS therefore scales it off
the O₂ its reference value was measured at,

```math
\Gamma^*(T, O) = \Gamma^*_{25}\,f_{\mathrm{Arr}}(T)\;\frac{O}{O_\mathrm{ref}},
\qquad O_\mathrm{ref} = 0.209 \qquad(1)
```

so `o2_mol_frac` reaches **both** places oxygen enters the C3 demand — the Michaelis term
$`K_c(1 + O/K_o)`$ and the compensation point. The shipped `gstar25` = 4.275 Pa is Bernacchi et al.
(2001), measured at 21 % O₂, so the factor is exactly 1 in any default run. It matters for
palaeo-atmosphere work (Carboniferous/Permian O₂ is reconstructed near 30–35 %) and for the low-O₂
sensitivity experiment that is the standard way to isolate photorespiration: measured on a
midday tropical leaf, halving O₂ raises $`A_\mathrm{net}`$ by 22 % and 35 % O₂ cuts it by 24 %.

**C4 (Collatz et al. 1992):** $A_c = V_{cmax}$, $A_j$ = a light-limited slope, $A_p = k_p C_i$
(PEP-case CO₂ limitation), with $`\Gamma^* \approx 0`$ (the CO₂-concentrating mechanism suppresses
photorespiration).

The three rates are combined into $A_g$ either as a sharp $\min(A_c,A_j,A_p)$ (`COLIM_MIN`) or as
two nested smoothing quadratics (`COLIM_QUADRATIC`, the same `quadratic_smaller_root` used for $J$).
Then $A = A_g - R_d$.

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

C3 used to pass `theta_j` into **both** of its smoothings while C4 already had its own pair. Because
$A_c$ and $A_j$ sit close together at ambient CO₂, that penalty landed exactly where the model spends
most of its time. Measured at PFT-1 kinetics, $C_i = 280$ µmol mol⁻¹, saturating light:

| curvature used | $A_g$ | vs $\min(A_c,A_j,A_p) = 14.40$ |
|---|---|---|
| `theta_j` = 0.85 (before) | 10.80 | **−25.0 %** |
| `theta_cj_c3` = 0.98, `theta_ip_c3` = 0.95 | 13.22 | −8.2 % |

The shortfall was nearly independent of $V_{cmax}$ — 30, 31, 31, 32 % at $V_{cmax,25}$ = 60, 90, 120,
150 — so it was a systematic offset rather than a regime effect, and **no measured $V_{cmax}$
reproduced a measured assimilation rate**. The shipped $V_{cmax,25}$ presets were checked against
this and left alone: they were set when the leaf module first landed and never revised, they sit
mid-range for their PFT descriptions, and the only other GPP-facing knob (`gpp_ref`) is the stub used
when the fast loop is off. The smoothing loss was therefore never compensated anywhere — so removing
it is a correction, not the unwinding of a hidden calibration.

Kernels: `assimilation_demand_c3`, `assimilation_demand_c4`, `electron_transport_j`, `combine_limits`.

---

### Thermal acclimation (opt-in)

With `[leaf_physiology].thermal_acclimation = true`, the peaked form's entropy term and the capacity
ratio follow the growth temperature (Kattge & Knorr 2007):

```math
\Delta S_{V} = 668.39 - 1.07\,T_g, \qquad
\Delta S_{J} = 659.70 - 0.75\,T_g, \qquad
\frac{J_{max,25}}{V_{cmax,25}} = 2.59 - 0.035\,T_g \qquad(2)
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
- **It is site-level, not per cohort**, for the same reason. Per-cohort acclimation is a different
  question — it would need a relation calibrated on leaf rather than air temperature — and is on the
  ROADMAP.

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
g_s = g_0 + \frac{g_1\,A}{(C_s - \Gamma^*)\,(1 + D/D_0)}
```

A semi-empirical law: conductance tracks assimilation, discounted by CO₂ drawdown
$`(C_s-\Gamma^*)`$ and by a hyperbolic humidity response $(1 + D/D_0)$. Parameters: `stomatal_g0`
($g_0$), `stomatal_g1` ($g_1$), `stomatal_d0` ($D_0$). Code: `stomata_gs_leuning`.

### 3.2 Medlyn et al. (2011) — unified stomatal optimization (USO)

```math
g_s = g_0 + 1.6\left(1 + \frac{g_1}{\sqrt{D}}\right)\frac{A}{C_s}, \qquad D \text{ in kPa}
```

Derived from the same $\max(A-\lambda E)$ optimization as Katul but solved analytically under a
linearized $A(C_i)$, giving a closed form. Here $g_1 \propto \lambda^{-1/2}$ — a key relation used
by the water-stress limbs (§4). Parameters: `stomatal_g0`, `stomatal_g1`. Code: `stomata_gs_medlyn`.

### 3.3 Katul et al. (2010) — marginal-WUE optimality

Stomata are assumed to **maximize carbon gain net of a water cost**:

```math
\max_{g_s}\bigl(A - \lambda E\bigr)
```

where $E = g_s D$ is transpiration and $\lambda$ is the **marginal water-use efficiency** — the
"price" of water (mol C per mol H₂O). The optimum satisfies $\partial A/\partial E = \lambda$.
Writing both fluxes in $C_i$-space (with $A = \tfrac{g_s}{1.6}(C_s - C_i)$, so
$`g_s = 1.6\,A/(C_s-C_i)`$ and $`E = 1.6\,D\,A/(C_s-C_i)`$) and setting
$`dA/dC_i = \lambda\,dE/dC_i`$ rearranges to the residual root-found on $C_i$:

```math
A'\,(C_s - C_i)^2 \;=\; 1.6\,D\,\lambda\,\bigl(A'(C_s - C_i) + A\bigr), \qquad A' \equiv \frac{dA}{dC_i}
```

$A'$ is taken by central difference (the co-limited $A(C_i)$ has no clean analytic slope). Unlike
the explicit models, $g_s$ is recovered *after* solving, as $`g_s = 1.6\,A/(C_s-C_i)`$; if that
would fall below the cuticular floor $g_0$, the solver re-solves g0-pinned (§5). Parameters:
`katul_lambda25` ($\lambda_{25}$), and the water-stress terms below. Code: `residual_optimality`,
`katul_lambda`.

---

## 4. Water stress — two limbs

MEDS follows the two-factor decomposition of Sabot et al. (2022) / Zhou et al. (2013): drought
acts on **capacity** (biochemistry) and on the **stomatal** aperture independently. Both are
per-PFT and tunable.

### 4.1 Non-stomatal (capacity) limb — $\beta_{ns}$

A linear ramp in **leaf** water potential, applied to $V_{cmax}$, $J_{max}$, and TPU for **all**
stomatal models (capacity downregulation is a biochemistry effect, scheme-independent):

```math
\beta_{ns} = \mathrm{clamp}\!\left(\frac{\psi_{leaf} - \psi_{close}}{\psi_{open} - \psi_{close}},\,0,\,1\right),
\qquad \{V_{cmax}, J_{max}, \mathrm{TPU}\} \mathrel{*}= \beta_{ns}
```

Parameters: `wstress_psi_open` ($\psi_{open}$, $\beta_{ns}=1$), `wstress_psi_close`
($\psi_{close}$, $\beta_{ns}=0$).

### 4.2 Stomatal limb — $\beta_s$ (Sabot 2022)

An exponential in **soil / predawn** water potential, downregulating the *aperture*:

```math
\beta_s = \min\!\bigl(1,\ \exp(s_{ref}\,\psi)\bigr)
```

- **Leuning / Medlyn:** scale the slope, $`g_{1,\text{eff}} = g_1\,\beta_s`$.
- **Katul:** raise the marginal price, $`\lambda = \lambda_{25}\,\beta_s^{-e}`$ (`katul_lambda`).
  Because Medlyn gives $g_1 \propto \lambda^{-1/2}$, the exponent $e = 2$ makes the Katul and
  Leuning/Medlyn responses *identical* for the same $\beta_s$.

Parameters: `wstress_sref_stomata` ($s_{ref}$, ~2 MPa⁻¹), `wstress_lambda_exp` ($e$).

> **Note.** The stomatal driver is the *supply-side* potential $\psi$ (via `leaf_env_t%psi`,
> default 0 = well-watered); the capacity driver is midday $\psi_{leaf}$, the *demand-side*
> tension. Sabot writes $\psi$ as soil/predawn potential; MEDS supplies the cohort's own **predawn
> (previous-day daily-max) leaf** potential — see below. The field was called `psi_soil` until
> issue #99; it was renamed because it never carried a soil potential. This mirrors ED2's Katul
> `stoma_beta` (with $`s_{ref}\,e \equiv`$ ED2's `stoma_beta`, and $\psi \equiv$ ED2's
> `dmax_leaf_psi` — the same quantity, same name); the divergences are tracked in **issue #47**.

> **What actually drives $`\beta_s`$ (issue #95).** `psi` is an *optional* argument, and until
> recently the fast-loop driver never passed it — so $`\psi=0`$ and $`\beta_s\equiv 1`$: there
> was **no stomatal water stress at all**. Over 8 dry days GPP stayed flat at ~98.5 µmol m⁻² s⁻¹ while
> the wood store sat empty. It is now fed the **previous day's daily-maximum leaf water potential** —
> the model's predawn potential, which is what the field measurement behind $`s_{ref}`$ actually is,
> and the plant's own integration over its rooted profile (so it already carries rooting depth,
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
> equilibrates overnight — which is why the kernel field is named plainly `psi` (issue #99).

### 4.3 Where this diverges from ED2 — deliberately (issue #47)

MEDS uses the Sabot two-limb scheme; ED2's `farq_katul.f90` uses a Manzoni-style single
`stoma_beta` plus a turgor-loss capacity term. **The decision is to keep Sabot and state the
difference**, rather than adopt ED2's forms. Three divergences, and why each one stands:

| | MEDS | ED2 (`farq_katul.f90`) |
|---|---|---|
| **capacity shape** | linear ramp in $\psi_{leaf}$ between $\psi_{open}$ and $\psi_{close}$ | $`1/\bigl(1+0.1(\psi_{leaf}/\psi_{tlp})^{6}\bigr)`$, a 6th-power turgor-loss curve |
| **capacity targets** | $V_{cmax}$, $J_{max}$ **and** TPU | $J_{max}$ and quantum yield $\alpha$ only ($V_{cmax}$ and TPU are commented out) |
| **stomatal scope** | Leuning/Medlyn $g_1$ **and** Katul $\lambda$ | Katul $\lambda$ only |
| **stomatal driver** | predawn (previous-day daily-max) leaf $\psi$ | `dmax_leaf_psi` — the same quantity |
| **parameter split** | $s_{ref}$ (shared with the $g_1$ limb) $\times$ $e$ (Katul-only) | one `stoma_beta` |

**Capacity shape and targets.** ED2's turgor-loss form is anchored on a measurable trait,
$\psi_{tlp}$, which is an advantage; but it is applied to $J_{max}$ and $\alpha$ while leaving
$V_{cmax}$ untouched, which makes the capacity limb act only on the light-limited branch. The Sabot
formulation downregulates all three potential rates together, which is what the underlying
observations describe, and it keeps the limb's effect independent of which limitation a leaf happens
to be under. MEDS also applies it for **every** stomatal model, not just Katul, because capacity
downregulation is biochemistry and should not depend on the aperture scheme selected.

The cost of the linear ramp is that it is a **linear amplifier on a $\psi_{leaf}$ that is itself
`dt_fast`-dependent** — at 900 s a 1 MPa $\psi_{leaf}$ error became a 33 % GPP shift. That is why
`wstress_nonstomata` ships **off** and the shipped configs say so. ED2's 6th power would amplify the
same error harder still, not less; this is an argument about the input, not the curve.

**Stomatal scope.** Restricting $\beta_s$ to Katul, as ED2 does, would leave the two explicit
schemes with no stomatal drought response at all. MEDS applies it to both, and the exponent
$e = 2$ makes them exactly equivalent because Medlyn gives $`g_1 \propto \lambda^{-1/2}`$ (§4.2).

**Parameter split.** MEDS's $s_{ref}\,e$ is precisely ED2's `stoma_beta` — the same exponential in
the same potential — so a value calibrated for one transfers to the other as
$`\texttt{stoma\_beta} = -s_{ref}\,e`$. For Katul alone only the product matters, so the split is
degenerate there; it earns its keep because $s_{ref}$ *alone* drives the Leuning/Medlyn limb, which
ED2 has no equivalent of.

> **The `psi_soil` wiring this issue also filed is already closed.** #47 recorded that the driver
> left the stomatal potential at 0 while `soil_psi_root` was computed only after the leaf call. That
> was fixed under issue #95: the field is `psi`, the driver passes `dmax_psi_leaf` on both call
> paths, and unset cohorts are seeded from the surface-layer soil potential. See the #95 note in
> §4.2 for what it was worth — $`\beta_s \equiv 1`$, i.e. no stomatal water stress at all.

### Arrestors: stopping a plant that has run out of water

$`\beta_s`$ scales $`g_1`$ only, so as it goes to zero the conductance falls to the residual
$`g_0`$ and never reaches it. Measured: ~2.6 mm day⁻¹ still leaving a plant whose wood store was
empty. An *arrestor* is therefore needed on top of the $`\beta_s`$ ramp. `[run].leaf_stress_arrestor`
selects it.

**`ARREST_GS_CLAMP` (default).** Below twice the turgor-loss point
$`\psi_{tlp}=\pi_0\varepsilon/(\pi_0+\varepsilon)`$ — the same PV curve the hydraulics solver uses —
the stomata shut completely: $`A_g=0`$, $`A_n=-R_d`$, $`g_s=0`$, $`E=0`$. It latches on the
*daily-max* potential, so it is a once-a-day decision on a slow integrated measure rather than a
per-step switch on a noisy sub-daily $`\psi`$.

Measured on the 8-day dry case ($`\theta=0.12`$, no rain), converged over six cadences:

| | wood store, day 8 | cumulative ET | error at $`dt_{fast}`$ = 900 s |
|---|---|---|---|
| `ARREST_NONE` | 1.229 kg plant⁻¹ | 27.5 mm | **54.9 %** |
| `ARREST_GS_CLAMP` | 1.940 kg plant⁻¹ | 15.8 mm | **2.3 %** |

The last column is the load-bearing one. Without an arrestor the drought trajectory converges only
*slowly*, because nothing bounds the drawdown until the state reaches a numerical floor — so the
day-8 answer at production cadence is 55 % wrong. The closure removes **96 % of that cadence error at
zero recurring cost**, which is why the structural fix proposed in issue #93 (making soil water
prognostic inside the ARK stages, at +14–26 % on every step of every run) was closed: an exact uptake
seam could at best address the residual 2.3 %.

The closure also bounds a defect it was not designed for. The frozen-uptake seam (issue #92) leaves
the wood store with a `dt`-dependent deficit that, in dry soil and **without** an arrestor, grows
without limit — −4.2 % on day 1 to **−65.5 %** by day 10 at 900 s, with $`\psi_{wood}`$ pinned on the
−10000 MPa clamp from day 5. With `ARREST_GS_CLAMP` it plateaus at **+2.3 %** from day 6 onward and
$`\psi_{wood}`$ recovers to −5.0 MPa. In moist soil the deficit is negligible either way (2e-4 %,
decaying) because the overnight wood↔soil relaxation resets it.

One caveat: days 1–2 are identical under both settings, because the feedback rolls over **once a
day** and cannot respond until the day after the potential drops. Drought *onset* therefore carries a
~9 % wood-water error at 900 s regardless of arrestor, relaxing to 2.3 % once the closure engages.

> ⚠️ **Earlier figures here were wrong.** This section previously quoted "the wood store refills from
> 0 to 2.29 kg plant⁻¹ and transpiration collapses from ~7 to 0.013 mm day⁻¹". Those came from a probe
> that never set `aenv%theta_atm`, leaving it at its 298.15 K default while the forcing drove the
> canopy over 282–294 K — Monin–Obukhov then saw a permanent stable layer and floored `ustar`,
> suppressing turbulent exchange ~44×. The conclusion held; the magnitudes were inflated. Fixed in
> issue #97, which also found that **no column test set `theta_atm` either**.

**Dynamic vapour pressure — built, measured, REMOVED (issue #96).** The substomatal air is in
equilibrium with leaf water at $`\psi_{leaf}`$, so its humidity is the Kelvin value
$`\mathrm{RH}=\exp\!\big(\psi/(\rho_w R_v T)\big)`$, not 1. Transpiration would then use
$`e_i=\mathrm{RH}\cdot e_{sat}`$, shrinking the gradient continuously as the leaf dries and
**reversing** it once $`e_i<e_a`$ — i.e. foliar water uptake — which arrests transpiration with no
threshold parameter at all, at $`\psi=\rho_w R_v T\ln(\mathrm{RH}_{air})`$. MEDS already applies this
exact relation to **soil** evaporation (`ground_evaporation`'s `alpha_soil`), so the model presently
assumes sub-saturated vapour over soil at −10 MPa but saturated vapour over leaf water at −100 MPa.

It was implemented and **removed**: it produced NaN by day 2 of the dry test. The physics is sound;
the wiring was not. `le_ref` in `surface_derivs` is the latent flux *linearised about* $`T_{cas}`$,
with the leaf–air temperature difference carried by `le_slope`·$`\Delta T_l`$, so scaling both terms
by RH moves the base point and the slope inconsistently — the symptom was −1.24 mm day⁻¹ of
*condensation* at midday on a nearly-turgid leaf. Flooring RH did **not** help, confirming the fault
is the linearisation rather than the magnitude. Doing it properly means rebuilding the flux *and* its
temperature derivative from $`\mathrm{RH}\cdot q_{sat}(T_{leaf})`$ inside the leaf energy balance.
Its likely future value is as a route to **foliar water uptake**, not as a stress arrestor.

> Note the Kelvin correction is deliberately kept out of the **stomatal** model: Medlyn/Leuning
> $`g_1`$ were calibrated against $`e_{sat}`$-based VPD, so feeding them a corrected VPD would
> silently re-tune $`g_1`$.

---

## 5. The coupled solver

`solve_leaf_gas_exchange` brackets $`C_i \in (\Gamma^*, C_a]`$ and bisects the residual to a
tolerance `ci_tol_ppm` (shared `meds_numerics%bisect_root`):

- **Night / closed branch** — if the best-case net rate $A(C_a) \le 0$, stomata sit at $g_0$ and
  the solve returns immediately.
- **Two residual providers** — `residual_explicit_gs` (Leuning/Medlyn diffusion identity, and the
  g0-pinned fallback with $g_s := g_0$) and `residual_optimality` (Katul, §3.3). Both are `pure`.
- **`force_g0` fallback** — if the model yields no sign change (e.g. Katul under strong stress
  wanting $g_s < g_0$), the solver retries g0-pinned, which always brackets when $A(C_a) > 0$.

Then $g_s$, $C_s$, $C_i$, and $E$ are assembled consistently. See
[`docs/dev_plans`](../dev_plans) for the design history if you need the *why* behind the structure.

---

## Parameters (config names)

| Symbol | Config key (`[pft]`) | Meaning |
|---|---|---|
| $g_0$ | `stomatal_g0` | cuticular / residual conductance |
| $g_1$ | `stomatal_g1` | stomatal slope (Leuning `--` / Medlyn kPa$^{0.5}$) |
| $D_0$ | `stomatal_d0` | Leuning humidity sensitivity |
| $\lambda_{25}$ | `katul_lambda25` | Katul marginal WUE (well-watered) |
| $`\psi_{open},\psi_{close}`$ | `wstress_psi_open`, `wstress_psi_close` | capacity-limb ramp bounds |
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
| C3 / C4 demand | `assimilation_demand_c3`, `assimilation_demand_c4` |
| electron transport $J$ | `electron_transport_j` |
| Leuning / Medlyn / Katul-$\lambda$ | `stomata_gs_leuning`, `stomata_gs_medlyn`, `katul_lambda` |
| coupled solve | `solve_leaf_gas_exchange` |
| residuals | `residual_explicit_gs`, `residual_optimality` |
| root finder | `meds_numerics%bisect_root` |
