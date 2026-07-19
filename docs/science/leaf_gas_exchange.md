# Leaf gas exchange

How MEDS computes leaf-level photosynthesis, stomatal conductance, and their coupling. The code
lives in [`src/plant/meds_leaf_gas_exchange.f90`](../../src/plant/meds_leaf_gas_exchange.f90); the
public entry point is `leaf_gas_exchange` in `meds_plant_interface`.

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

**C4 (Collatz et al. 1992):** $A_c = V_{cmax}$, $A_j$ = a light-limited slope, $A_p = k_p C_i$
(PEP-case CO₂ limitation), with $`\Gamma^* \approx 0`$ (the CO₂-concentrating mechanism suppresses
photorespiration).

The three rates are combined into $A_g$ either as a sharp $\min(A_c,A_j,A_p)$ (`COLIM_MIN`) or as
two nested smoothing quadratics (`COLIM_QUADRATIC`, the same `quadratic_smaller_root` used for $J$).
Then $A = A_g - R_d$.

Kernels: `assimilation_demand_c3`, `assimilation_demand_c4`, `electron_transport_j`, `combine_limits`.

---

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
\beta_s = \min\!\bigl(1,\ \exp(s_{ref}\,\psi_{soil})\bigr)
```

- **Leuning / Medlyn:** scale the slope, $`g_{1,\text{eff}} = g_1\,\beta_s`$.
- **Katul:** raise the marginal price, $`\lambda = \lambda_{25}\,\beta_s^{-e}`$ (`katul_lambda`).
  Because Medlyn gives $g_1 \propto \lambda^{-1/2}$, the exponent $e = 2$ makes the Katul and
  Leuning/Medlyn responses *identical* for the same $\beta_s$.

Parameters: `wstress_sref_stomata` ($s_{ref}$, ~2 MPa⁻¹), `wstress_lambda_exp` ($e$).

> **Note.** The stomatal driver is $\psi_{soil}$ (via `leaf_env_t%psi_soil`, default 0 =
> well-watered); the capacity driver is midday $\psi_{leaf}$. This mirrors ED2's Katul
> `stoma_beta` (with $`s_{ref}\,e \equiv`$ ED2's `stoma_beta`, $\psi_{soil}\approx$ ED2's
> `dmax_leaf_psi`); the divergences are tracked in **issue #47**.

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
