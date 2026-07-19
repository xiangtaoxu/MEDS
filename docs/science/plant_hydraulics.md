# Plant hydraulics

MEDS resolves per-individual water transport as a small network of water pools (nodes) linked by
xylem/rhizosphere conductances, driven by transpiration at the top and soil water at the bottom. The
prognostic state is the node water potentials $\psi$ (leaf, wood), carried in the cohort state and
advanced each fast sub-step. Everything is in **MPa** (never metres of head), matching
`leaf_env_t%psi_leaf` and the `[pft].wstress_psi_*` traits. The physical reference is Xu et al. 2016
(X16), revised for coupled nodes, nonlinear pressure–volume, and Kirchhoff-integrated conductance.

Total head is $`\Psi_i = \psi_i + \rho g\,z_i`$ with $`\rho g = 9.804\times10^{-3}\ \mathrm{MPa\,m^{-1}}`$
(`grav_head`). Signs: $\psi \le 0$ (tension), $z$ measured upward.

## 1. The node network

Per-node water mass balance (kg H₂O per plant):

```math
C_i(\psi_i)\,\frac{d\psi_i}{dt} \;=\; \sum_{j\sim i} K_{ij}(\psi)\,\big(\psi_j - \psi_i + g_{ij}\big) \;-\; S_i \qquad(1)
```

with $C_i$ the **capacitance** [kg MPa⁻¹] (§2), $K_{ij}$ the **edge conductance** [kg s⁻¹ MPa⁻¹] (§3),
$`g_{ij}=\rho g\,(z_j-z_i)`$ the gravity offset, and $S_i$ a sink. The default topology is **2-node**
(leaf L + lumped wood W):

```math
C_L\,\dot\psi_L = K_{LW}\big(\psi_W-\psi_L+g_{WL}\big) - E, \qquad
C_W\,\dot\psi_W = K_{LW}\big(\psi_L-\psi_W-g_{WL}\big) + Q_{\text{root}} \qquad(2)
```

where $`g_{WL}=-\rho g\,H`$ (leaf sits a height $H$ above the wood datum), $E$ is transpiration
[kg s⁻¹], and $Q_{\text{root}}$ is root water uptake (§4). The matrix $A$ of the linear system
$C\dot\psi = A\psi + b$ is a grounded, symmetric negative-definite Laplacian ⇒ all eigenvalues of
$M=C^{-1}A$ are **real $\le 0$** — the closed-form solver (§5) is always valid.

## 2. Pressure–volume curves and capacitance

Tissue water storage follows the nonlinear **Bartlett / Tyree–Hammel** pressure–volume relation, per
tissue (leaf, wood). Traits: osmotic potential at full turgor $\pi_0<0$, bulk elastic modulus
$\varepsilon>0$, apoplastic fraction $a_f\in[0,1)$, saturated water content $w_{sat}$ [kg H₂O / kgC].

**Turgor loss point** and the symplastic RWC there:

```math
\psi_{tlp}=\frac{\pi_0\,\varepsilon}{\pi_0+\varepsilon}, \qquad R_{tlp}=\frac{\pi_0+\varepsilon}{\varepsilon}
```

**Potential ↔ symplastic RWC** ($R$):

```math
\psi(R)=\begin{cases}\varepsilon\,(R-R_{tlp})+\pi_0/R & R\ge R_{tlp}\ \text{(turgid)}\\[2pt]
\pi_0/R & R< R_{tlp}\ \text{(flaccid)}\end{cases}
\qquad
R(\psi)=\begin{cases}\dfrac{b+\sqrt{b^2-4\varepsilon\pi_0}}{2\varepsilon},\ b=\psi+\varepsilon+\pi_0 & \psi\ge\psi_{tlp}\\[6pt]
\pi_0/\psi & \psi<\psi_{tlp}\end{cases}
```

**Tissue water and capacitance** ($`W_{sat}=w_{sat}\cdot\text{biomass}`$; the apoplast is a constant
reservoir):

```math
W(\psi)=(1-a_f)\,W_{sat}\,R(\psi)+a_f\,W_{sat}, \qquad
C(\psi)=\frac{dW}{d\psi}=(1-a_f)\,W_{sat}\,\frac{dR}{d\psi}
```

```math
\frac{dR}{d\psi}=\begin{cases}\big(\varepsilon-\pi_0/R^2\big)^{-1} & \text{turgid}\\[2pt]
-R^2/\pi_0 & \text{flaccid}\end{cases}
```

A pure-elastic (X16-style) linear proxy capacitance is available for calibration:
$`C_{\text{lin}}=(1-a_f)\,w_{sat}/(\varepsilon+|\pi_0|)`$.

## 3. Xylem vulnerability and the Kirchhoff conductance law

The retained conductance fraction (1 − PLC) is a Weibull-like curve in $r=\psi/\psi_{50}\ge 0$
($\psi_{50}<0$ is the potential at 50 % loss; $a=$ `wood_kexp` the shape):

```math
k(\psi)=\frac{1}{1+(\psi/\psi_{50})^{a}}=\frac{1}{1+r^{a}}
```

**Vulnerability enters through the Kirchhoff (matric flux) potential — never pointwise.** Define

```math
\Phi(\psi)=\int_0^{\psi} k(s)\,ds=\psi_{50}\!\int_0^{r}\frac{du}{1+u^{a}}
=\begin{cases}\psi_{50}\ln(1+r) & a=1\\ \psi_{50}\arctan(r) & a=2\\
\psi_{50}\cdot\text{(7-pt Gauss–Legendre)} & \text{else}\end{cases}
```

The **edge conductance** over a finite drop is the $\psi$-averaged retained fraction — the
finite-difference-consistent conductance that reproduces the exact Kirchhoff-integrated flux:

```math
K_{ij}=k_{\text{cond}}\,\frac{\Phi(\psi_{up})-\Phi(\psi_{down})}{\psi_{up}-\psi_{down}}
\;\xrightarrow[\Delta\psi\to0]{}\; k_{\text{cond}}\,k(\psi)
```

$k_{\text{cond}}$ is the maximum (fully-hydrated) per-plant conductance: whole-plant
$`k_{\text{cond}}=k_{plant\_max}\cdot\text{leaf area}`$ (default), or segment
$`k_{\text{cond}}=w_{kmax}\cdot A_{sap}/(H\cdot\text{vessel\_curl})`$ from stem allometry (Huber value
$`H_v=A_{sap}/A_{leaf}`$). For general $`a\notin\{1,2\}`$ the integral is precomputed once into a fixed
uniform-grid **lookup table** $G(r)$ and read by linear interpolation on the hot path (the closed
forms are kept for $`a\in\{1,2\}`$); the table stores $r$-normalized $G$, so $\psi_{50}$ is a runtime
scale.

## 4. Root water uptake

$Q_{\text{root}}$ enters the wood node from the soil. The **current** implementation uses a single
prescribed rhizosphere boundary, $`Q_{\text{root}}=\text{rhizo\_cond}\,(\psi_{soil}-\psi_W)`$.

The **multi-layer** formulation (opt-in `[hydraulics].multilayer_roots`; ED2-faithful; see
`MEDS_MULTILAYER_ROOTS_DESIGN.md`) couples to the prognostic soil column — per-layer soil ψ and the
**unsaturated** conductivity $K_{soil}(k)=K(\theta_k)$ (`soil_hydr_cond`) — and sums the soil layers the
roots reach, in parallel to the common wood node. Per-layer root+rhizosphere conductance
(Katul 2003; MEDS `rhizosphere_cond`, = ED2 `gw_cond`):

```math
g_k=\frac{K_{soil}(k)\,\sqrt{\text{RAI}_k}}{\pi\,\Delta z_k}\cdot\frac{1}{n_{plant}},\qquad
\text{RAI}_k=b_{root}\cdot\text{SRA}\cdot\text{root\_frac}(k)\cdot n_{plant}
```

with the ED2 cumulative-exponential root profile
$`\text{root\_frac}(k)=\beta^{\,d_{k-1}/D}-\beta^{\,d_k/D}`$ ($\beta=$ `root_beta`, $D$ the max rooting
depth). The parallel network collapses to an effective boundary,

```math
G_{\text{root}}=\sum_k g_k, \qquad
\psi_{\text{soil,eff}}=\frac{\sum_k g_k\,(\psi_{soil,k}+\rho g\,z_k)}{\sum_k g_k},\qquad
Q_{\text{root}}=G_{\text{root}}\,(\psi_{\text{soil,eff}}-\psi_W)
```

so the 2-node solver is unchanged; the converged total is distributed back per layer for the soil sink,
$`U_k = \text{supply}_k^{+}/\sum_j \text{supply}_j^{+}\cdot Q_{\text{root}}`$ with
$`\text{supply}_k^{+}=\max\!\big(g_k(\psi_{soil,k}-\psi_W+\rho g z_k),\,0\big)`$. **Hydraulic
redistribution is not enabled:** a dry layer that would give a negative supply (root→soil efflux) is
floored to 0, so every $U_k\ge 0$ and $`\sum_k U_k = Q_{\text{root}}`$. HR — the per-layer efflux and the
soil re-wetting it implies — is deferred to a future version (see
`docs/dev_plans/MEDS_MULTILAYER_ROOTS_DESIGN.md`).

## 5. The solver

Each fast step is integrated by freezing the linear system at the current state
($\dot\psi=M\psi+c$, coefficients from capacitance, the Kirchhoff edge conductance, and the soil
boundary) and advancing it **exactly** over a sub-step with the underflow-safe 2×2 matrix exponential
(sinh-c form; eigenvalues real $\le0$): $`\psi(h)=\psi^* + e^{Mh}(\psi_0-\psi^*)`$, where $`\psi^*`$ is
the Ohm's-law steady state $-M^{-1}c$. Adaptive **step-doubling** controls the sub-step via the shared
embedded-error controller (`meds_numerics%adaptive_step_update`); the explicit RHS
(`plant_water_tendency`) feeds the IMEX-ARK integrator. Boundary fluxes (sapflow, root uptake) close a
machine-precision water budget from the converged storage change $\Delta W$.

## Parameters (config names, `[hydraulics]` / `hydraulics_config_t`)

| Symbol | Config key | Meaning |
|---|---|---|
| $\pi_{0}$ | `leaf_pi0`, `wood_pi0` | osmotic potential at full turgor [MPa] |
| $\varepsilon$ | `leaf_elastic_mod`, `wood_elastic_mod` | bulk elastic modulus [MPa] |
| $a_f$ | `leaf_apoplast_frac`, `wood_apoplast_frac` | apoplastic fraction [–] |
| $w_{sat}$ | `leaf_water_sat`, `wood_water_sat` | saturated water content [kg H₂O / kgC] |
| $\psi_{50}$ | `wood_psi50` | xylem potential at 50 % loss of conductance [MPa] |
| $a$ | `wood_kexp` | vulnerability-curve shape [–] |
| $`k_{plant\_max}`$ | `k_plant_max` | max whole-plant conductance [kg s⁻¹ MPa⁻¹ m⁻²_leaf] |
| $K_s$ | `wood_kmax` | sapwood specific conductivity (segment mode) [kg m⁻¹ s⁻¹ MPa⁻¹] |
| — | `vessel_curl` | tortuosity / path-length factor [–] |
| — | `rhizo_cond` | rhizosphere conductance (single-BC) [kg s⁻¹ MPa⁻¹] |
| $\beta$ | `root_beta` | ED2 root-profile decay [–] (feeds `root_fraction_profile`) |
| $D$ | `root_depth` | maximum rooting depth [m] |
| SRA | `specific_root_area` | specific root area [m² kgC⁻¹] (multi-layer rhizosphere conductance) |
| — | `multilayer_roots` | opt-in soil↔plant per-layer coupling [bool, default `false`] |

## References
- Xu, Medvigy, Powers, Becknell & Guan (2016), *New Phytologist* — X16 hydraulics.
- Bartlett, Scoffoni & Sack (2012); Tyree & Hammel (1972) — pressure–volume theory.
- Katul, Leuning & Oren (2003) — rhizosphere conductance.
- ED2 `ED/src/dynamics/plant_hydro.f90`; `docs/dev_plans/MEDS_HYDRAULICS_DESIGN.md` (§4 governing
  equations, §16 per-layer roots); `docs/dev_plans/MEDS_MULTILAYER_ROOTS_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| PV curves / capacitance | `meds_hydro_curve`: `pv_psi_tlp`, `pv_rwc_tlp`, `psi_from_rwc`, `rwc_from_psi`, `water_content`, `capacitance` |
| vulnerability + Kirchhoff | `plc_retained`, `flux_potential`, `kirchhoff_edge` (+ table: `build_hydro_table`, `flux_potential_lin`, `kirchhoff_edge_tab`) |
| quadrature / root-find | `meds_numerics`: `gauss_legendre_7`, `bisect_root` |
| multi-layer root boundary | `meds_plant_hydraulics`: `rhizosphere_cond`, `root_fraction_profile`, `effective_root_boundary` |
| network solver | `meds_plant_hydraulics`: `solve_plant_water` (`freeze_coeffs` + `advance_exact_linear` + `exact_substep`), `plant_water_tendency` |
| config flatten / soil coupling | `meds_column_dynamics`: `apply_hydraulics_config`; opt-in per-layer soil↔plant in `column_fast_step` (`soil_hydr_cond` → K(θ)) |
