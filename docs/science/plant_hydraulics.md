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

$k_{\text{cond}}$ is the maximum (fully-hydrated) per-plant conductance. `[hydraulics].conductance`
picks its form: `"whole_plant"` (default), $`k_{\text{cond}}=k_{plant\_max}\cdot\text{leaf area}`$, or
`"segment"`, $`k_{\text{cond}}=w_{kmax}\cdot A_{sap}/(H\cdot\text{vessel\_curl})`$ from stem allometry
(Huber value $`H_v=A_{sap}/A_{leaf}`$). `wood_kmax` and `vessel_curl` are read only in segment
mode, and `k_plant_max` only in whole-plant mode. Before v0.3.1 there was no key for the mode, so
`wood_kmax` and `vessel_curl` were accepted and never used. For general $`a\notin\{1,2\}`$ the integral is precomputed once into a fixed
uniform-grid **lookup table** $G(r)$ and read by linear interpolation on the hot path (the closed
forms are kept for $`a\in\{1,2\}`$); the table stores $r$-normalized $G$, so $\psi_{50}$ is a runtime
scale.

## 4. Root water uptake

$Q_{\text{root}}$ enters the wood node from the soil, through a **per-layer** root boundary. The
single root-fraction-weighted rhizosphere boundary that used to be the default is **gone**, along with
its `[hydraulics].multilayer_roots` switch — the per-layer path is the only one.

The **multi-layer** formulation (unconditional since the switch was deleted; ED2-faithful; see
`MEDS_MULTILAYER_ROOTS_DESIGN.md`) couples to the prognostic soil column — per-layer soil ψ and the
**unsaturated** conductivity $K_{soil}(k)=K(\theta_k)$ (`soil_hydr_cond_from_theta`) — and sums the soil layers the
roots reach, in parallel to the common wood node.

**Each cohort has its own roots.** It roots to a depth set by its height (ED2's IALLOM 1 allometry,
Christoffersen 2013), capped at the soil column's bottom $z_{\text{bot}}$:

```math
D=\min\!\left(b_{1Rd}\,h^{\,b_{2Rd}},\ z_{\text{bot}}\right)
```

with $`b_{1Rd}`$ = `root_depth_b1` (1.114 m) and $`b_{2Rd}`$ = `root_depth_b2` (0.4223), which give 5 m
for a 35 m tree and 1.8 m for a 3 m sapling. Within $D$ its fine roots fall off exponentially, as in
ED2's `distrib_root`: the share above depth $d$ is

```math
Y(d)=\frac{1-\beta^{\,\min(d,D)/D}}{1-\beta},
```

so a layer between $d_{k-1}$ and $d_k$ holds $`Y(d_k)-Y(d_{k-1})`$, the shares sum to one, and a layer
below $D$ holds none ($\beta$ = `root_beta`, 0.1, the share the profile would leave below $D$ if it went
on; ED2's default). The patch's profile, which weights the root-zone temperature of fine-root
respiration and spreads the uptake when no layer supplies any, is its cohorts' profiles weighted by
their fine-root carbon.

**Soil to root, layer by layer: a single root** (Gardner 1960; Cowan 1965). Water flows radially to a
root of radius $`r_{\text{root}}`$ from a soil cylinder whose radius is half the distance between roots,
so per plant

```math
g_k=\frac{2\pi\,K_{soil}(k)\,L_k}{\ln\!\left(r_{\text{half},k}/r_{\text{root}}\right)},\qquad
L_k=b_{root}\cdot\text{SRL}\cdot Y_k,\qquad
r_{\text{half},k}=\left(\pi\,\frac{\sum_{\text{cohorts}} n_{plant}\,L_k}{\Delta z_k}\right)^{-1/2}
```

where $L_k$ is the plant's fine-root length in the layer (SRL = `specific_root_length`, m kgC⁻¹),
$`r_{\text{half},k}`$ comes from the patch's root-length density in the layer, and $`r_{\text{root}}`$ =
`fine_root_radius` (the logarithm is floored at ln 2 for very crowded roots). The conductance grows with
root length and falls only by the logarithm as roots crowd, so it does not depend on how the soil is
layered or how a stand is split into cohorts. ED2-hydro's per-layer form, $`K\sqrt{\text{RAI}_k}/(\pi\,
\Delta z_k)`$ after Katul et al. (2003), depends on both, and connects the thick deep layers one to two
orders of magnitude more weakly (#375). Whether the absorbing length really grows in proportion to fine-root carbon is open
(#377). The parallel network collapses to an effective boundary,

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
`docs/dev_plans/archive/MEDS_MULTILAYER_ROOTS_DESIGN.md`).

## 5. The solver

Each fast step is integrated by freezing the linear system at the current state
($\dot\psi=M\psi+c$, coefficients from capacitance, the Kirchhoff edge conductance, and the soil
boundary) and advancing it **exactly** over a sub-step with the underflow-safe 2×2 matrix exponential
(sinh-c form; eigenvalues real $\le0$): $`\psi(h)=\psi^* + e^{Mh}(\psi_0-\psi^*)`$, where $`\psi^*`$ is
the Ohm's-law steady state $-M^{-1}c$. Adaptive **step-doubling** controls the sub-step via the shared
embedded-error controller (`meds_numerics%adaptive_step_update`); the explicit RHS
(`advance_water_mass_full`, with the transpiration corrector) feeds the ESDIRK2 (`ark`) integrator. Boundary fluxes (sapflow, root uptake) close a
machine-precision water budget from the converged storage change $\Delta W$.

## Parameters (config names, `[hydraulics]` / `hydraulics_config_t`)

> **These are per-PFT** (#179). The `[hydraulics]` block is the shared base; any of the thirteen
> traits below can be given a per-PFT array in the `[pft]` table under the **same key name**, and
> that value wins for that PFT. An absent key falls back to the `[hydraulics]` scalar, so a config
> can make **one** trait per-PFT without restating the other twelve, and a config that mentions none
> is byte-identical to before.
>
> `apply_hydraulics_config` builds the per-PFT table and there is deliberately **no** PFT-uniform
> companion struct: a second, easier-to-reach copy is how a caller ends up silently running every
> PFT on the first one's hydraulics. Each entry carries its **own** Kirchhoff lookup, rebuilt from
> that PFT's `wood_kexp` — sharing one would have given every PFT the first one's vulnerability
> shape while every budget still closed.

| Symbol | Config key | Meaning |
|---|---|---|
| $\pi_{0}$ | `leaf_pi0`, `wood_pi0` | osmotic potential at full turgor [MPa] |
| $\varepsilon$ | `leaf_elastic_mod`, `wood_elastic_mod` | bulk elastic modulus [MPa] |
| $a_f$ | `leaf_apoplast_frac`, `wood_apoplast_frac` | apoplastic fraction [–] |
| $w_{sat}$ | `leaf_water_sat`, `wood_water_sat` | saturated water content [kg H₂O / kgC] |
| $\psi_{50}$ | `wood_psi50` | xylem potential at 50 % loss of conductance [MPa] |
| $a$ | `wood_kexp` | vulnerability-curve shape [–] |
| $`k_{plant\_max}`$ | `k_plant_max` | max whole-plant conductance [kg s⁻¹ MPa⁻¹ m⁻²_leaf] |
| — | `conductance` | `"whole_plant"` (default) or `"segment"`: which of the two rows below sets $k_{\text{cond}}$ |
| $K_s$ | `wood_kmax` | sapwood specific conductivity (segment mode) [kg m⁻¹ s⁻¹ MPa⁻¹] |
| — | `vessel_curl` | tortuosity / path-length factor (segment mode) [–] |
| — | `rhizo_cond` | rhizosphere conductance (single-BC) [kg s⁻¹ MPa⁻¹] |
| $\beta$ | `root_beta` | the root profile within the rooting depth [–], $0<\beta<1$ (ED2: 0.1) |
| $`b_{1Rd}, b_{2Rd}`$ | `root_depth_b1`, `root_depth_b2` | rooting depth $`b_{1Rd}\,h^{b_{2Rd}}`$ [m], capped at the soil column (ED2 IALLOM 1) |
| SRL | `specific_root_length` | fine-root length per unit fine-root carbon [m kgC⁻¹] |
| $`r_{\text{root}}`$ | `fine_root_radius` | absorbing-root radius [m] |
| $b_{1SA}, b_{2SA}$ | `sapwood_area_b1`, `sapwood_area_b2` | ED2 sapwood-area allometry; sets the sapwood ring, which is BOTH the hydraulic capacitance and the wood thermal store's internal water |

## References
- Xu, Medvigy, Powers, Becknell & Guan (2016), *New Phytologist* — X16 hydraulics.
- Bartlett, Scoffoni & Sack (2012); Tyree & Hammel (1972) — pressure–volume theory.
- Gardner (1960) *Soil Science* 89:63–73; Cowan (1965) *J. Appl. Ecol.* 2:221–239 — single-root rhizosphere conductance.
- Christoffersen (2013), PhD thesis, University of Arizona — the rooting-depth allometry (ED2 IALLOM 1).
- ED2 `ED/src/dynamics/plant_hydro.f90`; `docs/dev_plans/archive/MEDS_HYDRAULICS_DESIGN.md` (§4 governing
  equations, §16 per-layer roots); `docs/dev_plans/archive/MEDS_MULTILAYER_ROOTS_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| PV curves / capacitance | `meds_water_retention`: one `water_curve_t` record per tissue; `water_content`, `capacitance`, `psi_from_water_content`, `pv_psi_tlp`, `pv_rwc_tlp`, `psi_from_rwc`, `rwc_from_psi` |
| vulnerability + Kirchhoff | `meds_hydr_lib`: `plc_retained`, `flux_potential`, `kirchhoff_edge` (+ table: `build_hydro_table`, `flux_potential_lin`, `kirchhoff_edge_tab`) |
| Kirchhoff quadrature | `meds_hydr_lib`: `kirchhoff_integral` (7-point Gauss–Legendre, written out) |
| multi-layer root boundary | `meds_plant_hydraulics`: `cohort_root_depth`, `cohort_root_profile`, `rhizosphere_cond`, `effective_root_boundary` |
| network solver | `meds_plant_hydraulics`: `solve_plant_water` (`freeze_coeffs` + `advance_exact_linear` + `exact_substep`) |
| config flatten / soil coupling | `meds_fast_types`: `apply_hydraulics_config`; opt-in per-layer soil↔plant in `column_fast_step` (`soil_hydr_cond_from_theta` → K(θ)) |
