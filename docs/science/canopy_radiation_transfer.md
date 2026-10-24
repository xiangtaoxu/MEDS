# Canopy radiative transfer

How MEDS partitions incident shortwave and longwave radiation among the cohorts of a patch and the
ground. It is built on ED2's **two-stream** multi-layer canopy optics (the `icanrad=2` scheme in
`../ED2/ED/src/dynamics/twostream_rad.f90`; Liou 2002 ch.6, Longo et al. 2019), solved by the adding
method of §3 and generalized to run **one solver over every band** — the default set is VIS, NIR, and thermal LW. The
returned per-cohort **absorbed leaf VIS** is the field the leaf-photosynthesis kernel consumes as PAR,
and the absorbed leaf+wood SW/LW drive the leaf/wood energy balance. The library is state-free and
device-eligible: it takes plain arrays and value types, never a `polygon_t`.

Symbols: $\mu = \cos\theta_{sun}$ (solar-zenith cosine; the extinction $k$ guards $\mu\to0$); $\rho, \tau$ = leaf/wood
reflectance and transmittance per band; $\omega = \rho+\tau$ = single-scatter albedo; $g$ = leaf-angle
scattering asymmetry; the diffuse **backscatter** $\beta$ and the direct-beam **upscatter** $\beta_0$;
$k$ = direct-beam extinction; and the effective (clumping-corrected) area indices ELAI, EWAI. All
fluxes are **absolute W m⁻²** — there is no ED2-style normalize-to-unit-incidence.

The pipeline has four stages, top to bottom: (1) the leaf-angle distribution and its integrals
(`meds_optics_lib`); (2) the per-PFT, $\mu$-independent scattering table, built once at load
(`derive_rad_optics`, `meds_pft_derive`); (3) the per-cohort optics blend + solar geometry
(`blend_cohort_optics`); and (4) the single-band adding solver. Stages 3 and 4 run inside
`canopy_radiation` (`meds_canopy_radiation`). The pure optical-property kernels live in the
shared library `meds_optics_lib` (`src/share/functions/`); the RT assembly, the two-stream solver,
and its per-patch entry point live together in `meds_canopy_radiation`.

## 1. Leaf-angle distribution and the G-function

The leaf-inclination distribution (LIDF) is a two-parameter **Beta** density on the normalized
inclination $`t=\theta/(\pi/2)\in[0,1]`$ (Goel & Strebel 1984; the SCOPE / 4SAIL family, Verhoef 1984,
van der Tol et al. 2009). A single generic family subsumes the Verhoef archetypes (planophile,
erectophile, plagiophile, extremophile, and uniform $=\mathrm{Beta}(1,1)$). MEDS uses it in place of
ED2/CLM's $`\phi_1/\phi_2/\bar\mu`$ + $`(1+\chi_L)^2`$ backscatter approximation — a deliberate
departure (see `docs/dev_plans/archive/radiative_transfer_design.md`).

The class weights are integrals of the unnormalized Beta kernel over the SCOPE 13-class inclination
grid, renormalized to sum to one (the normalizing constant cancels):

```math
w_c \;\propto\; \int_{t_{c-1}}^{t_c} t^{\,p-1}(1-t)^{\,q-1}\,dt , \qquad \sum_c w_c = 1 \qquad(1)
```

evaluated by a 64-panel midpoint rule (`beta_lidf` over `beta_log_kernel`) — compiler-robust, no
special function, and it never touches the integrable endpoint singularities at $t=0,1$. The kernel
is summed relative to its largest value on the grid, so a narrow distribution, whose unnormalized
density underflows (a spread of a few degrees, or 0), still sums to one. An
interpretable $`(\overline{\theta},\,\sigma_\theta)`$ maps to $(p,q)$ by moment matching
(`moments_to_beta_params`).

Two quantities feed the two-stream optics. The **second moment** (mutually $\mu$-independent, driving
the back/forward scattering split) is

```math
\mathrm{bf} \;=\; \langle\cos^2\theta_{leaf}\rangle \;=\; \sum_c w_c\,\cos^2\theta_c \qquad(2)
```

($`\mathrm{bf}=1`$ horizontal leaves, $0$ vertical, $1/3$ spherical; `leaf_bf`). The **Ross
G-function** — the mean leaf-area projection toward the sun, $\mu$-dependent — uses the exact SCOPE
`volscatt` projection $`\chi_s(\theta_c,\mu)`$:

```math
G(\mu) \;=\; \sum_c w_c\,\chi_s(\theta_c,\mu), \qquad
k \;=\; \frac{G(\mu)}{\mu} \qquad(3)
```

`gfun_direct` returns $G(\mu)$; the direct-beam extinction coefficient is $`k=G(\mu)/\mu`$. Because
$G$ is band-independent, `canopy_radiation` evaluates it **once per PFT present** and reuses it for
every cohort of that PFT and every beam band.

## 2. Single-scatter albedo, asymmetry, and the scattering coefficients

For each PFT and band the $\mu$-independent optics are built once (`derive_rad_optics` over
`scatter_pair`):

```math
\omega = \rho + \tau, \qquad
g = \mathrm{bf}\,\frac{\rho-\tau}{\omega} \qquad(4)
```

($g=0$ for a perfect absorber, $\omega\to0$). For the **thermal band** the caller passes
$`\rho = 1-\varepsilon`$, $\tau=0$ (ED2/SCOPE convention), so $`\omega = 1-\varepsilon`$ and
$`g=\mathrm{bf}`$ — no special case is needed.

Per cohort, per band, `blend_cohort_optics` area-weights the leaf and wood optics by their
clumping-corrected areas $`\mathrm{ELAI}=\Omega_{leaf}\,\mathrm{LAI}`$,
$`\mathrm{EWAI}=\Omega_{wood}\,\mathrm{WAI}`$ (weights $`w_l=\mathrm{ELAI}/\mathrm{ETAI}`$,
$`w_w=\mathrm{EWAI}/\mathrm{ETAI}`$, $`\mathrm{ETAI}=\mathrm{ELAI}+\mathrm{EWAI}`$), then adds the
solar geometry:

```math
\omega = w_l\,\omega_{leaf} + w_w\,\omega_{wood}, \qquad
g = w_l\,g_{leaf} + w_w\,g_{wood} \qquad(5)
```

```math
\beta = \tfrac12\,(1 + g), \qquad
\beta_0 = \tfrac12\!\left(1 + \frac{g}{k}\right) \qquad(6)
```

$\beta$ is the **diffuse backscatter** fraction and $\beta_0$ the **direct-beam upscatter** fraction
(only meaningful, and only computed, for beam bands). The routine also returns the **leaf absorption
share** — the fraction of a cohort's absorbed radiation that lands on leaves rather than wood,
weighted by absorptivity $1-\omega$ and clumped area:

```math
\mathrm{leaf\_frac} = \frac{(1-\omega_{leaf})\,\mathrm{ELAI}}
{(1-\omega_{leaf})\,\mathrm{ELAI} + (1-\omega_{wood})\,\mathrm{EWAI}} \qquad(7)
```

## 3. The single-band two-stream (adding) solve

`solve_band` (`meds_canopy_radiation`) is **the** unified single-band solver over the vertical stack of
cohort layers, ordered **BOTTOM (index 1) to TOP (index n)** with a virtual sky interface above.
Shortwave bands carry a direct beam and no emission; the thermal band carries emission and no beam;
the diffuse multiple-scattering operator is identical across bands — only the source terms and boundary
conditions differ. It uses the **adding / layer-recursion** method: two $O(N)$ sweeps, bottom-up for
the reflectance and upward source of the stack below each layer, then top-down for the fluxes; each
layer covers the ground (crown area index 1). A layer's diffuse operator, the beam it scatters and its
emission are each the exact solution for a homogeneous layer, so the adding method is exact for the
stack, and a uniform canopy's light does not depend on how it is cut into cohorts.

**Per-layer diffuse reflectance/transmittance.** For a homogeneous layer of effective area index
ETAI, `layer_rt` gives the classical Stokes two-stream coefficients with backward/forward diffuse
scattering $`\sigma_b=\omega\beta`$, $`\sigma_f=\omega(1-\beta)`$, attenuation
$`\mathrm{att}=1-\sigma_f`$:

```math
m=\sqrt{\mathrm{att}^2-\sigma_b^2},\quad
r_\infty=\frac{\mathrm{att}-m}{\sigma_b},\quad
e=e^{-m\,\mathrm{ETAI}} \qquad(8)
```

```math
t_{dd}=\frac{(1-r_\infty^2)\,e}{1-r_\infty^2 e^2}, \qquad
r_{dd}=\frac{r_\infty(1-e^2)}{1-r_\infty^2 e^2} \qquad(9)
```

(The conservative-scattering singularity $\omega\to1$ is guarded; $`\sigma_b\to0`$ collapses to pure
Beer attenuation, $r_{dd}=0$.)

**Direct beam.** The collimated beam attenuates by Beer's law from the top down,

```math
I^\downarrow_0(i) = I^\downarrow_0(i{+}1)\,e^{-k_i\,\mathrm{ETAI}_i} \qquad(10)
```

and inside a layer, at depth $x$ in area index, it scatters $`\omega k S_0 e^{-kx}`$ into the diffuse
field, $\beta_0$ of it upward ($`S_0 = I^\downarrow_0(i{+}1)`$, the beam at the layer's top). The
layer's up- and down-sources are what of this leaves its top and its bottom, from the two-stream
equations solved through its depth with no diffuse entering either face (`layer_beam_sources`): with
the attenuation att and backscatter $`\sigma_b`$ of eq. (8), the particular solution is

```math
u_p = a\,e^{-kx},\quad d_p = b\,e^{-kx},\qquad
a = \frac{\omega k S_0\,[\beta_0(\mathrm{att}-k) + \sigma_b(1-\beta_0)]}{m^2-k^2},\quad
b = \frac{\omega k S_0\,[(1-\beta_0)(\mathrm{att}+k) + \sigma_b\beta_0]}{m^2-k^2} \qquad(10b)
```

and the diffuse modes $`e^{-mx}`$ and $`e^{m(x-\mathrm{ETAI})}`$, in the ratio $`r_\infty`$ of up to
down, take the boundary conditions. $k = m$ is a removable singularity, moved off by a relative
$`10^{-6}`$. A thin layer's sources tend to the single-scatter split, $`\beta_0\,\omega\,\Delta I`$
up and $`(1-\beta_0)\,\omega\,\Delta I`$ down of what it intercepts; a thick one's do not, because
light scattered deep in the layer is absorbed on its way out. Injecting the single-scatter split at the
faces would make the light depend on the cut: with the default optics, LAI 5, μ 0.87 and 70 % beam, the
VIS albedo would be 0.075 for one cohort, 0.053 for five and 0.043 in the limit, which the exact source
gives for any cut. For the thermal band each layer adds a blackbody source
$`(1-r_{dd}-t_{dd})\,\sigma T^4`$ to both faces. Each tissue emits at its absorptivity (Kirchhoff), so
a cohort's $`\sigma T^4`$ is its tissues' at their start-of-step temperatures, weighted as their
absorption is, by the thermal band's `leaf_frac` (eq. (7)):
$`\sigma T^4 = \mathrm{leaf\_frac}\,\sigma T_{leaf}^4 + (1-\mathrm{leaf\_frac})\,\sigma T_{wood}^4`$. The
longwave field thus sees both tissues.

**The recursion.** A bottom-up sweep builds the reflectance $`r_{bel}`$ and upward source $`s_{bel}`$
of the growing stack {ground + layers $1..i$}, seeded by the ground reflectance and its upward source
(thermal emission plus the reflected direct beam reaching the ground):

```math
r_{bel}(i) = r_{dd}(i) + \frac{t_{dd}(i)^2\,r_{bel}(i{-}1)}{1-r_{bel}(i{-}1)\,r_{dd}(i)},\qquad
s_{bel}(i) = s_{up}(i) + t_{dd}(i)\,\frac{s_{bel}(i{-}1)+r_{bel}(i{-}1)\,s_{dn}(i)}{1-r_{bel}(i{-}1)\,r_{dd}(i)}
\qquad(11)
```

A top-down sweep (seeded by the incident diffuse at canopy top) then recovers the up/down diffuse
fluxes at every interface. **Absorbed per layer** is the beam plus diffuse flux divergence — energy is
conserved by construction (in − out telescopes):

```math
A_i = \big(I^\downarrow_0(i{+}1)-I^\downarrow_0(i)\big)
      + \big(d^\downarrow_{in}(i)+u^\uparrow_{up}(i)\big)
      - \big(u^\uparrow_{top}(i)+d^\downarrow_{bot}(i)\big) \qquad(12)
```

The solver also returns the below-canopy downwelling $`d_{ground}`$, the upwelling from the ground
$`u_{ground}`$, and the patch **albedo** = upward flux leaving the top / total incident.

## 4. Ground optics and the returned fluxes

`ground_optics` closes the two-stream at the bottom (`meds_canopy_radiation`, surface block). For a shortwave
band the ground reflectance is the (per-band) soil albedo and its emission is zero; for the thermal
band the reflectance is $1-\varepsilon_{soil}$ and the emission is $`\varepsilon_{soil}\,\sigma\,
T_{soil}^4`$. Its input (`ground_optics_env_t`) is **bare soil plus snow**: the per-band soil albedo,
the emissivity and the skin temperature, which the driver blends toward the snowpack's albedo and
emissivity by the Niu-Yang snow-cover fraction, so a partial pack changes the ground continuously; the
snow's own albedo moves from its fresh to its aged value with the liquid fraction of the pack's surface
layer. The blend's temperature is the one at which the blended emissivity emits what the soil and the
pack emit at their own temperatures. What reaches the ground in each band is then shared between the
two at their own optics ([snow_biophysics](snow_biophysics.md), "The ground's radiation"). Litter,
standing water and soil moisture do not change the ground optics (see docs/ROADMAP.md).

`canopy_radiation` (the per-patch entry point, `meds_canopy_radiation`) loops the configured bands,
blends optics per cohort, and dispatches each band to `solve_band`.

**Co-dominant cohorts share a layer.** Cohorts of exactly equal height form one RT layer: its ETAI is
the members' sum, its $\omega$, $g$ (hence $\beta$, $\beta_0$) and $k$ are ETAI-weighted means of
theirs, and its $`\sigma T^4`$ is weighted, as within a cohort, by absorptivity-weighted area
$`(1-\omega_i)\,\mathrm{ETAI}_i`$. The solve's absorbed flux $`A_\ell`$ is net of the layer's emission.
Every element of the layer meets the same radiation, so what the layer takes in, $`A_\ell`$ plus the
emission $`E_\ell=\sum_{j\in\ell}E_j`$ of its members, goes back to them by absorptivity-weighted
area, the weight `leaf_frac` is built from, and each is charged its own emission:

```math
A_{b,i} = \big(A_{b,\ell}+E_{b,\ell}\big)\,\frac{(1-\omega_i)\,\mathrm{ETAI}_i}{\sum_{j\in\ell}(1-\omega_j)\,\mathrm{ETAI}_j}
          - E_{b,i}
\qquad(13)
```

A cohort's emission is its tissues', each from both faces of its area at its absorptivity,
$`E^{leaf}_{b,i} = 2\,(1-\omega_{leaf})\,\mathrm{ELAI}_i\,\sigma T_{leaf}^4`$ and
$`E^{wood}_{b,i} = 2\,(1-\omega_{wood})\,\mathrm{EWAI}_i\,\sigma T_{wood}^4`$, zero in a shortwave
band. Every member sees the light reaching the layer's top. A one-member layer is exactly the per-cohort
solve. Stacking tied cohorts instead would put whichever sorts first on top, so the light split between
PFTs would follow the PFT index rather than the ecology. The canopy total barely moves when they are
stacked, so no whole-canopy albedo or GPP check would show the error; only the split between the tied
cohorts changes. Ties arise at recruitment: every PFT recruits at the one smallest tracked height,
`min_cohort_height`.

`canopy_radiation` then splits each cohort's absorbed radiation into leaves and wood the same way:
they share what the cohort takes in by the `leaf_frac` weight the optics blend produced, and each is
charged its own emission:

```math
A^{leaf}_{b,i} = \big(A_{b,i}+E_{b,i}\big)\,\mathrm{leaf\_frac}_i - E^{leaf}_{b,i}, \qquad
A^{wood}_{b,i} = \big(A_{b,i}+E_{b,i}\big)\,(1-\mathrm{leaf\_frac}_i) - E^{wood}_{b,i} \qquad(14)
```

In a shortwave band this is the split by `leaf_frac`. In the thermal band it is exact for the
two-stream's well-mixed layer: leaves and wood meet the same radiation through its depth and absorb
it in the proportion `leaf_frac`, and the solve's net is what the layer absorbs less what it emits.
Wood warmer than its leaves warms them. With the sky, the ground and the leaves at 295 K and the wood
at 305 K, a cohort of LAI 3 and WAI 0.6 gives its leaves +19.4 W/m² and its wood −30.8 W/m², just
what leaves at 305 K give wood at 295 K. Split by `leaf_frac` alone, the cohort's net of −11.4 W/m²
would charge the cool leaves −10 W/m², a share of the warm wood's emission.

A patch with **no resolvable cohorts** short-circuits: all radiation reaches the ground, and the
albedo is the ground's alone.

> **Incident-equivalent PAR (load-bearing).** The absorbed leaf VIS $`A^{leaf}_{\mathrm{VIS},i}`$ is
> *absorbed* flux, but the leaf gas-exchange kernel expects *incident* PAR (it re-applies leaf
> absorptance internally when it computes $`I_2`$ for electron transport). The driver therefore passes
> $`A^{leaf}_{\mathrm{VIS}} / \alpha_{leaf}`$ as the incident-equivalent PAR — dividing by leaf
> absorptance. Skipping this double-counts absorptance and runs photosynthesis ~15 % light-starved
> (`apply_rt_forcing` in `meds_fast_forcing`). The absorbed leaf+wood SW/LW go to the energy balance
> *as absorbed* (no division).

## Parameters (config names / trait table)

| Symbol | Source | Meaning |
|---|---|---|
| $`\overline{\theta}, \sigma_\theta`$ | `[pft].leaf_angle_mean`, `leaf_angle_std` | mean / std leaf inclination [deg] |
| $p, q$ | derived from $`\overline{\theta}, \sigma_\theta`$ (`moments_to_beta_params`) | Beta-LIDF shape parameters |
| $`\rho_{leaf},\tau_{leaf}`$ | `[pft].leaf_reflect_*`, `leaf_transmit_*` (`_vis`, `_nir`) | leaf reflectance / transmittance |
| $`\rho_{wood},\tau_{wood}`$ | `[pft].wood_reflect_*`, `wood_transmit_*` (`_vis`, `_nir`) | wood reflectance / transmittance |
| $`\varepsilon_{leaf},\varepsilon_{wood}`$ | `[pft].leaf_emissivity`, `wood_emissivity` | LW: $`\rho=1-\varepsilon,\ \tau=0`$ |
| $`\Omega_{leaf},\Omega_{wood}`$ | `[pft].leaf_clumping`, `wood_clumping` | clumping factors $\in(0,1]$ |
| $`\varepsilon_{soil}`$ | `[soil].ground_emissivity` (0.95) | ground thermal emissivity |
| — | `[soil].ground_albedo_vis`, `ground_albedo_nir` (0.15, 0.30) | per-band bare-soil shortwave albedo |
| $`\alpha_{leaf}`$ | `[leaf_physiology].leaf_absorptance` | leaf VIS absorptance (PAR re-normalization) |

## References
- Longo et al. (2019), *GMD* 12:4309 — ED-2.2 technical description (the `icanrad=2` two-stream).
- Liou (2002), *An Introduction to Atmospheric Radiation*, ch. 6 — two-stream theory.
- Verhoef (1984), *Remote Sens. Environ.* 16:125–141 — SAIL / the diffuse two-stream layer operator.
- van der Tol et al. (2009), *Biogeosciences* 6:3109 — SCOPE (leaf-angle + volscatt).
- Goel & Strebel (1984), *Agron. J.* 76:800 — Beta leaf-inclination distribution.
- ED2 `../ED2/ED/src/dynamics/twostream_rad.f90`; `docs/dev_plans/archive/radiative_transfer_design.md`,
  `docs/dev_plans/archive/MEDS_ENERGY_BALANCE_DESIGN.md`.

## Surface radiative fluxes on the diagnostic path

The two-stream forms the top-of-canopy upward flux every step as `albedo(b) x incid(b)`, and the output
carries it, so a run compares directly against a radiometer or a satellite product. Five variables carry
it:

| variable | meaning |
|---|---|
| `mean_sw_in_vis_py`, `mean_sw_in_nir_py` | incident shortwave at canopy top, per band |
| `mean_sw_up_vis_py`, `mean_sw_up_nir_py` | upwelling shortwave leaving the canopy top |
| `mean_lw_up_py` | upwelling longwave, **surface emission included**, with the change in emission within the step that the tissues' and the ground's balances booked, so the reported budget closes on the booked net radiation ([vegetation_energy_dynamics](vegetation_energy_dynamics.md), "What the output reports") |

**These are fluxes, not a time-averaged albedo, deliberately.** A period-mean albedo is the mean of
a *ratio*, which is not the ratio of the means — and at night the shortwave ratio is 0/0. The albedo
a reader wants, and what a satellite product *is*, is

```math
\alpha_b = \frac{\sum_t \mathrm{SW}^\uparrow_b}{\sum_t \mathrm{SW}^\downarrow_b} \qquad(15)
```

over whatever period they choose. A mean albedo would be a number that looks right and is wrong.

Measured on a spun-up Ithaca stand, annual sums: **VIS 0.083, NIR 0.320, broadband 0.208**. The
VIS/NIR contrast is the vegetation red-edge — chlorophyll absorbs the visible and leaves scatter the
near-infrared — and it is the signal every vegetation index is built on. Winter months rise to 0.23
(VIS) and 0.45 (NIR) as the snow albedo ramps in, and `mean_lw_up_py` runs 305–468 W m⁻², consistent
with $`\sigma T^4`$ across the annual surface-temperature range.

## Code map

| Concept | Routine |
|---|---|
| Beta LIDF + params | `meds_optics_lib`: `beta_lidf`, `beta_log_kernel`; `meds_pft_derive`: `moments_to_beta_params` |
| moments $`\mathrm{bf}`$, $G(\mu)$ | `meds_optics_lib`: `leaf_bf`, `gfun_direct` |
| $\omega, g$ single-scatter pair (per PFT) | `meds_optics_lib`: `scatter_pair` |
| per-PFT optics table ($\omega, g$), eq. (4) | `meds_pft_derive`: `derive_rad_optics`, once at load onto the PFT table (`cfg%pft%optics`) |
| per-cohort blend ($\beta, \beta_0, k$, leaf_frac), eqs. (5)–(7) | `meds_canopy_radiation`: `blend_cohort_optics` |
| ground reflectance / emission | `meds_canopy_radiation`: `ground_optics` (`ground_optics_env_t`) |
| single-band adding solve | `meds_canopy_radiation`: `solve_band`, `layer_rt` |
| the per-patch entry point (+ tied-cohort layers) | `meds_canopy_radiation`: `canopy_radiation` |
| cohorts of one height are one layer (with the wind and the overtopping LAI) | `meds_numerics`: `equal_height_layers` |
| types | `meds_pft_params`: `rad_pft_optics_t`, the bands `RAD_*`; `meds_canopy_radiation`: `rad_env_t`, `rad_out_t`, `ground_optics_env_t` |
| fast-loop join (PAR renorm, the cohorts bottom to top) | `meds_fast_forcing`: `apply_rt_forcing` |
| each tissue's emission and the leaf/wood split, eqs. (13)–(14) | `meds_canopy_radiation`: `canopy_radiation` |
| the ground's optics under snow, and each tile's radiation | `meds_surface`: `snow_blend_ground_optics`, `ground_tile_radiation` |
