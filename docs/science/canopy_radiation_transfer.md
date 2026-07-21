# Canopy radiative transfer

How MEDS partitions incident shortwave and longwave radiation among the cohorts of a patch and the
ground. It is a faithful reimplementation of ED2's **two-stream** multi-layer canopy solver (the
`icanrad=2` scheme in `../ED2/ED/src/dynamics/twostream_rad.f90`; Liou 2002 ch.6, Longo et al. 2019),
generalized to run **one solver over every band** — the default set is VIS, NIR, and thermal LW. The
returned per-cohort **absorbed leaf VIS** is the field the leaf-photosynthesis kernel consumes as PAR,
and the absorbed leaf+wood SW/LW drive the leaf/wood energy balance. The library is state-free and
device-eligible: it takes plain arrays and value types, never a `site_t`.

Symbols: $\mu = \cos\theta_{sun}$ (solar-zenith cosine, floored $>0$); $\rho, \tau$ = leaf/wood
reflectance and transmittance per band; $\omega = \rho+\tau$ = single-scatter albedo; $g$ = leaf-angle
scattering asymmetry; the diffuse **backscatter** $\beta$ and the direct-beam **upscatter** $\beta_0$;
$k$ = direct-beam extinction; and the effective (clumping-corrected) area indices ELAI, EWAI. All
fluxes are **absolute W m⁻²** — there is no ED2-style normalize-to-unit-incidence.

The pipeline has four stages, top to bottom: (1) the leaf-angle distribution and its integrals
(`meds_optics_lib`, leaf-angle block); (2) the per-PFT, $\mu$-independent scattering table
(`meds_optics_lib`, canopy block); (3) the per-cohort optics blend + solar geometry
(`blend_cohort_optics`); and (4) the single-band adding solver — all three assembled, solved, and
sealed inside the public seam `meds_canopy_radiation`. The pure optical-property kernels live in the
shared library `meds_optics_lib` (`src/shared/functions/`); the RT assembly, the two-stream solver,
and the seam live together in `meds_canopy_radiation`.

## 1. Leaf-angle distribution and the G-function

The leaf-inclination distribution (LIDF) is a two-parameter **Beta** density on the normalized
inclination $`t=\theta/(\pi/2)\in[0,1]`$ (Goel & Strebel 1984; the SCOPE / 4SAIL family, Verhoef 1984,
van der Tol et al. 2009). A single generic family subsumes the Verhoef archetypes (planophile,
erectophile, plagiophile, extremophile, and uniform $=\mathrm{Beta}(1,1)$). MEDS uses it in place of
ED2/CLM's $`\phi_1/\phi_2/\bar\mu`$ + $`(1+\chi_L)^2`$ backscatter approximation — a deliberate
departure (see `docs/dev_plans/radiative_transfer_design.md`).

The class weights are integrals of the unnormalized Beta kernel over the SCOPE 13-class inclination
grid, renormalized to sum to one (the normalizing constant cancels):

```math
w_c \;\propto\; \int_{t_{c-1}}^{t_c} t^{\,p-1}(1-t)^{\,q-1}\,dt , \qquad \sum_c w_c = 1 \qquad(1)
```

evaluated by a 64-panel midpoint rule (`beta_lidf` over `beta_pdf_kernel`) — compiler-robust, no
special function, and it never touches the integrable endpoint singularities at $t=0,1$. An
interpretable $`(\overline{\theta},\,\sigma_\theta)`$ maps to $(p,q)$ by moment matching
(`beta_params_from_mean`).

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
$G$ is band-independent, `canopy_radiation` evaluates it **once per cohort** and reuses it for every
beam band.

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
conditions differ. It uses the **adding / layer-recursion** method: an $O(N)$, unconditionally stable
direct solve of the block-tridiagonal two-stream system (the crown fraction is $\mathrm{cai}=1$).

**Per-layer diffuse reflectance/transmittance.** For a homogeneous layer of effective area index
ETAI, `layer_rt` gives the classical Stokes two-stream coefficients with backward/forward diffuse
scattering $`\sigma_b=\omega\beta`$, $`\sigma_f=\omega(1-\beta)`$, attenuation
$`\mathrm{att}=1-\sigma_f`$:

```math
m=\sqrt{\mathrm{att}^2-\sigma_b^2},\quad
r_\infty=\frac{\mathrm{att}-m}{\sigma_b},\quad
e=e^{-m\,\mathrm{ETAI}}
```

```math
t_{dd}=\frac{(1-r_\infty^2)\,e}{1-r_\infty^2 e^2}, \qquad
r_{dd}=\frac{r_\infty(1-e^2)}{1-r_\infty^2 e^2} \qquad(8)
```

(The conservative-scattering singularity $\omega\to1$ is guarded; $`\sigma_b\to0`$ collapses to pure
Beer attenuation, $r_{dd}=0$.)

**Direct beam.** The collimated beam attenuates by Beer's law from the top down, and the fraction
intercepted by layer $i$ is single-scattered into the diffuse field:

```math
I^\downarrow_0(i) = I^\downarrow_0(i{+}1)\,e^{-k_i\,\mathrm{ETAI}_i}, \qquad
S_i = \omega_i\,\big(I^\downarrow_0(i{+}1)-I^\downarrow_0(i)\big) \qquad(9)
```

split into the layer's up-source $`\beta_{0,i}S_i`$ and down-source $`(1-\beta_{0,i})S_i`$. For the
thermal band each layer adds a blackbody source $`(1-r_{dd}-t_{dd})\,\sigma T^4`$ to both faces
(emission temperature = the cohort/canopy temperature; the driver sets it to the canopy-air temperature
so leaf emission is counted once).

**The recursion.** A bottom-up sweep builds the reflectance $`r_{bel}`$ and upward source $`s_{bel}`$
of the growing stack {ground + layers $1..i$}, seeded by the ground reflectance and its upward source
(thermal emission plus the reflected direct beam reaching the ground):

```math
r_{bel}(i) = r_{dd}(i) + \frac{t_{dd}(i)^2\,r_{bel}(i{-}1)}{1-r_{bel}(i{-}1)\,r_{dd}(i)},\qquad
s_{bel}(i) = s_{up}(i) + t_{dd}(i)\,\frac{s_{bel}(i{-}1)+r_{bel}(i{-}1)\,s_{dn}(i)}{1-r_{bel}(i{-}1)\,r_{dd}(i)}
\qquad(10)
```

A top-down sweep (seeded by the incident diffuse at canopy top) then recovers the up/down diffuse
fluxes at every interface. **Absorbed per layer** is the beam plus diffuse flux divergence — energy is
conserved by construction (in − out telescopes):

```math
A_i = \big(I^\downarrow_0(i{+}1)-I^\downarrow_0(i)\big)
      + \big(d^\downarrow_{in}(i)+u^\uparrow_{up}(i)\big)
      - \big(u^\uparrow_{top}(i)+d^\downarrow_{bot}(i)\big) \qquad(11)
```

The solver also returns the below-canopy downwelling $`d_{ground}`$, the upwelling from the ground
$`u_{ground}`$, and the patch **albedo** = upward flux leaving the top / total incident.

## 4. Ground optics and the returned fluxes

`ground_optics` closes the two-stream at the bottom (`meds_canopy_radiation`, surface block). For a shortwave
band the ground reflectance is the (per-band) soil albedo and its emission is zero; for the thermal
band the reflectance is $1-\varepsilon_{soil}$ and the emission is $`\varepsilon_{soil}\,\sigma\,
T_{soil}^4`$. This first cut models **bare soil only**; `surface_state_t` reserves the fields
(soil moisture, litter, standing water, snow) a full implementation will consult — and the snow model
already supplies a Niu-Yang-blended albedo through the driver.

`canopy_radiation` (the sealed public seam, `meds_canopy_radiation`) loops the configured bands,
blends optics per cohort, and dispatches each band to `solve_band`. It then splits each cohort's
absorbed radiation into leaves and wood by the same `leaf_frac` weight the optics blend produced:

```math
A^{leaf}_{b,i} = A_{b,i}\,\mathrm{leaf\_frac}_i, \qquad
A^{wood}_{b,i} = A_{b,i}\,(1-\mathrm{leaf\_frac}_i) \qquad(12)
```

A patch with **no resolvable cohorts** short-circuits: all radiation reaches the ground, and the
albedo is the ground's alone.

> **Incident-equivalent PAR (load-bearing).** The absorbed leaf VIS $`A^{leaf}_{\mathrm{VIS},i}`$ is
> *absorbed* flux, but the leaf gas-exchange kernel expects *incident* PAR (it re-applies leaf
> absorptance internally when it computes $`I_2`$ for electron transport). The driver therefore passes
> $`A^{leaf}_{\mathrm{VIS}} / \alpha_{leaf}`$ as the incident-equivalent PAR — dividing by leaf
> absorptance. Skipping this double-counts absorptance and runs photosynthesis ~15 % light-starved
> (`apply_rt_forcing` in `meds_fast_loop`; see also `project_meds_rt_join`). The absorbed leaf+wood
> SW/LW go to the energy balance *as absorbed* (no division).

## Parameters (config names / trait table)

| Symbol | Source | Meaning |
|---|---|---|
| $p, q$ | `beta_p`, `beta_q` (per PFT) | Beta-LIDF shape parameters |
| $`\overline{\theta}, \sigma_\theta`$ | (via `beta_params_from_mean`) | interpretable mean/std leaf inclination [deg] |
| $`\rho_{leaf},\tau_{leaf}`$ | per-PFT leaf spectra | leaf reflectance / transmittance per band |
| $`\rho_{wood},\tau_{wood}`$ | per-PFT wood spectra | wood reflectance / transmittance per band |
| $`\Omega_{leaf},\Omega_{wood}`$ | `clumping_leaf`, `clumping_wood` | clumping factors $\in(0,1]$ |
| $`\varepsilon_{soil}`$ | `soil_emiss` | ground thermal emissivity |
| — | `soil_albedo(band)` | per-band bare-soil shortwave albedo |
| $`\alpha_{leaf}`$ | `leaf_absorptance` (leaf kernel) | leaf VIS absorptance (PAR re-normalization) |

## References
- Longo et al. (2019), *GMD* 12:4309 — ED-2.2 technical description (the `icanrad=2` two-stream).
- Liou (2002), *An Introduction to Atmospheric Radiation*, ch. 6 — two-stream theory.
- Verhoef (1984), *Remote Sens. Environ.* 16:125–141 — SAIL / the diffuse two-stream layer operator.
- van der Tol et al. (2009), *Biogeosciences* 6:3109 — SCOPE (leaf-angle + volscatt).
- Goel & Strebel (1984), *Agron. J.* 76:800 — Beta leaf-inclination distribution.
- ED2 `../ED2/ED/src/dynamics/twostream_rad.f90`; `docs/dev_plans/radiative_transfer_design.md`,
  `docs/dev_plans/MEDS_ENERGY_BALANCE_DESIGN.md`.

## Code map

| Concept | Routine |
|---|---|
| Beta LIDF + params | `meds_optics_lib`: `beta_lidf`, `beta_pdf_kernel`, `beta_params_from_mean` |
| moments $`\mathrm{bf}`$, $G(\mu)$ | `meds_optics_lib`: `leaf_bf`, `gfun_direct` |
| $\omega, g$ single-scatter pair (per PFT) | `meds_optics_lib`: `scatter_pair` |
| $\omega, g$ optics table + per-cohort blend | `meds_canopy_radiation`: `derive_rad_optics`, `blend_cohort_optics` ($\beta,\beta_0,k,$ leaf_frac) |
| ground reflectance / emission | `meds_canopy_radiation`: `ground_optics` (`surface_state_t`) |
| single-band adding solve | `meds_canopy_radiation`: `solve_band`, `layer_rt` |
| public per-patch seam | `meds_canopy_radiation`: `canopy_radiation` |
| types (optics / forcing / flux / surface) | `meds_biophysics_types`: `rad_pft_optics_t`, `rad_forcing_t`, `rad_flux_t`, `surface_state_t` |
| fast-loop join (PAR renorm) | `meds_fast_loop`: `apply_rt_forcing` |
