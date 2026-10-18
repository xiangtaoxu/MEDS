# Canopy aerodynamics

How MEDS computes the turbulent transfer that couples the canopy air space (CAS) to the free
atmosphere above and to the leaf, wood, and ground surfaces within. The kernel is a **pure function** of
(free-atmosphere forcing, CAS state, canopy geometry) that produces every turbulence/conductance
quantity the other fast-loop kernels take as forced input: the friction velocity $`u_*`$, the scalar
profile factors that set the atm↔CAS conductances of all three CAS twins, the per-cohort in-canopy wind
and leaf/wood boundary-layer conductances, and the ground↔CAS conductance. Like ED2's
`canopy_turbulence8`, it carries **no integrated state** — it is recomputed each fast sub-step from the
current CAS + vegetation state. The physics reconciles ED2 and CLM (see
`docs/dev_plans/archive/MEDS_COLUMN_DYNAMICS_DESIGN.md`): CLM5's clean four-range Monin-Obukhov surface layer,
ED2's forced+free-convection Nusselt boundary layers, ED2's per-cohort wind extinction, and a CLM4-like
ground conductance.

Symbols: $`u_*`$ = friction velocity; $\zeta = z/L$ = Monin-Obukhov stability parameter ($L$ = Obukhov
length); $`z_{0m}, z_{0h}, z_{0q}`$ = roughness lengths for momentum / heat / vapour; $d$ = displacement
height; $`\kappa`$ = von Karman constant; $`\Psi_m,\Psi_h`$ = integrated stability functions;
$`\mathrm{Rib}`$ = bulk Richardson number. Everything is per unit ground area. All routines are `pure`
and write per-cohort outputs into caller-owned SoA slices — never an array-valued function result passed
straight into a call, which nvfortran miscompiles (see `CLAUDE.md`).

## 1. Roughness, displacement, and the reference geometry

Roughness and displacement scale with canopy top height $h$, blended toward the snow-surface roughness
by the burial fraction `snowfac`:

```math
z_{0m} = (1-\mathrm{snowfac})\,c_{z0}\,h + \mathrm{snowfac}\,z_{0,snow}, \qquad
d = c_d\,h\,(1-\mathrm{snowfac}) \qquad(1)
```

with $`c_{z0}=`$ `z0m_ratio` (0.13, ED2 `vh2vr`), $`c_d=`$ `d_ratio` (0.63, ED2 `vh2dh`) and
$`z_{0,snow}=`$ `snow_rough`, which also floors $`z_{0m}`$ (`canopy_roughness`). The driver passes an
unburied canopy ($`\mathrm{snowfac}=0`$), so the snow blend is inactive in a run (see docs/ROADMAP.md);
the open-canopy blend of §3 acts. The canopy displacement height
$`z_{ldis}=\max(z_{ref}-d,\,2z_{0m})`$ is the surface-layer thickness the Monin-Obukhov solve integrates
over, and the **prognostic CAS depth** is $`\max(\text{min\_canopy\_depth},\,h+\text{canopy\_freeboard})`$,
the tallest cohort plus a freeboard ([canopy_air_space_biophysics](canopy_air_space_biophysics.md)) — it
sets the CAS air mass per ground area, hence every twin's storage capacity.

## 2. The CLM5 Monin-Obukhov surface layer

`mo_surface_layer` returns $`u_*`$ and the scalar profile factor from the atmosphere-to-CAS gradients.
It is seeded by the analytic **MoninObukIni** bulk-Richardson guess and then does a **fixed** number of
sweeps (`n_iter_mo`, default 4) — no data-dependent exit, so it is warp-uniform on the GPU and needs no
root-find or I/O. Using virtual potential temperatures
$`\theta_v=\theta\,(1+0.61\,q)`$ and the CAS↔atmosphere differences $`\Delta\theta_v`$:

```math
\mathrm{Rib} = \frac{g\,z_{ldis}\,\Delta\theta_v}{\theta_{v,atm}\,u_m^2}, \qquad
u_m = \begin{cases}\bar u & \Delta\theta_v\ge0\ \text{(stable)}\\[2pt]
\sqrt{\bar u^2+w_c^2} & \Delta\theta_v<0\ \text{(convective)}\end{cases}, \qquad
\bar u = \max(u_{ref},\,u_{b,min}) \qquad(2)
```

($w_c$ = convective velocity scale `wc`; $`u_{b,min}`$ = `ubmin`, the reference-wind floor). Each sweep
evaluates the CLM5 momentum and scalar **profile denominators** $D_m(\zeta)$ and $D_h(\zeta)$ over four
stability ranges (very-unstable / unstable / stable / very-stable, `d_mom` and `d_heat`), then updates
$`u_*`$, the scalar factor, and $\zeta$:

```math
u_* = \max\!\Big(u_{*,min},\,\frac{\kappa\,u_m}{D_m(\zeta)}\Big), \qquad
\mathrm{temp1} = \frac{\kappa}{D_h(\zeta)}, \qquad
\zeta = \frac{z_{ldis}\,\kappa\,g\,\theta_{v*}}{u_*^2\,\theta_{v,atm}} \qquad(3)
```

The unstable-range stability integrals are the standard Businger-Dyer / Paulson forms with
$`x=(1-16\zeta)^{1/4}`$:

```math
\Psi_m(\zeta) = 2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan x + \frac{\pi}{2}, \qquad
\Psi_h(\zeta) = 2\ln\frac{1+x^2}{2} \qquad(4)
```

(`stabfunc1`/`stabfunc2`; the stable branch is $`-5\zeta`$). Careful placement — $`\zeta_{0m}=z_{0m}
\zeta/z_{ldis}=z_{0m}/L`$ appears without dividing by $\zeta$, and $\zeta$ only divides in the very-
unstable/very-stable branches where $|\zeta|>1$ — keeps the whole solve safe under `-fpe0`. The
neutral limit is $`u_*\to\kappa u/\ln(z_{ldis}/z_{0m})`$, $`\mathrm{temp1}\to\kappa/\ln(z_{ldis}/
z_{0m})`$.

### The scalar transfer factors that couple the CAS twins

The canopy uses $`z_{0h}=z_{0m}`$ and $`z_{0q}=z_{0h}`$, so the heat and vapour profile factors are equal:
$`\mathrm{temp2}=\mathrm{temp1}`$. These dimensionless Monin-Obukhov coefficients set the
**atm↔CAS conductances** every CAS twin shares (the single most load-bearing output of this kernel):

```math
g_{ah} = \rho\,u_*\,\mathrm{temp1}, \qquad
g_{aw} = \rho\,u_*\,\mathrm{temp2}, \qquad
g_{ac} = \rho_{mol}\,u_*\,\mathrm{temp2} \qquad(5)
```

$g_{ah}$ [kg m⁻² s⁻¹] drives the enthalpy twin, $g_{aw}$ the vapour twin, and $g_{ac}$ [mol m⁻² s⁻¹]
the CO₂ twin (molar, using dry-air molar density). Dropping the profile factor ($`\mathrm{temp1}=1`$)
over-couples the CAS to the free atmosphere by $`\sim1/\mathrm{temp1}`$ and damps the nocturnal
sub-canopy CO₂ build-up. The kernel also reports the scale stars
$`t_*=\mathrm{temp1}\,\Delta\theta`$, $`q_*=\mathrm{temp2}\,\Delta q`$, $`c_*=\mathrm{temp2}\,\Delta
C`$, and $`\mathrm{Rib}`$, $\zeta$, $L$ as diagnostics.

**The reference level.** $`z_{ref}`$ is the canopy-air top, where the forcing is moved to
([forcing](forcing.md) §8), and both potential temperatures are referenced to it: $`\theta_{atm}`$ is
the air temperature there and $`\theta_{CAS}`$ the canopy-air temperature, so
$`\Delta\theta = T_{atm} - T_{CAS}`$. That is the difference the canopy air's energy budget
exchanges, $`g_{ah}\,(h_{CAS}-h_{atm})`$ on actual temperatures, so the reported sensible heat
$`H = g_{ah}\,c_p\,(T_{CAS}-T_{atm})`$ is the flux the budget books. The two temperatures must share the
reference: referencing $`\theta_{atm}`$ to the ground, $`T_{atm}+(g/c_p)\,z_{ref}`$, while $`\theta_{CAS}`$
is not, would add $`(g/c_p)\,z_{ref}`$ (0.4 K at a 40 m canopy) of spurious stable stratification to every
solve and lower the reported H by $`g_{ah}\,g\,z_{ref}`$, about 100 W m⁻² at Barro Colorado Island.

### These three are re-solved at *every integrator stage*, and that is load-bearing

Almost everything the fast loop computes per-cohort — photosynthesis, stomatal conductance, the leaf
boundary layers, radiation, hydraulics — is **frozen once per `dt_fast`** (see
[numerical_scheme](numerical_scheme.md) §2). The conductances of eq. (5) are the exception:
`mo_surface_layer` is re-solved at each stage of the time integrator, against that stage's own canopy-air
temperature and humidity, and only $`g_{ah}`$, $`g_{aw}`$, $`g_{ac}`$ are updated
(`refresh_cas_conductances`). The per-cohort boundary layers around it stay frozen.

The reason is a feedback loop that runs entirely through this kernel. A warmer canopy air makes the
surface layer **more unstable**, which raises $`u_*`$ *and* $`\mathrm{temp1}`$, which vents the canopy
air harder and cools it. Measured sensitivity on a moderately-windy 18 m stand:
$`\mathrm{d}\ln g_{ah}/\mathrm{d}T \approx 2.2\ \mathrm{K^{-1}}`$ — a canopy air 0.5 K warmer triples
$`g_{ah}`$; 2 K warmer raises it twelvefold. That is a strong negative feedback, and evaluating it one
whole `dt_fast` behind the state it responds to turns it into a **numerical oscillation**: canopy-air
temperature alternating up to ~8 K between consecutive steps, with every conservation budget closing
to ~10⁻⁶ J throughout.

Re-solving per stage removes it. The measured one-step amplification factor of the canopy air falls
from −23.2 to −0.14 at `dt_fast = 900 s`, and the step-to-step swing from 7.7 K to 0.10 K. Holding
each of the other frozen coefficients fixed instead ($`g_{tr}`$, the leaf sensible coefficient, the
longwave emission base, the wetted fraction) changes the amplification by ≤ 0.7%, so this is the one
that matters. It is also cheap: `canopy_aerodynamics` is ~2% of the frozen pre-pass, against ~89% for
leaf gas exchange at 30 cohorts, so the re-solve costs a few percent of a sub-step and pays for itself
in fewer sub-steps.

There is **no switch** for this. ED2 does the same thing (`update_diagnostic_vars` →
`canopy_turbulence8` at every RK stage); the frozen alternative is numerically unstable at production
step sizes on most stand heights, so it is not a supported configuration. The full measurement record
is `docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §1g–1h.

## 3. Ground ↔ CAS conductance

The ground-to-CAS scalar conductance blends a bare-ground similarity term with CLM4's dense-canopy
form (`ggbare`/`ggveg`/`ggnet`), weighted by the patch's open fraction $W$ as CLM weighs them
(Zeng et al. 2005; the stability factor is Sakaguchi and Zeng 2009's):

```math
g_{g,bare} = u_*\,\mathrm{temp1}, \qquad
g_{g,veg} = \frac{c_s\,u_*}{1+\gamma_g\,\mathrm{stab}}, \qquad
g_{g,net} = W\,g_{g,bare} + (1-W)\,g_{g,veg}, \qquad
W = e^{-(\mathrm{LAI}+\mathrm{WAI})} \qquad(6)
```

where `stab` is a bounded buoyancy factor in the CAS↔ground temperature difference, $`c_s=`$
`cs_dense`, $`\gamma_g=`$ `gamma_g`, and LAI and WAI are the patch's leaf and wood area indices
(`open_canopy_fraction`; the kernel takes $W$ as `opencan_frac`). A patch with no cohorts ($W=1$)
couples its ground as bare soil; under a closed canopy the dense form sets the conductance ($W=0.0025$
at LAI + WAI = 6). CLM blends the transfer coefficients that multiply one velocity scale; with both
terms on $`u_*`$, that is this linear blend of the conductances. Snow does not bury the canopy here:
CLM's burial would raise $W$ through the exposed leaf and wood area, and waits for a snowy site (see
docs/ROADMAP.md). The **soil-evaporation aerodynamic resistance** is $`r_{aero}=1/g_{g,net}`$, and
$`g_{g,net}`$ is the ground sensible/latent conductance in the surface energy balance.

## 4. In-canopy wind and the leaf/wood boundary layers

**Canopy-top wind** comes from the log profile with the surface-layer stability correction
(`reduced_wind`):

```math
u_h = \frac{u_*}{\kappa}\big[\ln\tfrac{h-d}{z_{0m}} - \Psi_m(\zeta_h) + \Psi_m(\zeta_0)\big] \qquad(7)
```

floored at `ugbmin`. From the top cohort down, the wind attenuates through each cohort's crowns by
ED2's per-cohort extinction, here with **closed crowns**: every cohort's crowns cover its patch (crown
area 1, as in the two-stream radiation, `canopy_radiation_transfer.md` §3), so the wind meets all of the
cohort's leaf area. Cohorts of one height are one layer, as they are in the radiation
(`canopy_radiation_transfer.md` §4): LAI below is the layer's summed leaf area, every member meets the
layer's wind, and a cohort split into two of half the density meets the wind the whole one does. The
wind at a layer's mid-crown (also floored at `ugbmin`) and the value handed to the layer below are

```math
u_{ico} = u_h\,e^{-0.25\,\mathrm{LAI}}, \qquad
u_h \leftarrow u_h\,e^{-0.50\,\mathrm{LAI}} \qquad(8)
```

so the cascade walks **top → bottom**. The kernel's cohort arrays run bottom → top (the radiation
convention) while the cohort block is gathered tallest first, so `aero_bottom_to_top` reverses them at
the call site; the reversal is the identity for $n\le1$.

**Boundary-layer conductance** for each cohort's leaves (flat plate) and wood (cylinder) is ED2's
Nusselt kernel `boundary_gbh_mos`, the parallel sum of forced (Reynolds) and free (Grashof) convection
— the free term matters at low wind. With Reynolds $`\mathrm{Re}=u\,\ell/\alpha_{th}`$ and Grashof
$`\mathrm{Gr}=g\,|T_{elem}-T_{CAS}|\,\ell^3/(T_{CAS}\,\nu^2)`$ over a characteristic length $\ell$
(leaf width or branch diameter). $`T_{elem}`$ is the leaf temperature for the leaf layer and the wood
temperature for the wood layer. For the cylinder each forced term carries an offset,
$`o+a\,\mathrm{Re}^{\,n}`$ (`ocyli_*`).

```math
\mathrm{Nu}_{forced}=\max\!\big(a_{lam}\mathrm{Re}^{\,n_{lam}},\,a_{turb}\mathrm{Re}^{\,n_{turb}}\big),
\qquad
\mathrm{Nu}_{free}=\max\!\big(b_{lam}\mathrm{Gr}^{\,m_{lam}},\,b_{turb}\mathrm{Gr}^{\,m_{turb}}\big)
\qquad(9)
```

```math
g_{bh} = \max\!\Big(g_{bh,min},\ \big(\mathrm{Nu}_{forced}+\mathrm{Nu}_{free}\big)\,
\frac{\alpha_{th}}{\ell}\Big) \qquad(10)
```

The kinematic viscosity $\nu$ and thermal diffusivity $\alpha_{th}$ are linear-in-temperature about a
20 °C reference (`kin_visc0`, `dkin_visc`, `th_diff0`, `dth_diff`, `t_ref_air`), and
$`g_{bh,min}=`$ `gbhmos_min`. The vapour conductance is $`g_{bw}=1.075\,g_{bh}`$ (`gbh_2_gbw`).
$g_{bh}$ is a velocity [m s⁻¹] — the energy kernel consumes it directly, and the leaf gas-exchange seam
forms the molar $`g_{b,mol}=g_{bw}\,\rho_{mol}`$, with $`\rho_{mol}=p/(R\,T_{leaf})`$ the molar density of air at the leaf temperature [mol m⁻³].

## Parameters (config names, `[aerodynamics]` / `aero_cfg_t`)

| Symbol | Config key | Meaning |
|---|---|---|
| $`\kappa`$ | — (`vonkarman`, `meds_constants`) | von Karman constant (0.4) |
| $`c_{z0}`$ | `z0m_ratio` | roughness / canopy height (0.13) |
| $`c_d`$ | `d_ratio` | displacement / canopy height (0.63) |
| $`z_{0,snow}`$ | `snow_rough` | snow-surface roughness, also the $`z_{0m}`$ floor [m] |
| — | `n_iter_mo` | fixed Monin-Obukhov iterations (4, GPU-uniform) |
| $`\zeta_m,\zeta_t`$ | `zeta_m`, `zeta_t` | momentum / scalar range transitions (CLM5) |
| — | `zeta_max_stable` | cap on the stable initial guess of $\zeta$ (CLM5) |
| $w_c$ | `wc` | convective velocity scale (MoninObukIni) |
| $`u_{*,min}`$ | `ustmin` | friction-velocity floor |
| $`u_{b,min}`$ | `ubmin` | reference-wind floor in eq. (2) |
| — | `ugbmin` | floor on the canopy-top and in-canopy wind |
| $`c_s,\gamma_g`$ | `cs_dense`, `gamma_g` | CLM4 dense-canopy ground conductance |
| — | `gbh_2_gbw` | boundary-layer heat:vapour ratio (1.075) |
| $`g_{bh,min}`$ | `gbhmos_min` | boundary-layer conductance floor [m s⁻¹] |
| — | `kin_visc0`, `dkin_visc`, `th_diff0`, `dth_diff`, `t_ref_air` | linear-in-T air viscosity and diffusivity |
| — | `aflat_*/bflat_*` | flat-plate (leaf) Nusselt coefficients (`aero_cfg_t` only, no TOML key) |
| — | `ocyli_*/acyli_*/bcyli_*` | cylinder (wood) Nusselt coefficients (`aero_cfg_t` only, no TOML key) |
| — | `min_canopy_depth` | CAS-depth floor [m] |
| — | `canopy_freeboard` | CAS depth above the tallest cohort [m] |

## References
- ED2 `../ED2/ED/src/dynamics/canopy_struct_dynamics.f90` (`canopy_turbulence8`) — boundary layers +
  in-canopy wind.
- CLM5 Technical Note (Oleson et al. 2018) — Monin-Obukhov surface layer, stability functions.
- Businger et al. (1971); Paulson (1970) — surface-layer flux-profile relations.
- Zeng, X., et al. (2005), Treatment of undercanopy turbulence in land models, *J. Climate* 18,
  5086–5094 — the open-fraction blend of eq. (6); Sakaguchi, K., and X. Zeng (2009), *J. Geophys.
  Res.* 114, D01107 — its under-canopy stability factor.
- `docs/dev_plans/archive/MEDS_COLUMN_DYNAMICS_DESIGN.md` (Part I, ED2↔CLM reconciliation);
  `docs/dev_plans/archive/MEDS_ENERGY_BALANCE_DESIGN.md`.

## Code map

| Concept | Routine (`meds_canopy_aerodynamics`) |
|---|---|
| master per-patch seam | `canopy_aerodynamics` |
| roughness + displacement, eq. (1) | `canopy_roughness` |
| open fraction $W$, eq. (6) | `open_canopy_fraction` |
| Monin-Obukhov solve ($`u_*`$, temp1, $\zeta$) | `mo_surface_layer` |
| profile denominators $D_m, D_h$ | `d_mom`, `d_heat` |
| stability functions $`\Psi_m,\Psi_h`$ | `stabfunc1`, `stabfunc2`, `psim` |
| canopy-top / in-canopy wind | `reduced_wind` + closed-crown extinction (inline) |
| leaf/wood Nusselt boundary layer | `boundary_gbh_mos` |
| atm↔CAS conductances, eq. (5) | `cas_atm_conductances` |
| per-stage re-solve of eq. (5) | `meds_fast_time_derivs`: `refresh_cas_conductances` |
| types (cfg / env / geom / out) | `meds_biophysics_opts`: `aero_cfg_t`; `meds_canopy_types`: `aero_env_t`, `aero_geom_t`, `aero_out_t` |
| bottom→top cohort reversal | `meds_fast_prepass`: `aero_bottom_to_top` |
