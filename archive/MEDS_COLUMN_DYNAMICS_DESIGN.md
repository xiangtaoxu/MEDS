# MEDS Fast-Loop Coupling — Canopy Aerodynamics + Column Dynamics + Numerics Consolidation

**Status:** DESIGN ONLY (2026-07-06). This is the **P3 capstone** that the energy-balance, column-hydrology,
and column-CO₂ design docs each explicitly deferred to ("the coupled leaf↔CAS↔ground↔soil fixed point is
deferred to P3"; "canopy-air-space + leaf-boundary-layer are deferred sibling modules"). It specifies the
two files the user requested plus the numerics review:

1. **`src/biophysics/meds_canopy_aerodynamics.f90`** — a *stateless* kernel computing all in-canopy
   aerodynamics (friction velocity `ustar`, the similarity scales `tstar/qstar/cstar`, the in-canopy wind
   profile, per-cohort leaf/wood boundary-layer conductances, and the ground↔canopy-air conductance).
2. **`src/driver/meds_column_dynamics.f90`** — the fast-timescale integrator (timestep **`dt_fast`**, never
   `dt_lsm`) that couples every fast process into one sub-daily loop.
3. **`src/shared/meds_numerics.f90`** — the outcome of the solver-consolidation review (§4).

The scientific reference is ED2 (`../ED2/ED/src/dynamics/canopy_struct_dynamics.f90`, the `rk4_*` integrator
family); CLM5/CTSM (`FrictionVelocityMod`, `CanopyFluxesMod`) is the co-reference for the Monin–Obukhov
surface layer and the leaf-temperature fixed point. Modernization follows the MEDS rules (stateless pure
kernels, derived-type config passed by argument, `real(wp)`, ≤132 cols, GPU-portable, no hard-coded
parameters).

---

## 0. What this closes — the "missing supplier" framing

Every fast kernel MEDS already ships (radiation, leaf gas exchange, plant hydraulics, respiration, the three
canopy-air-space *twins* for energy/water/CO₂, soil thermal, soil hydrology) is **stateless** and takes its
turbulence/conductance inputs as **forced arguments with placeholder defaults**. The aerodynamics module's
entire job is to *replace those placeholders with physics*; the column-dynamics driver's entire job is to
*sequence the kernels and weave their state* each `dt_fast`. Nothing new is invented at the process level —
the two files close a wiring gap.

The single most important structural fact, confirmed by reading the shipped kernels: **each fast kernel is a
single implicit / linearized-backward-Euler step** — *mostly* L-stable. `veg_energy_balance` is a linearized
BE step; `canopy_air_update`'s **enthalpy** branch and `canopy_air_co2_update` are implicit-in-the-atmosphere
term (L-stable); `soil_energy_flux` and `column_hydrology_flux` are BE-Thomas with their own adaptive
substepping; the plant hydraulics kernel is a 2×2 matrix-exponential with step-doubling. So MEDS's native fast
integrator is a **first-order operator-split (Lie–Trotter) sweep of per-store implicit steps**, *not* ED2's
global adaptive RK4 over a monolithic state vector. This is cleaner and GPU-friendly — and it is what the
shipped kernels were built for.

**One exception the design must fix before it is uniformly stable:** the shipped `canopy_air_update`
**humidity** twin is *explicit* in the atmospheric exchange (`shv_new = cas_shv + dt·wci·(… + w_flux_ac)`),
so it is only conditionally stable — forward-Euler needs `dt·gaw/wcapcan ≲ 2`, i.e. `dt_fast ≲ can_depth/ustar
≈ 5 m / 0.3 m s⁻¹ ≈ 17–33 s`, far below the target `dt_fast = 300–900 s`. §3.5 makes this twin implicit (the
same form as the enthalpy/CO₂ twins) so the whole set is L-stable at the target step; that fix is **required in
the P3c MVP**, not deferred. Until then, the "unconditionally stable at large `dt_fast`" property holds only
for the enthalpy/CO₂/soil/veg stores. The column driver honours the split; §3.5 closes the last gap.

### 0.1 The forced-input → aerodynamics-output map (the master seam)

Every row is a field a shipped kernel currently receives as a placeholder; the right column is what
`meds_canopy_aerodynamics` must produce. This table *is* the module's requirements spec.

| Consumer kernel / field | Symbol | Units | Current placeholder | Aerodynamics supplies |
|---|---|---|---|---|
| `cas_atm_forcing_t%ustar`; `canopy_air_update(ustar)` | u⋆ | m/s | `0.0` | friction velocity (Monin–Obukhov) |
| atm↔CAS **scalar** conductance (all 3 twins) | g_ah / g_aw / g_ac | kg or mol /m²/s | placeholder `ρ·ustar` | `ρ·ustar·temp1` (heat), `ρ·ustar·temp2` (vap), `ρ_dmol·ustar·temp2` (CO₂) — the profile-factored form, §2.3/§3.5 |
| `canopy_air_update(w_flux_ac)` (vapour) | E_ac | kg/m²/s | `0.0` | `g_aw·(q_atm − can_shv) = ρ·ustar·qstar` (consistent; fold into implicit form, §3.5) |
| `column_co2_step(ustar)`; `gatm_co2 = c_dmol·ustar·temp2` | u⋆ / c⋆ | m/s | `0.0` | shared u⋆ × scalar profile factor; c⋆ sets the CO₂ exchange |
| `leaf_energy_env_t%gbh` (per cohort, leaf & wood) | g_bh | m/s | `0.0` | leaf/wood boundary-layer **heat** conductance |
| `leaf_energy_env_t%gbw` (per cohort, leaf & wood) | g_bw | m/s | `0.0` | leaf/wood boundary-layer **vapour** conductance |
| `leaf_env_t%gb` (leaf gas exchange) | g_b | **mol H₂O/m²/s** | `0.0` (⇒ Cs=Ca) | `gbw · ρ_mol` (molar bridge, §2.7) |
| `ground_surface_balance` `env%gbh` | ggnet | m/s | — | ground↔CAS aerodynamic conductance |
| `ground_surface_balance` `env%gbw` | ggnet | m/s | — | ground vapour conductance (= ggnet) |
| `chydro_forcing_t%r_aero` (soil-evap series) | r_aero | **s/m** | `100.0` | `1 / ggnet` |
| `cas_state_t%can_depth`; `canopy_air_update(can_depth)` | z_CAS | m | `20.0` | from canopy top height (floored, §2.5) |
| `leaf_env_t%vpd` (leaf gas exchange) | VPD | Pa | — | driver: `esat(leaf_temp) − e(can_shv,press)` |

Radiation forcing is **not** aerodynamics' job: `cosz` is computed from the time dimension (`solar_cosz`,
§3.7), the beam/diffuse split from `met%rshort + cosz`, and the ground reflectance/emission from
`meds_optics%ground_optics`. Rhizosphere conductance (`hydro_env_t%rhizo_cond`) is the soil/root module's, not
aerodynamics'.

---

# PART I — `meds_canopy_aerodynamics.f90`

## 1. Scientific scope and the ED2 ↔ CLM reconciliation

The kernel is a pure function of *(free-atmosphere forcing, canopy-air-space state, canopy geometry)* → the
"star family" + conductances. It carries **no integrated state** (ED2 recomputes it every RK stage; MEDS
recomputes it every `dt_fast` substep). The two references disagree on details; MEDS takes the best of each:

| Aspect | ED2 (`canopy_struct_dynamics`) | CLM5 (`FrictionVelocity`) | **MEDS decision** |
|---|---|---|---|
| ζ = z/L solve | Newton + bisection root-find on Rib (`zoobukhov8`), data-dependent trip count, `write`+`fatal_error` fallback | `MoninObukIni` analytic Rib guess + **fixed few-sweep** fixed-point on L | **Adopt CLM**: analytic init + fixed 3–4 iterations, no `while`, no I/O → GPU-portable, deterministic |
| ψ_m / ψ_h | 2 additive fns, 4-region CLM04 with awkward −1/6 exponents | Clean 4-region: unstable Businger–Dyer `x=(1−16ζ)¼`, stable linear `−5ζ`, very-stable log cap; fixed `ζm=1.574, ζt=0.465` | **Adopt CLM's piecewise** (§2.3) — monotone, well-behaved |
| z0h : z0m | ZD98 `exp(−a·Re^0.45)` (opt.) else 1 | canopy `z0h,v = z0m,v`; **ground** z0h via ZD98 `a=0.13, exp 0.45` | **Adopt CLM split**: canopy `z0h=z0m`; ground thermal-roughness = ZD98 |
| Leaf `gb` | full Nu/Re/Gr free+forced convection, per-PFT leaf width | `rb = (1/Cv)√(dleaf/Uav)`, Cv=0.01 | **Keep ED2's Nu/Re/Gr** — captures free convection at low wind (essential for leaf energy balance); expose CLM `rb` as a cheap fallback |
| In-canopy wind | Leuning per-cohort exponential, crown-area-limited | single-layer `Cs·Uav` | **Keep ED2** per-cohort extinction (needed for multi-cohort gas exchange) |
| Canopy-air point | **prognostic** CAS (θ,q,CO₂ + air-mass capacities) | **diagnostic** zero-storage junction | **Keep prognostic CAS** — MEDS already built it; CLM's junction is its `wcapcan→0` limit (§3.4) |

Net: **CLM math for the surface layer, ED2 math for the canopy interior, MEDS's prognostic CAS for the
reservoir.** `icanturb`-style layered closures (Massman 1997/1999, which need a 100-layer scratch grid) are
*out of scope for the first cut* — the seam is kept wide enough to add them later (§2.9).

## 2. Module structure — pure-function decomposition

Design rule: **per-patch quantities** are scalar `pure` functions/subroutines; **per-cohort** quantities are
`pure` kernels that write into caller-owned SoA slices (never return an array-valued function result into a
call — the nvfortran issue-#7 trap). One master driver-facing seam calls the pieces in order.

```
meds_canopy_aerodynamics            ! master, per patch → fills aero_out_t
 ├─ aero_roughness      (pure fn)   ! z0m, displacement d, z0h(ground) from height + LAI + snow blend
 ├─ aero_surface_layer  (pure sub)  ! CLM Monin–Obukhov: ustar, tstar, qstar, cstar, zeta, rib, obu, ggbare
 │    └─ psim / psih    (pure elem) ! CLM 4-range stability functions (§2.3)
 ├─ aero_reduced_wind   (pure fn)   ! wind at any height from (ustar, zeta) — used for canopy-top wind
 ├─ aero_wind_profile   (pure sub)  ! per-cohort in-canopy wind by crown-area extinction (SoA slice out)
 ├─ aero_leaf_boundary  (pure elem) ! per-cohort leaf gbh/gbw (flat-plate Nu/Re/Gr)   → SoA slice
 ├─ aero_wood_boundary  (pure elem) ! per-cohort wood gbh/gbw (cylinder Nu/Re/Gr)     → SoA slice
 ├─ aero_ground_cond    (pure fn)   ! ggveg (icanturb 0/4) + ggnet blend
 └─ aero_cas_capacity   (pure sub)  ! wcapcan/hcapcan/ccapcan (+ inverses) from ρ, ρ_mol, can_depth
```

### 2.1 The master seam

```fortran
pure subroutine meds_canopy_aerodynamics(cfg, env, geom, coh, n, out)
   type(aero_cfg_t),   intent(in)  :: cfg          ! run constants (device-constant; §2.8)
   type(aero_env_t),   intent(in)  :: env          ! free-atm forcing + CAS state + ground state
   type(aero_geom_t),  intent(in)  :: geom         ! per-patch canopy geometry (height, PAI, opencan, snowfac)
   integer(ik),        intent(in)  :: n            ! cohort count
   type(aero_cohort_t),intent(in)  :: coh          ! SoA slices (in): height, lai, wai, crown_area,
                                                    !   leaf_temp, wood_temp, pft, leaf_width, branch_diam
   type(aero_out_t),   intent(inout) :: out        ! preallocated; gbh/gbw/wind SoA slices written in place
end subroutine
```
`intent(inout)` on `out` (not `out`) so the caller owns the cohort arrays (`gbh(:)`, `gbw(:)`, `wind(:)`);
the routine writes them in place — the GPU-safe pattern.

### 2.2 Roughness, displacement, snow blend (per patch)

```
snowfac      = min(0.99, total_sfcw_depth / veg_height)            ! snow burial fraction
rough (z0m)  = snow_rough·snowfac                                                     &
             + (soil_rough·opencan_frac + veg_rough·(1−opencan_frac))·(1−snowfac)
veg_rough    = cfg%z0m_ratio · veg_height          ! ED2 vh2vr=0.13  (CLM Rz0m 0.055–0.12 per PFT, optional)
displace (d) = cfg%d_ratio   · rough / cfg%z0m_ratio ! ED2 vh2dh=0.63  →  d = 0.63·veg_height nominal
z0h_ground   = z0m_ground · exp(−cfg%zd98_a · (ustar·z0m_ground/nu)^cfg%zd98_b)   ! ZD98, a=0.13, b=0.45
```
Canopy `z0h = z0m` (CLM); only the **ground** gets a thermal-roughness reduction (ZD98). `nu` (kinematic
viscosity) is linear in `can_temp`.

### 2.3 The Monin–Obukhov surface layer (CLM piecewise, fixed-iteration)

`zeta = zldis/L`, `zldis = zref − d`; transition constants `ζm = 1.574` (momentum), `ζt = 0.465` (heat/vapour).
Unstable helpers (`x = (1−16ζ)^¼`):
```
StabilityFunc1(ζ) = 2·ln((1+x)/2) + ln((1+x²)/2) − 2·atan(x) + π/2      ! ψ_m, unstable
StabilityFunc2(ζ) = 2·ln((1+x²)/2)                                       ! ψ_h, unstable
```
Momentum profile denominator `D_m` (⇒ `ustar = vonk·um/D_m`), per range:

| Range | Condition | `D_m` |
|---|---|---|
| very unstable | ζ < −ζm | `ln(−ζm·L/z0m) − F1(−ζm) + F1(z0m/L) + 1.14·((−ζ)^⅓ − ζm^⅓)` |
| unstable | −ζm ≤ ζ < 0 | `ln(zldis/z0m) − F1(ζ) + F1(z0m/L)` |
| stable | 0 ≤ ζ ≤ 1 | `ln(zldis/z0m) + 5ζ − 5·z0m/L` |
| very stable | ζ > 1 | `ln(L/z0m) + 5 − 5·z0m/L + (5·ln ζ + ζ − 1)` |

Heat/vapour denominator `D_h` (⇒ `temp1 = vonk/D_h`, and `temp2` identically with `z0q`), transition `ζt`:

| Range | Condition | `D_h` |
|---|---|---|
| very unstable | ζ < −ζt | `ln(−ζt·L/z0h) − F2(−ζt) + F2(z0h/L) + 0.8·(ζt^−⅓ − (−ζ)^−⅓)` |
| unstable | −ζt ≤ ζ < 0 | `ln(zldis/z0h) − F2(ζ) + F2(z0h/L)` |
| stable | 0 ≤ ζ ≤ 1 | `ln(zldis/z0h) + 5ζ − 5·z0h/L` |
| very stable | ζ > 1 | `ln(L/z0h) + 5 − 5·z0h/L + (5·ln ζ + ζ − 1)` |

Initialization (`MoninObukIni`), then a **fixed `cfg%n_iter_mo` (=3–4) sweeps** — no data-dependent exit, for
warp uniformity:
```
ustar0 = 0.06 ; wc = 0.5
um  = (dthv≥0) ? max(ur,0.1) : sqrt(ur² + wc²)                 ! add convective velocity when unstable
rib = grav·zldis·dthv / (thv·um²)                              ! bulk Richardson first guess
ζ   = (rib≥0) ? clip(rib·ln(zldis/z0m)/(1−5·min(rib,0.19)), 0.01, cfg%zeta_max_stable)   ! CLM5 zetamaxstable = 0.5
              : clip(rib·ln(zldis/z0m),               −100.0, −0.01)
L   = zldis/ζ
do i = 1, n_iter_mo
   ustar = max(cfg%ustmin, vonk·um / D_m(ζ))
   temp1 = vonk / D_h(ζ) ;  temp2 = vonk / D_h_q(ζ)
   tstar = temp1·dth ;  qstar = temp2·dqh ;  cstar = temp2·dco2      ! CO₂ uses the vapour profile
   thvstar = tstar·(1+0.61·q_atm) + 0.61·θ_atm·qstar
   ζ = zldis·vonk·grav·thvstar / (ustar²·thv) ;  L = zldis/ζ
   um = (dthv<0) ? sqrt(ur² + (cfg%wc·wstar)²) : max(ur,0.1)          ! optional convective update
end do
ggbare = ustar·temp1                                                 ! ground↔CAS scalar conductance [m/s] (= 1/rah)
```
`dth = θ_atm − can_theta`, `dqh = q_atm − can_shv`, `dco2 = co2_atm − can_co2`, `dthv` the virtual-temp
difference. Clamps to port: `ustmin ≈ 0.10`, stable ζ∈[0.01, `zeta_max_stable`=0.5], unstable ζ∈[−100,−0.01],
Rib guess capped at 0.19. `zeta_max_stable` is CLM5's `zetamaxstable` (default 0.5, a config value — *not* the
CLM4-era 2.0); the clamp applies to the init guess, not re-clamped inside the iteration. **All in `real(dp)`**
(the Rib/Obukhov solve is ill-conditioned near ζ→0 in single precision).

**The atm↔CAS scalar conductances (the critical output for the CAS twins).** The four stars feed the CAS
budgets: momentum via `ustar`, sensible via `tstar`, latent via `qstar`, CO₂ via `cstar`. The atm↔CAS
exchange conductance is *not* bare `ρ·ustar` — it carries the same scalar profile factor `temp1`/`temp2` the
ground conductance already uses (this is ED2's `estar = c3·(h_atm−h_can)`, `eflxac = ρ·ustar·estar` form,
`c3 ≡ temp1`). The module therefore returns:
```
gah = ρ_air ·ustar·temp1    [kg/m²/s]   ! atm↔CAS ENTHALPY conductance
gaw = ρ_air ·ustar·temp2    [kg/m²/s]   ! atm↔CAS VAPOUR  conductance   (temp2 with z0q; ≈ temp1 when z0q=z0h)
gac = ρ_dmol·ustar·temp2    [mol/m²/s]  ! atm↔CAS CO₂     conductance (molar)
w_flux_ac = gaw·(q_atm − can_shv) = ρ_air·ustar·qstar                  ! explicit vapour flux (self-consistent)
```
The shipped `canopy_air_update`/`canopy_air_co2_update` hard-code `gatm = ρ·ustar` — a **placeholder that drops
`temp1`/`temp2` and so overstates the atm↔CAS coupling by ~5–10×** (for a typical `(zref−d)/z0h`,
`temp1 ≈ 0.1–0.2`), pinning the CAS far too tightly to the free atmosphere. **Correcting the twins to consume
the aerodynamics-supplied `gah`/`gaw`/`gac` (or equivalently `ggbare`) is a required part of this coupling**
(§3.5), landing in P3c alongside the implicit-vapour fix.

### 2.4 Leaf & wood boundary-layer conductances (ED2 Nu/Re/Gr)

Per cohort, from the in-canopy wind at crown mid-height `veg_wind` (§2.5). Monteith & Unsworth ch.10:
```
nu     = kin_visc0·(1 + dkin·(can_temp − T0)) ; α_th = th_diff0·(1 + dth·(can_temp − T0))
! FORCED (flat plate for leaf; cylinder for wood):
Re     = veg_wind·L_char / α_th
Nu_frc = max(a_lami·Re^n_lami , a_turb·Re^n_turb)                     ! flat: a_lami=0.60,n=0.5 / turb 0.032,0.8
! FREE (buoyant):
Gr     = (grav/(can_temp·nu²))·|elem_temp − can_temp|·L_char³
Nu_free= max(b_lami·Gr^m_lami , b_turb·Gr^m_turb)                     ! flat: b_lami=0.50,m=0.25 / turb 0.19,1/3
gbh_mos= max(gbhmos_min, (Nu_frc + Nu_free)·α_th / L_char)            ! [m/s]  the primitive conductance
gbh    = gbh_mos·can_rhos·can_cp    [J/K/m²/s]   ;   gbw = cfg%gbh_2_gbw·gbh_mos·can_rhos    [kg/m²/s]
```
`L_char = leaf_width(pft)` for leaves (flat plate), `branch_diam(pft)` for wood (cylinder Nu coefficients).
`gbh_2_gbw = 1.075` (Lewis-number heat:vapour ratio). **MEDS emits the primitive `gbh_mos [m/s]`** plus the
temperature/density it used, and lets each consumer apply its own convention (§2.7) — cleaner than shipping
three unit variants. The energy kernel wants `gbh, gbw` in **m/s** (its `leaf_energy_env_t%gbh/gbw` are
velocities: flux `= gbh·area·ρ·cp·ΔT`); the leaf gas-exchange kernel wants **molar** `gb = gbw·ρ_mol`.

**Resolvability gate:** ED2 skips the Nu/Re/Gr call for cohorts below an LAI/WAI threshold (or snow-buried)
and sets `gbw = f_bndlyr_init·gsw` (10×, so the boundary layer never limits). MEDS carries per-cohort
`leaf_resolvable`/`wood_resolvable` flags (or reproduces the fallback) — see §3.2.

### 2.5 In-canopy wind profile + ground conductance (ED2 icanturb 0/4)

Canopy-top wind from the log profile, then per-cohort exponential extinction, cohorts sorted tall→short
(a scalar `uh` carried down — a per-patch serial reduction, cheap, not vectorized across cohorts):
```
uh = aero_reduced_wind(ustar, ζ, rib, zref, d, veg_height, rough)    ! wind at canopy top, floored at ugbmin
do ico = tallest, shortest
   ext_half     = crown_area·exp(−0.25·lai/crown_area) + (1−crown_area)
   ext_full     = crown_area·exp(−0.50·lai/crown_area) + (1−crown_area)
   veg_wind(ico)= max(ugbmin, uh·ext_half)                           ! wind at crown mid-depth
   uh           = uh·ext_full                                        ! attenuate for the next cohort down
end do
```
Ground conductance (CLM-like `icanturb=4`, the recommended default):
```
stab   = clip( 2·grav·veg_height·(can_temp − t_ground) / ((can_temp + t_ground)·ustar²), 0, 10 )
ggveg  = cfg%cs_dense·ustar / (1 + cfg%gamma_g·stab)                 ! cs_dense=0.004, gamma_g=0.5
ggnet  = (opencan_frac > 0.999 .or. snowfac ≥ 0.9)                                                  &
         ? ggbare                                                                                    &
         : ggbare·ggveg / (ggveg + (1−opencan_frac)·ggbare)          ! series-parallel: 1/ggnet = 1/ggbare + (1−oc)/ggveg
```
`ggnet [m/s]` is the single ground↔CAS conductance that serves three consumers with reciprocal conventions:
`ground_surface_balance%gbh = ggnet` (m/s), `ground_surface_balance%gbw = ggnet` (m/s), and
`chydro_forcing_t%r_aero = 1/ggnet` (s/m). **This equivalence is a hard contract** (§3.6) — it is the fix for
the double-counted soil evaporation.

CAS depth: `can_depth = max(cfg%min_canopy_depth, veg_height)` (ED2 `minimum_canopy_depth = 5 m`).

### 2.6 CAS storage capacities

```
wcapcan = hcapcan = can_rhos·can_depth        [kg air/m²]           ! enthalpy & vapour capacity
ccapcan = can_dmol·can_depth                  [mol air/m²]          ! CO₂ (molar) capacity
```
The shipped `canopy_air_update`/`canopy_air_co2_update` currently recompute these internally from
`rho_air, can_depth` and `can_dmol = ρ(1−q)/M_dry`. The aerodynamics module can *also* return them (and their
inverses) so a future consolidated CAS update consumes them directly; for the MVP they stay internal to the
twins and aerodynamics just supplies `ustar` + `can_depth` + `rho_air`.

### 2.7 Units and the two conductance conventions (the bridge helpers)

Aerodynamics emits the primitive `gbh_mos [m/s]`. The driver bridges to each consumer:

| Consumer | Wants | Bridge (in the driver, via `meds_thermo`/`meds_constants`) |
|---|---|---|
| `leaf_energy_env_t%gbh` | m/s (heat) | `gbh_mos` directly |
| `leaf_energy_env_t%gbw` | m/s (vapour) | `gbh_2_gbw · gbh_mos` |
| `leaf_env_t%gb` | mol H₂O/m²/s | `gbw_ms · ρ_mol`, `ρ_mol = press/(R_univ·can_temp)` [mol/m³] |
| `chydro_forcing_t%r_aero` | s/m | `1 / ggnet` |
| leaf `vpd` | Pa | `sat_vapor_pressure(leaf_temp) − vapour_pressure(can_shv, press)` |

`meds_thermo` already provides `sat_vapor_pressure`, `sat_specific_humidity`, `air_density`,
`cas_temp_of_enthalpy`. Add one small `pure elemental` helper `molar_air_density(press, temp)` (or reuse an
existing constant). The molar bridge is the one genuinely new unit conversion the coupling introduces.

### 2.8 Parameters — config-sourced, never hard-coded

Per the MEDS "no hard-coded model parameters" rule, all constants live in `aero_cfg_t` (filled from the
`[aerodynamics]` TOML block at the aux/init layer; true physical constants like `vonk` may live in
`meds_constants`). Catalog (ED2 defaults / CLM values):

```
vonk 0.40  ·  z0m_ratio 0.13  ·  d_ratio 0.63  ·  exar 2.5 (wind extinction)
ustmin 0.10  ·  ubmin 0.65  ·  ugbmin 0.25  ·  gbhmos_min 1e-9  ·  min_canopy_depth 5.0
zeta_m 1.574  ·  zeta_t 0.465  ·  zeta_max_stable 0.5 (CLM5 zetamaxstable)  ·  n_iter_mo 4  ·  wc 0.5 (convective vel scale)
zd98_a 0.13  ·  zd98_b 0.45  ·  cs_dense 0.004  ·  gamma_g 0.5 (CLM ground)
gbh_2_gbw 1.075  ·  f_bndlyr_init 10.0
kin_visc0 1.33e-5 (slope 0.007)  ·  th_diff0 1.89e-5 (slope 0.007)
Nusselt flat: aflat_lami 0.60 nflat_lami 0.50 / aflat_turb 0.032 nflat_turb 0.80 /
              bflat_lami 0.50 mflat_lami 0.25 / bflat_turb 0.19 mflat_turb 0.333
Nusselt cyl:  ocyli_lami 0.32 acyli_lami 0.51 ncyli_lami 0.52 / ocyli_turb 0.0 acyli_turb 0.24 ncyli_turb 0.60 /
              bcyli_lami 0.48 mcyli_lami 0.25 / bcyli_turb 0.09 mcyli_turb 0.333
per-PFT: leaf_width, branch_diam   (add to meds_pft_params if absent)
```
`isfclyrm`/`icanturb` become an `integer(ik)` selector in `aero_cfg_t` with only the CLM-MO + icanturb-0/4
paths implemented at first (others `error stop 'not implemented'`).

### 2.9 GPU / portability notes specific to aerodynamics

- The whole kernel is `pure` and branch-light. The only non-vectorizable piece is the tall→short `uh` carry
  in `aero_wind_profile` — a per-patch serial scalar reduction, fine on host and device (cohorts within a
  patch are few).
- **No root-finder, no I/O**: the CLM fixed-iteration solve replaces ED2's `zoobukhov8` Newton+bisection —
  whose trip count is *data-dependent* (early-exit on convergence, bounded by `maxfpo`) and which prints
  diagnostics via `write` and calls `fatal_error` on failure. Those (data-dependent divergence + I/O + halt),
  not an unbounded loop, are the portability hazards; the fixed-iteration solve removes all three. The single
  biggest portability win.
- `real(dp)` for the surface-layer arithmetic; `real(wp)` elsewhere is fine.
- Do not pass an array-valued function result into a call (issue #7): `aero_leaf_boundary`/`aero_wood_boundary`
  are `pure elemental` writing into caller slices; `gbh_mos` etc. bind to named locals first.
- Deferred `icanturb 2/3` (Massman) would need a per-patch `ncanlyr` scratch grid — host-only, thread-local,
  add later behind the selector; do **not** let it into the device path.

---

# PART II — `meds_column_dynamics.f90`

## 3. The fast-timescale integrator

`meds_column_dynamics` lives in `src/driver/` (globbed into `meds_aux`, the ED2-`ed_model` analogue). It owns
the fast loop: for each patch, run `n_fast_per_slow` sub-steps of `dt_fast`, calling every fast kernel in
operator-split order, weaving the CAS twins as the shared coupling reservoir. It NEVER references `dt_lsm`.

### 3.1 Where it hooks — cadence and config

Today `meds_main` walks the calendar with a **minimum step of one day**; `meds_stepper%advance_one_step`
calls only the slow loop (`vegetation_dynamics`). The fast loop is a **new cadence layer below the day**:

```fortran
! meds_stepper.f90 — extended
subroutine advance_one_step(site, cfg, met, now, is_new_month, is_new_year)
   ...
   if (cfg%fast_biophysics_on) call column_dynamics(site, cfg, met, now)   ! NEW — fast, before slow
   call vegetation_dynamics(site, cfg, is_new_month, is_new_year)          ! slow, unchanged
end subroutine
```

**Config change — retire `ts_mode`, introduce `dt_slow`.** The current time config is a *string enum*
`[run] timestep = "daily"|"weekly"|"monthly"` → `ts_mode` (TS_DAILY/WEEKLY/MONTHLY) → `dt_years` derived
(`meds_config.f90:60,141` + the `TS_*` parameters). This design **replaces it with a single user-defined
slow-process resolution `dt_slow`, defaulting to 1 day** — the interval at which the slow (demography) loop
advances, and the outer bound of the fast sub-loop:
```fortran
! meds_config_t time-stepping block — REPLACES ts_mode/dt_years enum
real(wp)    :: dt_slow                 ! [s] slow-process timestep (user-defined resolution); DEFAULT 86400 (1 d)
! ── fast-loop additions ([fast] TOML block) ──
logical     :: fast_biophysics_on      ! master gate
real(wp)    :: dt_fast                 ! [s] fast biophysics timestep (ED2 DTLSM analogue), e.g. 300–900
integer(ik) :: integration_scheme      ! SCHEME_SPLIT_SEQUENTIAL (default) | SCHEME_PICARD_COUPLED (§3.4)
```
TOML: `[run] timestep = "daily"` → `[run] dt_slow = "1d"` (a duration string parsed to seconds — `"1d"`,
`"7d"`, `"6h"`; default `"1d"`). Derived in `derive_config` (using `yr_sec` from `meds_constants`):
```fortran
cfg%dt_years        = cfg%dt_slow / yr_sec                          ! demography rates stay per-year (unchanged)
cfg%n_fast_per_slow = max(1, nint(cfg%dt_slow / cfg%dt_fast))
```
`dt_years` remains the demography currency (rate integration, `growth_window_steps`) — only its *source*
changes (from the `ts_mode` case-switch to `dt_slow/yr_sec`), so the demography engine is untouched. The
`TS_DAILY/WEEKLY/MONTHLY` enums retire; the monthly/annual **fuse-fiss + disturbance cadence is unaffected** —
it is already driven by the calendar `is_new_month`/`is_new_year` rollovers (independent of `dt_slow`, §3.3).
`validate_config`: `dt_slow > 0`, `dt_fast > 0`, `dt_fast ≤ dt_slow`, scheme in range. The sub-daily cursor is
a driver-local `real(wp) :: t_sec` accumulator across the slow step — the calendar `now` need not change until
the day rolls, so `meds_time` needs no change (optionally add `time_advance_seconds` for sub-daily output
stamping; `time_to_decimal_year` already folds `sec_of_day`).

### 3.2 Prognostic state to add (and the lockstep obligation)

None of the fast prognostic stores exist in the SoA yet. This is real state-plumbing work, and every new
per-cohort field **must** be threaded through the *one* lockstep reorder machinery
(`cohort_reorder`/`cohort_compact`/`copy_cohort_slot`/`rebuild_csr`/`set_cohort_size` in
`meds_demography_types`), and every new per-patch field through `sort_patches` + `patch_compact` in
`meds_demography_fusefiss` — the fix for ED2's "forgot to reallocate" class.

**Per patch** (add to `patch_index`, or a parallel per-patch biophysics block referenced by CSR):
| field | type / units | note |
|---|---|---|
| `cas` | `cas_state_t` (can_enthalpy J/kg, can_shv kg/kg, can_co2 µmol/mol; can_temp, can_depth diagnosed/forced) | 3 prognostic twins |
| `soil_w` | `soil_column_t` (theta(1:n), w_surface, w_aquifer, z_wt) | prognostic moisture |
| `soil_e` | `soil_energy_column_t` (soil_energy(1:n) J/m³; soil_temp, soil_fliq diagnosed) | prognostic internal energy |
| `htry_fast` | `real(wp)` | warm-start substep, if any kernel keeps an adaptive `h` at the driver level |

`fast_soil_carbon [kgC/m²]` is read-only in the fast loop (written daily by biogeochem); it lives on the
existing per-patch carbon state.

**Per cohort** (add to `cohort_block`, thread through the reorder):
| field | units | note |
|---|---|---|
| `leaf_store_energy` | J/m² ground | `veg_energy_balance` prognostic (leaf) |
| `wood_store_energy` | J/m² ground | `veg_energy_balance` prognostic (wood) |
| `leaf_water` | kg/m² ground | interception film (`intercept_canopy_layer` inout; wetted fraction σ_w) |
| `psi(N_HYDRO=3)` | MPa | plant water potential (leaf/stem/root), hydraulics prognostic; store as `psi_leaf`,`psi_wood` (2-node) or a `(3,coh)` block |
| `leaf_resolvable`, `wood_resolvable` | logical | boundary-layer resolvability gate (§2.4) |

**Fusion/fission conserving averages** (a new obligation for these stores): energy and water are *extensive*
(J/m², kg/m²) → number/area-weighted sums, consistent with the AGB-conservation discipline; `psi` is
*intensive* (MPa) → carry the survivor's value or a leaf-mass-weighted mean. Every fuse/split must keep the
column's total internal energy and water within the 1% budget the engine already asserts.

**Diagnosed each step (never carried):** `soil_temp`/`soil_fliq` (from `soil_energy` via `uext_to_temp`),
`can_temp` (from enthalpy+shv via `cas_temp_of_enthalpy`), leaf/wood `temp`/`fliq`, `psi_soil` (from theta),
all conductances, stars, and residuals.

### 3.3 The operator-split sweep (per patch, per `dt_fast`)

The default `SCHEME_SPLIT_SEQUENTIAL`: a single Lie–Trotter forward sweep, each kernel seeing the latest
sibling state (Gauss–Seidel), the CAS twins integrated implicitly at the end. Radiation runs at a coarser
cadence (SW once per `dt_fast` or less; net LW updated analytically from evolving leaf temps).

```
── ONCE per dt_fast (or coarser): ──────────────────────────────────────────────
(0) canopy_radiation(opt, forcing, ...) → abs_leaf/abs_wood(band,coh), dn_ground   [B: meds_canopy_radiation]
    forcing.canopy_temp(coh) ← current leaf/wood store temps (the LW-band coupling seam)

── per dt_fast: ────────────────────────────────────────────────────────────────
(1) diagnose CAS thermo: can_temp = cas_temp_of_enthalpy(can_enthalpy, can_shv); ρ, ρ_mol, ρ_dmol
(2) meds_canopy_aerodynamics(cfg, env{CAS state, met wind, t_ground}, geom, coh, n, aero)
        → ustar, tstar/qstar/cstar, temp1/temp2 (scalar profile factors), ggnet, can_depth, veg_wind(:), gbh(:), gbw(:)
        set  gah = ρ·ustar·temp1 ;  gaw = ρ·ustar·temp2 ;  gac = ρ_dmol·ustar·temp2   ! atm↔CAS conductances
             r_aero = 1/ggnet ;  w_flux_ac = gaw·(q_atm − can_shv)  [= ρ·ustar·qstar]
(3) LEAF GAS EXCHANGE  per cohort:  leaf_env{ par←abs_leaf(VIS), leaf_temp←store, vpd←esat(Tleaf)−e(can_shv),
        ca←can_co2, press, psi_leaf←psi(LEAF), gb = gbw·ρ_mol } → leaf_gas_exchange → a_gross, a_net, gs, rd, E
(4) PLANT HYDRAULICS   per cohort:  hydro_env{ transp = E·leaf_area·M_h2o, soil_psi←last psi_soil (root-weighted),
        rhizo_cond, bleaf/bsap/broot/sap_area/height/leaf_area } → solve_plant_water(dt=dt_fast) → psi(:), root_uptake
(5) RESPIRATION        per cohort:  stem_maintenance_respiration, fine_root_maintenance_respiration → stem_resp, root_resp
(6) INTERCEPTION (capture) sweep tall→short: intercept_canopy_layer for the PRECIP/throughfall path only —
        capture q_intr, pass throughfall down; bottom throughfall = precip_ground. (The film-evap DEBIT of
        leaf_water is deferred to (7b), using the current-step energy w_flux — §3.6 contract 1b.)
(7) VEG ENERGY         per cohort, leaf & wood: veg_energy_balance(store_energy, env{abs_sw,abs_lw,can_temp,can_shv,
        gbh,gbw, gsw = gs·(m/s bridge), fs_open, area_index, leaf_water, wmass, dry_hcap, ρ, press}, dt=dt_fast)
        → new store temp, h_flux, qw_flux (film-evap enthalpy), q_transp, w_flux (film-evap mass = AUTHORITY), transp
(7b) FILM MASS UPDATE  leaf_water += (q_intr − w_flux)·dt, clamp [0,w_max]; clamp overflow → budget residual
        (same w_flux the CAS is credited at step 9 — closes the leaf_water↔CAS transfer)
(8) GROUND             ground_surface_balance(t_ground = soil_temp(1), env{gbh=gbw=ggnet, can_temp, can_shv,
        abs_sw, abs_lw}) → g_top, h_ground, le_ground     [and the DSL soil-evap authority, §3.6]
(9) CAS ENTHALPY+VAPOUR  sum cohort + ground fluxes → canopy_air_update(can_enthalpy, can_shv, can_temp,
        can_depth, Σcoh_h, Σcoh_qw+q_transp, Σcoh_w, Σcoh_transp, ground_h, ground_qw, ground_w, dew,
        ustar, enthalpy_atm, w_flux_ac, ρ, dt_fast, resid)      IMPLICIT in the atm term
(10) CAS CO₂           column_co2_step(can_co2, can_depth, can_shv, ustar, co2_atm, ρ, dt_fast, n, nplant,
        leaf_area, a_gross, rd, stem_resp, root_resp, growth_resp≈0, storage_resp≈0, fast_soil_carbon,
        soil_temp, theta, theta_dry, theta_sat, opts, budget)   IMPLICIT in the atm term
(11) SOIL THERMAL      soil_energy_flux(soil_e, forcing{ g_top←(8), geothermal, soil_water←theta, w_flux, root_heat_sink }, dt=dt_fast)
        → soil_energy, diagnose soil_temp(1) → next step t_ground
(12) SOIL HYDRO        column_hydrology_flux(soil_w, forcing{ precip_ground←(6), root_uptake←(4), t_ground,
        q_air = can_shv, ρ, r_aero = 1/ggnet }, dt=dt_fast) → theta, psi_soil (→ next step (4)), soil_evap, drainage
── accumulate for the slow loop: GPP, NPP, transpiration, ... (§3.8); check budget residuals ──
```

The one-step lag (`psi_soil`, `t_ground`, `q_air`, `e_canopy` lag by one `dt_fast`) is intentional and
physically benign at minute-scale `dt_fast` — it is exactly the "each kernel takes sibling temps as forced
inputs" contract the shipped kernels were designed around. Each kernel's *own* internal substepping
(soil BE step-doubling, hydraulics matrix-exp step-doubling) handles its private stiffness; the driver only
synchronizes at `dt_fast` boundaries. This *loose coupling* is deliberately GPU-friendlier than a global
lockstep (each kernel keeps its own stiffness control).

### 3.4 Coupling-architecture decision — prognostic-CAS split vs CLM inner fixed point

**Recommendation: OUTER = MEDS prognostic-CAS operator split (the default above); INNER = optional CLM-style
leaf fixed point.** Rationale, cross-checked against CLM's `CanopyFluxes`:

- MEDS already committed to **prognostic CAS** — three twins with implicit-atmosphere L-stable updates. CLM's
  *diagnostic* zero-storage canopy-air junction (`taf/qaf` as conductance-weighted means, re-solved every
  Newton iteration) is literally the `wcapcan→0` degenerate limit of MEDS's updates. Adopting CLM's junction
  would be a regression in generality and would throw away the closed-budget residuals (`resid`,
  `budget%resid ≈ 0`) MEDS relies on for verification. **Do not introduce a diagnostic `taf/qaf`** — the
  prognostic twin *is* the canopy-air point (`can_temp`, `q_air = can_shv`, `ca = can_co2` are its diagnostic
  reads).
- **But the leaf sub-problem is stiff** (leaf temp ↔ stomata ↔ leaf VPD/Cs ↔ boundary layer), exactly what
  CLM's ≤40-iteration Newton loop solves. MEDS's `veg_energy_balance` currently takes `can_temp/can_shv/gbh/
  gbw/gsw` as *forced* siblings (the P3-deferred fixed point). The opt-in `SCHEME_PICARD_COUPLED` wraps
  steps (3)→(7) in an inner iteration **with the CAS frozen**, re-linearizing leaf temp and re-solving
  stomata until convergence, borrowing CLM's constants verbatim: `itmax ≈ 40`, `itmin = 2`,
  `dtmin = 0.01 K` (|ΔT_veg|), `dlemin = 0.1 W/m²` (|ΔLE|), and the `nmozsgn ≥ 4` Monin–Obukhov sign-flip
  guard (freeze `L = zldis/(−0.01)`). This closes the *fast* leaf coupling within a substep while the CAS
  reservoir carries the *slow* coupling forward implicitly.
- `ustar` is computed **once per substep** and frozen through the CAS updates (matching how
  `canopy_air_update`/`canopy_air_co2_update` already take `ustar` as a scalar). CLM re-solves `ustar` every
  leaf iteration *because* its junction is algebraic; MEDS can freeze it because the reservoir carries the
  memory. If leaf temps swing hard within a `dt_fast`, the fix is to shrink `dt_fast` (or enable the Picard
  inner loop) — MEDS's replacement for CLM's inner `ustar` re-solve.

The `integration_scheme` selector is therefore **not** ED2's Euler/RK4/Heun/Hybrid global-integrator switch
(MEDS doesn't need one — each kernel is already implicit); it governs *how tightly the driver couples the
stores across a `dt_fast`*: `SPLIT_SEQUENTIAL` (one Gauss–Seidel sweep) or `PICARD_COUPLED` (outer/inner
iteration of the leaf↔CAS fixed point). A future high-accuracy `SCHEME_GLOBAL_RK` (pack the column into one
state vector, adaptive embedded RK45, re-evaluate aerodynamics per stage — ED2's production path) is
documented as a reserved option (§3.9) but is a large rewrite of the shipped single-step kernels into
RHS-derivative form and is *not* recommended for the first cut.

### 3.5 Unifying the three CAS twins + the correct atm↔CAS conductance (REQUIRED, in P3c)

Two coupled corrections to the shipped CAS twins, both enabled by aerodynamics and both **required for a stable,
physical MVP** — not optional refinements.

**(a) Use the profile-factored conductance, not bare `ρ·ustar`.** The shipped twins hard-code the atm↔CAS
conductance as `gatm = ρ·ustar` (enthalpy: `enth_new = (enth + dt·wci·(f + gatm·enth_atm))/(1 + dt·wci·gatm)`;
CO₂ analogously with `c_dmol·ustar`). As §2.3 shows, the physically faithful conductance carries the scalar
profile factor `temp1`/`temp2` (ED2's `c3`): `gah = ρ·ustar·temp1`, `gac = ρ_dmol·ustar·temp2`. Bare `ρ·ustar`
overstates the coupling by ~5–10×. The driver now passes the aerodynamics-supplied `gah`/`gac` in place of the
internally-reconstructed `ρ·ustar` — a small edit to both twins (accept the conductance, or accept `ggbare`
and multiply by ρ / ρ_dmol).

**(b) Make the vapour twin implicit like the other two.** The shipped `can_shv` update is *explicit* in the
atm term (`shv_new = shv + dt·wci·(… + w_flux_ac)`) → only conditionally stable (`dt_fast ≲ can_depth/ustar
≈ 17–33 s`, §0). Make it implicit, identical in form to enthalpy/CO₂:
```
shv_new = (can_shv + dt·wci·(Σsurface_vapour + gaw·q_atm)) / (1 + dt·wci·gaw),   gaw = ρ·ustar·temp2
```
This is a **signature change to `canopy_air_update`**, not a one-line internal edit: the routine currently
receives `enthalpy_atm` and `w_flux_ac` but **no `q_atm`**, which the implicit form needs. So the change is:
*add `q_atm` (and the conductances `gah`/`gaw`) to the argument list, drop `w_flux_ac`* — mirroring how
`enthalpy_atm` is already passed. After (a)+(b) all three twins share the same L-stable implicit-atm form
gated by `ustar × temp1/temp2`, the whole store set is L-stable at the target `dt_fast`, and one forced input
(`w_flux_ac`) disappears. **Sequenced in P3c** (the MVP diurnal cycle can blow up without it — §6).

### 3.6 Coupling contracts the driver must enforce

1. **Single evaporation authority — ground *and* leaf film (two instances of the same trap).**
   *(a) Ground.* Ground evaporation is computable *twice*: the hydrology kernel's `ground_evaporation()`
   (Swenson–Lawrence dry-surface-layer resistance + Philip pore-space RH `α_soil = exp(ψ₁·g/(R_wv·T_g))`,
   series with `r_aero`) and the energy kernel's `ground_surface_balance()` (plain conductance form
   `gbw·ρ·(q_sat(T_g) − can_shv)`, no DSL/α_soil). Pick the **DSL/α_soil hydrology form as authoritative**
   (physically richer), route its `flux%soil_evap [kg/m²/s]` into the CAS `ground_w_flux` twin, and make
   `ground_surface_balance` consume that same `LE_ground` (or drop its LE term) so the CAS never double-counts.
   *(b) Leaf film.* The **same double-computation exists for interception-film evaporation** and must be closed
   the same way: `veg_energy_balance`/`veg_surface_fluxes` *computes* the film-evap mass `w_flux = g_ev·ρ·grad`
   (credited to the CAS vapour twin, using the current `leaf_water` via σ_w), while `intercept_canopy_layer`
   *consumes* a supplied `e_canopy` that **debits** the film pool (`leaf_water += (q_intr − e_canopy)·dt`,
   clamped to `[0, w_max]`). If `e_canopy` lags by a step, the debit and the CAS credit are different numbers
   and the `leaf_water`↔CAS transfer does **not** close — breaking the machine-precision column-water budget
   §7 claims. Contract: the energy kernel's `w_flux` is the **film-evap authority**; feed that *same,
   current-step* value as `e_canopy` to the film-mass update (so run the film debit *after* `veg_energy_balance`,
   §3.3 step 7, not before at the interception step), and route any `[0, w_max]` clamp overflow as an explicit
   residual into the budget rather than silently dropping it. (The precip-throughfall part of interception
   stays early — it depends on rain, not on evap.)
2. **`r_aero ≡ 1/ggnet ≡ ground gbw`.** The hydrology soil-evap resistance and the energy ground-balance
   conductance are the *same* aerodynamic quantity in reciprocal units. The driver derives both from the one
   `ggnet` aerodynamics returns.
3. **Transpiration bookkeeping.** `leaf_flux_t%transpiration [mol H₂O/m²leaf/s]` from step (3) is the demand
   fed to hydraulics as `transp = E·leaf_area·M_h2o [kg/s/plant]`; the *supply-limited* flux from hydraulics
   (via `fs_open`/`psi_leaf` water stress) is what actually leaves through `veg_energy_balance`'s
   `g_series = gbw·gsw/(gbw+gsw)` and enters the CAS `coh_transp`. Root uptake `Σ nplant·root_uptake`,
   partitioned across layers by `root_frac`, is the hydrology `root_uptake(k)` sink. Set the hydrology
   `psi_wilt` *below* the plant's PLC/stomatal cutoff so `f_wilt ≈ 1` whenever the plant transpires —
   otherwise the soil sink double-counts the plant's own down-regulation and breaks the soil↔plant ΔW budget.
4. **Head↔MPa.** `flux%psi_soil [MPa] = grav_head · soil_psi_of_theta(...) [m]`, `grav_head = 9.804e-3 MPa/m`.
5. **LW-band coupling.** Radiation's `canopy_temp(coh)` is the leaf/wood temperature the fast loop is solving —
   a back-coupling. Refresh net LW from the evolving temps each `dt_fast` (analytically, ∝ T⁴) even when the
   full SW two-stream is run only once per `dt_fast`.

### 3.7 The meteorological forcing seam (a genuine prerequisite)

The fast loop needs a **sub-daily atmospheric boundary condition** that MEDS does not have — the carbon path
today uses a stub `cfg%gpp_ref`. This design adds a minimal met seam. **`cosz` is *not* a forcing field** — the
solar zenith cosine is fully determined by (calendar date, time-of-day, site latitude), so it is *computed*
from the time dimension, not read/stored. `met_forcing_t` therefore carries only the genuine atmospheric state:
```fortran
type :: met_forcing_t                    ! reference-height atmospheric state at time t
   real(wp) :: u_ref                      ! [m/s]   wind speed
   real(wp) :: t_air, theta_atm           ! [K]     air temp, potential temp
   real(wp) :: shv_atm                    ! [kg/kg] specific humidity
   real(wp) :: press                      ! [Pa]
   real(wp) :: co2_atm                    ! [µmol/mol]
   real(wp) :: rshort                      ! [W/m²]  total downwelling shortwave (beam+diffuse; split via cosz)
   real(wp) :: rlong                      ! [W/m²]  downwelling longwave
   real(wp) :: precip                     ! [kg/m²/s] rainfall
   real(wp) :: zref                       ! [m]     reference height
end type
```
Solar geometry is a small `pure` helper in `meds_time` (or a `meds_solar` sibling), keyed on the time cursor:
```fortran
pure real(wp) function solar_cosz(date, t_sec, latitude)   ! declination(day-of-year) + hour-angle(t_sec)
```
The driver computes `cosz = solar_cosz(now, t_sec, site%latitude)` each `dt_fast` and uses it to (a) build the
`rad_forcing_t` beam/diffuse split from `met%rshort` (e.g. Erbs/Spitters), and (b) drive the two-stream's
direct-beam geometry. So `cosz` is a derived quantity of the time dimension, consistent with the earlier note
that the stub SW/temperature cycle also comes from `solar_cosz(now, t_sec, latitude)` — one solar-geometry path
serves both the stub provider and the real met driver. (`site%latitude` is a polygon/site property, added to
state if absent — it is not a met input.)
Provider (behind an interface, so a real driver-file/netCDF reader drops in later): a **stub analytic diurnal
cycle** first (sinusoidal SW/temperature from the solar geometry, constant wind/CO₂, prescribed humidity) —
enough to exercise the full coupled loop and reproduce a diurnal GPP/energy cycle. **`solar_cosz` must be keyed
on the driver-local sub-daily cursor** `t_sec` (the hour-of-day accumulator, §3.1) — *not* on `now` alone:
since the calendar `now` is frozen at day resolution within the fast loop, a `now`-only `cosz` is constant
across all `n_fast_per_slow` substeps and there is no diurnal signal (the P3c
diurnal-cycle test would be vacuous). Make the `t_sec → solar-geometry → SW/temperature` path explicit. `enthalpy_atm` for the CAS is
`cas_enthalpy_of_temp(t_air, shv_atm)`. This seam is scoped here because the column driver cannot run without
it; the real met I/O is a follow-on (P3d).

### 3.8 Closing the loop — fast → slow handoff

The whole point of coupling the fast loop is to **replace the stub GPP** feeding demography. The driver
accumulates per-cohort daily carbon over the day's `n_fast_per_slow` steps (GPP = Σ `a_gross·leaf_area·dt_fast`,
leaf respiration, etc.), and at the day roll hands the daily totals to `vegetation_dynamics` in place of
`cfg%gpp_ref · leaf_area · dt` inside `carbon_growth`. The demography data-seam `update_demography` is
untouched — only the *source* of the carbon changes (stub → integrated fast-loop GPP−R). This is exactly the
last coupling the CLAUDE.md "Reserved follow-ups" names ("coupling the leaf-physiology module into the
demographic growth … needs canopy RT, leaf energy balance, a meteorological forcing source, and plant
hydraulics") — all four now exist, and this driver joins them.

### 3.9 Reserved: the global adaptive-RK path

Documented, not recommended for P3. ED2's production integrator packs `{can_enthalpy, can_shv, can_co2,
soil_energy(k), soil_water(k), sfcwater×3, virtual×3, per-cohort leaf/wood energy+water+internal-water+2 psi}`
into one `real(dp)` state vector and integrates with Cash–Karp RK45 (embedded 4th/5th-order error), warm-started
step `htry`, sanity-bounds rejection, and per-stage aerodynamics re-evaluation. If MEDS ever needs that
accuracy, the clean seam (from the ED2 integrator study) is: a scheme-generic `meds_integrator` harness (the
adaptive controller `h_new = clamp(safety·h·errmax^p, 0.1h, min(5h, dt_fast))`, the `yscal = |y|+|dydx·h|`
error policy, the Butcher tableau as data, the stage micro-kernel over a *flat* `real(wp) state(:)`) plus a
problem-specific `assemble_rhs` that unpacks state → diagnostics → aerodynamics → fluxes → `d/dt`. The stiff
CAS scalars would be handled IMEX (analytic exponential relaxation, `C(t+h) = C∞ + (C0−C∞)·exp(−b·h)`,
`b = κ·capci`), *not* ED2's leaky BDF2 hybrid (whose ×100 budget-tolerance relaxation is a red flag). This
requires rewriting the shipped single-step kernels into pure-derivative form and is a separate project.

### 3.10 Is the operator-split scheme faster than adaptive RK4?

The comparison is "operator-split-of-implicit-steps" (this design) vs ED2's monolithic **adaptive explicit**
Cash–Karp RK45 — a different *decomposition*, not the same system with a different stepper. Expected verdict:
**usually faster in the regimes that matter, decisively so on GPU — but it is a trade, and speed is the third
reason to prefer it, behind stability and portability.** Reasoning:

- **Stability, not stability-limited stepping.** RK45 is *explicit* → step size pinned by the stiffest mode.
  Here the stiff modes are severe: the CAS relaxation timescale is `wcapcan/(ρ·ustar·temp1)` (often seconds)
  and wet-soil diffusion. An explicit method must keep `h` below those or diverge, so it spends many small
  adaptive substeps (plus rejected retries). The split-implicit steps are L-stable, so `dt_fast` is limited by
  **accuracy**, not stability — far fewer, larger steps on exactly the parts that hurt RK4. (ED2 itself ships
  the Euler+BDF2 *hybrid* precisely because adaptive RK4 is the slow, stiffness-pinned option; MEDS takes the
  "implicit for speed" lesson but keeps conservation — each twin closes its budget by construction, avoiding
  the hybrid's leak.)
- **The implicit solves are cheap.** Every system is tiny and structured — tridiagonal `thomas_solve` (O(n),
  no pivot), 2×2 analytic matrix-exp, scalar-algebraic CAS updates, 2-evaluation linearized-BE leaf energy —
  so "implicit" carries no Newton/dense-factorization penalty. Meanwhile RK45 pays 6 stages × per-stage
  aerodynamics re-evaluation per substep. The expensive nonlinear physiology (Farquhar Ci solve) is frozen
  once per `dt_fast` in *both* schemes, so it does not separate them.
- **GPU throughput — the biggest win.** Adaptive RK45 across patches is **warp-divergent**: per-patch substep
  counts and step rejections differ, so the warp waits on the slowest lane. The split scheme in
  `*_SUBSTEP_FIXED` mode is fixed-step, branch-light, warp-uniform, root-finder-free (CLM fixed-iteration MO).
  For MEDS's GPU target this is a structural, not marginal, speedup.

Caveats, stated honestly: operator splitting is **1st-order** (Lie–Trotter, error ∝ `dt_fast`) vs RK45's
5th-order, so tight accuracy would force a smaller `dt_fast` and erode the gain — the split wins in practice
only because land-surface accuracy targets are loose (few-% fluxes), which is why CLM/Noah-MP/JULES all run
first-order splitting at 15–30 min steps. Per-kernel internal step-doubling still sub-cycles when a store is
accuracy-demanding, but (being L-stable) far less than an explicit method sub-cycles for *stability*. Net: the
speed advantage is expected but `dt_fast`-dependent — the **P3c diurnal benchmark should measure split vs a
reference RK4 at matched accuracy**, not assume it. The primary reasons for the choice remain unconditional
stability + GPU portability + behaviour-preservation (the kernels were already single-step implicit).

---

# PART III — Numerics consolidation review (`meds_numerics.f90`)

## 4. Should the solvers be consolidated? — Yes, one shared module.

A full grep sweep of `src/` found the solvers below. The DAG is
`shared ← {allometry, plant} ← state ← demography ← aux ← main`, with **biophysics linking `shared` only and
plant linking `shared` only** — so the *only* library reachable by both the biophysics and plant kernels is
`src/shared/`. Any solver shared across the biophysics/plant wall must therefore live in a new
**`src/shared/meds_numerics.f90`** at the root of the DAG (deps: `meds_kinds` only).

### 4.1 Inventory

| # | File / routine | Method | Shape | pure / GPU | Note |
|---|---|---|---|---|---|
| 1 | `biophysics/meds_soil_solver.f90` `thomas_solve` | Tridiagonal Thomas, no pivot | fixed `n_soil_layer_max` | **pure**, issue-#7-safe | canonical, 2 callers (#2, #3) |
| 2 | `biophysics/meds_column_energy` `soil_heat_be_step` | BE heat conduction | tridiag → #1 | pure | uses #1 |
| 3 | `biophysics/meds_column_hydrology` `soil_be_single_step` | BE + Celia Picard | tridiag → #1 | sub | uses #1 |
| 4 | `biophysics/meds_column_hydrology` `soil_water_advance` | adaptive step-doubling | wrapper | sub | controller ≈ #6 |
| 5 | `biophysics/meds_column_hydrology` `compute_psi_e` | 5-pt midpoint quadrature | 1-D quad | pure | domain-specific |
| 6 | `plant/meds_plant_hydraulics` `solve_plant_water`/`expm_step` | 2×2 matrix-exp + step-doubling | 2×2 analytic | sub | controller ≈ #4 |
| 7 | `plant/meds_plant_hydraulics` `phi_inverse` | bisection (60 it), monotone | scalar root | pure | ≈ #9 |
| 8 | `plant/meds_plant_hydraulics` `flux_potential` | 7-pt Gauss–Legendre | 1-D quad | pure | GL nodes reusable |
| 9 | `plant/meds_leaf_gas_exchange` `solve_leaf_gas_exchange` | bracketed bisection on Ci (100 it) | scalar root | sub | ≈ #7 |
| 10 | `plant/meds_leaf_gas_exchange` `smaller_root` | analytic quadratic smaller-root | scalar | **elem pure** | ≈ PV-curve quadratic |
| 11 | `biophysics/meds_twostream_band` `solve_band`/`layer_rt` | adding/matrix-operator block-tridiag | block-tridiag | sub | specialized RT |
| 12 | `shared/meds_thermo` `uext_to_temp` | closed-form phase inverter | scalar | elem pure | already shared |

### 4.2 Consolidation decisions

| Solver | Action | Target | Rationale |
|---|---|---|---|
| **`thomas_solve` (#1)** | **MOVE** to `shared/meds_numerics.f90`; delete `meds_soil_solver.f90` | shared | Used by biophysics energy (#2) + hydrology (#3); the CO₂ K-theory multi-layer-CAS forward-design and the future 3-node hydraulics BE (`HYDRO_SOLVER_BE`, declared but unimplemented in `plant`) will *also* need a tridiag — and plant **cannot link biophysics**. Only `shared` satisfies both walls. Keep the fixed-size signature (GPU-eligible). **Highest-value consolidation.** |
| **Bracketed bisection (#7, #9)** | **EXTRACT** `bracketed_root(...)` | shared | `phi_inverse` and the leaf-Ci solve are the same monotone-bracket bisection; a future hydrology Neumann↔Dirichlet switch wants it too. Offer a callback form (host) *and* keep the loop inlineable where device-resident (procedure-pointer args are awkward on GPU). |
| **Quadratic smaller-root (#10)** | **PROMOTE** `elemental pure quadratic_smaller_root(theta,a,b)` | shared | Co-limitation (leaf) + PV-curve inverse reuse the same guarded quadratic; tiny, elemental, GPU-clean. |
| **Adaptive step-doubling controller (#4, #6)** | **EXTRACT** the scalar `adaptive_step_update(err, h, accept, safety, fmin, fmax, h_new)` — constants as **arguments**, not hard-coded | shared | Same *shape* (`h·min(fmax, safety·err^−0.5)`) but **different constants and clamp policy** — NOT byte-for-byte: hydrology uses `safety=0.9, fmin=0.25, fmax=4` with no fmin clamp on the accept branch (`meds_column_hydrology.f90`), hydraulics uses `fmin=0.2` with an extra `max(fmin,·)` on accept (`meds_plant_hydraulics.f90`). The extracted helper must take `safety/fmin/fmax` + the clamp policy as args so each caller keeps its exact shipped behaviour (a hard-coded `fmin=0.25` would silently change the tested hydraulics kernel — violating the P3b bit-for-bit rule). Extract only the *scalar update*; the loop bodies stay local. |
| **Gauss–Legendre nodes (#8)** | **PROMOTE** `gl_x(7), gl_w(7)` parameter block | shared | any future quadrature reuses. |
| #2, #3, #4-body, #5, #6-body, #11 | **STAY** local | biophysics / plant | domain-specific integrands / assembly; they just `use meds_numerics` for the shared primitives. |
| #12 `uext_to_temp` | **STAY** in `meds_thermo` | shared | already shared, not iterative. |

`meds_numerics` slots next to `meds_kinds`/`meds_constants` at the DAG root. The two `use meds_soil_solver,
only: thomas_solve` sites (`meds_column_energy`, `meds_column_hydrology`) become
`use meds_numerics, only: thomas_solve`. Preserve the subroutine-with-`intent(out)`-array signature (never an
array-valued function result into a call — issue #7).

### 4.3 What this consolidation is *not*

It does **not** merge the process integrators (each kernel keeps its own BE/Picard/matrix-exp step — those are
physics, not shared numerics), and it does **not** introduce a global integrator harness (that is the reserved
§3.9 path). It is strictly the extraction of the reusable *primitives* (tridiag solve, bracket root, guarded
quadratic, GL nodes, step-doubling update) into one shared module — the minimum that removes the duplication
and unblocks the plant-side tridiag need.

---

# PART IV — Conservation checking (`meds_budget_check.f90`)

## 4A. Why a shared budget checker now

Coupling makes conservation a *cross-module* property. Today every fast kernel closes and checks its **own**
budget in isolation, with the *same logic re-implemented three times*: `soil_energy_flux` computes
`flux%energy_resid` and `error stop 'energy budget did not close'`; `column_co2_step` computes `budget%resid`
and guards `abs(resid) > rtol·scale + atol`; `column_hydrology_flux` computes `flux%mass_resid` and guards on
`atol`. Each is a per-kernel, per-store check with its own tolerance form. Once the driver sums cohort + ground
+ atmosphere fluxes into the CAS twins and threads soil↔plant water, the meaningful invariant becomes the
**column total**: does the change in *all* energy stores over `dt_fast` equal the net boundary energy flux (and
likewise for water and CO₂)? No single kernel can see that — only the driver can. So the coupling needs one
shared facility that both (a) unifies the per-kernel closure checks and (b) assembles the column-total budget.
This is the natural, real-code replacement for the removed empty `src/utils/` placeholder.

## 4B. `src/shared/meds_budget_check.f90`

Lives in `shared` (root of the DAG, deps `meds_kinds`/`meds_constants` only) so every process module *and* the
driver use the identical closure test. Pure/GPU-safe compute; the `error stop` path is a thin non-pure wrapper
used only in Debug, mirroring the existing `debug_error` discipline.

```fortran
module meds_budget_check
   ! A budget accumulator + a uniform closure test for energy [J/m2], water [kg/m2], CO2 [umol/m2 or mol/m2].
   type :: budget_t
      real(wp) :: store0   = 0.0_wp   ! store at t^n              (extensive, per m2 ground)
      real(wp) :: store1   = 0.0_wp   ! store at t^{n+1}
      real(wp) :: influx   = 0.0_wp   ! Σ boundary IN  over the accumulation window [store-unit/s]
      real(wp) :: outflux  = 0.0_wp   ! Σ boundary OUT over the window
      real(wp) :: resid    = 0.0_wp   ! (store1-store0) - dt*(influx-outflux)   ← the imbalance
      real(wp) :: worst    = 0.0_wp   ! running max |resid| over the window (for reporting)
      integer(ik) :: n_check = 0_ik, n_fail = 0_ik
   end type

   ! (1) The single closure predicate — the ONE place the mixed rtol/atol form lives.
   pure logical function closure_ok(resid, scale, rtol, atol)          ! |resid| <= rtol*scale + atol

   ! (2) Compute + record one store's residual (pure; fills b%resid, b%worst, counters).
   pure subroutine budget_residual(b, store0, store1, influx, outflux, dt)

   ! (3) Assert closure — status-code form (pure-friendly) OR Debug error-stop wrapper.
   pure subroutine budget_assert(b, scale, rtol, atol, ok)             ! ok=.false. if breached (no stop)
   subroutine budget_check_stop(b, scale, rtol, atol, label, debug)    ! error stop label in Debug; else no-op

   ! (4) Column-total assembly the DRIVER calls once per dt_fast per patch: fold every store's Δ and every
   !     boundary flux into one energy / water / CO2 budget_t, then assert. Sums are explicit do-loops
   !     (nvfortran array-temp trap, issue #7), not sum(...).
   subroutine column_energy_budget(...)   ! Σ Δ(soil_energy, leaf/wood store_energy, CAS enthalpy·wcapcan)
                                          !   vs (Rn_net − H_atm − LE_atm − drainage_enthalpy − ...)
   subroutine column_water_budget(...)    ! Σ Δ(theta·ρ_w·dz, w_surface, leaf_water, CAS shv·wcapcan)
                                          !   vs (precip − ET − drainage − runoff)
   subroutine column_co2_budget(...)      ! Δ(CAS co2·ccapcan) vs (Rh + Ra − GPP − loss2atm)
end module
```

**How it plugs in.** (i) The existing per-kernel guards (`soil_energy_flux`, `column_co2_step`,
`column_hydrology_flux`) migrate to call `closure_ok`/`budget_check_stop` instead of their bespoke inline
tests — one tolerance policy, one message path (a pure-refactor with unchanged behaviour, P3b). (ii) The driver
calls the three `column_*_budget` routines once per `dt_fast` per patch, accumulating `worst`/`n_fail` over the
slow step, and reports (or `error stop`s in Debug) if the **column** total drifts — the invariant that the
soil-evap and leaf-film single-authority contracts (§3.6) and the implicit-vapour fix (§3.5) exist to satisfy.
(iii) The P3c and P3e tests assert `worst → machine precision`; this module is what makes "every fast step
closes energy/water/CO₂ to machine precision" (§7) a *checked* invariant rather than a hope. It also subsumes
the removed `src/utils/` placeholder — the one general-purpose shared facility the coupling actually needs.

---

## 5. Files touched + CMake wiring

| Action | Path | Notes |
|---|---|---|
| NEW | `src/biophysics/meds_canopy_aerodynamics.f90` | stateless kernel (Part I) |
| NEW types | `src/biophysics/meds_biophysics_types.f90` | `aero_cfg_t`, `aero_env_t`, `aero_geom_t`, `aero_cohort_t`, `aero_out_t` + `AERO_*` selectors |
| NEW | `src/driver/meds_column_dynamics.f90` | fast integrator (Part II); globbed into `meds_aux` |
| NEW | `src/shared/meds_numerics.f90` | consolidated primitives (Part III) |
| NEW | `src/shared/meds_budget_check.f90` | conservation checker (Part IV); the real replacement for the removed `src/utils/` placeholder |
| DELETE | `src/biophysics/meds_soil_solver.f90` | `thomas_solve` → `meds_numerics`; repoint 2 `use` sites |
| **DELETE** | `src/utils/` | empty placeholder (README only, not in CMake); superseded by `meds_budget_check` in `shared`. Update the stale CLAUDE.md:186 line (which still calls `src/utils/` a placeholder) |
| NEW type | `src/driver/` or a met module | `met_forcing_t` (no `cosz` — derived) + stub diurnal provider (§3.7) |
| NEW fn | `src/shared/meds_time.f90` (or `meds_solar`) | `pure solar_cosz(date, t_sec, latitude)` — cosz from the time dimension (§3.7) |
| **EDIT** | `src/shared/meds_config.f90` | **retire `ts_mode`/`TS_*`/`dt_years`-from-enum; add `dt_slow` (default `"1d"`)**; derive `dt_years = dt_slow/yr_sec`, `n_fast_per_slow`; add `fast_biophysics_on`, `dt_fast`, `integration_scheme`; validation |
| EDIT | `src/io/meds_config_io.f90` | parse `[run] dt_slow` duration string (was `timestep` enum) + the `[fast]` block |
| EDIT | `src/driver/meds_main.f90` | calendar walk advances by `dt_slow` (was the `ts_mode` case-switch) |
| EDIT | `meds_config_main.toml` | `[run] timestep = "daily"` → `[run] dt_slow = "1d"`; add `[fast]` block (`fast_biophysics_on`, `dt_fast`, `integration_scheme`) |
| EDIT | `src/driver/meds_stepper.f90` | call `column_dynamics` before the slow loop |
| EDIT | `src/state/meds_demography_types.f90` | per-cohort `leaf_store_energy`, `wood_store_energy`, `leaf_water`, `psi(:)`, resolvable flags — through the lockstep reorder |
| EDIT | `src/demography/meds_demography_fusefiss.f90` | per-patch `cas`/`soil_w`/`soil_e` through `sort_patches` + `patch_compact`; conserving fusion averages for the new stores |
| EDIT | `src/shared/meds_pft_params.f90` | per-PFT `leaf_width`, `branch_diam` (if absent) |
| **EDIT** | `CMakeLists.txt` | **`meds_aux` must additionally link `meds_biophysics` + `meds_biogeochemistry`** (currently it links `meds_demography + meds_config_io + meds_plant` only) — the column driver calls all of them |

CMake auto-resolves `.mod` deps; the only manual edit is the `meds_aux` link set. `meds_numerics` compiles
standalone (deps: `meds_kinds`).

## 6. Phasing (independent, testable PRs)

- **P3a — aerodynamics kernel, standalone.** `meds_canopy_aerodynamics` + types + PFT traits. Unit tests: neutral
  log-law recovery (`ustar = vonk·u/ln(z/z0)`), stable/unstable stability-function limits vs analytic, leaf
  `gb` vs published Nu/Re/Gr values, monotonic in-canopy wind decay, `ggnet → ggbare` as `opencan → 1`. No
  driver, no state — a leaf kernel like RT. Ship a Python C-API demo (mirrors `meds.leaf`) if useful.
- **P3b — shared-facility consolidation.** `meds_numerics` (delete `meds_soil_solver`; repoint) **and**
  `meds_budget_check` (delete the empty `src/utils/`; migrate the per-kernel `soil_energy_flux`/
  `column_co2_step`/`column_hydrology_flux` closure guards to the shared `closure_ok`/`budget_check_stop`).
  Pure refactor, green suite unchanged (bit-for-bit).
- **P3c — column-dynamics MVP + config refactor.** The `dt_slow` config change (retire `ts_mode`; `[run]
  dt_slow = "1d"`, §3.1) + `solar_cosz` helper + `met_forcing_t`. State additions + reorder threading +
  `column_dynamics` (`SPLIT_SEQUENTIAL`) + stub met provider (SW/temperature from `solar_cosz(now, t_sec,
  latitude)`, §3.7). **Includes the correctness-critical fixes** — *not* deferrable, because the MVP diurnal
  budget test depends on them: (i) the profile-factored atm↔CAS conductance `gah/gaw/gac` replacing the
  placeholder `ρ·ustar` (§3.5a); (ii) the implicit vapour twin so the store set is L-stable at `dt_fast` =
  300–900 s (§3.5b — without it the diurnal water/energy cycle can blow up); (iii) the single evaporation
  authority for ground *and* leaf film so the column water budget closes (§3.6 contract 1). Single patch,
  single cohort, prescribed diurnal met → diurnal GPP/energy/water cycle; assert the per-kernel residuals
  *and* the `meds_budget_check` **column-total** energy/water/CO₂ budgets (Part IV) close to machine precision
  each step.
- **P3d — real met forcing.** netCDF/CSV driver reader behind the `met_forcing_t` interface.
- **P3e — Picard/inner-leaf coupling.** `SCHEME_PICARD_COUPLED` (the CLM-style leaf-temp/stomata inner fixed
  point, §3.4) — a pure accuracy/tight-coupling upgrade on top of the already-stable, already-conservative P3c.
- **P3f — fast→slow handoff.** Daily GPP/NPP accumulation replaces the stub `gpp_ref` in `carbon_growth`;
  end-to-end multi-year run with mechanistic carbon growth.

## 7. Testing strategy

- **Per-kernel unit tests** (CTest, ifx + nvfortran multicore — a green ifx run is *not* sufficient, per the
  issue-#7 portability trap): aerodynamics limits (above), `meds_numerics` primitives (tridiag vs dense LU
  reference, bisection convergence, GL quadrature exactness on polynomials, step-update monotonicity),
  `meds_budget_check` (`closure_ok` boundary cases; a hand-built store with a known injected imbalance is
  flagged, a closed one passes; `solar_cosz` vs an ephemeris at solstice/equinox).
- **Conservation** — the MEDS discipline, now a *checked* invariant via `meds_budget_check` (Part IV): every
  fast step closes the per-kernel residuals (`resid`, `budget%resid`, `mass_resid`) *and* the **column-total**
  energy/water/CO₂ budgets to machine precision; the driver accumulates `worst`/`n_fail` over the slow step and
  `error stop`s in Debug on drift. The column-total budget across all stores over a full day is the acceptance
  gate for the soil-evap/leaf-film single-authority (§3.6) and implicit-vapour (§3.5) fixes.
- **Physical sanity** — a bare-soil patch reproduces a Monin–Obukhov surface-flux diurnal cycle; a
  single-cohort patch reproduces a light-response GPP curve and a realistic leaf-to-air temperature
  difference; nighttime `ustar → ustmin`, GPP → 0, CAS CO₂ rises under respiration.
- **ED2 comparison** — drive `meds_canopy_aerodynamics` with an ED2 patch snapshot and compare `ustar`,
  `gbh/gbw`, `ggnet`, `veg_wind` against `canopy_turbulence8` outputs within tolerance (the EDTS analogue).

## 8. Open questions / risks

1. **Split accuracy at large `dt_fast`.** First-order Lie–Trotter splitting error scales with `dt_fast`;
   default 300–900 s should be safe (ED2's DTLSM is comparable) but the diurnal-cycle test must confirm the
   leaf-temp/CAS coupling doesn't oscillate. Mitigation: `SCHEME_PICARD_COUPLED` or a smaller `dt_fast`.
2. **Fusion of prognostic thermal/hydraulic stores.** Conserving averages for `psi` (intensive) on cohort
   fusion need care; energy/water (extensive) follow the AGB rule. Verify the 1% budget assertion holds
   across a fuse of two cohorts at different temperatures/potentials.
3. **The soil-evap double-count** (§3.6) is a real reconciliation, not just wiring — must be resolved before
   P3e or the CAS water budget will not close.
4. **`can_depth` as prognostic vs forced.** ED2 lets it float (constant-ρ,P); the MVP fixes it per step from
   canopy height (floored at 5 m). Revisit if canopy-air mass conservation matters.
5. **GPU path.** The fast loop is patch-parallel; keep aerodynamics/leaf kernels `pure` over SoA slices, force
   all `*_SUBSTEP_FIXED` selectors for warp uniformity, and keep the met provider and reductions host-visible.
   The `SPLIT_SEQUENTIAL` sweep (no global lockstep) is the GPU-friendly default; validate on the RTX 3050 Ti
   alongside ifx.

---

### Appendix A — proposed `aero_out_t` (the driver-facing contract)

```fortran
type :: aero_out_t
   ! --- per-patch scalars ---
   real(wp) :: ustar   = 0.0_wp     ! [m/s]   friction velocity (momentum)
   real(wp) :: temp1   = 0.0_wp     ! [-]     heat scalar profile factor vonk/D_h  → gah = ρ·ustar·temp1
   real(wp) :: temp2   = 0.0_wp     ! [-]     vapour/CO₂ profile factor vonk/D_h(z0q) → gaw, gac
   real(wp) :: tstar   = 0.0_wp     ! [K]     temperature scale (= temp1·dth)
   real(wp) :: qstar   = 0.0_wp     ! [kg/kg] humidity scale (= temp2·dqh)  → w_flux_ac = ρ·ustar·qstar
   real(wp) :: cstar   = 0.0_wp     ! [µmol/mol] CO₂ scale (= temp2·dco2)
   real(wp) :: zeta    = 0.0_wp     ! [-]     stability z/L (diagnostic)
   real(wp) :: rib     = 0.0_wp     ! [-]     bulk Richardson (diagnostic)
   real(wp) :: obu     = 0.0_wp     ! [m]     Obukhov length (diagnostic)
   real(wp) :: ggbare  = 0.0_wp     ! [m/s]   bare-ground conductance
   real(wp) :: ggnet   = 0.0_wp     ! [m/s]   net ground↔CAS conductance → r_aero = 1/ggnet; ground gbh=gbw
   real(wp) :: can_depth = 0.0_wp   ! [m]     CAS depth (from canopy height)
   ! --- per-cohort arrays (caller-owned; written in place) ---
   real(wp), allocatable :: wind(:)     ! [m/s]     in-canopy wind at crown mid-height
   real(wp), allocatable :: leaf_gbh(:) ! [m/s]     leaf boundary-layer heat conductance
   real(wp), allocatable :: leaf_gbw(:) ! [m/s]     leaf boundary-layer vapour conductance
   real(wp), allocatable :: wood_gbh(:) ! [m/s]     wood boundary-layer heat conductance
   real(wp), allocatable :: wood_gbw(:) ! [m/s]     wood boundary-layer vapour conductance
end type
```

### Appendix B — the reference call-order tension (why operator-split, one page)

ED2 wraps the *entire* fast column in an adaptive Cash–Karp RK45 over a monolithic `real(dp)` state vector,
re-evaluating aerodynamics every stage; the CAS scalars are prognostic ODEs with air-mass capacities. MEDS's
shipped kernels are instead each a *single implicit step* — `veg_energy_balance` (linearized BE),
`canopy_air_update`'s enthalpy branch / `canopy_air_co2_update` (implicit-in-atm algebraic), `soil_energy_flux`/
`column_hydrology_flux` (BE-Thomas with private step-doubling), hydraulics (2×2 matrix-exp) — all L-stable
*except* the shipped `can_shv` twin (explicit-in-atm), which §3.5 makes implicit (required in P3c). That design
*is* an operator-split-per-store integrator; the column driver is its sequencer, not a global ODE solver.
CLM confirms the middle path: keep the prognostic CAS (CLM's diagnostic junction is its degenerate limit) and,
if needed, borrow only CLM's *inner* leaf-temperature/stomata fixed point for the stiff leaf sub-problem. The
result — `SPLIT_SEQUENTIAL` default, `PICARD_COUPLED` opt-in, global-RK reserved — matches what MEDS already
built, is unconditionally stable at large `dt_fast` once §3.5 lands, and keeps the whole fast loop GPU-portable.
