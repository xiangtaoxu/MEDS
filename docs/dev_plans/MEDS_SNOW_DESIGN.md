# MEDS Surface-Water & Snow — Module Family Design

> **STATUS: P0 IMPLEMENTED** (2026-07-13, branch `feature/snow-model`, commits 12701c5 A / 1bebc93 B /
> ca75e32 C / 05a3ad4 D; ifx + nvfortran multicore 34/34, snow-off bit-identical). The stateless kernels
> (`meds_snow_energy` + `meds_snow_mass`), the per-patch `snow_column_t` state + lockstep, and the split-path
> fast-loop coupling all landed and **close the whole-column mass + energy budgets to machine precision**
> (|E|~5e-7 J/m², |W|~3e-13 kg/m²) through accumulation, sublimation, melt→infiltration, and the snow albedo
> ramp. **The `snowfac` continuity ramp (§4f/§4g/§6) IS wired** (commit 0a6b88b): a genuine sub-column area
> weight by the Niu-Yang07 fraction scales the snow store's boundary exchange, blends the surface
> (albedo/sensible/latent/base-conduction) as `(1−snowfac)·bare + snow`, and `snowfac→0` makes a thin pack a
> near no-op — so partial cover gives a partial (continuous) albedo and the old passive-thin `snow_stab_thresh`
> cliff is gone; snowfac=0 reduces to the bare-soil skin exactly. **Remaining MVP deviations (follow-ups):**
> the forced snow↔soil equilibrium for a genuinely melting thin pack (§6) is superseded by the snowfac no-op
> but a thin pack still can't fully self-melt from its own surface balance; ARK stays snow-free (split-path
> only). The dramatic clear-day 0 °C soil cap (§9 test 5) needs a deeply-frozen-soil + strong-sun
> + persistent-pack alignment the recycled forcing didn't cleanly provide, but the component physics
> (conservation, albedo reduction, `g_base` insulation damping the soil diurnal swing) is validated.



A **stateless** biophysics module family for MEDS (`src/biophysics/`) owning the **temporary surface-water / snow store**: the stacked mass+energy reservoir that sits *between* the canopy-air-space (CAS)/atmosphere and `soil_energy(1)`. It closes the last open surface seam — the winter ground skin — by inserting a snowpack that (a) caps the soil surface near `t_3ple` on clear cold days, (b) raises ground albedo, and (c) thermally insulates the soil from the atmosphere. It owns snow **accumulation** from frozen precip, the snow-surface **energy balance** (net SW/LW, sensible-to-CAS, sublimation/evaporation, base conduction), **freeze/thaw + melt** as an internal-energy read-off, **percolation** of meltwater into the soil, and the snow-base **conduction** that becomes the soil top BC.

The physics reference is **ED2**'s temporary-surface-water (`sfcwater`) model (`memory/ed_state_vars.F90:1242`, `dynamics/rk4_derivs.f90`, `dynamics/rk4_misc.f90:adjust_sfcw_properties`), **modernized** per MEDS conventions toward **CLM5/ClimaLand internal-energy state design**: the prognostic state is **internal energy + water-equivalent mass**, never a prognostic temperature — snow temperature and liquid fraction are a **read-off of the one shared inverter** `uext_to_temp` (`meds_thermo.f90:55`), exactly the machinery that already gives the soil column its freeze/thaw plateau under `ENERGY_PHASE_ON`.

**Per the user directive, this family is STATELESS-FIRST and OPERATOR-SPLIT.** The snowpack is a new reservoir advanced by a backward-Euler substep inside the split/Picard driver; it is **kept OUT of the ARK implicit (ESDIRK) set** in P0 (§6), mirroring how soil water and prognostic wood are already operator-split out of the ARK. The coupled surface fixed point remains a P3-orchestration concern of `column_fast_step`.

Units are SI throughout: mass/SWE `[kg/m²]`, internal energy `[J/m²]` (extensive, per unit ground area), depth `[m]`, temperature `[K]`, fluxes `[W/m²]` or `[kg/m²/s]` positive **upward**. Kinds `wp/ik`, `_wp` literals, `implicit none`, `pure`/`elemental` constitutive + inversion kernels, `error stop` (not `fatal_error`), ≤132 columns, nvfortran-safe.

**Deferred siblings (out of scope here).** Multi-layer snow with compaction/densification and aging albedo is P1. Fractional snow cover and canopy snow interception/unload are P2. Blowing/drifting snow, snow shortwave masking of vegetation beyond ground albedo, and glaciers/firn are out entirely (§8). Promotion of the snowpack into the ARK implicit set is deferred and gated (§6).

---

## 1. Scope, target processes & reachability [facts]

### 1.1 The store, the coordinate convention, and the enabling insight

**The governing decision (Bundle E, ED2, ClimaLand): the snowpack is a stacked surface store whose prognostic state is water-equivalent mass + internal energy, sharing the existing inverter, the BE-Thomas conduction, and the surface energy balance.** This is the `surface_layer(:)` array the energy doc already named as the deferred snow home (`MEDS_ENERGY_BALANCE_DESIGN.md:283`): the analogue of ED2's `sfcwater_energy(nzs)`/`sfcwater_mass(nzs)`. Snow, ponded water, and a future litter horizon are all **stacked surface layers** between the CAS/atmosphere and `soil_energy(1)`. MVP is a **single bulk layer**; the multi-layer stack is deferred to P1.

**Coordinate & sign convention (inherited).** Negative-z geometry (`MEDS_COLUMN_HYDROLOGY_DESIGN.md:17`): `z` is elevation `[m]`, positive up, `z=0` at the ground surface; soil nodes are negative, canopy positive. **Snow layers stack above `z=0`** (positive z), between the CAS and `soil_energy(1)` — which is precisely why the energy doc calls them a store "*between* the atmosphere/CAS and `soil_energy(1)`". Fluxes are **positive upward** (matching `ground_surface_balance` and `soil_energy_flux`'s `hf(0) = -g_top`).

**The enabling insight, stated once: snow is mostly a NEW STORE + wiring, not new numerics.** Everything the melt/refreeze physics needs already exists:

- The **internal-energy↔temperature inverter with liquid fraction** — `uext_to_temp(uext, wmass, dry_hcap, temp, fliq)` (`meds_thermo.f90:55-68`) — already pins `temp = t_3ple` on the plateau while `fliq` absorbs `wmass·latent_heat_fusion`. It is documented to work for **any per-area or per-volume store** (soil J/m³, leaf/wood J/m²). A snow layer reuses it verbatim with `dry_hcap ≈ 0` (snow is ~pure water/ice), `wmass = SWE`. Forward map `temp_to_uext` (`meds_thermo.f90:71`) seeds a fresh layer at first snowfall.
- The **surface energy seam** — `ground_surface_balance` (`meds_column_energy.f90:307`) plus its two live inlined copies (`meds_column_dynamics.f90:504-506`, `meds_column_derivs.f90:256-258`) — is the exact Rn−H−LE algebra the snow surface balance replaces.
- The **phase split** — `precip_phase(precip_total, tair_k, rainf, snowf)` (`meds_forcing_kernels.f90:291`) — already returns mass-conserving `snowf [kg/m²/s]` (ED2 Jin-1999 band, `PHASE_BAND_K=1.0`). Today `snowf` is computed and **dropped on the floor** (`meds_fast_loop.f90:482` threads only `met%rainf`); this is the primary dangling seam.
- The **soil thermal Thomas solver** — `soil_energy_flux` / `soil_heat_be_step` via `meds_soil_solver` — is reused unchanged for snow-column conduction (1-D heat on stacked layers is structurally identical to the soil BE solve).
- The **whole-column budget ledgers** — `budget_accumulate` / `track_resid` for `whole_energy`/`whole_water` (`meds_column_dynamics.f90:597-615`) — already have the exact insertion pattern used for prognostic wood/leaf (`wood_store0/1`, `leaf_store0/1`).

### 1.2 Target processes and MVP vs deferred

| Process | Governing choice | MVP (P0) | Deferred |
|---|---|---|---|
| Accumulation from snowfall | frozen-precip mass+ice-enthalpy add to top layer | ✔ | — |
| Snow-surface energy balance | net SW(snow albedo) + net LW − H_snow−LE_snow − G_base | ✔ | multi-layer top-skin (P1) |
| Freeze/thaw + melt | `uext_to_temp` plateau read-off, `dry_hcap≈0` | ✔ | — |
| Sublimation / evaporation | vapor twin: sublimation when frozen, evaporation when wet | ✔ | — |
| Meltwater → soil | route melt into `q_avail`, existing infiltration/ponding | ✔ | — |
| Snow albedo | fresh-vs-aged VIS/NIR on `fast_context_t`, into `ground_optics` | ✔ (fresh + simple age) | full spectral aging (P1) |
| Thermal insulation | snow-base conduction replaces bare-soil top BC | ✔ | — |
| Snow-cover / burial fraction | `snowfac` from SWE/depth → aero roughness | ✔ (Niu-Yang07) | — |
| Multi-layer stack | ED2 `nzs` layers + `thick` table redistribution | — | ✔ P1 |
| Compaction / densification | density evolves; `th_cond` density-dependent | — | ✔ P1 |
| Aging albedo | fresh→aged decay curve | — | ✔ P1 |
| Fractional snow cover | sub-grid patch coverage blend | — | ✔ P2 |
| Canopy snow interception / unload | intercepted snow on leaves, shed/unload | — | ✔ P2 |
| Blowing / drifting snow | wind redistribution, sublimation of blowing snow | — | out (§8) |

### 1.3 Reachability facts from the ED2 reference

- ED2 unifies snow and ponded liquid as **temporary surface water (TSW)**: up to `nzs` layers between the top soil layer and the CAS; the ice/liquid split is **diagnosed from internal energy every step**. Layer 1 = bottom (soil interface), `ksn` = top (atmosphere). A scalar **virtual layer** buffers sub-threshold precip/dew.
- **Unit gotcha (avoided by construction in MEDS):** ED2 stores `sfcwater_energy` as **J/kg** in the persistent csite state but **J/m²** inside the RK4 integrator (`rk4_copy_patch.f90:221,1301`). MEDS carries **extensive J/m²** everywhere — no dual-unit conversion (§13 S1).
- **Accumulation:** precip enters after canopy interception; the un-intercepted remainder always deposits onto the **top** layer (`rk4_derivs.f90:552`); if no layer exists, into the virtual pool, folded into layer 1 by `adjust_sfcw_properties`.
- **Precip partition** (`ed_met_driver.f90:2507-2562`, Jin 1999): rain-only above `t3ple+2.5`, snow-only below `t3ple`; mixed band in between; snow density `snden = 50 + 1.7·(T−258.15)^1.5` down to a `snden=50` floor. MEDS already reproduces the **liquid-fraction band** in `precip_phase` (`PHASE_BAND_K`), and computes the ice enthalpy internally (§4a).
- **Surface energy** (`rk4_derivs.f90:538`): `d(sfcwater_energy)(k) = h_flux_s(k) − h_flux_s(k+1) + rshort_s(k)`; top layer additionally gets `rlong_s`, committed SW of extinct layers, and latent/dew/precip enthalpy. Sensible to CAS `hflxsc = snowfac·ggnet·ρ·cp·(T_snow − can_temp)` (`canopy_derivs_two:1242`).
- **Sublimation is implicit in the enthalpy formulation** (`rk4_derivs.f90:1291`): vapor leaves at `tq2enthalpy8(T,1)`; removing vapor enthalpy from an ice layer automatically debits latent heat of **sublimation** (= vaporization + fusion) because the layer's own internal energy is ice-referenced. MEDS reproduces this with the vapor-enthalpy twin (§4d).
- **Melt/percolation** (`adjust_sfcw_properties`, top→bottom): "free" liquid above a holding capacity percolates down; at the soil interface it is capped by pore room; excess → runoff. Thresholds: `tiny_sfcwater_mass=1e-3`, `snowmin=5`, `water_stab_thresh=10 kg/m²`. A single layer thinner than `water_stab_thresh` is **forced into thermal equilibrium with the top soil layer** — the seed of the MVP thin-snow handling (§6).
- **Albedo** (`radiate_driver.f90:717`): TSW albedo interpolates damp-soil↔snow by top-layer liquid fraction; `snow_albedo_vis=0.518`, `snow_albedo_nir=0.435`, `snow_emiss_tir=0.970`. Snow-cover `snowfac` from Niu-Yang07 (`ny07_a=2.5`, `ny07_m=1.0`).
- **Insulation** (`rk4_derivs.f90:474`): soil/bottom-snow interface uses depth-weighted geometric-mean conductivity; low-density snow conductivity `≈0.1–0.3 W/m/K` sharply reduces the interface flux — the physical decoupling that caps the winter soil surface.

---

## 2. Where it lives (library DAG)

### 2.1 The state/process wall

```
shared ← {allometry, plant} ← state ← demography ← aux ← main
             ↑
        biophysics  (snow kernels: meds_snow_energy, meds_snow_mass — link meds_shared ONLY)
```

The snow kernels are **stateless per-store compute** in `src/biophysics/`, siblings of `plant`, linking `meds_shared` **only**. They **cannot** name `site_t`, `use meds_demography_types`, or read a global config; forcing arrives as value types (`snow_forcing_t`, `snow_params_t`, `snow_opts_t`).

**Controlling-precedent sentence.** A per-patch `snow_column_t` (SWE + internal energy [+ depth], per layer) is to `patch_index` what `soil_energy_column_t` and `soil_column_t` already are; the snow-surface skin temperature is a diagnosed read-off (like the current `t_ground = soil_temp(1)`), while the **SWE mass and column internal energy are the prognostic stores** advanced by an operator-split BE substep.

**Typed-vs-plain ownership rule (the recurring blocker).** `meds_config` is the DAG root in `src/shared` and must **not** `use meds_biophysics_types`; therefore `snow_params_t` is a **biophysics** type, never a `meds_config_t` component. `meds_config_t` carries only plain `real`/`integer` snow parameters (fresh-snow density, albedo endpoints, min-mass thresholds), assembled into `snow_params_t` by a biophysics `build_snow_params` routine called from the aux/init layer — the `build_soil_params`/`build_soil_thermal`/`derive_rad_optics` precedent.

**Selector enums** (`SNOW_ABSENT`/`SNOW_PRESENT`, `SNOW_PHASE_ON`) go in `meds_biophysics_types` beside the existing `SOIL_*`/`ENERGY_*` codes; TOML strings are mapped by `meds_config_io`. The snow melt plateau reuses the existing `ENERGY_PHASE_ON` precedent — snow ships with the plateau always on (there is no snow without phase change).

### 2.2 Files & CMake wiring

| File | Role | Analogue |
|---|---|---|
| `src/biophysics/meds_snow_energy.f90` | snow-surface balance + snow-column BE conduction + phase/melt read-off | `meds_column_energy.f90` (`ground_surface_balance`, `soil_energy_flux`) |
| `src/biophysics/meds_snow_mass.f90` | accumulation, sublimation sink, melt→free-liquid, percolation split, `snowfac` | `meds_column_hydrology.f90` (surface store + budget) |
| `src/shared/meds_column_state_types.f90` | **add** `snow_column_t` + `blend_snow` | `soil_energy_column_t`/`soil_column_t` + `blend_soil_e`/`blend_soil_w` |
| `src/shared/meds_thermo.f90` | **add** `internal_energy_ice` primitive (or reuse `uext_to_temp`) | `internal_energy_liquid` (`:87`) |
| `src/shared/meds_constants.f90` | snow phys constants (`fresh_snow_density`, `k_ice`, `k_air`, `snow_emiss`) | phase constants block |
| `src/biophysics/meds_biophysics_types.f90` | `snow_params_t`, `snow_opts_t`, `snow_forcing_t`, `snow_flux_t`; **add** `snow` to `patch_biophys_t` | `soil_params_t`/`energy_opts_t`/`energy_forcing_t` |
| `src/state/meds_demography_types.f90` | **add** `snow(:)` to `patch_index` + lockstep sites | `soil_e(:)`/`soil_w(:)` |
| `test/test_snow_energy.f90`, `test/test_snow_mass.f90` | standalone CTest per kernel + closure ledgers | `test_surface_energy.f90`, `test_column_hydrology.f90` |

`meds_soil_solver` (the hand-rolled `thomas_solve`) is reused **unchanged** for snow conduction.

### 2.3 The fast-loop gap and standalone testability

There is no general sub-daily loop wired for a kernel in isolation; each snow kernel ships on its own CTest with forcing passed as a value type and siblings lagged/forced (the energy/hydrology `§2.3` precedent). The snow-surface balance is tested against a prescribed CAS temp + radiation; the snow-mass kernel against a prescribed snowfall/melt schedule with a closed ledger; the coupled snow-in-column behavior is exercised through the split driver (§9, Jan diurnal run).

---

## 3. State variable & representation

### 3.1 The prognostic store and its read-offs

**Prognostic (persisted):** water-equivalent **mass** `swe [kg/m²]` and **internal energy** `snow_energy [J/m²]` (extensive, per ground area), per layer. MVP = one layer. Optionally **depth** `snow_depth [m]` is carried as a prognostic geometric state (needed for conduction thickness, `snowfac`, and albedo/optical depth); in the MVP single-layer path depth is derived from mass and a **fixed seasonal-mean bulk density `rho_snow ≈ 250 kg/m³`** (§4g — not fresh density; density only evolves via compaction in P1), so `snow_depth = swe/rho_snow` and there is no compaction in P0.

**Diagnosed each step (read-offs, NOT persisted authoritatively):** snow temperature `snow_temp [K]` and liquid fraction `snow_fliq [-]`, both from the shared inverter:

```
call uext_to_temp(snow_energy, swe, dry_hcap=0.0_wp, snow_temp, snow_fliq)
```

with the plateau logic (`meds_thermo.f90:59-67`): `u_freeze = swe·cp_ice·t_3ple`, `u_melt = u_freeze + swe·latent_heat_fusion`; below `u_freeze` all ice, above `u_melt` all liquid (a fully melted layer is meltwater to be drained), between them `snow_temp = t_3ple`, `snow_fliq = (snow_energy − u_freeze)/(swe·latent_heat_fusion)`. **Because `dry_hcap → 0`, the banded negligible-water guard (`abs(u_melt−u_freeze) < eps` when `swe → 0`) is load-bearing** — the last-flake transition must not divide by ~0 (§6, §13 S3).

**Why persist mass+energy, not temperature.** Internal energy is the conserved quantity; temperature is not (it pins at `t_3ple` through the entire melt). Persisting `(swe, snow_energy)` makes melt/refreeze a **smooth, monotone, invertible read-off** and makes conservation a machine-precision ledger. A prognostic-temperature snow model would have to bolt on a two-step phase adjustment (CLM5-style) — **rejected** in favor of reusing the one inverter, exactly as the soil column already does under `ENERGY_PHASE_ON`.

```fortran
!> src/shared/meds_column_state_types.f90  (new type, beside soil_energy_column_t)
integer(ik), parameter :: n_snow_layer_max = 1_ik    !< MVP single layer; P1 raises to ED2 nzs≈8

type :: snow_column_t
   real(wp) :: swe(n_snow_layer_max)         = 0.0_wp   !< [kg/m2] water-equivalent mass  (PROGNOSTIC)
   real(wp) :: snow_energy(n_snow_layer_max) = 0.0_wp   !< [J/m2]  extensive internal energy (PROGNOSTIC)
   real(wp) :: snow_depth(n_snow_layer_max)  = 0.0_wp   !< [m]     layer geometric depth (PROGNOSTIC, geom)
   real(wp) :: snow_temp(n_snow_layer_max)   = 0.0_wp   !< [K]     diagnosed each step (read-off)
   real(wp) :: snow_fliq(n_snow_layer_max)   = 0.0_wp   !< [-]     diagnosed liquid fraction (read-off)
   integer(ik) :: nlayer                     = 0_ik     !< active layer count (0 = no snow)
end type
```

**Each component carries its own default initializer** (Fortran initializes only the last name on a shared-initializer line) to stay clean under `-fpe0`/`-Ktrap=fp`. `nlayer = 0` is the "no snow" state — the store is present but empty, and the surface reverts to the bare-soil skin (§5).

**Conservation on fusion/disturbance.** `blend_snow` area-weight-mixes `swe` and `snow_energy` (extensive → additive-per-area, blend like `w_surface`/`soil_energy`); `snow_depth` blends by area weight; `snow_temp`/`snow_fliq` are **re-diagnosed** post-blend from the blended `(swe, snow_energy)`, never blended directly. A disturbance-gap patch seeds `nlayer=0`, `swe=0` (bare ground under a new gap; snow that was on the disturbed area is conserved into the surviving patches via the blend, or accounted as a mass sink if the gap model destroys it — MVP: conserve into survivors).

**One energy datum across all condensed-phase stores (load-bearing for conservation).** The snow store, the soil energy column, and the CAS vapor enthalpy MUST share ONE internal-energy reference, or every snow↔soil↔CAS transfer leaks silently per step. The convention (inherited from `meds_thermo`): the condensed phases are **ice-referenced** — `internal_energy_ice(T)=cp_ice·T` and `internal_energy_liquid(T)=cp_ice·t_3ple + latent_heat_fusion + cp_liq·(T−t_3ple)` share the same 0-K ice datum, and `uext_to_temp`/`temp_to_uext` are their exact inverse; the vapor enthalpy `enthalpy_vapor(T)` used to debit the pack for sublimation is the **same** function the CAS credits for the vapour it receives. Because every transfer in this design is a raw-additive move of `(mass, its enthalpy)` between two stores on this one datum — melt→soil (§4e), melt-out dump→soil (§4c), `g_base` conduction (§4g), sublimation→CAS (§4d) — conservation is a machine-precision round-trip. Test 8 (§9) is a raw-additive snow→soil→snow round-trip that asserts this.

**Initialization & restart.** Cold-start is the **empty store** (`nlayer=0`, `swe=snow_energy=snow_depth=0`) — no snow at model start; the pack builds from the first sub-freezing precip. Warm-start/restart serializes the three prognostic fields `(swe, snow_energy, snow_depth)` + `nlayer` per patch through the state I/O (the `soil_energy`/`soil_water` restart precedent); the two read-offs (`snow_temp`, `snow_fliq`) are re-diagnosed on read, never serialized authoritatively (like the soil column's diagnosed temperature).

---

## 4. Physics kernels

All kernels are `pure`/`elemental` where the store update is not tridiagonal; the conduction step is a `subroutine` reusing `thomas_solve`.

### 4a. Accumulation from snowfall (frozen-precip mass + ice enthalpy)

Snowfall `snowf [kg/m²/s]` arrives from `precip_phase` (already computed; §5 threads it). It deposits onto the top layer, or creates layer 1 if `nlayer=0`:

```
swe_new         = swe(top) + snowf·dt
snow_energy_new = snow_energy(top) + snowf·dt · internal_energy_ice(t_snowfall)
```

where `t_snowfall = min(t_3ple, tair_k)` (frozen precip cannot be warmer than freezing; ED2 uses `tl2uint(min(t3ple,T),0)`). **The one missing thermo primitive** is `internal_energy_ice(T) = cp_ice·T` (the ice-referenced specific internal energy `[J/kg]`), the frozen twin of the existing `internal_energy_liquid` (`meds_thermo.f90:87`); equivalently call `temp_to_uext(0, 1, T, fliq=0)`. Add it to `meds_thermo`. Depth is added at the P0 **seasonal-mean bulk density** `rho_snow` (§4g — NOT fresh density, to avoid over-deepening an aged pack): `snow_depth += snowf·dt / rho_snow`. **Fresh liquid precip on an existing pack (rain-on-snow, §4e) adds to the top layer's `swe` + `internal_energy_liquid(tair_k)` and refreezes via the inverter — it is not a separate branch here.**

**First-flake creation.** If `nlayer=0` and `snowf·dt ≥ min_new_snow_mass`, set `nlayer=1` and stamp `snow_energy` from `temp_to_uext(0, swe, t_snowfall, 0)`. Below `min_new_snow_mass`, snowfall is held/added to a scalar buffer (the MVP analogue of ED2's virtual layer; §6) or, simplest, folded directly into `soil_energy(1)` as a trace mass+enthalpy so no mass is lost. **Decision:** MVP folds sub-threshold new snow into the soil-top store (mass + ice enthalpy), avoiding a separate virtual-pool state; the buffer is a P1 refinement.

### 4b. The snow-surface ENERGY balance (the seam snow replaces)

When snow is present, the CAS-facing surface is the **snow surface** at temperature `t_snow_surf` (= `snow_temp(top)`), and the balance closes on the snow store, not the soil top. The algebra mirrors the inlined `ground_surface_balance` but with the snow-base conduction term:

```
rn_snow    = abs_sw_snow + abs_lw_snow                        ! net SW (snow albedo) + net LW (ε_snow)
h_snow     = ggnet · rho_air · cp_air · (t_snow_surf − can_temp)      ! sensible to CAS  [W/m2]
w_flux_snow= ggnet · rho_air · (sat_specific_humidity(t_snow_surf, press) − can_shv)  ! [kg/m2/s]
le_snow    = w_flux_snow · L_surface(t_snow_surf, snow_fliq_surf)     ! sublimation OR evaporation (§4d)
g_base     = k_snow_soil · (t_snow_base − t_ground_soil) / dz_interface   ! conduction into soil top [W/m2]
d(snow_energy)/dt = rn_snow − h_snow − le_snow − g_base
```

- `h_snow`/`le_snow` feed the CAS enthalpy source unchanged (`src_enth = coh_h + coh_qw + h_snow + le_snow`), and `w_flux_snow` feeds `src_vap`. The single ground-latent authority moves from soil evaporation to snow sublimation when snow is present (no double-count).
- **`g_base` is the new soil top BC.** It replaces `g_top = Rn − H − LE` (bare soil) with the snow-base conduction flux (§4g). `soil_energy_flux` receives `eforc%g_top = g_base` unchanged.
- **`abs_sw_snow`/`abs_lw_snow`** come from the RT path with snow albedo/emissivity (§4f), landing on the snow surface temperature for correct LW emission.

The whole-column update advances `snow_energy` by an operator-split backward-Euler substep. For the single-layer MVP the "conduction column" is degenerate (one interior node); the surface fluxes are quasi-steady on `t_snow_surf` and the store integrates `d(snow_energy)/dt` over `dt_fast` with the plateau read-off applied after (melt is diagnosed, not integrated as a stiff source — see §4c). This reuses the L-stable, bounded-flux discipline that fixed the wood store (`cap·(t_star − t_n)` + step-averaged flux, never unbounded `dt·r_star`).

### 4c. Freeze/thaw + MELT (internal-energy phase partition)

Melt is **not a source term in the ODE** — it is a read-off after the energy update, exactly the soil freeze/thaw pattern:

1. Advance `snow_energy` by `d(snow_energy)/dt·dt` (§4b).
2. Invert: `call uext_to_temp(snow_energy, swe, 0, snow_temp, snow_fliq)`.
3. **Melt = the liquid fraction that becomes "free" and drains** (§4e). If `snow_fliq > 0`, the layer is on the plateau (or fully liquid); the free-liquid fraction above the holding capacity becomes meltwater removed from `swe` and `snow_energy` (carrying `internal_energy_liquid(t_3ple)`), routed to soil infiltration.
4. **Refreeze** is automatic and needs no branch: if `d(snow_energy)/dt < 0` on a wet layer, `snow_energy` drops back below `u_melt`, `snow_fliq` falls, and the plateau reabsorbs latent heat — the zero-curtain, captured for free by the inverter.

**Full melt-out.** If, after draining free liquid, `swe ≤ tiny_snow_mass`, the residual mass+energy is dumped to `soil_energy(1)`/soil water (conserving both), `nlayer→0`, and the surface reverts to bare soil next step (§6 last-flake). **Documented bound:** when a layer melts fully mid-`dt_fast`, its surface fluxes (H/LE/LW) were evaluated at the snow skin for that substep (RT is lagged one substep, S6), so up to one sub-step of surface flux is mis-attributed on the melt-out step — an `O(dt_fast)` transient the adaptive step-doubling bounds; acceptable in P0, no intra-step re-solve.

### 4d. Sublimation / evaporation (the vapor twin)

The latent flux direction is set by the sign of `(sat_specific_humidity(t_snow_surf) − can_shv)`:

- **Sublimation** (frozen surface, `snow_fliq_surf < 1`, CAS undersaturated): vapor leaves at **sublimation enthalpy**. Following ED2's implicit-enthalpy trick, the vapor departs at vapor enthalpy `L_surface = enthalpy_vapor(t_snow_surf)`, and because the layer's own internal energy is **ice-referenced** (via `internal_energy_ice`), removing vapor enthalpy from an ice layer automatically debits latent heat of sublimation = vaporization + fusion. **The twin rule (bug E8):** move mass `w_flux_snow·dt` out of `swe` and enthalpy `w_flux_snow·dt·enthalpy_vapor(t_snow_surf)` out of `snow_energy`, **same T_source**. No separate `L_sublimation` constant is needed — the ice reference supplies the fusion component.
- **Evaporation** (wet surface, `snow_fliq_surf > 0`): identical form; the liquid reference (a wet layer's energy is above `u_freeze`) supplies only vaporization. `L_surface = enthalpy_vapor(t_snow_surf)` in both cases; the difference is which side of the plateau the layer sits on, handled automatically by the shared reference.
- **Deposition/frost** (`sat_ssh ≤ can_shv`): reverse sign — vapor deposits as frost, adding ice mass+enthalpy to `swe`/`snow_energy` **when a layer already exists** (`nlayer≥1`; MVP symmetric with sublimation, an explicit signed source in the ledger, tests 1–2). **Frost onto bare soil (`nlayer=0`) does NOT spawn a snow layer in P0** — it stays on the existing soil-dew path (`ground_evaporation` handles condensation onto the soil today); a deposition-initiated pack is a P1 refinement (alongside ED2's virtual pool).

`w_flux_snow` is capped to `swe/dt` (cannot sublime more than exists), FPE-safe; the deposition (negative) direction has no such cap but IS entered in the mass+energy ledgers as a signed source so it cannot leak.

### 4e. Percolation of meltwater into the soil

Free (percolatable) liquid uses the LEAF-3 / ED2 `ipercol=0` holding-capacity rule for the MVP:

```
wmass_free = max(0, swe·(snow_fliq − liquid_holding_frac)/(1 − liquid_holding_frac))
```

with `liquid_holding_frac ≈ 0.10` (ED2's 1:9 liquid:ice shed).

**Meltwater is a PAIRED (mass, enthalpy) transfer — the single most important conservation point.** `wmass_free` is removed from `swe`, **and the matching enthalpy `wmass_free·internal_energy_liquid(t_3ple)` is removed from `snow_energy` AND added to `soil_energy(1)`** in the same step. Routing only the mass (and dropping the enthalpy) would leak the meltwater's heat every melt event and, worse, would prevent the soil from **refreezing** percolating melt (the latent release that warms a cold soil under a melting pack). The transfer is a closed snow→soil handoff:

```
melt_mass   = wmass_free                                     ! [kg/m2] leaves swe
melt_enth   = wmass_free · internal_energy_liquid(t_3ple)    ! [J/m2]  leaves snow_energy, enters soil_energy(1)
swe        -= melt_mass ; snow_energy -= melt_enth
soil_energy(1) += melt_enth / dz(1)                          ! extensive->per-volume for the soil store
q_avail     = w_surface/dt + q_liq + melt_mass/dt            ! meds_column_hydrology.f90:137 insertion
```

The ledger **asserts snow-loss == soil-gain** (`melt_enth`) as a hard test (§9 test 2). The mass then flows through the **existing** conductivity-limited infiltration → ponding → runoff logic (`meds_column_hydrology.f90:135-197`) with **no new outflow path** — an internal store-to-store transfer (snow → soil water), like `clip_ex` routing to ponding, so melt **does not appear** in the water mass budget as an in/out term (§9). The Anderson-1976 capacity scheme (`ipercol=1`) is a P1 config option.

**Rain-on-snow (a dominant real melt/energy driver, must not be dropped).** When `nlayer ≥ 1`, un-intercepted `rainf` deposits onto the **top snow layer** (mass + `internal_energy_liquid(tair_k)`), NOT straight to the soil — matching ED2, which lands rain on the top TSW layer. Warm rain on a sub-freezing pack then **refreezes automatically** via the inverter (`snow_energy` rises but stays below `u_melt` → `snow_fliq` partial → the plateau reabsorbs the rain's sensible heat AND releases fusion latent heat that warms the pack toward melt). Only the resulting *free* liquid (§4e holding-capacity rule) drains to the soil. When `nlayer = 0`, rain bypasses to infiltration as today (bit-identical). This adds `rainf` to the forcing threaded in §5.

### 4f. Snow ALBEDO fed into ground radiation

The RT path (`apply_rt_forcing`, `meds_fast_loop.f90:496-576`) already assembles ground optics via `ground_optics(surf, ...)` from `ctx%soil_albedo(3)` and `ctx%soil_emiss`, then scatters net SW/LW back into `forc%abs_sw_ground`/`forc%abs_lw_ground`. **Snow raises exactly these inputs.** When snow is present:

- `surf%soil_albedo(VIS) = snow_albedo_vis`, `surf%soil_albedo(NIR) = snow_albedo_nir` (fresh values `≈0.85 VIS / 0.65 NIR` fresh, ED2's `0.518/0.435` are aged endpoints — MVP interpolates fresh↔aged by a simple age counter or by top-layer `snow_fliq`, ED2-style: all-ice→snow albedo, wet→lower).
- `surf%soil_emiss = snow_emiss (≈0.97)`.
- `surf%soil_temp = t_snow_surf` (not `soil_temp(1)`) — so LW emission comes off the snow surface.

Higher `grnd_refl` (snow) automatically lowers `abs_sw_ground`; `surf%soil_temp` adjusts net LW. Both flow straight into the §4b balance with no other RT change. Snow albedo/emiss live on `fast_context_t` (or a new per-patch snow-optics field) as a function of snow age/temperature; **all parameters from TOML**, none hard-coded (§8). The empty-canopy short-circuit in `meds_canopy_radiation.f90:51-59` already handles a bare snow patch (no cohorts).

**Ramp the optics by `snowfac`, do NOT flip binary at `nlayer≥1` (continuity, critical).** A trace flake (`min_new_snow_mass ~1e-3 kg/m²`) must not step ground albedo soil→0.85 and emissivity soil→0.97 in one substep — that is a discontinuity in net ground SW/LW. Blend by the same Niu-Yang07 burial fraction `snowfac ∈ [0,1]` (SWE/depth-driven) that the aerodynamic seam already uses:

```
surf%soil_albedo(b) = (1−snowfac)·soil_albedo(b) + snowfac·snow_albedo(b)      ! b = VIS, NIR
surf%soil_emiss     = (1−snowfac)·soil_emiss     + snowfac·snow_emiss
surf%soil_temp      = (1−snowfac)·soil_temp(1)   + snowfac·t_snow_surf          ! emitting-surface temp
```

so first-flake → full-cover → melt-out are all continuous in the radiation, matching the aero ramp. Using one `snowfac` for aero + radiation + the soil-BC blend (§4g) keeps a single cover authority.

**Cohort snow-burial (short-stature bias, documented).** P0 does **not** re-route the shortwave of cohorts shorter than the snow depth to the bright-ground path — buried cohorts still absorb full SW in the two-stream (§8). For this low-LAI winter stand most canopy is taller than a P0 pack, so the bias is small, but for cohorts with `height < snow_depth` it over-absorbs / under-reflects. P0 documents the bound; reassigning buried-cohort SW to `snow_albedo` is the first step of P2 fractional/burial optics.

### 4g. Thermal INSULATION (snow↔soil conduction replacing the bare-soil top BC)

Today `g_top = Rn − H − LE` conducts the full surface energy balance straight into `soil_energy(1)`. With snow, the soil top sees only the **snow-base conduction**:

```
k_snow      = snow_conductivity(rho_snow)          ! MVP: fixed k_snow from a SEASONAL-mean density (below)
dz_interface= 0.5·snow_depth(1) + 0.5·|z_node_soil(1)|
k_face      = depth-weighted geometric-mean(k_snow, soil_thermal_cond(1))    ! ED2 rk4_derivs.f90:474
g_base_snow = k_face · (t_snow_base − soil_temp(1)) / dz_interface
g_base      = (1−snowfac)·g_top_bare + snowfac·g_base_snow                    ! smooth blend, no threshold cliff
```

Because `k_snow ≪ k_soil`, a thick pack sharply throttles `g_base_snow` — the physical decoupling that **caps the winter soil surface near `t_3ple`** instead of letting it heat to 26 °C on a clear day.

**Blend the soil-top BC by `snowfac`, do not switch it at `snow_stab_thresh` (continuity, high).** The forced-equilibrium threshold (§6) governs the snow *layer's* stiffness, but the *soil top BC* must transition smoothly from the bare-soil `g_top_bare` to the insulated `g_base_snow` as cover grows — otherwise crossing `snow_stab_thresh` (~0.1 m) jumps the soil-top flux, and the melt-out step (`nlayer→0`) jumps it back. The `snowfac` blend (same authority as the §4f optics) makes first-flake, threshold, and last-flake all flux-continuous. The thin-persistent-snow regime (0 < SWE < `snow_stab_thresh`) then gets partial, continuous insulation rather than none.

**Sign convention (pin it, medium).** All fluxes are positive **upward**; `g_base` is the flux from the soil into the snow base. A **warm snow base** (`t_snow_base > soil_temp(1)`) gives `g_base > 0` upward-out-of-soil is wrong — so the soil sees heat flowing **down** into it: fed to `soil_energy_flux` as `eforc%g_top` with the existing `hf(0) = −g_top` (`meds_column_energy.f90`), a warm base **warms** the top soil node, and the identical `g_base` is a **loss** from the snow bottom layer. Worked case: `t_snow_base = 271 K`, `soil_temp(1) = 273 K` → `g_base < 0` → heat flows **up** from the warmer soil into the colder snow base (soil loses, snow gains) — the winter norm that keeps the pack from freezing solid to the ground. Tests 2 and 6 assert `|snow-loss| == |soil-gain|` and the correct direction. **Conservation by construction:** the same `g_base` debits the snow bottom layer and credits the soil top, counted once with opposite signs.

**Density in P0 (medium).** MVP uses a **seasonal-mean bulk density** (`rho_snow ≈ 250 kg/m³`), NOT fresh-snow density (~100 kg/m³): a season-old pack at fresh density would be 2–4× too deep → `g_base` too throttled and `snowfac` too high, so the January insulation test could pass for the *wrong* reason (excess insulation) before P1 lands. Depth then tracks `swe/rho_snow`. Density-dependent conductivity + prognostic compaction (ED2's Medvigy-2007 `ss` polynomial) is P1.

---

## 5. The seams to modify

**Who owns the surface: the snow-or-soil skin.** The rule is a single branch on `nlayer`:

- `nlayer == 0` (no snow): the surface skin is `soil_temp(1)`; the bare-soil `ground_surface_balance` runs unchanged; `g_top = Rn − H − LE_soil`. **Zero behavior change** from today when snow is absent (the Jan run before accumulation is bit-identical).
- `nlayer ≥ 1` (snow present): the surface skin is `t_snow_surf`; the snow-surface balance (§4b) runs; the soil top BC becomes the `snowfac`-blended `g_base` (§4g); ground optics and the latent authority (sublimation vs soil-evaporation) are **blended by `snowfac`** (§4f), not switched binary — so the snow↔bare-soil handoff is continuous at first-flake and melt-out. At full cover (`snowfac→1`) the surface is pure snow; at a trace flake (`snowfac→0`) it is essentially bare soil.

Exact files/routines:

| Seam | File:line | Change |
|---|---|---|
| **Forcing** — thread `snowf` | `meds_fast_loop.f90:482` (`apply_met_to_ctx`) → `fill_forcing:450` → `meds_column_dynamics.f90:477/896` | add `ctx%snow ← met%snowf`, `forc%snowf`, `hforc%snowf`; mirror the three existing `precip` assignments |
| **State** — working bundle | `meds_biophysics_types.f90:383` (`patch_biophys_t`) + `alloc_patch_biophys:442` | add `type(snow_column_t) :: snow`; gather `bio%snow = site%patch%snow(ip)` (`meds_fast_loop.f90:281`), scatter back (`:367`) |
| **State** — persistent + lockstep | `meds_demography_types.f90:130` | add `snow(:)`; update all 6 lockstep sites (alloc `:235`, dealloc `:195`, `patch_ensure_capacity:347/356`, `sort_patches`, `patch_compact`, fusion-blend `blend_snow`) + disturbance-gap seed |
| **Energy** — surface balance | `meds_column_dynamics.f90:503-506` (split, inlined) | branch on `nlayer`: snow-surface balance (§4b), `t_ground`→`t_snow_surf` for `:504-505`, `g_top`→`g_base` for `:526` |
| **Energy** — kernel | `meds_column_energy.f90:307` (`ground_surface_balance`) | add sibling `snow_surface_balance` (kernel + standalone test); the live inlined copy is the authority |
| **Soil top BC** | `soil_energy_flux` / `soil_heat_be_step` (`meds_column_energy.f90:44,201`) | **unchanged** — receives `g_base` as `g_top` transparently |
| **Hydrology** — melt + sublimation + accumulation | `meds_column_hydrology.f90:75` | add `snowf`/`melt` to `q_avail:137`; add `snow_sublimation` sink (pattern of `ground_evaporation:538`); snapshot+update `col%w_snow` OR the `snow_column_t` swe; extend `mass_resid:214` |
| **Radiation** — ground albedo | `meds_fast_loop.f90:553-559` (`apply_rt_forcing`) | when snow present set `surf%soil_albedo/emiss/temp` to snow optics (§4f) |
| **Aero coupling** — burial | `meds_fast_loop.f90:274` (`ageom%snowfac = 0.0`) | drive `snowfac` from SWE/depth (Niu-Yang07, §4f); the roughness blend (`meds_canopy_aerodynamics.f90:50-53,77`) already consumes it |
| **Coupling** — ARK RHS | `meds_column_derivs.f90:256-258` (`surface_derivs`) + `surface_frozen_t:60` | operator-split: snow committed once post-march (§6); `surface_frozen_t` gains `t_snow_surf`/`g_base`; `fro%t_ground`→snow-surface when present |
| **Budget** — whole-column | `meds_column_dynamics.f90:597-615` | add `swe`+`snow_energy` to `whole_water`/`whole_energy` store totals; add `snowf` to `e_in`/`w_in`, sublimation to `w_out`/`e_out` (mirror `wood_store0/1`) |

---

## 6. The thin-snow problem & the operator-split-vs-ARK decision

**First-flake / last-flake transitions.** The `nlayer` 0↔1 boundary is the delicate seam. Two guards:

- **Minimum new-snow mass** `min_new_snow_mass [kg/m²]`: below it, accumulating snowfall is folded into the soil-top store (mass + ice enthalpy) rather than creating a layer (§4a) — avoids a flickering micro-layer and a divide-by-tiny-`swe` in the inverter.
- **Melt-out threshold** `tiny_snow_mass [kg/m²]` (ED2 `tiny_sfcwater_mass=1e-3`): once `swe` drops below it, the residual mass+energy is dumped to soil and `nlayer→0`. The banded negligible-water guard in `uext_to_temp` (`abs(u_melt−u_freeze)<eps`) must fire before `swe` reaches machine-zero — **use a banded comparison, never `swe==0`** (FPE-safe, §13 S3).

**Thin-layer thermal stability (ED2's `water_stab_thresh=10 kg/m²`).** A snow layer thinner than a stability threshold has too little heat capacity to integrate independently against surface fluxes over `dt_fast` — it oscillates. **MVP rule (ED2-faithful):** when `0 < swe < snow_stab_thresh`, the thin layer is **forced into thermal equilibrium with the top soil layer** — combine `(swe·energy + soil-top energy)` into a single pool, invert once with `uext_to_temp`, share the resulting temperature/fliq. This is exactly ED2's `flag_sfcwater==1` path (`rk4_misc.f90:1299`) and it dissolves the stiffness: a thin cold film simply tracks the soil skin. Above `snow_stab_thresh`, the snow layer integrates on its own with the throttled `g_base` conduction (real insulation).

**Which surface balance drives the combined pool (ownership, continuity).** The forced-equilibrium pool is thin, so its *surface* fluxes must be continuous with BOTH the bare-soil endpoint (as `swe→0`) and the independent-snow endpoint (as `swe→snow_stab_thresh`). It is driven by the **`snowfac`-blended** surface balance, NOT full snow optics: absorbed radiation uses the `snowfac`-ramped albedo/emissivity (§4f), and the latent authority is `(1−snowfac)`·soil-evaporation + `snowfac`·sublimation. So a 0.001-kg film does not impose full snow albedo on the soil's heat capacity — it nudges the shared pool by an amount proportional to cover. This makes the whole `swe ∈ [0, snow_stab_thresh]` regime continuous at both ends, and the same `snowfac` authority governs the pool's radiation, its latent split, and (via §4g) its share of the soil-BC blend.

**Fractional cover (`snowfac`) is the single continuity authority.** MVP computes `snowfac` (Niu-Yang07, `ny07_a=2.5`, `ny07_m=1.0`) from SWE/depth and uses it to **ramp** three seams so the snow↔bare-soil handoff is everywhere continuous: the aerodynamic burial blend (roughness), the ground optics (§4f), and the soil-top-BC blend + latent split (§4g). This is a deliberate strengthening over the naive "binary at `nlayer≥1`": a trace flake nudges each seam by `snowfac≈0` (essentially bare soil), full cover gives `snowfac≈1` (pure snow), with no step at first-flake, `snow_stab_thresh`, or melt-out. What `snowfac` does NOT yet do in P0 is **sub-patch area partitioning** of the two-stream (running a separate snow-covered and bare-ground radiative sub-column and area-averaging) — that finer fractional-radiation treatment is P2; P0's single ground surface simply carries `snowfac`-weighted optics.

**Operator-split vs ARK — recommendation: operator-split the snowpack; keep it OUT of the ARK implicit (ESDIRK) set in P0.** Reasoning, all grounded in existing code:

1. **Precedent is explicit and unanimous.** Soil water is already operator-split out of the ARK and committed once (`y_out%theta = fro%theta1`, `meds_column_dynamics.f90:687`); prognostic wood is operator-split (`:448-466`); the ARK path **hard `error stop`s** on prognostic wood/leaf today (`meds_column_derivs.f90:242`). A snow SWE+energy reservoir is the same class of add-on.
2. **The melt plateau is a switch-like (zero-curtain) nonlinearity.** The ARK2 embedded-error step-doubling controller needs a **smooth** RHS — the code deliberately smoothed the condensation clamp into a `TAU_COND` relaxation for exactly this reason (`meds_column_derivs.f90:47,266`). A hard melt threshold in the implicit set would wreck the controller. Melt as a **post-step read-off** (§4c), outside the implicit solve, sidesteps this entirely.
3. **The snow surface skin is near-massless/stiff** (thin fresh snow has tiny heat capacity) and is best treated like the current ground skin — quasi-steady/diagnostic inside the per-substep Picard loop — while the SWE **mass** and **column internal energy** are prognostic stores advanced by an operator-split backward-Euler substep reusing the freeze/thaw plateau. This is the wood-store precedent (bounded `cap·(t_star−t_n)` + step-averaged flux, L-stable — the fix that stopped the wood blowup/hydrology hang).

**Concrete split insertion (P0):** add a `snow_surface_balance` + `snow_mass_energy_update` block between §3c (`:503-506`) and §3d′ (`:515`) in `column_fast_step`. When snow present, redefine `t_ground = t_snow_surf` for `:504-505`, set `le_ground = le_snow` (sublimation), replace the `g_top` handed to `soil_energy_flux` at `:526` with `g_base`. `src_enth`/`src_vap` (`:509-510`) carry `h_snow+le_snow`/sublimation into the CAS unchanged. In the ARK path, add `t_snow_surf`/`g_base` to `surface_frozen_t` and commit the SWE/snow-energy update **once** post-march (mirroring the `theta1` commit + the `build_column_frozen` scratch-solve pattern at `:900-923`). **ARK-implicit snow is deferred and gated** (§7 "later").

---

## 7. Phasing (MVP → full)

- **P0 — single-layer snow, split path (this doc's deliverable).** One `snow_column_t` layer. Accumulate from `snowf` (+ rain-on-snow refreeze); snow-surface energy balance (net SW with fresh/simple-aged albedo, net LW, sensible to CAS, sublimation/evaporation, base conduction); freeze/thaw + melt as `uext_to_temp` read-off; **meltwater as a paired (mass, enthalpy) transfer** into existing infiltration; snow albedo/emissivity/soil-BC/latent split all **`snowfac`-ramped** (continuous, not binary) into `ground_optics`/aero; snow-base conduction as the soil top BC (insulation); seasonal-mean density. Thin-snow forced-equilibrium below `snow_stab_thresh`. Split/Picard driver only; ARK guarded to snow-free. Whole-column mass+energy ledgers extended. **Tests:** §9 1–10.
- **P1 — multi-layer + compaction/densification + aging albedo.** Raise `n_snow_layer_max` to ED2's `nzs≈8`; port the `thick`-table redistribution/merge (`adjust_sfcw_properties:1833`) conserving mass/energy/depth; prognostic bulk density with compaction; density-dependent snow conductivity (Medvigy-2007 `ss` polynomial); fresh→aged albedo decay curve; Anderson-1976 percolation option. Snow conduction becomes a real BE-Thomas column via `meds_soil_solver`.
- **P2 — fractional snow cover + canopy snow interception/unload.** Sub-patch `snowfac` blending in radiation and turbulent exchange (soil path `(1−snowfac)`, snow path `snowfac`); canopy snow interception (intercepted frozen precip on leaves), sublimation of intercepted snow, and unload/shed to the ground pack — wires `intercept_canopy_layer` (`meds_column_hydrology.f90:48`, currently test-only) into the loop.
- **later (gated) — blowing/drifting snow; ARK-implicit snow.** Wind redistribution + blowing-snow sublimation. Promotion of the snow store into the ESDIRK implicit set (add to `column_state_t`/`column_tend_t`) — deferred, gated on a smoothed melt formulation compatible with the step-doubling controller.

**Reconciling the parents' phase labels.** The energy doc parked snow at **P2b** (gated on a hydrology sfcwater stack); the hydrology doc parked it at **P5**. This module unifies them: it needs both a mass/SWE store (hydrology's territory) and a thermal store (energy's), so it is delivered as its own family with **P0 = the unified single-layer MVP**, lifting the energy-P2b and hydrology-P5 restrictions together.

---

## 8. Out of scope

- **Canopy snow load / interception on P0.** Intercepted snow on leaves and its unload/sublimation are P2; P0 sends all frozen precip straight to the ground pack (interception deferred, as liquid precip already is).
- **Blowing / drifting snow.** Wind redistribution and blowing-snow sublimation are out.
- **Snow–vegetation shortwave masking beyond ground albedo.** Snow raises **ground** albedo/emissivity only; it does not mask cohort crowns in the two-stream (no "snow-buried canopy" optics). Aerodynamic burial (`snowfac` roughness blend) is the only vegetation-facing snow effect in P0.
- **Glaciers / firn / permanent snow.** No perennial pack, no firn densification, no ice-sheet coupling.
- **Sub-`t_3ple` supercooled liquid, wet-snow metamorphism microphysics, snow grain-size radiative aging** beyond the simple fresh↔aged albedo interpolation.

---

## 9. Validation & milestones

Every kernel closes a **machine-precision budget** as its strongest invariant; the residual is a diagnostic out-arg and `error stop`s in Debug if `|resid| > atol`. Numbered CTest assertions:

1. **Mass closure (snow store).** `Δ(swe + w_surface + Σθ_k·dz_k·ρ_w + w_aquifer) = dt·(snowf + rainf − sublimation − soil_evap − uptake − runoff − drainage)` to round-off across a mixed accumulate+melt schedule. **Melt is absent from the RHS** (internal snow→soil transfer, §4e); its presence in the residual is a bug.
2. **Energy closure (snow store).** `Δ(snow_energy + soil column energy + CAS enthalpy + leaf/wood stores) = dt·(net rad + snowfall enthalpy − atm exchange − sublimation enthalpy − drainage/runoff enthalpy)` to `atol`. Sublimation enthalpy uses the **same T_source** twin (§4d).
3. **Cold-snap accumulation test.** Drive constant `snowf` at `tair < t_3ple−5` with weak radiation: `swe` grows linearly, `snow_temp < t_3ple`, `snow_fliq = 0`, `g_base` throttles toward zero as depth grows, both ledgers closed.
4. **Melt-event test.** From a cold pack, apply a warm radiative forcing: `snow_energy` climbs to the plateau (`snow_temp` pins at `t_3ple`, `snow_fliq` rises 0→1), free liquid drains as **melt → infiltration** (soil `θ(nzg)` or ponding rises by the drained SWE), energy conserved through the plateau, `swe` monotonically decreases to melt-out, `nlayer→0`, surface reverts to bare soil.
5. **The January diurnal run — the headline acceptance test.** With a snowpack present, the model must **cap the soil surface near `t_3ple`** through a clear-sky diurnal cycle instead of the current bare-soil skin heating to ~26 °C. Assert `soil_temp(1)` stays within a few K of `t_3ple` under snow (high albedo cuts absorbed SW; `g_base` insulation throttles conduction), versus the bare-soil control run.
6. **Isothermal cold snowpack sanity.** A deep cold pack with no radiation/turbulent forcing holds constant `(swe, snow_energy, snow_temp)` across many substeps (no spurious drift); `g_base ≈ 0` when `t_snow_base = soil_temp(1)`.
7. **Diagnostic-vs-forcing air-temp relationship.** Under snow, `t_snow_surf` tracks between `can_temp` (CAS) and the radiative equilibrium, never exceeds `t_3ple` while ice remains, and the sensible flux `h_snow` has the correct sign (snow surface colder than CAS on a cold clear night → downward sensible). Cross-check `snowfac` burial lowers `ggnet` toward `ggbare` as depth grows.
8. **Energy-datum round-trip (single-datum invariant, §3.1).** Move a known `(mass, enthalpy)` parcel snow→soil→snow via the §4e melt + a synthetic refreeze path; assert the returned `(swe, snow_energy)` equals the original to round-off — proves the ice/liquid/vapor references are one datum and no transfer leaks.
9. **Snow↔bare-soil continuity (first-flake / threshold / melt-out).** Ramp `swe` slowly through `0 → min_new_snow_mass → snow_stab_thresh → full cover` and back to melt-out; assert `abs_sw_ground`, `soil_temp(1)`, and the soil-top flux (`g_base` blend) are **continuous** (no step at any threshold) because optics + BC + latent split all ramp by the single `snowfac` (§4f, §4g, §6). The bare-soil control (`nlayer=0`) branch stays bit-identical to today.
10. **Rain-on-snow refreeze (§4e).** Warm rain onto a sub-freezing pack: the pack's `snow_energy` rises and `snow_temp` climbs toward `t_3ple` as fusion latent heat is released by refreezing (mass added to `swe`, not passed to soil); only free liquid above the holding capacity drains; both ledgers close.

**Non-convergence contract.** On `nsub → max_substep`, set `flux%converged=.false.`, log, and `error stop` under `debug_error`; production never writes back an unconverged capped snow state unchecked.

**Phased milestones.**

| Phase | Deliverable | Tests |
|---|---|---|
| **P0** | single-layer snow: accumulate + surface balance + sublimation + melt→infiltration + `snowfac`-ramped albedo + insulation, split path | 1–10 |
| **P1** | multi-layer stack + compaction/densification + aging albedo + density-dependent conductivity | multi-layer mass/energy closure; redistribution conservation |
| **P2** | fractional snow cover + canopy snow interception/unload | sub-patch blend closure; interception mass balance |
| **later** | blowing snow; ARK-implicit snow (gated on smoothed melt) | controller stability under melt |

---

## 10. GPU / nvfortran portability

- **Two-layer structure.** Derived-type seam (host) wrapping a bare-array, `firstprivate`-scalar, device-eligible inner step — the `soil_heat_be_step`/`growth_step` precedent. Snow conduction (P1) reuses `thomas_solve` (sequential within a column); the parallel axis is **columns/patches, never within a snow column**.
- **Fixed compile-time layer count.** `n_snow_layer_max` (1 in P0, ≈8 in P1) — the `n_soil_layer_max=20` analogue of ED2's `nzs`; scratch as fixed-length automatics.
- **nvfortran issue #7 trap.** Never pass an array-valued function result straight into a call; `thomas_solve` is a `subroutine` with `intent(out) x(...)`; bind array results to named arrays first. **Build nvfortran multicore on every new module** — a green ifx suite is not sufficient.
- **Do NOT use `-stdpar=gpu`** (managed-allocator double-free); use OpenMP `target` + `-gpu=mem:separate`, `MEDS_GPU=multicore|gpu`.
- **FPE-safe under `-fpe0`/`-Ktrap=fp`.** Clamp `w_flux_snow` to `swe/dt`; floor densities/depths (`max(tiny_depth, ...)`); **banded `t_3ple` comparisons** and the banded negligible-water guard in the inverter (never exact-float equality, never `swe==0`) — the last-flake path is the exposure.
- Reuse the adaptive step-doubling `*_SUBSTEP_FIXED` wrapper so snow, heat, and water share substep control.

---

## 11. How this differs from ED2 / FATES / CLM / ClimaLand

| Aspect | ED2 (sfcwater) | CLM5 | ClimaLand | **MEDS (this doc)** |
|---|---|---|---|---|
| Prognostic state | mass + energy; **J/kg persistent, J/m² in RK4** (dual unit) | layer temperature + ice/liq mass, two-step phase adjust | internal energy + water content | **mass + extensive internal energy (J/m²), single unit** |
| Phase/melt | `uextcm2tl`/`uint2tl` inverter, diagnosed each step | temperature-prognostic + explicit phase-change adjust | energy read-off | **`uext_to_temp` read-off (shared with soil), plateau by construction** |
| Integrator | explicit adaptive RK4 | operator-split BE | implicit | **operator-split BE substep; melt post-step read-off; OUT of ARK implicit set** |
| Layers (MVP) | up to `nzs≈8`, dynamic | fixed 5 (+more) | flexible | **single bulk layer P0; multi-layer P1** |
| Conduction | density-dependent Medvigy-2007 | density-dependent | — | **fixed `k_snow` P0 (throttled `g_base`); density-dependent P1** |
| Sublimation | implicit ice-referenced enthalpy | explicit `L_subl` | enthalpy | **implicit ice-referenced enthalpy twin (no separate `L_subl`)** |
| Solver reuse | bespoke TSW RHS | bespoke | — | **reuses `meds_soil_solver` Thomas + `uext_to_temp` + `ground_optics` unchanged** |
| Conservation | `sfcw_toler=1e-6` end check | budget | budget | **machine-precision `budget_accumulate` ledger, `error stop` in Debug** |

**One line:** MEDS snow = ED2's `sfcwater` physics with ClimaLand's single-unit internal-energy state, melt as a post-step read-off of the *same* inverter the soil freeze/thaw already uses, operator-split out of the ARK (like soil water and wood), delivered as a single-layer MVP that closes the winter bare-soil skin gap.

---

## 12. Bugs / quirks found in the ED2 reference

- **S1 — dual-unit `sfcwater_energy` (J/kg persistent vs J/m² in RK4).** ED2 converts at `rk4_copy_patch.f90:221,1301`, an error-prone seam. **MEDS fixes by construction:** extensive J/m² everywhere; no conversion.
- **S2 — precip always onto the top layer + a scalar virtual pool.** ED2 buffers sub-threshold precip in a separate `virtual_water/energy/depth` scalar. **MEDS carries deliberately (simplified):** P0 folds sub-threshold new snow into the soil-top store (no virtual state); the buffer is a P1 refinement (§4a).
- **S3 — thin-mass divide risk.** ED2's inverters guard with `tiny`/`water_stab_thresh`. **MEDS carries deliberately with a fix:** banded negligible-water guard in `uext_to_temp` (`abs(u_melt−u_freeze)<eps`) + banded `t_3ple` comparisons; last-flake dumps to soil before `swe` reaches machine-zero (§6, §10).
- **S4 — thin-layer forced equilibrium (`flag_sfcwater==1`).** ED2 forces a sub-`water_stab_thresh` layer into equilibrium with soil-top in *two* places (`rk4_misc.f90:1299` and `update_diagnostic_vars:228`). **MEDS carries deliberately (single authority):** one forced-equilibrium branch in the snow-mass kernel, no duplicate in a diagnostic pass (§6, avoids ED2's double-site drift).
- **S5 — negative sub-threshold mass repaid by condensing CAS vapor.** ED2's `rk4min_sfcw_moist=−5e-4` negative floor and vapor-condensation repayment (`rk4_misc.f90:1071`) is a numerical band-aid for adaptive RK4 overshoot. **MEDS fixes by construction:** the operator-split BE substep with `swe/dt` sublimation cap cannot drive `swe` negative; no negative floor needed.
- **S6 — committed-radiation redirection on mid-step TSW vanishing.** ED2 redirects "committed" SW to soil when TSW disappears mid-`dtlsm` (`rk4_derivs.f90:215`). **MEDS carries deliberately:** melt-out at the end of a substep dumps residual energy to `soil_energy(1)` conserving the full budget (§4c); no intra-step radiation redirection because RT is lagged one substep (the existing MEDS RT-join convention).
- **S7 — sublimation latent heat is implicit in the ice reference.** Not a bug but a subtle correctness point: ED2 gets latent-heat-of-sublimation *for free* by referencing the ice layer's own internal energy. **MEDS reproduces deliberately** — `internal_energy_ice` reference + `enthalpy_vapor` twin, no separate `L_sublimation` constant (§4d, avoids the common double-count-fusion bug).
- **S8 — snowfall energy uses `min(t3ple, T)`.** Frozen precip enthalpy must clamp temperature to `≤ t_3ple` (`ed_met_driver.f90:2560`). **MEDS carries:** `t_snowfall = min(t_3ple, tair_k)` in `internal_energy_ice` (§4a) — omitting the clamp would add spurious sensible heat with the snow mass.

---

## 13. Parameters & config (no hard-coded values)

Every snow parameter is **required from TOML**, presence-mapped, `error stop` if missing; only genuine universal constants go in `meds_constants`.

```toml
[snow]
fresh_snow_density      = 100.0    # [kg/m3]  fsdns; fresh-fall density (used by P1 compaction/albedo aging)
rho_snow                = 250.0    # [kg/m3]  P0 fixed SEASONAL-mean bulk density: snow_depth = swe/rho_snow
snow_albedo_vis_fresh   = 0.85     # [-]      fresh-snow VIS albedo
snow_albedo_nir_fresh   = 0.65     # [-]      fresh-snow NIR albedo
snow_albedo_vis_aged    = 0.518    # [-]      aged endpoint (ED2)
snow_albedo_nir_aged    = 0.435    # [-]      aged endpoint (ED2)
snow_emiss              = 0.970    # [-]      thermal emissivity
k_snow                  = 0.15     # [W/m/K]  P0 fixed conductivity (density-dependent in P1)
liquid_holding_frac     = 0.10     # [-]      percolation holding capacity (LEAF-3 1:9)
snow_stab_thresh        = 10.0     # [kg/m2]  forced-equilibrium below this (ED2 water_stab_thresh)
min_new_snow_mass       = 1.0e-3   # [kg/m2]  first-flake / tiny threshold (ED2 tiny_sfcwater_mass)
snow_rough              = 0.001    # [m]      roughness floor when buried (already on aero_cfg_t)
ny07_a                  = 2.5      # [-]      Niu-Yang07 snow-cover shape
ny07_m                  = 1.0      # [-]      Niu-Yang07 snow-cover exponent
snow_phase              = "on"     # SNOW_PHASE_ON (melt plateau always active)
percolation_scheme      = "leaf3"  # leaf3 | anderson76 (P1)
```

Universal constants already in / added to `meds_constants` (SI): `cp_ice=2093`, `cp_liq=4186`, `latent_heat_fusion=3.34e5`, `t_3ple=273.16`, `tsupercool_liq/vap`, `k_ice=2.29`, `k_water=0.57`, `k_air=0.025 [W/m/K]`. The one **new thermo primitive**: `internal_energy_ice(T) = cp_ice·T [J/kg]` in `meds_thermo` (the frozen twin of `internal_energy_liquid`).

**Wiring.** `meds_config_t` carries only the plain `real`/`integer` snow fields above; a biophysics `build_snow_params(cfg) -> snow_params_t` (called from the aux/init layer) assembles the typed record; TOML strings (`snow_phase`, `percolation_scheme`) mapped to `SNOW_*` enums by `meds_config_io`.

---

## 14. Coupling contracts

| Field | Units | Location (when snow present) | Bare-soil value today | Closed by |
|---|---|---|---|---|
| `t_snow_surf` | K | surface skin seen by CAS + radiation | `soil_temp(1)` | `uext_to_temp(snow_energy, swe, 0, …)` |
| `abs_sw_snow` / `abs_lw_snow` | W/m² | net rad on snow surface | `abs_sw_ground`/`abs_lw_ground` on soil | `apply_rt_forcing` with snow albedo/emiss (§4f) |
| `h_snow` / `le_snow` | W/m² | CAS enthalpy source `src_enth` | `h_ground`/`le_ground` | snow-surface balance (§4b) |
| `w_flux_snow` (sublimation) | kg/m²/s | CAS vapor source `src_vap`; snow mass sink | `soil_evap` | `snow_sublimation` (§4d), capped `swe/dt` |
| `g_base` | W/m² | soil top Neumann BC `eforc%g_top` | `g_top = Rn−H−LE` | snow-base conduction (§4g) |
| `melt` | kg/m²/s | `q_avail` infiltration input | (n/a) | free-liquid drain (§4e); internal transfer |
| `swe`, `snow_energy` | kg/m², J/m² | `whole_water`/`whole_energy` store totals | (n/a) | `budget_accumulate` (§9 tests 1–2) |
| `snowfac` | – | aero roughness blend | `0.0` (hard-zeroed) | Niu-Yang07 from SWE/depth (§4f, §6) |

---

## 15. Open questions

1. **Depth as prognostic vs derived in P0.** *Resolved (§3.1/§4g):* carry `snow_depth` as a prognostic geometric state, but in P0 derive it from mass at a **fixed seasonal-mean `rho_snow ≈ 250 kg/m³`** (not fresh density — avoids over-deepening/over-insulating an aged pack) with **no compaction**. Compaction/densification (depth decoupling from mass) is P1. *Open:* the exact P0 `rho_snow` value to use for Ithaca (tune against the January insulation test, §9 test 5).
2. **`snowfac`-ramped vs sub-patch fractional cover for P0.** *Resolved (§4f/§4g/§6):* P0 ramps the ground optics, the soil-top-BC blend, and the latent split by a single `snowfac` (Niu-Yang07) so first-flake/threshold/melt-out are all continuous — NOT a binary `nlayer≥1` switch. *Deferred to P2:* true sub-patch fractional radiation (separate snow-covered and bare radiative sub-columns, area-averaged); P0 carries `snowfac`-weighted optics on one ground surface. *Open:* whether the `snowfac(ny07_a, ny07_m)` shape needs site calibration.
3. **Sub-threshold new-snow buffer.** *Recommend:* P0 folds sub-`min_new_snow_mass` snowfall into the soil-top store (mass+ice enthalpy), deferring a scalar virtual-pool state to P1 — one fewer state to lockstep, no mass loss.
4. **Aging albedo curve in P0.** *Resolved (§4f):* MVP interpolates fresh↔aged endpoints by top-layer `snow_fliq` (wet snow → lower albedo, ED2-style), plus an optional simple time-since-snowfall age counter. A physical grain-size decay is P1.
5. **`k_snow` fixed vs density-dependent.** *Resolved (§4g):* fixed `k_snow` in P0 (the throttled `g_base` already delivers the qualitative insulation that caps the soil skin); the Medvigy-2007 density polynomial is P1 with compaction.
6. **Does the ARK path need snow in P0?** *Resolved (§6):* no — P0 delivers the split/Picard path only; the ARK path stays snow-free (guarded, like prognostic wood/leaf today). ARK-integration is a later gated phase.

---

## 16. References

- **ED2 temporary-surface-water (sfcwater) physics.** State arrays and dual J/kg↔J/m² unit convention: `memory/ed_state_vars.F90:1242-1258`, `dynamics/rk4_copy_patch.f90:219-223,1297-1302`. Rain/snow partition + `pcpg/qpcpg/dpcpg`: `driver/ed_met_driver.f90:2507-2562` (Jin et al. 1999). TSW energy/mass RHS + conduction + sensible/latent/sublimation: `dynamics/rk4_derivs.f90` (`leaftw_derivs`, `canopy_derivs_two`). Accumulation/percolation/runoff/layer management + thresholds + `snowfac`: `dynamics/rk4_misc.f90` (`adjust_sfcw_properties`, `update_diagnostic_vars`). Albedo + snow-cover + SW/LW absorption: `dynamics/radiate_driver.f90:626-790`. Freeze-thaw inverters + enthalpy: `utils/therm_lib8.f90` (`uint2tl8`, `uextcm2tl8`, `tl2uint8`, `tq2enthalpy8`). Constants + `thick` table: `init/ed_params.f90`, `init/ed_init.F90:659-678`. **The physics reference (§4, §12).**
- **MEDS energy seams.** `ground_surface_balance` (`meds_column_energy.f90:307`), `soil_energy_flux`/`soil_heat_be_step` top BC (`:44,201`), `canopy_air_update` (`:324`); freeze/thaw inverter `uext_to_temp`/`temp_to_uext` (`meds_thermo.f90:55,71`), `ENERGY_PHASE_ON` (`meds_biophysics_types.f90:204,290`). **The surface seam + shared inverter (§1.1, §4).**
- **MEDS hydrology seams.** `column_hydrology_flux` + `q_avail`/`w_surface`/`mass_resid` (`meds_column_hydrology.f90:75,137,193,214`), `ground_evaporation` (`:538`); `precip_phase` (`meds_forcing_kernels.f90:291`) and its dangling `snowf` (`meds_fast_loop.f90:482`). **Accumulation + melt→infiltration + sublimation sink (§4a,§4d,§4e).**
- **MEDS coupling / state / radiation.** Split `column_fast_step` ground balance (`meds_column_dynamics.f90:503-506,526`), whole-column ledgers (`:597-615`); ARK `surface_derivs`/`surface_frozen_t` (`meds_column_derivs.f90:60,256-258`); `ground_optics`/`apply_rt_forcing` snow albedo hook (`meds_fast_loop.f90:553-575`); `ggnet` snow-aware branch + `snowfac`/`snow_rough` (`meds_canopy_aerodynamics.f90:76-82`, `meds_biophysics_types.f90:314,363`); reservoir types + blends (`meds_column_state_types.f90:31-90`); `patch_index` reservoirs + 6 lockstep sites (`meds_demography_types.f90:130`, `meds_demography_fusefiss.f90`). **The wiring (§5, §6, §3).**
- **MEDS internal:** `docs/dev_plans/MEDS_ENERGY_BALANCE_DESIGN.md` (internal-energy state decision, the deferred `surface_layer(:)` note at L283, BE-Thomas reuse, negative-z geometry, budget discipline — the primary parent); `docs/dev_plans/MEDS_COLUMN_HYDROLOGY_DESIGN.md` (stateless-kernel + state/process wall, surface store + mass ledger, `snowf` deferral — the primary parent); `CLAUDE.md` (the wall, no-hard-coded-params, spell-out-acronyms naming, nvfortran issue #7, `-stdpar` prohibition). This doc **reconciles** the two parents' phase labels (energy-P2b, hydrology-P5) into a unified snow family with P0 = single-layer MVP (§7).