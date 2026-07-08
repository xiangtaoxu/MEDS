# MEDS Meteorological Forcing — Source & Wiring Design

The prerequisite that closes MEDS's largest coupling gap. Today MEDS has **no meteorological
forcing source**: the fast biophysics loop is OFF by default, and *every* atmospheric boundary
condition enters the fast kernels as a **constant, horizontally-uniform, in-type-default value**
hand-assembled by `fast_context_t` (`meds_fast_loop.f90`) and translated to the kernel inputs by two
private routines, `build_forcing` (→ `column_forcing_t`) and `fill_aenv` (→ `aero_env_t`). The one
time-derived environmental quantity in the whole model is the solar cosine `solar_cosz` in
`meds_time` — and it is not yet wired to a site latitude, because the TOML has no `[site]` block. The
CLAUDE.md "Reserved follow-ups" name this precisely: *"a meteorological forcing source"* is one of the
four missing pieces (with canopy RT, leaf energy balance, hydraulics) blocking real GPP.

This document designs (1) a per-site **forcing derived type** (`met_forcing_t`) and a small
**buffer/reader state** (`met_driver_t`), living at the polygon/site level — "a location sharing one
meteorological forcing" (CLAUDE.md §Architecture-1); (2) a stateless **reader module**
`src/io/meds_met_driver.f90` reading a NetCDF/ALMA-like file **and** a CSV path for the BCI tower
record, aligned to the model calendar via `meds_time`; (3) **per-variable temporal
interpolation/disaggregation** from the forcing interval (30 min) down to `dt_fast` (900 s), with a
full BCI-column → MEDS-field mapping table and unit conversions; (4) the **exact wiring** of forcing
into the fast loop (the `build_forcing`/`fill_aenv` seam) and, via daily aggregation, into the slow
loop; (5) a **Python prep script** converting `BCI_v5.1.csv` into a MEDS sub-hourly forcing NetCDF;
(6) a **phasing** P0→P1→P2 and a **test plan** driving the leaf / RT kernels offline to reproduce a
diurnal cycle.

The physics reference is **ED2**'s `ed_met_driver.f90` + `radiate_utils.f90` (HDF5 read →
`cgrid%metinput` → temporal interpolation + radiation breakdown + thermodynamic closure →
`cgrid%met` → per-site `cpoly%met`). The community interchange reference is the **ALMA** convention
(8 core drivers in SI, positive-downward flux sign) as realised by FLUXNET/PLUMBER2. As always: MEDS
ports ED2's *algorithms* (cosine-zenith SW disaggregation, Weiss–Norman partition, step-constant
precip) as **stateless `pure`/`elemental` kernels**, keeps the reader's mutable buffer in a
driver/IO-owned state type, and threads config as a passed derived type — never module globals.

> **DATA NOTE (read first).** The task referenced `BCI_flux/BCI_v6.1.csv`, but **only
> `/home/xiangtao/projects/BCI_flux/BCI_v5.1.csv` exists** (90,528 data rows, **half-hourly**,
> 2012-07-03 onward). This design targets **v5.1**; if a v6.1 with different columns/units later
> lands, only the Python prep script's column map (§7) and the CSV header parse (§4.3) change — the
> Fortran reader consumes the derived NetCDF/canonical record, not the raw tower columns. The v5.1
> header is: `date, tair, RH, vpd, p_kpa, PPT, Rs, Rs_dn, Rl_dn, Rl_up, Rnet, LE, H, NEE, Par_tot,
> Par_diff, SWC, ubar, ustar, WD, gpp, FLAG`.

---

## 1. Scope, target variables & reachability

### 1.1 What the forcing source must drive

**Primary target: the fast biophysics loop.** `run_fast_biophysics` (§`meds_fast_loop`) advances each
patch by `n_fast_per_slow = nint(dt_slow/dt_fast)` operator-split sweeps of `column_fast_step`. With
`dt_slow = 1 d` and `dt_fast = 900 s` that is **96 sub-steps spanning one day** — so the *diurnal
cycle itself lives inside the fast sub-step loop*, and the forcing must be re-evaluated **per
sub-step**, not held constant per slow day (the single most important consequence for the wiring,
§6). The kernels that consume forcing, and the field/type each reads, are (consolidated from the
consumption-seam brief):

| Forcing quantity | Consumed by (kernel) | Field / type | Units |
|---|---|---|---|
| **Air temperature** `Tair` | aerodynamics | `aero_env_t%theta_atm` | [K] air temp at `zref` (θ≈T near surface; no Poisson at P0, §6.2) |
| | CAS energy (`canopy_air_update`) | `enthalpy_atm` (folded via `cas_enthalpy_of_temp`) | [J/kg] |
| | (driver) | `fast_context_t%air_temp` | [K] |
| | phenology (slow) | `pheno_env_t%temp_day` | [K] daily mean |
| **Specific humidity** `Qair` | aerodynamics | `aero_env_t%shv_atm` | [kg/kg] |
| | CAS energy (humidity twin) | `enthalpy_atm` + `w_flux_ac` | [J/kg],[kg/m²/s] |
| | soil evap (`ground_evaporation`) | `chydro_forcing_t%q_air` (= CAS shv) | [kg/kg] |
| | (driver) | `column_forcing_t%shv_atm`, `fast_context_t%shv_atm` | [kg/kg] |
| **Pressure** `Psurf` | leaf gas exchange | `leaf_env_t%pressure` | [Pa] |
| | aerodynamics | `aero_env_t%press` | [Pa] |
| | veg/ground energy | `leaf_energy_env_t%press` | [Pa] |
| | (driver) | `fast_context_t%press` | [Pa] |
| **Precipitation** `Rainf` | soil hydrology | `chydro_forcing_t%precip_ground` | [kg/m²/s] post-interception |
| | per-cohort interception | `intercept_canopy_layer(rain_above)` | [kg/m²/s] |
| | (driver) | `column_forcing_t%precip`, `fast_context_t%precip` | [kg/m²/s] |
| **Shortwave** (direct+diffuse, PAR+NIR) | canopy radiation | `rad_forcing_t%incid_beam(band)`, `%incid_diff(band)` | [W/m²] canopy top |
| | | `rad_forcing_t%cosz` | [–] (derived, `solar_cosz`) |
| | leaf gas exchange | `leaf_env_t%par` (= abs_sw/lai·par_per_w) | [µmol/m²/s] |
| | veg energy | `leaf_energy_env_t%abs_sw` | [W/m²] |
| | (driver) | `column_forcing_t%abs_sw(:)`, `%abs_sw_ground`; `fast_context_t%rad_sw_top`, `%rad_sw_ground` | [W/m²] |
| **Longwave down** `LWdown` | canopy radiation | `rad_forcing_t%incid_diff(RAD_LW)` | [W/m²] |
| | veg energy | `leaf_energy_env_t%abs_lw` | [W/m²] net |
| | (driver) | `column_forcing_t%abs_lw(:)`, `%abs_lw_ground` | [W/m²] |
| **Wind** `Wind` | aerodynamics | `aero_env_t%u_ref`, `%zref` | [m/s], [m] |
| | (driver) | `fast_context_t%u_ref`, `%zref` | [m/s], [m] |
| **CO₂** `CO2air` | leaf gas exchange | `leaf_env_t%ca` (reads **prognostic** `bio%cas%can_co2`) | [µmol/mol] |
| | aerodynamics (`cstar`) | `aero_env_t%co2_atm` | [µmol/mol] free-atm |
| | CAS CO₂ (`column_co2_step`) | arg `co2_atm` | [µmol/mol] free-atm |
| | biogeochem fallback | `co2_opts_t%co2_atm_ref` | [µmol/mol] |
| | (driver) | `column_forcing_t%co2_atm`, `fast_context_t%co2_atm` | [µmol/mol] |
| **Air density** `ρ_air` (property, derived) | aero / CAS / leaf-energy / soil-evap / CO₂ | `%rho_air` (pervasive) | [kg/m³] |

**Two seam subtleties the design must respect (concerns from the brief):**

1. **Some "forcing-looking" fields are prognostic, not atmospheric.** `leaf_env_t%ca` reads the
   *prognostic* CAS `bio%cas%can_co2`; `leaf_env_t%leaf_temp` is the *prognostic* leaf-energy
   temperature; `ustar` is an **output** of aerodynamics (`aero_out_t%ustar`) consumed by the CAS
   energy and CO₂ kernels. The met source must fill only the **free-atmosphere reference** values
   (`co2_atm`, `theta_atm`, `enthalpy_atm`, `u_ref`) and let the prognostic CAS/leaf/`ustar`
   machinery evolve — it must **not** overwrite `can_co2`, `leaf_temp`, or supply `ustar`.
2. **Air temperature never enters the energy kernels as raw [K].** It is folded into
   `enthalpy_atm [J/kg]` (via `cas_enthalpy_of_temp(air_temp, shv)`, already in `meds_thermo`) for the
   CAS budget and into `theta_atm [K]` (potential temperature) for aerodynamics. The met source
   produces `air_temp` + `shv_atm`; the *existing* `build_forcing` already computes
   `enthalpy_atm = cas_enthalpy_of_temp(air_temp, shv_atm)`. So the seam is: **the met source writes
   the raw met scalars into `fast_context_t`'s met fields; `build_forcing`/`fill_aenv` keep their
   present logic unchanged.** This is the minimal-diff wiring (§6.2).

### 1.2 The canonical variable set (ALMA + ED2 union, what the reader produces)

The reader's job is to fill one instantaneous per-site record. The **canonical 12-field set** (ALMA
core + ED2's 4-stream SW split + derived `cosz`/`ρ_air`), with the interpolation policy each carries:

| Canonical field | Symbol | Units | Interp policy (§5) |
|---|---|---|---|
| `tair_k` | Tair | [K] | linear |
| `qair` | Qair | [kg/kg] specific humidity | linear |
| `psurf_pa` | Psurf | [Pa] | linear |
| `rainf` | Rainf | [kg/m²/s] | **step-constant** (conserve accumulation) |
| `wind` | Wind | [m/s] | **energy-form** (√ of variance-weighted squares, §5.3) |
| `lwdown` | LWdown | [W/m²] | linear |
| `par_beam`, `par_diffuse` | — | [W/m²] | **cosz-weighted** |
| `nir_beam`, `nir_diffuse` | — | [W/m²] | **cosz-weighted** |
| `co2` | CO2air | [µmol/mol] | linear (or config constant) |
| `cosz` | — | [–] | **recomputed** (`solar_cosz`), never stored |
| `rho_air` | ρ_air | [kg/m³] | **recomputed** from `Tair,Psurf,Qair` |

`SWdown = par_beam + par_diffuse + nir_beam + nir_diffuse` and
`rshort_diffuse = par_diffuse + nir_diffuse` recover the aggregate ED2 diagnostics. Snow (`Snowf`) is
folded into `Rainf` with an air-temperature phase split at ingest (§5.4); the tropical BCI record is
rain-only, so P0 treats all precip as liquid.

### 1.3 Reachability facts (from the ED2 & MEDS references)

- **ED2 stores forcing in derived types passed by argument**, never global scalars: raw time-buffered
  `met_driv_data` (14 arrays) → instantaneous `met_driv_state` (25 scalars) → per-site `cpoly%met`.
  MEDS mirrors this two-tier split (`met_driver_t` buffer → `met_forcing_t` instantaneous, §3).
- **ED2 disaggregates each SW stream by a mean daytime zenith factor** (`mean_daysecz`), forming a
  beam-perpendicular flux and rescaling to the instantaneous geometry; night (`cosz ≤ cosz_min`) is
  forced to zero. **MEDS deliberately does NOT port ED2's `⟨sec z⟩ = ⟨1/cosz⟩`** — that quantity does
  *not* conserve the interval-mean flux (by Jensen `⟨1/cosz⟩ ≥ 1/⟨cosz⟩`, badly so near
  sunrise/sunset). MEDS instead uses the **reciprocal of the window-mean cosine**, `1/⟨cosz⟩`, which is
  the unique factor that makes the reconstructed flux integrate back to `F_avg` (§5.1, the exact
  identity test 1d checks). MEDS already has `solar_cosz(t, t_sec, latitude_deg)` in `meds_time` — the
  analogue of `ed_zen`. It needs a **site latitude**, absent from the TOML today: adding `[site]`
  (§6.4) is a hard prerequisite for any radiation forcing.
- **ED2 interpolates precipitation step-constant** and wind as a squared (energy) quantity. MEDS
  ports both (§5.3–§5.4).
- **The `fast_context_t` comment is explicit**: *"the CALLER builds this from TOML in the production
  path; the MVP holds constant, horizontally-uniform boundary conditions."* Today only
  `test/test_fast_loop.f90` builds it. This design supplies the production caller.
- **CRITICAL WIRING GAP (unchanged by this doc, closed by §6):** `meds_main` line 124 calls
  `advance_one_step(site, cfg, is_new_month, is_new_year)` **without** the optional `fast_ctx`, and
  never calls `init_fast_reservoirs`, so `run_fast_biophysics` never fires in production even with
  `fast_biophysics_on = true`. Wiring the met source is done **together with** passing `fast_ctx` and
  seeding the reservoirs.

---

## 2. Where it lives (library DAG)

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`, with `biophysics`/
`biogeochemistry` **stateless siblings** linking `shared` only. Meteorological forcing has two halves,
placed on opposite sides of the state/process wall:

```
shared ─┬─ meds_forcing_types   (NEW, src/shared) — met_forcing_t + interp POLICY codes + forcing_config_t
        │     the pure DATA record + config; the DAG root can name it, like meds_time / meds_config
        ├─ meds_forcing_kernels  (NEW, src/shared) — PURE/ELEMENTAL disaggregation + SW-partition +
        │     rh_to_specific_humidity + precip_phase. Depends ONLY on meds_thermo + meds_time (both
        │     shared), so io can call it with NO new io→biophysics edge. GPU-safe leaf math.
        │     REUSES meds_thermo%air_density and meds_thermo%sat_vapor_pressure — does NOT re-invent them.
        │
        ↓
   io (reader) ── meds_met_driver  (NEW, src/io) — the READER + BUFFER: opens CSV/NetCDF, holds the
        prev/next records (met_driver_t, MUTABLE), advances the cursor to the model date, and produces
        an instantaneous met_forcing_t via the shared kernels. Owns file I/O — belongs in src/io beside
        meds_config_io / meds_netcdf_c. Called by meds_main / meds_stepper. NEVER named by a physics lib.
```

**Rationale for the split (matches the RT / energy / hydrology precedent), and the DAG edge it avoids:**

- **`met_forcing_t` (the instantaneous per-site record) is a PURE DATA type in `src/shared`
  (`meds_forcing_types.f90`)** — like `meds_time_t` and `meds_config_t`, the DAG root may define it,
  and any layer may read it. It holds **no allocatables tied to cohort count** (it is per-site
  scalars + 4 SW-band reals), so it is trivially copyable and GPU-mappable.
- **The disaggregation / partition math is STATELESS in `src/shared`
  (`meds_forcing_kernels.f90`), NOT in `src/biophysics`.** The reader in `src/io` must call the
  disaggregation / partition / RH kernels; if those lived in `libmeds_biophysics` it would introduce a
  **new `libmeds_io → libmeds_biophysics` link edge** (today `libmeds_io` links only `config_io` +
  `netcdf`, and `config_io` pulls in `libmeds_shared`). Because the forcing math is *pure stateless
  math depending only on `meds_thermo` + `meds_time` (both already in `shared`)*, MEDS places it in
  `src/shared` — `libmeds_io` already reaches `libmeds_shared` transitively (add the direct
  `target_link_libraries(meds_io PUBLIC meds_shared)` edge to make it explicit), so **no new
  physics-library edge is created and the acyclic DAG is preserved**. The `pure`/`elemental` kernels
  are `disaggregate_shortwave`, `partition_shortwave`, `interpolate_forcing`, `cosz_reconstruct_factor`,
  `rh_to_specific_humidity`, `precip_phase`, `synthesize_lwdown`. `solar_cosz` **stays in `meds_time`**
  (calendar-derived, already present).
- **No duplicated shared code.** `meds_thermo` (DAG root, reachable everywhere) *already* provides
  `elemental air_density(t_k, p_pa, shv)` (virtual-temperature form) and the Bolton-1980
  `sat_vapor_pressure` / `sat_specific_humidity`. `met_instant` calls `meds_thermo%air_density`
  directly, and `rh_to_specific_humidity` is built on `meds_thermo%sat_vapor_pressure` — **not** a
  fresh `esat`. Only `rh_to_specific_humidity` and the SW disaggregation/partition are genuinely new
  math.
- **The reader + mutable buffer lives in `src/io` (`meds_met_driver.f90`), inside `libmeds_io`** — it
  does file I/O (NetCDF-C bindings via the existing `meds_netcdf_c`, or a Fortran CSV read), so it
  belongs beside `meds_config_io`. It owns the *only* mutable forcing state (`met_driver_t`: the two
  bracketing records + timestamps + cursor + file handle + `forcing_config_t`). **No physics library
  ever names it** — the acyclic DAG is preserved (physics siblings link `shared` only).

**The controlling analogy:** `met_forcing_t` is to the fast loop what `rad_forcing_t` /
`chydro_forcing_t` / `cas_atm_forcing_t` are to their kernels — a read-only boundary-condition value
type; `met_driver_t` is the *reader's* private buffer, the analogue of ED2 `cgrid%metinput`, owned by
the driver and threaded (like `meds_io_t`) explicitly, never a global.

### 2.1 Files & CMake wiring

| File | Role | Analogue |
|---|---|---|
| `src/shared/meds_forcing_types.f90` (new) | `met_forcing_t` (instantaneous record), `forcing_config_t` ([forcing]/[site] block), `INTERP_*`/`SWPART_*`/`METAVG_*` selector codes | `meds_biophysics_types` type block |
| `src/shared/meds_forcing_kernels.f90` (new) | `pure`/`elemental` `disaggregate_shortwave`, `partition_shortwave`, `interpolate_forcing`, `cosz_reconstruct_factor`, `rh_to_specific_humidity`, `precip_phase`, `synthesize_lwdown`. **Calls `meds_thermo%air_density` + `%sat_vapor_pressure` — no re-invented `air_density`/`esat`.** Lives in `shared` (not biophysics) so `io` calls it with NO new io→biophysics edge (§2 rationale) | `meds_thermo` |
| `src/shared/meds_time.f90` (**extend — prerequisite, not incidental**) | expose/implement `seconds_into_day(now)`, `time_advance_seconds(now, dsec)` (second-accurate JDN + sec-of-day arithmetic across day boundaries), `seconds_between(a, b)`. Today only day/month arithmetic is public and `sec_of_day` is **private** (line 164); the 30 min→900 s interpolation and window slide cannot compile without these | its `time_advance_days`/`days_between` block |
| `src/io/meds_met_driver.f90` (new) | `met_driver_t` buffer + `met_open`, `met_advance`, `met_instant`, `met_close`; CSV + NetCDF backends | `meds_config_io` + `meds_io` |
| `src/io/meds_netcdf_c.f90` (**extend, P1**) | add `nc_inq_dimid` + `nc_inq_dimlen` ISO-C bindings — the met reader needs the length of the unlimited `time` dimension (`nrec`); today only `nc_open`/`nc_inq_varid`/`nc_get_vara_{int,double}` are bound, and `io_read_state` sidesteps this by reading fixed sizes from metadata | its existing bindings block |
| `src/io/meds_config_io.f90` (extend) | `[forcing]` + `[site]` key loaders (`req_*`) | its `[fast]` block |
| `src/shared/meds_config.f90` (extend) | `forcing_config_t` + `[site]` scalars on `meds_config_t` | its `[fast]` fields |
| `src/driver/meds_fast_loop.f90` (edit) | per-sub-step met refresh in `run_fast_biophysics`; `build_forcing`/`fill_aenv` read a `met_forcing_t` | — |
| `src/driver/meds_main.f90` (edit) | open the driver, seed reservoirs, thread `[site].reference_height` into `ctx%zref`, pass `fast_ctx` (closes the wiring gap) | — |
| `scripts/prep_bci_forcing.py` (new) | `BCI_v5.1.csv` → MEDS sub-hourly forcing NetCDF | — |
| `test/test_met_driver.f90` (new) | CTest: read, interpolate, disaggregate, diurnal cycle | `test_fast_loop` |

CMake globs `src/shared/*.f90`, `src/io/*.f90` with `CONFIGURE_DEPENDS`, so
the new modules auto-add; append the `test_met_driver` executable following the `test_fast_loop`
block. **Build nvfortran multicore on every new module** (a green ifx suite is not sufficient —
CLAUDE.md portability trap #7; the SW-partition kernels return arrays and must not be passed
straight into a call).

---

## 3. Derived types & state model

Two tiers (ED2's `met_driv_data` / `met_driv_state` split), plus the reader's config.

### 3.1 `met_forcing_t` — the instantaneous per-site record (`meds_forcing_types`, `src/shared`)

```fortran
! SW band indices ALIAS the existing RT vocabulary (meds_biophysics_types: RAD_VIS=1, RAD_NIR=2,
! RAD_LW=3) so the forcing bands and the two-stream RT bands share ONE set of names at the seam:
integer(ik), parameter :: RAD_PAR = RAD_VIS   ! = 1_ik ; PAR ≡ the RT "visible" band (do NOT mint a new code)

type :: met_forcing_t                                        ! per-SITE instantaneous atmospheric state
   real(wp) :: tair_k       = 288.0_wp    ! [K]        air temperature at reference height
   real(wp) :: qair         = 0.008_wp    ! [kg/kg]    specific humidity
   real(wp) :: psurf_pa     = 101325.0_wp ! [Pa]       surface pressure
   real(wp) :: rainf        = 0.0_wp       ! [kg/m2/s] liquid precipitation rate (post phase-split)
   real(wp) :: snowf        = 0.0_wp       ! [kg/m2/s] frozen precip (0 for tropical BCI)
   real(wp) :: wind         = 2.0_wp       ! [m/s]     wind speed at reference height
   real(wp) :: lwdown       = 380.0_wp     ! [W/m2]    downwelling longwave (positive down)
   ! reference-climate SW: the four streams SUM to 400 W/m2 to reproduce fast_context_t%rad_sw_top=400
   ! (NOT 0 — a zeroed default gives zero canopy light and BREAKS the CONST no-op; see §6.2):
   real(wp) :: par_beam     = 180.0_wp     ! [W/m2]    direct-beam PAR at canopy top
   real(wp) :: par_diffuse  = 40.0_wp      ! [W/m2]    diffuse PAR
   real(wp) :: nir_beam     = 150.0_wp     ! [W/m2]    direct-beam NIR
   real(wp) :: nir_diffuse  = 30.0_wp      ! [W/m2]    diffuse NIR   (Σ = 400 W/m2)
   real(wp) :: co2          = 400.0_wp     ! [umol/mol] free-atmosphere CO2
   real(wp) :: cosz         = 0.0_wp       ! [-]       cos(solar zenith); DERIVED each substep, not stored
   real(wp) :: rho_air      = 1.2_wp       ! [kg/m3]   DERIVED from tair/psurf/qair
contains
   ! swdown() = par_beam+par_diffuse+nir_beam+nir_diffuse ; rshort_diffuse() = par_diffuse+nir_diffuse
end type met_forcing_t
```

All defaults are valid (clean under `-fpe0`), so `met_forcing_t()` is a usable "reference climate" box:
its four SW streams **sum to 400 W/m² = today's `fast_context_t%rad_sw_top`** (`meds_fast_loop.f90:41`;
`air_temp=288`, `shv_atm=0.008`, `co2=400`, `press=101325`, `u_ref=2.0`, `rho_air=1.2` likewise match).
The **ground SW default is not carried here** — today's ad-hoc `rad_sw_ground=60` (`:42`) is an
*independent* constant with no canopy-top relation; under CONST it is reproduced by threading a
reference `rad_sw_ground` (below) rather than by Beer's law, so the reference climate is an honest
reproduction of the current constants. See §6.2 for why this makes CONST a **near-no-op** (canopy-top
SW exact) rather than the bit-for-bit claim the earlier draft asserted.

### 3.2 `met_driver_t` — the reader buffer (MUTABLE, `meds_met_driver`, `src/io`)

The **only** mutable forcing state. Holds the two bracketing records read from file plus the read
cursor; it is the ED2 `cgrid%metinput` analogue, owned by the driver and threaded like `meds_io_t`.

```fortran
type :: met_record_t                       ! one raw timestamped record as read (pre-interpolation)
   type(meds_time_t) :: when               ! timestamp of this record (interval label per METAVG_*)
   real(wp) :: tair_k, qair, psurf_pa, rainf, wind, lwdown, co2
   real(wp) :: par_beam, par_diffuse, nir_beam, nir_diffuse   ! already 4-stream-split at ingest/prep
end type met_record_t

type :: met_driver_t                        ! MUTABLE per-SITE reader state
   type(forcing_config_t) :: fcfg           ! [forcing]/[site] config (path, format, dt, lat/lon, CO2, ...)
   integer(ik) :: backend    = MET_BACKEND_CONST   ! CSV | NETCDF | CONST(no file)
   integer(ik) :: ncid       = -1_ik        ! NetCDF handle (via meds_netcdf_c), -1 if CSV/const
   integer(ik) :: unit       = -1_ik        ! Fortran unit for CSV
   integer(ik) :: nrec       = 0_ik         ! total records in file (for cycling)
   integer(ik) :: irec_prev  = 0_ik         ! cursor: index of the "previous" bracketing record
   real(wp)    :: dt_forcing  = 1800.0_wp   ! [s] native forcing interval (30 min for BCI)
   type(met_record_t) :: rec_prev, rec_next ! the two records bracketing the model time
   type(meds_time_t)  :: base_time          ! file record #1 timestamp (calendar anchor)
end type met_driver_t
```

### 3.3 `forcing_config_t` — the [forcing]/[site] block (`meds_forcing_types`, plain scalars)

Lives on `meds_config_t` (a plain-scalar component — `meds_config` is the DAG root and must not
`use` a biophysics type, exactly the `soil_thermal_params_t` rule from the energy design):

```fortran
type :: forcing_config_t
   logical            :: forcing_on   = .false.    ! master gate (independent of fast_biophysics_on)
   integer(ik)        :: format       = MET_BACKEND_CONST   ! "csv"|"netcdf"|"const"
   character(len=256) :: path         = ''         ! forcing file path
   real(wp)           :: dt_forcing   = 1800.0_wp  ! [s] native interval ("1800s")
   integer(ik)        :: avg_convention = METAVG_INSTANT    ! instant|end|begin|center (timestamp semantics)
   integer(ik)        :: sw_partition = SWPART_PASSTHROUGH  ! passthrough|weiss_norman|sib|clearidx
   integer(ik)        :: lwdown_source = LW_FILE            ! file|synthesize (Brutsaert/Idso)
   integer(ik)        :: gap_policy   = GAP_HOLD            ! hold|error|climatology
   real(wp)           :: co2_const    = 400.0_wp   ! [umol/mol] used when file lacks CO2 (BCI)
   real(wp)           :: rad_sw_ground_const = 60.0_wp ! [W/m2] CONST-backend ground SW (reproduces fast_context_t:42)
   logical            :: recycle      = .true.     ! cycle the record when the run outruns the file
   integer(ik)        :: start_clamp  = CLAMP_ERROR ! model start < base_time: error | hold record #1 (§4.2)
   ! ---- site geolocation ([site]) — prerequisites for solar geometry & lapse ----
   real(wp)           :: latitude_deg  = 9.15_wp   ! [deg +N]  BCI ~9.15 N
   real(wp)           :: longitude_deg = -79.85_wp ! [deg +E]  BCI ~-79.85 E (used for the solar-phase note, §5.1)
   real(wp)           :: utc_offset_h  = -5.0_wp   ! [h]       file-clock timezone vs UTC (BCI local std = UTC-5)
   logical            :: apply_solar_longitude = .false. ! P0: false (clock time AS-IF apparent solar; §5.1 note)
   real(wp)           :: elevation_m   = 120.0_wp  ! [m]       reserved for lapse (P2)
   real(wp)           :: reference_height = 40.0_wp! [m]       measurement height; must exceed max hgt_max
end type forcing_config_t
```

### 3.4 State-model summary

No global mutable forcing state. The mutable buffer `met_driver_t` is owned by the driver
(`meds_main` constructs it, threads it read-write into a refresh call each fast sub-step). The
instantaneous `met_forcing_t` is a read-only value produced by the reader and consumed by the fast
loop; the disaggregation/partition kernels are pure (no saved state); `cosz`/`ρ_air` are recomputed,
never stored. This keeps the physics libraries free of hidden met state and GPU-offload-safe — the
identical discipline as `rad_forcing_t` / `energy_forcing_t`.

---

## 4. The reader module — `src/io/meds_met_driver.f90`

Public API (mirrors `meds_io`'s open/write/close cadence and `ED2`'s `init/read/update` split):

```fortran
public :: met_driver_t, met_forcing_t
public :: met_open      ! (drv, fcfg)                 open file, read header, load records #1..2
public :: met_advance   ! (drv, now)                  slide the bracketing window to contain `now`
public :: met_instant   ! (drv, now, t_sec, lat) -> met_forcing_t   interpolate+disaggregate to `now`
public :: met_close     ! (drv)                       close handle / free
```

### 4.1 `met_open(drv, fcfg)` — header & first records

```
select case (fcfg%format)
case (MET_BACKEND_CONST)  ! no file: drv%rec_prev = rec_next = reference climate (met_forcing_t defaults)
case (MET_BACKEND_CSV)    ! open unit; skip header line; parse the BCI column order (§4.3);
                          !   read record #1 -> base_time, record #2 -> rec_next; drv%dt_forcing from cfg
case (MET_BACKEND_NETCDF) ! nc_open via meds_netcdf_c; read time[] + the canonical vars; nrec = dim(time);
                          !   base_time from time units attribute; load records #1..2
end select
drv%irec_prev = 1
```

The NetCDF backend reuses the existing `meds_netcdf_c` ISO-C bindings (`nc_open`, `nc_inq_varid`,
`nc_get_vara_double`) — no netCDF-Fortran dependency (unavailable for ifx/nvfortran here). **It must
also add two small bindings** (`nc_inq_dimid` + `nc_inq_dimlen`, §2.1): the reader needs `nrec` = the
length of the unlimited `time` dimension for cycling/EOF, and `meds_netcdf_c` binds no `nc_inq_dim*`
today (`io_read_state` sidesteps this by reading fixed sizes from metadata; the met reader cannot). The
canonical NetCDF layout (what the Python prep writes, §7) is one unlimited `time` dimension with the
12 canonical variables named per ALMA (`Tair`, `Qair`, `PSurf`, `Rainf`, `Wind`, `LWdown`,
`SWdown_par_beam`, …) plus a `time` coordinate carrying `units = "seconds since <base>"`.

### 4.2 `met_advance(drv, now)` — window slide (calendar alignment via `meds_time`)

The bracketing records must satisfy `rec_prev%when ≤ now < rec_next%when`. Advance is a
`transfer`-style slide (ED2 `transfer_ol_month` analogue), reading forward until the window contains
`now`. **Start-before-`base_time` (spin-up) is handled explicitly** — BCI begins 2012-07-03, and a
plausible spin-up config sets `start_time` earlier, so `now < rec_prev%when` (record #1) is a real
case the earlier draft left unspecified (it only handled running off the *end*):

```
if ( time_lt(now, drv%base_time) ) then           ! model time precedes the first forcing record
   select case (drv%fcfg%start_clamp)
   case (CLAMP_HOLD)                                ! clamp to record #1 (persistence), warn once
      call warn_once('met_driver: start_time < base_time; holding record #1 until forcing begins')
      return                                        ! rec_prev = rec_next = record #1; w_next forced to 0 (§4.4)
   case default                                     ! CLAMP_ERROR (default): fail loud, symmetric with EOF
      error stop 'met_driver: model start precedes first forcing record and start_clamp=error'
   end select
end if
do while ( time_lt(drv%rec_next%when, now) .or. time_eq(drv%rec_next%when, now) )
   drv%rec_prev = drv%rec_next
   drv%irec_prev = drv%irec_prev + 1
   if (drv%irec_prev + 1 > drv%nrec) then          ! ran off the end
      if (drv%fcfg%recycle) then                    ! wrap to the first record (spin-up cycling)
         call rewind_to_first(drv) ; cycle
      else
         error stop 'met_driver: model time past end of forcing and recycle=off'
      end if
   end if
   call read_record(drv, drv%irec_prev + 1, drv%rec_next)
end do
```

Calendar alignment uses `meds_time` throughout: record timestamps are `meds_time_t`; `time_lt`/
`time_eq` bound the window; leap years are exact. **Feb-29 reconciliation** (ED2's leap quirk) is a P1
concern only when cycling a non-leap forcing year onto a leap model year — P0/BCI reads the record's
own real dates so there is no mismatch.

### 4.3 CSV backend column map (BCI v5.1 → canonical record)

The CSV backend parses the **v5.1 header order** and converts to SI on read (the Python prep, §7, does
the same conversions to produce a NetCDF; the CSV path exists so BCI can be driven **without** a prep
step). Missing/NaN sentinels are handled by `gap_policy` (§5.5). See the full mapping table with unit
conversions in §5.2.

**One canonical saturation formula on BOTH paths.** `rh_to_specific_humidity(RH, T, P)` (the RH→q
conversion) is built on **`meds_thermo%sat_vapor_pressure`, the Bolton-1980 form**
`e_sat = 611.2·exp(17.67·t_c/(t_c+243.5))` (`meds_thermo.f90:29`) — the reader does **not** carry its
own `esat`. The Python prep (§7.2) **must use the identical Bolton formula**, not a Buck-1981 `esat`:
otherwise the same `(RH, T, P)` yields a *different* `qair` between the NetCDF (Python-built) and the
on-the-fly CSV path, breaking both the "file == on-the-fly agree" goal (§7.2) and the Python↔Fortran
agreement test (§9.5). Bolton is the **canonical formula**; §7.2 is written against it.

### 4.4 `met_instant(drv, now, t_sec, lat) -> met_forcing_t`

The hot path (ED2 `update_met_drivers` analogue). Computes the interpolation weight from the
bracketing timestamps, applies the per-variable policy (§5), recomputes `cosz` and `ρ_air`, and
returns the instantaneous record:

`seconds_between` / `time_advance_seconds` / `seconds_into_day` are the **new second-accurate
`meds_time` helpers** (§2.1 — prerequisites, not incidental), and `air_density` is `meds_thermo`'s
existing virtual-temperature routine, not a re-implementation:

```
dt_win = real(seconds_between(rec_prev%when, rec_next%when), wp)   ! = dt_forcing  (NEW meds_time helper)
elapsed = real(seconds_between(rec_prev%when, now), wp)
w_next  = clamp(elapsed / max(dt_win, tiny_num), 0, 1) ; w_prev = 1 - w_next
cosz    = solar_cosz(now, t_sec, lat)                              ! meds_time, floored >= 0

met%tair_k  = interpolate_forcing(INTERP_LINEAR, prev%tair_k, next%tair_k, w_next)
met%qair    = interpolate_forcing(INTERP_LINEAR, prev%qair,   next%qair,   w_next)
met%psurf_pa= interpolate_forcing(INTERP_LINEAR, prev%psurf_pa,next%psurf_pa,w_next)
met%wind    = interpolate_forcing(INTERP_LINEAR, prev%wind,   next%wind,   w_next)   ! §5.3 energy form
met%lwdown  = interpolate_forcing(INTERP_LINEAR, prev%lwdown, next%lwdown, w_next)
met%co2     = merge(fcfg%co2_const, interpolate_forcing(INTERP_LINEAR,prev%co2,next%co2,w_next), no_co2_in_file)
met%rainf   = interpolate_forcing(INTERP_STEP,   prev%rainf,  next%rainf,  w_next)   ! step-constant
call disaggregate_sw_bands(prev, next, cosz, cosz_bar, w_prev, w_next, met)          ! §5.1 (reciprocal-mean-cosz)
met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)   ! REUSES meds_thermo%air_density (shared)
```

---

## 5. Temporal interpolation & disaggregation (30 min → 900 s)

The native BCI interval is **1800 s (half-hourly)**; `dt_fast = 900 s`. Each 30-min record is split
into two 15-min fast sub-steps (plus finer if `dt_fast` shrinks). Per-variable policy (community
practice; concerns from the literature brief):

### 5.1 Shortwave — solar-zenith-weighted (the load-bearing method)

Linear interpolation of SW smears the diurnal cycle and is badly wrong in the tropics where `cosz`
changes fast even at 30-min resolution. MEDS uses a **strictly interval-mean-conserving cosine-zenith
disaggregation**. The requirement is exact: reconstruct an instantaneous flux `F(t)` from the
interval-mean `F_avg` over `[rec_prev%when, rec_next%when]` such that its time-average over the SAME
window returns `F_avg`. If `F(t) = f_perp · cosz(t)` on the daytime part of the window (0 at night),
then

```
⟨F⟩_win = f_perp · (1/T_win) ∫_daytime cosz dt = f_perp · ⟨cosz⟩_win   ⟹   f_perp = F_avg / ⟨cosz⟩_win
F(now)  = F_avg · cosz(now) / ⟨cosz⟩_win        ! 0 at night (cosz(now) <= cosz_min)
```

where **`⟨cosz⟩_win` is the window-mean of `cosz` with night sub-samples clamped to 0** — NOT the mean
secant. This is the load-bearing correction over the earlier draft: ED2's `mean_daysecz` computes the
**mean of the secant `⟨sec z⟩ = ⟨1/cosz⟩`**, and by Jensen's inequality `⟨1/cosz⟩ ≥ 1/⟨cosz⟩` — the two
disagree, badly, near sunrise/sunset where `cosz → 0`. Only the **reciprocal of the mean cosine**
`1/⟨cosz⟩` conserves the interval mean (the identity test 1d checks); `⟨sec z⟩` reconstructs an SW that
does *not* integrate back to `F_avg` and is **biased high at low sun**. MEDS therefore deliberately
does **not** port `mean_daysecz`; the kernel is named `cosz_reconstruct_factor` (returning
`1/⟨cosz⟩_win`), not `mean_daytime_secant`.

```
cosz_bar = cosz_reconstruct_factor(lat, win_start, dt_sub, dt_win)   ! = 1 / ⟨cosz⟩_win ; ⟨⟩ over ALL subsamples
F(now)   = merge(F_avg * cosz(now) * cosz_bar, 0.0_wp, cosz_bar_valid .and. cosz(now) > cosz_min)
```

applied **independently to each of the four SW streams** (`par_beam, par_diffuse, nir_beam,
nir_diffuse`). `cosz_reconstruct_factor` subsamples `solar_cosz` over the window at `dt_sub` (midpoint
rule), forming `⟨cosz⟩_win = (1/N) Σ_i max(cosz_i, 0)` — night sub-samples contribute exactly 0, so the
normaliser is the mean over the **full window** (this is what makes the average close over a
sunrise/sunset window, where daytime seconds < window seconds). **The divide-by-small guard is on the
aggregate `⟨cosz⟩_win`, not per-step:** `cosz_bar_valid = ⟨cosz⟩_win > cosz_bar_min` (a fully-night
window has `⟨cosz⟩_win = 0` and `F_avg = 0`, so all SW routes to 0); there is no per-step secant to blow
up. **`avg_convention` (METAVG_*)** shifts the window reference (instant / end / begin / center) —
FLUXNET half-hourly is interval-ending-labelled, reanalysis differs; getting this wrong phase-shifts
the diurnal cycle by up to one window, on top of the local-time solar-phase offset (below).

**Solar-phase / timezone assumption at P0 (make it explicit, it is NOT "apparent solar time").**
`solar_cosz(t, t_sec, latitude_deg)` (`meds_time.f90:79`) treats `t_sec` as **local apparent solar
seconds** — `frac_day = t_sec/86400`, so `t_sec` at 0.5·86400 is solar noon — with **no longitude and
no equation-of-time (EoT) term**. But BCI timestamps are **local *clock* (standard) time** (UTC−5,
central meridian −75 °E; the site at −79.85 °E lies ~19 min west of that meridian), and EoT adds up to
±16 min. Feeding clock-seconds straight into `solar_cosz` therefore phase-shifts the modeled diurnal
SW by roughly **20–35 min**. The earlier draft's "longitude reserved; local apparent solar time only in
P0" is **backwards** — P0 in fact uses clock time *as if* it were apparent solar time. Two honest
options, both documented rather than hidden:
- **P0 default (`apply_solar_longitude=.false.`):** run in **local standard time**, accept the fixed
  ~20–35 min diurnal-phase error, and bound it in the test plan (test 2 cross-checks day-length /
  noon-time against the astronomical value and asserts the offset stays under one window).
- **Optional now (a few lines):** set `apply_solar_longitude=.true.` and pass
  `t_sec_solar = t_sec + (longitude_deg − 15·utc_offset_h)·240 s/deg + eot(doy)` into `solar_cosz`,
  removing the standing offset. Recommended once the longitude term is validated; it is cheap and makes
  P0 truly apparent-solar.

If the file supplies only *total* SWdown (no PAR/diffuse split — e.g. a reanalysis path), the
`sw_partition` kernel first splits it (§5.6) before disaggregation. BCI supplies `Par_tot` +
`Par_diff` directly (§5.2), so P0 needs only the NIR complement + diffuse-fraction transfer, not a
full partition.

### 5.2 BCI-column → MEDS-field mapping table (with unit conversions)

| BCI column | Units (raw) | → MEDS canonical field | Conversion | Policy / notes |
|---|---|---|---|---|
| `date` | `YYYY-MM-DD HH:MM` local | `met_record_t%when` | parse via `time_from_string` | interval label = **end** of 30-min avg (assume; set `avg_convention=end`) |
| `tair` | °C | `tair_k` | `+ 273.15` | linear |
| `RH` | % | `qair` | `rh_to_specific_humidity(RH/100, tair_k, psurf_pa)` | linear (needs `tair`,`p_kpa`) |
| `vpd` | kPa | (check) | — | not ingested; RH+T is authoritative. Cross-check only |
| `p_kpa` | kPa | `psurf_pa` | `× 1000` | linear |
| `PPT` | mm / 30 min | `rainf` | `× (rho_h2o/1000) / 1800 = PPT/1800` kg/m²/s | **step-constant**; all liquid (tropical) |
| `Rs` | W/m² (incoming SW) | `SWdown` total | as-is (clamp ≥ 0; night ≈ −0.4 → 0) | source for NIR complement |
| `Rs_dn` | W/m² | (alt SWdown) | — | **many NaN** early rows → not primary; `Rs` used |
| `Rl_dn` | W/m² (LW down) | `lwdown` | as-is | **NaN gap-fill / synthesize** (§5.5, §5.7) |
| `Rl_up` | W/m² (LW up) | — | — | not forcing (diagnostic) |
| `Rnet, LE, H, NEE, gpp` | W/m², µmol | — | — | **validation targets**, not forcing (§9) |
| `Par_tot` | µmol/m²/s | `par_beam+par_diffuse` (W) | `× 1/4.6` µmol→W (PAR); NaN at night→0 | split by `Par_diff` |
| `Par_diff` | µmol/m²/s | `par_diffuse` (W) | `× 1/4.6`; `par_beam = (Par_tot−Par_diff)/4.6` | cosz-weighted |
| `SWC` | m³/m³ | — | — | soil init only (reserved; not atmospheric forcing) |
| `ubar` | m/s | `wind` | as-is (floor at `u_min`) | linear (energy form) |
| `ustar` | m/s | (optional) | as-is | **not ingested** — `ustar` is an aero *output* (§1.1); reserved |
| `WD` | deg | — | — | wind direction, unused (single site) |
| `FLAG` | — | QC | — | gate gap-fill; `FLAG≠0` → treat suspect (§5.5) |

**NIR derivation.** BCI measures PAR (`Par_tot`) and total SW (`Rs`) but not NIR directly. Given
`PAR_W = Par_tot/4.6` and `SWdown = Rs`: `NIR_W = max(0, Rs − PAR_W)`. The diffuse *fraction* is
transferred from PAR (`fdiff = Par_diff/max(Par_tot, ε)`), giving `nir_diffuse = fdiff·NIR_W`,
`nir_beam = (1−fdiff)·NIR_W`. `par_beam = (Par_tot−Par_diff)/4.6`, `par_diffuse = Par_diff/4.6`. This
is done in the **Python prep** (§7) so the NetCDF carries the 4 streams; the CSV backend does the same
inline. (`1/4.6 ≈ 0.2174 W per µmol PAR`; the 4.6 µmol/J PAR conversion is standard.)

**Transferring PAR's diffuse fraction to NIR is a P0 approximation, and a known biased one.** Shorter
(PAR) wavelengths scatter more than NIR, so the *true* NIR diffuse fraction is **lower** than PAR's;
setting `fdiff_NIR = fdiff_PAR` therefore **overestimates `nir_diffuse` and underestimates `nir_beam`**.
This is acceptable at P0 (the two NIR streams still sum to the correct `NIR_W`, and P0 RT is a coarse
LAI-share split anyway), but it is exactly the reason the P1 **Weiss–Norman** partition (§5.6) computes
PAR and NIR diffuse fractions *separately* per band. P1 either replaces the transfer with WN85's
band-specific fractions or applies a simple NIR-diffuse reduction factor (`fdiff_NIR ≈ c·fdiff_PAR`,
`c < 1`).

**CO₂** — BCI has **no CO₂ column**; `co2 = fcfg%co2_const` (default 400 µmol/mol), the single
free-atmosphere authority that fans out to `aero_env_t%co2_atm`, `column_forcing_t%co2_atm`, and the
CAS-CO₂ reference (the brief flags this value is duplicated across ~5 types — one config authority
writes it once).

### 5.3 Wind — energy-conserving interpolation

Wind (and `ustar` if ever ingested) interpolate as **squared** quantities then square-rooted (ED2's
`vels` accumulates `u²`), so the kinetic-energy content is conserved across record boundaries:
`wind(now) = sqrt(w_prev·wind_prev² + w_next·wind_next²)`, floored at a `u_min` (numerical stability
of the M–O similarity in aerodynamics).

### 5.4 Precipitation — step-constant, phase-split at ingest

`Rainf` is held **step-constant** across the interval (never linearly interpolated — interpolation
smears intense events into physically wrong drizzle and breaks infiltration/runoff). Phase split
(`precip_phase`, ED2 Jin-1999 bands on `Tair`) partitions total precip into `rainf`/`snowf`; BCI is
tropical (all liquid), so P0 sets `snowf = 0`. Precip enters the fast loop as
`column_forcing_t%precip` (ground-reaching, post-interception) and `intercept_canopy_layer`'s
`rain_above`.

### 5.5 Gap policy (NaN handling)

BCI has NaN in `Rs_dn`, `Rl_dn`, `Par_tot/diff`, etc. in early rows. `gap_policy`:

- **`GAP_HOLD` (default P0):** carry the last valid value forward (persistence) for short gaps; for a
  gap at record #1, use the reference-climate default. Cheapest, adequate for P0.
- **`GAP_CLIMATOLOGY` (P1):** fill from a diurnal-mean climatology computed by the Python prep (a
  by-hour-of-day mean of valid values), the MDS-lite analogue.
- **`GAP_ERROR`:** hard stop on any NaN in a required field (strict validation runs).

The **Python prep does the heavy gap-fill** (§7), so the NetCDF the Fortran reader consumes is
gap-free; `gap_policy` in Fortran is the last-line guard for the raw-CSV path.

### 5.6 SW partition seam (`partition_shortwave`, for reanalysis / total-SW inputs)

A pluggable dispatch (ED2 `imetrad` analogue) over total SWdown + `psurf` + `cosz` → 4 streams:
`SWPART_PASSTHROUGH` (use the file's own split — the **BCI/P0 default**, since PAR beam/diffuse are
measured), `SWPART_WEISS_NORMAN` (WN85, needs pressure), `SWPART_SIB` (Sellers-86), `SWPART_CLEARIDX`
(Boland/Tsubo clearness index). Ported as `pure` kernels from `radiate_utils.f90`. P1 wires
`weiss_norman`; P2 exposes the full selector for gridded reanalysis that ships only total SW.

### 5.7 LWdown source / synthesis (`synthesize_lwdown`, `lwdown_source`)

When the file lacks LWdown or it is NaN (`LW_SYNTHESIZE`), reconstruct from `Tair`+`Qair` with a
clear-sky emissivity (Brutsaert 1975 / Idso-Jackson 1969) plus an optional cloud correction from the
SW clearness index: `LWdown = ε_clear·σ·Tair⁴·(1 + a·(1−kt))`. BCI *has* `Rl_dn` (with gaps), so P0
uses `LW_FILE` + gap-fill; `LW_SYNTHESIZE` is the fallback for records/datasets without longwave.

---

## 6. Wiring — how forcing enters the loops

### 6.1 The injection points (exact lines)

Two private routines in `meds_fast_loop` translate the constant `fast_context_t` into the kernel
inputs, and are the **exact seam** a real forcing source replaces:

- **`build_forcing(forc, coh, ctx, sum_lai)`** → `column_forcing_t`: today reads `ctx%air_temp`,
  `ctx%shv_atm`, `ctx%co2_atm`, `ctx%rad_sw_ground`, `ctx%precip`, and LAI-splits `ctx%rad_sw_top`
  across cohorts. It already computes `enthalpy_atm = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)`.
- **`fill_aenv(aenv, bio, ctx)`** → `aero_env_t`: reads `ctx%u_ref/zref/press/rho_air/air_temp/
  shv_atm/co2_atm` and the live CAS/ground state.

### 6.2 Minimal-diff wiring: refresh `fast_context_t`'s met fields per sub-step

The cleanest change that preserves `build_forcing`/`fill_aenv` verbatim: **make `fast_context_t`'s met
scalars time-varying**, refreshed from `met_instant` at the top of *each fast sub-step*. Because the
diurnal cycle lives across the `n_fast_per_slow` sub-steps (§1.1), the refresh must be **inside** the
`isub` loop of `run_fast_biophysics`, and the loop must know the sub-step wall-clock time.

Concretely, `run_fast_biophysics` gains the current date and the driver buffer, and the sub-step loop
rebuilds forcing each sweep:

```fortran
subroutine run_fast_biophysics(site, ctx, cfg, drv, now, worst_energy, worst_water, n_budget_fail)
   type(fast_context_t), intent(inout) :: ctx        ! met scalars refreshed each substep
   type(met_driver_t),   intent(inout) :: drv        ! reader buffer (mutable)
   type(meds_time_t),    intent(in)    :: now        ! start-of-slow-step date
   ...
   do isub = 1_ik, cfg%n_fast_per_slow
      t_sec  = seconds_into_day(now) + real(isub-1, wp)*cfg%dt_fast + 0.5_wp*cfg%dt_fast   ! substep midpoint
      t_sub  = time_advance_seconds(now, (isub-1)*nint(cfg%dt_fast))
      call met_advance(drv, t_sub)
      met = met_instant(drv, t_sub, t_sec, cfg%forcing%latitude_deg)     ! §4.4
      call apply_met_to_ctx(ctx, met)                                    ! write met scalars into ctx
      call build_forcing(forc, coh, ctx, sum_lai)                        ! UNCHANGED body
      call fill_aenv(aenv, bio, ctx)                                     ! UNCHANGED body
      call column_fast_step(cfg%dt_fast, cfg, ctx%ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      ... gpp_accum integration (unchanged) ...
   end do
```

`apply_met_to_ctx(ctx, met)` sets `ctx%air_temp = met%tair_k`, `ctx%shv_atm = met%qair`,
`ctx%press = met%psurf_pa`, `ctx%rho_air = met%rho_air`, `ctx%co2_atm = met%co2`,
`ctx%u_ref = met%wind`, `ctx%precip = met%rainf`, and `ctx%rad_sw_top = met%swdown()`. Four points the
earlier draft glossed:

- **Ground SW.** Today's `rad_sw_ground = 60` (`meds_fast_loop.f90:42`) is an *independent* ad-hoc
  constant, not a canopy-top derivative. To keep CONST a faithful reproduction, `apply_met_to_ctx`
  sets `ctx%rad_sw_ground = fcfg%rad_sw_ground_const` under the **CONST** backend (default 60), and
  `= met%swdown()·exp(−k·LAI)` (Beer's law) under CSV/NetCDF. At the P1 RT join (§6.3) this is
  replaced by the two-stream `abs_sw_ground`.
- **`theta_atm` is fed raw `T`, and the type label was misleading.** `fill_aenv` sets
  `aero_env_t%theta_atm = ctx%air_temp` (`meds_fast_loop.f90:220`) with **no Poisson
  `(p0/p)^(Rd/cp)` conversion**, so the "potential temp at zref" label in the §1.1 table is aspirational.
  MEDS uses the near-surface approximation **`θ ≈ T`** at P0 (the reference height is ~40 m, the factor
  is ≈1); the table label is corrected to "air temp at zref (θ≈T near surface)". Adding the true Poisson
  factor is a one-line P1 refinement in `fill_aenv`, not a forcing-reader concern.
- **`zref`.** `ctx%zref` defaults to **30 m** (`meds_fast_loop.f90:35`), but the configured BCI
  `[site].reference_height` is **40 m**. `meds_main` threads `cfg%site%reference_height → ctx%zref`
  when it builds the static context (§6.5) so the aerodynamics see the true measurement height; this is
  a build-time wiring line, not part of the per-substep refresh.
- **`lwdown` is carried but DEAD at P0.** `met%lwdown` (and `met%par_beam/par_diffuse/nir_*`,
  `met%cosz`) are handed to the RT seam **only once §6.3 lands**. At P0, `build_forcing` sets only
  `abs_sw`/`abs_sw_ground`; `column_forcing_t%abs_lw` is left unset until RT arrives, so **`lwdown`
  cannot affect any P0 result** and P0 tests must not assert on it.

Because the four SW streams sum to 400 and CONST holds them temporally flat (no cosz disaggregation —
a constant climate is diurnally flat by construction), the **canopy-top SW reproduces today's
`rad_sw_top = 400` exactly**, and air T / humidity / CO₂ / wind / pressure / ρ_air match the current
constants. This makes CONST a **near-no-op refactor** — but **not** the "bit-for-bit" claim the earlier
draft made: `cosz`/`rho_air` are recomputed each substep, and ground SW now flows through
`rad_sw_ground_const`. Test 4 (§9) is re-specified accordingly to compare against a documented CONST
baseline (canopy-top SW exact; budgets within round-off), **not** a byte-identical replay.

> **Design note — cleaner alternative (P1):** thread a `met_forcing_t` argument straight into
> `build_forcing`/`fill_aenv` (and later `column_forcing_t`), retiring the met scalars on
> `fast_context_t` entirely so `fast_context_t` carries only the *static* `column_config_t` +
> initial soil state. This removes the `apply_met_to_ctx` shim and the co2/temp duplication. P0 keeps
> the shim for minimal churn; P1 does the type-level cleanup once the RT-driven per-cohort SW lands.

### 6.3 Per-cohort shortwave — RT replaces the LAI-share split

Today `build_forcing` splits `rad_sw_top` across cohorts by **LAI share** (a placeholder). The real
per-cohort absorbed SW/PAR comes from **canopy RT** (`meds_canopy_radiation`), whose `rad_forcing_t`
consumes exactly `incid_beam(band)`, `incid_diff(band)`, `cosz` — which `met_forcing_t` now supplies
(`par_beam/par_diffuse/nir_beam/nir_diffuse` → the two bands, `cosz` from `solar_cosz`, `lwdown` →
`incid_diff(RAD_LW)`). The wiring: `run_fast_biophysics` calls `canopy_radiation(opt, rad_forcing_from
(met), coh, …)` to fill `column_forcing_t%abs_sw(:)`/`abs_lw(:)` per cohort and `abs_sw_ground`/
`abs_lw_ground`, **replacing** the LAI-share split. This is the RT↔forcing join; it is a P1 item
(P0 keeps the LAI split driven by the *time-varying* `rad_sw_top`, which already yields a diurnal GPP
cycle).

### 6.4 Slow-loop consumption — daily aggregation

The slow loop (`meds_vegetation_dynamics`) consumes forcing two ways:

1. **Via the fast→slow GPP handoff (primary).** `carbon_growth` reads `a_carbon =
   cohort%gpp_accum(j)` when `fast_biophysics_on`. `run_fast_biophysics` already integrates
   sub-daily GPP onto `gpp_accum`; with real time-varying SW forcing this becomes a **physically
   meaningful daily GPP** (today's constant SW gives a flat, unphysical value). No new slow-loop
   plumbing — the met source improves the *content* of `gpp_accum`, not the seam.
2. **Directly, for slow-loop drivers that need daily aggregates.** `pheno_env_t%temp_day` (daily-mean
   air T), `%avail_water`, `%daylength`, `%doy` are daily/seasonal cues. A small
   **daily accumulator** on the driver (sum `met%tair_k`, count, and `∫SWdown` over the day's fast
   sub-steps; `daylength` from `cosz > cosz_min` count; `doy = day_of_year(now)`) produces a
   `met_daily_t` the slow loop reads at day roll-over. P0 defers phenology forcing (the phenology
   module is a directional signal, already testable standalone); the accumulator is a P1 add when
   phenology is wired into the stepper.

### 6.5 Closing the production wiring gap (`meds_main`)

`meds_main` must (a) build a `met_driver_t` from `[forcing]`, (b) build a `fast_context_t` (static
`column_config_t` + initial soil from `[fast]`/`[site]`), **threading `cfg%site%reference_height` into
`ctx%zref`** so the aerodynamics use the configured 40 m rather than the 30 m field default
(`meds_fast_loop.f90:35`), (c) call `init_fast_reservoirs(site, ctx)` once, and (d) pass `fast_ctx`
**and** the driver/date into `advance_one_step`. Since
`advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx)` already has the optional `fast_ctx`
hook (and gates on `cfg%fast_biophysics_on .and. present(fast_ctx)`), it grows two optional args
`met_drv`/`now` threaded to `run_fast_biophysics`:

```fortran
! meds_main, replacing the line-124 call:
if (cfg%fast_biophysics_on) then
   call advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx=ctx, met_drv=drv, now=prev)
else
   call advance_one_step(site, cfg, is_new_month, is_new_year)
end if
```

**Latent-coupling hazard the wiring must avoid (brief concern):** if `fast_biophysics_on = true` but
`fast_ctx` is absent (today's `meds_main`), `carbon_growth` reads a `gpp_accum` that
`run_fast_biophysics` never filled/reset — stale/zero GPP. The fix is exactly this always-pass-when-on
wiring; `validate_config` additionally asserts `forcing_on` when `fast_biophysics_on` and
`growth_source == GS_CARBON` (no silent stub-GPP carbon run).

### 6.6 The `[forcing]` and `[site]` config blocks

New TOML sections, loaded by `meds_config_io` with the existing `req_*` helpers (`req_l`, `req_dur`
for durations, `req_r`, `req_i`, and a new `req_str` for paths / a mapped `req_enum`); every key
presence-mapped and hard-error on absence (the house rule — no defaults in the loader):

```toml
[site]
latitude         = 9.15        # [deg +N]  BCI (Barro Colorado Island, Panama)
longitude        = -79.85      # [deg +E]  used for the solar-longitude/EoT correction (§5.1)
utc_offset       = -5.0        # [h]       file clock timezone vs UTC (BCI local standard = UTC-5)
elevation        = 120.0       # [m]       reserved for lapse (P2)
reference_height = 40.0        # [m]       measurement height; MUST exceed every PFT hgt_max (asserted)

[forcing]
forcing_on       = true        # master gate; independent of [fast].fast_biophysics_on
format           = "csv"       # "csv" (BCI) | "netcdf" (ALMA/PLUMBER2) | "const" (reference climate)
path             = "BCI_flux/BCI_v5.1.csv"
timestep         = "1800s"     # native forcing interval (30 min)
avg_convention   = "end"       # timestamp semantics: "instant"|"end"|"begin"|"center" (BCI = interval-ending)
apply_solar_lon  = false       # P0: false = local standard time (clock-as-solar; ~20-35 min phase err, §5.1)
sw_partition     = "passthrough" # "passthrough"(BCI: measured PAR) | "weiss_norman" | "sib" | "clearidx"
lwdown_source    = "file"      # "file" (BCI Rl_dn, gap-filled) | "synthesize" (Brutsaert)
gap_policy       = "hold"      # "hold" | "climatology" | "error"
co2_const        = 400.0       # [umol/mol] free-atm CO2 (BCI has no CO2 column)
rad_sw_ground    = 60.0        # [W/m2] CONST-backend ground SW (reproduces today's constant)
start_clamp      = "error"     # model start < first record: "error" | "hold" (§4.2)
recycle          = true        # cycle the record when the run outruns the file (spin-up)
```

`req_dur("forcing.timestep")` reuses the duration parser; `format`/`avg_convention`/`sw_partition`/
`lwdown_source`/`gap_policy`/`start_clamp` map strings → the `MET_*`/`METAVG_*`/`SWPART_*`/`LW_*`/
`GAP_*`/`CLAMP_*` codes via small mappers (the `req_scheme` pattern). `validate_config` asserts:
`reference_height > max(hgt_max)` (ED2's abort rule), `latitude ∈ [−90,90]`, `dt_fast ≤ timestep`,
**`mod(timestep, dt_fast) == 0`** (the forcing interval must be an integer multiple of `dt_fast` for
clean sub-step windowing — 1800/900 = 2 for BCI), and — for the CSV/NetCDF backends — that
`start_time`/`end_time` lie **within the forcing record grid** (else the `start_clamp`/`recycle` paths
of §4.2 govern, and the assertion documents which).

**UTC vs local-standard-time policy (a P2 seam, stated now).** BCI timestamps are **local clock time**
(`utc_offset = −5`); P2 reanalysis products (CRUNCEP/GSWP3/WFDE5/ERA5) are **UTC**. `solar_cosz` needs
**local apparent solar** seconds. The single conversion seam is
`t_sec_solar = sec_of_day(now) − utc_offset·3600 + (longitude·240 + eot(doy))`, gated by
`apply_solar_lon`. At P0 (BCI, `apply_solar_lon=false`) the offset is left in as the documented
~20–35 min phase error (§5.1); at P2 the same one-line seam converts UTC forcing to solar time with no
new machinery. The earlier draft named none of this ("longitude reserved").

---

### 6.7 Forcing provenance — echo the assembled record to the diagnostic stream

Debugging the diurnal reconstruction (the cosz disaggregation, the avg-convention phase, the solar-time
offset) is nearly impossible from GPP alone. So the assembled instantaneous `met_forcing_t` is
**echoed to the diagnostic netCDF** (`meds_io`'s `-D-output.nc` stream) at the diagnostic cadence: the
13 scalar fields (`tair_k`, `qair`, `psurf_pa`, `rainf`, `wind`, `lwdown`, the four SW streams, `co2`,
`cosz`, `rho_air`) as site-level timeseries variables, plus `swdown()`/`rshort_diffuse()` diagnostics.
This is a small, additive `io_write_snapshot` extension (site-scalar variables, no new dimension) and
is the single most useful artefact for validating that the reconstructed SW tracks the measured `Rs`
and peaks at the right local time. It is **write-only provenance** — nothing consumes it — so it is
safe to add at P0 and costs one row of scalars per output record.

## 7. Python prep script — `scripts/prep_bci_forcing.py`

Converts `BCI_v5.1.csv` → a MEDS sub-hourly forcing NetCDF (`BCI_v5.1_forcing.nc`) the NetCDF backend
reads directly. Belongs in `scripts/` (beside the install helpers) or `python/` (the `meds` package);
a standalone script keeps it dependency-light (`pandas` + `numpy` + `netCDF4`). Rationale for a prep
step: the heavy gap-fill, NIR derivation, unit conversion, and optional 30→15 min disaggregation are
one-time, clearer in Python, and produce a clean ALMA-named NetCDF that any MEDS run (and other
models) can share.

### 7.1 Column table (raw → NetCDF variable)

| NetCDF var | Units | Built from | Conversion |
|---|---|---|---|
| `time` | s since `<base>` | `date` | `time_from_string`; seconds since row 0 |
| `Tair` | K | `tair` | `+ 273.15` |
| `Qair` | kg/kg | `RH`, `tair`, `p_kpa` | `rh_to_q(RH/100, Tair, PSurf)` via **Bolton `esat`** (= `meds_thermo`, §4.3) |
| `PSurf` | Pa | `p_kpa` | `× 1000` |
| `Rainf` | kg/m²/s | `PPT` | `PPT / 1800` (mm/30min → kg/m²/s) |
| `Wind` | m/s | `ubar` | as-is, floor `u_min = 0.1` |
| `LWdown` | W/m² | `Rl_dn` | as-is; gap-fill (§7.3) |
| `SWdown_par_beam` | W/m² | `Par_tot`, `Par_diff` | `(Par_tot − Par_diff) / 4.6` |
| `SWdown_par_diffuse` | W/m² | `Par_diff` | `Par_diff / 4.6` |
| `SWdown_nir_beam` | W/m² | `Rs`, `Par_tot`, `Par_diff` | `(1−fdiff)·max(0, Rs − Par_tot/4.6)` |
| `SWdown_nir_diffuse` | W/m² | ″ | `fdiff·max(0, Rs − Par_tot/4.6)`, `fdiff = Par_diff/Par_tot` |
| `CO2air` | µmol/mol | (none) | constant 400 (attribute-flagged synthetic) |

Global attributes: `latitude`, `longitude`, `elevation`, `timestep_seconds`, `avg_convention`,
`source = "BCI_v5.1.csv"`, provenance (fill methods, CO₂ constant).

### 7.2 Pseudo-code

```python
import pandas as pd, numpy as np, netCDF4

LAT, LON, ELEV = 9.15, -79.85, 120.0
DT_NATIVE = 1800.0        # s (half-hourly)
PAR_W_PER_UMOL = 1.0/4.6  # umol PAR -> W
CO2_CONST = 400.0

df = pd.read_csv("BCI_flux/BCI_v5.1.csv", parse_dates=["date"])

# 1. Unit conversions
Tair  = df.tair + 273.15
PSurf = df.p_kpa * 1000.0
def rh_to_q(rh, T_K, P_Pa):                       # Bolton-1980 esat — MUST match meds_thermo (§4.3)
    tc   = T_K - 273.15
    esat = 611.2*np.exp(17.67*tc/(tc + 243.5))    # identical to meds_thermo%sat_vapor_pressure
    e    = np.clip(rh,0,1)*esat
    return 0.622*e/(P_Pa - 0.378*e)
Qair  = rh_to_q(df.RH/100.0, Tair, PSurf)
Rainf = df.PPT.clip(lower=0)/DT_NATIVE
Wind  = df.ubar.clip(lower=0.1)

# 2. Shortwave: 4-stream split from measured PAR + total SW
Rs       = df.Rs.clip(lower=0)                      # night ~ -0.4 -> 0
par_tot_W= df.Par_tot.clip(lower=0)*PAR_W_PER_UMOL
par_dif_W= df.Par_diff.clip(lower=0)*PAR_W_PER_UMOL
fdiff    = np.where(df.Par_tot>1.0, df.Par_diff/df.Par_tot, 1.0).clip(0,1)
nir_W    = np.maximum(0.0, Rs - par_tot_W)
out = dict(SWdown_par_beam    = (par_tot_W - par_dif_W).clip(lower=0),
           SWdown_par_diffuse = par_dif_W,
           SWdown_nir_beam    = (1-fdiff)*nir_W,
           SWdown_nir_diffuse = fdiff*nir_W)

# 3. Gap-fill (LWdown etc.) — persistence + by-hour-of-day climatology
def gapfill(x, hour):
    x = x.copy(); clim = pd.Series(x).groupby(hour).transform(lambda s: s.mean())
    x = x.fillna(clim); return x.ffill().bfill()
LWdown = gapfill(df.Rl_dn, df.date.dt.hour)         # daytime SW streams are 0 at night -> no fill needed

# 4. (optional) 30 min -> 15 min: linear for state vars, step-hold for Rainf,
#    cosz-weighted for SW (mirror the Fortran kernel so file == on-the-fly agree). Off by default:
#    the Fortran reader disaggregates to dt_fast, so writing native 30 min is sufficient.

# 5. Write NetCDF (unlimited time; ALMA names; base-time units; global attrs lat/lon/elev/convention)
write_netcdf("BCI_v5.1_forcing.nc", time=..., Tair=Tair, Qair=Qair, PSurf=PSurf, Rainf=Rainf,
             Wind=Wind, LWdown=LWdown, CO2air=np.full(len(df),CO2_CONST), **out,
             attrs=dict(latitude=LAT, longitude=LON, elevation=ELEV,
                        timestep_seconds=DT_NATIVE, avg_convention="end"))
```

The optional 30→15 min disaggregation (step 4) is **off by default**: the Fortran reader already
disaggregates to `dt_fast` on the fly, so the NetCDF stays at native 30 min (smaller, and the
authoritative diurnal reconstruction lives in one place — the Fortran `cosz_reconstruct_factor` kernel).
It exists for validating the Fortran disaggregation against a Python reference (§9).

---

## 8. Phasing

**P0 — single-site reader wired to the fast loop (MVP).**
- `met_forcing_t`, `met_driver_t`, `forcing_config_t`; `[site]` + `[forcing]` config; `req_*` loaders.
- **`meds_time` second-level helpers** (`seconds_into_day`, `time_advance_seconds`, `seconds_between`)
  — a hard prerequisite; the interpolation and window slide do not compile without them (§2.1).
- `meds_met_driver`: **CSV backend** (BCI v5.1) + **CONST** backend; `met_open`/`advance`/`instant`/
  `close`; **linear** interpolation for state vars, **step-constant** precip, **reciprocal-mean-cosz**
  SW reconstruction (the existing `solar_cosz` + the new `cosz_reconstruct_factor` kernel — NOT
  ED2's `⟨sec z⟩`), `meds_thermo%air_density` (reused), `rh_to_specific_humidity` (Bolton esat).
- Wiring: per-sub-step met refresh in `run_fast_biophysics` (§6.2 shim); `meds_main` builds the
  driver, seeds reservoirs, passes `fast_ctx` + date (closes the gap). Per-cohort SW stays the
  LAI-share split of the **time-varying** `rad_sw_top` (real RT deferred).
- `scripts/prep_bci_forcing.py` (NetCDF is optional at P0 since the CSV backend works directly).
- Test: `test_met_driver` reproduces a diurnal GPP/temperature cycle offline (§9).

**P1 — physical fidelity.**
- **NetCDF/ALMA backend** (reads the prep output & PLUMBER2 files via `meds_netcdf_c`).
- **Canopy RT join** (§6.3): `rad_forcing_t` from `met_forcing_t`, per-cohort absorbed SW/PAR from
  `meds_canopy_radiation` replaces the LAI split; `cosz`/band SW/`lwdown` fully wired.
- **Gap-fill climatology**, **LWdown synthesis** (Brutsaert), **Weiss–Norman** SW partition,
  **multi-year cycling** (`recycle`, Feb-29 reconciliation), the **`met_forcing_t`-argument cleanup**
  (retire the `apply_met_to_ctx` shim), the **daily accumulator** feeding phenology (`temp_day`,
  `daylength`, `doy`).

**P2 — gridded / multi-polygon.**
- Multi-polygon forcing (a polygon dimension on the file; nearest-land grid match, ED2 `match_poly_
  grid`), elevation **lapse-rate** correction (`lapse.f90` analogue) between grid and site,
  climate-change intercept/slope perturbations, the full `sw_partition`/`avg_convention` selector for
  reanalysis products (CRUNCEP/GSWP3/WFDE5/ERA5), transient CO₂ (a CO₂-only stream that does not
  cycle, ED2's special case).

---

## 9. Test plan

Offline, driving the physics kernels with BCI to reproduce a **diurnal cycle** — the model-behaviour
analogue of ED2's EDTS, and the first end-to-end use of the fast loop in production.

1. **Reader unit tests (`test_met_driver`).** (a) `met_open` on a 3-row synthetic CSV loads records
   #1–2; (b) `met_advance` slides the window correctly across a record boundary, **recycles** at EOF,
   and takes the **start-before-`base_time`** path (`start_clamp=error` stops; `=hold` clamps to record
   #1 with `w_next=0`); (c) `met_instant` at the exact record time returns the record value; at the
   midpoint returns the linear mean (state vars) / the previous value (precip, step-constant);
   (d) **conservation** — the reciprocal-mean-cosz SW reconstruction integrated over a window returns
   the interval mean to round-off: `(1/T_win) ∫ F_avg·cosz(t)/⟨cosz⟩_win dt = F_avg` (this closes
   **because** the factor is `1/⟨cosz⟩_win`, not `⟨sec z⟩`; a companion assert shows the `⟨sec z⟩` form
   does **not** close and is biased high on a sunrise window), and step-constant precip conserves
   accumulation.
2. **Solar-geometry sanity.** `solar_cosz` at BCI (9.15 °N) peaks near local noon; `cosz_reconstruct_
   factor` returns `1/⟨cosz⟩_win` over daytime sub-samples (night clamped to 0); night SW is exactly 0
   (`cosz ≤ cosz_min`). Cross-check the day-length against the astronomical value, and **assert the
   modeled solar-noon offset stays under one window** — this quantifies (and bounds) the ~20–35 min
   local-standard-time phase error of the P0 clock-as-solar assumption (§5.1); with
   `apply_solar_lon=true` the offset should collapse to near-zero.
3. **Kernel offline drive (the headline).** Feed one day of BCI forcing (48 records → 96 fast
   sub-steps) into `column_fast_step` for a single on-allometry cohort and assert a **physical diurnal
   cycle**: GPP tracks PAR (zero at night, peak near noon), `leaf_temp` follows `Tair` + a radiative
   offset, transpiration tracks VPD, the CAS/soil budgets close each step (`budg%whole_*%worst` below
   tolerance). Compare GPP against the BCI `gpp` column (order-of-magnitude / diurnal-shape check —
   not a calibrated match, since RT/hydraulics are still coarse at P0, and **P0 leaf PAR is biased
   high**: `leaf_env%par = abs_sw/lai·par_per_w` converts the LAI-share of *total* SWdown, treating
   absorbed NIR as PAR-convertible, until the P1 RT join supplies true per-cohort absorbed PAR).
4. **CONST near-no-op refactor (re-specified — NOT bit-for-bit).** With `format = "const"`,
   `run_fast_biophysics` reproduces a **documented CONST baseline**: canopy-top SW is exactly 400 W/m²
   (the four reference streams sum to 400, held diurnally flat), and air T / humidity / CO₂ / wind /
   pressure match today's constants, so GPP and the CAS/soil budgets agree with the current
   `test_fast_loop` **within round-off**. It is **not** byte-identical (`cosz`/`rho_air` are recomputed
   each substep; ground SW flows through `rad_sw_ground_const=60`), and **`lwdown`/`abs_lw` is not
   asserted** (dead until the P1 RT join, §6.2).
5. **Python↔Fortran agreement.** The prep script's optional 30→15 min disaggregation matches the
   Fortran `cosz_reconstruct_factor` reconstruction within round-off on a shared day, **and the prep's
   `rh_to_q` (Bolton esat) matches the reader's `rh_to_specific_humidity` `qair` to round-off** — the
   single-saturation-formula check that guarantees file == on-the-fly (§4.3).
6. **Precip-phase / frozen branch (synthetic).** BCI is all-liquid, so the `precip_phase` frozen path
   never fires on real data; a **synthetic sub-freezing record** (`tair_k < t_3ple`) exercises it and
   asserts total precip splits conservatively into `rainf`+`snowf` (mass conserved, `snowf > 0` below
   the phase band) — otherwise the frozen branch ships untested.
7. **Fast→slow GPP handoff units.** With time-varying SW, assert the daily-integrated `gpp_accum`
   carries the correct units (`[kgC/plant/day]` at the seam `carbon_growth` reads) and is **reset once
   per slow day**, so a multi-day CONST run gives a constant daily GPP and a multi-day BCI run gives a
   day-varying one — pinning the reset/units the always-pass-`fast_ctx` wiring (§6.5) makes safe.
8. **Portability.** Build all new modules under **nvfortran multicore** (not just ifx): the
   `partition_shortwave` kernels return arrays — bind results to named arrays before any call (trap
   #7). Confirm the fast loop with real forcing runs on the multicore back end.

---

## 10. Open questions

1. **BCI unit ambiguities to confirm against the site metadata.** (a) `PPT` — assumed **mm per 30-min
   interval** (→ `/1800` for kg/m²/s); if it is already a rate (mm/s) or mm/hr the divisor changes.
   (b) `vpd` — assumed **kPa**; ingested only as a cross-check (RH+T is authoritative). (c) `Rs` vs
   `Rs_dn` — assumed `Rs` is incoming SWdown (night ≈ −0.4 → clamp 0) and `Rs_dn` a mostly-NaN
   alternate; confirm which is the calibrated downwelling stream. (d) `avg_convention` — assumed the
   30-min timestamp labels the **interval end** (FLUXNET convention); a mislabel phase-shifts the
   diurnal cycle by up to 30 min.
2. **PAR→W factor.** Used 4.6 µmol/J (standard for PAR); if the tower reports a different waveband or
   the site uses 4.57/4.6/2.02 conventions, the NIR complement `Rs − Par_tot/4.6` shifts. A per-site
   attribute in the prep output makes this explicit.
3. **Reference height vs `hgt_max`.** ED2 aborts if `zref ≤ hgt_max` of any PFT. BCI canopy ~40 m and
   tall-PFT `hgt_max` may approach the tower height — confirm `reference_height` clears every PFT cap,
   or the run hard-stops (by design).
4. **Where the daily phenology aggregation lives.** The `met_daily_t` accumulator (§6.4) could live on
   the driver (simplest) or on the site (visible to a future MPI/polygon layer). Deferred to P1 with
   phenology wiring; flag for the master-loop design.
5. **`met_forcing_t` on `fast_context_t` vs a first-class argument.** P0 uses the `apply_met_to_ctx`
   shim (minimal diff); P1 retires it. Decide at P1 whether `column_forcing_t`/`aero_env_t` should
   take `met_forcing_t` directly (cleaner, removes the co2/temp duplication the brief flags across ~5
   types) versus keeping the `fast_context_t` carrier for backward compatibility with `test_fast_loop`.
6. **Soil-moisture initialization from `SWC`.** BCI reports `SWC`; P0 ignores it (uses
   `ctx%theta_init`). Whether to seed `soil_w%theta` from the tower `SWC` (single-depth) is an
   init-side question, not a forcing one — but the reader is the natural place to surface it.
