# Meteorological forcing

Forcing is the **prescribed** half of the model's boundary: the atmospheric state above the canopy,
*read from a file* rather than evolved. Everything else in MEDS is state, advanced by a kernel that
owns it; nothing in the model writes forcing. One polygon / one site runs today, bound to one column
of a possibly multi-site file. The reader produces one **instantaneous per-site record**
(`met_forcing_t`) per `dt_fast` sub-step: air temperature, specific humidity, surface pressure, wind,
downwelling longwave, the four (beam/diffuse)×(PAR/NIR) shortwave streams, CO₂, plus the *derived*
$`\cos z`$ and $`\rho_{air}`$ — never stored in a file, always recomputed.

`libmeds_forcing` links the config leaf and the netCDF C bindings (`meds_netcdf_c`) and **nothing
else** — no state, no demography, no kernels. That placement is the point: a prescribed driver must sit
*low* in the dependency graph, because everything consumes it and it consumes nothing. It is also why
the disaggregation math is a separate `pure`/`elemental` kernel module from the reader that does I/O.

## 1. The forcing file

One format: an **unstructured multi-grid NetCDF** whose forcing variables are all shaped `(time, grid)`.
`time` is the record axis; `grid` is a plain *list* of locations (`grid = 1` for a single site), because
a polygon list is a list, not a raster — a regular reanalysis tile is flattened `y×x → grid`. Variables
marked optional below may be omitted; every other one is required.

```
netcdf ithaca_forcing {
dimensions:
    time = UNLIMITED ;  grid = 1 ;                     // grid = N for a multi-site file
variables:
    double time(time) ; time:calendar = "proleptic_gregorian" ;
        time:units = "seconds since 2024-01-01 01:00:00" ;      // the base-time anchor
    double latitude(grid) ["degrees_north"], longitude(grid) ["degrees_east"] ;
    double elevation(grid) ["m"] ;                              // optional, not read
    float Tair(time,grid) ["K"], Qair(time,grid) ["kg kg-1"],       // cell_methods = "time: point"
          PSurf(time,grid) ["Pa"], Wind(time,grid) ["m s-1"] ;
    float Rainf(time,grid) ["kg m-2 s-1"], LWdown(time,grid) ["W m-2"],  // cell_methods = "time: mean"
          SWdown(time,grid) ["W m-2"] ;                         // SWdown is the TOTAL; see sec. 6
    float CO2air(time,grid) ["umol mol-1"] ;                    // optional -> forcing.co2_const
// global attributes -- provenance only, never read by the model:
    :Conventions = "MEDS-forcing-1.0" ; :source = "ERA5-Land hourly" ; :time_zone = "UTC" ;
    :timestep_seconds = 3600 ; :avg_convention = "end" ; :sw_input_kind = "total" ;
}
```

With `sw_partition = "passthrough"`, `SWdown` is replaced by four pre-split streams `SWdown_par_beam`,
`SWdown_par_diffuse`, `SWdown_nir_beam`, `SWdown_nir_diffuse` [W m⁻²], all required.

**The time axis is the only metadata the reader parses:** it takes everything after `since` in
`time:units` as the base instant and the `time` values as seconds from it. Record interval, averaging
convention and shortwave scheme come from the TOML `[forcing]` block, *not* from the global attributes —
those are for a human, and a file whose attributes disagree with the config is not detected.

**Binding a site to a column.** `grid_match = "explicit"` uses `forcing.grid_index` verbatim (range checked
at open). `grid_match = "nearest"` instead reads the `latitude(grid)`/`longitude(grid)` vectors and picks
the cell minimising the great-circle distance to the `[site]` coordinates (`great_circle_distance`, ED2's
`dist_gc`), with strict `<` in the scan so ties keep the lowest index; only the ordering matters, so the
Earth radius is immaterial. Every value is read as one `(time, grid)` hyperslab of count `[1,1]` — the
reader never loads a variable's whole time series, only the cached `time` coordinate.

**Producing one from ERA5-Land.** `scripts/download_era5land.py` pulls the eight hourly variables (`t2m`,
`d2m`, `sp`, `u10`, `v10`, `tp`, `ssrd`, `strd`) for a small box around the site from the CDS;
`scripts/prep_era5land_forcing.py` writes the file above — `Tair = t2m`, `PSurf = sp`, `Qair` from the
dewpoint by (10) in the *same* Bolton form the model uses, $`\mathrm{Wind}=\sqrt{u_{10}^2+v_{10}^2}`$, and
the three **accumulated** fluxes de-accumulated then unit-converted:
$`\mathrm{Rainf}=\Delta tp\cdot 1000/3600`$ [kg m⁻² s⁻¹], $`\mathrm{SWdown}=\Delta ssrd/3600`$,
$`\mathrm{LWdown}=\Delta strd/3600`$ [W m⁻²].

*The 00Z trap.* ERA5-Land accumulations run from 00 UTC and reset daily, so the **00:00 stamp carries the
whole previous day's total** (step 24) — not zero, and not one hour. Ordered by valid time, the per-hour
amount is `raw(H) - raw(H-1)` everywhere *including* the 00:00 stamp (which then yields the prior day's
23Z–00Z hour), except at **01:00 UTC**, where it is `raw(01:00)` as-is (the first step of the period, with an
implicit 0 at 00:00). An off-by-one here doubles or zeroes the first hour of every day. Recovering day D's
last hour needs the 00:00 stamp of day D+1, so the download pads one extra trailing day; the single leading
sample that cannot be differenced is **dropped** and the time axis rebased, never filled, and GRIB-packing
negatives are clipped to zero.

## 2. The reader: a two-record window

`met_open` reads the dimensions, base time and time coordinate, validates the recycle window (§9) and
loads records #1–#2; `met_advance(drv, now)` slides the bracket so `rec_prev%when ≤ now < rec_next%when`,
marching the cursor incrementally; `met_close` releases the handle. `met_instant(drv, now)` is a **pure
function of (reader state, time)** — it interpolates and disaggregates the loaded bracket to the instant,
so it is safe to sample ahead of the threaded patch loop. A run starting before the first record either
hard-errors (`start_clamp = "error"`, the default) or holds record #1 (`"hold"`); a non-recycling run past
the last record clamps to the final interval. The no-file **`const`** backend returns the `met_forcing_t`
defaults — a reference climate whose four shortwave streams sum to 400 W m⁻², $`\cos z`$ and
$`\rho_{air}`$ still derived — so the fast loop can run with no forcing file at all.

## 3. Temporal interpolation, and why wind is different

Within the bracket the weight is $`w_{next}=(t-t_{prev})/(t_{next}-t_{prev})`$, clipped to $[0,1]$. Each
variable carries its own policy: **linear** for `Tair`, `Qair`, `PSurf`, `LWdown`, `CO2air` (smooth
atmospheric states); **step-constant** for `Rainf`, holding the previous value, because interpolating
precipitation smears intense events into physically wrong drizzle and breaks infiltration and runoff;
**cosz reconstruction** for the four shortwave streams (§5), because those records are interval *means*,
not point values; and for wind, the **energy form**:

```math
u(t) = \sqrt{\max\!\big[(1-w_{next})\,u_{prev}^2 + w_{next}\,u_{next}^2,\; u_{min}^2\big]},
\qquad u_{min} = 0.1\ \mathrm{m\,s^{-1}} \qquad(1)
```

Wind interpolates as a **squared** quantity, then square-rooted (ED2's `vels` convention), so it is the
kinetic-energy content, not the speed, that varies linearly across a record boundary; the floor keeps the
Monin-Obukhov solve in the aerodynamics kernel out of its degenerate zero-wind limit.

## 4. Solar geometry: apparent solar time

`meds_time%solar_cosz` expects **local apparent solar seconds** — it reads its time argument as a fraction
of day with noon at 0.5. A reanalysis file is stamped in UTC, so feeding clock seconds straight in would put
solar noon at 00:00 UTC everywhere. One transform fixes it:

```math
t_{solar} = t_{clock} + (\lambda - 15\,\Delta_{UTC})\cdot 240\ \mathrm{s\,deg^{-1}}
          + \mathrm{EoT}(\mathrm{doy}) \qquad(2)
```

```math
\mathrm{EoT} = 229.18\cdot 60\,\big(0.000075 + 0.001868\cos b - 0.032077\sin b
              - 0.014615\cos 2b - 0.040849\sin 2b\big)\ \mathrm{s},
\quad b = \frac{2\pi(\mathrm{doy}-1)}{365} \qquad(3)
```

Here $\lambda$ is the longitude (+E) and $`\Delta_{UTC}`$ the file clock's UTC offset in hours — zero for
a UTC file, which collapses (2) to the longitude plus the Spencer (1971) equation of time.

`apply_solar_longitude = false` passes the clock seconds through untouched, for a file already in local
apparent solar time. The zenith cosine is then the standard expression in declination $\delta$ and hour
angle $h$, floored at zero (night):

```math
\cos z = \sin\phi\sin\delta + \cos\phi\cos\delta\cos h, \quad
\delta = 23.45^\circ\sin\!\frac{2\pi(284+\mathrm{doy})}{365}, \quad
h = 2\pi\!\left(\frac{t_{solar}}{86400}-\tfrac12\right) \qquad(4)
```

`met_instant` recomputes $`\cos z`$ every sub-step from the **model** clock; it is never read from a file
and never interpolated.

## 5. Shortwave disaggregation: conserving the interval mean

This is the subtle piece. A shortwave record on an `avg_convention = "end"` file is the **mean over the
interval ending at its stamp**, not a value at the stamp. Redistributing that mean inside the interval
must satisfy an exact requirement: *the reconstructed flux, averaged back over the same window, must
return the recorded mean.* If the instantaneous flux is proportional to $`\cos z`$ by day and zero at
night, that fixes the normalizer uniquely:

```math
\langle F\rangle_{win} = f_\perp\,\langle\cos z\rangle_{win}
\quad\Longrightarrow\quad
F(t) = F_{avg}\,\frac{\cos z(t)}{\langle\cos z\rangle_{win}} \qquad(5)
```

`cosz_reconstruct_factor` returns **$`1/\langle\cos z\rangle_{win}`$ — the reciprocal of the mean cosine —
and emphatically not $`\langle\sec z\rangle`$.** ED2's `mean_daysecz` computes the mean secant; by Jensen's
inequality $`\langle 1/\cos z\rangle \ge 1/\langle\cos z\rangle`$, with the gap blowing up near sunrise and
sunset, so the secant form reconstructs a flux that does *not* integrate back to $`F_{avg}`$ and is biased
high at low sun. MEDS does not port it, and the test asserts both halves. The mean is a 10-sub-sample
midpoint rule over the **full** window, night sub-samples contributing exactly zero:

```math
\langle\cos z\rangle_{win} = \frac{1}{N}\sum_{i=1}^{N}\max\!\big(\cos z(t_i),\,0\big),
\qquad t_i = t_0 + \left(i-\tfrac12\right)\frac{\Delta t_{win}}{N} \qquad(6)
```

Averaging over the whole window rather than its daylit part is what makes a sunrise or sunset window close.
The guard is on the aggregate: $`\langle\cos z\rangle_{win}\le 10^{-3}`$ (a fully-night window, where
$`F_{avg}`$ is zero anyway) returns a factor of 0, and each stream is separately zeroed when
$`\cos z(t)\le 10^{-3}`$; there is no per-step secant to blow up.

**The reconstruction is anchored on the MODEL window, not the file's** — the window start is found by
stepping back from the model instant by the elapsed fraction of the bracket, so $`\langle\cos z\rangle_{win}`$
and $`\cos z(t)`$ see the same sun. Under calendar recycling (§9) the two calendars differ, and anchoring on
the model's is what keeps identity (5) true in the year actually simulated. The four streams are
disaggregated independently, each from the record carrying the interval mean — `rec_next` for
`avg_convention = "end"`, `rec_prev` for `"begin"`.

## 6. Splitting shortwave into four streams

Reanalysis ships total downwelling shortwave; the canopy two-stream needs beam and diffuse, in the visible
and the near-infrared. The split happens **at ingest**, once per record, on that record's interval-midpoint
$`\cos z`$, so `met_record_t` always holds four streams. Both schemes conserve energy exactly.

**`clearidx` (default)** — the Erbs et al. (1982) clearness-index correlation. With
$`k_t = \mathrm{SW}/(S_0\cos z)`$ clipped to $[0,1]$ and $`S_0 = 1361`$ W m⁻²,

```math
f_{diff} =
\begin{cases}
1 - 0.09\,k_t & k_t \le 0.22\\
0.9511 - 0.1604\,k_t + 4.388\,k_t^2 - 16.638\,k_t^3 + 12.336\,k_t^4 & 0.22 < k_t \le 0.80\\
0.165 & k_t > 0.80
\end{cases} \qquad(7)
```

Then $`F_{diff}=f_{diff}\,\mathrm{SW}`$ and $`F_{beam}=\mathrm{SW}-F_{diff}`$, each divided into PAR and
NIR by **fixed energy fractions** — 0.43 of the beam and 0.52 of the diffuse are PAR, the rest NIR; diffuse
light is PAR-enriched because Rayleigh scattering is strongest at short wavelengths.

**`weiss_norman`** — the band-specific Weiss & Norman (1985) scheme, a port of ED2's
`short_bdown_weissnorman`. It builds *potential* (clear-sky) band fluxes from air-mass optical depth, with
$`m=\sec z`$ and $`p_{rat}=P_{surf}/P_{std}`$ (so it needs surface pressure), then scales them to the
observed total. Writing $`S_{vis}=0.43\,S_0`$, $`S_{nir}=0.57\,S_0`$, and the NIR water-vapour absorption
$`w_{10}=S_0\,10^{\,-1.1950+0.4459\log_{10}m-0.0345(\log_{10}m)^2}`$:

```math
R^{pot}_{vis,b} = S_{vis}e^{-0.185\,p_{rat}m}\cos z, \qquad
R^{pot}_{vis,d} = 0.40\left(S_{vis}-R^{pot}_{vis,b}\right)\cos z
```

```math
R^{pot}_{nir,b} = \left(S_{nir}e^{-0.060\,p_{rat}m}-w_{10}\right)\cos z, \qquad
R^{pot}_{nir,d} = 0.60\left(S_{nir}-R^{pot}_{nir,b}-w_{10}\right)\cos z \qquad(8)
```

Each band's potential total is $`R^{pot}_{full}=R^{pot}_{b}+R^{pot}_{d}`$. With
$`r=\mathrm{SW}/(R^{pot}_{vis,full}+R^{pot}_{nir,full})`$ the observed-to-potential ratio (cloudier ⇒
smaller ⇒ more diffuse), the band total is $`r\,R^{pot}_{full}`$ and the actual beam fraction is

```math
f_{beam} = \frac{R^{pot}_{beam}}{R^{pot}_{full}}
\left[1-\left(\frac{a-\min(a,\max(0,r))}{c}\right)^{2/3}\right]_{0}^{1},
\quad (a,c)=(0.90,\,0.70)\ \mathrm{vis},\ (0.88,\,0.68)\ \mathrm{NIR} \qquad(9)
```

Below $`\cos z\le\cos 89^\circ`$ the secant is unstable and everything is routed to diffuse, 0.52 visible /
0.48 NIR. One ED2-faithful edge survives the port: in a narrow twilight band just above that cutoff
$`w_{10}`$ can exceed the attenuated NIR beam, so the scaled NIR diffuse goes slightly negative while the
four streams still sum to `SWdown` exactly.

## 7. Humidity and precipitation phase

Both humidity conversions go through the **same** Bolton (1980) saturation vapour pressure the rest of MEDS
uses (`meds_therm_lib%sat_vapor_pressure`, $`e_{sat}(T)=611.2\exp[17.67\,T_c/(T_c+243.5)]`$ Pa) — and so
does the ERA5-Land prep script, so a `Qair` built offline reconciles with any reader-side humidity math to
round-off. The dewpoint form is the identity that the actual vapour pressure *is* the saturation vapour
pressure evaluated at the dewpoint; RH is clipped to $[0,1]$ first.

```math
q = \frac{0.622\,e}{P-0.378\,e}, \qquad
e = e_{sat}(T_d)\ \text{(dewpoint)}, \qquad
e = \mathrm{RH}\cdot e_{sat}(T)\ \text{(relative humidity)} \qquad(10)
```

The file carries one total precipitation rate; the model needs rain and snow separately. The split is a
linear ramp across a 1 K half-width band centred on the triple point (ED2's Jin 1999 form, simplified),
mass-conserving by construction, on the *interpolated* sub-step air temperature rather than the record
value — so the phase follows the diurnal cycle within an interval straddling freezing.

```math
f_{liq} = \left[\frac{T-(T_3-1\,\mathrm{K})}{2\,\mathrm{K}}\right]_{0}^{1}, \qquad
\mathrm{rainf} = f_{liq}P, \qquad \mathrm{snowfall} = P-\mathrm{rainf} \qquad(11)
```

## 8. Wind height and elevation lapse

Two optional ingest-time corrections, both **off by default**, both applied per record before interpolation.
Reanalysis wind is diagnostic at 10 m, which at a temperate forest can sit below canopy top; a neutral-log
profile lifts it to the model reference height. The factor is independent of $u$, so it commutes with the
energy-form interpolation of §3; a degenerate roughness ($`z_0\le 0`$, or either height at or below
$`z_0`$) leaves the wind untouched.

```math
u_{ref} = u_{meas}\,\frac{\ln(z_{ref}/z_0)}{\ln(z_{meas}/z_0)} \qquad(12)
```

For elevation, with $`\Delta z = z_{site}-z_{grid}`$ and a positive environmental lapse rate $\Gamma$
(cooling upward), temperature and pressure move **together** so they stay ideal-gas consistent — the
pressure form is the hypsometric integral of $`dP/dz=-Pg/(R_dT)`$ under the *same* linear $T(z)$, not an
independent barometric guess:

```math
T_{site} = T_{grid}-\Gamma\,\Delta z, \qquad
P_{site} = P_{grid}\left(\frac{T_{site}}{T_{grid}}\right)^{g/(R_d\Gamma)} \qquad(13)
```

with the isothermal limit $`P_{site}=P_{grid}\exp[-g\,\Delta z/(R_dT_{grid})]`$ when
$`|\Gamma|\le 10^{-6}`$. Specific humidity is carried across unchanged; $`\rho_{air}`$ is re-derived.

**The reference height must clear the canopy.** `site.reference_height` is validated against every PFT's
`hgt_max` at config load and the run stops unless it exceeds all of them (ED2 aborts on the same condition):
a forcing height inside the canopy makes surface-layer similarity meaningless, silently.

## 9. Calendar recycling

A short forcing record can drive a long run by cycling. The rule: **the cycle window is DECLARED, never
inferred.** `forcing.recycle_start` and `forcing.recycle_end` are required whenever `forcing.recycle = true`,
and the window is half-open — `recycle_end` is the same instant one cycle later, so each record is covered
exactly once. It is validated three ways: **at config load**, both dates valid, end after start, and the span
an *exact whole number of calendar years*; **at `met_open`**, `recycle_start` must land exactly on a record
stamp (within 0.5 s — on an `"end"` file the first record of a calendar year is 01:00:00 for hourly data, not
00:00:00); and again at `met_open`, the file must cover the window, holding a record in the last interval
below `recycle_end`.

The whole-year requirement is the substantive one. A window that is not a whole number of calendar years
drifts **both** hour-of-day and day-of-year on every wrap — and because the shortwave reconstruction of §5
is mean-conserving, the daily mean stays correct while the sub-daily phase and the season slide. Nothing
downstream complains: the demography sees a plausible annual cycle over a scrambled one. Declaring the
window, and rejecting one that cannot wrap cleanly, is the only place that error is visible. Mapping a model
instant into the window is then a **year substitution**, anchor-relative so the cycle may begin anywhere in
the calendar:

```math
\mathrm{off} = \begin{cases}
0 & (\text{month},\text{day},\text{hh:mm:ss})_{model}\ \ge\ (\cdot)_{anchor}\\
1 & \text{otherwise}\end{cases}, \qquad
y_f = Y_1 + \mathrm{off} + \mathrm{mod}\!\left(Y_m-Y_1-\mathrm{off},\,N\right) \qquad(14)
```

with $`Y_1`$ the anchor year, $`Y_m`$ the model year and $N$ the cycle length in years. Month, day and **time
of day are preserved exactly** — that exactness is what lets the sub-daily phase and the day-of-year survive
an arbitrary number of wraps; for a Jan-1 00:00 anchor `off` is identically 0 and (14) reduces term for term
to plain year substitution. The one exception is the **leap day**: a Feb-29 model instant reads Feb-28 when
the target file year is not a leap year (ED2's `read_ol_file` repeats Feb 28), one-directionally — a non-leap
model year never asks for Feb-29 — so no record is skipped or double-counted. At the cycle edge the reader
loads a **seam bracket**: `rec_prev` is the last record inside the window, `rec_next` is the window's *first*
record, so the closing interval interpolates across the wrap instead of clamping. The window may be a
sub-range of a longer file; records past `recycle_end` belong to the next file year, not to this cycle.

## 10. MEDS never gap-fills

Every required field is checked for NaN as it is read, and a missing value **halts the run**, naming the
field, the time record and the grid index. There is no gap policy, no persistence hold, no climatology fill —
and the prep script applies the same rule, erroring rather than inventing a value. This is deliberate:
gap-filling is a modelling decision with real consequences for the fluxes it produces, and it belongs
upstream in the user's own preprocessing, where it is visible and documented. Buried in a reader, it turns a
silently degraded run into one that looks clean. What is *not* a gap: the deterministic derivations — the
total→four-stream partition, humidity from dewpoint, the de-accumulation boundary sample the prep script
drops rather than writes — which compute a variable the source does not store from variables it does.

## What is not here

- **No longwave synthesis.** `lwdown_source = "synthesize"` is a declared selector with no implementation,
  and the config validator **rejects** it rather than silently running the file path under a name that
  promises Brutsaert/Idso clear-sky synthesis — `LWdown` must be in the file.
- **No multi-polygon runtime.** The `(time, grid)` format, `grid_index` and nearest-cell matching are in
  place — file and reader are ready for N locations — but the model runs one site.
- **No transient CO₂ stream.** `CO2air` may be in the file; a run that omits it takes one constant
  (`forcing.co2_const`, default 420 µmol/mol) as the single free-atmosphere authority.
- **`avg_convention`** distinguishes `"end"` and `"begin"`. `"instant"` and `"center"` parse, but the
  disaggregation treats them as end-of-interval; only the partition's midpoint offset (§6) differs.

See [`docs/ROADMAP.md`](../ROADMAP.md) §8 for what is planned, and when.

## Where the code is

| Concept | Routine / file |
|---|---|
| interpolation, wind energy form | `meds_forcing_kernels`: `interpolate_forcing`, `interpolate_wind_energy` |
| apparent solar time, zenith | `meds_forcing_kernels`: `apparent_solar_seconds`, `equation_of_time`, `met_solar_cosz` over `meds_time%solar_cosz` |
| SW disaggregation | `meds_forcing_kernels`: `cosz_reconstruct_factor`, `disaggregate_sw` |
| SW partition | `meds_forcing_kernels`: `partition_shortwave`, `erbs_diffuse_fraction`, `weiss_norman_partition` |
| humidity, precip phase | `meds_forcing_kernels`: `dewpoint_to_specific_humidity`, `rh_to_specific_humidity`, `precip_phase` |
| grid match, wind, lapse | `meds_forcing_kernels`: `great_circle_distance`, `nearest_grid_index`, `wind_log_profile`, `lapse_air_temperature`, `lapse_pressure` |
| the reader | `meds_met_driver`: `met_open`, `met_advance`, `met_instant`, `met_close`; `read_record`, `assert_finite` |
| recycling | `meds_met_driver`: `validate_recycle_window`, `file_lookup_sec`, `recycle_model_to_file`, `load_wrap_bracket` |
| types | `meds_forcing_types`: `met_forcing_t`, `met_record_t`, `met_driver_t` |
| config + selectors | `meds_forcing_config`: `forcing_config_t`, `INTERP_*`, `SWPART_*`, `LW_*`, `CLAMP_*`, `METAVG_*`, `GRIDMATCH_*`; validated in `meds_config`, read by `meds_config_io` |
| TOML block | `[forcing]` + `[site]` (documented in `meds_config_main.toml`) |
| fast-loop join | `meds_fast_dynamics`: per-sub-step `met_advance`/`met_instant` sampling, `apply_met_to_ctx` |
| file production | `scripts/download_era5land.py`, `scripts/prep_era5land_forcing.py` |
| test | `test/test_met_driver.f90` — interpolation, humidity, phase, both SW schemes, the mean-conserving identity *and* the secant bias, CONST backend, NetCDF round-trip, clamp and recycle-window rejections, recycle phase over 29 years |

## References
- Erbs, Klein & Duffie (1982), *Solar Energy* 28:293 — diffuse fraction vs the clearness index.
- Weiss & Norman (1985), *Agric. For. Meteorol.* 34:205 — band-specific direct/diffuse partitioning;
  ED2 `../ED2/ED/src/utils/radiate_utils.f90` (`short_bdown_weissnorman`).
- Bolton (1980), *Mon. Weather Rev.* 108:1046 — saturation vapour pressure.
- Spencer (1971), *Search* 2:172 — Fourier series for the equation of time.
- Jin et al. (1999) — precipitation phase partitioning; ED2 `ed_met_driver.f90`.
- Longo et al. (2019), *GMD* 12:4309 — ED-2.2 technical description (the meteorological driver).
- Muñoz-Sabater et al. (2021), *ESSD* 13:4349 — ERA5-Land.
- Design doc: `docs/dev_plans/MEDS_FORCING_DESIGN.md`.
