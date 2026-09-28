# Meteorological forcing

Forcing is the **prescribed** half of the model's boundary: the atmospheric state above the canopy,
*read from a file* rather than evolved. Everything else in MEDS is state, advanced by a kernel that
owns it; nothing in the model writes forcing. A run is one site, or a region of polygons each at its
own forcing cell (`[run].mode = "region"`); either way a polygon reads one column of the forcing. The
reader produces one **instantaneous per-site record**
(`met_forcing_t`) per `dt_fast` sub-step: air temperature, specific humidity, surface pressure, wind,
downwelling longwave, the four (beam/diffuse)×(PAR/NIR) shortwave streams, CO₂ (prescribed on its own,
not read from the met file, §12), plus the *derived* $`\cos z`$ and $`\rho_{air}`$ — never stored in a
file, always recomputed.

`libmeds_forcing` links the config leaf and the netCDF C bindings (`meds_netcdf_c`) and **nothing
else** — no state, no demography, no kernels. That placement is the point: a prescribed driver must sit
*low* in the dependency graph, because everything consumes it and it consumes nothing. It is also why
the disaggregation math is a separate `pure`/`elemental` kernel module from the reader that does I/O.

## 1. The forcing file

Two sources, chosen by `[forcing].format`: a single **MEDS forcing file** (`"ED_default"`, below) and
the global **ED_ERA5land archive** (`"ED_ERA5land"`, at the end of this section). The single file is an
**unstructured multi-grid NetCDF** whose forcing variables are all shaped `(time, grid)`.
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
    double elevation(grid) ["m"] ;              // optional, not read: the cell's orography (archive only)
    float Tair(time,grid) ["K"], PSurf(time,grid) ["Pa"],          // cell_methods = "time: point"
          Wind(time,grid) ["m s-1"] ;
    float Tdew(time,grid) ["K"] ;              // the humidity: exactly one of Tdew, RHair ["1"]
                                               // and Qair ["kg kg-1"]; see sec. 7
    float u10(time,grid), v10(time,grid) ["m s-1"] ;            // optional wind vector; see below
    float Rainf(time,grid) ["kg m-2 s-1"], LWdown(time,grid) ["W m-2"],  // cell_methods = "time: mean"
          SWdown(time,grid) ["W m-2"] ;                         // SWdown is the TOTAL; see sec. 6
    // every float carries _FillValue = 1.e20, and MEDS stops on one (it never gap-fills)
// global attributes -- time_zone is required and must be "UTC"; avg_convention, sw_input_kind and,
// when present, the heights are checked against [forcing] at open; the rest is provenance:
    :Conventions = "MEDS-forcing-1.1" ; :source = "ERA5-Land hourly" ; :time_zone = "UTC" ;
    :timestep_seconds = 3600 ; :avg_convention = "end" ; :sw_input_kind = "total" ;
    :wind_meas_height_m = 10. ;          // = [forcing].wind_height; also tq_height_m, wind_height_m
                                         // and height_above ("zero_plane" | "ground") when present
}
```

With `sw_partition = "passthrough"`, `SWdown` is replaced by four pre-split streams `SWdown_par_beam`,
`SWdown_par_diffuse`, `SWdown_nir_beam`, `SWdown_nir_diffuse` [W m⁻²], all required.

There is no CO₂ variable: CO₂ comes from `[forcing].co2_source` (§12), and a file that carries
`CO2air` is rejected at open.

**The clock is UTC.** Every forcing file is on a UTC clock and says so (`time_zone = "UTC"`); a
missing or different value stops `met_open`. A local-time file read as UTC keeps its daily totals
and moves its sun by the offset, which nothing downstream notices, so the reader does not guess.
Local time belongs to the post-processing of the output.

**The humidity is the one the source measured.** A file carries exactly one of `Tdew` (a
reanalysis), `RHair` (a flux tower, a fraction) or `Qair` (a model-made source), and the reader
converts it to specific humidity with the model's own saturation curve (§7). None, two, or an
`RHair` above 1.5 (a percentage) stop `met_open`.

**The heights, when stated, are checked.** A file may say where it was measured: `tq_height_m`,
`wind_height_m` (or `wind_meas_height_m`) and `height_above`. When present, each must agree with
`[forcing]` (heights within 0.01 m), because a disagreement would move every sample from the wrong
height (§8).

**The time axis is the metadata the reader parses:** it takes everything after `since` in
`time:units` as the base instant and the `time` values as seconds from it. Record interval, averaging
convention and shortwave scheme come from the TOML `[forcing]` block; the reader checks the file's
record spacing and its `avg_convention` and `sw_input_kind` attributes against them and stops on a
disagreement (#185).

**The wind vector.** A file that carries `u10` and `v10` supplies the vector: the record keeps both
components and derives the speed $`\sqrt{u_{10}^2+v_{10}^2}`$ from them, ignoring `Wind`. A file with
`Wind` alone gives the speed only, and `met_forcing_t%has_wind_vector` stays false.

**Binding a site to a column.** `grid_match = "explicit"` uses `forcing.grid_index` verbatim (range checked
at open). `grid_match = "nearest"` instead reads the `latitude(grid)`/`longitude(grid)` vectors and picks
the cell minimising the great-circle distance to the `[site]` coordinates (`great_circle_distance`, ED2's
`dist_gc`), with strict `<` in the scan so ties keep the lowest index; only the ordering matters, so the
Earth radius is immaterial. At `met_open` the reader loads, for that column, every record the run can
reach — the recycle window, or the run period — into memory, one read per variable, so no step reads
the file; a file written in 1 × 1 chunks is read only over that range.

**Producing one from ERA5-Land.** `scripts/prepare_era5/make_forcing_file.py` writes the file above,
either by cutting the site's cell out of an ED_ERA5land archive (`--data-path`) or from a small download:
`download_era5land_cds.py` pulls the eight hourly variables (`t2m`, `d2m`, `sp`, `u10`, `v10`, `tp`,
`ssrd`, `strd`) for a box around the site from the CDS, `postprocess_era5land.py --split none` decodes them
into one box file per variable, and `make_forcing_file.py --box-dir` converts them. Either way,
`Tair = t2m`, `PSurf = sp`, `Tdew = d2m` (the reader makes q from it by (10), as it does for the
archive), the components `u10`, `v10` with their speed
$`\mathrm{Wind}=\sqrt{u_{10}^2+v_{10}^2}`$, unfloored (the reader floors every source's speed at
0.1 m s⁻¹ itself, for the Monin–Obukhov stability), and the three **accumulated** fluxes de-accumulated
then unit-converted:
$`\mathrm{Rainf}=\Delta tp\cdot 1000/3600`$ [kg m⁻² s⁻¹], $`\mathrm{SWdown}=\Delta ssrd/3600`$,
$`\mathrm{LWdown}=\Delta strd/3600`$ [W m⁻²]. Both inputs go through one table and one rule,
`ARCHIVE_VARIABLES` and `deaccumulate` in `scripts/prepare_era5/era5land_common.py`, the ones the archive
is built with, so the file names and describes each variable as the archive does. The location is always
named (`--lat` with `--lon`, `--cells`, or `--all-cells` for box files); only the archive knows the cell's
orography, so only a file cut from it carries `elevation`.

*The 00Z trap.* ERA5-Land accumulations run from 00 UTC and reset daily, so the **00:00 stamp carries the
whole previous day's total** (step 24) — not zero, and not one hour. Ordered by valid time, the per-hour
amount is `raw(H) - raw(H-1)` everywhere *including* the 00:00 stamp (which then yields the prior day's
23Z–00Z hour), except at **01:00 UTC**, where it is `raw(01:00)` as-is (the first step of the period, with an
implicit 0 at 00:00). An off-by-one here doubles or zeroes the first hour of every day. Recovering day D's
last hour needs the 00:00 stamp of day D+1, so the downloader fetches that closing stamp with one tiny
extra request; the box files therefore start at 01:00 on the first day and end at 00:00 after the
last, and nothing is ever filled.

**Producing one from flux-tower data.** `scripts/prepare_flux_tower/make_tower_forcing.py` reads
AmeriFlux BASE, FLUXNET (ONEFlux) or a plain CSV described by a site TOML, converts the clock to UTC,
writes `RHair` as measured, fills gaps explicitly with a per-variable `<Var>_qc` flag, and states the
tower's heights; `examples/example_flux_tower_bci/` builds Barro Colorado Island's file with it. Its
conversions and checks are in `docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md` §7–§9.

**The ED_ERA5land archive (`format = "ED_ERA5land"`).** The global archive built by
`scripts/prepare_era5/` holds one file per variable per month, `ED_ERA5land_<Var>_<YYYYMM>.nc`, each
a regular `(time, lat, lon)` grid at 0.1°: `Tair`, `Tdew`, `PSurf`, `u10`, `v10` as ERA5-Land
delivers them, and `Rainf`, `SWdown`, `LWdown` already de-accumulated to hourly means, end-stamped
from 01:00 on the 1st to 00:00 on the 1st of the next month. A static file carries the `valid`
mask and the orography (`MEDS_FORCING_DESIGN.md` §14). An archive run steps daily from midnight
(`[run].dt_slow = "1d"`, `start_time` at 00:00, checked by `validate_config`), because the reader
loads one month at a time before a step. Given `data_path`, the reader:

- **binds the site to a cell** by regular-grid arithmetic; a site on a no-data cell (a coast, a
  lake edge) takes the nearest valid cell by great-circle distance, lowest index on a tie, up to
  `max_distance_km`, and otherwise stops. The cell's orography becomes the grid elevation that the
  optional lapse correction (§8) starts from;
- **lays the months end to end** as one hourly axis: the recycle window when recycling, otherwise
  the run period. Every file the run needs must exist at open, so a gap stops the run before it
  starts;
- **reads a month at a time**, one chunk column per variable, into a buffer, checking each file's
  grid, stamps and units, and rejecting a missing value in the cell;
- **converts at each stamp:** $`q`$ from the stored dewpoint and pressure by (10), and the wind
  speed from the stored components, whose vector the record also carries.

Bracketing, recycling and the seam are the same code as for the single file: the record at 00:00 on
the 1st lives in the previous month's file, and the bracket that spans it reads both.

## 2. The reader: a two-record window

The reader state comes in two parts: one `met_source_t` per run, which `met_open` fills with the
dimensions, base time and time coordinate, the validated recycle window (§9) and the records in
memory, and one `met_cursor_t` per polygon, which `met_cursor_init` binds to the polygon's cell and
location and loads with its first bracket. `met_advance(src, cur, now)` slides the cursor's bracket so
`rec_prev%when ≤ now < rec_next%when`, marching incrementally; `met_close` releases the source.
`met_instant(src, cur, now)` is a **pure function of (reader state, time)** — it interpolates and
disaggregates the loaded bracket to the instant, so it is safe to sample ahead of the threaded patch
loop. A run starting before the first record either
hard-errors (`start_clamp = "error"`, the default) or holds record #1 (`"hold"`); a non-recycling run past
the last record clamps to the final interval. The no-file **`const`** backend returns the `met_forcing_t`
defaults — a reference climate whose four shortwave streams sum to 400 W m⁻², $`\cos z`$ and
$`\rho_{air}`$ still derived — so the fast loop can run with no forcing file at all.

**No step reads a file.** Everything a step's sub-samples need is in memory before the step starts
(`MEDS_POLYGON_RUNTIME_PLAN.md` §4): a MEDS forcing file's records for the recycle window or the run
period are read at `met_open`, and for the archive `met_prefetch` loads, before each step, the month
the step reads plus the one record before it — 00:00 on the 1st, which lives in the previous month's
file, or the window's last record at the recycle wrap. Moving into a new month that record comes from
the outgoing buffer, so each archive month is read once per pass through it.

## 3. Temporal interpolation, and why wind is different

Within the bracket the weight is $`w_{next}=(t-t_{prev})/(t_{next}-t_{prev})`$, clipped to $[0,1]$. Each
variable carries its own policy: **linear** for `Tair`, the humidity (as $q$), `PSurf`, `LWdown`
(smooth atmospheric states, read as values *at* the stamps); **step-constant** for `Rainf`, held over
the interval it is the mean of, because interpolating precipitation smears intense events into
physically wrong drizzle and breaks infiltration and runoff; **cosz reconstruction** for the four
shortwave streams (§5), because those records are interval *means*, not point values. Rain and
shortwave come from the same record, the one whose interval contains the instant: `rec_next` on an
`avg_convention = "end"` file, `rec_prev` on a `"begin"` one. Wind takes the **energy form**:

```math
u(t) = \sqrt{\max\!\big[(1-w_{next})\,u_{prev}^2 + w_{next}\,u_{next}^2,\; u_{min}^2\big]},
\qquad u_{min} = 0.1\ \mathrm{m\,s^{-1}} \qquad(1)
```

Wind interpolates as a **squared** quantity, then square-rooted (ED2's `vels` convention), so it is the
kinetic-energy content, not the speed, that varies linearly across a record boundary; the floor keeps the
Monin-Obukhov solve in the aerodynamics kernel out of its degenerate zero-wind limit.

## 4. Solar geometry: apparent solar time

`meds_time%solar_cosz` expects **local apparent solar seconds** — it reads its time argument as a fraction
of day with noon at 0.5. Every forcing file is stamped in UTC, so feeding clock seconds straight in would
put solar noon at 00:00 UTC everywhere. One transform fixes it:

```math
t_{solar} = t_{UTC} + \lambda\cdot 240\ \mathrm{s\,deg^{-1}} + \mathrm{EoT}(\mathrm{doy}) \qquad(2)
```

```math
\mathrm{EoT} = 229.18\cdot 60\,\big(0.000075 + 0.001868\cos b - 0.032077\sin b
              - 0.014615\cos 2b - 0.040849\sin 2b\big)\ \mathrm{s},
\quad b = \frac{2\pi(\mathrm{doy}-1)}{365} \qquad(3)
```

Here $\lambda$ is the longitude (+E), and EoT the Spencer (1971) equation of time. The zenith cosine is
then the standard expression in declination $\delta$ and hour angle $h$, floored at zero (night):

```math
\cos z = \sin\phi\sin\delta + \cos\phi\cos\delta\cos h, \quad
\delta = 23.45^\circ\sin\!\frac{2\pi(284+\mathrm{doy})}{365}, \quad
h = 2\pi\!\left(\frac{t_{solar}}{86400}-\tfrac12\right) \qquad(4)
```

`met_instant` recomputes $`\cos z`$ every sub-step from the **model** clock; it is never read from a file
and never interpolated.

That declination is `solar_declination(doy)` (Cooper 1969), and it is the **only** one in the model:
`daylength`, which drives the photoperiod phenology cue, reads the same function. It used to carry its
own White (1997) form $`-23.44^\circ\cos\frac{2\pi(\mathrm{doy}+9)}{365}`$. The two are the same
function offset by 1.25 days — $`-\cos x = \sin(x-\pi/2)`$ puts White's ascending zero crossing at
doy 82.25 against Cooper's 81 — and Cooper's is the closer to the true vernal equinox near doy 79–80.
Unifying moved Ithaca daylength by at most **3.7 minutes** (at the equinoxes; ~0 at the solstices) and
brought the 10.5 h autumn cue **2 days earlier**.

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

A file carries the humidity its source measured, and the reader converts it to specific humidity at
each stamp, with that record's own temperature and pressure: a dewpoint (`Tdew`, a reanalysis) or a
relative humidity (`RHair`, a flux tower) by (10), a specific humidity (`Qair`) as it is. Storing the
measured quantity is what keeps the provider's saturation curve out of the model: a tower's relative
humidity is a sensor reading, while its vapour-pressure deficit, or a $q$ made from it offline, carries
whichever curve the provider used (Barro Colorado Island's `vpd` column is the Alduchov–Eskridge
form, whose saturation pressure is 5.7 Pa below Bolton's at 25 °C). With `RHair`, 100 % at the tower is
saturation in the model.

Both conversions go through the **same** Bolton (1980) saturation vapour pressure the rest of MEDS
uses (`meds_therm_lib%sat_vapor_pressure`, $`e_{sat}(T)=611.2\exp[17.67\,T_c/(T_c+243.5)]`$ Pa) — over
**liquid**, deliberately: dewpoint is *defined* as the temperature at which the liquid saturation
vapour pressure equals the actual one, and a relative-humidity sensor reports against liquid water
(the WMO convention) even below freezing, so the ice branch that `sat_vapor_pressure` grew for frozen
surfaces (`snow_biophysics.md` §1) must not be applied here. Every output file of an ED_ERA5land run
names the formula in its `forcing_qair` global attribute. The dewpoint form is the identity that the
actual vapour pressure *is* the saturation vapour pressure evaluated at the dewpoint; RH is clipped
to $[0,1]$ first. The model's own inverse, `specific_humidity_to_vpd`, is the exact inverse of (10),
so a saturated record reads back as a vapour-pressure deficit of zero.

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

## 8. Vertical corrections: the terrain, then the canopy-air top

Every vertical correction lives in one module, `meds_lapse_rate`, and runs in two steps: the **terrain**
lapse per record as it is read, then the move to **each patch's canopy-air top** per patch and sub-step.

**What the forcing's heights mean.** For ERA5-Land, which is what MEDS reads today:
- heights are above the ground plus the roughness length; the ECMWF model has no displacement height
  (IFS Cy41r2 documentation, Part IV, §3.2);
- the **10 m wind is an open-terrain diagnostic**: the model's 40 m blending-height wind brought down to
  10 m with a fixed roughness of 0.03 m, whatever the land cover (Part IV, §3.10.2). ERA5-Land takes it
  from ERA5 by interpolation;
- over a forested cell, **2 m temperature and dewpoint are a clearing's**: they come from the cell's
  dominant low-vegetation tile, not from air above the canopy;
- ERA5-Land has already lapsed ERA5 to its own orography, holding relative humidity (Muñoz-Sabater
  et al. 2021, §2.3). **Its archived surface pressure does sit at that orography**, although
  ERA5-Land's documentation lists it as interpolated. Across 2,842 neighbouring-cell pairs over the
  Colorado Rockies (1,353–3,682 m, 1 July 2024), ln p differs by 0.9 of the hypsometric rate for the
  cells' elevation difference, with a correlation of 0.985. So the terrain lapse starts pressure, like
  temperature, from the cell's own elevation. What is left is the site's elevation, and the canopy's
  height.

### Step 1. Terrain: from the forcing cell's elevation to the site's

`[site].apply_elevation_lapse`, per record at ingest. With $`\Delta z = z_{site}-z_{grid}`$ (the archive's
static orography, or `[site].grid_elevation` for a single file) and the environmental lapse rate
$\Gamma$ (`[site].lapse_rate_tair`: one value, or twelve picked by the record's calendar month, such as
Kunkel's 1989 Northern Hemisphere rates tabulated by Liston & Elder 2006, from 4.4 K/km in January to
8.2 in June), temperature and pressure move together so they stay ideal-gas consistent — the pressure
form is the hypsometric integral of $`dP/dz=-Pg/(R_dT)`$ under the *same* linear $T(z)$:

```math
T_{site} = T_{grid}-\Gamma\,\Delta z, \qquad
P_{site} = P_{grid}\left(\frac{T_{site}}{T_{grid}}\right)^{g/(R_d\Gamma)} \qquad(12)
```

with the isothermal limit $`P_{site}=P_{grid}\exp[-g\,\Delta z/(R_dT_{grid})]`$ when $`|\Gamma|\le 10^{-6}`$.
Humidity and file longwave follow, as in NLDAS (Cosgrove et al. 2003):

```math
q_{site} = q\!\left(RH\,e_s(T_{site}),\,P_{site}\right),\quad RH=\frac{e(q_{grid},P_{grid})}{e_s(T_{grid})}, \qquad
L^{\downarrow}_{site} = L^{\downarrow}_{grid}\,\frac{\varepsilon(T_{site},q_{site},P_{site})\,T_{site}^4}{\varepsilon(T_{grid},q_{grid},P_{grid})\,T_{grid}^4} \qquad(13)
```

**Relative humidity is held**, not specific humidity: holding $q$ would dry a site below its cell by
about 6 % RH per K of warming, 20 % over 500 m. A file's `RHair` is used directly at the site; a
dewpoint or a specific humidity gives its relative humidity at the cell first. $\varepsilon$ is the clear-sky emissivity of §11
(`lw_clear_form`); a synthesized longwave is built afterwards from the lapsed $T$ and $q$ instead.
Wind, shortwave and rain are unchanged, and the rain/snow split (§7) follows the lapsed temperature.
A region's polygons sit at their cells' own elevations, so there $\Delta z = 0$.

### Step 2. To the top of each patch's canopy air space

The air the canopy exchanges with is the air at the top of its canopy air space, $`z_c`$ =
`can_depth` (the tallest cohort plus a freeboard, floored; the slow loop resizes it daily). **Each patch
has its own $`z_c`$, so each patch has its own forcing.** The move is a neutral surface layer: potential
temperature and specific humidity are conserved, and the wind follows the patch's own log profile —
roughness $`z_0 = 0.13\,h`$ and displacement $d = 0.63\,h$ from `canopy_roughness`, the profile the
aerodynamics itself starts from. The forcing's heights are declared in `[forcing]`: `tq_height` and
`wind_height`, measured above the **zero plane** (`height_above = "zero_plane"`: a reanalysis, whose model
has no displacement height, which is also how CLM and JULES read reanalysis heights) or above the
**ground** (`"ground"`: a flux tower above this canopy). An open-terrain wind (`wind_exposure =
"open_terrain"`, with `wind_exposure_z0` and `wind_blending_height`) is first returned to its blending
height $`z_b`$:

```math
u(z_c) = u_m\,\frac{\ln(z_b/z_{0e})}{\ln(z_m/z_{0e})}\;\frac{\ln(h_c/z_0)}{\ln(h_b/z_0)}, \qquad
T(z_c) = T_m-\frac{g}{c_p}\,(z_c-z_T) \qquad(16)
```

with $h$ a height above $d$ (a reanalysis's heights already are; a tower's minus $d$), floored at
$`2z_0`$ like the aerodynamics' own reference height, and $`z_T = d + `$`tq_height` above the zero plane
or `tq_height` above the ground. A local wind (`"local"`) skips the first factor. Because the
aerodynamics' potential temperature is $\theta = T + (g/c_p)\,z$ with the same constants, $\theta$ at
$`z_c`$ is exactly the forcing's. Humidity, pressure, radiation, rain and CO₂ are unchanged — pressure
stays at the ground, where the canopy air, ground and leaves use it — and $`\rho_{air}`$ is re-derived.
The aerodynamics then runs from $`z_c`$: there is no fixed reference height, and nothing has to clear the
canopy.

For a 25 m canopy ($`z_c`$ = 30 m, $d$ = 15.75 m, $`z_0`$ = 3.25 m) ERA5-Land's wind becomes 0.73·u10 and
its temperature 0.12 K cooler; over a 1 m regrowth patch the wind is 0.80·u10. A flux tower above that
canopy, measuring at 41 m above the ground (Barro Colorado Island), gives `height_above = "ground"`,
`wind_exposure = "local"` and both heights 41 m: its wind becomes 0.72 of the tower's at the same
canopy-air top, and its temperature 0.11 K warmer.

**Not corrected: stability.** The neutral move leaves the stability-dependent part of the 2 m → canopy-top
difference, about ±1 K between day and night. Correcting it needs fluxes ERA5-Land does not carry, and its
2 m value describes a clearing, so inverting similarity theory from it would add error.

**In the output**, the polygon's forcing echo (`air_temp_site`, `wind_site`, …) is the forcing after the
terrain lapse, at its own heights; each patch's forcing is `wind_cas_top_patch` and
`air_temp_cas_top_patch`, with the inputs that produced it — `cas_depth_patch` ($`z_c`$), `rough_patch`
and `displace_patch`.

**What it changes** (Ithaca, the r1 cases): moving the forcing to the canopy-air top raises the friction
velocity by 33 % over an established stand in the annual mean and by 21 % over regrowth in July; the
established stand's sensible heat falls by 2.3 W/m² and its GPP rises 0.6 %. The terrain lapse (the site
sits 47.5 m below its ERA5-Land cell) adds +0.31 K, +2.2 W/m² of longwave and +0.56 kPa.

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
model year never asks for Feb-29 — so no record is double-counted, and a leap file year's Feb-29 records
go unread in the model years that have none. At the cycle edge the reader
loads a **seam bracket**: `rec_prev` is the last record inside the window, `rec_next` is the window's *first*
record, so the closing interval interpolates across the wrap instead of clamping. The window may be a
sub-range of a longer file; records past `recycle_end` belong to the next file year, not to this cycle.
The window may start at any record stamp, so the seam may fall anywhere in a day. With the archive the
reader keeps the window's first day (its first record through the next midnight) in memory from open, so
a daily step that crosses the seam reads the window's last month and that day and loads nothing mid-step.
A region's window starts at 00:00 or 01:00 on the 1st, because a region loads each month once, before it.

## 10. MEDS never gap-fills

Every required field is checked for NaN as it is read, and a missing value **halts the run**, naming the
field, the time record and the grid index. There is no gap policy, no persistence hold, no climatology fill —
and the prep script applies the same rule, erroring rather than inventing a value. This is deliberate:
gap-filling is a modelling decision with real consequences for the fluxes it produces, and it belongs
upstream in the user's own preprocessing, where it is visible and documented. Buried in a reader, it turns a
silently degraded run into one that looks clean. What is *not* a gap: the deterministic derivations — the
total→four-stream partition, humidity from dewpoint, the de-accumulation boundary sample the prep script
drops rather than writes — which compute a variable the source does not store from variables it does.

## 11. Longwave synthesis, for a source without `LWdown`

`lwdown_source = "synthesize"` reconstructs downwelling longwave from temperature and humidity:

```math
\mathrm{LW}_\downarrow = \varepsilon_{clear}\,\sigma\,T_a^4\,\bigl[1 + a\,(1-k_t)\bigr] \qquad(11)
```

`lw_clear_form` picks the clear-sky emissivity — Brutsaert (1975), a power law in the screen-level
vapour pressure, or Idso & Jackson (1969), temperature alone (the fallback when a source carries
humidity you do not trust). $`k_t`$ is the **same** clearness index the shortwave partition uses;
`lw_cloud_a` (default 0.22) is the cloud coefficient.

**Why the cloud term is not optional in practice.** Against the Ithaca ERA5-Land file, clear-sky
Brutsaert alone underestimates `strd` by a mean of **−29.9 W m⁻²** (RMSE 43.6) on a file mean of
312.3 — clear-sky formulations miss the cloud enhancement, and a temperate site is cloudy most of the
time. Setting `lw_cloud_a = 0` gives a pure clear-sky sky and that bias back.

**At night $`k_t`$ is undefined** — there is no shortwave to divide — and `clearness_index` returns a
*negative sentinel* rather than zero, because a zero would read as fully overcast and apply the
maximum cloud correction to every night. The driver holds the **last daytime** $`k_t`$ through the
dark, which is dusk's cloudiness carried forward. That one scalar is updated in `met_advance`, not
`met_instant`: `met_instant` is `intent(in)` and runs per sub-step, and the advance sits outside the
patch loop, so the state cannot become a data race when the patch axis is threaded.

**It is a fallback, not a substitute.** Driving the Ithaca run from synthesis instead of the file's
`strd` leaves the soil surface **1.37 K cooler** in the annual mean. Use the file's longwave when the
file has it; this exists so a source that lacks the field can drive MEDS at all.

## 12. CO₂

Free-atmosphere CO₂ is prescribed like the met, but **not read from the met file**.
`[forcing].co2_source` gives it, the same way for every backend and every polygon:

- **`"const"`** (the default, so existing configs run unchanged): `co2_const` for the whole run.
- **`"file"`**: a **MEDS CO₂ file**, named by `co2_file` and read once at `met_open`. The repository
  ships one (below).

The key of the other mode is rejected, not ignored. CO₂ is kept out of the met file for two reasons:

- **It is looked up on model time**, not the met file's clock. A spin-up that recycles a few years
  of met gets the CO₂ of the model year, not of the met year; §9 maps only the met into its window.
- **A met file is not a CO₂ record.** ERA5-Land has no CO₂, so the `CO2air` column the prep script
  used to write was a constant, and it silently overrode `co2_const`. A met file that carries `CO2air`
  is therefore rejected at open. Drop it with `ncks -x -v CO2air in.nc out.nc`, or rebuild the file.

**The file, format 1.** Plain text, so a user can write one by hand:

```
# '#' starts a comment that runs to the end of the line; blank lines are ignored.
timestep  1 year
units     umol/mol
1850  284.297
1851  284.420
```

- **Two keyword lines** come before the first data row, each exactly once:
  - `timestep <n> <unit>`, with `n` a whole number ≥ 1 and `unit` one of `year`, `month`, `day`,
    `hour` or `minute`;
  - `units umol/mol` (dry-air mole fraction), the only unit accepted. The line makes a file state its
    unit, which catches one written in mol/mol.
- **Each data row** is `<period start> <value>`. The start is written to the precision of the unit,
  in model time (UTC for ERA5-Land):

  | Unit | Start | Rows |
  |---|---|---|
  | year | `YYYY` | `1850`, `1851`, … (`timestep 5 year`: `1000`, `1005`, …) |
  | month | `YYYY-MM` | `1850-01`, `1850-02`, … |
  | day | `YYYY-MM-DD` | `2024-02-28`, `2024-02-29`, … |
  | hour | `YYYY-MM-DDThh` | `2024-07-01T13`, `2024-07-01T14`, … |
  | minute | `YYYY-MM-DDThh:mm` | `2024-07-01T00:00`, `2024-07-01T00:30`, … (`timestep 30 minute`) |

- **The value** is the mean CO₂ over the period, from its start $`s_i`$ to $`s_{i+1}`$, the start plus
  `n` units.
- **Each row starts exactly `n` units after the one before it.** A missing period is an error, not an
  interpolation (§10). There are at least two rows, and every value is finite and positive.

A file that breaks a rule stops the run at open, naming the line and the rule.

**How it is read.** Each value sits at the middle of its period, $`m_i = (s_i + s_{i+1})/2`$ in
calendar seconds; a leap February's middle is the 15th at 12:00. At a model instant $t$ the CO₂ is
linear between middles:

```math
c(t) = c_i + \frac{t - m_i}{m_{i+1} - m_i}\,\bigl(c_{i+1} - c_i\bigr), \qquad m_i \le t < m_{i+1} \qquad(15)
```

Before the first middle and after the last, the end value is held. Those half-periods are still inside
the period the value is the mean of, so nothing is extrapolated. The run must lie between the start of
the first period and the end of the last; one the file does not cover stops at `met_open`, not
decades in.

Interpolating moves each period's mean slightly: over an interior period of equal length, the mean of
(15) is $`(c_{i-1} + 6c_i + c_{i+1})/8`$.

- On the shipped annual series, a year's mean moves by at most 0.016 µmol/mol before 1950 and
  0.19 µmol/mol after (1997, where the growth rate jumps between years). The last year, 2022, is
  lowered by 0.27, since its second half holds the end value. Over 1850–2022 the mean moves by 0.0015.
- A monthly seasonal cycle is damped by a factor $`(6 + 2\cos 30°)/8 = 0.97`$.

**The shipped series**, `data/co2/co2_cmip7_global_annual_1000-2022.txt`:
- **Contents.** Global-mean annual CO₂ for 1000–2022, from the CMIP7 input4MIPs greenhouse-gas
  concentrations (Nicholls et al., in prep.; source_id `CR-CMIP-1-0-0`, grid `gm`). It reads 282.437 in
  1000, 277.537 in 1700, 284.297 in 1850, 313.079 in 1950 and 417.320 in 2022.
- **How it was built.** `scripts/prepare_co2/make_co2_file.py` joins the two ESGF files, 1000–1749
  and 1750–2022, which meet without a step (278.007 on both sides).
- **It is its own template.** The file's header repeats the format above.
- **Licence.** The data is CC BY 4.0, so keep the attribution lines in any copy.

Why this series:
- **It is a true global mean.** It is built from ice cores (Law Dome and others) before the
  instrumental record, and from the NOAA networks after. Mauna Loa alone reads about 1.8 µmol/mol
  above the global mean today, and a Southern Hemisphere ice core spliced onto NOAA leaves a step at
  the join.
- **One global series serves every polygon.** CMIP7's 15° latitude bands stay within 0.9 µmol/mol of
  the global mean before 1950, and within −3.1 to +3.7 in 2022, highest at northern mid-latitudes.
- **It ends in 2022.** A run past 2022 is refused until a file that reaches further is given.

**A historical run** usually spins up at a fixed pre-industrial value: `co2_source = "const"`, with
`co2_const = 277.54` for 1700 or `284.30` for 1850. The transient then runs with `co2_source =
"file"` and the shipped series. With recycled met the CO₂ still rises; `test_met_driver` checks that.

## What is not here

- **Regions run serially, and cannot restart.** `[run].mode = "region"` runs every selected cell of a box
  as its own polygon, all sharing one reader. The polygons are stepped one after another, and a region
  writes no checkpoints. The OpenMP polygon loop and region restarts are R3 and R4 of
  `MEDS_POLYGON_RUNTIME_PLAN.md` (#183).
- **No forcing perturbations.** MEDS runs on the forcing it is given. Sensitivity offsets, scalings
  and delta-change climate scenarios belong upstream, in the forcing file; they are out of scope.
- **No latitude-resolved CO₂.** One global series drives every polygon (§12). CMIP7 also gives
  monthly 15° latitude bands, which a region spanning several bands would want for recent decades.
- **`avg_convention`** takes `"end"` or `"begin"`. `"instant"` and `"center"` parse but are rejected at
  config load: they have no disaggregation branch of their own, and running the end-of-interval one under
  their name would be silently wrong (#185).

See [`docs/ROADMAP.md`](../ROADMAP.md) §8 for what is planned, and when.

## Where the code is

| Concept | Routine / file |
|---|---|
| interpolation, wind energy form | `meds_forcing_kernels`: `interpolate_forcing`, `interpolate_wind_energy` |
| apparent solar time, zenith | `meds_forcing_kernels`: `apparent_solar_seconds`, `equation_of_time`, `met_solar_cosz` over `meds_time%solar_cosz` |
| SW disaggregation | `meds_forcing_kernels`: `cosz_reconstruct_factor`, `disaggregate_sw` |
| SW partition | `meds_forcing_kernels`: `partition_shortwave`, `erbs_diffuse_fraction`, `weiss_norman_partition` |
| humidity, precip phase | `meds_forcing_kernels`: `dewpoint_to_specific_humidity`, `rh_to_specific_humidity`, `precip_phase` |
| grid match | `meds_forcing_kernels`: `great_circle_distance`, `nearest_grid_index` |
| vertical corrections (§8) | `meds_lapse_rate`: terrain `lapse_air_temperature`, `lapse_pressure`, `lapse_specific_humidity`, `lapse_longwave`, `monthly_lapse_rate` (called by `read_record`); canopy-air top `cas_top_wind_factor`, `cas_top_air_temperature`, `met_to_cas_top` (called per patch by `fast_dynamics`, with `canopy_roughness` from `meds_canopy_aerodynamics`) |
| the reader | `meds_met_driver`: `met_open`, `met_cursor_init`, `met_prefetch`, `met_advance`, `met_instant`, `met_close`; `read_record`, `assert_finite`; for the archive `open_archive`, `load_axis_month`, `locate_record`, `keep_window_head` |
| the archive's files | `meds_era5land_reader`: `era5land_path`, `era5land_select_site`, `era5land_select_box`, `era5land_load_month` |
| recycling | `meds_met_driver`: `validate_recycle_window`, `file_lookup_sec`, `recycle_model_to_file`, `load_wrap_bracket` |
| CO₂ | `meds_co2_series`: `co2_series_read`, `co2_series_at`, `co2_series_covers`; `meds_met_driver`: `open_co2`, and `met_instant` sets `met%co2`; the `CO2air` rejection in `validate_file_against_config` |
| types | `meds_forcing_types`: `met_forcing_t`, `met_record_t`, `met_source_t`, `met_cursor_t`, `met_month_t`, `co2_series_t` |
| config + selectors | `meds_forcing_config`: `forcing_config_t`, `INTERP_*`, `SWPART_*`, `LW_*`, `CLAMP_*`, `METAVG_*`, `GRIDMATCH_*`, `CO2_SOURCE_*`; validated in `meds_config`, read by `meds_config_io` |
| TOML block | `[forcing]` + `[site]` (documented in `meds_config_main.toml`) |
| fast-loop join | `meds_fast_dynamics`: per-sub-step `met_advance`/`met_instant` sampling; `fill_forcing` and `fill_aenv` take the sampled `met_forcing_t` (`reference_met` without a forcing source) |
| file production | the archive: `scripts/prepare_era5/download_era5land_gdex.py` or `download_era5land_cds.py`, then `build_era5land_archive.py`; a single file: `make_forcing_file.py`, from the archive or from `download_era5land_cds.py` and `postprocess_era5land.py` box files, or `scripts/prepare_flux_tower/make_tower_forcing.py` from flux-tower data; the shared writer and the Python mirror of the model's conversions: `scripts/forcing_common/meds_forcing_file.py`; the shipped CO₂ series: `scripts/prepare_co2/make_co2_file.py` |
| test | `test/test_met_driver.f90` — interpolation, humidity, phase, both SW schemes, the mean-conserving identity *and* the secant bias, CONST backend, NetCDF round-trip, clamp and recycle-window rejections, recycle phase over 29 years, prescribed CO₂ (format 1 at five resolutions, every rejection, CO₂ not recycled with the met, the shipped series); `test/test_met_era5land.f90` — a synthetic archive: templates, site and box selection (across 180°), month loads, the NaN rejection, the month seam, recycling across months, dewpoint and wind-vector conversion, static elevation, rejections at open; `test/test_met_tower.f90` — the flux-tower contract of an `ED_default` file: the three humidity forms and their rejections, the UTC requirement, stated heights, rain and shortwave from the interval containing the instant under both stamp conventions, and the tower round trip (RH, VPD = 0 at saturation, the move to the canopy-air top) |

## References
- Erbs, Klein & Duffie (1982), *Solar Energy* 28:293 — diffuse fraction vs the clearness index.
- Weiss & Norman (1985), *Agric. For. Meteorol.* 34:205 — band-specific direct/diffuse partitioning;
  ED2 `../ED2/ED/src/utils/radiate_utils.f90` (`short_bdown_weissnorman`).
- Bolton (1980), *Mon. Weather Rev.* 108:1046 — saturation vapour pressure.
- Brutsaert (1975), *Water Resour. Res.* 11:742 — clear-sky emissivity.
- Idso & Jackson (1969), *J. Geophys. Res.* 74:5397 — temperature-only clear-sky emissivity.
- Spencer (1971), *Search* 2:172 — Fourier series for the equation of time.
- Jin et al. (1999) — precipitation phase partitioning; ED2 `ed_met_driver.f90`.
- Longo et al. (2019), *GMD* 12:4309 — ED-2.2 technical description (the meteorological driver).
- Muñoz-Sabater et al. (2021), *ESSD* 13:4349 — ERA5-Land; §2.3, its own lapse to the ERA5-Land orography.
- ECMWF (2016), *IFS Documentation Cy41r2, Part IV: Physical Processes*, §3.2 and §3.10 — the surface-layer
  heights and the open-terrain 10 m wind.
- Cosgrove et al. (2003), *J. Geophys. Res.* 108(D22):8842 — NLDAS elevation adjustment of T, P, q and LW.
- Kunkel (1989), *J. Climate* 2:656 — monthly temperature lapse rates; tabulated in Liston & Elder (2006),
  *J. Hydrometeor.* 7:217, Table 1 (MicroMet).
- Nicholls, Meinshausen, Lewis, Pflüger, Menking et al. (in prep., 2025) — CMIP7 greenhouse-gas
  concentrations, input4MIPs `CR-CMIP-1-0-0`, doi:10.5281/zenodo.14892947 (CC BY 4.0).
- Meinshausen et al. (2017), *GMD* 10:2057 — the CMIP6 greenhouse-gas concentrations CMIP7 succeeds.
- Design doc (archived): `docs/dev_plans/archive/MEDS_FORCING_DESIGN.md`.
