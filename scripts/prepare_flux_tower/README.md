# `prepare_flux_tower/` — MEDS forcing from flux-tower data

Builds an `ED_default` forcing file (`[forcing].format = "ED_default"`) from a flux tower's
meteorology: AmeriFlux BASE, FLUXNET/ONEFlux, or any CSV. A site TOML **declares** what the data
are. The tool **checks** each declaration against the sun and the data, converts everything to
MEDS's own conventions, fills gaps explicitly with a flag on every value, and states the tower's
heights in the file. MEDS then moves each sample from those heights to every patch's canopy-air top.
The worked example is [`examples/example_flux_tower_bci/`](../../examples/example_flux_tower_bci/);
the design record is
[`docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md`](../../docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md).

```bash
python make_tower_forcing.py --site my_site.toml --out my_forcing.nc
python compare_longwave_fill.py --site my_site.toml --out lw.json --figure lw.png
```

Needs numpy, pandas and netCDF4, plus tomli on Python < 3.11: the `meds-era5` environment of
[`../prepare_era5/environment.yml`](../prepare_era5/environment.yml) has them. The tests use pytest:
`python -m pytest scripts/prepare_flux_tower/tests`, which CTest also runs as `prepare_flux_tower`
when the Python it finds has the dependencies.

## The site TOML

```toml
[input]
format = "csv"                     # "csv" | "ameriflux_base" | "fluxnet"
path = "data/tower.csv"            # relative to this file
timestamp = "date"                 # csv: the time column (BASE/FLUXNET: TIMESTAMP_START by default)
timestamp_format = "%Y-%m-%d %H:%M"
missing = ["NaN"]                  # BASE/FLUXNET: -9999 by default

[site]
latitude = 9.1568
longitude = -79.8486
elevation = 150.0

[clock]
utc_offset = -5.0                  # [h] the data's clock (AmeriFlux: local standard time)
stamp = "begin"                    # the stamp marks the start ("begin") or end ("end") of each interval
timestep = 1800                    # [s]

[heights]                          # [m] above the ground at the tower base
tq_height = 41.0
wind_height = 41.0
pressure_height = 41.0             # 0 when the barometer is on the ground

[variables]                        # column and units of each variable
Tair   = { column = "tair",  units = "degC" }       # degC | K
RH     = { column = "RH",    units = "%" }          # % | 1
VPD    = { column = "vpd",   units = "kPa", curve = "alduchov_eskridge" }   # optional; kPa | hPa | Pa
PSurf  = { column = "p_kpa", units = "kPa" }        # kPa | hPa | Pa, at pressure_height
Rainf  = { column = "PPT",   units = "mm" }         # mm per interval | mm s-1 | kg m-2 s-1
SWdown = { column = "Rs",    units = "W m-2" }
LWdown = { column = "Rl_up", units = "W m-2" }      # optional: without it the longwave is all fill
Wind   = { column = "ubar",  units = "m s-1" }
PAR    = { column = "Par_tot", units = "umol m-2 s-1" }   # optional, reports only

[gapfill]
short_gap_max = 4                  # [records] interpolate gaps up to this long
```

`LWdown` points at `Rl_up` because the BCI file labels its two longwave columns the wrong way round
(the example's README shows how that was found). For AmeriFlux BASE and FLUXNET, point each
variable at its column, e.g. `TA_1_1_1` or `TA_F`.
With `timestamp = "TIMESTAMP_END"`, `stamp` must be `"end"`. FLUXNET's `_QC` columns mark the values
the provider filled.

## What it checks

The checks stop the build; V5 only reports.

| | Check |
|---|---|
| V1 | the time axis is uniform at `timestep`, sorted, without duplicates |
| V2 | under the declared clock and stamp, the shortwave fits the model's own sun within 10 min, and less than 0.1 % of it falls where the model sees night. This catches a local-time file declared UTC and a stamp at the wrong end of the interval, errors that keep daily totals right and scramble the sub-daily phase |
| V3 | with RH and VPD both given, the VPD is (1 − RH)·e_s(T) under the declared curve to 1 Pa. A wrong declaration is answered with the curve that fits |
| V4 | physical bounds in MEDS units. Out-of-bounds values become missing; more than 5 % of a variable is taken as a unit error |
| V5 | the JSON report beside the file: fill shares and rain by year, and the PAR/SW ratio by year (sensor drift) |

## What the file carries

- **Clock:** UTC. The source clock and stamp convention are converted, and the convention is kept
  (`avg_convention`).
- **`RHair`:** the measured relative humidity, a fraction. MEDS converts it with its own saturation
  curve. Where RH is missing, it is recovered from the provider's VPD under the declared curve.
- **`PSurf`:** at the ground, brought down from `pressure_height` hypsometrically.
- **States** (`Tair`, `RHair`, `PSurf`, `LWdown`, `Wind`): values at the stamps, as MEDS reads a
  state. Each is the mean of the two half-hour means that meet there; wind uses the root mean square.
- **Fluxes** (`Rainf`, `SWdown`): means over each interval.
- **`<Var>_qc`**, one per variable: 0 observed, 1 short gap, 3 synthesis regression or mean diurnal
  variation, 4 filled by the provider, 5 RH recovered from the provider's VPD. Code 2 is unused.
- **Global attributes:** `tq_height_m`, `wind_height_m` and `height_above = "ground"`, which MEDS
  checks against `[forcing]`, plus the fill methods and the source clock.

## Gap filling

MEDS never gap-fills, so the tool fills explicitly and flags each fill. Filling from another source,
such as ERA5-Land or a nearby station, is the user's to do in the tower file before the build.
- **Short gaps** of up to `short_gap_max` records are interpolated. Shortwave is interpolated through
  its clearness index, and wind in the energy form.
- **Long gaps** in the states and shortwave take the mean diurnal variation within ±7 days.
- **Rain** is never interpolated: a rain gap stops the build.
- **Longwave:** the model's own synthesis split into its clear-sky part εσT⁴ and its cloud part
  εσT⁴(1 − kt), regressed onto the tower by month and day/night, so the site sets the cloud
  coefficient instead of the model's 0.22. At Barro Colorado Island the fit gives 0.10, and the two
  parts beat one regression on the whole synthesis on every held-out draw. With fewer than 48
  observed longwave records there is nothing to fit, and the fill is the model's synthesis as MEDS
  computes it, with a warning.
- [`compare_longwave_fill.py`](compare_longwave_fill.py) scores the longwave fill against
  observations it never saw, beside a monthly climatology and the synthesis as MEDS computes it.

## Files

| File | What it does |
|---|---|
| `make_tower_forcing.py` | the build (the CLI) |
| `tower_inputs.py` | the site TOML and the three input formats, in MEDS units |
| `tower_checks.py` | V1–V5 |
| `tower_gapfill.py` | the fills, the longwave predictors, re-centring |
| `compare_longwave_fill.py` | the longwave comparison |
| `tests/test_tower_forcing.py` | synthetic towers from a known sun and known humidity |
| `tower_conversions.py` | the Python copy of the model's conversions (humidity, pressure, sun, longwave synthesis) |
| `../meds_forcing_file.py` | the writer of the forcing file, shared with `prepare_era5/make_forcing_file.py` |
