# MEDS flux-tower forcing plan

**Status:** written 2026-09-28, on branch `feat/flux-tower-forcing`. Phases P0–P5 below are open.
The worked example is Barro Colorado Island (BCI), Panama.

**Goal:** turn standard flux-tower meteorology (AmeriFlux BASE, FLUXNET/ONEFlux, or a plain CSV)
into a MEDS forcing file whose conversions agree with the model's own, and declare the tower's
measurement heights so the per-patch move to the canopy-air top (`forcing.md` §8, step 2) works from
tower data. The ERA5-Land variable conventions are not followed where a tower measures something
else; the model's saturation curve, pressure convention and solar geometry are.

## 1. Decisions

| # | Decision |
|---|---|
| D1 | **MEDS runs in UTC only.** `[site].utc_offset` and `[site].apply_solar_longitude` are removed, and a `ED_default` file must say `time_zone = "UTC"`. A shift to local time is done in post-processing of the output. |
| D2 | **The file carries the humidity the source measured, and the reader converts it** with the model's own saturation curve: `RHair` for a tower, `Tdew` for ERA5-Land (as the ED_ERA5land archive already does), `Qair` for a model-made source. Exactly one per file. |
| D3 | **Format names:** `[forcing].format = "ED_default"` (the single-file format, and the default) and `"ED_ERA5land"` (the archive). `"netcdf"` and `"era5land"` are refused with a message naming the new value. Names are matched exactly. |
| D4 | **Half-hour means are re-centred to point values at the stamps** for the state variables, because MEDS interpolates states as instants (§5.2). |
| D5 | **LW↓ gaps are filled two ways**, ERA5-Land regression and MEDS-synthesis regression, and the two are compared (§8). |
| D6 | **BCI heights:** temperature/humidity, wind and the barometer all at 41 m above the ground. |
| D7 | **No data in the repository.** The BCI CSV and every forcing file built from it stay out of git; the example downloads the CSV (CC0) from Zenodo. |
| D8 | **Recycle window** 2012-08-01 → 2017-08-01 (UTC). It includes the 2015–16 El Niño drought; the example README says so. |

## 2. What the BCI file shows (measured 2026-09-28)

`BCI_v5.1.csv` (Detto, Zenodo 6456527 / Dryad doi:10.5061/dryad.3tx95x6j5, CC0): 90,528
half-hours, 2012-07-03 00:00 → 2017-08-31 23:30, regular, no duplicates.

| Item | Result | Consequence |
|---|---|---|
| Clock and stamp | Local standard time (UTC−5); `date` stamps the **start** of each interval. Fitting the 95th-percentile Rs envelope against MEDS's window-mean cos z: RMSE 36 W m⁻² start-stamped, 61 end-stamped; best shift within 5 min of start-stamping in every year 2013–17 | convert to UTC, keep `avg_convention = "begin"` |
| `vpd` | Computed by the provider with the Alduchov–Eskridge Magnus curve (610.94, 17.625, 243.04): residual 0.000 Pa. Against Bolton the residual reaches 5.2 Pa | build humidity from RH, never from `vpd` (D2) |
| Complete columns | tair, RH, vpd, p_kpa, PPT, Rs: 0 % missing (filled upstream by an ANN, no flag) | provenance only |
| LW↓ | 61 % missing: none 2012–2014, 46 % of 2015, 88 % of 2016, 100 % of 2017 | long-gap fill (§8) |
| Wind (`ubar`) | 14 % NaN plus 212 negative values | negatives are QC failures |
| PAR | 2014–2016 only; daytime diffuse ≤ total; diffuse fraction 0.98 (kt<0.3) → 0.33 (kt>0.7); Par_tot/Rs 2.08–2.11 µmol J⁻¹ every year | evaluation of the SW partition, not forcing |
| MEDS LW synthesis | Brutsaert + cloud term (a = 0.22) vs observed: bias −12.7, RMSE 30.1 W m⁻²; clear-sky alone −65.5; refitting `a` only reaches −9.7 | synthesis alone cannot fill three years |
| PPT | 1,423–2,491 mm yr⁻¹ (2013–16) | report only |
| Pressure | mean 98.83 kPa ⇒ 182–208 m ASL for a sea-level pressure of 1009–1012 hPa; the ground is 140–150 m | the barometer is at the tower top (D6) |

## 3. BCI tower metadata

| Item | Value | Source |
|---|---|---|
| Network ID | AmeriFlux PA-Bar; BASE 2012–2017, CC-BY-4.0 | ameriflux.lbl.gov/sites/siteinfo/PA-Bar |
| Location | 9.1568 N, −79.8486 E; 150 m (papers: "about 140 m ASL", top plateau) | AmeriFlux; Detto et al. |
| Tower, 2012–2017 | 41 m "AVA tower"; replaced in 2018 by a 45 m tower | Detto et al.; STRI |
| EC height | 41 m, CSAT3 sonic + LI-7500A | JGR Atmospheres (OSTI 1811664) |
| Other sensors | HC2S3 T/RH, CMP11, CNR1, BF5 PAR, TB4 gauge; 5-min logging | Detto et al. (NSF PAR 10293706) |
| Canopy height | 24.9 m in the JGR analysis (d = 0.67 h, z₀ = 0.123 h) | JGR |
| T/RH and barometer heights | not published; in the PA-Bar BADM | assumed 41 m (D6) |

## 4. Findings in MEDS that this plan fixes

1. **Rain lags one record under `"end"`.** `met_instant` holds `rec_prev%rainf` for every
   convention, while shortwave takes `rec_next` under `"end"`. Every ERA5-Land run gets its rain one
   hour late; a begin-stamped file is read correctly.
2. **States are read as instants.** Tair, the humidity, PSurf, LWdown and Wind interpolate linearly
   between stamps as point values. A tower's half-hour means written as they are would lag (end) or
   lead (begin) the shortwave by 15 min. The tool re-centres them (D4); the reader is unchanged.
3. **Two humidity constants.** `rh_to_specific_humidity` uses 0.622/0.378;
   `specific_humidity_to_vpd` uses `EPS_MOL` = 0.621987. The round trip is off by 0.06 Pa at
   e = 3 kPa.

## 5. P0 — the MEDS forcing-format changes

One config-breaking series on this branch, with a CHANGELOG "Upgrading" paragraph.

1. **Format names** (D3). `req_met_backend` reads `"ED_default"` | `"ED_ERA5land"` | `"const"`;
   the old spellings stop with a message naming the replacement. The backend constants become
   `MET_BACKEND_ED_DEFAULT` and `MET_BACKEND_ED_ERA5LAND`.
2. **UTC only** (D1). `[site].utc_offset` and `[site].apply_solar_longitude` are refused. The solar
   kernels lose the offset and the switch: $`t_{solar} = t_{UTC} + 240\,\lambda + \mathrm{EoT}`$. An
   `ED_default` file must carry `time_zone = "UTC"`; a missing or different value stops `met_open`.
3. **Humidity** (D2). An `ED_default` file carries exactly one of `RHair` [1], `Tdew` [K] or `Qair`
   [kg kg⁻¹]; none or more than one stops `met_open`. `RHair` above 1.5 anywhere stops it too (a
   percentage written as a fraction). Each record is converted to specific humidity at ingest, with
   the record's own temperature and pressure, after the terrain lapse; the lapse holds RH, so a
   `RHair` or `Tdew` source lapses without first recovering RH from q.
4. **Rain from the interval-mean record** (§4 item 1): the same `mean_rec` shortwave uses.
5. **The file's heights are checked against the config** when present: `tq_height_m`,
   `wind_height_m` (within 0.01 m) and `height_above`. The existing provenance attribute
   `wind_meas_height_m` is checked the same way.
6. **One molar-mass ratio** (§4 item 3): `specific_humidity_to_vpd` uses 0.622/0.378, like every
   forward conversion and `sat_specific_humidity`, and `EPS_MOL`, whose only user it was, is gone.
   Unifying on the forward form keeps the model state bit-identical; only the VPD diagnostics move.
7. **`Conventions = "MEDS-forcing-1.1"`** in the files the scripts write. The reader does not
   require it.
8. **The terrain-lapse keys only when the lapse is on.** `[site].lapse_rate_tair` and
   `[site].grid_elevation` are required when `apply_elevation_lapse = true` and refused when it is
   false (a key that parses and does nothing is worse than an absent one).

**Consequence of D1 for users.** 00:00 UTC is 19:00 at BCI, so the daily slow step and daily
output means run 19:00 → 19:00 local. Sub-daily output is shifted in post-processing; a daily
comparison aggregates the tower on UTC days.

## 6. P1 — one Python module for the forcing file

`scripts/forcing_common/meds_forcing_file.py`, imported by both preparation tools:

- the writer: a variable list, `<Var>_qc` int8 flags, the global attributes of §5;
- the model's math, mirrored where a tool needs it: Bolton over liquid, the humidity conversions,
  the isothermal hypsometric pressure, the solar geometry (Cooper declination, Spencer equation of
  time, the 10-sub-sample window mean), the clearness index and `synthesize_lwdown`.

`make_forcing_file.py` imports it and writes `Tdew` instead of `Qair` (D2), with the §5 attributes.

## 7. P2 — `scripts/prepare_flux_tower/make_tower_forcing.py`

A **site TOML declares** what the data are; the tool **validates** each declaration against the
sun and the data (§9) and stops on disagreement. It declares: the input format and file, the
location, the source clock (UTC offset in hours, `stamp = "begin" | "end"`, timestep), the heights
(`tq_height`, `wind_height`, `pressure_height`, all above the ground), each variable's column and
units, the provider's VPD curve when VPD is used, and the gap-fill choices.

Input formats:

- `csv` — a column map (BCI);
- `ameriflux_base` — `TIMESTAMP_START`/`TIMESTAMP_END`, −9999 missing, `_H_V_R` qualifiers;
- `fluxnet` — ONEFlux `_F` variables, whose `_QC` maps onto the qc codes of §8.

Conversions:

| File variable | Rule |
|---|---|
| `time` | source clock → UTC; the stamp convention is kept (BCI: begin; first record 2012-07-03 05:00 UTC) |
| `Tair` [K] | °C + 273.15 |
| `RHair` [1] | % ÷ 100 as measured; the reader clips to [0, 1]. Where RH is missing and VPD is present, RH = 1 − VPD / e_provider(T) with the declared curve (checked by V3), flagged |
| `PSurf` [Pa] | kPa × 1000, from the barometer height to the ground: $`P_g = P_m \exp[g z_p/(R_d T)]`$ with $g$ = 9.80665 and $R_d$ = 287.04 (BCI, $z_p$ = 41 m: +0.46 %) |
| `Rainf` [kg m⁻² s⁻¹] | mm per interval ÷ interval seconds; negatives to 0 |
| `SWdown` [W m⁻²] | total; negatives to 0 (BCI night offsets reach −3.9) |
| `LWdown` [W m⁻²] | measured, gaps filled (§8) |
| `Wind` [m s⁻¹] | scalar speed; negatives are missing; no floor (the reader floors at 0.1) |

**Re-centring (D4).** For a state whose record k is the mean over interval k, the value written at
the record's stamp is the mean of the two intervals that meet there: for a begin-stamped file,
$`x(t_k) = \tfrac12(\bar x_{k-1} + \bar x_k)`$; for an end-stamped one, $`\tfrac12(\bar x_k + \bar x_{k+1})`$.
Wind uses the root-mean-square, matching the reader's energy-form interpolation. Rain and shortwave
stay interval means. The diurnal amplitude of a half-hourly series falls by cos²(π/48) = 0.996; the
record mean is kept.

## 8. P3 — gap filling

qc codes: 0 observed · 1 short-gap interpolation · 2 ERA5-Land regression · 3 synthesis
regression or mean diurnal variation. After filling, any missing value stops the tool (MEDS never
gap-fills, `forcing.md` §10).

1. **Short gaps** (≤ 2 h by default): linear for Tair, RH, PSurf and LWdown; energy form for wind;
   shortwave through the interpolated clearness index times the window-mean top-of-atmosphere flux.
   Rain is never interpolated.
2. **Long gaps, ERA5-Land** (`--lw-fill era5` for longwave; the same path for the other variables):
   the site's cell from box files or an ED_ERA5land archive, moved to UTC half-hours on the tower's
   stamp convention (states interpolated to the stamps, fluxes held over their hour), then a
   per-variable linear regression on the overlap, by calendar month and, for longwave, day and
   night (the FLUXNET2015 `_ERA` approach, Vuichard & Papale 2015).
3. **Long gaps, synthesis** (`--lw-fill synth`): MEDS's own Brutsaert + cloud-term longwave from the
   tower's T, humidity and clearness index (dusk's clearness held through the night, as the reader
   does), then the same regression.
4. **Long gaps, fallback for the other states**: mean diurnal variation over ±7 days (Falge et al.
   2001).

**The comparison (D5).** Offline: withhold 20 % of the observed longwave in 10-day blocks, fill it
both ways, and score bias, RMSE and the diurnal and seasonal cycles, day and night. In the model:
two forcing files identical except `LWdown`, the same runs, compared against the tower's
upwelling longwave, net radiation, H and LE. BCI's observed 2016–2017 act as a control, since the
two files agree there.

## 9. Validation gates

Errors, not warnings, except V5.

- **V1** — a uniform axis at the declared timestep, no duplicates.
- **V2** — the sun. Under the declared clock, the shift that best fits the shortwave envelope to
  MEDS's window-mean cos z lies within ±10 min, and less than 0.1 % of the shortwave falls in windows
  MEDS would zero. BCI: −5 min and 0.006 %. This is the check that catches a local-time file
  declared UTC, which keeps daily totals right and scrambles the sub-daily phase.
- **V3** — humidity. With RH and VPD both present, the VPD residual under the declared curve is
  below 1 Pa (BCI: 0.000 with Alduchov–Eskridge).
- **V4** — physical bounds per variable.
- **V5** — reports: PAR/SW per year, observed − filled longwave per year, annual rain, the
  barometer height implied by ERA5-Land pressure when ERA5-Land is given.
- **V6** — the round trip through the reader (a CTest): a tower-style `ED_default` file is read by
  `met_open`/`met_instant`, and the model's RH at the forcing temperature equals `RHair`, the rain
  and shortwave interval means are returned in the right interval, and the canopy-air-top
  temperature shift is $`-(g/c_p)(z_c - z_T)`$.

## 10. P4 — `examples/example_flux_tower_bci/`

A download script (Zenodo, pinned checksum), the site TOML, the two forcing builds, a spin-up on
the five-year cycle, an evaluation run with sub-daily output, and post-processing that shifts the
output to local time and compares it with the tower's LE, H, net radiation, upwelling longwave and
GPP (FLAG = 1 only). Figures: forcing QC, the longwave comparison, and the diffuse-PAR check of the
shortwave partition.

```toml
[site]
latitude  = 9.1568
longitude = -79.8486
elevation = 150.0
apply_elevation_lapse = false          # measured at the site
[forcing]
forcing_on = true
format = "ED_default"
path = "<user path>/bci_forcing_lw-era5.nc"
grid_index = 1
grid_match = "explicit"
tq_height   = 41.0                     # HC2S3, assumed at the EC height (D6)
wind_height = 41.0                     # CSAT3
height_above  = "ground"
wind_exposure = "local"
timestep = "1800s"
avg_convention = "begin"
sw_partition = "clearidx"
lwdown_source = "file"
co2_source = "file"
co2_file = "data/co2/co2_cmip7_global_annual_1000-2022.txt"
recycle = true
recycle_start = "2012-08-01 00:00:00"
recycle_end   = "2017-08-01 00:00:00"
```

## 11. P5 — tests

- The reader (Fortran, CTest): each new rejection — old format names, `utc_offset`, a non-UTC
  file, no or two humidity variables, `RHair` in percent, a height mismatch, lapse keys with the
  lapse off — and V6.
- The tool (Python, run by CTest when the interpreter has numpy and netCDF4): synthetic AmeriFlux
  BASE, FLUXNET and CSV inputs made inside the test from a known sun and known humidity; V2 and V3
  rejections (wrong clock, wrong stamp, wrong curve); re-centring and gap-fill flags.

## 12. Not in scope

- Stability corrections to the move to the canopy-air top (neutral, as in the model).
- Gap patches take the above-canopy tower wind through their own log profile, which is likely too
  windy for a gap inside a closed forest.
- Tower CO₂ as forcing: nocturnal build-up in the roughness sublayer makes it a poor free-atmosphere
  value. `co2_source` stays the CMIP7 series.
- Soil moisture from the tower's SWC as an initial condition.
