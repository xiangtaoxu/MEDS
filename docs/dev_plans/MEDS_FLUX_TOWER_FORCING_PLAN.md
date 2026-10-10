# MEDS flux-tower forcing plan

**Status:** written 2026-09-28; P0–P5 merged in #320, with the soil-water fix of §13 that the BCI
example needed. Revisited 2026-09-29 on branch `feat/bci-lw-fill`: the BCI file's two longwave
columns were found swapped (D9, §2), and the ERA5-Land fill, until then tested on synthetic
ERA5-Land only, was scored against the ED_ERA5land archive, found no better than the synthesis
regression at BCI, and removed (D5, §8).

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
| D5 | **LW↓ gaps are filled by the MEDS-synthesis regression, and by nothing else** (§8). An ERA5-Land regression was built beside it and scored equal at BCI; it was removed, with ERA5-Land filling of the other variables, because filling from another source is the user's to do in the tower file before the build. |
| D6 | **BCI heights:** temperature/humidity, wind and the barometer all at 41 m above the ground. |
| D7 | **No data in the repository.** The BCI CSV and every forcing file built from it stay out of git; the example downloads the CSV (CC0) from Zenodo. |
| D8 | **Recycle window** 2012-08-01 → 2017-08-01 (UTC). It includes the 2015–16 El Niño drought; the example README says so. |
| D9 | **BCI's longwave is the column `Rl_up`.** `BCI_v5.1.csv` labels its downwelling and upwelling longwave the wrong way round (§2), so the site TOML declares `LWdown = Rl_up`, and it and the example README say why. The tool does not test for a swap: no one pattern holds across towers. |

## 2. What the BCI file shows (measured 2026-09-28)

`BCI_v5.1.csv` (Detto, Zenodo 6456527 / Dryad doi:10.5061/dryad.3tx95x6j5, CC0): 90,528
half-hours, 2012-07-03 00:00 → 2017-08-31 23:30, regular, no duplicates.

| Item | Result | Consequence |
|---|---|---|
| Clock and stamp | Local standard time (UTC−5); `date` stamps the **start** of each interval. Fitting the 95th-percentile Rs envelope against MEDS's window-mean cos z: RMSE 36 W m⁻² start-stamped, 61 end-stamped; best shift within 5 min of start-stamping in every year 2013–17 | convert to UTC, keep `avg_convention = "begin"` |
| `vpd` | Computed by the provider with the Alduchov–Eskridge Magnus curve (610.94, 17.625, 243.04): residual 0.000 Pa. Against Bolton the residual reaches 5.2 Pa | build humidity from RH, never from `vpd` (D2) |
| Complete columns | tair, RH, vpd, p_kpa, PPT, Rs: 0 % missing (filled upstream by an ANN, no flag) | provenance only |
| LW↓ | 61 % missing: observed in none of 2012–2014, 46 % of 2015, 88 % of 2016, all of 2017 | long-gap fill (§8) |
| LW columns | **swapped** (D9). The provider's `Rnet` = Rs − Rs_dn + **Rl_up − Rl_dn** (RMS residual 1.7 W m⁻² after V4; 87.6 as labelled). At night `Rl_dn`/σT⁴ = 1.022, more than the air can emit, and `Rl_up`/σT⁴ = 0.948. ERA5-Land's daily longwave: r 0.86 with `Rl_up`, 0.12 with `Rl_dn`; bias −10 and −48 W m⁻² | `LWdown = Rl_up` (D9) |
| Wind (`ubar`) | 14 % NaN plus 212 negative values | negatives are QC failures |
| PAR | 2014–2016 only; daytime diffuse ≤ total; diffuse fraction 0.98 (kt<0.3) → 0.33 (kt>0.7); Par_tot/Rs 2.08–2.11 µmol J⁻¹ every year | evaluation of the SW partition, not forcing |
| MEDS LW synthesis | Brutsaert + cloud term (a = 0.22) vs observed (`Rl_up`): bias +21.5, RMSE 30.6 W m⁻², r 0.48; clear-sky alone −27.8; refitting `a` alone gives a = 0.12, bias −1.4, RMSE 17.3 | synthesis alone cannot fill three years; regress it (§8) |
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

`scripts/forcing_common/meds_forcing_file.py`, imported by both preparation tools (since 2026-09-30 the
writer is `scripts/meds_forcing_file.py` and the conversions are
`scripts/prepare_flux_tower/tower_conversions.py`; `MEDS_EFFICIENCY_SWEEP_PLAN.md` §5.1):

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

qc codes: 0 observed · 1 short-gap interpolation · 3 synthesis regression or mean diurnal
variation · 4 filled by the data provider (FLUXNET `_QC` > 0) · 5
relative humidity recovered from the provider's VPD through its declared curve. Code 2 is unused.
With too few
observed longwave records to fit (fewer than 48), the synthesis fill is the model's own synthesis,
cloud coefficient 0.22, with a warning. After
filling, any missing value stops the tool (MEDS never gap-fills, `forcing.md` §10).

1. **Short gaps** (≤ 2 h by default): linear for Tair, RH, PSurf and LWdown; energy form for wind;
   shortwave through the interpolated clearness index times the window-mean top-of-atmosphere flux.
   Rain is never interpolated, and a rain gap stops the build.
2. **Long gaps in the longwave, synthesis:** MEDS's own Brutsaert + cloud-term longwave from the
   tower's T, humidity and clearness index (dusk's clearness held through the night, as the reader
   does), regressed in its two parts, the clear-sky emission εσT⁴ and the cloud term εσT⁴(1 − kt), so
   the cloud coefficient is fitted at the site rather than taken as 0.22. At BCI the pooled fit
   gives 0.10 per unit of clear-sky emission, and the two parts beat a single regression on the
   whole synthesis on every held-out draw (13.7–14.5 against 14.5–16.7 W m⁻²). (Before D9 this
   paragraph reported that the observed longwave fell with daytime cloudiness; that was the
   canopy's emission, read from the mislabelled column.)
3. **Long gaps in the other states**: mean diurnal variation over ±7 days (Falge et al. 2001).

**The comparison (D5).** Offline: withhold 20 % of the observed longwave in 10-day blocks, fill it,
and score bias, RMSE and the diurnal and seasonal cycles, day and night. The ERA5-Land row is from
before its removal. At BCI, on 7,200
hidden half hours (seed 1; seeds 2 and 3 in brackets):

| fill | bias | RMSE | r | diurnal-cycle RMSE |
|---|---|---|---|---|
| synthesis regression | +2.3 | 14.5 (13.7, 13.9) | 0.79 | 4.3 |
| ERA5-Land regression (removed) | +2.4 | 14.8 (13.6, 13.8) | 0.77 | 3.9 |
| monthly day/night climatology | +4.4 | 18.6 (16.0, 16.3) | 0.63 | 7.1 |
| `lwdown_source = "synthesize"` as it is | +20.9 | 29.7 (30.0, 31.1) | 0.61 | 22.1 |

- **ERA5-Land**, before its removal, was the archive's cell at 9.2 N, −79.8 E, 7.2 km from the tower, 73 % land at 73 m.
  It was aligned in time with the tower (shortwave and air temperature correlate best at zero lag),
  its raw longwave was 10 W m⁻² low, and its daily means followed the tower's with r = 0.86. A fill
  on the synthesis parts and ERA5-Land together scored 13.6 against 14.3 for either alone, over
  five draws.
- **In the model:** two forcing files identical except `LWdown` (the synthesis and ERA5-Land fills),
  the same five-year census run.
  Every flux's mean agrees between them to 0.1 W m⁻², and in 2013–2014, when the tower measured
  net radiation but no longwave, the night-time net radiation is −18.9 (synthesis) and −19.7
  (ERA5-Land) against the tower's −28.4.
- **Not adopted: rebuilding LW↓ from the tower's own Rnet.** Where Rnet was measured and the
  longwave was not (48 % of 2013, 99 % of 2014, 42 % of 2015), LW↓ = Rnet − (Rs − Rs_dn) + LW↑, with
  LW↑ regressed on σT⁴ and the absorbed shortwave, scores RMSE 7.3–7.6 and r 0.94 on the same hidden
  records, half the error of either fill. It is not in the tool: the upwelling sensor's night-time
  ratio to σT⁴ drifts from 1.007 in 2015 to 1.038 in 2017, which a regression fitted on 2015–2017
  would carry into 2013–2014, and a forcing made from the tower's Rnet makes the model's net
  radiation partly circular against it.

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
- **V5** — reports: PAR/SW per year, observed − filled longwave per year, annual rain.
- **V6** — the round trip through the reader (a CTest): a tower-style `ED_default` file is read by
  `met_open`/`met_instant`, and the model's RH at the forcing temperature equals `RHair`, the rain
  and shortwave interval means are returned in the right interval, and the canopy-air-top
  temperature shift is $`-(g/c_p)(z_c - z_T)`$.

## 10. P4 — `examples/example04_column_biophysics/`

A download script (Zenodo, pinned checksum), the site TOML, the forcing build, a five-year run
from the 2010 census with sub-daily output (`MEDS_BCI_CENSUS_INIT_PLAN.md`; it replaced the
spin-up), and post-processing that shifts the
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
path = "data/bci_forcing.nc"
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
- The tool (Python, run by CTest as `prepare_flux_tower` when the interpreter has numpy, pandas,
  netCDF4 and pytest): synthetic AmeriFlux
  BASE, FLUXNET and CSV inputs made inside the test from a known sun and known humidity; V2 and V3
  rejections (wrong clock, wrong stamp, wrong curve); re-centring and gap-fill flags.

## 12. Not in scope

- Stability corrections to the move to the canopy-air top (neutral, as in the model).
- Gap patches take the above-canopy tower wind through their own log profile, which is likely too
  windy for a gap inside a closed forest.
- Tower CO₂ as forcing: nocturnal build-up in the roughness sublayer makes it a poor free-atmosphere
  value. `co2_source` stays the CMIP7 series.
- Soil moisture from the tower's SWC as an initial condition.

## 13. Found while running the example: dry soil does not re-wet

The BCI spin-up does not establish a stand. After 50 years from bare ground, LAI is 0.009 and AGB is
0.002 kgC m⁻². Seed-rain recruits never grow. On the same branch, the Ithaca spin-up takes off
around year 16, as its README describes, so the forcing changes are not the cause.

**What a three-year diagnostic run shows.** It had daily soil, water and carbon output:

| month | rain (mm) | ET (mm) | top-layer θ | top-layer ψ (MPa) |
|---|---|---|---|---|
| 1962-12 | 365 | 49 | 0.308 | −0.00 |
| 1963-01 | 2 | 47 | 0.128 | −6.2 |
| 1963-02 | 52 | 30 | 0.080 | −37.6 |
| 1963-05 | 180 | 24 | 0.081 | −25.6 |
| 1963-07 | 292 | 22 | 0.081 | −13.1 |

The first dry season dries the top layer to θ = 0.08. The next wet season's rain then runs off
without re-wetting it, and the stand stays under water stress all year: leaf water potential sits
at −3 to −70 MPa and the stomatal factor at 0.00–0.03.

**The mechanism.** `soil_water_step` limited infiltration by the top layer's own conductivity,
`q_inf_max = K(θ₁)·(1 + ψ₁/z₁)`. For the default van Genuchten loam, K at θ = 0.08 (just above
θ_res = 0.078) is 3×10⁻¹⁵ of K_sat. The capacity was therefore effectively zero, and Horton
overflow took the rain. Ithaca never dries that far, which is why no run found this.

**What the reference models do** (checked 2026-09-28):
- **ED2**, by default, has no conductivity limit at the surface. Water on the surface moves into the
  top soil layer until that layer's pore space is full (`rk4_misc.f90`). Its conductivity-limited
  alternative uses the top layer's own K, as MEDS did, and is disabled with `fatal_error`. Between
  layers, ED2 interpolates K log-linearly (`rk4_derivs.f90`).
- **CLM5** caps infiltration at (1 − f_sat)·Θ_ice·K_sat, independent of the top layer's moisture,
  and evaluates interior conductivity at the two layers' mean moisture.
- **FATES** leaves soil water to its host model (CLM or ELM).

**The options measured**, in isolated worktrees on the three-year BCI diagnostic:

| surface rule | top-layer θ, May 1963 | top-layer θ, July 1963 | 3-year ET / rain |
|---|---|---|---|
| top layer's K (old) | 0.081 | 0.081 | 16 % |
| geometric mean of K(θ₁) and K_sat | 0.254 | 0.326 | 26 % |
| K_sat floor | 0.254 | 0.318 | 26 % |

At θ = 0.08 the geometric mean's capacity is only 0.09 mm/h. It still recovers, because the first
~0.7 mm raises the top layer to θ = 0.10, where capacity is 1.8 mm/h, and the 5 mm pond holds the
storm's water meanwhile. The old rule's 5×10⁻⁹ mm/h never gets there. The integral (Kirchhoff) mean
across the half layer, the literature's accurate choice for steep fronts (Zaidel & Russo 1992),
gives about 39 mm/h at every dryness; it was not needed to fix BCI.

**Decision (the user, 2026-09-28): ED2's geometric rule on every soil-water face.**
- **Between layers:** log-linear interpolation of K to the face,
  K_face = K_k^(1−w)·K_(k+1)^w with w = dz(k)/(dz(k)+dz(k+1)). It replaces the upstream pick.
- **At the surface:** the geometric mean of K(θ₁) and K_sat, the same rule with the pond as thick
  as the top layer.
- **The aquifer bottom boundary stays upstream-weighted.** It is a boundary with the saturated zone,
  and changing it would change drainage and capillary rise.

**The cost it brings:** capillary rise into a dry profile is slower. In `test_column_hydrology`'s
aquifer case, the residual bottom flux passes 10⁻⁵ kg m⁻² s⁻¹ at about 1,000 h instead of about 450 h.
