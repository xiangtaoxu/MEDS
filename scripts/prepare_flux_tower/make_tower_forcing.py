#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""make_tower_forcing.py -- build a MEDS forcing file ([forcing].format = "ED_default") from
flux-tower meteorology: AmeriFlux BASE, FLUXNET/ONEFlux, or any CSV a site TOML describes
(docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md sec. 7-9).

What it does, in order:
  1. reads the site TOML, which DECLARES the file, the location, the clock (UTC offset and which
     end of each interval a stamp marks), the sensor heights and every variable's column and units;
  2. converts every variable to MEDS units and screens physical bounds (V4);
  3. checks the time axis (V1), the clock against the sun (V2) and a provider VPD against RH under
     its declared saturation curve (V3), and stops on any disagreement;
  4. moves the stamps to UTC, and brings the barometer's pressure down to the ground;
  5. fills every gap explicitly -- short gaps by interpolation, long ones from ERA5-Land or the
     mean diurnal variation, longwave from ERA5-Land or the model's own synthesis, each regressed
     onto the tower -- and flags each value with <Var>_qc;
  6. turns the states' interval means into values at the stamps, which is how MEDS reads a state;
  7. writes the file, with its heights, clock and fills stated, and a JSON report beside it (V5).

The file carries relative humidity as measured (RHair); MEDS converts it with its own saturation
curve, so the provider's never enters the model. The clock of the file is UTC; shift MEDS's output
back to local time in post-processing.

Usage:
  python make_tower_forcing.py --site examples/example_flux_tower_bci/bci_site.toml \\
      --out bci_forcing_lw-synth.nc --lw-fill synth
  python make_tower_forcing.py --site bci_site.toml --out bci_forcing_lw-era5.nc --lw-fill era5 \\
      --era5-file bci_era5land.nc        # an ED_default file from prepare_era5/make_forcing_file.py

Dependencies: numpy, pandas, netCDF4 (and tomli on Python < 3.11).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tower_checks as tc    # noqa: E402
import tower_gapfill as tg   # noqa: E402
import tower_inputs as ti    # noqa: E402

mff = ti.mff
STATES = ("Tair", "RH", "PSurf", "Wind", "SWdown")
FILE_NAMES = {"Tair": "Tair", "RH": "RHair", "PSurf": "PSurf", "Wind": "Wind", "Rainf": "Rainf",
              "SWdown": "SWdown", "LWdown": "LWdown"}
RECENTRED = ("Tair", "RH", "PSurf", "LWdown", "Wind")
BOUNDS_AFTER_FILL = {"RH": (0.0, 1.05), "Wind": (0.0, None), "SWdown": (0.0, None)}


def prepare(site, lw_fill, states_fill, era5_path=None, short_gap_max=4, lw_holdout=None):
    """Steps 1-5: the tower's filled interval means on UTC stamps, their qc flags and the checks'
    reports. `lw_holdout`, a boolean mask over the records, hides those longwave observations from
    the fill, so a caller can score a fill against them (compare_longwave_fill.py)."""
    data = ti.read_tower(site)
    report = {"site": site.name, "input": os.path.basename(site.input_path), "format": site.input_format}
    report["V1_axis"] = tc.check_axis(data.values.index.values, site.timestep)
    values, report["V4_bounds"] = tc.screen_bounds(data.values)
    stamps = ti.to_utc(values.index.values, site.utc_offset)
    starts, ends = ti.interval_bounds(stamps, site.stamp, site.timestep)
    mean_cosz = mff.window_mean_cosz(starts, site.timestep, site.latitude, site.longitude)
    report["V2_sun"] = tc.check_sun(stamps, values["SWdown"].to_numpy(), site)
    n = len(stamps)

    y, qc = {}, {}
    for name in ("Tair", "RH", "PSurf", "Rainf", "SWdown", "LWdown", "Wind"):
        if name in values:
            y[name] = values[name].to_numpy(dtype=float).copy()
            qc[name] = np.where(data.provider[name].to_numpy(), tg.QC_PROVIDER, tg.QC_OBSERVED).astype(np.int8)
        else:
            y[name] = np.full(n, np.nan)
            qc[name] = np.zeros(n, dtype=np.int8)

    # the humidity: RH as measured; where it is missing, from the provider's VPD by its own curve
    report["V3_humidity"] = dict(checked=False)
    if "VPD" in values:
        curve = site.variables["VPD"]["curve"]
        vpd = values["VPD"].to_numpy(dtype=float)
        if "RH" in values:
            report["V3_humidity"] = tc.check_vpd(y["Tair"], y["RH"], vpd, curve)
        from_vpd = ~np.isfinite(y["RH"]) & np.isfinite(vpd) & np.isfinite(y["Tair"])
        y["RH"][from_vpd] = tc.rh_from_vpd(y["Tair"][from_vpd], vpd[from_vpd], curve)
        qc["RH"][from_vpd] = np.where(data.provider["VPD"].to_numpy()[from_vpd], tg.QC_PROVIDER, tg.QC_OBSERVED)
        report["V3_humidity"]["rh_from_vpd_records"] = int(from_vpd.sum())

    era5 = None
    if era5_path:
        era5_file = tg.read_era5_file(era5_path, site)
        era5 = tg.era5_on_intervals(era5_file, starts, ends, site)
        report["era5land"] = dict(file=os.path.basename(era5_path), cell=era5_file["cell"],
                                  elevation=era5_file["elevation"])
    if (lw_fill == "era5" or states_fill == "era5") and era5 is None:
        raise SystemExit("ERROR: an ERA5-Land fill needs --era5-file (or gapfill.era5_file in the site TOML)")

    fills = {}
    groups_month = tg.regression_groups(stamps, mean_cosz, by_day_night=False)
    records_per_day = int(round(86400.0 / site.timestep))
    for name in STATES:
        kind = {"Wind": "energy", "SWdown": "shortwave"}.get(name, "linear")
        y[name], qc[name] = tg.fill_short(y[name], qc[name], short_gap_max, kind, mean_cosz)
        lower, upper = BOUNDS_AFTER_FILL.get(name, (None, None))
        if not np.isfinite(y[name]).all():
            if states_fill == "era5":
                y[name], qc[name], fills[name] = tg.fill_by_regression(
                    y[name], qc[name], era5[name], groups_month, tg.QC_ERA5, lower, upper)
            else:
                y[name], qc[name] = tg.fill_mean_diurnal(y[name], qc[name], records_per_day)
                fills[name] = dict(method="mean_diurnal_variation", filled=int((qc[name] == tg.QC_SYNTH_OR_MDV).sum()))
        if name == "SWdown":
            y[name] = np.where(mean_cosz > mff.COSZ_BAR_MIN, np.maximum(y[name], 0.0), 0.0)

    # rain is never interpolated or averaged: a gap takes ERA5-Land's rain, or stops the build
    missing_rain = ~np.isfinite(y["Rainf"])
    if missing_rain.any():
        if era5 is None:
            raise SystemExit(f"ERROR: Rainf has {int(missing_rain.sum())} missing records and there is no ERA5-Land "
                             f"file to fill them from; MEDS does not invent rain.")
        y["Rainf"][missing_rain] = era5["Rainf"][missing_rain]
        qc["Rainf"][missing_rain] = tg.QC_ERA5
        fills["Rainf"] = dict(method="era5land", filled=int(missing_rain.sum()))

    # the pressure at the ground, where MEDS keeps it (docs/science/forcing.md sec. 8)
    if era5 is not None and era5_file["elevation"] is not None:
        ok = (qc["PSurf"] == tg.QC_OBSERVED) & np.isfinite(era5["PSurf"])
        if ok.sum() > 100:
            dz = -mff.R_DRY * np.nanmean(y["Tair"][ok]) / mff.GRAV * np.mean(np.log(y["PSurf"][ok] / era5["PSurf"][ok]))
            report["V5_barometer_height_above_ground_m"] = round(float(dz + era5_file["elevation"] - site.elevation), 1)
    y["PSurf"] = mff.pressure_at_height(y["PSurf"], y["Tair"], -site.pressure_height)

    # the longwave, last, because the synthesis needs the filled states
    if lw_holdout is not None:
        y["LWdown"][lw_holdout] = np.nan
    y["LWdown"], qc["LWdown"] = tg.fill_short(y["LWdown"], qc["LWdown"], short_gap_max)
    groups_lw = tg.regression_groups(stamps, mean_cosz, by_day_night=True)
    if lw_fill == "era5":
        predictor, code = era5["LWdown"], tg.QC_ERA5
    else:
        predictor = tg.synthesis_predictors(y["Tair"], y["RH"], y["PSurf"], y["SWdown"], mean_cosz)
        code = tg.QC_SYNTH_OR_MDV
    y["LWdown"], qc["LWdown"], fills["LWdown"] = tg.fill_by_regression(y["LWdown"], qc["LWdown"], predictor,
                                                                        groups_lw, code, 30.0, 650.0)
    fills["LWdown"]["method"] = lw_fill
    report["fills"] = fills
    return dict(stamps=stamps, starts=starts, mean_cosz=mean_cosz, values=y, qc=qc, report=report,
                source=values)


def build(site, out, lw_fill, states_fill, era5_path=None, short_gap_max=4, start=None, end=None,
          report_path=None):
    """Steps 1-7."""
    p = prepare(site, lw_fill, states_fill, era5_path, short_gap_max)
    y, qc, stamps = p["values"], p["qc"], p["stamps"]
    for name in RECENTRED:
        y[name], qc[name] = tg.recentre(y[name], qc[name], site.stamp, energy=(name == "Wind"))
    keep = np.ones(len(stamps), dtype=bool)
    if start:
        keep &= stamps >= np.datetime64(start, "s")
    if end:
        keep &= stamps < np.datetime64(end, "s")
    if not keep.any():
        raise SystemExit("ERROR: --start/--end leave no records")
    stamps = stamps[keep]
    arrays = {FILE_NAMES[name]: y[name][keep][:, None] for name in FILE_NAMES}
    flags = {FILE_NAMES[name]: qc[name][keep][:, None] for name in FILE_NAMES}
    humidity_source = "RH as measured, as a fraction"
    if "VPD" in site.variables:
        humidity_source += (f"; where RH is missing, from the provider's VPD under its "
                            f"{site.variables['VPD']['curve']} saturation curve")
    attrs = dict(
        title="MEDS meteorological forcing from flux-tower data",
        source=f"{site.name}: {os.path.basename(site.input_path)} ({site.input_format})",
        history="make_tower_forcing.py (scripts/prepare_flux_tower)",
        timestep_seconds=int(round(site.timestep)),
        avg_convention=site.stamp,
        sw_input_kind="total",
        tq_height_m=site.tq_height,
        wind_height_m=site.wind_height,
        height_above="ground",
        pressure_height_m=site.pressure_height,
        pressure_note="PSurf is at the ground, brought down from the barometer height hypsometrically",
        source_utc_offset_hours=site.utc_offset,
        state_sampling=("Tair, RHair, PSurf, LWdown and Wind are values at the stamps, the mean of the two "
                        "intervals that meet there; Rainf and SWdown are means over each record's interval"),
        humidity_source=humidity_source,
        gapfill_longwave=lw_fill,
        gapfill_states=states_fill,
        gapfill_short_max_records=int(short_gap_max),
    )
    mff.write_forcing_file(out, stamps, [(site.latitude, site.longitude, site.elevation)], arrays, attrs, flags)
    report = p["report"]
    report["V5_years"] = tc.yearly_report(stamps, {k: v[:, 0] for k, v in arrays.items()},
                                          {k: v[:, 0] for k, v in flags.items()},
                                          {k: p["source"][k].to_numpy()[keep] for k in p["source"].columns})
    report["output"] = dict(file=os.path.basename(out), first=str(stamps[0]), last=str(stamps[-1]),
                            records=int(len(stamps)))
    report_path = report_path or os.path.splitext(out)[0] + ".report.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=1, default=float)
    print(f"wrote {out}: {len(stamps)} records, {stamps[0]} .. {stamps[-1]} UTC, avg_convention = {site.stamp!r}")
    print(f"  V2: best shift {report['V2_sun']['best_shift_min']:+.0f} min, night shortwave "
          f"{100 * report['V2_sun']['night_shortwave_fraction']:.3f}%")
    if report["V3_humidity"].get("checked"):
        print(f"  V3: VPD under {report['V3_humidity']['declared']} to {report['V3_humidity']['residual_p99_pa']:.3f} Pa")
    for name, f in report["fills"].items():
        print(f"  fill {name}: {f}")
    print(f"  report: {report_path}")
    return report


def parse_args(argv):
    ap = argparse.ArgumentParser(description="Flux-tower meteorology -> a MEDS ED_default forcing file.")
    ap.add_argument("--site", required=True, help="the site TOML that declares the data")
    ap.add_argument("--out", required=True, help="the forcing NetCDF to write")
    ap.add_argument("--lw-fill", choices=("synth", "era5"), help="longwave gaps: the model's synthesis or ERA5-Land, "
                    "regressed onto the tower (default: gapfill.longwave)")
    ap.add_argument("--states-fill", choices=("mdv", "era5"), help="long gaps in T, RH, P, wind, SW: the mean "
                    "diurnal variation or ERA5-Land (default: gapfill.states, else mdv)")
    ap.add_argument("--era5-file", help="an ED_default file cut from ERA5-Land for the site "
                    "(prepare_era5/make_forcing_file.py; default: gapfill.era5_file)")
    ap.add_argument("--short-gap-max", type=int, help="the longest gap interpolated, in records (default 4)")
    ap.add_argument("--start", help="first UTC stamp to write, e.g. 2012-08-01")
    ap.add_argument("--end", help="first UTC stamp NOT written")
    ap.add_argument("--report", help="the JSON report (default: <out>.report.json)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    site = ti.read_site(args.site)
    g = site.gapfill
    lw_fill = args.lw_fill or g.get("longwave")
    if lw_fill not in ("synth", "era5"):
        raise SystemExit("ERROR: choose the longwave fill: --lw-fill synth|era5 or gapfill.longwave")
    states_fill = args.states_fill or g.get("states", "mdv")
    era5_path = args.era5_file or g.get("era5_file")
    if era5_path and not os.path.isabs(era5_path) and not args.era5_file:
        era5_path = os.path.join(os.path.dirname(site.toml_path), era5_path)
    if era5_path and not os.path.exists(era5_path):
        if lw_fill == "era5" or states_fill == "era5":
            raise SystemExit(f"ERROR: the ERA5-Land file {era5_path} does not exist")
        era5_path = None
    short = args.short_gap_max if args.short_gap_max is not None else int(g.get("short_gap_max", 4))
    build(site, args.out, lw_fill, states_fill, era5_path, short, args.start, args.end, args.report)


if __name__ == "__main__":
    main()
