# SPDX-License-Identifier: Apache-2.0
"""tower_checks.py -- the validation gates of a tower forcing build (MEDS_FLUX_TOWER_FORCING_PLAN.md
sec. 9). Each gate checks a DECLARATION of the site TOML against the data and stops the build on a
disagreement; V5 only reports.

  V1  the time axis is uniform at the declared timestep, sorted, without duplicates
  V2  the sun: under the declared clock and stamp convention, the shortwave envelope lines up with
      the model's own window-mean cos z, and almost no shortwave falls where the model sees night
  V3  the humidity: a provider VPD agrees with RH under the declared saturation curve
  V4  physical bounds, after unit conversion (a wrong unit shows up here)
  V5  reports: coverage and fill by year, annual rain, PAR/SW by year, the longwave fills
"""
import numpy as np
import pandas as pd

import tower_inputs as ti

V2_MAX_SHIFT_MIN = 10.0        # [min] the envelope's best-fit shift may be at most this
V2_MAX_NIGHT_SW = 1.0e-3       # the shortwave fraction the model would zero may be at most this
V3_MAX_RESIDUAL_PA = 1.0       # [Pa] 99th percentile of |VPD - (1 - RH) e_s(T)|
V4_MAX_BAD_FRACTION = 0.05     # more than this out of bounds is a unit error, not a few spikes

BOUNDS = {                      # MEDS units
    "Tair": (170.0, 340.0), "RH": (0.0, 1.05), "VPD": (-50.0, 1.0e4), "PSurf": (3.0e4, 1.1e5),
    "Rainf": (0.0, 0.2), "SWdown": (0.0, 1500.0), "LWdown": (30.0, 650.0), "Wind": (0.0, 75.0),
    "PAR": (0.0, 3000.0),
}


def check_axis(stamps, timestep):
    """V1."""
    stamps = np.asarray(stamps, dtype="datetime64[s]")
    if len(stamps) < 2:
        raise SystemExit("ERROR (V1): the tower file has fewer than two records")
    step = np.diff(stamps).astype("timedelta64[s]").astype(float)
    if (step <= 0).any():
        k = int(np.argmax(step <= 0))
        raise SystemExit(f"ERROR (V1): the time axis is not strictly increasing at record {k + 2} "
                         f"({stamps[k + 1]}); duplicates or disorder")
    bad = np.abs(step - timestep) > 0.5
    if bad.any():
        k = int(np.argmax(bad))
        raise SystemExit(f"ERROR (V1): records {k + 1} and {k + 2} are {step[k]:.0f} s apart, not the declared "
                         f"clock.timestep = {timestep:.0f} s. MEDS reads a uniform axis: insert the missing "
                         f"stamps as missing values, which the gap filling then treats openly.")
    return dict(records=int(len(stamps)), first=str(stamps[0]), last=str(stamps[-1]))


def screen_bounds(values):
    """V4: values outside physical bounds become missing, and a variable with many of them stops the
    build, since that is a unit error. Negative shortwave (a night offset) and rain are set to 0 and
    negative wind speeds are missing, before the bounds apply."""
    report = {}
    values = values.copy()
    if "SWdown" in values:
        values["SWdown"] = values["SWdown"].where(~(values["SWdown"] < 0.0), 0.0)
    if "Rainf" in values:
        values["Rainf"] = values["Rainf"].where(~(values["Rainf"] < 0.0), 0.0)
    for name, (lo, hi) in BOUNDS.items():
        if name not in values:
            continue
        x = values[name]
        present = x.notna()
        bad = present & ((x < lo) | (x > hi))
        n_present = int(present.sum())
        report[name] = dict(present=n_present, out_of_bounds=int(bad.sum()))
        if n_present and bad.sum() / n_present > V4_MAX_BAD_FRACTION:
            raise SystemExit(f"ERROR (V4): {100.0 * bad.sum() / n_present:.1f}% of {name} lies outside "
                             f"[{lo}, {hi}] in MEDS units. Check its declared units.")
        values[name] = x.where(~bad)
    return values, report


def check_sun(stamps_utc, sw, site, shifts_min=range(-60, 61, 5)):
    """V2. Fits the 95th-percentile shortwave by (month, time of day) to a transmissivity times the
    model's window-mean cos z of each record's interval, shifted by each candidate offset; the best
    shift must be within V2_MAX_SHIFT_MIN of zero. Also the fraction of the shortwave that falls in
    intervals whose mean cos z the model treats as night (it would be zeroed)."""
    stamps_utc = np.asarray(stamps_utc, dtype="datetime64[s]")
    sw = np.asarray(sw, dtype=float)
    ok = np.isfinite(sw)
    if ok.sum() < 30 * 86400 / site.timestep:
        raise SystemExit("ERROR (V2): fewer than 30 days of shortwave; the clock cannot be checked")
    start, _ = ti.interval_bounds(stamps_utc, site.stamp, site.timestep)
    frame = pd.DataFrame({"sw": sw[ok]})
    t = pd.DatetimeIndex(stamps_utc[ok])
    frame["month"] = t.month
    frame["slot"] = ((t.hour * 3600 + t.minute * 60 + t.second) // int(site.timestep)).astype(int)
    envelope = frame.groupby(["month", "slot"])["sw"].quantile(0.95)

    def rmse(shift_min):
        cz = ti.conv.window_mean_cosz(start[ok] + np.timedelta64(int(shift_min * 60), "s"), site.timestep,
                                     site.latitude, site.longitude)
        model = pd.Series(cz, index=frame.index).groupby([frame["month"], frame["slot"]]).mean()
        k = float((envelope * model).sum() / max((model ** 2).sum(), 1e-30))
        return float(np.sqrt(((envelope - k * model) ** 2).mean()))

    scores = {s: rmse(s) for s in shifts_min}
    best = min(scores, key=scores.get)
    cz0 = ti.conv.window_mean_cosz(start, site.timestep, site.latitude, site.longitude)
    total = float(np.nansum(np.maximum(sw, 0.0)))
    night = float(np.nansum(np.where(cz0 <= ti.conv.COSZ_BAR_MIN, np.maximum(sw, 0.0), 0.0)))
    night_fraction = night / total if total > 0 else 0.0
    result = dict(best_shift_min=float(best), rmse_at_declared=scores[0], rmse_at_best=scores[best],
                  night_shortwave_fraction=night_fraction)
    if abs(best) > V2_MAX_SHIFT_MIN:
        raise SystemExit(f"ERROR (V2): the shortwave fits the model's sun best {best:+.0f} min from the declared "
                         f"clock (clock.utc_offset = {site.utc_offset:+g} h, clock.stamp = {site.stamp!r}); "
                         f"RMSE {scores[0]:.1f} W/m2 as declared, {scores[best]:.1f} at the best shift. A stamp "
                         f"at the other end of the interval is 30 min for half-hourly data; a clock in the "
                         f"wrong zone is whole hours.")
    if night_fraction > V2_MAX_NIGHT_SW:
        raise SystemExit(f"ERROR (V2): {100.0 * night_fraction:.2f}% of the shortwave falls in intervals the model "
                         f"treats as night under the declared clock; check clock.utc_offset and clock.stamp.")
    return result


def check_vpd(tair_k, rh, vpd_pa, curve):
    """V3: with both RH and VPD, the provider's VPD is (1 - RH) e_s(T) under the declared curve.
    Also names the curve of SATURATION_CURVES that fits best, so a wrong declaration is fixable."""
    tc = np.asarray(tair_k, dtype=float) - 273.15
    rh = np.asarray(rh, dtype=float)
    vpd_pa = np.asarray(vpd_pa, dtype=float)
    ok = np.isfinite(tc) & np.isfinite(rh) & np.isfinite(vpd_pa)
    if ok.sum() == 0:
        return dict(checked=False)
    residual = {name: float(np.percentile(np.abs(vpd_pa[ok] - (1.0 - rh[ok]) * f(tc[ok])), 99))
                for name, f in ti.SATURATION_CURVES.items()}
    best = min(residual, key=residual.get)
    result = dict(checked=True, declared=curve, residual_p99_pa=residual[curve], best_fit=best,
                  best_fit_residual_p99_pa=residual[best])
    if residual[curve] > V3_MAX_RESIDUAL_PA:
        if residual[best] > V3_MAX_RESIDUAL_PA:
            raise SystemExit(f"ERROR (V3): the provider's VPD differs from (1 - RH) e_s(T) by at least "
                             f"{residual[best]:.2f} Pa (99th percentile) under every known saturation curve, so "
                             f"the VPD and RH columns do not describe the same air (another sensor, or rounded "
                             f"values?). Declare RH alone, or VPD alone with its curve.")
        raise SystemExit(f"ERROR (V3): the provider's VPD differs from (1 - RH) e_s(T) under the declared curve "
                         f"{curve!r} by {residual[curve]:.2f} Pa (99th percentile). The best fit is {best!r} "
                         f"({residual[best]:.3f} Pa); set variables.VPD.curve to it.")
    return result


def rh_from_vpd(tair_k, vpd_pa, curve):
    """RH = 1 - VPD / e_s(T) under the provider's curve: the inverse of how the VPD was made."""
    tc = np.asarray(tair_k, dtype=float) - 273.15
    return 1.0 - np.asarray(vpd_pa, dtype=float) / ti.SATURATION_CURVES[curve](tc)


def yearly_report(stamps_utc, final, qc, source_values):
    """V5: per UTC year, each variable's share of each qc code, the rain total, and the PAR/SW ratio
    (a drifting pyranometer or quantum sensor shows as a trend in it)."""
    t = pd.DatetimeIndex(np.asarray(stamps_utc, dtype="datetime64[s]"))
    out = {}
    for year in sorted(set(t.year)):
        m = np.asarray(t.year == year)
        y = {"records": int(m.sum())}
        for name, flags in qc.items():
            codes, counts = np.unique(flags[m], return_counts=True)
            y[f"{name}_qc_fraction"] = {int(c): round(float(n) / m.sum(), 4) for c, n in zip(codes, counts)}
        if "Rainf" in final:
            dt = float(np.median(np.diff(np.asarray(stamps_utc, dtype="datetime64[s]")).astype(float)))
            y["rain_mm"] = round(float(np.sum(final["Rainf"][m]) * dt), 1)
        if "PAR" in source_values and "SWdown" in source_values:
            par, sw = source_values["PAR"][m], source_values["SWdown"][m]
            day = np.isfinite(par) & np.isfinite(sw) & (sw > 100.0)
            if day.sum() > 100:
                y["par_over_sw_umol_per_j"] = round(float(np.median(par[day] / sw[day])), 3)
        out[str(year)] = y
    return out
