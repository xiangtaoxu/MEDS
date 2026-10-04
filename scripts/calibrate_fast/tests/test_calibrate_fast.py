# SPDX-License-Identifier: Apache-2.0
"""Unit tests for scripts/calibrate_fast (MEDS_FAST_CALIBRATION_PLAN.md §7 P3 and the revision plan
§9): the registry, its priors and the key selection, the site settings against site_reference.toml,
the transforms, the trial writer, the closure correction, the targets' filters and weights, the
calibrated-config writer, the kernel models' anchor, the water stage's search, the filter
sensitivity, and the fit on a synthetic linear model with a known answer and covariance. No MEDS run."""
import datetime as dt
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import calibrate_fast as CF   # noqa: E402  (puts the source tree's meds package on the path)
import datarules as DR        # noqa: E402
import obsmodels as OM        # noqa: E402
import priors as PR           # noqa: E402
import fit as F               # noqa: E402
import residuals as R         # noqa: E402
import settings as SET        # noqa: E402
import stages as STG          # noqa: E402
import tower as TW            # noqa: E402
import trials as T            # noqa: E402
import json                   # noqa: E402
from meds.config import RunConfig, load_toml                     # noqa: E402

try:
    import tomllib
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib
from registry import H_U, SIGMA_U, Param, interval, load_registry, resolve_defaults, select   # noqa: E402

REGISTRY = Path(__file__).resolve().parents[1] / "parameters.toml"


# ----- registry and transforms --------------------------------------------------------------------
@pytest.mark.parametrize("transform,lo,hi,x", [("linear", 0.3, 0.55, 0.45), ("log", 2e-4, 2e-3, 6e-4),
                                               ("linear", -3.5, -1.0, -2.0)])
def test_transform_round_trip(transform, lo, hi, x):
    p = Param("k", "main", "a.k", lo, hi, transform, default=x)
    u = p.to_u(x)
    assert p.to_theta(u) == pytest.approx(x, rel=1e-12)
    for uu in (-30.0, -3.0, 0.0, 3.0, 30.0):
        assert lo <= p.to_theta(uu) <= hi
    eps = 1e-6
    assert p.dtheta_du(u) == pytest.approx((p.to_theta(u + eps) - p.to_theta(u - eps)) / (2 * eps), rel=1e-6)


def test_transform_refuses_a_bound():
    p = Param("k", "main", "a.k", 0.0, 1.0, "linear")
    with pytest.raises(ValueError):
        p.to_u(1.0)


def test_prior_band_is_the_range():
    p = Param("k", "main", "a.k", 0.0, 1.0, "linear", default=0.5)
    assert p.to_theta(p.u0 - 2 * SIGMA_U) == pytest.approx(0.025, abs=1e-9)
    assert p.to_theta(p.u0 + 2 * SIGMA_U) == pytest.approx(0.975, abs=1e-9)
    # the Jacobian step is about 1 % of the range at the centre
    assert p.to_theta(p.u0 + H_U) - 0.5 == pytest.approx(0.01, abs=1e-3)


def test_registry_loads_and_filters_variants():
    off = {p.name for p in load_registry(REGISTRY, "interception_off")}
    on = {p.name for p in load_registry(REGISTRY, "interception_on")}
    assert "leaf_surf_water_max" not in off and "leaf_surf_water_max" in on and "stomatal_g1" in off
    for p in load_registry(REGISTRY, "interception_on"):
        assert p.lo < p.hi and p.file in ("pft", "main", "obs")


def test_defaults_come_from_base_then_record():
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log"),
          Param("z0", "main", "aerodynamics.z0m_ratio", 0.05, 0.2)]
    for p in ps:
        if p.file == "pft":
            p.pft = 1
    record = {("main", "aerodynamics.z0m_ratio", 0): (False, 0.13)}
    resolve_defaults(ps, RunConfig({}, {"pft": {"stomatal_g1": [3.0]}}), record)
    assert ps[0].default == 3.0 and ps[1].default == 0.13
    with pytest.raises(KeyError):
        resolve_defaults([Param("x", "main", "nowhere.key", 0.0, 1.0)], RunConfig({}, {}), record)


# ----- the trial writer ------------------------------------------------------------------------
def test_build_trial_sets_keys(tmp_path):
    main = {"run": {"start_time": "2012-08-01"}, "init": {"init_mode": 1}, "output": {"fast": {}}}
    pft = {"pft": {"stomatal_g1": [3.0, 4.0], "vcmax25": [45.0, 50.0], "wood_density": [0.6, 0.7]}}
    base = RunConfig(main, pft)
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=2),
          Param("z0", "main", "aerodynamics.z0m_ratio", 0.05, 0.2)]
    w = T.Window("w1", dt.datetime(2016, 3, 5), 10, "cal", "cal")
    tdir = T.build_trial(base, ps, [5.0, 0.1], w, "/x/state.nc", tmp_path)
    m = load_toml(tdir / "main.toml")
    q = load_toml(tdir / "pft.toml")
    assert q["pft"]["stomatal_g1"] == [3.0, 5.0]
    assert m["aerodynamics"]["z0m_ratio"] == 0.1
    assert m["run"]["start_time"] == "2016-03-05 00:00:00" and m["run"]["end_time"] == "2016-03-15 00:00:00"
    assert m["init"]["init_mode"] == 2 and m["init"]["restart_file"] == "/x/state.nc"
    assert m["init"]["reacclimate_traits"] is True and m["run"]["slow_on"] is False
    assert m["init"]["pft_config"] == str(tdir / "pft.toml")
    # the same candidate is the same trial; another is another
    assert T.build_trial(base, ps, [5.0, 0.1], w, "/x/state.nc", tmp_path) == tdir
    assert T.build_trial(base, ps, [5.0, 0.11], w, "/x/state.nc", tmp_path) != tdir
    assert base.get("pft.stomatal_g1", file="pft") == [3.0, 4.0]      # the base is not changed


def test_command_runs_the_python_api_or_an_executable():
    assert T.command("python", Path("/t/main.toml")) == [sys.executable, "-m", "meds.model", "/t/main.toml"]
    assert T.command("/b/meds_main", Path("/t/main.toml")) == ["/b/meds_main", "/t/main.toml"]


def test_record_check(tmp_path):
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)]
    (tmp_path / "out").mkdir()
    rec = tmp_path / "out" / f"{T.PREFIX}_parameters.csv"
    main, pftf = tmp_path / "main.toml", tmp_path / "pft.toml"
    rec.write_text("source,key,index,present,value\n"
                   f'"{pftf}",pft.stomatal_g1,1,true,"4.0000000000000000E+000"\n')
    T.check_record(tmp_path, ps, [4.0])
    with pytest.raises(T.TrialError):
        T.check_record(tmp_path, ps, [4.5])
    rec.write_text("source,key,index,present,value\n")
    with pytest.raises(T.TrialError, match="not read"):
        T.check_record(tmp_path, ps, [4.0])


# ----- the closure correction --------------------------------------------------------------------
def synthetic_tower(days=40, closure=0.8):
    idx = pd.date_range("2016-01-01", periods=24 * days, freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    rnet = 0.7 * sw - 30.0
    h, le = 0.3 * closure * rnet, 0.7 * closure * rnet
    df = pd.DataFrame({"sw_in": sw, "sw_up": 0.13 * sw, "lw_up": 450.0, "rnet": rnet, "le": le, "h": h,
                       "nee": np.where(sw > 0, -10.0, 5.0), "gpp": np.where(sw > 0, 15.0, 0.0),
                       "ustar": 0.4}, index=idx)
    return df


def test_the_closure_factor_is_a_median_of_whole_days():
    df = synthetic_tower(closure=0.8)                               # H + LE = 0.8 (Rnet), every day
    values = pd.DataFrame({"Rnet": df["rnet"], "H": df["h"], "LE": df["le"]})
    measured = pd.DataFrame(True, index=df.index, columns=["H", "LE"])
    measured.iloc[24 * 10:24 * 11] = False                          # a day without measured turbulence
    cfg = {"window_days": 15, "min_measured": 0.7, "min_days": 5}
    f, rep = OM.closure_factor(values, measured, 0.0, cfg)
    assert np.allclose(f.dropna(), 1.25) and rep["valid_days"] == 39 and rep["daily_closure_median"] == pytest.approx(0.8)
    obs = df.assign(closure_f=f)
    h_c, le_c = OM.corrected(obs, 1.0)                              # the gap is H's
    assert np.allclose(le_c.dropna(), obs["le"][le_c.notna()]) and np.allclose((h_c + le_c).dropna(), 1.25 * (obs["h"] + obs["le"])[h_c.notna()])
    h_b, le_b = OM.corrected(obs, None)                             # Bowen
    d = h_b.notna() & (obs["le"] > 0)
    assert np.allclose(h_b[d] / le_b[d], obs["h"][d] / obs["le"][d])


@pytest.mark.parametrize("h,le,shares", [(0.5, -0.04, (1.0, 0.0)), (0.5, 0.3, (None, None)), (0.02, 0.3, (0.0, 1.0)),
                                          (0.02, -0.04, (0.0, 0.0))])
def test_the_attribution_test_gives_the_closure_shares(h, le, shares):
    s_h, s_le, why = OM.closure_shares("attribution", {"rises": {"h": h, "le": le}}, 0.10)
    assert (s_h, s_le) == shares and why
    assert OM.closure_shares("bowen", {}, 0.1)[:2] == (None, None) and OM.closure_shares("none", {}, 0.1)[:2] == (0.0, 0.0)


def test_the_attribution_test_finds_the_flux_that_rises_with_turbulence():
    rng = np.random.default_rng(2)
    n = 4000
    idx = pd.date_range("2016-01-01", periods=n, freq="30min")
    u, vpd, rnet = rng.uniform(0.1, 1.2, n), rng.uniform(100.0, 2000.0, n), rng.uniform(100.0, 700.0, n)
    obs = pd.DataFrame({"sw_in": 500.0, "ustar": u, "vpd": vpd, "rnet": rnet,
                        "h": rnet * (0.15 + 0.2 * u), "le": rnet * 0.45}, index=idx)
    test = OM.attribution(obs, {"vpd_classes": 4}, 50.0, 10.0)
    assert test["rises"]["h"] > 0.5 and abs(test["rises"]["le"]) < 1e-9
    assert OM.closure_shares("attribution", test, 0.10)[:2] == (1.0, 0.0)


def test_sigma_from_the_paired_days_at_a_smoothed_observation():
    rng = np.random.default_rng(3)
    idx = pd.date_range("2016-01-01", periods=48 * 200, freq="30min")
    hour = idx.hour.to_numpy() + idx.minute.to_numpy() / 60.0
    true = np.clip(300.0 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    sigma = 10.0 + 0.2 * true
    obs = pd.DataFrame({"le": true + sigma * rng.standard_normal(len(idx)), "sw_in": 2 * true, "par": np.nan,
                        "tair": 300.0, "vpd": 1000.0, "wind": 2.0}, index=idx)
    a, b, rep = OM.paired_sigma(obs, "le", SCFG, 1800.0)
    assert a == pytest.approx(10.0, abs=2.0) and b == pytest.approx(0.2, abs=0.03) and rep["pairs"] > 5000
    sm = OM.smoothed(obs["le"], 0.0, 7)
    noon = idx.hour == 12
    assert np.allclose(sm[noon], 300.0, rtol=0.2)                   # the same time of day over +-7 days (noise / sqrt(15))
    assert abs(sm[noon].mean() / 300.0 - 1.0) < 0.02
    targets = {"le": {"sigma_abs": 99.0, "sigma_rel": 0.0}, "albedo": {"sigma": 0.01}}
    rep = OM.set_sigmas(targets, obs, SCFG, 0.0, 1800.0)
    assert rep["le"]["source"] == "paired days" and targets["le"]["sigma_at"] == "le_smooth" and "albedo" not in rep
    obs["le_randunc"] = 5.0 + 0.1 * np.abs(obs["le"])
    rep = OM.set_sigmas(targets, obs, SCFG, 0.0, 1800.0)
    assert rep["le"]["source"] == "provider" and targets["le"]["sigma_abs"] == pytest.approx(5.0, abs=0.5)


SCFG = {"smooth_days": 7, "paired_dpar": 75.0, "paired_dsw": 35.0, "paired_dt": 3.0, "paired_dvpd": 200.0,
        "paired_dwind": 1.0, "min_pairs": 200, "bins": 10}


def test_kappa_is_an_observation_key_of_the_gpp_residual(tmp_path):
    idx = pd.date_range("2016-01-01", periods=4, freq="30min")
    t = R.TargetRows("gpp", np.arange(4), np.full(4, 10.0), np.ones(4), reco=np.full(4, 4.0))
    spec = R.WindowSpec("w", idx, [t], loss=("l2", 2.0))
    df = pd.DataFrame({"gpp_rate_fast": np.full(4, 12.0)}, index=idx)
    assert np.allclose(R.residual(spec, df, {"kappa": 1.0}), 2.0)
    assert np.allclose(R.residual(spec, df, {"kappa": 2.0 / 3.0}), 0.0)    # 10 + (1.5 - 1) 4 = 12
    #----- kappa is never written to a trial, and never fitted beside a shape key
    base = RunConfig({"run": {}, "init": {}, "output": {"fast": {}}}, {"pft": {"stomatal_g1": [3.0]}})
    kappa = Param("kappa", "obs", "kappa", 0.4, 1.0, prior={"centre": 0.65, "sd": 0.1}, kind="observation",
                  scope="observation")
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)
    resolve_defaults([kappa], base)
    assert kappa.default == 0.65
    cfg = T.with_params(base, [g1, kappa], [4.0, 0.7])
    assert cfg.pft["pft"]["stomatal_g1"] == [4.0] and "kappa" not in json.dumps(cfg.main)
    jv = Param("jv", "pft", "pft.jmax_vcmax_ratio", 1.4, 2.2, shape=True)
    with pytest.raises(ValueError, match="shape key"):
        select([kappa, jv], {"keys": ["kappa", "jv"]}, {})


def test_huber_is_undone_for_the_scores_and_sigma_is_scaled_by_the_misfit():
    idx = pd.date_range("2016-01-01", periods=100, freq="1h")
    spec = R.WindowSpec("w", idx, [R.TargetRows("h", np.arange(100), np.zeros(100), np.ones(100)),
                                   R.TargetRows("le", np.arange(100), np.zeros(100), np.ones(100))])
    df = pd.DataFrame({"h_flux_fast": np.r_[np.full(50, 3.0), np.full(50, -3.0)], "le_flux_fast": 0.5}, index=idx)
    r = R.residual(spec, df)
    assert np.allclose(np.abs(r[:100]), np.sqrt(2 * 2 * 3 - 4))      # Huber beyond c = 2
    assert np.allclose(R.raw([spec], r)[:100] ** 2, 9.0)              # undone for chi^2
    chi2 = R.chi2_per_row([spec], r)
    assert chi2["h"] == pytest.approx(9.0) and chi2["le"] == pytest.approx(0.25)
    scales = R.scale_sigma([spec], r, cap=2.5)
    assert scales == {"h": 2.5, "le": 1.0} and np.allclose(spec.targets[0].sigma, 2.5)   # 3 capped at 2.5; never below 1


SITE_TOML = """
[input]
format = "csv"
path = "t.csv"
timestamp = "date"
timestamp_format = "%Y-%m-%d %H:%M"
[site]
latitude = 9.15
longitude = -79.85
elevation = 100.0
[clock]
utc_offset = -5.0
stamp = "begin"
timestep = 1800
[heights]
tq_height = 30.0
wind_height = 30.0
pressure_height = 30.0
[variables]
Tair = { column = "tair", units = "degC" }
RH = { column = "RH", units = "%" }
PSurf = { column = "p_kpa", units = "kPa" }
Rainf = { column = "PPT", units = "mm" }
SWdown = { column = "Rs", units = "W m-2" }
Wind = { column = "ubar", units = "m s-1" }
[fluxes]
SW_out = { column = "Rs_dn", units = "W m-2" }
LW_out = { column = "Rl_up", units = "W m-2" }
Rnet = { column = "Rnet", units = "W m-2" }
LE = { column = "LE", units = "W m-2", measured = { column = "FLAG", equals = 1 } }
H = { column = "H", units = "W m-2", measured = { column = "FLAG", equals = 1 } }
NEE = { column = "NEE", units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
GPP = { column = "gpp", units = "umol m-2 s-1", measured = { column = "FLAG", equals = 1 } }
RECO = { sum = ["GPP", "NEE"] }
USTAR = { column = "ustar", units = "m s-1", measured = { column = "FLAG", equals = 1 } }
"""


def test_the_tower_is_its_site_toml_on_utc_starts(tmp_path):
    local = pd.date_range("2016-01-01", periods=96, freq="30min")
    raw = pd.DataFrame({"date": local.strftime("%Y-%m-%d %H:%M"), "tair": 25.0, "RH": 80.0, "p_kpa": 99.0,
                        "PPT": 0.0, "Rs": 100.0, "ubar": 2.0, "Rs_dn": 13.0, "Rl_up": 450.0, "Rnet": 60.0,
                        "LE": 40.0, "H": 10.0, "NEE": -5.0, "gpp": 10.0, "ustar": 0.3, "FLAG": 1})
    raw.loc[3, "FLAG"] = 0
    raw.loc[5, "Rnet"] = np.nan
    raw.to_csv(tmp_path / "t.csv", index=False)
    (tmp_path / "site.toml").write_text(SITE_TOML)
    spec = TW.TowerSpec.from_site(tmp_path / "site.toml")
    assert spec.step == 1800.0 and spec.utc_offset_h == -5.0 and spec.lat == 9.15
    obs = TW.observations(spec)
    assert len(obs) == 96 and obs.index[0] == pd.Timestamp("2016-01-01 05:00")     # the native interval, UTC starts
    assert np.isnan(obs["le"].iloc[3]) and np.isfinite(obs["le"].iloc[2])          # FLAG = 0: not measured
    assert np.isnan(obs["reco"].iloc[3]) and obs["reco"].iloc[2] == pytest.approx(5.0)   # RECO = GPP + NEE
    assert np.isnan(obs["rnet"].iloc[5]) and np.isfinite(obs["rnet"].iloc[3])      # radiation: wherever present
    assert obs["par"].isna().all()                                                  # not declared
    assert spec.report["fluxes"]["failures"] == []


# ----- the residuals ---------------------------------------------------------------------------
TARGETS = {"albedo": {"sigma": 0.01, "min_sw": 200.0}, "lw_up": {"sigma": 5.0},
           "le": {"sigma_abs": 10.0, "sigma_rel": 0.15}, "h": {"sigma_abs": 10.0, "sigma_rel": 0.15},
           "gpp": {"sigma_abs": 1.5, "sigma_rel": 0.15}, "ustar": {"sigma_abs": 0.1, "sigma_rel": 0.2}}


def with_closure(df, f=1.0, s_h=1.0):
    """A synthetic tower with its closure factor, corrected H and LE, and a respiration."""
    out = df.assign(closure_f=f, reco=4.0)
    out["h_c"], out["le_c"] = OM.corrected(out, s_h)
    return out


def model_like(obs, idx, kappa=1.0):
    o = obs.reindex(idx)
    return pd.DataFrame({"sw_in_fast": o["sw_in"], "sw_up_fast": o["sw_up"], "lw_up_fast": o["lw_up"],
                         "rnet_fast": o["rnet"], "le_flux_fast": o["le_c"], "h_flux_fast": o["h_c"],
                         "gpp_rate_fast": o["gpp"] + (1.0 / kappa - 1.0) * o["reco"], "nee_fast": o["nee"],
                         "ustar_fast": o["ustar"]}, index=idx)


def test_residual_is_zero_on_the_observations_and_rows_are_fixed():
    obs = with_closure(synthetic_tower(closure=0.8), f=1.25)
    fok = pd.Series(True, index=obs.index)
    idx = R.window_index(dt.datetime(2016, 1, 10), 10, 3600.0)
    spec = R.build_spec("w", idx, obs, fok, TARGETS, skip_hours=3)
    r = R.residual(spec, model_like(obs, idx, kappa=0.8), {"kappa": 0.8})
    assert np.allclose(r, 0.0)
    names = [t.name for t in spec.targets]
    assert names == list(R.TARGETS)
    day = obs.reindex(idx)["sw_in"].to_numpy() > 10.0
    assert all(day[t.rows].all() for t in spec.targets if t.name in ("le", "h", "gpp", "ustar"))
    # a model 10 % high in LE moves only the LE rows
    m = model_like(obs, idx, kappa=0.8)
    m["le_flux_fast"] *= 1.1
    r2 = R.residual(spec, m, {"kappa": 0.8})
    for name, s in spec.slices().items():
        assert np.any(r2[s] != 0) == (name == "le"), name


def test_ess_weights_follow_the_autocorrelation():
    rng = np.random.default_rng(0)
    x = np.zeros(4000)
    for i in range(1, len(x)):
        x[i] = 0.8 * x[i - 1] + rng.standard_normal()
    spec = R.WindowSpec("w", pd.date_range("2016-01-01", periods=len(x), freq="1h"),
                        [R.TargetRows("h", np.arange(len(x)), np.zeros(len(x)), np.ones(len(x)))])
    w = R.ess_weights([spec], x)
    assert w[0] == pytest.approx((1 - 0.8) / (1 + 0.8), abs=0.03)


# ----- the calibrated-config writer ------------------------------------------------------------
def test_set_toml_text_keeps_comments():
    text = "[pft]\nstomatal_g1 = [3.0]   # Medlyn\nvcmax25 = [45.0]\n\n[camac]\nx = 1\n"
    t = CF.set_toml_text(text, "pft.stomatal_g1", 4.25, index=0)
    assert "stomatal_g1 = [4.25]   # Medlyn" in t
    t = CF.set_toml_text(t, "pft.leaf_pi0", -1.8, index=0)
    assert tomllib.loads(t)["pft"]["leaf_pi0"] == [-1.8]
    t = CF.set_toml_text(t, "aerodynamics.z0m_ratio", 0.1)
    d = tomllib.loads(t)
    assert d["aerodynamics"]["z0m_ratio"] == 0.1 and d["camac"]["x"] == 1 and d["pft"]["vcmax25"] == [45.0]


# ----- the fit on a synthetic linear model -----------------------------------------------------
class LinearModel:
    """r(u) = A (u - u_true) + e: a linear least-squares problem with a known answer."""

    def __init__(self, params, A, u_true, e):
        self.params, self.A, self.u_true, self.e = params, A, u_true, e
        self.calls = 0

    def residuals(self, thetas, windows=None):
        self.calls += len(thetas)
        out = []
        for th in thetas:
            u = np.array([p.to_u(t) for p, t in zip(self.params, th)])
            out.append(self.A @ (u - self.u_true) + self.e)
        return out


def linear_problem(k=4, n=300, seed=3):
    rng = np.random.default_rng(seed)
    params = [Param(f"p{j}", "main", f"a.p{j}", 0.0, 10.0, default=5.0) for j in range(k)]
    A = rng.standard_normal((n, k)) * np.array([3.0, 1.0, 0.3, 0.005])[:k]
    u_true = rng.uniform(-1.0, 1.0, k)
    e = 0.3 * rng.standard_normal(n)
    model = LinearModel(params, A, u_true, e)
    prob = F.Problem(model, list(range(k)), np.array([p.default for p in params]))
    return prob, A, u_true, e


def test_lm_finds_the_linear_map_and_its_covariance():
    prob, A, u_true, e = linear_problem()
    k = A.shape[1]
    u0 = prob.u_prior
    # the analytic MAP of ||A(u - u_true) + e||^2 + ||(u - u0)/s||^2
    H = A.T @ A + np.eye(k) / SIGMA_U ** 2
    u_map = np.linalg.solve(H, A.T @ (A @ u_true - e) + u0 / SIGMA_U ** 2)
    out = F.lm(prob, u0, max_iter=30, rtol=1e-12, log=lambda *_: None)
    assert np.allclose(out["u"], u_map, atol=1e-5)
    assert np.allclose(out["J"], A, atol=1e-6)
    spec = R.WindowSpec("w", pd.date_range("2016-01-01", periods=A.shape[0], freq="1h"),
                        [R.TargetRows("h", np.arange(A.shape[0]), np.zeros(A.shape[0]), np.ones(A.shape[0]))])
    cov_u, _ = F.posterior(prob, out["J"], out["r"], [spec], weighted=False)
    assert np.allclose(cov_u, np.linalg.inv(H), rtol=1e-6)
    # the weakly constrained key keeps most of its prior sigma; the strong one does not
    ratio = np.sqrt(np.diag(cov_u)) / SIGMA_U
    assert ratio[3] > 0.8 and ratio[0] < 0.1


def test_screening_flags_dead_and_uninformed_keys():
    prob, A, u_true, e = linear_problem()
    A = A.copy()
    A[:, 2] = 0.0                                   # a dead key
    prob.model.A = A
    u0 = prob.u_prior
    r0 = prob.data([u0])[0]
    J, smooth, failed = F.jacobian(prob, u0, r0, log=lambda *_: None)
    spec = R.WindowSpec("w", pd.date_range("2016-01-01", periods=A.shape[0], freq="1h"),
                        [R.TargetRows("h", np.arange(A.shape[0]), np.zeros(A.shape[0]), np.ones(A.shape[0]))])
    keep, rep = F.screening(prob, J, r0, smooth, failed, [spec], max_free=20, mode="drop", weighted=False)
    assert rep["dead"] == ["p2"]
    assert "p0" in rep["fitted"] and "p2" in rep["fixed"]
    assert all(abs(s - 1.0) < 1e-6 for n, s in rep["smoothness"].items() if n != "p2")
    # the report mode (the default) keeps every requested key and says what the drop mode would fix
    keep, rep = F.screening(prob, J, r0, smooth, failed, [spec], max_free=20, mode="report", weighted=False)
    assert keep == list(range(4)) and rep["fixed"] == [] and "p2" in rep["would_drop"]


def test_linearity_is_exact_for_a_linear_model():
    prob, A, u_true, e = linear_problem()
    out = F.lm(prob, prob.u_prior, max_iter=30, rtol=1e-12, log=lambda *_: None)
    spec = R.WindowSpec("w", pd.date_range("2016-01-01", periods=A.shape[0], freq="1h"),
                        [R.TargetRows("h", np.arange(A.shape[0]), np.zeros(A.shape[0]), np.ones(A.shape[0]))])
    cov_u, w = F.posterior(prob, out["J"], out["r"], [spec], weighted=False)
    lin = F.linearity(prob, out["u"], cov_u, out["J"], log=lambda *_: None)
    for d in lin:
        assert d["delta_phi_plus"] == pytest.approx(1.0, rel=1e-4)
        assert d["ratio_minus"] == pytest.approx(1.0, rel=1e-4)
        assert not d["local_only"]


# ----- the registry's priors and the key selection ------------------------------------------------
def test_prior_sd_maps_into_u():
    p = Param("t", "pft", "pft.theta_j", 0.70, 0.90, default=0.9, prior={"centre": 0.8, "sd": 0.05})
    assert p.centre == 0.8 and p.u0 == pytest.approx(0.0)
    assert p.sigma_u == pytest.approx(0.05 / (0.25 * 0.2))           # d theta / d u at the centre is (b - a) / 4
    q = Param("g", "pft", "pft.stomatal_g1", 1.5, 8.0, "log", default=3.0, prior={"centre": 3.77, "log_sd": 0.35})
    eps = 1e-6
    dlog = (np.log(q.to_theta(q.u0 + eps)) - np.log(q.to_theta(q.u0 - eps))) / (2 * eps)
    assert q.sigma_u * dlog == pytest.approx(0.35, rel=1e-6)          # one prior sd in u is log_sd in log theta
    # a base value at a bound is fine; the prior's centre must be inside
    resolve_defaults([p], RunConfig({}, {}))
    with pytest.raises(ValueError):
        resolve_defaults([Param("t", "pft", "pft.theta_j", 0.70, 0.90, default=0.9)], RunConfig({}, {}))
    lo, hi = interval(q, q.u0, 1.0, 1.96)
    assert 1.5 < lo < q.centre < hi < 8.0 and (q.centre - lo) != pytest.approx(hi - q.centre)


def test_select_follows_the_site():
    menu = lambda: [Param("a", "main", "x.a", 0, 1, state="fit", default=0.5),
                    Param("b", "main", "x.b", 0, 1, state="optional", default=0.5),
                    Param("c", "main", "x.c", 0, 1, state="fixed", default=0.5)]
    assert [p.name for p in select(menu(), {}, {})] == ["a"]
    assert [p.name for p in select(menu(), {"add": ["b"]}, {})] == ["a", "b"]
    assert [p.name for p in select(menu(), {"keys": ["c"]}, {})] == ["c"]
    assert select(menu(), {"remove": ["a"]}, {}) == []
    ps = select(menu(), {}, {"a": {"centre": 0.3, "sd": 0.1, "range": [0.1, 0.9]}})
    assert ps[0].centre == 0.3 and (ps[0].lo, ps[0].hi) == (0.1, 0.9)
    for bad in ({"keys": ["nope"]}, {"keys": ["a"], "add": ["b"]}):
        with pytest.raises(ValueError):
            select(menu(), bad, {})
    with pytest.raises(ValueError):
        select(menu(), {}, {"a": {"width": 1}})


def test_the_registry_menu():
    ps = {p.name: p for p in load_registry(REGISTRY, "interception_off")}
    assert ps["theta_j"].state == "fixed" and ps["phi_psii"].state == "fixed" and ps["ea_vcmax"].state == "optional"
    assert ps["leaf_clumping"].state == "fixed" and ps["leaf_width"].state == "fixed" and ps["ds_jmax"].state == "fixed"
    assert ps["kappa"].state == "fit" and ps["kappa"].file == "obs"
    assert {n for n, p in ps.items() if p.shape} == {"theta_j", "jmax_vcmax_ratio", "phi_psii", "ds_vcmax", "ds_jmax",
                                                      "ea_vcmax", "ea_jmax"}
    assert ps["rd_vcmax_ratio"].state == "fixed" and ps["stomatal_g1"].prior["centre"] == "eeo"
    assert ps["stomatal_g1"].meta["tropical_evergreen_broadleaf"]["centre"] == 3.77
    assert {p.kind for p in ps.values()} <= {"trait", "effective", "observation"} and ps["z0m_ratio"].kind == "effective"
    assert {n for n, p in ps.items() if p.fixed_at} == {"jmax_vcmax_ratio", "ds_vcmax", "ds_jmax"}
    assert all(p.prior.get("sd") or p.prior.get("log_sd") for p in ps.values() if p.state == "fit")   # bounds apart
    assert {p.stage for p in ps.values()} <= {"optics", "photosynthesis", "energy", "water"}
    assert all(p.reason for p in ps.values() if p.state == "fixed")


# ----- the site settings against site_reference.toml ----------------------------------------------
BCI = Path(__file__).resolve().parents[3] / "examples/example_flux_tower_bci/calibration.toml"


def test_the_bci_declaration_is_complete_and_valid():
    d = SET.complete(load_toml(BCI))
    site = TW.ti.read_site(str(BCI.parent / d["tower"]["site"]))
    assert site.timestep == 1800 and site.utc_offset == -5.0
    assert site.fluxes["LW_out"]["column"] == "Rl_dn"                  # the file's longwave labels are swapped
    assert site.fluxes["RECO"] == {"sum": ["GPP", "NEE"]} and site.provider["ustar_threshold"] == 0.4
    assert site.leaf_on_months == list(range(1, 13))
    assert all(site.fluxes[q]["measured"] == {"column": "FLAG", "equals": 1} for q in ("LE", "H", "NEE", "GPP", "USTAR"))
    assert d["targets"]["gpp"]["ustar_min"] == "provider" and d["targets"]["gpp"]["sigma_abs"] == 1.5
    assert set(d["targets"]) == set(R.TARGETS) and d["closure"]["shares"] == "attribution"
    assert d["targets"]["le"]["ustar_min"] == "diagnostic" and d["targets"]["h"]["ustar_min"] == "diagnostic"
    assert d["targets"]["h"]["obs_model"] == "closure" and d["targets"]["gpp"]["obs_model"] == "respiration"
    assert d["priors"]["kappa"] == {"centre": 0.65, "sd": 0.10, "range": [0.4, 1.0], "source": d["priors"]["kappa"]["source"]}
    assert d["fit"]["loss"] == "huber"
    assert d["targets"]["albedo"]["min_solar_elevation"] == 20.0
    assert d["windows"]["list"] == [] and d["windows"]["seasonal"]["list"] == []       # chosen by the rule
    assert d["fit"]["stages"] == ["energy", "water", "polish"]
    assert d["stages"]["water"]["targets"] == ["le"]
    assert d["stages"]["polish"]["max_iter"] == 10                      # a default the site left out


def test_every_target_takes_the_filters():
    """A filter the reference does not list for a target is still accepted on it; a misspelt one is not."""
    good = {"base": {"main": "m.toml", "registry": "r.toml"}, "tower": {"site": "site.toml"},
            "windows": {"days": 10}, "targets": {"gpp": {"hours": [8, 17]}}}
    assert SET.complete(good)["targets"]["gpp"]["hours"] == [8, 17]
    good["targets"]["gpp"]["hourz"] = [8, 17]
    with pytest.raises(ValueError, match="hourz"):
        SET.complete(good)


def test_which_keys_each_stage_moves():
    ps = []
    for name, stage in (("rho", "optics"), ("vc", "photosynthesis"), ("g1", "energy"), ("sref", "water")):
        p = Param(name, "main", f"a.{name}", 0.0, 1.0)
        p.stage = stage
        ps.append(p)
    names = lambda idx: [ps[i].name for i in idx]                                      # noqa: E731
    full = list(CF.STAGE_ORDER)
    assert names(CF.free_of(ps, "energy", full)) == ["g1"]
    #----- a kernel stage not run has its keys fitted in the coupled stage
    assert names(CF.free_of(ps, "energy", ["energy", "water", "polish"])) == ["rho", "vc", "g1"]
    assert names(CF.free_of(ps, "energy", ["photosynthesis", "energy"])) == ["rho", "g1"]
    #----- the polish moves every key, the water keys only with seasonal runs
    assert names(CF.free_of(ps, "polish", ["polish"])) == ["rho", "vc", "g1", "sref"]
    assert names(CF.free_of(ps, "polish", ["polish"], seasonal=False)) == ["rho", "vc", "g1"]
    assert names(CF.free_of(ps, "polish", full, skipped=["optics"])) == ["vc", "g1", "sref"]
    #----- a fit without a polish: every run stage's keys together
    assert names(CF.fitted_keys(ps, ["energy"], [], seasonal=False)) == ["rho", "vc", "g1"]
    ps[2].stage = "dropped"
    assert names(CF.free_of(ps, "polish", ["polish"])) == ["rho", "vc", "sref"]


def test_settings_refuse_unknown_and_missing_keys():
    good = {"base": {"main": "m.toml", "registry": "r.toml"}, "tower": {"site": "site.toml"},
            "windows": {"days": 10}}
    d = SET.complete(good)
    assert d["targets"]["gpp"]["ustar_min"] == "provider" and d["windows"]["list"] == []
    assert d["windows"]["seasonal"]["list"] == []                        # the reference's entry is an example
    for path, bad in ((("targets", "gpp"), {"ustar_minn": 0.3}), (("fit",), {"stage": ["polish"]}),
                      (("windows",), {"list": [{"name": "w", "begin": "2015-01-01"}]}),
                      (("windows",), {"chains": {"cal": "2015-01-01"}}),
                      (("priors",), {"vcmax25": {"centre": 45, "width": 3}})):
        decl = {k: (dict(v) if isinstance(v, dict) else v) for k, v in good.items()}
        node = decl
        for k in path[:-1]:
            node = node.setdefault(k, {})
        node.setdefault(path[-1], {}).update(bad)
        with pytest.raises(ValueError):
            SET.complete(decl)
    with pytest.raises(ValueError, match="required"):
        SET.complete({k: v for k, v in good.items() if k != "tower"})


def test_every_setting_the_tool_reads_is_documented():
    """Each [fit], [stages.*], [uncertainty] and [targets.*] key calibrate_fast.py reads is in
    site_reference.toml (so the reference can neither miss a setting nor document a dead one)."""
    import re
    ref = SET.reference()
    src = (Path(CF.__file__)).read_text() + (Path(R.__file__)).read_text()
    for var, table in (("fc", ref["fit"]), ("pc", ref["stages"]["polish"]), ("site.uncertainty", ref["uncertainty"])):
        for key in re.findall(re.escape(var) + r'\["(\w+)"\]', src):
            assert key in table, f"{var}[{key!r}] is not documented"
    for stage, keys in (("optics", ["max_iter"]), ("photosynthesis", ["passes", "tol", "max_iter"]),
                        ("energy", ["max_iter", "rtol"]), ("water", ["rounds", "grid_sigma", "targets"])):
        assert set(keys) <= set(ref["stages"][stage]), stage
    for t in R.TARGETS:
        assert t in ref["targets"], t
    for f in R.FILTERS:
        assert f in (Path(SET.REFERENCE)).read_text(), f
    for key in re.findall(r'(?:ucfg|cfg|wcfg|scfg)\["(\w+)"\]', Path(DR.__file__).read_text()):
        assert key in ref["ustar"] or key in ref["windows"] or key in ref["windows"]["seasonal"], key


# ----- the targets' filters and weights -------------------------------------------------------------
def tower_with(**cols):
    idx = pd.date_range("2016-01-01", periods=48, freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    df = pd.DataFrame({"sw_in": sw, "sw_up": 0.13 * sw, "lw_up": 450.0, "rnet": 0.7 * sw, "le": 0.5 * sw,
                       "h": 0.2 * sw, "nee": -5.0, "gpp": np.where(sw > 0, 15.0, 0.0),
                       "ustar": np.where(hour < 10, 0.2, 0.6), "par": 2.0 * sw}, index=idx)
    for k, v in cols.items():
        df[k] = v
    return with_closure(df)


def test_filters_apply_in_order_and_are_counted():
    obs = tower_with()
    idx = obs.index
    fok = pd.Series(True, index=idx)
    cfg = {"gpp": {"sigma_abs": 2.5, "sigma_rel": 0.15, "ustar_min": 0.4, "par_min": 100.0, "hours": [8, 16]}}
    spec = R.build_spec("w", idx, obs, fok, cfg)
    t = spec.targets[0]
    hours = idx.hour.to_numpy()[t.rows]
    assert np.all(hours >= 10) and np.all(hours < 16)                 # u* >= 0.4 only from 10 h; the hours window
    assert np.all(obs["par"].to_numpy()[t.rows] >= 100.0)
    labels = [c[0] for c in t.counts]
    assert labels == ["measured", "forcing observed", "daytime", "u* >= 0.4", "PAR >= 100.0", "local hours [8, 16)"]
    assert all(a[1] >= b[1] for a, b in zip(t.counts, t.counts[1:]))
    assert np.allclose(t.sigma, 2.5 + 0.15 * 15.0)
    with pytest.raises(ValueError):
        R.build_spec("w", idx, obs, fok, {"ustar": {"sigma": 0.1, "ustar_min": 0.3}})
    elev = np.where(idx.hour.to_numpy() == 12, 80.0, 10.0)
    spec = R.build_spec("w", idx, obs, fok, {"albedo": {"sigma": 0.01, "min_sw": 50.0, "min_solar_elevation": 30.0}},
                        elev=elev)
    assert set(idx.hour.to_numpy()[spec.targets[0].rows]) == {12}


def test_windows_and_hours_follow_the_towers_interval_and_clock():
    idx = R.window_index(dt.datetime(2016, 1, 1, 5), 10, 1800.0)
    assert len(idx) == 480 and idx[1] - idx[0] == pd.Timedelta(minutes=30)
    obs = tower_with()                                               # an hourly tower on UTC starts
    spec = R.build_spec("w", obs.index, obs, pd.Series(True, index=obs.index),
                        {"le": {"sigma": 10.0, "hours": [12, 13]}}, utc_offset_h=-5.0)
    assert set(obs.index.hour.to_numpy()[spec.targets[0].rows]) == {17}   # local noon at UTC-5
    assert [CF.seconds(d) for d in ("900s", "15min", "1h", 900)] == [900.0, 900.0, 3600.0, 900.0]


def test_the_trial_output_must_be_on_the_towers_interval(tmp_path):
    from netCDF4 import Dataset
    when = pd.date_range("2016-01-01", periods=6, freq="1h")
    with Dataset(tmp_path / f"{T.PREFIX}-F-2016-01-01.nc", "w") as ds:
        ds.createDimension("time", len(when))
        for name, x in (("year", when.year), ("month", when.month), ("day", when.day), ("hour", when.hour),
                        ("minute", when.minute)):
            ds.createVariable(name, "i4", ("time",))[:] = np.asarray(x)
        for v in T.TRIAL_VARIABLES:
            ds.createVariable(v, "f8", ("time",))[:] = 1.0
    assert len(T.read_series(tmp_path, 3600.0)) == 6
    with pytest.raises(T.TrialError, match="not the tower's 1800 s"):
        T.read_series(tmp_path, 1800.0)


def test_weights_huber_and_sigma_scales():
    idx = pd.date_range("2016-01-01", periods=200, freq="1h")
    spec = R.WindowSpec("w", idx, [R.TargetRows("h", np.arange(200), np.zeros(200), np.ones(200))])
    rng = np.random.default_rng(0)
    x = np.zeros(200)
    for i in range(1, 200):
        x[i] = 0.7 * x[i - 1] + rng.standard_normal()
    w = R.set_weights([spec], x)
    assert spec.targets[0].weight == pytest.approx(np.sqrt(w["w/h"]))
    df = pd.DataFrame({"h_flux_fast": x}, index=idx)
    spec.loss = ("l2", 2.0)
    r = R.residual(spec, df)
    assert np.allclose(r, x * spec.targets[0].weight)
    assert np.allclose(R.raw([spec], r), x)
    assert R.chi2_per_row([spec], r)["h"] == pytest.approx(np.mean(x ** 2))
    z = np.array([-5.0, -1.0, 0.0, 1.5, 4.0])
    hz = R.huber(z, 2.0)
    assert np.allclose(hz[1:4], z[1:4]) and np.allclose(hz ** 2, np.where(np.abs(z) <= 2, z ** 2, 4 * np.abs(z) - 4))


UCFG = {"crit": 0.95, "edges": [0.0, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.8, 3.0], "driver_classes": 4,
        "min_driver": {"par": 100.0, "rnet": 50.0}, "min_class_records": 15, "min_records": 500,
        "n_boot": 30, "seed": 1}


def ustar_tower(shape, days=300):
    """A daytime tower whose flux / driver ratio depends on u* by `shape`: deficit below 0.4, flat,
    or rising to the top class."""
    idx = pd.date_range("2016-01-01", periods=24 * days, freq="1h")
    rng = np.random.default_rng(1)
    u = rng.uniform(0.05, 1.5, len(idx))
    par = rng.uniform(150.0, 1800.0, len(idx))
    eff = {"plateau": np.where(u < 0.4, 0.02, 0.03), "flat": np.full(len(idx), 0.03),
           "rising": 0.02 + 0.01 * u}[shape]
    return pd.DataFrame({"sw_in": par / 2.0, "gpp": par * eff, "par": par, "ustar": u,
                         "rnet": par / 3.0, "le": par / 3.0 * 0.5 * eff / 0.03, "h": par / 3.0 * eff}, index=idx)


@pytest.mark.parametrize("shape,outcome", [("plateau", "plateau"), ("flat", "flat"), ("rising", "rising")])
def test_the_ustar_diagnostic_tells_a_plateau_from_a_flat_or_rising_ratio(shape, outcome):
    d = DR.ustar_diagnostic(ustar_tower(shape), "gpp", UCFG, 10.0, 0.0)
    assert d["outcome"] == outcome
    if outcome == "plateau":
        assert d["threshold"] == 0.4 and d["bootstrap"]["outcome_share"]["plateau"] > 0.9
    assert DR.ustar_diagnostic(ustar_tower(shape).iloc[:300], "gpp", UCFG, 10.0, 0.0)["outcome"] == "too_few"


def test_each_target_gets_its_u_star_rule():
    obs = ustar_tower("plateau")
    cfg = {"gpp": {"ustar_min": "provider"}, "le": {"ustar_min": "diagnostic"}, "h": {"ustar_min": 0.6},
           "ustar": {"night": True}, "albedo": {}}
    out, rep = DR.resolve_ustar(cfg, obs, {2016: 0.35, 2017: 0.45}, UCFG, 10.0, 0.0)
    assert out["gpp"]["ustar_min"] == {2016: 0.35, 2017: 0.45} and out["h"]["ustar_min"] == 0.6
    assert out["le"]["ustar_min"] == 0.4 and rep["le"]["diagnostic"]["outcome"] == "plateau"
    assert out["ustar"]["night_ustar"] == {2016: 0.35, 2017: 0.45} and "ustar_min" not in out["albedo"]
    idx = pd.DatetimeIndex(["2016-05-01", "2017-05-01", "2018-05-01"])
    assert DR.threshold_at(idx, {2016: 0.35, 2017: 0.45}).tolist() == [0.35, 0.45, 0.4]   # a missing year: the median
    with pytest.raises(ValueError, match="provider"):
        DR.resolve_ustar({"gpp": {"ustar_min": "provider"}}, obs, None, UCFG, 10.0, 0.0)
    with pytest.raises(ValueError, match="must be a number"):
        DR.resolve_ustar({"gpp": {"ustar_min": "auto"}}, obs, 0.4, UCFG, 10.0, 0.0)


def test_turbulent_targets_are_daytime_unless_night_is_asked_for():
    obs = tower_with()
    idx = obs.index
    fok = pd.Series(True, index=idx)
    day = obs["sw_in"].to_numpy() > 10.0
    spec = R.build_spec("w", idx, obs, fok, {"le": {"sigma": 10.0}, "ustar": {"sigma": 0.1}})
    assert all(day[t.rows].all() for t in spec.targets)
    spec = R.build_spec("w", idx, obs, fok, {"le": {"sigma": 10.0, "night": True, "night_ustar": 0.5}})
    rows = spec.targets[0].rows
    assert (~day[rows]).any() and (obs["ustar"].to_numpy()[rows][~day[rows]] >= 0.5).all()


def test_the_albedo_waits_a_week_after_frost():
    obs = tower_with()
    t = pd.Series(280.0, index=obs.index)
    t.iloc[5] = 270.0                                                # a freezing night on the first day
    obs["days_since_frost"] = DR.days_since_frost(t)
    spec = R.build_spec("w", obs.index, obs, pd.Series(True, index=obs.index),
                        {"albedo": {"sigma": 0.01, "min_sw": 50.0, "snow_free_days": 1}})
    assert spec.targets == []                                         # day 2 is 1 day after the frost: still out
    obs["days_since_frost"] = 10.0
    spec = R.build_spec("w", obs.index, obs, pd.Series(True, index=obs.index),
                        {"albedo": {"sigma": 0.01, "min_sw": 50.0, "snow_free_days": 1}})
    assert set(obs.index.day[spec.targets[0].rows]) == {1, 2} and spec.targets[0].counts[-1][0] == "no frost in 1 d"


def test_keys_whose_process_the_data_never_sample_are_fixed():
    obs = tower_with()
    obs["rain"] = 0.0
    obs.iloc[12, obs.columns.get_loc("rain")] = 1e-4                  # rain at noon of day 1
    w = T.Window("w", dt.datetime(2016, 1, 1), 2, "cal", "w")
    specs = {"w": R.build_spec("w", obs.index, obs, pd.Series(True, index=obs.index), {"le": {"sigma": 10.0}})}
    cov = DR.process_coverage(specs, obs, [w], 10.0)
    assert cov["wet_canopy"] == 2 and cov["night"] == 0 and cov["drought"] == 0
    film = Param("film", "pft", "pft.leaf_surf_water_max", 0.05, 0.3, process="wet_canopy")
    sref = Param("sref", "pft", "pft.wstress_sref_stomata", 0.5, 5.0, process="drought")
    g1 = Param("g1", "pft", "pft.stomatal_g1", 1.5, 8.0)
    fixed = DR.fix_by_coverage([film, sref, g1], cov, 2)
    assert set(fixed) == {"sref"} and "drought" in fixed["sref"]


def two_year_tower():
    idx = pd.date_range("2015-01-01", "2016-12-31 23:00", freq="1h")
    hour = idx.hour.to_numpy()
    sw = np.clip(800 * np.sin(np.pi * (hour - 6) / 12), 0, None)
    obs = pd.DataFrame({"sw_in": sw, "h": 1.0, "le": 1.0, "gpp": 1.0}, index=idx)
    obs.loc[(obs.index.year == 2015) & (obs.index.month == 3), "h"] = np.nan    # March 2015 unmeasured
    return obs


def test_windows_one_per_slot_with_validation_in_another_year():
    obs = two_year_tower()
    fok = {"cal": pd.Series(True, index=obs.index), "val": pd.Series(True, index=obs.index)}
    wcfg = {"days": 10, "slots": 8, "chain_lead_days": 0, "min_score": 0.5}
    chosen, rep = DR.select_windows(obs, fok, wcfg, list(range(1, 13)), obs.index[0], 10.0, 0.0)
    cal = [c for c in chosen if c[2] == "cal"]
    val = [c for c in chosen if c[2] == "val"]
    assert len(cal) == 8 and len(val) == 8
    assert all(c[0].split("_")[1] != v[0].split("_")[1] for c, v in zip(cal, val))   # other years
    assert not any(c[1].year == 2015 and c[1].month == 3 for c in chosen)
    #----- leaf-on months: May to September only, and a chain lead that leaves out 2015
    chosen, _ = DR.select_windows(obs, fok, {**wcfg, "chain_lead_days": 400}, [5, 6, 7, 8, 9], obs.index[0], 10.0, 0.0)
    assert chosen and all(5 <= c[1].month <= 9 and c[1].year == 2016 for c in chosen)
    assert not [c for c in chosen if c[2] == "val"]


def test_seasonal_runs_end_at_the_deepest_water_deficit():
    days = pd.date_range("2015-01-01", "2017-12-31", freq="D")
    dry = days.month.isin([1, 2, 3, 4])
    daily = pd.DataFrame({"rain_mm": np.where(dry, 0.0, 8.0), "pet_mm": 4.0}, index=days)
    deficit = DR.water_deficit(daily)
    assert deficit.max() == 0.0 and deficit.loc["2016-04-30"] == pytest.approx(-4.0 * 121)
    scores = pd.Series(1.0, index=days)
    scores.loc["2017"] = 0.0                                          # 2017 unmeasured
    runs, rep = DR.seasonal_runs(deficit, {"days": 120, "max_runs": 2, "min_deficit_mm": 100.0},
                                 list(range(1, 13)), pd.Timestamp("2015-01-01"), scores, 0.5)
    assert [r[0] for r in runs] == ["dry2016"] and runs[0][1] == dt.datetime(2016, 1, 1)
    assert {y["year"]: y["usable"] for y in rep["years"]} == {2015: False, 2016: True, 2017: False}
    runs, rep = DR.seasonal_runs(deficit, {"days": 120, "max_runs": 2, "min_deficit_mm": 1000.0},
                                 list(range(1, 13)), pd.Timestamp("2015-01-01"), scores, 0.5)
    assert runs == [] and "drought keys are fixed" in rep["note"]


def test_priestley_taylor_and_the_range_coverage():
    pet = DR.priestley_taylor_mm(np.array([298.15]), np.array([200.0]), np.array([400.0]), np.array([1.0e5]))
    assert 3.0 < pet[0] < 6.0                                          # a tropical day: ~4-5 mm
    cov = DR.range_coverage({"t": np.arange(100.0)}, {"t": [10.0, 50.0]})
    assert cov["t"]["share"] == pytest.approx(40.0 / (98.01 - 0.99))   # the record's 1-99 % range


def test_the_share_of_area_whose_canopy_air_top_is_above_the_sensor(tmp_path):
    from netCDF4 import Dataset
    with Dataset(tmp_path / "s.nc", "w") as ds:
        ds.createDimension("c", 3)
        ds.createDimension("p", 2)
        ds.createVariable("height", "f8", ("c",))[:] = [30.0, 40.0, 10.0]
        ds.createVariable("owner_patch", "i4", ("c",))[:] = [1, 1, 2]
        ds.createVariable("patch_area", "f8", ("p",))[:] = [0.25, 0.75]
    out = DR.area_above_sensor(tmp_path / "s.nc", 5.0, 5.0, 41.0)
    assert out["share_above"] == 0.25 and out["top_range"] == [15.0, 45.0]


def test_solar_elevation():
    idx = pd.DatetimeIndex(["2016-03-20 11:45", "2016-03-20 23:45"])
    e = TW.solar_elevation(idx, 0.0, 0.0, 1800.0)                     # mid-record: 12:00 and 00:00 UTC
    assert e[0] > 85.0 and e[1] < -85.0


# ----- the kernel models' anchor, the water stage's search, the filter sensitivity -------------------
class ToyKernel(STG.KernelBase):
    column = "gpp_rate_fast"

    def model_value(self, name):
        return self.series[name]["gpp_rate_fast"].to_numpy()

    def kernel(self, cfg, name):
        a = float(cfg)
        return 0.9 * a * self.drivers[name]["light"]                   # the kernel is 10 % low everywhere


def test_the_anchor_makes_the_kernel_exact_there():
    idx = pd.date_range("2016-01-01", periods=24, freq="1h")
    light = np.linspace(0.0, 1.0, 24)
    ps = [Param("a", "main", "x.a", 0.1, 10.0, default=2.0)]
    obs = 2.0 * light
    spec = R.WindowSpec("w", idx, [R.TargetRows("gpp", np.arange(1, 24), obs[1:], np.ones(23))])
    w = T.Window("w", dt.datetime(2016, 1, 1), 1, "cal", "cal")
    km = ToyKernel(ps, [w], {"w": spec}, {"w": {"light": light}},
                   {"w": pd.DataFrame({"gpp_rate_fast": 2.0 * light}, index=idx)}, make_config=lambda th: th[0])
    km.set_anchor(np.array([2.0]))
    assert np.allclose(km.residuals([np.array([2.0])])[0], 0.0)    # exact at the anchor
    prob = F.Problem(km, [0], np.array([2.0]))
    out = F.lm(prob, np.array([ps[0].to_u(1.0)]), max_iter=30, rtol=1e-14, log=lambda *_: None)
    assert prob.theta(out["u"])[0] == pytest.approx(2.0, rel=1e-3)  # the prior is weak; the data win


def test_grid_search_finds_a_quadratic_minimum():
    prob, A, u_true, e = linear_problem(k=2, n=200)
    best_u, best_c, hist = STG.grid_search(prob, prob.u_prior, rounds=3, width=1.0, log=lambda *_: None)
    H = A.T @ A + np.eye(2) / SIGMA_U ** 2
    u_map = np.linalg.solve(H, A.T @ (A @ u_true - e) + prob.u_prior / SIGMA_U ** 2)
    assert np.allclose(best_u, u_map, atol=1e-3)                     # a quadratic's minimum, found exactly


def test_the_filter_shift_is_the_linear_maps_move():
    prob, A, u_true, e = linear_problem(k=3, n=300)
    out = F.lm(prob, prob.u_prior, max_iter=30, rtol=1e-14, log=lambda *_: None)
    keep = np.arange(0, 300, 2)                                         # another filter: every other row
    A2, e2 = A[keep], e[keep]
    H2 = A2.T @ A2 + np.eye(3) / SIGMA_U ** 2
    u_map2 = np.linalg.solve(H2, A2.T @ (A2 @ u_true - e2) + prob.u_prior / SIGMA_U ** 2)
    r1 = A @ (out["u"] - u_true) + e
    r2 = A2 @ (out["u"] - u_true) + e2
    du = F.shift(prob, out["u"], A, r1, A2, r2)
    assert np.allclose(out["u"] + du, u_map2, atol=1e-8)
    #----- from a point short of the optimum, the shift is still the difference of the two optima
    u_short = out["u"] + np.array([0.3, -0.2, 0.1])
    du = F.shift(prob, u_short, A, A @ (u_short - u_true) + e, A2, A2 @ (u_short - u_true) + e2)
    assert np.allclose(du, u_map2 - out["u"], atol=1e-8)


def test_a_driver_trial_is_its_own_directory_and_writes_the_drivers(tmp_path):
    base = RunConfig({"run": {}, "init": {}, "output": {"fast": {}}}, {"pft": {"stomatal_g1": [3.0]}})
    ps = [Param("g1", "pft", "pft.stomatal_g1", 1.5, 6.0, "log", pft=1)]
    w = T.Window("w1", dt.datetime(2016, 3, 5), 10, "cal", "cal")
    a = T.build_trial(base, ps, [3.0], w, "/x/state.nc", tmp_path)
    b = T.build_trial(base, ps, [3.0], w, "/x/state.nc", tmp_path, kind="driver")
    assert a != b
    text = (b / "output_variables.toml").read_text()
    assert all(v in text for v in T.DRIVER_COHORT + T.DRIVER_SITE) and 'wai_cohort = "D"' in text
    m = load_toml(b / "main.toml")
    assert m["state"]["write_state"] is True and m["output"]["daily"]["enabled"] is True


# ----- keys and priors (best-practice plan §3.5, §4) ------------------------------------------------
def test_kinds_scopes_and_plant_type_priors():
    g1 = Param("g1", "pft", "pft.stomatal_g1", 0.5, 15.0, "log", prior={"centre": "eeo", "log_sd": 0.5},
               meta={"c3_grass": {"centre": 5.25, "log_sd": 0.35}})
    g0 = Param("g0", "pft", "pft.stomatal_g0", 1e-4, 0.1, "log", meta={"c3_grass": {"centre": 0.02, "log_sd": 0.5}})
    out = select([g1, g0], {"plant_type": "c3_grass"}, {})
    assert out[0].prior["centre"] == "eeo" and out[1].prior["centre"] == 0.02     # EEO first, else the type's
    with pytest.raises(ValueError, match="not resolved"):
        out[0].centre
    with pytest.raises(ValueError, match="numerical"):
        select([Param("dt", "main", "fast.x", 0.1, 1.0, kind="numerical")], {}, {})
    with pytest.raises(ValueError, match="kind"):
        select([Param("x", "main", "a.x", 0.1, 1.0, kind="clever")], {}, {})
    with pytest.raises(ValueError, match="EEO centre"):
        select([Param("x", "main", "a.x", 0.1, 1.0, prior={"centre": "eeo"})], {}, {})


def test_kattge_knorr_and_the_least_cost_g1():
    kk = PR.kattge_knorr(25.5)
    assert kk["pft.jmax_vcmax_ratio"] == pytest.approx(1.6975) and kk["leaf_physiology.ds_vcmax"] == pytest.approx(641.105)
    assert kk["leaf_physiology.ds_jmax"] == pytest.approx(640.575) and PR.viscosity_ratio(298.15) == pytest.approx(1.0)
    clim = {"t_day_k": 299.4, "p_day_pa": 99300.0, "vpd_day_pa": 530.0, "par_day": 840.0}
    g1 = PR.eeo_g1(dict(PR.LEAF_DEFAULTS), clim)
    assert 2.6 < g1 < 3.0                                           # BCI's climate: 2.80
    hot = PR.eeo_g1(dict(PR.LEAF_DEFAULTS), {**clim, "t_day_k": 308.0})
    assert hot > g1                                                  # warmer: K + Gamma* up, viscosity down


def test_co2_from_the_run_setting(tmp_path):
    (tmp_path / "co2.txt").write_text("# a CO2 file\n2015  400.0\n2016  402.0\n2017  404.0\n")
    main = tmp_path / "main.toml"
    base = RunConfig({"forcing": {"co2_source": "file", "co2_file": "co2.txt"}}, {"pft": {}}, path=main)
    assert PR.co2_ppm(base, [2015, 2016]) == 401.0
    base = RunConfig({"forcing": {"co2_source": "const", "co2_const": 390.0}}, {"pft": {}}, path=main)
    assert PR.co2_ppm(base, [2015]) == 390.0


@pytest.mark.skipif(not os.environ.get("MEDS_LIB"), reason="the EEO vcmax25 uses MEDS's own leaf (libmeds.so)")
def test_the_coordination_vcmax25_balances_the_two_rates():
    from meds.plant import leaf
    lp = dict(PR.LEAF_DEFAULTS)
    clim = {"t_day_k": 299.4, "p_day_pa": 99300.0, "vpd_day_pa": 530.0, "par_day": 840.0}
    pft = {"theta_j": 0.7, "jmax_vcmax_ratio": 1.7}
    v25 = PR.eeo_vcmax25(lp, pft, clim, 400.0)
    assert 20.0 < v25 < 80.0
    brighter = PR.eeo_vcmax25(lp, pft, {**clim, "par_day": 1200.0}, 400.0)
    assert brighter > v25                                            # more light: more Rubisco to match it


def test_the_prior_z_and_who_pushes_a_key():
    prob, A, u_true, e = linear_problem(k=3, n=200)
    out = F.lm(prob, prob.u_prior, max_iter=30, rtol=1e-14, log=lambda *_: None)
    spec = R.WindowSpec("w", pd.date_range("2016-01-01", periods=200, freq="1h"),
                        [R.TargetRows("le", np.arange(100), np.zeros(100), np.ones(100)),
                         R.TargetRows("h", np.arange(100), np.zeros(100), np.ones(100))])
    z = F.prior_z(prob, out["u"], A, out["r"], [spec])
    for k, p in enumerate(prob.params):
        assert z[p.name]["z"] == pytest.approx((out["u"][k] - p.u0) / p.sigma_u)
        assert z[p.name]["pushed_by"] in ("le", "h")


def test_a_true_false_setting_from_the_parameter_record():
    """The record holds a logical as text: "false" must read as false (it once switched off the
    Kattge & Knorr values, as if the model's own acclimation were on)."""
    site = CF.Site.__new__(CF.Site)
    site.base = RunConfig({"leaf_physiology": {}}, {"pft": {}})
    rec = {("main", "leaf_physiology.thermal_acclimation", 0): (False, "false")}
    assert site.flag("leaf_physiology.thermal_acclimation", rec) is False
    rec = {("main", "leaf_physiology.thermal_acclimation", 0): (True, "true")}
    assert site.flag("leaf_physiology.thermal_acclimation", rec) is True
    site.base = RunConfig({"leaf_physiology": {"thermal_acclimation": True}}, {"pft": {}})
    assert site.flag("leaf_physiology.thermal_acclimation", None) is True
