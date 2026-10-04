# SPDX-License-Identifier: Apache-2.0
"""Unit tests for scripts/calibrate_fast (MEDS_FAST_CALIBRATION_PLAN.md §7 P3 and the revision plan
§9): the registry, its priors and the key selection, the site settings against site_reference.toml,
the transforms, the trial writer, the closure correction, the targets' filters and weights, the
calibrated-config writer, the kernel models' anchor, the water stage's search, the filter
sensitivity, and the fit on a synthetic linear model with a known answer and covariance. No MEDS run."""
import datetime as dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import calibrate_fast as CF   # noqa: E402  (puts the source tree's meds package on the path)
import datarules as DR        # noqa: E402
import fit as F               # noqa: E402
import residuals as R         # noqa: E402
import settings as SET        # noqa: E402
import stages as STG          # noqa: E402
import tower as TW            # noqa: E402
import trials as T            # noqa: E402
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
        assert p.lo < p.hi and p.file in ("pft", "main")


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


def test_closure_restores_the_balance_and_keeps_the_bowen_ratio():
    spec = TW.TowerSpec(site="", closure_days=31, daytime_sw=10.0)
    out = TW.add_closure(synthetic_tower(), spec)
    day = out["sw_in"] > 10.0
    mid = out.index[(out.index > "2016-01-16") & (out.index < "2016-01-25")]
    d = day & out.index.isin(mid)
    assert np.allclose(out.loc[d, "closure_f"], 1.25)
    assert np.allclose((out.loc[d, "h_c"] + out.loc[d, "le_c"]), out.loc[d, "rnet"])
    assert np.allclose(out.loc[d, "h_c"] / out.loc[d, "le_c"], out.loc[d, "h"] / out.loc[d, "le"])
    assert np.all(out.loc[~day, "closure_f"] == 1.0)


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
    spec = TW.TowerSpec.from_site(tmp_path / "site.toml", closure="none")
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
           "rnet": {"sigma_abs": 10.0, "sigma_rel": 0.05}, "le": {"sigma_abs": 10.0, "sigma_rel": 0.15},
           "h": {"sigma_abs": 10.0, "sigma_rel": 0.15}, "ef": {"sigma": 0.05},
           "gpp": {"sigma_abs": 1.5, "sigma_rel": 0.15}, "nee_night": {"sigma": 2.0, "ustar_min": 0.2},
           "ustar": {"sigma_abs": 0.1, "sigma_rel": 0.2}}


def model_like(obs, idx, gr=0.0):
    o = obs.reindex(idx)
    return pd.DataFrame({"sw_in_fast": o["sw_in"], "sw_up_fast": o["sw_up"], "lw_up_fast": o["lw_up"],
                         "rnet_fast": o["rnet"], "le_flux_fast": o["le_c"], "h_flux_fast": o["h_c"],
                         "gpp_rate_fast": o["gpp"], "nee_fast": o["nee"] - gr, "ustar_fast": o["ustar"]},
                        index=idx)


def test_residual_is_zero_on_the_observations_and_rows_are_fixed():
    spec_t = TW.TowerSpec(site="", closure_days=31)
    obs = TW.add_closure(synthetic_tower(closure=1.0), spec_t)
    fok = pd.Series(True, index=obs.index)
    idx = R.window_index(dt.datetime(2016, 1, 10), 10, 3600.0)
    spec = R.build_spec("w", idx, obs, fok, TARGETS, skip_hours=3)
    gr = {m: 1.5 for m in range(1, 13)}
    r = R.residual(spec, model_like(obs, idx, gr=1.5), gr)
    assert np.allclose(r, 0.0)
    names = [t.name for t in spec.targets]
    assert names == [n for n in R.TARGETS if n in names] and "ef" in names and "nee_night" in names
    ef = next(t for t in spec.targets if t.name == "ef")
    assert len(ef.obs) == 10 and np.allclose(ef.obs, 0.7)
    # a model 10 % high in LE moves only the LE and EF rows
    m = model_like(obs, idx, gr=1.5)
    m["le_flux_fast"] *= 1.1
    r2 = R.residual(spec, m, gr)
    sl = spec.slices()
    for name, s in sl.items():
        assert np.any(r2[s] != 0) == (name in ("le", "ef")), name


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
    assert ps["leaf_clumping"].state == "fixed" and ps["leaf_width"].state == "fixed" and ps["ds_jmax"].state == "fit"
    assert ps["rd_vcmax_ratio"].state == "fixed" and ps["stomatal_g1"].prior["centre"] == 3.77
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
    assert d["targets"]["gpp"]["ustar_min"] == "provider" and d["targets"]["gpp"]["sigma_abs"] == 2.5
    assert d["targets"]["nee_night"]["on"] is False
    assert d["tower"]["closure"] == "none" and not d["targets"]["rnet"]["on"] and not d["targets"]["ef"]["on"]
    assert d["targets"]["le"]["ustar_min"] == "diagnostic" and d["targets"]["h"]["ustar_min"] == "diagnostic"
    assert "hours" not in d["targets"]["h"] and d["targets"]["h"]["sigma_rel"] == 0.30
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
    return TW.add_closure(df, TW.TowerSpec(site="", closure_days=3))


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
    obs2 = obs.copy()
    obs2["closure_day"] = np.where(idx.day == 1, 0.5, 0.95)
    spec = R.build_spec("w", idx, obs2, fok, {"le": {"sigma": 10.0, "closure_range": [0.8, 1.2]}})
    assert set(idx.day[spec.targets[0].rows]) == {2}


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
    r = R.residual(spec, df)
    assert np.allclose(r, x * spec.targets[0].weight)
    assert np.allclose(R.unweighted([spec], r), x)
    assert R.sigma_scales([spec], r)["h"] == pytest.approx(np.mean(x ** 2))
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
