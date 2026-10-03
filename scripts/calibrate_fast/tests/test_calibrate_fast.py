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
    assert "dewmx" not in off and "dewmx" in on and "stomatal_g1" in off
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
    spec = TW.TowerSpec(path="", columns={}, closure_days=31, daytime_sw=10.0)
    out = TW.add_closure(synthetic_tower(), spec)
    day = out["sw_in"] > 10.0
    mid = out.index[(out.index > "2016-01-16") & (out.index < "2016-01-25")]
    d = day & out.index.isin(mid)
    assert np.allclose(out.loc[d, "closure_f"], 1.25)
    assert np.allclose((out.loc[d, "h_c"] + out.loc[d, "le_c"]), out.loc[d, "rnet"])
    assert np.allclose(out.loc[d, "h_c"] / out.loc[d, "le_c"], out.loc[d, "h"] / out.loc[d, "le"])
    assert np.all(out.loc[~day, "closure_f"] == 1.0)


def test_hourly_needs_both_half_hours(tmp_path):
    idx = pd.date_range("2016-01-01", periods=8, freq="30min")
    raw = pd.DataFrame({"date": idx, "Rs": 100.0, "Rs_dn": 13.0, "Rl_dn": 450.0, "Rnet": 60.0,
                        "LE": 40.0, "H": 10.0, "NEE": -5.0, "gpp": 10.0, "ustar": 0.3,
                        "FLAG": [1, 1, 1, 0, 1, 1, 1, 1]})
    raw.loc[5, "Rnet"] = np.nan
    raw.to_csv(tmp_path / "t.csv", index=False)
    cols = dict(sw_in="Rs", sw_up="Rs_dn", lw_up="Rl_dn", rnet="Rnet", le="LE", h="H", nee="NEE",
                gpp="gpp", ustar="ustar")
    obs = TW.load_tower(TW.TowerSpec(path=str(tmp_path / "t.csv"), columns=cols))
    assert len(obs) == 4
    assert np.isfinite(obs["le"].iloc[0]) and np.isnan(obs["le"].iloc[1])   # a FLAG = 0 half hour
    assert np.isfinite(obs["rnet"].iloc[1]) and np.isnan(obs["rnet"].iloc[2])  # a missing half hour


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
    spec_t = TW.TowerSpec(path="", columns={}, closure_days=31)
    obs = TW.add_closure(synthetic_tower(closure=1.0), spec_t)
    fok = pd.Series(True, index=obs.index)
    idx = R.window_index(dt.datetime(2016, 1, 10), 10, 0.0)
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
    assert ps["theta_j"].state == "optional" and ps["phi_psii"].state == "fixed"
    assert ps["rd_vcmax_ratio"].state == "fixed" and ps["stomatal_g1"].prior["centre"] == 3.77
    assert {p.stage for p in ps.values()} <= {"optics", "photosynthesis", "energy", "water"}
    assert all(p.reason for p in ps.values() if p.state == "fixed")


# ----- the site settings against site_reference.toml ----------------------------------------------
BCI = Path(__file__).resolve().parents[3] / "examples/example_flux_tower_bci/calibration.toml"


def test_the_bci_declaration_is_complete_and_valid():
    d = SET.complete(load_toml(BCI))
    assert d["targets"]["gpp"]["ustar_min"] == 0.4 and d["targets"]["gpp"]["sigma_abs"] == 2.5
    assert d["targets"]["nee_night"]["on"] is False
    assert d["fit"]["stages"] == ["optics", "photosynthesis", "energy", "water", "polish"]
    assert d["stages"]["polish"]["max_iter"] == 10                      # a default the site left out


def test_settings_refuse_unknown_and_missing_keys():
    good = {"base": {"main": "m.toml", "registry": "r.toml"}, "tower": {"path": "t.csv"},
            "windows": {"chains": {"cal": "2015-01-01"}, "list": []}}
    d = SET.complete(good)
    assert d["targets"]["gpp"]["ustar_min"] == 0.4 and d["windows"]["seasons"] == []
    assert d["stages"]["water"]["windows"] == []                         # the reference's entry is an example
    for path, bad in ((("targets", "gpp"), {"ustar_minn": 0.3}), (("fit",), {"stage": ["polish"]}),
                      (("windows",), {"list": [{"name": "w", "begin": "2015-01-01"}]}),
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
                        ("energy", ["max_iter", "rtol"]), ("water", ["rounds", "grid_sigma", "targets", "windows"])):
        assert set(keys) <= set(ref["stages"][stage]), stage
    for t in R.TARGETS:
        assert t in ref["targets"], t
    for f in R.FILTERS:
        assert f in (Path(SET.REFERENCE)).read_text(), f


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
    return TW.add_closure(df, TW.TowerSpec(path="", columns={}, closure_days=3))


def test_filters_apply_in_order_and_are_counted():
    obs = tower_with()
    idx = obs.index
    fok = pd.Series(True, index=idx)
    cfg = {"gpp": {"sigma_abs": 2.5, "sigma_rel": 0.15, "ustar_min": 0.4, "par_min": 100.0, "hours": [8, 16]}}
    spec = R.build_spec("w", idx, obs, fok, cfg)
    t = spec.targets[0]
    hours = idx.hour.to_numpy()[t.hours]
    assert np.all(hours >= 10) and np.all(hours < 16)                 # u* >= 0.4 only from 10 h; the hours window
    assert np.all(obs["par"].to_numpy()[t.hours] >= 100.0)
    labels = [c[0] for c in t.counts]
    assert labels == ["measured", "forcing observed", "daytime", "u* >= 0.4", "PAR >= 100.0", "local hours [8, 16)"]
    assert all(a[1] >= b[1] for a, b in zip(t.counts, t.counts[1:]))
    assert np.allclose(t.sigma, 2.5 + 0.15 * 15.0)
    with pytest.raises(ValueError):
        R.build_spec("w", idx, obs, fok, {"ustar": {"sigma": 0.1, "ustar_min": 0.3}})
    elev = np.where(idx.hour.to_numpy() == 12, 80.0, 10.0)
    spec = R.build_spec("w", idx, obs, fok, {"albedo": {"sigma": 0.01, "min_sw": 50.0, "min_solar_elevation": 30.0}},
                        elev=elev)
    assert set(idx.hour.to_numpy()[spec.targets[0].hours]) == {12}
    obs2 = obs.copy()
    obs2["closure_day"] = np.where(idx.day == 1, 0.5, 0.95)
    spec = R.build_spec("w", idx, obs2, fok, {"le": {"sigma": 10.0, "closure_range": [0.8, 1.2]}})
    assert set(idx.day[spec.targets[0].hours]) == {2}


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


def test_the_ustar_plateau_finds_the_threshold():
    idx = pd.date_range("2016-01-01", periods=24 * 400, freq="1h")
    rng = np.random.default_rng(1)
    u = rng.uniform(0.05, 1.0, len(idx))
    par = np.full(len(idx), 400.0)
    gpp = par * np.where(u < 0.5, 0.02, 0.03)                       # the deficit stops at u* 0.5
    obs = pd.DataFrame({"sw_in": 200.0, "gpp": gpp, "par": par, "ustar": u}, index=idx)
    rep = R.ustar_plateau(obs, (7, 9))
    assert rep["threshold"]["0.99"] == 0.5


def test_solar_elevation():
    idx = pd.DatetimeIndex(["2016-03-20 11:30", "2016-03-20 23:30"])
    e = TW.solar_elevation(idx, 0.0, 0.0, 0.0)                        # mid-hour: 12:00 and 00:00 UTC
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
