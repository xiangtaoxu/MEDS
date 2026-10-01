# SPDX-License-Identifier: Apache-2.0
"""Unit tests for scripts/calibrate_fast (MEDS_FAST_CALIBRATION_PLAN.md §7 P3): the registry and
transforms, the trial writer, the closure correction, the residuals, the calibrated-config writer,
and the fit on a synthetic linear model with a known answer and covariance. No MEDS run."""
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
import tower as TW            # noqa: E402
import trials as T            # noqa: E402
from meds.config import RunConfig, load_toml                     # noqa: E402

try:
    import tomllib
except ModuleNotFoundError:  # py3.10 and older
    import tomli as tomllib
from registry import H_U, SIGMA_U, Param, load_registry, resolve_defaults   # noqa: E402

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
    keep, rep = F.screening(prob, J, r0, smooth, failed, [spec], max_free=20)
    assert rep["dead"] == ["p2"]
    assert "p0" in rep["fitted"] and "p2" in rep["fixed"]
    assert all(abs(s - 1.0) < 1e-6 for n, s in rep["smoothness"].items() if n != "p2")


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
