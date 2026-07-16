"""Empirical vital-rate laws (the demography-only "example"), in numpy.

These are the phenomenological growth / mortality / recruitment relationships that
were REMOVED from the Fortran core in the reorg (the carbon path is the Fortran
model). They live here as an editable Python example: given the per-cohort stand
state read back through the C-API, they return the three rate arrays the engine's
law-free apply-primitives consume (via ``Site.apply_rates``).

A faithful port of the deleted ``meds_demography_rates`` + the pan-tropical
allometry; the PFT parameters mirror ``meds_config_pft.toml`` / the test config.
"""
import numpy as np

TINY = 1.0e-20
GROWTH_AVG_UNSET = -1.0

# --- pan-tropical allometry (from set_allometry in the shipped/test config) --------------
B1HT, B2HT = 1.139963, 0.564899
AGB_C1, AGB_C2 = 0.06080334, 1.0044785


def dbh_to_height(dbh, hgt_max):
    return np.minimum(np.exp(B1HT + B2HT * np.log(np.maximum(dbh, TINY))), hgt_max)


def height_to_dbh(h):
    return np.exp((np.log(np.maximum(h, TINY)) - B1HT) / B2HT)


def dbh_to_agb(dbh, h, rho):
    return AGB_C1 * rho ** AGB_C2 * (dbh * dbh * h) ** AGB_C2


class PFTParams:
    """3-PFT table matching meds_config_pft.toml / build_test_config."""
    wood_density = np.array([0.40, 0.60, 0.85])
    hgt_max = np.array([46.0, 46.0, 46.0])
    dbh_critical = np.array([100.0, 100.0, 100.0])
    growth_dbh_slope = np.array([0.25, 0.25, 0.25])
    growth_dbh_cap = np.array([100.0, 100.0, 100.0])
    growth_dbh_max = np.array([1.5, 1.0, 0.5])
    growth_lai_slope = np.array([-0.8, -0.7, -0.6])
    repro_invest = np.array([0.3, 0.3, 0.3])
    repro_carbon_eff = np.array([1.0e-3, 1.0e-3, 1.0e-3])
    seed_rain = np.array([0.01, 0.01, 0.01])
    include_pft = np.array([1, 1, 1])
    min_cohort_height = 2.0
    min_reproduction_height = 20.0
    # Camac mortality coefficients: a power law in wood density (derive_pft_rates):
    #   param = param_0 * (rho / rho_ref) ** exp
    _rho_ref = 0.6
    mort_gamma = 0.0094 * (wood_density / _rho_ref) ** (-1.8392)
    mort_alpha = 0.05 * (wood_density / _rho_ref) ** (-1.1493)
    mort_beta = 18.72 * (wood_density / _rho_ref) ** (0.2792)

    @classmethod
    def carbon_min(cls):
        d = height_to_dbh(cls.min_cohort_height)
        return dbh_to_agb(d, cls.min_cohort_height, cls.wood_density)


def growth_intrinsic(dbh, pf, P=PFTParams):
    return P.growth_dbh_max[pf] * np.exp(
        P.growth_dbh_slope[pf] * np.maximum(0.0, np.log(P.growth_dbh_cap[pf]) - np.log(dbh)))


def competition_suppression(over_lai, pf, P=PFTParams):
    return np.exp(P.growth_lai_slope[pf] * over_lai)


def empirical_rates(state, n_patch, P=PFTParams):
    """Return (growth[n], mortality[n], recruitment[n_pft, n_patch]) for the stand.

    ``state`` holds the per-cohort arrays read from the C-API: dbh, height,
    overtopping_lai, growth_avg, agb, nplant, pft (1-based), owner_patch (1-based).
    """
    n_pft = len(P.wood_density)
    recruitment = np.zeros((n_pft, max(n_patch, 1)))
    # Baseline external seed rain (independent of structure), gated by include_pft.
    for k in range(n_pft):
        if P.include_pft[k] == 1:
            recruitment[k, :n_patch] = P.seed_rain[k]

    n = len(state["dbh"])
    growth = np.zeros(n)
    mortality = np.zeros(n)
    if n == 0:
        return growth, mortality, recruitment

    dbh = state["dbh"]
    height = state["height"]
    over = state["overtopping_lai"]
    gavg = state["growth_avg"]
    agb = state["agb"]
    npl = state["nplant"]
    pf = state["pft"].astype(int) - 1          # -> 0-based
    op = state["owner_patch"].astype(int)      # 1-based patch

    gi = growth_intrinsic(dbh, pf, P)
    cs = competition_suppression(over, pf, P)
    mature = height >= P.min_reproduction_height
    repro_supp = np.where(mature, 1.0 - P.repro_invest[pf], 1.0)

    # growth = intrinsic x competition x reproductive-allocation suppression
    growth = gi * cs * repro_supp

    # mortality (Camac): tracked mean growth once seeded, else the instantaneous rate
    g_eff = np.where(gavg == GROWTH_AVG_UNSET, growth, gavg)
    mortality = P.mort_gamma[pf] + P.mort_alpha[pf] * np.exp(-P.mort_beta[pf] * g_eff)

    # reproduction -> recruits: the diverted growth becomes an AGB gain, converted to
    # min-size recruits over the recruit carbon; zero below the maturity height.
    cmin = P.carbon_min()
    repro_dbh = np.where(mature, gi * cs * P.repro_invest[pf], 0.0)
    new_dbh = dbh + repro_dbh
    dagb = dbh_to_agb(new_dbh, dbh_to_height(new_dbh, P.hgt_max[pf]), P.wood_density[pf]) - agb
    rec = np.where(mature, npl * dagb * P.repro_carbon_eff[pf] / cmin[pf], 0.0)
    for j in range(n):
        k = pf[j]
        if P.include_pft[k] == 1 and mature[j]:
            recruitment[k, op[j] - 1] += rec[j]

    return growth, mortality, recruitment
