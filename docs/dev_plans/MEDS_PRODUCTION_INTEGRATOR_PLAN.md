# MEDS production numerical scheme — the plan

**Status:** design, rewritten 2026-07-31 (second revision, same day). The first revision replaced the
retired default-choice plan with a production-efficiency plan; this one adds the **stability analysis
of §1**, which changes what the workstreams should be, and reorganises into two parts:

- **Part I — testing and comparing integrators** (§3–§4): the measurement program. Nothing in Part II
  is built before the relevant Part I number exists.
- **Part II — numerical improvement and multi-core** (§5–§8): the build program.

**Question this answers:** production runs were bounded by a `dt_fast` of 150 s set by *stability*,
not accuracy — 6× more fast steps than the model used to take, each paying a full frozen pre-pass.
What raises that bound, what makes each step cheaper, and what parallelises — in what order, and gated
on which measurement?

**Answered, in part, 2026-07-31.** The stability bound is **gone**: one frozen coefficient (the
CAS↔atmosphere conductance) carried 99% of it, and re-solving it per stage — 2% of the pre-pass —
removes it while *reducing* integrator work. `dt_fast` is now bounded by **accuracy**, and the binding
quantity is carbon, not temperature: GPP runs 33% low at 900 s through a second lagged feedback, this
one in the leaf water potential. That is N2b and it is the top of the queue. Multi-core (§7) is
untouched and still open.

**Companions.** `docs/science/numerical_scheme.md` is the user-facing description of the two schemes
and the freeze. `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–14 is the measurement record for the
oscillation. `MEDS_NUMERICS_SCOPING.md` §11 is the bare-array kernel convention. The archived parity
document is history and is not evidence — see §2b.

---

## 0. What is void from the previous version

The original subject of this file — which integrator family becomes the default — is **closed**.
`split` is retired, `ark` (ESDIRK2) is the default, `rk45` is the accuracy baseline, and its stiff
rescue redoes the step on `ark` (a within-family *stepper* swap, not a model swap). PR #88 closed the
last physics rows: `numerical_scheme.md` §4 now reads "there are no *physics* rows left."

| previous item | disposition |
|---|---|
| F3 — the `split` ↔ ARK/RK45 ~0.45 K gap; P1-a "attribute it" | **VOID.** The family was deleted, so there is nothing to attribute it to. |
| F4 — `ark @ 900 s` beats `split @ 225 s` on the efficiency frontier | **VOID** (one arm retired) and **stale** (both cells inside the oscillating regime). |
| P0-a…e — the five correctness blockers | **All landed.** `[fast].snow_on` and the `with_mass`/`with_theta` norm switches deleted. |
| P1-a′ — the three Class-1 physics rows | **CLOSED** (PR #88). |
| P2 — retire the cross-family fallback | **DONE** (`meds_fast_step.f90:97-113`). |
| P3 — a validator whose main rule was "reject `split` below 900 s" | **Rule void**; the cost-ceiling fragment survives into §6. |
| The layer-3 `dt_fast` → expected-error table | **Deleted** — every entry was measured inside the oscillating regime (§2b). |
| Workstream A variants A1/A2 (aerodynamics / surface coupling explicit) | **Refuted** before being built; retained only as §5 N6. |
| Workstream B3 as first drafted ("quasi-steady CAS, biggest stability win") | **Inverted by §1.** Naive `wcap → 0` is the *maximally unstable* configuration. Rewritten as N4. |

The one structural rule that survives intact — *adapt the step size within a run, choose the scheme at
configuration time* — is now honoured by construction, since the RK45 rescue is a stepper swap on one
semi-discretisation. §5 N5 proposes the first principled exception, at layer 3 rather than layer 2.

---

## 1. The governing diagnosis

### 1a. The cost model

```
  cost  =  n_patch  x  (dt_slow / dt_fast)  x  [ C_freeze  +  n_substep x C_stage x 2 ]
                        ^^^^^^^^^^^^^^^^^        ^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^
                        freeze cadence           pre-pass      the ESDIRK march
                        (stability-bounded)      (once per     (adaptive, ~2-3 substeps
                                                  dt_fast)      x 2 implicit stages)
```

Raising `dt_fast` is the only lever that cuts *both* bracketed terms, which is why §5 leads with it.
`C_freeze` is very probably dominant and **has never been measured** — T1.

### 1b. The canopy-air storage term is not what is failing

The canopy air is already implicit. `column_be_stage` commits it in flux form,

```
H_{n+1} = ( wcap*H_n + dt*(src_enth + gah*H_atm) ) / ( wcap + dt*gah )
```

inside an L-stable ESDIRK2, with the leaf↔CAS pair closed by a 2×2 Newton every stage
(`newton_surface_solve`). The fast mode is handled by an unconditionally stable method; small `wcap`
is not breaking the time discretisation.

Write the memory weight the old state carries:

```
a  =  wcap / (wcap + dt*gah)  =  1 / (1 + dt/tau),      tau = wcap/gah
```

**Measured** on the reference fixture at midday (§1g): `wcap = 24 kg/m2`, `gah = 1.68e-2 kg/m2/s`,
so `tau ≈ 1430 s` — an order of magnitude longer than a first estimate of 80–120 s assumed, because
`u_ref = 2 m/s` over an 18 m canopy gives `ustar = 0.167`, not the 0.3–0.6 of a windy day.

| `dt_fast` | `dt/tau` | `a` — fraction of `H_n` surviving |
|---|---|---|
| 900 s | 0.63 | **0.61** |
| 225 s | 0.16 | 0.86 |
| 150 s | 0.10 | 0.91 |
| 100 s | 0.07 | 0.93 |

So at production settings the canopy air is **not** nearly-diagnostic — it retains 60–90% of its state
across a step, and the storage term is doing real damping. The instability therefore cannot be blamed
on the storage treatment at all; it is entirely in the lag (§1g).

### 1c. The instability is an outer map on the frozen coefficients

`gah`, `g_tr_f`, `h_coeff_f` and `t_emit` are frozen at state `n`, so across one `dt_fast` the
canopy air obeys a map `H_{n+1} = Phi(H_n)`. Linearised, it splits into storage memory plus lag
feedback:

```
Phi'  =  a  +  (1 - a) * L        L = the quasi-steady map's coefficient-lag gain
```

`a` is positive and bounded by 1 — it is the **contraction**. `L` is the destabiliser.

**An earlier draft of this section fitted `L ≈ −2.3` to the amplitude scan. That estimate is
REFUTED — `Phi'` has now been measured directly (§1g) and the true lag gain is an order of magnitude
larger.** The structural decomposition survives; the number did not. Amplitude is a property of the
*nonlinear saturation*, not of the multiplier, and fitting one to the other was the error.

### 1d. Consequences that invert two earlier proposals

**Naive `wcap → 0` (a diagnostic canopy air with frozen coefficients) is the worst case, not the fix.**
It sets `a = 0`, deleting the damping the measured `a = 0.61–0.91` is currently supplying, and pins
the multiplier at `L ≈ −40 to −60`. This inverts the first draft's B3.

**Inner sub-stepping is very slightly harmful, which explains a previously unexplained observation.**
`n` sub-steps of `dt/n` give `a_eff = (1 + dt/(n·tau))^-n`, which *decreases* with `n`, driving the
canopy air closer to the frozen-coefficient equilibrium. `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9
recorded that "no amount of inner adaptivity fixes it" as an empirical surprise; it is a prediction of
this model.

**"Diagnostic" has two meanings and they land on opposite sides.** The second one works:

- `wcap → 0` with coefficients still frozen at state `n` — an explicit map with no damping. Above.
- An **index-1 DAE**: solve `S(H, c(H)) + gah(H)·(H_atm − H) = 0` with the coefficients evaluated at
  the *same* `H` as the balance. There is no lag, so `L` does not exist. This is what CLM5 does — its
  canopy air *is* diagnostic, and it works because the Monin–Obukhov and stomatal updates sit inside
  the surface-energy iteration. ED2 gets the same effect by recomputing `canopy_turbulence8` at every
  RK stage.

**Therefore "solve the CAS algebraically" and "put `∂c/∂y` into the Newton" are the same fix.** The
2×2 Newton already treats `H` and `q` implicitly; what it holds fixed are the Jacobian entries coming
through the coefficients. Adding them converts the lagged map into a genuine implicit solve. That is
the unifying statement of §5.

### 1e. Why the adaptive controller cannot find this by itself

A fair objection: MEDS already runs an embedded-error adaptive march — should it not discover that it
needs smaller steps? **No, and the reason is structural.**

The controller integrates `y' = f(y, c_frozen)` over `[0, dt_fast]` and controls the **local
truncation error of that ODE**. The frozen problem is smooth and the method L-stable, so the embedded
pair agrees closely, `err ≈ 0`, and the controller grows the step. The error actually being made is
`c_frozen ≠ c(y(t))` — a **semi-discretisation error, not a truncation error**. Both members of the
embedded pair solve the *same wrong right-hand side*, so their difference is blind to it by
construction. And the oscillation is a property of the map *between* `dt_fast` steps, which no
within-step estimator ever samples.

This is the same blindness as the conservation ledgers, for the same reason: every detector currently
in the code measures a within-step property of the frozen problem. Three independent instruments —
the ledgers, the embedded estimator, and the Newton's own convergence test — all report health while
the canopy air swings 8 K step to step.

**The corollary is the design:** give the freeze its own error estimator. The next step's pre-pass
computes `c(y_{n+1})` anyway, so the coefficient drift `||c(y_{n+1}) − c(y_n)||` is available at zero
extra cost, and it is the natural controlled quantity for an adaptive *freeze cadence* (N5).

### 1f. What this predicted, and how it scored

| # | prediction | outcome |
|---|---|---|
| 1 | `Phi'` crosses −1 between `dt_fast` 150 s and 225 s | **WRONG** — it crosses between **75 s and 100 s** (§1g T2a). The default 150 s is unstable, `Phi' = −2.74`. |
| 2 | `Phi' = a + (1−a)L` with `L` roughly constant | **PARTLY** — the form holds and `Phi'` is monotone in `dt_fast`, but the implied `L` drifts −37 → −62 across the sweep, so the one-mode model is a caricature. |
| 3 | under-relaxation with `θ < 0.5` stabilises 900 s | **WRONG** — with the measured `L ≈ −62`, `θ < 0.07` would be needed at 900 s. N1 is demoted. |
| 4 | a young stand is *stricter* than 150 s | **REFUTED** — the relationship is non-monotone. The worst stand is **mid-height (10 m, `Phi' = −8.2` at 150 s)**; the 2 m regenerating stand is the *most* stable of the five (§1h). |
| 5 | the tissue store makes 150 s pessimistic | **UNTESTED** — moot once N2a lands, since the CAS is then stable at every cadence measured. |

### 1h. MEASURED, T3 and T5

**T3 — the `ark` ↔ `rk45` gap is an artifact of the instability, not a scheme difference.**
5 cohorts, 24 h, scored against `ark` + refresh at `dt_fast = 12.5 s`, `rtol_all = 1e-9`,
`atol_scale = 1e-3`:

| | `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s |
|---|---|---|---|---|---|---|
| **frozen** | `\|ark − rk45\|` CAS-T, default tol | 0.0001 | 0.0002 | 0.0002 | 0.014 | **0.491** |
| **frozen** | ...tightened tol | 0.0001 | 0.0002 | 0.0002 | 0.010 | **0.483** |
| **refreshed** | `\|ark − rk45\|` CAS-T, default tol | 0.0009 | 0.0018 | 0.0027 | 0.0036 | 0.0074 |
| **refreshed** | ...tightened tol | 0.0005 | 0.0009 | 0.0012 | 0.0013 | **0.0013** |

With the refresh and a tight tolerance the two schemes — one implicit, one explicit — agree to
**1.3 mK at the production cadence**. Frozen, the same comparison reads 0.49 K, **370× larger**, and
tightening the tolerance does not shrink it. That 0.49 K is not a scheme difference: it is two
steppers sampling different phases of the same period-2 oscillation. This is the same trap that
produced the retired `split` gap, and it is worth stating as a rule: **a cross-scheme difference
measured in an unstable regime measures the instability, not the schemes.**

**The freeze dominates and the integrator tolerance is irrelevant — confirmed at a *stable* cadence.**
Refreshed `ark` at 150 s: 0.0174 K at default tolerance, 0.0180 K tightened. The time-stepping error
is ~0.6 mK against a 17 mK total, so ~97% of the error is the Category-0 freeze. (Frozen, the same
test is meaningless — tightening moves `ark` from 0.391 to 0.421 K, i.e. *worse*, because there is no
convergence to measure.)

**Convergence.** Refreshed `ark` errors 0.0006 / 0.0015 / 0.0024 / 0.0052 / 0.0174 K over
25 → 150 s: roughly **second order**. Frozen: 0.0011 / 0.0021 / 0.0032 / **0.153** / **0.391** — a
discontinuity between 75 and 100 s, which is the stability boundary, not a convergence curve. The
frozen scheme has no order across it.

**Carbon and water at 150 s** (against the same reference): GPP −3.8% refreshed vs −4.1% frozen; ET
**−2.6% refreshed vs −9.0% frozen**. ET improves 3.5×. GPP barely moves, because it is limited by the
*leaf gas-exchange* freeze, which N2a does not touch — that, and the ground fluxes, are what bound
`dt_fast` once the canopy air is stable.

`rk45_rescue = 0` in every cell, so every RK45 number above is genuinely RK45.

**T5 — the default is unsafe on most stands, and prediction #4 was backwards.**
`Phi'` at the probe state, `can_depth` following the stand (`max(min_canopy_depth, h + freeboard)`):

| stand | `wcap` | `Phi'` @ 50 s | @ 100 s | @ 150 s | @ 900 s | p2p @ 150 s |
|---|---|---|---|---|---|---|
| 2 m regen, LAI 1 | 8.4 | +0.16 | −0.44 | **−0.90** | −2.48 | 0.008 K |
| 5 m, LAI 2 | 12.0 | −0.93 | −2.51 | **−3.84** | −10.3 | 3.36 K |
| **10 m, LAI 3** | 18.0 | **−2.50** | −5.53 | **−8.22** | −25.9 | 2.07 K |
| 18 m, LAI 3 | 27.6 | −1.15 | −3.08 | −4.84 | −19.6 | 1.70 K |
| 30 m, LAI 4 | 42.0 | −1.10 | −3.05 | −4.88 | −23.3 | 0.79 K |

**Four of five stands are unstable at the default `dt_fast = 150 s`**, and the worst is *worse* than
the 18 m fixture the default was calibrated on (−8.2 vs −4.8). The frozen scheme would need
`dt_fast ≈ 50 s` — 3× more steps than today — to be stable on a 10 m stand. So N2a is not only an
optimisation: **without it, the current production default is numerically unsafe on ordinary stands.**

The prediction that a *shorter* stand is stricter is refuted. Capacity does fall with height, but so
does the flux driving it: over this height–LAI sequence the ratio LAI/`wcap` peaks in the middle
(0.17 at 5–10 m against 0.12 at 2 m and 0.10 at 30 m), and `Phi'` tracks that, not `wcap` alone. The
mechanism in §1c survives; the monotone corollary does not. **Caveat: LAI was varied with height, so
height and leaf area are confounded by construction — a 2 m stand carrying LAI 3 was not tested.**

**With the refresh, every stand is stable or marginal at every cadence measured**: `Phi'` at 150 s
runs −0.07 / −0.47 / −0.75 / −0.89 / −1.16 across the five stands, and peak-to-peak stays ≤ 0.04 K at
150 s and ≤ 0.44 K at 900 s. The 30 m stand is marginal (`Phi' = −1.16` at 150 s, −1.02 at 100 s) but
its amplitude is 0.017 K, so it is a weak mode rather than a live oscillation — worth re-checking on a
real tall-canopy configuration.

Also visible, and a caution for anyone reading amplitudes: the 2 m stand has the *mildest* multiplier
at 900 s (−2.48) and the *largest* excursion (20.4 K p2p). Multiplier and amplitude are decoupled by
the saturating nonlinearity. Rank stability by `Phi'`, not by peak-to-peak.

**Caveat on the probe.** The T5 base state was spun up at 50 s with the refresh off, which is itself
inside the unstable regime for the 5–30 m stands, so each probe is taken somewhere on a limit cycle
rather than at a fixed point. `Phi'` is a local derivative and the comparison across stands is
like-for-like, but the absolute values carry that.

### 1g. MEASURED, 2026-07-31 — the multiplier, its cause, and the cost split

Harness: the `test_column_dynamics` high-LAI sunlit column (LAI 3, 18 m canopy, single cohort),
spun up once to 12:00 solar and probed from that **one common base state**; central differences on
`bio%cas%can_enthalpy` applied *before* `build_column_frozen`, so what is measured is the lagged
outer map. `build_column_frozen` and `adaptive_ark_march` are public, so the coefficient
decomposition needed no source change.

**(i) The map is smooth, and there is no stable `dt_fast` near the default.**

| `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s | 225 s | 450 s | 900 s |
|---|---|---|---|---|---|---|---|---|
| `Phi'` | +0.36 | −0.26 | −0.88 | **−1.49** | **−2.74** | **−4.64** | **−10.6** | **−23.2** |

The stability boundary is at **~75–80 s**, not 150 s. `Phi'` converges cleanly as the perturbation
shrinks (−23.40, −23.37, −23.21 at δ = 10, 50, 100 J/kg, forward and backward differences agreeing to
1%), so **the map is differentiable** — a Jacobian-based fix is viable, and the switching-nonlinearity
worry is dismissed. Nonlinearity appears only above δ ≈ 0.5 K, where `|Phi'|` *falls* — a saturating
nonlinearity, which is what bounds the limit-cycle amplitude and is why amplitude is a poor proxy for
the multiplier.

**(ii) One coefficient carries essentially all of it: the CAS↔atmosphere conductance.** Holding each
frozen coefficient at its *unperturbed* value (i.e. "this coefficient did not see the perturbation"):

| pinned | `Phi'` @ 900 s | `Phi'` @ 150 s |
|---|---|---|
| nothing | −23.24 | −2.74 |
| **`gah`/`gaw`/`gac`** | **+0.80** | **+0.88** |
| `g_tr_f` (stomatal) | −23.33 | −2.77 |
| `h_coeff_f` (leaf sensible) | −23.29 | −2.76 |
| `abs_lw` (`t_emit`) | −23.24 | −2.74 |
| `f_wet_c` | −23.24 | −2.74 |
| all of the above | +0.66 | +0.83 |

Removing the aerodynamic conductance's lag alone takes the map from `Phi' = −23.2` to **+0.80 —
stable at 900 s**. Everything else contributes ≤ 0.7%. The mechanism is the Monin–Obukhov stability
feedback: a warmer canopy air is a more unstable surface layer, which raises `ustar` and the transfer
factor, which vents the canopy air harder — a negative feedback of very high gain (measured
`d ln gah/dT ≈ 2.2 K⁻¹`: +0.5 K gives `gah` ×3.1, +2 K gives ×12.2), lagged by exactly one `dt_fast`.

**This reconciles with, rather than contradicts, the earlier "aerodynamics only amplifies" result.**
That experiment pinned `ustar`. But `gah = rho·ustar·temp1`, and `temp1` — the MO scalar-transfer
factor — is *also* stability-dependent and was left free. At +0.5 K the measurement gives `ustar` ×1.75
against `gah` ×3.06, so pinning `ustar` alone removes only about half the log-sensitivity. The prior
test did not pin what it intended to pin.

**(iii) It is a light-to-moderate-wind phenomenon, which is the common case.** `Phi'` (full / with
`gah` pinned):

| `u_ref` | `ustar` | `Phi'` @ 150 s | `Phi'` @ 900 s |
|---|---|---|---|
| 1 m/s | 0.100 (at the floor) | +0.31 / +0.89 | −3.63 / +0.94 |
| **2 m/s** | 0.167 | **−2.74** / +0.88 | **−23.2** / +0.80 |
| 5 m/s | 0.630 | −0.48 / +0.48 | −2.42 / −0.02 |
| 10 m/s | 1.640 | −0.066 / +0.04 | −0.30 / −0.16 |

Strong winds are stable (high ventilation, weak stability sensitivity). Very light winds are stable
too, but only because `ustar` is pinned at its 0.1 floor, which incidentally zeroes the feedback — an
accidental stabiliser, not a designed one. The danger band is **light-to-moderate wind**, which is
the ordinary condition over a forest canopy.

**(iv) T1 — the pre-pass dominates, and the aerodynamics is a rounding error inside it.**
Copy-overhead subtracted, `dt_fast = 150 s`, 10 soil layers, 20 000 repeats:

| cohorts | `C_freeze` | whole sub-step | freeze share |
|---|---|---|---|
| 1 | 14.7 µs | 37.4 µs | 39% |
| 3 | 30.0 µs | 53.4 µs | **56%** |
| 10 | 86.9 µs | 137.2 µs | **63%** |
| 30 | 325.7 µs | 366.7 µs | **89%** |

The freeze scales with cohort count (leaf gas exchange, ~11 µs/cohort); the march does not. At any
realistic patch occupancy **`C_freeze` is the majority of the cost**, so cheapening the march is
low-value and §6 is confirmed as minor. And within the pre-pass, `canopy_aerodynamics` is **2.0%**
(0.28 µs of 13.9 µs) — so refreshing it at four stage points costs **~8% of one pre-pass**, roughly
3% of a sub-step.

**(v) The synthesis.** The single coefficient that causes the instability is also one of the cheapest
things in the pre-pass. Refreshing `gah`/`gaw`/`gac` per stage — and nothing else — is predicted to
move the stability limit from ~75 s to beyond 900 s for ~3% more work per sub-step, i.e. **up to a 6×
reduction in fast-loop cost.** That is now N2 and it is the plan's headline.

**Caveats.** `Phi'` is `∂H_{n+1}/∂H_n`, one diagonal entry of the map's Jacobian, not its spectral
radius; it is a strong indicator, not a proof of stability. One fixture, one time of day, one cohort,
one soil column, constant-forcing MVP (no live met, no RT join). The pinned-`gah` values (+0.80 at
900 s) are close enough to 1 that a different stand could still sit outside. Reproduce with
`scratchpad/meas_fastloop{2,3,4}.f90`.

---

## 2. Evidence

### 2a. What is measured and current

`dt_fast` is a stability boundary. Canopy-air temperature, peak-to-peak over *consecutive* steps,
high-LAI sunlit fixture (`MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–10):

| `dt_fast` | 900 s | 450 s | 300 s | 225 s | 150 s | 100 s |
|---|---|---|---|---|---|---|
| CAS-T p2p | ~8 K | ~8 K | ~5 K | ~4 K | ~1 K | smooth |

Every conservation ledger closed to ~1e-6 J at every row. Default is 150 s; `meds_config` warns above
300 s.

**Pinning the entire surface coupling** — aerodynamic conductances, `g_tr_f` and `h_coeff_f` together
— takes 9.5 K to 6.7 K at 900 s. Pinning is strictly more aggressive than any refresh (it removes the
lag *and* the response), so per-sub-step refresh of the surface coupling cannot recover `dt_fast`.
Also eliminated, each still oscillating at 900 s: plant hydraulics (10.2 K), soil thermal column
(9.1 K), soil water column (10.0 K).

**Untested suspects for the residual `L`:** the leaf gas-exchange pre-pass (`gsw` solved once per
`dt_fast` at state `n`; pinning `g_tr_f` removed its *variation* while the pre-pass still ran at a
lagged state), the longwave emission base `t_emit`, and `f_wet`.

**A methodological gap that T4 closes.** Every experiment so far *pinned* a coefficient, which asks
"is X **necessary**?" — answered *no* for everything tried. It never asked "would refreshing X **cure**
it?" Under §1c those are different questions: pinning sets that coefficient's contribution to `L` to
zero *and* removes its physical response, whereas the perturbation decomposition in T2/T4 measures its
contribution to `L` directly.

Both schemes close whole-column water and energy to machine precision, with snow, condensation and
the tissue store active.

### 2b. What is void as evidence

Anything measured at `dt_fast` ∈ {900, 1800} s in a sunlit high-LAI regime is a **phase-sample of an
oscillation**, not a convergence measurement. That voids: the cross-scheme accuracy tables in the
retired parity document; "ARK time-stepping error 0.039 K vs total 1.10 K at 1800 s" and the 8×
refinement ratios; the deleted `dt_fast` → error table; ARK's rejection-storm counts; and the 30-year
`runs/ithaca_ark30/` beds.

### 2c. Code facts this plan assumes

Verified on `main` at `0f95e67`:

- `ark` is a 2-solve **ESDIRK2**, `f_E ≡ 0`. There is no explicit tableau branch to populate — a true
  IMEX is new machinery, not a flag (`ark2_column_step`).
- Inside the tableau: CAS enthalpy/humidity/CO₂ + the soil thermal column. **Operator-split out:** soil
  water (committed once from the scratch `column_hydrology_flux`), plant internal water mass, canopy
  surface water, plant hydraulics (an exact matrix exponential — deliberately not a stage).
- The leaf↔CAS pair is a direct 2×2 Newton with a **numerical** Jacobian (`jac_surface`):
  `NEWT_MAX = 4`, `LS_MAX = 6`, `FEVAL_CAP = 24`. A numerical Jacobian is why enlarging the implicit
  block (N2) is mechanical rather than a derivation exercise.
- `TISSUE_STORE_SCALE = 1.0` — the tissue store is on, carried as a b-weighted time integral, adding
  ~1.5e4 J/m2/K against the canopy air's ~3.0e4.
- Every state combinator (`state_init`, `state_axpy`, `state_accum`, `state_extrap`, `state_sub`,
  `state_err_diff`) allocates four per-cohort arrays on entry; one sub-step calls several. §6.
- The patch loop (`fast_dynamics`, `meds_fast_dynamics.f90:322`) is serial, and its per-patch scratch
  was deliberately hoisted *out* of the loop (BB1 phase 1) — exactly the structure that blocks
  threading. §7.
- `t_sub` depends only on `isub`, yet `met_advance` / `met_instant` / `apply_met_to_ctx` run **inside**
  the patch loop, repeating site-uniform work `n_patch` times per sub-step. §7 C1.
- CMake gives OpenMP flags to NVHPC only, and only on `meds_core`. The fast loop is in `meds_aux`.
- **Dead split residue:** `picard_max_iter` / `picard_tol_temp` / `picard_tol_shv` / `picard_relax` /
  `picard_fixed_iter` are plumbed to no reader; `budg%picard_iters` / `picard_nonconv` are read and
  never written; `[fast].ark_relax` is vestigial; `[fast].ark_niter` is only tested as `np <= 1`, so
  it is a boolean in disguise (the real cap is `NEWT_MAX = 4`).

---

# PART I — TESTING AND COMPARING INTEGRATORS

The model in §1 is a hypothesis fitted to five points. Part I is what turns it into a measurement, and
what any candidate scheme is then judged against. **Nothing in Part II is built before its gating
number exists.**

## 3. The measurement program

**Progress, 2026-07-31: T1, T2 and T4 are DONE — results in §1g.** They were run with a standalone
harness (`scratchpad/meas_fastloop{2,3,4}.f90`) against the built libraries, exploiting the fact that
`build_column_frozen` and `adaptive_ark_march` are public, so no source instrumentation was needed.
T3 and T5 remain outstanding. Two protocol errors were made and corrected on the way, both worth
recording because they are generic: a first pass spun up **separately per `dt_fast`**, so each row
probed a different state (`tau` varied 40× across rows and `Phi'` was non-monotone); and it timed
`adaptive_ark_march` in isolation with a **cold start**, so the isolated march did more work than the
warm-started in-context call and freeze + march exceeded the total.

### T1 — per-phase cost profile of one `dt_fast`  *(gates all of Part II)* — **DONE, §1g(iv)**

**Answer: the pre-pass is 56–89% of a sub-step at realistic cohort counts, and `canopy_aerodynamics`
is 2% of the pre-pass.** So cheapening the march is low-value, and the one coefficient that must stop
being frozen is also one of the cheapest to refresh.

Still worth building as permanent instrumentation, behind `[output].numerics`: wall-clock and call
counts for `column_prepass` (split out:
leaf gas exchange, aerodynamics, respiration), `canopy_radiation`, `solve_plant_water_batch`,
`column_hydrology_flux`, the ESDIRK march as a whole, and within it `surface_derivs` evaluations and
`soil_energy_step_implicit` solves. Report per patch and per site, summer and winter month, at
`dt_fast = 150 s`.

*Why it gates everything:* if `C_freeze` is 80% of the step, then a cheaper march is irrelevant, §6 is
noise, and threading's speedup is set by how well the *pre-pass* threads. If it is 30%, priorities
invert. This is the single most valuable missing number in the document.

### T2 — the stability multiplier `Phi'`  *(gates §5)* — **DONE, §1g(i)–(iii)**

**Answer: the boundary is ~75–80 s, so the default `dt_fast = 150 s` is unstable (`Phi' = −2.74`);
the map is differentiable; and the CAS↔atmosphere conductance carries 99% of the lag gain.**

Protocol, for reproduction: on the high-LAI sunlit fixture, at a
fixed step, perturb `bio%cas%can_enthalpy` by δ **before** `column_fast_step` (so `build_column_frozen`
sees the perturbed state — that is what makes it the *lagged* map), advance one `dt_fast`, and form

```
Phi'  =  ( H_{n+1}(H_n + delta) - H_{n+1}(H_n) ) / delta
```

Repeat at `dt_fast` ∈ {100, 150, 225, 450, 900} s and test §1f predictions 1–2. Then **decompose `L`**
by repeating with each frozen coefficient (`gah`, `g_tr_f`, `h_coeff_f`, `t_emit`, `f_wet`) held at its
*unperturbed* value — that isolates each one's contribution to the lag gain, which is the sufficiency
question §2a says has never been asked.

Deliverables: `Phi'(dt_fast)`; the fitted `L`; a per-coefficient decomposition of `L`; and the
predicted stabilising under-relaxation factor `θ_crit`.

**Secondary deliverable worth keeping permanently:** a `test_cas_stability` fixture asserting
`|Phi'| < 1` at the default `dt_fast`. The project currently has **no detector** for this failure
mode (§1e) — the ledgers, the embedded estimator and the Newton convergence test are all blind. A
regression test on the multiplier is the missing instrument.

### T3 — the accuracy baseline at a stable cadence — **DONE, §1h**

**Answer: with the refresh on, `ark` and `rk45` agree to 1.3 mK at 150 s; frozen, the "gap" is 0.49 K
and is an artefact of the oscillation. The freeze is ~97% of the remaining error, so integrator
tolerance is irrelevant. Refreshed `ark` converges at roughly second order; the frozen scheme has no
order across its stability boundary.** Protocol, for reproduction:

`ark` vs `rk45` at `dt_fast` ∈ {150, 100, 75, 50} s with `rtol_all`/`atol_scale` tightened on the
reference, scoring CAS temperature, soil-surface temperature, GPP and ET. There is currently **no**
valid quantitative statement of how far apart the two schemes are, or of what any `dt_fast` costs in
accuracy — §2b deleted them all. Guard the reference cell as `numerics_sweep.py` already does (the
reference `dt_fast` must divide both `dt_slow` and 3600 s), and add the `parity_fidelity.py`
record-count / duplicate-key guard.

### T4 — the sufficiency question — **DONE, §1g(ii)**, and it answers differently than expected

Asked as a *derivative* decomposition rather than an amplitude sweep: build the frozen bundle at the
perturbed state, then overwrite one coefficient with its unperturbed value. **`gah`/`gaw`/`gac`:
`Phi'` −23.2 → +0.80. Everything else: ≤0.7%.** The earlier "pinning" experiments pinned `ustar`,
which leaves the stability-dependent transfer factor `temp1` free and therefore removes only about
half of `gah`'s sensitivity — which is why they concluded the aerodynamics merely amplifies.

**Still worth running as an amplitude confirmation** once N2 exists: the ED2 arm (refresh everything
per stage) as an upper bound, and a `dt_fast` amplitude scan with only `gah` refreshed, to confirm
that a multiplier below 1 actually produces a smooth trace rather than a differently-shaped one.

### T5 — the stability threshold where it is worst — **DONE (a), §1h**

**Answer: four of five stands are unstable at the default `dt_fast = 150 s`, and the worst is a
MID-HEIGHT canopy (10 m, LAI 3, `Phi' = −8.2`), not a regenerating one — prediction #4 was backwards.
The frozen scheme would need `dt_fast ≈ 50 s` to be safe there.** This is the correctness issue the
item was written to catch, and it outranks every optimisation in Part II: it means the *current*
default is unsafe, not merely slow. With the refresh, every stand is stable or marginal at every
cadence measured.

(b) — the tissue-store contribution — is **not run and is now moot**: with N2a the canopy air is
stable at every cadence measured, so separating the store's ~50% capacity contribution no longer
gates anything. Re-open only if N2a is rejected.

### T6 — evidence-base requirements for any cell reported

One site, one month per cell, one forcing year, no confidence intervals is the current standard and it
is too weak to rank anything within ~2×. Any claim promoted out of Part I carries: a dry-down window,
a cold season with snow, and a wet saturated window, in addition to the summer/winter pair; both
compilers; and `[output].numerics = true`.

## 4. How a candidate integrator or stabilisation is judged

A Part II option may become the default only when **all** hold:

1. It raises the measured stable `dt_fast` by ≥ 2× on **both** the high-LAI sunlit fixture and the
   young-stand fixture (T5), at CAS-T p2p ≤ 1 K — or, if it is a cost option rather than a stability
   option, it leaves `Phi'` unchanged within measurement error.
2. Net fast-loop cost per simulated day falls against today's 150 s configuration on the T1 profile.
   A stability win that costs more per step than it saves in steps is not a win.
3. **It does not change the limit.** The T3 refinement curve still converges, and to the same answer
   within the stated tolerance. A damping scheme that converges somewhere else is a *different model*
   — the lesson `split` cost this project, and it must not be repeated silently. This is the criterion
   most likely to fail for N1.
4. Every whole-column ledger still closes to machine precision with snow, condensation and the tissue
   store active, `debug_error = true`, on both schemes — **and** a consecutive-step CAS-T trace is
   inspected, because conservation is not stability.
5. Full `ctest` green on ifx Release, ifx Debug and nvfortran multicore (the issue-#7 rule: a green ifx
   run is not sufficient).
6. T6's scenario coverage.

A 30-year Ithaca run at the shipped configuration completes, conserves, and its AGB trajectory is
explicable. It is a **new baseline, not a comparison**: `runs/ithaca_ark30/` predates the recycle-phase
fix (PR #69), the snowfall routing fix, PR #88 and the `dt_fast` change.

---

# PART II — NUMERICAL IMPROVEMENT AND MULTI-CORE

## 5. Raising the stable `dt_fast` — the lever that cuts both cost terms

T2/T4 have run, so this section is no longer a menu of guesses — it has one indicated fix and a set
of alternatives that the measurement demoted.

**N2 — Unfreeze `gah`/`gaw`/`gac`. THE headline item.** *(indicated by §1g(ii)+(iv); do this first)*

The CAS↔atmosphere conductance carries 99% of the lag gain and costs 2% of the pre-pass. Refreshing
it — and nothing else — is predicted to take `Phi'` from −23.2 to +0.80 at `dt_fast = 900 s` for ~3%
more work per sub-step. If that holds end to end it is up to a **6× cut in fast-loop cost**, which is
larger than everything else in this document combined.

> ### N2a IS IMPLEMENTED AND MEASURED — 2026-07-31. It works, and by more than predicted.
>
> `[ccfg].aero_refresh` (default **0 = frozen**, bit-identical; `1` = refresh per ESDIRK stage).
> `refresh_cas_conductances` in `meds_fast_ark.f90` re-solves **only** `mo_surface_layer` — the bulk
> MO solve — at the stage's own `(T_cas, q_cas)`, leaving every per-cohort quantity frozen. The
> refreshed `gah/gaw/gac` are the *local* variables `column_be_stage` already used for the BE commit,
> the Newton and the boundary-flux ledger, so "one flux, both sides" holds with no extra work.
>
> **Stability — the oscillation is gone at every `dt_fast` tested.**
>
> | `dt_fast` | 25 s | 50 s | 75 s | 100 s | 150 s | 225 s | 450 s | 900 s |
> |---|---|---|---|---|---|---|---|---|
> | `Phi'` frozen | +0.36 | −0.26 | −0.88 | −1.49 | −2.74 | −4.64 | −10.6 | −23.2 |
> | `Phi'` refreshed | +0.35 | −0.21 | −0.48 | −0.74 | **−0.79** | **−0.67** | **−0.42** | **−0.14** |
> | CAS-T p2p frozen [K] | 0.12 | 0.03 | 0.25 | 0.73 | 2.01 | 4.04 | 6.60 | 7.69 |
> | CAS-T p2p refreshed [K] | 0.12 | 0.08 | 0.04 | 0.02 | **0.012** | **0.015** | **0.032** | **0.100** |
>
> `|Phi'| < 1` everywhere, and it *decreases* with `dt_fast` — the residual lag stops mattering as the
> storage term shrinks. Peak-to-peak at 900 s falls **77×**, from 7.7 K to 0.10 K.
>
> **Production economics — a full 24 h day, 5 cohorts, scored against `aero_refresh=1 @ 25 s`:**
>
> | cell | wall | RMSE CAS-T | RMSE soil-T | sub-steps | rejections |
> |---|---|---|---|---|---|
> | **frozen @ 150 s (today's default)** | 24.8 ms | 0.390 K | 0.079 K | 931 | 4 |
> | refreshed @ 150 s | 22.5 ms | **0.017 K** | **0.063 K** | 576 | **0** |
> | refreshed @ 225 s | 18.2 ms | 0.036 K | 0.145 K | 664 | 1 |
> | refreshed @ 450 s | 10.1 ms | 0.085 K | 0.372 K | 348 | 4 |
> | refreshed @ 900 s | **6.4 ms** | 0.160 K | 0.701 K | 239 | 39 |
>
> At **matched** `dt_fast = 150 s` the refresh is a strict improvement on every axis at once: 9%
> cheaper, 23× better on canopy-air temperature, better on soil temperature, and it removes the
> rejections. Against today's default, refreshed @ 900 s is **3.9× cheaper and 2.4× more accurate on
> CAS-T**; refreshed @ 450 s is 2.5× cheaper and 4.6× more accurate.
>
> **It does not move the limit (§4 criterion 3).** At `dt_fast = 50 s` the *frozen* scheme scores
> 0.0016 K against the refreshed reference — the two semi-discretisations converge to the same answer,
> so this is a consistency improvement, not a different model.
>
> **Conservation holds.** `n_fail = 0` across all six ledgers in every cell. Full `ctest` is green on
> ifx Release **and** nvfortran multicore both with the flag off *and* with it forced on (36/36 each).
>
> **The cost model needs correcting, in the plan's favour.** §1g(iv) predicted ~3% overhead from the
> refresh itself. Measured at a single midday probe it looked like +50%, but that was the *march*
> taking an extra sub-step at that one state. Over a whole day the refresh **reduces** integrator work
> at every `dt_fast ≥ 150 s` (576 vs 931 sub-steps at 150 s), because a non-oscillating column is an
> easier problem for the error controller. The refresh pays for itself before any `dt_fast` change.
>
> **Soil temperature is now the binding constraint.** CAS-T error is nearly `dt_fast`-independent once
> the feedback is live, but soil-T error still grows (0.06 → 0.70 K over 150 → 900 s) because the
> *ground* fluxes remain frozen. So the recommendation is **`dt_fast` 225–450 s**, not 900 s, until
> something equivalent is done for the ground surface — that is the natural next target, and §1g's
> decomposition method applies to it unchanged.
>
> ### RK45 WIRED — 2026-07-31. The two schemes now agree to 1%.
>
> The family asymmetry is closed. RK45's CAS tendency is built in `surface_derivs` (not
> `column_derivs`), so that is where its refresh lives; `stage_bnd` re-solves the same `pure` kernel
> at the same stage state so its boundary flux is charged at the conductance the tendency used.
>
> **ARK pays for one surface-layer solve per stage, not one per residual evaluation.** `surface_derivs`
> is called up to 24 times per stage by the Newton, to fill a CAS tendency ARK never reads — it commits
> the CAS through its own backward-Euler denominator. So `column_be_stage` refreshes once, writes the
> result into its stage-local `fs`, and **clears `fs%aero_refresh`**, after which `surface_derivs` reads
> the already-refreshed values. RK45, whose every RHS evaluation *is* one stage, keeps the flag set.
>
> **5 cohorts, 24 h day, scored against `ark`/refresh/25 s:**
>
> | scheme | refresh | `dt_fast` | wall | RMSE CAS-T | RMSE soil-T | sub-steps | rescue | `n_fail` |
> |---|---|---|---|---|---|---|---|---|
> | ark | off | 150 s | 39.8 ms | 0.390 K | 0.079 K | 931 | 0 | 0 |
> | ark | **on** | 150 s | 25.3 ms | **0.0171 K** | 0.063 K | 576 | 0 | 0 |
> | rk45 | off | 150 s | 32.6 ms | 0.375 K | 0.077 K | 576 | 0 | 0 |
> | rk45 | **on** | 150 s | 41.5 ms | **0.0169 K** | 0.061 K | 838 | 0 | 0 |
>
> **ARK and RK45 land within 1% of each other** (0.0171 vs 0.0169 K) — an implicit and an explicit
> method agreeing, which is the strongest available evidence that this is a consistency fix rather
> than a scheme-specific tweak. `n_fail = 0` on all six ledgers in every cell, which is also what
> validates the `stage_bnd` plumbing: had the ledger used a different `gah` than the tendency, whole-
> column energy would have failed immediately.
>
> **The refresh helps ARK's cost and hurts RK45's** — sub-steps 931 → 576 on ARK, 576 → 838 on RK45.
> Expected, and it sharpens the tier split: making `gah` state-dependent puts a strong negative
> feedback into the RHS, which an implicit method *absorbs* (the column stops oscillating, so the
> controller relaxes) and an explicit method must *resolve* (its stability limit tightens). At 900 s
> refreshed, `ark` is 2.3× cheaper than `rk45` and slightly more accurate. `rk45` is the accuracy
> baseline, not the production tier, so paying there to buy consistency is the right trade.
>
> **Bit-identity of the default path: verified, after initially losing it.** The first version reported
> the used conductances on `surface_tend_t` so `stage_bnd` could read them back — the more obviously
> safe design. Measured against a `git stash` of the pre-change tree it was **not** bit-identical:
> ~1e-12 relative, because three extra fields on that type perturb ifx's inlining/FMA decisions inside
> `surface_derivs`. Harmless numerically, but it silently shifted derived RMSEs by up to 8% in the
> cells where the frozen map is unstable and round-off gets amplified, and it would break the
> stash-and-`cmp` protocol the project verifies with. Removing the fields and re-solving in
> `stage_bnd` restores **bit-for-bit** identity against the pre-change tree, confirmed by `cmp` on a
> 288-sample trajectory across both schemes at three `dt_fast`.
>
> ### LANDED 2026-07-31 — the flag is deleted and the refresh is unconditional.
>
> User decision, on the T3/T5 evidence: **delete the switch, always refresh.** The frozen alternative
> is unstable at the production cadence on four of five stand heights (`Phi'` to −8.2), so it is
> known-wrong physics rather than a supported configuration — the same reasoning that deleted
> `snow_on` and the `with_theta` norm switch.
>
> `surface_frozen_t%mo_live` is what remains, and it is **not** that switch: it records whether the
> `mo_*` inputs are populated, and therefore whether a re-solve is possible at all. Three uses, none
> of them a user choice — default `.false.` ("use `gah/gaw/gac` as given", which is what a hand-built
> unit-test bundle means); `build_column_frozen` sets `.true.`; ARK's stage copy clears it after
> refreshing once, so the Newton's ≤24 residual evaluations per stage do not each re-solve.
>
> **That default direction was learned the hard way.** With the marker defaulting to "refresh",
> `test_column_derivs` — which builds a `surface_frozen_t` by hand and never populates the `mo_*`
> fields — called the surface-layer solve with zero roughness and zero wind and **hung**. The safe
> default is the one a bundle that knows nothing about this mechanism will get.
>
> `dt_fast` default moved 150 s → **900 s** with it, per the same decision. See the warning below.
>
> Validation: `ctest` 36/36 on ifx Release and nvfortran multicore.

### The carbon bias at long `dt_fast` — measured while landing the above, and it constrains the default

Moving `dt_fast` to 900 s costs far more in carbon than in temperature. Against a 12.5 s reference on
the high-LAI sunlit stand:

| `dt_fast` | 150 s | 225 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| GPP | −3.8% | −8.4% | −12.6% | −19.8% | **−33.1%** |
| ET | −2.6% | −5.8% | −8.7% | −14.5% | **−23.8%** |
| CAS-T RMSE | 0.017 K | 0.037 K | 0.055 K | 0.085 K | 0.16 K |
| soil-T RMSE | 0.06 K | 0.14 K | 0.22 K | 0.37 K | 0.70 K |

**Attributed in two steps, and the second changes the conclusion.**

*Step 1 — it is the plant hydraulics.* Freezing the plant-hydraulics store (`mask%hydraulics`) removes
the GPP bias entirely — **−33.1% → +0.93% at 900 s, and flat in `dt_fast`**.

*Step 2 — it is not a quadrature bug, and the stress limb is an amplifier rather than the source.* GPP
is accumulated as `Σ gpp(state n)·dt`, a left rectangle in the state, so an averaging artefact was the
obvious suspect. Disabling the **non-stomatal (capacity) water-stress limb** — a linear ramp
`beta = (ψ − ψ_close)/(ψ_open − ψ_close)` on Vcmax/Jmax/TPU — settles both questions at once:

| `dt_fast` | 12.5 s | 150 s | 300 s | 450 s | 900 s |
|---|---|---|---|---|---|
| GPP, limb on | 1.000 | 0.962 | 0.874 | 0.802 | **0.667** |
| GPP, limb **off** | 1.000 | 1.000 | 1.000 | 1.000 | **0.9995** |
| ET, limb off | 1.000 | 1.000 | 1.000 | 0.993 | 0.988 |
| daytime-mean ψ_leaf [MPa] | −0.23 | −0.50 | −0.72 | −0.88 | **−1.19** |
| max per-step \|Δψ\| [MPa] | 0.0001 | 0.010 | 0.034 | 0.071 | **0.84** |
| fraction of daylight steps with beta < 1 | 0.00 | 0.61 | 0.77 | 0.83 | 0.91 |
| fraction with beta clamped at 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |

So: **the quadrature is sound** (0.05% at 900 s with the limb off), **ψ_leaf itself is badly
non-converged** (−0.23 → −1.19 MPa), and **the limb converts that linearly into carbon** — slope
~0.5 MPa⁻¹ on the shipped PFT file, so ~1 MPa of ψ error ≈ 50% of Vcmax. The `beta = 0` clamp is never
reached, so this is straight linear amplification, not a saturation/Jensen artefact.

**Decision (issue #47, user, 2026-07-31): the capacity limb is OFF by default.**
`[leaf_physiology].wstress_nonstomatal`, default `.false.`. It is rarely measured directly and its two
parameters are weakly constrained, so wiring an unconverged ψ into carbon through it is not a trade
worth making. The stomatal limb (ψ_soil → `g1`/λ) is better constrained and stays on. The kernel keeps
the term and `test_leaf_physiology` still exercises it explicitly, plus a new assertion that ψ_leaf is
**inert** in the default configuration.

**Consequence for the `dt_fast = 900 s` default: it is now defensible.** With the limb off, 900 s costs
0.05% of GPP and 1.2% of ET. The `meds_config` warning is conditional on the limb accordingly.

**But the ψ error is unfixed — it is only unwired.** Anything keyed to leaf water potential still
inherits it, and it will come back the moment ψ feeds anything else (hydraulic mortality, ψ-driven
allocation, or the limb being re-enabled for a study). That is N2b, and it is still the top item.

### N2b — the leaf-water-potential lag. NEXT, and it is what unlocks a long `dt_fast`.

`psi_leaf` is diagnosed once per `dt_fast` from `leaf_water_mass` at state `n` (`column_prepass`) and
feeds `g_sw` through `leaf_gas_exchange_batch`. Unlike the conductances, re-solving the *whole* leaf
kernel per stage is not affordable — it is ~89% of the pre-pass at 30 cohorts. So the candidates are
the cheaper ones §5 already lists, now aimed at this seam:

- **tangent-linear `g_sw`** — evaluate `∂g_sw/∂ψ` (and `∂ψ/∂W`) once per `dt_fast` alongside the
  kernel call that is already happening, then correct at stage points. The mass update is already a
  closed-form Euler step, so ψ at a stage point is available for free.
- **advance ψ within the step** rather than only at its end, so the *frozen* value the next step
  inherits is a step-average rather than an endpoint.
- verify first, by the §1g method: perturb `leaf_water_mass`, measure `∂GPP/∂W`, and confirm the loop
  gain is what the −33% implies before building anything.

Three implementations, cheapest first:Three implementations, cheapest first:

- **N2a — per-stage refresh. ✅ IMPLEMENTED (see above).** Re-solve `mo_surface_layer` inside
  `column_be_stage` at the stage's own CAS state, leaving every other frozen quantity alone. This is
  what ED2 does (`update_diagnostic_vars` → `canopy_turbulence8` at every RK stage). Deliberately
  *not* a full `canopy_aerodynamics` call — the per-cohort boundary layers stay frozen, since they
  carry ≤0.7% of the lag gain and are what makes a full refresh expensive.
- **N2b — tangent-linear `gah`.** Evaluate `∂gah/∂H` once per `dt_fast` by one extra
  `canopy_aerodynamics` call and apply `gah(H) ≈ gah_n + (∂gah/∂H)(H − H_n)` at stage points. Cheaper
  than N2a, degrades to today's behaviour as `H → H_n`, and is exactly the "reuse the Jacobian, never
  the residual" discipline. But `d ln gah/dT ≈ 2.2 K⁻¹` is a *strong* nonlinearity, so a linear model
  over a multi-K excursion may be poor — measure `Phi'` under N2b before preferring it to N2a.
- **N2c — `gah` inside the Newton.** Add `∂gah/∂H` to `jac_surface`, making the ventilation term
  implicit rather than merely refreshed. Strictly the most robust (it cannot overshoot), and the
  Jacobian is already numerical so it costs one extra `surface_derivs` evaluation per iteration.

**Watch items for all three.** `gah` appears in the CAS commit denominator `(wcap + dt·gah)` and in
the boundary-flux ledger terms (`bf%cas_enth_in/out`, `whole_enth_out`), so a per-stage `gah` must be
the *same* value in the state update and in the ledger or conservation breaks — the "one flux, both
sides" discipline the ARK already relies on. And the `ustar` floor (0.1) is currently an accidental
stabiliser at very light wind (§1g(iii)); once the feedback is live, check that removing the lag does
not simply move the chattering onto the floor.

**N1 — Under-relax the frozen coefficients. DEMOTED.** `c_{n+1} = (1−θ)·c_n + θ·F(y_n)` gives
`Phi' = a + (1−a)θL`. With the measured `L ≈ −62` at 900 s this needs **`θ < 0.07`** — a coefficient
time constant of ~15 steps, which would lag sunrise and rain onset by hours and is nearly certain to
fail §4 criterion 3. At 150 s (`L ≈ −38`, `a = 0.91`) `θ ≈ 0.5` would suffice, so N1 remains viable as
a *cheap safety net at the current default* if N2 is delayed — but it cannot buy a longer `dt_fast`,
which was its entire attraction. Prediction 3 of §1f is withdrawn.

**N3 — Rosenbrock / W-method (ROS3P, RODAS3).** Linearly implicit: one Jacobian per step, linear
solves only, no nonlinear iteration. The relevant property is that **W-methods tolerate an inexact
Jacobian by design** — which legitimises carrying an approximate `∂c/∂y` (the exact thing N2 adds by
finite differences) without harming order. L-stable and stiffly accurate variants exist at order 3.
This is the best *formal* fit to the problem as §1 characterises it, and the natural landing place if
N2 works but its Newton cost is unattractive. Larger change than N2; do not start here.

**N4 — A self-consistent algebraic (index-1 DAE) canopy air.** The correct form of "make it
diagnostic": fast variables (CAS twins, leaf and wood temperature, ground skin) solved as an algebraic
block with coefficients evaluated at the same state, slow variables (soil heat, soil water, plant
water) integrated. This is CLM5's structure. It is N2 taken to completion, and it should be reached
*through* N2 rather than as a separate build. Two standing requirements: the CAS storage term appears
in the ledgers (`wcap*enth0`/`wcap*enth1` and the `cas_*` budgets in `column_fast_step_ark`), so those
must be re-derived rather than patched; and `wcap` stops being negligible at night and in calm
conditions, so validation must happen there, not at midday.

**N5 — An adaptive freeze cadence, with a real error estimator.** `tau` spans ~80 s (windy midday) to
~2 h (calm night), so the stability limit varies by more than an order of magnitude within a day, and
`dt_fast = 150 s` is set by the worst hour and paid for all 24. Per §1e the controller cannot find this
because it has no estimator for the freeze — so build one: the next step's pre-pass computes
`c(y_{n+1})` anyway, making the drift `||c(y_{n+1}) − c(y_n)||` (or the induced `a`, from `wcap` and
`gah`, both already on `fro%surf`) free. Choose `n_fast_per_slow` per slow step from the previous
step's indicator, keeping `dt_slow/dt_fast` an exact integer. Long steps at night and in calm periods,
short in windy sun.

This is the first proposed exception to the layer-3 rule that `dt_fast` is user-set, and it is a
narrow one: the *cadence* adapts, the *scheme* does not, and the quantity controlled is a measured
error indicator rather than a heuristic. It may well be the largest production win available, since it
recovers the cost of the worst hour across the other 23. Gated on T2 supplying an indicator that
actually tracks `Phi'`.

**N6 — Refuted or inverted; do not rebuild.**

- *A true IMEX tableau built for its own sake.* The general point stands: moving a term from *frozen*
  to *explicit* costs **more** evaluations per step, so a true IMEX can never cheapen the current step
  — its only route to a win is through `dt_fast`. **But note that N2a is the useful part of that idea,
  now targeted.** The first draft's A1/A2 proposed moving the *whole* surface coupling explicit, at a
  cost that scales with cohort count (leaf gas exchange, 89% of the pre-pass at 30 cohorts); §1g says
  only `gah`/`gaw`/`gac` need to move, at 2%. The refuted item is the *scope*, not the mechanism.
- *Partial refresh of the surface coupling per sub-step, as previously scoped* — the earlier
  refutation was measured with `ustar` pinned, which does not pin `gah` (§1g(ii)); it did not test
  what it claimed to test. What is refuted is refreshing the **expensive** parts.
- *Naive `wcap → 0`* — §1d: maximally unstable.
- *Coefficient extrapolation* (`2c_n − c_{n−1}`) — cancels the lag to first order but amplifies an
  already-oscillating state. Under-relaxation is the safe direction; the asymmetry is real. Revisit
  only after N1/N2 have restored contraction.
- *Inner sub-stepping as a stability measure* — §1d: slightly harmful.

## 6. Cheaper steps

Ordered by expected value once T1 says where the time goes:

1. **Allocation churn in the state combinators.** Six routines allocate four per-cohort arrays each on
   entry; one adaptive sub-step calls several. A pre-allocated workspace carried through the march (or
   fixed capacity sized to the site-wide max cohort count, as the driver's scratch already is) removes
   it. Also a §7 prerequisite — concurrent threads hammering the allocator is a classic scaling wall.
2. **`bflux_zero`'s per-call allocation** of the tissue integrals — same treatment.
3. **The RK45 rescue snapshot** copies whole `patch_biophys_t` / `column_budget_t` values per sub-step
   (`meds_fast_step.f90:101`) whether or not a rescue fires. Baseline tier only; measure first.
4. **A cost ceiling wired to a diagnostic, not a control-flow switch.** If a `dt_fast` needs more than
   *N* sub-steps, report it through the work counters and halt under L2. The machinery exists
   (`RK45_WORK_CAP`, the `dt_fast/64` floor); the reporting path does not.
5. **Delete the split residue** (§2c): the `picard_*` surface and its always-zero counters; give
   `ark_niter` real meaning (wire it to `NEWT_MAX`) or make it an honest boolean; resolve `ark_relax`
   (N1 may claim it). These make `[fast]` read as if it offers cost knobs it does not.

## 7. Multi-core

The only lever that costs no accuracy — **provided output is bit-identical regardless of thread
count**. That is a hard requirement, not a nicety: verification here is byte-for-byte netCDF
comparison (`MEDS_NUMERICS_SCOPING.md` §11.2), and a non-deterministic reduction destroys it.

**C0 — the axis.** Patches within a site: columns are independent within a `dt_fast`, coupled only
through the slow loop. Cohorts within a patch are too fine (and the batch kernels are the natural
vectorisation target); sub-steps are sequential; polygons/sites are the real long-run scaling axis but
the runtime is single-polygon today and that is out of scope. Speedup is bounded by `n_patch` and by
**load imbalance** — each patch runs its own adaptive march — so `schedule(dynamic, 1)`.

**C1 — hoist the met sampling out of the patch loop.** `t_sub` depends only on `isub`, yet
`met_advance` (file reader, `intent(inout)`, may reload a netCDF bracket), `met_instant` (solar
geometry, shortwave disaggregation, Weiss–Norman) and `apply_met_to_ctx` all run inside the patch
loop. Precompute the `n_fast_per_slow` samples once per slow step and index them. Correctness-neutral,
independently a win, **and it removes the file reader from the region about to become parallel.**
First, on its own commit, verified byte-identical.

**C2 — per-thread scratch.** BB1 phase 1 hoisted `coh`/`bio`/`aero`/`forc`/`budg`/`gpp_coh…` out of the
loop and sized them to the site-wide max cohort count, to cut `O(n_patch)` allocations per slow step.
That optimisation and threading are in direct conflict — those buffers are the loop-carried shared
state. Keep the hoist but make it per thread (a pool of `n_thread` copies, indexed by thread id), and
verify the allocation count does not regress; that was BB1's whole point.

**C3 — deterministic reductions.** Classify every write:

- *Disjoint, safe as-is:* `site%patch%{cas,soil_e,soil_w,snow}(ip)`, `site%patch%xi_accum(ip)`, and the
  per-cohort slices `site%cohort%*_accum(i0:i0+ncoh-1)`.
- *Reductions:* `site%et_accum`, `site%pheno_tair_sum/n`, the nine `site%work_*`, and local
  `we`/`ww`/`nfail`. OpenMP `reduction(+:…)` sums in **thread-arrival order**, so the answer moves with
  thread count. Accumulate into per-patch arrays and sum them **in patch order** after the loop.
- *Output staging:* `mgr%fast(isub)%…` and `mgr%fast_soil_*` are area-weighted over patches — same
  treatment. `mgr%fast_coh_*` is written by global cohort slot and is disjoint.
- *Not thread-safe:* `write_fast_probe` holds `save`d `unit`/`opened` and writes a shared file. Either
  serialize it or make `fast_probe` and threading mutually exclusive at config time, and say why.

**C4 — build wiring.** `-qopenmp` (IntelLLVM), `-fopenmp` (GNU), `-mp` (NVHPC multicore), applied to
the target owning `meds_fast_dynamics`, threading **off by default** so no existing result moves
without opt-in. Expose it once — a `[run].n_threads` key *or* the environment variable, not both
silently.

**C5 — validation.** Byte-identical diagnostic netCDF at 1/2/4/8 threads; allocation count not
regressed against BB1; `ctest` green on all three back ends; speedup reported against the measured
parallel fraction from T1, not against wall-clock alone.

## 8. Deliberately not in this plan

- **A higher-order tableau.** The freeze dominates the error at production `dt_fast`; an ARK4(3) over a
  first-order semi-discretisation buys nothing. Revisit only if §5 makes the tableau binding.
- **Folding plant hydraulics into the tableau.** Measured 4×–170× *worse*, and a category error: ψ is
  advanced by an exact matrix exponential and the coupled subsystem is ψ-independent within a step. ED2
  makes the same choice for the same reason.
- **Folding soil water into the ARK tableau.** Real work — it would retire the structurally-zero
  dilution in the ARK's error norm — but it buys correctness, not speed. Track separately.
- **GPU offload of the fast loop.** Gated on the bare-array kernel conversion
  (`MEDS_NUMERICS_SCOPING.md` §11). Note the ordering: §7's per-thread scratch and deterministic
  reductions are prerequisites for offload too, so C is not wasted if the GPU path resumes.
- **MPI / multi-polygon parallelism.** The right long-run scaling axis, explicitly not this document's.
- **Re-litigating `split`.** Archived with a tombstone; not evidence for anything here.

## 9. Where the code is

| piece | file |
|---|---|
| dispatch + the RK45→ARK stiff rescue | `src/driver/meds_fast_step.f90` (`column_fast_step`) |
| ESDIRK2 scheme (config name `ark`) | `src/driver/meds_fast_ark.f90` (`ark2_column_step`, `adaptive_ark_march`) |
| the frozen pre-pass — §5 acts here | `src/driver/meds_fast_ark.f90` (`column_prepass`, `build_column_frozen`) |
| the per-stage implicit surface solve — N2 acts here | `src/driver/meds_fast_ark.f90` (`column_be_stage`, `newton_surface_solve`, `jac_surface`) |
| the shared right-hand side | `src/driver/meds_fast_time_derivs.f90` (`surface_derivs`, `column_derivs`) |
| adaptive Cash–Karp scheme (baseline tier) | `src/driver/meds_fast_rk45.f90` |
| tolerances, error norm, step controller | `src/driver/meds_fast_control.f90` |
| the patch loop — §7 acts here | `src/driver/meds_fast_dynamics.f90` (`fast_dynamics`) |
| cadence owner + the `dt_fast` warning | `src/driver/meds_stepper.f90`, `src/shared/config/meds_config.f90` |
| benchmark + scoring harness | `scripts/numerics_sweep.py`, `scripts/parity_fidelity.py` |
