# The frozen-seam contract

**Status: a design note. No code.** Written 2026-09-13 (#201), to the specification in
`MEDS_CODE_STRUCTURE_DESIGN.md` §15.4. Its value is conceptual: the instruments it describes already
exist and are already reported. What was missing was the *rule* for reading them, so that classifying
the next seam somebody adds is a lookup rather than a re-derivation.

---

## 1. What a frozen seam is

The slow tier freezes a state across a whole slow step while the fast loop computes fluxes from it.
That pattern is everywhere in MEDS, and it is the same pattern the fast loop applies internally when
`build_column_frozen` holds photosynthesis, radiation and hydraulics fixed across a `dt_fast`.

It is an approximation, and it is a *good* one — but only under a condition that is checkable, and
which nothing in the code states.

## 2. The Λ criterion

For a frozen store $`S`$ consumed by a flux $`F`$ over a step $`dt`$, the governing number is the
fraction of the store the step withdraws:

```math
\Lambda = \frac{F\,dt}{S} \qquad(1)
```

Λ is dimensionless and scale-free. **Everything turns on what $`F`$ does as $`S \to 0`$**, and there
are exactly two cases:

**Case A — $`F = kS`$, linear in the store.** Then $`\Lambda = k\,dt`$, *independent of $`S`$*. The
freeze is sound at any store size: near-bare ground is no more dangerous than a mature soil. This is
the CENTURY case. `lambda_max = max_j (xi_int_j · K_j)` measures **3.6×10⁻³ per day** on the mature
Ithaca stand — a daily step withdraws about 0.4 % of the fastest pool — and it is *why bare-ground
spin-up works at all*.

**Case B — $`F`$ prescribed, or otherwise not vanishing with $`S`$.** Then Λ diverges as
$`S \to 0`$, and the freeze is not merely inaccurate, it is **unsound**: the step withdraws from a
store that is not there.

MEDS has a scar from exactly case B. The old `soil_carbon_on = .false.` branch respired a
*prescribed* 5 kgC m⁻² pool that did not exist. The seven pools were identically zero and received no
litter, so nothing was ever debited for the CO₂ leaving the column. It held the canopy air 47 ppm
above the 420 ppm datum, drove site NEE to **+30.7 µmol m⁻² s⁻¹** where the true answer is ~0, and
the elevated CO₂ then fertilised photosynthesis and inflated GPP by 5 %. None of it was visible,
because `rh_site` read the CENTURY matrix and reported 0 throughout.

**The rule:** a store seam is admissible when its flux is linear in the store, or when Λ is bounded
below 1 by an argument that does not assume the store stays large. A prescribed flux against a
prognostic store is never admissible.

## 3. The four seams, classified

| seam | what is frozen | flux form | verdict |
|---|---|---|---|
| `soil_carbon` → fast Rh | the seven CENTURY pools, held across the day | $`F = kS`$ | **Case A, sound.** Λ = 3.6×10⁻³/day, reported as `lambda_max` with the attaining pool in `lambda_pool`. |
| `xi_accum` (the return accumulator) | nothing — it is the *return* path | n/a | Not a store seam. Its check is `rh_seam_gap`, §4. |
| `shed_water_rate` | a frozen daily **rate** | n/a | **Rate seam, §4.** Λ is meaningless. |
| `slow_co2_rate` | a frozen daily **rate** | n/a | **Rate seam, §4.** Λ is meaningless. |

## 4. Rate seams are a different contract

Λ has no denominator for a frozen *rate*: there is no store in it. What makes a rate seam sound is a
sequencing invariant instead —

> **the slow tier debits its own store before it credits the rate the fast tier will consume.**

`shed_water_rate` and `slow_co2_rate` are both safe on that basis: the water and carbon have already
left the plant pools by the time the fast loop adds them to the ground. Today that is a **convention
stated in comments, not an assertion.** Making it one is the natural follow-up, and is deliberately
not claimed here.

The check that *does* exist for the return path is `rh_seam_gap` = `rh_out − rh_fast_accum`, which is
zero by construction on an ordinary day and non-zero when patch structure changes between the fast
window and the slow step (#192, documented at the field).

**Λ and the seam gap are not substitutes.** The seam gap catches a **broken contract** — it fires
after something is already wrong. Λ catches a **degrading approximation** *before* anything breaks.
Λ is the better instrument and the generalisation of the two.

## 5. Arbitration: when two consumers share one frozen store

Soil water has four consumers — transpiration, ground evaporation, drainage and runoff — and each was
historically clamped on its own. **Per-process clamping is order-dependent**: whichever process runs
first gets its full demand and the last one absorbs the shortfall, so the answer depends on call
order rather than on physics.

The rule is to **arbitrate, not clamp**: form the total demand $`D`$ over all consumers and scale
every demand by

```math
\min\!\left(1, \; \frac{S}{D}\right) \qquad(2)
```

so each consumer is reduced in proportion. This is order-independent and conserves the store exactly.

## 6. Why this note exists: two defect classes that came out of not having it

- **Borrowing one solve's flux while committing another's state.** Each fast-loop scheme freezes a
  scratch hydrology solve and then advances the column; wherever a scheme took a *flux* from the
  scratch while committing *state* from its own stages, the two disagreed by however far the
  trajectories had diverged — invisible while they nearly coincide, dominant once the column
  saturates. Three instances were found, all on RK45, two of matched magnitude and *opposite* sign so
  that each concealed the other. **A whole-column ledger cannot detect this class at all**, because
  the error is purely vertical: enthalpy in the wrong layer still sums correctly against the
  boundary.
- **The frozen conductance that drove a canopy-air oscillation.** Evaluating the CAS↔atmosphere
  conductances a whole `dt_fast` behind the state they respond to turned a stabilising feedback into
  a sustained period-2 oscillation, up to ~8 K step-to-step at 900 s, **while every conservation
  budget closed to ~10⁻⁶ J throughout**. See `docs/science/numerical_scheme.md` §2.

Both are seam questions that a stated rule turns into a review question rather than a discovery.
Neither is caught by conservation, which is the point: **conservation proves bookkeeping, not
plausibility.**

## 7. Explicitly not recommended

**Making the seam implicit.** `SOIL_LIN_PICARD` exists, but paying an outer iteration to fix a
3.6×10⁻³ error is the wrong trade.
