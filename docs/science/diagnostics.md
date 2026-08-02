# MEDS diagnostic output

How MEDS reports itself: what a variable means, how it is aggregated across the demographic
hierarchy and across time, and how to ask for the ones you want.

Implementation: `src/io/` (`meds_diagnostic_kernels`, `meds_diagnostic_reduce`,
`meds_output_{types,registry,integrate,stream,manager}`) plus the per-cohort / per-patch
accumulators in `src/core/meds_core_diag_types.f90`. Design and rationale:
`docs/dev_plans/MEDS_IO_V01_PLAN.md`; the temporal-aggregation engine underneath is
`docs/dev_plans/MEDS_IO_DESIGN.md`.

---

## 1. The five stages

```
  raw prognostic state (cohort SoA, patch reservoirs, site scalars)
      │
  [1] DERIVE      meds_diagnostic_kernels     quantities that are a closed-form function of state
      │                                       (LAI, gsc, WUE, soil psi/wetness, CAS VPD, DBH class)
      │
  [2] CAPTURE     meds_core_diag_types        dt-weighted accumulators for everything the fast loop
      │                                       computes per dt_fast and would otherwise discard
      │
  [3] REDUCE      meds_diagnostic_reduce      ONE weighted aggregation:
      │                                       cohort -> {patch, site, PFT, DBH class}; patch -> site
      │
  [4] INTEGRATE   meds_output_integrate       temporal folding per (variable, tier)
      │
  [5] SERIALIZE   meds_output_stream          per-tier, per-time-chunk netCDF
```

Stage [2] is why per-cohort ecophysiology is available at all. Sub-daily resolution exists only
inside the fast loop's sub-step; before it existed, `A_net`, `g_sw`, `C_i`, ψ_leaf, ψ_wood, PLC,
sapflow, root uptake, absorbed radiation and the turbulent fluxes were recomputed roughly 48 times
a day and thrown away.

---

## 2. Aggregation across scales — the rule that matters

Every variable declares three things at registration, beside its units: a **weight kind**, whether
the reduction is a **mean or a sum**, and a unit **scale**. Together these encode the
extensive/intensive distinction, which is a physical statement and the classic place a diagnostic
goes silently wrong.

| | reduction | typical weight | examples |
|---|---|---|---|
| **Extensive** (per plant) | weighted **SUM** → per unit ground area | `nplant` | `agb`, `leaf_area`, `gpp_accum` |
| **Intensive** (a state or a rate per unit leaf) | weighted **MEAN** | leaf area, basal area, `nplant` | `leaf_temp`, `gsw`, `psi_leaf`, `dbh` |

**Which weight is itself a physical statement.** `dbh` is reported basal-area-weighted — the
forestry convention, dominated by the trees that hold the stand — not stem-weighted, which a
regenerating understory would swamp. Canopy temperatures and conductances are leaf-area-weighted,
so a bare sapling does not pull the canopy mean as hard as a closed overstory. Demographic rates
are `nplant`-weighted.

The reduction chain, per patch `p` and cohort `i`, with per-cohort weight `w`:

```math
\text{patch (sum)} = \sum_{i \in p} w_i x_i
\qquad
\text{patch (mean)} = \frac{\sum_{i \in p} w_i x_i}{\sum_{i \in p} w_i}
```

```math
\text{site (sum)} = \sum_p a_p \sum_{i \in p} w_i x_i
\qquad
\text{site (mean)} = \frac{\sum_p a_p \sum_{i \in p} w_i x_i}{\sum_p a_p \sum_{i \in p} w_i}
```

Patch area `a_p` enters **only** at the patch → site step, never inside the weight, so the same `w`
serves every scale. A patch-axis value therefore carries no area factor: it is per m² of *that
patch's own* ground, which is what makes a gap-versus-closed-canopy comparison meaningful.

**Empty sets.** A patch with no cohorts, or a **mean** over a PFT or size class with no members,
emits `_FillValue` — never `0/0`, and never a bare 0 that a reader would take for a measurement.
A **sum** over an empty PFT or size class is a true `0`, because reporting fill there would break
the closure identity below the moment a PFT went locally extinct. The two conventions differ on
purpose.

**Closure identities** (asserted in `test_diagnostic_reduce`, and true to roundoff on real output):

```math
\sum_{\text{pft}} X_\text{pft} \;=\; \sum_{\text{class}} X_\text{size} \;=\; X_\text{site}
```

---

## 3. The axes

| netCDF dim | length | notes |
|---|---|---|
| `time` | UNLIMITED | period **start** stamp; `cell_methods` says how the period was reduced |
| `cohort` | live count, trimmed per file | slot order; `global_cohort_id` tracks a cohort across files |
| `patch` | live count, trimmed per file | ditto `global_patch_id` |
| `soil` | `n_soil_layer_max` | area-weighted site column |
| `pft` | **run-time** PFT count | carries a `pft` coordinate variable, so the file stays self-describing |
| `dbh_class` | from `[output].dbh_class_edges` | carries `dbh_lower` / `dbh_upper` coordinates |
| `(patch, soil)` | 2-D | per-patch soil profiles; `axes_soil_patch = true` |

**DBH classes are a size class of *plants*, ED2-style.** A cohort is assigned whole to the bin
containing its mean `dbh` — never split — and contributes with weight `area · nplant`. Bins are
half-open `[lo, hi)` except the last, which is closed at the top so the largest tree in the stand is
never dropped; a `dbh` below the first edge clamps into bin 1 for the same reason. Both choices
preserve the closure identity, which is what makes `nplant_size` a genuine stem-density
distribution comparable to a forest inventory.

**Cohort and patch axes may not appear on the annual stream.** A window longer than a month would
straddle the annual disturbance restructuring, so the slot set that was averaged would not be the
slot set present at flush. The registry rejects it at start-up.

---

## 4. Aggregation across time

| `agg` | `cell_methods` | meaning |
|---|---|---|
| `AGG_TMEAN` | `time: mean` | dt-weighted mean — the default for a physical state or a mean rate |
| `AGG_MEAN` | `time: mean` | equal-weight mean |
| `AGG_SUM` | `time: sum` | period total (accumulator variables and count tallies) |
| `AGG_LAST` | `time: point` | end-of-period snapshot (ids, CSR, PFT index) |
| `AGG_MIN` / `AGG_MAX` | `time: minimum` / `maximum` | period extremum |

Four tiers — `F` fast, `D` daily, `M` monthly, `Y` annual — each writing its own file family
`<prefix>-<letter>[-<stamp>].nc`. Each tier integrates raw state independently; for these operators
that is identical to chaining.

### A caveat worth stating plainly

**A FAST-tier per-cohort variable is a mean over the fast output window, not an instantaneous
sub-step value.** The per-cohort capture is one dt-weighted accumulator per cohort per slow step,
not a per-sub-step array (which would be `n_cohort × n_sub × n_var`). At
`[output.fast].interval_steps = 1` the two coincide; at the default `4` a FAST record is a
four-sub-step mean. This matters when comparing against a flux tower at sub-hourly resolution.

---

## 5. Choosing what to write

Resolution order — later wins:

1. registry defaults
2. `[output].axes_*` — suppress a whole trailing dimension
3. `[output].<group>` — suppress a whole variable group
4. `[output.<tier>].enabled` — suppress a whole tier
5. `meds_io_config.toml` — per-variable, the finest granularity

**Groups** (8): `structure`, `carbon`, `water`, `energy`, `biogeochem`, `numerics` on by default;
`radiation` and `ecophys` off. `ecophys` is the per-cohort leaf gas-exchange and hydraulics set —
by far the highest-volume group and the one a production run most often wants off.

`numerics` defaults **on** because it carries the energy and water budget residuals. A closure
nobody records is worse than one nobody looks at.

**Axes** are the biggest single lever on output volume: `axes_cohort = false` removes ~55
variables' worth of per-cohort slabs while leaving every site scalar intact.

**Per-variable control**, for debugging:

```bash
meds_main --dump-io-config          # writes meds_io_config.toml: every variable, ready to uncomment
```

```toml
[output]
io_config = "meds_io_config.toml"
```

```toml
[variables]
anet_cohort   = "F D"     # sub-daily leaf physiology for one diagnostic run
agb_size      = "M Y"
growth_avg_cohort = false
```

A name matching no registered variable is a hard error, not a silent no-op.

---

## 6. The minimum set — is this run sane?

Look at these first.

**Fast-scale (14).** `sw_in_site`, `rnet_site` — is the radiation forcing physically shaped?
`le_site`, `h_site` — the Bowen ratio, the single most diagnostic energy number.
`gpp_rate_site` — light-response shape and magnitude. `nee_site` — night `= +Reco`, day `= −uptake`;
the sign flip is the integration test. `et_rate_site` — closes against `le_site`, catches a
latent-heat unit error instantly. `cas_temp_site`, `cas_vpd_site` — canopy-air coupling, the first
thing to oscillate if `dt_fast` is too large. `cas_co2_site` — drawdown amplitude; a stuck value
means the CO₂ twin is not coupled. `soil_temp_site`, `soil_water_site` — the diurnal thermal wave
and the drydown shape. `ustar_site` — is the surface layer turbulent at all?
**`resid_energy_site`, `resid_water_site` — a nonzero value invalidates everything above.**

**Slow-scale (18).** `agb_site`, `lai_site`, `basal_area_site`, `nplant_site` — the four
stand-structure scalars. `agb_pft`, `lai_pft` — is PFT composition doing anything?
`agb_size`, `nplant_size` — the size distribution, the demographic core's actual output.
`npp_site`, `rh_site`, `nee_site` — the carbon balance; drift is the long-run sanity check.
`leaf_resp_site` + `stem_resp_site` + `root_resp_site` — the autotrophic fraction should be roughly
half of GPP, a strong parameter check. `agb_growth_site`, `agb_mort_site`,
`nplant_recruit_site` — the rates that *make* the AGB trajectory. `soilc_total_site` — spinning up,
equilibrating, or running away? `n_cohort_site`, `n_patch_site` — fuse/fission health; a
monotonically climbing count is a config bug. `canopy_height_site` — stand development.

**Budgets you can close from the file alone:**

```math
\frac{\mathrm{d}\,\text{agb}}{\mathrm{d}t} \approx \text{agb\_growth\_site} - \text{agb\_mort\_site}
\qquad
\frac{\mathrm{d}\,\text{soilc}}{\mathrm{d}t} \approx \sum \text{litter}_* - \text{rh\_site}
```

---

## 7. Variable inventory

203 registered variables. Run `meds_main --dump-io-config` for the authoritative list with units,
groups, axes and default streams — it is generated from the registry, so it cannot drift.

| group | count | | axis | count |
|---|---|---|---|---|
| structure | 62 | | site | 98 |
| energy | 38 | | cohort | 55 |
| carbon | 28 | | patch | 23 |
| ecophys | 24 | | pft | 9 |
| water | 18 | | soil | 7 |
| biogeochem | 15 | | dbh_class | 6 |
| numerics | 13 | | (patch, soil) | 5 |
| radiation | 5 | | | |

---

## 8. Adding a variable

Two coupled edits, as designed:

1. a `call add_variable(...)` line in the matching `register_*` routine
   (`meds_output_registry`), naming the field, units, dim, `agg`, group, default streams, and —
   if it is aggregated — the weight kind and mean/sum flag;
2. a `case` in the matching field accessor (`meds_output_integrate`) that copies the value out of
   live state.

Because the reduction is **data** (`dim` + `weight` + `mean` + `scale`) rather than code, a field
that already has an accessor case emits its patch, site, PFT and size-class twins for **one more
registry line each**, with no new extraction code.

For a quantity the fast loop computes and drops, add a row to `cohort_diag_block` or
`patch_diag_block` (`meds_core_diag_types`): a new index parameter, one fill line in the capture,
and its fusion kind. **Zero edits to the lockstep machinery** — the fields are rows of one 2-D
array, so every permutation is a single whole-array statement that cannot omit a field.

### The one trap

The per-cohort blocks are reset per slow step, but restructuring (fuse / split / cull / recruit /
disturb) happens **inside** the slow step, after the fast loop fills them and before the monthly
window closes. Anything read at the output tick must therefore ride the cohort lockstep.

`cohort_deriv_block` (`site%deriv`) does **not** — it is documented as transient and deliberately
unreordered, which is correct for its own consumer, since `update_cohort_states` applies it
immediately. Reading it from the output layer pairs tendency `i` with a different plant `i` on
exactly the month-boundary steps. That mistake was made during this work and was caught only by the
thread-invariance test, because thread count perturbs which cohorts fuse. Use `cohort%sdiag`.
