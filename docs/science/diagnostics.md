# MEDS diagnostic output

How MEDS reports itself: what a variable means, how it is aggregated across the demographic
hierarchy and across time, and how to ask for the ones you want.

Implementation: `src/io/` (`meds_diagnostic_kernels`, `meds_diagnostic_reduce`,
`meds_output_{types,registry,integrate,stream,manager}`) plus the per-cohort / per-patch
accumulators in `src/state/site/meds_site_diag_types.f90`. Design and rationale:
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
  [2] CAPTURE     meds_site_diag_types        dt-weighted accumulators for everything the fast loop
      │                                       computes per dt_fast and would otherwise discard
      │
  [3] REDUCE      meds_diagnostic_reduce      ONE weighted aggregation:
      │                                       cohort -> {patch, site, PFT, DBH class}; patch -> site
      │
  [4] INTEGRATE   meds_output_integrate       temporal folding per (variable, tier)
      │
  [5] SERIALIZE   meds_output_stream          per-tier, per-time-chunk netCDF
```

Stages [1]–[4] run inside the time step; a closed period's record is queued, and stage [5] writes the
queues when a calendar month closes and at the end of the run, so a step makes no netCDF call
(`MEDS_POLYGON_RUNTIME_PLAN.md` §4).

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

**Empty sets.** A patch with no cohorts, or a **mean** with nothing to average, emits
`_FillValue` — never `0/0`, and never a bare 0 that a reader would take for a measurement. "Nothing
to average" covers a PFT or size class with no members, and a whole site, as on bare ground. A
**sum** over an empty PFT or size class is a true `0`, because reporting fill there would break the
closure identity below the moment a PFT went locally extinct. The two conventions differ on purpose.

**A cohort with no samples.** A cohort recruited in the slow step has no fast-loop samples yet at
that step's output tick. It is fill on the cohort axis, and it is left out of every reduction:
- a mean is that of the cohorts that were sampled;
- a sum gains nothing from it.

Before this rule, it read 0, which pulled a canopy temperature towards 0 K.

**What a run does not simulate.** A run without the fast loop (`fast_biophysics_on = false`)
reports every fast-loop variable as `_FillValue`: the fluxes, the canopy air, the ground, the
forcing echo and the per-cohort fast diagnostics. The slow operators' own rates (disturbed area,
recruitment, mortality carbon, and litterfall when soil carbon is on) are weighted by the slow step
instead, and report their values (#299).

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

**Cohort and patch axes may not appear on the annual stream**, and neither may the `(patch, soil)`
profiles. A window longer than a month would straddle the annual disturbance restructuring, so the
slot set that was averaged would not be the slot set present at flush. The registry rejects it at
start-up.

---

## 4. Aggregation across time

| `agg` | `cell_methods` | meaning |
|---|---|---|
| `AGG_TMEAN` | `time: mean` | dt-weighted mean — the default for a physical state or a mean rate |
| `AGG_MEAN` | `time: mean` | equal-weight mean |
| `AGG_SUM` | `time: sum` | period total (accumulator variables and count tallies) |
| `AGG_LAST` | `time: point` | end-of-period snapshot (ids, CSR, PFT index) |
| `AGG_MIN` / `AGG_MAX` | `time: minimum` / `maximum` | period extremum |
| `AGG_FLUXSUM` | `time: sum` | dt-weighted integral of a rate (period total) |

#### Which steps a record holds

A slow step runs from `prev` to `now`. Its fluxes are accumulated over that interval, and its state is
read at `now`, after the step's dynamics. **A step belongs to the period it starts in**, so a record
stamped *d* holds the steps that start in period *d*:

- a daily flux is that day's: the daily record stamped 1 July is the mean of the fast records stamped
  1 July;
- a state is the value at the end of each of the period's steps, so the daily state stamped 31 July is
  the state at 1 August 00:00, and a monthly mean is the mean over the ends of the month's steps;
- no record is stamped at or after the run's end.

**The calendar's restructuring belongs to the period it opens.** At the turn of a month the stand's
cohorts are recruited, fused, split and culled; at the turn of a year its patches are disturbed and
fused. This restructuring runs **between** two slow steps: after the output has read the step that
ends on the boundary, and before the fast loop of the step that begins there. As a result:

- **every record describes one stand.** A monthly cohort or patch record is the mean of its daily
  records, slot by slot, and a file's cohort and patch axes are the same throughout;
- **the restructuring's events and the stand it leaves are recorded in the new period.** The year's
  disturbance area and disturbance mortality are in January's record, never December's, and the
  biomass the disturbance removes first shows in the 1 January daily state. The 31 December state is
  the stand before it.

**A checkpoint on a boundary holds the stand before the restructuring.** The global attribute
`restructure_pending` (`none`, `month` or `year`) records which restructuring is still owed. A run
resumed from the checkpoint performs it first, so its first record holds the same events as the
continuous run's. A state file without the attribute was written after its boundary's
restructuring, and owes none.

#### The skin temperature

`skin_temp_site` is the skin temperature in the land-model sense, as in CLM's `TSKIN`. It also comes
per patch (`skin_temp_patch`), sub-daily (`skin_temp_fast`) and both (`skin_temp_patch_fast`).

- **What it is.** The black-body temperature of the longwave leaving the canopy top, $`(L^\uparrow /
  \sigma)^{1/4}`$, taken per patch and sub-step and then averaged. It is what an infrared thermometer
  or a satellite land-surface temperature sees, and over a forest it is mostly the canopy.
- **What it is not.** It is not the soil. MEDS has no skin layer of its own: snow-free, the ground
  surface is the top soil layer, `soil_temp_top_site`. That layer's node sits 1.8 cm down, so its
  diurnal swing is damped and lagged against the air's.
- **Why it is formed inside the step.** Because of the fourth root, the skin temperature of the mean
  longwave is not the mean skin temperature.
- **The emissivity.** It is taken as 1, so the reflected sky longwave counts as emission. Under a sky
  colder than the surface, that puts the value a few tenths of a kelvin below the surface's own
  temperature, for an emissivity of 0.98.

Four tiers — `F` fast, `D` daily, `M` monthly, `Y` annual — each writing its own file family
`<prefix>-<letter>[-<stamp>].nc`. Each tier integrates raw state independently; for these operators
that is identical to chaining.

### The FAST tier by patch

A site mean over a closed canopy and a gap can describe neither. In the biophysics example stand,
the gap's surface soil ran 13 K above the air at midday while the closed patch's ran 0.7 K below it
(#270). So each FAST quantity of the patch block also has a per-patch twin:

- It is named after its coarse patch variable, plus `_fast`: `cas_temp_patch_fast`,
  `soil_temp_top_patch_fast`, `skin_temp_patch_fast`, `le_patch_fast`, `h_patch_fast`,
  `rnet_patch_fast`, `nee_patch_fast`.
- A quantity with no coarse patch variable keeps its site stem: `gpp_rate_patch_fast`,
  `sw_up_patch_fast`, `lw_up_patch_fast`, `ustar_patch_fast`, `npp_rate_patch_fast`,
  `reco_patch_fast`, `cas_co2_patch_fast`.
- The soil columns by layer and patch are `soil_temp_layer_patch_fast` and
  `soil_water_layer_patch_fast`.

The site `*_fast` value is the area-weighted sum of the patch values, taken in patch order, so the
two always agree. The twins follow `axes_patch` (on by default) and `axes_soil_patch` (off), as the
coarse patch variables do. Like every patch variable, they are not written for a region.

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

**Groups** (9): `structure`, `carbon`, `water`, `energy`, `biogeochem`, `numerics` and `forcing` on
by default; `radiation` and `ecophys` off. `forcing` is the atmospheric boundary the run used, after
the reader's shortwave partition, rain/snow split and optional corrections: site means at the daily,
monthly and yearly tiers (`air_temp_site`, `qair_site`, `psurf_site`, `wind_site`, `lwdown_site`,
the four shortwave streams `par_beam_site` … `nir_diffuse_site`, `snowfall_site`, `atm_co2_site`,
`cosz_site`, `rho_air_site`) and their sub-daily `*_fast` twins, and three patch rows at the daily and
monthly tiers: the forcing each patch saw at its canopy-air top, `wind_cas_top_patch` and
`air_temp_cas_top_patch`, and that top's height `cas_depth_patch` (with `rough_patch` and
`displace_patch`, in `energy`, the move can be rebuilt; [`forcing.md`](forcing.md) §8). With `sw_in_*`
and `precip_site`, it is what checks the sub-daily reconstruction against a tower. `ecophys` is the per-cohort leaf gas-exchange and hydraulics set —
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

252 registered variables. Run `meds_main --dump-io-config` for the authoritative list with units,
groups, axes and default streams — it is generated from the registry, so it cannot drift.

| group | count | | axis | count |
|---|---|---|---|---|
| structure | 67 | | site | 141 |
| energy | 48 | | cohort | 55 |
| carbon | 34 | | patch | 29 |
| forcing | 28 | | pft | 9 |
| ecophys | 24 | | soil | 7 |
| water | 18 | | dbh_class | 6 |
| biogeochem | 15 | | (patch, soil) | 5 |
| numerics | 13 | | | |
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
`patch_diag_block` (`meds_site_diag_types`): a new index parameter, one fill line in the capture,
and its fusion kind. The patch block's row is also the FAST tier's sample (`patch_diag_row`), so a
patch row reaches the sub-daily tier with one more registry line (`SRC_F_PD0 + PD_*`), per patch
or as a site mean. **Zero edits to the lockstep machinery** — the fields are rows of one 2-D
array, so every permutation is a single whole-array statement that cannot omit a field.

### The one trap

The per-cohort blocks run from one output tick to the next, and the stand is reordered while they
fill: `sort_cohorts` re-sorts it every slow step, after the fast loop has filled them. Restructuring
(fuse / split / cull / recruit / disturb) also runs between two steps. Anything read at the output
tick must therefore ride the cohort lockstep.

`cohort_deriv_block` (`site%deriv`) does **not** — it is documented as transient and deliberately
unreordered, which is correct for its own consumer, since `update_cohort_states` applies it
immediately. Reading it from the output layer pairs tendency `i` with a different plant `i` whenever
the stand has been reordered since it was formed. That mistake was made during this work and was caught only by the
thread-invariance test, because thread count perturbs which cohorts fuse. Use `cohort%sdiag`.
