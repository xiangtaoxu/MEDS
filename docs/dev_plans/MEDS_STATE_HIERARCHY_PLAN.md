# The state layer by model hierarchy — plan

**Revision 4, 2026-10-18. The owner agreed the design on 2026-10-18, all decisions are settled (§7),
and the owner asked for the implementation the same day.**

**Against:** `beta` at `b0d47a7` (after #447, which dissolved `meds_column_params`). Every line count,
file list and importer count below was read from that commit. Beta has not moved since.

**Asked for:** fold `meds_column_state_types` into the site state, then split the site state into
files by model hierarchy (cohort, patch, site). The outputs must not change: every step is a move, a
rename or a reorder, and each is checked bit for bit against `b0d47a7`.

**What changed since Revision 3** (the owner's review of the site and diagnostics files, 2026-10-18):
- **Site sections.** `site_t` is laid out in sections 0–6 like the cohort and patch blocks, and its
  routines run lifecycle, structure, identity, per step (§3.4).
- **Diagnostic rows in sections.** The rows of each level are regrouped by meaning, and the patch's
  slow rows come after its fast rows, so `patch_diag_slow_row` becomes a range check (§3.5).
- **Three follow-ups:** the unused `site_area`, the 11 output-only accumulators on `site_t`, and one
  way to store slow diagnostic rows (§9).
- **Decisions.** D14 and D15 are new (§7).

**What changed since Revision 2** (the owner's review of the patch file, 2026-10-18):
- **Patch sections.** The patch file is laid out like the cohort file: the store types and
  `patch_block` in one top-down order (canopy air, snow, soil water and heat, soil carbon), the helper
  fields grouped by lifecycle, and every list of patch fields in that order (§3.2, §6.5).
- **One rename.** `soil_column_t` becomes `soil_water_column_t`, to pair with `soil_energy_column_t`
  (§3.2).
- **Decisions.** D12 and D13 are new (§7).

**What changed since Revision 1** (the owner's review, 2026-10-18):
- **Names.** The files end in `_state_vars`, after ED2's `ed_state_vars.F90`. The diagnostics module
  becomes `meds_diag_vars` (§3.1).
- **Cohort sections.** The cohort file is laid out in eleven sections, and every list of cohort fields
  follows that order (§3.3, §6.5).
- **PFT table.** The 16 per-cohort `p_*` copies of PFT values go. Instead the cohort block carries one
  PFT table, as it already carries `allom` (§6.6).
- **Pool changes.** `carbon_flux_block` leaves the state layer. It moves to `meds_vegetation_dynamics`,
  its only user, and is renamed `pool_change_block`, because it holds each pool's net change rather
  than NPP (§6.7).
- **NPP by pool.** This is already a diagnostic (§3.5). Per-cohort and per-PFT outputs of it, and
  folding the pool changes into the tendency bundle, are follow-ups (§9).
- **Decisions.** D1–D7 are settled as recommended, and D8–D11 are new (§7).
- **Steps.** There are now eight steps instead of five (§8).

---

## 1. Summary

The state layer has two folders, each its own library:

| Now | Lines | Holds |
|---|---|---|
| `src/state/column/meds_column_state_types.f90` | 362 | one patch's stores: canopy air, soil water, soil heat, snow, soil carbon, the daily decomposition accumulator, the litter input; their fusion blends; `cas_set_depth` |
| `src/state/site/meds_site_state_types.f90` | 1,412 | `cohort_block`, `patch_block`, `site_t`, and every routine on them |
| `src/state/site/meds_site_diag_types.f90` | 630 | the diagnostic accumulators, per cohort, per patch and per polygon |

After the plan, `src/state/` is one folder and one library, with one file per level of the hierarchy
and one for the diagnostics:

| After | Lines (est.) | Holds |
|---|---|---|
| `src/state/meds_cohort_state_vars.f90` | ~930 | `cohort_block` in eleven sections, the tendency bundle, and everything that lists or sets the cohort fields |
| `src/state/meds_patch_state_vars.f90` | ~580 | one patch's stores and their blends (the whole of today's column module), `patch_block`, and everything that lists its fields, including the one patch reorder |
| `src/state/meds_site_state_vars.f90` | ~180 | `site_t`, which holds one cohort block and one patch block, and the routines that need both |
| `src/state/meds_diag_vars.f90` | 630 | today's `meds_site_diag_types.f90`, renamed, contents unchanged |

The cohort file is shorter than Revision 1's estimate (~1,070) because the 16 `p_*` fields leave every
field list.

There is no polygon file (§3.6). The rule "kernels never see `site_t`" stays. A new lint test enforces
it on every ctest run, replacing the library boundary, which only looked like it enforced it (§5).

---

## 2. Why the column folder can go now

The column/site split was decision #2 of `MEDS_CODE_STRUCTURE_DESIGN.md` (PR #125). It had two
reasons. #447 removed the first, and the second never worked as described.

1. **The build order.** `meds_config` read `n_soil_layer_max` and the soil parameter bundles from
   `state/column`, and the site state reads `meds_config`. One state library would therefore have
   been a cycle.
   - Since #447, config reads nothing from state: the order is shared → {config, state_column} →
     state_site.
   - One library linking shared and config is now acyclic.
2. **The kernel guard.** The kernel libraries link `meds_state_column`, not `meds_state_site`. This was
   meant to make a kernel that names `site_t` fail to build.
   - Tested on 2026-10-17 by adding `use meds_site_state_types, only : site_t` to a soil kernel.
     Building `meds_fast_kernels` alone from a clean tree fails, as intended. But a clean full build
     passes, and so does the kernel-only test `test_water_processes`.
   - The reason: every `.mod` file goes to one shared directory (`CMakeLists.txt:39-40`). Once
     `meds_site_state_types.mod` exists, any file can use it, and importing only a type needs nothing
     at link time.
   - No automated step builds a kernel library on its own: there are no CI workflows, and ctest
     doesn't. So today the rule is held by review.

So folding the column module costs no protection the build actually gives today. The guard moves to a
lint that runs on every ctest (§5), which is stronger than what it replaces.

---

## 3. Target layout

### 3.1 Names

- **`meds_cohort_state_vars`, `meds_patch_state_vars`, `meds_site_state_vars`.**
  - The names point to ED2's `memory/ed_state_vars.F90`. That module defines `patchtype`, `sitetype`
    and `polygontype` with their allocate and copy routines, which is the same job these files do, and
    `docs/ed2_comparison.md` already maps MEDS's state file to it.
  - "State variables" is also what an ecologist calls a cohort's `nplant` and carbon pools.
- **One difference, stated in each file's header.** ED2's module also holds the model's state itself,
  as the module-level global `edgrid_g` (line 3672). MEDS has no global state. These modules define the
  types and the routines on them, and the state is passed as arguments.
- **`meds_diag_vars`** for the diagnostics module.
  - It drops "site" because it holds blocks at three levels: `cohort_diag_block`, `patch_diag_block`
    and `polygon_diag_block`.
  - `meds_state_diag` was considered and not used. In MEDS's I/O vocabulary "state" means the restart
    stream (`[state]`, the `-S-` files), and these accumulators are exactly what never goes there.
- **Unchanged:**
  - the field-index names (`CD_*`, `CS_*`, `PD_*`, `PY_*`);
  - the kernel and work `_types` modules (`meds_plant_types`, `meds_soil_types`, `meds_fast_types`
    and the rest). They hold parameters, fluxes and work arrays, not the model's state.

So `src/state/` holds the state variables in three files and the diagnostic variables in one.

### 3.2 `meds_patch_state_vars.f90` — what one patch holds, and the block of all patches

**What moves in:**

| Item | From |
|---|---|
| `cas_state_t`, `soil_column_t` (renamed `soil_water_column_t`), `soil_energy_column_t`, `snow_column_t`, `soil_carbon_t`, `xi_accum_t`, `litter_input_t` | column module |
| `blend_cas`, `blend_soil_w`, `blend_soil_e`, `blend_snow`, `blend_soil_carbon`, `blend_litter_input`, `blend_xi_accum`, `cas_set_depth`, `tiny_mass` | column module |
| `patch_block` | site module, lines 255–322 |
| `patch_alloc`, `patch_ensure_capacity` | site module |
| `patch_free` (new name for the patch half of `site_free`'s deallocate list) | site module, `site_free` |
| `patch_birth_state` | site module, with a smaller argument list (§6.2) |
| `patch_reorder` | new, replacing the two patch field lists in demography (§6.3) |

The store types are part of the patch's state, so they belong in the patch file. Today they sit in a
separate folder only to keep them away from `site_t`.

**Sections.** The patch file follows the cohort file's rules (§3.3), plus one of its own:
- **The physical stores run top-down:** canopy air, snow, soil water and heat, soil carbon. That is the
  vertical order of the patch's column, and the order of the science pages
  (`canopy_air_space_biophysics`, `snow_biophysics`, `soil_biophysics`, `soil_carbon`).
  - Today there are two orders. The column module defines soil water, soil heat, snow, canopy air,
    litter, soil carbon, accumulator. `patch_block`, the sort and the fusion list canopy air, soil
    heat, soil water, snow.
  - The one order runs through the type definitions, the blends, `patch_block`, and every routine that
    lists all the patch fields: `patch_alloc`, `patch_ensure_capacity`, `patch_reorder`, `patch_free`,
    `fuse_2_patches`, and the new-gap code in `apply_patch_disturbance`.
- **The helper fields are grouped by lifecycle,** as in the cohort file. There are three kinds, and
  they behave differently on a restart:
  - daily rates passed from the slow loop to the fast loop (set once per slow step, held all day,
    restarted, 0 in a new gap);
  - accumulators within a step (not restarted, reset each step);
  - solver memory (numerical, not physical; a new gap starts it cold).
- **The cohort map is kept apart from identity.** `rebuild_csr` rebuilds it after every reorder, and it
  is never restarted.
- **Derived values stay in the store types, marked as derived.** The fast kernels read them together
  with the prognostic values, and the blends already handle them. `can_depth` is the patch's
  equivalent of cohort size: `refresh_canopy_depth` derives it from the tallest cohort at setup and at
  each slow step, and it is not in the restart file.

**How every field maps.** Kind: **P** = prognostic, **D** = derived, **M** = memory, **A** =
accumulator, **B** = bookkeeping.

| § | Field | Unit | Kind | In restart file | On fusion | New gap |
|---|---|---|---|---|---|---|
| **0. Block record** | `n`, `cap` | — | B | — | — | — |
| **1. Identity and area** | `global_id` | — | B | yes | receiver's | new id |
| | `dist_type` | — | B | yes | only patches of the same type fuse | treefall |
| | `area` | fraction of site | P (disturbance, fusion, termination) | yes | summed | its share of the disturbed area |
| | `age` | yr | M (advanced by `update_patch_states`) | yes | area-weighted | 0 |
| **2. Cohort map** | `cohort_offset`, `cohort_count` | — | D (`rebuild_csr`) | no | rebuilt | rebuilt |
| **3. Recruit bank** | `recruit_pool(pft, patch)` | plant m⁻² | P | yes | area-weighted | the donors' mean |
| **4. Canopy air** (`cas`) | `can_enthalpy`, `can_shv`, `can_co2` | J kg⁻¹, kg kg⁻¹, µmol mol⁻¹ | P | yes | weighted by air mass | the donors' mean |
| | `can_temp` | K | D (from enthalpy) | yes | re-diagnosed | re-diagnosed |
| | `can_depth` | m | D (tallest cohort) | no; refreshed at setup | area-weighted | the donors' mean |
| **5. Snow** (`snow`) | `swe`, `snow_energy` | kg m⁻², J m⁻² | P | yes | area-weighted | the donors' mean |
| | `snow_depth`, `snow_temp`, `snow_fliq`, `nlayer` | m, K, –, – | D | yes | re-diagnosed | re-diagnosed |
| **6. Soil water and heat** | `soil_w%theta(layer)` | m³ m⁻³ | P | yes | area-weighted | the donors' mean |
| | `soil_w%w_surface`, `w_surface_enth` (pond) | kg m⁻², J m⁻² | P | yes | area-weighted | the donors' mean |
| | `soil_e%soil_energy(layer)` | J m⁻³ | P | yes | area-weighted | the donors' mean |
| | `soil_e%soil_temp`, `soil_fliq` | K, – | D | yes | re-diagnosed | re-diagnosed |
| **7. Soil carbon** | `soil_carbon`: 7 pools and 2 lignin | kgC m⁻² | P | yes | area-weighted | the donors' mean |
| **8. Daily slow-to-fast rates** | `shed_water_rate` | kg m⁻² s⁻¹ | A (set once per slow step) | yes | area-weighted | 0 |
| | `slow_co2_rate` | µmol m⁻² s⁻¹ | A (set once per slow step) | yes | area-weighted | 0 |
| **9. Within-step accumulators** | `xi_accum`: 7 decomposition integrals and `rh_fast_accum` | day, kgC m⁻² | A (fast → slow) | no | area-weighted | the donors' mean |
| | `litter_in`: 6 pool streams and 6 by tissue | kgC m⁻² day⁻¹ | A (vegetation → soil) | no | area-weighted | the donors' mean |
| **10. Solver memory** | `adapt_dt_last` | s | M | yes | area-weighted | 0 (cold start) |
| **11. Output accumulators** | `diag` (from `meds_diag_vars`) | — | A | no | by fusion kind | — |

That is 17 per-patch arrays, the same as today: nothing is added or removed.

**The file from top to bottom:**
1. The store types in section order: `cas_state_t`, `snow_column_t`, `soil_water_column_t`,
   `soil_energy_column_t`, `soil_carbon_t`, then the accumulators `xi_accum_t` and `litter_input_t`.
2. `patch_block`, in sections 0–11.
3. The routines:
   - the blends, in the same order, with `cas_set_depth` after `blend_cas`;
   - `patch_alloc`, `patch_ensure_capacity`, `patch_reorder`, `patch_free`;
   - `patch_birth_state`.

**One rename: `soil_column_t` becomes `soil_water_column_t`.** It holds only water (layer moisture and
the ponded surface water), and the new name pairs it with `soil_energy_column_t`. The type name
appears 12 times in 4 source files (the two state modules, `meds_fast_types`,
`meds_biogeochem_types`) and in two science pages (`soil_biophysics.md`, `column_biophysics.md`). The
field names `soil_w` and `soil_e` stay: they appear 242 times.

### 3.3 `meds_cohort_state_vars.f90` — the cohort block and everything that lists its fields

**What moves in:**

| Item | From the site module |
|---|---|
| `GROWTH_AVG_UNSET`, `DMAX_PSI_LEAF_UNSET`, `DMAX_PSI_LEAF_ACCUM_RESET`, `PHENO_FLUSH_INIT`, `PHENO_SHED_INIT`; `N_COHORT_BOUNDED`, `COHORT_BOUNDED_NAME` | lines 49–68, 411–412 |
| `cohort_block`, `cohort_deriv_block` | lines 70–250, 324–335 |
| `cohort_deriv_alloc`, `cohort_alloc`, `cohort_ensure_capacity`, `move_alloc_block`, `cohort_reorder`, `cohort_compact`, `copy_cohort_slot` | the cohort lockstep |
| `cohort_free` (new name for the cohort half of `site_free`'s list) | `site_free` |
| `scale_cohort_ground_fields`, `fuse_cohort_fast_state`, `set_cohort_wood_geometry`, `set_cohort_size`, `set_cohort_size_from_carbon` | |
| `cohort_tissue_heat_capacity`, `cohort_tissue_water`, `cohort_wood_store_carbon`, `set_cohort_tissue_water`, `cohort_bad_field`, `cohort_field_name`, `init_cohort` | |

Two items from today's file do not move in:
- `gather_pft_params` is deleted (§6.6);
- `carbon_flux_block` leaves the state layer (§6.7).

This is where "add a per-cohort field" happens. It is the one file that lists every cohort field, in
each of its lockstep routines.

**Sections.** The type is laid out in sections. Every routine that lists all the cohort fields uses the
same order and the same section comments:
- `cohort_alloc`, both when allocating and when setting the initial values;
- `cohort_free`;
- `cohort_ensure_capacity`;
- `move_alloc_block`;
- `cohort_reorder`;
- `copy_cohort_slot`.

Today each list has its own order, so checking one against the type means hunting.

The sections follow three rules:
- **Prognostic state is kept apart from what is derived from it.** `nplant` and the four carbon pools
  are the only state the model integrates. Size and geometry are recomputed from the pools by one
  routine, `set_cohort_size_from_carbon`.
- **`cohort_block` holds no rates.** The tendencies are a separate, temporary type (below). What the
  block holds for demography is memory: the growth window that the mortality law reads.
- **Helper fields are grouped by lifecycle, not by currency.**
  - The phenology fields persist and are written to the restart file, so they are state, in a section
    of their own.
  - The fast-to-slow sums are added to on every fast sub-step, read by the slow step, then reset. The
    GPP and respiration sums sit there, not with the carbon pools.

**How every field maps.** Kind: **P** = prognostic (integrated), **D** = derived (recomputed),
**M** = memory (persists, not integrated), **A** = accumulator (reset each slow step), **B** =
bookkeeping.

| § | Field | Unit | Kind | In restart file | On fusion |
|---|---|---|---|---|---|
| **0. Block record** (one per block, not per cohort) | `n`, `cap` | — | B | — | — |
| | `allom` | — | parameters | no | — |
| | `pft_table` (new, §6.6) | — | parameters | no | — |
| **1. Identity and density** | `pft` | — | B | yes | same PFT only |
| | `global_id` | — | B | yes | survivor's |
| | `owner_patch` | — | B | yes | same patch only |
| | `nplant` | plant m⁻² | P | yes | summed |
| **2. Carbon pools** | `leaf_carbon`, `fineroot_carbon`, `wood_carbon`, `nonstructural_carbon` | kgC plant⁻¹ | P | yes | nplant-weighted, so totals are conserved |
| **3. Size, derived from carbon** | `dbh`, `height`, `basal_area`, `agb`, `leaf_area` | cm, m, cm², kgC, m² (per plant) | D | written for readers, recomputed on read | recomputed |
| | `wood_area`, `sapwood_carbon`, `sapwood_area` | m², kgC, cm² (per plant) | D | no | recomputed |
| | `overtopping_lai` (canopy position) | m² m⁻² | D, recomputed from the stand | yes | recomputed |
| **4. Water** | `leaf_water_mass`, `wood_water_mass` | kg plant⁻¹ | P | yes | nplant-weighted |
| | `leaf_surf_water`, `wood_surf_water` (intercepted water films) | kg m⁻² ground | P | yes | summed (per ground area) |
| | `dmax_psi_leaf` (yesterday's predawn) | MPa | M | yes | survivor's |
| **5. Energy** | `leaf_temp`, `wood_temp` | K | P | yes | weighted by heat capacity |
| **6. Plastic leaf traits** | `sla`, `vcmax25`, `rd25`, `leaf_turnover_rate` | m² kgC⁻¹, µmol m⁻² s⁻¹, yr⁻¹ | P (they acclimate) | yes | weighted by what each trait is per (SLA on leaf carbon, vcmax25 on leaf area) |
| **7. Leaf phenology** | `leaf_flush_tendency`, `leaf_shed_tendency` | 0–1 | M | yes | survivor's |
| | `growing_degree_days`, `cold_degree_days` | K day | M | yes | survivor's |
| | `dry_psi_sum`, `wet_psi_sum` | MPa day | M | yes | survivor's |
| | `light_hours_mean` | h day⁻¹ | M | yes | survivor's |
| **8. Growth history** (read by the mortality law) | `growth_avg`, `growth_accum`, `growth_count`, `growth_hist(window, cohort)` | cm yr⁻¹ | M | yes | survivor's |
| **9. Fast-to-slow accumulators** | `gpp_accum`, `leaf_resp_accum`, `stem_resp_accum`, `root_resp_accum` | kgC plant⁻¹ | A | no | nplant-weighted |
| | `light_hours_accum` | h | A | no | survivor's |
| | `dmax_psi_leaf_accum` (running daily maximum, reset at the daily rollover) | MPa | A | yes | survivor's |
| **10. Output accumulators** | `diag`, `sdiag` (blocks from `meds_diag_vars`) | — | A | no | by each field's fusion kind |
| **Removed** | the 16 `p_*` fields (§6.6) | — | copies of PFT values | no | — |

That leaves 45 per-cohort arrays plus the two diagnostic blocks, down from 61.

**Tendencies section** (after `cohort_block`, in the same file):
- **What it holds:** `cohort_deriv_block` and its allocator `cohort_deriv_alloc`. The type holds
  `d_leaf_carbon_dt`, `d_fineroot_carbon_dt`, `d_wood_carbon_dt`, `d_nonstructural_carbon_dt` and
  `dln_nplant_dt`: the demographic rates.
- **Why it stays a type of its own** rather than fields of `cohort_block`:
  - **Lifecycle.** The rates are valid only within one step, between being computed and applied, before
    any sort or fusion.
    - Inside `cohort_block`, they would either join all eight field lists for nothing, or sit in the
      block without following the cohorts.
    - Once the step re-sorts the cohorts, index `i` in the rates would no longer be the same plant.
  - **It is the applier's input.** The rates are `intent(in)` to `update_cohort_states`, so the
    compiler guarantees the applier can't change them (rule 4). The Python demography API fills the
    bundle with rates computed outside the model and calls the same applier.
  - **ED2 shows the cost of the alternative.** `ddbh_dt`, `dagb_dt` and `mort_rate` are fields of its
    cohort type (`ed_state_vars.F90` lines 164, 176 and 433). So every cohort copy and pack carries
    them (lines 11133 and 11878).
- **Why it stays in the state layer:** `site_t` carries it, so it is allocated once and reused every
  step, and the C API fills it.

### 3.4 `meds_site_state_vars.f90` — the site: one cohort block and one patch block

| Item | Note |
|---|---|
| `site_t` | holds `cohort_block`, `patch_block`, the polygon diagnostic block, the tendency bundle and the site scalars |
| `site_alloc`, `site_free` | `site_alloc` also sets the cohort block's PFT table (§6.6); `site_free` becomes two calls, `cohort_free` and `patch_free` (§6.1) |
| `rebuild_csr` | it needs both blocks: it counts cohorts per patch and reorders cohorts |
| `assign_cohort_id`, `assign_patch_id` | the id counters are site scalars |
| `reset_step_diagnostics` | resets all three diagnostic blocks |

**Sections of `site_t`.** The same order as the cohort and patch blocks: what it holds, bookkeeping,
memory, accumulators, scratch, output.

| § | Field | Unit | Kind | In restart file | Written by → read by |
|---|---|---|---|---|---|
| **0. The blocks it holds** | `cohort`, `patch` | — | — | their own fields | — |
| **1. Dimensions and id counters** | `n_pft` | — | B | yes | `site_alloc` → recruit pool, output, demography |
| | `next_cohort_id`, `next_patch_id` | — | B | yes | `assign_cohort_id`, `assign_patch_id` |
| | `site_area` | — | unused | yes | the restart only (§9) |
| **2. Shared cohort clock** | `growth_hist_pos` | — | M | yes | slow driver → the slot every cohort's `growth_hist` writes next |
| **3. Climate memory** | `t_growth_avg` | K | M | yes | `meds_polygon`, daily → thermal acclimation |
| **4. Fast-to-slow accumulators** (the model reads them) | `pheno_tair_sum`, `pheno_tair_n` | K, count | A | no | fast loop → phenology |
| **5. Tendencies** (scratch) | `deriv` | yr⁻¹ | — | no | slow driver, C API → `update_cohort_states` |
| **6. Output accumulators** (only output reads them) | `diag` (the polygon block) | — | A | no | fast loop → output |
| | `et_accum`, `runoff_surf_accum`, `drainage_accum` | kg m⁻² | A | no | fast loop → output |
| | 8 `work_*` counters | count | A | no | fast loop → output |

**Routines, in this order:** lifecycle (`site_alloc`, `site_free`); structure (`rebuild_csr`); identity
(`assign_cohort_id`, `assign_patch_id`); per step (`reset_step_diagnostics`).

### 3.5 `meds_diag_vars.f90` — the diagnostic accumulators

The diagnostics file is renamed. Its skeleton stays: the header with its six design decisions, the
fusion kinds, per level the rows and their fusion tables, the three block types, then the routines by
level in the same order for each (alloc, free, grow, reset, reorder, copy, clear, fuse or blend,
value).

**The rows are regrouped into sections (D15).** Today each level's rows are numbered in the order they
were added, so the patch's fast and slow rows alternate. Regrouping renumbers the rows. That is safe and
bit-identical: output, the tests and the C API read rows by name, nothing reads a row by number, and
the diagnostics are never restarted. Each fusion table is reordered with its rows.

- **Cohort, fast rows (`CD_*`, 30):**
  1. carbon: `ANET`, `AGROSS`, `RD`, `GPP_RATE`;
  2. stomata and leaf exchange: `GSW`, `GBW`, `CI`, `CS`, `TRANSP`, `BETA_STOM`, `BETA_NONSTOM`;
  3. plant water: `PSI_LEAF`, `PSI_WOOD`, `PLC`, `SAPFLOW`, `ROOT_UPTAKE`, `LEAF_WATER`, `WOOD_WATER`;
  4. energy and radiation: `LEAF_TEMP`, `WOOD_TEMP`, `ABS_PAR`, `ABS_SW`, `ABS_LW`;
  5. canopy environment: `LEAF_VPD`, `WIND`;
  6. the leaf solve's other inputs, kept together for the calibration: `LEAF_PAR`, `GB_MOL`,
     `PSI_PREDAWN`, `CA`, `PRESSURE`.
- **Cohort, slow rows (`CS_*`, 11):** the order stays; the rows get section comments (growth and
  mortality; carbon allocation; carbon losses).
- **Patch (`PD_*`, 42):** the rows the fast loop writes, top-down, then the rows the slow loop writes.
  - Fast rows (32): forcing at the canopy-air top (`WIND_CAS_TOP`, `TAIR_CAS_TOP`, `Z_CAS_TOP`);
    top-of-canopy radiation (`SW_IN_VIS`, `SW_IN_NIR`, `SW_UP_VIS`, `SW_UP_NIR`, `SW_UP`, `LW_UP`,
    `SKIN_TEMP`); turbulence (`USTAR`, `GGNET`, `ROUGH`, `DISPLACE`); surface energy balance (`RNET`,
    `H`, `LE`, `SW_GROUND`, `LW_GROUND`); carbon and water fluxes (`GPP`, `NPP`, `RECO`, `NEE`,
    `TRANSP`); canopy air and ground state (`CAS_TEMP`, `CAS_SHV`, `CAS_CO2`, `CAS_VPD`,
    `SOIL_TEMP_TOP`, `W_SURFACE`); budget residuals (`RESID_ENERGY`, `RESID_WATER`).
  - Slow rows (10): litter (`LITTER_LEAF`, `LITTER_FINEROOT`, `LITTER_WOOD`, `LITTER_REPRO`,
    `LITTER_STRUCT`); demography (`RECRUIT_NPLANT`, `DISTURB_AREA`, `MORT_C_BACKGROUND`,
    `MORT_C_CULL`, `MORT_C_DISTURB`).
  - With the slow rows last, `patch_diag_slow_row` loses its hand-kept list of the 10 slow rows and
    becomes a range check against a new `N_PDIAG_FAST`.
- **Polygon (`PY_*`, 16), the forcing as the model used it:** radiation (`SW_IN`, `PAR_BEAM`,
  `PAR_DIFFUSE`, `NIR_BEAM`, `NIR_DIFFUSE`, `LWDOWN`, `COSZ`); precipitation (`PRECIP`, `RAINF`,
  `SNOWFALL`); air (`TAIR`, `QAIR`, `PSURF`, `WIND`, `CO2`, `RHO_AIR`).

**It already holds NPP by pool.**
- **Rows.** The slow cohort block `sdiag` has the rows `CS_NPP_LEAF`, `CS_NPP_FINEROOT`, `CS_NPP_WOOD`,
  `CS_NPP_STORAGE` and `CS_NPP_REPRO`, together with `CS_GROWTH_RESP`, `CS_STORAGE_RESP` and
  `CS_EXUDATE` (lines 148–155), in kgC plant⁻¹ yr⁻¹.
- **Writer.** `compute_carbon_allocation` writes them (`meds_vegetation_dynamics.f90:709-719`).
- **Output.**
  - As site totals: `npp_leaf_site` through `root_exudate_site` (`meds_output_registry.f90:601-616`).
  - Per cohort: only through the C API (`meds_c_api_run.f90:281-286`).
- **What they record** is gross allocation, so they can't stand in for the model's own pool changes
  (§6.7).

### 3.6 No polygon file

`meds_polygon_t` (`src/main/meds_polygon.f90`) is not state. It holds the run's config copy, the
forcing cursor, the fast-loop context, the budgets and the slow ledger, so it needs the forcing, fast
and slow libraries, all of which sit above state. A polygon holds exactly one `site_t`, so the
polygon's ecological state already is `site_t` with its site scalars (`t_growth_avg`, the water and
work accumulators). `meds_polygon_t` stays orchestration, in `main/`.

### 3.7 Dependencies inside the library

```
meds_diag_vars           (shared only)
   ↑
meds_patch_state_vars    (diag; therm_lib; config + soil-column config for patch_birth_state)
meds_cohort_state_vars   (diag; allometry, water retention, numerics; config for init_cohort and pft_table_t)
   ↑
meds_site_state_vars     (cohort + patch + diag)
```

The cohort and patch files do not use each other. No module in the library uses a kernel module.

---

## 4. Library and build

One library replaces two:

```cmake
file(GLOB STATE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/state/*.f90)
add_library(meds_state STATIC ${STATE_SOURCES})
target_link_libraries(meds_state PUBLIC meds_shared meds_config)
```

| Target | Now links | After |
|---|---|---|
| `meds_fast_kernels`, `meds_slow_kernels` | `meds_state_column`, `meds_config` | `meds_state`, `meds_config` |
| `meds_demography` | `meds_state_site`, `meds_slow_kernels` | `meds_state`, `meds_slow_kernels` |
| `meds_fast` | … `meds_state_site` … | … `meds_state` … |
| `test_soil_column_config` | `meds_config`, `meds_state_column` | `meds_config` (it no longer needs state) |

The graph in `src/README.md` becomes:

```
  shared ──→ config ──→ state ──→ demography ──→ io_prep
                          │
                          ├──→ fast_kernels   (the kernels take only the per-patch store types:
                          └──→ slow_kernels    the kernels_no_site_state lint, §5)
```

That is 21 CMake libraries instead of 22. Each kernel library still builds on its own:
`cmake --build … --target meds_fast_kernels` from a clean tree builds `meds_shared`, `meds_config` and
`meds_state` first.

---

## 5. The kernel guard: a lint test

**Rule.** Code in the kernel folders (`fast_dynamics/{canopy,plant,soil}`, `slow_dynamics/{plant,soil}`):
- takes from the state layer only the per-patch store types and their blends, through
  `use meds_patch_state_vars, only : …`;
- never names `site_t`, `cohort_block` or `patch_block`;
- never uses `meds_site_state_vars`, `meds_cohort_state_vars` or `meds_diag_vars`.

A routine that needs more than that is driver code, as now.

**Check.** A new script, `scripts/lint/check_kernel_state.py`, run as ctest `kernels_no_site_state`
next to `clamps_named` and `no_issue_refs`. It:
- reads every `.f90` in the kernel folders, with comments stripped;
- fails on any of the three type names or the three module names;
- fails on a `use meds_patch_state_vars` whose only-list names anything outside the allowed list, or
  has no only-list at all.

Allowed list: `cas_state_t`, `soil_water_column_t`, `soil_energy_column_t`, `snow_column_t`,
`soil_carbon_t`, `xi_accum_t`, `litter_input_t`, the seven `blend_*`, and `cas_set_depth`.

Today the five kernel files that touch state already pass it. They import `snow_column_t`,
`soil_energy_column_t`, `soil_carbon_t` and `litter_input_t` (twice). The only mentions of the
forbidden names are two comments.

**Why a lint rather than the library boundary.**
- It runs on every ctest, in every build, and it can't be beaten by build order (§2).
- It is finer than any library split: it lets a kernel take `soil_carbon_t` but not `patch_block` from
  the same module.
- The alternative of giving each library its own `.mod` directory would make every library edge real.
  It still can't separate two types in one module, so it doesn't replace the lint here. It stays a
  separate option, not part of this plan (§9).

The implementation PR shows the lint failing once on a scratch copy with an injected `site_t` import,
and says so.

---

## 6. Changes that ride with the move (all bit-identical)

### 6.1 Each level frees its own fields

`site_free` deallocates every cohort and every patch array in one statement. It becomes:

```fortran
subroutine site_free(site)
   call cohort_free(site%cohort)
   call patch_free(site%patch)
end subroutine site_free
```

`cohort_free` and `patch_free` hold today's two halves of the list, unchanged. The point is that the
cohort file and the patch file each list their own fields, and the site file lists none.

### 6.2 `patch_birth_state` takes the patch, not the site

```fortran
! now
pure subroutine patch_birth_state(site, ip, cfg, soil, temp, psi)
! after
pure subroutine patch_birth_state(cas, soil_w, cfg, soil, temp, psi)
   type(cas_state_t),         intent(in) :: cas
   type(soil_water_column_t), intent(in) :: soil_w
```

- It reads only `site%patch%cas(ip)` and `site%patch%soil_w(ip)`, so its home is the patch file.
- Its callers pass the two elements:
  - `apply_recruitment`;
  - `add_cohort`;
  - the restart fallback in `meds_io`;
  - the two birth tests.

### 6.3 One reorder for the patch arrays

Every patch field is listed twice in demography, besides the lists in the state file:
- `sort_patches` permutes 15 arrays plus the diagnostics;
- `patch_compact` packs the same 15 with `pack`.

Because of this, the rules file has to state an exception: "patch arrays have no single reorder
routine". The plan adds `patch_reorder(patch, perm, m)` to the patch file, the twin of
`cohort_reorder`:
- `sort_patches` becomes four steps: compute the age order, `patch_reorder`, remap `owner_patch`,
  `rebuild_csr`;
- `patch_compact` builds the list of surviving indices, as `cohort_compact` does, and calls
  `patch_reorder`.

Moving data by index gives the same values as `pack`, so the result is bit-identical. The patch fields
are then listed in one file, as the cohort fields are, and the exception leaves the rules.

### 6.4 Comments that explain the old layout

Several comments justify a type's place by the library boundary that goes. Examples:
- "lives here (not in `meds_biogeochem_types`) so the site state's `patch_block` can carry it with no
  state → soil-kernel library edge";
- the CMake comments on the two state libraries;
- `meds_biogeochem_types.f90:47`.

They are rewritten to say why the type is in the patch file: it is what a patch holds, and the kernels
take it from there under the lint.

### 6.5 All four files in section order

This step reorders the declarations and the statements of each field list into the order of §3.3 (the
cohort file), §3.2 (the patch file), §3.4 (the site file) and §3.5 (the diagnostic rows), and adds the
section comments. For the patch, the order also
reaches the two field lists outside the file, `fuse_2_patches` and the new-gap code in
`apply_patch_disturbance`. The step also renames `soil_column_t` to `soil_water_column_t`. No
statement changes otherwise.
- **Bit-identical:** each field is moved or blended on its own, so statement order changes no values.
  The restart writer and reader in `meds_io` keep their own grouping by array shape and find variables
  by name, so this step doesn't touch them.
- **Its own commit,** so that the move in step 1 stays verbatim and easy to diff.

### 6.6 One PFT table on the cohort block

**What goes.** The 16 per-cohort copies of PFT values:
- `p_dbh_critical`, `p_wood_density`, `p_hgt_max`;
- `p_aboveground_frac`, `p_root_to_leaf_ratio`, `p_storage_cushion`;
- `p_wai_b1`, `p_wai_b2`, `p_sapwood_area_b1`, `p_sapwood_area_b2`;
- `p_leaf_width`, `p_branch_diameter`, `p_is_woody`;
- `p_stem_resp_factor25`, `p_root_resp_factor25`, `p_stem_resp_sapwood25`.

How they are handled today:
- They are copied from `cfg%pft` at birth (`init_cohort`) and at restart (`gather_pft_params`,
  `meds_io.f90:719`).
- Nothing else ever writes them, and they are not in the restart file.
- Each one is hand-written into every field list.

**Why they existed, and why that no longer holds.** The code comment says "kernels never index the PFT
table", meaning hot loops should not look a value up by PFT number. That doesn't match the code today:
- The hot fast-loop kernels don't read `cohort_block`. They read the column view, which is rebuilt once
  per patch per day.
- The same kernels already look up PFT tables by index:
  - the leaf photosynthesis table (`meds_leaf_gas_exchange.f90:547`);
  - the hydraulics table (`meds_fast_prepass.f90:243`, `meds_fast_step.f90:84`);
  - the phenology PAR threshold (`meds_fast_dynamics.f90:750`).

So the copied set is an arbitrary subset, and two styles live side by side.

**What replaces them.**
- `cohort_block` gets `type(pft_table_t) :: pft_table` in section 0, beside `allom`.
- `site_alloc` sets it from the config the site is built from, as it already sets `allom`. That means
  one new argument, `pft`.
- `column_cohort_init`, which builds a block outside `site_alloc`, sets it too.
- The table is never restarted (the restart path goes through `site_alloc`, as for `allom`), never
  reordered, and never written after it is set.
- The name is `pft_table`, not `traits`, so it can't be confused with section 6's plastic leaf traits,
  which belong to each cohort.

**How the readers change.**
- **Most reads.** `cohort%p_wood_density(i)` becomes `cohort%pft_table%wood_density(cohort%pft(i))`.
  - This covers the size routines, `copy_column_cohort` in `meds_column_view` (8 reads) and
    `meds_output_integrate` (1 read).
  - No signature changes: the 17 calls of `set_cohort_size_from_carbon` stay as they are.
- **Where a routine already has the config's table in hand,** it reads that instead, so no routine
  reads the same table two ways:
  - `meds_vegetation_dynamics` (7 reads), which already reads `pft%root_to_leaf_ratio(pf)` next to
    `cohort%p_wood_density(j)`;
  - `meds_c_api_demography` (5 reads).

**Other changes.**
- `gather_pft_params` is deleted, and so is its call in `meds_io`.
- `site_alloc` gains its argument in `meds_init` and `meds_io`, and in the three tests that call it:
  `test_output_integrate`, `test_output_roundtrip` and `test_diagnostic_reduce`.
- `test_carbon_growth` sets five `p_*` fields by hand; it sets the table instead.

**Results** are bit-identical: the values are the same and come from the same table.

**Memory.** `pft_table_t` holds about 124 per-PFT arrays, so the copy is 5–10 KB per site with 5–10
PFTs: a few MB for a 983-cell region. That is against 16 arrays per cohort removed.

### 6.7 `carbon_flux_block` becomes `pool_change_block`, in `meds_vegetation_dynamics`

**What it holds.** Each pool's net change over the step, in kgC per plant per step:

| Pool | `npp%…` in the block | The `sdiag` NPP row |
|---|---|---|
| leaf | growth − shed | growth |
| fine root | growth − shed | growth |
| wood | growth | growth |
| storage | allocation − storage maintenance + carbon resorbed from shed leaves and roots | allocation |

So it isn't NPP. A deciduous tree dropping its leaves has a gross leaf allocation of 0 and
`npp%leaf < 0`. The tissue-water shedding relies on that sign (`meds_vegetation_dynamics.f90:826`,
`:839`).

**Who uses it.** It is one local variable, `npp`, in `vegetation_dynamics`:
- `compute_carbon_allocation` fills it;
- `compute_vital_rates`, `shed_turnover_water` and `update_cohort_derivatives` read it;
- outside the module, only `test_slow_ledger` and `test_carbon_growth` use it.

Rule 3 says a type lives with whoever changes it.

**The change.**
- The type moves to `meds_vegetation_dynamics`, public for the two tests, renamed `pool_change_block`.
- The local and the dummy arguments named `npp` become `dpool`.
- `docs/science/plant_carbon_allocation.md` follows.
- Results are bit-identical.

**Why the NPP diagnostics can't replace it.**
- They record gross allocation, a different quantity (the table above).
- They exist only when output asks for them: the writes sit inside `if (cohort%sdiag%active)`. The
  model must never read something whose existence depends on output settings.

**What would remove it** is folding it into the tendency bundle. That is a follow-up (§9).

---

## 7. Decisions (settled 2026-10-18)

| | Question | Settled |
|---|---|---|
| **D1** | Folder layout | flat `src/state/` with four files |
| **D2** | File names | `meds_cohort_state_vars`, `meds_patch_state_vars`, `meds_site_state_vars`, after ED2's `ed_state_vars`; the diagnostics become `meds_diag_vars` (§3.1) |
| **D3** | The kernel guard | the lint of §5; per-library `.mod` directories stay a possible later change |
| **D4** | The diagnostic file | kept whole, because its three blocks share one design (a 2-D field array, the fusion kinds, alloc/grow/reorder); renamed only |
| **D5** | One reorder for the patch arrays | in this PR (§6.3) |
| **D6** | `patch_birth_state`'s arguments | element arguments, in the patch file (§6.2) |
| **D7** | The cohort file's hand-written field lists | a follow-up plan, after this lands (§9) |
| **D8** | Layout of the cohort file | sections 0–10 and a tendencies section; every field list in section order (§3.3, §6.5) |
| **D9** | The 16 `p_*` fields | one PFT table on the cohort block (option A, §6.6). Option B, passing `cfg%pft` to about 25 call sites, was not chosen |
| **D10** | `carbon_flux_block` | moves to `meds_vegetation_dynamics` as `pool_change_block` (§6.7); folding it into the tendency bundle is a follow-up |
| **D11** | NPP by pool | already in `meds_diag_vars` (§3.5); per-cohort and per-PFT outputs are a follow-up |
| **D12** | Layout of the patch file | sections 0–11, the physical stores top-down, the helpers by lifecycle; every patch field list in that order (§3.2, §6.5) |
| **D13** | `soil_column_t` | renamed `soil_water_column_t`; the field names `soil_w` and `soil_e` stay (§3.2) |
| **D14** | Layout of the site file | `site_t` in sections 0–6; routines by lifecycle, structure, identity, per step (§3.4) |
| **D15** | The diagnostic rows | regrouped into sections per level, the patch's slow rows last (§3.5) |

---

## 8. Steps and verification

One branch, one PR, one commit per step. Steps 1–7 are bit-identical. The work happens in the worktree
`state-hierarchy-plan` (branch `refactor/state-hierarchy`).

| Step | What | Files |
|---|---|---|
| **1. Move and rename** | Create the four files from verbatim blocks of the three modules, under the new module names. Delete `src/state/column/` and `src/state/site/`. One `meds_state` library. Rewrite every `use` with a script: the names in each statement go to the module they now live in. | 4 state files; `CMakeLists.txt`; **56 importers: 32 in `src/`, 24 tests** |
| **2. Each level owns its routines** | `cohort_free`/`patch_free` (§6.1); `patch_birth_state` in the patch file with element arguments (§6.2) | state files; `meds_demography_cohort_fusefiss`, `meds_init`, `meds_io`, 2 tests |
| **3. Sections** | the types and every field list in section order; `soil_column_t` → `soil_water_column_t`; the diagnostic rows regrouped (§6.5). Three commits: cohort, patch, site and diagnostics | the four state files; `meds_demography_patch_fusefiss`; `meds_fast_types`, `meds_biogeochem_types` |
| **4. PFT table** | `pft_table` on the block; the `p_*` fields and `gather_pft_params` go (§6.6) | cohort and site files; `meds_init`, `meds_io`, `meds_column_view`, `meds_vegetation_dynamics`, `meds_c_api_demography`, `meds_output_integrate`; 4 tests |
| **5. Pool changes** | `pool_change_block` in `meds_vegetation_dynamics` (§6.7) | cohort file; `meds_vegetation_dynamics`; 2 tests; `plant_carbon_allocation.md` |
| **6. The guard** | `scripts/lint/check_kernel_state.py`; ctest `kernels_no_site_state` | 1 script; `CMakeLists.txt` |
| **7. One patch reorder** | `patch_reorder`; `sort_patches` and `patch_compact` use it (§6.3) | patch file; `meds_demography_patch_fusefiss` |
| **8. Docs** | see §10 | |

**Importer churn by level** (from a scan of `b0d47a7`). All three module names change, so every
importer changes. In Revision 1, the 20 files that take only site-level names changed nothing; now
they change one word.

| They take | Files | Change |
|---|---|---|
| only site-level names (`site_t` and its routines) | 20 (7 in `src/`, 13 tests) | the module name |
| only patch store types | 13 (the 5 kernel files, 4 fast-loop numerics, 4 tests) | the module name |
| cohort names, alone or with one other level | 9 | the module name, and a second `use` where they take two levels |
| patch and site | 2 | a second `use` |
| all three levels | 8 (the demography operators, `meds_init`, `meds_io`, `meds_polygon`, `meds_fast_dynamics`, one test) | split into two or three `use` lines |
| only the diagnostics | 4 | the module name |

Of the 52 files that import state, 16 also import the diagnostics; that `use` changes name too. The
rewrite is mechanical, by the same kind of script #447 used for the column-params imports, keeping
each file's own alignment.

**Verification.**
- **Bit for bit:** the two-year census comparison of #447 (reference config, fast loop on, 51 netCDF
  files, 8,107 variables) matches `b0d47a7` exactly after each of steps 1–5 and 7.
- **The restart path:** `test_state_roundtrip` and `test_output_roundtrip` cover it. That is where
  `gather_pft_params` used to fill the copies, and where step 4 relies on `site_alloc` instead.
- **Test suites:**
  - ctest on ifx Release after every step;
  - at the end, ifx Debug, gfortran 11.5, and ifx with `MEDS_BUILD_PYLIB` plus `python/tests`.
- **Standalone builds:** `meds_state`, `meds_fast_kernels` and `meds_slow_kernels` each build alone
  from a clean tree.
- **The lint:** it passes on the tree, and fails once on a scratch copy with an injected `site_t`
  import.
- **Not possible here:** nvfortran is not installed. Steps 1–3 and 5–7 don't change any type's
  contents, but step 4 adds a component with allocatable arrays to `cohort_block`, so the PR will say
  the nvfortran check was not run.

---

## 9. Out of scope, and follow-ups

- **The cohort file's field lists (D7).** Every cohort field is written out by hand eight times:
  - the type;
  - `cohort_alloc` twice (allocate, then initial values);
  - the free list;
  - `cohort_ensure_capacity`'s copy;
  - `move_alloc_block`;
  - `cohort_reorder`;
  - `copy_cohort_slot`.

  Then come the creation, fusion and restart sites. After this plan, there are 45 fields instead of
  61, and every list is in the same order, which makes the next step easier to review. Possible
  consolidations:
  - one index-copy routine covering `cohort_reorder` and `copy_cohort_slot`;
  - a per-field grow helper replacing the capacity copy and `move_alloc_block`;
  - an `intent(out)` free.

  Together these would bring eight lists down to about four. Each touches compiler traps listed in
  `CLAUDE.md` (allocatable components, array temporaries), so it deserves its own plan and
  measurements.
- **Fold the pool changes into the tendency bundle (from D10).** Today the numbers make a round trip:
  1. the allocator produces amounts;
  2. `advance_pool` adds them to get tentative pools;
  3. `fill_cohort_deriv` converts back to rates, (new − old)/dt;
  4. `update_cohort_states` multiplies by dt again.

  If the allocator wrote the four `site%deriv%d_*_carbon_dt` rates directly:
  - `pool_change_block` would disappear;
  - the vital rates and the tissue-water shedding would read the rates;
  - `update_cohort_derivatives` would keep only the geometry update from the new pools, the growth
    history and the mortality.

  The C-API path, which takes new pools from Python, still goes through `fill_cohort_deriv` unchanged.
  This is not bit-identical: (x/dt)·dt can differ from x in the last bit, and demography can amplify
  that. So it needs its own issue and PR, compared within a tolerance.
- **Per-cohort and per-PFT NPP by pool (from D11).** These are registry entries on the existing `sdiag`
  rows, `npp_{leaf,fineroot,wood,storage,repro}_{cohort,pft}`, at one `add_variable` each. They add
  outputs, so they are not part of this PR.
- **"storage" or "nonstructural".** The diagnostics and outputs say storage (`CS_NPP_STORAGE`,
  `npp_storage_site`), and so do the PFT parameters (`storage_cushion`, `storage_turnover_rate`). The
  pool is `nonstructural_carbon`. One word could be chosen in a later naming pass; renaming outputs
  changes users' scripts.
- **Per-library `.mod` directories (D3).** These would make every library edge real. They are a
  possible later change and don't replace the lint.
- **`site%site_area` is never used.** Only the restart writes and reads it. With one site per polygon
  and patch areas as fractions that sum to 1, it is always 1. Deleting it changes the restart file's
  metadata layout, so it needs its own change that keeps old restart files readable.
- **11 output-only accumulators sit on `site_t`.** `et_accum`, `runoff_surf_accum`,
  `drainage_accum` and the 8 `work_*` counters are read only by output, so by this plan's rule they
  belong in `meds_diag_vars`, as rows of the polygon block. Today they are always accumulated (not only
  when output asks), reset every slow step by the fast loop rather than every output tick, and summed
  rather than weighted by time step. Moving them changes the output arithmetic, so it is its own change.
- **One way to store slow diagnostic rows.** The cohort keeps its slow rows in a second block, set once
  per step; the patch writes its slow rows into the fast block as a time-weighted contribution, which
  needs the `fast_rows` flag and `patch_diag_slow_row` for slow-only runs. A separate slow block per
  level, as the cohort has, needs no special case. It changes output arithmetic, so it is its own
  change.
- **Two placeholder defaults in `cas_state_t`.** `can_co2 = 400` and `can_depth = 20 m` are never
  seen by a run with the fast loop (the canopy air is seeded from the forcing, and its depth set from
  the stand, at setup). The type's other unseeded values are 0, and `cohort_birth_state` reads
  `can_temp > 0` as "seeded". Making the two 0 as well would change test fixtures and runs without the
  fast loop, so it is its own change.
- **The TOML site file** (the location, ground and starting state split from the model settings): the
  owner will revisit.
- **"column" elsewhere.** The fast loop's `column_cohort_t`, `meds_column_view`, `build_column_frozen`
  and the `column_*` tests keep the word. There it means one patch seen vertically, which is what the
  fast loop integrates; only the state folder of that name goes.

---

## 10. Documents to update in step 8

| File | What |
|---|---|
| `src/README.md` | the state tree and the four file names; rule 1 (kernels never see `site_t`, now checked by the lint); rule 3 (boundary state lives in the patch file); the library graph; the "where does a new file go" rows for a reservoir, a per-cohort field (its section, and every list in section order) and a per-patch field; the lockstep invariant (lines 98, 113, 131, 159–160, 178) |
| `CLAUDE.md` | rules 1 and 3 (lines 82, 88); the lockstep invariant |
| `.claude/rules/state-demography.md` | where the cohort lockstep lives (line 25), the section order, the PFT table on the block; the patch exception goes with D5 |
| `.claude/rules/fast-loop.md` | line 49: "canopy/, plant/ and soil/ link state_column and config only" |
| `CMakeLists.txt` | the state-layer comment block |
| `src/fast_dynamics/README.md` | line 24 (the stores live in `state/column`) |
| `docs/science/` code tables | `canopy_air_space_biophysics.md:101`, `snow_biophysics.md:89`, `soil_biophysics.md:400`, `column_biophysics.md:104` and `:150`, `soil_carbon.md:304`, `order_of_processes.md:244`, `plant_traits.md:105` and `:152`, `diagnostics.md:8` and its other two mentions of the diagnostics module; `soil_column_t` in `soil_biophysics.md` and `column_biophysics.md` |
| `docs/science/plant_carbon_allocation.md` | `carbon_flux_block` becomes `pool_change_block` (with step 5) |
| `docs/ed2_comparison.md:439` | the state files' paths, now named after `ed_state_vars` |
| `CHANGELOG.md` | one entry; no user-visible effect |

Archived plans keep the old paths, as the dev_plans rule says.

---

## 11. Implementation (2026-10-18)

Done on branch `refactor/state-hierarchy`, one commit per step. The order differs from §8 in one
way: the PFT table (step 4), the pool-change type (step 5) and the patch reorder (step 7) went in
before the section reorder of their files, so no list was written in section order only to be
rewritten in the next commit.

| Commit | Step |
|---|---|
| Plan the state layer by model hierarchy | this document |
| State: one library of four files by model hierarchy | 1 |
| State: each level frees its own fields; patch_birth_state takes the patch | 2 |
| Cohort block: one PFT table instead of 16 per-cohort copies | 4 |
| Move the pool-change type to its one user, named for what it holds | 5 |
| Cohort file in sections | 3, cohort |
| One reorder for the patch arrays | 7 |
| Patch file in sections; soil_column_t becomes soil_water_column_t | 3, patch |
| Site file and diagnostic rows in sections | 3, site and diagnostics |
| Lint: the kernels take only the per-patch store types | 6 |
| Docs | 8 |

**Verification.** After every step: ctest on ifx Release, and the two-year census comparison of #447
against `b0d47a7`, identical on all 51 files and 8,107 variables. After the diagnostic renumbering,
also a one-year run with every output group and the FAST tier on: identical on 26 files and 4,822
variables. The lint passes on the tree and fails on a copy with an injected `site_t` import and on
an injected `patch_block` in a kernel's only-list. ctest is 76 tests with the lint.

**Found on the way.** `cohort_ensure_capacity` did not carry `overtopping_lai` across a grow (the
only one of the 45 cohort arrays it left out), so a grow reset it to 0. The block grows in cohort
fission (`split_cohorts`), recruitment (`apply_recruitment`), the survivors copied into a new gap
(`apply_patch_disturbance`) and setup (`add_cohort`). The reset was never seen: each of those paths
ends with `update_overtopping_lai` (the end of `restructure_stand`, and the driver after setup),
nothing between reads the field, and the readers (trait plasticity, the output, the state file) all
come after. It is carried now, so a grow keeps every field, and the result is bit-identical.

**Constants, after review (2026-10-18).** The "unset" markers (`GROWTH_AVG_UNSET`,
`DMAX_PSI_LEAF_UNSET`, `DMAX_PSI_LEAF_ACCUM_RESET`) stay in the cohort file: they define what their
fields can hold, and `meds_constants` keeps only physical and calendar constants. The two newborn
phenology values are private to the cohort file, which alone uses them. The patch file's
`tiny_mass` duplicated `tiny_num` (both 1e-30) and is gone, and the snow store's default temperature
is written as `t_3ple`. Bit-identical.

**One rename after review (2026-10-18).** `patch_birth_state` is now `cohort_birth_state`: it reads a
patch's canopy air and top soil water and returns the temperature and water potential a new cohort
is born at, so the name says whose birth it describes. It stays in the patch file, since it reads
only patch data. The name `init_patch` is kept for a patch constructor the code does not have yet:
setup, the census reader and the new gap in `apply_patch_disturbance` each set a new patch's fields
themselves. Bit-identical.

**Cohort-file review (2026-10-18).** The header explains the CSR map and the comparison with ED2 moved
to `docs/ed2_comparison.md` §4.6 (from all three state headers). Field annotations are one line each,
with the detail in the science pages (the predawn double buffer: `leaf_gas_exchange.md` §4.2). The
"no value yet" markers, a newborn's phenology and the canopy-air placeholders are in a fifth state
file, `meds_state_defaults`. The block routines are `expand_cohort_block_size`,
`move_cohort_block`, `reorder_cohort_block` and `compact_cohort_block`, and
`set_cohort_wood_geometry` is `set_cohort_wood_allometry`. A tissue's heat capacity is
`tissue_heat_capacity` in `meds_therm_lib`; `cohort_tissue_heat_capacity` is its per-slot wrapper.
Bit-identical.
