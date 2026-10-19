# MEDS tissue energy plan: `leaf_energy` and `wood_energy`

> # 🗃️ ARCHIVED — 2026-10-19. Complete: every step shipped in PR #457.
>
> **What shipped**, one commit per step, with the departures recorded in §8:
> - one definition of the tissues' dry heat capacity (E1);
> - the canopy films per plant (E2);
> - the whole wood store in the wood's heat capacity (E3);
> - `leaf_energy` and `wood_energy` [J/plant] as each cohort's tissue state, with the temperatures derived from them (E4);
> - the shed water's enthalpy, `shed_water_enth` (E5).
>
> **The live description** is `docs/science/vegetation_energy_dynamics.md` (eq. (8) and "The tissue's energy is its state") and `docs/science/soil_plant_water.md` §2 and §5.
>
> **Deferred** to `docs/ROADMAP.md` §5: tissue freezing and cold-driven mortality (#456), and the films inside the tissue energy (#165).

**Revision 3, 2026-10-19.** Approved by the owner; implemented in PR #457 (log in §8).

- **Base:** `origin/beta` at f6dfa9e (after S5 of the multi-site plan).
- **Origin:** the last open item of the `src/state` review. Each cohort should carry its tissue energy, with the temperature diagnosed from it.

**Revision 2 records the owner's decisions of 2026-10-19:**

| Decision | Owner's choice |
|---|---|
| D1 | Energy per plant. The canopy films (`leaf_surf_water`, `wood_surf_water`) become per plant too: at cohort level a field is per plant, and it is scaled up with nplant or LAI where needed. |
| D3 | Temperatures stay stored, derived from the energy. |
| D4 | The whole wood store counts in the wood's heat capacity. The refactor should also remove the shed-water special cases. |
| D7 | The state file holds both the energy and the temperatures. |
| D8 | Tissue water stays supercooled liquid, because live plants keep their sap from freezing. The assumption, and a request for cold-driven mortality, is issue #456. |

The decisions the owner did not comment on (D2, D5, D6, D9, D10) are taken as recommended. Two new ones are listed in §6: D11 (§2.11) and D12 (§2.10).

## 0. Summary

- **New state.** Each cohort carries `leaf_energy` and `wood_energy` [J/plant] as its prognostic tissue state.
- **Temperatures become derived.** `leaf_temp` and `wood_temp` are derived from the energy, as `soil_temp` is derived from `soil_energy`.
- **One energy datum.** The energy uses the soil's convention, `temp_to_internal_energy`: dry heat capacity × T, plus the water's liquid internal energy. Every other store uses this datum, and so does every water flow that reaches the tissue.
- **The fast kernel is unchanged.** The step still relaxes each tissue exactly, at the heat capacity held over the step. Only the commit changes: it adds energy, then diagnoses the temperature.
- **The slow tier keeps its rule.** Matter a slow operator adds or removes does so at the tissue's own temperature. The energy follows at that held temperature, at the ledger brackets that already exist around each operator group.
- **Fusion sums energy,** as it sums water.
- **Wood heat capacity counts the whole wood-store water.** This includes the fine-root share, which the water flows already move at the wood's temperature.
- **Canopy films become per plant,** like the tissue water. The per-ground rescaling of films at splits, patch fusion and disturbance goes.
- **Shed water carries its own enthalpy to the ground.** Today it is valued as the precipitation it joins, or at the canopy-air temperature on a dry day.
- **§3 lists the special cases this removes.**
- **Steps:**

  | Step | Content | Numbers |
  |---|---|---|
  | E1 | one heat capacity | bit-identical |
  | E2 | per-plant films | round-off, plus the dead's film shed with them |
  | E3 | whole wood store | small physics change |
  | E4 | energy state | small numerical change |
  | E5 | shed enthalpy | small physics change |
  | E6 | docs | — |

## 1. Why

**Tissue is the one store held as a temperature.** The canopy air carries enthalpy (`can_enthalpy`). The soil, pond and snow carry internal energy (`soil_energy`, `w_surface_energy`, `snow_energy`). Leaf and wood carry only `leaf_temp` and `wood_temp`.

**The tissue store and the water that crosses it use two datums.**

- The tissue's heat is valued as cap·T:
  - `tissue_store0/1` in `meds_fast_ark.f90:429-469` and `meds_fast_rk45.f90:620-660`;
  - `slow_site_store` and `slow_tissue_heat` in `meds_slow_ledger.f90:367-396`.
- Water entering or leaving the tissue is valued with u_liq(T) = c_liq·(T − T_sc):
  - `qwflux_wl` and `q_wood_net` (soil_plant_water.md eq. (4));
  - `water_store_enth` (eq. (5)).
- **Size of the gap.** T_sc is `tsupercool_liq` = 56.8 K, so the two datums differ by c_liq·T_sc ≈ 2.4×10⁵ J per kg of tissue water.

**What follows from the two datums.**

- Each fast step closes its own energy budget.
- The store a step ends with, cap(W₀)·T₁ + ΔW·u(T₀), is not the store the next step starts with, cap(W₁)·T₁.
- The difference per step is ΔW·c_liq·(T₁ − T₀ + T_sc).
- No ledger compares one step's end with the next step's start, so nothing fails. But the tissue's energy is not a quantity the model can add up across steps.
- **Example:** at the BCI storm front, wet soil refills dry-season wood with about 7 kg m⁻² in a quarter hour (soil_plant_water.md §5). That one step re-values the store by about 1.7×10⁶ J m⁻².
- **The temperatures are right to second order.** This is a bookkeeping inconsistency, not a physics error.

**The wood's heat capacity and its water flows disagree on what the store is.**

- Eq. (8) counts only the sapwood share of the wood-store water (`wood_thermal_water`, i.e. W·B_sap/(B_sap + B_root)).
- Eqs. (4)–(5) move and book the whole store, fine roots included, at the wood's temperature.
- So the fine-root water's enthalpy is booked to the tissue, but its heat capacity is not.
- This needs special cases:
  - the fine-root share in `shed_turnover_water` (`meds_vegetation_dynamics.f90:848-862`);
  - fusion's non-additive remainder, declared as a thermal-mass change (`meds_cohort_state_vars.f90:812-817`);
  - the wood tolerance in `test_fusion_cohort` (line 109).

**The films are the one per-ground field on the cohort block.**

- `leaf_surf_water` and `wood_surf_water` are in kg per m² of ground, while every other cohort amount is per plant.
- So every operator that changes the ground a cohort's density refers to must rescale them (`scale_cohort_ground_fields`):
  - the split (`meds_demography_cohort_fusefiss.f90:303-304`);
  - patch fusion (`meds_demography_patch_fusefiss.f90:259, 266`);
  - disturbance survivors (`:578`).
- Two film-conservation error-stops guard that rescaling (split `:321-322`, patch fusion `:271-272`).
- Mortality needs its own rule: the dead's film stays on the survivors, because the film does not scale with nplant (`meds_vegetation_dynamics.f90:964-966`).

**The shed water has no enthalpy of its own.**

- Every slow path that hands water to the ground adds it to `patch%shed_water_rate`: turnover, excess water, mortality, the cull, and disturbance.
- In the fast tier that water joins the ground inflow (`meds_fast_frozen.f90:329-331`). It is valued at that inflow's temperature (`:343-349`):
  - the rain's, on a rainy day;
  - the melt's, under a pack;
  - the canopy air's, on a dry day.
- The slow ledger declares the water leaving at the tissue's temperature. The difference is a re-valuation at the handoff.

**There are two definitions of the heat capacity.**

| Where | Definition |
|---|---|
| `build_column_frozen` (`meds_fast_frozen.f90:216-246`) | `veg_thermal%c_leaf` and `c_wood_dry`, on the column view |
| `cohort_tissue_heat_capacity` (`meds_cohort_state_vars.f90:670`) | the constants |

They agree today only because the `veg_thermal` fields are not config keys and default to the constants.

**Freezing would need an energy state.** On the melt plateau the temperature is pinned while energy goes into the phase change. Under D8 the tissue does not freeze (§7, #456), but the energy state keeps that door open.

## 2. Design

### 2.1 The quantities

Per plant, for each tissue:

| Tissue | Dry heat capacity | Water |
|---|---|---|
| Leaf | D = 2·c_leaf·C_leaf | W = `leaf_water_mass` |
| Wood | D = 2·c_wood·C_wood | W = `wood_water_mass` (the whole store, §2.4) |

```text
E   = D·T + W·c_liq·(T − T_sc)              (= temp_to_internal_energy(D, W, T, 1))
T   = (E + W·c_liq·T_sc) / (D + W·c_liq)    (the liquid branch of internal_energy_to_temp)
cap = D + W·c_liq                           (dE/dT at fixed water: unchanged)
```

- **Units.** E is in J per plant, like the tissue water [kg/plant], the carbon pools [kgC/plant] and, from E2, the films [kg/plant]. So the energy:
  - fuses as an nplant-weighted amount;
  - needs no rescaling when a patch's area changes;
  - is copied whole when a cohort splits.
- **The fast tier** works per m² of ground (nplant × the per-plant value), as it does now.
- **Magnitudes.** A large tree's wood holds about 10⁹ J per plant, which is no problem in double precision.
- **A leafless leaf** has D = W = 0, so E = 0 and its temperature is that of its wood. This is the kernel's existing "no area" rule, and it becomes the diagnosis rule too. Wood always has carbon.
- **Below 0 °C** the liquid branch keeps the tissue water supercooled (D8, #456).

### 2.2 Where the state lives

The cohort file's section 5 (Energy) holds:

| Field | Units | Kind |
|---|---|---|
| `leaf_energy`, `wood_energy` | J plant⁻¹ | prognostic (restart, fusion, lockstep) |
| `leaf_temp`, `wood_temp` | K | derived: always exactly T(E, mass); saved as twins (D7) |

**Two helpers in the cohort file** (section "TISSUE WATER AND HEAT"):

- `set_cohort_tissue_temps(cohort, i)` sets the temperatures from the energy and the current mass. It is used after the fast commit and after fusion.
- `hold_cohort_tissue_temps(cohort, i)` re-derives the energy at the held temperature and the new mass, then re-derives the temperature from that energy. This means "the matter came or went at the tissue's own temperature". It is used by the slow brackets (§2.5) and at birth.

**Two kernels in `meds_therm_lib`:**

- `tissue_energy(D, W, T)`;
- `tissue_temp(E, W, D)`, with the leafless case handled by the caller.

**Why the temperature stays stored.**

- The slow operators that change mass need the temperature that matter comes or goes at.
- Inside one operator group, cohorts are reordered, merged, created and culled, so that temperature must travel with the cohort.
- The soil does the same: it stores `soil_temp` beside `soil_energy`.

**Invariant.** The stored temperature is always exactly `tissue_temp(E, mass)`: both helpers end by deriving T from E. So a resumed run matches a continuous one.

**Guard.** At the slow-to-fast seam, the fast loop derives T from E anyway. The guard error-stops if the derived T differs from the stored T by more than 10⁻⁶ K, which would mean some path changed a cohort's mass without holding its temperature. It costs two compares per cohort per day; that is the bounds-guard pattern of #392.

### 2.3 The fast tier

**At the seam** (copy-in, `meds_fast_dynamics.f90:594-609`):
- The working bundle gets the energy, the films and the tissue water per m² of ground (× nplant).
- It also gets the start-of-step temperatures derived from the energy.

**In `build_column_frozen`:**
- The capacities come from the one function (E1).
- Eq. (5) is kept per cohort, as `leaf_water_enth(i)` and `wood_water_enth(i)` [W m⁻²]. Today it is only the column sum, `water_store_enth`.
- Nothing else changes. The upwind temperatures of eqs. (4)–(5) are still the start-of-step ones.

**The kernel** (`veg_energy_balance`, `surface_derivs`) is unchanged.

**The commit.** One routine, shared by ARK and RK45, replaces the two copies at `ark:429-469` and `rk45:620-660`:

```text
E_end = E_0 + cap_0·(T* − T_0) + H_i·dt      T* = the march's time-average temperature, as now
T_end = tissue_temp(E_end, W_end, D)
store0 = Σ E_0,  store1 = Σ E_end            the same difference as today: the whole-column check is unchanged
```

- With the `veg_energy` mask off: T* = T₀, so E_end = E₀ + H·dt and T_end = T₀ exactly. Water that arrives at the tissue's temperature leaves the temperature unchanged.
- **What changes numerically:**
  - T_end = (cap₀·T* + ΔW·c_liq·T₀)/cap₁ instead of T*. The water a tissue gained is valued at T₀, as eq. (5) already books it. This is a second-order change.
  - The next step starts from E_end, so the cross-step re-valuation of §1 is gone.

**At the seam** (copy-out, `:813-821`): the per-plant values are the bundle's divided by nplant.

### 2.4 The wood's thermal mass: the whole wood store (D4)

```text
cap_wood = c_wood·B_wood + c_liq·n·W_wood
```

- This agrees with eqs. (4)–(5), which already move and book all of the store's water at T_wood.
- It removes `wood_thermal_water`, and makes the capacity additive across fusion and splitting.
- **Rough size of the change,** from example05's last state file (1780, 561 cohorts). Sapwood carbon there is estimated, because the state file does not carry it. The fine-root share of the wood-store water is:

  | Cohorts | Fine-root share (median) | Stand wood capacity rises by |
  |---|---|---|
  | all | 2 % | about 1 % |
  | dbh < 5 cm | 7 % | about 4 % |

- **What stays:** the fine roots' share of the water a fine-root loss sheds, W·B_root/(B_sap + B_root) in `shed_turnover_water`. That is hydraulics, not heat: the store's capacity is in proportion to its carbon. It is written directly from `cohort_wood_store_carbon`.

### 2.5 The slow tier

**Rule (unchanged):** matter that a slow operator adds or removes does so at its tissue's temperature. Growth makes tissue from CO₂ and water at the plant's temperature, and turnover sheds it at that temperature. The model tracks no heat for CO₂, so that heat crosses the system boundary as a declared thermal-mass exchange. This is the same mechanism as today.

**The brackets already exist.** `slow_tissue_heat(poly)` is taken before and after each operator group that changes mass:

| Bracket | Operators | Where |
|---|---|---|
| plasticity | `advance_plant_traits` | `meds_vegetation_dynamics.f90:114-118` |
| growth commit | `update_cohort_states` + `fill_new_leaf_water` | `:198-211` |
| cohort operators | fuse, cull, split, sort | `:339-354` |
| disturbance | `apply_patch_disturbance` | `:358-366` |
| patch operators | patch fusion, cohort fusion | `:368-380` |
| tissue water | `shed_excess_tissue_water` | `meds_slow_dynamics.f90:114-118` |

**What changes at the brackets:**

- The "after" call becomes `hold_tissue_temps(poly, heat)` in the slow ledger module. It applies `hold_cohort_tissue_temps` to every cohort and returns Σ area·n·(E_leaf + E_wood).
- It runs whether or not the ledger is on. The declaration stays under `present(ledger)`.
- `slow_tissue_heat` becomes a plain sum of the energy.

**`shed_turnover_water`** gets the same bracket. That replaces its two hand-written valuations, for the leaf (`:845-846`) and for the wood's sapwood share (`:856-862`).

**Per-operator effects.** E4 includes an audit of every path that changes a cohort's carbon, water or density, and the §2.2 guard catches one it misses.

| Operator | What happens |
|---|---|
| Fusion | Sums the energy, nplant-weighted beside the water. The temperatures are derived once the carbon pools are blended and the geometry re-derived. It is exact, with no declared remainder. The diagnostic twins (`CO_LEAF_TEMP`, `CO_WOOD_TEMP`) still fuse on the heat-capacity weights. |
| Split | Copies energy and temperature. The daughters' ±eps wood carbon is then held at the bracket's end. Both keep the parent's temperature, and the total is conserved exactly. |
| Mortality | Density falls while the per-plant energy is unchanged. The bracket measures the dead's heat, as today. |
| Birth | Temperature from `patch_physical_state` and water at the birth potential, as now; then the energy. The recruitment heat that `apply_recruitment` declares becomes Σ area·n·(E_leaf + E_wood). |
| New leaves' water (D6) | As today: the water taken from the wood is held at the leaf's temperature, inside the growth exchange. |

### 2.6 Heat capacity: one definition, computed on demand

- **Not stored.** It is a function of the mass.
- **One definition** for each half, the dry heat capacity and the water mass, per m² of ground, in `meds_therm_lib` beside `tissue_heat_capacity`:
  - `build_column_frozen` calls them on the column view.
  - `cohort_tissue_heat_capacity` calls them on the cohort block, and keeps its role for the diagnostics' fusion weights.
  - T(E) is the same per plant or per m², because numerator and denominator both scale with nplant.
- `veg_thermal%c_leaf` and `veg_thermal%c_wood_dry` are deleted. They are not config keys; the constants in `meds_constants` are the one definition.

### 2.7 Restart (D7: both)

- **Write** `leaf_energy` and `wood_energy` [J/plant], with `leaf_temp` and `wood_temp` beside them as diagnosed twins. That is the convention the state file already uses for `soil_temp`/`soil_fliq`, `snow_temp`/`snow_fliq` and `cas_can_temp` (`meds_io.f90:213-221`).
- **Read** both when the energy is present. By the §2.2 invariant, the twin read back is bit for bit what the run carried.
- **Older files are still readable:**
  - A file with temperatures but no energy is read and the energy derived (a hold).
  - A file with neither starts each cohort in its patch's state (`patch_physical_state`), as now.

### 2.8 Output, Python and the C API

- **Output: no change.** `leaf_temp_cohort`, `wood_temp_cohort`, `leaf_temp_site` and the sub-daily rows come from the cohort diagnostic rows, which the prepass fills from the working bundle's temperatures.
- **C API and Python: no change.** The kernels take temperatures, and no accessor exposes the cohort's tissue state or films.
- **The finiteness check** in `meds_diagnostic_reduce.f90:540-544` adds the energies.

### 2.9 Ledgers

- **Fast whole-column:** the tissue stores become Σ E, and the difference is unchanged.
- **Slow:** the store's tissue slice is Σ area·n·E, in the datum of the soil, snow and canopy air. The films become area·n·film.
- **New check:** across many fast steps, each step's Σ E at the start equals the previous step's at the end. Today this fails by ΔW·c_liq·(T₁ − T₀ + T_sc).

### 2.10 The canopy films per plant (D1)

`leaf_surf_water` and `wood_surf_water` become kg per plant; the capacity is `*_surf_water_max` × the plant's leaf or wood area. They become ordinary per-plant amounts:

| Operation | Today (per ground) | Per plant |
|---|---|---|
| Cohort fusion | added | nplant-weighted, beside the water |
| Split | each daughter halved | copied |
| Patch fusion | scaled by the area weight | unchanged |
| Disturbance survivor | scaled by the area factor | unchanged |
| Mortality | stays on the survivors (its own rule) | the dead's film goes to the ground with their tissue water, down the shed channel |
| Cull, disturbance kill | the cohort's film to the ground | n·film to the ground |

**What goes:**
- `scale_cohort_ground_fields` and its five calls;
- the two film-conservation error-stops;
- `FSK_GROUND` as a state kind, and its heading in `fuse_cohort_fast_state`. It stays as a diagnostic kind: the per-ground radiation rows use it.

**The fast tier is unchanged.**
- The seam converts: × nplant in, ÷ nplant out.
- The interception, film evaporation, the film clamp and the film ledger stay per m² of ground.

**What changes numerically:**
- the round-off of n·(x/n) at each seam;
- the dead's film, which now leaves with them. It is small, because a film holds at most 0.1 kg per m² of leaf.

**Restart (D12).**
- The film variables keep their names and get a `units` attribute, `kg/plant`.
- A file without the attribute holds per-ground films, which the reader divides by nplant.

### 2.11 Shed water carries its enthalpy (D11)

**A new patch field,** `shed_water_enth` [W m⁻²], beside `shed_water_rate` [kg m⁻² s⁻¹]. It has the same lifecycle:
- area-weighted in patch fusion;
- accumulated into a disturbance gap;
- in the patch lockstep and the restart.

**Every shed path adds its water's enthalpy:**
- tissue water at u(T) of its tissue;
- films at u(T) of the tissue they sit on.

**The fast tier:**
- It values the ground inflow as the mixture of the precipitation's enthalpy (or the melt's, under a pack) and the shed's, over their combined mass.
- The film's valuation temperature stays the rain's alone, since the film holds rain, not shed water.

**The slow ledger's declared water out and the fast tier's water in are then the same number in energy as well as mass.** What goes:
- the "valued like the precipitation it joins" rule;
- the dry-day fallback to the canopy-air temperature;
- the comment and ledger note that explain the re-valuation (`meds_vegetation_dynamics.f90:812-818, 840-844`).

**Old files** without `shed_water_enth` value the stored rate at the patch's canopy-air temperature. That is one day's shed water, once.

## 3. What this removes

This answers D4's question: yes, it removes special cases, the shed ones included.

| # | Removed | Step |
|---|---|---|
| 1 | The duplicate heat-capacity definition (`veg_thermal%c_leaf`, `c_wood_dry`) | E1 |
| 2 | `scale_cohort_ground_fields` and its 5 calls; the film-conservation error-stops in split and patch fusion | E2 |
| 3 | Mortality's rule that the dead's film stays on the survivors | E2 |
| 4 | `wood_thermal_water`: the sapwood-share thermal water, in the fast capacity, the cohort capacity and `shed_turnover_water` | E3 |
| 5 | Fusion's declared non-additive thermal remainder, and the wood tolerance in `test_fusion_cohort` | E3 |
| 6 | Heat-capacity weighting of the state temperatures in fusion (the diagnostic twins keep it) | E4 |
| 7 | The two hand-written enthalpy valuations in `shed_turnover_water` (leaf, wood share) | E4 |
| 8 | The duplicated tissue commit in ARK and RK45 | E4 |
| 9 | The column-sum water enthalpy that cap·T needed beside it (`water_store_enth`) | E4 |
| 10 | Shed water valued as the inflow it joins; the dry-day canopy-air fallback | E5 |

Added: two cohort fields (the energies), one patch field (`shed_water_enth`), two helpers, two kernels, one guard.

## 4. Files

| Area | Files | Steps |
|---|---|---|
| kernels | `meds_therm_lib` (`tissue_energy`, `tissue_temp`, the capacity halves); `meds_water_retention` (`wood_thermal_water` deleted) | E1, E3, E4 |
| cohort state | `meds_cohort_state_vars`: fields, the lockstep lists, birth, the fusion policy, the helpers, `scale_cohort_ground_fields` deleted | E1–E4 |
| patch state | `meds_patch_state_vars` (`shed_water_enth` and its lockstep) | E5 |
| fast | `meds_plant_types`; `meds_fast_types`; `meds_fast_frozen` (capacity, per-cohort eq. (5), the inflow mixture); `meds_fast_ark` and `meds_fast_rk45` (one shared commit); `meds_fast_dynamics` (seam, guard, `init_fast_reservoirs`) | E1–E5 |
| slow | `meds_slow_ledger`; `meds_vegetation_dynamics` (brackets, the shed paths); `meds_slow_dynamics`; `meds_demography_cohort_fusefiss` (split, cull, recruit heat); `meds_demography_patch_fusefiss` (patch fusion, disturbance) | E2–E5 |
| io | `meds_io` (energies, twins, film units, `shed_water_enth`, old files); `meds_diagnostic_reduce` | E2, E4, E5 |
| tests | `test_fusion_cohort`, `test_slow_ledger`, `test_state_roundtrip`, `test_surface_energy`, `test_fast_loop`, `test_several_sites`; new `test_tissue_energy` | E1–E5 |
| docs | vegetation_energy_dynamics.md; soil_plant_water.md §2 and §5; column_biophysics.md (state table, mixing); diagnostics.md (fusion kinds); numerical_scheme.md; ed2_comparison.md; ROADMAP (#456's assumption); CHANGELOG | E2–E6 |

## 5. Steps and verification

Each step is one verified commit on `feat/tissue-energy`, in a new worktree, delivered as one PR into `beta` (a merge commit).

**Every step runs:**
- ifx Release and Debug, and gfortran: 79/79.
- Python library 80/80, and pytest.
- The 8-case harness (`run_s1.sh`).

| Step | Content | Accepted when |
|---|---|---|
| E0 | This plan goes to `docs/dev_plans/MEDS_TISSUE_ENERGY_PLAN.md`. | — |
| E1 | One heat capacity: `build_column_frozen` takes its dry capacities and water masses from §2.6's functions; the `veg_thermal` duplicates go. | bit-identical on all 8 cases |
| E2 | Per-plant films (§2.10). | differences at round-off, apart from the dead's film; the film ledgers close as before |
| E3 | The whole wood store (§2.4). | the closure and change criteria below |
| E4 | The energy state (§2.2–§2.9). New `test_tissue_energy`: (a) with the mask off, a refilling tissue keeps its temperature; (b) no cross-step re-valuation; (c) fusion sums energy exactly, and a split keeps the parent's temperature; (d) restart round trip, and an old file read through its temperatures. | as below; a resumed run matches a continuous one bit for bit |
| E5 | Shed enthalpy (§2.11). A test checks that the energy handed off equals the energy the fast tier takes in. | as below |
| E6 | ed2_comparison.md, ROADMAP, CHANGELOG, this plan's log. | — |

**For E3 to E5** (any step whose numbers change):
- **Closure no worse.** The whole-column energy and water worst residuals stay within today's tolerances, and no slow-ledger phase's residual grows.
- **Results move little.** Site GPP, LE, H, NEE, AGB and mortality move by less than 0.1 % over the 2-year census case and the ex04 month; the mean |ΔT| of the tissue temperatures is reported.
- These thresholds are expectations, not results. If a case exceeds one, I stop and report before committing.

## 6. Decisions

| # | Question | Status |
|---|---|---|
| D1 | Energy per plant or per m² | **Owner: per plant; the films per plant too** (§2.10) |
| D2 | Datum | taken as recommended: the soil's internal energy, `temp_to_internal_energy` |
| D3 | Keep the temperatures stored? | **Owner: yes**, derived, always exactly T(E, mass) |
| D4 | Wood's thermal mass | **Owner: the whole wood store**, and remove the shed special cases (§3) |
| D5 | Slow operators | taken as recommended: hold the temperature; the energy is re-derived at the existing brackets |
| D6 | New leaves' water from the wood | taken as recommended: held at the leaf's temperature, as today |
| D7 | State file | **Owner: both**, the temperatures as diagnosed twins |
| D8 | Freezing | **Owner: supercooled liquid**; the assumption and cold-driven mortality are issue #456 |
| D9 | The seam guard | taken as recommended: always on, 10⁻⁶ K |
| D10 | Delivery | taken as recommended: one PR, commits E1–E6 |
| D11 | Shed water carries its enthalpy (§2.11) | **new**: recommended; the alternative keeps today's valuation and its special cases |
| D12 | Old per-ground films in state files | **new**: recommended, a `units` attribute on the same variable names; the alternative renames the variables (the `llspan` precedent) |

## 7. Later

**Tissue freeze/thaw is out of scope (D8, #456).**
- Tissue water stays liquid below 0 °C.
- #456 states the assumption and asks for cold-driven mortality: a PFT cold tolerance, frost leaf loss or death, and possibly freeze–thaw embolism.
- The energy state is what such work would build on.

**The films inside the tissue energy** (not proposed now).
- The films are valued at a per-step reference temperature (`t_film_valuation`) rather than carrying an energy, so they re-value between steps as the tissue did.
- Folding each film into its tissue's energy, at the tissue's temperature, would make rain cool or warm the leaves it wets and remove `t_film_valuation`.
- It would also raise a wet leaf's heat capacity by about 10–20 %.
- It is a physics change of its own, for a later plan if wanted.

## 8. Implementation log

**On `feat/tissue-energy`, 2026-10-19,** from `beta` at f6dfa9e.

- **Every step** was built with ifx Release and run on an 8-case harness against a reference built from f6dfa9e: a two-year census run, a one-year run with every output group on, example 04's evaluation month, a Palo Verde inventory site, a four-cell region, a site and a region restarted from state files written by f6dfa9e, and a site restarted from its own new state file.
- **The full matrix** (ifx Debug, gfortran, the standalone libraries, the Python library and pytest) ran on the final commit.
- **There are 80 tests, not 79:** E4 adds `test_tissue_energy`.

| Step | Commit | Content | Site means |
|---|---|---|---|
| E0 | 1d64c9d | This plan | — |
| E1 | 1d7e501 | `tissue_dry_heat_capacity`, used by the fast tier; `veg_thermal%c_leaf` and `c_wood_dry` gone | bit-identical on all 8 cases |
| E2 | ff17e13 | Films per plant: `scale_cohort_ground_fields` and the film error-stops gone, `cohort_film_water`, the dead's films shed with them, `units = "kg/plant"` in the state file | census cases bit-identical; elsewhere at most ~0.1 %, the BCI month 10⁻⁴ %; Palo Verde's shed water +2.8 % (the dead's films) |
| E3 | 6e5ec90 | The whole wood store in the wood's heat capacity; `wood_thermal_water` gone | from f6dfa9e at most 0.05 % (Palo Verde LE +0.04 %, H −0.04 %) |
| E4 | 518fb2c | `leaf_energy`, `wood_energy` as the state; `commit_tissue_energy` shared by ARK and RK45; the seam guard; fusion sums energy | at most 0.007 %; the cohort fuse/fiss phase's declared energy 1.1×10⁻² J → 0 |
| E5 | da3ff31 | `shed_water_enth`; the pond takes the mixture; the budget's `ground_inflow_enth` | below |
| E6 | — | ed2_comparison.md, ROADMAP (#165, #419, #456), CHANGELOG, this log; the plan archived | — |

**Departures from the plan:**

- **The temperatures are held inside each operator** that changes a tissue's mass, not at the ledger brackets (§2.5). This covers `advance_plant_traits`, `reacclimate_plant_traits`, `update_cohort_states`, `fill_new_leaf_water`, the shed paths, split, fusion and `init_cohort`. The seam guard showed why on the first restart: re-acclimating the traits after a restart changes leaf carbon outside any bracket.
- **Reading an old file through its temperatures** (E4 (d)) is tested by the harness, not a unit test. Removing `leaf_energy` from a written file failed in netCDF. The site and region restarts from f6dfa9e's state files take that path, and the guard never fired.
- **Tests that set masses or temperatures by hand** (`test_fast_loop`, `test_slow_ledger`, `test_fusion_cohort`) now hold the temperatures afterwards, as the operators do; without it the seam guard stops them.
- **The E5 test** is in `test_slow_ledger`: the enthalpy the turnover shed declares equals what the fast tier takes in over the day, `shed_water_enth`·dt_slow summed over the patches. `test_carbon_growth` checks the patch rate.

**From f6dfa9e to E5,** site means over each run:

| Case | GPP | LE | H | NEE | AGB | Mortality | Runoff |
|---|---|---|---|---|---|---|---|
| Census, 2 yr | 0.0000 % | −0.0000 % | 0.0000 % | −0.0005 % | −0.0000 % | +0.006 % | 0.0000 % |
| Palo Verde, 1 yr | −0.001 % | +0.038 % | −0.036 % | −0.001 % | 0.0000 % | −0.005 % | −0.017 % |
| Region, 4 cells | −0.0000 % | +0.0004 % | −0.0002 % | +0.002 % | +0.0005 % | −0.007 % | −0.0003 % |
| Region restart | +0.002 % | −0.0007 % | −0.0001 % | −0.004 % | +0.0001 % | −0.0003 % | 0.0000 % |

- **Example 04's month:** half-hourly H moves by −0.02 % on average (1.3 % at most), LE by +0.01 % and GPP by −0.003 %.
- **The tissue temperatures:**
  - Palo Verde's daily site means move by 0.9 mK (leaf) and 1.7 mK (wood) on average, and a cohort's by 0.16 K at most.
  - The census run's move by 10⁻⁶ K.
- **Closure is as before.**
  - Over a site run, the whole-column energy residual sums to at most 3.3×10⁻⁵ J m⁻² (3.7×10⁻⁵ before) and the water residual to at most 3.4×10⁻¹² kg m⁻² (4.0×10⁻¹² before).
  - The worst single check is about 1.5×10⁻⁶ J m⁻², as before, and no check fails.
  - Every slow-ledger phase closes at round-off.
