# MEDS leaf & wood energy — diagnostic/prognostic model option + separate wood store

**Status:** design / implementation plan only (no code yet). Adversarial review deliberately omitted (later pass).
**Goal:** make **leaf** energy and **wood** energy each a config-selectable **{diagnostic | prognostic}** model,
and give **wood its own temperature in every mode** (today wood is forced equal to leaf temperature). This lets us
test the diagnostic-vs-prognostic difference — most importantly the physically-lagging wood temperature that the
current `wood_temp = leaf_temp` slaving cannot represent.

---

## 1. The enabling insight — the infrastructure is ~80 % already here

Two facts, verified in-source, make this mostly a **wiring** job rather than new physics:

1. **The prognostic kernel already exists (unwired).** `veg_energy_balance(store_energy, env, tparams, dt, is_leaf,
   flux)` (`src/biophysics/meds_column_energy.f90:215`) is a linearized **backward-Euler, L-stable** step on a tissue
   internal-energy store, with `is_leaf` selecting **leaf** (transpiration on, `effarea_heat=2` flat plates) vs
   **wood** (no transpiration, `π` cylinders). It reads thermal mass `cap = dry_hcap + wmass·cp_liq`. The per-tissue
   thermal params (`c_leaf=3200`, `c_sapw=2700`, `c_dead`, `c_bark`, `veg_hcap_min`) already live in
   `veg_thermal_params_t`. Nothing in the driver calls it yet.

2. **Diagnostic is the `cap/dt → 0` limit of that same kernel.** The current diagnostic leaf solve
   `dtl = (net flux at tcas) / (h_coeff_f + le_slope + lw_slope)` (`meds_column_dynamics.f90:374`) is *exactly*
   `veg_energy_balance`'s `t_star = t_n + r_n/(cap/dt − drdt)` with `cap/dt → 0` (zero thermal mass / infinite
   conductance), residual evaluated at `t = tcas`, and `−drdt = h_coeff_f + le_slope + lw_slope`. So **diagnostic and
   prognostic are one parameterized code path** — diagnostic = drop the `cap/dt` term.

3. **The RT already splits absorbed radiation leaf-vs-wood.** `canopy_radiation` returns `flux%abs_wood(band,coh) =
   absorbed·(1−leaf_frac)` (`meds_canopy_radiation.f90:98`) alongside `abs_leaf` — MEDS just **discards** `abs_wood`
   in `apply_rt_forcing`. A separate wood store's radiation input is a one-line harvest, not a new two-stream.

What is genuinely missing: a **persistent per-cohort wood temperature**, the **config selectors**, the **wood↔CAS
coupling** (sensible + net-LW), and honoring the selectors in **both** integrator paths (legacy split + IMEX-ARK).

## 2. Config model

- **`leaf_energy_model`** (exists, `meds_config.f90:86`): `diagnostic (0, default)` | `prognostic (1)`. Today's leaf
  already *is* diagnostic, so no third value is needed.
- **`wood_energy_model`** (NEW): `slaved (0, default)` | `diagnostic (1)` | `prognostic (2)`.
  - `slaved` = today's `wood_temp = leaf_temp` with **zero** wood contribution to the CAS → **byte-identical**, the
    golden anchor `tc_split(54)=292.450065` is preserved and every existing config is unchanged.
  - `diagnostic` = separate steady-state wood temperature (own balance, own radiation, own boundary layer,
    sensible + net-LW to CAS).
  - `prognostic` = wood internal-energy store with real thermal mass (lagging diurnal wood temperature).

The tri-state on wood (vs binary on leaf) is *because* today wood≡leaf: `slaved` is the back-compat default; the two
new modes both yield a genuinely separate wood temperature. The two selectors are **independent** → the full test
matrix is `{leaf: diag | prog} × {wood: slaved | diag | prog}`. Both are read in `meds_config_io` (missing key →
default), validated in `validate_config` (range + integrator-support gates), and carried on `column_config_t`.

## 3. The stiffness asymmetry — it drives the phasing

From `veg_energy_balance`: `drdt = −8εσT³·AI − effarea_heat·AI·gbh·ρ·cp − Lv·ρ·g_vapor·dqsatdt` (all ≤ 0),
`τ = cap/|drdt|`.

- **LEAF is STIFF.** `cap_leaf ≈ veg_hcap_min(20) … c_leaf·bleaf`; `|drdt|` is tens–hundreds W/m²/K (transpiration
  term present, LAI-scaled) → `τ_leaf ~ O(1–100 s) ≪ dt_fast`. Prognostic leaf **must be implicit**. In the **split**
  path the per-sub-step linearized-BE (`cap/dt` added to the same implicit leaf↔CAS Picard denominator) is L-stable
  and absorbs it. In the **ARK** path it must enter the surface Newton arrowhead. An explicit/operator-split leaf
  would be unstable at `dt_fast` — this is precisely why leaf defaults to diagnostic.
- **WOOD is NON-STIFF.** `cap_wood ≈ c_sapw·bsap (+dead/bark) ~ 1e3–1e5 J/m²/K`, smaller `|drdt|` (no transpiration,
  WAI<LAI, cylinder area) → `τ_wood ~ minutes–hours ≫ dt_fast`. Prognostic wood is **operator-split** out of the ARK
  arrowhead (like plant hydraulics), advanced once over the full `dt` at the committed-CAS endpoint, and **excluded
  from the embedded-error/WRMS controller**. Cheap, stable, and it *is* the lagging wood temperature the feature is
  about.

**Consequence for phasing:** prognostic **wood** (non-stiff, operator-split) lands **before** prognostic **leaf**
(stiff, implicit); and prognostic-leaf-under-ARK (bordered arrowhead) is the last, hardest piece.

## 4. State design — persist temperature, reconstruct energy (minimal)

Persist only a per-cohort **`wood_temp`** (mirroring the existing `leaf_temp`), **not** a separate internal-energy
state field. For a prognostic store, reconstruct internal energy at the start of each sub-step from the persisted
temperature (`store_energy = tissue_internal_energy(wood_temp, wmass, dry_hcap)`), advance, and read the temperature
back off (`uext_to_temp`). This is legitimate because (a) fast-biophysics state is **warm-restarted, never
checkpointed** (`leaf_temp`/`psi`/`gpp_accum` aren't serialized by `meds_io`), and (b) the prognostic stores are
operator-split and excluded from the ARK error-controlled state vector, so they never need to live in
`column_state_t`. This avoids the heavy `column_state_t` + `state_axpy/accum/wrms/err_diff` expansion entirely.
(The only new `meds_thermo` primitive: `tissue_internal_energy(temp, wmass, dry_hcap)` — the exact forward inverse
of `uext_to_temp`'s liquid branch.)

`wood_temp` joins the cohort SoA lockstep exactly like `leaf_temp`: `cohort_block` decl + all 8 lockstep sites
(alloc/init, `site_free`, `cohort_ensure_capacity`, `move_alloc_block`, `cohort_reorder`, `copy_cohort_slot`), a
**WAI-weighted** (nplant-fallback) merge in `fuse_2_cohorts`, `patch_biophys_t%wood_temp(:)` + `alloc_patch_biophys`,
and the fast-loop gather/scatter (`meds_fast_loop.f90:284` / `:367`).

## 5. Phased implementation

Each phase is a green, testable increment; the default (`leaf=diagnostic`, `wood=slaved`) stays byte-identical on
**both** integrators throughout. Land each mode on the **split** path first, then the **ARK** path; error-stop any
(mode × integrator) combination not yet wired (mirroring the existing `LEAFEN_PROGNOSTIC` guard at
`meds_column_dynamics.f90:246`).

### P0 — Separate DIAGNOSTIC wood + `slaved` default + scaffolding (both integrators)

Deliverable: `wood_energy_model="diagnostic"` gives a separately-computed steady-state wood temperature coupled to
the CAS; default `"slaved"` reproduces today exactly. **This is the immediately-testable increment** ("separate wood
temperature").

- **Config:** add `wood_energy_model` to `meds_config_t` + `validate_config` range guard (gated on
  `fast_biophysics_on`); TOML read in `meds_config_io` (`{slaved→0, diagnostic→1, prognostic→2}`, missing→0);
  `WOODEN_SLAVED/DIAGNOSTIC/PROGNOSTIC` params + `column_config_t%wood_energy_model`; copy `cfg→ctx%ccfg`.
- **State:** `cohort_block%wood_temp(:)` through the 8 lockstep sites + WAI-weighted fusion merge;
  `patch_biophys_t%wood_temp(:)` + `alloc_patch_biophys`; fast-loop gather/scatter.
- **Wood radiation harvest:** add `abs_sw_wood(:)/abs_lw_wood(:)` to `column_forcing_t`; in `apply_rt_forcing` map
  `flux%abs_wood(VIS)+abs_wood(NIR)→abs_sw_wood`, `abs_wood(LW)→abs_lw_wood`; const/no-forcing `fill_forcing` sets
  wood absorption = 0 (MVP).
- **Aero at wood temperature:** give `aero_bottom_to_top` a `wood_temp` argument; pass `wt_bt` (lagged start-of-step
  `bio%wood_temp`) as `canopy_aerodynamics`' wood-temp argument at both call sites (today it passes `lt_bt,lt_bt`).
- **Split diagnostic wood balance:** in the `column_fast_step` cohort loop, when `wood_energy_model≥DIAGNOSTIC`,
  compute `h_coeff_w = π·wai·wood_gbh·ρ·cp_air`, `lw_slope_w = 4·leaf_emiss·σ·te_w³·wai`, `dtw = (abs_sw_wood +
  abs_lw_wood − lw_slope_w·(tcas−te_w)) / (h_coeff_w + lw_slope_w)`, `twood = tcas + dtw`; set `bio%wood_temp(i)`;
  accumulate `coh_h += h_coeff_w·dtw` and `coh_rnet += abs_sw_wood + abs_lw_wood − lw_slope_w·((tcas−te_w)+dtw)`.
  When `SLAVED`, keep `bio%wood_temp = bio%leaf_temp` and add **no** CAS term. Change `wenv%wood_temp` (respiration
  driver) from `leaf_temp` to `bio%wood_temp(i)`.
- **ARK/derivs twin:** carry `wood_gbh/wai/abs_sw_wood/abs_lw_wood` on `surface_frozen_t`, fill in
  `build_column_frozen`, set `wenv%wood_temp=bio%wood_temp`; add the same diagnostic wood branch to `surface_derivs`
  (wood sensible + net-LW into `src_enth`/`coh_rnet`); add `surface_tend_t%wood_temp(:)` and write it back on unpack.
- **Conservation:** the wood sensible+net-LW added to `src_enth` telescopes against `coh_rnet` (same fluxes on both
  sides) — diagnostic wood is conservative with no storage term. Ledger must still close.
- **Reuse note (optional cleanup):** the split (`:359-388`) and `surface_derivs` (`:220-236`) leaf algebra can be
  unified into one `veg_energy_diagnostic` (= `veg_energy_balance` at `cap=0`) that both the leaf and wood diagnostic
  branches call — nice, but must stay **bit-identical** for leaf, so treat as optional and verify against the golden
  anchor.

### P1 — Prognostic WOOD (operator-split, both integrators)

`wood_energy_model="prognostic"` advances a genuine wood internal-energy store via the existing
`veg_energy_balance(is_leaf=.false.)`; non-stiff → operator-split, warm-restarted through `wood_temp`.

- Add `tissue_internal_energy(temp, wmass, dry_hcap)` to `meds_thermo` (forward inverse of `uext_to_temp`).
- Add `pure veg_heat_capacity(is_leaf, biomass_c, nplant, water_kg, tparams) → (dry_hcap, wmass)`: wood `dry_hcap =
  (bsap·c_sapw [+dead·c_dead + bark·c_bark])·nplant·c2b` floored by `veg_hcap_min`, `wmass` from a fixed
  PFT-moisture-fraction × wood biomass (MVP; or 0); leaf `dry_hcap = bleaf·nplant·c2b·c_leaf`, `wmass` = leaf water.
- **Split:** after the CAS Picard commit, per cohort build a wood `leaf_energy_env_t` (is_leaf=.false., area=wai,
  gbh/gbw=aero%wood_*, abs_sw/lw=harvested, can_temp/shv=committed CAS, dry_hcap/wmass from `veg_heat_capacity`),
  seed the store from `bio%wood_temp` via `tissue_internal_energy`, call `veg_energy_balance`, set
  `bio%wood_temp=flux%temp`, and add `flux%h_flux + flux%qw_flux` into the CAS enthalpy source at the committed
  endpoint (single-flux telescoping). Feed stem respiration from the updated wood temperature.
- **ARK:** add `advance_wood_energy_full` mirroring `advance_hydraulics_full` — seed from `bio%wood_temp` at the
  committed-CAS endpoint, one `veg_energy_balance` step over full `dt`, write `bio%wood_temp`; pass `wood_temp`
  through `column_be_stage` untouched (like `psi`); **exclude** from `state_err_diff`/`state_wrms`.
- **Ledger:** add the wood storage delta (`cap_wood·ΔT_wood`) to the whole-energy budget (a slice of `coh_rnet` now
  becomes tissue storage rather than a CAS pass-through) so it closes.

### P2 — Prognostic LEAF (split path)

`leaf_energy_model="prognostic"` on the **split** integrator via the L-stable linearized-BE `cap/dt` term.

- Replace the steady-state `dtl` (`:374`) with the finite-cap form `dtl = (net)/(cap/dt + h_coeff_f + le_slope +
  lw_slope)` (or call `veg_energy_balance(is_leaf=.true.)` seeded from `bio%leaf_temp`); advance `bio%leaf_temp` from
  the returned temp; take `transp_i/coh_h/coh_qw` from the returned fluxes. Leaf stays **inside** the implicit
  leaf↔CAS Picard iterate (the `cap/dt` term stabilizes but does not decouple). Remove the `:246` error stop.
- **Gate** the stiff combination: error-stop at config load if `leaf_energy_model=prognostic` **and** the ARK
  integrator is selected (deferred to P3). The split path is fine.
- Verify `leaf_energy_model=diagnostic` (default) is byte-identical (the `cap/dt→0` limit).

### P3 — Prognostic LEAF under ARK — bordered arrowhead (stretch, deferred)

Make prognostic leaf work under IMEX-ARK by solving the stiff leaf implicitly inside the ESDIRK stage.

- In `surface_derivs`, when leaf prognostic, take `T_leaf` as an **input** (the stage iterate) rather than diagnosing
  it, returning H/LE at `T_leaf`.
- Extend `newton_surface_solve`/`jac_surface` to the **bordered arrowhead**: each cohort's leaf is a per-cohort
  1-DOF unknown coupled **only** to the two CAS scalars `(H, q)` — no leaf-leaf coupling — so per-cohort **Schur
  elimination** folds each leaf into a correction on the same **2×2** CAS Jacobian, then back-substitutes leaf temps.
  Wood stays operator-split (from P1). Remove the P2 ARK gate.
- Verify conservation telescopes and L-stability at `dt_fast` (the RK4 oracle blows up where this stays stable).
- Only here, *if* freeze/thaw or true error control on the leaf store is ever wanted, consider adding `leaf_energy`
  to `column_state_t` — otherwise keep the warm-restart.

## 6. Testing — the diagnostic-vs-prognostic difference

- **Unit (`test_surface_energy` / a new `test_veg_energy`):** diagnostic wood steady-state (net flux = 0, transp = 0,
  `twood ≠ tl`); prognostic wood relaxes on `τ ~ minutes–hours`; prognostic leaf relaxes to the diagnostic value
  within 1–2 sub-steps; the `cap→0` limit of prognostic equals diagnostic numerically; wood internal energy conserved
  across cohort fusion.
- **Golden-anchor regression:** default (`leaf=diag`, `wood=slaved`) is bit-identical on both integrators
  (`tc_split(54)=292.450065`; `test_fast_loop`).
- **Integration (the headline test the feature exists for):** a diurnal FAST-tier run (reuse the new `-F-` output)
  comparing the four modes — `wood=slaved` vs `diagnostic` vs `prognostic`, and `leaf=diagnostic` vs `prognostic`.
  Expect: **wood diurnal temperature amplitude damped and phase-lagged** under prognostic (large thermal mass) vs
  near-instant equilibration under diagnostic vs identical-to-leaf under slaved; **leaf** nearly unchanged
  diag-vs-prog (small thermal mass). Report the wood-temp diurnal amplitude/lag and the knock-on stem-respiration
  and CAS-temperature differences.
- **Portability:** ifx Debug + nvfortran multicore both build/pass; new `sum`/array kernels bound to named vars
  (issue #7); conservation budgets close on both integrators in every mode.

## 7. Key decisions

- **Wood selector is tri-state** (`slaved`/`diagnostic`/`prognostic`); leaf stays binary — reconciles "separate wood
  in both modes" with "existing configs unchanged" (today wood≡leaf).
- **Diagnostic = `cap/dt→0` of the existing prognostic kernel** — one shared code path; diagnostic is not a separate
  algebraic re-derivation.
- **Persist temperature, reconstruct energy per sub-step** — no `column_state_t` energy fields (warm-restart, matching
  `leaf_temp`/`psi`); only a new `meds_thermo` forward map.
- **Wood operator-split (non-stiff), leaf implicit (stiff)** — wood-prognostic lands first; ARK leaf-prognostic last.
- **Wood radiation = harvest existing `flux%abs_wood`**; keep RT LW *emission* blended at `tcas` (only split the
  *absorbed* net LW by `leaf_frac`) — a two-emission-temperature RT is out of scope.
- **Both integrators honor the selectors**; unimplemented (mode × integrator) combos error-stop at config load.

## 8. Out of scope (MVP)

Wood-water as a real prognostic pool (fixed moisture-fraction or 0); two-emission-temperature canopy RT (wood emits
`σT_wood⁴`); coupling wood `wmass` to the hydraulics PV curve; freeze/thaw (`ENERGY_PHASE_*`) gating of the veg store
(tropical scope — `uext_to_temp` read-off stays always-on); replacing the `0.20/0.10/0.05` WAI/bsap/sap_area MVP
proxies with real allometry (the wood heat capacity inherits that placeholder status); serializing tissue internal
energy (warm-restart, like `leaf_temp`/`psi`); the full arrowhead over *both* leaf and wood implicitly (wood stays
operator-split); adversarial review (a later pass).
