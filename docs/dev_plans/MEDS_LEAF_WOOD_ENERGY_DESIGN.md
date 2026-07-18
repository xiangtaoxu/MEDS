# MEDS leaf & wood energy — diagnostic/prognostic model option + separate wood store

**Status:** **P0–P3 IMPLEMENTED** on branch `feature/leaf-wood-energy` (commits 8dc56b9 P0, 1d3f40f P1,
84eac0b P2, e800d13 P3; ifx + nvfortran multicore 33/33, prognostic runs bit-identical across both compilers).
**P4 (prognostic leaf under ARK, bordered arrowhead) NOT yet implemented** — cleanly gated (error-stop) and
deferred; resume there. Adversarial review deliberately omitted (later pass).
**Goal:** make **leaf** energy and **wood** energy each a config-selectable **{diagnostic | prognostic}** model,
and give **wood its own temperature in every mode** (today wood is forced equal to leaf temperature). This lets us
test the diagnostic-vs-prognostic difference — most importantly the physically-lagging wood temperature that the
current `wood_temp = leaf_temp` slaving cannot represent.

### Implementation notes / findings (2026-07-13, P0–P3)
- **`veg_energy_balance` was NOT L-stable and was fixed in P2 (84eac0b).** The old store update applied the
  *endpoint* flux `store += dt*r_star`, which is unbounded when the store is stiff (`tau = cap/|drdt| << dt`) and
  `T^n` is far from equilibrium — a young-stand wood cohort seeded at 288 K overshot to 367 K then collapsed to
  −751 K, driving `t_ground` NaN and hanging the adaptive hydrology. It now applies the **bounded** energy change
  `cap*(t_star - t_n)` (t_star is one damped Newton step from t_n) and reports the **step-averaged** turbulent
  flux the store did not keep (`abs - delta_e/dt`). For a linear flux `delta_e == dt*r_star`, so the existing
  surface-energy conservation test is unchanged. This L-stable kernel is what makes both prognostic wood AND a
  future prognostic-leaf-under-ARK viable. MVP wood also gets `gbw=0` (dry bark, no film evap/dew — the dew branch
  made it a condensation surface near saturation) and `cap` floored by the **absolute** `veg_hcap_min`.
- **P3 prognostic leaf uses the BE `cap/dt` term added to the diagnostic linearization** (numerator gains
  `a_leaf*(t_emit - tcas)`, denominator gains `a_leaf = cap_leaf/dt`), so diagnostic (`cap_leaf=0`) is bit-identical
  and the leaf keeps the rich transpiration/emission flux definitions. **Requires `integration_scheme="picard"`**
  (error-stop otherwise): a single explicit split pass (niter=1) is marginally UNSTABLE with the storage term (a
  2·dt oscillation, ~1.7 K midday Tcas spikes on the census stand); the Picard iterate damps it to ~0.2 K. Wood has
  no transpiration feedback and is stable on the pure-split path, so P2 is NOT gated this way.
- **The leaf-inertia effect is negligible at production dt.** `tau_leaf ~ 14 s << dt_fast`, so prognostic ≈
  diagnostic: mean |ΔTcas| = 0.015 K, peak 0.2 K on the census stand. This is the physically-expected result and it
  is *why* diagnostic leaf is the correct default (and why ARK+diagnostic-leaf is fine for production). P4 would
  quantify this same negligible effect under ARK — architectural completeness, not new science.
- **Picard-correctness:** both `t_emit` (leaf) and `wood_emit` (wood) snapshot the start-of-sub-step temperature so
  the prognostic seed is correct across Picard passes (bit-identical at niter=1).

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
  already *is* diagnostic.
- **`wood_energy_model`** (NEW): `diagnostic (0, default)` | `prognostic (1)`.
  - `diagnostic` = separate steady-state wood temperature (own balance, own radiation, own boundary layer,
    sensible + net-LW to CAS).
  - `prognostic` = wood internal-energy store with real thermal mass (lagging diurnal wood temperature).

  There is **no `slaved` mode** — wood is *always* a genuine separate store that exchanges energy with the canopy
  air. Both selectors are binary and **independent** → the test matrix is `{leaf: diag | prog} × {wood: diag | prog}`
  (4 modes). Both are read in `meds_config_io` (missing key → default), range/integrator-support validated in
  `validate_config`, and carried on `column_config_t`.

**Behavior change (the default is NO longer byte-identical to today).** Because there is no `slaved` fallback, the
default (`wood=diagnostic`) replaces today's `wood_temp = leaf_temp`: wood now carries its own temperature and
contributes sensible + net-LW to the CAS. This is an intended **physics correction** — wood↔CAS energy exchange was
previously ignored — but it shifts the default fast-loop result for *every* run (the CAS, hence leaf temp too, moves
slightly; wood is a modest term, `π·WAI·gbh` sensible with `WAI≈0.2·LAI`). Concretely: the split golden anchor
`tc_split(54)=292.450065` must be **re-baselined** to the new default, and existing configs run the new physics with
no edit (missing key → `diagnostic`). If a strict reproduce-today path is ever needed it can be recovered by
`WAI→0` (zero wood area ⇒ zero wood coupling), but that is a test artifact, not a supported mode.

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

Each phase is a green, testable increment. The default becomes `leaf=diagnostic`, `wood=diagnostic` — **not**
byte-identical to today (the split anchor is re-baselined once, in P0), but stable thereafter. Land each mode on the
**split** path first, then the **ARK** path; error-stop any (mode × integrator) combination not yet wired (mirroring
the existing `LEAFEN_PROGNOSTIC` guard at `meds_column_dynamics.f90:246`).

### P0 — Separate DIAGNOSTIC wood (the new default) + scaffolding (both integrators)

Deliverable: wood gets a separately-computed steady-state temperature coupled to the CAS as the **new default**
(`wood_energy_model="diagnostic"`); the split golden anchor is re-baselined. **This is the immediately-testable
increment** ("separate wood temperature").

- **Config:** add `wood_energy_model` to `meds_config_t` + `validate_config` range guard (gated on
  `fast_biophysics_on`); TOML read in `meds_config_io` (`{diagnostic→0, prognostic→1}`, missing→0);
  `WOODEN_DIAGNOSTIC/PROGNOSTIC` params + `column_config_t%wood_energy_model`; copy `cfg→ctx%ccfg`.
- **State:** `cohort_block%wood_temp(:)` through the 8 lockstep sites + WAI-weighted fusion merge;
  `patch_biophys_t%wood_temp(:)` + `alloc_patch_biophys`; fast-loop gather/scatter.
- **Wood radiation harvest:** add `abs_sw_wood(:)/abs_lw_wood(:)` to `column_forcing_t`; in `apply_rt_forcing` map
  `flux%abs_wood(VIS)+abs_wood(NIR)→abs_sw_wood`, `abs_wood(LW)→abs_lw_wood`; const/no-forcing `fill_forcing` sets
  wood absorption = 0 (MVP).
- **Aero at wood temperature:** give `aero_bottom_to_top` a `wood_temp` argument; pass `wt_bt` (lagged start-of-step
  `bio%wood_temp`) as `canopy_aerodynamics`' wood-temp argument at both call sites (today it passes `lt_bt,lt_bt`).
- **Split diagnostic wood balance:** in the `column_fast_step` cohort loop, for `wood_energy_model=DIAGNOSTIC`
  (the default), compute `h_coeff_w = π·wai·wood_gbh·ρ·cp_air`, `lw_slope_w = 4·leaf_emiss·σ·te_w³·wai`,
  `dtw = (abs_sw_wood + abs_lw_wood − lw_slope_w·(tcas−te_w)) / (h_coeff_w + lw_slope_w)`, `twood = tcas + dtw`; set
  `bio%wood_temp(i)`; accumulate `coh_h += h_coeff_w·dtw` and `coh_rnet += abs_sw_wood + abs_lw_wood −
  lw_slope_w·((tcas−te_w)+dtw)`. Change `wenv%wood_temp` (respiration driver) from `leaf_temp` to `bio%wood_temp(i)`.
  (There is no `slaved` branch — wood always runs its own balance; `PROGNOSTIC` is dispatched in P2.)
- **ARK/derivs twin:** carry `wood_gbh/wai/abs_sw_wood/abs_lw_wood` on `surface_frozen_t`, fill in
  `build_column_frozen`, set `wenv%wood_temp=bio%wood_temp`; add the same diagnostic wood branch to `surface_derivs`
  (wood sensible + net-LW into `src_enth`/`coh_rnet`); add `surface_tend_t%wood_temp(:)` and write it back on unpack.
- **Conservation:** the wood sensible+net-LW added to `src_enth` telescopes against `coh_rnet` (same fluxes on both
  sides) — diagnostic wood is conservative with no storage term. Ledger must still close.
- **Reuse note (optional cleanup):** the split (`:359-388`) and `surface_derivs` (`:220-236`) leaf algebra can be
  unified into one `veg_energy_diagnostic` (= `veg_energy_balance` at `cap=0`) that both the leaf and wood diagnostic
  branches call — nice, but must stay **bit-identical** for leaf, so treat as optional and verify against the golden
  anchor.

### P1 — Separate leaf/wood longwave emission temperatures

**Goal:** leaf emits LW at `T_leaf` and wood at `T_wood`, so the canopy LW field — and each element's net LW —
reflects the now-separate tissue temperatures instead of the single blended `tcas` MEDS uses today. This is the
radiative *completion* of P0: P0 splits *absorbed* LW leaf-vs-wood but leaves *emission* blended, an inconsistency
largest on cold, clear, calm nights (when `T_wood`/`T_leaf` diverge most from `tcas`).

The two-stream already emits **per cohort** (`emission(i)=σ·canopy_temp(i)⁴`, `meds_canopy_radiation.f90:84`) — MEDS
just feeds it one blended temperature. So the core change is in the *caller*, not the solver:

- **Per-cohort effective emission temperature.** Build `canopy_temp(i)` from the **lagged** (start-of-sub-step) leaf
  and wood temperatures as the area-weighted radiative mean
  `canopy_temp(i)⁴ = leaf_frac·T_leaf,i⁴ + (1−leaf_frac)·T_wood,i⁴` (i.e. `(LAI·T_leaf⁴ + WAI·T_wood⁴)/(LAI+WAI)`).
  The cohort's total emission into the field is then *exactly* leaf-at-`T_leaf` + wood-at-`T_wood` (the T⁴ weights
  telescope with `leaf_frac`), fixing the inter-cohort / sky / ground LW coupling — the dominant effect. Wire this
  where `canopy_temp` is assembled for `apply_rt_forcing`; **the two-stream solver is unchanged**.
- **Emission base per element — DEFERRED (implementation finding).** The intent was to re-base each balance's
  `lw_slope·(T − te)` on the element's own lagged temperature so emission is "counted once." **This destabilizes the
  single-pass split integrator:** with `te = T_leaf_lag > tcas`, the `+lw_slope·(te − tcas)` term is a positive
  feedback (leaf warms → higher `te` next sub-step → warms more), which drove `test_fast_loop` to a NaN. (Picard/ARK
  are iterative/implicit and would tolerate it, but the split is single-pass.) So the balances **keep their local base
  at `tcas` (split) / `T_leaf_lag` (picard)** — the `tcan_bt` field change above is the stable, dominant win; the
  per-element base is a documented residual (this is the same residual noted below, made concrete). The `tcan_bt`
  change *is* stabilizing (higher emission temp → more negative `abs_lw` → cools → negative feedback).
- **Both integrators, both stores** inherit it (the change is in the shared RT-forcing path): split & ARK, diagnostic
  & prognostic leaf/wood.
- **Conservation & anchor.** Total canopy emission shifts from the blended `σ·tcas⁴` to the physically-correct
  area-weighted `σ·(LAI·T_leaf⁴ + WAI·T_wood⁴)` — a real flux redistribution — so the whole-column energy ledger still
  closes but the split anchor moves again; re-baseline `tc_split(54)` here (or land P0+P1 together and baseline once).
- **Residual approximation (honest).** The area (`leaf_frac`) split of *net* LW is exact for the absorbed field share
  but only approximate for *attributing emission* within a cohort when `T_leaf≠T_wood` (the effective temperature
  makes the cohort *total* emission exact; the leaf-vs-wood split of it stays area-weighted). The fully-exact
  alternative — RT returns *gross* absorbed and each energy balance subtracts its own full `σεT⁴·AI` — is a larger
  RT-interface change; adopt only if the split error proves material (high WAI/LAI with a large leaf–wood gap).
- **Test:** the diurnal FAST-tier comparison gains a separate-emission-vs-blended case — expected to matter most on
  cold clear January nights, lowering nighttime wood/leaf temperatures (more realistic radiative cooling).

### P2 — Prognostic WOOD (operator-split, both integrators)

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

### P3 — Prognostic LEAF (split path)

`leaf_energy_model="prognostic"` on the **split** integrator via the L-stable linearized-BE `cap/dt` term.

- Replace the steady-state `dtl` (`:374`) with the finite-cap form `dtl = (net)/(cap/dt + h_coeff_f + le_slope +
  lw_slope)` (or call `veg_energy_balance(is_leaf=.true.)` seeded from `bio%leaf_temp`); advance `bio%leaf_temp` from
  the returned temp; take `transp_i/coh_h/coh_qw` from the returned fluxes. Leaf stays **inside** the implicit
  leaf↔CAS Picard iterate (the `cap/dt` term stabilizes but does not decouple). Remove the `:246` error stop.
- **Gate** the stiff combination: error-stop at config load if `leaf_energy_model=prognostic` **and** the ARK
  integrator is selected (deferred to P4). The split path is fine.
- Verify `leaf_energy_model=diagnostic` (default) is byte-identical (the `cap/dt→0` limit).

### P4 — Prognostic LEAF under ARK — bordered arrowhead (stretch, deferred) — **RESUME POINT**

**NOT implemented. This is where to resume.** Currently gated: `column_fast_step` error-stops on
`leaf_energy_model="prognostic" .and. time_integrator==INTEG_ARK` ("P4 arrowhead"). ARK runs diagnostic leaf
(correct — leaf inertia is negligible, see findings above). The prognostic-wood-under-ARK combination is likewise
gated (P2 is split-only). Concrete resume steps below; the L-stable `veg_energy_balance` kernel (P2) is the enabler.
The current ARK leaf diagnosis to replace is `surface_derivs` at `meds_column_derivs.f90:230-233`
(`dtl = net/(h_coeff_f+le_slope+lw_slope)` → `f%leaf_temp(i)=tcas+dtl`).

Make prognostic leaf work under IMEX-ARK by solving the stiff leaf implicitly inside the ESDIRK stage.

- In `surface_derivs`, when leaf prognostic, take `T_leaf` as an **input** (the stage iterate) rather than diagnosing
  it, returning H/LE at `T_leaf`.
- Extend `newton_surface_solve`/`jac_surface` to the **bordered arrowhead**: each cohort's leaf is a per-cohort
  1-DOF unknown coupled **only** to the two CAS scalars `(H, q)` — no leaf-leaf coupling — so per-cohort **Schur
  elimination** folds each leaf into a correction on the same **2×2** CAS Jacobian, then back-substitutes leaf temps.
  Wood stays operator-split (from P2). Remove the P3 ARK gate.
- Verify conservation telescopes and L-stability at `dt_fast` (the RK4 oracle blows up where this stays stable).
- Only here, *if* freeze/thaw or true error control on the leaf store is ever wanted, consider adding `leaf_energy`
  to `column_state_t` — otherwise keep the warm-restart.

## 6. Testing — the diagnostic-vs-prognostic difference

- **Unit (`test_surface_energy` / a new `test_veg_energy`):** diagnostic wood steady-state (net flux = 0, transp = 0,
  `twood ≠ tl`); prognostic wood relaxes on `τ ~ minutes–hours`; prognostic leaf relaxes to the diagnostic value
  within 1–2 sub-steps; the `cap→0` limit of prognostic equals diagnostic numerically; wood internal energy conserved
  across cohort fusion.
- **Golden-anchor re-baseline (one-time):** the default is now `leaf=diag`, `wood=diagnostic` (wood↔CAS on), so the
  old split anchor `tc_split(54)=292.450065` no longer applies. Re-baseline `tc_split(54)` to the new value in P0 and
  assert it is (a) stable across ifx/nvfortran and (b) identical between the split and ARK default paths (both run the
  same new default physics). Document the shift magnitude.
- **Integration (the headline test the feature exists for):** a diurnal FAST-tier run (reuse the new `-F-` output)
  comparing the four modes — `wood=diagnostic` vs `prognostic`, crossed with `leaf=diagnostic` vs `prognostic`.
  Expect: **wood diurnal temperature amplitude damped and phase-lagged** under prognostic (large thermal mass) vs
  near-instant equilibration under diagnostic; **leaf** nearly unchanged diag-vs-prog (small thermal mass). Report the
  wood-temp diurnal amplitude/lag and the knock-on stem-respiration and CAS-temperature differences.
- **Portability:** ifx Debug + nvfortran multicore both build/pass; new `sum`/array kernels bound to named vars
  (issue #7); conservation budgets close on both integrators in every mode.

## 7. Key decisions

- **Wood selector is binary** (`diagnostic` default | `prognostic`), matching leaf — wood is *always* a genuine
  separate store; there is no back-compat `slaved` mode. Consequence: the default fast-loop physics changes (wood↔CAS
  coupling is now always on) and the split golden anchor is re-baselined once (a deliberate correctness improvement,
  not a regression). A strict reproduce-today path, if ever needed, is `WAI→0` (a test artifact, not a mode).
- **Diagnostic = `cap/dt→0` of the existing prognostic kernel** — one shared code path; diagnostic is not a separate
  algebraic re-derivation.
- **Persist temperature, reconstruct energy per sub-step** — no `column_state_t` energy fields (warm-restart, matching
  `leaf_temp`/`psi`); only a new `meds_thermo` forward map.
- **Wood operator-split (non-stiff), leaf implicit (stiff)** — wood-prognostic lands first; ARK leaf-prognostic last.
- **Wood radiation = harvest existing `flux%abs_wood`** (absorbed, P0), **plus separate leaf/wood LW emission** (P1)
  via a per-cohort area-weighted effective emission temperature — the two-stream already emits per cohort, so this is
  a caller change, not a solver change.
- **Both integrators honor the selectors**; unimplemented (mode × integrator) combos error-stop at config load.

## 8. Out of scope (MVP) — what is deliberately NOT included, and why

Each item below is a real simplification. **None blocks the diagnostic-vs-prognostic test**, but each bounds how
physically complete the wood/leaf store is. They are the natural follow-ons once the option is in and validated.

1. **Wood water as a prognostic pool.** Tissue thermal mass is `cap = dry_hcap + wmass·cp_liq`, so sapwood *water*
   content matters (it can dominate the wood heat capacity). The MVP sets `wmass_wood` from a **fixed** PFT
   moisture-fraction × wood biomass (or 0), not a water pool that fills/drains over time. *Consequence:* the wood heat
   capacity is static — it won't rise for a well-watered stem or fall as the stem dries. The diurnal *lag shape* is
   right; the *degree of damping* is fixed by the assumed moisture fraction.

2. **Coupling wood `wmass` to the hydraulics PV curve.** MEDS already carries a per-cohort wood water potential
   (`psi(NODE_WOOD)`), so in principle `wmass_wood` could be read from it via the pressure–volume curve, making the
   thermal mass track plant water status. The MVP keeps them **decoupled** (the fixed fraction in #1). *Consequence:*
   the wood thermal store and the hydraulic store evolve independently even though they share the same water; this is
   the obvious P-next if wood-water is wanted.

3. **Freeze/thaw of the veg store.** The soil kernel has an ice-aware freeze/thaw plateau (`ENERGY_PHASE_ON`); the same
   `veg_energy_balance` inverter (`uext_to_temp`) *could* freeze leaf/wood water, but the MVP keeps it
   **always-liquid** (no phase gating), consistent with MEDS's tropical scope. *Consequence — relevant to the cold
   January integrator test:* a sub-freezing wood store is modeled as super-cooled liquid. That is fine for the thermal
   lag and the numerics, but it omits the latent heat of fusion of tissue water (a zero-curtain plateau in the wood
   temperature during freeze/thaw). If Ithaca-January realism matters, this is the first item to revisit.

4. **Real WAI / sapwood allometry.** The wood geometry the heat capacity and boundary layer depend on uses the
   fast-loop **placeholder proxies** `WAI = 0.20·LAI`, `bsap = 0.10·wood_carbon`, `sap_area = 0.05·basal_area` (the MVP
   derived-geometry stubs, already flagged in the fast loop). *Consequence:* the *absolute* wood thermal mass (hence
   the *magnitude* of the prognostic damping) inherits these proxies — treat the prognostic-wood amplitude as
   indicative, not calibrated, until real allometry replaces the stubs. (This is a pre-existing MEDS gap, not new here.)

5. **Restart / serialization of tissue thermal state.** `wood_temp` / `leaf_temp` / `psi` / `gpp_accum` are
   **warm-restarted**, never written to the `-S-` checkpoint. *Consequence:* on restart, wood temperature re-seeds
   (from leaf/CAS) and re-equilibrates over its `τ` — **minutes–hours for wood** (leaf re-equilibrates in a sub-step).
   Harmless for a continuous run; a **frequently-restarted** run loses the wood thermal memory at each restart. Kept
   consistent with today's fast-state convention rather than plumbing new restart fields (which would be the fix).

6. **Full implicit arrowhead over *both* tissues.** In the ARK path, prognostic **leaf** enters the implicit surface
   Newton (P4, via per-cohort Schur), but **wood stays operator-split**. This is deliberate, not a gap: wood is
   non-stiff (`τ ≫ dt_fast`), so an implicit treatment buys no stability and only enlarges the arrowhead. *Consequence:*
   wood is never inside the coupled Newton — the physically-justified choice.

7. **State-vector residence for the stores.** Neither `leaf_energy` nor `wood_energy` is added to `column_state_t` or
   the adaptive-ARK `state_*` machinery (§4 — persist temperature, reconstruct energy per sub-step). *Consequence:* the
   stores are excluded from the ARK embedded-error controller and the WRMS norm (like `psi` today). This is only
   revisited if P4 wants true error control or freeze/thaw on the leaf store.

8. **Adversarial code review** of this plan — a later pass, per the standing preference.
