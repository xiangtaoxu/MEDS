# Bringing ARK and RK45 up to `split`'s physics

**Status:** 2026-07-30. **Phases 1 and 2 IMPLEMENTED** on `feature/integrator-physics-parity`
(plus an unplanned units fix found during Phase 1, recorded below). **Phases 0, 3, 4 and 5 are NOT
implemented** — see §7 for exactly where each stopped and what remains.

**Question this answers:** `docs/science/numerical_scheme.md` §4 lists three rows where the three
fast-loop integrators solve *different equations*: per-layer root-sink placement, prognostic
leaf/wood energy, and the non-free-drain bottom boundary. What does it take to remove all three, in
what order, and what does each one actually cost?

**Companions.** `MEDS_INTEGRATOR_PARITY.md` is the difference *inventory* (rows 3, 4 and 5 there are
these three, all marked "C6 — deferred"). `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` is the *default-choice*
plan. This file is the *closure* plan for the three remaining Class-1 rows, and it supersedes the "C6
deferred" disposition. `MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §5 P2/P4 is the pre-existing scoping for
row 4 and is folded in below rather than restated.

---

## 1. What the re-verification found

Verified by code reading against `main @ bf46030` (2026-07-30). The §4 table is accurate as far as it
goes, but four things about these rows are not recorded anywhere and they change the shape of the work.

**V1 — all three rows are inert at default settings, so none of them confounds the P1-a comparison.**
`multilayer_roots = .false.` (`meds_config.f90:124`), `leaf_energy_model = wood_energy_model = 0`
(diagnostic, `meds_config.f90:151-152`), `bottom_bc = SOIL_BC_FREE_DRAIN` and `zeng_decker = .false.`
(`meds_biophysics_opts.f90:52-55`). At defaults `split` takes the `root_frac` branch too, and neither
ARK nor RK45 can reach its `error stop`. **This work buys feature coverage under ARK/RK45; it does not
move the residual ~0.45 K family gap**, which lives in row 6 (the frozen soil-water solve) and is
`MEDS_PRODUCTION_INTEGRATOR_PLAN.md` P1-a. Nothing here should be sold as unblocking that.

**V2 — row 5 is asymmetric, and the ARK guard is nearly stale.** The ARK commits the scratch
`column_hydrology_flux` θ verbatim (`meds_fast_ark.f90:952`) together with the scratch's
`w_aquifer` / `z_wt` (`:1016-1017`, already plumbed as `fro%w_aquifer1` / `fro%z_wt1`). Its soil water
*is* the split's own solve, so the aquifer / bedrock / Zeng–Decker **mass** physics is already fully
present on that path — the guard at `meds_fast_ark.f90:899` mostly guards a capability the path has.
RK45 is the opposite case: it integrates its own θ but commits the *scratch's* `w_aquifer1`
(`meds_fast_rk45.f90:554`) — recharge from one trajectory, storage from another. That is the
borrow-a-flux-while-committing-own-state class (`MEDS_INTEGRATOR_PARITY.md` §7 item 5, issue #78
items 3 and 4) and it would go live the instant the guard is lifted. **The guard is load-bearing for
RK45 and near-vestigial for ARK.** *(Phase 0 dissolves this asymmetry rather than resolving it: with no
aquifer store there is nothing for RK45 to borrow.)*

**V3 — `SOIL_BC_AQUIFER` is not an aquifer boundary condition, and its bookkeeping is broken on `split`
too.** A real unconfined aquifer sets ψ = 0 at the water table, making the bottom flux head-driven and
**two-way** (capillary rise when the table is close and the soil above is dry). MEDS applies
`qbot = kk(n)` — the deep-water-table *limit* — unconditionally (`meds_soil_water.f90:419`, `:463`,
`:521`), identically to free-drain. `z_wt` reaches the column only through the Dunne `f_sat` split
(`:100`) and the Zeng–Decker **interior**-face gravity term (`:571`, `k = 1..n-1`, off by default); the
bottom face never sees it. So the water table can rise *into* the soil column —
`test_column_hydrology.f90:395` asserts exactly that — while the bottom layers keep free-draining as
though it were far below. **The missing head term is the negative feedback that bounds a water table**,
and without it the surviving feedbacks are both capped, so `z_wt` can pin at the surface with
`w_aquifer` growing without bound.

Three bookkeeping defects ride on top of the same subsystem, all of which Phase 0 deletes rather than
repairs: all three schemes debit the soil bottom face with `hflux%drainage` (`meds_fast_split.f90:922`,
`meds_fast_ark.f90:183`, `meds_fast_time_derivs.f90:338`), which under the aquifer BC is *baseflow*
(`:152`) while the mass leaving layer *n* is *recharge* (`:146`); the mass ledger counts the aquifer as
a store **inside** the column (`:261-262`, `:277`) while the energy side treats it as **outside**; and
`col%w_aquifer`'s `max(0, ...)` clip (`:149`) can fabricate water. **Reaching "the same physics as
`split`" on this row by copying `split` would propagate all of it**, which is why Phase 0 comes first.

**V4 — RK45 would honour a Zeng–Decker config silently wrongly.** `column_frozen_t%psi_e(:)` exists
(`meds_fast_types.f90:438`) and `column_derivs` passes it into `soil_water_time_deriv`
(`meds_fast_time_derivs.f90:279`), but `build_column_frozen` hardwires it to zero
(`meds_fast_ark.f90:1488`). Lifting the guard without populating it produces plausible numbers from a
scheme that is quietly not doing the thing the config asked for — the exact failure mode
`[energy].debug_error`'s missing TOML reader was (`MEDS_PRODUCTION_INTEGRATOR_PLAN.md` P0-e).

**V5 — row 3 is already half-built on the ARK path.** `build_column_frozen` computes
`root_uptake_layer_b` (`meds_fast_ark.f90:1546`), the per-(layer, cohort) breakdown it needs, and then
discards it, using the static `root_frac` at `:1572`. The state carrier `patch_biophys_t%root_sink_share`
already exists (`meds_biophysics_types.f90:438`). Row 3 is four call sites and one frozen field.

Two smaller items, both worth fixing in passing: the harness's `--parity` preset
(`scripts/numerics_sweep.py:92-109`) pins `soil.bottom_bc` but **not** `soil.zeng_decker`; and
`split`'s `root_sink_share` is applied one pass *late* (written at `meds_fast_split.f90:808-820`, read
at `:770` and `:973`) even though the current step's breakdown is available at `:710`, before either
read. The lag reads as historical rather than intentional.

---

## 2. Decisions taken

1. **Full feature parity on all three rows.** Any configuration valid on `split` becomes valid on
   `ark` and `rk45`, including prognostic leaf under the ARK bordered arrowhead.
2. **Fix the shared kernel before lifting any guard.** `SOIL_BC_AQUIFER` is not currently an aquifer
   boundary condition (V3), so lifting the ARK/RK45 guards onto it as it stands would extend a
   half-built BC from one path to three. Phase 0 rebuilds it in `meds_soil_water`, so all three schemes
   gain it at once, and **as a pure boundary condition with no storage**: a two-way head-driven flux
   against a saturated zone, the water table *defined* as the column base, and the whole aquifer
   bucket — storage, baseflow, Dunne `f_sat`, and their five parameters — deleted rather than adapted.
   *Corollary:* with no store there is no aquifer enthalpy to track. Inflow arrives at `T_aq = T_n`,
   which the retained `geothermal = 0` thermal BC (`dT/dz = 0` at the base) does not merely permit but
   **requires**.
3. **Current-step root shares on all three paths; the lag goes.** `split` changes too. This removes a
   state dependence, makes the three trivially identical by construction rather than by a matched
   convention. Combined with decision 4 this re-baselines the `split` anchor for **every** run.
4. **`[hydraulics].multilayer_roots` is deleted; per-layer root coupling is unconditional.** Row 3 is
   closed by removing the branch rather than implementing both sides of it on three schemes — the same
   move as `[fast].snow_on` (§3f F-1) and `with_mass`/`with_theta` (F-2), and for the same reason: a
   flag whose OFF path is the cruder physics gets plumbed everywhere and eventually reads as "this
   cannot happen". Consequence: Phase 1 changes the default answer for every run.
5. **Design first.** This document; no code until it is agreed. *(Superseded — implementation began
   2026-07-30; see the status line and §7.)*

---

## 3. Phases

Each phase is independently green and testable. Phase 0 unblocks 2 and 3; Phase 1 is independent and can
run in parallel.

**Two phases change physics; the rest must be bit-identical at default configuration.**

- **Phase 0** changes answers under `SOIL_BC_AQUIFER` only — a configuration nothing currently runs,
  which is exactly why the defect survived. Every default run stays bit-identical.
- **Phase 1** changes **every** default run: deleting `multilayer_roots` makes the per-layer root
  coupling unconditional. The `split` golden anchors are re-baselined once, with a documented
  before/after.

Each lands on its own PR, not folded into a parity change. For **every other phase**, a moved golden
anchor in a default run is a bug in the phase rather than a re-baselining — the features are all off by
default, so that is a strong and cheap check. State which case a PR is in.

### Phase 0 — `SOIL_BC_AQUIFER` becomes a real boundary condition, and nothing more

*Prerequisite for Phase 3. Replaces the earlier "aquifer store holds heat" and "aquifer outside the
ledger" drafts — both are moot once the aquifer stops being a store at all. **The only phase in this
plan that changes physics**, and only under `SOIL_BC_AQUIFER`, which is off by default, so acceptance
criterion 2 still holds. Land it on its own PR.*

**Design decisions (2026-07-30), settled.**

1. `free_drain` bottom water flux: unchanged (`q = K(θ_n)`, unit gradient).
2. `bedrock` bottom water flux: unchanged (`q = 0`).
3. `aquifer`: a **two-way head-driven flux against a saturated zone at the column base**. **No aquifer
   storage is tracked** — it is a boundary condition, not a reservoir.
4. **The water table is the column base by definition** (`z_wt ≡ z_bottom`): the soil column *is* the
   unsaturated zone above it. No water-table state, no parameter, no seasonality.
5. Bottom **thermal** flux: unchanged, `geothermal = 0` on every BC. A geothermal option is
   **deferred** to future development (see §5 for what it would need to address).
6. Boundary fluxes carry enthalpy into the ledger (below).

**The flux.** With ψ_wt = 0 at `z_wt = z_bottom` and `Δ = z_n − z_bottom = dz(n)/2`:

```math
q_{\mathrm{bot}} = K_{\mathrm{bot}}\left[\frac{\psi_n}{\Delta} + 1\right]
\qquad (\text{down-positive, } \psi \text{ in metres of head})
```

`ψ_n = ψ_wt` recovers free drainage `q = K`; `|ψ_n| > Δ` reverses it upward. Replaces `qbot` under
`SOIL_BC_AQUIFER` **only**, at all three sites (`meds_soil_water.f90:419`, `:463`, `:521`); free-drain
and bedrock branches untouched and bit-identical.

- **`K_bot` upstream-weighted**, matching the interior faces' own rule (`:564-568`): `K(θ_n)` downward,
  **`K_sat`** upward out of the saturated zone. `K(θ_n)` for upward flow would silently under-predict
  capillary rise, which is the point of the phase.
- **No degenerate case survives.** `Δ = dz(n)/2 > 0` always, so the `z_n − z_wt` sign flip that the
  storage-based `z_wt` made possible cannot occur.
- The two implicit sites need ψ_n in the Jacobian or the head term is a lagged source; RK45's explicit
  `soil_water_time_deriv` (`:521`) gets it naturally. State which was chosen.

**Enthalpy (decision 6), and why the inflow temperature is now forced rather than chosen.** Outflow
leaves upwind at `u_liq(T_n)` — as today, but now on the flux that actually crosses the face. Inflow
arrives at `T_aq`, and **decision 5 fixes `T_aq = T_n`**: `geothermal = 0` *is* the statement
`dT/dz = 0` at the base, whose only consistent companion is `T_aq = T_n`. Earlier drafts argued for
this on accuracy grounds; it is now required for the two boundary treatments to describe the same
sub-column body. Consequences: upward advection is thermally neutral to `T_n` (water arrives at the
layer's own temperature) but does raise its heat capacity, and the ledger books
`mass · u_liq(T_n)` as a boundary input. *This coupling is what a future geothermal option must
revisit — see §5.*

**What gets DELETED, not adapted.** Decisions 3 and 4 make a whole subsystem dead. Deleting it is the
phase's main risk reduction, since anything left behind will be read by something eventually:

- `soil_column_t%w_aquifer` and `%z_wt` (`meds_column_state_types.f90:59-60`), their `blend_soil_w`
  lines (`:207`), and their restart variables (`meds_io.f90:332`, `:401`, `:420`) — **a deliberate
  state-file change**; decide the older-`-S-` read path explicitly.
- `SPECIFIC_YIELD` (`meds_soil_water.f90:32`) and the whole baseflow bucket (`:145-153`).
- `q_drai_max`, `f_drai` — the bucket's parameters, now unreachable.
- **Dunne `f_sat` and its `f_max` / `f_over`** (`:99-104`), per the settled decision. The two-way flux
  lets the column saturate from below *explicitly*, and the existing saturation clip → pond → Horton
  overflow already turns that into runoff; keeping the parameterization on top would double-count the
  process it was standing in for. `f_sat` is computed only in the aquifer branch, so all three
  parameters die with it.
- `column_frozen_t%w_aquifer1`, `%z_wt1` (`meds_fast_types.f90:399-400`) and their commits
  (`meds_fast_ark.f90:1016-1017`, `meds_fast_rk45.f90:554`).

**Two unifications fall out.** `compute_psi_e` already used `z_wt = z_bottom` for every non-aquifer BC
(`:70-76`); with `z_wt ≡ z_bottom` that branch collapses and **Zeng–Decker becomes identical across all
three BCs**. And the aquifer BC now differs from free-drain in **exactly one place** — the bottom flux
— where today it differs in two side channels and not in the flux at all.

**Two consequences to accept knowingly.**

- **`aquifer` becomes the wet-site BC.** With the table pinned at the column base and
  `Δ = dz(n)/2 ≈ 0.27 m` on the default grid, the flux reverses whenever `|ψ_n| > 0.27 m` — i.e. for
  anything drier than near-saturation. The bottom layer will sit near saturation permanently, fed by an
  unlimited supply. That is *correct* for a site whose water table really is at the column base
  (riparian, floodplain, wetland) and wrong for a typical upland stand. Say so in the config docs; it is
  not a drop-in alternative to `free_drain`.
- **It may add a stiff mode.** The bottom-face conductance `K_bot/Δ` with a small `Δ` is a fast
  boundary term. The implicit paths absorb it; **RK45 integrates θ explicitly and may not**. Measure
  its sub-step count on a wet aquifer cell before Phase 3 concludes — do not assume.

**Test.** (i) A dry column over the saturated base shows **upward** flux — assert the sign; this
behaviour does not exist today. (ii) Bottom layer converges toward hydrostatic equilibrium with the base
(`ψ_n → −Δ`) under no other forcing. (iii) Mass **and** energy close to machine precision with the flux
in both directions. (iv) Free-drain and bedrock cells bit-identical. (v) Replace
`test_column_hydrology.f90`'s `test_aquifer` and `test_dunne` — both exercise the deleted bucket, and
`:395` (`'water table rises above bottom'`) asserts a symptom of the defect being removed.

### Phase 1 — row 3: delete `multilayer_roots`; per-layer root coupling becomes unconditional

*Independent of every other phase. **The second phase that changes default physics** (with Phase 0) —
see acceptance criterion 2 in §4.*

Per decision 4, row 3 is closed by deleting the branch rather than implementing both sides of it on
three schemes. Row 3 then leaves the `--parity` inventory not because both paths implement the option
but because **the option stops existing**, which is the cleaner end state.

**This changes the default answer.** The single root-fraction-weighted BC is replaced everywhere by
per-layer ψ_soil with rhizosphere conductances weighted by `K(θ_k)`, and the soil sink follows realized
uptake rather than the static `root_frac` profile. On a well-mixed column the difference is small; on a
drying one — where a dry surface layer's conductance collapses and uptake shifts downward — it is not.
That is exactly the P1-b dry-down regime, so **re-baseline the `split` golden anchors once, with a
documented before/after**, as the leaf/wood P0 landing already did.

**Delete, do not branch:**

- `[hydraulics].multilayer_roots` — the config field (`meds_config.f90:124`), its TOML reader
  (`meds_config_io.f90:737`), `ccfg%multilayer_roots` (`meds_fast_types.f90:94`) and its assignment
  (`meds_fast_dynamics.f90:112`).
- Every `if (ccfg%multilayer_roots)` guard, which becomes unconditional: `meds_fast_split.f90:678`,
  `:727`, `:770`, `:808-820`, `:973`; `meds_fast_ark.f90:1523`, `:1572`.
- The `multilayer_roots` argument to `solve_plant_water_batch` (`meds_plant_hydraulics.f90:361-411`) —
  the fast loop is its only caller. **Keep the per-cohort kernel's `n_root_layer <= 1` scalar path**
  (`:91-109`): `src/plant` is a standalone library that must compile and be usable without a soil
  column, so that branch is not dead, it is just no longer reached from the driver.
- Whatever that leaves unreferenced on the driver side — audit `ccfg%rhizo_cond` / `fro%rhizo_cond` /
  `[hydraulics].rhizo_cond` and the `soil_psi_scalar` argument, all of which feed only the retired
  single-BC path. Note `root_weighted_psi` itself **survives**: `meds_fast_ark.f90:1614` reuses it for
  the root-weighted soil *temperature*.

**Then make the three schemes agree, on the current step's shares:**

- **`split`.** Compute a local `root_share(1:nsl)` immediately after the hydraulics batch solve
  (`meds_fast_split.f90:710`) from this pass's `root_uptake_layer_b` — same non-negative floor, same
  `nplant` weighting, same normalisation as `:809-820` — and use it at `:770` (mass sink) and `:973`
  (heat sink). The lagged read of `bio%root_sink_share` goes, and so does the field itself
  (`meds_biophysics_types.f90:438`) once nothing lags — after confirming `meds_io` does not serialise
  it and no output-registry entry reads it.
- **ARK/RK45.** Add `root_share(n_soil_layer_max)` to `column_frozen_t`, populated in
  `build_column_frozen` from the `root_uptake_layer_b` it already computes (`meds_fast_ark.f90:1546`)
  by the identical floor/weight/normalise. Four sites then read it instead of `root_frac`:
  `meds_fast_ark.f90:1572` (scratch hydrology sink), `meds_fast_ark.f90:165` (`column_be_stage` root
  heat sink), `meds_fast_time_derivs.f90:277` (in-tableau soil-water sink) and `:304` (`column_derivs`
  root heat sink).
- **Keep the degenerate fallback.** When the accumulated share is ~0 (no uptake anywhere), fall back to
  `root_frac`. With the flag gone this is the only remaining branch, so it must not be dropped along
  with the others.

**Two things to confirm rather than assume.**

- **HR stays off.** Always-on per-layer coupling must not enable hydraulic redistribution. The
  non-negative floors on both sides (`meds_fast_split.f90:727`, `meds_plant_hydraulics.f90:251`) are
  what guarantee that; verify they survive the de-branching, because CLAUDE.md is explicit that HR is
  deferred.
- **Cost.** `rhizo_cond_all(nsl, n)` is now built every sub-step on the hot path, and the current loops
  recompute `k_theta` inside the cohort loop although it depends only on the layer
  (`meds_fast_split.f90:682-692`, `meds_fast_ark.f90:1524-1532`). Hoist it, and measure the fast-loop
  wall time before and after — this is the one place the deletion could cost something.

**Test.** From one identical state, one `dt_fast`, the three schemes must place the same per-layer sink
profile — compare layer-wise Δθ with the other stores masked off, not just the column total (a column
total agrees even when the vertical distribution does not, which is the whole failure mode). Add the
cell to `test_column_dynamics`, `test_column_ark` and `test_column_rk45`. Plus `sum(share) == 1` to
round-off, and a dry-down fixture where the surface layer is at wilting point and the deep layers are
wet: uptake must shift downward relative to `root_frac`, which is the behaviour being made default and
the one no existing test covers.
### Phase 2 — row 5a: bedrock and Zeng–Decker, which carry no prognostic state

*Depends on nothing; sequenced after Phase 0 only so the guards are narrowed once.*

- Have `column_hydrology_flux` **return** the `psi_e(:)` it computed (`meds_soil_water.f90:80`) rather
  than exporting `compute_psi_e` for a second, independent evaluation. One authority; a second
  `z_wt` diagnosis on the frozen path is exactly how two trajectories start to disagree. Copy it into
  `fro%psi_e`, replacing the hardwired zero at `meds_fast_ark.f90:1488` (fixes V4). No-op for the ARK,
  whose θ is passed through the tableau; load-bearing for RK45.
- Bedrock needs nothing further: `face_and_sink` and `soil_water_time_deriv` already honour it
  (`meds_soil_water.f90:463`, `:521`), and the ARK inherits it through the scratch solve.
- Narrow both guards from `zeng_decker .or. bottom_bc /= SOIL_BC_FREE_DRAIN` to
  `bottom_bc == SOIL_BC_AQUIFER` (`meds_fast_ark.f90:899`, `meds_fast_rk45.f90:449`). Keep the message
  shape and name the remaining restriction precisely.
- **Fix the mislabelled fixture** while here: `test_column_rk45`'s saturated case calls itself a
  "sealed bedrock column" and imports `SOIL_BC_BEDROCK` but never sets `bottom_bc`
  (`MEDS_INTEGRATOR_PARITY.md` §C2). It is a free-draining column. Wire the intent, and expect the
  numbers to move — that is the fixture starting to test what it claims.

**Test.** A genuinely sealed bedrock cell on all three schemes: zero bottom-face mass **and** zero
bottom-face enthalpy, whole-column ledgers closing, and the three agreeing within a few percent over
a wet month. A Zeng–Decker cell on RK45 that is measurably *different* from the same run with
`zeng_decker = false` — otherwise the phase has not proven `psi_e` is reaching the tendency.

### Phase 3 — row 5b: unlock the aquifer BC on ARK and RK45

*Depends on Phase 0 and Phase 2. Much smaller than the earlier drafts, because Phase 0 deleted the
subsystem those drafts were trying to make three schemes agree about: with no `w_aquifer`, no `z_wt`
state and no `f_sat`, there is nothing left to inherit, freeze, lag, or double-count.*

- **ARK — guard removal.** Its θ comes from the scratch `column_hydrology_flux` verbatim (V2), which
  now carries the head-driven flux, so the physics arrives for free. Narrow the guard at
  `meds_fast_ark.f90:899` to nothing and delete the now-dead `fro%w_aquifer1` / `fro%z_wt1` commits
  (`:1016-1017`). Run the aquifer cell with `[energy].debug_error` armed before declaring it done.
- **RK45 — guard removal plus one measurement.** The head term reaches its tendency through
  `soil_water_time_deriv` (`meds_soil_water.f90:521`) once Phase 0 lands, and the borrowed-`w_aquifer1`
  defect that made this guard load-bearing no longer has an object. Delete the guard
  (`meds_fast_rk45.f90:449`) and the inheritance (`:554`).
  **The one open risk:** the bottom-face conductance `K_bot/Δ` with `Δ = dz(n)/2` is a fast boundary
  term that the implicit paths absorb and an explicit march may not. **Measure RK45's sub-step count
  and `work_rk45_rescue_site` on a wet aquifer cell** before calling the BC supported there. If it
  forces a rescue storm, the honest outcome is a documented tier limit, not a new fallback.
- **Verify the sign on every path.** `flux%drainage` can now be negative on all three schemes. The
  bottom-face enthalpy debit (`meds_fast_split.f90:922`, `meds_fast_ark.f90:183`,
  `meds_fast_time_derivs.f90:338`) becomes a *credit* at `u_liq(T_n)` when it is. A sign-blind term
  will not fail loudly — it will quietly heat or cool layer *n*.

**Test.** Aquifer cells in `test_column_ark` and `test_column_rk45`: whole-column water and energy
close to machine precision with the bottom flux in **both** directions; the three schemes' bottom-layer
θ agree within a few percent over a wet month; `theta_ood_max` and `clamp_mass` stay at zero; and
RK45's sub-step count and `work_rk45_rescue_site` are recorded, not just checked. Run at least one cell
on the 29 mm h⁻¹ sealed fixture as well as on forced Ithaca — §9 of `numerical_scheme.md` is explicit
that the Ithaca cells never visit this corner.

### Phase 4 — row 4a: prognostic WOOD on ARK and RK45 (operator-split)

*Independent of Phases 0–3. This is `MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §5 P2's deferred half.*

Wood is operator-split out of both tableaux, exactly as plant hydraulics is. **The reason given in
`MEDS_LEAF_WOOD_ENERGY_DESIGN.md` §3 — "wood is non-stiff, `τ_wood` ≫ `dt_fast`" — is not the real
one and should not be relied on.** Reading the terms that kernel actually uses:

- `cap_wood = dbio_w·c_sapw + wmass_w·cp_liq` with `c_sapw = 2700` (`meds_biophysics_types.f90:249`),
  so at ~0.5 moisture fraction `cap_wood ≈ 4800·dbio_w` [J/m²/K];
- `|drdt_wood| ≈ π·WAI·gbh·ρcp + 4εσT³·WAI` (`meds_vegetation_biophysics.f90:150`, `:176`), convection
  dominating, `≈ 118·WAI` at `gbh = 0.03` m/s.

So `τ_wood ≈ 41·(dbio_w/WAI)` — sapwood mass per unit wood area, i.e. **a length scale**. For a
cylinder at `ρ_wood ≈ 500 kg/m³` that is `τ_wood ≈ 10⁴·r` seconds with `r` the stem radius in metres:
~50 s for a 1 cm sapling, ~500 s at 10 cm dbh, ~2500 s at 50 cm dbh. The crossover against
`dt_fast = 1800 s` sits near 35 cm dbh, and **most cohorts in a stand are below it**. §3's own range
(`cap_wood ~ 1e3–1e5` → "`τ` ~ minutes–hours") already contains the contradiction: minutes is not
≫ 1800 s. It read the top of its own range. A tableau is governed by the worst cohort in the patch,
not the mean.

The three reasons that do hold:

1. **`τ_wood` spans ~10¹–10⁴ s within one patch**, so an explicit in-tableau wood would set the whole
   column's step from the smallest sapling. RK45 averages ~3.4 sub-steps per `dt_fast` today (5033/month
   at 1800 s); a 50 s mode forces ≳36 on stability alone — roughly 10× the cost, in exactly the regime
   RK45 already rails and rescues in, to resolve a store whose entire point is that it lags slowly.
2. **The ARK has no explicit tableau branch to use.** `f_E == 0` and the scheme reduces to a 2-solve
   ESDIRK2 (`meds_fast_ark.f90:402-404`); `ark2_column_step` evaluates no explicit stage. "In the ARK
   tableau" therefore means "in the implicit block" — extending `newton_surface_solve`'s 2×2 arrowhead
   per cohort, which is Phase 5's work. Reviving `f_E ≠ 0` would not help, because by (1) wood belongs
   on the implicit side for most cohorts anyway.
3. **`split` operator-splits wood too** (`meds_fast_split.f90:624-639`: one BE step per cohort after the
   CAS Picard commit). Since parity with `split` is the goal, in-tableau wood on the adaptive schemes
   would close one cross-scheme difference by opening a new coupling-order one.

**The costs, stated rather than hidden.** The split is Lie–Trotter, so wood↔CAS coupling is
**first-order in `dt_fast`** — the lowest-order term in an otherwise second-order scheme. Two things
bound it: `split` carries the identical error, so no family gap opens; and the flux wood returns to the
CAS evolves on `τ_wood`, so the coupling error scales like `dt²/τ` rather than `dt`. Wood also sits
outside error control, mitigated only by `veg_energy_step_implicit` being L-stable — the failure mode
is lag, never blow-up. **Strang splitting would recover second order for one extra half-step** and is
the obvious refinement if the measurement below says the coupling error matters.

**Measure the premise before building on it.** The `τ_wood ≈ 10⁴·r` estimate above is
order-of-magnitude — `gbh` is assumed and the WAI↔stem-surface geometry is approximate. Emit the
per-cohort `cap_wood/|drdt|` distribution for a real census stand as the phase's first step, and record
it. If `min(τ_wood)` over a patch turns out to sit *above* `dt_fast` after all, reason (1) evaporates
and in-tableau wood on RK45 becomes the better option — cheaper to write and error-controlled. Do not
inherit this conclusion the way this plan inherited the last one.

Given the above:

- `advance_wood_energy_full`, mirroring the existing `advance_hydraulics_full` /
  `advance_water_mass_full` in `meds_fast_ark.f90`: seed the store from `bio%wood_temp` via
  `tissue_internal_energy`, take one `veg_energy_step_implicit(is_leaf = .false.)`
  (`meds_vegetation_biophysics.f90:191`) over the full `dt_fast` at the **committed** CAS endpoint,
  write `bio%wood_temp` back. `wood_temp` stays out of `column_state_t` and is warm-restarted, per
  that design's §4.
- RK45 calls the identical routine — it already shares `build_column_frozen`, and keeping the two
  adaptive schemes in one family is the whole point of `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §1.
- The wood storage delta joins the whole-energy ledger on both paths, or a slice of `coh_rnet` that
  used to pass through to the CAS now vanishes into a store the budget cannot see.
- Drop the guard at `meds_fast_split.f90:226`.

**A note that must go in the code, not just here.** Phase F-2 removed the `with_mass` / `with_theta`
norm opt-outs on the principle that "do not measure this state" is indistinguishable at the call site
from "this state cannot move". Excluding `wood_temp` from the WRMS norm is **not** a re-litigation of
that: wood is genuinely *outside* the tableau, like ψ, not inside-but-unmeasured. Write that
distinction where the exclusion happens so the next reader does not undo it — and write the
`τ_wood ≈ 10⁴·r` scaling next to it, since "wood is non-stiff" is the wrong reason and has already
propagated once.

**Test.** Prognostic-wood cells in `test_column_ark` / `test_column_rk45`: energy closes; wood
temperature demonstrably lags the CAS through a diel cycle; and the `cap_wood → 0` limit reproduces
the diagnostic-wood run bit-for-bit (the same limit check P3 used for the leaf on the split path).

### Phase 5 — row 4b: prognostic LEAF on ARK and RK45

*The hardest phase, and two genuinely different jobs. Depends on Phase 4 for the kernel.*

The leaf is stiff (`τ_leaf` ~ 1–100 s ≪ `dt_fast`), which is why it defaults to diagnostic and why
`split` requires `integration_scheme = "picard"` for it (`meds_fast_split.f90:312`).

- **ARK — the bordered arrowhead.** `surface_derivs` takes `T_leaf` as a *stage input* when the leaf is
  prognostic instead of diagnosing it, returning H/LE at that temperature.
  `newton_surface_solve` (`meds_fast_ark.f90:244`) and `jac_surface` (`:275`) extend from the current
  2×2 to a bordered arrowhead: each cohort's leaf is a 1-DOF unknown coupled **only** to the two CAS
  scalars — no leaf-leaf coupling — so per-cohort Schur elimination folds every leaf into a correction
  on the same 2×2, then back-substitutes the leaf temperatures. Cost stays O(ncoh) per Newton
  iteration.
  **Recommendation: keep leaf energy out of `column_state_t` for the MVP**, warm-restarted through
  `bio%leaf_temp`. `split`'s leaf rides the Picard iterate and is not error-controlled either, so
  adding it to the controlled state vector would make ARK *stricter* than the path we are matching —
  a new asymmetry in the name of removing one. Revisit only if leaf freeze/thaw is ever wanted.
- **RK45 — in the tableau, and expect it to hurt.** Architecturally the simpler of the two: add
  `leaf_energy(:)` to `column_state_t` plus its `state_axpy` / `state_accum` / `state_wrms` /
  `state_err_diff` entries and a `GRP_LEAF_E` tolerance group in `meds_fast_control.f90` (which
  currently has 8 groups, `GRP_ENTH`…`GRP_SOIL_T`). This is the ED2-faithful configuration — ED2's own
  `rk4` carries leaf internal energy explicitly — and it is exactly the stiff mode RK45 already rails
  on and rescues from. **Measure `work_rk45_rescue_site` before calling it supported.** If a
  vegetated winter month rescues most steps, the honest deliverable is a documented
  "supported, not recommended, here is the rescue rate" tier, not a quiet fallback. Do not add a new
  fallback path to make the number look better — §1 of the production plan is explicit that a
  cross-family rescue is a mid-run model switch.
- Drop the guard at `meds_fast_split.f90:228`.

**Test.** Prognostic-leaf cells on all three schemes: conservation telescopes; the ARK stays stable at
`dt_fast = 1800 s` where `meds_fast_rk4_oracle` blows up (the L-stability check P4's own scoping
names); `cap_leaf → 0` reproduces the diagnostic run bit-for-bit; and a three-scheme agreement cell
under `integration_scheme = "picard"` on the split side, since that is the only split configuration
prognostic leaf is legal in.

### Phase 6 — retire the `--parity` preset, and correct the docs

The completion criterion made mechanical: **`--parity` shrinks to empty as the plan lands.**

- `scripts/numerics_sweep.py:92-109`: add the missing `soil.zeng_decker` pin now (it is unpinned
  today), then delete each pin as its phase lands — `hydraulics.multilayer_roots` after Phase 1 (the
  key itself no longer exists),
  `soil.bottom_bc` and `soil.zeng_decker` after Phase 3, the two `*_energy_model` pins after Phase 5.
  When `PARITY` is `{}`, the Class-1 inventory is closed and the flag can go with it.
- Add the `parity_fidelity.py` record-count / duplicate-key guard that `MEDS_PRODUCTION_INTEGRATOR_PLAN.md`
  P0-d still lists as outstanding — every phase here is scored with that harness.
- `docs/science/numerical_scheme.md` §4: rewrite rows 3, 4 and 5 as they land.
- `MEDS_INTEGRATOR_PARITY.md`: rows 3/4/5 move off "C6 — deferred", and **row 5's "C5 DONE" framing is
  corrected** — C5 added a *guard* to RK45 so it failed the same way the ARK already did. It did not
  implement the boundary condition, and reading that row today suggests otherwise.

---

## 4. Acceptance criteria

The Class-1 inventory is closed when all of these hold:

1. Every phase green on ifx Release, ifx Debug (`-check all`) and nvfortran multicore. The nvfortran
   build is not optional here: several phases add array-valued expressions, and CLAUDE.md's
   array-temporary miscompilation trap is silent under a green ifx run.
2. **Every default-configuration run is bit-identical to `main`, except at the two phases that
   deliberately are not.** Phase 0 (aquifer-BC runs only) and Phase 1 (every run, since
   `multilayer_roots` becomes unconditional) each re-baseline once with a documented before/after.
   Everywhere else a moved golden anchor means a phase leaked.
3. For each of the three rows: a single-`dt_fast`-from-identical-state comparison shows `split`, `ark`
   and `rk45` producing the same physics to the level the row governs — per-layer sink profile for
   row 3 (now unconditional), store temperature and its ledger contribution for row 4 (feature on),
   bottom-face mass and enthalpy for row 5 (aquifer BC selected).
4. Whole-column water and energy ledgers close with `[energy].debug_error` armed, on all three
   schemes, in the regime each feature is *for* — a dry-down for row 3, a diel cycle for row 4, a
   saturating and a draining column for row 5. Verify each guard is genuinely armed by tightening its
   tolerance until it fires, rather than inferring closure from a clean exit (P0-e's lesson).
5. At least one cell per row runs on a **sealed unit fixture** at the corner, not only on forced
   Ithaca. Two of the defects in §7 of `numerical_scheme.md` were visible only at 29 mm h⁻¹ and are
   measurable no-ops in every Ithaca cell.
6. `PARITY` in `numerics_sweep.py` is empty, and `numerical_scheme.md` §4 has no `physics` rows left.

---

## 5. Risks and things deliberately not in this plan

- **Phase 5's RK45 half may not be usable.** A prognostic leaf inside an explicit tableau is the stiff
  mode that motivated the rescue path in the first place. The plan's deliverable is a *measurement*
  and an honest tier label, not a promise that it is fast. Do not let a bad rescue rate be papered
  over by a new fallback.
- **Phase 0 deletes a subsystem and changes a restart format.** `w_aquifer`, `z_wt`, `SPECIFIC_YIELD`,
  the baseflow bucket, Dunne `f_sat`, and the five parameters `f_max` / `f_over` / `q_drai_max` /
  `f_drai` all go. Two consequences to handle deliberately rather than discover: the `-S-` state file
  loses two variables, so decide the older-file read path explicitly; and `flux%drainage` can now be
  **negative** (inflow), which every consumer and every `[output]` diagnostic must tolerate. A
  sign-blind consumer will not fail loudly.
- **`aquifer` stops being a drop-in alternative to `free_drain`.** With the water table defined as the
  column base, an aquifer run is a permanently wet site with an unlimited supply from below. That is
  right for riparian, floodplain and wetland stands and wrong for a typical upland one. It needs to say
  so where the option is documented, or someone will select it as "the more realistic drainage".
- **The residual `split` ↔ ARK/RK45 gap (F3, ~0.45 K) is out of scope.** It is row 6, not rows 3–5,
  and it is `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` P1-a. Nothing in this plan should be expected to
  move it, and if something here *does* move it, that is a finding worth chasing rather than a win.
- **Folding soil water into the ARK tableau** stays out, for the reasons the production plan gives —
  it is P1-a's leading candidate and its scoping depends on that answer.
- **Hydraulic redistribution** stays disabled. Row 3 makes per-layer sink *placement* uniform across
  schemes; it does not enable root efflux, which is floored to zero in both the plant solver and the
  soil sink project-wide.
- **The soil column's bottom THERMAL boundary is out of scope, and it has a real problem.** Filed
  separately rather than folded in, because it is identical on all three schemes and mixing a
  model-depth question into a scheme-comparison question is exactly the confound this whole line of
  work exists to remove. Recorded here so it is not lost:

  `ENERGY_BC_GEOTHERMAL` is hardwired to `0.0` on all three paths (`meds_fast_split.f90:894`,
  `meds_fast_ark.f90:1488`) with **no TOML reader** — an insulated base nobody chose. Meanwhile
  `build_soil_hydr_params` is called with `soil_depth = 2.0` m, also hardwired and also with no TOML
  reader (`meds_fast_dynamics.f90:99`). The annual thermal damping depth for moist soil is
  `d = √(2κ/ρcω) ≈ 2.5 m`, so at the column base the **annual wave is still at ~44% of its surface
  amplitude, hitting a zero-flux wall**. The column is shallower than its own thermal skin depth. The
  diurnal wave is fully damped (`d_daily ≈ 0.13 m`), which is why no month-long fast-loop test can see
  this — it needs a multi-year run, i.e. the 30-yr Ithaca context.

  **Fix by deepening the column (8–10 m) and giving `soil_depth` a TOML reader**, or by a Dirichlet
  deep-temperature BC at the local annual mean if depth is too costly. Do **not** "continue the
  gradient" into a ghost cell: `T(n+1) = 2T(n) − T(n−1)` makes the bottom face flux exactly equal to
  the top face flux of layer *n*, so its divergence is identically zero and `dT(n)/dt = 0` — the
  deepest layer would be pinned at its initial temperature forever. That is a worse artefact than the
  reflection it is meant to cure.

  **The `bedrock` BC is separately inconsistent between its two fluxes.** It blocks water and heat
  equally, but rock conducts heat *better* than soil (κ ≈ 2–3 vs ~1.5 W/m/K). A bedrock base should be
  hydraulically sealed and thermally open; today the hydraulic BC and the thermal BC disagree about what
  lies beneath the column. Cheap to fix once `geothermal` has a config reader — and it belongs with that
  work, not with parity.

  Note the shape: a physically meaningful parameter, silently pinned at a default nobody selected,
  because it has no config reader. That is the same class as `[energy].debug_error`
  (`MEDS_PRODUCTION_INTEGRATOR_PLAN.md` P0-e). Worth a sweep for others.

  **A geothermal option is deferred (decision 5, 2026-07-30) — `geothermal = 0` stays on every BC.**
  When it is taken up, it must be taken up *together* with Phase 0's `T_aq = T_n`, because the two are
  the same assumption seen twice: `geothermal = 0` states `dT/dz = 0` at the base, and `T_aq = T_n` is
  its advective counterpart. Enabling a nonzero bottom heat flux without giving the inflowing water a
  matching temperature would leave the conductive and advective boundary treatments describing
  different bodies. The clean form when it is done is a single prescribed deep-boundary temperature
  `T_deep` driving both — a restoring BC, which also fixes the wave reflection above, rather than a
  prescribed flux. **Do not** extrapolate the interior gradient into a ghost cell (see the paragraph
  above for why that pins the deepest layer).


---

## 7. Implementation status (2026-07-30)

### Landed on `feature/integrator-physics-parity`

**Unplanned — soil ψ units at the hydraulics seam.** Found while preparing Phase 1.
`soil_psi_from_theta` returns **metres** of head (`meds_soil_water` adds `dz_node` to it directly and
exports `grav_head *` it), but `meds_fast_split.f90:668` / `meds_fast_ark.f90:1519` passed the raw
result into `hydro_env_t%soil_psi`, documented **MPa**. Both call sites replaced the older
`grav_head`-converted `hflux%psi_soil` when hydraulics moved ahead of the soil solve, and the
conversion went with it — so the plant solver was handed a potential **~102× too negative**. Inert in a
wet column (it stays above `wstress_psi_open = −0.5 MPa`), severe on drying, which is why no test
caught it and why a dry-down window would have. It had to land before Phase 1, because the multilayer
path adds a correct-MPa `grav_head·z_k` to an inflated `ψ_k` — inconsistent within one formula.

Measured on `test_column_dynamics` RUN 1 at noon: `ψ_leaf` −1.532 → −0.805 MPa, GPP 17.9 → 14.1,
`T_leaf` 306.4 → 295.7 K, `T_CAS` 292.7 → 304.4 K. Mechanism: with a spuriously dry soil the plant
could not take up water (uptake floored at 0), so the column stayed wet, ground evaporation held the
CAS down and the leaf ran hot because it could not transpire. GPP falls because the leaf cools 10 K
and the Arrhenius loss in Vcmax outweighs the relief in water stress.

**Phase 1 — DONE.** `multilayer_roots` deleted; per-layer coupling unconditional; current-step shares
on all three schemes via `column_frozen_t%root_share`; `root_sink_share`, `rhizo_cond` and the
single-BC batch arguments removed; `K(θ)` hoisted out of the cohort loop. New dry-down coverage in
`test_plant_hydraulics`. `test_picard_coupling` goldens rebased once with a documented before/after
(CAS +10.661 K, soil-surface −3.099 K, covering both this and the units fix).

**Phase 2 — DONE.** `column_hydrology_flux` exports `psi_e`; `build_column_frozen` copies it instead
of hardwiring zero (fixes V4 — RK45 silently ignored `zeng_decker`); both guards narrowed to
`bottom_bc == SOIL_BC_AQUIFER`. New `test_rk45_bedrock_and_zd`, whose ZD assertion is the mutation
proof that `psi_e` reaches the tendency.

37/37 on ifx Release and nvfortran multicore at each commit.

### NOT implemented

**Phase 0 (aquifer BC rebuild) and Phase 3 (unlock it on ARK/RK45).** Not started. The guards now
name the aquifer BC specifically, so the gap is non-silent and correctly scoped, but every item in
Phase 0 above — the head-driven flux, deleting the storage bucket / baseflow / Dunne `f_sat` and their
five parameters, the `recharge` vs `baseflow` split, `T_aq = T_n` inflow enthalpy — remains to do.

**Phase 4 (prognostic wood on ARK/RK45).** Not started. Note its first step is a *measurement*, not
code: emit the per-cohort `cap_wood/|drdt|` distribution before building on the `τ_wood ≈ 10⁴·r`
estimate.

**Phase 5 (prognostic leaf).** Not started. The largest phase by a wide margin — the ARK bordered
arrowhead with per-cohort Schur elimination, plus RK45's in-tableau leaf energy with its
`column_state_t` / `state_*` / tolerance-group expansion.

**Phase 6 — partial.** `--parity` lost its `multilayer_roots` pin (the key no longer exists) and
gained a note on `zeng_decker`; `soil.bottom_bc` is still pinned, correctly, until Phase 0/3.
`numerical_scheme.md` §4 updated for the two rows that closed. The `parity_fidelity.py`
record-count guard (P0-d's remaining half) is still outstanding.
