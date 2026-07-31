# MEDS design: multi-layer root water uptake (turn on `rhizosphere_cond`)

**Status:** design-only (no code). **Date:** 2026-07-18.

> **UPDATE 2026-07-30 — the `[hydraulics].multilayer_roots` opt-in flag is slated for DELETION and
> this path becomes unconditional.** Decided in `MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md` (decision 4,
> Phase 1), on the reasoning that retired `[fast].snow_on` and `with_mass`/`with_theta` before it: a
> flag whose OFF path is the cruder physics gets plumbed through every call site forever and
> eventually reads as "this cannot happen". Consequences for this document: the "bit-identical when
> off" property below stops being a design goal, the single-BC `(soil_psi, rhizo_cond)` path is
> removed from the *driver* (the per-cohort kernel keeps its `n_root_layer <= 1` branch, because
> `src/plant` must stand alone), and the default answer changes for every run. Hydraulic
> redistribution stays disabled — the non-negative floors on both sides are unaffected.
**Scope:** replace the single prescribed rhizosphere BC `(soil_psi, rhizo_cond)` with an ED2-faithful
**multi-layer root water uptake** — each soil layer the roots reach contributes to root uptake,
weighted by a root-distribution profile and the layer's soil water potential/conductivity — and
promote `rhizosphere_cond` from a test helper to the production per-layer conductance kernel.
Companion deliverable: a science reference at `docs/science/plant_hydraulics.md` (equations + params).

This realizes the *reduced* form of §16 of `MEDS_HYDRAULICS_DESIGN.md` ("per-soil-layer fine-root
water nodes"); §16's full per-layer **nodes** are the later Phase B here.

---

## 1. Motivation & the key reduction

Today `solve_plant_water` sees one scalar soil boundary: `env%soil_psi`, `env%rhizo_cond` (a prescribed
`5e-4`, currently from `cfg%hydraulics%rhizo_cond`). The wood-node soil term is
`rhizo_cond·(soil_psi − ψ_W)`. Roots don't "see" the soil profile: a wet deep layer and a dry shallow
layer are indistinguishable.

**The reduction that makes this cheap:** every soil layer connects to the *same* wood/collar node, so
the per-layer root+rhizosphere conductances are in **parallel to a common node**. A parallel network to
one node collapses exactly to an effective pair:

$$G_{\text{root}} = \sum_k g_k, \qquad
  \psi_{\text{soil,eff}} = \frac{\sum_k g_k\,(\psi_{\text{soil},k} + \rho g\,z_k)}{\sum_k g_k}$$

Substitute `(G_root, ψ_soil,eff)` for the scalar `(rhizo_cond, soil_psi)` in `freeze_coeffs` and **the
2×2 matrix-exponential solver is unchanged** — this is *pre-solve input aggregation*, not a solver
rewrite. After the solve, distribute the converged total uptake back to layers for the soil-water sink.

**This is exactly what ED2 does.** ED2's `plant_hydro.f90` computes `weighted_gw_cond = Σ_k gw_cond(k)`
and `weighted_soil_psi = Σ_k gw_cond(k)·soil_psi(k)`, advances the wood node with those, then sets
`wflux_gw_layer(k) = layer_water_supply(k)/Σ·wflux_gw`. So "multi-layer roots like ED2 hydro" **is** the
aggregation approach — not per-layer root nodes (that is FATES-HYDRO / Phase B). Better still, MEDS's
existing `rhizosphere_cond` **is** ED2's `gw_cond`:

```
ED2:  gw_cond(k) = soil_cond(k)·√RAI / (π·dslz(k)),     RAI = broot·SRA·root_frac(k)·nplant
MEDS: rhizosphere_cond = soil_cond·√RAI / (π·dz) / nplant   (per-plant; RAI same)
```

So the per-layer kernel already exists and is validated — it is only unused ("TEST helper"). Phase A is
mostly wiring + aggregation.

## 2. Root distribution profile (the one genuinely new science)

Per-layer root fraction from ED2's cumulative-exponential `root_beta` profile:

$$\text{root\_frac}(k) = \beta^{\,d_{k-1}/D} - \beta^{\,d_k/D}$$

with `β = root_beta` (per-PFT, ∈(0,1)), `d_k` the depth of layer-`k`'s bottom, and `D = −slz(krdepth)`
the maximum rooting depth. Roots exist only in layers `k ≥ krdepth`. New per-PFT traits: `root_beta`,
`root_depth` (→ `krdepth`), and `specific_root_area` (SRA); `broot` (fine-root C) comes from the carbon
state. These join `hydraulics_config_t` / the `[hydraulics]` block (and later become per-PFT).

## 3. The gate: per-layer soil state

The aggregation needs `ψ_soil,k` and `soil_cond(k)` **per layer** — real soil-column state, not one
prescribed scalar. That state lives in the biophysics soil-hydrology column
(`meds_column_hydrology` / `meds_soil_parameters`, which already compute per-layer `ψ(θ)` and `K(θ)`),
but it is **not yet coupled to plant hydraulics** — this is the deferred `soil_water_coupling` (P3) item.
So Phase A has a real dependency: either (a) land the soil→hydraulics `ψ_soil(:)` read first, or (b)
ship Phase A against a *prescribed per-layer* `soil_psi(:)`/`soil_cond(:)` profile (a vector
generalization of today's scalar BC) so the uptake machinery is exercised and tested ahead of the
coupling. **Recommend (b)** — it decouples this feature from the soil-coupling schedule and keeps each PR
small; the coupling then just fills the vectors from the real column.

## 4. Interface & code changes (Phase A)

| Piece | Change |
|---|---|
| `hydro_env_t` | `soil_psi` scalar → `soil_psi(NL)`; add `soil_cond(NL)`, `root_frac(NL)`, `layer_z(NL)`, `n_root_layer`. (Fixed `NL` = compile-time max soil layers, GPU/SoA-safe.) |
| `hydro_flux_t` | add `root_uptake_layer(NL)` (per-layer sink for the soil module); keep scalar `root_uptake` = Σ. |
| `rhizosphere_cond` | promote from "TEST helper" to production; call per layer in the new aggregator. |
| aggregator (new, in solver) | `aggregate_root_boundary(env) → (G_root, ψ_soil,eff)`; called once per `freeze_coeffs` (or once per sub-step, since `g_k` depends on `ψ_soil,k`, `soil_cond(k)`, not on the prognostic `ψ`). |
| `freeze_coeffs` | replace `rhz = env%rhizo_cond`, `soil_psi = env%soil_psi` with the aggregated pair. **Everything else identical** — the 2×2 assembly, the matrix-exp, the adaptive stepping are untouched. |
| post-solve | distribute total uptake: `U_k = supply_k / Σ supply · U`, `supply_k = g_k(ψ_soil,k − ψ_W + ρg z_k)` (ED2's exact-conservation split). |
| config | `root_beta`, `root_depth`, `specific_root_area` added to `hydraulics_config_t` + `[hydraulics]` TOML (defaulted; opt-in). |
| back-compat | `n_root_layer = 1` with the old scalar reproduces today's single-BC result **bit-identically** (`G_root = rhizo_cond`, `ψ_soil,eff = soil_psi`) — the regression guard. |

Numerically the wood ODE is unchanged in form:
`C_W ψ'_W = K_LW(ψ_L − ψ_W − g_WL) + G_root(ψ_soil,eff − ψ_W)`.

## 5. Phasing

- **Phase A — ED2-faithful aggregation (this design).** Vector soil BC + `root_beta` profile +
  per-layer `rhizosphere_cond` + parallel aggregation + per-layer uptake distribution. 2-node solver
  unchanged. **Hydraulic redistribution (HR) is intentionally NOT enabled** (deferred to a future
  version): the per-layer distribution floors each layer's supply at 0
  (`max(g_k(ψ_soil,k − ψ_W + ρg z_k), 0)`), so a dry layer that would efflux (`U_k < 0`) gets 0 — every
  `U_k ≥ 0` and `Σ U_k = Q_root`. HR (the root→soil efflux and the soil re-wetting it implies) is the
  clearest thing Phase B / a future pass adds.
- **Phase B — per-layer root nodes (§16, later).** A fine-root node `R_k` per layer with its own
  capacitance; `Leaf–Collar–{R_k}` arrow/bordered matrix; `N_HYDRO` bump; `BE` (Thomas/Schur) solver.
  Adds root storage + a true per-layer root potential (fuller HR). Phase A is a strict subset and a
  clean stepping stone — no rework, only added nodes/edges (§16 confirms compatibility point by point).
- **Coupling — DONE (soil→plant, opt-in `[hydraulics].multilayer_roots`, default off ⇒ bit-identical).**
  `column_fast_step` fills the multi-layer boundary from the prognostic soil column: `soil_psi_layer(k)
  = hflux%psi_soil(k)` (already [MPa], "EXPORTED to hydraulics"), `rhizo_cond_layer(k) =
  rhizosphere_cond(ρ·ksat(k)/grav_head, broot, SRA, root_frac(k), dz(k), nplant)`, `root_z_layer(k) =
  z_node(k)`. The plant then aggregates by conductance and returns `flux%root_uptake_layer`.
  `test_column_dynamics` RUN 3 asserts conservation still holds and the plant ψ shifts vs the single
  root-fraction-weighted BC.
  - **Plant→soil feedback — DONE (budget-safe).** The soil sink distributes the SAME total
    (`coh_transp`) by the previous fast step's normalized per-layer uptake shares
    (`bio%root_sink_share`, default 0 ⇒ root_frac fallback) instead of the static `root_frac` — so the
    soil dries where roots actually took water. Both share sets sum to 1, so the column-total water
    balance (and every budget) is unchanged; only the vertical distribution differs. The plant solver
    already floors per-layer efflux to 0 (HR disabled, above), so all shares are ≥ 0; shares are lagged
    one sub-step to avoid a plant↔soil ordering cycle.
  - **Unsaturated `K(θ)` — DONE.** The per-layer rhizosphere conductance now uses
    `soil_hydr_cond(θ)` (Mualem–van Genuchten / Campbell) instead of `ksat`, so dry layers are
    conductance-down-weighted, not just ψ-down-weighted.

## 6. Deliverables

1. **`docs/science/plant_hydraulics.md`** (this PR) — science reference: PV/capacitance curves,
   vulnerability + Kirchhoff conductance, the node network, the multi-layer root uptake (this design),
   the solver, and the parameter table. Matches `docs/science/leaf_gas_exchange.md` style.
2. This design doc.
3. (Implementation PR) Phase A per §4, with a `test_plant_hydraulics` case: an `n_root_layer > 1`
   profile whose aggregate reproduces the single-BC result, plus a wet/dry case asserting the drier
   layer is floored to 0 (no efflux; HR disabled) and the supplying layer takes all the uptake.

## 7. Validation

- **Back-compat:** `n_root_layer = 1` ⇒ bit-identical to today (the aggregate pair equals the scalar BC).
- **ED2 cross-check:** for a uniform-ψ profile, `Σ U_k = G_root(ψ̄ − ψ_W)`; `Σ U_k = U` exactly (the
  proportional split conserves).
- **No efflux (HR disabled):** wet/dry ⇒ the supplying layer takes all the uptake and the drier layer
  is floored to `U_k = 0` (never negative), `Σ U_k = Q_root`.
- Both compilers (ifx `-check all` + nvfortran multicore), per the repo's portability rule.

## 8. Risks / open questions

- **Fixed `NL`.** A compile-time max soil-layer count keeps `hydro_env_t`/`hydro_flux_t` fixed-shape
  (GPU-safe), mirroring `N_HYDRO`. Must match the soil column's `n_soil_layer_max`.
- **`g_k` recompute cadence.** `g_k` depends on `ψ_soil,k`/`soil_cond(k)` (frozen over `dt`), not the
  prognostic `ψ` — so aggregate **once** per `solve_plant_water` call, not per sub-step (cheaper, and
  matches the frozen-BC contract). Confirm with the adaptive path.
- **Units.** MPa throughout (per the module contract); `rhizosphere_cond`'s `soil_cond` is
  `[kg m⁻¹ s⁻¹ MPa⁻¹]` — the soil module must supply `K` in MPa-gradient units (a conversion at the seam).
- **Per-PFT vs uniform.** Ship uniform (matching current `hydraulics_config_t`); per-PFT root traits
  fold in when hydraulics becomes per-PFT.

## References
- ED2 `ED/src/dynamics/plant_hydro.f90` (`weighted_gw_cond`, `weighted_soil_psi`, `root_frac` via
  `root_beta`, `layer_water_supply`, `wflux_gw_layer`); `MEDS_HYDRAULICS_DESIGN.md` §4, §16.
- Katul et al. 2003 (rhizosphere conductance); Xu et al. 2016 (X16, the physical reference).
