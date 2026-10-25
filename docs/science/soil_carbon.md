# Soil carbon

The decomposable carbon of one patch is a seven-element vector $X$ advanced **once per day** by a
CENTURY-family **carbon matrix ODE**. It is a port of ED2's `update_C_and_N_pools` — the same pools,
decay rates, respired and transferred fractions — **reorganized as one linear-algebra object** instead
of a cascade of hand-written pool-to-pool transfers. That is not cosmetic: with the decay diagonal $K$,
the donor/transfer matrix $A$ and the environmental scalar $\xi$ assembled as matrices, the exact
matrix-exponential step, the semi-analytic steady state (SASU) and the traceability diagnostics all
fall out of the *same* operator instead of being written a second time.

## 1. The carbon matrix ODE

Per patch, with $X$ in $`[\mathrm{kgC\,m^{-2}}]`$, and element-wise the mass balance it encodes:

```math
\frac{dX}{dt} = B\,I(t) + A\,\xi(t)\,K\,X(t) \quad [\mathrm{kgC\,m^{-2}\,day^{-1}}],
\qquad
\frac{dX_i}{dt} = u_i + \sum_{j\neq i} a_{ij}\,\xi_j K_j X_j - \xi_i K_i X_i \qquad(1)
```

$`u = B\,I`$ is the litter input (§5); $K$ the diagonal of baseline decay rates, carried in
$`[\mathrm{yr^{-1}}]`$ by the config and divided by `yr_day` at assembly so the kernels run in
$`[\mathrm{day^{-1}}]`$; $\xi$ the diagonal per-pool environmental scalar (§4); $A$ the dimensionless
donor/transfer matrix (§3). The input is **litter carbon, never GPP** — autotrophic respiration is
taken out upstream, in the plant carbon budget.

## 2. The seven pools and the lignin sub-tracer

The state is `soil_carbon_state_t`, ordered litter → SOM → passive. The seven-pool ceiling matches ED2's
CENTURY set exactly, and the litter half is a **metabolic/structural × above/below** cross — chemistry
on one axis, position on the other (above/below is carried so a future fire class can burn the top).

| `IP_*` | pool | field | $K$ default $`[\mathrm{yr^{-1}}]`$ | gets litter | $\xi$ temperature |
|---|---|---|---|---|---|
| `IP_FAST_GRND` 1 | metabolic litter, above | `fast_grnd_carbon` | `k_fast` 11.0 | yes | surface |
| `IP_FAST_SOIL` 2 | metabolic litter, below | `fast_soil_carbon` | `k_fast` 11.0 | yes | active layer |
| `IP_STRUCT_GRND` 3 | structural litter + CWD, above | `struct_grnd_carbon` | `k_struct` 4.5 | yes | surface |
| `IP_STRUCT_SOIL` 4 | structural litter + CWD, below | `struct_soil_carbon` | `k_struct` 4.5 | yes | active layer |
| `IP_MICR` 5 | microbial SOM | `microbial_carbon` | `k_micr` 0.0 | no | active layer |
| `IP_SLOW` 6 | slow / humified SOM | `slow_carbon` | `k_slow` 0.2 | no | active layer |
| `IP_PASSIVE` 7 | passive SOM | `passive_carbon` | `k_pass` 0.0 | no | active layer |

A pool is **active** when $`K_j>0`$. Under the default `DECOMP_SCHEME_ED2` (ED2 schemes 0–4) microbial
and passive are inert, leaving **three active chemistries** — fast, structural, slow;
`DECOMP_SCHEME_CENTURY5` activates all seven for the **five-active** CENTURY topology. The selector
switches topology and temperature form, **not** the rate constants, which come from `[soil_carbon]` and
default to the scheme 0–4 values. Dead wood has no pool of its own: it enters the structural pools
through `f_labile_stem`, ED2's convention.

The structural pools carry a **lignin sub-tracer** `struct_grnd_lignin` / `struct_soil_lignin`, bounded
$`0\le L_s\le C_s`$. It is not an eighth pool — its carbon is already counted inside the structural
pools, so it stays outside the mass vector and keeps its own balance (§6). There is no nitrogen
(§10). `pack_pool_vector` / `unpack_pool_vector` marshal the named fields to and from the length-7
array the kernels use.

## 3. The transfer matrix, and the respired complement

`decomposition_transfer_matrix` builds $A$, $K$ and the respired fractions together. Each donor loses *all*
of its decayed carbon, so the diagonal is uniform; Rh is exactly what $A$ does *not* transfer:

```math
a_{jj} = -1 \ \ \forall j, \qquad e_{r,j} = 1 - \sum_{i\neq j} a_{ij} = -\sum_i a_{ij},
\qquad R_h = \sum_j e_{r,j}\,\xi_j K_j X_j = -\mathbf{1}^{\mathsf{T}} A\,\xi\,K\,X \qquad(2)
```

with $`a_{ij}\in[0,1]`$ the fraction of $j$'s loss entering pool $i$ and $`e_{r,j}`$ — the **respired
complement** — the fraction leaving the network as CO₂. Keeping $`a_{jj}=-1`$ on inert pools too is
deliberate: $`K_j=0`$ zeroes them anyway, and the column identity $`\sum_i a_{ij}=-e_{r,j}`$ then holds
for **every** pool, which the tests check.

**Default scheme.** Metabolic litter respires 100% and transfers nothing; structural sends its
complement to slow ($`a_{6,3}=1-e_{r,3}`$, $`a_{6,4}=1-e_{r,4}`$); slow is the sole SOM outlet and
respires 100%. The structural respired fraction is lignin-weighted,
$`e_{r,s}=(1-f_{\ell,s})r_{\mathrm{nonlig}}+f_{\ell,s}r_{\mathrm{lig}}`$ with $`f_{\ell,s}=L_s/C_s`$
(0 when $`C_s=0`$), so at the default 0.3 structural sends 70% to slow.

**Scheme 5** (Bolker/Pacala/Parton CENTURY). Fast litter → microbial; structural non-lignin →
microbial, lignin → slow; microbial ⇄ slow and slow ⇄ passive back-transfers, **texture-controlled**:
sand sets the microbial respired fraction
$`e_{r,\mathrm{micr}}=\texttt{er\_micr\_int}+\texttt{er\_micr\_slp}\cdot x_{\mathrm{sand}}`$, clay the
fractions routed to passive, $`\min(1-e_r,\ \texttt{int}+\texttt{slp}\cdot x_{\mathrm{clay}})`$ — and
those back-transfers sit above the diagonal, so scheme 5's active block is non-triangular (§8).

## 4. The environmental scalar, and where it sits in the product

`decomposition_env_scalar` returns one scalar per pool: temperature × moisture/oxygen × (structural only) a
lignin brake, over the relative saturation $`m`$ clipped to $`[0,1]`$.

```math
f_T = \begin{cases} q_{10}^{\,(T-T_{\mathrm{ref}})/10} & \text{scheme 5, UNCAPPED Q10}\\
\min\!\big(1,\ \exp[a_T(T-T_{\mathrm{sat}})]\big) & \text{schemes 0--4, ED2 capped exp}\end{cases}
\qquad
f_\theta = \begin{cases} \exp[(m-m_{\mathrm{opt}})\,b_{\mathrm{dry}}] & m \le m_{\mathrm{opt}}\\
\exp[(m_{\mathrm{opt}}-m)\,b_{\mathrm{wet}}] & m > m_{\mathrm{opt}}\end{cases}
\qquad(3)
```

```math
m = \frac{\theta-\theta_{\mathrm{dry}}}{\theta_{\mathrm{sat}}-\theta_{\mathrm{dry}}}, \qquad
\xi_j = f_T(T_j)\,f_\theta(\theta)\,\Lambda_j, \qquad
\Lambda_j = \begin{cases}\exp(-e_{\mathrm{lig}}f_{\ell,j})\cdot f_{\mathrm{decomp}} & j\in\{3,4\}\\
1 & \text{otherwise}\end{cases} \qquad(4)
```

One one-sided exponential covers both moisture limbs — dry-side substrate diffusion, wet-side anoxia —
and it is the **same functional form the fast Rh kernel uses**, so the two timescales reconcile at
matched inputs. $`T_j`$ is the **ground surface** temperature for the above-ground pools (1, 3) and the
**root-weighted soil** temperature for the rest (ED2's `A_decomp` / `B_decomp` split), $\theta$ the
depth-weighted column mean, $`\theta_{\mathrm{dry}}`$ and $`\theta_{\mathrm{sat}}`$ the top layer's residual and
saturated water contents, $`f_{\mathrm{decomp}}`$ the N immobilization brake (identically 1 while N
is off). The invariant is $`\xi>0`$, not $`\xi\le1`$: scheme 5's Q10 is uncapped.

**Placement is physics.** The scalar sits *inside* the transfer, $`A\,(\xi K X)`$, never
$`\xi\,(A K X)`$: **each donor decomposes at its own rate $`\xi_j K_j`$, and only then is the loss
split** into respired and transferred parts. $\xi$ and $K$ commute with each other but neither commutes
with $A$, so the receiver-scalar ordering silently breaks mass balance and every capacity/residence
formula. The tests assert it: the correct form's total loss equals $`R_h`$ of eq. (2); the other's does not.

## 5. Litter input

Dead plant carbon becomes litter in one place — `necromass_to_litter` (`meds_soil_carbon`), whose
arguments are plain carbon amounts and plain PFT traits, which is what lets both the slow driver and the
demographic operators call it. Each tissue goes where it is:

| Tissue | Labile share | Above / below ground |
|---|---|---|
| leaves; flowers, fruit and seed of reproduction | `f_labile_leaf` | all above |
| fine roots | `f_labile_leaf` | all below |
| wood (stem + coarse roots) | `f_labile_stem` | `aboveground_frac` / the rest |
| storage (non-structural carbon, in sapwood and coarse roots) | `f_labile_leaf` | `aboveground_frac` / the rest |

`aboveground_frac` is the woody split (ED2 `agf_bs`), so it applies to wood and storage only. Only the
structural streams carry lignin — `struct_lignin_frac` of each. The slow driver calls it through
`add_necromass`, which also records the carbon by the tissue it came from; those are the litter outputs
(`mean_litter_leaf_py`, `mean_litter_fineroot_py`, `mean_litter_wood_py`, `mean_litter_repro_py`, each after
resorption), and with `mean_root_exudate_py` they sum to the carbon entering pools 1–4.

`build_litter_input` maps those six numbers onto the matrix source term — four carbon streams into $u$,
two lignin fluxes alongside — and does **no** chemistry of its own, since the splits are per-cohort
decisions made before summing to the patch. So **only pools 1–4 ever receive litter**; microbial, slow
and passive are fed by transfers alone. Carbon arrives by **one route**: turnover and shedding,
background mortality, and the plants a calendar boundary's cull and treefall kill all accumulate into
`patch%litter_in` through `add_necromass`, and enter through $u$ at the next daily step. Once that step
has taken it up the accumulator is emptied, so a boundary's dead carbon reaches the pools with the
step after it.

## 6. The daily step, and why it consumes an accumulated integral

`soil_carbon_step` advances one patch by one day, and the decay it applies is **not** an instantaneous
$\xi$ — it is the fast loop's **accumulated per-pool loss integral** $`\xi^{\mathrm{int}}`$, with its
dimensionless decay fraction $`d`$. The default solver (`DECOMP_STEP_EULER`) is ED2-faithful forward
Euler on the frozen operator ($`a_{ii}=-1`$ removes donor $i$'s own loss, the off-diagonals transfer):

```math
\xi^{\mathrm{int}}_j = \int_{\mathrm{day}} \xi_j(t)\,dt \ [\mathrm{day}], \qquad
d_j = \xi^{\mathrm{int}}_j K_j \ [-], \qquad
X^{n+1}_i = X^n_i + u_i + \sum_j a_{ij}\,d_j\,X^n_j \qquad(5)
```

For **large accelerated steps** `DECOMP_STEP_EXPM` takes the *exact* exponential of that same frozen
operator, through an augmented matrix that carries the constant input without inverting $A$ — which
matters because $`A\,\mathrm{diag}(d)`$ is singular whenever a pool is inert:

```math
\begin{bmatrix} X^{n+1}\\ 1\end{bmatrix}
= \exp\!\left(\begin{bmatrix} A\,\mathrm{diag}(d) & u\\ 0 & 0\end{bmatrix}\right)
\begin{bmatrix} X^{n}\\ 1\end{bmatrix}, \qquad
R_h^{\mathrm{day}} = \sum_i u_i - \Big(\sum_i X^{n+1}_i - \sum_i X^{n}_i\Big) \qquad(6)
```

The day's respiration on the right is reported **solver-consistently** as the actual loss to the
atmosphere: inter-pool transfers cancel in the total, so it is what came in minus what stayed. Lignin
decays **in lockstep with its host structural carbon**, $`L_s\leftarrow L_s\sigma_s+\ell_s`$ with
$`\sigma_s=1-d_s`$ (Euler) or $`\exp(-d_s)`$ (EXPM), so $`f_\ell`$ is conserved under pure decay and
drifts only through inputs — the linear form under EXPM would drive $L$ negative at large $`d_s`$.

**Why the integral is the right object.** $\xi$ is a product of nonlinear functions of temperature and
moisture, both of which swing over a day, so by Jensen's inequality
$`\langle f_T(T) f_\theta(\theta)\rangle \neq f_T(\langle T\rangle) f_\theta(\langle\theta\rangle)`$: a
step that recomputed $\xi$ from daily-mean state would debit the pools by a *different* number than the
sub-daily flux actually delivered to the atmosphere, and ecosystem carbon would drift while every
individual budget still closed. Integrating $\xi$ as it happens removes the gap by construction rather
than bounding it; the tests run a diurnal cycle through both paths to show the daily-mean one does not.

The step also reports $`\lambda=\max_j d_j`$, the **largest fraction of any pool this step withdraws** —
the number that says whether freezing the pool across the slow step is sound. It is scale-free (the flux
is linear in the pool, so $\lambda$ does not grow as a pool empties, which is why bare-ground spin-up
works); a daily step on a mature stand withdraws a few tenths of a percent of the fastest pool, and past
1 forward Euler would drive that pool negative.

## 7. The fast/slow seam

Heterotrophic respiration must reach the canopy air **sub-daily**, but the pools may only be written
**daily** — the classic double-counting hazard. MEDS resolves it by having both ends read the same
object: once per day the fast driver copies the patch's pools into the fast state, and that **frozen**
copy is what every sub-step respires, through the very kernel the daily step's bookkeeping uses. Each
sub-step, `heterotrophic_respiration` evaluates eq. (2) on it, and the frozen build converts
$`R_h`$ to $`[\mathrm{\mu mol\,m^{-2}\,s^{-1}}]`$ with `kgcday_2_umols` for the canopy air's CO₂;
the driver integrates $`\xi_j`$ into the per-patch `decomp_xi%xi_days` (the $`\xi^{\mathrm{int}}`$ of
eq. (5)) and $`R_h`$, still in the pools' kgC m⁻² day⁻¹ so the seam check below compares like with
like, into `decomp_xi%rh_fast`, both reset once per slow step. Within a step the daily
`soil_carbon_step` is the **only** change to the real pools, debiting them with exactly that
$`\xi^{\mathrm{int}}`$ (it returns the new pools, and `meds_slow_update` writes them); at a calendar
boundary patch fusion and the treefall gaps blend them on area ([demography](demography.md) §4), and
the plants the boundary kills reach them through the litter (§5).

The day's total fast Rh therefore equals the daily pool debit **by construction**, and
$`\texttt{rh\_seam\_gap} = R_h^{\mathrm{day}} - \sum_k R_h(t_k)\Delta t_k \approx 0`$ is an
**assertion guard**, never a live correction; it closes to machine precision and the run reports its
worst value. Note what it is *not*: $\lambda$ (§6) catches the frozen-pool approximation degrading,
`rh_seam_gap` catches the contract being broken. This is the **single documented fast↔slow kernel seam
in the whole model** — one kernel called from both timescales, co-located with the pools it debits.

The gap is zero by construction only while the patch composition is fixed. On a day the patch
structure changes (the annual patch fusion, a termination or a disturbance), the fast loop has
accumulated `rh_fast` against one composition while the daily step debits another:
`blend_decomp_xi` and `blend_soil_carbon` area-weight the two ends separately, through a matrix that is
nonlinear in the lignin fraction. The gap is then nonzero and expected (8.4×10⁻⁴ kgC m⁻² at a year
rollover of a 50-year spin-up). A nonzero gap on an ordinary day is a broken contract; do not widen
the tolerance to absorb the structural-change days, which would blind the check on every other day.

## 8. Steady state (SASU) and traceability diagnostics

Under climatological-mean drivers the equilibrium is a direct linear solve rather than a
multi-millennial integration — and the same active-block machinery, evaluated at the *current* $\xi$
and $u$, gives the traceability trio:

```math
0 = \bar u + A\,\bar\xi K X_{ss} \Longrightarrow M X_{ss}[\mathcal{A}] = -\bar u[\mathcal{A}],\ \
M = \big(A\bar\xi K\big)[\mathcal{A},\mathcal{A}],\ \ \mathcal{A}=\{j: K_j>0\} \qquad(7)
```

```math
X_c = -\big(A\,\xi K\big)^{-1}u, \qquad X_p = X_c - X, \qquad
\tau_E = \frac{\sum_i X_{c,i}}{\sum_i u_i}\Big/\texttt{yr\_day} \qquad(8)
```

**The solve runs on the active sub-block only, because the full 7×7 system is singular** — an inert
pool's column of $`A\,\bar\xi K`$ is identically zero, so the inverse does not exist and a blind
back-substitution divides by zero on its diagonal. That is not a corner case; it is the *default*
scheme. Inert pools are held at their conserved value of zero. The reduced block is solved by a small
dense **Gaussian elimination with partial pivoting** — allocation- and LAPACK-free, and valid for
scheme 5's non-triangular back-transfers as well as the default block. A run can start from eq. (7)
(`[soil_carbon].spinup_steady`, a scalar climatological $`\bar\xi`$ and a constant litter estimate),
with two caveats: those means are **prescribed**, not fed back from a live canopy, and a mean-driver
equilibrium discards the covariance of periodic forcing.

In eq. (8), $`X_c`$ is the **storage capacity** — the equilibrium the patch is chasing under today's
climate and litter input, a moving attractor rather than a fixed point; $`X_p`$ is the **storage
potential**, the remaining sink (negative means a source); $`\tau_E`$ is the ecosystem **residence
time** in years, stock over throughput. The transient relaxes toward $`X_c`$ at a rate set by
$`A\,\xi K`$, making the trio a decomposition of *why* two runs store different carbon. The model
does not report the trio; it follows from the pools, the litter input and the scalars the output
holds.

## 9. Audits, and the difference between conservation and plausibility

Every step fills a `soil_carbon_out_t`: the new pools, `rh_seam_gap` (§7),
$\lambda$ with the pool that attained it, and the lignin passive-tracer residual
$`\max_s|\Delta L_s-(\ell_s-(1-\sigma_s)L_s)|`$, zero by construction for either solver. The carbon
check that can actually fail is the *independent* one, in the tests: recompute $`R_h`$ from
$`e_r\,\xi^{\mathrm{int}}K X^n`$ and compare against $`\sum u-\Delta X`$ — comparing eq. (6) against
itself cannot catch dropped litter or a misattributed transfer.

Conservation is not plausibility: a runaway that conserves carbon closes every budget while moving
absurd amounts of it. `soil_carbon_bad_pool` answers what a budget cannot — 0 when every pool is
physically possible, else the **index** of the first that is not (NaN, negative beyond a round-off
tolerance of $`10^{-9}\,\mathrm{kgC\,m^{-2}}`$, or above a divergence ceiling of
$`10^{4}\,\mathrm{kgC\,m^{-2}}`$, ~500× the richest real soil; a temperate forest's holds 20–40). It
returns an index, not a flag, so the caller's message can name the pool, and which pool went bad is
the first clue to which operator did it. Round-off negatives must pass — a daily
Euler step produces them routinely — and the ceiling catches divergence, not implausible science.

## 10. What is not here

- **Nitrogen.** The soil is carbon-only: no nitrogen pools, no C:N stoichiometry, no nitrogen
  limitation of decomposition or of plant growth. The nitrogen cycle and its coupling to carbon are
  planned after the move to ED-MORE.
- **DAMM, or any other respiration form.** The matrix form (§7) is the only heterotrophic
  respiration. Davidson's dual Arrhenius/Michaelis–Menten scheme (soluble-C and O₂ limitation, a
  unimodal moisture response *emergent* rather than imposed) is not available, nor are separate Q10
  or capped-exponential respiration kernels: those two temperature responses live inside $\xi$,
  eq. (3).
- **No vertical resolution.** Pools are per-patch scalars, as in ED2; only the environmental scalar
  sees layer-resolved temperature and moisture.
- **No fire**, so nothing yet consumes the above-ground pools; treefall is the only disturbance.

See `docs/ROADMAP.md` for where these stand.

## Parameters (`[soil_carbon]`)

| Key | Meaning |
|---|---|
| `soil_carbon_on`; `decomp_scheme` / `step_solver` | master switch (default on); `"ed2"` or `"century5"` / `"euler"` or `"expm"` |
| `k_fast`, `k_struct`, `k_micr`, `k_slow`, `k_pass` | baseline decay rates $`[\mathrm{yr^{-1}}]`$ |
| `er_fast`, `er_struct_nonlig`, `er_struct_lig`, `er_slow`, `er_pass`, `e_lignin` | respired fractions + the lignin brake |
| `er_micr_*`, `xsand`; `fx_micr_pass_*`, `fx_slow_pass_*`, `xclay` | scheme-5 texture control (sand → microbial respired fraction, clay → routing to passive) |
| `rh_q10`, `rh_t_ref`, `resp_temp_increase`, `resp_temp_ref`, `resp_opt_water`, `resp_water_{below,above}_opt` | the $\xi$ responses, eqs. (3)–(4) |
| `spinup_steady`, `spinup_xi`, `spinup_{labile,struct}_*` | SASU cold start, §8 |

The defaults are the `"ed2"` scheme's. `decomp_scheme = "century5"` changes the topology and the
temperature form but not the rates, so set them with it. Its reference values (ED2's scheme 5) are
`k_fast` 12.0, `k_struct` 1.5, `k_micr` 6.0, `k_slow` 0.15 and `k_pass` 0.012 yr⁻¹; `er_fast` 0.55,
`er_struct_nonlig` 0.50, `er_struct_lig` 0.20, `er_slow` 0.55 and `er_pass` 0.55; `er_micr_int` 0.60
and `er_micr_slp` 0.17.

**`soil_carbon_on = false` is no soil carbon, not a coarser one.** Litter is discarded at the slow
step and heterotrophic respiration is zero, so ecosystem respiration carries only its autotrophic
part. A prescribed pool respired through an empirical Q10 form would not be a coarser approximation
either: the pools are empty and take no litter, so its CO₂ would be carbon from nothing. Over a year
at Ithaca, annual-mean NEE then reads −4.657 µmol m⁻² s⁻¹ against −2.464 with soil carbon on, an
89 % stronger apparent sink (Rh 0 against 0.833 kgC m⁻² yr⁻¹), for no saving in wall clock.

The litter-partition traits are **per PFT**, not run constants: `f_labile_leaf`, `f_labile_stem`,
`aboveground_frac`, `struct_lignin_frac`.

## Where the code is

| Concept | Routine |
|---|---|
| $\xi$; and $A$, $K$, $`e_r`$ | `meds_respiration`: `decomposition_env_scalar`, `decomposition_transfer_matrix` |
| litter → source term $u$; necromass → destinations | `meds_soil_carbon`: `build_litter_input`, `necromass_to_litter`, `add_necromass` |
| daily pool advance (Euler / EXPM) + audits | `meds_soil_carbon`: `soil_carbon_step` |
| instantaneous $`R_h`$ (the fast↔slow seam kernel) | `meds_respiration`: `heterotrophic_respiration` |
| steady state (SASU) | `meds_soil_carbon`: `solve_soil_carbon_steady_state` |
| pool plausibility | `meds_soil_carbon`: `soil_carbon_bad_pool`, `soil_carbon_pool_name` |
| pools, indices, audit/diagnostic records | `meds_soil_carbon`: `soil_carbon_out_t`; the pool layout `N_SOIL_POOL`, `IP_*` in `meds_biogeochem_params`; `soil_carbon_state_t` / `decomp_xi_t` / `litter_input_t`, with the pool vector's `pack_pool_vector` / `unpack_pool_vector`, in `meds_patch_state_vars` |
| `[soil_carbon]` bundle + selectors | `meds_biogeochem_params`: `soil_carbon_params_t`, `DECOMP_SCHEME_*`, `DECOMP_STEP_*` |
| daily driver, seam guard, ledger | `meds_slow_state`: `advance_soil_carbon`; `meds_slow_update`: `commit_soil_carbon`; the run's worst seam checks in `meds_slow_ledger`: `slow_ledger_soil_check` |
| litter accumulation | `meds_slow_state`: `compute_carbon_allocation`, `accumulate_mortality_litter`; kills in `meds_slow_update`: `terminate_cohorts`, `apply_patch_disturbance` |
| sub-daily $`R_h`$ + $`\xi`$ accumulation | `meds_fast_integ_frozen`: `patch_heterotrophic_respiration`; `meds_fast_state`: `accumulate_fast_to_slow` |
| output; tests | `soilc_*_site`, `mean_soilc_total_pa`, `mean_rh_py` (`GRP_CARBON`); `test/test_soil_carbon.f90`, `test/test_slow_step.f90`, `test/test_column_co2.f90` |

## References
- ED2 `../ED2/ED/src/`: `dynamics/soil_respiration.f90` (environmental scalars, sub-daily
  accumulation), `dynamics/vegetation_dynamics.f90` (`update_C_and_N_pools`), `memory/decomp_coms.f90`
  + `init/ed_params.f90` (pools, decay rates, respired/transfer fractions).
- Parton et al. (1987), *SSSAJ* 51:1173 — the CENTURY pool structure. Bolker, Pacala & Parton (1998),
  *Ecol. Appl.* 8:425; Koven et al. (2013), *Biogeosciences* 10:7109 — the five-active topology.
- Luo et al. (2022), *JAMES* — the matrix formulation; Xia et al. (2012), *GMD* 5:1259 — semi-analytic
  spin-up; Luo et al. (2017), *Biogeosciences* 14:145 — capacity and potential; Davidson et al. (2012),
  *Glob. Change Biol.* 18:371 — DAMM.
