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
donor/transfer matrix (§3). The input is **litter carbon, never GPP** — autotrophic respiration was
removed upstream in the plant carbon budget.

## 2. The seven pools and the lignin sub-tracer

The state is `soil_carbon_t`, ordered litter → SOM → passive. The seven-pool ceiling matches ED2's
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
pools, so it stays outside the mass vector and keeps its own balance (§6). The nitrogen twin is shaped
into the type but inert (§10); `pack_pool_vector` / `unpack_pool_vector` marshal the named fields to
and from the length-7 array the kernels use.

## 3. The transfer matrix, and the respired complement

`assemble_transfer_matrix` builds $A$, $K$ and the respired fractions together. Each donor loses *all*
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

`assemble_env_scalar` returns one scalar per pool: temperature × moisture/oxygen × (structural only) a
lignin brake, over the relative saturation $`m`$ clipped to $`[0,1]`$.

```math
f_T = \begin{cases} q_{10}^{\,(T-T_{\mathrm{ref}})/10} & \text{scheme 5, UNCAPPED Q10}\\
\min\!\big(1,\ \exp[a_T(T-T_{\mathrm{sat}})]\big) & \text{schemes 0--4, ED2 capped exp}\end{cases}
\qquad
f_\theta = \begin{cases} \exp[(m-m_{\mathrm{opt}})\,b_{\mathrm{dry}}] & m \le m_{\mathrm{opt}}\\
\exp[(m_{\mathrm{opt}}-m)\,b_{\mathrm{wet}}] & m > m_{\mathrm{opt}}\end{cases}
```

```math
m = \frac{\theta-\theta_{\mathrm{dry}}}{\theta_{\mathrm{sat}}-\theta_{\mathrm{dry}}}, \qquad
\xi_j = f_T(T_j)\,f_\theta(\theta)\,\Lambda_j, \qquad
\Lambda_j = \begin{cases}\exp(-e_{\mathrm{lig}}f_{\ell,j})\cdot f_{\mathrm{decomp}} & j\in\{3,4\}\\
1 & \text{otherwise}\end{cases} \qquad(3)
```

One one-sided exponential covers both moisture limbs — dry-side substrate diffusion, wet-side anoxia —
and it is the **same functional form the fast Rh kernel uses**, so the two timescales reconcile at
matched inputs. $`T_j`$ is the **ground surface** temperature for the above-ground pools (1, 3) and the
**root-weighted soil** temperature for the rest (ED2's `A_decomp` / `B_decomp` split), $\theta$ the
depth-weighted column mean, $`f_{\mathrm{decomp}}`$ the N immobilization brake (identically 1 while N
is off). The invariant is $`\xi>0`$, not $`\xi\le1`$: scheme 5's Q10 is uncapped.

**Placement is physics.** The scalar sits *inside* the transfer, $`A\,(\xi K X)`$, never
$`\xi\,(A K X)`$: **each donor decomposes at its own rate $`\xi_j K_j`$, and only then is the loss
split** into respired and transferred parts. $\xi$ and $K$ commute with each other but neither commutes
with $A$, so the receiver-scalar ordering silently breaks mass balance and every capacity/residence
formula. The tests assert it: the correct form's total loss equals $`R_h`$ of eq (2); the other's does not.

## 5. Litter input

Dead plant carbon becomes litter in one place — `necromass_to_litter` (`meds_litter_partition`), whose
arguments are plain carbon amounts and plain PFT traits, which is what lets both the slow driver and the
demographic operators call it. Leaf necromass is pooled with **storage** necromass (both canopy-
associated and labile-eligible; dead seed and seedling arrive here) and split by `f_labile_leaf`;
fine-root necromass reuses `f_labile_leaf` but lands **entirely below-ground**; wood/CWD splits by
`f_labile_stem`. Each stream then splits above/below by `aboveground_frac`, and only the structural ones
carry lignin — `struct_lignin_frac` of each.

`build_litter_input` maps those six numbers onto the matrix source term — four carbon streams into $u$,
two lignin fluxes alongside — and does **no** chemistry of its own, since the splits are per-cohort
decisions made before summing to the patch. So **only pools 1–4 ever receive litter**; microbial, slow
and passive are fed by transfers alone. Carbon arrives by two routes: turnover/shed and continuous
background mortality accumulate into `patch%litter_in` and enter through $u$, while **cull-termination
and disturbance kills** add necromass *directly* onto `patch%soil_carbon` (the demographic operators
cannot link biogeochemistry). Hence the audit field `litter_in_matrix` — it is $`\sum_i u_i`$ only.

## 6. The daily step, and why it consumes an accumulated integral

`soil_carbon_step` advances one patch by one day, and the decay it applies is **not** an instantaneous
$\xi$ — it is the fast loop's **accumulated per-pool loss integral** $`\xi^{\mathrm{int}}`$, with its
dimensionless decay fraction $`d`$. The default solver (`DECOMP_STEP_EULER`) is ED2-faithful forward
Euler on the frozen operator ($`a_{ii}=-1`$ removes donor $i$'s own loss, the off-diagonals transfer):

```math
\xi^{\mathrm{int}}_j = \int_{\mathrm{day}} \xi_j(t)\,dt \ [\mathrm{day}], \qquad
d_j = \xi^{\mathrm{int}}_j K_j \ [-], \qquad
X^{n+1}_i = X^n_i + u_i + \sum_j a_{ij}\,d_j\,X^n_j \qquad(4)
```

For **large accelerated steps** `DECOMP_STEP_EXPM` takes the *exact* exponential of that same frozen
operator, through an augmented matrix that carries the constant input without inverting $A$ — which
matters because $`A\,\mathrm{diag}(d)`$ is singular whenever a pool is inert:

```math
\begin{bmatrix} X^{n+1}\\ 1\end{bmatrix}
= \exp\!\left(\begin{bmatrix} A\,\mathrm{diag}(d) & u\\ 0 & 0\end{bmatrix}\right)
\begin{bmatrix} X^{n}\\ 1\end{bmatrix}, \qquad
R_h^{\mathrm{day}} = \sum_i u_i - \Big(\sum_i X^{n+1}_i - \sum_i X^{n}_i\Big) \qquad(5)
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
sub-step, `heterotrophic_respiration_matrix` evaluates eq (2) on it and converts $`R_h`$ to
$`[\mathrm{\mu mol\,m^{-2}\,s^{-1}}]`$ with `kgCday_2_umols` for the canopy-air CO₂ twin; the driver
integrates $`\xi_j`$ into the per-patch `xi_accum` (the $`\xi^{\mathrm{int}}`$ of eq (4)) and $`R_h`$
into `rh_fast_accum`, both reset once per slow step; and the daily `soil_carbon_step` is the **sole
writer** of the real pools, debiting them with exactly that $`\xi^{\mathrm{int}}`$.

The day's total fast Rh therefore equals the daily pool debit **by construction**, and
$`\texttt{rh\_seam\_gap} = R_h^{\mathrm{day}} - \sum_k R_h(t_k)\Delta t_k \approx 0`$ is an
**assertion guard**, never a live correction; it closes to machine precision and the run reports its
worst value. Note what it is *not*: $\lambda$ (§6) catches the frozen-pool approximation degrading,
`rh_seam_gap` catches the contract being broken. This is the **single documented fast↔slow kernel seam
in the whole model** — one kernel called from both timescales, co-located with the pools it debits.

## 8. Steady state (SASU) and traceability diagnostics

Under climatological-mean drivers the equilibrium is a direct linear solve rather than a
multi-millennial integration — and the same active-block machinery, evaluated at the *current* $\xi$
and $u$, gives the traceability trio:

```math
0 = \bar u + A\,\bar\xi K X_{ss} \Longrightarrow M X_{ss}[\mathcal{A}] = -\bar u[\mathcal{A}],\ \
M = \big(A\bar\xi K\big)[\mathcal{A},\mathcal{A}],\ \ \mathcal{A}=\{j: K_j>0\} \qquad(6)
```

```math
X_c = -\big(A\,\xi K\big)^{-1}u, \qquad X_p = X_c - X, \qquad
\tau_E = \frac{\sum_i X_{c,i}}{\sum_i u_i}\Big/\texttt{yr\_day} \qquad(7)
```

**The solve runs on the active sub-block only, because the full 7×7 system is singular** — an inert
pool's column of $`A\,\bar\xi K`$ is identically zero, so the inverse does not exist and a blind
back-substitution divides by zero on its diagonal. That is not a corner case; it is the *default*
scheme. Inert pools are held at their conserved value of zero. The reduced block is solved by a small
dense **Gaussian elimination with partial pivoting** — allocation- and LAPACK-free, and valid for
scheme 5's non-triangular back-transfers as well as the default block. A run can start from eq (6)
(`[soil_carbon].spinup_steady`, a scalar climatological $`\bar\xi`$ and a constant litter estimate),
with two caveats: those means are **prescribed**, not fed back from a live canopy, and a mean-driver
equilibrium discards the covariance of periodic forcing.

In eq (7), $`X_c`$ is the **storage capacity** — the equilibrium the patch is chasing under today's
climate and litter input, a moving attractor rather than a fixed point; $`X_p`$ is the **storage
potential**, the remaining sink (negative means a source); $`\tau_E`$ is the ecosystem **residence
time** in years, stock over throughput. The transient relaxes toward $`X_c`$ at a rate set by
$`A\,\xi K`$, making the trio a decomposition of *why* two runs store different carbon.

## 9. Audits, and the difference between conservation and plausibility

Every step fills a `soilc_audit_t`: `litter_in_matrix`, `dC_pool`, `rh_out`, `rh_seam_gap` (§7),
$\lambda$ with the pool that attained it, and the lignin passive-tracer residual
$`\max_s|\Delta L_s-(\ell_s-(1-\sigma_s)L_s)|`$, zero by construction for either solver. The carbon
check that can actually fail is the *independent* one, in the tests: recompute $`R_h`$ from
$`e_r\,\xi^{\mathrm{int}}K X^n`$ and compare against $`\sum u-\Delta X`$ — comparing eq (5) against
itself cannot catch dropped litter or a misattributed transfer.

Conservation is not plausibility: a runaway that conserves carbon closes every budget while moving
absurd amounts of it. `soil_carbon_bad_pool` answers what a budget cannot — 0 when every pool is
physically possible, else the **index** of the first that is not (NaN, negative beyond a round-off
tolerance of $`10^{-9}\,\mathrm{kgC\,m^{-2}}`$, or above a divergence ceiling of
$`10^{4}\,\mathrm{kgC\,m^{-2}}`$, ~500× the richest real soil). Round-off negatives must pass — a daily
Euler step produces them routinely — and the ceiling catches divergence, not implausible science.

## 10. What is not here

- **Nitrogen.** `soil_carbon_t` and `litter_input_t` carry N fields and `decomp_opts_t` the C:N ratios
  and `n_cycle_on`, but **no kernel reads the flag**: $`f_{\mathrm{decomp}}\equiv1`$ and every N field
  stays zero. The cycle is shaped in, not implemented.
- **DAMM.** `heterotrophic_respiration_damm` implements Davidson's dual Arrhenius/Michaelis-Menten
  scheme (soluble-C and O₂ limitation, unimodal moisture response *emergent* rather than imposed), and
  `heterotrophic_respiration_flux` dispatches to it over `HR_Q10` / `HR_EXP_ED2`. Neither has a config
  key or a production caller — the matrix form (§7) is what the fast loop respires.
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
| `rh_q10`, `rh_t_ref`, `resp_temp_increase`, `resp_temp_ref`, `resp_opt_water`, `resp_water_{below,above}_opt` | the $\xi$ responses, eq (3) |
| `n_cycle_on`, `c2n_*`, `n_immobil_supply_scale` | nitrogen (shaped in, inert — §10) |
| `spinup_steady`, `spinup_xi`, `spinup_{labile,struct}_*` | SASU cold start, §8 |

The litter-partition traits are **per PFT**, not run constants: `f_labile_leaf`, `f_labile_stem`,
`aboveground_frac`, `struct_lignin_frac`.

## Where the code is

| Concept | Routine |
|---|---|
| $\xi$; and $A$, $K$, $`e_r`$ | `meds_soil_biogeochem`: `assemble_env_scalar`, `assemble_transfer_matrix` |
| litter → source term $u$; necromass → destinations | `meds_soil_biogeochem`: `build_litter_input`; `meds_litter_partition`: `necromass_to_litter` |
| daily pool advance (Euler / EXPM) + audits | `meds_soil_biogeochem`: `soil_carbon_step` |
| instantaneous $`R_h`$ (the fast↔slow seam kernel) | `meds_soil_biogeochem`: `heterotrophic_respiration_matrix` |
| steady state; capacity / potential / residence time | `meds_soil_biogeochem`: `solve_soil_carbon_steady_state`, `soil_carbon_diagnostics` |
| pool plausibility | `meds_soil_biogeochem`: `soil_carbon_bad_pool`, `soil_carbon_pool_name` |
| pools, indices, audit/diagnostic records | `meds_biogeochem_types`; `soil_carbon_t` / `xi_accum_t` / `litter_input_t` in `meds_column_state_types` |
| `[soil_carbon]` bundle + selectors | `meds_biogeochem_opts`: `decomp_opts_t`, `DECOMP_SCHEME_*`, `DECOMP_STEP_*` |
| daily driver, seam guard, ledger | `meds_biogeochem_dynamics`: `advance_biogeochem_dynamics` |
| litter accumulation | `meds_vegetation_dynamics`; kills in `meds_demography_{cohort,patch}_fusefiss` |
| sub-daily $`R_h`$ + $`\xi`$ accumulation | `meds_fast_prepass`: `patch_heterotrophic_respiration`; `meds_fast_dynamics` |
| output; tests | `soilc_*_site`, `soilc_total_patch`, `rh_site` (`GRP_CARBON`); `test/test_soil_biogeochem.f90`, `test/test_biogeochem_dynamics.f90`, `test/test_column_co2.f90` |

## References
- ED2 `../ED2/ED/src/`: `dynamics/soil_respiration.f90` (environmental scalars, sub-daily
  accumulation), `dynamics/vegetation_dynamics.f90` (`update_C_and_N_pools`), `memory/decomp_coms.f90`
  + `init/ed_params.f90` (pools, decay rates, respired/transfer fractions).
- Parton et al. (1987), *SSSAJ* 51:1173 — the CENTURY pool structure. Bolker, Pacala & Parton (1998),
  *Ecol. Appl.* 8:425; Koven et al. (2013), *Biogeosciences* 10:7109 — the five-active topology.
- Luo et al. (2022), *JAMES* — the matrix formulation; Xia et al. (2012), *GMD* 5:1259 — semi-analytic
  spin-up; Luo et al. (2017), *Biogeosciences* 14:145 — capacity and potential; Davidson et al. (2012),
  *Glob. Change Biol.* 18:371 — DAMM.
