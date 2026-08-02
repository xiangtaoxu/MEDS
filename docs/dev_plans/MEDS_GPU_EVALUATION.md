# MEDS on the GPU — an evaluation, MEASURED 2026-08-02

**Status: ❌ NOT VIABLE for the current model. Do not build BB2/BB3 as scoped.**

This document prices GPU offload for MEDS against measurement rather than intuition. It was written
because `MEDS_GPU=gpu` has been a supported build option since the port of `src/biophysics`, has never
been performance-evaluated, and the surrounding documentation implies a capability that does not exist.

Everything below was measured on the dev box (WSL2, i7-11370H 4c/8t + RTX 3050 Ti Laptop, ifx 2026 /
nvfortran 25.11 / CUDA 13.0) against `main` at `d06dac6`. Both nvfortran back ends were configured and
built fresh and pass **37/37** ctest.

Section 8 answers a separate question that arrived after the first pass — *does a stand with thousands
of cohorts change the answer?* — and section 9 answers *does FP32?* Those two, together, are the only
place where the verdict moves.

---

## 1. Verdict

| question | answer |
|---|---|
| Does the GPU build run faster today? | **No — 1.4× slower** than plain ifx (49.8 s vs 36.2 s on a 1-year run). |
| Is that a tuning problem? | **No.** Four independent walls, three of them quantitative. |
| Would BB2 (offload `column_fast_step`) fix it? | **No.** Launch latency alone exceeds the CPU cost of the step it replaces. |
| Does a 1000+ cohort stand change it? | **Partly** — it removes wall 1 and weakens wall 3. Walls 2 and 4 stand. §8. |
| Can the GPU be used as "20 slow cores", one patch per SM? | **No — 26× slower than 4 CPU cores at 20 patches**; one GPU lane is ~140× slower than one CPU core. §8d. |
| Does FP32 change it? | **Yes, wall 2 — decisively (56× on device).** But it is a mixed-precision program, not a flag. §9. |
| Is there *any* workload where a GPU wins? | **Yes: gridded/regional runs.** MEDS has no site axis. §10. |
| What should be done instead? | Thread/vectorise the **cohort** axis on the CPU, and remove the allocator traffic. §12. |

---

## 2. What is actually on the GPU today

Exactly **one** `!$omp target` region exists in the 24,569 lines of `src/`:

`src/core/meds_core_state_update.f90:117` — `update_cohort_states_kernel`, the slow-loop
`state += rate * dt` applier. It is offload-eligible for the reason CLAUDE.md gives: it takes bare
arrays, so the `map` clauses are clean.

Verified on device with `NVCOMPILER_ACC_NOTIFY=3` over 4 simulated days:

```
launch CUDA kernel  file=.../meds_core_state_update.f90  function=update_cohort_states_kernel  line=117
                    grid=<<<1,1,1>>>  block=<<<128,1,1>>>
4 launches | 80 uploads | 40 downloads | 920 B each | 110 kB total
```

`cuobjdump --list-elf` confirms it is the **only** kernel in the binary.

Three things follow:

- **Occupancy is 0.4%.** `cudaGetDeviceProperties` reports 20 SMs × 1536 threads = **30,720 resident
  threads to fill the device**. This kernel supplies one block of 128, of which 115 do work.
- **It runs once per simulated day** — it is a slow-loop kernel, not a fast-loop one.
- **The slow loop is 0.08% of runtime.** Same 1-year configuration, ifx Release:

  | | wall |
  |---|---|
  | full run | 36.24 s |
  | `fast_biophysics_on = false` (slow loop only) | **0.03 s** |

  Its 30 map transfers per call cost roughly as much as the entire slow loop does. The offload is a
  net negative, merely too small to measure.

**Two documentation statements are wrong and should be corrected.** `CMakeLists.txt:6` says *"The hot
kernels carry explicit OpenMP `target` regions"* (plural) and `CLAUDE.md:76` calls `update_cohort_states`
*"the hot kernel"*. There is one kernel and it is not hot.

---

## 3. Measured — three back ends, identical workload

Restart from the 50-year spin-up state (12 patches, 115 cohorts, 20 soil layers), 1 model year,
`dt_fast = 900 s`, `time_integrator = "ark"`, diagnostic output off.

| build | `n_threads = 1` | `n_threads = 4` |
|---|---|---|
| ifx Release (serial) | **36.2 s** | — |
| ifx Release `-DMEDS_OPENMP=ON` | 38.6 s | **19.1 s** |
| nvfortran `-mp` (multicore) | 49.3 s | 18.5 s |
| nvfortran `-mp=gpu -gpu=mem:separate` | **49.8 s** | 19.5 s |

GPU ≡ multicore to within noise. The 37% serial penalty is nvfortran-vs-ifx codegen, not the device.

The test suite shows the device tax directly: **4.63 s on the GPU build vs 1.95 s on multicore** — 2.4×,
entirely device initialisation and launch latency on tiny kernels.

---

## 4. Wall 1 — there is not enough parallel work (3 orders of magnitude short)

From the production state file `spinup-S-20740701000000.nc`: **12 patches, 115 cohorts, 20 soil layers,
`N_HYDRO = 3`.**

| axis | width | independent within a `dt_fast`? |
|---|---|---|
| **site** | **1** | no such axis exists — `site_t` is a single instance; there is no grid or ensemble dimension anywhere in `src/` |
| **patch** | 12 | yes — this is exactly what the §7 host threading already exploits |
| **cohort** | ~10/patch | partly — leaf gas exchange, hydraulics and respiration are per-cohort independent; the aerodynamic cascade is not |
| **soil layer** | 20 | no — Thomas recurrence |

The whole site's fast prognostic state is about **1,000 doubles (~8 kB)** — it fits in one SM's L1.
There is no problem here to spread across 2,560 CUDA cores.

The §7 patch threading already harvested this axis and hit the *hardware* ceiling at 4 cores (3.03×
measured aggregate throughput, 67% efficient). A GPU wants 30,720 concurrent work items; MEDS offers 12
coarse ones.

---

## 5. Wall 2 — this GPU is slower than this CPU in FP64

`wp = real64` throughout (`src/shared/base/meds_kinds.f90:13`). Consumer Ampere (GA107) runs FP64 at
1/64 of its FP32 rate. Measured with an 8-independent-chain FMA kernel (`scratchpad/fp64peak8.f90`):

| | measured GFLOP/s FP64 | spec peak |
|---|---|---|
| RTX 3050 Ti Laptop, `-mp=gpu` | **85.0** | ~119 |
| i7-11370H, ifx `-O3 -xHost`, 4 threads | 144.7 | ~256 |
| i7-11370H, ifx `-O3 -xHost`, 8 threads | **177.6** | — |

**The GPU delivers 0.48× the host's FP64 throughput.** A hypothetical perfect, zero-overhead offload of
100% of MEDS would run about half as fast on this machine.

This wall is hardware-specific — it disappears on a datacenter part (A100 FP64 ≈ 9.7 TFLOP/s) — and it
is the one wall FP32 removes outright (§9).

---

## 6. Wall 3 — one kernel launch costs more than an entire time step

| | cost |
|---|---|
| empty `!$omp target` region, data already resident | **59.3 µs** |
| `target` + `map(tofrom:)`, n = 12 … 1150 | 94.9 … 102.8 µs |
| raw CUDA empty kernel, synchronous round-trip | 64.1 µs |
| raw CUDA empty kernel, async / pipelined | 18.5 µs |
| **one whole `column_fast_step` on the CPU, 10 cohorts** | **58 µs** (`MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §5c(i)) |

A single empty GPU launch costs more than the entire column solve it would replace.

The raw-CUDA numbers confirm this is **WSL2 device round-trip latency**, not OpenMP-runtime bloat —
native Linux would be ~5 µs. But even at 5 µs the arithmetic does not close at production cohort counts.
A BB2-style offload of `column_fast_step` needs roughly 40–50 device round-trips per step (pre-pass +
hydraulics batch + scratch hydrology, then ~2 accepted sub-steps × 3 ESDIRK stages ×
{`surface_derivs`, `column_derivs`, Newton ≤ 8 iterations, Thomas, veg energy}). Against a CPU cost of
12 × 58 µs = 0.70 ms for the whole site per `dt_fast`:

- at WSL's 59 µs → 2.4–3.0 ms, **~4× slower**;
- at native Linux's 5 µs → 0.2–0.25 ms of pure launch overhead, *plus* kernels carrying 12–115 lanes on
  a device that needs 30,720. Every one is latency-bound, not throughput-bound.

This is the wall that a large cohort count weakens — see §8.

---

## 7. Wall 4 — the code is structurally hostile to `target`

These are what BB2 would have to fix. They are real, and the plans already name most of them.

**Allocatable-component derived types everywhere.** `map` cannot cleanly handle them, which is precisely
why `update_cohort_states` is the only offloaded kernel.

| type | allocatable components |
|---|---|
| `column_frozen_t` | 16 |
| `surface_frozen_t` | 14 |
| `column_cohort_t` | 8 |
| `patch_biophys_t` | 6 |
| `column_tend_t` | 5 |
| `column_state_t` | 4 |

**Heap allocation inside the hot path.** `build_column_frozen` issues ~26 `allocate` calls *per patch
per `dt_fast`*. There is no `allocate` on a device. This is also the measured **~24% allocator self-time**
in §5c(iv) — the single largest line item in the CPU profile is heap traffic, not physics.

**`error stop` in device-reachable code** — 4 in `meds_fast_ark`, 4 in `meds_plant_hydraulics`, 2 in
`meds_soil_water`, 1 in `meds_soil_energy`. Each blocks `pure`, and `pure` is the precondition for
offload. `MEDS_NUMERICS_SCOPING.md` §8a already flagged this for `solve_plant_water`.

**Serial recurrences at both inner widths.** `thomas_solve` (`meds_numerics.f90:74`) is a strict
forward/back substitution over the soil column. `canopy_aerodynamics` is a top-down wind cascade over
cohorts with an embedded Monin–Obukhov fixed-point iteration; §11.3's own inventory marks it
**"NOT batchable — cohort-coupled wind cascade."**

**Adaptivity is warp divergence by construction.** `adaptive_ark_march` runs a per-patch sub-step count
driven by a WRMS controller with rejections, and the plant-hydraulics sub-step count is measured at
**1.0–1.2 normally and 136.6 on a collapsed store** (§5c(v)). Under SIMT the whole warp pays the maximum.
BB3's warp-synchronous fixed-count profile is the acknowledged answer — and adopting it means giving up
the adaptive controller that PRs #65/#66/#90 were built around.

---

## 8. Does a stand with thousands of cohorts change the answer?

**Partly — it removes wall 1 and weakens wall 3. It does not touch walls 2 and 4, and it makes the
divergence problem qualitatively worse.**

Thousands of cohorts is reachable **by configuration alone**: `max_cohort = 0` disables fusion outright
(`meds_core_cohort_fusefiss.f90:102`). No code change is needed to get there — but note that it moves
MEDS from a size-and-age-structured model toward an individual-based one, which is a scientific scoping
decision, not just a performance knob.

### 8a. ⚠️ The synthetic cohort sweep DEGENERATES above n ≈ 30 — and that is itself the finding

The intent was to extend §5c(i)'s cost split to n = 10 … 10,000 (`scratchpad/meas_ncoh.f90`, adapted from
the T1-REDUX harness: `dt_fast = 900 s`, ARK, 10 soil layers, total LAI 5.1 and total plant density held
fixed and split over `n` cohorts). It reproduces §5c(i) exactly at the validated sizes and then falls off
a cliff:

| n | `build_column_frozen` | `adaptive_ark_march` | WHOLE | march `nsteps`/`nrej` | `hydro_nsub` per cohort | ψ_leaf(1) |
|---|---|---|---|---|---|---|
| 10 | 21.7 µs | 30.3 µs | **57.1 µs** | 2 / 0 | 1.0 | −0.20 MPa |
| 30 | 39.4 | 66.7 | **112.0** | 2 / 0 | 1.0 | −0.29 MPa |
| 60 | 3398 | 68298 | **70539** | 20 / 3 | **200 (the cap)** | **−15000 MPa** |
| 100 | 4997 | > 10⁶ | — | 28 / 1 | **200 (the cap)** | **−15000 MPa** |

The n = 10 and n = 30 rows agree with §5c(i)'s independent measurement to ~1% (57.1 vs 57.95; 112.0 vs
110.85), so the harness is sound. **The fixture is not.** Between n = 30 and n = 60 the synthetic stand
desiccates: ψ_leaf hits the `psi_from_rwc` clamp artefact (§6/E1 — RWC floored at `1e-4`, ψ returned as a
floor marker three to four orders outside anything physical), `hydro_nsub` pins at `max_substep = 200`
for every cohort, and the ARK march goes from 2 sub-steps to 28.

**Do not read a cost law off those rows.** This is the fixture flaw class that
`feedback_fixture_must_mirror_driver` and `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §5c's own "⚠️ A FIXTURE
FLAW invalidated part of the earlier record" both warn about — a synthetic stand that holds stand-level
totals fixed while subdividing does not stay a physically coherent stand.

**But it is direct evidence for the divergence argument in §8b.** It shows how *easily* a many-cohort
stand falls into §5c(v)'s collapse regime, and what happens when it does: **200× on hydraulics and 14× on
the march**, arriving as a binary switch rather than a trend, and applying to whole groups of cohorts at
once. Under SIMT, one collapsed cohort in a warp makes all 32 lanes pay the 200.

### 8a-bis. What can be said about scaling, from the validated fixture

Fitting §5c(i)'s trustworthy rows (n = 1, 3, 10, 30 at `dt_fast = 900 s`: 36.96, 52.97, 57.95,
110.85 µs) gives a marginal cost of **≈ 2.65 µs per cohort on a fixed per-patch base of ≈ 32 µs**.
Extrapolating that line — and it *is* an extrapolation, from a range 30× smaller than the target:

| site cohorts/patch | per-patch serial base | per-cohort work | per-cohort share |
|---|---|---|---|
| 10 (today) | 32 µs | 27 µs | 45% |
| 100 | 32 | 265 | 89% |
| 1,000 | 32 | 2,650 | **99%** |

So **Amdahl is favourable at high cohort count**: the per-patch serial remainder (CAS box, soil Thomas,
surface Newton, and the aerodynamic cascade's sequential part) becomes negligible, and essentially all the
work lands on an axis that is independent per cohort. That is the one genuinely encouraging fact for a
GPU — and it applies just as much to CPU vectorisation (§8c).

**Caveat that a real measurement would have to settle:** this linear fit assumes the per-cohort *states*
stay well-conditioned as `n` grows. §8a shows they need not, and the collapse cliff is worth 200× on the
affected cohorts. A trustworthy cost law needs a fusion-off run of the actual driver
(`max_cohort = 0`, long enough to accumulate thousands of cohorts naturally), not a subdivided fixture.

### 8b. What this means for a GPU

**Wall 1 dissolves — but only at the top of the range.** The relevant quantity is the site's *total*
cohort count, since the per-cohort kernels are independent across patches as well as within them. Against
the 30,720 resident threads the device needs:

| site total cohorts | device fill (per-cohort kernels) |
|---|---|
| 115 (today) | 0.4% |
| 1,000 | 3% |
| 10,000 | 33% |
| 30,000+ | saturated |

So "thousands of cohorts" is the threshold at which the device stops being empty, and tens of thousands
is where it is actually busy. For the per-cohort kernels — leaf gas exchange, hydraulics, respiration,
`veg_energy_diagnostic` — there would finally be work to spread. The per-*patch* kernels (CAS box, soil
Thomas, the surface Newton) do not benefit at all: they stay at 12 lanes no matter how many cohorts there
are, and by Amdahl they set the floor.

**Wall 3 becomes survivable.** At 1,000 cohorts per patch a single `column_fast_step` costs ~2.7 ms
(§8a-bis), so a whole site-step is ~32 ms against 40–50 launches at 59 µs = 2.4–3.0 ms — **~9% overhead**
instead of the 4× penalty it is today, and ~1% at native-Linux launch latency.

**Wall 2 is unchanged.** In FP64 this GPU is still 0.48× the host. Fixing this requires either FP32 (§9)
or different hardware.

**Wall 4 is unchanged and becomes the binding constraint** — the allocatable-component state, the
in-loop `allocate`s, and the `error stop`s all have to go before a single per-cohort kernel can be
mapped. Two further problems get *worse* at scale:

- **Divergence.** At 2,000 cohorts some will be in hydraulic collapse (136 sub-steps) while their
  warp-mates take 1. §5c(v) measured that gap as **binary, not gradual** — it coincides exactly with the
  wood store pinning on its `1e-30` floor. Under SIMT the warp pays 136× regardless. This is not a
  tail case at high cohort count; it is a certainty.
- **The automatic-array stack ceiling.** The fast loop carries ~740 B of stack arrays per cohort live
  across four nested frames — a hard SIGSEGV near 10–12k cohorts on an 8 MB stack, and it scales *down*
  with `OMP_STACKSIZE` once patches are threaded (`CMakeLists.txt`, the stack-ceiling note). Thousands
  of cohorts already stresses the *host* code before any GPU work starts.

### 8c. The cheaper move, and it is not the GPU

If cohort counts go to thousands, the cohort axis stops being a rounding error and becomes the dominant
dimension — but it is currently **serial within a patch**. The patch axis caps host parallelism at 12;
the cohort axis does not cap at anything.

The obvious first step is therefore to **thread and/or SIMD-vectorise the cohort axis on the CPU**:

- the `elemental pure` / `*_batch` convention (§11) already set this up — the per-cohort physiology set
  is complete, and `elemental pure` is exactly the compiler's licence to vectorise;
- AVX-512 gives 8 FP64 lanes for free on this host, and FP32 would give 16 (§9);
- the §7 order-preserving reduction machinery already exists and generalises from the patch axis to a
  cohort axis;
- it is worth roughly the same 4–8× as a GPU port would be *after* walls 2 and 4 were paid down, at a
  fraction of the cost and with no numerical-scheme change.

**Recommendation for the high-cohort scenario: do the CPU cohort axis first, and re-measure.** A GPU only
becomes interesting after that is exhausted *and* the FP32 question is settled.

---

### 8d. Related question — could the GPU just be used as *N slow cores*, one patch per SM?

**No. MEASURED: at 20 patches the GPU is 26× SLOWER than 4 CPU cores, and the one-patch-per-SM
formulation never wins at any size.**

This is the "forget vectorisation, the device has 20 SMs and I have 4 cores" idea, and it deserves its
own measurement because it needs none of BB2's kernel restructuring — just the existing patch loop
retargeted. `scratchpad/patchproxy.f90` runs one *proxy* column solve per work item: serial FP64
adaptive march, data-dependent accept/reject branches, `exp`/`sqrt` transcendentals — the shape of code
that has no ILP to extract and no vector width to use. Two device formulations are compared against the
host:

- **`1 thr/team`** — `num_teams(np) thread_limit(1)`, i.e. literally one patch per SM;
- **`warp-packed`** — `target teams distribute parallel do`, i.e. patches spread across threads, 32 to a
  warp. This is what a realistic patch-loop offload would compile to.

Wall time for `np` items, µs (min over 5 reps):

| np | HOST 1 thr | HOST 4 thr | HOST 8 thr | DEVICE 1 thr/team | DEVICE warp-packed |
|---|---|---|---|---|---|
| **20** | 212 | **53** | 30 | 1409 (**26× slower**) | 1392 (**26× slower**) |
| 80 | 745 | 241 | 121 | 2406 | 1891 (7.8× slower) |
| 320 | 3067 | 748 | 408 | 8292 | 2326 (3.1× slower) |
| 1280 | 13612 | 3167 | 2119 | 28161 | **1434 (2.2× faster)** |
| 5120 | 50307 | 12265 | 8177 | 72989 | **2969 (4.1× faster)** |
| 20480 | 200577 | 47206 | 29652 | 290128 | **10121 (4.7× faster)** |

Three conclusions:

1. **One GPU lane is ~140× slower than one CPU core for this code.** At np = 20 all items run
   concurrently, so the ~1400 µs wall is one item's *latency* — against 9.8 µs on a single CPU core. A
   CUDA core is a 1.5 GHz in-order scalar lane with a fraction of the register file and no branch
   prediction; a Tiger Lake core is a wide out-of-order machine. Twenty of the former do not beat four
   of the latter.
2. **The literal one-patch-per-SM formulation never wins.** It is 26× slower at np = 20 and still 6×
   slower at np = 20,480, because `thread_limit(1)` uses one lane of every 32-wide warp and throws away
   97% of the device.
3. **The warp-packed formulation crosses over at roughly np ≈ 600–700** against 4 cores (≈ 900–1000
   against 8 threads), and its asymptotic win on this GPU is only **~4–5×**.

**And this proxy is generous to the GPU in every respect that matters.** All work items have nearly
identical trip counts, so there is essentially no warp divergence — where the real fast loop measures
1 → 200 sub-steps per cohort (§8a). There is no allocatable state, no heap traffic, no derived types,
a two-array map, and a register footprint small enough for full occupancy — where `column_fast_step`
has all of those problems (§7). The real crossover is therefore well above 1,000 concurrent columns, and
the real asymptotic win well below 4×.

**So: raising `max_patch` from 12 to 20 does not open a GPU path.** It is a good move for the *CPU* —
20 patches over 4 threads measured a clean 4.0× here — and that is the win available at that size.

## 9. Does FP32 change it?

**Yes — it removes wall 2 outright, and it is the single largest lever available on consumer NVIDIA
hardware. But it is a mixed-precision engineering program, not a flag flip.**

### 9a. Measured

Same 8-chain FMA kernel, `real32` instead of `real64`:

| | FP64 | FP32 | ratio |
|---|---|---|---|
| RTX 3050 Ti Laptop, `-mp=gpu` | 85.0 GFLOP/s | **4,780 GFLOP/s** | **56×** |
| i7-11370H, ifx `-O3 -xHost`, 4–8 threads | 177.6 GFLOP/s | 303.0 GFLOP/s | 1.7× |
| **device / host** | **0.48×** | **16×** | — |

FP32 flips the GPU from *half* the host's throughput to **sixteen times** it. Nothing else in this
document moves the needle that far. It also roughly doubles CPU SIMD throughput and halves memory
traffic, so it is worth evaluating **independently of whether a GPU is ever used.**

### 9b. MEDS is unusually cheap to *test* this on

Precision is a single parameter — `wp` in `meds_kinds.f90:13` — followed by 3,053 sites across `src/`.
Only **67** hard-coded `real64`/`c_double` occurrences exist outside it, confined to 7 files:
`meds_netcdf_c.f90`, `meds_io.f90`, `meds_output_stream.f90`, `meds_met_driver.f90`,
`meds_demography_capi.f90`, `meds_plant_capi.f90`, `meds_core_state_types.f90` — i.e. the I/O, forcing,
and C-API boundaries, which should stay FP64 anyway.

**A `wp = real32` build is therefore roughly a day's work, and it is decisive.** Build it, run the
37-test suite, run a 50-year spin-up, and see what breaks. That experiment should be done before any
further GPU discussion.

### 9c. What will break, and why uniform FP32 is the wrong end state

FP32 carries ~7 decimal digits (eps ≈ 6e-8). Ranked by expected severity for MEDS specifically:

1. **Conservation ledgers keyed to machine precision.** The whole-column energy/water/carbon budgets
   close to ~1e-13 and `rh_seam_gap` is asserted at that level. In FP32 the *floor* is 6e-8. These
   assertions would need re-tuning by six orders of magnitude — and the L2 strict mode (§4) is built on
   them. This is a test-and-tolerance problem rather than a physics one, but it is the largest diff.
2. **Phase change near freezing.** The enthalpy formulation differences large numbers
   (`tsupercool_liq = 56.79 K`, `u_liq ≈ 1.0 MJ/kg` against sensible-heat increments). Catastrophic
   cancellation near the freezing point is the classic FP32 failure mode in land models, and MEDS's snow
   and soil-energy closure both live there.
3. **The hydraulics dynamic range.** `psi_from_rwc` floors RWC at `1e-4` and returns −1.0e4 MPa; water
   masses run down to `1e-30 kg` (representable in FP32, but barely — min normal is 1.2e-38). §5c(v)–(vi)
   showed the sub-step count and the resulting `W_leaf` are **binary** in exactly this region: FP32 noise
   would flip cohorts across that switch non-deterministically.
4. **`matrix_exp` scaling-and-squaring.** Repeated squaring amplifies rounding; the CENTURY biogeochem
   and the hydraulic stage map both use it.
5. **Long accumulations.** A 50-year spin-up makes 18,250 daily additions to the carbon pools and the
   patch-area bookkeeping. FP32 random-walk drift is ~1e-5 relative, worst case ~1e-3 — comparable to
   `conservation_tol = 0.001` itself.
6. **The adaptive controller.** At `ark_rtol = 1e-3` the embedded difference carries ~4.5 significant
   digits in FP32. Workable, but expect a noisier controller and more rejections.

Items 1, 2, 3 and 5 all point the same way: **the conserved accumulators, the patch-area bookkeeping,
and the enthalpy/phase-change path should stay FP64**, while the per-cohort physiology — photosynthesis,
respiration, boundary-layer conductances, radiative transfer — can move to FP32, since its *inputs* are
known to two or three digits at best. That is mixed precision, and it is the well-trodden path
(ECMWF's single-precision IFS, ICON's mixed-precision dynamics). It is also more work than a uniform
switch, and it should be justified by the CPU win (§9a) before the GPU is invoked.

---

## 10. Where a GPU could actually pay: the regional axis

The only viable GPU axis for a demographic land model is **many independent columns**, not one column
solved faster:

| workload | concurrent columns | verdict |
|---|---|---|
| single site (today) | 12 | hopeless — 0.04% device fill |
| single site, fusion off, 1000+ cohorts | ~10⁴ per-cohort items | fillable, but §8b's walls 2/4 remain |
| **gridded regional, 10⁴–10⁵ cells** | **10⁵–10⁶** | **the real target** |
| calibration / parameter ensembles, 10³–10⁴ members | 10⁴–10⁵ | viable, but also trivially MPI-parallel across CPU nodes |

MEDS has **no site or ensemble axis at all** — `site_t` is one instance, and the `pheno_tair_sum` comment
says as much (*"single-site forcing"*). Adding one is a large architectural change touching state, I/O,
forcing and the driver — and it would give a real speedup on plain CPU MPI *immediately*, before any
device work.

Even with the site axis, the prerequisites are the full BB2 + BB3 program **plus** FP32-or-datacenter
hardware. On an A100 (FP64 ≈ 9.7 TFLOP/s) against a 64-core server CPU (~2–4 TFLOP/s) the honest ceiling
is ~3× in FP64, or considerably more in FP32 — against a model whose single-site 50-year spin-up already
runs in **545 s on 4 cores**.

---

## 11. What survives from the BB1/BB2/BB3 track

`MEDS_NUMERICS_SCOPING.md` §7 committed BB1 and BB2 and staged BB3. Re-read against measurement:

- **BB1 — done, and was already largely moot.** Phase 1 (allocation-churn removal) shipped and was
  valuable. Phase 2's own note is the key finding and is correct: *"pointer-aliasing into the global SoA
  would not have fixed it… the real prerequisite is a bare-array calling convention."* The persistent
  reservoirs were **already** site-wide flat SoA before the pass.
- **§11's bare-array / `elemental pure` convention — keep, but re-justify.** It is a genuine win for the
  Python C-API and for host SIMD (and it is the enabler for §8c's cohort axis). It is *not* a GPU
  enabler on its own. Justify it on those grounds.
- **BB2 as scoped — ❌ refuted** by §4–§7. It targets `column_fast_step` on a 12-patch site; the launch
  arithmetic alone kills it. Revisit only under §8's high-cohort or §10's regional scenario, and only
  after FP32 is settled.
- **BB3 — deferred with BB2.** Its warp-uniform fixed-count profile is the right answer to §7's
  divergence problem, but it purchases warp-uniformity by abandoning the adaptive controller. That is a
  numerics regression to be paid only for a device that is actually winning.
- **`MEDS_GPU=multicore` — keep.** Its real value is as the issue-#7 nvfortran array-temp trap detector
  and a portability gate. It costs nothing and it has caught real bugs.

---

## 12. Recommendations

1. **Do not invest in GPU offload for the current single-site model.** Mark BB2/BB3 refuted-by-measurement
   in `MEDS_NUMERICS_SCOPING.md` §7, the way E1/E2/N2e were in the integrator plan.
2. **Correct the two overselling statements** — `CMakeLists.txt:6` and `CLAUDE.md:76` both imply a GPU
   path that exists as one cold kernel at 0.4% occupancy.
3. **Keep `MEDS_GPU=multicore`** as the portability gate. `MEDS_GPU=gpu` still passes 37/37, so there is
   no urgency to remove it; demote it in the docs to a curiosity rather than a supported path.
4. **The largest untaken CPU lever is the ~24% allocator self-time** (§5c(iv)), driven by the ~26
   `allocate` calls per patch per step in `build_column_frozen`. E2 was only ever priced against the ARK
   *combinators* (7%); the rest of the 24% lives in `build_column_frozen`'s own locals and has never been
   attacked. Bigger, cheaper, and lower-risk than anything on the device.
5. **If cohort counts are going to thousands:**
   a. First get a trustworthy cost law — a `max_cohort = 0` (fusion-off) run of the real driver, long
      enough to accumulate thousands of cohorts naturally. §8a shows a subdivided synthetic fixture
      cannot substitute for it.
   b. Then **build the CPU cohort axis** (§8c) — threading and/or `elemental` vectorisation over cohorts,
      which is where 99% of the work lands at that scale, and which needs no numerical-scheme change.
   c. Only then re-open the device question. And expect issue **#104** (the ψ clamp artefact and its
      200× sub-step cliff) to be a hard prerequisite: at thousands of cohorts it stops being a tail case.
6. **Run the `wp = real32` experiment** (§9b) — one line plus ~67 boundary sites, decisive either way,
   and valuable for the CPU regardless of the GPU verdict.
7. **If gridded/regional MEDS is on the roadmap, that decision — not GPU work — is the gate** (§10). The
   site axis pays off on CPU MPI immediately and is the only thing that ever makes a GPU worth revisiting.

---

## 13. Reproducing these measurements

```bash
export PATH=$HOME/opt/nvidia/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
cmake -S . -B build-nv-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release \
      -DMEDS_GPU=gpu -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-nv-gpu -j8

# prove what actually offloads, and at what occupancy
NVCOMPILER_ACC_NOTIFY=3 ./build-nv-gpu/meds_main <cfg>.toml 2>&1 | grep "launch CUDA"
cuobjdump --list-elf ./build-nv-gpu/meds_main
```

Probes used (scratchpad, not tracked):

| probe | measures |
|---|---|
| `fp64peak8.f90` / `fp32peak8.f90` | achievable FP64 / FP32 throughput, device vs host |
| `launchcost.f90` | `!$omp target` launch + `map` cost at MEDS-scale `n` |
| `cudalat.cu` | raw CUDA launch latency (separates WSL2 from OpenMP-runtime overhead) |
| `devq.cu` | SM count, resident-thread capacity |
| `meas_ncoh.f90` | fast-loop cost scaling vs cohort count — **valid only to n ≈ 30**, see §8a |
| `patchproxy.f90` | per-lane serial throughput, host core vs GPU lane; CPU/GPU crossover in work-item count (§8d) |

The single-year benchmark configuration restarts from `examples/example_biophysics/spinup-S-20740701000000.nc`
with `end_time = "2075-07-01"`, `dt_fast = "900s"`, and all diagnostic output disabled.
