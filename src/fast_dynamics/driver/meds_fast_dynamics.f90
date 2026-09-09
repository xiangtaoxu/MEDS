!==========================================================================================!
! meds_fast_dynamics -- the per-SITE fast-biophysics TIER ORCHESTRATOR (the fast-cadence peer   !
! of meds_vegetation_dynamics; advance_one_step fans out to both). It owns the ORCHESTRATION      !
! only: for each patch it gathers the demographic cohort slice into the column_cohort_t buffer,   !
! assembles a patch_biophys_t working bundle from the state-hub-owned per-patch reservoirs         !
! (cas/soil_e/soil_w), runs n_fast_per_slow operator-split sweeps of the per-patch kernel           !
! column_fast_step, and writes the evolved reservoirs back to the site. The static                 !
! column_config_t + base met arrive via fast_context_t (the caller builds them -- no model         !
! parameters are hard-coded here); per-cohort leaf_temp/wood_temp/leaf_water_mass/wood_water_mass/    !
! leaf_surf_water/wood_surf_water are PERSISTED on the cohort block and adopted here each slow step   !
! (no reseeding) -- same as the soil + CAS reservoirs, this is genuine cross-slow-step memory.        !
!                                                                                          !
! The stepper calls fast_dynamics before the slow loop when cfg%fast_biophysics_on. This is the    !
! fast->slow seam's fast half; the daily-GPP handoff into carbon growth lands in a later step.     !
!==========================================================================================!
module meds_fast_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, umol_2_kgC, grav, cp_air, latent_heat_vap, day_sec
   use meds_config,           only : meds_config_t
   use meds_budget_check,     only : budget_t, budget_merge
   use meds_biogeochem_types, only : IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   use meds_therm_lib,           only : cas_enthalpy_of_temp, cas_temp_of_enthalpy, temp_to_internal_energy
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_column_view, only : copy_column_cohort
   use meds_fast_reconcile,  only : reconcile_tissue_water_capacity
   use meds_time,             only : meds_time_t, time_advance_seconds, time_to_string
   use meds_output_types,     only : output_manager_t, fast_sample_t
   use meds_site_diag_types,  only : N_CDIAG, patch_diag_block,                                  &
                                     PD_LE, PD_H, PD_RNET, PD_SW_IN, PD_SW_GROUND, PD_LW_GROUND, &
                                     PD_USTAR, PD_GGNET, PD_ROUGH, PD_DISPLACE, PD_CAS_TEMP,     &
                                     PD_CAS_SHV, PD_CAS_CO2, PD_GPP, PD_NEE, PD_TRANSP,          &
                                     PD_ROOT_UPTAKE, PD_INFILTRATION, PD_DRAINAGE, PD_RUNOFF,    &
                                     PD_PRECIP, PD_GROUND_TEMP, PD_RESID_ENERGY, PD_RESID_WATER, &
                                     cohort_diag_grow, cohort_diag_reset, patch_diag_grow,        &
                                     patch_diag_reset
   use meds_column_params, only : n_soil_layer_max, PSI_INIT, build_soil_hydr_params, build_soil_therm_params
   use meds_column_state_types, only : xi_accum_t, snow_column_t
   use meds_forcing_types,    only : met_driver_t, met_forcing_t
   use meds_met_driver,       only : met_advance, met_instant
   use meds_site_state_types, only : site_t, DMAX_PSI_LEAF_UNSET, DMAX_PSI_LEAF_ACCUM_RESET
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, ensure_aero_out_capacity, rad_pft_optics_t, &
                                 rad_forcing_t, rad_flux_t, alloc_rad_forcing, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW, &
                                 set_aero_env_atm, set_aero_env_canopy
   use meds_fast_types, only : patch_biophys_t, ensure_patch_biophys_capacity
   use meds_hydr_lib, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : snow_params_t
   use meds_optics_lib,       only : beta_params_from_mean
   use meds_canopy_types, only : ground_optics_state_t
   use meds_canopy_radiation, only : canopy_radiation, derive_rad_optics, ground_optics
   use meds_ground_biophysics, only : snow_cover_fraction
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t,        &
                                     GRP_THETA, GRP_SOIL_T,                                       &
                                     column_budget_t,                                             &
                                     ensure_column_cohort_capacity, apply_hydraulics_config
   use meds_fast_step,       only : column_fast_step
   use meds_hydr_lib,         only : water_content, clamp_water_to_capacity
   !$ use omp_lib,            only : omp_get_thread_num
   implicit none
   private

   public :: fast_context_t, init_fast_reservoirs, fast_dynamics, build_fast_context

   !----- Absorbed-PAR (VIS) energy -> photon-flux conversion [umol photon / J], the 400-700 nm !
   !      value (~4.57). Used ONLY on the RT path (true absorbed PAR); the const path keeps the    !
   !      2.1 total-SW blend (column_forcing_t default) -- see fast_dynamics.                       !
   real(wp), parameter :: PAR_W_2_UMOL = 4.6_wp

   !----- §7 C3 (deterministic reductions). The site-level fast-loop accumulators are written once !
   !      per (patch, sub-step) and folded into site%... afterwards. Staging them in ONE array,      !
   !      indexed by the slots below, lets the patch loop write DISJOINT cells (no reduction clause,  !
   !      no atomics) and lets the fold run in the ORIGINAL (patch, sub-step) order. That order is     !
   !      the whole point: OpenMP `reduction(+:)` sums thread partials in thread-ARRIVAL order, so the  !
   !      last bits of every site diagnostic would move with the thread count -- and byte-identical     !
   !      netCDF across thread counts is §7's hard requirement, not a nicety. Folding in patch order    !
   !      reproduces the serial fold EXACTLY, so these results are also unchanged from before §7. ------!
   integer(ik), parameter :: RED_ET            =  1_ik   !< site ET [kg/m2], area-weighted
   integer(ik), parameter :: RED_PHENO_TAIR    =  2_ik   !< daily-mean air-temperature accumulator [K]
   integer(ik), parameter :: RED_INTEG_STEPS   =  3_ik   !< integrator work: accepted steps
   integer(ik), parameter :: RED_INTEG_REJ     =  4_ik   !< integrator work: rejected steps
   integer(ik), parameter :: RED_SOIL_NSUB     =  5_ik   !< soil-water sub-steps
   integer(ik), parameter :: RED_HYDRO_NSUB    =  6_ik   !< plant-hydraulics sub-steps
   integer(ik), parameter :: RED_HYDRO_THRASH  =  7_ik   !< pathological hydraulics sub-stepping (#104)
   integer(ik), parameter :: RED_NONCONV       =  8_ik   !< non-converged hydraulic solves
   integer(ik), parameter :: RED_RK45_RESCUE   =  9_ik   !< steps the RK45 path handed back to split
   integer(ik), parameter :: RED_CLAMP_STAGE   = 10_ik   !< stage-level state clamps
   integer(ik), parameter :: RED_CLAMP_COMMIT  = 11_ik   !< commit-level state clamps
   integer(ik), parameter :: RED_CLAMP_MASS    = 12_ik   !< water mass created by clamping [kg/m2]
   integer(ik), parameter :: RED_CLAMP_ENERGY  = 13_ik   !< energy created by clamping [J/m2]
   integer(ik), parameter :: N_RED             = 13_ik

   !----- Everything the fast driver needs beyond the site + cfg: the static column config plus !
   !      the reference met + initial soil state. The CALLER builds this (from TOML in the       !
   !      production path; the MVP holds constant, horizontally-uniform boundary conditions).    !
   type :: fast_context_t
      type(column_config_t) :: col_config                 !< soil/thermal/hydro/aero/resp column config
      real(wp) :: u_ref   = 2.0_wp, zref = 30.0_wp  !< [m/s],[m] reference wind + height
      real(wp) :: press   = 101325.0_wp             !< [Pa]
      real(wp) :: rho_air = 1.2_wp                  !< [kg/m3]
      real(wp) :: air_temp = 288.0_wp               !< [K]        reference-level air temperature
      real(wp) :: shv_atm  = 0.008_wp               !< [kg/kg]    reference-level specific humidity
      real(wp) :: co2_atm  = 400.0_wp               !< [umol/mol] free-atmosphere CO2
      real(wp) :: rad_sw_top    = 400.0_wp          !< [W/m2] shortwave into the canopy (leaves)
      real(wp) :: rad_sw_ground = 60.0_wp           !< [W/m2] shortwave reaching the ground
      real(wp) :: rainfall        = 0.0_wp            !< [kg/m2/s] ground-reaching rainfall
      real(wp) :: snowfall         = 0.0_wp            !< [kg/m2/s] frozen rainfall (snowfall)
      real(wp) :: theta_init      = 0.30_wp         !< [m3/m3] initial soil moisture (all layers)
      real(wp) :: soil_temp_init  = 288.0_wp        !< [K]     initial soil + CAS temperature
      real(wp) :: veg_height_bare = 1.0_wp          !< [m] canopy height for a cohort-free patch
      !----- Canopy-RT optics (per-PFT spectral/angle table + ground surface), for the RT-driven   !
      !      per-cohort absorbed SW/PAR on the forcing path (§6.3). MVP placeholders, built once.    !
      type(rad_pft_optics_t) :: rad_opt             !< per-PFT, per-band leaf/wood scattering + LIDF
      real(wp) :: soil_albedo(3) = [0.15_wp, 0.30_wp, 0.0_wp]  !< [VIS,NIR,LW] ground albedo (LW=0)
      real(wp) :: soil_emiss     = 0.95_wp          !< [-] ground longwave emissivity
   end type fast_context_t

contains

   !=======================================================================================!
   !  Build the fast-loop context for the production run. The scalar reference met + reservoir  !
   !  seeds keep the fast_context_t defaults (constant, horizontally-uniform MVP boundary          !
   !  conditions). The static column config `col_config` is assembled here from DOCUMENTED MVP            !
   !  PLACEHOLDER parameters (soil texture/thermal, stem/root maintenance-respiration factors,      !
   !  heterotrophic-Rh and plant-hydraulics constants) -- there is no per-site column_config TOML   !
   !  loader yet, so these mirror the validated test_fast_loop configuration (they pass the fast-   !
   !  loop energy/water/CO2 budget checks). The fast loop is OPT-IN (cfg%fast_biophysics_on,         !
   !  default .false.), so this leaves the default demographic run unchanged. Follow-up: source      !
   !  these from a [column]/[soil] TOML block + per-PFT respiration traits.                          !
   !=======================================================================================!
   subroutine build_fast_context(cfg, ctx)
      type(meds_config_t),  intent(in)  :: cfg
      type(fast_context_t), intent(out) :: ctx
      !----- The physical soil column, from [soil_column]. These were eleven hydraulic and three   !
      !      thermal LITERALS here, so the column's depth, layer count and texture could not be    !
      !      changed without a recompile -- and `depth` is the value the docs flag as shallower     !
      !      than the annual thermal damping depth, which made that defect unreachable from a       !
      !      config file. Defaults reproduce the old literals exactly. --------------------------!
      associate (sc => cfg%soil_column)
         call build_soil_hydr_params(sc%n_layer, sc%retention, sc%depth, sc%grid_growth,           &
                                     sc%theta_sat, sc%theta_res, sc%ksat, sc%curve_par_a,          &
                                     sc%curve_par_n, sc%root_beta, sc%psi_fc, ctx%col_config%soil)
         call build_soil_therm_params(sc%n_layer, sc%solid_conductivity, sc%dry_conductivity,      &
                                      sc%dry_heat_capacity, ctx%col_config%soil_thermal)
      end associate
      !----- Autotrophic maintenance-respiration + heterotrophic-Rh + prescribed soil-C pool. ---!
      ctx%col_config%wood%is_woody = .true.
      ctx%col_config%wood%stem_resp_factor25 = 0.06_wp
      ctx%col_config%root%root_resp_factor25 = 0.30_wp
      !----- Plant hydraulics: flatten the [hydraulics] config into hydraulics_params + rhizo_cond and build   !
      !       the vulnerability lookup table (dormant at kexp=2; consulted only if wood_kexp leaves    !
      !       {1,2}). Values come from cfg (MVP defaults unless a [hydraulics] block overrides). ------!
      call apply_hydraulics_config(cfg%hydraulics, ctx%col_config%hydraulics_params)
      call build_leaf_photo_table(cfg, ctx%col_config%leaf_photo)    ! per-PFT leaf parameters, once per run
      ctx%col_config%specific_root_area = cfg%hydraulics%specific_root_area
      !----- P3 coupled-surface (Picard) solver knobs + option selectors, from the [fast] block. --!
      ctx%col_config%canopy_water_on    = cfg%canopy_water_on
      !----- Fast-loop biophysics run-config from the [soil]/[energy]/[snow]/[aerodynamics] blocks   !
      !      (all opt-in; cfg carries the meds_biophysics_opts defaults unless a block overrides).    !
      !      Same types as the column config members, so a plain verbatim struct copy. --------------!
      ctx%col_config%soil_water_opts  = cfg%soil        ! [soil]         -> soil-water Richards solver opts
      ctx%col_config%energy = cfg%energy      ! [energy]       -> soil-thermal solver opts
      ctx%col_config%snow   = cfg%snow        ! [snow]         -> snow physical parameter table
      ctx%col_config%aero   = cfg%aero        ! [aerodynamics] -> canopy-aerodynamics constants

      !----- §8c Layer 1: ONE tolerance source drives the whole fast-loop hierarchy. build_tol_set     !
      !      SEEDS each group from the setting that governs it today, so these pushes are the           !
      !      IDENTITY by default (byte-identical); when [fast].rtol_all > 0 the single master dial       !
      !      propagates into every nested sub-solver as well as the ARK/RK45 march. hydraulics_opts (the plant-  !
      !      hydraulics sub-solver's OWN adaptive step-doubling tolerance) is NOT pushed from here any     !
      !      more (MEDS_ED2_RK45_DESIGN.md sec 4/6, P2): it still operates in PSI space internally         !
      !      (solve_plant_water's own matrix-exponential sub-stepping), and the retired GRP_PSI's outer     !
      !      WRMS group is now GRP_LEAF_W/GRP_WOOD_W in MASS units [kg/plant] -- feeding an MPa-space        !
      !      tolerance from a kg/plant-space group would be a unit mismatch, not a unification. hydraulics_opts    !
      !      keeps its own type default (rtol=atol=1e-3), unchanged from what it used implicitly before.  !
      !------------------------------------------------------------------------------------------------!
      ctx%col_config%integrator = build_integrator_opts(cfg)
      associate (tols => ctx%col_config%integrator%error_control%tols)
         ctx%col_config%soil_water_opts%rtol   = tols%rtol(GRP_THETA)
         ctx%col_config%soil_water_opts%atol   = tols%atol(GRP_THETA)
         ctx%col_config%energy%rtol  = tols%rtol(GRP_SOIL_T) ; ctx%col_config%energy%atol  = tols%atol(GRP_SOIL_T)
      end associate

      !----- §5.1 process mask: config logicals -> the mask the schemes honor. All-on = full column. --!
      ctx%col_config%mask%veg_energy = cfg%mask_veg_energy
      ctx%col_config%mask%cas_energy = cfg%mask_cas_energy
      ctx%col_config%mask%cas_vapour = cfg%mask_cas_vapour
      ctx%col_config%mask%cas_co2    = cfg%mask_cas_co2
      ctx%col_config%mask%soil_heat  = cfg%mask_soil_heat
      ctx%col_config%mask%soil_water = cfg%mask_soil_water
      ctx%col_config%mask%hydraulics = cfg%mask_hydraulics

      !----- Canopy-RT optics table (MVP placeholders; PFT-UNIFORM -- optics do not vary by PFT   !
      !      yet, that is the Phase-2 [radiation] PFT-TOML block). Values mirror                    !
      !      test_canopy_radiation.f90 build_optics. Built ONCE; read-only downstream.               !
      block
         integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT
         integer(ik) :: np
         real(wp), allocatable :: rl(:,:), tl(:,:), rw(:,:), tw(:,:), cl(:), cw(:), bp(:), bq(:)
         logical  :: hb(NB), he(NB)
         real(wp) :: bpp, bqq
         np = cfg%pft%n
         allocate(rl(NB,np), tl(NB,np), rw(NB,np), tw(NB,np), cl(np), cw(np), bp(np), bq(np))
         rl(RAD_VIS,:) = 0.10_wp ; tl(RAD_VIS,:) = 0.05_wp        ! leaf VIS reflect/transmit
         rl(RAD_NIR,:) = 0.45_wp ; tl(RAD_NIR,:) = 0.25_wp        ! leaf NIR
         rl(RAD_LW,:)  = 0.03_wp ; tl(RAD_LW,:)  = 0.0_wp         ! leaf_emiss = 0.97
         rw(RAD_VIS,:) = 0.11_wp ; tw(RAD_VIS,:) = 0.001_wp       ! wood VIS (near-opaque)
         rw(RAD_NIR,:) = 0.25_wp ; tw(RAD_NIR,:) = 0.001_wp       ! wood NIR
         rw(RAD_LW,:)  = 0.10_wp ; tw(RAD_LW,:)  = 0.0_wp         ! wood_emiss = 0.90
         cl = 0.80_wp ; cw = 0.50_wp                              ! leaf/wood clumping
         call beta_params_from_mean(45.0_wp, 20.0_wp, bpp, bqq)   ! mean 45deg, std 20deg leaf angle
         bp = bpp ; bq = bqq
         hb = [.true.,  .true.,  .false.]                         ! VIS/NIR have a beam; LW does not
         he = [.false., .false., .true. ]                         ! only LW emits
         call derive_rad_optics(NB, np, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, ctx%rad_opt)
      end block
   end subroutine build_fast_context

   !----- Seed every patch's fast reservoirs to a horizontally-uniform equilibrium. Called once !
   !      after the community is built (production) or before a fast-loop test. --------------!
   subroutine init_fast_reservoirs(site, ctx)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      integer(ik) :: ip, k, nsl
      nsl = ctx%col_config%soil%n_active
      do ip = 1_ik, site%patch%n
         associate (cas => site%patch%cas(ip), se => site%patch%soil_e(ip), sw => site%patch%soil_w(ip))
            sw%theta(1:nsl)  = ctx%theta_init ; sw%w_surface = 0.0_wp
            sw%w_surface_enth = 0.0_wp        ! dry pond -> zero enthalpy (issue #78 item 4)
            do k = 1_ik, nsl
               se%soil_energy(k) = temp_to_internal_energy(ctx%col_config%soil_thermal%soil_dry_heat_capacity(k),    &
                                   ctx%theta_init * rho_h2o, ctx%soil_temp_init, 1.0_wp)
               se%soil_temp(k)   = ctx%soil_temp_init ; se%soil_fliq(k) = 1.0_wp
            end do
            cas%can_shv      = ctx%shv_atm ; cas%can_co2 = ctx%co2_atm
            cas%can_enthalpy = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
            cas%can_temp     = ctx%air_temp
         end associate
      end do
   end subroutine init_fast_reservoirs

   !=======================================================================================!
   !  Advance every patch of the site by one slow-step's worth of fast biophysics: n_fast_per_ !
   !  slow operator-split sweeps over the state-hub reservoirs, on constant (MVP) forcing.     !
   !  Optional out-args report the worst whole-column budget residuals + the fail count so a    !
   !  caller/test can assert conservation.                                                     !
   !=======================================================================================!
   subroutine fast_dynamics(site, ctx, cfg, met_drv, step_start, worst_energy, worst_water, &
                            n_budget_fail, mgr, run_energy_budget, run_water_budget)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      type(meds_config_t),  intent(in)    :: cfg
      type(met_driver_t), optional, intent(inout) :: met_drv       !< live forcing reader (per-sub-step met)
      type(meds_time_t),  optional, intent(in)    :: step_start    !< calendar time at the START of this slow step
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail
      type(output_manager_t), optional, intent(inout) :: mgr       !< FAST-tier staging (filled when present + on)
      !----- RUN-level whole-column ledgers: this slow step's per-patch accumulators are area-      !
      !      weighted into a site accumulator and folded in here, so a caller that keeps them across !
      !      the whole run can report the SIGNED cumulative residual (a one-signed bias below the     !
      !      per-step tolerance is invisible to worst_* and n_budget_fail). -------------------------!
      type(budget_t), optional, intent(inout) :: run_energy_budget, run_water_budget

      !----- §7 C1: the n_fast_per_slow met samples + their sample TIMES, precomputed ONCE per slow  !
      !      step. `t_sub` depends only on `isub`, so met_advance (a FILE READER -- it may reload a  !
      !      netCDF bracket) and met_instant (solar geometry + Weiss-Norman disaggregation) were     !
      !      doing site-uniform work n_patch times over. Hoisting them is a win on its own AND is    !
      !      what takes the file reader out of the region §7 makes parallel. -------------------------!
      type(met_forcing_t), allocatable :: met_sample(:)
      type(meds_time_t),   allocatable :: t_sample(:)
      !----- §7 C3 reduction staging: one cell per (accumulator, sub-step, patch), written exactly   !
      !      once inside the patch loop and folded into site%... in patch order afterwards. Sized     !
      !      per call (O(1) allocations per slow step, like BB1's scratch); never zeroed, because      !
      !      every cell is ASSIGNED, not accumulated into. --------------------------------------------!
      real(wp),    allocatable :: red_site(:,:,:)                  !< (N_RED, sub-step, patch)
      real(wp),    allocatable :: red_worst_energy(:), red_worst_water(:)   !< (patch); max-folded
      type(budget_t), allocatable :: red_budget_energy(:), red_budget_water(:) !< (patch); area-merged
      type(budget_t) :: site_energy_budget, site_water_budget
      integer(ik), allocatable :: red_nfail(:)                     !< (patch) budget-failure counts
      type(fast_sample_t), allocatable :: red_fast(:,:)            !< (sub-step, patch) FAST-tier staging
      real(wp),    allocatable :: red_fast_soil_temp(:,:,:)        !< (layer, sub-step, patch)
      real(wp),    allocatable :: red_fast_soil_water(:,:,:)       !< (layer, sub-step, patch)
      real(wp)    :: f_ground
      integer(ik) :: ip, isub, npatch, ncoh_max, nsub, nl, n_thread
      logical     :: do_forcing, do_fast, do_cdiag, do_pdiag
      !----- Per-(cohort, sub-step) diagnostic scratch: filled by the pre-pass through cdiag, then   !
      !      folded (dt-weighted) into site%cohort%diag. PER THREAD, because the patch loop is         !
      !      parallel -- it rides the same per-thread pool discipline as the rest of the scratch.  ---!
      real(wp),    allocatable :: cdiag_pool(:,:,:)                !< (N_CDIAG, cohort, thread)
      !=========================================================================================!
      !  §7 C2 -- the PER-THREAD scratch POOL. One copy of every per-patch working buffer per       !
      !  thread, indexed by thread id, exactly as the plan specifies. The patch loop then aliases    !
      !  its own slice with `associate`, so the loop body below is textually unchanged from the      !
      !  serial version and nothing in it can accidentally reach another thread's buffer.            !
      !                                                                                             !
      !  A POOL, and not the obvious OpenMP data-sharing clauses, because BOTH of those are broken   !
      !  for these types on ifx 2026 and each fails SILENTLY in its own way:                          !
      !    - declaring the scratch in a BLOCK inside the region: ifx does not default-initialise      !
      !      block-scoped derived types in an outlined region, so their allocatable descriptors hold  !
      !      garbage, ensure_*_capacity reads them as already-allocated, and the first deallocate      !
      !      SIGSEGVs -- reproduced at n_threads = 1;                                                  !
      !    - `private(...)`: ifx constructs each thread's default-initialised copy through a           !
      !      compiler-generated STATIC mold (`AERO_OUT_T.omp.mold_ctor`) that every thread WRITES;     !
      !    - `firstprivate(...)`: same story one step over, in `FAST_CONTEXT_T.omp.copy_ctor`.          !
      !  The last two are real data races (found with ThreadSanitizer) and they produce exactly the     !
      !  failure mode §7 exists to prevent: plausible-looking answers that move with the thread count.  !
      !  A plain array indexed by thread id involves no compiler-generated constructor at all.          !
      !                                                                                                !
      !  BB1's hoist is preserved: ensure_*_capacity runs ONCE PER THREAD here (not once per patch), so !
      !  allocations are O(n_thread) per slow step -- identical to BB1's O(1) at the default 1 thread.  !
      !=========================================================================================!
      type(column_cohort_t),  allocatable :: coh_pool(:)
      type(column_forcing_t), allocatable :: forc_pool(:)
      type(aero_env_t),       allocatable :: aenv_pool(:)
      type(aero_geom_t),      allocatable :: ageom_pool(:)
      type(aero_out_t),       allocatable :: aero_pool(:)
      type(patch_biophys_t),  allocatable :: bio_pool(:)
      type(column_budget_t),  allocatable :: budg_pool(:)
      type(fast_context_t),   allocatable :: ctx_pool(:)     !< per-sub-step met overlay on ctx
      type(met_forcing_t),    allocatable :: met_pool(:)
      real(wp),               allocatable :: gpp_pool(:,:), leaf_resp_pool(:,:)
      real(wp),               allocatable :: stem_resp_pool(:,:), root_resp_pool(:,:), psi_leaf_pool(:,:)
      real(wp)    :: sum_lai, le_flux, h_flux, rnet, gpp_patch, npp_patch, w_area, dt_fast_days
      integer(ik) :: j, i, i0, ncoh, ith

      !----- Live forcing drives the fast loop only when it is ON and a reader + step time are    !
      !      supplied; otherwise ctx_now stays == ctx and the loop runs the CONSTANT-forcing MVP    !
      !      bit-identically (the diurnal cycle lives INSIDE the sub-step loop, design §1.1/§6.2).  !
      do_forcing = cfg%forcing%forcing_on .and. present(met_drv) .and. present(step_start)
      f_ground = ctx%rad_sw_ground / max(ctx%rad_sw_top, tiny_num)   ! ground/canopy-top SW transmittance
      npatch   = site%patch%n
      nsub     = cfg%n_fast_per_slow
      n_thread = max(1_ik, cfg%n_threads)

      !----- Reset the fast->slow carbon accumulators (gross GPP + maintenance-resp losses) BEFORE !
      !      the fast window (compute_carbon_allocation reads them after; it has site intent(in),   !
      !      cannot reset).                                                                          !
      site%cohort%gpp_accum(1:site%cohort%n)       = 0.0_wp
      site%cohort%leaf_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%stem_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%root_resp_accum(1:site%cohort%n) = 0.0_wp
      !----- ROLL OVER the predawn water status (#95): today's accumulated maximum becomes the value  !
      !      the leaf kernel uses tomorrow, then the accumulator restarts. A cohort whose max is still !
      !      the UNSET sentinel (a recruit born mid-day) keeps it, so column_prepass seeds it from the  !
      !      soil rather than inheriting a meaningless value. -----------------------------------------!
      !----- Roll over only where the day actually produced a sample (a cohort born mid-day, or one  !
      !      culled before its first step, has none). The guard is against the RESET value, not 0:    !
      !      psi_leaf is <= 0, so '<= 0' would also accept an untouched accumulator. -----------------!
      where (site%cohort%dmax_psi_leaf_accum(1:site%cohort%n) > 0.5_wp * DMAX_PSI_LEAF_ACCUM_RESET)          &
         site%cohort%dmax_psi_leaf(1:site%cohort%n) = site%cohort%dmax_psi_leaf_accum(1:site%cohort%n)
      site%cohort%dmax_psi_leaf_accum(1:site%cohort%n) = DMAX_PSI_LEAF_ACCUM_RESET
      !----- Reset the daily fast->slow soil-carbon accumulator (B2; opt-in [soil_carbon].            !
      !      soil_carbon_on -- harmless no-op accumulation when off, since column_prepass leaves        !
      !      budget%xi_step/rh_matrix_step at 0 in that case). ------------------------------------------!
      if (cfg%soil_carbon_on) site%patch%xi_accum(1:site%patch%n) = xi_accum_t()
      !----- The site daily-mean air-temperature accumulator (which the slow-loop phenology driver    !
      !      reads AFTER this fast window), site ET, and the §5.3 integrator WORK counters are all     !
      !      reset and refilled by the §7 C3 fold AFTER the patch loop -- see the bottom of this        !
      !      routine. Air temperature is site-uniform, so accumulating once per (patch, sub-step) and   !
      !      dividing by the count still yields the daily mean. ---------------------------------------!

      !----- FAST (sub-daily) output staging: fill mgr%fast(:) only when the tier is active and a       !
      !      diurnal signal exists (forcing on). Lazily allocate the per-sub-step buffers ONCE (sizes    !
      !      are run-constant), then zero the accumulators for this slow step. Serialize stays in main.  !
      do_fast = present(mgr) .and. do_forcing
      if (do_fast) do_fast = mgr%enabled .and. mgr%reg%nidx(1) > 0_ik
      !----- The per-cohort / per-patch diagnostic capture runs only when the registry actually      !
      !      asked for it (main sets `active` from the live variable list). Everything downstream is  !
      !      gated on these two flags, so a run that reports no ecophysiology takes the original      !
      !      code path exactly.  --------------------------------------------------------------!
      do_cdiag = site%cohort%diag%active
      do_pdiag = site%patch%diag%active
      !----- RESET the diagnostic accumulators for this slow step, BEFORE the patch loop folds into  !
      !      them. Same lifecycle as gpp_accum / et_accum / xi_accum: one window per slow step, read  !
      !      once at the output tick.  ---------------------------------------------------------!
      if (do_cdiag) then
         call cohort_diag_grow(site%cohort%diag, max(site%cohort%n, 1_ik))
         call cohort_diag_reset(site%cohort%diag)
         site%cohort%diag%n = site%cohort%n
      end if
      if (do_pdiag) then
         call patch_diag_grow(site%patch%diag, max(npatch, 1_ik))
         call patch_diag_reset(site%patch%diag)
         site%patch%diag%n = npatch
      end if
      if (do_fast) then
         nl = n_soil_layer_max
         if (.not. allocated(mgr%fast)) then
            allocate(mgr%fast(nsub), mgr%fast_time(nsub))
            allocate(mgr%fast_soil_temp(nl, nsub), mgr%fast_soil_water(nl, nsub))
            allocate(mgr%fast_coh_ltemp(max(mgr%cohort_max,1_ik), nsub),                            &
                     mgr%fast_coh_gpp(max(mgr%cohort_max,1_ik), nsub),                              &
                     mgr%fast_coh_height(max(mgr%cohort_max,1_ik), nsub))
         end if
         mgr%n_fast_sub    = nsub
         mgr%fast_n_soil   = nl
         mgr%fast_n_cohort = site%cohort%n
         do isub = 1_ik, nsub
            mgr%fast(isub) = fast_sample_t()
         end do
         mgr%fast_soil_temp = 0.0_wp ; mgr%fast_soil_water = 0.0_wp
         mgr%fast_coh_ltemp = 0.0_wp ; mgr%fast_coh_gpp = 0.0_wp ; mgr%fast_coh_height = 0.0_wp
      end if

      !----- BB1 phase 1 (MEDS_NUMERICS_SCOPING.md sec 7/10.2): size the per-patch fast-loop        !
      !      scratch (col_cohort/biophys/aero/forc + the per-cohort output accumulators) to the SITE-WIDE MAX   !
      !      cohort count ONCE here, instead of once PER PATCH inside the loop below. Every reader    !
      !      (column_fast_step and this driver) loops by the ACTIVE count (col_cohort%n / ncoh), never by     !
      !      size(...) (verified across src/test), so reusing a larger patch's leftover capacity for   !
      !      a smaller one is bit-identical -- this only cuts O(n_patch) heap allocations per slow      !
      !      step down to O(1). The persistent reservoirs (site%patch%cas/soil_e/soil_w/snow, site%     !
      !      cohort%leaf_temp/wood_temp/leaf_water_mass/wood_water_mass/leaf_surf_water/                !
      !      wood_surf_water) are UNCHANGED by this -- they were already site-wide flat SoA, not         !
      !      per-patch scratch (already true of MEDS's arch). ------------------------------------------!
      !----- Reconcile stored tissue water against the capacity today's biomass allows, BEFORE     !
      !      anything reads it. Once per call, outside every loop: this used to run per cohort per  !
      !      sub-step inside the gather, and it WRITES, so an unbooked mass edit sat in the          !
      !      integrator's inner loop where the whole-column ledger could not see it. ----------------!
      call reconcile_tissue_water_capacity(site, cfg)

      ncoh_max = 0_ik
      do ip = 1_ik, npatch
         ncoh_max = max(ncoh_max, site%patch%cohort_count(ip))
      end do

      !----- Per-thread scratch for the per-cohort diagnostic capture. Sized from ncoh_max, so it   !
      !      MUST follow the loop just above that computes it.  --------------------------------!
      allocate(cdiag_pool(N_CDIAG, max(ncoh_max,1_ik), max(n_thread,1_ik)))

      !----- §7 C3: the per-(sub-step, patch) reduction staging. Allocated (not zeroed -- every cell  !
      !      is assigned) per call, so this adds O(1) allocations per slow step, not O(n_patch). -----!
      allocate(red_site(N_RED, nsub, max(npatch,1_ik)))
      allocate(red_budget_energy(max(npatch,1_ik)), red_budget_water(max(npatch,1_ik)))
      allocate(red_worst_energy(max(npatch,1_ik)), red_worst_water(max(npatch,1_ik)),               &
               red_nfail(max(npatch,1_ik)))
      if (do_fast) then
         allocate(red_fast(nsub, max(npatch,1_ik)))
         allocate(red_fast_soil_temp(nl, nsub, max(npatch,1_ik)),                                   &
                  red_fast_soil_water(nl, nsub, max(npatch,1_ik)))
      end if

      !----- §7 C1: precompute this slow step's met samples. met_advance is called here in strictly   !
      !      increasing t exactly once per sub-step, instead of being rewound through the same t       !
      !      sequence once per patch; met_instant is a pure function of (reader state, t), so the      !
      !      per-sub-step values are unchanged and the loop below is byte-identical. The sample TIMES  !
      !      are kept too: they are what the FAST-tier output stamps and what the probe writes, and     !
      !      after the hoist a single scalar `t_sub` would leave both reading the LAST sub-step's time. !
      if (do_forcing) then
         allocate(met_sample(nsub), t_sample(nsub))
         do isub = 1_ik, nsub
            t_sample(isub) = time_advance_seconds(step_start,                                          &
                       (real(isub, wp) - 1.0_wp + cfg%forcing_sample_frac) * cfg%dt_fast)
            call met_advance(met_drv, t_sample(isub))
            met_sample(isub) = met_instant(met_drv, t_sample(isub))
         end do
      end if
      !----- Site-uniform, so it belongs OUT of the patch loop (where every patch used to rewrite it  !
      !      with the same value -- benign serially, a data race once threaded). ---------------------!
      if (do_fast) mgr%fast_time(1:nsub) = t_sample(1:nsub)

      !=========================================================================================!
      !  §7 C2 -- the parallel patch loop. Patch columns are independent within a dt_fast (they    !
      !  couple only through the slow loop), so patches are the parallel axis; `schedule(dynamic,1)` !
      !  because the adaptive march makes them WILDLY unequal (§5c(v): one collapsed patch runs 13x  !
      !  the others), which a static schedule would serialize behind.                                 !
      !                                                                                               !
      !  Without OpenMP flags every directive is a comment and the pool is one deep, so this is         !
      !  literally the serial code that preceded it.                                                    !
      !=========================================================================================!
      allocate(coh_pool(n_thread), forc_pool(n_thread), aenv_pool(n_thread), ageom_pool(n_thread),  &
               aero_pool(n_thread), bio_pool(n_thread), budg_pool(n_thread), ctx_pool(n_thread),    &
               met_pool(n_thread))
      allocate(gpp_pool(ncoh_max, n_thread), leaf_resp_pool(ncoh_max, n_thread),                    &
               stem_resp_pool(ncoh_max, n_thread), root_resp_pool(ncoh_max, n_thread),              &
               psi_leaf_pool(ncoh_max, n_thread))
      do ith = 1_ik, n_thread
         ctx_pool(ith) = ctx
         call ensure_column_cohort_capacity(coh_pool(ith), ncoh_max)
         call ensure_patch_biophys_capacity(bio_pool(ith), ncoh_max, ctx%air_temp, ctx%shv_atm,     &
                                            ctx%co2_atm, ctx%air_temp)
         call ensure_aero_out_capacity(aero_pool(ith), ncoh_max)
         call alloc_forcing(forc_pool(ith), ncoh_max)
      end do

      !$omp parallel do default(shared) schedule(dynamic, 1) num_threads(n_thread)                  &
      !$omp    private(ip, ith, isub, j, i, i0, ncoh,                                              &
      !$omp            sum_lai, le_flux, h_flux, rnet, gpp_patch, npp_patch, w_area, dt_fast_days)
      do ip = 1_ik, npatch
         !----- This thread's slot in the scratch pool. The `!$` sentinel keeps the non-OpenMP build  !
         !      on slot 1 with no dependence on omp_lib. ---------------------------------------------!
         ith = 1_ik
         !$ ith = int(omp_get_thread_num(), ik) + 1_ik
         associate (col_cohort           => coh_pool(ith),        forc          => forc_pool(ith),         &
                    aenv          => aenv_pool(ith),       ageom         => ageom_pool(ith),        &
                    aero          => aero_pool(ith),       biophys           => bio_pool(ith),          &
                    budget          => budg_pool(ith),       ctx_now       => ctx_pool(ith),          &
                    met           => met_pool(ith),        gpp_coh       => gpp_pool(:,ith),        &
                    leaf_resp_coh => leaf_resp_pool(:,ith), stem_resp_coh => stem_resp_pool(:,ith), &
                    root_resp_coh => root_resp_pool(:,ith), psi_leaf_coh => psi_leaf_pool(:,ith), &
                    cdiag_buf     => cdiag_pool(:,:,ith))
         ncoh = site%patch%cohort_count(ip)
         i0   = site%patch%cohort_offset(ip)

         !----- Gather the patch's cohort slice into the column buffer (+ MVP derived inputs).     !
         !      Capacity was ensured above (ncoh <= ncoh_max always); this just updates the ACTIVE   !
         !      count -- no allocation. -----------------------------------------------------------!
         call copy_column_cohort(col_cohort, site%cohort, i0, ncoh)
         sum_lai = 0.0_wp
         do j = 1_ik, ncoh
            sum_lai = sum_lai + col_cohort%lai(j)
         end do

         !----- Per-patch canopy geometry + constant forcing. -----------------------------!
         ageom%veg_height   = ctx%veg_height_bare
         do j = 1_ik, ncoh
            ageom%veg_height = max(ageom%veg_height, col_cohort%height(j))
         end do
         ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp

         call alloc_forcing(forc, ncoh)

         !----- Assemble the working bundle: adopt the owned per-patch reservoirs AND the        !
         !      PERSISTED per-cohort leaf_temp/leaf_water_mass carried on the cohort block (no     !
         !      reseeding). Capacity was ensured above; every field below is unconditionally        !
         !      (re)assigned from the site, so no alloc_patch_biophys seed call is needed here. -----!
         biophys%cas    = site%patch%cas(ip)
         biophys%soil_e = site%patch%soil_e(ip)
         biophys%soil_w = site%patch%soil_w(ip)
         biophys%snow   = site%patch%snow(ip)
         biophys%adapt_dt_last = site%patch%adapt_dt_last(ip)   ! issue #106: per-patch, not loop-carried
         !----- FROZEN slow soil-carbon pool (B2): a read-only snapshot for TODAY, held constant     !
         !      across the sub-step loop below (never written back -- the daily soil_carbon_step is   !
         !      the sole writer of the real site-level pool). site%patch%soil_carbon is ALWAYS         !
         !      allocated (default-initialised to 0 per patch at creation), so this copy is safe and   !
         !      bit-identical unconditionally: when soil_carbon_on=.false. nothing ever writes it, so   !
         !      it stays 0 -- the same value alloc_patch_biophys's intent(out) reset used to leave it   !
         !      at every patch (the OLD conditional skipped only a no-op copy of already-zero data).   !
         biophys%soil_carbon = site%patch%soil_carbon(ip)
         !----- FROZEN daily leaf/root-turnover shed-water rate (P4): same "read-only snapshot for  !
         !      TODAY, held constant across the sub-step loop" convention as soil_carbon just above. -!
         biophys%shed_water_rate = site%patch%shed_water_rate(ip)
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            biophys%leaf_temp(j) = site%cohort%leaf_temp(i)
            biophys%wood_temp(j) = site%cohort%wood_temp(i)
            !----- Tissue water is READ here, never written. The lazy PSI_INIT seed and the        !
            !      capacity clamp that used to live in this loop are a slow-loop concern -- capacity  !
            !      is a function of leaf/sapwood/root carbon, which only the slow loop changes -- and  !
            !      both branches WROTE, so an unbooked mass edit sat in the integrator's inner loop.   !
            !      They are now `reconcile_tissue_water_capacity`, run once per slow step before this  !
            !      loop (meds_stepper). Equivalent, not approximate: capacity is constant across a     !
            !      slow step, so the clamp is idempotent and the seed fires at most once. -------------!
            biophys%leaf_water_mass(j) = site%cohort%leaf_water_mass(i)
            biophys%wood_water_mass(j) = site%cohort%wood_water_mass(i)
            !----- Surface (interception film) water needs no lazy-init seed: 0 (bone dry) is a real  !
            !      initial condition here, not a placeholder -- a freshly-created cohort simply hasn't  !
            !      been rained on yet. -----------------------------------------------------------------!
            biophys%leaf_surf_water(j) = site%cohort%leaf_surf_water(i)
            biophys%wood_surf_water(j) = site%cohort%wood_surf_water(i)
         end do

         call ensure_aero_out_capacity(aero, ncoh)

         !----- n_fast_per_slow operator-split sweeps. Forcing is re-evaluated PER SUB-STEP (the   !
         !      diurnal cycle lives here): refresh the met overlay ctx_now, then fill_forcing +     !
         !      fill_aenv from it. CONSTANT path (do_forcing=.false.): ctx_now==ctx, so these        !
         !      reproduce the old build_forcing-once + fill_aenv sequence bit-identically.           !
         budget = column_budget_t()
         do isub = 1_ik, cfg%n_fast_per_slow
            !----- §8f: the met sample point within the sub-step (default 0.5 = midpoint) is applied   !
            !      when met_sample is built above. Only the cheap per-patch overlay write stays here --  !
            !      it targets ctx_now, which is per-patch (and will be per-THREAD) state. --------------!
            if (do_forcing) then
               met = met_sample(isub)
               call apply_met_to_ctx(ctx_now, met, f_ground)
            end if
            !----- Accumulate the sub-step air temperature for the daily-mean phenology driver. ------!
            red_site(RED_PHENO_TAIR, isub, ip) = ctx_now%air_temp
            call fill_forcing(forc, col_cohort, ctx_now, sum_lai)
            !----- RT join (§6.3): when forcing is on, REPLACE the LAI-share SW split with real     !
            !      per-cohort absorbed SW/PAR from the two-stream canopy radiation (ctx%rad_opt read !
            !      directly -- not the ctx_now overlay -- so the allocatable table is not deep-copied). !
            !----- LW emission base = the CAS temperature: the leaf energy balance linearizes leaf LW  !
            !      emission around tcas, so it needs abs_lw = NET LW AT tcas, and feeding the two-stream  !
            !      tcas as the canopy emission temperature makes abs_leaf(LW) exactly that. The           !
            !      prognostic CAS enthalpy is always valid here, unlike the lagged biophys%cas%can_temp. --!
            if (do_forcing) call apply_rt_forcing(forc, ncoh, col_cohort%pft, col_cohort%lai, col_cohort%wai,  &
                                 col_cohort%height, biophys%leaf_temp, biophys%wood_temp,                      &
                                 cas_temp_of_enthalpy(biophys%cas%can_enthalpy, biophys%cas%can_shv),          &
                                 biophys%soil_e%soil_temp(1), biophys%snow, ctx%col_config%snow,               &
                                 ctx%soil_albedo, ctx%soil_emiss, ctx%rad_opt, met, cfg%leaf_absorptance)
            call fill_aenv(aenv, biophys, ctx_now)
            !----- Slice to 1:ncoh (not the whole, possibly capacity-oversized backing array): the    !
            !      four accumulators are assumed-shape dummies in column_fast_step, so the ACTUAL      !
            !      argument's extent must equal col_cohort%n exactly, independent of the backing array's       !
            !      capacity (BB1 phase 1 pre-sizes it to the site-wide max, which can exceed ncoh). -----!
            !----- The per-cohort DIAGNOSTIC capture is passed only when the run reports per-cohort   !
            !      ecophysiology. Absent, the pre-pass never even asks the leaf kernel for the extra   !
            !      leaf_flux_t fields, so a production run pays nothing for a feature it is not using. !
            if (do_cdiag) then
               cdiag_buf(:, 1:ncoh) = 0.0_wp
               call column_fast_step(cfg%dt_fast, cfg, ctx_now%col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                     gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux, cdiag=cdiag_buf(:, 1:ncoh))
            else
               call column_fast_step(cfg%dt_fast, cfg, ctx_now%col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                     gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux)
            end if
            !----- Integrate the area-weighted CAS->atm latent flux -> site ET [kg/m2 = mm] over the step. !
            red_site(RED_ET, isub, ip) = site%patch%area(ip) * (le_flux / latent_heat_vap) * cfg%dt_fast
            !----- section 5.3 WORK: area-weight like every other site diagnostic, so a patch that     !
            !      needs more sub-steps is not double-counted by its area share. ------------------------!
            red_site(RED_INTEG_STEPS,  isub, ip) = site%patch%area(ip) * real(budget%integ_nsteps,   wp)
            red_site(RED_INTEG_REJ,    isub, ip) = site%patch%area(ip) * real(budget%integ_nrej,     wp)
            red_site(RED_SOIL_NSUB,    isub, ip) = site%patch%area(ip) * real(budget%soil_nsub,      wp)
            red_site(RED_HYDRO_NSUB,   isub, ip) = site%patch%area(ip) * real(budget%hydro_nsub,     wp)
            red_site(RED_HYDRO_THRASH, isub, ip) = site%patch%area(ip) * real(budget%hydro_thrash,   wp)
            red_site(RED_NONCONV,      isub, ip) = site%patch%area(ip) * real(budget%hydro_nonconv,  wp)
            !----- INTEGRATOR HEALTH, same area weighting. work_rk45_rescue is the one that decides       !
            !      whether an "RK45 run" was actually RK45: a nonzero total means some dt_fast steps      !
            !      were silently taken by the split path instead, which changes what the run measures.    !
            red_site(RED_RK45_RESCUE,  isub, ip) = site%patch%area(ip) * real(budget%rk45_rescue,    wp)
            red_site(RED_CLAMP_STAGE,  isub, ip) = site%patch%area(ip) * real(budget%clamp_stage_n,  wp)
            red_site(RED_CLAMP_COMMIT, isub, ip) = site%patch%area(ip) * real(budget%clamp_commit_n, wp)
            red_site(RED_CLAMP_MASS,   isub, ip) = site%patch%area(ip) * budget%clamp_mass
            red_site(RED_CLAMP_ENERGY, isub, ip) = site%patch%area(ip) * budget%clamp_energy
            !----- Integrate this sub-step's per-pool env scalar + matrix Rh into the day's totals    !
            !      (B2): dt_fast_days converts the instantaneous xi_step/rh_matrix_step (column_prepass !
            !      leaves both at 0 when soil_carbon_on=.false.) into the day-integral xi_int the daily  !
            !      soil_carbon_step consumes, and the accumulated Rh the audit cross-checks against.  ---!
            !      NOTE: dt_fast_days is declared at routine scope and listed `private`, not in a BLOCK.  !
            !      nvfortran 25.11 rejects a BLOCK construct anywhere inside a parallel region          !
            !      ("Unimplemented feature"), and MEDS must build on all three back ends. ---------------!
            if (cfg%soil_carbon_on) then
               dt_fast_days = cfg%dt_fast / day_sec
               site%patch%xi_accum(ip)%fast_grnd   = site%patch%xi_accum(ip)%fast_grnd               &
                                                     + budget%xi_step(IP_FAST_GRND)   * dt_fast_days
               site%patch%xi_accum(ip)%fast_soil   = site%patch%xi_accum(ip)%fast_soil               &
                                                     + budget%xi_step(IP_FAST_SOIL)   * dt_fast_days
               site%patch%xi_accum(ip)%struct_grnd = site%patch%xi_accum(ip)%struct_grnd             &
                                                     + budget%xi_step(IP_STRUCT_GRND) * dt_fast_days
               site%patch%xi_accum(ip)%struct_soil = site%patch%xi_accum(ip)%struct_soil             &
                                                     + budget%xi_step(IP_STRUCT_SOIL) * dt_fast_days
               site%patch%xi_accum(ip)%microbial   = site%patch%xi_accum(ip)%microbial               &
                                                     + budget%xi_step(IP_MICR)        * dt_fast_days
               site%patch%xi_accum(ip)%slow        = site%patch%xi_accum(ip)%slow                    &
                                                     + budget%xi_step(IP_SLOW)        * dt_fast_days
               site%patch%xi_accum(ip)%passive     = site%patch%xi_accum(ip)%passive                 &
                                                     + budget%xi_step(IP_PASSIVE)     * dt_fast_days
               site%patch%xi_accum(ip)%rh_fast_accum = site%patch%xi_accum(ip)%rh_fast_accum         &
                                                     + budget%rh_matrix_step * dt_fast_days
            end if
            !----- Sub-daily diagnostic PROBE (opt-in): per-(patch,sub-step) CAS temp / GPP / ET / soil. --!
            !      Thread-unsafe by construction (a `save`d unit on a shared file); load_meds_config     !
            !      rejects fast_probe together with n_threads > 1, so no guard is needed here. -----------!
            if (cfg%fast_probe .and. do_forcing)                                                    &
               call write_fast_probe(cfg, t_sample(isub), ip, ncoh, biophys, col_cohort, gpp_coh(1:ncoh),      &
                                     le_flux, ctx_now%rad_sw_top)
            !----- FAST (sub-daily) output staging: area-weight the LIVE per-sub-step site quantities   !
            !      onto the sub-step axis (patch areas sum to 1, so direct accumulation IS the site      !
            !      mean, mirroring et_accum). H and Rn are assembled here from the bulk conductances /    !
            !      absorbed radiation; per-cohort slabs are written by GLOBAL cohort slot (fixed within   !
            !      the <=1-day FAST file). main replays this into the FAST buffers (output_integrate_fast). !
            !----- Patch GPP and net absorbed radiation, once, for both the FAST staging and the patch  !
            !      diagnostics. h_flux is surfaced BY column_fast_step (like le_flux), computed from the   !
            !      in-call aero/aenv state -- computing it here from POST-call reads miscompiled to 0 on   !
            !      nvfortran (the same class as issue #7; le_flux/rnet read intent(in) forc, so safe). ----!
            gpp_patch = sum(gpp_coh(1:ncoh) * col_cohort%nplant(1:ncoh))    ! [umol/m2/s]
            rnet      = forc%abs_sw_ground + forc%abs_lw_ground                                          &
                        + sum(forc%abs_sw(1:ncoh)) + sum(forc%abs_lw(1:ncoh))
            if (do_fast) then
               w_area    = site%patch%area(ip)
               red_fast(isub,ip)%cas_temp      = w_area * biophys%cas%can_temp
               red_fast(isub,ip)%soil_temp_top = w_area * biophys%soil_e%soil_temp(1)
               red_fast(isub,ip)%gpp_rate      = w_area * gpp_patch
               red_fast(isub,ip)%le_flux       = w_area * le_flux
               red_fast(isub,ip)%h_flux        = w_area * h_flux
               red_fast(isub,ip)%rnet          = w_area * rnet
               red_fast(isub,ip)%sw_in         = w_area * ctx_now%rad_sw_top
               red_fast(isub,ip)%ustar         = w_area * aero%ustar
               red_fast(isub,ip)%air_temp      = w_area * ctx_now%air_temp
               !----- CARBON. budget%nee_last is the model's own NEE [umol/m2/s], sign-positive to     !
               !      the atmosphere -- the same number the CAS CO2 box is driven by, so the flux and  !
               !      the state it acts on cannot disagree. NPP is GPP net of the three MAINTENANCE    !
               !      respiration terms only (growth respiration is charged in the slow allocator, and !
               !      adding it here would double-count against npp_*_site). Reco then follows from    !
               !      the NEE identity, Reco = NEE + GPP, rather than being summed a second way.  -----!
               npp_patch = gpp_patch - sum((leaf_resp_coh(1:ncoh) + stem_resp_coh(1:ncoh)             &
                                            + root_resp_coh(1:ncoh)) * col_cohort%nplant(1:ncoh))
               red_fast(isub,ip)%nee_rate      = w_area * budget%nee_last
               red_fast(isub,ip)%npp_rate      = w_area * npp_patch
               red_fast(isub,ip)%reco_rate     = w_area * (budget%nee_last + gpp_patch)
               red_fast(isub,ip)%cas_co2       = w_area * biophys%cas%can_co2
               red_fast(isub,ip)%atm_co2       = w_area * ctx_now%co2_atm
               red_fast_soil_temp(1:nl,isub,ip)  = w_area * biophys%soil_e%soil_temp(1:nl)
               red_fast_soil_water(1:nl,isub,ip) = w_area * biophys%soil_w%theta(1:nl)
               !----- Per-cohort slabs are written by GLOBAL cohort slot, which is DISJOINT across      !
               !      patches (the CSR map partitions the flat SoA), so they go straight to mgr. -------!
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  mgr%fast_coh_ltemp(i,isub)  = biophys%leaf_temp(j)
                  mgr%fast_coh_gpp(i,isub)    = gpp_coh(j)
                  mgr%fast_coh_height(i,isub) = col_cohort%height(j)
               end do
            end if
            !----- FOLD the per-(cohort, sub-step) and per-patch DIAGNOSTICS into the site's        !
            !      dt-weighted accumulators (meds_site_diag_types). This is the whole point of the    !
            !      block: sub-daily resolution exists ONLY here, and before this everything but three !
            !      per-cohort quantities was recomputed ~48x/day and thrown away.                     !
            !                                                                                        !
            !      THREAD SAFETY: cohort slots [i0:i1] and patch slot ip are DISJOINT across patches  !
            !      (the CSR map partitions the flat SoA), so these writes need no reduction and stay  !
            !      byte-identical at any thread count -- the same argument the per-cohort state       !
            !      write-back below relies on.  ---------------------------------------------------!
            if (do_cdiag) then
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  site%cohort%diag%v(:, i) = site%cohort%diag%v(:, i) + cdiag_buf(:, j) * cfg%dt_fast
                  site%cohort%diag%w(i)    = site%cohort%diag%w(i)    + cfg%dt_fast
               end do
            end if
            if (do_pdiag) then
               call accumulate_patch_diag(site%patch%diag, ip, cfg%dt_fast, le_flux, h_flux, rnet,          &
                                          ctx_now%rad_sw_top, forc%abs_sw_ground, forc%abs_lw_ground,       &
                                          aero%ustar, aero%ggnet, aero%rough, aero%displace,               &
                                          biophys%cas%can_temp, biophys%cas%can_shv, biophys%cas%can_co2,   &
                                          gpp_patch, budget%nee_last, forc%rainfall + forc%snowfall,             &
                                          biophys%soil_e%soil_temp(1), budget%whole_energy%resid,           &
                                          budget%whole_water%resid)
            end if
            !----- Integrate GROSS GPP + maintenance-resp losses [umol/plant/s] -> [kgC/plant].  !
            !      Keep gross and loss terms SEPARATE (compute_carbon_allocation nets them; mirrors ED2). !
            do j = 1_ik, ncoh
               i = i0 + j - 1_ik
               site%cohort%gpp_accum(i)       = site%cohort%gpp_accum(i)       + gpp_coh(j)       * cfg%dt_fast * umol_2_kgC
               site%cohort%leaf_resp_accum(i) = site%cohort%leaf_resp_accum(i) + leaf_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%stem_resp_accum(i) = site%cohort%stem_resp_accum(i) + stem_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%root_resp_accum(i) = site%cohort%root_resp_accum(i) + root_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               !----- Running daily MAX of psi_leaf (#95). max(), not a sum: the daily maximum occurs !
               !      near dawn and IS the quantity that drives tomorrow's beta_stomata. --------------!
               site%cohort%dmax_psi_leaf_accum(i) = max(site%cohort%dmax_psi_leaf_accum(i), psi_leaf_coh(j))
            end do
         end do

         !----- Write the evolved state back to the site: per-patch reservoirs + per-cohort water. !
         site%patch%cas(ip)    = biophys%cas
         site%patch%soil_e(ip) = biophys%soil_e
         site%patch%soil_w(ip) = biophys%soil_w
         site%patch%snow(ip)   = biophys%snow
         site%patch%adapt_dt_last(ip) = biophys%adapt_dt_last
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            site%cohort%leaf_temp(i) = biophys%leaf_temp(j)
            site%cohort%wood_temp(i) = biophys%wood_temp(j)
            site%cohort%leaf_water_mass(i) = biophys%leaf_water_mass(j)
            site%cohort%wood_water_mass(i) = biophys%wood_water_mass(j)
            site%cohort%leaf_surf_water(i) = biophys%leaf_surf_water(j)
            site%cohort%wood_surf_water(i) = biophys%wood_surf_water(j)
         end do

         red_worst_energy(ip) = budget%whole_energy%worst
         red_worst_water(ip)  = budget%whole_water%worst
         red_budget_energy(ip) = budget%whole_energy
         red_budget_water(ip)  = budget%whole_water
         red_nfail(ip)        = budget%whole_energy%n_fail + budget%whole_water%n_fail
         end associate
      end do
      !$omp end parallel do

      !=========================================================================================!
      !  §7 C3 -- fold the staged per-(sub-step, patch) contributions back into the site, in the   !
      !  ORIGINAL serial order: patch outer, sub-step inner, exactly as the accumulate-in-place     !
      !  loop ran before threading. Floating-point addition is not associative, so this ordering IS  !
      !  the guarantee -- it makes every site diagnostic independent of the thread count, and equal   !
      !  bit-for-bit to what the serial driver produced. -----------------------------------------------!
      !=========================================================================================!
      site%pheno_tair_sum = 0.0_wp ; site%pheno_tair_n = 0_ik
      site%et_accum       = 0.0_wp
      site%work_integ_steps = 0.0_wp ; site%work_integ_rej  = 0.0_wp
      site%work_soil_nsub   = 0.0_wp ; site%work_hydro_nsub = 0.0_wp
      site%work_nonconv     = 0.0_wp ; site%work_hydro_thrash = 0.0_wp
      site%work_rk45_rescue = 0.0_wp ; site%work_clamp_stage  = 0.0_wp
      site%work_clamp_commit= 0.0_wp ; site%work_clamp_mass   = 0.0_wp
      site%work_clamp_energy= 0.0_wp
      do ip = 1_ik, npatch
         do isub = 1_ik, nsub
            site%et_accum          = site%et_accum          + red_site(RED_ET,            isub, ip)
            site%pheno_tair_sum    = site%pheno_tair_sum    + red_site(RED_PHENO_TAIR,    isub, ip)
            site%pheno_tair_n      = site%pheno_tair_n      + 1_ik
            site%work_integ_steps  = site%work_integ_steps  + red_site(RED_INTEG_STEPS,   isub, ip)
            site%work_integ_rej    = site%work_integ_rej    + red_site(RED_INTEG_REJ,     isub, ip)
            site%work_soil_nsub    = site%work_soil_nsub    + red_site(RED_SOIL_NSUB,     isub, ip)
            site%work_hydro_nsub   = site%work_hydro_nsub   + red_site(RED_HYDRO_NSUB,    isub, ip)
            site%work_hydro_thrash = site%work_hydro_thrash + red_site(RED_HYDRO_THRASH,  isub, ip)
            site%work_nonconv      = site%work_nonconv      + red_site(RED_NONCONV,       isub, ip)
            site%work_rk45_rescue  = site%work_rk45_rescue  + red_site(RED_RK45_RESCUE,   isub, ip)
            site%work_clamp_stage  = site%work_clamp_stage  + red_site(RED_CLAMP_STAGE,   isub, ip)
            site%work_clamp_commit = site%work_clamp_commit + red_site(RED_CLAMP_COMMIT,  isub, ip)
            site%work_clamp_mass   = site%work_clamp_mass   + red_site(RED_CLAMP_MASS,    isub, ip)
            site%work_clamp_energy = site%work_clamp_energy + red_site(RED_CLAMP_ENERGY,  isub, ip)
         end do
      end do
      !----- FAST-tier staging: each sub-step slot only ever took ONE contribution per patch, in       !
      !      increasing patch order, so folding over ip for fixed isub reproduces that order exactly. --!
      if (do_fast) then
         do isub = 1_ik, nsub
            do ip = 1_ik, npatch
               mgr%fast(isub)%cas_temp      = mgr%fast(isub)%cas_temp      + red_fast(isub,ip)%cas_temp
               mgr%fast(isub)%soil_temp_top = mgr%fast(isub)%soil_temp_top + red_fast(isub,ip)%soil_temp_top
               mgr%fast(isub)%gpp_rate      = mgr%fast(isub)%gpp_rate      + red_fast(isub,ip)%gpp_rate
               mgr%fast(isub)%le_flux       = mgr%fast(isub)%le_flux       + red_fast(isub,ip)%le_flux
               mgr%fast(isub)%h_flux        = mgr%fast(isub)%h_flux        + red_fast(isub,ip)%h_flux
               mgr%fast(isub)%rnet          = mgr%fast(isub)%rnet          + red_fast(isub,ip)%rnet
               mgr%fast(isub)%sw_in         = mgr%fast(isub)%sw_in         + red_fast(isub,ip)%sw_in
               mgr%fast(isub)%ustar         = mgr%fast(isub)%ustar         + red_fast(isub,ip)%ustar
               mgr%fast(isub)%air_temp      = mgr%fast(isub)%air_temp      + red_fast(isub,ip)%air_temp
               mgr%fast(isub)%nee_rate      = mgr%fast(isub)%nee_rate      + red_fast(isub,ip)%nee_rate
               mgr%fast(isub)%npp_rate      = mgr%fast(isub)%npp_rate      + red_fast(isub,ip)%npp_rate
               mgr%fast(isub)%reco_rate     = mgr%fast(isub)%reco_rate     + red_fast(isub,ip)%reco_rate
               mgr%fast(isub)%cas_co2       = mgr%fast(isub)%cas_co2       + red_fast(isub,ip)%cas_co2
               mgr%fast(isub)%atm_co2       = mgr%fast(isub)%atm_co2       + red_fast(isub,ip)%atm_co2
               mgr%fast_soil_temp(1:nl,isub)  = mgr%fast_soil_temp(1:nl,isub)                       &
                                              + red_fast_soil_temp(1:nl,isub,ip)
               mgr%fast_soil_water(1:nl,isub) = mgr%fast_soil_water(1:nl,isub)                      &
                                              + red_fast_soil_water(1:nl,isub,ip)
            end do
         end do
      end if

      if (present(worst_energy)) then
         worst_energy = 0.0_wp
         do ip = 1_ik, npatch ; worst_energy = max(worst_energy, red_worst_energy(ip)) ; end do
      end if
      if (present(worst_water)) then
         worst_water = 0.0_wp
         do ip = 1_ik, npatch ; worst_water = max(worst_water, red_worst_water(ip)) ; end do
      end if
      if (present(n_budget_fail)) n_budget_fail = sum(red_nfail(1:npatch))
      if (present(run_energy_budget) .or. present(run_water_budget)) then
         site_energy_budget = budget_t() ; site_water_budget = budget_t()
         do ip = 1_ik, npatch
            call budget_merge(site_energy_budget, red_budget_energy(ip), site%patch%area(ip))
            call budget_merge(site_water_budget,  red_budget_water(ip),  site%patch%area(ip))
         end do
         if (present(run_energy_budget)) call budget_merge(run_energy_budget, site_energy_budget, 1.0_wp)
         if (present(run_water_budget))  call budget_merge(run_water_budget,  site_water_budget,  1.0_wp)
      end if
      if (do_fast) mgr%fast_ready = .true.   ! signal main to replay + serialize the FAST tier
   end subroutine fast_dynamics

   !----- Sub-daily diagnostic PROBE. Opt-in ([fast].fast_probe): one CSV row per (patch, sub-step)  !
   !      with the fast-loop state the integrator/dt_fast evaluation (MEDS_INTEGRATOR_TEST.md §9)     !
   !      resolves on -- CAS temp, patch GPP, latent flux, top-soil temp, mean leaf temp. A saved     !
   !      unit opens the file on first use (header) and appends thereafter; the program-exit close     !
   !      flushes it. NOT part of the aggregation subsystem (that FAST tier stays deferred). ---------!
   subroutine write_fast_probe(cfg, t_sub, ip, ncoh, biophys, col_cohort, gpp_coh, le_flux, sw_in)
      type(meds_config_t),   intent(in) :: cfg
      type(meds_time_t),     intent(in) :: t_sub
      integer(ik),           intent(in) :: ip, ncoh
      type(patch_biophys_t), intent(in) :: biophys
      type(column_cohort_t), intent(in) :: col_cohort
      real(wp),              intent(in) :: gpp_coh(:), le_flux, sw_in
      integer, save :: unit           ! newunit returns a NEGATIVE handle -> track open state separately
      logical, save :: opened = .false.
      real(wp)      :: gpp_patch, leaf_temp_mean, wood_temp_mean

      if (.not. opened) then
         open(newunit=unit, file=trim(cfg%fast_probe_file), status='replace', action='write')
         write(unit,'(a)') 'datetime,patch,sw_in_W_m2,cas_temp_K,gpp_umol_m2_s,le_W_m2,soil_temp_top_K,leaf_temp_K,wood_temp_K'
         opened = .true.
      end if

      gpp_patch = 0.0_wp ; leaf_temp_mean = 0.0_wp ; wood_temp_mean = 0.0_wp
      if (ncoh > 0_ik) then
         gpp_patch      = sum(gpp_coh(1:ncoh) * col_cohort%nplant(1:ncoh))    ! per-plant [umol/plant/s] x nplant -> [umol/m2/s]
         leaf_temp_mean = sum(biophys%leaf_temp(1:ncoh)) / real(ncoh, wp)
         wood_temp_mean = sum(biophys%wood_temp(1:ncoh)) / real(ncoh, wp)
      end if

      write(unit,'(a,",",i0,7(",",es13.6))') trim(time_to_string(t_sub)), ip,   &
            sw_in, biophys%cas%can_temp, gpp_patch, le_flux, biophys%soil_e%soil_temp(1), leaf_temp_mean, wood_temp_mean
      flush(unit)
   end subroutine write_fast_probe

   !----- Grow-only capacity check for the per-patch forcing buffers (MEDS_NUMERICS_SCOPING.md BB1  !
   !      phase 1): `forc` is reused across the patch loop (and, since the caller now pre-sizes it   !
   !      to the site-wide max cohort count before the loop, across the WHOLE loop with zero          !
   !      reallocation). fill_forcing/apply_rt_forcing write indices 1..col_cohort%n (never size(forc%...)), !
   !      so reusing a larger patch's leftover capacity for a smaller one is bit-identical. -----------!
   subroutine alloc_forcing(forc, ncoh)
      type(column_forcing_t), intent(inout) :: forc
      integer(ik),            intent(in)    :: ncoh
      if (allocated(forc%abs_sw)) then
         if (size(forc%abs_sw) < ncoh) deallocate(forc%abs_sw, forc%abs_lw, forc%abs_par,        &
                                                  forc%abs_sw_wood, forc%abs_lw_wood)
      end if
      if (.not. allocated(forc%abs_sw)) allocate(forc%abs_sw(ncoh), forc%abs_lw(ncoh),           &
                                                 forc%abs_par(ncoh), forc%abs_sw_wood(ncoh), forc%abs_lw_wood(ncoh))
   end subroutine alloc_forcing

   !----- Fill the per-patch prescribed forcing from the (possibly per-sub-step) reference met. !
   subroutine fill_forcing(forc, col_cohort, ctx, sum_lai)
      type(column_forcing_t), intent(inout) :: forc
      type(column_cohort_t),  intent(in)    :: col_cohort
      type(fast_context_t),   intent(in)    :: ctx
      real(wp),               intent(in)    :: sum_lai
      integer(ik) :: j
      forc%enthalpy_atm  = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
      forc%shv_atm       = ctx%shv_atm
      forc%co2_atm       = ctx%co2_atm
      forc%abs_sw_ground = ctx%rad_sw_ground
      forc%abs_lw_ground = 0.0_wp
      forc%rainfall        = ctx%rainfall
      forc%snowfall         = ctx%snowfall                 ! frozen rainfall -> snow accumulation
      forc%air_temp          = ctx%air_temp              ! rainfall enthalpy reference (snow/rain-on-snow)
      forc%par_per_w     = 2.1_wp                    ! LAI-split path: total-SW->PAR blend (abs_par == abs_sw)
      !----- Split the canopy-top shortwave across cohorts by LAI share (MVP; the RT join (§6.3) !
      !      replaces this with real per-cohort absorbed SW/PAR when forcing is on).             !
      do j = 1_ik, col_cohort%n
         if (sum_lai > tiny_num) then
            forc%abs_sw(j) = ctx%rad_sw_top * col_cohort%lai(j) / sum_lai
         else
            forc%abs_sw(j) = 0.0_wp
         end if
         forc%abs_par(j) = forc%abs_sw(j)            ! no PAR/NIR split in the LAI path -> PAR==SW (biased high)
         forc%abs_lw(j)  = 0.0_wp
         forc%abs_sw_wood(j) = 0.0_wp                ! MVP const/no-forcing path: no wood absorption
         forc%abs_lw_wood(j) = 0.0_wp
      end do
   end subroutine fill_forcing

   !----- The met-source shim (design §6.2/§6.5, retired at P1): copy the instantaneous          !
   !      met_forcing_t's raw scalars into fast_context_t's met fields, so fill_forcing/fill_aenv  !
   !      keep their present logic unchanged. cosz/leaf_temp/ustar/can_co2 are NOT overwritten     !
   !      (they are prognostic or derived). Ground SW scales with canopy-top SW at the reference    !
   !      transmittance f_ground (reproduces the ad-hoc rad_sw_ground=60 at swdown=400).            !
   subroutine apply_met_to_ctx(ctx, met, f_ground)
      type(fast_context_t), intent(inout) :: ctx
      type(met_forcing_t),  intent(in)    :: met
      real(wp),             intent(in)    :: f_ground
      ctx%air_temp      = met%tair_k
      ctx%shv_atm       = met%qair
      ctx%press         = met%psurf_pa
      ctx%rho_air       = met%rho_air
      ctx%co2_atm       = met%co2
      ctx%u_ref         = met%wind
      ctx%rainfall        = met%rainf
      ctx%snowfall         = met%snowfall
      ctx%rad_sw_top    = met%swdown()
      ctx%rad_sw_ground = f_ground * met%swdown()
   end subroutine apply_met_to_ctx

   !----- RT join (§6.3): run the two-stream canopy radiation for this patch and OVERWRITE the      !
   !      LAI-share SW split in `forc` with real per-cohort absorbed SW (leaf energy) + PAR           !
   !      (photosynthesis) + NET longwave (leaf + ground) + below-canopy ground SW.                   !
   !      Cohorts are gathered height-DESCENDING (top=1) but the two-stream wants BOTTOM(1)->TOP, so  !
   !      we build an ascending-height permutation `perm` and inverse-scatter the outputs.            !
   !      NOTE on PAR: the leaf model's `env%par` is INCIDENT PAR (it re-applies cfg%leaf_absorptance !
   !      internally for electron transport), so we divide the two-stream ABSORBED VIS by that same   !
   !      absorptance to hand back an incident-equivalent PAR -- otherwise leaf absorptance would be  !
   !      applied twice. abs_sw stays true ABSORBED SW (the leaf energy balance wants absorbed).      !
   subroutine apply_rt_forcing(forc, ncoh, pft, lai, wai, height, leaf_temp, wood_temp, tcas,          &
                               soil_temp_top, snow, snow_params, soil_albedo, soil_emiss, rad_opt, met, &
                               leaf_absorptance)
      type(column_forcing_t),  intent(inout) :: forc
      integer(ik),             intent(in)    :: ncoh
      integer(ik),             intent(in)    :: pft(:)           !< per-cohort PFT (gather order, top first)
      real(wp),                intent(in)    :: lai(:), wai(:)   !< [m2/m2] leaf / wood area index
      real(wp),                intent(in)    :: height(:)        !< [m] cohort height
      real(wp),                intent(in)    :: leaf_temp(:), wood_temp(:)   !< [K] lagged tissue temperatures
      real(wp),                intent(in)    :: tcas             !< [K] CAS temperature (LW emission base)
      real(wp),                intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground emission)
      type(snow_column_t),     intent(in)    :: snow             !< pack store (albedo / emissivity ramp)
      type(snow_params_t),     intent(in)    :: snow_params
      real(wp),                intent(in)    :: soil_albedo(:)   !< [-] bare-soil albedo per band
      real(wp),                intent(in)    :: soil_emiss       !< [-] bare-soil emissivity
      type(rad_pft_optics_t),  intent(in)    :: rad_opt          !< per-PFT canopy optics (two-stream)
      type(met_forcing_t),     intent(in)    :: met
      real(wp),                intent(in)    :: leaf_absorptance !< [-] leaf PAR absorptance (incident-PAR conversion)
      integer(ik) :: j, k, ig, imin
      integer(ik) :: perm(ncoh), pft_bt(ncoh)
      real(wp)    :: lai_bt(ncoh), wai_bt(ncoh), tcan_bt(ncoh)
      logical     :: used(ncoh)
      real(wp)    :: hmin, lf_bt
      type(rad_forcing_t)   :: rf
      type(rad_flux_t)      :: flux
      type(ground_optics_state_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)
      real(wp) :: snow_fl, snow_fc

      !----- A bare patch (ncoh == 0) is NOT special-cased: the zero-trip perm/scatter loops fall     !
      !      through and canopy_radiation's own empty-canopy branch returns the correct NET ground SW  !
      !      (incident * (1 - soil albedo)), so a patch shedding its last cohort stays continuous.     !

      !----- perm: gather-indices in ASCENDING height (bottom -> top). Selection sort (ncoh small). !
      used = .false.
      do j = 1_ik, ncoh
         imin = 0_ik ; hmin = huge(1.0_wp)
         do k = 1_ik, ncoh
            if (.not. used(k) .and. height(k) <= hmin) then ; hmin = height(k) ; imin = k ; end if
         end do
         perm(j) = imin ; used(imin) = .true.
         pft_bt(j) = pft(imin) ; lai_bt(j) = lai(imin)
         wai_bt(j) = wai(imin)
         !----- LW emission temperature (P1): the cohort's AREA-WEIGHTED effective radiative temperature  !
         !      so it emits at leaf_temp over its LAI and wood_temp over its WAI (T^4 weights telescope    !
         !      with leaf_frac) -- so the RT FIELD (inter-cohort/sky/ground LW) reflects both tissue temps !
         !      instead of the single air temp. Lagged (start-of-sub-step). NOTE: the leaf/wood energy     !
         !      balances keep their LOCAL emission base at tcas (split)/leaf_temp (picard); re-basing the  !
         !      single-pass split on the lagged element temp is a positive-feedback instability, so the    !
         !      per-element "counted once" base is a documented residual (design §8/P1).                    !
         lf_bt      = lai(imin) / max(lai(imin) + wai(imin), tiny_num)
         tcan_bt(j) = (lf_bt * leaf_temp(imin) ** 4                                             &
                       + (1.0_wp - lf_bt) * wood_temp(imin) ** 4) ** 0.25_wp
      end do

      !----- rad_forcing_t from met (§6.3 mapping table; all W/m2, direct assignment). -----------!
      call alloc_rad_forcing(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = met%cosz
      rf%incid_beam(RAD_VIS) = met%par_beam ; rf%incid_diff(RAD_VIS) = met%par_diffuse
      rf%incid_beam(RAD_NIR) = met%nir_beam ; rf%incid_diff(RAD_NIR) = met%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp       ; rf%incid_diff(RAD_LW)  = met%lwdown

      !----- ground optics (MVP soil albedo/emiss from ctx; ground skin temp from the soil column). !
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = soil_albedo ; surf%soil_emiss = soil_emiss
      surf%soil_temp   = soil_temp_top
      !----- Snow raises the ground albedo/emissivity + emits off the snow surface (design §4f), RAMPED   !
      !      by the Niu-Yang07 snow-cover fraction so a partial pack gives a partial (continuous) albedo   !
      !      -- no threshold cliff. VIS/NIR fresh<->aged interpolated by the lagged surface liquid fraction. !
      if (snow%nlayer >= 1_ik .and. snow%swe(1) > snow_params%tiny_snow_mass) then
         associate (sp => snow_params)
            snow_fc = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), sp)
            snow_fl = snow%snow_fliq(1)
            surf%soil_albedo(RAD_VIS) = (1.0_wp - snow_fc) * soil_albedo(RAD_VIS)                 &
                 + snow_fc * ((1.0_wp - snow_fl) * sp%albedo_vis_fresh + snow_fl * sp%albedo_vis_aged)
            surf%soil_albedo(RAD_NIR) = (1.0_wp - snow_fc) * soil_albedo(RAD_NIR)                 &
                 + snow_fc * ((1.0_wp - snow_fl) * sp%albedo_nir_fresh + snow_fl * sp%albedo_nir_aged)
            surf%soil_emiss = (1.0_wp - snow_fc) * soil_emiss + snow_fc * sp%snow_emiss
            surf%soil_temp  = (1.0_wp - snow_fc) * soil_temp_top + snow_fc * snow%snow_temp(1)
         end associate
      end if
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, rf%grnd_refl, rf%grnd_emiss)

      call canopy_radiation(rad_opt, rf, ncoh, pft_bt, lai_bt, wai_bt, tcan_bt, flux)

      !----- inverse-scatter: RT index j (bottom->top) maps to gather index perm(j). --------------!
      do j = 1_ik, ncoh
         ig = perm(j)
         forc%abs_sw(ig)  = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! total ABSORBED leaf SW (energy)
         forc%abs_par(ig) = flux%abs_leaf(RAD_VIS, j) / max(leaf_absorptance, tiny_num)  ! -> INCIDENT-equiv PAR
         forc%abs_lw(ig)  = flux%abs_leaf(RAD_LW, j)                                ! NET leaf LW at tcas (emission incl.)
         forc%abs_sw_wood(ig) = flux%abs_wood(RAD_VIS, j) + flux%abs_wood(RAD_NIR, j)  ! ABSORBED wood SW (WAI share)
         forc%abs_lw_wood(ig) = flux%abs_wood(RAD_LW, j)                               ! NET wood LW
      end do
      forc%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                         + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forc%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)          ! NET ground LW (soil emission incl.)
      forc%par_per_w     = PAR_W_2_UMOL                        ! true VIS absorbed -> photon flux
   end subroutine apply_rt_forcing

   !----- Fill the aerodynamics env from the reference met + the patch's current CAS/ground. -!
   subroutine fill_aenv(aenv, biophys, ctx)
      type(aero_env_t),     intent(inout) :: aenv
      type(patch_biophys_t), intent(in)   :: biophys
      type(fast_context_t), intent(in)    :: ctx
      aenv%u_ref = ctx%u_ref ; aenv%zref = ctx%zref ; aenv%press = ctx%press ; aenv%rho_air = ctx%rho_air
      !----- The potential-temperature conversion and the CAS/ground refresh now live in            !
      !      meds_canopy_types/meds_soil_types (issue #97), so tests and probes assemble `aenv` through the SAME !
      !      routine this driver does instead of by a parallel hand-written copy -- which is how     !
      !      every column test ended up leaving `theta_atm` at its 298.15 K default. `zref` must be  !
      !      assigned before set_aero_env_atm, which reads it. -------------------------------------!
      call set_aero_env_atm(aenv, ctx%air_temp, ctx%shv_atm, ctx%co2_atm)
      call set_aero_env_canopy(aenv, biophys%cas%can_temp, biophys%cas%can_shv, biophys%cas%can_co2,           &
                               biophys%soil_e%soil_temp(1))
   end subroutine fill_aenv

   !=======================================================================================!
   !  Fold ONE (patch, sub-step) sample into the per-patch diagnostic accumulator.            !
   !                                                                                          !
   !  Every field is dt-weighted here and normalized on read, so the result is a correct        !
   !  time-mean even though the adaptive controller may have taken a different number of inner   !
   !  steps in different sub-steps. Values are the SAME numbers the physics just used -- nothing  !
   !  is recomputed, which is the point: before this they were computed and dropped.              !
   !=======================================================================================!
   subroutine accumulate_patch_diag(pd, ip, dt, le_flux, h_flux, rnet, sw_in, sw_ground, lw_ground,       &
                                    ustar, ggnet, rough, displace, cas_temp, cas_shv, cas_co2, gpp, nee,   &
                                    precip_total, ground_temp, resid_energy, resid_water)
      type(patch_diag_block), intent(inout) :: pd
      integer(ik),            intent(in)    :: ip
      real(wp),               intent(in)    :: dt                        !< [s]        sample weight
      real(wp),               intent(in)    :: le_flux, h_flux           !< [W/m2]     CAS -> atmosphere
      real(wp),               intent(in)    :: rnet                      !< [W/m2]     net all-wave radiation absorbed
      real(wp),               intent(in)    :: sw_in                     !< [W/m2]     incident shortwave, canopy top
      real(wp),               intent(in)    :: sw_ground, lw_ground      !< [W/m2]     ground SW / net LW
      real(wp),               intent(in)    :: ustar, ggnet              !< [m/s]      friction velocity, ground conductance
      real(wp),               intent(in)    :: rough, displace           !< [m]        roughness, displacement height
      real(wp),               intent(in)    :: cas_temp, cas_shv, cas_co2 !< [K],[kg/kg],[umol/mol] canopy air
      real(wp),               intent(in)    :: gpp, nee                  !< [umol/m2/s] gross uptake, net exchange (+ to atm)
      real(wp),               intent(in)    :: precip_total              !< [kg/m2/s]  rain + snow
      real(wp),               intent(in)    :: ground_temp               !< [K]        top soil-node temperature
      real(wp),               intent(in)    :: resid_energy, resid_water !< [J/m2],[kg/m2] this step's SIGNED ledger residuals
      pd%v(PD_LE,           ip) = pd%v(PD_LE,           ip) + le_flux                * dt
      pd%v(PD_H,            ip) = pd%v(PD_H,            ip) + h_flux                 * dt
      pd%v(PD_RNET,         ip) = pd%v(PD_RNET,         ip) + rnet                   * dt
      pd%v(PD_SW_IN,        ip) = pd%v(PD_SW_IN,        ip) + sw_in                  * dt
      pd%v(PD_SW_GROUND,    ip) = pd%v(PD_SW_GROUND,    ip) + sw_ground              * dt
      pd%v(PD_LW_GROUND,    ip) = pd%v(PD_LW_GROUND,    ip) + lw_ground              * dt
      pd%v(PD_USTAR,        ip) = pd%v(PD_USTAR,        ip) + ustar                  * dt
      pd%v(PD_GGNET,        ip) = pd%v(PD_GGNET,        ip) + ggnet                  * dt
      pd%v(PD_ROUGH,        ip) = pd%v(PD_ROUGH,        ip) + rough                  * dt
      pd%v(PD_DISPLACE,     ip) = pd%v(PD_DISPLACE,     ip) + displace               * dt
      pd%v(PD_CAS_TEMP,     ip) = pd%v(PD_CAS_TEMP,     ip) + cas_temp               * dt
      pd%v(PD_CAS_SHV,      ip) = pd%v(PD_CAS_SHV,      ip) + cas_shv                * dt
      pd%v(PD_CAS_CO2,      ip) = pd%v(PD_CAS_CO2,      ip) + cas_co2                * dt
      pd%v(PD_GPP,          ip) = pd%v(PD_GPP,          ip) + gpp                    * dt
      pd%v(PD_NEE,          ip) = pd%v(PD_NEE,          ip) + nee                    * dt
      !----- Transpiration as a WATER flux [kg/m2/s]: the latent flux is the canopy-air -> atmosphere  !
      !      total, so this is the evaporative flux the CAS actually shed, not a stomatal-only term.    !
      !      The stomatal share is available per cohort (CD_TRANSP) for anyone who needs the split.     !
      pd%v(PD_TRANSP,       ip) = pd%v(PD_TRANSP,       ip) + (le_flux/latent_heat_vap) * dt
      pd%v(PD_PRECIP,       ip) = pd%v(PD_PRECIP,       ip) + precip_total           * dt
      pd%v(PD_GROUND_TEMP,  ip) = pd%v(PD_GROUND_TEMP,  ip) + ground_temp            * dt
      !----- Whole-column budget residuals. These are the numbers that decide whether anything above  !
      !      this line can be believed, which is why they are captured on the same tick rather than    !
      !      left to an assertion nobody reads. budget%*%resid is THIS step's SIGNED imbalance [J/m2,    !
      !      kg/m2]; summed here and divided by the aggregation's sum(dt) it is the mean leak RATE      !
      !      [W/m2, kg/m2/s], sign-positive when store appears from nowhere. (It used to accumulate    !
      !      worst*dt -- a running max in J/m2 that the registry then labelled W/m2.) -----------------!
      pd%v(PD_RESID_ENERGY, ip) = pd%v(PD_RESID_ENERGY, ip) + resid_energy
      pd%v(PD_RESID_WATER,  ip) = pd%v(PD_RESID_WATER,  ip) + resid_water
      pd%w(ip)                  = pd%w(ip)                  + dt
   end subroutine accumulate_patch_diag

end module meds_fast_dynamics
