! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_dynamics -- the per-SITE fast-biophysics TIER ORCHESTRATOR (the fast-cadence   !
! peer of meds_vegetation_dynamics; advance_one_step fans out to both). It owns the        !
! ORCHESTRATION only (docs/science/order_of_processes.md §3): for each patch it gathers    !
! the demographic cohort slice into the column_cohort_t buffer, assembles a                !
! patch_biophys_t working bundle from the state-hub-owned per-patch reservoirs             !
! (cas/soil_energy/soil_water), runs n_fast_per_slow sub-steps of the per-patch kernel              !
! column_fast_step, and writes the evolved reservoirs back to the site. The static         !
! column_config_t + base met arrive via fast_context_t (the caller builds them -- no model !
! parameters are hard-coded here); per-cohort                                              !
! leaf_temp/wood_temp/leaf_water_mass/wood_water_mass/leaf_surf_water/wood_surf_water are  !
! PERSISTED on the cohort block and adopted here each slow step (no reseeding) -- same as  !
! the soil + CAS reservoirs, this is genuine cross-slow-step memory.                       !
!                                                                                          !
! The stepper calls fast_dynamics before the slow loop when cfg%fast_biophysics_on. This   !
! is the fast->slow seam's fast half: the slow loop's carbon allocation reads the GPP and  !
! maintenance respiration accumulated here (compute_carbon_allocation).                    !
!==========================================================================================!
module meds_fast_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, umol_2_kgC, grav, cp_air, latent_heat_vap, day_sec, p_std, &
                                     par_w_2_umol,                                                            &
                                     stefan
   use meds_numerics,         only : ascending_order
   use meds_config,           only : meds_config_t
   use meds_soil_column_config, only : soil_column_params, soil_column_start_theta, top_soil_psi
   use meds_budget_check,     only : budget_t, budget_merge
   use meds_biogeochem_types, only : IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   use meds_therm_lib,           only : cas_enthalpy_of_temp, temp_to_internal_energy, specific_humidity_to_vpd
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts, build_rad_optics
   use meds_column_view, only : copy_column_cohort
   use meds_time,             only : meds_time_t, time_advance_seconds, time_to_string
   use meds_output_types,     only : output_buffers_t
   use meds_diag_vars,        only : N_CDIAG, N_PDIAG, patch_diag_block, polygon_diag_block, PD_LE, PD_H, PD_RNET, PD_SW_GROUND, &
                                     PD_LW_GROUND, PD_SW_IN_VIS, PD_SW_IN_NIR, PD_SW_UP_VIS, PD_SW_UP_NIR, PD_LW_UP, PD_USTAR, &
                                     PD_GGNET, PD_ROUGH, PD_DISPLACE, PD_CAS_TEMP, PD_CAS_SHV, PD_CAS_CO2, PD_GPP, PD_NEE, &
                                     PD_TRANSP, PD_NPP, PD_RECO, PD_SW_UP, PD_SKIN_TEMP, PD_SOIL_TEMP_TOP, PD_CAS_VPD, &
                                     PD_W_SURFACE, PD_RESID_ENERGY, PD_RESID_WATER, PD_WIND_CAS_TOP, PD_TAIR_CAS_TOP, &
                                     PD_Z_CAS_TOP, PY_SW_IN, PY_PRECIP, PY_TAIR, PY_QAIR, PY_PSURF, PY_WIND, PY_LWDOWN, &
                                     PY_PAR_BEAM, PY_PAR_DIFFUSE, PY_NIR_BEAM, PY_NIR_DIFFUSE, PY_SNOWFALL, PY_CO2, PY_COSZ, &
                                     PY_RHO_AIR, PY_RAINF, N_PYDIAG, cohort_diag_grow, patch_diag_grow
   use meds_constants,     only : n_soil_layer_max
   use meds_patch_state_vars,   only : xi_accum_t, snow_column_t
   use meds_forcing_types,    only : met_source_t, met_cursor_t, met_forcing_t
   use meds_forcing_kernels,  only : precip_phase
   use meds_met_driver,       only : met_advance, met_instant
   use meds_lapse_rate,       only : met_to_cas_top
   use meds_canopy_aerodynamics, only : canopy_roughness, open_canopy_fraction
   use meds_cohort_state_vars, only : set_cohort_tissue_water
   use meds_state_defaults,    only : DMAX_PSI_LEAF_UNSET, DMAX_PSI_LEAF_ACCUM_RESET
   use meds_site_state_vars,  only : site_t
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, ensure_aero_out_capacity, rad_pft_optics_t, &
                                 rad_forcing_t, rad_flux_t, alloc_rad_forcing, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW, &
                                 set_aero_env_atm, set_aero_env_canopy
   use meds_fast_types, only : patch_biophys_t, ensure_patch_biophys_capacity, fast_workspace_t
   use meds_water_retention, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : snow_params_t
   use meds_canopy_types, only : ground_optics_state_t
   use meds_canopy_radiation, only : canopy_radiation, ground_optics
   use meds_ground_biophysics, only : snow_cover_fraction
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t,        &
                                     GRP_THETA,                                                   &
                                     column_budget_t,                                             &
                                     ensure_column_cohort_capacity, apply_hydraulics_config,      &
                                     apply_canopy_film_config
   use meds_fast_step,       only : column_fast_step
   use meds_water_retention,  only : water_content, clamp_water_to_capacity
   !$ use omp_lib,            only : omp_get_thread_num
   implicit none
   private

   public :: fast_context_t, init_fast_reservoirs, fast_dynamics, build_fast_context

   !----- PAR (VIS) energy -> photon flux: par_w_2_umol (meds_constants). Photosynthesis uses it  !
   !      ONLY on the RT path (true absorbed PAR; the no-forcing LAI-split path keeps the 2.1     !
   !      total-SW blend); the phenology light count uses it on the PAR reaching each cohort's    !
   !      top.                                                                                    !

   !----- DETERMINISTIC REDUCTIONS (docs/science/numerical_scheme.md §6a). The site-level fast-loop !
   !      accumulators are written once per (patch, sub-step) and folded into site%... afterwards.  !
   !      Staging them in ONE array, indexed by the slots below, lets the patch loop write DISJOINT !
   !      cells (no reduction clause, no atomics) and lets the fold run in the (patch, sub-step)    !
   !      order a serial loop uses. That order is the whole point: OpenMP `reduction(+:)` sums      !
   !      thread partials in thread-ARRIVAL order, so the last bits of every site diagnostic would  !
   !      move with the thread count -- and byte-identical netCDF across thread counts is a hard    !
   !      requirement, not a nicety. Folding in patch order reproduces the serial fold EXACTLY. ----!
   integer(ik), parameter :: RED_ET            =  1_ik   !< site ET [kg/m2], area-weighted
   integer(ik), parameter :: RED_PHENO_TAIR    =  2_ik   !< daily-mean air-temperature accumulator [K]
   integer(ik), parameter :: RED_INTEG_STEPS   =  3_ik   !< integrator work: accepted steps
   integer(ik), parameter :: RED_INTEG_REJ     =  4_ik   !< integrator work: rejected steps
   integer(ik), parameter :: RED_WATER_NSUB    =  5_ik   !< water network: accepted sub-steps
   integer(ik), parameter :: RED_WATER_NEWTON  =  6_ik   !< water network: Newton iterations
   integer(ik), parameter :: RED_WATER_REJECT  =  7_ik   !< water network: rejected sub-steps
   integer(ik), parameter :: RED_WATER_NFAIL   =  8_ik   !< water network: sub-steps taken unconverged
   integer(ik), parameter :: RED_RK45_RESCUE   =  9_ik   !< dt_fast steps RK45 handed to ARK (rescues)
   integer(ik), parameter :: RED_CLAMP_STAGE   = 10_ik   !< stage-level state clamps
   !----- The rain and snow each patch received, split at its canopy-air top (fill_forcing), area-  !
   !      weighted: the forcing outputs report the split the model used, not one at the forcing's     !
   !      own height. --------------------------------------------------------------------------------!
   integer(ik), parameter :: RED_RAINF         = 11_ik   !< rain [kg/m2/s]
   integer(ik), parameter :: RED_SNOWFALL      = 12_ik   !< snow [kg/m2/s]
   integer(ik), parameter :: RED_RUNOFF_SURF   = 13_ik   !< surface runoff [kg/m2], area-weighted
   integer(ik), parameter :: RED_DRAINAGE      = 14_ik   !< drainage at the column's bottom [kg/m2]
   integer(ik), parameter :: N_RED             = 14_ik

   !----- Everything the fast driver needs beyond the site + cfg: the static column config plus  !
   !      the reference met + initial soil state. The CALLER builds this (build_fast_context,    !
   !      from the run config); the scalar reference met is the constant, horizontally uniform   !
   !      climate a run without a forcing source uses.                                           !
   type :: fast_context_t
      type(column_config_t) :: col_config                 !< soil/thermal/hydro/aero/resp column config
      real(wp) :: u_ref   = 2.0_wp                  !< [m/s] reference wind (at the canopy-air top)
      real(wp) :: press   = 101325.0_wp             !< [Pa]
      real(wp) :: rho_air = 1.2_wp                  !< [kg/m3]
      real(wp) :: air_temp = 288.0_wp               !< [K]        reference-level air temperature
      real(wp) :: shv_atm  = 0.008_wp               !< [kg/kg]    reference-level specific humidity
      real(wp) :: co2_atm  = 400.0_wp               !< [umol/mol] free-atmosphere CO2
      real(wp) :: rad_sw_top    = 400.0_wp          !< [W/m2] shortwave into the canopy (leaves)
      real(wp) :: rad_sw_ground = 60.0_wp           !< [W/m2] shortwave reaching the ground
      real(wp) :: precip          = 0.0_wp            !< [kg/m2/s] total precipitation (split like any record)
      real(wp) :: theta_init(n_soil_layer_max) = 0.30_wp   !< [m3/m3] initial soil moisture, per layer
      real(wp) :: soil_temp_init  = 288.0_wp        !< [K]     initial soil + CAS temperature
      !----- Threads for this polygon's patch loop: [run].n_threads in a site run, 1 in a region,   !
      !      whose threads go to the polygons instead (polygon_prepare sets it). ------------------!
      integer(ik) :: patch_threads = 1_ik
      real(wp) :: veg_height_bare = 1.0_wp          !< [m] canopy height for a cohort-free patch
      !----- Canopy-RT optics (per-PFT spectral/angle table + ground surface), for the RT-driven   !
      !      per-cohort absorbed SW/PAR on the forcing path. Built once from the [pft] table.        !
      type(rad_pft_optics_t) :: rad_opt             !< per-PFT, per-band leaf/wood scattering + LIDF
      real(wp) :: soil_albedo(3) = [0.15_wp, 0.30_wp, 0.0_wp]  !< [VIS,NIR,LW] ground albedo (LW=0)
      real(wp) :: soil_emiss     = 0.95_wp          !< [-] ground longwave emissivity
   end type fast_context_t

contains

   !=======================================================================================!
   !  Build the fast-loop context for the production run. The scalar reference met +       !
   !  reservoir seeds keep the fast_context_t defaults, which are constant and             !
   !  horizontally uniform -- the reference climate used when no forcing file is           !
   !  configured. The static column config `col_config` is assembled here FROM THE RUN     !
   !  CONFIG: the soil column comes from [soil], the respiration factors and canopy        !
   !  optics are per-PFT traits, and the plant-hydraulics constants come from              !
   !  [hydraulics]. Nothing in this routine is a hard-coded model parameter.               !
   !=======================================================================================!
   subroutine build_fast_context(cfg, ctx)
      type(meds_config_t),  intent(in)  :: cfg
      type(fast_context_t), intent(out) :: ctx
      !----- The physical soil column, from [soil]: its layers, and each layer's texture,         !
      !      one for the whole column or one per layer from the soil maps (soil_column_params). -!
      call soil_column_params(cfg%soil_column, ctx%col_config%soil, ctx%col_config%soil_thermal)
      !----- Plant hydraulics: ONE seam builds the PER-PFT table. The shared [hydraulics] block is    !
      !      the base (cfg defaults unless a [hydraulics] block overrides); whatever per-PFT traits   !
      !      the [pft] table supplied are laid over it, and each entry gets its vulnerability lookup  !
      !      table (consulted only when wood_kexp is not 1 or 2). ------------------------------------!
      call apply_hydraulics_config(cfg%hydraulics, cfg%pft, ctx%col_config%hydraulics_table)
      call apply_canopy_film_config(cfg%pft, ctx%col_config%leaf_surf_water_max, ctx%col_config%wood_surf_water_max)
      call build_leaf_photo_table(cfg, ctx%col_config%leaf_photo)    ! per-PFT leaf parameters, once per run
      !----- The canopy interception film (+ film evaporation, dew) on/off, [fast] block. --------!
      ctx%col_config%canopy_water_on    = cfg%canopy_water_on
      !----- Fast-loop biophysics run-config from the [soil]/[snow]/[aerodynamics] blocks and [fast]  !
      !      (all opt-in; cfg carries the meds_biophysics_opts defaults unless a block overrides).    !
      !      Same types as the column config members, so a plain verbatim struct copy. --------------!
      ctx%col_config%soil_water_opts  = cfg%soil        ! [soil], [fast] -> soil water and the water network
      ctx%soil_albedo(1) = cfg%soil%ground_albedo_vis   ! [soil] ground optics -> the canopy radiation solver
      ctx%soil_albedo(2) = cfg%soil%ground_albedo_nir
      ctx%soil_emiss     = cfg%soil%ground_emissivity
      ctx%col_config%energy = cfg%energy      ! [soil], [fast] -> soil-thermal solver opts
      !----- The initial soil state, from [init]: what init_fast_reservoirs seeds every layer with. --!
      !      A column with a texture per layer starts each layer at its own field capacity. ----------!
      call soil_column_start_theta(cfg%soil_column, ctx%col_config%soil, cfg%init_soil_theta, ctx%theta_init)
      ctx%soil_temp_init = cfg%init_soil_temp
      ctx%col_config%snow   = cfg%snow        ! [snow]         -> snow physical parameter table
      ctx%col_config%aero   = cfg%aero        ! [aerodynamics] -> canopy-aerodynamics constants

      !----- ONE tolerance source drives the whole fast-loop hierarchy (build_tol_set). By default      !
      !      these pushes are the IDENTITY; when [fast].rtol_all > 0 the single master dial propagates  !
      !      into every nested sub-solver as well as the ARK/RK45 march: the soil-plant water network   !
      !      takes the GRP_THETA pair, in water content for soil and relative water content for plant   !
      !      tissues (docs/science/soil_plant_water.md §3). --------------------------------------------!
      ctx%col_config%integrator = build_integrator_opts(cfg)
      associate (tols => ctx%col_config%integrator%error_control%tols)
         ctx%col_config%soil_water_opts%rtol   = tols%rtol(GRP_THETA)
         ctx%col_config%soil_water_opts%atol   = tols%atol(GRP_THETA)
         !----- [fast].soil_energy_budget_atol is a budget-closure threshold, not an integrator   !
         !      tolerance, so atol_scale scales it DIRECTLY instead of through a tolerance group. ---!
         ctx%col_config%energy%atol  = cfg%energy%atol * cfg%atol_scale
      end associate

      !----- Process mask (MEDS_NUMERICS_SCOPING.md §5.1): config logicals -> the mask the schemes    !
      !      honor. All-on = full column. ------------------------------------------------------------!
      ctx%col_config%mask%veg_energy = cfg%mask_veg_energy
      ctx%col_config%mask%cas_energy = cfg%mask_cas_energy
      ctx%col_config%mask%cas_vapour = cfg%mask_cas_vapour
      ctx%col_config%mask%cas_co2    = cfg%mask_cas_co2
      ctx%col_config%mask%soil_heat  = cfg%mask_soil_heat
      ctx%col_config%mask%soil_water = cfg%mask_soil_water
      ctx%col_config%mask%hydraulics = cfg%mask_hydraulics

      !----- Canopy-RT optics table, built ONCE from the [pft] trait table and read-only           !
      !      downstream (meds_fast_config, which the calibration's canopy C API also calls). --------!
      call build_rad_optics(cfg, ctx%rad_opt)
   end subroutine build_fast_context

   !----- Seed every patch's fast reservoirs to a horizontally-uniform equilibrium, and start the  !
   !      patch's cohorts in equilibrium with it, as a recruit is born (patch_physical_state): at the !
   !      canopy-air temperature, holding the water their tissues hold at the top soil layer's     !
   !      potential. Called once after the community is built (production) or before a fast-loop  !
   !      test. ----------------------------------------------------------------------------------!
   subroutine init_fast_reservoirs(site, ctx)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      integer(ik) :: ip, k, nsl, i
      real(wp)    :: psi
      nsl = ctx%col_config%soil%n_active
      do ip = 1_ik, site%patch%n
         associate (cas => site%patch%cas(ip), se => site%patch%soil_energy(ip), sw => site%patch%soil_water(ip))
            sw%theta(1:nsl)  = ctx%theta_init(1:nsl) ; sw%w_surface = 0.0_wp
            sw%w_surface_energy = 0.0_wp        ! dry pond -> zero enthalpy
            do k = 1_ik, nsl
               se%soil_energy(k) = temp_to_internal_energy(ctx%col_config%soil_thermal%soil_dry_heat_capacity(k),    &
                                   ctx%theta_init(k) * rho_h2o, ctx%soil_temp_init, 1.0_wp)
               se%soil_temp(k)   = ctx%soil_temp_init ; se%soil_fliq(k) = 1.0_wp
            end do
            cas%can_shv      = ctx%shv_atm ; cas%can_co2 = ctx%co2_atm
            cas%can_enthalpy = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
            cas%can_temp     = ctx%air_temp
            psi = top_soil_psi(ctx%col_config%soil, sw%theta(1))
            do i = site%patch%cohort_offset(ip), site%patch%cohort_offset(ip) + site%patch%cohort_count(ip) - 1_ik
               site%cohort%leaf_temp(i) = cas%can_temp ; site%cohort%wood_temp(i) = cas%can_temp
               associate (hyd => ctx%col_config%hydraulics_table%pft(site%cohort%pft(i)))
                  call set_cohort_tissue_water(site%cohort, i, psi, hyd%leaf_curve, hyd%wood_curve)
               end associate
            end do
         end associate
      end do
   end subroutine init_fast_reservoirs

   !=======================================================================================!
   !  Advance every patch of the site by one slow-step's worth of fast biophysics:         !
   !  n_fast_per_slow sub-steps over the state-hub reservoirs, on the forcing source's     !
   !  per-sub-step samples (or the context's constant reference climate without one).      !
   !  Optional out-args report the worst whole-column budget residuals + the fail count so !
   !  a caller/test can assert conservation.                                               !
   !=======================================================================================!
   subroutine fast_dynamics(site, ctx, cfg, met_src, met_cur, step_start, worst_energy, worst_water, &
                            n_budget_fail, out_bufs, run_energy_budget, run_water_budget,             &
                            run_face_budget)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t), optional, intent(in)    :: met_src       !< the run's forcing source (read-only here)
      type(met_cursor_t), optional, intent(inout) :: met_cur       !< this polygon's cursor (per-sub-step met)
      type(meds_time_t),  optional, intent(in)    :: step_start    !< calendar time at the START of this slow step
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail
      type(output_buffers_t), optional, intent(inout) :: out_bufs  !< this polygon's FAST-tier staging (when on)
      !----- RUN-level whole-column ledgers: this slow step's per-patch accumulators are area-      !
      !      weighted into a site accumulator and folded in here, so a caller that keeps them across !
      !      the whole run can report the SIGNED cumulative residual (a one-signed bias below the     !
      !      per-step tolerance is invisible to worst_* and n_budget_fail). -------------------------!
      type(budget_t), optional, intent(inout) :: run_energy_budget, run_water_budget
      !----- The per-layer face-closure residual, folded the same way. It is the only one of the   !
      !      three that can see a purely VERTICAL misplacement. -----------------------------------!
      type(budget_t), optional, intent(inout) :: run_face_budget

      !----- The n_fast_per_slow met samples + their sample TIMES, precomputed ONCE per slow      !
      !      step. `t_sub` depends only on `isub`, so met_advance (a FILE READER -- it may reload !
      !      a netCDF bracket) and met_instant (solar geometry + Weiss-Norman disaggregation) are !
      !      site-uniform work; doing it outside the patch loop does it once, and keeps the file  !
      !      reader out of the parallel region. --------------------------------------------------!
      type(met_forcing_t), allocatable :: met_sample(:)
      type(meds_time_t),   allocatable :: t_sample(:)
      !----- Reduction staging: one cell per (accumulator, sub-step, patch), written exactly once  !
      !      inside the patch loop and folded into site%... in patch order afterwards. Sized per   !
      !      call (O(1) allocations per slow step); never zeroed, because every cell is ASSIGNED,  !
      !      not accumulated into. ----------------------------------------------------------------!
      real(wp),    allocatable :: red_site(:,:,:)                  !< (N_RED, sub-step, patch)
      real(wp),    allocatable :: red_worst_energy(:), red_worst_water(:)   !< (patch); max-folded
      type(budget_t), allocatable :: red_budget_energy(:), red_budget_water(:) !< (patch); area-merged
      type(budget_t), allocatable :: red_budget_face(:)                        !< (patch); area-merged
      type(budget_t) :: site_energy_budget, site_water_budget, site_face_budget
      integer(ik), allocatable :: red_nfail(:)                     !< (patch) budget-failure counts
      real(wp)    :: f_ground, sw_ground
      type(met_forcing_t) :: met_ref                  !< the context's reference climate (no forcing)
      integer(ik) :: ip, isub, npatch, ncoh_max, nsub, nl, n_thread
      logical     :: do_forcing, do_fast, do_cdiag, do_pdiag
      !----- Per-(cohort, sub-step) diagnostic scratch: filled by the pre-pass through cdiag, then   !
      !      folded (dt-weighted) into site%cohort%diag. PER THREAD, because the patch loop is         !
      !      parallel -- it rides the same per-thread pool discipline as the rest of the scratch.  ---!
      real(wp),    allocatable :: cdiag_pool(:,:,:)                !< (N_CDIAG, cohort, thread)
      !=========================================================================================!
      !  The PER-THREAD scratch POOL. One copy of every per-patch working buffer per thread,    !
      !  indexed by thread id. The patch loop aliases its own slice with `associate`, so the    !
      !  loop body reads as serial code and nothing in it can accidentally reach another        !
      !  thread's buffer.                                                                       !
      !                                                                                         !
      !  A POOL, and not the obvious OpenMP data-sharing clauses, because BOTH of those are     !
      !  broken for these types on ifx 2026 and each fails SILENTLY in its own way              !
      !  (docs/building.md, the portability traps):                                             !
      !    - a BLOCK inside the region: ifx does not default-initialise block-scoped derived    !
      !      types in an outlined region, so their allocatable descriptors hold garbage,        !
      !      ensure_*_capacity reads them as already-allocated, and the first deallocate        !
      !      SIGSEGVs -- even at n_threads = 1;                                                 !
      !    - `private(...)`: ifx constructs each thread's default-initialised copy through a    !
      !      compiler-generated STATIC mold (`AERO_OUT_T.omp.mold_ctor`) that every thread      !
      !      WRITES;                                                                            !
      !    - `firstprivate(...)`: the same, one step over, in `FAST_CONTEXT_T.omp.copy_ctor`.   !
      !  The last two are real data races (ThreadSanitizer finds them) and they produce exactly !
      !  the failure the deterministic fold exists to prevent: plausible-looking answers that   !
      !  move with the thread count. A plain array indexed by thread id involves no             !
      !  compiler-generated constructor at all.                                                 !
      !                                                                                         !
      !  Each entry is sized ONCE, to the site-wide max cohort count: ensure_*_capacity runs    !
      !  once per thread (not once per patch), so allocations are O(n_thread) per slow step,    !
      !  O(1) at the default 1 thread. Every reader loops by the ACTIVE count (col_cohort%n /   !
      !  ncoh), never by size(...), so reusing a larger patch's capacity for a smaller one is   !
      !  exact.                                                                                 !
      !=========================================================================================!
      type(column_cohort_t),  allocatable :: coh_pool(:)
      type(column_forcing_t), allocatable :: forc_pool(:)
      type(aero_env_t),       allocatable :: aenv_pool(:)
      type(aero_geom_t),      allocatable :: ageom_pool(:)
      type(aero_out_t),       allocatable :: aero_pool(:)
      type(patch_biophys_t),  allocatable :: bio_pool(:)
      type(column_budget_t),  allocatable :: budg_pool(:)
      type(fast_workspace_t), allocatable :: ws_pool(:)        !< the fast step's reused storage
      type(met_forcing_t),    allocatable :: met_pool(:)
      type(met_forcing_t),    allocatable :: met_top_pool(:)   !< the sample at the patch's canopy-air top
      real(wp),               allocatable :: gpp_pool(:,:), leaf_resp_pool(:,:)
      real(wp),               allocatable :: stem_resp_pool(:,:), root_resp_pool(:,:), psi_leaf_pool(:,:)
      real(wp)    :: sum_lai, sum_wai, le_flux, h_flux, rnet, gpp_patch, npp_patch, dt_fast_days
      real(wp)    :: prow(N_PDIAG)                !< the patch block's row for one (patch, sub-step)
      real(wp)    :: z_top, rough_p, displace_p   !< [m] this patch's canopy-air top, roughness, displacement
      integer(ik) :: j, i, i0, ncoh, ith, kb
      real(wp)    :: rain_s, snow_s, echo(N_PYDIAG)   !< the forcing echo's precipitation split, per sub-step

      !----- Live forcing drives the fast loop only when it is ON and a reader + step time are   !
      !      supplied; otherwise every sub-step reads the context's constant reference climate   !
      !      (met_ref). The diurnal cycle lives INSIDE the sub-step loop                         !
      !      (docs/science/order_of_processes.md §3).                                            !
      do_forcing = cfg%forcing%forcing_on .and. present(met_src) .and. present(met_cur)          &
                   .and. present(step_start)
      f_ground = ctx%rad_sw_ground / max(ctx%rad_sw_top, tiny_num)   ! ground/canopy-top SW transmittance
      met_ref  = reference_met(ctx)
      npatch   = site%patch%n
      nsub     = cfg%n_fast_per_slow
      n_thread = max(1_ik, ctx%patch_threads)

      !----- Reset the fast->slow carbon accumulators (gross GPP + maintenance-resp losses) BEFORE !
      !      the fast window (compute_carbon_allocation reads them after; it has site intent(in),   !
      !      cannot reset).                                                                          !
      site%cohort%gpp_accum(1:site%cohort%n)       = 0.0_wp
      site%cohort%leaf_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%stem_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%root_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%light_hours_accum(1:site%cohort%n) = 0.0_wp
      !----- ROLL OVER the predawn water status: today's accumulated maximum becomes the value the !
      !      leaf kernel's slow stomatal control uses tomorrow (docs/science/leaf_gas_exchange.md  !
      !      §4.2), then the accumulator restarts. A cohort whose max is still the UNSET sentinel  !
      !      (a recruit born mid-day) keeps it, so column_prepass seeds it from the soil rather    !
      !      than inheriting a meaningless value. -------------------------------------------------!
      !----- Roll over only where the day actually produced a sample (a cohort born mid-day, or one  !
      !      culled before its first step, has none). The guard is against the RESET value, not 0:    !
      !      psi_leaf is <= 0, so '<= 0' would also accept an untouched accumulator. -----------------!
      where (site%cohort%dmax_psi_leaf_accum(1:site%cohort%n) > 0.5_wp * DMAX_PSI_LEAF_ACCUM_RESET)          &
         site%cohort%dmax_psi_leaf(1:site%cohort%n) = site%cohort%dmax_psi_leaf_accum(1:site%cohort%n)
      site%cohort%dmax_psi_leaf_accum(1:site%cohort%n) = DMAX_PSI_LEAF_ACCUM_RESET
      !----- Reset the daily fast->slow soil-carbon accumulator (opt-in                             !
      !      [soil_carbon].soil_carbon_on -- a harmless no-op accumulation when off, since          !
      !      column_prepass leaves budget%xi_step/rh_matrix_step at 0 in that case). ---------------!
      if (cfg%soil_carbon_on) site%patch%xi_accum(1:site%patch%n) = xi_accum_t()
      !----- The site daily-mean air-temperature accumulator (which the slow-loop phenology driver !
      !      reads AFTER this fast window), site ET, and the integrator WORK counters are all      !
      !      reset and refilled by the patch-order fold AFTER the patch loop -- see the bottom of  !
      !      this routine. Air temperature is site-uniform, so accumulating once per (patch,       !
      !      sub-step) and dividing by the count still yields the daily mean. ---------------------!

      !----- FAST (sub-daily) output staging: fill it only when the tier is active and a diurnal    !
      !      signal exists (forcing on). Serialize stays in main.                                    !
      do_fast = present(out_bufs) .and. do_forcing
      if (do_fast) do_fast = out_bufs%fast_on
      !----- The per-cohort / per-patch diagnostic capture runs only when the registry actually   !
      !      asked for it (main sets `active` from the live variable list). Everything downstream !
      !      is gated on these two flags, so a run that reports no ecophysiology pays nothing for !
      !      it. ---------------------------------------------------------------------------------!
      do_cdiag = site%cohort%diag%active
      do_pdiag = site%patch%diag%active
      !----- Size the diagnostic accumulators to this step's stand before the patch loop folds into    !
      !      them. They are NOT zeroed here: they run from one output tick to the next, and the tick's !
      !      caller zeroes them after reading (reset_step_diagnostics), so a calendar boundary's       !
      !      events, written between the two steps, reach the next record. ----------------------------!
      if (do_cdiag) then
         call cohort_diag_grow(site%cohort%diag, max(site%cohort%n, 1_ik))
         site%cohort%diag%n = site%cohort%n
      end if
      if (do_pdiag) then
         call patch_diag_grow(site%patch%diag, max(npatch, 1_ik))
         site%patch%diag%n = npatch
      end if
      nl = n_soil_layer_max
      if (do_fast) then
         call size_fast_staging(out_bufs, nsub, npatch, site%cohort%n, do_cdiag)
         out_bufs%n_fast_sub    = nsub
         !----- The slabs are sized to the ceiling, but only the ACTIVE layers are data: the layers !
         !      past n_active stay unwritten and read back as the fill value, as on the coarse      !
         !      tiers. -----------------------------------------------------------------------------!
         out_bufs%fast_n_soil   = min(ctx%col_config%soil%n_active, nl)
         out_bufs%fast_n_cohort = site%cohort%n
         out_bufs%fast_n_patch  = npatch
         out_bufs%fast_site = 0.0_wp
         out_bufs%fast_soil_temp = 0.0_wp ; out_bufs%fast_soil_water = 0.0_wp
         out_bufs%fast_coh_ltemp = 0.0_wp ; out_bufs%fast_coh_gpp = 0.0_wp ; out_bufs%fast_coh_height = 0.0_wp
         if (do_cdiag) out_bufs%fast_coh_cdiag = 0.0_wp
      end if

      ncoh_max = 0_ik
      do ip = 1_ik, npatch
         ncoh_max = max(ncoh_max, site%patch%cohort_count(ip))
      end do

      !----- Per-thread scratch for the per-cohort diagnostic capture. Sized from ncoh_max, so it   !
      !      MUST follow the loop just above that computes it.  --------------------------------!
      allocate(cdiag_pool(N_CDIAG, max(ncoh_max,1_ik), max(n_thread,1_ik)))

      !----- The per-(sub-step, patch) reduction staging. Allocated (not zeroed -- every cell is   !
      !      assigned) per call, so this adds O(1) allocations per slow step, not O(n_patch). -----!
      allocate(red_site(N_RED, nsub, max(npatch,1_ik)))
      allocate(red_budget_energy(max(npatch,1_ik)), red_budget_water(max(npatch,1_ik)))
      allocate(red_budget_face(max(npatch,1_ik)))
      allocate(red_worst_energy(max(npatch,1_ik)), red_worst_water(max(npatch,1_ik)),               &
               red_nfail(max(npatch,1_ik)))

      !----- Precompute this slow step's met samples. met_advance is called here in strictly        !
      !      increasing t exactly once per sub-step, and met_instant is a pure function of (reader  !
      !      state, t). The sample TIMES are kept per sub-step too: the probe writes them.          !
      if (do_forcing) then
         allocate(met_sample(nsub), t_sample(nsub))
         do isub = 1_ik, nsub
            t_sample(isub) = time_advance_seconds(step_start,                                          &
                       (real(isub, wp) - 1.0_wp + cfg%forcing_sample_frac) * cfg%dt_fast)
            call met_advance(met_src, met_cur, t_sample(isub))
            met_sample(isub) = met_instant(met_src, met_cur, t_sample(isub))
         end do
      end if
      !----- The FAST tier stamps a record by its window's start, as the slower tiers do: each      !
      !      sub-step's start, not the instant its forcing is sampled at. Site-uniform, so out of   !
      !      the patch loop. -----------------------------------------------------------------------!
      if (do_fast) then
         do isub = 1_ik, nsub
            out_bufs%fast_time(isub) = time_advance_seconds(step_start, real(isub - 1_ik, wp) * cfg%dt_fast)
         end do
      end if

      !=========================================================================================!
      !  The parallel patch loop (docs/science/numerical_scheme.md §6a). Patch columns are      !
      !  independent within a dt_fast (they couple only through the slow loop), so patches are  !
      !  the parallel axis; `schedule(dynamic,1)` because the adaptive march makes them WILDLY  !
      !  unequal (one collapsed patch can run 13x the others), which a static schedule would    !
      !  serialize behind.                                                                      !
      !                                                                                         !
      !  Without OpenMP flags every directive is a comment and the pool is one deep, so this is !
      !  the serial code.                                                                       !
      !=========================================================================================!
      allocate(coh_pool(n_thread), forc_pool(n_thread), aenv_pool(n_thread), ageom_pool(n_thread),  &
               aero_pool(n_thread), bio_pool(n_thread), budg_pool(n_thread), ws_pool(n_thread),     &
               met_pool(n_thread), met_top_pool(n_thread))
      allocate(gpp_pool(ncoh_max, n_thread), leaf_resp_pool(ncoh_max, n_thread),                    &
               stem_resp_pool(ncoh_max, n_thread), root_resp_pool(ncoh_max, n_thread),              &
               psi_leaf_pool(ncoh_max, n_thread))
      do ith = 1_ik, n_thread
         met_pool(ith) = met_ref                            ! overwritten per sub-step when forcing is on
         met_top_pool(ith) = met_ref
         call ensure_column_cohort_capacity(coh_pool(ith), ncoh_max)
         call ensure_patch_biophys_capacity(bio_pool(ith), ncoh_max, ctx%air_temp, ctx%shv_atm,     &
                                            ctx%co2_atm, ctx%air_temp)
         call ensure_aero_out_capacity(aero_pool(ith), ncoh_max)
         call alloc_forcing(forc_pool(ith), ncoh_max)
      end do

      !$omp parallel do default(shared) schedule(dynamic, 1) num_threads(n_thread)                  &
      !$omp    private(ip, ith, isub, j, i, i0, ncoh,                                              &
      !$omp            sum_lai, sum_wai, le_flux, h_flux, rnet, gpp_patch, npp_patch, dt_fast_days,   &
      !$omp            prow, kb,                                                                    &
      !$omp            sw_ground, z_top, rough_p, displace_p)
      do ip = 1_ik, npatch
         !----- This thread's slot in the scratch pool. The `!$` sentinel keeps the non-OpenMP build  !
         !      on slot 1 with no dependence on omp_lib. ---------------------------------------------!
         ith = 1_ik
         !$ ith = int(omp_get_thread_num(), ik) + 1_ik
         associate (col_cohort           => coh_pool(ith),        forc          => forc_pool(ith),         &
                    aenv          => aenv_pool(ith),       ageom         => ageom_pool(ith),        &
                    aero          => aero_pool(ith),       biophys           => bio_pool(ith),          &
                    budget          => budg_pool(ith),       fast_ws       => ws_pool(ith),          &
                    met           => met_pool(ith),        met_top       => met_top_pool(ith),      &
                    gpp_coh       => gpp_pool(:,ith),                                               &
                    leaf_resp_coh => leaf_resp_pool(:,ith), stem_resp_coh => stem_resp_pool(:,ith), &
                    root_resp_coh => root_resp_pool(:,ith), psi_leaf_coh => psi_leaf_pool(:,ith), &
                    cdiag_buf     => cdiag_pool(:,:,ith))
         ncoh = site%patch%cohort_count(ip)
         i0   = site%patch%cohort_offset(ip)

         !----- Gather the patch's cohort slice into the column buffer (+ the per-ground indices it !
         !      forms). Capacity was ensured above (ncoh <= ncoh_max always); this just updates the !
         !      ACTIVE count -- no allocation. -----------------------------------------------------!
         call copy_column_cohort(col_cohort, site%cohort, i0, ncoh)
         sum_lai = 0.0_wp ; sum_wai = 0.0_wp
         do j = 1_ik, ncoh
            sum_lai = sum_lai + col_cohort%lai(j)
            sum_wai = sum_wai + col_cohort%wai(j)
         end do

         !----- Per-patch canopy geometry. The ground conductance weighs bare ground by the patch's  !
         !      open fraction exp(-(LAI + WAI)) (docs/science/canopy_aerodynamics.md eq. (6)), so a   !
         !      patch without cohorts couples its ground as bare soil. The canopy is never buried in  !
         !      snow (snowfac = 0; see docs/ROADMAP.md). ------------------------------------------------!
         ageom%veg_height   = ctx%veg_height_bare
         do j = 1_ik, ncoh
            ageom%veg_height = max(ageom%veg_height, col_cohort%height(j))
         end do
         ageom%opencan_frac = open_canopy_fraction(sum_lai, sum_wai)
         ageom%snowfac      = 0.0_wp

         call alloc_forcing(forc, ncoh)

         !----- Assemble the working bundle: adopt the owned per-patch reservoirs AND the        !
         !      PERSISTED per-cohort leaf_temp/leaf_water_mass carried on the cohort block (no     !
         !      reseeding). Capacity was ensured above; every field below is unconditionally        !
         !      (re)assigned from the site, so no alloc_patch_biophys seed call is needed here. -----!
         biophys%cas    = site%patch%cas(ip)
         biophys%soil_energy = site%patch%soil_energy(ip)
         biophys%soil_water = site%patch%soil_water(ip)
         biophys%snow   = site%patch%snow(ip)
         biophys%adapt_dt_last = site%patch%adapt_dt_last(ip)   ! per-patch, not loop-carried
         !----- FROZEN slow soil-carbon pool: a read-only snapshot for TODAY, held constant across  !
         !      the sub-step loop below (never written back -- the daily soil_carbon_step is the    !
         !      sole writer of the real site-level pool). The pool is ALWAYS allocated (0 per patch !
         !      at creation), so this copy is safe unconditionally: with soil_carbon_on = .false.   !
         !      nothing ever writes it, and it stays 0.                                             !
         biophys%soil_carbon = site%patch%soil_carbon(ip)
         !----- FROZEN daily leaf/root-turnover shed-water rate: the same "read-only snapshot for   !
         !      TODAY, held constant across the sub-step loop" convention as soil_carbon just       !
         !      above. -----------------------------------------------------------------------------!
         biophys%shed_water_rate = site%patch%shed_water_rate(ip)
         biophys%slow_co2_rate   = site%patch%slow_co2_rate(ip)
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            biophys%leaf_temp(j) = site%cohort%leaf_temp(i)
            biophys%wood_temp(j) = site%cohort%wood_temp(i)
            !----- Tissue water is READ here, never written. Fitting it to the capacity the biomass  !
            !      allows is the slow tier's: a cohort is born with its water, new leaves take theirs   !
            !      from the wood, and water above the capacity is shed before the fast loop runs       !
            !      (docs/science/soil_plant_water.md §2), all booked by the slow ledger. ---------------!
            biophys%leaf_water_mass(j) = site%cohort%leaf_water_mass(i)
            biophys%wood_water_mass(j) = site%cohort%wood_water_mass(i)
            !----- Surface (interception film) water needs no lazy-init seed: 0 (bone dry) is a real  !
            !      initial condition here, not a placeholder -- a freshly-created cohort simply hasn't  !
            !      been rained on yet. -----------------------------------------------------------------!
            biophys%leaf_surf_water(j) = site%cohort%leaf_surf_water(i)
            biophys%wood_surf_water(j) = site%cohort%wood_surf_water(i)
         end do

         call ensure_aero_out_capacity(aero, ncoh)

         !----- This patch's canopy-air top and roughness. Every forcing sample is moved to that top  !
         !      along the patch's own log profile (meds_lapse_rate, docs/science/forcing.md §8), and    !
         !      the aerodynamics runs from there. Both are fixed for the day: the canopy height and the !
         !      canopy-air depth change only on a slow step. -------------------------------------------!
         z_top = biophys%cas%can_depth
         call canopy_roughness(ctx%col_config%aero, ageom, rough_p, displace_p)

         !----- n_fast_per_slow sub-steps. Forcing is re-evaluated PER SUB-STEP (the diurnal     !
         !      cycle lives here): take the sub-step's forcing record, then fill_forcing +       !
         !      fill_aenv from it. CONSTANT path (do_forcing=.false.): the record is the         !
         !      context's reference climate (met_ref) and the ground shortwave its own reference !
         !      value.                                                                           !
         budget = column_budget_t()
         do isub = 1_ik, cfg%n_fast_per_slow
            !----- The met sample point within the sub-step ([fast].forcing_sample_frac, default 0.5 = !
            !      midpoint) is applied when met_sample is built above. The ground shortwave scales    !
            !      with the canopy-top shortwave at the reference transmittance f_ground (the constant !
            !      path keeps the reference value itself). `met` and `sw_ground` are per-THREAD. ------!
            if (do_forcing) then
               met = met_sample(isub)
               met_top = met_to_cas_top(met, cfg%forcing, z_top, displace_p, rough_p)
               sw_ground = f_ground * met%swdown()
            else
               met_top = met                   ! the reference climate IS the air at the canopy-air top
               sw_ground = ctx%rad_sw_ground
            end if
            !----- Accumulate the sub-step air temperature for the daily-mean phenology driver. ------!
            red_site(RED_PHENO_TAIR, isub, ip) = met%tair_k
            call fill_forcing(forc, col_cohort, met_top, sw_ground, sum_lai)
            red_site(RED_RAINF,    isub, ip) = site%patch%area(ip) * forc%rainfall
            red_site(RED_SNOWFALL, isub, ip) = site%patch%area(ip) * forc%snowfall
            !----- RT join (docs/science/canopy_radiation_transfer.md §4): when forcing is on,        !
            !      REPLACE the LAI-share SW split with real per-cohort absorbed SW/PAR and net LW     !
            !      from the two-stream canopy radiation, which emits each cohort at its radiative     !
            !      temperature from the start-of-sub-step tissue temperatures (see                    !
            !      apply_rt_forcing). ----------------------------------------------------------------!
            if (do_forcing) call apply_rt_forcing(forc, ncoh, col_cohort%pft, col_cohort%lai, col_cohort%wai,  &
                                 col_cohort%height, biophys%leaf_temp, biophys%wood_temp,                      &
                                 biophys%soil_energy%soil_temp(1), biophys%snow, ctx%col_config%snow,               &
                                 ctx%soil_albedo, ctx%soil_emiss, ctx%rad_opt, met, cfg%leaf_absorptance)
            call fill_aenv(aenv, biophys, met_top, z_top)
            !----- Slice to 1:ncoh (not the whole, possibly capacity-oversized backing array): the    !
            !      five per-cohort accumulators are assumed-shape dummies in column_fast_step, so the !
            !      ACTUAL argument's extent must equal col_cohort%n exactly, whatever the backing     !
            !      array's capacity (the pool sizes it to the site-wide max). ------------------------!
            !----- The per-cohort DIAGNOSTIC capture is passed only when the run reports per-cohort   !
            !      ecophysiology. Absent, the pre-pass never even asks the leaf kernel for the extra   !
            !      leaf_flux_t fields, so a production run pays nothing for a feature it is not using. !
            if (do_cdiag) then
               cdiag_buf(:, 1:ncoh) = 0.0_wp
               call column_fast_step(cfg%dt_fast, cfg, ctx%col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                     fast_ws, gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux, cdiag=cdiag_buf(:, 1:ncoh))
            else
               call column_fast_step(cfg%dt_fast, cfg, ctx%col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                     fast_ws, gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux)
            end if
            !----- Integrate the area-weighted CAS->atm latent flux -> site ET [kg/m2 = mm] over the step. !
            red_site(RED_ET, isub, ip) = site%patch%area(ip) * (le_flux / latent_heat_vap) * cfg%dt_fast
            !----- And the water leaving its soil column: surface runoff and bottom drainage [kg/m2]. ----!
            red_site(RED_RUNOFF_SURF, isub, ip) = site%patch%area(ip) * budget%runoff_surf_export
            red_site(RED_DRAINAGE,    isub, ip) = site%patch%area(ip) * budget%drainage_export
            !----- Integrator WORK counters: area-weight like every other site diagnostic, so a patch  !
            !      that needs more sub-steps is not double-counted by its area share. -----------------!
            red_site(RED_INTEG_STEPS,  isub, ip) = site%patch%area(ip) * real(budget%integ_nsteps,   wp)
            red_site(RED_INTEG_REJ,    isub, ip) = site%patch%area(ip) * real(budget%integ_nrej,     wp)
            red_site(RED_WATER_NSUB,   isub, ip) = site%patch%area(ip) * real(budget%water_nsub,   wp)
            red_site(RED_WATER_NEWTON, isub, ip) = site%patch%area(ip) * real(budget%water_newton, wp)
            red_site(RED_WATER_REJECT, isub, ip) = site%patch%area(ip) * real(budget%water_reject, wp)
            red_site(RED_WATER_NFAIL,  isub, ip) = site%patch%area(ip) * real(budget%water_nfail,  wp)
            !----- INTEGRATOR HEALTH, same area weighting. work_rk45_rescue is the one that decides     !
            !      whether an "RK45 run" was actually RK45: a nonzero total means some dt_fast steps    !
            !      were rescued by ARK instead, which changes what the run measures.                    !
            red_site(RED_RK45_RESCUE,  isub, ip) = site%patch%area(ip) * real(budget%rk45_rescue,    wp)
            red_site(RED_CLAMP_STAGE,  isub, ip) = site%patch%area(ip) * real(budget%clamp_stage_n,  wp)
            !----- Integrate this sub-step's per-pool env scalar + matrix Rh into the day's totals:   !
            !      dt_fast_days converts the instantaneous xi_step/rh_matrix_step (column_prepass     !
            !      leaves both at 0 when soil_carbon_on=.false.) into the day-integral xi_int the     !
            !      daily soil_carbon_step consumes (docs/science/soil_carbon.md §6), and the          !
            !      accumulated Rh the audit cross-checks against.                                     !
            !                                                                                         !
            !      NOTE: dt_fast_days is declared at routine scope and listed `private`, not in a     !
            !      BLOCK. nvfortran 25.11 rejects a BLOCK construct anywhere inside a parallel region !
            !      ("Unimplemented feature"), and MEDS must build on all three back ends. ------------!
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
                                     le_flux, met%swdown())
            !----- Patch GPP and net absorbed radiation, once, for the patch block's row. h_flux is    !
            !      surfaced BY column_fast_step (like le_flux), computed from the in-call aero/aenv    !
            !      state: computing it here from POST-call reads miscompiles to 0 on nvfortran (the    !
            !      same class as the array-valued-result trap in docs/building.md). le_flux and rnet   !
            !      read intent(in) forc, so they are safe. --------------------------------------------!
            gpp_patch = sum(gpp_coh(1:ncoh) * col_cohort%nplant(1:ncoh))    ! [umol/m2/s]
            rnet      = forc%abs_sw_ground + forc%abs_lw_ground                                          &
                        + sum(forc%abs_sw(1:ncoh)) + sum(forc%abs_lw(1:ncoh))
            !----- The patch block's row for this sub-step: the patch diagnostics accumulate it, and the  !
            !      FAST tier stages it as it is, per patch; its site means are the area-weighted sums     !
            !      after the patch loop. NPP is GPP net of the three MAINTENANCE respiration terms only    !
            !      (growth respiration is charged in the slow allocator, and adding it here would double-  !
            !      count against npp_*_site). --------------------------------------------------------!
            if (do_pdiag .or. do_fast) then
               npp_patch = gpp_patch - sum((leaf_resp_coh(1:ncoh) + stem_resp_coh(1:ncoh)             &
                                            + root_resp_coh(1:ncoh)) * col_cohort%nplant(1:ncoh))
               call patch_diag_row(prow, cfg%dt_fast, le_flux, h_flux, rnet,                         &
                                   forc%abs_sw_ground, forc%abs_lw_ground,                           &
                                   aero%ustar, aero%ggnet, aero%rough, aero%displace,                &
                                   biophys%cas%can_temp, biophys%cas%can_shv, biophys%cas%can_co2,   &
                                   gpp_patch, budget%nee_last, npp_patch,                            &
                                   biophys%soil_energy%soil_temp(1), budget%whole_energy%resid,           &
                                   budget%whole_water%resid, met_top%wind, met_top%tair_k, z_top,    &
                                   forc%sw_in_vis, forc%sw_in_nir, forc%sw_up_vis, forc%sw_up_nir,   &
                                   forc%lw_up, biophys%soil_water%w_surface)
            end if
            !----- FAST (sub-daily) output staging, by the patch's own slot (and the per-cohort slabs by  !
            !      the site's cohort slot), both DISJOINT across patches, so no reduction is needed here.  !
            !      main replays this into the FAST buffers (output_integrate_fast). ----------------------!
            if (do_fast) then
               out_bufs%fast_patch(:, ip, isub) = prow
               kb = (ip - 1_ik) * nl
               out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)  = biophys%soil_energy%soil_temp(1:nl)
               out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub) = biophys%soil_water%theta(1:nl)
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  out_bufs%fast_coh_ltemp(i,isub)  = biophys%leaf_temp(j)
                  out_bufs%fast_coh_gpp(i,isub)    = gpp_coh(j)
                  out_bufs%fast_coh_height(i,isub) = col_cohort%height(j)
                  if (do_cdiag) out_bufs%fast_coh_cdiag(i, :, isub) = cdiag_buf(:, j)
               end do
            end if
            !----- FOLD the per-(cohort, sub-step) and per-patch DIAGNOSTICS into the site's         !
            !      dt-weighted accumulators (meds_diag_vars). Sub-daily resolution exists ONLY       !
            !      here, so this is where it is captured.                                            !
            !                                                                                        !
            !      THREAD SAFETY: cohort slots [i0:i1] and patch slot ip are DISJOINT across patches !
            !      (the CSR map partitions the flat SoA), so these writes need no reduction and stay !
            !      byte-identical at any thread count -- the same argument the per-cohort state      !
            !      write-back below relies on. ------------------------------------------------------!
            if (do_cdiag) then
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  site%cohort%diag%v(:, i) = site%cohort%diag%v(:, i) + cdiag_buf(:, j) * cfg%dt_fast
                  site%cohort%diag%w(i)    = site%cohort%diag%w(i)    + cfg%dt_fast
               end do
            end if
            if (do_pdiag) call accumulate_patch_diag(site%patch%diag, ip, cfg%dt_fast, prow)
            !----- Integrate GROSS GPP + maintenance-resp losses [umol/plant/s] -> [kgC/plant].  !
            !      Keep gross and loss terms SEPARATE (compute_carbon_allocation nets them; mirrors ED2). !
            do j = 1_ik, ncoh
               i = i0 + j - 1_ik
               site%cohort%gpp_accum(i)       = site%cohort%gpp_accum(i)       + gpp_coh(j)       * cfg%dt_fast * umol_2_kgC
               site%cohort%leaf_resp_accum(i) = site%cohort%leaf_resp_accum(i) + leaf_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%stem_resp_accum(i) = site%cohort%stem_resp_accum(i) + stem_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%root_resp_accum(i) = site%cohort%root_resp_accum(i) + root_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               !----- The hours of light at the cohort's top over the slow step, the phenology light   !
               !      cue's driver: this step counts when its PAR exceeds the PFT's par_min, so each      !
               !      cohort sees its own clouds and its own place in the canopy. ----------------------!
               if (forc%par_top(j) * par_w_2_umol > cfg%pft%pheno_par_min(site%cohort%pft(i)))           &
                  site%cohort%light_hours_accum(i) = site%cohort%light_hours_accum(i) + cfg%dt_fast / 3600.0_wp
               !----- Running daily MAX of psi_leaf. max(), not a sum: the daily maximum occurs near !
               !      dawn and IS the quantity that drives tomorrow's beta_stomata                   !
               !      (docs/science/leaf_gas_exchange.md §4.2). -------------------------------------!
               site%cohort%dmax_psi_leaf_accum(i) = max(site%cohort%dmax_psi_leaf_accum(i), psi_leaf_coh(j))
            end do
         end do

         !----- Write the evolved state back to the site: per-patch reservoirs + per-cohort water. !
         site%patch%cas(ip)    = biophys%cas
         site%patch%soil_energy(ip) = biophys%soil_energy
         site%patch%soil_water(ip) = biophys%soil_water
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
         red_budget_face(ip)   = budget%soil_face_mass
         red_nfail(ip)        = budget%whole_energy%n_fail + budget%whole_water%n_fail
         end associate
      end do
      !$omp end parallel do

      !=========================================================================================!
      !  Fold the staged per-(sub-step, patch) contributions back into the site in patch order: !
      !  patch outer, sub-step inner, the order a serial loop accumulates in. Floating-point    !
      !  addition is not associative, so this ordering IS the guarantee: it makes every site    !
      !  diagnostic independent of the thread count, and equal bit-for-bit to the serial fold   !
      !  (docs/science/numerical_scheme.md §6a).                                                !
      !=========================================================================================!
      site%pheno_tair_sum = 0.0_wp ; site%pheno_tair_n = 0_ik
      site%et_accum       = 0.0_wp ; site%runoff_surf_accum = 0.0_wp ; site%drainage_accum = 0.0_wp
      site%work_integ_steps = 0.0_wp ; site%work_integ_rej  = 0.0_wp
      site%work_water_nsub  = 0.0_wp ; site%work_water_newton = 0.0_wp
      site%work_water_nfail = 0.0_wp ; site%work_water_reject = 0.0_wp
      site%work_rk45_rescue = 0.0_wp ; site%work_clamp_stage  = 0.0_wp
      do ip = 1_ik, npatch
         do isub = 1_ik, nsub
            site%et_accum          = site%et_accum          + red_site(RED_ET,            isub, ip)
            site%runoff_surf_accum = site%runoff_surf_accum + red_site(RED_RUNOFF_SURF,   isub, ip)
            site%drainage_accum    = site%drainage_accum    + red_site(RED_DRAINAGE,      isub, ip)
            site%pheno_tair_sum    = site%pheno_tair_sum    + red_site(RED_PHENO_TAIR,    isub, ip)
            site%pheno_tair_n      = site%pheno_tair_n      + 1_ik
            site%work_integ_steps  = site%work_integ_steps  + red_site(RED_INTEG_STEPS,   isub, ip)
            site%work_integ_rej    = site%work_integ_rej    + red_site(RED_INTEG_REJ,     isub, ip)
            site%work_water_nsub   = site%work_water_nsub   + red_site(RED_WATER_NSUB,   isub, ip)
            site%work_water_newton = site%work_water_newton + red_site(RED_WATER_NEWTON, isub, ip)
            site%work_water_reject = site%work_water_reject + red_site(RED_WATER_REJECT, isub, ip)
            site%work_water_nfail  = site%work_water_nfail  + red_site(RED_WATER_NFAIL,  isub, ip)
            site%work_rk45_rescue  = site%work_rk45_rescue  + red_site(RED_RK45_RESCUE,   isub, ip)
            site%work_clamp_stage  = site%work_clamp_stage  + red_site(RED_CLAMP_STAGE,   isub, ip)
         end do
      end do
      !----- FAST-tier site means: each patch's staged row and soil column, area-weighted, folded over  !
      !      the patches in increasing order for each sub-step. -------------------------------------!
      if (do_fast) then
         do isub = 1_ik, nsub
            do ip = 1_ik, npatch
               kb = (ip - 1_ik) * nl
               out_bufs%fast_site(:, isub) = out_bufs%fast_site(:, isub)                                  &
                                             + site%patch%area(ip) * out_bufs%fast_patch(:, ip, isub)
               out_bufs%fast_soil_temp(1:nl, isub)  = out_bufs%fast_soil_temp(1:nl, isub)                 &
                                             + site%patch%area(ip) * out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)
               out_bufs%fast_soil_water(1:nl, isub) = out_bufs%fast_soil_water(1:nl, isub)                &
                                             + site%patch%area(ip) * out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub)
            end do
         end do
      end if
      !----- The forcing echo (PY_*, MEDS_FORCING_DESIGN.md §6.7), once per sub-step: the          !
      !      sub-step's sample (the reference climate without a source) with its precipitation     !
      !      split as the patches received it, at their canopy-air tops, area-summed in patch      !
      !      order. The FAST tier stages it as is and the polygon block accumulates it for the     !
      !      daily, monthly and yearly tiers. -----------------------------------------------------!
      if (do_fast .or. site%diag%active) then
         do isub = 1_ik, nsub
            rain_s = 0.0_wp ; snow_s = 0.0_wp
            do ip = 1_ik, npatch
               rain_s = rain_s + red_site(RED_RAINF,    isub, ip)
               snow_s = snow_s + red_site(RED_SNOWFALL, isub, ip)
            end do
            if (do_forcing) then
               call forcing_echo(met_sample(isub), rain_s, snow_s, echo)
            else
               call forcing_echo(met_ref, rain_s, snow_s, echo)
            end if
            if (do_fast) out_bufs%fast_forcing(:, isub) = echo
            if (site%diag%active) call accumulate_polygon_diag(site%diag, cfg%dt_fast, echo)
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
      if (present(run_energy_budget) .or. present(run_water_budget) .or. present(run_face_budget)) then
         site_energy_budget = budget_t() ; site_water_budget = budget_t() ; site_face_budget = budget_t()
         do ip = 1_ik, npatch
            call budget_merge(site_energy_budget, red_budget_energy(ip), site%patch%area(ip))
            call budget_merge(site_water_budget,  red_budget_water(ip),  site%patch%area(ip))
            call budget_merge(site_face_budget,   red_budget_face(ip),   site%patch%area(ip))
         end do
         if (present(run_energy_budget)) call budget_merge(run_energy_budget, site_energy_budget, 1.0_wp)
         if (present(run_water_budget))  call budget_merge(run_water_budget,  site_water_budget,  1.0_wp)
         if (present(run_face_budget))   call budget_merge(run_face_budget,   site_face_budget,   1.0_wp)
      end if
      if (do_fast) out_bufs%fast_ready = .true.   ! signal main to replay + serialize the FAST tier
   end subroutine fast_dynamics

   !----- Sub-daily diagnostic PROBE. Opt-in ([fast].fast_probe): one CSV row per (patch, sub-step) !
   !      with the fast-loop state an integrator or dt_fast evaluation resolves on -- shortwave in, !
   !      CAS temp, patch GPP, latent flux, top-soil temp, mean leaf and wood temp. A saved unit    !
   !      opens the file on first use (header) and appends thereafter; the program-exit close       !
   !      flushes it. NOT part of the output subsystem, whose FAST tier is the supported sub-daily  !
   !      output. ----------------------------------------------------------------------------------!
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
            sw_in, biophys%cas%can_temp, gpp_patch, le_flux, biophys%soil_energy%soil_temp(1), leaf_temp_mean, wood_temp_mean
      flush(unit)
   end subroutine write_fast_probe

   !----- Grow-only capacity check for the per-patch forcing buffers. `forc` is reused across the    !
   !      patch loop, and the caller pre-sizes it to the site-wide max cohort count, so the loop     !
   !      never reallocates. fill_forcing and apply_rt_forcing write indices 1..col_cohort%n (never  !
   !      size(forc%...)), so reusing a larger patch's capacity for a smaller one is exact. ---------!
   subroutine alloc_forcing(forc, ncoh)
      type(column_forcing_t), intent(inout) :: forc
      integer(ik),            intent(in)    :: ncoh
      if (allocated(forc%abs_sw)) then
         if (size(forc%abs_sw) < ncoh) deallocate(forc%abs_sw, forc%abs_lw, forc%abs_par,        &
                                                  forc%abs_sw_wood, forc%abs_lw_wood, forc%par_top)
      end if
      if (.not. allocated(forc%abs_sw)) allocate(forc%abs_sw(ncoh), forc%abs_lw(ncoh),           &
                                                 forc%abs_par(ncoh), forc%abs_sw_wood(ncoh), forc%abs_lw_wood(ncoh), &
                                                 forc%par_top(ncoh))
   end subroutine alloc_forcing

   !----- Fill the per-patch prescribed forcing from the sub-step's forcing record and the ground  !
   !      shortwave (the caller's: f_ground x canopy-top with forcing, the reference value without). !
   subroutine fill_forcing(forc, col_cohort, met, sw_ground, sum_lai)
      type(column_forcing_t), intent(inout) :: forc
      type(column_cohort_t),  intent(in)    :: col_cohort
      type(met_forcing_t),    intent(in)    :: met
      real(wp),               intent(in)    :: sw_ground, sum_lai
      integer(ik) :: j
      forc%enthalpy_atm  = cas_enthalpy_of_temp(met%tair_k, met%qair)
      forc%shv_atm       = met%qair
      forc%co2_atm       = met%co2
      forc%abs_sw_ground = sw_ground
      forc%abs_lw_ground = 0.0_wp
      !----- The ONE rain/snow split: the record's total, at the record's own temperature -- with a    !
      !      forcing source, the patch's canopy-air top (met_to_cas_top), where the precipitation's      !
      !      enthalpy is valued too (precip_enthalpy, forc%air_temp below). -----------------------------!
      call precip_phase(met%precip, met%tair_k, forc%rainfall, forc%snowfall)
      forc%air_temp          = met%tair_k                ! rainfall enthalpy reference (snow/rain-on-snow)
      forc%par_per_w     = 2.1_wp                    ! LAI-split path: total-SW->PAR blend (abs_par == abs_sw)
      !----- NO-FORCING FALLBACK: split the canopy-top shortwave across cohorts by LAI share.   !
      !      When forcing is on, apply_rt_forcing overrides this with the real two-stream solve.  !
      do j = 1_ik, col_cohort%n
         if (sum_lai > tiny_num) then
            forc%abs_sw(j) = met%swdown() * col_cohort%lai(j) / sum_lai
         else
            forc%abs_sw(j) = 0.0_wp
         end if
         forc%abs_par(j) = forc%abs_sw(j)            ! no PAR/NIR split in the LAI path -> PAR==SW (biased high)
         forc%par_top(j) = met%par_beam + met%par_diffuse   ! no light gradient in the LAI path
         forc%abs_lw(j)  = 0.0_wp
         forc%abs_sw_wood(j) = 0.0_wp                ! no-forcing path: no wood absorption
         forc%abs_lw_wood(j) = 0.0_wp
      end do
   end subroutine fill_forcing

   !----- The context's constant reference climate as a forcing record, for a run without a forcing  !
   !      source. Only the TOTAL shortwave matters there -- fill_forcing's LAI-share split; the       !
   !      two-stream RT runs only with forcing -- so it sits in one stream and the total is exact.    !
   pure function reference_met(ctx) result(met)
      type(fast_context_t), intent(in) :: ctx
      type(met_forcing_t) :: met
      met%tair_k   = ctx%air_temp ; met%qair = ctx%shv_atm ; met%psurf_pa = ctx%press
      met%rho_air  = ctx%rho_air  ; met%co2  = ctx%co2_atm ; met%wind     = ctx%u_ref
      met%precip   = ctx%precip
      met%par_beam = ctx%rad_sw_top
      met%par_diffuse = 0.0_wp ; met%nir_beam = 0.0_wp ; met%nir_diffuse = 0.0_wp
   end function reference_met

   !----- RT join (docs/science/canopy_radiation_transfer.md §4): run the two-stream canopy         !
   !      radiation for this patch and OVERWRITE the LAI-share SW split in `forc` with real         !
   !      per-cohort absorbed SW (leaf energy) + PAR (photosynthesis) + NET longwave (leaf +        !
   !      ground) + below-canopy ground SW. Cohorts are gathered height-DESCENDING (top=1) but the  !
   !      two-stream wants BOTTOM(1)->TOP, so we build an ascending-height permutation `perm` and   !
   !      inverse-scatter the outputs. PAR: the leaf model's `env%par` is INCIDENT PAR (it          !
   !      re-applies cfg%leaf_absorptance internally for electron transport), so we divide the      !
   !      two-stream ABSORBED VIS by that same absorptance to hand back an incident-equivalent PAR  !
   !      -- otherwise leaf absorptance would be applied twice. abs_sw stays true ABSORBED SW (the  !
   !      leaf energy balance wants absorbed).                                                      !
   subroutine apply_rt_forcing(forc, ncoh, pft, lai, wai, height, leaf_temp, wood_temp,                &
                               soil_temp_top, snow, snow_params, soil_albedo, soil_emiss, rad_opt, met, &
                               leaf_absorptance)
      type(column_forcing_t),  intent(inout) :: forc
      integer(ik),             intent(in)    :: ncoh
      integer(ik),             intent(in)    :: pft(:)           !< per-cohort PFT (gather order, top first)
      real(wp),                intent(in)    :: lai(:), wai(:)   !< [m2/m2] leaf / wood area index
      real(wp),                intent(in)    :: height(:)        !< [m] cohort height
      real(wp),                intent(in)    :: leaf_temp(:), wood_temp(:)   !< [K] lagged tissue temperatures
      real(wp),                intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground emission)
      type(snow_column_t),     intent(in)    :: snow             !< pack store (albedo / emissivity ramp)
      type(snow_params_t),     intent(in)    :: snow_params
      real(wp),                intent(in)    :: soil_albedo(:)   !< [-] bare-soil albedo per band
      real(wp),                intent(in)    :: soil_emiss       !< [-] bare-soil emissivity
      type(rad_pft_optics_t),  intent(in)    :: rad_opt          !< per-PFT canopy optics (two-stream)
      type(met_forcing_t),     intent(in)    :: met
      real(wp),                intent(in)    :: leaf_absorptance !< [-] leaf PAR absorptance (incident-PAR conversion)
      integer(ik) :: j, i, ig
      integer(ik) :: perm(ncoh), pft_bt(ncoh)
      real(wp)    :: lai_bt(ncoh), wai_bt(ncoh), tcan_bt(ncoh), hgt_bt(ncoh)
      real(wp)    :: lf_bt
      type(rad_forcing_t)   :: rf
      type(rad_flux_t)      :: flux
      type(ground_optics_state_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)
      real(wp) :: snow_fl, snow_fc

      !----- A bare patch (ncoh == 0) is NOT special-cased: the zero-trip perm/scatter loops fall     !
      !      through and canopy_radiation's own empty-canopy branch returns the correct NET ground SW  !
      !      (incident * (1 - soil albedo)), so a patch shedding its last cohort stays continuous.     !

      !----- perm: gather indices in ASCENDING height (bottom -> top). ------------------------------!
      call ascending_order(height, ncoh, perm)
      do j = 1_ik, ncoh
         i = perm(j)
         pft_bt(j) = pft(i) ; lai_bt(j) = lai(i) ; hgt_bt(j) = height(i)
         wai_bt(j) = wai(i)
         !----- LW emission temperature: the cohort's AREA-WEIGHTED effective radiative temperature     !
         !      T_rad = (f_l*T_leaf^4 + (1-f_l)*T_wood^4)^(1/4), f_l = LAI/(LAI+WAI), from the lagged   !
         !      (start-of-sub-step) tissue temperatures (docs/science/canopy_radiation_transfer.md §3), !
         !      so the longwave field between cohorts, sky and ground sees both tissues. abs_lw below   !
         !      is therefore the NET longwave with the tissue emitting at T_rad, while                  !
         !      veg_energy_balance linearizes the emission change about T_cas                           !
         !      (docs/science/vegetation_energy_dynamics.md eq. (1)): emission is not counted exactly   !
         !      once at the tissue's own temperature, a known residual ("The emission base" there). ----!
         lf_bt      = lai(i) / max(lai(i) + wai(i), tiny_num)
         tcan_bt(j) = (lf_bt * leaf_temp(i) ** 4                                             &
                       + (1.0_wp - lf_bt) * wood_temp(i) ** 4) ** 0.25_wp
      end do

      !----- rad_forcing_t from met (all W/m2, direct assignment). --------------------------------!
      call alloc_rad_forcing(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = met%cosz
      rf%incid_beam(RAD_VIS) = met%par_beam ; rf%incid_diff(RAD_VIS) = met%par_diffuse
      rf%incid_beam(RAD_NIR) = met%nir_beam ; rf%incid_diff(RAD_NIR) = met%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp       ; rf%incid_diff(RAD_LW)  = met%lwdown

      !----- Ground optics: bare-soil albedo/emissivity from [soil], skin temperature from the top soil node. --!
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = soil_albedo ; surf%soil_emiss = soil_emiss
      surf%soil_temp   = soil_temp_top
      !----- Snow raises the ground albedo/emissivity + emits off the snow surface                     !
      !      (docs/science/canopy_radiation_transfer.md §4), RAMPED by the Niu-Yang07 snow-cover       !
      !      fraction so a partial pack gives a partial (continuous) albedo -- no threshold cliff.     !
      !      VIS/NIR fresh<->aged interpolated by the lagged surface liquid fraction. -----------------!
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

      call canopy_radiation(rad_opt, rf, ncoh, pft_bt, hgt_bt, lai_bt, wai_bt, tcan_bt, flux)

      !----- inverse-scatter: RT index j (bottom->top) maps to gather index perm(j). --------------!
      do j = 1_ik, ncoh
         ig = perm(j)
         forc%abs_sw(ig)  = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! total ABSORBED leaf SW (energy)
         forc%abs_par(ig) = flux%abs_leaf(RAD_VIS, j) / max(leaf_absorptance, tiny_num)  ! -> INCIDENT-equiv PAR
         forc%abs_lw(ig)  = flux%abs_leaf(RAD_LW, j)                                ! NET leaf LW (emission at T_rad incl.)
         forc%abs_sw_wood(ig) = flux%abs_wood(RAD_VIS, j) + flux%abs_wood(RAD_NIR, j)  ! ABSORBED wood SW (WAI share)
         forc%abs_lw_wood(ig) = flux%abs_wood(RAD_LW, j)                               ! NET wood LW
         forc%par_top(ig)     = flux%incid_top(RAD_VIS, j)                             ! PAR at the cohort's top
      end do
      forc%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                         + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forc%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)          ! NET ground LW (soil emission incl.)
      forc%par_per_w     = par_w_2_umol                        ! true VIS absorbed -> photon flux
      !----- TOP-OF-CANOPY fluxes (docs/science/canopy_radiation_transfer.md, "Surface radiative    !
      !      fluxes on the diagnostic path"). canopy_radiation already forms albedo(b) =            !
      !      utop/incid, so the upwelling is that ratio times the incident it was formed from -- no !
      !      second solve. Kept as FLUXES rather than as a time-averaged albedo: a period-mean      !
      !      albedo is the mean of a RATIO, and at night the SW ratio is 0/0; a reader forms the    !
      !      albedo as sum(up)/sum(down) over the period, eq. (15), which is what a satellite       !
      !      product is. ---------------------------------------------------------------------------!
      forc%sw_in_vis = rf%incid_beam(RAD_VIS) + rf%incid_diff(RAD_VIS)
      forc%sw_in_nir = rf%incid_beam(RAD_NIR) + rf%incid_diff(RAD_NIR)
      forc%sw_up_vis = flux%albedo(RAD_VIS) * forc%sw_in_vis
      forc%sw_up_nir = flux%albedo(RAD_NIR) * forc%sw_in_nir
      forc%lw_up     = flux%albedo(RAD_LW)  * (rf%incid_beam(RAD_LW) + rf%incid_diff(RAD_LW))
   end subroutine apply_rt_forcing

   !----- Fill the aerodynamics env from the sub-step's forcing + the patch's current CAS/ground. -!
   subroutine fill_aenv(aenv, biophys, met, zref)
      type(aero_env_t),      intent(inout) :: aenv
      type(patch_biophys_t), intent(in)    :: biophys
      type(met_forcing_t),   intent(in)    :: met
      real(wp),              intent(in)    :: zref     !< [m] the reference height the forcing is applied at
      aenv%u_ref = met%wind ; aenv%zref = zref ; aenv%press = met%psurf_pa ; aenv%rho_air = met%rho_air
      !----- The atmospheric reference state and the CAS/ground refresh go through               !
      !      set_aero_env_atm and set_aero_env_canopy (meds_canopy_types), the routines the      !
      !      column tests use too, so no caller assembles `aenv` by a parallel hand-written copy !
      !      that could leave `theta_atm` at its 298.15 K default. ------------------------------!
      call set_aero_env_atm(aenv, met%tair_k, met%qair, met%co2)
      call set_aero_env_canopy(aenv, biophys%cas%can_temp, biophys%cas%can_shv, biophys%cas%can_co2,           &
                               biophys%soil_energy%soil_temp(1))
   end subroutine fill_aenv

   !=======================================================================================!
   !  The patch block's row for ONE (patch, sub-step): the SAME numbers the physics just   !
   !  used, nothing recomputed. The slow rows stay 0 (the slow driver writes them). Every  !
   !  entry is a rate or a state, so the patch diagnostics can dt-weight the row whole     !
   !  (accumulate_patch_diag) and the FAST tier can stage it as it is.                     !
   !=======================================================================================!
   pure subroutine patch_diag_row(row, dt, le_flux, h_flux, rnet, sw_ground, lw_ground,             &
                                  ustar, ggnet, rough, displace, cas_temp, cas_shv, cas_co2, gpp, nee, &
                                  npp, soil_temp_top, resid_energy, resid_water, wind_top, tair_top,  &
                                  z_top, sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up, w_surface)
      real(wp), intent(out) :: row(N_PDIAG)
      real(wp), intent(in)  :: dt                        !< [s]        the sub-step
      real(wp), intent(in)  :: le_flux, h_flux           !< [W/m2]     CAS -> atmosphere
      real(wp), intent(in)  :: rnet                      !< [W/m2]     net all-wave radiation absorbed
      real(wp), intent(in)  :: sw_ground, lw_ground      !< [W/m2]     ground SW / net LW
      real(wp), intent(in)  :: ustar, ggnet              !< [m/s]      friction velocity, ground conductance
      real(wp), intent(in)  :: rough, displace           !< [m]        roughness, displacement height
      real(wp), intent(in)  :: cas_temp, cas_shv, cas_co2 !< [K],[kg/kg],[umol/mol] canopy air
      real(wp), intent(in)  :: gpp, nee, npp             !< [umol/m2/s] gross uptake, net exchange (+ to atm), NPP
      real(wp), intent(in)  :: soil_temp_top             !< [K]        top soil-node temperature
      real(wp), intent(in)  :: resid_energy, resid_water !< [J/m2],[kg/m2] this step's SIGNED ledger residuals
      real(wp), intent(in)  :: wind_top, tair_top        !< [m/s],[K] the forcing at the canopy-air top
      real(wp), intent(in)  :: z_top                     !< [m]   that top (the canopy-air depth)
      !----- Top-of-canopy radiative fluxes per band. -----------------------------------------!
      real(wp), intent(in)  :: sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up
      real(wp), intent(in)  :: w_surface                 !< [kg/m2]   ponded surface water
      row = 0.0_wp
      row(PD_LE)           = le_flux
      row(PD_H)            = h_flux
      row(PD_RNET)         = rnet
      row(PD_SW_GROUND)    = sw_ground
      row(PD_LW_GROUND)    = lw_ground
      row(PD_SW_IN_VIS)    = sw_in_vis
      row(PD_SW_IN_NIR)    = sw_in_nir
      row(PD_SW_UP_VIS)    = sw_up_vis
      row(PD_SW_UP_NIR)    = sw_up_nir
      row(PD_SW_UP)        = sw_up_vis + sw_up_nir
      row(PD_LW_UP)        = lw_up
      row(PD_SKIN_TEMP)    = sqrt(sqrt(max(lw_up, 0.0_wp) / stefan))
      row(PD_USTAR)        = ustar
      row(PD_GGNET)        = ggnet
      row(PD_ROUGH)        = rough
      row(PD_DISPLACE)     = displace
      row(PD_CAS_TEMP)     = cas_temp
      row(PD_CAS_SHV)      = cas_shv
      row(PD_CAS_CO2)      = cas_co2
      !----- VPD is NONLINEAR in temperature, so the dt-weighted mean of VPD is not the VPD of    !
      !      the mean twins -- it has to be formed here, per sub-step, not derived at read        !
      !      time. -------------------------------------------------------------------------------!
      row(PD_CAS_VPD)      = specific_humidity_to_vpd(cas_temp, cas_shv, p_std)
      row(PD_W_SURFACE)    = w_surface
      !----- The patch's own forcing: the sample moved to its canopy-air top                        !
      !      (docs/science/forcing.md §8), the only forcing values that differ between patches. ----!
      row(PD_WIND_CAS_TOP) = wind_top
      row(PD_TAIR_CAS_TOP) = tair_top
      row(PD_Z_CAS_TOP)    = z_top
      !----- CARBON. nee is the model's own NEE, sign-positive to the atmosphere -- the same number  !
      !      the CAS CO2 box is driven by, so the flux and the state it acts on cannot disagree.      !
      !      Reco follows from the NEE identity, Reco = NEE + GPP, rather than being summed a second  !
      !      way. -----------------------------------------------------------------------------------!
      row(PD_GPP)          = gpp
      row(PD_NEE)          = nee
      row(PD_NPP)          = npp
      row(PD_RECO)         = nee + gpp
      !----- Transpiration as a WATER flux [kg/m2/s]: the latent flux is the canopy-air -> atmosphere  !
      !      total, so this is the evaporative flux the CAS actually shed, not a stomatal-only term.    !
      !      Each cohort's own transpiration is available per cohort (CD_TRANSP) for the split.         !
      row(PD_TRANSP)       = le_flux / latent_heat_vap
      row(PD_SOIL_TEMP_TOP) = soil_temp_top
      !----- Whole-column budget residuals. These are the numbers that decide whether anything above  !
      !      this line can be believed, which is why they are captured on the same tick rather than    !
      !      left to an assertion nobody reads. The ledger gives this step's SIGNED imbalance [J/m2,    !
      !      kg/m2]; over the step it is the leak RATE [W/m2, kg/m2/s], sign-positive when store       !
      !      appears from nowhere, and the time mean of that rate is what the registry reports.        !
      row(PD_RESID_ENERGY) = resid_energy / dt
      row(PD_RESID_WATER)  = resid_water / dt
   end subroutine patch_diag_row

   !----- Fold ONE (patch, sub-step) row into the per-patch diagnostic accumulator, dt-weighted and  !
   !      normalized on read, so the result is a correct time-mean even though the adaptive          !
   !      controller may have taken a different number of inner steps in different sub-steps. -----!
   subroutine accumulate_patch_diag(pd, ip, dt, row)
      type(patch_diag_block), intent(inout) :: pd
      integer(ik),            intent(in)    :: ip
      real(wp),               intent(in)    :: dt        !< [s] sample weight
      real(wp),               intent(in)    :: row(N_PDIAG)
      pd%v(:, ip) = pd%v(:, ip) + row * dt
      pd%w(ip)    = pd%w(ip)    + dt
   end subroutine accumulate_patch_diag

   !----- Size the FAST tier's staging for this slow step. The sub-step count is run-constant; the  !
   !      patch and cohort counts are this step's, and the buffers only grow (to at least twice    !
   !      their size, so a growing stand reallocates rarely). The output layer checks the live      !
   !      counts against output.patch_max and cohort_max, as it does for the coarse tiers. --------!
   subroutine size_fast_staging(b, nsub, npatch, ncohort, with_cdiag)
      type(output_buffers_t), intent(inout) :: b
      integer(ik),            intent(in)    :: nsub, npatch, ncohort
      logical,                intent(in)    :: with_cdiag   !< the per-cohort diagnostic block is on
      integer(ik) :: nl, np, nc
      nl = n_soil_layer_max
      if (.not. allocated(b%fast_time))                                                            &
         allocate(b%fast_time(nsub), b%fast_forcing(N_PYDIAG, nsub), b%fast_site(N_PDIAG, nsub),    &
                  b%fast_soil_temp(nl, nsub), b%fast_soil_water(nl, nsub))
      np = max(npatch, 1_ik)
      if (allocated(b%fast_patch)) then
         if (size(b%fast_patch, 2, kind=ik) < np) then
            np = max(np, 2_ik * size(b%fast_patch, 2, kind=ik))
            deallocate(b%fast_patch, b%fast_soil_temp_patch, b%fast_soil_water_patch)
         end if
      end if
      if (.not. allocated(b%fast_patch))                                                           &
         allocate(b%fast_patch(N_PDIAG, np, nsub), b%fast_soil_temp_patch(nl * np, nsub),           &
                  b%fast_soil_water_patch(nl * np, nsub))
      nc = max(ncohort, 1_ik)
      if (allocated(b%fast_coh_ltemp)) then
         if (size(b%fast_coh_ltemp, 1, kind=ik) < nc) then
            nc = max(nc, 2_ik * size(b%fast_coh_ltemp, 1, kind=ik))
            deallocate(b%fast_coh_ltemp, b%fast_coh_gpp, b%fast_coh_height)
         end if
      end if
      if (.not. allocated(b%fast_coh_ltemp))                                                       &
         allocate(b%fast_coh_ltemp(nc, nsub), b%fast_coh_gpp(nc, nsub), b%fast_coh_height(nc, nsub))
      !----- The diagnostic block's rows, staged only while the block is on: its FAST-tier variables  !
      !      (SRC_F_CD0 + CD_*) are what turn it on, so a run that reports none allocates nothing. ----!
      if (with_cdiag) then
         nc = size(b%fast_coh_ltemp, 1, kind=ik)
         if (allocated(b%fast_coh_cdiag)) then
            if (size(b%fast_coh_cdiag, 1, kind=ik) < nc) deallocate(b%fast_coh_cdiag)
         end if
         if (.not. allocated(b%fast_coh_cdiag)) allocate(b%fast_coh_cdiag(nc, N_CDIAG, nsub))
      end if
   end subroutine size_fast_staging

   !----- The polygon's forcing for one sub-step (MEDS_FORCING_DESIGN.md §6.7), as the fast loop  !
   !      used it: after the reader's shortwave partition, phase split and optional corrections. --!
   subroutine accumulate_polygon_diag(d, dt, f)
      type(polygon_diag_block), intent(inout) :: d
      real(wp),                 intent(in)    :: dt     !< [s] sample weight
      real(wp),                 intent(in)    :: f(N_PYDIAG)   !< the sub-step's forcing echo
      d%v = d%v + f * dt
      d%w = d%w + dt
   end subroutine accumulate_polygon_diag

   !----- The forcing one sub-step used, in the polygon block's table (PY_*): the one list of the  !
   !      forcing echo, which the coarse tiers accumulate (accumulate_polygon_diag) and the FAST    !
   !      tier stages as it is. -------------------------------------------------------------------!
   pure subroutine forcing_echo(met, rainf, snowfall, f)
      type(met_forcing_t), intent(in)  :: met
      real(wp),            intent(in)  :: rainf, snowfall   !< [kg/m2/s] the split the patches used, area-weighted
      real(wp),            intent(out) :: f(N_PYDIAG)
      f(PY_SW_IN)       = met%swdown()
      f(PY_PRECIP)      = met%precip
      f(PY_TAIR)        = met%tair_k
      f(PY_QAIR)        = met%qair
      f(PY_PSURF)       = met%psurf_pa
      f(PY_WIND)        = met%wind
      f(PY_LWDOWN)      = met%lwdown
      f(PY_PAR_BEAM)    = met%par_beam
      f(PY_PAR_DIFFUSE) = met%par_diffuse
      f(PY_NIR_BEAM)    = met%nir_beam
      f(PY_NIR_DIFFUSE) = met%nir_diffuse
      f(PY_SNOWFALL)    = snowfall
      f(PY_CO2)         = met%co2
      f(PY_COSZ)        = met%cosz
      f(PY_RHO_AIR)     = met%rho_air
      f(PY_RAINF)       = rainf
   end subroutine forcing_echo

end module meds_fast_dynamics
