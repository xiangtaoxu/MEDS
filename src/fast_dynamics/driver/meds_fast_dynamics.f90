! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_dynamics -- the polygon's fast loop (docs/science/order_of_processes.md §3).   !
! fast_dynamics gathers each patch's cohorts into the column view, adopts the patch's      !
! stores and the cohorts' tissue energy, water and films, runs n_fast_per_slow sub-steps   !
! of column_fast_step, and writes the evolved state back. It accumulates the GPP and       !
! maintenance respiration the slow step's carbon allocation reads, and stages the output.  !
! build_fast_context and init_fast_reservoirs prepare a polygon for it.                    !
!==========================================================================================!
module meds_fast_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, umol_2_kgc, latent_heat_vap, day_sec, p_std,              &
                                     par_w_2_umol, stefan_boltzmann, n_soil_layer_max
   use meds_numerics,         only : ascending_order
   use meds_config,           only : meds_config_t, top_soil_psi
   use meds_budget_check,     only : budget_t, budget_merge
   use meds_biogeochem_types, only : IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   use meds_therm_lib,           only : cas_temp_to_enthalpy, temp_to_internal_energy,                            &
                                        specific_humidity_to_vpd, precip_phase
   use meds_fast_config, only : fast_context_t, site_column_configs
   use meds_column_view, only : copy_column_cohort
   use meds_time,             only : meds_time_t, time_advance_seconds, time_to_string, seconds_into_day,         &
                                     slope_incidence_cosz
   use meds_output_types,     only : output_buffers_t
   use meds_diag_vars,        only : N_CO_FAST_ROWS, N_PA_ROWS, patch_diag_block, polygon_diag_block, PA_LE, PA_H, &
                                     PA_RNET, PA_SW_GROUND, &
                                     PA_LW_GROUND, PA_SW_IN_VIS, PA_SW_IN_NIR, PA_SW_UP_VIS, PA_SW_UP_NIR, PA_LW_UP, PA_USTAR, &
                                     PA_GGNET, PA_ROUGH, PA_DISPLACE, PA_CAS_TEMP, PA_CAS_SHV, PA_CAS_CO2, PA_GPP, PA_NEE, &
                                     PA_TRANSP, PA_NPP, PA_RECO, PA_SW_UP, PA_SKIN_TEMP, PA_SOIL_TEMP_TOP, PA_CAS_VPD, &
                                     PA_W_SURFACE, PA_RESID_ENERGY, PA_RESID_WATER, PA_WIND_CAS_TOP, PA_TAIR_CAS_TOP, &
                                     PA_Z_CAS_TOP, PY_SW_IN, PY_PRECIP, PY_TAIR, PY_QAIR, PY_PSURF, PY_WIND, PY_LWDOWN, &
                                     PY_PAR_BEAM, PY_PAR_DIFFUSE, PY_NIR_BEAM, PY_NIR_DIFFUSE, PY_SNOWFALL, PY_CO2, PY_COSZ, &
                                     PY_RHO_AIR, PY_RAINF, N_PY_ROWS, expand_cohort_diag, expand_patch_diag, &
                                     site_diag_block, SI_SW_IN, SI_PAR_BEAM, SI_NIR_BEAM, SI_LWDOWN, SI_TAIR,  &
                                     SI_QAIR, SI_PSURF, SI_RHO_AIR, N_SI_ROWS
   use meds_patch_state_vars,   only : xi_accum_t, snow_column_t
   use meds_met_types,    only : met_source_t, met_cursor_t, met_forcing_t
   use meds_radiation_lib,    only : slope_beam_factor
   use meds_met,       only : met_advance, met_instant, met_to_cas_top, met_at_site
   use meds_canopy_aerodynamics, only : canopy_roughness, open_canopy_fraction
   use meds_cohort_state_vars, only : set_cohort_tissue_water, hold_cohort_tissue_temps, cohort_tissue_temps, &
                                      set_cohort_tissue_temps
   use meds_state_defaults,    only : DMAX_PSI_LEAF_ACCUM_RESET
   use meds_polygon_state_vars, only : polygon_t, patch_polygon_area
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, ensure_aero_out_capacity,                   &
                                 rad_pft_optics_t, rad_forcing_t, rad_flux_t, alloc_rad_forcing,                  &
                                 N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW, set_aero_env_atm,                  &
                                 set_aero_env_canopy, ground_optics_state_t
   use meds_fast_types, only : patch_biophys_t, ensure_patch_biophys_capacity, fast_workspace_t,                  &
                               column_config_t, column_cohort_t, column_forcing_t, column_budget_t,               &
                               ensure_column_cohort_capacity, column_site_t
   use meds_biophysics_params, only : snow_params_t
   use meds_canopy_radiation, only : canopy_radiation, ground_optics
   use meds_ground_biophysics, only : snow_cover_fraction
   use meds_fast_step,       only : column_fast_step
   !$ use omp_lib,            only : omp_get_thread_num
   implicit none
   private

   public :: init_fast_reservoirs, fast_dynamics

   !----- Slots of the reduction staging: each (patch, sub-step) writes its own cells, folded in patch !
   !      order after the loop, so the sums do not move with the thread count                       !
   !      (docs/science/numerical_scheme.md §6a). ---------------------------------------------------!
   integer(ik), parameter :: RED_ET            =  1_ik   !< polygon ET [kg/m2], area-weighted
   integer(ik), parameter :: RED_PHENO_TAIR    =  2_ik   !< daily-mean air-temperature accumulator [K]
   integer(ik), parameter :: RED_INTEG_STEPS   =  3_ik   !< integrator work: accepted steps
   integer(ik), parameter :: RED_INTEG_REJ     =  4_ik   !< integrator work: rejected steps
   integer(ik), parameter :: RED_WATER_NSUB    =  5_ik   !< water network: accepted sub-steps
   integer(ik), parameter :: RED_WATER_NEWTON  =  6_ik   !< water network: Newton iterations
   integer(ik), parameter :: RED_WATER_REJECT  =  7_ik   !< water network: rejected sub-steps
   integer(ik), parameter :: RED_WATER_NFAIL   =  8_ik   !< water network: sub-steps taken unconverged
   integer(ik), parameter :: RED_RK45_RESCUE   =  9_ik   !< dt_fast steps RK45 handed to ARK (rescues)
   integer(ik), parameter :: RED_CLAMP_STAGE   = 10_ik   !< stage-level state clamps
   !----- Rain and snow as split at each patch's canopy-air top, area-weighted (docs/science/forcing.md §7). !
   integer(ik), parameter :: RED_RAINF         = 11_ik   !< rain [kg/m2/s]
   integer(ik), parameter :: RED_SNOWFALL      = 12_ik   !< snow [kg/m2/s]
   integer(ik), parameter :: RED_RUNOFF_SURF   = 13_ik   !< surface runoff [kg/m2], area-weighted
   integer(ik), parameter :: RED_DRAINAGE      = 14_ik   !< drainage at the column's bottom [kg/m2]
   integer(ik), parameter :: N_RED             = 14_ik
   !----- [K] round-off allowed between a cohort's stored tissue temperatures and those its energy gives. !
   real(wp),    parameter :: TISSUE_TEMP_TOL   = 1.0e-6_wp

contains

   !----- Seed every patch's fast stores to a horizontally uniform equilibrium, and its cohorts with  !
   !      them as a recruit is born: at the canopy-air temperature, holding the water their tissues  !
   !      hold at the top soil layer's potential. Called once, after the community is built. -------!
   subroutine init_fast_reservoirs(poly, ctx)
      type(polygon_t),      intent(inout) :: poly
      type(fast_context_t), intent(in)    :: ctx
      integer(ik) :: ip, k, nsl, i
      real(wp)    :: psi
      do ip = 1_ik, poly%patch%n
         associate (cas => poly%patch%cas(ip), se => poly%patch%soil_energy(ip), sw => poly%patch%soil_water(ip), &
                    surf => poly%site%surface(poly%patch%owner_site(ip)))
            nsl = surf%soil%n_active
            sw%theta(1:nsl)  = surf%theta_start(1:nsl)
            sw%w_surface = 0.0_wp
            sw%w_surface_energy = 0.0_wp        ! dry pond -> zero enthalpy
            do k = 1_ik, nsl
               se%soil_energy(k) = temp_to_internal_energy(surf%soil_thermal%soil_dry_heat_capacity(k),             &
                                   surf%theta_start(k) * rho_h2o, ctx%soil_temp_init, 1.0_wp)
               se%soil_temp(k)   = ctx%soil_temp_init
               se%soil_fliq(k) = 1.0_wp
            end do
            cas%can_shv      = ctx%shv_atm
            cas%can_co2 = ctx%co2_atm
            cas%can_enthalpy = cas_temp_to_enthalpy(ctx%air_temp, ctx%shv_atm)
            cas%can_temp     = ctx%air_temp
            psi = top_soil_psi(surf%soil, sw%theta(1))
            do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
               poly%cohort%leaf_temp(i) = cas%can_temp
               poly%cohort%wood_temp(i) = cas%can_temp
               associate (hyd => ctx%col_config%hydraulics_table%pft(poly%cohort%pft(i)))
                  call set_cohort_tissue_water(poly%cohort, i, psi, hyd%leaf_curve, hyd%wood_curve)
               end associate
               call hold_cohort_tissue_temps(poly%cohort, i)
            end do
         end associate
      end do
   end subroutine init_fast_reservoirs

   !=======================================================================================!
   !  Advance every patch of the polygon by one slow step of fast biophysics:              !
   !  n_fast_per_slow sub-steps on the forcing source's samples, or on the context's       !
   !  reference climate without one. The optional outputs report the worst whole-column    !
   !  budget residuals and the failure count, for a caller that asserts conservation.      !
   !=======================================================================================!
   subroutine fast_dynamics(poly, ctx, cfg, met_src, met_cur, step_start, worst_energy, worst_water, &
                            n_budget_fail, out_bufs, run_energy_budget, run_water_budget,             &
                            run_face_budget)
      type(polygon_t),      intent(inout) :: poly
      type(fast_context_t), intent(in)    :: ctx
      type(meds_config_t),  intent(in)    :: cfg
      type(met_source_t), optional, intent(in)    :: met_src       !< the run's forcing source (read-only here)
      type(met_cursor_t), optional, intent(inout) :: met_cur       !< this polygon's cursor (per-sub-step met)
      type(meds_time_t),  optional, intent(in)    :: step_start    !< calendar time at the start of this slow step
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail
      type(output_buffers_t), optional, intent(inout) :: out_bufs  !< this polygon's FAST-tier staging (when on)
      !----- Run-level whole-column ledgers, which take the step's area-weighted ledgers: their signed  !
      !      cumulative residual shows a one-signed bias below the per-step tolerance. ----------------!
      type(budget_t), optional, intent(inout) :: run_energy_budget, run_water_budget
      !----- The per-layer face-closure ledger, folded the same way: it sees a vertical misplacement. -!
      type(budget_t), optional, intent(inout) :: run_face_budget

      !----- The sub-step met samples and their times, made once per slow step outside the patch loop: !
      !      they are the same for every patch, and the file reader stays out of the parallel region. -!
      type(met_forcing_t), allocatable :: met_sample(:)
      type(meds_time_t),   allocatable :: t_sample(:)
      !----- Each site's forcing per sub-step: the polygon's sample after the site's terrain. -----!
      type(met_forcing_t), allocatable :: met_site(:,:)   !< (sub-step, site)
      !----- Reduction staging: one cell per (slot, sub-step, patch), assigned once in the patch loop  !
      !      (so never zeroed) and folded into the polygon in patch order after it. ----------------!
      real(wp),    allocatable :: red_site(:,:,:)                  !< (N_RED, sub-step, patch)
      real(wp),    allocatable :: red_worst_energy(:), red_worst_water(:)   !< (patch); max-folded
      type(budget_t), allocatable :: red_budget_energy(:), red_budget_water(:) !< (patch); area-merged
      type(budget_t), allocatable :: red_budget_face(:)                        !< (patch); area-merged
      type(budget_t) :: site_energy_budget, site_water_budget, site_face_budget
      integer(ik), allocatable :: red_nfail(:)                     !< (patch) budget-failure counts
      real(wp)    :: f_ground, sw_ground
      type(met_forcing_t) :: met_ref                  !< the context's reference climate (no forcing)
      integer(ik) :: ip, isub, npatch, ncoh_max, nsub, nl, n_thread, is, ks
      logical     :: do_forcing, do_fast, do_cdiag, do_pdiag
      !----- Per-thread scratch for the per-(cohort, sub-step) diagnostics, folded into poly%cohort%diag. !
      real(wp),    allocatable :: cdiag_pool(:,:,:)                !< (N_CO_FAST_ROWS, cohort, thread)
      !----- Each site's part of the column configuration (site_column_configs), read-only in the   !
      !      patch loop. ----------------------------------------------------------------------------!
      type(column_site_t), allocatable :: site_col(:)
      !=========================================================================================!
      !  The per-thread scratch pool: one copy of each per-patch working buffer per thread,     !
      !  indexed by thread number and aliased with `associate` in the loop. A pool, not OpenMP  !
      !  private copies or BLOCK locals, which ifx builds through shared compiler-made storage  !
      !  (CLAUDE.md, compiler traps; docs/building.md). Each entry is sized once, to the        !
      !  largest patch; every reader loops to the active count, never to size().                !
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
      real(wp)    :: prow(N_PA_ROWS)                !< the patch block's row for one (patch, sub-step)
      real(wp)    :: z_top, rough_p, displace_p   !< [m] this patch's canopy-air top, roughness, displacement
      integer(ik) :: j, i, i0, ncoh, ith, kb
      real(wp)    :: rain_s, snow_s, echo(N_PY_ROWS)   !< the forcing echo's precipitation split, per sub-step
      real(wp)    :: t_leaf_e, t_wood_e            !< [K] the tissue temperatures a cohort's energy gives

      !----- Live forcing drives the sub-steps when it is on and a source, cursor and step time are   !
      !      given; otherwise every sub-step reads the context's reference climate (met_ref). -------!
      do_forcing = cfg%forcing%forcing_on .and. present(met_src) .and. present(met_cur)          &
                   .and. present(step_start)
      f_ground = ctx%rad_sw_ground / max(ctx%rad_sw_top, tiny_num)   ! ground/canopy-top SW transmittance
      met_ref  = reference_met(ctx)
      npatch   = poly%patch%n
      nsub     = cfg%n_fast_per_slow
      n_thread = max(1_ik, ctx%patch_threads)

      !----- Reset the fast->slow accumulators the slow step reads after this window. -------------!
      poly%cohort%gpp_accum(1:poly%cohort%n)       = 0.0_wp
      poly%cohort%leaf_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%stem_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%root_resp_accum(1:poly%cohort%n) = 0.0_wp
      poly%cohort%light_hours_accum(1:poly%cohort%n) = 0.0_wp
      !----- Roll over the predawn water status: today's maximum leaf potential becomes the value the  !
      !      slow stomatal control uses tomorrow (docs/science/leaf_gas_exchange.md §4.2). A cohort   !
      !      with no sample today keeps its value; psi_leaf <= 0, so the guard tests the reset value. -!
      where (poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n) > 0.5_wp * DMAX_PSI_LEAF_ACCUM_RESET)          &
         poly%cohort%dmax_psi_leaf(1:poly%cohort%n) = poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n)
      poly%cohort%dmax_psi_leaf_accum(1:poly%cohort%n) = DMAX_PSI_LEAF_ACCUM_RESET
      !----- Reset the daily soil-carbon accumulator ([soil_carbon].soil_carbon_on). --------------!
      if (cfg%soil_carbon%soil_carbon_on) poly%patch%xi_accum(1:poly%patch%n) = xi_accum_t()
      !----- The phenology air-temperature sums, ET and the solver work counters are reset and     !
      !      refilled by the patch-order fold after the patch loop. A site's patches share its air   !
      !      temperature, so the mean over (patch, sub-step) samples is the daily mean. -----------!

      !----- Each site's column configuration, for the patch loop. ---------------------------------!
      allocate(site_col(poly%site%n))
      call site_column_configs(ctx, poly, site_col)

      !----- FAST-tier staging, filled when the F tier or a snapshot of a fast-loop value reads it. --!
      do_fast = present(out_bufs)
      if (do_fast) do_fast = out_bufs%fast_on
      !----- The per-cohort and per-patch diagnostic capture runs only when the registry asks for it. -!
      do_cdiag = poly%cohort%diag%active
      do_pdiag = poly%patch%diag%active
      !----- Size the diagnostic accumulators to the stand. Not zeroed: they run from one output tick  !
      !      to the next, and reset_step_diagnostics zeroes them (docs/science/order_of_processes.md §6). !
      if (do_cdiag) then
         call expand_cohort_diag(poly%cohort%diag, max(poly%cohort%n, 1_ik))
         poly%cohort%diag%n = poly%cohort%n
      end if
      if (do_pdiag) then
         call expand_patch_diag(poly%patch%diag, max(npatch, 1_ik))
         poly%patch%diag%n = npatch
      end if
      nl = n_soil_layer_max
      if (do_fast) then
         call size_fast_staging(out_bufs, nsub, npatch, poly%site%n, poly%cohort%n, do_cdiag)
         out_bufs%n_fast_sub    = nsub
         !----- Only the active layers are data: a site's column and each of its patches' have the  !
         !      site's layers, and the rest of the slab reads back as fill. ----------------------!
         out_bufs%fast_n_site = poly%site%n
         out_bufs%fast_soil_live = .false.
         out_bufs%fast_soil_patch_live = .false.
         do is = 1_ik, poly%site%n
            ks = (is - 1_ik) * nl
            out_bufs%fast_soil_live(ks + 1_ik:ks + poly%site%surface(is)%soil%n_active) = .true.
         end do
         do ip = 1_ik, npatch
            ks = (ip - 1_ik) * nl
            out_bufs%fast_soil_patch_live(ks + 1_ik:ks + poly%site%surface(poly%patch%owner_site(ip))%soil%n_active) &
               = .true.
         end do
         out_bufs%fast_n_cohort = poly%cohort%n
         out_bufs%fast_n_patch  = npatch
         out_bufs%fast_site = 0.0_wp
         out_bufs%fast_soil_temp = 0.0_wp
         out_bufs%fast_soil_water = 0.0_wp
         out_bufs%fast_cohort_gpp = 0.0_wp
         out_bufs%fast_cohort_height = 0.0_wp
         if (do_cdiag) out_bufs%fast_cohort_diag = 0.0_wp
      end if

      ncoh_max = 0_ik
      do ip = 1_ik, npatch
         ncoh_max = max(ncoh_max, poly%patch%cohort_count(ip))
      end do

      !----- Per-thread diagnostic scratch, sized from ncoh_max just above. -----------------------!
      allocate(cdiag_pool(N_CO_FAST_ROWS, max(ncoh_max,1_ik), max(n_thread,1_ik)))

      !----- The per-(sub-step, patch) reduction staging; every cell is assigned, so none is zeroed. -!
      allocate(red_site(N_RED, nsub, max(npatch,1_ik)))
      allocate(red_budget_energy(max(npatch,1_ik)), red_budget_water(max(npatch,1_ik)))
      allocate(red_budget_face(max(npatch,1_ik)))
      allocate(red_worst_energy(max(npatch,1_ik)), red_worst_water(max(npatch,1_ik)),               &
               red_nfail(max(npatch,1_ik)))

      !----- This slow step's met samples: met_advance once per sub-step in increasing time, and     !
      !      met_instant, a pure function of (reader state, time). The probe writes the sample times. !
      if (do_forcing) then
         allocate(met_sample(nsub), t_sample(nsub))
         met_cur%kt_last_day = poly%kt_last_day          ! the polygon's memory of the sky, carried over
         do isub = 1_ik, nsub
            t_sample(isub) = time_advance_seconds(step_start,                                          &
                       (real(isub, wp) - 1.0_wp + cfg%forcing_sample_frac) * cfg%dt_fast)
            call met_advance(met_src, met_cur, t_sample(isub))
            met_sample(isub) = met_instant(met_src, met_cur, t_sample(isub))
         end do
         poly%kt_last_day = met_cur%kt_last_day
         !----- Each site's sample, moved to its elevation and its beam to its slope (met_at_site,   !
         !      docs/science/forcing.md §8); a flat site at the polygon's elevation takes it as is. -!
         allocate(met_site(nsub, poly%site%n))
         do is = 1_ik, poly%site%n
            associate (surf => poly%site%surface(is))
               do isub = 1_ik, nsub
                  met_site(isub, is) = met_sample(isub)
                  if (surf%elevation_m == poly%elevation_m .and. surf%slope_deg == 0.0_wp) cycle
                  associate (m => met_sample(isub))
                     met_site(isub, is) = met_at_site(m, cfg%forcing, t_sample(isub)%month,             &
                        surf%elevation_m - poly%elevation_m,                                          &
                        slope_beam_factor(m%par_beam + m%nir_beam, m%cosz,                               &
                                          slope_incidence_cosz(t_sample(isub), seconds_into_day(t_sample(isub)), &
                                                               met_cur%latitude_deg, met_cur%longitude_deg,   &
                                                               surf%slope_deg, surf%aspect_deg)))
                  end associate
               end do
            end associate
         end do
      end if
      !----- A FAST record is stamped by its sub-step's start, as the slower tiers stamp theirs. ---!
      if (do_fast) then
         do isub = 1_ik, nsub
            out_bufs%fast_time(isub) = time_advance_seconds(step_start, real(isub - 1_ik, wp) * cfg%dt_fast)
         end do
      end if

      !=========================================================================================!
      !  The parallel patch loop (docs/science/numerical_scheme.md §6a): patches are            !
      !  independent within a dt_fast. schedule(dynamic, 1), because the adaptive march makes   !
      !  their costs very unequal. Without OpenMP the directives are comments and the pool is   !
      !  one deep.                                                                              !
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
      !$omp            sw_ground, z_top, rough_p, displace_p, t_leaf_e, t_wood_e)
      do ip = 1_ik, npatch
         !----- This thread's slot in the pool; without OpenMP the `!$` line is a comment: slot 1. --!
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
         ncoh = poly%patch%cohort_count(ip)
         i0   = poly%patch%cohort_offset(ip)

         !----- Gather the patch's cohorts into the column view (capacity ensured above). ----------!
         call copy_column_cohort(col_cohort, poly%cohort, i0, ncoh)
         sum_lai = 0.0_wp
         sum_wai = 0.0_wp
         do j = 1_ik, ncoh
            sum_lai = sum_lai + col_cohort%lai(j)
            sum_wai = sum_wai + col_cohort%wai(j)
         end do

         !----- Canopy geometry. The open fraction exp(-(LAI + WAI)) weighs bare ground in the ground !
         !      conductance (docs/science/canopy_aerodynamics.md eq. (6)), so a patch without cohorts !
         !      couples as bare soil. The canopy is never buried in snow (snowfac = 0). ------------!
         ageom%veg_height   = ctx%veg_height_bare
         do j = 1_ik, ncoh
            ageom%veg_height = max(ageom%veg_height, col_cohort%height(j))
         end do
         ageom%opencan_frac = open_canopy_fraction(sum_lai, sum_wai)
         ageom%snowfac      = 0.0_wp

         call alloc_forcing(forc, ncoh)

         !----- Adopt the patch's stores and each cohort's tissue state into the working bundle;    !
         !      every field is assigned here, so the bundle needs no seeding. --------------------!
         biophys%cas    = poly%patch%cas(ip)
         biophys%soil_energy = poly%patch%soil_energy(ip)
         biophys%soil_water = poly%patch%soil_water(ip)
         biophys%snow   = poly%patch%snow(ip)
         biophys%adapt_dt_last = poly%patch%adapt_dt_last(ip)   ! per-patch, not loop-carried
         !----- The soil-carbon pools, frozen for the day and never written back: soil_carbon_step is !
         !      their only writer (they stay 0 with soil carbon off). ----------------------------!
         biophys%soil_carbon = poly%patch%soil_carbon(ip)
         !----- The day's turnover shed water and slow-step CO2 rate, frozen the same way. ---------!
         biophys%shed_water_rate = poly%patch%shed_water_rate(ip)
         biophys%shed_water_enth = poly%patch%shed_water_enth(ip)
         biophys%slow_co2_rate   = poly%patch%slow_co2_rate(ip)
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            !----- Tissue energy is the state the fast loop advances (per m2 of ground here). A slow  !
            !      operator that changed a tissue's mass without hold_cohort_tissue_temps leaves the  !
            !      stored temperatures off the ones the energy gives, and stops the run here. -------!
            call cohort_tissue_temps(poly%cohort, i, t_leaf_e, t_wood_e)
            if (abs(t_leaf_e - poly%cohort%leaf_temp(i)) > TISSUE_TEMP_TOL .or.                          &
                abs(t_wood_e - poly%cohort%wood_temp(i)) > TISSUE_TEMP_TOL) then
               write(*, '(a,i0,4(a,es12.5))') ' fast_dynamics: cohort ', i, ' leaf ', poly%cohort%leaf_temp(i),    &
                     ' K stored vs ', t_leaf_e, ' K from its energy; wood ', poly%cohort%wood_temp(i), ' vs ', t_wood_e
               error stop 'fast_dynamics: a cohort''s tissue temperatures are not the ones its energy gives'
            end if
            biophys%leaf_energy(j) = poly%cohort%leaf_energy(i) * poly%cohort%nplant(i)
            biophys%wood_energy(j) = poly%cohort%wood_energy(i) * poly%cohort%nplant(i)
            biophys%leaf_temp(j) = poly%cohort%leaf_temp(i)
            biophys%wood_temp(j) = poly%cohort%wood_temp(i)
            !----- Tissue water is adopted as it is: fitting it to the capacity the biomass allows is    !
            !      the slow step's (docs/science/soil_plant_water.md §2). -------------------------------!
            biophys%leaf_water_mass(j) = poly%cohort%leaf_water_mass(i)
            biophys%wood_water_mass(j) = poly%cohort%wood_water_mass(i)
            !----- The canopy films, per plant on the cohort and per m2 of ground here (0 is dry). ---!
            biophys%leaf_surf_water(j) = poly%cohort%leaf_surf_water(i) * poly%cohort%nplant(i)
            biophys%wood_surf_water(j) = poly%cohort%wood_surf_water(i) * poly%cohort%nplant(i)
         end do

         call ensure_aero_out_capacity(aero, ncoh)

         !----- The patch's canopy-air top and roughness, fixed for the day. Each forcing sample is   !
         !      moved to that top along the patch's log profile (met_to_cas_top,                     !
         !      docs/science/forcing.md §8), and the aerodynamics runs from there. -----------------!
         z_top = biophys%cas%can_depth
         call canopy_roughness(ctx%col_config%aero, ageom, rough_p, displace_p)

         !----- The sub-steps, each on its own forcing sample (the reference climate without forcing). !
         budget = column_budget_t()
         do isub = 1_ik, cfg%n_fast_per_slow
            !----- With forcing, the ground shortwave set here is replaced by the two-stream's in      !
            !      apply_rt_forcing; without, it is the reference value. --------------------------!
            if (do_forcing) then
               met = met_site(isub, poly%patch%owner_site(ip))
               met_top = met_to_cas_top(met, cfg%forcing, z_top, displace_p, rough_p)
               sw_ground = f_ground * met%swdown()
            else
               met_top = met                   ! the reference climate IS the air at the canopy-air top
               sw_ground = ctx%rad_sw_ground
            end if
            !----- Accumulate the sub-step air temperature for the daily-mean phenology driver. ------!
            red_site(RED_PHENO_TAIR, isub, ip) = met%tair_k
            call fill_forcing(forc, col_cohort, met_top, sw_ground, sum_lai)
            red_site(RED_RAINF,    isub, ip) = patch_polygon_area(poly, ip) * forc%rainfall
            red_site(RED_SNOWFALL, isub, ip) = patch_polygon_area(poly, ip) * forc%snowfall
            !----- With forcing, the two-stream canopy radiation replaces the LAI-share split          !
            !      (docs/science/canopy_radiation_transfer.md §4). ------------------------------------!
            if (do_forcing) call apply_rt_forcing(forc, ncoh, col_cohort%pft, col_cohort%lai, col_cohort%wai,  &
                                 col_cohort%height, biophys%leaf_temp, biophys%wood_temp,                      &
                                 biophys%soil_energy%soil_temp(1), biophys%snow, ctx%col_config%snow,               &
                                 ctx%col_config%aero%snow_rough, ctx%soil_albedo, ctx%soil_emiss, ctx%rad_opt, met, &
                                 cfg%leaf_absorptance)
            call fill_aenv(aenv, biophys, met_top, z_top)
            !----- The per-cohort outputs go as 1:ncoh slices of the pool, as column_fast_step takes them  !
            !      assumed-shape. The diagnostic capture is passed only when the run reports it. ---------!
            if (do_cdiag) then
               cdiag_buf(:, 1:ncoh) = 0.0_wp
               call column_fast_step(cfg%dt_fast, ctx%col_config, site_col(poly%patch%owner_site(ip)), aenv, ageom, col_cohort, &
                                     forc, biophys, aero, budget,                                                     &
                                     fast_ws, gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux, cdiag=cdiag_buf(:, 1:ncoh))
            else
               call column_fast_step(cfg%dt_fast, ctx%col_config, site_col(poly%patch%owner_site(ip)), aenv, ageom, col_cohort, &
                                     forc, biophys, aero, budget,                                                     &
                                     fast_ws, gpp_coh=gpp_coh(1:ncoh), leaf_resp_coh=leaf_resp_coh(1:ncoh),            &
                                     psi_leaf_coh=psi_leaf_coh(1:ncoh),                                        &
                                     stem_resp_coh=stem_resp_coh(1:ncoh), root_resp_coh=root_resp_coh(1:ncoh), &
                                     le_flux=le_flux, h_flux=h_flux)
            end if
            !----- The area-weighted CAS -> atmosphere latent flux, integrated into ET [kg/m2 = mm]. -!
            red_site(RED_ET, isub, ip) = patch_polygon_area(poly, ip) * (le_flux / latent_heat_vap) * cfg%dt_fast
            !----- And the water leaving its soil column: surface runoff and bottom drainage [kg/m2]. ----!
            red_site(RED_RUNOFF_SURF, isub, ip) = patch_polygon_area(poly, ip) * budget%runoff_surf_export
            red_site(RED_DRAINAGE,    isub, ip) = patch_polygon_area(poly, ip) * budget%drainage_export
            !----- The solvers' work counters, area-weighted like every other polygon diagnostic. -----!
            red_site(RED_INTEG_STEPS,  isub, ip) = patch_polygon_area(poly, ip) * real(budget%integ_nsteps,   wp)
            red_site(RED_INTEG_REJ,    isub, ip) = patch_polygon_area(poly, ip) * real(budget%integ_nrej,     wp)
            red_site(RED_WATER_NSUB,   isub, ip) = patch_polygon_area(poly, ip) * real(budget%water_nsub,   wp)
            red_site(RED_WATER_NEWTON, isub, ip) = patch_polygon_area(poly, ip) * real(budget%water_newton, wp)
            red_site(RED_WATER_REJECT, isub, ip) = patch_polygon_area(poly, ip) * real(budget%water_reject, wp)
            red_site(RED_WATER_NFAIL,  isub, ip) = patch_polygon_area(poly, ip) * real(budget%water_nfail,  wp)
            !----- Integrator health: a nonzero RK45 rescue count means ARK took some of an RK45 run's   !
            !      dt_fast steps. --------------------------------------------------------------------!
            red_site(RED_RK45_RESCUE,  isub, ip) = patch_polygon_area(poly, ip) * real(budget%rk45_rescue,    wp)
            red_site(RED_CLAMP_STAGE,  isub, ip) = patch_polygon_area(poly, ip) * real(budget%clamp_stage_n,  wp)
            !----- The day's integrals of the per-pool environment scalar and the matrix Rh, which the   !
            !      daily soil_carbon_step consumes (docs/science/soil_carbon.md §6). dt_fast_days is a   !
            !      private routine variable: nvfortran rejects a BLOCK here (CLAUDE.md compiler traps). -!
            if (cfg%soil_carbon%soil_carbon_on) then
               dt_fast_days = cfg%dt_fast / day_sec
               poly%patch%xi_accum(ip)%fast_grnd   = poly%patch%xi_accum(ip)%fast_grnd               &
                                                     + budget%xi_step(IP_FAST_GRND)   * dt_fast_days
               poly%patch%xi_accum(ip)%fast_soil   = poly%patch%xi_accum(ip)%fast_soil               &
                                                     + budget%xi_step(IP_FAST_SOIL)   * dt_fast_days
               poly%patch%xi_accum(ip)%struct_grnd = poly%patch%xi_accum(ip)%struct_grnd             &
                                                     + budget%xi_step(IP_STRUCT_GRND) * dt_fast_days
               poly%patch%xi_accum(ip)%struct_soil = poly%patch%xi_accum(ip)%struct_soil             &
                                                     + budget%xi_step(IP_STRUCT_SOIL) * dt_fast_days
               poly%patch%xi_accum(ip)%microbial   = poly%patch%xi_accum(ip)%microbial               &
                                                     + budget%xi_step(IP_MICR)        * dt_fast_days
               poly%patch%xi_accum(ip)%slow        = poly%patch%xi_accum(ip)%slow                    &
                                                     + budget%xi_step(IP_SLOW)        * dt_fast_days
               poly%patch%xi_accum(ip)%passive     = poly%patch%xi_accum(ip)%passive                 &
                                                     + budget%xi_step(IP_PASSIVE)     * dt_fast_days
               poly%patch%xi_accum(ip)%rh_fast_accum = poly%patch%xi_accum(ip)%rh_fast_accum         &
                                                     + budget%rh_matrix_step * dt_fast_days
            end if
            !----- The opt-in sub-daily probe. It writes one shared file, so the configuration check    !
            !      rejects it with n_threads > 1. -------------------------------------------------!
            if (cfg%fast_probe .and. do_forcing)                                                    &
               call write_fast_probe(cfg, t_sample(isub), ip, ncoh, biophys, col_cohort, gpp_coh(1:ncoh),      &
                                     le_flux, met%swdown())
            !----- Patch GPP and net absorbed radiation for the patch row. h_flux comes from          !
            !      column_fast_step, as le_flux does: formed here from post-call reads, nvfortran      !
            !      miscompiles it to 0 (the array-valued-result class of trap, docs/building.md). -----!
            gpp_patch = sum(gpp_coh(1:ncoh) * col_cohort%nplant(1:ncoh))    ! [umol/m2/s]
            rnet      = forc%abs_sw_ground + forc%abs_lw_ground                                          &
                        + sum(forc%abs_sw(1:ncoh)) + sum(forc%abs_lw(1:ncoh))
            !----- The patch row for this sub-step, accumulated and staged as it is. NPP is GPP net of    !
            !      maintenance respiration only: the slow allocator charges growth respiration. ---------!
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
            !----- FAST-tier staging in the patch's and its cohorts' own slots, disjoint across patches. -!
            if (do_fast) then
               out_bufs%fast_patch(:, ip, isub) = prow
               kb = (ip - 1_ik) * nl
               out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)  = biophys%soil_energy%soil_temp(1:nl)
               out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub) = biophys%soil_water%theta(1:nl)
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  out_bufs%fast_cohort_gpp(i,isub)    = gpp_coh(j)
                  out_bufs%fast_cohort_height(i,isub) = col_cohort%height(j)
                  if (do_cdiag) out_bufs%fast_cohort_diag(i, :, isub) = cdiag_buf(:, j)
               end do
            end if
            !----- Fold the sub-step's diagnostics into the dt-weighted accumulators (meds_diag_vars).  !
            !      The cohort and patch slots are disjoint across patches, so no reduction is needed. -!
            if (do_cdiag) then
               do j = 1_ik, ncoh
                  i = i0 + j - 1_ik
                  poly%cohort%diag%v(:, i) = poly%cohort%diag%v(:, i) + cdiag_buf(:, j) * cfg%dt_fast
                  poly%cohort%diag%w(i)    = poly%cohort%diag%w(i)    + cfg%dt_fast
               end do
            end if
            if (do_pdiag) call accumulate_patch_diag(poly%patch%diag, ip, cfg%dt_fast, prow)
            !----- Gross GPP and the maintenance losses [umol/plant/s -> kgC/plant], kept apart: the    !
            !      slow step's compute_carbon_allocation nets them. ---------------------------------!
            do j = 1_ik, ncoh
               i = i0 + j - 1_ik
               poly%cohort%gpp_accum(i)       = poly%cohort%gpp_accum(i)       + gpp_coh(j)       * cfg%dt_fast * umol_2_kgc
               poly%cohort%leaf_resp_accum(i) = poly%cohort%leaf_resp_accum(i) + leaf_resp_coh(j) * cfg%dt_fast * umol_2_kgc
               poly%cohort%stem_resp_accum(i) = poly%cohort%stem_resp_accum(i) + stem_resp_coh(j) * cfg%dt_fast * umol_2_kgc
               poly%cohort%root_resp_accum(i) = poly%cohort%root_resp_accum(i) + root_resp_coh(j) * cfg%dt_fast * umol_2_kgc
               !----- Hours of light at the cohort's top, the phenology light cue: a step counts when  !
               !      the PAR there exceeds the PFT's par_min. ---------------------------------------!
               if (forc%par_top(j) * par_w_2_umol > cfg%pft%pheno(poly%cohort%pft(i))%par_min)           &
                  poly%cohort%light_hours_accum(i) = poly%cohort%light_hours_accum(i) + cfg%dt_fast / 3600.0_wp
               !----- Running daily maximum of psi_leaf, the predawn value that drives tomorrow's    !
               !      stomatal control (docs/science/leaf_gas_exchange.md §4.2). -------------------!
               poly%cohort%dmax_psi_leaf_accum(i) = max(poly%cohort%dmax_psi_leaf_accum(i), psi_leaf_coh(j))
            end do
         end do

         !----- Write the evolved state back: the patch's stores and each cohort's tissue state. --!
         poly%patch%cas(ip)    = biophys%cas
         poly%patch%soil_energy(ip) = biophys%soil_energy
         poly%patch%soil_water(ip) = biophys%soil_water
         poly%patch%snow(ip)   = biophys%snow
         poly%patch%adapt_dt_last(ip) = biophys%adapt_dt_last
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            poly%cohort%leaf_energy(i) = biophys%leaf_energy(j) / poly%cohort%nplant(i)
            poly%cohort%wood_energy(i) = biophys%wood_energy(j) / poly%cohort%nplant(i)
            poly%cohort%leaf_water_mass(i) = biophys%leaf_water_mass(j)
            poly%cohort%wood_water_mass(i) = biophys%wood_water_mass(j)
            poly%cohort%leaf_surf_water(i) = biophys%leaf_surf_water(j) / poly%cohort%nplant(i)
            poly%cohort%wood_surf_water(i) = biophys%wood_surf_water(j) / poly%cohort%nplant(i)
            call set_cohort_tissue_temps(poly%cohort, i)
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
      !  Fold the staged cells into the polygon in patch order, patch outer and sub-step inner, !
      !  as a serial loop sums them, so the totals do not move with the thread count            !
      !  (docs/science/numerical_scheme.md §6a).                                                !
      !=========================================================================================!
      poly%site%pheno_tair_sum = 0.0_wp
      poly%site%pheno_tair_n = 0_ik
      poly%et_accum       = 0.0_wp
      poly%runoff_surf_accum = 0.0_wp
      poly%drainage_accum = 0.0_wp
      poly%work_integ_steps = 0.0_wp
      poly%work_integ_rej  = 0.0_wp
      poly%work_water_nsub  = 0.0_wp
      poly%work_water_newton = 0.0_wp
      poly%work_water_nfail = 0.0_wp
      poly%work_water_reject = 0.0_wp
      poly%work_rk45_rescue = 0.0_wp
      poly%work_clamp_stage  = 0.0_wp
      do ip = 1_ik, npatch
         do isub = 1_ik, nsub
            poly%et_accum          = poly%et_accum          + red_site(RED_ET,            isub, ip)
            poly%runoff_surf_accum = poly%runoff_surf_accum + red_site(RED_RUNOFF_SURF,   isub, ip)
            poly%drainage_accum    = poly%drainage_accum    + red_site(RED_DRAINAGE,      isub, ip)
            associate (is_p => poly%patch%owner_site(ip))
               poly%site%pheno_tair_sum(is_p) = poly%site%pheno_tair_sum(is_p) + red_site(RED_PHENO_TAIR, isub, ip)
               poly%site%pheno_tair_n(is_p)   = poly%site%pheno_tair_n(is_p) + 1_ik
            end associate
            poly%work_integ_steps  = poly%work_integ_steps  + red_site(RED_INTEG_STEPS,   isub, ip)
            poly%work_integ_rej    = poly%work_integ_rej    + red_site(RED_INTEG_REJ,     isub, ip)
            poly%work_water_nsub   = poly%work_water_nsub   + red_site(RED_WATER_NSUB,   isub, ip)
            poly%work_water_newton = poly%work_water_newton + red_site(RED_WATER_NEWTON, isub, ip)
            poly%work_water_reject = poly%work_water_reject + red_site(RED_WATER_REJECT, isub, ip)
            poly%work_water_nfail  = poly%work_water_nfail  + red_site(RED_WATER_NFAIL,  isub, ip)
            poly%work_rk45_rescue  = poly%work_rk45_rescue  + red_site(RED_RK45_RESCUE,   isub, ip)
            poly%work_clamp_stage  = poly%work_clamp_stage  + red_site(RED_CLAMP_STAGE,   isub, ip)
         end do
      end do
      !----- FAST-tier means, folded over the patches in increasing order for each sub-step: the       !
      !      polygon's patch rows, each patch weighted by its share of the polygon, and each site's    !
      !      soil column, each of its patches weighted by its area in the site. -----------------------!
      if (do_fast) then
         do isub = 1_ik, nsub
            do ip = 1_ik, npatch
               kb = (ip - 1_ik) * nl
               ks = (poly%patch%owner_site(ip) - 1_ik) * nl
               out_bufs%fast_site(:, isub) = out_bufs%fast_site(:, isub)                                  &
                                             + patch_polygon_area(poly, ip) * out_bufs%fast_patch(:, ip, isub)
               out_bufs%fast_soil_temp(ks+1_ik:ks+nl, isub)  = out_bufs%fast_soil_temp(ks+1_ik:ks+nl, isub)  &
                                             + poly%patch%area(ip) * out_bufs%fast_soil_temp_patch(kb+1_ik:kb+nl, isub)
               out_bufs%fast_soil_water(ks+1_ik:ks+nl, isub) = out_bufs%fast_soil_water(ks+1_ik:ks+nl, isub) &
                                             + poly%patch%area(ip) * out_bufs%fast_soil_water_patch(kb+1_ik:kb+nl, isub)
            end do
         end do
      end if
      !----- The forcing echo (PY_*) per sub-step: the sample, with the precipitation split as the    !
      !      patches received it (docs/science/forcing.md §7), staged and accumulated. --------------!
      if (do_fast .or. poly%diag%active) then
         do isub = 1_ik, nsub
            rain_s = 0.0_wp
            snow_s = 0.0_wp
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
            if (poly%diag%active) call accumulate_polygon_diag(poly%diag, cfg%dt_fast, echo)
         end do
      end if
      !----- Each site's forcing (SI_*), for the coarse tiers. ---------------------------------------!
      if (poly%site%diag%active) then
         do isub = 1_ik, nsub
            do is = 1_ik, poly%site%n
               if (do_forcing) then
                  call accumulate_site_diag(poly%site%diag, is, cfg%dt_fast, met_site(isub, is))
               else
                  call accumulate_site_diag(poly%site%diag, is, cfg%dt_fast, met_ref)
               end if
            end do
         end do
      end if

      if (present(worst_energy)) then
         worst_energy = 0.0_wp
         do ip = 1_ik, npatch
            worst_energy = max(worst_energy, red_worst_energy(ip))
         end do
      end if
      if (present(worst_water)) then
         worst_water = 0.0_wp
         do ip = 1_ik, npatch
            worst_water = max(worst_water, red_worst_water(ip))
         end do
      end if
      if (present(n_budget_fail)) n_budget_fail = sum(red_nfail(1:npatch))
      if (present(run_energy_budget) .or. present(run_water_budget) .or. present(run_face_budget)) then
         site_energy_budget = budget_t()
         site_water_budget = budget_t()
         site_face_budget = budget_t()
         do ip = 1_ik, npatch
            call budget_merge(site_energy_budget, red_budget_energy(ip), patch_polygon_area(poly, ip))
            call budget_merge(site_water_budget,  red_budget_water(ip),  patch_polygon_area(poly, ip))
            call budget_merge(site_face_budget,   red_budget_face(ip),   patch_polygon_area(poly, ip))
         end do
         if (present(run_energy_budget)) call budget_merge(run_energy_budget, site_energy_budget, 1.0_wp)
         if (present(run_water_budget))  call budget_merge(run_water_budget,  site_water_budget,  1.0_wp)
         if (present(run_face_budget))   call budget_merge(run_face_budget,   site_face_budget,   1.0_wp)
      end if
      if (do_fast) out_bufs%fast_ready = .true.   ! signal main: the staging holds this slow step
   end subroutine fast_dynamics

   !----- The opt-in sub-daily probe ([fast].fast_probe): one CSV row per (patch, sub-step), flushed  !
   !      as written, from a saved unit opened on first use. The FAST tier is the supported sub-daily  !
   !      output; the probe serves integrator and dt_fast checks. -----------------------------------!
   subroutine write_fast_probe(cfg, t_sub, ip, ncoh, biophys, col_cohort, gpp_coh, le_flux, sw_in)
      type(meds_config_t),   intent(in) :: cfg
      type(meds_time_t),     intent(in) :: t_sub
      integer(ik),           intent(in) :: ip, ncoh
      type(patch_biophys_t), intent(in) :: biophys
      type(column_cohort_t), intent(in) :: col_cohort
      real(wp),              intent(in) :: gpp_coh(:), le_flux, sw_in
      integer, save :: unit           ! newunit returns a negative handle: track open state apart
      logical, save :: opened = .false.
      real(wp)      :: gpp_patch, leaf_temp_mean, wood_temp_mean

      if (.not. opened) then
         open(newunit=unit, file=trim(cfg%fast_probe_file), status='replace', action='write')
         write(unit,'(a)') 'datetime,patch,sw_in_W_m2,cas_temp_K,gpp_umol_m2_s,le_W_m2,soil_temp_top_K,leaf_temp_K,wood_temp_K'
         opened = .true.
      end if

      gpp_patch = 0.0_wp
      leaf_temp_mean = 0.0_wp
      wood_temp_mean = 0.0_wp
      if (ncoh > 0_ik) then
         gpp_patch      = sum(gpp_coh(1:ncoh) * col_cohort%nplant(1:ncoh))    ! per-plant [umol/plant/s] x nplant -> [umol/m2/s]
         leaf_temp_mean = sum(biophys%leaf_temp(1:ncoh)) / real(ncoh, wp)
         wood_temp_mean = sum(biophys%wood_temp(1:ncoh)) / real(ncoh, wp)
      end if

      write(unit,'(a,",",i0,7(",",es13.6))') trim(time_to_string(t_sub)), ip,   &
            sw_in, biophys%cas%can_temp, gpp_patch, le_flux, biophys%soil_energy%soil_temp(1), leaf_temp_mean, wood_temp_mean
      flush(unit)
   end subroutine write_fast_probe

   !----- Grow-only capacity of the per-patch forcing buffers, pre-sized to the largest patch; their  !
   !      writers loop to col_cohort%n, so a larger capacity is exact. ------------------------------!
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

   !----- Fill the patch's column forcing from the sub-step's forcing record and ground shortwave. --!
   subroutine fill_forcing(forc, col_cohort, met, sw_ground, sum_lai)
      type(column_forcing_t), intent(inout) :: forc
      type(column_cohort_t),  intent(in)    :: col_cohort
      type(met_forcing_t),    intent(in)    :: met
      real(wp),               intent(in)    :: sw_ground, sum_lai
      integer(ik) :: j
      forc%enthalpy_atm  = cas_temp_to_enthalpy(met%tair_k, met%qair)
      forc%shv_atm       = met%qair
      forc%co2_atm       = met%co2
      forc%abs_sw_ground = sw_ground
      forc%abs_lw_ground = 0.0_wp
      !----- The rain/snow split at the record's temperature: with forcing, the patch's canopy-air top,  !
      !      where the precipitation's enthalpy is valued too (docs/science/forcing.md §7). --------------!
      call precip_phase(met%precip, met%tair_k, forc%rainfall, forc%snowfall)
      forc%air_temp          = met%tair_k                ! rainfall enthalpy reference (snow/rain-on-snow)
      forc%par_per_w     = 2.1_wp                    ! LAI-split path: total-SW->PAR blend (abs_par == abs_sw)
      !----- Split the canopy-top shortwave across the cohorts by LAI share; with forcing on,     !
      !      apply_rt_forcing replaces it with the two-stream solve. --------------------------!
      do j = 1_ik, col_cohort%n
         if (sum_lai > tiny_num) then
            forc%abs_sw(j) = met%swdown() * col_cohort%lai(j) / sum_lai
         else
            forc%abs_sw(j) = 0.0_wp
         end if
         forc%abs_par(j) = forc%abs_sw(j)            ! LAI path: no PAR/NIR split, PAR is the whole SW
         forc%par_top(j) = met%par_beam + met%par_diffuse   ! no light gradient in the LAI path
         forc%abs_lw(j)  = 0.0_wp
         forc%abs_sw_wood(j) = 0.0_wp                ! no-forcing path: no wood absorption
         forc%abs_lw_wood(j) = 0.0_wp
      end do
   end subroutine fill_forcing

   !----- The context's reference climate as a forcing record. Without forcing only the total       !
   !      shortwave is read (fill_forcing's LAI split), so it all sits in one stream. ---------------!
   pure function reference_met(ctx) result(met)
      type(fast_context_t), intent(in) :: ctx
      type(met_forcing_t) :: met
      met%tair_k   = ctx%air_temp
      met%qair = ctx%shv_atm
      met%psurf_pa = ctx%press
      met%rho_air  = ctx%rho_air
      met%co2  = ctx%co2_atm
      met%wind     = ctx%u_ref
      met%precip   = ctx%precip
      met%par_beam = ctx%rad_sw_top
      met%par_diffuse = 0.0_wp
      met%nir_beam = 0.0_wp
      met%nir_diffuse = 0.0_wp
   end function reference_met

   !----- Run the patch's two-stream canopy radiation and replace the LAI split in `forc` with each   !
   !      cohort's absorbed SW, incident-equivalent PAR and net LW, the ground's SW and LW, and the  !
   !      top-of-canopy fluxes (docs/science/canopy_radiation_transfer.md §4). Cohorts arrive       !
   !      tallest first and the solver takes them bottom to top, so outputs go back through `perm`. -!
   subroutine apply_rt_forcing(forc, ncoh, pft, lai, wai, height, leaf_temp, wood_temp,                &
                               soil_temp_top, snow, snow_params, snow_rough, soil_albedo, soil_emiss,   &
                               rad_opt, met, leaf_absorptance)
      type(column_forcing_t),  intent(inout) :: forc
      integer(ik),             intent(in)    :: ncoh
      integer(ik),             intent(in)    :: pft(:)           !< per-cohort PFT (gather order, top first)
      real(wp),                intent(in)    :: lai(:), wai(:)   !< [m2/m2] leaf / wood area index
      real(wp),                intent(in)    :: height(:)        !< [m] cohort height
      real(wp),                intent(in)    :: leaf_temp(:), wood_temp(:)   !< [K] lagged tissue temperatures
      real(wp),                intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground emission)
      type(snow_column_t),     intent(in)    :: snow             !< pack store (albedo / emissivity ramp)
      type(snow_params_t),     intent(in)    :: snow_params
      real(wp),                intent(in)    :: snow_rough       !< [m] snow-surface roughness ([aerodynamics])
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

      !----- A bare patch needs no special case: canopy_radiation's empty canopy gives the ground's. -!

      !----- perm: gather indices in ascending height (bottom -> top). ------------------------------!
      call ascending_order(height, ncoh, perm)
      do j = 1_ik, ncoh
         i = perm(j)
         pft_bt(j) = pft(i)
         lai_bt(j) = lai(i)
         hgt_bt(j) = height(i)
         wai_bt(j) = wai(i)
         !----- LW emission temperature: T_rad = (f_l*T_leaf^4 + (1-f_l)*T_wood^4)^(1/4), f_l =        !
         !      LAI/(LAI+WAI), at the start of the sub-step (docs/science/canopy_radiation_transfer.md  !
         !      §3); abs_lw is net of that emission (docs/science/vegetation_energy_dynamics.md, "The   !
         !      emission base"). -----------------------------------------------------------------------!
         lf_bt      = lai(i) / max(lai(i) + wai(i), tiny_num)
         tcan_bt(j) = (lf_bt * leaf_temp(i) ** 4                                             &
                       + (1.0_wp - lf_bt) * wood_temp(i) ** 4) ** 0.25_wp
      end do

      !----- rad_forcing_t from met (all W/m2, direct assignment). --------------------------------!
      call alloc_rad_forcing(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = met%cosz
      rf%incid_beam(RAD_VIS) = met%par_beam
      rf%incid_diff(RAD_VIS) = met%par_diffuse
      rf%incid_beam(RAD_NIR) = met%nir_beam
      rf%incid_diff(RAD_NIR) = met%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp
      rf%incid_diff(RAD_LW)  = met%lwdown

      !----- Ground optics: bare-soil albedo/emissivity from [soil], skin temperature from the top soil node. --!
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = soil_albedo
      surf%soil_emiss = soil_emiss
      surf%soil_temp   = soil_temp_top
      !----- Snow blends the ground's albedo, emissivity and temperature toward the pack's by the      !
      !      Niu-Yang 2007 snow-cover fraction, so a partial pack is continuous; fresh and aged albedos !
      !      are weighted by the surface liquid fraction (docs/science/canopy_radiation_transfer.md §4). !
      if (snow%nlayer >= 1_ik .and. snow%swe(1) > snow_params%tiny_snow_mass) then
         associate (sp => snow_params)
            snow_fc = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), sp, snow_rough)
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
         forc%abs_sw(ig)  = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! absorbed leaf SW (energy balance)
         forc%abs_par(ig) = flux%abs_leaf(RAD_VIS, j) / max(leaf_absorptance, tiny_num)  ! incident-equivalent PAR
         forc%abs_lw(ig)  = flux%abs_leaf(RAD_LW, j)                                ! net leaf LW, emission at T_rad included
         forc%abs_sw_wood(ig) = flux%abs_wood(RAD_VIS, j) + flux%abs_wood(RAD_NIR, j)  ! absorbed wood SW
         forc%abs_lw_wood(ig) = flux%abs_wood(RAD_LW, j)                               ! net wood LW
         forc%par_top(ig)     = flux%incid_top(RAD_VIS, j)                             ! PAR at the cohort's top
      end do
      forc%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                         + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forc%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)          ! net ground LW, emission included
      forc%par_per_w     = par_w_2_umol                        ! true VIS absorbed -> photon flux
      !----- Top-of-canopy fluxes: each band's albedo times the incident it was formed from, kept as   !
      !      fluxes (docs/science/canopy_radiation_transfer.md, Surface radiative fluxes on the       !
      !      diagnostic path). ------------------------------------------------------------------!
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
      aenv%u_ref = met%wind
      aenv%zref = zref
      aenv%press = met%psurf_pa
      aenv%rho_air = met%rho_air
      !----- Both halves go through set_aero_env_atm and set_aero_env_canopy, as in the column tests,  !
      !      so `theta_atm` is never left at its default. ----------------------------------------!
      call set_aero_env_atm(aenv, met%tair_k, met%qair, met%co2)
      call set_aero_env_canopy(aenv, biophys%cas%can_temp, biophys%cas%can_shv, biophys%cas%can_co2,           &
                               biophys%soil_energy%soil_temp(1))
   end subroutine fill_aenv

   !=======================================================================================!
   !  The patch block's row for one (patch, sub-step): the numbers the physics used, none  !
   !  recomputed (docs/science/diagnostics.md §4, The patch row). The slow rows stay 0.    !
   !=======================================================================================!
   pure subroutine patch_diag_row(row, dt, le_flux, h_flux, rnet, sw_ground, lw_ground,             &
                                  ustar, ggnet, rough, displace, cas_temp, cas_shv, cas_co2, gpp, nee, &
                                  npp, soil_temp_top, resid_energy, resid_water, wind_top, tair_top,  &
                                  z_top, sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up, w_surface)
      real(wp), intent(out) :: row(N_PA_ROWS)
      real(wp), intent(in)  :: dt                        !< [s]        the sub-step
      real(wp), intent(in)  :: le_flux, h_flux           !< [W/m2]     CAS -> atmosphere
      real(wp), intent(in)  :: rnet                      !< [W/m2]     net all-wave radiation absorbed
      real(wp), intent(in)  :: sw_ground, lw_ground      !< [W/m2]     ground SW / net LW
      real(wp), intent(in)  :: ustar, ggnet              !< [m/s]      friction velocity, ground conductance
      real(wp), intent(in)  :: rough, displace           !< [m]        roughness, displacement height
      real(wp), intent(in)  :: cas_temp, cas_shv, cas_co2 !< [K],[kg/kg],[umol/mol] canopy air
      real(wp), intent(in)  :: gpp, nee, npp             !< [umol/m2/s] gross uptake, net exchange (+ to atm), NPP
      real(wp), intent(in)  :: soil_temp_top             !< [K]        top soil-node temperature
      real(wp), intent(in)  :: resid_energy, resid_water !< [J/m2],[kg/m2] this step's signed ledger residuals
      real(wp), intent(in)  :: wind_top, tair_top        !< [m/s],[K] the forcing at the canopy-air top
      real(wp), intent(in)  :: z_top                     !< [m]   that top (the canopy-air depth)
      !----- Top-of-canopy radiative fluxes per band. -----------------------------------------!
      real(wp), intent(in)  :: sw_in_vis, sw_in_nir, sw_up_vis, sw_up_nir, lw_up
      real(wp), intent(in)  :: w_surface                 !< [kg/m2]   ponded surface water
      row = 0.0_wp
      row(PA_LE)           = le_flux
      row(PA_H)            = h_flux
      row(PA_RNET)         = rnet
      row(PA_SW_GROUND)    = sw_ground
      row(PA_LW_GROUND)    = lw_ground
      row(PA_SW_IN_VIS)    = sw_in_vis
      row(PA_SW_IN_NIR)    = sw_in_nir
      row(PA_SW_UP_VIS)    = sw_up_vis
      row(PA_SW_UP_NIR)    = sw_up_nir
      row(PA_SW_UP)        = sw_up_vis + sw_up_nir
      row(PA_LW_UP)        = lw_up
      row(PA_SKIN_TEMP)    = sqrt(sqrt(max(lw_up, 0.0_wp) / stefan_boltzmann))
      row(PA_USTAR)        = ustar
      row(PA_GGNET)        = ggnet
      row(PA_ROUGH)        = rough
      row(PA_DISPLACE)     = displace
      row(PA_CAS_TEMP)     = cas_temp
      row(PA_CAS_SHV)      = cas_shv
      row(PA_CAS_CO2)      = cas_co2
      !----- VPD is nonlinear in temperature, so it is formed per sub-step, not from mean states. --!
      row(PA_CAS_VPD)      = specific_humidity_to_vpd(cas_temp, cas_shv, p_std)
      row(PA_W_SURFACE)    = w_surface
      !----- The patch's own forcing, at its canopy-air top (docs/science/forcing.md §8). ---------!
      row(PA_WIND_CAS_TOP) = wind_top
      row(PA_TAIR_CAS_TOP) = tair_top
      row(PA_Z_CAS_TOP)    = z_top
      !----- nee drives the CAS CO2 box (positive to the atmosphere); Reco = NEE + GPP. -----------!
      row(PA_GPP)          = gpp
      row(PA_NEE)          = nee
      row(PA_NPP)          = npp
      row(PA_RECO)         = nee + gpp
      !----- The canopy air's whole evaporative flux to the atmosphere [kg/m2/s], not transpiration  !
      !      alone; each cohort's own transpiration is CO_TRANSP. -------------------------------!
      row(PA_TRANSP)       = le_flux / latent_heat_vap
      row(PA_SOIL_TEMP_TOP) = soil_temp_top
      !----- The whole-column ledgers' signed residuals as leak rates [W/m2, kg/m2/s], positive when  !
      !      store appears from nowhere. --------------------------------------------------------!
      row(PA_RESID_ENERGY) = resid_energy / dt
      row(PA_RESID_WATER)  = resid_water / dt
   end subroutine patch_diag_row

   !----- Fold one (patch, sub-step) row into the patch accumulator, dt-weighted; read as v / w. ---!
   subroutine accumulate_patch_diag(pd, ip, dt, row)
      type(patch_diag_block), intent(inout) :: pd
      integer(ik),            intent(in)    :: ip
      real(wp),               intent(in)    :: dt        !< [s] sample weight
      real(wp),               intent(in)    :: row(N_PA_ROWS)
      pd%v(:, ip) = pd%v(:, ip) + row * dt
      pd%w(ip)    = pd%w(ip)    + dt
   end subroutine accumulate_patch_diag

   !----- Size the FAST tier's staging. Sub-step and site counts are run-constant; the patch and     !
   !      cohort buffers grow to at least twice their size when outgrown. The output layer checks   !
   !      the counts against output.patch_max and cohort_max. ---------------------------------------!
   subroutine size_fast_staging(b, nsub, npatch, nsite, ncohort, with_cdiag)
      type(output_buffers_t), intent(inout) :: b
      integer(ik),            intent(in)    :: nsub, npatch, nsite, ncohort
      logical,                intent(in)    :: with_cdiag   !< the per-cohort diagnostic block is on
      integer(ik) :: nl, np, nc, ns
      nl = n_soil_layer_max
      ns = max(nsite, 1_ik)
      if (.not. allocated(b%fast_time))                                                            &
         allocate(b%fast_time(nsub), b%fast_forcing(N_PY_ROWS, nsub), b%fast_site(N_PA_ROWS, nsub),    &
                  b%fast_soil_temp(nl * ns, nsub), b%fast_soil_water(nl * ns, nsub), b%fast_soil_live(nl * ns))
      np = max(npatch, 1_ik)
      if (allocated(b%fast_patch)) then
         if (size(b%fast_patch, 2, kind=ik) < np) then
            np = max(np, 2_ik * size(b%fast_patch, 2, kind=ik))
            deallocate(b%fast_patch, b%fast_soil_temp_patch, b%fast_soil_water_patch, b%fast_soil_patch_live)
         end if
      end if
      if (.not. allocated(b%fast_patch))                                                           &
         allocate(b%fast_patch(N_PA_ROWS, np, nsub), b%fast_soil_temp_patch(nl * np, nsub),           &
                  b%fast_soil_water_patch(nl * np, nsub), b%fast_soil_patch_live(nl * np))
      nc = max(ncohort, 1_ik)
      if (allocated(b%fast_cohort_gpp)) then
         if (size(b%fast_cohort_gpp, 1, kind=ik) < nc) then
            nc = max(nc, 2_ik * size(b%fast_cohort_gpp, 1, kind=ik))
            deallocate(b%fast_cohort_gpp, b%fast_cohort_height)
         end if
      end if
      if (.not. allocated(b%fast_cohort_gpp))                                                      &
         allocate(b%fast_cohort_gpp(nc, nsub), b%fast_cohort_height(nc, nsub))
      !----- The per-cohort diagnostic rows, staged only while their FAST-tier variables are on. ---!
      if (with_cdiag) then
         nc = size(b%fast_cohort_gpp, 1, kind=ik)
         if (allocated(b%fast_cohort_diag)) then
            if (size(b%fast_cohort_diag, 1, kind=ik) < nc) deallocate(b%fast_cohort_diag)
         end if
         if (.not. allocated(b%fast_cohort_diag)) allocate(b%fast_cohort_diag(nc, N_CO_FAST_ROWS, nsub))
      end if
   end subroutine size_fast_staging

   !----- The polygon's forcing for one sub-step as the fast loop used it (docs/science/forcing.md). -!
   subroutine accumulate_polygon_diag(d, dt, f)
      type(polygon_diag_block), intent(inout) :: d
      real(wp),                 intent(in)    :: dt     !< [s] sample weight
      real(wp),                 intent(in)    :: f(N_PY_ROWS)   !< the sub-step's forcing echo
      d%v = d%v + f * dt
      d%w = d%w + dt
   end subroutine accumulate_polygon_diag

   !----- One sub-step of a site's forcing, in the site block's table (SI_*). ---------------------!
   subroutine accumulate_site_diag(d, is, dt, met)
      type(site_diag_block), intent(inout) :: d
      integer(ik),           intent(in)    :: is
      real(wp),              intent(in)    :: dt     !< [s] sample weight
      type(met_forcing_t),   intent(in)    :: met    !< the site's forcing
      real(wp) :: f(N_SI_ROWS)
      f(SI_SW_IN)    = met%swdown()
      f(SI_PAR_BEAM) = met%par_beam
      f(SI_NIR_BEAM) = met%nir_beam
      f(SI_LWDOWN)   = met%lwdown
      f(SI_TAIR)     = met%tair_k
      f(SI_QAIR)     = met%qair
      f(SI_PSURF)    = met%psurf_pa
      f(SI_RHO_AIR)  = met%rho_air
      d%v(:, is) = d%v(:, is) + f * dt
      d%w(is)    = d%w(is) + dt
   end subroutine accumulate_site_diag

   !----- The forcing one sub-step used, as the polygon block's row (PY_*): the one list of the     !
   !      forcing echo, accumulated by the coarse tiers and staged by the FAST tier. --------------!
   pure subroutine forcing_echo(met, rainf, snowfall, f)
      type(met_forcing_t), intent(in)  :: met
      real(wp),            intent(in)  :: rainf, snowfall   !< [kg/m2/s] the split the patches used, area-weighted
      real(wp),            intent(out) :: f(N_PY_ROWS)
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
