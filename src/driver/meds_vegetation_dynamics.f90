!==========================================================================================!
! meds_vegetation_dynamics -- MEDS's slow-loop VEGETATION-DYNAMICS driver (the analogue of      !
! ED2's veg_dynamics_driver). It lives in src/driver/ (compiled into meds_aux) and is THE         !
! ORCHESTRATOR: it (1) assembles the per-cohort carbon NPP from the plant seam, (2) turns that     !
! into the demographic RATES by calling the per-individual plant vital-rate KERNELS               !
! (meds_plant_vital_rates: carbon growth rate, Camac mortality, recruitment), and (3) sequences    !
! the demography engine's APPLY-PRIMITIVES (grow / die / age / recruit / fuse-fiss) folding the    !
! calendar cadence into the structural triggers.                                                  !
!                                                                                          !
! This is THE single place plant and demography meet: the driver calls the plant kernels; the      !
! demography engine only APPLIES supplied arrays, so `demography _|_ plant` holds. (The former      !
! update_demography seam was dissolved into this routine. The empirical growth/recruitment LAWS    !
! were removed from Fortran -- they live in the Python example -- so the Fortran model is the       !
! carbon path.)                                                                                    !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : day_sec
   use meds_config,               only : meds_config_t, growth_window_steps
   use meds_allometry,            only : size2leaf_carbon, carbon_to_structure
   use meds_time,                 only : daylength
   use meds_core_state_types,      only : carbon_flux_block, cohort_deriv_alloc, GROWTH_AVG_UNSET
   use meds_core_interface, only : site_t, update_cohort_states, fill_cohort_deriv,            &
                                         update_patch_states, apply_recruitment,                &
                                         apply_patch_disturbance, new_fuse_cohorts,              &
                                         terminate_cohorts, split_cohorts, new_fuse_patches,     &
                                         terminate_patches, sort_cohorts, sort_patches,          &
                                         update_overtopping_lai
   use meds_plant_vital_rates,    only : carbon_growth_rate, camac_mortality, min_cohort_carbon, &
                                         recruitment_contribution
   use meds_plant_interface,      only : get_plant_flux_slow, growth_respiration,                &
                                         carbon_env_t, carbon_demand_t, carbon_npp_t,            &
                                         pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,&
                                         update_phenology, pheno_drives_to_rates
   implicit none
   private

   public :: vegetation_dynamics, advance_leaf_phenology

   !----- Wood is the residual carbon sink: a demand large enough to take all remaining NPP. --!
   real(wp), parameter :: WOOD_DEMAND_BIG   = 1.0e6_wp    !< [kgC/plant]
   real(wp), parameter :: STUB_TISSUE_TEMP  = 298.15_wp   !< [K] 25 degC (no met forcing yet)
   !----- Phenology-OFF flush rate [1/day]: a large relative rate so the flush cap is non-binding !
   !       (the leaf deficit fills as it did before the refactor -- bit-identical non-phenology path).!
   real(wp), parameter :: PHENOLOGY_OFF_FLUSH = 1.0e6_wp  !< [1/day]
   !----- Patch structural dynamics (fusion/fission/disturbance) are an ANNUAL process, so       !
   !       disturbance integrates its yearly rate over this interval, not the per-step dt.        !
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !< [yr]

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the carbon NPP, compute the carbon !
   ! vital rates via the plant kernels, and sequence the demography apply-primitives + cadence. !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, is_new_month, is_new_year, doy)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      integer(ik),         intent(in), optional :: doy   !< day-of-year at the step start (needed if phenology_on)
      real(wp), allocatable    :: mortality(:), recruitment(:,:), npp_repro(:)
      type(carbon_flux_block)  :: npp
      logical                  :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse
      integer(ik)              :: n_window

      !----- 0. Leaf phenology (opt-in): advance the per-cohort governor drives ONE daily step,   !
      !         BEFORE carbon_growth reads them (the folded phenology_driver -- ED2 calls it inside !
      !         the vegetation-dynamics slow loop). Needs the step-start day-of-year. --------------!
      if (cfg%phenology_on) then
         if (.not. present(doy))                                                                 &
            error stop 'vegetation_dynamics: phenology_on=.true. but no doy supplied'
         call advance_leaf_phenology(site, cfg, doy)
      end if

      !----- 1. Carbon NPP from the plant seam (the ONLY plant call). -----------------------!
      call carbon_growth(site, cfg, cfg%dt_years, npp, npp_repro)

      !----- 2. Carbon vital RATES via the plant kernels (PRE-apply, so mortality sees the same !
      !         growth_avg the former carbon_vital_rates did -- behaviour preserved).            !
      call carbon_rates(site, cfg, npp%wood, npp_repro, cfg%dt_years, mortality, recruitment)

      !----- 3. Fold the calendar cadence + demography on/off + fiss/fuse switches. ---------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- 4. Sequence the apply-primitives (was update_demography). ----------------------!
      n_window = growth_window_steps(cfg)
      site%growth_hist_pos = mod(site%growth_hist_pos, n_window) + 1_ik

      !----- Growth + mortality via the TENDENCY seam: compute the per-cohort time-derivatives    !
      !      (carbon flip -> geometry/pool rates + the realized-dbh-rate moving-average ring       !
      !      buffer; mortality -> log-space rate) HERE in the driver, into the site-carried scratch !
      !      bundle, then let the core engine's pure applier advance the state. The ring buffer is  !
      !      refreshed inside compute_slow_derivatives, so mortality's growth_avg dependency holds. !
      call cohort_deriv_alloc(site%deriv, site%cohort%n)
      call compute_slow_derivatives(site, npp, mortality, cfg%dt_years, n_window, site%growth_hist_pos)
      call update_cohort_states(site%cohort, site%deriv, cfg%dt_years, cfg%negligible_nplant)

      !----- Re-sort every step: growth changed heights, so re-establish the tallest-first order  !
      !      the overtopping-LAI sweep + the patch-light profiles depend on (fuse/fiss also sorts, !
      !      but only on the monthly/annual cadence). -------------------------------------------!
      call sort_cohorts(site)

      !----- Slow per-patch state (patch ageing now; the soil-carbon step will join here). -!
      call update_patch_states(site%patch, cfg%dt_years)

      !----- Cohort restructuring (monthly): recruit + fuse/split + sort. -------------------!
      if (do_cohort_fissfuse) then
         call apply_recruitment(site, cfg, recruitment)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg)
         call split_cohorts(site, cfg)
         call sort_cohorts(site)
      end if

      !----- Patch disturbance (annual) then patch restructuring (annual, independent). ----!
      if (do_patch_disturbance) then
         call apply_patch_disturbance(site, cfg, PATCH_DYNAMICS_INTERVAL)
      end if
      if (do_patch_fissfuse) then
         call sort_patches(site)
         call new_fuse_patches(site, cfg)
         call terminate_patches(site, cfg)
         call new_fuse_cohorts(site, cfg)
         call terminate_cohorts(site, cfg)
         call sort_cohorts(site)
      end if

      !----- 5. Refresh the overtopping-LAI competition diagnostic (the stand is sorted -- either  !
      !         by the per-step sort above or by the cadence fuse/fiss). ------------------------!
      call update_overtopping_lai(site)
   end subroutine vegetation_dynamics

   !---------------------------------------------------------------------------------------!
   ! CARBON-mode TENDENCY computer (the carbon analogue of the former apply_growth, split so    !
   ! the ENGINE only applies). Per cohort it adds this step's NPP to the four pools (tentatively, !
   ! in locals -- state is NOT committed here), runs the allometric FLIP wood_carbon->dbh via     !
   ! carbon_to_structure to get the new geometry, and backs out every field's time-derivative     !
   ! d_X_dt = (X_new - X_old)/dt for the core applier. It ALSO advances the moving-average ring    !
   ! buffer with the realized dbh-increment RATE (persistent state, identical eviction logic to    !
   ! the former apply_growth) so Camac mortality sees carbon growth, and records the log-space     !
   ! mortality rate dln_nplant_dt = -mortality. Host-only (the flip calls the branchy wood_to_dbh).!
   !---------------------------------------------------------------------------------------!
   subroutine compute_slow_derivatives(site, npp, mortality, dt_yr, n_window, hist_pos)
      type(site_t),             intent(inout) :: site       ! inout: fills site%deriv + refreshes the ring buffer
      type(carbon_flux_block),  intent(in)    :: npp
      real(wp),                 intent(in)    :: mortality(:)
      real(wp),                 intent(in)    :: dt_yr
      integer(ik),              intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: lc_new, fc_new, wc_new, nc_new
      real(wp)    :: dbh_new, height_new, ba_new, agb_new, la_new, dbh_rate

      associate (cohort => site%cohort)
         do i = 1_ik, cohort%n
            !----- Tentative new pools (never let a pool go negative), in locals. -------------!
            lc_new = max(cohort%leaf_carbon(i)          + npp%leaf(i),          0.0_wp)
            fc_new = max(cohort%fineroot_carbon(i)      + npp%fineroot(i),      0.0_wp)
            wc_new = max(cohort%wood_carbon(i)          + npp%wood(i),          0.0_wp)
            nc_new = max(cohort%nonstructural_carbon(i) + npp%nonstructural(i), 0.0_wp)
            !----- FLIP the geometry from the tentative pools (dbh from wood_carbon, leaf_area    !
            !      from leaf_carbon) into locals -- state is committed by update_cohort_states. --!
            call carbon_to_structure(wc_new, lc_new, cohort%p_wood_density(i), cohort%p_hgt_max(i), &
                                     cohort%p_aboveground_frac(i), cohort%p_sla(i),                 &
                                     dbh_new, height_new, ba_new, agb_new, la_new)
            !----- Back out the tendencies + advance the ring buffer with the REALIZED dbh rate  !
            !      (the shared core builder; same fill/evict logic the empirical path uses). ----!
            dbh_rate = (dbh_new - cohort%dbh(i)) / dt_yr
            call fill_cohort_deriv(cohort, i, site%deriv, dt_yr, mortality(i), dbh_rate,            &
                                   dbh_new, height_new, ba_new, agb_new, la_new,                    &
                                   lc_new, fc_new, wc_new, nc_new, n_window, hist_pos)
         end do
      end associate
   end subroutine compute_slow_derivatives

   !---------------------------------------------------------------------------------------!
   ! CARBON vital-RATE assembler: turn the per-cohort carbon fluxes (npp_wood for the dbh rate, !
   ! npp_repro for recruits) into the mortality [1/yr] + recruitment [plant/m2/yr] arrays, by     !
   ! calling the per-individual plant kernels. The site-level sweep (patch reduction + baseline   !
   ! seed rain + the GROWTH_AVG_UNSET test) stays here in the DRIVER, so the plant kernels remain  !
   ! state-free. Mirrors the former meds_demography_rates%carbon_vital_rates exactly.             !
   !---------------------------------------------------------------------------------------!
   subroutine carbon_rates(site, cfg, npp_wood, npp_repro, dt_yr, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp),              intent(in)  :: npp_wood(:), npp_repro(:), dt_yr
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr] per cohort
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, j, pf, ip
      real(wp)    :: dbh_rate
      real(wp), allocatable :: carbon_min(:)

      n    = site%cohort%n
      np   = site%patch%n
      npft = cfg%pft%n
      allocate(mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp

      associate (cohort => site%cohort, pft => cfg%pft)
         allocate(carbon_min(npft))
         do pf = 1_ik, npft
            carbon_min(pf) = min_cohort_carbon(pft%min_cohort_height, pft%wood_density(pf))
         end do
         !----- Baseline external seed rain (independent of structure). --------------------!
         do ip = 1_ik, np
            do pf = 1_ik, npft
               if (pft%include_pft(pf) == 1_ik) recruitment(pf, ip) = pft%seed_rain_recruits(pf)
            end do
         end do
         do j = 1_ik, n
            pf = cohort%pft(j)
            ip = cohort%owner_patch(j)
            !----- Prospective carbon dbh-increment rate (the mortality instantaneous term). --!
            dbh_rate = carbon_growth_rate(cohort%wood_carbon(j), npp_wood(j), cohort%dbh(j), dt_yr, &
                                          cohort%p_wood_density(j), cohort%p_hgt_max(j),            &
                                          cohort%p_aboveground_frac(j))
            mortality(j) = camac_mortality(cohort%growth_avg(j), dbh_rate,                          &
                                           cohort%growth_avg(j) /= GROWTH_AVG_UNSET,                &
                                           pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))
            !----- Reproduction carbon -> recruits (annual rate; monthly share applied later),    !
            !      gated by include_pft (a disabled PFT must not re-establish).  -----------------!
            if (pft%include_pft(pf) == 1_ik)                                                       &
            recruitment(pf, ip) = recruitment(pf, ip)                                              &
                 + recruitment_contribution(cohort%nplant(j), npp_repro(j), dt_yr,                 &
                                            pft%repro_carbon_efficiency(pf), carbon_min(pf))
         end do
      end associate
   end subroutine carbon_rates

   !---------------------------------------------------------------------------------------!
   ! CARBON NPP assembler. Per cohort it builds the allometric demands (target - pool) and the  !
   ! net carbon (fast-loop GPP - maintenance resp, or the stub gpp_ref when fast is off), minus   !
   ! growth respiration, calls the plant slow-flux seam get_plant_flux_slow, and returns the      !
   ! per-pool NPP block + the reproduction carbon. This is the ONLY place the driver calls plant  !
   ! for fluxes (the vital-rate kernels above are the other plant call).                          !
   !---------------------------------------------------------------------------------------!
   subroutine carbon_growth(site, cfg, dt_yr, npp, npp_repro)
      type(site_t),            intent(in)  :: site
      type(meds_config_t),     intent(in)  :: cfg
      real(wp),                intent(in)  :: dt_yr
      type(carbon_flux_block), intent(out) :: npp
      real(wp), allocatable,   intent(out) :: npp_repro(:)
      integer(ik) :: n, j, pf
      real(wp)    :: leaf_target, a_carbon, gross_gpp, resp_maint, dt_day
      type(carbon_env_t)    :: env
      type(carbon_demand_t) :: demand
      type(carbon_npp_t)    :: out

      n = site%cohort%n
      dt_day = cfg%dt_slow / day_sec            ! days in the slow step (phenology rates are [1/day])
      allocate(npp%leaf(n), npp%fineroot(n), npp%wood(n), npp%nonstructural(n), npp_repro(n))
      associate (cohort => site%cohort, pft => cfg%pft)
         do j = 1_ik, n
            pf = cohort%pft(j)
            !----- Allometric demands (target - current pool); wood is the residual sink. ----!
            leaf_target      = size2leaf_carbon(cohort%dbh(j), cohort%height(j), pft%sla(pf))
            demand%leaf      = max(0.0_wp, leaf_target                              - cohort%leaf_carbon(j))
            demand%fineroot  = max(0.0_wp, pft%root_to_leaf_ratio(pf) * leaf_target - cohort%fineroot_carbon(j))
            demand%storage   = max(0.0_wp, pft%storage_cushion(pf)    * leaf_target - cohort%nonstructural_carbon(j))
            demand%wood      = WOOD_DEMAND_BIG
            demand%reproduction_fraction = merge(pft%reproduction_investment_fraction(pf), 0.0_wp, &
                                                 cohort%height(j) >= pft%min_reproduction_height)
            if (cfg%fast_biophysics_on) then
               gross_gpp  = cohort%gpp_accum(j)
               resp_maint = cohort%leaf_resp_accum(j) + cohort%stem_resp_accum(j)                  &
                          + cohort%root_resp_accum(j)
            else
               gross_gpp  = cfg%gpp_ref * cohort%leaf_area(j) * dt_yr
               resp_maint = 0.0_wp
            end if
            a_carbon       = gross_gpp - resp_maint
            env%net_carbon = a_carbon - growth_respiration(a_carbon, pft%growth_resp_factor(pf))
            env%nonstructural    = cohort%nonstructural_carbon(j)
            env%leaf_carbon      = cohort%leaf_carbon(j)
            env%fineroot_carbon  = cohort%fineroot_carbon(j)
            env%leaf_carbon_full = leaf_target
            env%tissue_temp      = STUB_TISSUE_TEMP
            env%dt_yr            = dt_yr
            env%dt_day           = dt_day
            !----- Phenology (opt-in): derive the two RELATIVE rates from the stored governor      !
            !      drives; OFF => large flush (uncapped fill) + zero active shed (bit-identical).  !
            if (cfg%phenology_on) then
               call pheno_drives_to_rates(cohort%pheno_flush_drive(j), cohort%pheno_shed_drive(j), &
                                          pft%pheno_k_flush_max(pf), pft%pheno_k_shed_max(pf),      &
                                          env%leaf_flush_rate, env%leaf_shed_rate)
            else
               env%leaf_flush_rate = PHENOLOGY_OFF_FLUSH
               env%leaf_shed_rate  = 0.0_wp
            end if
            call get_plant_flux_slow(env, cfg, pf, demand, out)
            npp%leaf(j)          = out%leaf
            npp%fineroot(j)      = out%fineroot
            npp%wood(j)          = out%wood
            npp%nonstructural(j) = out%nonstructural
            npp_repro(j)         = out%repro
         end do
      end associate
   end subroutine carbon_growth

   !---------------------------------------------------------------------------------------!
   ! Advance the leaf phenology of every cohort over one slow step (the folded phenology driver, !
   ! ED2 phenology_driv analogue -- see docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md). It  !
   ! flattens the per-PFT cue params, builds the daily env from the site daily-mean air temperature !
   ! the fast loop accumulated + latitude + day-of-year, advances the two governor drives + thermal !
   ! memory via update_phenology, and writes the drives back to the cohort. It touches NO leaf/      !
   ! storage carbon; carbon_growth derives the two rates from the stored drives. A no-temperature    !
   ! step (no fast sub-steps ran) is skipped so the memory is never advanced on a bogus 0/0 mean.    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_leaf_phenology(site, cfg, doy)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: doy
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: i, pf
      real(wp)    :: dt_days, temp_day, dlen
      logical     :: north

      if (site%pheno_tair_n < 1_ik) return             ! no fast sub-steps this slow step -> no drivers
      dt_days  = cfg%dt_slow / day_sec
      temp_day = site%pheno_tair_sum / real(site%pheno_tair_n, wp)
      north    = cfg%forcing%latitude_deg >= 0.0_wp
      dlen     = daylength(cfg%forcing%latitude_deg, doy)

      do i = 1_ik, site%cohort%n
         pf = site%cohort%pft(i)
         call flatten_pheno_params(cfg, pf, params)

         !----- Daily environment (P1-P2: air temp drives GDD/chilling AND stands in for the cold-  !
         !      drop soil-temp trigger; WATER/HYDRO/LIGHT drivers are unused because those cue bits  !
         !      are rejected by validate_config until P3). ---------------------------------------!
         env%temp_day      = temp_day
         env%soil_temp     = temp_day
         env%daylength     = dlen
         env%doy           = doy
         env%hemis_north   = north
         env%avail_water   = 0.0_wp
         env%dmax_leaf_psi = 0.0_wp
         env%rad           = 0.0_wp

         !----- Pack the cohort's governor + thermal memory, advance, unpack. (Reset the whole state !
         !      each cohort so the unused WATER/HYDRO/LIGHT accumulators cannot leak across cohorts.) !
         state             = pheno_state_t()
         state%flush_drive = site%cohort%pheno_flush_drive(i)
         state%shed_drive  = site%cohort%pheno_shed_drive(i)
         state%gdd         = site%cohort%pheno_gdd(i)
         state%chill       = site%cohort%pheno_chill(i)
         call update_phenology(env, params, dt_days, state, out)
         site%cohort%pheno_flush_drive(i) = state%flush_drive
         site%cohort%pheno_shed_drive(i)  = state%shed_drive
         site%cohort%pheno_gdd(i)         = state%gdd
         site%cohort%pheno_chill(i)       = state%chill
      end do
   end subroutine advance_leaf_phenology

   !---------------------------------------------------------------------------------------!
   ! Flatten the per-PFT phenology traits (cfg%pft%pheno_*) into the self-contained kernel param  !
   ! set. The WATER/HYDRO param fields keep their pheno_params_t defaults (their cues are rejected  !
   ! in P1-P2). Mirrors meds_plant_interface's leaf-trait flattening.                              !
   !---------------------------------------------------------------------------------------!
   subroutine flatten_pheno_params(cfg, ipft, p)
      type(meds_config_t),  intent(in)  :: cfg
      integer(ik),          intent(in)  :: ipft
      type(pheno_params_t), intent(out) :: p
      associate (t => cfg%pft)
         p%flush_cue_mask      = t%pheno_flush_cue_mask(ipft)
         p%shed_cue_mask       = t%pheno_shed_cue_mask(ipft)
         p%cue_sharpness       = t%pheno_cue_sharpness(ipft)
         p%k_flush_max         = t%pheno_k_flush_max(ipft)
         p%k_shed_max          = t%pheno_k_shed_max(ipft)
         p%tau_flush           = t%pheno_tau_flush(ipft)
         p%tau_shed            = t%pheno_tau_shed(ipft)
         p%gdd_base_temp       = t%pheno_gdd_base_temp(ipft)
         p%chill_base_temp     = t%pheno_chill_base_temp(ipft)
         p%phen_a              = t%pheno_phen_a(ipft)
         p%phen_b              = t%pheno_phen_b(ipft)
         p%phen_c              = t%pheno_phen_c(ipft)
         p%cold_drop_daylength = t%pheno_cold_drop_daylength(ipft)
         p%cold_drop_soiltemp1 = t%pheno_cold_drop_soiltemp1(ipft)
         p%cold_drop_soiltemp2 = t%pheno_cold_drop_soiltemp2(ipft)
         p%water_width         = t%pheno_water_width(ipft)
         p%photo_crit          = t%pheno_photo_crit(ipft)
         p%photo_slope         = t%pheno_photo_slope(ipft)
         p%light_on_threshold  = t%pheno_light_on_threshold(ipft)
         p%light_width         = t%pheno_light_width(ipft)
         p%light_window        = t%pheno_light_window(ipft)
      end associate
   end subroutine flatten_pheno_params

end module meds_vegetation_dynamics
