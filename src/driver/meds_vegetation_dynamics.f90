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
   use meds_config,               only : meds_config_t, growth_window_steps
   use meds_allometry,            only : size2leaf_carbon
   use meds_core_state_types,      only : carbon_flux_block, GROWTH_AVG_UNSET
   use meds_core_interface, only : site_t, apply_growth, mortality_step, apply_recruitment,&
                                         apply_patch_disturbance, new_fuse_cohorts,              &
                                         terminate_cohorts, split_cohorts, new_fuse_patches,     &
                                         terminate_patches, sort_cohorts, sort_patches,          &
                                         update_overtopping_lai
   use meds_plant_vital_rates,    only : carbon_growth_rate, camac_mortality, min_cohort_carbon, &
                                         recruitment_contribution
   use meds_plant_interface,      only : get_plant_flux_slow, growth_respiration,                &
                                         carbon_env_t, carbon_demand_t, carbon_npp_t, PHEN_ON
   implicit none
   private

   public :: vegetation_dynamics

   !----- Wood is the residual carbon sink: a demand large enough to take all remaining NPP. --!
   real(wp), parameter :: WOOD_DEMAND_BIG   = 1.0e6_wp    !< [kgC/plant]
   real(wp), parameter :: STUB_TISSUE_TEMP  = 298.15_wp   !< [K] 25 degC (no met forcing yet)
   !----- Patch structural dynamics (fusion/fission/disturbance) are an ANNUAL process, so       !
   !       disturbance integrates its yearly rate over this interval, not the per-step dt.        !
   real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp   !< [yr]

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the carbon NPP, compute the carbon !
   ! vital rates via the plant kernels, and sequence the demography apply-primitives + cadence. !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, is_new_month, is_new_year)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable   :: mortality(:), recruitment(:,:), npp_repro(:)
      type(carbon_flux_block) :: npp
      logical                 :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse
      integer(ik)             :: ip, n_window

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

      !----- Growth (npp -> pools -> geometry + the moving-average ring buffer), then mortality. -!
      call apply_growth(site%cohort, npp, cfg%dt_years, n_window, site%growth_hist_pos)
      call mortality_step(site%cohort%n, site%cohort%nplant, mortality, cfg%dt_years, cfg%negligible_nplant)

      !----- Re-sort every step: growth changed heights, so re-establish the tallest-first order  !
      !      the overtopping-LAI sweep + the patch-light profiles depend on (fuse/fiss also sorts, !
      !      but only on the monthly/annual cadence). -------------------------------------------!
      call sort_cohorts(site)

      !----- Patch ageing (every step). ----------------------------------------------------!
      do ip = 1_ik, site%patch%n
         site%patch%age(ip) = site%patch%age(ip) + cfg%dt_years
      end do

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
      real(wp)    :: leaf_target, a_carbon, gross_gpp, resp_maint
      type(carbon_env_t)    :: env
      type(carbon_demand_t) :: demand
      type(carbon_npp_t)    :: out

      n = site%cohort%n
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
            env%tissue_temp      = STUB_TISSUE_TEMP
            env%dt_yr            = dt_yr
            env%phenology_status = merge(cohort%phenology_status(j), PHEN_ON, cfg%phenology_on)
            call get_plant_flux_slow(env, cfg, pf, demand, out)
            npp%leaf(j)          = out%leaf
            npp%fineroot(j)      = out%fineroot
            npp%wood(j)          = out%wood
            npp%nonstructural(j) = out%nonstructural
            npp_repro(j)         = out%repro
         end do
      end associate
   end subroutine carbon_growth

end module meds_vegetation_dynamics
