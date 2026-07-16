!==========================================================================================!
! meds_vegetation_dynamics -- MEDS's slow-loop VEGETATION-DYNAMICS driver (the analogue of      !
! ED2's veg_dynamics_driver). It lives in src/driver/ (compiled into meds_aux) and drives the   !
! slow loop: it assembles the per-cohort demographic rates -- EMPIRICAL (meds_demography_rates)  !
! or CARBON-sourced (the plant slow-flux seam) per cfg%growth_source -- folds the calendar        !
! cadence into the structural triggers, and applies everything through update_demography.         !
!                                                                                          !
! The carbon path is the ONLY place plant and demography meet: this driver computes the           !
! allometric demands + stub GPP, calls the stateless plant carbon seam get_plant_flux_slow to     !
! get the per-pool NPP, and hands it to the (plant-free) demography engine. `demography ⊥ plant`.  !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t, GS_CARBON
   use meds_allometry,            only : size2leaf_carbon
   use meds_demography_interface, only : site_t, update_demography
   use meds_ecosystem_state,     only : carbon_flux_block
   use meds_demography_rates,     only : empirical_vital_rates, carbon_vital_rates
   use meds_plant_interface,      only : get_plant_flux_slow, growth_respiration,               &
                                         carbon_env_t, carbon_demand_t, carbon_npp_t, PHEN_ON
   implicit none
   private

   public :: vegetation_dynamics

   !----- Wood is the residual carbon sink: a demand large enough to take all remaining NPP. --!
   real(wp), parameter :: WOOD_DEMAND_BIG   = 1.0e6_wp    !< [kgC/plant]
   real(wp), parameter :: STUB_TISSUE_TEMP  = 298.15_wp   !< [K] 25 degC (no met forcing yet)

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the demographic rates (empirical  !
   ! or carbon), fold the calendar cadence + switches into the structural triggers, and apply  !
   ! everything through the demography data seam (which dispatches growth_step | apply_carbon_npp).!
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, is_new_month, is_new_year)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable   :: growth(:), mortality(:), recruitment(:,:), npp_repro(:)
      type(carbon_flux_block) :: npp
      logical                 :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      !----- 1. Assemble the demographic rates: empirical OR carbon-sourced. ----------------!
      if (cfg%growth_source == GS_CARBON) then
         call carbon_growth(site, cfg, cfg%dt_years, npp, npp_repro)
         call carbon_vital_rates(site, cfg, npp%wood, npp_repro, cfg%dt_years, mortality, recruitment)
         allocate(growth(0))                              ! growth[] unused in carbon mode
      else
         call empirical_vital_rates(site, cfg, growth, mortality, recruitment)
      end if

      !----- 2. Fold the calendar cadence + demography on/off + fiss/fuse switches. ---------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- 3. Apply the rates through the demography engine's data seam. --------------------!
      call update_demography(site, growth, npp, mortality, recruitment, cfg, cfg%dt_years,     &
                             cfg%growth_source, do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse)
   end subroutine vegetation_dynamics

   !---------------------------------------------------------------------------------------!
   ! CARBON growth assembler (the carbon twin of empirical_vital_rates). Per cohort it builds  !
   ! the allometric demands (target - pool) and the net carbon (stub GPP per leaf area minus    !
   ! growth respiration), calls the plant slow-flux seam get_plant_flux_slow, and returns the   !
   ! per-pool NPP block + the reproduction carbon. This is the ONLY place the driver calls plant.!
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
            !----- Net carbon over the slow step. Fast-on: accumulated GROSS GPP minus the         !
            !      accumulated autotrophic MAINTENANCE respiration (leaf Rd + stem + fine root),    !
            !      all [kgC/plant] from the fast loop (BUG3 fix -- maintenance was previously        !
            !      dropped). Fast-off stub: gpp_ref per leaf area, no mechanistic maintenance term.  !
            !      GROWTH respiration is then deducted ONCE from the maintenance-adjusted carbon.    !
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
            !----- Leaf-flush gate: obey the phenology driver's per-cohort status when phenology is  !
            !      ON; otherwise hold the legacy always-ON behaviour (bit-identical to before).  -----!
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
