!==========================================================================================!
! meds_vegetation_dynamics -- MEDS's slow-loop VEGETATION-DYNAMICS driver (the analogue of      !
! ED2's veg_dynamics_driver). It lives in src/driver/ (compiled into meds_aux), the only layer  !
! above BOTH the plant-ecophysiology library and the demography engine, and it is where they    !
! meet: it ASSEMBLES the per-cohort demographic rates from the stateless plant rate kernels and  !
! APPLIES them through the demography engine's data seam (update_demography). Neither library     !
! depends on the other; this driver weaves them.                                                 !
!                                                                                          !
! Rate assembly (empirical path): the height-sorted per-patch sweep that gives each cohort its   !
! overtopping LAI is a STAND-structure operation (it reads the cohort/patch CSR), so it lives     !
! here; the per-cohort LAWS (growth/mortality/recruitment) are stateless kernels in              !
! meds_plant_vital_rates. A mechanistic carbon path will be selected here later (growth_source),  !
! reusing the same apply-through-update_demography tail.                                          !
!                                                                                          !
! NB: distinct from meds_demography_dynamics -- that is the ENGINE (the growth/mortality/         !
! disturbance kernels that MUTATE state); this is the DRIVER that assembles rates and drives it.  !
!==========================================================================================!
module meds_vegetation_dynamics
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t, update_demography
   use meds_plant_interface,      only : growth_rate_empirical, mortality_rate,               &
                                         recruitment_rate, min_cohort_carbon
   implicit none
   private

   !----- vegetation_dynamics is the driver; empirical_vital_rates is the rate PROVIDER it     !
   !       calls (also exposed so it can be unit-tested and, later, A/B-compared with a          !
   !       mechanistic carbon provider that produces the same three arrays).                     !
   public :: vegetation_dynamics, empirical_vital_rates

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the vegetation dynamics for one step: assemble the demographic rates for the    !
   ! current stand, fold the calendar cadence + switches into the structural triggers, and    !
   ! apply everything through the demography data seam. (A future master step would call this  !
   ! on the slow cadence alongside the fast-loop biophysics processes.)                       !
   !---------------------------------------------------------------------------------------!
   subroutine vegetation_dynamics(site, cfg, is_new_month, is_new_year)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      real(wp), allocatable :: growth(:), mortality(:), recruitment(:,:)
      logical               :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse

      !----- 1. Assemble demographic rates for the current stand (empirical path for now; a   !
      !         mechanistic carbon path will be selected here later via growth_source). -------!
      call empirical_vital_rates(site, cfg, growth, mortality, recruitment)

      !----- 2. Fold the calendar cadence + demography on/off + fiss/fuse switches into the   !
      !         structural triggers the engine consumes. -------------------------------------!
      do_cohort_fissfuse   = is_new_month .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
      do_patch_disturbance = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_disturbance
      do_patch_fissfuse    = is_new_year  .and. cfg%demography_on .and. cfg%do_patch_fissfuse

      !----- 3. Apply the rates through the demography engine's data seam. --------------------!
      call update_demography(site, growth, mortality, recruitment, cfg, cfg%dt_years,         &
                             do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse)
   end subroutine vegetation_dynamics

   !---------------------------------------------------------------------------------------!
   ! Evaluate the EMPIRICAL (phenomenological) vital rates for the current site_t in ONE       !
   ! height-sorted, per-patch overtopping-LAI sweep: growth and mortality per cohort, and --   !
   ! in the same pass -- the reproduction part of recruitment, reduced into the per-(PFT,patch) !
   ! array. The per-cohort laws are the stateless kernels in meds_plant_vital_rates; THIS       !
   ! routine owns the stand sweep and the reduction (both need the cohort/patch structure).    !
   ! Produces the three arrays update_demography consumes; a mechanistic provider produces the  !
   ! same arrays with no engine change.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine empirical_vital_rates(site, cfg, growth, mortality, recruitment)
      type(site_t),          intent(in)  :: site
      type(meds_config_t),   intent(in)  :: cfg
      real(wp), allocatable, intent(out) :: growth(:)         !< [cm/yr]  per cohort
      real(wp), allocatable, intent(out) :: mortality(:)      !< [1/yr]   per cohort, total
      real(wp), allocatable, intent(out) :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: n, np, npft, i, j, k, ip, pf, i0, i1
      real(wp)    :: cum, over_lai, layer_lai
      real(wp), allocatable :: carbon_min(:)     ! [kgC/plant] carbon of one min-size cohort, per PFT

      n    = site%cohort%n
      np   = site%patch%n
      npft = cfg%pft%n
      allocate(growth(n), mortality(n), recruitment(npft, max(np, 1_ik)))
      recruitment = 0.0_wp

      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         !----- Carbon of one minimum-size cohort (the recruit "unit"), per PFT. -----------!
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

         !----- Per-patch top-down sweep over the overtopping LAI. -------------------------!
         do ip = 1_ik, np
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            i   = i0
            do while (i <= i1)
               !----- Equal-height cohorts are co-dominant (share the overtopping LAI). -----!
               k = i
               do while (k < i1)
                  if (cohort%height(k + 1_ik) /= cohort%height(i)) exit
                  k = k + 1_ik
               end do
               over_lai  = cum
               layer_lai = 0.0_wp
               do j = i, k
                  pf = cohort%pft(j)
                  !----- Growth: empirical law (intrinsic x competition x repro allocation). --!
                  growth(j) = growth_rate_empirical(cohort%dbh(j), over_lai, cohort%height(j),  &
                                 pft%growth_dbh_max(pf), pft%growth_dbh_slope(pf),              &
                                 pft%growth_dbh_cap(pf), pft%growth_lai_slope(pf),              &
                                 pft%min_reproduction_height, pft%reproduction_investment_fraction(pf))

                  !----- Mortality (Camac additive hazard) from the tracked running-mean       !
                  !      growth; a not-yet-seeded cohort uses its instantaneous growth.        !
                  mortality(j) = mortality_rate(cohort%growth_avg(j), growth(j),               &
                                 pft%mort_gamma(pf), pft%mort_alpha(pf), pft%mort_beta(pf))

                  !----- Reproduction flux -> recruits (zero below the maturity height). -------!
                  recruitment(pf, ip) = recruitment(pf, ip)                                    &
                       + recruitment_rate(cohort%dbh(j), cohort%height(j), over_lai,           &
                                          cohort%agb(j), cohort%nplant(j), cohort%p_hgt_max(j),&
                                          cohort%p_wood_density(j), pft%growth_dbh_max(pf),    &
                                          pft%growth_dbh_slope(pf), pft%growth_dbh_cap(pf),    &
                                          pft%growth_lai_slope(pf), pft%min_reproduction_height,&
                                          pft%reproduction_investment_fraction(pf),            &
                                          pft%repro_carbon_efficiency(pf), carbon_min(pf))

                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine empirical_vital_rates

end module meds_vegetation_dynamics
