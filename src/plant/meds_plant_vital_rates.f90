!==========================================================================================!
! meds_plant_vital_rates -- the per-INDIVIDUAL vital-rate KERNELS, as pure scalar functions.    !
! These are the physiological/phenomenological laws that turn one plant's carbon + tracked        !
! growth into its demographic rates: a carbon-driven diameter growth rate, the Camac (2018)       !
! additive mortality hazard, and the reproduction-carbon -> recruit conversion. They are          !
! `elemental pure`, take ONLY scalars + PFT traits (NO site_t, NO cohort SoA), and use the         !
! shared allometry -- so they live in the stateless plant library and are CALLED BY THE DRIVER     !
! (meds_vegetation_dynamics), never by the demography engine. That keeps `demography _|_ plant`:   !
! the engine only APPLIES the resulting rate arrays.                                               !
!                                                                                          !
! (Lifted from the deleted meds_demography_rates carbon path. The empirical growth/recruitment     !
! LAWS are NOT here -- they moved to the Python example. Camac mortality stays in Fortran.)         !
!==========================================================================================!
module meds_plant_vital_rates
   use meds_kinds,     only : wp
   use meds_allometry, only : wood_to_dbh
   implicit none
   private

   public :: npp_to_growth, camac_mortality, npp_to_recruitment

contains

   !----- Prospective carbon diameter-growth rate [cm/yr]: the dbh the cohort WOULD reach if    !
   !       this step's wood NPP were added to its size-anchor wood_carbon (floored >= 0), minus   !
   !       its current dbh, over the step. Feeds both the applied growth and the mortality hazard.!
   elemental pure function npp_to_growth(wood_carbon, npp_wood, dbh, dt_yr,                     &
                                         rho, hgt_max, aboveground_frac) result(dbh_rate)
      real(wp), intent(in) :: wood_carbon, npp_wood, dbh, dt_yr, rho, hgt_max, aboveground_frac
      real(wp)             :: dbh_rate, dbh_new
      dbh_new  = wood_to_dbh(max(wood_carbon + npp_wood, 0.0_wp), rho, hgt_max, aboveground_frac)
      dbh_rate = (dbh_new - dbh) / dt_yr
   end function npp_to_growth

   !----- Camac et al. (2018) additive mortality hazard [1/yr] on the cohort's EFFECTIVE growth:  !
   !       its tracked running-mean growth (growth_avg) once seeded, else the instantaneous rate.  !
   !       `avg_is_set` is the caller's test of the moving-average sentinel (kept out of this      !
   !       state-free kernel so plant _|_ state holds).                                             !
   elemental pure function camac_mortality(growth_avg, growth_inst, avg_is_set,                 &
                                           mort_gamma, mort_alpha, mort_beta) result(mort)
      real(wp), intent(in) :: growth_avg, growth_inst, mort_gamma, mort_alpha, mort_beta
      logical,  intent(in) :: avg_is_set
      real(wp)             :: mort, g_eff
      g_eff = merge(growth_avg, growth_inst, avg_is_set)
      mort  = mort_gamma + mort_alpha * exp(-mort_beta * g_eff)
   end function camac_mortality

   !----- One cohort's reproduction-carbon recruit contribution [plant/m2/yr]: the per-plant       !
   !       reproduction NPP over the step, as an annual rate, times the density and the            !
   !       establishment efficiency, over the min-size recruit carbon (min_cohort_carbon, from      !
   !       meds_allometry). Caller reduces these into the (PFT, patch) array and gates by include_pft.!
   elemental pure function npp_to_recruitment(nplant, npp_repro, dt_yr,                         &
                                              repro_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: nplant, npp_repro, dt_yr, repro_carbon_efficiency, carbon_min
      real(wp)             :: rec
      rec = nplant * (npp_repro / dt_yr) * repro_carbon_efficiency / carbon_min
   end function npp_to_recruitment

end module meds_plant_vital_rates
