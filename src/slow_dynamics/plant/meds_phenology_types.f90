! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_phenology_types -- the derived types of the leaf-phenology kernel (meds_phenology): !
! its daily drivers, per-PFT traits, per-cohort state and outputs, and the cue bits        !
! (docs/science/plant_phenology.md). They are kept apart from the fast leaf types, so the  !
! slow plant kernels and the fast ones do not depend on each other. Data only.             !
!==========================================================================================!
module meds_phenology_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pheno_env_t, pheno_args_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_LIGHT, CUE_WATER, CUE_ALL

   !----- Cue bits: flush_cue_mask and shed_cue_mask are sums of them, selecting the cues of !
   !      each side (docs/science/plant_phenology.md §2-§3). ------------------------------!
   integer(ik), parameter :: CUE_NONE  = 0_ik   !< no cues (always flushing / no senescence)
   integer(ik), parameter :: CUE_TEMP  = 1_ik   !< temperature: warmth sum (flush), cold sum (shed)
   integer(ik), parameter :: CUE_LIGHT = 2_ik   !< light: running-mean hours a day of PAR above par_min
   integer(ik), parameter :: CUE_WATER = 4_ik   !< water: predawn leaf psi summed against the TLP
   integer(ik), parameter :: CUE_ALL   = CUE_TEMP + CUE_LIGHT + CUE_WATER

   !----- Daily environmental drivers (read-only). ------------------------------------------!
   type :: pheno_env_t
      real(wp)    :: temp_day         = 0.0_wp    !< [K]   daily-mean air temperature
      real(wp)    :: par_hours        = 12.0_wp   !< [h/day] hours with PAR at the cohort's top above par_min
      real(wp)    :: predawn_leaf_psi = 0.0_wp    !< [MPa, <=0] predawn (daily-max) leaf water potential
      integer(ik) :: doy              = 1_ik      !< [-]   day of year
      logical     :: hemis_north      = .true.    !< northern hemisphere (season windows)
   end type pheno_env_t

   !----- The prognostic memory: two smoothed tendencies + the cue accumulators. -----------!
   type :: pheno_state_t
      real(wp) :: leaf_flush_tendency = 1.0_wp   !< [-] smoothed flush signal; a cohort is born flushing
      real(wp) :: leaf_shed_tendency  = 0.0_wp   !< [-] smoothed senescence signal
      real(wp) :: growing_degree_days = 0.0_wp   !< [K day] warmth sum since midwinter     (TEMP flush)
      real(wp) :: cold_degree_days    = 0.0_wp   !< [K day] cold sum since midsummer       (TEMP shed)
      real(wp) :: wet_psi_sum         = 0.0_wp   !< [MPa day] sum of psi above the TLP     (WATER flush)
      real(wp) :: dry_psi_sum         = 0.0_wp   !< [MPa day] sum of psi below the TLP     (WATER shed)
      real(wp) :: light_hours_mean    = -1.0_wp  !< [h/day] running-mean par_hours (LIGHT); < 0: unset, the first day sets it
   end type pheno_state_t

   !----- Per-PFT traits, filled by the driver from cfg%pft. Each switch is sigma(s (x - x*)), !
   !      with a centre x* and a signed sharpness s (plant_phenology.md eq. (2)). The defaults !
   !      are a PFT with no [phenology] section: no cues, so only the rates and timescales act. !
   type :: pheno_args_t
      integer(ik) :: flush_cue_mask      = CUE_NONE         !< cues of the flush side
      integer(ik) :: shed_cue_mask       = CUE_NONE         !< cues of the shed side
      real(wp)    :: flush_cue_timescale = 5.0_wp           !< [day]  smoothing of the flush tendency
      real(wp)    :: shed_cue_timescale  = 5.0_wp           !< [day]  smoothing of the shed tendency
      real(wp)    :: flush_rate_max      = 0.06667_wp       !< [1/day] leaf growth at full tendency
      real(wp)    :: shed_rate_max       = 0.05_wp          !< [1/day] senescence at full tendency
      !----- Temperature: warmth above flush_base_temp from midwinter, cold below          !
      !      shed_base_temp from midsummer.                                                !
      real(wp)    :: flush_base_temp       = 278.15_wp      !< [K]
      real(wp)    :: flush_degree_days     = 100.0_wp       !< [K day] warmth requirement
      real(wp)    :: flush_temp_sharpness  = 0.04_wp        !< [1/(K day)]
      real(wp)    :: shed_base_temp        = 290.15_wp      !< [K]
      real(wp)    :: shed_degree_days      = 50.0_wp        !< [K day] cold requirement
      real(wp)    :: shed_temp_sharpness   = 0.1_wp         !< [1/(K day)]
      !----- Light: the hours a day PAR at the cohort's top exceeds par_min, as a running mean. !
      real(wp)    :: par_min               = 5.0_wp         !< [umol/m2/s] light that counts as an hour of light
      real(wp)    :: flush_light_hours     = 12.0_wp        !< [h/day]
      real(wp)    :: flush_light_sharpness = 1.0_wp         !< [1/h] > 0: long or bright days permit flushing
      real(wp)    :: shed_light_hours      = 11.0_wp        !< [h/day]
      real(wp)    :: shed_light_sharpness  = -1.0_wp        !< [1/h] < 0: short or dim days trigger senescence; > 0: bright days
      real(wp)    :: light_window          = 10.0_wp        !< [day] running-mean window
      !----- Water: predawn leaf psi summed against the turgor-loss point. -----------------!
      real(wp)    :: leaf_psi_tlp          = -2.0_wp        !< [MPa] turgor-loss point
      real(wp)    :: flush_water_sum       = 10.0_wp        !< [MPa day] wet sum that permits flushing
      real(wp)    :: flush_water_sharpness = 0.5_wp         !< [1/(MPa day)]
      real(wp)    :: shed_water_sum        = 10.0_wp        !< [MPa day] dry sum that triggers shedding
      real(wp)    :: shed_water_sharpness  = 0.5_wp         !< [1/(MPa day)]
   end type pheno_args_t

   !----- Outputs: the potential relative rates the tendencies allow. ----------------------!
   type :: pheno_out_t
      real(wp) :: leaf_flush_potential = 0.0_wp   !< [1/day] flush_rate_max * leaf_flush_tendency
      real(wp) :: leaf_shed_potential  = 0.0_wp   !< [1/day] shed_rate_max  * leaf_shed_tendency
   end type pheno_out_t

end module meds_phenology_types
