! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_phenology_types -- the derived types of the PHENOLOGY kernel.                            !
!                                                                                          !
! Split out of meds_plant_types when the plant library split by timescale: phenology runs on !
! the DAILY tier, everything else that module held (leaf gas exchange, hydraulics, non-leaf   !
! maintenance respiration) runs sub-daily. Keeping one types module would have made the slow   !
! plant kernels depend on the fast ones and back again -- a genuine cycle, not a style point.  !
!                                                                                          !
! Pure DATA: no methods, no hidden state.                                                    !
!==========================================================================================!
module meds_phenology_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_LIGHT, CUE_WATER, CUE_ALL
   public :: LIGHT_DAYLENGTH, LIGHT_RADIATION

   !=======================================================================================!
   !     PHENOLOGY -- a SIGNAL kernel: daily cues + per-PFT traits -> two smoothed tendencies  !
   !     in [0,1] (leaf_flush_tendency, leaf_shed_tendency). It touches no carbon; the carbon  !
   !     layer turns the tendencies into leaf growth and leaf loss (leaf_turnover_step).       !
   !     docs/science/plant_phenology.md.                                                     !
   !=======================================================================================!
   !----- Cue-enable bits. flush_cue_mask and shed_cue_mask select the cues of each side: the !
   !      flush signal is the PRODUCT of its cues' switches, the shed signal the larger of the !
   !      seasonal trigger (TEMP x LIGHT) and the water trigger.                               !
   integer(ik), parameter :: CUE_NONE  = 0_ik    !< no cues (always flushing / no senescence)
   integer(ik), parameter :: CUE_TEMP  = 1_ik    !< temperature: warmth sum (flush), cold sum (shed)
   integer(ik), parameter :: CUE_LIGHT = 2_ik    !< light: day length or running-mean radiation
   integer(ik), parameter :: CUE_WATER = 4_ik    !< water: predawn leaf psi summed against the TLP
   integer(ik), parameter :: CUE_ALL   = CUE_TEMP + CUE_LIGHT + CUE_WATER

   !----- Which variable the light cue reads (per PFT). -------------------------------------!
   integer(ik), parameter :: LIGHT_DAYLENGTH = 1_ik   !< day length [h]
   integer(ik), parameter :: LIGHT_RADIATION = 2_ik   !< running-mean incident shortwave [W/m2]

   !----- Daily environmental drivers (read-only). ------------------------------------------!
   type :: pheno_env_t
      real(wp)    :: temp_day         = 0.0_wp    !< [K]   daily-mean air temperature
      real(wp)    :: daylength        = 12.0_wp   !< [h]   day length
      real(wp)    :: rad              = 0.0_wp    !< [W/m2] daily-mean incident shortwave
      real(wp)    :: predawn_leaf_psi = 0.0_wp    !< [MPa, <=0] predawn (daily-max) leaf water potential
      integer(ik) :: doy              = 1_ik      !< [-]   day of year
      logical     :: hemis_north      = .true.    !< northern hemisphere (season windows)
   end type pheno_env_t

   !----- The prognostic memory: two smoothed tendencies + the cue accumulators. -----------!
   type :: pheno_state_t
      real(wp) :: leaf_flush_tendency = 1.0_wp   !< [-] smoothed flush signal; born flushing (the
                                                 !<     always-flushing, never-senescing fixed point)
      real(wp) :: leaf_shed_tendency  = 0.0_wp   !< [-] smoothed senescence signal
      real(wp) :: growing_degree_days = 0.0_wp   !< [K day] warmth sum since midwinter     (TEMP flush)
      real(wp) :: cold_degree_days    = 0.0_wp   !< [K day] cold sum since midsummer       (TEMP shed)
      real(wp) :: wet_psi_sum         = 0.0_wp   !< [MPa day] sum of psi above the TLP     (WATER flush)
      real(wp) :: dry_psi_sum         = 0.0_wp   !< [MPa day] sum of psi below the TLP     (WATER shed)
      real(wp) :: shortwave_mean      = 0.0_wp   !< [W/m2] running-mean shortwave  (LIGHT, radiation)
   end type pheno_state_t

   !----- Flat per-PFT trait set (filled by the driver from cfg%pft). ----------------------!
   !      Every switch is sigma(s (x - x*)): a centre x* and a signed sharpness s [1/units of x];
   !      its 0.12-0.88 transition spans x* +- 2/|s|. The defaults below are the ones a PFT     !
   !      without a [phenology] section keeps (both masks CUE_NONE: always flushing, no       !
   !      senescence), so only the rates and the timescales act.                              !
   type :: pheno_params_t
      integer(ik) :: flush_cue_mask      = CUE_NONE         !< cues of the flush side
      integer(ik) :: shed_cue_mask       = CUE_NONE         !< cues of the shed side
      real(wp)    :: flush_cue_timescale = 5.0_wp           !< [day]  smoothing of the flush tendency
      real(wp)    :: shed_cue_timescale  = 5.0_wp           !< [day]  smoothing of the shed tendency
      real(wp)    :: flush_rate_max      = 0.06667_wp       !< [1/day] leaf growth at full tendency
      real(wp)    :: shed_rate_max       = 0.05_wp          !< [1/day] senescence at full tendency
      !----- Temperature: warmth above flush_base_temp from midwinter, cold below            !
      !      shed_base_temp from midsummer.                                                   !
      real(wp)    :: flush_base_temp       = 278.15_wp      !< [K]
      real(wp)    :: flush_degree_days     = 100.0_wp       !< [K day] warmth requirement
      real(wp)    :: flush_temp_sharpness  = 0.04_wp        !< [1/(K day)]
      real(wp)    :: shed_base_temp        = 290.15_wp      !< [K]
      real(wp)    :: shed_degree_days      = 50.0_wp        !< [K day] cold requirement
      real(wp)    :: shed_temp_sharpness   = 0.1_wp         !< [1/(K day)]
      !----- Light: day length [h] or running-mean shortwave [W/m2], per light_variable. ---!
      integer(ik) :: light_variable        = LIGHT_DAYLENGTH
      real(wp)    :: flush_light_threshold = 12.0_wp        !< [h | W/m2]
      real(wp)    :: flush_light_sharpness = 1.0_wp         !< [1/h | m2/W] > 0: more light permits
      real(wp)    :: shed_light_threshold  = 11.0_wp        !< [h | W/m2]
      real(wp)    :: shed_light_sharpness  = -1.0_wp        !< < 0: short days trigger senescence;
                                                            !< > 0: bright light does (leaf exchange)
      real(wp)    :: light_window          = 10.0_wp        !< [day] running mean (LIGHT_RADIATION)
      !----- Water: predawn leaf psi summed against the turgor-loss point. -----------------!
      real(wp)    :: leaf_psi_tlp          = -2.0_wp        !< [MPa] turgor-loss point
      real(wp)    :: flush_water_sum       = 10.0_wp        !< [MPa day] wet sum that permits flushing
      real(wp)    :: flush_water_sharpness = 0.5_wp         !< [1/(MPa day)]
      real(wp)    :: shed_water_sum        = 10.0_wp        !< [MPa day] dry sum that triggers shedding
      real(wp)    :: shed_water_sharpness  = 0.5_wp         !< [1/(MPa day)]
   end type pheno_params_t

   !----- Outputs: the potential relative rates the tendencies allow. ----------------------!
   type :: pheno_out_t
      real(wp) :: leaf_flush_potential = 0.0_wp   !< [1/day] flush_rate_max * leaf_flush_tendency
      real(wp) :: leaf_shed_potential  = 0.0_wp   !< [1/day] shed_rate_max  * leaf_shed_tendency
   end type pheno_out_t

end module meds_phenology_types
