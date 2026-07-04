!==========================================================================================!
! meds_pheno_types -- data structures of the leaf-phenology interface seam.                 !
!                                                                                          !
! Phenology is a PURE DIRECTIONAL signal generator: given plant traits and environmental      !
! conditions it predicts the phenological STATUS -- a direction of leaf-display change        !
! (ON = seek growth, OFF = drop leaves, DORMANT = hold) -- and nothing else. It does NOT       !
! compute leaf growth, leaf drop, a target leaf level, or carbon; the displayed leaf fraction   !
! EMERGES downstream (a separate leaf-dynamics module) from these states + the leaf growth/     !
! shed rates. Four public types describe the calculation:                                       !
!   * pheno_env_t    -- the EXPLICIT daily environmental drivers (read-only).                    !
!   * pheno_params_t -- a flat, self-contained per-PFT trait set (cue-enable mask + selectors +   !
!                       per-cue parameters), so the engine never reaches back into meds_config.    !
!   * pheno_state_t  -- the prognostic cue ACCUMULATORS (phenological memory: degree-day sums,      !
!                       chilling count, running-mean water, leaf-psi day counters). NOT leaf mass.   !
!   * pheno_out_t    -- the phenological status {ON,OFF,DORMANT} + which cue governed it.             !
!                                                                                          !
! The only per-PFT "strategy" is the cue-enable MASK; every ED2 habit is a mask + parameter       !
! special-case (evergreen = CUE_NONE). These types are pure DATA -- no methods, no hidden state.   !
!==========================================================================================!
module meds_pheno_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO
   public :: PHEN_ON, PHEN_DORMANT, PHEN_OFF

   !----- Cue-enable mask bits (per-PFT, OR-combined). The ONLY strategy selector. ---------!
   integer(ik), parameter :: CUE_NONE  = 0_ik   !< no cues => evergreen (perpetually ON)
   integer(ik), parameter :: CUE_TEMP  = 1_ik   !< temperature (GDD/CDD cold-deciduous)
   integer(ik), parameter :: CUE_WATER = 2_ik   !< soil-water drought (running-mean available water)
   integer(ik), parameter :: CUE_HYDRO = 4_ik   !< leaf water potential (hydraulic)
   integer(ik), parameter :: CUE_PHOTO = 8_ik   !< photoperiod (daylength gate on the temperature cue)

   !----- Phenological status: the direction of leaf-display change. -----------------------!
   integer(ik), parameter :: PHEN_OFF     = -1_ik  !< unfavorable -- actively dropping leaves
   integer(ik), parameter :: PHEN_DORMANT =  0_ik  !< neutral deadband -- hold the current state
   integer(ik), parameter :: PHEN_ON      =  1_ik  !< favorable -- actively seeking leaf growth

   !----- Raw daily environmental drivers (read-only; all a caller boundary condition). ----!
   type :: pheno_env_t
      real(wp)    :: temp_day    = 0.0_wp    !< [K]   daily-mean air/canopy temperature (thermal sums)
      real(wp)    :: soil_temp   = 0.0_wp    !< [K]   shallow-layer soil temperature (cold-drop trigger)
      real(wp)    :: avail_water = 0.0_wp    !< [-] fraction OR [MPa] soil-water potential (CUE_WATER)
      real(wp)    :: psi_leaf    = 0.0_wp    !< [MPa, <=0] leaf water potential (CUE_HYDRO)
      real(wp)    :: daylength   = 12.0_wp   !< [h]   photoperiod (CUE_PHOTO; caller-supplied)
      integer(ik) :: doy         = 1_ik      !< [-]   day-of-year (thermal-sum season gating)
      logical     :: hemis_north = .true.    !< northern hemisphere (season gating)
   end type pheno_env_t

   !----- Cue accumulators: the prognostic phenological MEMORY (not leaf-mass state). ------!
   type :: pheno_state_t
      real(wp) :: gdd           = 0.0_wp   !< [K day] growing-degree-day sum          (CUE_TEMP)
      real(wp) :: chill         = 0.0_wp   !< [day]   chilling-day count              (CUE_TEMP)
      real(wp) :: water_avg     = 0.0_wp   !< [-]|[MPa] running-mean available water  (CUE_WATER)
      real(wp) :: low_psi_days  = 0.0_wp   !< [day]   consecutive days psi_leaf < psi_tlp   (CUE_HYDRO)
      real(wp) :: high_psi_days = 0.0_wp   !< [day]   consecutive days psi_leaf >= 0.5 psi_tlp (CUE_HYDRO)
   end type pheno_state_t

   !----- Flat per-PFT trait set (self-contained; filled by the seam from cfg%pft). --------!
   type :: pheno_params_t
      !----- Selector: the cue mask, the shared logistic sharpness, the on/off band. --------!
      integer(ik) :: cue_mask           = CUE_NONE   !< OR of CUE_* (evergreen = CUE_NONE)
      real(wp)    :: cue_sharpness       = 2.0_wp     !< [-] dimensionless logistic slope (large => ED2-sharp)
      real(wp)    :: phen_on_threshold   = 0.6_wp     !< favorability above which status = ON
      real(wp)    :: phen_off_threshold  = 0.4_wp     !< favorability below which status = OFF (< on => DORMANT band)
      !----- Thermal (CUE_TEMP): GDD flush threshold a+b*exp(c*chill) + autumn cold drop. ---!
      real(wp)    :: gdd_base_temp       = 278.15_wp  !< [K] GDD accumulation base (5 degC)
      real(wp)    :: chill_base_temp     = 278.15_wp  !< [K] chilling-day base
      real(wp)    :: phen_a              = -68.0_wp   !< [K day] GDD threshold intercept  (Botta 2000)
      real(wp)    :: phen_b              = 638.0_wp   !< [K day] GDD threshold amplitude
      real(wp)    :: phen_c              = -0.01_wp   !< [1/day] chilling exponent (more chill => lower GDD need)
      real(wp)    :: cold_drop_daylength = 10.9_wp    !< [h] autumn short-day drop trigger  (White 1997)
      real(wp)    :: cold_drop_soiltemp1 = 284.3_wp   !< [K] cool-soil drop (with short days)
      real(wp)    :: cold_drop_soiltemp2 = 275.15_wp  !< [K] very-cold-soil drop (unconditional)
      !----- Water (CUE_WATER): running-mean ramp between off/on thresholds. ----------------!
      logical     :: water_use_potential = .false.    !< .false.: moisture fraction; .true.: soil-psi [MPa]
      real(wp)    :: water_off_threshold = 0.2_wp      !< available water at which favorability = 0
      real(wp)    :: water_on_threshold  = 0.5_wp      !< available water at which favorability = 1 (> off)
      real(wp)    :: water_window        = 10.0_wp     !< [day] running-mean window
      !----- Hydraulic (CUE_HYDRO): leaf-psi consecutive-day counters vs the TLP. -----------!
      real(wp)    :: leaf_psi_tlp        = -2.0_wp     !< [MPa] turgor-loss point (Xu 2016)
      real(wp)    :: low_psi_threshold   = 10.0_wp     !< [day] dry days to full unfavorable
      real(wp)    :: high_psi_threshold  = 10.0_wp     !< [day] wet days to full favorable
      !----- Photoperiod (CUE_PHOTO): daylength logistic gate. ------------------------------!
      real(wp)    :: photo_crit          = 11.0_wp     !< [h] critical daylength
      real(wp)    :: photo_slope         = 2.0_wp      !< [1/h] daylength logistic slope
   end type pheno_params_t

   !----- Outputs: the phenological status + the governing (most-limiting) cue. ------------!
   type :: pheno_out_t
      integer(ik) :: phenology_status = PHEN_DORMANT   !< PHEN_ON | PHEN_OFF | PHEN_DORMANT
      integer(ik) :: cue_limiting     = CUE_NONE       !< the CUE_* bit with the lowest favorability
   end type pheno_out_t

end module meds_pheno_types
