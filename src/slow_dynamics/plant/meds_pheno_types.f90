!==========================================================================================!
! meds_pheno_types -- the derived types of the PHENOLOGY kernel.                            !
!                                                                                          !
! Split out of meds_plant_types when the plant library split by timescale: phenology runs on !
! the DAILY tier, everything else that module held (leaf gas exchange, hydraulics, non-leaf   !
! maintenance respiration) runs sub-daily. Keeping one types module would have made the slow   !
! plant kernels depend on the fast ones and back again -- a genuine cycle, not a style point.  !
!                                                                                          !
! Pure DATA: no methods, no hidden state.                                                    !
!==========================================================================================!
module meds_pheno_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   public :: CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO, CUE_LIGHT

   !=======================================================================================!
   !     PHENOLOGY -- pure SIGNAL kernel: env cues + traits -> two RELATIVE rate tendencies. !
   !     Emits leaf_flush_rate + leaf_shed_rate [1/day]; touches NO carbon, NO leaf/storage   !
   !     state, NO elongf. All leaf/storage carbon update lives in meds_plant_carbon_dynamics.!
   !     The kernel carries TWO governor accumulators (flush_drive, shed_drive) as its memory.!
   !     See docs/dev_plans/MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN.md.                            !
   !=======================================================================================!
   !----- Cue-enable bits. flush_cue_mask and shed_cue_mask select which cues drive each side  !
   !      (min over flush cues, max over shed cues) -- the mechanism for the four target        !
   !      patterns (evergreen flushes on TEMP but never sheds; a leaf-exchanger flushes          !
   !      permissively but sheds on LIGHT).                                                       !
   integer(ik), parameter :: CUE_NONE  = 0_ik    !< no cues (permissive flush / no active shed)
   integer(ik), parameter :: CUE_TEMP  = 1_ik    !< temperature (GDD flush + autumn cold-drop shed)
   integer(ik), parameter :: CUE_WATER = 2_ik    !< soil-water running mean
   integer(ik), parameter :: CUE_HYDRO = 4_ik    !< daily-max leaf water potential (dmax_leaf_psi)
   integer(ik), parameter :: CUE_PHOTO = 8_ik    !< photoperiod (gates the temperature flush)
   integer(ik), parameter :: CUE_LIGHT = 16_ik   !< light: active shed rises with running-mean radiation

   !----- Raw daily environmental drivers (read-only; NO leaf status). ---------------------!
   type :: pheno_env_t
      real(wp)    :: temp_day      = 0.0_wp    !< [K]   daily-mean air/canopy temperature (thermal sums)
      real(wp)    :: soil_temp     = 0.0_wp    !< [K]   shallow-layer soil temperature (cold-drop trigger)
      real(wp)    :: avail_water   = 0.0_wp    !< [-] fraction OR [MPa] soil-water potential (CUE_WATER)
      real(wp)    :: dmax_leaf_psi = 0.0_wp    !< [MPa, <=0] DAILY-MAX leaf water potential (CUE_HYDRO)
      real(wp)    :: rad           = 0.0_wp    !< [W/m2] daily-mean radiation (CUE_LIGHT)
      real(wp)    :: daylength     = 12.0_wp   !< [h]   photoperiod (CUE_PHOTO; caller-supplied)
      integer(ik) :: doy           = 1_ik      !< [-]   day-of-year (thermal-sum season gating)
      logical     :: hemis_north   = .true.    !< northern hemisphere (season gating)
   end type pheno_env_t

   !----- The TWO governor accumulators (the prognostic memory) + cue sub-accumulators. ----!
   type :: pheno_state_t
      real(wp) :: flush_drive   = 1.0_wp   !< [-] smoothed flush permission in [0,1]; 0 => dormant
      real(wp) :: shed_drive    = 0.0_wp   !< [-] smoothed active-shed pressure in [0,1]; 0 => none
      real(wp) :: gdd           = 0.0_wp   !< [K day] growing-degree-day sum            (CUE_TEMP)
      real(wp) :: chill         = 0.0_wp   !< [day]   chilling-day count                (CUE_TEMP)
      real(wp) :: water_avg     = 0.0_wp   !< [-]|[MPa] running-mean available water    (CUE_WATER)
      real(wp) :: low_psi_days  = 0.0_wp   !< [day]   consecutive dry days (dmax<tlp)   (CUE_HYDRO)
      real(wp) :: high_psi_days = 0.0_wp   !< [day]   consecutive wet days (dmax>=.5tlp)(CUE_HYDRO)
      real(wp) :: light_avg     = 0.0_wp   !< [W/m2]  running-mean radiation            (CUE_LIGHT)
   end type pheno_state_t

   !----- Flat per-PFT trait set (self-contained; filled by the driver from cfg%pft). ------!
   type :: pheno_params_t
      !----- Selectors: the two cue masks + shared logistic sharpness. --------------------!
      integer(ik) :: flush_cue_mask = CUE_NONE   !< OR of CUE_* driving the flush side (min)
      integer(ik) :: shed_cue_mask  = CUE_NONE   !< OR of CUE_* driving the shed side (max)
      real(wp)    :: cue_sharpness   = 2.0_wp     !< [-] dimensionless logistic slope (large => ED2-sharp)
      !----- Per-cue transition WIDTHS (normalize each driver so cue_sharpness is a shared slope). !
      real(wp)    :: gdd_width       = 50.0_wp    !< [K day] GDD flush transition width
      real(wp)    :: daylen_width    = 1.0_wp     !< [h]     autumn daylength transition width
      real(wp)    :: soiltemp_width  = 2.0_wp     !< [K]     autumn soil-temperature transition width
      !----- Relative rate scales [1/day] + governor smoothing timescales [day]. -----------!
      real(wp)    :: k_flush_max = 0.06667_wp    !< [1/day] max relative flush rate (~full in 15 d)
      real(wp)    :: k_shed_max  = 0.05_wp       !< [1/day] max relative active-shed rate (~bare in 20 d)
      real(wp)    :: tau_flush   = 5.0_wp        !< [day]   flush governor low-pass timescale
      real(wp)    :: tau_shed    = 5.0_wp        !< [day]   shed governor low-pass timescale
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
      real(wp)    :: water_off_threshold = 0.2_wp      !< available water at which shed = 1 / flush = 0
      real(wp)    :: water_on_threshold  = 0.5_wp      !< available water at which flush = 1 / shed = 0 (> off)
      real(wp)    :: water_window        = 10.0_wp     !< [day] running-mean window
      real(wp)    :: water_width         = 0.1_wp      !< transition width for the water logistics
      !----- Hydraulic (CUE_HYDRO): dmax_leaf_psi consecutive-day counters vs the TLP. ------!
      real(wp)    :: leaf_psi_tlp        = -2.0_wp     !< [MPa] turgor-loss point (Xu 2016)
      real(wp)    :: low_psi_threshold   = 10.0_wp     !< [day] dry days to full shed
      real(wp)    :: high_psi_threshold  = 10.0_wp     !< [day] wet days to full flush
      !----- Photoperiod (CUE_PHOTO): daylength logistic gate (multiplies the flush side). --!
      real(wp)    :: photo_crit          = 11.0_wp     !< [h] critical daylength
      real(wp)    :: photo_slope         = 2.0_wp      !< [1/h] daylength logistic slope
      !----- Light (CUE_LIGHT): active shed rises with running-mean radiation. --------------!
      real(wp)    :: light_on_threshold  = 200.0_wp    !< [W/m2] radiation at which the light shed onsets
      real(wp)    :: light_width         = 50.0_wp     !< [W/m2] transition width for the light shed logistic
      real(wp)    :: light_window        = 10.0_wp     !< [day] running-mean window (ED2 rad_avg)
   end type pheno_params_t

   !----- Outputs: the two RELATIVE rate tendencies + the governing shed cue (diagnostic). --!
   type :: pheno_out_t
      real(wp)    :: leaf_flush_rate = 0.0_wp    !< [1/day] relative flush tendency; 0 == dormancy
      real(wp)    :: leaf_shed_rate  = 0.0_wp    !< [1/day] relative active-shed tendency; 0 == no active shed
      integer(ik) :: cue_limiting    = CUE_NONE  !< strongest active shed cue (argmax over shed_cue_mask)
   end type pheno_out_t

end module meds_pheno_types
