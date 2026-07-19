!==========================================================================================!
! meds_phenology -- the leaf-phenology SIGNAL kernel.                                       !
!                                                                                          !
! A pure signal generator: environment cues + per-PFT traits -> TWO relative rate tendencies  !
! (leaf_flush_rate, leaf_shed_rate, both [1/day]). It touches NO carbon, NO leaf/storage state, !
! and NO elongf -- all leaf/storage carbon update lives downstream in meds_plant_carbon_dynamics.!
!                                                                                          !
! One daily update, four steps:                                                             !
!   (1) ACCUMULATE the cue memory -- season-gated GDD/chilling sums, the running-mean available  !
!       water, the dmax-leaf-psi consecutive dry/wet-day counters, and the running-mean radiation. !
!   (2) PER-CUE (s_flush, s_shed) in [0,1] -- each cue supplies a flush signal and a shed signal.  !
!   (3) COMBINE over the TWO masks: s_flush = MIN over flush_cue_mask (build only when every flush   !
!       cue is clear); s_shed = MAX over shed_cue_mask (any shed cue commands senescence). Splitting  !
!       the mask is what expresses the four target patterns (evergreen flushes on TEMP but never       !
!       sheds; a leaf-exchanger flushes permissively but sheds on LIGHT).                               !
!   (4) LOW-PASS the two governors (flush_drive, shed_drive -- the prognostic memory) and map to rates: !
!       leaf_flush_rate = k_flush_max*flush_drive; leaf_shed_rate = k_shed_max*shed_drive.               !
!                                                                                          !
! Everything is pure/scalar, arithmetic + intrinsics only, fixed-size flat data -- GPU/SIMD-friendly    !
! and reentrant. logistic/clamp01 come from meds_numerics; doy_effective from meds_time.                !
!==========================================================================================!
module meds_phenology
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : safe_exp
   use meds_numerics,    only : logistic, clamp01
   use meds_time,        only : doy_effective
   use meds_plant_types
   implicit none
   private

   public :: phenology_kernel, pheno_drives_to_rates

   !----- Per-cue transition WIDTHS: normalize each driver so cue_sharpness is a shared,    !
   !      dimensionless slope (larger => sharper, approaching an ED2 hard threshold). -------!
   real(wp), parameter :: gdd_width      = 50.0_wp   !< [K day] GDD flush transition width
   real(wp), parameter :: daylen_width   = 1.0_wp    !< [h]     autumn daylength transition width
   real(wp), parameter :: soiltemp_width = 2.0_wp    !< [K]     autumn soil-temperature transition width

contains

   !=======================================================================================!
   !  Advance one cohort's phenology over one daily step: accumulate the cue memory, form the !
   !  per-cue signals, combine over the two masks, low-pass the governors, emit two rates.    !
   !=======================================================================================!
   pure subroutine phenology_kernel(env, params, dt, state, out)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      type(pheno_out_t),    intent(out)   :: out
      integer(ik) :: both, lim
      real(wp)    :: k, s_flush, s_shed, f, f_photo, w_f, w_s
      real(wp)    :: t_flush, t_shed, gdd_thresh, g_day, g_st1, g_st2

      both = ior(params%flush_cue_mask, params%shed_cue_mask)

      !----- (1) update the cue memory (only for cues enabled on either side). ------------!
      call accumulate(env, params, dt, both, state)

      !----- (2) temperature signals (computed once; used by whichever side lists TEMP). --!
      k = params%cue_sharpness
      t_flush = 0.0_wp ; t_shed = 0.0_wp
      if (iand(both, CUE_TEMP) /= 0_ik) then
         gdd_thresh = params%phen_a + params%phen_b * safe_exp(params%phen_c * state%chill)
         t_flush    = logistic(k * (state%gdd - gdd_thresh) / gdd_width)
         g_day      = logistic(k * (params%cold_drop_daylength - env%daylength) / daylen_width)
         g_st1      = logistic(k * (params%cold_drop_soiltemp1 - env%soil_temp) / soiltemp_width)
         g_st2      = logistic(k * (params%cold_drop_soiltemp2 - env%soil_temp) / soiltemp_width)
         t_shed     = max(g_day * g_st1, g_st2)
      end if

      !----- (3a) FLUSH signal = MIN over the flush cues (CUE_NONE => permissive 1). -------!
      s_flush = 1.0_wp
      f_photo = 1.0_wp
      if (iand(params%flush_cue_mask, CUE_PHOTO) /= 0_ik)                                       &
         f_photo = logistic(params%photo_slope * (env%daylength - params%photo_crit))
      if (iand(params%flush_cue_mask, CUE_TEMP) /= 0_ik) then
         f = t_flush
         if (iand(params%flush_cue_mask, CUE_PHOTO) /= 0_ik) f = f * f_photo   ! photoperiod gates temp flush
         s_flush = min(s_flush, f)
      else if (iand(params%flush_cue_mask, CUE_PHOTO) /= 0_ik) then
         s_flush = min(s_flush, f_photo)                                        ! photoperiod acting alone
      end if
      if (iand(params%flush_cue_mask, CUE_WATER) /= 0_ik)                                       &
         s_flush = min(s_flush, logistic(k * (state%water_avg - params%water_on_threshold)      &
                                           / max(params%water_width, tiny(1.0_wp))))
      if (iand(params%flush_cue_mask, CUE_HYDRO) /= 0_ik)                                       &
         s_flush = min(s_flush, clamp01(state%high_psi_days / max(params%high_psi_threshold, tiny(1.0_wp))))
      !----- CUE_LIGHT contributes a flush signal of 1 (non-limiting) -- no term needed. ---!

      !----- (3b) SHED signal = MAX over the shed cues (CUE_NONE => no active shed 0). -----!
      s_shed = 0.0_wp ; lim = CUE_NONE
      if (iand(params%shed_cue_mask, CUE_TEMP) /= 0_ik)  call consider(t_shed, CUE_TEMP,  s_shed, lim)
      if (iand(params%shed_cue_mask, CUE_WATER) /= 0_ik)                                        &
         call consider(logistic(k * (params%water_off_threshold - state%water_avg)              &
                                  / max(params%water_width, tiny(1.0_wp))), CUE_WATER, s_shed, lim)
      if (iand(params%shed_cue_mask, CUE_HYDRO) /= 0_ik)                                        &
         call consider(clamp01(state%low_psi_days / max(params%low_psi_threshold, tiny(1.0_wp))), &
                       CUE_HYDRO, s_shed, lim)
      if (iand(params%shed_cue_mask, CUE_LIGHT) /= 0_ik)                                        &
         call consider(logistic(k * (state%light_avg - params%light_on_threshold)              &
                                  / max(params%light_width, tiny(1.0_wp))), CUE_LIGHT, s_shed, lim)

      !----- (4) low-pass the two governors (guarded weight caps at 1; FPE-safe), map to rates. -!
      w_f = dt / max(params%tau_flush, dt)
      w_s = dt / max(params%tau_shed,  dt)
      state%flush_drive = clamp01(state%flush_drive + w_f * (s_flush - state%flush_drive))
      state%shed_drive  = clamp01(state%shed_drive  + w_s * (s_shed  - state%shed_drive))

      out%leaf_flush_rate = params%k_flush_max * state%flush_drive
      out%leaf_shed_rate  = params%k_shed_max  * state%shed_drive
      out%cue_limiting    = lim
   end subroutine phenology_kernel

   !----- Track the running MAXIMUM shed signal and the cue that produced it. --------------!
   pure subroutine consider(f, cue, s, lim)
      real(wp),    intent(in)    :: f
      integer(ik), intent(in)    :: cue
      real(wp),    intent(inout) :: s
      integer(ik), intent(inout) :: lim
      if (f > s) then
         s   = f
         lim = cue
      end if
   end subroutine consider

   !=======================================================================================!
   !  Accumulate the cue memory from today's drivers (only for the enabled cues -- `both` is  !
   !  the OR of the two masks). Identical thermal/water/hydro logic to the tri-state kernel,  !
   !  plus the light running mean; the hydraulic counters key on the DAILY-MAX leaf psi.      !
   !=======================================================================================!
   pure subroutine accumulate(env, params, dt, both, state)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      integer(ik),          intent(in)    :: both
      type(pheno_state_t),  intent(inout) :: state
      integer(ik) :: de
      logical     :: growing, chilling, warm
      real(wp)    :: w

      !----- Thermal sums (season-gated), CUE_TEMP. --------------------------------------!
      if (iand(both, CUE_TEMP) /= 0_ik) then
         de       = doy_effective(env%doy, env%hemis_north)
         growing  = de <= 244_ik
         chilling = de >= 305_ik .or. de <= 181_ik
         warm     = env%temp_day > params%gdd_base_temp
         if (growing) then
            if (warm) state%gdd = state%gdd + (env%temp_day - params%gdd_base_temp) * dt
         else
            state%gdd = 0.0_wp
         end if
         if (chilling) then
            if (env%temp_day < params%chill_base_temp) state%chill = state%chill + dt
         else
            state%chill = 0.0_wp
         end if
      end if

      !----- Soil-water running mean (exponential), CUE_WATER. ----------------------------!
      if (iand(both, CUE_WATER) /= 0_ik) then
         w = min(1.0_wp, dt / max(params%water_window, dt))
         state%water_avg = state%water_avg + w * (env%avail_water - state%water_avg)
      end if

      !----- Daily-max-leaf-psi consecutive-day counters, CUE_HYDRO. ----------------------!
      if (iand(both, CUE_HYDRO) /= 0_ik) then
         if (env%dmax_leaf_psi < params%leaf_psi_tlp) then
            state%low_psi_days = state%low_psi_days + dt
         else
            state%low_psi_days = 0.0_wp
         end if
         if (env%dmax_leaf_psi >= 0.5_wp * params%leaf_psi_tlp) then
            state%high_psi_days = state%high_psi_days + dt
         else
            state%high_psi_days = 0.0_wp
         end if
      end if

      !----- Radiation running mean (exponential), CUE_LIGHT (ED2 rad_avg). ---------------!
      if (iand(both, CUE_LIGHT) /= 0_ik) then
         w = min(1.0_wp, dt / max(params%light_window, dt))
         state%light_avg = state%light_avg + w * (env%rad - state%light_avg)
      end if
   end subroutine accumulate

   !=======================================================================================!
   !  Map the two governor drives to the two RELATIVE rate tendencies [1/day]. A subroutine  !
   !  with scalar intent(out) args (issue-#7 safe: never an array-valued function result fed  !
   !  into a call). The consumer (carbon_growth) calls this from the stored cohort drives.    !
   !=======================================================================================!
   pure subroutine pheno_drives_to_rates(flush_drive, shed_drive, k_flush_max, k_shed_max,     &
                                         leaf_flush_rate, leaf_shed_rate)
      real(wp), intent(in)  :: flush_drive, shed_drive, k_flush_max, k_shed_max
      real(wp), intent(out) :: leaf_flush_rate, leaf_shed_rate
      leaf_flush_rate = k_flush_max * flush_drive
      leaf_shed_rate  = k_shed_max  * shed_drive
   end subroutine pheno_drives_to_rates

end module meds_phenology
