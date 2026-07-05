!==========================================================================================!
! meds_phenology -- the leaf-phenology compute kernel.                                   !
!                                                                                          !
! One daily update, two steps:                                                             !
!   (1) ACCUMULATE the cue memory (state) from today's drivers -- season-gated growing-       !
!       degree-day + chilling sums, the running-mean available water, and the leaf-psi         !
!       consecutive dry/wet-day counters.                                                       !
!   (2) FAVORABILITY -> STATUS -- each active cue maps its (accumulated) driver to a             !
!       favorability f_i in [0,1]; the MOST-LIMITING active cue governs (phi = min_i f_i;        !
!       cue_limiting = argmin_i), and phi is banded into the tri-state via the on/off              !
!       thresholds (the DORMANT deadband is the hysteresis, so no latch is needed).                 !
!                                                                                          !
! Everything is pure/elemental, arithmetic + intrinsics only, with fixed-size flat data --      !
! GPU/SIMD-friendly and reentrant. Combining several simultaneously-active cues more richly     !
! than the limiting-factor min is deferred (a future development).                                !
!==========================================================================================!
module meds_phenology
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi, safe_exp
   use meds_plant_types
   implicit none
   private

   public :: phenology_kernel
   public :: logistic, daylength           ! pure helpers, exposed for tests

   !----- Per-cue transition WIDTHS: normalize each driver so cue_sharpness is a shared,    !
   !      dimensionless slope (larger => sharper, approaching an ED2 hard threshold). -------!
   real(wp), parameter :: gdd_width      = 50.0_wp   !< [K day] GDD flush transition width
   real(wp), parameter :: daylen_width   = 1.0_wp    !< [h]     autumn daylength transition width
   real(wp), parameter :: soiltemp_width = 2.0_wp    !< [K]     autumn soil-temperature transition width

contains

   !=======================================================================================!
   !  Advance one cohort's phenology over one daily step: accumulate, then band the         !
   !  most-limiting favorability into the ON/OFF/DORMANT status.                             !
   !=======================================================================================!
   pure subroutine phenology_kernel(env, params, dt, state, out)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      type(pheno_out_t),    intent(out)   :: out
      real(wp)    :: phi, f, f_photo
      integer(ik) :: lim
      logical     :: have_temp, have_photo

      !----- (1) update the cue memory. --------------------------------------------------!
      call accumulate(env, params, dt, state)

      !----- (2) per-cue favorability -> most-limiting phi. ------------------------------!
      phi = 1.0_wp                              ! CUE_NONE (evergreen) stays fully favorable
      lim = CUE_NONE
      have_temp  = iand(params%cue_mask, CUE_TEMP)  /= 0_ik
      have_photo = iand(params%cue_mask, CUE_PHOTO) /= 0_ik

      f_photo = 1.0_wp
      if (have_photo) f_photo = favor_photo(env, params)

      if (have_temp) then
         f = favor_temp(env, params, state)
         if (have_photo) f = f * f_photo        ! photoperiod gates the temperature cue
         call consider(f, CUE_TEMP, phi, lim)
      else if (have_photo) then
         call consider(f_photo, CUE_PHOTO, phi, lim)   ! photoperiod acting alone
      end if

      if (iand(params%cue_mask, CUE_WATER) /= 0_ik) then
         f = favor_water(state, params)
         call consider(f, CUE_WATER, phi, lim)
      end if

      if (iand(params%cue_mask, CUE_HYDRO) /= 0_ik) then
         f = favor_hydro(state, params)
         call consider(f, CUE_HYDRO, phi, lim)
      end if

      !----- (3) band phi into the directional tri-state (deadband = hysteresis). ---------!
      if (phi > params%phen_on_threshold) then
         out%phenology_status = PHEN_ON
      else if (phi < params%phen_off_threshold) then
         out%phenology_status = PHEN_OFF
      else
         out%phenology_status = PHEN_DORMANT
      end if
      out%cue_limiting = lim
   end subroutine phenology_kernel

   !----- Track the running minimum favorability and the cue that produced it. ------------!
   pure subroutine consider(f, cue, phi, lim)
      real(wp),    intent(in)    :: f
      integer(ik), intent(in)    :: cue
      real(wp),    intent(inout) :: phi
      integer(ik), intent(inout) :: lim
      if (f < phi) then
         phi = f
         lim = cue
      end if
   end subroutine consider

   !=======================================================================================!
   !  Accumulate the cue memory from today's drivers (only for the enabled cues).           !
   !=======================================================================================!
   pure subroutine accumulate(env, params, dt, state)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      integer(ik) :: de
      logical     :: growing, chilling, warm
      real(wp)    :: w

      !----- Thermal sums (season-gated), CUE_TEMP. GDD accumulates on warm days in the      !
      !      growing half of the year and resets outside it; chilling counts cold days in the  !
      !      cool half (so a fresh GDD sum drives each spring flush, ED2 update_thermal_sums).  !
      if (iand(params%cue_mask, CUE_TEMP) /= 0_ik) then
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
      if (iand(params%cue_mask, CUE_WATER) /= 0_ik) then
         w = min(1.0_wp, dt / max(params%water_window, dt))
         state%water_avg = state%water_avg + w * (env%avail_water - state%water_avg)
      end if

      !----- Leaf-psi consecutive-day counters, CUE_HYDRO. --------------------------------!
      if (iand(params%cue_mask, CUE_HYDRO) /= 0_ik) then
         if (env%psi_leaf < params%leaf_psi_tlp) then
            state%low_psi_days = state%low_psi_days + dt
         else
            state%low_psi_days = 0.0_wp
         end if
         if (env%psi_leaf >= 0.5_wp * params%leaf_psi_tlp) then
            state%high_psi_days = state%high_psi_days + dt
         else
            state%high_psi_days = 0.0_wp
         end if
      end if
   end subroutine accumulate

   !=======================================================================================!
   !  Per-cue favorabilities (all in [0,1]; high = favors leaves).                          !
   !=======================================================================================!

   !----- Temperature: chilling-adaptive GDD flush, cut by the autumn cold-drop predicate. -!
   pure real(wp) function favor_temp(env, params, state) result(f)
      type(pheno_env_t),    intent(in) :: env
      type(pheno_params_t), intent(in) :: params
      type(pheno_state_t),  intent(in) :: state
      real(wp) :: k, gdd_thresh, f_gdd, g_day, g_st1, g_st2, drop
      k          = params%cue_sharpness
      gdd_thresh = params%phen_a + params%phen_b * safe_exp(params%phen_c * state%chill)
      f_gdd      = logistic(k * (state%gdd - gdd_thresh) / gdd_width)
      !----- Autumn drop: (short day AND cool soil) OR very cold soil (each smooth in [0,1]).-!
      g_day = logistic(k * (params%cold_drop_daylength - env%daylength) / daylen_width)
      g_st1 = logistic(k * (params%cold_drop_soiltemp1 - env%soil_temp) / soiltemp_width)
      g_st2 = logistic(k * (params%cold_drop_soiltemp2 - env%soil_temp) / soiltemp_width)
      drop  = max(g_day * g_st1, g_st2)
      f     = f_gdd * (1.0_wp - drop)
   end function favor_temp

   !----- Water: linear ramp of the running-mean available water between off and on. -------!
   pure real(wp) function favor_water(state, params) result(f)
      type(pheno_state_t),  intent(in) :: state
      type(pheno_params_t), intent(in) :: params
      real(wp) :: denom
      denom = params%water_on_threshold - params%water_off_threshold
      if (abs(denom) < tiny(1.0_wp)) then
         if (state%water_avg >= params%water_on_threshold) then
            f = 1.0_wp
         else
            f = 0.0_wp
         end if
      else
         f = clamp01((state%water_avg - params%water_off_threshold) / denom)
      end if
   end function favor_water

   !----- Hydraulic: rises with sustained wet days, falls with sustained dry days;          !
   !      neither sustained => 0.5 (the DORMANT deadband).                                   !
   pure real(wp) function favor_hydro(state, params) result(f)
      type(pheno_state_t),  intent(in) :: state
      type(pheno_params_t), intent(in) :: params
      real(wp) :: hi, lo
      hi = clamp01(state%high_psi_days / max(params%high_psi_threshold, tiny(1.0_wp)))
      lo = clamp01(state%low_psi_days  / max(params%low_psi_threshold,  tiny(1.0_wp)))
      f  = clamp01(0.5_wp + 0.5_wp * hi - 0.5_wp * lo)
   end function favor_hydro

   !----- Photoperiod: logistic gate on daylength (multiplies the temperature cue). --------!
   pure real(wp) function favor_photo(env, params) result(f)
      type(pheno_env_t),    intent(in) :: env
      type(pheno_params_t), intent(in) :: params
      f = logistic(params%photo_slope * (env%daylength - params%photo_crit))
   end function favor_photo

   !=======================================================================================!
   !  Pure helpers.                                                                         !
   !=======================================================================================!

   !----- Numerically safe logistic (smooth 0->1 step). -----------------------------------!
   elemental real(wp) function logistic(z) result(f)
      real(wp), intent(in) :: z
      f = 1.0_wp / (1.0_wp + safe_exp(-z))
   end function logistic

   !----- Clamp to [0,1]. -----------------------------------------------------------------!
   elemental real(wp) function clamp01(x) result(y)
      real(wp), intent(in) :: x
      y = min(1.0_wp, max(0.0_wp, x))
   end function clamp01

   !----- Effective (northern-hemisphere-equivalent) day-of-year: shift SH by half a year. -!
   elemental integer(ik) function doy_effective(doy, hemis_north) result(de)
      integer(ik), intent(in) :: doy
      logical,     intent(in) :: hemis_north
      if (hemis_north) then
         de = doy
      else
         de = modulo(doy - 1_ik + 182_ik, 365_ik) + 1_ik
      end if
   end function doy_effective

   !----- Daylength [h] from latitude [deg] + day-of-year (White et al. 1997 form). The ED2  !
   !      polar branch is FIXED here: polar DAY is arg <= -1 (ED2 wrote arg <= 1, a bug). ---!
   elemental real(wp) function daylength(lat_deg, doy) result(dl)
      real(wp),    intent(in) :: lat_deg
      integer(ik), intent(in) :: doy
      real(wp) :: latr, decl, arg
      latr = lat_deg * pi / 180.0_wp
      decl = -23.44_wp * pi / 180.0_wp * cos(2.0_wp * pi / 365.0_wp * (real(doy, wp) + 9.0_wp))
      arg  = -tan(latr) * tan(decl)
      if (arg >= 1.0_wp) then
         dl = 0.0_wp                 ! polar night
      else if (arg <= -1.0_wp) then
         dl = 24.0_wp                ! polar day (FIX: ED2's branch read arg <= 1)
      else
         dl = 24.0_wp / pi * acos(arg)
      end if
   end function daylength

end module meds_phenology
