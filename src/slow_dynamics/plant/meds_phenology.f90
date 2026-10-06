! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_phenology -- the leaf-phenology kernel.                                               !
!                                                                                          !
! A signal generator: daily cues + per-PFT traits -> two smoothed tendencies in [0,1]       !
! (leaf_flush_tendency, leaf_shed_tendency). It touches no carbon; leaf_turnover_step, here  !
! too, is the one place the tendencies become leaf growth and leaf loss, for the carbon      !
! layer and for the Python mirror alike.                                                    !
!                                                                                          !
! One daily update (docs/science/plant_phenology.md):                                        !
!   (1) ACCUMULATE the cue memory: warmth above flush_base_temp from midwinter, cold below    !
!       shed_base_temp from midsummer, predawn leaf psi above and below the turgor-loss point, !
!       and the running-mean hours of light (PAR at the cohort's top above par_min).           !
!   (2) SWITCH each cue through sigma(s (x - x*)), one centre and one signed sharpness each.   !
!   (3) COMBINE: flush = the product of its cues' switches (every cue must permit flushing);   !
!       shed = the larger of the seasonal trigger (temperature x light) and the water        !
!       trigger.                                                                               !
!   (4) SMOOTH each signal into its tendency over a timescale; the potential relative rates    !
!       are rate_max * tendency.                                                               !
! Pure, scalar, arithmetic only: device- and SIMD-friendly, and reentrant.                     !
!==========================================================================================!
module meds_phenology
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : yr_day
   use meds_numerics,    only : logistic, clamp01
   use meds_time,        only : doy_effective
   use meds_phenology_types
   implicit none
   private

   public :: phenology_kernel, leaf_turnover_step

   !----- Northern-equivalent day of the summer solstice: the cold sum counts from here, while !
   !      days shorten (doy_effective shifts the southern hemisphere by half a year).          !
   integer(ik), parameter :: MIDSUMMER_DOY     = 172_ik
   !----- A flush rate at or below this is dormant: a senescing canopy then snaps to bare     !
   !      rather than leaving an exponential tail (a numerical "off", not a trait).            !
   real(wp),    parameter :: DORMANT_FLUSH_EPS = 1.0e-6_wp    !< [1/day]

contains

   !=======================================================================================!
   !  Advance one cohort's phenology over one step of dt days.                               !
   !=======================================================================================!
   pure subroutine phenology_kernel(env, params, dt, state, out)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      type(pheno_state_t),  intent(inout) :: state
      type(pheno_out_t),    intent(out)   :: out
      integer(ik) :: both
      real(wp)    :: wet_before, dry_before, wet_switch, dry_switch, s_flush, seasonal, s_shed, w

      both = ior(params%flush_cue_mask, params%shed_cue_mask)
      wet_before = wet_water_switch(params, state)
      dry_before = dry_water_switch(params, state)
      call accumulate(env, params, dt, both, state)

      !----- Water has no calendar: each sum resets when the OTHER side's switch crosses 0.5  !
      !      upward today (the event, not the level). A brief rain adds wet credit without     !
      !      wiping a drought, and a long wet season cannot block the next drought.            !
      if (iand(both, CUE_WATER) /= 0_ik) then
         if (wet_before <= 0.5_wp .and. wet_water_switch(params, state) > 0.5_wp) state%dry_psi_sum = 0.0_wp
         if (dry_before <= 0.5_wp .and. dry_water_switch(params, state) > 0.5_wp) state%wet_psi_sum = 0.0_wp
      end if
      wet_switch = wet_water_switch(params, state)
      dry_switch = dry_water_switch(params, state)

      !----- (3a) Flush: every enabled cue must permit it (an empty mask always flushes). ----!
      s_flush = 1.0_wp
      if (iand(params%flush_cue_mask, CUE_TEMP) /= 0_ik) s_flush = s_flush *                     &
         logistic(params%flush_temp_sharpness * (state%growing_degree_days - params%flush_degree_days))
      if (iand(params%flush_cue_mask, CUE_LIGHT) /= 0_ik) s_flush = s_flush *                    &
         logistic(params%flush_light_sharpness * (state%light_hours_mean - params%flush_light_hours))
      if (iand(params%flush_cue_mask, CUE_WATER) /= 0_ik) s_flush = s_flush * wet_switch

      !----- (3b) Shed: the seasonal trigger needs all of its enabled cues (cold AND the     !
      !      light condition); the water trigger acts on its own. An empty mask never         !
      !      senesces.                                                                        !
      seasonal = 0.0_wp
      if (iand(params%shed_cue_mask, CUE_TEMP + CUE_LIGHT) /= 0_ik) then
         seasonal = 1.0_wp
         if (iand(params%shed_cue_mask, CUE_TEMP) /= 0_ik) seasonal = seasonal *                 &
            logistic(params%shed_temp_sharpness * (state%cold_degree_days - params%shed_degree_days))
         if (iand(params%shed_cue_mask, CUE_LIGHT) /= 0_ik) seasonal = seasonal *                &
            logistic(params%shed_light_sharpness * (state%light_hours_mean - params%shed_light_hours))
      end if
      s_shed = seasonal
      if (iand(params%shed_cue_mask, CUE_WATER) /= 0_ik) s_shed = max(s_shed, dry_switch)

      !----- (4) Smooth into the tendencies (the weight caps at 1 when dt exceeds the       !
      !      timescale) and map to the potential rates.                                     !
      w = dt / max(params%flush_cue_timescale, dt)
      state%leaf_flush_tendency = clamp01(state%leaf_flush_tendency + w * (s_flush - state%leaf_flush_tendency))
      w = dt / max(params%shed_cue_timescale, dt)
      state%leaf_shed_tendency  = clamp01(state%leaf_shed_tendency  + w * (s_shed  - state%leaf_shed_tendency))
      out%leaf_flush_potential = params%flush_rate_max * state%leaf_flush_tendency
      out%leaf_shed_potential  = params%shed_rate_max  * state%leaf_shed_tendency
   end subroutine phenology_kernel

   !----- The two water switches, from the sums of predawn leaf psi above and below the TLP. -!
   pure real(wp) function wet_water_switch(params, state) result(sw)
      type(pheno_params_t), intent(in) :: params
      type(pheno_state_t),  intent(in) :: state
      sw = logistic(params%flush_water_sharpness * (state%wet_psi_sum - params%flush_water_sum))
   end function wet_water_switch

   pure real(wp) function dry_water_switch(params, state) result(sw)
      type(pheno_params_t), intent(in) :: params
      type(pheno_state_t),  intent(in) :: state
      sw = logistic(params%shed_water_sharpness * (state%dry_psi_sum - params%shed_water_sum))
   end function dry_water_switch

   !=======================================================================================!
   !  Accumulate the cue memory from today's drivers, for the cues either side uses.         !
   !=======================================================================================!
   pure subroutine accumulate(env, params, dt, both, state)
      type(pheno_env_t),    intent(in)    :: env
      type(pheno_params_t), intent(in)    :: params
      real(wp),             intent(in)    :: dt
      integer(ik),          intent(in)    :: both
      type(pheno_state_t),  intent(inout) :: state
      integer(ik) :: de
      real(wp)    :: w

      !----- Temperature: both sums restart at midwinter (the first step of the northern-   !
      !      equivalent year); warmth counts all year, cold only once days shorten.          !
      if (iand(both, CUE_TEMP) /= 0_ik) then
         de = doy_effective(env%doy, env%hemis_north)
         if (real(de, wp) < 1.0_wp + dt) then
            state%growing_degree_days = 0.0_wp
            state%cold_degree_days    = 0.0_wp
         end if
         state%growing_degree_days = state%growing_degree_days + max(0.0_wp, env%temp_day - params%flush_base_temp) * dt
         if (de >= MIDSUMMER_DOY)                                                                &
            state%cold_degree_days = state%cold_degree_days + max(0.0_wp, params%shed_base_temp - env%temp_day) * dt
      end if

      !----- Water: predawn leaf psi above and below the turgor-loss point. ----------------!
      if (iand(both, CUE_WATER) /= 0_ik) then
         state%wet_psi_sum = state%wet_psi_sum + max(0.0_wp, env%predawn_leaf_psi - params%leaf_psi_tlp) * dt
         state%dry_psi_sum = state%dry_psi_sum + max(0.0_wp, params%leaf_psi_tlp - env%predawn_leaf_psi) * dt
      end if

      !----- Light: an exponential running mean of the daily hours of light. A cohort with  !
      !      no light memory yet (a cold start or a recruit) starts from today's hours, so its  !
      !      first days do not read as darkness.                                             !
      if (iand(both, CUE_LIGHT) /= 0_ik) then
         if (state%light_hours_mean < 0.0_wp) then
            state%light_hours_mean = env%par_hours
         else
            w = min(1.0_wp, dt / max(params%light_window, dt))
            state%light_hours_mean = state%light_hours_mean + w * (env%par_hours - state%light_hours_mean)
         end if
      end if
   end subroutine accumulate

   !=======================================================================================!
   !  One step of leaf loss and the flush cap, in whatever carbon unit `leaf` is in. The      !
   !  carbon layer calls it per cohort; the C-API exposes it so Python applies the same rule. !
   !                                                                                         !
   !    background = leaf_turnover_rate * leaf_flush_tendency * leaf * dt   (old leaves turned !
   !                 over while new ones grow; none in dormancy)                               !
   !    senescence = shed_rate_max * leaf_shed_tendency * leaf * dt, stopping where leaf_cover !
   !                 (leaf / leaf_full) reaches min_leaf_cover                                 !
   !    flush_cap  = flush_rate_max * leaf_flush_tendency * leaf_full * dt                     !
   !                                                                                         !
   !  A dormant canopy (flush ~ 0) that senescence would leave below bare_leaf_cover snaps to  !
   !  bare, when its floor lies below that threshold (a deciduous PFT).                        !
   !=======================================================================================!
   pure subroutine leaf_turnover_step(leaf, leaf_full, leaf_flush_tendency, leaf_shed_tendency,   &
                                      flush_rate_max, shed_rate_max, leaf_turnover_rate,          &
                                      min_leaf_cover, bare_leaf_cover, dt,                        &
                                      senescence, background, flush_cap)
      real(wp), intent(in)  :: leaf, leaf_full                          !< [carbon] now and full canopy
      real(wp), intent(in)  :: leaf_flush_tendency, leaf_shed_tendency  !< [-]
      real(wp), intent(in)  :: flush_rate_max, shed_rate_max            !< [1/day]
      real(wp), intent(in)  :: leaf_turnover_rate                       !< [1/yr] background loss
      real(wp), intent(in)  :: min_leaf_cover, bare_leaf_cover          !< [-] fractions of leaf_full
      real(wp), intent(in)  :: dt                                       !< [day]
      real(wp), intent(out) :: senescence, background, flush_cap        !< [carbon] this step
      real(wp) :: pool

      pool       = max(leaf, 0.0_wp)
      background = min(max(leaf_turnover_rate, 0.0_wp) / yr_day * leaf_flush_tendency * pool * dt, pool)
      senescence = min(max(shed_rate_max, 0.0_wp) * leaf_shed_tendency * pool * dt,                  &
                       max(0.0_wp, pool - background - min_leaf_cover * leaf_full))
      if (flush_rate_max * leaf_flush_tendency <= DORMANT_FLUSH_EPS .and. senescence > 0.0_wp      &
          .and. min_leaf_cover < bare_leaf_cover                                                  &
          .and. pool - background - senescence < bare_leaf_cover * leaf_full)                     &
         senescence = pool - background
      flush_cap  = max(flush_rate_max, 0.0_wp) * leaf_flush_tendency * leaf_full * dt
   end subroutine leaf_turnover_step

end module meds_phenology
