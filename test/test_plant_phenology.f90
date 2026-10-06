! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_phenology -- unit tests for the leaf-phenology kernel and leaf_turnover_step.  !
!                                                                                          !
! The kernel turns daily cues into two smoothed tendencies (leaf_flush_tendency,            !
! leaf_shed_tendency) and the potential rates rate_max * tendency. The tests:                !
!                                                                                          !
!   1. NO CUES             : both masks empty => flush potential = flush_rate_max, shed 0.    !
!   2. TEMPERATE DECIDUOUS : warmth x day-length flush, cold x short-day shed; the sums reset  !
!                            at midwinter and the cold sum waits for midsummer.              !
!   3. SOUTHERN HEMISPHERE : the same PFT half a year later.                                 !
!   4. DROUGHT DECIDUOUS   : wet/dry psi sums; a long wet season does not block a drought, a  !
!                            brief rain does not wipe one.                                    !
!   5. LIGHT EXCHANGING    : bright running-mean shortwave triggers senescence, flush stays on.!
!   6. RATE MAPPING        : potential == rate_max * tendency, tendencies in [0,1].            !
!   7. DEGENERATE / FPE    : extreme drivers keep the potentials finite and >= 0.             !
!   8. DAYLENGTH           : the polar branches of meds_time daylength.                       !
!   9. LEAF TURNOVER STEP  : background, senescence, the leaf-cover floor (an emergent        !
!                            evergreen), the dormant snap to bare, the pool clamp.            !
!==========================================================================================!
program test_plant_phenology
   use meds_test_assert, only : check_close, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : yr_day
   use meds_time,            only : daylength
   use meds_phenology_types, only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,     &
                                    CUE_NONE, CUE_TEMP, CUE_LIGHT, CUE_WATER, CUE_ALL,           &
                                    LIGHT_RADIATION
   use meds_phenology,       only : phenology_kernel, leaf_turnover_step
   implicit none

   real(wp), parameter :: twopi = 6.283185307179586_wp
   real(wp), parameter :: LAT   = 42.5_wp                 ! Harvard Forest

   call test_no_cues()
   call test_temperate_deciduous()
   call test_southern_hemisphere()
   call test_drought_deciduous()
   call test_light_exchanging()
   call test_rate_mapping()
   call test_degenerate()
   call test_daylength_polar()
   call test_leaf_turnover_step()

   call test_report('test_plant_phenology')

contains

   !----- Annual air temperature (northern hemisphere; warmest ~ doy 201). -----------------!
   pure real(wp) function annual_temp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 12.0_wp * cos(twopi * real(doy - 201_ik, wp) / 365.0_wp)
   end function annual_temp

   !----- A temperate deciduous PFT (warmth x day length, cold x short days). ---------------!
   pure function temperate_deciduous() result(p)
      type(pheno_params_t) :: p
      p%flush_cue_mask        = CUE_TEMP + CUE_LIGHT
      p%shed_cue_mask         = CUE_TEMP + CUE_LIGHT
      p%flush_degree_days     = 92.0_wp
      p%shed_base_temp        = 290.37_wp
      p%shed_degree_days      = 48.0_wp
      p%flush_light_threshold = 10.35_wp
      p%shed_light_threshold  = 9.83_wp
   end function temperate_deciduous

   !----- 1. No cues: always flushing, never senescing, under ANY drivers. -----------------!
   subroutine test_no_cues()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: flush_ok, shed_ok
      print '(a)', '-- 1. no cues (flush={}, shed={}) --'
      params%flush_cue_mask = CUE_NONE ; params%shed_cue_mask = CUE_NONE
      flush_ok = .true. ; shed_ok = .true.
      do d = 1_ik, 400_ik
         env%doy              = modulo(d - 1_ik, 365_ik) + 1_ik
         env%temp_day         = annual_temp(env%doy)
         env%predawn_leaf_psi = -8.0_wp                  ! cavitated ...
         env%rad              = 800.0_wp                 ! ... and blazing -- none of it is a cue
         env%daylength        = daylength(LAT, env%doy)
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (abs(out%leaf_flush_potential - params%flush_rate_max) > 1.0e-12_wp) flush_ok = .false.
         if (out%leaf_shed_potential /= 0.0_wp) shed_ok = .false.
      end do
      call check_true('no cues: flush potential = flush_rate_max regardless of drivers', flush_ok)
      call check_true('no cues: shed potential = 0 regardless of drivers',               shed_ok)
   end subroutine test_no_cues

   !----- Run a temperate deciduous PFT for two years and sample the second. ----------------!
   subroutine run_two_years(north, offset, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      logical,     intent(in)  :: north
      integer(ik), intent(in)  :: offset          !< calendar lag of the climate [day] (183: doy_effective)
      real(wp),    intent(out) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      type(pheno_env_t)    :: env
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      type(pheno_params_t) :: params
      integer(ik) :: d, de
      real(wp)    :: lat_signed
      params = temperate_deciduous()
      lat_signed = merge(LAT, -LAT, north)
      env%hemis_north = north
      do d = 1_ik, 730_ik
         env%doy       = modulo(d - 1_ik + offset, 365_ik) + 1_ik
         de            = modulo(d - 1_ik, 365_ik) + 1_ik          ! the season, northern-equivalent
         env%temp_day  = annual_temp(de)
         env%daylength = daylength(lat_signed, env%doy)
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (d == 366_ik) gdd_1   = state%growing_degree_days
         if (d == 536_ik) cdd_171 = state%cold_degree_days
         if (d == 565_ik) then
            fl_200 = state%leaf_flush_tendency ; sh_200 = state%leaf_shed_tendency
            gdd_200 = state%growing_degree_days
         end if
         if (d == 705_ik) then
            fl_340 = state%leaf_flush_tendency ; sh_340 = state%leaf_shed_tendency
         end if
      end do
   end subroutine run_two_years

   !----- 2. Temperate deciduous: flush in summer, senescence in late autumn. --------------!
   subroutine test_temperate_deciduous()
      real(wp) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      print '(a)', '-- 2. temperate deciduous (flush={TEMP,LIGHT}, shed={TEMP,LIGHT}) --'
      call run_two_years(.true., 0_ik, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      call check_true('deciduous: flush tendency HIGH in mid-summer', fl_200 > 0.9_wp, fl_200)
      call check_true('deciduous: shed tendency  LOW  in mid-summer', sh_200 < 0.05_wp, sh_200)
      call check_true('deciduous: shed tendency  HIGH in late autumn', sh_340 > 0.5_wp, sh_340)
      !----- The warmth sum is still high in autumn: it is the day-length gate that closes. --!
      call check_true('deciduous: flush tendency LOW in late autumn (day-length gate)', fl_340 < 0.5_wp, fl_340)
      call check_true('deciduous: warmth sum restarts at midwinter',   gdd_1 < 1.0_wp, gdd_1)
      call check_true('deciduous: warmth sum accumulated by summer',   gdd_200 > 92.0_wp, gdd_200)
      call check_close('deciduous: no cold sum before midsummer',      cdd_171, 0.0_wp, 0.0_wp)
   end subroutine test_temperate_deciduous

   !----- 3. Southern hemisphere: the same seasons half a year later on the calendar. -----!
   subroutine test_southern_hemisphere()
      real(wp) :: fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171
      print '(a)', '-- 3. southern hemisphere --'
      call run_two_years(.false., 183_ik, fl_200, sh_200, fl_340, sh_340, gdd_1, gdd_200, cdd_171)
      call check_true('south: flush tendency HIGH in its mid-summer', fl_200 > 0.9_wp, fl_200)
      call check_true('south: shed tendency  LOW  in its mid-summer', sh_200 < 0.05_wp, sh_200)
      call check_true('south: shed tendency  HIGH in its late autumn', sh_340 > 0.5_wp, sh_340)
      call check_true('south: warmth sum restarts at its midwinter',  gdd_1 < 1.0_wp, gdd_1)
   end subroutine test_southern_hemisphere

   !----- 4. Drought deciduous: the wet and dry sums of predawn psi against the TLP. -------!
   subroutine test_drought_deciduous()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      real(wp)    :: dry_before_rain
      print '(a)', '-- 4. drought deciduous (flush={WATER}, shed={WATER}) --'
      params%flush_cue_mask = CUE_WATER ; params%shed_cue_mask = CUE_WATER   ! tlp -2, sums 10
      !----- (a) A long wet season: the wet sum grows far past its centre. -----------------!
      do d = 1_ik, 60_ik
         env%predawn_leaf_psi = -0.5_wp                     ! 1.5 MPa above the TLP each day
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('drought-decid: watered => no senescence', out%leaf_shed_potential < 0.05_wp * params%shed_rate_max)
      call check_true('drought-decid: watered => flushing',      out%leaf_flush_potential > 0.9_wp * params%flush_rate_max)
      !----- (b) Drought: the dry sum crosses its centre DESPITE the large wet sum (the     !
      !      crossing resets the wet sum; a level test would wipe the dry sum every day). -!
      do d = 1_ik, 30_ik
         env%predawn_leaf_psi = -3.0_wp                     ! 1 MPa below the TLP each day
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('drought-decid: drought after a wet season => senescence', &
                      out%leaf_shed_potential > 0.9_wp * params%shed_rate_max, out%leaf_shed_potential)
      call check_true('drought-decid: drought => flush falls',   out%leaf_flush_potential < 0.1_wp * params%flush_rate_max)
      !----- (c) A brief rain adds wet credit but does not wipe the drought. ---------------!
      dry_before_rain = state%dry_psi_sum
      do d = 1_ik, 2_ik
         env%predawn_leaf_psi = -0.5_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('drought-decid: brief rain keeps the dry sum', state%dry_psi_sum, dry_before_rain, 0.0_wp)
      call check_true('drought-decid: brief rain keeps senescing', out%leaf_shed_potential > 0.5_wp * params%shed_rate_max)
      !----- (d) Rewetting: the wet sum crosses, the dry sum resets, flushing resumes. -----!
      do d = 1_ik, 30_ik
         env%predawn_leaf_psi = -0.5_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('drought-decid: rewet resets the dry sum', state%dry_psi_sum, 0.0_wp, 0.0_wp)
      call check_true('drought-decid: rewet => flush recovers',  out%leaf_flush_potential > 0.9_wp * params%flush_rate_max)
      call check_true('drought-decid: rewet => senescence stops', out%leaf_shed_potential < 0.05_wp * params%shed_rate_max)
   end subroutine test_drought_deciduous

   !----- 5. Light-driven leaf exchange: flush stays on, senescence rises with light. -------!
   subroutine test_light_exchanging()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 5. light-driven leaf exchange (flush={}, shed={LIGHT} on shortwave) --'
      params%flush_cue_mask       = CUE_NONE ; params%shed_cue_mask = CUE_LIGHT
      params%light_variable       = LIGHT_RADIATION
      params%shed_light_threshold = 200.0_wp
      params%shed_light_sharpness = 0.05_wp              ! > 0: bright light triggers senescence
      do d = 1_ik, 40_ik
         env%rad = 50.0_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('leaf-exch: running mean tracks the shortwave', state%shortwave_mean, 50.0_wp, 1.0_wp)
      call check_true('leaf-exch: dim => little senescence', out%leaf_shed_potential < 0.01_wp * params%shed_rate_max)
      call check_close('leaf-exch: flush potential = flush_rate_max', out%leaf_flush_potential, params%flush_rate_max, 1.0e-12_wp)
      do d = 1_ik, 40_ik
         env%rad = 500.0_wp
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_true('leaf-exch: bright => senescence', out%leaf_shed_potential > 0.9_wp * params%shed_rate_max)
      call check_close('leaf-exch: flush unchanged by light', out%leaf_flush_potential, params%flush_rate_max, 1.0e-12_wp)
   end subroutine test_light_exchanging

   !----- 6. Rate mapping: the potentials are exactly rate_max times the tendencies. ------!
   subroutine test_rate_mapping()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 6. rate mapping (potential = rate_max * tendency) --'
      params%flush_cue_mask = CUE_WATER ; params%shed_cue_mask = CUE_LIGHT
      params%light_variable = LIGHT_RADIATION ; params%shed_light_threshold = 250.0_wp
      params%shed_light_sharpness = 0.05_wp
      do d = 1_ik, 25_ik
         env%predawn_leaf_psi = -1.7_wp ; env%rad = 260.0_wp   ! partial flush + partial shed
         call phenology_kernel(env, params, 1.0_wp, state, out)
      end do
      call check_close('flush potential == flush_rate_max * flush tendency', out%leaf_flush_potential, &
                       params%flush_rate_max * state%leaf_flush_tendency, 1.0e-12_wp)
      call check_close('shed potential  == shed_rate_max  * shed tendency',  out%leaf_shed_potential,  &
                       params%shed_rate_max  * state%leaf_shed_tendency,  1.0e-12_wp)
      call check_true('tendencies are partial and in [0,1]',                                           &
                      state%leaf_flush_tendency > 0.0_wp .and. state%leaf_flush_tendency < 1.0_wp .and.  &
                      state%leaf_shed_tendency  > 0.0_wp .and. state%leaf_shed_tendency  < 1.0_wp)
   end subroutine test_rate_mapping

   !----- 7. Degenerate drivers: potentials stay finite and non-negative, no FP trap. -----!
   subroutine test_degenerate()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: ok
      print '(a)', '-- 7. degenerate drivers (every cue on both sides) --'
      params%flush_cue_mask = CUE_ALL ; params%shed_cue_mask = CUE_ALL
      params%light_variable = LIGHT_RADIATION
      ok = .true.
      do d = 1_ik, 60_ik
         env%doy = modulo(d - 1_ik, 365_ik) + 1_ik
         if (mod(d, 2_ik) == 0_ik) then
            env%temp_day = 350.0_wp ; env%predawn_leaf_psi = 5.0_wp ; env%daylength = 30.0_wp ; env%rad = 5000.0_wp
         else
            env%temp_day = 200.0_wp ; env%predawn_leaf_psi = -100.0_wp ; env%daylength = -5.0_wp ; env%rad = -50.0_wp
         end if
         call phenology_kernel(env, params, 1.0_wp, state, out)
         if (out%leaf_flush_potential < 0.0_wp .or. out%leaf_shed_potential < 0.0_wp) ok = .false.
         if (out%leaf_flush_potential /= out%leaf_flush_potential .or.                            &
             out%leaf_shed_potential  /= out%leaf_shed_potential) ok = .false.        ! NaN
      end do
      call check_true('degenerate: potentials finite and >= 0, no trap', ok)
   end subroutine test_degenerate

   !----- 8. Daylength: polar day/night + equator. ------------------------------------------!
   subroutine test_daylength_polar()
      print '(a)', '-- 8. daylength polar branches (meds_time) --'
      call check_true('polar day ~ 24 h',  daylength(80.0_wp, 172_ik) > 23.5_wp)
      call check_true('polar night ~ 0 h', daylength(80.0_wp, 355_ik) < 0.5_wp)
      call check_close('equator ~ 12 h',   daylength(0.0_wp, 172_ik), 12.0_wp, 0.5_wp)
   end subroutine test_daylength_polar

   !----- 9. leaf_turnover_step: the one place the tendencies become leaf loss. ------------!
   subroutine test_leaf_turnover_step()
      real(wp) :: sen, bg, cap
      real(wp), parameter :: TURN = 0.5_wp            ! [1/yr]
      print '(a)', '-- 9. leaf_turnover_step --'
      !----- Flushing, not senescing: background turnover only, and the flush cap. ----------!
      call leaf_turnover_step(1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.06_wp, 0.3_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('flushing: no senescence',             sen, 0.0_wp, 0.0_wp)
      call check_close('flushing: background = turnover/yr',  bg,  TURN / yr_day, 1.0e-15_wp)
      call check_close('flushing: flush cap = rate*full*dt',  cap, 0.06_wp, 1.0e-15_wp)
      !----- Background turnover follows the flush tendency: none in dormancy. --------------!
      call leaf_turnover_step(1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.06_wp, 0.3_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('dormant: no background turnover', bg,  0.0_wp, 0.0_wp)
      call check_close('dormant: no flush cap',           cap, 0.0_wp, 0.0_wp)
      !----- Senescence = shed_rate_max * tendency * pool * dt; background adds to it. -------!
      call leaf_turnover_step(0.5_wp, 1.0_wp, 0.5_wp, 1.0_wp, 0.06_wp, 0.1_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('senescing: senescence = rate*pool*dt', sen, 0.05_wp, 1.0e-15_wp)
      call check_close('senescing: background still acts',    bg,  TURN / yr_day * 0.5_wp * 0.5_wp, 1.0e-15_wp)
      !----- The leaf-cover floor: senescence stops at min_leaf_cover of the full canopy, ---!
      !      which is what makes a PFT evergreen (an emergent habit, not a flag).          --!
      call leaf_turnover_step(0.82_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.3_wp, 0.0_wp, 0.8_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('floor: senescence stops at min_leaf_cover', 0.82_wp - sen, 0.8_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.7_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.3_wp, 0.0_wp, 0.8_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('floor: none below min_leaf_cover', sen, 0.0_wp, 0.0_wp)
      !----- Dormant snap: a dormant canopy left below bare_leaf_cover goes bare. -----------!
      call leaf_turnover_step(0.021_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, TURN, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: a dormant canopy goes bare', sen + bg, 0.021_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.021_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: none while flushing', sen, 0.0021_wp, 1.0e-15_wp)
      call leaf_turnover_step(0.021_wp, 1.0_wp, 0.0_wp, 1.0_wp, 0.06_wp, 0.1_wp, 0.0_wp, 0.02_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('snap: none when the floor is not below bare', sen, 0.001_wp, 1.0e-15_wp)
      !----- A rate * dt above 1 removes at most the pool. -----------------------------------!
      call leaf_turnover_step(0.4_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.06_wp, 5.0_wp, 0.0_wp, 0.0_wp, 0.02_wp, &
                              1.0_wp, sen, bg, cap)
      call check_close('clamp: senescence removes at most the pool', sen, 0.4_wp, 1.0e-15_wp)
   end subroutine test_leaf_turnover_step

end program test_plant_phenology
