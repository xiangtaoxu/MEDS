!==========================================================================================!
! test_plant_phenology -- unit tests for the stateless leaf-phenology SIGNAL kernel.         !
!                                                                                          !
! The kernel emits TWO relative rate tendencies (leaf_flush_rate, leaf_shed_rate [1/day]) from  !
! two governor accumulators. The tests exercise the FOUR target patterns (design §1a) plus the   !
! rate-mapping boundary values, hysteresis, and FPE safety:                                       !
!                                                                                          !
!   1. EVERGREEN            : flush={} shed={} => flush_rate = k_flush_max, shed_rate = 0 always.  !
!   2. TEMPERATE DECIDUOUS  : flush={TEMP} shed={TEMP} => flush pulse in summer, shed pulse autumn. !
!   3. FACULTATIVE DROUGHT  : flush={HYDRO} shed={HYDRO} => sheds under sustained low dmax_leaf_psi, !
!                             flushes on rewet; stays flushed when never droughted.                 !
!   4. LIGHT LEAF-EXCHANGING: flush={} shed={LIGHT} => flush stays HIGH while shed RISES with light  !
!                             (canopy full while turning over).                                     !
!   5. RATE MAPPING         : leaf_flush_rate == k_flush_max*flush_drive; shed likewise (boundaries).!
!   6. DEGENERATE / FPE     : extreme drivers keep both rates finite and >= 0 (no trap).             !
!   7. DAYLENGTH            : the relocated meds_time daylength has the polar branches right.         !
!==========================================================================================!
program test_plant_phenology
   use meds_kinds,           only : wp, ik
   use meds_time,            only : daylength
   use meds_plant_interface, only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,   &
                                    update_phenology,                                          &
                                    CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO, CUE_LIGHT
   implicit none

   real(wp),    parameter :: twopi = 6.283185307179586_wp
   integer(ik) :: nfail
   nfail = 0_ik

   call test_evergreen()
   call test_temperate_deciduous()
   call test_drought_deciduous()
   call test_light_exchanging()
   call test_rate_mapping()
   call test_degenerate()
   call test_daylength_polar()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_phenology: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_phenology: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check_true(name, cond)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a)', '  FAIL : ', name
      end if
   end subroutine check_true

   subroutine check_close(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,f10.5,a,f10.5,a)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,f10.5,a,f10.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check_close

   !----- Annual forcing (northern hemisphere; summer solstice ~ doy 201). -----------------!
   pure real(wp) function annual_temp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 12.0_wp * cos(twopi * real(doy - 201_ik, wp) / 365.0_wp)
   end function annual_temp

   pure real(wp) function annual_soiltemp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 8.0_wp * cos(twopi * real(doy - 215_ik, wp) / 365.0_wp)
   end function annual_soiltemp

   !----- 1. Evergreen: flush stays at k_flush_max, no active shed, under ANY drivers. ------!
   subroutine test_evergreen()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: flush_ok, shed_ok
      print '(a)', '-- 1. evergreen (flush={}, shed={}) --'
      params%flush_cue_mask = CUE_NONE ; params%shed_cue_mask = CUE_NONE
      state    = pheno_state_t()
      flush_ok = .true. ; shed_ok = .true.
      do d = 1_ik, 400_ik
         env%doy           = modulo(d - 1_ik, 365_ik) + 1_ik
         env%temp_day      = annual_temp(env%doy)
         env%soil_temp     = annual_soiltemp(env%doy)
         env%avail_water   = 0.05_wp                     ! bone dry ...
         env%dmax_leaf_psi = -8.0_wp                     ! ... and cavitated ...
         env%rad           = 800.0_wp                    ! ... and blazing -- none of it is a cue
         env%daylength     = daylength(45.0_wp, env%doy)
         call update_phenology(env, params, 1.0_wp, state, out)
         if (abs(out%leaf_flush_rate - params%k_flush_max) > 1.0e-9_wp) flush_ok = .false.
         if (out%leaf_shed_rate /= 0.0_wp) shed_ok = .false.
      end do
      call check_true('evergreen: flush_rate = k_flush_max regardless of drivers', flush_ok)
      call check_true('evergreen: shed_rate = 0 regardless of drivers',            shed_ok)
   end subroutine test_evergreen

   !----- 2. Temperate deciduous: flush pulse in summer, shed pulse in autumn. --------------!
   subroutine test_temperate_deciduous()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d, doy
      real(wp)    :: fl_200, sh_200, fl_340, sh_340, gdd_200
      print '(a)', '-- 2. temperate deciduous (flush={TEMP}, shed={TEMP}) --'
      params%flush_cue_mask = CUE_TEMP ; params%shed_cue_mask = CUE_TEMP
      state           = pheno_state_t()
      env%hemis_north = .true.
      fl_200 = -9.0_wp ; sh_200 = -9.0_wp ; fl_340 = -9.0_wp ; sh_340 = -9.0_wp ; gdd_200 = 0.0_wp
      do d = 1_ik, 730_ik                                 ! two years (year 2 has prior-winter chill)
         doy           = modulo(d - 1_ik, 365_ik) + 1_ik
         env%doy       = doy
         env%temp_day  = annual_temp(doy)
         env%soil_temp = annual_soiltemp(doy)
         env%daylength = daylength(45.0_wp, doy)
         call update_phenology(env, params, 1.0_wp, state, out)
         if (d == 565_ik) then                            ! year-2 mid-summer (doy 200)
            fl_200 = state%flush_drive ; sh_200 = state%shed_drive ; gdd_200 = state%gdd
         end if
         if (d == 705_ik) then                            ! year-2 late autumn (doy 340)
            fl_340 = state%flush_drive ; sh_340 = state%shed_drive
         end if
      end do
      call check_true('deciduous: flush_drive HIGH in mid-summer',   fl_200 > 0.5_wp)
      call check_true('deciduous: shed_drive  LOW  in mid-summer',   sh_200 < 0.2_wp)
      call check_true('deciduous: shed_drive  HIGH in late autumn',  sh_340 > 0.4_wp)
      call check_true('deciduous: flush_drive LOW  in late autumn',  fl_340 < 0.5_wp)
      call check_true('deciduous: GDD accumulated by summer',        gdd_200 > 100.0_wp)
   end subroutine test_temperate_deciduous

   !----- 3. Facultative drought-deciduous: sheds on drought, flushes on rewet. -------------!
   subroutine test_drought_deciduous()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      real(wp)    :: shed_watered, shed_drought, flush_rewet
      print '(a)', '-- 3. facultative drought-deciduous (flush={HYDRO}, shed={HYDRO}) --'
      params%flush_cue_mask = CUE_HYDRO ; params%shed_cue_mask = CUE_HYDRO   ! tlp=-2, thresholds=10 d
      state = pheno_state_t()
      env%doy = 1_ik
      !----- (a) well-watered: sustained high dmax_leaf_psi -> no active shed (facultative). --!
      do d = 1_ik, 30_ik
         env%dmax_leaf_psi = -0.5_wp                      ! >= 0.5*tlp (-1): a wet day
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      shed_watered = out%leaf_shed_rate
      call check_true('drought-decid: watered => no active shed', shed_watered < 0.05_wp * params%k_shed_max)
      call check_true('drought-decid: watered => flush is on',    out%leaf_flush_rate > 0.5_wp * params%k_flush_max)
      !----- (b) drought: sustained low dmax_leaf_psi -> shed rises, flush falls. ------------!
      do d = 1_ik, 30_ik
         env%dmax_leaf_psi = -3.0_wp                      ! < tlp (-2): a dry day
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      shed_drought = out%leaf_shed_rate
      call check_true('drought-decid: drought => active shed rises', shed_drought > 0.5_wp * params%k_shed_max)
      call check_true('drought-decid: drought => flush falls',       out%leaf_flush_rate < 0.2_wp * params%k_flush_max)
      !----- (c) rewet: flush recovers. ----------------------------------------------------!
      do d = 1_ik, 30_ik
         env%dmax_leaf_psi = -0.5_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      flush_rewet = out%leaf_flush_rate
      call check_true('drought-decid: rewet => flush recovers',      flush_rewet > 0.5_wp * params%k_flush_max)
   end subroutine test_drought_deciduous

   !----- 4. Light-driven leaf-exchanging: flush stays HIGH, shed rises with light. ---------!
   subroutine test_light_exchanging()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      real(wp)    :: flush_lowlight, shed_lowlight, flush_highlight, shed_highlight
      print '(a)', '-- 4. light-driven leaf-exchanging (flush={}, shed={LIGHT}) --'
      params%flush_cue_mask = CUE_NONE ; params%shed_cue_mask = CUE_LIGHT   ! on=200, width=50, window=10
      state = pheno_state_t()
      !----- (a) low light: little active shed; flush stays at k_flush_max. -----------------!
      do d = 1_ik, 40_ik
         env%rad = 50.0_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      flush_lowlight = out%leaf_flush_rate ; shed_lowlight = out%leaf_shed_rate
      call check_true('leaf-exch: low light => shed small',           shed_lowlight < 0.1_wp * params%k_shed_max)
      call check_close('leaf-exch: flush = k_flush_max (permissive)', flush_lowlight, params%k_flush_max, 1.0e-9_wp)
      !----- (b) high light: shed rises; flush UNCHANGED (canopy stays full while exchanging). !
      do d = 1_ik, 40_ik
         env%rad = 500.0_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      flush_highlight = out%leaf_flush_rate ; shed_highlight = out%leaf_shed_rate
      call check_true('leaf-exch: high light => active shed rises',   shed_highlight > 0.5_wp * params%k_shed_max)
      call check_close('leaf-exch: flush unchanged by light',         flush_highlight, params%k_flush_max, 1.0e-9_wp)
      call check_true('leaf-exch: BOTH rates > 0 under high light',   flush_highlight > 0.0_wp .and. shed_highlight > 0.0_wp)
   end subroutine test_light_exchanging

   !----- 5. Rate mapping: the outputs are exactly k_*_max times the governor drives. -------!
   subroutine test_rate_mapping()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 5. rate mapping (leaf_rate = k*drive) --'
      params%flush_cue_mask = CUE_HYDRO ; params%shed_cue_mask = CUE_LIGHT
      state = pheno_state_t()
      env%doy = 1_ik
      do d = 1_ik, 25_ik
         env%dmax_leaf_psi = -1.5_wp ; env%rad = 260.0_wp     ! partial flush + partial shed
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_close('flush_rate == k_flush_max * flush_drive', out%leaf_flush_rate, &
                       params%k_flush_max * state%flush_drive, 1.0e-12_wp)
      call check_close('shed_rate  == k_shed_max  * shed_drive',  out%leaf_shed_rate,  &
                       params%k_shed_max  * state%shed_drive,  1.0e-12_wp)
      call check_true ('drives are in [0,1]', state%flush_drive >= 0.0_wp .and. state%flush_drive <= 1.0_wp &
                       .and. state%shed_drive >= 0.0_wp .and. state%shed_drive <= 1.0_wp)
   end subroutine test_rate_mapping

   !----- 6. Degenerate drivers: rates stay finite and non-negative, no FP trap. ------------!
   subroutine test_degenerate()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: ok
      print '(a)', '-- 6. degenerate drivers (all cues both masks) --'
      params%flush_cue_mask = ior(ior(CUE_TEMP, CUE_WATER), ior(CUE_HYDRO, CUE_PHOTO))
      params%shed_cue_mask  = ior(ior(CUE_TEMP, CUE_WATER), ior(CUE_HYDRO, CUE_LIGHT))
      state = pheno_state_t()
      ok    = .true.
      do d = 1_ik, 60_ik
         env%doy = modulo(d - 1_ik, 365_ik) + 1_ik
         if (mod(d, 2_ik) == 0_ik) then
            env%temp_day = 350.0_wp ; env%soil_temp = 330.0_wp ; env%avail_water = 10.0_wp
            env%dmax_leaf_psi = 5.0_wp ; env%daylength = 30.0_wp ; env%rad = 5000.0_wp
         else
            env%temp_day = 200.0_wp ; env%soil_temp = 210.0_wp ; env%avail_water = -5.0_wp
            env%dmax_leaf_psi = -100.0_wp ; env%daylength = -5.0_wp ; env%rad = -50.0_wp
         end if
         call update_phenology(env, params, 1.0_wp, state, out)
         if (out%leaf_flush_rate < 0.0_wp .or. out%leaf_shed_rate < 0.0_wp)          ok = .false.
         if (out%leaf_flush_rate /= out%leaf_flush_rate .or. out%leaf_shed_rate /= out%leaf_shed_rate) ok = .false.  ! NaN
      end do
      call check_true('degenerate: rates finite and >= 0, no trap', ok)
   end subroutine test_degenerate

   !----- 7. Daylength: polar day/night + equator (relocated to meds_time; ED2 branch fixed). -!
   subroutine test_daylength_polar()
      real(wp) :: dl_summer, dl_winter, dl_eq
      print '(a)', '-- 7. daylength polar branches (meds_time) --'
      dl_summer = daylength(80.0_wp, 172_ik)             ! high N latitude, near summer solstice
      dl_winter = daylength(80.0_wp, 355_ik)             ! near winter solstice
      dl_eq     = daylength(0.0_wp, 172_ik)              ! equator
      call check_true('polar day ~ 24 h (ED2 bug fixed)', dl_summer > 23.5_wp)
      call check_true('polar night ~ 0 h', dl_winter < 0.5_wp)
      call check_close('equator ~ 12 h', dl_eq, 12.0_wp, 0.5_wp)
   end subroutine test_daylength_polar

end program test_plant_phenology
