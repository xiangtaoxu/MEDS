!==========================================================================================!
! test_plant_phenology -- unit tests for the stateless leaf-phenology library.              !
!                                                                                          !
!   1. EVERGREEN      : CUE_NONE stays PHEN_ON under any drivers.                             !
!   2. COLD-DECIDUOUS : a seasonal temperature cycle gives ON in summer, OFF in winter.        !
!   3. DROUGHT        : the running-mean water ramp gives ON when wet, OFF when dry.            !
!   4. HYDRAULIC      : sustained high leaf-psi -> ON; sustained low leaf-psi -> OFF.            !
!   5. DEADBAND       : a marginal driver sits in DORMANT between OFF and ON (hysteresis).        !
!   6. MULTI-CUE      : the most-limiting cue governs (a summer drought forces OFF; WATER limits). !
!   7. DEGENERATE     : extreme drivers keep the status in {ON,OFF,DORMANT} with no FPE.            !
!   8. DAYLENGTH      : the ported daylength has the polar-day/night branches right (ED2 bug fixed). !
!==========================================================================================!
program test_plant_phenology
   use meds_kinds,           only : wp, ik
   use meds_plant_phenology, only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,   &
                                    update_phenology, daylength,                                &
                                    CUE_NONE, CUE_TEMP, CUE_WATER, CUE_HYDRO, CUE_PHOTO,        &
                                    PHEN_ON, PHEN_DORMANT, PHEN_OFF
   implicit none

   real(wp),    parameter :: twopi = 6.283185307179586_wp
   integer(ik) :: nfail
   nfail = 0_ik

   call test_evergreen()
   call test_cold_deciduous()
   call test_drought()
   call test_hydraulic()
   call test_deadband()
   call test_multicue()
   call test_degenerate()
   call test_daylength_polar()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_phenology: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_phenology: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- Assertions. ----------------------------------------------------------------------!
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

   subroutine check_int(name, got, expect)
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: got, expect
      if (got == expect) then
         print '(a,a,a,i0,a,i0,a)', '  ok   : ', name, '  (', got, ' == ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,i0,a,i0)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check_int

   subroutine check_real(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,f8.3,a,f8.3)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,f8.3,a,f8.3)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check_real

   !----- Annual forcing (northern hemisphere; summer solstice ~ doy 201). -----------------!
   pure real(wp) function annual_temp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 12.0_wp * cos(twopi * real(doy - 201_ik, wp) / 365.0_wp)
   end function annual_temp

   pure real(wp) function annual_soiltemp(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 8.0_wp * cos(twopi * real(doy - 215_ik, wp) / 365.0_wp)
   end function annual_soiltemp

   !----- 1. Evergreen: perpetually ON. ----------------------------------------------------!
   subroutine test_evergreen()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      logical     :: all_on
      print '(a)', '-- 1. evergreen (CUE_NONE) --'
      params%cue_mask = CUE_NONE
      state  = pheno_state_t()
      all_on = .true.
      do d = 1_ik, 400_ik
         env%doy         = modulo(d - 1_ik, 365_ik) + 1_ik
         env%temp_day    = annual_temp(env%doy)
         env%soil_temp   = annual_soiltemp(env%doy)
         env%avail_water = 0.05_wp                        ! bone dry ...
         env%psi_leaf    = -8.0_wp                        ! ... and cavitated ...
         env%daylength   = daylength(45.0_wp, env%doy)
         call update_phenology(env, params, 1.0_wp, state, out)
         if (out%phenology_status /= PHEN_ON) all_on = .false.
      end do
      call check_true('evergreen stays ON regardless of drivers', all_on)
   end subroutine test_evergreen

   !----- 2. Cold-deciduous: ON in summer, OFF in winter. ----------------------------------!
   subroutine test_cold_deciduous()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d, doy, status_doy(365), n_on_summer, n_off_winter, i
      print '(a)', '-- 2. cold-deciduous (CUE_TEMP) --'
      params%cue_mask = CUE_TEMP
      state           = pheno_state_t()
      status_doy      = PHEN_DORMANT
      env%hemis_north = .true.
      do d = 1_ik, 730_ik                                 ! two years (year 2 has prior-winter chill)
         doy           = modulo(d - 1_ik, 365_ik) + 1_ik
         env%doy       = doy
         env%temp_day  = annual_temp(doy)
         env%soil_temp = annual_soiltemp(doy)
         env%daylength = daylength(45.0_wp, doy)
         call update_phenology(env, params, 1.0_wp, state, out)
         status_doy(doy) = out%phenology_status           ! last write per doy = year 2
      end do
      call check_int('midsummer (doy 200) is ON', status_doy(200), PHEN_ON)
      call check_int('midwinter (doy 20) is OFF', status_doy(20), PHEN_OFF)
      n_on_summer  = 0_ik
      n_off_winter = 0_ik
      do i = 150_ik, 240_ik
         if (status_doy(i) == PHEN_ON) n_on_summer = n_on_summer + 1_ik
      end do
      do i = 1_ik, 45_ik
         if (status_doy(i) == PHEN_OFF) n_off_winter = n_off_winter + 1_ik
      end do
      call check_true('most of summer (doy 150-240) is ON', n_on_summer >= 60_ik)
      call check_true('deep winter (doy 1-45) is OFF', n_off_winter >= 30_ik)
   end subroutine test_cold_deciduous

   !----- 3. Drought: ON when wet, OFF when dry, ON again on recovery. ----------------------!
   subroutine test_drought()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 3. drought (CUE_WATER) --'
      params%cue_mask = CUE_WATER
      state           = pheno_state_t()
      state%water_avg = 0.8_wp
      env%doy         = 1_ik
      do d = 1_ik, 40_ik
         env%avail_water = 0.8_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('sustained wet -> ON', out%phenology_status, PHEN_ON)
      do d = 1_ik, 60_ik
         env%avail_water = 0.02_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('sustained dry -> OFF', out%phenology_status, PHEN_OFF)
      do d = 1_ik, 60_ik
         env%avail_water = 0.9_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('re-wet -> ON', out%phenology_status, PHEN_ON)
   end subroutine test_drought

   !----- 4. Hydraulic: sustained wet/dry leaf-psi. ----------------------------------------!
   subroutine test_hydraulic()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 4. hydraulic (CUE_HYDRO) --'
      params%cue_mask = CUE_HYDRO             ! psi_tlp = -2 MPa, low/high thresholds = 10 days
      state           = pheno_state_t()
      env%doy         = 1_ik
      do d = 1_ik, 20_ik
         env%psi_leaf = -0.5_wp               ! well-hydrated (>= 0.5*psi_tlp)
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('sustained high psi -> ON', out%phenology_status, PHEN_ON)
      do d = 1_ik, 20_ik
         env%psi_leaf = -3.0_wp               ! past the turgor-loss point
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('sustained low psi -> OFF', out%phenology_status, PHEN_OFF)
   end subroutine test_hydraulic

   !----- 5. Deadband: a marginal driver holds DORMANT between OFF and ON. ------------------!
   subroutine test_deadband()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 5. deadband hysteresis (CUE_WATER) --'
      params%cue_mask = CUE_WATER             ! off = 0.2, on = 0.5 => mid 0.35 is DORMANT
      state           = pheno_state_t()
      state%water_avg = 0.02_wp
      env%doy         = 1_ik
      do d = 1_ik, 40_ik
         env%avail_water = 0.02_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('dry start -> OFF', out%phenology_status, PHEN_OFF)
      do d = 1_ik, 60_ik
         env%avail_water = 0.35_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('marginal water -> DORMANT (deadband)', out%phenology_status, PHEN_DORMANT)
      do d = 1_ik, 60_ik
         env%avail_water = 0.9_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('re-wet -> ON', out%phenology_status, PHEN_ON)
   end subroutine test_deadband

   !----- 6. Multi-cue: the most-limiting cue governs. -------------------------------------!
   subroutine test_multicue()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d
      print '(a)', '-- 6. multi-cue limiting factor (CUE_TEMP | CUE_WATER) --'
      params%cue_mask = ior(CUE_TEMP, CUE_WATER)
      state           = pheno_state_t()
      env%hemis_north = .true.
      env%doy         = 200_ik                          ! warm, long-day summer
      env%temp_day    = annual_temp(200_ik)
      env%soil_temp   = annual_soiltemp(200_ik)
      env%daylength   = daylength(45.0_wp, 200_ik)
      do d = 1_ik, 80_ik                                 ! thermally favorable, but keep it dry
         env%avail_water = 0.02_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('warm summer but dry -> OFF', out%phenology_status, PHEN_OFF)
      call check_int('  governing cue is WATER', out%cue_limiting, CUE_WATER)
      do d = 1_ik, 60_ik                                 ! now also wet
         env%avail_water = 0.9_wp
         call update_phenology(env, params, 1.0_wp, state, out)
      end do
      call check_int('warm summer and wet -> ON', out%phenology_status, PHEN_ON)
   end subroutine test_multicue

   !----- 7. Degenerate drivers: status stays valid, no floating-point trap. ----------------!
   subroutine test_degenerate()
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: d, s
      logical     :: ok
      print '(a)', '-- 7. degenerate drivers (all cues) --'
      params%cue_mask = ior(ior(CUE_TEMP, CUE_WATER), ior(CUE_HYDRO, CUE_PHOTO))
      state           = pheno_state_t()
      ok              = .true.
      do d = 1_ik, 50_ik
         env%doy = modulo(d - 1_ik, 365_ik) + 1_ik
         if (mod(d, 2_ik) == 0_ik) then
            env%temp_day = 350.0_wp ; env%soil_temp = 330.0_wp ; env%avail_water = 10.0_wp
            env%psi_leaf = 5.0_wp   ; env%daylength = 30.0_wp
         else
            env%temp_day = 200.0_wp ; env%soil_temp = 210.0_wp ; env%avail_water = -5.0_wp
            env%psi_leaf = -100.0_wp ; env%daylength = -5.0_wp
         end if
         call update_phenology(env, params, 1.0_wp, state, out)
         s = out%phenology_status
         if (s /= PHEN_ON .and. s /= PHEN_OFF .and. s /= PHEN_DORMANT) ok = .false.
      end do
      call check_true('extreme drivers -> status in {ON,OFF,DORMANT}, no trap', ok)
   end subroutine test_degenerate

   !----- 8. Daylength: polar day/night + equator (the ED2 polar branch is fixed). ----------!
   subroutine test_daylength_polar()
      real(wp) :: dl_summer, dl_winter, dl_eq
      print '(a)', '-- 8. daylength polar branches --'
      dl_summer = daylength(80.0_wp, 172_ik)             ! high N latitude, near summer solstice
      dl_winter = daylength(80.0_wp, 355_ik)             ! near winter solstice
      dl_eq     = daylength(0.0_wp, 172_ik)              ! equator
      call check_true('polar day ~ 24 h (ED2 bug fixed)', dl_summer > 23.5_wp)
      call check_true('polar night ~ 0 h', dl_winter < 0.5_wp)
      call check_real('equator ~ 12 h', dl_eq, 12.0_wp, 0.5_wp)
   end subroutine test_daylength_polar

end program test_plant_phenology
