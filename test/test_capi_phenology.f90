!==========================================================================================!
! test_capi_phenology -- COVERAGE FOR THE PHENOLOGY C-API SHIM.                               !
!                                                                                          !
! WHY THIS TEST EXISTS. The same reason test_capi_leaf does, one subsystem over. A C-API shim  !
! that only the optional `-DMEDS_BUILD_PYLIB=ON` library compiles is invisible to ctest on both  !
! back ends, so a field inserted into `pheno_params_t` or a renamed component can leave the C     !
! API unable to compile while the whole suite stays green. That is not hypothetical: it happened  !
! to the leaf shim (#95 -> #100) and, while this file was being written, to the DEMOGRAPHY shim,   !
! which did not compile against main at all because a signature moved under it a day earlier.     !
!                                                                                          !
! Registering the shim in a ctest target COMPILES it in every build, which is the substance of    !
! the fix. The assertions then cover what a compile cannot: that the bind(c) mirrors carry values  !
! through to the kernel and back, rather than silently reading the wrong member.                    !
!                                                                                          !
! Needs no .so and no Python -- it calls the bind(c) procedure directly as Fortran.                 !
!==========================================================================================!
program test_capi_phenology
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use meds_kinds,           only : wp, ik
   use meds_capi_phenology,  only : pheno_env_c, pheno_params_c, pheno_state_c, pheno_out_c,      &
                                    meds_phenology_step
   use meds_phenology_types, only : CUE_NONE, CUE_TEMP
   use meds_test_support,    only : check, check_close, banner
   implicit none

   type(pheno_env_c)    :: env
   type(pheno_params_c) :: p
   type(pheno_state_c)  :: st, st0
   type(pheno_out_c)    :: out
   real(c_double), parameter :: DT = 1.0_c_double
   real(wp) :: gdd_gain

   call banner('phenology C-API shim')

   !=== 1. A PERMISSIVE (no-cue) canopy flushes and does not shed. ==========================!
   !      CUE_NONE is the evergreen fixed point, so this pins the plainest path through the   !
   !      shim: the masks reach the kernel and the two rates come back through pheno_out_c.   !
   env = warm_day()
   p   = base_params(CUE_NONE, CUE_NONE)
   st  = zero_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check(out%leaf_flush_rate > 0.0_wp, 'no-cue canopy flushes (flush rate reaches pheno_out_c)')
   call check_close(out%leaf_shed_rate, 0.0_wp, 1.0e-12_wp, 'no-cue canopy has no ACTIVE shed')

   !=== 2. The STATE struct is in/out: thermal memory must come back changed. ================!
   !      pheno_state_c is the one mirror the caller owns across steps, so a shim that passed  !
   !      it by value, or unpacked it into the wrong member, would look fine on a single call. !
   env = warm_day()
   p   = base_params(CUE_TEMP, CUE_TEMP)
   st  = zero_state()
   st0 = st
   call meds_phenology_step(env, p, DT, st, out)
   gdd_gain = st%gdd - st0%gdd
   call check(gdd_gain > 0.0_wp, 'a warm day accumulates GDD back into pheno_state_c')
   call check_close(gdd_gain, env%temp_day - p%gdd_base_temp, 1.0e-9_wp,                          &
                    'the GDD increment is (temp_day - gdd_base_temp), so BOTH fields mapped')

   !=== 3. A COLD day must not accumulate GDD -- the value reaches the kernel, not just the   !
   !       struct. A mirror that read a neighbouring member would still "work" in test 2.     !
   env = warm_day() ; env%temp_day = p%gdd_base_temp - 5.0_c_double
   st  = zero_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check_close(st%gdd, 0.0_wp, 1.0e-12_wp, 'a sub-base day accumulates no GDD')

   print '(a)', 'test_capi_phenology: ALL PASSED'

contains

   type(pheno_env_c) function warm_day() result(e)
      e%temp_day      = 20.0_c_double
      e%soil_temp     = 18.0_c_double
      e%avail_water   = 0.8_c_double
      e%dmax_leaf_psi = -0.2_c_double
      e%rad           = 400.0_c_double
      e%daylength     = 14.0_c_double
      e%doy           = 150_c_int
      e%hemis_north   = 1_c_int
   end function warm_day

   !----- Parameters with the two cue masks under test and everything else benign. -----------!
   type(pheno_params_c) function base_params(flush_mask, shed_mask) result(q)
      integer(ik), intent(in) :: flush_mask, shed_mask
      q%flush_cue_mask      = int(flush_mask, c_int)
      q%shed_cue_mask       = int(shed_mask,  c_int)
      q%cue_sharpness       = 1.0_c_double
      q%k_flush_max         = 0.05_c_double
      q%k_shed_max          = 0.05_c_double
      q%tau_flush           = 10.0_c_double
      q%tau_shed            = 10.0_c_double
      q%gdd_base_temp       = 5.0_c_double
      q%chill_base_temp     = 5.0_c_double
      q%phen_a              = -68.0_c_double
      q%phen_b              = 638.0_c_double
      q%phen_c              = -0.01_c_double
      q%cold_drop_daylength = 10.5_c_double
      q%cold_drop_soiltemp1 = 8.0_c_double
      q%cold_drop_soiltemp2 = 2.0_c_double
      q%water_use_potential = 0_c_int
      q%water_off_threshold = 0.2_c_double
      q%water_on_threshold  = 0.4_c_double
      q%water_window        = 10.0_c_double
      q%water_width         = 0.05_c_double
      q%leaf_psi_tlp        = -2.0_c_double
      q%low_psi_threshold   = 5.0_c_double
      q%high_psi_threshold  = 5.0_c_double
      q%photo_crit          = 11.0_c_double
      q%photo_slope         = 1.0_c_double
      q%light_on_threshold  = 100.0_c_double
      q%light_width         = 20.0_c_double
      q%light_window        = 10.0_c_double
   end function base_params

   type(pheno_state_c) function zero_state() result(s)
      s%flush_drive  = 0.0_c_double ; s%shed_drive    = 0.0_c_double
      s%gdd          = 0.0_c_double ; s%chill         = 0.0_c_double
      s%water_avg    = 0.5_c_double ; s%low_psi_days  = 0.0_c_double
      s%high_psi_days= 0.0_c_double ; s%light_avg     = 300.0_c_double
   end function zero_state

end program test_capi_phenology
