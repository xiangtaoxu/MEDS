! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_c_api_phenology -- COVERAGE FOR THE PHENOLOGY C-API SHIM.                              !
!                                                                                          !
! WHY THIS TEST EXISTS. The same reason test_c_api_leaf does, one subsystem over. A C-API shim !
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
! Needs no .so and no Python -- it calls the bind(c) procedures directly as Fortran.                !
!==========================================================================================!
program test_c_api_phenology
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : yr_day
   use meds_c_api_phenology, only : pheno_env_c, pheno_params_c, pheno_state_c, pheno_out_c,      &
                                    meds_phenology_step, meds_leaf_turnover_step, meds_daylength
   use meds_phenology_types, only : CUE_NONE, CUE_TEMP, CUE_WATER
   use meds_test_support, only : banner, check, check_close
   implicit none

   type(pheno_env_c)    :: env
   type(pheno_params_c) :: p
   type(pheno_state_c)  :: st
   type(pheno_out_c)    :: out
   real(c_double), parameter :: DT = 1.0_c_double
   real(c_double) :: sen, bg, cap

   call banner('phenology C-API shim')

   !=== 1. A no-cue canopy flushes at flush_rate_max and does not senesce. =================!
   !      This pins the plainest path through the shim: the masks and the rate maxima reach  !
   !      the kernel and the two potentials come back through pheno_out_c.                   !
   env = warm_day()
   p   = base_params(CUE_NONE, CUE_NONE)
   st  = birth_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check_close(out%leaf_flush_potential, p%flush_rate_max, 1.0e-12_wp,                         &
                    'no-cue canopy flushes at flush_rate_max (reaches pheno_out_c)')
   call check_close(out%leaf_shed_potential, 0.0_wp, 1.0e-12_wp, 'no-cue canopy does not senesce')

   !=== 2. The STATE struct is in/out: the warmth sum must come back changed by exactly      !
   !       (temp_day - flush_base_temp), so BOTH fields mapped. A shim that passed the state    !
   !       by value, or unpacked it into the wrong member, would look fine on a single call.   !
   p  = base_params(CUE_TEMP, CUE_TEMP)
   st = birth_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check_close(st%growing_degree_days, env%temp_day - p%flush_base_temp, 1.0e-9_wp,            &
                    'the warmth sum grows by (temp_day - flush_base_temp)')

   !=== 3. A COLD day must not add warmth -- the value reaches the kernel, not just the      !
   !       struct. A mirror that read a neighbouring member would still "work" in test 2.     !
   env%temp_day = p%flush_base_temp - 5.0_c_double
   st = birth_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check_close(st%growing_degree_days, 0.0_wp, 1.0e-12_wp, 'a sub-base day adds no warmth')

   !=== 4. The water fields: the dry sum grows by (leaf_psi_tlp - predawn_leaf_psi). ========!
   env = warm_day() ; env%predawn_leaf_psi = -3.5_c_double
   p   = base_params(CUE_WATER, CUE_WATER)
   st  = birth_state()
   call meds_phenology_step(env, p, DT, st, out)
   call check_close(st%dry_psi_sum, p%leaf_psi_tlp - env%predawn_leaf_psi, 1.0e-12_wp,               &
                    'the dry sum grows by (leaf_psi_tlp - predawn_leaf_psi)')
   call check_close(st%wet_psi_sum, 0.0_wp, 0.0_wp, 'a dry day adds no wet credit')

   !=== 5. meds_leaf_turnover_step carries every argument through (background turnover,     !
   !       senescence, the flush cap), each a distinct product of its inputs.                !
   call meds_leaf_turnover_step(0.5_c_double, 1.0_c_double, 0.5_c_double, 1.0_c_double,              &
                                0.06_c_double, 0.1_c_double, 0.5_c_double, 0.0_c_double,              &
                                0.02_c_double, 1.0_c_double, sen, bg, cap)
   call check_close(sen, 0.05_wp, 1.0e-15_wp, 'leaf turnover: senescence = shed_rate_max*tendency*leaf*dt')
   call check_close(bg, 0.5_wp / yr_day * 0.5_wp * 0.5_wp, 1.0e-15_wp,                              &
                    'leaf turnover: background = rate/yr * flush tendency * leaf * dt')
   call check_close(cap, 0.03_wp, 1.0e-15_wp, 'leaf turnover: flush cap = flush_rate_max*tendency*full*dt')

   !=== 6. meds_daylength passes latitude and day of year through: 12 h at the equator,      !
   !       a long day at 60 N in June, polar night at 80 N in December.                      !
   call check_close(meds_daylength(0.0_c_double, 172_c_int), 12.0_wp, 0.5_wp, 'daylength: equator ~12 h')
   call check(meds_daylength(60.0_c_double, 172_c_int) > 18.0_wp, 'daylength: 60 N in June > 18 h')
   call check(meds_daylength(80.0_c_double, 355_c_int) < 0.5_wp, 'daylength: 80 N in December ~0 h')

   print '(a)', 'test_c_api_phenology: ALL PASSED'

contains

   type(pheno_env_c) function warm_day() result(e)
      e%temp_day         = 293.15_c_double
      e%daylength        = 14.0_c_double
      e%rad              = 400.0_c_double
      e%predawn_leaf_psi = -0.2_c_double
      e%doy              = 150_c_int
      e%hemis_north      = 1_c_int
   end function warm_day

   !----- Parameters with the two cue masks under test and everything else benign. -----------!
   type(pheno_params_c) function base_params(flush_mask, shed_mask) result(q)
      integer(ik), intent(in) :: flush_mask, shed_mask
      q%flush_cue_mask        = int(flush_mask, c_int)
      q%shed_cue_mask         = int(shed_mask,  c_int)
      q%flush_cue_timescale   = 5.0_c_double
      q%shed_cue_timescale    = 5.0_c_double
      q%flush_rate_max        = 0.05_c_double
      q%shed_rate_max         = 0.04_c_double
      q%flush_base_temp       = 278.15_c_double
      q%flush_degree_days     = 100.0_c_double
      q%flush_temp_sharpness  = 0.04_c_double
      q%shed_base_temp        = 290.15_c_double
      q%shed_degree_days      = 50.0_c_double
      q%shed_temp_sharpness   = 0.1_c_double
      q%light_variable        = 1_c_int
      q%flush_light_threshold = 12.0_c_double
      q%flush_light_sharpness = 1.0_c_double
      q%shed_light_threshold  = 11.0_c_double
      q%shed_light_sharpness  = -1.0_c_double
      q%light_window          = 10.0_c_double
      q%leaf_psi_tlp          = -2.0_c_double
      q%flush_water_sum       = 10.0_c_double
      q%flush_water_sharpness = 0.5_c_double
      q%shed_water_sum        = 10.0_c_double
      q%shed_water_sharpness  = 0.5_c_double
   end function base_params

   type(pheno_state_c) function birth_state() result(s)
      s%leaf_flush_tendency = 1.0_c_double ; s%leaf_shed_tendency = 0.0_c_double
      s%growing_degree_days = 0.0_c_double ; s%cold_degree_days   = 0.0_c_double
      s%wet_psi_sum         = 0.0_c_double ; s%dry_psi_sum        = 0.0_c_double
      s%shortwave_mean      = 0.0_c_double
   end function birth_state

end program test_c_api_phenology
