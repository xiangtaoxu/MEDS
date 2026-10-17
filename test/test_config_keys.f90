! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_config_keys -- a config may hold only the keys its reference lists.                  !
!                                                                                          !
! The check (meds_config_keys) reports a retired key with what replaced it, and any other     !
! unknown key with the listed key it most likely meant: a misspelling, or a key that moved     !
! section. A key the reference lists only commented out, at its default, is a listed key.     !
!==========================================================================================!
program test_config_keys
   use meds_kinds,        only : ik
   use meds_toml,         only : toml_table_t, toml_parse_file
   use meds_config_keys,  only : key_report_t, check_config_keys, nearest_key, MAIN_KEYS, PFT_KEYS
   use meds_test_support, only : banner, check
   implicit none
   type(toml_table_t) :: t
   type(key_report_t) :: report
   logical :: ok
   integer :: u

   call banner('config_keys')

   !----- The suggestion: the same name in another section first, then the closest spelling. -!
   call check(nearest_key('state.cohort_max', MAIN_KEYS) == 'output.cohort_max',                    &
              'a key that moved section is matched by its name')
   call check(nearest_key('forcing.timestpe', MAIN_KEYS) == 'forcing.timestep',                    &
              'a misspelt key is matched by its spelling')
   call check(len_trim(nearest_key('zzzz.qqqq', MAIN_KEYS)) == 0,                                   &
              'nothing close gives no suggestion')

   !----- A config with one key of each kind. -----!
   open(newunit=u, file='test_config_keys.toml', status='replace', action='write')
   write(u,'(a)') '[run]'
   write(u,'(a)') 'dt_slow   = "1d"'                ! set in the reference
   write(u,'(a)') 'n_threads = 4'                   ! listed there commented out, at its default
   write(u,'(a)') '[forcing]'
   write(u,'(a)') 'timestpe  = "3600s"'             ! misspelt
   write(u,'(a)') '[state]'
   write(u,'(a)') 'cohort_max = 2048'               ! moved to [output]
   write(u,'(a)') '[output]'
   write(u,'(a)') 'strict_caps = true'              ! retired
   close(u)
   call toml_parse_file('test_config_keys.toml', t, ok)
   call check(ok, 'the test config parses')
   call check_config_keys(t, MAIN_KEYS, report)
   call check(report%n == 3_ik, 'three keys are reported, and the two listed ones are not')
   call check(index(report%line(1), 'forcing.timestpe') > 0 .and.                                  &
              index(report%line(1), 'did you mean forcing.timestep?') > 0, 'the misspelling is named')
   call check(index(report%line(2), 'did you mean output.cohort_max?') > 0, 'the moved key is named')
   call check(index(report%line(3), 'output.strict_caps') > 0 .and. index(report%line(3), 'retired') > 0 &
              .and. index(report%line(3), 'raise the cap') > 0, 'the retired key says what replaced it')

   !----- The PFT file is checked against its own list: a main-file key is unknown there. -----!
   report = key_report_t()
   call check_config_keys(t, PFT_KEYS, report)
   call check(report%n == 5_ik, 'every main-file key is unknown to the PFT reference')

   !----- Keys whose feature is gone are refused with what replaced it: the inert above-ground  !
   !      litter shares and the nitrogen keys of a cycle that is not built. -------------------!
   open(newunit=u, file='test_config_keys.toml', status='replace', action='write')
   write(u,'(a)') '[soil_carbon]'
   write(u,'(a)') 'n_cycle_on = true'
   write(u,'(a)') 'agf_fast   = 0.5'
   close(u)
   call toml_parse_file('test_config_keys.toml', t, ok)
   call check(ok, 'the retired-key config parses')
   report = key_report_t()
   call check_config_keys(t, MAIN_KEYS, report)
   call check(report%n == 2_ik, 'both retired soil-carbon keys are reported')
   call check(index(report%line(1), 'soil_carbon.n_cycle_on') > 0 .and. index(report%line(1), 'retired') > 0 &
              .and. index(report%line(1), 'ED-MORE') > 0, 'the nitrogen key says the cycle is not built')
   call check(index(report%line(2), 'soil_carbon.agf_fast') > 0 .and. index(report%line(2), 'retired') > 0 &
              .and. index(report%line(2), 'f_labile_leaf') > 0, 'the litter share says the litter splits by tissue')

   !----- In the PFT file: the crown-area allometry and the crown fraction, gone with the       !
   !      closed crown. ----------------------------------------------------------------------!
   open(newunit=u, file='test_config_keys.toml', status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'crown_area_frac = [1.0]'
   write(u,'(a)') '[allometry]'
   write(u,'(a)') 'ca_b1 = 0.370'
   close(u)
   call toml_parse_file('test_config_keys.toml', t, ok)
   call check(ok, 'the retired PFT-key config parses')
   report = key_report_t()
   call check_config_keys(t, PFT_KEYS, report)
   call check(report%n == 2_ik, 'both retired crown-area keys are reported')
   call check(index(report%line(1), 'pft.crown_area_frac') > 0 .and. index(report%line(1), 'retired') > 0 &
              .and. index(report%line(1), 'closed-crown') > 0, 'the crown fraction says the crowns close')
   call check(index(report%line(2), 'allometry.ca_b1') > 0 .and. index(report%line(2), 'retired') > 0  &
              .and. index(report%line(2), 'closed-crown') > 0, 'the crown-area allometry says the crowns close')

   open(newunit=u, file='test_config_keys.toml', status='old')
   close(u, status='delete')
   write(*,'(a)') 'test_config_keys: ALL PASSED'
end program test_config_keys
