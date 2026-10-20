! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_config_keys -- a config may hold only the keys its reference lists, and only the    !
! keys its settings use.                                                                   !
!                                                                                          !
! check_config_keys reports an unknown key with the listed key it most likely meant: a     !
! misspelling, or a key that moved section. A key the reference lists only commented out,  !
! at its default, is a listed key. check_key_use reports a listed key the loaded settings  !
! give nothing to do.                                                                      !
!==========================================================================================!
program test_config_keys
   use meds_kinds,          only : ik
   use meds_toml,           only : toml_table_t, toml_parse_file
   use meds_config,         only : meds_config_t
   use meds_forcing_config, only : MET_BACKEND_ED_ERA5LAND
   use meds_config_check,   only : key_report_t, check_config_keys, check_key_use, nearest_key,     &
                                   MAIN_KEYS, PFT_KEYS
   use meds_test_support,   only : banner, check, build_test_config
   implicit none
   type(toml_table_t)  :: t, tp
   type(key_report_t)  :: report
   type(meds_config_t) :: cfg
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
   write(u,'(a)') 'strict_caps = true'              ! not a key
   close(u)
   call toml_parse_file('test_config_keys.toml', t, ok)
   call check(ok, 'the test config parses')
   call check_config_keys(t, MAIN_KEYS, report)
   call check(report%n == 3_ik, 'three keys are reported, and the two listed ones are not')
   call check(index(report%line(1), 'forcing.timestpe') > 0 .and.                                  &
              index(report%line(1), 'did you mean forcing.timestep?') > 0, 'the misspelling is named')
   call check(index(report%line(2), 'did you mean output.cohort_max?') > 0, 'the moved key is named')
   call check(index(report%line(3), 'output.strict_caps') > 0 .and.                                &
              index(report%line(3), 'is not a key MEDS reads') > 0, 'an unknown key is named')

   !----- The PFT file is checked against its own list: a main-file key is unknown there. -----!
   report = key_report_t()
   call check_config_keys(t, PFT_KEYS, report)
   call check(report%n == 5_ik, 'every main-file key is unknown to the PFT reference')

   !----- A listed key the settings give nothing to do: the forcing file's name with the archive, !
   !      and a [region] block in a site run. --------------------------------------------------!
   open(newunit=u, file='test_config_keys.toml', status='replace', action='write')
   write(u,'(a)') '[forcing]'
   write(u,'(a)') 'path = "met.nc"'
   write(u,'(a)') '[region]'
   write(u,'(a)') 'land_fraction_min = 0.5'
   close(u)
   call toml_parse_file('test_config_keys.toml', t, ok)
   open(newunit=u, file='test_config_keys_pft.toml', status='replace', action='write')
   write(u,'(a)') '[pft]'
   close(u)
   call toml_parse_file('test_config_keys_pft.toml', tp, ok)
   cfg = build_test_config()
   cfg%forcing%forcing_on = .true.
   cfg%forcing%backend    = MET_BACKEND_ED_ERA5LAND
   report = key_report_t()
   call check_key_use(t, tp, cfg, report)
   call check(report%n == 2_ik, 'both inapplicable keys are reported')
   call check(index(report%line(1), '[region]') > 0 .and. index(report%line(1), 'mode = "region"') > 0, &
              'a [region] block in a site run is named')
   call check(index(report%line(2), 'forcing.path') > 0 .and. index(report%line(2), 'ED_ERA5land') > 0, &
              'the file key under the archive is named')

   open(newunit=u, file='test_config_keys.toml', status='old')
   close(u, status='delete')
   open(newunit=u, file='test_config_keys_pft.toml', status='old')
   close(u, status='delete')
   write(*,'(a)') 'test_config_keys: ALL PASSED'
end program test_config_keys
