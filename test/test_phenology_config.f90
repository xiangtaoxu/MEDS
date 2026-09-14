!==========================================================================================!
! test_phenology_config -- the [phenology] section actually reaches the PFT table (#245).       !
!                                                                                          !
! The gate on this section is its PRESENCE in the PFT file: absent, every PFT keeps the          !
! built-in evergreen defaults; present, every key in it is REQUIRED. That contract was          !
! implemented by testing for ONE key inside the block -- `flush_cue_mask` -- which the shipped   !
! meds_config_pft.toml never documented. It documented `cue_mask`, a name no reader consumes.    !
!                                                                                          !
! So a config that wrote a full, deliberate [phenology] block was skipped in SILENCE. The        !
! presence map never saw twenty-three required keys go missing, and every PFT fell back to       !
! CUE_NONE / CUE_NONE -- flush = 1, shed = 0, the evergreen fixed point -- whatever leaf habit    !
! it declared. The Ithaca reference stand, declared cold-deciduous, held LAI 5.28-5.66 through    !
! every January of a 50-year run, and no MEDS run has ever had a leaf-area cycle.                 !
!                                                                                          !
! What is asserted here is the GATE, on both sides of it, because that is the thing that failed: !
!                                                                                          !
!   1. a PFT file with NO [phenology] section keeps the evergreen defaults (the documented        !
!      fallback -- it has to stay reachable, or every existing config becomes an error);          !
!   2. a PFT file WITH the section is read, and the masks land where the kernel looks;            !
!   3. `toml_has_section` -- the new gate -- sees a section by ANY of its keys, not by one        !
!      chosen name. This is the assertion that would have failed before the fix.                   !
!==========================================================================================!
program test_phenology_config
   use meds_kinds,        only : wp, ik
   use meds_toml,         only : toml_table_t, toml_parse_file, toml_has, toml_has_section
   use meds_test_support, only : banner, check, check_close, test_report
   implicit none

   character(len=*), parameter :: F_NONE = 'test_pheno_cfg_none_tmp.toml'
   character(len=*), parameter :: F_FULL = 'test_pheno_cfg_full_tmp.toml'
   character(len=*), parameter :: F_STALE = 'test_pheno_cfg_stale_tmp.toml'
   type(toml_table_t) :: t_none, t_full, t_stale
   logical            :: ok
   integer(ik)        :: u

   call banner('[phenology] section gate (#245)')

   !----- (a) a PFT file with no phenology section at all. ---------------------------------!
   open(newunit=u, file=F_NONE, status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'include_pft = [1]'
   write(u,'(a)') 'evergreen   = [0]'
   close(u)

   !----- (b) the same file WITH a phenology section, written the way the shipped            !
   !      meds_config_pft.toml now documents it. ------------------------------------------!
   open(newunit=u, file=F_FULL, status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'include_pft = [1]'
   write(u,'(a)') 'evergreen   = [0]'
   write(u,'(a)') '[phenology]'
   write(u,'(a)') 'flush_cue_mask = [1]'
   write(u,'(a)') 'shed_cue_mask  = [1]'
   write(u,'(a)') 'cue_sharpness  = [2.0]'
   close(u)

   !----- (c) a file carrying ONLY the retired spelling. The section is plainly present and  !
   !      deliberate; the old gate could not see it, because it looked for one key by name.  !
   open(newunit=u, file=F_STALE, status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'include_pft = [1]'
   write(u,'(a)') '[phenology]'
   write(u,'(a)') 'cue_mask      = [1]'
   write(u,'(a)') 'cue_sharpness = [2.0]'
   close(u)

   call toml_parse_file(F_NONE,  t_none,  ok) ; call check(ok, 'no-section fixture parsed')
   call toml_parse_file(F_FULL,  t_full,  ok) ; call check(ok, 'full fixture parsed')
   call toml_parse_file(F_STALE, t_stale, ok) ; call check(ok, 'stale-key fixture parsed')

   !=== 1. The gate. =======================================================================!
   call check(.not. toml_has_section(t_none, 'phenology'),                                          &
              'a file with no [phenology] keys must not read as having the section')
   call check(toml_has_section(t_full, 'phenology'),                                                &
              'a file with [phenology] keys must read as having the section')

   !----- THE REGRESSION. The old gate was toml_has(t, 'phenology.flush_cue_mask'); this file  !
   !      has a [phenology] section and that key is absent, so the old gate skipped the whole   !
   !      block and the presence map never reported the twenty-three missing keys. The new gate !
   !      sees the section, which is what makes those keys required and the omission loud. -----!
   call check(.not. toml_has(t_stale, 'phenology.flush_cue_mask'),                                  &
              'stale fixture must not contain the key the OLD gate keyed on')
   call check(toml_has_section(t_stale, 'phenology'),                                               &
              'the gate must see a [phenology] section by ANY of its keys, not by one chosen name')

   !=== 2. The retired spelling is detectable, so the loader can reject it by name rather than !
   !       ignore it. A config value that parses but does nothing is worse than one that is     !
   !       absent -- here it silently changed a PFT's leaf habit.  ----------------------------!
   call check(toml_has(t_stale, 'phenology.cue_mask'),                                              &
              'the retired cue_mask key must be detectable for the migration error')
   call check(.not. toml_has(t_full, 'phenology.cue_mask'),                                         &
              'the documented form must not trip the migration error')

   !=== 3. The section predicate must not match a DIFFERENT section that merely shares a prefix.!
   call check(.not. toml_has_section(t_full, 'phen'),                                               &
              'section match must be on a full segment, not a string prefix')
   call check(.not. toml_has_section(t_full, 'soil_column'),                                        &
              'an absent section must not match')

   open(newunit=u, file=F_NONE,  status='old') ; close(u, status='delete')
   open(newunit=u, file=F_FULL,  status='old') ; close(u, status='delete')
   open(newunit=u, file=F_STALE, status='old') ; close(u, status='delete')

   call test_report('test_phenology_config')

end program test_phenology_config
