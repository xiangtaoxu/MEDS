! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_phenology_config -- the [phenology] section actually reaches the PFT table.         !
!                                                                                          !
! The gate on this section is its PRESENCE in the PFT file: absent, every PFT keeps the    !
! built-in no-cue defaults; present, every key in it is REQUIRED. A gate that tests for    !
! ONE key inside the block skips, in SILENCE, a block written with any other spelling: the !
! presence map never sees the required keys go missing, and every PFT falls back to        !
! CUE_NONE / CUE_NONE -- always flushing, never senescing -- whatever leaf habit it        !
! declares. A cold-deciduous stand then keeps its summer leaf area through every January,  !
! and no budget notices.                                                                   !
!                                                                                          !
! What is asserted here is the GATE, on both sides of it:                                  !
!                                                                                          !
!   1. a file with no [phenology] keys reads as having no section, a file with them as     !
!      having it, and a file carrying only a retired spelling still reads as having it:    !
!      `toml_has_section` sees a section by ANY of its keys, not by one chosen name;       !
!   2. the retired spelling `cue_mask` is detectable, so the loader can reject it by name; !
!   3. the section predicate matches a full segment, not a string prefix.                  !
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

   call banner('[phenology] section gate')

   !----- (a) a PFT file with no phenology section at all. ---------------------------------!
   open(newunit=u, file=F_NONE, status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'include_pft = [1]'
   close(u)

   !----- (b) the same file WITH a phenology section, written the way the shipped            !
   !      meds_config_pft.toml documents it. ----------------------------------------------!
   open(newunit=u, file=F_FULL, status='replace', action='write')
   write(u,'(a)') '[pft]'
   write(u,'(a)') 'include_pft = [1]'
   write(u,'(a)') '[phenology]'
   write(u,'(a)') 'flush_cue_mask = [1]'
   write(u,'(a)') 'shed_cue_mask  = [1]'
   write(u,'(a)') 'flush_rate_max = [0.05]'
   close(u)

   !----- (c) a file carrying ONLY the retired spelling. The section is plainly present and  !
   !      deliberate; a gate that looks for one key by name cannot see it.                   !
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

   !----- The defect this guards: a gate of toml_has(t, 'phenology.flush_cue_mask') skips this !
   !      file's [phenology] section, because that key is absent, and the presence map never   !
   !      reports the missing keys. toml_has_section sees the section, which is what makes     !
   !      those keys required and the omission loud. ------------------------------------------!
   call check(.not. toml_has(t_stale, 'phenology.flush_cue_mask'),                                  &
              'stale fixture must not contain the key the OLD gate keyed on')
   call check(toml_has_section(t_stale, 'phenology'),                                               &
              'the gate must see a [phenology] section by ANY of its keys, not by one chosen name')

   !=== 2. The retired spelling is detectable, so the loader can reject it by name rather than !
   !       ignore it. A config value that parses but does nothing is worse than one that is     !
   !       absent -- here it would silently change a PFT's leaf habit. ------------------------!
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
