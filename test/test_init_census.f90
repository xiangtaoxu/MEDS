! SPDX-License-Identifier: Apache-2.0
!----- Census restart: build a site from a cohort CSV (site filter, patch map, allometry). --!
program test_init_census
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_init,             only : init_from_census
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik)         :: u, i
   logical             :: found

   call banner('init: from census CSV')
   cfg = build_test_config()

   !----- Write a small census file (with a comment + blank line to exercise skipping). ----!
   open(newunit=u, file='test_census_tmp.csv', status='replace', action='write')
   write(u,'(a)') '# pseudo cohort census (two patches in site 1, one row for site 2)'
   write(u,'(a)') 'site_id,patch_id,cohort_id,dbh,height,pft,nplant'
   write(u,'(a)') ''
   write(u,'(a)') '1,1,1,30.0,21.0,3,0.05'
   write(u,'(a)') '1,1,2,10.0,12.0,2,0.20'
   write(u,'(a)') '1,2,3, 5.0, 8.0,1,0.40'
   write(u,'(a)') '2,1,1,40.0,27.0,3,0.02'
   close(u)

   !=== Default: first site_id (1) selected; site 2 row excluded. ==========================!
   call init_from_census(site, cfg, 'test_census_tmp.csv', found)
   call check(found, 'census file should be usable')
   call check(site%patch%n  == 2_ik, 'site 1 should yield 2 patches')
   call check(site%cohort%n == 3_ik, 'site 1 should yield 3 cohorts (site 2 row excluded)')
   call check(sum(site%patch%cohort_count(1:site%patch%n)) == site%cohort%n, 'CSR count mismatch')
   !----- dbh drove the allometry: every cohort has a positive derived height/agb/leaf area. !
   call check(all(site%cohort%height(1:site%cohort%n)    > 0.0_wp), 'derived height not positive')
   call check(all(site%cohort%agb(1:site%cohort%n)       > 0.0_wp), 'derived agb not positive')
   call check(all(site%cohort%leaf_area(1:site%cohort%n) > 0.0_wp), 'derived leaf area not positive')

   !=== Explicit site selection: site 2 -> one patch, one (climax, pft 3) cohort. ==========!
   call init_from_census(site, cfg, 'test_census_tmp.csv', found, use_site_id = 2_ik)
   call check(found, 'site 2 should be usable')
   call check(site%patch%n  == 1_ik, 'site 2 should yield 1 patch')
   call check(site%cohort%n == 1_ik, 'site 2 should yield 1 cohort')
   call check(site%cohort%pft(1) == 3_ik, 'site 2 cohort should be pft 3')

   !=== Header-LESS CSV: the first data line must be KEPT, not silently dropped as a header. ==!
   open(newunit=u, file='test_census_nohdr.csv', status='replace', action='write')
   write(u,'(a)') '1,1,1,30.0,21.0,3,0.05'      ! no header row -> the first line is DATA
   write(u,'(a)') '1,1,2,10.0,12.0,2,0.20'
   write(u,'(a)') '1,2,3, 5.0, 8.0,1,0.40'
   close(u)
   call init_from_census(site, cfg, 'test_census_nohdr.csv', found)
   call check(found, 'header-less census should be usable')
   call check(site%cohort%n == 3_ik, 'header-less census must keep all 3 cohorts (first not dropped)')
   open(newunit=u, file='test_census_nohdr.csv', status='old', action='read')
   close(u, status='delete')

   !=== Columns by NAME, in any order, with patch_area and without site_id. =================!
   !    Patch 7 has three times the area of patch 9, so the site fractions are 0.75 and 0.25.  !
   open(newunit=u, file='test_census_named.csv', status='replace', action='write')
   write(u,'(a)') '# columns reordered, areas in m2'
   write(u,'(a)') 'nplant,pft,dbh,patch_area,patch_id'
   write(u,'(a)') '0.05,1,30.0,300.0,7'
   write(u,'(a)') '0.20,1,10.0,300.0,7'
   write(u,'(a)') '0.40,1,5.0,100.0,9'
   close(u)
   call init_from_census(site, cfg, 'test_census_named.csv', found)
   call check(found, 'named-column census should be usable')
   call check(site%patch%n  == 2_ik, 'named-column census should yield 2 patches')
   call check(site%cohort%n == 3_ik, 'named-column census should yield 3 cohorts')
   call check_close(sum(site%patch%area(1:2)), 1.0_wp, 1.0e-12_wp, 'patch areas must sum to the site')
   call check_close(site%patch%area(1), 0.75_wp, 1.0e-12_wp, 'patch_area not normalized to the site (7)')
   call check_close(site%patch%area(2), 0.25_wp, 1.0e-12_wp, 'patch_area not normalized to the site (9)')
   do i = 1_ik, site%cohort%n
      if (abs(site%cohort%dbh(i) - 5.0_wp) < 1.0e-9_wp) then
         call check(site%cohort%owner_patch(i) == 2_ik, 'the 5 cm cohort must sit in patch 9')
         call check_close(site%cohort%nplant(i), 0.40_wp, 1.0e-12_wp, 'nplant not read by name')
      end if
   end do
   open(newunit=u, file='test_census_named.csv', status='old', action='read')
   close(u, status='delete')
   !=== A missing file is reported (found=.false.), not a crash. ===========================!
   call init_from_census(site, cfg, 'no_such_census_file.csv', found)
   call check(.not. found, 'missing census file should report found=.false.')

   open(newunit=u, file='test_census_tmp.csv', status='old', action='read')
   close(u, status='delete')

   write(*,'(a)') '   PASS'
end program test_init_census
