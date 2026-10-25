! SPDX-License-Identifier: Apache-2.0
!----- CSR integrity, cohort sort order, and termination compaction. ----------------------!
program test_containers
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_polygon_state_vars,         only : polygon_t
   use meds_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use meds_slow_update, only : terminate_cohorts
   use meds_diag_vars,       only : alloc_cohort_diag, cohort_diag_value
   use meds_test_support, only : banner, build_test_config, check
   use meds_slow_state,  only : kill_report_t
   implicit none
   type(kill_report_t)  :: kills   !< what the cull and the treefall killed
   real(wp)            :: shed_water   !< an operator's shed water, unread here

   type(meds_config_t) :: cfg
   type(polygon_t)  :: poly
   integer(ik)         :: ip, k, i0, i1, nbefore, a, b, gid_dropped
   integer(ik)         :: nread
   real(wp)            :: xbuf(8), xshort(2)

   call banner('containers: CSR + sort + termination')
   cfg = build_test_config()
   call init_bare_ground(poly, cfg, 2_ik)

   !----- Add cohorts (same PFT so height order == DBH order), unsorted. ------------------!
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.2_wp,  5.0_wp)
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.2_wp, 20.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 15.0_wp)
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.15_wp, 3.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.1_wp,  9.0_wp)
   call finalize_init(poly)

   !----- CSR: counts sum to n, offsets contiguous, owners match their slice. -------------!
   call check(sum(poly%patch%cohort_count(1:poly%patch%n)) == poly%cohort%n, 'cohort_count sum /= n')
   call check(poly%patch%cohort_offset(1) == 1_ik, 'cohort_offset(1) /= 1')
   call check(poly%patch%cohort_count(1) == 3_ik, 'patch 1 should have 3 cohorts')
   call check(poly%patch%cohort_count(2) == 2_ik, 'patch 2 should have 2 cohorts')
   do ip = 1_ik, poly%patch%n
      i0 = poly%patch%cohort_offset(ip)
      i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
      if (ip < poly%patch%n) call check(poly%patch%cohort_offset(ip+1) == i1 + 1_ik, 'cohort_offset not contiguous')
      do k = i0, i1
         call check(poly%cohort%owner_patch(k) == ip, 'owner_patch mismatch in slice')
      end do
      !----- Sorted height-descending within the slice. ---------------------------------!
      do k = i0, i1 - 1_ik
         call check(poly%cohort%height(k) >= poly%cohort%height(k+1), 'cohorts not height-descending')
      end do
   end do

   !----- Persistent global ids: assigned (>0) and unique across the live set, even after   !
   !       the sort permuted the cohort order. --------------------------------------------!
   call check(all(poly%cohort%global_id(1:poly%cohort%n) > 0_ik), 'global_id not assigned')
   do a = 1_ik, poly%cohort%n
      do b = a + 1_ik, poly%cohort%n
         call check(poly%cohort%global_id(a) /= poly%cohort%global_id(b), 'global_id not unique')
      end do
   end do

   !----- Termination removes a sub-threshold cohort and keeps CSR consistent. ------------!
   !                                                                                          !
   !      THE DIAGNOSTIC BLOCKS RIDE THE SAME COMPACTION, COUNT INCLUDED. They are armed   !
   !      here because a cull is the operator that SHRINKS the cohort array, and a block   !
   !      whose `n` stays at the pre-cull value claims more live slots than the array has. !
   !      extract_variable sizes its scratch buffer from cohort%n and cohort_diag_value    !
   !      fills 1..d%n into it, so a stale count is an out-of-bounds WRITE -- silent in    !
   !      Release, a SIGSEGV some steps later. A stand that never shrinks cannot show it   !
   !      (the Ithaca reference goes 114 -> 125 cohorts over three years). ----------------!
   call alloc_cohort_diag(poly%cohort%diag,  poly%cohort%cap, .true.)
   call alloc_cohort_diag(poly%cohort%sdiag, poly%cohort%cap, .true., slow_rows=.true.)
   poly%cohort%diag%n  = poly%cohort%n
   poly%cohort%sdiag%n = poly%cohort%n

   nbefore = poly%cohort%n
   gid_dropped = poly%cohort%global_id(poly%patch%cohort_offset(1))
   poly%cohort%nplant(poly%patch%cohort_offset(1)) = 1.0e-12_wp     ! below negligible_nplant
   call terminate_cohorts(poly, cfg, shed_water, kills)
   call check(poly%cohort%n == nbefore - 1_ik, 'termination did not drop exactly one cohort')
   call check(sum(poly%patch%cohort_count(1:poly%patch%n)) == poly%cohort%n, 'CSR inconsistent after terminate')
   !----- The culled cohort's id never reappears (ids are not reused). --------------------!
   call check(all(poly%cohort%global_id(1:poly%cohort%n) /= gid_dropped), 'dropped global_id reused')

   call check(poly%cohort%diag%n  == poly%cohort%n, 'fast diag block kept a stale n across the cull')
   call check(poly%cohort%sdiag%n == poly%cohort%n, 'slow diag block kept a stale n across the cull')

   !----- ...and independently, the extractor must never write past the caller's array, however  !
   !      the counts got out of step. `xshort` is deliberately smaller than the live cohort       !
   !      count: the reader has to report what it actually wrote, not what the block holds. ------!
   call cohort_diag_value(poly%cohort%diag, 1_ik, xbuf, nread, 0.0_wp)
   call check(nread == poly%cohort%n, 'extractor did not report the live cohort count')
   call cohort_diag_value(poly%cohort%diag, 1_ik, xshort, nread, 0.0_wp)
   call check(nread <= int(size(xshort), ik), 'extractor reported more than the caller''s array holds')

   write(*,'(a)') '   PASS'
end program test_containers
