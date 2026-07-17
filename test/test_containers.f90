!----- CSR integrity, cohort sort order, and termination compaction. ----------------------!
program test_containers
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_core_state_types,           only : site_t
   use meds_init,            only : init_bare_ground, add_cohort, finalize_init
   use meds_core_cohort_fusefiss, only : terminate_cohorts
   use meds_test_support, only : build_test_config, check, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   integer(ik)         :: ip, k, i0, i1, nbefore, a, b, gid_dropped

   call banner('containers: CSR + sort + termination')
   cfg = build_test_config()
   call init_bare_ground(site, cfg, 2_ik)

   !----- Add cohorts (same PFT so height order == DBH order), unsorted. ------------------!
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.2_wp,  5.0_wp)
   call add_cohort(site, cfg, 2_ik, 1_ik, 0.2_wp, 20.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 15.0_wp)
   call add_cohort(site, cfg, 2_ik, 1_ik, 0.15_wp, 3.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.1_wp,  9.0_wp)
   call finalize_init(site)

   !----- CSR: counts sum to n, offsets contiguous, owners match their slice. -------------!
   call check(sum(site%patch%cohort_count(1:site%patch%n)) == site%cohort%n, 'cohort_count sum /= n')
   call check(site%patch%cohort_offset(1) == 1_ik, 'cohort_offset(1) /= 1')
   call check(site%patch%cohort_count(1) == 3_ik, 'patch 1 should have 3 cohorts')
   call check(site%patch%cohort_count(2) == 2_ik, 'patch 2 should have 2 cohorts')
   do ip = 1_ik, site%patch%n
      i0 = site%patch%cohort_offset(ip)
      i1 = i0 + site%patch%cohort_count(ip) - 1_ik
      if (ip < site%patch%n) call check(site%patch%cohort_offset(ip+1) == i1 + 1_ik, 'cohort_offset not contiguous')
      do k = i0, i1
         call check(site%cohort%owner_patch(k) == ip, 'owner_patch mismatch in slice')
      end do
      !----- Sorted height-descending within the slice. ---------------------------------!
      do k = i0, i1 - 1_ik
         call check(site%cohort%height(k) >= site%cohort%height(k+1), 'cohorts not height-descending')
      end do
   end do

   !----- Persistent global ids: assigned (>0) and unique across the live set, even after   !
   !       the sort permuted the cohort order. --------------------------------------------!
   call check(all(site%cohort%global_id(1:site%cohort%n) > 0_ik), 'global_id not assigned')
   do a = 1_ik, site%cohort%n
      do b = a + 1_ik, site%cohort%n
         call check(site%cohort%global_id(a) /= site%cohort%global_id(b), 'global_id not unique')
      end do
   end do

   !----- Termination removes a sub-threshold cohort and keeps CSR consistent. ------------!
   nbefore = site%cohort%n
   gid_dropped = site%cohort%global_id(site%patch%cohort_offset(1))
   site%cohort%nplant(site%patch%cohort_offset(1)) = 1.0e-12_wp     ! below negligible_nplant
   call terminate_cohorts(site, cfg)
   call check(site%cohort%n == nbefore - 1_ik, 'termination did not drop exactly one cohort')
   call check(sum(site%patch%cohort_count(1:site%patch%n)) == site%cohort%n, 'CSR inconsistent after terminate')
   !----- The culled cohort's id never reappears (ids are not reused). --------------------!
   call check(all(site%cohort%global_id(1:site%cohort%n) /= gid_dropped), 'dropped global_id reused')

   write(*,'(a)') '   PASS'
end program test_containers
