!----- CSR integrity, cohort sort order, and termination compaction. ----------------------!
program test_containers
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t, build_config
   use meds_demography_types,           only : community
   use meds_setup,           only : init_bare_ground, add_cohort, finalize_init
   use meds_cohort_dynamics, only : terminate_cohorts
   use meds_test_support,    only : check, banner
   implicit none

   type(meds_config_t) :: cfg
   type(community)     :: comm
   integer(ik)         :: ip, k, i0, i1, nbefore

   call banner('containers: CSR + sort + termination')
   cfg = build_config()
   call init_bare_ground(comm, cfg, 2_ik, 298.15_wp, 295.15_wp)

   !----- Add cohorts (same PFT so height order == DBH order), unsorted. ------------------!
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.2_wp,  5.0_wp)
   call add_cohort(comm, cfg, 2_ik, 1_ik, 0.2_wp, 20.0_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.3_wp, 15.0_wp)
   call add_cohort(comm, cfg, 2_ik, 1_ik, 0.15_wp, 3.0_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.1_wp,  9.0_wp)
   call finalize_init(comm)

   !----- CSR: counts sum to n, offsets contiguous, owners match their slice. -------------!
   call check(sum(comm%pat%ccount(1:comm%pat%n)) == comm%coh%n, 'ccount sum /= n')
   call check(comm%pat%coff(1) == 1_ik, 'coff(1) /= 1')
   call check(comm%pat%ccount(1) == 3_ik, 'patch 1 should have 3 cohorts')
   call check(comm%pat%ccount(2) == 2_ik, 'patch 2 should have 2 cohorts')
   do ip = 1_ik, comm%pat%n
      i0 = comm%pat%coff(ip)
      i1 = i0 + comm%pat%ccount(ip) - 1_ik
      if (ip < comm%pat%n) call check(comm%pat%coff(ip+1) == i1 + 1_ik, 'coff not contiguous')
      do k = i0, i1
         call check(comm%coh%owner_patch(k) == ip, 'owner_patch mismatch in slice')
      end do
      !----- Sorted height-descending within the slice. ---------------------------------!
      do k = i0, i1 - 1_ik
         call check(comm%coh%hite(k) >= comm%coh%hite(k+1), 'cohorts not height-descending')
      end do
   end do

   !----- Termination removes a sub-threshold cohort and keeps CSR consistent. ------------!
   nbefore = comm%coh%n
   comm%coh%nplant(comm%pat%coff(1)) = 1.0e-12_wp     ! below negligible_nplant
   call terminate_cohorts(comm, cfg)
   call check(comm%coh%n == nbefore - 1_ik, 'termination did not drop exactly one cohort')
   call check(sum(comm%pat%ccount(1:comm%pat%n)) == comm%coh%n, 'CSR inconsistent after terminate')

   write(*,'(a)') '   PASS'
end program test_containers
