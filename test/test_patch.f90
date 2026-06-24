!----- Patch fusion conserves site plant number & area; termination renormalizes area. -----!
program test_patch
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t, build_config
   use meds_demography_types,          only : community
   use meds_setup,          only : init_bare_ground, add_cohort, finalize_init
   use meds_patch_dynamics, only : new_fuse_patches, terminate_patches
   use meds_diagnostics,    only : total_nplant, total_area
   use meds_test_support,   only : check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(community)     :: comm
   real(wp)            :: n0

   call banner('patch fusion & termination conservation')
   cfg = build_config()

   !=== Two identical patches must fuse, conserving site N and total area. ================!
   call init_bare_ground(comm, cfg, 2_ik, 298.15_wp, 295.15_wp)   ! area 0.5 each
   call add_cohort(comm, cfg, 1_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.7_wp,  6.0_wp)
   call add_cohort(comm, cfg, 2_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(comm, cfg, 2_ik, 1_ik, 0.7_wp,  6.0_wp)
   call finalize_init(comm)
   n0 = total_nplant(comm)
   call new_fuse_patches(comm, cfg)
   call check(comm%pat%n == 1_ik, 'identical patches should fuse to one')
   call check_close(total_area(comm),   1.0_wp, 1.0e-9_wp, 'total area not conserved by fusion')
   call check_close(total_nplant(comm), n0,     1.0e-9_wp, 'site plant number not conserved by fusion')

   !=== A negligible-area patch is removed and the remaining areas renormalize to 1. ======!
   call init_bare_ground(comm, cfg, 3_ik, 298.15_wp, 295.15_wp)
   comm%pat%area(1) = 0.6_wp
   comm%pat%area(2) = 0.4_wp - 1.0e-6_wp
   comm%pat%area(3) = 1.0e-6_wp                  ! below min_patch_area
   call add_cohort(comm, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call finalize_init(comm)
   call terminate_patches(comm, cfg)
   call check(comm%pat%n == 2_ik, 'sub-threshold patch not removed')
   call check_close(total_area(comm), 1.0_wp, 1.0e-9_wp, 'areas not renormalized to 1')

   write(*,'(a)') '   PASS'
end program test_patch
