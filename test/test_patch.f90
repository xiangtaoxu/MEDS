!----- Patch fusion conserves site_t plant number & area; termination renormalizes area. -----!
program test_patch
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_demography_types,          only : site_t
   use meds_init,           only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_structure, only : new_fuse_patches, terminate_patches
   use meds_demography_diagnostics, only : total_nplant, total_area
   use meds_test_support, only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   real(wp)            :: n0

   call banner('patch fusion & termination conservation')
   cfg = build_test_config()

   !=== Two identical patches must fuse, conserving site_t N and total area. ================!
   call init_bare_ground(site, cfg, 2_ik)   ! area 0.5 each
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.7_wp,  6.0_wp)
   call add_cohort(site, cfg, 2_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(site, cfg, 2_ik, 1_ik, 0.7_wp,  6.0_wp)
   call finalize_init(site)
   n0 = total_nplant(site)
   call new_fuse_patches(site, cfg)
   call check(site%patch%n == 1_ik, 'identical patches should fuse to one')
   call check_close(total_area(site),   1.0_wp, 1.0e-9_wp, 'total area not conserved by fusion')
   call check_close(total_nplant(site), n0,     1.0e-9_wp, 'site_t plant number not conserved by fusion')

   !=== A negligible-area patch is removed and the remaining areas renormalize to 1. ======!
   call init_bare_ground(site, cfg, 3_ik)
   site%patch%area(1) = 0.6_wp
   site%patch%area(2) = 0.4_wp - 1.0e-6_wp
   site%patch%area(3) = 1.0e-6_wp                  ! below min_patch_area
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call finalize_init(site)
   call terminate_patches(site, cfg)
   call check(site%patch%n == 2_ik, 'sub-threshold patch not removed')
   call check_close(total_area(site), 1.0_wp, 1.0e-9_wp, 'areas not renormalized to 1')

   write(*,'(a)') '   PASS'
end program test_patch
