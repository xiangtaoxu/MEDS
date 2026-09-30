! SPDX-License-Identifier: Apache-2.0
!----- Patch fusion conserves site_t plant number & area; termination renormalizes area. -----!
program test_patch
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t, derive_parameters
   use meds_site_state_types,          only : site_t
   use meds_init,           only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_patch_fusefiss, only : new_fuse_patches, terminate_patches, sort_patches,  &
                                              patch_light_profile
   use meds_diagnostic_reduce, only : total_nplant, total_area
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)     :: site
   real(wp)            :: n0, w_before, w_after, cas_before, cas_after, davg
   real(wp), allocatable :: light1(:), light2(:)
   integer(ik)         :: ip

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
   !----- REVIEW 2026-09 (item 1B #2): films are per m2 of THEIR patch; on fusion they must dilute   !
   !      by the area weights like nplant. Site total before = 0.5*(0.2+0.1) + 0.5*(0.4+0.0) = 0.35. --!
   site%cohort%leaf_surf_water(1) = 0.20_wp ; site%cohort%wood_surf_water(2) = 0.10_wp   ! patch 1
   site%cohort%leaf_surf_water(3) = 0.40_wp                                               ! patch 2
   call new_fuse_patches(site, cfg)
   call check(site%patch%n == 1_ik, 'identical patches should fuse to one')
   call check_close(total_area(site),   1.0_wp, 1.0e-9_wp, 'total area not conserved by fusion')
   call check_close(total_nplant(site), n0,     1.0e-9_wp, 'site_t plant number not conserved by fusion')
   call check_close(site%patch%area(1) * sum(site%cohort%leaf_surf_water(1:site%cohort%n)               &
                                             + site%cohort%wood_surf_water(1:site%cohort%n)),            &
                    0.35_wp, 1.0e-12_wp, 'patch fusion broke canopy film water conservation')

   !=== A negligible-area patch is removed and the remaining areas renormalize to 1. ======!
   !----- TWO doomed patches, not one: the merge loop reads the survivor's CSR slice, so a       !
   !      second sliver is the case a stale map would corrupt.  ----------------------------------!
   call init_bare_ground(site, cfg, 4_ik)
   site%patch%area(1) = 0.6_wp
   site%patch%area(2) = 0.4_wp - 2.0e-6_wp
   site%patch%area(3) = 1.0e-6_wp                  ! below min_patch_area
   site%patch%area(4) = 1.0e-6_wp                  ! and another
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 3_ik, 1_ik, 0.4_wp,  8.0_wp)
   call add_cohort(site, cfg, 4_ik, 1_ik, 0.3_wp,  7.0_wp)
   call finalize_init(site)
   !----- Give the doomed patch content UNLIKE its neighbours', which is the realistic case (a    !
   !      sliver is usually a fresh gap). Dropping it and renormalizing the survivors used to     !
   !      destroy that content outright; it is merged into the largest survivor now, so the site  !
   !      total must come through untouched (review item 1B #9).  --------------------------------!
   site%patch%soil_w(1)%theta(1) = 0.20_wp
   site%patch%soil_w(2)%theta(1) = 0.20_wp
   site%patch%soil_w(3)%theta(1) = 0.45_wp        ! the slivers are much wetter
   site%patch%soil_w(4)%theta(1) = 0.50_wp
   n0 = total_nplant(site)
   w_before = sum(site%patch%area(1:4) * [(site%patch%soil_w(ip)%theta(1), ip = 1, 4)])
   call terminate_patches(site, cfg)
   call check(site%patch%n == 2_ik, 'sub-threshold patch not removed')
   call check_close(total_area(site), 1.0_wp, 1.0e-9_wp, 'areas not renormalized to 1')
   w_after = sum(site%patch%area(1:site%patch%n)                                                   &
                 * [(site%patch%soil_w(ip)%theta(1), ip = 1, site%patch%n)])
   call check_close(w_after, w_before, 1.0e-12_wp,                                                 &
                    'terminate_patches must MERGE the sliver, not discard its content')
   call check_close(total_nplant(site), n0, 1.0e-12_wp,                                            &
                    'terminate_patches must carry the slivers'' plants into the survivor')

   !=== Fast-biophysics reservoirs ride the patch lockstep. ================================!
   !    (a) area-weighted MERGE on fusion conserves the area-extensive store; (b) SORT keeps  !
   !    each reservoir bound to its patch. Both are silent-wrong if the lockstep misses them. !
   call init_bare_ground(site, cfg, 2_ik)
   site%patch%area(1) = 0.75_wp ; site%patch%area(2) = 0.25_wp
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 12.0_wp)   ! identical cohort in each patch
   call add_cohort(site, cfg, 2_ik, 1_ik, 0.5_wp, 12.0_wp)   ! -> identical light profiles -> fuse
   call finalize_init(site)
   site%patch%soil_w(1)%theta(1)       = 0.20_wp   ; site%patch%soil_w(2)%theta(1)       = 0.40_wp
   site%patch%cas(1)%can_enthalpy      = 1000.0_wp ; site%patch%cas(2)%can_enthalpy      = 2000.0_wp
   site%patch%soil_e(1)%soil_energy(3) = 1.0e6_wp  ; site%patch%soil_e(2)%soil_energy(3) = 3.0e6_wp
   !----- DIFFERENT canopy depths, or this fixture cannot tell the two weightings apart. The      !
   !      canopy air is intensive per kg of AIR, and the air mass is area x depth, so the merged  !
   !      value is mass-weighted: (0.75*30*1000 + 0.25*10*2000)/(0.75*30 + 0.25*10) = 1100, not   !
   !      the 1250 an area weighting gives. Both patches sat at the default 20 m until now, which !
   !      is exactly the case where the wrong answer and the right one agree.  -------------------!
   site%patch%cas(1)%can_depth = 30.0_wp ; site%patch%cas(2)%can_depth = 10.0_wp
   cas_before = site%patch%area(1) * site%patch%cas(1)%can_depth * site%patch%cas(1)%can_enthalpy  &
              + site%patch%area(2) * site%patch%cas(2)%can_depth * site%patch%cas(2)%can_enthalpy
   call new_fuse_patches(site, cfg)
   call check(site%patch%n == 1_ik, 'reservoir patches should fuse to one')
   call check_close(site%patch%soil_w(1)%theta(1),       0.25_wp,   1.0e-9_wp, 'soil water not area-weighted on fusion')
   call check_close(site%patch%soil_e(1)%soil_energy(3), 1.5e6_wp,  1.0e-3_wp, 'soil energy not area-weighted on fusion')
   call check_close(site%patch%cas(1)%can_depth,         25.0_wp,   1.0e-9_wp, 'CAS depth not area-weighted on fusion')
   call check_close(site%patch%cas(1)%can_enthalpy,      1100.0_wp, 1.0e-9_wp, 'CAS enthalpy not AIR-MASS-weighted on fusion')
   !----- and the point of that weighting: the extensive content survives the merge. -------------!
   cas_after = site%patch%area(1) * site%patch%cas(1)%can_depth * site%patch%cas(1)%can_enthalpy
   call check_close(cas_after, cas_before, 1.0e-12_wp, 'patch fusion must conserve canopy-air energy')

   call init_bare_ground(site, cfg, 2_ik)
   site%patch%age(1) = 1.0_wp ; site%patch%age(2) = 5.0_wp             ! patch 2 is older
   site%patch%soil_w(1)%theta(1) = 0.11_wp ; site%patch%soil_w(2)%theta(1) = 0.55_wp
   call finalize_init(site)
   call sort_patches(site)                                             ! sorts age-descending: older first
   call check_close(site%patch%age(1), 5.0_wp, 1.0e-12_wp, 'sort did not put the older patch first')
   call check_close(site%patch%soil_w(1)%theta(1), 0.55_wp, 1.0e-12_wp, 'reservoir did not follow its patch through sort')

   !=== The tolerance ceiling: dissimilar patches stay apart even above max_patch. ===========!
   !    A closed stand and a sparse understory, and max_patch = 1 asking for one patch. Under the !
   !    shipped ceiling (0.15) they stay two; under the old schedule, 0.10 x 1.5 per pass with no  !
   !    ceiling (0.76 by the sixth pass), they fuse.                                               !
   cfg = build_test_config()
   cfg%max_patch = 1_ik
   call check_close(cfg%patch_light_tol * cfg%patch_light_tol_mult**(cfg%n_patch_fusion_iter - 1_ik),  &
                    cfg%patch_light_tol_max, 1.0e-12_wp, 'the last pass must use exactly the ceiling')
   call two_unlike_patches()
   allocate(light1(cfg%n_height_layers), light2(cfg%n_height_layers))
   call patch_light_profile(site, cfg, 1_ik, light1)
   call patch_light_profile(site, cfg, 2_ik, light2)
   davg = sum(abs(light1 - light2), mask=(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp)) &
        / real(count(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp), wp)
   call check(davg > 0.3_wp .and. davg < 0.7_wp,                                                  &
              'fixture: the two light profiles must differ by between 0.3 and 0.7 on average')
   call new_fuse_patches(site, cfg)
   call check(site%patch%n == 2_ik, 'patches more different than the ceiling must stay apart above max_patch')
   !----- The old schedule, spelled as a ceiling: 0.10 * 1.5**5 over six passes is a step of 1.5. --!
   cfg%patch_light_tol_max = 0.10_wp * 1.5_wp**5
   call derive_parameters(cfg)
   call check_close(cfg%patch_light_tol_mult, 1.5_wp, 1.0e-12_wp, 'ceiling 0.759375 must give the old 1.5 step')
   call two_unlike_patches()
   call new_fuse_patches(site, cfg)
   call check(site%patch%n == 1_ik, 'without the ceiling the same two patches fuse down to max_patch')

   write(*,'(a)') '   PASS'

contains

   !----- Two equal-area patches: a closed stand of 30 cm trees, and 5 cm saplings alone. ------!
   subroutine two_unlike_patches()
      call init_bare_ground(site, cfg, 2_ik)
      call add_cohort(site, cfg, 1_ik, 1_ik, 0.012_wp, 30.0_wp)
      call add_cohort(site, cfg, 2_ik, 1_ik, 0.050_wp,  5.0_wp)
      call finalize_init(site)
   end subroutine two_unlike_patches
end program test_patch
