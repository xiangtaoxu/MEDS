! SPDX-License-Identifier: Apache-2.0
!----- Patch fusion conserves polygon_t plant number & area; termination renormalizes area. --!
program test_patch
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t, derive_parameters
   use meds_polygon_state_vars,        only : polygon_t
   use meds_init,           only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_patch_fusefiss, only : new_fuse_patches, terminate_patches, sort_patches,  &
                                              patch_light_profile
   use meds_output_reduce, only : total_nplant, total_area
   use meds_cohort_state_vars, only : cohort_film_water
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(polygon_t)  :: poly
   real(wp)            :: n0, w_before, w_after, cas_before, cas_after, davg, film0
   real(wp), allocatable :: light1(:), light2(:)
   integer(ik)         :: ip

   call banner('patch fusion & termination conservation')
   cfg = build_test_config()

   !=== Two identical patches must fuse, conserving the polygon's N and total area. =========!
   call init_bare_ground(poly, cfg, 2_ik)   ! area 0.5 each
   call add_cohort(poly, cfg, 1_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.7_wp,  6.0_wp)
   call add_cohort(poly, cfg, 2_ik, 2_ik, 0.4_wp, 18.0_wp)
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.7_wp,  6.0_wp)
   call finalize_init(poly)
   n0 = total_nplant(poly)
   !----- Films are per plant: fusion rescales the densities and leaves the films alone, and the   !
   !      site's film water, the area- and nplant-weighted sum, is what it was. -------------------!
   poly%cohort%leaf_surf_water(1) = 0.20_wp ; poly%cohort%wood_surf_water(2) = 0.10_wp   ! patch 1
   poly%cohort%leaf_surf_water(3) = 0.40_wp                                               ! patch 2
   film0 = site_film(poly)
   call new_fuse_patches(poly, cfg)
   call check(poly%patch%n == 1_ik, 'identical patches should fuse to one')
   call check_close(total_area(poly),   1.0_wp, 1.0e-9_wp, 'total area not conserved by fusion')
   call check_close(total_nplant(poly), n0,     1.0e-9_wp, 'the polygon''s plant number not conserved by fusion')
   call check(poly%site%patch_offset(1) == 1_ik .and. poly%site%patch_count(1) == poly%patch%n,            &
              'the site map covers the fused patch block')
   call check_close(site_film(poly), film0, 1.0e-12_wp, 'patch fusion broke canopy film water conservation')

   !=== A negligible-area patch is removed and the remaining areas renormalize to 1. ======!
   !----- TWO doomed patches, not one: the merge loop reads the survivor's CSR slice, so a       !
   !      second sliver is the case a stale map would corrupt.  ----------------------------------!
   call init_bare_ground(poly, cfg, 4_ik)
   poly%patch%area(1) = 0.6_wp
   poly%patch%area(2) = 0.4_wp - 2.0e-6_wp
   poly%patch%area(3) = 1.0e-6_wp                  ! below min_patch_area
   poly%patch%area(4) = 1.0e-6_wp                  ! and another
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(poly, cfg, 3_ik, 1_ik, 0.4_wp,  8.0_wp)
   call add_cohort(poly, cfg, 4_ik, 1_ik, 0.3_wp,  7.0_wp)
   call finalize_init(poly)
   !----- Give the doomed patch content UNLIKE its neighbours', which is the realistic case (a    !
   !      sliver is usually a fresh gap). Dropping it and renormalizing the survivors would       !
   !      destroy that content outright; it is merged into the largest survivor, so the site      !
   !      total must come through untouched. ------------------------------------------------------!
   poly%patch%soil_water(1)%theta(1) = 0.20_wp
   poly%patch%soil_water(2)%theta(1) = 0.20_wp
   poly%patch%soil_water(3)%theta(1) = 0.45_wp        ! the slivers are much wetter
   poly%patch%soil_water(4)%theta(1) = 0.50_wp
   n0 = total_nplant(poly)
   w_before = sum(poly%patch%area(1:4) * [(poly%patch%soil_water(ip)%theta(1), ip = 1, 4)])
   call terminate_patches(poly, cfg)
   call check(poly%patch%n == 2_ik, 'sub-threshold patch not removed')
   call check_close(total_area(poly), 1.0_wp, 1.0e-9_wp, 'areas not renormalized to 1')
   w_after = sum(poly%patch%area(1:poly%patch%n)                                                   &
                 * [(poly%patch%soil_water(ip)%theta(1), ip = 1, poly%patch%n)])
   call check_close(w_after, w_before, 1.0e-12_wp,                                                 &
                    'terminate_patches must MERGE the sliver, not discard its content')
   call check_close(total_nplant(poly), n0, 1.0e-12_wp,                                            &
                    'terminate_patches must carry the slivers'' plants into the survivor')

   !=== Fast-biophysics reservoirs ride the patch lockstep. ================================!
   !    (a) area-weighted MERGE on fusion conserves the area-extensive store; (b) SORT keeps  !
   !    each reservoir bound to its patch. Both are silent-wrong if the lockstep misses them. !
   call init_bare_ground(poly, cfg, 2_ik)
   poly%patch%area(1) = 0.75_wp ; poly%patch%area(2) = 0.25_wp
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 12.0_wp)   ! identical cohort in each patch
   call add_cohort(poly, cfg, 2_ik, 1_ik, 0.5_wp, 12.0_wp)   ! -> identical light profiles -> fuse
   call finalize_init(poly)
   poly%patch%soil_water(1)%theta(1)       = 0.20_wp   ; poly%patch%soil_water(2)%theta(1)       = 0.40_wp
   poly%patch%cas(1)%can_enthalpy      = 1000.0_wp ; poly%patch%cas(2)%can_enthalpy      = 2000.0_wp
   poly%patch%soil_energy(1)%soil_energy(3) = 1.0e6_wp  ; poly%patch%soil_energy(2)%soil_energy(3) = 3.0e6_wp
   !----- DIFFERENT canopy depths, or this fixture cannot tell the two weightings apart. The      !
   !      canopy air is intensive per kg of AIR, and the air mass is area x depth, so the merged  !
   !      value is mass-weighted: (0.75*30*1000 + 0.25*10*2000)/(0.75*30 + 0.25*10) = 1100, not   !
   !      the 1250 an area weighting gives. Equal depths are exactly the case where the wrong     !
   !      answer and the right one agree. ---------------------------------------------------------!
   poly%patch%cas(1)%can_depth = 30.0_wp ; poly%patch%cas(2)%can_depth = 10.0_wp
   cas_before = poly%patch%area(1) * poly%patch%cas(1)%can_depth * poly%patch%cas(1)%can_enthalpy  &
              + poly%patch%area(2) * poly%patch%cas(2)%can_depth * poly%patch%cas(2)%can_enthalpy
   call new_fuse_patches(poly, cfg)
   call check(poly%patch%n == 1_ik, 'reservoir patches should fuse to one')
   call check_close(poly%patch%soil_water(1)%theta(1),       0.25_wp,   1.0e-9_wp, 'soil water not area-weighted on fusion')
   call check_close(poly%patch%soil_energy(1)%soil_energy(3), 1.5e6_wp,  1.0e-3_wp, 'soil energy not area-weighted on fusion')
   call check_close(poly%patch%cas(1)%can_depth,         25.0_wp,   1.0e-9_wp, 'CAS depth not area-weighted on fusion')
   call check_close(poly%patch%cas(1)%can_enthalpy,      1100.0_wp, 1.0e-9_wp, 'CAS enthalpy not AIR-MASS-weighted on fusion')
   !----- and the point of that weighting: the extensive content survives the merge. -------------!
   cas_after = poly%patch%area(1) * poly%patch%cas(1)%can_depth * poly%patch%cas(1)%can_enthalpy
   call check_close(cas_after, cas_before, 1.0e-12_wp, 'patch fusion must conserve canopy-air energy')

   call init_bare_ground(poly, cfg, 2_ik)
   poly%patch%age(1) = 1.0_wp ; poly%patch%age(2) = 5.0_wp             ! patch 2 is older
   poly%patch%soil_water(1)%theta(1) = 0.11_wp ; poly%patch%soil_water(2)%theta(1) = 0.55_wp
   call finalize_init(poly)
   call sort_patches(poly)                                             ! sorts age-descending: older first
   call check_close(poly%patch%age(1), 5.0_wp, 1.0e-12_wp, 'sort did not put the older patch first')
   call check_close(poly%patch%soil_water(1)%theta(1), 0.55_wp, 1.0e-12_wp, 'reservoir did not follow its patch through sort')

   !=== The tolerance ceiling: dissimilar patches stay apart even above max_patch. ===========!
   !    A closed stand and a sparse understory, and max_patch = 1 asking for one patch. Under the !
   !    shipped ceiling (0.15) they stay two; under an uncapped schedule, 0.10 x 1.5 per pass     !
   !    (0.76 by the sixth pass), they fuse.                                                      !
   cfg = build_test_config()
   cfg%max_patch = 1_ik
   call check_close(cfg%patch_light_tol * cfg%patch_light_tol_mult**(cfg%n_patch_fusion_iter - 1_ik),  &
                    cfg%patch_light_tol_max, 1.0e-12_wp, 'the last pass must use exactly the ceiling')
   call two_unlike_patches()
   allocate(light1(cfg%n_height_layers), light2(cfg%n_height_layers))
   call patch_light_profile(poly, cfg, 1_ik, light1)
   call patch_light_profile(poly, cfg, 2_ik, light2)
   davg = sum(abs(light1 - light2), mask=(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp)) &
        / real(count(light1 < 1.0_wp - 1.0e-7_wp .or. light2 < 1.0_wp - 1.0e-7_wp), wp)
   call check(davg > 0.3_wp .and. davg < 0.7_wp,                                                  &
              'fixture: the two light profiles must differ by between 0.3 and 0.7 on average')
   call new_fuse_patches(poly, cfg)
   call check(poly%patch%n == 2_ik, 'patches more different than the ceiling must stay apart above max_patch')
   !----- The uncapped schedule, spelled as a ceiling: 0.10 * 1.5**5 over six passes is a step of 1.5. -!
   cfg%patch_light_tol_max = 0.10_wp * 1.5_wp**5
   call derive_parameters(cfg)
   call check_close(cfg%patch_light_tol_mult, 1.5_wp, 1.0e-12_wp, 'ceiling 0.759375 must give the old 1.5 step')
   call two_unlike_patches()
   call new_fuse_patches(poly, cfg)
   call check(poly%patch%n == 1_ik, 'without the ceiling the same two patches fuse down to max_patch')

   write(*,'(a)') '   PASS'

contains

   !----- Two equal-area patches: a closed stand of 30 cm trees, and 5 cm saplings alone. ------!
   subroutine two_unlike_patches()
      call init_bare_ground(poly, cfg, 2_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.012_wp, 30.0_wp)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.050_wp,  5.0_wp)
      call finalize_init(poly)
   end subroutine two_unlike_patches

   !----- The site's canopy film water [kg/m2 of site]: each cohort's, area-weighted by its patch. -!
   real(wp) function site_film(p) result(w)
      type(polygon_t), intent(in) :: p
      integer(ik) :: i
      w = 0.0_wp
      do i = 1_ik, p%cohort%n
         w = w + p%patch%area(p%cohort%owner_patch(i)) * cohort_film_water(p%cohort, i)
      end do
   end function site_film
end program test_patch
