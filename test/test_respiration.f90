! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_respiration -- the plants' maintenance respiration (meds_respiration;               !
! docs/science/plant_respiration.md eqs. (3)-(5)).                                         !
!   1. 25 C (root_zone_temp_scale, fine_root_maintenance_respiration,                      !
!      stem_maintenance_respiration): the temperature response is 1 at 25 C, so the roots  !
!      respire their rate times their carbon and the sapwood its rate times its volume.    !
!   2. LAYERS (root_zone_temp_scale, eq. (4)): the response summed over the root layers is !
!      the mean temperature's on a uniform profile, and higher on a profile spread about  !
!      the same mean, the response being convex below its optimum.                         !
!   3. SAPWOOD (stem_maintenance_respiration, eq. (3)): the sapwood term adds to the       !
!      surface term.                                                                       !
!   4. GRASS (stem_maintenance_respiration): a cohort that is not woody respires no stem.  !
!   5. ABOVEGROUND FRACTION (stem_maintenance_respiration): the stem area is the plant's   !
!      above ground over its aboveground fraction, so halving the fraction doubles it.     !
! Growth respiration is the daily allocator's (test_plant_carbon_allocation).              !
!==========================================================================================!
program test_respiration
   use meds_test_assert, only : check_close, check_true, test_report
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : t_ref_photo
   use meds_respiration, only : stem_maintenance_respiration, fine_root_maintenance_respiration, root_zone_temp_scale
   implicit none
   !----- Leaf Rd's peaked response, as the reference configuration sets it ([leaf_physiology]). -!
   real(wp), parameter :: EA_RD = 46390.0_wp, HD_RD = 200000.0_wp, DS_RD = 490.0_wp
   real(wp), parameter :: AGF  = 0.7_wp    !< [-] the PFT aboveground fraction, passed per cohort
   real(wp), parameter :: SRF  = 0.06_wp   !< [umol/m2 stem/s @25C] the PFT stem baseline
   real(wp), parameter :: RRF  = 0.30_wp   !< [umol/kgC root/s @25C] the PFT root baseline
   real(wp), parameter :: RSAP = 40.0_wp   !< [umol/m3 sapwood/s @25C] the PFT sapwood rate
   real(wp), parameter :: VSAP = 0.05_wp   !< [m3/plant] a cohort's sapwood volume
   logical,  parameter :: WOODY = .true.

   call test_reference_temperature()
   call test_root_layers()
   call test_sapwood_adds()
   call test_grass()
   call test_aboveground_frac()

   call test_report('test_respiration')

contains

   !----- 1. At 25 C the peaked response is 1: the roots respire their 25 C rate times their      !
   !      carbon, and the sapwood, with no surface rate, its rate times its volume. --------------!
   subroutine test_reference_temperature()
      real(wp) :: r_root, r_sap
      call fine_root_maintenance_respiration(root_zone_temp_scale([t_ref_photo], [1.0_wp], 1_ik, EA_RD, HD_RD, DS_RD), &
                                             2.0_wp, RRF, r_root)
      call check_close('root @25C = factor25 * broot', r_root, RRF * 2.0_wp)
      call stem_maintenance_respiration(t_ref_photo, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,   &
                                        0.0_wp, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_sap)
      call check_close('sapwood @25C = rate x volume', r_sap, RSAP * VSAP)
   end subroutine test_reference_temperature

   !----- 2. The response is summed OVER LAYERS (eq. (4)), not taken at a mean temperature. A      !
   !      UNIFORM profile must give exactly the mean-temperature answer (the layer sum is inert     !
   !      where there is no gradient), and a SPLIT profile with the same root-weighted mean (285 K, !
   !      over a 40 K range) must NOT: below the optimum the response is convex, so the layered sum !
   !      is higher, which a mean-temperature formulation cannot express. --------------------------!
   subroutine test_root_layers()
      real(wp), parameter :: FRAC(4) = [0.4_wp, 0.3_wp, 0.2_wp, 0.1_wp]
      real(wp), parameter :: T_UNIFORM(4) = [285.0_wp, 285.0_wp, 285.0_wp, 285.0_wp]
      real(wp), parameter :: T_SPLIT(4) = [275.0_wp, 295.0_wp, 275.0_wp, 315.0_wp]
      real(wp) :: f_mean, f_layer
      f_mean  = root_zone_temp_scale([285.0_wp], [1.0_wp], 1_ik, EA_RD, HD_RD, DS_RD)
      f_layer = root_zone_temp_scale(T_UNIFORM, FRAC, 4_ik, EA_RD, HD_RD, DS_RD)
      call check_close('uniform profile: layered == mean-temperature', f_layer, f_mean)
      f_layer = root_zone_temp_scale(T_SPLIT, FRAC, 4_ik, EA_RD, HD_RD, DS_RD)
      call check_true('a temperature GRADIENT changes the response at fixed mean',              &
                      abs(f_layer - f_mean) > 1.0e-3_wp * f_mean, f_layer / f_mean)
      call check_true('convex below the optimum => layered exceeds mean-temperature',           &
                      f_layer > f_mean, f_layer / f_mean)
   end subroutine test_root_layers

   !----- 3. The surface and sapwood terms share the temperature response and add. ----------------!
   subroutine test_sapwood_adds()
      real(wp) :: r_sap, r_surf, r_both
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        0.0_wp, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_sap)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        SRF, 0.0_wp, 0.0_wp, EA_RD, HD_RD, DS_RD, r_surf)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        SRF, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_both)
      call check_close('surface and sapwood terms add', r_both, r_surf + r_sap)
   end subroutine test_sapwood_adds

   !----- 4. is_woody is a per-PFT trait: a grass cohort respires no stem, surface or sapwood. ----!
   subroutine test_grass()
      real(wp) :: r_grass
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, .false.,    &
                                        SRF, VSAP, RSAP, EA_RD, HD_RD, DS_RD, r_grass)
      call check_true('a grass cohort respires no stem', r_grass == 0.0_wp, r_grass)
   end subroutine test_grass

   !----- 5. The aboveground fraction is a per-cohort trait and must scale the stem area. A        !
   !      run-uniform value would make two PFTs with different allocation respire identically      !
   !      while demography used their real values. The kernel divides the stem area by it, so      !
   !      halving the fraction must double the respiration exactly. -------------------------------!
   subroutine test_aboveground_frac()
      real(wp) :: r_high, r_low
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.7_wp,  WOODY, SRF, 0.0_wp,      &
                                        0.0_wp, EA_RD, HD_RD, DS_RD, r_high)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.35_wp, WOODY, SRF, 0.0_wp,      &
                                        0.0_wp, EA_RD, HD_RD, DS_RD, r_low)
      call check_close('aboveground_frac scales stem respiration inversely', r_low, 2.0_wp * r_high)
   end subroutine test_aboveground_frac

end program test_respiration
