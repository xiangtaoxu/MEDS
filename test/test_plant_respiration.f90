!==========================================================================================!
! test_plant_respiration -- unit tests for the non-leaf respiration kernels.                 !
!                                                                                          !
!   1. GRASS      : is_woody = .false. => stem respiration is identically zero.               !
!   2. STEM @25C  : at 25 degC the peaked scale is 1, so R = factor25 * stem_area (cylinder).   !
!   3. T-RESP     : at T /= 25 degC, R = factor25 * peaked_scale(T) * stem_area.                    !
!   4. WAI        : the WAI branch term adds pi*wai/nplant/aboveground_frac to the stem area.       !
!   5. SIZE       : the DBH size scaler multiplies the baseline by 10^(scaler*dbh).                  !
!   6. ROOT       : R = factor25 * broot at 25 degC; broot = 0 => 0; monotone up to the optimum.     !
!   (Growth/construction respiration moved to meds_plant_carbon_allocation; see test_plant_carbon_allocation.)!
! Oracle = a reference implementation of the chosen form, evaluated in-test (no ecosystem consumer   !
! yet, so numerical agreement -- not "it runs" -- is the correctness check).                           !
!==========================================================================================!
program test_plant_respiration
   use meds_test_assert, only : check_close, check_true, test_report
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : pi, t_ref_photo
   use meds_temp_response, only : peaked_arrhenius_scale
   use meds_plant_types, only : wood_params_t, root_params_t
   use meds_plant_respiration, only : stem_maintenance_respiration, fine_root_maintenance_respiration,     &
                                     root_zone_temp_scale
   implicit none
   real(wp), parameter :: AGF  = 0.7_wp   !< the PFT aboveground fraction, once (issue #128)
   real(wp), parameter :: SRF  = 0.06_wp  !< [umol/m2 stem/s @25C] the PFT stem baseline
   real(wp), parameter :: RRF  = 0.30_wp  !< [umol/kgC root/s @25C] the PFT root baseline
   logical,  parameter :: WOODY = .true.


   call test_grass_zero()
   call test_stem_identity_25c()
   call test_tresponse()
   call test_wai_branch()
   call test_aboveground_frac_scales()
   call test_size_scaler()
   call test_root()
   call test_pft_respiration_traits()

   call test_report('test_plant_respiration')

contains



   !----- Reference per-plant stem surface area (the ED2 form, written independently). -----!
   pure real(wp) function stem_area_ref(dbh, height, wai, nplant, aboveground_frac) result(a)
      real(wp), intent(in) :: dbh, height, wai, nplant, aboveground_frac
      a = ( pi * (dbh * 1.0e-2_wp) * height + pi * wai / nplant ) / aboveground_frac
   end function stem_area_ref

   !----- 1. Grass has no stem. ------------------------------------------------------------!
   subroutine test_grass_zero()
      type(wood_params_t) :: p
      real(wp)            :: stem_resp
      !----- is_woody is a per-PFT TRAIT now, passed per cohort, not a field on wood_params_t. --!
      call stem_maintenance_respiration(300.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, .false., SRF, &
                                        p, stem_resp)
      call check_true('grass: stem_resp == 0', stem_resp == 0.0_wp)
   end subroutine test_grass_zero

   !----- 2. At 25 degC the peaked scale is 1 => R = factor25 * stem_area. ------------------!
   subroutine test_stem_identity_25c()
      type(wood_params_t) :: p
      real(wp) :: expect, stem_resp
      p%stem_resp_size_scaler = 0.0_wp
      call stem_maintenance_respiration(t_ref_photo, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, AGF, WOODY, SRF, p, stem_resp)
      expect = 0.06_wp * stem_area_ref(20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, 0.7_wp)
      call check_close('stem @25C = factor25 * stem_area', stem_resp, expect)
   end subroutine test_stem_identity_25c

   !----- 3. Temperature response: at T /= 25C, R = factor25 * peaked_scale(T) * stem_area. -!
   subroutine test_tresponse()
      type(wood_params_t) :: p
      real(wp) :: expect, tscale, stem_resp
      p%stem_resp_size_scaler = 0.0_wp
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, AGF, WOODY, SRF, p, stem_resp)
      tscale = peaked_arrhenius_scale(1.0_wp, p%ea, p%hd, p%ds, 305.0_wp)
      expect = 0.06_wp * tscale * stem_area_ref(20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, 0.7_wp)
      call check_close('stem T-response = factor25 * peaked(T) * area', stem_resp, expect)
   end subroutine test_tresponse

   !----- 4. The WAI branch term adds pi*wai/nplant/aboveground_frac to the per-plant stem area. ------!
   subroutine test_wai_branch()
      type(wood_params_t) :: p
      real(wp) :: tscale, expect_delta, r0, r1
      p%stem_resp_size_scaler = 0.0_wp
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, AGF, WOODY, SRF, p, r0)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY, SRF, p, r1)
      tscale = peaked_arrhenius_scale(1.0_wp, p%ea, p%hd, p%ds, 305.0_wp)
      expect_delta = 0.06_wp * tscale * ( pi * 1.0_wp / 0.1_wp / 0.7_wp )
      call check_close('WAI adds the branch-area term', r1 - r0, expect_delta)
   end subroutine test_wai_branch

   !----- 4b. The aboveground fraction is a PER-COHORT trait and must scale the stem area.    !
   !          Before issue #128 this was a run-uniform 0.7 on wood_params_t, so two PFTs with  !
   !          different allocation respired identically while demography used their real       !
   !          values. The kernel divides the stem area by it, so halving the fraction must      !
   !          double the respiration EXACTLY -- an inverse-proportionality, not just a change.   !
   subroutine test_aboveground_frac_scales()
      type(wood_params_t) :: p
      real(wp) :: r_high, r_low
      p%stem_resp_size_scaler = 0.0_wp
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.7_wp,  WOODY, SRF, p, r_high)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, 0.35_wp, WOODY, SRF, p, r_low)
      call check_close('aboveground_frac scales stem respiration inversely', r_low, 2.0_wp * r_high)
   end subroutine test_aboveground_frac_scales

   !----- 5. The DBH size scaler multiplies the baseline by 10^(scaler*dbh). ----------------!
   subroutine test_size_scaler()
      type(wood_params_t) :: p0, p1
      real(wp) :: r0, r1
      p0%stem_resp_size_scaler = 0.0_wp
      p1%stem_resp_size_scaler = 0.0041_wp
      call stem_maintenance_respiration(t_ref_photo, 50.0_wp, 25.0_wp, 0.0_wp, 0.05_wp, AGF, WOODY, SRF, p0, r0)
      call stem_maintenance_respiration(t_ref_photo, 50.0_wp, 25.0_wp, 0.0_wp, 0.05_wp, AGF, WOODY, SRF, p1, r1)
      call check_close('size scaler = x 10^(scaler*dbh)', r1, r0 * 10.0_wp**(0.0041_wp*50.0_wp))
   end subroutine test_size_scaler

   !----- 6. Fine root: identity at 25C, zero at broot=0, monotone below the optimum. The rate  !
   !         is the PFT's `root_resp_factor25`, passed per cohort now, so the expected value    !
   !         is written against the same constant the call uses -- not a literal that can drift. !
   subroutine test_root()
      type(root_params_t) :: p
      real(wp) :: r_hot, r_cold, r_zero, f_hot, f_cold
      real(wp) :: t_uniform(4), t_split(4), frac(4), f_mean, f_layer
      f_hot  = root_zone_temp_scale([t_ref_photo], [1.0_wp], 1_ik, p)
      f_cold = root_zone_temp_scale([t_ref_photo - 10.0_wp], [1.0_wp], 1_ik, p)
      call fine_root_maintenance_respiration(f_hot, 2.0_wp, RRF, r_hot)
      call check_close('root @25C = factor25 * broot', r_hot, RRF * 2.0_wp)
      call fine_root_maintenance_respiration(f_hot, 0.0_wp, RRF, r_zero)
      call check_true('root: broot == 0 => 0', r_zero == 0.0_wp)
      call fine_root_maintenance_respiration(f_cold, 2.0_wp, RRF, r_cold)
      call check_true('root: warmer (below optimum) respires more', r_hot > r_cold)

      !----- #178: the response is summed OVER LAYERS, not taken at a mean temperature. A       !
      !      UNIFORM profile must give exactly the old answer (so the change is inert where      !
      !      there is no gradient), and a SPLIT profile with the same mean must NOT -- which is  !
      !      the whole point, and what a mean-temperature formulation cannot express.            !
      frac      = [0.4_wp, 0.3_wp, 0.2_wp, 0.1_wp]
      t_uniform = [285.0_wp, 285.0_wp, 285.0_wp, 285.0_wp]
      f_layer = root_zone_temp_scale(t_uniform, frac, 4_ik, p)
      f_mean  = root_zone_temp_scale([sum(frac*t_uniform)], [1.0_wp], 1_ik, p)
      call check_close('uniform profile: layered == mean-temperature', f_layer, f_mean)
      !----- Same root-weighted MEAN (285 K), spread over a 20 K range. Below the optimum the    !
      !      peaked response is convex, so the layered sum must come out HIGHER. ----------------!
      t_split = [275.0_wp, 295.0_wp, 275.0_wp, 315.0_wp]
      call check_close('the split profile has the same weighted mean', sum(frac*t_split), 285.0_wp)
      f_layer = root_zone_temp_scale(t_split, frac, 4_ik, p)
      call check_true('a temperature GRADIENT changes the response at fixed mean',              &
                      abs(f_layer - f_mean) > 1.0e-3_wp * f_mean, f_layer / f_mean)
      call check_true('convex below the optimum => layered exceeds mean-temperature',           &
                      f_layer > f_mean, f_layer / f_mean)
   end subroutine test_root

   !----- 7. The three non-leaf respiration traits are PER-PFT and must differentiate cohorts.  !
   !         They were run-uniform literals in the fast driver until now, so two PFTs respired   !
   !         their stems and roots identically however they differed. Same shape as the          !
   !         aboveground_frac check above: assert the trait CHANGES the answer, not merely that   !
   !         the argument compiles.                                                                !
   subroutine test_pft_respiration_traits()
      type(wood_params_t) :: pw
      type(root_params_t) :: pr
      real(wp) :: r_a, r_b
      pw%stem_resp_size_scaler = 0.0_wp
      !----- stem baseline is a plain multiplier: double it, double the flux. ------------------!
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        SRF, pw, r_a)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, WOODY,      &
                                        2.0_wp * SRF, pw, r_b)
      call check_close('stem_resp_factor25 is per-PFT and scales the flux', r_b, 2.0_wp * r_a)
      !----- and a NON-woody PFT respires no stem at all, whatever its baseline. ----------------!
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, AGF, .false.,    &
                                        2.0_wp * SRF, pw, r_b)
      call check_true('is_woody is per-PFT: a grass cohort respires no stem', r_b == 0.0_wp)
      !----- root baseline likewise. -----------------------------------------------------------!
      call fine_root_maintenance_respiration(root_zone_temp_scale([t_ref_photo], [1.0_wp], 1_ik, pr),  &
                                             2.0_wp, RRF,          r_a)
      call fine_root_maintenance_respiration(root_zone_temp_scale([t_ref_photo], [1.0_wp], 1_ik, pr),  &
                                             2.0_wp, 2.0_wp * RRF, r_b)
      call check_close('root_resp_factor25 is per-PFT and scales the flux', r_b, 2.0_wp * r_a)
   end subroutine test_pft_respiration_traits

end program test_plant_respiration
