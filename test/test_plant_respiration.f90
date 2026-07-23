!==========================================================================================!
! test_plant_respiration -- unit tests for the non-leaf respiration kernels.                 !
!                                                                                          !
!   1. GRASS      : is_woody = .false. => stem respiration is identically zero.               !
!   2. STEM @25C  : at 25 degC the peaked scale is 1, so R = factor25 * stem_area (cylinder).   !
!   3. T-RESP     : at T /= 25 degC, R = factor25 * peaked_scale(T) * stem_area.                    !
!   4. WAI        : the WAI branch term adds pi*wai/nplant/agf_bs to the per-plant stem area.       !
!   5. SIZE       : the DBH size scaler multiplies the baseline by 10^(scaler*dbh).                  !
!   6. ROOT       : R = factor25 * broot at 25 degC; broot = 0 => 0; monotone up to the optimum.     !
!   (Growth/construction respiration moved to meds_plant_carbon_allocation; see test_plant_carbon_allocation.)!
! Oracle = a reference implementation of the chosen form, evaluated in-test (no ecosystem consumer   !
! yet, so numerical agreement -- not "it runs" -- is the correctness check).                           !
!==========================================================================================!
program test_plant_respiration
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : pi, t_ref_photo
   use meds_temp_response, only : peaked_arrhenius_scale
   use meds_plant_interface, only : wood_params_t, root_params_t,                               &
                                    stem_maintenance_respiration,                               &
                                    fine_root_maintenance_respiration
   implicit none

   integer(ik) :: nfail
   nfail = 0_ik

   call test_grass_zero()
   call test_stem_identity_25c()
   call test_tresponse()
   call test_wai_branch()
   call test_size_scaler()
   call test_root()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_respiration: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_respiration: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check_true(name, cond)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a)', '  FAIL : ', name
      end if
   end subroutine check_true

   subroutine check_close(name, got, expect)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect
      real(wp)                     :: tol
      tol = 1.0e-9_wp * max(1.0_wp, abs(expect))
      if (abs(got - expect) <= tol) then
         print '(a,a,a,es13.6,a,es13.6,a)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.6,a,es13.6)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check_close

   !----- Reference per-plant stem surface area (the ED2 form, written independently). -----!
   pure real(wp) function stem_area_ref(dbh, height, wai, nplant, agf_bs) result(a)
      real(wp), intent(in) :: dbh, height, wai, nplant, agf_bs
      a = ( pi * (dbh * 1.0e-2_wp) * height + pi * wai / nplant ) / agf_bs
   end function stem_area_ref

   !----- 1. Grass has no stem. ------------------------------------------------------------!
   subroutine test_grass_zero()
      type(wood_params_t) :: p
      real(wp)            :: stem_resp
      p%is_woody = .false. ; p%stem_resp_factor25 = 0.06_wp
      call stem_maintenance_respiration(300.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, p, stem_resp)
      call check_true('grass: stem_resp == 0', stem_resp == 0.0_wp)
   end subroutine test_grass_zero

   !----- 2. At 25 degC the peaked scale is 1 => R = factor25 * stem_area. ------------------!
   subroutine test_stem_identity_25c()
      type(wood_params_t) :: p
      real(wp) :: expect, stem_resp
      p%is_woody = .true. ; p%stem_resp_factor25 = 0.06_wp ; p%stem_resp_size_scaler = 0.0_wp
      p%agf_bs = 0.7_wp
      call stem_maintenance_respiration(t_ref_photo, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, p, stem_resp)
      expect = 0.06_wp * stem_area_ref(20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, 0.7_wp)
      call check_close('stem @25C = factor25 * stem_area', stem_resp, expect)
   end subroutine test_stem_identity_25c

   !----- 3. Temperature response: at T /= 25C, R = factor25 * peaked_scale(T) * stem_area. -!
   subroutine test_tresponse()
      type(wood_params_t) :: p
      real(wp) :: expect, tscale, stem_resp
      p%is_woody = .true. ; p%stem_resp_factor25 = 0.06_wp ; p%stem_resp_size_scaler = 0.0_wp ; p%agf_bs = 0.7_wp
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, p, stem_resp)
      tscale = peaked_arrhenius_scale(1.0_wp, p%ea, p%hd, p%ds, 305.0_wp)
      expect = 0.06_wp * tscale * stem_area_ref(20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, 0.7_wp)
      call check_close('stem T-response = factor25 * peaked(T) * area', stem_resp, expect)
   end subroutine test_tresponse

   !----- 4. The WAI branch term adds pi*wai/nplant/agf_bs to the per-plant stem area. ------!
   subroutine test_wai_branch()
      type(wood_params_t) :: p
      real(wp) :: tscale, expect_delta, r0, r1
      p%is_woody = .true. ; p%stem_resp_factor25 = 0.06_wp ; p%stem_resp_size_scaler = 0.0_wp
      p%agf_bs = 0.7_wp
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 0.0_wp, 0.1_wp, p, r0)
      call stem_maintenance_respiration(305.0_wp, 20.0_wp, 15.0_wp, 1.0_wp, 0.1_wp, p, r1)
      tscale = peaked_arrhenius_scale(1.0_wp, p%ea, p%hd, p%ds, 305.0_wp)
      expect_delta = 0.06_wp * tscale * ( pi * 1.0_wp / 0.1_wp / 0.7_wp )
      call check_close('WAI adds the branch-area term', r1 - r0, expect_delta)
   end subroutine test_wai_branch

   !----- 5. The DBH size scaler multiplies the baseline by 10^(scaler*dbh). ----------------!
   subroutine test_size_scaler()
      type(wood_params_t) :: p0, p1
      real(wp) :: r0, r1
      p0%is_woody = .true. ; p0%stem_resp_factor25 = 0.06_wp ; p0%stem_resp_size_scaler = 0.0_wp    ; p0%agf_bs = 0.7_wp
      p1%is_woody = .true. ; p1%stem_resp_factor25 = 0.06_wp ; p1%stem_resp_size_scaler = 0.0041_wp ; p1%agf_bs = 0.7_wp
      call stem_maintenance_respiration(t_ref_photo, 50.0_wp, 25.0_wp, 0.0_wp, 0.05_wp, p0, r0)
      call stem_maintenance_respiration(t_ref_photo, 50.0_wp, 25.0_wp, 0.0_wp, 0.05_wp, p1, r1)
      call check_close('size scaler = x 10^(scaler*dbh)', r1, r0 * 10.0_wp**(0.0041_wp*50.0_wp))
   end subroutine test_size_scaler

   !----- 6. Fine root: identity at 25C, zero at broot=0, monotone below the optimum. -------!
   subroutine test_root()
      type(root_params_t) :: p
      real(wp) :: r_hot, r_cold, r_zero
      p%root_resp_factor25 = 0.25_wp
      call fine_root_maintenance_respiration(t_ref_photo, 2.0_wp, p, r_hot)
      call check_close('root @25C = factor25 * broot', r_hot, 0.25_wp * 2.0_wp)
      call fine_root_maintenance_respiration(t_ref_photo, 0.0_wp, p, r_zero)
      call check_true('root: broot == 0 => 0', r_zero == 0.0_wp)
      call fine_root_maintenance_respiration(t_ref_photo,          2.0_wp, p, r_hot)
      call fine_root_maintenance_respiration(t_ref_photo - 10.0_wp, 2.0_wp, p, r_cold)
      call check_true('root: warmer (below optimum) respires more', r_hot > r_cold)
   end subroutine test_root

end program test_plant_respiration
