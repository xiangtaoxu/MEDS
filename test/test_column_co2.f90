!==========================================================================================!
! test_column_co2 -- unit tests for the P0 column CO2 balance (the well-mixed canopy-air-space  !
! CO2 box, meds_cas_biophysics). Mirrors the design's test plan (section 7):                          !
!   1. RESID ~ 0 : the closed CO2 budget residual vanishes for a spread of inputs.                 !
!   2. STEADY STATE : f_bio = 0 and can_co2 = co2_atm  =>  no change, loss2atm = 0.                  !
!   3. ATM RELAXATION : f_bio = 0, can_co2 /= co2_atm  =>  L-stable relaxation toward co2_atm.        !
!   4. STEADY Ca : constant source  =>  can_co2 -> co2_atm + f_bio/gatm_co2 (analytic fixed point).    !
!   5. SIGN discipline : GPP-only lowers can_co2 (nep>0); respiration-only raises it (nee>0, vents).    !
!   6. CONSERVATION identity : d(storage) = dt*(nee - loss2atm), assembled independently.               !
!   7. AGGREGATE units : a hand-built 3-cohort patch matches the analytic umol/m2/s sums.               !
!   8. HETEROTROPHIC response : Q10 doubling, ED2 cap, moisture hump, zero pool.                         !
!   9. NEP identity : nep = -nee exactly.                                                                 !
!  10-13. DAMM (P1) : Harvard-Forest hand value ~2.15, moisture unimodality, Arrhenius Vmax, anoxia limit. !
!==========================================================================================!
program test_column_co2
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, kgCday_2_umols, r_gas_kj
   use meds_biophysics_types, only : column_co2_budget_t, cohort_co2_flux_t
   use meds_biogeochem_types, only : co2_opts_t, HR_Q10, HR_EXP_ED2, HR_DAMM
   use meds_cas_biophysics,   only : canopy_air_co2_update, aggregate_cohort_co2_fluxes
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_resid_zero()
   call test_steady_state()
   call test_atm_relaxation()
   call test_steady_ca()
   call test_sign_discipline()
   call test_conservation_identity()
   call test_aggregate_units()
   call test_heterotrophic()
   call test_nep_identity()
   call test_damm_hand_value()
   call test_damm_moisture_unimodality()
   call test_damm_arrhenius()
   call test_damm_anoxia_limit()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_co2: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_co2: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (', val, ')'
      else
         nfail = nfail + 1_ik ; print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

   !----- Dry-air molar CAS capacity, recomputed independently of the kernel. ----------------!
   pure function ccapcan_of(rho_air, can_shv, can_depth) result(ccapcan)
      real(wp), intent(in) :: rho_air, can_shv, can_depth
      real(wp) :: ccapcan
      ccapcan = rho_air * (1.0_wp - can_shv) / mmdry * can_depth
   end function ccapcan_of

   !----- 1. Closed-budget residual ~ 0 over a spread of states/fluxes/steps. ------------------!
   subroutine test_resid_zero()
      type(column_co2_budget_t) :: b
      real(wp), parameter :: cco2(6) = [400.0_wp, 380.0_wp, 450.0_wp, 500.0_wp, 350.0_wp, 420.0_wp]
      real(wp), parameter :: gpp(6)  = [ 20.0_wp,   5.0_wp,  30.0_wp,   0.0_wp,  12.0_wp,  18.0_wp]
      real(wp), parameter :: pr(6)   = [  6.0_wp,   4.0_wp,   9.0_wp,   2.0_wp,   5.0_wp,   7.0_wp]
      real(wp), parameter :: hr(6)   = [  3.0_wp,   2.0_wp,   5.0_wp,   1.0_wp,   4.0_wp,   2.5_wp]
      real(wp), parameter :: us(6)   = [ 0.30_wp,  0.10_wp,  0.50_wp,  0.05_wp,  0.20_wp,  0.40_wp]
      real(wp), parameter :: rho(6)  = [  1.2_wp,  1.15_wp,  1.25_wp,   1.1_wp,   1.3_wp,  1.18_wp]
      real(wp), parameter :: shv(6)  = [ 0.00_wp,  0.01_wp,  0.02_wp, 0.005_wp, 0.015_wp,  0.00_wp]
      real(wp), parameter :: dts(6)  = [ 60.0_wp, 300.0_wp,  30.0_wp,1800.0_wp, 120.0_wp, 900.0_wp]
      real(wp)    :: cc, worst
      integer(ik) :: i
      print '(a)', 'test_resid_zero:'
      worst = 0.0_wp
      do i = 1_ik, 6_ik
         cc = cco2(i)
         call canopy_air_co2_update(cc, 20.0_wp, shv(i), gpp(i), pr(i), hr(i), us(i), 1.0_wp,          &
                                    400.0_wp, rho(i), dts(i), b)
         worst = max(worst, abs(b%resid) / max(abs(b%storage), 1.0_wp))
      end do
      call check_true('CO2 budget residual ~ 0 (relative)', worst < 1.0e-9_wp, worst)
   end subroutine test_resid_zero

   !----- 2. f_bio = 0 and can_co2 = co2_atm  =>  no drift, loss2atm = 0. ----------------------!
   subroutine test_steady_state()
      type(column_co2_budget_t) :: b
      real(wp) :: cc
      print '(a)', 'test_steady_state:'
      cc = 400.0_wp
      call canopy_air_co2_update(cc, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,          &
                                 400.0_wp, 1.2_wp, 1800.0_wp, b)
      call check('steady can_co2 unchanged', cc, 400.0_wp, 1.0e-9_wp)
      call check('steady loss2atm = 0', b%loss2atm, 0.0_wp, 1.0e-9_wp)
   end subroutine test_steady_state

   !----- 3. f_bio = 0, can_co2 /= co2_atm  =>  relaxes toward co2_atm; L-stable (no overshoot). !
   subroutine test_atm_relaxation()
      type(column_co2_budget_t) :: b
      real(wp) :: cc_mod, cc_big
      print '(a)', 'test_atm_relaxation:'
      cc_mod = 500.0_wp
      call canopy_air_co2_update(cc_mod, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,      &
                                 400.0_wp, 1.2_wp, 60.0_wp, b)
      call check_true('moderate dt: relaxes toward atm (400 < cc < 500)',                       &
                      cc_mod < 500.0_wp .and. cc_mod > 400.0_wp, cc_mod)
      cc_big = 500.0_wp
      call canopy_air_co2_update(cc_big, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,      &
                                 400.0_wp, 1.2_wp, 1.0e9_wp, b)
      call check_true('huge dt: L-stable, no overshoot below atm', cc_big >= 400.0_wp, cc_big)
      call check('huge dt: converges to atm', cc_big, 400.0_wp, 1.0e-3_wp)
      !----- BUG8: the M-O scalar-transfer coefficient temp2 (c3) scales the atm<->CAS conductance. !
      !      Halving it must SLOW the approach to co2_atm, leaving the CAS farther from the atm      !
      !      after one step. On the old hardcoded c3=1 both runs would be identical.  ---------------!
      block
         real(wp) :: cc_full, cc_half
         cc_full = 500.0_wp
         call canopy_air_co2_update(cc_full, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp, &
                                    400.0_wp, 1.2_wp, 60.0_wp, b)
         cc_half = 500.0_wp
         call canopy_air_co2_update(cc_half, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 0.5_wp, &
                                    400.0_wp, 1.2_wp, 60.0_wp, b)
         call check_true('smaller temp2 (c3) slows atm coupling (CAS stays farther from atm)',        &
                         (cc_half - 400.0_wp) > (cc_full - 400.0_wp), cc_half - cc_full)
      end block
   end subroutine test_atm_relaxation

   !----- 4. Constant source  =>  can_co2 -> co2_atm + f_bio/gatm_co2 (analytic fixed point). ---!
   subroutine test_steady_ca()
      type(column_co2_budget_t) :: b
      real(wp), parameter :: rho_air = 1.2_wp, ustar = 0.30_wp, hetero = 5.0_wp, co2_atm = 400.0_wp
      real(wp) :: cc, gatm_co2, expect
      integer(ik) :: step
      print '(a)', 'test_steady_ca:'
      gatm_co2 = rho_air / mmdry * ustar                     ! can_shv = 0 => can_dmol = rho/mmdry
      expect   = co2_atm + hetero / gatm_co2                 ! f_bio = hetero (gpp = plant_resp = 0)
      cc = co2_atm
      do step = 1_ik, 4000_ik
         call canopy_air_co2_update(cc, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, hetero, ustar, 1.0_wp,         &
                                    co2_atm, rho_air, 60.0_wp, b)
      end do
      call check('steady canopy CO2 = atm + f_bio/gatm', cc, expect, 1.0e-4_wp)
      call check('steady loss2atm balances the source', b%loss2atm, hetero, 1.0e-4_wp)
   end subroutine test_steady_ca

   !----- 5. GPP-only draws CO2 down (nep>0); respiration-only pushes it up (nee>0, vents up). --!
   subroutine test_sign_discipline()
      type(column_co2_budget_t) :: b
      real(wp) :: cc_gpp, cc_resp
      print '(a)', 'test_sign_discipline:'
      cc_gpp = 400.0_wp
      call canopy_air_co2_update(cc_gpp, 20.0_wp, 0.0_wp, 15.0_wp, 0.0_wp, 0.0_wp, 0.10_wp, 1.0_wp,     &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check_true('GPP-only pulls can_co2 below atm', cc_gpp < 400.0_wp, cc_gpp)
      call check_true('GPP-only: nep > 0 (uptake)', b%nep > 0.0_wp, b%nep)
      cc_resp = 400.0_wp
      call canopy_air_co2_update(cc_resp, 20.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 6.0_wp, 0.10_wp, 1.0_wp,     &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check_true('respiration-only pushes can_co2 above atm', cc_resp > 400.0_wp, cc_resp)
      call check_true('respiration-only: nee > 0 (source)', b%nee > 0.0_wp, b%nee)
      call check_true('respiration-only: loss2atm > 0 (vents up)', b%loss2atm > 0.0_wp, b%loss2atm)
   end subroutine test_sign_discipline

   !----- 6. d(storage) = dt*(nee - loss2atm), assembled independently of the kernel's resid. --!
   subroutine test_conservation_identity()
      type(column_co2_budget_t) :: b
      real(wp), parameter :: rho_air = 1.2_wp, can_shv = 0.01_wp, can_depth = 20.0_wp, dt = 600.0_wp
      real(wp) :: cc, ccapcan, storage_before, actual_delta, expect_delta
      print '(a)', 'test_conservation_identity:'
      cc = 430.0_wp
      ccapcan        = ccapcan_of(rho_air, can_shv, can_depth)
      storage_before = ccapcan * cc
      call canopy_air_co2_update(cc, can_depth, can_shv, 22.0_wp, 7.0_wp, 4.0_wp, 0.25_wp, 1.0_wp,      &
                                 400.0_wp, rho_air, dt, b)
      actual_delta = b%storage - storage_before
      expect_delta = dt * (b%nee - b%loss2atm)
      call check('d(storage) = dt*(nee - loss2atm)', actual_delta, expect_delta, 1.0e-6_wp)
   end subroutine test_conservation_identity

   !----- 7. Hand-built 3-cohort patch: leaf-area-weighted GPP/Rd, nplant-weighted stem/root. --!
   subroutine test_aggregate_units()
      type(cohort_co2_flux_t) :: coh
      real(wp), parameter :: np(3)  = [0.10_wp, 0.20_wp, 0.05_wp]   ! plant/m2
      real(wp), parameter :: la(3)  = [5.00_wp, 3.00_wp, 8.00_wp]   ! m2 leaf/plant
      real(wp), parameter :: ag(3)  = [10.0_wp, 8.00_wp, 12.0_wp]   ! umol/m2 leaf/s
      real(wp), parameter :: rd(3)  = [1.00_wp, 0.80_wp, 1.20_wp]   ! umol/m2 leaf/s
      real(wp), parameter :: sr(3)  = [2.00_wp, 1.50_wp, 3.00_wp]   ! umol/plant/s
      real(wp), parameter :: rr(3)  = [1.00_wp, 0.70_wp, 1.50_wp]   ! umol/plant/s
      print '(a)', 'test_aggregate_units:'
      call aggregate_cohort_co2_fluxes(3_ik, np, la, ag, rd, sr, rr, coh)
      call check('GPP  = SUM a_gross*leaf_area*nplant', coh%gross_primary_prod, 14.60_wp, 1.0e-10_wp)
      call check('Rd   = SUM rd*leaf_area*nplant',      coh%leaf_respiration,    1.46_wp, 1.0e-10_wp)
      call check('stem = SUM stem_resp*nplant',         coh%stem_respiration,    0.65_wp, 1.0e-10_wp)
      call check('root = SUM root_resp*nplant',         coh%root_respiration,   0.315_wp, 1.0e-10_wp)
   end subroutine test_aggregate_units

   !----- 8. Rh: Q10 doubling per 10 K, ED2 cap <= 1, moisture hump, zero pool. ----------------!
   subroutine test_heterotrophic()
      type(co2_opts_t) :: opts
      real(wp) :: rh_lo, rh_hi, rh_zero, fw_dry, fw_opt, fw_wet, ft_cap
      real(wp), parameter :: th_opt = 0.8938_wp     ! rel = opt (theta_dry=0, theta_sat=1) => f_water = 1
      print '(a)', 'test_heterotrophic:'
      !----- Q10 doubling: at rh_q10 = 2, +10 K doubles Rh (f_water = 1 at rel = opt). --------!
      opts%hr_model = HR_Q10 ; opts%rh_q10 = 2.0_wp ; opts%rh_t_ref = 288.15_wp ; opts%rh_k_base = 1.0_wp
      rh_lo = heterotrophic_respiration_flux(1.0_wp, 288.15_wp,          th_opt, 0.0_wp, 1.0_wp, opts)
      rh_hi = heterotrophic_respiration_flux(1.0_wp, 288.15_wp + 10.0_wp, th_opt, 0.0_wp, 1.0_wp, opts)
      call check('Q10=2: Rh doubles per +10 K', rh_hi / rh_lo, 2.0_wp, 1.0e-9_wp)
      call check('Q10 baseline Rh = pool*k*kgCday_2_umols', rh_lo, kgCday_2_umols, 1.0e-6_wp)
      !----- Zero pool => zero flux. ----------------------------------------------------------!
      rh_zero = heterotrophic_respiration_flux(0.0_wp, 300.0_wp, th_opt, 0.0_wp, 1.0_wp, opts)
      call check('zero soil-C pool => Rh = 0', rh_zero, 0.0_wp, 1.0e-12_wp)
      !----- ED2 capped exponential: f_temp <= 1 even above the saturation temperature. --------!
      opts%hr_model = HR_EXP_ED2
      ft_cap = heterotrophic_respiration_flux(1.0_wp, 330.0_wp, th_opt, 0.0_wp, 1.0_wp, opts)   ! T > resp_temp_ref
      call check_true('ED2 f_temp capped: Rh <= pool*k*kgCday_2_umols', ft_cap <= kgCday_2_umols + 1.0e-9_wp, ft_cap)
      !----- Moisture hump: f_water peaks at the optimum, falls on both sides. -----------------!
      opts%hr_model = HR_Q10 ; opts%rh_q10 = 1.0_wp    ! q10=1 => f_temp=1, so Rh tracks f_water
      fw_dry = heterotrophic_respiration_flux(1.0_wp, 288.15_wp, 0.20_wp, 0.0_wp, 1.0_wp, opts)
      fw_opt = heterotrophic_respiration_flux(1.0_wp, 288.15_wp, th_opt, 0.0_wp, 1.0_wp, opts)
      fw_wet = heterotrophic_respiration_flux(1.0_wp, 288.15_wp, 0.99_wp, 0.0_wp, 1.0_wp, opts)
      call check_true('moisture hump: Rh(opt) > Rh(dry)', fw_opt > fw_dry, fw_opt - fw_dry)
      call check_true('moisture hump: Rh(opt) > Rh(wet)', fw_opt > fw_wet, fw_opt - fw_wet)
   end subroutine test_heterotrophic

   !----- 9. nep = -nee exactly. --------------------------------------------------------------!
   subroutine test_nep_identity()
      type(column_co2_budget_t) :: b
      real(wp) :: cc
      print '(a)', 'test_nep_identity:'
      cc = 410.0_wp
      call canopy_air_co2_update(cc, 20.0_wp, 0.0_wp, 17.0_wp, 6.0_wp, 3.0_wp, 0.20_wp, 1.0_wp,         &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check('nep = -nee', b%nep, -b%nee, 1.0e-12_wp)
   end subroutine test_nep_identity

   !----- 10. DAMM at the Harvard-Forest point reproduces the published/hand-checked Rh. -------!
   !      (soil_temp=15C, theta=0.229, porosity=0.6825, pool=4.8 kgC/m2, depth=10 cm) => ~2.15.  !
   subroutine test_damm_hand_value()
      type(co2_opts_t) :: opts
      real(wp) :: rh
      print '(a)', 'test_damm_hand_value:'
      opts%hr_model = HR_DAMM                              ! default damm params = Davidson-2012
      rh = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.229_wp, 0.0_wp, 0.6825_wp, opts)
      call check('DAMM Harvard-Forest Rh ~ 2.15 umol/m2/s', rh, 2.147_wp, 0.05_wp)
   end subroutine test_damm_hand_value

   !----- 11. DAMM moisture response is UNIMODAL: rises then falls, ~0 at both ends. -----------!
   subroutine test_damm_moisture_unimodality()
      type(co2_opts_t) :: opts
      real(wp), parameter :: ts = 0.6825_wp
      real(wp) :: rh_dry, rh_peak, rh_wet, rh_sat
      print '(a)', 'test_damm_moisture_unimodality:'
      opts%hr_model = HR_DAMM
      rh_dry  = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.05_wp,   0.0_wp, ts, opts)
      rh_peak = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.35_wp,   0.0_wp, ts, opts)
      rh_wet  = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.60_wp,   0.0_wp, ts, opts)
      rh_sat  = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, ts,        0.0_wp, ts, opts)
      call check_true('DAMM hump: Rh(mid) > Rh(dry) [substrate-limited]', rh_peak > rh_dry, rh_peak - rh_dry)
      call check_true('DAMM hump: Rh(mid) > Rh(wet) [O2-limited]',        rh_peak > rh_wet, rh_peak - rh_wet)
      call check('DAMM: Rh -> 0 at saturation (anoxia, finite)', rh_sat, 0.0_wp, 1.0e-12_wp)
   end subroutine test_damm_moisture_unimodality

   !----- 12. DAMM temperature response is exactly Arrhenius in Vmax (Ea/R pairing check). -----!
   subroutine test_damm_arrhenius()
      type(co2_opts_t) :: opts
      real(wp), parameter :: t1 = 288.15_wp, t2 = 298.15_wp, ea = 72.26_wp
      real(wp) :: rh1, rh2, ratio_expect
      print '(a)', 'test_damm_arrhenius:'
      opts%hr_model = HR_DAMM                              ! default ea_sx = 72.26 kJ/mol
      rh1 = heterotrophic_respiration_flux(4.8_wp, t1, 0.40_wp, 0.0_wp, 0.6825_wp, opts)
      rh2 = heterotrophic_respiration_flux(4.8_wp, t2, 0.40_wp, 0.0_wp, 0.6825_wp, opts)
      ratio_expect = exp((ea / r_gas_kj) * (1.0_wp / t1 - 1.0_wp / t2))    ! only Vmax depends on T
      call check('DAMM Rh ratio = Arrhenius factor', rh2 / rh1, ratio_expect, 1.0e-9_wp * ratio_expect)
   end subroutine test_damm_arrhenius

   !----- 13. DAMM O2 anoxia limit: Rh -> 0 as theta -> porosity, and stays 0 (no NaN) above. --!
   subroutine test_damm_anoxia_limit()
      type(co2_opts_t) :: opts
      real(wp), parameter :: ts = 0.6825_wp
      real(wp) :: rh_near, rh_at, rh_over
      print '(a)', 'test_damm_anoxia_limit:'
      opts%hr_model = HR_DAMM
      rh_near = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.66_wp,      0.0_wp, ts, opts)
      rh_at   = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, ts,          0.0_wp, ts, opts)
      rh_over = heterotrophic_respiration_flux(4.8_wp, 288.15_wp, 0.70_wp,      0.0_wp, ts, opts)
      call check_true('anoxia: Rh(near-sat) < Rh at optimum band', rh_near < 0.5_wp, rh_near)
      call check('anoxia: Rh = 0 at saturation', rh_at, 0.0_wp, 1.0e-12_wp)
      call check_true('anoxia: theta > porosity clamps to Rh = 0 (finite, not NaN)',            &
                      rh_over == 0.0_wp .and. rh_over == rh_over, rh_over)
   end subroutine test_damm_anoxia_limit

end program test_column_co2
