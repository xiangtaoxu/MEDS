! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_column_co2 -- unit tests for the column CO2 balance (the well-mixed                 !
! canopy-air-space CO2 box, docs/science/canopy_air_space_biophysics.md eq. (1)). The box  !
! is advanced by the production kernel meds_canopy_air_space::cas_box_step_implicit; the   !
! local co2_step helper feeds it the raw scalar inputs and assembles the budget            !
! diagnostics (nee/nep/loss2atm/storage/resid) from its result. Test plan:                 !
!   1. RESID ~ 0 : the closed CO2 budget residual vanishes for a spread of inputs.                 !
!   2. STEADY STATE : f_bio = 0 and can_co2 = co2_atm  =>  no change, loss2atm = 0.                  !
!   3. ATM RELAXATION : f_bio = 0, can_co2 /= co2_atm  =>  L-stable relaxation toward co2_atm.        !
!   4. STEADY Ca : constant source  =>  can_co2 -> co2_atm + f_bio/gatm_co2 (analytic fixed point).    !
!   5. SIGN discipline : GPP-only lowers can_co2 (nep>0); respiration-only raises it (nee>0, vents).    !
!   6. CONSERVATION identity : d(storage) = dt*(nee - loss2atm), assembled independently.               !
!   7. NEP identity : nep = -nee exactly.                                                                 !
! Heterotrophic respiration is the CENTURY soil-carbon matrix, tested in test_soil_biogeochem and   !
! test_biogeochem_dynamics.                                                                         !
!==========================================================================================!
program test_column_co2
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, kgcday_2_umols
   use meds_canopy_air_space,   only : cas_box_t, cas_source_t, cas_box_step_implicit
   implicit none

   !----- Local CO2 budget diagnostics, derived in co2_step from the production                  !
   !      cas_box_step_implicit result (see below).                                              !
   type :: co2_budget_t
      real(wp) :: nee = 0.0_wp, nep = 0.0_wp, loss2atm = 0.0_wp, storage = 0.0_wp, resid = 0.0_wp
   end type co2_budget_t


   call test_resid_zero()
   call test_steady_state()
   call test_atm_relaxation()
   call test_steady_ca()
   call test_sign_discipline()
   call test_conservation_identity()
   call test_nep_identity()

   call test_report('test_column_co2')

contains



   !----- Dry-air molar CAS capacity, recomputed independently of the kernel. ----------------!
   pure function ccapcan_of(rho_air, can_shv, can_depth) result(ccapcan)
      real(wp), intent(in) :: rho_air, can_shv, can_depth
      real(wp) :: ccapcan
      ccapcan = rho_air * (1.0_wp - can_shv) / mmdry * can_depth
   end function ccapcan_of

   !----- Advance the CO2 twin via the PRODUCTION kernel cas_box_step_implicit, then assemble  !
   !      the budget diagnostics from its result.                                                   !
   !      The enthalpy/vapour twins are held inert (nonzero air_mass_capacity avoids a 0/0).        !
   subroutine co2_step(can_co2, can_depth, can_shv, gpp, plant_resp, hetero,                    &
                       ustar, vapour_transfer_factor, co2_atm, rho_air, dt, budget)
      real(wp), intent(inout) :: can_co2
      real(wp), intent(in)    :: can_depth, can_shv, gpp, plant_resp, hetero
      real(wp), intent(in)    :: ustar, vapour_transfer_factor, co2_atm, rho_air, dt
      type(co2_budget_t), intent(out) :: budget
      type(cas_box_t) :: column
      type(cas_source_t) :: source
      real(wp) :: can_dmol, ccapcan, gatm_co2, f_bio, enth_dummy, shv_dummy, co2_new
      can_dmol = rho_air * (1.0_wp - can_shv) / mmdry        ! dry-air molar density [mol/m3]
      ccapcan  = can_dmol * can_depth                        ! molar CAS capacity   [mol/m2]
      gatm_co2 = can_dmol * ustar * vapour_transfer_factor                    ! atm<->CAS molar conductance
      f_bio    = hetero + plant_resp - gpp                   ! Reco - GPP  [umol/m2/s]
      column%air_mass_capacity   = 1.0_wp                    ! inert enthalpy/vapour twins (avoid 0/0)
      column%air_molar_capacity  = ccapcan
      column%atm_conductance_co2 = gatm_co2
      column%atm_co2             = co2_atm
      source%biotic_co2_source   = f_bio
      call cas_box_step_implicit(0.0_wp, 0.0_wp, can_co2, source, column, dt,                 &
                                    enth_dummy, shv_dummy, co2_new)
      budget%resid    = ccapcan * (co2_new - can_co2) - dt * (f_bio + gatm_co2 * (co2_atm - co2_new))
      budget%nee      = f_bio
      budget%nep      = -f_bio
      budget%loss2atm = gatm_co2 * (co2_new - co2_atm)
      budget%storage  = ccapcan * co2_new
      can_co2 = co2_new
   end subroutine co2_step

   !----- 1. Closed-budget residual ~ 0 over a spread of states/fluxes/steps. ------------------!
   subroutine test_resid_zero()
      type(co2_budget_t) :: b
      real(wp), parameter :: CCO2(6) = [400.0_wp, 380.0_wp, 450.0_wp, 500.0_wp, 350.0_wp, 420.0_wp]
      real(wp), parameter :: GPP(6)  = [ 20.0_wp,   5.0_wp,  30.0_wp,   0.0_wp,  12.0_wp,  18.0_wp]
      real(wp), parameter :: PR(6)   = [  6.0_wp,   4.0_wp,   9.0_wp,   2.0_wp,   5.0_wp,   7.0_wp]
      real(wp), parameter :: HR(6)   = [  3.0_wp,   2.0_wp,   5.0_wp,   1.0_wp,   4.0_wp,   2.5_wp]
      real(wp), parameter :: US(6)   = [ 0.30_wp,  0.10_wp,  0.50_wp,  0.05_wp,  0.20_wp,  0.40_wp]
      real(wp), parameter :: RHO(6)  = [  1.2_wp,  1.15_wp,  1.25_wp,   1.1_wp,   1.3_wp,  1.18_wp]
      real(wp), parameter :: SHV(6)  = [ 0.00_wp,  0.01_wp,  0.02_wp, 0.005_wp, 0.015_wp,  0.00_wp]
      real(wp), parameter :: DTS(6)  = [ 60.0_wp, 300.0_wp,  30.0_wp,1800.0_wp, 120.0_wp, 900.0_wp]
      real(wp)    :: cc, worst
      integer(ik) :: i
      print '(a)', 'test_resid_zero:'
      worst = 0.0_wp
      do i = 1_ik, 6_ik
         cc = CCO2(i)
         call co2_step(cc, 20.0_wp, SHV(i), GPP(i), PR(i), HR(i), US(i), 1.0_wp,          &
                                    400.0_wp, RHO(i), DTS(i), b)
         worst = max(worst, abs(b%resid) / max(abs(b%storage), 1.0_wp))
      end do
      call check_true('CO2 budget residual ~ 0 (relative)', worst < 1.0e-9_wp, worst)
   end subroutine test_resid_zero

   !----- 2. f_bio = 0 and can_co2 = co2_atm  =>  no drift, loss2atm = 0. ----------------------!
   subroutine test_steady_state()
      type(co2_budget_t) :: b
      real(wp) :: cc
      print '(a)', 'test_steady_state:'
      cc = 400.0_wp
      call co2_step(cc, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,          &
                                 400.0_wp, 1.2_wp, 1800.0_wp, b)
      call check('steady can_co2 unchanged', cc, 400.0_wp, 1.0e-9_wp)
      call check('steady loss2atm = 0', b%loss2atm, 0.0_wp, 1.0e-9_wp)
   end subroutine test_steady_state

   !----- 3. f_bio = 0, can_co2 /= co2_atm  =>  relaxes toward co2_atm; L-stable (no overshoot). !
   subroutine test_atm_relaxation()
      type(co2_budget_t) :: b
      real(wp) :: cc_mod, cc_big
      print '(a)', 'test_atm_relaxation:'
      cc_mod = 500.0_wp
      call co2_step(cc_mod, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,      &
                                 400.0_wp, 1.2_wp, 60.0_wp, b)
      call check_true('moderate dt: relaxes toward atm (400 < cc < 500)',                       &
                      cc_mod < 500.0_wp .and. cc_mod > 400.0_wp, cc_mod)
      cc_big = 500.0_wp
      call co2_step(cc_big, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp,      &
                                 400.0_wp, 1.2_wp, 1.0e9_wp, b)
      call check_true('huge dt: L-stable, no overshoot below atm', cc_big >= 400.0_wp, cc_big)
      call check('huge dt: converges to atm', cc_big, 400.0_wp, 1.0e-3_wp)
      !----- The M-O scalar-transfer coefficient vapour_transfer_factor (c3) scales the atm<->CAS conductance        !
      !      (docs/science/canopy_aerodynamics.md eq. (5)). Halving it must SLOW the approach to      !
      !      co2_atm, leaving the CAS farther from the atm after one step; a kernel that ignored     !
      !      vapour_transfer_factor would give two identical runs. ---------------------------------------------------!
      block
         real(wp) :: cc_full, cc_half
         cc_full = 500.0_wp
         call co2_step(cc_full, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 1.0_wp, &
                                    400.0_wp, 1.2_wp, 60.0_wp, b)
         cc_half = 500.0_wp
         call co2_step(cc_half, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.30_wp, 0.5_wp, &
                                    400.0_wp, 1.2_wp, 60.0_wp, b)
         call check_true('smaller vapour_transfer_factor (c3) slows atm coupling (CAS stays farther from atm)',        &
                         (cc_half - 400.0_wp) > (cc_full - 400.0_wp), cc_half - cc_full)
      end block
   end subroutine test_atm_relaxation

   !----- 4. Constant source  =>  can_co2 -> co2_atm + f_bio/gatm_co2 (analytic fixed point). ---!
   subroutine test_steady_ca()
      type(co2_budget_t) :: b
      real(wp), parameter :: RHO_AIR = 1.2_wp, USTAR = 0.30_wp, HETERO = 5.0_wp, CO2_ATM = 400.0_wp
      real(wp) :: cc, gatm_co2, expect
      integer(ik) :: step
      print '(a)', 'test_steady_ca:'
      gatm_co2 = RHO_AIR / mmdry * USTAR                     ! can_shv = 0 => can_dmol = rho/mmdry
      expect   = CO2_ATM + HETERO / gatm_co2                 ! f_bio = hetero (gpp = plant_resp = 0)
      cc = CO2_ATM
      do step = 1_ik, 4000_ik
         call co2_step(cc, 20.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, HETERO, USTAR, 1.0_wp,         &
                                    CO2_ATM, RHO_AIR, 60.0_wp, b)
      end do
      call check('steady canopy CO2 = atm + f_bio/gatm', cc, expect, 1.0e-4_wp)
      call check('steady loss2atm balances the source', b%loss2atm, HETERO, 1.0e-4_wp)
   end subroutine test_steady_ca

   !----- 5. GPP-only draws CO2 down (nep>0); respiration-only pushes it up (nee>0, vents up). --!
   subroutine test_sign_discipline()
      type(co2_budget_t) :: b
      real(wp) :: cc_gpp, cc_resp
      print '(a)', 'test_sign_discipline:'
      cc_gpp = 400.0_wp
      call co2_step(cc_gpp, 20.0_wp, 0.0_wp, 15.0_wp, 0.0_wp, 0.0_wp, 0.10_wp, 1.0_wp,     &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check_true('GPP-only pulls can_co2 below atm', cc_gpp < 400.0_wp, cc_gpp)
      call check_true('GPP-only: nep > 0 (uptake)', b%nep > 0.0_wp, b%nep)
      cc_resp = 400.0_wp
      call co2_step(cc_resp, 20.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 6.0_wp, 0.10_wp, 1.0_wp,     &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check_true('respiration-only pushes can_co2 above atm', cc_resp > 400.0_wp, cc_resp)
      call check_true('respiration-only: nee > 0 (source)', b%nee > 0.0_wp, b%nee)
      call check_true('respiration-only: loss2atm > 0 (vents up)', b%loss2atm > 0.0_wp, b%loss2atm)
   end subroutine test_sign_discipline

   !----- 6. d(storage) = dt*(nee - loss2atm), assembled independently of the kernel's resid. --!
   subroutine test_conservation_identity()
      type(co2_budget_t) :: b
      real(wp), parameter :: RHO_AIR = 1.2_wp, CAN_SHV = 0.01_wp, CAN_DEPTH = 20.0_wp, DT = 600.0_wp
      real(wp) :: cc, ccapcan, storage_before, actual_delta, expect_delta
      print '(a)', 'test_conservation_identity:'
      cc = 430.0_wp
      ccapcan        = ccapcan_of(RHO_AIR, CAN_SHV, CAN_DEPTH)
      storage_before = ccapcan * cc
      call co2_step(cc, CAN_DEPTH, CAN_SHV, 22.0_wp, 7.0_wp, 4.0_wp, 0.25_wp, 1.0_wp,      &
                                 400.0_wp, RHO_AIR, DT, b)
      actual_delta = b%storage - storage_before
      expect_delta = DT * (b%nee - b%loss2atm)
      call check('d(storage) = dt*(nee - loss2atm)', actual_delta, expect_delta, 1.0e-6_wp)
   end subroutine test_conservation_identity


   !----- 7. nep = -nee exactly. --------------------------------------------------------------!
   subroutine test_nep_identity()
      type(co2_budget_t) :: b
      real(wp) :: cc
      print '(a)', 'test_nep_identity:'
      cc = 410.0_wp
      call co2_step(cc, 20.0_wp, 0.0_wp, 17.0_wp, 6.0_wp, 3.0_wp, 0.20_wp, 1.0_wp,         &
                                 400.0_wp, 1.2_wp, 300.0_wp, b)
      call check('nep = -nee', b%nep, -b%nee, 1.0e-12_wp)
   end subroutine test_nep_identity





end program test_column_co2
