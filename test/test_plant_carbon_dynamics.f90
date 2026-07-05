!==========================================================================================!
! test_plant_carbon_dynamics -- unit tests for the carbon budget + allocation kernels.       !
!                                                                                          !
!   1. TURNOVER       : non-evergreen returns the flat baseline rates; evergreen cold-        !
!                       suppresses them (factor = 0.5 at 5 degC, ~1 when warm, monotone).       !
!   2. ON-ALLOMETRY   : pools at target (zero deficits), gain > turnover => turnover is         !
!                       replaced (npp_leaf/fineroot = 0) and the surplus grows wood.            !
!   3. WOOD RESIDUAL  : gain only covers turnover + a leaf deficit => npp_wood = 0.             !
!   4. FLUSH (ON)     : a leaf deficit with short NPP draws storage when PHEN_ON               !
!                       (npp_leaf > 0, npp_nonstructural < 0).                                  !
!   5. NO-FLUSH (OFF) : the same case with PHEN_DORMANT does NOT draw storage for growth.        !
!   6. STORAGE REFILL : surplus fills the storage deficit before wood.                          !
!   7. NEGATIVE NPP   : respiration debt paid from storage (npp_nonstructural = net_carbon),     !
!                       no growth, not starving; storage too small => starving + deficit > 0.     !
!   8. CLOSURE        : (leaf+fineroot+wood+nonstructural) - deficit == net_carbon - turnover,    !
!                       asserted for every scenario above.                                        !
! Oracle = the closure identity + the documented priority ordering (no ecosystem consumer yet).   !
!==========================================================================================!
program test_plant_carbon_dynamics
   use meds_kinds,           only : wp, ik
   use meds_plant_interface, only : turnover_env_t, turnover_params_t, turnover_rates_t,        &
                                    carbon_gain_t, carbon_loss_t, carbon_demand_t, carbon_npp_t, &
                                    tissue_turnover_rates, plant_carbon_allocation
   implicit none

   !----- Phenology status flags (mirror meds_plant_types; re-declared here to keep the test   !
   !       self-contained, exactly as the flags are used at call sites). --------------------!
   integer(ik), parameter :: PHEN_OFF     = -1_ik
   integer(ik), parameter :: PHEN_DORMANT =  0_ik
   integer(ik), parameter :: PHEN_ON      =  1_ik

   integer(ik) :: nfail
   nfail = 0_ik

   call test_turnover()
   call test_on_allometry()
   call test_wood_residual()
   call test_flush_on()
   call test_no_flush_dormant()
   call test_storage_refill()
   call test_negative_npp()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_carbon_dynamics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_carbon_dynamics: ', nfail, ' FAILED'
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

   !----- The carbon-closure identity the kernel guarantees on every call. -----------------!
   subroutine check_closure(name, gain, loss, npp)
      character(len=*),    intent(in) :: name
      type(carbon_gain_t), intent(in) :: gain
      type(carbon_loss_t), intent(in) :: loss
      type(carbon_npp_t),  intent(in) :: npp
      real(wp) :: lhs, rhs
      lhs = npp%leaf + npp%fineroot + npp%wood + npp%nonstructural - npp%deficit
      rhs = gain%net_carbon - (loss%leaf + loss%fineroot)
      call check_close(name, lhs, rhs)
   end subroutine check_closure

   !----- 1. Turnover rates: flat baseline vs evergreen cold suppression. -------------------!
   subroutine test_turnover()
      type(turnover_env_t)    :: env
      type(turnover_params_t) :: pd, pe
      type(turnover_rates_t)  :: rd, r_cold, r_ref, r_warm
      pd%leaf_turnover_rate = 0.5_wp ; pd%fineroot_turnover_rate = 1.2_wp ; pd%evergreen = .false.
      pe%leaf_turnover_rate = 0.5_wp ; pe%fineroot_turnover_rate = 1.2_wp ; pe%evergreen = .true.

      env%tissue_temp = 300.0_wp
      call tissue_turnover_rates(env, pd, rd)
      call check_close('deciduous: leaf rate = baseline',     rd%leaf,     0.5_wp)
      call check_close('deciduous: fineroot rate = baseline', rd%fineroot, 1.2_wp)

      env%tissue_temp = 278.15_wp                         ! 5 degC reference => factor 0.5
      call tissue_turnover_rates(env, pe, r_ref)
      call check_close('evergreen @5C: factor = 0.5 (leaf)', r_ref%leaf, 0.5_wp * 0.5_wp)

      env%tissue_temp = 273.15_wp ; call tissue_turnover_rates(env, pe, r_cold)
      env%tissue_temp = 303.15_wp ; call tissue_turnover_rates(env, pe, r_warm)
      call check_true('evergreen: warmer => more turnover', r_warm%leaf > r_cold%leaf)
      call check_true('evergreen: warm approaches baseline', r_warm%leaf < 0.5_wp .and. r_warm%leaf > 0.45_wp)
   end subroutine test_turnover

   !----- 2. On allometry (zero deficits): turnover replaced, surplus -> wood. --------------!
   subroutine test_on_allometry()
      type(carbon_gain_t)   :: gain
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: npp
      gain = carbon_gain_t(net_carbon=1.0_wp, storage=5.0_wp)
      loss = carbon_loss_t(leaf=0.1_wp, fineroot=0.2_wp)
      dem  = carbon_demand_t(leaf=0.0_wp, fineroot=0.0_wp, storage=0.0_wp, wood=1.0e6_wp)  ! wood = residual sink
      call plant_carbon_allocation(gain, loss, dem, PHEN_ON, npp)
      call check_close('on-allom: npp_leaf = 0 (turnover replaced)',     npp%leaf,     0.0_wp)
      call check_close('on-allom: npp_fineroot = 0',                     npp%fineroot, 0.0_wp)
      call check_close('on-allom: npp_nonstructural = 0',                npp%nonstructural, 0.0_wp)
      call check_close('on-allom: npp_wood = gain - turnover',           npp%wood, 1.0_wp - 0.3_wp)
      call check_true ('on-allom: not starving',                    .not. npp%starving)
      call check_closure('on-allom: carbon closes', gain, loss, npp)
   end subroutine test_on_allometry

   !----- 3. Wood is a residual sink: no surplus => npp_wood = 0. ---------------------------!
   subroutine test_wood_residual()
      type(carbon_gain_t)   :: gain
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: npp
      gain = carbon_gain_t(net_carbon=0.5_wp, storage=0.0_wp)
      loss = carbon_loss_t(leaf=0.1_wp, fineroot=0.0_wp)
      dem  = carbon_demand_t(leaf=0.4_wp, fineroot=0.0_wp, storage=0.0_wp, wood=1.0e6_wp)
      call plant_carbon_allocation(gain, loss, dem, PHEN_ON, npp)
      ! 0.5 NPP: 0.1 replaces leaf turnover, 0.4 fills the leaf deficit -> nothing left for wood.
      call check_close('residual: npp_wood = 0',            npp%wood, 0.0_wp)
      call check_close('residual: npp_leaf = leaf deficit', npp%leaf, 0.4_wp)
      call check_closure('residual: carbon closes', gain, loss, npp)
   end subroutine test_wood_residual

   !----- 4. Flush from storage when PHEN_ON and NPP is short. ------------------------------!
   subroutine test_flush_on()
      type(carbon_gain_t)   :: gain
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: npp
      gain = carbon_gain_t(net_carbon=0.2_wp, storage=3.0_wp)
      loss = carbon_loss_t(leaf=0.0_wp, fineroot=0.0_wp)
      dem  = carbon_demand_t(leaf=1.0_wp, fineroot=0.0_wp, storage=0.0_wp, wood=1.0e6_wp)
      call plant_carbon_allocation(gain, loss, dem, PHEN_ON, npp)
      ! leaf deficit 1.0: 0.2 from NPP + 0.8 drawn from storage.
      call check_close('flush(ON): npp_leaf = full deficit',         npp%leaf, 1.0_wp)
      call check_close('flush(ON): npp_nonstructural = -0.8 (drawn)', npp%nonstructural, -0.8_wp)
      call check_close('flush(ON): npp_wood = 0',                     npp%wood, 0.0_wp)
      call check_closure('flush(ON): carbon closes', gain, loss, npp)
   end subroutine test_flush_on

   !----- 5. No storage draw for growth when dormant. --------------------------------------!
   subroutine test_no_flush_dormant()
      type(carbon_gain_t)   :: gain
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: npp
      gain = carbon_gain_t(net_carbon=0.2_wp, storage=3.0_wp)
      loss = carbon_loss_t(leaf=0.0_wp, fineroot=0.0_wp)
      dem  = carbon_demand_t(leaf=1.0_wp, fineroot=0.0_wp, storage=0.0_wp, wood=1.0e6_wp)
      call plant_carbon_allocation(gain, loss, dem, PHEN_DORMANT, npp)
      ! dormant: growth funded from NPP only (0.2); storage untouched by the flush.
      call check_close('no-flush(DORMANT): npp_leaf = NPP only',       npp%leaf, 0.2_wp)
      call check_close('no-flush(DORMANT): storage untouched',         npp%nonstructural, 0.0_wp)
      call check_closure('no-flush(DORMANT): carbon closes', gain, loss, npp)
   end subroutine test_no_flush_dormant

   !----- 6. Storage refill takes priority over wood. --------------------------------------!
   subroutine test_storage_refill()
      type(carbon_gain_t)   :: gain
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: npp
      gain = carbon_gain_t(net_carbon=1.0_wp, storage=0.0_wp)
      loss = carbon_loss_t(leaf=0.0_wp, fineroot=0.0_wp)
      dem  = carbon_demand_t(leaf=0.0_wp, fineroot=0.0_wp, storage=0.3_wp, wood=1.0e6_wp)
      call plant_carbon_allocation(gain, loss, dem, PHEN_ON, npp)
      call check_close('refill: npp_nonstructural = storage deficit', npp%nonstructural, 0.3_wp)
      call check_close('refill: npp_wood = remainder',                npp%wood, 0.7_wp)
      call check_closure('refill: carbon closes', gain, loss, npp)
   end subroutine test_storage_refill

   !----- 7. Negative NPP: storage pays the debt; exhausted storage => starving + deficit. --!
   subroutine test_negative_npp()
      type(carbon_gain_t)   :: g1, g2
      type(carbon_loss_t)   :: loss
      type(carbon_demand_t) :: dem
      type(carbon_npp_t)    :: n1, n2
      loss = carbon_loss_t(leaf=0.05_wp, fineroot=0.05_wp)
      dem  = carbon_demand_t(leaf=0.0_wp, fineroot=0.0_wp, storage=0.0_wp, wood=1.0e6_wp)

      !----- Storage covers the debt. ----------------------------------------------------!
      g1 = carbon_gain_t(net_carbon=-0.4_wp, storage=2.0_wp)
      call plant_carbon_allocation(g1, loss, dem, PHEN_ON, n1)
      call check_true ('neg NPP (covered): not starving',      .not. n1%starving)
      call check_close('neg NPP (covered): deficit = 0',        n1%deficit, 0.0_wp)
      call check_close('neg NPP (covered): storage pays debt',  n1%nonstructural, -0.4_wp)
      call check_close('neg NPP (covered): no wood growth',     n1%wood, 0.0_wp)
      call check_closure('neg NPP (covered): carbon closes', g1, loss, n1)

      !----- Storage too small: starving, unpaid respiration reported as deficit. ---------!
      g2 = carbon_gain_t(net_carbon=-1.0_wp, storage=0.3_wp)
      call plant_carbon_allocation(g2, loss, dem, PHEN_ON, n2)
      call check_true ('neg NPP (starving): starving flag set', n2%starving)
      call check_close('neg NPP (starving): deficit = debt - storage', n2%deficit, 0.7_wp)
      call check_close('neg NPP (starving): storage drained',   n2%nonstructural, -0.3_wp)
      call check_closure('neg NPP (starving): carbon closes', g2, loss, n2)
   end subroutine test_negative_npp

end program test_plant_carbon_dynamics
