!==========================================================================================!
! test_plant_carbon_allocation -- unit tests for the elemental, GROWTH-ONLY daily carbon      !
! allocation kernel + growth_respiration + the turnover baseline rates. (Tissue shed/turnover  !
! is applied by the driver's update_biomass_turnover -- tested in test_carbon_growth.)         !
!                                                                                          !
!   1. GROWTH RESP     : growth_respiration = g * max(0, npp_growth).                            !
!   2. CLOSURE         : (growth pools + npp_store) - deficit = (gpp - resp) - growth_resp.       !
!   3. GROWTH-RESP CHG : charged (1+g) on REALIZED growth only; storage refill is exempt.         !
!   4. WOOD RESIDUAL   : a budget that only covers the leaf demand => wood 0.                     !
!   5. STORAGE GROWTH  : storage funds leaf growth EVEN when net < 0 (spring leaf-out).            !
!   6. STARVING        : maintenance debt beyond storage => starving + deficit, no growth.         !
!   7. TURNOVER FLOOR  : baseline turnover as a [1/day] shed rate; max(active, baseline).          !
!==========================================================================================!
program test_plant_carbon_allocation
   use meds_kinds,           only : wp, ik
   use meds_plant_interface, only : plant_carbon_allocation, growth_respiration,                 &
                                    pheno_drives_to_rates, turnover_shed_rates
   implicit none

   real(wp), parameter :: YR_DAY = 365.2425_wp
   integer(ik) :: nfail
   nfail = 0_ik

   call test_growth_respiration()
   call test_closure()
   call test_growth_resp_charge()
   call test_wood_residual()
   call test_storage_growth()
   call test_starving()
   call test_turnover_floor()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_carbon_allocation: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_carbon_allocation: ', nfail, ' FAILED'
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

   !----- The growth-side carbon-closure identity the kernel guarantees on every call. -----!
   subroutine check_closure(name, gpp, resp, gl, gf, gw, gr, gs, gresp, def)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: gpp, resp, gl, gf, gw, gr, gs, gresp, def
      call check_close(name, (gl + gf + gw + gr + gs) - def, (gpp - resp) - gresp)
   end subroutine check_closure

   !----- 1. growth_respiration = g * max(0, npp_growth) (relocated from the respiration module). !
   subroutine test_growth_respiration()
      call check_close('growth_resp: 0.3 * 10 = 3', growth_respiration(10.0_wp, 0.3_wp), 3.0_wp)
      call check_close('growth_resp: npp<=0 => 0',  growth_respiration(-5.0_wp, 0.3_wp), 0.0_wp)
   end subroutine test_growth_respiration

   !----- 2. Carbon closes (g = 0 for the simplest arithmetic): surplus -> wood. ------------!
   subroutine test_closure()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def
      logical  :: starv
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.0_wp, 5.0_wp, 0.3_wp, 0.2_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, starv)
      call check_close('closure: growth_leaf = demand', gl, 0.3_wp)
      call check_close('closure: growth_wood = residual', gw, 0.5_wp)
      call check_true ('closure: not starving', .not. starv)
      call check_closure('closure: g=0', 1.0_wp, 0.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_closure

   !----- 3. Growth respiration charged (1+g) on realized growth; storage refill exempt. ----!
   subroutine test_growth_resp_charge()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def
      logical  :: starv
      ! g = 0.5: leaf demand 0.3 built at cost 1.5 (gresp += 0.15); storage 0.2 at 1:1 (exempt);
      ! residual 0.35 -> wood 0.35/1.5.
      call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.2_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, starv)
      call check_close('growth resp: on realized growth', gresp, 0.15_wp + 0.5_wp*(0.35_wp/1.5_wp))
      call check_close('growth resp: storage exempt (npp_store = 0.2)', gs, 0.2_wp)
      call check_close('growth resp: wood tissue = residual/(1+g)', gw, 0.35_wp/1.5_wp)
      call check_closure('growth resp: closes', 1.0_wp, 0.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_growth_resp_charge

   !----- 4. Wood is the residual sink: a budget that only covers the leaf demand => wood 0. !
   subroutine test_wood_residual()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def
      logical  :: starv
      call plant_carbon_allocation(0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, starv)
      call check_close('residual: growth_leaf = demand', gl, 0.5_wp)
      call check_close('residual: growth_wood = 0',      gw, 0.0_wp)
   end subroutine test_wood_residual

   !----- 5. Storage funds leaf growth even when net < 0 (spring leaf-out from reserves). ---!
   subroutine test_storage_growth()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def
      logical  :: starv
      ! net = -0.05 (paid from storage), then 0.1 leaf built from the remaining reserves (g=0).
      call plant_carbon_allocation(0.0_wp, 0.05_wp, 0.0_wp, 1.0_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, starv)
      call check_true ('storage growth: leaf grows despite net<0', gl > 0.0_wp)
      call check_close('storage growth: growth_leaf = demand',     gl, 0.1_wp)
      call check_close('storage growth: storage draw = debt+growth', gs, -0.15_wp)
      call check_true ('storage growth: not starving',             .not. starv)
      call check_closure('storage growth: closes', 0.0_wp, 0.05_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_storage_growth

   !----- 6. Maintenance debt beyond storage => starving + deficit; no growth. --------------!
   subroutine test_starving()
      real(wp) :: gl, gf, gw, gs, gr, gresp, def
      logical  :: starv
      call plant_carbon_allocation(0.0_wp, 1.0_wp, 0.0_wp, 0.3_wp, 0.1_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
           gl, gf, gw, gs, gr, gresp, def, starv)
      call check_true ('starving: flag set',           starv)
      call check_close('starving: deficit = debt - storage', def, 0.7_wp)
      call check_close('starving: no leaf growth',      gl, 0.0_wp)
      call check_close('starving: storage drained',     gs, -0.3_wp)
      call check_closure('starving: closes', 0.0_wp, 1.0_wp, gl, gf, gw, gr, gs, gresp, def)
   end subroutine test_starving

   !----- 7. Turnover as a baseline [1/day] shed rate, and the max(active, baseline) floor. --!
   subroutine test_turnover_floor()
      real(wp) :: lb, rb, flush, shed, root_shed
      call turnover_shed_rates(0.5_wp, 1.2_wp, .false., 278.15_wp, 0.4_wp, 300.0_wp, lb, rb)
      call check_close('turnover: deciduous leaf base', lb, 0.5_wp / YR_DAY)
      call check_close('turnover: deciduous root base', rb, 1.2_wp / YR_DAY)
      call turnover_shed_rates(0.5_wp, 1.2_wp, .true., 278.15_wp, 0.4_wp, 278.15_wp, lb, rb)
      call check_close('turnover: evergreen @ref => x0.5', lb, 0.5_wp * (0.5_wp / YR_DAY))
      call pheno_drives_to_rates(1.0_wp, 0.0_wp, 0.0667_wp, 0.05_wp, 0.5_wp, 1.2_wp, .false.,     &
           278.15_wp, 0.4_wp, 300.0_wp, flush, shed, root_shed)
      call check_close('drives: flush = k_flush*drive',   flush, 0.0667_wp)
      call check_close('drives: shed floored to baseline', shed, 0.5_wp / YR_DAY)
      call pheno_drives_to_rates(1.0_wp, 1.0_wp, 0.0667_wp, 0.05_wp, 0.5_wp, 1.2_wp, .false.,     &
           278.15_wp, 0.4_wp, 300.0_wp, flush, shed, root_shed)
      call check_close('drives: active shed wins when large', shed, 0.05_wp)
   end subroutine test_turnover_floor

end program test_plant_carbon_allocation
