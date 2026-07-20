!==========================================================================================!
! test_plant_trait_dynamics -- unit tests for the light-plasticity kernel.                  !
!                                                                                          !
!   1. NO SHADE      : cum_lai_above = 0 => every target equals its top-of-canopy value.       !
!   2. SHADE SIGNS   : shaded => SLA up, Vcmax/Rd down, leaf lifespan up (ED2 directions).      !
!   3. GRADIENT      : target = trait_toc * exp(kplastic * cum_lai_above).                      !
!   4. RELAX RATE    : relax_trait moves by f = 1 - exp(-dt/llspan); dt->0 barely moves,        !
!                      a long-lived (shaded) canopy acclimates slowly, dt>>llspan reaches target.!
!==========================================================================================!
program test_plant_trait_dynamics
   use meds_kinds,                only : wp, ik
   use meds_plant_trait_dynamics, only : light_plastic_traits, relax_trait
   implicit none

   integer(ik) :: nfail
   real(wp)    :: sla_t, vc_t, rd_t, ll_t
   nfail = 0_ik

   !----- 1. No shade: targets == top-of-canopy. ------------------------------------------!
   call light_plastic_traits(0.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 2.0_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, 0.05_wp, sla_t, vc_t, rd_t, ll_t)
   call cc('no shade: sla = toc',    sla_t, 15.0_wp)
   call cc('no shade: vcmax = toc',  vc_t,  50.0_wp)
   call cc('no shade: rd = toc',     rd_t,  1.0_wp)
   call cc('no shade: llspan = toc', ll_t,  2.0_wp)

   !----- 2 + 3. Shaded (cum_lai_above = 3): signs + exact exponential gradient. -----------!
   call light_plastic_traits(3.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 2.0_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, 0.05_wp, sla_t, vc_t, rd_t, ll_t)
   call ct('shade: SLA increases',        sla_t > 15.0_wp)
   call ct('shade: Vcmax decreases',      vc_t  < 50.0_wp)
   call ct('shade: Rd decreases',         rd_t  < 1.0_wp)
   call ct('shade: leaf lifespan increases', ll_t > 2.0_wp)
   call cc('shade: SLA gradient',   sla_t, 15.0_wp * exp( 0.1_wp * 3.0_wp))
   call cc('shade: Vcmax gradient', vc_t,  50.0_wp * exp(-0.1_wp * 3.0_wp))
   call cc('shade: llspan gradient', ll_t, 2.0_wp  * exp( 0.05_wp * 3.0_wp))

   !----- 4. Replacement-weighted relaxation. ---------------------------------------------!
   call cc('relax: dt = llspan => f = 1-exp(-1)', relax_trait(15.0_wp, 20.0_wp, 2.0_wp, 2.0_wp), &
           15.0_wp + (1.0_wp - exp(-1.0_wp)) * 5.0_wp)
   call ct('relax: tiny dt barely moves',   abs(relax_trait(15.0_wp, 20.0_wp, 2.0_wp, 1.0e-3_wp) - 15.0_wp) < 0.01_wp)
   call ct('relax: long-lived => slow',     relax_trait(15.0_wp, 20.0_wp, 100.0_wp, 1.0_wp) < 15.1_wp)
   call ct('relax: dt >> llspan => target', abs(relax_trait(15.0_wp, 20.0_wp, 0.1_wp, 10.0_wp) - 20.0_wp) < 1.0e-6_wp)

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_trait_dynamics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_trait_dynamics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine ct(name, cond)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a)', '  FAIL : ', name
      end if
   end subroutine ct

   subroutine cc(name, got, expect)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect
      if (abs(got - expect) <= 1.0e-9_wp * max(1.0_wp, abs(expect))) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a,es13.6,a,es13.6)', '  FAIL : ', name, got, ' expected ', expect
      end if
   end subroutine cc

end program test_plant_trait_dynamics
