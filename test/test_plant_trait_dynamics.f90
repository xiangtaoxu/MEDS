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
   use meds_test_assert, only : check_close, check_true, test_report
   use meds_kinds,                only : wp, ik
   use meds_plant_trait_dynamics, only : light_plastic_traits, update_plastic_trait
   implicit none

   real(wp)    :: sla_t, vc_t, rd_t, ll_t

   !----- 1. No shade: targets == top-of-canopy. ------------------------------------------!
   call light_plastic_traits(0.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 2.0_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, 0.05_wp, sla_t, vc_t, rd_t, ll_t)
   call check_close('no shade: sla = toc', sla_t, 15.0_wp)
   call check_close('no shade: vcmax = toc', vc_t, 50.0_wp)
   call check_close('no shade: rd = toc', rd_t, 1.0_wp)
   call check_close('no shade: llspan = toc', ll_t, 2.0_wp)

   !----- 2 + 3. Shaded (cum_lai_above = 3): signs + exact exponential gradient. -----------!
   call light_plastic_traits(3.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 2.0_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, 0.05_wp, sla_t, vc_t, rd_t, ll_t)
   call check_true('shade: SLA increases', sla_t > 15.0_wp)
   call check_true('shade: Vcmax decreases', vc_t  < 50.0_wp)
   call check_true('shade: Rd decreases', rd_t  < 1.0_wp)
   call check_true('shade: leaf lifespan increases', ll_t > 2.0_wp)
   call check_close('shade: SLA gradient', sla_t, 15.0_wp * exp( 0.1_wp * 3.0_wp))
   call check_close('shade: Vcmax gradient', vc_t, 50.0_wp * exp(-0.1_wp * 3.0_wp))
   call check_close('shade: llspan gradient', ll_t, 2.0_wp  * exp( 0.05_wp * 3.0_wp))

   !----- 4. Replacement-weighted update (gradual) + the instantaneous (census) mode. ------!
   call check_close('update: dt = llspan => f = 1-exp(-1)',                                                           &
           update_plastic_trait(15.0_wp, 20.0_wp, 2.0_wp, 2.0_wp, .false.),                                           &
           15.0_wp + (1.0_wp - exp(-1.0_wp)) * 5.0_wp)
   call check_true('update: tiny dt barely moves',                                                                    &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 2.0_wp, 1.0e-3_wp, .false.) - 15.0_wp) < 0.01_wp)
   call check_true('update: long-lived => slow',                                                                      &
           update_plastic_trait(15.0_wp, 20.0_wp, 100.0_wp, 1.0_wp, .false.) < 15.1_wp)
   call check_true('update: dt >> llspan => target',                                                                  &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 0.1_wp, 10.0_wp, .false.) - 20.0_wp) < 1.0e-6_wp)
   call check_close('update: instant => target regardless of dt/llspan',                                              &
           update_plastic_trait(15.0_wp, 20.0_wp, 100.0_wp, 1.0e-6_wp, .true.), 20.0_wp)

   call test_report('test_plant_trait_dynamics')

contains



end program test_plant_trait_dynamics
