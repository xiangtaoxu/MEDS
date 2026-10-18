! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_trait_dynamics -- unit tests for the light-plasticity kernel                   !
! (docs/science/plant_traits.md eqs. (1)-(2)).                                             !
!                                                                                          !
!   1. NO SHADE      : cum_lai_above = 0 => every target equals its top-of-canopy value.       !
!   2. SHADE SIGNS   : shaded => SLA up, Vcmax/Rd and leaf turnover down (ED2 directions).    !
!   3. GRADIENT      : target = trait_toc * exp(kplastic * cum_lai_above); the turnover rate's  !
!                      target is the reciprocal of a lifespan that lengthens as exp(-k L).      !
!   4. RELAX RATE    : update_plastic_trait moves by f = 1 - exp(-rate dt); dt->0 barely moves, !
!                      a slow-turnover (shaded) canopy acclimates slowly, dt >> 1/rate reaches  !
!                      target, and the instantaneous (census) mode reaches it at once.            !
!==========================================================================================!
program test_plant_trait_dynamics
   use meds_test_assert, only : check_close, check_true, test_report
   use meds_kinds,                only : wp, ik
   use meds_plant_trait_dynamics, only : light_plastic_traits, update_plastic_trait
   implicit none

   real(wp)    :: sla_t, vc_t, rd_t, turnover_t

   !----- 1. No shade: targets == top-of-canopy. ------------------------------------------!
   call light_plastic_traits(0.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 0.5_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, -0.05_wp, sla_t, vc_t, rd_t, turnover_t)
   call check_close('no shade: sla = toc', sla_t, 15.0_wp)
   call check_close('no shade: vcmax = toc', vc_t, 50.0_wp)
   call check_close('no shade: rd = toc', rd_t, 1.0_wp)
   call check_close('no shade: leaf turnover rate = toc', turnover_t, 0.5_wp)

   !----- 2 + 3. Shaded (cum_lai_above = 3): signs + exact exponential gradient. -----------!
   call light_plastic_traits(3.0_wp, 15.0_wp, 50.0_wp, 1.0_wp, 0.5_wp,                          &
                             0.1_wp, -0.1_wp, -0.1_wp, -0.05_wp, sla_t, vc_t, rd_t, turnover_t)
   call check_true('shade: SLA increases', sla_t > 15.0_wp)
   call check_true('shade: Vcmax decreases', vc_t  < 50.0_wp)
   call check_true('shade: Rd decreases', rd_t  < 1.0_wp)
   call check_true('shade: leaf turnover slows', turnover_t < 0.5_wp)
   call check_close('shade: SLA gradient', sla_t, 15.0_wp * exp( 0.1_wp * 3.0_wp))
   call check_close('shade: Vcmax gradient', vc_t, 50.0_wp * exp(-0.1_wp * 3.0_wp))
   call check_close('shade: turnover gradient', turnover_t, 0.5_wp * exp(-0.05_wp * 3.0_wp))
   !----- A lifespan of 2 yr lengthening as exp(0.05 L) in shade is this rate, exactly. ---!
   call check_close('shade: turnover = 1 / the lifespan form', turnover_t, 1.0_wp / (2.0_wp * exp(0.05_wp * 3.0_wp)))

   !----- 4. Replacement-weighted update (gradual) + the instantaneous (census) mode. ------!
   call check_close('update: rate dt = 1 => f = 1-exp(-1)',                                                           &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.5_wp, 2.0_wp, .false.),                                           &
           15.0_wp + (1.0_wp - exp(-1.0_wp)) * 5.0_wp)
   call check_true('update: tiny dt barely moves',                                                                    &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 0.5_wp, 1.0e-3_wp, .false.) - 15.0_wp) < 0.01_wp)
   call check_true('update: slow turnover => slow',                                                                   &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.01_wp, 1.0_wp, .false.) < 15.1_wp)
   call check_true('update: dt >> 1/rate => target',                                                                  &
           abs(update_plastic_trait(15.0_wp, 20.0_wp, 10.0_wp, 10.0_wp, .false.) - 20.0_wp) < 1.0e-6_wp)
   call check_close('update: instant => target regardless of rate dt',                                                &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.01_wp, 1.0e-6_wp, .true.), 20.0_wp)
   call check_close('update: a rate of 0 replaces no leaves',                                                         &
           update_plastic_trait(15.0_wp, 20.0_wp, 0.0_wp, 1.0_wp, .false.), 15.0_wp)

   call test_report('test_plant_trait_dynamics')

contains



end program test_plant_trait_dynamics
