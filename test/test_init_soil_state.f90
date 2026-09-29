! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_init_soil_state -- [init].soil_temp and soil_theta reach every soil layer of every patch. !
!                                                                                          !
! The initial soil state was two constants in the fast context. A run without a spin-up starts   !
! from it, so the keys matter only if they reach the column: this asserts the path from the      !
! config through build_fast_context to init_fast_reservoirs, and pins the defaults to the old     !
! constants so that no existing config changes.                                                 !
!==========================================================================================!
program test_init_soil_state
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_init,             only : init_bare_ground, finalize_init
   use meds_fast_dynamics,    only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_test_support,     only : banner, build_test_config, check, check_close
   implicit none
   type(meds_config_t)  :: cfg
   type(site_t)         :: site
   type(fast_context_t) :: ctx
   integer(ik)          :: ip, nsl

   call banner('[init].soil_temp and soil_theta reach the soil column')
   cfg = build_test_config()

   !=== 1. The defaults are the constants the fast context used to carry. ====================!
   call check_close(cfg%init_soil_temp,  288.0_wp, 1.0e-14_wp, 'default soil_temp changed (was 288 K)')
   call check_close(cfg%init_soil_theta, 0.30_wp,  1.0e-14_wp, 'default soil_theta changed (was 0.30)')

   !=== 2. Set values reach the context, then every layer of every patch. ===================!
   cfg%init_soil_temp  = 298.65_wp
   cfg%init_soil_theta = 0.25_wp
   call build_fast_context(cfg, ctx)
   call check_close(ctx%soil_temp_init, 298.65_wp, 1.0e-14_wp, 'soil_temp did not reach the fast context')
   call check_close(ctx%theta_init,     0.25_wp,   1.0e-14_wp, 'soil_theta did not reach the fast context')
   call init_bare_ground(site, cfg, 3_ik)
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)
   nsl = ctx%col_config%soil%n_active
   do ip = 1_ik, site%patch%n
      call check(all(abs(site%patch%soil_e(ip)%soil_temp(1:nsl) - 298.65_wp) < 1.0e-12_wp),         &
                 'a soil layer did not start at [init].soil_temp')
      call check(all(abs(site%patch%soil_w(ip)%theta(1:nsl) - 0.25_wp) < 1.0e-12_wp),              &
                 'a soil layer did not start at [init].soil_theta')
   end do
   write(*,'(a)') '   PASS'
end program test_init_soil_state
