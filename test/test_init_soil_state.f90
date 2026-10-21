! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_init_soil_state -- [init].soil_temp and soil_theta reach every soil layer of every patch. !
!                                                                                          !
! A run without a spin-up starts from the initial soil state, so the keys matter only if   !
! they reach the column: this asserts the path from the config through build_fast_context  !
! to init_fast_reservoirs, and pins the defaults (288 K, theta = 0.30). A column with a    !
! texture per layer (the soil maps) builds each layer's own curve and starts each layer at !
! its own field capacity. A bare-ground start is empty unless a PFT's nbg_init_density     !
! seeds it, at min_cohort_height.                                                          !
!==========================================================================================!
program test_init_soil_state
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t, soil_column_layers
   use meds_polygon_state_vars, only : polygon_t
   use meds_init_bare_ground, only : init_bare_ground, seed_bare_ground, finalize_init
   use meds_constants,        only : n_soil_layer_max
   use meds_fast_dynamics,    only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_test_support,     only : banner, build_test_config, check, check_close
   implicit none
   type(meds_config_t)  :: cfg
   type(polygon_t)      :: poly
   type(fast_context_t) :: ctx
   integer(ik)          :: ip, nsl, k, n
   real(wp)             :: bottoms(n_soil_layer_max)

   call banner('[init].soil_temp and soil_theta reach the soil column')
   cfg = build_test_config()

   !=== 1. The defaults: 288 K and theta = 0.30. ===========================================!
   call check_close(cfg%init%soil_temp,  288.0_wp, 1.0e-14_wp, 'default soil_temp changed (was 288 K)')
   call check_close(cfg%init%soil_theta, 0.30_wp,  1.0e-14_wp, 'default soil_theta changed (was 0.30)')

   !=== 2. Set values reach the context, then every layer of every patch. ===================!
   cfg%init%soil_temp  = 298.65_wp
   cfg%init%soil_theta = 0.25_wp
   call build_fast_context(cfg, ctx)
   call check_close(ctx%soil_temp_init, 298.65_wp, 1.0e-14_wp, 'soil_temp did not reach the fast context')
   call init_bare_ground(poly, cfg, 3_ik)
   call check_close(poly%site%surface(1)%theta_start(1), 0.25_wp, 1.0e-14_wp, 'soil_theta did not reach the site')
   call finalize_init(poly)
   call init_fast_reservoirs(poly, ctx)
   nsl = ctx%col_config%soil%n_active
   do ip = 1_ik, poly%patch%n
      call check(all(abs(poly%patch%soil_energy(ip)%soil_temp(1:nsl) - 298.65_wp) < 1.0e-12_wp),         &
                 'a soil layer did not start at [init].soil_temp')
      call check(all(abs(poly%patch%soil_water(ip)%theta(1:nsl) - 0.25_wp) < 1.0e-12_wp),              &
                 'a soil layer did not start at [init].soil_theta')
   end do

   !=== 3. A texture per layer: each layer's own curve, each starting at its field capacity. ======!
   call soil_column_layers(cfg%soil_column, bottoms, n)
   associate (sc => cfg%soil_column)
      sc%n_layer_values = n
      do k = 1_ik, n
         sc%layer_theta_sat(k)   = 0.40_wp + 0.01_wp * real(k, wp)
         sc%layer_theta_res(k)   = 0.05_wp
         sc%layer_ksat(k)        = 1.0e-6_wp / real(k, wp)
         sc%layer_curve_par_a(k) = 2.0_wp + 0.1_wp * real(k, wp)
         sc%layer_curve_par_n(k) = 1.3_wp + 0.02_wp * real(k, wp)
         sc%layer_solid_conductivity(k) = 3.0_wp ; sc%layer_dry_conductivity(k) = 0.2_wp
         sc%layer_dry_heat_capacity(k)  = 1.5e6_wp + 1.0e4_wp * real(k, wp)
      end do
   end associate
   call build_fast_context(cfg, ctx)
   associate (soil => ctx%col_config%soil)
      call check(soil%n_active == n .and. all(soil%theta_sat(1:n) == cfg%soil_column%layer_theta_sat(1:n))  &
                 .and. all(soil%curve_par_n(1:n) == cfg%soil_column%layer_curve_par_n(1:n)),                       &
                 'a layer did not take its own curve')
      call check(all(ctx%col_config%soil_thermal%soil_dry_heat_capacity(1:n) ==                             &
                     cfg%soil_column%layer_dry_heat_capacity(1:n)), 'a layer did not take its own heat capacity')
   end associate
   call init_bare_ground(poly, cfg, 3_ik)
   associate (surf => poly%site%surface(1))
      call check(all(surf%theta_start(1:n) == surf%soil%theta_fc(1:n)) .and.                                 &
                 surf%soil%theta_fc(1) /= surf%soil%theta_fc(n), 'a layer did not start at its own field capacity')
   end associate
   call finalize_init(poly)
   call init_fast_reservoirs(poly, ctx)
   call check(all(poly%patch%soil_water(2)%theta(1:n) == ctx%col_config%soil%theta_fc(1:n)),                    &
              'a patch''s soil did not start at field capacity')

   !=== 4. Bare ground is empty by default; nbg_init_density seeds each patch with each PFT. ======!
   call init_bare_ground(poly, cfg, 3_ik)
   call seed_bare_ground(poly, cfg)
   call check(poly%cohort%n == 0_ik, 'bare ground is not empty with nbg_init_density = 0')
   cfg%pft%nbg_init_density = [0.5_wp, 0.0_wp, 0.2_wp]
   call init_bare_ground(poly, cfg, 3_ik)
   call seed_bare_ground(poly, cfg)
   call check(poly%cohort%n == 6_ik, 'not one cohort per seeded PFT per patch')
   call check(all(abs(poly%cohort%height(1:6) - cfg%pft%min_cohort_height) < 1.0e-9_wp),                    &
              'a seeded cohort is not at min_cohort_height')
   call check(all(poly%cohort%nplant(1:6) == 0.5_wp .or. poly%cohort%nplant(1:6) == 0.2_wp) .and.             &
              count(poly%cohort%pft(1:6) == 1_ik) == 3 .and. count(poly%cohort%pft(1:6) == 3_ik) == 3,       &
              'a seeded cohort has the wrong PFT or density')
   write(*,'(a)') '   PASS'
end program test_init_soil_state
