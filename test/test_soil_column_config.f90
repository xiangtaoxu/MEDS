! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_soil_column_config -- the [soil] column keys reach the column they describe.         !
!                                                                                          !
! The soil column's grid, texture and thermal properties are config keys, worth anything   !
! only if they actually change the column, and that wiring is the kind that breaks         !
! silently: a mis-threaded argument leaves the model running on a default, every budget    !
! still closes, and the run looks healthy. So this asserts the knobs turn, not merely that !
! the code compiles.                                                                       !
!                                                                                          !
! It also pins the default column: SOIL_LAYER_TEMPLATE (SoilGrids' depths, then deeper) cut at  !
! the 2 m default depth, and the cut and merge rule of soil_layer_grid                          !
! (docs/science/soil_biophysics.md §The soil layers).                                      !
!==========================================================================================!
program test_soil_column_config
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : n_soil_layer_max
   use meds_soil_column_config, only : soil_column_config_t, soil_column_layers, soil_params_t, soil_thermal_params_t, &
                                       build_soil_hydr_params, build_soil_therm_params
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use meds_test_support, only : banner, check, check_close
   implicit none

   type(soil_column_config_t)  :: sc
   type(soil_params_t)         :: p, p2
   type(soil_thermal_params_t) :: th
   real(wp)                    :: bottoms(n_soil_layer_max)
   integer(ik)                 :: n_layer, k
   real(wp), parameter         :: SOILGRIDS_BOTTOMS(5) = [0.05_wp, 0.15_wp, 0.30_wp, 0.60_wp, 1.00_wp]

   call banner('[soil]: the column keys reach the column')

   !=== 1. The defaults: a 2 m column on the default layers. =================================!
   call check(sc%n_layer_bottoms == 0_ik,                'default layers are not the template')
   call check_close(sc%depth,       2.0_wp,   1.0e-14_wp, 'default depth changed (was the hard-coded 2.0 m)')
   call check_close(sc%theta_sat,   0.43_wp,  1.0e-14_wp, 'default theta_sat changed')
   call check_close(sc%ksat,        2.89e-6_wp, 1.0e-14_wp, 'default ksat changed')
   call check(sc%retention == SOIL_RETENTION_VG,         'default retention is not van Genuchten')

   !=== 2. The default column is the template cut at 2 m: seven layers, a 5 cm top layer,    !
   !       SoilGrids' bottoms to 1 m, and the base at -depth exactly. ========================!
   call build_column(sc, p)
   call check(p%n_active == 7_ik,                        'the default 2 m column is not 7 layers')
   call check_close(p%dz(1), 0.05_wp, 1.0e-12_wp,        'the top layer is not 5 cm')
   do k = 1_ik, 5_ik
      call check_close(-p%soil_layer_z(k + 1_ik), SOILGRIDS_BOTTOMS(k), 1.0e-12_wp,                    &
                       'a layer bottom is off SoilGrids'' depths')
   end do
   call check_close(p%soil_layer_z(p%n_active + 1_ik), -sc%depth, 1.0e-12_wp, 'column base is not at -depth')

   !=== 3. depth cuts the same list: 0.3 m is 5, 10 and 15 cm; 4 m is nine layers whose top  !
   !       is still 5 cm, the last running from 3 to 4 m. ===================================!
   sc%depth = 0.3_wp
   call build_column(sc, p2)
   call check(p2%n_active == 3_ik .and. abs(p2%dz(3) - 0.15_wp) < 1.0e-12_wp,                       &
              'a 0.3 m column is not 5, 10 and 15 cm')
   sc%depth = 4.0_wp
   call build_column(sc, p2)
   call check(p2%n_active == 9_ik,                       'a 4 m column is not 9 layers')
   call check_close(p2%dz(1), 0.05_wp, 1.0e-12_wp,       'a deeper column changed the top layer')
   call check_close(p2%dz(9), 1.0_wp,  1.0e-12_wp,       'the 4 m column does not end in a 3-4 m layer')

   !=== 4. A remainder under half of the list layer it falls in joins the layer above: 0.65 m !
   !       ends in a 30-65 cm layer, while 0.85 m keeps a 60-85 cm layer of its own. ==========!
   sc%depth = 0.65_wp
   call soil_column_layers(sc, bottoms, n_layer)
   call check(n_layer == 4_ik .and. abs(bottoms(4) - 0.65_wp) < 1.0e-12_wp .and.                    &
              abs(bottoms(3) - 0.30_wp) < 1.0e-12_wp, 'a 5 cm remainder was not merged upward')
   sc%depth = 0.85_wp
   call soil_column_layers(sc, bottoms, n_layer)
   call check(n_layer == 5_ik .and. abs(bottoms(5) - 0.85_wp) < 1.0e-12_wp,                         &
              'a 25 cm remainder was merged although it is over half its layer')

   !=== 5. layer_bottoms replaces the list, and depth still cuts it. ========================!
   sc%n_layer_bottoms = 3_ik ; sc%layer_bottoms(1:3) = [0.2_wp, 0.5_wp, 1.0_wp] ; sc%depth = 1.0_wp
   call build_column(sc, p2)
   call check(p2%n_active == 3_ik .and. abs(p2%dz(1) - 0.2_wp) < 1.0e-12_wp,                         &
              'layer_bottoms did not reach the column')
   sc%depth = 0.8_wp
   call build_column(sc, p2)
   call check(p2%n_active == 3_ik .and. abs(p2%soil_layer_z(4) + 0.8_wp) < 1.0e-12_wp,              &
              'depth did not cut the layer_bottoms list')
   sc%n_layer_bottoms = 0_ik ; sc%layer_bottoms = 0.0_wp ; sc%depth = 2.0_wp

   !=== 6. Texture reaches the layers, and the DERIVED thresholds follow it. ================!
   sc%theta_sat = 0.55_wp ; sc%theta_res = 0.02_wp
   call build_column(sc, p2)
   call check_close(p2%theta_sat(1), 0.55_wp, 1.0e-14_wp, 'theta_sat did not reach the layers')
   call check_close(p2%theta_res(1), 0.02_wp, 1.0e-14_wp, 'theta_res did not reach the layers')
   call check(p2%theta_fc(1) > p2%theta_res(1),          'field capacity must exceed residual content')
   call check(p2%theta_fc(1) < p2%theta_sat(1),          'field capacity must be below saturation')
   sc%theta_sat = 0.43_wp ; sc%theta_res = 0.078_wp

   !=== 7. The retention family selects which pair curve_par_a/n mean: for Campbell, the air-entry !
   !       suction psi_sat [m] (< 0) and the exponent b, as validate_config requires. ================!
   sc%retention = SOIL_RETENTION_CAMPBELL
   sc%curve_par_a = -0.26_wp ; sc%curve_par_n = 5.65_wp
   call build_column(sc, p2)
   call check(p2%retention == SOIL_RETENTION_CAMPBELL,   'retention family did not reach the column')
   call check(p2%theta_fc(1) > p2%theta_res(1) .and. p2%theta_fc(1) < p2%theta_sat(1),           &
              'Campbell field capacity lies between residual content and saturation')
   sc%retention = SOIL_RETENTION_VG
   sc%curve_par_a = 3.6_wp ; sc%curve_par_n = 1.56_wp                  ! the van Genuchten defaults

   !=== 8. The thermal keys reach the thermal column, on the same layers. ====================!
   sc%solid_conductivity = 2.5_wp ; sc%dry_conductivity = 0.20_wp ; sc%dry_heat_capacity = 1.5e6_wp
   call soil_column_layers(sc, bottoms, n_layer)
   call build_soil_therm_params(n_layer, sc%solid_conductivity, sc%dry_conductivity,                &
                                sc%dry_heat_capacity, th)
   call check(th%nzg_active == n_layer,                  'the layer count did not reach the thermal column')
   call check_close(th%soil_solid_conductivity(1), 2.5_wp,  1.0e-14_wp, 'solid_conductivity did not reach it')
   call check_close(th%soil_dry_conductivity(1),   0.20_wp, 1.0e-14_wp, 'dry_conductivity did not reach it')
   call check_close(th%soil_dry_heat_capacity(1),  1.5e6_wp, 1.0e-8_wp, 'dry_heat_capacity did not reach it')

   print '(a)', 'test_soil_column_config: all checks passed'

contains

   !----- Exactly the call build_fast_context makes, so the test cannot drift from it. -------!
   subroutine build_column(c, out)
      type(soil_column_config_t), intent(in)  :: c
      type(soil_params_t),        intent(out) :: out
      real(wp)    :: column_bottoms(n_soil_layer_max)
      integer(ik) :: n_column
      call soil_column_layers(c, column_bottoms, n_column)
      call build_soil_hydr_params(n_column, c%retention, column_bottoms, c%theta_sat, c%theta_res, &
                                  c%ksat, c%curve_par_a, c%curve_par_n, c%psi_fc, out)
   end subroutine build_column

end program test_soil_column_config
