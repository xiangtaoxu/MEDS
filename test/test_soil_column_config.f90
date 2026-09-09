!==========================================================================================!
! test_soil_column_config -- the [soil_column] keys reach the column they describe.         !
!                                                                                          !
! The soil column's grid, texture and thermal properties used to be literals inside          !
! build_fast_context. Moving them to config is only worth anything if the keys actually       !
! change the column, and that wiring is the kind that breaks silently: a mis-threaded          !
! argument leaves the model running the old placeholder, every budget still closes, and the    !
! run looks healthy. So this asserts the knobs turn, not merely that the code compiles.        !
!                                                                                          !
! It also pins the DEFAULTS to the literals they replaced, because the whole migration rests   !
! on an absent block reproducing the old behaviour.                                            !
!==========================================================================================!
program test_soil_column_config
   use meds_kinds,         only : wp, ik
   use meds_config,        only : soil_column_config_t
   use meds_column_params, only : soil_params_t, soil_thermal_params_t,                            &
                                  build_soil_hydr_params, build_soil_therm_params
   use meds_hydr_lib,      only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use meds_test_support,  only : check, check_close, banner
   implicit none

   type(soil_column_config_t)  :: sc
   type(soil_params_t)         :: p, p2
   type(soil_thermal_params_t) :: th

   call banner('[soil_column]: the keys reach the column')

   !=== 1. The DEFAULTS are the literals build_fast_context used to carry. ==================!
   call check(sc%n_layer == 10_ik,                       'default n_layer changed (was the hard-coded 10)')
   call check_close(sc%depth,       2.0_wp,   1.0e-14_wp, 'default depth changed (was the hard-coded 2.0 m)')
   call check_close(sc%grid_growth, 3.0_wp,   1.0e-14_wp, 'default grid_growth changed')
   call check_close(sc%theta_sat,   0.43_wp,  1.0e-14_wp, 'default theta_sat changed')
   call check_close(sc%ksat,        2.89e-6_wp, 1.0e-14_wp, 'default ksat changed')
   call check(sc%retention == SOIL_RETENTION_VG,         'default retention is not van Genuchten')

   !=== 2. depth reaches the grid: the deepest interface is -depth, exactly. ================!
   call build_column(sc, p)
   call check_close(p%soil_layer_z(sc%n_layer + 1_ik), -sc%depth, 1.0e-12_wp,                      &
                    'column base is not at -depth')
   sc%depth = 4.0_wp
   call build_column(sc, p2)
   call check_close(p2%soil_layer_z(sc%n_layer + 1_ik), -4.0_wp, 1.0e-12_wp,                       &
                    'changing depth did not move the column base')
   call check(abs(p2%dz(1) - p%dz(1)) > 1.0e-9_wp,                                                 &
              'a deeper column must have thicker layers on the same grid growth')
   sc%depth = 2.0_wp

   !=== 3. n_layer reaches the active count, and the grid still spans the full depth. =======!
   sc%n_layer = 4_ik
   call build_column(sc, p2)
   call check(p2%n_active == 4_ik,                       'n_layer did not reach soil_params_t%n_active')
   call check_close(p2%soil_layer_z(5_ik), -2.0_wp, 1.0e-12_wp,                                    &
                    'a coarser grid must still reach the full depth')
   sc%n_layer = 10_ik

   !=== 4. Texture reaches the layers, and the DERIVED thresholds follow it. ================!
   sc%theta_sat = 0.55_wp ; sc%theta_res = 0.02_wp
   call build_column(sc, p2)
   call check_close(p2%theta_sat(1), 0.55_wp, 1.0e-14_wp, 'theta_sat did not reach the layers')
   call check_close(p2%theta_res(1), 0.02_wp, 1.0e-14_wp, 'theta_res did not reach the layers')
   call check(p2%theta_fc(1) > p2%theta_wp(1),           'field capacity must exceed wilting point')
   call check(p2%theta_fc(1) < p2%theta_sat(1),          'field capacity must be below saturation')
   sc%theta_sat = 0.43_wp ; sc%theta_res = 0.078_wp

   !=== 5. The retention family selects which pair curve_par_a/n mean. ======================!
   sc%retention = SOIL_RETENTION_CAMPBELL
   call build_column(sc, p2)
   call check(p2%retention == SOIL_RETENTION_CAMPBELL,   'retention family did not reach the column')
   sc%retention = SOIL_RETENTION_VG

   !=== 6. The thermal keys reach the thermal column. =======================================!
   sc%solid_conductivity = 2.5_wp ; sc%dry_conductivity = 0.20_wp ; sc%dry_heat_capacity = 1.5e6_wp
   call build_soil_therm_params(sc%n_layer, sc%solid_conductivity, sc%dry_conductivity,            &
                                sc%dry_heat_capacity, th)
   call check(th%nzg_active == sc%n_layer,               'n_layer did not reach the thermal column')
   call check_close(th%soil_solid_conductivity(1), 2.5_wp,  1.0e-14_wp, 'solid_conductivity did not reach it')
   call check_close(th%soil_dry_conductivity(1),   0.20_wp, 1.0e-14_wp, 'dry_conductivity did not reach it')
   call check_close(th%soil_dry_heat_capacity(1),  1.5e6_wp, 1.0e-8_wp, 'dry_heat_capacity did not reach it')

   print '(a)', 'test_soil_column_config: all checks passed'

contains

   !----- Exactly the call build_fast_context makes, so the test cannot drift from it. -------!
   subroutine build_column(c, out)
      type(soil_column_config_t), intent(in)  :: c
      type(soil_params_t),        intent(out) :: out
      call build_soil_hydr_params(c%n_layer, c%retention, c%depth, c%grid_growth,                  &
                                  c%theta_sat, c%theta_res, c%ksat, c%curve_par_a,                 &
                                  c%curve_par_n, c%root_beta, c%psi_fc, out)
   end subroutine build_column

end program test_soil_column_config
