! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_pft_optics_config -- the [pft] canopy optics reach the two-stream's per-PFT table,    !
! and DIFFERENTIATE the PFTs.                                                                !
!                                                                                          !
! These fourteen traits are what lets two plant functional types differ in how they        !
! intercept or scatter light -- the competition the two-stream solver exists to resolve.   !
! They are only worth anything if the values arrive per PFT, and that wiring fails         !
! silently: a mis-threaded index leaves every PFT on the first one's optics, the radiation !
! still solves, the energy budget still closes, and the canopy is quietly one optical      !
! type. So this asserts the values arrive AND that PFTs differ, not merely that the code   !
! compiles. The same holds for the per-PFT hydraulics, down to a cohort's water at birth   !
! and the slow tier's shed of water above capacity.                                       !
!==========================================================================================!
program test_pft_optics_config
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t, pft_leaf_psi_tlp, soil_column_params
   use meds_config_io,      only : write_pft_params_csv
   use meds_plant_types,    only : hydro_args_table_t
   use meds_fast_types,     only : apply_hydraulics_config
   use meds_plant_hydr_lib, only : conductance_retained_fraction, tissue_psi_tlp,  &
                                   tissue_psi_to_water
   use meds_polygon_state_vars, only : polygon_t
   use meds_patch_state_vars, only : patch_physical_state
   use meds_soil_params,    only : soil_params_t
   use meds_init,           only : init_bare_ground, add_cohort, finalize_init
   use meds_vegetation_dynamics, only : shed_excess_tissue_water
   use meds_fast_dynamics,  only : fast_context_t, build_fast_context
   use meds_canopy_types,   only : RAD_VIS, RAD_NIR, RAD_LW
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t)  :: cfg
   type(fast_context_t) :: ctx
   real(wp) :: om1, om2

   call banner('[pft] canopy optics: per PFT, and reaching the two-stream')
   cfg = build_test_config()

   !=== 1. The test configuration's DEFAULT optics. ========================================!
   call check_close(cfg%pft%leaf_reflect_vis(1),  0.10_wp,  1.0e-14_wp, 'default leaf VIS reflectance changed')
   call check_close(cfg%pft%leaf_transmit_nir(1), 0.25_wp,  1.0e-14_wp, 'default leaf NIR transmittance changed')
   call check_close(cfg%pft%leaf_emissivity(1),   0.97_wp,  1.0e-14_wp, 'default leaf emissivity changed')
   call check_close(cfg%pft%wood_emissivity(1),   0.90_wp,  1.0e-14_wp, 'default wood emissivity changed')
   call check_close(cfg%pft%leaf_clumping(1),     0.80_wp,  1.0e-14_wp, 'default leaf clumping changed')
   call check_close(cfg%pft%leaf_angle_mean(1),   45.0_wp,  1.0e-14_wp, 'default leaf-angle mean changed')

   !=== 2. Uniform traits => a uniform table. The control for test 3. =======================!
   call build_fast_context(cfg, ctx)
   call check(ctx%rad_opt%n_pft >= 2_ik, 'fixture must have at least two PFTs to differentiate')
   call check_close(ctx%rad_opt%omega_leaf(RAD_VIS,1), ctx%rad_opt%omega_leaf(RAD_VIS,2), 1.0e-14_wp,       &
                    'identical traits must give identical optics')

   !=== 3. Differentiate ONE PFT and the table must follow it -- only that one. ==============!
   !       Single-scatter albedo is reflectance + transmittance, so raising the reflectance    !
   !       of PFT 2 must raise ITS omega and leave PFT 1 alone.                                !
   om1 = ctx%rad_opt%omega_leaf(RAD_VIS,1)
   cfg%pft%leaf_reflect_vis(2) = 0.30_wp
   call build_fast_context(cfg, ctx)
   om2 = ctx%rad_opt%omega_leaf(RAD_VIS,2)
   call check_close(ctx%rad_opt%omega_leaf(RAD_VIS,1), om1, 1.0e-14_wp,                                    &
                    'changing PFT 2 must not move PFT 1 (per-PFT indexing)')
   call check_close(om2, 0.30_wp + cfg%pft%leaf_transmit_vis(2), 1.0e-12_wp,                               &
                    'leaf VIS omega must be reflectance + transmittance of THAT PFT')
   call check(om2 > om1, 'a more reflective PFT must have the larger single-scatter albedo')

   !=== 4. Longwave is configured as EMISSIVITY: omega = 1 - emissivity, transmittance 0. ====!
   cfg%pft%leaf_emissivity(2) = 0.80_wp
   call build_fast_context(cfg, ctx)
   call check_close(ctx%rad_opt%omega_leaf(RAD_LW,2), 1.0_wp - 0.80_wp, 1.0e-12_wp,                        &
                    'leaf LW omega must be 1 - emissivity (a leaf is opaque in the thermal band)')

   !=== 4b. The ground optics come from [soil], not a constant in the fast context. =============!
   cfg%soil%ground_albedo_vis = 0.11_wp ; cfg%soil%ground_albedo_nir = 0.22_wp ; cfg%soil%ground_emissivity = 0.93_wp
   call build_fast_context(cfg, ctx)
   call check_close(ctx%soil_albedo(1), 0.11_wp, 0.0_wp, 'soil.ground_albedo_vis did not reach the radiation solver')
   call check_close(ctx%soil_albedo(2), 0.22_wp, 0.0_wp, 'soil.ground_albedo_nir did not reach the radiation solver')
   call check_close(ctx%soil_emiss,     0.93_wp, 0.0_wp, 'soil.ground_emissivity did not reach the radiation solver')
   cfg%soil%ground_albedo_vis = 0.15_wp ; cfg%soil%ground_albedo_nir = 0.30_wp ; cfg%soil%ground_emissivity = 0.95_wp

   !=== 5. Clumping and the leaf-angle distribution are per PFT too. ========================!
   cfg%pft%leaf_clumping(2)   = 0.40_wp
   cfg%pft%leaf_angle_mean(2) = 20.0_wp        ! erectophile against PFT 1's 45 deg (spherical is 57.3)
   call build_fast_context(cfg, ctx)
   call check_close(ctx%rad_opt%clumping_leaf(2), 0.40_wp, 1.0e-14_wp, 'leaf clumping did not reach the table')
   call check_close(ctx%rad_opt%clumping_leaf(1), 0.80_wp, 1.0e-14_wp, 'PFT 1 clumping must be untouched')
   call check(abs(ctx%rad_opt%bf(2) - ctx%rad_opt%bf(1)) > 1.0e-6_wp,                            &
              'a different leaf-angle mean must give a different angle distribution')

   !=== 5b. PER-PFT HYDRAULIC TRAITS. Same shape as the optics checks above: each PFT's table   !
   !        entry carries its own [hydraulics] record, and a different record CHANGES the answer,
   !        not merely that it compiles.
   block
      type(hydro_args_table_t)   :: tab
      real(wp) :: plc1, plc2
      cfg%pft%hydro(2)%wood_psi50     = -0.8_wp   ! vs PFT 1's -2.0
      cfg%pft%hydro(2)%leaf_curve%pi0 = -2.5_wp   ! vs PFT 1's -1.5
      cfg%pft%hydro(2)%wood_kmax      =  2.0_wp
      call apply_hydraulics_config(cfg%pft, tab)
      call check_close(tab%pft(1)%wood_psi50, cfg%pft%hydro(1)%wood_psi50, 0.0_wp,                &
                       'PFT 1 keeps its own psi50')
      call check_close(tab%pft(2)%wood_psi50, -0.8_wp,   0.0_wp, 'PFT 2 takes its own psi50')
      call check_close(tab%pft(2)%leaf_curve%pi0, -2.5_wp,   0.0_wp, 'PFT 2 takes its own leaf pi0')
      call check_close(tab%pft(2)%wood_kmax, 2.0_wp, 0.0_wp, 'PFT 2 takes its own conductivity')
      call check_close(tab%pft(1)%wood_kmax, cfg%pft%hydro(1)%wood_kmax, 0.0_wp,                  &
                       'PFT 1 keeps its own conductivity')
      call check_close(tab%pft(2)%vessel_curl, cfg%pft%hydro(2)%vessel_curl, 0.0_wp,              &
                       'every trait of the table is the PFT''s own')
      !----- And it MATTERS: at the same potential the more vulnerable xylem has lost more. -----!
      plc1 = 1.0_wp - conductance_retained_fraction(-1.5_wp, tab%pft(1)%wood_psi50, tab%pft(1)%wood_kexp)
      plc2 = 1.0_wp - conductance_retained_fraction(-1.5_wp, tab%pft(2)%wood_psi50, tab%pft(2)%wood_kexp)
      call check(plc2 > plc1, 'a more vulnerable PFT loses more conductance at the same psi')
      !----- The fast water-stress limb (the turgor-loss point to twice it) follows the PFT's OWN   !
      !      pressure-volume curve, the one the plant-water solver stores its leaf water on. -------!
      call build_fast_context(cfg, ctx)
      call check_close(ctx%col_config%leaf_photo%pft(2)%psi_open,                                 &
                       tissue_psi_tlp(-2.5_wp, cfg%pft%hydro(2)%leaf_curve%elastic_mod), 1.0e-12_wp,   &
                       'PFT 2 loses capacity from its own turgor-loss point')
      call check_close(ctx%col_config%leaf_photo%pft(2)%psi_close,                                &
                       2.0_wp * tissue_psi_tlp(-2.5_wp, cfg%pft%hydro(2)%leaf_curve%elastic_mod),        &
                       1.0e-12_wp, 'PFT 2 has lost it at twice its turgor-loss point')
      call check_close(ctx%col_config%leaf_photo%pft(1)%psi_open,                                 &
                       tissue_psi_tlp(cfg%pft%hydro(1)%leaf_curve%pi0, cfg%pft%hydro(1)%leaf_curve%elastic_mod), &
                       1.0e-12_wp, 'PFT 1 keeps its own turgor-loss point')
      call check_close(pft_leaf_psi_tlp(cfg, 2_ik), ctx%col_config%leaf_photo%pft(2)%psi_open, 0.0_wp, &
                       'the phenology reads the same per-PFT turgor-loss point')
   end block

   !=== A cohort is born with its tissue water on ITS PFT's curve, the one the fast loop reads    !
   !    the water back on, and the slow tier caps a store on that curve too. ======================!
   block
      type(hydro_args_table_t)   :: tab
      type(polygon_t) :: poly
      real(wp) :: w_expect, w_shed, t_birth, psi_birth
      integer(ik) :: i2
      type(soil_params_t) :: soil
      cfg%pft%hydro(2)%leaf_curve%water_sat = 3.0_wp   ! vs PFT 1's 2.0
      call apply_hydraulics_config(cfg%pft, tab)
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 20.0_wp)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 20.0_wp)
      call finalize_init(poly)
      i2 = merge(1_ik, 2_ik, poly%cohort%pft(1) == 2_ik)        ! the PFT-2 cohort, after the sort
      call patch_physical_state(poly%patch%cas(1), poly%patch%soil_water(1), cfg, poly%site%surface(1),     &
                                t_birth, psi_birth)
      w_expect = tissue_psi_to_water(psi_birth, tab%pft(2)%leaf_curve, poly%cohort%leaf_carbon(i2))
      call check_close(poly%cohort%leaf_water_mass(i2), w_expect, 1.0e-12_wp,                      &
                       'PFT 2 leaf water at birth is on its own curve')
      call check(abs(w_expect - tissue_psi_to_water(psi_birth, tab%pft(1)%leaf_curve, poly%cohort%leaf_carbon(i2))) &
                 > 1.0e-6_wp * w_expect, 'the two PFTs'' curves give different birth water (the fixture bites)')
      !----- Capped on its own ceiling too: a supersaturated store falls to 3.0 x leaf carbon. ---!
      poly%cohort%leaf_water_mass(i2) = 10.0_wp * poly%cohort%leaf_carbon(i2)
      call shed_excess_tissue_water(poly, cfg, w_shed)
      call check_close(poly%cohort%leaf_water_mass(i2), 3.0_wp * poly%cohort%leaf_carbon(i2), 1.0e-12_wp, &
                       'PFT 2 leaf water is capped at its own saturation')
   end block

   !=== 6. The PFT-parameter CSV dump has as many values as it has column headers. ==========!
   !        A Fortran format SHORTER than its output list does not fail: it REVERTS to the last  !
   !        repeat group and keeps going, so an integer edit descriptor silently receives a real  !
   !        and prints its bit pattern, and the trailing columns disappear. Adding PFT traits     !
   !        without widening the format does exactly that, and nothing else in the suite notices, !
   !        because the run completes and every budget still closes. Counting the two lines is   !
   !        the whole test.                                                                       !
   block
      character(len=*), parameter :: CSVFILE = 'test_pft_params_tmp.csv'
      character(len=8192) :: head_line, val_line
      integer(ik)         :: u, ios, nhead, nval
      call build_fast_context(cfg, ctx)          ! keep cfg consistent with the edits above
      call write_pft_params_csv(cfg, CSVFILE)
      open(newunit=u, file=CSVFILE, status='old', action='read', iostat=ios)
      call check(ios == 0, 'PFT parameter CSV must be written')
      read(u,'(a)') head_line
      read(u,'(a)') val_line
      close(u, status='delete')
      nhead = count_fields(head_line)
      nval  = count_fields(val_line)
      call check(nhead == nval, 'PFT CSV: header column count must equal the value count')
      if (nhead /= nval) print '(a,i0,a,i0)', '   header fields = ', nhead, ', value fields = ', nval
      !----- And no field may be a bit pattern: every value is either an integer or es15.8, so a  !
      !      19-digit run of digits is the format-reversion signature, not a real number. --------!
      call check(index(val_line, '4605380978949069210') == 0,                                     &
                 'PFT CSV: a real printed through an integer descriptor (format reversion)')
   end block

   print '(a)', 'test_pft_optics_config: all checks passed'

contains

   pure integer(ik) function count_fields(line) result(n)
      character(len=*), intent(in) :: line
      integer(ik) :: i
      n = 1_ik
      do i = 1_ik, len_trim(line)
         if (line(i:i) == ',') n = n + 1_ik
      end do
   end function count_fields

end program test_pft_optics_config
