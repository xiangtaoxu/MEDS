!==========================================================================================!
! test_pft_optics_config -- the [pft] canopy optics reach the two-stream's per-PFT table,    !
! and DIFFERENTIATE the PFTs.                                                                !
!                                                                                          !
! These fourteen traits were literals broadcast over every PFT, which meant two plant        !
! functional types could not differ in how they intercept or scatter light -- in a model      !
! whose two-stream solver exists to resolve that competition. Moving them to the table is      !
! only worth anything if the values arrive per PFT, and that wiring fails silently: a          !
! mis-threaded index leaves every PFT on the first one's optics, the radiation still solves,   !
! the energy budget still closes, and the canopy is quietly one optical type again. So this    !
! asserts the values arrive AND that PFTs differ, not merely that the code compiles.           !
!==========================================================================================!
program test_pft_optics_config
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_fast_dynamics,  only : fast_context_t, build_fast_context
   use meds_canopy_types,   only : RAD_VIS, RAD_NIR, RAD_LW
   use meds_test_support,   only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t)  :: cfg
   type(fast_context_t) :: ctx
   real(wp) :: om1, om2

   call banner('[pft] canopy optics: per PFT, and reaching the two-stream')
   cfg = build_test_config()

   !=== 1. The DEFAULTS are the literals the fast driver used to carry. =====================!
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

   !=== 5. Clumping and the leaf-angle distribution are per PFT too. ========================!
   cfg%pft%leaf_clumping(2)   = 0.40_wp
   cfg%pft%leaf_angle_mean(2) = 20.0_wp        ! erectophile against PFT 1's spherical 45 deg
   call build_fast_context(cfg, ctx)
   call check_close(ctx%rad_opt%clumping_leaf(2), 0.40_wp, 1.0e-14_wp, 'leaf clumping did not reach the table')
   call check_close(ctx%rad_opt%clumping_leaf(1), 0.80_wp, 1.0e-14_wp, 'PFT 1 clumping must be untouched')
   call check(abs(ctx%rad_opt%bf(2) - ctx%rad_opt%bf(1)) > 1.0e-6_wp,                            &
              'a different leaf-angle mean must give a different angle distribution')

   print '(a)', 'test_pft_optics_config: all checks passed'

end program test_pft_optics_config
