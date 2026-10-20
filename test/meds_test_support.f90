! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_test_support -- the in-code test CONFIGURATION fixture, plus a re-export of the      !
! assertion helpers so a test needs one `use` rather than two. The assertions themselves live !
! in meds_test_assert, which depends on nothing but a kind; see that module for why, and for  !
! the two assertion families (fatal and accumulating).                                     !
!==========================================================================================!
module meds_test_support
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : day_sec
   use meds_time,       only : meds_time_t
   use meds_pft_params, only : pft_table_t, alloc_pft_table, PATH_C3, PATH_C4
   use meds_allometry_lib,  only : allometry_t
   use meds_plant_hydr_lib, only : plant_hydraulics_t, water_curve_t
   use meds_config,     only : meds_config_t, derive_parameters, BK_SERIAL, INIT_BARE, INTEG_ARK
   use meds_leaf_opts,     only : SM_MEDLYN, COLIM_QUADRATIC
   use meds_therm_lib,     only : temp_response_peaked
   use meds_test_assert, only : check, check_close, check_true, check_int, banner,             &
                                test_report, test_reset, test_failures
   implicit none
   private

   !----- Re-exported, so a test that links the full support library needs one `use`. -----!
   public :: check, check_close, check_true, check_int, banner
   public :: test_report, test_reset, test_failures
   public :: build_test_config, test_pft_table

contains

   !---------------------------------------------------------------------------------------!
   ! Build a COMPLETE test configuration in code. The model itself has no built-in defaults  !
   ! (every parameter comes from TOML), so the test fixtures live here, explicitly. Mirrors    !
   ! the canonical shipped config; pass dt_slow [s] to vary the slow-process timestep.          !
   !---------------------------------------------------------------------------------------!
   function build_test_config(dt_slow) result(cfg)
      real(wp), intent(in), optional :: dt_slow
      type(meds_config_t)            :: cfg

      cfg%dt_slow    = day_sec ; if (present(dt_slow)) cfg%dt_slow = dt_slow
      cfg%fast_biophysics_on = .false.
      cfg%dt_fast            = 150.0_wp
      cfg%time_integrator    = INTEG_ARK
      cfg%backend    = BK_SERIAL
      cfg%start_time = meds_time_t(2000_ik, 1_ik, 1_ik)
      cfg%end_time   = meds_time_t(2100_ik, 1_ik, 1_ik)
      cfg%demography_on = .true.
      cfg%do_cohort_fissfuse = .true. ; cfg%do_patch_fissfuse = .true. ; cfg%do_patch_disturbance = .true.

      cfg%max_cohort = 60_ik ; cfg%n_cohort_fusion_iter = 6_ik
      cfg%cohort_size_tol_min = 0.02_wp ; cfg%cohort_size_tol_max = 0.10_wp
      cfg%cohort_lai_cap = 1.0_wp ; cfg%min_cohort_agb = 1.0e-6_wp
      cfg%negligible_nplant = 1.0e-8_wp ; cfg%split_eps = 1.0e-4_wp ; cfg%enable_cohort_fission = .true.

      cfg%n_height_layers = 16_ik ; cfg%max_patch = 12_ik ; cfg%n_patch_fusion_iter = 6_ik
      cfg%patch_light_tol = 0.10_wp ; cfg%patch_light_tol_max = 0.15_wp ; cfg%patch_light_maxdev_factor = 1.5_wp
      cfg%patch_diff_age_tol = 1.0_wp ; cfg%min_patch_area = 1.0e-4_wp
      cfg%patch_min_area_remain = 0.99_wp ; cfg%enable_patch_fission = .false.

      cfg%patch_disturbance_rate = 0.014_wp ; cfg%disturbance_survive_height = 10.0_wp
      cfg%growth_memory_days = 90.0_wp ; cfg%min_recruit_size = 1.0e-2_wp ; cfg%conservation_tol = 0.001_wp

      cfg%init_mode = INIT_BARE ; cfg%init_restart_file = '' ; cfg%init_census_file = '' ; cfg%pft_config = ''
      cfg%state_output_dir = '.' ; cfg%state_output_prefix = 'test'
      cfg%state_write_state = .false. ; cfg%state_interval_years_cfg = 50_ik
      cfg%gpp_ref = 0.3_wp

      !----- Leaf physiology: model selectors + shared biochemistry (mirrors the TOML). ---!
      cfg%stomatal_model = SM_MEDLYN ; cfg%temp_response_form = temp_response_peaked
      cfg%colimitation = COLIM_QUADRATIC ; cfg%leaf_use_boundary_layer = .false.
      cfg%kc25 = 40.49_wp ; cfg%ko25 = 27840.0_wp ; cfg%gstar25 = 4.275_wp
      cfg%ea_kc = 79430.0_wp ; cfg%ea_ko = 36380.0_wp ; cfg%ea_gstar = 37830.0_wp
      cfg%ea_vcmax = 65330.0_wp ; cfg%ea_jmax = 43540.0_wp ; cfg%ea_rd = 46390.0_wp
      cfg%hd_vcmax = 200000.0_wp ; cfg%hd_jmax = 200000.0_wp ; cfg%hd_rd = 200000.0_wp
      cfg%ds_vcmax = 650.0_wp ; cfg%ds_jmax = 640.0_wp ; cfg%ds_rd = 490.0_wp
      cfg%o2_mol_frac = 0.209_wp ; cfg%leaf_absorptance = 0.85_wp ; cfg%phi_psii = 0.85_wp

      call alloc_pft_table(cfg%pft, 3_ik)
      associate (p => cfg%pft)
         p%wood_density = [ 0.40_wp, 0.60_wp, 0.85_wp ]
         p%reproduction_investment_fraction = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%recruit_carbon_efficiency = [ 1.0e-3_wp, 1.0e-3_wp, 1.0e-3_wp ]
         p%seed_rain_recruits = [ 0.01_wp, 0.01_wp, 0.01_wp ]
         p%include_pft = [ 1_ik, 1_ik, 1_ik ]
         p%min_cohort_height = 2.0_wp
         p%min_repro_height  = [ 20.0_wp, 20.0_wp, 20.0_wp ]
         !----- The mortality hazard and the plasticity slopes (mirrors meds_config_pft.toml). --!
         p%mort_gamma = [ 0.019815037439806053_wp, 0.0094_wp, 0.004953547651645506_wp ]
         p%mort_alpha = [ 0.07968043353213256_wp, 0.05_wp, 0.0335056513529311_wp ]
         p%mort_beta  = [ 16.71633880800162_wp, 18.72_wp, 20.631926375047502_wp ]
         p%kplastic_vcmax = [ -0.1459341188305866_wp, -0.11760191025318557_wp, -0.10943768264434321_wp ]
         p%kplastic_rd    = p%kplastic_vcmax
         p%kplastic_sla   = [ 0.15654750929099287_wp, 0.12615477645341722_wp, 0.11739678683665906_wp ]
         p%kplastic_leaf_turnover = [ -0.058535787713143994_wp, -0.015560662518427415_wp, 0.0_wp ]
         !----- Leaf-photosynthesis per-PFT traits (PFT 3 is C4). ---------------------------!
         p%photosynthetic_pathway = [ PATH_C3, PATH_C3, PATH_C4 ]
         p%vcmax25          = [ 60.0_wp, 45.0_wp, 40.0_wp ]
         p%jmax_vcmax_ratio = [ 1.8_wp, 1.7_wp, 4.0_wp ]
         p%tpu_vcmax_ratio  = [ 0.167_wp, 0.167_wp, 0.167_wp ]
         p%rd_vcmax_ratio   = [ 0.015_wp, 0.015_wp, 0.025_wp ]
         p%kp25             = [ 1.0e-9_wp, 1.0e-9_wp, 0.7_wp ]
         p%stomatal_g0      = [ 0.01_wp, 0.01_wp, 0.04_wp ]
         p%stomatal_g1      = [ 4.0_wp, 3.0_wp, 1.6_wp ]
         p%stomatal_d0      = [ 1500.0_wp, 1500.0_wp, 1500.0_wp ]
         p%quantum_yield_c4 = [ 0.0_wp, 0.0_wp, 0.04_wp ]
         p%theta_j          = [ 0.85_wp, 0.90_wp, 0.85_wp ]
         p%theta_cj_c3      = [ 0.98_wp, 0.98_wp, 0.98_wp ]
         p%theta_ip_c3      = [ 0.95_wp, 0.95_wp, 0.95_wp ]
         p%theta_cj_c4      = [ 0.80_wp, 0.80_wp, 0.80_wp ]
         p%theta_ic_c4      = [ 0.95_wp, 0.95_wp, 0.95_wp ]
         p%katul_lambda25   = [ 600.0_wp, 800.0_wp, 350.0_wp ]
         p%wstress_lambda_exp = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%wstress_sref_stomata = [ 2.0_wp, 2.0_wp, 2.0_wp ]
         !----- Carbon-dynamics traits (mirrors meds_config_pft.toml). -----------------------!
         p%sla                    = [ 16.0_wp, 13.0_wp, 10.0_wp ]
         !----- Canopy-element geometry. --------------------------------------------------!
         p%leaf_width             = [ 0.04_wp, 0.04_wp, 0.04_wp ]
         p%branch_diameter        = [ 0.02_wp, 0.02_wp, 0.02_wp ]
         p%root_to_leaf_ratio     = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%huber_value            = [ 1.0e-4_wp, 1.5e-4_wp, 2.0e-4_wp ]
         !----- Non-leaf maintenance respiration. -----------------------------------------!
         p%is_woody               = [ .true., .true., .true. ]
         p%stem_resp_factor25     = [ 0.06_wp, 0.06_wp, 0.06_wp ]
         p%root_resp_factor25     = [ 0.30_wp, 0.30_wp, 0.30_wp ]
         !----- Canopy optics. ------------------------------------------------------------!
         p%leaf_reflect_vis  = [ 0.10_wp, 0.10_wp, 0.10_wp ]
         p%leaf_transmit_vis = [ 0.05_wp, 0.05_wp, 0.05_wp ]
         p%leaf_reflect_nir  = [ 0.45_wp, 0.45_wp, 0.45_wp ]
         p%leaf_transmit_nir = [ 0.25_wp, 0.25_wp, 0.25_wp ]
         p%leaf_emissivity   = [ 0.97_wp, 0.97_wp, 0.97_wp ]
         p%wood_reflect_vis  = [ 0.11_wp, 0.11_wp, 0.11_wp ]
         p%wood_transmit_vis = [ 0.001_wp, 0.001_wp, 0.001_wp ]
         p%wood_reflect_nir  = [ 0.25_wp, 0.25_wp, 0.25_wp ]
         p%wood_transmit_nir = [ 0.001_wp, 0.001_wp, 0.001_wp ]
         p%wood_emissivity   = [ 0.90_wp, 0.90_wp, 0.90_wp ]
         p%leaf_clumping     = [ 0.80_wp, 0.80_wp, 0.80_wp ]
         p%wood_clumping     = [ 0.50_wp, 0.50_wp, 0.50_wp ]
         p%leaf_angle_mean   = [ 45.0_wp, 45.0_wp, 45.0_wp ]
         p%leaf_angle_std    = [ 20.0_wp, 20.0_wp, 20.0_wp ]
         p%storage_cushion        = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%growth_resp_factor     = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%leaf_turnover_rate_toc = [ 1.0_wp, 0.5_wp, 0.33_wp ]          ! [1/yr]
         p%fineroot_turnover_rate = [ 1.0_wp, 0.8_wp, 0.6_wp ]
         p%wood_carbon_density    = [ 200.0_wp, 300.0_wp, 425.0_wp ]
         p%f_labile_leaf          = [ 0.80_wp, 0.70_wp, 0.60_wp ]
         p%f_labile_stem          = [ 0.05_wp, 0.05_wp, 0.05_wp ]
         p%struct_lignin_frac     = [ 0.20_wp, 0.25_wp, 0.30_wp ]
      end associate

      !----- Every PFT on the shipped allometry (mirrors meds_config_pft.toml). ---------------!
      cfg%pft%allom = allometry_t(hgt_max = 46.0_wp, b1Ht = 1.139963_wp, b2Ht = 0.564899_wp,          &
                                  agb_c1 = 0.03365_wp, agb_c2 = 0.976_wp, lai_b1 = 0.23384770_wp,     &
                                  lai_b2 = 0.6410495_wp, wai_b1 = 0.0096_wp, wai_b2 = 2.0947_wp,       &
                                  sapwood_area_b1 = 1.582_wp, sapwood_area_b2 = 1.764_wp,             &
                                  aboveground_frac = 0.7_wp)
      !----- ... and on the shipped hydraulics (mirrors meds_config_pft.toml). ------------------!
      cfg%pft%hydro = plant_hydraulics_t(                                                            &
         leaf_curve = water_curve_t(pi0 = -1.5_wp, elastic_mod = 12.0_wp, apoplast_frac = 0.30_wp,   &
                                    water_sat = 2.0_wp),                                             &
         wood_curve = water_curve_t(pi0 = -1.0_wp, elastic_mod = 8.0_wp, apoplast_frac = 0.20_wp,    &
                                    water_sat = 1.0_wp),                                             &
         wood_psi50 = -2.0_wp, wood_kexp = 2.0_wp, wood_kmax = 0.0_wp, vessel_curl = 1.5_wp,          &
         root_beta = 0.1_wp, root_depth_b1 = 1.1140580_wp, root_depth_b2 = 0.4223014_wp,             &
         specific_root_length = 2.0e4_wp, fine_root_radius = 2.5e-4_wp)
      cfg%pft%hydro%wood_kmax = [ 5.340397027887704_wp, 3.2683981418222547_wp, 1.7692398658114543_wp ]
      cfg%light_ext_coeff_fusefiss = 0.5_wp
      call derive_parameters(cfg)
   end function build_test_config

   !----- The test configuration's PFT table, for a fixture site that has no config of its own. --!
   function test_pft_table() result(pft)
      type(pft_table_t)   :: pft
      type(meds_config_t) :: cfg
      cfg = build_test_config()
      pft = cfg%pft
   end function test_pft_table

end module meds_test_support
