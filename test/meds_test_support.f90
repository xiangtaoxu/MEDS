!==========================================================================================!
! meds_test_support -- shared assert helpers for the CTest programs.                       !
!==========================================================================================!
module meds_test_support
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : day_sec
   use meds_time,       only : meds_time_t
   use meds_allometry,  only : set_allometry
   use meds_pft_params, only : alloc_pft_table, derive_pft_rates, derive_leaf_params,          &
                               PATH_C3, PATH_C4
   use meds_config,     only : meds_config_t, derive_config, BK_SERIAL, INIT_BARE,             &
                               SM_MEDLYN, TRESP_PEAKED, COLIM_QUADRATIC,                     &
                               SCHEME_SPLIT_SEQUENTIAL, INTEG_SPLIT
   implicit none
   private

   public :: check, check_close, banner, build_test_config

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
      cfg%dt_fast            = 900.0_wp
      cfg%integration_scheme = SCHEME_SPLIT_SEQUENTIAL
      cfg%time_integrator    = INTEG_SPLIT
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
      cfg%patch_light_tol = 0.10_wp ; cfg%patch_light_maxdev_factor = 1.5_wp
      cfg%patch_diff_age_tol = 1.0_wp ; cfg%min_patch_area = 1.0e-4_wp
      cfg%patch_min_area_remain = 0.99_wp ; cfg%enable_patch_fission = .false.

      cfg%patch_disturbance_rate = 0.014_wp ; cfg%disturbance_survive_height = 10.0_wp
      cfg%growth_memory_days = 90.0_wp ; cfg%min_recruit_size = 1.0e-2_wp ; cfg%conservation_tol = 0.001_wp

      cfg%init_mode = INIT_BARE ; cfg%init_restart_file = '' ; cfg%init_census_file = '' ; cfg%pft_config = ''
      cfg%io_write_output = .false. ; cfg%io_output_dir = '.' ; cfg%io_output_prefix = 'test'
      cfg%io_output_interval_years = 1_ik ; cfg%io_cohort_max = 2048_ik ; cfg%io_patch_max = 64_ik
      cfg%io_write_state = .false. ; cfg%io_state_interval_years = 50_ik
      cfg%override_derived = .false.
      cfg%gpp_ref = 0.3_wp

      !----- Leaf physiology: model selectors + shared biochemistry (mirrors the TOML). ---!
      cfg%stomatal_model = SM_MEDLYN ; cfg%temp_response_form = TRESP_PEAKED
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
         p%dbh_critical = [ 100.0_wp, 100.0_wp, 100.0_wp ]
         p%hgt_max      = [ 46.0_wp, 46.0_wp, 46.0_wp ]
         p%reproduction_investment_fraction = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%repro_carbon_efficiency = [ 1.0e-3_wp, 1.0e-3_wp, 1.0e-3_wp ]
         p%seed_rain_recruits = [ 0.01_wp, 0.01_wp, 0.01_wp ]
         p%include_pft = [ 1_ik, 1_ik, 1_ik ]
         p%min_cohort_height = 2.0_wp ; p%min_reproduction_height = 20.0_wp
         p%mort_rho_ref = 0.6_wp
         p%mort_gamma_0 = 0.0094_wp ; p%mort_gamma_exp = -1.8392_wp
         p%mort_alpha_0 = 0.05_wp   ; p%mort_alpha_exp = -1.1493_wp
         p%mort_beta_0  = 18.72_wp  ; p%mort_beta_exp  =  0.2792_wp
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
         p%theta_cj_c4      = [ 0.80_wp, 0.80_wp, 0.80_wp ]
         p%theta_ic_c4      = [ 0.95_wp, 0.95_wp, 0.95_wp ]
         p%katul_lambda25   = [ 600.0_wp, 800.0_wp, 350.0_wp ]
         p%wstress_psi_open = [ -0.5_wp, -0.5_wp, -0.5_wp ]
         p%wstress_psi_close= [ -2.5_wp, -3.0_wp, -2.0_wp ]
         p%wstress_lambda_exp = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         !----- Carbon-dynamics traits (mirrors meds_config_pft.toml). -----------------------!
         p%sla                    = [ 16.0_wp, 13.0_wp, 10.0_wp ]
         p%root_to_leaf_ratio     = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%huber_value            = [ 1.0e-4_wp, 1.5e-4_wp, 2.0e-4_wp ]
         p%aboveground_frac       = [ 0.7_wp, 0.7_wp, 0.7_wp ]
         p%storage_cushion        = [ 1.0_wp, 1.0_wp, 1.0_wp ]
         p%growth_resp_factor     = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%leaf_turnover_rate     = [ 1.0_wp, 0.5_wp, 0.33_wp ]
         p%fineroot_turnover_rate = [ 1.0_wp, 0.8_wp, 0.6_wp ]
         p%wood_carbon_density    = [ 200.0_wp, 300.0_wp, 425.0_wp ]
         p%evergreen              = [ 1_ik, 1_ik, 1_ik ]
      end associate

      call set_allometry(1.139963_wp, 0.564899_wp, 0.06080334_wp, 1.0044785_wp,        &
                         0.370_wp, 0.464_wp, 0.46769540_wp, 0.6410495_wp, 0.5_wp)
      call derive_pft_rates(cfg%pft)
      call derive_leaf_params(cfg%pft)
      call derive_config(cfg)
   end function build_test_config

   subroutine check(cond, msg)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: msg
      if (.not. cond) then
         write(*,'(2a)') 'FAIL: ', msg
         error stop 1
      end if
   end subroutine check

   subroutine check_close(a, b, rtol, msg)
      real(wp),         intent(in) :: a, b, rtol
      character(len=*), intent(in) :: msg
      real(wp) :: tol
      tol = rtol * max(abs(b), 1.0e-30_wp) + 1.0e-12_wp
      if (abs(a - b) > tol) then
         write(*,'(2a)') 'FAIL: ', msg
         write(*,'(a,es16.8,a,es16.8,a,es10.2)') '   got=', a, ' expected=', b, ' rtol=', rtol
         error stop 1
      end if
   end subroutine check_close

   subroutine banner(name)
      character(len=*), intent(in) :: name
      write(*,'(2a)') '[test] ', name
   end subroutine banner

end module meds_test_support
