!==========================================================================================!
! meds_test_support -- shared assert helpers for the CTest programs.                       !
!==========================================================================================!
module meds_test_support
   use meds_kinds,      only : wp, ik
   use meds_time,       only : meds_time_t
   use meds_allometry,  only : set_allometry
   use meds_pft_params, only : alloc_pft_table, derive_pft_rates
   use meds_config,     only : meds_config_t, derive_config, TS_DAILY, BK_SERIAL, INIT_BARE
   implicit none
   private

   public :: check, check_close, banner, build_test_config

contains

   !---------------------------------------------------------------------------------------!
   ! Build a COMPLETE test configuration in code. The model itself has no built-in defaults  !
   ! (every parameter comes from TOML), so the test fixtures live here, explicitly. Mirrors    !
   ! the canonical shipped config; pass ts_mode to vary the timestep.                          !
   !---------------------------------------------------------------------------------------!
   function build_test_config(ts_mode) result(cfg)
      integer(ik), intent(in), optional :: ts_mode
      type(meds_config_t)               :: cfg

      cfg%ts_mode    = TS_DAILY ; if (present(ts_mode)) cfg%ts_mode = ts_mode
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

      call alloc_pft_table(cfg%pft, 3_ik)
      associate (p => cfg%pft)
         p%wood_density = [ 0.40_wp, 0.60_wp, 0.85_wp ]
         p%dbh_critical = [ 100.0_wp, 100.0_wp, 100.0_wp ]
         p%growth_dbh_slope = [ 0.25_wp, 0.25_wp, 0.25_wp ]
         p%growth_dbh_cap   = [ 100.0_wp, 100.0_wp, 100.0_wp ]
         p%growth_dbh_max   = [ 1.5_wp, 1.0_wp, 0.5_wp ]
         p%growth_lai_slope = [ -0.8_wp, -0.7_wp, -0.6_wp ]
         p%reproduction_investment_fraction = [ 0.3_wp, 0.3_wp, 0.3_wp ]
         p%repro_carbon_efficiency = [ 1.0e-3_wp, 1.0e-3_wp, 1.0e-3_wp ]
         p%seed_rain_recruits = [ 0.01_wp, 0.01_wp, 0.01_wp ]
         p%include_pft = [ 1_ik, 1_ik, 1_ik ]
         p%min_cohort_height = 2.0_wp ; p%min_reproduction_height = 20.0_wp
         p%mort_rho_ref = 0.6_wp
         p%mort_gamma_0 = 0.0094_wp ; p%mort_gamma_exp = -1.8392_wp
         p%mort_alpha_0 = 0.05_wp   ; p%mort_alpha_exp = -1.1493_wp
         p%mort_beta_0  = 18.72_wp  ; p%mort_beta_exp  =  0.2792_wp
      end associate

      call set_allometry(1.139963_wp, 0.564899_wp, 46.0_wp, 0.06080334_wp, 1.0044785_wp,        &
                         0.370_wp, 0.464_wp, 0.46769540_wp, 0.6410495_wp, 0.5_wp)
      call derive_pft_rates(cfg%pft)
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
