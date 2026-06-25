!==========================================================================================!
! meds_config_io -- build a run configuration, optionally overridden by a meds_config.toml   !
! file. Lives in src/io/ (config input) but has NO external dependency, so it is always       !
! compiled (unlike the netCDF output layer, which is gated behind MEDS_ENABLE_IO).            !
!                                                                                          !
! Flow: read the time-step from [run] first (build_config needs it to set dt and tolerances), !
! build the defaults, then overlay every recognized key. PFT arrays ([pft].wood_density,      !
! dbh_critical, ...) overlay the default trait table; if wood_density changes, the            !
! growth/mortality traits are re-derived (derive_pft_rates) so the trade-off axis stays        !
! consistent. Unknown keys are ignored; a missing file just keeps the built-in defaults.      !
!==========================================================================================!
module meds_config_io
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, build_config, validate_config,                  &
                               TS_DAILY, TS_WEEKLY, TS_MONTHLY
   use meds_pft_params, only : derive_pft_rates
   use meds_toml,       only : toml_table_t, toml_parse_file, toml_has, toml_int, toml_real,  &
                               toml_logical, toml_string, toml_real_array
   implicit none
   private

   public :: load_meds_config

contains

   !---------------------------------------------------------------------------------------!
   ! Load configuration from `path`. If the file is absent, return the built-in defaults     !
   ! (found=.false.). `run_years` is the simulation length read from [run].years.            !
   !---------------------------------------------------------------------------------------!
   subroutine load_meds_config(path, cfg, run_years, found)
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(out) :: cfg
      integer(ik),         intent(out) :: run_years
      logical,             intent(out) :: found
      type(toml_table_t)    :: t
      character(len=256)    :: ts
      integer(ik)           :: ts_mode, nout, npft
      real(wp), allocatable :: arr(:)

      run_years = 60_ik

      call toml_parse_file(path, t, found)
      if (.not. found) then
         cfg = build_config()
         call validate_config(cfg)
         return
      end if

      !----- Time step first (build_config derives dt and the tolerance schedule from it). -!
      ts = toml_string(t, 'run.timestep', 'daily')
      select case (trim(ts))
      case ('monthly') ; ts_mode = TS_MONTHLY
      case ('weekly')  ; ts_mode = TS_WEEKLY
      case default     ; ts_mode = TS_DAILY
      end select
      cfg = build_config(ts_mode = ts_mode)

      run_years = toml_int(t, 'run.years', run_years)

      !----- Cohort/patch structural tunables. --------------------------------------------!
      cfg%vegetation_dynamics_on = toml_logical(t, 'demography.vegetation_dynamics_on', cfg%vegetation_dynamics_on)
      cfg%do_cohort_fissfuse     = toml_logical(t, 'demography.do_cohort_fissfuse',     cfg%do_cohort_fissfuse)
      cfg%do_patch_fissfuse      = toml_logical(t, 'demography.do_patch_fissfuse',      cfg%do_patch_fissfuse)
      cfg%do_patch_disturbance   = toml_logical(t, 'demography.do_patch_disturbance',   cfg%do_patch_disturbance)
      cfg%max_cohort             = toml_int (t, 'demography.max_cohort',     cfg%max_cohort)
      cfg%cohort_lai_cap         = toml_real(t, 'demography.cohort_lai_cap', cfg%cohort_lai_cap)
      cfg%min_cohort_agb         = toml_real(t, 'demography.min_cohort_agb', cfg%min_cohort_agb)
      cfg%max_patch              = toml_int (t, 'demography.max_patch',      cfg%max_patch)
      cfg%patch_light_tol        = toml_real(t, 'demography.patch_light_tol', cfg%patch_light_tol)
      cfg%patch_light_maxdev_factor = toml_real(t, 'demography.patch_light_maxdev_factor', cfg%patch_light_maxdev_factor)
      cfg%conservation_tol       = toml_real(t, 'demography.conservation_tol', cfg%conservation_tol)

      !----- Disturbance. -----------------------------------------------------------------!
      cfg%patch_disturbance_rate     = toml_real(t, 'disturbance.patch_disturbance_rate',     cfg%patch_disturbance_rate)
      cfg%disturbance_survive_height = toml_real(t, 'disturbance.disturbance_survive_height', cfg%disturbance_survive_height)

      !----- Recruitment (scalar recruit_dens applies to all PFTs unless [pft] overrides). -!
      cfg%min_recruit_size            = toml_real(t, 'recruitment.min_recruit_size', cfg%min_recruit_size)
      cfg%pft%min_reproduction_height = toml_real(t, 'recruitment.min_reproduction_height',   &
                                                  cfg%pft%min_reproduction_height)
      if (toml_has(t, 'recruitment.recruit_dens'))                                            &
         cfg%pft%recruit_dens = toml_real(t, 'recruitment.recruit_dens', cfg%pft%recruit_dens(1))

      !----- Initial conditions (empty => bare ground; a path => restart from a census CSV).-!
      cfg%init_census_file = toml_string(t, 'init.census_file', cfg%init_census_file)

      !----- netCDF output. ---------------------------------------------------------------!
      cfg%io_output_interval_years = toml_int(t, 'io.output_interval_years', cfg%io_output_interval_years)
      cfg%io_cohort_max            = toml_int(t, 'io.cohort_max', cfg%io_cohort_max)
      cfg%io_patch_max             = toml_int(t, 'io.patch_max',  cfg%io_patch_max)

      !----- Per-PFT arrays (length cfg%pft%n; extra entries ignored, short arrays partial).-!
      npft = cfg%pft%n
      allocate(arr(npft))
      call toml_real_array(t, 'pft.wood_density', arr, nout)
      if (nout > 0_ik) then
         cfg%pft%wood_density(1:min(nout, npft)) = arr(1:min(nout, npft))
         call derive_pft_rates(cfg%pft)                 ! keep gr_max/mort_* consistent
      end if
      call toml_real_array(t, 'pft.dbh_critical', arr, nout)
      if (nout > 0_ik) cfg%pft%dbh_critical(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.recruit_dens', arr, nout)
      if (nout > 0_ik) cfg%pft%recruit_dens(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.recruit_min_temp', arr, nout)
      if (nout > 0_ik) cfg%pft%recruit_min_temp(1:min(nout, npft)) = arr(1:min(nout, npft))

      call validate_config(cfg)
   end subroutine load_meds_config

end module meds_config_io
