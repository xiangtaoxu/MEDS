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
   use meds_time,       only : meds_time_t, time_from_string
   use meds_pft_params, only : derive_pft_rates
   use meds_toml,       only : toml_table_t, toml_parse_file, toml_has, toml_int, toml_real,  &
                               toml_logical, toml_string, toml_real_array
   implicit none
   private

   public :: load_meds_config, write_pft_params_csv

contains

   !---------------------------------------------------------------------------------------!
   ! Write every per-PFT model parameter (one row per PFT) to a CSV at `path`. netCDF-free,  !
   ! so it works in every build. The two shared height thresholds are repeated on each row    !
   ! so the file is self-contained. A failed open warns but does not stop the run.            !
   !---------------------------------------------------------------------------------------!
   subroutine write_pft_params_csv(cfg, path)
      type(meds_config_t), intent(in) :: cfg
      character(len=*),    intent(in) :: path
      integer     :: u, ios
      integer(ik) :: pf

      open(newunit=u, file=path, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
         write(*,'(3a)') ' warning: could not write PFT parameters CSV "', trim(path), '"'
         return
      end if
      write(u,'(a)') 'pft,wood_density,dbh_critical,growth_dbh_slope,growth_dbh_cap,growth_dbh_max,'   &
           //'growth_lai_slope,reproduction_investment_fraction,repro_carbon_efficiency,'              &
           //'mort_gamma,mort_alpha,mort_beta,seed_rain_recruits,include_pft,'                         &
           //'min_cohort_height,min_reproduction_height'
      associate (p => cfg%pft)
         do pf = 1_ik, p%n
            write(u,'(i0,12(",",es15.8),",",i0,2(",",es15.8))')                                        &
                 pf, p%wood_density(pf), p%dbh_critical(pf), p%growth_dbh_slope(pf),                    &
                 p%growth_dbh_cap(pf), p%growth_dbh_max(pf), p%growth_lai_slope(pf),                    &
                 p%reproduction_investment_fraction(pf), p%repro_carbon_efficiency(pf),                &
                 p%mort_gamma(pf), p%mort_alpha(pf), p%mort_beta(pf), p%seed_rain_recruits(pf),         &
                 p%include_pft(pf), p%min_cohort_height, p%min_reproduction_height
         end do
      end associate
      close(u)
      write(*,'(2a)') ' pft   : ', trim(path)
   end subroutine write_pft_params_csv

   !---------------------------------------------------------------------------------------!
   ! Load configuration from `path`. If the file is absent, return the built-in defaults     !
   ! (found=.false.). The run is bounded by [run].start_time / [run].end_time (calendar       !
   ! dates "YYYY-MM-DD[ HH:MM:SS]"); an unparseable date keeps the built-in default.          !
   !---------------------------------------------------------------------------------------!
   subroutine load_meds_config(path, cfg, found)
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(out) :: cfg
      logical,             intent(out) :: found
      type(toml_table_t)    :: t
      character(len=256)    :: ts, sval
      integer(ik)           :: ts_mode, nout, npft
      real(wp), allocatable :: arr(:)
      type(meds_time_t)     :: tt
      logical               :: tok

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

      !----- Run span as calendar dates (default kept on a missing/unparseable value). -----!
      sval = toml_string(t, 'run.start_time', '')
      if (len_trim(sval) > 0) then
         call time_from_string(trim(sval), tt, tok)
         if (tok) cfg%start_time = tt
      end if
      sval = toml_string(t, 'run.end_time', '')
      if (len_trim(sval) > 0) then
         call time_from_string(trim(sval), tt, tok)
         if (tok) cfg%end_time = tt
      end if

      !----- Cohort/patch structural tunables. --------------------------------------------!
      cfg%demography_on          = toml_logical(t, 'demography.demography_on', cfg%demography_on)
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
      cfg%growth_memory_days     = toml_real(t, 'demography.growth_memory_days', cfg%growth_memory_days)

      !----- Disturbance. -----------------------------------------------------------------!
      cfg%patch_disturbance_rate     = toml_real(t, 'disturbance.patch_disturbance_rate',     cfg%patch_disturbance_rate)
      cfg%disturbance_survive_height = toml_real(t, 'disturbance.disturbance_survive_height', cfg%disturbance_survive_height)

      !----- Recruitment + the two height thresholds (scalar seed_rain_recruits applies to ALL!
      !       PFTs unless the [pft] array overrides it). ------------------------------------!
      cfg%min_recruit_size            = toml_real(t, 'recruitment.min_recruit_size', cfg%min_recruit_size)
      cfg%pft%min_cohort_height       = toml_real(t, 'recruitment.min_cohort_height',       &
                                                  cfg%pft%min_cohort_height)
      cfg%pft%min_reproduction_height = toml_real(t, 'recruitment.min_reproduction_height',   &
                                                  cfg%pft%min_reproduction_height)
      if (toml_has(t, 'recruitment.seed_rain_recruits'))                                      &
         cfg%pft%seed_rain_recruits = toml_real(t, 'recruitment.seed_rain_recruits',          &
                                                cfg%pft%seed_rain_recruits(1))

      !----- Initial conditions: init_mode selects the source; both files are read but only -!
      !       the one matching init_mode is used. ------------------------------------------!
      cfg%init_mode         = toml_int   (t, 'init.init_mode',    cfg%init_mode)
      cfg%init_restart_file = toml_string(t, 'init.restart_file', cfg%init_restart_file)
      cfg%init_census_file  = toml_string(t, 'init.census_file',  cfg%init_census_file)

      !----- netCDF output (path/prefix + what & how often to write). ---------------------!
      cfg%io_write_output          = toml_logical(t, 'io.write_output', cfg%io_write_output)
      cfg%io_output_dir            = toml_string(t, 'io.output_dir',    cfg%io_output_dir)
      cfg%io_output_prefix         = toml_string(t, 'io.output_prefix', cfg%io_output_prefix)
      cfg%io_output_interval_years = toml_int(t, 'io.output_interval_years', cfg%io_output_interval_years)
      cfg%io_cohort_max            = toml_int(t, 'io.cohort_max', cfg%io_cohort_max)
      cfg%io_patch_max             = toml_int(t, 'io.patch_max',  cfg%io_patch_max)
      cfg%io_write_state           = toml_logical(t, 'io.write_state', cfg%io_write_state)
      cfg%io_state_interval_years  = toml_int(t, 'io.state_interval_years', cfg%io_state_interval_years)

      !----- Per-PFT arrays (length cfg%pft%n; extra entries ignored, short arrays partial).-!
      npft = cfg%pft%n
      allocate(arr(npft))
      call toml_real_array(t, 'pft.wood_density', arr, nout)
      if (nout > 0_ik) then
         cfg%pft%wood_density(1:min(nout, npft)) = arr(1:min(nout, npft))
         call derive_pft_rates(cfg%pft)                 ! keep mort_gamma/alpha/beta consistent
      end if
      call toml_real_array(t, 'pft.dbh_critical', arr, nout)
      if (nout > 0_ik) cfg%pft%dbh_critical(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.seed_rain_recruits', arr, nout)
      if (nout > 0_ik) cfg%pft%seed_rain_recruits(1:min(nout, npft)) = arr(1:min(nout, npft))
      !----- Phenomenological growth / reproduction parameters (per PFT). ------------------!
      call toml_real_array(t, 'pft.growth_dbh_slope', arr, nout)
      if (nout > 0_ik) cfg%pft%growth_dbh_slope(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.growth_dbh_cap', arr, nout)
      if (nout > 0_ik) cfg%pft%growth_dbh_cap(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.growth_dbh_max', arr, nout)
      if (nout > 0_ik) cfg%pft%growth_dbh_max(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.growth_lai_slope', arr, nout)
      if (nout > 0_ik) cfg%pft%growth_lai_slope(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.reproduction_investment_fraction', arr, nout)
      if (nout > 0_ik) cfg%pft%reproduction_investment_fraction(1:min(nout, npft)) = arr(1:min(nout, npft))
      call toml_real_array(t, 'pft.repro_carbon_efficiency', arr, nout)
      if (nout > 0_ik) cfg%pft%repro_carbon_efficiency(1:min(nout, npft)) = arr(1:min(nout, npft))

      call validate_config(cfg)
   end subroutine load_meds_config

end module meds_config_io
