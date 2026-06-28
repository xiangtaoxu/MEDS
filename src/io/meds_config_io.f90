!==========================================================================================!
! meds_config_io -- load a run configuration from TWO mandatory TOML files and a presence-map !
! check. Lives in src/io/ (config input) with NO external dependency, so it is always compiled.!
!                                                                                          !
! Philosophy: NO hard-coded model parameters. Every parameter is supplied by the user via TOML !
! (or DERIVED from supplied ones); a missing parameter or a missing file is a HARD ERROR.      !
! The MAIN file (passed on the CLI) holds all non-PFT parameters and names the PFT file; the    !
! PFT file holds all PFT-specific traits plus the allometry and mortality-hazard coefficients.  !
! Presence is tracked while reading (a "missing key" list); after both files are read, any      !
! missing required key aborts with the full list. Derived quantities are then computed          !
! (derive_config + derive_pft_rates); with [options].override_derived a [derived] block in the  !
! PFT file overwrites the computed mortality-hazard parameters.                                 !
!==========================================================================================!
module meds_config_io
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, derive_config, validate_config,                  &
                               TS_DAILY, TS_WEEKLY, TS_MONTHLY, BK_SERIAL
   use meds_time,       only : meds_time_t, time_from_string
   use meds_allometry,  only : set_allometry
   use meds_pft_params, only : alloc_pft_table, derive_pft_rates
   use meds_toml,       only : toml_table_t, toml_parse_file, toml_has, toml_int, toml_real,  &
                               toml_logical, toml_string, toml_real_array
   implicit none
   private

   public :: load_meds_config, write_pft_params_csv

   integer(ik), parameter :: MAXPFT  = 64_ik    !< largest PFT count the array buffers handle
   integer,     parameter :: MAXMISS = 256       !< largest missing-key list

   !----- Accumulator of required keys that were absent (the presence map). ----------------!
   type :: keymiss_t
      integer(ik)       :: n = 0_ik
      character(len=72) :: key(MAXMISS) = ''
   end type keymiss_t

contains

   !=======================================================================================!
   !  Presence-map helpers: read a REQUIRED key, or record it as missing.                   !
   !=======================================================================================!
   subroutine note_missing(m, key)
      type(keymiss_t),  intent(inout) :: m
      character(len=*), intent(in)    :: key
      if (m%n < MAXMISS) then ; m%n = m%n + 1_ik ; m%key(m%n) = key ; end if
   end subroutine note_missing

   subroutine req_r(t, key, out, m)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      real(wp),           intent(out) :: out
      type(keymiss_t),    intent(inout) :: m
      out = 0.0_wp
      if (toml_has(t, key)) then ; out = toml_real(t, key, 0.0_wp) ; else ; call note_missing(m, key) ; end if
   end subroutine req_r

   subroutine req_i(t, key, out, m)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      integer(ik),        intent(out) :: out
      type(keymiss_t),    intent(inout) :: m
      out = 0_ik
      if (toml_has(t, key)) then ; out = toml_int(t, key, 0_ik) ; else ; call note_missing(m, key) ; end if
   end subroutine req_i

   subroutine req_l(t, key, out, m)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      logical,            intent(out) :: out
      type(keymiss_t),    intent(inout) :: m
      out = .false.
      if (toml_has(t, key)) then ; out = toml_logical(t, key, .false.) ; else ; call note_missing(m, key) ; end if
   end subroutine req_l

   subroutine req_s(t, key, out, m)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      character(len=*),   intent(out) :: out
      type(keymiss_t),    intent(inout) :: m
      out = ''
      if (toml_has(t, key)) then ; out = toml_string(t, key, '') ; else ; call note_missing(m, key) ; end if
   end subroutine req_s

   subroutine req_ts(t, key, mode, m)        ! timestep string -> TS_* mode
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      integer(ik),        intent(out) :: mode
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: ts
      mode = TS_DAILY
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      ts = toml_string(t, key, 'daily')
      select case (trim(ts))
      case ('monthly') ; mode = TS_MONTHLY
      case ('weekly')  ; mode = TS_WEEKLY
      case default     ; mode = TS_DAILY
      end select
   end subroutine req_ts

   subroutine req_date(t, key, out, m)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      type(meds_time_t),  intent(out) :: out
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: sval
      logical :: ok
      out = meds_time_t()
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      sval = toml_string(t, key, '')
      call time_from_string(trim(sval), out, ok)
      if (.not. ok) call note_missing(m, key)         ! present but unparseable
   end subroutine req_date

   !----- Required per-PFT real array of length npft (record missing / length-mismatch). ---!
   subroutine req_pa(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      real(wp)    :: buf(MAXPFT)
      integer(ik) :: nout
      call toml_real_array(t, key, buf, nout)
      if (nout == npft) then ; out(1:npft) = buf(1:npft) ; else ; call note_missing(m, key) ; end if
   end subroutine req_pa

   subroutine req_pa_int(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      integer(ik),        intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      real(wp)    :: buf(MAXPFT)
      integer(ik) :: nout
      call toml_real_array(t, key, buf, nout)
      if (nout == npft) then ; out(1:npft) = nint(buf(1:npft), ik) ; else ; call note_missing(m, key) ; end if
   end subroutine req_pa_int

   !=======================================================================================!
   !  Load configuration from `path` (the MAIN file). Both files are mandatory; a missing    !
   !  file or any missing required key is a hard error.                                       !
   !=======================================================================================!
   subroutine load_meds_config(path, cfg)
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(out) :: cfg
      type(toml_table_t) :: tm, tp
      type(keymiss_t)    :: miss
      logical            :: found
      integer(ik)        :: npft, nout, i
      real(wp)           :: buf(MAXPFT)
      real(wp) :: b1Ht, b2Ht, height_max, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2, light_ext

      !----- MAIN file. -------------------------------------------------------------------!
      call toml_parse_file(path, tm, found)
      if (.not. found) error stop 'meds_config: main config file not found: '//trim(path)
      cfg%backend = BK_SERIAL                        ! reporting only, not a config parameter

      call req_ts  (tm, 'run.timestep',   cfg%ts_mode,    miss)
      call req_date(tm, 'run.start_time', cfg%start_time, miss)
      call req_date(tm, 'run.end_time',   cfg%end_time,   miss)

      call req_l(tm, 'demography.demography_on',          cfg%demography_on,          miss)
      call req_l(tm, 'demography.do_cohort_fissfuse',     cfg%do_cohort_fissfuse,     miss)
      call req_l(tm, 'demography.do_patch_fissfuse',      cfg%do_patch_fissfuse,      miss)
      call req_l(tm, 'demography.do_patch_disturbance',   cfg%do_patch_disturbance,   miss)
      call req_i(tm, 'demography.max_cohort',             cfg%max_cohort,             miss)
      call req_i(tm, 'demography.n_cohort_fusion_iter',   cfg%n_cohort_fusion_iter,   miss)
      call req_r(tm, 'demography.cohort_size_tol_min',    cfg%cohort_size_tol_min,    miss)
      call req_r(tm, 'demography.cohort_size_tol_max',    cfg%cohort_size_tol_max,    miss)
      call req_r(tm, 'demography.cohort_lai_cap',         cfg%cohort_lai_cap,         miss)
      call req_r(tm, 'demography.min_cohort_agb',         cfg%min_cohort_agb,         miss)
      call req_r(tm, 'demography.negligible_nplant',      cfg%negligible_nplant,      miss)
      call req_r(tm, 'demography.split_eps',              cfg%split_eps,              miss)
      call req_l(tm, 'demography.enable_cohort_fission',  cfg%enable_cohort_fission,  miss)
      call req_i(tm, 'demography.n_height_layers',        cfg%n_height_layers,        miss)
      call req_i(tm, 'demography.max_patch',              cfg%max_patch,              miss)
      call req_i(tm, 'demography.n_patch_fusion_iter',    cfg%n_patch_fusion_iter,    miss)
      call req_r(tm, 'demography.patch_light_tol',        cfg%patch_light_tol,        miss)
      call req_r(tm, 'demography.patch_light_maxdev_factor', cfg%patch_light_maxdev_factor, miss)
      call req_r(tm, 'demography.patch_diff_age_tol',     cfg%patch_diff_age_tol,     miss)
      call req_r(tm, 'demography.min_patch_area',         cfg%min_patch_area,         miss)
      call req_r(tm, 'demography.patch_min_area_remain',  cfg%patch_min_area_remain,  miss)
      call req_l(tm, 'demography.enable_patch_fission',   cfg%enable_patch_fission,   miss)
      call req_r(tm, 'demography.conservation_tol',       cfg%conservation_tol,       miss)
      call req_r(tm, 'demography.growth_memory_days',     cfg%growth_memory_days,     miss)

      call req_r(tm, 'disturbance.patch_disturbance_rate',     cfg%patch_disturbance_rate,     miss)
      call req_r(tm, 'disturbance.disturbance_survive_height', cfg%disturbance_survive_height, miss)

      call req_r(tm, 'recruitment.min_recruit_size',      cfg%min_recruit_size,       miss)

      call req_i(tm, 'init.init_mode',     cfg%init_mode,         miss)
      call req_s(tm, 'init.restart_file',  cfg%init_restart_file, miss)
      call req_s(tm, 'init.census_file',   cfg%init_census_file,  miss)
      call req_s(tm, 'init.pft_config',    cfg%pft_config,        miss)

      call req_l(tm, 'io.write_output',          cfg%io_write_output,          miss)
      call req_s(tm, 'io.output_dir',            cfg%io_output_dir,            miss)
      call req_s(tm, 'io.output_prefix',         cfg%io_output_prefix,         miss)
      call req_i(tm, 'io.output_interval_years', cfg%io_output_interval_years, miss)
      call req_i(tm, 'io.cohort_max',            cfg%io_cohort_max,            miss)
      call req_i(tm, 'io.patch_max',             cfg%io_patch_max,             miss)
      call req_l(tm, 'io.write_state',           cfg%io_write_state,           miss)
      call req_i(tm, 'io.state_interval_years',  cfg%io_state_interval_years,  miss)

      call req_l(tm, 'options.override_derived', cfg%override_derived,         miss)

      !----- PFT file (named in the main file). -------------------------------------------!
      call toml_parse_file(trim(cfg%pft_config), tp, found)
      if (.not. found) error stop 'meds_config: PFT config file not found: '//trim(cfg%pft_config)

      !----- PFT count comes from the length of the wood_density array. -------------------!
      call toml_real_array(tp, 'pft.wood_density', buf, nout)
      if (nout < 1_ik) error stop 'meds_config: pft.wood_density missing/empty in '//trim(cfg%pft_config)
      npft = nout
      call alloc_pft_table(cfg%pft, npft)
      cfg%pft%wood_density = buf(1:npft)

      call req_pa(tp, 'pft.dbh_critical',                     cfg%pft%dbh_critical,                     npft, miss)
      call req_pa(tp, 'pft.growth_dbh_slope',                 cfg%pft%growth_dbh_slope,                 npft, miss)
      call req_pa(tp, 'pft.growth_dbh_cap',                   cfg%pft%growth_dbh_cap,                   npft, miss)
      call req_pa(tp, 'pft.growth_dbh_max',                   cfg%pft%growth_dbh_max,                   npft, miss)
      call req_pa(tp, 'pft.growth_lai_slope',                 cfg%pft%growth_lai_slope,                 npft, miss)
      call req_pa(tp, 'pft.reproduction_investment_fraction', cfg%pft%reproduction_investment_fraction, npft, miss)
      call req_pa(tp, 'pft.repro_carbon_efficiency',          cfg%pft%repro_carbon_efficiency,          npft, miss)
      call req_pa(tp, 'pft.seed_rain_recruits',               cfg%pft%seed_rain_recruits,               npft, miss)
      call req_pa_int(tp, 'pft.include_pft',                  cfg%pft%include_pft,                      npft, miss)
      call req_r(tp, 'pft.min_cohort_height',       cfg%pft%min_cohort_height,       miss)
      call req_r(tp, 'pft.min_reproduction_height', cfg%pft%min_reproduction_height, miss)

      call req_r(tp, 'camac.mort_rho_ref',   cfg%pft%mort_rho_ref,   miss)
      call req_r(tp, 'camac.mort_gamma_0',   cfg%pft%mort_gamma_0,   miss)
      call req_r(tp, 'camac.mort_gamma_exp', cfg%pft%mort_gamma_exp, miss)
      call req_r(tp, 'camac.mort_alpha_0',   cfg%pft%mort_alpha_0,   miss)
      call req_r(tp, 'camac.mort_alpha_exp', cfg%pft%mort_alpha_exp, miss)
      call req_r(tp, 'camac.mort_beta_0',    cfg%pft%mort_beta_0,    miss)
      call req_r(tp, 'camac.mort_beta_exp',  cfg%pft%mort_beta_exp,  miss)

      call req_r(tp, 'allometry.b1Ht',       b1Ht,       miss)
      call req_r(tp, 'allometry.b2Ht',       b2Ht,       miss)
      call req_r(tp, 'allometry.height_max', height_max, miss)
      call req_r(tp, 'allometry.agb_c1',     agb_c1,     miss)
      call req_r(tp, 'allometry.agb_c2',     agb_c2,     miss)
      call req_r(tp, 'allometry.ca_b1',      ca_b1,      miss)
      call req_r(tp, 'allometry.ca_b2',      ca_b2,      miss)
      call req_r(tp, 'allometry.lai_b1',     lai_b1,     miss)
      call req_r(tp, 'allometry.lai_b2',     lai_b2,     miss)
      call req_r(tp, 'allometry.light_ext',  light_ext,  miss)

      !----- Presence check: abort with the full list of any missing required keys. -------!
      if (miss%n > 0_ik) then
         write(*,'(a,i0,a)') ' meds_config: ', miss%n, ' required parameter(s) missing:'
         do i = 1_ik, miss%n
            write(*,'(3a)') '   - ', trim(miss%key(i)), ''
         end do
         error stop 'meds_config: incomplete configuration (see missing keys above)'
      end if

      !----- Install allometry, then compute every derived quantity. ----------------------!
      call set_allometry(b1Ht, b2Ht, height_max, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2, light_ext)
      call derive_config(cfg)
      call derive_pft_rates(cfg%pft)

      !----- Global override: a [derived] block in the PFT file pins the mortality params. -!
      if (cfg%override_derived) then
         call toml_real_array(tp, 'derived.mort_gamma', buf, nout)
         if (nout == npft) cfg%pft%mort_gamma = buf(1:npft)
         call toml_real_array(tp, 'derived.mort_alpha', buf, nout)
         if (nout == npft) cfg%pft%mort_alpha = buf(1:npft)
         call toml_real_array(tp, 'derived.mort_beta', buf, nout)
         if (nout == npft) cfg%pft%mort_beta = buf(1:npft)
      end if

      call validate_config(cfg)
   end subroutine load_meds_config

   !---------------------------------------------------------------------------------------!
   ! Write every per-PFT model parameter (one row per PFT) to a CSV at `path`. netCDF-free,  !
   ! a provenance record of what the run actually used. A failed open warns, does not stop.   !
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

end module meds_config_io
