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
   use meds_config,     only : meds_config_t, derive_parameters, validate_config,               &
                               BK_SERIAL,                                                       &
                               SM_LEUNING, SM_MEDLYN, SM_KATUL,                                 &
                               TRESP_ARRHENIUS, TRESP_PEAKED, COLIM_MIN, COLIM_QUADRATIC,        &
                               SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED, INTEG_SPLIT, INTEG_ARK
   use meds_forcing_config, only : forcing_config_t,                                            &
                                   MET_BACKEND_CONST, MET_BACKEND_NETCDF,                       &
                                   METAVG_INSTANT, METAVG_END, METAVG_BEGIN, METAVG_CENTER,      &
                                   SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN,      &
                                   LW_FILE, LW_SYNTHESIZE, CLAMP_ERROR, CLAMP_HOLD,              &
                                   GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   use meds_output_config, only : output_config_t, GRP_STRUCTURE, GRP_CARBON, GRP_WATER,          &
                                   GRP_ENERGY, SYNC_FLUSH, SYNC_NEVER, FC_DAY, FC_MONTH, FC_YEAR, &
                                   FC_RUN
   use meds_time,       only : meds_time_t, time_from_string
   use meds_pft_params, only : alloc_pft_table
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

   subroutine req_dur(t, key, secs, m)       ! duration string ("1d","7d","12h","1800s") -> seconds
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(out)   :: secs
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: s
      character         :: u
      integer(ik)       :: n
      integer           :: ios
      real(wp)          :: v, mult
      secs = 86400.0_wp                        ! 1 day (used only if the key is missing)
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = adjustl(toml_string(t, key, '1d'))
      n = int(len_trim(s), ik)
      if (n == 0_ik) then ; call note_missing(m, key) ; return ; end if   ! present but empty -> invalid
      u = s(n:n)
      select case (u)
      case ('s', 'S') ; mult = 1.0_wp
      case ('m', 'M') ; mult = 60.0_wp
      case ('h', 'H') ; mult = 3600.0_wp
      case ('d', 'D') ; mult = 86400.0_wp
      case ('w', 'W') ; mult = 604800.0_wp
      case default    ; mult = 0.0_wp          ! no unit suffix -> the whole string is seconds
      end select
      if (mult == 0.0_wp) then
         read(s(1:n), *, iostat=ios) v ; mult = 1.0_wp
      else
         read(s(1:n-1_ik), *, iostat=ios) v
      end if
      if (ios == 0) then
         secs = v * mult
      else
         call note_missing(m, key)          ! present but unparseable -> aggregated required-key error
      end if
   end subroutine req_dur

   subroutine req_stomatal_model(t, key, mode, m)   ! stomatal-model string -> SM_* mode
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      integer(ik),        intent(out)   :: mode
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: s
      mode = SM_MEDLYN
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'medlyn')
      select case (trim(s))
      case ('leuning') ; mode = SM_LEUNING
      case ('medlyn')  ; mode = SM_MEDLYN
      case ('katul')   ; mode = SM_KATUL
      case default     ; call note_missing(m, key)   ! present but unrecognized -> hard error
      end select
   end subroutine req_stomatal_model

   subroutine req_scheme(t, key, mode, m)           ! fast-loop coupling scheme string -> SCHEME_* mode
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      integer(ik),        intent(out)   :: mode
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: s
      mode = SCHEME_SPLIT_SEQUENTIAL
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'split')
      select case (trim(s))
      case ('split')  ; mode = SCHEME_SPLIT_SEQUENTIAL
      case ('picard') ; mode = SCHEME_PICARD_COUPLED
      case default    ; call note_missing(m, key)      ! present but unrecognized -> hard error
      end select
   end subroutine req_scheme

   subroutine req_temp_response(t, key, mode, m)    ! temperature-response string -> TRESP_* mode
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      integer(ik),        intent(out)   :: mode
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: s
      mode = TRESP_PEAKED
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'peaked')
      select case (trim(s))
      case ('arrhenius') ; mode = TRESP_ARRHENIUS
      case ('peaked')    ; mode = TRESP_PEAKED
      case default       ; call note_missing(m, key)   ! present but unrecognized -> hard error
      end select
   end subroutine req_temp_response

   subroutine req_colimitation(t, key, mode, m)     ! co-limitation string -> COLIM_* mode
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      integer(ik),        intent(out)   :: mode
      type(keymiss_t),    intent(inout) :: m
      character(len=64) :: s
      mode = COLIM_QUADRATIC
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'quadratic')
      select case (trim(s))
      case ('min')       ; mode = COLIM_MIN
      case ('quadratic') ; mode = COLIM_QUADRATIC
      case default       ; call note_missing(m, key)   ! present but unrecognized -> hard error
      end select
   end subroutine req_colimitation

   !----- [forcing] string-enum mappers (the req_scheme pattern). --------------------------!
   subroutine req_met_backend(t, key, mode, m)      ! "netcdf" | "const"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = MET_BACKEND_NETCDF
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'netcdf')
      select case (trim(s))
      case ('netcdf') ; mode = MET_BACKEND_NETCDF
      case ('const')  ; mode = MET_BACKEND_CONST
      case default    ; call note_missing(m, key)
      end select
   end subroutine req_met_backend

   subroutine req_avg_convention(t, key, mode, m)   ! "instant" | "end" | "begin" | "center"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = METAVG_END
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'end')
      select case (trim(s))
      case ('instant') ; mode = METAVG_INSTANT
      case ('end')     ; mode = METAVG_END
      case ('begin')   ; mode = METAVG_BEGIN
      case ('center')  ; mode = METAVG_CENTER
      case default     ; call note_missing(m, key)
      end select
   end subroutine req_avg_convention

   subroutine req_sw_partition(t, key, mode, m)     ! "passthrough" | "clearidx" | "weiss_norman"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = SWPART_CLEARIDX
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'clearidx')
      select case (trim(s))
      case ('passthrough')  ; mode = SWPART_PASSTHROUGH
      case ('clearidx')     ; mode = SWPART_CLEARIDX
      case ('weiss_norman') ; mode = SWPART_WEISS_NORMAN
      case default          ; call note_missing(m, key)
      end select
   end subroutine req_sw_partition

   subroutine req_lwdown_source(t, key, mode, m)    ! "file" | "synthesize"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = LW_FILE
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'file')
      select case (trim(s))
      case ('file')       ; mode = LW_FILE
      case ('synthesize') ; mode = LW_SYNTHESIZE
      case default        ; call note_missing(m, key)
      end select
   end subroutine req_lwdown_source

   subroutine req_start_clamp(t, key, mode, m)      ! "error" | "hold"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = CLAMP_ERROR
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'error')
      select case (trim(s))
      case ('error') ; mode = CLAMP_ERROR
      case ('hold')  ; mode = CLAMP_HOLD
      case default   ; call note_missing(m, key)
      end select
   end subroutine req_start_clamp

   subroutine req_grid_match(t, key, mode, m)       ! "explicit" | "nearest"
      type(toml_table_t), intent(in) :: t ; character(len=*), intent(in) :: key
      integer(ik), intent(out) :: mode ; type(keymiss_t), intent(inout) :: m
      character(len=64) :: s
      mode = GRIDMATCH_EXPLICIT
      if (.not. toml_has(t, key)) then ; call note_missing(m, key) ; return ; end if
      s = toml_string(t, key, 'explicit')
      select case (trim(s))
      case ('explicit') ; mode = GRIDMATCH_EXPLICIT
      case ('nearest')  ; mode = GRIDMATCH_NEAREST
      case default      ; call note_missing(m, key)
      end select
   end subroutine req_grid_match

   !----- Load the [forcing]/[site] block. GATED on forcing.forcing_on (a DEFAULTED read, so a  !
   !      config with no [forcing] block runs the constant-forcing MVP unchanged). When ON, the   !
   !      rest are REQUIRED (the no-silent-defaults rule) -- design MEDS_FORCING_DESIGN.md §6.6.   !
   subroutine load_forcing_config(t, cfg, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      type(keymiss_t),     intent(inout) :: m
      cfg%forcing = forcing_config_t()                                  ! Ithaca/ERA5-Land defaults
      cfg%forcing%forcing_on = toml_logical(t, 'forcing.forcing_on', .false.)   ! opt-in gate (defaulted)
      if (.not. cfg%forcing%forcing_on) return
      call req_s            (t, 'forcing.path',           cfg%forcing%path,                  m)
      call req_i            (t, 'forcing.grid_index',     cfg%forcing%grid_index,            m)
      call req_grid_match   (t, 'forcing.grid_match',     cfg%forcing%grid_match,            m)
      call req_dur          (t, 'forcing.timestep',       cfg%forcing%dt_forcing,            m)
      call req_met_backend  (t, 'forcing.format',         cfg%forcing%backend,               m)
      call req_avg_convention(t, 'forcing.avg_convention', cfg%forcing%avg_convention,       m)
      call req_sw_partition (t, 'forcing.sw_partition',   cfg%forcing%sw_partition,          m)
      call req_lwdown_source(t, 'forcing.lwdown_source',  cfg%forcing%lwdown_source,         m)
      call req_start_clamp  (t, 'forcing.start_clamp',    cfg%forcing%start_clamp,           m)
      call req_l            (t, 'forcing.recycle',        cfg%forcing%recycle,               m)
      call req_r            (t, 'forcing.co2_const',      cfg%forcing%co2_const,             m)
      call req_r            (t, 'site.latitude',          cfg%forcing%latitude_deg,          m)
      call req_r            (t, 'site.longitude',         cfg%forcing%longitude_deg,         m)
      call req_r            (t, 'site.utc_offset',        cfg%forcing%utc_offset_h,          m)
      call req_l            (t, 'site.apply_solar_longitude', cfg%forcing%apply_solar_longitude, m)
      call req_r            (t, 'site.reference_height',  cfg%forcing%reference_height,      m)
      call req_r            (t, 'site.wind_meas_height',  cfg%forcing%wind_meas_height,      m)
      call req_r            (t, 'site.elevation',         cfg%forcing%elevation_m,           m)
      !----- wind-height + elevation-lapse corrections (§5.2/§10-Q2). -------------------------!
      call req_l            (t, 'site.apply_wind_profile',    cfg%forcing%apply_wind_profile,    m)
      call req_l            (t, 'site.apply_elevation_lapse', cfg%forcing%apply_elevation_lapse, m)
      call req_r            (t, 'site.wind_roughness_z0',     cfg%forcing%wind_roughness_z0,     m)
      call req_r            (t, 'site.lapse_rate_tair',       cfg%forcing%lapse_rate_tair,       m)
      call req_r            (t, 'site.grid_elevation',        cfg%forcing%grid_elevation_m,      m)
   end subroutine load_forcing_config

   !----- Load the [output] diagnostic-aggregation block. OPT-IN: gated on output.enabled (a       !
   !      DEFAULTED read, so a config with no [output] block leaves the new path OFF and the legacy  !
   !      [io] path runs unchanged). ALL keys are optional-with-default (§6.1 softening), so no       !
   !      keymiss is recorded. Per-variable overrides live in the optional meds_io_config.toml,        !
   !      applied by the driver against the built registry (§6, MEDS_IO_DESIGN.md).                    !
   subroutine load_output_config(t, cfg)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      cfg%output = output_config_t()                                     ! defaults
      cfg%output%enabled = toml_logical(t, 'output.enabled', .false.)    ! opt-in gate (defaulted)
      if (.not. cfg%output%enabled) return
      cfg%output%dir                 = toml_string (t, 'output.dir',                 '.')
      cfg%output%prefix              = toml_string (t, 'output.prefix',              'meds')
      cfg%output%io_config           = toml_string (t, 'output.io_config',           '')
      cfg%output%cohort_max          = toml_int    (t, 'output.cohort_max',          4096_ik)
      cfg%output%patch_max           = toml_int    (t, 'output.patch_max',           256_ik)
      cfg%output%strict_caps         = toml_logical(t, 'output.strict_caps',         .false.)
      cfg%output%fast_interval_steps = toml_int    (t, 'output.fast_interval_steps', 4_ik)
      cfg%output%sync_every = merge(SYNC_NEVER, SYNC_FLUSH,                                        &
                                    trim(toml_string(t, 'output.sync_every', 'flush')) == 'never')
      !----- high-level flux-group toggles (§6). ---!
      cfg%output%grp_on(GRP_STRUCTURE) = toml_logical(t, 'output.structure',     .true.)
      cfg%output%grp_on(GRP_CARBON)    = toml_logical(t, 'output.carbon_fluxes', .true.)
      cfg%output%grp_on(GRP_WATER)     = toml_logical(t, 'output.water_fluxes',  .false.)
      cfg%output%grp_on(GRP_ENERGY)    = toml_logical(t, 'output.energy_fluxes', .false.)
      !----- per-tier enable + file-chunk (index order FAST/DAILY/MONTHLY/ANNUAL). ---!
      cfg%output%freq_on(1) = toml_logical(t, 'output.fast.enabled',    .false.)
      cfg%output%freq_on(2) = toml_logical(t, 'output.daily.enabled',   .true.)
      cfg%output%freq_on(3) = toml_logical(t, 'output.monthly.enabled', .true.)
      cfg%output%freq_on(4) = toml_logical(t, 'output.annual.enabled',  .true.)
      cfg%output%file_chunk(1) = parse_file_chunk(toml_string(t, 'output.fast.file_chunk',    'day'))
      cfg%output%file_chunk(2) = parse_file_chunk(toml_string(t, 'output.daily.file_chunk',   'month'))
      cfg%output%file_chunk(3) = parse_file_chunk(toml_string(t, 'output.monthly.file_chunk', 'year'))
      cfg%output%file_chunk(4) = parse_file_chunk(toml_string(t, 'output.annual.file_chunk',  'run'))
   end subroutine load_output_config

   !----- Load the per-PFT [phenology] cue params from the PFT file. GATED on cfg%phenology_on (read !
   !      from the MAIN file just above): OFF -> keep the alloc_pft_table literature defaults (nothing  !
   !      required). ON -> every per-PFT array is REQUIRED (the no-silent-defaults rule). The arrays     !
   !      are flattened per cohort into a meds_plant pheno_params_t by the slow-loop phenology driver.   !
   subroutine load_phenology_pft(t, cfg, npft, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      integer(ik),         intent(in)    :: npft
      type(keymiss_t),     intent(inout) :: m
      if (.not. cfg%phenology_on) return
      call req_pa_int(t, 'phenology.cue_mask',        cfg%pft%pheno_cue_mask,            npft, m)
      call req_pa(t, 'phenology.cue_sharpness',       cfg%pft%pheno_cue_sharpness,       npft, m)
      call req_pa(t, 'phenology.on_threshold',        cfg%pft%pheno_on_threshold,        npft, m)
      call req_pa(t, 'phenology.off_threshold',       cfg%pft%pheno_off_threshold,       npft, m)
      call req_pa(t, 'phenology.gdd_base_temp',       cfg%pft%pheno_gdd_base_temp,       npft, m)
      call req_pa(t, 'phenology.chill_base_temp',     cfg%pft%pheno_chill_base_temp,     npft, m)
      call req_pa(t, 'phenology.phen_a',              cfg%pft%pheno_phen_a,              npft, m)
      call req_pa(t, 'phenology.phen_b',              cfg%pft%pheno_phen_b,              npft, m)
      call req_pa(t, 'phenology.phen_c',              cfg%pft%pheno_phen_c,              npft, m)
      call req_pa(t, 'phenology.cold_drop_daylength', cfg%pft%pheno_cold_drop_daylength, npft, m)
      call req_pa(t, 'phenology.cold_drop_soiltemp1', cfg%pft%pheno_cold_drop_soiltemp1, npft, m)
      call req_pa(t, 'phenology.cold_drop_soiltemp2', cfg%pft%pheno_cold_drop_soiltemp2, npft, m)
      call req_pa(t, 'phenology.photo_crit',          cfg%pft%pheno_photo_crit,          npft, m)
      call req_pa(t, 'phenology.photo_slope',         cfg%pft%pheno_photo_slope,         npft, m)
   end subroutine load_phenology_pft

   !----- "day" | "month" | "year" | "run" -> FC_* (unknown -> error stop naming the offender). ---!
   pure integer(ik) function parse_file_chunk(s) result(fc)
      character(len=*), intent(in) :: s
      select case (trim(adjustl(s)))
      case ('day')   ; fc = FC_DAY
      case ('month') ; fc = FC_MONTH
      case ('year')  ; fc = FC_YEAR
      case ('run')   ; fc = FC_RUN
      case default   ; error stop 'meds_config_io: unknown output file_chunk (expected day|month|year|run)'
      end select
   end function parse_file_chunk

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

      !----- MAIN file. -------------------------------------------------------------------!
      call toml_parse_file(path, tm, found)
      if (.not. found) error stop 'meds_config: main config file not found: '//trim(path)
      cfg%backend = BK_SERIAL                        ! reporting only, not a config parameter

      call req_dur (tm, 'run.dt_slow',    cfg%dt_slow,    miss)
      call req_date(tm, 'run.start_time', cfg%start_time, miss)
      call req_date(tm, 'run.end_time',   cfg%end_time,   miss)

      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      call req_l     (tm, 'fast.fast_biophysics_on', cfg%fast_biophysics_on, miss)
      call req_dur   (tm, 'fast.dt_fast',            cfg%dt_fast,            miss)
      call req_scheme(tm, 'fast.integration_scheme', cfg%integration_scheme, miss)
      !----- P3 coupled-surface (Picard) knobs + option selectors: DEFAULTED reads (solver tuning, !
      !      not physical params), so a config without them runs the documented defaults. ---------!
      cfg%picard_max_iter     = toml_int    (tm, 'fast.picard_max_iter',    20_ik)
      cfg%picard_tol_temp     = toml_real   (tm, 'fast.picard_tol_temp',    1.0e-3_wp)
      cfg%picard_tol_shv      = toml_real   (tm, 'fast.picard_tol_shv',     1.0e-6_wp)
      cfg%picard_relax        = toml_real   (tm, 'fast.picard_relax',       0.5_wp)
      cfg%picard_fixed_iter   = toml_logical(tm, 'fast.picard_fixed_iter',  .false.)
      cfg%leaf_energy_model   = merge(1_ik, 0_ik,                                                  &
                                trim(toml_string(tm, 'fast.leaf_energy_model',   'diagnostic')) == 'prognostic')
      cfg%wood_energy_model   = merge(1_ik, 0_ik,                                                  &
                                trim(toml_string(tm, 'fast.wood_energy_model',   'diagnostic')) == 'prognostic')
      cfg%soil_water_coupling = merge(1_ik, 0_ik,                                                  &
                                trim(toml_string(tm, 'fast.soil_water_coupling', 'lagged')) == 'coupled')
      cfg%snow_on             = toml_logical(tm, 'fast.snow_on', .false.)
      cfg%snow_init_swe       = toml_real(tm, 'fast.snow_init_swe',  0.0_wp)
      cfg%snow_init_temp      = toml_real(tm, 'fast.snow_init_temp', 270.0_wp)
      !----- Fast-loop TIME integrator selector + ARK knobs: DEFAULTED reads (absent -> INTEG_SPLIT,   !
      !      so every existing config + the golden anchor stay byte-identical). ------------------------!
      cfg%time_integrator   = merge(INTEG_ARK, INTEG_SPLIT,                                          &
                              trim(toml_string(tm, 'fast.time_integrator', 'split')) == 'ark')
      cfg%ark_adaptive      = toml_logical(tm, 'fast.ark_adaptive',      .true.)
      cfg%ark_rtol          = toml_real   (tm, 'fast.ark_rtol',          1.0e-3_wp)
      cfg%ark_dt_init       = toml_real   (tm, 'fast.ark_dt_init',       0.0_wp)
      cfg%ark_fixed_substep = toml_int    (tm, 'fast.ark_fixed_substep', 4_ik)
      cfg%ark_niter         = toml_int    (tm, 'fast.ark_niter',         8_ik)
      cfg%ark_relax         = toml_real   (tm, 'fast.ark_relax',         0.6_wp)
      cfg%fast_probe        = toml_logical(tm, 'fast.fast_probe',        .false.)
      cfg%fast_probe_file   = toml_string (tm, 'fast.fast_probe_file',   'fast_probe.csv')

      !----- [hydraulics] (opt-in; each key DEFAULTED to the placeholder in hydraulics_config_t, so an  !
      !      absent block is byte-identical to the former hardcoded fast-loop stub). Consumed only by    !
      !      the fast loop (fast_biophysics_on); PFT-uniform MVP. -----------------------------------!
      cfg%hydraulics%leaf_pi0       = toml_real(tm, 'hydraulics.leaf_pi0',       cfg%hydraulics%leaf_pi0)
      cfg%hydraulics%leaf_elastic_mod       = toml_real(tm, 'hydraulics.leaf_elastic_mod',       cfg%hydraulics%leaf_elastic_mod)
      cfg%hydraulics%leaf_apoplast_frac        = toml_real(tm, 'hydraulics.leaf_apoplast_frac',        cfg%hydraulics%leaf_apoplast_frac)
      cfg%hydraulics%leaf_water_sat = toml_real(tm, 'hydraulics.leaf_water_sat', cfg%hydraulics%leaf_water_sat)
      cfg%hydraulics%wood_pi0       = toml_real(tm, 'hydraulics.wood_pi0',       cfg%hydraulics%wood_pi0)
      cfg%hydraulics%wood_elastic_mod       = toml_real(tm, 'hydraulics.wood_elastic_mod',       cfg%hydraulics%wood_elastic_mod)
      cfg%hydraulics%wood_apoplast_frac        = toml_real(tm, 'hydraulics.wood_apoplast_frac',        cfg%hydraulics%wood_apoplast_frac)
      cfg%hydraulics%wood_water_sat = toml_real(tm, 'hydraulics.wood_water_sat', cfg%hydraulics%wood_water_sat)
      cfg%hydraulics%wood_psi50     = toml_real(tm, 'hydraulics.wood_psi50',     cfg%hydraulics%wood_psi50)
      cfg%hydraulics%wood_kexp      = toml_real(tm, 'hydraulics.wood_kexp',      cfg%hydraulics%wood_kexp)
      cfg%hydraulics%k_plant_max    = toml_real(tm, 'hydraulics.k_plant_max',    cfg%hydraulics%k_plant_max)
      cfg%hydraulics%wood_kmax      = toml_real(tm, 'hydraulics.wood_kmax',      cfg%hydraulics%wood_kmax)
      cfg%hydraulics%vessel_curl    = toml_real(tm, 'hydraulics.vessel_curl',    cfg%hydraulics%vessel_curl)
      cfg%hydraulics%rhizo_cond     = toml_real(tm, 'hydraulics.rhizo_cond',     cfg%hydraulics%rhizo_cond)
      cfg%hydraulics%root_beta          = toml_real(tm, 'hydraulics.root_beta',          cfg%hydraulics%root_beta)
      cfg%hydraulics%root_depth         = toml_real(tm, 'hydraulics.root_depth',         cfg%hydraulics%root_depth)
      cfg%hydraulics%specific_root_area = toml_real(tm, 'hydraulics.specific_root_area', cfg%hydraulics%specific_root_area)
      cfg%hydraulics%multilayer_roots   = toml_logical(tm, 'hydraulics.multilayer_roots', cfg%hydraulics%multilayer_roots)

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

      !----- Carbon growth: stub GPP (used when the fast loop is off). ---------------------!
      call req_r            (tm, 'carbon.gpp_ref',       cfg%gpp_ref,       miss)

      !----- Leaf phenology (opt-in; gated on phenology.phenology_on, a DEFAULTED read, so a config !
      !      with no [phenology] block leaves phenology OFF and the leaf-flush gate hard-wired ON). --!
      cfg%phenology_on = toml_logical(tm, 'phenology.phenology_on', .false.)

      !----- Meteorological forcing (opt-in; gated on forcing.forcing_on, defaulted false). !
      call load_forcing_config(tm, cfg, miss)

      !----- Diagnostic-aggregation output (opt-in; gated on output.enabled, defaulted false; !
      !      a config with no [output] block runs the legacy [io] path unchanged).  ----------!
      call load_output_config(tm, cfg)

      !----- Leaf physiology: model selectors + shared biochemistry (non-PFT). ------------!
      call req_stomatal_model(tm, 'leaf_physiology.stomatal_model',     cfg%stomatal_model,     miss)
      call req_temp_response (tm, 'leaf_physiology.temp_response_form', cfg%temp_response_form, miss)
      call req_colimitation  (tm, 'leaf_physiology.colimitation',       cfg%colimitation,       miss)
      call req_l(tm, 'leaf_physiology.use_boundary_layer', cfg%leaf_use_boundary_layer, miss)
      call req_r(tm, 'leaf_physiology.kc25',     cfg%kc25,     miss)
      call req_r(tm, 'leaf_physiology.ko25',     cfg%ko25,     miss)
      call req_r(tm, 'leaf_physiology.gstar25',  cfg%gstar25,  miss)
      call req_r(tm, 'leaf_physiology.ea_kc',    cfg%ea_kc,    miss)
      call req_r(tm, 'leaf_physiology.ea_ko',    cfg%ea_ko,    miss)
      call req_r(tm, 'leaf_physiology.ea_gstar', cfg%ea_gstar, miss)
      call req_r(tm, 'leaf_physiology.ea_vcmax', cfg%ea_vcmax, miss)
      call req_r(tm, 'leaf_physiology.ea_jmax',  cfg%ea_jmax,  miss)
      call req_r(tm, 'leaf_physiology.ea_rd',    cfg%ea_rd,    miss)
      call req_r(tm, 'leaf_physiology.hd_vcmax', cfg%hd_vcmax, miss)
      call req_r(tm, 'leaf_physiology.hd_jmax',  cfg%hd_jmax,  miss)
      call req_r(tm, 'leaf_physiology.hd_rd',    cfg%hd_rd,    miss)
      call req_r(tm, 'leaf_physiology.ds_vcmax', cfg%ds_vcmax, miss)
      call req_r(tm, 'leaf_physiology.ds_jmax',  cfg%ds_jmax,  miss)
      call req_r(tm, 'leaf_physiology.ds_rd',    cfg%ds_rd,    miss)
      call req_r(tm, 'leaf_physiology.o2_mol_frac',      cfg%o2_mol_frac,      miss)
      call req_r(tm, 'leaf_physiology.leaf_absorptance', cfg%leaf_absorptance, miss)
      call req_r(tm, 'leaf_physiology.phi_psii',         cfg%phi_psii,         miss)

      !----- PFT file (named in the main file). -------------------------------------------!
      call toml_parse_file(trim(cfg%pft_config), tp, found)
      if (.not. found) error stop 'meds_config: PFT config file not found: '//trim(cfg%pft_config)

      !----- PFT count comes from the length of the wood_density array. -------------------!
      call toml_real_array(tp, 'pft.wood_density', buf, nout)
      if (nout < 1_ik) error stop 'meds_config: pft.wood_density missing/empty in '//trim(cfg%pft_config)
      if (nout > MAXPFT) error stop 'meds_config: pft.wood_density lists more than MAXPFT PFTs in '//trim(cfg%pft_config)
      npft = nout
      call alloc_pft_table(cfg%pft, npft)
      cfg%pft%wood_density = buf(1:npft)

      call req_pa(tp, 'pft.dbh_critical',                     cfg%pft%dbh_critical,                     npft, miss)
      call req_pa(tp, 'pft.hgt_max',                          cfg%pft%hgt_max,                          npft, miss)
      call req_pa(tp, 'pft.reproduction_investment_fraction', cfg%pft%reproduction_investment_fraction, npft, miss)
      call req_pa(tp, 'pft.repro_carbon_efficiency',          cfg%pft%repro_carbon_efficiency,          npft, miss)
      call req_pa(tp, 'pft.seed_rain_recruits',               cfg%pft%seed_rain_recruits,               npft, miss)
      call req_pa_int(tp, 'pft.include_pft',                  cfg%pft%include_pft,                      npft, miss)
      call req_r(tp, 'pft.min_cohort_height',       cfg%pft%min_cohort_height,       miss)
      call req_r(tp, 'pft.min_reproduction_height', cfg%pft%min_reproduction_height, miss)

      !----- Leaf-photosynthesis per-PFT traits. ------------------------------------------!
      call req_pa_int(tp, 'pft.photosynthetic_pathway', cfg%pft%photosynthetic_pathway, npft, miss)
      call req_pa(tp, 'pft.vcmax25',          cfg%pft%vcmax25,          npft, miss)
      call req_pa(tp, 'pft.jmax_vcmax_ratio', cfg%pft%jmax_vcmax_ratio, npft, miss)
      call req_pa(tp, 'pft.tpu_vcmax_ratio',  cfg%pft%tpu_vcmax_ratio,  npft, miss)
      call req_pa(tp, 'pft.rd_vcmax_ratio',   cfg%pft%rd_vcmax_ratio,   npft, miss)
      call req_pa(tp, 'pft.kp25',             cfg%pft%kp25,             npft, miss)
      call req_pa(tp, 'pft.stomatal_g0',      cfg%pft%stomatal_g0,      npft, miss)
      call req_pa(tp, 'pft.stomatal_g1',      cfg%pft%stomatal_g1,      npft, miss)
      call req_pa(tp, 'pft.stomatal_d0',      cfg%pft%stomatal_d0,      npft, miss)
      call req_pa(tp, 'pft.quantum_yield_c4', cfg%pft%quantum_yield_c4, npft, miss)
      call req_pa(tp, 'pft.theta_j',          cfg%pft%theta_j,          npft, miss)
      call req_pa(tp, 'pft.theta_cj_c4',      cfg%pft%theta_cj_c4,      npft, miss)
      call req_pa(tp, 'pft.theta_ic_c4',      cfg%pft%theta_ic_c4,      npft, miss)
      call req_pa(tp, 'pft.katul_lambda25',   cfg%pft%katul_lambda25,   npft, miss)
      call req_pa(tp, 'pft.wstress_psi_open', cfg%pft%wstress_psi_open, npft, miss)
      call req_pa(tp, 'pft.wstress_psi_close',cfg%pft%wstress_psi_close,npft, miss)
      call req_pa(tp, 'pft.wstress_lambda_exp',cfg%pft%wstress_lambda_exp,npft, miss)
      call req_pa(tp, 'pft.wstress_sref_stomata',cfg%pft%wstress_sref_stomata,npft, miss)

      !----- Carbon-dynamics per-PFT traits (meds_plant_carbon_dynamics). ------------------!
      call req_pa(tp, 'pft.sla',                    cfg%pft%sla,                    npft, miss)
      call req_pa(tp, 'pft.root_to_leaf_ratio',     cfg%pft%root_to_leaf_ratio,     npft, miss)
      call req_pa(tp, 'pft.huber_value',            cfg%pft%huber_value,            npft, miss)
      call req_pa(tp, 'pft.aboveground_frac',       cfg%pft%aboveground_frac,       npft, miss)
      call req_pa(tp, 'pft.storage_cushion',        cfg%pft%storage_cushion,        npft, miss)
      call req_pa(tp, 'pft.growth_resp_factor',     cfg%pft%growth_resp_factor,     npft, miss)
      call req_pa(tp, 'pft.leaf_turnover_rate',     cfg%pft%leaf_turnover_rate,     npft, miss)
      call req_pa(tp, 'pft.fineroot_turnover_rate', cfg%pft%fineroot_turnover_rate, npft, miss)
      call req_pa(tp, 'pft.wood_carbon_density',    cfg%pft%wood_carbon_density,    npft, miss)
      call req_pa_int(tp, 'pft.evergreen',          cfg%pft%evergreen,              npft, miss)

      !----- Leaf-phenology per-PFT cue params (opt-in; only required when phenology_on). ---!
      call load_phenology_pft(tp, cfg, npft, miss)

      call req_r(tp, 'camac.mort_rho_ref',   cfg%pft%mort_rho_ref,   miss)
      call req_r(tp, 'camac.mort_gamma_0',   cfg%pft%mort_gamma_0,   miss)
      call req_r(tp, 'camac.mort_gamma_exp', cfg%pft%mort_gamma_exp, miss)
      call req_r(tp, 'camac.mort_alpha_0',   cfg%pft%mort_alpha_0,   miss)
      call req_r(tp, 'camac.mort_alpha_exp', cfg%pft%mort_alpha_exp, miss)
      call req_r(tp, 'camac.mort_beta_0',    cfg%pft%mort_beta_0,    miss)
      call req_r(tp, 'camac.mort_beta_exp',  cfg%pft%mort_beta_exp,  miss)

      call req_r(tp, 'allometry.b1Ht',       cfg%allom%b1Ht,      miss)
      call req_r(tp, 'allometry.b2Ht',       cfg%allom%b2Ht,      miss)
      call req_r(tp, 'allometry.agb_c1',     cfg%allom%agb_c1,    miss)
      call req_r(tp, 'allometry.agb_c2',     cfg%allom%agb_c2,    miss)
      call req_r(tp, 'allometry.ca_b1',      cfg%allom%ca_b1,     miss)
      call req_r(tp, 'allometry.ca_b2',      cfg%allom%ca_b2,     miss)
      call req_r(tp, 'allometry.lai_b1',     cfg%allom%lai_b1,    miss)
      call req_r(tp, 'allometry.lai_b2',     cfg%allom%lai_b2,    miss)
      call req_r(tp, 'allometry.light_ext',  cfg%allom%light_ext, miss)

      !----- Presence check: abort with the full list of any missing required keys. -------!
      if (miss%n > 0_ik) then
         write(*,'(a,i0,a)') ' meds_config: ', miss%n, ' required parameter(s) missing:'
         do i = 1_ik, miss%n
            write(*,'(3a)') '   - ', trim(miss%key(i)), ''
         end do
         error stop 'meds_config: incomplete configuration (see missing keys above)'
      end if

      !----- Install allometry + compute every derived quantity (one consolidation point). -----!
      call derive_parameters(cfg)

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
      write(u,'(a)') 'pft,wood_density,dbh_critical,hgt_max,'                                       &
           //'reproduction_investment_fraction,repro_carbon_efficiency,'                            &
           //'mort_gamma,mort_alpha,mort_beta,seed_rain_recruits,include_pft,'                         &
           //'min_cohort_height,min_reproduction_height,'                                              &
           //'photosynthetic_pathway,vcmax25,jmax25,tpu25,rd25,kp25,'                                  &
           //'stomatal_g0,stomatal_g1,stomatal_d0,quantum_yield_c4,theta_j,theta_cj_c4,theta_ic_c4,'   &
           //'katul_lambda25,wstress_psi_open,wstress_psi_close,wstress_lambda_exp,wstress_sref_stomata,' &
           //'sla,root_to_leaf_ratio,huber_value,aboveground_frac,storage_cushion,growth_resp_factor,' &
           //'leaf_turnover_rate,fineroot_turnover_rate,wood_carbon_density,evergreen'
      associate (p => cfg%pft)
         do pf = 1_ik, p%n
            write(u,'(i0,9(",",es15.8),",",i0,2(",",es15.8),",",i0,17(",",es15.8),9(",",es15.8),",",i0)') &
                 pf, p%wood_density(pf), p%dbh_critical(pf), p%hgt_max(pf),                             &
                 p%reproduction_investment_fraction(pf), p%repro_carbon_efficiency(pf),                &
                 p%mort_gamma(pf), p%mort_alpha(pf), p%mort_beta(pf), p%seed_rain_recruits(pf),         &
                 p%include_pft(pf), p%min_cohort_height, p%min_reproduction_height,                     &
                 p%photosynthetic_pathway(pf), p%vcmax25(pf), p%jmax25(pf), p%tpu25(pf),                &
                 p%rd25(pf), p%kp25(pf), p%stomatal_g0(pf), p%stomatal_g1(pf), p%stomatal_d0(pf),       &
                 p%quantum_yield_c4(pf), p%theta_j(pf), p%theta_cj_c4(pf), p%theta_ic_c4(pf),           &
                 p%katul_lambda25(pf), p%wstress_psi_open(pf), p%wstress_psi_close(pf),                 &
                 p%wstress_lambda_exp(pf), p%wstress_sref_stomata(pf),                                  &
                 p%sla(pf), p%root_to_leaf_ratio(pf), p%huber_value(pf), p%aboveground_frac(pf),        &
                 p%storage_cushion(pf), p%growth_resp_factor(pf), p%leaf_turnover_rate(pf),             &
                 p%fineroot_turnover_rate(pf), p%wood_carbon_density(pf), p%evergreen(pf)
         end do
      end associate
      close(u)
      write(*,'(2a)') ' pft   : ', trim(path)
   end subroutine write_pft_params_csv

end module meds_config_io
