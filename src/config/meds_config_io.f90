! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config_io -- load a run configuration from TWO mandatory TOML files and a presence- !
! map check (docs/configuration.md). In src/config (config input) with NO external         !
! dependency, so it is always compiled.                                                    !
!                                                                                          !
! Philosophy: NO hard-coded model parameters. Every parameter is supplied by the user via  !
! TOML (or DERIVED from supplied ones); a missing parameter or a missing file is a HARD    !
! ERROR. The exceptions are the optional keys -- solver settings, process switches and a   !
! few blocks ([soil], [snow], [aerodynamics], [soil_carbon]) -- which fall back to their   !
! defaults. The MAIN file (passed on the CLI) holds all non-PFT parameters and names the   !
! PFT file; the PFT file holds every PFT trait, its allometry and its hydraulics among     !
! them. Presence is tracked while reading (a "missing key" list); after both files are     !
! read, any missing required key aborts with the full list. Derived quantities are then    !
! computed (derive_parameters).                                                            !
!==========================================================================================!
module meds_config_io
   use meds_kinds,      only : wp, ik
   use meds_allometry_lib, only : HEIGHT_POWER, HEIGHT_CHAVE
   use meds_config,     only : meds_config_t, derive_parameters, INIT_CENSUS, INIT_RESTART
   use meds_fast_opts,  only : INTEG_ARK, INTEG_RK45, CTRL_L0_FIXED, CTRL_L1_ADAPTIVE,                &
                               CTRL_L2_STRICT, CTRL_I, CTRL_PI
   use meds_config,     only : pft_stomata_psi_onset, pft_wstress_psi_open, pft_wstress_psi_close
   use meds_soil_column_config, only : soil_column_config_t
   use meds_region_opts, only : RUN_MODE_SITE, RUN_MODE_REGION
   use meds_maps_opts, only : soil_maps_on
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   use meds_leaf_opts,     only : SM_LEUNING, SM_MEDLYN, SM_KATUL, COLIM_MIN, COLIM_QUADRATIC
   use meds_therm_lib,     only : TEMP_RESPONSE_ARRHENIUS, TEMP_RESPONSE_PEAKED
   use meds_forcing_config, only : LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO, CO2_SOURCE_CONST, CO2_SOURCE_FILE
   use meds_forcing_config, only : HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND,                     &
                                   WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL
   use meds_forcing_config, only : forcing_config_t,                                            &
                                   MET_BACKEND_CONST, MET_BACKEND_ED_DEFAULT, MET_BACKEND_ED_ERA5LAND, &
                                   METAVG_END, METAVG_BEGIN,                                     &
                                   SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN,      &
                                   LW_FILE, LW_SYNTHESIZE, CLAMP_ERROR, CLAMP_HOLD,              &
                                   GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   use meds_output_config, only : output_config_t, GRP_STRUCTURE, GRP_CARBON, GRP_WATER,          &
                                 GRP_RADIATION, GRP_ECOPHYS, GRP_BIOGEOCHEM, GRP_FORCING,        &
                                 AXIS_COHORT, AXIS_PATCH, AXIS_PFT, AXIS_SIZE, AXIS_SOIL_PATCH,  &
                                 MAX_DBH_EDGE,                                                   &
                                   GRP_ENERGY, GRP_NUMERICS, SYNC_FLUSH, SYNC_NEVER, FC_DAY, FC_MONTH, FC_YEAR, &
                                   FC_RUN
   use meds_time,       only : meds_time_t, string_to_time
   use meds_pft_params, only : alloc_pft_table
   use meds_biophysics_opts, only : ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t,        &
                                    SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use meds_biogeochem_opts, only : decomp_opts_t, DECOMP_STEP_EULER, DECOMP_STEP_EXPM,            &
                                    DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   use meds_toml,       only : toml_table_t, toml_parse_file, toml_has, toml_has_section,    &
                              toml_int, toml_real, toml_table_count,                        &
                               toml_logical, toml_string, toml_real_array, toml_record_reset,  &
                               toml_is_array, toml_string_array
   use meds_config_check, only : key_report_t, check_config_keys, check_key_use, validate_config,     &
                                 MAIN_KEYS, PFT_KEYS, SITE_GROUND_KEYS
   implicit none
   private

   public :: load_meds_config, write_pft_params_csv
   public :: load_soil_opts, load_energy_opts, load_snow_params, load_aero_cfg   !< exposed for unit tests

   integer(ik), parameter :: MAXPFT  = 64_ik    !< largest PFT count the array buffers handle
   integer,     parameter :: MAXMISS = 256       !< largest missing-key list
   integer,     parameter :: NAME_LEN = 16       !< longest selector name

   !----- Accumulator of required keys that were absent or unreadable (the presence map). -----!
   type :: keymiss_t
      integer(ik)        :: n = 0_ik
      character(len=200) :: key(MAXMISS) = ''
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

   !----- A selector: a string key naming one of `names`, read into the matching code of `codes`.  !
   !      Absent, a required one is recorded missing and an optional one keeps `out`. A name not in   !
   !      the list is recorded with the names it may take, or, with no presence map to record it in   !
   !      (the stand-alone block readers), stops the run. ------------------------------------------!
   subroutine req_choice(t, key, names, codes, out, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key, names(:)
      integer(ik),        intent(in)    :: codes(:)
      integer(ik),        intent(inout) :: out
      type(keymiss_t),    intent(inout) :: m
      if (.not. toml_has(t, key)) then
         call note_missing(m, key)
         return
      end if
      call read_choice(t, key, names, codes, out, m)
   end subroutine req_choice

   subroutine opt_choice(t, key, names, codes, out, m)
      type(toml_table_t),        intent(in)    :: t
      character(len=*),          intent(in)    :: key, names(:)
      integer(ik),               intent(in)    :: codes(:)
      integer(ik),               intent(inout) :: out
      type(keymiss_t), optional, intent(inout) :: m
      character(len=:), allocatable :: kept
      integer :: i
      if (.not. toml_has(t, key)) then
         !----- Absent: `out` stands, and the parameter record lists its name as the default. -----!
         kept = ''
         do i = 1, size(codes)
            if (codes(i) == out) kept = trim(names(i))
         end do
         kept = toml_string(t, key, kept)
         return
      end if
      call read_choice(t, key, names, codes, out, m)
   end subroutine opt_choice

   subroutine read_choice(t, key, names, codes, out, m)
      type(toml_table_t),        intent(in)    :: t
      character(len=*),          intent(in)    :: key, names(:)
      integer(ik),               intent(in)    :: codes(:)
      integer(ik),               intent(inout) :: out
      type(keymiss_t), optional, intent(inout) :: m
      character(len=:), allocatable :: s, allowed
      integer :: i
      s = trim(toml_string(t, key, ''))
      do i = 1, size(names)
         if (s == trim(names(i))) then
            out = codes(i)
            return
         end if
      end do
      allowed = trim(names(1))
      do i = 2, size(names)
         allowed = allowed//' | '//trim(names(i))
      end do
      if (.not. present(m)) error stop 'load_meds_config: '//key//' = "'//s//'" is not one of '//allowed
      call note_missing(m, key//' = "'//s//'" is not one of '//allowed)
   end subroutine read_choice


   !=======================================================================================!
   !  Fast-loop biophysics run-config loaders ([soil]/[snow]/[aerodynamics], and the solver !
   !  settings in [fast]). All OPTIONAL: every key is read with its current (default)       !
   !  value as the getter fallback, so an absent key or block keeps the type's default.     !
   !  String selectors keep the default when absent and HARD-ERROR when present but         !
   !  unrecognized (like the [forcing] enum mappers, but optional -- no note_missing).      !
   !  Public for unit testing.                                                              !
   !                                                                                       !
   !  [soil] is the ground, the one block that describes it: the column's layers and        !
   !  texture (soil_column_config_t), its surface (soil_opts_t: pond, dry surface layer,    !
   !  optics) and both bottom boundaries (soil_opts_t for the water, energy_opts_t for the  !
   !  heat). The Fortran types are grouped by the solver that reads them, so one block      !
   !  fills three types.                                                                   !
   !=======================================================================================!
   !---------------------------------------------------------------------------------------!
   ! [soil] -> the soil column: its layers, hydraulic texture and thermal texture. Every    !
   ! key is optional and defaults to its soil_column_config_t value.                       !
   !---------------------------------------------------------------------------------------!
   subroutine load_soil_column(tm, c)
      type(toml_table_t),         intent(in)    :: tm
      type(soil_column_config_t), intent(inout) :: c
      !----- A layer_bottoms list without a depth sets the column to end at its last bottom. ---!
      call toml_real_array(tm, 'soil.layer_bottoms', c%layer_bottoms, c%n_layer_bottoms)
      if (c%n_layer_bottoms > 0_ik) c%depth = c%layer_bottoms(c%n_layer_bottoms)
      c%depth       = toml_real(tm, 'soil.depth',       c%depth)
      call opt_choice(tm, 'soil.retention', [character(len=NAME_LEN) :: 'van_genuchten', 'campbell'],   &
                      [SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL], c%retention)
      c%theta_sat   = toml_real(tm, 'soil.theta_sat',   c%theta_sat)
      c%theta_res   = toml_real(tm, 'soil.theta_res',   c%theta_res)
      c%ksat        = toml_real(tm, 'soil.ksat',        c%ksat)
      c%curve_par_a = toml_real(tm, 'soil.curve_par_a', c%curve_par_a)
      c%curve_par_n = toml_real(tm, 'soil.curve_par_n', c%curve_par_n)
      c%psi_fc      = toml_real(tm, 'soil.psi_fc',      c%psi_fc)
      c%solid_conductivity = toml_real(tm, 'soil.solid_conductivity', c%solid_conductivity)
      c%dry_conductivity   = toml_real(tm, 'soil.dry_conductivity',   c%dry_conductivity)
      c%dry_heat_capacity  = toml_real(tm, 'soil.dry_heat_capacity',  c%dry_heat_capacity)
   end subroutine load_soil_column

   !----- The soil water's options: the ground's surface and water boundary from [soil], the   !
   !      water network's tolerances and the canopy interception from [fast]. ------------------!
   subroutine load_soil_opts(tm, s)
      type(toml_table_t), intent(in)    :: tm
      type(soil_opts_t),  intent(inout) :: s
      call load_water_bottom_bc(tm, s%water_bottom_bc, s%drain_factor)
      s%water_atol           = toml_real   (tm, 'fast.water_atol',      s%water_atol)
      s%w_pond_max     = toml_real   (tm, 'soil.w_pond_max',      s%w_pond_max)
      s%intercept_alpha= toml_real   (tm, 'fast.intercept_alpha', s%intercept_alpha)
      s%intercept_k    = toml_real   (tm, 'fast.intercept_k',     s%intercept_k)
      s%dsl_dmax       = toml_real   (tm, 'soil.dsl_dmax',        s%dsl_dmax)
      s%dsl_theta_init = toml_real   (tm, 'soil.dsl_theta_init',  s%dsl_theta_init)
      s%ground_albedo_vis = toml_real(tm, 'soil.ground_albedo_vis', s%ground_albedo_vis)
      s%ground_albedo_nir = toml_real(tm, 'soil.ground_albedo_nir', s%ground_albedo_nir)
      s%ground_emissivity = toml_real(tm, 'soil.ground_emissivity', s%ground_emissivity)
   end subroutine load_soil_opts

   !----- [soil].water_bottom_bc: the column's water boundary and its drainage factor, kept when   !
   !      the key is absent. ----------------------------------------------------------------------!
   subroutine load_water_bottom_bc(tm, bc, drain_factor)
      type(toml_table_t), intent(in)    :: tm
      integer(ik),        intent(inout) :: bc
      real(wp),           intent(inout) :: drain_factor
      integer(ik), parameter :: FREE_DRAIN = 1_ik, BEDROCK = 2_ik, AQUIFER = 3_ik
      integer(ik) :: choice
      if (.not. toml_has(tm, 'soil.water_bottom_bc')) return
      choice = FREE_DRAIN
      call opt_choice(tm, 'soil.water_bottom_bc', [character(len=NAME_LEN) :: 'free_drain', 'bedrock', 'aquifer'], &
                      [FREE_DRAIN, BEDROCK, AQUIFER], choice)
      !----- free_drain and bedrock are the drainage boundary's two ends, f = 1 and f = 0 !
      !      (docs/science/soil_biophysics.md §Boundary conditions). ---------------------!
      bc = merge(SOIL_BC_AQUIFER, SOIL_BC_DRAIN, choice == AQUIFER)
      if (choice == FREE_DRAIN) drain_factor = 1.0_wp
      if (choice == BEDROCK)    drain_factor = 0.0_wp
   end subroutine load_water_bottom_bc

   !----- The soil heat's options: the ground's thermal bottom boundary from [soil], the budget  !
   !      check from [fast]. --------------------------------------------------------------------!
   subroutine load_energy_opts(tm, e)
      type(toml_table_t),  intent(in)    :: tm
      type(energy_opts_t), intent(inout) :: e
      call opt_choice(tm, 'soil.thermal_bottom_bc', [character(len=NAME_LEN) :: 'geothermal', 'dirichlet'], &
                      [ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET], e%thermal_bottom_bc)
      e%deep_temp   = toml_real(tm, 'soil.deep_temp',   e%deep_temp)
      e%deep_depth  = toml_real(tm, 'soil.deep_depth',  e%deep_depth)
      e%soil_energy_budget_atol        = toml_real(tm, 'fast.soil_energy_budget_atol', e%soil_energy_budget_atol)
      !----- HARD-STOP on a non-closing budget. Gates the stop on the fast loop's budget checks      !
      !      (meds_fast_ark, meds_fast_rk45) and the kernels' own closure checks. Default            !
      !      .false., so a production run reports a breach rather than halting on it; set it to      !
      !      halt when you are bisecting one.                                                        !
      !                                                                                              !
      !      This key MUST keep its reader. Without one the flag can be set only from Fortran,       !
      !      every configured run sits at .false., and the protocol "run with fast.debug_error       !
      !      = .true. so a breach HALTS" is silently vacuous: a run would complete because           !
      !      nothing was checking, not because nothing was wrong. -----------------------------------!
      e%debug_error = toml_logical(tm, 'fast.debug_error', e%debug_error)
   end subroutine load_energy_opts

   subroutine load_snow_params(tm, s)                   ! [snow] -> snow physical params
      type(toml_table_t), intent(in)    :: tm
      type(snow_params_t), intent(inout) :: s
      s%rho_snow            = toml_real(tm, 'snow.rho_snow',            s%rho_snow)
      s%k_snow              = toml_real(tm, 'snow.k_snow',              s%k_snow)
      s%albedo_vis_fresh    = toml_real(tm, 'snow.albedo_vis_fresh',    s%albedo_vis_fresh)
      s%albedo_nir_fresh    = toml_real(tm, 'snow.albedo_nir_fresh',    s%albedo_nir_fresh)
      s%albedo_vis_aged     = toml_real(tm, 'snow.albedo_vis_aged',     s%albedo_vis_aged)
      s%albedo_nir_aged     = toml_real(tm, 'snow.albedo_nir_aged',     s%albedo_nir_aged)
      s%snow_emiss          = toml_real(tm, 'snow.snow_emiss',          s%snow_emiss)
      s%liquid_holding_frac = toml_real(tm, 'snow.liquid_holding_frac', s%liquid_holding_frac)
      s%min_new_snow_mass   = toml_real(tm, 'snow.min_new_snow_mass',   s%min_new_snow_mass)
      s%tiny_snow_mass      = toml_real(tm, 'snow.tiny_snow_mass',      s%tiny_snow_mass)
      s%ny07_a              = toml_real(tm, 'snow.ny07_a',              s%ny07_a)
      s%ny07_m              = toml_real(tm, 'snow.ny07_m',              s%ny07_m)
   end subroutine load_snow_params

   subroutine load_aero_cfg(tm, a)                      ! [aerodynamics] -> canopy-aero constants
      type(toml_table_t), intent(in)    :: tm            !  (the Nusselt-coefficient tables keep
      type(aero_cfg_t),   intent(inout) :: a             !   their ED2/CLM defaults; not exposed)
      a%z0m_ratio        = toml_real(tm, 'aerodynamics.z0m_ratio',        a%z0m_ratio)
      a%d_ratio          = toml_real(tm, 'aerodynamics.d_ratio',          a%d_ratio)
      a%snow_rough       = toml_real(tm, 'aerodynamics.snow_rough',       a%snow_rough)
      a%zeta_m           = toml_real(tm, 'aerodynamics.zeta_m',           a%zeta_m)
      a%zeta_t           = toml_real(tm, 'aerodynamics.zeta_t',           a%zeta_t)
      a%zeta_max_stable  = toml_real(tm, 'aerodynamics.zeta_max_stable',  a%zeta_max_stable)
      a%n_iter_mo        = toml_int (tm, 'aerodynamics.n_iter_mo',        a%n_iter_mo)
      a%wc               = toml_real(tm, 'aerodynamics.wc',               a%wc)
      a%ustmin           = toml_real(tm, 'aerodynamics.ustmin',           a%ustmin)
      a%ubmin            = toml_real(tm, 'aerodynamics.ubmin',            a%ubmin)
      a%ugbmin           = toml_real(tm, 'aerodynamics.ugbmin',           a%ugbmin)
      a%gbhmos_min       = toml_real(tm, 'aerodynamics.gbhmos_min',       a%gbhmos_min)
      a%min_canopy_depth = toml_real(tm, 'aerodynamics.min_canopy_depth', a%min_canopy_depth)
      a%canopy_freeboard = toml_real(tm, 'aerodynamics.canopy_freeboard', a%canopy_freeboard)
      a%cs_dense         = toml_real(tm, 'aerodynamics.cs_dense',         a%cs_dense)
      a%gamma_g          = toml_real(tm, 'aerodynamics.gamma_g',          a%gamma_g)
      a%gbh_2_gbw        = toml_real(tm, 'aerodynamics.gbh_2_gbw',        a%gbh_2_gbw)
      a%kin_visc0        = toml_real(tm, 'aerodynamics.kin_visc0',        a%kin_visc0)
      a%dkin_visc        = toml_real(tm, 'aerodynamics.dkin_visc',        a%dkin_visc)
      a%th_diff0         = toml_real(tm, 'aerodynamics.th_diff0',         a%th_diff0)
      a%dth_diff         = toml_real(tm, 'aerodynamics.dth_diff',         a%dth_diff)
      a%t_ref_air        = toml_real(tm, 'aerodynamics.t_ref_air',        a%t_ref_air)
   end subroutine load_aero_cfg

   !----- [soil_carbon]: the switch, the decomposition selectors and rates, and the cold start. Each   !
   !      absent key keeps the record's default, so the record, not the reader, decides it. ----------!
   subroutine load_soil_carbon_opts(tm, s, m)
      type(toml_table_t),  intent(in)    :: tm
      type(decomp_opts_t), intent(inout) :: s
      type(keymiss_t),     intent(inout) :: m
      s%soil_carbon_on     = toml_logical(tm, 'soil_carbon.soil_carbon_on',     s%soil_carbon_on)
      s%spinup_steady      = toml_logical(tm, 'soil_carbon.spinup_steady',      s%spinup_steady)
      s%spinup_xi          = toml_real   (tm, 'soil_carbon.spinup_xi',          s%spinup_xi)
      s%spinup_labile_grnd = toml_real   (tm, 'soil_carbon.spinup_labile_grnd', s%spinup_labile_grnd)
      s%spinup_labile_soil = toml_real   (tm, 'soil_carbon.spinup_labile_soil', s%spinup_labile_soil)
      s%spinup_struct_grnd = toml_real   (tm, 'soil_carbon.spinup_struct_grnd', s%spinup_struct_grnd)
      s%spinup_struct_soil = toml_real   (tm, 'soil_carbon.spinup_struct_soil', s%spinup_struct_soil)
      call opt_choice(tm, 'soil_carbon.decomp_scheme', [character(len=NAME_LEN) :: 'ed2', 'century5'],     &
                      [DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5], s%decomp_scheme, m)
      call opt_choice(tm, 'soil_carbon.step_solver', [character(len=NAME_LEN) :: 'euler', 'expm'],          &
                      [DECOMP_STEP_EULER, DECOMP_STEP_EXPM], s%step_solver, m)
      s%k_fast   = toml_real(tm, 'soil_carbon.k_fast',   s%k_fast)
      s%k_struct = toml_real(tm, 'soil_carbon.k_struct', s%k_struct)
      s%k_micr   = toml_real(tm, 'soil_carbon.k_micr',   s%k_micr)
      s%k_slow   = toml_real(tm, 'soil_carbon.k_slow',   s%k_slow)
      s%k_pass   = toml_real(tm, 'soil_carbon.k_pass',   s%k_pass)
      s%er_fast          = toml_real(tm, 'soil_carbon.er_fast',          s%er_fast)
      s%er_struct_nonlig = toml_real(tm, 'soil_carbon.er_struct_nonlig', s%er_struct_nonlig)
      s%er_struct_lig    = toml_real(tm, 'soil_carbon.er_struct_lig',    s%er_struct_lig)
      s%er_slow          = toml_real(tm, 'soil_carbon.er_slow',          s%er_slow)
      s%er_pass          = toml_real(tm, 'soil_carbon.er_pass',          s%er_pass)
      s%er_micr_int      = toml_real(tm, 'soil_carbon.er_micr_int',      s%er_micr_int)
      s%er_micr_slp      = toml_real(tm, 'soil_carbon.er_micr_slp',      s%er_micr_slp)
      s%e_lignin         = toml_real(tm, 'soil_carbon.e_lignin',         s%e_lignin)
      s%xsand = toml_real(tm, 'soil_carbon.xsand', s%xsand)
      s%xclay = toml_real(tm, 'soil_carbon.xclay', s%xclay)
      s%fx_micr_pass_int = toml_real(tm, 'soil_carbon.fx_micr_pass_int', s%fx_micr_pass_int)
      s%fx_micr_pass_slp = toml_real(tm, 'soil_carbon.fx_micr_pass_slp', s%fx_micr_pass_slp)
      s%fx_slow_pass_int = toml_real(tm, 'soil_carbon.fx_slow_pass_int', s%fx_slow_pass_int)
      s%fx_slow_pass_slp = toml_real(tm, 'soil_carbon.fx_slow_pass_slp', s%fx_slow_pass_slp)
      s%rh_q10                = toml_real(tm, 'soil_carbon.rh_q10',                s%rh_q10)
      s%rh_t_ref              = toml_real(tm, 'soil_carbon.rh_t_ref',              s%rh_t_ref)
      s%resp_temp_increase    = toml_real(tm, 'soil_carbon.resp_temp_increase',    s%resp_temp_increase)
      s%resp_temp_ref         = toml_real(tm, 'soil_carbon.resp_temp_ref',         s%resp_temp_ref)
      s%resp_opt_water        = toml_real(tm, 'soil_carbon.resp_opt_water',        s%resp_opt_water)
      s%resp_water_below_opt  = toml_real(tm, 'soil_carbon.resp_water_below_opt',  s%resp_water_below_opt)
      s%resp_water_above_opt  = toml_real(tm, 'soil_carbon.resp_water_above_opt',  s%resp_water_above_opt)
   end subroutine load_soil_carbon_opts


   !----- [allometry].height_allometry, one for every PFT or one per PFT: "power" (the        !
   !      default) or "chave2014", the log-quadratic, whose PFTs take Chave's coefficients        !
   !      unless given their own fit, and the site's stress index chave_e. ------------------------!
   subroutine opt_height_allometry(t, cfg, npft, m)
      type(toml_table_t), intent(in) :: t ; type(meds_config_t), intent(inout) :: cfg
      integer(ik), intent(in) :: npft ; type(keymiss_t), intent(inout) :: m
      character(len=16) :: forms(MAXPFT)
      integer(ik) :: nform, pf, k
      real(wp)    :: chave_e
      forms = 'power'
      if (toml_has(t, 'allometry.height_allometry')) then
         call toml_string_array(t, 'allometry.height_allometry', forms, nform)
         if (nform == 1_ik) forms(1:npft) = forms(1)
         if (nform /= 1_ik .and. nform /= npft) then
            call note_missing(m, 'allometry.height_allometry')
            return
         end if
      end if
      associate (allom => cfg%pft%allom)
         do pf = 1_ik, npft
            k = int(findloc([character(len=16) :: 'power', 'chave2014'], forms(pf), dim=1), ik)
            if (k == 0_ik) then
               call note_missing(m, 'allometry.height_allometry = "'//trim(forms(pf))//                 &
                                    '" is not one of power | chave2014')
               return
            end if
            allom(pf)%height_allometry = merge(HEIGHT_POWER, HEIGHT_CHAVE, k == 1_ik)
         end do
         if (any(allom(1:npft)%height_allometry == HEIGHT_CHAVE)) then
            call opt_pa_bcast(t, 'allometry.chave_a',   allom%chave_a,   npft, m)
            call opt_pa_bcast(t, 'allometry.chave_b',   allom%chave_b,   npft, m)
            call opt_pa_bcast(t, 'allometry.chave_c',   allom%chave_c,   npft, m)
            call opt_pa_bcast(t, 'allometry.chave_rse', allom%chave_rse, npft, m)
            !----- E is the site's, the same for every PFT: from [maps].chave_e_map at each       !
            !      polygon's cell, or set here. ------------------------------------------------!
            if (len_trim(cfg%maps%chave_e_map) == 0) then
               call req_r(t, 'allometry.chave_e', chave_e, m)
               allom(1:npft)%chave_e = chave_e
            end if
         end if
      end associate
   end subroutine opt_height_allometry


   !----- Load the [forcing]/[site] block (docs/science/forcing.md). GATED on                     !
   !      forcing.forcing_on (a DEFAULTED read, so a config with no [forcing] block runs on       !
   !      constant forcing). When ON, the rest are REQUIRED (the no-silent-defaults rule). The     !
   !      source, the location and the switches pick which keys apply; check_key_use refuses the   !
   !      rest. ------------------------------------------------------------------------------------!
   subroutine load_forcing_config(t, cfg, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      type(keymiss_t),     intent(inout) :: m
      real(wp)    :: gamma_month(12)
      integer(ik) :: nout
      cfg%forcing = forcing_config_t()                                  ! Ithaca/ERA5-Land defaults
      cfg%forcing%forcing_on = toml_logical(t, 'forcing.forcing_on', .false.)   ! opt-in gate (defaulted)
      if (.not. cfg%forcing%forcing_on) return
      call req_choice(t, 'forcing.format', [character(len=NAME_LEN) :: 'ED_default', 'ED_ERA5land', 'const'],    &
                      [MET_BACKEND_ED_DEFAULT, MET_BACKEND_ED_ERA5LAND, MET_BACKEND_CONST], cfg%forcing%format, m)
      !----- The archive (format = "ED_ERA5land"; docs/science/forcing.md §1) finds the site's cell  !
      !      and its orography itself; a forcing file is named, with its grid slot. -----------------!
      if (cfg%forcing%format == MET_BACKEND_ED_ERA5LAND) then
         call req_s(t, 'forcing.data_path', cfg%forcing%data_path, m)
         !----- A region's polygons sit on their cells, so no site-to-cell distance exists. ------!
         if (cfg%mode /= RUN_MODE_REGION)                                                     &
            call req_r(t, 'forcing.max_distance_km', cfg%forcing%max_distance_km, m)
         cfg%forcing%file_template = toml_string(t, 'forcing.file_template', '')
         cfg%forcing%static_file   = toml_string(t, 'forcing.static_file',   '')
      else
         call req_s(t, 'forcing.path',       cfg%forcing%path,       m)
         call req_i(t, 'forcing.grid_index', cfg%forcing%grid_index, m)
         call req_choice(t, 'forcing.grid_match', [character(len=NAME_LEN) :: 'explicit', 'nearest'],           &
                         [GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST], cfg%forcing%grid_match, m)
      end if
      call req_dur(t, 'forcing.timestep', cfg%forcing%timestep, m)
      !----- Only the end- and begin-stamped conventions have a disaggregation (forcing.md §5). ----!
      call req_choice(t, 'forcing.avg_convention', [character(len=NAME_LEN) :: 'end', 'begin'],                &
                      [METAVG_END, METAVG_BEGIN], cfg%forcing%avg_convention, m)
      call req_choice(t, 'forcing.sw_partition',                                                     &
                      [character(len=NAME_LEN) :: 'passthrough', 'clearidx', 'weiss_norman'],         &
                      [SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN], cfg%forcing%sw_partition, m)
      call req_choice(t, 'forcing.lwdown_source', [character(len=NAME_LEN) :: 'file', 'synthesize'],           &
                      [LW_FILE, LW_SYNTHESIZE], cfg%forcing%lwdown_source, m)
      !----- Longwave-synthesis parameters (docs/science/forcing.md §11), optional. ---------------!
      call opt_choice(t, 'forcing.lw_clear_form', [character(len=NAME_LEN) :: 'brutsaert', 'idso'],            &
                      [LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO], cfg%forcing%lw_clear_form, m)
      cfg%forcing%lw_cloud_a = toml_real(t, 'forcing.lw_cloud_a', cfg%forcing%lw_cloud_a)
      call req_choice(t, 'forcing.start_clamp', [character(len=NAME_LEN) :: 'error', 'hold'],                  &
                      [CLAMP_ERROR, CLAMP_HOLD], cfg%forcing%start_clamp, m)
      call req_l(t, 'forcing.recycle', cfg%forcing%recycle, m)
      !----- The recycle WINDOW is declared, not inferred (docs/science/forcing.md §9). -----------!
      if (cfg%forcing%recycle) then
         call req_date(t, 'forcing.recycle_start', cfg%forcing%recycle_start, m)
         call req_date(t, 'forcing.recycle_end',   cfg%forcing%recycle_end,   m)
      end if
      !----- CO2 (docs/science/forcing.md §12): one source for every backend, never the met file. ---!
      call opt_choice(t, 'forcing.co2_source', [character(len=NAME_LEN) :: 'const', 'file'],                   &
                      [CO2_SOURCE_CONST, CO2_SOURCE_FILE], cfg%forcing%co2_source, m)
      if (cfg%forcing%co2_source == CO2_SOURCE_CONST) then
         call req_r(t, 'forcing.co2_const', cfg%forcing%co2_const, m)
      else
         call req_s(t, 'forcing.co2_file', cfg%forcing%co2_file, m)
      end if
      !----- The location: a region's polygons each take theirs from their cell. -----------------!
      if (cfg%mode /= RUN_MODE_REGION) then
         call req_r(t, 'site.latitude',  cfg%forcing%latitude,  m)
         call req_r(t, 'site.longitude', cfg%forcing%longitude, m)
         call req_r(t, 'site.elevation', cfg%forcing%elevation,   m)
      end if
      !----- The forcing's own vertical frame (docs/science/forcing.md §8). ---------------------!
      call req_r(t, 'forcing.tq_height',   cfg%forcing%tq_height,   m)
      call req_r(t, 'forcing.wind_height', cfg%forcing%wind_height, m)
      call req_choice(t, 'forcing.height_above', [character(len=NAME_LEN) :: 'zero_plane', 'ground'],          &
                      [HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND], cfg%forcing%height_above, m)
      call req_choice(t, 'forcing.wind_exposure', [character(len=NAME_LEN) :: 'open_terrain', 'local'],        &
                      [WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL], cfg%forcing%wind_exposure, m)
      if (cfg%forcing%wind_exposure == WIND_EXPOSURE_OPEN_TERRAIN) then
         call req_r(t, 'forcing.wind_exposure_z0',     cfg%forcing%wind_exposure_z0,     m)
         call req_r(t, 'forcing.wind_blending_height', cfg%forcing%wind_blending_height, m)
      end if
      !----- The terrain lapse (docs/science/forcing.md §8, eqs. (13)-(14)): one lapse rate for     !
      !      the year, or twelve (January .. December), read only when the lapse is on. -----------!
      call req_l(t, 'site.apply_elevation_lapse', cfg%forcing%apply_elevation_lapse, m)
      if (.not. cfg%forcing%apply_elevation_lapse) return
      if (toml_is_array(t, 'site.lapse_rate_tair')) then
         call toml_real_array(t, 'site.lapse_rate_tair', gamma_month, nout)
         if (nout == 12_ik) then
            cfg%forcing%lapse_rate_tair = gamma_month
         else
            call note_missing(m, 'site.lapse_rate_tair takes one value, or a list of twelve (January .. December)')
         end if
      else if (toml_has(t, 'site.lapse_rate_tair')) then
         cfg%forcing%lapse_rate_tair = toml_real(t, 'site.lapse_rate_tair', 0.0_wp)
      else
         call note_missing(m, 'site.lapse_rate_tair')
      end if
      if (cfg%forcing%format /= MET_BACKEND_ED_ERA5LAND)                                        &
         call req_r(t, 'site.grid_elevation', cfg%forcing%grid_elevation, m)
   end subroutine load_forcing_config

   !----- [maps] (docs/configuration.md, "Input maps"): the paths of the static maps. The three     !
   !      soil maps go together: the depth map gives the column, the texture and chemistry maps     !
   !      its layers. Each key the soil maps provide is refused when they are set (check_key_use),  !
   !      so a run never holds two values for one quantity. --------------------------------------!
   subroutine load_maps_config(t, cfg)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      integer(ik) :: n_set
      cfg%maps%soil_depth_map     = toml_string(t, 'maps.soil_depth_map',     '')
      cfg%maps%soil_texture_map   = toml_string(t, 'maps.soil_texture_map',   '')
      cfg%maps%soil_chemistry_map = toml_string(t, 'maps.soil_chemistry_map', '')
      cfg%maps%chave_e_map        = toml_string(t, 'maps.chave_e_map',        '')
      n_set = count([len_trim(cfg%maps%soil_depth_map) > 0, len_trim(cfg%maps%soil_texture_map) > 0, &
                     len_trim(cfg%maps%soil_chemistry_map) > 0])
      if (n_set /= 0_ik .and. n_set /= 3_ik)                                                         &
         error stop 'load_meds_config: [maps] needs soil_depth_map, soil_texture_map and '//       &
                    'soil_chemistry_map together'
   end subroutine load_maps_config

   !----- [[sites]] (docs/configuration.md, "Several sites"): the polygon's sites, one table each,   !
   !      after [soil] and [maps]. A table that sets any of [soil]'s column keys gives its site its  !
   !      own ground, [soil]'s values for each key it leaves out; with the soil maps every site      !
   !      takes the cell's soil (check_key_use). Each table is read through a view under            !
   !      [soil]'s names (site_view), so the [soil] readers serve it and the parameter record lists  !
   !      its values with the table as their source. ------------------------------------------------!
   subroutine load_sites(t, cfg, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      type(keymiss_t),     intent(inout) :: m
      type(toml_table_t) :: v
      integer(ik) :: ns, is
      ns = toml_table_count(t, 'sites')
      if (ns == 0_ik) return
      allocate(cfg%sites(ns))
      do is = 1_ik, ns
         call site_view(t, is, v)
         associate (s => cfg%sites(is))
            call req_r(v, 'sites.area_fraction', s%area_fraction, m)
            s%has_elevation  = toml_has(v, 'sites.elevation')
            s%elevation    = toml_real(v, 'sites.elevation', cfg%forcing%elevation)
            s%slope      = toml_real(v, 'sites.slope',  s%slope)
            s%aspect     = toml_real(v, 'sites.aspect', s%aspect)
            s%census_site_id = toml_int (v, 'sites.census_site_id', s%census_site_id)
            s%own_ground     = toml_has_section(v, 'soil')
            if (s%own_ground) then
               s%soil_column = cfg%soil_column
               call load_soil_column(v, s%soil_column)
               s%water_bottom_bc = cfg%soil%water_bottom_bc
               s%drain_factor    = cfg%soil%drain_factor
               call load_water_bottom_bc(v, s%water_bottom_bc, s%drain_factor)
               s%deep_temp  = toml_real(v, 'soil.deep_temp',  cfg%energy%deep_temp)
               s%deep_depth = toml_real(v, 'soil.deep_depth', cfg%energy%deep_depth)
            end if
         end associate
      end do
   end subroutine load_sites

   !----- Table `is` of [[sites]] as its own table: its [soil] keys (SITE_GROUND_KEYS) as            !
   !      soil.<key>, the rest as sites.<key>, its source the file and the table's place. -----------!
   subroutine site_view(t, is, v)
      type(toml_table_t), intent(in)  :: t
      integer(ik),        intent(in)  :: is
      type(toml_table_t), intent(out) :: v
      character(len=32) :: pre
      character(len=64) :: key
      integer(ik) :: i
      write(pre, '(a,i0,a)') 'sites.', is, '.'
      write(v%source, '(2a,i0)') trim(t%source), ' [[sites]] ', is
      allocate(v%key(t%n), v%val(t%n))
      v%n = 0_ik
      do i = 1_ik, t%n
         if (index(t%key(i), trim(pre)) /= 1) cycle
         key = t%key(i)(len_trim(pre) + 1:)
         v%n = v%n + 1_ik
         if (any(SITE_GROUND_KEYS == key)) then
            v%key(v%n) = 'soil.'//key
         else
            v%key(v%n) = 'sites.'//key
         end if
         v%val(v%n) = t%val(i)
      end do
   end subroutine site_view

   !----- [run].mode and the [region] block (docs/configuration.md, "Regional runs"). The mode      !
   !      is a DEFAULTED read ("site"). In region mode the box is required; the selection rule      !
   !      and the cell mask default. ----------------------------------------------------------------!
   subroutine load_region_config(t, cfg, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      type(keymiss_t),     intent(inout) :: m
      real(wp)    :: buf(5)
      integer(ik) :: nout
      call opt_choice(t, 'run.mode', [character(len=NAME_LEN) :: 'site', 'region'],                       &
                      [RUN_MODE_SITE, RUN_MODE_REGION], cfg%mode, m)
      if (cfg%mode /= RUN_MODE_REGION) return
      call toml_real_array(t, 'region.box_nwse', buf, nout)
      if (nout == 4_ik) then
         cfg%region%box_nwse = buf(1:4)
      else
         call note_missing(m, 'region.box_nwse')
      end if
      cfg%region%land_fraction_min = toml_real(t, 'region.land_fraction_min', 0.5_wp)
      cfg%region%cell_mask         = toml_string(t, 'region.cell_mask', '')
      cfg%region%cell_mask_var     = toml_string(t, 'region.cell_mask_var', 'mask')
      cfg%region%cell_mask_min     = toml_real(t, 'region.cell_mask_min', 0.5_wp)
   end subroutine load_region_config

   !----- Load the [output] diagnostic-aggregation block (docs/science/diagnostics.md §5).        !
   !      OPT-IN: gated on output.enabled (a DEFAULTED read, so a config with no [output]         !
   !      block writes no diagnostic stream). ALL keys are optional-with-default, so no           !
   !      keymiss is recorded. Per-variable overrides live in the optional                        !
   !      meds_io_config.toml, applied by the driver against the built registry. -----------------!
   subroutine load_output_config(t, cfg, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      type(keymiss_t),     intent(inout) :: m
      character(len=NAME_LEN), parameter :: CHUNKS(4) = [character(len=NAME_LEN) :: 'day', 'month', 'year', 'run']
      cfg%output = output_config_t()                                     ! defaults
      cfg%output%enabled = toml_logical(t, 'output.enabled', .false.)    ! opt-in gate (defaulted)
      if (.not. cfg%output%enabled) return
      cfg%output%dir                 = toml_string (t, 'output.dir',                 '.')
      cfg%output%prefix              = toml_string (t, 'output.prefix',              'meds')
      cfg%output%io_config           = toml_string (t, 'output.io_config',           '')
      cfg%output%cohort_max          = toml_int    (t, 'output.cohort_max',          4096_ik)
      cfg%output%patch_max           = toml_int    (t, 'output.patch_max',           256_ik)
      cfg%output%fast_interval_steps = toml_int    (t, 'output.fast_interval_steps', 4_ik)
      call opt_choice(t, 'output.sync_every', [character(len=NAME_LEN) :: 'flush', 'never'],               &
                      [SYNC_FLUSH, SYNC_NEVER], cfg%output%sync_every, m)
      !----- Variable-GROUP toggles: the groups that make a run judgeable are on, and the two      !
      !      heavy evaluation groups (radiation, ecophys) off. ------------------------------------!
      cfg%output%grp_on(GRP_STRUCTURE)  = toml_logical(t, 'output.structure',  .true.)
      cfg%output%grp_on(GRP_CARBON)     = toml_logical(t, 'output.carbon',     .true.)
      cfg%output%grp_on(GRP_WATER)      = toml_logical(t, 'output.water',      .true.)
      cfg%output%grp_on(GRP_ENERGY)     = toml_logical(t, 'output.energy',     .true.)
      cfg%output%grp_on(GRP_RADIATION)  = toml_logical(t, 'output.radiation',  .false.)
      cfg%output%grp_on(GRP_ECOPHYS)    = toml_logical(t, 'output.ecophys',    .false.)
      cfg%output%grp_on(GRP_BIOGEOCHEM) = toml_logical(t, 'output.biogeochem', .true.)
      cfg%output%grp_on(GRP_FORCING)    = toml_logical(t, 'output.forcing',    .true.)
      !----- AXIS toggles: suppress a whole trailing dimension without naming variables. The 2-D    !
      !      (soil layer x patch) axis is the highest-volume non-cohort axis and stays off. --------!
      cfg%output%axis_on(AXIS_COHORT)     = toml_logical(t, 'output.axes_cohort',     .true.)
      cfg%output%axis_on(AXIS_PATCH)      = toml_logical(t, 'output.axes_patch',      .true.)
      cfg%output%axis_on(AXIS_PFT)        = toml_logical(t, 'output.axes_pft',        .true.)
      cfg%output%axis_on(AXIS_SIZE)       = toml_logical(t, 'output.axes_size',       .true.)
      cfg%output%axis_on(AXIS_SOIL_PATCH) = toml_logical(t, 'output.axes_soil_patch', .false.)
      !----- DBH size-class edges (DIM_SIZE), ascending, n_class+1 of them. Absent -> manager_setup !
      !      installs the built-in inventory default. ----------------------------------------------!
      call load_dbh_edges(t, cfg%output)
      !----- Integrator WORK counters and the budget residuals; on by default. -------------------!
      cfg%output%grp_on(GRP_NUMERICS)  = toml_logical(t, 'output.numerics',      .true.)
      !----- per-tier enable + file-chunk (index order FAST/DAILY/MONTHLY/ANNUAL). ---!
      cfg%output%freq_on(1) = toml_logical(t, 'output.fast.enabled',    .false.)
      cfg%output%freq_on(2) = toml_logical(t, 'output.daily.enabled',   .true.)
      cfg%output%freq_on(3) = toml_logical(t, 'output.monthly.enabled', .true.)
      cfg%output%freq_on(4) = toml_logical(t, 'output.annual.enabled',  .true.)
      call opt_choice(t, 'output.fast.file_chunk',    CHUNKS, [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN],              &
                      cfg%output%file_chunk(1), m)
      call opt_choice(t, 'output.daily.file_chunk',   CHUNKS, [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN],              &
                      cfg%output%file_chunk(2), m)
      call opt_choice(t, 'output.monthly.file_chunk', CHUNKS, [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN],              &
                      cfg%output%file_chunk(3), m)
      call opt_choice(t, 'output.annual.file_chunk',  CHUNKS, [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN],              &
                      cfg%output%file_chunk(4), m)
   end subroutine load_output_config

   !----- Load the per-PFT [phenology] cue params from the PFT file                               !
   !      (docs/science/plant_phenology.md). Leaf phenology always runs; what is gated is the     !
   !      PRESENCE of a [phenology] section in THIS PFT file: absent -> keep the                  !
   !      alloc_pft_table defaults (no cues: always flushing, never senescing, losing leaves      !
   !      to background turnover alone, with a ~15-day flush cap). Present -> every per-PFT       !
   !      array is REQUIRED except the optional water keys (the no-silent-defaults rule still     !
   !      holds for a deliberate override). The arrays are flattened per cohort into a            !
   !      pheno_params_t by the slow-loop phenology driver.                                       !
   !                                                                                              !
   !      THE GATE IS THE SECTION, NOT A KEY INSIDE IT (toml_has_section). Gated on one key,      !
   !      a full, deliberate [phenology] block that misspells that key, or spells it under a      !
   !      retired name, is skipped in SILENCE: the presence map never sees the required keys      !
   !      go missing, and every PFT falls back to CUE_NONE/CUE_NONE -- always flushing, never     !
   !      senescing -- whatever leaf habit it declared. Gating on the section asks whether        !
   !      the author meant to configure phenology at all, which is the question. -----------------!
   subroutine load_phenology_pft(t, cfg, npft, m)
      type(toml_table_t),  intent(in)    :: t
      type(meds_config_t), intent(inout) :: cfg
      integer(ik),         intent(in)    :: npft
      type(keymiss_t),     intent(inout) :: m
      if (.not. toml_has_section(t, 'phenology')) return
      call req_pa_int(t, 'phenology.flush_cue_mask',      cfg%pft%pheno%flush_cue_mask,        npft, m)
      call req_pa_int(t, 'phenology.shed_cue_mask',       cfg%pft%pheno%shed_cue_mask,         npft, m)
      call req_pa(t, 'phenology.flush_cue_timescale',     cfg%pft%pheno%flush_cue_timescale,   npft, m)
      call req_pa(t, 'phenology.shed_cue_timescale',      cfg%pft%pheno%shed_cue_timescale,    npft, m)
      call req_pa(t, 'phenology.flush_rate_max',          cfg%pft%pheno%flush_rate_max,        npft, m)
      call req_pa(t, 'phenology.shed_rate_max',           cfg%pft%pheno%shed_rate_max,         npft, m)
      call req_pa(t, 'phenology.min_leaf_cover',          cfg%pft%pheno%min_leaf_cover,        npft, m)
      call req_pa(t, 'phenology.bare_leaf_cover',         cfg%pft%pheno%bare_leaf_cover,       npft, m)
      call req_pa(t, 'phenology.flush_base_temp',         cfg%pft%pheno%flush_base_temp,       npft, m)
      call req_pa(t, 'phenology.flush_degree_days',       cfg%pft%pheno%flush_degree_days,     npft, m)
      call req_pa(t, 'phenology.flush_temp_sharpness',    cfg%pft%pheno%flush_temp_sharpness,  npft, m)
      call req_pa(t, 'phenology.shed_base_temp',          cfg%pft%pheno%shed_base_temp,        npft, m)
      call req_pa(t, 'phenology.shed_degree_days',        cfg%pft%pheno%shed_degree_days,      npft, m)
      call req_pa(t, 'phenology.shed_temp_sharpness',     cfg%pft%pheno%shed_temp_sharpness,   npft, m)
      call req_pa(t, 'phenology.par_min',                 cfg%pft%pheno%par_min,               npft, m)
      call req_pa(t, 'phenology.flush_light_hours',       cfg%pft%pheno%flush_light_hours,     npft, m)
      call req_pa(t, 'phenology.flush_light_sharpness',   cfg%pft%pheno%flush_light_sharpness, npft, m)
      call req_pa(t, 'phenology.shed_light_hours',        cfg%pft%pheno%shed_light_hours,      npft, m)
      call req_pa(t, 'phenology.shed_light_sharpness',    cfg%pft%pheno%shed_light_sharpness,  npft, m)
      call req_pa(t, 'phenology.light_window',            cfg%pft%pheno%light_window,          npft, m)
      !----- OPTIONAL: the water cue, which a PFT turns on through its masks. An absent key keeps  !
      !      the alloc_pft_table default, which equals the pheno_params_t default. --------------!
      call opt_pa(t, 'phenology.flush_water_sum',         cfg%pft%pheno%flush_water_sum,       npft, m)
      call opt_pa(t, 'phenology.flush_water_sharpness',   cfg%pft%pheno%flush_water_sharpness, npft, m)
      call opt_pa(t, 'phenology.shed_water_sum',          cfg%pft%pheno%shed_water_sum,        npft, m)
      call opt_pa(t, 'phenology.shed_water_sharpness',    cfg%pft%pheno%shed_water_sharpness,  npft, m)
   end subroutine load_phenology_pft

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
      call string_to_time(trim(sval), out, ok)
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

   !----- OPTIONAL per-PFT array: keeps whatever `out` already holds when the key is absent,      !
   !      so a config that never mentions it gets the alloc_pft_table default. Used for           !
   !      parameters that only bite when an OPT-IN feature is selected (the phenology water       !
   !      cue, the per-PFT hydraulic overrides, storage maintenance) or that have a               !
   !      defensible default. A WRONG-LENGTH array is still an error: silently ignoring           !
   !      `[phenology] flush_water_sum = [10.0]` in a three-PFT run would be exactly the typo     !
   !      this table's required-key discipline exists to catch. ----------------------------------!
   subroutine opt_pa(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      real(wp)    :: buf(MAXPFT)
      integer(ik) :: nout
      if (.not. toml_has(t, key)) return                  ! absent -> keep the default
      call toml_real_array(t, key, buf, nout)
      if (nout == npft) then ; out(1:npft) = buf(1:npft) ; else ; call note_missing(m, key) ; end if
   end subroutine opt_pa

   !----- A required per-PFT real given as one value, for every PFT, or as an array of npft. --!
   subroutine req_pa_bcast(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      if (.not. toml_has(t, key)) then
         call note_missing(m, key)
         return
      end if
      call read_pa_bcast(t, key, out, npft, m)
   end subroutine req_pa_bcast

   !----- An optional req_pa_bcast: absent, `out` keeps what it holds. ---------------------!
   subroutine opt_pa_bcast(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      if (.not. toml_has(t, key)) return
      call read_pa_bcast(t, key, out, npft, m)
   end subroutine opt_pa_bcast

   subroutine read_pa_bcast(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      real(wp),           intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      real(wp)    :: buf(MAXPFT)
      integer(ik) :: nout
      if (toml_is_array(t, key)) then
         call toml_real_array(t, key, buf, nout)
         if (nout == npft) then
            out(1:npft) = buf(1:npft)
         else
            call note_missing(m, key)
         end if
      else
         out(1:npft) = toml_real(t, key, 0.0_wp)
      end if
   end subroutine read_pa_bcast

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

   !----- A per-PFT LOGICAL array. The minimal TOML reader has no boolean-array primitive, and    !
   !      adding one would mean a second array parser; a 0/1 real array reuses the one that exists  !
   !      and reads acceptably in the config (`is_woody = [1, 1, 0]`). Anything other than 0 is     !
   !      true, so a typo cannot silently read as false.  ------------------------------------------!
   subroutine req_pa_log(t, key, out, npft, m)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: key
      logical,            intent(inout) :: out(:)
      integer(ik),        intent(in)    :: npft
      type(keymiss_t),    intent(inout) :: m
      real(wp)    :: buf(MAXPFT)
      integer(ik) :: nout
      call toml_real_array(t, key, buf, nout)
      if (nout == npft) then
         out(1:npft) = nint(buf(1:npft), ik) /= 0_ik
      else
         call note_missing(m, key)
      end if
   end subroutine req_pa_log

   !=======================================================================================!
   !  Load configuration from `path` (the MAIN file). Both files are mandatory; a missing    !
   !  file or any missing required key is a hard error.                                       !
   !=======================================================================================!
   subroutine load_meds_config(path, cfg)
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(out) :: cfg
      type(toml_table_t) :: tm, tp
      type(keymiss_t)    :: miss
      type(key_report_t) :: unknown
      logical            :: found
      integer(ik)        :: npft, nout, i
      real(wp)           :: buf(MAXPFT)

      !----- The parameter record is this config's alone: rows from a config loaded earlier in    !
      !      the same process (the Python API runs one after another) would otherwise stay in it. !
      call toml_record_reset()

      !----- MAIN file. Every key it holds must be one its reference lists (meds_config_keys);   !
      !      the unknown and retired ones are reported with the missing required keys below. ---!
      call toml_parse_file(path, tm, found)
      if (.not. found) error stop 'meds_config: main config file not found: '//trim(path)
      call check_config_keys(tm, MAIN_KEYS, unknown)

      call req_dur (tm, 'run.dt_slow',    cfg%dt_slow,    miss)
      call req_date(tm, 'run.start_time', cfg%start_time, miss)
      call req_date(tm, 'run.end_time',   cfg%end_time,   miss)
      !----- Master slow-tier freeze and the slow ledger (both DEFAULTED true). ------------------!
      cfg%slow_on = toml_logical(tm, 'run.slow_on', .true.)
      cfg%slow_ledger_on = toml_logical(tm, 'run.slow_ledger_on', .true.)
      !----- Host threads (DEFAULTED 1 = serial). Exposed ONCE, here -- deliberately NOT also       !
      !      from OMP_NUM_THREADS, so the run's thread count is recorded in the config that         !
      !      produced it. --------------------------------------------------------------------------!
      cfg%n_threads = toml_int(tm, 'run.n_threads', 1_ik)
      !----- [run].mode (DEFAULTED "site") and, for a region, the [region] block. ----------------!
      call load_region_config(tm, cfg, miss)

      !----- Fast (sub-daily) biophysics loop. --------------------------------------------!
      call req_l     (tm, 'fast.fast_biophysics_on', cfg%fast_biophysics_on, miss)
      call req_dur   (tm, 'fast.dt_fast',            cfg%dt_fast,            miss)
      !----- Optional fast-loop settings: DEFAULTED reads, so a config without them runs the      !
      !      documented defaults. ----------------------------------------------------------------!
      cfg%snow_init_swe       = toml_real(tm, 'fast.snow_init_swe',  0.0_wp)
      cfg%snow_init_temp      = toml_real(tm, 'fast.snow_init_temp', 270.0_wp)
      cfg%canopy_water_on     = toml_logical(tm, 'fast.canopy_water_on', .false.)
      !----- Fast-loop TIME integrator: ARK (ESDIRK2) when absent (docs/science/numerical_scheme.md §3). !
      call opt_choice(tm, 'fast.time_integrator', [character(len=NAME_LEN) :: 'ark', 'rk45'],              &
                      [INTEG_ARK, INTEG_RK45], cfg%time_integrator, miss)
      cfg%ark_adaptive      = toml_logical(tm, 'fast.ark_adaptive',      .true.)
      cfg%rtol              = toml_real   (tm, 'fast.rtol',              cfg%rtol)
      !----- Error-control selectors (DEFAULTED). ---------------------------------------------!
      cfg%atol_scale        = toml_real   (tm, 'fast.atol_scale',        1.0_wp)
      cfg%forcing_sample_frac = toml_real (tm, 'fast.forcing_sample_frac', 0.5_wp)
      cfg%cas_condensation  = toml_logical(tm, 'fast.cas_condensation', .true.)
      !----- Process mask: all-on default => the full column. -----------------------------------!
      cfg%mask%veg_energy   = toml_logical(tm, 'fast.mask.veg_energy', .true.)
      cfg%mask%cas_energy   = toml_logical(tm, 'fast.mask.cas_energy', .true.)
      cfg%mask%cas_vapour   = toml_logical(tm, 'fast.mask.cas_vapour', .true.)
      cfg%mask%cas_co2      = toml_logical(tm, 'fast.mask.cas_co2',    .true.)
      cfg%mask%soil_heat    = toml_logical(tm, 'fast.mask.soil_heat',  .true.)
      cfg%mask%soil_water   = toml_logical(tm, 'fast.mask.soil_water', .true.)
      cfg%mask%hydraulics   = toml_logical(tm, 'fast.mask.hydraulics', .true.)
      call opt_choice(tm, 'fast.step_controller', [character(len=NAME_LEN) :: 'i', 'pi'],                  &
                      [CTRL_I, CTRL_PI], cfg%step_controller, miss)
      call opt_choice(tm, 'fast.error_level', [character(len=NAME_LEN) :: 'fixed', 'adaptive', 'strict'],  &
                      [CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT], cfg%error_level, miss)
      cfg%ark_dt_init       = toml_real   (tm, 'fast.ark_dt_init',       0.0_wp)
      cfg%ark_fixed_substep = toml_int    (tm, 'fast.ark_fixed_substep', 4_ik)
      cfg%ark_coupled       = toml_logical(tm, 'fast.ark_coupled',       .true.)
      cfg%fast_probe        = toml_logical(tm, 'fast.fast_probe',        .false.)
      cfg%fast_probe_file   = toml_string (tm, 'fast.fast_probe_file',   'fast_probe.csv')

      !----- [soil]/[snow]/[aerodynamics] fast-loop biophysics run-config, and its solver         !
      !      settings in [fast] (every key DEFAULTED to its type's value, so an absent block      !
      !      keeps the defaults). build_fast_context copies these verbatim into the column        !
      !      config. -----------------------------------------------------------------------------!
      call load_soil_opts  (tm, cfg%soil)
      call load_soil_column(tm, cfg%soil_column)
      call load_energy_opts(tm, cfg%energy)
      !----- ONE debug switch: [fast].debug_error arms every Debug-only hard stop, including the       !
      !      soil-water kernel's own mass-closure stops, which have no key of their own. --------------!
      cfg%soil%debug_error = cfg%energy%debug_error
      call load_snow_params(tm, cfg%snow)
      call load_aero_cfg   (tm, cfg%aero)

      !----- [soil_carbon] (docs/science/soil_carbon.md): every key optional, at its record default. -!
      call load_soil_carbon_opts(tm, cfg%soil_carbon, miss)

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
      !----- The tolerance ceiling is optional. Absent, it is 0.15, or patch_light_tol when         !
      !      that is larger, which keeps a looser config valid. ------------------------------------!
      cfg%patch_light_tol_max = 0.15_wp
      if (toml_has(tm, 'demography.patch_light_tol'))                                               &
         cfg%patch_light_tol_max = max(0.15_wp, cfg%patch_light_tol)   ! clamp-ok: the setting's default
      cfg%patch_light_tol_max = toml_real(tm, 'demography.patch_light_tol_max', cfg%patch_light_tol_max)
      cfg%light_ext_coeff_fusefiss = toml_real(tm, 'demography.light_ext_coeff_fusefiss', 0.5_wp)
      call req_r(tm, 'demography.min_patch_area',         cfg%min_patch_area,         miss)
      call req_r(tm, 'demography.conservation_tol',       cfg%conservation_tol,       miss)
      call req_r(tm, 'demography.growth_memory_days',     cfg%growth_memory_days,     miss)

      call req_r(tm, 'disturbance.patch_disturbance_rate',     cfg%patch_disturbance_rate,     miss)
      call req_r(tm, 'disturbance.disturbance_survive_height', cfg%disturbance_survive_height, miss)

      call req_r(tm, 'recruitment.min_recruit_size',      cfg%min_recruit_size,       miss)

      call req_i(tm, 'init.init_mode',     cfg%init%init_mode,         miss)
      !----- The state file and the census file are required by the mode that starts from them,     !
      !      and optional otherwise: the Python census start reads census_file in any mode. ---------!
      if (cfg%init%init_mode == INIT_RESTART) then
         call req_s(tm, 'init.restart_file', cfg%init%restart_file, miss)
      else
         cfg%init%restart_file = toml_string(tm, 'init.restart_file', '')
      end if
      if (cfg%init%init_mode == INIT_CENSUS) then
         call req_s(tm, 'init.census_file', cfg%init%census_file, miss)
      else
         cfg%init%census_file = toml_string(tm, 'init.census_file', '')
      end if
      call req_s(tm, 'init.pft_config',    cfg%init%pft_config,        miss)
      cfg%init%soil_temp  = toml_real(tm, 'init.soil_temp',  cfg%init%soil_temp)
      cfg%init%soil_theta = toml_real(tm, 'init.soil_theta', cfg%init%soil_theta)
      cfg%init%reacclimate_traits = toml_logical(tm, 'init.reacclimate_traits', cfg%init%reacclimate_traits)

      !----- [state]: the restart stream; every diagnostic is [output]'s. ------------------------!
      call req_s(tm, 'state.output_dir',     cfg%state%output_dir,         miss)
      call req_s(tm, 'state.output_prefix',  cfg%state%output_prefix,      miss)
      call req_l(tm, 'state.write_state',    cfg%state%write_state,        miss)
      call req_i(tm, 'state.interval_years', cfg%state%interval_years, miss)


      !----- Carbon growth: the stub GPP of a run without the fast loop, required there only. ---!
      if (.not. cfg%fast_biophysics_on) then
         call req_r(tm, 'carbon.gpp_ref', cfg%gpp_ref, miss)
      else
         cfg%gpp_ref = toml_real(tm, 'carbon.gpp_ref', cfg%gpp_ref)
      end if

      !----- Light trait plasticity (opt-in; DEFAULTED false). Leaf phenology always runs;        !
      !      its per-PFT cue params are read by load_phenology_pft. ------------------------------!
      cfg%trait_plasticity_on = toml_logical(tm, 'trait_dynamics.trait_plasticity_on', .false.)
      cfg%fineroot_turnover_plastic = toml_logical(tm, 'trait_dynamics.fineroot_turnover_plastic', .false.)

      !----- Meteorological forcing (opt-in; gated on forcing.forcing_on, defaulted false). !
      call load_forcing_config(tm, cfg, miss)
      !----- [maps]: the maps a run takes its soil (and its Chave site index) from. ----------!
      call load_maps_config(tm, cfg)
      !----- [[sites]]: the polygon's sites, after the [soil] they inherit and the maps that    !
      !      refuse their soil. ----------------------------------------------------------------!
      call load_sites(tm, cfg, miss)

      !----- Diagnostic-aggregation output (opt-in; gated on output.enabled, defaulted false; !
      !      a config with no [output] block emits no diagnostic stream).  ------------------!
      call load_output_config(tm, cfg, miss)

      !----- Leaf physiology: model selectors + shared biochemistry (non-PFT). ------------!
      call req_choice(tm, 'leaf_physiology.stomatal_model', [character(len=NAME_LEN) :: 'leuning', 'medlyn', 'katul'], &
                      [SM_LEUNING, SM_MEDLYN, SM_KATUL], cfg%stomatal_model, miss)
      call req_choice(tm, 'leaf_physiology.temp_response_form', [character(len=NAME_LEN) :: 'arrhenius', 'peaked'],   &
                      [TEMP_RESPONSE_ARRHENIUS, TEMP_RESPONSE_PEAKED], cfg%temp_response_form, miss)
      call req_choice(tm, 'leaf_physiology.colimitation', [character(len=NAME_LEN) :: 'min', 'quadratic'],            &
                      [COLIM_MIN, COLIM_QUADRATIC], cfg%colimitation, miss)
      call req_l(tm, 'leaf_physiology.use_boundary_layer', cfg%use_boundary_layer, miss)
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
      !----- Thermal acclimation (docs/science/leaf_gas_exchange.md eq. (4)): OPT-IN, every       !
      !      key optional, the six Kattge & Knorr coefficients defaulting to their published      !
      !      values. -----------------------------------------------------------------------------!
      cfg%thermal_acclimation = toml_logical(tm, 'leaf_physiology.thermal_acclimation',      &
                                                  cfg%thermal_acclimation)
      cfg%acclim_ds_vcmax_a  = toml_real(tm, 'leaf_physiology.acclim_ds_vcmax_a',  cfg%acclim_ds_vcmax_a)
      cfg%acclim_ds_vcmax_b  = toml_real(tm, 'leaf_physiology.acclim_ds_vcmax_b',  cfg%acclim_ds_vcmax_b)
      cfg%acclim_ds_jmax_a   = toml_real(tm, 'leaf_physiology.acclim_ds_jmax_a',   cfg%acclim_ds_jmax_a)
      cfg%acclim_ds_jmax_b   = toml_real(tm, 'leaf_physiology.acclim_ds_jmax_b',   cfg%acclim_ds_jmax_b)
      cfg%acclim_jv_a        = toml_real(tm, 'leaf_physiology.acclim_jv_a',        cfg%acclim_jv_a)
      cfg%acclim_jv_b        = toml_real(tm, 'leaf_physiology.acclim_jv_b',        cfg%acclim_jv_b)
      cfg%acclim_window_days = toml_real(tm, 'leaf_physiology.acclim_window_days', cfg%acclim_window_days)
      call req_r(tm, 'leaf_physiology.ds_vcmax', cfg%ds_vcmax, miss)
      call req_r(tm, 'leaf_physiology.ds_jmax',  cfg%ds_jmax,  miss)
      call req_r(tm, 'leaf_physiology.ds_rd',    cfg%ds_rd,    miss)
      call req_r(tm, 'leaf_physiology.o2_mol_frac',      cfg%o2_mol_frac,      miss)
      call req_r(tm, 'leaf_physiology.leaf_absorptance', cfg%leaf_absorptance, miss)
      call req_r(tm, 'leaf_physiology.phi_psii',         cfg%phi_psii,         miss)
      cfg%medlyn_vpd_min = toml_real(tm, 'leaf_physiology.medlyn_vpd_min', cfg%medlyn_vpd_min)

      !----- PFT file (named in the main file). -------------------------------------------!
      call toml_parse_file(trim(cfg%init%pft_config), tp, found)
      if (.not. found) error stop 'meds_config: PFT config file not found: '//trim(cfg%init%pft_config)
      call check_config_keys(tp, PFT_KEYS, unknown)

      !----- PFT count comes from the length of the wood_density array. -------------------!
      call toml_real_array(tp, 'pft.wood_density', buf, nout)
      if (nout < 1_ik) error stop 'meds_config: pft.wood_density missing/empty in '//trim(cfg%init%pft_config)
      if (nout > MAXPFT) error stop 'meds_config: pft.wood_density lists more than MAXPFT PFTs in '//trim(cfg%init%pft_config)
      npft = nout
      call alloc_pft_table(cfg%pft, npft)
      cfg%pft%wood_density = buf(1:npft)

      call req_pa(tp, 'pft.reproduction_investment_fraction', cfg%pft%reproduction_investment_fraction, npft, miss)
      call req_pa(tp, 'pft.recruit_carbon_efficiency',        cfg%pft%recruit_carbon_efficiency,        npft, miss)
      call opt_pa(tp, 'pft.recruit_shade_decay',              cfg%pft%recruit_shade_decay,              npft, miss)
      call opt_pa(tp, 'pft.nbg_init_density',                 cfg%pft%nbg_init_density,                 npft, miss)
      call req_pa(tp, 'pft.seed_rain_recruits',               cfg%pft%seed_rain_recruits,               npft, miss)
      call req_pa_int(tp, 'pft.include_pft',                  cfg%pft%include_pft,                      npft, miss)
      call req_r(tp, 'pft.min_cohort_height',       cfg%pft%min_cohort_height,       miss)
      call req_pa(tp, 'pft.min_repro_height',                 cfg%pft%min_repro_height,                 npft, miss)

      !----- Leaf-photosynthesis per-PFT traits. ------------------------------------------!
      call req_pa_int(tp, 'pft.photosynthetic_pathway', cfg%pft%photosynthetic_pathway, npft, miss)
      call req_pa_log(tp, 'pft.is_woody',           cfg%pft%is_woody,           npft, miss)
      call req_pa(tp,     'pft.stem_resp_factor25', cfg%pft%stem_resp_factor25, npft, miss)
      call opt_pa(tp,     'pft.stem_resp_sapwood25', cfg%pft%stem_resp_sapwood25, npft, miss)
      call opt_pa(tp,     'pft.max_relative_growth_rate', cfg%pft%max_relative_growth_rate, npft, miss)
      call opt_pa(tp,     'pft.max_absolute_growth_rate', cfg%pft%max_absolute_growth_rate, npft, miss)
      call opt_pa(tp,     'pft.max_absolute_growth_exponent', cfg%pft%max_absolute_growth_exponent, npft, miss)
      call req_pa(tp,     'pft.root_resp_factor25', cfg%pft%root_resp_factor25, npft, miss)
      !----- Canopy optics: shortwave as reflect/transmit per band, longwave as EMISSIVITY. -----!
      call req_pa(tp, 'pft.leaf_reflect_vis',  cfg%pft%leaf_reflect_vis,  npft, miss)
      call req_pa(tp, 'pft.leaf_transmit_vis', cfg%pft%leaf_transmit_vis, npft, miss)
      call req_pa(tp, 'pft.leaf_reflect_nir',  cfg%pft%leaf_reflect_nir,  npft, miss)
      call req_pa(tp, 'pft.leaf_transmit_nir', cfg%pft%leaf_transmit_nir, npft, miss)
      call req_pa(tp, 'pft.leaf_emissivity',   cfg%pft%leaf_emissivity,   npft, miss)
      call req_pa(tp, 'pft.wood_reflect_vis',  cfg%pft%wood_reflect_vis,  npft, miss)
      call req_pa(tp, 'pft.wood_transmit_vis', cfg%pft%wood_transmit_vis, npft, miss)
      call req_pa(tp, 'pft.wood_reflect_nir',  cfg%pft%wood_reflect_nir,  npft, miss)
      call req_pa(tp, 'pft.wood_transmit_nir', cfg%pft%wood_transmit_nir, npft, miss)
      call req_pa(tp, 'pft.wood_emissivity',   cfg%pft%wood_emissivity,   npft, miss)
      call req_pa(tp, 'pft.leaf_clumping',     cfg%pft%leaf_clumping,     npft, miss)
      call req_pa(tp, 'pft.wood_clumping',     cfg%pft%wood_clumping,     npft, miss)
      call req_pa(tp, 'pft.leaf_angle_mean',   cfg%pft%leaf_angle_mean,   npft, miss)
      call req_pa(tp, 'pft.leaf_angle_std',    cfg%pft%leaf_angle_std,    npft, miss)
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
      call req_pa(tp, 'pft.theta_cj_c3',      cfg%pft%theta_cj_c3,      npft, miss)
      call req_pa(tp, 'pft.theta_ip_c3',      cfg%pft%theta_ip_c3,      npft, miss)
      call req_pa(tp, 'pft.theta_cj_c4',      cfg%pft%theta_cj_c4,      npft, miss)
      call req_pa(tp, 'pft.theta_ic_c4',      cfg%pft%theta_ic_c4,      npft, miss)
      call req_pa(tp, 'pft.katul_lambda25',   cfg%pft%katul_lambda25,   npft, miss)
      call req_pa(tp, 'pft.wstress_lambda_exp',cfg%pft%wstress_lambda_exp,npft, miss)
      call req_pa(tp, 'pft.wstress_sref_stomata',cfg%pft%wstress_sref_stomata,npft, miss)
      call opt_pa(tp, 'pft.stomata_psi_onset',   cfg%pft%stomata_psi_onset,   npft, miss)
      call opt_pa(tp, 'pft.wstress_psi_open',    cfg%pft%wstress_psi_open,    npft, miss)
      call opt_pa(tp, 'pft.wstress_psi_close',   cfg%pft%wstress_psi_close,   npft, miss)
      call opt_pa(tp, 'pft.leaf_surf_water_max', cfg%pft%leaf_surf_water_max, npft, miss)
      call opt_pa(tp, 'pft.wood_surf_water_max', cfg%pft%wood_surf_water_max, npft, miss)

      !----- Carbon-allocation per-PFT traits (meds_plant_carbon_allocation). --------------!
      call req_pa(tp, 'pft.sla',                    cfg%pft%sla,                    npft, miss)
      call req_pa(tp, 'pft.leaf_width',              cfg%pft%leaf_width,              npft, miss)
      call req_pa(tp, 'pft.branch_diameter',         cfg%pft%branch_diameter,         npft, miss)
      call req_pa(tp, 'pft.root_to_leaf_ratio',     cfg%pft%root_to_leaf_ratio,     npft, miss)
      call req_pa(tp, 'pft.storage_cushion',        cfg%pft%storage_cushion,        npft, miss)
      call req_pa(tp, 'pft.growth_resp_factor',     cfg%pft%growth_resp_factor,     npft, miss)
      !----- Storage maintenance: OPTIONAL, default 0 (off). -----------------------------!
      call opt_pa(tp, 'pft.storage_turnover_rate', cfg%pft%storage_turnover_rate,  npft, miss)
      !----- Leaf and fine-root resorption: OPTIONAL, default 0.3. ----------------------!
      call opt_pa(tp, 'pft.retained_carbon_fraction', cfg%pft%retained_carbon_fraction, npft, miss)
      call req_pa(tp, 'pft.leaf_turnover_rate_toc', cfg%pft%leaf_turnover_rate_toc, npft, miss)
      call req_pa(tp, 'pft.fineroot_turnover_rate', cfg%pft%fineroot_turnover_rate, npft, miss)
      call req_pa(tp, 'pft.f_labile_leaf',          cfg%pft%f_labile_leaf,          npft, miss)
      call req_pa(tp, 'pft.f_labile_stem',          cfg%pft%f_labile_stem,          npft, miss)
      call req_pa(tp, 'pft.struct_lignin_frac',     cfg%pft%struct_lignin_frac,     npft, miss)

      !----- Leaf-phenology per-PFT cue params (required only if this PFT file has a          !
      !      [phenology] section; see load_phenology_pft). ------------------------------------!
      call load_phenology_pft(tp, cfg, npft, miss)

      !----- The mortality hazard and the light-plasticity slopes, per PFT. -------------------!
      call req_pa(tp, 'pft.mort_gamma',             cfg%pft%mort_gamma,             npft, miss)
      call req_pa(tp, 'pft.mort_alpha',             cfg%pft%mort_alpha,             npft, miss)
      call req_pa(tp, 'pft.mort_beta',              cfg%pft%mort_beta,              npft, miss)
      call req_pa(tp, 'pft.kplastic_vcmax',         cfg%pft%kplastic_vcmax,         npft, miss)
      call req_pa(tp, 'pft.kplastic_rd',            cfg%pft%kplastic_rd,            npft, miss)
      call req_pa(tp, 'pft.kplastic_sla',           cfg%pft%kplastic_sla,           npft, miss)
      call req_pa(tp, 'pft.kplastic_leaf_turnover', cfg%pft%kplastic_leaf_turnover, npft, miss)

      !----- Each PFT's hydraulic traits: every key one value for all PFTs, or one per PFT. ---!
      associate (hydro => cfg%pft%hydro)
         call req_pa_bcast(tp, 'hydraulics.leaf_pi0',             hydro%leaf_curve%pi0,           npft, miss)
         call req_pa_bcast(tp, 'hydraulics.leaf_elastic_mod',     hydro%leaf_curve%elastic_mod,   npft, miss)
         call req_pa_bcast(tp, 'hydraulics.leaf_apoplast_frac',   hydro%leaf_curve%apoplast_frac, npft, miss)
         call req_pa_bcast(tp, 'hydraulics.leaf_water_sat',       hydro%leaf_curve%water_sat,     npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_pi0',             hydro%wood_curve%pi0,           npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_elastic_mod',     hydro%wood_curve%elastic_mod,   npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_apoplast_frac',   hydro%wood_curve%apoplast_frac, npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_water_sat',       hydro%wood_curve%water_sat,     npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_psi50',           hydro%wood_psi50,               npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_kexp',            hydro%wood_kexp,                npft, miss)
         call req_pa_bcast(tp, 'hydraulics.wood_kmax',            hydro%wood_kmax,                npft, miss)
         call req_pa_bcast(tp, 'hydraulics.vessel_curl',          hydro%vessel_curl,              npft, miss)
         call req_pa_bcast(tp, 'hydraulics.root_beta',            hydro%root_beta,                npft, miss)
         call req_pa_bcast(tp, 'hydraulics.root_depth_b1',        hydro%root_depth_b1,            npft, miss)
         call req_pa_bcast(tp, 'hydraulics.root_depth_b2',        hydro%root_depth_b2,            npft, miss)
         call req_pa_bcast(tp, 'hydraulics.specific_root_length', hydro%specific_root_length,     npft, miss)
         call req_pa_bcast(tp, 'hydraulics.fine_root_radius',     hydro%fine_root_radius,         npft, miss)
      end associate

      !----- Each PFT's allometry: every key one value for all PFTs, or one per PFT. ---------!
      associate (allom => cfg%pft%allom)
         call req_pa_bcast(tp, 'allometry.hgt_max',          allom%hgt_max,          npft, miss)
         call req_pa_bcast(tp, 'allometry.b1Ht',             allom%b1Ht,             npft, miss)
         call req_pa_bcast(tp, 'allometry.b2Ht',             allom%b2Ht,             npft, miss)
         call req_pa_bcast(tp, 'allometry.agb_c1',           allom%agb_c1,           npft, miss)
         call req_pa_bcast(tp, 'allometry.agb_c2',           allom%agb_c2,           npft, miss)
         call req_pa_bcast(tp, 'allometry.lai_b1',           allom%lai_b1,           npft, miss)
         call req_pa_bcast(tp, 'allometry.lai_b2',           allom%lai_b2,           npft, miss)
         call req_pa_bcast(tp, 'allometry.wai_b1',           allom%wai_b1,           npft, miss)
         call req_pa_bcast(tp, 'allometry.wai_b2',           allom%wai_b2,           npft, miss)
         call req_pa_bcast(tp, 'allometry.sapwood_area_b1',  allom%sapwood_area_b1,  npft, miss)
         call req_pa_bcast(tp, 'allometry.sapwood_area_b2',  allom%sapwood_area_b2,  npft, miss)
         call req_pa_bcast(tp, 'allometry.aboveground_frac', allom%aboveground_frac, npft, miss)
      end associate
      call opt_height_allometry(tp, cfg, npft, miss)

      !----- Key check: abort with the full list of the missing required keys and of the keys    !
      !      MEDS does not read (unknown or retired), from both files at once. -----------------!
      if (miss%n > 0_ik) then
         write(*,'(a,i0,a)') ' meds_config: ', miss%n, ' required parameter(s) missing or unreadable:'
         do i = 1_ik, miss%n
            write(*,'(3a)') '   - ', trim(miss%key(i)), ''
         end do
      end if
      call check_key_use(tm, tp, cfg, unknown)
      if (unknown%n > 0_ik) then
         write(*,'(a,i0,a)') ' meds_config: ', unknown%n, ' key(s) MEDS does not read here:'
         do i = 1_ik, unknown%n
            write(*,'(2a)') '   - ', trim(unknown%line(i))
         end do
         write(*,'(a)') '   Every key MEDS reads is listed in meds_config_main.toml and meds_config_pft.toml;'
         write(*,'(a)') '   CHANGELOG.md lists the keys that were renamed or retired.'
      end if
      if (miss%n > 0_ik .or. unknown%n > 0_ik)                                                     &
         error stop 'meds_config: the configuration has missing, unknown or inapplicable keys (see above)'

      !----- Install allometry + compute every derived quantity (one consolidation point). -----!
      call derive_parameters(cfg)

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
      write(u,'(a)') 'pft,wood_density,hgt_max,'                                                    &
           //'reproduction_investment_fraction,recruit_carbon_efficiency,'                            &
           //'mort_gamma,mort_alpha,mort_beta,seed_rain_recruits,include_pft,'                         &
           //'min_cohort_height,min_repro_height,'                                         &
           //'photosynthetic_pathway,vcmax25,jmax25,tpu25,rd25,kp25,'                                  &
           //'stomatal_g0,stomatal_g1,stomatal_d0,quantum_yield_c4,theta_j,theta_cj_c3,theta_ip_c3,'   &
           //'theta_cj_c4,theta_ic_c4,'                                                               &
           //'katul_lambda25,wstress_psi_open,wstress_psi_close,wstress_lambda_exp,wstress_sref_stomata,' &
           //'stomata_psi_onset,leaf_surf_water_max,wood_surf_water_max,'                                &
           //'sla,root_to_leaf_ratio,aboveground_frac,storage_cushion,growth_resp_factor,'      &
           //'storage_turnover_rate,retained_carbon_fraction,'                                     &
           //'leaf_turnover_rate_toc,fineroot_turnover_rate,'                                &
           //'f_labile_leaf,f_labile_stem,struct_lignin_frac,recruit_shade_decay,'                  &
           //'kplastic_vcmax,kplastic_rd,kplastic_sla,kplastic_leaf_turnover,'                        &
           //'stem_resp_factor25,stem_resp_sapwood25,root_resp_factor25,max_relative_growth_rate,'      &
           //'max_absolute_growth_rate,max_absolute_growth_exponent,nbg_init_density,'                &
           //'height_allometry,b1Ht,b2Ht,chave_a,chave_b,chave_c,chave_rse,chave_e,'                          &
           //'agb_c1,agb_c2,lai_b1,lai_b2,wai_b1,wai_b2,sapwood_area_b1,sapwood_area_b2,'            &
           //'leaf_pi0,leaf_elastic_mod,leaf_apoplast_frac,leaf_water_sat,wood_pi0,wood_elastic_mod,'   &
           //'wood_apoplast_frac,wood_water_sat,wood_psi50,wood_kexp,wood_kmax,vessel_curl,root_beta,'  &
           //'root_depth_b1,root_depth_b2,specific_root_length,fine_root_radius'
      associate (p => cfg%pft)
         do pf = 1_ik, p%n
            !----- 92 ITEMS: i0 + 8 + i0 + 2 + i0 + 46 + i0 + 32. A format SHORTER than the value    !
            !      list does not fail -- Fortran reverts and re-uses the last repeat group, so an    !
            !      integer slot silently receives a real and prints its bit pattern, and the trailing !
            !      columns vanish. Keep the count here in step with both the header and the list. ---!
            write(u,'(i0,8(",",es15.8),",",i0,2(",",es15.8),",",i0,46(",",es15.8),",",i0,32(",",es15.8))') &
                 pf, p%wood_density(pf), p%allom(pf)%hgt_max,                                           &
                 p%reproduction_investment_fraction(pf), p%recruit_carbon_efficiency(pf),                &
                 p%mort_gamma(pf), p%mort_alpha(pf), p%mort_beta(pf), p%seed_rain_recruits(pf),         &
                 p%include_pft(pf), p%min_cohort_height, p%min_repro_height(pf),                        &
                 p%photosynthetic_pathway(pf), p%vcmax25(pf), p%jmax25(pf), p%tpu25(pf),                &
                 p%rd25(pf), p%kp25(pf), p%stomatal_g0(pf), p%stomatal_g1(pf), p%stomatal_d0(pf),       &
                 p%quantum_yield_c4(pf), p%theta_j(pf), p%theta_cj_c3(pf), p%theta_ip_c3(pf),           &
                 p%theta_cj_c4(pf), p%theta_ic_c4(pf),                                                  &
                 p%katul_lambda25(pf), pft_wstress_psi_open(cfg, pf), pft_wstress_psi_close(cfg, pf),   &
                 p%wstress_lambda_exp(pf), p%wstress_sref_stomata(pf), pft_stomata_psi_onset(cfg, pf),   &
                 p%leaf_surf_water_max(pf), p%wood_surf_water_max(pf),                                &
                 p%sla(pf), p%root_to_leaf_ratio(pf), p%allom(pf)%aboveground_frac,                     &
                 p%storage_cushion(pf), p%growth_resp_factor(pf), p%storage_turnover_rate(pf),        &
                 p%retained_carbon_fraction(pf),                                                     &
                 p%leaf_turnover_rate_toc(pf),                                                       &
                 p%fineroot_turnover_rate(pf),                                                        &
                 p%f_labile_leaf(pf), p%f_labile_stem(pf), p%struct_lignin_frac(pf),                  &
                 p%recruit_shade_decay(pf), p%kplastic_vcmax(pf), p%kplastic_rd(pf),                   &
                 p%kplastic_sla(pf), p%kplastic_leaf_turnover(pf),                                   &
                 p%stem_resp_factor25(pf), p%stem_resp_sapwood25(pf), p%root_resp_factor25(pf),    &
                 p%max_relative_growth_rate(pf), p%max_absolute_growth_rate(pf),                  &
                 p%max_absolute_growth_exponent(pf), p%nbg_init_density(pf),                        &
                 p%allom(pf)%height_allometry, p%allom(pf)%b1Ht, p%allom(pf)%b2Ht,                         &
                 p%allom(pf)%chave_a, p%allom(pf)%chave_b, p%allom(pf)%chave_c,                       &
                 p%allom(pf)%chave_rse, p%allom(pf)%chave_e,                                          &
                 p%allom(pf)%agb_c1, p%allom(pf)%agb_c2, p%allom(pf)%lai_b1, p%allom(pf)%lai_b2,      &
                 p%allom(pf)%wai_b1, p%allom(pf)%wai_b2,                                              &
                 p%allom(pf)%sapwood_area_b1, p%allom(pf)%sapwood_area_b2,                            &
                 p%hydro(pf)%leaf_curve%pi0, p%hydro(pf)%leaf_curve%elastic_mod,                      &
                 p%hydro(pf)%leaf_curve%apoplast_frac, p%hydro(pf)%leaf_curve%water_sat,              &
                 p%hydro(pf)%wood_curve%pi0, p%hydro(pf)%wood_curve%elastic_mod,                      &
                 p%hydro(pf)%wood_curve%apoplast_frac, p%hydro(pf)%wood_curve%water_sat,              &
                 p%hydro(pf)%wood_psi50, p%hydro(pf)%wood_kexp, p%hydro(pf)%wood_kmax,                &
                 p%hydro(pf)%vessel_curl, p%hydro(pf)%root_beta, p%hydro(pf)%root_depth_b1,           &
                 p%hydro(pf)%root_depth_b2, p%hydro(pf)%specific_root_length,                         &
                 p%hydro(pf)%fine_root_radius
         end do
      end associate
      close(u)
      write(*,'(2a)') ' pft   : ', trim(path)
   end subroutine write_pft_params_csv

   !----- Read [output].dbh_class_edges into the config. An explicit list must have at least two  !
   !      entries (one class); anything shorter is a config error rather than a silent fallback,    !
   !      because a user who typed the key meant to control the binning.  ------------------------!
   subroutine load_dbh_edges(t, out_cfg)
      type(toml_table_t),    intent(in)    :: t
      type(output_config_t), intent(inout) :: out_cfg
      real(wp)    :: tmp(MAX_DBH_EDGE)
      integer(ik) :: nout
      out_cfg%n_dbh_class = 0_ik
      if (.not. toml_has(t, 'output.dbh_class_edges')) return
      call toml_real_array(t, 'output.dbh_class_edges', tmp, nout)
      if (nout < 2_ik)                                                                           &
         error stop 'meds_config_io: [output].dbh_class_edges needs at least 2 edges (1 class)'
      out_cfg%n_dbh_class        = nout - 1_ik
      out_cfg%dbh_class_edges(1:nout)  = tmp(1:nout)
   end subroutine load_dbh_edges

end module meds_config_io
