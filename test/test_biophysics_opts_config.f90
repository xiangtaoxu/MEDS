! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_biophysics_opts_config -- unit tests for the [soil]/[fast]/[snow]/[aerodynamics] TOML    !
! wiring of the fast-loop biophysics run-config (meds_config_io load_* helpers -> the            !
! meds_biophysics_opts bundles). Writes a small TOML fixture, parses it, and checks: (1) string    !
! selectors map to the right enum codes, (2) scalar keys override the defaults, (3) keys ABSENT    !
! from the block keep their meds_biophysics_opts default (opt-in, no clobber).                    !
!==========================================================================================!
program test_biophysics_opts_config
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t,        &
                                    SOIL_BC_AQUIFER, SOIL_BC_DRAIN,                                &
                                    ENERGY_BC_GEOTHERMAL, ENERGY_BC_DIRICHLET
   use meds_toml,            only : toml_table_t, toml_parse_file, toml_write_record
   use meds_config_io,       only : load_soil_opts, load_energy_opts, load_snow_params, load_aero_cfg
   implicit none

   character(len=*), parameter :: TOMLFILE = 'test_biophysics_opts_tmp.toml'
   type(toml_table_t) :: tm
   type(soil_opts_t)    :: s
   type(energy_opts_t)  :: e
   type(snow_params_t)  :: sn
   type(aero_cfg_t)     :: a
   integer(ik) :: u
   logical     :: ok


   !----- Write a fixture that overrides a SUBSET of each block (leaving others to default). ----!
   open(newunit=u, file=TOMLFILE, status='replace', action='write')
   write(u,'(a)') '[soil]'
   write(u,'(a)') 'water_bottom_bc   = "aquifer"'
   write(u,'(a)') 'w_pond_max        = 9.0'
   write(u,'(a)') 'ground_albedo_nir = 0.25'
   write(u,'(a)') 'thermal_bottom_bc = "dirichlet"'
   write(u,'(a)') 'deep_temp         = 283.24'
   write(u,'(a)') ''
   write(u,'(a)') '[fast]'
   write(u,'(a)') 'water_rtol              = 2.5e-4'
   write(u,'(a)') 'soil_energy_budget_atol = 0.05'
   write(u,'(a)') ''
   write(u,'(a)') '[snow]'
   write(u,'(a)') 'rho_snow   = 300.0'
   write(u,'(a)') 'snow_emiss = 0.95'
   write(u,'(a)') ''
   write(u,'(a)') '[aerodynamics]'
   write(u,'(a)') 'z0m_ratio = 0.10'
   write(u,'(a)') 'n_iter_mo = 6'
   write(u,'(a)') 'z0m_ration = 0.20'             ! a misspelling: the parameter record must not list it
   close(u)

   call toml_parse_file(TOMLFILE, tm, ok)
   call check_true('TOML fixture parsed', ok, 0.0_wp)

   call load_soil_opts(tm, s)
   call load_energy_opts(tm, e)
   call load_snow_params(tm, sn)
   call load_aero_cfg(tm, a)

   print '(a)', 'test_biophysics_opts_config:'

   !----- [soil] and [fast]: string selectors -> enums, scalar overrides, and absent-keeps-default. -!
   call check_true('soil.water_bottom_bc "aquifer" -> SOIL_BC_AQUIFER', s%bottom_bc == SOIL_BC_AQUIFER,   &
                   real(s%bottom_bc, wp))
   call check('fast.water_rtol overridden', s%rtol,     2.5e-4_wp, 1.0e-12_wp)
   call check('soil.w_pond_max overridden', s%w_pond_max, 9.0_wp, 0.0_wp)
   call check('soil.ground_albedo_nir overridden', s%ground_albedo_nir, 0.25_wp, 1.0e-12_wp)
   call check('soil.ground_albedo_vis absent -> default 0.15', s%ground_albedo_vis, 0.15_wp, 0.0_wp)
   call check('soil.ground_emissivity absent -> default 0.95', s%ground_emissivity, 0.95_wp, 0.0_wp)
   call check('soil.dsl_dmax absent -> default 0.015', s%dsl_dmax, 0.015_wp, 0.0_wp)
   call check('fast.water_atol absent -> default 1e-4',  s%atol,   1.0e-4_wp, 1.0e-12_wp)

   !----- The soil heat's budget-closure threshold, for the debug halt. The soil-thermal solve is   !
   !      one backward-Euler step, so there are no substep keys. ---------------------------------------!
   call check('fast.soil_energy_budget_atol overridden', e%atol, 0.05_wp, 1.0e-12_wp)
   !----- The [soil] thermal bottom BC: selector -> enum, anchor temperature read, and the         !
   !      anchor DEPTH left at its derived default when the key is absent. -------------------------!
   call check_true('soil.thermal_bottom_bc "dirichlet" -> ENERGY_BC_DIRICHLET',                    &
                   e%bottom_bc == ENERGY_BC_DIRICHLET, real(e%bottom_bc, wp))
   call check('soil.deep_temp overridden', e%deep_temp, 283.24_wp, 1.0e-9_wp)
   call check('soil.deep_depth absent -> derived default 3.145', e%deep_depth, 3.145_wp, 1.0e-12_wp)

   !----- [snow]: overrides + default. ---------------------------------------------------------!
   call check('snow.rho_snow overridden',   sn%rho_snow,   300.0_wp, 1.0e-9_wp)
   call check('snow.snow_emiss overridden', sn%snow_emiss, 0.95_wp,  1.0e-12_wp)
   call check('snow.k_snow absent -> default', sn%k_snow,  0.15_wp,  1.0e-12_wp)

   !----- [aerodynamics]: override (real + int) + default. -------------------------------------!
   call check('aero.z0m_ratio overridden',   a%z0m_ratio, 0.10_wp, 1.0e-12_wp)
   call check('aero.n_iter_mo overridden',   real(a%n_iter_mo, wp), 6.0_wp, 0.0_wp)
   call check('aero.d_ratio absent -> default', a%d_ratio, 0.63_wp, 1.0e-12_wp)

   !----- Opt-in no-op: loading on a table WITHOUT the blocks leaves a fresh struct at defaults. !
   block
      type(soil_opts_t) :: s0
      type(energy_opts_t) :: e0
      call load_soil_opts(tm, s0)     ! same table, fresh struct: same result as s (idempotent)
      call check('reload fast.water_rtol matches', s0%rtol, s%rtol, 0.0_wp)
      call load_energy_opts(tm, e0)
      call check_true('reload soil.thermal_bottom_bc matches', e0%bottom_bc == e%bottom_bc, real(e0%bottom_bc, wp))
   end block

   !----- A table with NO [soil] block keeps the compiled default bottom boundary, geothermal. -----!
   block
      type(toml_table_t)  :: tm2
      type(energy_opts_t) :: e2
      logical             :: ok2
      open(newunit=u, file='test_biophysics_opts_empty.toml', status='replace', action='write')
      write(u,'(a)') '[run]'
      write(u,'(a)') 'dt_slow = "1d"'
      close(u)
      call toml_parse_file('test_biophysics_opts_empty.toml', tm2, ok2)
      call check_true('empty-block TOML parsed', ok2, 0.0_wp)
      call load_energy_opts(tm2, e2)
      call check_true('no [soil] block -> default ENERGY_BC_GEOTHERMAL',                           &
                      e2%bottom_bc == ENERGY_BC_GEOTHERMAL, real(e2%bottom_bc, wp))
      open(newunit=u, file='test_biophysics_opts_empty.toml', status='old', action='write')
      close(u, status='delete')
   end block

   !----- The drainage boundary's two named ends (docs/science/soil_biophysics.md §Boundary         !
   !      conditions): free_drain is drainage at f = 1, the default, and bedrock at f = 0. ------------!
   block
      type(toml_table_t) :: tb
      type(soil_opts_t)  :: s_free, s_rock, s_default
      logical            :: okb
      character(len=*), parameter :: NAMES(2) = [character(len=10) :: 'free_drain', 'bedrock']
      integer(ik) :: i
      do i = 1_ik, 2_ik
         open(newunit=u, file='test_biophysics_opts_bc.toml', status='replace', action='write')
         write(u,'(a)') '[soil]'
         write(u,'(3a)') 'water_bottom_bc = "', trim(NAMES(i)), '"'
         close(u)
         call toml_parse_file('test_biophysics_opts_bc.toml', tb, okb)
         if (i == 1_ik) call load_soil_opts(tb, s_free)
         if (i == 2_ik) call load_soil_opts(tb, s_rock)
      end do
      open(newunit=u, file='test_biophysics_opts_bc.toml', status='old', action='write')
      close(u, status='delete')
      call check_true('soil.water_bottom_bc "free_drain" -> drainage at f = 1, the default',          &
                      s_free%bottom_bc == SOIL_BC_DRAIN .and. s_free%drain_factor == 1.0_wp .and.     &
                      s_default%bottom_bc == SOIL_BC_DRAIN .and. s_default%drain_factor == 1.0_wp)
      call check_true('soil.water_bottom_bc "bedrock" -> drainage at f = 0',                          &
                      s_rock%bottom_bc == SOIL_BC_DRAIN .and. s_rock%drain_factor == 0.0_wp)
   end block

   !----- The parameter record: every key the loaders read, set or defaulted, and no key they did   !
   !      not read. A calibration trial checks the keys it set against it. ------------------------!
   block
      character(len=*), parameter :: RECFILE = 'test_parameter_record_tmp.csv'
      character(len=1024) :: line, z0m_line
      integer :: ios
      logical :: set_seen, default_seen, typo_seen, header_ok
      call toml_write_record(RECFILE)
      set_seen = .false. ; default_seen = .false. ; typo_seen = .false. ; z0m_line = '(none)'
      open(newunit=u, file=RECFILE, status='old', action='read')
      read(u, '(a)') line
      header_ok = trim(line) == 'source,key,index,present,value'
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, ',aerodynamics.z0m_ratio,') > 0) z0m_line = line
         if (index(line, ',aerodynamics.z0m_ratio,0,true,') > 0 .and. index(line, '1.0000000000000001E-001') > 0) &
            set_seen = .true.
         if (index(line, ',aerodynamics.d_ratio,0,false,') > 0) default_seen = .true.
         if (index(line, 'z0m_ration') > 0) typo_seen = .true.
      end do
      close(u, status='delete')
      call check_true('record: header', header_ok, 0.0_wp)
      if (.not. set_seen) write(*,'(2a)') '   z0m_ratio row: ', trim(z0m_line)
      call check_true('record: a key set in the file, with the value used to the last digit', set_seen, 0.0_wp)
      call check_true('record: a key left to its default is listed as a default', default_seen, 0.0_wp)
      call check_true('record: a misspelt key that nothing reads is absent', .not. typo_seen, 0.0_wp)
   end block

   !----- Clean up the fixture. ----------------------------------------------------------------!
   open(newunit=u, file=TOMLFILE, status='old', action='write')
   close(u, status='delete')

   call test_report('test_biophysics_opts_config')

contains



end program test_biophysics_opts_config
