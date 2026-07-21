!==========================================================================================!
! test_biophysics_opts_config -- unit tests for the [soil]/[energy]/[snow]/[aerodynamics] TOML  !
! wiring of the fast-loop biophysics run-config (meds_config_io load_* helpers -> the            !
! meds_biophysics_opts bundles). Writes a small TOML fixture, parses it, and checks: (1) string    !
! selectors map to the right enum codes, (2) scalar keys override the defaults, (3) keys ABSENT    !
! from the block keep their meds_biophysics_opts default (opt-in, no clobber). House style.        !
!==========================================================================================!
program test_biophysics_opts_config
   use meds_kinds,           only : wp, ik
   use meds_biophysics_opts, only : soil_opts_t, energy_opts_t, snow_params_t, aero_cfg_t,        &
                                    SOIL_BC_AQUIFER, SOIL_LIN_PICARD, SOIL_SUBSTEP_ADAPTIVE,       &
                                    ENERGY_PHASE_ON, ENERGY_PHASE_OFF
   use meds_toml,            only : toml_table_t, toml_parse_file
   use meds_config_io,       only : load_soil_opts, load_energy_opts, load_snow_params, load_aero_cfg
   implicit none

   character(len=*), parameter :: TOMLFILE = 'test_biophysics_opts_tmp.toml'
   type(toml_table_t) :: tm
   type(soil_opts_t)    :: s
   type(energy_opts_t)  :: e
   type(snow_params_t)  :: sn
   type(aero_cfg_t)     :: a
   integer(ik) :: nfail, u
   logical     :: ok

   nfail = 0_ik

   !----- Write a fixture that overrides a SUBSET of each block (leaving others to default). ----!
   open(newunit=u, file=TOMLFILE, status='replace', action='write')
   write(u,'(a)') '[soil]'
   write(u,'(a)') 'bottom_bc  = "aquifer"'
   write(u,'(a)') 'linearize  = "picard"'
   write(u,'(a)') 'rtol       = 2.5e-4'
   write(u,'(a)') 'max_picard = 9'
   write(u,'(a)') 'psi_wilt   = -160.0'
   write(u,'(a)') ''
   write(u,'(a)') '[energy]'
   write(u,'(a)') 'phase_change = "on"'
   write(u,'(a)') 'atol         = 0.05'
   write(u,'(a)') ''
   write(u,'(a)') '[snow]'
   write(u,'(a)') 'rho_snow   = 300.0'
   write(u,'(a)') 'snow_emiss = 0.95'
   write(u,'(a)') ''
   write(u,'(a)') '[aerodynamics]'
   write(u,'(a)') 'z0m_ratio = 0.10'
   write(u,'(a)') 'n_iter_mo = 6'
   close(u)

   call toml_parse_file(TOMLFILE, tm, ok)
   call check_true('TOML fixture parsed', ok, 0.0_wp)

   call load_soil_opts(tm, s)
   call load_energy_opts(tm, e)
   call load_snow_params(tm, sn)
   call load_aero_cfg(tm, a)

   print '(a)', 'test_biophysics_opts_config:'

   !----- [soil]: string selectors -> enums, scalar overrides, and absent-keeps-default. -------!
   call check_true('soil.bottom_bc "aquifer" -> SOIL_BC_AQUIFER', s%bottom_bc == SOIL_BC_AQUIFER, real(s%bottom_bc, wp))
   call check_true('soil.linearize "picard" -> SOIL_LIN_PICARD',  s%linearize == SOIL_LIN_PICARD, real(s%linearize, wp))
   call check('soil.rtol overridden',       s%rtol,     2.5e-4_wp, 1.0e-12_wp)
   call check('soil.max_picard overridden', real(s%max_picard, wp), 9.0_wp, 0.0_wp)
   call check('soil.psi_wilt overridden',   s%psi_wilt, -160.0_wp, 1.0e-9_wp)
   call check_true('soil.substep absent -> default ADAPTIVE', s%substep == SOIL_SUBSTEP_ADAPTIVE, real(s%substep, wp))
   call check('soil.h_init absent -> default 900', s%h_init, 900.0_wp, 1.0e-9_wp)
   call check('soil.atol absent -> default 1e-4',  s%atol,   1.0e-4_wp, 1.0e-12_wp)

   !----- [energy]: phase_change toggle + override + default. -----------------------------------!
   call check_true('energy.phase_change "on" -> ENERGY_PHASE_ON', e%phase_change == ENERGY_PHASE_ON, real(e%phase_change, wp))
   call check('energy.atol overridden',       e%atol, 0.05_wp,    1.0e-12_wp)
   call check('energy.rtol absent -> default', e%rtol, 1.0e-3_wp, 1.0e-12_wp)

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
      call check('reload soil.rtol matches', s0%rtol, s%rtol, 0.0_wp)
      call load_energy_opts(tm, e0)
      call check_true('reload energy.phase_change matches', e0%phase_change == e%phase_change, 0.0_wp)
   end block

   !----- Clean up the fixture. ----------------------------------------------------------------!
   open(newunit=u, file=TOMLFILE, status='old', action='write')
   close(u, status='delete')

   if (nfail == 0_ik) then
      print '(a)', 'test_biophysics_opts_config: ALL PASSED'
   else
      print '(a,i0,a)', 'test_biophysics_opts_config: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik ; print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

end program test_biophysics_opts_config
