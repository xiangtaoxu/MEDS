!==========================================================================================!
! test_plant_capi -- COVERAGE FOR THE C-API SHIM (issue #100).                                    !
!                                                                                          !
! WHY THIS TEST EXISTS. `meds_plant_capi.f90` used to be compiled ONLY into the optional shared    !
! library (`-DMEDS_BUILD_PYLIB=ON`), which is in neither the default build nor ctest. So a change   !
! to `leaf_env_t` or `leaf_photo_params_t` could pass the whole suite on BOTH back ends while       !
! leaving the C API unable to compile -- and that is exactly what happened: #95 inserted `psi_tlp`  !
! and `stress_arrestor` into the middle of `leaf_photo_params_t`, the POSITIONAL structure          !
! constructor in `to_params` shifted every trailing value by three, and the file stopped compiling  !
! ("Omitted component is not initialized" for o2_mol_frac / absorptance / phi_psii). It sat on main !
! through a fully green suite and was found only by hand-running the Python API.                     !
!                                                                                          !
! Registering this test COMPILES the shim in every build, which is the fix -- the constructors are  !
! keyword-form now, so the remaining failure mode is a compile error, and this makes that error     !
! impossible to miss. The assertions below then cover the part a compile cannot: that the bind(c)   !
! structs actually map onto the Fortran types field-for-field.                                       !
!                                                                                          !
! Deliberately does NOT need the .so, Python, or MEDS_BUILD_PYLIB -- it calls the bind(c) procedures !
! directly as Fortran.                                                                              !
!==========================================================================================!
program test_plant_capi
   use, intrinsic :: iso_c_binding, only : c_double, c_int
   use meds_kinds,        only : wp, ik
   use meds_plant_capi,   only : leaf_env_c, leaf_params_c, leaf_flux_c, meds_leaf_solve,        &
                                 meds_arrhenius, meds_peaked_arrhenius
   use meds_plant_types,  only : PATH_C3
   use meds_config,       only : SM_MEDLYN, TRESP_PEAKED, COLIM_QUADRATIC
   use meds_test_support, only : check, check_close, banner
   implicit none

   type(leaf_env_c)    :: env
   type(leaf_params_c) :: p
   type(leaf_flux_c)   :: f0, f1, f2, fb
   real(wp)            :: k25, arr

   call banner('plant C-API shim (issue #100)')

   p = c3_defaults()

   !=== 1. The shim solves at all -- i.e. to_env/to_params bound every component. ===========!
   env = make_env(0.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TRESP_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f0)
   call check(f0%converged /= 0_c_int, 'C-API solve converges on well-watered C3 defaults')
   call check(f0%A_net > 0.0_wp,  'C-API A_net > 0 in full light')
   call check(f0%gs    > 0.0_wp,  'C-API gs > 0 when unstressed')
   call check(f0%ci > 0.0_wp .and. f0%ci < 400.0_wp, 'C-API ci is between 0 and ca')

   !=== 2. `psi` (renamed from psi_soil, #99) actually reaches beta_stomata. ================!
   !      This is the ABI assertion a compile cannot make: if leaf_env_c's components were         !
   !      mis-bound, `psi` would land on some other field and gs would not respond.                !
   env = make_env(-1.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TRESP_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f1)
   env = make_env(-3.0_c_double)
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TRESP_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 0_c_int, f2)
   call check(f1%gs < f0%gs, 'C-API gs falls when psi drops 0 -> -1 MPa (beta_stomata is live)')
   call check(f2%gs < f1%gs, 'C-API gs falls further at psi = -3 MPa (monotone in psi)')
   call check(f2%A_net < f0%A_net, 'C-API A_net falls with water stress')

   !=== 3. A scalar entry point, to catch a by-value/by-reference regression. ===============!
   !      meds_arrhenius takes its arguments BY VALUE; a binding slip here shows up as a wrong      !
   !      number rather than a link error, so assert against the closed form at the reference        !
   !      temperature, where the response must be exactly 1x k25. ---------------------------------!
   k25 = 60.0_wp
   arr = real(meds_arrhenius(real(k25, c_double), 65330.0_c_double, 298.15_c_double), wp)
   call check_close(arr, k25, 1.0e-9_wp, 'C-API meds_arrhenius returns k25 at the 25 C reference')

   !=== 4. The boundary-layer switch is an int flag, not a logical -- assert it is honoured. =!
   env = make_env(0.0_c_double) ; env%gb = 0.5_c_double
   call meds_leaf_solve(env, p, int(SM_MEDLYN, c_int), int(TRESP_PEAKED, c_int),                 &
                        int(COLIM_QUADRATIC, c_int), 1_c_int, fb)
   call check(fb%converged /= 0_c_int, 'C-API solve converges with the boundary layer ON')
   call check(fb%cs < 400.0_wp, 'C-API boundary layer draws cs below ca (flag reached the kernel)')

   write(*,'(a)') '   PASS'

contains

   !----- Env at full light / 25 C, parameterised only by the stress potential under test. ---!
   pure function make_env(psi) result(e)
      real(c_double), intent(in) :: psi
      type(leaf_env_c) :: e
      e%par = 1500.0_c_double ; e%leaf_temp = 298.0_c_double ; e%vpd = 1500.0_c_double
      e%ca  = 400.0_c_double  ; e%pressure  = 101325.0_c_double
      e%psi_leaf = 0.0_c_double          ! demand-side limb off, so only `psi` moves gs
      e%gb  = 0.0_c_double
      e%psi = psi                         ! supply-side: what beta_stomata keys on
   end function make_env

   !----- The same C3 trait set the Python wrapper ships as its defaults (meds/plant/leaf.py), !
   !      so this test and the Python path exercise the same numbers. -------------------------!
   pure function c3_defaults() result(q)
      type(leaf_params_c) :: q
      q%pathway = int(PATH_C3, c_int)
      q%vcmax25 = 60.0_c_double  ; q%jmax25 = 108.0_c_double ; q%tpu25 = 1.0e6_c_double
      q%rd25    = 0.9_c_double   ; q%kp25   = 0.0_c_double
      q%g0      = 0.01_c_double  ; q%g1     = 4.0_c_double   ; q%d0    = 1500.0_c_double
      q%quantum_yield = 0.0_c_double
      q%theta_j = 0.85_c_double  ; q%theta_cj = 0.80_c_double ; q%theta_ic = 0.95_c_double
      q%lambda25 = 600.0_c_double
      q%psi_open = -0.5_c_double ; q%psi_close = -2.5_c_double ; q%lambda_psi_exp = 1.0_c_double
      q%sref_stomata = 2.0_c_double
      q%kc25 = 40.49_c_double ; q%ko25 = 27840.0_c_double ; q%gstar25 = 4.275_c_double
      q%ea_kc = 79430.0_c_double ; q%ea_ko = 36380.0_c_double ; q%ea_gstar = 37830.0_c_double
      q%ea_vcmax = 65330.0_c_double ; q%ea_jmax = 43540.0_c_double ; q%ea_rd = 46390.0_c_double
      q%hd_vcmax = 200000.0_c_double ; q%hd_jmax = 200000.0_c_double ; q%hd_rd = 1.0e9_c_double
      q%ds_vcmax = 650.0_c_double ; q%ds_jmax = 640.0_c_double ; q%ds_rd = 490.0_c_double
      q%o2_mol_frac = 0.209_c_double ; q%absorptance = 0.85_c_double ; q%phi_psii = 0.85_c_double
   end function c3_defaults

end program test_plant_capi
