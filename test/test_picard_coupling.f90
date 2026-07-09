!==========================================================================================!
! test_picard_coupling -- the P3 coupled-surface (Picard) fixed point in column_fast_step.      !
!   A. SPLIT UNCHANGED: SCHEME_SPLIT_SEQUENTIAL reproduces a golden CAS/soil anchor from the pre-  !
!      P3 operator-split sweep (all P3 changes are gated on the picard scheme; test_column_dynamics !
!      validates the full split trajectory + tight budgets). NOTE split != picard@1 once P3c re-    !
!      bases the leaf-LW emission temperature under picard only.                                    !
!   B. CONVERGENCE: PICARD @ max_iter=20 converges every sub-step (picard_nonconv==0) in a few     !
!      iterations (a contraction), and ALL conservation budgets still close (n_fail==0).           !
!   C. CONSISTENCY: the converged Picard trajectory is physical + bounded and MEASURABLY damps the !
!      explicit-split midday overshoot (the value of the fixed point), without blowing up.         !
!   D. FIXED-ITER: the GPU-uniform picard_fixed_iter path (no early exit) still closes every budget !
!      (guards the t_ground/t_bot handling on the non-early-exit path).                             !
!==========================================================================================!
program test_picard_coupling
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t, SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_thermo,              only : cas_enthalpy_of_temp, temp_to_uext
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_column_dynamics,     only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, column_fast_step
   use meds_plant_interface,     only : NODE_LEAF
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik, nstep = 96_ik
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: ccfg
   type(column_cohort_t)  :: coh
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: bio
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budg
   type(meds_time_t)      :: sim_date
   real(wp) :: tc_split(nstep), sh_split(nstep), lf_split(nstep), ss_split(nstep)
   real(wp) :: tc_p20(nstep),   sh_p20(nstep),   lf_p20(nstep),   ss_p20(nstep)
   real(wp) :: dmax_tc, dmax_lf
   integer(ik) :: nfail, istep

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp

   call alloc_column_cohort(coh, n)
   coh%pft(1) = 1_ik ; coh%lai(1) = 3.0_wp ; coh%wai(1) = 0.5_wp
   coh%height(1) = 16.0_wp ; coh%crown(1) = 0.9_wp
   coh%leaf_width(1) = 0.04_wp ; coh%branch_diam(1) = 0.02_wp
   coh%leaf_area(1) = 10.0_wp ; coh%nplant(1) = 0.3_wp ; coh%dbh(1) = 20.0_wp ; coh%broot(1) = 0.5_wp
   coh%bleaf(1) = 0.5_wp ; coh%bsap(1) = 5.0_wp ; coh%sap_area(1) = 0.01_wp

   call build_soil_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ccfg%soil)
   call build_soil_thermal(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ccfg%soil_thermal)
   ccfg%wood%is_woody = .true. ; ccfg%wood%stem_resp_factor25 = 0.06_wp ; ccfg%wood%agf_bs = 0.7_wp
   ccfg%root%root_resp_factor25 = 0.30_wp
   ccfg%co2%rh_k_base = 0.01_wp
   ccfg%fast_soil_carbon = 5.0_wp
   ccfg%hydro_p%leaf_pi0 = -1.5_wp ; ccfg%hydro_p%leaf_eps = 12.0_wp
   ccfg%hydro_p%leaf_af  = 0.30_wp ; ccfg%hydro_p%leaf_water_sat = 2.0_wp
   ccfg%hydro_p%wood_pi0 = -1.0_wp ; ccfg%hydro_p%wood_eps = 8.0_wp
   ccfg%hydro_p%wood_af  = 0.20_wp ; ccfg%hydro_p%wood_water_sat = 1.0_wp
   ccfg%hydro_p%wood_psi50 = -2.0_wp ; ccfg%hydro_p%wood_kexp = 2.0_wp
   ccfg%hydro_p%k_plant_max = 6.0e-4_wp ; ccfg%hydro_p%wood_kmax = 8.0_wp
   ccfg%hydro_p%vessel_curl = 1.5_wp
   ccfg%rhizo_cond = 5.0e-4_wp

   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n))

   print '(a)', '[test_picard_coupling]'

   !=== A. The SPLIT path is UNCHANGED by P3 (all P3 changes are gated on the picard scheme). ===!
   !       Guard it against a golden anchor from the pre-P3 operator-split sweep (test_column_       !
   !       dynamics validates the full split trajectory + tight budgets; this pins two exact points). !
   call run_day(SCHEME_SPLIT_SEQUENTIAL, 1_ik, tc_split, sh_split, lf_split, ss_split, budg)
   call ck(abs(tc_split(54) - 292.450065_wp) < 1.0e-3_wp .and.                                  &
           abs(ss_split(54) - 296.218258_wp) < 1.0e-3_wp,                                       &
           'split path unchanged (golden CAS + soil-surface temp at noon)',                     &
           abs(tc_split(54) - 292.450065_wp))

   !=== B. picard@20 -- CONVERGES + conserves. ===========================================!
   call run_day(SCHEME_PICARD_COUPLED, 20_ik, tc_p20, sh_p20, lf_p20, ss_p20, budg)
   call ck(budg%picard_nonconv == 0_ik, 'picard@20 converges every sub-step (nonconv==0)',      &
           real(budg%picard_nonconv, wp))
   call ck(budg%picard_iters >= 2_ik .and. budg%picard_iters <= 20_ik,                          &
           'picard actually iterated + converged within the cap (2 <= worst iters <= 20)',      &
           real(budg%picard_iters, wp))
   call ck(budg%cas_energy%n_fail == 0_ik .and. budg%cas_water%n_fail == 0_ik .and.             &
           budg%cas_co2%n_fail == 0_ik .and. budg%soil_energy%n_fail == 0_ik .and.              &
           budg%soil_water%n_fail == 0_ik .and. budg%whole_energy%n_fail == 0_ik .and.          &
           budg%whole_water%n_fail == 0_ik, 'picard@20: all conservation budgets close',        &
           real(budg%whole_energy%n_fail, wp))

   !=== C. picard is the IMPLICIT solution: physical + bounded, and MEASURABLY different from   !
   !       the explicit split sweep (it removes the operator-split midday overshoot -- the whole   !
   !       point of P3), without blowing up. ===================================================!
   dmax_tc = maxval(abs(tc_p20 - tc_split)) ; dmax_lf = maxval(abs(lf_p20 - lf_split))
   call ck(minval(tc_p20) > 280.0_wp .and. maxval(tc_p20) < 310.0_wp .and.                      &
           minval(lf_p20) > 280.0_wp .and. maxval(lf_p20) < 315.0_wp,                           &
           'picard trajectory is physical (bounded diurnal CAS + leaf temperature)', maxval(tc_p20))
   call ck(dmax_tc > 0.3_wp .and. dmax_tc < 12.0_wp,                                            &
           'picard differs from split (removes splitting-error overshoot) but stays bounded', dmax_tc)
   call ck(maxval(tc_split) > maxval(tc_p20),                                                   &
           'implicit picard damps the explicit-split midday overshoot', maxval(tc_split) - maxval(tc_p20))

   !=== D. picard_fixed_iter (GPU warp-uniform, NO early exit): must still close every budget. ==!
   !       Exercises the non-early-exit loop tail, where t_ground/t_bot must stay at the values the  !
   !       last pass's advection used (not the freshly diagnosed ones) or the energy budget leaks.   !
   ccfg%picard_fixed_iter = .true.
   call run_day(SCHEME_PICARD_COUPLED, 8_ik, tc_p20, sh_p20, lf_p20, ss_p20, budg)
   ccfg%picard_fixed_iter = .false.
   call ck(budg%cas_energy%n_fail == 0_ik .and. budg%cas_water%n_fail == 0_ik .and.              &
           budg%soil_energy%n_fail == 0_ik .and. budg%soil_water%n_fail == 0_ik .and.            &
           budg%whole_energy%n_fail == 0_ik .and. budg%whole_water%n_fail == 0_ik,               &
           'picard_fixed_iter (no early exit): all conservation budgets close',                  &
           real(budg%whole_energy%n_fail, wp))

   if (nfail == 0_ik) then
      print '(a)', 'test_picard_coupling: ALL PASSED'
      print '(a,i0,a,es10.3,a,es10.3,a)', '   (picard@20 worst iters=', budg%picard_iters,       &
            '  max|dTcas|=', dmax_tc, ' K  max|dTleaf|=', dmax_lf, ' K vs split)'
   else
      print '(a,i0,a)', 'test_picard_coupling: ', nfail, ' FAILED' ; error stop 1
   end if

contains

   subroutine ck(cond, name, val)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es12.4,a)', '  ok   : ', name, '  (', val, ')'
      else
         print '(a,a,a,es12.4)', '  FAIL : ', name, '  val = ', val ; nfail = nfail + 1_ik
      end if
   end subroutine ck

   !----- One 24 h diurnal integration under `scheme` with `max_iter`; captures per-step CAS   !
   !      temp / humidity, leaf temp, soil-surface temp; leaves the final budg in the module var.!
   subroutine run_day(scheme, max_iter, tc, sh, lf, ss, b)
      integer(ik),           intent(in)  :: scheme, max_iter
      real(wp),              intent(out) :: tc(:), sh(:), lf(:), ss(:)
      type(column_budget_t), intent(out) :: b
      real(wp)    :: t_sec, cosz, t_air
      integer(ik) :: is, k
      cfg%integration_scheme = scheme
      ccfg%picard_max_iter   = max_iter
      if (allocated(bio%leaf_temp)) deallocate(bio%leaf_temp)
      if (allocated(bio%psi))       deallocate(bio%psi)
      call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
      b = column_budget_t()
      bio%soil_w%theta(1:nsl) = theta0
      do k = 1_ik, nsl
         bio%soil_e%soil_energy(k) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(k),  &
                                     theta0 * rho_h2o, t0, 1.0_wp)
         bio%soil_e%soil_temp(k)   = t0
      end do
      do is = 1_ik, nstep
         t_sec = (real(is, wp) - 0.5_wp) * dt_fast
         cosz  = solar_cosz(sim_date, t_sec, lat)
         t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
         forc%abs_sw   = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
         forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
         forc%precip   = 0.0_wp
         if (is >= 12_ik .and. is <= 28_ik) forc%precip = 1.5e-4_wp
         forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
         forc%shv_atm      = 0.008_wp ; forc%co2_atm = 400.0_wp
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, b)
         tc(is) = bio%cas%can_temp ; sh(is) = bio%cas%can_shv
         lf(is) = bio%leaf_temp(1) ; ss(is) = bio%soil_e%soil_temp(1)
      end do
   end subroutine run_day

end program test_picard_coupling
