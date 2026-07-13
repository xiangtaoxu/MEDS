!==========================================================================================!
! test_column_ark -- the INTEG_ARK wiring of the coupled IMEX-ARK integrator into column_fast_step.  !
!   A. GPP PARITY: build_column_frozen runs the SAME frozen pre-pass as the split, so gpp_coh from     !
!      the ARK path is bit-identical to the split's on step 1 (identical initial state).               !
!   B. PHYSICAL: a 24-step dry-window march under INTEG_ARK keeps the CAS / soil / leaf / psi state    !
!      finite, bounded, and sub-saturated (the integrator + unpack are correct end to end).            !
!   C. The split path is validated unchanged by test_picard_coupling (the golden anchor); here we only !
!      exercise the opt-in ARK path (free-drain, no Zeng-Decker).                                      !
!   F. PRECIP>0 guard-lift: a wet march runs (no guard error stop), wets the soil, and closes the      !
!      budgets (ENERGY machine via the frozen boundary advection, WATER to the lagged-ponding split).  !
!==========================================================================================!
program test_column_ark
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t, INTEG_SPLIT, INTEG_ARK
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_thermo,              only : cas_enthalpy_of_temp, temp_to_uext, cas_temp_of_enthalpy, &
                                        sat_specific_humidity
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_column_dynamics,     only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, column_fast_step
   use meds_plant_interface,     only : NODE_LEAF, NODE_WOOD
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik
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
   real(wp)    :: gpp_split(n), gpp_ark(n), gpp_coh(n), tcas, qsat, worst_super, tcas_1, tcas_8
   integer(ik) :: nfail, is, k
   logical     :: physical

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- column setup (mirrors test_picard_coupling). ---------------------------------------!
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
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   print '(a)', '[test_column_ark]'

   !=== A. GPP parity on step 1 (identical initial state; ARK pre-pass == split pre-pass). =====!
   call set_noon_forcing()
   call reset_state()
   cfg%time_integrator = INTEG_SPLIT
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   gpp_split = gpp_coh
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   gpp_ark = gpp_coh
   call ck(abs(gpp_ark(1) - gpp_split(1)) < 1.0e-12_wp,                                          &
           'ARK pre-pass gpp bit-identical to the split (build_column_frozen)', abs(gpp_ark(1) - gpp_split(1)))
   call ck(gpp_ark(1) > 0.0_wp, 'ARK midday gpp > 0', gpp_ark(1))

   !=== B. A dry-window march under INTEG_ARK stays physical + bounded + sub-saturated. ========!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
   physical = .true. ; worst_super = -1.0_wp
   do is = 1_ik, 24_ik
      call set_diurnal_forcing(is)
      call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      tcas = bio%cas%can_temp
      qsat = sat_specific_humidity(tcas, aenv%press)
      worst_super = max(worst_super, bio%cas%can_shv - qsat)
      physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. bio%cas%can_shv > 0.0_wp
      physical = physical .and. bio%leaf_temp(1) > 250.0_wp .and. bio%leaf_temp(1) < 340.0_wp
      do k = 1_ik, nsl
         physical = physical .and. bio%soil_w%theta(k) > 0.0_wp .and. bio%soil_w%theta(k) < 0.6_wp
         physical = physical .and. bio%soil_e%soil_temp(k) > 260.0_wp .and. bio%soil_e%soil_temp(k) < 340.0_wp
      end do
      physical = physical .and. bio%psi(NODE_LEAF,1) < 0.5_wp .and. bio%psi(NODE_LEAF,1) > -12.0_wp
   end do
   call ck(physical, 'INTEG_ARK dry-window march stays physical + bounded (24 steps)', bio%cas%can_temp)
   call ck(worst_super <= 1.0e-4_wp, 'INTEG_ARK CAS stays sub-saturated', worst_super)

   !=== C. Fixed-substep (GPU-lockstep) path also runs + stays physical. =======================!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .false. ; cfg%ark_fixed_substep = 4_ik
   call set_noon_forcing()
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   call ck(bio%cas%can_temp > 270.0_wp .and. bio%cas%can_temp < 325.0_wp,                        &
           'INTEG_ARK fixed-substep path physical', bio%cas%can_temp)

   !=== D. ark_niter reaches the inner solver (np<=1 baseline vs np>1 Newton must differ). ======!
   call reset_state()
   cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true. ; cfg%ark_rtol = 1.0e-4_wp
   call set_noon_forcing()
   cfg%ark_niter = 1_ik
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   tcas_1 = bio%cas%can_temp
   call reset_state() ; call set_noon_forcing()
   cfg%ark_niter = 8_ik
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   tcas_8 = bio%cas%can_temp
   call ck(abs(tcas_1 - tcas_8) > 1.0e-4_wp, 'ARK ark_niter reaches ark2 (baseline vs Newton differ)', abs(tcas_1-tcas_8))
   call ck(tcas_1 > 270.0_wp .and. tcas_8 < 325.0_wp, 'both niter paths physical', tcas_8)

   !=== E. WHOLE-COLUMN CONSERVATION LEDGER: the 7 budgets close over a 24 h dry diurnal march. ==!
   call test_ark_budgets(.true.)      ! adaptive path (accept/reject accumulation)
   call test_ark_budgets(.false.)     ! fixed-substep path (uniform accumulation)

   !=== F. PRECIP>0 GUARD-LIFT: a WET diurnal march runs (no guard error stop), wets the soil, and     !
   !       still closes the budgets -- ENERGY to machine precision (the rain/runoff/drainage advection  !
   !       is a frozen fixed source), WATER to the lagged-ponding operator-split tolerance. ============!
   call test_ark_budgets_wet()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_ark: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_ark: ', nfail, ' FAILED' ; error stop 1
   end if

contains

   !----- march 96 sub-steps (24 h) of INTEG_ARK over MOIST free-draining soil (src_frac==1, no clamp, !
   !       no psi-limit, precip==0) and assert the 7 conservation budgets close. -------------------!
   subroutine test_ark_budgets(adaptive)
      logical, intent(in) :: adaptive
      integer(ik) :: istep
      character(len=6) :: tag
      tag = merge('adapt ', 'fixed ', adaptive)
      call reset_state()                               ! theta0=0.30 (moist), budg zeroed
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = adaptive
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_fixed_substep = 4_ik ; cfg%ark_niter = 8_ik
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)               ! precip==0 always (dry); diurnal SW
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      call ck(budg%cas_energy%n_fail == 0_ik .and. budg%cas_water%n_fail == 0_ik .and.            &
              budg%cas_co2%n_fail == 0_ik .and. budg%soil_energy%n_fail == 0_ik .and.             &
              budg%soil_water%n_fail == 0_ik .and. budg%whole_energy%n_fail == 0_ik .and.         &
              budg%whole_water%n_fail == 0_ik, 'ARK '//trim(tag)//': all 7 budgets close (n_fail==0)', &
              real(budg%whole_energy%n_fail, wp))
      call ck(budg%cas_energy%n_check == 96_ik, 'ARK '//trim(tag)//': ledger fired every dispatched dt_fast', &
              real(budg%cas_energy%n_check, wp))
      call ck(budg%cas_energy%worst < 1.0e-3_wp, 'ARK '//trim(tag)//': CAS energy machine-precision closure', &
              budg%cas_energy%worst)
      call ck(budg%whole_energy%worst < 1.0_wp, 'ARK '//trim(tag)//': whole-energy closes < 1 J', &
              budg%whole_energy%worst)
   end subroutine test_ark_budgets

   !----- march 96 sub-steps (24 h) of INTEG_ARK over free-draining soil WITH continuous rain (precip>0 !
   !       -> the guard-lift path: boundary water-enthalpy advection + persisted ponding). Assert the    !
   !       run completes, the soil wets, and the budgets close (ENERGY machine, WATER split-tol). ------!
   subroutine test_ark_budgets_wet()
      integer(ik) :: istep
      real(wp)    :: theta_col0, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_niter = 8_ik
      theta_col0 = sum(bio%soil_w%theta(1:nsl))
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         forc%precip = 8.0e-5_wp                         ! ~0.29 mm/hr continuous rain (precip>0)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(bio%soil_w%theta(1:nsl))
      call ck(budg%cas_energy%n_check == 96_ik, 'ARK wet: ran 96 wet steps (no guard error stop)',      &
              real(budg%cas_energy%n_check, wp))
      call ck(theta_col1 > theta_col0, 'ARK wet: rain wetted the soil column (theta rose)',             &
              theta_col1 - theta_col0)
      call ck(budg%cas_energy%n_fail == 0_ik .and. budg%soil_energy%n_fail == 0_ik .and.                &
              budg%whole_energy%n_fail == 0_ik, 'ARK wet: ENERGY budgets close (incl. advection)',      &
              budg%whole_energy%worst)
      call ck(budg%soil_water%n_fail == 0_ik .and. budg%whole_water%n_fail == 0_ik,                     &
              'ARK wet: WATER budgets close (lagged-ponding split tolerance)', budg%whole_water%worst)
      call ck(budg%soil_energy%worst < 1.0e-3_wp, 'ARK wet: soil energy machine-precision closure',     &
              budg%soil_energy%worst)
   end subroutine test_ark_budgets_wet

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

   subroutine reset_state()
      integer(ik) :: kk
      if (allocated(bio%leaf_temp)) deallocate(bio%leaf_temp)
      if (allocated(bio%psi))       deallocate(bio%psi)
      call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
      budg = column_budget_t()
      bio%soil_w%theta(1:nsl) = theta0
      do kk = 1_ik, nsl
         bio%soil_e%soil_energy(kk) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(kk), &
                                      theta0 * rho_h2o, t0, 1.0_wp)
         bio%soil_e%soil_temp(kk)   = t0
      end do
   end subroutine reset_state

   subroutine set_noon_forcing()
      forc%abs_sw = 450.0_wp ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 70.0_wp ; forc%abs_lw_ground = 0.0_wp
      forc%precip = 0.0_wp                                    ! inert-hydrology regime (ARK MVP)
      forc%enthalpy_atm = cas_enthalpy_of_temp(295.0_wp, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
   end subroutine set_noon_forcing

   subroutine set_diurnal_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
      forc%abs_sw = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%precip = 0.0_wp                                    ! dry window (no precip -> inert hydrology)
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
   end subroutine set_diurnal_forcing

end program test_column_ark
