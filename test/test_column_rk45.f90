!==========================================================================================!
! test_column_rk45 -- the INTEG_RK4 wiring of the ED2-faithful adaptive Cash-Karp integrator    !
! into column_fast_step (MEDS_ED2_RK45_DESIGN.md P2).                                          !
!   A. GPP PARITY: build_column_frozen (shared with ARK) runs the SAME frozen pre-pass as the     !
!      split, so gpp_coh from the RK45 path is bit-identical to the split's on step 1.            !
!   B. PHYSICAL: a 24-step dry-window march under INTEG_RK4 keeps the CAS / soil / leaf / plant    !
!      water-mass state finite, bounded, and sub-saturated.                                        !
!   C. WHOLE-COLUMN CONSERVATION LEDGER: the water AND energy budgets close over a 24 h dry          !
!      diurnal march, machine-precision (the design doc's headline gates 2/3) -- unlike ARK,         !
!      RK45 has no operator split at all (mass and soil water are genuinely integrated), so there   !
!      is no tolerance inflation to remove; it is gated at split's own tight tolerance directly.     !
!==========================================================================================!
program test_column_rk45
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t, INTEG_SPLIT, INTEG_RK4
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_uext, cas_temp_of_enthalpy, &
                                        sat_specific_humidity
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_column_state_types, only : build_soil_hydr_params, PSI_INIT
   use meds_column_state_types, only : build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config
   use meds_fast_split,          only : column_fast_step
   use meds_hydr_lib,            only : psi_from_water_content, water_content
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
   real(wp)    :: gpp_split(n), gpp_rk45(n), gpp_coh(n), tcas, qsat, worst_super
   real(wp)    :: psi_leaf_diag
   integer(ik) :: nfail, is, k
   logical     :: physical

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- column setup (mirrors test_column_ark). ---------------------------------------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   call alloc_column_cohort(coh, n)
   coh%pft(1) = 1_ik ; coh%lai(1) = 3.0_wp ; coh%wai(1) = 0.5_wp
   coh%vcmax25(1) = cfg%pft%vcmax25(1) ; coh%rd25(1) = cfg%pft%rd25(1)
   coh%height(1) = 16.0_wp ; coh%crown(1) = 0.9_wp
   coh%leaf_width(1) = 0.04_wp ; coh%branch_diam(1) = 0.02_wp
   coh%leaf_area(1) = 10.0_wp ; coh%nplant(1) = 0.3_wp ; coh%dbh(1) = 20.0_wp ; coh%broot(1) = 0.5_wp
   coh%bleaf(1) = 0.5_wp ; coh%bsap(1) = 5.0_wp ; coh%sap_area(1) = 0.01_wp
   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ccfg%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ccfg%soil_thermal)
   ccfg%wood%is_woody = .true. ; ccfg%wood%stem_resp_factor25 = 0.06_wp ; ccfg%wood%agf_bs = 0.7_wp
   ccfg%root%root_resp_factor25 = 0.30_wp
   ccfg%co2%rh_k_base = 0.01_wp
   ccfg%fast_soil_carbon = 5.0_wp
   call apply_hydraulics_config(cfg%hydraulics, ccfg%hydro_p, ccfg%rhizo_cond)
   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   print '(a)', '[test_column_rk45]'

   !=== A. GPP parity on step 1 (identical initial state; RK45 pre-pass == split pre-pass). ====!
   call set_noon_forcing()
   call reset_state()
   cfg%time_integrator = INTEG_SPLIT
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   gpp_split = gpp_coh
   call reset_state()
   cfg%time_integrator = INTEG_RK4
   call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
   gpp_rk45 = gpp_coh
   call ck(abs(gpp_rk45(1) - gpp_split(1)) < 1.0e-12_wp,                                          &
           'RK45 pre-pass gpp bit-identical to the split (build_column_frozen)', abs(gpp_rk45(1) - gpp_split(1)))
   call ck(gpp_rk45(1) > 0.0_wp, 'RK45 midday gpp > 0', gpp_rk45(1))

   !=== B. A dry-window march under INTEG_RK4 stays physical + bounded + sub-saturated. ========!
   call reset_state()
   cfg%time_integrator = INTEG_RK4
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
      physical = physical .and. bio%leaf_water_mass(1) > 0.0_wp .and. bio%wood_water_mass(1) > 0.0_wp
      psi_leaf_diag = psi_from_water_content(bio%leaf_water_mass(1), ccfg%hydro_p%leaf_pi0,       &
           ccfg%hydro_p%leaf_elastic_mod, ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, &
           coh%bleaf(1))
      physical = physical .and. psi_leaf_diag < 0.5_wp .and. psi_leaf_diag > -12.0_wp
   end do
   call ck(physical, 'INTEG_RK4 dry-window march stays physical + bounded (24 steps)', bio%cas%can_temp)
   call ck(worst_super <= 1.0e-4_wp, 'INTEG_RK4 CAS stays sub-saturated', worst_super)

   !=== C. WHOLE-COLUMN CONSERVATION LEDGER: the water AND energy budgets close over a 24 h dry   !
   !       diurnal march -- RK45 has no operator split (mass + soil water are genuinely             !
   !       integrated), so this is gated at split's own tight tolerance (1e-4 water, 1 J energy),    !
   !       not ARK's inflated one. =================================================================!
   call test_rk45_budgets()

   !=== D. PRECIP>0: a WET diurnal march runs and still closes the whole-column water AND energy   !
   !       budgets to the SAME tight (split/gate 2-3) tolerance -- unlike ARK's wet test (which        !
   !       needs the lagged-ponding operator-split water tolerance), RK45 genuinely integrates soil     !
   !       water, so there is no split-vs-continuous mismatch to tolerate; this is the test that        !
   !       exercises the P2-part-3 boundary water-enthalpy advection fix (column_derivs's e_infil/       !
   !       e_runof/e_drain -> root_heat_sink) under a REAL, substantial precip rate, not just the         !
   !       incidental background drainage that first caught it. ===================================!
   call test_rk45_budgets_wet()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_rk45: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_rk45: ', nfail, ' FAILED' ; error stop 1
   end if

contains

   !----- march 96 sub-steps (24 h) of INTEG_RK4 over MOIST free-draining soil and assert the      !
   !      whole-column water + energy budgets close (design doc sec 8 gates 2/3). -----------------!
   subroutine test_rk45_budgets()
      integer(ik) :: istep
      call reset_state()
      cfg%time_integrator = INTEG_RK4
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      call ck(budg%whole_water%n_fail == 0_ik, 'RK45: whole_water closes (n_fail==0)',            &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik, 'RK45: whole_energy closes (n_fail==0)',          &
              real(budg%whole_energy%n_fail, wp))
      call ck(budg%whole_water%n_check == 96_ik, 'RK45: ledger fired every dispatched dt_fast',   &
              real(budg%whole_water%n_check, wp))
      call ck(budg%whole_energy%worst < 1.0_wp, 'RK45: whole-energy closes < 1 J',                &
              budg%whole_energy%worst)
      call ck(budg%integ_nsteps >= 1_ik .and. budg%integ_nsteps <= 64_ik,                          &
              'RK45: adaptive substeps bounded (sec 6: ~3 expected, cap 64)', real(budg%integ_nsteps, wp))
      print '(a,i0,a,i0)', '   RK45 last dt_fast: substeps = ', budg%integ_nsteps,                 &
            ' , rejects = ', budg%integ_nrej
      print '(a,es10.3,a,es10.3,a)', '   (worst whole-column resid: energy= ', budg%whole_energy%worst, &
            ' J/m2  water= ', budg%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets

   !----- march 96 sub-steps (24 h) of INTEG_RK4 over free-draining soil WITH continuous rain       !
   !      (precip>0), mirroring test_column_ark's test_ark_budgets_wet. Asserts the run completes,   !
   !      the soil wets, and BOTH whole-column budgets close at the tight (non-inflated) tolerance.   !
   subroutine test_rk45_budgets_wet()
      integer(ik) :: istep
      real(wp)    :: theta_col0, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_RK4
      theta_col0 = sum(bio%soil_w%theta(1:nsl))
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         forc%precip = 8.0e-5_wp                         ! ~0.29 mm/hr continuous rain (precip>0)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(bio%soil_w%theta(1:nsl))
      call ck(budg%whole_water%n_check == 96_ik, 'RK45 wet: ran 96 wet steps (no guard error stop)',   &
              real(budg%whole_water%n_check, wp))
      call ck(theta_col1 > theta_col0, 'RK45 wet: rain wetted the soil column (theta rose)',           &
              theta_col1 - theta_col0)
      call ck(budg%whole_energy%n_fail == 0_ik, 'RK45 wet: whole_energy closes (n_fail==0)',           &
              real(budg%whole_energy%n_fail, wp))
      call ck(budg%whole_water%n_fail == 0_ik, 'RK45 wet: whole_water closes (n_fail==0)',             &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%worst < 1.0_wp, 'RK45 wet: whole-energy closes < 1 J',                 &
              budg%whole_energy%worst)
      print '(a,es10.3,a,es10.3,a)', '   (RK45 wet worst whole-column resid: energy= ',                &
            budg%whole_energy%worst, ' J/m2  water= ', budg%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets_wet

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
      call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
      bio%leaf_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%leaf_pi0, ccfg%hydro_p%leaf_elastic_mod, &
           ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, coh%bleaf(1:n))
      bio%wood_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%wood_pi0, ccfg%hydro_p%wood_elastic_mod, &
           ccfg%hydro_p%wood_apoplast_frac, ccfg%hydro_p%wood_water_sat, coh%bsap(1:n) + coh%broot(1:n))
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
      forc%precip = 0.0_wp
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
      forc%precip = 0.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
   end subroutine set_diurnal_forcing

end program test_column_rk45
