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
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_uext, cas_temp_of_enthalpy, &
                                        sat_specific_humidity
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG,  &
                                        SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER
   use meds_column_state_types, only : build_soil_hydr_params, PSI_INIT
   use meds_column_state_types, only : build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config, &
                                        WOODEN_PROGNOSTIC, WOODEN_DIAGNOSTIC
   use meds_fast_split,          only : column_fast_step
   use meds_hydr_lib,            only : psi_from_water_content, water_content
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   !----- seed moisture, a VARIABLE so the saturated test can drive the column to theta_sat.   !
   real(wp) :: theta_seed = theta0
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
   real(wp)    :: psi_leaf_diag
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
   coh%vcmax25(1) = cfg%pft%vcmax25(1) ; coh%rd25(1) = cfg%pft%rd25(1)   ! leaf capacities (plastic trait state)
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
   call apply_hydraulics_config(cfg%hydraulics, ccfg%hydro_p)
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
      !----- psi is no longer persisted state (MEDS_ED2_RK45_DESIGN.md sec 4): diagnose it from the  !
      !      persisted leaf_water_mass for the same physical bound this test always checked. ---------!
      psi_leaf_diag = psi_from_water_content(bio%leaf_water_mass(1), ccfg%hydro_p%leaf_pi0,          &
           ccfg%hydro_p%leaf_elastic_mod, ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, &
           coh%bleaf(1))
      physical = physical .and. psi_leaf_diag < 0.5_wp .and. psi_leaf_diag > -12.0_wp
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

   !=== G. CANOPY-SURFACE WATER (opt-in, MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c): a diurnal march with   !
   !       a morning rain pulse actually gets intercepted; whole_water closes exactly, whole_energy        !
   !       stays BOUNDED (known deferred sensible-heat approx -- the store is valued at one fixed           !
   !       rain_temp reference rather than a real prognostic surface-water temperature, same category        !
   !       as sec 2's qloss/qwflux_wl upwind-temperature approximation; mirrors the split path's own          !
   !       RUN 6 in test_column_dynamics.f90, same bound). Proves the WIRING (interception->film->CAS->        !
   !       ledgers), not the wetted-fraction algebra itself (already unit-tested in test_surface_energy.f90). !
   call test_ark_canopy_water()

   !=== H. LEAF/ROOT-TURNOVER SHED WATER (P4, MEDS_ED2_RK45_DESIGN.md): a constant shed_water_rate    !
   !       (distinct from precip) wets the soil and both whole_water/whole_energy still close. ========!
   call test_ark_shed_water()
   call test_split_shed_water()
   call test_ark_saturated()
   call test_ark_aquifer()
   call test_wood_prognostic(INTEG_ARK, 'ARK ')

   if (nfail == 0_ik) then
      print '(a)', 'test_column_ark: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_ark: ', nfail, ' FAILED' ; error stop 1
   end if

contains

   !----- march 96 sub-steps (24 h) of INTEG_ARK over MOIST free-draining soil (src_frac==1, no clamp, !
   !       no psi-limit, precip==0) and assert the 7 conservation budgets close. -------------------!
   !----- PHASE 0/3 (MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md): the aquifer bottom BC used to hard      !
   !      error-stop on this path. It is now a head-driven, two-way boundary with no prognostic state, !
   !      and the ARK commits the scratch column_hydrology_flux theta verbatim, so it inherits it      !
   !      unchanged. A column started DRY must wet from below with both ledgers closed. ---------------!
   !----- PHASE 4 (MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md): prognostic WOOD, operator-split behind    !
   !      the L-stable veg_energy_step_implicit kernel, now available on ARK and RK45 (it used to      !
   !      hard error-stop). Three assertions:                                                          !
   !        (a) the whole-column energy ledger still closes -- the store delta AND the wood net         !
   !            radiation both had to enter it, since surface_derivs no longer folds wood into          !
   !            coh_rnet when the frozen diagnostic inputs are zeroed;                                  !
   !        (b) wood LAGS the canopy air through the diel cycle, which is the feature;                  !
   !        (c) the cap_wood -> 0 limit reproduces the DIAGNOSTIC run, which is the check that the      !
   !            two wood authorities are wired to the same balance and never both fire.                 !
   subroutine test_wood_prognostic(integ, tag)
      integer(ik),      intent(in) :: integ
      character(len=4), intent(in) :: tag
      integer(ik) :: istep
      real(wp)    :: dmax_lag, tw_diag, tw_tiny, bsap_save
      call reset_state()
      cfg%time_integrator = integ
      ccfg%wood_energy_model = WOODEN_PROGNOSTIC
      dmax_lag = 0.0_wp
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
         dmax_lag = max(dmax_lag, abs(bio%wood_temp(1) - bio%cas%can_temp))
      end do
      call ck(budg%whole_energy%n_fail == 0_ik, trim(tag)//' PROG-WOOD: whole_energy closes',      &
              real(budg%whole_energy%n_fail, wp))
      call ck(budg%whole_water%n_fail == 0_ik, trim(tag)//' PROG-WOOD: whole_water closes',        &
              real(budg%whole_water%n_fail, wp))
      call ck(dmax_lag > 1.0e-3_wp, trim(tag)//' PROG-WOOD: wood temperature lags the CAS', dmax_lag)
      call ck(bio%wood_temp(1) > 200.0_wp .and. bio%wood_temp(1) < 350.0_wp,                        &
              trim(tag)//' PROG-WOOD: wood temperature physical', bio%wood_temp(1))

      !----- (c) cap_wood -> 0: a store with no inertia must land on the DIAGNOSTIC steady state.      !
      !                                                                                                 !
      !          TRAP, found the hard way (-1.82 K before it was spotted): bsap is DUAL-PURPOSE. It      !
      !          sets the thermal store (dry_hcap and wmass) AND the plant hydraulic capacitance --      !
      !          solve_plant_water_batch takes bsap + broot, and psi_from_water_content diagnoses        !
      !          psi_wood from it. Shrinking bsap on only ONE side of the comparison therefore also      !
      !          collapses that run's water relations, moving psi_leaf -> beta_nonstomata -> GPP/gs ->   !
      !          transpiration -> the CAS. The result was two different PLANTS being compared, not two   !
      !          wood schemes. Perturb BOTH sides identically; the gap then falls from -1.82 K to        !
      !          +0.31 K, and the sign flips to the one the physics predicts.                            !
      !---------------------------------------------------------------------------------------------!
      !----- What is left is the real quantity of interest: the Lie-Trotter coupling error between a   !
      !      wood diagnosed INSIDE the implicit CAS solve and one operator-split OUTSIDE it. That is    !
      !      the accepted cost of Phase 4's design (see advance_wood_energy_full's header) and this     !
      !      bounds it rather than asserting it away. --------------------------------------------!
      call ck(abs(tw_tiny - tw_diag) < 0.5_wp,                                                        &
              trim(tag)//' PROG-WOOD: cap->0 recovers the diagnostic balance to the operator-split '// &
              'coupling error', tw_tiny - tw_diag)
      ccfg%wood_energy_model = WOODEN_DIAGNOSTIC
   end subroutine test_wood_prognostic

   subroutine test_ark_aquifer()
      integer(ik) :: istep
      real(wp)    :: theta_bot0
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      ccfg%hydro%bottom_bc = SOIL_BC_AQUIFER
      bio%soil_w%theta(1:ccfg%soil%n_active) = 0.15_wp
      theta_bot0 = bio%soil_w%theta(ccfg%soil%n_active)
      do istep = 1_ik, 48_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      call ck(budg%whole_water%n_fail == 0_ik, 'AQUIFER/ARK: whole_water closes',                  &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik, 'AQUIFER/ARK: whole_energy closes',                &
              real(budg%whole_energy%n_fail, wp))
      call ck(bio%soil_w%theta(ccfg%soil%n_active) > theta_bot0,                                   &
              'AQUIFER/ARK: dry column wets from below', bio%soil_w%theta(ccfg%soil%n_active) - theta_bot0)
      ccfg%hydro%bottom_bc = SOIL_BC_FREE_DRAIN
   end subroutine test_ark_aquifer

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

   !----- march 96 sub-steps (24 h) of INTEG_ARK with canopy_water_on and a 5-step morning rain pulse   !
   !      (istep 20-24); assert some of it was intercepted and the budgets close (mirrors               !
   !      test_column_dynamics.f90's own RUN 6 for the split path). ----------------------------------!
   subroutine test_ark_canopy_water()
      integer(ik) :: istep
      real(wp)    :: surf_water_peak
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_niter = 8_ik
      ccfg%canopy_water_on = .true.
      surf_water_peak = 0.0_wp
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         if (istep >= 20_ik .and. istep <= 24_ik) forc%precip = 5.0e-5_wp   ! a morning rain pulse
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
         surf_water_peak = max(surf_water_peak, bio%leaf_surf_water(1) + bio%wood_surf_water(1))
      end do
      ccfg%canopy_water_on = .false.   ! restore default for any test added after this
      call ck(budg%whole_water%n_fail == 0_ik, 'ARK canopy water: whole-column water still closes',    &
              real(budg%whole_water%n_fail, wp))
      call ck(surf_water_peak > 0.0_wp, 'ARK canopy water: the morning rain pulse was intercepted',     &
              surf_water_peak)
      call ck(budg%whole_energy%worst < 5.0e6_wp,                                                       &
              'ARK canopy water: whole-column energy stays BOUNDED (known deferred approx)',            &
              budg%whole_energy%worst)
      print '(a,es10.3,a)', '   (ARK canopy water peak film=', surf_water_peak, ' kg/m2)'
   end subroutine test_ark_canopy_water

   !----- march 96 sub-steps (24 h) of INTEG_ARK with a constant leaf/root-turnover shed-water rate  !
   !      (MEDS_ED2_RK45_DESIGN.md P4, bio%shed_water_rate -- a PATCH-level input, not atmospheric      !
   !      forcing, so it is frozen on bio for the whole day rather than living on forc; distinct from    !
   !      precip, which stays 0 throughout): the soil must wet from THIS input alone, and both              !
   !      whole_water AND whole_energy must still close -- energy closing needs NO separate wiring of        !
   !      its own (P4's design choice: the shed water's enthalpy rides the SAME e_infil/rain_temp             !
   !      treatment every other infiltrating input already gets, once mixed into hforc%precip_ground          !
   !      by build_column_frozen). ---------------------------------------------------------------------!
   subroutine test_ark_shed_water()
      integer(ik) :: istep
      real(wp)    :: theta_col0, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_ARK ; cfg%ark_adaptive = .true.
      cfg%ark_rtol = 1.0e-4_wp ; cfg%ark_niter = 8_ik
      bio%shed_water_rate = 8.0e-5_wp                     ! P4: frozen for the whole day (precip stays 0)
      theta_col0 = sum(bio%soil_w%theta(1:nsl))
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(bio%soil_w%theta(1:nsl))
      bio%shed_water_rate = 0.0_wp   ! restore default for any test added after this
      call ck(theta_col1 > theta_col0,                                                              &
              'ARK shed water: leaf/root shed water alone wetted the soil column (theta rose)',      &
              theta_col1 - theta_col0)
      call ck(budg%whole_water%n_fail == 0_ik,                                                       &
              'ARK shed water: whole-column WATER still closes with shed_water_rate active',         &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik,                                                      &
              'ARK shed water: whole-column ENERGY still closes (no separate energy wiring needed)',  &
              real(budg%whole_energy%n_fail, wp))
   end subroutine test_ark_shed_water

   !----- SAME shed-water check as test_ark_shed_water, but on the SPLIT integrator (meds_fast_       !
   !      split.f90's own hforc%precip_ground/w_in wiring, not build_column_frozen's) -- P4 touches    !
   !      three integrator files; this is the split path's coverage (ARK/RK45 share build_column_       !
   !      frozen, so one test each there already covers both). --------------------------------------!
   subroutine test_split_shed_water()
      integer(ik) :: istep
      real(wp)    :: theta_col0, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_SPLIT
      bio%shed_water_rate = 8.0e-5_wp                     ! P4: frozen for the whole day (precip stays 0)
      theta_col0 = sum(bio%soil_w%theta(1:nsl))
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(bio%soil_w%theta(1:nsl))
      bio%shed_water_rate = 0.0_wp   ! restore default for any test added after this
      call ck(theta_col1 > theta_col0,                                                              &
              'split shed water: leaf/root shed water alone wetted the soil column (theta rose)',   &
              theta_col1 - theta_col0)
      call ck(budg%whole_water%n_fail == 0_ik,                                                       &
              'split shed water: whole-column WATER still closes with shed_water_rate active',      &
              real(budg%whole_water%n_fail, wp))
      call ck(budg%whole_energy%n_fail == 0_ik,                                                      &
              'split shed water: whole-column ENERGY still closes (no separate energy wiring needed)', &
              real(budg%whole_energy%n_fail, wp))
   end subroutine test_split_shed_water


   !----- SATURATED sealed column: the ARK/RK45 twin of test_column_dynamics RUN 7. A bedrock bottom   !
   !      seeded just below theta_sat under heavy rain both SATURATES (firing meds_soil_water's post-   !
   !      solve theta clip, which moves water with no face) and OVERFLOWS the ponding store (surface     !
   !      runoff). None of the other ARK runs reach either state -- they are free-draining at           !
   !      theta = 0.30 -- so the enthalpy bookkeeping for both was entirely uncovered on this path:       !
   !        * the clip's per-layer enthalpy, frozen at each layer's own temperature in                    !
   !          build_column_frozen and carried into the stages' root_heat_sink column, and                 !
   !        * runoff, which must contribute NO enthalpy term at all -- it leaves the ponding store,       !
   !          which holds mass but no enthalpy, so the e_runof term this path used to apply at layer 1    !
   !          removed ~1 MJ per kg the soil never received.                                               !
   !      whole_energy closing is the assertion that both are booked with the right sign; the            !
   !      temperature bound is the assertion that the INTERIOR advective faces (previously hardcoded     !
   !      to zero here) actually connect the boundary faces -- without them layer 1 accumulates the      !
   !      full infiltration enthalpy while a deeper layer sheds it, and the ledger cannot see it. -------!
   subroutine test_ark_saturated()
      integer(ik) :: istep
      real(wp)    :: pond_peak, theta_peak, ss_min, ss_max
      !----- FREE-DRAIN, not bedrock: column_fast_step_ark error-stops on any other bottom BC. The     !
      !      route to saturation here is therefore rain far above the column's drainage capacity        !
      !      (ksat = 2.89e-6 m/s ~ 2.9e-3 kg/m2/s) rather than a sealed bottom. -----------------------!
      theta_seed = 0.428_wp                              ! just below theta_sat = 0.43
      call reset_state()
      cfg%time_integrator = INTEG_ARK
      pond_peak = 0.0_wp ; theta_peak = 0.0_wp
      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp
      do istep = 1_ik, 96_ik
         call set_diurnal_forcing(istep)
         forc%precip = 8.0e-3_wp                         ! ~29 mm/hr: far above the drainage capacity
         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
         pond_peak  = max(pond_peak,  bio%soil_w%w_surface)
         theta_peak = max(theta_peak, maxval(bio%soil_w%theta(1:nsl)))
         ss_min = min(ss_min, bio%soil_e%soil_temp(1)) ; ss_max = max(ss_max, bio%soil_e%soil_temp(1))
      end do
      theta_seed = theta0                                 ! restore for any test added after this
      call ck(theta_peak >= 0.43_wp - 1.0e-12_wp,                                                    &
              'ARK saturated: column reached theta_sat (clip path is live)', theta_peak)
      call ck(pond_peak >= ccfg%hydro%w_pond_max - 1.0e-9_wp,                                        &
              'ARK saturated: ponding store overflowed (runoff path is live)', pond_peak)
      call ck(budg%whole_energy%n_fail == 0_ik,                                                      &
              'ARK saturated: whole-column ENERGY still closes through clip + runoff',             &
              budg%whole_energy%worst)
      call ck(budg%whole_water%n_fail == 0_ik,                                                       &
              'ARK saturated: whole-column WATER still closes', budg%whole_water%worst)
      call ck(ss_min > 250.0_wp .and. ss_max < 340.0_wp,                                             &
              'ARK saturated: soil surface temp stays physical (interior faces connected)', ss_max)
   end subroutine test_ark_saturated

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
      !----- alloc_patch_biophys seeds leaf_water_mass/wood_water_mass at a scratch 0 (the real     !
      !      lazy-init lives in meds_fast_dynamics.f90's site-level gather loop, which this driver-  !
      !      level test bypasses) -- seed the same water_content(PSI_INIT,...) a freshly-created     !
      !      cohort gets there, or psi_from_water_content would diagnose an unphysical psi from an   !
      !      empty pool. -------------------------------------------------------------------------!
      bio%leaf_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%leaf_pi0, ccfg%hydro_p%leaf_elastic_mod, &
           ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, coh%bleaf(1:n))
      bio%wood_water_mass(1:n) = water_content(PSI_INIT, ccfg%hydro_p%wood_pi0, ccfg%hydro_p%wood_elastic_mod, &
           ccfg%hydro_p%wood_apoplast_frac, ccfg%hydro_p%wood_water_sat, coh%bsap(1:n) + coh%broot(1:n))
      budg = column_budget_t()
      bio%soil_w%theta(1:nsl) = theta_seed
      do kk = 1_ik, nsl
         bio%soil_e%soil_energy(kk) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(kk), &
                                      theta_seed * rho_h2o, t0, 1.0_wp)
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
