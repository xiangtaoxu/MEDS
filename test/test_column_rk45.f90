!==========================================================================================!
! test_column_rk45 -- the INTEG_RK45 wiring of the ED2-faithful adaptive Cash-Karp integrator    !
! into column_fast_step (MEDS_ED2_RK45_DESIGN.md P2).                                          !
!   A. GPP PARITY: build_column_frozen (shared with ARK) runs the SAME frozen pre-pass as the     !
!      split, so gpp_coh from the RK45 path is bit-identical to the split's on step 1.            !
!   B. PHYSICAL: a 24-step dry-window march under INTEG_RK45 keeps the CAS / soil / leaf / plant    !
!      water-mass state finite, bounded, and sub-saturated.                                        !
!   C. WHOLE-COLUMN CONSERVATION LEDGER: the water AND energy budgets close over a 24 h dry          !
!      diurnal march, machine-precision (the design doc's headline gates 2/3) -- unlike ARK,         !
!      RK45 has no operator split at all (mass and soil water are genuinely integrated), so there   !
!      is no tolerance inflation to remove; it is gated at split's own tight tolerance directly.     !
!==========================================================================================!
program test_column_rk45
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t, INTEG_ARK, INTEG_RK45
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_therm_lib,              only : cas_enthalpy_of_temp, temp_to_internal_energy, cas_temp_of_enthalpy, &
                                        sat_specific_humidity
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out
   use meds_fast_types, only : patch_biophys_t, alloc_patch_biophys
   use meds_hydr_lib, only : SOIL_RETENTION_VG
   use meds_biophysics_opts, only : SOIL_BC_BEDROCK, SOIL_BC_FREE_DRAIN, SOIL_BC_AQUIFER
   use meds_canopy_types, only : set_aero_env_atm
   use meds_column_constants, only : PSI_INIT
   use meds_column_params, only : build_soil_hydr_params
   use meds_column_params, only : build_soil_therm_params
   use meds_fast_types,          only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, apply_hydraulics_config
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_fast_step,          only : column_fast_step
   use meds_hydr_lib,            only : psi_from_water_content, water_content, soil_psi_from_theta
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik
   real(wp),    parameter :: dt_fast = 150.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   !----- seed moisture, a VARIABLE so the saturated test can drive the column to theta_sat.   !
   real(wp) :: theta_seed = theta0
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: col_config
   type(column_cohort_t)  :: col_cohort
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: biophys
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budget
   type(meds_time_t)      :: sim_date
   real(wp)    :: gpp_split(n), gpp_rk45(n), gpp_coh(n), tcas, qsat, worst_super
   real(wp)    :: psi_leaf_diag, psi_leaf_probe(n)
   integer(ik) :: nfail, is, k
   logical     :: physical

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- column setup (mirrors test_column_ark). ---------------------------------------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   call alloc_column_cohort(col_cohort, n)
   col_cohort%pft(1) = 1_ik ; col_cohort%lai(1) = 3.0_wp ; col_cohort%wai(1) = 0.5_wp
   col_cohort%vcmax25(1) = cfg%pft%vcmax25(1) ; col_cohort%rd25(1) = cfg%pft%rd25(1)
   col_cohort%height(1) = 16.0_wp ; col_cohort%crown(1) = 0.9_wp
   col_cohort%leaf_width(1) = 0.04_wp ; col_cohort%branch_diam(1) = 0.02_wp
   col_cohort%leaf_area(1) = 10.0_wp ; col_cohort%nplant(1) = 0.3_wp ; col_cohort%dbh(1) = 20.0_wp ; col_cohort%broot(1) = 0.5_wp
   col_cohort%bleaf(1) = 0.5_wp ; col_cohort%bsap(1) = 5.0_wp ; col_cohort%sap_area(1) = 0.01_wp
   call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, col_config%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, col_config%soil_thermal)
   col_config%wood%is_woody = .true. ; col_config%wood%stem_resp_factor25 = 0.06_wp ; col_config%wood%agf_bs = 0.7_wp
   col_config%root%root_resp_factor25 = 0.30_wp
   col_config%co2%rh_k_base = 0.01_wp
   col_config%fast_soil_carbon = 5.0_wp
   call apply_hydraulics_config(cfg%hydraulics, col_config%hydraulics_params)
   call build_leaf_photo_table(cfg, col_config%leaf_photo)
   col_config%integrator = build_integrator_opts(cfg)
   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   print '(a)', '[test_column_rk45]'

   !=== A. GPP parity on step 1 (identical initial state; RK45 pre-pass == split pre-pass). ====!
   call set_noon_forcing()
   call reset_state()
   cfg%time_integrator = INTEG_ARK
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
   gpp_split = gpp_coh
   call reset_state()
   cfg%time_integrator = INTEG_RK45
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
   gpp_rk45 = gpp_coh
   call ck(abs(gpp_rk45(1) - gpp_split(1)) < 1.0e-12_wp,                                          &
           'RK45 pre-pass gpp bit-identical to the split (build_column_frozen)', abs(gpp_rk45(1) - gpp_split(1)))
   call ck(gpp_rk45(1) > 0.0_wp, 'RK45 midday gpp > 0', gpp_rk45(1))

   !=== B. A dry-window march under INTEG_RK45 stays physical + bounded + sub-saturated. ========!
   call reset_state()
   cfg%time_integrator = INTEG_RK45
   col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
   physical = .true. ; worst_super = -1.0_wp
   do is = 1_ik, 24_ik
      call set_diurnal_forcing(is)
      call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      tcas = biophys%cas%can_temp
      qsat = sat_specific_humidity(tcas, aenv%press)
      worst_super = max(worst_super, biophys%cas%can_shv - qsat)
      physical = physical .and. tcas > 270.0_wp .and. tcas < 325.0_wp .and. biophys%cas%can_shv > 0.0_wp
      physical = physical .and. biophys%leaf_temp(1) > 250.0_wp .and. biophys%leaf_temp(1) < 340.0_wp
      do k = 1_ik, nsl
         physical = physical .and. biophys%soil_w%theta(k) > 0.0_wp .and. biophys%soil_w%theta(k) < 0.6_wp
         physical = physical .and. biophys%soil_e%soil_temp(k) > 260.0_wp .and. biophys%soil_e%soil_temp(k) < 340.0_wp
      end do
      physical = physical .and. biophys%leaf_water_mass(1) > 0.0_wp .and. biophys%wood_water_mass(1) > 0.0_wp
      psi_leaf_diag = psi_from_water_content(biophys%leaf_water_mass(1), col_config%hydraulics_params%leaf_pi0,       &
           col_config%hydraulics_params%leaf_elastic_mod, col_config%hydraulics_params%leaf_apoplast_frac, &
                col_config%hydraulics_params%leaf_water_sat, &
           col_cohort%bleaf(1))
      physical = physical .and. psi_leaf_diag < 0.5_wp .and. psi_leaf_diag > -12.0_wp
   end do
   call ck(physical, 'INTEG_RK45 dry-window march stays physical + bounded (24 steps)', biophys%cas%can_temp)
   call ck(worst_super <= 1.0e-4_wp, 'INTEG_RK45 CAS stays sub-saturated', worst_super)

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
   !       e_runof/e_drain -> root_heat_sink) under a REAL, substantial rainfall rate, not just the         !
   !       incidental background drainage that first caught it. ===================================!
   call test_rk45_budgets_wet()

   !=== E. CANOPY-SURFACE WATER (opt-in, MEDS_ED2_RK45_DESIGN.md sec 3.4, P2c): a diurnal march with   !
   !       a morning rain pulse actually gets intercepted; whole_water closes exactly, whole_energy        !
   !       stays BOUNDED (known deferred sensible-heat approx, same category as sec 2's qloss/qwflux_wl     !
   !       upwind-temperature approximation; mirrors test_column_dynamics.f90's own RUN 6 for the split      !
   !       path and test_column_ark.f90's Test G for ARK). Proves the WIRING end to end for RK45's           !
   !       genuinely-integrated surf_water ODE (not the wetted-fraction algebra itself, already unit-         !
   !       tested in test_surface_energy.f90). =========================================================!
   call test_rk45_canopy_water()

   !=== F. LEAF/ROOT-TURNOVER SHED WATER (P4, MEDS_ED2_RK45_DESIGN.md): a constant shed_water_rate    !
   !       (distinct from rainfall) wets the soil and both whole_water/whole_energy still close. ========!
   call test_rk45_shed_water()
   call test_rk45_saturated()

   !=== F2. DRYDOWN to the RESIDUAL floor: the OTHER edge of the constitutive domain. RK45 removed the  !
   !       scratch solve's floor ENTHALPY compensation (#78 item 3) because that mass never moves on    !
   !       this trajectory, which left theta_res unguarded on the committed state -- so the commit      !
   !       guard now applies the floor itself and books the created water as a boundary input. This is  !
   !       the only scenario in the file that pushes a layer down onto that floor. =====================!
   call test_rk45_drydown()

   !=== G. TINY NEWLY-RECRUITED COHORT: LAI/WAI at a just-recruited cohort's real magnitude (an        !
   !       actual 30-yr Ithaca run's first recruitment event) reproduces two RK45-only failures a        !
   !       mature-canopy cohort (LAI=3/WAI=0.5, tests A-F above) never exercises: (1) leaf/wood_temp       !
   !       overflowing via a near-zero veg_energy_balance denominator (fixed by the coupling floor,     !
   !       meds_plant_biophysics.f90), and (2) soil theta/soil_energy escaping their physical          !
   !       domain within a single explicit stage for a near-bare patch, crashing ground_evaporation's       !
   !       fractional pow() (fixed by clamp_theta/clamp_cas/clamp_soil_energy in every RK45 stage, plus      !
   !       the committed y_out, meds_fast_rk45.f90). Neither fix is exercised by a mature canopy, whose      !
   !       coupling/heat-capacity terms sit far above the floor and whose surface<->soil coupling is far     !
   !       from stiff -- this test is the ONLY one in the suite with a cohort small enough to trip either.  !
   call test_rk45_tiny_cohort()

   !=== H. DENSE COLD CANOPY (P6, MEDS_ED2_RK45_DESIGN.md): a mature canopy (LAI=3) under a hard cold    !
   !       snap makes the diagnostic leaf<->CAS coupling stiff enough that the FULLY-EXPLICIT RK45        !
   !       surface oscillates and rails the CAS/soil to the clamp bounds within one dt_fast -- the        !
   !       instability that collapsed the 30-yr Ithaca run at LAI~1 in winter 2050. The dispatcher's      !
   !       hybrid rescue must detect the railed commit and REDO the step on the stable implicit-CAS       !
   !       split path, keeping the state physical. Asserts: the rescue actually fires (budget%rk45_rescue   !
   !       > 0, i.e. the test exercises the fix, not a no-op), and CAS/soil/leaf stay in [250,340] K      !
   !       with GPP finite throughout. ================================================================!
   call test_rk45_dense_cold_canopy()
   call test_rk45_bedrock_and_aquifer()
   call test_rk45_prognostic_wood()

   !=== I. REVIEW 2026-09 (item 2 #1): the dispatcher must report psi_leaf_coh on the RK45 success !
   !       path too. It used to fill it only after the RK45 block, whose success branch returns      !
   !       early, so under time_integrator=rk45 the daily-max accumulator read whatever the caller's  !
   !       per-thread buffer held (another patch, or uninitialised memory) and beta_stomata was       !
   !       computed from garbage. Seed the buffer with an impossible POSITIVE potential and assert     !
   !       the step overwrote it with a physical (<= 0) value without an ARK rescue having run. ======!
   call test_rk45_reports_psi_leaf()

   !=== J. REVIEW 2026-09 (item 2 #3): on moderately DRY soil the wood<->soil interface must cancel !
   !       to machine precision. frozen%roots%uptake is advance_soil_water_column's realized supply, which already  !
   !       carries the psi-wilting ramp; the RK45 RHS used to pass it back through the ramp, so the    !
   !       soil lost fwilt*uptake while wood gained uptake -- water created from nothing whenever       !
   !       psi_open > psi_soil > psi_wilt. The existing dry-down (theta ~ theta_res, fwilt ~ 0) and     !
   !       moist tests (fwilt = 1) both sit where the product uptake*(1-fwilt) vanishes; this one sits  !
   !       in the middle of the ramp, where it does not, and asks for closure at the ledger's own       !
   !       round-off rather than at its 1e-4 relative tolerance. ======================================!
   call test_rk45_dry_uptake_seam()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_rk45: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_rk45: ', nfail, ' FAILED' ; error stop 1
   end if

contains

   !----- march 96 sub-steps (24 h) of INTEG_RK45 over MOIST free-draining soil and assert the      !
   !      whole-column water + energy budgets close (design doc sec 8 gates 2/3). -----------------!
   !----- PHASE 0/2/3 (MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN.md): the ARK/RK45 bottom-BC guard is    !
   !      GONE. All three BCs are now pure boundary conditions with no prognostic state behind them, !
   !      so RK45 -- which integrates its own theta -- has nothing left to borrow from the scratch   !
   !      solve. Bedrock seals the face; the aquifer is head-driven and TWO-WAY, and RK45 gets it     !
   !      through soil_water_time_deriv like any other face. ---------------------------------------!
   !----- PHASE 4: prognostic WOOD on RK45, via the SAME shared advance_wood_energy_full the ARK uses !
   !      (test_column_ark's test_wood_prognostic carries the full rationale and the bsap dual-purpose !
   !      trap). Keeping both adaptive schemes on one wood model is the point of the phase. -----------!
   subroutine test_rk45_prognostic_wood()
      integer(ik) :: istep
      real(wp)    :: dmax_lag
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      dmax_lag = 0.0_wp
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         dmax_lag = max(dmax_lag, abs(biophys%wood_temp(1) - biophys%cas%can_temp))
      end do
      call ck(budget%whole_energy%n_fail == 0_ik, 'RK45 PROG-WOOD: whole_energy closes',              &
              real(budget%whole_energy%n_fail, wp))
      call ck(budget%whole_water%n_fail == 0_ik, 'RK45 PROG-WOOD: whole_water closes',                &
              real(budget%whole_water%n_fail, wp))
      call ck(dmax_lag > 1.0e-3_wp, 'RK45 PROG-WOOD: wood temperature lags the CAS', dmax_lag)
      call ck(biophys%wood_temp(1) > 200.0_wp .and. biophys%wood_temp(1) < 350.0_wp,                        &
              'RK45 PROG-WOOD: wood temperature physical', biophys%wood_temp(1))
      print '(a,i0,a,i0)', '   RK45 PROG-WOOD last dt_fast: substeps = ', budget%integ_nsteps,        &
            ' , rescues = ', budget%rk45_rescue
   end subroutine test_rk45_prognostic_wood

   subroutine test_rk45_bedrock_and_aquifer()
      integer(ik) :: istep
      real(wp)    :: theta_bot0
      !----- (a) sealed bedrock column on RK45. -----------------------------------------------!
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%soil_water_opts%bottom_bc = SOIL_BC_BEDROCK
      do istep = 1_ik, 48_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      end do
      call ck(budget%whole_water%n_fail == 0_ik, 'BEDROCK/RK45: whole_water closes',                 &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_energy%n_fail == 0_ik, 'BEDROCK/RK45: whole_energy closes',               &
              real(budget%whole_energy%n_fail, wp))

      !----- (b) aquifer BC on RK45: it used to hard error-stop here. A column started DRY must wet  !
      !          from below, and both ledgers must close with the bottom flux running UPWARD -- the   !
      !          direction that did not exist before Phase 0. -----------------------------------!
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%soil_water_opts%bottom_bc = SOIL_BC_AQUIFER
      biophys%soil_w%theta(1:col_config%soil%n_active) = 0.15_wp
      theta_bot0 = biophys%soil_w%theta(col_config%soil%n_active)
      do istep = 1_ik, 48_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      end do
      call ck(budget%whole_water%n_fail == 0_ik, 'AQUIFER/RK45: whole_water closes',                 &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_energy%n_fail == 0_ik, 'AQUIFER/RK45: whole_energy closes',               &
              real(budget%whole_energy%n_fail, wp))
      call ck(biophys%soil_w%theta(col_config%soil%n_active) > theta_bot0,                                    &
              'AQUIFER/RK45: dry column wets from below', biophys%soil_w%theta(col_config%soil%n_active) - theta_bot0)
      !----- The plan's open risk: K_bot/Delta with Delta = dz(n)/2 is a fast boundary term the       !
      !      implicit paths absorb and an explicit march may not. Report the cost rather than assume. !
      print '(a,i0,a,i0,a,i0)', '   AQUIFER/RK45 last dt_fast: substeps = ', budget%integ_nsteps,     &
            ' , rejects = ', budget%integ_nrej, ' , rescues = ', budget%rk45_rescue
      col_config%soil_water_opts%bottom_bc = SOIL_BC_FREE_DRAIN
   end subroutine test_rk45_bedrock_and_aquifer

   subroutine test_rk45_budgets()
      integer(ik) :: istep
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      end do
      call ck(budget%whole_water%n_fail == 0_ik, 'RK45: whole_water closes (n_fail==0)',            &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_energy%n_fail == 0_ik, 'RK45: whole_energy closes (n_fail==0)',          &
              real(budget%whole_energy%n_fail, wp))
      call ck(budget%whole_water%n_check == 576_ik, 'RK45: ledger fired every dispatched dt_fast',   &
              real(budget%whole_water%n_check, wp))
      call ck(budget%whole_energy%worst < 1.0_wp, 'RK45: whole-energy closes < 1 J',                &
              budget%whole_energy%worst)
      call ck(budget%integ_nsteps >= 1_ik .and. budget%integ_nsteps <= 64_ik,                          &
              'RK45: adaptive substeps bounded (sec 6: ~3 expected, cap 64)', real(budget%integ_nsteps, wp))
      print '(a,i0,a,i0)', '   RK45 last dt_fast: substeps = ', budget%integ_nsteps,                 &
            ' , rejects = ', budget%integ_nrej
      print '(a,es10.3,a,es10.3,a)', '   (worst whole-column resid: energy= ', budget%whole_energy%worst, &
            ' J/m2  water= ', budget%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets

   !----- march 96 sub-steps (24 h) of INTEG_RK45 over free-draining soil WITH continuous rain       !
   !      (rainfall>0), mirroring test_column_ark's test_ark_budgets_wet. Asserts the run completes,   !
   !      the soil wets, and BOTH whole-column budgets close at the tight (non-inflated) tolerance.   !
   subroutine test_rk45_budgets_wet()
      integer(ik) :: istep, commit_n
      real(wp)    :: theta_col0, theta_col1, commit_mass, commit_energy
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      theta_col0 = sum(biophys%soil_w%theta(1:nsl))
      commit_n = 0_ik ; commit_mass = 0.0_wp ; commit_energy = 0.0_wp
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-5_wp                         ! ~0.29 mm/hr continuous rain (rainfall>0)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         !----- the clamp counters are per-sub-step (each stepper zeroes them on entry), so a test that  !
         !      wants a window total has to accumulate, exactly as the site-level driver does. ----------!
         commit_n      = commit_n      + budget%clamp_commit_n
         commit_mass   = commit_mass   + budget%clamp_mass
         commit_energy = commit_energy + budget%clamp_energy
      end do
      theta_col1 = sum(biophys%soil_w%theta(1:nsl))
      call ck(budget%whole_water%n_check == 576_ik, 'RK45 wet: ran 96 wet steps (no guard error stop)',   &
              real(budget%whole_water%n_check, wp))
      call ck(theta_col1 > theta_col0, 'RK45 wet: rain wetted the soil column (theta rose)',           &
              theta_col1 - theta_col0)
      call ck(budget%whole_energy%n_fail == 0_ik, 'RK45 wet: whole_energy closes (n_fail==0)',           &
              real(budget%whole_energy%n_fail, wp))
      call ck(budget%whole_water%n_fail == 0_ik, 'RK45 wet: whole_water closes (n_fail==0)',             &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_energy%worst < 1.0_wp, 'RK45 wet: whole-energy closes < 1 J',                 &
              budget%whole_energy%worst)
      !----- WHY the books close here, pinned as a mechanism rather than an outcome. Measuring the      !
      !      clamps established something the closure assertions alone cannot show, and that was NOT    !
      !      the expected result: the committed-state clamp fires on essentially EVERY sub-step in      !
      !      this benign window too (~1250 activations over 96 steps), because the Richards solve       !
      !      routinely lands a whisker outside [theta_res, theta_sat]. So an activation COUNT does not  !
      !      separate a healthy window from a broken one -- both fire constantly. The MAGNITUDE does,   !
      !      by ~7 orders of magnitude: ~3e-5 kg/m2 cumulative here against ~1.4e2 kg/m2 in the         !
      !      saturated twin. Assert on the magnitude, therefore, and let the count be telemetry. -------!
      call ck(commit_mass + commit_energy < 1.0e-3_wp,                                                &
              'RK45 wet: commit-clamp corrections stay negligible (closure is not luck)',             &
              commit_mass + commit_energy)
      print '(a,i0,a,es10.3,a,es10.3,a)', '   (RK45 wet commit clamps: n= ', commit_n,                 &
            '  unbookkept mass= ', commit_mass, ' kg/m2  energy= ', commit_energy, ' J/m2)'
      print '(a,es10.3,a,es10.3,a)', '   (RK45 wet worst whole-column resid: energy= ',                &
            budget%whole_energy%worst, ' J/m2  water= ', budget%whole_water%worst, ' kg/m2)'
   end subroutine test_rk45_budgets_wet

   !----- march 96 sub-steps (24 h) of INTEG_RK45 with canopy_water_on and a 5-step morning rain pulse   !
   !      (istep 20-24); assert some of it was intercepted and the budgets close (mirrors                !
   !      test_column_dynamics.f90's own RUN 6 / test_column_ark.f90's Test G). -----------------------!
   subroutine test_rk45_canopy_water()
      integer(ik) :: istep
      real(wp)    :: surf_water_peak
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      col_config%canopy_water_on = .true.
      surf_water_peak = 0.0_wp
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         if (istep >= 20_ik .and. istep <= 24_ik) forc%rainfall = 5.0e-5_wp   ! a morning rain pulse
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         surf_water_peak = max(surf_water_peak, biophys%leaf_surf_water(1) + biophys%wood_surf_water(1))
      end do
      col_config%canopy_water_on = .false.   ! restore default for any test added after this
      call ck(budget%whole_water%n_fail == 0_ik, 'RK45 canopy water: whole-column water still closes',   &
              real(budget%whole_water%n_fail, wp))
      call ck(surf_water_peak > 0.0_wp, 'RK45 canopy water: the morning rain pulse was intercepted',    &
              surf_water_peak)
      call ck(budget%whole_energy%worst < 5.0e6_wp,                                                       &
              'RK45 canopy water: whole-column energy stays BOUNDED (known deferred approx)',           &
              budget%whole_energy%worst)
      print '(a,es10.3,a)', '   (RK45 canopy water peak film=', surf_water_peak, ' kg/m2)'
   end subroutine test_rk45_canopy_water

   !----- march 96 sub-steps (24 h) of INTEG_RK45 with a constant leaf/root-turnover shed-water rate  !
   !      (MEDS_ED2_RK45_DESIGN.md P4, biophys%shed_water_rate -- a PATCH-level input, not atmospheric      !
   !      forcing, so it is frozen on biophys for the whole day rather than living on forc; distinct from    !
   !      rainfall, which stays 0 throughout), mirroring test_column_ark's test_ark_shed_water: the soil    !
   !      must wet from THIS input alone, and both whole_water AND whole_energy must still close at        !
   !      RK45's own tight (non-split-inflated) tolerance -- energy closing needs no separate wiring         !
   !      (rides the SAME e_infil/t_film_valuation treatment every other infiltrating input already gets). ---------!
   subroutine test_rk45_shed_water()
      integer(ik) :: istep
      real(wp)    :: theta_col0, theta_col1
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      biophys%shed_water_rate = 8.0e-5_wp                     ! P4: frozen for the whole day (rainfall stays 0)
      theta_col0 = sum(biophys%soil_w%theta(1:nsl))
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      end do
      theta_col1 = sum(biophys%soil_w%theta(1:nsl))
      biophys%shed_water_rate = 0.0_wp   ! restore default for any test added after this
      call ck(theta_col1 > theta_col0,                                                              &
              'RK45 shed water: leaf/root shed water alone wetted the soil column (theta rose)',     &
              theta_col1 - theta_col0)
      call ck(budget%whole_water%n_fail == 0_ik,                                                       &
              'RK45 shed water: whole-column WATER still closes with shed_water_rate active',        &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_energy%n_fail == 0_ik,                                                      &
              'RK45 shed water: whole-column ENERGY still closes (no separate energy wiring needed)', &
              real(budget%whole_energy%n_fail, wp))
      call ck(budget%whole_energy%worst < 1.0_wp, 'RK45 shed water: whole-energy closes < 1 J',        &
              budget%whole_energy%worst)
   end subroutine test_rk45_shed_water

   !----- march 48 sub-steps (12 h, spanning full daylight) of INTEG_RK45 with the cohort shrunk to a   !
   !      just-recruited size (LAI/WAI/nplant/bleaf/bsap/broot at the actual magnitude the 30-yr Ithaca   !
   !      run's first recruitment event produced, MEDS_ED2_RK45_DESIGN.md bug report). Before the two      !
   !      fixes above, this reproduced BOTH failures directly: leaf/wood_temp overflowing (T**4 in the     !
   !      RT forcing) and, once that was patched alone, a SEPARATE ground_evaporation crash from theta/     !
   !      soil_energy escaping their domain. Restores the shared col_cohort/biophys state on exit. -------------------!
   subroutine test_rk45_tiny_cohort()
      integer(ik) :: istep, k
      real(wp) :: lai0, wai0, leaf_area0, bleaf0, bsap0, broot0, sap_area0, nplant0
      logical  :: physical
      lai0 = col_cohort%lai(1) ; wai0 = col_cohort%wai(1) ; leaf_area0 = col_cohort%leaf_area(1)
      bleaf0 = col_cohort%bleaf(1) ; bsap0 = col_cohort%bsap(1) ; broot0 = col_cohort%broot(1)
      sap_area0 = col_cohort%sap_area(1) ; nplant0 = col_cohort%nplant(1)

      col_cohort%lai(1) = 0.00265_wp ; col_cohort%wai(1) = 0.000529_wp ; col_cohort%leaf_area(1) = 0.00265_wp
      col_cohort%bleaf(1) = 0.0203507_wp ; col_cohort%bsap(1) = 0.00212943_wp ; col_cohort%broot(1) = 0.0203507_wp
      col_cohort%sap_area(1) = 1.0e-4_wp ; col_cohort%nplant(1) = 0.01_wp

      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      physical = .true.
      do istep = 1_ik, 48_ik
         call set_diurnal_forcing(istep)
         !----- The shared fixture hands EVERY cohort a full-canopy absorbed load; for a near-zero-LAI  !
         !      cohort that is unphysical (the real two-stream RT, apply_rt_forcing, scales absorption   !
         !      by leaf area, so a 0.0027-LAI seedling intercepts ~0.1% of the beam). Apply the same     !
         !      Beer-law scaling here, else the test asserts on a leaf absorbing 500 W/m2 through        !
         !      essentially zero area -- which genuinely MUST run hot, and which the pre-P6 discontinuous !
         !      floor only hid by zeroing the whole balance. -------------------------------------------!
         block
            real(wp) :: rad_frac
            rad_frac      = 1.0_wp - exp(-0.5_wp * col_cohort%lai(1))
            forc%abs_sw   = forc%abs_sw  * rad_frac
            forc%abs_par  = forc%abs_sw
            forc%abs_lw   = forc%abs_lw  * rad_frac
         end block
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         physical = physical .and. biophys%leaf_temp(1) > 250.0_wp .and. biophys%leaf_temp(1) < 340.0_wp
         physical = physical .and. biophys%wood_temp(1) > 250.0_wp .and. biophys%wood_temp(1) < 340.0_wp
         physical = physical .and. biophys%cas%can_temp > 250.0_wp .and. biophys%cas%can_temp < 340.0_wp
         do k = 1_ik, nsl
            physical = physical .and. biophys%soil_w%theta(k) >= 0.0_wp .and. biophys%soil_w%theta(k) <= 0.6_wp
            physical = physical .and. biophys%soil_e%soil_temp(k) > 250.0_wp .and. biophys%soil_e%soil_temp(k) < 340.0_wp
         end do
         physical = physical .and. gpp_coh(1) == gpp_coh(1)   ! not NaN
      end do
      call ck(physical, 'RK45 tiny (just-recruited) cohort stays physical + bounded (48 steps)',   &
              biophys%leaf_temp(1))
      call ck(gpp_coh(1) == gpp_coh(1), 'RK45 tiny cohort: gpp is not NaN', gpp_coh(1))

      col_cohort%lai(1) = lai0 ; col_cohort%wai(1) = wai0 ; col_cohort%leaf_area(1) = leaf_area0
      col_cohort%bleaf(1) = bleaf0 ; col_cohort%bsap(1) = bsap0 ; col_cohort%broot(1) = broot0
      col_cohort%sap_area(1) = sap_area0 ; col_cohort%nplant(1) = nplant0
   end subroutine test_rk45_tiny_cohort

   !----- ROBUSTNESS STRESS (P6, MEDS_ED2_RK45_DESIGN.md): a very dense canopy (LAI=8) under a hard      !
   !      cold snap (dry −20 C air) at the production macro-step (2*dt_fast = 1800 s, the Ithaca step)    !
   !      -- the operating-point FAMILY that rails the explicit RK45 surface at scale. Whether the        !
   !      hybrid rescue fires on this exact synthetic point or the P5 stage/y_out clamps alone hold it,   !
   !      the committed state MUST stay physical (CAS/leaf/soil in [250,340] K, GPP finite): a regression !
   !      that lost surface stability here would rail. (The rescue path itself is proven at integration    !
   !      scale by the 30-yr Ithaca run, which collapses at LAI~1 in winter 2050 WITHOUT the rescue and    !
   !      completes healthy WITH it -- a delicate decades-evolved point no short march reproduces; the     !
   !      rescue count here is reported for information, not asserted.) Restores col_cohort geometry on exit. ----!
   subroutine test_rk45_dense_cold_canopy()
      integer(ik) :: istep, k, total_rescue
      real(wp) :: lai0, wai0, la0
      logical  :: physical
      lai0 = col_cohort%lai(1) ; wai0 = col_cohort%wai(1) ; la0 = col_cohort%leaf_area(1)
      ! very dense -> high coupling gain
      col_cohort%lai(1) = 8.0_wp
      col_cohort%wai(1) = 1.5_wp
      col_cohort%leaf_area(1) = 26.0_wp
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      physical = .true. ; total_rescue = 0_ik
      do istep = 1_ik, 576_ik
         call set_coldsnap_forcing(istep)
         call column_fast_step(2.0_wp*dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
               gpp_coh=gpp_coh)
         total_rescue = total_rescue + budget%rk45_rescue ; budget%rk45_rescue = 0_ik
         physical = physical .and. biophys%cas%can_temp > 250.0_wp .and. biophys%cas%can_temp < 340.0_wp
         physical = physical .and. biophys%leaf_temp(1)  > 250.0_wp .and. biophys%leaf_temp(1)  < 340.0_wp
         do k = 1_ik, nsl
            physical = physical .and. biophys%soil_e%soil_temp(k) > 250.0_wp .and. biophys%soil_e%soil_temp(k) < 340.0_wp
         end do
         physical = physical .and. gpp_coh(1) == gpp_coh(1)   ! not NaN
      end do
      call ck(physical, 'RK45 dense cold-snap canopy stays physical (P6 surface stability)',           &
              biophys%cas%can_temp)
      print '(a,i0,a)', '   (dense cold-snap: ', total_rescue, ' RK45->split rescues over 96 steps)'
      col_cohort%lai(1) = lai0 ; col_cohort%wai(1) = wai0 ; col_cohort%leaf_area(1) = la0
   end subroutine test_rk45_dense_cold_canopy


   !----- SATURATED sealed column: the ARK/RK45 twin of test_column_dynamics RUN 7. A bedrock bottom   !
   !      seeded just below theta_sat under heavy rain both SATURATES (firing meds_soil_water's post-   !
   !      solve theta clip, which moves water with no face) and OVERFLOWS the ponding store (surface     !
   !      runoff). None of the other RK45 runs reach either state -- they are free-draining at           !
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
   !----- Drive the column down toward theta_res: seed just above it, no rainfall, and a dry atmosphere   !
   !      so ground evaporation and transpiration both pull hard. The assertion is the DOMAIN INVARIANT  !
   !      -- theta never commits below theta_res, so soil_psi_from_theta is never evaluated at Se < 0    !
   !      on the following step -- plus closure of both whole-column ledgers.                            !
   !                                                                                                    !
   !      MEASURED, and worth recording rather than implying otherwise: over 48 h this bottoms out at    !
   !      theta = 0.0807 against theta_res = 0.078, so the commit floor never actually fires. Both       !
   !      sinks self-limit before it: the root sink is psi-limited (fwilt shuts it off) and ground        !
   !      evaporation is capped at (theta(1)-theta_res)*dz*rho_w/dt inside the scratch solve, so the      !
   !      frozen q_top cannot over-extract either. The floor is therefore a DEFENSIVE guard on an        !
   !      invariant that nothing else on this path enforces (RK45's old frozen%hydrology%floor_enth compensation was   !
   !      an enthalpy term, never a state edit), and its ledger booking is UNEXERCISED by this suite --   !
   !      a mutation that unbooks floor_mass_rk/floor_enth_rk passes. The invariant below is the right    !
   !      assertion regardless of which mechanism happens to enforce it, and this drydown regime is      !
   !      covered by nothing else in the file. -------------------------------------------------------!
   subroutine test_rk45_drydown()
      integer(ik) :: istep
      real(wp)    :: theta_min, res_min, e_worst, w_worst
      theta_seed = 0.085_wp                              ! just above theta_res = 0.078
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      theta_min = 1.0e9_wp
      res_min   = minval(col_config%soil%theta_res(1:nsl))
      do istep = 1_ik, 192_ik
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         forc%shv_atm = 0.5_wp * forc%shv_atm            ! halve the atmospheric humidity: dry the air
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         theta_min = min(theta_min, minval(biophys%soil_w%theta(1:nsl)))
      end do
      e_worst = budget%whole_energy%worst ; w_worst = budget%whole_water%worst
      theta_seed = theta0                                 ! restore for any test added after this
      call ck(theta_min >= res_min - 1.0e-12_wp,                                                       &
              'RK45 drydown: theta never commits below theta_res (Se >= 0 guaranteed)', theta_min)
      call ck(w_worst < 1.0e-6_wp,                                                                     &
              'RK45 drydown: whole-column WATER closes with the residual floor booked', w_worst)
      call ck(e_worst < 1.0e-3_wp,                                                                     &
              'RK45 drydown: whole-column ENERGY closes with the residual floor booked', e_worst)
      print '(a,es10.3,a,es10.3,a)', '   (RK45 drydown: theta_min= ', theta_min, '  theta_res= ',      &
            res_min, ')'
   end subroutine test_rk45_drydown

   subroutine test_rk45_saturated()
      integer(ik) :: istep, commit_n
      real(wp)    :: pond_peak, theta_peak, ss_min, ss_max, commit_mass, commit_energy, ood_peak
      theta_seed = 0.428_wp                              ! just below theta_sat = 0.43
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      pond_peak = 0.0_wp ; theta_peak = 0.0_wp
      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp
      commit_n = 0_ik ; commit_mass = 0.0_wp ; commit_energy = 0.0_wp ; ood_peak = 0.0_wp
      do istep = 1_ik, 576_ik
         call set_diurnal_forcing(istep)
         forc%rainfall = 8.0e-3_wp                         ! ~29 mm/hr: far above the drainage capacity
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
         pond_peak  = max(pond_peak,  biophys%soil_w%w_surface)
         theta_peak = max(theta_peak, maxval(biophys%soil_w%theta(1:nsl)))
         ss_min = min(ss_min, biophys%soil_e%soil_temp(1)) ; ss_max = max(ss_max, biophys%soil_e%soil_temp(1))
         ood_peak      = max(ood_peak,   budget%theta_ood_max)
         commit_n      = commit_n      + budget%clamp_commit_n
         commit_mass   = commit_mass   + budget%clamp_mass
         commit_energy = commit_energy + budget%clamp_energy
      end do
      theta_seed = theta0                                 ! restore for any test added after this
      call ck(theta_peak >= 0.43_wp - 1.0e-12_wp,                                                    &
              'RK45 saturated: column reached theta_sat (clip path is live)', theta_peak)
      call ck(pond_peak >= col_config%soil_water_opts%w_pond_max - 1.0e-9_wp,                                        &
              'RK45 saturated: ponding store overflowed (runoff path is live)', pond_peak)
      !----- ENERGY is asserted as a BOUND, not exact closure, and the reason is a VERIFIED defect in  !
      !      RK45's own stability guard rather than in the water-enthalpy wiring this test covers.     !
      !      clamp_soil_energy rebuilds soil_energy at a clamped temperature with NO ledger term, so    !
      !      when it bites on the COMMITTED state it silently injects energy. Measured: disabling that  !
      !      one call drops the residual here from 1.5e5 to 5.5e-7 J/m2 and the soil-surface peak from  !
      !      329.5 K to 303.4 K. It is trajectory-dependent -- ifx never triggers it in this scenario   !
      !      (8.8e-7 J/m2) while nvfortran does -- which is exactly why a bound, not n_fail, belongs     !
      !      here. clamp_theta has the same unbookkept character on the water side (~0.85 kg/m2). Both   !
      !      are the same "correction with no bookkeeping" family as the clip/floor fixes, but they      !
      !      live in the explicit integrator's guards and belong with the deferred RK45 work. The        !
      !      UNSATURATED path closes exactly on both compilers (test_rk45_budgets_wet, n_fail == 0),     !
      !      which is what shows the enthalpy plumbing itself is right.                                 !
      !      (#78 item 3 note: with the interior faces wired to RK45's own theta the guard no longer     !
      !      bites at all here -- the residual is ~6e-7 J/m2 on ifx -- but the bound stays, because      !
      !      whether it bites is trajectory- and compiler-dependent, and that dependence is the finding.)!
      call ck(budget%whole_energy%worst < 1.0e-3_wp,                                                    &
              'RK45 saturated: whole-column ENERGY closes (C1 removed the unbookkept clamp)',      &
              budget%whole_energy%worst)
      !----- WATER now closes to MACHINE PRECISION here, which it did not for most of this file's    !
      !      history. The trajectory is worth keeping because each step was a distinct defect:         !
      !        4.35  -- RK45 committed its own theta but took drainage, ponding and runoff from the    !
      !                 frozen Act-1 scratch solve                                                     !
      !        5.86  -- C2 made DRAINAGE state-consistent while runoff was not, so the imbalance MOVED !
      !                 rather than shrank (the expected cost of fixing half a pair)                   !
      !        3.33  -- the residual saturation clip closed that pair, routing RK45's own excess to    !
      !                 the pond with paired enthalpy                                                  !
      !        1e-13 -- the pond is rebuilt from RK45's OWN trajectory. frozen%hydrology%w_surface1 already held    !
      !                 the SCRATCH solve's clip, mass this theta never shed, so adding RK45's clip    !
      !                 on top counted that water twice (issue #75).                                   !
      !      Assert closure, not a bound: there is no longer a known gap to tolerate. -----------------!
      call ck(budget%whole_water%n_fail == 0_ik,                                                       &
              'RK45 saturated: whole-column WATER closes (n_fail==0)',                                &
              real(budget%whole_water%n_fail, wp))
      call ck(budget%whole_water%worst < 1.0e-6_wp,                                                    &
              'RK45 saturated: whole-column WATER closes to machine precision',                       &
              budget%whole_water%worst)
      !----- SOIL-SURFACE TEMPERATURE is the one assertion that names #78 item 3's mechanism, so it is  !
      !      bounded at 310 K rather than the 340 K that used to be needed. The history is worth the      !
      !      lines, because two defects of matched magnitude and OPPOSITE sign lived here:                !
      !        285.3 K -- column_derivs advected interior soil enthalpy on the SCRATCH solve's frozen     !
      !                   faces while committing RK45's own theta (a ~2.6e6 J/m2/step vertical            !
      !                   misplacement, invisible to any whole-column ledger), and separately took the    !
      !                   scratch's clip ENTHALPY as a stage sink for mass this trajectory never shed     !
      !                   (~2.6e6 J/m2/step of spurious cooling). The two very nearly cancelled.          !
      !        345.0 K -- removing ONLY the borrowed clip enthalpy, which is what exposed the faces.      !
      !        295.1 K -- both fixed: the faces now carry the flux this stage's own theta is using.       !
      !      A 29 mm/hr rain event cannot warm a soil surface past ~300 K, so 310 K is a bound with       !
      !      real diagnostic power; 340 K only ever tolerated the defect. -------------------------------!
      call ck(ss_min > 250.0_wp .and. ss_max < 310.0_wp,                                             &
              'RK45 saturated: soil surface temp stays physical (interior faces on OWN theta)', ss_max)
      !----- C1 left NOTHING corrected off-ledger here; assert the mechanism is gone rather than only  !
      !      that the residuals are small, since a bound alone would still pass if the clamp were       !
      !      silently replaced by some other unbookkept                                                 !
      !      correction, whereas this names the mechanism. The magnitudes are printed rather than       !
      !      pinned because they are trajectory- (and therefore compiler-) dependent -- that            !
      !      dependence is itself the finding, and Phase C removes the commit clamp outright, at which  !
      !      point commit_n here must fall to 0 and these residuals collapse. -------------------------!
      !----- C1: the committed state is no longer clamped, so NOTHING is corrected outside the ledger. !
      !      Before the fix this window logged 142.7 kg/m2 of unbookkept water and (on nvfortran)      !
      !      1.5e5 J/m2 of unbookkept energy, and the soil surface peaked at 329.4 K; it now peaks at  !
      !      297.0 K. Assert the mechanism is gone, not just that the residuals are small. ------------!
      call ck(commit_n == 0_ik .and. commit_mass == 0.0_wp .and. commit_energy == 0.0_wp,            &
              'RK45 saturated: no committed-state clamp (C1 -- nothing corrected off-ledger)',        &
              commit_mass + commit_energy)
      !----- HONEST CONSEQUENCE of C1, asserted so it cannot drift unnoticed: with the commit clamp    !
      !      gone, theta is now committed slightly ABOVE theta_sat (0.4536 vs 0.43) on this sealed     !
      !      saturated column. That is the SAME error the clamp was hiding, now visible and on the     !
      !      books instead of silently paid for in fabricated mass. It is C2's to remove -- RK45 takes !
      !      ponding/drainage/runoff from the frozen scratch solve while integrating its own theta, so !
      !      the water that should have left as runoff has nowhere to go. Bound it so a REGRESSION     !
      !      still trips, and tighten this once C2 lands. ------------------------------------------!
      !----- C2 removed the overshoot entirely: the residual saturation clip now routes RK45's own     !
      !      excess to the ponding store with paired enthalpy, so theta commits AT theta_sat rather     !
      !      than above it (0.4382 -> 0.43000). Assert equality-to-tolerance, not just a bound. -------!
      call ck(theta_peak <= 0.43_wp + 1.0e-9_wp,                                                     &
              'RK45 saturated: theta commits AT theta_sat, no overshoot (C2)', theta_peak)
      !----- ISSUE #78 ITEM 2, bounded rather than merely inspected. RK45's stage-1 RHS reads the      !
      !      previous SUB-step's committed theta, which C1 deliberately stopped clamping, so this is    !
      !      the one place a constitutive kernel sees theta outside [theta_res, theta_sat]. It is       !
      !      HARMLESS to evaluate -- every kernel clamps its own effective saturation, proven by        !
      !      test_column_derivs' test_rhs_domain_safety, and an oversaturated cell is then treated as   !
      !      exactly saturated, which is the physically right answer for one. What is asserted here is  !
      !      that the excursion stays SMALL, because "harmless" rests on the curves being evaluated     !
      !      AT their endpoint: measured 7.9e-3 in theta (Se = 1.022) on this, the wettest fixture in   !
      !      the suite, and identically 0 on a month-long forced Ithaca cell. Bound at 2e-2 -- ~2.5x    !
      !      the measurement, tight enough that a real regression trips it. --------------------------!
      call ck(ood_peak < 2.0e-2_wp,                                                                    &
              'RK45 saturated: constitutive-domain excursion stays small (#78 item 2)', ood_peak)
      print '(a,es10.3)', '   (RK45 saturated: worst constitutive-domain excursion in theta= ', ood_peak
      print '(a,i0,a,es10.3,a,es10.3,a)', '   (RK45 saturated commit clamps: n= ', commit_n,          &
            '  unbookkept mass= ', commit_mass, ' kg/m2  energy= ', commit_energy, ' J/m2)'
   end subroutine test_rk45_saturated

   subroutine test_rk45_reports_psi_leaf()
      call set_noon_forcing()
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      psi_leaf_probe = 1.0_wp                            ! impossible: psi_leaf is <= 0 by construction
      call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,           &
                            gpp_coh=gpp_coh, psi_leaf_coh=psi_leaf_probe)
      call ck(budget%rk45_rescue == 0_ik, 'RK45 psi report: the RK45 path itself ran (no ARK rescue)', &
              real(budget%rk45_rescue, wp))
      call ck(psi_leaf_probe(1) <= 0.0_wp .and. psi_leaf_probe(1) > -50.0_wp,                     &
              'RK45 psi report: psi_leaf_coh is filled on the RK45 success path', psi_leaf_probe(1))
   end subroutine test_rk45_reports_psi_leaf

   subroutine test_rk45_dry_uptake_seam()
      integer(ik) :: istep
      real(wp)    :: w_worst, psi_top
      theta_seed = 0.14_wp                               ! psi ~ -10 m: inside the f_wilt ramp
      call reset_state()
      cfg%time_integrator = INTEG_RK45
      col_config%integrator = build_integrator_opts(cfg)   ! the schemes read the record, not cfg
      do istep = 1_ik, 48_ik
         call set_diurnal_forcing(istep)
         forc%rainfall = 0.0_wp
         call column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, gpp_coh=gpp_coh)
      end do
      w_worst = budget%whole_water%worst
      psi_top = soil_psi_from_theta(col_config%soil%retention, biophys%soil_w%theta(1), col_config%soil%theta_sat(1), &
                                    col_config%soil%theta_res(1), col_config%soil%vg_alpha(1), col_config%soil%vg_n(1))
      theta_seed = theta0
      call ck(psi_top < col_config%soil_water_opts%psi_open .and. psi_top > col_config%soil_water_opts%psi_wilt,        &
              'RK45 dry seam: the fixture really sits inside the wilting ramp (psi_open > psi > psi_wilt)', psi_top)
      call ck(budget%rk45_rescue == 0_ik, 'RK45 dry seam: the RK45 path itself ran (no ARK rescue)',  &
              real(budget%rk45_rescue, wp))
      call ck(w_worst < 1.0e-8_wp,                                                                  &
              'RK45 dry seam: whole-column WATER closes to round-off (soil debit == wood credit)', w_worst)
   end subroutine test_rk45_dry_uptake_seam

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
      if (allocated(biophys%leaf_temp)) deallocate(biophys%leaf_temp)
      call alloc_patch_biophys(biophys, n, t0, 0.008_wp, 400.0_wp, t0)
      biophys%leaf_water_mass(1:n) = water_content(PSI_INIT, col_config%hydraulics_params%leaf_pi0, &
                              col_config%hydraulics_params%leaf_elastic_mod, &
           col_config%hydraulics_params%leaf_apoplast_frac, col_config%hydraulics_params%leaf_water_sat, col_cohort%bleaf(1:n))
      biophys%wood_water_mass(1:n) = water_content(PSI_INIT, col_config%hydraulics_params%wood_pi0, &
                              col_config%hydraulics_params%wood_elastic_mod, &
           col_config%hydraulics_params%wood_apoplast_frac, col_config%hydraulics_params%wood_water_sat, &
                col_cohort%bsap(1:n) + col_cohort%broot(1:n))
      budget = column_budget_t()
      biophys%soil_w%theta(1:nsl) = theta_seed
      do kk = 1_ik, nsl
         biophys%soil_e%soil_energy(kk) = temp_to_internal_energy(col_config%soil_thermal%soil_dry_heat_capacity(kk), &
                                      theta_seed * rho_h2o, t0, 1.0_wp)
         biophys%soil_e%soil_temp(kk)   = t0
      end do
   end subroutine reset_state

   subroutine set_noon_forcing()
      forc%abs_sw = 450.0_wp ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 70.0_wp ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(295.0_wp, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, 295.0_wp, forc%shv_atm, forc%co2_atm)   ! #97: else MO sees a fixed 298.15 K
   end subroutine set_noon_forcing

   subroutine set_diurnal_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)
      forc%abs_sw = 500.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm = 0.008_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! #97: else MO sees a fixed 298.15 K
   end subroutine set_diurnal_forcing

   !----- HARD cold snap: very cold, very dry reference air (~248-256 K) with a weak winter sun over the !
   !      warm (t0=288 K) soil the reset_state seeds -- a large canopy-through gradient that makes the    !
   !      dense-canopy leaf<->CAS coupling stiff enough to rail the explicit RK45 surface (P6). ----------!
   subroutine set_coldsnap_forcing(istep)
      integer(ik), intent(in) :: istep
      real(wp) :: t_sec, cosz, t_air
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 252.0_wp + 4.0_wp * (cosz - 0.3_wp)                 ! ~248-256 K (a −20 C cold snap)
      forc%abs_sw = 120.0_wp * cosz ; forc%abs_par = forc%abs_sw ; forc%abs_lw = 0.0_wp
      forc%abs_sw_ground = 20.0_wp * cosz ; forc%abs_lw_ground = 0.0_wp
      forc%rainfall = 0.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.001_wp)   ! very dry cold air
      forc%shv_atm = 0.001_wp ; forc%co2_atm = 400.0_wp
      call set_aero_env_atm(aenv, t_air, forc%shv_atm, forc%co2_atm)   ! #97: else MO sees a fixed 298.15 K
   end subroutine set_coldsnap_forcing

end program test_column_rk45
