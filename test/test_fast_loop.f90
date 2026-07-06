!==========================================================================================!
! test_fast_loop -- the SITE-LEVEL fast-biophysics driver over the state-hub-owned per-patch     !
! reservoirs. Builds a two-patch site (one vegetated, one bare ground), seeds the reservoirs,     !
! and runs the fast loop:                                                                        !
!   1. CONSERVATION: every patch's whole-column energy + water budgets close (n_fail == 0) --     !
!      the fast loop conserves over the reservoirs the state hub owns, for a VEGETATED and a      !
!      BARE (zero-cohort) patch alike.                                                            !
!   2. ACTIVITY: the reservoirs evolve (CAS warms under the absorbed shortwave; soil responds).   !
!   3. STEPPER WIRING + GATE: advance_one_step runs the fast loop only when fast_biophysics_on    !
!      AND a fast context is supplied; with the gate off the reservoirs are untouched.            !
!==========================================================================================!
program test_fast_loop
   use meds_kinds,               only : wp, ik
   use meds_config,              only : meds_config_t, GS_CARBON
   use meds_demography_types,    only : site_t
   use meds_init,                only : init_bare_ground, add_cohort, finalize_init
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_biophysics_types,    only : SOIL_RETENTION_VG
   use meds_fast_loop,           only : fast_context_t, init_fast_reservoirs, run_fast_biophysics
   use meds_stepper,             only : advance_one_step
   use meds_test_support,        only : build_test_config, check, check_close, banner
   implicit none

   integer(ik), parameter :: nsl = 10_ik
   type(meds_config_t)  :: cfg
   type(site_t)         :: site
   type(fast_context_t) :: ctx
   real(wp)    :: we, ww, t_cas0, t_cas1, theta0_1, theta1_1, psi0_leaf, psi1_leaf
   real(wp)    :: cbal0, cbal1
   integer(ik) :: nfail

   call banner('site-level fast-biophysics loop (state-hub reservoirs)')

   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg%dt_fast            = 900.0_wp
   cfg%n_fast_per_slow    = 8_ik            ! 8 x 900 s = 2 h of fast integration per slow step

   !----- Build the (MVP) column config + reference met inside the fast context. -----------!
   call build_soil_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,            &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ctx%ccfg%soil)
   call build_soil_thermal(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%ccfg%soil_thermal)
   ctx%ccfg%wood%is_woody = .true. ; ctx%ccfg%wood%stem_resp_factor25 = 0.06_wp ; ctx%ccfg%wood%agf_bs = 0.7_wp
   ctx%ccfg%root%root_resp_factor25 = 0.30_wp
   ctx%ccfg%co2%rh_k_base = 0.01_wp
   ctx%ccfg%fast_soil_carbon = 5.0_wp
   ctx%ccfg%hydro_p%leaf_pi0 = -1.5_wp ; ctx%ccfg%hydro_p%leaf_eps = 12.0_wp
   ctx%ccfg%hydro_p%leaf_af  = 0.30_wp ; ctx%ccfg%hydro_p%leaf_water_sat = 2.0_wp
   ctx%ccfg%hydro_p%wood_pi0 = -1.0_wp ; ctx%ccfg%hydro_p%wood_eps = 8.0_wp
   ctx%ccfg%hydro_p%wood_af  = 0.20_wp ; ctx%ccfg%hydro_p%wood_water_sat = 1.0_wp
   ctx%ccfg%hydro_p%wood_psi50 = -2.0_wp ; ctx%ccfg%hydro_p%wood_kexp = 2.0_wp
   ctx%ccfg%hydro_p%k_plant_max = 6.0e-4_wp ; ctx%ccfg%hydro_p%wood_kmax = 8.0_wp
   ctx%ccfg%hydro_p%vessel_curl = 1.5_wp
   ctx%ccfg%rhizo_cond = 5.0e-4_wp
   ctx%air_temp = 290.0_wp ; ctx%rad_sw_top = 500.0_wp ; ctx%rad_sw_ground = 75.0_wp
   ctx%theta_init = 0.30_wp ; ctx%soil_temp_init = 288.0_wp

   !----- Two patches: patch 1 vegetated, patch 2 BARE (zero cohorts). ---------------------!
   call init_bare_ground(site, cfg, 2_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)     ! patch 1: one tree cohort
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)

   t_cas0    = site%patch%cas(1)%can_temp
   theta0_1  = site%patch%soil_w(1)%theta(1)
   psi0_leaf = site%cohort%psi(1, 1)                 ! patch-1 cohort-1 leaf-node psi (== PSI_INIT)

   !=== 1+2. Run the fast loop directly; conservation + activity + per-cohort persistence. =!
   call run_fast_biophysics(site, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=nfail)
   t_cas1    = site%patch%cas(1)%can_temp
   theta1_1  = site%patch%soil_w(1)%theta(1)
   psi1_leaf = site%cohort%psi(1, 1)

   call check(nfail == 0_ik, 'whole-column budgets closed on every patch (n_fail == 0)')
   call check(we < 1.0e-3_wp, 'whole-column energy residual tiny')
   call check(ww < 1.0e-8_wp, 'whole-column water residual tiny')
   call check(abs(t_cas1 - t_cas0) > 0.05_wp, 'CAS temperature evolved under the fast loop')
   call check(site%patch%cas(2)%can_temp > 200.0_wp, 'bare (zero-cohort) patch fast step stayed physical')
   !----- The fast loop READ psi from the cohort block, evolved it, and WROTE it back (persist). !
   call check(psi1_leaf < psi0_leaf - 1.0e-4_wp, 'per-cohort leaf psi evolved + persisted on the cohort block')
   !----- Fast->slow carbon bridge: the vegetated cohort accumulated positive GROSS GPP. -------!
   call check(site%cohort%gpp_accum(1) > 0.0_wp, 'fast loop accumulated positive gross GPP (fast->slow bridge)')

   !=== 3. Stepper wiring + gate. =========================================================!
   block
      real(wp) :: t_before
      !----- Gate OFF (no fast context): reservoirs must NOT change. --------------------!
      call init_fast_reservoirs(site, ctx)
      t_before = site%patch%cas(1)%can_temp
      call advance_one_step(site, cfg, .false., .false.)          ! no fast_ctx -> fast loop skipped
      call check_close(site%patch%cas(1)%can_temp, t_before, 1.0e-12_wp, &
                       'fast loop must NOT run without a fast context')
      !----- Gate ON (context supplied): the hook fires, reservoirs change. -------------!
      call init_fast_reservoirs(site, ctx)
      t_before = site%patch%cas(1)%can_temp
      call advance_one_step(site, cfg, .false., .false., ctx)
      call check(abs(site%patch%cas(1)%can_temp - t_before) > 0.05_wp, &
                 'advance_one_step ran the fast loop when gated on')
   end block

   !=== 4. END-TO-END fast->slow handoff via the stepper. The SAME carbon-mode site is advanced  !
   !    with the fast loop OFF then ON, with gpp_ref = 0 so the stub adds nothing. Turnover is    !
   !    identical between the two, so more carbon under fast-ON is EXACTLY the fast GPP flowing    !
   !    through gpp_accum into carbon growth -- the fast->slow bridge, isolated from turnover.  ===!
   cfg%growth_source = GS_CARBON
   cfg%gpp_ref       = 0.0_wp

   cfg%fast_biophysics_on = .false.                    ! stub GPP = 0; carbon change is pure turnover
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(site)
   call advance_one_step(site, cfg, .false., .false., ctx)
   cbal0 = site%cohort%leaf_carbon(1) + site%cohort%fineroot_carbon(1)                          &
         + site%cohort%wood_carbon(1) + site%cohort%nonstructural_carbon(1)

   cfg%fast_biophysics_on = .true.                     ! fast loop -> gpp_accum -> carbon growth
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)
   call advance_one_step(site, cfg, .false., .false., ctx)
   cbal1 = site%cohort%leaf_carbon(1) + site%cohort%fineroot_carbon(1)                          &
         + site%cohort%wood_carbon(1) + site%cohort%nonstructural_carbon(1)

   call check(cbal1 > cbal0, 'fast GPP raised carbon vs fast-off (gpp_ref=0): the fast->slow bridge is live')

   write(*,'(a)')          '   PASS'
   write(*,'(a,f7.2,a,f7.2,a)') '   (patch-1 CAS temp ', t_cas0, ' -> ', t_cas1, ' K)'
   write(*,'(a,es10.3,a,es10.3,a)') '   (worst whole-column resid: energy=', we, ' J/m2  water=', ww, ' kg/m2)'
   write(*,'(a,es10.3,a)') '   (2 h accumulated gross GPP = ', site%cohort%gpp_accum(1), ' kgC/plant)'
end program test_fast_loop
