! SPDX-License-Identifier: Apache-2.0
!----- Leaf physiology: temperature response, FvCB C3 + Collatz C4 demand, the three stomatal !
!      models, the coupled A-gs-Ci solver, water stress, and the night/closed branch          !
!      (docs/science/leaf_gas_exchange.md). ----------------------------------------------------!
program test_leaf_physiology
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_kelvin
   use meds_config,             only : meds_config_t, pft_leaf_psi_tlp, pft_stomata_psi_onset,  &
                                       pft_wstress_psi_open, pft_wstress_psi_close
   use meds_pft_params,         only : HYD_UNSET
   use meds_leaf_params,        only : SM_LEUNING, SM_MEDLYN, SM_KATUL, COLIM_MIN, COLIM_QUADRATIC
   use meds_therm_lib,     only : arrhenius_scale, peaked_arrhenius_scale, kattge_knorr_entropy,  &
                                  kattge_knorr_jv_ratio
   use meds_leaf_gas_exchange,only : assimilation_demand_c3, rubisco_kinetics
   use meds_leaf_gas_exchange,       only : stomata_gs_medlyn
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_args_t, LIM_NONE, LIM_RUBISCO, LIM_RUBP,   &
                                LIM_C4_PEP
   use meds_fast_config, only : leaf_gas_exchange, leaf_photo_params_for_pft
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t) :: cfg
   type(leaf_env_t)    :: env
   type(leaf_flux_t)   :: flux, flux2
   real(wp) :: a_gross, ac, aj, ap, an0, an1, an2, prev
   real(wp) :: psi_tlp, jump, gs_top
   integer(ik) :: i, j
   integer(ik), dimension(3) :: sms = [ SM_LEUNING, SM_MEDLYN, SM_KATUL ]

   call banner('leaf physiology (photosynthesis + stomata)')
   cfg = build_test_config_local()

   !=== 1. Temperature response: identity at the reference, monotonicity, peaked optimum. ===!
   call check_close(arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 25.0_wp), 40.49_wp,        &
                    1.0e-12_wp, 'Arrhenius must equal k25 at the reference temperature')
   call check(arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 5.0_wp) <                        &
              arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 35.0_wp),                        &
              'Arrhenius must increase with temperature')
   call check_close(peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, 650.0_wp,          &
                    t_kelvin + 25.0_wp), 60.0_wp, 1.0e-12_wp, 'peaked must equal k25 at reference')
   !----- Peaked form: a single interior maximum near ~32 degC (rise then fall). ------------!
   an0 = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, 650.0_wp, t_kelvin + 25.0_wp)
   an1 = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, 650.0_wp, t_kelvin + 32.0_wp)
   an2 = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, 650.0_wp, t_kelvin + 45.0_wp)
   call check(an1 > an0 .and. an1 > an2, 'peaked Vcmax must have an interior thermal optimum')

   !=== 2. C3 demand: hand-computed Ac/Aj known-answer (vcmax=100, j=120, ci=300, gstar=40, ==!
   !       kc=400, ko=275000, o2=210000, tpu=15). Ac=100*260/(300+400*(1+210000/275000))      !
   !       = 25.8593; Aj=120*260/(4*300+8*40)=20.5263; Ap=3*15=45. =============================!
   call assimilation_demand_c3(300.0_wp, 100.0_wp, 120.0_wp, 15.0_wp, 40.0_wp, 400.0_wp, 275000.0_wp, &
                        210000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check_close(ac, 25.85930_wp, 1.0e-4_wp, 'C3 Ac known-answer')
   call check_close(aj, 20.52632_wp, 1.0e-4_wp, 'C3 Aj known-answer')
   call check_close(ap, 45.0_wp,     1.0e-6_wp, 'C3 Ap known-answer')
   call check_close(a_gross, 20.52632_wp, 1.0e-4_wp, 'C3 min co-limitation picks Aj')

   !=== 1b. THERMAL ACCLIMATION, Kattge & Knorr (2007): docs/science/leaf_gas_exchange.md eq. (4). =!
   !        dS = a - b*T_growth[degC] is an arithmetic identity, so what these assert is that
   !        the CODE computes the published relation -- and then that the shift does what it is
   !        for: moving where the response peaks.
   block
      real(wp) :: ds_cold, ds_warm, jv_cold, jv_warm, k_cold, k_warm, topt_cold, topt_warm
      real(wp), parameter :: T_COLD = t_kelvin + 10.0_wp, T_WARM = t_kelvin + 30.0_wp
      print '(a)', 'test_thermal_acclimation:'
      !----- Hand-computed: dS_vcmax = 668.39 - 1.07*T. At 10 C -> 657.69; at 30 C -> 636.29. -!
      ds_cold = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_COLD)
      ds_warm = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_WARM)
      call check_close(ds_cold, 657.69_wp, 1.0e-10_wp, 'dS_vcmax at a 10 C growth temperature')
      call check_close(ds_warm, 636.29_wp, 1.0e-10_wp, 'dS_vcmax at a 30 C growth temperature')
      !----- Jmax:Vcmax = 2.59 - 0.035*T. At 10 C -> 2.24; at 30 C -> 1.54. -------------------!
      jv_cold = kattge_knorr_jv_ratio(2.59_wp, 0.035_wp, T_COLD)
      jv_warm = kattge_knorr_jv_ratio(2.59_wp, 0.035_wp, T_WARM)
      call check_close(jv_cold, 2.24_wp, 1.0e-10_wp, 'Jmax:Vcmax at a 10 C growth temperature')
      call check_close(jv_warm, 1.54_wp, 1.0e-10_wp, 'Jmax:Vcmax at a 30 C growth temperature')
      call check(jv_warm < jv_cold, 'a warm-grown plant invests relatively less in electron transport')
      !----- The POINT of the dS shift: at a hot leaf temperature the warm-acclimated plant must !
      !      retain more capacity than the cold-acclimated one.  --------------------------------!
      k_cold = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, ds_cold, t_kelvin + 38.0_wp)
      k_warm = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, ds_warm, t_kelvin + 38.0_wp)
      call check(k_warm > k_cold, 'warm acclimation retains more Vcmax at 38 C than cold acclimation')
      topt_cold = peak_temperature(ds_cold)
      topt_warm = peak_temperature(ds_warm)
      call check(topt_warm > topt_cold + 1.0_wp,                                                 &
                 'the thermal optimum moves UP with growth temperature (> 1 K over a 20 K range)')
      print '(a,2(f7.2,a))', '   [acclimation] Vcmax optimum: cold-grown ', topt_cold - t_kelvin,       &
            ' C, warm-grown ', topt_warm - t_kelvin, ' C'
   end block

   !=== 2b. CO-LIMITATION strength. The two C3 smoothings must use the CO-LIMITATION ========!
   !        curvatures, not theta_j (docs/science/leaf_gas_exchange.md §2). This is a REGIME     !
   !        assertion: at ambient CO2 Ac and Aj sit close together, so the smoothing penalty     !
   !        lands exactly where the model spends most of its time, and borrowing theta_j = 0.85  !
   !        costs ~29 % of assimilation against the sharp min(Ac,Aj,Ap). Measured here at PFT-1  !
   !        kinetics, Ci = 280, saturating light.                                                !
   block
      real(wp) :: a_min, a_new, a_old
      real(wp), parameter :: CI = 280.0_wp, VC = 60.0_wp, JR = 108.0_wp, TP = 15.0_wp
      real(wp), parameter :: GS = 42.75_wp, KC = 404.9_wp, KO = 278400.0_wp, O2 = 209000.0_wp
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_MIN,                     &
                           0.98_wp, 0.95_wp, a_min, ac, aj, ap)
      !----- The shipped curvatures. --------------------------------------------------------!
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_QUADRATIC,               &
                           0.98_wp, 0.95_wp, a_new, ac, aj, ap)
      !----- What borrowing theta_j would do: the size of the error it makes. ----------------!
      call assimilation_demand_c3(CI, VC, JR, TP, GS, KC, KO, O2, COLIM_QUADRATIC,               &
                           0.85_wp, 0.85_wp, a_old, ac, aj, ap)
      call check(a_old < 0.80_wp * a_min,                                                        &
                 'theta_j = 0.85 as a co-limitation curvature must cost > 20 % vs min()')
      call check(a_new > 0.88_wp * a_min,                                                        &
                 'co-limitation curvatures must keep the smoothed rate within 12 % of min()')
      call check(a_new < a_min,                                                                  &
                 'smoothing must still sit below the sharp min(), not above it')
      print '(a,3(f8.4,a))', '   [co-limitation] min(Ac,Aj,Ap) = ', a_min, ' | theta_*_c3 = ', a_new,     &
            ' | theta_j borrowed = ', a_old, ''
   end block

   !=== 3. C3 limitation regimes (raw demand, COLIM_MIN). ===================================!
   !----- High light, low Ci -> Rubisco-limited. ---------------------------------------------!
   call assimilation_demand_c3(120.0_wp, 60.0_wp, 104.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,  &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(ac < aj .and. ac < ap, 'low Ci, high light should be Rubisco-limited')
   !----- High Ci, low light (small J) -> RuBP/light-limited. ---------------------------------!
   call assimilation_demand_c3(300.0_wp, 60.0_wp, 20.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(aj < ac .and. aj < ap, 'high Ci, low light should be RuBP-limited')
   !----- High Ci, high light, low TPU -> product-limited. ------------------------------------!
   call assimilation_demand_c3(600.0_wp, 60.0_wp, 104.0_wp, 2.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.98_wp, 0.95_wp, a_gross, ac, aj, ap)
   call check(ap < ac .and. ap < aj, 'high Ci, low TPU should be product-limited')

   !=== 4. Full solve (PFT 1, C3, Medlyn default): convergence, bounds, diffusion closure. ==!
   env = std_env()
   call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check(flux%converged, 'C3 solve must converge')
   call check(flux%A_net > 0.0_wp, 'C3 net assimilation must be positive at midday')
   call check(flux%gs >= cfg%pft%stomatal_g0(1), 'gs must be at least the cuticular g0')
   call check(flux%ci > 0.0_wp .and. flux%ci < flux%cs, 'Ci must lie between 0 and Cs')
   call check(flux%limitation == LIM_RUBISCO .or. flux%limitation == LIM_RUBP,                 &
              'midday C3 leaf should be Rubisco- or RuBP-limited')
   !----- Diffusion identity A_net = gs*(Cs-Ci)/1.6 (back-computation consistency). ----------!
   call check_close(flux%A_net, flux%gs * (flux%cs - flux%ci) / 1.6_wp, 1.0e-6_wp,             &
                    'A_net must satisfy the CO2 diffusion identity')
   !----- The converged gs reproduces the Medlyn law (root really solved the coupled system). -!
   call check_close(flux%gs, stomata_gs_medlyn(flux%A_net, flux%cs, env%vpd,                   &
                    cfg%pft%stomatal_g0(1), cfg%pft%stomatal_g1(1), cfg%medlyn_vpd_min), 1.0e-4_wp, &
                    'back-computed gs must match the Medlyn model at the solution')

   !=== 5. All three stomatal models converge; gs decreases as VPD rises. ===================!
   do i = 1_ik, 3_ik
      cfg%stomatal_model = sms(i)
      env = std_env() ; env%vpd = 800.0_wp
      call leaf_gas_exchange(env, cfg, 1_ik, flux)
      env%vpd = 2500.0_wp
      call leaf_gas_exchange(env, cfg, 1_ik, flux2)
      call check(flux%converged .and. flux2%converged, 'each stomatal model must converge')
      call check(flux%gs >= cfg%pft%stomatal_g0(1) .and. flux%A_net > 0.0_wp,                  &
                 'each stomatal model must give a positive, open-stomata solution')
      call check(flux2%gs < flux%gs, 'gs must decrease as VPD increases')
   end do
   cfg%stomatal_model = SM_MEDLYN

   !=== 6. C3 vs C4 contrast (PFT 1 vs PFT 3): C4 runs and draws Ci lower (concentrating). ==!
   env = std_env()
   call leaf_gas_exchange(env, cfg, 1_ik, flux)    ! C3
   call leaf_gas_exchange(env, cfg, 3_ik, flux2)   ! C4
   call check(flux2%converged .and. flux2%A_net > 0.0_wp, 'C4 leaf must converge and assimilate')
   call check(flux2%limitation == LIM_RUBISCO .or. flux2%limitation == LIM_RUBP .or.           &
              flux2%limitation == LIM_C4_PEP, 'C4 limitation flag must be a C4-valid regime')
   call check(flux2%ci / flux2%cs < flux%ci / flux%cs, 'C4 should operate at a lower Ci/Cs than C3')

   !=== 6b. C4 kp temperature response: at a CO2/PEP-limiting Ci the PEP-limited gross rate     !
   !     must RISE with leaf temperature, which it cannot with kp held at its 25 degC value. =====!
   block
      real(wp) :: a_cold, a_warm
      env = std_env() ; env%ca = 40.0_wp                     ! very low CO2 -> C4 PEP-limited
      env%leaf_temp = t_kelvin + 18.0_wp
      call leaf_gas_exchange(env, cfg, 3_ik, flux)           ! C4 (PFT 3), cold
      a_cold = flux%A_gross
      env%leaf_temp = t_kelvin + 34.0_wp
      call leaf_gas_exchange(env, cfg, 3_ik, flux2)          ! C4 (PFT 3), warm
      a_warm = flux2%A_gross
      call check(flux%converged .and. flux2%converged, 'C4 low-CO2 solve must converge at both temperatures')
      call check(flux%limitation == LIM_C4_PEP, 'low-CO2 C4 leaf should be PEP-limited (isolates kp)')
      call check(a_warm > a_cold, 'C4 PEP-limited gross assimilation must rise with temperature (kp now temp-scaled)')
   end block

   !=== 6c. Transpiration puts stomata gs and the boundary layer gb in SERIES when use_boundary_layer is on   !
   !     (gs alone would overestimate the water flux). ===========================================!
   block
      real(wp) :: e_series
      cfg%use_boundary_layer = .true.
      env = std_env() ; env%gb = 0.6_wp                      ! finite water-vapour boundary conductance [mol/m2/s]
      call leaf_gas_exchange(env, cfg, 1_ik, flux)
      e_series = flux%gs * env%gb / (flux%gs + env%gb) * env%vpd / env%pressure
      call check_close(flux%transpiration, e_series, 1.0e-9_wp, 'transpiration uses the gs-gb series conductance')
      call check(flux%transpiration < flux%gs * env%vpd / env%pressure,                            &
                 'series transpiration is below the gs-only value (gb resistance applied)')
      cfg%use_boundary_layer = .false.
   end block

   !=== 7-. The fast limb is idle while the leaf holds turgor: psi_leaf above the turgor-loss   !
   !     point is inert. ====================================================================!
   psi_tlp = pft_leaf_psi_tlp(cfg, 1_ik)
   env = std_env() ; env%psi_leaf =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an0 = flux%A_net
   env%psi_leaf = 0.9_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an1 = flux%A_net
   call check_close(an1, an0, 1.0e-12_wp, 'psi_leaf above the turgor-loss point is inert')

   !=== 7. Water stress (PFT 1, Medlyn): A_net falls monotonically as the leaf wilts past its      !
   !     turgor-loss point, to -Rd once the fast limb has shut the biochemistry. ===================!
   env = std_env() ; env%psi_leaf =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an0 = flux%A_net
   env%psi_leaf = 1.5_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an1 = flux%A_net
   env%psi_leaf = 2.5_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an2 = flux%A_net
   call check(an0 > an1 .and. an1 > an2, 'A_net must decrease as the leaf dries')
   call check_close(an2, -flux%rd, 1.0e-6_wp, 'fully stressed C3 leaf nets -Rd (beta_nonstomata = 0)')

   !=== 7b. Katul under full water stress must CLOSE (g0 fallback), not return an open flux. =!
   cfg%stomatal_model = SM_KATUL
   env = std_env() ; env%psi_leaf =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an0 = flux%A_net
   env%psi_leaf = -5.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux2)
   call check(flux2%converged, 'Katul full-stress solve must converge via the g0 fallback')
   call check(flux2%A_net < an0, 'Katul full water stress must reduce assimilation')
   call check_close(flux2%gs, cfg%pft%stomatal_g0(1), 1.0e-2_wp, 'Katul full-stress leaf must close to ~g0')

   !=== 7c. Every OPEN Katul solve must be diffusion-consistent: gs*(cs-ci)/1.6 == A_net across a  !
   !     water-stress sweep. The g0-pinned re-solve keeps A/gs/Ci mutually consistent even when the  !
   !     Katul optimum would otherwise fall below the cuticular floor g0. ==========================!
   do i = 0_ik, 8_ik
      env = std_env() ; env%psi_leaf = -0.5_wp * real(i, wp)                 ! 0 .. -4 MPa
      call leaf_gas_exchange(env, cfg, 1_ik, flux)
      call check(flux%converged, 'Katul stress sweep must converge')
      if (flux%gs > cfg%pft%stomatal_g0(1) * (1.0_wp + 1.0e-6_wp)) then      ! open stomata only
         call check_close(flux%gs * (flux%cs - flux%ci) / 1.6_wp, flux%A_net, 1.0e-6_wp,          &
                          'Katul open solve is diffusion-consistent (gs, cs, ci, A match)')
      end if
   end do
   cfg%stomatal_model = SM_MEDLYN

   !=== 7d. Stomatal-limb stress (Sabot beta_stomata via psi): with psi_leaf = 0 the fast limb  !
   !     is idle, so a drop in the predawn potential must close stomata (gs falls) through g1. ===!
   env = std_env() ; env%psi_leaf = 0.0_wp
   env%psi =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an0 = flux%gs
   env%psi = -1.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an1 = flux%gs
   env%psi = -3.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux2) ; an2 = flux2%gs
   call check(flux2%converged, 'stomatal-limb (psi) solve must converge')
   call check(an0 > an1 .and. an1 > an2, 'gs must fall as soil water potential drops (beta_stomata on g1)')

   !=== 7e. The fast limb (beta_nonstomata) on the leaf's current potential: 1 down to the     !
   !     turgor-loss point, falling linearly to 0 at twice it (the defaults; the PFT file may set  !
   !     wstress_psi_open/_close). Every stomatal model ties gs to A, so as the leaf wilts gs falls !
   !     continuously to the cuticular g0, and a wilted leaf nets -Rd by day. ======================!
   call check_close(pft_wstress_psi_open(cfg, 1_ik),  psi_tlp,          1.0e-14_wp,                 &
                    'fast limb starts at the turgor-loss point')
   call check_close(pft_wstress_psi_close(cfg, 1_ik), 2.0_wp * psi_tlp, 1.0e-14_wp,                 &
                    'fast limb completes at twice the turgor-loss point')
   env = std_env() ; env%psi_leaf = 0.9_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 1.0_wp, 0.0_wp, 'fast limb: 1 while the leaf holds turgor')
   env%psi_leaf = 1.5_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 0.5_wp, 1.0e-12_wp, 'fast limb: 0.5 half way')
   env%psi_leaf = 2.0_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check_close(flux%beta_nonstomata, 0.0_wp, 0.0_wp, 'fast limb: 0 at twice the turgor-loss point')
   do j = 1_ik, size(sms)
      cfg%stomatal_model = sms(j)
      env = std_env()
      env%psi_leaf = psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; gs_top = flux%gs
      jump = 0.0_wp ; prev = huge(1.0_wp)
      do i = 0_ik, 280_ik
         env%psi_leaf = psi_tlp * (0.8_wp + 1.4_wp * real(i, wp) / 280.0_wp)
         call leaf_gas_exchange(env, cfg, 1_ik, flux)
         call check(flux%converged .and. flux%gs >= cfg%pft%stomatal_g0(1) * (1.0_wp - 1.0e-12_wp),   &
                    'fast-limb sweep: converged, gs >= g0')
         !----- to the Ci solve's tolerance (1e-3 ppm): near the g0 floor Cs - Ci is ~30 ppm, so the  !
         !      back-computed gs carries ~1e-5 relative noise (Katul switches to its g0 pass there). !
         if (i > 0_ik) call check(flux%gs <= prev * (1.0_wp + 1.0e-4_wp) + 1.0e-12_wp,            &
                                  'fast-limb sweep: gs never rises as the leaf dries')
         if (i > 0_ik) jump = max(jump, prev - flux%gs)
         if (flux%gs > cfg%pft%stomatal_g0(1) * (1.0_wp + 1.0e-6_wp)) then
            call check_close(flux%gs * (flux%cs - flux%ci) / 1.6_wp, flux%A_net, 1.0e-6_wp,           &
                             'fast-limb sweep: A, gs, Ci and Cs stay diffusion-consistent')
         end if
         prev = flux%gs
      end do
      call check(jump < 0.05_wp * gs_top, 'fast-limb sweep: no step larger than 5% of gs at the turgor-loss point')
      env%psi_leaf = 2.1_wp * psi_tlp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
      call check_close(flux%gs, cfg%pft%stomatal_g0(1), 1.0e-9_wp, 'wilted leaf by day: gs is the cuticular g0')
      call check_close(flux%A_gross, 0.0_wp, 1.0e-12_wp, 'wilted leaf by day: no photosynthesis')
      call check_close(flux%A_net, -flux%rd, 1.0e-12_wp, 'wilted leaf by day: the leaf respires')
   end do
   cfg%stomatal_model = SM_MEDLYN

   !=== 7f. The stomatal stress begins at an onset: half the turgor-loss point unless the PFT     !
   !     file sets stomata_psi_onset. Above it beta_stomata is 1, so the gravity head a tall tree   !
   !     has in wet soil does not close its stomata; below it beta = exp(sref*(psi - onset)). =====!
   psi_tlp = pft_leaf_psi_tlp(cfg, 1_ik)
   call check_close(pft_stomata_psi_onset(cfg, 1_ik), 0.5_wp * psi_tlp, 1.0e-14_wp,                 &
                    'onset defaults to half the turgor-loss point')
   env = std_env() ; env%psi_leaf = 0.0_wp
   env%psi = 0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; gs_top = flux%gs
   env%psi = 0.99_wp * pft_stomata_psi_onset(cfg, 1_ik) ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check(flux%beta_stomata == 1.0_wp .and. flux%gs == gs_top, 'onset: no stress just above it')
   env%psi = pft_stomata_psi_onset(cfg, 1_ik) - 0.5_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check_close(flux%beta_stomata, exp(-0.5_wp * cfg%pft%wstress_sref_stomata(1)), 1.0e-14_wp,     &
                    'onset: beta = exp(sref*(psi - onset)) below it')
   cfg%pft%stomata_psi_onset(1) = -0.3_wp
   call check_close(pft_stomata_psi_onset(cfg, 1_ik), -0.3_wp, 0.0_wp, 'onset: a set value is used')
   env%psi = -0.5_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check_close(flux%beta_stomata, exp(-0.2_wp * cfg%pft%wstress_sref_stomata(1)), 1.0e-14_wp,     &
                    'onset: a set value moves the threshold')
   cfg%pft%stomata_psi_onset(1) = HYD_UNSET

   !=== 8. PAR sweep: night/closed branch at PAR=0, monotone rise, no NaNs, all converge. ===!
   env = std_env() ; env%par = 0.0_wp
   call leaf_gas_exchange(env, cfg, 1_ik, flux)
   call check(flux%limitation == LIM_NONE, 'PAR=0 must hit the closed/night branch')
   call check_close(flux%A_net, -flux%rd, 1.0e-6_wp, 'night A_net must equal -Rd')
   call check_close(flux%gs, cfg%pft%stomatal_g0(1), 1.0e-12_wp, 'night gs must be g0')
   prev = -1.0e30_wp
   do i = 0_ik, 20_ik
      env%par = 100.0_wp * real(i, wp)
      call leaf_gas_exchange(env, cfg, 1_ik, flux)
      call check(flux%converged, 'every PAR level must converge')
      call check(flux%A_net == flux%A_net, 'A_net must never be NaN')   ! NaN /= itself
      call check(flux%A_net >= prev - 1.0e-9_wp, 'A_net must rise monotonically with PAR')
      prev = flux%A_net
   end do

   !=== 9. O2 propagates to BOTH places oxygen enters the C3 demand. =======================!
   !      Gamma* = 0.5*O/S_(c/o) is proportional to the O2 partial pressure
   !      (docs/science/leaf_gas_exchange.md eq. (3)). If o2_mol_frac reached only the Michaelis
   !      term Kc(1+O/Ko), raising O2 would inhibit carboxylation while leaving the
   !      photorespiratory penalty on Aj untouched.
   env = std_env()
   cfg%o2_mol_frac = 0.209_wp
   call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an0 = flux%A_net       ! the calibration point
   cfg%o2_mol_frac = 0.105_wp
   call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an1 = flux%A_net       ! half O2
   cfg%o2_mol_frac = 0.350_wp
   call leaf_gas_exchange(env, cfg, 1_ik, flux2) ; an2 = flux2%A_net      ! paleo-high O2
   call check(flux%converged .and. flux2%converged, 'the O2 sweep must converge')
   !----- Less O2 => less photorespiration => MORE assimilation, and vice versa. Monotone in O2. -!
   call check(an1 > an0, 'halving O2 must RAISE A_net (less photorespiration)')
   call check(an2 < an0, 'raising O2 to 35% must LOWER A_net (more photorespiration)')
   !----- The response must run through Gamma*, NOT just the Michaelis term. This is the assertion -!
   !      a frozen Gamma* fails, and the thresholds are set from measuring BOTH on this exact       !
   !      fixture rather than guessed:                                                              !
   !                                                                                          !
   !                        O2 = 10.5%      O2 = 35%                                                !
   !        with Gamma*(O2)   +22.0 %       -23.7 %                                                  !
   !        frozen Gamma*      +8.2 %       -10.5 %   <- Kc(1+O/Ko) acting alone                      !
   !                                                                                          !
   !      A >5% gate would pass on BOTH; 15% separates them with room on either side.                 !
   call check((an1 - an0) / abs(an0) > 0.15_wp,                                                  &
              'halving O2 must raise A_net > 15% -- the response runs through Gamma*, not just Kc')
   call check((an2 - an0) / abs(an0) < -0.15_wp,                                                 &
              'raising O2 to 35% must cut A_net > 15% -- likewise')
   cfg%o2_mol_frac = 0.209_wp                 ! back to the shipped default

   !=== 10. Gamma* as a mole fraction does not depend on pressure (docs/science/leaf_gas_exchange.md!
   !       eq. (3)): it is proportional to the O2 partial pressure, which scales with pressure as   !
   !       the CO2 partial pressure does. At 25 C and 21 % O2 it is Bernacchi et al. (2001)'s       !
   !       42.75 umol/mol at sea level and at 3000 m (70 kPa) alike. Kc and Ko are partial          !
   !       pressures, so as mole fractions they scale as 1/P. ======================================!
   block
      type(leaf_photo_args_t)   :: pp
      real(wp) :: kc_low, ko_low, gstar_low, o2_low, kc_high, ko_high, gstar_high, o2_high
      call leaf_photo_params_for_pft(cfg, 1_ik, pp)
      call rubisco_kinetics(pp, t_kelvin + 25.0_wp, 101325.0_wp, kc_low, ko_low, gstar_low, o2_low)
      call rubisco_kinetics(pp, t_kelvin + 25.0_wp, 70000.0_wp, kc_high, ko_high, gstar_high, o2_high)
      call check(gstar_high == gstar_low, 'Gamma* in umol/mol does not depend on pressure')
      call check_close(gstar_low, 42.75_wp, 1.0e-12_wp, 'Gamma* at 25 C and 21 % O2 is 42.75 umol/mol')
      call check(o2_high == o2_low, 'O2 in umol/mol does not depend on pressure')
      call check_close(kc_high * 70000.0_wp, kc_low * 101325.0_wp, 1.0e-12_wp, 'Kc is a partial pressure')
      call check_close(ko_high * 70000.0_wp, ko_low * 101325.0_wp, 1.0e-12_wp, 'Ko is a partial pressure')
      !----- The same Gamma* at 35 C, where the Arrhenius factor is not 1. -----------------------!
      call rubisco_kinetics(pp, t_kelvin + 35.0_wp, 101325.0_wp, kc_low, ko_low, gstar_low, o2_low)
      call rubisco_kinetics(pp, t_kelvin + 35.0_wp, 70000.0_wp, kc_high, ko_high, gstar_high, o2_high)
      call check(gstar_high == gstar_low, 'Gamma* in umol/mol does not depend on pressure at 35 C')
   end block

   write(*,'(a)') '   PASS'

contains

   !----- Locate the peak of the peaked-Arrhenius response by a coarse scan, so the optimum is   !
   !      MEASURED from the function rather than re-derived from a formula the test would then    !
   !      be checking against itself. ------------------------------------------------------------!
   pure real(wp) function peak_temperature(ds) result(t_opt)
      real(wp), intent(in) :: ds
      real(wp)    :: t, k, k_best
      integer(ik) :: j
      k_best = -1.0_wp ; t_opt = 0.0_wp
      do j = 0_ik, 600_ik
         t = t_kelvin + 5.0_wp + 0.05_wp * real(j, wp)
         k = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, ds, t)
         if (k > k_best) then ; k_best = k ; t_opt = t ; end if
      end do
   end function peak_temperature

   !----- A representative midday tropical leaf environment. ------------------------------!
   function std_env() result(e)
      type(leaf_env_t) :: e
      e%par = 1500.0_wp ; e%leaf_temp = t_kelvin + 25.0_wp ; e%vpd = 1500.0_wp
      e%ca = 400.0_wp ; e%pressure = 101325.0_wp ; e%psi_leaf = 0.0_wp ; e%gb = 0.0_wp ; e%psi = 0.0_wp
   end function std_env

   !----- build_test_config lives in meds_test_support; wrap it for clarity. --------------!
   function build_test_config_local() result(c)
      use meds_test_support, only : build_test_config
      type(meds_config_t) :: c
      c = build_test_config()
   end function build_test_config_local

end program test_leaf_physiology
