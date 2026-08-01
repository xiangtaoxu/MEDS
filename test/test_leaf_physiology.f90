!----- Leaf physiology: temperature response, FvCB C3 + Collatz C4 demand, the three stomatal !
!      models, the coupled A-gs-Ci solver, water stress, and the night/closed branch.         !
program test_leaf_physiology
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_kelvin
   use meds_config,             only : meds_config_t, SM_LEUNING, SM_MEDLYN, SM_KATUL,         &
                                       COLIM_MIN
   use meds_temp_response, only : arrhenius_scale, peaked_arrhenius_scale
   use meds_leaf_gas_exchange,only : assimilation_demand_c3
   use meds_leaf_gas_exchange,       only : stomata_gs_medlyn
   use meds_plant_types,         only : leaf_env_t, leaf_flux_t, LIM_NONE, LIM_RUBISCO,         &
                                       LIM_RUBP, LIM_C4_PEP
   use meds_plant_interface,    only : leaf_gas_exchange
   use meds_test_support,       only : check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(leaf_env_t)    :: env
   type(leaf_flux_t)   :: flux, flux2
   real(wp) :: a_gross, ac, aj, ap, an0, an1, an2, prev
   integer(ik) :: i
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
                        210000.0_wp, COLIM_MIN, 0.85_wp, a_gross, ac, aj, ap)
   call check_close(ac, 25.85930_wp, 1.0e-4_wp, 'C3 Ac known-answer')
   call check_close(aj, 20.52632_wp, 1.0e-4_wp, 'C3 Aj known-answer')
   call check_close(ap, 45.0_wp,     1.0e-6_wp, 'C3 Ap known-answer')
   call check_close(a_gross, 20.52632_wp, 1.0e-4_wp, 'C3 min co-limitation picks Aj')

   !=== 3. C3 limitation regimes (raw demand, COLIM_MIN). ===================================!
   !----- High light, low Ci -> Rubisco-limited. ---------------------------------------------!
   call assimilation_demand_c3(120.0_wp, 60.0_wp, 104.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,  &
                        209000.0_wp, COLIM_MIN, 0.85_wp, a_gross, ac, aj, ap)
   call check(ac < aj .and. ac < ap, 'low Ci, high light should be Rubisco-limited')
   !----- High Ci, low light (small J) -> RuBP/light-limited. ---------------------------------!
   call assimilation_demand_c3(300.0_wp, 60.0_wp, 20.0_wp, 10.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.85_wp, a_gross, ac, aj, ap)
   call check(aj < ac .and. aj < ap, 'high Ci, low light should be RuBP-limited')
   !----- High Ci, high light, low TPU -> product-limited. ------------------------------------!
   call assimilation_demand_c3(600.0_wp, 60.0_wp, 104.0_wp, 2.0_wp, 42.0_wp, 400.0_wp, 275000.0_wp,   &
                        209000.0_wp, COLIM_MIN, 0.85_wp, a_gross, ac, aj, ap)
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
                    cfg%pft%stomatal_g0(1), cfg%pft%stomatal_g1(1)), 1.0e-4_wp,                 &
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

   !=== 6b. C4 kp temperature response (BUG4): at a CO2/PEP-limiting Ci the PEP-limited gross    !
   !     rate must RISE with leaf temperature -- kp was previously FROZEN at its 25 degC value. ==!
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
   !     (was gs alone, overestimating the water flux). ==========================================!
   block
      real(wp) :: e_series
      cfg%leaf_use_boundary_layer = .true.
      env = std_env() ; env%gb = 0.6_wp                      ! finite water-vapour boundary conductance [mol/m2/s]
      call leaf_gas_exchange(env, cfg, 1_ik, flux)
      e_series = flux%gs * env%gb / (flux%gs + env%gb) * env%vpd / env%pressure
      call check_close(flux%transpiration, e_series, 1.0e-9_wp, 'transpiration uses the gs-gb series conductance')
      call check(flux%transpiration < flux%gs * env%vpd / env%pressure,                            &
                 'series transpiration is below the gs-only value (gb resistance applied)')
      cfg%leaf_use_boundary_layer = .false.
   end block

   !=== 7-. The capacity limb is OFF by default (issue #47), so psi_leaf must be INERT. ======!
   !     Guards the default: if someone re-enables it silently, this fires before the sweeps below. !
   env = std_env() ; env%psi_leaf =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an0 = flux%A_net
   env%psi_leaf = -5.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an1 = flux%A_net
   call check_close(an1, an0, 1.0e-12_wp, 'psi_leaf is inert when wstress_nonstomatal is off (default)')

   !=== 7. Water stress (PFT 1, Medlyn): A_net falls monotonically as psi_leaf drops. =======!
   !     The capacity limb still EXISTS and is still tested -- it is opt-in, not deleted, so every  !
   !     assertion below enables it explicitly. ====================================================!
   cfg%leaf_wstress_nonstomatal = .true.
   env = std_env() ; env%psi_leaf =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an0 = flux%A_net
   env%psi_leaf = -1.5_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an1 = flux%A_net
   env%psi_leaf = -5.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux) ; an2 = flux%A_net
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

   !=== 7d. Stomatal-limb stress (Sabot beta_stomata via psi_soil): with psi_leaf = 0 the capacity !
   !     limb is OFF, so a drop in SOIL water potential must close stomata (gs falls) through g1. ==!
   env = std_env() ; env%psi_leaf = 0.0_wp
   env%psi_soil =  0.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an0 = flux%gs
   env%psi_soil = -1.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux)  ; an1 = flux%gs
   env%psi_soil = -3.0_wp ; call leaf_gas_exchange(env, cfg, 1_ik, flux2) ; an2 = flux2%gs
   call check(flux2%converged, 'stomatal-limb (psi_soil) solve must converge')
   call check(an0 > an1 .and. an1 > an2, 'gs must fall as soil water potential drops (beta_stomata on g1)')

   cfg%leaf_wstress_nonstomatal = .false.     ! back to the shipped default for the rest of the suite

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

   write(*,'(a)') '   PASS'

contains

   !----- A representative midday tropical leaf environment. ------------------------------!
   function std_env() result(e)
      type(leaf_env_t) :: e
      e%par = 1500.0_wp ; e%leaf_temp = t_kelvin + 25.0_wp ; e%vpd = 1500.0_wp
      e%ca = 400.0_wp ; e%pressure = 101325.0_wp ; e%psi_leaf = 0.0_wp ; e%gb = 0.0_wp ; e%psi_soil = 0.0_wp
   end function std_env

   !----- build_test_config lives in meds_test_support; wrap it for clarity. --------------!
   function build_test_config_local() result(c)
      use meds_test_support, only : build_test_config
      type(meds_config_t) :: c
      c = build_test_config()
   end function build_test_config_local

end program test_leaf_physiology
