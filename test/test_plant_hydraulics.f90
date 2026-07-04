!==========================================================================================!
! test_plant_hydraulics -- unit tests for the stateless plant-hydraulics library.            !
!                                                                                          !
!   1. PV CURVE      : turgor-loss point matches Bartlett eqn 1; psi<->RWC round-trips;        !
!                      capacitance C matches a finite-difference dW/dpsi.                        !
!   2. KIRCHHOFF EDGE : the integrated conductance uses <plc> in (plc(leaf), plc(wood)); the      !
!                       |dpsi|->0 limit reduces to the pointwise value.                            !
!   3. NO-FLOW EQUIL. : E=0 at the hydrostatic equilibrium (psi_L = psi_W - G, psi_W = psi_soil)    !
!                       gives dpsi ~ 0 and zero boundary fluxes (conservation).                      !
!   4. STEADY / OHM   : under constant E the wood node converges to psi_soil - E/rhizo_cond and      !
!                       root_uptake -> E (mass balance); leaf is more tensioned than wood.             !
!   5. vs REFERENCE   : the adaptive matrix-exponential matches a fine forward-Euler integration       !
!                       of the same ODE over one step.                                                  !
!   6. DIURNAL        : a sinusoidal E gives a midday psi_leaf minimum and nocturnal recovery, bounded.  !
!   7. DEGENERATE     : a near-leafless cohort returns finite values (no NaN/Inf).                        !
!==========================================================================================!
program test_plant_hydraulics
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : grav_head
   use meds_hydro_pv,          only : pv_psi_tlp, rwc_from_psi, psi_from_rwc, water_content,   &
                                      capacitance
   use meds_hydro_conductance, only : plc_retained, flux_potential, kirchhoff_edge
   use meds_plant_hydraulics,  only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t, &
                                      plant_water_flux, N_HYDRO, NODE_LEAF, NODE_WOOD,          &
                                      HYDRO_SUBSTEP_FIXED
   implicit none

   integer(ik) :: nfail
   nfail = 0_ik

   call test_pv_curve()
   call test_kirchhoff_edge()
   call test_no_flow_equilibrium()
   call test_steady_ohm()
   call test_vs_reference()
   call test_diurnal()
   call test_degenerate()

   if (nfail == 0_ik) then
      print '(a)', 'test_plant_hydraulics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_plant_hydraulics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- Assertion. -----------------------------------------------------------------------!
   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5,a,es10.2)', '  FAIL : ', name, '  got ', got,          &
               ' expected ', expect, '  |diff|>', atol
      end if
   end subroutine check

   subroutine check_true(name, cond)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a)', '  FAIL : ', name
      end if
   end subroutine check_true

   !----- A pioneer-like default trait set + cohort environment. ---------------------------!
   subroutine defaults(p, env, o)
      type(hydro_params_t), intent(out) :: p
      type(hydro_env_t),    intent(out) :: env
      type(hydro_opts_t),   intent(out) :: o
      p%leaf_pi0 = -1.5_wp ; p%leaf_eps = 12.0_wp ; p%leaf_af = 0.30_wp ; p%leaf_water_sat = 2.0_wp
      p%wood_pi0 = -1.0_wp ; p%wood_eps =  8.0_wp ; p%wood_af = 0.20_wp ; p%wood_water_sat = 1.0_wp
      p%wood_psi50 = -2.0_wp ; p%wood_kexp = 2.0_wp
      p%k_plant_max = 6.0e-4_wp ; p%wood_kmax = 8.0_wp ; p%vessel_curl = 1.5_wp
      env%transp = 1.5e-4_wp ; env%soil_psi = -0.3_wp ; env%rhizo_cond = 5.0e-4_wp
      env%bleaf = 0.5_wp ; env%bsap = 5.0_wp ; env%broot = 2.0_wp
      env%sap_area = 0.01_wp ; env%height = 20.0_wp ; env%leaf_area = 5.0_wp
      o = hydro_opts_t()   ! defaults: 2-node, EXPM, KPLANT, adaptive, gravity on
   end subroutine defaults

   !=======================================================================================!
   subroutine test_pv_curve()
      real(wp), parameter :: pi0 = -1.5_wp, eps = 12.0_wp, af = 0.3_wp, ws = 2.0_wp, bm = 0.5_wp
      real(wp) :: psi_tlp, r, psi, cnum, cana, dp, w1, w2
      print '(a)', '-- PV curve --'
      !----- Turgor loss point (Bartlett eqn 1). --------------------------------------!
      psi_tlp = pv_psi_tlp(pi0, eps)
      call check('psi_tlp = pi0*eps/(pi0+eps)', psi_tlp, pi0*eps/(pi0+eps), 1.0e-12_wp)
      !----- psi<->RWC round trip in the turgid region. -------------------------------!
      r   = 0.9_wp
      psi = psi_from_rwc(r, pi0, eps)
      call check('rwc round-trip (turgid)', rwc_from_psi(psi, pi0, eps), r, 1.0e-9_wp)
      !----- ... and in the flaccid region. -------------------------------------------!
      r   = 0.5_wp*((pi0+eps)/eps)
      psi = psi_from_rwc(r, pi0, eps)
      call check('rwc round-trip (flaccid)', rwc_from_psi(psi, pi0, eps), r, 1.0e-9_wp)
      !----- Full-turgor water = saturated water. -------------------------------------!
      call check('W(psi=0) = water_sat*biomass', water_content(0.0_wp, pi0, eps, af, ws, bm),  &
                 ws*bm, 1.0e-9_wp)
      !----- Capacitance = finite-difference dW/dpsi. ---------------------------------!
      psi  = -0.8_wp ; dp = 1.0e-5_wp
      w1   = water_content(psi - dp, pi0, eps, af, ws, bm)
      w2   = water_content(psi + dp, pi0, eps, af, ws, bm)
      cnum = (w2 - w1)/(2.0_wp*dp)
      cana = capacitance(psi, pi0, eps, af, ws, bm)
      call check('C(psi) ~ dW/dpsi (turgid)', cana, cnum, 1.0e-6_wp)
   end subroutine test_pv_curve

   !=======================================================================================!
   subroutine test_kirchhoff_edge()
      real(wp), parameter :: psi50 = -2.0_wp, kexp = 2.0_wp, kcond = 3.0e-3_wp
      real(wp) :: psi_l, psi_w, avg_plc, keff, keff0, pmid
      print '(a)', '-- Kirchhoff edge --'
      psi_l = -1.2_wp ; psi_w = -0.4_wp
      keff  = kirchhoff_edge(psi_w, psi_l, kcond, psi50, kexp)
      avg_plc = keff/kcond
      !----- <plc> lies between the two endpoint retained fractions. -------------------!
      call check_true('<plc> in (plc(leaf), plc(wood))',                                       &
           avg_plc > plc_retained(psi_l,psi50,kexp) .and. avg_plc < plc_retained(psi_w,psi50,kexp))
      !----- <plc> = (Phi(up)-Phi(down))/(up-down). -----------------------------------!
      call check('<plc> = DeltaPhi/Deltapsi', avg_plc,                                          &
           (flux_potential(psi_w,psi50,kexp)-flux_potential(psi_l,psi50,kexp))/(psi_w-psi_l),  &
           1.0e-12_wp)
      !----- |dpsi|->0 limit is the pointwise value at the mean. -----------------------!
      pmid  = -0.9_wp
      keff0 = kirchhoff_edge(pmid, pmid - 1.0e-9_wp, kcond, psi50, kexp)
      call check('Kirchhoff dpsi->0 limit', keff0/kcond, plc_retained(pmid,psi50,kexp), 1.0e-6_wp)
   end subroutine test_kirchhoff_edge

   !=======================================================================================!
   subroutine test_no_flow_equilibrium()
      type(hydro_params_t) :: p ; type(hydro_env_t) :: env ; type(hydro_opts_t) :: o
      type(hydro_flux_t) :: flux
      real(wp) :: psi(N_HYDRO), g, psi_l0, psi_w0
      print '(a)', '-- No-flow equilibrium (conservation) --'
      call defaults(p, env, o)
      env%transp = 0.0_wp
      g = grav_head*env%height
      psi_w0 = env%soil_psi ; psi_l0 = env%soil_psi - g       ! hydrostatic equilibrium
      psi(:) = 0.0_wp ; psi(NODE_LEAF) = psi_l0 ; psi(NODE_WOOD) = psi_w0
      call plant_water_flux(env, p, o, 900.0_wp, psi, flux)
      call check('psi_leaf unchanged at equilibrium', psi(NODE_LEAF), psi_l0, 1.0e-6_wp)
      call check('psi_wood unchanged at equilibrium', psi(NODE_WOOD), psi_w0, 1.0e-6_wp)
      call check('root_uptake ~ 0 at equilibrium',    flux%root_uptake, 0.0_wp, 1.0e-9_wp)
      call check('sapflow ~ 0 at equilibrium',        flux%sapflow,     0.0_wp, 1.0e-9_wp)
   end subroutine test_no_flow_equilibrium

   !=======================================================================================!
   subroutine test_steady_ohm()
      type(hydro_params_t) :: p ; type(hydro_env_t) :: env ; type(hydro_opts_t) :: o
      type(hydro_flux_t) :: flux
      real(wp) :: psi(N_HYDRO)
      integer(ik) :: it
      print '(a)', '-- Steady state / Ohm law --'
      call defaults(p, env, o)
      psi(:) = 0.0_wp ; psi(NODE_LEAF) = env%soil_psi ; psi(NODE_WOOD) = env%soil_psi
      do it = 1_ik, 400_ik                                     ! march to steady state
         call plant_water_flux(env, p, o, 900.0_wp, psi, flux)
      end do
      !----- Wood balance is exact at steady state (rhizo_cond constant): psi_W = psi_soil - E/K_S.
      call check('psi_wood = psi_soil - E/rhizo_cond', psi(NODE_WOOD),                          &
                 env%soil_psi - env%transp/env%rhizo_cond, 1.0e-4_wp)
      call check('root_uptake -> E at steady state', flux%root_uptake, env%transp, 1.0e-8_wp)
      call check_true('leaf more tensioned than wood', psi(NODE_LEAF) < psi(NODE_WOOD))
      call check_true('wood more tensioned than soil', psi(NODE_WOOD) < env%soil_psi)
   end subroutine test_steady_ohm

   !=======================================================================================!
   subroutine test_vs_reference()
      type(hydro_params_t) :: p ; type(hydro_env_t) :: env ; type(hydro_opts_t) :: o
      type(hydro_flux_t) :: flux
      real(wp) :: psi(N_HYDRO), refL, refW
      print '(a)', '-- vs fine explicit reference --'
      call defaults(p, env, o)
      o%rtol = 1.0e-5_wp ; o%atol = 1.0e-5_wp
      psi(:) = 0.0_wp ; psi(NODE_LEAF) = env%soil_psi ; psi(NODE_WOOD) = env%soil_psi
      call plant_water_flux(env, p, o, 900.0_wp, psi, flux)
      call reference_2node(p, env, env%soil_psi, env%soil_psi, 900.0_wp, 200000_ik, refL, refW)
      call check('psi_leaf vs reference', psi(NODE_LEAF), refL, 2.0e-3_wp)
      call check('psi_wood vs reference', psi(NODE_WOOD), refW, 2.0e-3_wp)
   end subroutine test_vs_reference

   !----- Fine forward-Euler reference for the 2-node ODE (coefficients updated each step). ---!
   subroutine reference_2node(p, env, psiL0, psiW0, dt, nstep, psiL, psiW)
      type(hydro_params_t), intent(in)  :: p
      type(hydro_env_t),    intent(in)  :: env
      real(wp),             intent(in)  :: psiL0, psiW0, dt
      integer(ik),          intent(in)  :: nstep
      real(wp),             intent(out) :: psiL, psiW
      real(wp) :: dtr, g, kc, bwood, cL, cW, keff, dL, dW
      integer(ik) :: i
      dtr   = dt/real(nstep, wp)
      g     = grav_head*env%height
      kc    = p%k_plant_max*env%leaf_area
      bwood = env%bsap + env%broot
      psiL  = psiL0 ; psiW = psiW0
      do i = 1_ik, nstep
         cL   = capacitance(psiL, p%leaf_pi0, p%leaf_eps, p%leaf_af, p%leaf_water_sat, env%bleaf)
         cW   = capacitance(psiW, p%wood_pi0, p%wood_eps, p%wood_af, p%wood_water_sat, bwood)
         keff = kirchhoff_edge(psiW, psiL, kc, p%wood_psi50, p%wood_kexp)
         dL   = (keff*(psiW - psiL - g) - env%transp)/cL
         dW   = (keff*(psiL - psiW + g) + env%rhizo_cond*(env%soil_psi - psiW))/cW
         psiL = psiL + dtr*dL
         psiW = psiW + dtr*dW
      end do
   end subroutine reference_2node

   !=======================================================================================!
   subroutine test_diurnal()
      type(hydro_params_t) :: p ; type(hydro_env_t) :: env ; type(hydro_opts_t) :: o
      type(hydro_flux_t) :: flux
      real(wp) :: psi(N_HYDRO), tsec, emax, psi_midday, psi_predawn, g
      real(wp), parameter :: pi = 3.14159265358979_wp, day = 86400.0_wp
      integer(ik) :: it
      print '(a)', '-- Diurnal cycle --'
      call defaults(p, env, o)
      emax = 2.0e-4_wp
      g    = grav_head*env%height
      psi(:) = 0.0_wp ; psi(NODE_LEAF) = env%soil_psi - g ; psi(NODE_WOOD) = env%soil_psi
      psi_midday = 0.0_wp ; psi_predawn = -1.0e30_wp
      do it = 1_ik, 96_ik                                     ! 2 days at 30-min steps
         tsec = real(mod(it-1_ik, 48_ik), wp)*1800.0_wp        ! seconds into the day
         env%transp = emax*max(0.0_wp, -cos(2.0_wp*pi*tsec/day))   ! 0 at night, peak midday
         call plant_water_flux(env, p, o, 1800.0_wp, psi, flux)
         if (mod(it-1_ik, 48_ik) == 24_ik) psi_midday  = psi(NODE_LEAF)   ! noon of each day
         if (mod(it-1_ik, 48_ik) == 0_ik)  psi_predawn = psi(NODE_LEAF)   ! predawn
      end do
      call check_true('midday psi_leaf below predawn', psi_midday < psi_predawn)
      call check_true('psi_leaf stays bounded (<= 0)', psi(NODE_LEAF) <= 0.0_wp)
      call check_true('psi_leaf finite', abs(psi(NODE_LEAF)) < 1.0e3_wp)
   end subroutine test_diurnal

   !=======================================================================================!
   subroutine test_degenerate()
      type(hydro_params_t) :: p ; type(hydro_env_t) :: env ; type(hydro_opts_t) :: o
      type(hydro_flux_t) :: flux
      real(wp) :: psi(N_HYDRO)
      print '(a)', '-- Degenerate (near-leafless) --'
      call defaults(p, env, o)
      env%bleaf = 1.0e-8_wp ; env%leaf_area = 1.0e-8_wp ; env%transp = 0.0_wp
      psi(:) = 0.0_wp ; psi(NODE_LEAF) = -0.5_wp ; psi(NODE_WOOD) = -0.4_wp
      call plant_water_flux(env, p, o, 900.0_wp, psi, flux)
      call check_true('leafless psi_leaf finite', abs(psi(NODE_LEAF)) < 1.0e3_wp)
      call check_true('leafless psi_wood finite', abs(psi(NODE_WOOD)) < 1.0e3_wp)
      call check_true('leafless fluxes finite',                                                 &
           abs(flux%root_uptake) < 1.0e3_wp .and. abs(flux%sapflow) < 1.0e3_wp)
   end subroutine test_degenerate

end program test_plant_hydraulics
