! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_plant_energy -- a tissue's energy balance (meds_plant_energy: tissue_energy_balance; !
! docs/science/vegetation_energy_dynamics.md eqs. (1)-(8)).                                 !
!   1. EXACT RELAXATION (eqs. (3)-(7)): the closed-form tissue step, its diagnostic and    !
!      frozen limits, its energy identity and small-x branch, marched exactly against the  !
!      analytic relaxation; a store with no coupling warms by its forcing.                 !
!   2. LATENT PATHWAYS: the stomata's and the wetted film's; f_wet = 0 leaves no trace of  !
!      the film, f_wet = 1 is the film alone, and partial wetting closes the balance.      !
!   3. THIN TISSUE (eqs. (1), (8)): energy closes at every area index, halving the area    !
!      halves every flux and keeps the temperature, and a tissue with no area takes the    !
!      temperature of what it grows from.                                                  !
!==========================================================================================!
program test_plant_energy
   use meds_test_assert,  only : check, check_true, test_report
   use meds_kinds,        only : wp, ik
   use meds_constants,    only : cp_air, latent_heat_vap, stefan_boltzmann, pi
   use meds_plant_energy, only : tissue_energy_balance, latent_pathway_t
   implicit none

   call test_veg_exponential()
   call test_veg_energy_diagnostic_wetted()
   call test_veg_energy_thin_tissue()

   call test_report('test_plant_energy')

contains

   !=======================================================================================!
   ! EXACT EXPONENTIAL RELAXATION in tissue_energy_balance.                               !
   !                                                                                          !
   ! Under the Category-0 freeze the tissue equation is linear with tau = cap/denom, so the     !
   ! step has a closed form and the kernel uses it: endpoint weight exp(-x), flux weight        !
   ! (1-exp(-x))/x, x = denom/store_hcap_per_dt (docs/science/vegetation_energy_dynamics.md     !
   ! eqs. (3)-(7)). These checks pin the properties that make it usable -- the two limits, the  !
   ! conservation identity, the small-x branch, and (the one that matters) that marching it is  !
   ! EXACT, not merely convergent.                                                              !
   !=======================================================================================!
   subroutine test_veg_exponential()
      real(wp) :: dtt, ts, tr, dh, drn, dt_diag_ref, ts_ref
      real(wp) :: t_cas, t0, sw, lw, h, les, lws, ler, cap, dt, denom, tau, a_st
      real(wp) :: t_exact, t_march, worst, fe
      type(latent_pathway_t) :: dry, none   !< the stomata's pathway; none: no film
      integer(ik) :: k, nstep
      print '(a)', 'test_veg_exponential:'
      t_cas = 295.0_wp
      sw = 300.0_wp
      lw = -40.0_wp
      h = 120.0_wp
      les = 30.0_wp
      lws = 5.0_wp
      ler = 25.0_wp
      dry = latent_pathway_t(slope=les, ref=ler)
      denom = h + les + lws
      dt    = 1800.0_wp

      !----- (a) store_hcap_per_dt -> 0 is the DIAGNOSTIC limit, exactly. This is why there is no     !
      !          diagnostic/prognostic selector: diagnostic is a limit of one formula, not a second   !
      !          branch. -----------------------------------------------------------------------------!
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, 0.0_wp, t_cas, &
                                 dtt, ts, tr, fe, dh, drn)
      dt_diag_ref = (sw + lw - ler) / denom
      call check('store_hcap_per_dt = 0 gives the exact diagnostic offset', dtt, dt_diag_ref, 1.0e-12_wp)
      ts_ref = t_cas + dt_diag_ref

      !----- (b) store_hcap_per_dt -> huge FREEZES the store at its entry temperature. -----------------------!
      t0 = 288.0_wp
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, 1.0e12_wp, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('store_hcap_per_dt -> infinity holds the store at t_store0', ts, t0, 1.0e-6_wp)

      !----- (c) THE CONSERVATION IDENTITY. The endpoint and step-average weights are different      !
      !          numbers, and pairing them correctly is exactly what makes the balance close:        !
      !            store_hcap_per_dt*(dt_end - dt_prev) + denom*dt_avg == numer.                     !
      !          dt_avg is not returned, so recover it from dh = h*dt_avg. ------------------------!
      cap  = 4.0e4_wp                       ! a big cohort: tau = cap/denom is order dt
      a_st = cap / dt
      tau  = cap / denom
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('energy balance closes with the endpoint/average pair',                            &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-9_wp)

      !----- (d) THE POINT: marching the kernel reproduces the analytic relaxation EXACTLY, not      !
      !          approximately. A backward-Euler storage term would drift here; this does not.       !
      !          T(t) = T_eq + (T0 - T_eq)*exp(-t/tau).                                              !
      nstep = 20_ik
      t_march = t0
      worst = 0.0_wp
      do k = 1_ik, nstep
         call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t_march, &
                                    dtt, ts, tr, fe, dh, drn)
         t_march = ts
         t_exact = ts_ref + (t0 - ts_ref) * exp(-real(k, wp) * dt / tau)
         worst   = max(worst, abs(t_march - t_exact))
      end do
      call check_true('marching 20 steps is EXACT against the analytic relaxation', worst < 1.0e-10_wp, worst)

      !----- (e) The small-x series branch (a very large capacity) must not lose precision where     !
      !          1 - exp(-x) would cancel. tau = 1e6 s against dt = 1800 s puts x ~ 2.7e-4. ---------!
      a_st = 1.0e6_wp * denom / dt
      call tissue_energy_balance(sw, lw, 0.0_wp, h, lws, dry, none, 0.0_wp, t_cas, t_cas, a_st, t0, &
                                 dtt, ts, tr, fe, dh, drn)
      call check('large-capacity limit stays accurate (series branch)',                             &
           a_st*(dtt - (t0 - t_cas)) + denom*(dh/h), sw + lw - ler, 1.0e-6_wp)

      !----- (f) A store with NO coupling (x = 0) warms by its forcing over the step, numer*dt/cap,  !
      !          and passes nothing on: the relaxation needs no division by the coupling. ---------!
      a_st = cap / dt
      call tissue_energy_balance(sw, lw, 0.0_wp, 0.0_wp, 0.0_wp, latent_pathway_t(ref=ler), none, 0.0_wp, t_cas, &
                                 t_cas, a_st, t0, dtt, ts, tr, fe, dh, drn)
      call check('no coupling: the store warms by numer*dt/cap', ts, t0 + (sw + lw - ler) / a_st, 1.0e-12_wp)
      call check_true('no coupling: no sensible heat leaves', dh == 0.0_wp, dh)
   end subroutine test_veg_exponential

   subroutine test_veg_energy_diagnostic_wetted()
      !----- Direct unit test of tissue_energy_balance's two latent pathways (the column tests exercise  !
      !      them only through the full integration). f_wet splits the latent flux into a (1-f_wet)     !
      !      DRY (stomatal) share and an f_wet WET (boundary-layer film) share; DRY and WET carry two    !
      !      DIFFERENT conductances, so the three properties below are non-trivial (a bug that swapped  !
      !      or mixed them would still "look" plausible).                                               !
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_temp2, t_store2, transp2, dh2, drnet2, film_evap2
      real(wp), parameter :: ABS_SW = 400.0_wp, ABS_LW = -50.0_wp, H_COEFF = 12.0_wp
      real(wp), parameter :: LW_SLOPE = 6.0_wp, T_CAS = 298.0_wp, T_PARENT = 298.0_wp
      type(latent_pathway_t), parameter :: DRY = latent_pathway_t(slope=8.0_wp,  ref=15.0_wp)   ! stomatal
      type(latent_pathway_t), parameter :: WET = latent_pathway_t(slope=20.0_wp, ref=40.0_wp)   ! film
      type(latent_pathway_t) :: none
      print '(a)', 'test_veg_energy_diagnostic_wetted:'

      !----- 1. f_wet = 0 leaves no trace of the film's pathway: the same result, bit for bit, as no   !
      !      film pathway at all. -----------------------------------------------------------------!
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, none, 0.0_wp,           &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 0.0_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp2, t_store2, transp2, film_evap2, dh2, &
                                 drnet2)
      call check('f_wet=0 == no film pathway (dt_temp)', dt_temp2, dt_temp, 1.0e-14_wp)
      call check('f_wet=0 == no film pathway (transp)',  transp2,  transp,  1.0e-14_wp)
      call check_true('f_wet=0 -> film_evap = 0 exactly', film_evap2 == 0.0_wp, film_evap2)

      !----- 2. FULLY WET (f_wet=1): the dry (stomatal) pathway must vanish -- a saturated film          !
      !      transpires ~0 and evaporates its whole latent capacity through the wet pathway instead.      !
      !      SWAP-EQUIVALENCE: this must be algebraically IDENTICAL to calling the kernel with the film's  !
      !      pathway as the DRY one and no film -- not just approximately equal (both reduce to the same   !
      !      2-unknown linear solve). ---------------------------------------------------------------------!
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 1.0_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, WET, none, 0.0_wp,           &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp2, t_store2, transp2, film_evap2, dh2, &
                                 drnet2)
      call check_true('f_wet=1: dry pathway (transp) vanishes', abs(transp) < 1.0e-14_wp, transp)
      call check('f_wet=1: film_evap == swap-equivalent transp', film_evap, transp2, 1.0e-12_wp)
      call check('f_wet=1: dt_temp == swap-equivalent dt_temp',  dt_temp,   dt_temp2, 1.0e-12_wp)

      !----- 3. PARTIAL wetting (f_wet=0.4): the energy balance must still close -- absorbed radiation     !
      !      net of LW emission (drnet) splits EXACTLY into sensible + BOTH latent pathways, each at the    !
      !      pathway's h_vap (here the default, latent_heat_vap) -- the identity the surface path relies    !
      !      on to make the CAS's temperature-dependent vapour-enthalpy credit balance against a store's     !
      !      liquid-enthalpy debit (for transpiration the leaf pays the full vapour enthalpy itself; the     !
      !      surface-water accounting does the same for film_evap). -----------------------------------------!
      call tissue_energy_balance(ABS_SW, ABS_LW, 0.0_wp, H_COEFF, LW_SLOPE, DRY, WET, 0.4_wp,            &
                                 T_CAS, T_PARENT, 0.0_wp, T_CAS, dt_temp, t_store, transp, film_evap, dh, drnet)
      call check('partial f_wet: energy balance closes (drnet = dh + latent)', drnet,                  &
                 dh + (transp + film_evap) * latent_heat_vap, 1.0e-9_wp)
      call check_true('partial f_wet: film_evap > 0 (wet pathway active)', film_evap > 0.0_wp, film_evap)
      call check_true('partial f_wet: transp > 0 (dry pathway still active)', transp > 0.0_wp, transp)
   end subroutine test_veg_energy_diagnostic_wetted


   subroutine test_veg_energy_thin_tissue()
      !----- Every term of a tissue's balance -- absorbed radiation, the sapflow's enthalpy, the        !
      !      conductances -- and its heat capacity scale with the tissue's own area                     !
      !      (docs/science/vegetation_energy_dynamics.md eqs. (1), (8)), so a thin tissue is simply a    !
      !      small one. Three properties, swept over an area index from 1e-6 to 10 for leaf and wood,    !
      !      with and without a heat store:                                                             !
      !                                                                                          !
      !        * CONSERVATION: within the step everything that enters leaves or is stored,               !
      !            drnet + q_extra == dh + Lv*(transp + film_evap) + store_hcap_per_dt*(t_store - t0),   !
      !          to round-off. surface_derivs relies on it: the whole-column ledger books coh_rnet on    !
      !          the input side and coh_h/coh_qw as they reach the CAS, and nothing re-checks the        !
      !          tissue in between.                                                                      !
      !        * SPLIT INVARIANCE: halving the area (one cohort of two at half density) keeps the        !
      !          temperature and halves every flux. A floor in ground units on the coupling or the       !
      !          capacity breaks this for a thin tissue, pulling it toward the canopy air.               !
      !        * NO AREA: a tissue with no coupling and no heat capacity takes the temperature of what   !
      !          it grows from (t_parent) and carries no flux. -------------------------------------------!
      real(wp) :: dt_temp, t_store, transp, dh, drnet, film_evap
      real(wp) :: dt_h, t_h, transp_h, dh_h, drnet_h, film_h
      real(wp) :: ai, ast, resid, worst, worst_split
      integer(ik) :: ia, il, is
      real(wp), parameter :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp), parameter :: CAP_PER_AREA = 1500.0_wp   !< [J/m2 tissue/K] a leaf's dry tissue and water
      logical :: is_leaf
      print '(a)', 'test_veg_energy_thin_tissue:'
      worst = 0.0_wp
      worst_split = 0.0_wp
      do is = 0_ik, 1_ik
         do il = 1_ik, 2_ik
            is_leaf = (il == 1_ik)
            do ia = 1_ik, 15_ik
               ai = 1.0e-6_wp * (10.0_wp ** (0.5_wp * real(ia - 1_ik, wp)))
               ast = 0.0_wp
               if (is == 1_ik) ast = CAP_PER_AREA * ai / 900.0_wp
               call thin_tissue(ai, ast, is_leaf, dt_temp, t_store, transp, dh, drnet, film_evap)
               resid = drnet + 25.0_wp*ai - dh - latent_heat_vap*(transp + film_evap) - ast * (t_store - T0)
               worst = max(worst, abs(resid) / (abs(drnet) + 25.0_wp*ai))
               call thin_tissue(0.5_wp*ai, 0.5_wp*ast, is_leaf, dt_h, t_h, transp_h, dh_h, drnet_h, film_h)
               worst_split = max(worst_split, abs(t_h - t_store), abs(2.0_wp*dh_h - dh) / abs(dh),          &
                                 abs(2.0_wp*drnet_h - drnet) / abs(drnet))
            end do
         end do
      end do
      call check_true('a thin tissue conserves energy at every area index (60 cases)', worst < 1.0e-12_wp, worst)
      call check_true('half the area: the same temperature, half of every flux', worst_split < 1.0e-12_wp,       &
                      worst_split)
      !----- No area: no coupling and no heat capacity. -----------------------------------------!
      call thin_tissue(0.0_wp, 0.0_wp, .true., dt_temp, t_store, transp, dh, drnet, film_evap)
      call check('a tissue with no area takes the temperature of what it grows from', t_store, T_PARENT, 0.0_wp)
      call check_true('a tissue with no area carries no flux',                                                  &
                      dh == 0.0_wp .and. drnet == 0.0_wp .and. transp == 0.0_wp .and. film_evap == 0.0_wp, dh)
   end subroutine test_veg_energy_thin_tissue

   !----- One tissue of area index a at 298 K canopy air, starting 0.5 K below it, whose parent sits at  !
   !      301 K: every coupling and forcing term proportional to a. -----------------------------------!
   subroutine thin_tissue(a, a_store, leaf, dtt, ts, tr, h_out, rn, fe)
      real(wp), intent(in)  :: a, a_store
      logical,  intent(in)  :: leaf
      real(wp), intent(out) :: dtt, ts, tr, h_out, rn, fe
      real(wp), parameter   :: GBH = 0.03_wp, RHO = 1.2_wp, EMISS = 0.97_wp
      real(wp), parameter   :: T_CAS = 298.0_wp, T0 = T_CAS - 0.5_wp, T_PARENT = 301.0_wp
      real(wp) :: hc, lws, les, ler
      hc  = merge(2.0_wp, pi, leaf) * a * GBH * RHO * cp_air
      lws = 4.0_wp * EMISS * stefan_boltzmann * T_CAS**3 * a
      les = merge( 8.0_wp * a, 0.0_wp, leaf)    ! stomatal (wood does not transpire)
      ler = merge(15.0_wp * a, 0.0_wp, leaf)
      call tissue_energy_balance(400.0_wp*a, -50.0_wp*a, 25.0_wp*a, hc, lws, latent_pathway_t(slope=les, ref=ler), &
                                 latent_pathway_t(slope=20.0_wp*a, ref=40.0_wp*a), 0.4_wp, T_CAS, T_PARENT,      &
                                 a_store, T0, dtt, ts, tr, fe, h_out, rn)
   end subroutine thin_tissue

end program test_plant_energy
