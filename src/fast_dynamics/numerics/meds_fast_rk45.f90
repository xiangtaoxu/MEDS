! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_rk45 -- the explicit fast-loop integrator a run can choose ("rk45"): an        !
! adaptive Cash-Karp 5(4) march over the whole column state (canopy air, soil energy and   !
! water, plant water, films) through the shared right-hand side column_derivs, from the    !
! same frozen record as ARK. It is also the explicit reference the tests compare ARK       !
! against. A step it cannot resolve is redone on ARK by meds_fast_step. Main routines:     !
! column_fast_step_rk45, rk45_column_step (docs/science/numerical_scheme.md §3 `rk45`, §5b).!
!==========================================================================================!
module meds_fast_rk45
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, cp_air, n_soil_layer_max
   use meds_therm_lib,           only : cas_enthalpy_to_temp, internal_energy_liquid
   use meds_fast_time_derivs, only : column_derivs, cas_conductances
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t, error_control_t,             &
                                     surface_tend_t, column_config_t, column_cohort_t, column_forcing_t,          &
                                     column_budget_t, patch_biophys_t, fast_workspace_t, column_site_t
   use meds_fast_frozen,      only : build_column_frozen
   use meds_column_state_ops, only : state_init, state_accum, state_sub, zero_like, clamp_cas,                    &
                                     soil_water_store, soil_energy_store, plant_water_store, canopy_film_store,   &
                                     deposit_condensate, clamp_canopy_film, unpack_column_state,                  &
                                     diagnose_soil_temps, transp_per_leaf_area, commit_tissue_energy
   use meds_diag_vars,        only : CO_TRANSP
   use meds_fast_control,     only : state_wrms_grouped,                                         &
                                     step_control_factor
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_budget_check,     only : budget_check, last_check_closed, BUDGET_ENERGY_RATE_FLOOR,      &
                                     BUDGET_WATER_RATE_FLOOR
   implicit none
   private

   public :: column_fast_step_rk45
   public :: rk45_state_railed, rk45_ledgers_stop
   public :: rk45_column_step   ! used only by the tests

   !----- Cash-Karp 5(4) tableau (Cash & Karp 1990, ACM TOMS 16:201; Numerical Recipes' rkck). The    !
   !      frozen record has no explicit time dependence, so the stage times never enter. ------------!
   real(wp), parameter :: A21 = 0.2_wp
   real(wp), parameter :: A31 = 3.0_wp/40.0_wp,      A32 = 9.0_wp/40.0_wp
   real(wp), parameter :: A41 = 0.3_wp,               A42 = -0.9_wp,             A43 = 1.2_wp
   real(wp), parameter :: A51 = -11.0_wp/54.0_wp,     A52 = 2.5_wp,              A53 = -70.0_wp/27.0_wp, &
                          A54 = 35.0_wp/27.0_wp
   real(wp), parameter :: A61 = 1631.0_wp/55296.0_wp, A62 = 175.0_wp/512.0_wp,   A63 = 575.0_wp/13824.0_wp, &
                          A64 = 44275.0_wp/110592.0_wp, A65 = 253.0_wp/4096.0_wp
   !----- 5th-order (committed) and embedded 4th-order weights; b2 = 0 in both, b5 = 0 at 5th order. !
   real(wp), parameter :: B1 = 37.0_wp/378.0_wp,    B3 = 250.0_wp/621.0_wp, B4 = 125.0_wp/594.0_wp, &
                          B6 = 512.0_wp/1771.0_wp
   real(wp), parameter :: BS1 = 2825.0_wp/27648.0_wp, BS3 = 18575.0_wp/48384.0_wp,                &
                          BS4 = 13525.0_wp/55296.0_wp, BS5 = 277.0_wp/14336.0_wp, BS6 = 0.25_wp
   !----- The embedded estimate's order, p in docs/science/numerical_scheme.md eq. (4). ----------------!
   integer(ik), parameter :: RK45_P_ORDER = 4_ik

   !----- Work cap per dt_fast. A healthy step takes a few sub-steps (docs/science/numerical_scheme.md   !
   !      §3 `rk45`); on a dense, cold canopy the explicit march instead thrashes and still rails, and   !
   !      the dispatcher redoes the step on ARK, so bailing at the cap bounds the wasted work. ---------!
   integer(ik), parameter :: RK45_WORK_CAP = 64_ik

contains

   !---------------------------------------------------------------------------------------!
   ! stage_bnd -- the whole-column boundary-flux rates at one stage, the explicit counterpart of      !
   ! column_be_stage's stage_bflux_t: the state is the stage input y_stage, and surf_tend is          !
   ! column_derivs' surface diagnostic there, condensation included at the stage's own CAS temperature. !
   pure subroutine stage_bnd(y_stage, frozen, surf_tend, rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i, &
                             atm_heat_i)
      type(column_state_t),  intent(in)  :: y_stage
      type(column_frozen_t), intent(in)  :: frozen
      type(surface_tend_t),  intent(in)  :: surf_tend
      real(wp),               intent(out) :: rnet_i, atm_enth_i, atm_vap_i, cond_i, cond_enth_i
      real(wp),               intent(out) :: atm_heat_i   !< [W/m2] sensible export, the reported H
      real(wp) :: gah_i, gaw_i, gac_i
      rnet_i      = surf_tend%coh_rnet
      !----- Charge the boundary at the conductances the tendency used, re-solved at the stage's state, !
      !      not the state^n ones (docs/science/numerical_scheme.md §5b). cas_conductances is pure, so   !
      !      this matches column_derivs; carrying them on surface_tend_t instead would perturb ifx's     !
      !      inlining and FMA choices in surface_derivs (~1e-12 in every result). ----------------------!
      call cas_conductances(frozen%cas, y_stage%cas_enthalpy, y_stage%cas_shv, gah_i, gaw_i, gac_i)
      atm_enth_i  = gah_i * (y_stage%cas_enthalpy - frozen%cas%enthalpy_atm)
      atm_vap_i   = gaw_i * (y_stage%cas_shv      - frozen%cas%shv_atm)
      atm_heat_i  = gah_i * cp_air * (cas_enthalpy_to_temp(y_stage%cas_enthalpy, y_stage%cas_shv)          &
                                       - frozen%cas%mo_theta_atm)
      cond_i      = surf_tend%cond
      cond_enth_i = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
   end subroutine stage_bnd

   !---------------------------------------------------------------------------------------!
   ! rk45_column_step -- one Cash-Karp 5(4) step over column_derivs, committing the 5th-order    !
   ! solution; y_err is the 5th - 4th difference. w_out/e_in/e_out are the whole-column boundary  !
   ! amounts over dt, b-weighted like the commit; infiltration, a pond -> soil transfer, is not    !
   ! one. Each stage input's canopy air is clamped (clamp_cas) so that an oversized trial gives a  !
   ! finite derivative, which the embedded error then rejects. Soil energy needs no clamp: any     !
   ! finite energy has a temperature, and a stage reads it only through terms linear in it.       !
   !---------------------------------------------------------------------------------------!
   pure subroutine rk45_column_step(y, frozen, col_site, n, nsl, dt, y_out, y_err, w_out, e_in, e_out,          &
                                    clamp_stage_n, cond_out,                                         &
                                    tissue_leaf_int, tissue_wood_int, transp_int, cond_enth_out,       &
                                    atm_heat_out, atm_vap_out, face_resid_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_site_t),   intent(in)  :: col_site
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      real(wp),               intent(out) :: w_out, e_in, e_out
      !----- Stage-clamp count (column_budget_t%clamp_stage_n), expected on an oversized trial. ------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n
      real(wp),    optional,  intent(out)   :: cond_out   !< [kg/m2] condensate to deposit into soil layer 1
      real(wp),    optional,  intent(out)   :: cond_enth_out !< [J/m2] its liquid enthalpy, b-weighted at the stage CAS temps
      real(wp), optional,     intent(out) :: atm_heat_out, atm_vap_out   !< [J/m2],[kg/m2] net CAS->atm sensible/vapour over dt
      real(wp),    optional,  intent(out)   :: tissue_leaf_int(n), tissue_wood_int(n)  !< [K*s]
      real(wp),    optional,  intent(out)   :: transp_int(n)   !< [kg/m2] each cohort's transpiration over dt
      !----- Per-layer face closure: summed over interior layers, |the water the committed theta moved - !
      !      the water the soil-energy faces were charged for| (docs/science/numerical_scheme.md §5b). ---!
      real(wp),    optional,  intent(out)   :: face_resid_out   !< [kg/m2]

      type(column_tend_t)  :: k1, k2, k3, k4, k5, k6
      type(column_state_t) :: y_stage, y_4th
      type(surface_tend_t) :: surf_tend
      real(wp) :: rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6)
      real(wp) :: bw_rnet, bw_atm_enth, bw_atm_vap, bw_cond, bw_cond_enth, bw_drain, bw_atm_heat, bw_g_bottom
      !----- Per-stage tissue temperatures, for the time integral the tissues' energy is set from. ---!
      real(wp) :: tleaf_s(n,6), twood_s(n,6), transp_s(n,6)
      real(wp) :: bw_face(n_soil_layer_max), f_in, f_out, uptk, face_resid
      integer(ik) :: kf

      call column_derivs(y, frozen, col_site, n, nsl, k1, sf_out=surf_tend)
      call stage_bnd(y, frozen, surf_tend, rnet(1), atm_enth(1), atm_vap(1), cond(1), cond_enth(1), atm_heat(1))
      tleaf_s(1:n,1) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,1) = surf_tend%wood_temp(1:n)
      transp_s(1:n,1) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A21, k1, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_site, n, nsl, k2, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(2), atm_enth(2), atm_vap(2), cond(2), cond_enth(2), atm_heat(2))
      tleaf_s(1:n,2) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,2) = surf_tend%wood_temp(1:n)
      transp_s(1:n,2) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A31, k1, n, nsl)
      call state_accum(y_stage, dt*A32, k2, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_site, n, nsl, k3, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(3), atm_enth(3), atm_vap(3), cond(3), cond_enth(3), atm_heat(3))
      tleaf_s(1:n,3) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,3) = surf_tend%wood_temp(1:n)
      transp_s(1:n,3) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A41, k1, n, nsl)
      call state_accum(y_stage, dt*A42, k2, n, nsl)
      call state_accum(y_stage, dt*A43, k3, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_site, n, nsl, k4, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(4), atm_enth(4), atm_vap(4), cond(4), cond_enth(4), atm_heat(4))
      tleaf_s(1:n,4) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,4) = surf_tend%wood_temp(1:n)
      transp_s(1:n,4) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A51, k1, n, nsl)
      call state_accum(y_stage, dt*A52, k2, n, nsl)
      call state_accum(y_stage, dt*A53, k3, n, nsl)
      call state_accum(y_stage, dt*A54, k4, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_site, n, nsl, k5, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(5), atm_enth(5), atm_vap(5), cond(5), cond_enth(5), atm_heat(5))
      tleaf_s(1:n,5) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,5) = surf_tend%wood_temp(1:n)
      transp_s(1:n,5) = surf_tend%transp_c(1:n)

      call state_init(y, n, nsl, y_stage)
      call state_accum(y_stage, dt*A61, k1, n, nsl)
      call state_accum(y_stage, dt*A62, k2, n, nsl)
      call state_accum(y_stage, dt*A63, k3, n, nsl)
      call state_accum(y_stage, dt*A64, k4, n, nsl)
      call state_accum(y_stage, dt*A65, k5, n, nsl)
      call clamp_cas(y_stage, nfire=clamp_stage_n)
      call column_derivs(y_stage, frozen, col_site, n, nsl, k6, sf_out=surf_tend)
      call stage_bnd(y_stage, frozen, surf_tend, rnet(6), atm_enth(6), atm_vap(6), cond(6), cond_enth(6), atm_heat(6))
      tleaf_s(1:n,6) = surf_tend%leaf_temp(1:n)
      twood_s(1:n,6) = surf_tend%wood_temp(1:n)
      transp_s(1:n,6) = surf_tend%transp_c(1:n)

      !----- y_out = y + dt*(B1*k1 + B3*k3 + B4*k4 + B6*k6), the 5th-order solution. -----------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt*B1, k1, n, nsl)
      call state_accum(y_out, dt*B3, k3, n, nsl)
      call state_accum(y_out, dt*B4, k4, n, nsl)
      call state_accum(y_out, dt*B6, k6, n, nsl)
      !----- Per-layer face closure on the committed step, interior layers only: the water column's   !
      !      top face nets out ground evaporation, whose enthalpy is paid at the surface rather than   !
      !      advected, and the bottom face is the drainage; the whole-column ledgers cover both. -----!
      if (present(face_resid_out)) then
         do kf = 1_ik, nsl
            bw_face(kf) = B1*k1%soil_face(kf) + B3*k3%soil_face(kf)                                &
                        + B4*k4%soil_face(kf) + B6*k6%soil_face(kf)
         end do
         face_resid = 0.0_wp
         do kf = 2_ik, nsl - 1_ik
            f_in  = bw_face(kf-1_ik) * rho_h2o * dt
            f_out = bw_face(kf)      * rho_h2o * dt
            uptk  = frozen%roots%uptake * frozen%roots%root_share(kf) * dt
            face_resid = face_resid                                                                &
               + abs((y_out%theta(kf) - y%theta(kf)) * col_site%soil_profile%dz(kf) * rho_h2o           &
                     - (f_in - f_out - uptk))
         end do
         face_resid_out = face_resid
      end if

      !----- The committed state is not clamped: at the floor a step is committed with the controller !
      !      out of moves, and a clamp there would move mass or energy with no ledger entry. A railed  !
      !      commit is caught by rk45_state_railed, and the dispatcher redoes the dt_fast on ARK. -----!

      !----- y_4th, the embedded 4th-order solution, is its own temporary: passing y_err as both      !
      !      state_sub's b and out would alias two dummies, one modified, which Fortran forbids. -----!
      call state_init(y, n, nsl, y_4th)
      call state_accum(y_4th, dt*BS1, k1, n, nsl)
      call state_accum(y_4th, dt*BS3, k3, n, nsl)
      call state_accum(y_4th, dt*BS4, k4, n, nsl)
      call state_accum(y_4th, dt*BS5, k5, n, nsl)
      call state_accum(y_4th, dt*BS6, k6, n, nsl)
      call state_sub(y_out, y_4th, n, nsl, y_err)   ! y_err := y_5th - y_4th (the embedded estimate)

      !----- Boundary amounts over dt: the stage rates b-weighted like the commit, plus the frozen     !
      !      constants added once (sum(b) = 1). Drainage is the stages' too: the network's rate. ------!
      bw_drain     = B1*k1%drainage_rate + B3*k3%drainage_rate + B4*k4%drainage_rate                  &
                     + B6*k6%drainage_rate
      !----- Heat across the column's bottom face: 0 on the adiabatic base, the conduction to         !
      !      deep_temp under the Dirichlet anchor. -----------------------------------------------------!
      bw_g_bottom  = B1*k1%g_bottom + B3*k3%g_bottom + B4*k4%g_bottom + B6*k6%g_bottom
      bw_rnet      = B1*rnet(1)      + B3*rnet(3)      + B4*rnet(4)      + B6*rnet(6)
      bw_atm_enth  = B1*atm_enth(1)  + B3*atm_enth(3)  + B4*atm_enth(4)  + B6*atm_enth(6)
      bw_atm_vap   = B1*atm_vap(1)   + B3*atm_vap(3)   + B4*atm_vap(4)   + B6*atm_vap(6)
      bw_atm_heat  = B1*atm_heat(1)  + B3*atm_heat(3)  + B4*atm_heat(4)  + B6*atm_heat(6)
      bw_cond      = B1*cond(1)      + B3*cond(3)      + B4*cond(4)      + B6*cond(6)
      bw_cond_enth = B1*cond_enth(1) + B3*cond_enth(3) + B4*cond_enth(4) + B6*cond_enth(6)
      !----- Tissue-temperature time integral over the step, b-weighted like the commit. -----------!
      if (present(tissue_leaf_int)) then
         tissue_leaf_int(1:n) = dt * (B1*tleaf_s(1:n,1) + B3*tleaf_s(1:n,3)                       &
                                      + B4*tleaf_s(1:n,4) + B6*tleaf_s(1:n,6))
         tissue_wood_int(1:n) = dt * (B1*twood_s(1:n,1) + B3*twood_s(1:n,3)                       &
                                      + B4*twood_s(1:n,4) + B6*twood_s(1:n,6))
      end if
      !----- ...and each cohort's transpiration, as its leaf water was debited. ---------------------!
      if (present(transp_int)) transp_int(1:n) = dt * (B1*transp_s(1:n,1) + B3*transp_s(1:n,3)       &
                                                     + B4*transp_s(1:n,4) + B6*transp_s(1:n,6))

      !----- No clip term: the water network keeps every layer inside its curve                       !
      !      (docs/science/numerical_scheme.md §5b). ground_rad is the snow-blended ground radiation.  !
      !      Infiltration is a pond -> soil transfer between two stores, not a boundary input. -------!
      e_in  = (bw_rnet + frozen%snow%ground_rad + bw_g_bottom) * dt
      e_out = bw_atm_enth * dt                                                                    &
              + bw_drain * dt * internal_energy_liquid(frozen%hydrology%t_bot)
      !----- Condensate is not in w_out: the caller deposits it into soil layer 1 (cond_out). The     !
      !      caller adds the network's runoff once per dt_fast. --------------------------------------!
      w_out = bw_atm_vap * dt + bw_drain * dt
      if (present(cond_out))      cond_out      = bw_cond * dt
      if (present(cond_enth_out)) cond_enth_out = bw_cond_enth * dt
      if (present(atm_heat_out))  atm_heat_out  = bw_atm_heat * dt
      if (present(atm_vap_out))   atm_vap_out   = bw_atm_vap  * dt
   end subroutine rk45_column_step

   !---------------------------------------------------------------------------------------!
   ! adaptive_rk45_march -- march to t_end over rk45_column_step with the same error norm,       !
   ! controller, warm start and accept/reject rule as adaptive_ark_march                         !
   ! (docs/science/numerical_scheme.md eqs. (3)-(4), p = 4). Every integrated state but the films !
   ! is scored; soil and plant water's embedded difference is round-off, their rates constant.   !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_rk45_march(y0, frozen, col_site, n, nsl, t_end, ec, dt_init, y_out, nsteps, nrej,       &
                                  w_out_acc, e_in_acc, e_out_acc, dt_warm_out,                    &
                                  clamp_stage_n, cond_acc,                                            &
                                  cond_enth_acc, atm_heat_acc, atm_vap_acc, face_resid_acc,           &
                                  ood_max, tissue_leaf_acc, tissue_wood_acc, transp_acc)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      type(column_site_t),   intent(in)  :: col_site
      integer(ik),            intent(in)  :: n, nsl
      real(wp),               intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec
      type(column_state_t),  intent(out) :: y_out
      integer(ik),            intent(out) :: nsteps, nrej
      real(wp),               intent(out) :: w_out_acc, e_in_acc, e_out_acc
      real(wp),    optional,  intent(out) :: dt_warm_out
      !----- Stage-clamp count over the march, rejected trials included, as on ARK. ------------------!
      integer(ik), optional,  intent(inout) :: clamp_stage_n
      !----- Per-layer face closure over the accepted sub-steps only. ---------------------------------!
      real(wp),    optional,  intent(inout) :: face_resid_acc
      real(wp),    optional,  intent(inout) :: cond_acc   !< [kg/m2] accumulated condensate
      real(wp),    optional,  intent(inout) :: cond_enth_acc !< [J/m2] ...and its stage-valued liquid enthalpy
      real(wp), optional,     intent(inout) :: atm_heat_acc, atm_vap_acc   !< [J/m2],[kg/m2] net CAS->atm export, accepted sub-steps
      !----- [K*s] per-cohort tissue-temperature time integrals over the accepted march. Caller zeroes.!
      real(wp),    optional,  intent(inout) :: tissue_leaf_acc(n), tissue_wood_acc(n)
      real(wp),    optional,  intent(inout) :: transp_acc(n)   !< [kg/m2] per-cohort transpiration, accepted march
      !----- [m3/m3] running max of theta's excursion outside its domain at a stage-1 input, every    !
      !      trial included, since each is a place the right-hand side was evaluated (theta_ood_max). -!
      real(wp),    optional,  intent(inout) :: ood_max

      type(column_state_t) :: y, y_new, y_err, y_zero
      real(wp) :: t, dt, err, err_prev, fac, dt_floor
      real(wp) :: w_out, e_in, e_out, dt_try, dt_warm
      real(wp) :: cond_i, cond_enth_i, face_i
      real(wp)    :: atm_heat_i, atm_vap_i
      real(wp) :: tl_int_i(n), tw_int_i(n), tr_int_i(n)
      integer(ik) :: kood
      logical  :: clamped

      w_out_acc = 0.0_wp
      e_in_acc = 0.0_wp
      e_out_acc = 0.0_wp
      !----- Sub-step floor t_end/64, as on ARK. A floor step is not guaranteed stable here, but a    !
      !      healthy step takes a few sub-steps, so reaching the floor marks a pathological one. ------!
      dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)   ! clamp-ok: the smallest sub-step [s] (above)

      call state_init(y0, n, nsl, y)
      t = 0.0_wp
      dt = min(dt_init, t_end)
      nsteps = 0_ik
      nrej = 0_ik
      err_prev = -1.0_wp
      dt_warm = dt
      do
         if (t >= t_end - tiny_num) exit
         dt_try = dt
         dt = min(dt, t_end - t)
         clamped = dt < dt_try - tiny_num
         !----- Record how far the step's starting theta lies outside its domain (never clamped). ----!
         if (present(ood_max)) then
            do kood = 1_ik, nsl
               ood_max = max(ood_max, y%theta(kood) - col_site%soil_profile%theta_sat(kood),                         &
                             col_site%soil_profile%theta_res(kood) - y%theta(kood))
            end do
         end if
         call rk45_column_step(y, frozen, col_site, n, nsl, dt, y_new, y_err, w_out, e_in, e_out,               &
                               clamp_stage_n=clamp_stage_n, cond_out=cond_i,                      &
                               tissue_leaf_int=tl_int_i, tissue_wood_int=tw_int_i, transp_int=tr_int_i, &
                               cond_enth_out=cond_enth_i, atm_heat_out=atm_heat_i, atm_vap_out=atm_vap_i, &
                               face_resid_out=face_i)
         !----- A named temporary, not a function result passed to a call (docs/building.md traps). -!
         y_zero = zero_like(y_err, n, nsl)
         err = state_wrms_grouped(y_err, y_zero, y_new, n, nsl, ec%tols)
         if (err /= err .or. dt /= dt) then    ! non-finite: too big a step, reject deterministically
            if (dt <= dt_floor) then
               !----- A floor-forced commit of a non-finite trial, tallied like any accept. ----------!
               if (present(face_resid_acc)) face_resid_acc = face_resid_acc + face_i
               if (present(cond_acc))       cond_acc       = cond_acc       + cond_i
               if (present(cond_enth_acc))  cond_enth_acc  = cond_enth_acc  + cond_enth_i
            if (present(atm_heat_acc))   atm_heat_acc   = atm_heat_acc   + atm_heat_i
            if (present(atm_vap_acc))    atm_vap_acc    = atm_vap_acc    + atm_vap_i
               call state_init(y_new, n, nsl, y)
               t = t + dt_floor
               nsteps = nsteps + 1_ik
               exit
            end if
            nrej = nrej + 1_ik
            dt = max(dt * ec%fmin, dt_floor)
            cycle
         end if
         fac = step_control_factor(err, err_prev, ec)
         if (err <= 1.0_wp .or. dt <= dt_floor) then
            call state_init(y_new, n, nsl, y)
            if (present(face_resid_acc)) face_resid_acc = face_resid_acc + face_i
            if (present(cond_acc))       cond_acc       = cond_acc       + cond_i
            if (present(cond_enth_acc))  cond_enth_acc  = cond_enth_acc  + cond_enth_i
            if (present(atm_heat_acc))   atm_heat_acc   = atm_heat_acc   + atm_heat_i
            if (present(atm_vap_acc))    atm_vap_acc    = atm_vap_acc    + atm_vap_i
            w_out_acc = w_out_acc + w_out
            e_in_acc = e_in_acc + e_in
            e_out_acc = e_out_acc + e_out
            !----- Over accepted sub-steps only, like every amount above. --------------------------!
            if (present(tissue_leaf_acc)) then
               tissue_leaf_acc(1:n) = tissue_leaf_acc(1:n) + tl_int_i(1:n)
               tissue_wood_acc(1:n) = tissue_wood_acc(1:n) + tw_int_i(1:n)
            end if
            if (present(transp_acc)) transp_acc(1:n) = transp_acc(1:n) + tr_int_i(1:n)
            t = t + dt
            nsteps = nsteps + 1_ik
            err_prev = err
            if (.not. clamped) dt_warm = dt
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         !----- Work cap: stop on a step the dispatcher will redo on ARK. The partial-time y_out is    !
         !      never committed (stiff_bail). ---------------------------------------------------------!
         if (nsteps + nrej >= RK45_WORK_CAP) exit
      end do
      call state_init(y, n, nsl, y_out)
      if (present(dt_warm_out)) dt_warm_out = dt_warm
   end subroutine adaptive_rk45_march

   !=======================================================================================!
   !  column_fast_step_rk45 -- the INTEG_RK45 fast step: build the frozen record as ARK does,     !
   !  march one dt_fast, commit, and close the whole-column water and energy ledgers.            !
   !=======================================================================================!
   subroutine column_fast_step_rk45(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget,  &
                                    ws, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
                                    stiff_bail, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(column_config_t),   intent(in)    :: col_config
      type(column_site_t),     intent(in)    :: col_site
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      type(fast_workspace_t),  intent(inout) :: ws       !< this thread's reusable step storage
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters
      !----- True when the march hit RK45_WORK_CAP without resolving the step: biophys is left at     !
      !      state^n, and the caller redoes this dt_fast on ARK. --------------------------------------!
      logical,     optional,   intent(out)   :: stiff_bail
      real(wp),    optional,   intent(inout) :: cdiag(:,:)   !< (N_CO_FAST_ROWS, ncoh) per-cohort diagnostic capture

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out
      real(wp)    :: dt0, cas_mass_capacity, enth0, shv0, enth1, shv1
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1, w_surface0
      real(wp)    :: w_out_acc, e_in_acc, e_out_acc, w_in, w_out, e_in, e_out
      real(wp)    :: dt_warm_next, cond_dep, cond_dep_enth, face_acc
      real(wp)    :: w_pond_rk, runoff_rk
      real(wp)    :: e_pond0, e_pond_rk, over_enth_rk   ! [J/m2] pond enthalpy (start, end); runoff enthalpy
      !----- Canopy-film ledger scratch. --------------------------------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, intercept_total
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: tl_int_acc(col_cohort%n), tw_int_acc(col_cohort%n), tr_int_acc(col_cohort%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, nsteps, nrej

      n = col_cohort%n
      nsl = col_site%soil_profile%n_active
      !----- The stage-clamp count accumulates down the call chain, so this step's tally starts here. !
      budget%clamp_stage_n = 0_ik
      budget%theta_ood_max = 0.0_wp
      !----- State^n pond, taken before the commit below overwrites it. -----------------------------!
      w_surface0 = biophys%soil_water%w_surface
      e_pond0    = biophys%soil_water%w_surface_energy
      if (present(stiff_bail)) stiff_bail = .false.

      cond_dep = 0.0_wp
      cond_dep_enth = 0.0_wp
      face_acc = 0.0_wp
      budget%atm_heat_export = 0.0_wp
      budget%atm_vap_export = 0.0_wp
      call build_column_frozen(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws%spac, &
                               n, nsl, frozen, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)

      dt0 = dt_fast
      if (biophys%adapt_dt_last > tiny_num) dt0 = min(biophys%adapt_dt_last, dt_fast)
      ec = col_config%integrator%error_control
      ec%p_order = RK45_P_ORDER
      tl_int_acc(1:n) = 0.0_wp
      tw_int_acc(1:n) = 0.0_wp
      tr_int_acc(1:n) = 0.0_wp
      call adaptive_rk45_march(y, frozen, col_site, n, nsl, dt_fast, ec, dt0, y_out, nsteps, nrej,             &
                              w_out_acc, e_in_acc, e_out_acc, dt_warm_out=dt_warm_next,           &
                              clamp_stage_n=budget%clamp_stage_n,                                  &
                              cond_acc=cond_dep, cond_enth_acc=cond_dep_enth,                     &
                              atm_heat_acc=budget%atm_heat_export, atm_vap_acc=budget%atm_vap_export, &
                              ood_max=budget%theta_ood_max, face_resid_acc=face_acc,                 &
                              tissue_leaf_acc=tl_int_acc, tissue_wood_acc=tw_int_acc, transp_acc=tr_int_acc)
      biophys%adapt_dt_last = dt_warm_next
      budget%integ_nsteps = nsteps
      budget%integ_nrej = nrej
      !----- Per-layer face closure: the committed theta against the faces and uptake it was given. --!
      budget%soil_face_mass%resid   = face_acc
      budget%soil_face_mass%worst   = max(budget%soil_face_mass%worst, abs(face_acc))
      budget%soil_face_mass%abs_sum = budget%soil_face_mass%abs_sum + abs(face_acc)
      budget%soil_face_mass%n_check = budget%soil_face_mass%n_check + 1_ik

      !----- The march gave up (stiff regime): return before committing anything, since y_out is a    !
      !      partial-time state, so the caller redoes this dt_fast on ARK from state^n. ---------------!
      if (nsteps + nrej >= RK45_WORK_CAP) then
         if (present(stiff_bail)) stiff_bail = .true.
         if (present(converged))  converged  = .false.
         if (present(iters))      iters      = nsteps
         return
      end if

      !----- Canopy films: capacity clamp and overflow/deficit at the store's boundary after the march, !
      !      as on ARK: a film can overdraw below 0, so it is floored and the deficit booked as        !
      !      negative outflow. ----------------------------------------------------------------------!
      surf_overflow = 0.0_wp
      surf_deficit = 0.0_wp
      if (col_config%canopy_water_on) then
         call clamp_canopy_film(y_out, col_cohort%lai, col_cohort%wai,                                  &
                                col_config%leaf_surf_water_max(col_cohort%pft(1:n)),                    &
                                col_config%wood_surf_water_max(col_cohort%pft(1:n)), n, surf_overflow, surf_deficit)
      end if

      !----- Unpack into biophys. -----------------------------------------------------------------!
      call unpack_column_state(y_out, n, nsl, biophys)
      !----- Deposit the condensate into soil layer 1, as on ARK (docs/science/numerical_scheme.md §7, !
      !      item 4), at the b-weighted stage enthalpy the canopy air was debited. --------------------!
      if (cond_dep > 0.0_wp) then
         call deposit_condensate(y_out, col_site%soil_profile%dz(1), cond_dep, cond_dep_enth)
         biophys%soil_energy%soil_energy(1) = y_out%soil_energy(1)
         biophys%soil_water%theta(1) = y_out%theta(1)
      end if
      call diagnose_soil_temps(y_out, col_site%soil_profile%dry_heat_capacity, nsl, biophys%soil_energy%soil_temp, &
            biophys%soil_energy%soil_fliq)
      !----- The pond, its runoff and their enthalpy are the water network's, as theta is            !
      !      (docs/science/soil_plant_water.md). ---------------------------------------------------!
      w_pond_rk    = frozen%hydrology%w_surface1
      e_pond_rk    = frozen%hydrology%w_surface_energy1
      runoff_rk    = frozen%hydrology%runoff_surf * dt_fast
      over_enth_rk = frozen%hydrology%runoff_enth * dt_fast
      biophys%soil_water%w_surface      = w_pond_rk
      biophys%soil_water%w_surface_energy = e_pond_rk

      !----- Commit the tissues' energy from their time-average temperature, as on ARK. -------------!
      call commit_tissue_energy(frozen%tissue, col_cohort%nplant, dt_fast, tl_int_acc(1:n) / dt_fast,         &
                                tw_int_acc(1:n) / dt_fast, n, biophys, tissue_store0, tissue_store1)
      !----- The transpiration each cohort's leaves lost over the step, as on the ARK path. ---------!
      if (present(cdiag)) cdiag(CO_TRANSP, 1:n) = transp_per_leaf_area(tr_int_acc(1:n), dt_fast,          &
                                                                       col_cohort%nplant(1:n), col_cohort%leaf_area(1:n))

      !----- Whole-column water and energy ledgers only: every store rides the same right-hand side,  !
      !      so these are the per-store ledgers here. --------------------------------------------------!
      cas_mass_capacity = frozen%cas%cas_mass_capacity
      enth0 = y%cas_enthalpy
      shv0 = y%cas_shv
      enth1 = y_out%cas_enthalpy
      shv1 = y_out%cas_shv
      e_soil0 = soil_energy_store(y%soil_energy,     col_site%soil_profile%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, col_site%soil_profile%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     col_site%soil_profile%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, col_site%soil_profile%dz, nsl)
      w_plant0 = plant_water_store(col_cohort%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(col_cohort%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- The films, valued as on ARK. -------------------------------------------------------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))

      !----- In: rainfall, snowfall and shed water; for energy, the stage inputs (e_in_acc), the      !
      !      intercepted rain at the film's valuation, the rainfall enthalpy the pack took, and what   !
      !      the pond's inflow brings, shed water included (ground_inflow_enth). The _acc terms and    !
      !      surf_overflow/surf_deficit are amounts; the rates take dt_fast. -------------------------!
      w_in  = (forc%rainfall + forc%snowfall + biophys%shed_water_rate) * dt_fast
      !----- Out: the stages' vapour and drainage, the film overflow net of deficit, the pond's runoff. !
      w_out = w_out_acc + surf_overflow - surf_deficit + runoff_rk
      e_in  = e_in_acc + intercept_total * dt_fast * internal_energy_liquid(frozen%hydrology%t_film_valuation)        &
              + frozen%snow%acc_enth                                                             &
              + frozen%hydrology%ground_inflow_enth * dt_fast
      !----- Runoff leaves at the pond's own temperature (docs/science/numerical_scheme.md §5b). -----!
      e_out = e_out_acc + (surf_overflow - surf_deficit) * internal_energy_liquid(frozen%hydrology%t_film_valuation)     &
              + over_enth_rk

      !----- These checks never stop the run: column_fast_step may still discard the step for the    !
      !      ARK rescue, and calls rk45_ledgers_stop only for a step it keeps. -----------------------!
      call budget_check(budget%whole_water,                                                            &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0                     &
                        + frozen%snow%swe0,                                                          &
                        w_soil1 + cas_mass_capacity*shv1 + w_pond_rk + w_plant1 + surf_water1                        &
                        + frozen%snow%swe1,                                                          &
                        w_in, w_out, dt_fast, BUDGET_WATER_RATE_FLOOR, 'whole_water (rk45)', .false.)
      !----- The snow store and its accumulated rainfall enthalpy join the ledger; 0 without snow. ----!
      call budget_check(budget%whole_energy,                                                           &
                        !----- Meltwater passes pack -> pond, both stores, so it telescopes. ----------!
                        e_soil0                           + cas_mass_capacity*enth0 + surf_enth0                    &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                              &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond_rk            &
                        + tissue_store1,                                                                &
                        e_in, e_out, dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'whole_energy (rk45)', .false.)
      !----- No canopy-air CO2 ledger here: rk45_column_step does not accumulate the stages' CO2     !
      !      exchange with the atmosphere. ----------------------------------------------------------!

      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_rk45

   !----- Under [fast].debug_error, stop the run if a step column_fast_step keeps breached a        !
   !      whole-column ledger; a step handed to the ARK rescue is discarded and must not stop it. -!
   subroutine rk45_ledgers_stop(budget, dt_fast, col_config, col_site)
      type(column_budget_t), intent(in) :: budget
      real(wp),              intent(in) :: dt_fast
      type(column_config_t), intent(in) :: col_config
      type(column_site_t),   intent(in) :: col_site
      if (.not. col_site%debug_error) return
      if (.not. last_check_closed(budget%whole_water, dt_fast, BUDGET_WATER_RATE_FLOOR)) then
         write(*, '(a,es13.5)') 'meds_fast_rk45: whole_water (rk45) did not close, resid = ',               &
                                budget%whole_water%resid
         error stop 'meds_fast_rk45: a kept RK45 step breached its whole-column water ledger'
      end if
      if (.not. last_check_closed(budget%whole_energy, dt_fast, BUDGET_ENERGY_RATE_FLOOR)) then
         write(*, '(a,es13.5)') 'meds_fast_rk45: whole_energy (rk45) did not close, resid = ',              &
                                budget%whole_energy%resid
         error stop 'meds_fast_rk45: a kept RK45 step breached its whole-column energy ledger'
      end if
   end subroutine rk45_ledgers_stop

   !----- Did the explicit march commit a railed state? Its 5th- and 4th-order solutions can rail  !
   !      together, so the controller sees err ~ 0 and accepts an impossible temperature; this is  !
   !      the dispatcher's second rescue trigger. -------------------------------------------------!
   pure function rk45_state_railed(biophys, nsl) result(railed)
      type(patch_biophys_t), intent(in) :: biophys
      integer(ik),           intent(in) :: nsl
      logical :: railed
      real(wp), parameter :: T_LO = 185.0_wp, T_HI = 345.0_wp
      railed =      biophys%cas%can_temp <= T_LO .or. biophys%cas%can_temp >= T_HI                          &
               .or. any(biophys%soil_energy%soil_temp(1:nsl) <= T_LO)                                        &
               .or. any(biophys%soil_energy%soil_temp(1:nsl) >= T_HI)
   end function rk45_state_railed

end module meds_fast_rk45
