!==========================================================================================!
! meds_leaf_solver -- solve the coupled assimilation / stomatal-conductance / intercellular- !
! CO2 (A-gs-Ci) system for one leaf.                                                        !
!                                                                                          !
! The biochemistry is temperature-scaled ONCE (Vcmax, Jmax, Rd, TPU, Kc, Ko, Gamma*, the     !
! electron-transport rate J or the C4 light slope), water stress is applied, and then a       !
! single intercellular-CO2 Ci is found by a bracketed bisection on [Gamma*+eps, Ca]. ONE       !
! residual covers all stomatal models and both pathways:                                       !
!   * Leuning / Medlyn -- a diffusion residual  f(Ci) = Ci - (Cs - 1.6 A/gs(A)),             !
!   * Katul            -- the optimality residual  A'(Cs-Ci)^2 - 1.6 D lambda (A'(Cs-Ci)+A), !
!     so Katul sets Ci directly (it has no gs(A) law) and gs is back-computed.                !
! A respiring/night leaf (no positive-assimilation root, A(Ca) <= 0) is handled by a closed-   !
! stomata branch (gs = g0) before bracketing. Working concentration unit is the mole fraction   !
! [umol/mol]; the Pa Michaelis constants are converted with the air pressure.                   !
!==========================================================================================!
module meds_leaf_solver
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : p_std, tiny_num
   use meds_config,             only : SM_LEUNING, SM_MEDLYN, SM_KATUL
   use meds_leaf_types,         only : leaf_env_t, leaf_flux_t, leaf_photo_params_t,           &
                                       PATH_C3, PATH_C4, LIM_NONE, LIM_RUBISCO, LIM_RUBP,      &
                                       LIM_PRODUCT, LIM_C4_PEP
   use meds_leaf_temp_response, only : temp_response, arrhenius_scale
   use meds_leaf_photosynthesis,only : assim_demand_c3, assim_demand_c4, electron_transport_j
   use meds_leaf_stomata,       only : stomata_gs_leuning, stomata_gs_medlyn, katul_lambda
   implicit none
   private

   public :: solve_leaf_gas_exchange

   real(wp),    parameter :: ci_tol_ppm = 1.0e-3_wp    !< [umol/mol] Ci convergence tolerance (~1e-4 Pa)
   real(wp),    parameter :: lo_eps_ppm = 1.0e-3_wp    !< [umol/mol] offset of the lower bracket above Gamma*
   integer(ik), parameter :: max_iter   = 100_ik       !< bisection iteration cap (safety net)

contains

   !---------------------------------------------------------------------------------------!
   ! Solve the leaf A-gs-Ci system. The selectors (sm, tresp, colim, use_bl) come from the  !
   ! run config; p carries every per-PFT and shared parameter so this routine is self-       !
   ! contained and unit-testable without a full meds_config_t.                              !
   !---------------------------------------------------------------------------------------!
   subroutine solve_leaf_gas_exchange(env, p, sm, tresp, colim, use_bl, flux)
      type(leaf_env_t),          intent(in)  :: env
      type(leaf_photo_params_t), intent(in)  :: p
      integer(ik),               intent(in)  :: sm, tresp, colim
      logical,                   intent(in)  :: use_bl
      type(leaf_flux_t),         intent(out) :: flux

      real(wp) :: t_leaf, pressure, ca_ppm, o2_ppm, ddef, beta, stress
      real(wp) :: vcmax, jmax, jrate, tpu, rd, kc_ppm, ko_ppm, gstar_ppm
      real(wp) :: aj_light, kp_eff, lambda_eff
      real(wp) :: lo, hi, mid, flo, fhi, fmid, ci_sol, an_open
      real(wp) :: a_gross, ac, aj, ap, an, cs_sol, gs_sol
      integer(ik) :: it
      logical     :: converged, do_bl, force_g0

      t_leaf   = env%leaf_temp
      pressure = env%pressure
      ca_ppm   = env%ca
      o2_ppm   = p%o2_mol_frac * 1.0e6_wp
      ddef     = env%vpd / pressure                       ! mole-fraction water deficit D
      do_bl    = use_bl .and. env%gb > 0.0_wp

      !----- Temperature-scale the biochemistry (Kc/Ko/Gamma* always Arrhenius; Pa -> ppm). -!
      kc_ppm    = arrhenius_scale(p%kc25,    p%ea_kc,    t_leaf) / pressure * 1.0e6_wp
      ko_ppm    = arrhenius_scale(p%ko25,    p%ea_ko,    t_leaf) / pressure * 1.0e6_wp
      gstar_ppm = arrhenius_scale(p%gstar25, p%ea_gstar, t_leaf) / pressure * 1.0e6_wp
      vcmax = temp_response(tresp, p%vcmax25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      jmax  = temp_response(tresp, p%jmax25,  p%ea_jmax,  p%hd_jmax,  p%ds_jmax,  t_leaf)
      tpu   = temp_response(tresp, p%tpu25,   p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      rd    = temp_response(tresp, p%rd25,    p%ea_rd,    p%hd_rd,    p%ds_rd,    t_leaf)

      !----- Water stress beta(psi_leaf): downregulate capacity (Leuning/Medlyn) or lambda  !
      !       (Katul). -------------------------------------------------------------------!
      beta = (env%psi_leaf - p%psi_close) / (p%psi_open - p%psi_close)
      beta = min(max(beta, 0.0_wp), 1.0_wp)
      lambda_eff = katul_lambda(p%lambda25, beta, p%lambda_psi_exp)
      if (sm /= SM_KATUL) then
         stress = beta ; vcmax = vcmax * stress ; jmax = jmax * stress ; tpu = tpu * stress
      end if

      !----- Light: C3 non-rectangular hyperbola J; C4 linear light-limited slope. --------!
      if (p%pathway == PATH_C4) then
         gstar_ppm = 0.0_wp
         aj_light  = p%quantum_yield * p%absorptance * env%par
         kp_eff    = p%kp25 * pressure / p_std
         jrate     = 0.0_wp
      else
         jrate     = electron_transport_j(env%par, p%absorptance, p%phi_psii, jmax, p%theta_j)
         aj_light  = 0.0_wp ; kp_eff = 0.0_wp
      end if

      !----- Closed/night branch: no positive-assimilation root (best-case net <= 0). ------!
      an_open = net_at(ca_ppm)
      if (an_open <= 0.0_wp) then
         gs_sol = p%g0
         an     = an_open
         cs_sol = ca_ppm ; if (do_bl) cs_sol = ca_ppm - an * 1.4_wp / env%gb
         ci_sol = cs_sol - 1.6_wp * an / max(gs_sol, tiny_num)
         call fill_flux(an + rd, an, gs_sol, ci_sol, cs_sol, rd, LIM_NONE, .true.)
         return
      end if

      !----- Bracket Ci in (Gamma*, Ca] and bisect the residual to ci_tol_ppm. If the chosen  !
      !       stomatal model yields no consistent open solution (no sign change -- e.g. Katul  !
      !       under strong water stress where lambda -> large), fall back to a g0-pinned       !
      !       (closed-stomata) diffusion solve, which always brackets when net A(Ca) > 0. ----!
      force_g0 = .false.
      lo = gstar_ppm + lo_eps_ppm ; hi = ca_ppm
      flo = residual(lo) ; fhi = residual(hi)
      if (flo * fhi > 0.0_wp) then
         force_g0 = .true.
         lo = gstar_ppm + lo_eps_ppm ; hi = ca_ppm
         flo = residual(lo) ; fhi = residual(hi)
      end if
      converged = .false. ; ci_sol = 0.5_wp * (lo + hi)
      if (flo * fhi <= 0.0_wp) then
         do it = 1_ik, max_iter
            mid  = 0.5_wp * (lo + hi)
            fmid = residual(mid)
            if (flo * fmid <= 0.0_wp) then ; hi = mid ; fhi = fmid ; else ; lo = mid ; flo = fmid ; end if
            if (hi - lo < ci_tol_ppm) exit
         end do
         ci_sol = 0.5_wp * (lo + hi) ; converged = (hi - lo < ci_tol_ppm)
      end if

      !----- Assemble the solution: net A, surface CO2, back-computed gs, transpiration. ---!
      call eval_demand(ci_sol, a_gross, ac, aj, ap) ; an = a_gross - rd
      cs_sol = ca_ppm ; if (do_bl) cs_sol = ca_ppm - an * 1.4_wp / env%gb
      !----- Back-compute gs from the diffusion identity; if the boundary layer pushed Cs at  !
      !       or below Ci (degenerate), pin gs to g0 and report Cs as the surface CO2. -------!
      if (cs_sol - ci_sol > tiny_num) then
         gs_sol = max(1.6_wp * an / (cs_sol - ci_sol), p%g0)
      else
         gs_sol = p%g0 ; ci_sol = cs_sol
      end if
      call fill_flux(a_gross, an, gs_sol, ci_sol, cs_sol, rd, pick_limit(ac, aj, ap, an), converged)

   contains

      !----- Gross + raw limitation rates at a trial Ci (dispatch on pathway). -----------!
      pure subroutine eval_demand(ci, ag, rac, raj, rap)
         real(wp), intent(in)  :: ci
         real(wp), intent(out) :: ag, rac, raj, rap
         if (p%pathway == PATH_C4) then
            call assim_demand_c4(ci, vcmax, aj_light, kp_eff, colim, p%theta_cj, p%theta_ic,   &
                                 ag, rac, raj, rap)
         else
            call assim_demand_c3(ci, vcmax, jrate, tpu, gstar_ppm, kc_ppm, ko_ppm, o2_ppm,     &
                                 colim, p%theta_j, ag, rac, raj, rap)
         end if
      end subroutine eval_demand

      !----- Net assimilation at a trial Ci. ---------------------------------------------!
      pure function net_at(ci) result(an_loc)
         real(wp), intent(in) :: ci
         real(wp)             :: an_loc, ag, rac, raj, rap
         call eval_demand(ci, ag, rac, raj, rap) ; an_loc = ag - rd
      end function net_at

      !----- Root residual: diffusion (Leuning/Medlyn) or Katul optimality. --------------!
      function residual(ci) result(r)
         real(wp), intent(in) :: ci
         real(wp)             :: r, an_loc, csl, gsl, ci_pred, anp, dh
         an_loc = net_at(ci)
         csl    = ca_ppm ; if (do_bl) csl = ca_ppm - an_loc * 1.4_wp / env%gb
         if (force_g0) then                              ! closed-stomata fallback: gs pinned to g0
            ci_pred = csl - 1.6_wp * an_loc / max(p%g0, tiny_num)
            r       = ci - ci_pred ; return
         end if
         if (sm == SM_KATUL) then
            dh  = max(1.0e-3_wp * abs(ci), 1.0e-2_wp)
            anp = (net_at(ci + dh) - net_at(ci - dh)) / (2.0_wp * dh)
            r   = anp * (csl - ci)**2 - 1.6_wp * ddef * lambda_eff * (anp * (csl - ci) + an_loc)
         else
            if (sm == SM_LEUNING) then
               gsl = stomata_gs_leuning(an_loc, csl, gstar_ppm, env%vpd, p%g0, p%g1, p%d0)
            else
               gsl = stomata_gs_medlyn(an_loc, csl, env%vpd, p%g0, p%g1)
            end if
            ci_pred = csl - 1.6_wp * an_loc / max(gsl, tiny_num)
            r       = ci - ci_pred
         end if
      end function residual

      !----- Map the binding gross rate to a limitation flag. ----------------------------!
      pure function pick_limit(rac, raj, rap, an_loc) result(lim)
         real(wp), intent(in) :: rac, raj, rap, an_loc
         integer(ik)          :: lim
         if (an_loc <= 0.0_wp) then
            lim = LIM_NONE
         else if (rac <= raj .and. rac <= rap) then
            lim = LIM_RUBISCO
         else if (raj <= rac .and. raj <= rap) then
            lim = LIM_RUBP
         else if (p%pathway == PATH_C4) then
            lim = LIM_C4_PEP
         else
            lim = LIM_PRODUCT
         end if
      end function pick_limit

      !----- Pack the output flux record. ------------------------------------------------!
      subroutine fill_flux(ag, an_loc, gs, ci, cs, rd_loc, lim, conv)
         real(wp),    intent(in) :: ag, an_loc, gs, ci, cs, rd_loc
         integer(ik), intent(in) :: lim
         logical,     intent(in) :: conv
         flux%a_gross = ag ; flux%a_net = an_loc ; flux%gs = gs ; flux%ci = ci ; flux%cs = cs
         flux%transpiration = gs * env%vpd / pressure
         flux%rd = rd_loc ; flux%limitation = lim ; flux%converged = conv
      end subroutine fill_flux

   end subroutine solve_leaf_gas_exchange

end module meds_leaf_solver
