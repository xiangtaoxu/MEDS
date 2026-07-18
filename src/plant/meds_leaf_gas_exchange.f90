!==========================================================================================!
! meds_leaf_gas_exchange -- leaf-level gas-exchange COMPUTE kernels (photosynthesis + stomata !
! + coupled Ci solver), merged into one module. FvCB C3 / Collatz C4 demand, the electron-     !
! transport hyperbola, the Leuning / Medlyn / Katul stomatal models, and the bracketed Ci      !
! root-find (solve_leaf_gas_exchange). The public seam leaf_gas_exchange lives in               !
! meds_plant_interface; the raw kernels here are also called directly by meds_plant_capi.        !
!==========================================================================================!
module meds_leaf_gas_exchange
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : p_std, tiny_num, gsw_2_gsc, gbw_2_gbc, mol_2_umol
   use meds_config,        only : COLIM_MIN, COLIM_QUADRATIC, SM_LEUNING, SM_MEDLYN, SM_KATUL
   use meds_plant_types,   only : leaf_env_t, leaf_flux_t, leaf_photo_params_t,                 &
                                  PATH_C3, PATH_C4, LIM_NONE, LIM_RUBISCO, LIM_RUBP,            &
                                  LIM_PRODUCT, LIM_C4_PEP
   use meds_temp_response, only : temp_response, arrhenius_scale
   use meds_numerics,      only : quadratic_smaller_root, bisect_root
   implicit none
   private

   !----- from meds_leaf_photosynthesis.f90 ----------------------------------------------!

   public :: assimilation_demand_c3, assimilation_demand_c4, electron_transport_j

   !----- from meds_leaf_stomata.f90 -----------------------------------------------------!

   public :: stomata_gs_leuning, stomata_gs_medlyn, katul_lambda

   real(wp), parameter :: vpd_floor_pa = 50.0_wp     !< [Pa] VPD floor (avoid 1/sqrt(0) in Medlyn)
   real(wp), parameter :: beta_floor   = 1.0e-4_wp   !< [--] water-stress floor (bound lambda as beta->0)

   !----- from meds_leaf_solver.f90 ------------------------------------------------------!

   public :: solve_leaf_gas_exchange

   real(wp),    parameter :: ci_tol_ppm = 1.0e-3_wp    !< [umol/mol] Ci convergence tolerance (~1e-4 Pa)
   real(wp),    parameter :: lo_eps_ppm = 1.0e-3_wp    !< [umol/mol] offset of the lower bracket above Gamma*
   integer(ik), parameter :: max_iter   = 100_ik       !< bisection iteration cap (safety net)


contains

   !========== meds_leaf_photosynthesis.f90 =============================================!

   !---------------------------------------------------------------------------------------!
   ! Actual electron transport rate J from absorbed PAR via the non-rectangular hyperbola   !
   ! theta J^2 - (I2 + Jmax) J + I2 Jmax = 0 (smaller root). I2 = 0.5 phi_psii absorptance PAR.!
   !---------------------------------------------------------------------------------------!
   elemental pure function electron_transport_j(par, absorptance, phi_psii, jmax, theta) result(j)
      real(wp), intent(in) :: par         !< [umol photon/m2/s] incident PAR
      real(wp), intent(in) :: absorptance !< [--] leaf PAR absorptance
      real(wp), intent(in) :: phi_psii    !< [--] PSII quantum yield (electrons/photon)
      real(wp), intent(in) :: jmax        !< [umol/m2/s] electron-transport capacity (T-scaled)
      real(wp), intent(in) :: theta       !< [--] curvature (0 < theta < 1)
      real(wp)             :: j, i2
      i2 = 0.5_wp * phi_psii * absorptance * par
      j  = quadratic_smaller_root(theta, i2, jmax)
   end function electron_transport_j

   !---------------------------------------------------------------------------------------!
   ! C3 demand: gross assimilation A_gross and the three raw limitation rates (Ac/Aj/Ap).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colim, theta,        &
                                   A_gross, Ac, Aj, Ap)
      real(wp),    intent(in)  :: ci, vcmax, j, tpu, gstar, kc, ko, o2, theta
      integer(ik), intent(in)  :: colim
      real(wp),    intent(out) :: A_gross, Ac, Aj, Ap
      Ac = vcmax * (ci - gstar) / (ci + kc * (1.0_wp + o2 / ko))
      Aj = j     * (ci - gstar) / (4.0_wp * ci + 8.0_wp * gstar)
      Ap = 3.0_wp * tpu
      A_gross = combine_limits(Ac, Aj, Ap, colim, theta)
   end subroutine assimilation_demand_c3

   !---------------------------------------------------------------------------------------!
   ! C4 demand (Collatz 1992): Ac = Vcmax, Aj = light-limited slope (passed in), Ap = PEPcase !
   ! CO2 limitation kp_eff*Ci. Gamma* ~ 0 (CO2-concentrating mechanism suppresses photoresp). !
   !---------------------------------------------------------------------------------------!
   pure subroutine assimilation_demand_c4(ci, vcmax, Aj_light, kp_eff, colim, theta_cj, theta_ic,    &
                                   A_gross, Ac, Aj, Ap)
      real(wp),    intent(in)  :: ci, vcmax, Aj_light, kp_eff, theta_cj, theta_ic
      integer(ik), intent(in)  :: colim
      real(wp),    intent(out) :: A_gross, Ac, Aj, Ap
      real(wp) :: Ai
      Ac = vcmax
      Aj = Aj_light
      Ap = kp_eff * ci
      if (colim == COLIM_QUADRATIC) then
         Ai      = quadratic_smaller_root(theta_cj, Ac, Aj)   ! co-limit Rubisco & light
         A_gross = quadratic_smaller_root(theta_ic, Ai, Ap)   ! co-limit with PEPcase CO2
      else
         A_gross = min(Ac, Aj, Ap)
      end if
   end subroutine assimilation_demand_c4

   !---------------------------------------------------------------------------------------!
   ! Combine the three C3 limitation rates: sharp min, or two nested smoothing quadratics.  !
   ! (The smaller root of the co-limitation quadratic is the shared                         !
   ! meds_numerics%quadratic_smaller_root.)                                                 !
   !---------------------------------------------------------------------------------------!
   pure function combine_limits(Ac, Aj, Ap, colim, theta) result(A)
      real(wp),    intent(in) :: Ac, Aj, Ap, theta
      integer(ik), intent(in) :: colim
      real(wp)                :: A, Ai
      if (colim == COLIM_QUADRATIC) then
         Ai = quadratic_smaller_root(theta, Ac, Aj)      ! co-limit Rubisco & RuBP
         A  = quadratic_smaller_root(theta, Ai, Ap)      ! co-limit with product (TPU)
      else
         A = min(Ac, Aj, Ap)
      end if
   end function combine_limits


   !========== meds_leaf_stomata.f90 ====================================================!

   !---------------------------------------------------------------------------------------!
   ! Leuning (1995) semi-empirical stomatal conductance.                                   !
   !---------------------------------------------------------------------------------------!
   pure function stomata_gs_leuning(a_net, cs, gstar, vpd, g0, g1, d0) result(gs)
      real(wp), intent(in) :: a_net   !< [umol/m2/s] net assimilation
      real(wp), intent(in) :: cs      !< [umol/mol]  leaf-surface CO2
      real(wp), intent(in) :: gstar   !< [umol/mol]  CO2 compensation point
      real(wp), intent(in) :: vpd     !< [Pa]        leaf-to-air VPD
      real(wp), intent(in) :: g0, g1  !< [mol/m2/s], [--]
      real(wp), intent(in) :: d0      !< [Pa]        humidity sensitivity
      real(wp)             :: gs, denom
      if (a_net <= 0.0_wp) then
         gs = g0
         return
      end if
      denom = max(cs - gstar, tiny_num) * (1.0_wp + vpd / d0)
      gs = g0 + g1 * a_net / denom
   end function stomata_gs_leuning

   !---------------------------------------------------------------------------------------!
   ! Medlyn et al. (2011) unified stomatal optimization (USO) conductance.                 !
   !---------------------------------------------------------------------------------------!
   pure function stomata_gs_medlyn(a_net, cs, vpd, g0, g1) result(gs)
      real(wp), intent(in) :: a_net   !< [umol/m2/s] net assimilation
      real(wp), intent(in) :: cs      !< [umol/mol]  leaf-surface CO2
      real(wp), intent(in) :: vpd     !< [Pa]        leaf-to-air VPD
      real(wp), intent(in) :: g0, g1  !< [mol/m2/s], [kPa^0.5]
      real(wp)             :: gs, vpd_kpa
      if (a_net <= 0.0_wp) then
         gs = g0
         return
      end if
      vpd_kpa = max(vpd, vpd_floor_pa) * 1.0e-3_wp
      gs = g0 + gsw_2_gsc * (1.0_wp + g1 / sqrt(vpd_kpa)) * a_net / max(cs, tiny_num)
   end function stomata_gs_medlyn

   !---------------------------------------------------------------------------------------!
   ! Katul marginal water-use efficiency lambda [umol CO2/mol H2O], scaled by the water-     !
   ! stress factor beta in (0,1]: drier leaves carry a larger marginal water cost, so lambda  !
   ! rises (lambda = lambda25 * beta^(-exp)) and the optimal stomata close. exp = 0 disables.  !
   !---------------------------------------------------------------------------------------!
   pure function katul_lambda(lambda25, beta, lambda_psi_exp) result(lambda)
      real(wp), intent(in) :: lambda25        !< [umol CO2/mol H2O] well-watered marginal WUE
      real(wp), intent(in) :: beta            !< [--] water-stress factor (0 = closed, 1 = open)
      real(wp), intent(in) :: lambda_psi_exp  !< [--] water-stress exponent
      real(wp)             :: lambda
      lambda = lambda25 * max(beta, beta_floor) ** (-lambda_psi_exp)
   end function katul_lambda


   !========== meds_leaf_solver.f90 =====================================================!

   !---------------------------------------------------------------------------------------!
   ! Solve the leaf A-gs-Ci system. The selectors (sm, tresp, colim, use_boundary_layer) come from the  !
   ! run config; p carries every per-PFT and shared parameter so this routine is self-       !
   ! contained and unit-testable without a full meds_config_t.                              !
   !   use_boundary_layer -- account for the leaf boundary layer (env%gb). When .true. (and gb > 0),     !
   !     leaf-surface CO2 is drawn down from ambient (Cs = Ca - 1.4*A/gb) and transpiration  !
   !     puts stomata gs in SERIES with gb; when .false. the leaf is well-coupled (Cs = Ca,  !
   !     E = gs*VPD/pressure). The 1.4 / 1.6 factors are the boundary-layer / stomatal       !
   !     H2O:CO2 diffusivity ratios.                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine solve_leaf_gas_exchange(env, p, sm, tresp, colim, use_boundary_layer, flux)
      type(leaf_env_t),          intent(in)  :: env
      type(leaf_photo_params_t), intent(in)  :: p
      integer(ik),               intent(in)  :: sm, tresp, colim
      logical,                   intent(in)  :: use_boundary_layer
      type(leaf_flux_t),         intent(out) :: flux

      real(wp) :: t_leaf, pressure, ca_ppm, o2_ppm, ddef, beta_nonstomata, beta_stomata, g1_eff
      real(wp) :: vcmax, jmax, jrate, tpu, rd, kc_ppm, ko_ppm, gstar_ppm
      real(wp) :: Aj_light, kp_eff, lambda_eff
      real(wp) :: lo0, hi0, ci_sol, An_open
      real(wp) :: A_gross, Ac, Aj, Ap, An, cs_sol, gs_sol
      logical  :: converged, do_boundary_layer, force_g0

      t_leaf   = env%leaf_temp
      pressure = env%pressure
      ca_ppm   = env%ca
      o2_ppm   = p%o2_mol_frac * mol_2_umol
      ddef     = env%vpd / pressure                       ! mole-fraction water deficit D
      do_boundary_layer    = use_boundary_layer .and. env%gb > 0.0_wp

      !----- Temperature-scale the biochemistry (Kc/Ko/Gamma* always Arrhenius; Pa -> ppm). -!
      kc_ppm    = arrhenius_scale(p%kc25,    p%ea_kc,    t_leaf) / pressure * mol_2_umol
      ko_ppm    = arrhenius_scale(p%ko25,    p%ea_ko,    t_leaf) / pressure * mol_2_umol
      gstar_ppm = arrhenius_scale(p%gstar25, p%ea_gstar, t_leaf) / pressure * mol_2_umol
      vcmax = temp_response(tresp, p%vcmax25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      jmax  = temp_response(tresp, p%jmax25,  p%ea_jmax,  p%hd_jmax,  p%ds_jmax,  t_leaf)
      tpu   = temp_response(tresp, p%tpu25,   p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      rd    = temp_response(tresp, p%rd25,    p%ea_rd,    p%hd_rd,    p%ds_rd,    t_leaf)

      !----- Water stress, split into two independently-tunable limbs (Sabot 2022 / Zhou 2013): !
      !   beta_nonstomata -- capacity limb: a linear psi_LEAF ramp downregulating Vcmax/Jmax/TPU, !
      !     applied to ALL stomatal models (a leaf-biochemistry effect, scheme-independent).      !
      !   beta_stomata    -- stomatal limb: min(1, exp(sref*psi_SOIL)) downregulating the         !
      !     Leuning/Medlyn slope g1 and the Katul marginal WUE lambda (lambda ~                   !
      !     beta_stomata^(-lambda_psi_exp); lambda_psi_exp = 2 recovers Sabot's g1<->lambda).     !
      beta_nonstomata = (env%psi_leaf - p%psi_close) / (p%psi_open - p%psi_close)
      beta_nonstomata = min(max(beta_nonstomata, 0.0_wp), 1.0_wp)
      vcmax = vcmax * beta_nonstomata
      jmax  = jmax  * beta_nonstomata
      tpu   = tpu   * beta_nonstomata
      beta_stomata = min(1.0_wp, exp(p%sref_stomata * env%psi_soil))
      g1_eff       = p%g1 * beta_stomata
      lambda_eff   = katul_lambda(p%lambda25, beta_stomata, p%lambda_psi_exp)

      !----- Light: C3 non-rectangular hyperbola J; C4 linear light-limited slope. --------!
      if (p%pathway == PATH_C4) then
         gstar_ppm = 0.0_wp
         Aj_light  = p%quantum_yield * p%absorptance * env%par
         !----- C4 PEP/CO2-limited slope: temperature-scale kp like every other C4 term (BUG4).   !
         !      ED2 sets kp = klowco2*vm, so kp inherits Vcmax's temperature response (and, under   !
         !      the peaked form, its high-T deactivation) -- reuse the Vcmax Ea/Hd/dS set, exactly  !
         !      as the tpu term above does. temp_response == 1 at 25 degC, so 25 degC is unchanged. !
         kp_eff    = temp_response(tresp, p%kp25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)       &
                     * pressure / p_std
         jrate     = 0.0_wp
      else
         jrate     = electron_transport_j(env%par, p%absorptance, p%phi_psii, jmax, p%theta_j)
         Aj_light  = 0.0_wp
         kp_eff    = 0.0_wp
      end if

      !----- Closed/night branch: no positive-assimilation root (best-case net <= 0). ------!
      An_open = Anet_at_ci(ca_ppm)
      if (An_open <= 0.0_wp) then
         gs_sol = p%g0
         An     = An_open
         cs_sol = ca_ppm
         if (do_boundary_layer) cs_sol = ca_ppm - gbw_2_gbc * An / env%gb
         ci_sol = cs_sol - gsw_2_gsc * An / max(gs_sol, tiny_num)
         call fill_flux(An + rd, An, gs_sol, ci_sol, cs_sol, rd, LIM_NONE, .true.)
         return
      end if

      !----- Bracket Ci in (Gamma*, Ca] and bisect the residual (meds_numerics%bisect_root)   !
      !       to ci_tol_ppm. If the chosen stomatal model yields no consistent open solution   !
      !       (no sign change -- e.g. Katul under strong water stress where lambda -> large),   !
      !       fall back to a g0-pinned (closed-stomata) diffusion solve, which always brackets  !
      !       when net A(Ca) > 0. The explicit-gs provider serves BOTH the Leuning/Medlyn model !
      !       and the g0 fallback (gsl := g0 when force_g0); Katul uses the optimality provider. !
      lo0 = gstar_ppm + lo_eps_ppm
      hi0 = ca_ppm
      force_g0 = .false.
      do                                                 ! at most twice: retry g0-pinned
         if (.not. force_g0 .and. sm == SM_KATUL) then
            call bisect_root(residual_optimality,  lo0, hi0, ci_tol_ppm, max_iter, ci_sol, converged)
         else
            call bisect_root(residual_explicit_gs, lo0, hi0, ci_tol_ppm, max_iter, ci_sol, converged)
         end if
         !----- No sign change on the first (model) attempt -> retry g0-pinned. --------------!
         if (.not. converged .and. .not. force_g0) then
            force_g0 = .true.
            cycle
         end if

         !----- Assemble the solution: net A, surface CO2, back-computed gs, transpiration. ---!
         call eval_assimilation_demand(ci_sol, A_gross, Ac, Aj, Ap)
         An     = A_gross - rd
         cs_sol = ca_ppm
         if (do_boundary_layer) cs_sol = ca_ppm - gbw_2_gbc * An / env%gb
         !----- Katul optimum can land below the cuticular floor g0; re-solve once g0-pinned so   !
         !       A/gs/Ci/E stay mutually consistent (Leuning/Medlyn already return gs >= g0). ----!
         if (sm == SM_KATUL .and. .not. force_g0 .and. cs_sol - ci_sol > tiny_num) then
            if (gsw_2_gsc * An / (cs_sol - ci_sol) < p%g0) then
               force_g0 = .true.
               cycle
            end if
         end if
         exit
      end do
      !----- Back-compute gs from the diffusion identity; if the boundary layer pushed Cs at  !
      !       or below Ci (degenerate), pin gs to g0 and report Cs as the surface CO2. -------!
      if (cs_sol - ci_sol > tiny_num) then
         gs_sol = max(gsw_2_gsc * An / (cs_sol - ci_sol), p%g0)
      else
         gs_sol = p%g0
         ci_sol = cs_sol
      end if
      call fill_flux(A_gross, An, gs_sol, ci_sol, cs_sol, rd, pick_limit(Ac, Aj, Ap, An), converged)

   contains

      !----- Gross + raw limitation rates at a trial Ci (dispatch on pathway). -----------!
      pure subroutine eval_assimilation_demand(ci, Ag, Rac, Raj, Rap)
         real(wp), intent(in)  :: ci
         real(wp), intent(out) :: Ag, Rac, Raj, Rap
         if (p%pathway == PATH_C4) then
            call assimilation_demand_c4(ci, vcmax, Aj_light, kp_eff, colim, p%theta_cj, p%theta_ic,   &
                                 Ag, Rac, Raj, Rap)
         else
            call assimilation_demand_c3(ci, vcmax, jrate, tpu, gstar_ppm, kc_ppm, ko_ppm, o2_ppm,     &
                                 colim, p%theta_j, Ag, Rac, Raj, Rap)
         end if
      end subroutine eval_assimilation_demand

      !----- Net assimilation at a trial Ci. ---------------------------------------------!
      pure function Anet_at_ci(ci) result(An_loc)
         real(wp), intent(in) :: ci
         real(wp)             :: An_loc, Ag, Rac, Raj, Rap
         call eval_assimilation_demand(ci, Ag, Rac, Raj, Rap)
         An_loc = Ag - rd
      end function Anet_at_ci

      !----- Explicit-conductance residual for Leuning / Medlyn (and the g0-pinned fallback):    !
      !       the trial Ci must satisfy the CO2 diffusion identity Ci = Cs - A / gs_co2, where    !
      !       the stomatal model supplies gs. Returns Ci - Ci_predicted, whose ROOT is the        !
      !       solution. PURE so it feeds bisect_root. -----------------------------------------!
      pure function residual_explicit_gs(ci) result(r)
         real(wp), intent(in) :: ci
         real(wp)             :: r, An_loc, cs_surf, gs, ci_pred
         An_loc  = Anet_at_ci(ci)                       ! net assimilation A at this trial Ci
         !----- Leaf-surface CO2: ambient, less the boundary-layer drawdown by the CO2 flux. ---!
         cs_surf = ca_ppm
         if (do_boundary_layer) cs_surf = ca_ppm - gbw_2_gbc * An_loc / env%gb
         !----- Stomatal conductance from the chosen model (or the cuticular floor g0). --------!
         if (force_g0) then
            gs = p%g0                                   ! closed-stomata fallback: gs pinned to g0
         else if (sm == SM_LEUNING) then
            gs = stomata_gs_leuning(An_loc, cs_surf, gstar_ppm, env%vpd, p%g0, g1_eff, p%d0)
         else
            gs = stomata_gs_medlyn(An_loc, cs_surf, env%vpd, p%g0, g1_eff)
         end if
         !----- Ci predicted by CO2 diffusion through the stomata (gs is a WATER conductance, so  !
         !       the CO2 conductance is gs / gsw_2_gsc); residual = trial Ci minus predicted Ci.  !
         ci_pred = cs_surf - gsw_2_gsc * An_loc / max(gs, tiny_num)
         r       = ci - ci_pred
      end function residual_explicit_gs

      !----- Katul optimality residual (no explicit gs). Stomata maximize A - lambda*E; the      !
      !       first-order condition dA/dCi = lambda * dE/dCi, with CO2 supply                     !
      !       A = gs/gsw_2_gsc * (Cs - Ci) and water loss E = gs * D (D = VPD/P = ddef),          !
      !       rearranges to the residual below, whose ROOT is the optimal Ci. PURE (feeds         !
      !       bisect_root). --------------------------------------------------------------------!
      pure function residual_optimality(ci) result(r)
         real(wp), intent(in) :: ci
         real(wp)             :: r, An_loc, cs_surf, dAn_dci, dci
         An_loc  = Anet_at_ci(ci)                       ! net assimilation A at this trial Ci
         !----- Leaf-surface CO2 (ambient less boundary-layer drawdown), as for the gs models. -!
         cs_surf = ca_ppm
         if (do_boundary_layer) cs_surf = ca_ppm - gbw_2_gbc * An_loc / env%gb
         !----- Marginal demand A' = dA/dCi by central difference (A(Ci) is the co-limited FvCB   !
         !       envelope, so the slope is taken numerically; dci is a relative step, abs-floored).!
         dci     = max(1.0e-3_wp * abs(ci), 1.0e-2_wp)
         dAn_dci = (Anet_at_ci(ci + dci) - Anet_at_ci(ci - dci)) / (2.0_wp * dci)
         !----- First-order optimality: A'(Cs-Ci)^2 = gsw_2_gsc * D * lambda * (A'(Cs-Ci) + A). --!
         r = dAn_dci * (cs_surf - ci)**2                                                          &
             - gsw_2_gsc * ddef * lambda_eff * (dAn_dci * (cs_surf - ci) + An_loc)
      end function residual_optimality

      !----- Map the binding gross rate to a limitation flag. ----------------------------!
      pure function pick_limit(Rac, Raj, Rap, An_loc) result(lim)
         real(wp), intent(in) :: Rac, Raj, Rap, An_loc
         integer(ik)          :: lim
         if (An_loc <= 0.0_wp) then
            lim = LIM_NONE
         else if (Rac <= Raj .and. Rac <= Rap) then
            lim = LIM_RUBISCO
         else if (Raj <= Rac .and. Raj <= Rap) then
            lim = LIM_RUBP
         else if (p%pathway == PATH_C4) then
            lim = LIM_C4_PEP
         else
            lim = LIM_PRODUCT
         end if
      end function pick_limit

      !----- Pack the output flux record. ------------------------------------------------!
      subroutine fill_flux(Ag, An_loc, gs, ci, cs, rd_loc, lim, conv)
         real(wp),    intent(in) :: Ag, An_loc, gs, ci, cs, rd_loc
         integer(ik), intent(in) :: lim
         logical,     intent(in) :: conv
         flux%A_gross = Ag
         flux%A_net   = An_loc
         flux%gs      = gs
         flux%ci      = ci
         flux%cs      = cs
         !----- Transpiration uses the TOTAL leaf-to-air water conductance: stomata gs in SERIES    !
         !      with the boundary layer gb (env%gb), consistent with the CO2 solve (which puts       !
         !      1.4*gb in series). Without gb the water flux is overestimated; gated by use_boundary_layer.  ----!
         if (do_boundary_layer) then
            flux%transpiration = gs * env%gb / (gs + env%gb) * env%vpd / pressure
         else
            flux%transpiration = gs * env%vpd / pressure
         end if
         flux%rd         = rd_loc
         flux%limitation = lim
         flux%converged  = conv
      end subroutine fill_flux

   end subroutine solve_leaf_gas_exchange

end module meds_leaf_gas_exchange
