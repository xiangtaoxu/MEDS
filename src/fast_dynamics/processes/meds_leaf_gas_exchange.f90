! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_leaf_gas_exchange -- leaf photosynthesis, stomatal conductance and the coupled Ci   !
! solve, docs/science/leaf_gas_exchange.md §1-§5, eqs. (1)-(10): FvCB C3 / Collatz C4      !
! demand, the Leuning / Medlyn / Katul stomatal models, the two water-stress controls, and !
! solve_leaf_gas_exchange (a bracketed root on Ci). leaf_gas_exchange_batch runs it     !
! over n leaves.                                                                          !
!==========================================================================================!
module meds_leaf_gas_exchange
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : p_std, tiny_num, gsw_2_gsc, gbw_2_gbc, mol_2_umol
   use meds_leaf_params,   only : COLIM_QUADRATIC, SM_LEUNING, SM_KATUL
   use meds_pft_params,    only : PATH_C4
   use meds_pft_params,  only : leaf_photo_pft_t, leaf_photo_table_t
   use meds_therm_lib,     only : temp_response, arrhenius_scale
   use meds_numerics,      only : quadratic_smaller_root, root_bracket_t, regula_falsi_point, illinois_narrow
   implicit none

   !----- The O2 mole fraction and pressure gstar25 is stated at (Bernacchi et al. 2001); Gamma* scales !
   !      with O2 and not with pressure (docs/science/leaf_gas_exchange.md eq. (3)). --------------------!
   real(wp), parameter :: O2_REF_GSTAR = 0.209_wp   !< [mol/mol]
   real(wp), parameter :: P_REF_GSTAR  = 1.0e5_wp   !< [Pa]
   private

   !----- The leaf's environment and fluxes, and the limitation flags. ------------------!
   public :: leaf_env_t, leaf_out_t
   public :: LIM_NONE, LIM_RUBISCO, LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP

   !----- Photosynthetic demand. ---------------------------------------------------------!
   public :: assimilation_demand_c3, assimilation_demand_c4, electron_transport_j

   !----- Stomatal conductance. ----------------------------------------------------------!
   public :: stomata_gs_leuning, stomata_gs_medlyn

   real(wp), parameter :: BETA_FLOOR   = 1.0e-4_wp   !< [--] water-stress floor (bound lambda as beta->0)

   !----- The coupled Ci solver. ---------------------------------------------------------!
   public :: solve_leaf_gas_exchange, rubisco_kinetics
   public :: leaf_gas_exchange_batch, leaf_params_at_capacity

   real(wp),    parameter :: CI_TOL_PPM = 1.0e-3_wp    !< [umol/mol] Ci convergence tolerance (~1e-4 Pa)
   real(wp),    parameter :: LO_EPS_PPM = 1.0e-3_wp    !< [umol/mol] offset of the lower bracket above Gamma*
   integer(ik), parameter :: MAX_ITER   = 100_ik       !< root iteration cap (safety net)

   !----- Which stomatal conductance a pass of the Ci solve uses (solve_leaf_gas_exchange). -----!
   integer(ik), parameter :: GS_FROM_MODEL = 1_ik      !< the configured stomatal model
   integer(ik), parameter :: GS_CUTICULAR  = 2_ik      !< the cuticular floor g0 (the fallback)

   !----- Limitation-regime flag of the binding term at the solution. ---------------------!
   integer(ik), parameter :: LIM_NONE    = 0_ik     !< degenerate (night / net <= 0)
   integer(ik), parameter :: LIM_RUBISCO = 1_ik     !< Rubisco-limited (C3 Ac) / Vcmax-limited (C4)
   integer(ik), parameter :: LIM_RUBP    = 2_ik     !< RuBP / light-limited (Aj)
   integer(ik), parameter :: LIM_PRODUCT = 3_ik     !< triose-phosphate-use limited (C3 Ap)
   integer(ik), parameter :: LIM_C4_PEP  = 4_ik     !< C4 PEPcase CO2 limitation (Ap)

   !----- Leaf-level environmental drivers (no canopy RT / energy balance / hydraulics). ---!
   type :: leaf_env_t
      real(wp) :: par        !< [umol photon/m2/s] incident PAR (leaf absorptance applied internally)
      real(wp) :: leaf_temp  !< [K]    leaf temperature
      real(wp) :: vpd        !< [Pa]   leaf-to-air vapour-pressure deficit
      real(wp) :: ca         !< [umol/mol] reference (canopy-air) CO2 mole fraction
      real(wp) :: pressure   !< [Pa]   air pressure
      real(wp) :: psi_leaf   !< [MPa]  current leaf water potential (<= 0); drives beta_nonstomata (fast)
      real(wp) :: gb         !< [mol H2O/m2/s] boundary-layer conductance (<= 0 => skip, Cs = Ca)
      !----- The supply-side potential driving the slow stomatal control (default 0, well-watered): the !
      !      fast loop passes the cohort's predawn (previous-day maximum) leaf potential, dmax_psi_leaf. !
      !      Named neutrally: it equals the soil's only in wet soil (leaf_gas_exchange.md §4.2). --------!
      real(wp) :: psi = 0.0_wp  !< [MPa] supply-side water potential driving beta_stomata (Sabot 2022)
   end type leaf_env_t

   !----- Leaf-level fluxes returned by the solver. ---------------------------------------!
   type :: leaf_out_t
      real(wp)    :: A_net        !< [umol CO2/m2 leaf/s] net assimilation (A_gross - Rd)
      real(wp)    :: A_gross      !< [umol CO2/m2 leaf/s] gross assimilation
      real(wp)    :: gs           !< [mol H2O/m2 leaf/s]  stomatal conductance to water vapour
      real(wp)    :: ci           !< [umol/mol]      intercellular CO2 mole fraction
      real(wp)    :: cs           !< [umol/mol]      leaf-surface CO2 mole fraction (= ca if no BL)
      real(wp)    :: transpiration!< [mol H2O/m2 leaf/s]  E = g VPD / pressure: g = gs, or gs and gb in series
      real(wp)    :: rd           !< [umol/m2 leaf/s]     leaf respiration used
      integer(ik) :: limitation   !< LIM_* of the binding term
      logical     :: converged    !< .true. if the Ci solve met tolerance
      !----- The two water-stress limbs the kernel applied, reported so an inert control shows in the  !
      !      output (env%psi left at 0 makes beta_stomata identically 1). ---------------------------!
      real(wp)    :: beta_stomata    = 1.0_wp  !< [-] stomatal limb min(1, exp(sref*(psi-psi_onset))) (1 = unstressed)
      real(wp)    :: beta_nonstomata = 1.0_wp  !< [-] fast limb, the psi_leaf ramp on Vcmax/Jmax/TPU
   end type leaf_out_t

   !----- Everything the Ci residual reads, passed explicitly so the residual is a module function    !
   !      (compiler trap: no contained function handed to a solver, docs/building.md). -------------!
   type :: ci_problem_t
      integer(ik) :: pathway, colimitation, stomatal_model
      integer(ik) :: gs_rule = GS_FROM_MODEL
      real(wp)    :: vcmax, jrate, tpu, rd                    !< [umol/m2/s] at leaf temperature
      real(wp)    :: gstar, kc, ko, o2                        !< [umol/mol]
      real(wp)    :: aj_light, kp_eff                         !< C4 light-limited rate, PEP slope
      real(wp)    :: theta_cj_c3, theta_ip_c3, theta_cj_c4, theta_ic_c4
      real(wp)    :: ca                                       !< [umol/mol] ambient CO2
      real(wp)    :: gb                                       !< [mol/m2/s] boundary-layer conductance
      logical     :: boundary_layer                           !< draw leaf-surface CO2 down through gb
      real(wp)    :: vpd, ddef                                !< [Pa], [mol/mol] water deficit
      real(wp)    :: g0, g1, d0, lambda                       !< the stomatal model's parameters
      real(wp)    :: vpd_min                                  !< [kPa] the Medlyn model's least VPD
   end type ci_problem_t


contains

   !========== Photosynthetic demand (docs/science/leaf_gas_exchange.md §2) ===============!

   !---------------------------------------------------------------------------------------!
   ! Electron transport J: the smaller root of the non-rectangular hyperbola in absorbed   !
   ! PAR, I2 = 0.5 phi_psii absorptance PAR (§2).                                          !
   !---------------------------------------------------------------------------------------!
   elemental pure function electron_transport_j(par, absorptance, phi_psii, jmax, theta) result(j)
      real(wp), intent(in) :: par         !< [umol photon/m2/s] incident PAR
      real(wp), intent(in) :: absorptance !< [--] leaf PAR absorptance
      real(wp), intent(in) :: phi_psii    !< [--] electron yield of linear transport in low light
      real(wp), intent(in) :: jmax        !< [umol/m2/s] electron-transport capacity (T-scaled)
      real(wp), intent(in) :: theta       !< [--] curvature (0 < theta < 1)
      real(wp)             :: j, i2
      i2 = 0.5_wp * phi_psii * absorptance * par
      j  = quadratic_smaller_root(theta, i2, jmax)
   end function electron_transport_j

   !---------------------------------------------------------------------------------------!
   ! C3 demand: gross assimilation and the three raw limitation rates Ac/Aj/Ap, eq. (2).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colim,               &
                                   theta_cj, theta_ip, A_gross, Ac, Aj, Ap)
      real(wp),    intent(in)  :: ci, vcmax, j, tpu, gstar, kc, ko, o2, theta_cj, theta_ip
      integer(ik), intent(in)  :: colim
      real(wp),    intent(out) :: A_gross, Ac, Aj, Ap
      Ac = vcmax * (ci - gstar) / (ci + kc * (1.0_wp + o2 / ko))
      Aj = j     * (ci - gstar) / (4.0_wp * ci + 8.0_wp * gstar)
      Ap = 3.0_wp * tpu
      A_gross = combine_limits(Ac, Aj, Ap, colim, theta_cj, theta_ip)
   end subroutine assimilation_demand_c3

   !---------------------------------------------------------------------------------------!
   ! C4 demand (Collatz 1992): Ac = Vcmax, Aj = the light-limited rate, Ap = kp_eff*Ci     !
   ! (PEPcase); Gamma* ~ 0. The co-limitation curvatures are C4's own pair.                !
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
   ! Combine the three C3 rates: sharp min, or two nested smoothing quadratics with one    !
   ! curvature per transition, neither of them theta_j (docs/science/leaf_gas_exchange.md  !
   ! §2, "The co-limitation curvatures are not theta_J").                                  !
   !---------------------------------------------------------------------------------------!
   pure function combine_limits(Ac, Aj, Ap, colim, theta_cj, theta_ip) result(A)
      real(wp),    intent(in) :: Ac, Aj, Ap, theta_cj, theta_ip
      integer(ik), intent(in) :: colim
      real(wp)                :: A, Ai
      if (colim == COLIM_QUADRATIC) then
         Ai = quadratic_smaller_root(theta_cj, Ac, Aj)   ! co-limit Rubisco & RuBP
         A  = quadratic_smaller_root(theta_ip, Ai, Ap)   ! co-limit with product (TPU)
      else
         A = min(Ac, Aj, Ap)
      end if
   end function combine_limits


   !========== Stomatal conductance (docs/science/leaf_gas_exchange.md §3) ================!

   !---------------------------------------------------------------------------------------!
   ! Leuning (1995) semi-empirical stomatal conductance, eq. (5).                          !
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
   ! Medlyn et al. (2011) unified stomatal optimization (USO) conductance, eq. (6).        !
   !---------------------------------------------------------------------------------------!
   pure function stomata_gs_medlyn(a_net, cs, vpd, g0, g1, vpd_min) result(gs)
      real(wp), intent(in) :: a_net   !< [umol/m2/s] net assimilation
      real(wp), intent(in) :: cs      !< [umol/mol]  leaf-surface CO2
      real(wp), intent(in) :: vpd     !< [Pa]        leaf-to-air VPD
      real(wp), intent(in) :: g0, g1  !< [mol/m2/s], [kPa^0.5]
      real(wp), intent(in) :: vpd_min !< [kPa]       the VPD used at least (g1/sqrt(D) is undefined at 0)
      real(wp)             :: gs, vpd_kpa
      if (a_net <= 0.0_wp) then
         gs = g0
         return
      end if
      vpd_kpa = max(vpd * 1.0e-3_wp, vpd_min)
      gs = g0 + gsw_2_gsc * (1.0_wp + g1 / sqrt(vpd_kpa)) * a_net / max(cs, tiny_num)
   end function stomata_gs_medlyn

   !---------------------------------------------------------------------------------------!
   ! Katul marginal water-use efficiency lambda = lambda25 * beta^(-exp) (§4.2): a drier   !
   ! leaf pays more for water, so its optimal stomata close. exp = 0 disables it.          !
   !---------------------------------------------------------------------------------------!
   pure function katul_lambda(lambda25, beta, lambda_psi_exp) result(lambda)
      real(wp), intent(in) :: lambda25        !< [umol CO2/mol H2O] well-watered marginal WUE
      real(wp), intent(in) :: beta            !< [--] water-stress factor (0 = closed, 1 = open)
      real(wp), intent(in) :: lambda_psi_exp  !< [--] water-stress exponent
      real(wp)             :: lambda
      lambda = lambda25 * max(beta, BETA_FLOOR) ** (-lambda_psi_exp)
   end function katul_lambda


   !---------------------------------------------------------------------------------------!
   ! C3 Rubisco kinetics at leaf temperature, as mole fractions [umol/mol]. Kc and Ko are  !
   ! partial pressures, divided by the local pressure; Gamma* scales with the O2 mole      !
   ! fraction and not with pressure (docs/science/leaf_gas_exchange.md eqs. (2)-(3)).      !
   !---------------------------------------------------------------------------------------!
   pure subroutine rubisco_kinetics(p, t_leaf, pressure, kc, ko, gstar, o2)
      type(leaf_photo_pft_t),    intent(in)  :: p
      real(wp),                  intent(in)  :: t_leaf      !< [K]
      real(wp),                  intent(in)  :: pressure    !< [Pa] the air's pressure at the leaf
      real(wp),                  intent(out) :: kc, ko, gstar, o2   !< [umol/mol]
      kc    = arrhenius_scale(p%kc25, p%ea_kc, t_leaf) / pressure * mol_2_umol
      ko    = arrhenius_scale(p%ko25, p%ea_ko, t_leaf) / pressure * mol_2_umol
      gstar = arrhenius_scale(p%gstar25, p%ea_gstar, t_leaf) / P_REF_GSTAR * mol_2_umol                &
              * (p%o2_mol_frac / O2_REF_GSTAR)
      o2    = p%o2_mol_frac * mol_2_umol
   end subroutine rubisco_kinetics


   !========== The coupled Ci solver (docs/science/leaf_gas_exchange.md §5) ===============!

   !---------------------------------------------------------------------------------------!
   ! Solve the leaf A-gs-Ci system (docs/science/leaf_gas_exchange.md §1, §5). p carries   !
   ! every parameter, so the routine is testable without a meds_config_t. With            !
   ! use_boundary_layer (and gb > 0) the boundary layer draws Cs down from Ca (eq. (1))    !
   ! and sits in series with gs for transpiration; without it Cs = Ca.                     !
   !---------------------------------------------------------------------------------------!
   subroutine solve_leaf_gas_exchange(env, p, sm, temp_response_form, colim, use_boundary_layer, flux)
      type(leaf_env_t),          intent(in)  :: env
      type(leaf_photo_pft_t),    intent(in)  :: p
      integer(ik),               intent(in)  :: sm, temp_response_form, colim
      logical,                   intent(in)  :: use_boundary_layer
      type(leaf_out_t),         intent(out) :: flux

      real(wp) :: t_leaf, pressure, ca_ppm, o2_ppm, ddef, beta_nonstomata, beta_stomata, g1_eff
      real(wp) :: vcmax, jmax, jrate, tpu, rd, kc_ppm, ko_ppm, gstar_ppm
      real(wp) :: Aj_light, kp_eff, lambda_eff
      real(wp) :: lo0, hi0, ci_sol, An_open
      real(wp) :: A_gross, Ac, Aj, Ap, An, cs_sol, gs_sol
      logical  :: converged, do_boundary_layer
      type(ci_problem_t) :: prob

      t_leaf   = env%leaf_temp
      pressure = env%pressure
      ca_ppm   = env%ca
      ddef     = env%vpd / pressure                       ! mole-fraction water deficit D
      do_boundary_layer    = use_boundary_layer .and. env%gb > 0.0_wp

      call rubisco_kinetics(p, t_leaf, pressure, kc_ppm, ko_ppm, gstar_ppm, o2_ppm)
      vcmax = temp_response(temp_response_form, p%vcmax25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      jmax  = temp_response(temp_response_form, p%jmax25,  p%ea_jmax,  p%hd_jmax,  p%ds_jmax,  t_leaf)
      tpu   = temp_response(temp_response_form, p%tpu25,   p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      rd    = temp_response(temp_response_form, p%rd25,    p%ea_rd,    p%hd_rd,    p%ds_rd,    t_leaf)

      !----- Water stress, two controls on two time scales (docs/science/leaf_gas_exchange.md §4). The fast !
      !      one, eq. (9), ramps Vcmax/Jmax/TPU on the step-start leaf potential (Rd is not scaled); the     !
      !      slow one, eq. (10), scales g1 and the Katul lambda on the previous day's predawn potential. ---!
      beta_nonstomata = (env%psi_leaf - p%psi_close) / (p%psi_open - p%psi_close)
      beta_nonstomata = min(max(beta_nonstomata, 0.0_wp), 1.0_wp)
      vcmax = vcmax * beta_nonstomata
      jmax  = jmax  * beta_nonstomata
      tpu   = tpu   * beta_nonstomata
      !----- No stress above the onset, so a tall tree in wet soil keeps its full g1 (§4.2). ---------------!
      beta_stomata = min(1.0_wp, exp(p%sref_stomata * (env%psi - p%psi_onset)))
      g1_eff       = p%g1 * beta_stomata
      lambda_eff   = katul_lambda(p%lambda25, beta_stomata, p%lambda_psi_exp)

      !----- Light: C3 non-rectangular hyperbola J; C4 linear light-limited slope. --------!
      if (p%pathway == PATH_C4) then
         gstar_ppm = 0.0_wp
         Aj_light  = p%quantum_yield * p%absorptance * env%par
         !----- The PEPcase slope kp takes Vcmax's temperature response, as tpu does (§2). ------!
         kp_eff    = temp_response(temp_response_form, p%kp25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)       &
                     * pressure / p_std
         jrate     = 0.0_wp
      else
         jrate     = electron_transport_j(env%par, p%absorptance, p%phi_psii, jmax, p%theta_j)
         Aj_light  = 0.0_wp
         kp_eff    = 0.0_wp
      end if

      prob = ci_problem_t(pathway = p%pathway, colimitation = colim, stomatal_model = sm,            &
                          vcmax = vcmax, jrate = jrate, tpu = tpu, rd = rd,                           &
                          gstar = gstar_ppm, kc = kc_ppm, ko = ko_ppm, o2 = o2_ppm,                   &
                          aj_light = Aj_light, kp_eff = kp_eff,                                       &
                          theta_cj_c3 = p%theta_cj_c3, theta_ip_c3 = p%theta_ip_c3,                   &
                          theta_cj_c4 = p%theta_cj_c4, theta_ic_c4 = p%theta_ic_c4,                   &
                          ca = ca_ppm, gb = env%gb, boundary_layer = do_boundary_layer,               &
                          vpd = env%vpd, ddef = ddef, g0 = p%g0, g1 = g1_eff, d0 = p%d0,              &
                          vpd_min = p%medlyn_vpd_min,                                                 &
                          lambda = lambda_eff)

      !----- Closed/night branch: no positive-assimilation root (best-case net <= 0). ------!
      An_open = ci_net_assimilation(prob, ca_ppm)
      if (An_open <= 0.0_wp) then
         gs_sol = p%g0
         An     = An_open
         cs_sol = ca_ppm
         if (do_boundary_layer) cs_sol = ca_ppm - gbw_2_gbc * An / env%gb
         if (gs_sol > tiny_num) then
            ci_sol = cs_sol - gsw_2_gsc * An / gs_sol
         else
            ci_sol = cs_sol                               ! stomata shut: no exchange to set Ci by
         end if
         call fill_flux(An + rd, An, gs_sol, ci_sol, cs_sol, rd, LIM_NONE, .true.)
         return
      end if

      !----- Bracket Ci in (Gamma*, Ca] and find the root. With no sign change (Katul under strong stress) !
      !       fall back to the cuticular g0, which always brackets when net A(Ca) > 0 (§5). ---------!
      lo0 = gstar_ppm + LO_EPS_PPM
      hi0 = ca_ppm
      prob%gs_rule = GS_FROM_MODEL
      do                                                 ! at most two passes: model, then g0-pinned
         call root_ci(prob, lo0, hi0, ci_sol, converged)
         !----- No sign change -> retry with the cuticular floor g0. ----------------------------!
         if (.not. converged .and. prob%gs_rule /= GS_CUTICULAR) then
            prob%gs_rule = GS_CUTICULAR
            cycle
         end if

         !----- Assemble the solution: net A, surface CO2, back-computed gs, transpiration. ---!
         call ci_assimilation_demand(prob, ci_sol, A_gross, Ac, Aj, Ap)
         An     = A_gross - rd
         cs_sol = ca_ppm
         if (do_boundary_layer) cs_sol = ca_ppm - gbw_2_gbc * An / env%gb
         !----- Katul optimum can land below the cuticular floor g0; re-solve once g0-pinned so   !
         !       A/gs/Ci/E stay mutually consistent (Leuning/Medlyn already return gs >= g0). ----!
         if (sm == SM_KATUL .and. prob%gs_rule == GS_FROM_MODEL .and. cs_sol - ci_sol > tiny_num) then
            if (gsw_2_gsc * An / (cs_sol - ci_sol) < p%g0) then
               prob%gs_rule = GS_CUTICULAR
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
         !----- Transpiration on the total leaf-to-air conductance, gs in series with gb, as the  !
         !      CO2 solve does. ------------------------------------------------------------------!
         if (do_boundary_layer) then
            flux%transpiration = gs * env%gb / (gs + env%gb) * env%vpd / pressure
         else
            flux%transpiration = gs * env%vpd / pressure
         end if
         flux%rd         = rd_loc
         flux%limitation = lim
         flux%converged  = conv
         !----- Report the two water-stress limbs the kernel just applied (diagnostic only). ---!
         flux%beta_stomata    = beta_stomata
         flux%beta_nonstomata = beta_nonstomata
      end subroutine fill_flux

   end subroutine solve_leaf_gas_exchange

   !----- Gross assimilation and the three raw limitation rates at a trial Ci. ----------------!
   pure subroutine ci_assimilation_demand(prob, ci, Ag, Ac, Aj, Ap)
      type(ci_problem_t), intent(in)  :: prob
      real(wp),           intent(in)  :: ci
      real(wp),           intent(out) :: Ag, Ac, Aj, Ap
      if (prob%pathway == PATH_C4) then
         call assimilation_demand_c4(ci, prob%vcmax, prob%aj_light, prob%kp_eff, prob%colimitation,  &
                                     prob%theta_cj_c4, prob%theta_ic_c4, Ag, Ac, Aj, Ap)
      else
         call assimilation_demand_c3(ci, prob%vcmax, prob%jrate, prob%tpu, prob%gstar, prob%kc,      &
                                     prob%ko, prob%o2, prob%colimitation, prob%theta_cj_c3,          &
                                     prob%theta_ip_c3, Ag, Ac, Aj, Ap)
      end if
   end subroutine ci_assimilation_demand

   !----- Net assimilation at a trial Ci. -----------------------------------------------------!
   pure real(wp) function ci_net_assimilation(prob, ci) result(an)
      type(ci_problem_t), intent(in) :: prob
      real(wp),           intent(in) :: ci
      real(wp) :: Ag, Ac, Aj, Ap
      call ci_assimilation_demand(prob, ci, Ag, Ac, Aj, Ap)
      an = Ag - prob%rd
   end function ci_net_assimilation

   !----- The residual whose root is the solved Ci: Katul's optimality condition on its own     !
   !      pass, and otherwise the CO2 diffusion identity with an explicit conductance. -----------!
   pure real(wp) function ci_residual(prob, ci) result(r)
      type(ci_problem_t), intent(in) :: prob
      real(wp),           intent(in) :: ci
      if (prob%gs_rule == GS_FROM_MODEL .and. prob%stomatal_model == SM_KATUL) then
         r = residual_optimality(prob, ci)
      else
         r = residual_explicit_gs(prob, ci)
      end if
   end function ci_residual

   !----- Explicit-conductance residual for Leuning / Medlyn and the cuticular fallback: the     !
   !       trial Ci must satisfy the CO2 diffusion identity Ci = Cs - A / gs_co2 (eq. (1)),      !
   !       where the current conductance rule supplies gs. Returns Ci - Ci_predicted. -----------!
   pure real(wp) function residual_explicit_gs(prob, ci) result(r)
      type(ci_problem_t), intent(in) :: prob
      real(wp),           intent(in) :: ci
      real(wp) :: An_loc, cs_surf, gs, ci_pred
      An_loc  = ci_net_assimilation(prob, ci)            ! net assimilation A at this trial Ci
      !----- Leaf-surface CO2: ambient, less the boundary-layer drawdown by the CO2 flux. ------!
      cs_surf = prob%ca
      if (prob%boundary_layer) cs_surf = prob%ca - gbw_2_gbc * An_loc / prob%gb
      !----- Stomatal conductance from the current rule. ---------------------------------------!
      select case (prob%gs_rule)
      case (GS_CUTICULAR)
         gs = prob%g0                                    ! closed-stomata fallback: gs pinned to g0
      case default
         if (prob%stomatal_model == SM_LEUNING) then
            gs = stomata_gs_leuning(An_loc, cs_surf, prob%gstar, prob%vpd, prob%g0, prob%g1, prob%d0)
         else
            gs = stomata_gs_medlyn(An_loc, cs_surf, prob%vpd, prob%g0, prob%g1, prob%vpd_min)
         end if
      end select
      !----- Ci predicted by CO2 diffusion through the stomata (gs is a water conductance, so     !
      !       the CO2 conductance is gs / gsw_2_gsc); residual = trial Ci minus predicted Ci. ---!
      ci_pred = cs_surf - gsw_2_gsc * An_loc / max(gs, tiny_num)
      r       = ci - ci_pred
   end function residual_explicit_gs

   !----- Katul optimality residual (no explicit gs), eqs. (7)-(8): the first-order condition of    !
   !       max(A - lambda*E) in Ci-space, whose root is the optimal Ci. --------------------------!
   pure real(wp) function residual_optimality(prob, ci) result(r)
      type(ci_problem_t), intent(in) :: prob
      real(wp),           intent(in) :: ci
      real(wp) :: An_loc, cs_surf, dAn_dci, dci
      An_loc  = ci_net_assimilation(prob, ci)            ! net assimilation A at this trial Ci
      !----- Leaf-surface CO2 (ambient less boundary-layer drawdown), as for the gs models. -----!
      cs_surf = prob%ca
      if (prob%boundary_layer) cs_surf = prob%ca - gbw_2_gbc * An_loc / prob%gb
      !----- Marginal demand A' = dA/dCi by central difference (A(Ci) is the co-limited FvCB      !
      !       envelope, so the slope is taken numerically; dci is a relative step, abs-floored). --!
      dci     = max(1.0e-3_wp * abs(ci), 1.0e-2_wp)   ! clamp-ok: the finite-difference step [umol/mol]
      dAn_dci = (ci_net_assimilation(prob, ci + dci) - ci_net_assimilation(prob, ci - dci))          &
                / (2.0_wp * dci)
      !----- First-order optimality: A'(Cs-Ci)^2 = gsw_2_gsc * D * lambda * (A'(Cs-Ci) + A). ----!
      r = dAn_dci * (cs_surf - ci)**2                                                              &
          - gsw_2_gsc * prob%ddef * prob%lambda * (dAn_dci * (cs_surf - ci) + An_loc)
   end function residual_optimality

   !----- The root of ci_residual on [lo, hi], to CI_TOL_PPM, by the Illinois method            !
   !      (meds_numerics: superlinear, far fewer residuals than bisection's ~21). It stops when   !
   !      the bracket is within CI_TOL_PPM: a small move alone is not enough, as regula falsi can !
   !      creep. A same-sign bracket returns the midpoint with converged = .false., so the caller !
   !      can change the conductance rule and retry. ---------------------------------------------!
   pure subroutine root_ci(prob, lo, hi, ci, converged)
      type(ci_problem_t), intent(in)  :: prob
      real(wp),           intent(in)  :: lo, hi
      real(wp),           intent(out) :: ci
      logical,            intent(out) :: converged
      type(root_bracket_t) :: br
      real(wp)    :: fc
      integer(ik) :: it
      br%a  = lo
      br%b  = hi
      br%fa = ci_residual(prob, br%a)
      br%fb = ci_residual(prob, br%b)
      converged = .false.
      ci = 0.5_wp * (br%a + br%b)
      if (br%fa * br%fb > 0.0_wp) return
      converged = .true.
      if (br%fa == 0.0_wp .or. br%fb == 0.0_wp) then
         ci = merge(br%a, br%b, br%fa == 0.0_wp)
         return
      end if
      do it = 1_ik, MAX_ITER
         ci = regula_falsi_point(br)
         fc = ci_residual(prob, ci)
         if (fc == 0.0_wp) exit
         call illinois_narrow(br, ci, fc)
         if (br%b - br%a < CI_TOL_PPM) exit
      end do
      converged = fc == 0.0_wp .or. br%b - br%a < CI_TOL_PPM
   end subroutine root_ci

   !----- leaf_params_at_capacity -- a PFT's table entry with one leaf's plastic capacities on top:   !
   !      Jmax25 and TPU25 scale with the overriding Vcmax25. The one place a cohort's capacities meet !
   !      the table: leaf_gas_exchange_batch and the canopy C API (meds_c_api_canopy) both take it. --!
   pure function leaf_params_at_capacity(table, ipft, vcmax25, rd25) result(p)
      type(leaf_photo_table_t), intent(in) :: table
      integer(ik),              intent(in) :: ipft
      real(wp),                 intent(in) :: vcmax25, rd25
      type(leaf_photo_pft_t)               :: p
      p         = table%pft(ipft)
      p%vcmax25 = vcmax25
      p%jmax25  = table%jmax_vcmax_ratio(ipft) * vcmax25
      p%tpu25   = table%tpu_vcmax_ratio(ipft)  * vcmax25
      p%rd25    = rd25
   end function leaf_params_at_capacity

   !---------------------------------------------------------------------------------------!
   ! leaf_gas_exchange_batch -- solve_leaf_gas_exchange over n leaves on bare arrays (the  !
   ! fast-loop prepass). Per-leaf inputs are arrays; ca, pressure and the leaf table are   !
   ! shared. psi absent means 0 (the slow stomatal control inert). A_gross, gs and rd are  !
   ! what the fast loop reads; the optional outputs are diagnostics.                       !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_gas_exchange_batch(n, par, leaf_temp, vpd, ca, pressure, psi_leaf, gb,      &
                                      table, pft, vcmax25, rd25, a_gross, gs, rd, psi,        &
                                      a_net, ci, cs, transp, limitation, beta_stom, beta_nonstom)
      integer(ik),         intent(in)  :: n
      real(wp),            intent(in)  :: par(n), leaf_temp(n), vpd(n), psi_leaf(n), gb(n)  !< per-leaf env
      real(wp),            intent(in)  :: ca, pressure                                      !< patch-uniform (broadcast)
      type(leaf_photo_table_t), intent(in) :: table                                       !< per-PFT parameters (once per run)
      integer(ik),         intent(in)  :: pft(n)
      real(wp),            intent(in)  :: vcmax25(n), rd25(n)                               !< per-leaf plastic capacities
      real(wp),            intent(out) :: a_gross(n), gs(n), rd(n)
      real(wp), optional,  intent(in)  :: psi(n)                                       !< absent => 0 (well-watered)
      !----- Diagnostic-only outputs (the per-cohort ecophysiology block). --------------------!
      real(wp),    optional, intent(out) :: a_net(n), ci(n), cs(n), transp(n)
      real(wp),    optional, intent(out) :: beta_stom(n), beta_nonstom(n)
      integer(ik), optional, intent(out) :: limitation(n)
      type(leaf_env_t)          :: env
      type(leaf_out_t)         :: flux
      type(leaf_photo_pft_t)    :: p
      integer(ik) :: i
      do i = 1_ik, n
         env%par = par(i)
         env%leaf_temp = leaf_temp(i)
         env%vpd = vpd(i)
         env%ca = ca
         env%pressure = pressure
         env%psi_leaf = psi_leaf(i)
         env%gb = gb(i)
         if (present(psi)) then
            env%psi = psi(i)
         else
            env%psi = 0.0_wp
         end if
         p = leaf_params_at_capacity(table, pft(i), vcmax25(i), rd25(i))
         call solve_leaf_gas_exchange(env, p, table%stomatal_model, table%temp_response_form,   &
                                      table%colimitation, table%use_boundary_layer, flux)
         a_gross(i) = flux%A_gross
         gs(i) = flux%gs
         rd(i) = flux%rd
         if (present(a_net))        a_net(i)        = flux%A_net
         if (present(ci))           ci(i)           = flux%ci
         if (present(cs))           cs(i)           = flux%cs
         if (present(transp))       transp(i)       = flux%transpiration
         if (present(limitation))   limitation(i)   = flux%limitation
         if (present(beta_stom))    beta_stom(i)    = flux%beta_stomata
         if (present(beta_nonstom)) beta_nonstom(i) = flux%beta_nonstomata
      end do
   end subroutine leaf_gas_exchange_batch

end module meds_leaf_gas_exchange
