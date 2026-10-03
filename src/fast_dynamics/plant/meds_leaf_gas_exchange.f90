! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_leaf_gas_exchange -- leaf-level gas-exchange COMPUTE kernels (photosynthesis + stomata !
! + coupled Ci solver), merged into one module. FvCB C3 / Collatz C4 demand, the electron-     !
! transport hyperbola, the Leuning / Medlyn / Katul stomatal models, and the bracketed Ci      !
! root-find (solve_leaf_gas_exchange). The public seam leaf_gas_exchange lives in               !
! meds_fast_config; the raw kernels here are also called directly by meds_c_api_leaf.        !
!==========================================================================================!
module meds_leaf_gas_exchange
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : p_std, tiny_num, gsw_2_gsc, gbw_2_gbc, mol_2_umol
   use meds_leaf_opts,     only : COLIM_MIN, COLIM_QUADRATIC, SM_LEUNING, SM_MEDLYN, SM_KATUL
   use meds_plant_types, only : leaf_env_t, leaf_flux_t, leaf_photo_params_t, leaf_photo_table_t, PATH_C3, PATH_C4, &
                                LIM_NONE, LIM_RUBISCO, &
                                LIM_RUBP, LIM_PRODUCT, LIM_C4_PEP
   use meds_temp_response, only : temp_response, arrhenius_scale
   use meds_numerics,      only : quadratic_smaller_root
   implicit none

   !----- O2 mole fraction the shipped gstar25 was MEASURED at (Bernacchi et al. 2001). Gamma* is  !
   !      proportional to O2, so this is the reference the config's o2_mol_frac is scaled against. !
   real(wp), parameter :: O2_REF_GSTAR = 0.209_wp   !< [mol/mol]
   private

   !----- from meds_leaf_photosynthesis.f90 ----------------------------------------------!

   public :: assimilation_demand_c3, assimilation_demand_c4, electron_transport_j

   !----- from meds_leaf_stomata.f90 -----------------------------------------------------!

   public :: stomata_gs_leuning, stomata_gs_medlyn, katul_lambda, low_psi_gs_factor

   real(wp), parameter :: beta_floor   = 1.0e-4_wp   !< [--] water-stress floor (bound lambda as beta->0)

   !----- from meds_leaf_solver.f90 ------------------------------------------------------!

   public :: solve_leaf_gas_exchange
   public :: leaf_gas_exchange_batch, leaf_params_at_capacity

   real(wp),    parameter :: ci_tol_ppm = 1.0e-3_wp    !< [umol/mol] Ci convergence tolerance (~1e-4 Pa)
   real(wp),    parameter :: lo_eps_ppm = 1.0e-3_wp    !< [umol/mol] offset of the lower bracket above Gamma*
   integer(ik), parameter :: max_iter   = 100_ik       !< bisection iteration cap (safety net)

   !----- Which stomatal conductance a pass of the Ci solve uses (solve_leaf_gas_exchange). -----!
   integer(ik), parameter :: GS_FROM_MODEL = 1_ik      !< the configured stomatal model
   integer(ik), parameter :: GS_CUTICULAR  = 2_ik      !< the cuticular floor g0 (the fallback)
   integer(ik), parameter :: GS_PINNED     = 3_ik      !< a fixed conductance (Katul, low potential)

   !----- Everything the Ci residual reads: the leaf's biochemistry at its temperature, its air,  !
   !      its stomatal model, and the conductance rule of the current pass. It is passed to the   !
   !      residual explicitly, so the residual is an ordinary module function: under ifx, handing !
   !      a contained function to a solver costs a lock-guarded allocation on every call, which   !
   !      made more than four threads slower than four (#325). ------------------------------------!
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
      real(wp)    :: f_lwp                                    !< low-water-potential factor on gs
      real(wp)    :: gs_pin = 0.0_wp                          !< [mol/m2/s] the GS_PINNED conductance
   end type ci_problem_t


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
   !                                                                                          !
   ! TWO CURVATURES, ONE PER TRANSITION, and neither of them is theta_j (#118). A co-limitation !
   ! curvature says how sharply the leaf switches between two LIMITING PROCESSES; theta_j is the !
   ! curvature of the electron-transport hyperbola, a different quantity that happens to live in !
   ! the same units. C3 used theta_j for both smoothings here while C4 already had its own pair. !
   ! Because Ac and Aj sit close together at ambient CO2, the penalty landed exactly where the   !
   ! model spends most of its time: ~29 % of assimilation against min(Ac,Aj,Ap), and the shortfall!
   ! is nearly independent of Vcmax, so no measured Vcmax reproduced a measured rate.            !
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
   !----- The low-water-potential factor on the calculated stomatal conductance (#332): 1 at or above  !
   !      the turgor-loss point, falling linearly to 0 at twice it, 0 below. psi is the previous      !
   !      day's daily-max (predawn) leaf potential [MPa]. -----------------------------------------!
   elemental pure function low_psi_gs_factor(psi, psi_tlp) result(f)
      real(wp), intent(in) :: psi, psi_tlp      !< [MPa]; psi_tlp < 0
      real(wp) :: f
      if (psi_tlp >= 0.0_wp) then
         f = 1.0_wp
      else
         f = min(max((psi - 2.0_wp * psi_tlp) / (-psi_tlp), 0.0_wp), 1.0_wp)
      end if
   end function low_psi_gs_factor

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
      real(wp) :: f_lwp             !< low-water-potential factor on the calculated gs (1 = none)
      real(wp) :: gs_pin            !< Katul: the conductance the factor pins the solve to
      real(wp) :: vcmax, jmax, jrate, tpu, rd, kc_ppm, ko_ppm, gstar_ppm
      real(wp) :: Aj_light, kp_eff, lambda_eff
      real(wp) :: lo0, hi0, ci_sol, An_open
      real(wp) :: A_gross, Ac, Aj, Ap, An, cs_sol, gs_sol
      logical  :: converged, do_boundary_layer
      type(ci_problem_t) :: prob

      t_leaf   = env%leaf_temp
      pressure = env%pressure
      ca_ppm   = env%ca
      o2_ppm   = p%o2_mol_frac * mol_2_umol
      ddef     = env%vpd / pressure                       ! mole-fraction water deficit D
      do_boundary_layer    = use_boundary_layer .and. env%gb > 0.0_wp

      !----- Temperature-scale the biochemistry (Kc/Ko/Gamma* always Arrhenius; Pa -> ppm). -!
      kc_ppm    = arrhenius_scale(p%kc25,    p%ea_kc,    t_leaf) / pressure * mol_2_umol
      ko_ppm    = arrhenius_scale(p%ko25,    p%ea_ko,    t_leaf) / pressure * mol_2_umol
      !----- Gamma* is set by Rubisco's CO2/O2 specificity, Gamma* = 0.5*O/S_(c/o), so it is       !
      !      PROPORTIONAL to the O2 partial pressure. Scale it off the O2 the shipped gstar25 was  !
      !      measured at (Bernacchi et al. 2001, 21% O2), so the parameter's provenance is stated  !
      !      in code rather than implied. Without this, o2_mol_frac reached only Kc(1+O/Ko) and    !
      !      not the compensation point, so raising O2 inhibited carboxylation but left the whole  !
      !      photorespiratory penalty on Aj untouched -- an error of exactly zero at 20.9% growing !
      !      monotonically in both directions (-13% of A at 10% O2, +20% at 35%), which is the     !
      !      signature of a missing scaling rather than anything else (#117).                      !
      !      At the shipped 0.209 the factor is EXACTLY 1, so the default path is bit-identical.   !
      gstar_ppm = arrhenius_scale(p%gstar25, p%ea_gstar, t_leaf) / pressure * mol_2_umol           &
                  * (p%o2_mol_frac / O2_REF_GSTAR)
      vcmax = temp_response(tresp, p%vcmax25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      jmax  = temp_response(tresp, p%jmax25,  p%ea_jmax,  p%hd_jmax,  p%ds_jmax,  t_leaf)
      tpu   = temp_response(tresp, p%tpu25,   p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, t_leaf)
      rd    = temp_response(tresp, p%rd25,    p%ea_rd,    p%hd_rd,    p%ds_rd,    t_leaf)

      !----- Water stress, split into two independently-tunable limbs (Sabot 2022 / Zhou 2013): !
      !   beta_nonstomata -- capacity limb: a linear psi_LEAF ramp downregulating Vcmax/Jmax/TPU, !
      !     applied to ALL stomatal models (a leaf-biochemistry effect, scheme-independent).      !
      !   beta_stomata    -- stomatal limb: 1 above the onset potential psi_onset and             !
      !     exp(sref*(psi - psi_onset)) below it, downregulating the                               !
      !     Leuning/Medlyn slope g1 and the Katul marginal WUE lambda (lambda ~                   !
      !     beta_stomata^(-lambda_psi_exp); lambda_psi_exp = 2 recovers Sabot's g1<->lambda).     !
      !----- The capacity limb is OFF by default (issue #47): rarely measured directly, weakly  !
      !      constrained, and a LINEAR AMPLIFIER on psi_leaf -- see meds_config_t%leaf_wstress_  !
      !      nonstomatal. beta = 1 when off, so Vcmax/Jmax/TPU pass through untouched. ----------!
      beta_nonstomata = 1.0_wp
      if (p%wstress_nonstomatal) then
         beta_nonstomata = (env%psi_leaf - p%psi_close) / (p%psi_open - p%psi_close)
         beta_nonstomata = min(max(beta_nonstomata, 0.0_wp), 1.0_wp)
         vcmax = vcmax * beta_nonstomata
         jmax  = jmax  * beta_nonstomata
         tpu   = tpu   * beta_nonstomata
      end if
      !----- No stress above the onset (Sabot et al. 2022 Eq. 5 has none while the soil is at field     !
      !      capacity), so a tall tree in wet soil, whose predawn potential is its gravity head, keeps   !
      !      its full g1. Below the onset the decline is exponential, at the rate sref. ---------------!
      beta_stomata = min(1.0_wp, exp(p%sref_stomata * (env%psi - p%psi_onset)))
      !----- LOW-WATER-POTENTIAL CONTROL (#332). The Sabot beta above scales g1 only, so as it goes to  !
      !      0 the conductance falls to the RESIDUAL g0 and never reaches zero -- measured at ~2.6       !
      !      mm/day of transpiration still leaving a plant whose wood store was empty and whose predawn  !
      !      potential was -116 MPa (#95). So the conductance the stomatal model calculates, g0          !
      !      included, is multiplied by f_lwp, which falls LINEARLY from 1 at the turgor-loss point to  !
      !      0 at twice it. The solve below stays coupled: A and Ci are those of the reduced gs, so a    !
      !      leaf with f_lwp = 0 exchanges no CO2 or water. It replaced a hard shutdown at 2*psi_tlp, a  !
      !      step no calibration could see past.                                                       !
      !                                                                                          !
      !      CARBON NOTE (#332, option kept by the owner). Rd is untouched: it is computed as before  !
      !      and charged in full as leaf respiration. What the factor changes is GROSS A. The solve   !
      !      ties net A to the stomata, net A = gs*(Cs - Ci)/1.6, so as gs -> 0 net A -> 0 and the     !
      !      leaf photosynthesises at the internal compensation point on the CO2 it respires: gross   !
      !      A -> Rd, its respiration refixed. The canopy sums gross A as GPP and Rd as leaf          !
      !      respiration (canopy_leaf_gas_exchange, meds_fast_prepass), so a fully closed leaf is     !
      !      carbon-neutral by day -- credited GPP = Rd, charged Rd -- where the former clamp set     !
      !      gross A = 0 and lost Rd. At night gross A = 0 and the leaf loses Rd, as before. This is  !
      !      the standard coupled leaf model, and the same refixation already happens whenever gs is  !
      !      small. A tower's GPP, partitioned from NEE, cannot see refixed CO2. ---------------------!
      !                                                                                          !
      !      The driver feeds env%psi the previous day's daily-MAX leaf water potential -- the   !
      !      model's predawn potential -- so f_lwp is set once a day from a slow, integrated measure,  !
      !      not from a noisy sub-daily psi. Phenology compares dmax against the same tlp            !
      !      (pheno_state_t%low_psi_days). A thermodynamic limit on transpiration (issue #96) would    !
      !      make this control matter less. ----------------------------------------------------!
      f_lwp = 1.0_wp
      if (p%low_psi_control == 1_ik) f_lwp = low_psi_gs_factor(env%psi, p%psi_tlp)   ! linear decline
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

      prob = ci_problem_t(pathway = p%pathway, colimitation = colim, stomatal_model = sm,            &
                          vcmax = vcmax, jrate = jrate, tpu = tpu, rd = rd,                           &
                          gstar = gstar_ppm, kc = kc_ppm, ko = ko_ppm, o2 = o2_ppm,                   &
                          aj_light = Aj_light, kp_eff = kp_eff,                                       &
                          theta_cj_c3 = p%theta_cj_c3, theta_ip_c3 = p%theta_ip_c3,                   &
                          theta_cj_c4 = p%theta_cj_c4, theta_ic_c4 = p%theta_ic_c4,                   &
                          ca = ca_ppm, gb = env%gb, boundary_layer = do_boundary_layer,               &
                          vpd = env%vpd, ddef = ddef, g0 = p%g0, g1 = g1_eff, d0 = p%d0,              &
                          vpd_min = p%medlyn_vpd_min,                                                 &
                          lambda = lambda_eff, f_lwp = f_lwp)

      !----- Closed/night branch: no positive-assimilation root (best-case net <= 0). ------!
      An_open = ci_net_assimilation(prob, ca_ppm)
      if (An_open <= 0.0_wp) then
         gs_sol = f_lwp * p%g0
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

      !----- Bracket Ci in (Gamma*, Ca] and bisect the residual to ci_tol_ppm. If the chosen    !
      !       stomatal model yields no consistent open solution (no sign change -- e.g. Katul      !
      !       under strong water stress where lambda -> large), fall back to the cuticular floor   !
      !       g0 (closed stomata), which always brackets when net A(Ca) > 0. The explicit-gs       !
      !       residual serves Leuning/Medlyn, the g0 fallback and the pinned pass; Katul's own     !
      !       pass uses the optimality residual (ci_residual). ------------------------------------!
      lo0 = gstar_ppm + lo_eps_ppm
      hi0 = ca_ppm
      prob%gs_rule = GS_FROM_MODEL
      do                                                 ! at most three passes: model, g0-pinned, f_lwp-pinned
         call bisect_ci(prob, lo0, hi0, ci_sol, converged)
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
         !----- Katul: the optimum is the calculated gs; scale it by f_lwp and re-solve with gs      !
         !       pinned there, so A/gs/Ci/E stay consistent (Leuning/Medlyn scale inside the solve). -!
         if (sm == SM_KATUL .and. prob%gs_rule == GS_FROM_MODEL .and. f_lwp < 1.0_wp) then
            gs_pin = p%g0
            if (cs_sol - ci_sol > tiny_num) gs_pin = max(gsw_2_gsc * An / (cs_sol - ci_sol), p%g0)
            prob%gs_pin  = f_lwp * gs_pin
            prob%gs_rule = GS_PINNED
            cycle
         end if
         exit
      end do
      !----- f_lwp = 0: the solve's limit, written exactly. With no conductance the leaf exchanges   !
      !       nothing -- net A = 0 and gross A = Rd, its respiration refixed at the Ci the solve found  !
      !       (see the CARBON NOTE above: the canopy counts that gross A as GPP). ---------------------!
      if (f_lwp <= 0.0_wp) then
         call fill_flux(rd, 0.0_wp, 0.0_wp, ci_sol, ca_ppm, rd, LIM_NONE, converged)
         return
      end if
      !----- Back-compute gs from the diffusion identity; if the boundary layer pushed Cs at  !
      !       or below Ci (degenerate), pin gs to g0 and report Cs as the surface CO2. -------!
      if (cs_sol - ci_sol > tiny_num) then
         gs_sol = max(gsw_2_gsc * An / (cs_sol - ci_sol), f_lwp * p%g0)
      else
         gs_sol = f_lwp * p%g0
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

   !----- Explicit-conductance residual for Leuning / Medlyn, the cuticular fallback and the     !
   !       pinned pass: the trial Ci must satisfy the CO2 diffusion identity Ci = Cs - A / gs_co2, !
   !       where the current conductance rule supplies gs. Returns Ci - Ci_predicted. ----------!
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
         gs = prob%f_lwp * prob%g0                       ! closed-stomata fallback: gs pinned to g0
      case (GS_PINNED)
         gs = prob%gs_pin                                ! Katul, scaled by the low-psi factor
      case default
         if (prob%stomatal_model == SM_LEUNING) then
            gs = prob%f_lwp * stomata_gs_leuning(An_loc, cs_surf, prob%gstar, prob%vpd, prob%g0,       &
                                                 prob%g1, prob%d0)
         else
            gs = prob%f_lwp * stomata_gs_medlyn(An_loc, cs_surf, prob%vpd, prob%g0, prob%g1, prob%vpd_min)
         end if
      end select
      !----- Ci predicted by CO2 diffusion through the stomata (gs is a WATER conductance, so     !
      !       the CO2 conductance is gs / gsw_2_gsc); residual = trial Ci minus predicted Ci. ---!
      ci_pred = cs_surf - gsw_2_gsc * An_loc / max(gs, tiny_num)
      r       = ci - ci_pred
   end function residual_explicit_gs

   !----- Katul optimality residual (no explicit gs). Stomata maximize A - lambda*E; the         !
   !       first-order condition dA/dCi = lambda * dE/dCi, with CO2 supply                        !
   !       A = gs/gsw_2_gsc * (Cs - Ci) and water loss E = gs * D (D = VPD/P), rearranges to the   !
   !       residual below, whose root is the optimal Ci. ---------------------------------------!
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
      dci     = max(1.0e-3_wp * abs(ci), 1.0e-2_wp)
      dAn_dci = (ci_net_assimilation(prob, ci + dci) - ci_net_assimilation(prob, ci - dci))          &
                / (2.0_wp * dci)
      !----- First-order optimality: A'(Cs-Ci)^2 = gsw_2_gsc * D * lambda * (A'(Cs-Ci) + A). ----!
      r = dAn_dci * (cs_surf - ci)**2                                                              &
          - gsw_2_gsc * prob%ddef * prob%lambda * (dAn_dci * (cs_surf - ci) + An_loc)
   end function residual_optimality

   !----- Bisect ci_residual on [lo, hi] to ci_tol_ppm. A same-sign bracket returns the midpoint  !
   !      with converged = .false., so the caller can change the conductance rule and retry. ----!
   pure subroutine bisect_ci(prob, lo, hi, ci, converged)
      type(ci_problem_t), intent(in)  :: prob
      real(wp),           intent(in)  :: lo, hi
      real(wp),           intent(out) :: ci
      logical,            intent(out) :: converged
      real(wp)    :: a, b, flo, fhi, mid, fmid
      integer(ik) :: it
      a = lo ; b = hi
      flo = ci_residual(prob, a) ; fhi = ci_residual(prob, b)
      converged = .false. ; ci = 0.5_wp * (a + b)
      if (flo * fhi <= 0.0_wp) then
         do it = 1_ik, max_iter
            mid = 0.5_wp * (a + b) ; fmid = ci_residual(prob, mid)
            if (flo * fmid <= 0.0_wp) then ; b = mid ; else ; a = mid ; flo = fmid ; end if
            if (b - a < ci_tol_ppm) exit
         end do
         ci = 0.5_wp * (a + b) ; converged = (b - a < ci_tol_ppm)
      end if
   end subroutine bisect_ci

   !---------------------------------------------------------------------------------------!
   ! leaf_gas_exchange_batch -- BARE-ARRAY entry point over n leaves (MEDS_NUMERICS_SCOPING.md    !
   ! "bare-array process kernels": one call solves a whole array of leaf environments, so the      !
   ! per-cohort fast-loop driver loops need no longer thread the leaf_env_t/leaf_flux_t derived     !
   ! types, and a Python/ctypes wrapper (meds_c_api_leaf) can vectorise over numpy arrays instead   !
   ! of calling one leaf at a time). The per-leaf PHYSICS is UNCHANGED: this loops `do i=1,n`        !
   ! calling the SAME leaf_gas_exchange above, so it is bit-identical to an inline caller loop.      !
   !                                                                                          !
   ! CONVENTION (shared by every *_batch kernel): genuinely PER-ELEMENT quantities are bare arrays   !
   ! of length n; PATCH/RUN-UNIFORM quantities (here ca, pressure, and the whole cfg trait table)    !
   ! are passed as scalars/one config object (broadcast to every element); outputs are bare arrays.  !
   ! psi is OPTIONAL (absent => 0, the leaf_env_t default = the drought-stomata limb inert,     !
   ! matching the current fast-loop wiring). A_gross / gs / rd are what the fast loop CONSUMES and   !
   ! are mandatory; the remaining leaf_flux_t fields are OPTIONAL outputs for DIAGNOSTICS only --    !
   ! absent means the caller does not report per-cohort ecophysiology, and nothing extra is copied.  !
   !---------------------------------------------------------------------------------------!
   !----- leaf_params_at_capacity -- a PFT's table entry with one leaf's plastic capacities on top:   !
   !      Jmax25 and TPU25 scale with the overriding Vcmax25. The one place a cohort's capacities meet !
   !      the table: leaf_gas_exchange_batch and the canopy C API (meds_c_api_canopy) both take it. --!
   pure function leaf_params_at_capacity(table, ipft, vcmax25, rd25) result(p)
      type(leaf_photo_table_t), intent(in) :: table
      integer(ik),              intent(in) :: ipft
      real(wp),                 intent(in) :: vcmax25, rd25
      type(leaf_photo_params_t)            :: p
      p         = table%pft(ipft)
      p%vcmax25 = vcmax25
      p%jmax25  = table%jmax_vcmax_ratio(ipft) * vcmax25
      p%tpu25   = table%tpu_vcmax_ratio(ipft)  * vcmax25
      p%rd25    = rd25
   end function leaf_params_at_capacity

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
      !----- DIAGNOSTIC-only outputs (MEDS_IO_V01_PLAN.md section 4.7). ----------------------!
      real(wp),    optional, intent(out) :: a_net(n), ci(n), cs(n), transp(n)
      real(wp),    optional, intent(out) :: beta_stom(n), beta_nonstom(n)
      integer(ik), optional, intent(out) :: limitation(n)
      type(leaf_env_t)          :: env
      type(leaf_flux_t)         :: flux
      type(leaf_photo_params_t) :: p
      integer(ik) :: i
      do i = 1_ik, n
         env%par = par(i) ; env%leaf_temp = leaf_temp(i) ; env%vpd = vpd(i)
         env%ca = ca ; env%pressure = pressure ; env%psi_leaf = psi_leaf(i) ; env%gb = gb(i)
         if (present(psi)) then ; env%psi = psi(i) ; else ; env%psi = 0.0_wp ; end if
         p = leaf_params_at_capacity(table, pft(i), vcmax25(i), rd25(i))
         call solve_leaf_gas_exchange(env, p, table%stomatal_model, table%temp_response_form,   &
                                      table%colimitation, table%use_boundary_layer, flux)
         a_gross(i) = flux%A_gross ; gs(i) = flux%gs ; rd(i) = flux%rd
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
