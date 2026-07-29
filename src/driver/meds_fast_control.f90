!==========================================================================================!
! meds_fast_control -- the shared FAST-LOOP ERROR-CONTROL facility (MEDS_NUMERICS_SCOPING.md      !
! goal (a) / §4 / §9.3). One place owns: (1) the per-STATE-GROUP tolerance set that weights the    !
! embedded-error WRMS norm, (2) the adaptive step-size CONTROLLER (an I-controller or a PI          !
! controller over the embedded estimate), and (3) the STRICTNESS level (fixed / adaptive /          !
! strict-conservation). Previously each adaptive path carried a single scalar `rtol` + its own      !
! hard-coded per-field atols + an I-controller inline; this module unifies them so every scheme     !
! measures error and adapts identically, and a "target accuracy" is well-defined across the         !
! physically-heterogeneous column state.                                                            !
!                                                                                          !
! GROUPS. The fast column state mixes incomparable units/scales, so a single tolerance is           !
! meaningless: the WRMS normalizes each state by its GROUP's (rtol, atol). The 5 groups are the      !
! 5 physical field classes the ARK integrates -- CAS specific enthalpy, CAS specific humidity, CAS   !
! CO2, soil internal energy, and plant water potential (soil moisture theta is operator-split OUT    !
! of the ESDIRK stages, so its stage-difference is identically zero and it carries no group).        !
!                                                                                          !
! BIT-IDENTICAL DEFAULT. default_error_control(rtol) broadcasts the one `rtol` to every group and    !
! keeps the historical per-field atols + the I-controller, so state_wrms_grouped + step_control_     !
! factor reproduce the old state_wrms(rtol) + adaptive_step_update EXACTLY. Per-group tolerances,    !
! the PI controller, and the L2 strict mode are opt-in (config), so turning them off is byte-        !
! identical to before this module existed.                                                           !
!==========================================================================================!
module meds_fast_control
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : tiny_num
   use meds_numerics,    only : adaptive_step_update, clamp
   use meds_config,      only : CTRL_L0_FIXED, CTRL_L1_ADAPTIVE, CTRL_L2_STRICT, CTRL_I, CTRL_PI, &
                                meds_config_t
   use meds_fast_types,  only : column_state_t
   implicit none
   private

   public :: GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, GRP_SOIL_T, N_TOL_GROUP
   public :: tol_set_t, error_control_t
   public :: default_tol_set, default_error_control, state_wrms_grouped, step_control_factor
   public :: build_tol_set, build_error_control

   !----- The tolerance GROUPS -- one per physical field class in the fast-loop state. Groups 1-5 are   !
   !      the ARK-INTEGRATED state (what the embedded-error WRMS measures); groups 6-7 belong to the     !
   !      nested SUB-SOLVERS (soil-water Richards on theta, soil-energy on temperature) that the driver   !
   !      drives from this same set -- so ONE tolerance source governs the whole hierarchy (§8c Layer 1). !
   integer(ik), parameter :: GRP_ENTH    = 1_ik   !< CAS specific enthalpy   [J/kg]
   integer(ik), parameter :: GRP_SHV     = 2_ik   !< CAS specific humidity   [kg/kg]
   integer(ik), parameter :: GRP_CO2     = 3_ik   !< CAS CO2 mole fraction   [umol/mol]
   integer(ik), parameter :: GRP_SE      = 4_ik   !< soil internal energy    [J/m3]
   !----- Internal leaf/wood water MASS (MEDS_ED2_RK45_DESIGN.md sec 6, P2) -- replaces the retired      !
   !      GRP_PSI (plant water potential): mass, not psi, is the fast-loop prognostic state now, and      !
   !      leaf/wood get SEPARATE groups (unlike psi, which shared one) since their natural capacities      !
   !      (leaf vs. sapwood+fine-root biomass) differ enough to want independent tolerances. -------------!
   integer(ik), parameter :: GRP_LEAF_W  = 5_ik   !< leaf internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_WOOD_W  = 6_ik   !< wood internal water mass [kg/plant] (RK45 WRMS)
   integer(ik), parameter :: GRP_THETA   = 7_ik   !< soil moisture           [m3/m3] (soil-water sub-solver)
   integer(ik), parameter :: GRP_SOIL_T  = 8_ik   !< soil temperature        [K]     (soil-energy sub-solver)
   integer(ik), parameter :: N_TOL_GROUP = 8_ik

   !----- Historical per-field absolute tolerances (were `parameter`s in meds_fast_ark / the sub-solver !
   !      opts defaults; used as the group defaults so every path is byte-identical unless overridden). !
   real(wp), parameter :: ATOL_ENTH_DEF   = 5.0e1_wp    !< [J/kg]      (~0.05 K in enthalpy)
   real(wp), parameter :: ATOL_SHV_DEF    = 1.0e-6_wp   !< [kg/kg]
   real(wp), parameter :: ATOL_CO2_DEF    = 1.0e-1_wp   !< [umol/mol]
   real(wp), parameter :: ATOL_SE_DEF     = 1.0e3_wp    !< [J/m3]
   !----- [kg/plant]: a small fraction of a typical saturated tissue capacity (MVP defaults -- rtol       !
   !      does the scale-appropriate work for larger/smaller cohorts, same spirit as every other group). !
   real(wp), parameter :: ATOL_LEAF_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_WOOD_W_DEF = 1.0e-4_wp   !< [kg/plant]
   real(wp), parameter :: ATOL_THETA_DEF  = 1.0e-4_wp   !< [m3/m3] (== soil_opts_t's own default)
   real(wp), parameter :: ATOL_SOIL_T_DEF = 1.0e-2_wp   !< [K]     (== energy_opts_t's own default)

   !----- Default PI gains for a 1st-order embedded pair (accepted solution order p+1 = 2): a = 0.7/2, !
   !      b = 0.4/2 (Gustafsson 1988 / Soderlind). fac = safety*err^-a*err_prev^b. b = 0 recovers a    !
   !      pure I-controller with exponent a. -----------------------------------------------------------!
   real(wp), parameter :: PI_ALPHA_DEF = 0.35_wp
   real(wp), parameter :: PI_BETA_DEF  = 0.20_wp

   !----- Per-group (rtol, atol). The WRMS normalizes state group g by atol(g) + rtol(g)*|y|. ---------!
   type :: tol_set_t
      real(wp) :: rtol(N_TOL_GROUP) = 1.0e-3_wp
      real(wp) :: atol(N_TOL_GROUP) = [ATOL_ENTH_DEF, ATOL_SHV_DEF, ATOL_CO2_DEF, ATOL_SE_DEF,   &
                                       ATOL_LEAF_W_DEF, ATOL_WOOD_W_DEF, ATOL_THETA_DEF, ATOL_SOIL_T_DEF]
   end type tol_set_t

   !----- The bundle threaded into an adaptive march: strictness + controller + step-clamp knobs +     !
   !      PI gains + the tolerance set. -----------------------------------------------------------------!
   type :: error_control_t
      integer(ik)      :: level      = CTRL_L1_ADAPTIVE
      integer(ik)      :: controller = CTRL_I
      real(wp)         :: safety     = 0.9_wp
      real(wp)         :: fmin       = 0.2_wp
      real(wp)         :: fmax       = 5.0_wp
      real(wp)         :: pi_alpha   = PI_ALPHA_DEF
      real(wp)         :: pi_beta    = PI_BETA_DEF
      !----- Embedded-pair LOWER order (MEDS_ED2_RK45_DESIGN.md sec 6): default 1 matches ARK's       !
      !      ARS(2,2,2) 1st-order embedded estimate (byte-identical to every march before this field   !
      !      existed); Cash-Karp's RK45 sets this to 4 so step_control_factor's I-controller uses the   !
      !      correct -1/5 exponent instead of silently reusing ARK's -1/2. -------------------------!
      integer(ik)      :: p_order    = 1_ik
      type(tol_set_t)  :: tols
   end type error_control_t

contains

   !---------------------------------------------------------------------------------------!
   ! Build an error_control_t that is BYTE-IDENTICAL to the pre-module behaviour: one `rtol`     !
   ! broadcast to every group, the historical atols, the I-controller, L1. The caller overrides   !
   ! %controller / %level / per-group %tols from config to opt into PI / L2 / per-group tolerances. !
   !---------------------------------------------------------------------------------------!
   pure function default_error_control(rtol) result(ec)
      real(wp), intent(in)  :: rtol
      type(error_control_t) :: ec
      ec%tols = default_tol_set(rtol)
   end function default_error_control

   !----- The default tolerance set: one `rtol` broadcast to every group + the historical per-field    !
   !      atols. state_wrms_grouped with this set is byte-identical to the pre-module state_wrms(rtol).  !
   pure function default_tol_set(rtol) result(tols)
      real(wp), intent(in) :: rtol
      type(tol_set_t)      :: tols
      tols%rtol = rtol              ! broadcast to all groups (legacy single-rtol behaviour)
      ! tols%atol keeps its default (the historical per-field constants)
   end function default_tol_set

   !---------------------------------------------------------------------------------------!
   ! build_tol_set -- THE single tolerance source for the whole fast loop (§8c Layer 1). Each group   !
   ! is SEEDED from the setting that governs it today, so the result is byte-identical to the         !
   ! pre-unification behaviour:                                                                        !
   !   * ARK/RK45-integrated groups (enthalpy/shv/CO2/soil-energy/leaf_w/wood_w) <- [fast].ark_rtol +   !
   !     historical atols;                                                                              !
   !   * GRP_THETA   <- the [soil]   sub-solver's own (rtol, atol)  -- soil-water Richards step-doubling; !
   !   * GRP_SOIL_T  <- the [energy] sub-solver's own (rtol, atol)  -- soil-energy substepping;           !
   !   * GRP_LEAF_W/GRP_WOOD_W (MEDS_ED2_RK45_DESIGN.md sec 6, P2, replaces the retired GRP_PSI): only     !
   !     RK45 actually folds these into its embedded-error WRMS (mass is operator-split out of ARK's      !
   !     ESDIRK tableau, like psi was, via with_mass=.false.) -- seeded here regardless so the group        !
   !     exists uniformly. --------------------------------------------------------------------------------!
   !                                                                                          !
   ! ONE MASTER DIAL: when cfg%rtol_all > 0 it OVERRIDES every group's rtol, so a single number sets the  !
   ! relative accuracy of the entire hierarchy (the "target accuracy" axis goals (b)/(c) need). Left at   !
   ! its 0 default, each group keeps its own per-sub-solver value => byte-identical, even for a config    !
   ! that already customised [soil]/[energy] tolerances.                                                  !
   !---------------------------------------------------------------------------------------!
   pure function build_tol_set(cfg) result(tols)
      type(meds_config_t), intent(in) :: cfg
      type(tol_set_t)                 :: tols
      !----- ARK/RK45-integrated groups: the single ark_rtol, historical atols (atol defaults kept). --!
      tols%rtol(GRP_ENTH)   = cfg%ark_rtol
      tols%rtol(GRP_SHV)    = cfg%ark_rtol
      tols%rtol(GRP_CO2)    = cfg%ark_rtol
      tols%rtol(GRP_SE)     = cfg%ark_rtol
      tols%rtol(GRP_LEAF_W) = cfg%ark_rtol
      tols%rtol(GRP_WOOD_W) = cfg%ark_rtol
      !----- Sub-solver groups: seed from the opts that drive them today. ----------------------------!
      tols%rtol(GRP_THETA)  = cfg%soil%rtol   ; tols%atol(GRP_THETA)  = cfg%soil%atol
      tols%rtol(GRP_SOIL_T) = cfg%energy%rtol ; tols%atol(GRP_SOIL_T) = cfg%energy%atol
      !----- The one master accuracy dial (0 => unset => keep the per-group values above). ------------!
      if (cfg%rtol_all > 0.0_wp) tols%rtol = cfg%rtol_all
      !----- ...and its ABSOLUTE companion. The WRMS denominator is atol + rtol*|y|, so rtol_all alone   !
      !      SATURATES once atol dominates: measured on the split path, rtol 1e-3 -> 1e-6 raises the      !
      !      soil-water error estimate only ~4x (4e-4 -> 1e-4 denominator) -- never enough to force a     !
      !      substep, while scaling BOTH does. atol is dimensional and differs per group, so it scales    !
      !      rather than broadcasts. The default 1.0 is an exact IEEE identity => byte-identical. --------!
      tols%atol = tols%atol * cfg%atol_scale
   end function build_tol_set

   !----- The full error-control bundle from config: unified tolerances + controller + strictness. ----!
   pure function build_error_control(cfg) result(ec)
      type(meds_config_t), intent(in) :: cfg
      type(error_control_t)           :: ec
      ec%tols       = build_tol_set(cfg)
      ec%controller = cfg%step_controller
      ec%level      = cfg%error_level
   end function build_error_control

   !---------------------------------------------------------------------------------------!
   ! Grouped WRMS error norm of (a - b), each state normalized by its group's atol + rtol*|y_ref|.  !
   !                                                                                          !
   ! ONE NORM OVER THE WHOLE COLUMN STATE -- no per-caller opt-outs. Every prognostic field of      !
   ! column_state_t contributes: the three CAS twins, soil internal energy, soil MOISTURE, and the   !
   ! leaf/wood internal water MASS. A scheme does not get to declare a state uninteresting.          !
   !                                                                                          !
   ! This used to carry `with_mass` / `with_theta` switches, because on the ARK path both theta and   !
   ! mass are operator-split OUT of the ESDIRK tableau -- passed through column_be_stage unchanged,   !
   ! then advanced once over the full step -- so their embedded differences are STRUCTURALLY ZERO      !
   ! there, and summing them adds (nsl + 2n) exact zeros that still increment cnt and divide the norm  !
   ! down. That dilution is real: it makes the ARK march run looser than its stated ark_rtol.           !
   !                                                                                          !
   ! It is nonetheless the right trade, for two reasons. (1) A switch that says "do not measure this   !
   ! state" is indistinguishable, at the call site, from a switch that says "this state cannot move",   !
   ! and the two were in fact confused: RK45 integrates theta inside its tableau and the norm silently  !
   ! discarded it, leaving an integrated state with no error control at all. (2) Measured, the cost is  !
   ! not there: on b4_stand_summer the ARK's CAS-T RMSE moves by 0-1% at every dt_fast while its        !
   ! sub-step count FALLS 9-15%, because at production dt_fast the error is dominated by the Category-0 !
   ! coefficient freeze, not by the stepper (MEDS_INTEGRATOR_PARITY.md sec 3d/3e). A dilution the       !
   ! measurement cannot find is not worth a configuration axis.                                         !
   !                                                                                          !
   ! Note that soil moisture is error-CONTROLLED on every path regardless of this norm: split and the   !
   ! ARK take theta wholly from column_hydrology_flux, whose own adaptive step-doubling is driven by    !
   ! the SAME GRP_THETA tolerances build_tol_set seeds here. What this norm adds is control for the one !
   ! scheme (RK45) that took theta out of that solver and into its own stages.                          !
   !---------------------------------------------------------------------------------------!
   pure function state_wrms_grouped(a, b, y_ref, n, nsl, tols) result(err)
      type(column_state_t), intent(in) :: a, b, y_ref
      integer(ik),          intent(in) :: n, nsl
      type(tol_set_t),      intent(in) :: tols
      real(wp)    :: err, s
      integer(ik) :: k, i, cnt
      s = 0.0_wp ; cnt = 0_ik
      s = s + ((a%cas_enthalpy - b%cas_enthalpy)                                                &
               / (tols%atol(GRP_ENTH) + tols%rtol(GRP_ENTH)*abs(y_ref%cas_enthalpy)))**2 ; cnt = cnt + 1_ik
      s = s + ((a%cas_shv - b%cas_shv)                                                          &
               / (tols%atol(GRP_SHV) + tols%rtol(GRP_SHV)*abs(y_ref%cas_shv)))**2 ; cnt = cnt + 1_ik
      s = s + ((a%cas_co2 - b%cas_co2)                                                          &
               / (tols%atol(GRP_CO2) + tols%rtol(GRP_CO2)*abs(y_ref%cas_co2)))**2 ; cnt = cnt + 1_ik
      do k = 1_ik, nsl
         s = s + ((a%soil_energy(k) - b%soil_energy(k))                                         &
                  / (tols%atol(GRP_SE) + tols%rtol(GRP_SE)*abs(y_ref%soil_energy(k))))**2
         cnt = cnt + 1_ik
      end do
      do k = 1_ik, nsl
         s = s + ((a%theta(k) - b%theta(k))                                                      &
                  / (tols%atol(GRP_THETA) + tols%rtol(GRP_THETA)*abs(y_ref%theta(k))))**2
         cnt = cnt + 1_ik
      end do
      do i = 1_ik, n
         s = s + ((a%leaf_water_mass(i) - b%leaf_water_mass(i))                                  &
                  / (tols%atol(GRP_LEAF_W) + tols%rtol(GRP_LEAF_W)*abs(y_ref%leaf_water_mass(i))))**2
         s = s + ((a%wood_water_mass(i) - b%wood_water_mass(i))                                  &
                  / (tols%atol(GRP_WOOD_W) + tols%rtol(GRP_WOOD_W)*abs(y_ref%wood_water_mass(i))))**2
         cnt = cnt + 2_ik
      end do
      err = sqrt(s / real(cnt, wp))
   end function state_wrms_grouped

   !---------------------------------------------------------------------------------------!
   ! Step-size FACTOR from the (normalized) embedded error `err` and the PREVIOUS accepted step's   !
   ! error `err_prev`. err is the WRMS with tolerance baked into the denominators, so err <= 1 means !
   ! "within tolerance". Dispatch on the controller:                                                  !
   !   * CTRL_I  (or the very first step, err_prev <= 0): the elementary integral controller          !
   !     `adaptive_step_update` = clamp(safety*err^-1/2, fmin, fmax) -- BYTE-IDENTICAL to the legacy   !
   !     inline controller.                                                                            !
   !   * CTRL_PI: Gustafsson's PI, fac = clamp(safety*err^-alpha * err_prev^beta, fmin, fmax). Using   !
   !     the previous error damps the step oscillation the I-controller shows on stiff transients.     !
   ! err/err_prev are floored at tiny_num so a perfectly-converged step gives fmax, never 0^-a = Inf.  !
   !---------------------------------------------------------------------------------------!
   pure function step_control_factor(err, err_prev, ec) result(fac)
      real(wp),              intent(in) :: err, err_prev
      type(error_control_t), intent(in) :: ec
      real(wp) :: fac, e, ep
      if (ec%controller /= CTRL_PI .or. err_prev <= 0.0_wp) then
         fac = adaptive_step_update(max(err, tiny_num), ec%safety, ec%fmin, ec%fmax, ec%p_order)
      else
         e  = max(err,      tiny_num)
         ep = max(err_prev, tiny_num)
         fac = clamp(ec%safety * e**(-ec%pi_alpha) * ep**(ec%pi_beta), ec%fmin, ec%fmax)
      end if
   end function step_control_factor

end module meds_fast_control
