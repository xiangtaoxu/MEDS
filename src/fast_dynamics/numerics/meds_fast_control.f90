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
   use meds_config,      only : CTRL_I, CTRL_PI
   use meds_fast_types,  only : column_state_t, tol_set_t, error_control_t, integrator_opts_t,     &
                                GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, &
                                GRP_SOIL_T, N_TOL_GROUP
   implicit none
   private

   public :: default_tol_set, default_error_control, state_wrms_grouped, step_control_factor


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
   ! coefficient freeze, not by the stepper (MEDS_INTEGRATOR_PARITY.md [RETIRED] sec 3d/3e). A dilution the       !
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
