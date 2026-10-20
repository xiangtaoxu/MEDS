! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_control -- the shared FAST-LOOP ERROR-CONTROL facility. One place owns: (1) the      !
! per-STATE-GROUP tolerance set that weights the embedded-error WRMS norm, (2) the adaptive      !
! step-size CONTROLLER (an I-controller or a PI controller over the embedded estimate), and (3)  !
! the STRICTNESS level (fixed / adaptive / strict-conservation). Every scheme measures error and !
! adapts through it, so a "target accuracy" means the same thing on both and across the          !
! physically heterogeneous column state: docs/science/numerical_scheme.md eqs. (3)-(4).          !
!                                                                                                !
! GROUPS. The fast column state mixes incomparable units and scales, so a single tolerance is    !
! meaningless: the WRMS normalizes each state by its GROUP's (rtol, atol). The seven groups      !
! (meds_fast_types) are CAS specific enthalpy, CAS specific humidity, CAS CO2, soil internal     !
! energy, leaf and wood internal water mass, and soil moisture theta; the canopy films and the   !
! pond carry no group (meds_column_state_ops%state_entry_rules).                                 !
!                                                                                                !
! DEFAULTS. default_error_control(rtol) gives every group the one `rtol`, the per-field default  !
! atols and the I-controller. Per-group tolerances, the PI controller and the L2 strict mode are !
! opt-in through config (meds_fast_config%build_error_control).                                  !
!==========================================================================================!
module meds_fast_control
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : tiny_num
   use meds_numerics,    only : adaptive_step_update, clamp
   use meds_config,      only : CTRL_I, CTRL_PI
   use meds_column_state_ops, only : state_to_array, state_entry_rules
   use meds_fast_types,  only : column_state_t, tol_set_t, error_control_t, integrator_opts_t,     &
                                GRP_ENTH, GRP_SHV, GRP_CO2, GRP_SE, GRP_LEAF_W, GRP_WOOD_W, GRP_THETA, &
                                N_TOL_GROUP
   implicit none
   private

   public :: default_tol_set, default_error_control, state_wrms_grouped, step_control_factor


contains

   !---------------------------------------------------------------------------------------!
   ! The default error_control_t: one `rtol` for every group, the per-field default atols, the       !
   ! I-controller, level L1. The caller overrides %controller / %level / per-group %tols from config !
   ! to opt into PI / L2 / per-group tolerances.                                                     !
   !---------------------------------------------------------------------------------------!
   pure function default_error_control(rtol) result(ec)
      real(wp), intent(in)  :: rtol
      type(error_control_t) :: ec
      ec%tols = default_tol_set(rtol)
   end function default_error_control

   !----- The default tolerance set: one `rtol` for every group and the per-field default atols    !
   !      (tol_set_t's initial values, meds_fast_types). ------------------------------------------!
   pure function default_tol_set(rtol) result(tols)
      real(wp), intent(in) :: rtol
      type(tol_set_t)      :: tols
      tols%rtol = rtol              ! the same rtol for every group
      ! tols%atol keeps its default (the per-field constants)
   end function default_tol_set

   !---------------------------------------------------------------------------------------!
   ! Grouped WRMS error norm of (a - b), each state normalized by its group's atol + rtol*|y_ref|:    !
   ! docs/science/numerical_scheme.md eq. (3).                                                        !
   !                                                                                                  !
   ! ONE NORM OVER THE WHOLE COLUMN STATE -- no per-caller opt-outs. Every prognostic field with a    !
   ! group contributes: the three CAS twins, soil internal energy, soil MOISTURE, and the leaf/wood   !
   ! internal water MASS. A scheme does not get to declare a state uninteresting.                     !
   !                                                                                                  !
   ! On the ARK path theta and the plant water ride operator-split maps outside the ESDIRK tableau,   !
   ! so their embedded differences are STRUCTURALLY ZERO there, and summing them adds (nsl + 2n)      !
   ! exact zeros that still increment cnt and divide the norm down: the ARK march runs slightly       !
   ! looser than its stated rtol. That is the right trade, for two reasons. (1) A switch that     !
   ! says "do not measure this state" is indistinguishable, at the call site, from a switch that      !
   ! says "this state cannot move": RK45 integrates theta inside its tableau, and such a switch       !
   ! would leave an integrated state with no error control at all. (2) Measured on a summer stand,    !
   ! the accuracy cost is under 1% while the ARK's sub-step count FALLS 9-15%, because at production  !
   ! dt_fast the error is dominated by the Category-0 coefficient freeze, not by the stepper          !
   ! (docs/science/numerical_scheme.md §5 and §7 item 3).                                             !
   !                                                                                                  !
   ! Soil moisture is error-CONTROLLED on every path regardless of this norm: both schemes take theta !
   ! from the soil-plant water network, whose own sub-step control is driven by the SAME GRP_THETA    !
   ! tolerances build_tol_set seeds here; within a step theta moves at the network's steady rate,     !
   ! which every tableau integrates exactly.                                                          !
   !---------------------------------------------------------------------------------------!
   pure function state_wrms_grouped(a, b, y_ref, n, nsl, tols) result(err)
      type(column_state_t), intent(in) :: a, b, y_ref
      integer(ik),          intent(in) :: n, nsl
      type(tol_set_t),      intent(in) :: tols
      real(wp)    :: err, s
      real(wp)    :: va(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n), vr(5_ik + 2_ik*nsl + 4_ik*n)
      integer(ik) :: j, grp, cnt
      logical     :: counted
      call state_to_array(a, n, nsl, va)
      call state_to_array(b, n, nsl, vb)
      call state_to_array(y_ref, n, nsl, vr)
      s = 0.0_wp ; cnt = 0_ik
      do j = 1_ik, size(va, kind=ik)            ! the layout's order is the order the terms are added
         call state_entry_rules(j, n, nsl, grp, counted)
         if (grp == 0_ik) cycle
         s = s + ((va(j) - vb(j)) / (tols%atol(grp) + tols%rtol(grp)*abs(vr(j))))**2
         cnt = cnt + 1_ik
      end do
      err = sqrt(s / real(cnt, wp))
   end function state_wrms_grouped

   !---------------------------------------------------------------------------------------!
   ! Step-size FACTOR from the (normalized) embedded error `err` and the PREVIOUS accepted step's     !
   ! error `err_prev`: docs/science/numerical_scheme.md eq. (4). err is the WRMS with tolerance       !
   ! baked into the denominators, so err <= 1 means "within tolerance". Dispatch on the controller:   !
   !   * CTRL_I  (or the very first step, err_prev <= 0): the elementary integral controller          !
   !     `adaptive_step_update` = clamp(safety*err^(-1/(p+1)), fmin, fmax), p = ec%p_order the        !
   !     order of the embedded estimate.                                                              !
   !   * CTRL_PI: Gustafsson's PI, fac = clamp(safety*err^-alpha * err_prev^beta, fmin, fmax). Using  !
   !     the previous error damps the step oscillation the I-controller shows on stiff transients.    !
   ! err/err_prev are floored at tiny_num so a perfectly-converged step gives fmax, never 0^-a = Inf. !
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
