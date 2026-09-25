! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_rk4_oracle -- TEST-SUPPORT ONLY: the independent cross-validation oracle for the       !
! production ARK path (meds_fast_ark) -- NOT production, reachable only from test/test_column_derivs.f90. !
! Two tiers, both retained as verification infrastructure (do not delete; see                       !
! docs/dev_plans/archive/MEDS_DRIVER_REORG_DESIGN.md §3.3):                                                  !
!                                                                                          !
!   * rk4_column_step -- the fully-EXPLICIT classical RK4 reference integrator over the pure RHS     !
!     meds_fast_time_derivs%column_derivs. It shares NO code with the split's/ARK''s backward-Euler   !
!     machinery, so agreement between it and any implicit integrator rules out a shared-bug false     !
!     pass. It is NOT a production integrator -- the fast loop is stiff (stiffness ratio ~6.5e5), so   !
!     RK4 is stable only for dt below ~2.785*tau_fast (the ~17 s plant-hydraulic mode among the        !
!     INTEGRATED reservoirs; leaf energy is diagnostic, not integrated). Used as an oracle at small dt.!
!                                                                                          !
! There is ONE oracle here, and that is the point. The first-order gamma=1 "IMEX-Euler tier" that     !
! used to sit beside it was retired (#198): it returned no reference trajectory, so it was not an      !
! oracle, and it composed column_be_stage -- the production scheme's own kernel -- so its              !
! independence was in the tableau, not in the machinery it was supposed to check. What it was          !
! genuinely useful for is degrading the production scheme in one known way and watching what fails,    !
! and that probe now lives in test_column_derivs (be_euler_step), which is the code that degrades it.  !
!==========================================================================================!
module meds_fast_rk4_oracle
   !----- The import list IS the independence claim. Nothing here reaches the implicit machinery:  !
   !      no column_be_stage, no newton_surface_solve, no advance_water_mass_full. The oracle sees   !
   !      the pure right-hand side and the state algebra, and nothing else, which is what makes       !
   !      agreement between it and an implicit scheme rule out a shared-bug false pass. Retiring the   !
   !      IMEX-Euler tier (#198) is what let meds_fast_be_stage leave this list.                        !
   use meds_kinds,            only : wp, ik
   use meds_fast_time_derivs, only : column_derivs
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t
   use meds_column_state_ops, only : state_init, state_axpy, state_accum
   implicit none
   private

   public :: rk4_column_step

contains

   !---------------------------------------------------------------------------------------!
   ! rk4_column_step -- one classical 4th-order Runge-Kutta step of the whole column state over the  !
   ! pure RHS column_derivs, with the frozen forcing `frozen` held constant across the four stages (the  !
   ! explicit part of the additive split). Commits into y_out; y is unchanged.                        !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_column_step(y, frozen, n, nsl, dt, y_out, freeze_theta)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      logical, optional,     intent(in)  :: freeze_theta   !< zero the soil-water tendency (theta held fixed):
                                                           !< makes the oracle solve the SAME reduced system the
                                                           !< ARK stepper does (soil water operator-split OUT).

      type(column_tend_t)  :: k1, k2, k3, k4
      type(column_state_t) :: y_stage
      logical              :: frz

      frz = .false. ; if (present(freeze_theta)) frz = freeze_theta

      call column_derivs(y, frozen, n, nsl, k1) ; if (frz) k1%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k1, n, nsl, y_stage)
      call column_derivs(y_stage, frozen, n, nsl, k2)
      if (frz) k2%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k2, n, nsl, y_stage)
      call column_derivs(y_stage, frozen, n, nsl, k3)
      if (frz) k3%dtheta_dt = 0.0_wp
      call state_axpy(y,          dt, k3, n, nsl, y_stage)
      call column_derivs(y_stage, frozen, n, nsl, k4)
      if (frz) k4%dtheta_dt = 0.0_wp

      !----- y_out = y + dt/6 (k1 + 2 k2 + 2 k3 + k4). -------------------------------------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt / 6.0_wp, k1, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k2, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k3, n, nsl)
      call state_accum(y_out, dt / 6.0_wp, k4, n, nsl)
   end subroutine rk4_column_step

end module meds_fast_rk4_oracle
