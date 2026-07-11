!==========================================================================================!
! meds_ark_stepper -- time-integrators for the fast-loop column state over the pure RHS            !
! meds_column_derivs (design archive/MEDS_IMEX_ARK_DESIGN.md). The overhaul replaces the split +    !
! Picard fast step with ONE additive Runge-Kutta advance; this module is where the tableaux live.   !
!                                                                                          !
! PHASE P1 (here): the fully-EXPLICIT classical RK4 reference integrator `rk4_column_step`. It is    !
! the plan's verification ORACLE (roadmap INTEG_RK4): it shares no code with the split's backward-   !
! Euler machinery, so agreement between it and any implicit integrator rules out a shared-bug false  !
! pass. It is NOT a production integrator -- the fast loop is stiff (stiffness ratio ~6.5e5), so RK4  !
! is stable only for dt below ~2.785*tau_fast (the ~17 s plant-hydraulic mode among the INTEGRATED   !
! reservoirs; leaf energy is diagnostic, not integrated). Used as an oracle at small dt.             !
!                                                                                          !
! DEFERRED here: the L-stable IMEX-ARK stage solve (P2 arrowhead + Thomas + 2x2) and the adaptive    !
! embedded-error controller (P3). This module is their home; `rk4_column_step` establishes the       !
! state-arithmetic + stage-evaluation scaffold they reuse.                                           !
!==========================================================================================!
module meds_ark_stepper
   use meds_kinds,            only : wp, ik
   use meds_biophysics_types, only : n_soil_layer_max
   use meds_plant_types,      only : N_HYDRO
   use meds_column_derivs,    only : column_state_t, column_frozen_t, column_tend_t, column_derivs
   implicit none
   private

   public :: rk4_column_step

contains

   !---------------------------------------------------------------------------------------!
   ! rk4_column_step -- one classical 4th-order Runge-Kutta step of the whole column state over the  !
   ! pure RHS column_derivs, with the frozen forcing `fro` held constant across the four stages (the  !
   ! explicit part of the additive split). Commits into y_out; y is unchanged.                        !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_column_step(y, fro, n, nsl, dt, y_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out

      type(column_tend_t)  :: k1, k2, k3, k4
      type(column_state_t) :: ys

      call column_derivs(y, fro, n, nsl, k1)
      call state_axpy(y, 0.5_wp * dt, k1, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k2)
      call state_axpy(y, 0.5_wp * dt, k2, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k3)
      call state_axpy(y,          dt, k3, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k4)

      !----- y_out = y + dt/6 (k1 + 2 k2 + 2 k3 + k4). -------------------------------------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt / 6.0_wp, k1, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k2, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k3, n, nsl)
      call state_accum(y_out, dt / 6.0_wp, k4, n, nsl)
   end subroutine rk4_column_step

   !----- ys = y + a*k  (state + a * tendency). ------------------------------------------!
   pure subroutine state_axpy(y, a, k, n, nsl, ys)
      type(column_state_t), intent(in)  :: y
      real(wp),             intent(in)  :: a
      type(column_tend_t),  intent(in)  :: k
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: ys
      integer(ik) :: j, i
      ys%cas_enthalpy = y%cas_enthalpy + a * k%d_cas_enthalpy
      ys%cas_shv      = y%cas_shv      + a * k%d_cas_shv
      ys%cas_co2      = y%cas_co2      + a * k%d_cas_co2
      ys%soil_energy  = y%soil_energy
      ys%theta        = y%theta
      do j = 1_ik, nsl
         ys%soil_energy(j) = y%soil_energy(j) + a * k%dedt(j)
         ys%theta(j)       = y%theta(j)       + a * k%dtheta_dt(j)
      end do
      allocate(ys%psi(N_HYDRO, n))
      do i = 1_ik, n
         ys%psi(:, i) = y%psi(:, i) + a * k%dpsi_dt(:, i)
      end do
   end subroutine state_axpy

   !----- copy the prognostic state (used to seed the RK combination). --------------------!
   pure subroutine state_init(y, n, nsl, ys)
      type(column_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: ys
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      ys%soil_energy  = y%soil_energy  ; ys%theta   = y%theta
      allocate(ys%psi(N_HYDRO, n))
      ys%psi(:, 1:n) = y%psi(:, 1:n)
   end subroutine state_init

   !----- ys += a*k  (accumulate a weighted tendency into a state). -----------------------!
   pure subroutine state_accum(ys, a, k, n, nsl)
      type(column_state_t), intent(inout) :: ys
      real(wp),             intent(in)    :: a
      type(column_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      integer(ik) :: j, i
      ys%cas_enthalpy = ys%cas_enthalpy + a * k%d_cas_enthalpy
      ys%cas_shv      = ys%cas_shv      + a * k%d_cas_shv
      ys%cas_co2      = ys%cas_co2      + a * k%d_cas_co2
      do j = 1_ik, nsl
         ys%soil_energy(j) = ys%soil_energy(j) + a * k%dedt(j)
         ys%theta(j)       = ys%theta(j)       + a * k%dtheta_dt(j)
      end do
      do i = 1_ik, n
         ys%psi(:, i) = ys%psi(:, i) + a * k%dpsi_dt(:, i)
      end do
   end subroutine state_accum

end module meds_ark_stepper
