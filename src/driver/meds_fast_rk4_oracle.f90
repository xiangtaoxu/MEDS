!==========================================================================================!
! meds_fast_rk4_oracle -- TEST-SUPPORT ONLY: the independent cross-validation oracle for the       !
! production ARK path (meds_fast_ark) -- NOT production, reachable only from test/test_column_derivs.f90. !
! Two tiers, both retained as verification infrastructure (do not delete; see                       !
! docs/dev_plans/MEDS_DRIVER_REORG_DESIGN.md §3.3):                                                  !
!                                                                                          !
!   * rk4_column_step -- the fully-EXPLICIT classical RK4 reference integrator over the pure RHS     !
!     meds_fast_time_derivs%column_derivs. It shares NO code with the split's/ARK''s backward-Euler   !
!     machinery, so agreement between it and any implicit integrator rules out a shared-bug false     !
!     pass. It is NOT a production integrator -- the fast loop is stiff (stiffness ratio ~6.5e5), so   !
!     RK4 is stable only for dt below ~2.785*tau_fast (the ~17 s plant-hydraulic mode among the        !
!     INTEGRATED reservoirs; leaf energy is diagnostic, not integrated). Used as an oracle at small dt.!
!   * imex_euler_column_step / adaptive_imex_march -- the first-order gamma=1 IMEX-Euler tier         !
!     (the P2 baseline), SUPERSEDED in production by the embedded-error ark2_column_step /             !
!     adaptive_ark_march (2 solves/step vs step-doubling''s 3x cost) but kept here as a simpler,        !
!     independently-derived cross-check on the shared column_be_stage/advance_hydraulics_full          !
!     building blocks (imported cross-module from meds_fast_ark, the production home).                 !
!==========================================================================================!
module meds_fast_rk4_oracle
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num
   use meds_numerics,         only : adaptive_step_update
   use meds_fast_time_derivs, only : column_derivs
   use meds_fast_types,       only : column_state_t, column_frozen_t, column_tend_t, surface_tend_t
   use meds_fast_ark,         only : column_be_stage, advance_water_mass_full, state_init
   use meds_fast_control,     only : state_wrms_grouped, default_tol_set, tol_set_t
   implicit none
   private

   public :: rk4_column_step, imex_euler_column_step, adaptive_imex_march

contains

   !---------------------------------------------------------------------------------------!
   ! adaptive_imex_march -- integrate the column from t=0 to t_end with a STEP-DOUBLING adaptive     !
   ! controller over the coupled IMEX-Euler step (the plan's P3 adaptive machinery, using the shipped !
   ! adaptive_step_update primitive; p=1 -> exponent -1/(p+1) = -1/2). Each trial compares one step   !
   ! of dt against two of dt/2; the local extrapolation (the two-half-step result) is committed and   !
   ! the WRMS of their difference drives accept/reject + the next dt. Reports the step + reject count. !
   ! The embedded-estimate ARK tableau (2nd order, one solve/stage) is the follow-on that replaces    !
   ! step-doubling's 3x cost; this establishes the adaptive-controller contract.                      !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_imex_march(y0, fro, n, nsl, t_end, rtol, dt_init, y_out, nsteps, nrej)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, rtol, dt_init
      type(column_state_t),  intent(out) :: y_out
      integer(ik),           intent(out) :: nsteps, nrej

      type(column_state_t) :: y, y_big, y_h, y_small
      type(tol_set_t)      :: tols
      real(wp)             :: t, dt, err, fac
      real(wp), parameter  :: DT_FLOOR = 1.0e-2_wp, SAFETY = 0.9_wp, FMIN = 0.2_wp, FMAX = 5.0_wp
      integer(ik), parameter :: NP = 8_ik

      tols = default_tol_set(rtol)                ! single rtol broadcast to all groups (legacy WRMS)
      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik

      do
         if (t >= t_end - tiny_num) exit
         dt = min(dt, t_end - t)
         call imex_euler_column_step(y,   fro, n, nsl, dt,          y_big,   niter=NP)
         call imex_euler_column_step(y,   fro, n, nsl, 0.5_wp*dt,   y_h,     niter=NP)
         call imex_euler_column_step(y_h, fro, n, nsl, 0.5_wp*dt,   y_small, niter=NP)
         err = state_wrms_grouped(y_big, y_small, y, n, nsl, tols)
         fac = adaptive_step_update(max(err, tiny_num), SAFETY, FMIN, FMAX)
         if (err <= 1.0_wp .or. dt <= DT_FLOOR) then
            call state_init(y_small, n, nsl, y)       ! local extrapolation: commit the finer result
            t = t + dt ; nsteps = nsteps + 1_ik
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac                        ! reject + shrink (fac < 1 since err > 1)
         end if
      end do
      call state_init(y, n, nsl, y_out)
   end subroutine adaptive_imex_march
   !---------------------------------------------------------------------------------------!
   ! imex_euler_column_step -- the COMPLETE gamma=1, first-order IMEX-Euler column integrator = one    !
   ! column_be_stage (CAS + soil) composed with the closed-form plant water-mass operator split over    !
   ! the full dt (advance_water_mass_full). It is L-stable at the full dt_fast where the explicit RK4   !
   ! oracle blows up, and is the gamma=1 member of the ARK family / the P2 baseline. The mass split is   !
   ! EXPLICIT (not folded into the BE stage) because it needs the SAME per-cohort transp_c the single    !
   ! BE stage already evaluated (weight 1.0 -- there is only one stage at gamma=1), not a separately-    !
   ! evaluated endpoint value. ark2_column_step composes the analogous 2nd-order pair (two column_be_    !
   ! stage calls, b-weighted transp), so there is ONE mass-update path. --------------------------------!
   subroutine imex_euler_column_step(y, fro, n, nsl, dt, y_out, niter)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter
      type(surface_tend_t) :: sf

      call column_be_stage(y, fro, n, nsl, dt, y_out, niter, sf_out=sf)     ! CAS + soil (mass passed through)
      call advance_water_mass_full(y, fro, n, nsl, dt, sf%transp_c(1:n), y_out)  ! closed-form mass Euler, this stage's transp
   end subroutine imex_euler_column_step
   !---------------------------------------------------------------------------------------!
   ! rk4_column_step -- one classical 4th-order Runge-Kutta step of the whole column state over the  !
   ! pure RHS column_derivs, with the frozen forcing `fro` held constant across the four stages (the  !
   ! explicit part of the additive split). Commits into y_out; y is unchanged.                        !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_column_step(y, fro, n, nsl, dt, y_out, freeze_theta)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      logical, optional,     intent(in)  :: freeze_theta   !< zero the soil-water tendency (theta held fixed):
                                                           !< makes the oracle solve the SAME reduced system the
                                                           !< ARK stepper does (soil water operator-split OUT).

      type(column_tend_t)  :: k1, k2, k3, k4
      type(column_state_t) :: ys
      logical              :: frz

      frz = .false. ; if (present(freeze_theta)) frz = freeze_theta

      call column_derivs(y, fro, n, nsl, k1) ; if (frz) k1%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k1, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k2) ; if (frz) k2%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k2, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k3) ; if (frz) k3%dtheta_dt = 0.0_wp
      call state_axpy(y,          dt, k3, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k4) ; if (frz) k4%dtheta_dt = 0.0_wp

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
      allocate(ys%leaf_water_mass(n), ys%wood_water_mass(n))
      do i = 1_ik, n
         ys%leaf_water_mass(i) = y%leaf_water_mass(i) + a * k%d_leaf_water_mass(i)
         ys%wood_water_mass(i) = y%wood_water_mass(i) + a * k%d_wood_water_mass(i)
      end do
   end subroutine state_axpy
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
         ys%leaf_water_mass(i) = ys%leaf_water_mass(i) + a * k%d_leaf_water_mass(i)
         ys%wood_water_mass(i) = ys%wood_water_mass(i) + a * k%d_wood_water_mass(i)
      end do
   end subroutine state_accum

end module meds_fast_rk4_oracle
