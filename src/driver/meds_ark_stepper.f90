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
   use meds_constants,        only : rho_h2o, tiny_num
   use meds_thermo,           only : uext_to_temp
   use meds_biophysics_types, only : n_soil_layer_max, soil_energy_column_t, energy_forcing_t, energy_flux_t
   use meds_plant_types,      only : N_HYDRO, hydro_env_t, hydro_flux_t
   use meds_column_energy,    only : soil_energy_flux
   use meds_column_hydrology, only : soil_be_single_step
   use meds_plant_hydraulics, only : solve_plant_water
   use meds_column_derivs,    only : column_state_t, column_frozen_t, column_tend_t, column_derivs, &
                                     surface_state_t, surface_frozen_t, surface_tend_t, surface_derivs
   implicit none
   private

   public :: rk4_column_step, imex_euler_column_step

contains

   !---------------------------------------------------------------------------------------!
   ! imex_euler_column_step -- one first-order IMEX-Euler step of the whole column: the stiff        !
   ! reservoirs are advanced by their UNCONDITIONALLY-STABLE backward-Euler / exact kernels (CAS     !
   ! twins BE-in-atm, soil heat + water BE-Thomas, plant hydraulics exact 2x2), driven by the        !
   ! surface sources frozen at the start-of-step state. This is the gamma = 1, single-stage member   !
   ! of the ARK family (P2 baseline): it is STABLE at the full dt_fast where the explicit RK4 oracle  !
   ! blows up, and reduces to the current operator-split coupling. The coupled-surface arrowhead      !
   ! (which additionally removes the Lie-Trotter coupling error) and the higher-order tableau are     !
   ! the P2/P3 work that builds on this. Reuses the validated production kernels -- no new numerics.  !
   !---------------------------------------------------------------------------------------!
   subroutine imex_euler_column_step(y, fro, n, nsl, dt, y_out, niter, relax)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter    !< CAS<->leaf Picard passes (default 1 = uncoupled baseline)
      real(wp),    optional, intent(in)  :: relax    !< under-relaxation of the CAS iterate (default 1)

      type(surface_state_t)      :: ys
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: sf
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      type(hydro_env_t)          :: henv
      type(hydro_flux_t)         :: hfx
      real(wp)    :: t_ground, fliq1, wmass1, wcap, ccap, gah, gaw, gac
      real(wp)    :: root_uptake(n_soil_layer_max), psi_e(n_soil_layer_max)
      real(wp)    :: theta_out(n_soil_layer_max), wface(n_soil_layer_max), drain, uptk, psi_i(N_HYDRO)
      real(wp)    :: enth1, shv1, enth_it, shv_it, rlx
      integer(ik) :: k, i, rc, np, it
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)
      rlx = 1.0_wp ; if (present(relax)) rlx = relax

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, fro%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      wcap = fro%surf%wcap ; ccap = fro%surf%ccap
      gah  = fro%surf%gah  ; gaw  = fro%surf%gaw ; gac = fro%surf%gac
      fs = fro%surf ; fs%t_ground = t_ground

      !----- CAS enthalpy + humidity: a leaf<->CAS Picard fixed point (np passes). The leaf source   !
      !      is re-evaluated at the CURRENT CAS iterate so the VPD self-limits the transpiration and  !
      !      the pair can approach saturation without over-shooting -- the coupled-surface solve that !
      !      removes the Lie-Trotter coupling error. np = 1 is the uncoupled BE baseline exactly (one !
      !      source evaluation at the start state, committed). The final-pass sf drives the soil +     !
      !      hydraulics sinks, so the water the CAS gains equals the water the soil/plant released.    !
      enth1 = y%cas_enthalpy ; shv1 = y%cas_shv ; enth_it = y%cas_enthalpy ; shv_it = y%cas_shv
      do it = 1_ik, np
         ys%cas_enthalpy = enth_it ; ys%cas_shv = shv_it ; ys%cas_co2 = y%cas_co2
         call surface_derivs(ys, fs, n, sf)
         enth1 = (wcap*y%cas_enthalpy + dt*(sf%src_enth + gah*fro%surf%enth_atm)) / (wcap + dt*gah)
         shv1  = (wcap*y%cas_shv      + dt*(sf%src_vap  + gaw*fro%surf%shv_atm )) / (wcap + dt*gaw)
         if (it < np) then
            enth_it = rlx*enth1 + (1.0_wp - rlx)*enth_it
            shv_it  = rlx*shv1  + (1.0_wp - rlx)*shv_it
         end if
      end do
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = (ccap*y%cas_co2 + dt*(fro%surf%nee_biotic + gac*fro%surf%co2_atm)) / (ccap + dt*gac)

      !----- soil-heat column: implicit BE-Thomas (soil_energy_flux). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top = sf%g_top ; eforc%geothermal = fro%geothermal
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%root_heat_sink(k) = sf%coh_qsoil * fro%soil%root_frac(k)
         eforc%w_flux(k)         = 0.0_wp
      end do
      call soil_energy_flux(se, eforc, fro%therm, fro%soil, fro%energy_opts, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- soil-water column: implicit BE-Thomas Richards (soil_be_single_step). --------------!
      rc = fro%soil%retention
      psi_e(1:nsl) = fro%psi_e(1:nsl)
      do k = 1_ik, nsl
         root_uptake(k) = sf%coh_transp * fro%soil%root_frac(k)
      end do
      call soil_be_single_step(y%theta, fro%soil, fro%hydro_opts, rc, nsl, dt, fro%q_top, psi_e,  &
                               root_uptake, theta_out, drain, uptk, wface, ok)
      y_out%theta(1:nsl) = theta_out(1:nsl)

      !----- plant hydraulics: exact frozen-2x2 matrix exponential per cohort. -------------------!
      do i = 1_ik, n
         henv%transp     = sf%transp_c(i) * fro%surf%src_frac / max(fro%nplant(i), tiny_num)
         henv%soil_psi   = fro%soil_psi_root ; henv%rhizo_cond = fro%rhizo_cond
         henv%bleaf      = fro%bleaf(i) ; henv%bsap = fro%bsap(i) ; henv%broot = fro%broot(i)
         henv%sap_area   = fro%sap_area(i) ; henv%height = fro%height(i) ; henv%leaf_area = fro%leaf_area(i)
         psi_i = y%psi(:, i)
         call solve_plant_water(henv, fro%hydro_p, fro%hydro_o, dt, psi_i, hfx)
         y_out%psi(:, i) = psi_i
      end do
   end subroutine imex_euler_column_step

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
