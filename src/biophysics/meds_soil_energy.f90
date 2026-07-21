!==========================================================================================!
! meds_soil_energy -- the stateless soil-thermal-column kernels (design                        !
! MEDS_ENERGY_BALANCE_DESIGN.md). Prognostic INTERNAL ENERGY / enthalpy per layer (phase-safe);  !
! temperature diagnosed via the meds_therm_lib inverter. Exposes the store in two forms:          !
!   * soil_energy_time_deriv   -- the EXPLICIT RHS dE_k/dt [W/m3] at the current state (for the    !
!                                 IMEX-ARK integrator; pure, commits nothing).                      !
!   * soil_energy_step_implicit -- the IMPLICIT backward-Euler advance over dt (mutates the store,   !
!                                  re-diagnoses temp/fliq, closes the energy budget).                 !
! Both share the same conductive-face + water-enthalpy-advection physics (evaluated at T^n for the    !
! tendency, at T^{n+1} for the step). The implicit Thomas solve reuses meds_numerics + the negative-z  !
! geometry. STATELESS: sibling temps arrive as forced inputs (the coupled fixed point is the driver's). !
!==========================================================================================!
module meds_soil_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o
   use meds_biophysics_types, only : soil_energy_column_t, energy_forcing_t, soil_thermal_params_t, &
                                     soil_params_t, energy_opts_t, energy_flux_t, n_soil_layer_max, &
                                     ENERGY_PHASE_OFF
   use meds_therm_lib,        only : uext_to_temp, internal_energy_liquid,                        &
                                     soil_thermal_cond, soil_heat_cap_vol
   use meds_numerics,         only : thomas_solve
   implicit none
   private

   public :: soil_energy_step_implicit, soil_energy_time_deriv

contains

   !---------------------------------------------------------------------------------------!
   ! The soil seam: advance one patch's soil thermal column over dt (implicit backward-Euler).!
   !---------------------------------------------------------------------------------------!
   subroutine soil_energy_step_implicit(col, forcing, therm, soil, opts, dt, flux)
      type(soil_energy_column_t),  intent(inout) :: col
      type(energy_forcing_t),      intent(in)    :: forcing
      type(soil_thermal_params_t), intent(in)    :: therm
      type(soil_params_t),         intent(in)    :: soil
      type(energy_opts_t),         intent(in)    :: opts
      real(wp),                    intent(in)    :: dt
      type(energy_flux_t),         intent(out)   :: flux

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max) :: t_n, fl_n, kappa, c_eff, q_src, t_new, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: fliq_use, wmass, e0, e1, div

      n = soil%n_active

      !----- Invert energy -> (temp, fliq) and node-wise thermal properties at state^n. ----!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o                          ! [kg/m3] water mass per volume
         call uext_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k), t_n(k), fl_n(k))
         fliq_use = fl_n(k)
         if (opts%phase_change == ENERGY_PHASE_OFF) fliq_use = 1.0_wp     ! liquid-only in P1
         kappa(k) = soil_thermal_cond(forcing%soil_water(k), fliq_use, soil%theta_sat(k),     &
                                      therm%soil_solid_conductivity(k), therm%soil_dry_conductivity(k))
         c_eff(k) = soil_heat_cap_vol(forcing%soil_water(k), fliq_use, therm%soil_dry_heat_capacity(k))
         q_src(k) = -forcing%root_heat_sink(k) / soil%dz(k)               ! [W/m3] source (sink is negative)
      end do

      !----- Implicit BE conduction solve for temperature^{n+1}. ----------------------------!
      call soil_heat_be_solve(t_n, soil%dz, soil%dz_node, kappa, c_eff, q_src, forcing%g_top,   &
                              forcing%geothermal, dt, n, t_new, kf)

      !----- Conservative energy update: conductive faces from T^{n+1}, advective upwind. ---!
      hf(0)  = -forcing%g_top                                             ! top face (positive up)
      qwf(0) = 0.0_wp
      do k = 1_ik, n - 1_ik
         hf(k)  = -kf(k) * (t_new(k) - t_new(k+1)) / soil%dz_node(k)      ! kf reused from soil_heat_be_solve
         !----- Upwind the liquid enthalpy on the SOURCE layer. w_flux is UPWARD-positive, so     !
         !      <= 0 (downward flow) draws from layer k (above face k); > 0 (upward flow) from      !
         !      k+1 (below). The sign is carried by w_flux (no extra minus), matching hf and ED2    !
         !      rk4_derivs qw_flux_g. (This convention was verified correct -- BUG6 was refuted.)   !
         if (forcing%w_flux(k) <= 0.0_wp) then
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k+1))
         end if
      end do
      hf(n)  = forcing%geothermal                                        ! bottom geothermal (positive up)
      qwf(n) = 0.0_wp

      e0 = 0.0_wp
      e1 = 0.0_wp
      do k = 1_ik, n
         e0  = e0 + col%soil_energy(k) * soil%dz(k)
         div = (hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))
         col%soil_energy(k) = col%soil_energy(k) + dt / soil%dz(k) * div + dt * q_src(k)
         e1  = e1 + col%soil_energy(k) * soil%dz(k)
      end do

      !----- Re-diagnose temperature + liquid fraction from the committed energy. ----------!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o
         call uext_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k),         &
                           col%soil_temp(k), col%soil_fliq(k))
      end do

      !----- Diagnostics + closed energy budget. -------------------------------------------!
      flux%ground_heat  = forcing%g_top
      flux%bottom_heat  = -forcing%geothermal
      flux%energy_resid = (e1 - e0) - dt * (forcing%g_top - flux%bottom_heat                   &
                          - sum(forcing%root_heat_sink(1:n)))
      flux%nsub      = 1_ik
      flux%converged = .true.
      if (opts%debug_error .and. abs(flux%energy_resid) > opts%atol * sum(c_eff(1:n)*soil%dz(1:n))) then
         error stop 'soil_energy_step_implicit: energy budget did not close'
      end if
   end subroutine soil_energy_step_implicit

   !---------------------------------------------------------------------------------------!
   ! soil_energy_time_deriv -- the EXPLICIT soil-heat RHS dE_k/dt [W/m3] at the current state, for !
   ! the IMEX-ARK fast integrator (docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md). Same flux-divergence form as  !
   ! the implicit step (conductive faces + optional water-enthalpy advection + root heat sink), but !
   ! the faces are evaluated at the CURRENT temperature T^n (the step uses the implicit T^{n+1}); as !
   ! dt -> 0 the committed BE update / dt converges to this tendency. Commits nothing.               !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_energy_time_deriv(col, forcing, therm, soil, opts, dedt)
      type(soil_energy_column_t),  intent(in)  :: col
      type(energy_forcing_t),      intent(in)  :: forcing
      type(soil_thermal_params_t), intent(in)  :: therm
      type(soil_params_t),         intent(in)  :: soil
      type(energy_opts_t),         intent(in)  :: opts
      real(wp),                    intent(out) :: dedt(n_soil_layer_max)   !< [W/m3] dE/dt per layer (0 for k>n)

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max)   :: t_n, fl_n, kappa, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: fliq_use, wmass

      dedt = 0.0_wp
      n = soil%n_active

      !----- Diagnose (temp, fliq) + node conductivity at the current internal energy. --------!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o
         call uext_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k), t_n(k), fl_n(k))
         fliq_use = fl_n(k)
         if (opts%phase_change == ENERGY_PHASE_OFF) fliq_use = 1.0_wp
         kappa(k) = soil_thermal_cond(forcing%soil_water(k), fliq_use, soil%theta_sat(k),        &
                                      therm%soil_solid_conductivity(k), therm%soil_dry_conductivity(k))
      end do
      do k = 1_ik, n - 1_ik
         kf(k) = (soil%dz(k) + soil%dz(k+1)) / (soil%dz(k) / kappa(k) + soil%dz(k+1) / kappa(k+1))
      end do

      !----- Conductive faces at the CURRENT temperature (explicit) + upwind liquid enthalpy. --!
      hf(0)  = -forcing%g_top
      qwf(0) = 0.0_wp
      do k = 1_ik, n - 1_ik
         hf(k) = -kf(k) * (t_n(k) - t_n(k+1)) / soil%dz_node(k)
         if (forcing%w_flux(k) <= 0.0_wp) then
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_n(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_n(k+1))
         end if
      end do
      hf(n)  = forcing%geothermal
      qwf(n) = 0.0_wp

      !----- dE_k/dt = flux divergence + source (q_src = -root_heat_sink/dz). ---------------!
      do k = 1_ik, n
         dedt(k) = ((hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))) / soil%dz(k)                        &
                   - forcing%root_heat_sink(k) / soil%dz(k)
      end do
   end subroutine soil_energy_time_deriv

   !---------------------------------------------------------------------------------------!
   ! Bare-array inner BE solve: implicit conduction (dz-weighted harmonic face kappa) +      !
   ! explicit source q_src, top Neumann g_top, bottom geothermal geo. Returns temperature^{n+1}. !
   ! (Advection is applied in the conservative energy update, not here.)                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_heat_be_solve(t_n, dz, dz_node, kappa, c_eff, q_src, g_top, geo, dt, nzg, t_new, kf)
      integer(ik), intent(in)  :: nzg
      real(wp),    intent(in)  :: t_n(n_soil_layer_max), dz(n_soil_layer_max)
      real(wp),    intent(in)  :: dz_node(n_soil_layer_max), kappa(n_soil_layer_max)
      real(wp),    intent(in)  :: c_eff(n_soil_layer_max), q_src(n_soil_layer_max)
      real(wp),    intent(in)  :: g_top, geo, dt
      real(wp),    intent(out) :: t_new(n_soil_layer_max)
      real(wp),    intent(out) :: kf(n_soil_layer_max)     ! series-resistor face conductivity (reused by caller)
      real(wp) :: a(n_soil_layer_max), b(n_soil_layer_max)
      real(wp) :: c(n_soil_layer_max), r(n_soil_layer_max)
      integer(ik) :: k
      do k = 1_ik, nzg - 1_ik
         kf(k) = (dz(k) + dz(k+1)) / (dz(k) / kappa(k) + dz(k+1) / kappa(k+1))     ! series-resistor face
      end do
      do k = 1_ik, nzg
         if (k >= 2_ik) then
            a(k) = -kf(k-1) / dz_node(k-1)
         else
            a(k) = 0.0_wp
         end if
         if (k <= nzg-1_ik) then
            c(k) = -kf(k) / dz_node(k)
         else
            c(k) = 0.0_wp
         end if
         b(k) = c_eff(k) * dz(k) / dt - a(k) - c(k)
         r(k) = c_eff(k) * dz(k) / dt * t_n(k) + dz(k) * q_src(k)
      end do
      r(1)   = r(1)   + g_top                                            ! top Neumann surface flux
      r(nzg) = r(nzg) + geo                                              ! bottom geothermal source
      call thomas_solve(a, b, c, r, t_new, nzg)
   end subroutine soil_heat_be_solve

end module meds_soil_energy
