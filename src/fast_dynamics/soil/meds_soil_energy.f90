! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_energy -- the stateless soil thermal-column kernels, prognostic in internal    !
! energy per layer, temperature diagnosed by the inverter                                  !
! (docs/science/soil_biophysics.md, "The soil thermal column", eqs. (10)-(13)). Two forms  !
! share the conductive and water-enthalpy faces:                                           !
!   * soil_energy_time_deriv    -- the explicit tendency dE_k/dt [W/m3], faces at T^n;     !
!   * soil_energy_step_implicit -- the backward-Euler advance over dt, faces at T^{n+1}.   !
!==========================================================================================!
module meds_soil_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, n_soil_layer_max
   use meds_soil_types, only : energy_forcing_t, energy_flux_t
   use meds_soil_params,      only : soil_thermal_params_t, soil_params_t, energy_opts_t,             &
                                     ENERGY_BC_DIRICHLET
   use meds_patch_state_vars,   only : soil_energy_column_t
   use meds_therm_lib,        only : internal_energy_to_temp, internal_energy_liquid, soil_thermal_conductivity,  &
                                     soil_volumetric_heat_capacity
   use meds_linear_algebra,   only : thomas_solve
   implicit none
   private

   public :: soil_energy_step_implicit, soil_energy_time_deriv

contains

   !---------------------------------------------------------------------------------------!
   ! water_enthalpy_faces -- the liquid enthalpy moving water carries across each soil     !
   ! face [W/m2, positive up], upwind on the given layer temperatures (the qwf terms of    !
   ! docs/science/soil_biophysics.md eq. (10)): one rule for the step and the tendency.    !
   ! Evaporation is not a face: its enthalpy leaves in g_top's latent flux.                !
   !---------------------------------------------------------------------------------------!
   pure subroutine water_enthalpy_faces(forcing, temp, n, qwf)
      type(energy_forcing_t), intent(in)  :: forcing
      real(wp),               intent(in)  :: temp(:)   !< layer temps to upwind on (T^n or T^{n+1})
      integer(ik),            intent(in)  :: n
      real(wp),               intent(out) :: qwf(0:)   !< [W/m2] face 0 = top ... face n = bottom
      integer(ik) :: k
      !----- Top face: infiltration (flow down) carries the pond's water at t_water_top; flow up   !
      !      draws from layer 1. -----------------------------------------------------------------!
      if (forcing%w_flux_top <= 0.0_wp) then
         qwf(0) = forcing%w_flux_top * rho_h2o * internal_energy_liquid(forcing%t_water_top)
      else
         qwf(0) = forcing%w_flux_top * rho_h2o * internal_energy_liquid(temp(1))
      end if
      do k = 1_ik, n - 1_ik
         if (forcing%w_flux(k) <= 0.0_wp) then
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(temp(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(temp(k+1))
         end if
      end do
      !----- Bottom face: no advective term; drainage enthalpy leaves through the bottom layer's   !
      !      root_heat_sink (assemble_soil_energy_forcing). ---------------------------------------!
      qwf(n) = 0.0_wp
   end subroutine water_enthalpy_faces

   !---------------------------------------------------------------------------------------!
   ! deep_conductance -- the Dirichlet anchor's conductance [W/m2/K], eq. (12): the bottom !
   ! node conducts to deep_temp, deep_depth - |z_node(n)| below it, through the bottom     !
   ! layer's conductivity; 0 for the Neumann (geothermal) boundary. The anchor sits below  !
   ! the column so the annual wave reflects least (eq. (13)).                              !
   !---------------------------------------------------------------------------------------!
   pure function deep_conductance(opts, soil, kappa_bot, n) result(g_deep)
      type(energy_opts_t), intent(in) :: opts
      type(soil_params_t), intent(in) :: soil
      real(wp),            intent(in) :: kappa_bot        !< [W/m/K] bottom-layer conductivity
      integer(ik),         intent(in) :: n
      real(wp) :: g_deep, l_deep
      if (opts%thermal_bottom_bc == ENERGY_BC_DIRICHLET) then
         l_deep = opts%deep_depth - abs(soil%z_node(n))
         g_deep = kappa_bot / l_deep
      else
         g_deep = 0.0_wp
      end if
   end function deep_conductance

   !---------------------------------------------------------------------------------------!
   ! bottom_heat_face -- the heat flux across the bottom face [W/m2, positive up]: the     !
   ! Dirichlet anchor, eq. (12), or the prescribed geothermal flux. One expression for     !
   ! both time levels, so the implicit and explicit paths agree.                           !
   !---------------------------------------------------------------------------------------!
   pure function bottom_heat_face(opts, forcing, g_deep, t_bot) result(hf_bot)
      type(energy_opts_t),    intent(in) :: opts
      type(energy_forcing_t), intent(in) :: forcing
      real(wp),               intent(in) :: g_deep      !< [W/m2/K] from deep_conductance
      real(wp),               intent(in) :: t_bot       !< [K] bottom-layer temperature
      real(wp) :: hf_bot
      if (opts%thermal_bottom_bc == ENERGY_BC_DIRICHLET) then
         hf_bot = -g_deep * (t_bot - opts%deep_temp)
      else
         hf_bot = forcing%geothermal
      end if
   end function bottom_heat_face

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
      real(wp) :: wmass, e0, e1, div, g_deep

      n = soil%n_active

      !----- Invert energy -> (temp, fliq) and node-wise thermal properties at state^n. ----!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o                          ! [kg/m3] water mass per volume
         call internal_energy_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k), t_n(k), fl_n(k))
         !----- Ice-aware conductivity and heat capacity, on the liquid fraction the inverter returns. !
         kappa(k) = soil_thermal_conductivity(forcing%soil_water(k), fl_n(k), soil%theta_sat(k),      &
                                      therm%soil_solid_conductivity(k), therm%soil_dry_conductivity(k))
         c_eff(k) = soil_volumetric_heat_capacity(forcing%soil_water(k), fl_n(k), therm%soil_dry_heat_capacity(k))
         q_src(k) = -forcing%root_heat_sink(k) / soil%dz(k)               ! [W/m3] source (sink is negative)
      end do

      !----- The bottom boundary: one conductance, shared by the matrix row and the flux. ---!
      g_deep = deep_conductance(opts, soil, kappa(n), n)

      !----- Implicit BE conduction solve for temperature^{n+1}. The Dirichlet anchor enters   !
      !      the matrix (not the residual) so it stays unconditionally stable. ---------------!
      call soil_heat_be_solve(t_n, soil%dz, soil%dz_node, kappa, c_eff, q_src, forcing%g_top,   &
                              forcing%geothermal, dt, n, t_new, kf, g_deep, opts%deep_temp)

      !----- Conservative energy update: conductive faces from T^{n+1}, advective upwind. ---!
      hf(0)  = -forcing%g_top                                             ! top face (positive up)
      do k = 1_ik, n - 1_ik
         hf(k)  = -kf(k) * (t_new(k) - t_new(k+1)) / soil%dz_node(k)      ! kf reused from soil_heat_be_solve
      end do
      hf(n)  = bottom_heat_face(opts, forcing, g_deep, t_new(n))         ! bottom boundary (positive up)
      !----- Water-enthalpy faces at T^{n+1}: the shared seam, upwind rule single-sourced. ------!
      call water_enthalpy_faces(forcing, t_new, n, qwf)

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
         call internal_energy_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k),         &
                           col%soil_temp(k), col%soil_fliq(k))
      end do

      !----- Diagnostics + closed energy budget. -------------------------------------------!
      flux%bottom_heat  = -hf(n)
      !----- The column budget carries the boundary water enthalpy too; the interior faces        !
      !      telescope out of the sum exactly, as they must. -------------------------------------!
      flux%energy_resid = (e1 - e0) - dt * (forcing%g_top - flux%bottom_heat                   &
                          - sum(forcing%root_heat_sink(1:n)) + (qwf(n) - qwf(0)))
      if (opts%debug_error .and. abs(flux%energy_resid) > opts%soil_energy_budget_atol * sum(c_eff(1:n)*soil%dz(1:n))) then
         error stop 'soil_energy_step_implicit: energy budget did not close'
      end if
   end subroutine soil_energy_step_implicit

   !---------------------------------------------------------------------------------------!
   ! soil_energy_time_deriv -- the explicit soil-heat tendency dE_k/dt [W/m3] RK45         !
   ! integrates: the implicit step's flux divergence with the faces at the current T^n.    !
   ! Commits nothing; bottom_face returns the bottom face's flux for the caller's ledger.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_energy_time_deriv(col, forcing, therm, soil, opts, dedt, bottom_face)
      type(soil_energy_column_t),  intent(in)  :: col
      type(energy_forcing_t),      intent(in)  :: forcing
      type(soil_thermal_params_t), intent(in)  :: therm
      type(soil_params_t),         intent(in)  :: soil
      type(energy_opts_t),         intent(in)  :: opts
      real(wp),                    intent(out) :: dedt(n_soil_layer_max)   !< [W/m3] dE/dt per layer (0 for k>n)
      real(wp), optional,          intent(out) :: bottom_face   !< [W/m2] bottom-face heat flux, positive UP

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max)   :: t_n, fl_n, kappa, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: wmass, g_deep

      dedt = 0.0_wp
      n = soil%n_active

      !----- Diagnose (temp, fliq) + node conductivity at the current internal energy. --------!
      do k = 1_ik, n
         wmass = forcing%soil_water(k) * rho_h2o
         call internal_energy_to_temp(col%soil_energy(k), wmass, therm%soil_dry_heat_capacity(k), t_n(k), fl_n(k))
         kappa(k) = soil_thermal_conductivity(forcing%soil_water(k), fl_n(k), soil%theta_sat(k),         &
                                      therm%soil_solid_conductivity(k), therm%soil_dry_conductivity(k))
      end do
      do k = 1_ik, n - 1_ik
         kf(k) = (soil%dz(k) + soil%dz(k+1)) / (soil%dz(k) / kappa(k) + soil%dz(k+1) / kappa(k+1))
      end do

      !----- Conductive faces at the current temperature (explicit) + upwind liquid enthalpy. --!
      hf(0)  = -forcing%g_top
      do k = 1_ik, n - 1_ik
         hf(k) = -kf(k) * (t_n(k) - t_n(k+1)) / soil%dz_node(k)
      end do
      g_deep = deep_conductance(opts, soil, kappa(n), n)
      hf(n)  = bottom_heat_face(opts, forcing, g_deep, t_n(n))
      if (present(bottom_face)) bottom_face = hf(n)
      !----- Water-enthalpy faces at T^n, the same rule the implicit step applies at T^{n+1}. ------!
      call water_enthalpy_faces(forcing, t_n, n, qwf)

      !----- dE_k/dt = flux divergence + source (q_src = -root_heat_sink/dz). ---------------!
      do k = 1_ik, n
         dedt(k) = ((hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))) / soil%dz(k)                        &
                   - forcing%root_heat_sink(k) / soil%dz(k)
      end do
   end subroutine soil_energy_time_deriv

   !---------------------------------------------------------------------------------------!
   ! Bare-array inner BE solve: implicit conduction (dz-weighted harmonic face kappa) +      !
   ! explicit source q_src, top Neumann g_top, and a bottom boundary that is either a prescribed  !
   ! flux `geo` (g_deep = 0) or conduction to a fixed deep temperature (g_deep > 0). Returns       !
   ! temperature^{n+1}; advection is applied in the conservative energy update.                   !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_heat_be_solve(t_n, dz, dz_node, kappa, c_eff, q_src, g_top, geo, dt, nzg, t_new, kf,  &
                                      g_deep, t_deep)
      integer(ik), intent(in)  :: nzg
      real(wp),    intent(in)  :: t_n(n_soil_layer_max), dz(n_soil_layer_max)
      real(wp),    intent(in)  :: dz_node(n_soil_layer_max), kappa(n_soil_layer_max)
      real(wp),    intent(in)  :: c_eff(n_soil_layer_max), q_src(n_soil_layer_max)
      real(wp),    intent(in)  :: g_top, geo, dt
      real(wp),    intent(in)  :: g_deep     !< [W/m2/K] bottom Dirichlet conductance (0 = Neumann)
      real(wp),    intent(in)  :: t_deep     !< [K] bottom Dirichlet anchor (unread when g_deep = 0)
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
      !----- Bottom boundary. With g_deep = 0 (Neumann) the prescribed flux is a pure source term   !
      !      and the diagonal is untouched; the Dirichlet anchor adds g_deep to the diagonal and    !
      !      g_deep*t_deep to the source (eq. (12)). -----------------------------------------------!
      b(nzg) = b(nzg) + g_deep
      r(nzg) = r(nzg) + geo + g_deep * t_deep
      call thomas_solve(a, b, c, r, t_new, nzg)
   end subroutine soil_heat_be_solve

end module meds_soil_energy
