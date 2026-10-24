! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_heat -- the stateless soil thermal-column kernels, prognostic in internal      !
! energy per layer, temperature diagnosed by the inverter                                  !
! (docs/science/soil_biophysics.md, "The soil thermal column", eqs. (10)-(13)). Two forms  !
! share the conductive and water-enthalpy faces:                                           !
!   * soil_energy_time_deriv    -- the explicit tendency dE_k/dt [W/m3], faces at T^n;     !
!   * soil_energy_step_implicit -- the backward-Euler advance over dt, faces at T^{n+1}.   !
! set_soil_heat_env assembles their forcing from the ground heat flux and the water        !
! network's flows; soil_layer_temp is one layer's temperature from its energy.             !
!==========================================================================================!
module meds_soil_heat
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, n_soil_layer_max
   use meds_soil_params,      only : soil_profile_t, soil_params_t,                                    &
                                     ENERGY_BC_DIRICHLET
   use meds_patch_state_vars,   only : soil_energy_state_t
   use meds_therm_lib,        only : internal_energy_to_temp, internal_energy_liquid, soil_thermal_conductivity,  &
                                     soil_volumetric_heat_capacity
   use meds_linear_algebra,   only : thomas_solve
   implicit none
   private

   public :: soil_energy_step_implicit, soil_energy_time_deriv, set_soil_heat_env, soil_layer_temp
   public :: soil_heat_env_t, soil_heat_out_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: soil_heat_env_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] prescribed bottom flux, positive up (held at 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux, upward-positive like hf
                                                            !<        (set_soil_heat_env flips the network's)
      real(wp) :: water_heat_out(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy leaving each layer with water other
                                                            !<        than across a face (root uptake, drainage)
      !----- The top face's water, applied by the same upwind rule and time level as the interior      !
      !      faces (docs/science/soil_biophysics.md, "Water-thermal coupling"); no bottom face term. ---!
      real(wp) :: w_flux_top = 0.0_wp                       !< [m/s] water across the top face (<0 = infiltration in)
      real(wp) :: t_water_top = 0.0_wp                      !< [K] temperature of the infiltrating pond water
   end type soil_heat_env_t

   !----- Soil-column energy outputs + diagnostics. ----------------------------------------!
   type :: soil_heat_out_t
      real(wp) :: bottom_heat  = 0.0_wp                     !< [W/m2] heat conducted out through the bottom face
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
      real(wp) :: heat_capacity = 0.0_wp                    !< [J/m2/K] the column's, at the step's start
   end type soil_heat_out_t

contains

   !---------------------------------------------------------------------------------------!
   ! set_soil_heat_env -- the soil column's forcing from the ground heat flux and the      !
   ! water network's flows (docs/science/soil_plant_water.md §5): the enthalpy the water   !
   ! takes out of each layer other than across a face (root uptake, and the top layer's    !
   ! exfiltration, at the layer's temperature), and with the drainage at t_drain, which    !
   ! leaves through the bottom layer, so the bottom face carries none; the faces' flows,   !
   ! which the network gives down-positive (the column's are up-positive, hence the flip); !
   ! and the infiltration across the top face at t_infil, upwind like the interior faces.  !
   ! No runoff term: runoff leaves the pond. water_heat and e_drain return the two parts   !
   ! of what leaves, for the caller's ledger.                                              !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_soil_heat_env(nsl, g_top, theta, root_heat, exfil_enth, face_flux, infiltration,       &
                                     t_infil, drainage, t_drain, env, water_heat, e_drain)
      integer(ik),           intent(in)  :: nsl
      real(wp),              intent(in)  :: g_top          !< [W/m2]    net ground heat flux into the soil top
      real(wp),              intent(in)  :: theta(:)       !< [m3/m3]   soil moisture (thermal properties)
      real(wp),              intent(in)  :: root_heat(:)   !< [W/m2]    enthalpy leaving each layer with root uptake
      real(wp),              intent(in)  :: exfil_enth     !< [W/m2]    enthalpy leaving the top layer by exfiltration
      real(wp),              intent(in)  :: face_flux(:)   !< [m/s]     DOWN-positive water flux at the layer faces
      real(wp),              intent(in)  :: infiltration   !< [kg/m2/s] top-face infiltration (pond -> soil)
      real(wp),              intent(in)  :: t_infil        !< [K]       temperature of the infiltrating water
      real(wp),              intent(in)  :: drainage       !< [kg/m2/s] out of the bottom layer
      real(wp),              intent(in)  :: t_drain        !< [K]       temperature it leaves at
      type(soil_heat_env_t), intent(out) :: env
      real(wp),              intent(out) :: water_heat(:)  !< [W/m2]    root uptake and exfiltration, by layer
      real(wp),              intent(out) :: e_drain        !< [W/m2]    the drainage's
      integer(ik) :: k
      water_heat(1:nsl) = root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + exfil_enth
      e_drain = drainage * internal_energy_liquid(t_drain)
      env%g_top = g_top
      do k = 1_ik, nsl
         env%soil_water(k)     = theta(k)
         env%water_heat_out(k) = water_heat(k)
         env%w_flux(k)         = -face_flux(k)
      end do
      env%w_flux_top  = -infiltration / rho_h2o
      env%t_water_top = t_infil
      env%water_heat_out(nsl) = env%water_heat_out(nsl) + e_drain
   end subroutine set_soil_heat_env

   !----- One layer's temperature and liquid fraction from its internal energy and water content: the  !
   !      shared inverter on the layer's water mass per volume (docs/science/soil_biophysics.md, "The   !
   !      soil thermal column"). --------------------------------------------------------------------!
   elemental pure subroutine soil_layer_temp(soil_energy, theta, dry_heat_capacity, temp, fliq)
      real(wp), intent(in)  :: soil_energy         !< [J/m3]
      real(wp), intent(in)  :: theta               !< [m3/m3]
      real(wp), intent(in)  :: dry_heat_capacity   !< [J/m3/K]
      real(wp), intent(out) :: temp                !< [K]
      real(wp), intent(out) :: fliq                !< [-]
      call internal_energy_to_temp(soil_energy, theta * rho_h2o, dry_heat_capacity, temp, fliq)
   end subroutine soil_layer_temp

   !----- Each layer's temperature, liquid fraction and ice-aware conductivity at its current energy, !
   !      for the step and the tendency alike. ------------------------------------------------------!
   pure subroutine layer_thermal_state(col, forcing, profile, n, temp, fliq, kappa)
      type(soil_energy_state_t), intent(in)  :: col
      type(soil_heat_env_t),     intent(in)  :: forcing
      type(soil_profile_t),      intent(in)  :: profile
      integer(ik),               intent(in)  :: n
      real(wp),                  intent(out) :: temp(:), fliq(:)   !< [K], [-]
      real(wp),                  intent(out) :: kappa(:)           !< [W/m/K]
      integer(ik) :: k
      do k = 1_ik, n
         call soil_layer_temp(col%soil_energy(k), forcing%soil_water(k), profile%dry_heat_capacity(k), temp(k), fliq(k))
         kappa(k) = soil_thermal_conductivity(forcing%soil_water(k), fliq(k), profile%theta_sat(k),            &
                                              profile%solid_conductivity(k), profile%dry_conductivity(k))
      end do
   end subroutine layer_thermal_state

   !----- Each interior face's conductivity [W/m/K]: the two half-layers in series, the dz-weighted  !
   !      harmonic mean of the layers' conductivities. --------------------------------------------!
   pure subroutine face_conductivity(dz, kappa, n, kf)
      real(wp),    intent(in)  :: dz(:), kappa(:)
      integer(ik), intent(in)  :: n
      real(wp),    intent(out) :: kf(:)   !< face k between layers k and k+1, k < n
      integer(ik) :: k
      do k = 1_ik, n - 1_ik
         kf(k) = (dz(k) + dz(k+1)) / (dz(k) / kappa(k) + dz(k+1) / kappa(k+1))
      end do
   end subroutine face_conductivity

   !---------------------------------------------------------------------------------------!
   ! water_enthalpy_faces -- the liquid enthalpy moving water carries across each soil     !
   ! face [W/m2, positive up], upwind on the given layer temperatures (the qwf terms of    !
   ! docs/science/soil_biophysics.md eq. (10)): one rule for the step and the tendency.    !
   ! Evaporation is not a face: its enthalpy leaves in g_top's latent flux.                !
   !---------------------------------------------------------------------------------------!
   pure subroutine water_enthalpy_faces(forcing, temp, n, qwf)
      type(soil_heat_env_t), intent(in)  :: forcing
      real(wp),               intent(in)  :: temp(:)   !< layer temps to upwind on (T^n or T^{n+1})
      integer(ik),            intent(in)  :: n
      real(wp),               intent(out) :: qwf(0:)   !< [W/m2] face 0 = top ... face n = bottom
      integer(ik) :: k
      !----- Top face: infiltration carries the pond's water at t_water_top. Exfiltration leaves   !
      !      through water_heat_out, so the face never flows up; the other branch keeps the rule.  !
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
      !      water_heat_out (set_soil_heat_env). --------------------------------------------------!
      qwf(n) = 0.0_wp
   end subroutine water_enthalpy_faces

   !---------------------------------------------------------------------------------------!
   ! deep_conductance -- the Dirichlet anchor's conductance [W/m2/K], eq. (12): the bottom !
   ! node conducts to deep_temp, deep_depth - |z_node(n)| below it, through the bottom     !
   ! layer's conductivity; 0 for the Neumann (geothermal) boundary. The anchor sits below  !
   ! the column so the annual wave reflects least (eq. (13)).                              !
   !---------------------------------------------------------------------------------------!
   pure function deep_conductance(soil, profile, kappa_bot, n) result(g_deep)
      type(soil_params_t),  intent(in) :: soil
      type(soil_profile_t), intent(in) :: profile
      real(wp),             intent(in) :: kappa_bot       !< [W/m/K] bottom-layer conductivity
      integer(ik),          intent(in) :: n
      real(wp) :: g_deep, l_deep
      if (soil%thermal_bottom_bc == ENERGY_BC_DIRICHLET) then
         l_deep = soil%deep_depth - abs(profile%z_node(n))
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
   pure function bottom_heat_face(soil, forcing, g_deep, t_bot) result(hf_bot)
      type(soil_params_t),    intent(in) :: soil
      type(soil_heat_env_t), intent(in) :: forcing
      real(wp),               intent(in) :: g_deep      !< [W/m2/K] from deep_conductance
      real(wp),               intent(in) :: t_bot       !< [K] bottom-layer temperature
      real(wp) :: hf_bot
      if (soil%thermal_bottom_bc == ENERGY_BC_DIRICHLET) then
         hf_bot = -g_deep * (t_bot - soil%deep_temp)
      else
         hf_bot = forcing%geothermal
      end if
   end function bottom_heat_face

   !---------------------------------------------------------------------------------------!
   ! The soil seam: advance one patch's soil thermal column over dt (implicit backward-Euler).!
   !---------------------------------------------------------------------------------------!
   subroutine soil_energy_step_implicit(col, forcing, profile, soil, dt, flux)
      type(soil_energy_state_t),   intent(inout) :: col
      type(soil_heat_env_t),      intent(in)    :: forcing
      type(soil_profile_t),        intent(in)    :: profile
      type(soil_params_t),         intent(in)    :: soil
      real(wp),                    intent(in)    :: dt
      type(soil_heat_out_t),         intent(out)   :: flux

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max) :: t_n, fl_n, kappa, c_eff, q_src, t_new, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: e0, e1, div, g_deep

      n = profile%n_active

      !----- Invert energy -> (temp, fliq) and node-wise thermal properties at state^n: the ice-aware  !
      !      conductivity and heat capacity, on the liquid fraction the inverter returns. -----------!
      call layer_thermal_state(col, forcing, profile, n, t_n, fl_n, kappa)
      do k = 1_ik, n
         c_eff(k) = soil_volumetric_heat_capacity(forcing%soil_water(k), fl_n(k), profile%dry_heat_capacity(k))
         q_src(k) = -forcing%water_heat_out(k) / profile%dz(k)               ! [W/m3] source (sink is negative)
      end do

      !----- The bottom boundary: one conductance, shared by the matrix row and the flux. ---!
      g_deep = deep_conductance(soil, profile, kappa(n), n)

      !----- Implicit BE conduction solve for temperature^{n+1}. The Dirichlet anchor enters   !
      !      the matrix (not the residual) so it stays unconditionally stable. ---------------!
      call soil_heat_be_solve(t_n, profile%dz, profile%dz_node, kappa, c_eff, q_src, forcing%g_top,   &
                              forcing%geothermal, dt, n, t_new, kf, g_deep, soil%deep_temp)

      !----- Conservative energy update: conductive faces from T^{n+1}, advective upwind. ---!
      hf(0)  = -forcing%g_top                                             ! top face (positive up)
      do k = 1_ik, n - 1_ik
         hf(k)  = -kf(k) * (t_new(k) - t_new(k+1)) / profile%dz_node(k)      ! kf reused from soil_heat_be_solve
      end do
      hf(n)  = bottom_heat_face(soil, forcing, g_deep, t_new(n))         ! bottom boundary (positive up)
      !----- Water-enthalpy faces at T^{n+1}: the shared seam, upwind rule single-sourced. ------!
      call water_enthalpy_faces(forcing, t_new, n, qwf)

      e0 = 0.0_wp
      e1 = 0.0_wp
      do k = 1_ik, n
         e0  = e0 + col%soil_energy(k) * profile%dz(k)
         div = (hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))
         col%soil_energy(k) = col%soil_energy(k) + dt / profile%dz(k) * div + dt * q_src(k)
         e1  = e1 + col%soil_energy(k) * profile%dz(k)
      end do

      !----- Re-diagnose temperature + liquid fraction from the committed energy. ----------!
      call soil_layer_temp(col%soil_energy(1:n), forcing%soil_water(1:n), profile%dry_heat_capacity(1:n),        &
                           col%soil_temp(1:n), col%soil_fliq(1:n))

      !----- Diagnostics + closed energy budget. -------------------------------------------!
      flux%bottom_heat  = -hf(n)
      !----- The column budget carries the boundary water enthalpy too; the interior faces        !
      !      telescope out of the sum exactly, as they must. -------------------------------------!
      flux%energy_resid = (e1 - e0) - dt * (forcing%g_top - flux%bottom_heat                   &
                          - sum(forcing%water_heat_out(1:n)) + (qwf(n) - qwf(0)))
      flux%heat_capacity = sum(c_eff(1:n) * profile%dz(1:n))
   end subroutine soil_energy_step_implicit

   !---------------------------------------------------------------------------------------!
   ! soil_energy_time_deriv -- the explicit soil-heat tendency dE_k/dt [W/m3] RK45         !
   ! integrates: the implicit step's flux divergence with the faces at the current T^n.    !
   ! Commits nothing; bottom_face returns the bottom face's flux for the caller's ledger.  !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_energy_time_deriv(col, forcing, profile, soil, dedt, bottom_face)
      type(soil_energy_state_t),   intent(in)  :: col
      type(soil_heat_env_t),      intent(in)  :: forcing
      type(soil_profile_t),        intent(in)  :: profile
      type(soil_params_t),         intent(in)  :: soil
      real(wp),                    intent(out) :: dedt(n_soil_layer_max)   !< [W/m3] dE/dt per layer (0 for k>n)
      real(wp), optional,          intent(out) :: bottom_face   !< [W/m2] bottom-face heat flux, positive UP

      integer(ik) :: n, k
      real(wp), dimension(n_soil_layer_max)   :: t_n, fl_n, kappa, kf
      real(wp), dimension(0:n_soil_layer_max) :: hf, qwf
      real(wp) :: g_deep

      dedt = 0.0_wp
      n = profile%n_active

      !----- Diagnose (temp, fliq) + node conductivity at the current internal energy. --------!
      call layer_thermal_state(col, forcing, profile, n, t_n, fl_n, kappa)
      call face_conductivity(profile%dz, kappa, n, kf)

      !----- Conductive faces at the current temperature (explicit) + upwind liquid enthalpy. --!
      hf(0)  = -forcing%g_top
      do k = 1_ik, n - 1_ik
         hf(k) = -kf(k) * (t_n(k) - t_n(k+1)) / profile%dz_node(k)
      end do
      g_deep = deep_conductance(soil, profile, kappa(n), n)
      hf(n)  = bottom_heat_face(soil, forcing, g_deep, t_n(n))
      if (present(bottom_face)) bottom_face = hf(n)
      !----- Water-enthalpy faces at T^n, the same rule the implicit step applies at T^{n+1}. ------!
      call water_enthalpy_faces(forcing, t_n, n, qwf)

      !----- dE_k/dt = flux divergence + source (q_src = -water_heat_out/dz). ---------------!
      do k = 1_ik, n
         dedt(k) = ((hf(k) - hf(k-1)) + (qwf(k) - qwf(k-1))) / profile%dz(k)                        &
                   - forcing%water_heat_out(k) / profile%dz(k)
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
      call face_conductivity(dz, kappa, nzg, kf)
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

end module meds_soil_heat
