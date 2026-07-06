!==========================================================================================!
! meds_column_energy -- THE stateless seam of the 1-D soil THERMAL column (design              !
! MEDS_ENERGY_BALANCE_DESIGN.md, 4c/6). The heat twin of meds_column_hydrology: same negative-z  !
! geometry (soil_params_t), same meds_soil_solver Thomas sweep. It advances one patch's           !
! prognostic soil INTERNAL ENERGY over a step dt by implicit backward-Euler heat diffusion         !
! (dz-weighted harmonic-mean interface conductivity), an explicit upwind advective-heat term       !
! carried by the hydrology's inter-layer water flux, and a root-uptake enthalpy sink; it returns    !
! the diagnosed soil temperature + a machine-precision energy budget. Freeze/thaw is automatic in   !
! the meds_thermo inverter (P1 runs liquid-only, phase_change OFF).                                  !
!                                                                                          !
! P1: one implicit step per call, no adaptive substepping (deferred with the coupling). The inner    !
! `soil_heat_be_step` is bare-array + device-eligible (the growth_step precedent).                    !
!==========================================================================================!
module meds_column_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num
   use meds_biophysics_types, only : soil_energy_column_t, energy_forcing_t, soil_thermal_params_t, &
                                     soil_params_t, energy_opts_t, energy_flux_t, n_soil_layer_max, &
                                     ENERGY_PHASE_OFF
   use meds_thermo,           only : uext_to_temp, internal_energy_liquid
   use meds_soil_thermal,     only : soil_thermal_cond, soil_heat_cap_vol
   use meds_soil_solver,      only : thomas_solve
   implicit none
   private

   public :: soil_energy_flux

contains

   !---------------------------------------------------------------------------------------!
   ! The soil seam: advance one patch's soil thermal column over dt.                        !
   !---------------------------------------------------------------------------------------!
   subroutine soil_energy_flux(col, forcing, therm, soil, opts, dt, flux)
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
      call soil_heat_be_step(t_n, soil%dz, soil%dz_node, kappa, c_eff, q_src, forcing%g_top,   &
                             forcing%geothermal, dt, n, t_new)

      !----- Conservative energy update: conductive faces from T^{n+1}, advective upwind. ---!
      hf(0)  = -forcing%g_top                                             ! top face (positive up)
      qwf(0) = 0.0_wp
      do k = 1_ik, n - 1_ik
         kf(k)  = (soil%dz(k) + soil%dz(k+1)) / (soil%dz(k) / kappa(k) + soil%dz(k+1) / kappa(k+1))
         hf(k)  = -kf(k) * (t_new(k) - t_new(k+1)) / soil%dz_node(k)
         if (forcing%w_flux(k) <= 0.0_wp) then                           ! upwind on the source layer
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k))
         else
            qwf(k) = forcing%w_flux(k) * rho_h2o * internal_energy_liquid(t_new(k+1))
         end if
      end do
      hf(n)  = forcing%geothermal                                        ! bottom geothermal (positive up)
      qwf(n) = 0.0_wp

      e0 = 0.0_wp ; e1 = 0.0_wp
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
         error stop 'soil_energy_flux: energy budget did not close'
      end if
   end subroutine soil_energy_flux

   !---------------------------------------------------------------------------------------!
   ! Bare-array device-eligible inner step: implicit BE conduction (dz-weighted harmonic     !
   ! face kappa) + explicit source q_src, top Neumann g_top, bottom geothermal geo. Returns   !
   ! temperature^{n+1}. (Advection is applied in the conservative energy update, not here.)   !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_heat_be_step(t_n, dz, dz_node, kappa, c_eff, q_src, g_top, geo, dt, nzg, t_new)
      integer(ik), intent(in)  :: nzg
      real(wp),    intent(in)  :: t_n(n_soil_layer_max), dz(n_soil_layer_max)
      real(wp),    intent(in)  :: dz_node(n_soil_layer_max), kappa(n_soil_layer_max)
      real(wp),    intent(in)  :: c_eff(n_soil_layer_max), q_src(n_soil_layer_max)
      real(wp),    intent(in)  :: g_top, geo, dt
      real(wp),    intent(out) :: t_new(n_soil_layer_max)
      real(wp) :: kf(n_soil_layer_max), a(n_soil_layer_max), b(n_soil_layer_max)
      real(wp) :: c(n_soil_layer_max), r(n_soil_layer_max)
      integer(ik) :: k
      do k = 1_ik, nzg - 1_ik
         kf(k) = (dz(k) + dz(k+1)) / (dz(k) / kappa(k) + dz(k+1) / kappa(k+1))     ! series-resistor face
      end do
      do k = 1_ik, nzg
         if (k >= 2_ik) then ; a(k) = -kf(k-1) / dz_node(k-1) ; else ; a(k) = 0.0_wp ; end if
         if (k <= nzg-1_ik) then ; c(k) = -kf(k) / dz_node(k) ; else ; c(k) = 0.0_wp ; end if
         b(k) = c_eff(k) * dz(k) / dt - a(k) - c(k)
         r(k) = c_eff(k) * dz(k) / dt * t_n(k) + dz(k) * q_src(k)
      end do
      r(1)   = r(1)   + g_top                                            ! top Neumann surface flux
      r(nzg) = r(nzg) + geo                                              ! bottom geothermal source
      call thomas_solve(a, b, c, r, t_new, nzg)
   end subroutine soil_heat_be_step

end module meds_column_energy
