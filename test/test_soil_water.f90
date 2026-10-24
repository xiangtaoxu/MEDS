! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_soil_water -- the soil's water flows (meds_soil_water;                              !
! docs/science/soil_biophysics.md, docs/science/soil_plant_water.md §3).                   !
!   1. DARCY (soil_water_flux, soil_biophysics.md eq. (3)): its four derivatives against   !
!      central differences; the pond face at zero head is the infiltration capacity.       !
!   2. EVAPORATION (ground_evaporation_rate, ground_evaporation_setup): its derivatives     !
!      against central differences, and the closed form where no dry surface layer forms. !
!==========================================================================================!
program test_soil_water
   use meds_test_assert,   only : check_true, check_close, test_report
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : grav, grav_head, rho_h2o, r_wv, p_std
   use meds_soil_hydr_lib, only : SOIL_RETENTION_VG
   use meds_soil_water,    only : ground_evaporation_rate, ground_evaporation_setup, soil_water_flux, ground_evap_env_t
   use meds_therm_lib,     only : sat_specific_humidity
   use meds_soil_params,   only : soil_profile_t, soil_params_t
   implicit none

   call test_soil_water_flux()
   call test_evaporation()

   call test_report('test_soil_water')

contains

   !----- Relative agreement of an analytic derivative with a central difference. -----------!
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   !----- 1. Darcy flow across the face between two layers 0.4 m apart (the geometric-mean face),  !
   !      and the pond's face at zero head. ----------------------------------------------------!
   subroutine test_soil_water_flux()
      real(wp) :: q, qa, qb, dqa, dqb, dka, dkb, dum, h, ka, kb, pa, pb, psi1, k1, ks, z1, q_ref
      logical  :: ok
      ka = 3.0e-7_wp
      kb = 8.0e-9_wp
      pa = -0.02_wp
      pb = -0.6_wp
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, q, dqa, dqb, dka, dkb)
      ok = .true.
      h = 1.0e-7_wp
      call soil_water_flux(ka, kb, 0.4_wp, pa + h, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb, 0.4_wp, pa - h, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dqa, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb + h, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb, 0.4_wp, pa, pb - h, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dqb, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      h = 1.0e-6_wp * ka
      call soil_water_flux(ka + h, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka - h, kb, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dka, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      h = 1.0e-6_wp * kb
      call soil_water_flux(ka, kb + h, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qa, dum, dum, dum, dum)
      call soil_water_flux(ka, kb - h, 0.4_wp, pa, pb, 0.07_wp, 0.07_wp, qb, dum, dum, dum, dum)
      if (.not. fd_agrees(dkb, (qa - qb) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call check_true('darcy: derivatives match finite differences', ok)
      !----- The pond face at zero head: the infiltration capacity,                               !
      !      sqrt(K1 Ksat) (1 + (-psi1)/(-z1)), times the water density. -------------------------!
      psi1 = -3.0_wp
      k1 = 2.0e-9_wp
      ks = 4.5e-6_wp
      z1 = -0.015_wp
      call soil_water_flux(ks, k1, 0.5_wp, 0.0_wp, grav_head * psi1, -z1, -z1, q, dum, dum, dum, dum)
      q_ref = sqrt(k1 * ks) * (1.0_wp + (-psi1) / (-z1)) * rho_h2o
      call check_close('darcy: pond face at zero head = infiltration capacity', q, q_ref, 1.0e-12_wp * q_ref)
   end subroutine test_soil_water_flux

   !----- 2. Ground evaporation from BCI's clay top layer, over water contents from 0.11 to 0.43, at  !
   !      the potential its retention curve gives. ---------------------------------------------!
   subroutine test_evaporation()
      type(soil_profile_t)   :: sp
      type(ground_evap_env_t) :: f
      type(soil_params_t)    :: o
      real(wp) :: theta, head, psi, e, det, dep, ep, em, dum, ref, h
      integer  :: i
      logical  :: ok_d, ok_v
      sp%n_active = 3
      sp%retention = SOIL_RETENTION_VG
      sp%dz(1:3) = 0.03_wp
      sp%theta_sat(1:3) = 0.4314_wp
      sp%theta_res(1:3) = 0.0981_wp
      sp%curve_par_a(1:3) = 0.7607_wp
      sp%curve_par_n(1:3) = 1.3554_wp
      sp%ksat(1:3) = 4.54e-6_wp
      f%t_ground = 300.0_wp
      f%q_air = 0.012_wp
      f%rho_air = 1.15_wp
      f%r_aero = 60.0_wp
      call ground_evaporation_setup(sp, o, f)
      ok_d = .true.
      ok_v = .true.
      do i = 1, 30
         theta = 0.10_wp + 0.011_wp * real(i, wp)
         head  = -(1.0_wp / 0.7607_wp) * (((theta - 0.0981_wp) / (0.4314_wp - 0.0981_wp))              &
                  ** (-1.0_wp / (1.0_wp - 1.0_wp / 1.3554_wp)) - 1.0_wp) ** (1.0_wp / 1.3554_wp)
         psi   = head * grav_head                      ! [MPa], the potential the kernel takes
         call ground_evaporation_rate(theta, psi, f, o, e, det, dep)
         h = 1.0e-7_wp
         call ground_evaporation_rate(theta + h, psi, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta - h, psi, f, o, em, dum, dum)
         if (.not. fd_agrees(det, (ep - em) / (2*h), 1.0e-4_wp, 1.0e-14_wp)) ok_d = .false.
         h = 1.0e-5_wp * max(abs(psi), grav_head)
         call ground_evaporation_rate(theta, psi + h, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta, psi - h, f, o, em, dum, dum)
         if (.not. fd_agrees(dep, (ep - em) / (2*h), 1.0e-4_wp, 1.0e-14_wp)) ok_d = .false.
         !----- Moist enough that no dry surface layer forms: the closed form rho (alpha q_sat - q_air) !
         !      / r_aero, with Philip's alpha = exp(psi g / (R_v T)). --------------------------------!
         if (theta > o%dsl_theta_init * 0.4314_wp + 0.01_wp) then
            ref = f%rho_air * (exp(head * grav / (r_wv * f%t_ground)) * sat_specific_humidity(f%t_ground, p_std, 1.0_wp) &
                               - f%q_air) / f%r_aero
            if (abs(e - ref) > 1.0e-10_wp * abs(ref)) ok_v = .false.
         end if
      end do
      call check_true('evaporation: derivatives match finite differences', ok_d)
      call check_true('evaporation: the closed form where no dry surface layer forms', ok_v)
   end subroutine test_evaporation

end program test_soil_water
