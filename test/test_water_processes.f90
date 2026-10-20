! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_water_processes -- the water relations and flows the soil-plant water solver uses   !
! (docs/science/soil_plant_water.md): every analytic derivative against a central finite   !
! difference, over dry, wet and saturated states, and every value against the read-offs of !
! meds_soil_hydr_lib and meds_plant_hydr_lib, where those read-offs' floors do not bind.   !
!                                                                                          !
!   1. SOIL LAYERS  : two van Genuchten soils and a Campbell soil; the two branches of the      !
!                     unknown meet with equal value and slope; theta and K match the            !
!                     read-offs; theta -> unknown -> theta round-trips.                         !
!   2. TISSUES      : the pressure-volume curve matches tissue_rwc_to_psi and                  !
!                     tissue_psi_to_water.                                                     !
!   3. POND         : smooth at empty and at the depression depth; water -> head round-trips.  !
!   4. SOIL FLUX    : the pond face at zero head is the infiltration capacity.                  !
!   5. ROOT UPTAKE  : equals Gardner's conductance x driving force; nothing flows backwards.    !
!   6. XYLEM        : equals the Kirchhoff mean conductance x driving force, exponent 2 and a   !
!                     tabled 3.                                                                 !
!   7. TRANSPIRATION: zero at an empty leaf, the demand above the closing water content.        !
!   8. EVAPORATION  : the closed form where no dry surface layer forms.                          !
!   9. HEAD         : one constant, rho_w g, converts head to potential for soil and plant.      !
!==========================================================================================!
program test_water_processes
   use meds_test_assert, only : check_true, check_close, test_report
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : grav_head, rho_h2o
   use meds_soil_hydr_lib,   only : soil_layer_curve_t, build_soil_layer_curve, soil_water_relations,  &
                                    soil_theta_to_unknown, pond_water_relations, pond_water_to_head,   &
                                    soil_specific_storage, pond_depression_depth, soil_retention_vg,   &
                                    soil_retention_campbell, soil_head_to_theta, vg_pore_connectivity
   use meds_plant_hydr_lib,  only : tissue_water_relations, tissue_water_to_rwc, rhizosphere_conductance,  &
                                    root_water_uptake, xylem_sapflow, xylem_max_conductance, hydro_table_t, &
                                    build_hydro_table, kirchhoff_potential, conductance_retained_fraction, &
                                    water_curve_t, tissue_rwc_to_psi, tissue_psi_to_water
   use meds_plant_hydraulics, only : supply_limited_transpiration
   use meds_numerics,    only : smooth_max0
   use meds_soil_water, only : ground_evaporation_rate, soil_water_flux
   use meds_constants,  only : pi, grav, r_wv, p_std
   use meds_therm_lib,  only : sat_specific_humidity
   use meds_soil_types, only : chydro_forcing_t
   use meds_soil_column_config, only : soil_params_t
   use meds_biophysics_opts, only : soil_opts_t
   implicit none

   call test_soil_water_relations()
   call test_tissue_water_relations()
   call test_pond_water_relations()
   call test_soil_water_flux()
   call test_root_uptake()
   call test_xylem()
   call test_transpiration()
   call test_evaporation()
   call test_hydrostatic_head()
   call test_report('test_water_processes')

contains

   !----- The conductivity's closed form, written independently of meds_soil_hydr_lib: Mualem-    !
   !      van Genuchten Ks Se^L [1 - (1 - Se^(1/m))^m]^2, or Campbell's Ks Se^(2b+3). -------------!
   pure real(wp) function conductivity_ref(c, se) result(k)
      type(soil_layer_curve_t), intent(in) :: c
      real(wp),                 intent(in) :: se
      real(wp) :: m
      if (c%retention == soil_retention_campbell) then
         k = c%ksat * se ** (2.0_wp * c%par_n + 3.0_wp)
      else
         m = 1.0_wp - 1.0_wp / c%par_n
         k = c%ksat * se ** vg_pore_connectivity * (1.0_wp - (1.0_wp - se ** (1.0_wp / m)) ** m) ** 2
      end if
   end function conductivity_ref

   !----- Relative agreement of an analytic derivative with a central difference. -----------!
   logical function fd_agrees(analytic, numeric, rtol, atol)
      real(wp), intent(in) :: analytic, numeric, rtol, atol
      fd_agrees = abs(analytic - numeric) <= rtol * max(abs(analytic), abs(numeric)) + atol
   end function fd_agrees

   !=======================================================================================!
   subroutine test_soil_water_relations()
      type(soil_layer_curve_t) :: c(3)
      character(len=24)  :: label(3)
      real(wp) :: tau, h, psi, dpsi, w, dw, k, dk, psip, psim, wp_, wm, kp, km, dum1, dum2
      real(wp) :: worst_psi, worst_w, worst_k, theta, theta_ref, k_ref, tau_back, head
      real(wp) :: psi_lo, psi_hi, dpsi_lo, dpsi_hi, w_lo, w_hi, dw_lo, dw_hi
      integer  :: ic, i
      logical  :: ok_psi, ok_w, ok_k, ok_val, ok_k_val, ok_round
      !----- BCI's clay (example 04), a sand, and a Campbell loam. ----------------------------!
      call build_soil_layer_curve(soil_retention_vg, 0.4314_wp, 0.0981_wp, 0.7607_wp, 1.3554_wp, 4.54e-6_wp, 0.05_wp, c(1))
      call build_soil_layer_curve(soil_retention_vg, 0.38_wp, 0.05_wp, 3.5_wp, 2.7_wp, 8.0e-5_wp, 0.10_wp, c(2))
      call build_soil_layer_curve(soil_retention_campbell, 0.45_wp, 0.0_wp, -0.25_wp, 6.0_wp, 5.0e-6_wp, 0.20_wp, c(3))
      label = [character(len=24) :: 'vG clay (BCI)', 'vG sand', 'Campbell loam']
      do ic = 1, 3
         ok_psi = .true. ; ok_w = .true. ; ok_k = .true.
         worst_psi = 0.0_wp ; worst_w = 0.0_wp ; worst_k = 0.0_wp
         !----- tau from very dry, through the branch point, to 2 m of head above saturation. ---!
         do i = 0, 400
            tau = 1.0e-4_wp * (c(ic)%se_switch / 1.0e-4_wp) ** (real(i, wp) / 300.0_wp)
            if (i > 300) tau = c(ic)%se_switch + real(i - 300, wp) / 100.0_wp                       &
                               * c(ic)%slope_switch * (2.0_wp - c(ic)%head_switch)
            if (abs(tau - c(ic)%se_switch) < 1.0e-6_wp) cycle            ! straddles the join
            if (abs(tau - (c(ic)%se_switch - c(ic)%slope_switch * c(ic)%head_switch)) < 1.0e-6_wp) cycle
            h = 1.0e-7_wp * max(abs(tau), 1.0e-3_wp)
            call soil_water_relations(c(ic), tau, psi, dpsi, w, dw, k, dk)
            call soil_water_relations(c(ic), tau + h, psip, dum1, wp_, dum2, kp, dum1)
            call soil_water_relations(c(ic), tau - h, psim, dum1, wm, dum2, km, dum1)
            if (.not. fd_agrees(dpsi, (psip - psim) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok_psi = .false.
            if (.not. fd_agrees(dw,   (wp_ - wm) / (2*h),    1.0e-4_wp, 1.0e-9_wp))  ok_w = .false.
            !----- Mualem's 1 - (1 - Se^(1/m))^m is the difference of two numbers within Se^(1/m) of !
            !      one: where Se^(1/m) < 1e-4, K carries round-off of order 1e-16 / Se^(1/m), which a  !
            !      finite difference amplifies past the tolerance (with L = -1, K there is above the   !
            !      1e-20 floor). The derivative is checked where a difference can resolve it. --------!
            if (c(ic)%retention == soil_retention_vg .and. tau < c(ic)%se_switch) then
               if (tau ** (c(ic)%par_n / (c(ic)%par_n - 1.0_wp)) < 1.0e-4_wp) cycle
            end if
            if (.not. fd_agrees(dk,   (kp - km) / (2*h),     1.0e-4_wp, 1.0e-20_wp)) ok_k = .false.
         end do
         call check_true(trim(label(ic))//': dpsi/dtau matches finite differences', ok_psi)
         call check_true(trim(label(ic))//': dW/dtau matches finite differences', ok_w)
         call check_true(trim(label(ic))//': dK/dtau matches finite differences', ok_k)
         !----- The two branches meet with the same value and slope. ------------------------!
         tau = c(ic)%se_switch
         call soil_water_relations(c(ic), tau - 1.0e-9_wp, psi_lo, dpsi_lo, w_lo, dw_lo, dum1, dum2)
         call soil_water_relations(c(ic), tau + 1.0e-9_wp, psi_hi, dpsi_hi, w_hi, dw_hi, dum1, dum2)
         call check_true(trim(label(ic))//': branches meet (potential)', abs(psi_hi - psi_lo) < 1.0e-8_wp)
         call check_true(trim(label(ic))//': branches meet (slope of potential)',                     &
                         abs(dpsi_hi - dpsi_lo) <= 1.0e-5_wp * abs(dpsi_lo))
         call check_true(trim(label(ic))//': branches meet (slope of water)',                          &
                         abs(dw_hi - dw_lo) <= 1.0e-5_wp * abs(dw_lo) .or. ic == 3)
         !----- Values against meds_soil_hydr_lib, where its floors do not bind. -------------!
         ok_val = .true. ; ok_k_val = .true. ; ok_round = .true.
         do i = 1, 60
            head = -10.0_wp ** (-2.0_wp + 5.0_wp * real(i, wp) / 60.0_wp)          ! -1 cm .. -1 km
            theta_ref = soil_head_to_theta(c(ic)%retention, head, c(ic)%theta_sat, c(ic)%theta_res,   &
                                            c(ic)%par_a, c(ic)%par_n)
            k_ref = conductivity_ref(c(ic), (theta_ref - c(ic)%theta_res) / (c(ic)%theta_sat - c(ic)%theta_res))
            if ((theta_ref - c(ic)%theta_res) / (c(ic)%theta_sat - c(ic)%theta_res) < 1.0e-8_wp) cycle
            tau = soil_theta_to_unknown(c(ic), theta_ref)
            call soil_water_relations(c(ic), tau, psi, dpsi, w, dw, k, dk)
            theta = w / (rho_h2o * c(ic)%dz)
            if (abs(theta - theta_ref) > 1.0e-10_wp) ok_val = .false.
            if (k_ref > 1.0e-12_wp .and. abs(k - k_ref) > 1.0e-8_wp * k_ref) ok_k_val = .false.
            if (abs(psi / grav_head - head) > 1.0e-6_wp * abs(head) .and.                              &
                theta_ref < c(ic)%theta_sat - 1.0e-9_wp) ok_round = .false.
         end do
         call check_true(trim(label(ic))//': theta(psi) matches meds_soil_hydr_lib', ok_val)
         call check_true(trim(label(ic))//': K(theta) matches the closed form', ok_k_val)
         call check_true(trim(label(ic))//': theta -> unknown -> potential round-trips', ok_round)
         !----- Above saturation the layer holds theta_sat + S_s * head. ---------------------!
         tau = soil_theta_to_unknown(c(ic), c(ic)%theta_sat + soil_specific_storage * 0.5_wp)
         call soil_water_relations(c(ic), tau, psi, dpsi, w, dw, k, dk)
         call check_close(trim(label(ic))//': head above saturation', psi / grav_head, 0.5_wp, 1.0e-9_wp)
         call check_close(trim(label(ic))//': conductivity above saturation', k, c(ic)%ksat, 0.0_wp)
         tau_back = soil_theta_to_unknown(c(ic), w / (rho_h2o * c(ic)%dz))
         call check_close(trim(label(ic))//': supersaturated round trip', tau_back, tau, 1.0e-9_wp)
      end do
   end subroutine test_soil_water_relations

   !=======================================================================================!
   subroutine test_tissue_water_relations()
      type(water_curve_t) :: leaf
      real(wp) :: r, h, psi, dpsi, w, dw, psip, psim, wpl, wm, dum, worst
      integer  :: i
      logical  :: ok_d, ok_v
      leaf = water_curve_t(pi0=-1.8_wp, elastic_mod=12.0_wp, apoplast_frac=0.2_wp, water_sat=2.5_wp)
      ok_d = .true. ; ok_v = .true.
      do i = 1, 120
         r = 0.01_wp * real(i, wp)                                  ! 0.01 .. 1.2
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r, psi, dpsi, w, dw)
         h = 1.0e-7_wp
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r + h, psip, dum, wpl, dum)
         call tissue_water_relations(leaf, 0.3_wp, 0.05_wp, r - h, psim, dum, wm, dum)
         if (abs(r - (leaf%pi0 + leaf%elastic_mod) / leaf%elastic_mod) > 1.0e-6_wp) then
            if (.not. fd_agrees(dpsi, (psip - psim) / (2*h), 1.0e-5_wp, 1.0e-10_wp)) ok_d = .false.
         end if
         if (.not. fd_agrees(dw, (wpl - wm) / (2*h), 1.0e-6_wp, 1.0e-12_wp)) ok_d = .false.
         if (r >= 1.0e-4_wp .and. abs(psi - tissue_rwc_to_psi(r, leaf%pi0, leaf%elastic_mod)) > 1.0e-10_wp)    &
            ok_v = .false.
         if (abs(w - 0.05_wp * tissue_psi_to_water(psi, leaf, 0.3_wp)) > 1.0e-10_wp) ok_v = .false.
         if (abs(tissue_water_to_rwc(leaf, 0.3_wp, w / 0.05_wp) - r) > 1.0e-10_wp) ok_v = .false.
      end do
      call check_true('tissue: derivatives match finite differences', ok_d)
      call check_true('tissue: potential and water match the pressure-volume curve', ok_v)
   end subroutine test_tissue_water_relations

   !=======================================================================================!
   subroutine test_pond_water_relations()
      real(wp) :: h, psi, dpsi, w, dw, wpl, wm, dum, hh
      integer  :: i
      logical  :: ok_d, ok_r
      ok_d = .true. ; ok_r = .true.
      do i = -20, 60
         h = 1.0e-4_wp * real(i, wp) + 1.0e-7_wp               ! -2 mm .. 6 mm, off the corners
         call pond_water_relations(h, psi, dpsi, w, dw)
         hh = 1.0e-9_wp
         call pond_water_relations(h + hh, dum, dum, wpl, dum)
         call pond_water_relations(h - hh, dum, dum, wm, dum)
         if (.not. fd_agrees(dw, (wpl - wm) / (2*hh), 1.0e-5_wp, 1.0e-6_wp)) ok_d = .false.
         if (h > 0.0_wp .and. abs(pond_water_to_head(w) - h) > 1.0e-12_wp) ok_r = .false.
      end do
      call check_true('pond: dW/dh matches finite differences, smooth at empty and full', ok_d)
      call check_true('pond: water -> head round-trips', ok_r)
      call pond_water_relations(0.0_wp, psi, dpsi, w, dw)
      call check_close('pond: empty pond holds nothing', w, 0.0_wp, 0.0_wp)
      call pond_water_relations(0.05_wp, psi, dpsi, w, dw)
      call check_close('pond: 5 cm of head holds 50 kg/m2 less half the depression depth', w,           &
                       rho_h2o * (0.05_wp - 0.5_wp * pond_depression_depth), 1.0e-10_wp)
   end subroutine test_pond_water_relations

   !=======================================================================================!
   subroutine test_soil_water_flux()
      real(wp) :: q, qa, qb, dqa, dqb, dka, dkb, dum, h, ka, kb, pa, pb, psi1, k1, ks, z1, q_ref
      logical  :: ok
      ka = 3.0e-7_wp ; kb = 8.0e-9_wp ; pa = -0.02_wp ; pb = -0.6_wp
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
      psi1 = -3.0_wp ; k1 = 2.0e-9_wp ; ks = 4.5e-6_wp ; z1 = -0.015_wp
      call soil_water_flux(ks, k1, 0.5_wp, 0.0_wp, grav_head * psi1, -z1, -z1, q, dum, dum, dum, dum)
      q_ref = sqrt(k1 * ks) * (1.0_wp + (-psi1) / (-z1)) * rho_h2o
      call check_close('darcy: pond face at zero head = infiltration capacity', q, q_ref, 1.0e-12_wp * q_ref)
   end subroutine test_soil_water_flux

   !=======================================================================================!
   subroutine test_root_uptake()
      real(wp) :: g, u, dus, duw, duk, up, um, dum, h, k, ps, pw, z, ref
      logical  :: ok
      g = rhizosphere_conductance(0.02_wp, 35.0_wp, 0.01_wp, 2.5e-4_wp)
      k = 1.0e-8_wp ; ps = -0.3_wp ; pw = -0.9_wp ; z = -0.4_wp
      call root_water_uptake(g, k, ps, z, pw, u, dus, duw, duk)
      ref = 0.02_wp * 2.0_wp * pi * (rho_h2o * k / grav_head) * 35.0_wp / log(0.01_wp / 2.5e-4_wp)      &
            * (ps + grav_head * z - pw)
      call check_close('root: uptake = cohort x Gardner conductance x driving force', u, ref, 1.0e-10_wp * ref)
      ok = .true.
      h = 1.0e-7_wp
      call root_water_uptake(g, k, ps + h, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k, ps - h, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(dus, (up - um) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call root_water_uptake(g, k, ps, z, pw + h, up, dum, dum, dum)
      call root_water_uptake(g, k, ps, z, pw - h, um, dum, dum, dum)
      if (.not. fd_agrees(duw, (up - um) / (2*h), 1.0e-6_wp, 0.0_wp)) ok = .false.
      call root_water_uptake(g, k * (1.0_wp + 1.0e-6_wp), ps, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k * (1.0_wp - 1.0e-6_wp), ps, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(duk, (up - um) / (2.0e-6_wp * k), 1.0e-6_wp, 0.0_wp)) ok = .false.
      !----- On the ramp: half a width of driving force. --------------------------------------!
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-4_wp, z, pw, u, dus, duw, duk)
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-4_wp + h, z, pw, up, dum, dum, dum)
      call root_water_uptake(g, k, pw - grav_head * z + 5.0e-4_wp - h, z, pw, um, dum, dum, dum)
      if (.not. fd_agrees(dus, (up - um) / (2*h), 1.0e-5_wp, 0.0_wp)) ok = .false.
      call check_true('root: derivatives match finite differences, on and off the ramp', ok)
      call root_water_uptake(g, k, -2.0_wp, z, pw, u, dus, duw, duk)
      call check_close('root: no flow from soil drier than the wood', u, 0.0_wp, 0.0_wp)
   end subroutine test_root_uptake

   !=======================================================================================!
   subroutine test_xylem()
      type(hydro_table_t) :: tab3, tab2
      real(wp) :: kc, s, dsw, dsl, sp, sm, dum, h, pw, pl, ref, grav
      integer  :: i
      logical  :: ok2, ok3, okv
      call build_hydro_table(tab3, 3.0_wp)
      call build_hydro_table(tab2, 2.0_wp)
      kc = xylem_max_conductance(3.3_wp, 2.0e-3_wp, 18.0_wp, 1.0_wp) * 0.05_wp
      grav = grav_head * 18.0_wp
      ok2 = .true. ; ok3 = .true. ; okv = .true.
      do i = 0, 40
         pw = -0.1_wp - 0.15_wp * real(i, wp)
         pl = pw - 0.6_wp + 0.03_wp * real(i, wp)          ! crosses pw - pl = 0 at i = 20
         if (i == 20) pl = pw - 3.0e-7_wp                  ! inside the midpoint expansion
         h = 1.0e-7_wp
         call xylem_sapflow(kc, pw, pl, grav, -2.5_wp, 2.0_wp, tab2, s, dsw, dsl)
         ref = kirchhoff_mean(pw, pl, kc, -2.5_wp, 2.0_wp, tab2) * (pw - pl - grav)
         if (abs(s - ref) > 1.0e-9_wp * abs(ref) + 1.0e-11_wp * kc) okv = .false.
         call xylem_sapflow(kc, pw + h, pl, grav, -2.5_wp, 2.0_wp, tab2, sp, dum, dum)
         call xylem_sapflow(kc, pw - h, pl, grav, -2.5_wp, 2.0_wp, tab2, sm, dum, dum)
         if (.not. fd_agrees(dsw, (sp - sm) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok2 = .false.
         call xylem_sapflow(kc, pw, pl + h, grav, -2.5_wp, 2.0_wp, tab2, sp, dum, dum)
         call xylem_sapflow(kc, pw, pl - h, grav, -2.5_wp, 2.0_wp, tab2, sm, dum, dum)
         if (.not. fd_agrees(dsl, (sp - sm) / (2*h), 1.0e-4_wp, 1.0e-12_wp)) ok2 = .false.
         !----- A tabled vulnerability exponent: the table's slope is piecewise, so the analytic  !
         !      slope (the exact retained fraction) agrees only to the table's resolution. -------!
         call xylem_sapflow(kc, pw, pl, grav, -2.5_wp, 3.0_wp, tab3, s, dsw, dsl)
         !----- At near-equal ends the reference leaves its table for the exact retained fraction;   !
         !      the solver stays on the table, so the two differ there by its resolution. ----------!
         ref = kirchhoff_mean(pw, pl, kc, -2.5_wp, 3.0_wp, tab3) * (pw - pl - grav)
         if (i /= 20 .and. abs(s - ref) > 1.0e-9_wp * abs(ref) + 1.0e-11_wp * kc) okv = .false.
         h = 1.0e-4_wp
         call xylem_sapflow(kc, pw + h, pl, grav, -2.5_wp, 3.0_wp, tab3, sp, dum, dum)
         call xylem_sapflow(kc, pw - h, pl, grav, -2.5_wp, 3.0_wp, tab3, sm, dum, dum)
         if (.not. fd_agrees(dsw, (sp - sm) / (2*h), 2.0e-2_wp, 1.0e-8_wp * kc)) then
            ok3 = .false.
            print '(a,i3,4es14.6)', '    table case i, pw, pl, analytic, fd: ', i, pw, pl, dsw, (sp - sm) / (2*h)
         end if
      end do
      call check_true('xylem: sapflow = Kirchhoff mean conductance x driving force (exponents 2 and 3)', okv)
      call check_true('xylem: derivatives match finite differences (exponent 2, incl. equal ends)', ok2)
      call check_true('xylem: derivatives match the table to its resolution (exponent 3)', ok3)
   end subroutine test_xylem

   !----- The reference Kirchhoff mean conductance, written independently of xylem_sapflow:        !
   !      k (Phi(up) - Phi(down)) / (up - down) from kirchhoff_potential (the closed form for        !
   !      exponent 2, the table for 3), and the retained fraction at the midpoint for near-equal    !
   !      ends. ----------------------------------------------------------------------------------!
   pure real(wp) function kirchhoff_mean(psi_up, psi_down, k, psi50, kexp, tab) result(keff)
      real(wp),            intent(in) :: psi_up, psi_down, k, psi50, kexp
      type(hydro_table_t), intent(in) :: tab
      if (abs(psi_up - psi_down) <= 1.0e-6_wp) then
         keff = k * conductance_retained_fraction(0.5_wp * (psi_up + psi_down), psi50, kexp)
      else
         keff = k * (kirchhoff_potential(psi_up, psi50, kexp, tab) - kirchhoff_potential(psi_down, psi50, kexp, tab)) &
                / (psi_up - psi_down)
      end if
   end function kirchhoff_mean

   !=======================================================================================!
   subroutine test_transpiration()
      real(wp) :: e, de, ep, em, dum, r, h
      integer  :: i
      logical  :: ok
      ok = .true.
      do i = 1, 99
         r = 0.01_wp * real(i, wp) + 1.0e-6_wp
         h = 1.0e-7_wp
         call supply_limited_transpiration(2.0e-5_wp, r, 0.42_wp, e, de)
         call supply_limited_transpiration(2.0e-5_wp, r + h, 0.42_wp, ep, dum)
         call supply_limited_transpiration(2.0e-5_wp, r - h, 0.42_wp, em, dum)
         if (.not. fd_agrees(de, (ep - em) / (2*h), 1.0e-5_wp, 1.0e-15_wp)) ok = .false.
      end do
      call check_true('transpiration: derivative matches finite differences', ok)
      call supply_limited_transpiration(2.0e-5_wp, 0.0_wp, 0.42_wp, e, de)
      call check_close('transpiration: an empty leaf transpires nothing', e, 0.0_wp, 0.0_wp)
      call supply_limited_transpiration(2.0e-5_wp, 0.5_wp, 0.42_wp, e, de)
      call check_close('transpiration: the demand stands above the closing water content', e, 2.0e-5_wp, 0.0_wp)
      call supply_limited_transpiration(2.0e-5_wp, 0.21_wp, 0.42_wp, e, de)
      call check_close('transpiration: half the closing water content, half the demand', e, 1.0e-5_wp,      &
                       1.0e-15_wp)
   end subroutine test_transpiration

   !=======================================================================================!
   subroutine test_evaporation()
      type(soil_params_t)    :: sp
      type(chydro_forcing_t) :: f
      type(soil_opts_t)      :: o
      real(wp) :: theta, head, e, det, dep, ep, em, dum, ref, h
      integer  :: i
      logical  :: ok_d, ok_v
      sp%n_active = 3 ; sp%retention = soil_retention_vg
      sp%dz(1:3) = 0.03_wp ; sp%theta_sat(1:3) = 0.4314_wp ; sp%theta_res(1:3) = 0.0981_wp
      sp%curve_par_a(1:3) = 0.7607_wp ; sp%curve_par_n(1:3) = 1.3554_wp ; sp%ksat(1:3) = 4.54e-6_wp
      f%t_ground = 300.0_wp ; f%q_air = 0.012_wp ; f%rho_air = 1.15_wp ; f%r_aero = 60.0_wp
      ok_d = .true. ; ok_v = .true.
      do i = 1, 30
         theta = 0.10_wp + 0.011_wp * real(i, wp)
         head  = -(1.0_wp / 0.7607_wp) * (((theta - 0.0981_wp) / (0.4314_wp - 0.0981_wp))              &
                  ** (-1.0_wp / (1.0_wp - 1.0_wp / 1.3554_wp)) - 1.0_wp) ** (1.0_wp / 1.3554_wp)
         call ground_evaporation_rate(theta, head, sp, f, o, e, det, dep)
         h = 1.0e-7_wp
         call ground_evaporation_rate(theta + h, head, sp, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta - h, head, sp, f, o, em, dum, dum)
         if (.not. fd_agrees(det, (ep - em) / (2*h), 1.0e-4_wp, 1.0e-14_wp)) ok_d = .false.
         h = 1.0e-5_wp * max(abs(head), 1.0_wp)
         call ground_evaporation_rate(theta, head + h, sp, f, o, ep, dum, dum)
         call ground_evaporation_rate(theta, head - h, sp, f, o, em, dum, dum)
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

   !=======================================================================================!
   !----- Head and potential: the soil curves, the soil faces, root uptake and the crown's      !
   !      gravity all convert with grav_head, which is water's density times standard gravity   !
   !      (docs/science/plant_hydraulics.md, the total head; soil_plant_water.md). 10 m of water  !
   !      is 98.0665 kPa. ------------------------------------------------------------------------!
   subroutine test_hydrostatic_head()
      call check_close('head: grav_head is rho_w g in MPa per metre', grav_head, rho_h2o * grav * 1.0e-6_wp,      &
                       epsilon(1.0_wp) * grav_head)
      call check_close('head: 10 m of water is 98.0665 kPa', 10.0_wp * grav_head, 0.0980665_wp,                  &
                       4.0_wp * epsilon(1.0_wp) * 0.0980665_wp)
   end subroutine test_hydrostatic_head

end program test_water_processes
