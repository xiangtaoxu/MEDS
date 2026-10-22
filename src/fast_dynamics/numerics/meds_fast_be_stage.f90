! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_be_stage -- the implicit stage of the fast-loop march, and the water advance that  !
! rides with it. One ESDIRK stage (docs/science/numerical_scheme.md §3 `ark`, eq. (1)) is a     !
! backward-Euler solve of the canopy air against the leaf, wood and ground temperatures, by a   !
! small Newton on a forward-difference Jacobian. column_be_stage is the stage,                  !
! newton_surface_solve and jac_surface its inner solve, cas_box_commit the box commit;          !
! advance_plant_water and advance_surf_water_full move plant and film water once per step.      !
! meds_fast_ark keeps the tableau and the march.                                                !
!==========================================================================================!
module meds_fast_be_stage
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air, rho_h2o
   use meds_soil_types, only : energy_forcing_t, energy_flux_t
   use meds_patch_state_vars,   only : soil_energy_column_t
   use meds_fast_time_derivs, only : surface_derivs, cas_conductances
   use meds_cas_biophysics, only : cas_column_step_implicit, cas_column_t, cas_source_t
   use meds_column_state_ops, only : state_init, assemble_soil_energy_forcing
   use meds_fast_types, only : column_config_t, column_state_t, column_frozen_t, surface_state_t,                 &
                               cas_boundary_t, surface_tend_t, stage_bflux_t
   use meds_soil_energy, only : soil_energy_step_implicit
   use meds_therm_lib, only : cas_enthalpy_to_temp, internal_energy_liquid, internal_energy_to_temp
   use meds_numerics,     only : sufficient_decrease
   implicit none
   private

   public :: column_be_stage, advance_plant_water, advance_surf_water_full

contains

   !---------------------------------------------------------------------------------------!
   ! column_be_stage -- one backward-Euler stage Y = base + dt*f(Y) of the stiff block (each ESDIRK2   !
   ! stage, dt = gamma*h, is one call): the CAS twins, implicit in the atmospheric exchange (niter > 1:  !
   ! the coupled 2x2 Newton), and soil heat (implicit, Thomas). Soil water moves at the water network's !
   ! steady rate; plant water is passed through and advanced once (advance_plant_water).               !
   !---------------------------------------------------------------------------------------!
   subroutine column_be_stage(y, frozen, col_config, n, nsl, dt, y_out, surf_tend, niter, bf)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      !----- The stage's state and its own surface tendencies, the ones its CAS source and bf used; the   !
      !      caller b-weights transp_c into the plant water, so the leaves' water debit and the canopy     !
      !      air's vapour credit agree (docs/science/numerical_scheme.md §5a). Both filled in place. -----!
      type(column_state_t),  intent(inout) :: y_out
      type(surface_tend_t),  intent(inout) :: surf_tend
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux rates for the ARK ledger

      type(surface_state_t)      :: y_stage
      type(cas_boundary_t)       :: cas_stage   !< the CAS boundary with this stage's conductances
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp)    :: enth1, shv1, e_infil, e_drain, t_cas1, water_heat(nsl)
      real(wp)    :: co21
      integer(ik) :: np, nfeval
      logical     :: ok

      np = 1_ik
      if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, col_config%soil_thermal%soil_dry_heat_capacity(1), t_ground, fliq1)

      cas_mass_capacity = frozen%cas%cas_mass_capacity
      cas_molar_capacity = frozen%cas%cas_molar_capacity
      cas_stage = frozen%cas              ! plain scalars: a cheap copy, overridden with this stage's conductances below

      !----- Re-solve the surface layer at this stage's canopy-air state, so the stage pays for the       !
      !      ventilation its own temperature earns (docs/science/numerical_scheme.md §2); the box commit,  !
      !      the Newton and the ledger all read these values. cas_stage carries them with mo_live cleared,  !
      !      so surface_derivs does not re-solve on each Newton residual evaluation. ---------------------!
      call cas_conductances(frozen%cas, y%cas_enthalpy, y%cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      cas_stage%g_atm_heat = g_atm_heat
      cas_stage%g_atm_vapour = g_atm_vapour
      cas_stage%g_atm_co2 = g_atm_co2
      cas_stage%mo_live = .false.

      !----- np == 1: one evaluation at the stage start; np > 1: the coupled 2x2 Newton. The final      !
      !      surf_tend drives the box commit and the soil sinks, one flux per interface. ------------!
      if (np <= 1_ik) then
         y_stage%cas_enthalpy = y%cas_enthalpy
         y_stage%cas_shv = y%cas_shv
         y_stage%cas_co2 = y%cas_co2
         call surface_derivs(y_stage, cas_stage, frozen%tissue, frozen%film, frozen%ground, frozen%snow,    &
                             t_ground, n, surf_tend)
      else
         call newton_surface_solve(y, cas_stage, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                   surf_tend, nfeval, ok)
      end if
      !----- One box commit for all three twins, in flux form, so the budgets close for any surf_tend: !
      !      CO2 rides the same box as enthalpy and humidity. -----------------------------------------!
      call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, cas_stage, cas_mass_capacity, cas_molar_capacity, &
                          g_atm_heat, g_atm_vapour, g_atm_co2, dt, enth1, shv1, co21)
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = co21

      !----- Soil heat: implicit backward Euler, Thomas solve. -------------------------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !----- The enthalpy of the water the network moved: across the faces at its time-mean flows, in     !
      !      from the pond at t_infil, out with root uptake and exfiltration at each layer's start-of-step  !
      !      temperature (what the wood and the pond were credited), and out with drainage                !
      !      (docs/science/soil_plant_water.md §5). The leaf pays its transpiration's full vapour         !
      !      enthalpy in surface_derivs. ---------------------------------------------------------------!
      water_heat(1:nsl) = frozen%roots%root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + frozen%hydrology%exfil_enth
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      e_drain = frozen%hydrology%drainage     * internal_energy_liquid(frozen%hydrology%t_bot)
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, y%theta,          &
                                        water_heat, frozen%hydrology%w_flux_frozen,                                &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain)
      call soil_energy_step_implicit(se, eforc, col_config%soil_thermal, col_config%soil, col_config%energy, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- Soil water is the network's, its end-of-step theta committed once in column_fast_step_ark     !
      !      (docs/science/numerical_scheme.md §5b). The stages move theta at its steady rate, so each      !
      !      layer holds the water whose enthalpy the faces brought in; a steady rate is integrated exactly !
      !      and adds nothing to the error estimate. ---------------------------------------------------!
      y_out%theta(1:nsl) = y%theta(1:nsl) + dt * frozen%hydrology%theta_rate(1:nsl)

      !----- The pond is passed through: column_fast_step_ark commits the network's end-of-step pond,  !
      !      and the pond is outside the error norm. --------------------------------------------------!
      y_out%w_surface      = y%w_surface
      y_out%w_surface_energy = y%w_surface_energy

      !----- Plant water mass is passed through (advance_plant_water advances it). -----------!
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !----- This stage's boundary-flux rates for the ledger, all at the committed state, so the amounts   !
      !      ark2_column_step and adaptive_ark_march b-weight and accumulate telescope to closure. -------!
      if (present(bf)) then
         t_cas1 = cas_enthalpy_to_temp(enth1, shv1)         ! committed CAS temperature, for the reported H
            bf%cas_enth_in  = surf_tend%src_enth + g_atm_heat*frozen%cas%enthalpy_atm
            bf%cas_enth_out = g_atm_heat*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + g_atm_vapour*frozen%cas%shv_atm
            bf%cas_vap_out  = g_atm_vapour*shv1
            bf%cas_co2_in   = frozen%cas%nee_biotic + g_atm_co2*frozen%cas%co2_atm
            bf%cas_co2_out  = g_atm_co2*y_out%cas_co2
            !----- The bottom face the soil solve committed (zero on the adiabatic base, the conduction  !
            !      to deep_temp under the Dirichlet anchor) crosses both the soil's and the column's book. !
            bf%soil_enth_in = surf_tend%g_top - eflux%bottom_heat + e_infil
            bf%soil_enth_out= sum(water_heat(1:nsl)) + e_drain
            !----- Soil water's storage change and flows are the network's, booked once per step in     !
            !      column_fast_step_ark: the per-stage soil-water amounts are zero, and the whole-column   !
            !      water here is only the CAS-vapour exchange. ------------------------------------------!
            bf%soil_wat_in  = 0.0_wp
            bf%soil_wat_out = 0.0_wp
            !----- ground_rad is the snowfac-blended radiation (abs_sw + abs_lw when bare). Infiltration,  !
            !      exfiltration and uptake move enthalpy between tracked stores and cancel; rain and runoff  !
            !      are added once at the outer level (docs/science/numerical_scheme.md §5b). -------------!
            bf%whole_enth_in= surf_tend%coh_rnet + frozen%snow%ground_rad - eflux%bottom_heat
            !----- The condensate is not a boundary loss: the caller deposits it in soil layer 1 with the  !
            !      enthalpy surface_derivs debited from the CAS (whole_cond, whole_cond_enth). -----------!
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%enthalpy_atm) + e_drain
            bf%whole_wat_in = 0.0_wp
            bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%shv_atm)
            bf%whole_cond   = surf_tend%cond                     ! deposited into soil layer 1, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! exactly what surface_derivs debited from the CAS
            !----- The reported H, against theta_atm (T at zref), the ledger's basis. ----------------------!
            bf%atm_heat_out = g_atm_heat*cp_air*(t_cas1 - frozen%cas%mo_theta_atm)
            bf%atm_vap_out  = g_atm_vapour*(shv1  - frozen%cas%shv_atm)
      end if
   end subroutine column_be_stage

   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- a direct 2x2 Newton solve of the coupled backward-Euler canopy-air block   !
   ! {R_H, R_q} = 0 in (cas_enthalpy, cas_shv), whose surface sources depend on both through tcas,    !
   ! qcas and qsat; robust near saturation, where a leaf <-> CAS Picard iteration is not. The forward- !
   ! difference Jacobian (jac_surface) captures the strong VPD self-limiting d src_vap/dq. A singular- !
   ! Jacobian fallback, a line search and an evaluation cap; it never stops the run. Returns surf_tend !
   ! at the solution, from which column_be_stage commits the box.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, cas, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, &
                                   surf_tend, nfeval, ok)
      type(column_state_t),   intent(in)    :: y
      type(cas_boundary_t),   intent(in)    :: cas      !< the stage's CAS boundary (its own conductances)
      type(column_frozen_t),  intent(in)    :: frozen   !< tissue / film / ground / snow for surface_derivs
      real(wp),               intent(in)    :: t_ground
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, cas_mass_capacity, g_atm_heat, g_atm_vapour
      type(surface_tend_t),   intent(inout) :: surf_tend   !< the evaluation at the solution, on return
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: y_stage
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy
      q0 = y%cas_shv
      Hk = H0
      qk = q0
      nfeval = 0_ik
      ok = .false.
      y_stage%cas_co2 = y%cas_co2
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, surf_tend)
      nfeval = nfeval + 1_ik
      R_H = cas_mass_capacity*(Hk - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%enthalpy_atm - Hk)
      R_q = cas_mass_capacity*(qk - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%shv_atm  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/cas_mass_capacity <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/cas_mass_capacity <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true.
            exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, &
                          J11, J12, J21, J22, nfeval)
         detJ = J11*J22 - J12*J21
         if (detJ <= DETEPS*abs(J11*J22) .or. detJ <= 0.0_wp) then       ! singular / sign-flipped guard
            delH = -R_H / max(J11, tiny_num)                             ! damped-diagonal (Picard-like) fallback
            delq = -R_q / max(J22, tiny_num)
         else
            delH = (-R_H*J22 + R_q*J12)/detJ
            delq = (-R_q*J11 + R_H*J21)/detJ
         end if
         lam = 1.0_wp
         rn0 = R_H*R_H + R_q*R_q
         Ht = Hk
         qt = qk
         RHt = R_H
         Rqt = R_q
         do ls = 1_ik, LS_MAX
            Ht = Hk + lam*delH
            qt = qk + lam*delq
            !----- No supersaturation clamp on the state: a jump there inflates the embedded error near  !
            !      RH = 1. Transient supersaturation is relaxed by surface_derivs' condensation sink. ----!
            y_stage%cas_enthalpy = Ht
            y_stage%cas_shv = qt
            call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,         &
               t_ground, n, surf_tend)
            nfeval = nfeval + 1_ik
            RHt = cas_mass_capacity*(Ht - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%enthalpy_atm - Ht)
            Rqt = cas_mass_capacity*(qt - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%shv_atm  - qt)
            if (sufficient_decrease(RHt*RHt + Rqt*Rqt, rn0, lam)) exit
            lam = 0.5_wp*lam
         end do
         Hk = Ht
         qk = qt
         R_H = RHt
         R_q = Rqt
         if (nfeval >= FEVAL_CAP) exit
      end do
      !----- surf_tend holds the last evaluation, which is at (Hk, qk): the start when the first test  !
      !      passes, else the accepted line-search trial. The caller commits the box from it. --------!
   end subroutine newton_surface_solve

   !---------------------------------------------------------------------------------------!
   ! The backward-Euler canopy-air box commit, through the shared box kernel                           !
   ! meds_cas_biophysics%cas_column_step_implicit (docs/science/canopy_air_space_biophysics.md         !
   ! eq. (2)): one implementation of the box, not a copy of its formulas here.                         !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_box_commit(h0, q0, c0, surf_tend, cas, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
                                  g_atm_co2, dt, h1, q1, c1)
      real(wp),               intent(in)  :: h0, q0, c0, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
           g_atm_co2, dt
      type(surface_tend_t),   intent(in)  :: surf_tend
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(out) :: h1, q1, c1
      type(cas_column_t) :: box
      type(cas_source_t) :: src
      box%air_mass_capacity        = cas_mass_capacity
      box%air_molar_capacity     = cas_molar_capacity
      box%atm_conductance_enthalpy = g_atm_heat
      box%atm_conductance_vapor  = g_atm_vapour
      box%atm_conductance_co2 = g_atm_co2
      box%atm_enthalpy             = cas%enthalpy_atm
      box%atm_specific_humidity = cas%shv_atm
      box%atm_co2 = cas%co2_atm
      src%surface_enthalpy_source  = surf_tend%src_enth
      src%surface_vapor_source = surf_tend%src_vap
      src%biotic_co2_source        = cas%nee_biotic
      call cas_column_step_implicit(h0, q0, c0, src, box, dt, h1, q1, c1)
   end subroutine cas_box_commit

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, J11, &
                          J12, J21, J22, &
                          nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt
      type(cas_boundary_t),   intent(in)    :: cas
      type(column_frozen_t),  intent(in)    :: frozen
      real(wp),               intent(in)    :: t_ground
      type(surface_tend_t),   intent(in)    :: surf_tend         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      y_stage%cas_co2 = co2
      y_stage%cas_enthalpy = Hk + dH
      y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, sfp)
      nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - surf_tend%src_enth)/dH
      dsv_dH = (sfp%src_vap - surf_tend%src_vap)/dH
      y_stage%cas_enthalpy = Hk
      y_stage%cas_shv = qk + dq
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, sfp)
      nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - surf_tend%src_enth)/dq
      dsv_dq = (sfp%src_vap - surf_tend%src_vap)/dq
      J11 = cas_mass_capacity/dt + g_atm_heat - dse_dH
      J12 =              - dse_dq
      J21 =              - dsv_dH
      J22 = cas_mass_capacity/dt + g_atm_vapour - dsv_dq
   end subroutine jac_surface

   !---------------------------------------------------------------------------------------!
   ! advance_plant_water -- each cohort's leaf and wood water over an ARK step at the network's      !
   ! time-mean rates: sapflow and uptake in; sapflow and the b-weighted transpiration the stages     !
   ! charged the canopy air out. The stages hold transpiration at the network's, so a dt_fast step  !
   ! ends where the network did. Nothing is floored: the network keeps every store inside its curve. !
   !---------------------------------------------------------------------------------------!
   pure subroutine advance_plant_water(y, frozen, n, dt, transp_c_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  !< [kg/m2 ground/s] per cohort, b-weighted over the step
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      do i = 1_ik, n
         y_out%leaf_water_mass(i) = y%leaf_water_mass(i)                                                  &
              + dt * (frozen%plant%sapflow_frozen(i) - transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num))
         y_out%wood_water_mass(i) = y%wood_water_mass(i)                                                  &
              + dt * (frozen%plant%uptake_frozen(i) - frozen%plant%sapflow_frozen(i))
      end do
   end subroutine advance_plant_water

   !---------------------------------------------------------------------------------------!
   ! advance_surf_water_full -- the canopy film over the full dt, the film analogue of                 !
   ! advance_plant_water: frozen interception in, b-weighted film evaporation out (the number the CAS  !
   ! was credited). Unclamped: column_fast_step_ark clamps it and books the overflow and deficit.     !
   !---------------------------------------------------------------------------------------!
   pure subroutine advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !----- Not floored: the caller's deficit bookkeeping must see a negative film to credit the       !
      !      shortfall to the whole-column outflow; a floor here would create that water. ------------!
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%film%intercept_leaf(i) - film_evap_leaf_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%film%intercept_wood(i) - film_evap_wood_bw(i))
      end do
   end subroutine advance_surf_water_full

end module meds_fast_be_stage
