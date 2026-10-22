! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_update -- the end of a fast step, the same for both schemes.                    !
! update_dynamic_state commits the march's end state, the water network's soil water and    !
! pond, the tissues' energy, the condensate and the snow pack: the one place a step writes  !
! the dynamic state. close_ledgers closes the seven conservation ledgers on the march's      !
! boundary-flux amounts (docs/science/numerical_scheme.md §5b).                              !
!==========================================================================================!
module meds_fast_update
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o
   use meds_config,           only : meds_config_t
   use meds_site_state_vars,  only : site_surface_t
   use meds_fast_types,       only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,     &
                                     step_report_t, column_state_t, column_frozen_t, column_bflux_t
   use meds_column_state_ops, only : clamp_canopy_film, unpack_column_state, diagnose_soil_temps,                 &
                                     commit_tissue_energy, transp_per_leaf_area, deposit_condensate,              &
                                     soil_water_store, soil_energy_store, plant_water_store, canopy_film_store
   use meds_diag_vars,        only : CO_TRANSP
   use meds_therm_lib,        only : internal_energy_liquid
   use meds_budget_check,     only : budget_check, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR,           &
                                     BUDGET_CO2_RATE_FLOOR
   implicit none
   private

   public :: commit_amounts_t, update_dynamic_state, close_ledgers

   !----- What the update moved besides the march, which the ledgers book: the films' overflow at   !
   !      their capacity and deficit below zero, the condensate deposited in the top soil layer with  !
   !      its enthalpy, and the tissue store before and after. -----------------------------------!
   type :: commit_amounts_t
      real(wp) :: surf_overflow = 0.0_wp, surf_deficit = 0.0_wp   !< [kg/m2]
      real(wp) :: cond_mass     = 0.0_wp                          !< [kg/m2]
      real(wp) :: cond_enth     = 0.0_wp                          !< [J/m2]
      real(wp) :: tissue_store0 = 0.0_wp, tissue_store1 = 0.0_wp  !< [J/m2]
   end type commit_amounts_t

contains

   !---------------------------------------------------------------------------------------!
   ! update_dynamic_state -- commit one fast step. Soil water and the pond are the water      !
   ! network's: the march moved theta at the network's steady rate, and its end-of-step theta !
   ! and pond are committed once, here. A film evaporating at its frozen conductance can      !
   ! overdraw below zero, which would create water in the next step's interception, so the    !
   ! films are floored at zero and capped at capacity, and the update reports both. Each      !
   ! tissue's energy is committed from its time-average temperature over the march, the      !
   ! complement of what the canopy air received (docs/science/vegetation_energy_dynamics.md   !
   ! §One formulation). The condensate, dew and fog, lands inside the column, so it goes,      !
   ! mass and enthalpy paired, into soil layer 1 (docs/science/numerical_scheme.md §7, item 4).!
   !---------------------------------------------------------------------------------------!
   subroutine update_dynamic_state(dt_fast, cfg, site, static, frozen, acc, dt_warm, y_out, dynamic, out, moved)
      real(wp),                   intent(in)    :: dt_fast
      type(meds_config_t),        intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(column_frozen_t),      intent(in)    :: frozen
      type(column_bflux_t),       intent(in)    :: acc       !< the march's boundary-flux amounts
      real(wp),                   intent(in)    :: dt_warm   !< [s] the step the next march starts from
      type(column_state_t),       intent(inout) :: y_out     !< the march's end state, committed as stored
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(commit_amounts_t),     intent(out)   :: moved
      real(wp)    :: t_leaf_avg(static%n), t_wood_avg(static%n)
      integer(ik) :: n, nsl

      n   = static%n
      nsl = site%soil_profile%n_active
      y_out%theta(1:nsl)     = frozen%hydrology%theta1(1:nsl)
      y_out%w_surface        = frozen%hydrology%w_surface1
      y_out%w_surface_energy = frozen%hydrology%w_surface_energy1
      if (cfg%canopy_water_on) then
         call clamp_canopy_film(y_out, static%lai, static%wai, cfg%pft%leaf_surf_water_max(static%pft(1:n)),    &
                                cfg%pft%wood_surf_water_max(static%pft(1:n)), n, moved%surf_overflow,          &
                                moved%surf_deficit)
      end if

      call unpack_column_state(y_out, n, nsl, dynamic)
      dynamic%soil_water%w_surface        = y_out%w_surface
      dynamic%soil_water%w_surface_energy = y_out%w_surface_energy
      call diagnose_soil_temps(y_out, site%soil_profile%dry_heat_capacity, nsl, dynamic%soil_energy%soil_temp,    &
                               dynamic%soil_energy%soil_fliq)
      t_leaf_avg(1:n) = acc%tissue_leaf_int(1:n) / dt_fast
      t_wood_avg(1:n) = acc%tissue_wood_int(1:n) / dt_fast
      call commit_tissue_energy(frozen%tissue, static%nplant, dt_fast, t_leaf_avg, t_wood_avg, n, dynamic,        &
                                moved%tissue_store0, moved%tissue_store1)
      !----- Each cohort's transpiration over the step, as the leaf water was debited. -------------!
      if (out%capture .and. allocated(acc%transp_coh))                                                       &
         out%cdiag(CO_TRANSP, 1:n) = transp_per_leaf_area(acc%transp_coh(1:n), dt_fast, static%nplant(1:n),     &
                                                          static%leaf_area(1:n))
      if (acc%whole_cond > 0.0_wp) then
         moved%cond_mass = acc%whole_cond
         moved%cond_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (column_bflux_t)
         call deposit_condensate(y_out, site%soil_profile%dz(1), moved%cond_mass, moved%cond_enth)
         dynamic%soil_water%theta(1)        = y_out%theta(1)
         dynamic%soil_energy%soil_energy(1) = y_out%soil_energy(1)
      end if
      !----- The pack the snow stage advanced, and the warm start of the next march. --------------!
      dynamic%snow          = frozen%snow%pack
      dynamic%adapt_dt_last = dt_warm
   end subroutine update_dynamic_state

   !---------------------------------------------------------------------------------------!
   ! close_ledgers -- the seven conservation ledgers, from the march's b-weighted boundary-flux !
   ! amounts, which close to round-off (docs/science/numerical_scheme.md §5b): the canopy air's  !
   ! energy, water and CO2, the soil's energy and water, and the whole column's water and      !
   ! energy; and the step's canopy-air -> atmosphere exports. y is the march's start state and  !
   ! y_out the committed end. With halt ([fast].debug_error) a budget that does not close stops !
   ! the run; the tolerance scales with the boundary flux, not the store.                       !
   !---------------------------------------------------------------------------------------!
   subroutine close_ledgers(dt_fast, site, forcing, static, frozen, y, y_out, acc, moved, halt, report)
      real(wp),                  intent(in)    :: dt_fast
      type(site_surface_t),      intent(in)    :: site
      type(fast_forcing_t),      intent(in)    :: forcing
      type(fast_static_state_t), intent(in)    :: static
      type(column_frozen_t),     intent(in)    :: frozen
      type(column_state_t),      intent(in)    :: y, y_out
      type(column_bflux_t),      intent(in)    :: acc
      type(commit_amounts_t),    intent(in)    :: moved
      logical,                   intent(in)    :: halt
      type(step_report_t),       intent(inout) :: report
      real(wp)    :: cas_mass_capacity, cas_molar_capacity, enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, w_plant0, w_plant1, w_surface0, e_pond0, e_pond1
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1, intercept_total
      integer(ik) :: n, nsl

      n   = static%n
      nsl = site%soil_profile%n_active
      cas_mass_capacity  = frozen%cas%box%air_mass_capacity
      cas_molar_capacity = frozen%cas%box%air_molar_capacity
      enth0 = y%cas_enthalpy
      shv0  = y%cas_shv
      co20  = y%cas_co2
      enth1 = y_out%cas_enthalpy
      shv1  = y_out%cas_shv
      co21  = y_out%cas_co2
      !----- The pond at the start, as the march started from it, and at the end, the network's. -----!
      w_surface0 = y%w_surface
      e_pond0    = y%w_surface_energy
      e_pond1    = frozen%hydrology%w_surface_energy1
      e_soil0 = soil_energy_store(y%soil_energy,     site%soil_profile%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, site%soil_profile%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     site%soil_profile%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, site%soil_profile%dz, nsl)
      !----- Plant water is a store: it absorbs the step's mismatch between transpiration and uptake  !
      !      (docs/science/plant_hydraulics.md eqs. (1)-(2)). ----------------------------------------!
      w_plant0 = plant_water_store(static%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(static%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- The films, per m2 of ground, valued at the liquid enthalpy the intercepted water arrived  !
      !      with, as surface_derivs charges the tissue (0 with canopy_water_on off). ----------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))

      report%exports%atm_heat_export = acc%atm_heat_out
      report%exports%atm_vap_export  = acc%atm_vap_out
      call budget_check(report%ledgers%cas_energy,                                                           &
                        cas_mass_capacity*enth0, cas_mass_capacity*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'cas_energy', halt)
      call budget_check(report%ledgers%cas_water,                                                            &
                        cas_mass_capacity*shv0,  cas_mass_capacity*shv1,  acc%cas_vap_in,  acc%cas_vap_out,  &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'cas_water', halt)
      call budget_check(report%ledgers%cas_co2,                                                              &
                        cas_molar_capacity*co20, cas_molar_capacity*co21, acc%cas_co2_in,  acc%cas_co2_out,  &
                        dt_fast, BUDGET_CO2_RATE_FLOOR, 'cas_co2', halt)
      !----- The condensate deposit is an inflow to the soil budgets, though the whole-column ledger   !
      !      sees it as a transfer between two of its stores. -----------------------------------------!
      call budget_check(report%ledgers%soil_energy, e_soil0, e_soil1, acc%soil_enth_in + moved%cond_enth,          &
                        acc%soil_enth_out, dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'soil_energy', halt)
      !----- Soil water: theta^n -> theta1, in through the top face net of evaporation and            !
      !      exfiltration (q_top), out through the bottom and to the roots, plus the condensate. ------!
      call budget_check(report%ledgers%soil_water,  w_soil0, w_soil1,                                              &
                        frozen%hydrology%q_top*rho_h2o*dt_fast + moved%cond_mass,                          &
                        (frozen%hydrology%drainage + frozen%roots%uptake)*dt_fast,                         &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'soil_water', halt)
      !----- Whole water: rainfall, snowfall and shed water in; the CAS vapour export, drainage and    !
      !      runoff out; soil, pond, plant water, films, snow pack and canopy-air vapour stored.       !
      !      Sublimation leaves through the vapour export and meltwater moves pack -> pond, so        !
      !      neither needs a term. -------------------------------------------------------------------!
      call budget_check(report%ledgers%whole_water,                                                                &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%snow%swe0,  &
                        w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 + surf_water1 &
                        + frozen%snow%swe1,                                                                &
                        acc%whole_wat_in + (forcing%rainfall + forcing%snowfall + static%shed_water_rate)*dt_fast, &
                        acc%whole_wat_out + (frozen%hydrology%runoff_surf + frozen%hydrology%drainage)*dt_fast  &
                        + moved%surf_overflow - moved%surf_deficit,                                        &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'whole_water', halt)
      !----- Whole energy: meltwater passes pack -> pond, both stores, so it telescopes. -------------!
      call budget_check(report%ledgers%whole_energy,                                                               &
                        e_soil0 + cas_mass_capacity*enth0 + surf_enth0 + frozen%snow%enth0 + e_pond0               &
                        + moved%tissue_store0,                                                             &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond1               &
                        + moved%tissue_store1,                                                             &
                        acc%whole_enth_in + intercept_total*dt_fast*internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                          + frozen%snow%acc_enth + frozen%hydrology%ground_inflow_enth*dt_fast, &
                        acc%whole_enth_out + (moved%surf_overflow - moved%surf_deficit)                     &
                        * internal_energy_liquid(frozen%hydrology%t_film_valuation)                         &
                                           + frozen%hydrology%runoff_enth*dt_fast,                          &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'whole_energy', halt)
   end subroutine close_ledgers

end module meds_fast_update
