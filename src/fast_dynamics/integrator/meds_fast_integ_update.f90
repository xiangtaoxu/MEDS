! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_integ_update -- the end of a fast step, the same for both schemes.             !
! update_dynamic_state commits the march's end state, the water network's soil water and   !
! pond, the tissues' energy, the condensate and the snow pack: the one place a step writes !
! the dynamic state. close_ledgers closes the seven conservation ledgers on the march's    !
! boundary-flux amounts (docs/science/numerical_scheme.md §5b). Also the commit's parts    !
! and the store totals the ledgers read.                                                   !
!==========================================================================================!
module meds_fast_integ_update
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num, mmh2o, cp_liq
   use meds_config,           only : meds_config_t
   use meds_site_state_vars,  only : site_surface_t
   use meds_fast_integ_types,       only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,     &
                                     step_report_t, march_state_t, frozen_t, boundary_flux_t,              &
                                     tissue_coefficients_t, diagnose_soil_temps
   use meds_diag_vars,        only : CO_TRANSP
   use meds_therm_lib,        only : internal_energy_liquid, cas_enthalpy_to_temp, tissue_energy_to_temp
   use meds_budget_check,     only : budget_check, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR,           &
                                     BUDGET_CO2_RATE_FLOOR
   implicit none
   private

   public :: commit_amounts_t, update_dynamic_state, close_ledgers
   public :: commit_tissue_energy, unpack_column_state, transp_per_leaf_area   ! used only by the tests

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
   ! mass and enthalpy paired, into soil layer 1 (docs/science/numerical_scheme.md §7, item 4) !
   ! before the soil's temperatures are diagnosed.                                            !
   !---------------------------------------------------------------------------------------!
   subroutine update_dynamic_state(dt_fast, cfg, site, static, frozen, acc, dt_warm, y_out, dynamic, out, moved)
      real(wp),                   intent(in)    :: dt_fast
      type(meds_config_t),        intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(frozen_t),      intent(in)    :: frozen
      type(boundary_flux_t),       intent(in)    :: acc       !< the march's boundary-flux amounts
      real(wp),                   intent(in)    :: dt_warm   !< [s] the step the next march starts from
      type(march_state_t),       intent(inout) :: y_out     !< the march's end state, committed as stored
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
      !----- The condensate lands before the soil's temperatures are diagnosed, so the top layer's  !
      !      temperature includes its water and enthalpy. ------------------------------------------!
      if (acc%whole_cond > 0.0_wp) then
         moved%cond_mass = acc%whole_cond
         moved%cond_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (boundary_flux_t)
         call deposit_condensate(y_out, site%soil_profile%dz(1), moved%cond_mass, moved%cond_enth)
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
      type(frozen_t),     intent(in)    :: frozen
      type(march_state_t),      intent(in)    :: y, y_out
      type(boundary_flux_t),      intent(in)    :: acc
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

   !=======================================================================================!
   !  The commit's parts and the store totals the ledgers read.                              !
   !=======================================================================================!

   !----- A cohort's transpiration over a step [kg/m2 ground] as the rate per leaf area the            !
   !      diagnostics report [mol H2O/m2 leaf/s]; zero for a leafless cohort. -------------------------!
   elemental pure function transp_per_leaf_area(amount, dt, nplant, leaf_area) result(rate)
      real(wp), intent(in) :: amount      !< [kg/m2 ground] over dt
      real(wp), intent(in) :: dt          !< [s]
      real(wp), intent(in) :: nplant      !< [plant/m2]
      real(wp), intent(in) :: leaf_area   !< [m2 leaf/plant]
      real(wp) :: rate
      rate = 0.0_wp
      if (nplant * leaf_area > tiny_num) rate = amount / dt / (nplant * leaf_area) / mmh2o
   end function transp_per_leaf_area

   !---------------------------------------------------------------------------------------!
   ! Whole-column store totals for the conservation ledgers [per m2 ground]. Plain ascending loops, !
   ! shared by both schemes, so a ledger built from these is the same arithmetic on both.           !
   !---------------------------------------------------------------------------------------!
   pure function soil_water_store(theta, dz, nsl) result(w)     ! [kg/m2]
      real(wp),    intent(in) :: theta(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: w
      integer(ik) :: k
      w = 0.0_wp
      do k = 1_ik, nsl
         w = w + theta(k) * dz(k) * rho_h2o
      end do
   end function soil_water_store

   pure function soil_energy_store(soil_energy, dz, nsl) result(e)   ! [J/m2] from [J/m3]
      real(wp),    intent(in) :: soil_energy(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: e
      integer(ik) :: k
      e = 0.0_wp
      do k = 1_ik, nsl
         e = e + soil_energy(k) * dz(k)
      end do
   end function soil_energy_store

   pure function plant_water_store(nplant, leaf_mass, wood_mass, n) result(w)   ! [kg/m2] from [kg/plant]
      real(wp),    intent(in) :: nplant(:), leaf_mass(:), wood_mass(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(nplant(1:n) * (leaf_mass(1:n) + wood_mass(1:n)))
   end function plant_water_store

   pure function canopy_film_store(leaf_surf, wood_surf, n) result(w)   ! [kg/m2 ground], already ground-referenced
      real(wp),    intent(in) :: leaf_surf(:), wood_surf(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(leaf_surf(1:n) + wood_surf(1:n))
   end function canopy_film_store

   !---------------------------------------------------------------------------------------!
   ! Canopy-film capacity clamp on the committed state: cap each cohort's leaf/wood film at        !
   ! leaf_max*LAI / wood_max*WAI and floor it at 0; the clipped excess (overflow) and the floored    !
   ! shortfall (deficit) are returned so the caller books them in the whole-column ledger.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine clamp_canopy_film(y_out, lai, wai, leaf_max, wood_max, n, surf_overflow, surf_deficit)
      type(march_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: lai(:), wai(:)
      real(wp),             intent(in)    :: leaf_max(:), wood_max(:)   !< [kg/m2 leaf], [kg/m2 wood] per cohort
      integer(ik),          intent(in)    :: n
      real(wp),             intent(inout) :: surf_overflow, surf_deficit   !< [kg/m2] accumulated
      real(wp)    :: leaf_cap_i, wood_cap_i
      integer(ik) :: i
      do i = 1_ik, n
         leaf_cap_i = leaf_max(i) * lai(i)
         wood_cap_i = wood_max(i) * wai(i)
         surf_overflow = surf_overflow + max(0.0_wp, y_out%leaf_surf_water(i) - leaf_cap_i)          &
                                        + max(0.0_wp, y_out%wood_surf_water(i) - wood_cap_i)
         surf_deficit  = surf_deficit  + max(0.0_wp, -y_out%leaf_surf_water(i))                      &
                                        + max(0.0_wp, -y_out%wood_surf_water(i))
         y_out%leaf_surf_water(i) = min(max(y_out%leaf_surf_water(i), 0.0_wp), leaf_cap_i)
         y_out%wood_surf_water(i) = min(max(y_out%wood_surf_water(i), 0.0_wp), wood_cap_i)
      end do
   end subroutine clamp_canopy_film

   !---------------------------------------------------------------------------------------!
   ! Deposit condensed canopy-air vapour into soil layer 1 as paired mass [kg/m2] + liquid         !
   ! enthalpy [J/m2] (an internal CAS -> soil transfer; the ledger needs no boundary term).         !
   !---------------------------------------------------------------------------------------!
   pure subroutine deposit_condensate(y_out, dz1, mass, enth)
      type(march_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: dz1, mass, enth
      y_out%theta(1)       = y_out%theta(1) + mass / (rho_h2o * dz1)
      y_out%soil_energy(1) = y_out%soil_energy(1) + enth / dz1
   end subroutine deposit_condensate

   !---------------------------------------------------------------------------------------!
   ! Unpack the committed column state into the dynamic state: CAS twins (and their                !
   ! temperature), soil energy and moisture, plant internal water, canopy films; not the pond.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine unpack_column_state(y_out, n, nsl, dynamic)
      type(march_state_t),       intent(in)    :: y_out
      integer(ik),                intent(in)    :: n, nsl
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      dynamic%cas%can_enthalpy = y_out%cas_enthalpy
      dynamic%cas%can_shv = y_out%cas_shv
      dynamic%cas%can_co2 = y_out%cas_co2
      dynamic%cas%can_temp = cas_enthalpy_to_temp(y_out%cas_enthalpy, y_out%cas_shv)
      dynamic%soil_energy%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      dynamic%soil_water%theta(1:nsl)       = y_out%theta(1:nsl)
      dynamic%leaf_water_mass(1:n) = y_out%leaf_water_mass(1:n)
      dynamic%wood_water_mass(1:n) = y_out%wood_water_mass(1:n)
      dynamic%leaf_surf_water(1:n) = y_out%leaf_surf_water(1:n)
      dynamic%wood_surf_water(1:n) = y_out%wood_surf_water(1:n)
      !----- The pond (w_surface, w_surface_energy) is left out on purpose: each scheme commits it from  !
      !      the water network, so writing it here would commit it twice (test_state_combinators). ---!
   end subroutine unpack_column_state

   !------------------------------------------------------------------------------------------------!
   ! commit_tissue_energy -- each tissue's end-of-step energy and temperature, and the tissue store !
   ! the whole-column ledger counts (docs/science/vegetation_energy_dynamics.md eqs. (9)-(10)). The  !
   ! energy gains the heat cap*(t_avg - t0), t_avg being the march's time-average of where the      !
   ! tissue relaxed to, and the enthalpy its water gained (docs/science/soil_plant_water.md eq. (5)); !
   ! the temperature is read back with the end-of-step water. A masked-off process adds nothing; a  !
   ! leaf with no heat capacity takes its wood's temperature.                                        !
   !------------------------------------------------------------------------------------------------!
   pure subroutine commit_tissue_energy(tissue, nplant, dt, t_leaf_avg, t_wood_avg, n,                          &
                                        dynamic, store0, store1)
      type(tissue_coefficients_t), intent(in)    :: tissue
      real(wp),                    intent(in)    :: nplant(:)                     !< [plant/m2]
      real(wp),                    intent(in)    :: dt                            !< [s] dt_fast
      real(wp),                    intent(in)    :: t_leaf_avg(:), t_wood_avg(:)  !< [K]
      integer(ik),                 intent(in)    :: n
      type(fast_dynamic_state_t),  intent(inout) :: dynamic
      real(wp),                    intent(out)   :: store0, store1                !< [J/m2] the tissue store
      real(wp)    :: heat_leaf, heat_wood, enth_leaf, enth_wood, w_leaf, w_wood
      integer(ik) :: i
      store0 = 0.0_wp
      store1 = 0.0_wp
      do i = 1_ik, n
         heat_leaf = tissue%leaf_hcap_per_dt(i) * dt * (t_leaf_avg(i) - tissue%t_leaf0(i))
         heat_wood = tissue%wood_hcap_per_dt(i) * dt * (t_wood_avg(i) - tissue%t_wood0(i))
         enth_leaf = tissue%leaf_water_enth(i) * dt
         enth_wood = tissue%wood_water_enth(i) * dt
         store0 = store0 + dynamic%leaf_energy(i) + dynamic%wood_energy(i)
         dynamic%leaf_energy(i) = dynamic%leaf_energy(i) + heat_leaf + enth_leaf
         dynamic%wood_energy(i) = dynamic%wood_energy(i) + heat_wood + enth_wood
         store1 = store1 + dynamic%leaf_energy(i) + dynamic%wood_energy(i)
         w_leaf = max(dynamic%leaf_water_mass(i), 0.0_wp) * nplant(i)
         w_wood = max(dynamic%wood_water_mass(i), 0.0_wp) * nplant(i)
         if (tissue%wood_dry_hcap(i) + w_wood * cp_liq > 0.0_wp)                                        &
            dynamic%wood_temp(i) = tissue_energy_to_temp(dynamic%wood_energy(i), w_wood, tissue%wood_dry_hcap(i))
         dynamic%leaf_temp(i) = dynamic%wood_temp(i)
         if (tissue%leaf_dry_hcap(i) + w_leaf * cp_liq > 0.0_wp)                                        &
            dynamic%leaf_temp(i) = tissue_energy_to_temp(dynamic%leaf_energy(i), w_leaf, tissue%leaf_dry_hcap(i))
      end do
   end subroutine commit_tissue_energy

end module meds_fast_integ_update
