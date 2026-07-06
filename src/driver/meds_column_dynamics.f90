!==========================================================================================!
! meds_column_dynamics -- the fast-timescale (dt_fast) integrator that couples the fast          !
! biophysics processes into one sub-daily operator-split sweep (design Part II). The canopy-air-  !
! space (CAS) twins are the shared coupling reservoir: each dt_fast the aerodynamics kernel sets   !
! the conductances, the leaf and GROUND surface energy balances emit sensible + latent fluxes into  !
! the CAS, the soil thermal column is advanced from the net ground heat flux, and the three CAS      !
! twins (enthalpy / specific humidity / CO2) are advanced IMPLICITLY in the atmospheric-exchange     !
! term, gated by the friction velocity.                                                             !
!                                                                                          !
! Coupled so far: aerodynamics -> {per-cohort diagnostic leaf balance, ground surface balance} ->    !
! CAS enthalpy/vapour/CO2 (implicit) -> soil thermal column (implicit BE-Thomas). Carries the         !
! design's §3.5 fix: the atm<->CAS scalar conductance is the profile-factored gah = rho*ustar*temp1   !
! (heat), gaw = rho*ustar*temp2 (vapour), gac = rho_dmol*ustar*temp2 (CO2) -- not the bare rho*ustar.  !
!                                                                                          !
! The leaf temperature is DIAGNOSED from a linearized steady-state balance (a near-massless leaf      !
! equilibrates within a dt_fast). Prescribed for now (next layers): absorbed radiation, stomatal      !
! conductance, net biotic CO2, and the soil-moisture profile (from the soil-water column). Ground      !
! evaporation uses the energy-kernel conductance form; when the hydrology column lands, the §3.6        !
! single-authority reconciliation picks the DSL form. Soil hydrology, photosynthesis/hydraulics, the    !
! stepper hook, and cross-demography state persistence remain to wire.                                  !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     leaf_energy_env_t, soil_params_t, soil_thermal_params_t,  &
                                     energy_forcing_t, energy_opts_t, energy_flux_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_column_energy,    only : ground_surface_balance, soil_energy_flux
   use meds_thermo,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     d_sat_vapor_pressure_dt, enthalpy_vapor
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
   implicit none
   private

   public :: column_config_t, column_forcing_t, column_fast_step

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero          !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal  !< leaf/wood thermal params
      type(soil_params_t)         :: soil         !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal !< soil thermal texture
      type(energy_opts_t)         :: energy       !< soil-thermal solver options
   end type column_config_t

   !----- Prescribed per-step forcing the higher layers (RT, photosynthesis, met, hydrology). --!
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: nee_biotic    = 0.0_wp   !< [umol/m2/s] net biotic CO2 source (Reco - GPP)
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW / net LW per cohort (leaf)
      real(wp), allocatable :: gsw(:), fs_open(:)       !< [m/s], [-] stomatal conductance, open fraction
      real(wp), allocatable :: soil_water(:)            !< [m3/m3] soil moisture profile (from the water column)
   end type column_forcing_t

contains

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch: aerodynamics -> per-cohort    !
   !  DIAGNOSTIC leaf balance + GROUND surface balance -> CAS three-twin implicit update ->      !
   !  soil thermal column. Budgets accumulate the closed residuals. Cohort arrays BOTTOM(1)->TOP.!
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, ccfg, aenv, ageom, n, lai, wai, height, crown,         &
                               leaf_width, branch_diam, forc, bio, aero,                       &
                               be_energy, be_water, be_co2, be_soil)
      real(wp),                intent(in)    :: dt_fast
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv        !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      integer(ik),             intent(in)    :: n
      real(wp),                intent(in)    :: lai(n), wai(n), height(n), crown(n)
      real(wp),                intent(in)    :: leaf_width(n), branch_diam(n)
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero        !< preallocated (alloc_aero_out)
      type(budget_t),          intent(inout) :: be_energy, be_water, be_co2, be_soil

      type(leaf_energy_env_t) :: env_g
      type(energy_forcing_t)  :: eforc
      type(energy_flux_t)     :: sflux
      real(wp)    :: tcas, qcas, press, rho, h_coeff, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: g_tr, le_ref, dtl, tl, h_i, le_i, transp_i
      real(wp)    :: coh_h, coh_le, coh_transp, t_ground, g_top, h_ground, le_ground, ground_w
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol, src_enth, src_vap
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      integer(ik) :: i, nsl

      !----- 1. Refresh the aerodynamics env from the current CAS state, then solve. ---------!
      bio%cas%can_temp = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      tcas = bio%cas%can_temp ; qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = bio%soil_e%soil_temp(1)
      call canopy_aerodynamics(ccfg%aero, aenv, ageom, n, height, lai, crown, bio%leaf_temp,   &
                               bio%leaf_temp, leaf_width, branch_diam, aero)

      !----- 2. Per-cohort DIAGNOSTIC leaf steady-state -> sensible + latent flux into CAS. ---!
      qsat_c = sat_specific_humidity(tcas, press)
      dqdt   = 0.622_wp * press / max((press - 0.378_wp * sat_e(tcas))**2, tiny_num)           &
               * d_sat_vapor_pressure_dt(tcas)
      coh_h = 0.0_wp ; coh_le = 0.0_wp ; coh_transp = 0.0_wp
      do i = 1_ik, n
         h_coeff = ccfg%veg_thermal%effarea_heat * lai(i) * aero%leaf_gbh(i) * rho * cp_air
         g_tr    = 0.0_wp
         if (aero%leaf_gbw(i) + forc%gsw(i) > tiny_num) then
            g_tr = ccfg%veg_thermal%effarea_transp * lai(i) * forc%fs_open(i)                  &
                   * aero%leaf_gbw(i) * forc%gsw(i) / (aero%leaf_gbw(i) + forc%gsw(i))
         end if
         lw_slope = 4.0_wp * ccfg%veg_thermal%leaf_emiss * stefan * tcas**3 * lai(i)
         le_slope = latent_heat_vap * rho * g_tr * dqdt
         le_ref   = latent_heat_vap * rho * g_tr * (qsat_c - qcas)
         dtl = (forc%abs_sw(i) + forc%abs_lw(i) - le_ref) / max(h_coeff + le_slope + lw_slope, tiny_num)
         tl  = tcas + dtl
         bio%leaf_temp(i) = tl
         h_i      = h_coeff * dtl
         le_i     = le_ref + le_slope * dtl
         transp_i = le_i / latent_heat_vap
         coh_h      = coh_h      + h_i
         coh_le     = coh_le     + le_i
         coh_transp = coh_transp + transp_i
      end do

      !----- 3. GROUND surface balance (t_ground from the soil column; ggnet conductance). ----!
      t_ground = bio%soil_e%soil_temp(1)
      env_g%gbh = aero%ggnet ; env_g%gbw = aero%ggnet
      env_g%can_temp = tcas ; env_g%can_shv = qcas ; env_g%rho_air = rho ; env_g%press = press
      env_g%abs_sw = forc%abs_sw_ground ; env_g%abs_lw = forc%abs_lw_ground
      call ground_surface_balance(t_ground, env_g, g_top, h_ground, le_ground)
      ground_w = le_ground / max(enthalpy_vapor(t_ground), tiny_num)      ! ground evap mass [kg/m2/s]

      !----- 4. CAS three-twin update: IMPLICIT in the profile-factored atm exchange (§3.5). --!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      wcap     = rho      * bio%cas%can_depth
      ccap     = can_dmol * bio%cas%can_depth
      gah      = rho      * aero%ustar * aero%temp1
      gaw      = rho      * aero%ustar * aero%temp2
      gac      = can_dmol * aero%ustar * aero%temp2
      src_enth = coh_h + coh_le + h_ground + le_ground                    ! [W/m2]  sensible + latent
      src_vap  = coh_transp + ground_w                                    ! [kg/m2/s] vapour mass

      enth0 = bio%cas%can_enthalpy ; shv0 = qcas ; co20 = bio%cas%can_co2
      enth1 = (wcap*enth0 + dt_fast*(src_enth + gah*forc%enthalpy_atm)) / (wcap + dt_fast*gah)
      shv1  = (wcap*shv0  + dt_fast*(src_vap  + gaw*forc%shv_atm))       / (wcap + dt_fast*gaw)
      co21  = (ccap*co20  + dt_fast*(forc%nee_biotic + gac*forc%co2_atm)) / (ccap + dt_fast*gac)

      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- 5. Soil thermal column: advance from the net ground heat flux (top Neumann BC). --!
      nsl = ccfg%soil%n_active
      eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
      eforc%soil_water(1:nsl)     = forc%soil_water(1:nsl)
      eforc%w_flux(1:nsl)         = 0.0_wp
      eforc%root_heat_sink(1:nsl) = 0.0_wp
      call soil_energy_flux(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)

      !----- 6. Closed-budget accumulation. --------------------------------------------------!
      call budget_accumulate(be_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm,  &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(be_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,          &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(be_co2,    ccap*co20, ccap*co21, forc%nee_biotic + gac*forc%co2_atm,  &
                             gac*co21, dt_fast, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      !----- Soil column closes its OWN energy budget (kernel residual); track it directly. --!
      be_soil%n_check = be_soil%n_check + 1_ik
      be_soil%resid   = sflux%energy_resid
      be_soil%worst   = max(be_soil%worst, abs(sflux%energy_resid))
      if (.not. closure_ok(sflux%energy_resid, abs(g_top) + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp))       &
         be_soil%n_fail = be_soil%n_fail + 1_ik
   end subroutine column_fast_step

   !----- Saturation vapour pressure [Pa] (Bolton; local mirror of the thermo form for dqdt). -!
   pure real(wp) function sat_e(t_k) result(esat)
      real(wp), intent(in) :: t_k
      real(wp) :: tc
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
   end function sat_e

end module meds_column_dynamics
