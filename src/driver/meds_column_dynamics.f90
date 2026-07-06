!==========================================================================================!
! meds_column_dynamics -- the fast-timescale (dt_fast) integrator that couples the fast          !
! biophysics processes into one sub-daily operator-split sweep (design Part II). The canopy-air-  !
! space (CAS) twins are the shared coupling reservoir: each dt_fast the aerodynamics kernel sets   !
! the conductances, the leaf/wood surface energy balance emits sensible + latent fluxes into the   !
! CAS, and the three CAS twins (enthalpy / specific humidity / CO2) are advanced IMPLICITLY in     !
! the atmospheric-exchange term, gated by the friction velocity.                                    !
!                                                                                          !
! This is the MVP core: aerodynamics -> per-cohort surface energy balance -> CAS enthalpy/vapour/   !
! CO2. It carries the design's §3.5 correctness fix relative to the shipped canopy_air_update:       !
! the atm<->CAS scalar conductance is the PROFILE-FACTORED gah = rho*ustar*temp1 (heat/enthalpy),     !
! gaw = rho*ustar*temp2 (vapour), gac = rho_dmol*ustar*temp2 (CO2) -- NOT the bare rho*ustar that     !
! overstates the coupling ~5-10x -- and the VAPOUR twin is advanced IMPLICITLY like the other two.    !
!                                                                                          !
! The leaf temperature is DIAGNOSED from a linearized steady-state balance (a near-massless leaf     !
! equilibrates within a dt_fast; this is the physically-correct treatment and avoids the small-       !
! thermal-mass stiffness of time-integrating a prognostic leaf-energy store at large dt_fast). The     !
! prognostic leaf/wood-energy store, the soil column, and the photosynthesis/hydraulics fluxes are     !
! the next layer; absorbed radiation, stomatal conductance and the net biotic CO2 flux are prescribed. !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, cas_state_t,        &
                                     patch_biophys_t, alloc_patch_biophys
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_thermo,           only : cas_temp_of_enthalpy,                                     &
                                     sat_specific_humidity, d_sat_vapor_pressure_dt
   use meds_budget_check,     only : budget_t, budget_accumulate
   implicit none
   private

   public :: column_forcing_t, column_fast_step

   !----- Prescribed per-step forcing the higher layers (RT, photosynthesis, met) will supply. !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm = 0.0_wp  !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm      = 0.0_wp  !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm      = 400.0_wp!< [umol/mol] free-atmosphere CO2
      real(wp)              :: nee_biotic   = 0.0_wp  !< [umol/m2/s] net biotic CO2 source (Reco - GPP)
      real(wp), allocatable :: abs_sw(:), abs_lw(:)   !< [W/m2] absorbed SW / net LW per cohort (leaf)
      real(wp), allocatable :: gsw(:), fs_open(:)     !< [m/s], [-] stomatal conductance, open fraction
   end type column_forcing_t

contains

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch: aerodynamics -> per-cohort    !
   !  DIAGNOSTIC leaf energy balance -> CAS three-twin implicit update. Budgets accumulate the   !
   !  closed residuals. Cohort arrays are ordered BOTTOM(1) -> TOP(n).                           !
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, acfg, aenv, ageom, n, lai, wai, height, crown,         &
                               leaf_width, branch_diam, tparams, forc, bio, aero,              &
                               be_energy, be_water, be_co2)
      real(wp),                   intent(in)    :: dt_fast
      type(aero_cfg_t),           intent(in)    :: acfg
      type(aero_env_t),           intent(inout) :: aenv       !< can_* fields refreshed from CAS state
      type(aero_geom_t),          intent(in)    :: ageom
      integer(ik),                intent(in)    :: n
      real(wp),                   intent(in)    :: lai(n), wai(n), height(n), crown(n)
      real(wp),                   intent(in)    :: leaf_width(n), branch_diam(n)
      type(veg_thermal_params_t), intent(in)    :: tparams
      type(column_forcing_t),     intent(in)    :: forc
      type(patch_biophys_t),      intent(inout) :: bio
      type(aero_out_t),           intent(inout) :: aero       !< preallocated (alloc_aero_out)
      type(budget_t),             intent(inout) :: be_energy, be_water, be_co2

      real(wp)    :: tcas, qcas, press, rho, h_coeff, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: g_tr, le_ref, dtl, tl, h_i, le_i, transp_i
      real(wp)    :: coh_h, coh_le, coh_transp
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      integer(ik) :: i

      !----- 1. Refresh the aerodynamics env from the current CAS state, then solve. ---------!
      bio%cas%can_temp = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      tcas = bio%cas%can_temp ; qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      call canopy_aerodynamics(acfg, aenv, ageom, n, height, lai, crown, bio%leaf_temp,        &
                               bio%leaf_temp, leaf_width, branch_diam, aero)

      !----- 2. Per-cohort DIAGNOSTIC leaf steady-state -> sensible + latent flux into CAS. ---!
      qsat_c = sat_specific_humidity(tcas, press)
      dqdt   = 0.622_wp * press / max((press - 0.378_wp * sat_e(tcas))**2, tiny_num)           &
               * d_sat_vapor_pressure_dt(tcas)
      coh_h = 0.0_wp ; coh_le = 0.0_wp ; coh_transp = 0.0_wp
      do i = 1_ik, n
         !----- Linearized conductances (heat: 2 sides; transpiration in series w/ stomata). -!
         h_coeff  = tparams%effarea_heat * lai(i) * aero%leaf_gbh(i) * rho * cp_air            ! [W/m2/K]
         g_tr     = 0.0_wp
         if (aero%leaf_gbw(i) + forc%gsw(i) > tiny_num) then
            g_tr = tparams%effarea_transp * lai(i) * forc%fs_open(i)                           &
                   * aero%leaf_gbw(i) * forc%gsw(i) / (aero%leaf_gbw(i) + forc%gsw(i))         ! [m/s]
         end if
         lw_slope = 4.0_wp * tparams%leaf_emiss * stefan * tcas**3 * lai(i)                    ! [W/m2/K]
         le_slope = latent_heat_vap * rho * g_tr * dqdt                                        ! [W/m2/K]
         le_ref   = latent_heat_vap * rho * g_tr * (qsat_c - qcas)                             ! [W/m2] at tcas
         !----- Steady-state leaf temperature (near-massless leaf: solve, don't integrate). --!
         dtl = (forc%abs_sw(i) + forc%abs_lw(i) - le_ref) / max(h_coeff + le_slope + lw_slope, tiny_num)
         tl  = tcas + dtl
         bio%leaf_temp(i) = tl
         h_i      = h_coeff * dtl
         le_i     = le_ref + le_slope * dtl
         transp_i = le_i / latent_heat_vap                                                     ! [kg/m2/s]
         coh_h      = coh_h      + h_i
         coh_le     = coh_le     + le_i
         coh_transp = coh_transp + transp_i
      end do

      !----- 3. CAS three-twin update: IMPLICIT in the profile-factored atm exchange (§3.5). --!
      can_dmol = rho * (1.0_wp - qcas) / mmdry                       ! [mol_dryair/m3]
      wcap     = rho      * bio%cas%can_depth                        ! [kg/m2]  enthalpy+vapour capacity
      ccap     = can_dmol * bio%cas%can_depth                        ! [mol/m2] CO2 capacity
      gah      = rho      * aero%ustar * aero%temp1                  ! [kg/m2/s]  atm<->CAS enthalpy conductance
      gaw      = rho      * aero%ustar * aero%temp2                  ! [kg/m2/s]  atm<->CAS vapour conductance
      gac      = can_dmol * aero%ustar * aero%temp2                  ! [mol/m2/s] atm<->CAS CO2 conductance

      enth0 = bio%cas%can_enthalpy ; shv0 = qcas ; co20 = bio%cas%can_co2
      enth1 = (wcap*enth0 + dt_fast*(coh_h + coh_le + gah*forc%enthalpy_atm)) / (wcap + dt_fast*gah)
      shv1  = (wcap*shv0  + dt_fast*(coh_transp    + gaw*forc%shv_atm))       / (wcap + dt_fast*gaw)
      co21  = (ccap*co20  + dt_fast*(forc%nee_biotic + gac*forc%co2_atm))     / (ccap + dt_fast*gac)

      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- 4. Closed-budget accumulation (each twin: store change vs boundary fluxes). ------!
      call budget_accumulate(be_energy, wcap*enth0, wcap*enth1,                                &
                             coh_h + coh_le + gah*forc%enthalpy_atm, gah*enth1, dt_fast,        &
                             abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(be_water,  wcap*shv0, wcap*shv1,                                   &
                             coh_transp + gaw*forc%shv_atm, gaw*shv1, dt_fast,                  &
                             max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(be_co2,    ccap*co20, ccap*co21,                                   &
                             forc%nee_biotic + gac*forc%co2_atm, gac*co21, dt_fast,             &
                             abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
   end subroutine column_fast_step

   !----- Saturation vapour pressure [Pa] (Bolton; local mirror of the thermo form for dqdt). -!
   pure real(wp) function sat_e(t_k) result(esat)
      real(wp), intent(in) :: t_k
      real(wp) :: tc
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
   end function sat_e

end module meds_column_dynamics
