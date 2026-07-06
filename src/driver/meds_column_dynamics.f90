!==========================================================================================!
! meds_column_dynamics -- the fast-timescale (dt_fast) integrator that couples the fast          !
! biophysics processes into one sub-daily operator-split sweep (design Part II). The canopy-air-  !
! space (CAS) twins are the shared coupling reservoir: each dt_fast the aerodynamics kernel sets   !
! the conductances; the leaf and GROUND surface fluxes feed the CAS; the soil WATER column         !
! (Richards) and soil THERMAL column (heat) are advanced; and the three CAS twins (enthalpy /       !
! specific humidity / CO2) are advanced IMPLICITLY in the atmospheric-exchange term, gated by u*.    !
!                                                                                          !
! Coupled so far: aerodynamics -> {diagnostic leaf balance, ground balance} -> soil WATER column     !
! (infiltration / DSL soil-evap / root uptake / drainage) -> CAS enthalpy/vapour/CO2 (implicit) ->    !
! soil THERMAL column (implicit BE-Thomas), which reads the just-updated soil moisture. Carries the    !
! §3.5 fix (atm<->CAS conductance = profile-factored rho*ustar*temp1/temp2) AND the §3.6 GROUND-        !
! evaporation single-authority: the hydrology kernel's DSL/alpha_soil soil_evap is THE ground latent    !
! flux -- it drives the CAS vapour twin AND the ground energy balance's LE (no double-count).           !
!                                                                                          !
! Leaf temperature is DIAGNOSED from a linearized steady-state balance. Prescribed for now (next        !
! layers): absorbed radiation, stomatal conductance, net biotic CO2, precip. Canopy interception (leaf  !
! film), photosynthesis/hydraulics (real GPP/E), the stepper hook, and cross-demography persistence     !
! remain to wire.                                                                                        !
!                                                                                          !
! WHOLE-COLUMN CONSERVATION -- verified by budg%whole_energy / budg%whole_water (Δ of ALL stores vs the     !
! true boundary fluxes; these CATCH cross-seam leaks the per-kernel budgets miss). Water-borne enthalpy    !
! is transported consistently: the CAS latent uses enthalpy_vapor(tl) (matching the CAS inverter + ground); !
! the soil sheds the transpiration water's liquid enthalpy via root_heat_sink; infiltration/drainage water  !
! carry internal_energy_liquid across the soil boundaries. Remaining approximations (small): infiltration   !
! enthalpy is lumped at the top layer (the hydrology kernel does not yet expose per-layer water fluxes for   !
! exact inter-layer advective heat), and a SUPPLY-limited leaf is not re-solved for the extra sensible      !
! (closes with the plant-hydraulics layer) -- so the whole-column energy budget is exact for well-watered    !
! soil and carries a small residual under drought. The stepper hook + cross-demography persistence remain.   !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap, rho_h2o
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, chydro_forcing_t, chydro_flux_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_column_energy,    only : soil_energy_flux
   use meds_column_hydrology, only : column_hydrology_flux
   use meds_thermo,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     d_sat_vapor_pressure_dt, enthalpy_vapor, internal_energy_liquid
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
   implicit none
   private

   public :: column_config_t, column_forcing_t, column_budget_t, column_fast_step

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero          !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal  !< leaf/wood thermal params
      type(soil_params_t)         :: soil         !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal !< soil thermal texture
      type(energy_opts_t)         :: energy       !< soil-thermal solver options
      type(soil_opts_t)           :: hydro        !< soil-water (Richards) solver options
   end type column_config_t

   !----- Prescribed per-step forcing the higher layers (RT, photosynthesis, met) will supply. !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: nee_biotic    = 0.0_wp   !< [umol/m2/s] net biotic CO2 source (Reco - GPP)
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: precip        = 0.0_wp   !< [kg/m2/s] ground-reaching rainfall (interception deferred)
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW / net LW per cohort (leaf)
      real(wp), allocatable :: gsw(:), fs_open(:)       !< [m/s], [-] stomatal conductance, open fraction
   end type column_forcing_t

   !----- The per-patch conservation budgets (one place; the driver accumulates the closed resids).!
   !      The per-kernel budgets close BY CONSTRUCTION; whole_energy/whole_water are the CROSS-      !
   !      seam column totals (Δ all stores vs the true boundary fluxes) that actually catch leaks.   !
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
   end type column_budget_t

contains

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch. Cohort arrays BOTTOM(1)->TOP. !
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, ccfg, aenv, ageom, n, lai, wai, height, crown,         &
                               leaf_width, branch_diam, forc, bio, aero, budg)
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
      type(column_budget_t),   intent(inout) :: budg

      type(chydro_forcing_t) :: hforc
      type(chydro_flux_t)    :: hflux
      type(energy_forcing_t) :: eforc
      type(energy_flux_t)    :: sflux
      real(wp)    :: tcas, qcas, press, rho, h_coeff, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: g_tr, le_ref, dtl, tl, h_i, le_i, transp_i
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol, src_enth, src_vap, src_frac
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, e_in, e_out, w_in, w_out
      integer(ik) :: i, nsl, k

      nsl = ccfg%soil%n_active

      !----- Snapshot start-of-step SOIL stores (for the whole-column budgets). --------------!
      e_soil0 = 0.0_wp ; w_soil0 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil0 = e_soil0 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil0 = w_soil0 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do

      !----- 1. Refresh the aerodynamics env from the current CAS state, then solve. ---------!
      bio%cas%can_temp = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      tcas = bio%cas%can_temp ; qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air
      t_ground = bio%soil_e%soil_temp(1) ; t_bot = bio%soil_e%soil_temp(nsl) ; rain_temp = tcas
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      call canopy_aerodynamics(ccfg%aero, aenv, ageom, n, height, lai, crown, bio%leaf_temp,   &
                               bio%leaf_temp, leaf_width, branch_diam, aero)

      !----- 2. Per-cohort DIAGNOSTIC leaf steady-state -> sensible + latent flux into CAS. ---!
      !         The CAS latent is credited with the VAPOUR enthalpy enthalpy_vapor(tl) (matching  !
      !         the CAS inverter + the ground term); the leaf balance keeps latent_heat_vap. The   !
      !         difference (coh_qsoil) is the water's liquid enthalpy, shed by the soil (§3 fix).   !
      qsat_c = sat_specific_humidity(tcas, press)
      dqdt   = 0.622_wp * press / max((press - 0.378_wp * sat_e(tcas))**2, tiny_num)           &
               * d_sat_vapor_pressure_dt(tcas)
      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
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
         coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)                          ! CAS latent (vapour enthalpy)
         coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)      ! liquid enthalpy the soil sheds
         coh_transp = coh_transp + transp_i
         coh_rnet   = coh_rnet   + (forc%abs_sw(i) + forc%abs_lw(i) - lw_slope * dtl)     ! NET leaf radiation (post LW emission)
      end do

      !----- 3. Soil WATER column: infiltration + DSL soil-evap + root uptake + drainage. -----!
      hforc%precip_ground          = forc%precip
      hforc%root_uptake(1:nsl)     = coh_transp * ccfg%soil%root_frac(1:nsl)
      hforc%t_ground               = t_ground
      hforc%q_air                  = qcas
      hforc%rho_air                = rho
      hforc%r_aero                 = 1.0_wp / max(aero%ggnet, tiny_num)     ! §3.6: r_aero = 1/ggnet
      call column_hydrology_flux(bio%soil_w, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
      soil_evap = hflux%soil_evap                                          ! §3.6: THE ground latent authority

      !----- Reconcile transpiration SUPPLY vs demand: the CAS gains only the water the soil    !
      !      actually gave up (no water creation under stress); latent enthalpy tracks the water. !
      src_frac = 1.0_wp
      if (coh_transp > tiny_num) src_frac = min(1.0_wp, hflux%uptake_total / coh_transp)
      coh_qw     = coh_qw     * src_frac
      coh_qsoil  = coh_qsoil  * src_frac
      coh_transp = coh_transp * src_frac

      !----- 4. GROUND surface energy: sensible + the authoritative (hydrology) latent. -------!
      h_ground  = aero%ggnet * rho * cp_air * (t_ground - tcas)
      le_ground = soil_evap * enthalpy_vapor(t_ground)
      g_top     = forc%abs_sw_ground + forc%abs_lw_ground - h_ground - le_ground

      !----- 5. CAS three-twin update: IMPLICIT in the profile-factored atm exchange (§3.5). --!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      wcap     = rho      * bio%cas%can_depth
      ccap     = can_dmol * bio%cas%can_depth
      gah      = rho      * aero%ustar * aero%temp1
      gaw      = rho      * aero%ustar * aero%temp2
      gac      = can_dmol * aero%ustar * aero%temp2
      src_enth = coh_h + coh_qw + h_ground + le_ground                    ! [W/m2]  sensible + latent (vapour enthalpy)
      src_vap  = coh_transp + soil_evap                                   ! [kg/m2/s] leaf transp + soil evap

      enth0 = bio%cas%can_enthalpy ; shv0 = qcas ; co20 = bio%cas%can_co2
      enth1 = (wcap*enth0 + dt_fast*(src_enth + gah*forc%enthalpy_atm)) / (wcap + dt_fast*gah)
      shv1  = (wcap*shv0  + dt_fast*(src_vap  + gaw*forc%shv_atm))       / (wcap + dt_fast*gaw)
      co21  = (ccap*co20  + dt_fast*(forc%nee_biotic + gac*forc%co2_atm)) / (ccap + dt_fast*gac)

      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- 6. Advective water enthalpy across the soil boundaries (top infiltration + rain temp; !
      !         bottom drainage; surface runoff), THEN the soil thermal solve reads the corrected   !
      !         energy + the just-updated moisture. Root-uptake liquid enthalpy is shed via         !
      !         eforc%root_heat_sink. (Per-layer infiltration/inter-layer enthalpy is lumped at the  !
      !         top -- an approximation pending the hydrology kernel exposing per-layer water fluxes.)!
      bio%soil_e%soil_energy(1)   = bio%soil_e%soil_energy(1)                                       &
           + (hflux%infiltration * internal_energy_liquid(rain_temp)                                &
              - hflux%runoff_surf * internal_energy_liquid(t_ground)) * dt_fast / ccfg%soil%dz(1)
      bio%soil_e%soil_energy(nsl) = bio%soil_e%soil_energy(nsl)                                     &
           - hflux%drainage * internal_energy_liquid(t_bot) * dt_fast / ccfg%soil%dz(nsl)
      eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
      eforc%soil_water(1:nsl)     = bio%soil_w%theta(1:nsl)
      eforc%w_flux(1:nsl)         = 0.0_wp
      eforc%root_heat_sink(1:nsl) = coh_qsoil * ccfg%soil%root_frac(1:nsl)     ! §3 fix: shed transpiration-water enthalpy
      call soil_energy_flux(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)

      !----- 7. Per-kernel closed budgets (each closes by construction). ----------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm, &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(budg%cas_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,        &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(budg%cas_co2,    ccap*co20, ccap*co21, forc%nee_biotic + gac*forc%co2_atm,&
                             gac*co21, dt_fast, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_energy, sflux%energy_resid, abs(g_top)*dt_fast + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_water,  hflux%mass_resid,   1.0_wp,             1.0e-6_wp, 1.0e-4_wp)

      !----- 7b. WHOLE-COLUMN budgets: Δ(all stores) vs the TRUE boundary fluxes (catches leaks). !
      e_soil1 = 0.0_wp ; w_soil1 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil1 = e_soil1 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil1 = w_soil1 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- Water: precip IN; drainage + runoff + atm-vapour OUT. --------------------------!
      w_in  = forc%precip
      w_out = hflux%drainage + hflux%runoff_surf + gaw * (shv1 - forc%shv_atm)
      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0, w_soil1 + wcap*shv1, w_in, w_out, &
                             dt_fast, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      !----- Energy: net radiation + rain enthalpy IN; atm exchange + drainage/runoff enthalpy OUT.!
      e_in  = coh_rnet + forc%abs_sw_ground + forc%abs_lw_ground                                     &
              + hflux%infiltration * internal_energy_liquid(rain_temp)   ! energy that reached the SOIL store
      e_out = gah * (enth1 - forc%enthalpy_atm) + hflux%drainage * internal_energy_liquid(t_bot)      &
              + hflux%runoff_surf * internal_energy_liquid(t_ground)
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0, e_soil1 + wcap*enth1, e_in, e_out, &
                             dt_fast, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)
   end subroutine column_fast_step

   !----- Track a kernel's own closed-budget residual into a budget_t (worst + fail count). ---!
   pure subroutine track_resid(b, resid, scale, rtol, atol)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: resid, scale, rtol, atol
      b%n_check = b%n_check + 1_ik
      b%resid   = resid
      b%worst   = max(b%worst, abs(resid))
      if (.not. closure_ok(resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine track_resid

   !----- Saturation vapour pressure [Pa] (Bolton; local mirror of the thermo form for dqdt). -!
   pure real(wp) function sat_e(t_k) result(esat)
      real(wp), intent(in) :: t_k
      real(wp) :: tc
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
   end function sat_e

end module meds_column_dynamics
