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
! Leaf temperature is DIAGNOSED from a linearized steady-state balance. CARBON is now REAL: leaf gas     !
! exchange gives GPP + stomatal conductance (driving transpiration) + leaf Rd; stem + fine-root           !
! maintenance respiration (autotrophic) and heterotrophic Rh assemble NEE = (Rd+stem+root) + Rh - GPP     !
! for the CAS CO2 twin. Plant HYDRAULICS is now REAL too: each cohort carries prognostic node water       !
! potentials psi(NODE_LEAF/WOOD) in patch_biophys_t; plant_water_flux advances them from the realized     !
! (supply-limited) transpiration demand and the root-weighted soil psi, and the updated psi_leaf feeds     !
! NEXT step's leaf gas exchange -- the soil -> plant -> stomata drought feedback (lagged one dt_fast).     !
! Still prescribed (next layers): absorbed radiation, precip. Canopy interception (leaf film), the         !
! stepper hook, and cross-demography persistence remain to wire.                                            !
!                                                                                          !
! WHOLE-COLUMN CONSERVATION -- verified by budg%whole_energy / budg%whole_water (Δ of ALL stores vs the     !
! true boundary fluxes; these CATCH cross-seam leaks the per-kernel budgets miss). Water-borne enthalpy    !
! is transported consistently: the CAS latent uses enthalpy_vapor(tl) (matching the CAS inverter + ground); !
! the soil sheds the transpiration water's liquid enthalpy via root_heat_sink; infiltration/drainage water  !
! carry internal_energy_liquid across the soil boundaries. INTER-LAYER advective heat: the hydrology kernel  !
! now EXPOSES the time-mean per-face Darcy flux (hflux%w_flux), and soil_energy_flux can advect the liquid   !
! enthalpy on it -- an OPT-IN coupling (cfg%advect_soil_heat, default OFF). It reconciles the soil moisture   !
! <-> energy coupling that the standalone kernels never exercised (their unit tests forced moisture constant):!
! with the internal-energy soil store referenced at tsupercool_liq (~193 K), a per-step Δθ carries a large    !
! liquid enthalpy, so turning the advection ON shifts the coupled soil/CAS temperatures NOTABLY. It still      !
! conserves the whole-column energy to machine precision and keeps soil temps bounded, but the MAGNITUDE      !
! awaits validation against a reference (ED2) with real met forcing -- hence gated OFF by default. A          !
! SUPPLY-limited leaf is not re-solved for the extra sensible (small drought residual). The stepper hook +    !
! cross-demography persistence remain to wire.                                                                !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap, rho_h2o, r_gas
   use meds_config,           only : meds_config_t
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, chydro_forcing_t, chydro_flux_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_column_energy,    only : soil_energy_flux
   use meds_column_hydrology, only : column_hydrology_flux
   use meds_plant_interface,  only : leaf_env_t, leaf_flux_t, leaf_gas_exchange,               &
                                     wood_env_t, wood_params_t, wood_flux_t, stem_maintenance_respiration, &
                                     root_env_t, root_params_t, root_flux_t, fine_root_maintenance_respiration, &
                                     hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     plant_water_flux, N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_column_co2,       only : heterotrophic_respiration_flux
   use meds_biogeochem_types, only : co2_opts_t
   use meds_thermo,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     d_sat_vapor_pressure_dt, enthalpy_vapor, internal_energy_liquid
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
   implicit none
   private

   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: alloc_column_cohort, column_fast_step

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal    !< leaf/wood thermal params
      type(soil_params_t)         :: soil           !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal   !< soil thermal texture
      type(energy_opts_t)         :: energy         !< soil-thermal solver options
      type(soil_opts_t)           :: hydro          !< soil-water (Richards) solver options
      type(wood_params_t)         :: wood           !< stem-respiration parameters
      type(root_params_t)         :: root           !< fine-root-respiration parameters
      type(co2_opts_t)            :: co2            !< heterotrophic-respiration options
      type(hydro_params_t)        :: hydro_p        !< plant-hydraulics parameters (PV curves, vulnerability)
      type(hydro_opts_t)          :: hydro_o        !< plant-hydraulics solver options
      real(wp)                    :: fast_soil_carbon = 5.0_wp   !< [kgC/m2] decomposable soil-C pool (prescribed, MVP)
      real(wp)                    :: rhizo_cond       = 5.0e-4_wp !< [kg/s/MPa] soil->root conductance (prescribed, MVP)
      logical                     :: advect_soil_heat = .false.  !< opt-in: advect liquid enthalpy on the interior
                                                                 !< per-layer Darcy flux (moisture<->energy coupling)
   end type column_config_t

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply; photosynthesis/    !
   !      respiration/NEE are now computed from the plant kernels (no longer prescribed).     !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: precip        = 0.0_wp   !< [kg/m2/s] ground-reaching rainfall (interception deferred)
      real(wp)              :: par_per_w     = 2.1_wp   !< [umol photon / (W SW absorbed)] SW->incident-PAR factor (MVP)
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW / net LW per cohort (leaf)
   end type column_forcing_t

   !----- The per-patch conservation budgets (one place; the driver accumulates the closed resids).!
   !      The per-kernel budgets close BY CONSTRUCTION; whole_energy/whole_water are the CROSS-      !
   !      seam column totals (Δ all stores vs the true boundary fluxes) that actually catch leaks.   !
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   !< [umol/m2/s] last-step diagnostics
   end type column_budget_t

contains

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(coh, n)
      type(column_cohort_t), intent(out) :: coh
      integer(ik),           intent(in)  :: n
      coh%n = n
      allocate(coh%pft(n), coh%lai(n), coh%wai(n), coh%height(n), coh%crown(n),                &
               coh%leaf_width(n), coh%branch_diam(n), coh%leaf_area(n), coh%nplant(n),         &
               coh%dbh(n), coh%broot(n), coh%bleaf(n), coh%bsap(n), coh%sap_area(n))
      coh%pft = 1_ik
      coh%lai = 0.0_wp ; coh%wai = 0.0_wp ; coh%height = 0.0_wp ; coh%crown = 1.0_wp
      coh%leaf_width = 0.04_wp ; coh%branch_diam = 0.02_wp
      coh%leaf_area = 0.0_wp ; coh%nplant = 0.0_wp ; coh%dbh = 0.0_wp ; coh%broot = 0.0_wp
      coh%bleaf = 0.0_wp ; coh%bsap = 0.0_wp ; coh%sap_area = 0.0_wp
   end subroutine alloc_column_cohort

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch. Cohort arrays BOTTOM(1)->TOP. !
   !  Leaf gas exchange (real GPP + stomata + leaf Rd), stem/root maintenance respiration and    !
   !  heterotrophic Rh feed a physically-decomposed NEE = (Rd_leaf + stem + root) + Rh - GPP.    !
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh, &
                               leaf_resp_coh, stem_resp_coh, root_resp_coh)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg          !< PFT traits for leaf gas exchange
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv         !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero         !< preallocated (alloc_aero_out)
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:)   !< [umol CO2/plant/s] per-cohort GROSS GPP (fast->slow)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration (fast->slow)
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp (fast->slow)
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp (fast->slow)

      type(chydro_forcing_t) :: hforc
      type(chydro_flux_t)    :: hflux
      type(energy_forcing_t) :: eforc
      type(energy_flux_t)    :: sflux
      type(leaf_env_t)       :: lenv
      type(leaf_flux_t)      :: lf
      type(wood_env_t)       :: wenv
      type(wood_flux_t)      :: wf
      type(root_env_t)       :: renv
      type(root_flux_t)      :: rf
      type(hydro_env_t)      :: henv
      type(hydro_flux_t)     :: hfx
      real(wp)               :: transp_c(coh%n)     !< [kg/m2 ground/s] per-cohort transpiration demand (automatic)
      real(wp)    :: soil_psi_root
      real(wp)    :: tcas, qcas, press, rho, h_coeff, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: g_tr, le_ref, dtl, tl, h_i, le_i, transp_i, gsw_ms, e_air, rho_mol
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      real(wp)    :: gpp, ra_leaf, ra_stem, ra_root, rh, nee_biotic, soil_temp_root, theta_mean
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol, src_enth, src_vap, src_frac
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, e_in, e_out, w_in, w_out
      integer(ik) :: i, n, nsl, k

      n = coh%n ; nsl = ccfg%soil%n_active

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
      call canopy_aerodynamics(ccfg%aero, aenv, ageom, n, coh%height, coh%lai, coh%crown,      &
                               bio%leaf_temp, bio%leaf_temp, coh%leaf_width, coh%branch_diam, aero)

      !----- Root-weighted soil temperature + column-mean moisture (root / heterotrophic resp). !
      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- 2. Per cohort: LEAF gas exchange (real GPP / gs / Rd) + DIAGNOSTIC leaf energy      !
      !         balance (transpiration via aero-gb in series with the REAL gs) + stem/root        !
      !         maintenance respiration. CAS latent uses enthalpy_vapor(tl); the soil sheds the    !
      !         water's liquid enthalpy (coh_qsoil). NEE is assembled after the loop.              !
      qsat_c = sat_specific_humidity(tcas, press)
      dqdt   = 0.622_wp * press / max((press - 0.378_wp * sat_e(tcas))**2, tiny_num)           &
               * d_sat_vapor_pressure_dt(tcas)
      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      do i = 1_ik, n
         !----- Leaf gas exchange: incident PAR, leaf-to-air VPD, CAS CO2, molar boundary gb. --!
         rho_mol       = press / (r_gas * bio%leaf_temp(i))                     ! [mol/m3] molar air density
         e_air         = qcas * press / (0.622_wp + 0.378_wp * qcas)            ! [Pa] canopy-air vapour pressure
         lenv%par      = forc%abs_sw(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w
         lenv%leaf_temp = bio%leaf_temp(i)
         lenv%vpd      = max(sat_e(bio%leaf_temp(i)) - e_air, 0.0_wp)
         lenv%ca       = bio%cas%can_co2
         lenv%pressure = press
         lenv%psi_leaf = bio%psi(NODE_LEAF, i)                                  ! lagged plant water status (hydraulics)
         lenv%gb       = aero%leaf_gbw(i) * rho_mol                             ! m/s -> mol H2O/m2/s
         call leaf_gas_exchange(lenv, cfg, coh%pft(i), lf)
         gsw_ms        = lf%gs / max(rho_mol, tiny_num)                         ! mol/m2/s -> m/s
         gpp           = gpp     + lf%a_gross * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = lf%a_gross * coh%leaf_area(i)          ! [umol/plant/s] per-plant gross

         ra_leaf       = ra_leaf + lf%rd      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)   ! [umol/plant/s] per-plant leaf Rd
         !----- Diagnostic leaf energy balance (transpiration = aero-gb in series with real gs). !
         h_coeff = ccfg%veg_thermal%effarea_heat * coh%lai(i) * aero%leaf_gbh(i) * rho * cp_air
         g_tr    = 0.0_wp
         if (aero%leaf_gbw(i) + gsw_ms > tiny_num) then
            g_tr = ccfg%veg_thermal%effarea_transp * coh%lai(i)                                &
                   * aero%leaf_gbw(i) * gsw_ms / (aero%leaf_gbw(i) + gsw_ms)
         end if
         lw_slope = 4.0_wp * ccfg%veg_thermal%leaf_emiss * stefan * tcas**3 * coh%lai(i)
         le_slope = latent_heat_vap * rho * g_tr * dqdt
         le_ref   = latent_heat_vap * rho * g_tr * (qsat_c - qcas)
         dtl = (forc%abs_sw(i) + forc%abs_lw(i) - le_ref) / max(h_coeff + le_slope + lw_slope, tiny_num)
         tl  = tcas + dtl
         bio%leaf_temp(i) = tl
         h_i      = h_coeff * dtl
         le_i     = le_ref + le_slope * dtl
         transp_i = le_i / latent_heat_vap
         transp_c(i) = transp_i                                                          ! per-cohort demand (for hydraulics)
         coh_h      = coh_h      + h_i
         coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)                          ! CAS latent (vapour enthalpy)
         coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)      ! liquid enthalpy the soil sheds
         coh_transp = coh_transp + transp_i
         coh_rnet   = coh_rnet   + (forc%abs_sw(i) + forc%abs_lw(i) - lw_slope * dtl)     ! NET leaf radiation (post LW emission)
         !----- Autotrophic maintenance respiration: stem + fine root (per plant -> per m2). ---!
         wenv%wood_temp = bio%leaf_temp(i) ; wenv%dbh = coh%dbh(i) ; wenv%height = coh%height(i)
         wenv%wai = coh%wai(i) ; wenv%nplant = coh%nplant(i)
         call stem_maintenance_respiration(wenv, ccfg%wood, wf)
         renv%soil_temp = soil_temp_root ; renv%broot = coh%broot(i)
         call fine_root_maintenance_respiration(renv, ccfg%root, rf)
         ra_stem = ra_stem + wf%stem_resp * coh%nplant(i)
         ra_root = ra_root + rf%root_resp * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp   ! [umol/plant/s] already per-plant
         if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp   ! [umol/plant/s] already per-plant
      end do

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP. --------------!
      rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,   &
                                          ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

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

      !----- 3b. Plant HYDRAULICS: advance each cohort's psi from the realized transpiration     !
      !          demand + the root-weighted soil psi. The updated psi_leaf feeds NEXT step's      !
      !          leaf gas exchange (the soil -> plant -> stomata water-stress feedback).          !
      soil_psi_root = 0.0_wp
      do k = 1_ik, nsl
         soil_psi_root = soil_psi_root + hflux%psi_soil(k) * ccfg%soil%root_frac(k)
      end do
      do i = 1_ik, n
         henv%transp     = transp_c(i) * src_frac / max(coh%nplant(i), tiny_num)   ! [kg/plant/s]
         henv%soil_psi   = soil_psi_root
         henv%rhizo_cond = ccfg%rhizo_cond
         henv%bleaf      = coh%bleaf(i) ; henv%bsap = coh%bsap(i) ; henv%broot = coh%broot(i)
         henv%sap_area   = coh%sap_area(i) ; henv%height = coh%height(i) ; henv%leaf_area = coh%leaf_area(i)
         call plant_water_flux(henv, ccfg%hydro_p, ccfg%hydro_o, dt_fast, bio%psi(:,i), hfx)
      end do

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
      co21  = (ccap*co20  + dt_fast*(nee_biotic + gac*forc%co2_atm)) / (ccap + dt_fast*gac)

      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- 6. Advective water enthalpy across the soil BOUNDARIES (top infiltration + rain temp; !
      !         bottom drainage; surface runoff) is applied here; the INTERIOR inter-layer advection  !
      !         is handed to soil_energy_flux via eforc%w_flux (the hydrology kernel now exposes the   !
      !         time-mean per-face Darcy flux). Sign flip: hydrology is downward-positive, the energy   !
      !         kernel is upward-positive. Root-uptake liquid enthalpy is shed via root_heat_sink.       !
      bio%soil_e%soil_energy(1)   = bio%soil_e%soil_energy(1)                                       &
           + (hflux%infiltration * internal_energy_liquid(rain_temp)                                &
              - hflux%runoff_surf * internal_energy_liquid(t_ground)) * dt_fast / ccfg%soil%dz(1)
      bio%soil_e%soil_energy(nsl) = bio%soil_e%soil_energy(nsl)                                     &
           - hflux%drainage * internal_energy_liquid(t_bot) * dt_fast / ccfg%soil%dz(nsl)
      eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
      eforc%soil_water(1:nsl)     = bio%soil_w%theta(1:nsl)
      if (ccfg%advect_soil_heat) then
         eforc%w_flux(1:nsl)      = -hflux%w_flux(1:nsl)      ! down-positive (hydro) -> up-positive (energy)
      else
         eforc%w_flux(1:nsl)      = 0.0_wp                    ! interior advection lumped (validated baseline)
      end if
      eforc%root_heat_sink(1:nsl) = coh_qsoil * ccfg%soil%root_frac(1:nsl)     ! §3 fix: shed transpiration-water enthalpy
      call soil_energy_flux(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)

      !----- 7. Per-kernel closed budgets (each closes by construction). ----------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm, &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(budg%cas_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,        &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(budg%cas_co2,    ccap*co20, ccap*co21, nee_biotic + gac*forc%co2_atm,&
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
