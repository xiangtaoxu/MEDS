!==========================================================================================!
! meds_biogeochem_types -- shared derived types + selector codes for the biogeochemistry     !
! domain: the CARBON / nutrient cycle of the ecosystem column, spanning the FAST canopy-air    !
! CO2 exchange and (later) the SLOW soil-carbon pools. Pure DATA, no methods -- the biogeochem   !
! analogue of meds_biophysics_types. Links src/shared (meds_kinds) ONLY.                          !
!                                                                                          !
! P0 (this cut) is the well-mixed canopy-air-space CO2 box (design MEDS_COLUMN_CO2_BALANCE_DESIGN !
! .md, section 4). Types:                                                                          !
!   * co2_opts_t           -- pre-extracted heterotrophic-respiration model selector + params.      !
!   * soil_carbon_t        -- the slow, stateful decomposable pool (prescribed at P0; daily at P2).   !
!   * cohort_co2_flux_t    -- the per-patch, ground-area cohort CO2 fluxes (from aggregate).           !
!   * column_co2_budget_t  -- NEE/NEP/loss2atm/storage + the closed-budget residual.                    !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: HR_Q10, HR_EXP_ED2, HR_DAMM
   public :: co2_opts_t, soil_carbon_t, column_co2_budget_t, cohort_co2_flux_t, damm_params_t

   !----- Heterotrophic-respiration MODEL selector codes (the co2_opts_t%hr_model field). --------!
   integer(ik), parameter :: HR_Q10     = 1_ik   !< Collatz/K13 Q10 q10**((T-T_ref)/10) x f_water x pool
   integer(ik), parameter :: HR_EXP_ED2 = 2_ik   !< ED2 scheme-0 capped exp min(1,exp(a*(T-T_sat))) x f_water x pool
   integer(ik), parameter :: HR_DAMM    = 3_ik   !< Davidson 2012 dual-Arrhenius/Michaelis-Menten (mechanistic moisture)

   !----- Slow, stateful per-patch soil-carbon pool (written DAILY, read-only in the fast loop). --!
   !      MVP: one lumped decomposable pool, prescribed at P0. Shaped so the CENTURY expansion is    !
   !      additive (structural/slow/microbial/passive land in a later slow-loop phase).               !
   type :: soil_carbon_t
      real(wp) :: fast_soil_carbon = 0.0_wp   !< [kgC/m2] lumped decomposable pool (MVP; prescribed at P0)
   end type soil_carbon_t

   !----- Pre-summed patch-ground cohort CO2 fluxes (filled by aggregate_cohort_co2_fluxes). ------!
   type :: cohort_co2_flux_t
      real(wp) :: gross_primary_prod  = 0.0_wp   !< [umol/m2/s] SUM a_gross*leaf_area*nplant (uptake, >=0)
      real(wp) :: leaf_respiration    = 0.0_wp   !< [umol/m2/s] SUM rd*leaf_area*nplant       (source, >=0)
      real(wp) :: stem_respiration    = 0.0_wp   !< [umol/m2/s] SUM stem_resp*nplant          (source, >=0)
      real(wp) :: root_respiration    = 0.0_wp   !< [umol/m2/s] SUM root_resp*nplant          (source, >=0)
      real(wp) :: growth_respiration  = 0.0_wp   !< [umol/m2/s] committed daily growth resp, amortized; MVP=0
      real(wp) :: storage_respiration = 0.0_wp   !< [umol/m2/s] committed daily storage resp, amortized; MVP=0
   end type cohort_co2_flux_t

   !----- Column CO2 budget + diagnostics (mirrors energy_flux_t / chydro_flux_t). ----------------!
   type :: column_co2_budget_t
      real(wp) :: nee      = 0.0_wp   !< [umol/m2/s] Reco - GPP  (atmospheric sign: >0 = source to atm)
      real(wp) :: nep      = 0.0_wp   !< [umol/m2/s] GPP - Reco  (= -nee; >0 = ecosystem uptake; DERIVED)
      real(wp) :: loss2atm = 0.0_wp   !< [umol/m2/s] CAS->free-atm venting = gatm_co2*(co2_new - co2_atm)
      real(wp) :: storage  = 0.0_wp   !< [umol/m2]   ccapcan*can_co2 (physical CAS CO2 inventory)
      real(wp) :: resid    = 0.0_wp   !< [umol/m2]   closed-budget residual (~0 by construction)
   end type column_co2_budget_t

   !----- DAMM diffusion parameters (Davidson et al. 2012 GCB 18:371, Table 2, Harvard Forest). ---!
   !      Used only when hr_model = HR_DAMM. bd/pd are provenance -- they enter the runtime SOLELY   !
   !      through theta_sat (porosity = 1 - bd/pd), so the kernel needs only soil_temp/theta/         !
   !      theta_sat/pool. p_soluble is 4.14e-4 (confirmed; NOT 0.0414).                                !
   type :: damm_params_t
      real(wp) :: alpha_sx  = 5.38e10_wp   !< [mgC cm-3 soil h-1] Arrhenius pre-exponential (calibrated; depth via depth_cm)
      real(wp) :: ea_sx     = 72.26_wp     !< [kJ/mol]            activation energy          (calibrated)
      real(wp) :: km_sx     = 9.95e-7_wp   !< [gC cm-3 soil]      soluble-C half-saturation  (weak prior)
      real(wp) :: km_o2     = 0.121_wp     !< [cm3 O2 cm-3 air]   O2 half-saturation         (weak prior)
      real(wp) :: p_soluble = 4.14e-4_wp   !< [-] soluble fraction of total soil C  (physically fixed)
      real(wp) :: d_liq     = 3.17_wp      !< [-] liquid-diffusion coefficient      (physically fixed)
      real(wp) :: d_gas     = 1.67_wp      !< [-] gas-diffusion coefficient         (physically fixed)
      real(wp) :: depth_cm  = 10.0_wp      !< [cm] effective respiring depth (SOC->conc AND flux depth-integral)
      real(wp) :: bd        = 0.80_wp      !< [g cm-3] bulk density   (provenance; enters only via theta_sat)
      real(wp) :: pd        = 2.52_wp      !< [g cm-3] particle density(provenance; 1 - bd/pd = 0.6825 = porosity)
   end type damm_params_t

   !----- Pre-extracted CO2/decomposition selectors + parameters (NOT the whole config). ---------!
   !      In-type defaults so the standalone kernels/tests compile pre-P3 (like soil_opts_t).       !
   type :: co2_opts_t
      integer(ik) :: hr_model              = HR_Q10        !< {HR_Q10, HR_EXP_ED2, HR_DAMM}
      real(wp)    :: rh_k_base             = 0.0_wp        !< [1/day]  effective decomposition rate (x pool)
      real(wp)    :: rh_q10                = 1.5_wp        !< [-]      Collatz/K13 Q10            (HR_Q10)
      real(wp)    :: rh_t_ref              = 288.15_wp     !< [K]      15 C Q10 reference         (HR_Q10)
      real(wp)    :: resp_temp_increase    = 0.0757_wp     !< [1/K]    ED2 scheme-0 slope        (HR_EXP_ED2)
      real(wp)    :: resp_temp_ref         = 318.15_wp     !< [K]      45 C saturation           (HR_EXP_ED2)
      real(wp)    :: resp_opt_water        = 0.8938_wp     !< [-]      moisture optimum (relative)
      real(wp)    :: resp_water_below_opt  = 5.0786_wp     !< [-]      dry-side exponential slope
      real(wp)    :: resp_water_above_opt  = 4.5139_wp     !< [-]      wet-side (anoxia) exponential slope
      real(wp)    :: co2_atm_ref           = 400.0_wp      !< [umol/mol] fixed atm CO2 when not met-forced
      real(wp)    :: rtol = 1.0e-8_wp                      !< [-]       relative closure tolerance
      real(wp)    :: atol = 1.0e-3_wp                      !< [umol/m2] absolute closure floor
      logical     :: debug_error = .false.
      type(damm_params_t) :: damm                         !< DAMM diffusion parameters (HR_DAMM only)
   end type co2_opts_t

end module meds_biogeochem_types
