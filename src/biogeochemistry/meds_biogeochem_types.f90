!==========================================================================================!
! meds_biogeochem_types -- shared derived types + selector codes for the biogeochemistry     !
! domain: the SLOW soil-carbon / nutrient cycle of the ecosystem column. Pure DATA +           !
! parameters, no methods -- the biogeochem analogue of meds_biophysics_types. Links              !
! src/shared (meds_kinds, meds_column_state_types) ONLY.                                          !
!                                                                                          !
! (The FAST canopy-air-space CO2 exchange -- a sub-daily biophysical diffusion/venting process --  !
! moved to meds_biophysics_types + meds_cas_biophysics under src/biophysics; this module now holds      !
! only the SLOW soil-carbon pools.)                                                                 !
!                                                                                          !
! SLOW soil carbon (design MEDS_BIOGEOCHEMISTRY_DESIGN.md): the CENTURY-family multi-pool soil-    !
! carbon state advanced daily by the carbon matrix ODE dX/dt = B*I + A*xi*K*X. `soil_carbon_t` (the   !
! 7-pool vector + lignin sub-state + optional N) is DEFINED in meds_column_state_types (shared/state, !
! docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part I §8.1) and RE-EXPORTED here; field                !
! `fast_soil_carbon` KEEPS its name/index (2) so meds_cas_biophysics -- which reads it as a BARE       !
! SCALAR in the fast loop -- compiles unchanged. `decomp_opts_t` + its selector codes are similarly    !
! DEFINED in meds_biogeochem_opts (shared/config, MEDS_SLOW_DYNAMICS_DESIGN.md Part II B0) so           !
! meds_config can carry it with no shared->biogeochemistry edge, and RE-EXPORTED here. Slow types      !
! owned HERE: litter_input_t / soilc_audit_t / soilc_diag_t, and the fixed pool count `n_soil_pool` +   !
! `IP_*`.                                                                                          !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   use meds_column_state_types, only : soil_carbon_t
   use meds_biogeochem_opts, only : decomp_opts_t, DECOMP_STEP_EULER, DECOMP_STEP_EXPM,             &
                                    DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   implicit none
   private

   public :: soil_carbon_t
   !----- Slow soil-carbon matrix additions (P0). ------------------------------------------------!
   public :: n_soil_pool
   public :: IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   public :: DECOMP_STEP_EULER, DECOMP_STEP_EXPM
   public :: DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   public :: decomp_opts_t, litter_input_t, soilc_audit_t, soilc_diag_t
   !----- Fast heterotrophic-respiration selectors + params (kernels in meds_soil_biogeochem). ----!
   public :: HR_Q10, HR_EXP_ED2, HR_DAMM
   public :: co2_opts_t, damm_params_t

   !----- Heterotrophic-respiration MODEL selector codes (the co2_opts_t%hr_model field). --------!
   integer(ik), parameter :: HR_Q10     = 1_ik   !< Collatz/K13 Q10 q10**((T-T_ref)/10) x f_water x pool
   integer(ik), parameter :: HR_EXP_ED2 = 2_ik   !< ED2 scheme-0 capped exp min(1,exp(a*(T-T_sat))) x f_water x pool
   integer(ik), parameter :: HR_DAMM    = 3_ik   !< Davidson 2012 dual-Arrhenius/Michaelis-Menten (mechanistic moisture)

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

   !----- Pre-extracted fast-Rh selectors + parameters. In-type defaults so the standalone kernels/ !
   !      tests compile; the env-scalar params share semantics with decomp_opts_t (matched chemistry).!
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

   !==========================================================================================!
   !  SLOW soil-carbon matrix: fixed pool count + index parameters (the matrix state vector X).   !
   !  Ordered litter -> SOM -> passive so the default-scheme ACTIVE block of A is lower-triangular  !
   !  (back-substitution / GE solvable). The 7-pool ceiling matches ED2's CENTURY set exactly.      !
   !==========================================================================================!
   integer(ik), parameter :: n_soil_pool = 7_ik   !< length of the carbon matrix state vector X

   integer(ik), parameter :: IP_FAST_GRND   = 1_ik   !< metabolic litter, above-ground (flammable)
   integer(ik), parameter :: IP_FAST_SOIL   = 2_ik   !< metabolic litter, below-ground  (the P0 lumped pool)
   integer(ik), parameter :: IP_STRUCT_GRND = 3_ik   !< structural litter + CWD, above-ground
   integer(ik), parameter :: IP_STRUCT_SOIL = 4_ik   !< structural litter + CWD, below-ground
   integer(ik), parameter :: IP_MICR        = 5_ik   !< microbial SOM        (active only in scheme 5)
   integer(ik), parameter :: IP_SLOW        = 6_ik   !< slow / humified SOM
   integer(ik), parameter :: IP_PASSIVE     = 7_ik   !< passive SOM          (active only in scheme 5)

   !==========================================================================================!
   !  Slow, stateful per-patch soil-carbon pools (written DAILY, read-only in the fast loop).       !
   !  NOW DEFINED in meds_column_state_types (src/shared/state) -- re-exported here (its conceptual  !
   !  home) so meds_soil_biogeochem's kernels + meds_cas_biophysics compile unchanged. It lives in    !
   !  shared/state, not here, so `patch_block` (src/core, links `shared` ONLY) can carry a per-patch  !
   !  soil_carbon_t with no core->biogeochemistry library edge (the same reason cas_state_t/           !
   !  soil_column_t/soil_energy_column_t/snow_column_t live there). The ACTIVE pool count is set by     !
   !  decomp_scheme (default 3-active: idx 5,7 inert).                                                  !
   !==========================================================================================!

   !----- Per-patch, per-day litter input -- ALREADY PARTITIONED to pool destinations (§2.2, §5.4). !
   !      The driver sums each cohort's leaf/fineroot/wood/storage necromass into these bins USING    !
   !      that cohort's PFT f_labile_leaf/f_labile_stem and agf; the kernel maps them straight onto u. !
   !      Units [kgC/m2/day].                                                                          !
   type :: litter_input_t
      real(wp) :: labile_grnd  = 0.0_wp      !< -> X(1) fast_grnd  (labile leaf/storage, above)
      real(wp) :: labile_soil  = 0.0_wp      !< -> X(2) fast_soil  (labile fineroot/storage, below)
      real(wp) :: struct_grnd  = 0.0_wp      !< -> X(3) struct_grnd (structural leaf + CWD, above)
      real(wp) :: struct_soil  = 0.0_wp      !< -> X(4) struct_soil (structural fineroot + belowground CWD)
      real(wp) :: lignin_grnd  = 0.0_wp      !< lignin flux to struct_grnd (sets f_lignin of incoming litter)
      real(wp) :: lignin_soil  = 0.0_wp      !< lignin flux to struct_soil
      ! optional N twin (driver-split by tissue C:N), present only when n_cycle_on:
      real(wp) :: n_labile_grnd = 0.0_wp, n_labile_soil = 0.0_wp
      real(wp) :: n_struct_grnd = 0.0_wp, n_struct_soil = 0.0_wp
   end type litter_input_t

   !----- Daily carbon-mass conservation guard (the fast/slow contract). --------------------------!
   type :: soilc_audit_t
      real(wp) :: litter_in    = 0.0_wp   !< [kgC/m2/day] sum(u)
      real(wp) :: rh_out       = 0.0_wp   !< [kgC/m2/day] Rh reported by soil_carbon_step (= litter_in - dC_pool)
      real(wp) :: rh_fast_accum= 0.0_wp   !< [kgC/m2/day] fast loop's accumulated today_rh (the CAS-fed flux)
      real(wp) :: dC_pool      = 0.0_wp   !< [kgC/m2/day] net pool change
      real(wp) :: resid        = 0.0_wp   !< [kgC/m2/day] dC_pool - (litter_in - rh_out); ~0 by construction
      real(wp) :: rh_seam_gap  = 0.0_wp   !< [kgC/m2/day] rh_out - rh_fast_accum; fast/slow reconciliation check (~0)
      real(wp) :: lignin_resid = 0.0_wp   !< [kgC/m2/day] max_s |dL_s - (lignin_in_s - d_s*L_s)|; passive-tracer check
   end type soilc_audit_t

   !----- Traceability diagnostics (pure post-processing off the assembled matrices). -------------!
   type :: soilc_diag_t
      real(wp) :: x_c(n_soil_pool) = 0.0_wp   !< [kgC/m2] storage CAPACITY (equilibrium under current climate+input)
      real(wp) :: x_p(n_soil_pool) = 0.0_wp   !< [kgC/m2] storage POTENTIAL = x_c - X (remaining sink; <0 = source)
      real(wp) :: tau_e            = 0.0_wp   !< [yr]     ecosystem soil-carbon residence time = sum(X_c)/sum(u)
      real(wp) :: rh               = 0.0_wp   !< [kgC/m2/day] instantaneous Rh = -1^T*A*xi*K*X
   end type soilc_diag_t

end module meds_biogeochem_types
