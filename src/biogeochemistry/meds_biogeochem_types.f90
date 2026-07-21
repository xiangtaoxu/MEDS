!==========================================================================================!
! meds_biogeochem_types -- shared derived types + selector codes for the biogeochemistry     !
! domain: the SLOW soil-carbon / nutrient cycle of the ecosystem column. Pure DATA +           !
! parameters, no methods -- the biogeochem analogue of meds_biophysics_types. Links              !
! src/shared (meds_kinds) ONLY.                                                                   !
!                                                                                          !
! (The FAST canopy-air-space CO2 exchange -- a sub-daily biophysical diffusion/venting process --  !
! moved to meds_biophysics_types + meds_cas_biophysics under src/biophysics; this module now holds      !
! only the SLOW soil-carbon pools.)                                                                 !
!                                                                                          !
! SLOW soil carbon (design MEDS_BIOGEOCHEMISTRY_DESIGN.md): the CENTURY-family multi-pool soil-    !
! carbon state advanced daily by the carbon matrix ODE dX/dt = B*I + A*xi*K*X. The one-pool MVP     !
! `soil_carbon_t` is expanded HERE (P0) to the 7-pool vector + lignin sub-state + optional N,        !
! ADDITIVELY: field `fast_soil_carbon` KEEPS its name/index (2) so meds_cas_biophysics -- which reads     !
! it as a BARE SCALAR in the fast loop -- compiles unchanged. Slow types: decomp_opts_t /              !
! litter_input_t / soilc_audit_t / soilc_diag_t, and the fixed pool count `n_soil_pool` + `IP_*`.       !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
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

   !----- Daily-step solver selector (decomp_opts_t%step_solver). --------------------------------!
   integer(ik), parameter :: DECOMP_STEP_EULER = 1_ik  !< ED2-faithful forward Euler (the coupled daily step)
   integer(ik), parameter :: DECOMP_STEP_EXPM  = 2_ik  !< exact matrix exponential (LARGE accelerated/spin-up steps)

   !----- Decomposition SCHEME family (decomp_opts_t%decomp_scheme). The full enum is {0..5}; the  !
   !      two that matter for the MVP are named: 0-4 collapse to 3 active pools (microbial/passive    !
   !      inert, K=0), 5 is the full 5-active CENTURY (Bolker 1998 / Koven 2013).                     !
   integer(ik), parameter :: DECOMP_SCHEME_ED2      = 0_ik  !< 3-active, ED2 capped-exp temperature
   integer(ik), parameter :: DECOMP_SCHEME_CENTURY5 = 5_ik  !< 5-active, Collatz Q10 temperature

   !==========================================================================================!
   !  Slow, stateful per-patch soil-carbon pools (written DAILY, read-only in the fast loop).       !
   !  Expanded from the P0 one-pool shape; fast_soil_carbon KEEPS its name/role (index 2) so           !
   !  meds_cas_biophysics compiles unchanged. The ACTIVE pool count is set by decomp_scheme (default        !
   !  3-active: idx 5,7 inert). Lignin is a passive tracer of the structural pools (0 <= L <= C),         !
   !  OUTSIDE the n_soil_pool carbon mass vector (its carbon is already counted in struct_*_carbon).       !
   !==========================================================================================!
   type :: soil_carbon_t
      ! carbon pools [kgC/m2] -- the matrix state vector X, ordered litter -> SOM -> passive
      real(wp) :: fast_grnd_carbon    = 0.0_wp   !< X(1) metabolic litter, above-ground (flammable)
      real(wp) :: fast_soil_carbon    = 0.0_wp   !< X(2) metabolic litter, below-ground  (the P0 pool)
      real(wp) :: struct_grnd_carbon  = 0.0_wp   !< X(3) structural litter + CWD, above
      real(wp) :: struct_soil_carbon  = 0.0_wp   !< X(4) structural litter + CWD, below
      real(wp) :: microbial_carbon    = 0.0_wp   !< X(5) microbial SOM        (scheme-5 only)
      real(wp) :: slow_carbon         = 0.0_wp   !< X(6) slow / humified SOM
      real(wp) :: passive_carbon      = 0.0_wp   !< X(7) passive SOM          (scheme-5 only)
      ! lignin sub-state of the structural pools [kgC/m2] (fraction f_lignin = L/C brakes decomposition)
      real(wp) :: struct_grnd_lignin  = 0.0_wp
      real(wp) :: struct_soil_lignin  = 0.0_wp
      ! optional nitrogen twin [kgN/m2] -- present only when opts%n_cycle_on (P1/P2); C-only default
      real(wp) :: fast_grnd_n = 0.0_wp, fast_soil_n = 0.0_wp
      real(wp) :: struct_grnd_n = 0.0_wp, struct_soil_n = 0.0_wp
      real(wp) :: mineralized_n = 0.0_wp
   end type soil_carbon_t

   !==========================================================================================!
   !  Pre-extracted decomposition selectors + parameters for the SLOW soil-carbon matrix (the      !
   !  decomp_coms.f90 / ed_params.f90 provenance; verified against ED2 at config load, §6). Baseline  !
   !  decay rates are carried in [1/yr] -- the kernel converts to [1/day] at the step seam (yr_day).    !
   !  In-type defaults so the standalone kernels/tests compile pre-P3. Env-scalar params MIRROR         !
   !  co2_opts_t so the daily xi and the fast Rh use matched chemistry (the fast/slow reconciliation).   !
   !==========================================================================================!
   type :: decomp_opts_t
      integer(ik) :: decomp_scheme = DECOMP_SCHEME_ED2    !< {0..5}; 0-4 = 3-active, 5 = full 5-pool CENTURY
      integer(ik) :: step_solver   = DECOMP_STEP_EULER    !< {EULER, EXPM}
      logical     :: n_cycle_on    = .false.              !< optional nitrogen (P1/P2); C-only default (f_decomp = 1)
      !----- baseline decay rates [1/yr] (verified vs ED2 ed_params.f90 decay_rate_* at load). -----!
      real(wp)    :: k_fast   = 11.0_wp    !< decay_rate_fsc  (scheme 0-4); 12.0 (scheme 5)
      real(wp)    :: k_struct = 4.5_wp     !< decay_rate_stsc (scheme 0-4); 1.5  (scheme 5)
      real(wp)    :: k_micr   = 0.0_wp     !< decay_rate_msc  = 0 (inert) schemes 0-4; 6.0 (scheme 5)
      real(wp)    :: k_slow   = 0.2_wp     !< decay_rate_ssc  scheme-dependent (100.2 / 0.2 branches); 0.15 (scheme 5)
      real(wp)    :: k_pass   = 0.0_wp     !< decay_rate_psc  = 0 (inert) schemes 0-4; 0.012 (scheme 5)
      !----- respired fractions [-] (er = 1 - transferred; = decomp_coms r_*). ---------------------!
      real(wp)    :: er_fast         = 1.0_wp   !< r_fsc   ; 0.55 (scheme 5)
      real(wp)    :: er_struct_nonlig= 0.3_wp   !< r_stsc_o; 0.50 (scheme 5)
      real(wp)    :: er_struct_lig   = 0.3_wp   !< r_stsc_l; 0.20 (scheme 5)
      real(wp)    :: er_slow         = 1.0_wp   !< r_ssc   ; 0.55 (scheme 5)
      real(wp)    :: er_pass         = 1.0_wp   !< r_psc   ; 0.55 (scheme 5)
      real(wp)    :: er_micr_int     = 0.0_wp   !< r_msc_int; 0.60 (scheme 5)  -> er_micr = int + slp*xsand
      real(wp)    :: er_micr_slp     = 0.0_wp   !< r_msc_slp; 0.17 (scheme 5)
      real(wp)    :: e_lignin        = 3.0_wp   !< lignin brake exponent  L = exp(-e_lignin*f_lignin)
      !----- CENTURY texture controls (scheme 5): sand -> er_micr, clay -> transfers to passive. ---!
      real(wp)    :: xsand = 0.35_wp, xclay = 0.25_wp
      real(wp)    :: fx_micr_pass_int = 0.003_wp, fx_micr_pass_slp = 0.032_wp   !< ex_micr->pass = int + slp*xclay
      real(wp)    :: fx_slow_pass_int = 0.003_wp, fx_slow_pass_slp = 0.009_wp   !< ex_slow->pass = int + slp*xclay
      !----- above/below split FALLBACKS (ED2 agf_fsc/agf_stsc); the driver splits per-COHORT with   !
      !      per-PFT agf/f_labile at P3, so these are NOT used by the kernels. -----------------------!
      real(wp)    :: agf_fast = 0.5_wp, agf_struct = 0.7_wp
      !----- environmental-scalar params (SHARED semantics with co2_opts_t; matched fast/slow). -----!
      real(wp)    :: rh_q10               = 1.5_wp
      real(wp)    :: rh_t_ref             = 288.15_wp     !< [K] 15 C Q10 reference (scheme 5)
      real(wp)    :: resp_temp_increase   = 0.0757_wp     !< [1/K] ED2 scheme-0 slope
      real(wp)    :: resp_temp_ref        = 318.15_wp     !< [K] 45 C saturation
      real(wp)    :: resp_opt_water       = 0.8938_wp
      real(wp)    :: resp_water_below_opt = 5.0786_wp
      real(wp)    :: resp_water_above_opt = 4.5139_wp
      !----- N stoichiometry (used only when n_cycle_on). ------------------------------------------!
      real(wp)    :: c2n_structural = 150.0_wp, c2n_slow = 10.0_wp, c2n_fast = 30.0_wp
      real(wp)    :: n_immobil_supply_scale = 0.0_wp      !< [1/yr] 40.0 in ED2 (converted by yr_day)
   end type decomp_opts_t

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
