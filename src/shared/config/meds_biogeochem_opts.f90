!==========================================================================================!
! meds_biogeochem_opts -- the [soil_carbon] run-configuration bundle for the SLOW soil-carbon !
! matrix (decomposition selectors + rate/respired-fraction/environmental-scalar parameters),    !
! plus its selector codes.                                                                       !
!                                                                                          !
! Placed in src/shared (NOT src/biogeochemistry) so meds_config -- the DAG ROOT -- can carry     !
! decomp_opts_t as a plain component with NO backward `shared -> biogeochemistry` edge (the       !
! same rule meds_forcing_config/meds_biophysics_opts already follow for forcing_config_t/          !
! soil_opts_t etc.). Pure DATA + parameters; links meds_kinds only. meds_biogeochem_types           !
! re-exports every name below so `use meds_biogeochem_types, only : decomp_opts_t` (every        !
! existing kernel/test call site) compiles unchanged.                                             !
!==========================================================================================!
module meds_biogeochem_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: DECOMP_STEP_EULER, DECOMP_STEP_EXPM
   public :: DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   public :: decomp_opts_t

   !----- Daily-step solver selector (decomp_opts_t%step_solver). --------------------------------!
   integer(ik), parameter :: DECOMP_STEP_EULER = 1_ik  !< ED2-faithful forward Euler (the coupled daily step)
   integer(ik), parameter :: DECOMP_STEP_EXPM  = 2_ik  !< exact matrix exponential (LARGE accelerated/spin-up steps)

   !----- Decomposition SCHEME family (decomp_opts_t%decomp_scheme). The full enum is {0..5}; the  !
   !      two that matter for the MVP are named: 0-4 collapse to 3 active pools (microbial/passive    !
   !      inert, K=0), 5 is the full 5-active CENTURY (Bolker 1998 / Koven 2013).                     !
   integer(ik), parameter :: DECOMP_SCHEME_ED2      = 0_ik  !< 3-active, ED2 capped-exp temperature
   integer(ik), parameter :: DECOMP_SCHEME_CENTURY5 = 5_ik  !< 5-active, Collatz Q10 temperature

   !==========================================================================================!
   !  Pre-extracted decomposition selectors + parameters for the SLOW soil-carbon matrix (the      !
   !  decomp_coms.f90 / ed_params.f90 provenance; verified against ED2 at config load, §6). Baseline  !
   !  decay rates are carried in [1/yr] -- the kernel converts to [1/day] at the step seam (yr_day).    !
   !  In-type defaults so the standalone kernels/tests compile pre-B2, and so the DEFAULTED [soil_carbon] !
   !  TOML reads (meds_config_io) fall back to these ED2-verified values when a key is absent -- this   !
   !  block is NOT required-when-on (unlike PFT traits): soil_carbon_on gates whether the FEATURE runs,  !
   !  not whether these params are supplied. Env-scalar params MIRROR co2_opts_t so the daily xi and the !
   !  fast Rh use matched chemistry (the fast/slow reconciliation, wired at B2).                          !
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

end module meds_biogeochem_opts
