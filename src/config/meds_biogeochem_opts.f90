! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_opts -- the [soil_carbon] block for the slow soil-carbon matrix          !
! (docs/science/soil_carbon.md): the decomposition selectors and their codes, the decay    !
! rates, the respired fractions and the environmental-scalar parameters.                   !
!                                                                                          !
! A config leaf, so meds_config can carry decomp_opts_t as a plain component with no edge  !
! back to the soil-carbon kernels, which import it from here. Pure data; links meds_kinds  !
! only.                                                                                    !
!==========================================================================================!
module meds_biogeochem_opts
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: DECOMP_STEP_EULER, DECOMP_STEP_EXPM
   public :: DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   public :: decomp_opts_t

   !----- Daily-step solver (decomp_opts_t%step_solver; docs/science/soil_carbon.md §6). -------!
   integer(ik), parameter :: DECOMP_STEP_EULER = 1_ik  !< ED2-faithful forward Euler (the coupled daily step)
   integer(ik), parameter :: DECOMP_STEP_EXPM  = 2_ik  !< exact matrix exponential (LARGE accelerated/spin-up steps)

   !----- Decomposition scheme (decomp_opts_t%decomp_scheme), in ED2's numbering; the TOML key     !
   !      sets one of these two. ED2's schemes 0-4 share the 3-active topology (microbial and      !
   !      passive inert, K = 0); 5 is the full 5-active CENTURY (Bolker 1998 / Koven 2013);        !
   !      docs/science/soil_carbon.md §2. ---------------------------------------------------------!
   integer(ik), parameter :: DECOMP_SCHEME_ED2      = 0_ik  !< 3-active, ED2 capped-exp temperature
   integer(ik), parameter :: DECOMP_SCHEME_CENTURY5 = 5_ik  !< 5-active, Collatz Q10 temperature

   !==========================================================================================!
   !  The decomposition selectors and parameters of the slow soil-carbon matrix, with ED2's   !
   !  values (decomp_coms.f90, ed_params.f90) as in-type defaults. Every [soil_carbon] key is !
   !  optional: meds_config_io falls back to these defaults when a key is absent, so the      !
   !  block may be omitted (soil_carbon_on switches the feature, not whether these are        !
   !  supplied). Decay rates are carried in [1/yr]; the kernel converts them to [1/day]       !
   !  (yr_day). The one set of environmental-scalar parameters serves both the sub-daily Rh   !
   !  and the daily step, which is what lets the two reconcile (docs/science/soil_carbon.md   !
   !  §4, §7).                                                                                !
   !==========================================================================================!
   type :: decomp_opts_t
      integer(ik) :: decomp_scheme = DECOMP_SCHEME_ED2    !< DECOMP_SCHEME_ED2 (3-active) or DECOMP_SCHEME_CENTURY5 (5-active)
      integer(ik) :: step_solver   = DECOMP_STEP_EULER    !< {EULER, EXPM}
      logical     :: n_cycle_on    = .false.              !< nitrogen switch: read, but no kernel reads it (f_decomp = 1)
      !----- Baseline decay rates [1/yr] (ED2 ed_params.f90 decay_rate_*). ------------------------!
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
      !----- Above-ground shares of the fast and structural litter (ED2 agf_fsc/agf_stsc): read,    !
      !      but used by no kernel. necromass_to_litter splits each cohort's litter with its PFT's  !
      !      f_labile_leaf, f_labile_stem and aboveground_frac (docs/science/soil_carbon.md §5). ---!
      real(wp)    :: agf_fast = 0.5_wp, agf_struct = 0.7_wp
      !----- Environmental scalar xi (docs/science/soil_carbon.md eqs. (3)-(4)). ------------------!
      real(wp)    :: rh_q10               = 1.5_wp
      real(wp)    :: rh_t_ref             = 288.15_wp     !< [K] 15 C Q10 reference (scheme 5)
      real(wp)    :: resp_temp_increase   = 0.0757_wp     !< [1/K] ED2 scheme-0 slope
      real(wp)    :: resp_temp_ref        = 318.15_wp     !< [K] 45 C saturation
      real(wp)    :: resp_opt_water       = 0.8938_wp
      real(wp)    :: resp_water_below_opt = 5.0786_wp
      real(wp)    :: resp_water_above_opt = 4.5139_wp
      !----- N stoichiometry: read, but no kernel reads it (docs/science/soil_carbon.md §10). ------!
      real(wp)    :: c2n_structural = 150.0_wp, c2n_slow = 10.0_wp, c2n_fast = 30.0_wp
      real(wp)    :: n_immobil_supply_scale = 0.0_wp      !< [1/yr] 40.0 in ED2
   end type decomp_opts_t

end module meds_biogeochem_opts
