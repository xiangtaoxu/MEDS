! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_diag_vars -- the DIAGNOSTIC accumulators, per cohort, per patch and per polygon.    !
!                                                                                          !
! Stage [2] of the diagnostic wall, CAPTURE (docs/science/diagnostics.md §1). The problem  !
! it solves: the fast loop computes ~40 per-cohort and per-patch quantities every dt_fast  !
! -- leaf assimilation, stomatal conductance, tissue water potentials, sapflow, absorbed   !
! radiation, turbulent fluxes, runoff, budget residuals -- and sub-daily resolution exists !
! only inside the fast loop's sub-step, so without a home they would be recomputed and     !
! discarded ~48 times a day. These blocks are where they land instead.                     !
!                                                                                          !
! SIX DESIGN DECISIONS:                                                                    !
!                                                                                          !
!  1. dt-WEIGHTED SUMS, normalized on read. Every field accumulates x*dt and `w`           !
!     accumulates dt; the reader divides. That is correct under the adaptive integrator    !
!     (unequal sub-steps) and under AGG_TMEAN chaining for free, where a plain running     !
!     mean would not be.                                                                   !
!                                                                                          !
!  2. RESET ONCE PER SLOW STEP -- the same lifecycle as gpp_accum / et_accum / xi_accum.   !
!     No second reset cadence is invented.                                                 !
!                                                                                          !
!  3. TRANSIENT: never written to the restart file. A checkpoint is prognostic state at an !
!     instant; a time-averaged diagnostic is not state.                                    !
!                                                                                          !
!  4. THEY STILL RIDE THE LOCKSTEP. The reset is per slow step, but the slow step re-sorts !
!     the cohorts after the fast loop has filled these and before the output tick reads    !
!     them, and the calendar restructuring (fuse / split / cull / recruit / disturb)       !
!     permutes the axis at a boundary. So slot i must keep meaning cohort i across every   !
!     permutation. That obligation is REAL and is the main cost of this design.            !
!                                                                                          !
!     It is discharged by STORAGE LAYOUT rather than by discipline: the fields are rows of !
!     ONE 2-D array v(N_CDIAG, cap), so every lockstep operation is a single whole-array   !
!     statement that cannot omit a field. Adding a diagnostic is a new index parameter     !
!     plus a fill line -- ZERO edits to the reorder machinery. The alternative, one        !
!     allocatable per field, would need a new line per field in each of alloc_cohort_block /     !
!     expand_cohort_block / move_cohort_block / reorder_cohort_block / copy_cohort --               !
!     precisely the "forgot to reallocate" class CLAUDE.md names as the reason the         !
!     centralized reorder exists.                                                          !
!                                                                                          !
!  5. ALLOCATED ONLY ON DEMAND. A run that reports no per-cohort ecophysiology pays        !
!     nothing: main sets `active` from the registry, and every entry point is a no-op      !
!     while it is .false.                                                                  !
!                                                                                          !
!  6. FUSION BLENDS BY THE FIELD'S OWN KIND. Each field declares the weight that carries   !
!     it: leaf area (INTENSIVE -- a conductance, a potential, a per-leaf-area flux), plant !
!     number (EXTENSIVE -- a per-plant flux), the tissue's heat capacity (a leaf or wood   !
!     temperature, as its prognostic twin fuses), or none (GROUND-REFERENCED: summed,      !
!     already per m2 of patch). Blending a temperature by nplant, or an extensive flux by  !
!     leaf area, is the classic silent fusion bug; the kind is DATA here, so the blend     !
!     cannot disagree with the quantity.                                                   !
!==========================================================================================!
module meds_diag_vars
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_numerics,  only : pair_mean
   implicit none
   private

   public :: cohort_diag_block, patch_diag_block, polygon_diag_block
   public :: polygon_diag_reset, polygon_diag_value
   public :: cohort_diag_alloc, cohort_diag_free, cohort_diag_grow, cohort_diag_reset
   public :: cohort_diag_reorder, cohort_diag_copy_slot, cohort_diag_fuse, cohort_diag_clear_slot
   public :: patch_diag_alloc, patch_diag_free, patch_diag_grow, patch_diag_reset
   public :: patch_diag_reorder, patch_diag_copy_slot, patch_diag_blend, patch_diag_clear_slot
   public :: patch_diag_inherit
   public :: cohort_diag_value, patch_diag_value, patch_diag_slow_row

   !----- Fusion kinds (decision 6 above). ---------------------------------------------------!
   integer(ik), parameter, public :: FK_INTENSIVE = 1_ik  !< leaf-area-weighted mean (gs, psi, per-leaf-area fluxes)
   integer(ik), parameter, public :: FK_EXTENSIVE = 2_ik  !< nplant-weighted mean (per-plant fluxes)
   integer(ik), parameter, public :: FK_GROUND    = 3_ik  !< summed (already per m2 of patch ground)
   integer(ik), parameter, public :: FK_LEAF_HEAT = 4_ik  !< leaf-heat-capacity-weighted mean (leaf temperature)
   integer(ik), parameter, public :: FK_WOOD_HEAT = 5_ik  !< wood-heat-capacity-weighted mean (wood temperature)

   !==========================================================================================!
   !  PER-COHORT diagnostic fields. Order is free; the parameters are the only thing that binds  !
   !  a row to a meaning.                                                                        !
   !==========================================================================================!
   !----- 1. Carbon -----------------------------------------------------------------------------!
   integer(ik), parameter, public :: CD_ANET         =  1_ik  !< [umol CO2/m2 leaf/s] net assimilation
   integer(ik), parameter, public :: CD_AGROSS       =  2_ik  !< [umol CO2/m2 leaf/s] gross assimilation
   integer(ik), parameter, public :: CD_RD           =  3_ik  !< [umol/m2 leaf/s] leaf dark respiration
   integer(ik), parameter, public :: CD_GPP_RATE     =  4_ik  !< [umol CO2/plant/s]
   !----- 2. Stomata and leaf exchange ----------------------------------------------------------!
   integer(ik), parameter, public :: CD_GSW          =  5_ik  !< [mol H2O/m2 leaf/s]  stomatal conductance
   integer(ik), parameter, public :: CD_GBW          =  6_ik  !< [m/s] leaf boundary-layer conductance
   integer(ik), parameter, public :: CD_CI           =  7_ik  !< [umol/mol] intercellular CO2
   integer(ik), parameter, public :: CD_CS           =  8_ik  !< [umol/mol] leaf-surface CO2
   integer(ik), parameter, public :: CD_TRANSP       =  9_ik  !< [mol H2O/m2 leaf/s] leaf transpiration
   integer(ik), parameter, public :: CD_BETA_STOM    = 10_ik  !< [-] stomatal water-stress limb
   integer(ik), parameter, public :: CD_BETA_NONSTOM = 11_ik  !< [-] non-stomatal (capacity) limb
   !----- 3. Plant water ------------------------------------------------------------------------!
   integer(ik), parameter, public :: CD_PSI_LEAF     = 12_ik  !< [MPa]
   integer(ik), parameter, public :: CD_PSI_WOOD     = 13_ik  !< [MPa]
   integer(ik), parameter, public :: CD_PLC          = 14_ik  !< [-] percent loss of conductance
   integer(ik), parameter, public :: CD_SAPFLOW      = 15_ik  !< [kg/plant/s] wood -> leaf
   integer(ik), parameter, public :: CD_ROOT_UPTAKE  = 16_ik  !< [kg/plant/s] soil -> root
   integer(ik), parameter, public :: CD_LEAF_WATER   = 17_ik  !< [kg/plant] internal leaf water
   integer(ik), parameter, public :: CD_WOOD_WATER   = 18_ik  !< [kg/plant] internal wood water
   !----- 4. Energy and radiation ---------------------------------------------------------------!
   integer(ik), parameter, public :: CD_LEAF_TEMP    = 19_ik  !< [K]
   integer(ik), parameter, public :: CD_WOOD_TEMP    = 20_ik  !< [K]
   integer(ik), parameter, public :: CD_ABS_PAR      = 21_ik  !< [W/m2 ground] absorbed PAR
   integer(ik), parameter, public :: CD_ABS_SW       = 22_ik  !< [W/m2 ground] absorbed shortwave
   integer(ik), parameter, public :: CD_ABS_LW       = 23_ik  !< [W/m2 ground] net longwave
   !----- 5. Canopy environment -----------------------------------------------------------------!
   integer(ik), parameter, public :: CD_LEAF_VPD     = 24_ik  !< [Pa] leaf-to-air VPD
   integer(ik), parameter, public :: CD_WIND         = 25_ik  !< [m/s] in-canopy wind
   !----- 6. The leaf solve's other inputs ------------------------------------------------------!
   !----- The rest of the leaf solve's own inputs, as it saw them, so the solve can be repeated     !
   !      outside the model (the calibration's canopy of leaf kernels): with CD_LEAF_TEMP,           !
   !      CD_LEAF_VPD and CD_PSI_LEAF they are every driver of leaf_gas_exchange_batch. ------------!
   integer(ik), parameter, public :: CD_LEAF_PAR     = 26_ik  !< [umol photon/m2 leaf/s] incident-equivalent PAR per leaf area
   integer(ik), parameter, public :: CD_GB_MOL       = 27_ik  !< [mol H2O/m2 leaf/s] boundary-layer conductance
   integer(ik), parameter, public :: CD_PSI_PREDAWN  = 28_ik  !< [MPa] the stomatal limb's predawn leaf potential
   integer(ik), parameter, public :: CD_CA           = 29_ik  !< [umol/mol] canopy-air CO2 the leaf sees
   integer(ik), parameter, public :: CD_PRESSURE     = 30_ik  !< [Pa] air pressure the leaf sees
   integer(ik), parameter, public :: N_CDIAG         = 30_ik

   !----- Fusion kind per cohort field. The leaf and wood temperatures fuse on the heat capacity  !
   !      of their tissue, as their prognostic twins do (fuse_cohort_fast_state), so a fused       !
   !      cohort's mean temperature over the output window is the one its tissue energy implies.   !
   !      Conductances, potentials and the leaf-area-normalized fluxes are INTENSIVE (leaf-area-   !
   !      weighted). The per-PLANT quantities (GPP rate, sapflow, uptake, tissue water mass) are   !
   !      EXTENSIVE (nplant-weighted), matching how their prognostic twins fuse. The absorbed-     !
   !      radiation terms are per m2 of GROUND, so two cohorts' contributions to the same patch    !
   !      simply add -- the same convention leaf_surf_water follows. -----------------------------!
   integer(ik), parameter, public :: CDIAG_FUSE(N_CDIAG) = [                                              &
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_EXTENSIVE,                                           &   ! 1. Carbon
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, &   ! 2. Stomata and leaf exchange
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, &   ! 3. Plant water
        FK_LEAF_HEAT, FK_WOOD_HEAT, FK_GROUND, FK_GROUND, FK_GROUND,                                      &   ! 4. Energy and radiation
        FK_INTENSIVE, FK_INTENSIVE,                                                                       &   ! 5. Canopy environment
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE                              ]   ! 6. The leaf solve's other inputs

   !==========================================================================================!
   !  PER-COHORT SLOW-LOOP diagnostic fields -- a SECOND block of the same shape, written once per  !
   !  slow step by the vegetation-dynamics driver.                                                  !
   !                                                                                          !
   !  WHY A SEPARATE BLOCK AND NOT MORE ROWS ABOVE. The two have different sample semantics: the     !
   !  fast rows accumulate ~48 sub-step samples per slow step, the slow rows get exactly one. One     !
   !  shared dt weight cannot normalize both.                                                         !
   !                                                                                          !
   !  WHY NOT READ site%deriv DIRECTLY. `cohort_deriv_block` is TRANSIENT and deliberately NOT        !
   !  lockstep-reordered -- which is fine for its own consumer, since update_cohort_states applies it  !
   !  immediately. But the step re-sorts the cohorts after that, so at the output tick `deriv(i)` and  !
   !  `cohort(i)` need not be the same plant, and an error of that kind shows only when the ordering   !
   !  changes, as thread count can make it. Storing the tendencies in a block that DOES ride the       !
   !  lockstep removes the failure mode rather than timing around it.                                  !
   !==========================================================================================!
   !----- 1. Growth and mortality ---------------------------------------------------------------!
   integer(ik), parameter, public :: CS_DDBH_DT      =  1_ik  !< [cm/yr]         diameter growth
   integer(ik), parameter, public :: CS_DAGB_DT      =  2_ik  !< [kgC/plant/yr]  AGB growth
   integer(ik), parameter, public :: CS_MORT_RATE    =  3_ik  !< [1/yr]          mortality hazard (positive)
   !----- 2. Carbon allocation ------------------------------------------------------------------!
   integer(ik), parameter, public :: CS_NPP_LEAF     =  4_ik  !< [kgC/plant/yr]  NPP to leaf
   integer(ik), parameter, public :: CS_NPP_FINEROOT =  5_ik  !< [kgC/plant/yr]  NPP to fine root
   integer(ik), parameter, public :: CS_NPP_WOOD     =  6_ik  !< [kgC/plant/yr]  NPP to wood
   integer(ik), parameter, public :: CS_NPP_STORAGE  =  7_ik  !< [kgC/plant/yr]  NPP to non-structural storage
   integer(ik), parameter, public :: CS_NPP_REPRO    =  8_ik  !< [kgC/plant/yr]  NPP to reproduction
   !----- 3. Carbon losses ----------------------------------------------------------------------!
   integer(ik), parameter, public :: CS_GROWTH_RESP  =  9_ik  !< [kgC/plant/yr]  growth respiration
   integer(ik), parameter, public :: CS_STORAGE_RESP = 10_ik  !< [kgC/plant/yr]  storage maintenance respiration
   integer(ik), parameter, public :: CS_EXUDATE      = 11_ik  !< [kgC/plant/yr]  root exudate (carbon the growth sink could not use)
   integer(ik), parameter, public :: N_CSDIAG        = 11_ik

   !----- Fusion kinds for the slow rows. All of them are per-PLANT rates, so all are EXTENSIVE    !
   !      (nplant-weighted) -- the same convention gpp_accum and the maintenance-respiration        !
   !      accumulators already follow in fuse_2_cohorts.  ---------------------------------------!
   integer(ik), parameter, public :: CSDIAG_FUSE(N_CSDIAG) = [                                   &
        FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE,                    &
        FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE,                    &
        FK_EXTENSIVE ]

   !==========================================================================================!
   !  PER-PATCH diagnostic fields. The fast loop fills the fast rows once per (patch, sub-step)  !
   !  (patch_diag_row): the coarse tiers accumulate that row, and the FAST tier stages it as it  !
   !  is, so the two read one list.                                                             !
   !==========================================================================================!
   !----- 1. Forcing at the canopy-air top ------------------------------------------------------!
   !----- The patch's own forcing: the sample moved to its canopy-air top (docs/science/forcing.md  !
   !      §8). Only the wind and the air temperature differ between patches.                       !
   integer(ik), parameter, public :: PD_WIND_CAS_TOP      =  1_ik  !< [m/s] wind speed at the canopy-air top
   integer(ik), parameter, public :: PD_TAIR_CAS_TOP      =  2_ik  !< [K]   air temperature at the canopy-air top
   integer(ik), parameter, public :: PD_Z_CAS_TOP         =  3_ik  !< [m]   that top: the canopy-air depth the fast loop used
   !----- 2. Top-of-canopy radiation ------------------------------------------------------------!
   !----- TOP-OF-CANOPY radiative fluxes. Kept as FLUXES, not as a time-averaged albedo: a        !
   !      period-mean albedo is the mean of a RATIO, which is not the ratio of the means, and     !
   !      at night the shortwave ratio is 0/0. The albedo a reader wants is sum(up)/sum(down)     !
   !      over the period (docs/science/canopy_radiation_transfer.md eq. (15)), which these       !
   !      make computable -- and which is what a satellite product is. ---------------------------!
   integer(ik), parameter, public :: PD_SW_IN_VIS         =  4_ik  !< [W/m2] incident VIS at canopy top
   integer(ik), parameter, public :: PD_SW_IN_NIR         =  5_ik  !< [W/m2] incident NIR at canopy top
   integer(ik), parameter, public :: PD_SW_UP_VIS         =  6_ik  !< [W/m2] upwelling VIS leaving the canopy top
   integer(ik), parameter, public :: PD_SW_UP_NIR         =  7_ik  !< [W/m2] upwelling NIR leaving the canopy top
   integer(ik), parameter, public :: PD_SW_UP             =  8_ik  !< [W/m2] upwelling shortwave, both bands (VIS + NIR)
   integer(ik), parameter, public :: PD_LW_UP             =  9_ik  !< [W/m2] upwelling LW (surface emission included)
   !----- The SKIN temperature, in the land-model sense (CLM's TSKIN): the black-body             !
   !      temperature of the longwave leaving the canopy top, what an infrared thermometer or     !
   !      a satellite land-surface temperature sees (docs/science/diagnostics.md §4, "The         !
   !      skin temperature"). MEDS has no skin layer of its own -- snow-free, the ground          !
   !      surface is the top soil layer (PD_SOIL_TEMP_TOP) -- and over a forest this is           !
   !      mostly the canopy. It is a fourth root, so it is formed per sub-step, like the VPD. ----!
   integer(ik), parameter, public :: PD_SKIN_TEMP         = 10_ik  !< [K] (longwave up / sigma)^(1/4)
   !----- 3. Turbulence -------------------------------------------------------------------------!
   integer(ik), parameter, public :: PD_USTAR             = 11_ik  !< [m/s]
   integer(ik), parameter, public :: PD_GGNET             = 12_ik  !< [m/s] ground conductance
   integer(ik), parameter, public :: PD_ROUGH             = 13_ik  !< [m] roughness length
   integer(ik), parameter, public :: PD_DISPLACE          = 14_ik  !< [m] displacement height
   !----- 4. Surface energy balance -------------------------------------------------------------!
   integer(ik), parameter, public :: PD_RNET              = 15_ik  !< [W/m2] net all-wave radiation
   integer(ik), parameter, public :: PD_H                 = 16_ik  !< [W/m2] sensible heat (CAS -> atm)
   integer(ik), parameter, public :: PD_LE                = 17_ik  !< [W/m2] latent heat (CAS -> atm)
   integer(ik), parameter, public :: PD_SW_GROUND         = 18_ik  !< [W/m2] shortwave reaching the ground
   integer(ik), parameter, public :: PD_LW_GROUND         = 19_ik  !< [W/m2] net longwave at the ground
   !----- 5. Carbon and water fluxes ------------------------------------------------------------!
   integer(ik), parameter, public :: PD_GPP               = 20_ik  !< [umol CO2/m2/s] patch GPP rate
   integer(ik), parameter, public :: PD_NPP               = 21_ik  !< [umol CO2/m2/s] GPP net of maintenance respiration
   integer(ik), parameter, public :: PD_RECO              = 22_ik  !< [umol CO2/m2/s] ecosystem respiration, NEE + GPP
   integer(ik), parameter, public :: PD_NEE               = 23_ik  !< [umol CO2/m2/s] net ecosystem exchange
   integer(ik), parameter, public :: PD_TRANSP            = 24_ik  !< [kg/m2/s] canopy transpiration
   !----- 6. Canopy air and ground state --------------------------------------------------------!
   !----- Fast-varying CAS / surface state, accumulated dt-weighted rather than read at the      !
   !      output tick. A state read at the tick is ONE instantaneous sample per output           !
   !      window, taken at whatever local time the boundary falls on, so its bias is a           !
   !      function of the site's longitude -- measured at +1.4 K for the top soil layer at       !
   !      Ithaca. -------------------------------------------------------------------------------!
   integer(ik), parameter, public :: PD_CAS_TEMP          = 25_ik  !< [K]
   integer(ik), parameter, public :: PD_CAS_SHV           = 26_ik  !< [kg/kg]
   integer(ik), parameter, public :: PD_CAS_CO2           = 27_ik  !< [umol/mol]
   integer(ik), parameter, public :: PD_CAS_VPD           = 28_ik  !< [Pa]
   integer(ik), parameter, public :: PD_SOIL_TEMP_TOP     = 29_ik  !< [K] top soil layer (the ground surface when snow-free)
   integer(ik), parameter, public :: PD_W_SURFACE         = 30_ik  !< [kg/m2] ponded surface water
   !----- 7. Budget residuals -------------------------------------------------------------------!
   integer(ik), parameter, public :: PD_RESID_ENERGY      = 31_ik  !< [W/m2] whole-column energy residual
   integer(ik), parameter, public :: PD_RESID_WATER       = 32_ik  !< [kg/m2/s] whole-column water residual
   integer(ik), parameter, public :: N_PDIAG_FAST         = 32_ik  !< the fast loop's rows are 1..N_PDIAG_FAST; the slow driver's follow
   !----- 8. SLOW rows: litter ------------------------------------------------------------------!
   !----- SLOW-loop patch diagnostics. They share this block (and its weight) with the fast rows    !
   !      above because the reader normalizes by w either way, and the slow driver writes them as   !
   !      a dt-weighted contribution over the same slow step the fast rows span. -------------------!
   integer(ik), parameter, public :: PD_LITTER_LEAF       = 33_ik  !< [kgC/m2/yr] leaf litterfall (after resorption)
   integer(ik), parameter, public :: PD_LITTER_FINEROOT   = 34_ik  !< [kgC/m2/yr] fine-root litter (after resorption)
   !----- Litter by the tissue it came from, completing PD_LITTER_LEAF / PD_LITTER_FINEROOT:     !
   !      with root_exudate_site they sum to the carbon entering the soil pools                  !
   !      (docs/science/soil_carbon.md §5). -----------------------------------------------------!
   integer(ik), parameter, public :: PD_LITTER_WOOD       = 35_ik  !< [kgC/m2/yr] wood + storage of dead plants
   integer(ik), parameter, public :: PD_LITTER_REPRO      = 36_ik  !< [kgC/m2/yr] fruit and seed that made no recruit
   integer(ik), parameter, public :: PD_LITTER_STRUCT     = 37_ik  !< [kgC/m2/yr] structural litter + CWD
   !----- 9. SLOW rows: demography --------------------------------------------------------------!
   integer(ik), parameter, public :: PD_RECRUIT_NPLANT    = 38_ik  !< [plant/m2/yr] recruitment density flux
   integer(ik), parameter, public :: PD_DISTURB_AREA      = 39_ik  !< [1/yr] area fraction disturbed
   !----- MORTALITY CARBON BY PATHWAY: the three ways a MEDS plant can die, kept apart so a       !
   !      user can tell a stand thinning continuously from one being knocked over. The three      !
   !      are mutually exclusive and jointly exhaustive -- every kilogram of carbon that          !
   !      leaves the live pools as a whole-individual death passes through exactly one of         !
   !      them -- so they SUM to the total, which is the property test_mortality_pathways         !
   !      asserts.                                                                                !
   !                                                                                              !
   !      Written UNCONDITIONALLY, unlike the litter each pathway also generates: two of the      !
   !      three litter sites sit behind [soil_carbon].soil_carbon_on, and how much biomass        !
   !      died is a DEMOGRAPHIC question that does not stop being asked when the soil pools       !
   !      are switched off. Whole-individual death carries every pool (leaf + fine root +         !
   !      wood + storage), as the litter partitioning does. --------------------------------------!
   integer(ik), parameter, public :: PD_MORT_C_BACKGROUND = 40_ik  !< [kgC/m2/yr] continuous hazard mortality
   integer(ik), parameter, public :: PD_MORT_C_CULL       = 41_ik  !< [kgC/m2/yr] cohorts culled below the tracking floor
   integer(ik), parameter, public :: PD_MORT_C_DISTURB    = 42_ik  !< [kgC/m2/yr] canopy killed by patch disturbance
   integer(ik), parameter, public :: N_PDIAG              = 42_ik

   !==========================================================================================!
   !  PER-POLYGON diagnostic fields: quantities that are the same everywhere in the polygon,  !
   !  so they are accumulated ONCE, not per patch. The prefixes: CD_ cohort (fast loop), CS_  !
   !  cohort (slow loop), PD_ patch, PY_ polygon. The forcing echo (the `forcing` group,      !
   !  docs/science/diagnostics.md §5) is the forcing the fast loop used, after the reader's   !
   !  shortwave partition, rain/snow split and optional height and lapse corrections, dt-     !
   !  weighted over the sub-steps.                                                            !
   !==========================================================================================!
   !----- 1. Radiation --------------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_SW_IN       =  1_ik  !< [W/m2] incident shortwave at canopy top
   integer(ik), parameter, public :: PY_PAR_BEAM    =  2_ik  !< [W/m2] direct-beam PAR
   integer(ik), parameter, public :: PY_PAR_DIFFUSE =  3_ik  !< [W/m2] diffuse PAR
   integer(ik), parameter, public :: PY_NIR_BEAM    =  4_ik  !< [W/m2] direct-beam NIR
   integer(ik), parameter, public :: PY_NIR_DIFFUSE =  5_ik  !< [W/m2] diffuse NIR
   integer(ik), parameter, public :: PY_LWDOWN      =  6_ik  !< [W/m2] downward longwave
   integer(ik), parameter, public :: PY_COSZ        =  7_ik  !< [-] cosine of the solar zenith angle (0 with the sun down)
   !----- 2. Precipitation ----------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_PRECIP      =  8_ik  !< [kg/m2/s] total precipitation (rain + snow)
   integer(ik), parameter, public :: PY_RAINF       =  9_ik  !< [kg/m2/s] liquid precipitation
   integer(ik), parameter, public :: PY_SNOWFALL    = 10_ik  !< [kg/m2/s] frozen precipitation
   !----- 3. Air --------------------------------------------------------------------------------!
   integer(ik), parameter, public :: PY_TAIR        = 11_ik  !< [K] reference-height air temperature
   integer(ik), parameter, public :: PY_QAIR        = 12_ik  !< [kg/kg] specific humidity
   integer(ik), parameter, public :: PY_PSURF       = 13_ik  !< [Pa] surface pressure
   integer(ik), parameter, public :: PY_WIND        = 14_ik  !< [m/s] wind speed at the reference height
   integer(ik), parameter, public :: PY_CO2         = 15_ik  !< [umol/mol] free-atmosphere CO2
   integer(ik), parameter, public :: PY_RHO_AIR     = 16_ik  !< [kg/m3] air density
   integer(ik), parameter, public :: N_PYDIAG       = 16_ik

   !==========================================================================================!
   !  The blocks themselves. `v` is (field, slot): field-major so a lockstep permutation of the  !
   !  SLOT axis is one contiguous whole-array statement (decision 4). Extraction reads one field  !
   !  across slots, which is strided -- acceptable, because this is the once-per-step diagnostic   !
   !  path, not the inner integrator.                                                              !
   !==========================================================================================!
   type :: cohort_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      integer(ik) :: nfield = 0_ik           !< N_CDIAG (fast block) or N_CSDIAG (slow block)
      logical     :: active = .false.        !< .false. => every entry point is a no-op, nothing allocated
      real(wp), allocatable :: v(:,:)        !< (nfield, cap) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (cap)         running Sum(dt)
   end type cohort_diag_block

   type :: patch_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      logical     :: active = .false.
      !----- Whether the fast loop fills the block. A slow-only run fills only the slow operators'  !
      !      rows (patch_diag_slow_row); the rest are not simulated and read as missing. -----------!
      logical     :: fast_rows = .true.
      real(wp), allocatable :: v(:,:)        !< (N_PDIAG, cap) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (cap)          running Sum(dt)
   end type patch_diag_block

   !----- One slot: the polygon (the site). Fixed size, so nothing to allocate or blend. -------!
   type :: polygon_diag_block
      logical  :: active = .false.
      real(wp) :: v(N_PYDIAG) = 0.0_wp      !< running Sum(x*dt)
      real(wp) :: w = 0.0_wp                !< running Sum(dt)
   end type polygon_diag_block

contains

   !=======================================================================================!
   !  COHORT block lifecycle. Every routine no-ops while `active` is .false., so a run that   !
   !  reports no per-cohort ecophysiology never allocates and never pays.                     !
   !=======================================================================================!

   subroutine cohort_diag_alloc(d, cap, active, nfield)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: cap
      logical,                 intent(in)    :: active
      integer(ik), optional,   intent(in)    :: nfield   !< default N_CDIAG (the fast block)
      integer(ik) :: nf
      nf = N_CDIAG
      if (present(nfield)) nf = nfield
      d%active = active
      d%nfield = nf
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(nf, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%v = 0.0_wp
      d%w = 0.0_wp          ! the first window starts empty: nothing zeroes it later
      d%cap = max(cap, 1_ik)
      d%n = 0_ik
      call cohort_diag_reset(d)
   end subroutine cohort_diag_alloc

   subroutine cohort_diag_free(d)
      type(cohort_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik
      d%n = 0_ik
      d%active = .false.
   end subroutine cohort_diag_free

   !----- Grow to at least `need` slots, preserving the live prefix (the ensure_capacity twin). !
   subroutine cohort_diag_grow(d, need)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      !----- Keep every existing slot, not just the first n: a slow operator may have written rows    !
      !      ahead of the count (a boundary's restructuring before any fast loop has set n). ---------!
      m = d%cap
      allocate(v2(int(size(d%v,1),ik), need), w2(need))
      v2 = 0.0_wp
      w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m)
         w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v)
      call move_alloc(w2, d%w)
      d%cap = need
   end subroutine cohort_diag_grow

   !----- Zero the whole block for a fresh slow step. ---------------------------------------!
   subroutine cohort_diag_reset(d)
      type(cohort_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine cohort_diag_reset

   !----- Permute the live prefix in lockstep with the cohort SoA. ONE statement -- it cannot   !
   !      omit a field, which is the whole point of the 2-D layout (decision 4).                 !
   !                                                                                          !
   !      THE COUNT RIDES THE LOCKSTEP TOO. `m` is the post-operator cohort count, and         !
   !      reorder_cohort_block sets cohort%n = m right after calling this; a block that kept its     !
   !      old d%n would claim more live slots than the cohort array has. Everything            !
   !      downstream believes d%n: cohort_diag_value fills 1..d%n into a caller buffer         !
   !      that extract_variable sizes from cohort%n, so after any operator that SHRINKS        !
   !      the array -- a cull is the common one -- that write would run off the end: an        !
   !      out-of-bounds WRITE, silent in Release. Assigned BEFORE the m <= 0 bail, so          !
   !      emptying the block records the emptying. --------------------------------------------!
   subroutine cohort_diag_reorder(d, perm, m)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: perm(:), m
      if (.not. d%active) return
      d%n = max(m, 0_ik)
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine cohort_diag_reorder

   subroutine cohort_diag_copy_slot(d, dst, src)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine cohort_diag_copy_slot

   !----- Zero one slot (a fresh birth, or a reused stale cull slot). -----------------------!
   subroutine cohort_diag_clear_slot(d, i)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp
      d%w(i) = 0.0_wp
   end subroutine cohort_diag_clear_slot

   !----- Blend the donor into the survivor on cohort fusion, each field by its own kind.      !
   !                                                                                           !
   !      `li_*` (leaf area) and `np_*` (nplant) are the weights the caller,                   !
   !      fuse_2_cohorts, computes. The extensive rows share their weight with their           !
   !      prognostic twins; the intensive rows, temperatures included, are leaf-area           !
   !      means, while the prognostic temperatures fuse on heat capacity                       !
   !      (fuse_cohort_fast_state). The dt weight `w` is a property of the window, not of      !
   !      either cohort (see below). ----------------------------------------------------------!
   subroutine cohort_diag_fuse(d, recc, donc, li_r, li_d, np_r, np_d, hl_r, hl_d, hw_r, hw_d, kinds)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: recc, donc
      real(wp),                intent(in)    :: li_r, li_d   !< leaf-area weights (nplant*leaf_area)
      real(wp),                intent(in)    :: np_r, np_d   !< nplant weights
      real(wp),                intent(in)    :: hl_r, hl_d   !< [J/m2/K] leaf heat capacities, before the merge
      real(wp),                intent(in)    :: hw_r, hw_d   !< [J/m2/K] wood heat capacities, before the merge
      integer(ik),             intent(in)    :: kinds(:)     !< CDIAG_FUSE or CSDIAG_FUSE
      integer(ik) :: f
      if (.not. d%active) return
      do f = 1_ik, min(int(size(kinds), ik), int(size(d%v, 1), ik))
         select case (kinds(f))
         case (FK_INTENSIVE)
            d%v(f, recc) = pair_mean(d%v(f, recc), li_r, d%v(f, donc), li_d)
         case (FK_EXTENSIVE)
            d%v(f, recc) = pair_mean(d%v(f, recc), np_r, d%v(f, donc), np_d)
         case (FK_LEAF_HEAT)
            d%v(f, recc) = pair_mean(d%v(f, recc), hl_r, d%v(f, donc), hl_d)
         case (FK_WOOD_HEAT)
            d%v(f, recc) = pair_mean(d%v(f, recc), hw_r, d%v(f, donc), hw_d)
         case default   ! FK_GROUND -- already per m2 of the SAME patch ground, so they add
            d%v(f, recc) = d%v(f, recc) + d%v(f, donc)
         end select
      end do
      !----- The dt weight is a property of the WINDOW, not of either cohort, and both slots saw  !
      !      the same sub-steps -- so the survivor keeps its own rather than summing (summing would !
      !      halve every normalized mean).  -------------------------------------------------!
   end subroutine cohort_diag_fuse

   !----- Extract one field, normalized by its own weight. `n` is what the caller may          !
   !      READ, and it is bounded by the caller's array as well as by the block: a kernel      !
   !      that writes past a caller's buffer on a count mismatch should not be possible,       !
   !      and this costs nothing on the once-per-step diagnostic path. With the count          !
   !      riding the lockstep (cohort_diag_reorder) the clamp never bites; it is what          !
   !      turns the next instance of that class into a short read instead of memory            !
   !      corruption. -------------------------------------------------------------------------!
   pure subroutine cohort_diag_value(d, field, x, n, empty)
      type(cohort_diag_block), intent(in)  :: d
      integer(ik),             intent(in)  :: field
      real(wp),                intent(out) :: x(:)
      integer(ik),             intent(out) :: n
      real(wp),                intent(in)  :: empty   !< what a slot with no samples reads
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      !----- A cohort with no samples reads `empty`: one recruited in the slow step, after the      !
      !      fast loop, or a slot past the block's count. Reading 0 instead would make a            !
      !      recruit's leaf temperature 0 K and pull every canopy mean towards 0. ------------------!
      x = empty
      n = min(d%n, int(size(x), ik))
      do i = 1_ik, n
         if (d%w(i) > tiny_num) x(i) = d%v(field, i) / d%w(i)
      end do
   end subroutine cohort_diag_value

   !=======================================================================================!
   !  PATCH block lifecycle -- the same shape, one axis smaller.                             !
   !=======================================================================================!

   subroutine patch_diag_alloc(d, cap, active)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: cap
      logical,                intent(in)    :: active
      d%active = active
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(N_PDIAG, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%v = 0.0_wp
      d%w = 0.0_wp          ! as cohort_diag_alloc
      d%cap = max(cap, 1_ik)
      d%n = 0_ik
      call patch_diag_reset(d)
   end subroutine patch_diag_alloc

   subroutine patch_diag_free(d)
      type(patch_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik
      d%n = 0_ik
      d%active = .false.
   end subroutine patch_diag_free

   subroutine patch_diag_grow(d, need)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      m = d%cap                              ! every existing slot (see cohort_diag_grow)
      allocate(v2(N_PDIAG, need), w2(need))
      v2 = 0.0_wp
      w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m)
         w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v)
      call move_alloc(w2, d%w)
      d%cap = need
   end subroutine patch_diag_grow

   subroutine patch_diag_reset(d)
      type(patch_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine patch_diag_reset

   !----- The patch analogue, and the count rides it for the same reason. ------------------!
   subroutine patch_diag_reorder(d, perm, m)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: perm(:), m
      if (.not. d%active) return
      d%n = max(m, 0_ik)
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine patch_diag_reorder

   subroutine patch_diag_copy_slot(d, dst, src)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine patch_diag_copy_slot

   subroutine patch_diag_clear_slot(d, i)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp
      d%w(i) = 0.0_wp
   end subroutine patch_diag_clear_slot

   !----- Area-weighted blend on patch fusion. Every patch diagnostic is per m2 of ITS OWN     !
   !      ground, so all of them blend the same way -- there is no extensive/intensive split at  !
   !      this level, which is why there is no kind array here.                                  !
   subroutine patch_diag_blend(d, recp, donp, area_r, area_d)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: recp, donp
      real(wp),               intent(in)    :: area_r, area_d
      real(wp) :: wtot
      if (.not. d%active) return
      wtot = area_r + area_d
      if (wtot <= tiny_num) return
      d%v(:, recp) = (area_r * d%v(:, recp) + area_d * d%v(:, donp)) / wtot
   end subroutine patch_diag_blend

   !----- A patch carved out of others inherits their slots, sums and weight, blended by the area  !
   !      each donor gave. Every donor loses the same fraction, so blending by the donors' areas is  !
   !      exact and the stand total Sum(area * v) is unchanged by the carving; a cleared slot would  !
   !      lose the carved area's share. At a calendar boundary the donors hold only that boundary's !
   !      events (weight 0: the step's sums were read and reset); inside a step (the C API) they    !
   !      hold the step so far, and the site mean Sum(area * v/w) is unchanged as well. `n` is      !
   !      raised to cover the new slot, so a reader sees it before the next reorder.                !
   subroutine patch_diag_inherit(d, dst, donor_area)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: dst
      real(wp),               intent(in)    :: donor_area(:)   !< donors are slots 1..size(donor_area)
      real(wp)    :: atot, wd
      integer(ik) :: i
      if (.not. d%active) return
      if (dst < 1_ik .or. dst > d%cap) return
      d%v(:, dst) = 0.0_wp
      d%w(dst) = 0.0_wp
      d%n = max(d%n, dst)
      atot = sum(donor_area)
      if (atot <= tiny_num) return
      do i = 1_ik, size(donor_area, kind=ik)
         wd = donor_area(i) / atot
         d%v(:, dst) = d%v(:, dst) + wd * d%v(:, i)
         d%w(dst)    = d%w(dst)    + wd * d%w(i)
      end do
   end subroutine patch_diag_inherit

   !----- The rows the slow operators write (litterfall, recruitment, disturbance, mortality      !
   !      carbon): the ones a slow-only run fills. Every other row is the fast loop's. ------------!
   pure logical function patch_diag_slow_row(field) result(slow)
      integer(ik), intent(in) :: field
      !----- The slow driver's rows follow the fast loop's (see N_PDIAG_FAST). ------------------!
      slow = field > N_PDIAG_FAST
   end function patch_diag_slow_row

   pure subroutine patch_diag_value(d, field, x, n)
      type(patch_diag_block), intent(in)  :: d
      integer(ik),            intent(in)  :: field
      real(wp),               intent(out) :: x(:)
      integer(ik),            intent(out) :: n
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      x = 0.0_wp                         ! slots with no entry read 0
      n = min(d%n, int(size(x), ik))     ! bounded by the caller's array too (see cohort_diag_value)
      do i = 1_ik, n
         if (d%w(i) > tiny_num) x(i) = d%v(field, i) / d%w(i)
      end do
   end subroutine patch_diag_value

   !=======================================================================================!
   !  POLYGON block. Reset at the start of each slow step's fast loop, read as the step mean.   !
   !=======================================================================================!
   subroutine polygon_diag_reset(d)
      type(polygon_diag_block), intent(inout) :: d
      if (.not. d%active) return
      d%v = 0.0_wp
      d%w = 0.0_wp
   end subroutine polygon_diag_reset

   !----- The step mean of one field; 0 if nothing was accumulated, as the patch block reads. --!
   pure real(wp) function polygon_diag_value(d, field) result(x)
      type(polygon_diag_block), intent(in) :: d
      integer(ik),              intent(in) :: field
      x = 0.0_wp
      if (d%active .and. d%w > tiny_num) x = d%v(field) / d%w
   end function polygon_diag_value

end module meds_diag_vars
