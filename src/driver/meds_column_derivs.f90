!==========================================================================================!
! meds_column_derivs -- the pure whole-column tendency (RHS) for the fast-loop IMEX-ARK          !
! integrator (design docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md, phase P0). The overhaul replaces the       !
! operator-split + Picard fast step with ONE additive Runge-Kutta advance; an ARK needs a         !
! side-effect-free f(y) = dy/dt for the whole column, whereas today every kernel is               !
! *advance-and-commit* (forms a flux, hides it inside a backward-Euler denominator, and mutates   !
! its `intent(inout)` store). This module gives the ARK stepper that RHS.                          !
!                                                                                          !
! TWO ENTRY POINTS:                                                                               !
!   * surface_derivs -- the CAS surface block (leaf-energy diagnostic + ground skin + the three    !
!     CAS twins). Depends only on meds_therm_lib; validated bit-for-bit against the split.            !
!   * column_derivs  -- the WHOLE column: surface_derivs + the soil-heat column, the soil-water    !
!     (Richards) column, and the per-cohort plant-hydraulics 2x2, assembled from the tendency-      !
!     exposing kernels soil_energy_time_deriv / soil_water_time_deriv / plant_water_tendency (each a    !
!     side-effect-free sibling of its BE kernel, reusing the SAME flux helpers -- no re-derivation).!
!                                                                                          !
! Faithfulness: each reservoir's tendency is the EXPLICIT RHS whose backward-Euler advance over    !
! dt reproduces the corresponding split kernel as dt -> 0 (soil heat/water) or exactly (the CAS     !
! implicit-in-atm twins; the 2x2 matrix-exponential). Verified per-reservoir in test_column_derivs.!
! The surface hydrology BCs (q_top, psi_e, soil_psi_root) and the leaf gas-exchange pre-pass are    !
! the frozen, explicit part of the additive split (held constant across the ARK macro-step).       !
!==========================================================================================!
module meds_column_derivs
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap, stefan, cp_air, tiny_num, rho_h2o
   use meds_therm_lib,           only : cas_temp_of_enthalpy, sat_specific_humidity,                    &
                                     sat_specific_humidity_temp_deriv, enthalpy_vapor, uext_to_temp,       &
                                     internal_energy_liquid
   use meds_biophysics_types, only : n_soil_layer_max, soil_energy_column_t, energy_forcing_t,   &
                                     soil_thermal_params_t, soil_params_t, soil_opts_t, energy_opts_t
   use meds_soil_energy,      only : soil_energy_time_deriv
   use meds_soil_water,       only : soil_water_time_deriv
   use meds_cas_biophysics,   only : cas_column_t, cas_source_t, cas_column_time_deriv
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_plant_types,      only : hydro_env_t, hydro_params_t, hydro_opts_t, N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_plant_hydraulics, only : plant_water_tendency
   implicit none
   private

   public :: surface_state_t, surface_frozen_t, surface_tend_t, surface_derivs
   public :: column_state_t, column_frozen_t, column_tend_t, column_derivs
   public :: stage_bflux_t, column_bflux_t

   !----- CAS supersaturation relaxation timescale [s] for the smooth condensation sink: q relaxes to    !
   !      qsat with this e-folding time (dq/dt|_cond = -(q-qsat)/TAU_COND). Short vs dt_fast=1800s so a    !
   !      supersaturated canopy air condenses within the step, but SMOOTH (no clamp) so the ARK2 embedded  !
   !      error stays small. Physical dew/fog forms fast; the exact value is not critical (any strong     !
   !      restoring rate pins q near qsat). ----------------------------------------------------------!
   real(wp), parameter :: TAU_COND = 300.0_wp

   !----- The prognostic CAS surface state advanced by the fast loop. ---------------------------!
   type :: surface_state_t
      real(wp) :: cas_enthalpy = 0.0_wp    !< [J/kg]      canopy-air specific enthalpy
      real(wp) :: cas_shv      = 0.0_wp    !< [kg/kg]     canopy-air specific humidity
      real(wp) :: cas_co2      = 0.0_wp    !< [umol/mol]  canopy-air CO2 mixing ratio
   end type surface_state_t

   !----- Frozen-per-substep inputs to the surface block (pre-pass coefficients, aerodynamic       !
   !      capacities/conductances, atmospheric BCs, lagged radiation + ground-latent forcing,       !
   !      and the soil-water supply fraction). Mirrors what column_fast_step freezes once per        !
   !      sub-step (meds_column_dynamics.f90:256-405).                                               !
   type :: surface_frozen_t
      real(wp), allocatable :: h_coeff_f(:)   !< [W/m2/K]  frozen sensible coefficient
      real(wp), allocatable :: g_tr_f(:)      !< [m/s]     frozen leaf transpiration series conductance
      real(wp), allocatable :: abs_sw(:)      !< [W/m2]    absorbed shortwave (frozen source)
      real(wp), allocatable :: abs_lw(:)      !< [W/m2]    net longwave at the emission base (frozen source)
      real(wp), allocatable :: lai(:)         !< [m2/m2]   cohort leaf area index
      real(wp), allocatable :: h_coeff_w(:)   !< [W/m2/K]  frozen WOOD sensible coefficient (pi*wai*wood_gbh*rho*cp)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] frozen absorbed SW / net LW on wood
      real(wp), allocatable :: wai(:)         !< [m2/m2]   cohort wood area index
      real(wp) :: leaf_emiss    = 0.95_wp     !< [-]       leaf LW emissivity
      real(wp) :: wcap          = 0.0_wp      !< [kg/m2]   CAS mass capacity  -> enthalpy & vapour
      real(wp) :: ccap          = 0.0_wp      !< [mol/m2]  CAS molar capacity -> CO2
      real(wp) :: gah           = 0.0_wp      !< [kg/m2/s] CAS<->atm enthalpy conductance
      real(wp) :: gaw           = 0.0_wp      !< [kg/m2/s] CAS<->atm vapour   conductance
      real(wp) :: gac           = 0.0_wp      !< [mol/m2/s]CAS<->atm CO2      conductance
      real(wp) :: enth_atm      = 0.0_wp      !< [J/kg]    reference-level specific enthalpy
      real(wp) :: shv_atm       = 0.0_wp      !< [kg/kg]   reference-level specific humidity
      real(wp) :: co2_atm       = 400.0_wp    !< [umol/mol]free-atmosphere CO2
      real(wp) :: nee_biotic    = 0.0_wp      !< [umol/m2/s] frozen biotic CO2 source (Ra+Rh-GPP)
      real(wp) :: abs_sw_ground = 0.0_wp      !< [W/m2]    shortwave reaching the ground (frozen source)
      real(wp) :: abs_lw_ground = 0.0_wp      !< [W/m2]    net longwave at the ground (frozen source)
      real(wp) :: ggnet         = 0.0_wp      !< [m/s]     ground<->CAS aerodynamic conductance
      real(wp) :: soil_evap     = 0.0_wp      !< [kg/m2/s] ground latent flux (frozen hydrology authority)
      real(wp) :: rho           = 0.0_wp      !< [kg/m3]   canopy-air density
      real(wp) :: press         = 0.0_wp      !< [Pa]      canopy-air pressure
      real(wp) :: src_frac      = 1.0_wp      !< [-]       soil-water supply fraction (uptake / demand)
      real(wp) :: t_ground      = 0.0_wp      !< [K]       soil-top temperature (diagnosed from the state in column_derivs)
   end type surface_frozen_t

   !----- Surface-block tendencies + the diagnostics the ARK ledger and the soil/hydraulics         !
   !      tendencies consume (coh_qsoil -> soil-heat sink; coh_transp -> soil-water sink; transp_c   !
   !      -> per-cohort hydraulic demand). ------------------------------------------------------!
   type :: surface_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp     !< [J/kg/s]     dH/dt
      real(wp) :: d_cas_shv      = 0.0_wp     !< [kg/kg/s]    dq/dt
      real(wp) :: d_cas_co2      = 0.0_wp     !< [umol/mol/s] dC/dt
      real(wp) :: src_enth       = 0.0_wp     !< [W/m2]     summed surface enthalpy source into the CAS
      real(wp) :: src_vap        = 0.0_wp     !< [kg/m2/s]  summed surface vapour source into the CAS
      real(wp) :: g_top          = 0.0_wp     !< [W/m2]     net energy into the soil-top store
      real(wp) :: h_ground       = 0.0_wp     !< [W/m2]     ground sensible flux to the CAS
      real(wp) :: le_ground      = 0.0_wp     !< [W/m2]     ground latent flux to the CAS
      real(wp) :: coh_rnet       = 0.0_wp     !< [W/m2]     net radiation absorbed by the canopy
      real(wp) :: coh_qsoil      = 0.0_wp     !< [W/m2]     liquid enthalpy the soil sheds (post src_frac)
      real(wp) :: coh_transp     = 0.0_wp     !< [kg/m2/s]  total realized transpiration (post src_frac)
      real(wp) :: cond           = 0.0_wp     !< [kg/m2/s]  smooth condensation sink (dew) draining CAS supersat
      real(wp), allocatable :: leaf_temp(:)   !< [K]        diagnosed per-cohort leaf temperature
      real(wp), allocatable :: wood_temp(:)   !< [K]        diagnosed per-cohort wood temperature
      real(wp), allocatable :: transp_c(:)    !< [kg/m2/s]  per-cohort transpiration DEMAND (pre src_frac)
   end type surface_tend_t

   !----- The full prognostic column state advanced per dt_fast. --------------------------------!
   type :: column_state_t
      real(wp) :: cas_enthalpy = 0.0_wp                    !< [J/kg]
      real(wp) :: cas_shv      = 0.0_wp                    !< [kg/kg]
      real(wp) :: cas_co2      = 0.0_wp                    !< [umol/mol]
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp   !< [J/m3]   per soil layer
      real(wp) :: theta(n_soil_layer_max)       = 0.0_wp   !< [m3/m3]  per soil layer
      real(wp), allocatable :: psi(:,:)                    !< [MPa]    (N_HYDRO, ncoh) leaf/wood water potentials
   end type column_state_t

   !----- Frozen inputs for the whole column: the surface pre-pass + the soil/hydraulics params +   !
   !      the frozen hydrology surface BCs + per-cohort geometry the hydraulics kernel needs.        !
   type :: column_frozen_t
      type(surface_frozen_t)      :: surf         !< the surface-block frozen inputs (t_ground overwritten per call)
      type(soil_params_t)         :: soil         !< soil geometry + texture (dz, root_frac, ...)
      type(soil_thermal_params_t) :: therm        !< soil thermal texture
      type(energy_opts_t)         :: energy_opts  !< soil-thermal options (phase change)
      type(soil_opts_t)           :: hydro_opts   !< soil-water (Richards) options
      type(hydro_params_t)        :: hydro_p      !< plant-hydraulics parameters
      type(hydro_opts_t)          :: hydro_o      !< plant-hydraulics solver options
      real(wp) :: geothermal    = 0.0_wp          !< [W/m2]    bottom heat flux BC
      real(wp) :: q_top         = 0.0_wp          !< [m/s]     Richards top water flux (infiltration - evaporation)
      real(wp) :: soil_psi_root = 0.0_wp          !< [MPa]     root-zone soil water potential (hydraulics BC)
      real(wp) :: rhizo_cond    = 0.0_wp          !< [kg/s/MPa]soil->root conductance (hydraulics BC)
      !----- frozen boundary hydrology for the precip>0 guard-lift: the throughfall/drainage/runoff    !
      !      water carries internal_energy_liquid across the soil boundaries (matches the split's       !
      !      :436-439,518-520 advection), and the scratch column_hydrology_flux's end-of-step ponding/  !
      !      aquifer/water-table is persisted (column_state_t does NOT carry these surface stores). ----!
      real(wp) :: infiltration  = 0.0_wp          !< [kg/m2/s] throughfall reaching the soil top face
      real(wp) :: drainage      = 0.0_wp          !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf   = 0.0_wp          !< [kg/m2/s] surface runoff
      real(wp) :: rain_temp     = 0.0_wp          !< [K]       rain temperature (CAS temp @ state^n)
      real(wp) :: t_bot         = 0.0_wp          !< [K]       bottom-layer soil temperature @ state^n
      real(wp) :: w_surface1    = 0.0_wp          !< [kg/m2]   end-of-step ponded surface water
      real(wp) :: w_aquifer1    = 0.0_wp          !< [kg/m2]   end-of-step aquifer store
      real(wp) :: z_wt1         = 0.0_wp          !< [m]       end-of-step water-table elevation
      real(wp) :: uptake        = 0.0_wp          !< [kg/m2/s] realized root uptake (soil_wat_out ledger term)
      !----- the AUTHORITATIVE end-of-step soil moisture from the scratch column_hydrology_flux (the robust  !
      !      ponding/runoff/free-drain Richards solve). The ARK COMMITS this instead of re-solving theta in   !
      !      the ESDIRK stages (soil water is fully operator-split out; see column_fast_step_ark).            !
      real(wp), allocatable :: theta1(:)          !< [m3/m3]   committed post-step soil moisture (per layer)
      real(wp), allocatable :: psi_e(:)           !< [m]       Zeng-Decker equilibrium potential per layer (frozen)
      !----- per-cohort geometry the hydraulics kernel reads (frozen over the step). ------------!
      real(wp), allocatable :: nplant(:), bleaf(:), bsap(:), broot(:), sap_area(:), height(:), leaf_area(:)
   end type column_frozen_t

   !----- The whole-column tendency vector + diagnostics. ---------------------------------------!
   type :: column_tend_t
      real(wp) :: d_cas_enthalpy = 0.0_wp
      real(wp) :: d_cas_shv      = 0.0_wp
      real(wp) :: d_cas_co2      = 0.0_wp
      real(wp) :: dedt(n_soil_layer_max)   = 0.0_wp   !< [W/m3] dsoil_energy/dt
      real(wp) :: dtheta_dt(n_soil_layer_max) = 0.0_wp!< [1/s]  dtheta/dt
      real(wp), allocatable :: dpsi_dt(:,:)           !< [MPa/s] (N_HYDRO, ncoh)
      real(wp) :: g_top = 0.0_wp, drainage_rate = 0.0_wp, uptake_rate = 0.0_wp
      real(wp), allocatable :: leaf_temp(:)
   end type column_tend_t

   !----- ARK conservation ledger: per-stage boundary-flux RATES (emitted by column_be_stage) and     !
   !      the b-weighted, cross-substep-ACCUMULATED amounts (the time-integral of the true boundary    !
   !      fluxes). Because Y3 - y = (1-gamma)*h*K2 + gamma*h*K3 exactly (BETA*gamma = 1-gamma), the     !
   !      accumulated in/out amounts telescope against the committed store change to machine precision  !
   !      for the flux-form CAS twins + the (energy_resid=0) soil-heat column -- the ARK path can then  !
   !      close the same 7 budgets the split closes. Reflects the CURRENT inert ARK (no soil-boundary   !
   !      water-enthalpy advection); the deferred precip>0 guard-lift adds those terms.                 !
   type :: stage_bflux_t                                    !< per-stage RATES
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp    !< [W/m2]
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp    !< [kg/m2/s]
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp    !< [umol/m2/s]
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp  !< [W/m2]
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp  !< [kg/m2/s]
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp!< [W/m2]
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp!< [kg/m2/s]
   end type stage_bflux_t

   type :: column_bflux_t                                  !< accumulated AMOUNTS (J/m2, kg/m2, umol/m2)
      real(wp) :: cas_enth_in = 0.0_wp, cas_enth_out = 0.0_wp
      real(wp) :: cas_vap_in  = 0.0_wp, cas_vap_out  = 0.0_wp
      real(wp) :: cas_co2_in  = 0.0_wp, cas_co2_out  = 0.0_wp
      real(wp) :: soil_enth_in = 0.0_wp, soil_enth_out = 0.0_wp
      real(wp) :: soil_wat_in  = 0.0_wp, soil_wat_out  = 0.0_wp
      real(wp) :: whole_enth_in = 0.0_wp, whole_enth_out = 0.0_wp
      real(wp) :: whole_wat_in  = 0.0_wp, whole_wat_out  = 0.0_wp
   end type column_bflux_t

contains

   !---------------------------------------------------------------------------------------!
   ! surface_derivs -- the CAS surface-block RHS (leaf-energy diagnostic + ground skin + the three  !
   ! CAS twins). Faithful, side-effect-free transcription of column_fast_step's surface path         !
   ! (meds_column_dynamics.f90:349-412,462). The longwave emission base is the current tcas, so the   !
   ! split's `tcas - te` term is identically zero. Integrating d_cas_* with a backward-Euler-in-the-  !
   ! atmosphere step reproduces the split's committed enth1/shv1/co21 exactly.                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine surface_derivs(y, fro, n, f)
      type(surface_state_t),  intent(in)  :: y
      type(surface_frozen_t), intent(in)  :: fro
      integer(ik),            intent(in)  :: n
      type(surface_tend_t),   intent(out) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat
      real(wp)    :: lw_slope, le_slope, le_ref, dtl, tl, transp_i
      real(wp)    :: lw_slope_w, dtw                                  !< diagnostic WOOD balance (own store)
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      type(cas_source_t) :: cas_src
      type(cas_column_t) :: cas_col
      integer(ik) :: i

      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n))

      tcas   = cas_temp_of_enthalpy(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, fro%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, fro%press)

      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
      do i = 1_ik, n
         lw_slope = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%lai(i)
         le_slope = latent_heat_vap * fro%rho * fro%g_tr_f(i) * dqdt
         le_ref   = latent_heat_vap * fro%rho * fro%g_tr_f(i) * (qsat_c - qcas)
         dtl      = (fro%abs_sw(i) + fro%abs_lw(i) - le_ref)                                     &
                    / max(fro%h_coeff_f(i) + le_slope + lw_slope, tiny_num)
         tl       = tcas + dtl
         f%leaf_temp(i) = tl
         transp_i       = (le_ref + le_slope * dtl) / latent_heat_vap
         f%transp_c(i)  = transp_i                                          ! per-cohort demand (pre src_frac)
         coh_h      = coh_h      + fro%h_coeff_f(i) * dtl
         coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)
         coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)
         coh_transp = coh_transp + transp_i
         coh_rnet   = coh_rnet   + (fro%abs_sw(i) + fro%abs_lw(i) - lw_slope * dtl)
         !----- Diagnostic WOOD balance (own store; emission base = tcas, no transpiration). Wood        !
         !      sensible + net-LW join coh_h / coh_rnet; a diagnostic wood has no storage so the two     !
         !      wood terms are equal (h_coeff_w*dtw) and telescope in the ledger. Frozen wood inputs are !
         !      zero when wood is not diagnostic (build_column_frozen), making this a no-op then.        !
         lw_slope_w = 4.0_wp * fro%leaf_emiss * stefan * tcas ** 3 * fro%wai(i)
         dtw        = (fro%abs_sw_wood(i) + fro%abs_lw_wood(i))                                    &
                      / max(fro%h_coeff_w(i) + lw_slope_w, tiny_num)
         f%wood_temp(i) = tcas + dtw
         coh_h    = coh_h    + fro%h_coeff_w(i) * dtw
         coh_rnet = coh_rnet + (fro%abs_sw_wood(i) + fro%abs_lw_wood(i) - lw_slope_w * dtw)
      end do
      coh_qw     = coh_qw     * fro%src_frac
      coh_qsoil  = coh_qsoil  * fro%src_frac
      coh_transp = coh_transp * fro%src_frac

      call ground_surface_fluxes(fro%t_ground, tcas, fro%ggnet, fro%rho, fro%soil_evap, f%h_ground, f%le_ground)
      f%g_top     = fro%abs_sw_ground + fro%abs_lw_ground - f%h_ground - f%le_ground

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + fro%soil_evap
      f%coh_rnet   = coh_rnet
      f%coh_qsoil  = coh_qsoil
      f%coh_transp = coh_transp

      !----- SMOOTH condensation sink (dew/fog): when the CAS is supersaturated, condense the excess as a  !
      !      continuous relaxation of q toward qsat with timescale TAU_COND -- C1-smooth except a mild kink !
      !      at RH=1 (vs the removed HARD state clamp). The condensed vapour leaves as liquid at Tcas: the  !
      !      CAS loses vapour mass (src_vap) AND only the LIQUID enthalpy (src_enth -= cond*u_liq); the     !
      !      latent heat (enthalpy_vapor - u_liq) STAYS in the air (warming it, as real condensation does). !
      !      Replaces the discontinuous clamp -> keeps the ARK2 embedded error smooth through saturation.   !
      f%cond     = (fro%wcap / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
      f%src_vap  = f%src_vap  - f%cond
      f%src_enth = f%src_enth - f%cond * internal_energy_liquid(tcas)

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_src%surface_enthalpy_source = f%src_enth
      cas_src%surface_vapor_source    = f%src_vap
      cas_src%biotic_co2_source       = fro%nee_biotic
      cas_col%air_mass_capacity        = fro%wcap
      cas_col%air_molar_capacity       = fro%ccap
      cas_col%atm_conductance_enthalpy = fro%gah
      cas_col%atm_conductance_vapor    = fro%gaw
      cas_col%atm_conductance_co2      = fro%gac
      cas_col%atm_enthalpy             = fro%enth_atm
      cas_col%atm_specific_humidity    = fro%shv_atm
      cas_col%atm_co2                  = fro%co2_atm
      call cas_column_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_src, cas_col,        &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !---------------------------------------------------------------------------------------!
   ! column_derivs -- the WHOLE-column RHS. Diagnoses the soil-top temperature from the state (so    !
   ! the ground skin couples to the current soil-top energy), runs the surface block, then assembles !
   ! the soil-heat, soil-water, and per-cohort hydraulics tendencies from the frozen forcing + the   !
   ! surface couplings (coh_qsoil -> root heat sink, coh_transp -> root water sink, transp_c -> the   !
   ! per-cohort transpiration demand). Commits nothing. This mirrors column_fast_step's operator      !
   ! sequence WITHOUT the backward-Euler denominators -- it is the tendency an IMEX-ARK stage evaluates.!
   !---------------------------------------------------------------------------------------!
   pure subroutine column_derivs(y, fro, n, nsl, f)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(column_tend_t),   intent(out) :: f

      type(surface_state_t)      :: ys
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: sf
      type(soil_energy_column_t) :: soil_e
      type(energy_forcing_t)     :: eforc
      type(hydro_env_t)          :: henv
      real(wp)                   :: t_ground, fliq1, wmass1, root_uptake(n_soil_layer_max)
      real(wp)                   :: dpsi(N_HYDRO)
      integer(ik)                :: k, i

      allocate(f%dpsi_dt(N_HYDRO, n), f%leaf_temp(n))

      !----- Diagnose the soil-top temperature from the current state so the ground skin sees the   !
      !      prognostic soil-top energy (the coupling the surface block needs). ---------------------!
      wmass1   = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, fro%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      fs = fro%surf ; fs%t_ground = t_ground
      call surface_derivs(ys, fs, n, sf)
      f%d_cas_enthalpy = sf%d_cas_enthalpy ; f%d_cas_shv = sf%d_cas_shv ; f%d_cas_co2 = sf%d_cas_co2
      f%g_top = sf%g_top ; f%leaf_temp(1:n) = sf%leaf_temp(1:n)

      !----- 2. Soil-heat column: g_top from the surface, root heat sink from the shed enthalpy. --!
      soil_e%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top      = sf%g_top
      eforc%geothermal = fro%geothermal
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%root_heat_sink(k) = sf%coh_qsoil * fro%soil%root_frac(k)
         eforc%w_flux(k)         = 0.0_wp                       ! interior advection lumped (baseline)
      end do
      call soil_energy_time_deriv(soil_e, eforc, fro%therm, fro%soil, fro%energy_opts, f%dedt)

      !----- 3. Soil-water (Richards) column: frozen top flux + psi-limited root sink. ---------!
      do k = 1_ik, nsl
         root_uptake(k) = sf%coh_transp * fro%soil%root_frac(k)
      end do
      call soil_water_time_deriv(y%theta, fro%soil, fro%hydro_opts, nsl, fro%q_top, fro%psi_e,     &
                               root_uptake, f%dtheta_dt, f%drainage_rate, f%uptake_rate)

      !----- 4. Per-cohort plant hydraulics (2x2): the transpiration demand drives leaf<->wood<->soil. !
      do i = 1_ik, n
         henv%transp     = sf%transp_c(i) * fro%surf%src_frac / max(fro%nplant(i), tiny_num)
         henv%soil_psi   = fro%soil_psi_root
         henv%rhizo_cond = fro%rhizo_cond
         henv%bleaf      = fro%bleaf(i) ; henv%bsap = fro%bsap(i) ; henv%broot = fro%broot(i)
         henv%sap_area   = fro%sap_area(i) ; henv%height = fro%height(i) ; henv%leaf_area = fro%leaf_area(i)
         call plant_water_tendency(henv, fro%hydro_p, fro%hydro_o, y%psi(:,i), dpsi)
         f%dpsi_dt(:,i) = dpsi
      end do
   end subroutine column_derivs

end module meds_column_derivs
