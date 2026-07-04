!==========================================================================================!
! meds_hydro_types -- data structures of the plant-hydraulics interface seam.               !
!                                                                                          !
! Four public types describe the stateless per-individual water-transport calculation, plus  !
! the topology / solver / conductance-mode flags:                                            !
!   * hydro_env_t     -- the EXPLICIT boundary conditions + geometry (transpiration demand,    !
!                        soil water potential, rhizosphere conductance, per-plant biomass and   !
!                        canopy geometry). Everything here is read-only.                         !
!   * hydro_params_t  -- a flat, self-contained per-PFT trait set (pressure-volume, xylem         !
!                        vulnerability, whole-plant / segment conductance) so the solver never     !
!                        reaches back into meds_config.                                             !
!   * hydro_opts_t    -- the run selectors + numerical controls (topology, solver, conductance      !
!                        mode, sub-step mode, gravity switch, tolerances).                            !
!   * hydro_flux_t    -- the computed fluxes + updated diagnostics (sapflow, root uptake, leaf/wood   !
!                        water potential, plant loss of conductance, sub-step count, convergence).      !
!                                                                                          !
! The prognostic state (node water potentials psi, [MPa]) is NOT stored here -- it lives in the  !
! cohort SoA and is passed by argument (the FATES *Mem/compute split). These types are pure DATA. !
!==========================================================================================!
module meds_hydro_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t
   public :: N_HYDRO, NODE_LEAF, NODE_STEM, NODE_ROOT, NODE_WOOD
   public :: HYDRO_NODES_2, HYDRO_NODES_3
   public :: HYDRO_SOLVER_EXPM, HYDRO_SOLVER_BE
   public :: HYDRO_COND_KPLANT, HYDRO_COND_SEGMENT
   public :: HYDRO_SUBSTEP_ADAPTIVE, HYDRO_SUBSTEP_FIXED

   !----- Fixed compile-time node dimension (max over supported topologies). The active count  !
   !      is a run option; the state array shape never changes, so it is GPU/SoA-friendly.      !
   integer(ik), parameter :: N_HYDRO   = 3_ik
   integer(ik), parameter :: NODE_LEAF = 1_ik   !< leaf water pool
   integer(ik), parameter :: NODE_STEM = 2_ik   !< stem (3-node) / lumped wood (2-node)
   integer(ik), parameter :: NODE_WOOD = 2_ik   !< alias of NODE_STEM in the 2-node run
   integer(ik), parameter :: NODE_ROOT = 3_ik   !< root (3-node only)

   !----- Topology (number of active nodes). -----------------------------------------------!
   integer(ik), parameter :: HYDRO_NODES_2 = 2_ik   !< leaf + lumped wood (default)
   integer(ik), parameter :: HYDRO_NODES_3 = 3_ik   !< leaf + stem + root (opt-in)

   !----- Sub-step integrator. -------------------------------------------------------------!
   integer(ik), parameter :: HYDRO_SOLVER_EXPM = 1_ik  !< frozen-coefficient matrix exponential (default)
   integer(ik), parameter :: HYDRO_SOLVER_BE   = 2_ik  !< linearly-implicit backward Euler (3-node / stiff)

   !----- Internal-conductance parameterization. -------------------------------------------!
   integer(ik), parameter :: HYDRO_COND_KPLANT  = 1_ik !< leaf-area-specific whole-plant conductance (default)
   integer(ik), parameter :: HYDRO_COND_SEGMENT = 2_ik !< X16 stem allometry: kmax*sap_area/(height*curl)

   !----- Sub-step control. ----------------------------------------------------------------!
   integer(ik), parameter :: HYDRO_SUBSTEP_ADAPTIVE = 1_ik !< step-doubling error control (default)
   integer(ik), parameter :: HYDRO_SUBSTEP_FIXED    = 2_ik !< fixed n_sub equal steps (GPU lockstep)

   !----- Boundary conditions + per-plant geometry (all read-only; per plant). --------------!
   type :: hydro_env_t
      real(wp) :: transp     = 0.0_wp   !< [kg/s]  transpiration demand E (per plant)
      real(wp) :: soil_psi   = 0.0_wp   !< [MPa]   aggregated soil (rhizosphere) water potential
      real(wp) :: rhizo_cond = 0.0_wp   !< [kg/s/MPa] soil->root conductance (per plant, given BC)
      real(wp) :: bleaf      = 0.0_wp   !< [kgC]   leaf biomass  (sets leaf capacitance)
      real(wp) :: bsap       = 0.0_wp   !< [kgC]   sapwood biomass
      real(wp) :: broot      = 0.0_wp   !< [kgC]   fine-root biomass (bsap+broot set wood capacitance)
      real(wp) :: sap_area   = 0.0_wp   !< [m2]    sapwood cross-sectional area (segment cond. mode)
      real(wp) :: height     = 0.0_wp   !< [m]     plant height (gravity head + segment path length)
      real(wp) :: leaf_area  = 0.0_wp   !< [m2]    leaf area (scales whole-plant conductance)
   end type hydro_env_t

   !----- Flat per-PFT hydraulic trait set (self-contained, filled by the seam from cfg%pft). !
   type :: hydro_params_t
      !----- Pressure-volume (Bartlett/Tyree-Hammel), per tissue. --------------------------!
      real(wp) :: leaf_pi0 = 0.0_wp, leaf_eps = 0.0_wp, leaf_af = 0.0_wp  !< [MPa],[MPa],[-]
      real(wp) :: wood_pi0 = 0.0_wp, wood_eps = 0.0_wp, wood_af = 0.0_wp
      real(wp) :: leaf_water_sat = 0.0_wp, wood_water_sat = 0.0_wp        !< [kg H2O / kgC] at saturation
      !----- Xylem vulnerability (loss of conductance). -----------------------------------!
      real(wp) :: wood_psi50 = 0.0_wp   !< [MPa, <0] potential at 50% loss
      real(wp) :: wood_kexp  = 0.0_wp   !< [-]  vulnerability shape (a)
      !----- Conductance parameterization. ------------------------------------------------!
      real(wp) :: k_plant_max = 0.0_wp  !< [kg/s/MPa/m2_leaf] whole-plant (HYDRO_COND_KPLANT)
      real(wp) :: wood_kmax   = 0.0_wp  !< [kg/m/s/MPa] sapwood specific conductivity (HYDRO_COND_SEGMENT)
      real(wp) :: vessel_curl = 1.0_wp  !< [-] tortuosity / path-length factor (HYDRO_COND_SEGMENT)
   end type hydro_params_t

   !----- Run selectors + numerical controls. ----------------------------------------------!
   type :: hydro_opts_t
      integer(ik) :: topology     = HYDRO_NODES_2
      integer(ik) :: solver       = HYDRO_SOLVER_EXPM
      integer(ik) :: cond_mode    = HYDRO_COND_KPLANT
      integer(ik) :: substep_mode = HYDRO_SUBSTEP_ADAPTIVE
      logical     :: gravity_on   = .true.
      real(wp)    :: rtol         = 1.0e-3_wp   !< [-]  relative sub-step tolerance
      real(wp)    :: atol         = 1.0e-3_wp   !< [MPa] absolute sub-step tolerance
      real(wp)    :: h_init       = 0.0_wp      !< [s]  initial sub-step (0 => start at dt)
      integer(ik) :: max_substep  = 200_ik      !< sub-step cap / fixed count
   end type hydro_opts_t

   !----- Outputs (per plant). -------------------------------------------------------------!
   type :: hydro_flux_t
      real(wp)    :: sapflow     = 0.0_wp   !< [kg/s]  wood->leaf sapflow (time-mean over dt)
      real(wp)    :: root_uptake = 0.0_wp   !< [kg/s]  soil->root uptake (time-mean; the budget term)
      real(wp)    :: psi_leaf    = 0.0_wp   !< [MPa]   leaf water potential at end of step
      real(wp)    :: psi_wood    = 0.0_wp   !< [MPa]   wood water potential at end of step
      real(wp)    :: plc         = 0.0_wp   !< [-]     plant loss of conductance (1 - retained)
      integer(ik) :: nsub        = 0_ik     !< sub-steps taken
      logical     :: converged   = .false.  !< .true. if the sub-step cap was not hit
   end type hydro_flux_t

end module meds_hydro_types
