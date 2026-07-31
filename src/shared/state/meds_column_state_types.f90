!==========================================================================================!
! meds_column_state_types -- the per-store column DESCRIPTORS of the fast (sub-daily) biophysics !
! loop: the PROGNOSTIC state (canopy-air-space thermal twins + the two soil columns + the snow    !
! store) AND the static per-column PARAMETER bundles (soil_params_t / soil_thermal_params_t) that  !
! describe those same stores, with their `pure` assemblers. They live in `shared` -- the ONE place  !
! reachable by BOTH the biophysics kernels (which mutate the state) AND the demographic state hub    !
! (which will OWN them, per-patch, and thread them through the lockstep reorder). The param builders !
! are state-free constructors (they call only the shared retention curve), so this stays below the   !
! biophysics kernels; meds_biophysics_types re-exports these names so the fast kernels compile        !
! unchanged.                                                                                          !
!                                                                                          !
! Fixed-size (n_soil_layer_max) so the soil columns stay allocatable-free and GPU-eligible.     !
!==========================================================================================!
module meds_column_state_types
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_hydr_lib, only : soil_theta_from_psi, SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   implicit none
   private

   public :: n_soil_layer_max, n_snow_layer_max, N_HYDRO_NODE, LEAF_TEMP_INIT, PSI_INIT
   public :: cas_state_t, soil_column_t, soil_energy_column_t, snow_column_t, soil_carbon_t
   public :: xi_accum_t   !< daily fast->slow accumulator for the soil-carbon matrix (B2)
   public :: blend_cas, blend_soil_w, blend_soil_e, blend_snow, blend_soil_carbon, blend_xi_accum
                                                                    !< area-weighted mix (patch fusion / disturbance seed)
   public :: necromass_to_litter  !< necromass -> litter-destination split (B1; DAG-safe, plain scalars)
   !----- Per-column soil PARAMETER bundles (static geometry/texture) + their pure assemblers.  !
   !      These are NOT prognostic state, but they describe the SAME per-column stores the state !
   !      types below hold, and their builders are stateless parameter constructors -- so they    !
   !      live beside the column-state types (the builders call the shared retention curve).       !
   public :: soil_params_t, soil_thermal_params_t
   public :: build_soil_hydr_params, build_soil_therm_params

   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time soil-column-depth ceiling

   real(wp),    parameter :: PSI_WP = -152.96_wp           !< [m] -1.5 MPa head (wilting-point derivation)

   !----- Per-cohort fast state carried on cohort_block (rides the cohort lockstep). N_HYDRO_NODE !
   !      MUST equal meds_plant_types%N_HYDRO (the psi node count); a fresh cohort starts at a     !
   !      mild tension + a neutral leaf temperature (both relax within one fast step).            !
   integer(ik), parameter :: N_HYDRO_NODE   = 3_ik          !< == N_HYDRO (leaf/wood/root psi nodes)
   real(wp),    parameter :: LEAF_TEMP_INIT = 288.15_wp     !< [K]   fresh-cohort leaf temperature
   real(wp),    parameter :: PSI_INIT       = -0.1_wp       !< [MPa] fresh-cohort node water potential

   !----- Prognostic per-patch soil WATER column (the value the hydrology kernel updates). --!
   type :: soil_column_t
      real(wp) :: theta(n_soil_layer_max) = 0.0_wp   !< [m3/m3] volumetric soil moisture (PROGNOSTIC)
      real(wp) :: w_surface = 0.0_wp                 !< [kg/m2] ponded surface water
      !----- ENTHALPY of the ponded water (issue #78 item 4). The pond used to be a MASS buffer with no  !
      !      thermal state, so water crossing into it shed its enthalpy at the source layer's           !
      !      temperature and that energy left the whole-column ledger -- while the water itself sat on   !
      !      the surface still holding its heat. The books closed (both sides agreed) but the column     !
      !      lost energy that had not gone anywhere. Making this prognostic lets every pond seam be a    !
      !      PAIRED (mass, enthalpy) transfer, the same discipline the snow pack already follows.        !
      !      EXTENSIVE [J/m2], not volumetric, so the transfers are plain additions; temperature is a    !
      !      read-off of uext_to_temp with dry_hcap = 0, exactly as for snow -- which means freeze/thaw  !
      !      of ponded water comes for free. ---------------------------------------------------------!
      real(wp) :: w_surface_enth = 0.0_wp            !< [J/m2] ponded-water internal energy (PROGNOSTIC)
   end type soil_column_t

   !----- Prognostic per-patch soil THERMAL column (internal energy; temp/fliq diagnosed). --!
   type :: soil_energy_column_t
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp    !< [J/m3] volumetric internal energy (PROGNOSTIC)
      real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp    !< [K]    diagnosed each step
      real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp    !< [-]    diagnosed liquid fraction
   end type soil_energy_column_t

   integer(ik), parameter :: n_snow_layer_max = 1_ik        !< MVP single bulk layer; P1 raises to ED2 nzs~8

   !----- Temporary-surface-water / SNOW store: a stacked mass+energy reservoir between the CAS      !
   !      and soil_energy(1). PROGNOSTIC = water-equivalent mass (swe) + EXTENSIVE internal energy    !
   !      (J/m2, unlike soil_energy's J/m3); temperature + liquid fraction are read-offs of           !
   !      uext_to_temp (dry_hcap=0). nlayer=0 is the "no snow" state (store present but empty).        !
   type :: snow_column_t
      real(wp) :: swe(n_snow_layer_max)         = 0.0_wp   !< [kg/m2] water-equivalent mass  (PROGNOSTIC)
      real(wp) :: snow_energy(n_snow_layer_max) = 0.0_wp   !< [J/m2]  extensive internal energy (PROGNOSTIC)
      real(wp) :: snow_depth(n_snow_layer_max)  = 0.0_wp   !< [m]     geometric depth (PROGNOSTIC, = swe/rho_snow in P0)
      real(wp) :: snow_temp(n_snow_layer_max)   = 273.16_wp !< [K]    diagnosed each step (read-off)
      real(wp) :: snow_fliq(n_snow_layer_max)   = 0.0_wp   !< [-]     diagnosed liquid fraction (read-off)
      integer(ik) :: nlayer                     = 0_ik     !< active layer count (0 = no snow)
   end type snow_column_t

   !----- Prognostic per-patch canopy-air-space thermal state (three implicit twins). -------!
   type :: cas_state_t
      real(wp) :: can_enthalpy = 0.0_wp                     !< [J/kg] specific enthalpy (PROGNOSTIC)
      real(wp) :: can_shv      = 0.0_wp                     !< [kg/kg] specific humidity (PROGNOSTIC twin)
      real(wp) :: can_co2      = 400.0_wp                   !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC third twin)
      real(wp) :: can_temp     = 0.0_wp                     !< [K]    diagnosed
      real(wp) :: can_depth    = 20.0_wp                    !< [m]    CAS depth (from canopy height; forcing)
   end type cas_state_t

   !----- Per-column geometry + texture (assembled once per site; ED2 negative-z convention:   !
   !      interface elevations z <= 0 below ground; dz, dz_node are positive magnitudes). Fixed- !
   !      size (n_soil_layer_max) so the hydrology kernel stays allocatable-free and GPU-eligible.!
   type :: soil_params_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0, ED2 slz)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      real(wp) :: vg_alpha(n_soil_layer_max)  = 0.0_wp       !< [1/m] van Genuchten inverse air-entry
      real(wp) :: vg_n(n_soil_layer_max)      = 0.0_wp       !< [-] van Genuchten pore-size index (> 1)
      real(wp) :: psi_sat(n_soil_layer_max)   = 0.0_wp       !< [m] Campbell air-entry potential (option)
      real(wp) :: b_camp(n_soil_layer_max)    = 0.0_wp       !< [-] Campbell exponent (option)
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
      real(wp) :: theta_wp(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] wilting point (DERIVED)
      real(wp) :: root_frac(n_soil_layer_max) = 0.0_wp       !< [-] normalized root fraction (sum = 1)
   end type soil_params_t

   !----- Per-column soil THERMAL texture (geometry + porosity arrive via soil_params_t). ----!
   type :: soil_thermal_params_t
      integer(ik) :: nzg_active = n_soil_layer_max
      real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_thermal_params_t

   !==========================================================================================!
   !  Slow, stateful per-patch soil-carbon pools (written DAILY by meds_soil_biogeochem%          !
   !  soil_carbon_step; READ-ONLY, frozen across the day, by the fast loop's heterotrophic         !
   !  respiration). Lives here (not in meds_biogeochem_types, its conceptual home) so `patch_block` !
   !  (src/core, links `shared` ONLY) can carry it with no core->biogeochemistry library edge --    !
   !  the SAME reason cas_state_t/soil_column_t/soil_energy_column_t/snow_column_t live here.       !
   !  meds_biogeochem_types re-exports this name so the biogeochemistry kernels compile unchanged.  !
   !  The matrix state vector X, ordered litter -> SOM -> passive (n_soil_pool=7 in                 !
   !  meds_biogeochem_types); fast_soil_carbon KEEPS its name/index(2) so meds_fast_ark's frozen-    !
   !  pool Rh reads it as a bare scalar. Lignin is a passive tracer of the structural pools           !
   !  (0 <= L <= C), outside the carbon-mass vector. The optional N twin is present only when         !
   !  opts%n_cycle_on (P1/P2); C-only default (all N fields 0).                                        !
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
   !  Daily fast->slow accumulator for the soil-carbon matrix (MEDS_SLOW_DYNAMICS_DESIGN.md      !
   !  Part II, B2): the per-pool day-INTEGRAL of the environmental decomposition scalar,          !
   !  `xi_int_j = INT_day xi_j dt` [day] (annotated per the author's request -- kept as `xi_int`,  !
   !  not renamed), accumulated once per (patch, fast sub-step) by column_prepass over the SAME     !
   !  frozen pool the fast loop's heterotrophic_respiration_matrix respires against. Reset to 0      !
   !  at the start of each slow step (mirrors cohort%gpp_accum's reset in fast_dynamics), consumed    !
   !  by the daily soil_carbon_step, which decrements each donor pool by dvec_j = xi_int_j*K_j.       !
   !  `rh_fast_accum` is the day's ACCUMULATED fast-loop Rh (audit-only cross-check against            !
   !  soil_carbon_step's own rh_today -- design section 9's rh_seam_gap). Named fields (not an        !
   !  n_soil_pool-sized array) for the SAME reason soil_carbon_t uses named fields: this lives in       !
   !  shared/state so patch_block (core, links shared only) can carry it with no core->biogeochem      !
   !  edge; meds_soil_biogeochem's pack/unpack_pool_vector marshal to/from the array form kernels use.  !
   !==========================================================================================!
   type :: xi_accum_t
      real(wp) :: fast_grnd   = 0.0_wp   !< INT_day xi_1 dt [day]
      real(wp) :: fast_soil   = 0.0_wp   !< INT_day xi_2 dt [day]
      real(wp) :: struct_grnd = 0.0_wp   !< INT_day xi_3 dt [day]
      real(wp) :: struct_soil = 0.0_wp   !< INT_day xi_4 dt [day]
      real(wp) :: microbial   = 0.0_wp   !< INT_day xi_5 dt [day]
      real(wp) :: slow        = 0.0_wp   !< INT_day xi_6 dt [day]
      real(wp) :: passive     = 0.0_wp   !< INT_day xi_7 dt [day]
      real(wp) :: rh_fast_accum = 0.0_wp !< [kgC/m2] today's accumulated fast-loop Rh (audit-only)
   end type xi_accum_t

contains

   !=======================================================================================!
   !  Area-weighted linear mixes: result = w1*a + w2*b. The caller passes NORMALIZED weights !
   !  (w1 = a1/(a1+a2), w2 = a2/(a1+a2)) so an intensive quantity (theta, enthalpy, [J/m3])   !
   !  is conserved on an area basis when two patches fuse or a disturbance gap is carved from  !
   !  its donors. Diagnosed fields (temp/fliq) mix too and are re-diagnosed next fast step.    !
   !=======================================================================================!
   pure function blend_cas(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(cas_state_t), intent(in) :: a, b
      type(cas_state_t)             :: c
      c%can_enthalpy = w1 * a%can_enthalpy + w2 * b%can_enthalpy
      c%can_shv      = w1 * a%can_shv      + w2 * b%can_shv
      c%can_co2      = w1 * a%can_co2      + w2 * b%can_co2
      c%can_temp     = w1 * a%can_temp     + w2 * b%can_temp
      c%can_depth    = w1 * a%can_depth    + w2 * b%can_depth
   end function blend_cas

   pure function blend_soil_w(w1, a, w2, b) result(c)
      real(wp),           intent(in) :: w1, w2
      type(soil_column_t), intent(in) :: a, b
      type(soil_column_t)             :: c
      c%theta     = w1 * a%theta     + w2 * b%theta
      c%w_surface = w1 * a%w_surface + w2 * b%w_surface
      !----- EXTENSIVE, so it blends additively per area exactly like the mass it belongs to. ---------!
      c%w_surface_enth = w1 * a%w_surface_enth + w2 * b%w_surface_enth
   end function blend_soil_w

   pure function blend_soil_e(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_energy_column_t), intent(in) :: a, b
      type(soil_energy_column_t)             :: c
      c%soil_energy = w1 * a%soil_energy + w2 * b%soil_energy
      c%soil_temp   = w1 * a%soil_temp   + w2 * b%soil_temp
      c%soil_fliq   = w1 * a%soil_fliq   + w2 * b%soil_fliq
   end function blend_soil_e

   !----- Area-weighted mix of two snow stores (patch fusion / disturbance seed). The extensive       !
   !      mass + energy blend additively-per-area (like w_surface); depth blends by area weight;       !
   !      temp/fliq are RE-DIAGNOSED by the caller from the blended (swe, snow_energy), never blended.  !
   !      nlayer of the mix = 1 if any blended swe remains, else 0 (caller re-diagnoses / cleans up).   !
   pure function blend_snow(w1, a, w2, b) result(c)
      real(wp),            intent(in) :: w1, w2
      type(snow_column_t), intent(in) :: a, b
      type(snow_column_t)             :: c
      c%swe         = w1 * a%swe         + w2 * b%swe
      c%snow_energy = w1 * a%snow_energy + w2 * b%snow_energy
      c%snow_depth  = w1 * a%snow_depth  + w2 * b%snow_depth
      c%snow_temp   = w1 * a%snow_temp   + w2 * b%snow_temp   ! provisional; caller re-diagnoses from (swe, energy)
      c%snow_fliq   = w1 * a%snow_fliq   + w2 * b%snow_fliq   ! provisional; caller re-diagnoses
      c%nlayer      = max(a%nlayer, b%nlayer)                 ! caller collapses to 0 if blended swe < tiny
   end function blend_snow

   !----- Area-weighted mix of two soil-carbon columns (patch fusion / disturbance seed). Every    !
   !      field is a plain per-area density [kgC/m2] or [kgN/m2] (no diagnosed/re-derived fields,   !
   !      unlike temp/fliq elsewhere), so a straight area-weighted average of all 12 fields          !
   !      conserves total site-wide soil carbon/N exactly, mirroring blend_soil_w/blend_soil_e.       !
   pure function blend_soil_carbon(w1, a, w2, b) result(c)
      real(wp),             intent(in) :: w1, w2
      type(soil_carbon_t),  intent(in) :: a, b
      type(soil_carbon_t)              :: c
      c%fast_grnd_carbon   = w1 * a%fast_grnd_carbon   + w2 * b%fast_grnd_carbon
      c%fast_soil_carbon   = w1 * a%fast_soil_carbon   + w2 * b%fast_soil_carbon
      c%struct_grnd_carbon = w1 * a%struct_grnd_carbon + w2 * b%struct_grnd_carbon
      c%struct_soil_carbon = w1 * a%struct_soil_carbon + w2 * b%struct_soil_carbon
      c%microbial_carbon   = w1 * a%microbial_carbon   + w2 * b%microbial_carbon
      c%slow_carbon        = w1 * a%slow_carbon        + w2 * b%slow_carbon
      c%passive_carbon     = w1 * a%passive_carbon     + w2 * b%passive_carbon
      c%struct_grnd_lignin = w1 * a%struct_grnd_lignin + w2 * b%struct_grnd_lignin
      c%struct_soil_lignin = w1 * a%struct_soil_lignin + w2 * b%struct_soil_lignin
      c%fast_grnd_n        = w1 * a%fast_grnd_n        + w2 * b%fast_grnd_n
      c%fast_soil_n        = w1 * a%fast_soil_n        + w2 * b%fast_soil_n
      c%struct_grnd_n      = w1 * a%struct_grnd_n      + w2 * b%struct_grnd_n
      c%struct_soil_n      = w1 * a%struct_soil_n      + w2 * b%struct_soil_n
      c%mineralized_n      = w1 * a%mineralized_n      + w2 * b%mineralized_n
   end function blend_soil_carbon

   !----- Area-weighted mix of two daily xi accumulators (patch fusion / disturbance seed) --   !
   !      same rationale as blend_cas/blend_soil_w: an intra-day fusion should blend the         !
   !      partial-day accumulation exactly like the other per-patch fast reservoirs. -----------!
   pure function blend_xi_accum(w1, a, w2, b) result(c)
      real(wp),         intent(in) :: w1, w2
      type(xi_accum_t), intent(in) :: a, b
      type(xi_accum_t)             :: c
      c%fast_grnd     = w1 * a%fast_grnd     + w2 * b%fast_grnd
      c%fast_soil     = w1 * a%fast_soil     + w2 * b%fast_soil
      c%struct_grnd   = w1 * a%struct_grnd   + w2 * b%struct_grnd
      c%struct_soil   = w1 * a%struct_soil   + w2 * b%struct_soil
      c%microbial     = w1 * a%microbial     + w2 * b%microbial
      c%slow          = w1 * a%slow          + w2 * b%slow
      c%passive       = w1 * a%passive       + w2 * b%passive
      c%rh_fast_accum = w1 * a%rh_fast_accum + w2 * b%rh_fast_accum
   end function blend_xi_accum

   !=======================================================================================!
   !  NECROMASS -> LITTER-DESTINATION split (MEDS_SLOW_DYNAMICS_DESIGN.md Part II, B1). Every input !
   !  is a plain per-area carbon AMOUNT [kgC/m2] (leaf/storage necromass, fine-root necromass, wood  !
   !  necromass) already converted from the per-plant cohort pools via nplant; every PFT input is a   !
   !  plain scalar trait, so this has NO derived-type dependency (callable from BOTH the driver,      !
   !  which routes the result through litter_input_t/build_litter_input, and the core engine's         !
   !  termination/disturbance kills, which cannot link biogeochemistry -- it adds these six outputs    !
   !  straight onto a soil_carbon_t's carbon+lignin fields). Leaf necromass (+ storage, bundled since  !
   !  both are canopy-associated labile-eligible pools) splits labile/structural by f_labile_leaf,      !
   !  then each of those splits above/below by aboveground_frac; fine-root necromass reuses            !
   !  f_labile_leaf (no separate root trait) but lands ENTIRELY below-ground (no agf split -- fine      !
   !  roots are definitionally below-ground); wood/CWD necromass splits labile/structural by            !
   !  f_labile_stem, then above/below by aboveground_frac. Only the STRUCTURAL streams carry lignin     !
   !  (struct_lignin_frac of each). Pass 0 for any necromass channel that does not apply (e.g. the       !
   !  turnover/shed call site has no wood or storage component).                                         !
   !=======================================================================================!
   pure subroutine necromass_to_litter(leaf_c, root_c, wood_c, storage_c,                            &
                                       f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac,   &
                                       labile_grnd, labile_soil, struct_grnd, struct_soil,             &
                                       lignin_grnd, lignin_soil)
      real(wp), intent(in)  :: leaf_c, root_c, wood_c, storage_c
      real(wp), intent(in)  :: f_labile_leaf, f_labile_stem, aboveground_frac, lignin_frac
      real(wp), intent(out) :: labile_grnd, labile_soil, struct_grnd, struct_soil
      real(wp), intent(out) :: lignin_grnd, lignin_soil
      real(wp) :: canopy_c, canopy_lab, canopy_str, root_lab, root_str, wood_lab, wood_str

      canopy_c   = leaf_c + storage_c
      canopy_lab = canopy_c * f_labile_leaf
      canopy_str = canopy_c * (1.0_wp - f_labile_leaf)
      root_lab   = root_c * f_labile_leaf
      root_str   = root_c * (1.0_wp - f_labile_leaf)
      wood_lab   = wood_c * f_labile_stem
      wood_str   = wood_c * (1.0_wp - f_labile_stem)

      labile_grnd = (canopy_lab + wood_lab) * aboveground_frac
      labile_soil = (canopy_lab + wood_lab) * (1.0_wp - aboveground_frac) + root_lab
      struct_grnd = (canopy_str + wood_str) * aboveground_frac
      struct_soil = (canopy_str + wood_str) * (1.0_wp - aboveground_frac) + root_str
      lignin_grnd = lignin_frac * struct_grnd
      lignin_soil = lignin_frac * struct_soil
   end subroutine necromass_to_litter

   !=======================================================================================!
   !  Per-column soil PARAMETER assemblers. Both are `pure` -- they take plain scalar texture/    !
   !  geometry inputs and fill a parameter struct; they depend on NO soil STATE (theta/energy),  !
   !  so they are state-free constructors, not part of the fast loop. Per-layer texture + TOML    !
   !  wiring land at P3; these standalone builders keep the kernels testable now.                  !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Assemble a per-column soil_params_t: the ED2 negative-z geometry (exponential generator !
   ! OR an explicit soil_layer_z interface array), uniform texture broadcast, an exponential  !
   ! root profile, and the DERIVED thresholds theta_fc/theta_wp (from the retention curve).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_hydr_params(n_active, retention, soil_depth, grid_growth, theta_sat,  &
                                theta_res, ksat, par_a, par_n, root_beta, psi_fc_m, params,    &
                                soil_layer_z_in)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: soil_depth, grid_growth, theta_sat, theta_res
      real(wp),            intent(in)  :: ksat, par_a, par_n, root_beta, psi_fc_m
      type(soil_params_t), intent(out) :: params
      real(wp), optional,  intent(in)  :: soil_layer_z_in(:)
      integer(ik) :: k
      real(wp)    :: denom, rsum

      params%n_active  = n_active
      params%retention = retention

      !----- Interface elevations soil_layer_z(k) <= 0 (k=1 surface, deeper = more negative). -!
      if (present(soil_layer_z_in)) then
         params%soil_layer_z(1:n_active+1) = soil_layer_z_in(1:n_active+1)
      else
         if (abs(grid_growth) < tiny_num) then          ! uniform grid: the exact 0/0 limit of below
            do k = 1_ik, n_active + 1_ik
               params%soil_layer_z(k) = -soil_depth * real(k - 1_ik, wp) / real(n_active, wp)
            end do
         else
            denom = exp(grid_growth) - 1.0_wp
            do k = 1_ik, n_active + 1_ik
               params%soil_layer_z(k) = -soil_depth                                         &
                  * (exp(grid_growth * real(k - 1_ik, wp) / real(n_active, wp)) - 1.0_wp) / denom
            end do
         end if
      end if

      !----- Thicknesses, node elevations, internode spacings (dz, dz_node > 0). ----------!
      do k = 1_ik, n_active
         params%dz(k)     = params%soil_layer_z(k) - params%soil_layer_z(k+1)
         params%z_node(k) = 0.5_wp * (params%soil_layer_z(k) + params%soil_layer_z(k+1))
      end do
      do k = 1_ik, n_active - 1_ik
         params%dz_node(k) = params%z_node(k) - params%z_node(k+1)
      end do
      params%dz_node(n_active) = params%dz(n_active)      ! unused at the bottom face; kept > 0

      !----- Uniform texture broadcast to every active layer. -----------------------------!
      params%theta_sat(1:n_active) = theta_sat
      params%theta_res(1:n_active) = theta_res
      params%ksat(1:n_active)      = ksat
      if (retention == SOIL_RETENTION_CAMPBELL) then
         params%psi_sat(1:n_active) = par_a
         params%b_camp(1:n_active)  = par_n
      else
         params%vg_alpha(1:n_active) = par_a
         params%vg_n(1:n_active)     = par_n
      end if

      !----- Derived thresholds (from the retention curve). -------------------------------!
      params%theta_fc(1:n_active) =                                                          &
         soil_theta_from_psi(retention, psi_fc_m, theta_sat, theta_res, par_a, par_n)
      params%theta_wp(1:n_active) =                                                          &
         soil_theta_from_psi(retention, PSI_WP, theta_sat, theta_res, par_a, par_n)

      !----- Exponential root profile (z_node <= 0 => decays with depth), normalized. -----!
      rsum = 0.0_wp
      do k = 1_ik, n_active
         params%root_frac(k) = exp(root_beta * params%z_node(k)) * params%dz(k)
         rsum = rsum + params%root_frac(k)
      end do
      if (rsum > tiny_num) params%root_frac(1:n_active) = params%root_frac(1:n_active) / rsum
   end subroutine build_soil_hydr_params

   !---------------------------------------------------------------------------------------!
   ! Assemble a soil_thermal_params_t from plain (uniform) inputs -- the thermal twin of      !
   ! build_soil_hydr_params (per-layer texture + TOML land at P3, with hydrology).            !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_therm_params(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid, k_dry, dry_cvol
      type(soil_thermal_params_t), intent(out) :: therm
      therm%nzg_active = n_active
      therm%soil_solid_conductivity(1:n_active) = k_solid
      therm%soil_dry_conductivity(1:n_active)   = k_dry
      therm%soil_dry_heat_capacity(1:n_active)  = dry_cvol
   end subroutine build_soil_therm_params

end module meds_column_state_types
