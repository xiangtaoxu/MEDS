!==========================================================================================!
! meds_soil_types -- the argument records of the ground-column kernels: soil hydrology, soil  !
! thermal energy, and the snow store. Pure DATA, no methods, no hidden state.                 !
!                                                                                          !
! Each domain follows the same shape: a *_forcing_t / *_env_t the kernel reads, and a *_flux_t !
! it returns, with the PROGNOSTIC stores themselves in state/column and the run-config          !
! selectors in meds_biophysics_opts. Split out of the former meds_canopy_types/meds_soil_types.             !
!==========================================================================================!
module meds_soil_types
   use meds_kinds,             only : wp, ik
   use meds_column_params,  only : n_soil_layer_max, n_snow_layer_max
   implicit none
   private

   public :: chydro_forcing_t, chydro_flux_t
   public :: energy_forcing_t, energy_flux_t
   public :: snow_env_t, snow_flux_t, snow_melt_t

   !=======================================================================================!
   !  Soil-column hydrology types + selector codes (meds_soil_water, design §4).  !
   !  Fixed-size (n_soil_layer_max) so the kernel stays allocatable-free and GPU-eligible.   !
   !  ED2 negative-z convention: elevation z <= 0 below ground; dz, dz_node are positive       !
   !  magnitudes.                                                                               !
   !=======================================================================================!
   !----- n_soil_layer_max + the prognostic column-state types live in meds_column_state_types    !
   !      (src/shared); the SOIL_* solver selectors + soil_opts_t live in meds_biophysics_opts     !
   !      (shared/config); the constitutive SOIL_RETENTION_* live in meds_hydr_lib. All re-exported !
   !      below so the fast kernels + callers keep `use meds_canopy_types/meds_soil_types` unchanged. ----------!

   !----- Soil-column boundary conditions (read-only). ------------------------------------!
   type :: chydro_forcing_t
      real(wp) :: precip_ground = 0.0_wp                  !< [kg/m2/s] ground-reaching liquid (post interception)
      real(wp) :: root_uptake(n_soil_layer_max) = 0.0_wp  !< [kg/m2/s] per-layer transpiration DEMAND (x nplant)
      real(wp) :: t_ground = 298.15_wp                    !< [K] ground skin temp (FORCED = T_air until soil energy)
      !----- PER-LAYER soil temperature and the rainfall temperature, needed because this kernel now owns  !
      !      the ponding store's ENTHALPY as well as its mass (issue #78 item 4). Valuing the saturation  !
      !      clip requires layer k's own temperature, and valuing the rain that ponds requires the        !
      !      rainfall temperature. Keeping mass here and enthalpy in the callers is what produced the two   !
      !      defects fixed in PR #81 (the ARK condensate deposit and the RK45 double-clip): a store whose  !
      !      two halves are owned in different places drifts. Defaults make the enthalpy terms 0, so a     !
      !      caller that does not set them gets the pre-#78 mass-only behaviour. -------------------------!
      !----- Defaulted to a PHYSICAL temperature, not 0: internal_energy_liquid is referenced to
      !      tsupercool_liq (~57 K), so a 0 K default would value every pond transfer at a large
      !      NEGATIVE enthalpy rather than at zero. A caller that leaves these alone gets a
      !      self-consistent (if arbitrary) thermal treatment and unchanged mass behaviour; a caller
      !      that wants meaningful pond enthalpy must set both.
      real(wp) :: soil_temp(n_soil_layer_max) = 298.15_wp !< [K] per-layer soil temperature (clip enthalpy)
      real(wp) :: t_pond_inflow = 298.15_wp                    !< [K] temperature of precip_ground (pond inflow)
      real(wp) :: q_air    = 0.0_wp                       !< [kg/kg] canopy-air specific humidity (soil evap)
      real(wp) :: rho_air  = 1.2_wp                       !< [kg/m3] canopy-air density (soil evap)
      real(wp) :: r_aero   = 100.0_wp                     !< [s/m] aerodynamic resistance of the BARE-SOIL tile
                                                          !<       (1/ggnet) -- NOT snow-inflated; see below.
      real(wp) :: snow_free_frac = 1.0_wp                 !< [-] snow-free AREA fraction (1-snowfac). Ground
                                                          !<     evaporation is an AREA-weighted tile flux:
                                                          !<     E = snow_free_frac * rho*(q_g-q_air)/(r_aero+r_soil).
                                                          !<     Throttling partial cover by inflating r_aero
                                                          !<     instead is NOT equivalent whenever the dry-surface
                                                          !<     -layer resistance r_soil > 0 -- it divides only the
                                                          !<     aerodynamic leg by (1-snowfac), so it over-predicts
                                                          !<     bare-soil evaporation at intermediate cover, worst
                                                          !<     where a dry surface makes r_soil dominant. The two
                                                          !<     agree exactly at snowfac = 0 and snowfac = 1.
   end type chydro_forcing_t

   !----- soil_params_t (per-column geometry + texture) is defined in (and re-exported from)    !
   !      meds_column_state_types, beside the prognostic soil columns it describes. -------------!

   !----- soil_opts_t (soil-water solver selectors + tolerances) lives in meds_biophysics_opts    !
   !      (shared/config); re-exported above.                                                     !

   !----- Soil-column outputs + diagnostics. ----------------------------------------------!
   type :: chydro_flux_t
      real(wp) :: infiltration   = 0.0_wp                !< [kg/m2/s] top-face infiltration
      real(wp) :: drainage       = 0.0_wp                !< [kg/m2/s] bottom-face drainage
      real(wp) :: runoff_surf    = 0.0_wp                !< [kg/m2/s] surface runoff
      real(wp) :: soil_evap      = 0.0_wp                !< [kg/m2/s] ground evaporation
      real(wp) :: uptake_total   = 0.0_wp                !< [kg/m2/s] realized root uptake (after theta_wp cap)
      real(wp) :: uptake_deficit = 0.0_wp                !< [kg/m2/s] capped (unmet) sink
      real(wp) :: clip_excess    = 0.0_wp                !< [kg/m2/s] theta-clip water routed to ponding
      !----- POST-SOLVE mass corrections, PER LAYER. These move water with NO face, so the soil ENERGY  !
      !      column cannot see them: it consumes the corrected theta but advects enthalpy only on the    !
      !      faces. Because internal_energy_liquid carries the tsupercool_liq datum (~1.0 MJ/kg in       !
      !      ABSOLUTE terms), an uncompensated correction lands entirely in the diagnosed temperature    !
      !      -- of order 3 K per kg/m2 in a 0.1 m layer. Exported per layer so the caller can debit /    !
      !      credit each layer's enthalpy at ITS OWN temperature, leaving T unchanged (the correction    !
      !      is numerical, so it must be temperature-NEUTRAL). ------------------------------------------!
      real(wp) :: clip_layer(n_soil_layer_max)  = 0.0_wp !< [kg/m2/s] saturation-clip water LEAVING layer k
                                                         !<           for the ponding store (>= 0)
      real(wp) :: floor_layer(n_soil_layer_max) = 0.0_wp !< [kg/m2/s] theta_res hard-floor water CREATED in
                                                         !<           layer k (>= 0; the last-resort guard
                                                         !<           that also shows up in mass_resid)
      real(wp) :: face_mass_resid = 0.0_wp              !< [kg/m2] |per-layer mass change - net face flux|, summed.
                                                        !<   The interior-face contract the soil ENERGY column
                                                        !<   relies on: w_flux must carry the mass that actually
                                                        !<   moved, or the enthalpy advected on it is fiction.
                                                        !<   mass_resid cannot see this -- interior face errors
                                                        !<   cancel in a column-vs-boundary sum.
      !----- PONDING-STORE ENTHALPY exports (issue #78 item 4). The pond is now a real thermal store, so   !
      !      its seams stop being boundary losses and become paired transfers:                             !
      !        * t_infil  -- the temperature of the water that infiltrates. It comes OUT OF THE POND        !
      !          (rain enters the pond first, then infiltration draws from the mixture), so the soil's      !
      !          top-face advection must use this, NOT t_film_valuation. With a dry pond it IS the rainfall           !
          !          temperature, so the common case is unchanged.                                          !
      !        * runoff_enth -- runoff is a genuine boundary energy OUTPUT now that the water it carries     !
      !          had a temperature. Covers both the Dunne share (at t_pond_inflow, never entered the pond) and    !
      !          the pond overflow (at the pond temperature).                                                !
      !        * clip_layer's enthalpy is NO LONGER a boundary loss: it moves layer k -> pond, and both       !
      !          ends are tracked stores, so it telescopes out of the whole-column ledger entirely. ---------!
      real(wp) :: t_infil     = 0.0_wp                   !< [K] temperature of the infiltrating water
      real(wp) :: runoff_enth = 0.0_wp                   !< [W/m2] enthalpy leaving with surface runoff
      real(wp) :: psi_soil(n_soil_layer_max) = 0.0_wp    !< [MPa] per-layer matric potential (EXPORTED to hydraulics)
      real(wp) :: w_flux(n_soil_layer_max)   = 0.0_wp    !< [m/s] time-mean DOWNWARD Darcy flux BELOW node k (k=1..n-1);
                                                         !<       interior interfaces only (EXPORTED for advective heat)
      real(wp) :: mass_resid     = 0.0_wp                !< [kg/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik                         !< sub-steps taken
      logical  :: converged = .true.                     !< .false. on any cap-hit
   end type chydro_flux_t

   !----- Soil-column thermal boundary conditions (read-only). ------------------------------!
   type :: energy_forcing_t
      real(wp) :: g_top      = 0.0_wp                       !< [W/m2] net ground heat flux (Rn-H-LE), top Neumann
      real(wp) :: geothermal = 0.0_wp                       !< [W/m2] bottom flux (default 0)
      real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp     !< [m3/m3] theta from hydrology (kappa, C_eff)
      real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp     !< [m/s]  inter-layer water flux for advective heat,
                                                            !<        UPWARD-positive (matches the hf face convention;
                                                            !<        the caller flips the DOWNWARD-positive hydrology
                                                            !<        flux -- see meds_fast_ark.f90:
                                                            !<        eforc%w_flux = -hflux%w_flux). Populating it with
                                                            !<        the raw downward flux would reverse the advection.
      real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp !< [W/m2] enthalpy removed with root uptake
      !----- BOUNDARY water fluxes, same UPWARD-positive convention as w_flux above, with the       !
      !      temperature of water arriving from OUTSIDE the column. These close the water-enthalpy    !
      !      advection at the two boundary faces so ALL of it -- top, interior, bottom -- is applied   !
      !      by ONE upwind rule at ONE time level. Previously the driver added the top/bottom terms    !
      !      to soil_energy directly, BEFORE the step: the inflow was evaluated on state^n while the   !
      !      interior outflow used the post-conduction T^{n+1}. With internal_energy_liquid ~1.0       !
      !      MJ/kg (the tsupercool_liq datum), each face term reaches ~1300 W/m2 for a few mm/h of     !
      !      percolation, so the physical signal is the small DIFFERENCE of two very large numbers --   !
      !      and a 5 K time-level inconsistency in one of them is worth ~30 W/m2, i.e. the entire       !
      !      signal. Evaluating both consistently is what makes the scheme well posed. ------------------!
      real(wp) :: w_flux_top = 0.0_wp                       !< [m/s] water across the TOP face (<0 = infiltration in)
      real(wp) :: w_flux_bot = 0.0_wp                       !< [m/s] water across the BOTTOM face (<0 = drainage out)
      real(wp) :: t_water_top = 0.0_wp                      !< [K] temperature of water entering from above
      real(wp) :: t_water_bot = 0.0_wp                      !< [K] temperature of water entering from below
   end type energy_forcing_t

   !----- Soil-column energy outputs + diagnostics. ----------------------------------------!
   type :: energy_flux_t
      real(wp) :: ground_heat  = 0.0_wp                     !< [W/m2] conductive flux into layer 1
      real(wp) :: bottom_heat  = 0.0_wp                     !< [W/m2] advective+geothermal bottom loss
      real(wp) :: energy_resid = 0.0_wp                     !< [J/m2] closed-budget residual (~0)
      integer(ik) :: nsub = 0_ik
      logical  :: converged = .true.
   end type energy_flux_t

   !=======================================================================================!
   !  Canopy-air-space CO2 balance: the prognostic third twin is carried in cas_state_t and       !
   !  advanced by meds_cas_biophysics's cas_column_* kernels; the Rh selectors (HR_*) + co2_opts_t !
   !  / damm_params_t live in meds_biogeochem_types (heterotrophic respiration is decomposition).  !
   !=======================================================================================!

   !=======================================================================================!
   !  SNOW / temporary-surface-water types (meds_ground_biophysics; design MEDS_SNOW_DESIGN.md P0). !
   !  STATELESS per-store kernels; forcing arrives as value types. (snow_params_t -- the physical  !
   !  parameter table -- lives in meds_biophysics_opts, shared/config; re-exported above.)          !
   !=======================================================================================!

   !----- Surface boundary conditions for the snow-surface energy balance (read-only, FORCED). --!
   type :: snow_env_t
      real(wp) :: abs_sw   = 0.0_wp, abs_lw = 0.0_wp    !< [W/m2] absorbed SW (snow albedo), NET LW at the snow surface
      real(wp) :: can_temp = 273.16_wp, can_shv = 0.0_wp !< [K],[kg/kg] CAS state (forced sibling)
      real(wp) :: ggnet    = 0.0_wp                     !< [m/s] ground<->CAS conductance (heat = vapour)
      real(wp) :: rho_air  = 1.2_wp, press = 101325.0_wp !< [kg/m3],[Pa] CAS air
      real(wp) :: t_soil_top = 273.16_wp                !< [K] top soil-node temperature (snow-base conduction target)
      real(wp) :: k_soil_top = 1.5_wp                   !< [W/m/K] top soil thermal conductivity (k_face geom-mean)
      real(wp) :: dz_soil_top = 0.05_wp                 !< [m] top soil node depth |z_node(1)| (interface spacing)
   end type snow_env_t

   !----- Snow-surface energy-balance outputs (all positive UPWARD; per unit ground area). ------!
   type :: snow_flux_t
      real(wp) :: t_surf  = 273.16_wp, fliq = 0.0_wp    !< [K],[-] diagnosed snow-surface temperature + liquid fraction
      real(wp) :: h_snow  = 0.0_wp, le_snow = 0.0_wp    !< [W/m2] sensible + latent (sublimation/evap) to the CAS
      real(wp) :: w_flux  = 0.0_wp                      !< [kg/m2/s] vapour mass to CAS (>0 sublimation, <0 deposition)
      real(wp) :: g_base  = 0.0_wp                      !< [W/m2] conduction snow-base -> soil top (the new soil top BC)
      real(wp) :: rnet    = 0.0_wp                      !< [W/m2] net all-wave radiation absorbed by the snow surface
      real(wp) :: snowfac = 0.0_wp                      !< [-] Niu-Yang07 snow-cover / burial fraction
      real(wp) :: energy_resid = 0.0_wp                 !< [J/m2] closed-budget residual (~0 by construction)
   end type snow_flux_t

   !----- Snow-mass outputs from accumulation + melt drainage (a store->store handoff to soil). -!
   type :: snow_melt_t
      real(wp) :: melt_mass = 0.0_wp                    !< [kg/m2] free liquid drained this step -> soil infiltration
      real(wp) :: melt_enth = 0.0_wp                    !< [J/m2]  matching enthalpy (internal_energy_liquid(t_3ple))
      real(wp) :: dump_mass = 0.0_wp                    !< [kg/m2] full-melt-out residual water dumped to the soil
      real(wp) :: dump_enth = 0.0_wp                    !< [J/m2]  full-melt-out residual energy dumped to the soil
      logical  :: melted_out = .false.                  !< .true. when the layer vanished (nlayer -> 0)
      real(wp) :: mass_resid = 0.0_wp                   !< [kg/m2] closed-budget residual (~0 by construction)
   end type snow_melt_t

end module meds_soil_types
