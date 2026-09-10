!==========================================================================================!
! meds_column_state_types -- the PROGNOSTIC per-patch stores of the fast (sub-daily) loop: the !
! canopy-air-space thermal twins, the soil water and soil energy columns, the snow store, the  !
! slow soil-carbon pools, and the daily fast->slow decomposition accumulator -- plus the        !
! area-weighted blends that mix two of each when patches fuse or a disturbance gap is carved.    !
!                                                                                          !
! This is BOUNDARY state: the biophysics kernels mutate it and the demographic state hub owns  !
! it per patch and threads it through the lockstep reorder. That is why it is its own layer     !
! rather than a member of either side. Parameters that describe these same stores live next     !
! door in meds_column_params -- a parameter is not state (see the plan's placement rule 3).     !
!==========================================================================================!
module meds_column_state_types
   use meds_kinds,            only : wp, ik
   use meds_column_params, only : n_soil_layer_max, n_snow_layer_max
   use meds_therm_lib,        only : cas_molar_density
   implicit none
   private

   !----- Below this total air mass a blend has nothing to conserve (see blend_cas). ------!
   real(wp), parameter :: tiny_mass = 1.0e-30_wp

   public :: cas_state_t, cas_set_depth
   public :: soil_column_t, soil_energy_column_t, snow_column_t, soil_carbon_t
   public :: xi_accum_t   !< daily fast->slow accumulator for the soil-carbon matrix (B2)
   public :: blend_cas, blend_soil_w, blend_soil_e, blend_snow, blend_soil_carbon, blend_xi_accum
                                                                    !< area-weighted mix (patch fusion / disturbance seed)

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
      !      read-off of internal_energy_to_temp with dry_hcap = 0, exactly as for snow -- which means freeze/thaw  !
      !      of ponded water comes for free. ---------------------------------------------------------!
      real(wp) :: w_surface_enth = 0.0_wp            !< [J/m2] ponded-water internal energy (PROGNOSTIC)
   end type soil_column_t

   !----- Prognostic per-patch soil THERMAL column (internal energy; temp/fliq diagnosed). --!
   type :: soil_energy_column_t
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp    !< [J/m3] volumetric internal energy (PROGNOSTIC)
      real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp    !< [K]    diagnosed each step
      real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp    !< [-]    diagnosed liquid fraction
   end type soil_energy_column_t

   !----- Temporary-surface-water / SNOW store: a stacked mass+energy reservoir between the CAS      !
   !      and soil_energy(1). PROGNOSTIC = water-equivalent mass (swe) + EXTENSIVE internal energy    !
   !      (J/m2, unlike soil_energy's J/m3); temperature + liquid fraction are read-offs of           !
   !      internal_energy_to_temp (dry_hcap=0). nlayer=0 is the "no snow" state (store present but empty).        !
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
      !----- CAS DEPTH is PROGNOSTIC-ish: a per-patch geometry state the SLOW loop owns, set from   !
      !      the tallest cohort plus a freeboard. It used to be this hardcoded 20 m and NOTHING ever !
      !      assigned it -- the aerodynamics computed the right value and threw it away -- so cas_mass_capacity   !
      !      and cas_molar_capacity were a fixed 24 kg/m2 and 0.83 mol/m2 for every stand, 4x too much air over a  !
      !      1 m regenerating gap and ~1.8x too little over a 35 m tropical canopy. -----------------!
      real(wp) :: can_depth    = 20.0_wp                    !< [m]    CAS depth (slow loop owns it)
   end type cas_state_t

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
   !=========================================================================================!
   ! blend_cas -- mix two canopy-air control volumes on their AIR MASS, not their ground area.   !
   !                                                                                          !
   ! Enthalpy, specific humidity and CO2 mixing ratio are INTENSIVE -- per kg of air. What the   !
   ! merged volume must conserve is the EXTENSIVE content, area x rho x depth x value, and the   !
   ! merged depth is the area-weighted mean (the volumes add over the summed ground). Putting     !
   ! those together, the weight that conserves is area x depth:                                   !
   !                                                                                          !
   !     (a1+a2) . d_new . v_new = a1.d1.v1 + a2.d2.v2 ,   (a1+a2).d_new = a1.d1 + a2.d2          !
   !     =>  v_new = (a1.d1.v1 + a2.d2.v2) / (a1.d1 + a2.d2)                                      !
   !                                                                                          !
   ! Weighting on AREA alone -- what this did until now -- drops the covariance term              !
   ! w1.w2.(d1-d2).(v1-v2), so fusing a tall patch with a gap corrupted the canopy air's energy,   !
   ! humidity AND CO2 together, all three by the same relative error. It was exact only when the   !
   ! two depths matched, which is precisely the case that does not need a blend. Measured at        !
   ! ~1e4 J per patch-fusion event before this (plan §10.2.4, review item 1B #9).                   !
   !                                                                                          !
   ! rho cancels: it is site-uniform, so it never has to be passed in. `can_depth` stays AREA-      !
   ! weighted -- it is the merged volume over the merged ground, which is what a depth is -- and    !
   ! that value equals the total mass weight, so the two are computed once.                         !
   !                                                                                          !
   ! `can_temp` is re-diagnosed from the blended enthalpy by the caller; it is mixed on the same    !
   ! weights only so an un-refreshed read is not wildly wrong.                                      !
   !=========================================================================================!
   pure function blend_cas(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(cas_state_t), intent(in) :: a, b
      type(cas_state_t)             :: c
      real(wp) :: m1, m2, mt
      m1 = w1 * a%can_depth ; m2 = w2 * b%can_depth ; mt = m1 + m2
      c%can_depth = mt                                   ! = w1*d1 + w2*d2, the area-weighted depth
      if (mt > tiny_mass) then
         c%can_enthalpy = (m1 * a%can_enthalpy + m2 * b%can_enthalpy) / mt
         c%can_shv      = (m1 * a%can_shv      + m2 * b%can_shv     ) / mt
         c%can_co2      = (m1 * a%can_co2      + m2 * b%can_co2     ) / mt
         c%can_temp     = (m1 * a%can_temp     + m2 * b%can_temp    ) / mt
      else
         !----- Both volumes are empty: nothing to conserve, so fall back to the area weights     !
         !      rather than divide by zero. The result is carried only until the next refresh.    !
         c%can_enthalpy = w1 * a%can_enthalpy + w2 * b%can_enthalpy
         c%can_shv      = w1 * a%can_shv      + w2 * b%can_shv
         c%can_co2      = w1 * a%can_co2      + w2 * b%can_co2
         c%can_temp     = w1 * a%can_temp     + w2 * b%can_temp
      end if
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

   !=========================================================================================!
   ! cas_set_depth -- resize the canopy-air control volume when the canopy height changes.       !
   !                                                                                          !
   ! The CAS is an OPEN control volume. Growing the canopy does not create air: it enlarges the  !
   ! well-mixed layer by entraining air from just above, and shrinking it detrains air back. The !
   ! entrained/detrained air is at essentially the canopy-air state (it is the air the CAS was    !
   ! already exchanging with, one aerodynamic timescale away), so the INTENSIVE state -- specific !
   ! enthalpy, specific humidity, CO2 mixing ratio -- is the invariant, and it is left untouched. !
   !                                                                                          !
   ! The EXTENSIVE content therefore changes: mass by rho*d(depth), energy by rho*d(depth)*enth. !
   ! That is a real exchange with the atmosphere, not a leak, and the three *_open outputs report !
   ! it -- one per twin, because all three intensive quantities are invariant across the resize    !
   ! and all three extensive contents therefore move together -- so a caller can book it rather    !
   ! than discover it as an unexplained jump. Two reasons it is done HERE and                      !
   ! on the slow step rather than inside the fast loop:                                          !
   !                                                                                          !
   !   * canopy height only changes on a slow step, so the fast ledger never has to carry a      !
   !     moving control volume -- its cas_mass_capacity is constant across every sub-step of a day;           !
   !   * the jumps that matter are not growth (a 0.003 m/day increment is ~7 J/m2, negligible)    !
   !     but DISTURBANCE and FUSION, where a 20 m canopy can become a 1 m gap in one step. Those  !
   !     happen in the slow loop, so this is where the term belongs.                             !
   !                                                                                          !
   ! Conserving TOTAL energy instead (rescaling can_enthalpy by the mass ratio) would be wrong:   !
   ! it would cool the canopy air simply because the trees grew.                                 !
   !=========================================================================================!
   pure subroutine cas_set_depth(cas, depth_new, rho_air, de_open, dw_open, dc_open)
      type(cas_state_t),  intent(inout) :: cas
      real(wp),           intent(in)    :: depth_new  !< [m] new CAS depth
      !----- All optional, and the *_open outputs need `rho_air` to mean anything. A caller that   !
      !      just wants the geometry updated omits them; one booking the exchange asks for the      !
      !      twins it tracks. ------------------------------------------------------------------!
      real(wp), optional, intent(in)    :: rho_air    !< [kg/m3]
      real(wp), optional, intent(out)   :: de_open    !< [J/m2]     energy    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dw_open    !< [kg/m2]    vapour    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dc_open    !< [umol/m2]  CO2       entrained (+) / detrained (-)
      real(wp) :: dmass, dmol
      dmass = 0.0_wp ; dmol = 0.0_wp
      if (present(rho_air)) then
         dmass = rho_air * (depth_new - cas%can_depth)
         dmol  = cas_molar_density(rho_air, cas%can_shv) * (depth_new - cas%can_depth)
      end if
      if (present(de_open)) de_open = dmass * cas%can_enthalpy
      if (present(dw_open)) dw_open = dmass * cas%can_shv
      if (present(dc_open)) dc_open = dmol  * cas%can_co2
      cas%can_depth = depth_new
   end subroutine cas_set_depth

end module meds_column_state_types
