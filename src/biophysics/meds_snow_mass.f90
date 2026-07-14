!==========================================================================================!
! meds_snow_mass -- the MASS side of the MEDS temporary-surface-water / snow store (design      !
! MEDS_SNOW_DESIGN.md P0, single bulk layer). Stateless per-store kernels: snowfall + rain-on-    !
! snow ACCUMULATION (mass + its enthalpy on the shared ice/liquid datum), the Niu-Yang07 snow-    !
! COVER fraction, and meltwater PERCOLATION (drain free liquid above the holding capacity, as a    !
! PAIRED (mass, enthalpy) hand-off to the soil, plus full-melt-out residual dump). Links           !
! meds_shared ONLY. The energy side (surface balance + phase read-off + base conduction) lives in   !
! meds_snow_energy; the fast-loop driver orchestrates accumulate -> energy step -> drain.           !
!==========================================================================================!
module meds_snow_mass
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : t_3ple, tiny_num
   use meds_thermo,            only : uext_to_temp, internal_energy_ice, internal_energy_liquid
   use meds_column_state_types, only : snow_column_t
   use meds_biophysics_types,  only : snow_params_t, snow_melt_t
   implicit none
   private

   public :: snow_cover_fraction, snow_accumulate, snow_drain_meltwater

contains

   !----- Niu-Yang (2007) snow-cover / burial fraction from geometric depth. Monotone tanh in      !
   !      [0,1]; the (rho_snow/rho_fresh) normalization folds into the tuned scale. Used to RAMP     !
   !      the ground optics, the soil-top-BC blend, and the aerodynamic roughness (design §4f/§4g/§6). !
   pure function snow_cover_fraction(swe, snow_depth, params) result(snowfac)
      real(wp),            intent(in) :: swe, snow_depth
      type(snow_params_t), intent(in) :: params
      real(wp) :: snowfac, scale
      if (swe <= params%tiny_snow_mass) then
         snowfac = 0.0_wp
      else
         scale   = max(params%ny07_a * params%z0_snow * (params%rho_snow / 100.0_wp) ** params%ny07_m, tiny_num)
         snowfac = tanh(snow_depth / scale)
      end if
   end function snow_cover_fraction

   !----- Accumulation: frozen precip (snowf) + rain-on-snow (rainf) onto the single bulk layer.    !
   !      snowf lands as ICE at min(t_3ple, tair); rain-on-snow lands as LIQUID at tair and refreezes !
   !      automatically via the inverter downstream. A layer is CREATED only when the total new snow  !
   !      mass reaches min_new_snow_mass; below that (and on bare ground, nlayer=0) NOTHING is added   !
   !      here -- the caller folds sub-threshold snowfall into the soil-top store and routes rain to    !
   !      infiltration (design §4a). Mass + enthalpy are conserved on the shared datum.                 !
   pure subroutine snow_accumulate(snow, snowf, rainf, tair, dt, params)
      type(snow_column_t), intent(inout) :: snow
      real(wp),            intent(in)    :: snowf, rainf, tair, dt
      type(snow_params_t), intent(in)    :: params
      real(wp) :: dm_snow, dm_rain
      logical  :: has_layer
      dm_snow = max(0.0_wp, snowf) * dt
      dm_rain = max(0.0_wp, rainf) * dt
      has_layer = (snow%nlayer >= 1_ik) .or. (dm_snow >= params%min_new_snow_mass)
      if (.not. has_layer) return                          ! bare ground + sub-threshold snow: caller handles it
      snow%swe(1)         = snow%swe(1) + dm_snow + dm_rain
      snow%snow_energy(1) = snow%snow_energy(1)                                                     &
                          + dm_snow * internal_energy_ice(min(t_3ple, tair))                        &
                          + dm_rain * internal_energy_liquid(tair)           ! rain-on-snow: liquid enthalpy, refreezes later
      snow%snow_depth(1)  = snow%snow_depth(1) + dm_snow / params%rho_snow   ! rain fills pores / refreezes -> no bulk depth
      snow%nlayer         = 1_ik
   end subroutine snow_accumulate

   !----- Percolation: after the energy step, drain the FREE liquid (above the holding capacity) as a !
   !      PAIRED (mass, enthalpy) transfer to the soil (melt%melt_*), and handle full melt-out (residual !
   !      dump%). The caller adds melt_mass to infiltration availability AND melt_enth to soil_energy(1)  !
   !      (design §4e). temp/fliq are re-diagnosed from the drained (swe, energy). Mass conserved.        !
   pure subroutine snow_drain_meltwater(snow, params, melt)
      type(snow_column_t), intent(inout) :: snow
      type(snow_params_t), intent(in)    :: params
      type(snow_melt_t),   intent(out)   :: melt
      real(wp) :: wmass_free, hold
      melt = snow_melt_t()
      if (snow%nlayer == 0_ik) return

      !----- Vanished / trace layer: dump the whole residual to the soil, revert to bare ground. -----!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         melt%melted_out = .true.
         call clear_layer(snow)
         return
      end if

      call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))

      !----- Drain the free liquid above the holding capacity (LEAF-3 1:9). ---------------------------!
      hold       = max(1.0_wp - params%liquid_holding_frac, tiny_num)
      wmass_free = max(0.0_wp, snow%swe(1) * (snow%snow_fliq(1) - params%liquid_holding_frac) / hold)
      wmass_free = min(wmass_free, snow%swe(1))
      if (wmass_free > 0.0_wp) then
         melt%melt_mass      = wmass_free
         melt%melt_enth      = wmass_free * internal_energy_liquid(t_3ple)   ! drains at the plateau
         snow%swe(1)         = snow%swe(1) - wmass_free
         snow%snow_energy(1) = snow%snow_energy(1) - melt%melt_enth
         snow%snow_depth(1)  = snow%swe(1) / params%rho_snow
      end if

      !----- If draining took the layer below the melt-out threshold, dump the remainder to soil. -----!
      if (snow%swe(1) <= params%tiny_snow_mass) then
         melt%dump_mass  = max(0.0_wp, snow%swe(1))
         melt%dump_enth  = snow%snow_energy(1)
         melt%melted_out = .true.
         call clear_layer(snow)
      else
         call uext_to_temp(snow%snow_energy(1), snow%swe(1), 0.0_wp, snow%snow_temp(1), snow%snow_fliq(1))
      end if
   end subroutine snow_drain_meltwater

   !----- Revert to the empty "no snow" state (nlayer = 0). ----------------------------------------!
   pure subroutine clear_layer(snow)
      type(snow_column_t), intent(inout) :: snow
      snow%swe(1)         = 0.0_wp
      snow%snow_energy(1) = 0.0_wp
      snow%snow_depth(1)  = 0.0_wp
      snow%snow_fliq(1)   = 0.0_wp
      snow%snow_temp(1)   = t_3ple
      snow%nlayer         = 0_ik
   end subroutine clear_layer

end module meds_snow_mass
