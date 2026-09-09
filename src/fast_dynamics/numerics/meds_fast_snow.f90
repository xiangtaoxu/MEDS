!==========================================================================================!
! meds_fast_snow -- the SHARED pre-column snow stage (MEDS_INTEGRATOR_PARITY.md [RETIRED] row 2 / C4,     !
! GitHub issue #76). Snow used to live inside the retired operator-split driver, which meant      !
! ARK and RK45 imported the snow kernels and never called them: a run with [fast].snow_on = true   !
! silently had no snow under either. This module is the one authority all three dispatch branches   !
! call, so the snow physics is identical by construction rather than by three-way maintenance.      !
!                                                                                          !
! WHY A PRE-COLUMN STAGE, not a member of the column ODE. The pack is advanced over the whole      !
! dt_fast on its own, ahead of the state^n snapshot, and its results become FROZEN boundary        !
! conditions for whatever stepper advances the column afterwards -- a Lie-Trotter operator split.  !
! That is the position snow already occupied on the split path, so this is a relocation of the      !
! approximation rather than a new one, and it is what lets the three schemes share it at all.       !
!                                                                                          !
! What is frozen for the step, snow -> column:                                                      !
!   * snowfac        -- the Niu-Yang cover fraction that area-weights snow vs bare-soil exchange     !
!   * g_base         -- the snow-base series conductance, i.e. the soil's top BC under the pack      !
!   * subl_rate      -- sublimation vapour into the CAS                                              !
!   * melt_rate/t_melt -- meltwater to the PONDING store, with the temperature that values it       !
!   * ground_rad     -- the blended radiative input the whole-column ledgers book                    !
! and column -> snow, LAGGED to state^n: CAS temperature/humidity, the aerodynamic conductance,      !
! and absorbed SW/LW.                                                                                !
!                                                                                          !
! The cost is that snow<->CAS coupling is first order for EVERY scheme, RK45 included -- a           !
! deliberate order reduction on that one coupling. Justified because the pack's timescale is long    !
! relative to dt_fast except for a thin pack (which the snowfac ramp already renders a near-no-op),  !
! and because ED2's own IMEX makes the same call, restricting its implicit matrix to [T_can,         !
! T_veg...] and forward-Eulering hydrology since "some of these processes are nearly binary".        !
! Putting snow INSIDE the tableau needs swe + pack internal energy as column_state_t members with a  !
! column_derivs contribution, and the melt plateau is a non-smooth RHS (T pinned at 0 C) of exactly  !
! the class an embedded-error adaptive controller cannot resolve. That is a separate project.        !
!==========================================================================================!
module meds_fast_snow
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : t_3ple
   use meds_therm_lib,        only : temp_of_liquid_enthalpy
   use meds_biophysics_types, only : snow_env_t, snow_flux_t, snow_melt_t
   use meds_biophysics_opts, only : snow_params_t
   use meds_column_reservoirs, only : snow_column_t
   use meds_ground_biophysics, only : snow_energy_step, snow_accumulate, snow_drain_meltwater,    &
                                      snow_cover_fraction
   implicit none
   private

   public :: snow_stage_t, advance_snow_stage

   !----- The frozen outcome of one pre-column snow advance. Every field is 0/.false. when snow is  !
   !      off or no pack exists, and the consumers are written so that those values reduce their     !
   !      arithmetic EXACTLY to the pre-C4 snow-free form -- which is what makes "snow-off            !
   !      bit-identical" a structural property rather than something to re-verify per scheme. -------!
   type :: snow_stage_t
      logical  :: exists     = .false.   !< a pack is present (drives precip routing + rain_temp)
      real(wp) :: snowfac    = 0.0_wp    !< [-]        Niu-Yang cover fraction actually used
      real(wp) :: h_snow     = 0.0_wp    !< [W/m2]     snowfac-weighted sensible flux to the CAS
      real(wp) :: le_snow    = 0.0_wp    !< [W/m2]     snowfac-weighted latent (sublimation) flux
      real(wp) :: g_base     = 0.0_wp    !< [W/m2]     throttled base conduction into the soil top
      real(wp) :: subl_rate  = 0.0_wp    !< [kg/m2/s]  sublimation vapour source for the CAS
      real(wp) :: melt_rate  = 0.0_wp    !< [kg/m2/s]  meltwater to the ponding store (see t_melt)
      real(wp) :: ground_rad = 0.0_wp    !< [W/m2]     blended ground radiative input for the ledgers
      real(wp) :: acc_enth   = 0.0_wp    !< [J/m2]     precip enthalpy that entered the pack (boundary in)
      real(wp) :: swe0       = 0.0_wp    !< [kg/m2]    pack mass BEFORE the stage (ledger store term)
      real(wp) :: swe1       = 0.0_wp    !< [kg/m2]    pack mass AFTER  the stage (ledger store term)
      real(wp) :: enth0      = 0.0_wp    !< [J/m2]     pack internal energy BEFORE (ledger store term)
      real(wp) :: enth1      = 0.0_wp    !< [J/m2]     pack internal energy AFTER  (ledger store term)
      !----- enthalpy the melt transfer moved pack -> soil layer 1. Needed by any caller whose soil    !
      !      baseline is snapshotted AFTER this stage runs: that snapshot already contains the melt    !
      !      energy while enth0 still contains it too, so the pair double-counts it by exactly this    !
      !      amount. Split snapshots BEFORE the stage and needs no correction. ---------------------!
      real(wp) :: melt_enth  = 0.0_wp    !< [J/m2] melt enthalpy leaving the pack with the meltwater
      !----- Temperature that VALUES the meltwater, i.e. the T with u_liq(T)*melt_mass == melt_enth.     !
      !      The caller hands this to the hydrology kernel as chydro_forcing_t%t_precip so the pond      !
      !      receives exactly melt_enth when it receives melt_rate*dt of mass -- one number, both        !
      !      sides. Falls back to t_3ple when there is no melt mass to value. ------------------------!
      real(wp) :: t_melt     = 0.0_wp    !< [K] effective temperature of the meltwater
   end type snow_stage_t

contains

   !---------------------------------------------------------------------------------------!
   ! advance_snow_stage -- accumulate snowfall + rain-on-snow, advance the snow-surface energy     !
   ! balance at the LAGGED CAS, and drain meltwater to the PONDING store as a PAIRED (mass,        !
   ! enthalpy) transfer. Mutates the pack store only; everything else it reports through st, so the  !
   ! caller decides how the frozen results reach its own stepper. Inputs are the physical boundary   !
   ! quantities, not the driver's aggregates (2026-09 review, item 4 #7).                           !
   !                                                                                          !
   ! The meltwater's enthalpy is NOT handed to the soil here (it was, before issue #78 item 4 gave   !
   ! the pond a thermal state). It leaves the pack via snow_energy and is reported as melt_enth      !
   ! together with t_melt, the temperature that values it; the caller passes t_melt to the           !
   ! hydrology kernel as chydro_forcing_t%t_precip, and the ONE pond inflow carries both halves.     !
   ! Pack and pond are both tracked stores, so the transfer telescopes out of the whole-column       !
   ! ledger rather than needing a boundary term -- and no consumer has to rebase a soil baseline.    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_snow_stage(snow, snow_params, dz_soil_top, abs_sw_ground, abs_lw_ground,        &
                                 snowfall, rainfall, t_air, ggnet, t_soil_top, dt_fast, tcas, qcas,    &
                                 rho, press, st)
      type(snow_column_t),  intent(inout) :: snow           !< the pack store (the ONLY state mutated here)
      type(snow_params_t),  intent(in)    :: snow_params
      real(wp),             intent(in)    :: dz_soil_top    !< [m]       top soil-node depth |z_node(1)|
      real(wp),             intent(in)    :: abs_sw_ground  !< [W/m2]    shortwave reaching the ground
      real(wp),             intent(in)    :: abs_lw_ground  !< [W/m2]    net longwave at the ground
      real(wp),             intent(in)    :: snowfall       !< [kg/m2/s] frozen precipitation
      real(wp),             intent(in)    :: rainfall       !< [kg/m2/s] liquid precipitation (rain-on-snow)
      real(wp),             intent(in)    :: t_air          !< [K]       reference-level air temperature
      real(wp),             intent(in)    :: ggnet          !< [m/s]     ground <-> CAS conductance
      real(wp),             intent(in)    :: t_soil_top     !< [K]       top soil-node temperature
      real(wp),             intent(in)    :: dt_fast, tcas, qcas, rho, press
      type(snow_stage_t),   intent(out)   :: st

      type(snow_env_t)  :: senv
      type(snow_flux_t) :: sfx
      type(snow_melt_t) :: smelt
      real(wp)          :: snow_e0

      !----- default = the bare-ground boundary the snow-free column expects (snowfac = 0). --------!
      !                                                                                             !
      !      ALWAYS-ON. There is no `snow_on` switch any more: snowfall is a boundary water input    !
      !      like rain, and a model that receives it must have somewhere to put it. The flag existed  !
      !      because snow was split-only (C4 shared the stage across all three integrators), and      !
      !      while it existed the DEFAULT (.false.) silently discarded frozen precipitation on the     !
      !      ARK/RK45 paths -- precip_phase splits rain from snow without consulting it, so `off`      !
      !      never meant "no snow", it meant "snow with nowhere to go".                                !
      !                                                                                                !
      !      Always-on costs nothing on a snow-free column: snow_accumulate returns immediately unless  !
      !      a pack exists or the snowfall clears params%min_new_snow_mass, so st stays at the bare-     !
      !      ground defaults set just above, snowfac = 0, and surface_derivs' snow blend reduces         !
      !      EXACTLY to its pre-C4 form. Sub-threshold snowfall onto bare ground still reaches the       !
      !      soil as liquid via the caller's throughfall routing -- nothing is dropped either way. ------!
      st%ground_rad = abs_sw_ground + abs_lw_ground
      st%swe0       = snow%swe(1)         ; st%swe1  = snow%swe(1)
      st%enth0      = snow%snow_energy(1) ; st%enth1 = snow%snow_energy(1)
      snow_e0 = snow%snow_energy(1)
      call snow_accumulate(snow, snowfall, rainfall, t_air, dt_fast, snow_params)
      st%acc_enth = snow%snow_energy(1) - snow_e0   ! precip enthalpy into the pack (boundary in)
      st%exists   = snow%nlayer >= 1_ik             ! accumulate took snow+rain -> precip routing

      if (st%exists .and. snow%swe(1) > snow_params%tiny_snow_mass) then
         !----- SUB-COLUMN: snowfac is snow, (1-snowfac) is bare soil. The pack's boundary exchange   !
         !      is SCALED by snowfac inside snow_energy_step, so a thin/patchy pack barely exchanges  !
         !      -- continuous and stable, with no threshold cliff -- and its returned fluxes are      !
         !      already snowfac-weighted. The bare-soil share is blended by the consumer. -----------!
         st%snowfac       = snow_cover_fraction(snow%swe(1), snow%snow_depth(1), snow_params)
         senv%abs_sw      = abs_sw_ground ; senv%abs_lw = abs_lw_ground
         senv%can_temp    = tcas ; senv%can_shv = qcas ; senv%ggnet = ggnet
         senv%rho_air     = rho ; senv%press = press
         senv%t_soil_top  = t_soil_top
         senv%dz_soil_top = dz_soil_top
         call snow_energy_step(snow, senv, snow_params, dt_fast, st%snowfac, sfx)
         call snow_drain_meltwater(snow, snow_params, smelt)
         st%h_snow  = sfx%h_snow ; st%le_snow = sfx%le_snow ; st%g_base = sfx%g_base
         st%snowfac = sfx%snowfac                         ! the clamped fraction the kernel actually used
         !----- ground radiation boundary in = snow's snowfac-weighted net + bare's (1-snowfac) share !
         st%ground_rad = sfx%rnet + (1.0_wp - st%snowfac) * (abs_sw_ground + abs_lw_ground)
         st%subl_rate  = sfx%w_flux
         st%melt_rate  = (smelt%melt_mass + smelt%dump_mass) / dt_fast
         !----- PAIRED enthalpy: snow store -> soil top (extensive J/m2 -> volumetric J/m3). The mass !
         !      half rides melt_rate into infiltration, and the caller MUST infiltrate it at zero     !
         !      enthalpy (rain_temp = tsupercool_liq) or this enthalpy is counted twice. ------------!
         !----- MELTWATER GOES TO THE POND, not straight into soil layer 1 (issue #78 item 4).           !
         !                                                                                              !
         !      The pack used to hand its melt enthalpy directly to soil_energy(1), and the caller then  !
         !      set rain_temp = tsupercool_liq so the meltwater MASS infiltrated carrying zero enthalpy  !
         !      -- "the energy already moved, paired, here". That worked while the soil was the only     !
         !      place surface water could go. Once the ponding store has a real thermal state the        !
         !      meltwater ponds FIRST and infiltrates from the pond, so the direct transfer would be     !
         !      counted twice: once here, and again in the pond->layer-1 advection at the mixed pond     !
         !      temperature. (That is the C4 double-count in a new guise.)                               !
         !                                                                                              !
         !      So report the enthalpy and the temperature that values it, and let the ONE pond inflow   !
         !      carry both. The pack still loses it (snow_energy was already debited in                  !
         !      snow_drain_meltwater), the pond gains it, and both are tracked stores -- so it           !
         !      telescopes out of the whole-column ledger instead of needing a boundary term. ----------!
         st%melt_enth = smelt%melt_enth + smelt%dump_enth
         st%t_melt    = t_3ple
         if (smelt%melt_mass + smelt%dump_mass > 0.0_wp)                                            &
            st%t_melt = temp_of_liquid_enthalpy(st%melt_enth / (smelt%melt_mass + smelt%dump_mass))
      end if
      st%swe1  = snow%swe(1)
      st%enth1 = snow%snow_energy(1)
   end subroutine advance_snow_stage

end module meds_fast_snow
