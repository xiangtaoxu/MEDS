!==========================================================================================!
! meds_fast_snow -- the SHARED pre-column snow stage (MEDS_INTEGRATOR_PARITY.md row 2 / C4,     !
! GitHub issue #76). Snow used to live inside meds_fast_split's column_fast_step, which meant     !
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
!   * melt_rate      -- meltwater to infiltration, its enthalpy ALREADY handed to the soil top       !
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
   use meds_constants,        only : tiny_num
   use meds_biophysics_types, only : aero_out_t, patch_biophys_t, snow_env_t, snow_flux_t, snow_melt_t
   use meds_fast_types,       only : column_config_t, column_forcing_t
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
      real(wp) :: melt_rate  = 0.0_wp    !< [kg/m2/s]  meltwater to infiltration (enthalpy already paid)
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
      real(wp) :: melt_enth  = 0.0_wp    !< [J/m2] melt enthalpy transferred to soil layer 1
   end type snow_stage_t

contains

   !---------------------------------------------------------------------------------------!
   ! advance_snow_stage -- accumulate snowfall + rain-on-snow, advance the snow-surface energy     !
   ! balance at the LAGGED CAS, and drain meltwater INTO the soil top as a PAIRED (mass, enthalpy)  !
   ! transfer. Mutates bio%snow and bio%soil_e%soil_energy(1); everything else it reports through   !
   ! st, so the caller decides how the frozen results reach its own stepper.                        !
   !                                                                                          !
   ! MUST be called BEFORE the caller's state^n snapshot, so the melt enthalpy is inside the        !
   ! snapshot and the whole-column energy ledger sees a consistent starting soil column. That       !
   ! ordering was load-bearing on the split path and is load-bearing here.                          !
   !---------------------------------------------------------------------------------------!
   subroutine advance_snow_stage(ccfg, forc, aero, bio, dt_fast, tcas, qcas, rho, press, st)
      type(column_config_t),  intent(in)    :: ccfg
      type(column_forcing_t), intent(in)    :: forc
      type(aero_out_t),       intent(in)    :: aero
      type(patch_biophys_t),  intent(inout) :: bio
      real(wp),               intent(in)    :: dt_fast, tcas, qcas, rho, press
      type(snow_stage_t),     intent(out)   :: st

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
      st%ground_rad = forc%abs_sw_ground + forc%abs_lw_ground
      st%swe0       = bio%snow%swe(1)        ; st%swe1  = bio%snow%swe(1)
      st%enth0      = bio%snow%snow_energy(1) ; st%enth1 = bio%snow%snow_energy(1)
      snow_e0 = bio%snow%snow_energy(1)
      call snow_accumulate(bio%snow, forc%snowf, forc%precip, forc%tair, dt_fast, ccfg%snow)
      st%acc_enth = bio%snow%snow_energy(1) - snow_e0   ! precip enthalpy into the pack (boundary in)
      st%exists   = bio%snow%nlayer >= 1_ik             ! accumulate took snow+rain -> precip routing

      if (st%exists .and. bio%snow%swe(1) > ccfg%snow%tiny_snow_mass) then
         !----- SUB-COLUMN: snowfac is snow, (1-snowfac) is bare soil. The pack's boundary exchange   !
         !      is SCALED by snowfac inside snow_energy_step, so a thin/patchy pack barely exchanges  !
         !      -- continuous and stable, with no threshold cliff -- and its returned fluxes are      !
         !      already snowfac-weighted. The bare-soil share is blended by the consumer. -----------!
         st%snowfac       = snow_cover_fraction(bio%snow%swe(1), bio%snow%snow_depth(1), ccfg%snow)
         senv%abs_sw      = forc%abs_sw_ground ; senv%abs_lw = forc%abs_lw_ground
         senv%can_temp    = tcas ; senv%can_shv = qcas ; senv%ggnet = aero%ggnet
         senv%rho_air     = rho ; senv%press = press
         senv%t_soil_top  = bio%soil_e%soil_temp(1)
         senv%dz_soil_top = max(-ccfg%soil%z_node(1), tiny_num)   ! |z_node(1)| = top-node depth
         call snow_energy_step(bio%snow, senv, ccfg%snow, dt_fast, st%snowfac, sfx)
         call snow_drain_meltwater(bio%snow, ccfg%snow, smelt)
         st%h_snow  = sfx%h_snow ; st%le_snow = sfx%le_snow ; st%g_base = sfx%g_base
         st%snowfac = sfx%snowfac                         ! the clamped fraction the kernel actually used
         !----- ground radiation boundary in = snow's snowfac-weighted net + bare's (1-snowfac) share !
         st%ground_rad = sfx%rnet + (1.0_wp - st%snowfac) * (forc%abs_sw_ground + forc%abs_lw_ground)
         st%subl_rate  = sfx%w_flux
         st%melt_rate  = (smelt%melt_mass + smelt%dump_mass) / dt_fast
         !----- PAIRED enthalpy: snow store -> soil top (extensive J/m2 -> volumetric J/m3). The mass !
         !      half rides melt_rate into infiltration, and the caller MUST infiltrate it at zero     !
         !      enthalpy (rain_temp = tsupercool_liq) or this enthalpy is counted twice. ------------!
         st%melt_enth = smelt%melt_enth + smelt%dump_enth
         bio%soil_e%soil_energy(1) = bio%soil_e%soil_energy(1) + st%melt_enth / ccfg%soil%dz(1)
      end if
      st%swe1  = bio%snow%swe(1)
      st%enth1 = bio%snow%snow_energy(1)
   end subroutine advance_snow_stage

end module meds_fast_snow
