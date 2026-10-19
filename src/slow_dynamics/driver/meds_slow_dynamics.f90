! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_slow_dynamics -- the THIN slow-tier coordinator (docs/science/order_of_processes.md       !
! §4). Sequences the slow-tier domain drivers as PEERS -- vegetation dynamics, then              !
! biogeochemistry -- and owns the shared per-patch slow-state application                        !
! (update_patch_states) between them. It stays a coordinator, not a numerical stepper: every     !
! actual computation lives in the domain drivers it calls.                                       !
!==========================================================================================!
module meds_slow_dynamics
   use meds_kinds,                 only : wp, ik
   use meds_config,                only : meds_config_t
   use meds_time,                  only : meds_time_t, day_of_year
   use meds_polygon_state_vars, only : polygon_t
   use meds_demography_update, only : update_patch_states
   use meds_vegetation_dynamics,   only : vegetation_dynamics, restructure_stand, shed_excess_tissue_water
   use meds_biogeochem_dynamics,   only : advance_biogeochem_dynamics
   use meds_biogeochem_types, only : litter_input_t, soilc_seam_t
   use meds_patch_state_vars,   only : cas_set_depth
   use meds_slow_ledger,           only : slow_ledger_t, slow_ledger_open, slow_ledger_mark,       &
                                          slow_ledger_rebase, slow_tissue_heat,                   &
                                          slow_ledger_declare, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC, &
                                          SLOW_PHASE_TISSUE_WATER, KGC_PER_UMOL_C
   implicit none
   private

   public :: advance_slow_dynamics, advance_boundary_dynamics, refresh_canopy_depth, settle_tissue_water

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step: vegetation dynamics (demographic rates + apply-primitives; returns  !
   ! this step's litter accumulator), the shared patch-ageing applier, then -- when soil_carbon_  !
   ! on -- the daily soil-carbon matrix step consuming that litter + the fast loop's day-          !
   ! integrated environmental scalar. `step_start` gives leaf phenology its day of year; without  !
   ! it (a bare test call) the phenology is not advanced.                                         !
   !---------------------------------------------------------------------------------------!
   subroutine advance_slow_dynamics(poly, cfg, step_start, seam, ledger, rho_air, latitude_deg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(meds_time_t),   intent(in), optional :: step_start   !< the slow step's start (its day of year)
      type(soilc_seam_t),  intent(inout), optional :: seam   !< per-run seam diagnostics (see the biogeochem driver)
      !----- The site conservation ledger (meds_slow_ledger). Optional so a test, a probe or the   !
      !      C-API can advance a slow step without one; every entry point is a no-op when it is    !
      !      absent or inactive. `rho_air` values the canopy-air store and comes from the fast     !
      !      context.                                                                              !
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp),            intent(in),    optional :: rho_air
      real(wp),            intent(in),    optional :: latitude_deg   !< the polygon's (default [site])

      if (present(ledger)) call slow_ledger_open(ledger, poly, cfg, rho_air)

      if (present(step_start)) then
         call vegetation_dynamics(poly, cfg, day_of_year(step_start), ledger, latitude_deg)
      else
         call vegetation_dynamics(poly, cfg, ledger=ledger, latitude_deg=latitude_deg)
      end if

      call update_patch_states(poly%patch, cfg%dt_years)

      !----- Resize each patch's canopy-air control volume to the stand that now exists. Done HERE, !
      !      after growth and mortality have settled, so the fast loop sees a can_depth that is     !
      !      constant across every sub-step of the coming day (advance_boundary_dynamics resizes    !
      !      again after a restructuring).                                                          !
      call refresh_canopy_depth(poly, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)

      !----- `seam` is intent(inout) and accumulates over the run, so when soil carbon is off it     !
      !      simply stays as it was. ----------------------------------------------------------------!
      if (cfg%soil_carbon_on) call advance_biogeochem_dynamics(poly, cfg, seam, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_SOILC)
   end subroutine advance_slow_dynamics

   !---------------------------------------------------------------------------------------!
   ! advance_boundary_dynamics -- the calendar boundary between two slow steps: restructure the !
   ! stand (restructure_stand), then resize each canopy-air volume to the stand it left, so the !
   ! next fast loop sees a constant can_depth. Called after the output has read the step that   !
   ! ended on the boundary, and before the next one's fast loop; a no-op off a month boundary.  !
   !                                                                                            !
   ! The ledger is NOT reopened: nothing moves the site between the last step's final mark and  !
   ! this call, so the marks below continue that step's window and each operator's residual is  !
   ! still attributed to its own phase. A resumed run's first call precedes any step, so there  !
   ! it only takes the reference store (slow_ledger_rebase).                                    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_boundary_dynamics(poly, cfg, is_new_month, is_new_year, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      type(slow_ledger_t), intent(inout), optional :: ledger
      if (.not. is_new_month) return
      if (present(ledger)) then
         if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
      end if
      call restructure_stand(poly, cfg, is_new_month, is_new_year, ledger)
      call refresh_canopy_depth(poly, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)
   end subroutine advance_boundary_dynamics

   !---------------------------------------------------------------------------------------!
   ! settle_tissue_water -- before the fast loop reads the plants' water, shed what lies above  !
   ! the capacity their biomass now allows (shed_excess_tissue_water) and book it: the water as !
   ! a handoff to the ground, the heat it takes as a change in the tissue's thermal mass. It    !
   ! runs once per step, after the slow step and any calendar boundary have changed the stand,  !
   ! so it continues their ledger window in a phase of its own. Before a run's first step       !
   ! there is no window yet, and it takes the reference store, as a resumed run's boundary does. !
   !---------------------------------------------------------------------------------------!
   subroutine settle_tissue_water(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp) :: water_shed, heat0
      heat0 = 0.0_wp
      if (present(ledger)) then
         if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
         heat0 = slow_tissue_heat(poly)
      end if
      call shed_excess_tissue_water(poly, cfg, water_shed)
      if (present(ledger)) then
         call slow_ledger_declare(ledger, water_out = water_shed, energy_in = slow_tissue_heat(poly) - heat0)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_TISSUE_WATER)
      end if
   end subroutine settle_tissue_water

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_depth -- per-patch CAS depth = tallest cohort height + freeboard, floored.  !
   !                                                                                          !
   ! The ONLY writer of cas%can_depth, called by the slow step, at a calendar boundary and once   !
   ! when a polygon is prepared. The depth sets the canopy air's mass and heat capacity,          !
   ! so it sets how fast the canopy air responds and the freeze-cadence stability limit on         !
   ! dt_fast, and the forcing is moved to its top. A fixed depth would be 4x too much canopy air   !
   ! over a 1 m regenerating gap and ~1.8x too little under a 35 m tropical canopy.                !
   !                                                                                          !
   ! An empty patch falls back to the floor alone (no cohorts, no canopy to sit under).            !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_depth(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp)    :: h_top, depth_new, de, dw, dc, rho
      integer(ik) :: ip, i0, i1, i
      logical     :: book
      if (.not. allocated(poly%patch%cas)) return
      !----- Booking the open-volume exchange needs the air density the CAS is valued at; without  !
      !      a ledger there is nothing to book it into and the resize is a plain geometry update.  !
      book = .false. ; rho = 0.0_wp
      if (present(ledger)) then
         book = ledger%active
         rho  = ledger%rho_air
      end if
      do ip = 1_ik, poly%patch%n
         i0 = poly%patch%cohort_offset(ip)
         i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
         h_top = 0.0_wp
         do i = i0, i1
            h_top = max(h_top, poly%cohort%height(i))
         end do
         depth_new = max(cfg%aero%min_canopy_depth, h_top + cfg%aero%canopy_freeboard)
         if (book) then
            !----- Growing the canopy ENTRAINS air from just above; shrinking DETRAINS it. The     !
            !      intensive state is the invariant, so the extensive content moves, and that is   !
            !      a genuine exchange with the atmosphere -- a DECLARED boundary term, not a leak. ----!
            call cas_set_depth(poly%patch%cas(ip), depth_new, rho_air=rho,                         &
                               de_open=de, dw_open=dw, dc_open=dc)
            call slow_ledger_declare(ledger,                                                       &
                     energy_in=poly%patch%area(ip) * de, water_in=poly%patch%area(ip) * dw,        &
                     carbon_in=poly%patch%area(ip) * dc * KGC_PER_UMOL_C)
         else
            call cas_set_depth(poly%patch%cas(ip), depth_new)
         end if
      end do
   end subroutine refresh_canopy_depth

end module meds_slow_dynamics
