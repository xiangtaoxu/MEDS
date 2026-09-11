!==========================================================================================!
! meds_slow_dynamics -- the THIN slow-tier coordinator (MEDS_SLOW_DYNAMICS_DESIGN.md Part II,   !
! section 10a). Sequences the slow-tier domain drivers as PEERS -- vegetation dynamics, then     !
! biogeochemistry -- and owns the shared per-patch slow-state application (update_patch_states,  !
! hoisted out of meds_vegetation_dynamics) between them. It stays a coordinator, not a numerical   !
! stepper: every actual computation lives in the domain drivers it calls.                          !
!==========================================================================================!
module meds_slow_dynamics
   use meds_kinds,                 only : wp, ik
   use meds_config,                only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_demography_update, only : update_patch_states
   use meds_vegetation_dynamics,   only : vegetation_dynamics
   use meds_biogeochem_dynamics,   only : advance_biogeochem_dynamics
   use meds_biogeochem_types, only : litter_input_t
   use meds_column_state_types, only : cas_set_depth
   use meds_slow_ledger,           only : slow_ledger_t, slow_ledger_open, slow_ledger_mark,       &
                                          slow_ledger_declare, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC, &
                                          KGC_PER_UMOL_C
   implicit none
   private

   public :: advance_slow_dynamics

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step: vegetation dynamics (demographic rates + apply-primitives; returns  !
   ! this step's litter accumulator), the shared patch-ageing applier, then -- when soil_carbon_  !
   ! on -- the daily soil-carbon matrix step consuming that litter + the fast loop's day-          !
   ! integrated environmental scalar. Same signature as vegetation_dynamics (doy optional).        !
   !---------------------------------------------------------------------------------------!
   subroutine advance_slow_dynamics(site, cfg, is_new_month, is_new_year, doy, worst_rh_seam_gap,  &
                                    ledger, rho_air)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      integer(ik),         intent(in), optional :: doy
      real(wp),            intent(out),  optional :: worst_rh_seam_gap
      !----- The site conservation ledger (plan §10.2). Optional so a test, a probe or the C-API   !
      !      can advance a slow step without one; every entry point is a no-op when it is absent   !
      !      or inactive. `rho_air` values the canopy-air store and comes from the fast context.   !
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp),            intent(in),    optional :: rho_air

      if (present(ledger)) call slow_ledger_open(ledger, site, cfg, rho_air)

      if (present(doy)) then
         call vegetation_dynamics(site, cfg, is_new_month, is_new_year, doy, ledger)
      else
         call vegetation_dynamics(site, cfg, is_new_month, is_new_year, ledger=ledger)
      end if

      call update_patch_states(site%patch, cfg%dt_years)

      !----- Resize each patch's canopy-air control volume to the stand that now exists. Done HERE, !
      !      after growth / mortality / recruitment / fusion / disturbance have all settled, so the !
      !      fast loop sees a can_depth that is constant across every sub-step of the coming day.   !
      call refresh_canopy_depth(site, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_CANOPY)

      !----- Define it unconditionally: an intent(out) optional the callee never sets is undefined,  !
      !      and with soil carbon off there is no seam to report a gap for.  --------------------------!
      if (present(worst_rh_seam_gap)) worst_rh_seam_gap = 0.0_wp
      if (cfg%soil_carbon_on) call advance_biogeochem_dynamics(site, cfg, worst_rh_seam_gap, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, site, cfg, SLOW_PHASE_SOILC)
   end subroutine advance_slow_dynamics


   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_depth -- per-patch CAS depth = tallest cohort height + freeboard, floored.  !
   !                                                                                          !
   ! This is the ONLY writer of cas%can_depth. Before this existed the field was a hardcoded 20 m !
   ! that nothing ever assigned: the aerodynamics computed max(floor, veg_height) and the value    !
   ! was discarded, so every stand got cas_mass_capacity = 24 kg/m2 and cas_molar_capacity = 0.83 mol/m2 regardless of its    !
   ! height -- 4x too much canopy air over a 1 m regenerating gap, ~1.8x too little over a 35 m    !
   ! tropical canopy. Since cas_mass_capacity is the canopy air's heat capacity, that error goes straight into  !
   ! how fast the canopy air responds, and into the freeze-cadence stability limit on dt_fast.     !
   !                                                                                          !
   ! An empty patch falls back to the floor alone (no cohorts, no canopy to sit under).            !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_depth(site, cfg, ledger)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp)    :: h_top, depth_new, de, dw, dc, rho
      integer(ik) :: ip, i0, i1, i
      logical     :: book
      if (.not. allocated(site%patch%cas)) return
      !----- Booking the open-volume exchange needs the air density the CAS is valued at; without  !
      !      a ledger there is nothing to book it into and the resize is a plain geometry update.  !
      book = .false. ; rho = 0.0_wp
      if (present(ledger)) then
         book = ledger%active
         rho  = ledger%rho_air
      end if
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip)
         i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         h_top = 0.0_wp
         do i = i0, i1
            h_top = max(h_top, site%cohort%height(i))
         end do
         depth_new = max(cfg%aero%min_canopy_depth, h_top + cfg%aero%canopy_freeboard)
         if (book) then
            !----- Growing the canopy ENTRAINS air from just above; shrinking DETRAINS it. The     !
            !      intensive state is the invariant, so the extensive content moves, and that is   !
            !      a genuine exchange with the atmosphere -- a DECLARED boundary term, not a leak. !
            !      The hook has existed on cas_set_depth since the field was made prognostic; this !
            !      is the first caller to ask for it (plan §10.2.4). ----------------------------!
            call cas_set_depth(site%patch%cas(ip), depth_new, rho_air=rho,                         &
                               de_open=de, dw_open=dw, dc_open=dc)
            call slow_ledger_declare(ledger,                                                       &
                     energy_in=site%patch%area(ip) * de, water_in=site%patch%area(ip) * dw,        &
                     carbon_in=site%patch%area(ip) * dc * KGC_PER_UMOL_C)
         else
            call cas_set_depth(site%patch%cas(ip), depth_new)
         end if
      end do
   end subroutine refresh_canopy_depth

end module meds_slow_dynamics
