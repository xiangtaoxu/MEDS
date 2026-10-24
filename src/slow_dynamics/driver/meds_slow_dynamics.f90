! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_slow_dynamics -- the slow-tier coordinator (docs/science/order_of_processes.md      !
! §4-§5). advance_slow_dynamics runs vegetation dynamics, the patch ageing, the canopy-air !
! depth and the soil biogeochemistry in order; advance_boundary_dynamics restructures the  !
! stand at a calendar boundary; settle_tissue_water and acclimate_site_leaf_tables run     !
! once per step. It sequences the drivers and marks the slow ledger's phases               !
! (docs/science/conservation.md §3); every computation lives in the drivers it calls.      !
!==========================================================================================!
module meds_slow_dynamics
   use meds_kinds,                 only : wp, ik
   use meds_config,                only : meds_config_t
   use meds_leaf_params,           only : leaf_physiology_params_t
   use meds_time,                  only : meds_time_t, day_of_year
   use meds_polygon_state_vars, only : polygon_t, patch_polygon_area
   use meds_demography_update, only : update_patch_states
   use meds_vegetation_dynamics,   only : vegetation_dynamics, restructure_stand, shed_excess_tissue_water
   use meds_biogeochem_dynamics,   only : advance_biogeochem_dynamics
   use meds_biogeochem_types, only : litter_input_t, soilc_seam_t
   use meds_patch_state_vars,   only : cas_set_depth
   use meds_state_defaults,     only : is_unset_positive
   use meds_plant_trait_dynamics, only : acclimate_leaf_photo_table
   use meds_slow_ledger,           only : slow_ledger_t, slow_ledger_open, slow_ledger_mark,       &
                                          slow_ledger_rebase, slow_tissue_energy,                   &
                                          slow_ledger_declare, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC, &
                                          SLOW_PHASE_TISSUE_WATER, KGC_PER_UMOL_C
   implicit none
   private

   public :: advance_slow_dynamics, advance_boundary_dynamics, refresh_canopy_depth, settle_tissue_water
   public :: acclimate_site_leaf_tables

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step (order_of_processes.md §4): vegetation dynamics, patch ageing,  !
   ! the canopy-air depth, then the daily soil-carbon step. Without `step_start` (a bare   !
   ! test call) phenology is not advanced.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine advance_slow_dynamics(poly, cfg, step_start, seam, ledger, rho_air, latitude_deg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(meds_time_t),   intent(in), optional :: step_start   !< the slow step's start (its day of year)
      type(soilc_seam_t),  intent(inout), optional :: seam   !< run-worst fast/slow seam checks (soil_carbon.md §7)
      !----- The slow ledger (optional: a test may step without one) and its air density;  !
      !      the stepper passes the constant rho_air_ref (conservation.md §2). ------------!
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

      !----- Resize each canopy-air volume to the stand growth and mortality left, so the  !
      !      fast loop sees one depth all day. --------------------------------------------!
      call refresh_canopy_depth(poly, cfg, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_CANOPY)

      !----- `seam` accumulates over the run; with soil carbon off it stays as it was. ----!
      if (cfg%soil_carbon%soil_carbon_on) call advance_biogeochem_dynamics(poly, cfg, seam, ledger)
      if (present(ledger)) call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_SOILC)
   end subroutine advance_slow_dynamics

   !---------------------------------------------------------------------------------------!
   ! advance_boundary_dynamics -- the calendar boundary (order_of_processes.md §5):        !
   ! restructure the stand, then resize the canopy-air volumes; a no-op off a month        !
   ! boundary. Its marks continue the last step's ledger window.                           !
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
   ! settle_tissue_water -- before the fast loop, shed the tissue water above what the     !
   ! biomass holds and declare it: the water as a handoff to the ground, the heat as a     !
   ! thermal-mass change (conservation.md §3-§4), in a phase of its own.                   !
   !---------------------------------------------------------------------------------------!
   subroutine settle_tissue_water(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp) :: water_shed, heat0
      heat0 = 0.0_wp
      if (present(ledger)) then
         if (.not. ledger%has_first) call slow_ledger_rebase(ledger, poly, cfg)
         heat0 = slow_tissue_energy(poly)
      end if
      call shed_excess_tissue_water(poly, cfg, water_shed)
      if (present(ledger)) then
         call slow_ledger_declare(ledger, water_out = water_shed, energy_in = slow_tissue_energy(poly) - heat0)
         call slow_ledger_mark(ledger, poly, cfg, SLOW_PHASE_TISSUE_WATER)
      end if
   end subroutine settle_tissue_water

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_depth -- each patch's canopy-air depth, the tallest cohort plus the    !
   ! freeboard, floored. The only writer of can_depth, which sets the canopy air's         !
   ! capacity and the height the forcing is moved to (canopy_air_space_biophysics.md;      !
   ! forcing.md Step 2).                                                                   !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_depth(poly, cfg, ledger)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(slow_ledger_t), intent(inout), optional :: ledger
      real(wp)    :: h_top, depth_new, de, dw, dc, rho
      integer(ik) :: ip, i0, i1, i
      logical     :: book
      if (.not. allocated(poly%patch%cas)) return
      !----- Without a ledger the resize is a plain geometry update. ----------------------!
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
            !----- Air entrained or detrained is declared (conservation.md §3). -----------!
            call cas_set_depth(poly%patch%cas(ip), depth_new, rho_air=rho,                         &
                               de_open=de, dw_open=dw, dc_open=dc)
            call slow_ledger_declare(ledger,                                                       &
                     energy_in=patch_polygon_area(poly, ip) * de, water_in=patch_polygon_area(poly, ip) * dw, &
                     carbon_in=patch_polygon_area(poly, ip) * dc * KGC_PER_UMOL_C)
         else
            call cas_set_depth(poly%patch%cas(ip), depth_new)
         end if
      end do
   end subroutine refresh_canopy_depth

   !---------------------------------------------------------------------------------------!
   ! acclimate_site_leaf_tables -- advance each site's growth-temperature running mean     !
   ! (weight dt/window, on the last fast loop's daily mean) and acclimate its leaf table   !
   ! (leaf_gas_exchange.md §2 "Thermal acclimation", eq. (4)). An unset mean takes the     !
   ! first day's value, so the optimum starts from the site's own climate.                 !
   !---------------------------------------------------------------------------------------!
   subroutine acclimate_site_leaf_tables(poly, lp, step_days)
      type(polygon_t),                intent(inout) :: poly
      type(leaf_physiology_params_t), intent(in)    :: lp          !< [leaf_physiology]
      integer(ik),                    intent(in)    :: step_days
      real(wp)    :: t_day, w
      integer(ik) :: is
      do is = 1_ik, poly%site%n
         if (poly%site%pheno_tair_n(is) < 1_ik) cycle      ! no fast sub-steps ran: nothing to average
         t_day = poly%site%pheno_tair_sum(is) / real(poly%site%pheno_tair_n(is), wp)
         associate (t_growth => poly%site%t_growth_avg(is))
            if (is_unset_positive(t_growth)) then
               t_growth = t_day
            else
               w = min(1.0_wp, real(step_days, wp) / max(lp%acclim_window_days, real(step_days, wp)))
               t_growth = t_growth + w * (t_day - t_growth)
            end if
            call acclimate_leaf_photo_table(t_growth, lp%acclim_ds_vcmax_a, lp%acclim_ds_vcmax_b,            &
                                            lp%acclim_ds_jmax_a, lp%acclim_ds_jmax_b, lp%acclim_jv_a,        &
                                            lp%acclim_jv_b, poly%site%leaf_photo(is))
         end associate
      end do
   end subroutine acclimate_site_leaf_tables

end module meds_slow_dynamics
