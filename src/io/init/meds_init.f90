! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init -- the start of a polygon: init_polygon runs the start [init].init_mode names   !
! (bare ground, a census or a restart file; docs/configuration.md, How a run starts) and   !
! that start's own finishing: a census stand is restructured and its plastic traits        !
! acclimated at once, a restart's traits re-acclimated when [init].reacclimate_traits asks !
! (docs/science/order_of_processes.md §1, setup step 2).                                    !
!==========================================================================================!
module meds_init
   use meds_kinds,               only : wp, ik
   use meds_time,                only : meds_time_t, time_to_string
   use meds_config,              only : meds_config_t, INIT_CENSUS, INIT_RESTART
   use meds_polygon_state_vars,  only : polygon_t
   use meds_demography_update,   only : update_overtopping_lai
   use meds_demography_cohort_fusefiss, only : max_cohort_count
   use meds_vegetation_dynamics, only : advance_plant_traits, reacclimate_plant_traits
   use meds_init_bare_ground,    only : init_bare_ground, init_bare_sites, seed_bare_ground, N_PATCH_INIT
   use meds_init_census,         only : init_from_census, restructure_census_stand
   use meds_init_restart,        only : read_restart_file
   implicit none
   private

   public :: init_polygon

contains

   !---------------------------------------------------------------------------------------!
   ! Start `poly` as cfg%init%init_mode says. A restart reads row `restart_row` of the      !
   ! restart file (0: a file without the polygon axis) and sets `now` to the file's time; a !
   ! census reads the census file's rows for the polygon's site(s); bare ground lays out    !
   ! N_PATCH_INIT patches per site and seeds them. site_cfg and site_area give each site of !
   ! a polygon of several. ok = .false. (with the reason printed, naming `label`) when the  !
   ! start's file cannot serve it.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine init_polygon(poly, cfg, label, verbose, restart_row, now, fast_found, ok, site_cfg, site_area)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      character(len=*),    intent(in)    :: label        !< the polygon's name, for messages
      logical,             intent(in)    :: verbose      !< report the start
      integer(ik),         intent(in)    :: restart_row  !< the polygon's row of the restart file (0: no axis)
      type(meds_time_t),   intent(inout) :: now          !< the run's start; a restart's time on return
      logical,             intent(out)   :: fast_found   !< a restart file held the fast reservoirs
      logical,             intent(out)   :: ok
      type(meds_config_t), intent(in), optional :: site_cfg(:)   !< each site's configuration ...
      real(wp),            intent(in), optional :: site_area(:)  !< ... and share of the polygon
      character(len=:), allocatable :: msg
      type(meds_time_t) :: restart_time
      logical           :: found

      ok = .false.
      fast_found = .false.
      select case (cfg%init%init_mode)
      case (INIT_RESTART)
         if (restart_row > 0_ik) then
            call read_restart_file(poly, cfg, trim(cfg%init%restart_file), restart_time, found,          &
                                   fast_found=fast_found, polygon=restart_row, site_cfg=site_cfg, why=msg)
         else
            call read_restart_file(poly, cfg, trim(cfg%init%restart_file), restart_time, found,          &
                                   fast_found=fast_found, site_cfg=site_cfg, why=msg)
         end if
         if (.not. found) then
            write(*,'(7a)') ' init  : restart (mode 2): cannot read ', trim(label),                       &
                            ' from the state file ', trim(cfg%init%restart_file), ', which ', msg
            return
         end if
         now = restart_time
         if (verbose) call report_restart(cfg, now, fast_found, poly%restructure_pending)
      case (INIT_CENSUS)
         if (present(site_cfg)) then
            call init_from_census(poly, cfg, trim(cfg%init%census_file), found, why=msg,                   &
                                  site_ids=cfg%sites%census_site_id, site_cfg=site_cfg, site_area=site_area)
         else
            call init_from_census(poly, cfg, trim(cfg%init%census_file), found, why=msg)
         end if
         if (.not. found) then
            write(*,'(4a)') ' init  : census (mode 1): the census file ', trim(cfg%init%census_file), ' ', msg
            return
         end if
         if (verbose) then
            write(*,'(2a)') ' init  : census (mode 1) ', trim(cfg%init%census_file)
            write(*,'(a,i0,a,i0,a)') '         read ', poly%patch%n, ' patches, ', poly%cohort%n, ' cohorts'
         end if
         !----- A census stand is restructured by the slow step's own operators before the first  !
         !      step: it arrives with a cohort per measured size and a patch per plot cell.  -------!
         call restructure_census_stand(poly, cfg)
         if (verbose) write(*,'(a,i0,a,i0,a,i0,a)') '         restructured to ', poly%patch%n,          &
                                                     ' patches, ', poly%cohort%n, ' cohorts (at most ',   &
                                                     max_cohort_count(poly), ' in a patch)'
      case default
         if (present(site_cfg)) then
            call init_bare_sites(poly, cfg, spread(N_PATCH_INIT, 1, size(site_cfg)), site_cfg, site_area)
         else
            call init_bare_ground(poly, cfg, N_PATCH_INIT)
         end if
         call seed_bare_ground(poly, cfg)
         if (verbose) write(*,'(a)') ' init  : bare ground (mode 0)'
      end select

      !----- Census cohorts sit in an established stand but carry no trait history, so with      !
      !      plasticity on their leaf traits acclimate to their light at once. Bare ground starts !
      !      at top of canopy; a restart keeps the traits it read, unless                         !
      !      [init].reacclimate_traits asks for this run's PFT file instead.                      !
      if (cfg%trait_plasticity_on .and. cfg%init%init_mode == INIT_CENSUS) then
         call update_overtopping_lai(poly)
         call advance_plant_traits(poly, cfg, cfg%dt_years, instantaneous=.true.)
      end if
      if (cfg%init%reacclimate_traits .and. cfg%init%init_mode == INIT_RESTART) then
         call reacclimate_plant_traits(poly, cfg)
         if (verbose) write(*,'(a)') ' init  : plastic traits re-acclimated to this run''s PFT file'
      end if
      ok = .true.
   end subroutine init_polygon

   !----- What a restart read: its file and time, whether the fast reservoirs came with it, and a   !
   !      restructuring the stand still owes. --------------------------------------------------------!
   subroutine report_restart(cfg, now, fast_found, pending)
      type(meds_config_t), intent(in) :: cfg
      type(meds_time_t),   intent(in) :: now
      logical,             intent(in) :: fast_found, pending
      write(*,'(4a)') ' init  : restart (mode 2) ', trim(cfg%init%restart_file), ' @ ', time_to_string(now)
      if (fast_found) then
         write(*,'(a)') '         (CAS/soil/snow fast reservoirs restored -- exact restart)'
      else
         write(*,'(a)') '         (a file without the fast reservoirs -- re-seeding)'
      end if
      if (pending) write(*,'(a)') '         (the boundary''s restructuring is pending: it runs before the first step)'
   end subroutine report_restart

end module meds_init
