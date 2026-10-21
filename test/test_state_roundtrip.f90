! SPDX-License-Identifier: Apache-2.0
!----- State restart round-trip: the per-cohort PLASTIC leaf traits (sla/vcmax25/rd25 and the ----!
!----- leaf turnover rate) must be written to and recovered from a state checkpoint, not reset ----!
!----- to top-of-canopy; a file holding the rate as a leaf lifespan (`llspan`) still reads, and ---!
!----- one holding the canopy films per m2 of ground (no `units`) reads them per plant. ----------!
program test_state_roundtrip
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_init,             only : init_bare_ground, add_cohort, finalize_init
   use meds_cohort_state_vars, only : hold_cohort_tissue_temps
   use meds_init_restart,               only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                     read_restart_file
   use meds_time,             only : meds_time_t
   use meds_netcdf_c,     only : NC_NOERR, nc_open_f, nc_inq_varid_f, nc_put_vara_double, nc_close, cstr
   use meds_test_support, only : banner, build_test_config, check, check_close
   use iso_c_binding,     only : c_int, c_size_t, c_double, c_char
   implicit none

   integer(c_int), parameter :: NC_WRITE = 1_c_int        ! netcdf.h: open for writing
   interface
      integer(c_int) function nc_rename_var(ncid, varid, name) bind(c, name="nc_rename_var")
         import :: c_int, c_char
         integer(c_int), value, intent(in) :: ncid, varid
         character(kind=c_char), intent(in) :: name(*)
      end function nc_rename_var
      integer(c_int) function nc_del_att(ncid, varid, name) bind(c, name="nc_del_att")
         import :: c_int, c_char
         integer(c_int), value, intent(in) :: ncid, varid
         character(kind=c_char), intent(in) :: name(*)
      end function nc_del_att
   end interface

   type(meds_config_t) :: cfg
   type(polygon_t)     :: poly, poly2
   type(meds_time_t)   :: now, rt
   logical             :: found, fast_ok
   real(wp)            :: sla_set, vc_set, rd_set, turnover_set
   integer(ik)         :: i

   call banner('state roundtrip: plastic traits')
   cfg = build_test_config()
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 0.2_wp, 20.0_wp)
   call finalize_init(poly)
   call check(poly%cohort%n >= 1_ik, 'have at least one cohort to checkpoint')

   !----- Stamp distinct, NON-top-of-canopy plastic traits on every cohort. ----------------!
   sla_set = 99.0_wp ; vc_set = 42.0_wp ; rd_set = 3.5_wp ; turnover_set = 0.7_wp
   do i = 1_ik, poly%cohort%n
      poly%cohort%sla(i)     = sla_set ; poly%cohort%vcmax25(i) = vc_set
      poly%cohort%rd25(i)    = rd_set  ; poly%cohort%leaf_turnover_rate(i) = turnover_set
      !----- Per-cohort hydraulics/temperature, stamped to distinct NON-sentinel values so a           !
      !      round-trip that fails to persist them (leaving the lazy-init 0 / LEAF_TEMP_INIT) is caught. !
      poly%cohort%leaf_water_mass(i) = 0.0123_wp + 0.001_wp * real(i, wp)
      poly%cohort%wood_water_mass(i) = 0.0456_wp + 0.001_wp * real(i, wp)
      poly%cohort%leaf_temp(i)       = 291.5_wp + real(i, wp)
      poly%cohort%wood_temp(i)       = 289.5_wp + real(i, wp)
      !----- The predawn-psi stomatal feedback. Stamped NEGATIVE and distinct so a round-trip that     !
      !      drops them is caught -- both sentinels are the opposite sign (UNSET_NEGATIVE = +1,            !
      !      DMAX_PSI_LEAF_ACCUM_RESET = -1e30), so a silent fallback cannot masquerade as a recovered value.  !
      poly%cohort%dmax_psi_leaf(i)   = -0.75_wp - 0.01_wp * real(i, wp)
      poly%cohort%dmax_psi_leaf_accum(i) = -0.40_wp - 0.01_wp * real(i, wp)
      !----- PHENOLOGY MEMORY. Set every column to a distinct, non-birth value: the point is that  !
      !      a restart must not resurrect the cohort at flush tendency 1 / shed 0 / all             !
      !      accumulators zero, which is what a restart that does not write them does. -------------!
      poly%cohort%leaf_flush_tendency(i) = 0.30_wp + 0.01_wp * real(i, wp)
      poly%cohort%leaf_shed_tendency(i)  = 0.70_wp - 0.01_wp * real(i, wp)
      poly%cohort%growing_degree_days(i) = 123.0_wp + real(i, wp)
      poly%cohort%cold_degree_days(i)    = 17.0_wp + real(i, wp)
      poly%cohort%dry_psi_sum(i)         = 5.0_wp + real(i, wp)
      poly%cohort%wet_psi_sum(i)         = 2.0_wp + real(i, wp)
      poly%cohort%light_hours_mean(i) = 10.0_wp + 0.1_wp * real(i, wp)
      !----- The canopy films [kg/plant]. ---------------------------------------------------------!
      poly%cohort%leaf_surf_water(i) = 0.020_wp + 0.001_wp * real(i, wp)
      poly%cohort%wood_surf_water(i) = 0.007_wp + 0.001_wp * real(i, wp)
      !----- The energies the stamped temperatures and waters give, the state. -------------------!
      call hold_cohort_tissue_temps(poly%cohort, i)
   end do

   !----- Stamp distinct, non-default FAST reservoir values on every patch -- a bare-ground site     !
   !      never calls init_fast_reservoirs (that is meds_main's own, fast_biophysics_on-gated step),  !
   !      so these still sit at alloc_patch_block's raw defaults (theta=0, soil_energy=0 -> 0 K) unless     !
   !      stamped here, and a restart that dropped them would restart from those defaults. ----------!
   call check(poly%patch%n >= 1_ik, 'have at least one patch to checkpoint')
   do i = 1_ik, poly%patch%n
      poly%patch%cas(i)%can_enthalpy      = 3.1e5_wp ; poly%patch%cas(i)%can_shv = 0.011_wp
      poly%patch%cas(i)%can_co2           = 415.0_wp ; poly%patch%cas(i)%can_temp = 291.0_wp
      poly%patch%soil_water(i)%theta(1:3)     = [0.28_wp, 0.31_wp, 0.34_wp]
      poly%patch%soil_energy(i)%soil_energy(1:3) = [2.1e8_wp, 2.2e8_wp, 2.3e8_wp]
      poly%patch%soil_energy(i)%soil_temp(1:3)   = [284.0_wp, 285.0_wp, 286.0_wp]
      poly%patch%soil_energy(i)%soil_fliq(1:3)   = [0.9_wp, 0.95_wp, 1.0_wp]
      poly%patch%snow(i)%swe(1)           = 4.0_wp   ; poly%patch%snow(i)%snow_energy(1) = -1.0e5_wp
      poly%patch%snow(i)%snow_depth(1)    = 0.02_wp  ; poly%patch%snow(i)%snow_temp(1)   = 270.0_wp
      poly%patch%snow(i)%snow_fliq(1)     = 0.0_wp   ; poly%patch%snow(i)%nlayer         = 1_ik
      !----- A ponded store WITH enthalpy, and the daily shed-water seam rate: the writer must put   !
      !      the pond's own enthalpy column, and the shed rate and its enthalpy, in the file. -------!
      poly%patch%soil_water(i)%w_surface      = 1.5_wp   ; poly%patch%soil_water(i)%w_surface_energy = 1.5_wp * 9.0e5_wp
      poly%patch%shed_water_rate(i)       = 3.0e-8_wp
      poly%patch%shed_water_enth(i)       = 3.0e-8_wp * 1.1e5_wp
   end do

   now = meds_time_t(year=2000, month=1, day=1, hour=0, minute=0, second=0)
   !----- The file is written on a year boundary with its restructuring still owed, as the driver's   !
   !      yearly checkpoint is (polygon_t%restructure_pending); the reader must say so. --------------!
   call write_state(poly, now, .true., .true.)
   call read_restart_file(poly2, cfg, './test_sr-S-20000101000000.nc', rt, found, polygon=1_ik, fast_found=fast_ok)
   call check(poly2%restructure_pending .and. poly2%restructure_new_year,                                 &
              'a checkpoint records the year''s restructuring as still owed')

   call check(found, 'state file read back')
   call check(fast_ok, 'fast reservoirs found in state file just written')
   call check(poly2%cohort%n == poly%cohort%n, 'cohort count preserved')
   do i = 1_ik, poly2%cohort%n
      call check_close(poly2%cohort%sla(i),     sla_set, 1.0e-9_wp, 'sla recovered from state')
      call check_close(poly2%cohort%vcmax25(i), vc_set,  1.0e-9_wp, 'vcmax25 recovered from state')
      call check_close(poly2%cohort%rd25(i),    rd_set,  1.0e-9_wp, 'rd25 recovered from state')
      call check_close(poly2%cohort%leaf_turnover_rate(i), turnover_set, 1.0e-9_wp,              &
                       'leaf_turnover_rate recovered from state')
      call check_close(poly2%cohort%leaf_water_mass(i), poly%cohort%leaf_water_mass(i), 1.0e-12_wp, &
                        'leaf_water_mass recovered from state (not re-seeded)')
      call check_close(poly2%cohort%wood_water_mass(i), poly%cohort%wood_water_mass(i), 1.0e-12_wp, &
                        'wood_water_mass recovered from state (not re-seeded)')
      call check_close(poly2%cohort%leaf_temp(i), poly%cohort%leaf_temp(i), 1.0e-9_wp, &
                        'leaf_temp recovered from state (not reset to LEAF_TEMP_INIT)')
      call check_close(poly2%cohort%wood_temp(i), poly%cohort%wood_temp(i), 1.0e-9_wp, &
                        'wood_temp recovered from state')
      !----- Without these a mid-drought restart re-seeds predawn psi from surface soil (reads as       !
      !      UNSTRESSED, silently re-opening stomata the run had closed), and a restart at any time but  !
      !      midnight discards the part of today's maximum already seen. --------------------------------!
      call check_close(poly2%cohort%dmax_psi_leaf(i), poly%cohort%dmax_psi_leaf(i), 1.0e-12_wp, &
                        'dmax_psi_leaf recovered from state (not re-seeded from soil)')
      call check_close(poly2%cohort%leaf_flush_tendency(i), poly%cohort%leaf_flush_tendency(i), 1.0e-12_wp, &
                        'leaf_flush_tendency recovered (not reset to the birth value)')
      call check_close(poly2%cohort%leaf_shed_tendency(i), poly%cohort%leaf_shed_tendency(i), 1.0e-12_wp, &
                        'leaf_shed_tendency recovered')
      call check_close(poly2%cohort%growing_degree_days(i), poly%cohort%growing_degree_days(i), 1.0e-12_wp, &
                        'growing_degree_days recovered (a season of thermal memory, not re-accumulated)')
      call check_close(poly2%cohort%cold_degree_days(i), poly%cohort%cold_degree_days(i), 1.0e-12_wp, &
                        'cold_degree_days recovered')
      call check_close(poly2%cohort%dry_psi_sum(i), poly%cohort%dry_psi_sum(i), 1.0e-12_wp, &
                        'dry_psi_sum recovered')
      call check_close(poly2%cohort%wet_psi_sum(i), poly%cohort%wet_psi_sum(i), 1.0e-12_wp, &
                        'wet_psi_sum recovered')
      call check_close(poly2%cohort%light_hours_mean(i), poly%cohort%light_hours_mean(i), 1.0e-12_wp, &
                        'light_hours_mean recovered')
      call check_close(poly2%cohort%dmax_psi_leaf_accum(i), poly%cohort%dmax_psi_leaf_accum(i), 1.0e-12_wp, &
                        'dmax_psi_leaf_accum recovered from state (partial day not discarded)')
      call check(poly2%cohort%leaf_energy(i) == poly%cohort%leaf_energy(i) .and.                  &
                 poly2%cohort%wood_energy(i) == poly%cohort%wood_energy(i),                       &
                 'the tissue energies recovered bit for bit')
      call check(poly2%cohort%leaf_temp(i) == poly%cohort%leaf_temp(i) .and.                      &
                 poly2%cohort%wood_temp(i) == poly%cohort%wood_temp(i),                           &
                 'and their temperatures, read as twins, bit for bit')
      call check(poly2%cohort%leaf_surf_water(i) == poly%cohort%leaf_surf_water(i) .and.          &
                 poly2%cohort%wood_surf_water(i) == poly%cohort%wood_surf_water(i),               &
                 'the canopy films [kg/plant] recovered bit for bit')
   end do
   do i = 1_ik, poly2%patch%n
      call check_close(poly2%patch%cas(i)%can_enthalpy, poly%patch%cas(i)%can_enthalpy, 1.0e-6_wp, &
                        'cas can_enthalpy recovered from state')
      call check_close(poly2%patch%cas(i)%can_shv,      poly%patch%cas(i)%can_shv,      1.0e-9_wp, &
                        'cas can_shv recovered from state')
      call check_close(poly2%patch%cas(i)%can_temp,     poly%patch%cas(i)%can_temp,     1.0e-9_wp, &
                        'cas can_temp recovered from state (not re-seeded to a generic value)')
      call check_close(poly2%patch%soil_water(i)%theta(2),  poly%patch%soil_water(i)%theta(2),  1.0e-9_wp, &
                        'soil theta recovered from state (not reset to 0 = bone dry)')
      call check_close(poly2%patch%soil_energy(i)%soil_temp(2), poly%patch%soil_energy(i)%soil_temp(2), 1.0e-9_wp, &
                        'soil temp recovered from state (not reset to 0 K)')
      call check_close(poly2%patch%snow(i)%swe(1),      poly%patch%snow(i)%swe(1),      1.0e-9_wp, &
                        'snow swe recovered from state')
      call check(poly2%patch%snow(i)%nlayer == poly%patch%snow(i)%nlayer, 'snow nlayer recovered from state')
      call check_close(poly2%patch%soil_water(i)%w_surface, poly%patch%soil_water(i)%w_surface, 1.0e-12_wp,       &
                       'pond mass recovered from state')
      call check_close(poly2%patch%soil_water(i)%w_surface_energy, poly%patch%soil_water(i)%w_surface_energy, 1.0e-3_wp, &
                       'pond ENTHALPY recovered from state (the writer used to put an unassigned column)')
      call check_close(poly2%patch%shed_water_rate(i), poly%patch%shed_water_rate(i), 1.0e-18_wp,          &
                       'shed_water_rate recovered from state (was not persisted)')
      call check_close(poly2%patch%shed_water_enth(i), poly%patch%shed_water_enth(i), 1.0e-15_wp,          &
                       'shed_water_enth recovered from state')
   end do

   !----- Nothing owed: a checkpoint between boundaries, or one written before the flag existed.  !
   call write_state(poly, now, .false., .false.)
   call read_restart_file(poly2, cfg, './test_sr-S-20000101000000.nc', rt, found, polygon=1_ik)
   call check(found .and. .not. poly2%restructure_pending .and. .not. poly2%restructure_new_year,          &
              'a checkpoint with nothing owed reads so')

   !----- A state file written while the trait was a leaf lifespan holds `llspan` [yr] in place of  !
   !      `leaf_turnover_rate`: the reader takes the rate as its reciprocal. The file just written is  !
   !      made into one by renaming the variable and writing a 4-yr lifespan into it. And a file      !
   !      written while the films were per m2 of ground has no `units` on them: removing it makes    !
   !      the reader take the same numbers as per m2 of ground, so per plant they are 1/nplant as   !
   !      large. ----------------------------------------------------------------------------------!
   block
      integer(c_int)              :: ncid, vid
      integer(c_size_t)           :: start(2), count(2)
      real(c_double), allocatable :: lifespan(:)
      character(kind=c_char)      :: new_name(7)
      call check(nc_open_f('./test_sr-S-20000101000000.nc', NC_WRITE, ncid) == NC_NOERR,             &
                 'the state file opens for writing')
      call check(nc_inq_varid_f(ncid, 'leaf_turnover_rate', vid) == NC_NOERR,                        &
                 'the state file holds leaf_turnover_rate')
      new_name = cstr('llspan')
      call check(nc_rename_var(ncid, vid, new_name) == NC_NOERR, 'the variable is renamed llspan')
      allocate(lifespan(poly%cohort%n)) ; lifespan = 4.0_c_double
      start = 0_c_size_t
      count = [1_c_size_t, int(poly%cohort%n, c_size_t)]                 ! the polygon's row
      call check(nc_put_vara_double(ncid, vid, start, count, lifespan) == NC_NOERR,                  &
                 'a 4-yr lifespan is written into it')
      call check(nc_inq_varid_f(ncid, 'leaf_surf_water', vid) == NC_NOERR .and.                      &
                 nc_del_att(ncid, vid, cstr('units')) == NC_NOERR, 'the leaf film loses its units')
      call check(nc_close(ncid) == NC_NOERR, 'the state file closes')
   end block
   call read_restart_file(poly2, cfg, './test_sr-S-20000101000000.nc', rt, found, polygon=1_ik)
   call check(found, 'a state file with llspan reads')
   do i = 1_ik, poly2%cohort%n
      call check_close(poly2%cohort%leaf_turnover_rate(i), 0.25_wp, 1.0e-15_wp,                     &
                       'llspan [yr] is read as the leaf turnover rate 1/llspan [1/yr]')
      call check_close(poly2%cohort%sla(i), sla_set, 1.0e-9_wp, 'the other traits read as before')
      call check_close(poly2%cohort%leaf_surf_water(i), poly%cohort%leaf_surf_water(i) / poly%cohort%nplant(i), &
                       1.0e-15_wp, 'a film without units is read per m2 of ground')
   end do

   write(*,'(a)') '   PASS'

contains

   !----- The site's state, as the one polygon of a state file (id 0, at 0 N 0 E). --------------!
   subroutine write_state(s, t, pending, new_year)
      type(polygon_t),   intent(in) :: s
      type(meds_time_t), intent(in) :: t
      logical,           intent(in) :: pending, new_year
      type(restart_file_t) :: sf
      call restart_file_open(sf, '.', 'test_sr', t, s%n_pft, int(size(s%cohort%growth_hist, 1), ik),      &
                           s%cohort%n, s%patch%n, s%site%n, polygon_id=[0_ik], polygon_lat=[0.0_wp],   &
                           polygon_lon=[0.0_wp], restructure_pending=pending, restructure_new_year=new_year)
      call restart_file_write(sf, s, t, 1_ik)
      call restart_file_close(sf)
   end subroutine write_state

end program test_state_roundtrip
