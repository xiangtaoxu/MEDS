! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init_restart -- the state (restart) file: every polygon's prognostic state at an    !
! instant, in a netCDF-4 checkpoint a resumed run continues from exactly                   !
! (docs/science/order_of_processes.md §6, The restart file). restart_file_open, _write and !
! _close write a file a polygon's row at a time; read_restart_file reads a row, or the     !
! whole of a file without the polygon axis (one site's); restart_file_slots checks a file  !
! before a run starts from it. Each level's variables are defined, written and read in one !
! place, in its state module's section order. Diagnostic output is meds_output_*.          !
!==========================================================================================!
module meds_init_restart
   use iso_c_binding,       only : c_int, c_size_t
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : n_soil_layer_max, n_snow_layer_max
   use meds_state_defaults, only : is_unset_positive
   use meds_therm_lib,      only : internal_energy_liquid
   use meds_config,         only : meds_config_t, growth_window_steps
   use meds_time,           only : meds_time_t, time_to_stamp, time_to_decimal_year
   use meds_netcdf_c,       only : nc_check, nc_close, nc_create_f, nc_def_dim_f, nc_enddef, nc_get_att_text_f,   &
                                   nc_get_vara_double, nc_get_vara_int, nc_inq_dimlen, nc_inq_dimlen_f,           &
                                   nc_inq_varid_f, nc_open_f, nc_put_vara_double, nc_put_vara_int, NC_CLOBBER,    &
                                   NC_DOUBLE, NC_GLOBAL, NC_INT, NC_NETCDF4, NC_NOERR, NC_NOWRITE
   use meds_netcdf_write,   only : define_var, put_global_text, define_polygon_axis, write_polygon_axis,             &
                                   RESTART_DEFLATE
   use meds_polygon_state_vars,         only : polygon_t, alloc_polygon, rebuild_patch_map, rebuild_cohort_map
   use meds_site_state_vars,            only : make_site_surface
   use meds_patch_state_vars,           only : patch_physical_state
   use meds_cohort_state_vars,          only : set_cohort_size, set_cohort_size_from_carbon, set_cohort_tissue_water, &
                                               set_cohort_tissue_temps, hold_cohort_tissue_temps
   use meds_demography_update,          only : update_overtopping_lai
   use meds_demography_cohort_fusefiss, only : sort_cohorts
   implicit none
   private

   public :: restart_file_t, restart_file_open, restart_file_write, restart_file_close
   public :: read_restart_file, restart_file_slots

   character(len=*), parameter :: STATE_TITLE = 'MEDS instantaneous state (restart) file'
   integer,          parameter :: N_META_INT = 11, N_META_REAL = 2
   character(len=*), parameter :: FILM_UNITS  = 'kg/plant'   !< the canopy films' units attribute

   !----- An open state file, its polygon axis ahead of every other. --------------------------------!
   type :: restart_file_t
      integer(c_int)     :: ncid   = -1_c_int
      character(len=512) :: path   = ''
   end type restart_file_t

   !----- The dimensions of a file being defined. --------------------------------------------------!
   type :: file_dims_t
      integer(c_int) :: ncid = -1_c_int
      integer(c_int) :: polygon, cohort, patch, site, pft, meta_int, meta_real, soil_layer, snow_layer
      integer(c_int) :: growth_window
   end type file_dims_t

   !----- One polygon's row of an open file: row `row` (1-based) of a file with the polygon axis, or !
   !      the whole of a file without it (one site's). -----------------------------------------------!
   type :: file_row_t
      integer(c_int) :: ncid     = -1_c_int
      integer(ik)    :: row      = 1_ik
      logical        :: has_axis = .true.
   end type file_row_t

   !----- What a file held of the state that older files lack; section 4 derives the rest. ----------!
   type :: file_held_t
      logical :: pools           = .false.   !< the carbon pools (else each cohort is placed from its dbh)
      logical :: overtopping_lai = .false.
      logical :: growth_rings    = .false.   !< the growth rings, of this run's window length
      logical :: films_per_plant = .true.    !< the canopy films in kg/plant (else per m2 of ground)
      logical :: tissue_water    = .true.
      logical :: tissue_temps    = .true.
      logical :: tissue_energy   = .true.
      logical :: fast            = .false.   !< the fast reservoirs (else the caller re-seeds them)
      logical :: pond_energy     = .false.
      logical :: shed_energy     = .false.
   end type file_held_t

contains

   !=======================================================================================!
   !  0. THE FILE: open, write, close and read, the checks before a run starts from one,   !
   !     and one variable defined, written or read at a polygon's row                      !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Create <dir>/<prefix>-S-YYYYMMDDHHMMSS.nc and define every dimension and variable.    !
   ! The cohort and patch axes are n_cohort and n_patch long (at least 1), the largest of  !
   ! the polygons'. A `polygon` axis leads every variable, `polygon_id` names its rows,    !
   ! and the variables are compressed, since most polygons fill only part of the padded    !
   ! axes. Why each group persists: docs/science/order_of_processes.md §6, The restart     !
   ! file.                                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine restart_file_open(sf, dir, prefix, now, n_pft, n_window, n_cohort, n_patch, n_site,     &
                              polygon_id, polygon_lat, polygon_lon, restructure_pending, restructure_new_year)
      type(restart_file_t),  intent(out) :: sf
      character(len=*),      intent(in)  :: dir, prefix
      type(meds_time_t),     intent(in)  :: now
      integer(ik),           intent(in)  :: n_pft, n_window, n_cohort, n_patch, n_site
      integer(ik),           intent(in)  :: polygon_id(:)                  !< the polygons' ids ...
      real(wp),              intent(in)  :: polygon_lat(:), polygon_lon(:)  !< ... and locations [deg]
      logical,     optional, intent(in)  :: restructure_pending, restructure_new_year
      type(file_dims_t) :: d
      character(len=5)  :: pending
      integer(c_int)    :: ncid

      !----- The stamp is the simulated date: <prefix>-S-YYYYMMDDHHMMSS.nc. -----------------!
      sf%path = trim(dir)//'/'//trim(prefix)//'-S-'//time_to_stamp(now)//'.nc'
      call nc_check(nc_create_f(trim(sf%path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'state nc_create')
      sf%ncid = ncid
      d%ncid = ncid
      call define_polygon_axis(ncid, size(polygon_id, kind=ik), d%polygon)
      call nc_check(nc_def_dim_f(ncid, 'cohort', int(max(n_cohort,1_ik), c_size_t), d%cohort), 'state dim cohort')
      call nc_check(nc_def_dim_f(ncid, 'patch',  int(max(n_patch,1_ik), c_size_t), d%patch),   'state dim patch')
      call nc_check(nc_def_dim_f(ncid, 'site',   int(max(n_site,1_ik), c_size_t),  d%site),    'state dim site')
      call nc_check(nc_def_dim_f(ncid, 'pft',    int(n_pft, c_size_t),             d%pft),     'state dim pft')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_int',  int(N_META_INT, c_size_t), d%meta_int), 'state dim mi')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_real', int(N_META_REAL, c_size_t), d%meta_real), 'state dim mr')
      call nc_check(nc_def_dim_f(ncid, 'soil_layer', int(n_soil_layer_max, c_size_t), d%soil_layer), 'state dim soill')
      call nc_check(nc_def_dim_f(ncid, 'snow_layer', int(n_snow_layer_max, c_size_t), d%snow_layer), 'state dim snowl')
      call nc_check(nc_def_dim_f(ncid, 'growth_window', int(max(n_window,1_ik), c_size_t), d%growth_window),     &
                    'state dim gwin')
      call define_polygon_vars(d)
      call define_cohort_vars(d)
      call define_patch_vars(d)
      call put_global_text(ncid, 'title', STATE_TITLE)
      pending = 'none'
      if (present(restructure_pending)) then
         if (restructure_pending) pending = 'month'
         if (restructure_pending .and. present(restructure_new_year)) then
            if (restructure_new_year) pending = 'year'
         end if
      end if
      call put_global_text(ncid, 'restructure_pending', trim(pending))
      call nc_check(nc_enddef(ncid), 'state enddef')
      call write_polygon_axis(ncid, polygon_id, polygon_lat, polygon_lon)
   end subroutine restart_file_open

   !---------------------------------------------------------------------------------------!
   ! Write one polygon's state into its row (1-based) of an open file.                     !
   !---------------------------------------------------------------------------------------!
   subroutine restart_file_write(sf, poly, now, polygon)
      type(restart_file_t), intent(in) :: sf
      type(polygon_t),      intent(in) :: poly
      type(meds_time_t),    intent(in) :: now
      integer(ik),          intent(in) :: polygon
      type(file_row_t) :: r
      r = file_row_t(ncid=sf%ncid, row=polygon, has_axis=.true.)
      call put_polygon_vars(r, poly, now)
      if (poly%cohort%n > 0_ik) call put_cohort_vars(r, poly)
      if (poly%patch%n > 0_ik)  call put_patch_vars(r, poly)
   end subroutine restart_file_write

   subroutine restart_file_close(sf)
      type(restart_file_t), intent(inout) :: sf
      call nc_check(nc_close(sf%ncid), 'state nc_close')
      sf%ncid = -1_c_int
      write(*,'(2a)') ' state : ', trim(sf%path)
   end subroutine restart_file_close

   !---------------------------------------------------------------------------------------!
   ! Reconstruct a polygon from a state file written by restart_file_write: row `polygon`  !
   ! (1-based) of a file with the polygon axis, or the whole of a file without it. Each    !
   ! site's surface is made from `site_cfg` (else `cfg`) and the geometry from the carbon  !
   ! pools; the maps and sort order are rebuilt, and the restructuring owed is the file's  !
   ! restructure_pending (docs/science/order_of_processes.md §5). found = .false. if the   !
   ! file cannot be opened or has another layout (restart_file_slots says why); a PFT,     !
   ! site or soil-layer count unlike the run's stops the run.                              !
   !---------------------------------------------------------------------------------------!
   subroutine read_restart_file(poly, cfg, path, restart_time, found, fast_found, polygon, site_cfg)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      type(meds_time_t),   intent(out) :: restart_time
      logical,             intent(out) :: found
      !----- .true. when the file's fast reservoirs were read; otherwise the caller re-seeds them   !
      !      (init_fast_reservoirs). -----------------------------------------------------------------!
      logical, optional,   intent(out) :: fast_found
      integer(ik), optional, intent(in) :: polygon    !< its row (1-based); absent for a file without the axis
      type(meds_config_t), optional, intent(in) :: site_cfg(:)   !< each site's configuration
      type(file_row_t)  :: r
      type(file_held_t) :: held
      integer(c_int)    :: ncid, st

      found = .false.
      restart_time = meds_time_t()
      if (present(fast_found)) fast_found = .false.
      if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) return
      r = file_row_t(ncid=ncid, has_axis=present(polygon))
      if (present(polygon)) r%row = polygon
      !----- The polygon axis, where a file has it, leads every variable: reading the other kind   !
      !      of file would take the wrong elements. --------------------------------------------------!
      if (has_polygon_axis(ncid) .neqv. r%has_axis) then
         st = nc_close(ncid)
         return
      end if
      call get_polygon_vars(r, poly, cfg, restart_time, site_cfg)
      if (poly%cohort%n > 0_ik) call get_cohort_vars(r, poly, cfg, held)
      if (poly%patch%n > 0_ik)  call get_patch_vars(r, poly, held)
      call nc_check(nc_close(ncid), 'state nc_close (read)')
      call derive_after_read(poly, cfg, held)
      if (present(fast_found)) fast_found = held%fast
      found = .true.
   end subroutine read_restart_file

   !----- Each polygon's row of the restart file at `path`, by its id: 0 for every polygon of a file  !
   !      without the polygon axis (one site's). why = '' when the file can start a run of these     !
   !      polygons, else the reason it cannot (restart_file_refusal, or a polygon it does not hold). -!
   subroutine restart_file_slots(path, ids, slot, why)
      character(len=*),              intent(in)  :: path
      integer(ik),                   intent(in)  :: ids(:)
      integer(ik), allocatable,      intent(out) :: slot(:)
      character(len=:), allocatable, intent(out) :: why
      integer(ik), allocatable :: file_ids(:)
      integer(ik)       :: p
      logical           :: found
      character(len=24) :: id_str
      allocate(slot(size(ids)), source=0_ik)
      why = restart_file_refusal(path, size(ids, kind=ik))
      if (len(why) > 0) return
      call restart_file_polygons(path, file_ids, found)
      if (.not. found) return
      do p = 1_ik, size(ids, kind=ik)
         slot(p) = int(findloc(file_ids, ids(p), dim=1), ik)
         if (slot(p) == 0_ik) then
            write(id_str, '(i0)') ids(p)
            why = 'does not hold polygon '//trim(id_str)//', and must hold every polygon of the run'
            return
         end if
      end do
   end subroutine restart_file_slots

   !----- A state file's polygon ids, in its rows' order; found = .false. for a file without the  !
   !      polygon axis (one site's) or one that cannot be opened. -------------------------------!
   subroutine restart_file_polygons(path, polygon_id, found)
      character(len=*),         intent(in)  :: path
      integer(ik), allocatable, intent(out) :: polygon_id(:)
      logical,                  intent(out) :: found
      integer(c_int)    :: ncid, vid, st
      integer(c_size_t) :: n, start1(1), count1(1)
      found = .false.
      if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) return
      st = nc_inq_varid_f(ncid, 'polygon_id', vid)
      if (st == NC_NOERR) st = nc_inq_dimlen_f(ncid, 'polygon', n)
      if (st == NC_NOERR) then
         allocate(polygon_id(n))
         start1 = 0_c_size_t
         count1 = n
         st = nc_get_vara_int(ncid, vid, start1, count1, polygon_id)
      end if
      found = st == NC_NOERR
      st = nc_close(ncid)
   end subroutine restart_file_polygons

   !---------------------------------------------------------------------------------------!
   ! Why the state file at `path` cannot start a run of `n_polygon` polygons, or '' when   !
   ! it can: it cannot be opened, or it holds one site without a `polygon` axis and the    !
   ! run has more than one polygon.                                                        !
   !---------------------------------------------------------------------------------------!
   function restart_file_refusal(path, n_polygon) result(why)
      character(len=*), intent(in)  :: path
      integer(ik),      intent(in)  :: n_polygon
      character(len=:), allocatable :: why
      integer(c_int) :: ncid, st
      logical        :: has_axis
      why = ''
      if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) then
         why = 'cannot be opened'
         return
      end if
      has_axis = has_polygon_axis(ncid)
      st = nc_close(ncid)
      if (n_polygon > 1_ik .and. .not. has_axis)                                                     &
         why = 'holds one site without a polygon axis; a run of several polygons needs a row for each'
   end function restart_file_refusal

   !----- Does the open state file have a `polygon` axis ahead of every other? -------------------!
   logical function has_polygon_axis(ncid)
      integer(c_int), intent(in) :: ncid
      integer(c_size_t) :: n
      has_polygon_axis = nc_inq_dimlen_f(ncid, 'polygon', n) == NC_NOERR
   end function has_polygon_axis

   !----- A per-polygon variable: one value per row, uncompressed, with its long_name. ------------!
   subroutine def_polygon_var(d, name, xtype, lname)
      type(file_dims_t), intent(in) :: d
      character(len=*),  intent(in) :: name, lname
      integer(c_int),    intent(in) :: xtype
      integer(c_int) :: v
      call define_var(d%ncid, name, xtype, [d%polygon], v, long_name=lname)
   end subroutine def_polygon_var

   !----- A state variable on its own dimensions, behind the polygon axis, with its long_name and,  !
   !      where the reader needs it, its units; one polygon's row per chunk, since each polygon   !
   !      writes its own. ----------------------------------------------------------------------!
   subroutine def_var(d, name, xtype, dimids, lname, units)
      type(file_dims_t), intent(in) :: d
      character(len=*),  intent(in) :: name, lname
      integer(c_int),    intent(in) :: xtype, dimids(:)
      character(len=*),  intent(in), optional :: units
      integer(c_int)    :: v, dims(size(dimids) + 1)
      integer(c_size_t) :: chunk(size(dimids) + 1)
      integer           :: k
      dims(1) = d%polygon
      dims(2:) = dimids
      chunk(1) = 1_c_size_t
      do k = 2, size(dims)
         call nc_check(nc_inq_dimlen(d%ncid, dims(k), chunk(k)), 'state dimlen '//name)
         chunk(k) = max(chunk(k), 1_c_size_t)
      end do
      call define_var(d%ncid, name, xtype, dims, v, units=units, long_name=lname, chunk=chunk,     &
                      deflate=RESTART_DEFLATE)
   end subroutine def_var

   !----- A per-polygon value (n_cohort, n_patch, n_site, n_window), at the polygon's row. --------!
   subroutine put_polygon_int(r, name, n)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: n
      integer(c_size_t) :: st(1), ct(1)
      integer(ik)       :: x(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      x = n
      call nc_check(nc_put_vara_int(r%ncid, var_id(r%ncid, name), st, ct, x), 'put '//name)
   end subroutine put_polygon_int

   integer(ik) function get_polygon_int(r, name) result(n)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(c_size_t) :: st(1), ct(1)
      integer(ik)       :: x(1)
      st = int(r%row - 1_ik, c_size_t)
      ct = 1_c_size_t
      call nc_check(nc_get_vara_int(r%ncid, var_id(r%ncid, name), st, ct, x), 'get '//name)
      n = x(1)
   end function get_polygon_int

   !----- One vector into, or out of, the variable `name` at the polygon's row: its first entries  !
   !      (a per-cohort, per-patch, per-site or meta vector), or with `item`, the layers of that    !
   !      cohort or patch (1-based) of a variable with a second axis. ---------------------------!
   subroutine put_dbl(r, name, x, item)
      type(file_row_t),      intent(in) :: r
      character(len=*),      intent(in) :: name
      real(wp),              intent(in) :: x(:)
      integer(ik), optional, intent(in) :: item
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), item, st, ct, nd)
      call nc_check(nc_put_vara_double(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'put '//name)
   end subroutine put_dbl

   subroutine put_int(r, name, x)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: x(:)
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), st=st, ct=ct, nd=nd)
      call nc_check(nc_put_vara_int(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'put '//name)
   end subroutine put_int

   subroutine get_dbl(r, name, x, item)
      type(file_row_t),      intent(in)  :: r
      character(len=*),      intent(in)  :: name
      real(wp),              intent(out) :: x(:)
      integer(ik), optional, intent(in)  :: item
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), item, st, ct, nd)
      call nc_check(nc_get_vara_double(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'get '//name)
   end subroutine get_dbl

   subroutine get_int(r, name, x)
      type(file_row_t), intent(in)  :: r
      character(len=*), intent(in)  :: name
      integer(ik),      intent(out) :: x(:)
      integer(c_size_t) :: st(3), ct(3)
      integer           :: nd
      call where_in_file(r, size(x, kind=ik), st=st, ct=ct, nd=nd)
      call nc_check(nc_get_vara_int(r%ncid, var_id(r%ncid, name), st(1:nd), ct(1:nd), x), 'get '//name)
   end subroutine get_int

   !----- get_dbl where the file holds the variable; `x` keeps what it had where it does not. -----!
   subroutine get_dbl_if_held(r, name, x)
      type(file_row_t), intent(in)    :: r
      character(len=*), intent(in)    :: name
      real(wp),         intent(inout) :: x(:)
      if (file_has(r, name)) call get_dbl(r, name, x)
   end subroutine get_dbl_if_held

   logical function file_has(r, name)
      type(file_row_t), intent(in) :: r
      character(len=*), intent(in) :: name
      integer(c_int) :: v
      file_has = nc_inq_varid_f(r%ncid, name, v) == NC_NOERR
   end function file_has

   integer(c_int) function var_id(ncid, name)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name
      call nc_check(nc_inq_varid_f(ncid, name, var_id), 'state inq '//name)
   end function var_id

   !----- Where a polygon's n values sit: its row, when the file has the polygon axis; then entry  !
   !      `item` (1-based) of the cohort or patch axis, when given; then the first n entries. -----!
   pure subroutine where_in_file(r, n, item, st, ct, nd)
      type(file_row_t),      intent(in)  :: r
      integer(ik),           intent(in)  :: n
      integer(ik), optional, intent(in)  :: item
      integer(c_size_t),     intent(out) :: st(3), ct(3)
      integer,               intent(out) :: nd
      nd = 0
      if (r%has_axis) then
         nd = 1
         st(1) = int(r%row - 1_ik, c_size_t)
         ct(1) = 1_c_size_t
      end if
      if (present(item)) then
         nd = nd + 1
         st(nd) = int(item - 1_ik, c_size_t)
         ct(nd) = 1_c_size_t
      end if
      nd = nd + 1
      st(nd) = 0_c_size_t
      ct(nd) = int(n, c_size_t)
   end subroutine where_in_file

   !=======================================================================================!
   !  1. THE POLYGON AND ITS SITES, in meds_polygon_state_vars' and meds_site_state_vars'  !
   !     section order                                                                     !
   !=======================================================================================!

   subroutine define_polygon_vars(d)
      type(file_dims_t), intent(in) :: d
      !----- 1. Identity and location: the polygon axis (restart_file_open). ---------------------!
      !----- 2. Dimensions and id counters, with the file's time in meta_int. -----------------!
      call def_polygon_var(d, 'n_cohort',   NC_INT,    'the polygon''s live cohorts (the rest is fill)')
      call def_polygon_var(d, 'n_patch',    NC_INT,    'the polygon''s live patches (the rest is fill)')
      call def_polygon_var(d, 'n_site',     NC_INT,    'the polygon''s sites (the rest is fill)')
      call def_polygon_var(d, 'n_window',   NC_INT,    'the polygon''s growth-window length')
      call def_var(d, 'meta_int',  NC_INT,    [d%meta_int],                                  &
                   '[n_cohort,n_patch,n_pft,next_cohort_id,next_patch_id,year,month,day,hour,minute,second]')
      call def_var(d, 'meta_real', NC_DOUBLE, [d%meta_real], '[area of the sites (1), decimal_year]')
      !----- 3. Calendar: the cohorts' ring-buffer slot (the restructuring owed is a global      !
      !      attribute). ----------------------------------------------------------------------!
      call def_var(d, 'growth_hist_pos', NC_INT, [d%meta_real], 'ring-buffer write slot (element 1)')
      !----- The sites: 1. each one's area; 3. the surface it was made with (its layer count is   !
      !      checked on read); 4. thermal acclimation's running mean                              !
      !      (docs/science/leaf_gas_exchange.md §2). ------------------------------------------!
      call def_var(d, 'site_area',       NC_DOUBLE, [d%site], 'site area fraction of the polygon')
      call def_var(d, 'site_n_soil',     NC_INT,    [d%site], 'the site''s active soil layers')
      call def_var(d, 'site_soil_depth', NC_DOUBLE, [d%site], 'the site''s soil column depth [m]')
      call def_var(d, 'site_elevation',  NC_DOUBLE, [d%site], 'the site''s elevation [m]')
      call def_var(d, 'site_slope',      NC_DOUBLE, [d%site], 'the site''s slope [degrees]')
      call def_var(d, 'site_aspect',     NC_DOUBLE, [d%site], 'the site''s aspect [degrees from north]')
      call def_var(d, 't_growth_avg',    NC_DOUBLE, [d%site], 'the site''s running-mean growth temperature [K]')
   end subroutine define_polygon_vars

   subroutine put_polygon_vars(r, poly, now)
      type(file_row_t),  intent(in) :: r
      type(polygon_t),   intent(in) :: poly
      type(meds_time_t), intent(in) :: now
      integer(ik) :: ns, is, nsoil(poly%site%n)
      real(wp)    :: depth(poly%site%n), elev(poly%site%n), slope(poly%site%n), aspect(poly%site%n)
      ns = poly%site%n
      !----- 2. Dimensions and id counters, and the time. ------------------------------------!
      call put_polygon_int(r, 'n_cohort', poly%cohort%n)
      call put_polygon_int(r, 'n_patch',  poly%patch%n)
      call put_polygon_int(r, 'n_site',   ns)
      call put_polygon_int(r, 'n_window', int(size(poly%cohort%growth_hist, 1), ik))
      call put_int(r, 'meta_int', [poly%cohort%n, poly%patch%n, poly%n_pft, poly%next_cohort_id,              &
                                   poly%next_patch_id, now%year, now%month, now%day, now%hour, now%minute,     &
                                   now%second])
      call put_dbl(r, 'meta_real', [sum(poly%site%area(1:ns)), time_to_decimal_year(now)])
      !----- 3. Calendar. --------------------------------------------------------------------!
      call put_int(r, 'growth_hist_pos', [poly%growth_hist_pos])
      !----- The sites, the surface's numbers through contiguous buffers. --------------------!
      do is = 1_ik, ns
         nsoil(is)  = poly%site%surface(is)%soil%n_active
         depth(is)  = poly%site%surface(is)%soil_column%depth
         elev(is)   = poly%site%surface(is)%elevation_m
         slope(is)  = poly%site%surface(is)%slope_deg
         aspect(is) = poly%site%surface(is)%aspect_deg
      end do
      call put_dbl(r, 'site_area',       poly%site%area(1:ns))
      call put_int(r, 'site_n_soil',     nsoil)
      call put_dbl(r, 'site_soil_depth', depth)
      call put_dbl(r, 'site_elevation',  elev)
      call put_dbl(r, 'site_slope',      slope)
      call put_dbl(r, 'site_aspect',     aspect)
      call put_dbl(r, 't_growth_avg',    poly%site%t_growth_avg(1:ns))
   end subroutine put_polygon_vars

   !----- The polygon, allocated to the file's counts, and its sites, each surface made from     !
   !      site_cfg (else cfg). A file without the site axis holds one site. ---------------------!
   subroutine get_polygon_vars(r, poly, cfg, restart_time, site_cfg)
      type(file_row_t),    intent(in)  :: r
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      type(meds_time_t),   intent(out) :: restart_time
      type(meds_config_t), intent(in), optional :: site_cfg(:)
      character(len=16)  :: pending
      character(len=160) :: msg
      integer(ik) :: meta_i(N_META_INT), ns, is
      real(wp)    :: meta_r(N_META_REAL)

      !----- 2. Dimensions and id counters, and the file's time. -----------------------------!
      call get_int(r, 'meta_int',  meta_i)
      call get_dbl(r, 'meta_real', meta_r)
      if (meta_i(3) /= cfg%pft%n) error stop 'read_restart_file: PFT count in state file /= config PFT count'
      ns = 1_ik
      if (file_has(r, 'n_site')) ns = get_polygon_int(r, 'n_site')
      if (present(site_cfg)) then
         if (size(site_cfg, kind=ik) /= ns) then
            write(msg, '(a,i0,a,i0,a)') 'read_restart_file: the state file holds ', ns,                     &
                                        ' sites and the configuration ', size(site_cfg), ' for the polygon'
            error stop trim(msg)
         end if
      end if
      call alloc_polygon(poly, cfg%pft%n, max(meta_i(1), 1_ik), max(meta_i(2), 1_ik), growth_window_steps(cfg), &
                         cfg%pft, n_site=ns)
      poly%cohort%n = meta_i(1)
      poly%patch%n = meta_i(2)
      poly%next_cohort_id = meta_i(4)
      poly%next_patch_id = meta_i(5)
      restart_time = meds_time_t(year=meta_i(6), month=meta_i(7), day=meta_i(8),               &
                                 hour=meta_i(9), minute=meta_i(10), second=meta_i(11))
      !----- 3. Calendar: the restructuring owed (the ring's slot comes with the rings, section 2). !
      pending = 'none'
      if (nc_get_att_text_f(r%ncid, NC_GLOBAL, 'restructure_pending', pending) /= NC_NOERR) pending = 'none'
      poly%restructure_pending  = trim(pending) == 'month' .or. trim(pending) == 'year'
      poly%restructure_new_year = trim(pending) == 'year'
      !----- The sites: 3. each surface; 1. each area, and each site's soil state is per layer of   !
      !      the column it was written on; a file without them has one site, every patch in it. ---!
      do is = 1_ik, ns
         if (present(site_cfg)) then
            call make_site_surface(poly%site%surface(is), site_cfg(is))
         else
            call make_site_surface(poly%site%surface(is), cfg)
         end if
      end do
      if (file_has(r, 'site_area')) then
         call get_dbl(r, 'site_area', poly%site%area(1:ns))
         block
            integer(ik) :: nsoil(ns)
            call get_int(r, 'site_n_soil', nsoil)
            do is = 1_ik, ns
               if (nsoil(is) /= poly%site%surface(is)%soil%n_active) then
                  write(msg, '(a,i0,a,i0,a,i0,a)') 'read_restart_file: site ', is, ' has ', nsoil(is),       &
                        ' soil layers in the state file and ', poly%site%surface(is)%soil%n_active,     &
                        ' in the configuration'
                  error stop trim(msg)
               end if
            end do
         end block
      else
         poly%site%area(1) = meta_r(1)
      end if
      !----- 4. Growth temperature: a file without it keeps alloc_polygon's unseeded sentinel. --!
      call get_dbl_if_held(r, 't_growth_avg', poly%site%t_growth_avg(1:ns))
   end subroutine get_polygon_vars

   !=======================================================================================!
   !  2. THE COHORTS, in meds_cohort_state_vars' section order                             !
   !=======================================================================================!

   subroutine define_cohort_vars(d)
      type(file_dims_t), intent(in) :: d
      !----- 1. Identity and density. --------------------------------------------------------!
      call def_var(d, 'pft',              NC_INT,    [d%cohort], 'plant functional type index')
      call def_var(d, 'global_cohort_id', NC_INT,    [d%cohort], 'persistent cohort id')
      call def_var(d, 'owner_patch',      NC_INT,    [d%cohort], 'owning patch index (1-based)')
      call def_var(d, 'nplant',           NC_DOUBLE, [d%cohort], 'plant number density [plant/m2]')
      !----- 2. The carbon pools, as the run holds them. -------------------------------------!
      call def_var(d, 'leaf_carbon',          NC_DOUBLE, [d%cohort], 'leaf carbon per plant [kgC]')
      call def_var(d, 'fineroot_carbon',      NC_DOUBLE, [d%cohort], 'fine-root carbon per plant [kgC]')
      call def_var(d, 'wood_carbon',          NC_DOUBLE, [d%cohort], 'wood carbon per plant [kgC]')
      call def_var(d, 'nonstructural_carbon', NC_DOUBLE, [d%cohort], 'nonstructural carbon per plant [kgC]')
      !----- 3. Size: the reader derives it from the pools, so dbh, height, basal_area, agb and   !
      !      leaf_area are for readers of the file; overtopping_lai is the cohort's light. -------!
      call def_var(d, 'dbh',             NC_DOUBLE, [d%cohort], 'diameter at breast height [cm]')
      call def_var(d, 'height',          NC_DOUBLE, [d%cohort], 'height [m]')
      call def_var(d, 'basal_area',      NC_DOUBLE, [d%cohort], 'basal area per plant [cm2]')
      call def_var(d, 'agb',             NC_DOUBLE, [d%cohort], 'aboveground biomass per plant [kgC]')
      call def_var(d, 'leaf_area',       NC_DOUBLE, [d%cohort], 'leaf area per plant [m2]')
      call def_var(d, 'overtopping_lai', NC_DOUBLE, [d%cohort], 'LAI of the taller cohorts in the patch [m2/m2]')
      !----- 4. Water: the tissue water; the canopy interception films ([fast].canopy_water_on),  !
      !      per plant and marked so; the stomatal feedback's predawn potential                    !
      !      (docs/science/leaf_gas_exchange.md §4.2). ----------------------------------------!
      call def_var(d, 'leaf_water_mass', NC_DOUBLE, [d%cohort], 'per-cohort internal leaf water [kg/plant]')
      call def_var(d, 'wood_water_mass', NC_DOUBLE, [d%cohort], 'per-cohort internal wood water [kg/plant]')
      call def_var(d, 'leaf_surf_water', NC_DOUBLE, [d%cohort], 'leaf interception film [kg/plant]', units=FILM_UNITS)
      call def_var(d, 'wood_surf_water', NC_DOUBLE, [d%cohort], 'wood interception film [kg/plant]', units=FILM_UNITS)
      call def_var(d, 'dmax_psi_leaf',   NC_DOUBLE, [d%cohort],                                                 &
                   'yesterday daily-max leaf water potential [MPa] (drives beta_stomata)')
      !----- 5. Energy: the tissue energies, with their temperatures as twins. ----------------!
      call def_var(d, 'leaf_energy', NC_DOUBLE, [d%cohort], 'per-cohort leaf internal energy [J/plant] (prognostic)')
      call def_var(d, 'wood_energy', NC_DOUBLE, [d%cohort], 'per-cohort wood internal energy [J/plant] (prognostic)')
      call def_var(d, 'leaf_temp',   NC_DOUBLE, [d%cohort], 'per-cohort leaf temperature [K] (diagnosed twin)')
      call def_var(d, 'wood_temp',   NC_DOUBLE, [d%cohort], 'per-cohort wood temperature [K] (diagnosed twin)')
      !----- 6. Plastic leaf traits. ---------------------------------------------------------!
      call def_var(d, 'sla',                NC_DOUBLE, [d%cohort], 'specific leaf area [m2/kgC] (plastic)')
      call def_var(d, 'vcmax25',            NC_DOUBLE, [d%cohort], 'max carboxylation @25C [umol/m2/s] (plastic)')
      call def_var(d, 'rd25',               NC_DOUBLE, [d%cohort], 'leaf dark respiration @25C [umol/m2/s] (plastic)')
      call def_var(d, 'leaf_turnover_rate', NC_DOUBLE, [d%cohort], 'baseline leaf turnover rate [1/yr] (plastic)')
      !----- 7. Leaf phenology: the season's cue sums and running light mean. -----------------!
      call def_var(d, 'leaf_flush_tendency', NC_DOUBLE, [d%cohort], 'smoothed flush signal [0,1]')
      call def_var(d, 'leaf_shed_tendency',  NC_DOUBLE, [d%cohort], 'smoothed senescence signal [0,1]')
      call def_var(d, 'growing_degree_days', NC_DOUBLE, [d%cohort], 'warmth sum since midwinter [K day]')
      call def_var(d, 'cold_degree_days',    NC_DOUBLE, [d%cohort], 'cold sum since midsummer [K day]')
      call def_var(d, 'dry_psi_sum',         NC_DOUBLE, [d%cohort], 'predawn leaf psi below the TLP, summed [MPa day]')
      call def_var(d, 'wet_psi_sum',         NC_DOUBLE, [d%cohort], 'predawn leaf psi above the TLP, summed [MPa day]')
      call def_var(d, 'light_hours_mean',    NC_DOUBLE, [d%cohort],                                             &
                   'running-mean hours of light at the cohort top [h/day]')
      !----- 8. Growth history: growth_avg (the mortality predictor) and the ring buffer behind it. !
      call def_var(d, 'growth_avg',   NC_DOUBLE, [d%cohort], 'simple moving-average growth [cm/yr] (mortality predictor)')
      call def_var(d, 'growth_accum', NC_DOUBLE, [d%cohort], 'running sum of the growth window [cm/yr]')
      call def_var(d, 'growth_count', NC_INT,    [d%cohort], 'samples in the growth window')
      call def_var(d, 'growth_hist',  NC_DOUBLE, [d%cohort, d%growth_window], 'growth-rate ring buffer [cm/yr]')
      !----- 9. Fast-to-slow accumulators: the day's running maximum leaf potential. -----------!
      call def_var(d, 'dmax_psi_leaf_accum', NC_DOUBLE, [d%cohort],                                             &
                   'running daily-max leaf water potential accumulator [MPa]')
   end subroutine define_cohort_vars

   subroutine put_cohort_vars(r, poly)
      type(file_row_t), intent(in) :: r
      type(polygon_t),  intent(in) :: poly
      integer(ik) :: ncoh, nwin, ic
      ncoh = poly%cohort%n
      nwin = int(size(poly%cohort%growth_hist, 1), ik)
      associate (c => poly%cohort)
         !----- 1. Identity and density. -----------------------------------------------------!
         call put_int(r, 'pft',              c%pft(1:ncoh))
         call put_int(r, 'global_cohort_id', c%global_id(1:ncoh))
         call put_int(r, 'owner_patch',      c%owner_patch(1:ncoh))
         call put_dbl(r, 'nplant',           c%nplant(1:ncoh))
         !----- 2. Carbon pools. --------------------------------------------------------------!
         call put_dbl(r, 'leaf_carbon',          c%leaf_carbon(1:ncoh))
         call put_dbl(r, 'fineroot_carbon',      c%fineroot_carbon(1:ncoh))
         call put_dbl(r, 'wood_carbon',          c%wood_carbon(1:ncoh))
         call put_dbl(r, 'nonstructural_carbon', c%nonstructural_carbon(1:ncoh))
         !----- 3. Size. ----------------------------------------------------------------------!
         call put_dbl(r, 'dbh',             c%dbh(1:ncoh))
         call put_dbl(r, 'height',          c%height(1:ncoh))
         call put_dbl(r, 'basal_area',      c%basal_area(1:ncoh))
         call put_dbl(r, 'agb',             c%agb(1:ncoh))
         call put_dbl(r, 'leaf_area',       c%leaf_area(1:ncoh))
         call put_dbl(r, 'overtopping_lai', c%overtopping_lai(1:ncoh))
         !----- 4. Water. ---------------------------------------------------------------------!
         call put_dbl(r, 'leaf_water_mass', c%leaf_water_mass(1:ncoh))
         call put_dbl(r, 'wood_water_mass', c%wood_water_mass(1:ncoh))
         call put_dbl(r, 'leaf_surf_water', c%leaf_surf_water(1:ncoh))
         call put_dbl(r, 'wood_surf_water', c%wood_surf_water(1:ncoh))
         call put_dbl(r, 'dmax_psi_leaf',   c%dmax_psi_leaf(1:ncoh))
         !----- 5. Energy. --------------------------------------------------------------------!
         call put_dbl(r, 'leaf_energy', c%leaf_energy(1:ncoh))
         call put_dbl(r, 'wood_energy', c%wood_energy(1:ncoh))
         call put_dbl(r, 'leaf_temp',   c%leaf_temp(1:ncoh))
         call put_dbl(r, 'wood_temp',   c%wood_temp(1:ncoh))
         !----- 6. Plastic leaf traits. -------------------------------------------------------!
         call put_dbl(r, 'sla',                c%sla(1:ncoh))
         call put_dbl(r, 'vcmax25',            c%vcmax25(1:ncoh))
         call put_dbl(r, 'rd25',               c%rd25(1:ncoh))
         call put_dbl(r, 'leaf_turnover_rate', c%leaf_turnover_rate(1:ncoh))
         !----- 7. Leaf phenology. ------------------------------------------------------------!
         call put_dbl(r, 'leaf_flush_tendency', c%leaf_flush_tendency(1:ncoh))
         call put_dbl(r, 'leaf_shed_tendency',  c%leaf_shed_tendency(1:ncoh))
         call put_dbl(r, 'growing_degree_days', c%growing_degree_days(1:ncoh))
         call put_dbl(r, 'cold_degree_days',    c%cold_degree_days(1:ncoh))
         call put_dbl(r, 'dry_psi_sum',         c%dry_psi_sum(1:ncoh))
         call put_dbl(r, 'wet_psi_sum',         c%wet_psi_sum(1:ncoh))
         call put_dbl(r, 'light_hours_mean',    c%light_hours_mean(1:ncoh))
         !----- 8. Growth history. ------------------------------------------------------------!
         call put_dbl(r, 'growth_avg',   c%growth_avg(1:ncoh))
         call put_dbl(r, 'growth_accum', c%growth_accum(1:ncoh))
         call put_int(r, 'growth_count', c%growth_count(1:ncoh))
         if (nwin > 0_ik) then
            do ic = 1_ik, ncoh
               call put_dbl(r, 'growth_hist', c%growth_hist(1:nwin, ic), item=ic)
            end do
         end if
         !----- 9. Fast-to-slow accumulators. -------------------------------------------------!
         call put_dbl(r, 'dmax_psi_leaf_accum', c%dmax_psi_leaf_accum(1:ncoh))
      end associate
   end subroutine put_cohort_vars

   !----- The cohorts as the file holds them; what an older file lacks keeps alloc_polygon's    !
   !      value, or is set in section 4 as `held` records. --------------------------------------!
   subroutine get_cohort_vars(r, poly, cfg, held)
      type(file_row_t),    intent(in)    :: r
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(file_held_t),   intent(inout) :: held
      character(len=16) :: units
      integer(c_int)    :: vid
      integer(c_size_t) :: nwin_file
      integer(ik)       :: ncoh, nwin, i, gpos(1)
      ncoh = poly%cohort%n
      associate (c => poly%cohort)
         !----- 1. Identity and density. -----------------------------------------------------!
         call get_int(r, 'pft',              c%pft(1:ncoh))
         call get_int(r, 'global_cohort_id', c%global_id(1:ncoh))
         call get_int(r, 'owner_patch',      c%owner_patch(1:ncoh))
         call get_dbl(r, 'nplant',           c%nplant(1:ncoh))
         !----- 2. Carbon pools: a file without them places each cohort from its dbh. ----------!
         held%pools = file_has(r, 'leaf_carbon')
         if (held%pools) then
            call get_dbl(r, 'leaf_carbon',          c%leaf_carbon(1:ncoh))
            call get_dbl(r, 'fineroot_carbon',      c%fineroot_carbon(1:ncoh))
            call get_dbl(r, 'wood_carbon',          c%wood_carbon(1:ncoh))
            call get_dbl(r, 'nonstructural_carbon', c%nonstructural_carbon(1:ncoh))
         end if
         !----- 3. Size: derived from the pools in section 4, so only dbh (a file without pools  !
         !      places the cohort from it) and the overtopping LAI are read. ---------------------!
         call get_dbl(r, 'dbh', c%dbh(1:ncoh))
         held%overtopping_lai = file_has(r, 'overtopping_lai')
         if (held%overtopping_lai) call get_dbl(r, 'overtopping_lai', c%overtopping_lai(1:ncoh))
         !----- 4. Water. A file whose films carry no units holds them per m2 of ground; a file  !
         !      without the predawn potential keeps alloc_polygon's unset value (seeded from the !
         !      soil). ---------------------------------------------------------------------------!
         held%tissue_water = file_has(r, 'leaf_water_mass')
         call get_dbl_if_held(r, 'leaf_water_mass', c%leaf_water_mass(1:ncoh))
         call get_dbl_if_held(r, 'wood_water_mass', c%wood_water_mass(1:ncoh))
         call get_dbl_if_held(r, 'leaf_surf_water', c%leaf_surf_water(1:ncoh))
         call get_dbl_if_held(r, 'wood_surf_water', c%wood_surf_water(1:ncoh))
         if (nc_inq_varid_f(r%ncid, 'leaf_surf_water', vid) == NC_NOERR) then
            if (nc_get_att_text_f(r%ncid, vid, 'units', units) /= NC_NOERR) units = ''
            held%films_per_plant = trim(units) == FILM_UNITS
         end if
         call get_dbl_if_held(r, 'dmax_psi_leaf', c%dmax_psi_leaf(1:ncoh))
         !----- 5. Energy. --------------------------------------------------------------------!
         held%tissue_temps  = file_has(r, 'leaf_temp')
         held%tissue_energy = file_has(r, 'leaf_energy')
         call get_dbl_if_held(r, 'leaf_energy', c%leaf_energy(1:ncoh))
         call get_dbl_if_held(r, 'wood_energy', c%wood_energy(1:ncoh))
         call get_dbl_if_held(r, 'leaf_temp',   c%leaf_temp(1:ncoh))
         call get_dbl_if_held(r, 'wood_temp',   c%wood_temp(1:ncoh))
         !----- 6. Plastic leaf traits: the PFT's top-of-canopy values, then the file's where it  !
         !      has them; a file that holds a leaf lifespan [yr], `llspan`, gives the rate as its  !
         !      reciprocal. ----------------------------------------------------------------------!
         do i = 1_ik, ncoh
            c%sla(i)     = cfg%pft%sla(c%pft(i))
            c%vcmax25(i) = cfg%pft%vcmax25(c%pft(i))
            c%rd25(i)    = cfg%pft%rd25(c%pft(i))
            c%leaf_turnover_rate(i) = cfg%pft%leaf_turnover_rate_toc(c%pft(i))
         end do
         call get_dbl_if_held(r, 'sla',                c%sla(1:ncoh))
         call get_dbl_if_held(r, 'vcmax25',            c%vcmax25(1:ncoh))
         call get_dbl_if_held(r, 'rd25',               c%rd25(1:ncoh))
         call get_dbl_if_held(r, 'leaf_turnover_rate', c%leaf_turnover_rate(1:ncoh))
         if (file_has(r, 'llspan')) then
            call get_dbl(r, 'llspan', c%leaf_turnover_rate(1:ncoh))
            c%leaf_turnover_rate(1:ncoh) = 1.0_wp / c%leaf_turnover_rate(1:ncoh)
         end if
         !----- 7. Leaf phenology: a file without it keeps the birth values alloc_polygon set. ---!
         call get_dbl_if_held(r, 'leaf_flush_tendency', c%leaf_flush_tendency(1:ncoh))
         call get_dbl_if_held(r, 'leaf_shed_tendency',  c%leaf_shed_tendency(1:ncoh))
         call get_dbl_if_held(r, 'growing_degree_days', c%growing_degree_days(1:ncoh))
         call get_dbl_if_held(r, 'cold_degree_days',    c%cold_degree_days(1:ncoh))
         call get_dbl_if_held(r, 'dry_psi_sum',         c%dry_psi_sum(1:ncoh))
         call get_dbl_if_held(r, 'wet_psi_sum',         c%wet_psi_sum(1:ncoh))
         call get_dbl_if_held(r, 'light_hours_mean',    c%light_hours_mean(1:ncoh))
         !----- 8. Growth history: the ring buffer, with the polygon's write slot, when the file   !
         !      has one of this run's window length; section 4 otherwise seeds it full of the     !
         !      saved growth_avg. ----------------------------------------------------------------!
         call get_dbl(r, 'growth_avg', c%growth_avg(1:ncoh))
         nwin = growth_window_steps(cfg)
         held%growth_rings = file_has(r, 'growth_hist')
         if (held%growth_rings) then
            call nc_check(nc_inq_dimlen_f(r%ncid, 'growth_window', nwin_file), 'inq growth_window')
            held%growth_rings = int(nwin_file, ik) == nwin .and. int(size(c%growth_hist, 1), ik) == nwin
         end if
         if (held%growth_rings) then
            call get_dbl(r, 'growth_accum', c%growth_accum(1:ncoh))
            call get_int(r, 'growth_count', c%growth_count(1:ncoh))
            do i = 1_ik, ncoh
               call get_dbl(r, 'growth_hist', c%growth_hist(1:nwin, i), item=i)
            end do
            call get_int(r, 'growth_hist_pos', gpos)
            poly%growth_hist_pos = gpos(1)
         end if
         !----- 9. Fast-to-slow accumulators. -------------------------------------------------!
         call get_dbl_if_held(r, 'dmax_psi_leaf_accum', c%dmax_psi_leaf_accum(1:ncoh))
      end associate
   end subroutine get_cohort_vars

   !=======================================================================================!
   !  3. THE PATCHES, in meds_patch_state_vars' section order                              !
   !=======================================================================================!

   subroutine define_patch_vars(d)
      type(file_dims_t), intent(in) :: d
      !----- 1. Identity and area. -----------------------------------------------------------!
      call def_var(d, 'global_patch_id', NC_INT,    [d%patch], 'persistent patch id')
      call def_var(d, 'dist_type',       NC_INT,    [d%patch], 'disturbance type (1=primary,2=treefall)')
      call def_var(d, 'owner_site',      NC_INT,    [d%patch], 'owning site index (1-based)')
      call def_var(d, 'patch_area',      NC_DOUBLE, [d%patch], 'patch area fraction of its site')
      call def_var(d, 'patch_age',       NC_DOUBLE, [d%patch], 'time since last disturbance [yr]')
      !----- 3. Recruit bank. ----------------------------------------------------------------!
      call def_var(d, 'recruit_pool', NC_DOUBLE, [d%patch, d%pft], 'carry-forward recruit pool [plant/m2]')
      !----- 4-6. The fast reservoirs (canopy air, snow, soil, pond), with their diagnosed twins. !
      call def_var(d, 'cas_can_enthalpy', NC_DOUBLE, [d%patch], 'CAS specific enthalpy [J/kg] (prognostic)')
      call def_var(d, 'cas_can_shv',      NC_DOUBLE, [d%patch], 'CAS specific humidity [kg/kg] (prognostic)')
      call def_var(d, 'cas_can_co2',      NC_DOUBLE, [d%patch], 'CAS CO2 mixing ratio [umol/mol] (prognostic)')
      call def_var(d, 'cas_can_temp',     NC_DOUBLE, [d%patch], 'CAS temperature [K] (diagnosed twin)')
      call def_var(d, 'snow_nlayer', NC_INT,    [d%patch], 'active snow layer count (0 = no snow)')
      call def_var(d, 'snow_swe',    NC_DOUBLE, [d%patch, d%snow_layer], 'snow water-equivalent mass [kg/m2] (prognostic)')
      call def_var(d, 'snow_energy', NC_DOUBLE, [d%patch, d%snow_layer],                                   &
                   'snow extensive internal energy [J/m2] (prognostic)')
      call def_var(d, 'snow_depth',  NC_DOUBLE, [d%patch, d%snow_layer], 'snow geometric depth [m] (prognostic)')
      call def_var(d, 'snow_temp',   NC_DOUBLE, [d%patch, d%snow_layer], 'snow temperature [K] (diagnosed twin)')
      call def_var(d, 'snow_fliq',   NC_DOUBLE, [d%patch, d%snow_layer], 'snow liquid fraction [-] (diagnosed twin)')
      call def_var(d, 'soil_theta',  NC_DOUBLE, [d%patch, d%soil_layer], 'soil volumetric moisture [m3/m3] (prognostic)')
      call def_var(d, 'soil_energy', NC_DOUBLE, [d%patch, d%soil_layer],                                   &
                   'soil volumetric internal energy [J/m3] (prognostic)')
      call def_var(d, 'soil_temp',   NC_DOUBLE, [d%patch, d%soil_layer], 'soil temperature [K] (diagnosed twin)')
      call def_var(d, 'soil_fliq',   NC_DOUBLE, [d%patch, d%soil_layer], 'soil liquid-water fraction [-] (diagnosed twin)')
      call def_var(d, 'soil_w_surface',      NC_DOUBLE, [d%patch], 'ponded surface water [kg/m2]')
      call def_var(d, 'soil_w_surface_enth', NC_DOUBLE, [d%patch], 'ponded-water internal energy [J/m2]')
      !----- 7. Soil carbon (docs/science/soil_carbon.md §2). --------------------------------!
      call def_var(d, 'soilc_fast_grnd',   NC_DOUBLE, [d%patch], 'fast/metabolic litter carbon, above-ground [kgC/m2]')
      call def_var(d, 'soilc_fast_soil',   NC_DOUBLE, [d%patch], 'fast/metabolic litter carbon, below-ground [kgC/m2]')
      call def_var(d, 'soilc_struct_grnd', NC_DOUBLE, [d%patch], 'structural litter + CWD carbon, above-ground [kgC/m2]')
      call def_var(d, 'soilc_struct_soil', NC_DOUBLE, [d%patch], 'structural litter + CWD carbon, below-ground [kgC/m2]')
      call def_var(d, 'soilc_microbial',   NC_DOUBLE, [d%patch], 'microbial SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_slow',        NC_DOUBLE, [d%patch], 'slow/humified SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_passive',     NC_DOUBLE, [d%patch], 'passive SOM carbon [kgC/m2]')
      call def_var(d, 'soilc_lignin_grnd', NC_DOUBLE, [d%patch], 'structural litter lignin, above-ground [kgC/m2]')
      call def_var(d, 'soilc_lignin_soil', NC_DOUBLE, [d%patch], 'structural litter lignin, below-ground [kgC/m2]')
      !----- 8. The daily slow-to-fast rates: the turnover shed water and the slow loop's CO2. --!
      call def_var(d, 'shed_water_rate', NC_DOUBLE, [d%patch], 'turnover shed water to the ground [kg/m2/s]')
      call def_var(d, 'shed_water_enth', NC_DOUBLE, [d%patch], 'the enthalpy that water carries [W/m2]')
      call def_var(d, 'slow_co2_rate',   NC_DOUBLE, [d%patch], 'slow-loop CO2 flux the fast NEE adds [umol/m2/s]')
      !----- 10. Solver memory: the step the fast controller last accepted, its warm start. ----!
      call def_var(d, 'adapt_dt_last', NC_DOUBLE, [d%patch], 'fast-controller warm-start step [s] (0 = cold)')
   end subroutine define_patch_vars

   !----- The patches. Each per-patch number of a store goes through a contiguous buffer: a     !
   !      component section is never an actual argument (CLAUDE.md, compiler traps). ----------!
   subroutine put_patch_vars(r, poly)
      type(file_row_t), intent(in) :: r
      type(polygon_t),  intent(in) :: poly
      real(wp)    :: buf(poly%patch%n)
      integer(ik) :: ibuf(poly%patch%n), npat, npft, ip
      npat = poly%patch%n
      npft = poly%n_pft
      associate (p => poly%patch)
         !----- 1. Identity and area. ---------------------------------------------------------!
         call put_int(r, 'global_patch_id', p%global_id(1:npat))
         call put_int(r, 'dist_type',       p%dist_type(1:npat))
         call put_int(r, 'owner_site',      p%owner_site(1:npat))
         call put_dbl(r, 'patch_area',      p%area(1:npat))
         call put_dbl(r, 'patch_age',       p%age(1:npat))
         !----- 3. Recruit bank. --------------------------------------------------------------!
         do ip = 1_ik, npat
            call put_dbl(r, 'recruit_pool', p%recruit_pool(1:npft, ip), item=ip)
         end do
         !----- 4. Canopy air. ----------------------------------------------------------------!
         buf = p%cas(1:npat)%can_enthalpy
         call put_dbl(r, 'cas_can_enthalpy', buf)
         buf = p%cas(1:npat)%can_shv
         call put_dbl(r, 'cas_can_shv', buf)
         buf = p%cas(1:npat)%can_co2
         call put_dbl(r, 'cas_can_co2', buf)
         buf = p%cas(1:npat)%can_temp
         call put_dbl(r, 'cas_can_temp', buf)
         !----- 5. Snow. ----------------------------------------------------------------------!
         ibuf = p%snow(1:npat)%nlayer
         call put_int(r, 'snow_nlayer', ibuf)
         do ip = 1_ik, npat
            call put_dbl(r, 'snow_swe',    p%snow(ip)%swe(1:n_snow_layer_max),         item=ip)
            call put_dbl(r, 'snow_energy', p%snow(ip)%snow_energy(1:n_snow_layer_max), item=ip)
            call put_dbl(r, 'snow_depth',  p%snow(ip)%snow_depth(1:n_snow_layer_max),  item=ip)
            call put_dbl(r, 'snow_temp',   p%snow(ip)%snow_temp(1:n_snow_layer_max),   item=ip)
            call put_dbl(r, 'snow_fliq',   p%snow(ip)%snow_fliq(1:n_snow_layer_max),   item=ip)
         end do
         !----- 6. Soil water and heat, and the pond. -----------------------------------------!
         do ip = 1_ik, npat
            call put_dbl(r, 'soil_theta',  p%soil_water(ip)%theta(1:n_soil_layer_max),        item=ip)
            call put_dbl(r, 'soil_energy', p%soil_energy(ip)%soil_energy(1:n_soil_layer_max), item=ip)
            call put_dbl(r, 'soil_temp',   p%soil_energy(ip)%soil_temp(1:n_soil_layer_max),   item=ip)
            call put_dbl(r, 'soil_fliq',   p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max),   item=ip)
         end do
         buf = p%soil_water(1:npat)%w_surface
         call put_dbl(r, 'soil_w_surface', buf)
         buf = p%soil_water(1:npat)%w_surface_energy
         call put_dbl(r, 'soil_w_surface_enth', buf)
         !----- 7. Soil carbon. ---------------------------------------------------------------!
         buf = p%soil_carbon(1:npat)%fast_grnd_carbon
         call put_dbl(r, 'soilc_fast_grnd', buf)
         buf = p%soil_carbon(1:npat)%fast_soil_carbon
         call put_dbl(r, 'soilc_fast_soil', buf)
         buf = p%soil_carbon(1:npat)%struct_grnd_carbon
         call put_dbl(r, 'soilc_struct_grnd', buf)
         buf = p%soil_carbon(1:npat)%struct_soil_carbon
         call put_dbl(r, 'soilc_struct_soil', buf)
         buf = p%soil_carbon(1:npat)%microbial_carbon
         call put_dbl(r, 'soilc_microbial', buf)
         buf = p%soil_carbon(1:npat)%slow_carbon
         call put_dbl(r, 'soilc_slow', buf)
         buf = p%soil_carbon(1:npat)%passive_carbon
         call put_dbl(r, 'soilc_passive', buf)
         buf = p%soil_carbon(1:npat)%struct_grnd_lignin
         call put_dbl(r, 'soilc_lignin_grnd', buf)
         buf = p%soil_carbon(1:npat)%struct_soil_lignin
         call put_dbl(r, 'soilc_lignin_soil', buf)
         !----- 8. Daily slow-to-fast rates. --------------------------------------------------!
         call put_dbl(r, 'shed_water_rate', p%shed_water_rate(1:npat))
         call put_dbl(r, 'shed_water_enth', p%shed_water_enth(1:npat))
         call put_dbl(r, 'slow_co2_rate',   p%slow_co2_rate(1:npat))
         !----- 10. Solver memory. ------------------------------------------------------------!
         call put_dbl(r, 'adapt_dt_last', p%adapt_dt_last(1:npat))
      end associate
   end subroutine put_patch_vars

   !----- The patches as the file holds them. The fast reservoirs were written together, so the  !
   !      first stands for all, and the shed water and the controller's warm start come with     !
   !      them: a file without them leaves alloc_polygon's values for the caller to re-seed. ----!
   subroutine get_patch_vars(r, poly, held)
      type(file_row_t),  intent(in)    :: r
      type(polygon_t),   intent(inout) :: poly
      type(file_held_t), intent(inout) :: held
      real(wp)    :: buf(poly%patch%n), sc(poly%patch%n, 9)
      integer(ik) :: ibuf(poly%patch%n), npat, npft, ip
      npat = poly%patch%n
      npft = poly%n_pft
      held%fast = file_has(r, 'cas_can_enthalpy')
      associate (p => poly%patch)
         !----- 1. Identity and area; a file without owner_site has one site, alloc_polygon's   !
         !      default. -------------------------------------------------------------------------!
         call get_int(r, 'global_patch_id', p%global_id(1:npat))
         call get_int(r, 'dist_type',       p%dist_type(1:npat))
         if (file_has(r, 'owner_site')) call get_int(r, 'owner_site', p%owner_site(1:npat))
         call get_dbl(r, 'patch_area',      p%area(1:npat))
         call get_dbl(r, 'patch_age',       p%age(1:npat))
         !----- 3. Recruit bank. --------------------------------------------------------------!
         do ip = 1_ik, npat
            call get_dbl(r, 'recruit_pool', p%recruit_pool(1:npft, ip), item=ip)
         end do
         if (held%fast) then
            !----- 4. Canopy air. -------------------------------------------------------------!
            call get_dbl(r, 'cas_can_enthalpy', buf)
            p%cas(1:npat)%can_enthalpy = buf
            call get_dbl(r, 'cas_can_shv', buf)
            p%cas(1:npat)%can_shv = buf
            call get_dbl(r, 'cas_can_co2', buf)
            p%cas(1:npat)%can_co2 = buf
            call get_dbl(r, 'cas_can_temp', buf)
            p%cas(1:npat)%can_temp = buf
            !----- 5. Snow. -------------------------------------------------------------------!
            call get_int(r, 'snow_nlayer', ibuf)
            p%snow(1:npat)%nlayer = ibuf
            do ip = 1_ik, npat
               call get_dbl(r, 'snow_swe',    p%snow(ip)%swe(1:n_snow_layer_max),         item=ip)
               call get_dbl(r, 'snow_energy', p%snow(ip)%snow_energy(1:n_snow_layer_max), item=ip)
               call get_dbl(r, 'snow_depth',  p%snow(ip)%snow_depth(1:n_snow_layer_max),  item=ip)
               call get_dbl(r, 'snow_temp',   p%snow(ip)%snow_temp(1:n_snow_layer_max),   item=ip)
               call get_dbl(r, 'snow_fliq',   p%snow(ip)%snow_fliq(1:n_snow_layer_max),   item=ip)
            end do
            !----- 6. Soil water and heat, and the pond; a file without the pond's energy has it  !
            !      derived in section 4. ----------------------------------------------------------!
            do ip = 1_ik, npat
               call get_dbl(r, 'soil_theta',  p%soil_water(ip)%theta(1:n_soil_layer_max),        item=ip)
               call get_dbl(r, 'soil_energy', p%soil_energy(ip)%soil_energy(1:n_soil_layer_max), item=ip)
               call get_dbl(r, 'soil_temp',   p%soil_energy(ip)%soil_temp(1:n_soil_layer_max),   item=ip)
               call get_dbl(r, 'soil_fliq',   p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max),   item=ip)
            end do
            call get_dbl(r, 'soil_w_surface', buf)
            p%soil_water(1:npat)%w_surface = buf
            held%pond_energy = file_has(r, 'soil_w_surface_enth')
            if (held%pond_energy) then
               call get_dbl(r, 'soil_w_surface_enth', buf)
               p%soil_water(1:npat)%w_surface_energy = buf
            end if
         end if
         !----- 7. Soil carbon: a file without it restarts the pools at zero. -----------------!
         sc = 0.0_wp
         call get_dbl_if_held(r, 'soilc_fast_grnd',   sc(:,1))
         call get_dbl_if_held(r, 'soilc_fast_soil',   sc(:,2))
         call get_dbl_if_held(r, 'soilc_struct_grnd', sc(:,3))
         call get_dbl_if_held(r, 'soilc_struct_soil', sc(:,4))
         call get_dbl_if_held(r, 'soilc_microbial',   sc(:,5))
         call get_dbl_if_held(r, 'soilc_slow',        sc(:,6))
         call get_dbl_if_held(r, 'soilc_passive',     sc(:,7))
         call get_dbl_if_held(r, 'soilc_lignin_grnd', sc(:,8))
         call get_dbl_if_held(r, 'soilc_lignin_soil', sc(:,9))
         p%soil_carbon(1:npat)%fast_grnd_carbon   = sc(:,1)
         p%soil_carbon(1:npat)%fast_soil_carbon   = sc(:,2)
         p%soil_carbon(1:npat)%struct_grnd_carbon = sc(:,3)
         p%soil_carbon(1:npat)%struct_soil_carbon = sc(:,4)
         p%soil_carbon(1:npat)%microbial_carbon   = sc(:,5)
         p%soil_carbon(1:npat)%slow_carbon        = sc(:,6)
         p%soil_carbon(1:npat)%passive_carbon     = sc(:,7)
         p%soil_carbon(1:npat)%struct_grnd_lignin = sc(:,8)
         p%soil_carbon(1:npat)%struct_soil_lignin = sc(:,9)
         !----- 8. Daily slow-to-fast rates: the shed water, with the fast reservoirs (a file     !
         !      without its enthalpy has it derived in section 4); the slow loop's CO2. ----------!
         if (held%fast) then
            call get_dbl_if_held(r, 'shed_water_rate', p%shed_water_rate(1:npat))
            held%shed_energy = file_has(r, 'shed_water_enth')
            if (held%shed_energy) call get_dbl(r, 'shed_water_enth', p%shed_water_enth(1:npat))
         end if
         call get_dbl_if_held(r, 'slow_co2_rate', p%slow_co2_rate(1:npat))
         !----- 10. Solver memory: a file without the warm start cold-starts the controller. --!
         if (held%fast) call get_dbl_if_held(r, 'adapt_dt_last', p%adapt_dt_last(1:npat))
      end associate
   end subroutine get_patch_vars

   !=======================================================================================!
   !  4. AFTER THE READ: the geometry from the pools, what an older file lacked from what  !
   !     it held, and the maps, sort order and overtopping LAI                             !
   !=======================================================================================!

   subroutine derive_after_read(poly, cfg, held)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      type(file_held_t),   intent(in)    :: held
      real(wp)    :: t_birth, psi_birth
      integer(ik) :: ncoh, npat, nwin, i, ip
      ncoh = poly%cohort%n
      npat = poly%patch%n

      !----- Canopy films held per m2 of ground, to per plant. -------------------------------!
      if (.not. held%films_per_plant) then
         poly%cohort%leaf_surf_water(1:ncoh) = poly%cohort%leaf_surf_water(1:ncoh) / poly%cohort%nplant(1:ncoh)
         poly%cohort%wood_surf_water(1:ncoh) = poly%cohort%wood_surf_water(1:ncoh) / poly%cohort%nplant(1:ncoh)
      end if
      !----- The geometry from the pools, as the run derives it; a file without pools places    !
      !      each cohort on the allometry from its dbh. ----------------------------------------!
      do i = 1_ik, ncoh
         if (held%pools) then
            call set_cohort_size_from_carbon(poly%cohort, i)
         else
            call set_cohort_size(poly%cohort, i)
         end if
      end do
      !----- A growth ring buffer the file lacked, seeded full of the saved growth_avg. -------!
      if (.not. held%growth_rings) then
         nwin = growth_window_steps(cfg)
         do i = 1_ik, ncoh
            if (.not. is_unset_positive(poly%cohort%growth_avg(i))) then
               poly%cohort%growth_count(i)   = nwin
               poly%cohort%growth_accum(i)   = poly%cohort%growth_avg(i) * real(nwin, wp)
               poly%cohort%growth_hist(:, i) = poly%cohort%growth_avg(i)
            end if
         end do
      end if
      !----- The shed water's enthalpy, valued at the canopy air's temperature, and the pond's,  !
      !      its mass at the top soil layer's temperature (docs/science/order_of_processes.md §6, !
      !      The restart file). -------------------------------------------------------------!
      if (held%fast .and. .not. held%shed_energy) then
         do ip = 1_ik, npat
            poly%patch%shed_water_enth(ip) = poly%patch%shed_water_rate(ip)                                   &
                                           * internal_energy_liquid(poly%patch%cas(ip)%can_temp)
         end do
      end if
      if (held%fast .and. .not. held%pond_energy) then
         do ip = 1_ik, npat
            poly%patch%soil_water(ip)%w_surface_energy = poly%patch%soil_water(ip)%w_surface                 &
                  * internal_energy_liquid(poly%patch%soil_energy(ip)%soil_temp(1))
         end do
      end if
      !----- A file without the tissue water or temperatures starts each cohort in its patch's state,  !
      !      as a recruit is born: the canopy air and top soil layer just read, or, without those, the !
      !      run's starting state, from which init_fast_reservoirs re-seeds the patch and cohorts. ----!
      if (ncoh > 0_ik .and. npat > 0_ik .and. .not. (held%tissue_water .and. held%tissue_temps)) then
         do i = 1_ik, ncoh
            ip = poly%cohort%owner_patch(i)
            call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip), cfg,                &
                                      poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
            if (.not. held%tissue_temps) then
               poly%cohort%leaf_temp(i) = t_birth
               poly%cohort%wood_temp(i) = t_birth
            end if
            if (.not. held%tissue_water) then
               associate (hydro => cfg%pft%hydro(poly%cohort%pft(i)))
                  call set_cohort_tissue_water(poly%cohort, i, psi_birth, hydro%leaf_curve, hydro%wood_curve)
               end associate
            end if
         end do
      end if
      !----- The tissue energies are the state and the temperatures their twins. A file without the  !
      !      energies, or whose tissue water was set above, derives them from the temperatures at the  !
      !      masses; one without only the temperatures derives those from the energies. --------------!
      if (ncoh > 0_ik .and. .not. (held%tissue_energy .and. held%tissue_water .and. held%tissue_temps)) then
         do i = 1_ik, ncoh
            if (held%tissue_energy .and. held%tissue_water) then
               call set_cohort_tissue_temps(poly%cohort, i)
            else
               call hold_cohort_tissue_temps(poly%cohort, i)
            end if
         end do
      end if

      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
      !----- A file without overtopping_lai: compute it from the stand as read, so the traits'     !
      !      light environment is the stand's from the first step rather than an open sky. --------!
      if (ncoh > 0_ik .and. .not. held%overtopping_lai) call update_overtopping_lai(poly)
   end subroutine derive_after_read

end module meds_init_restart
