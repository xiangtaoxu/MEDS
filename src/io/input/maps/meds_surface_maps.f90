! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_surface_maps -- every static per-cell map a run reads, in one place: no other module opens  !
! a map (MEDS_POLYGON_RUNTIME_PLAN.md §14, RP2).                                                    !
!                                                                                          !
! A map is a variable in a netCDF file on the forcing grid: the file's lat and lon axes are the     !
! ED_ERA5land static file's, so a cell's value sits at the cell's own (row, col) and MEDS never     !
! regrids. A 2-D map is (lat, lon); a layered one is (layer, lat, lon). A value the map does not    !
! hold -- the variable's _FillValue or missing_value, or a NaN -- reads as NaN: filling gaps is the  !
! product's job, not MEDS's (§14.6, OR17), and the caller decides whether a missing value is an     !
! error. Packed values (scale_factor, add_offset) are unpacked.                                     !
!                                                                                          !
! A read takes the band of rows the cells span, every column, as one hyperslab: one read per       !
! variable for a whole region, which also serves a box that crosses 180 degrees.                   !
!==========================================================================================!
module meds_surface_maps
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use meds_kinds,            only : wp, ik
   use meds_netcdf_c
   use meds_forcing_types,    only : met_cells_t
   use meds_era5land_reader,  only : era5land_same_grid, era5land_keep_cells, ERA_OK, ERA_ERR_OPEN
   implicit none
   private

   public :: map_values_at_cells, map_layers_at_cells, apply_cell_mask, map_status_message
   public :: MAP_OK, MAP_ERR_OPEN, MAP_ERR_GRID, MAP_ERR_VAR

   integer(ik), parameter :: MAP_OK       = 0_ik
   integer(ik), parameter :: MAP_ERR_OPEN = 1_ik   !< the file (or the static file) cannot be opened
   integer(ik), parameter :: MAP_ERR_GRID = 2_ik   !< the file is not on the forcing grid
   integer(ik), parameter :: MAP_ERR_VAR  = 3_ik   !< the variable (or its layer dimension) is not there

contains

   !----- A 2-D (lat, lon) map variable at every cell: values(cell), NaN where the map is empty. -----!
   subroutine map_values_at_cells(path, var, static_path, cells, values, stat)
      character(len=*),      intent(in)  :: path, var, static_path
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:)
      integer(ik),           intent(out) :: stat
      real(wp), allocatable :: v(:,:)
      call read_cells(path, var, '', static_path, cells, v, stat)
      if (stat == MAP_OK) values = v(1, :)
   end subroutine map_values_at_cells

   !----- A layered (layer, lat, lon) map variable at every cell: values(layer, cell), NaN where the  !
   !      map is empty. `layer_dim` names the variable's first dimension. ------------------------------!
   subroutine map_layers_at_cells(path, var, layer_dim, static_path, cells, values, stat)
      character(len=*),      intent(in)  :: path, var, layer_dim, static_path
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:,:)
      integer(ik),           intent(out) :: stat
      call read_cells(path, var, layer_dim, static_path, cells, values, stat)
   end subroutine map_layers_at_cells

   !----- The cell mask of a region (§14 RP2, OR11): keep the cells whose mask value is at least      !
   !      mask_min. A cell the mask leaves empty is not kept. n_dropped counts the cells it removes. ----!
   subroutine apply_cell_mask(cells, mask_path, mask_var, mask_min, static_path, n_dropped, stat)
      type(met_cells_t), intent(inout) :: cells
      character(len=*),  intent(in)    :: mask_path, mask_var, static_path
      real(wp),          intent(in)    :: mask_min
      integer(ik),       intent(out)   :: n_dropped, stat
      real(wp), allocatable :: v(:)
      logical,  allocatable :: keep(:)
      n_dropped = 0_ik
      call map_values_at_cells(mask_path, mask_var, static_path, cells, v, stat)
      if (stat /= MAP_OK) return
      keep = .not. ieee_is_nan(v)
      where (keep) keep = v >= mask_min
      n_dropped = cells%ncell - int(count(keep), ik)
      call era5land_keep_cells(cells, keep)
   end subroutine apply_cell_mask

   !----- What a status code means, for the caller's message. ------------------------------------!
   function map_status_message(stat) result(msg)
      integer(ik), intent(in)       :: stat
      character(len=:), allocatable :: msg
      select case (stat)
      case (MAP_OK)       ; msg = 'read'
      case (MAP_ERR_OPEN) ; msg = 'cannot be opened'
      case (MAP_ERR_GRID) ; msg = 'is not on the forcing grid (its lat and lon must be the static file''s)'
      case (MAP_ERR_VAR)  ; msg = 'does not hold the variable (or its layer dimension)'
      case default        ; msg = 'cannot be read'
      end select
   end function map_status_message

   !=======================================================================================!
   !  The read: check the grid, read the band of rows the cells span (every column, every     !
   !  layer), then pick each cell's values. layer_dim = '' reads a 2-D variable as one layer.   !
   !=======================================================================================!
   subroutine read_cells(path, var, layer_dim, static_path, cells, values, stat)
      character(len=*),      intent(in)  :: path, var, layer_dim, static_path
      type(met_cells_t),     intent(in)  :: cells
      real(wp), allocatable, intent(out) :: values(:,:)
      integer(ik),           intent(out) :: stat
      integer(c_int)    :: ncid, vid, st
      integer(c_size_t) :: nlayer_c, start3(3), count3(3), start2(2), count2(2)
      integer(ik)       :: nlayer, r0, nr, k, l
      real(c_double), allocatable :: band(:,:,:)
      real(wp) :: fill, missing, scale, offset

      allocate(values(1, cells%ncell)) ; values = ieee_value(1.0_wp, ieee_quiet_nan)
      call era5land_same_grid(path, static_path, stat)
      if (stat /= ERA_OK) then
         stat = merge(MAP_ERR_OPEN, MAP_ERR_GRID, stat == ERA_ERR_OPEN) ; return
      end if
      stat = MAP_ERR_OPEN
      st = nc_open_f(trim(path), NC_NOWRITE, ncid) ; if (st /= NC_NOERR) return
      stat = MAP_ERR_VAR
      st = nc_inq_varid_f(ncid, trim(var), vid)
      if (st /= NC_NOERR) then ; st = nc_close(ncid) ; return ; end if
      nlayer = 1_ik
      if (len_trim(layer_dim) > 0) then
         st = nc_inq_dimlen_f(ncid, trim(layer_dim), nlayer_c)
         if (st /= NC_NOERR) then ; st = nc_close(ncid) ; return ; end if
         nlayer = int(nlayer_c, ik)
      end if
      if (cells%ncell == 0_ik) then
         st = nc_close(ncid) ; deallocate(values) ; allocate(values(nlayer, 0)) ; stat = MAP_OK ; return
      end if

      !----- The band of rows the cells span, every column: one hyperslab, x(lon, lat, layer). ------!
      r0 = minval(cells%row) ; nr = maxval(cells%row) - r0 + 1_ik
      allocate(band(cells%nlon, nr, nlayer))
      if (len_trim(layer_dim) > 0) then
         start3 = [0_c_size_t, int(r0, c_size_t), 0_c_size_t]
         count3 = [int(nlayer, c_size_t), int(nr, c_size_t), int(cells%nlon, c_size_t)]
         st = nc_get_vara_double(ncid, vid, start3, count3, band)
      else
         start2 = [int(r0, c_size_t), 0_c_size_t]
         count2 = [int(nr, c_size_t), int(cells%nlon, c_size_t)]
         st = nc_get_vara_double(ncid, vid, start2, count2, band)
      end if
      if (st /= NC_NOERR) then ; st = nc_close(ncid) ; return ; end if
      fill    = attribute(ncid, vid, '_FillValue')
      missing = attribute(ncid, vid, 'missing_value')
      scale   = attribute(ncid, vid, 'scale_factor')
      offset  = attribute(ncid, vid, 'add_offset')
      st = nc_close(ncid)

      deallocate(values) ; allocate(values(nlayer, cells%ncell))
      do k = 1_ik, cells%ncell
         do l = 1_ik, nlayer
            values(l, k) = real(band(cells%col(k) + 1_ik, cells%row(k) - r0 + 1_ik, l), wp)
         end do
      end do
      !----- Empty where the map says so; then unpack what is there. -------------------------------!
      if (.not. ieee_is_nan(fill))    where (values == fill)    values = ieee_value(1.0_wp, ieee_quiet_nan)
      if (.not. ieee_is_nan(missing)) where (values == missing) values = ieee_value(1.0_wp, ieee_quiet_nan)
      if (.not. ieee_is_nan(scale))  values = values * scale
      if (.not. ieee_is_nan(offset)) values = values + offset
      stat = MAP_OK

   contains

      !----- A numeric attribute of the variable, NaN when it has none. ---------------------------!
      real(wp) function attribute(ncid, vid, name) result(x)
         integer(c_int),   intent(in) :: ncid, vid
         character(len=*), intent(in) :: name
         real(c_double) :: a
         x = ieee_value(1.0_wp, ieee_quiet_nan)
         if (nc_get_att_double_f(ncid, vid, name, a) == NC_NOERR) x = real(a, wp)
      end function attribute

   end subroutine read_cells

end module meds_surface_maps
