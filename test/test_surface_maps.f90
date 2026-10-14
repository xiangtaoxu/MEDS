! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_surface_maps -- the one module that reads static per-cell maps (MEDS_POLYGON_RUNTIME_PLAN.md  !
! §14, RP2), on the synthetic archive's 4 x 18 grid.                                                  !
!                                                                                          !
! A 2-D map is read at a box's cells, each at its own (row, col); a cell the map leaves empty reads  !
! as NaN. A layered map gives every layer at every cell. A map off the forcing grid, a missing       !
! variable and a missing file are refused, each with its own status. The cell mask keeps the cells   !
! at or above its threshold, drops the empty one, and leaves the kept cells grouped by chunk the way !
! the forcing reader needs them. The work directory is the first argument.                          !
!==========================================================================================!
program test_surface_maps
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use meds_kinds,                 only : wp, ik
   use meds_test_assert,           only : check_true, test_report
   use meds_test_era5land_archive, only : NLAT, NLON, write_static, write_map
   use meds_forcing_types,         only : met_cells_t
   use meds_era5land_reader,       only : era5land_select_box, ERA_OK
   use meds_surface_maps,          only : map_values_at_cells, map_layers_at_cells, apply_cell_mask,  &
                                          MAP_OK, MAP_ERR_OPEN, MAP_ERR_GRID, MAP_ERR_VAR
   implicit none

   real(wp), parameter :: FILL = -9999.0_wp
   integer(ik), parameter :: NL = 3_ik
   character(len=512) :: work
   character(len=:), allocatable :: static
   type(met_cells_t)  :: cells, masked
   real(wp) :: grid2(NLON, NLAT, 1), grid3(NLON, NLAT, NL)
   real(wp), allocatable :: v(:), v3(:,:)
   integer(ik) :: st, r, c, l, k, n_dropped
   logical :: ok

   call get_command_argument(1, work)
   if (len_trim(work) == 0) error stop 'test_surface_maps: pass the work directory'
   call execute_command_line('rm -rf "'//trim(work)//'" && mkdir -p "'//trim(work)//'"')
   static = trim(work)//'/static.nc'
   call write_static(static, .true.)

   !----- The box: rows 0-3 (45 .. 15 N) and columns 0-3 (180 .. 120 W), every cell valid. ---------!
   call era5land_select_box(static, [45.0_wp, -180.0_wp, 15.0_wp, -120.0_wp], cells, st)
   call check_true('the box selects its 16 cells', st == ERA_OK .and. cells%ncell == 16_ik)

   !----- A 2-D map: 100 row + col, with (row 1, col 2) left empty. ----------------------------------!
   do r = 0_ik, NLAT - 1_ik
      do c = 0_ik, NLON - 1_ik
         grid2(c + 1_ik, r + 1_ik, 1) = 100.0_wp * real(r, wp) + real(c, wp)
      end do
   end do
   grid2(3, 2, 1) = FILL
   call write_map(trim(work)//'/map2.nc', 'value', grid2, FILL)
   call map_values_at_cells(trim(work)//'/map2.nc', 'value', static, cells, v, st)
   call check_true('a 2-D map is read', st == MAP_OK .and. size(v) == cells%ncell)
   ok = .true.
   do k = 1_ik, cells%ncell
      if (cells%row(k) == 1_ik .and. cells%col(k) == 2_ik) then
         ok = ok .and. ieee_is_nan(v(k))
      else
         ok = ok .and. v(k) == 100.0_wp * real(cells%row(k), wp) + real(cells%col(k), wp)
      end if
   end do
   call check_true('each cell reads its own (row, col); the empty cell reads NaN', ok)

   !----- A layered map: 1000 layer + 100 row + col. -----------------------------------------------!
   do l = 1_ik, NL
      grid3(:, :, l) = grid2(:, :, 1) + 1000.0_wp * real(l, wp)
   end do
   grid3(3, 2, :) = FILL
   call write_map(trim(work)//'/map3.nc', 'value', grid3, FILL, layer_dim='depth')
   call map_layers_at_cells(trim(work)//'/map3.nc', 'value', 'depth', static, cells, v3, st)
   call check_true('a layered map is read', st == MAP_OK .and. size(v3, 1) == NL .and. size(v3, 2) == cells%ncell)
   ok = .true.
   do k = 1_ik, cells%ncell
      do l = 1_ik, NL
         if (cells%row(k) == 1_ik .and. cells%col(k) == 2_ik) then
            ok = ok .and. ieee_is_nan(v3(l, k))
         else
            ok = ok .and. v3(l, k) == 1000.0_wp * real(l, wp) + 100.0_wp * real(cells%row(k), wp)       &
                                    + real(cells%col(k), wp)
         end if
      end do
   end do
   call check_true('every layer of every cell, NaN where empty', ok)

   !----- Refusals. ---------------------------------------------------------------------------------!
   call write_map(trim(work)//'/offgrid.nc', 'value', grid2, FILL, lat_shift=0.05_wp)
   call map_values_at_cells(trim(work)//'/offgrid.nc', 'value', static, cells, v, st)
   call check_true('a map off the forcing grid is refused', st == MAP_ERR_GRID)
   call map_values_at_cells(trim(work)//'/map2.nc', 'other', static, cells, v, st)
   call check_true('a missing variable is refused', st == MAP_ERR_VAR)
   call map_layers_at_cells(trim(work)//'/map3.nc', 'value', 'layer', static, cells, v3, st)
   call check_true('a missing layer dimension is refused', st == MAP_ERR_VAR)
   call map_values_at_cells(trim(work)//'/none.nc', 'value', static, cells, v, st)
   call check_true('a missing file is refused', st == MAP_ERR_OPEN)

   !----- The cell mask: keep the even columns at a 0.5 threshold; the empty cell is dropped. --------!
   do c = 0_ik, NLON - 1_ik
      grid2(c + 1_ik, :, 1) = merge(1.0_wp, 0.0_wp, modulo(c, 2_ik) == 0_ik)
   end do
   grid2(3, 2, 1) = FILL                                    ! (row 1, col 2): an even column, but empty
   call write_map(trim(work)//'/mask.nc', 'keep', grid2, FILL)
   masked = cells
   call apply_cell_mask(masked, trim(work)//'/mask.nc', 'keep', 0.5_wp, static, n_dropped, st)
   call check_true('the mask is read', st == MAP_OK)
   call check_true('it keeps the 7 even-column cells that are not empty', masked%ncell == 7_ik            &
                   .and. n_dropped == 9_ik)
   call check_true('every kept cell is in an even column, none is the empty one',                         &
                   all(modulo(masked%col, 2_ik) == 0_ik) .and.                                           &
                   .not. any(masked%row == 1_ik .and. masked%col == 2_ik))
   call check_true('the kept cells keep their order', all(masked%row(2:) >= masked%row(:masked%ncell - 1_ik)))
   call check_true('the kept cells are regrouped by chunk', masked%chunk_first(masked%nchunk + 1_ik) - 1_ik  &
                   == masked%ncell .and. size(masked%by_chunk) == masked%ncell)
   !----- A fraction and a threshold: 0.1 col, kept at >= 0.25, so columns 3 and up. -----------------!
   do c = 0_ik, NLON - 1_ik
      grid2(c + 1_ik, :, 1) = 0.1_wp * real(c, wp)
   end do
   call write_map(trim(work)//'/frac.nc', 'natural', grid2, FILL)
   masked = cells
   call apply_cell_mask(masked, trim(work)//'/frac.nc', 'natural', 0.25_wp, static, n_dropped, st)
   call check_true('a fraction map with a threshold keeps column 3 only', st == MAP_OK .and.               &
                   masked%ncell == 4_ik .and. all(masked%col == 3_ik))

   call test_report('test_surface_maps')
end program test_surface_maps
