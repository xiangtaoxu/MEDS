! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_surface_maps -- every static per-cell map a run reads, in one place: no other module opens  !
! a map (docs/configuration.md, Input maps).                                                        !
!                                                                                          !
! A map is a variable in a netCDF file on the forcing grid: the file's lat and lon axes are the     !
! ED_ERA5land static file's, so a cell's value sits at the cell's own (row, col) and MEDS never     !
! regrids. A 2-D map is (lat, lon); a layered one is (layer, lat, lon). A value the map does not    !
! hold -- the variable's _FillValue or missing_value, or a NaN -- reads as NaN: filling gaps is the  !
! product's job, not MEDS's, and the caller decides whether a missing value is an error. Packed     !
! values (scale_factor, add_offset) are unpacked.                                                   !
!                                                                                          !
! A read takes the band of rows the cells span, every column, as one hyperslab: one read per       !
! variable for a whole region, which also serves a box that crosses 180 degrees.                   !
!                                                                                          !
! The surface inputs ([maps]; docs/science/soil_biophysics.md §The column from the soil maps):  !
! read_map_inputs reads the soil maps and the Chave site index at every cell once, and          !
! apply_map_inputs turns one cell's values into that polygon's configuration -- its soil column,  !
! layer by layer, its bottom boundaries and its height allometry's site index. A run cell the maps !
! leave empty is an input error naming the cell.                                                    !
!==========================================================================================!
module meds_surface_maps
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use meds_kinds,            only : wp, ik
   use meds_netcdf_c
   use meds_forcing_types,    only : met_cells_t
   use meds_era5land_reader,  only : era5land_same_grid, era5land_keep_cells, era5land_select_site,   &
                                     era5land_static_file, ERA_OK, ERA_ERR_OPEN, ERA_ERR_NO_CELL
   use meds_config,           only : meds_config_t, soil_column_layers, soil_column_params
   use meds_config_check,     only : validate_config
   use meds_maps_opts,     only : maps_opts_t, soil_maps_on
   use meds_constants,        only : n_soil_layer_max
   use meds_soil_column_config, only : soil_params_t
   use meds_soil_hydr_lib,    only : soil_texture_t, soil_texture_to_hydraulics, SOIL_RETENTION_VG
   use meds_therm_lib,        only : soil_texture_to_thermal
   use meds_biophysics_opts,  only : SOIL_BC_DRAIN, SOIL_BC_AQUIFER, ANNUAL_DAMPING_DEPTH
   implicit none
   private

   public :: map_values_at_cells, map_layers_at_cells, apply_cell_mask, map_status_message
   public :: MAP_OK, MAP_ERR_OPEN, MAP_ERR_GRID, MAP_ERR_VAR
   public :: map_inputs_t, read_map_inputs, apply_map_inputs, site_map_inputs

   integer(ik), parameter :: MAP_OK       = 0_ik
   integer(ik), parameter :: MAP_ERR_OPEN = 1_ik   !< the file (or the static file) cannot be opened
   integer(ik), parameter :: MAP_ERR_GRID = 2_ik   !< the file is not on the forcing grid
   integer(ik), parameter :: MAP_ERR_VAR  = 3_ik   !< the variable (or its layer dimension) is not there

   !----- The soil maps' layers (ED-SoilMaps v2: SoilGrids 2.0's six depths, to 2 m). ------------!
   character(len=*), parameter :: SOIL_LAYER_DIM = 'depth', SOIL_LAYER_BOUNDS = 'depth_bnds'
   !----- The bottom water boundary of the soil depth map's bc_class. ---------------------------!
   integer(ik), parameter :: BC_FREE_DRAIN = 1_ik   !< q = K(theta)
   integer(ik), parameter :: BC_LEAKY      = 2_ik   !< q = f_leak K(theta)
   integer(ik), parameter :: BC_SEALED     = 3_ik   !< q = 0
   integer(ik), parameter :: BC_AQUIFER    = 4_ik   !< a water table at the column base

   !----- What the [maps] maps hold at each cell of a run, read once (read_map_inputs). --!
   type :: map_inputs_t
      logical :: soil  = .false.                          !< the soil maps were read
      logical :: chave = .false.                          !< the chave_e map was read
      real(wp), allocatable :: layer_top(:), layer_bottom(:)   !< [m] the soil maps' layers, positive down
      real(wp), allocatable :: soil_depth(:)              !< [m] (cell) the soil column's depth
      real(wp), allocatable :: bc_class(:)                !< [-] (cell) its bottom water boundary, BC_*
      real(wp), allocatable :: f_leak(:)                  !< [-] (cell) a leaky base's drainage factor
      real(wp), allocatable :: deep_temp(:)               !< [K] (cell) the deep soil temperature
      type(soil_texture_t), allocatable :: texture(:,:)   !< (map layer, cell) each layer's soil
      real(wp), allocatable :: chave_e(:)                 !< [-] (cell) Chave et al. (2014)'s E
   end type map_inputs_t

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

   !----- The cell mask of a region (docs/configuration.md, Regional runs): keep the cells whose    !
   !      mask value is at least mask_min. A cell the mask leaves empty is not kept. n_dropped counts !
   !      the cells it removes. -----------------------------------------------------------------------!
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
   !  read_map_inputs -- every [maps] map at every cell: the soil depth map's depth,      !
   !  bottom boundary, leaky drainage factor and deep temperature, the texture and chemistry      !
   !  maps' layers, and the Chave map's E. A map that cannot be read, or that holds no value at   !
   !  a cell, is an error (so is a leaky base whose f_leak is outside (0, 1]):                    !
   !  msg says which map, which variable and which cell, and is '' when every value is there.     !
   !=======================================================================================!
   subroutine read_map_inputs(opts, static_path, cells, inp, msg)
      type(maps_opts_t),          intent(in)  :: opts
      character(len=*),              intent(in)  :: static_path
      type(met_cells_t),             intent(in)  :: cells
      type(map_inputs_t),        intent(out) :: inp
      character(len=:), allocatable, intent(out) :: msg
      real(wp), allocatable :: v(:,:), top_chem(:), bottom_chem(:)
      integer(ik) :: k
      msg = ''
      if (soil_maps_on(opts)) then
         inp%soil = .true.
         call need(opts%soil_depth_map, 'soil_depth', '', v) ; if (len(msg) > 0) return
         inp%soil_depth = v(1, :)
         call need(opts%soil_depth_map, 'bc_class', '', v) ; if (len(msg) > 0) return
         inp%bc_class = v(1, :)
         do k = 1_ik, cells%ncell
            if (all(nint(inp%bc_class(k), ik) /= [BC_FREE_DRAIN, BC_LEAKY, BC_SEALED, BC_AQUIFER])) then
               msg = 'soil map '//trim(opts%soil_depth_map)//': bc_class is not 1-4 at '//cell_name(k)
               return
            end if
         end do
         call need(opts%soil_depth_map, 'f_leak', '', v) ; if (len(msg) > 0) return
         inp%f_leak = v(1, :)
         do k = 1_ik, cells%ncell
            if (nint(inp%bc_class(k), ik) == BC_LEAKY .and. .not. (inp%f_leak(k) > 0.0_wp .and.           &
                                                                  inp%f_leak(k) <= 1.0_wp)) then
               msg = 'soil map '//trim(opts%soil_depth_map)//': f_leak is not in (0, 1] at '//cell_name(k)// &
                     ', whose base is leaky'
               return
            end if
         end do
         call need(opts%soil_depth_map, 'deep_temp', '', v) ; if (len(msg) > 0) return
         inp%deep_temp = v(1, :)

         !----- The layers: the texture and the chemistry maps must share them. ------------------!
         call layer_bounds(opts%soil_texture_map, inp%layer_top, inp%layer_bottom)
         if (len(msg) > 0) return
         call layer_bounds(opts%soil_chemistry_map, top_chem, bottom_chem)
         if (len(msg) > 0) return
         if (size(top_chem) /= size(inp%layer_top)) then
            msg = 'soil maps '//trim(opts%soil_texture_map)//' and '//trim(opts%soil_chemistry_map)// &
                  ' do not have the same layers' ; return
         end if
         if (any(top_chem /= inp%layer_top) .or. any(bottom_chem /= inp%layer_bottom)) then
            msg = 'soil maps '//trim(opts%soil_texture_map)//' and '//trim(opts%soil_chemistry_map)// &
                  ' do not have the same layers' ; return
         end if
         allocate(inp%texture(size(inp%layer_top), cells%ncell))
         call need(opts%soil_texture_map, 'sand', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%sand = v
         call need(opts%soil_texture_map, 'silt', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%silt = v
         call need(opts%soil_texture_map, 'clay', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%clay = v
         call need(opts%soil_texture_map, 'bulk_density', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%bulk_density = v
         call need(opts%soil_texture_map, 'coarse_fragments', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%coarse_frag = v
         call need(opts%soil_chemistry_map, 'organic_carbon', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%org_carbon = v
         call need(opts%soil_chemistry_map, 'cec', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%cec = v
         call need(opts%soil_chemistry_map, 'ph_h2o', SOIL_LAYER_DIM, v) ; if (len(msg) > 0) return
         inp%texture%ph = v
      end if
      if (len_trim(opts%chave_e_map) > 0) then
         inp%chave = .true.
         call need(opts%chave_e_map, 'chave_e', '', v) ; if (len(msg) > 0) return
         inp%chave_e = v(1, :)
      end if

   contains

      !----- One variable at every cell, and every layer of it; msg on a failed read or a gap. ----!
      subroutine need(path, var, layer_dim, values)
         character(len=*),      intent(in)  :: path, var, layer_dim
         real(wp), allocatable, intent(out) :: values(:,:)
         integer(ik) :: st, c
         call read_cells(trim(path), var, layer_dim, static_path, cells, values, st)
         if (st /= MAP_OK) then
            msg = 'map '//trim(path)//' '//map_status_message(st)//': '//var ; return
         end if
         do c = 1_ik, cells%ncell
            if (any(ieee_is_nan(values(:, c)))) then
               msg = 'map '//trim(path)//': '//var//' has no value at '//cell_name(c) ; return
            end if
         end do
      end subroutine need

      !----- A soil map's layers: the tops and bottoms [m] of its depth_bnds(depth, nv). ----------!
      subroutine layer_bounds(path, top, bottom)
         character(len=*),      intent(in)  :: path
         real(wp), allocatable, intent(out) :: top(:), bottom(:)
         integer(c_int)    :: ncid, vid, st, st_close
         integer(c_size_t) :: n, start2(2), count2(2)
         real(c_double), allocatable :: b(:,:)
         msg = 'map '//trim(path)//' '//map_status_message(MAP_ERR_OPEN)
         if (nc_open_f(trim(path), NC_NOWRITE, ncid) /= NC_NOERR) return
         msg = 'map '//trim(path)//' '//map_status_message(MAP_ERR_VAR)//': '//SOIL_LAYER_BOUNDS
         st = nc_inq_dimlen_f(ncid, SOIL_LAYER_DIM, n)
         if (st == NC_NOERR) st = nc_inq_varid_f(ncid, SOIL_LAYER_BOUNDS, vid)
         if (st == NC_NOERR) then
            allocate(b(2, n))
            start2 = 0_c_size_t ; count2 = [n, 2_c_size_t]
            st = nc_get_vara_double(ncid, vid, start2, count2, b)
         end if
         st_close = nc_close(ncid)
         if (st /= NC_NOERR .or. st_close /= NC_NOERR) return
         top = real(b(1, :), wp) ; bottom = real(b(2, :), wp)
         msg = ''
      end subroutine layer_bounds

      !----- A cell, as a message names it. ----------------------------------------------------!
      function cell_name(c) result(name)
         integer(ik), intent(in) :: c
         character(len=:), allocatable :: name
         character(len=96) :: buf
         write(buf, '(a,f0.2,a,f0.2,a,i0,a,i0,a)') 'the cell (', cells%lat(c), ', ', cells%lon(c),      &
                                                  ') (row ', cells%row(c), ', col ', cells%col(c), ')'
         name = trim(buf)
      end function cell_name

   end subroutine read_map_inputs

   !=======================================================================================!
   !  apply_map_inputs -- cell k's inputs into the configuration its polygon runs.         !
   !                                                                                          !
   !  * The soil column ends at the map's soil depth, on the default layers (soil_column_layers); !
   !    each layer takes the map layers it overlaps, weighted by the overlap, and below the      !
   !    deepest map layer (2 m) the deepest one. Its curve and Ks come from the pedotransfer     !
   !    function, its thermal properties from the same texture (docs/science/soil_biophysics.md  !
   !    §The column from the soil maps).                                                         !
   !  * The bottom water boundary is the map's: drainage with f = 1 on a free-draining base, the  !
   !    map's f_leak on a leaky one and 0 on a sealed one, or the aquifer                          !
   !    (docs/science/soil_biophysics.md §Boundary conditions).                                   !
   !  * The Dirichlet thermal anchor sits at the map's deep temperature, ANNUAL_DAMPING_DEPTH /  !
   !    sqrt(2) below the bottom node, where it reflects the annual wave least                    !
   !    (docs/science/soil_biophysics.md eq. (13)).                                               !
   !  * The Chave height allometry takes the cell's E.                                          !
   !=======================================================================================!
   subroutine apply_map_inputs(inp, k, cfg)
      type(map_inputs_t), intent(in)    :: inp
      integer(ik),            intent(in)    :: k
      type(meds_config_t),    intent(inout) :: cfg
      type(soil_texture_t) :: tex
      type(soil_params_t)  :: soil
      real(wp)    :: bottoms(n_soil_layer_max), top
      integer(ik) :: n, l

      if (inp%soil) then
         associate (sc => cfg%soil_column)
            sc%depth = inp%soil_depth(k) ; sc%n_layer_bottoms = 0_ik
            sc%retention = SOIL_RETENTION_VG
            call soil_column_layers(sc, bottoms, n)
            top = 0.0_wp
            do l = 1_ik, n
               tex = layer_texture(inp, k, top, bottoms(l))
               call soil_texture_to_hydraulics(tex, sc%layer_theta_sat(l), sc%layer_theta_res(l),      &
                                                 sc%layer_curve_par_a(l), sc%layer_curve_par_n(l),       &
                                                 sc%layer_ksat(l))
               call soil_texture_to_thermal(tex%sand, tex%clay, sc%layer_theta_sat(l),                     &
                                              sc%layer_solid_conductivity(l), sc%layer_dry_conductivity(l),  &
                                              sc%layer_dry_heat_capacity(l))
               top = bottoms(l)
            end do
            sc%n_layer_values = n
         end associate
         select case (nint(inp%bc_class(k), ik))
         case (BC_AQUIFER) ; cfg%soil%water_bottom_bc = SOIL_BC_AQUIFER
         case (BC_LEAKY)   ; cfg%soil%water_bottom_bc = SOIL_BC_DRAIN ; cfg%soil%drain_factor = inp%f_leak(k)
         case (BC_SEALED)  ; cfg%soil%water_bottom_bc = SOIL_BC_DRAIN ; cfg%soil%drain_factor = 0.0_wp
         case default      ; cfg%soil%water_bottom_bc = SOIL_BC_DRAIN ; cfg%soil%drain_factor = 1.0_wp
         end select
         call soil_column_params(cfg%soil_column, soil)
         cfg%energy%deep_temp  = inp%deep_temp(k)
         cfg%energy%deep_depth = abs(soil%z_node(soil%n_active)) + ANNUAL_DAMPING_DEPTH / sqrt(2.0_wp)
      end if
      if (inp%chave) cfg%pft%allom(:)%chave_e = inp%chave_e(k)
   end subroutine apply_map_inputs

   !----- The soil of the column's layer [top, bottom]: the map layers it overlaps, each weighted  !
   !      by the overlap. The deepest map layer reaches down without end. ----------------------------!
   pure function layer_texture(inp, k, top, bottom) result(tex)
      type(map_inputs_t), intent(in) :: inp
      integer(ik),            intent(in) :: k
      real(wp),               intent(in) :: top, bottom
      type(soil_texture_t) :: tex
      real(wp)    :: w, wsum, hi
      integer(ik) :: m, nm
      nm = size(inp%layer_top, kind=ik)
      tex%sand = 0.0_wp ; tex%silt = 0.0_wp ; tex%clay = 0.0_wp ; tex%org_carbon = 0.0_wp
      tex%cec  = 0.0_wp ; tex%ph   = 0.0_wp ; tex%bulk_density = 0.0_wp ; tex%coarse_frag = 0.0_wp
      wsum = 0.0_wp
      do m = 1_ik, nm
         hi = bottom
         if (m < nm) hi = min(bottom, inp%layer_bottom(m))
         w = hi - max(top, inp%layer_top(m))
         if (w <= 0.0_wp) cycle
         wsum = wsum + w
         associate (t => inp%texture(m, k))
            tex%sand         = tex%sand         + w * t%sand
            tex%silt         = tex%silt         + w * t%silt
            tex%clay         = tex%clay         + w * t%clay
            tex%org_carbon   = tex%org_carbon   + w * t%org_carbon
            tex%cec          = tex%cec          + w * t%cec
            tex%ph           = tex%ph           + w * t%ph
            tex%bulk_density = tex%bulk_density + w * t%bulk_density
            tex%coarse_frag  = tex%coarse_frag  + w * t%coarse_frag
         end associate
      end do
      tex%sand = tex%sand / wsum ; tex%silt = tex%silt / wsum ; tex%clay = tex%clay / wsum
      tex%org_carbon = tex%org_carbon / wsum ; tex%cec = tex%cec / wsum ; tex%ph = tex%ph / wsum
      tex%bulk_density = tex%bulk_density / wsum ; tex%coarse_frag = tex%coarse_frag / wsum
   end function layer_texture

   !=======================================================================================!
   !  site_map_inputs -- a site run's [maps] inputs, at the forcing cell its [site]       !
   !  location reads (the archive's own choice, era5land_select_site), into cfg, which is then   !
   !  validated again. msg is '' on success, else why the inputs could not be had.               !
   !=======================================================================================!
   subroutine site_map_inputs(cfg, msg)
      type(meds_config_t),           intent(inout) :: cfg
      character(len=:), allocatable, intent(out)   :: msg
      type(met_cells_t)      :: cell
      type(map_inputs_t) :: inp
      character(len=:), allocatable :: static
      real(wp)    :: distance_km
      integer(ik) :: st
      static = trim(era5land_static_file(cfg%forcing))
      call era5land_select_site(static, cfg%forcing%latitude, cfg%forcing%longitude,              &
                                cfg%forcing%max_distance_km, cell, distance_km, st)
      if (st /= ERA_OK) then
         if (st == ERA_ERR_NO_CELL) then
            msg = 'no valid forcing cell near the site for the [maps] maps'
         else
            msg = 'cannot read the forcing''s static file '//static//' for the [maps] maps'
         end if
         return
      end if
      call read_map_inputs(cfg%maps, static, cell, inp, msg)
      if (len(msg) > 0) return
      call apply_map_inputs(inp, 1_ik, cfg)
      call validate_config(cfg)
   end subroutine site_map_inputs

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
