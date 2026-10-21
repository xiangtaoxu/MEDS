! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init_restart -- the state (restart) file: every polygon's prognostic state at an instant, in  !
! a netCDF-4 checkpoint a resumed run continues from exactly, through the netCDF C library !
! (docs/science/order_of_processes.md §6, The restart file). restart_file_open, _put and     !
! _close write a file a polygon's row at a time; read_restart_file reads a row, or the whole   !
! of a file without the polygon axis (one site's); restart_file_polygons and                 !
! restart_file_refusal check a file before a run starts from it. Restart only: diagnostic    !
! output is meds_output_*.                                                                 !
!==========================================================================================!
module meds_init_restart
   use iso_c_binding, only : c_int, c_size_t, c_double
   use meds_kinds,    only : wp, ik
   use meds_state_defaults, only : is_unset_positive
   use meds_config,   only : meds_config_t, growth_window_steps
   use meds_time,     only : meds_time_t, time_to_string, time_to_stamp, time_to_decimal_year
   use meds_netcdf_c, only : nc_check, nc_close, nc_create_f, nc_def_dim_f, nc_def_var_f, nc_def_var_chunking,   &
                             nc_def_var_deflate, nc_enddef, nc_get_att_text_f, nc_get_vara_double,               &
                             nc_get_vara_int, nc_inq_dimlen, nc_inq_dimlen_f, nc_inq_varid_f, nc_open_f,          &
                             nc_put_att_text_f, nc_put_vara_double, nc_put_vara_int, NC_CHUNKED, NC_CLOBBER,      &
                             NC_DOUBLE, NC_GLOBAL, NC_INT, NC_NETCDF4, NC_NOERR, NC_NOWRITE
   use meds_polygon_state_vars, only : polygon_t
   use meds_patch_state_vars, only : patch_physical_state
   use meds_cohort_state_vars,      only : set_cohort_size, set_cohort_size_from_carbon, set_cohort_tissue_water, &
                                           set_cohort_tissue_temps, hold_cohort_tissue_temps
   use meds_site_state_vars,        only : make_site_surface
   use meds_polygon_state_vars,     only : alloc_polygon, rebuild_patch_map, rebuild_cohort_map
   use meds_demography_update,      only : update_overtopping_lai
   use meds_constants,     only : n_soil_layer_max, n_snow_layer_max
   use meds_therm_lib,              only : internal_energy_liquid
   use meds_demography_cohort_fusefiss,   only : sort_cohorts
   implicit none
   private

   public :: read_restart_file
   public :: restart_file_t, restart_file_open, restart_file_write, restart_file_close, restart_file_slots

   character(len=*), parameter :: STATE_TITLE = 'MEDS instantaneous state (restart) file'
   integer,          parameter :: N_META_INT = 11, N_META_REAL = 2
   character(len=*), parameter :: FILM_UNITS  = 'kg/plant'   !< the canopy films' units attribute
   integer(c_int),   parameter :: DEFLATE_LEVEL = 1_c_int    !< the variables' zlib level

   !----- An open state file, its polygon axis ahead of every other. --------------------------------!
   type :: restart_file_t
      integer(c_int)     :: ncid   = -1_c_int
      character(len=512) :: path   = ''
   end type restart_file_t

contains

   !=======================================================================================!
   !  State (restart) files                                                                 !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Create <dir>/<prefix>-S-YYYYMMDDHHMMSS.nc and define every dimension and variable.    !
   ! The cohort and patch axes are n_cohort and n_patch long (at least 1), the largest of  !
   ! the polygons'. A `polygon` axis leads every variable, `polygon_id` names its rows,    !
   ! and the variables are compressed, since most polygons fill only part of the padded    !
   ! axes. Why each group below persists: docs/science/order_of_processes.md §6, The       !
   ! restart file.                                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine restart_file_open(sf, dir, prefix, now, n_pft, n_window, n_cohort, n_patch, n_site,     &
                              polygon_id, polygon_lat, polygon_lon, restructure_pending, restructure_new_year)
      type(restart_file_t),    intent(out) :: sf
      character(len=*),      intent(in)  :: dir, prefix
      type(meds_time_t),     intent(in)  :: now
      integer(ik),           intent(in)  :: n_pft, n_window, n_cohort, n_patch, n_site
      integer(ik),           intent(in)  :: polygon_id(:)                  !< the polygons' ids ...
      real(wp),              intent(in)  :: polygon_lat(:), polygon_lon(:)  !< ... and locations [deg]
      logical,     optional, intent(in)  :: restructure_pending, restructure_new_year
      character(len=5) :: pending
      integer(c_int) :: ncid, d_poly, d_cohort, d_patch, d_site, d_pft, d_mi, d_mr, d_soill, d_snowl, d_gwin
      integer(c_int) :: v_id, v_lat, v_lon
      integer(c_size_t) :: start1(1), count1(1)
      real(c_double), allocatable :: coord(:)

      !----- The stamp is the simulated date: <prefix>-S-YYYYMMDDHHMMSS.nc. -----------------!
      sf%path = trim(dir)//'/'//trim(prefix)//'-S-'//time_to_stamp(now)//'.nc'
      call nc_check(nc_create_f(trim(sf%path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'state nc_create')
      sf%ncid = ncid
      call nc_check(nc_def_dim_f(ncid, 'polygon', int(size(polygon_id), c_size_t), d_poly), 'state dim polygon')
      call nc_check(nc_def_dim_f(ncid, 'cohort', int(max(n_cohort,1_ik), c_size_t), d_cohort), 'state dim cohort')
      call nc_check(nc_def_dim_f(ncid, 'patch',  int(max(n_patch,1_ik), c_size_t), d_patch),   'state dim patch')
      call nc_check(nc_def_dim_f(ncid, 'site',   int(max(n_site,1_ik), c_size_t),  d_site),    'state dim site')
      call nc_check(nc_def_dim_f(ncid, 'pft',    int(n_pft, c_size_t),             d_pft),     'state dim pft')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_int',  int(N_META_INT, c_size_t), d_mi), 'state dim mi')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_real', int(N_META_REAL, c_size_t), d_mr), 'state dim mr')
      call nc_check(nc_def_dim_f(ncid, 'soil_layer', int(n_soil_layer_max, c_size_t), d_soill), 'state dim soill')
      call nc_check(nc_def_dim_f(ncid, 'snow_layer', int(n_snow_layer_max, c_size_t), d_snowl), 'state dim snowl')
      call nc_check(nc_def_dim_f(ncid, 'growth_window', int(max(n_window,1_ik), c_size_t), d_gwin), 'state dim gwin')
      call coord_var('polygon_id', NC_INT,    'forcing cell row-major index: row * n_lon + col', v_id)
      call coord_var('lat',        NC_DOUBLE, 'the polygon''s latitude [degrees_north]',          v_lat)
      call coord_var('lon',        NC_DOUBLE, 'the polygon''s longitude [degrees_east]',          v_lon)
      call coord_var('n_cohort',   NC_INT,    'the polygon''s live cohorts (the rest is fill)',    v_id)
      call coord_var('n_patch',    NC_INT,    'the polygon''s live patches (the rest is fill)',    v_id)
      call coord_var('n_site',     NC_INT,    'the polygon''s sites (the rest is fill)',           v_id)
      call coord_var('n_window',   NC_INT,    'the polygon''s growth-window length',              v_id)
      call nc_check(nc_inq_varid_f(ncid, 'polygon_id', v_id), 'state inq polygon_id')

      call dv('meta_int',         NC_INT,    [d_mi],                                   &
              '[n_cohort,n_patch,n_pft,next_cohort_id,next_patch_id,year,month,day,hour,minute,second]')
      call dv('meta_real',        NC_DOUBLE, [d_mr],     '[area of the sites (1), decimal_year]')
      !----- The sites: each one's area and the surface it was made with (its layer count is checked on read). !
      call dv('site_area',        NC_DOUBLE, [d_site],   'site area fraction of the polygon')
      call dv('site_n_soil',      NC_INT,    [d_site],   'the site''s active soil layers')
      call dv('site_soil_depth',  NC_DOUBLE, [d_site],   'the site''s soil column depth [m]')
      call dv('site_elevation',   NC_DOUBLE, [d_site],   'the site''s elevation [m]')
      call dv('site_slope',       NC_DOUBLE, [d_site],   'the site''s slope [degrees]')
      call dv('site_aspect',      NC_DOUBLE, [d_site],   'the site''s aspect [degrees from north]')
      call dv('pft',              NC_INT,    [d_cohort], 'plant functional type index')
      call dv('nplant',           NC_DOUBLE, [d_cohort], 'plant number density [plant/m2]')
      call dv('dbh',              NC_DOUBLE, [d_cohort], 'diameter at breast height [cm]')
      call dv('growth_avg',       NC_DOUBLE, [d_cohort], 'simple moving-average growth [cm/yr] (mortality predictor)')
      call dv('owner_patch',      NC_INT,    [d_cohort], 'owning patch index (1-based)')
      call dv('global_cohort_id', NC_INT,    [d_cohort], 'persistent cohort id')
      call dv('sla',              NC_DOUBLE, [d_cohort], 'specific leaf area [m2/kgC] (plastic)')
      call dv('vcmax25',          NC_DOUBLE, [d_cohort], 'max carboxylation @25C [umol/m2/s] (plastic)')
      call dv('rd25',             NC_DOUBLE, [d_cohort], 'leaf dark respiration @25C [umol/m2/s] (plastic)')
      call dv('leaf_turnover_rate', NC_DOUBLE, [d_cohort], 'baseline leaf turnover rate [1/yr] (plastic)')
      !----- Each cohort's tissue water, and its energies with their temperatures as twins. ---------!
      call dv('leaf_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal leaf water [kg/plant]')
      call dv('wood_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal wood water [kg/plant]')
      call dv('leaf_energy',      NC_DOUBLE, [d_cohort], 'per-cohort leaf internal energy [J/plant] (prognostic)')
      call dv('wood_energy',      NC_DOUBLE, [d_cohort], 'per-cohort wood internal energy [J/plant] (prognostic)')
      call dv('leaf_temp',        NC_DOUBLE, [d_cohort], 'per-cohort leaf temperature [K] (diagnosed twin)')
      call dv('wood_temp',        NC_DOUBLE, [d_cohort], 'per-cohort wood temperature [K] (diagnosed twin)')
      !----- The stomatal feedback's predawn potential and accumulator (docs/science/leaf_gas_exchange.md §4.2). !
      call dv('dmax_psi_leaf',   NC_DOUBLE, [d_cohort], &
            'yesterday daily-max leaf water potential [MPa] (drives beta_stomata)')
      call dv('dmax_psi_leaf_accum', NC_DOUBLE, [d_cohort], 'running daily-max leaf water potential accumulator [MPa]')
      !----- The phenology memory: the season's cue sums and running light mean. ---------------------!
      call dv('leaf_flush_tendency', NC_DOUBLE, [d_cohort], 'smoothed flush signal [0,1]')
      call dv('leaf_shed_tendency',  NC_DOUBLE, [d_cohort], 'smoothed senescence signal [0,1]')
      call dv('growing_degree_days', NC_DOUBLE, [d_cohort], 'warmth sum since midwinter [K day]')
      call dv('cold_degree_days',    NC_DOUBLE, [d_cohort], 'cold sum since midsummer [K day]')
      call dv('dry_psi_sum',         NC_DOUBLE, [d_cohort], 'predawn leaf psi below the TLP, summed [MPa day]')
      call dv('wet_psi_sum',         NC_DOUBLE, [d_cohort], 'predawn leaf psi above the TLP, summed [MPa day]')
      call dv('light_hours_mean', NC_DOUBLE, [d_cohort], 'running-mean hours of light at the cohort top [h/day]')
      !----- The carbon pools as the run holds them; height, basal_area, agb, leaf_area for readers. !
      call dv('height',               NC_DOUBLE, [d_cohort], 'height [m]')
      call dv('basal_area',           NC_DOUBLE, [d_cohort], 'basal area per plant [cm2]')
      call dv('agb',                  NC_DOUBLE, [d_cohort], 'aboveground biomass per plant [kgC]')
      call dv('leaf_area',            NC_DOUBLE, [d_cohort], 'leaf area per plant [m2]')
      call dv('leaf_carbon',          NC_DOUBLE, [d_cohort], 'leaf carbon per plant [kgC]')
      call dv('fineroot_carbon',      NC_DOUBLE, [d_cohort], 'fine-root carbon per plant [kgC]')
      call dv('wood_carbon',          NC_DOUBLE, [d_cohort], 'wood carbon per plant [kgC]')
      call dv('nonstructural_carbon', NC_DOUBLE, [d_cohort], 'nonstructural carbon per plant [kgC]')
      call dv('overtopping_lai',      NC_DOUBLE, [d_cohort], 'LAI of the taller cohorts in the patch [m2/m2]')
      !----- The canopy interception films ([fast].canopy_water_on), per plant and marked so. -------!
      call dv('leaf_surf_water',      NC_DOUBLE, [d_cohort], 'leaf interception film [kg/plant]', units=FILM_UNITS)
      call dv('wood_surf_water',      NC_DOUBLE, [d_cohort], 'wood interception film [kg/plant]', units=FILM_UNITS)
      !----- The growth-rate ring buffer behind growth_avg (the mortality predictor), its write slot. !
      call dv('growth_accum',         NC_DOUBLE, [d_cohort], 'running sum of the growth window [cm/yr]')
      call dv('growth_count',         NC_INT,    [d_cohort], 'samples in the growth window')
      call dv('growth_hist',          NC_DOUBLE, [d_cohort, d_gwin], 'growth-rate ring buffer [cm/yr]')
      call dv('growth_hist_pos',      NC_INT,    [d_mr], 'ring-buffer write slot (element 1)')
      !----- Thermal acclimation's running mean, one per site (docs/science/leaf_gas_exchange.md §2). !
      call dv('t_growth_avg', NC_DOUBLE, [d_site], 'the site''s running-mean growth temperature [K]')
      call dv('patch_area',       NC_DOUBLE, [d_patch],  'patch area fraction of its site')
      call dv('owner_site',       NC_INT,    [d_patch],  'owning site index (1-based)')
      call dv('patch_age',        NC_DOUBLE, [d_patch],  'time since last disturbance [yr]')
      call dv('dist_type',        NC_INT,    [d_patch],  'disturbance type (1=primary,2=treefall)')
      call dv('global_patch_id',  NC_INT,    [d_patch],  'persistent patch id')
      call dv('recruit_pool',     NC_DOUBLE, [d_patch, d_pft], 'carry-forward recruit pool [plant/m2]')
      call dv('slow_co2_rate',    NC_DOUBLE, [d_patch], 'slow-loop CO2 flux the fast NEE adds [umol/m2/s]')
      !----- The fast reservoirs (canopy air, soil, pond, snow), with their diagnosed twins. ----------!
      call dv('cas_can_enthalpy', NC_DOUBLE, [d_patch], 'CAS specific enthalpy [J/kg] (prognostic)')
      call dv('cas_can_shv',      NC_DOUBLE, [d_patch], 'CAS specific humidity [kg/kg] (prognostic)')
      call dv('cas_can_co2',      NC_DOUBLE, [d_patch], 'CAS CO2 mixing ratio [umol/mol] (prognostic)')
      call dv('cas_can_temp',     NC_DOUBLE, [d_patch], 'CAS temperature [K] (diagnosed twin)')
      call dv('soil_theta',       NC_DOUBLE, [d_patch, d_soill], 'soil volumetric moisture [m3/m3] (prognostic)')
      call dv('soil_w_surface',   NC_DOUBLE, [d_patch], 'ponded surface water [kg/m2]')
      call dv('soil_w_surface_enth', NC_DOUBLE, [d_patch], 'ponded-water internal energy [J/m2]')
      call dv('soil_energy',      NC_DOUBLE, [d_patch, d_soill], 'soil volumetric internal energy [J/m3] (prognostic)')
      call dv('soil_temp',        NC_DOUBLE, [d_patch, d_soill], 'soil temperature [K] (diagnosed twin)')
      call dv('soil_fliq',        NC_DOUBLE, [d_patch, d_soill], 'soil liquid-water fraction [-] (diagnosed twin)')
      call dv('snow_swe',         NC_DOUBLE, [d_patch, d_snowl], 'snow water-equivalent mass [kg/m2] (prognostic)')
      call dv('snow_energy',      NC_DOUBLE, [d_patch, d_snowl], 'snow extensive internal energy [J/m2] (prognostic)')
      call dv('snow_depth',       NC_DOUBLE, [d_patch, d_snowl], 'snow geometric depth [m] (prognostic)')
      call dv('snow_temp',        NC_DOUBLE, [d_patch, d_snowl], 'snow temperature [K] (diagnosed twin)')
      call dv('snow_fliq',        NC_DOUBLE, [d_patch, d_snowl], 'snow liquid fraction [-] (diagnosed twin)')
      call dv('snow_nlayer',      NC_INT,    [d_patch], 'active snow layer count (0 = no snow)')
      !----- The turnover shed water the slow step hands the fast loop for the day. ------------------!
      call dv('shed_water_rate',  NC_DOUBLE, [d_patch], 'turnover shed water to the ground [kg/m2/s]')
      call dv('shed_water_enth',  NC_DOUBLE, [d_patch], 'the enthalpy that water carries [W/m2]')
      !----- The fast controller's warm start: the step it last accepted in each patch. --------------!
      call dv('adapt_dt_last',    NC_DOUBLE, [d_patch], 'fast-controller warm-start step [s] (0 = cold)')
      !----- The soil-carbon pools (docs/science/soil_carbon.md §2). ---------------------------------!
      call dv('soilc_fast_grnd',   NC_DOUBLE, [d_patch], 'fast/metabolic litter carbon, above-ground [kgC/m2]')
      call dv('soilc_fast_soil',   NC_DOUBLE, [d_patch], 'fast/metabolic litter carbon, below-ground [kgC/m2]')
      call dv('soilc_struct_grnd', NC_DOUBLE, [d_patch], 'structural litter + CWD carbon, above-ground [kgC/m2]')
      call dv('soilc_struct_soil', NC_DOUBLE, [d_patch], 'structural litter + CWD carbon, below-ground [kgC/m2]')
      call dv('soilc_microbial',   NC_DOUBLE, [d_patch], 'microbial SOM carbon [kgC/m2]')
      call dv('soilc_slow',        NC_DOUBLE, [d_patch], 'slow/humified SOM carbon [kgC/m2]')
      call dv('soilc_passive',     NC_DOUBLE, [d_patch], 'passive SOM carbon [kgC/m2]')
      call dv('soilc_lignin_grnd', NC_DOUBLE, [d_patch], 'structural litter lignin, above-ground [kgC/m2]')
      call dv('soilc_lignin_soil', NC_DOUBLE, [d_patch], 'structural litter lignin, below-ground [kgC/m2]')
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, 'title',                            &
                    int(len_trim(STATE_TITLE), c_size_t), STATE_TITLE), 'state title')
      pending = 'none'
      if (present(restructure_pending)) then
         if (restructure_pending) pending = 'month'
         if (restructure_pending .and. present(restructure_new_year)) then
            if (restructure_new_year) pending = 'year'
         end if
      end if
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, 'restructure_pending',              &
                    int(len_trim(pending), c_size_t), trim(pending)), 'state restructure_pending')
      call nc_check(nc_enddef(ncid), 'state enddef')
      start1 = 0_c_size_t
      count1 = int(size(polygon_id), c_size_t)
      call nc_check(nc_put_vara_int(ncid, v_id, start1, count1, polygon_id), 'put polygon_id')
      coord = real(polygon_lat, c_double)
      call nc_check(nc_put_vara_double(ncid, v_lat, start1, count1, coord), 'put lat')
      coord = real(polygon_lon, c_double)
      call nc_check(nc_put_vara_double(ncid, v_lon, start1, count1, coord), 'put lon')

   contains
      !----- A per-polygon variable: one value per row, with its long_name. ------------------------!
      subroutine coord_var(name, xtype, lname, v)
         character(len=*), intent(in)  :: name, lname
         integer(c_int),   intent(in)  :: xtype
         integer(c_int),   intent(out) :: v
         integer(c_int) :: dims(1)
         dims = d_poly
         call nc_check(nc_def_var_f(ncid, name, xtype, 1_c_int, dims, v), 'state def '//name)
         call nc_check(nc_put_att_text_f(ncid, v, 'long_name', int(len_trim(lname), c_size_t), lname), &
                       'state long_name '//name)
      end subroutine coord_var
      !----- Define a state variable on its own dimensions, behind the polygon axis, with its      !
      !      long_name and, where the reader needs it, its units. ----------------------------------!
      subroutine dv(name, xtype, dimids, lname, units)
         character(len=*), intent(in) :: name, lname
         integer(c_int),   intent(in) :: xtype, dimids(:)
         character(len=*), intent(in), optional :: units
         integer(c_int)    :: v, dims(size(dimids) + 1)
         integer(c_size_t) :: chunk(size(dimids) + 1)
         integer           :: nd, k
         dims(1) = d_poly
         dims(2:) = dimids
         nd = size(dimids) + 1
         call nc_check(nc_def_var_f(ncid, name, xtype, int(nd, c_int), dims(1:nd), v), 'state def '//name)
         call nc_check(nc_put_att_text_f(ncid, v, 'long_name',                              &
                       int(len_trim(lname), c_size_t), lname), 'state long_name '//name)
         if (present(units)) call nc_check(nc_put_att_text_f(ncid, v, 'units',              &
                                           int(len_trim(units), c_size_t), units), 'state units '//name)
         !----- One polygon's row per chunk: each polygon writes its own (docs/science/order_of_processes.md §6). !
         if (nd > 1) then
            chunk(1) = 1_c_size_t
            do k = 2, nd
               call nc_check(nc_inq_dimlen(ncid, dims(k), chunk(k)), 'state dimlen '//name)
               chunk(k) = max(chunk(k), 1_c_size_t)
            end do
            call nc_check(nc_def_var_chunking(ncid, v, NC_CHUNKED, chunk(1:nd)), 'state chunk '//name)
         end if
         call nc_check(nc_def_var_deflate(ncid, v, 1_c_int, 1_c_int, DEFLATE_LEVEL), 'state deflate '//name)
      end subroutine dv
   end subroutine restart_file_open

   !---------------------------------------------------------------------------------------!
   ! Write one polygon's state into its row (1-based) of an open file.                         !
   !---------------------------------------------------------------------------------------!
   subroutine restart_file_write(sf, poly, now, polygon)
      type(restart_file_t),    intent(in) :: sf
      type(polygon_t),       intent(in) :: poly
      type(meds_time_t),     intent(in) :: now
      integer(ik),           intent(in) :: polygon
      integer(ik) :: ncoh, npat, npft, ip, ic, nwin, slot, ns, meta_i(N_META_INT)
      real(wp)    :: meta_r(N_META_REAL)

      ncoh = poly%cohort%n
      npat = poly%patch%n
      npft = poly%n_pft
      ns = poly%site%n
      nwin = int(size(poly%cohort%growth_hist, 1), ik)
      slot = polygon - 1_ik

      meta_i = [ncoh, npat, npft, poly%next_cohort_id, poly%next_patch_id,                     &
                now%year, now%month, now%day, now%hour, now%minute, now%second]
      call put_count('n_cohort', ncoh)
      call put_count('n_patch', npat)
      call put_count('n_window', nwin)
      call put_count('n_site', ns)
      meta_r = [sum(poly%site%area(1:ns)), time_to_decimal_year(now)]
      call put_int('meta_int', meta_i)
      call put_dbl('meta_real', meta_r)
      block                                 ! the sites, from contiguous buffers as below
         integer(ik) :: nsoil(ns), is
         real(wp)    :: depth(ns), elev(ns), slope(ns), aspect(ns)
         do is = 1_ik, ns
            nsoil(is)  = poly%site%surface(is)%soil%n_active
            depth(is)  = poly%site%surface(is)%soil_column%depth
            elev(is)   = poly%site%surface(is)%elevation_m
            slope(is)  = poly%site%surface(is)%slope_deg
            aspect(is) = poly%site%surface(is)%aspect_deg
         end do
         call put_dbl('site_area',       poly%site%area(1:ns))
         call put_int('site_n_soil',     nsoil)
         call put_dbl('site_soil_depth', depth)
         call put_dbl('site_elevation',  elev)
         call put_dbl('site_slope',      slope)
         call put_dbl('site_aspect',     aspect)
      end block
      call put_dbl('t_growth_avg', poly%site%t_growth_avg(1:ns))
      if (ncoh > 0_ik) then
         associate (c => poly%cohort)
            call put_int('pft',         c%pft(1:ncoh))
            call put_dbl('nplant',      c%nplant(1:ncoh))
            call put_dbl('dbh',         c%dbh(1:ncoh))
            call put_dbl('growth_avg',  c%growth_avg(1:ncoh))
            call put_int('owner_patch', c%owner_patch(1:ncoh))
            call put_int('global_cohort_id', c%global_id(1:ncoh))
            call put_dbl('sla',         c%sla(1:ncoh))
            call put_dbl('vcmax25',     c%vcmax25(1:ncoh))
            call put_dbl('rd25',        c%rd25(1:ncoh))
            call put_dbl('leaf_turnover_rate', c%leaf_turnover_rate(1:ncoh))
            call put_dbl('leaf_water_mass', c%leaf_water_mass(1:ncoh))
            call put_dbl('wood_water_mass', c%wood_water_mass(1:ncoh))
            call put_dbl('leaf_energy', c%leaf_energy(1:ncoh))
            call put_dbl('wood_energy', c%wood_energy(1:ncoh))
            call put_dbl('leaf_temp',   c%leaf_temp(1:ncoh))
            call put_dbl('wood_temp',   c%wood_temp(1:ncoh))
            call put_dbl('dmax_psi_leaf',       c%dmax_psi_leaf(1:ncoh))
            call put_dbl('dmax_psi_leaf_accum', c%dmax_psi_leaf_accum(1:ncoh))
            call put_dbl('leaf_flush_tendency', c%leaf_flush_tendency(1:ncoh))
            call put_dbl('leaf_shed_tendency',  c%leaf_shed_tendency(1:ncoh))
            call put_dbl('growing_degree_days', c%growing_degree_days(1:ncoh))
            call put_dbl('cold_degree_days',    c%cold_degree_days(1:ncoh))
            call put_dbl('dry_psi_sum',         c%dry_psi_sum(1:ncoh))
            call put_dbl('wet_psi_sum',         c%wet_psi_sum(1:ncoh))
            call put_dbl('light_hours_mean',    c%light_hours_mean(1:ncoh))
            call put_dbl('height',               c%height(1:ncoh))
            call put_dbl('basal_area',           c%basal_area(1:ncoh))
            call put_dbl('agb',                  c%agb(1:ncoh))
            call put_dbl('leaf_area',            c%leaf_area(1:ncoh))
            call put_dbl('leaf_carbon',          c%leaf_carbon(1:ncoh))
            call put_dbl('fineroot_carbon',      c%fineroot_carbon(1:ncoh))
            call put_dbl('wood_carbon',          c%wood_carbon(1:ncoh))
            call put_dbl('nonstructural_carbon', c%nonstructural_carbon(1:ncoh))
            call put_dbl('overtopping_lai',      c%overtopping_lai(1:ncoh))
            call put_dbl('leaf_surf_water',      c%leaf_surf_water(1:ncoh))
            call put_dbl('wood_surf_water',      c%wood_surf_water(1:ncoh))
            call put_dbl('growth_accum',         c%growth_accum(1:ncoh))
            call put_int('growth_count',         c%growth_count(1:ncoh))
            if (nwin > 0_ik) then
               do ic = 1_ik, ncoh
                  call put_dbl('growth_hist', c%growth_hist(1:nwin, ic), ic - 1_ik)
               end do
            end if
         end associate
      end if
      call put_int('growth_hist_pos', [poly%growth_hist_pos])
      if (npat > 0_ik) then
         associate (p => poly%patch)
            call put_dbl('patch_area',      p%area(1:npat))
            call put_int('owner_site',      p%owner_site(1:npat))
            call put_dbl('patch_age',       p%age(1:npat))
            call put_int('dist_type',       p%dist_type(1:npat))
            call put_int('global_patch_id', p%global_id(1:npat))
            call put_dbl('slow_co2_rate',   p%slow_co2_rate(1:npat))
            do ip = 1_ik, npat
               call put_dbl('recruit_pool', p%recruit_pool(1:npft, ip), ip - 1_ik)
            end do
            block                                 ! the fast reservoirs, through contiguous buffers:
               real(wp) :: fc(npat, 4_ik)          ! no component section into a C call (CLAUDE.md, compiler traps)
               real(wp) :: fw(npat, 5_ik)
               real(wp) :: fsoil(npat, n_soil_layer_max, 4_ik)
               real(wp) :: fsnow(npat, n_snow_layer_max, 5_ik)
               integer(ik) :: fnl(npat)
               do ip = 1_ik, npat
                  fc(ip,1) = p%cas(ip)%can_enthalpy
                  fc(ip,2) = p%cas(ip)%can_shv
                  fc(ip,3) = p%cas(ip)%can_co2
                  fc(ip,4) = p%cas(ip)%can_temp
                  fw(ip,1) = p%soil_water(ip)%w_surface
                  fw(ip,2) = p%soil_water(ip)%w_surface_energy
                  fw(ip,3) = p%shed_water_rate(ip)
                  fw(ip,4) = p%adapt_dt_last(ip)
                  fw(ip,5) = p%shed_water_enth(ip)
                  fsoil(ip,:,1) = p%soil_water(ip)%theta(1:n_soil_layer_max)
                  fsoil(ip,:,2) = p%soil_energy(ip)%soil_energy(1:n_soil_layer_max)
                  fsoil(ip,:,3) = p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max)
                  fsoil(ip,:,4) = p%soil_energy(ip)%soil_temp(1:n_soil_layer_max)
                  fsnow(ip,:,1) = p%snow(ip)%swe(1:n_snow_layer_max)
                  fsnow(ip,:,2) = p%snow(ip)%snow_energy(1:n_snow_layer_max)
                  fsnow(ip,:,3) = p%snow(ip)%snow_depth(1:n_snow_layer_max)
                  fsnow(ip,:,4) = p%snow(ip)%snow_temp(1:n_snow_layer_max)
                  fsnow(ip,:,5) = p%snow(ip)%snow_fliq(1:n_snow_layer_max)
                  fnl(ip) = p%snow(ip)%nlayer
               end do
               call put_dbl('cas_can_enthalpy', fc(:,1))
               call put_dbl('cas_can_shv',      fc(:,2))
               call put_dbl('cas_can_co2',      fc(:,3))
               call put_dbl('cas_can_temp',     fc(:,4))
               call put_dbl('soil_w_surface',   fw(:,1))
               call put_dbl('soil_w_surface_enth', fw(:,2))
               call put_dbl('shed_water_rate',  fw(:,3))
               call put_dbl('shed_water_enth',  fw(:,5))
               call put_dbl('adapt_dt_last',    fw(:,4))
               call put_int('snow_nlayer',      fnl)
               do ip = 1_ik, npat
                  call put_dbl('soil_theta',  fsoil(ip,:,1), ip - 1_ik)
                  call put_dbl('soil_energy', fsoil(ip,:,2), ip - 1_ik)
                  call put_dbl('soil_fliq',   fsoil(ip,:,3), ip - 1_ik)
                  call put_dbl('soil_temp',   fsoil(ip,:,4), ip - 1_ik)
                  call put_dbl('snow_swe',    fsnow(ip,:,1), ip - 1_ik)
                  call put_dbl('snow_energy', fsnow(ip,:,2), ip - 1_ik)
                  call put_dbl('snow_depth',  fsnow(ip,:,3), ip - 1_ik)
                  call put_dbl('snow_temp',   fsnow(ip,:,4), ip - 1_ik)
                  call put_dbl('snow_fliq',   fsnow(ip,:,5), ip - 1_ik)
               end do
            end block
            block                              ! the soil carbon, through a contiguous buffer as above
               real(wp) :: sc(npat, 9_ik)
               do ip = 1_ik, npat
                  sc(ip,1) = p%soil_carbon(ip)%fast_grnd_carbon
                  sc(ip,2) = p%soil_carbon(ip)%fast_soil_carbon
                  sc(ip,3) = p%soil_carbon(ip)%struct_grnd_carbon
                  sc(ip,4) = p%soil_carbon(ip)%struct_soil_carbon
                  sc(ip,5) = p%soil_carbon(ip)%microbial_carbon
                  sc(ip,6) = p%soil_carbon(ip)%slow_carbon
                  sc(ip,7) = p%soil_carbon(ip)%passive_carbon
                  sc(ip,8) = p%soil_carbon(ip)%struct_grnd_lignin
                  sc(ip,9) = p%soil_carbon(ip)%struct_soil_lignin
               end do
               call put_dbl('soilc_fast_grnd',   sc(:,1))
               call put_dbl('soilc_fast_soil',   sc(:,2))
               call put_dbl('soilc_struct_grnd', sc(:,3))
               call put_dbl('soilc_struct_soil', sc(:,4))
               call put_dbl('soilc_microbial',   sc(:,5))
               call put_dbl('soilc_slow',        sc(:,6))
               call put_dbl('soilc_passive',     sc(:,7))
               call put_dbl('soilc_lignin_grnd', sc(:,8))
               call put_dbl('soilc_lignin_soil', sc(:,9))
            end block
         end associate
      end if

   contains
      !----- A polygon's count, at its row. ----------------------------------------------------------!
      subroutine put_count(name, n)
         character(len=*), intent(in) :: name
         integer(ik),      intent(in) :: n
         integer(c_size_t) :: st(1), ct(1)
         integer(ik)       :: x(1)
         st = int(slot, c_size_t)
         ct = 1_c_size_t
         x = n
         call nc_check(nc_put_vara_int(sf%ncid, varid(name), st, ct, x), 'put '//name)
      end subroutine put_count
      !----- One double vector into the variable `name`: its first entries (a per-cohort, per-patch !
      !      or meta vector), or with `row`, row `row` (0-based) of a 2-D (cohort or patch, layer)     !
      !      variable; behind the polygon's slot. ---------------------------------------------------!
      subroutine put_dbl(name, x, row)
         character(len=*),      intent(in) :: name
         real(wp),              intent(in) :: x(:)
         integer(ik), optional, intent(in) :: row
         integer(c_size_t) :: st(3), ct(3)
         integer :: n
         call where_to(size(x), row, st, ct, n)
         call nc_check(nc_put_vara_double(sf%ncid, varid(name), st(1:n), ct(1:n), x), 'put '//name)
      end subroutine put_dbl
      subroutine put_int(name, x, row)
         character(len=*),      intent(in) :: name
         integer(ik),           intent(in) :: x(:)
         integer(ik), optional, intent(in) :: row
         integer(c_size_t) :: st(3), ct(3)
         integer :: n
         call where_to(size(x), row, st, ct, n)
         call nc_check(nc_put_vara_int(sf%ncid, varid(name), st(1:n), ct(1:n), x), 'put '//name)
      end subroutine put_int
      subroutine where_to(nx, row, st, ct, n)
         integer,               intent(in)  :: nx
         integer(ik), optional, intent(in)  :: row
         integer(c_size_t),     intent(out) :: st(3), ct(3)
         integer,               intent(out) :: n
         n = 1
         st(1) = int(slot, c_size_t)
         ct(1) = 1_c_size_t
         if (present(row)) then
            n = n + 1
            st(n) = int(row, c_size_t)
            ct(n) = 1_c_size_t
         end if
         n = n + 1
         st(n) = 0_c_size_t
         ct(n) = int(nx, c_size_t)
      end subroutine where_to
      integer(c_int) function varid(name)
         character(len=*), intent(in) :: name
         call nc_check(nc_inq_varid_f(sf%ncid, name, varid), 'state inq '//name)
      end function varid
   end subroutine restart_file_write

   subroutine restart_file_close(sf)
      type(restart_file_t), intent(inout) :: sf
      call nc_check(nc_close(sf%ncid), 'state nc_close')
      sf%ncid = -1_c_int
      write(*,'(2a)') ' state : ', trim(sf%path)
   end subroutine restart_file_close

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

   !---------------------------------------------------------------------------------------!
   ! Reconstruct a polygon from a state file written by restart_file_write: row `polygon`      !
   ! (1-based) of a file with the polygon axis, or the whole of a file without it. Each    !
   ! site's surface is made from `site_cfg` (else `cfg`) and the geometry from the carbon  !
   ! pools; the maps and sort order are rebuilt, and the restructuring owed is the file's  !
   ! restructure_pending (docs/science/order_of_processes.md §5). found = .false. if the   !
   ! file cannot be opened or has another layout (restart_file_refusal says why); a PFT,     !
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
      character(len=16) :: pending, units_attr
      character(len=160) :: msg
      integer(c_int) :: ncid, vid, st
      integer(ik)    :: ncoh, npat, npft, ip, i, nwin, slot, ns, is, meta_i(N_META_INT)
      real(wp)       :: meta_r(N_META_REAL)
      logical        :: fast_ok, pools_ok, olai_ok, ghist_ok, region, tissue_water_ok, tissue_temp_ok
      logical        :: tissue_energy_ok
      integer(c_size_t) :: nwin_file

      found = .false.
      restart_time = meds_time_t()
      tissue_water_ok = .true.
      tissue_temp_ok = .true.
      tissue_energy_ok = .true.
      if (present(fast_found)) fast_found = .false.
      region = present(polygon)
      slot = 0_ik
      if (region) slot = polygon - 1_ik
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) return
      !----- The polygon axis, where a file has it, leads every variable: reading the other kind   !
      !      of file would take the wrong elements. --------------------------------------------------!
      if (has_polygon_axis(ncid) .neqv. region) then
         st = nc_close(ncid)
         return
      end if
      pending = 'none'
      if (nc_get_att_text_f(ncid, NC_GLOBAL, 'restructure_pending', pending) /= NC_NOERR) pending = 'none'

      call gv_int(ncid, 'meta_int',  int(N_META_INT, ik),  meta_i)
      call gv_dbl(ncid, 'meta_real', int(N_META_REAL, ik), meta_r)
      ncoh = meta_i(1)
      npat = meta_i(2)
      npft = meta_i(3)
      if (npft /= cfg%pft%n) error stop 'read_restart_file: PFT count in state file /= config PFT count'

      !----- The sites: a file without the site axis holds one. ---------------------------------!
      ns = 1_ik
      if (nc_inq_varid_f(ncid, 'n_site', vid) == NC_NOERR) ns = polygon_count(vid)
      if (present(site_cfg)) then
         if (size(site_cfg, kind=ik) /= ns) then
            write(msg, '(a,i0,a,i0,a)') 'read_restart_file: the state file holds ', ns,                     &
                                        ' sites and the configuration ', size(site_cfg), ' for the polygon'
            error stop trim(msg)
         end if
      end if
      call alloc_polygon(poly, cfg%pft%n, max(ncoh, 1_ik), max(npat, 1_ik), growth_window_steps(cfg), &
                      cfg%pft, n_site=ns)
      do is = 1_ik, ns
         if (present(site_cfg)) then
            call make_site_surface(poly%site%surface(is), site_cfg(is))
         else
            call make_site_surface(poly%site%surface(is), cfg)
         end if
      end do
      poly%cohort%n = ncoh
      poly%patch%n = npat
      poly%restructure_pending  = trim(pending) == 'month' .or. trim(pending) == 'year'
      poly%restructure_new_year = trim(pending) == 'year'
      !----- Growth temperature: a file without it keeps alloc_polygon's unseeded sentinel. -----!
      call gv_dbl_opt(ncid, 't_growth_avg', ns, poly%site%t_growth_avg(1:ns))
      poly%next_cohort_id = meta_i(4)
      poly%next_patch_id = meta_i(5)
      if (nc_inq_varid_f(ncid, 'site_area', vid) == NC_NOERR) then
         call gv_dbl(ncid, 'site_area', ns, poly%site%area(1:ns))
         !----- Each site's soil state is per layer of the column it was written on. -------------!
         block
            integer(ik) :: nsoil(ns)
            call gv_int(ncid, 'site_n_soil', ns, nsoil)
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
         poly%site%area(1) = meta_r(1)              ! one site, every patch in it
      end if
      restart_time = meds_time_t(year=meta_i(6), month=meta_i(7), day=meta_i(8),               &
                                 hour=meta_i(9), minute=meta_i(10), second=meta_i(11))

      if (ncoh > 0_ik) then
         associate (c => poly%cohort)
            call gv_int (ncid, 'pft',              ncoh, c%pft(1:ncoh))
            call gv_dbl (ncid, 'nplant',           ncoh, c%nplant(1:ncoh))
            call gv_dbl (ncid, 'dbh',              ncoh, c%dbh(1:ncoh))
            call gv_dbl (ncid, 'growth_avg',       ncoh, c%growth_avg(1:ncoh))  ! ring buffer reseeded below
            call gv_int (ncid, 'owner_patch',      ncoh, c%owner_patch(1:ncoh))
            call gv_int (ncid, 'global_cohort_id', ncoh, c%global_id(1:ncoh))
            !----- Plastic leaf traits: the PFT's top-of-canopy values, then the file's where it has !
            !      them. -------------------------------------------------------------------------!
            do i = 1_ik, ncoh
               c%sla(i)     = cfg%pft%sla(c%pft(i))
               c%vcmax25(i) = cfg%pft%vcmax25(c%pft(i))
               c%rd25(i)    = cfg%pft%rd25(c%pft(i))
               c%leaf_turnover_rate(i) = cfg%pft%leaf_turnover_rate_toc(c%pft(i))
            end do
            call gv_dbl_opt(ncid, 'sla',     ncoh, c%sla(1:ncoh))
            call gv_dbl_opt(ncid, 'vcmax25', ncoh, c%vcmax25(1:ncoh))
            call gv_dbl_opt(ncid, 'rd25',    ncoh, c%rd25(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_turnover_rate', ncoh, c%leaf_turnover_rate(1:ncoh))
            !----- A file that holds a leaf lifespan [yr], `llspan`, gives the rate as its reciprocal. -!
            if (nc_inq_varid_f(ncid, 'llspan', vid) == NC_NOERR) then
               call gv_dbl(ncid, 'llspan', ncoh, c%leaf_turnover_rate(1:ncoh))
               c%leaf_turnover_rate(1:ncoh) = 1.0_wp / c%leaf_turnover_rate(1:ncoh)
            end if
            !----- Tissue water, temperatures, energies: what a file lacks is set below. --------------!
            tissue_water_ok = nc_inq_varid_f(ncid, 'leaf_water_mass', vid) == NC_NOERR
            tissue_temp_ok  = nc_inq_varid_f(ncid, 'leaf_temp', vid) == NC_NOERR
            tissue_energy_ok = nc_inq_varid_f(ncid, 'leaf_energy', vid) == NC_NOERR
            call gv_dbl_opt(ncid, 'leaf_water_mass', ncoh, c%leaf_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_water_mass', ncoh, c%wood_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_temp',       ncoh, c%leaf_temp(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_temp',       ncoh, c%wood_temp(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_energy',     ncoh, c%leaf_energy(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_energy',     ncoh, c%wood_energy(1:ncoh))
            !----- Stomatal feedback: a file without it keeps alloc_polygon's unset values (seeded   !
            !      from the soil). ---------------------------------------------------------------!
            call gv_dbl_opt(ncid, 'dmax_psi_leaf',   ncoh, c%dmax_psi_leaf(1:ncoh))
            call gv_dbl_opt(ncid, 'dmax_psi_leaf_accum', ncoh, c%dmax_psi_leaf_accum(1:ncoh))
            !----- Phenology memory: a file without it keeps the birth values alloc_polygon set. ---!
            call gv_dbl_opt(ncid, 'leaf_flush_tendency', ncoh, c%leaf_flush_tendency(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_shed_tendency',  ncoh, c%leaf_shed_tendency(1:ncoh))
            call gv_dbl_opt(ncid, 'growing_degree_days', ncoh, c%growing_degree_days(1:ncoh))
            call gv_dbl_opt(ncid, 'cold_degree_days',    ncoh, c%cold_degree_days(1:ncoh))
            call gv_dbl_opt(ncid, 'dry_psi_sum',         ncoh, c%dry_psi_sum(1:ncoh))
            call gv_dbl_opt(ncid, 'wet_psi_sum',         ncoh, c%wet_psi_sum(1:ncoh))
            call gv_dbl_opt(ncid, 'light_hours_mean', ncoh, c%light_hours_mean(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_surf_water',    ncoh, c%leaf_surf_water(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_surf_water',    ncoh, c%wood_surf_water(1:ncoh))
            !----- A file whose films carry no units holds them per m2 of ground. -------------------!
            if (nc_inq_varid_f(ncid, 'leaf_surf_water', vid) == NC_NOERR) then
               if (nc_get_att_text_f(ncid, vid, 'units', units_attr) /= NC_NOERR) units_attr = ''
               if (trim(units_attr) /= FILM_UNITS) then
                  c%leaf_surf_water(1:ncoh) = c%leaf_surf_water(1:ncoh) / c%nplant(1:ncoh)
                  c%wood_surf_water(1:ncoh) = c%wood_surf_water(1:ncoh) / c%nplant(1:ncoh)
               end if
            end if
         end associate
         !----- The carbon pools, and the geometry derived from them as the run derives it; a file !
         !      without pools places each cohort on the allometry from its dbh. ------------------!
         pools_ok = nc_inq_varid_f(ncid, 'leaf_carbon', vid) == NC_NOERR
         if (pools_ok) then
            associate (c => poly%cohort)
               call gv_dbl(ncid, 'leaf_carbon',          ncoh, c%leaf_carbon(1:ncoh))
               call gv_dbl(ncid, 'fineroot_carbon',      ncoh, c%fineroot_carbon(1:ncoh))
               call gv_dbl(ncid, 'wood_carbon',          ncoh, c%wood_carbon(1:ncoh))
               call gv_dbl(ncid, 'nonstructural_carbon', ncoh, c%nonstructural_carbon(1:ncoh))
            end associate
            do i = 1_ik, ncoh
               call set_cohort_size_from_carbon(poly%cohort, i)
            end do
         else
            do i = 1_ik, ncoh
               call set_cohort_size(poly%cohort, i)          ! pools from dbh, then the geometry
            end do
         end if
         olai_ok = nc_inq_varid_f(ncid, 'overtopping_lai', vid) == NC_NOERR
         if (olai_ok) call gv_dbl(ncid, 'overtopping_lai', ncoh, poly%cohort%overtopping_lai(1:ncoh))
         !----- The growth-rate ring buffer, when the file has one of this run's window length; !
         !      otherwise it is seeded full of the saved growth_avg. ----------------------------!
         nwin = growth_window_steps(cfg)
         ghist_ok = nc_inq_varid_f(ncid, 'growth_hist', vid) == NC_NOERR
         if (ghist_ok) then
            call nc_check(nc_inq_dimlen_f(ncid, 'growth_window', nwin_file), 'inq growth_window')
            ghist_ok = int(nwin_file, ik) == nwin .and. int(size(poly%cohort%growth_hist, 1), ik) == nwin
         end if
         if (ghist_ok) then
            call gv_dbl(ncid, 'growth_accum', ncoh, poly%cohort%growth_accum(1:ncoh))
            call gv_int(ncid, 'growth_count', ncoh, poly%cohort%growth_count(1:ncoh))
            do i = 1_ik, ncoh
               call gv_dbl2_row(ncid, 'growth_hist', i, nwin, poly%cohort%growth_hist(1:nwin, i))
            end do
            block
               integer(ik) :: gpos(1)
               call gv_int(ncid, 'growth_hist_pos', 1_ik, gpos)
               poly%growth_hist_pos = gpos(1)
            end block
         else
            do i = 1_ik, ncoh
               if (.not. is_unset_positive(poly%cohort%growth_avg(i))) then
                  poly%cohort%growth_count(i)   = nwin
                  poly%cohort%growth_accum(i)   = poly%cohort%growth_avg(i) * real(nwin, wp)
                  poly%cohort%growth_hist(:, i) = poly%cohort%growth_avg(i)
               end if
            end do
         end if
      end if
      if (npat > 0_ik) then
         associate (p => poly%patch)
            call gv_dbl (ncid, 'patch_area',      npat, p%area(1:npat))
            !----- Each patch's site; a file without them has one site, alloc_polygon's default. ---!
            if (nc_inq_varid_f(ncid, 'owner_site', vid) == NC_NOERR)                                 &
               call gv_int(ncid, 'owner_site', npat, p%owner_site(1:npat))
            call gv_dbl (ncid, 'patch_age',       npat, p%age(1:npat))
            call gv_int (ncid, 'dist_type',       npat, p%dist_type(1:npat))
            call gv_int (ncid, 'global_patch_id', npat, p%global_id(1:npat))
            call gv_dbl_opt(ncid, 'slow_co2_rate', npat, p%slow_co2_rate(1:npat))
            do ip = 1_ik, npat
               call gv_dbl2_row(ncid, 'recruit_pool', ip, npft, p%recruit_pool(1:npft, ip))
            end do
            !----- The fast reservoirs, written together, so the first stands for all: a file without !
            !      them leaves alloc_polygon's defaults for the caller to re-seed. -------------------!
            fast_ok = nc_inq_varid_f(ncid, 'cas_can_enthalpy', vid) == NC_NOERR
            if (fast_ok) then
               block
                  real(wp) :: fc(npat), fw1(npat), fw2(npat), fw3(npat), fw4(npat), fsh(npat)
                  logical  :: pond_enth_ok, shed_ok
                  integer(ik) :: fnl(npat)
                  call gv_dbl(ncid, 'cas_can_enthalpy', npat, fc)
                  do ip = 1_ik, npat
                     p%cas(ip)%can_enthalpy = fc(ip)
                  end do
                  call gv_dbl(ncid, 'cas_can_shv',       npat, fc)
                  do ip = 1_ik, npat
                     p%cas(ip)%can_shv      = fc(ip)
                  end do
                  call gv_dbl(ncid, 'cas_can_co2',       npat, fc)
                  do ip = 1_ik, npat
                     p%cas(ip)%can_co2      = fc(ip)
                  end do
                  call gv_dbl(ncid, 'cas_can_temp',      npat, fc)
                  do ip = 1_ik, npat
                     p%cas(ip)%can_temp     = fc(ip)
                  end do
                  call gv_dbl(ncid, 'soil_w_surface',    npat, fw1)
                  do ip = 1_ik, npat
                     p%soil_water(ip)%w_surface = fw1(ip)
                  end do
                  !----- Pond enthalpy, if the file has it; otherwise it is rebuilt below, from soil_temp. -!
                  pond_enth_ok = nc_inq_varid_f(ncid, 'soil_w_surface_enth', vid) == NC_NOERR
                  fw4 = 0.0_wp
                  if (pond_enth_ok) call gv_dbl(ncid, 'soil_w_surface_enth', npat, fw4)
                  !----- Shed water: a file without it keeps zero until the slow step sets it. -----------!
                  shed_ok = nc_inq_varid_f(ncid, 'shed_water_rate', vid) == NC_NOERR
                  fsh = 0.0_wp
                  if (shed_ok) call gv_dbl(ncid, 'shed_water_rate', npat, fsh)
                  if (shed_ok) then
                     do ip = 1_ik, npat
                        p%shed_water_rate(ip) = fsh(ip)
                     end do
                  end if
                  !----- Its enthalpy: a file without it values the rate at the canopy air's temperature. -!
                  if (nc_inq_varid_f(ncid, 'shed_water_enth', vid) == NC_NOERR) then
                     call gv_dbl(ncid, 'shed_water_enth', npat, p%shed_water_enth(1:npat))
                  else
                     do ip = 1_ik, npat
                        p%shed_water_enth(ip) = p%shed_water_rate(ip) * internal_energy_liquid(p%cas(ip)%can_temp)
                     end do
                  end if
                  !----- The controller's warm start: a file without it cold-starts it. ----------------!
                  call gv_dbl_opt(ncid, 'adapt_dt_last', npat, p%adapt_dt_last(1:npat))
                  call gv_int(ncid, 'snow_nlayer', npat, fnl)
                  do ip = 1_ik, npat
                     p%snow(ip)%nlayer = fnl(ip)
                  end do
                  do ip = 1_ik, npat
                     call gv_dbl2_row(ncid, 'soil_theta',  ip, n_soil_layer_max, p%soil_water(ip)%theta(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_energy', ip, n_soil_layer_max, p%soil_energy(ip)%soil_energy(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_temp',   ip, n_soil_layer_max, p%soil_energy(ip)%soil_temp(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_fliq',   ip, n_soil_layer_max, p%soil_energy(ip)%soil_fliq(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'snow_swe',    ip, n_snow_layer_max, p%snow(ip)%swe(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_energy', ip, n_snow_layer_max, p%snow(ip)%snow_energy(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_depth',  ip, n_snow_layer_max, p%snow(ip)%snow_depth(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_temp',   ip, n_snow_layer_max, p%snow(ip)%snow_temp(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_fliq',   ip, n_snow_layer_max, p%snow(ip)%snow_fliq(1:n_snow_layer_max))
                  end do
                  !----- Pond enthalpy: as read, or the pond's mass at the top soil layer's temperature   !
                  !      (docs/science/order_of_processes.md §6, The restart file). ---------------------!
                  do ip = 1_ik, npat
                     if (pond_enth_ok) then
                        p%soil_water(ip)%w_surface_energy = fw4(ip)
                     else
                        p%soil_water(ip)%w_surface_energy = p%soil_water(ip)%w_surface                            &
                              * internal_energy_liquid(p%soil_energy(ip)%soil_temp(1))
                     end if
                  end do
               end block
               if (present(fast_found)) fast_found = .true.
            end if
            !----- Soil carbon: a file without it restarts the pools at zero. ----------------------!
            block
               real(wp) :: sc(npat, 9_ik)
               sc = 0.0_wp
               call gv_dbl_opt(ncid, 'soilc_fast_grnd',   npat, sc(:,1))
               call gv_dbl_opt(ncid, 'soilc_fast_soil',   npat, sc(:,2))
               call gv_dbl_opt(ncid, 'soilc_struct_grnd', npat, sc(:,3))
               call gv_dbl_opt(ncid, 'soilc_struct_soil', npat, sc(:,4))
               call gv_dbl_opt(ncid, 'soilc_microbial',   npat, sc(:,5))
               call gv_dbl_opt(ncid, 'soilc_slow',        npat, sc(:,6))
               call gv_dbl_opt(ncid, 'soilc_passive',     npat, sc(:,7))
               call gv_dbl_opt(ncid, 'soilc_lignin_grnd', npat, sc(:,8))
               call gv_dbl_opt(ncid, 'soilc_lignin_soil', npat, sc(:,9))
               do ip = 1_ik, npat
                  p%soil_carbon(ip)%fast_grnd_carbon   = sc(ip,1)
                  p%soil_carbon(ip)%fast_soil_carbon   = sc(ip,2)
                  p%soil_carbon(ip)%struct_grnd_carbon = sc(ip,3)
                  p%soil_carbon(ip)%struct_soil_carbon = sc(ip,4)
                  p%soil_carbon(ip)%microbial_carbon   = sc(ip,5)
                  p%soil_carbon(ip)%slow_carbon        = sc(ip,6)
                  p%soil_carbon(ip)%passive_carbon     = sc(ip,7)
                  p%soil_carbon(ip)%struct_grnd_lignin = sc(ip,8)
                  p%soil_carbon(ip)%struct_soil_lignin = sc(ip,9)
               end do
            end block
         end associate
      end if
      !----- A file without the tissue water or temperatures starts each cohort in its patch's state,  !
      !      as a recruit is born: the canopy air and top soil layer just read, or, without those, the !
      !      run's starting state, from which init_fast_reservoirs re-seeds the patch and cohorts. ----!
      if (ncoh > 0_ik .and. npat > 0_ik .and. .not. (tissue_water_ok .and. tissue_temp_ok)) then
         block
            real(wp) :: t_birth, psi_birth
            integer(ik) :: ip
            do i = 1_ik, ncoh
               ip = poly%cohort%owner_patch(i)
               call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip), cfg,                &
                                         poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
               if (.not. tissue_temp_ok) then
                  poly%cohort%leaf_temp(i) = t_birth
                  poly%cohort%wood_temp(i) = t_birth
               end if
               if (.not. tissue_water_ok) then
                  associate (hydro => cfg%pft%hydro(poly%cohort%pft(i)))
                     call set_cohort_tissue_water(poly%cohort, i, psi_birth, hydro%leaf_curve, hydro%wood_curve)
                  end associate
               end if
            end do
         end block
      end if
      !----- The tissue energies are the state and the temperatures their twins. A file without the  !
      !      energies, or whose tissue water was set above, derives them from the temperatures at the  !
      !      masses; one without only the temperatures derives those from the energies. --------------!
      if (ncoh > 0_ik .and. .not. (tissue_energy_ok .and. tissue_water_ok .and. tissue_temp_ok)) then
         do i = 1_ik, ncoh
            if (tissue_energy_ok .and. tissue_water_ok) then
               call set_cohort_tissue_temps(poly%cohort, i)
            else
               call hold_cohort_tissue_temps(poly%cohort, i)
            end if
         end do
      end if
      call nc_check(nc_close(ncid), 'state nc_close (read)')

      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
      !----- A file without overtopping_lai: compute it from the stand as read, so the traits'     !
      !      light environment is the stand's from the first step rather than an open sky. --------!
      if (ncoh > 0_ik .and. .not. olai_ok) call update_overtopping_lai(poly)
      found = .true.

   contains
      !----- A per-polygon count (n_site), at the polygon's row. ----------------------------------!
      integer(ik) function polygon_count(v) result(n)
         integer(c_int), intent(in) :: v
         integer(c_size_t) :: st1(1), ct1(1)
         integer(ik)       :: x(1)
         st1 = int(slot, c_size_t)
         ct1 = 1_c_size_t
         call nc_check(nc_get_vara_int(ncid, v, st1, ct1, x), 'get polygon count')
         n = x(1)
      end function polygon_count
      !----- Where a read starts: the polygon's row in a file with that axis, then `row` (1-based) of a !
      !      2-D variable when given, then the first n entries. ------------------------------------!
      subroutine where_from(n, row, st, ct, nd)
         integer(ik),           intent(in)  :: n
         integer(ik), optional, intent(in)  :: row
         integer(c_size_t),     intent(out) :: st(3), ct(3)
         integer,               intent(out) :: nd
         nd = 0
         if (region) then
            nd = 1
            st(1) = int(slot, c_size_t)
            ct(1) = 1_c_size_t
         end if
         if (present(row)) then
            nd = nd + 1
            st(nd) = int(row - 1_ik, c_size_t)
            ct(nd) = 1_c_size_t
         end if
         nd = nd + 1
         st(nd) = 0_c_size_t
         ct(nd) = int(n, c_size_t)
      end subroutine where_from
      subroutine gv_int(nc, name, n, out)        ! read [1:n] of an int variable by name
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: n
         integer(ik),      intent(out) :: out(:)
         integer(c_int)    :: v
         integer(c_size_t) :: st(3), ct(3)
         integer           :: nd
         call where_from(n, st=st, ct=ct, nd=nd)
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_int(nc, v, st(1:nd), ct(1:nd), out), 'get '//name)
      end subroutine gv_int
      subroutine gv_dbl(nc, name, n, out)        ! read [1:n] of a double variable by name
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: n
         real(wp),         intent(out) :: out(:)
         integer(c_int)    :: v
         integer(c_size_t) :: st(3), ct(3)
         integer           :: nd
         call where_from(n, st=st, ct=ct, nd=nd)
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_double(nc, v, st(1:nd), ct(1:nd), out), 'get '//name)
      end subroutine gv_dbl
      subroutine gv_dbl_opt(nc, name, n, out)    ! read if the var exists; leave `out` unchanged otherwise
         integer(c_int),   intent(in)    :: nc
         character(len=*), intent(in)    :: name
         integer(ik),      intent(in)    :: n
         real(wp),         intent(inout) :: out(:)
         integer(c_int)    :: v
         integer(c_size_t) :: st(3), ct(3)
         integer           :: nd
         if (nc_inq_varid_f(nc, name, v) /= NC_NOERR) return
         call where_from(n, st=st, ct=ct, nd=nd)
         call nc_check(nc_get_vara_double(nc, v, st(1:nd), ct(1:nd), out), 'get '//name)
      end subroutine gv_dbl_opt
      subroutine gv_dbl2_row(nc, name, row, nlayer, out)   ! read row [row,1:nlayer] of a [patch,layer] var
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: row, nlayer
         real(wp),         intent(out) :: out(:)
         integer(c_int)    :: v
         integer(c_size_t) :: st(3), ct(3)
         integer           :: nd
         call where_from(nlayer, row, st, ct, nd)
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_double(nc, v, st(1:nd), ct(1:nd), out), 'get '//name)
      end subroutine gv_dbl2_row
   end subroutine read_restart_file

end module meds_init_restart
