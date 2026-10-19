! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_io -- the STATE (restart) stream: write the full prognostic demographic state to a     !
! netCDF-4 checkpoint and read it back, via the netCDF C library (meds_netcdf_c).             !
!                                                                                          !
! RESTART ONLY. Diagnostic output is the [output] aggregation subsystem (meds_output_*), which !
! is registry-driven -- adding a variable never means editing this file. Do not add a second     !
! diagnostic writer here: it would collide with [output] on the `-D-` filename prefix.           !
!                                                                                          !
! The restart stream stays deliberately orthogonal to the diagnostic streams: a checkpoint must  !
! be the raw prognostic state at an INSTANT, never a time-average, and a run resumed from it     !
! must continue exactly as the run that wrote it would have. So the cohort carbon pools are       !
! stored as the running model holds them, not re-derived from dbh: fusion leaves a cohort off the !
! allometry, so re-deriving them would move the leaf area and reset the fine-root carbon. The     !
! geometry is a function of the pools, and the reader re-derives it from them.                    !
!                                                                                          !
! A file holds every polygon of the run (docs/science/order_of_processes.md §6): each      !
! variable has a leading `polygon` axis, and the cohort, patch and site axes are the       !
! largest any polygon has. n_cohort, n_patch, n_site and n_window mark each polygon's live !
! part (the first two also in its row of meta_int); the rest holds the fill value.         !
! polygon_id, lat and lon name the rows. Each site's area and surface, and each patch's    !
! site, are stored; a file without them is one site, its area meta_real(1). It is written  !
! by state_file_open / state_file_put / state_file_close and read by io_read_state, a      !
! polygon's row at a time; io_read_state also reads a file without the axis, one site's,   !
! as every one-polygon run wrote one before runs became regions.                           !
!==========================================================================================!
module meds_io
   use iso_c_binding, only : c_int, c_size_t, c_double
   use meds_kinds,    only : wp, ik
   use meds_state_defaults, only : is_unset_positive
   use meds_config,   only : meds_config_t, growth_window_steps, pft_tissue_curves
   use meds_time,     only : meds_time_t, time_to_string, time_to_stamp, time_to_decimal_year
   use meds_netcdf_c
   use meds_polygon_state_vars, only : polygon_t
   use meds_patch_state_vars, only : patch_physical_state
   use meds_cohort_state_vars,      only : set_cohort_size, set_cohort_size_from_carbon, set_cohort_tissue_water
   use meds_site_state_vars,        only : make_site_surface
   use meds_polygon_state_vars,     only : alloc_polygon, rebuild_patch_map, rebuild_cohort_map
   use meds_demography_update,      only : update_overtopping_lai
   use meds_constants,     only : n_soil_layer_max, n_snow_layer_max
   use meds_water_retention, only : water_curve_t
   use meds_therm_lib,              only : internal_energy_liquid
   use meds_demography_cohort_fusefiss,   only : sort_cohorts
   implicit none
   private

   public :: io_read_state
   public :: state_file_t, state_file_open, state_file_put, state_file_close, state_file_polygons
   public :: state_file_refusal

   character(len=*), parameter :: STATE_TITLE = 'MEDS instantaneous state (restart) file'
   integer,          parameter :: N_META_INT = 11, N_META_REAL = 2

   !----- An open state file, its polygon axis ahead of every other. --------------------------------!
   type :: state_file_t
      integer(c_int)     :: ncid   = -1_c_int
      character(len=512) :: path   = ''
   end type state_file_t

contains

   !=======================================================================================!
   !  STATE (restart) files                                                                 !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Create <dir>/<prefix>-S-YYYYMMDDHHMMSS.nc and define every dimension and variable. The    !
   ! cohort and patch axes are n_cohort and n_patch long (at least 1), the largest of the       !
   ! polygons'. A `polygon` axis leads every variable, `polygon_id` names its rows, and the       !
   ! variables are compressed, since most polygons fill only part of the padded axes.            !
   !---------------------------------------------------------------------------------------!
   subroutine state_file_open(sf, dir, prefix, now, n_pft, n_window, n_cohort, n_patch, n_site,     &
                              polygon_id, polygon_lat, polygon_lon, restructure_pending, restructure_new_year)
      type(state_file_t),    intent(out) :: sf
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

      !----- The timestamp IS the simulated calendar date -> <prefix>-S-YYYYMMDDHHMMSS.nc. --!
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
      !----- The sites: each one's area and the surface it was made with. The surface is made again   !
      !      from the configuration on read; its layer count must match, or the soil state would be   !
      !      read onto another column. ---------------------------------------------------------------!
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
      !----- Per-cohort plant-hydraulics internal water + tissue temperatures. Persisting these makes  !
      !      a restart EXACT for the fast loop -- without them the restart starts the tissues in their   !
      !      patch's state (patch_physical_state), whose step from the evolved state makes the             !
      !      plant-hydraulics solve grind for a few days on a healthy high-LAI restart. OPTIONAL on read !
      !      (gv_dbl_opt), so a state file without them still restarts, on those values. The            !
      !      leaf-surface (interception) film is stored below with the pools. ---------------------------!
      call dv('leaf_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal leaf water [kg/plant]')
      call dv('wood_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal wood water [kg/plant]')
      call dv('leaf_temp',        NC_DOUBLE, [d_cohort], 'per-cohort leaf temperature [K]')
      call dv('wood_temp',        NC_DOUBLE, [d_cohort], 'per-cohort wood temperature [K]')
      !----- Per-cohort stomatal water-stress feedback (docs/science/leaf_gas_exchange.md §4.2).      !
      !      `dmax_psi_leaf` is YESTERDAY's completed daily-max leaf water potential and drives TODAY's !
      !      beta_stomata; `dmax_psi_leaf_accum` is the running accumulator the fast loop fills,        !
      !      rolled over at day end. Both must persist:                                                 !
      !                                                                                          !
      !        * without predawn, a restart falls back to the UNSET sentinel and re-seeds from surface-   !
      !          layer soil psi -- which reads as unstressed, so a restart mid-drought silently re-opens  !
      !          stomata that the running model had closed.                                                !
      !        * without the accumulator, a restart at any time other than midnight discards the part of  !
      !          today's maximum already seen, so the value handed to tomorrow is biased LOW (too         !
      !          negative => too much stress) by however much of the day was replayed.                     !
      !                                                                                          !
      !      OPTIONAL on read like the block above, so a state file without them restarts on the        !
      !      sentinels. NOT derivable from leaf_water_mass: predawn is a time-MAXIMUM over the previous   !
      !      day, not a function of the instantaneous state. ---------------------------------------------!
      call dv('dmax_psi_leaf',   NC_DOUBLE, [d_cohort], &
            'yesterday daily-max leaf water potential [MPa] (drives beta_stomata)')
      call dv('dmax_psi_leaf_accum', NC_DOUBLE, [d_cohort], 'running daily-max leaf water potential accumulator [MPa]')
      !----- PHENOLOGY MEMORY. Without it a restart would resurrect every cohort at its birth state  !
      !      (flushing, no senescence, empty sums): a deciduous stand restarted in January would leaf  !
      !      out in midwinter. None of it is derivable from the rest of the file: these are time      !
      !      integrals over the preceding season (degree-day sums, psi sums, a running mean), not     !
      !      functions of the instantaneous state. OPTIONAL on read, so a file without them restarts  !
      !      on the birth values.  ---------------------------------------------------------------------!
      call dv('leaf_flush_tendency', NC_DOUBLE, [d_cohort], 'smoothed flush signal [0,1]')
      call dv('leaf_shed_tendency',  NC_DOUBLE, [d_cohort], 'smoothed senescence signal [0,1]')
      call dv('growing_degree_days', NC_DOUBLE, [d_cohort], 'warmth sum since midwinter [K day]')
      call dv('cold_degree_days',    NC_DOUBLE, [d_cohort], 'cold sum since midsummer [K day]')
      call dv('dry_psi_sum',         NC_DOUBLE, [d_cohort], 'predawn leaf psi below the TLP, summed [MPa day]')
      call dv('wet_psi_sum',         NC_DOUBLE, [d_cohort], 'predawn leaf psi above the TLP, summed [MPa day]')
      call dv('light_hours_mean', NC_DOUBLE, [d_cohort], 'running-mean hours of light at the cohort top [h/day]')
      !----- The carbon pools as the running model holds them. A cohort's pools sit off the        !
      !      allometry for its dbh (fusion, phenology), so re-deriving them from dbh on read would   !
      !      move the leaf area (BCI: LAI 5.6380 written would read back as 5.6389). The reader      !
      !      derives the geometry from the pools, as the run does; height/basal_area/agb/leaf_area   !
      !      are written for the file's readers (the calibration's stand check). overtopping_lai is  !
      !      what the plastic traits were acclimated to. OPTIONAL on read (a file without the pools  !
      !      re-derives them from dbh). ----------------------------------------------------------!
      call dv('height',               NC_DOUBLE, [d_cohort], 'height [m]')
      call dv('basal_area',           NC_DOUBLE, [d_cohort], 'basal area per plant [cm2]')
      call dv('agb',                  NC_DOUBLE, [d_cohort], 'aboveground biomass per plant [kgC]')
      call dv('leaf_area',            NC_DOUBLE, [d_cohort], 'leaf area per plant [m2]')
      call dv('leaf_carbon',          NC_DOUBLE, [d_cohort], 'leaf carbon per plant [kgC]')
      call dv('fineroot_carbon',      NC_DOUBLE, [d_cohort], 'fine-root carbon per plant [kgC]')
      call dv('wood_carbon',          NC_DOUBLE, [d_cohort], 'wood carbon per plant [kgC]')
      call dv('nonstructural_carbon', NC_DOUBLE, [d_cohort], 'nonstructural carbon per plant [kgC]')
      call dv('overtopping_lai',      NC_DOUBLE, [d_cohort], 'LAI of the taller cohorts in the patch [m2/m2]')
      !----- The canopy interception film ([fast].canopy_water_on). A dry canopy is a valid start,   !
      !      but not the state the writing run was in. OPTIONAL on read (a file without it starts   !
      !      dry). ---------------------------------------------------------------------------------!
      call dv('leaf_surf_water',      NC_DOUBLE, [d_cohort], 'leaf interception film [kg/m2 ground]')
      call dv('wood_surf_water',      NC_DOUBLE, [d_cohort], 'wood interception film [kg/m2 ground]')
      !----- The growth-rate ring buffer behind growth_avg (the mortality predictor) and its site-wide !
      !      write slot. OPTIONAL on read: a file without it, or one written with another window       !
      !      length, reseeds the buffer from growth_avg. ---------------------------------------------!
      call dv('growth_accum',         NC_DOUBLE, [d_cohort], 'running sum of the growth window [cm/yr]')
      call dv('growth_count',         NC_INT,    [d_cohort], 'samples in the growth window')
      call dv('growth_hist',          NC_DOUBLE, [d_cohort, d_gwin], 'growth-rate ring buffer [cm/yr]')
      call dv('growth_hist_pos',      NC_INT,    [d_mr], 'ring-buffer write slot (element 1)')
      !----- THERMAL ACCLIMATION (docs/science/leaf_gas_exchange.md §2, Thermal acclimation). A SITE !
      !      scalar, but written as its own variable rather than appended to meta_real, because      !
      !      meta_real is read at a fixed length of 2 and growing it would stop older state files    !
      !      loading. Optional on read, so a file without it restarts on the unseeded sentinel --    !
      !      which re-adopts the first day's mean rather than relaxing from an arbitrary origin over !
      !      a month.  ------------------------------------------------------------------------------!
      call dv('t_growth_avg', NC_DOUBLE, [d_mr], 'running-mean growth temperature [K]')
      call dv('patch_area',       NC_DOUBLE, [d_patch],  'patch area fraction of its site')
      call dv('owner_site',       NC_INT,    [d_patch],  'owning site index (1-based)')
      call dv('patch_age',        NC_DOUBLE, [d_patch],  'time since last disturbance [yr]')
      call dv('dist_type',        NC_INT,    [d_patch],  'disturbance type (1=primary,2=treefall)')
      call dv('global_patch_id',  NC_INT,    [d_patch],  'persistent patch id')
      call dv('recruit_pool',     NC_DOUBLE, [d_patch, d_pft], 'carry-forward recruit pool [plant/m2]')
      call dv('slow_co2_rate',    NC_DOUBLE, [d_patch], 'slow-loop CO2 flux the fast NEE adds [umol/m2/s]')
      !----- FAST reservoirs: the true evolved CAS/soil/snow state, not just its cached geometry. A     !
      !      restart that re-seeds these generically hands an already-mature canopy a discontinuous jump  !
      !      on day 1, which the adaptive integrators resolve only by grinding through many tiny          !
      !      sub-steps (RK45's explicit stepper worst of all). Diagnosed twins (soil_temp/fliq,           !
      !      snow_temp/fliq, cas_temp) are saved directly alongside their prognostic partner rather than   !
      !      recomputed on read, so this reader needs no outside PFT/soil config (matching every other     !
      !      read here). Always written; reading them back is optional as a whole, so a state file         !
      !      without them falls back to the generic seed (meds_polygon, init_fast_reservoirs) with no      !
      !      format-version bump. --------------------------------------------------------------------------!
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
      !----- Daily leaf/root-turnover shed water handed to the ground: a slow->fast seam rate the    !
      !      fast loop reads all day, so a restart without it would lose that day's shed water. -----!
      call dv('shed_water_rate',  NC_DOUBLE, [d_patch], 'turnover shed water to the ground [kg/m2/s]')
      !----- The fast integrators' warm start, per patch: the step the adaptive controller last        !
      !      accepted, where the next fast step starts. Without it a resumed run cold-starts the        !
      !      controller and takes different sub-steps from the continuous run. -------------------------!
      call dv('adapt_dt_last',    NC_DOUBLE, [d_patch], 'fast-controller warm-start step [s] (0 = cold)')
      !----- Slow soil-carbon pools ([soil_carbon].soil_carbon_on, on by default;                      !
      !      docs/science/soil_carbon.md §2). There are no nitrogen fields: the nitrogen cycle is not   !
      !      built (see docs/ROADMAP.md). ----------------------------------------------------------------!
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
      !----- Define a state variable on its site dimensions, behind the polygon axis, with its     !
      !      long_name. ----------------------------------------------------------------------------!
      subroutine dv(name, xtype, dimids, lname)
         character(len=*), intent(in) :: name, lname
         integer(c_int),   intent(in) :: xtype, dimids(:)
         integer(c_int)    :: v, dims(size(dimids) + 1)
         integer(c_size_t) :: chunk(size(dimids) + 1)
         integer           :: nd, k
         dims(1) = d_poly
         dims(2:) = dimids
         nd = size(dimids) + 1
         call nc_check(nc_def_var_f(ncid, name, xtype, int(nd, c_int), dims(1:nd), v), 'state def '//name)
         call nc_check(nc_put_att_text_f(ncid, v, 'long_name',                              &
                       int(len_trim(lname), c_size_t), lname), 'state long_name '//name)
         !----- One polygon's row per chunk. Each polygon writes its own rows (state_file_put), so a   !
         !      chunk shared by many polygons is decompressed and recompressed at each of their writes   !
         !      once the chunks of a row outgrow HDF5's chunk cache: a 606-polygon region (example 06)   !
         !      then takes hours per checkpoint instead of minutes. A polygon-only variable is one       !
         !      small chunk, which the cache holds. -----------------------------------------------------!
         if (nd > 1) then
            chunk(1) = 1_c_size_t
            do k = 2, nd
               call nc_check(nc_inq_dimlen(ncid, dims(k), chunk(k)), 'state dimlen '//name)
               chunk(k) = max(chunk(k), 1_c_size_t)
            end do
            call nc_check(nc_def_var_chunking(ncid, v, NC_CHUNKED, chunk(1:nd)), 'state chunk '//name)
         end if
         call nc_check(nc_def_var_deflate(ncid, v, 1_c_int, 1_c_int, 1_c_int), 'state deflate '//name)
      end subroutine dv
   end subroutine state_file_open

   !---------------------------------------------------------------------------------------!
   ! Write one polygon's state into its row (1-based) of an open file.                         !
   !---------------------------------------------------------------------------------------!
   subroutine state_file_put(sf, poly, now, polygon)
      type(state_file_t),    intent(in) :: sf
      type(polygon_t),       intent(in) :: poly
      type(meds_time_t),     intent(in) :: now
      integer(ik),           intent(in) :: polygon
      integer(ik) :: ncoh, npat, npft, ip, ic, nwin, slot, ns, meta_i(N_META_INT)
      real(wp)    :: meta_r(N_META_REAL)

      ncoh = poly%cohort%n ; npat = poly%patch%n ; npft = poly%n_pft
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
      call put_dbl('t_growth_avg', [poly%t_growth_avg])
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
            block                                 ! FAST reservoirs: contiguous buffers, same
               real(wp) :: fc(npat, 4_ik)          ! discipline as the soil-carbon block below (never a
               real(wp) :: fw(npat, 4_ik)          ! bare derived-type component-array-section straight
               real(wp) :: fsoil(npat, n_soil_layer_max, 4_ik)     ! into a C-bound call)
               real(wp) :: fsnow(npat, n_snow_layer_max, 5_ik)
               integer(ik) :: fnl(npat)
               do ip = 1_ik, npat
                  fc(ip,1) = p%cas(ip)%can_enthalpy ; fc(ip,2) = p%cas(ip)%can_shv
                  fc(ip,3) = p%cas(ip)%can_co2      ; fc(ip,4) = p%cas(ip)%can_temp
                  fw(ip,1) = p%soil_water(ip)%w_surface ; fw(ip,2) = p%soil_water(ip)%w_surface_energy
                  fw(ip,3) = p%shed_water_rate(ip) ; fw(ip,4) = p%adapt_dt_last(ip)
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
            block                              ! explicit contiguous buffer (never a bare derived-type
               real(wp) :: sc(npat, 9_ik)       ! component-array-section straight into a C-bound call)
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
         st = int(slot, c_size_t) ; ct = 1_c_size_t ; x = n
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
         if (present(row)) then ; n = n + 1 ; st(n) = int(row, c_size_t) ; ct(n) = 1_c_size_t ; end if
         n = n + 1 ; st(n) = 0_c_size_t ; ct(n) = int(nx, c_size_t)
      end subroutine where_to
      integer(c_int) function varid(name)
         character(len=*), intent(in) :: name
         call nc_check(nc_inq_varid_f(sf%ncid, name, varid), 'state inq '//name)
      end function varid
   end subroutine state_file_put

   subroutine state_file_close(sf)
      type(state_file_t), intent(inout) :: sf
      call nc_check(nc_close(sf%ncid), 'state nc_close')
      sf%ncid = -1_c_int
      write(*,'(2a)') ' state : ', trim(sf%path)
   end subroutine state_file_close

   !----- A region state file's polygon ids, in its rows' order; found = .false. for a file      !
   !      without them (a site's) or one that cannot be opened. ---------------------------------!
   subroutine state_file_polygons(path, polygon_id, found)
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
         start1 = 0_c_size_t ; count1 = n
         st = nc_get_vara_int(ncid, vid, start1, count1, polygon_id)
      end if
      found = st == NC_NOERR
      st = nc_close(ncid)
   end subroutine state_file_polygons

   !---------------------------------------------------------------------------------------!
   ! Why the state file at `path` cannot start a run of `n_polygon` polygons, or '' when it     !
   ! can: it cannot be opened, or it holds one site without a `polygon` axis (a file written    !
   ! before every run had one) and the run has more than one polygon.                           !
   !---------------------------------------------------------------------------------------!
   function state_file_refusal(path, n_polygon) result(why)
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
   end function state_file_refusal

   !----- Does the open state file hold a region's state, a `polygon` axis ahead of every other? ---!
   logical function has_polygon_axis(ncid)
      integer(c_int), intent(in) :: ncid
      integer(c_size_t) :: n
      has_polygon_axis = nc_inq_dimlen_f(ncid, 'polygon', n) == NC_NOERR
   end function has_polygon_axis

   !---------------------------------------------------------------------------------------!
   ! Reconstruct a polygon from a state file written by state_file_put: with `polygon`        !
   ! (1-based) that row of a file with the polygon axis, or the whole of an older file without !
   ! it. Each site's surface is made from `site_cfg` (the run's `cfg` for every site when       !
   ! absent). The geometry is re-derived from the carbon pools (set_cohort_size_from_carbon);   !
   ! the maps and within-patch sort order are rebuilt, and the restructuring the stand owes at  !
   ! restart_time is set on the polygon (state_file_open; a file without the attribute is an    !
   ! older format, written AFTER its boundary's restructuring, so it owes nothing).             !
   ! found=.false. (no error stop) if the file cannot be opened or its layout is not the one    !
   ! `polygon` asks for (state_file_refusal says which). Errors out only on a genuine           !
   ! inconsistency: the PFT count, the site count or a site's soil layer count differs from the !
   ! configuration's.                                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine io_read_state(poly, cfg, path, restart_time, found, fast_found, polygon, site_cfg)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      type(meds_time_t),   intent(out) :: restart_time
      logical,             intent(out) :: found
      !----- FAST reservoirs: .true. only when THIS file's cas/soil/snow state was actually read     !
      !      back (always written together, so one variable's presence implies all of them); .false.  !
      !      means a file without them -- the caller then re-seeds them generically                    !
      !      (init_fast_reservoirs, from meds_polygon). ----------------------------------------------!
      logical, optional,   intent(out) :: fast_found
      integer(ik), optional, intent(in) :: polygon    !< its row (1-based); absent for a file without the axis
      type(meds_config_t), optional, intent(in) :: site_cfg(:)   !< each site's configuration
      character(len=16) :: pending
      character(len=160) :: msg
      integer(c_int) :: ncid, vid, st
      integer(ik)    :: ncoh, npat, npft, ip, i, nwin, slot, ns, is, meta_i(N_META_INT)
      real(wp)       :: meta_r(N_META_REAL)
      logical        :: fast_ok, pools_ok, olai_ok, ghist_ok, region, tissue_water_ok, tissue_temp_ok
      integer(c_size_t) :: nwin_file

      found = .false. ; restart_time = meds_time_t()
      tissue_water_ok = .true. ; tissue_temp_ok = .true.
      if (present(fast_found)) fast_found = .false.
      region = present(polygon)
      slot = 0_ik
      if (region) slot = polygon - 1_ik
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) return
      !----- A region's variables have a polygon axis ahead of the others, so reading the other     !
      !      kind's file would take the wrong elements. ----------------------------------------------!
      if (has_polygon_axis(ncid) .neqv. region) then
         st = nc_close(ncid)
         return
      end if
      pending = 'none'
      if (nc_get_att_text_f(ncid, NC_GLOBAL, 'restructure_pending', pending) /= NC_NOERR) pending = 'none'

      call gv_int(ncid, 'meta_int',  int(N_META_INT, ik),  meta_i)
      call gv_dbl(ncid, 'meta_real', int(N_META_REAL, ik), meta_r)
      ncoh = meta_i(1) ; npat = meta_i(2) ; npft = meta_i(3)
      if (npft /= cfg%pft%n) error stop 'io_read_state: PFT count in state file /= config PFT count'

      !----- The sites: a file without the site axis holds one. ---------------------------------!
      ns = 1_ik
      if (nc_inq_varid_f(ncid, 'n_site', vid) == NC_NOERR) ns = polygon_count(vid)
      if (present(site_cfg)) then
         if (size(site_cfg, kind=ik) /= ns) then
            write(msg, '(a,i0,a,i0,a)') 'io_read_state: the state file holds ', ns,                     &
                                        ' sites and the configuration ', size(site_cfg), ' for the polygon'
            error stop trim(msg)
         end if
      end if
      call alloc_polygon(poly, cfg%pft%n, max(ncoh, 1_ik), max(npat, 1_ik), growth_window_steps(cfg), &
                      cfg%allom, cfg%pft, n_site=ns)
      do is = 1_ik, ns
         if (present(site_cfg)) then
            call make_site_surface(poly%site%surface(is), site_cfg(is))
         else
            call make_site_surface(poly%site%surface(is), cfg)
         end if
      end do
      poly%cohort%n = ncoh ; poly%patch%n = npat
      poly%restructure_pending  = trim(pending) == 'month' .or. trim(pending) == 'year'
      poly%restructure_new_year = trim(pending) == 'year'
      !----- Growth temperature: OPTIONAL, so a state file without it keeps the unseeded          !
      !      sentinel alloc_polygon set and the first slow step adopts that day's mean.  ------------!
      block
         real(wp) :: tg(1)
         tg(1) = poly%t_growth_avg
         call gv_dbl_opt(ncid, 't_growth_avg', 1_ik, tg)
         poly%t_growth_avg = tg(1)
      end block
      poly%next_cohort_id = meta_i(4) ; poly%next_patch_id = meta_i(5)
      if (nc_inq_varid_f(ncid, 'site_area', vid) == NC_NOERR) then
         call gv_dbl(ncid, 'site_area', ns, poly%site%area(1:ns))
         !----- Each site's soil state is per layer of the column it was written on. -------------!
         block
            integer(ik) :: nsoil(ns)
            call gv_int(ncid, 'site_n_soil', ns, nsoil)
            do is = 1_ik, ns
               if (nsoil(is) /= poly%site%surface(is)%soil%n_active) then
                  write(msg, '(a,i0,a,i0,a,i0,a)') 'io_read_state: site ', is, ' has ', nsoil(is),       &
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
            !----- Plastic leaf traits: default to PFT top-of-canopy, then overwrite from the state  !
            !       file where present (a file without them restarts at top-of-canopy). ----------!
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
            !----- A state file written while this trait was a leaf lifespan [yr] holds that instead, -!
            !      as `llspan`: the turnover rate is its reciprocal. -----------------------------------!
            if (nc_inq_varid_f(ncid, 'llspan', vid) == NC_NOERR) then
               call gv_dbl(ncid, 'llspan', ncoh, c%leaf_turnover_rate(1:ncoh))
               c%leaf_turnover_rate(1:ncoh) = 1.0_wp / c%leaf_turnover_rate(1:ncoh)
            end if
            !----- Per-cohort hydraulics/temperature (OPTIONAL): a state file with them restores the   !
            !      true state; without them each cohort starts in its patch's state, set once the patches !
            !      and the pools are read (below).                                                       !
            tissue_water_ok = nc_inq_varid_f(ncid, 'leaf_water_mass', vid) == NC_NOERR
            tissue_temp_ok  = nc_inq_varid_f(ncid, 'leaf_temp', vid) == NC_NOERR
            call gv_dbl_opt(ncid, 'leaf_water_mass', ncoh, c%leaf_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_water_mass', ncoh, c%wood_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_temp',       ncoh, c%leaf_temp(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_temp',       ncoh, c%wood_temp(1:ncoh))
            !----- Stomatal feedback (OPTIONAL): alloc_polygon has already set UNSET_NEGATIVE /        !
            !      DMAX_PSI_LEAF_ACCUM_RESET, so a state file without them seeds from the soil; one with  !
            !      them restores the true feedback state. ------------------------------------------------!
            call gv_dbl_opt(ncid, 'dmax_psi_leaf',   ncoh, c%dmax_psi_leaf(1:ncoh))
            call gv_dbl_opt(ncid, 'dmax_psi_leaf_accum', ncoh, c%dmax_psi_leaf_accum(1:ncoh))
            !----- Phenology memory (OPTIONAL): a file without it restarts on the birth values      !
            !      alloc_polygon has already set; a file with it restores the season. --------------!
            call gv_dbl_opt(ncid, 'leaf_flush_tendency', ncoh, c%leaf_flush_tendency(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_shed_tendency',  ncoh, c%leaf_shed_tendency(1:ncoh))
            call gv_dbl_opt(ncid, 'growing_degree_days', ncoh, c%growing_degree_days(1:ncoh))
            call gv_dbl_opt(ncid, 'cold_degree_days',    ncoh, c%cold_degree_days(1:ncoh))
            call gv_dbl_opt(ncid, 'dry_psi_sum',         ncoh, c%dry_psi_sum(1:ncoh))
            call gv_dbl_opt(ncid, 'wet_psi_sum',         ncoh, c%wet_psi_sum(1:ncoh))
            call gv_dbl_opt(ncid, 'light_hours_mean', ncoh, c%light_hours_mean(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_surf_water',    ncoh, c%leaf_surf_water(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_surf_water',    ncoh, c%wood_surf_water(1:ncoh))
         end associate
         !----- The carbon pools the writing run held, and the geometry derived from them as the  !
         !      run derives it (the file's height/agb/leaf_area are for its readers, not read back). !
         !      A file without pools places each cohort on the allometry from its dbh. -----------!
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
         !----- The growth-rate ring buffer, when the file has one of this run's window length.  !
         !      Otherwise seed it as if it were full of the saved growth_avg; it is overwritten     !
         !      within a window. -----------------------------------------------------------------!
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
            !----- FAST reservoirs: OPTIONAL as a WHOLE (existence of the first variable stands for     !
            !      all of them, always written together) -- a file without them leaves p%cas/soil_energy/      !
            !      soil_water/snow at alloc_polygon's defaults, and the caller re-seeds them via              !
            !      init_fast_reservoirs. -----------------------------------------------------------------!
            fast_ok = nc_inq_varid_f(ncid, 'cas_can_enthalpy', vid) == NC_NOERR
            if (fast_ok) then
               block
                  real(wp) :: fc(npat), fw1(npat), fw2(npat), fw3(npat), fw4(npat), fsh(npat)
                  logical  :: pond_enth_ok, shed_ok
                  integer(ik) :: fnl(npat)
                  call gv_dbl(ncid, 'cas_can_enthalpy', npat, fc)
                  do ip = 1_ik, npat ; p%cas(ip)%can_enthalpy = fc(ip) ; end do
                  call gv_dbl(ncid, 'cas_can_shv',       npat, fc)
                  do ip = 1_ik, npat ; p%cas(ip)%can_shv      = fc(ip) ; end do
                  call gv_dbl(ncid, 'cas_can_co2',       npat, fc)
                  do ip = 1_ik, npat ; p%cas(ip)%can_co2      = fc(ip) ; end do
                  call gv_dbl(ncid, 'cas_can_temp',      npat, fc)
                  do ip = 1_ik, npat ; p%cas(ip)%can_temp     = fc(ip) ; end do
                  call gv_dbl(ncid, 'soil_w_surface',    npat, fw1)
                  do ip = 1_ik, npat
                     p%soil_water(ip)%w_surface = fw1(ip)
                  end do
                  !----- Pond ENTHALPY: OPTIONAL, because older state files do not carry it. Read the
                  !      value here if present; the RECONSTRUCTION for a file without it waits until
                  !      soil_temp is in (below).
                  pond_enth_ok = nc_inq_varid_f(ncid, 'soil_w_surface_enth', vid) == NC_NOERR
                  fw4 = 0.0_wp
                  if (pond_enth_ok) call gv_dbl(ncid, 'soil_w_surface_enth', npat, fw4)
                  !----- shed_water_rate: OPTIONAL (older files lack it; the slow loop rewrites it daily). -!
                  shed_ok = nc_inq_varid_f(ncid, 'shed_water_rate', vid) == NC_NOERR
                  fsh = 0.0_wp
                  if (shed_ok) call gv_dbl(ncid, 'shed_water_rate', npat, fsh)
                  if (shed_ok) then
                     do ip = 1_ik, npat ; p%shed_water_rate(ip) = fsh(ip) ; end do
                  end if
                  !----- The controller's warm start: OPTIONAL, a file without it cold-starts it. -------!
                  call gv_dbl_opt(ncid, 'adapt_dt_last', npat, p%adapt_dt_last(1:npat))
                  call gv_int(ncid, 'snow_nlayer', npat, fnl)
                  do ip = 1_ik, npat ; p%snow(ip)%nlayer = fnl(ip) ; end do
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
                  !----- Pond ENTHALPY, committed now that soil_temp is in. When the file carries it, use  !
                  !      it verbatim. When it does not, RECONSTRUCT from the pond mass at the top soil      !
                  !      layer's temperature rather than defaulting to 0: such a file records a pond        !
                  !      whose enthalpy is UNKNOWN, not empty, and seeding 0 J/m2 would restart the run     !
                  !      ~1 MJ short per kg of ponded water -- which the energy budgets would               !
                  !      immediately (and correctly) report as a leak. The top layer's temperature is the   !
                  !      one the pond exchanges with, so it is the defensible equilibrium guess. A dry      !
                  !      pond (the common case) reconstructs to exactly 0 either way. --------------------!
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
            !----- Slow soil-carbon pools ([soil_carbon].soil_carbon_on): OPTIONAL variables         !
            !      (gv_dbl_opt), so a state file without them restarts the pools at the just-allocated  !
            !      zero (bare ground), as the plastic leaf traits restart at top-of-canopy. -----------!
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
      !----- A state file without the cohorts' tissue water or temperatures starts each cohort in its  !
      !      patch's state, as a recruit is born (patch_physical_state): the canopy air and top soil layer !
      !      just read, or, for a file without them either, the run's starting state, from which       !
      !      init_fast_reservoirs then re-seeds the patch and its cohorts. -----------------------------!
      if (ncoh > 0_ik .and. npat > 0_ik .and. .not. (tissue_water_ok .and. tissue_temp_ok)) then
         block
            type(water_curve_t) :: leaf_curve, wood_curve
            real(wp) :: t_birth, psi_birth
            integer(ik) :: ip
            do i = 1_ik, ncoh
               ip = poly%cohort%owner_patch(i)
               call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip), cfg,                &
                                         poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
               if (.not. tissue_temp_ok) then
                  poly%cohort%leaf_temp(i) = t_birth ; poly%cohort%wood_temp(i) = t_birth
               end if
               if (.not. tissue_water_ok) then
                  call pft_tissue_curves(cfg%hydraulics, cfg%pft, poly%cohort%pft(i), leaf_curve, wood_curve)
                  call set_cohort_tissue_water(poly%cohort, i, psi_birth, leaf_curve, wood_curve)
               end if
            end do
         end block
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
      !----- Where a read starts: the polygon's row in a region's file, then `row` (1-based) of a   !
      !      2-D variable when given, then the first n entries. ------------------------------------!
      subroutine where_from(n, row, st, ct, nd)
         integer(ik),           intent(in)  :: n
         integer(ik), optional, intent(in)  :: row
         integer(c_size_t),     intent(out) :: st(3), ct(3)
         integer,               intent(out) :: nd
         nd = 0
         if (region) then ; nd = 1 ; st(1) = int(slot, c_size_t) ; ct(1) = 1_c_size_t ; end if
         if (present(row)) then ; nd = nd + 1 ; st(nd) = int(row - 1_ik, c_size_t) ; ct(nd) = 1_c_size_t ; end if
         nd = nd + 1 ; st(nd) = 0_c_size_t ; ct(nd) = int(n, c_size_t)
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
   end subroutine io_read_state

end module meds_io
