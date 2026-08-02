!==========================================================================================!
! meds_io -- write the full demographic state (cohort + patch + site levels) to a netCDF-4   !
! file over time, via the netCDF C library (meds_netcdf_c).                                  !
!                                                                                          !
! Layout (single-site, ragged per record): an UNLIMITED `time` dimension and fixed `cohort`  !
! and `patch` dimensions (sized from cfg caps). Each call to io_write_snapshot appends ONE    !
! time record, writing only the valid prefix [1:n] of each per-cohort / per-patch array and   !
! recording n_cohort(time) / n_patch(time). A reader reconstructs patch->cohort membership    !
! per record from cohort_offset / cohort_count (1-based) within n_patch / n_cohort. Per-record !
! cohort/patch slabs are chunked + deflated so the unused fill tail compresses away.          !
!                                                                                          !
! Writes are INFREQUENT (default once a year, cfg%io_output_interval_years) relative to the   !
! daily kernels, so I/O never throttles compute and no buffering beyond the live state is      !
! needed. The file is opened once (io_create) and closed once (io_close).                     !
!==========================================================================================!
module meds_io
   use iso_c_binding, only : c_int, c_size_t, c_double
   use meds_kinds,    only : wp, ik
   use meds_config,   only : meds_config_t, growth_window_steps
   use meds_time,     only : meds_time_t, time_to_string, time_to_stamp, time_to_decimal_year
   use meds_netcdf_c
   use meds_core_interface,   only : site_t
   use meds_core_state_types,       only : site_alloc, gather_pft_params, set_cohort_size, rebuild_csr
   use meds_column_state_types,     only : n_soil_layer_max, n_snow_layer_max
   use meds_therm_lib,              only : internal_energy_liquid
   use meds_core_cohort_fusefiss,   only : sort_cohorts
   use meds_diagnostic_reduce, only : total_nplant, total_basal_area, total_agb,         &
                                           total_lai, mean_dbh
   implicit none
   private

   public :: meds_io_t, io_create, io_write_snapshot, io_close
   public :: io_write_state, io_read_state

   character(len=*), parameter :: TITLE       = 'MEDS demographic diagnostic output'
   character(len=*), parameter :: STATE_TITLE = 'MEDS instantaneous state (restart) file'

   type :: meds_io_t
      integer(c_int) :: ncid = -1_c_int
      integer(ik)    :: nrec = 0_ik
      integer(ik)    :: cohort_max = 0_ik, patch_max = 0_ik
      integer(c_int) :: d_time, d_cohort, d_patch
      integer(c_int) :: v_time, v_year, v_month, v_day, v_ncohort, v_npatch
      integer(c_int) :: v_c_pft, v_c_nplant, v_c_dbh, v_c_height, v_c_ba, v_c_agb, v_c_la, v_c_owner
      integer(c_int) :: v_c_gid, v_c_gavg
      integer(c_int) :: v_p_area, v_p_age, v_p_dist, v_p_coff, v_p_ccount
      integer(c_int) :: v_p_gid
      integer(c_int) :: v_s_nplant, v_s_ba, v_s_agb, v_s_lai, v_s_dmean
   end type meds_io_t

contains

   !---------------------------------------------------------------------------------------!
   ! Create the file and define all dimensions, variables, attributes (then leave define    !
   ! mode). The cohort/patch dimensions are fixed at the config caps.                        !
   !---------------------------------------------------------------------------------------!
   subroutine io_create(io, path, cfg)
      type(meds_io_t),     intent(out) :: io
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(in)  :: cfg

      io%cohort_max = cfg%io_cohort_max
      io%patch_max  = cfg%io_patch_max
      io%nrec       = 0_ik

      call nc_check(nc_create_f(path, ior(NC_NETCDF4, NC_CLOBBER), io%ncid), 'nc_create')
      call nc_check(nc_def_dim_f(io%ncid, 'time',   NC_UNLIMITED,                  io%d_time),   'dim time')
      call nc_check(nc_def_dim_f(io%ncid, 'cohort', int(io%cohort_max, c_size_t),  io%d_cohort), 'dim cohort')
      call nc_check(nc_def_dim_f(io%ncid, 'patch',  int(io%patch_max,  c_size_t),  io%d_patch),  'dim patch')

      !----- time coordinate (real calendar) + per-record counts. -------------------------!
      call def_t(io%v_time,    'time',     NC_DOUBLE, 'year', 'decimal calendar year')
      call def_t(io%v_year,    'year',     NC_INT,    '--',   'calendar year')
      call def_t(io%v_month,   'month',    NC_INT,    '--',   'calendar month (1-12)')
      call def_t(io%v_day,     'day',      NC_INT,    '--',   'calendar day of month')
      call def_t(io%v_ncohort, 'n_cohort', NC_INT,    '--',   'number of cohorts in use this record')
      call def_t(io%v_npatch,  'n_patch',  NC_INT,    '--',   'number of patches in use this record')

      !----- per-cohort state (time, cohort). ---------------------------------------------!
      call def_c(io%v_c_pft,    'pft',         NC_INT,    '--',        'plant functional type index')
      call def_c(io%v_c_nplant, 'nplant',      NC_DOUBLE, 'plant/m2',  'plant number density')
      call def_c(io%v_c_dbh,    'dbh',         NC_DOUBLE, 'cm',        'diameter at breast height')
      call def_c(io%v_c_height, 'height',      NC_DOUBLE, 'm',         'height')
      call def_c(io%v_c_ba,     'basal_area',  NC_DOUBLE, 'cm2/plant', 'basal area per plant')
      call def_c(io%v_c_agb,    'agb',         NC_DOUBLE, 'kgC/plant', 'aboveground biomass per plant')
      call def_c(io%v_c_la,     'leaf_area',   NC_DOUBLE, 'm2/plant',  'leaf area per plant')
      call def_c(io%v_c_gavg,   'growth_avg',  NC_DOUBLE, 'cm/yr',     'simple moving-average growth (mortality predictor)')
      call def_c(io%v_c_owner,  'owner_patch', NC_INT,    '--',        'owning patch index (1-based)')
      call def_c(io%v_c_gid,    'global_cohort_id', NC_INT, '--',                                &
                 'persistent cohort id (stable across records until fused/culled)')

      !----- per-patch state (time, patch). -----------------------------------------------!
      call def_p(io%v_p_area,   'patch_area',     NC_DOUBLE, '--', 'patch area fraction')
      call def_p(io%v_p_age,    'patch_age',      NC_DOUBLE, 'yr', 'time since last disturbance')
      call def_p(io%v_p_dist,   'dist_type',      NC_INT,    '--', 'disturbance type (1=primary,2=treefall)')
      call def_p(io%v_p_coff,   'cohort_offset',  NC_INT,    '--', 'first cohort index of patch (1-based)')
      call def_p(io%v_p_ccount, 'cohort_count',   NC_INT,    '--', 'number of cohorts in patch')
      call def_p(io%v_p_gid,    'global_patch_id', NC_INT,   '--',                               &
                 'persistent patch id (stable across records until fused/culled)')

      !----- per-site totals (time). ------------------------------------------------------!
      call def_t(io%v_s_nplant, 'total_nplant',     NC_DOUBLE, 'plant/m2', 'site total plant number')
      call def_t(io%v_s_ba,     'total_basal_area', NC_DOUBLE, 'm2/m2',    'site total basal area')
      call def_t(io%v_s_agb,    'total_agb',        NC_DOUBLE, 'kgC/m2',   'site total aboveground biomass')
      call def_t(io%v_s_lai,    'total_lai',        NC_DOUBLE, 'm2/m2',    'site total leaf area index')
      call def_t(io%v_s_dmean,  'mean_dbh',         NC_DOUBLE, 'cm',       'basal-area-weighted mean DBH')

      call nc_check(nc_put_att_text_f(io%ncid, NC_GLOBAL, 'title',                         &
                    int(len_trim(TITLE), c_size_t), TITLE), 'global title')
      call put_global_text('start_time', time_to_string(cfg%start_time))
      call put_global_text('end_time',   time_to_string(cfg%end_time))
      call nc_check(nc_enddef(io%ncid), 'enddef')
      write(*,'(2a)') ' output: ', trim(path)

   contains
      !----- 1-D record variable (time). -----------------------------------------------!
      subroutine def_t(vid, name, xtype, units, lname)
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, units, lname
         integer(c_int),   intent(in)  :: xtype
         call nc_check(nc_def_var_f(io%ncid, name, xtype, 1_c_int, [io%d_time], vid), 'def '//name)
         call put_var_attrs(vid, units, lname)
      end subroutine def_t

      !----- 2-D per-cohort variable (time, cohort), chunked + deflated. ----------------!
      subroutine def_c(vid, name, xtype, units, lname)
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, units, lname
         integer(c_int),   intent(in)  :: xtype
         call nc_check(nc_def_var_f(io%ncid, name, xtype, 2_c_int,                          &
                       [io%d_time, io%d_cohort], vid), 'def '//name)
         call nc_check(nc_def_var_chunking(io%ncid, vid, NC_CHUNKED,                            &
                       [1_c_size_t, int(io%cohort_max, c_size_t)]), 'chunk '//name)
         call nc_check(nc_def_var_deflate(io%ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//name)
         call put_var_attrs(vid, units, lname)
      end subroutine def_c

      !----- 2-D per-patch variable (time, patch), chunked + deflated. ------------------!
      subroutine def_p(vid, name, xtype, units, lname)
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, units, lname
         integer(c_int),   intent(in)  :: xtype
         call nc_check(nc_def_var_f(io%ncid, name, xtype, 2_c_int,                          &
                       [io%d_time, io%d_patch], vid), 'def '//name)
         call nc_check(nc_def_var_chunking(io%ncid, vid, NC_CHUNKED,                            &
                       [1_c_size_t, int(io%patch_max, c_size_t)]), 'chunk '//name)
         call nc_check(nc_def_var_deflate(io%ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//name)
         call put_var_attrs(vid, units, lname)
      end subroutine def_p

      subroutine put_var_attrs(vid, units, lname)
         integer(c_int),   intent(in) :: vid
         character(len=*), intent(in) :: units, lname
         call nc_check(nc_put_att_text_f(io%ncid, vid, 'units',                             &
                       int(len_trim(units), c_size_t), units), 'units')
         call nc_check(nc_put_att_text_f(io%ncid, vid, 'long_name',                         &
                       int(len_trim(lname), c_size_t), lname), 'long_name')
      end subroutine put_var_attrs

      subroutine put_global_text(name, text)
         character(len=*), intent(in) :: name, text
         call nc_check(nc_put_att_text_f(io%ncid, NC_GLOBAL, name,                          &
                       int(len_trim(text), c_size_t), text), 'global '//name)
      end subroutine put_global_text
   end subroutine io_create

   !---------------------------------------------------------------------------------------!
   ! Append one time record with the current live state.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine io_write_snapshot(io, site, now)
      type(meds_io_t),   intent(inout) :: io
      type(site_t),      intent(in)    :: site
      type(meds_time_t), intent(in)    :: now
      integer(ik)       :: ncoh, npat
      integer(c_size_t) :: t0, sc(2), cc(2), sp(2), cp(2), i1(1)

      ncoh = site%cohort%n
      npat = site%patch%n
      if (ncoh > io%cohort_max) error stop 'meds_io: cohort count exceeds io_cohort_max'
      if (npat > io%patch_max)  error stop 'meds_io: patch count exceeds io_patch_max'

      t0 = int(io%nrec, c_size_t)                 ! 0-based record index
      i1 = [t0]
      sc = [t0, 0_c_size_t] ; cc = [1_c_size_t, int(ncoh, c_size_t)]   ! cohort slab
      sp = [t0, 0_c_size_t] ; cp = [1_c_size_t, int(npat, c_size_t)]   ! patch slab

      !----- calendar date + counts. ------------------------------------------------------!
      call nc_check(nc_put_var1_double(io%ncid, io%v_time, i1, time_to_decimal_year(now)), 'put time')
      call nc_check(nc_put_vara_int(io%ncid, io%v_year,  i1, [1_c_size_t], [int(now%year,  c_int)]), 'put year')
      call nc_check(nc_put_vara_int(io%ncid, io%v_month, i1, [1_c_size_t], [int(now%month, c_int)]), 'put month')
      call nc_check(nc_put_vara_int(io%ncid, io%v_day,   i1, [1_c_size_t], [int(now%day,   c_int)]), 'put day')
      call nc_check(nc_put_vara_int(io%ncid, io%v_ncohort, i1, [1_c_size_t], [int(ncoh, c_int)]), 'put n_cohort')
      call nc_check(nc_put_vara_int(io%ncid, io%v_npatch,  i1, [1_c_size_t], [int(npat, c_int)]), 'put n_patch')

      !----- per-cohort slabs (only the valid prefix 1:ncoh). -----------------------------!
      if (ncoh > 0_ik) then
         associate (c => site%cohort)
            call nc_check(nc_put_vara_int   (io%ncid, io%v_c_pft,    sc, cc, c%pft(1:ncoh)),        'put pft')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_nplant, sc, cc, c%nplant(1:ncoh)),     'put nplant')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_dbh,    sc, cc, c%dbh(1:ncoh)),        'put dbh')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_height, sc, cc, c%height(1:ncoh)),     'put height')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_ba,     sc, cc, c%basal_area(1:ncoh)), 'put basal_area')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_agb,    sc, cc, c%agb(1:ncoh)),        'put agb')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_la,     sc, cc, c%leaf_area(1:ncoh)),  'put leaf_area')
            call nc_check(nc_put_vara_double(io%ncid, io%v_c_gavg,   sc, cc, c%growth_avg(1:ncoh)), 'put growth_avg')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_c_owner,  sc, cc, c%owner_patch(1:ncoh)),'put owner_patch')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_c_gid,    sc, cc, c%global_id(1:ncoh)),  'put global_cohort_id')
         end associate
      end if

      !----- per-patch slabs. -------------------------------------------------------------!
      if (npat > 0_ik) then
         associate (p => site%patch)
            call nc_check(nc_put_vara_double(io%ncid, io%v_p_area,  sp, cp, p%area(1:npat)),           'put patch_area')
            call nc_check(nc_put_vara_double(io%ncid, io%v_p_age,   sp, cp, p%age(1:npat)),            'put patch_age')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_p_dist,  sp, cp, p%dist_type(1:npat)),      'put dist_type')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_p_coff,  sp, cp, p%cohort_offset(1:npat)),  'put cohort_offset')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_p_ccount,sp, cp, p%cohort_count(1:npat)),   'put cohort_count')
            call nc_check(nc_put_vara_int   (io%ncid, io%v_p_gid,   sp, cp, p%global_id(1:npat)),      'put global_patch_id')
         end associate
      end if

      !----- per-site totals. -------------------------------------------------------------!
      call nc_check(nc_put_var1_double(io%ncid, io%v_s_nplant, i1, total_nplant(site)),     'put total_nplant')
      call nc_check(nc_put_var1_double(io%ncid, io%v_s_ba,     i1, total_basal_area(site)), 'put total_basal_area')
      call nc_check(nc_put_var1_double(io%ncid, io%v_s_agb,    i1, total_agb(site)),        'put total_agb')
      call nc_check(nc_put_var1_double(io%ncid, io%v_s_lai,    i1, total_lai(site)),        'put total_lai')
      call nc_check(nc_put_var1_double(io%ncid, io%v_s_dmean,  i1, mean_dbh(site)),         'put mean_dbh')

      io%nrec = io%nrec + 1_ik
   end subroutine io_write_snapshot

   subroutine io_close(io)
      type(meds_io_t), intent(inout) :: io
      if (io%ncid >= 0_c_int) call nc_check(nc_close(io%ncid), 'nc_close')
      io%ncid = -1_c_int
   end subroutine io_close

   !=======================================================================================!
   !  STATE (restart) files                                                                 !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Write the full INSTANTANEOUS prognostic state to a self-contained, timestamped file     !
   ! <dir>/<prefix>-S-YYYYMMDDHHMMSS.nc -- everything needed to restart, and NO diagnostics.  !
   ! Cached geometry (height/basal_area/agb/leaf_area) is omitted: it is re-derived from dbh  !
   ! on restart. Scalars travel in a small meta_int/meta_real vector.                         !
   !---------------------------------------------------------------------------------------!
   subroutine io_write_state(site, cfg, dir, prefix, now)
      type(site_t),        intent(in) :: site
      type(meds_config_t), intent(in) :: cfg
      character(len=*),    intent(in) :: dir, prefix
      type(meds_time_t),   intent(in) :: now
      integer(c_int) :: ncid, d_cohort, d_patch, d_pft, d_mi, d_mr, d_soill, d_snowl
      integer(c_int) :: vmi, vmr, vc_pft, vc_np, vc_dbh, vc_own, vc_gid, vc_gavg
      integer(c_int) :: vc_sla, vc_vc, vc_rd, vc_ll        ! plastic leaf traits
      integer(c_int) :: vc_lwm, vc_wwm, vc_lt, vc_wt       ! P6: per-cohort hydraulics/temperature state
      integer(c_int) :: vc_dmax, vc_dmax_acc                     ! #95: per-cohort predawn-psi stomatal feedback
      integer(c_int) :: vp_area, vp_age, vp_dist, vp_gid, vp_rec
      integer(c_int) :: vp_sc1, vp_sc2, vp_sc3, vp_sc4, vp_sc5, vp_sc6, vp_sc7, vp_lig1, vp_lig2
      !----- FAST reservoirs (P5 restart-completeness fix, MEDS_ED2_RK45_DESIGN.md): persisted so a  !
      !      restart resumes the true evolved CAS/soil/snow state instead of a generic re-seed --      !
      !      the mismatch between a freshly-reset reservoir and an already-mature canopy's demand       !
      !      is exactly the kind of stiff transient RK45's explicit stepper has no L-stable defense     !
      !      against (unlike split's implicit CAS box or ARK's Newton surface solve). ------------------!
      integer(c_int) :: vf_centh, vf_cshv, vf_cco2, vf_ctemp
      integer(c_int) :: vf_theta, vf_wsurf, vf_wsenth
      integer(c_int) :: vf_se, vf_stemp, vf_sfliq
      integer(c_int) :: vf_swe, vf_sneng, vf_sdep, vf_stmp, vf_sfl, vf_snl
      integer(ik)    :: ncoh, npat, npft, ip, meta_i(11)
      real(wp)       :: meta_r(2)
      character(len=512) :: fname

      ncoh = site%cohort%n ; npat = site%patch%n ; npft = site%n_pft
      !----- The timestamp IS the simulated calendar date -> <prefix>-S-YYYYMMDDHHMMSS.nc. --!
      fname = trim(dir)//'/'//trim(prefix)//'-S-'//time_to_stamp(now)//'.nc'

      call nc_check(nc_create_f(trim(fname), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'state nc_create')
      call nc_check(nc_def_dim_f(ncid, 'cohort', int(max(ncoh,1_ik), c_size_t), d_cohort), 'state dim cohort')
      call nc_check(nc_def_dim_f(ncid, 'patch',  int(max(npat,1_ik), c_size_t), d_patch),  'state dim patch')
      call nc_check(nc_def_dim_f(ncid, 'pft',    int(npft, c_size_t),           d_pft),    'state dim pft')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_int',  11_c_size_t, d_mi), 'state dim mi')
      call nc_check(nc_def_dim_f(ncid, 'nmeta_real',  2_c_size_t, d_mr), 'state dim mr')
      call nc_check(nc_def_dim_f(ncid, 'soil_layer', int(n_soil_layer_max, c_size_t), d_soill), 'state dim soill')
      call nc_check(nc_def_dim_f(ncid, 'snow_layer', int(n_snow_layer_max, c_size_t), d_snowl), 'state dim snowl')

      call dv(vmi,    'meta_int',         NC_INT,    [d_mi],                                   &
              '[n_cohort,n_patch,n_pft,next_cohort_id,next_patch_id,year,month,day,hour,minute,second]')
      call dv(vmr,    'meta_real',        NC_DOUBLE, [d_mr],     '[site_area, decimal_year]')
      call dv(vc_pft, 'pft',              NC_INT,    [d_cohort], 'plant functional type index')
      call dv(vc_np,  'nplant',           NC_DOUBLE, [d_cohort], 'plant number density [plant/m2]')
      call dv(vc_dbh, 'dbh',              NC_DOUBLE, [d_cohort], 'diameter at breast height [cm]')
      call dv(vc_gavg,'growth_avg',       NC_DOUBLE, [d_cohort], 'simple moving-average growth [cm/yr] (mortality predictor)')
      call dv(vc_own, 'owner_patch',      NC_INT,    [d_cohort], 'owning patch index (1-based)')
      call dv(vc_gid, 'global_cohort_id', NC_INT,    [d_cohort], 'persistent cohort id')
      call dv(vc_sla, 'sla',              NC_DOUBLE, [d_cohort], 'specific leaf area [m2/kgC] (plastic)')
      call dv(vc_vc,  'vcmax25',          NC_DOUBLE, [d_cohort], 'max carboxylation @25C [umol/m2/s] (plastic)')
      call dv(vc_rd,  'rd25',             NC_DOUBLE, [d_cohort], 'leaf dark respiration @25C [umol/m2/s] (plastic)')
      call dv(vc_ll,  'llspan',           NC_DOUBLE, [d_cohort], 'leaf lifespan [yr] (plastic)')
      !----- P6 (MEDS_ED2_RK45_DESIGN.md): per-cohort plant-hydraulics internal water + tissue temps.  !
      !      Persisting these makes a restart EXACT for the fast loop -- without them the restart        !
      !      re-seeds leaf_water_mass to near-saturated (lazy-init) and leaf/wood_temp to LEAF_TEMP_INIT, !
      !      whose psi/temperature discontinuity makes the plant-hydraulics sub-stepper grind for a few   !
      !      days on a healthy high-LAI restart. OPTIONAL on read (gv_dbl_opt), so an older-format state   !
      !      file still restarts (re-seeding as before). Leaf-surface (interception) water stays          !
      !      unpersisted on purpose: 0 (bone-dry) is a real, discontinuity-free initial condition. -------!
      call dv(vc_lwm, 'leaf_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal leaf water [kg/plant]')
      call dv(vc_wwm, 'wood_water_mass',  NC_DOUBLE, [d_cohort], 'per-cohort internal wood water [kg/plant]')
      call dv(vc_lt,  'leaf_temp',        NC_DOUBLE, [d_cohort], 'per-cohort leaf temperature [K]')
      call dv(vc_wt,  'wood_temp',        NC_DOUBLE, [d_cohort], 'per-cohort wood temperature [K]')
      !----- #95 per-cohort stomatal water-stress feedback. `dmax_psi_leaf` is YESTERDAY's completed !
      !      daily-max leaf water potential and drives TODAY's beta_stomata; `dmax_psi_leaf_accum` is the !
      !      running   !
      !      accumulator the fast loop fills, rolled over at day end. Both must persist:                  !
      !                                                                                          !
      !        * without predawn, a restart falls back to the UNSET sentinel and re-seeds from surface-   !
      !          layer soil psi -- which reads as unstressed, so a restart mid-drought silently re-opens  !
      !          stomata that the running model had closed.                                                !
      !        * without the accumulator, a restart at any time other than midnight discards the part of  !
      !          today's maximum already seen, so the value handed to tomorrow is biased LOW (too         !
      !          negative => too much stress) by however much of the day was replayed.                     !
      !                                                                                          !
      !      OPTIONAL on read like the P6 block above, so an older state file still restarts on the       !
      !      sentinels. NOT derivable from leaf_water_mass: predawn is a time-MAXIMUM over the previous   !
      !      day, not a function of the instantaneous state. ---------------------------------------------!
      call dv(vc_dmax, 'dmax_psi_leaf',   NC_DOUBLE, [d_cohort], 'yesterday daily-max leaf water potential [MPa] (drives beta_stomata)')
      call dv(vc_dmax_acc, 'dmax_psi_leaf_accum', NC_DOUBLE, [d_cohort], 'running daily-max leaf water potential accumulator [MPa]')
      call dv(vp_area,'patch_area',       NC_DOUBLE, [d_patch],  'patch area fraction')
      call dv(vp_age, 'patch_age',        NC_DOUBLE, [d_patch],  'time since last disturbance [yr]')
      call dv(vp_dist,'dist_type',        NC_INT,    [d_patch],  'disturbance type (1=primary,2=treefall)')
      call dv(vp_gid, 'global_patch_id',  NC_INT,    [d_patch],  'persistent patch id')
      call dv(vp_rec, 'recruit_pool',     NC_DOUBLE, [d_patch, d_pft], 'carry-forward recruit pool [plant/m2]')
      !----- FAST reservoirs (P5, MEDS_ED2_RK45_DESIGN.md): the true evolved CAS/soil/snow state, not    !
      !      just its cached geometry -- a restart that reset these to a generic seed (the ONLY prior      !
      !      behavior) handed an already-mature canopy a discontinuous jump on day 1, which RK45's fully   !
      !      explicit stepper (no L-stable defense, unlike split's implicit CAS box or ARK's Newton         !
      !      surface solve) could only resolve by grinding through many tiny adaptive substeps. Diagnosed    !
      !      twins (soil_temp/fliq, snow_temp/fliq, cas_temp) are saved directly alongside their             !
      !      prognostic partner rather than recomputed on read, so this reader needs no outside PFT/soil     !
      !      config (matching every other read here). Always written (not opt-in): reading them back is     !
      !      itself optional (gv_*_opt), so an OLDER state file lacking them still falls back to the          !
      !      previous behavior (meds_main re-seeds via init_fast_reservoirs) with no format-version bump. ---!
      call dv(vf_centh, 'cas_can_enthalpy', NC_DOUBLE, [d_patch], 'CAS specific enthalpy [J/kg] (prognostic)')
      call dv(vf_cshv,  'cas_can_shv',      NC_DOUBLE, [d_patch], 'CAS specific humidity [kg/kg] (prognostic)')
      call dv(vf_cco2,  'cas_can_co2',      NC_DOUBLE, [d_patch], 'CAS CO2 mixing ratio [umol/mol] (prognostic)')
      call dv(vf_ctemp, 'cas_can_temp',     NC_DOUBLE, [d_patch], 'CAS temperature [K] (diagnosed twin)')
      call dv(vf_theta, 'soil_theta',       NC_DOUBLE, [d_patch, d_soill], 'soil volumetric moisture [m3/m3] (prognostic)')
      call dv(vf_wsurf, 'soil_w_surface',   NC_DOUBLE, [d_patch], 'ponded surface water [kg/m2]')
      call dv(vf_wsenth,'soil_w_surface_enth', NC_DOUBLE, [d_patch], 'ponded-water internal energy [J/m2]')
      call dv(vf_se,    'soil_energy',      NC_DOUBLE, [d_patch, d_soill], 'soil volumetric internal energy [J/m3] (prognostic)')
      call dv(vf_stemp, 'soil_temp',        NC_DOUBLE, [d_patch, d_soill], 'soil temperature [K] (diagnosed twin)')
      call dv(vf_sfliq, 'soil_fliq',        NC_DOUBLE, [d_patch, d_soill], 'soil liquid-water fraction [-] (diagnosed twin)')
      call dv(vf_swe,   'snow_swe',         NC_DOUBLE, [d_patch, d_snowl], 'snow water-equivalent mass [kg/m2] (prognostic)')
      call dv(vf_sneng, 'snow_energy',      NC_DOUBLE, [d_patch, d_snowl], 'snow extensive internal energy [J/m2] (prognostic)')
      call dv(vf_sdep,  'snow_depth',       NC_DOUBLE, [d_patch, d_snowl], 'snow geometric depth [m] (prognostic)')
      call dv(vf_stmp,  'snow_temp',        NC_DOUBLE, [d_patch, d_snowl], 'snow temperature [K] (diagnosed twin)')
      call dv(vf_sfl,   'snow_fliq',        NC_DOUBLE, [d_patch, d_snowl], 'snow liquid fraction [-] (diagnosed twin)')
      call dv(vf_snl,   'snow_nlayer',      NC_INT,    [d_patch], 'active snow layer count (0 = no snow)')
      !----- Slow soil-carbon pools (opt-in, [soil_carbon].soil_carbon_on; MEDS_SLOW_DYNAMICS_DESIGN.md !
      !      Part II B0). N-cycle fields are skipped (n_cycle_on defaults false; C-only MVP). -------------!
      call dv(vp_sc1, 'soilc_fast_grnd',   NC_DOUBLE, [d_patch], 'fast/metabolic litter carbon, above-ground [kgC/m2]')
      call dv(vp_sc2, 'soilc_fast_soil',   NC_DOUBLE, [d_patch], 'fast/metabolic litter carbon, below-ground [kgC/m2]')
      call dv(vp_sc3, 'soilc_struct_grnd', NC_DOUBLE, [d_patch], 'structural litter + CWD carbon, above-ground [kgC/m2]')
      call dv(vp_sc4, 'soilc_struct_soil', NC_DOUBLE, [d_patch], 'structural litter + CWD carbon, below-ground [kgC/m2]')
      call dv(vp_sc5, 'soilc_microbial',   NC_DOUBLE, [d_patch], 'microbial SOM carbon [kgC/m2]')
      call dv(vp_sc6, 'soilc_slow',        NC_DOUBLE, [d_patch], 'slow/humified SOM carbon [kgC/m2]')
      call dv(vp_sc7, 'soilc_passive',     NC_DOUBLE, [d_patch], 'passive SOM carbon [kgC/m2]')
      call dv(vp_lig1,'soilc_lignin_grnd', NC_DOUBLE, [d_patch], 'structural litter lignin, above-ground [kgC/m2]')
      call dv(vp_lig2,'soilc_lignin_soil', NC_DOUBLE, [d_patch], 'structural litter lignin, below-ground [kgC/m2]')
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, 'title',                            &
                    int(len_trim(STATE_TITLE), c_size_t), STATE_TITLE), 'state title')
      call nc_check(nc_enddef(ncid), 'state enddef')

      meta_i = [ncoh, npat, npft, site%next_cohort_id, site%next_patch_id,                     &
                now%year, now%month, now%day, now%hour, now%minute, now%second]
      meta_r = [site%site_area, time_to_decimal_year(now)]
      call nc_check(nc_put_vara_int   (ncid, vmi, [0_c_size_t], [11_c_size_t], meta_i), 'put meta_int')
      call nc_check(nc_put_vara_double(ncid, vmr, [0_c_size_t], [2_c_size_t], meta_r), 'put meta_real')
      if (ncoh > 0_ik) then
         associate (c => site%cohort)
            call nc_check(nc_put_vara_int   (ncid, vc_pft, [0_c_size_t], [int(ncoh,c_size_t)], c%pft(1:ncoh)),        'put pft')
            call nc_check(nc_put_vara_double(ncid, vc_np,  [0_c_size_t], [int(ncoh,c_size_t)], c%nplant(1:ncoh)),     'put nplant')
            call nc_check(nc_put_vara_double(ncid, vc_dbh, [0_c_size_t], [int(ncoh,c_size_t)], c%dbh(1:ncoh)),        'put dbh')
            call nc_check(nc_put_vara_double(ncid, vc_gavg, [0_c_size_t], [int(ncoh,c_size_t)],          &
                          c%growth_avg(1:ncoh)), 'put growth_avg')
            call nc_check(nc_put_vara_int   (ncid, vc_own, [0_c_size_t], [int(ncoh,c_size_t)], c%owner_patch(1:ncoh)),'put owner')
            call nc_check(nc_put_vara_int   (ncid, vc_gid, [0_c_size_t], [int(ncoh,c_size_t)], c%global_id(1:ncoh)),  'put cgid')
            call nc_check(nc_put_vara_double(ncid, vc_sla, [0_c_size_t], [int(ncoh,c_size_t)], c%sla(1:ncoh)),       'put sla')
            call nc_check(nc_put_vara_double(ncid, vc_vc,  [0_c_size_t], [int(ncoh,c_size_t)], c%vcmax25(1:ncoh)),   'put vcmax25')
            call nc_check(nc_put_vara_double(ncid, vc_rd,  [0_c_size_t], [int(ncoh,c_size_t)], c%rd25(1:ncoh)),      'put rd25')
            call nc_check(nc_put_vara_double(ncid, vc_ll,  [0_c_size_t], [int(ncoh,c_size_t)], c%llspan(1:ncoh)),    'put llspan')
            call nc_check(nc_put_vara_double(ncid, vc_lwm, [0_c_size_t], [int(ncoh,c_size_t)], c%leaf_water_mass(1:ncoh)), 'put leaf_water_mass')
            call nc_check(nc_put_vara_double(ncid, vc_wwm, [0_c_size_t], [int(ncoh,c_size_t)], c%wood_water_mass(1:ncoh)), 'put wood_water_mass')
            call nc_check(nc_put_vara_double(ncid, vc_lt,  [0_c_size_t], [int(ncoh,c_size_t)], c%leaf_temp(1:ncoh)), 'put leaf_temp')
            call nc_check(nc_put_vara_double(ncid, vc_wt,  [0_c_size_t], [int(ncoh,c_size_t)], c%wood_temp(1:ncoh)), 'put wood_temp')
            call nc_check(nc_put_vara_double(ncid, vc_dmax, [0_c_size_t], [int(ncoh,c_size_t)], c%dmax_psi_leaf(1:ncoh)), 'put dmax_psi_leaf')
            call nc_check(nc_put_vara_double(ncid, vc_dmax_acc, [0_c_size_t], [int(ncoh,c_size_t)], c%dmax_psi_leaf_accum(1:ncoh)), 'put dmax_psi_leaf_accum')
         end associate
      end if
      if (npat > 0_ik) then
         associate (p => site%patch)
            call nc_check(nc_put_vara_double(ncid, vp_area, [0_c_size_t], [int(npat,c_size_t)], p%area(1:npat)),      'put area')
            call nc_check(nc_put_vara_double(ncid, vp_age,  [0_c_size_t], [int(npat,c_size_t)], p%age(1:npat)),       'put age')
            call nc_check(nc_put_vara_int   (ncid, vp_dist, [0_c_size_t], [int(npat,c_size_t)], p%dist_type(1:npat)), 'put dist')
            call nc_check(nc_put_vara_int   (ncid, vp_gid,  [0_c_size_t], [int(npat,c_size_t)], p%global_id(1:npat)), 'put pgid')
            do ip = 1_ik, npat
               call nc_check(nc_put_vara_double(ncid, vp_rec, [int(ip-1_ik,c_size_t), 0_c_size_t],    &
                             [1_c_size_t, int(npft,c_size_t)], p%recruit_pool(1:npft, ip)), 'put recruit_pool')
            end do
            block                                 ! FAST reservoirs (P5): contiguous buffers, same
               real(wp) :: fc(npat, 4_ik)          ! discipline as the soil-carbon block below (never a
               real(wp) :: fw(npat, 4_ik)          ! bare derived-type component-array-section straight
               real(wp) :: fsoil(npat, n_soil_layer_max, 3_ik)     ! into a C-bound call)
               real(wp) :: fsnow(npat, n_snow_layer_max, 5_ik)
               integer(ik) :: fnl(npat)
               do ip = 1_ik, npat
                  fc(ip,1) = p%cas(ip)%can_enthalpy ; fc(ip,2) = p%cas(ip)%can_shv
                  fc(ip,3) = p%cas(ip)%can_co2      ; fc(ip,4) = p%cas(ip)%can_temp
                  fw(ip,1) = p%soil_w(ip)%w_surface ; fw(ip,2) = p%soil_w(ip)%w_surface_enth
                  fsoil(ip,:,1) = p%soil_w(ip)%theta(1:n_soil_layer_max)
                  fsoil(ip,:,2) = p%soil_e(ip)%soil_energy(1:n_soil_layer_max)
                  fsoil(ip,:,3) = p%soil_e(ip)%soil_fliq(1:n_soil_layer_max)
                  !----- soil_temp shares fsoil's slab shape but comes from a separate diagnosed array; !
                  !      written as its own call below (a 4th soil field would not fit the (:,3) shape). !
                  fsnow(ip,:,1) = p%snow(ip)%swe(1:n_snow_layer_max)
                  fsnow(ip,:,2) = p%snow(ip)%snow_energy(1:n_snow_layer_max)
                  fsnow(ip,:,3) = p%snow(ip)%snow_depth(1:n_snow_layer_max)
                  fsnow(ip,:,4) = p%snow(ip)%snow_temp(1:n_snow_layer_max)
                  fsnow(ip,:,5) = p%snow(ip)%snow_fliq(1:n_snow_layer_max)
                  fnl(ip) = p%snow(ip)%nlayer
               end do
               call nc_check(nc_put_vara_double(ncid, vf_centh, [0_c_size_t], [int(npat,c_size_t)], fc(:,1)), 'put cas_can_enthalpy')
               call nc_check(nc_put_vara_double(ncid, vf_cshv,  [0_c_size_t], [int(npat,c_size_t)], fc(:,2)), 'put cas_can_shv')
               call nc_check(nc_put_vara_double(ncid, vf_cco2,  [0_c_size_t], [int(npat,c_size_t)], fc(:,3)), 'put cas_can_co2')
               call nc_check(nc_put_vara_double(ncid, vf_ctemp, [0_c_size_t], [int(npat,c_size_t)], fc(:,4)), 'put cas_can_temp')
               call nc_check(nc_put_vara_double(ncid, vf_wsurf, [0_c_size_t], [int(npat,c_size_t)], fw(:,1)), 'put soil_w_surface')
               call nc_check(nc_put_vara_double(ncid, vf_wsenth,[0_c_size_t], [int(npat,c_size_t)], fw(:,4)), 'put soil_w_surface_enth')
               call nc_check(nc_put_vara_int   (ncid, vf_snl,   [0_c_size_t], [int(npat,c_size_t)], fnl),     'put snow_nlayer')
               do ip = 1_ik, npat
                  call nc_check(nc_put_vara_double(ncid, vf_theta, [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_soil_layer_max,c_size_t)], fsoil(ip,:,1)), 'put soil_theta')
                  call nc_check(nc_put_vara_double(ncid, vf_se,    [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_soil_layer_max,c_size_t)], fsoil(ip,:,2)), 'put soil_energy')
                  call nc_check(nc_put_vara_double(ncid, vf_sfliq, [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_soil_layer_max,c_size_t)], fsoil(ip,:,3)), 'put soil_fliq')
                  call nc_check(nc_put_vara_double(ncid, vf_stemp, [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_soil_layer_max,c_size_t)],                              &
                                p%soil_e(ip)%soil_temp(1:n_soil_layer_max)), 'put soil_temp')
                  call nc_check(nc_put_vara_double(ncid, vf_swe,   [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_snow_layer_max,c_size_t)], fsnow(ip,:,1)), 'put snow_swe')
                  call nc_check(nc_put_vara_double(ncid, vf_sneng, [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_snow_layer_max,c_size_t)], fsnow(ip,:,2)), 'put snow_energy')
                  call nc_check(nc_put_vara_double(ncid, vf_sdep,  [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_snow_layer_max,c_size_t)], fsnow(ip,:,3)), 'put snow_depth')
                  call nc_check(nc_put_vara_double(ncid, vf_stmp,  [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_snow_layer_max,c_size_t)], fsnow(ip,:,4)), 'put snow_temp')
                  call nc_check(nc_put_vara_double(ncid, vf_sfl,   [int(ip-1_ik,c_size_t), 0_c_size_t],   &
                                [1_c_size_t, int(n_snow_layer_max,c_size_t)], fsnow(ip,:,5)), 'put snow_fliq')
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
               call nc_check(nc_put_vara_double(ncid, vp_sc1,  [0_c_size_t], [int(npat,c_size_t)], sc(:,1)), 'put soilc_fast_grnd')
               call nc_check(nc_put_vara_double(ncid, vp_sc2,  [0_c_size_t], [int(npat,c_size_t)], sc(:,2)), 'put soilc_fast_soil')
               call nc_check(nc_put_vara_double(ncid, vp_sc3,  [0_c_size_t], [int(npat,c_size_t)], sc(:,3)), 'put soilc_struct_grnd')
               call nc_check(nc_put_vara_double(ncid, vp_sc4,  [0_c_size_t], [int(npat,c_size_t)], sc(:,4)), 'put soilc_struct_soil')
               call nc_check(nc_put_vara_double(ncid, vp_sc5,  [0_c_size_t], [int(npat,c_size_t)], sc(:,5)), 'put soilc_microbial')
               call nc_check(nc_put_vara_double(ncid, vp_sc6,  [0_c_size_t], [int(npat,c_size_t)], sc(:,6)), 'put soilc_slow')
               call nc_check(nc_put_vara_double(ncid, vp_sc7,  [0_c_size_t], [int(npat,c_size_t)], sc(:,7)), 'put soilc_passive')
               call nc_check(nc_put_vara_double(ncid, vp_lig1, [0_c_size_t], [int(npat,c_size_t)], sc(:,8)), 'put soilc_lignin_grnd')
               call nc_check(nc_put_vara_double(ncid, vp_lig2, [0_c_size_t], [int(npat,c_size_t)], sc(:,9)), 'put soilc_lignin_soil')
            end block
         end associate
      end if
      call nc_check(nc_close(ncid), 'state nc_close')
      write(*,'(2a)') ' state : ', trim(fname)

   contains
      subroutine dv(vid, name, xtype, dimids, lname)        ! define a state variable + long_name
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, lname
         integer(c_int),   intent(in)  :: xtype, dimids(:)
         call nc_check(nc_def_var_f(ncid, name, xtype, int(size(dimids), c_int), dimids, vid), 'state def '//name)
         call nc_check(nc_put_att_text_f(ncid, vid, 'long_name',                              &
                       int(len_trim(lname), c_size_t), lname), 'state long_name '//name)
      end subroutine dv
   end subroutine io_write_state

   !---------------------------------------------------------------------------------------!
   ! Reconstruct a site from a state file written by io_write_state. The cached geometry is   !
   ! re-derived from dbh (gather_pft_params + set_cohort_size); the CSR map and within-patch   !
   ! sort order are rebuilt. found=.false. (no error stop) if the file cannot be opened, so    !
   ! the caller can fall back. Errors out only on a genuine inconsistency (PFT-count mismatch).!
   !---------------------------------------------------------------------------------------!
   subroutine io_read_state(site, cfg, path, restart_time, found, fast_found)
      type(site_t),        intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      type(meds_time_t),   intent(out) :: restart_time
      logical,             intent(out) :: found
      !----- FAST reservoirs (P5): .true. only when THIS file's cas/soil/snow state was actually    !
      !      read back (always written together, so one variable's presence implies all of them);     !
      !      absent/.false. means an older-format file -- the caller (meds_main) must then fall back    !
      !      to init_fast_reservoirs's generic seed, exactly the pre-P5 behavior. --------------------!
      logical, optional,   intent(out) :: fast_found
      integer(c_int) :: ncid, vid, vrec, st
      integer(ik)    :: ncoh, npat, npft, ip, i, nwin, meta_i(11)
      real(wp)       :: meta_r(2)
      logical        :: fast_ok

      found = .false. ; restart_time = meds_time_t()
      if (present(fast_found)) fast_found = .false.
      st = nc_open_f(trim(path), NC_NOWRITE, ncid)
      if (st /= NC_NOERR) return

      call nc_check(nc_inq_varid_f(ncid, 'meta_int',  vid), 'inq meta_int')
      call nc_check(nc_get_vara_int(ncid, vid, [0_c_size_t], [11_c_size_t], meta_i), 'get meta_int')
      call nc_check(nc_inq_varid_f(ncid, 'meta_real', vid), 'inq meta_real')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t], [2_c_size_t], meta_r), 'get meta_real')
      ncoh = meta_i(1) ; npat = meta_i(2) ; npft = meta_i(3)
      if (npft /= cfg%pft%n) error stop 'io_read_state: PFT count in state file /= config PFT count'

      call site_alloc(site, cfg%pft%n, max(ncoh, 1_ik), max(npat, 1_ik), growth_window_steps(cfg))
      site%cohort%n = ncoh ; site%patch%n = npat
      site%next_cohort_id = meta_i(4) ; site%next_patch_id = meta_i(5)
      site%site_area = meta_r(1)
      restart_time = meds_time_t(year=meta_i(6), month=meta_i(7), day=meta_i(8),               &
                                 hour=meta_i(9), minute=meta_i(10), second=meta_i(11))

      if (ncoh > 0_ik) then
         associate (c => site%cohort)
            call gv_int (ncid, 'pft',              ncoh, c%pft(1:ncoh))
            call gv_dbl (ncid, 'nplant',           ncoh, c%nplant(1:ncoh))
            call gv_dbl (ncid, 'dbh',              ncoh, c%dbh(1:ncoh))
            call gv_dbl (ncid, 'growth_avg',       ncoh, c%growth_avg(1:ncoh))  ! ring buffer reseeded below
            call gv_int (ncid, 'owner_patch',      ncoh, c%owner_patch(1:ncoh))
            call gv_int (ncid, 'global_cohort_id', ncoh, c%global_id(1:ncoh))
            !----- Plastic leaf traits: default to PFT top-of-canopy, then overwrite from the state  !
            !       file where present (states written before this feature restart at top-of-canopy).!
            do i = 1_ik, ncoh
               c%sla(i)     = cfg%pft%sla(c%pft(i))
               c%vcmax25(i) = cfg%pft%vcmax25(c%pft(i))
               c%rd25(i)    = cfg%pft%rd25(c%pft(i))
               c%llspan(i)  = cfg%pft%leaf_lifespan_toc(c%pft(i))
            end do
            call gv_dbl_opt(ncid, 'sla',     ncoh, c%sla(1:ncoh))
            call gv_dbl_opt(ncid, 'vcmax25', ncoh, c%vcmax25(1:ncoh))
            call gv_dbl_opt(ncid, 'rd25',    ncoh, c%rd25(1:ncoh))
            call gv_dbl_opt(ncid, 'llspan',  ncoh, c%llspan(1:ncoh))
            !----- P6 per-cohort hydraulics/temperature (OPTIONAL): site_alloc has already set the      !
            !      lazy-init sentinels (leaf/wood_water_mass = 0 -> re-seed near-saturated on first       !
            !      touch; leaf/wood_temp = LEAF_TEMP_INIT), so an older state file lacking these keeps    !
            !      exactly the pre-P6 re-seed behaviour; a current file overwrites with the true state.   !
            call gv_dbl_opt(ncid, 'leaf_water_mass', ncoh, c%leaf_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_water_mass', ncoh, c%wood_water_mass(1:ncoh))
            call gv_dbl_opt(ncid, 'leaf_temp',       ncoh, c%leaf_temp(1:ncoh))
            call gv_dbl_opt(ncid, 'wood_temp',       ncoh, c%wood_temp(1:ncoh))
            !----- #95 stomatal feedback (OPTIONAL): alloc_column_cohort has already set               !
            !      DMAX_PSI_LEAF_UNSET / DMAX_PSI_LEAF_ACCUM_RESET, so an older state file keeps the seed-from- !
            !      soil behaviour; a current file restores the true feedback state. -------------------!
            call gv_dbl_opt(ncid, 'dmax_psi_leaf',   ncoh, c%dmax_psi_leaf(1:ncoh))
            call gv_dbl_opt(ncid, 'dmax_psi_leaf_accum', ncoh, c%dmax_psi_leaf_accum(1:ncoh))
         end associate
         call gather_pft_params(site%cohort, cfg%pft)        ! p_dbh_critical / p_wood_density
         do i = 1_ik, ncoh
            call set_cohort_size(site%cohort, i)             ! height/basal_area/agb/leaf_area from dbh
         end do
         !----- Seed the moving-average ring buffer as if it were full of the saved growth_avg !
         !       (the per-step history itself is not stored); it is overwritten within a window.!
         nwin = growth_window_steps(cfg)
         do i = 1_ik, ncoh
            if (site%cohort%growth_avg(i) >= 0.0_wp) then
               site%cohort%growth_count(i)   = nwin
               site%cohort%growth_accum(i)   = site%cohort%growth_avg(i) * real(nwin, wp)
               site%cohort%growth_hist(:, i) = site%cohort%growth_avg(i)
            end if
         end do
      end if
      if (npat > 0_ik) then
         associate (p => site%patch)
            call gv_dbl (ncid, 'patch_area',      npat, p%area(1:npat))
            call gv_dbl (ncid, 'patch_age',       npat, p%age(1:npat))
            call gv_int (ncid, 'dist_type',       npat, p%dist_type(1:npat))
            call gv_int (ncid, 'global_patch_id', npat, p%global_id(1:npat))
            call nc_check(nc_inq_varid_f(ncid, 'recruit_pool', vrec), 'inq recruit_pool')
            do ip = 1_ik, npat
               call nc_check(nc_get_vara_double(ncid, vrec, [int(ip-1_ik,c_size_t), 0_c_size_t],  &
                             [1_c_size_t, int(npft,c_size_t)], p%recruit_pool(1:npft, ip)), 'get recruit_pool')
            end do
            !----- FAST reservoirs (P5): OPTIONAL as a WHOLE (existence of the first variable stands   !
            !      for all of them, always written together) -- an older-format file leaves p%cas/       !
            !      soil_e/soil_w/snow at whatever site_alloc/patch_alloc's own defaults are, and the      !
            !      caller (meds_main) re-seeds via init_fast_reservoirs exactly as it did pre-P5. ---------!
            fast_ok = nc_inq_varid_f(ncid, 'cas_can_enthalpy', vid) == NC_NOERR
            if (fast_ok) then
               block
                  real(wp) :: fc(npat), fw1(npat), fw2(npat), fw3(npat), fw4(npat)
                  logical  :: pond_enth_ok
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
                     p%soil_w(ip)%w_surface = fw1(ip)
                  end do
                  !----- Pond ENTHALPY (issue #78 item 4): OPTIONAL, because state files written before
                  !      the pond had a thermal state do not carry it. Read the value here if present;
                  !      the RECONSTRUCTION for an older file has to wait until soil_temp is in (below).
                  pond_enth_ok = nc_inq_varid_f(ncid, 'soil_w_surface_enth', vid) == NC_NOERR
                  fw4 = 0.0_wp
                  if (pond_enth_ok) call gv_dbl(ncid, 'soil_w_surface_enth', npat, fw4)
                  call gv_int(ncid, 'snow_nlayer', npat, fnl)
                  do ip = 1_ik, npat ; p%snow(ip)%nlayer = fnl(ip) ; end do
                  do ip = 1_ik, npat
                     call gv_dbl2_row(ncid, 'soil_theta',  ip, n_soil_layer_max, p%soil_w(ip)%theta(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_energy', ip, n_soil_layer_max, p%soil_e(ip)%soil_energy(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_temp',   ip, n_soil_layer_max, p%soil_e(ip)%soil_temp(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'soil_fliq',   ip, n_soil_layer_max, p%soil_e(ip)%soil_fliq(1:n_soil_layer_max))
                     call gv_dbl2_row(ncid, 'snow_swe',    ip, n_snow_layer_max, p%snow(ip)%swe(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_energy', ip, n_snow_layer_max, p%snow(ip)%snow_energy(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_depth',  ip, n_snow_layer_max, p%snow(ip)%snow_depth(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_temp',   ip, n_snow_layer_max, p%snow(ip)%snow_temp(1:n_snow_layer_max))
                     call gv_dbl2_row(ncid, 'snow_fliq',   ip, n_snow_layer_max, p%snow(ip)%snow_fliq(1:n_snow_layer_max))
                  end do
                  !----- Pond ENTHALPY, committed now that soil_temp is in. When the file carries it, use  !
                  !      it verbatim. When it does not, RECONSTRUCT from the pond mass at the top soil      !
                  !      layer's temperature rather than defaulting to 0: an older file records a pond      !
                  !      whose enthalpy is UNKNOWN, not empty, and seeding 0 J/m2 would restart the run     !
                  !      ~1 MJ short per kg of ponded water -- which the newly-armed budgets would          !
                  !      immediately (and correctly) report as a leak. The top layer's temperature is the   !
                  !      one the pond exchanges with, so it is the defensible equilibrium guess. A dry      !
                  !      pond (the common case) reconstructs to exactly 0 either way. --------------------!
                  do ip = 1_ik, npat
                     if (pond_enth_ok) then
                        p%soil_w(ip)%w_surface_enth = fw4(ip)
                     else
                        p%soil_w(ip)%w_surface_enth = p%soil_w(ip)%w_surface                            &
                              * internal_energy_liquid(p%soil_e(ip)%soil_temp(1))
                     end if
                  end do
               end block
               if (present(fast_found)) fast_found = .true.
            end if
            !----- Slow soil-carbon pools (opt-in, [soil_carbon].soil_carbon_on): OPTIONAL variables --!
            !      (gv_dbl_opt), so a state file written before this feature restarts the pools at the  !
            !      just-allocated zero (bare-ground philosophy), matching the plastic-leaf-trait pattern.!
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
      call nc_check(nc_close(ncid), 'state nc_close (read)')

      call rebuild_csr(site)
      call sort_cohorts(site)
      found = .true.

   contains
      subroutine gv_int(nc, name, n, out)        ! read [1:n] of an int variable by name
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: n
         integer(ik),      intent(out) :: out(:)
         integer(c_int) :: v
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_int(nc, v, [0_c_size_t], [int(n,c_size_t)], out), 'get '//name)
      end subroutine gv_int
      subroutine gv_dbl(nc, name, n, out)        ! read [1:n] of a double variable by name
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: n
         real(wp),         intent(out) :: out(:)
         integer(c_int) :: v
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_double(nc, v, [0_c_size_t], [int(n,c_size_t)], out), 'get '//name)
      end subroutine gv_dbl
      subroutine gv_dbl_opt(nc, name, n, out)    ! read if the var exists; leave `out` unchanged otherwise
         integer(c_int),   intent(in)    :: nc
         character(len=*), intent(in)    :: name
         integer(ik),      intent(in)    :: n
         real(wp),         intent(inout) :: out(:)
         integer(c_int) :: v
         if (nc_inq_varid_f(nc, name, v) /= NC_NOERR) return
         call nc_check(nc_get_vara_double(nc, v, [0_c_size_t], [int(n,c_size_t)], out), 'get '//name)
      end subroutine gv_dbl_opt
      subroutine gv_dbl2_row(nc, name, row, nlayer, out)   ! read row [row,1:nlayer] of a [patch,layer] var
         integer(c_int),   intent(in)  :: nc
         character(len=*), intent(in)  :: name
         integer(ik),      intent(in)  :: row, nlayer
         real(wp),         intent(out) :: out(:)
         integer(c_int) :: v
         call nc_check(nc_inq_varid_f(nc, name, v), 'inq '//name)
         call nc_check(nc_get_vara_double(nc, v, [int(row-1_ik,c_size_t), 0_c_size_t],              &
                       [1_c_size_t, int(nlayer,c_size_t)], out), 'get '//name)
      end subroutine gv_dbl2_row
   end subroutine io_read_state

end module meds_io
