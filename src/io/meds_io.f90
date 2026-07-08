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
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_alloc, gather_pft_params, set_cohort_size, rebuild_csr
   use meds_demography_fusefiss,   only : sort_cohorts
   use meds_demography_diagnostics, only : total_nplant, total_basal_area, total_agb,         &
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
      integer(c_int) :: ncid, d_cohort, d_patch, d_pft, d_mi, d_mr
      integer(c_int) :: vmi, vmr, vc_pft, vc_np, vc_dbh, vc_own, vc_gid, vc_gavg
      integer(c_int) :: vp_area, vp_age, vp_dist, vp_gid, vp_rec
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

      call dv(vmi,    'meta_int',         NC_INT,    [d_mi],                                   &
              '[n_cohort,n_patch,n_pft,next_cohort_id,next_patch_id,year,month,day,hour,minute,second]')
      call dv(vmr,    'meta_real',        NC_DOUBLE, [d_mr],     '[site_area, decimal_year]')
      call dv(vc_pft, 'pft',              NC_INT,    [d_cohort], 'plant functional type index')
      call dv(vc_np,  'nplant',           NC_DOUBLE, [d_cohort], 'plant number density [plant/m2]')
      call dv(vc_dbh, 'dbh',              NC_DOUBLE, [d_cohort], 'diameter at breast height [cm]')
      call dv(vc_gavg,'growth_avg',       NC_DOUBLE, [d_cohort], 'simple moving-average growth [cm/yr] (mortality predictor)')
      call dv(vc_own, 'owner_patch',      NC_INT,    [d_cohort], 'owning patch index (1-based)')
      call dv(vc_gid, 'global_cohort_id', NC_INT,    [d_cohort], 'persistent cohort id')
      call dv(vp_area,'patch_area',       NC_DOUBLE, [d_patch],  'patch area fraction')
      call dv(vp_age, 'patch_age',        NC_DOUBLE, [d_patch],  'time since last disturbance [yr]')
      call dv(vp_dist,'dist_type',        NC_INT,    [d_patch],  'disturbance type (1=primary,2=treefall)')
      call dv(vp_gid, 'global_patch_id',  NC_INT,    [d_patch],  'persistent patch id')
      call dv(vp_rec, 'recruit_pool',     NC_DOUBLE, [d_patch, d_pft], 'carry-forward recruit pool [plant/m2]')
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
            call nc_check(nc_put_vara_double(ncid, vc_gavg,[0_c_size_t], [int(ncoh,c_size_t)], c%growth_avg(1:ncoh)),'put growth_avg')
            call nc_check(nc_put_vara_int   (ncid, vc_own, [0_c_size_t], [int(ncoh,c_size_t)], c%owner_patch(1:ncoh)),'put owner')
            call nc_check(nc_put_vara_int   (ncid, vc_gid, [0_c_size_t], [int(ncoh,c_size_t)], c%global_id(1:ncoh)),  'put cgid')
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
   subroutine io_read_state(site, cfg, path, restart_time, found)
      type(site_t),        intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      type(meds_time_t),   intent(out) :: restart_time
      logical,             intent(out) :: found
      integer(c_int) :: ncid, vid, vrec, st
      integer(ik)    :: ncoh, npat, npft, ip, i, nwin, meta_i(11)
      real(wp)       :: meta_r(2)

      found = .false. ; restart_time = meds_time_t()
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
   end subroutine io_read_state

end module meds_io
