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
   use meds_config,   only : meds_config_t
   use meds_netcdf_c
   use meds_demography_interface,   only : site_t
   use meds_demography_diagnostics, only : total_nplant, total_basal_area, total_agb,         &
                                           total_lai, mean_dbh
   implicit none
   private

   public :: meds_io_t, io_create, io_write_snapshot, io_close

   character(len=*), parameter :: TITLE = 'MEDS demographic state output'

   type :: meds_io_t
      integer(c_int) :: ncid = -1_c_int
      integer(ik)    :: nrec = 0_ik
      integer(ik)    :: cohort_max = 0_ik, patch_max = 0_ik
      integer(c_int) :: d_time, d_cohort, d_patch
      integer(c_int) :: v_time, v_ncohort, v_npatch
      integer(c_int) :: v_c_pft, v_c_nplant, v_c_dbh, v_c_height, v_c_ba, v_c_agb, v_c_la, v_c_owner
      integer(c_int) :: v_c_gid
      integer(c_int) :: v_p_area, v_p_age, v_p_dist, v_p_atemp, v_p_mtemp, v_p_coff, v_p_ccount
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

      call nc_check(nc_create(cstr(path), ior(NC_NETCDF4, NC_CLOBBER), io%ncid), 'nc_create')
      call nc_check(nc_def_dim(io%ncid, cstr('time'),   NC_UNLIMITED,                  io%d_time),   'dim time')
      call nc_check(nc_def_dim(io%ncid, cstr('cohort'), int(io%cohort_max, c_size_t),  io%d_cohort), 'dim cohort')
      call nc_check(nc_def_dim(io%ncid, cstr('patch'),  int(io%patch_max,  c_size_t),  io%d_patch),  'dim patch')

      !----- time coordinate + per-record counts. -----------------------------------------!
      call def_t(io%v_time,    'time',     NC_DOUBLE, 'yr', 'time since start of run')
      call def_t(io%v_ncohort, 'n_cohort', NC_INT,    '--', 'number of cohorts in use this record')
      call def_t(io%v_npatch,  'n_patch',  NC_INT,    '--', 'number of patches in use this record')

      !----- per-cohort state (time, cohort). ---------------------------------------------!
      call def_c(io%v_c_pft,    'pft',         NC_INT,    '--',        'plant functional type index')
      call def_c(io%v_c_nplant, 'nplant',      NC_DOUBLE, 'plant/m2',  'plant number density')
      call def_c(io%v_c_dbh,    'dbh',         NC_DOUBLE, 'cm',        'diameter at breast height')
      call def_c(io%v_c_height, 'height',      NC_DOUBLE, 'm',         'height')
      call def_c(io%v_c_ba,     'basal_area',  NC_DOUBLE, 'cm2/plant', 'basal area per plant')
      call def_c(io%v_c_agb,    'agb',         NC_DOUBLE, 'kgC/plant', 'aboveground biomass per plant')
      call def_c(io%v_c_la,     'leaf_area',   NC_DOUBLE, 'm2/plant',  'leaf area per plant')
      call def_c(io%v_c_owner,  'owner_patch', NC_INT,    '--',        'owning patch index (1-based)')
      call def_c(io%v_c_gid,    'global_cohort_id', NC_INT, '--',                                &
                 'persistent cohort id (stable across records until fused/culled)')

      !----- per-patch state (time, patch). -----------------------------------------------!
      call def_p(io%v_p_area,   'patch_area',     NC_DOUBLE, '--', 'patch area fraction')
      call def_p(io%v_p_age,    'patch_age',      NC_DOUBLE, 'yr', 'time since last disturbance')
      call def_p(io%v_p_dist,   'dist_type',      NC_INT,    '--', 'disturbance type (1=primary,2=treefall)')
      call def_p(io%v_p_atemp,  'avg_daily_temp', NC_DOUBLE, 'K',  'mean daily temperature')
      call def_p(io%v_p_mtemp,  'min_month_temp', NC_DOUBLE, 'K',  'minimum monthly temperature')
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

      call nc_check(nc_put_att_text(io%ncid, NC_GLOBAL, cstr('title'),                         &
                    int(len_trim(TITLE), c_size_t), cstr(TITLE)), 'global title')
      call nc_check(nc_enddef(io%ncid), 'enddef')

   contains
      !----- 1-D record variable (time). -----------------------------------------------!
      subroutine def_t(vid, name, xtype, units, lname)
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, units, lname
         integer(c_int),   intent(in)  :: xtype
         call nc_check(nc_def_var(io%ncid, cstr(name), xtype, 1_c_int, [io%d_time], vid), 'def '//name)
         call put_var_attrs(vid, units, lname)
      end subroutine def_t

      !----- 2-D per-cohort variable (time, cohort), chunked + deflated. ----------------!
      subroutine def_c(vid, name, xtype, units, lname)
         integer(c_int),   intent(out) :: vid
         character(len=*), intent(in)  :: name, units, lname
         integer(c_int),   intent(in)  :: xtype
         call nc_check(nc_def_var(io%ncid, cstr(name), xtype, 2_c_int,                          &
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
         call nc_check(nc_def_var(io%ncid, cstr(name), xtype, 2_c_int,                          &
                       [io%d_time, io%d_patch], vid), 'def '//name)
         call nc_check(nc_def_var_chunking(io%ncid, vid, NC_CHUNKED,                            &
                       [1_c_size_t, int(io%patch_max, c_size_t)]), 'chunk '//name)
         call nc_check(nc_def_var_deflate(io%ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//name)
         call put_var_attrs(vid, units, lname)
      end subroutine def_p

      subroutine put_var_attrs(vid, units, lname)
         integer(c_int),   intent(in) :: vid
         character(len=*), intent(in) :: units, lname
         call nc_check(nc_put_att_text(io%ncid, vid, cstr('units'),                             &
                       int(len_trim(units), c_size_t), cstr(units)), 'units')
         call nc_check(nc_put_att_text(io%ncid, vid, cstr('long_name'),                         &
                       int(len_trim(lname), c_size_t), cstr(lname)), 'long_name')
      end subroutine put_var_attrs
   end subroutine io_create

   !---------------------------------------------------------------------------------------!
   ! Append one time record with the current live state.                                    !
   !---------------------------------------------------------------------------------------!
   subroutine io_write_snapshot(io, site, time_years)
      type(meds_io_t), intent(inout) :: io
      type(site_t),    intent(in)    :: site
      real(wp),        intent(in)    :: time_years
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

      !----- time + counts. ---------------------------------------------------------------!
      call nc_check(nc_put_var1_double(io%ncid, io%v_time, i1, time_years), 'put time')
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
            call nc_check(nc_put_vara_double(io%ncid, io%v_p_atemp, sp, cp, p%avg_daily_temp(1:npat)), 'put avg_daily_temp')
            call nc_check(nc_put_vara_double(io%ncid, io%v_p_mtemp, sp, cp, p%min_month_temp(1:npat)), 'put min_month_temp')
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

end module meds_io
