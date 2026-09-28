! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_stream -- the netCDF serializer for the diagnostic-aggregation streams: create a  !
! per-tier, per-time-chunk file (dims + registry-driven variable defs + CF metadata), append one  !
! averaged record, and roll to a new file at the chunk boundary (§5). Reuses meds_netcdf_c.       !
!                                                                                          !
! netCDF wall: this is the ONLY diagnostic module that touches C; referenced from the manager      !
! (main-only), never from the stepper. Generalizes meds_io's io_create/io_write_snapshot from a     !
! hard-coded field list to a loop over the registry's per-tier live-variable index.                 !
!==========================================================================================!
module meds_output_stream
   use iso_c_binding, only : c_int, c_size_t, c_double
   use meds_kinds,    only : wp, ik
   use meds_time,     only : meds_time_t, time_to_decimal_year
   use meds_column_params, only : n_soil_layer_max
   use meds_netcdf_c
   use meds_output_config, only : FC_DAY, FC_MONTH, FC_YEAR, FC_RUN, SYNC_FLUSH, freq_letter,     &
                                  freq_tier_index
   use meds_output_types,  only : output_registry_t, stream_file_t, pending_record_t, var_desc_t, slab_col, &
                                  output_files_t, output_buffers_t,                                  &
                                  DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,            &
                                  DIM_SIZE, DIM_SOIL_PATCH, diag_params_t,                        &
                                  XTYPE_DOUBLE, XTYPE_INT, AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX,    &
                                  AGG_LAST, AGG_VARIANCE, AGG_TMEAN, AGG_FLUXSUM,                  &
                                  MISSING_VALUE, MISSING_INT
   implicit none
   private

   public :: stream_write_record, stream_close_file, region_write_record

   character(len=*), parameter :: TITLE = 'MEDS diagnostic aggregation output'

contains

   !----- Write one closed-period record, opening / rolling the per-tier file as needed. ------!
   subroutine stream_write_record(stream, reg, dg, pr, dir, prefix, file_chunk, cohort_max,      &
                                  patch_max, sync_every, forcing_qair)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(diag_params_t),     intent(in)    :: dg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, cohort_max, patch_max, sync_every
      character(len=*),        intent(in)    :: forcing_qair   !< provenance attribute ('' -> none)
      integer(ik) :: tier, bucket, fc
      tier = freq_tier_index(pr%freq)
      !----- Cohort/patch counts are invariant only WITHIN a month (§4.4); to trim the cohort/patch !
      !      dimension to the live count the file must not span more than a month. Cap the effective  !
      !      file_chunk to FC_MONTH for any tier that carries a cohort/patch variable (site-only tiers !
      !      keep their configured chunk, e.g. the annual run-file).  ------------------------------!
      fc = file_chunk
      if (tier_has_cohort_or_patch(reg, tier)) then
         fc = min(fc, FC_MONTH)
         !----- FAST (tier 1) cohort/patch: force one file PER DAY. Cohort/patch counts are invariant   !
         !      within a day (fusion/fission is monthly/annual), so the fixed-slot ≤1-day window keeps    !
         !      the count-grew guard from firing at sub-daily cadence (§4.4; no global_id keying).  ------!
         if (tier == 1_ik) fc = FC_DAY
      end if
      bucket = bucket_key(pr%t_open, fc)
      if (stream%ncid < 0_ik .or. bucket /= stream%chunk_bucket) then
         call stream_close_file(stream)
         call stream_open_file(stream, reg, dg, pr, dir, prefix, fc, bucket, cohort_max,          &
                               patch_max, tier, forcing_qair)
      end if
      call write_one_record(stream, reg, pr, tier)
      !----- Skip the per-record nc_sync for the FAST tier: ~n_fast_per_slow records/day would else    !
      !      each fsync. The chunk-boundary nc_close still flushes the file. Coarse tiers honor sync.    !
      if (sync_every == SYNC_FLUSH .and. tier /= 1_ik) call nc_check(nc_sync(int(stream%ncid, c_int)), 'nc_sync')
   end subroutine stream_write_record

   subroutine stream_close_file(stream)
      type(stream_file_t), intent(inout) :: stream
      if (stream%ncid >= 0_ik) call nc_check(nc_close(int(stream%ncid, c_int)), 'nc_close stream')
      stream%ncid = -1_ik ; stream%nrec = 0_ik ; stream%chunk_bucket = -1_ik
      stream%cohort_dim = 0_ik ; stream%patch_dim = 0_ik
   end subroutine stream_close_file

   !----- .true. if the tier defines a cohort- or patch-dimensioned variable (drives the ≤1-month  !
   !      file-chunk cap, since those axes are invariant only within a month, §4.4). ---------------!
   pure logical function tier_has_cohort_or_patch(reg, tier) result(yes)
      type(output_registry_t), intent(in) :: reg
      integer(ik),             intent(in) :: tier
      integer(ik) :: j, k
      yes = .false.
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         if (reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) then
            yes = .true. ; return
         end if
      end do
   end function tier_has_cohort_or_patch

   !----- Integer bucket key identifying the time-chunk a period belongs to (§5.2). ----------!
   pure integer(ik) function bucket_key(t, file_chunk) result(b)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: file_chunk
      select case (file_chunk)
      case (FC_DAY)   ; b = t%year*10000_ik + t%month*100_ik + t%day
      case (FC_MONTH) ; b = t%year*100_ik + t%month
      case (FC_YEAR)  ; b = t%year
      case default    ; b = 0_ik                    ! FC_RUN: single bucket
      end select
   end function bucket_key

   !----- Filename stamp for the time-chunk (empty for FC_RUN -> stampless -Y.nc, §5.1). ------!
   pure function chunk_stamp(t, file_chunk) result(s)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: file_chunk
      character(len=16) :: s
      s = ''
      select case (file_chunk)
      case (FC_DAY)   ; write(s,'(i4.4,i2.2,i2.2)') t%year, t%month, t%day
      case (FC_MONTH) ; write(s,'(i4.4,i2.2)')      t%year, t%month
      case (FC_YEAR)  ; write(s,'(i4.4)')           t%year
      end select
   end function chunk_stamp

   pure function cell_methods_of(agg) result(cm)
      integer(ik), intent(in) :: agg
      character(len=16) :: cm
      select case (agg)
      case (AGG_MEAN, AGG_TMEAN) ; cm = 'time: mean'
      case (AGG_VARIANCE)        ; cm = 'time: variance'
      case (AGG_SUM, AGG_FLUXSUM)            ; cm = 'time: sum'
      case (AGG_MIN)                         ; cm = 'time: minimum'
      case (AGG_MAX)                         ; cm = 'time: maximum'
      case default                           ; cm = 'time: point'   ! AGG_LAST
      end select
   end function cell_methods_of

   !=======================================================================================!
   !  Create the file + define dims / registry-driven variables / CF metadata (§5.3).        !
   !=======================================================================================!
   subroutine stream_open_file(stream, reg, dg, pr, dir, prefix, file_chunk, bucket, cohort_max,  &
                               patch_max, tier, forcing_qair)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(diag_params_t),     intent(in)    :: dg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, bucket, cohort_max, patch_max, tier
      character(len=*),        intent(in)    :: forcing_qair
      character(len=512) :: path
      character(len=16)  :: stamp
      character(len=1)   :: letter
      integer(c_int)     :: ncid, dt, dc, dp, ds, dpf, dsz, vid, dims1(1)
      integer(ik)        :: j, k, cohort_dim, patch_dim
      logical            :: hasc, hasp, hass, haspf, hassz

      !----- which trailing dims does this tier need? A 2-D (soil layer x patch) variable needs   !
      !      BOTH the patch and soil dims, which is why DIM_SOIL_PATCH sets two flags.  ----------!
      hasc = .false. ; hasp = .false. ; hass = .false. ; haspf = .false. ; hassz = .false.
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         select case (reg%var(k)%dim)
         case (DIM_COHORT)     ; hasc = .true.
         case (DIM_PATCH)      ; hasp = .true.
         case (DIM_SOIL)       ; hass = .true.
         case (DIM_PFT)        ; haspf = .true.
         case (DIM_SIZE)       ; hassz = .true.
         case (DIM_SOIL_PATCH) ; hasp = .true. ; hass = .true.
         end select
      end do

      !----- Trim the cohort/patch axes to the live count of the record that opens this file. The     !
      !      count is invariant across a ≤1-month file (the caller capped file_chunk to FC_MONTH for   !
      !      cohort/patch tiers), so this dim fits every record in the file; write_one_record asserts   !
      !      it. max(.,1) avoids a zero-length dim on an empty (bare-ground) window.  ------------------!
      if (pr%n_cohort > cohort_max) error stop 'meds_output_stream: n_cohort exceeds cohort_max buffer'
      if (pr%n_patch  > patch_max)  error stop 'meds_output_stream: n_patch exceeds patch_max buffer'
      cohort_dim = max(pr%n_cohort, 1_ik)
      patch_dim  = max(pr%n_patch,  1_ik)

      letter = freq_letter(pr%freq)
      stamp  = chunk_stamp(pr%t_open, file_chunk)
      if (len_trim(stamp) > 0) then
         path = trim(dir)//'/'//trim(prefix)//'-'//letter//'-'//trim(stamp)//'.nc'
      else
         path = trim(dir)//'/'//trim(prefix)//'-'//letter//'.nc'
      end if

      call nc_check(nc_create_f(trim(path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'stream nc_create')
      call nc_check(nc_def_dim_f(ncid, 'time', NC_UNLIMITED, dt), 'dim time')
      dc = -1_c_int ; dp = -1_c_int ; ds = -1_c_int ; dpf = -1_c_int ; dsz = -1_c_int
      if (hasc) call nc_check(nc_def_dim_f(ncid, 'cohort', int(cohort_dim, c_size_t), dc), 'dim cohort')
      if (hasp) call nc_check(nc_def_dim_f(ncid, 'patch',  int(patch_dim,  c_size_t), dp), 'dim patch')
      if (hass) call nc_check(nc_def_dim_f(ncid, 'soil',   int(n_soil_layer_max, c_size_t), ds), 'dim soil')
      if (haspf) call nc_check(nc_def_dim_f(ncid, 'pft',   int(max(dg%n_pft,1_ik), c_size_t), dpf), 'dim pft')
      if (hassz) call nc_check(nc_def_dim_f(ncid, 'dbh_class',                                       &
                               int(max(dg%n_dbh_class,1_ik), c_size_t), dsz), 'dim dbh_class')

      !----- time coordinate + calendar companions (period-start stamp). ---!
      stream%v_time  = int(def_scalar_var(ncid, dt, 'time',  NC_DOUBLE, 'year', 'decimal calendar year (period start)'), ik)
      stream%v_year  = int(def_scalar_var(ncid, dt, 'year',  NC_INT,    '1', 'calendar year (period start)'), ik)
      stream%v_month = int(def_scalar_var(ncid, dt, 'month', NC_INT,    '1', 'calendar month (period start)'), ik)
      stream%v_day   = int(def_scalar_var(ncid, dt, 'day',   NC_INT,    '1', 'calendar day (period start)'), ik)
      !----- FAST (tier 1): human-readable sub-daily companions so a reader can group-by-hour without    !
      !      decoding the decimal `time`. Period-start stamp, matching the calendar companions above.  --!
      stream%v_hour = -1_ik ; stream%v_minute = -1_ik ; stream%v_second = -1_ik
      if (tier == 1_ik) then
         stream%v_hour   = int(def_scalar_var(ncid, dt, 'hour',   NC_INT, '1', 'clock hour (period start)'), ik)
         stream%v_minute = int(def_scalar_var(ncid, dt, 'minute', NC_INT, '1', 'clock minute (period start)'), ik)
         stream%v_second = int(def_scalar_var(ncid, dt, 'second', NC_INT, '1', 'clock second (period start)'), ik)
      end if
      if (hasc) stream%v_ncohort = int(def_scalar_var(ncid, dt, 'n_cohort', NC_INT, '1', 'live cohorts this record'), ik)
      if (hasp) stream%v_npatch  = int(def_scalar_var(ncid, dt, 'n_patch',  NC_INT, '1', 'live patches this record'), ik)

      !----- SELF-DESCRIBING AXIS COORDINATES. Both the pft and dbh_class axes have RUN-DEPENDENT   !
      !      lengths (the PFT count comes from the PFT table, the class edges from TOML), so a file  !
      !      that carried only the data would be un-interpretable next to a file from another run.   !
      !      Writing the coordinates makes each file stand on its own.                                !
      stream%v_pft = -1_ik ; stream%v_dbh_lower = -1_ik ; stream%v_dbh_upper = -1_ik
      stream%v_soil_z = -1_ik
      if (haspf) then
         dims1 = [dpf]
         call nc_check(nc_def_var_f(ncid, 'pft', NC_INT, 1_c_int, dims1, vid), 'def pft coord')
         call put_var_text(ncid, vid, 'long_name', 'plant functional type index')
         stream%v_pft = int(vid, ik)
      end if
      if (hass) then
         dims1 = [ds]
         call nc_check(nc_def_var_f(ncid, 'soil_z', NC_DOUBLE, 1_c_int, dims1, vid), 'def soil_z')
         call put_var_text(ncid, vid, 'units', 'm')
         call put_var_text(ncid, vid, 'long_name', 'soil layer node depth (negative downward)')
         stream%v_soil_z = int(vid, ik)
      end if
      if (hassz) then
         dims1 = [dsz]
         call nc_check(nc_def_var_f(ncid, 'dbh_lower', NC_DOUBLE, 1_c_int, dims1, vid), 'def dbh_lower')
         call put_var_text(ncid, vid, 'units', 'cm')
         call put_var_text(ncid, vid, 'long_name', 'DBH class lower edge (inclusive)')
         stream%v_dbh_lower = int(vid, ik)
         call nc_check(nc_def_var_f(ncid, 'dbh_upper', NC_DOUBLE, 1_c_int, dims1, vid), 'def dbh_upper')
         call put_var_text(ncid, vid, 'units', 'cm')
         call put_var_text(ncid, vid, 'long_name',                                                  &
              'DBH class upper edge (exclusive, except the last class which is closed)')
         stream%v_dbh_upper = int(vid, ik)
      end if

      !----- registry-driven variable definitions (dims per DIM_*, chunk+deflate + CF attrs). ---!
      stream%vid = -1_ik
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         vid = def_registry_var(ncid, dt, dc, dp, ds, dpf, dsz, reg%var(k), cohort_dim, patch_dim, dg)
         stream%vid(k) = int(vid, ik)
      end do

      call put_global(ncid, 'title', TITLE)
      call put_global(ncid, 'Conventions', 'CF-1.10')
      if (len_trim(forcing_qair) > 0) call put_global(ncid, 'forcing_qair', trim(forcing_qair))
      call nc_check(nc_enddef(ncid), 'stream enddef')

      !----- Write the axis coordinates once, right after enddef (they do not vary by record). --!
      if (hass)  call write_soil_coord(ncid, stream%v_soil_z, dg)
      if (haspf) call write_pft_coord(ncid, stream%v_pft, dg%n_pft)
      if (hassz) call write_size_coord(ncid, stream%v_dbh_lower, stream%v_dbh_upper, dg)

      stream%ncid = int(ncid, ik) ; stream%nrec = 0_ik ; stream%chunk_bucket = bucket
      stream%has_cohort = hasc ; stream%has_patch = hasp ; stream%has_soil = hass
      stream%has_pft = haspf ; stream%has_size = hassz
      stream%cohort_dim = cohort_dim ; stream%patch_dim = patch_dim
      stream%d_time = int(dt, ik) ; stream%d_cohort = int(dc, ik)
      stream%d_patch = int(dp, ik) ; stream%d_soil = int(ds, ik)
      stream%d_pft = int(dpf, ik) ; stream%d_size = int(dsz, ik)
      write(*,'(2a)') ' output: ', trim(path)
   end subroutine stream_open_file

   !----- The soil coordinate: layer node depths [m], negative downward. Inactive layers past    !
   !      n_soil are left at 0 and are also where the data slabs carry _FillValue, so a reader can   !
   !      mask on either.  ----------------------------------------------------------------------!
   subroutine write_soil_coord(ncid, vid_ik, dg)
      integer(c_int),      intent(in) :: ncid
      integer(ik),         intent(in) :: vid_ik
      type(diag_params_t), intent(in) :: dg
      real(c_double)    :: z(n_soil_layer_max)
      integer(c_size_t) :: st(1), cn(1)
      integer(ik)       :: k, nz
      !----- Only the ACTIVE layers (#246). Writing the full ceiling put a run of 0.0 m node      !
      !      depths on the end of the coordinate, which reads as ten more layers all at the        !
      !      surface; the unwritten tail is left for netCDF to fill, like every other padded axis. !
      nz = min(max(dg%n_soil, 1_ik), n_soil_layer_max)
      do k = 1_ik, nz ; z(k) = real(dg%soil_z(k), c_double) ; end do
      st = [0_c_size_t] ; cn = [int(nz, c_size_t)]
      call nc_check(nc_put_vara_double(ncid, int(vid_ik, c_int), st, cn, z(1:nz)), 'put soil_z')
   end subroutine write_soil_coord

   !----- The pft coordinate: the 1-based PFT indices this run carries. ----------------------!
   subroutine write_pft_coord(ncid, vid_ik, n_pft)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: vid_ik, n_pft
      integer(c_int)    :: iarr(max(n_pft,1_ik))
      integer(c_size_t) :: st(1), cn(1)
      integer(ik)       :: i
      do i = 1_ik, max(n_pft, 1_ik) ; iarr(i) = int(i, c_int) ; end do
      st = [0_c_size_t] ; cn = [int(max(n_pft,1_ik), c_size_t)]
      call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), st, cn, iarr), 'put pft coord')
   end subroutine write_pft_coord

   !----- The dbh_class coordinates: the lower and upper edge of each class, so a reader can    !
   !      label the axis without knowing the run's TOML.                                         !
   subroutine write_size_coord(ncid, vlo_ik, vhi_ik, dg)
      integer(c_int),      intent(in) :: ncid
      integer(ik),         intent(in) :: vlo_ik, vhi_ik
      type(diag_params_t), intent(in) :: dg
      real(c_double)    :: lo(max(dg%n_dbh_class,1_ik)), hi(max(dg%n_dbh_class,1_ik))
      integer(c_size_t) :: st(1), cn(1)
      integer(ik)       :: i, nc
      nc = max(dg%n_dbh_class, 1_ik)
      do i = 1_ik, nc
         lo(i) = real(dg%dbh_edges(i),        c_double)
         hi(i) = real(dg%dbh_edges(i + 1_ik), c_double)
      end do
      st = [0_c_size_t] ; cn = [int(nc, c_size_t)]
      call nc_check(nc_put_vara_double(ncid, int(vlo_ik, c_int), st, cn, lo), 'put dbh_lower')
      call nc_check(nc_put_vara_double(ncid, int(vhi_ik, c_int), st, cn, hi), 'put dbh_upper')
   end subroutine write_size_coord

   !----- One text attribute on an already-defined variable. --------------------------------!
   subroutine put_var_text(ncid, vid, name, text)
      integer(c_int),   intent(in) :: ncid, vid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, vid, name, int(len_trim(text), c_size_t), text), 'att '//name)
   end subroutine put_var_text

   !----- Define a 1-D record variable (time) + units/long_name; returns its varid. ----------!
   integer(c_int) function def_scalar_var(ncid, dt, name, xtype, units, lname) result(vid)
      integer(c_int),   intent(in) :: ncid, dt, xtype
      character(len=*), intent(in) :: name, units, lname
      integer(c_int) :: dims(1)
      dims = [dt]
      call nc_check(nc_def_var_f(ncid, name, xtype, 1_c_int, dims, vid), 'def '//name)
      call nc_check(nc_put_att_text_f(ncid, vid, 'units', int(len_trim(units), c_size_t), units), 'units '//name)
      call nc_check(nc_put_att_text_f(ncid, vid, 'long_name', int(len_trim(lname), c_size_t), lname), 'lname '//name)
   end function def_scalar_var

   !----- Define one registry variable at its dim, chunk+deflate slabs, CF attrs + _FillValue. -!
   integer(c_int) function def_registry_var(ncid, dt, dc, dp, ds, dpf, dsz, v, cohort_dim,       &
                                            patch_dim, dg) result(vid)
      integer(c_int),      intent(in) :: ncid, dt, dc, dp, ds, dpf, dsz
      type(var_desc_t),    intent(in) :: v
      integer(ik),         intent(in) :: cohort_dim, patch_dim  !< the file's TRIMMED (live-count) axis lengths
      type(diag_params_t), intent(in) :: dg
      integer(c_int)    :: xt, dims(2), dims1(1), dims3(3), axislen
      integer(c_size_t) :: chunk2(2), chunk3(3)
      character(len=16) :: cm
      xt = merge(NC_INT, NC_DOUBLE, v%xtype == XTYPE_INT)
      if (v%dim == DIM_SCALAR) then
         dims1 = [dt]                       ! named local (avoids the nvfortran/ifx arg-temp trap, issue #7)
         call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 1_c_int, dims1, vid), 'def '//trim(v%name))
      else if (v%dim == DIM_SOIL_PATCH) then
         !----- The 2-D axis: RANK-3 on disk, (time, patch, soil). The buffer side is a flat 1-D   !
         !      slab strided by n_soil_layer_max (see DIM_SOIL_PATCH in meds_output_types), so the  !
         !      layout here and the flattening there must agree: patch-major, layer-minor. netCDF   !
         !      is row-major in C order, and nc_put_vara takes count = [1, n_patch, n_soil], which  !
         !      walks layer fastest -- exactly the flattened order. No new C binding is needed:     !
         !      nc_put_vara_double's startp/countp are assumed-size.  ------------------------------!
         dims3  = [dt, dp, ds]
         chunk3 = [1_c_size_t, int(patch_dim, c_size_t), int(n_soil_layer_max, c_size_t)]
         call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 3_c_int, dims3, vid), 'def '//trim(v%name))
         call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk3), 'chunk '//trim(v%name))
         call nc_check(nc_def_var_deflate(ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//trim(v%name))
      else
         select case (v%dim)
         case (DIM_COHORT) ; dims = [dt, dc]  ; axislen = int(cohort_dim, c_int)
         case (DIM_PATCH)  ; dims = [dt, dp]  ; axislen = int(patch_dim,  c_int)
         case (DIM_SOIL)   ; dims = [dt, ds]  ; axislen = int(n_soil_layer_max, c_int)
         case (DIM_PFT)    ; dims = [dt, dpf] ; axislen = int(max(dg%n_pft, 1_ik), c_int)
         case (DIM_SIZE)   ; dims = [dt, dsz] ; axislen = int(max(dg%n_dbh_class, 1_ik), c_int)
         case default      ; dims = [dt, dc]  ; axislen = int(cohort_dim, c_int)
         end select
         chunk2 = [1_c_size_t, int(axislen, c_size_t)]
         call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 2_c_int, dims, vid), 'def '//trim(v%name))
         call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk2), 'chunk '//trim(v%name))
         call nc_check(nc_def_var_deflate(ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//trim(v%name))
      end if
      call nc_check(nc_put_att_text_f(ncid, vid, 'units', int(len_trim(v%units), c_size_t), trim(v%units)), 'units')
      call nc_check(nc_put_att_text_f(ncid, vid, 'long_name', int(len_trim(v%long_name), c_size_t), trim(v%long_name)), 'lname')
      cm = cell_methods_of(v%agg)
      call nc_check(nc_put_att_text_f(ncid, vid, 'cell_methods', int(len_trim(cm), c_size_t), trim(cm)), 'cell_methods')
      if (v%xtype == XTYPE_INT) then
         call nc_check(nc_put_att_int_f(ncid, vid, '_FillValue', NC_INT, int(MISSING_INT, c_int)), 'fill')
      else
         call nc_check(nc_put_att_double_f(ncid, vid, '_FillValue', NC_DOUBLE, real(MISSING_VALUE, c_double)), 'fill')
      end if
   end function def_registry_var

   subroutine put_global(ncid, name, text)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, name, int(len_trim(text), c_size_t), text), 'global '//name)
   end subroutine put_global

   !=======================================================================================!
   !  Append one record's data (calendar + counts + each live variable's scalar / slab).     !
   !=======================================================================================!
   subroutine write_one_record(stream, reg, pr, tier)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(pending_record_t),  intent(in)    :: pr
      integer(ik),             intent(in)    :: tier
      integer(c_int)    :: ncid
      integer(c_size_t) :: t0, i1(1), s2(2), c2(2), s3(3), c3(3)
      integer(ik)       :: j, k, ns, nsp
      ncid = int(stream%ncid, c_int)
      t0   = int(stream%nrec, c_size_t)
      i1   = [t0]
      !----- The cohort/patch axes were trimmed to the opening record's live count; every record in a  !
      !      ≤1-month file shares that count (§4.4). Assert it rather than silently overflow the slab.  !
      if (stream%has_cohort .and. pr%n_cohort > stream%cohort_dim)                                    &
         error stop 'meds_output_stream: cohort count grew within a file (cohort output needs file_chunk <= month)'
      if (stream%has_patch .and. pr%n_patch > stream%patch_dim)                                       &
         error stop 'meds_output_stream: patch count grew within a file (patch output needs file_chunk <= month)'
      !----- calendar (period start) + live counts. ---!
      call nc_check(nc_put_var1_double(ncid, int(stream%v_time, c_int), i1, time_to_decimal_year(pr%t_open)), 'put time')
      call put_int_rec(ncid, stream%v_year,  t0, int(pr%t_open%year,  c_int))
      call put_int_rec(ncid, stream%v_month, t0, int(pr%t_open%month, c_int))
      call put_int_rec(ncid, stream%v_day,   t0, int(pr%t_open%day,   c_int))
      if (stream%v_hour   >= 0_ik) call put_int_rec(ncid, stream%v_hour,   t0, int(pr%t_open%hour,   c_int))
      if (stream%v_minute >= 0_ik) call put_int_rec(ncid, stream%v_minute, t0, int(pr%t_open%minute, c_int))
      if (stream%v_second >= 0_ik) call put_int_rec(ncid, stream%v_second, t0, int(pr%t_open%second, c_int))
      if (stream%has_cohort) call put_int_rec(ncid, stream%v_ncohort, t0, int(pr%n_cohort, c_int))
      if (stream%has_patch)  call put_int_rec(ncid, stream%v_npatch,  t0, int(pr%n_patch,  c_int))

      do j = 1_ik, reg%nidx(tier)
         k  = reg%idx_freq(j, tier)
         ns = pr%nslab(k)
         if (reg%var(k)%dim == DIM_SCALAR) then
            if (reg%var(k)%xtype == XTYPE_INT) then
               call put_int_rec(ncid, stream%vid(k), t0, real_to_int(pr%sval(k), pr%svalid(k)))
            else
               call nc_check(nc_put_var1_double(ncid, int(stream%vid(k), c_int), i1, pr%sval(k)), 'put '//trim(reg%var(k)%name))
            end if
         else if (reg%var(k)%dim == DIM_SOIL_PATCH) then
            !----- Rank-3 (time, patch, soil). The flat slab already holds patch-major /          !
            !      layer-minor data strided by n_soil_layer_max, so the hyperslab count is just    !
            !      [1, n_patch, n_soil] over the same contiguous memory.  ------------------------!
            if (pr%n_patch > 0_ik) then
               s3 = [t0, 0_c_size_t, 0_c_size_t]
               c3 = [1_c_size_t, int(pr%n_patch, c_size_t), int(n_soil_layer_max, c_size_t)]
               nsp = pr%n_patch * n_soil_layer_max
               call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), s3, c3,          &
                             pr%slab(1:nsp, slab_col(pr, k))), 'put '//trim(reg%var(k)%name))
            end if
         else if (ns > 0_ik) then
            s2 = [t0, 0_c_size_t] ; c2 = [1_c_size_t, int(ns, c_size_t)]
            if (reg%var(k)%xtype == XTYPE_INT) then
               call put_int_slab(ncid, stream%vid(k), s2, c2, pr%slab(1:ns, slab_col(pr, k)),              &
                                 pr%slabvalid(1:ns, slab_col(pr, k)), ns)
            else
               call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), s2, c2,          &
                             pr%slab(1:ns, slab_col(pr, k))), 'put '//trim(reg%var(k)%name))
            end if
         end if
      end do
      stream%nrec = stream%nrec + 1_ik
   end subroutine write_one_record

   !=======================================================================================!
   !  REGION files (MEDS_POLYGON_RUNTIME_PLAN.md §6): one file per tier per time chunk for all     !
   !  polygons, with a `polygon` dimension after `time`. Record i of tier t is the same closed      !
   !  period in the buffers of every polygon that holds it (all polygons step the same calendar);   !
   !  each variable is packed from them into one region-wide array and written with ONE hyperslab.  !
   !  A polygon's slice holds exactly what its single-site file would: the slab entries past its    !
   !  live length, which the site writer leaves unwritten, carry the fill value here too, and so    !
   !  does the whole slice of a polygon that does not hold record i (it failed earlier this month). !
   !=======================================================================================!
   subroutine region_write_record(files, bufs, t, i)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(in)    :: bufs(:)
      integer(ik),           intent(in)    :: t, i
      integer(ik) :: bucket, np, p, p1
      np = size(bufs, kind=ik)
      if (np /= files%n_polygon) error stop 'region_write_record: one set of buffers per polygon of the region'
      p1 = 1_ik
      do while (bufs(p1)%queue(t)%n < i)
         p1 = p1 + 1_ik
         if (p1 > np) error stop 'region_write_record: no polygon holds the record'
      end do
      associate (r1 => bufs(p1)%queue(t)%rec(i))
         do p = p1 + 1_ik, np
            if (bufs(p)%queue(t)%n < i) cycle
            if (.not. same_time(bufs(p)%queue(t)%rec(i)%t_open, r1%t_open))                     &
               error stop 'region_write_record: the polygons'' records are not the same period'
         end do
         bucket = bucket_key(r1%t_open, files%file_chunk(t))
         if (files%stream(t)%ncid < 0_ik .or. bucket /= files%stream(t)%chunk_bucket) then
            call stream_close_file(files%stream(t))
            call region_open_file(files, t, r1, bucket)
         end if
      end associate
      call region_write_one(files, bufs, t, i, p1)
      if (files%sync_every == SYNC_FLUSH) call nc_check(nc_sync(int(files%stream(t)%ncid, c_int)), 'nc_sync')
   end subroutine region_write_record

   pure logical function same_time(a, b)
      type(meds_time_t), intent(in) :: a, b
      same_time = a%year == b%year .and. a%month == b%month .and. a%day == b%day .and.            &
                  a%hour == b%hour .and. a%minute == b%minute .and. a%second == b%second
   end function same_time

   !----- Create a region file: dims, polygon coordinates, the tier's variables, CF metadata. ----!
   subroutine region_open_file(files, tier, pr, bucket)
      type(output_files_t),   intent(inout) :: files
      integer(ik),            intent(in)    :: tier, bucket
      type(pending_record_t), intent(in)    :: pr
      character(len=512) :: path
      character(len=16)  :: stamp
      integer(c_int)     :: ncid, dt, dpo, ds, dpf, dsz, vid, dims1(1), dims2(2), dims3(3), xt, axislen
      integer(c_int)     :: v_id, v_lat, v_lon, v_row, v_col
      integer(c_size_t)  :: chunk2(2), chunk3(3), st1(1), cn1(1)
      integer(ik)        :: j, k, np
      logical            :: hass, haspf, hassz
      character(len=16)  :: cm
      associate (stream => files%stream(tier), reg => files%reg, dg => files%diag)
      np = files%n_polygon
      hass = .false. ; haspf = .false. ; hassz = .false.
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         select case (reg%var(k)%dim)
         case (DIM_SOIL) ; hass  = .true.
         case (DIM_PFT)  ; haspf = .true.
         case (DIM_SIZE) ; hassz = .true.
         case (DIM_COHORT, DIM_PATCH, DIM_SOIL_PATCH)
            error stop 'region_open_file: a region file holds only fixed-shape variables'
         end select
      end do

      stamp = chunk_stamp(pr%t_open, files%file_chunk(tier))
      if (len_trim(stamp) > 0) then
         path = trim(files%dir)//'/'//trim(files%prefix)//'-'//freq_letter(pr%freq)//'-'//trim(stamp)//'.nc'
      else
         path = trim(files%dir)//'/'//trim(files%prefix)//'-'//freq_letter(pr%freq)//'.nc'
      end if
      call nc_check(nc_create_f(trim(path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'region nc_create')
      call nc_check(nc_def_dim_f(ncid, 'time', NC_UNLIMITED, dt), 'dim time')
      call nc_check(nc_def_dim_f(ncid, 'polygon', int(np, c_size_t), dpo), 'dim polygon')
      ds = -1_c_int ; dpf = -1_c_int ; dsz = -1_c_int
      if (hass)  call nc_check(nc_def_dim_f(ncid, 'soil', int(n_soil_layer_max, c_size_t), ds), 'dim soil')
      if (haspf) call nc_check(nc_def_dim_f(ncid, 'pft', int(max(dg%n_pft,1_ik), c_size_t), dpf), 'dim pft')
      if (hassz) call nc_check(nc_def_dim_f(ncid, 'dbh_class',                                      &
                               int(max(dg%n_dbh_class,1_ik), c_size_t), dsz), 'dim dbh_class')

      stream%v_time  = int(def_scalar_var(ncid, dt, 'time',  NC_DOUBLE, 'year', 'decimal calendar year (period start)'), ik)
      stream%v_year  = int(def_scalar_var(ncid, dt, 'year',  NC_INT,    '1', 'calendar year (period start)'), ik)
      stream%v_month = int(def_scalar_var(ncid, dt, 'month', NC_INT,    '1', 'calendar month (period start)'), ik)
      stream%v_day   = int(def_scalar_var(ncid, dt, 'day',   NC_INT,    '1', 'calendar day (period start)'), ik)
      stream%v_hour = -1_ik ; stream%v_minute = -1_ik ; stream%v_second = -1_ik

      !----- The polygon axis: its id, its cell centre and its grid indices, so a region maps back  !
      !      onto the lat/lon grid in one step (OR7). ---------------------------------------------!
      dims1 = [dpo]
      call nc_check(nc_def_var_f(ncid, 'polygon_id', NC_INT, 1_c_int, dims1, v_id), 'def polygon_id')
      call put_var_text(ncid, v_id, 'long_name', 'forcing cell row-major index: row * n_lon + col')
      call nc_check(nc_def_var_f(ncid, 'lat', NC_DOUBLE, 1_c_int, dims1, v_lat), 'def lat')
      call put_var_text(ncid, v_lat, 'units', 'degrees_north')
      call put_var_text(ncid, v_lat, 'standard_name', 'latitude')
      call nc_check(nc_def_var_f(ncid, 'lon', NC_DOUBLE, 1_c_int, dims1, v_lon), 'def lon')
      call put_var_text(ncid, v_lon, 'units', 'degrees_east')
      call put_var_text(ncid, v_lon, 'standard_name', 'longitude')
      call nc_check(nc_def_var_f(ncid, 'row', NC_INT, 1_c_int, dims1, v_row), 'def row')
      call put_var_text(ncid, v_row, 'long_name', '0-based latitude index on the forcing grid')
      call nc_check(nc_def_var_f(ncid, 'col', NC_INT, 1_c_int, dims1, v_col), 'def col')
      call put_var_text(ncid, v_col, 'long_name', '0-based longitude index on the forcing grid')

      stream%v_pft = -1_ik ; stream%v_dbh_lower = -1_ik ; stream%v_dbh_upper = -1_ik ; stream%v_soil_z = -1_ik
      if (haspf) then
         dims1 = [dpf]
         call nc_check(nc_def_var_f(ncid, 'pft', NC_INT, 1_c_int, dims1, vid), 'def pft coord')
         call put_var_text(ncid, vid, 'long_name', 'plant functional type index')
         stream%v_pft = int(vid, ik)
      end if
      if (hass) then
         dims1 = [ds]
         call nc_check(nc_def_var_f(ncid, 'soil_z', NC_DOUBLE, 1_c_int, dims1, vid), 'def soil_z')
         call put_var_text(ncid, vid, 'units', 'm')
         call put_var_text(ncid, vid, 'long_name', 'soil layer node depth (negative downward)')
         stream%v_soil_z = int(vid, ik)
      end if
      if (hassz) then
         dims1 = [dsz]
         call nc_check(nc_def_var_f(ncid, 'dbh_lower', NC_DOUBLE, 1_c_int, dims1, vid), 'def dbh_lower')
         call put_var_text(ncid, vid, 'units', 'cm')
         call put_var_text(ncid, vid, 'long_name', 'DBH class lower edge (inclusive)')
         stream%v_dbh_lower = int(vid, ik)
         call nc_check(nc_def_var_f(ncid, 'dbh_upper', NC_DOUBLE, 1_c_int, dims1, vid), 'def dbh_upper')
         call put_var_text(ncid, vid, 'units', 'cm')
         call put_var_text(ncid, vid, 'long_name',                                                  &
              'DBH class upper edge (exclusive, except the last class which is closed)')
         stream%v_dbh_upper = int(vid, ik)
      end if

      !----- The tier's variables: (time, polygon) or (time, polygon, axis), one record's polygons  !
      !      per chunk, deflated like the site files. ------------------------------------------------!
      stream%vid = -1_ik
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         associate (v => reg%var(k))
            xt = merge(NC_INT, NC_DOUBLE, v%xtype == XTYPE_INT)
            if (v%dim == DIM_SCALAR) then
               dims2 = [dt, dpo] ; chunk2 = [1_c_size_t, int(np, c_size_t)]
               call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 2_c_int, dims2, vid), 'def '//trim(v%name))
               call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk2), 'chunk '//trim(v%name))
            else
               select case (v%dim)
               case (DIM_SOIL) ; dims3 = [dt, dpo, ds]  ; axislen = int(n_soil_layer_max, c_int)
               case (DIM_PFT)  ; dims3 = [dt, dpo, dpf] ; axislen = int(max(dg%n_pft, 1_ik), c_int)
               case default    ; dims3 = [dt, dpo, dsz] ; axislen = int(max(dg%n_dbh_class, 1_ik), c_int)
               end select
               chunk3 = [1_c_size_t, int(np, c_size_t), int(axislen, c_size_t)]
               call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 3_c_int, dims3, vid), 'def '//trim(v%name))
               call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk3), 'chunk '//trim(v%name))
            end if
            call nc_check(nc_def_var_deflate(ncid, vid, 1_c_int, 1_c_int, 4_c_int), 'deflate '//trim(v%name))
            call put_var_text(ncid, vid, 'units', trim(v%units))
            call put_var_text(ncid, vid, 'long_name', trim(v%long_name))
            cm = cell_methods_of(v%agg)
            call put_var_text(ncid, vid, 'cell_methods', trim(cm))
            if (v%xtype == XTYPE_INT) then
               call nc_check(nc_put_att_int_f(ncid, vid, '_FillValue', NC_INT, int(MISSING_INT, c_int)), 'fill')
            else
               call nc_check(nc_put_att_double_f(ncid, vid, '_FillValue', NC_DOUBLE,              &
                                                  real(MISSING_VALUE, c_double)), 'fill')
            end if
         end associate
         stream%vid(k) = int(vid, ik)
      end do

      call put_global(ncid, 'title', TITLE)
      call put_global(ncid, 'Conventions', 'CF-1.10')
      if (len_trim(files%forcing_qair) > 0) call put_global(ncid, 'forcing_qair', trim(files%forcing_qair))
      call nc_check(nc_enddef(ncid), 'region enddef')

      st1 = [0_c_size_t] ; cn1 = [int(np, c_size_t)]
      call nc_check(nc_put_vara_int(ncid, v_id, st1, cn1, int(files%polygon_id, c_int)), 'put polygon_id')
      call nc_check(nc_put_vara_double(ncid, v_lat, st1, cn1, real(files%polygon_lat, c_double)), 'put lat')
      call nc_check(nc_put_vara_double(ncid, v_lon, st1, cn1, real(files%polygon_lon, c_double)), 'put lon')
      call nc_check(nc_put_vara_int(ncid, v_row, st1, cn1, int(files%polygon_row, c_int)), 'put row')
      call nc_check(nc_put_vara_int(ncid, v_col, st1, cn1, int(files%polygon_col, c_int)), 'put col')
      if (hass)  call write_soil_coord(ncid, stream%v_soil_z, dg)
      if (haspf) call write_pft_coord(ncid, stream%v_pft, dg%n_pft)
      if (hassz) call write_size_coord(ncid, stream%v_dbh_lower, stream%v_dbh_upper, dg)

      stream%ncid = int(ncid, ik) ; stream%nrec = 0_ik ; stream%chunk_bucket = bucket
      stream%has_cohort = .false. ; stream%has_patch = .false. ; stream%has_soil = hass
      stream%has_pft = haspf ; stream%has_size = hassz
      stream%d_time = int(dt, ik) ; stream%d_polygon = int(dpo, ik)
      stream%d_soil = int(ds, ik) ; stream%d_pft = int(dpf, ik) ; stream%d_size = int(dsz, ik)
      write(*,'(2a)') ' output: ', trim(path)
      end associate
   end subroutine region_open_file

   !----- Append record i of tier t: the calendar, from polygon p1, the first that holds the record, !
   !      then each variable packed over the polygons, the fill value for one that does not. -------!
   subroutine region_write_one(files, bufs, t, i, p1)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(in)    :: bufs(:)
      integer(ik),           intent(in)    :: t, i, p1
      integer(c_int)    :: ncid
      integer(c_size_t) :: t0, i1(1), s2(2), c2(2), s3(3), c3(3)
      integer(ik)       :: j, k, np, p, n, ns, c, m
      real(c_double),  allocatable :: x1(:), x2(:,:)
      integer(c_int),  allocatable :: k1(:), k2(:,:)
      associate (stream => files%stream(t), reg => files%reg, r1 => bufs(p1)%queue(t)%rec(i))
      ncid = int(stream%ncid, c_int)
      t0 = int(stream%nrec, c_size_t) ; i1 = [t0]
      np = size(bufs, kind=ik)
      call nc_check(nc_put_var1_double(ncid, int(stream%v_time, c_int), i1, time_to_decimal_year(r1%t_open)), 'put time')
      call put_int_rec(ncid, stream%v_year,  t0, int(r1%t_open%year,  c_int))
      call put_int_rec(ncid, stream%v_month, t0, int(r1%t_open%month, c_int))
      call put_int_rec(ncid, stream%v_day,   t0, int(r1%t_open%day,   c_int))
      do j = 1_ik, reg%nidx(t)
         k = reg%idx_freq(j, t)
         if (reg%var(k)%dim == DIM_SCALAR) then
            s2 = [t0, 0_c_size_t] ; c2 = [1_c_size_t, int(np, c_size_t)]
            if (reg%var(k)%xtype == XTYPE_INT) then
               allocate(k1(np)) ; k1 = int(MISSING_INT, c_int)
               do p = 1_ik, np
                  if (bufs(p)%queue(t)%n < i) cycle
                  k1(p) = real_to_int(bufs(p)%queue(t)%rec(i)%sval(k), bufs(p)%queue(t)%rec(i)%svalid(k))
               end do
               call nc_check(nc_put_vara_int(ncid, int(stream%vid(k), c_int), s2, c2, k1), 'put '//trim(reg%var(k)%name))
               deallocate(k1)
            else
               allocate(x1(np)) ; x1 = real(MISSING_VALUE, c_double)
               do p = 1_ik, np
                  if (bufs(p)%queue(t)%n < i) cycle
                  x1(p) = bufs(p)%queue(t)%rec(i)%sval(k)
               end do
               call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), s2, c2, x1), 'put '//trim(reg%var(k)%name))
               deallocate(x1)
            end if
         else
            select case (reg%var(k)%dim)
            case (DIM_SOIL) ; n = n_soil_layer_max
            case (DIM_PFT)  ; n = max(files%diag%n_pft, 1_ik)
            case default    ; n = max(files%diag%n_dbh_class, 1_ik)
            end select
            s3 = [t0, 0_c_size_t, 0_c_size_t] ; c3 = [1_c_size_t, int(np, c_size_t), int(n, c_size_t)]
            if (reg%var(k)%xtype == XTYPE_INT) then
               allocate(k2(n, np)) ; k2 = int(MISSING_INT, c_int)
               do p = 1_ik, np
                  if (bufs(p)%queue(t)%n < i) cycle
                  associate (r => bufs(p)%queue(t)%rec(i))
                     ns = min(r%nslab(k), n) ; c = slab_col(r, k)
                     do m = 1_ik, ns
                        k2(m, p) = real_to_int(r%slab(m, c), r%slabvalid(m, c))
                     end do
                  end associate
               end do
               call nc_check(nc_put_vara_int(ncid, int(stream%vid(k), c_int), s3, c3, k2), 'put '//trim(reg%var(k)%name))
               deallocate(k2)
            else
               allocate(x2(n, np)) ; x2 = real(MISSING_VALUE, c_double)
               do p = 1_ik, np
                  if (bufs(p)%queue(t)%n < i) cycle
                  associate (r => bufs(p)%queue(t)%rec(i))
                     ns = min(r%nslab(k), n) ; c = slab_col(r, k)
                     if (ns > 0_ik) x2(1:ns, p) = r%slab(1:ns, c)
                  end associate
               end do
               call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), s3, c3, x2), 'put '//trim(reg%var(k)%name))
               deallocate(x2)
            end if
         end if
      end do
      stream%nrec = stream%nrec + 1_ik
      end associate
   end subroutine region_write_one

   !----- Map a normalized real (with validity) to an int record value (MISSING_INT if invalid). !
   pure integer(c_int) function real_to_int(x, valid) result(iv)
      real(wp), intent(in) :: x
      logical,  intent(in) :: valid
      if (valid) then ; iv = int(nint(x, ik), c_int) ; else ; iv = int(MISSING_INT, c_int) ; end if
   end function real_to_int

   subroutine put_int_rec(ncid, vid_ik, t0, val)
      integer(c_int),    intent(in) :: ncid
      integer(ik),       intent(in) :: vid_ik
      integer(c_size_t), intent(in) :: t0
      integer(c_int),    intent(in) :: val
      integer(c_size_t) :: st(1), cn(1)
      integer(c_int)    :: vv(1)
      st = [t0] ; cn = [1_c_size_t] ; vv = [val]
      call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), st, cn, vv), 'put int rec')
   end subroutine put_int_rec

   subroutine put_int_slab(ncid, vid_ik, s2, c2, x, valid, ns)
      integer(c_int),    intent(in) :: ncid
      integer(ik),       intent(in) :: vid_ik, ns
      integer(c_size_t), intent(in) :: s2(2), c2(2)
      real(wp),          intent(in) :: x(:)
      logical,           intent(in) :: valid(:)
      integer(c_int) :: iarr(ns)
      integer(ik)    :: i
      do i = 1_ik, ns
         if (valid(i)) then ; iarr(i) = int(nint(x(i), ik), c_int) ; else ; iarr(i) = int(MISSING_INT, c_int) ; end if
      end do
      call nc_check(nc_put_vara_int(ncid, int(vid_ik, c_int), s2, c2, iarr), 'put int slab')
   end subroutine put_int_slab

end module meds_output_stream
