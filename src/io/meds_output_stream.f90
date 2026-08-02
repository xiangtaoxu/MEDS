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
   use meds_column_state_types, only : n_soil_layer_max
   use meds_netcdf_c
   use meds_output_config, only : FC_DAY, FC_MONTH, FC_YEAR, FC_RUN, SYNC_FLUSH, freq_letter,     &
                                  freq_tier_index
   use meds_output_types,  only : output_registry_t, stream_file_t, pending_record_t, var_desc_t, &
                                  DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,            &
                                  DIM_SIZE, DIM_SOIL_PATCH, diag_params_t,                        &
                                  XTYPE_DOUBLE, XTYPE_INT, AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX,    &
                                  AGG_LAST, AGG_MEANSQ, AGG_TMEAN, AGG_FLUXSUM,                    &
                                  MISSING_VALUE, MISSING_INT
   implicit none
   private

   public :: stream_write_record, stream_close_file

   character(len=*), parameter :: TITLE = 'MEDS diagnostic aggregation output'

contains

   !----- Write one closed-period record, opening / rolling the per-tier file as needed. ------!
   subroutine stream_write_record(stream, reg, dg, pr, dir, prefix, file_chunk, cohort_max,      &
                                  patch_max, sync_every)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(diag_params_t),     intent(in)    :: dg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, cohort_max, patch_max, sync_every
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
                               patch_max, tier)
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
      case (AGG_MEAN, AGG_TMEAN, AGG_MEANSQ) ; cm = 'time: mean'
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
                               patch_max, tier)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(diag_params_t),     intent(in)    :: dg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, bucket, cohort_max, patch_max, tier
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
      if (haspf) then
         dims1 = [dpf]
         call nc_check(nc_def_var_f(ncid, 'pft', NC_INT, 1_c_int, dims1, vid), 'def pft coord')
         call put_var_text(ncid, vid, 'long_name', 'plant functional type index')
         stream%v_pft = int(vid, ik)
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
      call nc_check(nc_enddef(ncid), 'stream enddef')

      !----- Write the axis coordinates once, right after enddef (they do not vary by record). --!
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
                             pr%slab(1:nsp, k)), 'put '//trim(reg%var(k)%name))
            end if
         else if (ns > 0_ik) then
            s2 = [t0, 0_c_size_t] ; c2 = [1_c_size_t, int(ns, c_size_t)]
            if (reg%var(k)%xtype == XTYPE_INT) then
               call put_int_slab(ncid, stream%vid(k), s2, c2, pr%slab(1:ns, k), pr%slabvalid(1:ns, k), ns)
            else
               call nc_check(nc_put_vara_double(ncid, int(stream%vid(k), c_int), s2, c2,          &
                             pr%slab(1:ns, k)), 'put '//trim(reg%var(k)%name))
            end if
         end if
      end do
      stream%nrec = stream%nrec + 1_ik
   end subroutine write_one_record

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
