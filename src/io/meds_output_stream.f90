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
   use meds_netcdf_c
   use meds_output_config, only : FC_DAY, FC_MONTH, FC_YEAR, FC_RUN, SYNC_FLUSH, freq_letter,     &
                                  freq_tier_index
   use meds_output_types,  only : output_registry_t, stream_file_t, pending_record_t, var_desc_t, &
                                  DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT,            &
                                  XTYPE_DOUBLE, XTYPE_INT, AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX,    &
                                  AGG_LAST, AGG_MEANSQ, AGG_TMEAN, AGG_FLUXSUM,                    &
                                  MISSING_VALUE, MISSING_INT
   implicit none
   private

   public :: stream_write_record, stream_close_file

   character(len=*), parameter :: TITLE = 'MEDS diagnostic aggregation output'

contains

   !----- Write one closed-period record, opening / rolling the per-tier file as needed. ------!
   subroutine stream_write_record(stream, reg, pr, dir, prefix, file_chunk, cohort_max,          &
                                  patch_max, sync_every)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, cohort_max, patch_max, sync_every
      integer(ik) :: tier, bucket
      tier   = freq_tier_index(pr%freq)
      bucket = bucket_key(pr%t_open, file_chunk)
      if (stream%ncid < 0_ik .or. bucket /= stream%chunk_bucket) then
         call stream_close_file(stream)
         call stream_open_file(stream, reg, pr, dir, prefix, file_chunk, bucket, cohort_max,      &
                               patch_max, tier)
      end if
      call write_one_record(stream, reg, pr, tier)
      if (sync_every == SYNC_FLUSH) call nc_check(nc_sync(int(stream%ncid, c_int)), 'nc_sync')
   end subroutine stream_write_record

   subroutine stream_close_file(stream)
      type(stream_file_t), intent(inout) :: stream
      if (stream%ncid >= 0_ik) call nc_check(nc_close(int(stream%ncid, c_int)), 'nc_close stream')
      stream%ncid = -1_ik ; stream%nrec = 0_ik ; stream%chunk_bucket = -1_ik
   end subroutine stream_close_file

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
   subroutine stream_open_file(stream, reg, pr, dir, prefix, file_chunk, bucket, cohort_max,      &
                               patch_max, tier)
      type(stream_file_t),     intent(inout) :: stream
      type(output_registry_t), intent(in)    :: reg
      type(pending_record_t),  intent(in)    :: pr
      character(len=*),        intent(in)    :: dir, prefix
      integer(ik),             intent(in)    :: file_chunk, bucket, cohort_max, patch_max, tier
      character(len=512) :: path
      character(len=16)  :: stamp
      character(len=1)   :: letter
      integer(c_int)     :: ncid, dt, dc, dp, vid
      integer(ik)        :: j, k
      logical            :: hasc, hasp

      !----- which ragged dims does this tier need? ---!
      hasc = .false. ; hasp = .false.
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         if (reg%var(k)%dim == DIM_COHORT) hasc = .true.
         if (reg%var(k)%dim == DIM_PATCH)  hasp = .true.
      end do

      letter = freq_letter(pr%freq)
      stamp  = chunk_stamp(pr%t_open, file_chunk)
      if (len_trim(stamp) > 0) then
         path = trim(dir)//'/'//trim(prefix)//'-'//letter//'-'//trim(stamp)//'.nc'
      else
         path = trim(dir)//'/'//trim(prefix)//'-'//letter//'.nc'
      end if

      call nc_check(nc_create_f(trim(path), ior(NC_NETCDF4, NC_CLOBBER), ncid), 'stream nc_create')
      call nc_check(nc_def_dim_f(ncid, 'time', NC_UNLIMITED, dt), 'dim time')
      dc = -1_c_int ; dp = -1_c_int
      if (hasc) call nc_check(nc_def_dim_f(ncid, 'cohort', int(cohort_max, c_size_t), dc), 'dim cohort')
      if (hasp) call nc_check(nc_def_dim_f(ncid, 'patch',  int(patch_max,  c_size_t), dp), 'dim patch')

      !----- time coordinate + calendar companions (period-start stamp). ---!
      stream%v_time  = int(def_scalar_var(ncid, dt, 'time',  NC_DOUBLE, 'year', 'decimal calendar year (period start)'), ik)
      stream%v_year  = int(def_scalar_var(ncid, dt, 'year',  NC_INT,    '1', 'calendar year (period start)'), ik)
      stream%v_month = int(def_scalar_var(ncid, dt, 'month', NC_INT,    '1', 'calendar month (period start)'), ik)
      stream%v_day   = int(def_scalar_var(ncid, dt, 'day',   NC_INT,    '1', 'calendar day (period start)'), ik)
      if (hasc) stream%v_ncohort = int(def_scalar_var(ncid, dt, 'n_cohort', NC_INT, '1', 'live cohorts this record'), ik)
      if (hasp) stream%v_npatch  = int(def_scalar_var(ncid, dt, 'n_patch',  NC_INT, '1', 'live patches this record'), ik)

      !----- registry-driven variable definitions (dims per DIM_*, chunk+deflate + CF attrs). ---!
      stream%vid = -1_ik
      do j = 1_ik, reg%nidx(tier)
         k = reg%idx_freq(j, tier)
         vid = def_registry_var(ncid, dt, dc, dp, reg%var(k), cohort_max, patch_max)
         stream%vid(k) = int(vid, ik)
      end do

      call put_global(ncid, 'title', TITLE)
      call put_global(ncid, 'Conventions', 'CF-1.10')
      call nc_check(nc_enddef(ncid), 'stream enddef')

      stream%ncid = int(ncid, ik) ; stream%nrec = 0_ik ; stream%chunk_bucket = bucket
      stream%has_cohort = hasc ; stream%has_patch = hasp
      stream%d_time = int(dt, ik) ; stream%d_cohort = int(dc, ik) ; stream%d_patch = int(dp, ik)
      write(*,'(2a)') ' output: ', trim(path)
   end subroutine stream_open_file

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
   integer(c_int) function def_registry_var(ncid, dt, dc, dp, v, cohort_max, patch_max) result(vid)
      integer(c_int),   intent(in) :: ncid, dt, dc, dp
      type(var_desc_t), intent(in) :: v
      integer(ik),      intent(in) :: cohort_max, patch_max
      integer(c_int)    :: xt, dims(2), dims1(1), axislen
      character(len=16) :: cm
      xt = merge(NC_INT, NC_DOUBLE, v%xtype == XTYPE_INT)
      if (v%dim == DIM_SCALAR) then
         dims1 = [dt]                       ! named local (avoids the nvfortran/ifx arg-temp trap, issue #7)
         call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 1_c_int, dims1, vid), 'def '//trim(v%name))
      else
         select case (v%dim)
         case (DIM_COHORT) ; dims = [dt, dc] ; axislen = int(cohort_max, c_int)
         case (DIM_PATCH)  ; dims = [dt, dp] ; axislen = int(patch_max,  c_int)
         case default      ; dims = [dt, dc] ; axislen = int(cohort_max, c_int)
         end select
         call nc_check(nc_def_var_f(ncid, trim(v%name), xt, 2_c_int, dims, vid), 'def '//trim(v%name))
         call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED,                                 &
                       [1_c_size_t, int(axislen, c_size_t)]), 'chunk '//trim(v%name))
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
      integer(c_size_t) :: t0, i1(1), s2(2), c2(2)
      integer(ik)       :: j, k, ns
      ncid = int(stream%ncid, c_int)
      t0   = int(stream%nrec, c_size_t)
      i1   = [t0]
      !----- calendar (period start) + live counts. ---!
      call nc_check(nc_put_var1_double(ncid, int(stream%v_time, c_int), i1, time_to_decimal_year(pr%t_open)), 'put time')
      call put_int_rec(ncid, stream%v_year,  t0, int(pr%t_open%year,  c_int))
      call put_int_rec(ncid, stream%v_month, t0, int(pr%t_open%month, c_int))
      call put_int_rec(ncid, stream%v_day,   t0, int(pr%t_open%day,   c_int))
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
