!==========================================================================================!
! meds_netcdf_c -- a thin ISO_C_BINDING layer over the netCDF *C* library.                  !
!                                                                                          !
! netCDF-Fortran is not available for the ifx/nvfortran toolchains on this machine (its     !
! .mod files are gfortran-only and compiler-specific), but the netCDF C library is. C        !
! linkage is compiler-agnostic, so we bind the handful of C entry points we need directly.   !
! Only the subset used by meds_io is declared. Higher-level logic lives in meds_io.          !
!                                                                                          !
! Pitfalls baked in here: C is 0-based (callers pass 0-based start); nc_def_var lists dimIDs !
! slowest-varying FIRST (so [time, cohort]); sizes are C size_t (c_size_t); C strings are     !
! NUL-terminated (use cstr()).                                                               !
!==========================================================================================!
module meds_netcdf_c
   use iso_c_binding, only : c_int, c_size_t, c_double, c_char, c_ptr, c_null_char,           &
                             c_f_pointer, c_associated
   use meds_kinds,    only : ik
   implicit none
   public

   !----- netCDF constants (values from netcdf.h, v4.9.2). --------------------------------!
   integer(c_int), parameter :: NC_NOERR     = 0_c_int
   integer(c_int), parameter :: NC_NETCDF4   = 4096_c_int   ! 0x1000: HDF5-backed format
   integer(c_int), parameter :: NC_CLOBBER   = 0_c_int
   integer(c_int), parameter :: NC_NOWRITE   = 0_c_int      ! read-only open mode
   integer(c_int), parameter :: NC_GLOBAL    = -1_c_int
   integer(c_int), parameter :: NC_INT       = 4_c_int
   integer(c_int), parameter :: NC_DOUBLE    = 6_c_int
   integer(c_int), parameter :: NC_CHUNKED   = 0_c_int      ! storage mode for nc_def_var_chunking
   integer(c_size_t), parameter :: NC_UNLIMITED = 0_c_size_t

   interface
      integer(c_int) function nc_create(path, cmode, ncidp) bind(c, name="nc_create")
         import :: c_char, c_int
         character(kind=c_char), intent(in)  :: path(*)
         integer(c_int), value, intent(in)   :: cmode
         integer(c_int),        intent(out)  :: ncidp
      end function nc_create

      integer(c_int) function nc_def_dim(ncid, name, length, idp) bind(c, name="nc_def_dim")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,    intent(in)  :: ncid
         character(kind=c_char),   intent(in)  :: name(*)
         integer(c_size_t), value, intent(in)  :: length
         integer(c_int),           intent(out) :: idp
      end function nc_def_dim

      integer(c_int) function nc_def_var(ncid, name, xtype, ndims, dimids, varidp)           &
                              bind(c, name="nc_def_var")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_int), value,  intent(in)  :: xtype
         integer(c_int), value,  intent(in)  :: ndims
         integer(c_int),         intent(in)  :: dimids(*)
         integer(c_int),         intent(out) :: varidp
      end function nc_def_var

      integer(c_int) function nc_put_att_text(ncid, varid, name, length, tp)                 &
                              bind(c, name="nc_put_att_text")
         import :: c_char, c_int, c_size_t
         integer(c_int), value,    intent(in) :: ncid, varid
         character(kind=c_char),   intent(in) :: name(*)
         integer(c_size_t), value, intent(in) :: length
         character(kind=c_char),   intent(in) :: tp(*)
      end function nc_put_att_text

      integer(c_int) function nc_def_var_deflate(ncid, varid, shuffle, deflate, level)       &
                              bind(c, name="nc_def_var_deflate")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid, varid, shuffle, deflate, level
      end function nc_def_var_deflate

      integer(c_int) function nc_def_var_chunking(ncid, varid, storage, chunks)              &
                              bind(c, name="nc_def_var_chunking")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in) :: ncid, varid, storage
         integer(c_size_t),     intent(in) :: chunks(*)
      end function nc_def_var_chunking

      integer(c_int) function nc_enddef(ncid) bind(c, name="nc_enddef")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_enddef

      integer(c_int) function nc_put_vara_double(ncid, varid, startp, countp, op)            &
                              bind(c, name="nc_put_vara_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: startp(*), countp(*)
         real(c_double),        intent(in) :: op(*)
      end function nc_put_vara_double

      integer(c_int) function nc_put_vara_int(ncid, varid, startp, countp, op)               &
                              bind(c, name="nc_put_vara_int")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: startp(*), countp(*)
         integer(c_int),        intent(in) :: op(*)
      end function nc_put_vara_int

      integer(c_int) function nc_put_var1_double(ncid, varid, indexp, op)                     &
                              bind(c, name="nc_put_var1_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in) :: ncid, varid
         integer(c_size_t),     intent(in) :: indexp(*)
         real(c_double),        intent(in) :: op
      end function nc_put_var1_double

      integer(c_int) function nc_close(ncid) bind(c, name="nc_close")
         import :: c_int
         integer(c_int), value, intent(in) :: ncid
      end function nc_close

      !----- Read side (used by the state-restart reader, meds_io::io_read_state). ----------!
      integer(c_int) function nc_open(path, omode, ncidp) bind(c, name="nc_open")
         import :: c_char, c_int
         character(kind=c_char), intent(in)  :: path(*)
         integer(c_int), value,  intent(in)  :: omode
         integer(c_int),         intent(out) :: ncidp
      end function nc_open

      integer(c_int) function nc_inq_varid(ncid, name, varidp) bind(c, name="nc_inq_varid")
         import :: c_char, c_int
         integer(c_int), value,  intent(in)  :: ncid
         character(kind=c_char), intent(in)  :: name(*)
         integer(c_int),         intent(out) :: varidp
      end function nc_inq_varid

      integer(c_int) function nc_get_vara_int(ncid, varid, startp, countp, vals)             &
                              bind(c, name="nc_get_vara_int")
         import :: c_int, c_size_t
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_size_t),     intent(in)  :: startp(*), countp(*)
         integer(c_int),        intent(out) :: vals(*)
      end function nc_get_vara_int

      integer(c_int) function nc_get_vara_double(ncid, varid, startp, countp, vals)          &
                              bind(c, name="nc_get_vara_double")
         import :: c_int, c_size_t, c_double
         integer(c_int), value, intent(in)  :: ncid, varid
         integer(c_size_t),     intent(in)  :: startp(*), countp(*)
         real(c_double),        intent(out) :: vals(*)
      end function nc_get_vara_double

      type(c_ptr) function nc_strerror(ncerr) bind(c, name="nc_strerror")
         import :: c_int, c_ptr
         integer(c_int), value, intent(in) :: ncerr
      end function nc_strerror
   end interface

contains

   !----- Build a NUL-terminated C string from a Fortran string. --------------------------!
   pure function cstr(f) result(c)
      character(len=*), intent(in)  :: f
      character(kind=c_char) :: c(len_trim(f) + 1)
      integer(ik) :: i, n
      n = len_trim(f)
      do i = 1_ik, n
         c(i) = f(i:i)
      end do
      c(n + 1) = c_null_char
   end function cstr

   !---------------------------------------------------------------------------------------!
   ! String-argument wrappers. They bind cstr() to a NAMED local before the bind(c) call, so !
   ! an array-valued FUNCTION RESULT is never passed straight into a C call -- the nvfortran   !
   ! array-temp miscompile trap (CLAUDE.md issue #7: silently wrong at -O2). meds_io calls      !
   ! these `_f` forms with plain character(len=*) arguments; the raw bind(c) forms stay for      !
   ! the numeric calls (which already pass named array sections).                              !
   !---------------------------------------------------------------------------------------!
   integer(c_int) function nc_create_f(path, cmode, ncidp) result(st)
      character(len=*), intent(in)  :: path
      integer(c_int),   intent(in)  :: cmode
      integer(c_int),   intent(out) :: ncidp
      character(kind=c_char) :: cpath(len_trim(path) + 1)
      cpath = cstr(path)
      st = nc_create(cpath, cmode, ncidp)
   end function nc_create_f

   integer(c_int) function nc_def_dim_f(ncid, name, length, idp) result(st)
      integer(c_int),    intent(in)  :: ncid
      character(len=*),  intent(in)  :: name
      integer(c_size_t), intent(in)  :: length
      integer(c_int),    intent(out) :: idp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_def_dim(ncid, cname, length, idp)
   end function nc_def_dim_f

   integer(c_int) function nc_def_var_f(ncid, name, xtype, ndims, dimids, varidp) result(st)
      integer(c_int),   intent(in)  :: ncid, xtype, ndims
      character(len=*), intent(in)  :: name
      integer(c_int),   intent(in)  :: dimids(*)
      integer(c_int),   intent(out) :: varidp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_def_var(ncid, cname, xtype, ndims, dimids, varidp)
   end function nc_def_var_f

   integer(c_int) function nc_put_att_text_f(ncid, varid, name, length, tp) result(st)
      integer(c_int),    intent(in) :: ncid, varid
      character(len=*),  intent(in) :: name, tp
      integer(c_size_t), intent(in) :: length
      character(kind=c_char) :: cname(len_trim(name) + 1)
      character(kind=c_char) :: ctp(len_trim(tp) + 1)
      cname = cstr(name) ; ctp = cstr(tp)
      st = nc_put_att_text(ncid, varid, cname, length, ctp)
   end function nc_put_att_text_f

   integer(c_int) function nc_open_f(path, omode, ncidp) result(st)
      character(len=*), intent(in)  :: path
      integer(c_int),   intent(in)  :: omode
      integer(c_int),   intent(out) :: ncidp
      character(kind=c_char) :: cpath(len_trim(path) + 1)
      cpath = cstr(path)
      st = nc_open(cpath, omode, ncidp)
   end function nc_open_f

   integer(c_int) function nc_inq_varid_f(ncid, name, varidp) result(st)
      integer(c_int),   intent(in)  :: ncid
      character(len=*), intent(in)  :: name
      integer(c_int),   intent(out) :: varidp
      character(kind=c_char) :: cname(len_trim(name) + 1)
      cname = cstr(name)
      st = nc_inq_varid(ncid, cname, varidp)
   end function nc_inq_varid_f

   !----- Abort with the netCDF error message if `status` is not NC_NOERR. ----------------!
   subroutine nc_check(status, context)
      integer(c_int),   intent(in) :: status
      character(len=*), intent(in) :: context
      character(kind=c_char), pointer :: cbuf(:)
      type(c_ptr) :: cp
      integer     :: i
      if (status == NC_NOERR) return
      cp = nc_strerror(status)
      write(*,'(3a,i0)') 'meds_netcdf_c: netCDF error in ', trim(context), ', status=', status
      if (c_associated(cp)) then
         call c_f_pointer(cp, cbuf, [256])
         block
            character(len=256) :: text
            text = ''
            do i = 1, 256
               if (cbuf(i) == c_null_char) exit
               text(i:i) = cbuf(i)
            end do
            write(*,'(2a)') '   ', trim(text)
         end block
      end if
      error stop 'meds_netcdf_c: netCDF call failed'
   end subroutine nc_check

end module meds_netcdf_c
