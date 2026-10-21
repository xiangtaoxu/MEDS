! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_netcdf_write -- what the two netCDF writers share, the restart file                 !
! (meds_init_restart) and the diagnostic files (meds_output_stream): a variable defined    !
! with its attributes, chunks and compression in one call, text and integer attributes,    !
! the zlib levels, and the polygon axis both files lead with, its id and cell centre.      !
!==========================================================================================!
module meds_netcdf_write
   use iso_c_binding, only : c_int, c_size_t, c_double, c_char
   use meds_kinds,    only : wp, ik
   use meds_netcdf_c, only : cstr, nc_check, nc_def_dim_f, nc_def_var_f, nc_def_var_chunking, nc_def_var_deflate, &
                             nc_inq_varid_f, nc_put_att_int, nc_put_att_text_f, nc_put_vara_double,             &
                             nc_put_vara_int, NC_CHUNKED, NC_DOUBLE, NC_GLOBAL, NC_INT
   implicit none
   private

   public :: define_var, put_text_att, put_int_atts, put_global_text
   public :: define_polygon_axis, write_polygon_axis
   public :: RESTART_DEFLATE, OUTPUT_DEFLATE

   !----- The zlib levels: the restart file's variables, and the diagnostic files'. ------------!
   integer(c_int), parameter :: RESTART_DEFLATE = 1_c_int
   integer(c_int), parameter :: OUTPUT_DEFLATE  = 4_c_int

contains

   !----- Define a variable on `dims` (slowest first), with its units and long_name when given,  !
   !      its chunks, and its zlib level. -----------------------------------------------------!
   subroutine define_var(ncid, name, xtype, dims, vid, units, long_name, chunk, deflate)
      integer(c_int),    intent(in)  :: ncid, xtype, dims(:)
      character(len=*),  intent(in)  :: name
      integer(c_int),    intent(out) :: vid
      character(len=*),  intent(in), optional :: units, long_name
      integer(c_size_t), intent(in), optional :: chunk(:)
      integer(c_int),    intent(in), optional :: deflate
      call nc_check(nc_def_var_f(ncid, name, xtype, size(dims, kind=c_int), dims, vid), 'def '//name)
      if (present(chunk))   call nc_check(nc_def_var_chunking(ncid, vid, NC_CHUNKED, chunk), 'chunk '//name)
      if (present(deflate)) call nc_check(nc_def_var_deflate(ncid, vid, 1_c_int, 1_c_int, deflate), 'deflate '//name)
      if (present(units))     call put_text_att(ncid, vid, 'units', units)
      if (present(long_name)) call put_text_att(ncid, vid, 'long_name', long_name)
   end subroutine define_var

   !----- One text attribute on a defined variable. -------------------------------------------!
   subroutine put_text_att(ncid, vid, name, text)
      integer(c_int),   intent(in) :: ncid, vid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, vid, name, int(len_trim(text), c_size_t), text), 'att '//name)
   end subroutine put_text_att

   !----- One integer-array attribute on a defined variable. ----------------------------------!
   subroutine put_int_atts(ncid, vid, name, vals)
      integer(c_int),   intent(in) :: ncid, vid
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: vals(:)
      character(kind=c_char) :: cname(len_trim(name) + 1)
      integer(c_int)         :: iv(size(vals))
      cname = cstr(name)
      iv = int(vals, c_int)
      call nc_check(nc_put_att_int(ncid, vid, cname, NC_INT, int(size(iv), c_size_t), iv), 'att '//name)
   end subroutine put_int_atts

   !----- One global text attribute. ------------------------------------------------------------!
   subroutine put_global_text(ncid, name, text)
      integer(c_int),   intent(in) :: ncid
      character(len=*), intent(in) :: name, text
      call nc_check(nc_put_att_text_f(ncid, NC_GLOBAL, name, int(len_trim(text), c_size_t), text), 'global '//name)
   end subroutine put_global_text

   !----- The polygon axis: the `polygon` dimension, n_polygon long, and its id and cell centre, one !
   !      value per polygon; d_polygon is the dimension. ---------------------------------------!
   subroutine define_polygon_axis(ncid, n_polygon, d_polygon)
      integer(c_int), intent(in)  :: ncid
      integer(ik),    intent(in)  :: n_polygon
      integer(c_int), intent(out) :: d_polygon
      integer(c_int) :: v
      call nc_check(nc_def_dim_f(ncid, 'polygon', int(n_polygon, c_size_t), d_polygon), 'dim polygon')
      call define_var(ncid, 'polygon_id', NC_INT, [d_polygon], v,                                  &
                      long_name='forcing cell row-major index: row * n_lon + col')
      call define_var(ncid, 'lat', NC_DOUBLE, [d_polygon], v, units='degrees_north')
      call put_text_att(ncid, v, 'standard_name', 'latitude')
      call define_var(ncid, 'lon', NC_DOUBLE, [d_polygon], v, units='degrees_east')
      call put_text_att(ncid, v, 'standard_name', 'longitude')
   end subroutine define_polygon_axis

   !----- The polygon axis's values, every polygon at once, after the definitions end. ----------!
   subroutine write_polygon_axis(ncid, polygon_id, lat, lon)
      integer(c_int), intent(in) :: ncid
      integer(ik),    intent(in) :: polygon_id(:)
      real(wp),       intent(in) :: lat(:), lon(:)   !< [deg]
      integer(c_size_t) :: start1(1), count1(1)
      integer(c_int)    :: v, ids(size(polygon_id))
      real(c_double)    :: x(size(lat))
      start1 = 0_c_size_t
      count1 = int(size(polygon_id), c_size_t)
      !----- Each through a named array: no function result goes straight into a call (CLAUDE.md). -!
      ids = int(polygon_id, c_int)
      call nc_check(nc_inq_varid_f(ncid, 'polygon_id', v), 'inq polygon_id')
      call nc_check(nc_put_vara_int(ncid, v, start1, count1, ids), 'put polygon_id')
      x = real(lat, c_double)
      call nc_check(nc_inq_varid_f(ncid, 'lat', v), 'inq lat')
      call nc_check(nc_put_vara_double(ncid, v, start1, count1, x), 'put lat')
      x = real(lon, c_double)
      call nc_check(nc_inq_varid_f(ncid, 'lon', v), 'inq lon')
      call nc_check(nc_put_vara_double(ncid, v, start1, count1, x), 'put lon')
   end subroutine write_polygon_axis

end module meds_netcdf_write
