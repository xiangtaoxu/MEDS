!==========================================================================================!
! meds_io (STUB) -- compiled in place of src/io/meds_io.f90 when MEDS is built WITHOUT         !
! -DMEDS_ENABLE_IO (i.e. without the netCDF C library). It satisfies the SAME public interface !
! as the real writer (meds_io_t + io_create/io_write_snapshot/io_close) so the single entry     !
! point meds_main always compiles and runs; the routines simply do nothing (no file is          !
! written). io_create notes once that output is unavailable.                                   !
!                                                                                          !
! KEEP THIS PUBLIC INTERFACE IN SYNC with src/io/meds_io.f90 -- exactly one of the two is        !
! compiled per build (selected by MEDS_ENABLE_IO in CMakeLists.txt), so they must agree.        !
!==========================================================================================!
module meds_io
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t
   implicit none
   private

   public :: meds_io_t, io_create, io_write_snapshot, io_close

   type :: meds_io_t
      integer(ik) :: nrec = 0_ik
   end type meds_io_t

contains

   subroutine io_create(io, path, cfg)
      type(meds_io_t),     intent(out) :: io
      character(len=*),    intent(in)  :: path
      type(meds_config_t), intent(in)  :: cfg
      io%nrec = 0_ik
      write(*,'(3a)') ' note: built without netCDF output (-DMEDS_ENABLE_IO=OFF) -- "',         &
                      trim(path), '" will NOT be written.'
   end subroutine io_create

   subroutine io_write_snapshot(io, site, time_years)
      type(meds_io_t), intent(inout) :: io
      type(site_t),    intent(in)    :: site
      real(wp),        intent(in)    :: time_years
      io%nrec = io%nrec + 1_ik          ! count requested snapshots; nothing is written
   end subroutine io_write_snapshot

   subroutine io_close(io)
      type(meds_io_t), intent(inout) :: io
      io%nrec = 0_ik                    ! nothing to flush
   end subroutine io_close

end module meds_io
