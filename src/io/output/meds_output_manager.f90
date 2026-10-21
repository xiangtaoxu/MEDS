! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_manager -- the serializer-side glue: drain queued records into the file set's       !
! streams (output_serialize, called only from the run's I/O phase), and close them at run end      !
! (flushing any final partial period). netCDF via meds_output_stream. The netCDF-free half of the  !
! manager -- allocation (setup_output_files / finalize_output_files / alloc_output_buffers in               !
! meds_output_registry) and the per-step tick (output_integrate in meds_output_integrate) -- is    !
! deliberately in meds_io_prep, so the stepper stays off netCDF (docs/science/diagnostics.md §1).   !
!==========================================================================================!
module meds_output_manager
   use meds_kinds,            only : ik
   use meds_output_params,    only : N_FREQ
   use meds_output_types,     only : output_files_t, output_buffers_t
   use meds_output_integrate, only : close_tier
   use meds_output_stream,    only : write_record, stream_close_file
   implicit none
   private

   public :: output_serialize, output_close

contains

   !----- Drain every queued record to its per-tier file, in closing order (the flush wall;       !
   !      main-only), then empty the queues. Called from the I/O phase at month boundaries and at    !
   !      the end of the run. `bufs` holds every polygon's buffers. Record i of tier t is one period !
   !      across the polygons; a polygon that failed earlier in the month closed fewer, so the phase !
   !      writes as many records as the longest queue holds, with the fill value where a polygon   !
   !      has none, and a failure loses no polygon's month. ------------------------------------------!
   subroutine output_serialize(files, bufs)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      integer(ik) :: t, i, p, nrec
      if (.not. files%enabled) return
      do t = 1_ik, N_FREQ
         nrec = 0_ik
         do p = 1_ik, size(bufs, kind=ik)
            nrec = max(nrec, bufs(p)%queue(t)%n)
         end do
         do i = 1_ik, nrec
            call write_record(files, bufs, t, i)
         end do
         do p = 1_ik, size(bufs, kind=ik)
            bufs(p)%queue(t)%n = 0_ik
         end do
      end do
   end subroutine output_serialize

   !----- End of run: optionally close every polygon's final PARTIAL windows (incl. FAST) and     !
   !      write them, then close the files. ----------------------------------------------------!
   subroutine output_close(files, bufs, flush_partial)
      type(output_files_t),   intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      logical, optional,      intent(in)    :: flush_partial
      logical     :: fp
      integer(ik) :: t, p
      if (.not. files%enabled) return
      fp = .true. ; if (present(flush_partial)) fp = flush_partial
      if (fp) then
         do p = 1_ik, size(bufs, kind=ik)
            do t = 1_ik, N_FREQ
               if (bufs(p)%has_data(t)) call close_tier(files, bufs(p), t)
            end do
         end do
         call output_serialize(files, bufs)
      end if
      do t = 1_ik, N_FREQ
         call stream_close_file(files%stream(t))
      end do
   end subroutine output_close

end module meds_output_manager
