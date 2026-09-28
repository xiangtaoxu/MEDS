! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_manager -- the serializer-side glue: drain a polygon's queued records into its file   !
! set's streams (output_serialize_pending, the ONLY flush -- called by main), and close them at     !
! run end (flushing any final partial period). netCDF via meds_output_stream. The netCDF-free half  !
! of the manager (allocation = manager_setup / manager_finalize / manager_alloc_buffers in              !
! meds_output_registry; the per-step tick = output_integrate in meds_output_integrate) is            !
! deliberately in the core library so the stepper stays off netCDF (§2, §4.5).                        !
!==========================================================================================!
module meds_output_manager
   use meds_kinds,            only : ik
   use meds_output_config,    only : N_FREQ
   use meds_output_types,     only : output_files_t, output_buffers_t
   use meds_output_integrate, only : close_tier
   use meds_output_stream,    only : stream_write_record, stream_close_file, region_write_record
   implicit none
   private

   public :: output_serialize_pending, output_manager_close
   public :: output_serialize_region, output_region_close

contains

   !----- Drain every queued record to its per-tier file, in closing order (the flush wall;     !
   !      main-only). Called from the I/O phase at month boundaries and at the end of the run. ----!
   subroutine output_serialize_pending(files, bufs)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik) :: t, i
      if (.not. files%enabled) return
      do t = 1_ik, N_FREQ
         do i = 1_ik, bufs%queue(t)%n
            call stream_write_record(files%stream(t), files%reg, files%diag, bufs%queue(t)%rec(i), files%dir,  &
                                     files%prefix,                                                  &
                                     files%file_chunk(t), files%cohort_max, files%patch_max, files%sync_every, &
                                     files%forcing_qair)
         end do
         bufs%queue(t)%n = 0_ik
      end do
   end subroutine output_serialize_pending

   !----- End of run: optionally flush each tier's final PARTIAL window, then close files. ----!
   subroutine output_manager_close(files, bufs, flush_partial)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs
      logical, optional,     intent(in)    :: flush_partial
      logical     :: fp
      integer(ik) :: t
      if (.not. files%enabled) return
      fp = .true. ; if (present(flush_partial)) fp = flush_partial
      if (fp) then
         do t = 1_ik, N_FREQ                 ! close every tier's final partial window (incl. FAST)
            if (bufs%has_data(t)) call close_tier(files, bufs, t)
         end do
         call output_serialize_pending(files, bufs)
      end if
      do t = 1_ik, N_FREQ
         call stream_close_file(files%stream(t))
      end do
   end subroutine output_manager_close

   !----- A REGION's I/O phase: record i of tier t is one period across the polygons; write it to the !
   !      region file as one hyperslab per variable (MEDS_POLYGON_RUNTIME_PLAN.md §6.1), in closing   !
   !      order, then empty every queue. The polygons close the same periods, except in a month one   !
   !      of them failed: it closed fewer, and the polygons after it none. The phase writes what       !
   !      closed -- as many records as the longest queue holds, the fill value where a polygon has    !
   !      none -- so a failure loses no polygon's month. ----------------------------------------------!
   subroutine output_serialize_region(files, bufs)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      integer(ik) :: t, i, p, nrec
      if (.not. files%enabled) return
      do t = 1_ik, N_FREQ
         nrec = 0_ik
         do p = 1_ik, size(bufs, kind=ik)
            nrec = max(nrec, bufs(p)%queue(t)%n)
         end do
         do i = 1_ik, nrec
            call region_write_record(files, bufs, t, i)
         end do
         do p = 1_ik, size(bufs, kind=ik)
            bufs(p)%queue(t)%n = 0_ik
         end do
      end do
   end subroutine output_serialize_region

   !----- End of a region run: close every polygon's final partial windows, write them, and close  !
   !      the region files. --------------------------------------------------------------------!
   subroutine output_region_close(files, bufs)
      type(output_files_t),  intent(inout) :: files
      type(output_buffers_t), intent(inout) :: bufs(:)
      integer(ik) :: t, p
      if (.not. files%enabled) return
      do p = 1_ik, size(bufs, kind=ik)
         do t = 1_ik, N_FREQ
            if (bufs(p)%has_data(t)) call close_tier(files, bufs(p), t)
         end do
      end do
      call output_serialize_region(files, bufs)
      do t = 1_ik, N_FREQ
         call stream_close_file(files%stream(t))
      end do
   end subroutine output_region_close

end module meds_output_manager
