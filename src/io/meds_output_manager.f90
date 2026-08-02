!==========================================================================================!
! meds_output_manager -- the serializer-side glue: drain the manager's pending stage to disk      !
! (output_serialize_pending, the ONLY flush -- called by main), and close the streams at run end   !
! (flushing any final partial period). netCDF via meds_output_stream. The netCDF-free half of the   !
! manager (buffer alloc = manager_alloc in meds_output_registry; the per-step tick = output_        !
! integrate in meds_output_integrate) is deliberately in the core library so the stepper stays off   !
! netCDF (§2, §4.5).                                                                                  !
!==========================================================================================!
module meds_output_manager
   use meds_kinds,            only : ik
   use meds_output_config,    only : N_FREQ
   use meds_output_types,     only : output_manager_t
   use meds_output_integrate, only : close_tier
   use meds_output_stream,    only : stream_write_record, stream_close_file
   implicit none
   private

   public :: output_serialize_pending, output_manager_close

contains

   !----- Drain every staged record to its per-tier file (the flush wall; main-only). --------!
   subroutine output_serialize_pending(mgr)
      type(output_manager_t), intent(inout) :: mgr
      integer(ik) :: t
      if (.not. mgr%enabled) return
      do t = 1_ik, N_FREQ
         if (.not. mgr%pending(t)%used) cycle
         call stream_write_record(mgr%stream(t), mgr%reg, mgr%diag, mgr%pending(t), mgr%dir,       &
                                  mgr%prefix,                                                    &
                                  mgr%file_chunk(t), mgr%cohort_max, mgr%patch_max, mgr%sync_every)
         mgr%pending(t)%used = .false.
      end do
   end subroutine output_serialize_pending

   !----- End of run: optionally flush each tier's final PARTIAL window, then close files. ----!
   subroutine output_manager_close(mgr, flush_partial)
      type(output_manager_t), intent(inout) :: mgr
      logical, optional,      intent(in)    :: flush_partial
      logical     :: fp
      integer(ik) :: t
      if (.not. mgr%enabled) return
      fp = .true. ; if (present(flush_partial)) fp = flush_partial
      if (fp) then
         do t = 1_ik, N_FREQ                 ! close every tier's final partial window (incl. FAST)
            if (mgr%has_data(t)) call close_tier(mgr, t)
         end do
         call output_serialize_pending(mgr)
      end if
      do t = 1_ik, N_FREQ
         call stream_close_file(mgr%stream(t))
      end do
   end subroutine output_manager_close

end module meds_output_manager
