!==========================================================================================!
! test_output_integrate -- the netCDF-free temporal-reduction kernels + the extract switchboard. !
! Covers MEDS_IO_DESIGN.md test 1 (integrator arithmetic + zero-sample guard) and test 3         !
! (fixed-within-window slab integration + extract_variable over a synthetic site).                !
!                                                                                          !
! nvfortran note: the slab test mutates the site's agb through the OPAQUE setter below (a          !
! separate module scope) between steps -- mirroring production, where state changes happen deep    !
! inside vegetation_dynamics/update_demography. Mutating an allocatable component INLINE between    !
! two extract_variable() calls lets nvfortran -O2 wrongly CSE the reads (the CLAUDE.md issue #7     !
! allocatable-component aliasing blind spot); an opaque mutation defeats it, as production does.     !
!==========================================================================================!
module test_output_integrate_support
   use meds_kinds,            only : wp, ik
   use meds_demography_types, only : site_t
   implicit none
contains
   !----- Opaque per-cohort agb setter (separate compilation scope; nvfortran cannot CSE across it). !
   subroutine set_cohort_agb(site, vals, n)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: vals(:)
      integer(ik),  intent(in)    :: n
      site%cohort%agb(1:n) = vals(1:n)
   end subroutine set_cohort_agb
end module test_output_integrate_support

program test_output_integrate
   use test_output_integrate_support, only : set_cohort_agb
   use meds_kinds,            only : wp, ik
   use meds_demography_types, only : site_t, site_alloc, site_free
   use meds_output_types,     only : var_desc_t, integ_buffer_t, MISSING_VALUE,                  &
                                     AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST, AGG_MEANSQ,   &
                                     AGG_TMEAN, AGG_FLUXSUM, DIM_SCALAR, DIM_COHORT
   use meds_output_integrate, only : alloc_integ_buffer, reset_buffer, integrate_scalar,         &
                                     integrate_slab, normalize_scalar, normalize_slab,           &
                                     extract_variable, SRC_C_AGB, SRC_S_AGB
   use meds_output_config,    only : FREQ_MONTHLY
   use meds_test_support,     only : check, check_close, banner
   implicit none

   call banner('output_integrate')
   call test_scalar_operators()
   call test_zero_sample_guard()
   call test_slab_and_extract()
   write(*,'(a)') 'test_output_integrate: ALL PASSED'

contains

   !----- Build a scalar buffer for one operator, fold a sequence, normalize. ---------------!
   subroutine run_scalar(agg, x, dt, out, valid, out2, has2)
      integer(ik), intent(in)  :: agg
      real(wp),    intent(in)  :: x(:), dt(:)
      real(wp),    intent(out) :: out
      logical,     intent(out) :: valid
      real(wp), optional, intent(out) :: out2
      logical,  optional, intent(out) :: has2
      type(var_desc_t)     :: v
      type(integ_buffer_t) :: buf
      integer(ik) :: i
      v%dim = DIM_SCALAR ; v%agg = agg
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      do i = 1_ik, int(size(x), ik)
         call integrate_scalar(buf, x(i), dt(i))
      end do
      call normalize_scalar(buf, out, valid, out2, has2)
   end subroutine run_scalar

   subroutine test_scalar_operators()
      real(wp) :: out, out2
      logical  :: valid, has2
      !----- MEAN (equal weight). -----!
      call run_scalar(AGG_MEAN, [2.0_wp,4.0_wp,6.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check(valid, 'MEAN valid'); call check_close(out, 4.0_wp, 1.0e-12_wp, 'AGG_MEAN')
      !----- TMEAN on NON-uniform dt (17.5, not the plain mean 15). -----!
      call run_scalar(AGG_TMEAN, [10.0_wp,20.0_wp], [1.0_wp,3.0_wp], out, valid)
      call check_close(out, 17.5_wp, 1.0e-12_wp, 'AGG_TMEAN dt-weighted')
      !----- FLUXSUM: constant rate 5 over 6 s -> 30. -----!
      call run_scalar(AGG_FLUXSUM, [5.0_wp,5.0_wp,5.0_wp], [2.0_wp,2.0_wp,2.0_wp], out, valid)
      call check_close(out, 30.0_wp, 1.0e-12_wp, 'AGG_FLUXSUM integral')
      !----- SUM (count-like). -----!
      call run_scalar(AGG_SUM, [1.0_wp,2.0_wp,3.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check_close(out, 6.0_wp, 1.0e-12_wp, 'AGG_SUM')
      !----- MIN / MAX. -----!
      call run_scalar(AGG_MIN, [3.0_wp,1.0_wp,2.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check_close(out, 1.0_wp, 1.0e-12_wp, 'AGG_MIN')
      call run_scalar(AGG_MAX, [3.0_wp,1.0_wp,2.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check_close(out, 3.0_wp, 1.0e-12_wp, 'AGG_MAX')
      !----- LAST. -----!
      call run_scalar(AGG_LAST, [7.0_wp,8.0_wp,9.0_wp], [1.0_wp,1.0_wp,1.0_wp], out, valid)
      call check_close(out, 9.0_wp, 1.0e-12_wp, 'AGG_LAST')
      !----- MEANSQ-derived variance (x=[2,4] -> mean 3, var 1). -----!
      call run_scalar(AGG_MEANSQ, [2.0_wp,4.0_wp], [1.0_wp,1.0_wp], out, valid, out2, has2)
      call check(has2, 'MEANSQ has variance'); call check_close(out, 3.0_wp, 1.0e-12_wp, 'MEANSQ mean')
      call check_close(out2, 1.0_wp, 1.0e-12_wp, 'MEANSQ variance')
   end subroutine test_scalar_operators

   subroutine test_zero_sample_guard()
      type(var_desc_t)     :: v
      type(integ_buffer_t) :: buf
      real(wp) :: out
      logical  :: valid
      !----- Never-touched MEAN buffer -> invalid, MISSING, no NaN / no huge leak. -----!
      v%dim = DIM_SCALAR ; v%agg = AGG_MEAN
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      call normalize_scalar(buf, out, valid)
      call check(.not. valid, 'zero-sample -> invalid')
      call check(out == MISSING_VALUE, 'zero-sample -> MISSING')
      !----- MIN re-seed after reset is +huge (no leak). -----!
      v%agg = AGG_MIN ; call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      call reset_buffer(buf)
      call check(buf%scal == huge(1.0_wp), 'MIN re-seeds to +huge')
      call normalize_scalar(buf, out, valid)
      call check(.not. valid .and. out == MISSING_VALUE, 'still-seeded MIN -> MISSING')
   end subroutine test_zero_sample_guard

   subroutine test_slab_and_extract()
      type(site_t)         :: site
      type(var_desc_t)     :: v_agb_c, v_agb_s
      type(integ_buffer_t) :: buf
      real(wp)             :: slab(64), out(64), scal
      logical              :: valid(64)
      integer(ik)          :: n_out, i
      !----- A fixed 3-cohort / 1-patch site (slot set constant across the window, §4.4). -----!
      call site_alloc(site, 2_ik, 64_ik, 8_ik, 4_ik)
      site%cohort%n = 3_ik
      site%cohort%nplant(1:3)     = [1.0_wp, 1.0_wp, 1.0_wp]
      call set_cohort_agb(site, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)   ! opaque (see header note)
      site%cohort%pft(1:3)        = [1_ik, 1_ik, 2_ik]
      site%cohort%owner_patch(1:3)= [1_ik, 1_ik, 1_ik]
      site%patch%n = 1_ik
      site%patch%area(1)          = 1.0_wp
      site%patch%cohort_offset(1) = 1_ik
      site%patch%cohort_count(1)  = 3_ik

      !----- extract_variable returns the live cohort slab + count. -----!
      v_agb_c%name = 'agb_cohort' ; v_agb_c%dim = DIM_COHORT ; v_agb_c%agg = AGG_TMEAN
      v_agb_c%source_id = SRC_C_AGB
      call extract_variable(site, v_agb_c, scal, slab, n_out)
      call check(n_out == 3_ik, 'extract cohort n')
      call check_close(slab(2), 2.0_wp, 1.0e-12_wp, 'extract cohort agb(2)')

      !----- Integrate the cohort slab over a fixed-slot window: TMEAN per slot, dt-weighted. -!
      !      step1 agb=[1,2,3] dt=1 ; step2 agb=[3,4,5] dt=3 -> slot1 = (1+9)/4 = 2.5. The agb   !
      !      mutation goes through the OPAQUE set_cohort_agb (production mutates via veg dynamics). !
      call alloc_integ_buffer(buf, v_agb_c, FREQ_MONTHLY, 64_ik)
      call set_cohort_agb(site, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)
      call extract_variable(site, v_agb_c, scal, slab, n_out)
      call integrate_slab(buf, slab, n_out, 1.0_wp)
      call set_cohort_agb(site, [3.0_wp, 4.0_wp, 5.0_wp], 3_ik)
      call extract_variable(site, v_agb_c, scal, slab, n_out)
      call integrate_slab(buf, slab, n_out, 3.0_wp)
      call normalize_slab(buf, out, valid, n_out)
      call check(n_out == 3_ik, 'slab n_out')
      call check(all(valid(1:3)), 'slab slots valid')
      call check_close(out(1), 2.5_wp, 1.0e-12_wp, 'slab TMEAN slot 1')
      call check_close(out(2), 3.5_wp, 1.0e-12_wp, 'slab TMEAN slot 2')
      call check_close(out(3), 4.5_wp, 1.0e-12_wp, 'slab TMEAN slot 3')
      !----- Slots past the live count are invalid (fill tail). -----!
      call check(.not. valid(4), 'fill-tail slot invalid')

      !----- Site scalar reduction via extract (total_agb = sum area*nplant*agb = 1*(1+2+3)). -!
      call set_cohort_agb(site, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)
      v_agb_s%name = 'agb_site' ; v_agb_s%dim = DIM_SCALAR ; v_agb_s%agg = AGG_MEAN
      v_agb_s%source_id = SRC_S_AGB
      call extract_variable(site, v_agb_s, scal, slab, n_out)
      call check(n_out == 0_ik, 'scalar extract n_out=0')
      call check_close(scal, 6.0_wp, 1.0e-12_wp, 'extract site agb')

      call site_free(site)
      if (.false.) i = 0_ik   ! silence unused
   end subroutine test_slab_and_extract

end program test_output_integrate
