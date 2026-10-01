! SPDX-License-Identifier: Apache-2.0
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
   use meds_site_state_types, only : site_t
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
   use meds_config,           only : meds_config_t
   use meds_site_state_types, only : site_t, site_alloc, site_free
   use meds_output_types,     only : var_desc_t, integ_buffer_t, output_files_t, output_buffers_t,  &
                                     diag_params_t, slab_col,                                     &
                                     MISSING_VALUE,                                               &
                                     AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST,                 &
                                     AGG_TMEAN, AGG_FLUXSUM, DIM_SCALAR, DIM_COHORT
   use meds_output_integrate, only : alloc_integ_buffer, reset_buffer, integrate_scalar,         &
                                     integrate_slab, normalize_scalar, normalize_slab,           &
                                     extract_variable, output_integrate_fast, close_tier,        &
                                     extract_fast_scalar, FLD_C_AGB, SRC_F_PD0, SRC_F_PY0,       &
                                     output_integrate
   use meds_site_diag_types,  only : N_PYDIAG, PY_TAIR, PY_CO2, N_PDIAG, PD_CAS_TEMP, PD_LE, PD_H, &
                                     PD_GPP, PD_SW_UP, PD_LW_UP
   use meds_time,             only : meds_time_t, time_advance_days
   use meds_output_registry,  only : manager_alloc, manager_alloc_buffers, find_var_index,        &
                                     manager_setup, manager_finalize, build_freq_index,          &
                                     apply_variable_override, OVR_MASK
   use meds_diagnostic_reduce, only : W_NPLANT
   use meds_output_config,    only : FREQ_MONTHLY, FREQ_FAST, FREQ_NONE
   use meds_column_params,    only : n_soil_layer_max
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   call banner('output_integrate')
   call test_scalar_operators()
   call test_zero_sample_guard()
   call test_slab_and_extract()
   call test_fast_tier()
   call test_slab_sized_after_overrides()
   call test_soil_patch_record()
   call test_two_buffers()
   call test_boundary_step()
   write(*,'(a)') 'test_output_integrate: ALL PASSED'

contains

   !----- Build a scalar buffer for one operator, fold a sequence, normalize. ---------------!
   subroutine run_scalar(agg, x, dt, out, valid)
      integer(ik), intent(in)  :: agg
      real(wp),    intent(in)  :: x(:), dt(:)
      real(wp),    intent(out) :: out
      logical,     intent(out) :: valid
      type(var_desc_t)     :: v
      type(integ_buffer_t) :: buf
      integer(ik) :: i
      v%dim = DIM_SCALAR ; v%agg = agg
      call alloc_integ_buffer(buf, v, FREQ_MONTHLY, 1_ik)
      do i = 1_ik, int(size(x), ik)
         call integrate_scalar(buf, x(i), dt(i))
      end do
      call normalize_scalar(buf, out, valid)
   end subroutine run_scalar

   subroutine test_scalar_operators()
      real(wp) :: out
      logical  :: valid
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
      logical              :: valid(64), vslab(64)
      type(diag_params_t)  :: dp
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
      v_agb_c%source_id = FLD_C_AGB
      call extract_variable(site, dp, v_agb_c, scal, slab, vslab, n_out)
      call check(n_out == 3_ik, 'extract cohort n')
      call check_close(slab(2), 2.0_wp, 1.0e-12_wp, 'extract cohort agb(2)')

      !----- Integrate the cohort slab over a fixed-slot window: TMEAN per slot, dt-weighted. -!
      !      step1 agb=[1,2,3] dt=1 ; step2 agb=[3,4,5] dt=3 -> slot1 = (1+9)/4 = 2.5. The agb   !
      !      mutation goes through the OPAQUE set_cohort_agb (production mutates via veg dynamics). !
      call alloc_integ_buffer(buf, v_agb_c, FREQ_MONTHLY, 64_ik)
      call set_cohort_agb(site, [1.0_wp, 2.0_wp, 3.0_wp], 3_ik)
      call extract_variable(site, dp, v_agb_c, scal, slab, vslab, n_out)
      call integrate_slab(buf, slab, n_out, 1.0_wp)
      call set_cohort_agb(site, [3.0_wp, 4.0_wp, 5.0_wp], 3_ik)
      call extract_variable(site, dp, v_agb_c, scal, slab, vslab, n_out)
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
      !----- The site twin is the SAME field with a different dim + weight: agb is EXTENSIVE, so   !
      !      it reduces as a weighted SUM over nplant (W_NPLANT), landing in kgC/m2 ground.  --------!
      v_agb_s%name = 'agb_site' ; v_agb_s%dim = DIM_SCALAR ; v_agb_s%agg = AGG_MEAN
      v_agb_s%source_id = FLD_C_AGB ; v_agb_s%weight = W_NPLANT ; v_agb_s%mean = .false.
      call extract_variable(site, dp, v_agb_s, scal, slab, vslab, n_out)
      call check(n_out == 0_ik, 'scalar extract n_out=0')
      call check_close(scal, 6.0_wp, 1.0e-12_wp, 'extract site agb')

      call site_free(site)
      if (.false.) i = 0_ik   ! silence unused
   end subroutine test_slab_and_extract

   !----- The FAST (sub-daily) tier: extract_fast_scalar id mapping, then a 2-sub-step window folded  !
   !      by output_integrate_fast + close_tier(1) -> TMEAN in pending(1), across the scalar, soil-     !
   !      column, and per-cohort slab paths. netCDF-free; mirrors what main's replay loop does.  -------!
   subroutine test_fast_tier()
      type(meds_config_t)    :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      real(wp)    :: f0(N_PYDIAG), s(N_PDIAG)
      integer(ik) :: k_cas, k_soil, k_leaf, k_air, k_hgt, k_cas_p, k_soil_p, nl
      real(wp), parameter :: DT = 900.0_wp    ! uniform sub-step -> TMEAN == plain mean

      !----- extract_fast_scalar: a patch-block id reads its row of the sub-step's site mean, and a  !
      !      forcing id its row of the sub-step's forcing table. -----------------------------------!
      f0 = 0.0_wp ; f0(PY_CO2) = 415.0_wp
      s = 0.0_wp
      s(PD_CAS_TEMP) = 290.0_wp ; s(PD_LE) = 100.0_wp ; s(PD_H) = 50.0_wp ; s(PD_GPP) = 12.0_wp
      s(PD_SW_UP) = 83.0_wp ; s(PD_LW_UP) = 455.0_wp
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_CAS_TEMP, s, f0), 290.0_wp, 1.0e-12_wp, 'fast extract cas_temp')
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_LE,       s, f0), 100.0_wp, 1.0e-12_wp, 'fast extract le')
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_H,        s, f0),  50.0_wp, 1.0e-12_wp, 'fast extract h')
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_GPP,      s, f0),  12.0_wp, 1.0e-12_wp, 'fast extract gpp_rate')
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_SW_UP,    s, f0),  83.0_wp, 1.0e-12_wp, 'fast extract sw_up')
      call check_close(extract_fast_scalar(SRC_F_PD0 + PD_LW_UP,    s, f0), 455.0_wp, 1.0e-12_wp, 'fast extract lw_up')
      call check_close(extract_fast_scalar(SRC_F_PY0 + PY_CO2, s, f0), 415.0_wp, 1.0e-12_wp, 'fast extract forcing CO2')

      !----- Build a manager with the FAST tier + all groups ON (so the energy/water/carbon FAST vars   !
      !      register), then fold a 2-sub-step window and close it. -----!
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      cfg%output%grp_on     = .true.
      cfg%output%axis_on    = .true.               ! the soil-by-patch axis too
      cfg%output%cohort_max = 8_ik
      call manager_alloc(files, bufs, cfg)
      call check(files%reg%nidx(1) > 0_ik, 'FAST tier has live variables')

      !----- Stage 2 sub-steps of two patches with known values (as fast_dynamics would). -----!
      nl = n_soil_layer_max
      call alloc_fast_staging(bufs, 2_ik, 2_ik, 8_ik)
      bufs%fast_n_soil = 2_ik ; bufs%fast_n_cohort = 2_ik
      bufs%fast_site(PD_CAS_TEMP, :) = [290.0_wp, 294.0_wp]                 ! mean 292
      bufs%fast_site(PD_LE, :)       = [100.0_wp, 200.0_wp]
      bufs%fast_patch(PD_CAS_TEMP, 1, :) = [290.0_wp, 294.0_wp]             ! patch 1 mean 292
      bufs%fast_patch(PD_CAS_TEMP, 2, :) = [300.0_wp, 304.0_wp]             ! patch 2 mean 302
      bufs%fast_forcing(PY_TAIR, 1) = 300.0_wp ; bufs%fast_forcing(PY_TAIR, 2) = 302.0_wp    ! mean 301
      bufs%fast_soil_temp(1:2,1) = [280.0_wp, 281.0_wp]                     ! slot1 mean 281
      bufs%fast_soil_temp(1:2,2) = [282.0_wp, 283.0_wp]                     ! slot2 mean 282
      bufs%fast_soil_temp_patch(1:2,1)       = [280.0_wp, 281.0_wp]         ! patch 1: 281, 282
      bufs%fast_soil_temp_patch(1:2,2)       = [282.0_wp, 283.0_wp]
      bufs%fast_soil_temp_patch(nl+1:nl+2,1) = [290.0_wp, 291.0_wp]         ! patch 2: 291, 292
      bufs%fast_soil_temp_patch(nl+1:nl+2,2) = [292.0_wp, 293.0_wp]
      bufs%fast_coh_ltemp(1:2,1) = [288.0_wp, 289.0_wp]                     ! slot1 mean 289
      bufs%fast_coh_ltemp(1:2,2) = [290.0_wp, 291.0_wp]                     ! slot2 mean 290
      bufs%fast_coh_height(1:2,1) = [10.0_wp, 12.0_wp]                      ! slot1=10, slot2=12 (constant)
      bufs%fast_coh_height(1:2,2) = [10.0_wp, 12.0_wp]

      call output_integrate_fast(files, bufs, 1_ik, DT)
      call output_integrate_fast(files, bufs, 2_ik, DT)
      call close_tier(files, bufs, 1_ik)

      !----- Scalar path: cas_temp_fast = TMEAN(290,294) = 292. The FAST-tier variables are        !
      !      DISTINCT registry entries from their coarse-tier namesakes (cas_temp_site), because    !
      !      they read the staged sub-step samples rather than live site state.  -------------------!
      k_cas = find_var_index(files%reg, 'cas_temp_fast')
      call check(k_cas > 0_ik, 'cas_temp_fast registered')
      call check(bufs%pending(1)%svalid(k_cas), 'FAST cas_temp valid after close')
      call check_close(bufs%pending(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'FAST cas_temp TMEAN')
      !----- Soil-column slab path: soil_temp_site slot means. -----!
      k_soil = find_var_index(files%reg, 'soil_temp_site_fast')
      call check(bufs%pending(1)%nslab(k_soil) == 2_ik, 'FAST soil slab length 2')
      call check_close(bufs%pending(1)%slab(1,k_soil), 281.0_wp, 1.0e-10_wp, 'FAST soil_temp slot 1')
      call check_close(bufs%pending(1)%slab(2,k_soil), 282.0_wp, 1.0e-10_wp, 'FAST soil_temp slot 2')
      !----- Per-cohort slab path (P2): leaf_temp_cohort_fast slot means + n_cohort. -----!
      k_leaf = find_var_index(files%reg, 'leaf_temp_cohort_fast')
      call check(bufs%pending(1)%n_cohort == 2_ik, 'FAST n_cohort = 2')
      call check_close(bufs%pending(1)%slab(1,k_leaf), 289.0_wp, 1.0e-10_wp, 'FAST leaf_temp cohort 1')
      call check_close(bufs%pending(1)%slab(2,k_leaf), 290.0_wp, 1.0e-10_wp, 'FAST leaf_temp cohort 2')
      !----- Forcing air-temp scalar path + per-cohort HEIGHT slab path (tallest-cohort post-proc). -----!
      k_air = find_var_index(files%reg, 'air_temp_fast')
      call check(k_air > 0_ik, 'air_temp_fast registered')
      call check_close(bufs%pending(1)%sval(k_air), 301.0_wp, 1.0e-10_wp, 'FAST air_temp TMEAN')
      k_hgt = find_var_index(files%reg, 'height_cohort_fast')
      call check(k_hgt > 0_ik, 'height_cohort_fast registered')
      call check_close(bufs%pending(1)%slab(1,k_hgt), 10.0_wp, 1.0e-10_wp, 'FAST height cohort 1')
      call check_close(bufs%pending(1)%slab(2,k_hgt), 12.0_wp, 1.0e-10_wp, 'FAST height cohort 2')
      !----- The patch axis (#270): each patch's own mean, and each patch's soil column in the      !
      !      soil-by-patch layout, with the inactive layers as fill. -----!
      k_cas_p = find_var_index(files%reg, 'cas_temp_patch_fast')
      call check(k_cas_p > 0_ik, 'cas_temp_patch_fast registered')
      call check(bufs%pending(1)%n_patch == 2_ik, 'FAST n_patch = 2')
      call check_close(bufs%pending(1)%slab(1,k_cas_p), 292.0_wp, 1.0e-10_wp, 'FAST cas_temp patch 1')
      call check_close(bufs%pending(1)%slab(2,k_cas_p), 302.0_wp, 1.0e-10_wp, 'FAST cas_temp patch 2')
      k_soil_p = find_var_index(files%reg, 'soil_temp_layer_patch_fast')
      call check(k_soil_p > 0_ik, 'soil_temp_layer_patch_fast registered')
      call check_close(bufs%pending(1)%slab(1,k_soil_p),    281.0_wp, 1.0e-10_wp, 'FAST soil patch 1 layer 1')
      call check_close(bufs%pending(1)%slab(2,k_soil_p),    282.0_wp, 1.0e-10_wp, 'FAST soil patch 1 layer 2')
      call check_close(bufs%pending(1)%slab(nl+1,k_soil_p), 291.0_wp, 1.0e-10_wp, 'FAST soil patch 2 layer 1')
      call check_close(bufs%pending(1)%slab(nl+2,k_soil_p), 292.0_wp, 1.0e-10_wp, 'FAST soil patch 2 layer 2')
      call check(.not. bufs%pending(1)%slabvalid(3,k_soil_p), 'FAST soil patch 1 layer 3 is fill (inactive)')
   end subroutine test_fast_tier

   !----- The driver's order: manager_setup, then an [output].io_config override that switches on a  !
   !      slab variable setup left off, then manager_finalize and manager_alloc_buffers. The slab must !
   !      hold every soil layer, read from the QUEUED record through slab_col, as the writer reads it. !
   subroutine test_slab_sized_after_overrides()
      type(meds_config_t)    :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      integer(ik) :: k, c, i, nl, n_live_slab
      logical     :: found, every_layer
      nl = n_soil_layer_max
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      cfg%output%grp_on     = .false.              ! every group off: no slab variable is live yet
      call manager_setup(files, cfg)
      n_live_slab = 0_ik
      do k = 1_ik, files%reg%nvar
         if (files%reg%var(k)%dim /= DIM_SCALAR .and. files%reg%var(k)%streams /= FREQ_NONE)        &
            n_live_slab = n_live_slab + 1_ik
      end do
      call check(n_live_slab == 0_ik, 'override test premise: setup leaves no slab variable live')

      call apply_variable_override(files%reg, 'soil_temp_site_fast', OVR_MASK, FREQ_FAST, found)
      call check(found, 'soil_temp_site_fast is a registry variable')
      call build_freq_index(files%reg)             ! what apply_io_overrides does after its overrides
      call manager_finalize(files)
      call manager_alloc_buffers(files, bufs)
      call check(files%max_slab >= nl, 'max_slab covers the soil axis an override switched on')

      call alloc_fast_staging(bufs, 1_ik, 0_ik, 1_ik)
      bufs%fast_n_soil = nl ; bufs%fast_n_cohort = 0_ik
      bufs%fast_soil_temp(:,1) = [(270.0_wp + real(i, wp), i = 1_ik, nl)]
      call output_integrate_fast(files, bufs, 1_ik, 900.0_wp)
      call close_tier(files, bufs, 1_ik)

      k = find_var_index(files%reg, 'soil_temp_site_fast')
      associate (r => bufs%queue(1)%rec(bufs%queue(1)%n))
         c = slab_col(r, k)
         call check(r%nslab(k) == nl, 'the queued soil slab has every layer')
         call check(size(r%slab, 1) >= nl, 'the queued record has a row for every layer')
         every_layer = size(r%slab, 1) >= nl
         if (every_layer) every_layer = all(r%slab(1:nl, c) == [(270.0_wp + real(i, wp), i = 1_ik, nl)])
         call check(every_layer, 'every soil layer reads back through slab_col, as the writer reads it')
      end associate
   end subroutine test_slab_sized_after_overrides

   !----- A tier whose only patch output is soil by patch still records its patch count: the slab  !
   !      is patch-major, so it holds nslab / n_soil_layer_max patches. Before, only a DIM_PATCH    !
   !      variable set the count, and such a record wrote with one patch and every value as fill. -!
   subroutine test_soil_patch_record()
      type(meds_config_t)    :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      integer(ik) :: k, nl
      logical     :: found
      nl = n_soil_layer_max
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.               ! FAST tier on
      cfg%output%grp_on     = .false.              ! every group off: no slab variable is live yet
      call manager_setup(files, cfg)
      call apply_variable_override(files%reg, 'soil_temp_layer_patch_fast', OVR_MASK, FREQ_FAST, found)
      call check(found, 'soil_temp_layer_patch_fast is a registry variable')
      call build_freq_index(files%reg)
      call manager_finalize(files)
      call manager_alloc_buffers(files, bufs)
      call alloc_fast_staging(bufs, 1_ik, 3_ik, 1_ik)
      bufs%fast_n_soil = 2_ik
      bufs%fast_soil_temp_patch(:,1) = 285.0_wp
      call output_integrate_fast(files, bufs, 1_ik, 900.0_wp)
      call close_tier(files, bufs, 1_ik)
      k = find_var_index(files%reg, 'soil_temp_layer_patch_fast')
      associate (r => bufs%queue(1)%rec(bufs%queue(1)%n))
         call check(r%n_patch == 3_ik, 'a soil-by-patch-only record counts its patches')
         call check(r%nslab(k) == 3_ik * nl, 'the soil-by-patch slab holds every patch''s column')
      end associate
   end subroutine test_soil_patch_record

   !----- Two polygons' buffers for one file set (MEDS_POLYGON_RUNTIME_PLAN.md R2): each          !
   !      reduces only its own samples, and closing one leaves the other's window open. Folds are  !
   !      interleaved the way a region's month loop would interleave its polygons.  ----------------!
   subroutine test_two_buffers()
      type(meds_config_t)   :: cfg
      type(output_files_t)  :: files
      type(output_buffers_t) :: a, b
      integer(ik) :: k_cas
      real(wp), parameter :: DT = 900.0_wp
      cfg = build_test_config(86400.0_wp)
      cfg%output%enabled    = .true.
      cfg%output%freq_on(1) = .true.
      cfg%output%grp_on     = .true.
      cfg%output%cohort_max = 8_ik
      call manager_alloc(files, a, cfg)
      call manager_alloc_buffers(files, b)
      call check(a%fast_on .and. b%fast_on, 'both polygons'' buffers stage the FAST tier')
      call stage_cas(a, [290.0_wp, 294.0_wp])                              ! mean 292
      call stage_cas(b, [270.0_wp, 280.0_wp])                              ! mean 275
      call output_integrate_fast(files, a, 1_ik, DT)
      call output_integrate_fast(files, b, 1_ik, DT)
      call output_integrate_fast(files, a, 2_ik, DT)
      call close_tier(files, a, 1_ik)
      k_cas = find_var_index(files%reg, 'cas_temp_fast')
      call check_close(a%pending(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'polygon a: its own mean')
      call check(a%queue(1)%n == 1_ik .and. b%queue(1)%n == 0_ik, 'closing polygon a''s window queues only a''s record')
      call check(b%has_data(1) .and. .not. a%has_data(1), 'polygon b''s window stays open')
      call output_integrate_fast(files, b, 2_ik, DT)
      call close_tier(files, b, 1_ik)
      call check_close(b%pending(1)%sval(k_cas), 275.0_wp, 1.0e-10_wp, 'polygon b: its own mean')
      call check_close(a%queue(1)%rec(1)%sval(k_cas), 292.0_wp, 1.0e-10_wp, 'polygon a''s queued record is untouched')
   end subroutine test_two_buffers

   !----- The slow tick across a month boundary (#294). The step from 31 January to 1 February     !
   !      is January's: it folds into January, which then closes. The stand is restructured only   !
   !      after that tick (advance_boundary), so January's cohort slab holds one slot set, and     !
   !      February opens on 1 February with the new one. ------------------------------------------!
   subroutine test_boundary_step()
      type(meds_config_t)    :: cfg
      type(output_files_t)   :: files
      type(output_buffers_t) :: bufs
      type(site_t)           :: site
      type(meds_time_t)      :: jan30, jan31, feb1
      integer(ik) :: k_s, k_c, k_n
      real(wp), parameter :: DT = 86400.0_wp
      cfg = build_test_config(DT)
      cfg%output%enabled    = .true.
      cfg%output%freq_on    = [.false., .true., .true., .false.]    ! daily + monthly
      cfg%output%grp_on     = .false.
      cfg%output%grp_on(1)  = .true.                                 ! structure
      cfg%output%cohort_max = 8_ik
      call manager_alloc(files, bufs, cfg)
      k_s = find_var_index(files%reg, 'agb_site')
      k_c = find_var_index(files%reg, 'agb_cohort')
      k_n = find_var_index(files%reg, 'n_cohort_site')
      call check(k_s > 0_ik .and. k_c > 0_ik .and. k_n > 0_ik, 'boundary fixture registers its variables')

      call site_alloc(site, 2_ik, 8_ik, 8_ik, 4_ik)
      site%patch%n = 1_ik ; site%patch%area(1) = 1.0_wp ; site%patch%cohort_offset(1) = 1_ik
      jan30 = meds_time_t(year=2000_ik, month=1_ik, day=30_ik)
      jan31 = time_advance_days(jan30, 1_ik)
      feb1  = time_advance_days(jan31, 1_ik)

      !----- 30 and 31 January: one cohort, agb 10 then 12 at the ends of the two steps. The      !
      !      second step ends on the month boundary. ---------------------------------------------!
      call set_cohorts(site, [10.0_wp])
      call output_integrate(files, bufs, site, jan30, DT, .true., .false., .false.)
      call set_cohorts(site, [12.0_wp])
      call output_integrate(files, bufs, site, jan31, DT, .true., .true., .false.)

      call check(bufs%pending(2)%t_open%month == 1_ik .and. bufs%pending(2)%t_open%day == 31_ik,     &
                 'the step from 31 January to 1 February is the 31 January daily record')
      call check(bufs%queue(3)%n == 1_ik, 'the month turning closes January')
      call check(bufs%pending(3)%t_open%month == 1_ik .and. bufs%pending(3)%t_open%day == 30_ik,     &
                 'January opens at the start of its first step')
      call check_close(bufs%pending(3)%sval(k_s), 11.0_wp, 1.0e-12_wp,                             &
                       'January''s site mean includes its last step: (10 + 12) / 2')
      call check(bufs%pending(3)%nslab(k_c) == 1_ik, 'January''s cohort slab: its one slot set')
      call check_close(bufs%pending(3)%slab(1,k_c), 11.0_wp, 1.0e-12_wp,                           &
                       'January''s cohort slab includes its last step too')
      call check(.not. bufs%has_data(3), 'nothing of the boundary step is left for February')

      !----- The boundary: a recruit joins (restructure_stand runs here in the driver). Then 1    !
      !      February's own step; close the month by hand. ---------------------------------------!
      call set_cohorts(site, [20.0_wp, 5.0_wp])
      call output_integrate(files, bufs, site, feb1, DT, .true., .false., .false.)
      call check(bufs%t_open(3)%month == 2_ik .and. bufs%t_open(3)%day == 1_ik, 'February opens on 1 February')
      call close_tier(files, bufs, 3_ik)
      call check_close(bufs%pending(3)%sval(k_s), 25.0_wp, 1.0e-12_wp, 'February''s site mean: its own step')
      call check(bufs%pending(3)%nslab(k_c) == 2_ik, 'February''s cohort slab: the new slot set')
      call check_close(bufs%pending(3)%slab(2,k_c), 5.0_wp, 1.0e-12_wp, 'February''s recruit')
      call check_close(bufs%pending(3)%sval(k_n), 2.0_wp, 1.0e-12_wp, 'February''s cohort count')
      call site_free(site)
   end subroutine test_boundary_step

   !----- One patch holding `agb`'s cohorts, one plant each. ------------------------------------!
   subroutine set_cohorts(site, agb)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: agb(:)
      integer(ik) :: n
      n = size(agb, kind=ik)
      site%cohort%n = n ; site%patch%cohort_count(1) = n
      site%cohort%nplant(1:n) = 1.0_wp ; site%cohort%pft(1:n) = 1_ik ; site%cohort%owner_patch(1:n) = 1_ik
      call set_cohort_agb(site, agb, n)
   end subroutine set_cohorts

   !----- Stage two FAST sub-steps with the given CAS temperatures and zero everything else. -----!
   subroutine stage_cas(bufs, cas)
      type(output_buffers_t), intent(inout) :: bufs
      real(wp),            intent(in)    :: cas(2)
      call alloc_fast_staging(bufs, 2_ik, 1_ik, 8_ik)
      bufs%fast_n_soil = 2_ik ; bufs%fast_n_cohort = 1_ik
      bufs%fast_site(PD_CAS_TEMP, :) = cas
      bufs%fast_patch(PD_CAS_TEMP, 1, :) = cas
      bufs%fast_soil_temp = 280.0_wp ; bufs%fast_soil_temp_patch = 280.0_wp
      bufs%fast_coh_ltemp = 285.0_wp ; bufs%fast_coh_height = 10.0_wp
   end subroutine stage_cas

   !----- One polygon's FAST staging as fast_dynamics sizes it (size_fast_staging): nsub sub-steps, !
   !      np patches and nc cohort slots, all zero. ---------------------------------------------!
   subroutine alloc_fast_staging(bufs, nsub, np, nc)
      type(output_buffers_t), intent(inout) :: bufs
      integer(ik),            intent(in)    :: nsub, np, nc
      integer(ik) :: nl
      nl = n_soil_layer_max
      allocate(bufs%fast_time(nsub), bufs%fast_forcing(N_PYDIAG, nsub), bufs%fast_site(N_PDIAG, nsub),  &
               bufs%fast_patch(N_PDIAG, max(np, 1_ik), nsub))
      allocate(bufs%fast_soil_temp(nl, nsub), bufs%fast_soil_water(nl, nsub),                        &
               bufs%fast_soil_temp_patch(nl * max(np, 1_ik), nsub),                                  &
               bufs%fast_soil_water_patch(nl * max(np, 1_ik), nsub))
      allocate(bufs%fast_coh_ltemp(nc, nsub), bufs%fast_coh_gpp(nc, nsub), bufs%fast_coh_height(nc, nsub))
      bufs%fast_forcing = 0.0_wp ; bufs%fast_site = 0.0_wp ; bufs%fast_patch = 0.0_wp
      bufs%fast_soil_temp = 0.0_wp ; bufs%fast_soil_water = 0.0_wp
      bufs%fast_soil_temp_patch = 0.0_wp ; bufs%fast_soil_water_patch = 0.0_wp
      bufs%fast_coh_ltemp = 0.0_wp ; bufs%fast_coh_gpp = 0.0_wp ; bufs%fast_coh_height = 0.0_wp
      bufs%n_fast_sub = nsub ; bufs%fast_n_patch = np
   end subroutine alloc_fast_staging

end program test_output_integrate
