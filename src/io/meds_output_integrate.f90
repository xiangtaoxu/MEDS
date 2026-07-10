!==========================================================================================!
! meds_output_integrate -- the netCDF-FREE temporal-reduction engine (ED2 integrate_/normalize_/  !
! zero_ family): pull the current value of a variable out of live state (extract_variable), fold    !
! it into a per-(variable,tier) running reduction (integrate_scalar / integrate_slab), close a       !
! period (normalize_*), and re-seed for the next (reset_buffer).                                      !
!                                                                                          !
! Links meds_output_types + meds_output_config + the demographic state/reductions ONLY -- no C        !
! bindings, so the stepper edge that calls these pulls no netCDF (the §2 DAG-hygiene wall).            !
!                                                                                          !
! P0 SIMPLIFICATION (honest, and exact for the P0 operators): each ACTIVE tier integrates raw state    !
! independently every step; there is NO inter-tier chaining/feeder map. For AGG_MEAN/TMEAN/SUM/MIN/    !
! MAX/LAST a direct dt-weighted reduction over the coarse period equals the chained roll-up            !
! (§4.1 "chaining exactness"), so the result is identical. Chaining is needed only for AGG_MEANSQ      !
! (variance) and for a sub-dt_fast fast tier feeding daily -- both DEFERRED to P1 (§9). Design §3-§4.    !
!==========================================================================================!
module meds_output_integrate
   use meds_kinds,          only : wp, ik
   use meds_output_types,   only : var_desc_t, integ_buffer_t,                                   &
                                   AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST, AGG_MEANSQ,     &
                                   AGG_TMEAN, AGG_FLUXSUM, DIM_SCALAR, DIM_COHORT, DIM_PATCH,     &
                                   DIM_SOIL, DIM_PFT, MISSING_VALUE
   use meds_demography_types, only : site_t
   use meds_demography_diagnostics, only : total_nplant, total_basal_area, total_agb, total_lai
   implicit none
   private

   public :: alloc_integ_buffer, reset_buffer, integrate_scalar, integrate_slab
   public :: normalize_scalar, normalize_slab, slab_capacity
   public :: extract_variable, extract_scalar_source
   !----- SRC_* : which SoA field / reduction supplies a variable (paired with extract_variable). !
   public :: SRC_C_NPLANT, SRC_C_DBH, SRC_C_HEIGHT, SRC_C_BASAL_AREA, SRC_C_AGB, SRC_C_LEAF_AREA
   public :: SRC_C_GROWTH_AVG, SRC_C_PFT, SRC_C_OWNER_PATCH, SRC_C_GLOBAL_ID
   public :: SRC_P_AREA, SRC_P_AGE, SRC_P_DIST_TYPE, SRC_P_COHORT_OFFSET, SRC_P_COHORT_COUNT
   public :: SRC_P_GLOBAL_ID
   public :: SRC_S_NPLANT, SRC_S_BASAL_AREA, SRC_S_AGB, SRC_S_LAI, SRC_S_N_COHORT, SRC_S_N_PATCH

   !----- Cohort-dimensioned sources (DIM_COHORT). ------------------------------------------!
   integer(ik), parameter :: SRC_C_NPLANT     = 101_ik
   integer(ik), parameter :: SRC_C_DBH        = 102_ik
   integer(ik), parameter :: SRC_C_HEIGHT     = 103_ik
   integer(ik), parameter :: SRC_C_BASAL_AREA = 104_ik
   integer(ik), parameter :: SRC_C_AGB        = 105_ik
   integer(ik), parameter :: SRC_C_LEAF_AREA  = 106_ik
   integer(ik), parameter :: SRC_C_GROWTH_AVG = 107_ik
   integer(ik), parameter :: SRC_C_PFT        = 108_ik
   integer(ik), parameter :: SRC_C_OWNER_PATCH= 109_ik
   integer(ik), parameter :: SRC_C_GLOBAL_ID  = 110_ik
   !----- Patch-dimensioned sources (DIM_PATCH). --------------------------------------------!
   integer(ik), parameter :: SRC_P_AREA          = 201_ik
   integer(ik), parameter :: SRC_P_AGE           = 202_ik
   integer(ik), parameter :: SRC_P_DIST_TYPE     = 203_ik
   integer(ik), parameter :: SRC_P_COHORT_OFFSET = 204_ik
   integer(ik), parameter :: SRC_P_COHORT_COUNT  = 205_ik
   integer(ik), parameter :: SRC_P_GLOBAL_ID     = 206_ik
   !----- Site-scalar reductions (DIM_SCALAR). ----------------------------------------------!
   integer(ik), parameter :: SRC_S_NPLANT      = 301_ik
   integer(ik), parameter :: SRC_S_BASAL_AREA  = 302_ik
   integer(ik), parameter :: SRC_S_AGB         = 303_ik
   integer(ik), parameter :: SRC_S_LAI         = 304_ik
   integer(ik), parameter :: SRC_S_N_COHORT    = 305_ik
   integer(ik), parameter :: SRC_S_N_PATCH     = 306_ik

contains

   !=======================================================================================!
   !  Buffer lifecycle                                                                      !
   !=======================================================================================!

   !----- The re-seed value for an operator: MIN=+huge, MAX=-huge, everything else 0. ------!
   pure real(wp) function seed_of(agg) result(s)
      integer(ik), intent(in) :: agg
      select case (agg)
      case (AGG_MIN) ; s =  huge(1.0_wp)
      case (AGG_MAX) ; s = -huge(1.0_wp)
      case default   ; s =  0.0_wp
      end select
   end function seed_of

   !----- Allocate + seed one (variable, tier) buffer. slab arrays sized to `cap` for a       !
   !      non-scalar dim; scalar dims allocate none.                                            !
   subroutine alloc_integ_buffer(buf, v, freq, cap)
      type(integ_buffer_t), intent(out) :: buf
      type(var_desc_t),     intent(in)  :: v
      integer(ik),          intent(in)  :: freq, cap
      buf%var_id = 0_ik      ! set by the caller (registry index)
      buf%freq   = freq
      buf%agg    = v%agg
      buf%dim    = v%dim
      buf%active = .true.
      buf%seed   = seed_of(v%agg)
      if (v%dim /= DIM_SCALAR) then
         allocate(buf%slab(cap), buf%slab2(cap), buf%wsum_slab(cap), buf%hits(cap))
      end if
      call reset_buffer(buf)
   end subroutine alloc_integ_buffer

   !----- Re-seed a buffer for the next period (ED2 zero_*). --------------------------------!
   subroutine reset_buffer(buf)
      type(integ_buffer_t), intent(inout) :: buf
      buf%scal  = buf%seed
      buf%scal2 = 0.0_wp
      buf%wsum  = 0.0_wp
      buf%nsamp = 0_ik
      buf%n_slab = 0_ik
      if (allocated(buf%slab)) then
         buf%slab(:)      = buf%seed
         buf%slab2(:)     = 0.0_wp
         buf%wsum_slab(:) = 0.0_wp
         buf%hits(:)      = 0_ik
      end if
   end subroutine reset_buffer

   pure integer(ik) function slab_capacity(buf) result(cap)
      type(integ_buffer_t), intent(in) :: buf
      if (allocated(buf%slab)) then ; cap = int(size(buf%slab), ik) ; else ; cap = 0_ik ; end if
   end function slab_capacity

   !=======================================================================================!
   !  Integrate (fold one step)                                                             !
   !=======================================================================================!

   !----- Fold one scalar sample x over a step of length dt [s] into a scalar buffer. ------!
   elemental subroutine integrate_scalar(buf, x, dt)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x, dt
      select case (buf%agg)
      case (AGG_MEAN)   ; buf%scal = buf%scal + x       ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_SUM)    ; buf%scal = buf%scal + x       ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_MIN)    ; buf%scal = min(buf%scal, x)   ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_MAX)    ; buf%scal = max(buf%scal, x)   ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_LAST)   ; buf%scal = x                  ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_TMEAN)  ; buf%scal = buf%scal + x*dt    ; buf%wsum  = buf%wsum + dt
      case (AGG_FLUXSUM); buf%scal = buf%scal + x*dt    ; buf%wsum  = buf%wsum + dt ; buf%nsamp = buf%nsamp + 1_ik
      case (AGG_MEANSQ) ; buf%scal = buf%scal + x*dt    ; buf%scal2 = buf%scal2 + x*x*dt ; buf%wsum = buf%wsum + dt
      end select
   end subroutine integrate_scalar

   !----- Fold a slab of n samples (direct slot index; slot set fixed within a window, §4.4). !
   subroutine integrate_slab(buf, x, n, dt)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x(:)
      integer(ik),          intent(in)    :: n
      real(wp),             intent(in)    :: dt
      integer(ik) :: i
      buf%n_slab = max(buf%n_slab, n)
      select case (buf%agg)
      case (AGG_MEAN, AGG_SUM)
         do i = 1_ik, n
            buf%slab(i) = buf%slab(i) + x(i) ; buf%hits(i) = buf%hits(i) + 1_ik
         end do
      case (AGG_MIN)
         do i = 1_ik, n
            buf%slab(i) = min(buf%slab(i), x(i)) ; buf%hits(i) = buf%hits(i) + 1_ik
         end do
      case (AGG_MAX)
         do i = 1_ik, n
            buf%slab(i) = max(buf%slab(i), x(i)) ; buf%hits(i) = buf%hits(i) + 1_ik
         end do
      case (AGG_LAST)
         do i = 1_ik, n
            buf%slab(i) = x(i) ; buf%hits(i) = buf%hits(i) + 1_ik
         end do
      case (AGG_TMEAN, AGG_FLUXSUM)
         do i = 1_ik, n
            buf%slab(i) = buf%slab(i) + x(i)*dt ; buf%wsum_slab(i) = buf%wsum_slab(i) + dt
            buf%hits(i) = buf%hits(i) + 1_ik
         end do
      case (AGG_MEANSQ)
         do i = 1_ik, n
            buf%slab(i)  = buf%slab(i)  + x(i)*dt
            buf%slab2(i) = buf%slab2(i) + x(i)*x(i)*dt
            buf%wsum_slab(i) = buf%wsum_slab(i) + dt ; buf%hits(i) = buf%hits(i) + 1_ik
         end do
      end select
   end subroutine integrate_slab

   !=======================================================================================!
   !  Normalize (close a period). valid=.false. -> emit MISSING (no NaN, no ±huge leak).     !
   !=======================================================================================!
   subroutine normalize_scalar(buf, out, valid, out2, has2)
      type(integ_buffer_t), intent(in)  :: buf
      real(wp),             intent(out) :: out
      logical,              intent(out) :: valid
      real(wp), optional,   intent(out) :: out2   !< variance companion (AGG_MEANSQ)
      logical,  optional,   intent(out) :: has2
      real(wp) :: mean
      valid = .true. ; out = MISSING_VALUE
      if (present(has2)) has2 = .false.
      select case (buf%agg)
      case (AGG_MEAN)
         if (buf%nsamp > 0_ik) then ; out = buf%scal / real(buf%nsamp, wp) ; else ; valid = .false. ; end if
      case (AGG_TMEAN)
         if (buf%wsum > 0.0_wp) then ; out = buf%scal / buf%wsum ; else ; valid = .false. ; end if
      case (AGG_FLUXSUM)
         if (buf%nsamp > 0_ik) then ; out = buf%scal ; else ; valid = .false. ; end if
      case (AGG_SUM)
         if (buf%nsamp > 0_ik) then ; out = buf%scal ; else ; valid = .false. ; end if
      case (AGG_MIN, AGG_MAX, AGG_LAST)
         if (buf%nsamp > 0_ik) then ; out = buf%scal ; else ; valid = .false. ; end if
      case (AGG_MEANSQ)
         if (buf%wsum > 0.0_wp) then
            mean = buf%scal / buf%wsum ; out = mean
            if (present(out2)) out2 = max(0.0_wp, buf%scal2 / buf%wsum - mean*mean)
            if (present(has2)) has2 = .true.
         else
            valid = .false.
         end if
      end select
   end subroutine normalize_scalar

   !----- Slab normalize: one value + validity per live slot [1:n_out]. --------------------!
   subroutine normalize_slab(buf, out, valid, n_out)
      type(integ_buffer_t), intent(in)  :: buf
      real(wp),             intent(out) :: out(:)
      logical,              intent(out) :: valid(:)
      integer(ik),          intent(out) :: n_out
      integer(ik) :: i
      n_out = buf%n_slab
      out(:)   = MISSING_VALUE
      valid(:) = .false.
      do i = 1_ik, n_out
         select case (buf%agg)
         case (AGG_MEAN)
            if (buf%hits(i) > 0_ik) then ; out(i) = buf%slab(i) / real(buf%hits(i), wp) ; valid(i) = .true. ; end if
         case (AGG_TMEAN)
            if (buf%wsum_slab(i) > 0.0_wp) then ; out(i) = buf%slab(i) / buf%wsum_slab(i) ; valid(i) = .true. ; end if
         case (AGG_FLUXSUM, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST)
            if (buf%hits(i) > 0_ik) then ; out(i) = buf%slab(i) ; valid(i) = .true. ; end if
         case (AGG_MEANSQ)
            if (buf%wsum_slab(i) > 0.0_wp) then ; out(i) = buf%slab(i) / buf%wsum_slab(i) ; valid(i) = .true. ; end if
         end select
      end do
   end subroutine normalize_slab

   !=======================================================================================!
   !  Extract (copy the current value of one variable out of live state; NO pointers, §3.3). !
   !                                                                                        !
   !  nvfortran caveat (CLAUDE.md issue #7 family): nvfortran -O2 does not always see a store  !
   !  to site%<allocatable component> as affecting a subsequent extract_variable(site,...)     !
   !  read, so it can wrongly CSE two extracts across an INLINE mutation. Production is safe --  !
   !  state changes happen through opaque routines (vegetation_dynamics/update_demography) that !
   !  the optimizer cannot see through, and output_integrate extracts each variable ONCE per     !
   !  step. Do NOT extract twice across an inline component store (see test_output_integrate).    !
   !=======================================================================================!
   subroutine extract_variable(site, v, scalar_out, slab_out, n_out)
      type(site_t),     intent(in)  :: site
      type(var_desc_t), intent(in)  :: v
      real(wp),         intent(out) :: scalar_out
      real(wp),         intent(out) :: slab_out(:)
      integer(ik),      intent(out) :: n_out
      integer(ik) :: nc, np
      scalar_out = 0.0_wp ; n_out = 0_ik
      nc = site%cohort%n ; np = site%patch%n
      select case (v%source_id)
      !----- Cohort slabs (DIM_COHORT). -----------------------------------------------------!
      case (SRC_C_NPLANT)     ; slab_out(1:nc) = site%cohort%nplant(1:nc)              ; n_out = nc
      case (SRC_C_DBH)        ; slab_out(1:nc) = site%cohort%dbh(1:nc)                 ; n_out = nc
      case (SRC_C_HEIGHT)     ; slab_out(1:nc) = site%cohort%height(1:nc)              ; n_out = nc
      case (SRC_C_BASAL_AREA) ; slab_out(1:nc) = site%cohort%basal_area(1:nc)          ; n_out = nc
      case (SRC_C_AGB)        ; slab_out(1:nc) = site%cohort%agb(1:nc)                 ; n_out = nc
      case (SRC_C_LEAF_AREA)  ; slab_out(1:nc) = site%cohort%leaf_area(1:nc)           ; n_out = nc
      case (SRC_C_GROWTH_AVG) ; slab_out(1:nc) = site%cohort%growth_avg(1:nc)          ; n_out = nc
      case (SRC_C_PFT)        ; slab_out(1:nc) = real(site%cohort%pft(1:nc), wp)         ; n_out = nc
      case (SRC_C_OWNER_PATCH); slab_out(1:nc) = real(site%cohort%owner_patch(1:nc), wp) ; n_out = nc
      case (SRC_C_GLOBAL_ID)  ; slab_out(1:nc) = real(site%cohort%global_id(1:nc), wp)   ; n_out = nc
      !----- Patch slabs (DIM_PATCH). -------------------------------------------------------!
      case (SRC_P_AREA)          ; slab_out(1:np) = site%patch%area(1:np)                       ; n_out = np
      case (SRC_P_AGE)           ; slab_out(1:np) = site%patch%age(1:np)                        ; n_out = np
      case (SRC_P_DIST_TYPE)     ; slab_out(1:np) = real(site%patch%dist_type(1:np), wp)          ; n_out = np
      case (SRC_P_COHORT_OFFSET) ; slab_out(1:np) = real(site%patch%cohort_offset(1:np), wp)      ; n_out = np
      case (SRC_P_COHORT_COUNT)  ; slab_out(1:np) = real(site%patch%cohort_count(1:np), wp)       ; n_out = np
      case (SRC_P_GLOBAL_ID)     ; slab_out(1:np) = real(site%patch%global_id(1:np), wp)          ; n_out = np
      !----- Site scalars (DIM_SCALAR). -----------------------------------------------------!
      case default
         scalar_out = extract_scalar_source(site, v%source_id)
      end select
   end subroutine extract_variable

   !----- The scalar reductions, split out so a scalar-only caller (and the manager) can use it.-!
   pure real(wp) function extract_scalar_source(site, source_id) result(val)
      type(site_t), intent(in) :: site
      integer(ik),  intent(in) :: source_id
      select case (source_id)
      case (SRC_S_NPLANT)     ; val = total_nplant(site)
      case (SRC_S_BASAL_AREA) ; val = total_basal_area(site)
      case (SRC_S_AGB)        ; val = total_agb(site)
      case (SRC_S_LAI)        ; val = total_lai(site)
      case (SRC_S_N_COHORT)   ; val = real(site%cohort%n, wp)
      case (SRC_S_N_PATCH)    ; val = real(site%patch%n, wp)
      case default            ; val = MISSING_VALUE
      end select
   end function extract_scalar_source

end module meds_output_integrate
