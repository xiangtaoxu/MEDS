!==========================================================================================!
! meds_output_registry -- the ONE source-defined registration list of diagnostic variables, and !
! the config resolution that turns it into per-tier live-variable index lists (§3.2, §3.4, §6.1). !
!                                                                                          !
! netCDF-FREE (links types + config + demography state + the integrate SRC_* codes). Adding a      !
! variable is TWO coupled edits: one add_variable() line here + one case in extract_variable        !
! (meds_output_integrate). Resolution order (§6.1): source defaults -> group toggles -> per-tier      !
! enables -> per-variable overrides (meds_io_config.toml, applied separately) -> annual guard +        !
! freq index. Design: archive/MEDS_IO_DESIGN.md §3, §6.                                               !
!==========================================================================================!
module meds_output_registry
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_output_config,  only : output_config_t, N_FREQ, N_GRP,                               &
                                   FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE
   use meds_output_types,   only : var_desc_t, output_registry_t, MAX_OUTPUT_VARS,               &
                                   AGG_MEAN, AGG_LAST, AGG_TMEAN, DIM_SCALAR, DIM_COHORT,         &
                                   DIM_PATCH, DIM_SOIL, DIM_PFT, XTYPE_DOUBLE, XTYPE_INT
   use meds_output_integrate, only :                                                            &
        SRC_C_NPLANT, SRC_C_DBH, SRC_C_HEIGHT, SRC_C_BASAL_AREA, SRC_C_AGB, SRC_C_LEAF_AREA,      &
        SRC_C_GROWTH_AVG, SRC_C_PFT, SRC_C_OWNER_PATCH, SRC_C_GLOBAL_ID,                          &
        SRC_P_AREA, SRC_P_AGE, SRC_P_DIST_TYPE, SRC_P_COHORT_OFFSET, SRC_P_COHORT_COUNT,          &
        SRC_P_GLOBAL_ID, SRC_S_NPLANT, SRC_S_BASAL_AREA, SRC_S_AGB, SRC_S_LAI,                    &
        SRC_S_N_COHORT, SRC_S_N_PATCH
   implicit none
   private

   public :: build_output_registry, build_freq_index, find_var_index, parse_stream_mask
   public :: apply_group_toggles, apply_freq_enables, apply_variable_override
   public :: freq_bit, OVR_TRUE, OVR_FALSE, OVR_MASK

   !----- Per-variable override kinds (meds_io_config.toml value grammar, §6.1). --------------!
   integer(ik), parameter :: OVR_TRUE  = 1_ik  !< bool true  -> restore registry-default streams
   integer(ik), parameter :: OVR_FALSE = 2_ik  !< bool false -> disable everywhere
   integer(ik), parameter :: OVR_MASK  = 3_ik  !< stream-string -> replace the stream mask

contains

   !----- FREQ_* bit for a 1-based tier index (1=FAST,2=DAILY,3=MONTHLY,4=ANNUAL). ------------!
   pure integer(ik) function freq_bit(tier) result(b)
      integer(ik), intent(in) :: tier
      b = ishft(1_ik, tier - 1_ik)
   end function freq_bit

   !=======================================================================================!
   !  Build the registry: the source list, then apply the main-config group/tier resolution. !
   !=======================================================================================!
   subroutine build_output_registry(reg, cfg)
      type(output_registry_t), intent(out) :: reg
      type(meds_config_t),     intent(in)  :: cfg
      integer(ik), parameter :: DAY_MON    = ior(FREQ_DAILY, FREQ_MONTHLY)
      integer(ik), parameter :: MON        = FREQ_MONTHLY
      integer(ik), parameter :: DAY_MON_YR = ior(ior(FREQ_DAILY, FREQ_MONTHLY), FREQ_ANNUAL)

      allocate(reg%var(MAX_OUTPUT_VARS)) ; reg%nvar = 0_ik
      allocate(reg%idx_freq(MAX_OUTPUT_VARS, N_FREQ)) ; reg%idx_freq = 0_ik ; reg%nidx = 0_ik

      !----- P0 registration list (§3.5, structure group). Args: name, long_name, units, dim, agg,  !
      !      group, default streams, source [, xtype]. Cohort & patch rows carry NO annual bit;       !
      !      GRP_STRUCTURE = 1. (Wrapped to <=132 cols; scale suffixes per §3.5.)                      !
      !--- cohort (DIM_COHORT) ---!
      call add_variable(reg, 'nplant_cohort', 'plant number density', 'plant/m2',                &
                        DIM_COHORT, AGG_LAST, 1_ik, DAY_MON, SRC_C_NPLANT)
      call add_variable(reg, 'dbh_cohort', 'diameter at breast height', 'cm',                    &
                        DIM_COHORT, AGG_MEAN, 1_ik, DAY_MON, SRC_C_DBH)
      call add_variable(reg, 'height_cohort', 'height', 'm',                                     &
                        DIM_COHORT, AGG_MEAN, 1_ik, DAY_MON, SRC_C_HEIGHT)
      call add_variable(reg, 'basal_area_cohort', 'basal area per plant', 'cm2/plant',           &
                        DIM_COHORT, AGG_MEAN, 1_ik, MON, SRC_C_BASAL_AREA)
      call add_variable(reg, 'agb_cohort', 'aboveground biomass per plant', 'kgC/plant',         &
                        DIM_COHORT, AGG_MEAN, 1_ik, DAY_MON, SRC_C_AGB)
      call add_variable(reg, 'leaf_area_cohort', 'leaf area per plant', 'm2/plant',              &
                        DIM_COHORT, AGG_MEAN, 1_ik, MON, SRC_C_LEAF_AREA)
      call add_variable(reg, 'growth_avg_cohort', 'moving-average growth', 'cm/yr',              &
                        DIM_COHORT, AGG_MEAN, 1_ik, MON, SRC_C_GROWTH_AVG)
      call add_variable(reg, 'pft_cohort', 'plant functional type index', '-',                  &
                        DIM_COHORT, AGG_LAST, 1_ik, DAY_MON, SRC_C_PFT, XTYPE_INT)
      call add_variable(reg, 'owner_patch', 'owning patch (1-based)', '-',                       &
                        DIM_COHORT, AGG_LAST, 1_ik, DAY_MON, SRC_C_OWNER_PATCH, XTYPE_INT)
      call add_variable(reg, 'global_cohort_id', 'persistent cohort id', '-',                    &
                        DIM_COHORT, AGG_LAST, 1_ik, DAY_MON, SRC_C_GLOBAL_ID, XTYPE_INT)
      !--- patch (DIM_PATCH) ---!
      call add_variable(reg, 'area_patch', 'patch area fraction', '-',                           &
                        DIM_PATCH, AGG_LAST, 1_ik, DAY_MON, SRC_P_AREA)
      call add_variable(reg, 'age_patch', 'time since disturbance', 'yr',                        &
                        DIM_PATCH, AGG_LAST, 1_ik, MON, SRC_P_AGE)
      call add_variable(reg, 'dist_type_patch', 'disturbance type', '-',                         &
                        DIM_PATCH, AGG_LAST, 1_ik, MON, SRC_P_DIST_TYPE, XTYPE_INT)
      call add_variable(reg, 'cohort_offset', 'first cohort of patch (CSR)', '-',                &
                        DIM_PATCH, AGG_LAST, 1_ik, DAY_MON, SRC_P_COHORT_OFFSET, XTYPE_INT)
      call add_variable(reg, 'cohort_count', 'cohorts in patch (CSR)', '-',                      &
                        DIM_PATCH, AGG_LAST, 1_ik, DAY_MON, SRC_P_COHORT_COUNT, XTYPE_INT)
      call add_variable(reg, 'global_patch_id', 'persistent patch id', '-',                      &
                        DIM_PATCH, AGG_LAST, 1_ik, DAY_MON, SRC_P_GLOBAL_ID, XTYPE_INT)
      !--- site scalars (DIM_SCALAR) ---!
      call add_variable(reg, 'nplant_site', 'site total plant number', 'plant/m2',               &
                        DIM_SCALAR, AGG_MEAN, 1_ik, DAY_MON_YR, SRC_S_NPLANT)
      call add_variable(reg, 'basal_area_site', 'site total basal area', 'm2/m2',                &
                        DIM_SCALAR, AGG_MEAN, 1_ik, DAY_MON_YR, SRC_S_BASAL_AREA)
      call add_variable(reg, 'agb_site', 'site aboveground biomass', 'kgC/m2',                   &
                        DIM_SCALAR, AGG_MEAN, 1_ik, DAY_MON_YR, SRC_S_AGB)
      call add_variable(reg, 'lai_site', 'site leaf area index', 'm2/m2',                        &
                        DIM_SCALAR, AGG_MEAN, 1_ik, DAY_MON_YR, SRC_S_LAI)
      call add_variable(reg, 'n_cohort', 'cohorts in use', '-',                                  &
                        DIM_SCALAR, AGG_LAST, 1_ik, DAY_MON_YR, SRC_S_N_COHORT, XTYPE_INT)
      call add_variable(reg, 'n_patch', 'patches in use', '-',                                   &
                        DIM_SCALAR, AGG_LAST, 1_ik, DAY_MON_YR, SRC_S_N_PATCH, XTYPE_INT)

      call enforce_annual_guard(reg)              ! cohort/patch var MUST NOT declare FREQ_ANNUAL (§3.1)
      call apply_group_toggles(reg, cfg%output)   ! high-level flux-group toggles (§6.1 step 2)
      call apply_freq_enables(reg, cfg%output)    ! per-tier enable (§6.1 step 3)
      call build_freq_index(reg)                  ! precompute idx_freq (per-variable overrides applied later)
   end subroutine build_output_registry

   !----- Append one descriptor (the single "add a variable" registry edit). -----------------!
   subroutine add_variable(reg, nm, ln, un, dm, ag, gr, st, src, xt)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: nm, ln, un
      integer(ik),             intent(in)    :: dm, ag, gr, st, src
      integer(ik), optional,   intent(in)    :: xt
      integer(ik) :: k
      if (reg%nvar >= MAX_OUTPUT_VARS) error stop 'meds_output_registry: MAX_OUTPUT_VARS exceeded'
      reg%nvar = reg%nvar + 1_ik ; k = reg%nvar
      reg%var(k)%name            = nm
      reg%var(k)%long_name       = ln
      reg%var(k)%units           = un
      reg%var(k)%dim             = dm
      reg%var(k)%agg             = ag
      reg%var(k)%group           = gr
      reg%var(k)%streams         = st
      reg%var(k)%streams_default = st
      reg%var(k)%enabled         = .true.
      reg%var(k)%source_id       = src
      reg%var(k)%xtype           = XTYPE_DOUBLE ; if (present(xt)) reg%var(k)%xtype = xt
   end subroutine add_variable

   !----- Guard: a cohort/patch-dimensioned variable may NOT be on the annual stream (§3.1). --!
   subroutine enforce_annual_guard(reg)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: k
      do k = 1_ik, reg%nvar
         if ((reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) .and.               &
             iand(reg%var(k)%streams_default, FREQ_ANNUAL) /= 0_ik)                              &
            error stop 'meds_output_registry: cohort/patch variable on the annual stream ('//trim(reg%var(k)%name)//')'
      end do
   end subroutine enforce_annual_guard

   !----- Step 2: disable every variable in a group whose main-config toggle is .false. -------!
   subroutine apply_group_toggles(reg, out_cfg)
      type(output_registry_t), intent(inout) :: reg
      type(output_config_t),   intent(in)    :: out_cfg
      integer(ik) :: k, g
      do k = 1_ik, reg%nvar
         g = reg%var(k)%group
         if (g >= 1_ik .and. g <= N_GRP) then
            if (.not. out_cfg%grp_on(g)) reg%var(k)%streams = FREQ_NONE
         end if
      end do
   end subroutine apply_group_toggles

   !----- Step 3: a disabled tier clears its FREQ_* bit from every variable's mask. -----------!
   subroutine apply_freq_enables(reg, out_cfg)
      type(output_registry_t), intent(inout) :: reg
      type(output_config_t),   intent(in)    :: out_cfg
      integer(ik) :: k, t, bit
      do t = 1_ik, N_FREQ
         if (.not. out_cfg%freq_on(t)) then
            bit = freq_bit(t)
            do k = 1_ik, reg%nvar
               reg%var(k)%streams = iand(reg%var(k)%streams, not(bit))
            end do
         end if
      end do
   end subroutine apply_freq_enables

   !----- Step 4 primitive: apply one per-variable override (from meds_io_config.toml, §6.1).  !
   !      found=.false. => the name matched no registry variable (the unknown-key trap caller). !
   subroutine apply_variable_override(reg, name, kind, mask, found)
      type(output_registry_t), intent(inout) :: reg
      character(len=*),        intent(in)    :: name
      integer(ik),             intent(in)    :: kind   !< OVR_TRUE | OVR_FALSE | OVR_MASK
      integer(ik),             intent(in)    :: mask   !< used for OVR_MASK
      logical,                 intent(out)   :: found
      integer(ik) :: k, m
      k = find_var_index(reg, name)
      found = (k > 0_ik)
      if (.not. found) return
      select case (kind)
      case (OVR_TRUE)  ; reg%var(k)%streams = reg%var(k)%streams_default
      case (OVR_FALSE) ; reg%var(k)%streams = FREQ_NONE
      case (OVR_MASK)
         m = mask
         if ((reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) .and.               &
             iand(m, FREQ_ANNUAL) /= 0_ik)                                                       &
            error stop 'meds_output_registry: annual (Y) stream forbidden on cohort/patch variable ('//trim(name)//')'
         reg%var(k)%streams = m
      end select
   end subroutine apply_variable_override

   !----- 1-based registry index of a variable by name; 0 if absent. -------------------------!
   pure integer(ik) function find_var_index(reg, name) result(idx)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik) :: k
      idx = 0_ik
      do k = 1_ik, reg%nvar
         if (trim(reg%var(k)%name) == trim(name)) then ; idx = k ; return ; end if
      end do
   end function find_var_index

   !----- Parse a "F D M Y" stream-string into an ior(FREQ_*) mask (§6.1 value grammar). ------!
   !      status = 0 ok; status = 1 unrecognized token (the offender returned in bad_token).    !
   subroutine parse_stream_mask(s, mask, status, bad_token)
      character(len=*), intent(in)  :: s
      integer(ik),      intent(out) :: mask, status
      character(len=*), intent(out) :: bad_token
      integer(ik) :: i
      character(len=1) :: c
      mask = FREQ_NONE ; status = 0_ik ; bad_token = ''
      do i = 1_ik, int(len_trim(s), ik)
         c = s(i:i)
         select case (c)
         case ('F','f') ; mask = ior(mask, FREQ_FAST)
         case ('D','d') ; mask = ior(mask, FREQ_DAILY)
         case ('M','m') ; mask = ior(mask, FREQ_MONTHLY)
         case ('Y','y') ; mask = ior(mask, FREQ_ANNUAL)
         case (' ')     ! separator
         case default   ; status = 1_ik ; bad_token = c ; return
         end select
      end do
   end subroutine parse_stream_mask

   !----- Step 5: precompute, per tier, the list of live (enabled + in-tier) variable indices.-!
   subroutine build_freq_index(reg)
      type(output_registry_t), intent(inout) :: reg
      integer(ik) :: t, k, bit
      reg%idx_freq = 0_ik ; reg%nidx = 0_ik
      do t = 1_ik, N_FREQ
         bit = freq_bit(t)
         do k = 1_ik, reg%nvar
            if (reg%var(k)%enabled .and. iand(reg%var(k)%streams, bit) /= 0_ik) then
               reg%nidx(t) = reg%nidx(t) + 1_ik
               reg%idx_freq(reg%nidx(t), t) = k
            end if
         end do
      end do
   end subroutine build_freq_index

end module meds_output_registry
