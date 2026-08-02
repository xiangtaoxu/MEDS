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
   use meds_time,           only : meds_time_t
   use meds_output_config,  only : N_FREQ
   use meds_output_types,   only : var_desc_t, integ_buffer_t, output_manager_t, fast_sample_t,   &
                                   diag_params_t,                                                  &
                                   AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST, AGG_MEANSQ,     &
                                   AGG_TMEAN, AGG_FLUXSUM, DIM_SCALAR, DIM_COHORT, DIM_PATCH,     &
                                   DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH, MISSING_VALUE
   use meds_core_state_types,   only : site_t
   use meds_column_state_types, only : n_soil_layer_max
   use meds_diagnostic_kernels, only : cohort_lai, cohort_npp_per_plant, soil_wetness,            &
                                       soil_matric_potential, specific_humidity_to_vpd
   use meds_diagnostic_reduce,  only : reduce_cohort_to_site, reduce_cohort_to_patch,             &
                                       reduce_cohort_to_pft, reduce_cohort_to_size,               &
                                       reduce_patch_to_site, reduce_patch_column_to_site,         &
                                       gather_patch_columns, total_et
   implicit none
   private

   public :: alloc_integ_buffer, reset_buffer, integrate_scalar, integrate_slab
   public :: normalize_scalar, normalize_slab, slab_capacity
   public :: extract_variable, extract_scalar_source, extract_fast_scalar
   public :: output_integrate, output_integrate_fast, close_tier
   public :: src_class, SRCK_COHORT, SRCK_PATCH, SRCK_LAYER, SRCK_SITE, SRCK_FAST
   !----- Cohort FIELDS (1000-1999). ---------------------------------------------------------!
   public :: FLD_C_NPLANT, FLD_C_DBH, FLD_C_HEIGHT, FLD_C_BASAL_AREA, FLD_C_AGB, FLD_C_LEAF_AREA
   public :: FLD_C_GROWTH_AVG, FLD_C_PFT, FLD_C_OWNER_PATCH, FLD_C_GLOBAL_ID
   public :: FLD_C_LAI, FLD_C_LEAF_CARBON, FLD_C_FINEROOT_CARBON, FLD_C_WOOD_CARBON
   public :: FLD_C_STORAGE_CARBON, FLD_C_BGB, FLD_C_VEG_CARBON, FLD_C_ONE
   public :: FLD_C_SLA, FLD_C_VCMAX25, FLD_C_RD25, FLD_C_LLSPAN, FLD_C_OVERTOP_LAI
   public :: FLD_C_GPP_ACCUM, FLD_C_NPP_ACCUM, FLD_C_LEAF_RESP, FLD_C_STEM_RESP, FLD_C_ROOT_RESP
   public :: FLD_C_DMAX_PSI_LEAF, FLD_C_PHENO_FLUSH, FLD_C_PHENO_SHED
   public :: FLD_C_LEAF_TEMP, FLD_C_WOOD_TEMP
   public :: FLD_C_DDBH_DT, FLD_C_DAGB_DT, FLD_C_MORT_RATE
   !----- Patch FIELDS (2000-2999). ----------------------------------------------------------!
   public :: FLD_P_AREA, FLD_P_AGE, FLD_P_DIST_TYPE, FLD_P_COHORT_OFFSET, FLD_P_COHORT_COUNT
   public :: FLD_P_GLOBAL_ID, FLD_P_CAS_TEMP, FLD_P_CAS_SHV, FLD_P_CAS_CO2, FLD_P_CAS_VPD
   public :: FLD_P_CAS_DEPTH, FLD_P_SOIL_TEMP_TOP, FLD_P_SWE, FLD_P_SNOW_DEPTH, FLD_P_W_SURFACE
   public :: FLD_P_SOILC_FAST_GRND, FLD_P_SOILC_FAST_SOIL, FLD_P_SOILC_STRUCT_GRND
   public :: FLD_P_SOILC_STRUCT_SOIL, FLD_P_SOILC_MICROBIAL, FLD_P_SOILC_SLOW
   public :: FLD_P_SOILC_PASSIVE, FLD_P_SOILC_TOTAL, FLD_P_RH
   !----- Soil LAYER fields (3000-3999). -----------------------------------------------------!
   public :: FLD_L_SOIL_TEMP, FLD_L_SOIL_WATER, FLD_L_SOIL_PSI, FLD_L_SOIL_WETNESS, FLD_L_SOIL_FLIQ
   !----- Site-only scalars (4000-4999). -----------------------------------------------------!
   public :: SRC_S_ET, SRC_S_N_COHORT, SRC_S_N_PATCH, SRC_S_CANOPY_HEIGHT
   public :: SRC_S_WORK_STEPS, SRC_S_WORK_REJ, SRC_S_WORK_SOIL_NSUB,                              &
             SRC_S_WORK_HYDRO_NSUB, SRC_S_WORK_NONCONV, SRC_S_WORK_HYDRO_THRASH
   public :: SRC_S_WORK_RK45_RESCUE, SRC_S_WORK_CLAMP_STAGE, SRC_S_WORK_CLAMP_COMMIT,             &
             SRC_S_WORK_CLAMP_MASS, SRC_S_WORK_CLAMP_ENERGY
   !----- FAST-tier instantaneous sources (5000-5999): resolved against the live fast_sample_t.  !
   public :: SRC_F_GPP_RATE, SRC_F_LE, SRC_F_H, SRC_F_RNET, SRC_F_SW_IN, SRC_F_USTAR, SRC_F_AIR_TEMP
   public :: SRC_F_CAS_TEMP, SRC_F_SOIL_TEMP_TOP, SRC_F_SOIL_TEMP, SRC_F_SOIL_WATER
   public :: SRC_F_COH_LEAF_TEMP, SRC_F_COH_GPP, SRC_F_COH_HEIGHT

   !==========================================================================================!
   !  SOURCE CODE SPACE. Each source id names a FIELD, and its NUMERIC RANGE says which entity   !
   !  the field lives on. The dispatcher reads the range to decide which reduction applies, so a   !
   !  variable's scale (cohort / patch / site / PFT / size) is a property of its `dim`, not of a    !
   !  separate source id -- which is what lets one field emit every scale from one registry line.   !
   !                                                                                          !
   !  THE RANGES ARE LOAD-BEARING, NOT COSMETIC. The previous flat numbering had SRC_F_GPP_RATE      !
   !  == SRC_S_SOILC_FAST_GRND == 312 (and five more collisions in the same block). Those were       !
   !  invisible only because the two ids were consumed by different switchboards -- but a user who    !
   !  moved a FAST variable onto the daily tier via meds_io_config.toml (`gpp_rate_fast = "F D"`)      !
   !  would have had the daily record silently filled with a soil-carbon pool. Disjoint ranges +       !
   !  the src_class() dispatcher make that class of collision impossible to reintroduce.               !
   !==========================================================================================!
   integer(ik), parameter :: SRCK_NONE   = 0_ik
   integer(ik), parameter :: SRCK_COHORT = 1_ik   !< 1000-1999 : per-cohort field
   integer(ik), parameter :: SRCK_PATCH  = 2_ik   !< 2000-2999 : per-patch field
   integer(ik), parameter :: SRCK_LAYER  = 3_ik   !< 3000-3999 : per-soil-layer field (per patch)
   integer(ik), parameter :: SRCK_SITE   = 4_ik   !< 4000-4999 : site-only scalar (no reduction)
   integer(ik), parameter :: SRCK_FAST   = 5_ik   !< 5000-5999 : FAST-tier staged sample

   !----- Cohort FIELDS (SRCK_COHORT). Raw SoA columns first, then derived ones. --------------!
   integer(ik), parameter :: FLD_C_NPLANT          = 1001_ik
   integer(ik), parameter :: FLD_C_DBH             = 1002_ik
   integer(ik), parameter :: FLD_C_HEIGHT          = 1003_ik
   integer(ik), parameter :: FLD_C_BASAL_AREA      = 1004_ik
   integer(ik), parameter :: FLD_C_AGB             = 1005_ik
   integer(ik), parameter :: FLD_C_LEAF_AREA       = 1006_ik
   integer(ik), parameter :: FLD_C_GROWTH_AVG      = 1007_ik
   integer(ik), parameter :: FLD_C_PFT             = 1008_ik
   integer(ik), parameter :: FLD_C_OWNER_PATCH     = 1009_ik
   integer(ik), parameter :: FLD_C_GLOBAL_ID       = 1010_ik
   integer(ik), parameter :: FLD_C_LEAF_CARBON     = 1011_ik
   integer(ik), parameter :: FLD_C_FINEROOT_CARBON = 1012_ik
   integer(ik), parameter :: FLD_C_WOOD_CARBON     = 1013_ik
   integer(ik), parameter :: FLD_C_STORAGE_CARBON  = 1014_ik
   integer(ik), parameter :: FLD_C_SLA             = 1015_ik
   integer(ik), parameter :: FLD_C_VCMAX25         = 1016_ik
   integer(ik), parameter :: FLD_C_RD25            = 1017_ik
   integer(ik), parameter :: FLD_C_LLSPAN          = 1018_ik
   integer(ik), parameter :: FLD_C_OVERTOP_LAI     = 1019_ik
   integer(ik), parameter :: FLD_C_GPP_ACCUM       = 1020_ik
   integer(ik), parameter :: FLD_C_LEAF_RESP       = 1021_ik
   integer(ik), parameter :: FLD_C_STEM_RESP       = 1022_ik
   integer(ik), parameter :: FLD_C_ROOT_RESP       = 1023_ik
   integer(ik), parameter :: FLD_C_DMAX_PSI_LEAF   = 1024_ik
   integer(ik), parameter :: FLD_C_PHENO_FLUSH     = 1025_ik
   integer(ik), parameter :: FLD_C_PHENO_SHED      = 1026_ik
   integer(ik), parameter :: FLD_C_LEAF_TEMP       = 1027_ik
   integer(ik), parameter :: FLD_C_WOOD_TEMP       = 1028_ik
   !----- DERIVED cohort fields (computed by meds_diagnostic_kernels / from the deriv bundle). --!
   integer(ik), parameter :: FLD_C_ONE             = 1050_ik !< constant 1 (stem counts on any axis)
   integer(ik), parameter :: FLD_C_LAI             = 1051_ik !< nplant*leaf_area [m2/m2]
   integer(ik), parameter :: FLD_C_NPP_ACCUM       = 1052_ik !< gpp - (leaf+stem+root) maintenance resp
   integer(ik), parameter :: FLD_C_BGB             = 1053_ik !< belowground biomass per plant
   integer(ik), parameter :: FLD_C_VEG_CARBON      = 1054_ik !< leaf+fineroot+wood+storage per plant
   !----- Slow-loop TENDENCIES, read straight off site%deriv (still live at the output tick). ---!
   integer(ik), parameter :: FLD_C_DDBH_DT         = 1060_ik !< [cm/yr]
   integer(ik), parameter :: FLD_C_DAGB_DT         = 1061_ik !< [kgC/plant/yr]
   integer(ik), parameter :: FLD_C_MORT_RATE       = 1062_ik !< [1/yr] = -dln_nplant_dt

   !----- Patch FIELDS (SRCK_PATCH). ---------------------------------------------------------!
   integer(ik), parameter :: FLD_P_AREA            = 2001_ik
   integer(ik), parameter :: FLD_P_AGE             = 2002_ik
   integer(ik), parameter :: FLD_P_DIST_TYPE       = 2003_ik
   integer(ik), parameter :: FLD_P_COHORT_OFFSET   = 2004_ik
   integer(ik), parameter :: FLD_P_COHORT_COUNT    = 2005_ik
   integer(ik), parameter :: FLD_P_GLOBAL_ID       = 2006_ik
   integer(ik), parameter :: FLD_P_CAS_TEMP        = 2010_ik
   integer(ik), parameter :: FLD_P_CAS_SHV         = 2011_ik
   integer(ik), parameter :: FLD_P_CAS_CO2         = 2012_ik
   integer(ik), parameter :: FLD_P_CAS_VPD         = 2013_ik
   integer(ik), parameter :: FLD_P_CAS_DEPTH       = 2014_ik
   integer(ik), parameter :: FLD_P_SOIL_TEMP_TOP   = 2015_ik
   integer(ik), parameter :: FLD_P_SWE             = 2016_ik
   integer(ik), parameter :: FLD_P_SNOW_DEPTH      = 2017_ik
   integer(ik), parameter :: FLD_P_W_SURFACE       = 2018_ik
   integer(ik), parameter :: FLD_P_SOILC_FAST_GRND   = 2020_ik
   integer(ik), parameter :: FLD_P_SOILC_FAST_SOIL   = 2021_ik
   integer(ik), parameter :: FLD_P_SOILC_STRUCT_GRND = 2022_ik
   integer(ik), parameter :: FLD_P_SOILC_STRUCT_SOIL = 2023_ik
   integer(ik), parameter :: FLD_P_SOILC_MICROBIAL   = 2024_ik
   integer(ik), parameter :: FLD_P_SOILC_SLOW        = 2025_ik
   integer(ik), parameter :: FLD_P_SOILC_PASSIVE     = 2026_ik
   integer(ik), parameter :: FLD_P_SOILC_TOTAL       = 2027_ik
   integer(ik), parameter :: FLD_P_RH                = 2028_ik

   !----- Soil LAYER fields (SRCK_LAYER): per (layer, patch); reduced to a site column for      !
   !      DIM_SOIL, or flattened for DIM_SOIL_PATCH. ------------------------------------------!
   integer(ik), parameter :: FLD_L_SOIL_TEMP    = 3001_ik
   integer(ik), parameter :: FLD_L_SOIL_WATER   = 3002_ik
   integer(ik), parameter :: FLD_L_SOIL_PSI     = 3003_ik
   integer(ik), parameter :: FLD_L_SOIL_WETNESS = 3004_ik
   integer(ik), parameter :: FLD_L_SOIL_FLIQ    = 3005_ik

   !----- Site-only scalars (SRCK_SITE): no reduction -- these ARE site quantities already. ----!
   integer(ik), parameter :: SRC_S_ET               = 4001_ik !< [kg/m2] slow-step ET accumulator
   integer(ik), parameter :: SRC_S_N_COHORT         = 4002_ik
   integer(ik), parameter :: SRC_S_N_PATCH          = 4003_ik
   integer(ik), parameter :: SRC_S_CANOPY_HEIGHT    = 4004_ik !< [m] tallest live cohort
   integer(ik), parameter :: SRC_S_WORK_STEPS       = 4030_ik
   integer(ik), parameter :: SRC_S_WORK_REJ         = 4031_ik
   integer(ik), parameter :: SRC_S_WORK_SOIL_NSUB   = 4032_ik
   integer(ik), parameter :: SRC_S_WORK_HYDRO_NSUB  = 4033_ik
   integer(ik), parameter :: SRC_S_WORK_NONCONV     = 4034_ik
   integer(ik), parameter :: SRC_S_WORK_RK45_RESCUE = 4035_ik
   integer(ik), parameter :: SRC_S_WORK_CLAMP_STAGE = 4036_ik
   integer(ik), parameter :: SRC_S_WORK_CLAMP_COMMIT= 4037_ik
   integer(ik), parameter :: SRC_S_WORK_CLAMP_MASS  = 4038_ik
   integer(ik), parameter :: SRC_S_WORK_CLAMP_ENERGY= 4039_ik
   integer(ik), parameter :: SRC_S_WORK_HYDRO_THRASH= 4040_ik

   !----- FAST-tier staged sources (SRCK_FAST): read from the live fast_sample_t / manager       !
   !      slabs, never from site state (which at replay time is the end-of-slow-step snapshot).   !
   integer(ik), parameter :: SRC_F_CAS_TEMP      = 5001_ik
   integer(ik), parameter :: SRC_F_SOIL_TEMP_TOP = 5002_ik
   integer(ik), parameter :: SRC_F_GPP_RATE      = 5003_ik
   integer(ik), parameter :: SRC_F_LE            = 5004_ik
   integer(ik), parameter :: SRC_F_H             = 5005_ik
   integer(ik), parameter :: SRC_F_RNET          = 5006_ik
   integer(ik), parameter :: SRC_F_SW_IN         = 5007_ik
   integer(ik), parameter :: SRC_F_USTAR         = 5008_ik
   integer(ik), parameter :: SRC_F_AIR_TEMP      = 5009_ik
   integer(ik), parameter :: SRC_F_SOIL_TEMP     = 5010_ik  !< DIM_SOIL, from the fast soil slab
   integer(ik), parameter :: SRC_F_SOIL_WATER    = 5011_ik  !< DIM_SOIL, from the fast soil slab
   integer(ik), parameter :: SRC_F_COH_LEAF_TEMP = 5020_ik  !< DIM_COHORT
   integer(ik), parameter :: SRC_F_COH_GPP       = 5021_ik  !< DIM_COHORT
   integer(ik), parameter :: SRC_F_COH_HEIGHT    = 5022_ik  !< DIM_COHORT

   !----- Reference pressure [Pa] for the DIAGNOSTIC canopy-air VPD read-off. The CAS box       !
   !      carries no prognostic pressure, so a VPD from its two twins needs one supplied. Using   !
   !      the standard atmosphere makes this a diagnostic-grade signal (right shape, right         !
   !      magnitude for canopy coupling) rather than a thermodynamic state variable, and that      !
   !      limitation is stated here rather than left for a reader to discover from the numbers.    !
   real(wp), parameter :: PRSS_REF = 101325.0_wp

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

   !----- Fold a slab of n samples (direct slot index; slot set fixed within a window, §4.4).  !
   !                                                                                          !
   !      `valid` (optional) marks slots that had NOTHING to report this step -- an empty patch, !
   !      a PFT with no members, a size class with no stems, a soil column with no texture. Those !
   !      slots are SKIPPED, so their hits/wsum stay 0 and normalize emits _FillValue. Folding a   !
   !      0 into them instead would be the classic diagnostics lie: an empty bin reported as a      !
   !      real zero, which then drags a period mean toward 0 and is indistinguishable in the file    !
   !      from a genuine measurement of zero. Absent `valid` => every slot in [1:n] is folded.       !
   subroutine integrate_slab(buf, x, n, dt, valid)
      type(integ_buffer_t), intent(inout) :: buf
      real(wp),             intent(in)    :: x(:)
      integer(ik),          intent(in)    :: n
      real(wp),             intent(in)    :: dt
      logical, optional,    intent(in)    :: valid(:)
      integer(ik) :: i
      logical     :: use_v
      buf%n_slab = max(buf%n_slab, n)
      use_v = present(valid)
      do i = 1_ik, n
         if (use_v) then
            if (.not. valid(i)) cycle
         end if
         select case (buf%agg)
         case (AGG_MEAN, AGG_SUM)
            buf%slab(i) = buf%slab(i) + x(i) ; buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_MIN)
            buf%slab(i) = min(buf%slab(i), x(i)) ; buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_MAX)
            buf%slab(i) = max(buf%slab(i), x(i)) ; buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_LAST)
            buf%slab(i) = x(i) ; buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_TMEAN, AGG_FLUXSUM)
            buf%slab(i) = buf%slab(i) + x(i)*dt ; buf%wsum_slab(i) = buf%wsum_slab(i) + dt
            buf%hits(i) = buf%hits(i) + 1_ik
         case (AGG_MEANSQ)
            buf%slab(i)  = buf%slab(i)  + x(i)*dt
            buf%slab2(i) = buf%slab2(i) + x(i)*x(i)*dt
            buf%wsum_slab(i) = buf%wsum_slab(i) + dt ; buf%hits(i) = buf%hits(i) + 1_ik
         end select
      end do
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
   !----- Which entity a source id names, from its numeric range (see the SOURCE CODE SPACE     !
   !      block at the top). This one function is what makes the id ranges enforceable rather      !
   !      than a convention -- every dispatch below asks it, so a mis-ranged id fails loudly here   !
   !      instead of silently resolving to a neighbour's field.                                     !
   pure integer(ik) function src_class(src) result(k)
      integer(ik), intent(in) :: src
      if      (src >= 1000_ik .and. src < 2000_ik) then ; k = SRCK_COHORT
      else if (src >= 2000_ik .and. src < 3000_ik) then ; k = SRCK_PATCH
      else if (src >= 3000_ik .and. src < 4000_ik) then ; k = SRCK_LAYER
      else if (src >= 4000_ik .and. src < 5000_ik) then ; k = SRCK_SITE
      else if (src >= 5000_ik .and. src < 6000_ik) then ; k = SRCK_FAST
      else                                              ; k = SRCK_NONE
      end if
   end function src_class

   !=======================================================================================!
   !  FIELD ACCESSORS -- copy one entity's raw column out of live state (stage [1] + the       !
   !  derived kernels). NO reduction here; the reduction is chosen by the descriptor's dim.    !
   !=======================================================================================!

   !----- Per-cohort field over the whole site SoA [1:site%cohort%n]. -----------------------!
   pure subroutine cohort_source_field(site, src, x, n)
      type(site_t), intent(in)  :: site
      integer(ik),  intent(in)  :: src
      real(wp),     intent(out) :: x(:)
      integer(ik),  intent(out) :: n
      n = site%cohort%n
      if (n <= 0_ik) return
      select case (src)
      case (FLD_C_NPLANT)          ; x(1:n) = site%cohort%nplant(1:n)
      case (FLD_C_DBH)             ; x(1:n) = site%cohort%dbh(1:n)
      case (FLD_C_HEIGHT)          ; x(1:n) = site%cohort%height(1:n)
      case (FLD_C_BASAL_AREA)      ; x(1:n) = site%cohort%basal_area(1:n)
      case (FLD_C_AGB)             ; x(1:n) = site%cohort%agb(1:n)
      case (FLD_C_LEAF_AREA)       ; x(1:n) = site%cohort%leaf_area(1:n)
      case (FLD_C_GROWTH_AVG)      ; x(1:n) = site%cohort%growth_avg(1:n)
      case (FLD_C_PFT)             ; x(1:n) = real(site%cohort%pft(1:n), wp)
      case (FLD_C_OWNER_PATCH)     ; x(1:n) = real(site%cohort%owner_patch(1:n), wp)
      case (FLD_C_GLOBAL_ID)       ; x(1:n) = real(site%cohort%global_id(1:n), wp)
      case (FLD_C_LEAF_CARBON)     ; x(1:n) = site%cohort%leaf_carbon(1:n)
      case (FLD_C_FINEROOT_CARBON) ; x(1:n) = site%cohort%fineroot_carbon(1:n)
      case (FLD_C_WOOD_CARBON)     ; x(1:n) = site%cohort%wood_carbon(1:n)
      case (FLD_C_STORAGE_CARBON)  ; x(1:n) = site%cohort%nonstructural_carbon(1:n)
      case (FLD_C_SLA)             ; x(1:n) = site%cohort%sla(1:n)
      case (FLD_C_VCMAX25)         ; x(1:n) = site%cohort%vcmax25(1:n)
      case (FLD_C_RD25)            ; x(1:n) = site%cohort%rd25(1:n)
      case (FLD_C_LLSPAN)          ; x(1:n) = site%cohort%llspan(1:n)
      case (FLD_C_OVERTOP_LAI)     ; x(1:n) = site%cohort%overtopping_lai(1:n)
      case (FLD_C_GPP_ACCUM)       ; x(1:n) = site%cohort%gpp_accum(1:n)
      case (FLD_C_LEAF_RESP)       ; x(1:n) = site%cohort%leaf_resp_accum(1:n)
      case (FLD_C_STEM_RESP)       ; x(1:n) = site%cohort%stem_resp_accum(1:n)
      case (FLD_C_ROOT_RESP)       ; x(1:n) = site%cohort%root_resp_accum(1:n)
      case (FLD_C_DMAX_PSI_LEAF)   ; x(1:n) = site%cohort%dmax_psi_leaf(1:n)
      case (FLD_C_PHENO_FLUSH)     ; x(1:n) = site%cohort%pheno_flush_drive(1:n)
      case (FLD_C_PHENO_SHED)      ; x(1:n) = site%cohort%pheno_shed_drive(1:n)
      case (FLD_C_LEAF_TEMP)       ; x(1:n) = site%cohort%leaf_temp(1:n)
      case (FLD_C_WOOD_TEMP)       ; x(1:n) = site%cohort%wood_temp(1:n)
      !----- derived ------------------------------------------------------------------------!
      case (FLD_C_ONE)             ; x(1:n) = 1.0_wp
      case (FLD_C_LAI)             ; x(1:n) = cohort_lai(site%cohort%nplant(1:n),               &
                                                         site%cohort%leaf_area(1:n))
      case (FLD_C_NPP_ACCUM)       ; x(1:n) = cohort_npp_per_plant(site%cohort%gpp_accum(1:n),  &
                                                site%cohort%leaf_resp_accum(1:n),              &
                                                site%cohort%stem_resp_accum(1:n),              &
                                                site%cohort%root_resp_accum(1:n))
      case (FLD_C_BGB)             ; x(1:n) = site%cohort%fineroot_carbon(1:n)                  &
                                    + site%cohort%wood_carbon(1:n)                             &
                                      * (1.0_wp - site%cohort%p_aboveground_frac(1:n))
      case (FLD_C_VEG_CARBON)      ; x(1:n) = site%cohort%leaf_carbon(1:n)                      &
                                    + site%cohort%fineroot_carbon(1:n)                         &
                                    + site%cohort%wood_carbon(1:n)                             &
                                    + site%cohort%nonstructural_carbon(1:n)
      !----- slow-loop tendencies, read off the still-live per-step deriv bundle. ------------!
      !      site%deriv is refilled every slow step and is NOT lockstep-reordered, so it is      !
      !      meaningful only at the output tick, which is exactly when this runs. A run whose     !
      !      slow tier is frozen ([run].slow_on = false) leaves it at 0, which is the truth.      !
      case (FLD_C_DDBH_DT)
         if (allocated(site%deriv%d_dbh_dt)) then
            x(1:n) = site%deriv%d_dbh_dt(1:n)
         else ; x(1:n) = 0.0_wp ; end if
      case (FLD_C_DAGB_DT)
         if (allocated(site%deriv%d_agb_dt)) then
            x(1:n) = site%deriv%d_agb_dt(1:n)
         else ; x(1:n) = 0.0_wp ; end if
      case (FLD_C_MORT_RATE)
         !----- Reported as a POSITIVE hazard [1/yr]: the engine stores the log-space density    !
         !      tendency dln(nplant)/dt, which is <= 0 under mortality, so the diagnostic is its   !
         !      negation. Sign errors here are the classic way a mortality plot comes out upside   !
         !      down, hence the explicit note.                                                     !
         if (allocated(site%deriv%dln_nplant_dt)) then
            x(1:n) = -site%deriv%dln_nplant_dt(1:n)
         else ; x(1:n) = 0.0_wp ; end if
      case default                 ; x(1:n) = MISSING_VALUE
      end select
   end subroutine cohort_source_field

   !----- Per-patch field [1:site%patch%n]. -------------------------------------------------!
   pure subroutine patch_source_field(site, src, x, n)
      type(site_t), intent(in)  :: site
      integer(ik),  intent(in)  :: src
      real(wp),     intent(out) :: x(:)
      integer(ik),  intent(out) :: n
      integer(ik) :: ip
      n = site%patch%n
      if (n <= 0_ik) return
      select case (src)
      case (FLD_P_AREA)          ; x(1:n) = site%patch%area(1:n)
      case (FLD_P_AGE)           ; x(1:n) = site%patch%age(1:n)
      case (FLD_P_DIST_TYPE)     ; x(1:n) = real(site%patch%dist_type(1:n), wp)
      case (FLD_P_COHORT_OFFSET) ; x(1:n) = real(site%patch%cohort_offset(1:n), wp)
      case (FLD_P_COHORT_COUNT)  ; x(1:n) = real(site%patch%cohort_count(1:n), wp)
      case (FLD_P_GLOBAL_ID)     ; x(1:n) = real(site%patch%global_id(1:n), wp)
      case (FLD_P_CAS_TEMP)      ; do ip = 1_ik, n ; x(ip) = site%patch%cas(ip)%can_temp     ; end do
      case (FLD_P_CAS_SHV)       ; do ip = 1_ik, n ; x(ip) = site%patch%cas(ip)%can_shv      ; end do
      case (FLD_P_CAS_CO2)       ; do ip = 1_ik, n ; x(ip) = site%patch%cas(ip)%can_co2      ; end do
      case (FLD_P_CAS_DEPTH)     ; do ip = 1_ik, n ; x(ip) = site%patch%cas(ip)%can_depth    ; end do
      case (FLD_P_SOIL_TEMP_TOP) ; do ip = 1_ik, n ; x(ip) = site%patch%soil_e(ip)%soil_temp(1) ; end do
      case (FLD_P_W_SURFACE)     ; do ip = 1_ik, n ; x(ip) = site%patch%soil_w(ip)%w_surface ; end do
      case (FLD_P_SWE)
         do ip = 1_ik, n ; x(ip) = sum(site%patch%snow(ip)%swe(:))        ; end do
      case (FLD_P_SNOW_DEPTH)
         do ip = 1_ik, n ; x(ip) = sum(site%patch%snow(ip)%snow_depth(:)) ; end do
      !----- CAS vapour-pressure deficit: a DERIVED read-off of the two prognostic twins.      !
      !      Pressure is the standard-atmosphere reference (the CAS box carries no prognostic   !
      !      pressure), so this is a diagnostic-grade VPD, adequate for the canopy-coupling       !
      !      signal it exists to show and not to be mistaken for a thermodynamic state variable.  !
      case (FLD_P_CAS_VPD)
         do ip = 1_ik, n
            x(ip) = specific_humidity_to_vpd(site%patch%cas(ip)%can_temp,                      &
                                             site%patch%cas(ip)%can_shv, PRSS_REF)
         end do
      case (FLD_P_SOILC_FAST_GRND)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%fast_grnd_carbon   ; end do
      case (FLD_P_SOILC_FAST_SOIL)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%fast_soil_carbon   ; end do
      case (FLD_P_SOILC_STRUCT_GRND)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%struct_grnd_carbon ; end do
      case (FLD_P_SOILC_STRUCT_SOIL)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%struct_soil_carbon ; end do
      case (FLD_P_SOILC_MICROBIAL)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%microbial_carbon   ; end do
      case (FLD_P_SOILC_SLOW)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%slow_carbon        ; end do
      case (FLD_P_SOILC_PASSIVE)
         do ip = 1_ik, n ; x(ip) = site%patch%soil_carbon(ip)%passive_carbon     ; end do
      case (FLD_P_SOILC_TOTAL)
         do ip = 1_ik, n
            x(ip) = site%patch%soil_carbon(ip)%fast_grnd_carbon                                &
                  + site%patch%soil_carbon(ip)%fast_soil_carbon                                &
                  + site%patch%soil_carbon(ip)%struct_grnd_carbon                              &
                  + site%patch%soil_carbon(ip)%struct_soil_carbon                              &
                  + site%patch%soil_carbon(ip)%microbial_carbon                                &
                  + site%patch%soil_carbon(ip)%slow_carbon                                     &
                  + site%patch%soil_carbon(ip)%passive_carbon
         end do
      case (FLD_P_RH)
         do ip = 1_ik, n ; x(ip) = site%patch%xi_accum(ip)%rh_fast_accum ; end do
      case default               ; x(1:n) = MISSING_VALUE
      end select
   end subroutine patch_source_field

   !----- Per-(soil layer, patch) field. `col` is (n_soil_layer_max, n_patch); nlayer is the   !
   !      ACTIVE layer count. psi/wetness need the retention curve, so they emit MISSING when    !
   !      the manager has no soil parameters wired (dp%soil_ready = .false.) rather than         !
   !      inventing a texture -- a plausible wrong psi is worse than an honest fill value.        !
   pure subroutine layer_source_field(site, dp, src, col, nlayer, ok)
      type(site_t),        intent(in)  :: site
      type(diag_params_t), intent(in)  :: dp
      integer(ik),         intent(in)  :: src
      real(wp),            intent(out) :: col(:,:)
      integer(ik),         intent(out) :: nlayer
      logical,             intent(out) :: ok
      integer(ik) :: ip, np
      np = site%patch%n
      nlayer = n_soil_layer_max ; ok = .true.
      if (np <= 0_ik) then ; ok = .false. ; return ; end if
      select case (src)
      case (FLD_L_SOIL_TEMP)
         do ip = 1_ik, np
            col(1:nlayer, ip) = site%patch%soil_e(ip)%soil_temp(1:nlayer)
         end do
      case (FLD_L_SOIL_WATER)
         do ip = 1_ik, np
            col(1:nlayer, ip) = site%patch%soil_w(ip)%theta(1:nlayer)
         end do
      case (FLD_L_SOIL_FLIQ)
         do ip = 1_ik, np
            col(1:nlayer, ip) = site%patch%soil_e(ip)%soil_fliq(1:nlayer)
         end do
      case (FLD_L_SOIL_PSI)
         if (.not. dp%soil_ready) then ; ok = .false. ; return ; end if
         do ip = 1_ik, np
            col(1:nlayer, ip) = soil_matric_potential(dp%retention,                            &
                                   site%patch%soil_w(ip)%theta(1:nlayer),                      &
                                   dp%theta_sat(1:nlayer), dp%theta_res(1:nlayer),             &
                                   dp%par_a(1:nlayer), dp%par_n(1:nlayer))
         end do
      case (FLD_L_SOIL_WETNESS)
         if (.not. dp%soil_ready) then ; ok = .false. ; return ; end if
         do ip = 1_ik, np
            col(1:nlayer, ip) = soil_wetness(site%patch%soil_w(ip)%theta(1:nlayer),            &
                                             dp%theta_res(1:nlayer), dp%theta_sat(1:nlayer))
         end do
      case default ; ok = .false.
      end select
   end subroutine layer_source_field

   !=======================================================================================!
   !  THE DISPATCHER. One variable -> one value (scalar) or one slab, by reducing its FIELD    !
   !  to the scale its `dim` names, with the weight/mean/scale its descriptor declares.         !
   !                                                                                          !
   !  This is the routine that makes "one registry line emits every scale" true: the same        !
   !  FLD_C_AGB serves agb_cohort (raw), agb_patch, agb_site, agb_pft and agb_size, and the       !
   !  ONLY difference between those five registry rows is dim + weight + mean.                    !
   !                                                                                          !
   !  nvfortran caveat (CLAUDE.md issue #7 family): nvfortran -O2 does not always see a store     !
   !  to site%<allocatable component> as affecting a subsequent extract_variable(site,...) read,   !
   !  so it can wrongly CSE two extracts across an INLINE mutation. Production is safe -- state     !
   !  changes happen through opaque routines the optimizer cannot see through, and output_integrate !
   !  extracts each variable ONCE per step. Do NOT extract twice across an inline component store.   !
   !=======================================================================================!
   subroutine extract_variable(site, dp, v, scalar_out, slab_out, valid_out, n_out)
      type(site_t),        intent(in)  :: site
      type(diag_params_t), intent(in)  :: dp
      type(var_desc_t),    intent(in)  :: v
      real(wp),            intent(out) :: scalar_out
      real(wp),            intent(out) :: slab_out(:)
      logical,             intent(out) :: valid_out(:)
      integer(ik),         intent(out) :: n_out
      real(wp)    :: x(max(site%cohort%n, site%patch%n, 1_ik))
      real(wp)    :: col(n_soil_layer_max, max(site%patch%n, 1_ik))
      integer(ik) :: nf, nlayer, kcls
      logical     :: ok
      scalar_out = 0.0_wp ; n_out = 0_ik
      kcls = src_class(v%source_id)

      select case (v%dim)

      !----- RAW cohort slab: the field itself, no reduction. -------------------------------!
      case (DIM_COHORT)
         call cohort_source_field(site, v%source_id, slab_out, n_out)
         valid_out(1:max(n_out,1_ik)) = .true.

      !----- Patch axis: either a per-patch field verbatim, or a cohort field reduced to it. -!
      case (DIM_PATCH)
         if (kcls == SRCK_COHORT) then
            call cohort_source_field(site, v%source_id, x, nf)
            call reduce_cohort_to_patch(site, x, v%weight, v%mean, slab_out, valid_out, n_out, &
                                        scale=v%scale)
         else
            call patch_source_field(site, v%source_id, slab_out, n_out)
            valid_out(1:max(n_out,1_ik)) = .true.
         end if

      !----- PFT axis (run-time length; empty PFTs are a true 0 for sums, see the reducer). --!
      case (DIM_PFT)
         call cohort_source_field(site, v%source_id, x, nf)
         call reduce_cohort_to_pft(site, x, v%weight, v%mean, dp%n_pft, slab_out, valid_out,   &
                                   n_out, scale=v%scale)

      !----- DBH size class (ED2-style: bin by the cohort's mean dbh, never split a cohort). -!
      case (DIM_SIZE)
         call cohort_source_field(site, v%source_id, x, nf)
         call reduce_cohort_to_size(site, x, v%weight, v%mean, dp%dbh_edges, dp%n_dbh_class,   &
                                    slab_out, valid_out, n_out, scale=v%scale)

      !----- Area-weighted SITE soil column. ------------------------------------------------!
      case (DIM_SOIL)
         call layer_source_field(site, dp, v%source_id, col, nlayer, ok)
         if (ok) then
            call reduce_patch_column_to_site(site, col, nlayer, slab_out, n_out)
            valid_out(1:max(n_out,1_ik)) = .true.
         else
            n_out = 0_ik
         end if

      !----- 2-D (soil layer x patch), flattened at (ip-1)*n_soil_layer_max + k. -------------!
      case (DIM_SOIL_PATCH)
         call layer_source_field(site, dp, v%source_id, col, nlayer, ok)
         if (ok) then
            call gather_patch_columns(site, col, nlayer, slab_out, valid_out, n_out)
         else
            n_out = 0_ik
         end if

      !----- Site scalar: a site-only quantity, or a cohort/patch field reduced to the site. -!
      case default
         select case (kcls)
         case (SRCK_COHORT)
            call cohort_source_field(site, v%source_id, x, nf)
            call reduce_cohort_to_site(site, x, v%weight, v%mean, scalar_out, ok, scale=v%scale)
         case (SRCK_PATCH)
            call patch_source_field(site, v%source_id, x, nf)
            call reduce_patch_to_site(site, x, v%mean, scalar_out, ok)
            scalar_out = scalar_out * v%scale
         case default
            scalar_out = extract_scalar_source(site, v%source_id)
         end select
      end select
   end subroutine extract_variable

   !----- Site-ONLY scalars: quantities that are already site-level and need no reduction. ---!
   pure real(wp) function extract_scalar_source(site, source_id) result(val)
      type(site_t), intent(in) :: site
      integer(ik),  intent(in) :: source_id
      select case (source_id)
      case (SRC_S_ET)         ; val = total_et(site)
      case (SRC_S_N_COHORT)   ; val = real(site%cohort%n, wp)
      case (SRC_S_N_PATCH)    ; val = real(site%patch%n, wp)
      case (SRC_S_CANOPY_HEIGHT)
         if (site%cohort%n > 0_ik) then
            val = maxval(site%cohort%height(1:site%cohort%n))
         else
            val = 0.0_wp
         end if
      case (SRC_S_WORK_STEPS)        ; val = site%work_integ_steps
      case (SRC_S_WORK_REJ)          ; val = site%work_integ_rej
      case (SRC_S_WORK_SOIL_NSUB)    ; val = site%work_soil_nsub
      case (SRC_S_WORK_HYDRO_NSUB)   ; val = site%work_hydro_nsub
      case (SRC_S_WORK_NONCONV)      ; val = site%work_nonconv
      case (SRC_S_WORK_HYDRO_THRASH) ; val = site%work_hydro_thrash
      case (SRC_S_WORK_RK45_RESCUE)  ; val = site%work_rk45_rescue
      case (SRC_S_WORK_CLAMP_STAGE)  ; val = site%work_clamp_stage
      case (SRC_S_WORK_CLAMP_COMMIT) ; val = site%work_clamp_commit
      case (SRC_S_WORK_CLAMP_MASS)   ; val = site%work_clamp_mass
      case (SRC_S_WORK_CLAMP_ENERGY) ; val = site%work_clamp_energy
      case default            ; val = MISSING_VALUE
      end select
   end function extract_scalar_source

   !=======================================================================================!
   !  The per-step tick (netCDF-FREE): close any period that ended, then fold the current    !
   !  step into every active tier. Called by the stepper (aux); stages closed records into    !
   !  mgr%pending for main to serialize. FAST tier (index 1) is DEFERRED to P1 (§9); P0 ticks   !
   !  DAILY/MONTHLY/ANNUAL at the slow step. The MONTHLY cohort/patch flush must run BEFORE that  !
   !  month boundary's fiss/fuse (§4.4/§4.5) -- the caller orders this by where it calls the tick. !
   !=======================================================================================!
   subroutine output_integrate(mgr, site, now, dt, is_new_day, is_new_month, is_new_year)
      type(output_manager_t), intent(inout) :: mgr
      type(site_t),           intent(in)    :: site
      type(meds_time_t),      intent(in)    :: now
      real(wp),               intent(in)    :: dt
      logical,                intent(in)    :: is_new_day, is_new_month, is_new_year
      integer(ik) :: t, j, k, n_out
      real(wp)    :: scal
      real(wp)    :: slab(max(mgr%max_slab, 1_ik))
      logical     :: vslab(max(mgr%max_slab, 1_ik))
      if (.not. mgr%enabled) return
      !----- close periods that just ended (before folding the new step; tiers independent). ---!
      if (is_new_year  .and. mgr%has_data(4_ik)) call close_tier(mgr, 4_ik)
      if (is_new_month .and. mgr%has_data(3_ik)) call close_tier(mgr, 3_ik)
      if (is_new_day   .and. mgr%has_data(2_ik)) call close_tier(mgr, 2_ik)
      !----- fold the current step into DAILY/MONTHLY/ANNUAL. The FAST tier is fed separately    !
      !      from the staged sub-step samples (output_integrate_fast), because sub-daily          !
      !      resolution exists only inside the fast loop.                                         !
      do t = 2_ik, N_FREQ
         if (mgr%reg%nidx(t) == 0_ik) cycle
         if (.not. mgr%has_data(t)) mgr%t_open(t) = now
         do j = 1_ik, mgr%reg%nidx(t)
            k = mgr%reg%idx_freq(j, t)
            call extract_variable(site, mgr%diag, mgr%reg%var(k), scal, slab, vslab, n_out)
            if (mgr%reg%var(k)%dim == DIM_SCALAR) then
               call integrate_scalar(mgr%buf(k,t), scal, dt)
            else
               call integrate_slab(mgr%buf(k,t), slab, n_out, dt, vslab)
            end if
         end do
         mgr%has_data(t) = .true.
      end do
   end subroutine output_integrate

   !----- Normalize a tier's buffers into its pending record + reset them (staging, §4.5). -----!
   subroutine close_tier(mgr, t)
      type(output_manager_t), intent(inout) :: mgr
      integer(ik),            intent(in)    :: t
      integer(ik) :: j, k, ns
      mgr%pending(t)%used   = .true.
      mgr%pending(t)%freq   = ishft(1_ik, t - 1_ik)
      mgr%pending(t)%t_open = mgr%t_open(t)
      mgr%pending(t)%n_cohort = 0_ik ; mgr%pending(t)%n_patch = 0_ik
      mgr%pending(t)%sval(:)   = MISSING_VALUE
      mgr%pending(t)%svalid(:) = .false.
      mgr%pending(t)%nslab(:)  = 0_ik
      do j = 1_ik, mgr%reg%nidx(t)
         k = mgr%reg%idx_freq(j, t)
         if (mgr%reg%var(k)%dim == DIM_SCALAR) then
            call normalize_scalar(mgr%buf(k,t), mgr%pending(t)%sval(k), mgr%pending(t)%svalid(k))
         else
            call normalize_slab(mgr%buf(k,t), mgr%pending(t)%slab(:,k),                          &
                                mgr%pending(t)%slabvalid(:,k), ns)
            mgr%pending(t)%nslab(k) = ns
            if (mgr%reg%var(k)%dim == DIM_COHORT)                                                &
               mgr%pending(t)%n_cohort = max(mgr%pending(t)%n_cohort, ns)
            if (mgr%reg%var(k)%dim == DIM_PATCH)                                                 &
               mgr%pending(t)%n_patch  = max(mgr%pending(t)%n_patch,  ns)
         end if
         call reset_buffer(mgr%buf(k,t))
      end do
      mgr%has_data(t) = .false.
   end subroutine close_tier

   !=======================================================================================!
   !  FAST tier (sub-daily). Resolve one DIM_SCALAR variable's instantaneous value out of a  !
   !  live fast_sample_t (the temps reuse the SRC_S_* site ids; the fluxes use SRC_F_*). The   !
   !  slab (soil / cohort) sources are resolved directly in output_integrate_fast.             !
   !=======================================================================================!
   pure real(wp) function extract_fast_scalar(source_id, s) result(val)
      integer(ik),        intent(in) :: source_id
      type(fast_sample_t), intent(in) :: s
      select case (source_id)
      case (SRC_F_CAS_TEMP)     ; val = s%cas_temp
      case (SRC_F_SOIL_TEMP_TOP); val = s%soil_temp_top
      case (SRC_F_GPP_RATE)     ; val = s%gpp_rate
      case (SRC_F_LE)           ; val = s%le_flux
      case (SRC_F_H)            ; val = s%h_flux
      case (SRC_F_RNET)         ; val = s%rnet
      case (SRC_F_SW_IN)        ; val = s%sw_in
      case (SRC_F_USTAR)        ; val = s%ustar
      case (SRC_F_AIR_TEMP)     ; val = s%air_temp
      case default              ; val = MISSING_VALUE
      end select
   end function extract_fast_scalar

   !=======================================================================================!
   !  Fold ONE staged sub-step (index isub) into the FAST tier buffers (tier 1). Called from   !
   !  main once per staged sub-step; main closes the tier every fast_interval_steps folds.      !
   !  netCDF-FREE. Reads only the manager's own fast(:) sample + soil/cohort slabs -- never live !
   !  site state (which by then is the post-fast-loop, end-of-slow-step snapshot).               !
   !=======================================================================================!
   subroutine output_integrate_fast(mgr, isub, dt)
      type(output_manager_t), intent(inout) :: mgr
      integer(ik),            intent(in)    :: isub
      real(wp),               intent(in)    :: dt
      integer(ik) :: j, k, src
      if (.not. mgr%enabled) return
      if (mgr%reg%nidx(1) == 0_ik) return
      if (.not. mgr%has_data(1)) mgr%t_open(1) = mgr%fast_time(isub)
      do j = 1_ik, mgr%reg%nidx(1)
         k   = mgr%reg%idx_freq(j, 1_ik)
         src = mgr%reg%var(k)%source_id
         select case (mgr%reg%var(k)%dim)
         case (DIM_SCALAR)
            call integrate_scalar(mgr%buf(k,1), extract_fast_scalar(src, mgr%fast(isub)), dt)
         case (DIM_SOIL)
            if (src == SRC_F_SOIL_WATER) then
               call integrate_slab(mgr%buf(k,1), mgr%fast_soil_water(:,isub), mgr%fast_n_soil, dt)
            else
               call integrate_slab(mgr%buf(k,1), mgr%fast_soil_temp(:,isub), mgr%fast_n_soil, dt)
            end if
         case (DIM_COHORT)
            select case (src)
            case (SRC_F_COH_GPP)
               call integrate_slab(mgr%buf(k,1), mgr%fast_coh_gpp(:,isub), mgr%fast_n_cohort, dt)
            case (SRC_F_COH_HEIGHT)
               call integrate_slab(mgr%buf(k,1), mgr%fast_coh_height(:,isub), mgr%fast_n_cohort, dt)
            case default   ! SRC_F_COH_LEAF_TEMP
               call integrate_slab(mgr%buf(k,1), mgr%fast_coh_ltemp(:,isub), mgr%fast_n_cohort, dt)
            end select
         end select
      end do
      mgr%has_data(1) = .true.
   end subroutine output_integrate_fast

end module meds_output_integrate
