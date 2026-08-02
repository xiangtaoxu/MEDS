!==========================================================================================!
! meds_core_diag_types -- the per-COHORT and per-PATCH fast-loop DIAGNOSTIC accumulators.    !
!                                                                                          !
! Stage [2] of the diagnostic wall (MEDS_IO_V01_PLAN.md sections 3.1, 3.4). The problem it solves:  !
! the fast loop computes ~40 per-cohort and per-patch quantities every dt_fast -- leaf                !
! assimilation, stomatal conductance, tissue water potentials, sapflow, absorbed radiation,           !
! turbulent fluxes, runoff, budget residuals -- and, before this, kept three of them. Sub-daily        !
! resolution existed ONLY inside the fast loop's sub-step, so everything else was recomputed and       !
! discarded ~48 times a day. These blocks are where it lands instead.                                   !
!                                                                                          !
! SIX DESIGN DECISIONS, stated because this is where the change could have gone wrong:                  !
!                                                                                          !
!  1. dt-WEIGHTED SUMS, normalized on read. Every field accumulates x*dt and `w` accumulates dt; the    !
!     reader divides. That is correct under the adaptive integrator (unequal sub-steps) and under        !
!     AGG_TMEAN chaining for free, where a plain running mean would not be.                              !
!                                                                                          !
!  2. RESET ONCE PER SLOW STEP -- the same lifecycle as gpp_accum / et_accum / xi_accum. No second        !
!     reset cadence is invented.                                                                          !
!                                                                                          !
!  3. TRANSIENT: never written to the restart file. A checkpoint is prognostic state at an instant;       !
!     a time-averaged diagnostic is not state. This is what keeps the change out of io_write_state.       !
!                                                                                          !
!  4. THEY STILL RIDE THE LOCKSTEP. The reset is per slow step, but restructuring (fuse / split /         !
!     cull / recruit / disturb) happens INSIDE the slow step, after the fast loop has filled these        !
!     and before the monthly window closes. So slot i must keep meaning cohort i across every             !
!     permutation. That obligation is REAL and is the main cost of this design.                            !
!                                                                                          !
!     It is discharged by STORAGE LAYOUT rather than by discipline: the fields are rows of ONE 2-D        !
!     array v(N_CDIAG, cap), so every lockstep operation is a single whole-array statement that            !
!     cannot omit a field. Adding a diagnostic is a new index parameter plus a fill line -- ZERO           !
!     edits to the reorder machinery. Compare the alternative (25 separate allocatables), which would      !
!     have needed 25 new lines in each of cohort_alloc / ensure_capacity / move_alloc_block /              !
!     cohort_reorder / copy_cohort_slot -- precisely the "forgot to reallocate" class CLAUDE.md names      !
!     as the reason the centralized reorder exists.                                                         !
!                                                                                          !
!  5. ALLOCATED ONLY ON DEMAND. A run that reports no per-cohort ecophysiology pays nothing: main          !
!     sets `active` from the registry, and every entry point is a no-op while it is .false.                 !
!                                                                                          !
!  6. FUSION BLENDS BY THE FIELD'S OWN KIND. Each field declares whether it is INTENSIVE (leaf-area-        !
!     weighted -- a temperature, a conductance, a potential) or EXTENSIVE (nplant-weighted -- a per-plant   !
!     flux) or GROUND-REFERENCED (summed -- already per m2 of patch). Blending a temperature by nplant, or  !
!     an extensive flux by leaf area, is the classic silent fusion bug; the kind is DATA here, so the       !
!     blend cannot disagree with the quantity.                                                               !
!==========================================================================================!
module meds_core_diag_types
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: cohort_diag_block, patch_diag_block
   public :: cohort_diag_alloc, cohort_diag_free, cohort_diag_grow, cohort_diag_reset
   public :: cohort_diag_reorder, cohort_diag_copy_slot, cohort_diag_fuse, cohort_diag_clear_slot
   public :: patch_diag_alloc, patch_diag_free, patch_diag_grow, patch_diag_reset
   public :: patch_diag_reorder, patch_diag_copy_slot, patch_diag_blend, patch_diag_clear_slot
   public :: cohort_diag_value, patch_diag_value

   !----- Fusion kinds (decision 6 above). ---------------------------------------------------!
   integer(ik), parameter, public :: FK_INTENSIVE = 1_ik  !< leaf-area-weighted mean (temps, gs, psi)
   integer(ik), parameter, public :: FK_EXTENSIVE = 2_ik  !< nplant-weighted mean (per-plant fluxes)
   integer(ik), parameter, public :: FK_GROUND    = 3_ik  !< summed (already per m2 of patch ground)

   !==========================================================================================!
   !  PER-COHORT diagnostic fields. Order is free; the parameters are the only thing that binds  !
   !  a row to a meaning.                                                                        !
   !==========================================================================================!
   integer(ik), parameter, public :: CD_ANET         =  1_ik  !< [umol CO2/m2 leaf/s] net assimilation
   integer(ik), parameter, public :: CD_AGROSS       =  2_ik  !< [umol CO2/m2 leaf/s] gross assimilation
   integer(ik), parameter, public :: CD_GSW          =  3_ik  !< [mol H2O/m2 leaf/s]  stomatal conductance
   integer(ik), parameter, public :: CD_GBW          =  4_ik  !< [m/s] leaf boundary-layer conductance
   integer(ik), parameter, public :: CD_CI           =  5_ik  !< [umol/mol] intercellular CO2
   integer(ik), parameter, public :: CD_CS           =  6_ik  !< [umol/mol] leaf-surface CO2
   integer(ik), parameter, public :: CD_RD           =  7_ik  !< [umol/m2 leaf/s] leaf dark respiration
   integer(ik), parameter, public :: CD_TRANSP       =  8_ik  !< [mol H2O/m2 leaf/s] leaf transpiration
   integer(ik), parameter, public :: CD_BETA_STOM    =  9_ik  !< [-] stomatal water-stress limb
   integer(ik), parameter, public :: CD_BETA_NONSTOM = 10_ik  !< [-] non-stomatal (capacity) limb
   integer(ik), parameter, public :: CD_LEAF_TEMP    = 11_ik  !< [K]
   integer(ik), parameter, public :: CD_WOOD_TEMP    = 12_ik  !< [K]
   integer(ik), parameter, public :: CD_LEAF_VPD     = 13_ik  !< [Pa] leaf-to-air VPD
   integer(ik), parameter, public :: CD_PSI_LEAF     = 14_ik  !< [MPa]
   integer(ik), parameter, public :: CD_PSI_WOOD     = 15_ik  !< [MPa]
   integer(ik), parameter, public :: CD_PLC          = 16_ik  !< [-] percent loss of conductance
   integer(ik), parameter, public :: CD_SAPFLOW      = 17_ik  !< [kg/plant/s] wood -> leaf
   integer(ik), parameter, public :: CD_ROOT_UPTAKE  = 18_ik  !< [kg/plant/s] soil -> root
   integer(ik), parameter, public :: CD_ABS_PAR      = 19_ik  !< [W/m2 ground] absorbed PAR
   integer(ik), parameter, public :: CD_ABS_SW       = 20_ik  !< [W/m2 ground] absorbed shortwave
   integer(ik), parameter, public :: CD_ABS_LW       = 21_ik  !< [W/m2 ground] net longwave
   integer(ik), parameter, public :: CD_WIND         = 22_ik  !< [m/s] in-canopy wind
   integer(ik), parameter, public :: CD_GPP_RATE     = 23_ik  !< [umol CO2/plant/s]
   integer(ik), parameter, public :: CD_LEAF_WATER   = 24_ik  !< [kg/plant] internal leaf water
   integer(ik), parameter, public :: CD_WOOD_WATER   = 25_ik  !< [kg/plant] internal wood water
   integer(ik), parameter, public :: N_CDIAG         = 25_ik

   !----- Fusion kind per cohort field. Temperatures, conductances, potentials and the leaf-area- !
   !      normalized fluxes are INTENSIVE (leaf-area-weighted). The per-PLANT quantities (GPP rate, !
   !      sapflow, uptake, tissue water mass) are EXTENSIVE (nplant-weighted), matching how the      !
   !      prognostic twins of exactly these quantities are already fused. The absorbed-radiation      !
   !      terms are per m2 of GROUND, so two cohorts' contributions to the same patch simply add --   !
   !      the same convention leaf_surf_water follows.  -----------------------------------------!
   integer(ik), parameter, public :: CDIAG_FUSE(N_CDIAG) = [                                     &
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE,                    &
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE,                    &
        FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE, FK_INTENSIVE,                    &
        FK_INTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_GROUND,    FK_GROUND,                       &
        FK_GROUND,    FK_INTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE ]

   !==========================================================================================!
   !  PER-COHORT SLOW-LOOP diagnostic fields -- a SECOND block of the same shape, written once per  !
   !  slow step by the vegetation-dynamics driver.                                                  !
   !                                                                                          !
   !  WHY A SEPARATE BLOCK AND NOT MORE ROWS ABOVE. The two have different sample semantics: the     !
   !  fast rows accumulate ~48 sub-step samples per slow step, the slow rows get exactly one. One     !
   !  shared dt weight cannot normalize both.                                                         !
   !                                                                                          !
   !  WHY NOT READ site%deriv DIRECTLY. That was the first implementation and it was WRONG, in a      !
   !  way only the thread-invariance test exposed. `cohort_deriv_block` is documented as TRANSIENT     !
   !  and deliberately NOT lockstep-reordered -- which is fine for its own consumer, since            !
   !  update_cohort_states applies it immediately. But the output tick runs at the END of the step,    !
   !  AFTER the monthly fiss/fuse has permuted the cohort axis, so `deriv(i)` and `cohort(i)` refer     !
   !  to different plants on exactly the boundary steps. The error was invisible at one thread and     !
   !  appeared at four only because thread count perturbs which cohorts fuse. Storing the tendencies    !
   !  in a block that DOES ride the lockstep removes the failure mode rather than timing around it.     !
   !==========================================================================================!
   integer(ik), parameter, public :: CS_DDBH_DT      = 1_ik  !< [cm/yr]         diameter growth
   integer(ik), parameter, public :: CS_DAGB_DT      = 2_ik  !< [kgC/plant/yr]  AGB growth
   integer(ik), parameter, public :: CS_MORT_RATE    = 3_ik  !< [1/yr]          mortality hazard (positive)
   integer(ik), parameter, public :: CS_NPP_LEAF     = 4_ik  !< [kgC/plant/yr]  NPP to leaf
   integer(ik), parameter, public :: CS_NPP_FINEROOT = 5_ik  !< [kgC/plant/yr]  NPP to fine root
   integer(ik), parameter, public :: CS_NPP_WOOD     = 6_ik  !< [kgC/plant/yr]  NPP to wood
   integer(ik), parameter, public :: CS_NPP_STORAGE  = 7_ik  !< [kgC/plant/yr]  NPP to non-structural storage
   integer(ik), parameter, public :: CS_NPP_REPRO    = 8_ik  !< [kgC/plant/yr]  NPP to reproduction
   integer(ik), parameter, public :: CS_GROWTH_RESP  = 9_ik  !< [kgC/plant/yr]  growth respiration
   integer(ik), parameter, public :: N_CSDIAG        = 9_ik

   !----- Fusion kinds for the slow rows. All of them are per-PLANT rates, so all are EXTENSIVE    !
   !      (nplant-weighted) -- the same convention gpp_accum and the maintenance-respiration        !
   !      accumulators already follow in fuse_2_cohorts.  ---------------------------------------!
   integer(ik), parameter, public :: CSDIAG_FUSE(N_CSDIAG) = [                                   &
        FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE,                    &
        FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE, FK_EXTENSIVE ]

   !==========================================================================================!
   !  PER-PATCH diagnostic fields.                                                              !
   !==========================================================================================!
   integer(ik), parameter, public :: PD_LE            =  1_ik  !< [W/m2] latent heat (CAS -> atm)
   integer(ik), parameter, public :: PD_H             =  2_ik  !< [W/m2] sensible heat (CAS -> atm)
   integer(ik), parameter, public :: PD_RNET          =  3_ik  !< [W/m2] net all-wave radiation
   integer(ik), parameter, public :: PD_SW_IN         =  4_ik  !< [W/m2] incident shortwave at canopy top
   integer(ik), parameter, public :: PD_SW_GROUND     =  5_ik  !< [W/m2] shortwave reaching the ground
   integer(ik), parameter, public :: PD_LW_GROUND     =  6_ik  !< [W/m2] net longwave at the ground
   integer(ik), parameter, public :: PD_USTAR         =  7_ik  !< [m/s]
   integer(ik), parameter, public :: PD_GGNET         =  8_ik  !< [m/s] ground conductance
   integer(ik), parameter, public :: PD_ROUGH         =  9_ik  !< [m] roughness length
   integer(ik), parameter, public :: PD_DISPLACE      = 10_ik  !< [m] displacement height
   integer(ik), parameter, public :: PD_CAS_TEMP      = 11_ik  !< [K]
   integer(ik), parameter, public :: PD_CAS_SHV       = 12_ik  !< [kg/kg]
   integer(ik), parameter, public :: PD_CAS_CO2       = 13_ik  !< [umol/mol]
   integer(ik), parameter, public :: PD_GPP           = 14_ik  !< [umol CO2/m2/s] patch GPP rate
   integer(ik), parameter, public :: PD_NEE           = 15_ik  !< [umol CO2/m2/s] net ecosystem exchange
   integer(ik), parameter, public :: PD_TRANSP        = 16_ik  !< [kg/m2/s] canopy transpiration
   integer(ik), parameter, public :: PD_ROOT_UPTAKE   = 17_ik  !< [kg/m2/s] realized root uptake
   integer(ik), parameter, public :: PD_INFILTRATION  = 18_ik  !< [kg/m2/s]
   integer(ik), parameter, public :: PD_DRAINAGE      = 19_ik  !< [kg/m2/s] bottom-face drainage
   integer(ik), parameter, public :: PD_RUNOFF        = 20_ik  !< [kg/m2/s] surface runoff
   integer(ik), parameter, public :: PD_PRECIP        = 21_ik  !< [kg/m2/s] total precipitation
   integer(ik), parameter, public :: PD_GROUND_TEMP   = 22_ik  !< [K] ground/skin temperature
   integer(ik), parameter, public :: PD_RESID_ENERGY  = 23_ik  !< [W/m2] whole-column energy residual
   integer(ik), parameter, public :: PD_RESID_WATER   = 24_ik  !< [kg/m2/s] whole-column water residual
   integer(ik), parameter, public :: N_PDIAG          = 24_ik

   !==========================================================================================!
   !  The blocks themselves. `v` is (field, slot): field-major so a lockstep permutation of the  !
   !  SLOT axis is one contiguous whole-array statement (decision 4). Extraction reads one field  !
   !  across slots, which is strided -- acceptable, because this is the once-per-step diagnostic   !
   !  path, not the inner integrator.                                                              !
   !==========================================================================================!
   type :: cohort_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      integer(ik) :: nfield = 0_ik           !< N_CDIAG (fast block) or N_CSDIAG (slow block)
      logical     :: active = .false.        !< .false. => every entry point is a no-op, nothing allocated
      real(wp), allocatable :: v(:,:)        !< (nfield, cap) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (cap)         running Sum(dt)
   end type cohort_diag_block

   type :: patch_diag_block
      integer(ik) :: n = 0_ik, cap = 0_ik
      logical     :: active = .false.
      real(wp), allocatable :: v(:,:)        !< (N_PDIAG, cap) running Sum(x*dt)
      real(wp), allocatable :: w(:)          !< (cap)          running Sum(dt)
   end type patch_diag_block

contains

   !=======================================================================================!
   !  COHORT block lifecycle. Every routine no-ops while `active` is .false., so a run that   !
   !  reports no per-cohort ecophysiology never allocates and never pays.                     !
   !=======================================================================================!

   subroutine cohort_diag_alloc(d, cap, active, nfield)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: cap
      logical,                 intent(in)    :: active
      integer(ik), optional,   intent(in)    :: nfield   !< default N_CDIAG (the fast block)
      integer(ik) :: nf
      nf = N_CDIAG ; if (present(nfield)) nf = nfield
      d%active = active
      d%nfield = nf
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(nf, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%cap = max(cap, 1_ik) ; d%n = 0_ik
      call cohort_diag_reset(d)
   end subroutine cohort_diag_alloc

   subroutine cohort_diag_free(d)
      type(cohort_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik ; d%n = 0_ik ; d%active = .false.
   end subroutine cohort_diag_free

   !----- Grow to at least `need` slots, preserving the live prefix (the ensure_capacity twin). !
   subroutine cohort_diag_grow(d, need)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      m = d%n
      allocate(v2(int(size(d%v,1),ik), need), w2(need))
      v2 = 0.0_wp ; w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m) ; w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v) ; call move_alloc(w2, d%w)
      d%cap = need
   end subroutine cohort_diag_grow

   !----- Zero the whole block for a fresh slow step. ---------------------------------------!
   subroutine cohort_diag_reset(d)
      type(cohort_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine cohort_diag_reset

   !----- Permute the live prefix in lockstep with the cohort SoA. ONE statement -- it cannot   !
   !      omit a field, which is the whole point of the 2-D layout (decision 4).                 !
   subroutine cohort_diag_reorder(d, perm, m)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: perm(:), m
      if (.not. d%active) return
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine cohort_diag_reorder

   subroutine cohort_diag_copy_slot(d, dst, src)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine cohort_diag_copy_slot

   !----- Zero one slot (a fresh birth, or a reused stale cull slot). -----------------------!
   subroutine cohort_diag_clear_slot(d, i)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp ; d%w(i) = 0.0_wp
   end subroutine cohort_diag_clear_slot

   !----- Blend the donor into the survivor on cohort fusion, each field by its own kind.     !
   !                                                                                          !
   !      `w_int` / `w_ext` are the SAME weights fuse_2_cohorts uses for the prognostic twins of !
   !      these quantities (leaf area for the intensive ones, nplant for the extensive ones), so  !
   !      a diagnostic and its prognostic counterpart can never be fused on different weights.    !
   !      The dt weight `w` is carried across too, so a fused cohort's period mean stays correct.  !
   subroutine cohort_diag_fuse(d, recc, donc, li_r, li_d, np_r, np_d, kinds)
      type(cohort_diag_block), intent(inout) :: d
      integer(ik),             intent(in)    :: recc, donc
      real(wp),                intent(in)    :: li_r, li_d   !< leaf-area weights (nplant*leaf_area)
      real(wp),                intent(in)    :: np_r, np_d   !< nplant weights
      integer(ik),             intent(in)    :: kinds(:)     !< CDIAG_FUSE or CSDIAG_FUSE
      real(wp)    :: wi, we
      integer(ik) :: f
      if (.not. d%active) return
      wi = li_r + li_d ; we = np_r + np_d
      do f = 1_ik, min(int(size(kinds), ik), int(size(d%v, 1), ik))
         select case (kinds(f))
         case (FK_INTENSIVE)
            if (wi > tiny_num) d%v(f, recc) = (li_r * d%v(f, recc) + li_d * d%v(f, donc)) / wi
         case (FK_EXTENSIVE)
            if (we > tiny_num) d%v(f, recc) = (np_r * d%v(f, recc) + np_d * d%v(f, donc)) / we
         case default   ! FK_GROUND -- already per m2 of the SAME patch ground, so they add
            d%v(f, recc) = d%v(f, recc) + d%v(f, donc)
         end select
      end do
      !----- The dt weight is a property of the WINDOW, not of either cohort, and both slots saw  !
      !      the same sub-steps -- so the survivor keeps its own rather than summing (summing would !
      !      halve every normalized mean).  -------------------------------------------------!
   end subroutine cohort_diag_fuse

   !----- Read one field's live prefix, normalized by the dt weight. Slots with no samples     !
   !      return 0 (the caller's `valid` mask, not this, decides _FillValue).                   !
   pure subroutine cohort_diag_value(d, field, x, n)
      type(cohort_diag_block), intent(in)  :: d
      integer(ik),             intent(in)  :: field
      real(wp),                intent(out) :: x(:)
      integer(ik),             intent(out) :: n
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      n = d%n
      do i = 1_ik, n
         if (d%w(i) > tiny_num) then ; x(i) = d%v(field, i) / d%w(i) ; else ; x(i) = 0.0_wp ; end if
      end do
   end subroutine cohort_diag_value

   !=======================================================================================!
   !  PATCH block lifecycle -- the same shape, one axis smaller.                             !
   !=======================================================================================!

   subroutine patch_diag_alloc(d, cap, active)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: cap
      logical,                intent(in)    :: active
      d%active = active
      if (.not. active) return
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      allocate(d%v(N_PDIAG, max(cap,1_ik)), d%w(max(cap,1_ik)))
      d%cap = max(cap, 1_ik) ; d%n = 0_ik
      call patch_diag_reset(d)
   end subroutine patch_diag_alloc

   subroutine patch_diag_free(d)
      type(patch_diag_block), intent(inout) :: d
      if (allocated(d%v)) deallocate(d%v)
      if (allocated(d%w)) deallocate(d%w)
      d%cap = 0_ik ; d%n = 0_ik ; d%active = .false.
   end subroutine patch_diag_free

   subroutine patch_diag_grow(d, need)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: need
      real(wp), allocatable :: v2(:,:), w2(:)
      integer(ik) :: m
      if (.not. d%active) return
      if (need <= d%cap) return
      m = d%n
      allocate(v2(N_PDIAG, need), w2(need))
      v2 = 0.0_wp ; w2 = 0.0_wp
      if (m > 0_ik) then
         v2(:, 1:m) = d%v(:, 1:m) ; w2(1:m) = d%w(1:m)
      end if
      call move_alloc(v2, d%v) ; call move_alloc(w2, d%w)
      d%cap = need
   end subroutine patch_diag_grow

   subroutine patch_diag_reset(d)
      type(patch_diag_block), intent(inout) :: d
      if (.not. d%active) return
      if (allocated(d%v)) d%v = 0.0_wp
      if (allocated(d%w)) d%w = 0.0_wp
   end subroutine patch_diag_reset

   subroutine patch_diag_reorder(d, perm, m)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: perm(:), m
      if (.not. d%active) return
      if (m <= 0_ik) return
      d%v(:, 1:m) = d%v(:, perm(1:m))
      d%w(1:m)    = d%w(perm(1:m))
   end subroutine patch_diag_reorder

   subroutine patch_diag_copy_slot(d, dst, src)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: dst, src
      if (.not. d%active) return
      d%v(:, dst) = d%v(:, src)
      d%w(dst)    = d%w(src)
   end subroutine patch_diag_copy_slot

   subroutine patch_diag_clear_slot(d, i)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: i
      if (.not. d%active) return
      if (i < 1_ik .or. i > d%cap) return
      d%v(:, i) = 0.0_wp ; d%w(i) = 0.0_wp
   end subroutine patch_diag_clear_slot

   !----- Area-weighted blend on patch fusion. Every patch diagnostic is per m2 of ITS OWN     !
   !      ground, so all of them blend the same way -- there is no extensive/intensive split at  !
   !      this level, which is why there is no kind array here.                                  !
   subroutine patch_diag_blend(d, recp, donp, area_r, area_d)
      type(patch_diag_block), intent(inout) :: d
      integer(ik),            intent(in)    :: recp, donp
      real(wp),               intent(in)    :: area_r, area_d
      real(wp) :: wtot
      if (.not. d%active) return
      wtot = area_r + area_d
      if (wtot <= tiny_num) return
      d%v(:, recp) = (area_r * d%v(:, recp) + area_d * d%v(:, donp)) / wtot
   end subroutine patch_diag_blend

   pure subroutine patch_diag_value(d, field, x, n)
      type(patch_diag_block), intent(in)  :: d
      integer(ik),            intent(in)  :: field
      real(wp),               intent(out) :: x(:)
      integer(ik),            intent(out) :: n
      integer(ik) :: i
      n = 0_ik
      if (.not. d%active) return
      n = d%n
      do i = 1_ik, n
         if (d%w(i) > tiny_num) then ; x(i) = d%v(field, i) / d%w(i) ; else ; x(i) = 0.0_wp ; end if
      end do
   end subroutine patch_diag_value

end module meds_core_diag_types
