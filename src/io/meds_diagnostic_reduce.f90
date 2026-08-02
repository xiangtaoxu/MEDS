!==========================================================================================!
! meds_diagnostic_reduce -- generic WEIGHTED aggregation of the demographic state across scales. !
! Replaces meds_output_diagnostics (the bag of ~20 hand-written total_* patch loops).             !
!                                                                                          !
! Stage [3] of the diagnostic wall (MEDS_IO_V01_PLAN.md section 3.1). The idea it replaces: each    !
! site-level number used to be its own subroutine that re-walked the patch/CSR loop with its own    !
! weighting inlined, so a new one cost ~10 lines of boilerplate and a per-PATCH or per-PFT twin      !
! of the same quantity was simply not expressible. Here there is ONE reduction, parameterized by     !
! an explicit WEIGHT KIND, and every scale is a different entry point into it -- which is what lets  !
! one registry line emit a quantity's cohort / patch / site / PFT / size-class variants.             !
!                                                                                          !
! THE TWO RULES THAT MAKE AGGREGATION CORRECT (getting either wrong is the classic diagnostics bug): !
!                                                                                          !
!   EXTENSIVE vs INTENSIVE. An extensive per-plant quantity (agb [kgC/plant], leaf_area [m2/plant])  !
!   aggregates as a weighted SUM and lands in per-ground-area units. An intensive one (dbh,          !
!   leaf_temp, psi_leaf, gsw) aggregates as a weighted MEAN, and WHICH weight is a physical           !
!   statement, not a detail: basal area for dbh, leaf area for canopy fluxes and temperatures,        !
!   nplant for demographic rates. Both the weight kind and the sum/mean choice are carried as DATA     !
!   on the variable descriptor, so they are declared once at registration and cannot drift.            !
!                                                                                          !
!   EMPTY SETS. A patch with no cohorts, a PFT with no members, a size class with no stems -- every    !
!   one of those has a zero denominator. They return `valid = .false.` and the caller emits            !
!   _FillValue; they never return 0/0, and a MEAN never silently reports 0 for "nothing here".         !
!                                                                                          !
! ASSOCIATION ORDER IS PART OF THE CONTRACT. The reductions below deliberately use the same           !
! `sum(array-expression)` construct, in the same order, as the total_* routines they replace, so the   !
! port is bit-identical rather than merely close (the P0 acceptance test asserts exactly that). Do     !
! not "tidy" a sum() into an explicit loop here -- the compiler is free to associate the two            !
! differently, and the last bit is what the gate tests.                                                !
!                                                                                          !
! netCDF-FREE: links meds_core (site_t) + meds_shared. Part of the meds_io_prep target.                !
!==========================================================================================!
module meds_diagnostic_reduce
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_core_state_types,   only : site_t
   use meds_column_state_types, only : n_soil_layer_max
   use meds_diagnostic_kernels, only : dbh_class_index
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private

   !----- The generic API (stage 3). -------------------------------------------------------!
   public :: W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA, W_AGB
   public :: cohort_weight, reduce_cohort_to_site, reduce_cohort_to_patch
   public :: reduce_cohort_to_pft, reduce_cohort_to_size
   public :: reduce_patch_to_site, reduce_patch_column_to_site, gather_patch_columns
   !----- Convenience site totals, kept as thin wrappers so existing callers (meds_main, the   !
   !      legacy meds_io writer, the conservation tests) do not change.                         !
   public :: total_nplant, total_basal_area, total_agb, total_lai, total_area, mean_dbh
   public :: total_gpp, total_npp, site_soil_temp_column, site_soil_water_column
   public :: mean_can_temp, mean_soil_temp_top, total_et
   public :: total_soilc_fast_grnd, total_soilc_fast_soil, total_soilc_struct_grnd,               &
             total_soilc_struct_soil, total_soilc_microbial, total_soilc_slow,                    &
             total_soilc_passive, total_rh
   public :: count_cohorts, has_nan, print_summary

   !----- cm2 -> m2, for basal area (stored per plant in cm2, reported per ground area in m2). !
   real(wp), parameter, public :: cm2_to_m2 = 1.0e-4_wp

   !=======================================================================================!
   !  WEIGHT KINDS. A weight is a PER-COHORT quantity w_i; patch area is applied separately    !
   !  at the patch -> site step, never folded in here, so the same w works at every scale.      !
   !                                                                                          !
   !    patch value (sum)  = SUM_{i in p} w_i * x_i                                             !
   !    patch value (mean) = SUM_{i in p} w_i * x_i / SUM_{i in p} w_i                          !
   !    site  value (sum)  = SUM_p area_p * (patch sum)_p                                       !
   !    site  value (mean) = SUM_p area_p * (patch num)_p / SUM_p area_p * (patch den)_p         !
   !=======================================================================================!
   integer(ik), parameter :: W_NONE       = 0_ik  !< w = 1        (count-like; x already per ground area)
   integer(ik), parameter :: W_NPLANT     = 1_ik  !< w = nplant   (the standard EXTENSIVE per-plant weight)
   integer(ik), parameter :: W_LEAF_AREA  = 2_ik  !< w = nplant*leaf_area (canopy-flux / leaf-state weighting)
   integer(ik), parameter :: W_BASAL_AREA = 3_ik  !< w = nplant*basal_area (structural weighting; dbh)
   integer(ik), parameter :: W_AGB        = 4_ik  !< w = nplant*agb (biomass weighting)

contains

   !=======================================================================================!
   !  Weights                                                                               !
   !=======================================================================================!

   !----- Fill w(1:n) with the per-cohort weight of `wkind` over the WHOLE site SoA. --------!
   pure subroutine cohort_weight(site, wkind, w, n)
      type(site_t), intent(in)  :: site
      integer(ik),  intent(in)  :: wkind
      real(wp),     intent(out) :: w(:)
      integer(ik),  intent(out) :: n
      n = site%cohort%n
      if (n <= 0_ik) return
      select case (wkind)
      case (W_NPLANT)     ; w(1:n) = site%cohort%nplant(1:n)
      case (W_LEAF_AREA)  ; w(1:n) = site%cohort%nplant(1:n) * site%cohort%leaf_area(1:n)
      case (W_BASAL_AREA) ; w(1:n) = site%cohort%nplant(1:n) * site%cohort%basal_area(1:n)
      case (W_AGB)        ; w(1:n) = site%cohort%nplant(1:n) * site%cohort%agb(1:n)
      case default        ; w(1:n) = 1.0_wp                       ! W_NONE
      end select
   end subroutine cohort_weight

   !=======================================================================================!
   !  cohort -> site                                                                        !
   !=======================================================================================!

   !----- Reduce a per-cohort field x(1:ncohort) to one site scalar.                          !
   !                                                                                          !
   !      `scale` multiplies the per-patch contribution BEFORE it is accumulated, matching the  !
   !      `area * sum(...) * cm2_to_m2` order the total_* routines used -- applying it once at   !
   !      the end instead would change the rounding and break the bit-identity gate.            !
   pure subroutine reduce_cohort_to_site(site, x, wkind, mean, val, valid, scale)
      type(site_t),       intent(in)  :: site
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: val
      logical,            intent(out) :: valid
      real(wp), optional, intent(in)  :: scale
      real(wp)    :: w(max(site%cohort%n, 1_ik)), num, den, sc
      integer(ik) :: ip, i0, i1, nw
      sc = 1.0_wp ; if (present(scale)) sc = scale
      call cohort_weight(site, wkind, w, nw)
      num = 0.0_wp ; den = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 < i0) cycle
         num = num + site%patch%area(ip) * sum(w(i0:i1) * x(i0:i1)) * sc
         if (mean) den = den + site%patch%area(ip) * sum(w(i0:i1))
      end do
      if (mean) then
         valid = (abs(den) > tiny_num)
         val   = 0.0_wp ; if (valid) val = num / den
      else
         valid = .true. ; val = num
      end if
   end subroutine reduce_cohort_to_site

   !=======================================================================================!
   !  cohort -> patch                                                                       !
   !=======================================================================================!

   !----- Reduce a per-cohort field to one value per patch. NO area factor (a patch value is  !
   !      already per unit of ITS OWN ground). An empty patch is valid=.false. -> _FillValue.   !
   pure subroutine reduce_cohort_to_patch(site, x, wkind, mean, out, valid, n_out, scale)
      type(site_t),       intent(in)  :: site
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp)    :: w(max(site%cohort%n, 1_ik)), num, den, sc
      integer(ik) :: ip, i0, i1, nw
      sc = 1.0_wp ; if (present(scale)) sc = scale
      call cohort_weight(site, wkind, w, nw)
      n_out = site%patch%n
      out(1:max(n_out,1_ik)) = 0.0_wp ; valid(1:max(n_out,1_ik)) = .false.
      do ip = 1_ik, n_out
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 < i0) cycle                                   ! empty patch -> _FillValue
         num = sum(w(i0:i1) * x(i0:i1)) * sc
         if (mean) then
            den = sum(w(i0:i1))
            if (abs(den) > tiny_num) then ; out(ip) = num / den ; valid(ip) = .true. ; end if
         else
            out(ip) = num ; valid(ip) = .true.
         end if
      end do
   end subroutine reduce_cohort_to_patch

   !=======================================================================================!
   !  cohort -> PFT   (DIM_PFT)                                                             !
   !=======================================================================================!

   !----- Reduce a per-cohort field onto the PFT axis, area-weighted across patches so the    !
   !      result is per m2 of SITE ground (hence sum_pft == the site total, exactly).           !
   pure subroutine reduce_cohort_to_pft(site, x, wkind, mean, n_pft, out, valid, n_out, scale)
      type(site_t),       intent(in)  :: site
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind, n_pft
      logical,            intent(in)  :: mean
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp)    :: w(max(site%cohort%n, 1_ik)), sc, aw
      real(wp)    :: num(max(n_pft,1_ik)), den(max(n_pft,1_ik))
      integer(ik) :: ip, i0, i1, i, k, nw
      sc = 1.0_wp ; if (present(scale)) sc = scale
      call cohort_weight(site, wkind, w, nw)
      n_out = n_pft
      num = 0.0_wp ; den = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         aw = site%patch%area(ip)
         do i = i0, i1
            k = site%cohort%pft(i)
            if (k < 1_ik .or. k > n_pft) cycle
            num(k) = num(k) + aw * w(i) * x(i) * sc
            den(k) = den(k) + aw * w(i)
         end do
      end do
      do k = 1_ik, n_pft
         if (mean) then
            valid(k) = (abs(den(k)) > tiny_num)
            out(k)   = 0.0_wp ; if (valid(k)) out(k) = num(k) / den(k)
         else
            !----- A SUM over a PFT with no members is a true 0 (that PFT contributes nothing),  !
            !      not a missing value -- reporting _FillValue there would break the closure       !
            !      sum_pft(agb_pft) == agb_site the moment a PFT goes locally extinct.             !
            valid(k) = .true. ; out(k) = num(k)
         end if
      end do
   end subroutine reduce_cohort_to_pft

   !=======================================================================================!
   !  cohort -> DBH size class   (DIM_SIZE)                                                 !
   !=======================================================================================!

   !----- Reduce a per-cohort field onto the DBH-class axis, ED2-style: bin by the cohort's   !
   !      MEAN dbh (never split a cohort across bins), area-weight across patches. Same          !
   !      empty-bin convention as the PFT axis above, for the same closure reason.               !
   pure subroutine reduce_cohort_to_size(site, x, wkind, mean, edges, n_class, out, valid,     &
                                         n_out, scale)
      type(site_t),       intent(in)  :: site
      real(wp),           intent(in)  :: x(:)
      integer(ik),        intent(in)  :: wkind, n_class
      logical,            intent(in)  :: mean
      real(wp),           intent(in)  :: edges(:)
      real(wp),           intent(out) :: out(:)
      logical,            intent(out) :: valid(:)
      integer(ik),        intent(out) :: n_out
      real(wp), optional, intent(in)  :: scale
      real(wp)    :: w(max(site%cohort%n, 1_ik)), sc, aw
      real(wp)    :: num(max(n_class,1_ik)), den(max(n_class,1_ik))
      integer(ik) :: ip, i0, i1, i, k, nw
      sc = 1.0_wp ; if (present(scale)) sc = scale
      call cohort_weight(site, wkind, w, nw)
      n_out = n_class
      num = 0.0_wp ; den = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         aw = site%patch%area(ip)
         do i = i0, i1
            k = dbh_class_index(site%cohort%dbh(i), edges, n_class)
            num(k) = num(k) + aw * w(i) * x(i) * sc
            den(k) = den(k) + aw * w(i)
         end do
      end do
      do k = 1_ik, n_class
         if (mean) then
            valid(k) = (abs(den(k)) > tiny_num)
            out(k)   = 0.0_wp ; if (valid(k)) out(k) = num(k) / den(k)
         else
            valid(k) = .true. ; out(k) = num(k)
         end if
      end do
   end subroutine reduce_cohort_to_size

   !=======================================================================================!
   !  patch -> site                                                                         !
   !=======================================================================================!

   !----- Area-weight a per-patch field to one site scalar.                                   !
   !                                                                                          !
   !      `mean = .false.` gives the plain SUM_p area_p * x_p. For a per-ground-area quantity   !
   !      that IS the area-weighted mean, because MEDS renormalizes patch area to 1 on every     !
   !      restructuring (meds_main asserts |sum-1| < 1e-5 at the end of a run). `mean = .true.`  !
   !      divides by SUM area anyway, which differs only in the last bit but is robust to a      !
   !      caller that has not renormalized. Both forms exist because the routines this module    !
   !      replaces used both, and the port is bit-identical.                                     !
   pure subroutine reduce_patch_to_site(site, x, mean, val, valid)
      type(site_t), intent(in)  :: site
      real(wp),     intent(in)  :: x(:)
      logical,      intent(in)  :: mean
      real(wp),     intent(out) :: val
      logical,      intent(out) :: valid
      real(wp)    :: num, den
      integer(ik) :: ip
      num = 0.0_wp ; den = 0.0_wp
      do ip = 1_ik, site%patch%n
         num = num + site%patch%area(ip) * x(ip)
         den = den + site%patch%area(ip)
      end do
      if (mean) then
         valid = (den > tiny_num)
         val   = 0.0_wp ; if (valid) val = num / den
      else
         valid = (site%patch%n > 0_ik) ; val = num
      end if
   end subroutine reduce_patch_to_site

   !----- Area-weight a per-patch COLUMN (soil / snow layers) into one site column. ---------!
   pure subroutine reduce_patch_column_to_site(site, col, nlayer, out, n_out)
      type(site_t), intent(in)  :: site
      real(wp),     intent(in)  :: col(:,:)   !< (layer, patch)
      integer(ik),  intent(in)  :: nlayer
      real(wp),     intent(out) :: out(:)
      integer(ik),  intent(out) :: n_out
      real(wp)    :: wsum
      integer(ik) :: ip
      n_out = nlayer
      out(1:n_out) = 0.0_wp ; wsum = 0.0_wp
      do ip = 1_ik, site%patch%n
         out(1:n_out) = out(1:n_out) + site%patch%area(ip) * col(1:n_out, ip)
         wsum = wsum + site%patch%area(ip)
      end do
      if (wsum > tiny_num) out(1:n_out) = out(1:n_out) / wsum
   end subroutine reduce_patch_column_to_site

   !----- FLATTEN per-patch columns into the 2-D (layer, patch) slab, no reduction at all      !
   !      (DIM_SOIL_PATCH, MEDS_IO_V01_PLAN.md section 3.6.1).                                   !
   !                                                                                          !
   !      The stride is the COMPILE-TIME n_soil_layer_max, never the live layer count. That is   !
   !      the load-bearing detail: striding by a live count would silently re-map every profile   !
   !      the moment the active layer count differed from the one the reader assumed, producing   !
   !      plausible shifted columns rather than an error.                                         !
   pure subroutine gather_patch_columns(site, col, nlayer, out, valid, n_out)
      type(site_t), intent(in)  :: site
      real(wp),     intent(in)  :: col(:,:)   !< (layer, patch)
      integer(ik),  intent(in)  :: nlayer
      real(wp),     intent(out) :: out(:)
      logical,      intent(out) :: valid(:)
      integer(ik),  intent(out) :: n_out
      integer(ik) :: ip, k, base
      n_out = site%patch%n * n_soil_layer_max
      out(1:max(n_out,1_ik))   = 0.0_wp
      valid(1:max(n_out,1_ik)) = .false.
      do ip = 1_ik, site%patch%n
         base = (ip - 1_ik) * n_soil_layer_max
         do k = 1_ik, nlayer
            out(base + k)   = col(k, ip)
            valid(base + k) = .true.
         end do
      end do
   end subroutine gather_patch_columns

   !=======================================================================================!
   !  Convenience site totals (thin wrappers over the generic reduction).                    !
   !                                                                                        !
   !  These are the names meds_main, the legacy meds_io writer and the conservation tests     !
   !  already call. They are kept -- and kept EXACTLY equivalent, association order included --  !
   !  so the reduction refactor is invisible to every existing caller.                          !
   !=======================================================================================!

   !----- Site plant number [plant m-2 ground]. -------------------------------------------!
   pure real(wp) function total_nplant(site) result(tot)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%nplant, W_NONE, .false., tot, ok)
   end function total_nplant

   !----- Site basal area [m2 m-2]: per-plant cm2 -> per-ground m2. ------------------------!
   pure real(wp) function total_basal_area(site) result(tot)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%basal_area, W_NPLANT, .false., tot, ok,    &
                                 scale=cm2_to_m2)
   end function total_basal_area

   !----- Site aboveground biomass [kgC m-2 ground]. --------------------------------------!
   pure real(wp) function total_agb(site) result(tot)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%agb, W_NPLANT, .false., tot, ok)
   end function total_agb

   !----- Site leaf area index [m2 m-2 ground]. -------------------------------------------!
   pure real(wp) function total_lai(site) result(tot)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%leaf_area, W_NPLANT, .false., tot, ok)
   end function total_lai

   pure real(wp) function total_area(site) result(tot)
      type(site_t), intent(in) :: site
      tot = sum(site%patch%area(1:site%patch%n))
   end function total_area

   !----- Basal-area-weighted mean diameter [cm] (a simple structural index). --------------!
   pure real(wp) function mean_dbh(site) result(dm)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%dbh, W_BASAL_AREA, .true., dm, ok)
   end function mean_dbh

   !----- Site GROSS GPP accumulated over the slow step [kgC m-2 ground]. The per-cohort     !
   !      gpp_accum is populated by the fast biophysics loop (0 when it is off).              !
   pure real(wp) function total_gpp(site) result(tot)
      type(site_t), intent(in) :: site
      logical :: ok
      call reduce_cohort_to_site(site, site%cohort%gpp_accum, W_NPLANT, .false., tot, ok)
   end function total_gpp

   !----- Site NPP over the slow step [kgC m-2 ground] = GPP - autotrophic MAINTENANCE resp   !
   !      (leaf + stem + fine-root), all per-slow-step accumulators from the fast loop.        !
   pure real(wp) function total_npp(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp)    :: x(max(site%cohort%n, 1_ik))
      integer(ik) :: n
      logical     :: ok
      n = site%cohort%n
      if (n > 0_ik) x(1:n) = site%cohort%gpp_accum(1:n) - site%cohort%leaf_resp_accum(1:n)    &
                             - site%cohort%stem_resp_accum(1:n) - site%cohort%root_resp_accum(1:n)
      call reduce_cohort_to_site(site, x, W_NPLANT, .false., tot, ok)
   end function total_npp

   !----- Site canopy-air-space temperature [K], area-weighted mean over patches. -----------!
   pure real(wp) function mean_can_temp(site) result(tbar)
      type(site_t), intent(in) :: site
      real(wp)    :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical     :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%cas(ip)%can_temp ; end do
      call reduce_patch_to_site(site, x, .false., tbar, ok)
   end function mean_can_temp

   !----- Site soil-top (layer 1) temperature [K], area-weighted mean over patches. ---------!
   pure real(wp) function mean_soil_temp_top(site) result(tbar)
      type(site_t), intent(in) :: site
      real(wp)    :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical     :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_e(ip)%soil_temp(1) ; end do
      call reduce_patch_to_site(site, x, .false., tbar, ok)
   end function mean_soil_temp_top

   !----- Site evapotranspiration over the slow step [kg/m2 = mm] (fast-loop accumulator). --!
   pure real(wp) function total_et(site) result(tot)
      type(site_t), intent(in) :: site
      tot = site%et_accum
   end function total_et

   !----- Site slow soil-carbon pools [kgC m-2 ground]: area-weighted patch sums. All 0 when  !
   !      soil_carbon_on is off (the pools stay at their allocation-time default).             !
   pure real(wp) function total_soilc_fast_grnd(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%fast_grnd_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_fast_grnd

   pure real(wp) function total_soilc_fast_soil(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%fast_soil_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_fast_soil

   pure real(wp) function total_soilc_struct_grnd(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%struct_grnd_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_struct_grnd

   pure real(wp) function total_soilc_struct_soil(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%struct_soil_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_struct_soil

   pure real(wp) function total_soilc_microbial(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%microbial_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_microbial

   pure real(wp) function total_soilc_slow(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%slow_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_slow

   pure real(wp) function total_soilc_passive(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%soil_carbon(ip)%passive_carbon ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_soilc_passive

   !----- Site heterotrophic respiration over the slow step [kgC m-2 ground]: the fast loop's  !
   !      independently-accumulated matrix Rh. Equals the daily soil_carbon_step's rh_today by  !
   !      construction (the B2 double-count identity), so this is exact, not an approximation.  !
   pure real(wp) function total_rh(site) result(tot)
      type(site_t), intent(in) :: site
      real(wp) :: x(max(site%patch%n, 1_ik))
      integer(ik) :: ip
      logical :: ok
      do ip = 1_ik, site%patch%n ; x(ip) = site%patch%xi_accum(ip)%rh_fast_accum ; end do
      call reduce_patch_to_site(site, x, .false., tot, ok)
   end function total_rh

   !----- Area-weighted SITE soil temperature column [K]. ----------------------------------!
   pure subroutine site_soil_temp_column(site, out, n)
      type(site_t), intent(in)  :: site
      real(wp),     intent(out) :: out(:)
      integer(ik),  intent(out) :: n
      real(wp)    :: col(n_soil_layer_max, max(site%patch%n, 1_ik))
      integer(ik) :: ip
      do ip = 1_ik, site%patch%n
         col(1:n_soil_layer_max, ip) = site%patch%soil_e(ip)%soil_temp(1:n_soil_layer_max)
      end do
      call reduce_patch_column_to_site(site, col, n_soil_layer_max, out, n)
   end subroutine site_soil_temp_column

   !----- Area-weighted SITE volumetric soil-moisture column [m3/m3]. ----------------------!
   pure subroutine site_soil_water_column(site, out, n)
      type(site_t), intent(in)  :: site
      real(wp),     intent(out) :: out(:)
      integer(ik),  intent(out) :: n
      real(wp)    :: col(n_soil_layer_max, max(site%patch%n, 1_ik))
      integer(ik) :: ip
      do ip = 1_ik, site%patch%n
         col(1:n_soil_layer_max, ip) = site%patch%soil_w(ip)%theta(1:n_soil_layer_max)
      end do
      call reduce_patch_column_to_site(site, col, n_soil_layer_max, out, n)
   end subroutine site_soil_water_column

   !=======================================================================================!
   !  Health checks / reporting (unchanged behaviour, moved here with their neighbours).     !
   !=======================================================================================!

   pure integer(ik) function count_cohorts(site) result(nc)
      type(site_t), intent(in) :: site
      nc = site%cohort%n
   end function count_cohorts

   !----- True if any diameter, density, AGB, or wood carbon is NaN. -----------------------!
   pure logical function has_nan(site) result(bad)
      type(site_t), intent(in) :: site
      integer(ik) :: i
      bad = .false.
      do i = 1_ik, site%cohort%n
         if (ieee_is_nan(site%cohort%dbh(i))    .or. ieee_is_nan(site%cohort%nplant(i)) .or.   &
             ieee_is_nan(site%cohort%agb(i))    .or. ieee_is_nan(site%cohort%wood_carbon(i)))  &
            bad = .true.
      end do
   end function has_nan

   !----- Human-readable one-line summary. ------------------------------------------------!
   subroutine print_summary(site, label)
      type(site_t),  intent(in) :: site
      character(len=*), intent(in) :: label
      write(*,'(a,t18,a,i6,a,i4,a,f9.4,a,f8.4,a,f8.3,a,f7.2)')                             &
         trim(label), 'cohorts=', site%cohort%n, '  patches=', site%patch%n,                    &
         '  N=', total_nplant(site), '  LAI=', total_lai(site),                            &
         '  AGB=', total_agb(site), '  Dmean=', mean_dbh(site)
   end subroutine print_summary

end module meds_diagnostic_reduce
