!==========================================================================================!
! meds_demography_diagnostics -- pure, side-effect-free reductions over the site_t state.      !
!                                                                                          !
! Site-level totals weight per-patch quantities by patch area. These are used by the demo,   !
! the I/O writer, and the test suite to check conservation, detect NaNs, and report the      !
! evolving stand structure. Lives in the demography library (reads the state types directly).!
!==========================================================================================!
module meds_demography_diagnostics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_demography_types,   only : site_t
   use meds_column_state_types, only : n_soil_layer_max
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private

   public :: total_nplant, total_basal_area, total_agb, total_lai, total_area, mean_dbh
   public :: total_gpp, total_npp, site_soil_temp_column, site_soil_water_column
   public :: count_cohorts, has_nan, print_summary

   real(wp), parameter :: cm2_to_m2 = 1.0e-4_wp

contains

   !----- Site plant number [plant m-2 ground]. -------------------------------------------!
   pure real(wp) function total_nplant(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip) * sum(site%cohort%nplant(i0:i1))
      end do
   end function total_nplant

   !----- Site basal area [m2 m-2]. -------------------------------------------------------!
   pure real(wp) function total_basal_area(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip)                                        &
                            * sum(site%cohort%nplant(i0:i1) * site%cohort%basal_area(i0:i1)) * cm2_to_m2
      end do
   end function total_basal_area

   !----- Site aboveground biomass [kgC m-2 ground]. --------------------------------------!
   pure real(wp) function total_agb(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip)                                        &
                            * sum(site%cohort%nplant(i0:i1) * site%cohort%agb(i0:i1))
      end do
   end function total_agb

   !----- Site leaf area index [m2 m-2 ground]. -------------------------------------------!
   pure real(wp) function total_lai(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip)                                        &
                            * sum(site%cohort%nplant(i0:i1) * site%cohort%leaf_area(i0:i1))
      end do
   end function total_lai

   pure real(wp) function total_area(site) result(tot)
      type(site_t), intent(in) :: site
      tot = sum(site%patch%area(1:site%patch%n))
   end function total_area

   !----- Basal-area-weighted mean diameter [cm] (a simple structural index). -------------!
   pure real(wp) function mean_dbh(site) result(dm)
      type(site_t), intent(in) :: site
      real(wp)    :: wsum, dsum, w
      integer(ik) :: ip, i0, i1, i
      wsum = 0.0_wp ; dsum = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            w    = site%patch%area(ip) * site%cohort%nplant(i) * site%cohort%basal_area(i)
            wsum = wsum + w
            dsum = dsum + w * site%cohort%dbh(i)
         end do
      end do
      dm = dsum / max(wsum, tiny_num)
   end function mean_dbh

   !----- Site GROSS GPP accumulated over the slow step [kgC m-2 ground]. The per-cohort           !
   !      gpp_accum is populated by the fast biophysics loop (0 when it is off); summing it here      !
   !      over the period (AGG_SUM in the registry) gives a period carbon TOTAL, not a rate.          !
   pure real(wp) function total_gpp(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip)                                             &
                            * sum(site%cohort%nplant(i0:i1) * site%cohort%gpp_accum(i0:i1))
      end do
   end function total_gpp

   !----- Site NPP over the slow step [kgC m-2 ground] = GPP - autotrophic maintenance resp        !
   !      (leaf + stem + fine-root), all per-slow-step accumulators from the fast loop.              !
   pure real(wp) function total_npp(site) result(tot)
      type(site_t), intent(in) :: site
      integer(ik) :: ip, i0, i1
      tot = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip) ; i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         if (i1 >= i0) tot = tot + site%patch%area(ip)                                             &
              * sum(site%cohort%nplant(i0:i1) * (site%cohort%gpp_accum(i0:i1)                      &
                    - site%cohort%leaf_resp_accum(i0:i1) - site%cohort%stem_resp_accum(i0:i1)      &
                    - site%cohort%root_resp_accum(i0:i1)))
      end do
   end function total_npp

   !----- Area-weighted SITE soil temperature column [K] (average of the per-patch columns). ------!
   !      Written on the fixed DIM_SOIL axis (n_soil_layer_max); meaningful once the fast loop runs. !
   pure subroutine site_soil_temp_column(site, out, n)
      type(site_t), intent(in)  :: site
      real(wp),     intent(out) :: out(:)
      integer(ik),  intent(out) :: n
      integer(ik) :: ip
      real(wp)    :: wsum
      n = n_soil_layer_max
      out(1:n) = 0.0_wp ; wsum = 0.0_wp
      do ip = 1_ik, site%patch%n
         out(1:n) = out(1:n) + site%patch%area(ip) * site%patch%soil_e(ip)%soil_temp(1:n)
         wsum = wsum + site%patch%area(ip)
      end do
      if (wsum > tiny_num) out(1:n) = out(1:n) / wsum
   end subroutine site_soil_temp_column

   !----- Area-weighted SITE volumetric soil-moisture column [m3/m3]. -----------------------------!
   pure subroutine site_soil_water_column(site, out, n)
      type(site_t), intent(in)  :: site
      real(wp),     intent(out) :: out(:)
      integer(ik),  intent(out) :: n
      integer(ik) :: ip
      real(wp)    :: wsum
      n = n_soil_layer_max
      out(1:n) = 0.0_wp ; wsum = 0.0_wp
      do ip = 1_ik, site%patch%n
         out(1:n) = out(1:n) + site%patch%area(ip) * site%patch%soil_w(ip)%theta(1:n)
         wsum = wsum + site%patch%area(ip)
      end do
      if (wsum > tiny_num) out(1:n) = out(1:n) / wsum
   end subroutine site_soil_water_column

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

end module meds_demography_diagnostics
