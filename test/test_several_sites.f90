! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_several_sites -- a polygon of several sites.                                        !
!                                                                                          !
! 1. Two identical sites of areas 0.3 and 0.7 each evolve as a polygon of one site does,    !
!    over 60 days of fast and slow steps on a short calendar -- a month boundary every 15    !
!    days, a year boundary every 30, each with a disturbance and patch fusion: every site's !
!    stand is the one-site stand, and the polygon totals are the one site's to round-off.   !
! 2. Two sites on different soils each evolve as a one-site polygon on its own soil, and   !
!    their soil columns differ.                                                             !
! 3. The patch loop gives the same polygon on one thread and on four.                       !
!==========================================================================================!
program test_several_sites
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : rho_h2o
   use meds_config,               only : meds_config_t
   use meds_polygon_state_vars,   only : polygon_t
   use meds_init,                 only : init_bare_ground, init_bare_sites, add_cohort, finalize_init
   use meds_fast_dynamics,        only : fast_context_t, build_fast_context, init_fast_reservoirs
   use meds_stepper,              only : advance_one_step, advance_boundary
   use meds_diagnostic_reduce,    only : total_agb, total_nplant, total_area
   use meds_test_support,         only : banner, build_test_config, check, check_close, test_report
   implicit none

   integer(ik), parameter :: N_DAY = 60_ik            ! two short years of two 15-day months
   real(wp),    parameter :: RTOL  = 1.0e-12_wp
   type(meds_config_t)  :: cfg, cfg_shallow
   type(fast_context_t) :: ctx, ctx4
   type(polygon_t)      :: one, two, two_mt, one_shallow, mixed
   integer(ik)          :: is

   call banner('several sites per polygon')
   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   call build_fast_context(cfg, ctx)
   ctx4 = ctx
   ctx4%patch_threads = 4_ik
   cfg_shallow = cfg
   cfg_shallow%soil_column%depth = 0.6_wp

   !=== 1. Two identical sites, 0.3 and 0.7 of the polygon. ================================!
   call stand(one, cfg, [cfg], [1.0_wp])
   call stand(two, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call stand(two_mt, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call run(one, ctx)
   call run(two, ctx)
   call run(two_mt, ctx4)
   call check(one%patch%n > 2_ik, 'the year boundary added a disturbance gap to the one site')
   do is = 1_ik, 2_ik
      call check(two%site%patch_count(is) == one%patch%n,                                            &
                 'each identical site has the one site''s patches after the disturbance and fusion')
      call check_close(site_carbon(two, is), site_carbon(one, 1_ik), RTOL * site_carbon(one, 1_ik),   &
                       'each identical site holds the one site''s plant carbon')
      call check_close(site_water(two, is), site_water(one, 1_ik), RTOL * site_water(one, 1_ik),      &
                       'each identical site holds the one site''s soil water')
   end do
   call check(all(two%patch%owner_site(1:two%patch%n) == [(1_ik, is = 1_ik, two%site%patch_count(1)),  &
                                                            (2_ik, is = 1_ik, two%site%patch_count(2))]), &
              'the patch block stays sorted by site through the disturbance and fusion')
   call check_close(total_agb(two), total_agb(one), RTOL * total_agb(one),                          &
                    'the polygon''s AGB is the one site''s to round-off')
   call check_close(total_nplant(two), total_nplant(one), RTOL * total_nplant(one),                 &
                    'the polygon''s stem density is the one site''s to round-off')
   call check_close(total_area(two), 1.0_wp, 1.0e-12_wp, 'the patches still cover the polygon')
   do is = 1_ik, 2_ik
      associate (o => two%site%patch_offset(is), c => two%site%patch_count(is))
         call check_close(sum(two%patch%area(o:o + c - 1_ik)), 1.0_wp, 1.0e-12_wp,                  &
                          'each site''s patch areas still sum to 1')
      end associate
   end do

   !=== 3. One thread or four: the same polygon. ===========================================!
   call check(all(two_mt%patch%area(1:two%patch%n) == two%patch%area(1:two%patch%n)) .and.        &
              site_carbon(two_mt, 1_ik) == site_carbon(two, 1_ik) .and.                             &
              site_carbon(two_mt, 2_ik) == site_carbon(two, 2_ik) .and.                             &
              site_water(two_mt, 2_ik) == site_water(two, 2_ik),                                    &
              'four patch threads give the polygon one thread gives, bit for bit')

   !=== 2. Two soils: each site evolves on its own. =========================================!
   call stand(one_shallow, cfg_shallow, [cfg_shallow], [1.0_wp])
   call stand(mixed, cfg, [cfg, cfg_shallow], [0.5_wp, 0.5_wp])
   call run(one_shallow, ctx)
   call run(mixed, ctx)
   call check(mixed%site%surface(1)%soil%n_active /= mixed%site%surface(2)%soil%n_active,           &
              'the two sites have different soil columns')
   call check(abs(site_water(mixed, 1_ik) - site_water(mixed, 2_ik)) > 1.0_wp,                      &
              'the two soils hold different water')
   call check_close(site_carbon(mixed, 1_ik), site_carbon(one, 1_ik), RTOL * site_carbon(one, 1_ik),  &
                    'the deep site grows the one-site stand on the deep soil')
   call check_close(site_carbon(mixed, 2_ik), site_carbon(one_shallow, 1_ik),                      &
                    RTOL * site_carbon(one_shallow, 1_ik), 'the shallow site grows the one-site stand on its soil')
   call check_close(site_water(mixed, 2_ik), site_water(one_shallow, 1_ik),                        &
                    RTOL * site_water(one_shallow, 1_ik), 'the shallow site keeps the one-site water on its soil')

   call test_report('test_several_sites')

contains

   !----- A two-patch stand in each site: a canopy tree and an understorey in one patch, a     !
   !      sapling in the other. ------------------------------------------------------------!
   subroutine stand(poly, run_cfg, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: run_cfg, site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)
      integer(ik) :: k, ns
      ns = int(size(site_cfg), ik)
      call init_bare_sites(poly, run_cfg, [(2_ik, k = 1_ik, ns)], site_cfg, site_area)
      do k = 0_ik, ns - 1_ik
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 1_ik, 0.05_wp, 30.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 2_ik, 0.20_wp,  4.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 2_ik, 3_ik, 0.50_wp,  2.0_wp)
      end do
      call finalize_init(poly)
      call init_fast_reservoirs(poly, ctx)
   end subroutine stand

   !----- N_DAY daily steps, a month's restructuring every 15 days and a year's every 30. -----!
   subroutine run(poly, fctx)
      type(polygon_t),      intent(inout) :: poly
      type(fast_context_t), intent(in)    :: fctx
      integer(ik) :: day
      do day = 1_ik, N_DAY
         call advance_one_step(poly, cfg, fctx)
         if (mod(day, 15_ik) == 0_ik) call advance_boundary(poly, cfg, .true., mod(day, 30_ik) == 0_ik)
      end do
   end subroutine run

   !----- Plant carbon [kgC/m2 of the site]: its patches' cohorts, weighted by patch area. ----!
   real(wp) function site_carbon(poly, is) result(c)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, i
      c = 0.0_wp
      do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
         do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
            c = c + poly%patch%area(ip) * poly%cohort%nplant(i) * (poly%cohort%leaf_carbon(i)           &
                  + poly%cohort%fineroot_carbon(i) + poly%cohort%wood_carbon(i)                       &
                  + poly%cohort%nonstructural_carbon(i))
         end do
      end do
   end function site_carbon

   !----- Soil water [kg/m2 of the site] over its patches' active layers. ---------------------!
   real(wp) function site_water(poly, is) result(w)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, k
      w = 0.0_wp
      associate (soil => poly%site%surface(is)%soil)
         do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
            do k = 1_ik, soil%n_active
               w = w + poly%patch%area(ip) * poly%patch%soil_water(ip)%theta(k) * soil%dz(k) * rho_h2o
            end do
         end do
      end associate
   end function site_water

end program test_several_sites
