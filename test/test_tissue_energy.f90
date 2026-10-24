! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_tissue_energy -- a cohort's tissue energy is its state, and its temperatures are    !
! derived from it (docs/science/vegetation_energy_dynamics.md eqs. (9)-(10)).              !
!                                                                                          !
! 1. The datum: tissue_temp_to_energy is temp_to_internal_energy with liquid water,        !
!    tissue_energy_to_temp is its inverse, and water that comes or goes at the tissue's    !
!    temperature leaves it there.                                                          !
! 2. Holding: wood that grows at its temperature keeps it, its energy rising by the new    !
!    dry matter's heat; a leafless crown has no energy and takes its wood's temperature.   !
! 3. A split keeps the parent's temperatures in both daughters and its energy in total.    !
! 4. The fast commit: with the vegetation-energy mask off, water gained at the tissue's    !
!    temperature leaves the temperature where it was; with it on, the energy gains exactly !
!    the heat the tissue relaxed through and its water's enthalpy; and the temperature it  !
!    commits gives back the energy, so the next step starts where this one ended.          !
!==========================================================================================!
program test_tissue_energy
   use meds_test_assert,                only : check_true, test_report
   use meds_kinds,                      only : wp, ik
   use meds_constants,                  only : cp_liq, tsupercool_liq, tissue_c_wood_dry
   use meds_therm_lib,                  only : tissue_temp_to_energy, tissue_energy_to_temp,     &
                                               temp_to_internal_energy, internal_energy_liquid,  &
                                               tissue_dry_heat_capacity
   use meds_config,                     only : meds_config_t
   use meds_polygon_state_vars,         only : polygon_t
   use meds_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use meds_cohort_state_vars,          only : hold_cohort_tissue_temps, cohort_tissue_temps
   use meds_demography_cohort_fusefiss, only : split_cohorts
   use meds_fast_integ_types,                 only : fast_dynamic_state_t, alloc_dynamic_state
   use meds_plant_energy,                     only : plant_energy_env_t, commit_tissue_energy
   use meds_test_support,               only : build_test_config
   implicit none

   real(wp), parameter :: DT = 900.0_wp
   type(meds_config_t)  :: cfg
   type(polygon_t)      :: poly
   real(wp) :: e, dw, tl, tw, t_wood, e_wood, c_wood, e_total, t_leaf_parent, t_wood_parent
   integer(ik) :: i

   cfg = build_test_config()

   !=== 1. The datum. ===========================================================================!
   e  = tissue_temp_to_energy(5.0e4_wp, 12.0_wp, 297.3_wp)
   dw = 0.7_wp
   call check_true('the datum is the soil''s: temp_to_internal_energy with liquid water',           &
                   abs(e - temp_to_internal_energy(5.0e4_wp, 12.0_wp, 297.3_wp, 1.0_wp)) <= 1.0e-12_wp * abs(e))
   call check_true('tissue_energy_to_temp inverts tissue_temp_to_energy',                           &
                   abs(tissue_energy_to_temp(e, 12.0_wp, 5.0e4_wp) - 297.3_wp) < 1.0e-10_wp)
   call check_true('water that arrives at the tissue''s temperature leaves it there',               &
                   abs(tissue_energy_to_temp(e + dw * internal_energy_liquid(297.3_wp), 12.0_wp + dw, 5.0e4_wp)  &
                       - 297.3_wp) < 1.0e-10_wp)
   call check_true('and so does water that leaves at it',                                           &
                   abs(tissue_energy_to_temp(e - dw * internal_energy_liquid(297.3_wp), 12.0_wp - dw, 5.0e4_wp)  &
                       - 297.3_wp) < 1.0e-10_wp)

   !=== 2. Holding the temperatures through a change of mass. ==================================!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
   call finalize_init(poly)
   call cohort_tissue_temps(poly%cohort, 1_ik, tl, tw)
   call check_true('a cohort is born with the temperatures its energies give',                      &
                   tl == poly%cohort%leaf_temp(1) .and. tw == poly%cohort%wood_temp(1))
   t_wood = poly%cohort%wood_temp(1)
   e_wood = poly%cohort%wood_energy(1)
   c_wood = poly%cohort%wood_carbon(1)
   poly%cohort%wood_carbon(1) = 1.2_wp * c_wood
   call hold_cohort_tissue_temps(poly%cohort, 1_ik)
   call check_true('wood that grows at its temperature keeps it', abs(poly%cohort%wood_temp(1) - t_wood) < 1.0e-10_wp)
   call check_true('and its energy rises by the new dry wood''s heat at that temperature',          &
                   abs(poly%cohort%wood_energy(1) - e_wood                                          &
                       - tissue_dry_heat_capacity(1.0_wp, 0.2_wp * c_wood, tissue_c_wood_dry) * t_wood) &
                   <= 1.0e-9_wp * e_wood)
   poly%cohort%leaf_carbon(1)     = 0.0_wp
   poly%cohort%leaf_water_mass(1) = 0.0_wp
   call hold_cohort_tissue_temps(poly%cohort, 1_ik)
   call check_true('a leafless crown holds no energy', poly%cohort%leaf_energy(1) == 0.0_wp)
   call check_true('and takes its wood''s temperature', poly%cohort%leaf_temp(1) == poly%cohort%wood_temp(1))

   !=== 3. A split. ==============================================================================!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 2.0_wp, 50.0_wp)      ! nplant*leaf_area well above the cap
   call finalize_init(poly)
   poly%cohort%leaf_temp(1) = 301.0_wp
   poly%cohort%wood_temp(1) = 296.0_wp
   call hold_cohort_tissue_temps(poly%cohort, 1_ik)
   t_leaf_parent = poly%cohort%leaf_temp(1)
   t_wood_parent = poly%cohort%wood_temp(1)
   e_total = site_energy(poly)
   call split_cohorts(poly, cfg)
   call check_true('the split made daughters', poly%cohort%n >= 2_ik)
   do i = 1_ik, poly%cohort%n
      call check_true('each daughter keeps its parent''s temperatures',                             &
                      abs(poly%cohort%leaf_temp(i) - t_leaf_parent) < 1.0e-9_wp .and.                 &
                      abs(poly%cohort%wood_temp(i) - t_wood_parent) < 1.0e-9_wp)
   end do
   call check_true('and together they hold the parent''s energy',                                   &
                   abs(site_energy(poly) - e_total) <= 1.0e-12_wp * e_total)

   !=== 4. The fast commit. =====================================================================!
   call commit_case('commit')

   call test_report('test_tissue_energy')

contains

   !----- The tissues' energy of the site's cohorts [J/m2 of the one patch]. ---------------------!
   real(wp) function site_energy(p) result(s)
      type(polygon_t), intent(in) :: p
      integer(ik) :: k
      s = 0.0_wp
      do k = 1_ik, p%cohort%n
         s = s + p%cohort%nplant(k) * (p%cohort%leaf_energy(k) + p%cohort%wood_energy(k))
      end do
   end function site_energy

   !----- One cohort, one plant per m2, whose leaf and wood gain water over the step and relax from  !
   !      300 K and 296 K towards 301.5 K and 296.8 K. --------------------------------------------!
   subroutine commit_case(tag)
      character(len=*), intent(in) :: tag
      real(wp), parameter :: DL = 2.0e3_wp, DW = 4.0e4_wp          ! [J/K] dry heat capacities
      real(wp), parameter :: WL0 = 0.8_wp, WW0 = 9.0_wp            ! [kg] start-of-step water
      real(wp), parameter :: DWL = 0.05_wp, DWW = 0.4_wp           ! [kg] the water each gains
      real(wp), parameter :: TL0 = 300.0_wp, TW0 = 296.0_wp, TL_AVG = 301.5_wp, TW_AVG = 296.8_wp
      type(plant_energy_env_t)    :: tissue
      type(fast_dynamic_state_t)  :: dynamic
      real(wp) :: store0, store1, heat, enth, tl_expect, tw_expect, gap
      allocate(tissue%leaf_hcap_per_dt(1), tissue%wood_hcap_per_dt(1), tissue%t_leaf0(1), tissue%t_wood0(1))
      allocate(tissue%leaf_dry_hcap(1), tissue%wood_dry_hcap(1), tissue%leaf_water_enth(1), tissue%wood_water_enth(1))
      tissue%leaf_dry_hcap    = DL
      tissue%wood_dry_hcap    = DW
      tissue%leaf_hcap_per_dt = (DL + WL0 * cp_liq) / DT
      tissue%wood_hcap_per_dt = (DW + WW0 * cp_liq) / DT
      tissue%t_leaf0          = TL0
      tissue%t_wood0          = TW0
      tissue%leaf_water_enth  = DWL / DT * internal_energy_liquid(TL0)
      tissue%wood_water_enth  = DWW / DT * internal_energy_liquid(TW0)
      call alloc_dynamic_state(dynamic, 1_ik, 295.0_wp, 0.008_wp, 400.0_wp, TL0)
      dynamic%leaf_energy     = tissue_temp_to_energy(DL, WL0, TL0)
      dynamic%wood_energy     = tissue_temp_to_energy(DW, WW0, TW0)
      dynamic%leaf_water_mass = WL0 + DWL                          ! the network's end-of-step water
      dynamic%wood_water_mass = WW0 + DWW
      call commit_tissue_energy(tissue, [1.0_wp], DT, [TL_AVG], [TW_AVG], 1_ik, dynamic%leaf_water_mass,       &
                                dynamic%wood_water_mass, dynamic%leaf_energy, dynamic%wood_energy,              &
                                dynamic%leaf_temp, dynamic%wood_temp, store0, store1)
      heat = (DL + WL0 * cp_liq) * (TL_AVG - TL0) + (DW + WW0 * cp_liq) * (TW_AVG - TW0)
      tl_expect = ((DL + WL0 * cp_liq) * TL_AVG + DWL * cp_liq * TL0) / (DL + (WL0 + DWL) * cp_liq)
      tw_expect = ((DW + WW0 * cp_liq) * TW_AVG + DWW * cp_liq * TW0) / (DW + (WW0 + DWW) * cp_liq)
      enth = DWL * internal_energy_liquid(TL0) + DWW * internal_energy_liquid(TW0)
      call check_true(tag//': the tissues gain the heat they relaxed through and their water''s enthalpy', &
                      abs((store1 - store0) - (heat + enth)) <= 1.0e-12_wp * abs(heat + enth))
      call check_true(tag//': the gained water arrived at the tissues'' start-of-step temperatures',        &
                      abs(dynamic%leaf_temp(1) - tl_expect) < 1.0e-9_wp .and.                            &
                      abs(dynamic%wood_temp(1) - tw_expect) < 1.0e-9_wp)
      call check_true(tag//': the committed temperatures give back the energies, so the next step '//      &
                      'starts where this one ended',                                                      &
                      abs(tissue_temp_to_energy(DL, WL0 + DWL, dynamic%leaf_temp(1)) - dynamic%leaf_energy(1))    &
                      <= 1.0e-12_wp * abs(dynamic%leaf_energy(1)) .and.                                   &
                      abs(tissue_temp_to_energy(DW, WW0 + DWW, dynamic%wood_temp(1)) - dynamic%wood_energy(1))    &
                      <= 1.0e-12_wp * abs(dynamic%wood_energy(1)))
      !----- What heat capacity x temperature, valued with the water at its liquid enthalpy, would  !
      !      re-value between the two steps: dW*c_liq*(T1 - T0 + T_sc). --------------------------!
      gap = DWW * cp_liq * (dynamic%wood_temp(1) - TW0 + tsupercool_liq)
      print '(a,a,es10.3,a)', '   (', tag, gap, ' J the cap*T bookkeeping would re-value the wood by between steps)'
   end subroutine commit_case

end program test_tissue_energy
