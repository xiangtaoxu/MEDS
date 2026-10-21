! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_output_registry -- the registration list + the [output] config resolution (group    !
! toggles, per-tier enables, per-variable overrides, unknown-key trap;                     !
! docs/science/diagnostics.md §5).                                                         !
!==========================================================================================!
program test_output_registry
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_output_params,   only : FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, N_FREQ, N_GRP, &
                                    GRP_STRUCTURE
   use meds_output_types,    only : output_registry_t, MAX_OUTPUT_VARS, DIM_COHORT, DIM_PATCH
   use meds_output_registry, only : build_output_registry, build_freq_index, find_var_index,                      &
                                    apply_variable_override, apply_group_toggles, parse_stream_mask, freq_bit,    &
                                    OVR_TRUE, OVR_FALSE, OVR_MASK
   use meds_test_support, only : banner, build_test_config, check
   implicit none

   type(meds_config_t)     :: cfg
   type(output_registry_t) :: reg
   integer(ik) :: mask, status, t, k
   character(len=8) :: bad
   logical :: found

   call banner('output_registry')
   cfg = build_test_config(86400.0_wp)     ! 1-day slow step; cfg%output holds the [output] defaults

   !----- Build with the defaults (radiation + ecophys groups off, freq_on=[F,T,T,T]). -----!
   call build_output_registry(reg, cfg)
   !----- Size sanity, not an exact count. A hard-coded total has to be edited every time a       !
   !      variable is added, which trains the reader to update the number without checking what     !
   !      changed -- the assertion then guards nothing. Bound it instead, and let the NAMED          !
   !      presence checks below carry the real meaning.  ---------------------------------------!
   call check(reg%nvar > 80_ik .and. reg%nvar <= MAX_OUTPUT_VARS,                                 &
              'registry populated and within MAX_OUTPUT_VARS')
   call check(find_var_index(reg, 'mean_agb_co') > 0_ik, 'mean_agb_co present')
   call check(find_var_index(reg, 'mean_agb_py')   > 0_ik, 'mean_agb_py present')
   call check(find_var_index(reg, 'mean_gpp_py')   > 0_ik, 'mean_gpp_py (carbon) present')
   call check(find_var_index(reg, 'mean_soilc_slow_py') > 0_ik, 'mean_soilc_slow_py (carbon) present')
   call check(find_var_index(reg, 'mean_rh_py')         > 0_ik, 'mean_rh_py (carbon) present')
   call check(find_var_index(reg, 'mean_soil_temp_layer_si') > 0_ik, 'mean_soil_temp_layer_si (energy) present')
   call check(find_var_index(reg, 'nope')       == 0_ik, 'unknown name -> 0')
   !----- The default groups (docs/science/diagnostics.md §5): all ON but radiation + ecophys,    !
   !      the two heavy evaluation groups. Assert BOTH directions of the toggle, so a default     !
   !      flip cannot pass by accident. ----------------------------------------------------------!
   call check(in_tier(reg, 'mean_soil_temp_layer_si', 3_ik), 'energy ON by default -> mean_soil_temp_layer_si in MONTHLY')
   cfg%output%grp_on(4) = .false.       ! GRP_ENERGY off
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'mean_soil_temp_layer_si', 3_ik), 'energy off -> mean_soil_temp_layer_si in no tier')
   cfg%output%grp_on(4) = .true.
   call build_output_registry(reg, cfg)

   !----- FAST tier off by default. The ANNUAL tier must contain NO cohort/patch-dimensioned      !
   !      variable (docs/science/diagnostics.md §3) -- that invariant is what the check below     !
   !      states, rather than a count that has to be re-derived on every registry edit. ----------!
   call check(reg%nidx(1) == 0_ik, 'FAST tier empty by default (freq_on(FAST)=false)')
   call check(reg%nidx(2) > 0_ik,  'DAILY tier populated')
   call check(reg%nidx(4) > 0_ik,  'ANNUAL tier populated')
   call check(annual_is_site_level(reg), 'ANNUAL tier carries no cohort/patch variable')
   !----- The AXIS toggles suppress a whole trailing dimension without naming variables. -----!
   cfg%output%axis_on(1) = .false.                     ! AXIS_COHORT off
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'mean_agb_co', 2_ik), 'axes_cohort=false -> mean_agb_co in no tier')
   call check(in_tier(reg, 'mean_agb_py', 2_ik),         'axes_cohort=false leaves site scalars alone')
   cfg%output%axis_on(1) = .true.
   call build_output_registry(reg, cfg)

   !----- PFT + DBH-size axes are registered and live at the monthly cadence. -----!
   call check(find_var_index(reg, 'mean_agb_pft_py')  > 0_ik, 'mean_agb_pft_py present')
   call check(find_var_index(reg, 'mean_agb_size_py') > 0_ik, 'mean_agb_size_py present')
   call check(in_tier(reg, 'mean_agb_pft_py',  3_ik), 'mean_agb_pft_py in MONTHLY')
   call check(in_tier(reg, 'mean_agb_size_py', 3_ik), 'mean_agb_size_py in MONTHLY')
   call check(in_tier(reg, 'mean_agb_pft_py',  4_ik), 'mean_agb_pft_py in ANNUAL (site-level axis, allowed)')

   call check(in_tier(reg, 'mean_agb_co', 2_ik),        'mean_agb_co in DAILY')
   call check(.not. in_tier(reg, 'mean_agb_co', 4_ik),  'mean_agb_co NOT in ANNUAL (no cohort dim annually)')
   call check(in_tier(reg, 'mean_agb_py', 4_ik),          'mean_agb_py in ANNUAL')
   call check(in_tier(reg, 'mean_gpp_py', 4_ik),          'mean_gpp_py in ANNUAL')

   !----- Per-tier enable: turn DAILY off -> DAILY empties, MONTHLY unaffected. -----!
   cfg%output%freq_on(2) = .false.
   call build_output_registry(reg, cfg)
   call check(reg%nidx(2) == 0_ik, 'DAILY off -> empty')
   call check(reg%nidx(3) > 0_ik,  'MONTHLY still populated')
   cfg%output%freq_on(2) = .true.

   !----- Group toggle: CARBON off -> mean_gpp_py/mean_npp_py vanish, structure stays. -----!
   cfg%output%grp_on(2) = .false.          ! GRP_CARBON
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'mean_gpp_py', 3_ik), 'carbon off -> mean_gpp_py absent from MONTHLY')
   call check(in_tier(reg, 'mean_agb_co', 2_ik),     'carbon off -> structure mean_agb_co still present')
   cfg%output%grp_on(2) = .true.
   !----- Structure OFF -> no GRP_STRUCTURE variable survives anywhere, and the other groups     !
   !      are untouched. Stated as an invariant over the group field rather than a count.  -----!
   cfg%output%grp_on(1) = .false.
   call build_output_registry(reg, cfg)
   call check(group_absent(reg, GRP_STRUCTURE), 'structure off -> no structure var in any tier')
   call check(in_tier(reg, 'mean_gpp_py', 4_ik),   'structure off leaves carbon alone')
   cfg%output%grp_on(1) = .true.

   !----- Per-variable override: OFF removes it from every tier. -----!
   call build_output_registry(reg, cfg)
   call apply_variable_override(reg, 'mean_growth_avg_co', OVR_FALSE, 0_ik, found)
   call check(found, 'override found mean_growth_avg_co')
   call build_freq_index(reg)
   call check(.not. in_tier(reg, 'mean_growth_avg_co', 3_ik), 'mean_growth_avg_co OFF -> absent from MONTHLY')

   !----- Per-variable override: MASK re-targets streams (mean_basal_area_co M -> D M). -----!
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'mean_basal_area_co', 2_ik), 'mean_basal_area_co default NOT in DAILY')
   call parse_stream_mask('D M', mask, status, bad)
   call check(status == 0_ik, 'parse_stream_mask("D M") ok')
   call apply_variable_override(reg, 'mean_basal_area_co', OVR_MASK, mask, found)
   call build_freq_index(reg)
   call check(in_tier(reg, 'mean_basal_area_co', 2_ik), 'mean_basal_area_co override -> now in DAILY')

   !----- Unknown-key trap: a typo matches no registry variable. -----!
   call apply_variable_override(reg, 'grwoth_avg_cohort', OVR_FALSE, 0_ik, found)
   call check(.not. found, 'typo key -> not found (unknown-key trap)')

   !----- parse_stream_mask flags an unrecognized token. -----!
   call parse_stream_mask('D X', mask, status, bad)
   call check(status == 1_ik .and. trim(bad) == 'X', 'bad stream token flagged')

   write(*,'(a)') 'test_output_registry: ALL PASSED'

contains

   !----- .true. if no variable of `grp` is live in any tier. --------------------------------!
   logical function group_absent(reg, grp) result(ok)
      type(output_registry_t), intent(in) :: reg
      integer(ik),             intent(in) :: grp
      integer(ik) :: t, j, k
      ok = .true.
      do t = 1_ik, N_FREQ
         do j = 1_ik, reg%nidx(t)
            k = reg%idx_freq(j, t)
            if (reg%var(k)%group == grp) ok = .false.
         end do
      end do
   end function group_absent

   !----- The annual stream is SITE-LEVEL ONLY (docs/science/diagnostics.md §3). A cohort/patch   !
   !      window longer than a month would straddle the annual disturbance restructuring, so the   !
   !      slot set it averaged would not be the slot set present at flush.  ----------------------!
   logical function annual_is_site_level(reg) result(ok)
      type(output_registry_t), intent(in) :: reg
      integer(ik) :: j, k
      ok = .true.
      do j = 1_ik, reg%nidx(4_ik)
         k = reg%idx_freq(j, 4_ik)
         if (reg%var(k)%dim == DIM_COHORT .or. reg%var(k)%dim == DIM_PATCH) ok = .false.
      end do
   end function annual_is_site_level

   logical function in_tier(reg, name, tier) result(yes)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik),             intent(in) :: tier
      integer(ik) :: idx, j
      idx = find_var_index(reg, name) ; yes = .false.
      do j = 1_ik, reg%nidx(tier)
         if (reg%idx_freq(j, tier) == idx) then ; yes = .true. ; return ; end if
      end do
   end function in_tier

end program test_output_registry
