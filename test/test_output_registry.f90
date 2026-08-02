!==========================================================================================!
! test_output_registry -- the registration list + the §6.1 config resolution (group toggles,     !
! per-tier enables, per-variable overrides, unknown-key trap). MEDS_IO_DESIGN.md test 4.          !
!==========================================================================================!
program test_output_registry
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_output_config,   only : FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, N_FREQ, N_GRP, &
                                    GRP_STRUCTURE
   use meds_output_types,    only : output_registry_t, MAX_OUTPUT_VARS, DIM_COHORT, DIM_PATCH
   use meds_output_registry, only : build_output_registry, build_freq_index, find_var_index,     &
                                    apply_variable_override, apply_group_toggles, parse_stream_mask, &
                                    freq_bit, OVR_TRUE, OVR_FALSE, OVR_MASK
   use meds_test_support,    only : check, banner, build_test_config
   implicit none

   type(meds_config_t)     :: cfg
   type(output_registry_t) :: reg
   integer(ik) :: mask, status, t
   character(len=8) :: bad
   logical :: found

   call banner('output_registry')
   cfg = build_test_config(86400.0_wp)     ! 1-day slow step; cfg%output holds the [output] defaults

   !----- Build with the defaults (grp_on=[T,T,F,F,F], freq_on=[F,T,T,T]). -----!
   call build_output_registry(reg, cfg)
   !----- Size sanity, not an exact count. A hard-coded total has to be edited every time a       !
   !      variable is added, which trains the reader to update the number without checking what     !
   !      changed -- the assertion then guards nothing. Bound it instead, and let the NAMED          !
   !      presence checks below carry the real meaning.  ---------------------------------------!
   call check(reg%nvar > 80_ik .and. reg%nvar <= MAX_OUTPUT_VARS,                                 &
              'registry populated and within MAX_OUTPUT_VARS')
   call check(find_var_index(reg, 'agb_cohort') > 0_ik, 'agb_cohort present')
   call check(find_var_index(reg, 'agb_site')   > 0_ik, 'agb_site present')
   call check(find_var_index(reg, 'gpp_site')   > 0_ik, 'gpp_site (carbon) present')
   call check(find_var_index(reg, 'soilc_slow_site') > 0_ik, 'soilc_slow_site (carbon, B3) present')
   call check(find_var_index(reg, 'rh_site')         > 0_ik, 'rh_site (carbon, B3) present')
   call check(find_var_index(reg, 'soil_temp_site') > 0_ik, 'soil_temp_site (energy) present')
   call check(find_var_index(reg, 'nope')       == 0_ik, 'unknown name -> 0')
   !----- v0.1 defaults (section 8 D6): structure/carbon/water/energy/numerics/biogeochem ON,      !
   !      radiation + ecophys OFF -- the two heavy evaluation groups. Assert BOTH directions of the  !
   !      toggle, so a default flip cannot pass by accident.  -----------------------------------!
   call check(in_tier(reg, 'soil_temp_site', 3_ik), 'energy ON by default -> soil_temp_site in MONTHLY')
   cfg%output%grp_on(4) = .false.       ! GRP_ENERGY off
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'soil_temp_site', 3_ik), 'energy off -> soil_temp_site in no tier')
   cfg%output%grp_on(4) = .true.
   call build_output_registry(reg, cfg)

   !----- FAST tier off by default. The ANNUAL tier must contain NO cohort/patch-dimensioned      !
   !      variable (the section 3.1 guard) -- that invariant is what the check below states, rather  !
   !      than a count that has to be re-derived on every registry edit.  ----------------------!
   call check(reg%nidx(1) == 0_ik, 'FAST tier empty by default (freq_on(FAST)=false)')
   call check(reg%nidx(2) > 0_ik,  'DAILY tier populated')
   call check(reg%nidx(4) > 0_ik,  'ANNUAL tier populated')
   call check(annual_is_site_level(reg), 'ANNUAL tier carries no cohort/patch variable')
   !----- The AXIS toggles suppress a whole trailing dimension without naming variables. -----!
   cfg%output%axis_on(1) = .false.                     ! AXIS_COHORT off
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'agb_cohort', 2_ik), 'axes_cohort=false -> agb_cohort in no tier')
   call check(in_tier(reg, 'agb_site', 2_ik),         'axes_cohort=false leaves site scalars alone')
   cfg%output%axis_on(1) = .true.
   call build_output_registry(reg, cfg)

   !----- PFT + DBH-size axes are registered and live at the monthly cadence. -----!
   call check(find_var_index(reg, 'agb_pft')  > 0_ik, 'agb_pft present')
   call check(find_var_index(reg, 'agb_size') > 0_ik, 'agb_size present')
   call check(in_tier(reg, 'agb_pft',  3_ik), 'agb_pft in MONTHLY')
   call check(in_tier(reg, 'agb_size', 3_ik), 'agb_size in MONTHLY')
   call check(in_tier(reg, 'agb_pft',  4_ik), 'agb_pft in ANNUAL (site-level axis, allowed)')

   call check(in_tier(reg, 'agb_cohort', 2_ik),        'agb_cohort in DAILY')
   call check(.not. in_tier(reg, 'agb_cohort', 4_ik),  'agb_cohort NOT in ANNUAL (no cohort dim annually)')
   call check(in_tier(reg, 'agb_site', 4_ik),          'agb_site in ANNUAL')
   call check(in_tier(reg, 'gpp_site', 4_ik),          'gpp_site in ANNUAL')

   !----- Per-tier enable: turn DAILY off -> DAILY empties, MONTHLY unaffected. -----!
   cfg%output%freq_on(2) = .false.
   call build_output_registry(reg, cfg)
   call check(reg%nidx(2) == 0_ik, 'DAILY off -> empty')
   call check(reg%nidx(3) > 0_ik,  'MONTHLY still populated')
   cfg%output%freq_on(2) = .true.

   !----- Group toggle: CARBON off -> gpp_site/npp_site vanish, structure stays. -----!
   cfg%output%grp_on(2) = .false.          ! GRP_CARBON
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'gpp_site', 3_ik), 'carbon off -> gpp_site absent from MONTHLY')
   call check(in_tier(reg, 'agb_cohort', 2_ik),     'carbon off -> structure agb_cohort still present')
   cfg%output%grp_on(2) = .true.
   !----- Structure OFF -> no GRP_STRUCTURE variable survives anywhere, and the other groups     !
   !      are untouched. Stated as an invariant over the group field rather than a count.  -----!
   cfg%output%grp_on(1) = .false.
   call build_output_registry(reg, cfg)
   call check(group_absent(reg, GRP_STRUCTURE), 'structure off -> no structure var in any tier')
   call check(in_tier(reg, 'gpp_site', 4_ik),   'structure off leaves carbon alone')
   cfg%output%grp_on(1) = .true.

   !----- Per-variable override: OFF removes it from every tier. -----!
   call build_output_registry(reg, cfg)
   call apply_variable_override(reg, 'growth_avg_cohort', OVR_FALSE, 0_ik, found)
   call check(found, 'override found growth_avg_cohort')
   call build_freq_index(reg)
   call check(.not. in_tier(reg, 'growth_avg_cohort', 3_ik), 'growth_avg_cohort OFF -> absent from MONTHLY')

   !----- Per-variable override: MASK re-targets streams (basal_area_cohort M -> D M). -----!
   call build_output_registry(reg, cfg)
   call check(.not. in_tier(reg, 'basal_area_cohort', 2_ik), 'basal_area_cohort default NOT in DAILY')
   call parse_stream_mask('D M', mask, status, bad)
   call check(status == 0_ik, 'parse_stream_mask("D M") ok')
   call apply_variable_override(reg, 'basal_area_cohort', OVR_MASK, mask, found)
   call build_freq_index(reg)
   call check(in_tier(reg, 'basal_area_cohort', 2_ik), 'basal_area_cohort override -> now in DAILY')

   !----- Unknown-key trap: a typo matches no registry variable. -----!
   call apply_variable_override(reg, 'grwoth_avg_cohort', OVR_FALSE, 0_ik, found)
   call check(.not. found, 'typo key -> not found (unknown-key trap)')

   !----- parse_stream_mask flags an unrecognized token. -----!
   call parse_stream_mask('D X', mask, status, bad)
   call check(status == 1_ik .and. trim(bad) == 'X', 'bad stream token flagged')

   write(*,'(a)') 'test_output_registry: ALL PASSED'

contains

   !----- The section 3.1 invariant: the annual stream is SITE-LEVEL ONLY. A cohort/patch window   !
   !      longer than a month would straddle the annual disturbance restructuring, so the slot set   !
   !      it averaged would not be the slot set present at flush.  ------------------------------!
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
