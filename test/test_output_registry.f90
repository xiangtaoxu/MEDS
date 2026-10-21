! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_output_registry -- the registration list + the selection: the defaults, the io_config  !
! lines (names and patterns, the most specific winning), the axes' bookkeeping, the tier     !
! switches, the unknown-key trap (docs/science/diagnostics.md §5).                           !
!==========================================================================================!
program test_output_registry
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_output_params,   only : FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE, N_FREQ
   use meds_output_types,    only : output_registry_t, MAX_OUTPUT_VARS, DIM_COHORT, DIM_PATCH
   use meds_output_registry, only : build_output_registry, build_freq_index, find_var_index,                      &
                                    apply_variable_override, parse_stream_mask, freq_bit, selection_t,            &
                                    apply_selections, resolve_streams, glob_match, OVR_TRUE, OVR_FALSE, OVR_MASK
   use meds_test_support, only : banner, build_test_config, check
   implicit none

   type(meds_config_t)     :: cfg
   type(output_registry_t) :: reg
   integer(ik) :: mask, status, t, k
   character(len=8) :: bad
   logical :: found

   call banner('output_registry')
   cfg = build_test_config(86400.0_wp)     ! 1-day slow step; cfg%output holds the [output] defaults

   !----- Build with the defaults (freq_on=[F,T,T,T]). -----!
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
   !----- The defaults (docs/science/diagnostics.md §5): the radiation and ecophys groups, the    !
   !      soil-by-patch columns and the snapshots are off until a selection names them; the        !
   !      axes' bookkeeping is written with the default cohort and patch variables. ---------------!
   call check(in_tier(reg, 'mean_soil_temp_layer_si', 3_ik), 'energy on by default -> mean_soil_temp_layer_si in MONTHLY')
   call check(.not. live(reg, 'mean_gsw_co'), 'the ecophys group off by default')
   call check(.not. live(reg, 'mean_abs_sw_co'), 'the radiation group off by default')
   call check(.not. live(reg, 'mean_soil_water_layer_pa'), 'the soil by patch off by default')
   call check(.not. live(reg, 'agb_py'), 'a snapshot off by default')
   call check(in_tier(reg, 'global_cohort_id', 2_ik) .and. in_tier(reg, 'global_patch_id', 3_ik),         &
              'the bookkeeping with the default cohort and patch variables')

   !----- FAST tier off by default. The ANNUAL tier must contain NO cohort/patch-dimensioned      !
   !      variable (docs/science/diagnostics.md §3) -- that invariant is what the check below     !
   !      states, rather than a count that has to be re-derived on every registry edit. ----------!
   call check(reg%nidx(1) == 0_ik, 'FAST tier empty by default (freq_on(FAST)=false)')
   call check(reg%nidx(2) > 0_ik,  'DAILY tier populated')
   call check(reg%nidx(4) > 0_ik,  'ANNUAL tier populated')
   call check(annual_is_site_level(reg), 'ANNUAL tier carries no cohort/patch variable')
   !----- A pattern: "*_co" = false drops the cohort axis, its bookkeeping included, and nothing  !
   !      else. -----------------------------------------------------------------------------------!
   call choose([selection_t('*_co', OVR_FALSE, FREQ_NONE)])
   call check(.not. live(reg, 'mean_agb_co'), '"*_co" = false -> mean_agb_co in no tier')
   call check(.not. live(reg, 'global_cohort_id'), '"*_co" = false -> no cohort bookkeeping')
   call check(in_tier(reg, 'mean_agb_py', 2_ik), '"*_co" = false leaves the polygon variables alone')
   call check(in_tier(reg, 'global_patch_id', 2_ik), '"*_co" = false leaves the patch axis alone')
   !----- The most specific line wins: a name over a pattern, more literal characters over fewer. !
   call choose([selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('mean_*_py', OVR_MASK, FREQ_MONTHLY),      &
                selection_t('mean_gpp_py', OVR_MASK, ior(FREQ_DAILY, FREQ_ANNUAL))])
   call check(in_tier(reg, 'mean_agb_py', 3_ik) .and. .not. in_tier(reg, 'mean_agb_py', 2_ik),              &
              '"mean_*_py" wins over "*"')
   call check(in_tier(reg, 'mean_gpp_py', 2_ik) .and. in_tier(reg, 'mean_gpp_py', 4_ik) .and.               &
              .not. in_tier(reg, 'mean_gpp_py', 3_ik), 'a name wins over a pattern')
   call check(.not. live(reg, 'nplant_co') .and. .not. live(reg, 'owner_patch'),                             &
              '"*" = false: no other variable, so no bookkeeping')
   !----- A pattern gives each variable the tiers it asks that the variable can be on. -----!
   call choose([selection_t('mean_agb_*', OVR_MASK, ior(FREQ_MONTHLY, FREQ_ANNUAL))])
   call check(in_tier(reg, 'mean_agb_py', 4_ik) .and. in_tier(reg, 'mean_agb_co', 3_ik) .and.               &
              .not. in_tier(reg, 'mean_agb_co', 4_ik), 'a pattern gives Y only off the cohort axis')
   call choose([selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('mean_*_co', OVR_MASK, FREQ_FAST)])
   call check(.not. live(reg, 'mean_gsw_co'), 'the F tier switched off writes nothing')
   !----- The bookkeeping follows its axis, tier by tier. -----!
   call choose([selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('mean_agb_co', OVR_MASK, FREQ_MONTHLY)])
   call check(in_tier(reg, 'owner_patch', 3_ik) .and. .not. in_tier(reg, 'owner_patch', 2_ik),              &
              'the cohort bookkeeping at the cohort variables'' tiers only')
   call check(.not. live(reg, 'cohort_offset'), 'no patch variable, no patch bookkeeping')
   call choose([selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('mean_soil_water_layer_pa', OVR_MASK,     &
                FREQ_DAILY)])
   call check(in_tier(reg, 'cohort_offset', 2_ik), 'the soil by patch is on the patch axis')
   !----- `true`: a variable's default; a pattern's true leaves one off by default off. -----!
   call choose([selection_t('*', OVR_FALSE, FREQ_NONE), selection_t('*_py', OVR_TRUE, FREQ_NONE)])
   call check(in_tier(reg, 'mean_agb_py', 4_ik) .and. .not. live(reg, 'agb_py'),                            &
              'a pattern''s true restores each default, off included')
   !----- The patterns. -----!
   call check(glob_match('*_py', 'mean_gpp_py') .and. glob_match('mean_*_py', 'mean_gpp_pft_py') .and.      &
              glob_match('*', 'x') .and. glob_match('a*b*c', 'axxbyyc') .and. glob_match('a*b', 'ab'),       &
              'glob_match: a * matches any run of characters, none included')
   call check(.not. glob_match('mean_*_co', 'mean_gpp_py') .and. .not. glob_match('a*b', 'abc') .and.       &
              .not. glob_match('gpp_py', 'mean_gpp_py'), 'glob_match: the literal parts must match')
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

   !----- Per-tier enable: turn DAILY off -> DAILY empties, MONTHLY unaffected, and a selection   !
   !      cannot put a variable back on it. -----!
   cfg%output%freq_on(2) = .false.
   call build_output_registry(reg, cfg)
   call check(reg%nidx(2) == 0_ik, 'DAILY off -> empty')
   call check(reg%nidx(3) > 0_ik,  'MONTHLY still populated')
   call choose([selection_t('mean_agb_py', OVR_MASK, ior(FREQ_DAILY, FREQ_MONTHLY))])
   call check(reg%nidx(2) == 0_ik, 'DAILY off -> a selection writes nothing there')
   cfg%output%freq_on(2) = .true.

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

   !----- The registry from the config, with these io_config lines applied. -------------------!
   subroutine choose(sel)
      type(selection_t), intent(in) :: sel(:)
      call build_output_registry(reg, cfg)
      call apply_selections(reg, sel)
      call resolve_streams(reg)
   end subroutine choose

   !----- .true. if the variable is live in some tier. ---------------------------------------!
   logical function live(reg, name) result(yes)
      type(output_registry_t), intent(in) :: reg
      character(len=*),        intent(in) :: name
      integer(ik) :: t
      yes = .false.
      do t = 1_ik, N_FREQ
         if (in_tier(reg, name, t)) yes = .true.
      end do
   end function live

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
