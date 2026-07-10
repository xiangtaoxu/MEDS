!==========================================================================================!
! test_output_registry -- the registration list + the §6.1 config resolution (group toggles,     !
! per-tier enables, per-variable overrides, unknown-key trap). MEDS_IO_DESIGN.md test 4.          !
!==========================================================================================!
program test_output_registry
   use meds_kinds,           only : wp, ik
   use meds_config,          only : meds_config_t
   use meds_output_config,   only : FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, N_FREQ, N_GRP
   use meds_output_types,    only : output_registry_t
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

   !----- Build with the defaults (grp_on=[T,T,F,F], freq_on=[F,T,T,T]). -----!
   call build_output_registry(reg, cfg)
   call check(reg%nvar == 20_ik, 'P0 registry has 20 variables (n_cohort/n_patch are structural)')
   call check(find_var_index(reg, 'agb_cohort') > 0_ik, 'agb_cohort present')
   call check(find_var_index(reg, 'agb_site')   > 0_ik, 'agb_site present')
   call check(find_var_index(reg, 'nope')       == 0_ik, 'unknown name -> 0')

   !----- FAST tier off by default -> no variables; annual = the 4 site scalars only. -----!
   call check(reg%nidx(1) == 0_ik, 'FAST tier empty by default (freq_on(FAST)=false)')
   call check(reg%nidx(2) > 0_ik,  'DAILY tier populated')
   call check(reg%nidx(4) == 4_ik, 'ANNUAL tier = 4 site scalars only')
   call check(in_tier(reg, 'agb_cohort', 2_ik),        'agb_cohort in DAILY')
   call check(.not. in_tier(reg, 'agb_cohort', 4_ik),  'agb_cohort NOT in ANNUAL (no cohort dim annually)')
   call check(in_tier(reg, 'agb_site', 4_ik),          'agb_site in ANNUAL')

   !----- Per-tier enable: turn DAILY off -> DAILY empties, MONTHLY unaffected. -----!
   cfg%output%freq_on(2) = .false.
   call build_output_registry(reg, cfg)
   call check(reg%nidx(2) == 0_ik, 'DAILY off -> empty')
   call check(reg%nidx(3) > 0_ik,  'MONTHLY still populated')
   cfg%output%freq_on(2) = .true.

   !----- Group toggle: structure OFF -> every P0 var vanishes (all are GRP_STRUCTURE). -----!
   cfg%output%grp_on(1) = .false.
   call build_output_registry(reg, cfg)
   call check(all(reg%nidx(1:N_FREQ) == 0_ik), 'structure group off -> all tiers empty')
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
