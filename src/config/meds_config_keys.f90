! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config_keys -- which keys a config may hold (N-10).                                    !
!                                                                                          !
! The two reference configs, meds_config_main.toml and meds_config_pft.toml, list every key   !
! MEDS reads, set or commented out at its default. A run checks each of its two files against !
! its reference: a key the reference does not list stops the run, because a key that parses   !
! and does nothing is worse than one that is absent -- `dt_forcing = 7200` (the key is         !
! `timestep`) used to run silently at the file's own spacing. A retired key is named with what  !
! replaced it; any other unknown key with the listed key it most likely meant.                  !
!                                                                                          !
! The lists come from the references at build time (cmake/meds_config_keys.cmake), and the     !
! test check_config_keys.cmake holds the references equal to the keys the loader reads.         !
!==========================================================================================!
module meds_config_keys
   use meds_kinds, only : ik
   use meds_toml,  only : toml_table_t
   implicit none
   private

   public :: key_report_t, check_config_keys, nearest_key, MAIN_KEYS, PFT_KEYS

   integer, parameter :: KEY_LEN = 64            !< meds_toml's key length
   include 'meds_config_keys.inc'                ! MAIN_KEYS, PFT_KEYS

   !----- Keys MEDS used to read. Each is refused with what replaced it, so a config written    !
   !      against an older MEDS says what to change rather than only that something is wrong. ---!
   character(len=*), parameter :: GROWTH_CURVE_GONE = 'growth is carbon-driven, by the NPP the '// &
                                                      'cohort allocates; the empirical growth curve is gone'
   character(len=*), parameter :: PHENO_SHARP = 'every phenology switch is sigma(s (x - x*)) with one '// &
      'signed sharpness: flush_temp_sharpness, shed_temp_sharpness, flush_daylength_sharpness, '//     &
      'shed_daylength_sharpness, flush_par_sharpness, shed_par_sharpness, flush_water_sharpness, '//   &
      'shed_water_sharpness'
   character(len=*), parameter :: PHENO_CHILL = 'chilling is gone: the warmth requirement is the '//   &
      'constant flush_degree_days, warmth above flush_base_temp summed from midwinter'
   character(len=*), parameter :: PHENO_COLD = 'senescence follows a cold sum below shed_base_temp '//  &
      'from midsummer (shed_degree_days) and short days (shed_daylength_threshold, '//             &
      'shed_daylength_sharpness)'
   character(len=*), parameter :: PHENO_WATER = 'the water cue is predawn leaf psi summed against the '// &
      'turgor-loss point: flush_water_sum and shed_water_sum [MPa day]; no soil-water cue'
   character(len=*), parameter :: PHENO_HABIT = 'leaf habit comes from the phenology parameters: '//    &
      'min_leaf_cover (the leaf cover senescence stops at) and the cue masks; background turnover '//  &
      'is not cold-suppressed'
   type :: retired_key_t
      character(len=KEY_LEN) :: key
      character(len=320)     :: why
   end type retired_key_t
   type(retired_key_t), parameter :: RETIRED(*) = [                                          &
      retired_key_t('soil_column.root_beta', 'the root profile is a plant trait: set '//          &
                    '[hydraulics].root_beta (0 < beta < 1) and root_depth; root_beta = '//        &
                    'exp(-b*root_depth) gives the old exponential decay b per metre'),           &
      retired_key_t('soil.dewmx', 'the canopy films are plant traits now: set pft.leaf_surf_water_max '// &
                    '[kg/m2 leaf] and pft.wood_surf_water_max [kg/m2 wood] in the PFT file; '//   &
                    'both default to 0.1, the old dewmx'),                                       &
      retired_key_t('energy.phase_change', 'ice-aware soil conductivity and heat capacity are '//  &
                    'always on now (the freeze/thaw plateau always was); delete the key'),       &
      retired_key_t('output.strict_caps', 'a run always stops at the step its live cohort or '//  &
                    'patch count exceeds output.cohort_max or patch_max; raise the cap'),        &
      retired_key_t('output.carbon_fluxes', 'renamed [output].carbon'),                          &
      retired_key_t('output.water_fluxes', 'renamed [output].water'),                            &
      retired_key_t('output.energy_fluxes', 'renamed [output].energy'),                          &
      retired_key_t('output.fast.interval_steps', 'this spelling was never read; the key is '//   &
                    '[output].fast_interval_steps'),                                             &
      retired_key_t('fast.ark_niter', 'use [fast].ark_coupled: true for the coupled Newton '//    &
                    'solve (ark_niter > 1), false for one pass'),                                &
      retired_key_t('site.utc_offset', 'every forcing file is in UTC (its time_zone attribute '//  &
                    'says so); convert a local-time source to UTC when you build the file'),     &
      retired_key_t('site.apply_solar_longitude', 'every forcing file is in UTC (its time_zone '// &
                    'attribute says so); convert a local-time source to UTC when you build it'), &
      retired_key_t('site.reference_height', 'the forcing is moved to each patch''s canopy-air '// &
                    'top; declare its own heights in [forcing]: tq_height, wind_height, '//       &
                    'height_above, wind_exposure'),                                              &
      retired_key_t('site.wind_meas_height', 'the forcing is moved to each patch''s canopy-air '// &
                    'top; declare its own heights in [forcing]: tq_height, wind_height, '//       &
                    'height_above, wind_exposure'),                                              &
      retired_key_t('site.apply_wind_profile', 'the forcing is moved to each patch''s canopy-air '//&
                    'top; declare its own heights in [forcing]: tq_height, wind_height, '//       &
                    'height_above, wind_exposure'),                                              &
      retired_key_t('site.wind_roughness_z0', 'the forcing is moved to each patch''s canopy-air '// &
                    'top; declare its own heights in [forcing]: tq_height, wind_height, '//       &
                    'height_above, wind_exposure'),                                              &
      retired_key_t('phenology.cue_mask', 'the flush and shed cues are selected independently: '// &
                    'flush_cue_mask and shed_cue_mask (TEMP 1, LIGHT 2, WATER 4); a cold-'//     &
                    'deciduous PFT is flush_cue_mask = [3], shed_cue_mask = [3]'),               &
      retired_key_t('phenology.k_flush_max', 'renamed flush_rate_max'),                          &
      retired_key_t('phenology.k_shed_max', 'renamed shed_rate_max'),                            &
      retired_key_t('phenology.tau_flush', 'renamed flush_cue_timescale'),                       &
      retired_key_t('phenology.tau_shed', 'renamed shed_cue_timescale'),                         &
      retired_key_t('phenology.gdd_base_temp', 'renamed flush_base_temp'),                       &
      retired_key_t('phenology.phen_a', PHENO_CHILL),                                            &
      retired_key_t('phenology.phen_b', PHENO_CHILL),                                            &
      retired_key_t('phenology.phen_c', PHENO_CHILL),                                            &
      retired_key_t('phenology.chill_base_temp', PHENO_CHILL),                                   &
      retired_key_t('phenology.cold_drop_daylength', PHENO_COLD),                                &
      retired_key_t('phenology.cold_drop_soiltemp1', PHENO_COLD),                                &
      retired_key_t('phenology.cold_drop_soiltemp2', PHENO_COLD),                                &
      retired_key_t('phenology.photo_crit', 'renamed flush_daylength_threshold'),                &
      retired_key_t('phenology.photo_slope', 'renamed flush_daylength_sharpness'),               &
      retired_key_t('phenology.cue_sharpness', PHENO_SHARP),                                     &
      retired_key_t('phenology.gdd_width', PHENO_SHARP),                                         &
      retired_key_t('phenology.daylen_width', PHENO_SHARP),                                      &
      retired_key_t('phenology.soiltemp_width', PHENO_SHARP),                                    &
      retired_key_t('phenology.light_width', PHENO_SHARP),                                       &
      retired_key_t('phenology.water_width', PHENO_WATER),                                       &
      retired_key_t('phenology.water_off_threshold', PHENO_WATER),                               &
      retired_key_t('phenology.water_on_threshold', PHENO_WATER),                                &
      retired_key_t('phenology.water_window', PHENO_WATER),                                      &
      retired_key_t('phenology.low_psi_threshold', 'renamed shed_water_sum, now a sum of psi '//   &
                    'below the turgor-loss point [MPa day]'),                                    &
      retired_key_t('phenology.high_psi_threshold', 'renamed flush_water_sum, now a sum of psi '// &
                    'above the turgor-loss point [MPa day]'),                                    &
      retired_key_t('phenology.light_on_threshold', 'renamed shed_par_threshold [umol/m2/s], '//  &
                    'with a positive shed_par_sharpness (the PAR cue, mask bit 8)'),             &
      retired_key_t('phenology.light_window', 'renamed par_window'),                             &
      retired_key_t('phenology.evg_ref_temp', PHENO_HABIT),                                      &
      retired_key_t('phenology.evg_slope', PHENO_HABIT),                                         &
      retired_key_t('pft.evergreen', PHENO_HABIT),                                               &
      retired_key_t('phenology.bare_snap_frac', 'renamed bare_leaf_cover'),                      &
      retired_key_t('phenology.phenology_on', 'leaf phenology always runs; the per-PFT cues are '//&
                    'the PFT file''s [phenology] block'),                                        &
      retired_key_t('carbon.growth_source', 'growth is always carbon-driven, by the NPP the '//    &
                    'cohort allocates'),                                                         &
      retired_key_t('pft.growth_dbh_slope', GROWTH_CURVE_GONE),                                  &
      retired_key_t('pft.growth_dbh_cap',   GROWTH_CURVE_GONE),                                  &
      retired_key_t('pft.growth_dbh_max',   GROWTH_CURVE_GONE),                                  &
      retired_key_t('pft.growth_lai_slope', GROWTH_CURVE_GONE) ]

   !----- What the check found, one line per key, reported with the missing required keys. -----!
   integer, parameter :: MAX_REPORT = 256
   type :: key_report_t
      integer(ik)        :: n = 0_ik
      character(len=480) :: line(MAX_REPORT) = ''
   end type key_report_t

contains

   !----- Check every key of `t` against `known`, the keys its reference lists (MAIN_KEYS or    !
   !      PFT_KEYS), adding a line to `report` for each retired or unknown one. ------------------!
   subroutine check_config_keys(t, known, report)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: known(:)
      type(key_report_t), intent(inout) :: report
      character(len=480)     :: line
      character(len=KEY_LEN) :: hint
      integer(ik) :: i, j
      do i = 1_ik, t%n
         if (any(known == t%key(i))) cycle
         line = trim(t%key(i))//' (in '//trim(t%source)//')'
         j = retired_index(t%key(i))
         if (j > 0_ik) then
            line = trim(line)//' is retired: '//trim(RETIRED(j)%why)
         else
            line = trim(line)//' is not a key MEDS reads'
            hint = nearest_key(t%key(i), known)
            if (len_trim(hint) > 0) line = trim(line)//'; did you mean '//trim(hint)//'?'
         end if
         if (report%n < MAX_REPORT) then
            report%n = report%n + 1_ik
            report%line(report%n) = line
         end if
      end do
   end subroutine check_config_keys

   pure integer(ik) function retired_index(key) result(j)
      character(len=*), intent(in) :: key
      do j = 1_ik, size(RETIRED, kind=ik)
         if (RETIRED(j)%key == key) return
      end do
      j = 0_ik
   end function retired_index

   !----- The listed key an unknown one most likely meant. First the same name under another    !
   !      section (a key that moved, like state.cohort_max to output.cohort_max), the closest     !
   !      in spelling if there are several; else the closest spelling within a third of the      !
   !      key's length (a typo). Blank when nothing is that close. ------------------------------!
   pure function nearest_key(key, known) result(best)
      character(len=*), intent(in) :: key, known(:)
      character(len=KEY_LEN) :: best
      integer :: i, d, d_best
      best = ''
      d_best = huge(1)
      do i = 1, size(known)
         if (leaf(known(i)) /= leaf(key)) cycle
         d = edit_distance(trim(key), trim(known(i)))
         if (d < d_best) then ; d_best = d ; best = known(i) ; end if
      end do
      if (len_trim(best) > 0) return
      d_best = max(2, len_trim(key) / 3) + 1
      do i = 1, size(known)
         d = edit_distance(trim(key), trim(known(i)))
         if (d < d_best) then ; d_best = d ; best = known(i) ; end if
      end do
   end function nearest_key

   !----- A key's name without its section: what follows the last dot. -------------------------!
   pure function leaf(key) result(name)
      character(len=*), intent(in) :: key
      character(len=len(key)) :: name
      name = key(index(trim(key), '.', back=.true.) + 1:)
   end function leaf

   !----- The number of single-character insertions, deletions and substitutions that turn a   !
   !      into b (Levenshtein), by the two-row recurrence. ------------------------------------!
   pure integer function edit_distance(a, b) result(d)
      character(len=*), intent(in) :: a, b
      integer :: prev(0:len(b)), cur(0:len(b)), i, j
      prev = [(j, j = 0, len(b))]
      do i = 1, len(a)
         cur(0) = i
         do j = 1, len(b)
            cur(j) = min(prev(j) + 1, cur(j-1) + 1, prev(j-1) + merge(0, 1, a(i:i) == b(j:j)))
         end do
         prev = cur
      end do
      d = prev(len(b))
   end function edit_distance

end module meds_config_keys
