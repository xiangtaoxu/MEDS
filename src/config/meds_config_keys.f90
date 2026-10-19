! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_config_keys -- which keys a config may hold (docs/configuration.md).                !
!                                                                                          !
! The two reference configs, meds_config_main.toml and meds_config_pft.toml, list every    !
! key MEDS reads, set or commented out at its default. A run checks each of its two files  !
! against its reference: a key the reference does not list stops the run, because a key    !
! that parses and does nothing is worse than one that is absent -- `dt_forcing = 7200`     !
! (the key is `timestep`) would otherwise run silently at the file's own spacing. A        !
! retired key is named with what replaced it; any other unknown key with the listed key it !
! most likely meant.                                                                       !
!                                                                                          !
! The lists come from the references at build time (cmake/meds_config_keys.cmake), and the !
! test check_config_keys.cmake holds the references equal to the keys the loader reads.    !
!==========================================================================================!
module meds_config_keys
   use meds_kinds, only : ik
   use meds_toml,  only : toml_table_t
   implicit none
   private

   public :: key_report_t, check_config_keys, nearest_key, MAIN_KEYS, PFT_KEYS, SITE_GROUND_KEYS

   integer, parameter :: KEY_LEN = 64            !< meds_toml's key length
   include 'meds_config_keys.inc'                ! MAIN_KEYS, PFT_KEYS

   !----- Retired keys. Each is refused with what replaced it, so a config written against       !
   !      an older MEDS says what to change rather than only that something is wrong. -----------!
   character(len=*), parameter :: GROWTH_CURVE_GONE = 'growth is carbon-driven, by the NPP the '// &
                                                      'cohort allocates; the empirical growth curve is gone'
   character(len=*), parameter :: PHENO_SHARP = 'every phenology switch is sigma(s (x - x*)) with one '// &
      'signed sharpness: flush_temp_sharpness, shed_temp_sharpness, flush_light_sharpness, '//         &
      'shed_light_sharpness, flush_water_sharpness, shed_water_sharpness'
   character(len=*), parameter :: PHENO_CHILL = 'chilling is gone: the warmth requirement is the '//   &
      'constant flush_degree_days, warmth above flush_base_temp summed from midwinter'
   character(len=*), parameter :: PHENO_COLD = 'senescence follows a cold sum below shed_base_temp '//  &
      'from midsummer (shed_degree_days) and short days (shed_light_hours, '//                     &
      'shed_light_sharpness, with a low par_min)'
   character(len=*), parameter :: PHENO_WATER = 'the water cue is predawn leaf psi summed against the '// &
      'turgor-loss point: flush_water_sum and shed_water_sum [MPa day]; no soil-water cue'
   character(len=*), parameter :: PHENO_HABIT = 'leaf habit comes from the phenology parameters: '//    &
      'min_leaf_cover (the leaf cover senescence stops at) and the cue masks; background turnover '//  &
      'is not cold-suppressed'
   character(len=*), parameter :: SOIL_NETWORK = 'soil water is solved with the plants'' water as '// &
      'one network (implicit Euler, Newton, its own sub-steps), whose error is set by '//         &
      '[fast].water_rtol and [fast].water_atol; the Richards solver''s linearization and sub-step '// &
      'controls are gone'
   !----- [soil] is the one block for the ground: the column, its surface and both bottom        !
   !      boundaries. Solver settings are in [fast]. -------------------------------------------------!
   character(len=*), parameter :: GROUND = ' ([soil] holds the whole ground, [fast] the solver '// &
      'settings)'
   character(len=*), parameter :: NITROGEN = 'the nitrogen cycle is not built; it is planned after '// &
      'the move to ED-MORE, and MEDS is carbon-only until then: delete the key'
   character(len=*), parameter :: LITTER_SPLIT = 'the litter is split by the tissue it comes from, with '// &
      'each PFT''s f_labile_leaf, f_labile_stem and aboveground_frac; delete the key'
   character(len=*), parameter :: CROWN_AREA = 'every crown closes over its patch (crown area 1, '// &
      'ED2''s closed-crown path), so the crown-area allometry and the crown fraction are gone; the '// &
      'in-canopy wind attenuates by leaf area alone: delete the key'
   type :: retired_key_t
      character(len=KEY_LEN) :: key
      character(len=320)     :: why
   end type retired_key_t
   type(retired_key_t), parameter :: RETIRED(*) = [                                          &
      retired_key_t('soil_column.root_beta', 'the root profile is a plant trait: set '//          &
                    '[hydraulics].root_beta (0 < beta < 1), and the rooting depth through '//    &
                    'root_depth_b1 and root_depth_b2'),                                          &
      retired_key_t('soil_column.n_layer', 'the soil layers are SoilGrids'' depths (5, 15, 30, '// &
                    '60, 100, 150, 200 cm, then deeper) cut at soil.depth; set '//               &
                    'soil.layer_bottoms [m] for other layers'),                                  &
      retired_key_t('soil_column.grid_growth', 'the layers no longer thicken by a growth '//     &
                    'factor: see soil_column.n_layer'),                                          &
      retired_key_t('soil_column.depth', 'moved to soil.depth'//GROUND),                         &
      retired_key_t('soil_column.layer_bottoms', 'moved to soil.layer_bottoms'//GROUND),         &
      retired_key_t('soil_column.retention', 'moved to soil.retention'//GROUND),                 &
      retired_key_t('soil_column.theta_sat', 'moved to soil.theta_sat'//GROUND),                 &
      retired_key_t('soil_column.theta_res', 'moved to soil.theta_res'//GROUND),                 &
      retired_key_t('soil_column.ksat', 'moved to soil.ksat'//GROUND),                           &
      retired_key_t('soil_column.curve_par_a', 'moved to soil.curve_par_a'//GROUND),             &
      retired_key_t('soil_column.curve_par_n', 'moved to soil.curve_par_n'//GROUND),             &
      retired_key_t('soil_column.psi_fc', 'moved to soil.psi_fc'//GROUND),                       &
      retired_key_t('soil_column.solid_conductivity', 'moved to soil.solid_conductivity'//GROUND), &
      retired_key_t('soil_column.dry_conductivity', 'moved to soil.dry_conductivity'//GROUND),   &
      retired_key_t('soil_column.dry_heat_capacity', 'moved to soil.dry_heat_capacity'//GROUND), &
      retired_key_t('soil.bottom_bc', 'renamed soil.water_bottom_bc (the heat''s is '//          &
                    'soil.thermal_bottom_bc)'//GROUND),                                        &
      retired_key_t('soil.rtol', 'moved to fast.water_rtol'//GROUND),                            &
      retired_key_t('soil.atol', 'moved to fast.water_atol'//GROUND),                            &
      retired_key_t('soil.intercept_alpha', 'moved to fast.intercept_alpha, beside canopy_water_on'), &
      retired_key_t('soil.intercept_k', 'moved to fast.intercept_k, beside canopy_water_on'),    &
      retired_key_t('energy.bottom_bc', 'renamed soil.thermal_bottom_bc'//GROUND),               &
      retired_key_t('energy.deep_temp', 'moved to soil.deep_temp'//GROUND),                      &
      retired_key_t('energy.deep_depth', 'moved to soil.deep_depth'//GROUND),                    &
      retired_key_t('energy.atol', 'renamed fast.soil_energy_budget_atol'//GROUND),              &
      retired_key_t('energy.debug_error', 'moved to fast.debug_error'//GROUND),                  &
      retired_key_t('hydraulics.root_depth', 'each cohort roots to root_depth_b1 * height^'//     &
                    'root_depth_b2, capped at the soil column: set [hydraulics].root_depth_b1 '// &
                    'and root_depth_b2 (defaults 1.114 and 0.4223, ED2), or pft.root_depth_b1/b2'), &
      retired_key_t('hydraulics.specific_root_area', 'the soil->root conductance takes root '//   &
                    'length: set [hydraulics].specific_root_length [m/kgC] and fine_root_radius '// &
                    '[m], or their pft.* values'),                                               &
      retired_key_t('hydraulics.k_plant_max', 'the plant''s conductance is its sapwood''s, '//      &
                    'wood_kmax * sapwood area / (height * vessel_curl), with wood_kmax from the '// &
                    'PFT''s wood density (Xu et al. 2016); set pft.wood_kmax to override it'),   &
      retired_key_t('pft.k_plant_max', 'the plant''s conductance is its sapwood''s, wood_kmax '//    &
                    '* sapwood area / (height * vessel_curl): set pft.wood_kmax [kg/m/s/MPa], '// &
                    'which defaults to the PFT''s wood density (Xu et al. 2016)'),             &
      retired_key_t('hydraulics.conductance', 'the sapwood (segment) conductance is the only '//   &
                    'form: wood_kmax * sapwood area / (height * vessel_curl)'),                  &
      retired_key_t('hydraulics.wood_kmax', 'wood_kmax follows each PFT''s wood density (Xu et '// &
                    'al. 2016: 3.3 kg/m/s/MPa at 0.6 g/cm3); set pft.wood_kmax to override it'), &
      retired_key_t('leaf_physiology.wstress_nonstomatal', 'the leaf-potential limb on '//      &
                    'Vcmax/Jmax/TPU is always on: it is the fast water-stress control, from the '// &
                    'turgor-loss point to twice it (pft.wstress_psi_open/_close override)'),     &
      retired_key_t('leaf_physiology.low_water_potential_control', 'the linear decline of gs '// &
                    'on the daily predawn potential is gone: the always-on leaf-potential limb '// &
                    'on Vcmax/Jmax/TPU closes the stomata within the day instead'),              &
      retired_key_t('soil.psi_wilt', 'there is no soil-side wilting point: the plant hydraulics '// &
                    'limits root uptake (the soil->root conductance falls with K(theta), and a '//   &
                    'layer drier than the wood gives nothing), so a wilting ramp cut the same '//    &
                    'water twice'),                                                              &
      retired_key_t('soil.psi_open', 'there is no soil-side wilting ramp: the plant hydraulics '// &
                    'limits root uptake (see soil.psi_wilt)'),                                   &
      retired_key_t('soil.linearize', SOIL_NETWORK),                                             &
      retired_key_t('soil.max_picard', SOIL_NETWORK),                                            &
      retired_key_t('soil.substep', SOIL_NETWORK),                                               &
      retired_key_t('soil.h_init', SOIL_NETWORK),                                                &
      retired_key_t('soil.max_substep', SOIL_NETWORK),                                           &
      retired_key_t('soil.dewmx', 'the canopy films are plant traits now: set pft.leaf_surf_water_max '// &
                    '[kg/m2 leaf] and pft.wood_surf_water_max [kg/m2 wood] in the PFT file; '//   &
                    'both default to 0.1, the old dewmx'),                                       &
      retired_key_t('energy.phase_change', 'ice-aware soil conductivity and heat capacity are '//  &
                    'always on now (the freeze/thaw plateau always was); delete the key'),       &
      retired_key_t('output.strict_caps', 'a run always stops at the step its live cohort or '//  &
                    'patch count exceeds output.cohort_max or patch_max; raise the cap'),        &
      retired_key_t('pft.leaf_lifespan_toc', 'leaves are parameterized by their turnover rate: '// &
                    'set pft.leaf_turnover_rate_toc [1/yr] = 1 / the lifespan [yr] (a 2-yr '//      &
                    'lifespan is 0.5)'),                                                         &
      retired_key_t('derived.kplastic_llspan', 'the leaf turnover rate is plastic in its place: '// &
                    'set derived.kplastic_leaf_turnover = -kplastic_llspan (<= 0: rate = '//       &
                    'leaf_turnover_rate_toc * exp(k * LAI above), slower in shade)'),             &
      retired_key_t('trait_dynamics.fineroot_lifespan_plastic', 'renamed '//                  &
                    'trait_dynamics.fineroot_turnover_plastic: fine-root turnover x (the '//       &
                    'cohort''s leaf turnover rate / the top-of-canopy one)'),                    &
      retired_key_t('pft.repro_carbon_efficiency', 'renamed pft.recruit_carbon_efficiency: the '// &
                    'share of reproduction carbon that becomes recruits (seed set, '//         &
                    'germination and seedling survival together)'),                            &
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
      retired_key_t('phenology.photo_crit', 'the light cue counts hours of light: use '//         &
                    'flush_light_hours [h/day] with a low par_min (then ~ the day length)'),     &
      retired_key_t('phenology.photo_slope', 'renamed flush_light_sharpness [1/h]'),             &
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
      retired_key_t('phenology.light_on_threshold', 'the light cue counts hours of light: '//    &
                    'set par_min to that PAR [umol/m2/s] and shed_light_hours with a positive '// &
                    'shed_light_sharpness (LIGHT, mask bit 2)'),                                 &
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
      retired_key_t('pft.growth_lai_slope', GROWTH_CURVE_GONE),                                  &
      retired_key_t('soil_carbon.agf_fast', LITTER_SPLIT),                                      &
      retired_key_t('soil_carbon.agf_struct', LITTER_SPLIT),                                    &
      retired_key_t('soil_carbon.n_cycle_on', NITROGEN),                                        &
      retired_key_t('soil_carbon.c2n_structural', NITROGEN),                                    &
      retired_key_t('soil_carbon.c2n_slow', NITROGEN),                                          &
      retired_key_t('soil_carbon.c2n_fast', NITROGEN),                                          &
      retired_key_t('soil_carbon.n_immobil_supply_scale', NITROGEN),                            &
      retired_key_t('allometry.ca_b1', CROWN_AREA),                                             &
      retired_key_t('allometry.ca_b2', CROWN_AREA),                                             &
      retired_key_t('pft.crown_area_frac', CROWN_AREA),                                         &
      retired_key_t('region.detail_polygons', 'the region files hold every variable and tier for '// &
                    'every polygon, cohort, patch and sub-daily output included; choose the '//      &
                    'variables with [output].io_config') ]

   !----- The [soil] keys a [[sites]] table may set, under the same names: the site's soil column !
   !      and its bottom. The ground's surface and its thermal boundary are the run's. -----------!
   character(len=*), parameter :: SITE_GROUND_KEYS(15) = [character(len=24) ::                     &
      'depth', 'layer_bottoms', 'retention', 'theta_sat', 'theta_res', 'ksat', 'curve_par_a',        &
      'curve_par_n', 'psi_fc', 'solid_conductivity', 'dry_conductivity', 'dry_heat_capacity',         &
      'water_bottom_bc', 'deep_temp', 'deep_depth']

   !----- What the check found, one line per key, reported with the missing required keys. -----!
   integer, parameter :: MAX_REPORT = 256
   type :: key_report_t
      integer(ik)        :: n = 0_ik
      character(len=480) :: line(MAX_REPORT) = ''
   end type key_report_t

contains

   !----- Check every key of `t` against `known`, the keys its reference lists (MAIN_KEYS or    !
   !      PFT_KEYS), adding a line to `report` for each retired or unknown one. A key of an     !
   !      array's table, `sites.2.depth`, is checked as `sites.depth`; a [[sites]] table also   !
   !      takes SITE_GROUND_KEYS. --------------------------------------------------------------!
   subroutine check_config_keys(t, known, report)
      type(toml_table_t), intent(in)    :: t
      character(len=*),   intent(in)    :: known(:)
      type(key_report_t), intent(inout) :: report
      character(len=480)     :: line
      character(len=KEY_LEN) :: hint, key
      integer(ik) :: i, j
      do i = 1_ik, t%n
         key = table_key(t%key(i))
         if (any(known == key)) cycle
         if (key(1:6) == 'sites.' .and. any(SITE_GROUND_KEYS == key(7:))) cycle
         line = trim(t%key(i))//' (in '//trim(t%source)//')'
         j = retired_index(key)
         if (j > 0_ik) then
            line = trim(line)//' is retired: '//trim(RETIRED(j)%why)
         else
            line = trim(line)//' is not a key MEDS reads'
            hint = nearest_key(key, known)
            if (len_trim(hint) > 0) line = trim(line)//'; did you mean '//trim(hint)//'?'
         end if
         if (report%n < MAX_REPORT) then
            report%n = report%n + 1_ik
            report%line(report%n) = line
         end if
      end do
   end subroutine check_config_keys

   !----- A key without its table's number: `sites.2.depth` is `sites.depth`. ------------------!
   pure function table_key(key) result(k)
      character(len=*), intent(in) :: key
      character(len=KEY_LEN) :: k
      integer :: d1, d2
      k = key
      d1 = index(key, '.')
      if (d1 == 0) return
      d2 = index(key(d1 + 1:), '.')
      if (d2 <= 1) return
      if (verify(key(d1 + 1:d1 + d2 - 1), '0123456789') == 0) k = key(1:d1)//key(d1 + d2 + 1:)
   end function table_key

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
