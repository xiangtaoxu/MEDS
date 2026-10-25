# SPDX-License-Identifier: Apache-2.0
#----- meds_main refuses a run it cannot start, and names the cause. ------------------------------#
#
# Each case edits the reference config (meds_config_main.toml), or a copy of the PFT file it names
# (meds_config_pft.toml), into a run that must stop before its first step, runs meds_main on it, and
# asserts that it stopped and that its report names the cause:
#   config_unknown_keys        a misspelt key, a key that moved block, and a [region] block in a site run
#   config_rk45_loose_rtol     RK45 at a tolerance looser than it can integrate (fast.rtol > 1e-4)
#   config_malformed_value     a value the reader cannot parse as the number its key takes
#   hydraulics_pv_traits       a leaf elastic modulus below |pi0|, whose turgor-loss point is not
#                              negative, in a PFT file that sets its stress points by hand
#   hydraulics_vulnerability   a xylem psi50 above zero
#   soil_nan_layer             a NaN soil conductivity, as a map layer with no sand or clay gives; NaN
#                              fails every "<= 0" test, so the check must state what holds (Release
#                              builds only: a Debug build traps the NaN at that comparison)
#   toml_unclosed_string       a string value with no closing quote, which would otherwise read as ""
#   init_census_missing        a census file that is not there (init mode 1)
#   init_census_empty          a census file with no row for the site
#   init_restart_missing       a state file that is not there (init mode 2)
# Every case runs in its own directory of the build tree, on a one-day run whose input paths are
# absolute: a refusal that stops working then costs a day's run there, never a long run in the
# source tree.
#
# Usage: cmake -DMEDS_MAIN=<meds_main> -DSOURCE_DIR=<repo> -DWORK_DIR=<build-tree dir> -DCASE=<case>
#              -P check_main_refuses.cmake
#--------------------------------------------------------------------------------------------------#
cmake_minimum_required(VERSION 3.20)

#----- The edits the reference config must still accept; a key that moved stops the test here. ---#
function(replace_or_fail var pattern replacement what)
   string(REGEX MATCH "${pattern}" found "${${var}}")
   if(NOT found)
      message(FATAL_ERROR "the reference configuration no longer holds ${what}: the test cannot build its run")
   endif()
   string(REGEX REPLACE "${pattern}" "${replacement}" edited "${${var}}")
   set(${var} "${edited}" PARENT_SCOPE)
endfunction()

file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
file(READ "${SOURCE_DIR}/meds_config_main.toml" text)

#----- A one-day run, its input files named by absolute path. ------------------------------------#
replace_or_fail(text "\nstart_time *= *\"[^\"]*\"" "\nstart_time = \"2000-01-01\"" "start_time")
replace_or_fail(text "\nend_time *= *\"[^\"]*\"" "\nend_time = \"2000-01-02\"" "end_time")
replace_or_fail(text "\npft_config *= *\"([^/\"][^\"]*)\"" "\npft_config = \"${SOURCE_DIR}/\\1\"" "pft_config")
replace_or_fail(text "\ncensus_file *= *\"([^/\"][^\"]*)\"" "\ncensus_file = \"${SOURCE_DIR}/\\1\"" "census_file")

if(CASE STREQUAL "config_unknown_keys")
   string(APPEND text "\n[forcing]\ntimestpe = \"3600s\"\n[state]\ncohort_max = 2048\n[region]\nland_fraction_min = 0.5\n")
   set(expects "forcing.timestpe .*did you mean forcing.timestep\\?"
               "state.cohort_max .*did you mean output.cohort_max\\?"
               "\\[region\\] .*needs \\[run\\]\\.mode = \"region\""
               "the configuration has missing, unknown or inapplicable keys")
elseif(CASE STREQUAL "config_rk45_loose_rtol")
   #----- The fast loop on, RK45, rtol 1e-3: any line naming either key goes, then both are set. -#
   replace_or_fail(text "\nfast_biophysics_on *= *false" "\nfast_biophysics_on = true" "fast_biophysics_on = false")
   string(REGEX REPLACE "\n#? *time_integrator *=[^\n]*" "" text "${text}")
   string(REGEX REPLACE "\n#? *rtol *=[^\n]*" "" text "${text}")
   replace_or_fail(text "\n\\[fast\\][^\n]*\n" "\n[fast]\ntime_integrator = \"rk45\"\nrtol = 1.0e-3\n" "a [fast] block")
   set(expects "fast\\.rtol > 1e-4 with time_integrator = \"rk45\"")
elseif(CASE STREQUAL "config_malformed_value")
   replace_or_fail(text "\ndepth *= *[0-9.]+" "\ndepth = two" "[soil] depth")
   set(expects "cannot parse real for key \"soil\\.depth\"" "raw value: \"two\"" "malformed configuration value")
elseif(CASE STREQUAL "soil_nan_layer")
   replace_or_fail(text "\nsolid_conductivity *= *[0-9.eE+-]+" "\nsolid_conductivity = nan" "[soil] solid_conductivity")
   set(expects "soil conductivities must be positive \\(a map layer needs sand or clay\\)")
elseif(CASE STREQUAL "toml_unclosed_string")
   #----- The diagnostics on, and the file of their variables named with no closing quote. --------#
   replace_or_fail(text "\nenabled *= *false" "\nenabled = true" "[output] enabled")
   replace_or_fail(text "\nio_config *= *\"\"[^\n]*" "\nio_config = \"my_variables.txt" "[output] io_config")
   set(expects "cannot parse string[^\n]*for key \"output\\.io_config\"")
elseif(CASE MATCHES "^hydraulics_")
   #----- The PFT file's copy beside the run, with one trait made impossible for every PFT. --------#
   file(READ "${SOURCE_DIR}/meds_config_pft.toml" pft)
   if(CASE STREQUAL "hydraulics_pv_traits")
      #----- The stress points set by hand, as a PFT file may: left to their default, the positive    #
      #      turgor-loss point would already be refused through them. -------------------------------#
      foreach(key wstress_psi_open wstress_psi_close stomata_psi_onset)
         replace_or_fail(pft "\n# *${key}( *=)" "\n${key}\\1" "[pft] ${key}, commented out")
      endforeach()
      replace_or_fail(pft "\nleaf_elastic_mod *= *[^ #\n]+" "\nleaf_elastic_mod = 1.0" "[hydraulics] leaf_elastic_mod")
      set(expects "hydraulics\\.leaf_elastic_mod and wood_elastic_mod must exceed \\|pi0\\|")
   else()
      replace_or_fail(pft "\nwood_psi50 *= *[^ #\n]+" "\nwood_psi50 = 0.5" "[hydraulics] wood_psi50")
      set(expects "hydraulics\\.wood_psi50 must be < 0")
   endif()
   file(WRITE "${WORK_DIR}/meds_config_pft.toml" "${pft}")
   replace_or_fail(text "\npft_config *= *\"[^\"]*\"" "\npft_config = \"${WORK_DIR}/meds_config_pft.toml\"" "pft_config")
elseif(CASE MATCHES "^init_")
   if(CASE STREQUAL "init_restart_missing")
      set(mode 2)
      set(key "restart_file")
      set(file "${WORK_DIR}/no_such_state.nc")
   else()
      set(mode 1)
      set(key "census_file")
      if(CASE STREQUAL "init_census_empty")
         set(file "${WORK_DIR}/empty_census.csv")
         file(WRITE "${file}" "site_id,patch_id,cohort_id,dbh,height,pft,nplant\n")
      else()
         set(file "${WORK_DIR}/no_such_census.csv")
      endif()
   endif()
   string(REGEX REPLACE "([][+.*()^$?|\\\\])" "\\\\\\1" file_re "${file}")   # the path, matched literally
   if(CASE STREQUAL "init_restart_missing")
      set(expects "restart \\(mode 2\\): the state file ${file_re} cannot be opened")
   elseif(CASE STREQUAL "init_census_empty")
      set(expects "census \\(mode 1\\): the census file ${file_re} holds no cohort for the site")
   else()
      set(expects "census \\(mode 1\\): the census file ${file_re} cannot be opened")
   endif()
   string(REGEX REPLACE "\n${key} *= *\"[^\"]*\"[^\n]*" "" text "${text}")
   replace_or_fail(text "\ninit_mode *= *[0-9]" "\ninit_mode = ${mode}\n${key} = \"${file}\"" "init_mode")
else()
   message(FATAL_ERROR "check_main_refuses: unknown CASE '${CASE}'")
endif()

file(WRITE "${WORK_DIR}/meds_config_main.toml" "${text}")
execute_process(COMMAND "${MEDS_MAIN}" "${WORK_DIR}/meds_config_main.toml" WORKING_DIRECTORY "${WORK_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
set(problems "")
if(rc EQUAL 0)
   string(APPEND problems "\n  the run did not stop")
endif()
if(CASE MATCHES "^init_" AND log MATCHES "bare ground")
   string(APPEND problems "\n  the run started from bare ground")
endif()
foreach(expect IN LISTS expects)
   if(NOT log MATCHES "${expect}")
      string(APPEND problems "\n  the report lacks: ${expect}")
   endif()
endforeach()
if(problems)
   message(FATAL_ERROR "meds_main did not refuse case ${CASE} as expected:${problems}\n--- output:\n${log}")
endif()
message(STATUS "meds_main refused case ${CASE}, naming the cause")
