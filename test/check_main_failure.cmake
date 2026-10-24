# SPDX-License-Identifier: Apache-2.0
# Runs meds_main on a three-day run whose first step fails, and fails unless meds_main stops with an
# error after finalizing the run. The failure: a soil-carbon spin-up fed so much litter that its
# steady pools pass the plausibility ceiling, which the step's bounds guard stops. The end-of-run
# report in the log shows that region_finalize ran, and region_finalize is what closes the output
# files after a failure (test_region reads them).
#
# The run is the reference config under a prepended override block (meds_toml returns the first
# match of a key), and it runs in WORK_DIR, never in the source tree.
#
# Usage: cmake -DMEDS_MAIN=<meds_main> -DSOURCE_DIR=<repo> -DWORK_DIR=<scratch dir>
#              -P check_main_failure.cmake
cmake_minimum_required(VERSION 3.20)
file(REMOVE_RECURSE "${WORK_DIR}")
file(MAKE_DIRECTORY "${WORK_DIR}")
set(overrides "[run]
start_time = \"2000-01-01\"
end_time   = \"2000-01-04\"
[init]
init_mode  = 0
pft_config = \"${SOURCE_DIR}/meds_config_pft.toml\"
[soil_carbon]
soil_carbon_on     = true
spinup_steady      = true
spinup_labile_grnd = 1000.0
[output]
enabled = true
dir     = \"${WORK_DIR}/out\"
")
file(READ "${SOURCE_DIR}/meds_config_main.toml" text)
file(WRITE "${WORK_DIR}/main.toml" "${overrides}${text}")
execute_process(COMMAND "${MEDS_MAIN}" "${WORK_DIR}/main.toml" WORKING_DIRECTORY "${WORK_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
set(problems "")
if(rc EQUAL 0)
   string(APPEND problems "\n  the run did not stop with an error")
endif()
foreach(expect "a cohort or soil-carbon state out of bounds" "site area start=")
   string(FIND "${log}" "${expect}" at)
   if(at EQUAL -1)
      string(APPEND problems "\n  the log lacks: ${expect}")
   endif()
endforeach()
if(problems)
   message(FATAL_ERROR "meds_main did not finalize the run before stopping on a failed step:${problems}\n--- output:\n${log}")
endif()
message(STATUS "meds_main finalized the run, then stopped on the failed step")
