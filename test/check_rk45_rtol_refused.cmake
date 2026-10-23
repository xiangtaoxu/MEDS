# SPDX-License-Identifier: Apache-2.0
# Runs meds_main on the reference config with the fast loop on and switched to the explicit integrator at
# an rtol looser than it accepts, and fails unless the run stops before simulating and says why. (With the
# fast loop off the [fast] keys are not read, so the reference config as shipped would not be refused.)
#
# Usage: cmake -DMEDS_MAIN=<meds_main> -DSOURCE_DIR=<repo> -DWORK_FILE=<scratch toml>
#              -P check_rk45_rtol_refused.cmake
cmake_minimum_required(VERSION 3.20)
file(READ "${SOURCE_DIR}/meds_config_main.toml" text)
foreach(line "fast_biophysics_on  = false" "# time_integrator     = \"esdirk2\"" "# rtol                = 2.0e-3")
   string(FIND "${text}" "${line}" at)
   if(at EQUAL -1)
      message(FATAL_ERROR "meds_config_main.toml no longer holds the line this test edits: ${line}")
   endif()
endforeach()
string(REPLACE "fast_biophysics_on  = false" "fast_biophysics_on  = true" text "${text}")
string(REPLACE "# time_integrator     = \"esdirk2\"" "time_integrator = \"rk45\"" text "${text}")
string(REPLACE "# rtol                = 2.0e-3" "rtol = 1.0e-3" text "${text}")
file(WRITE "${WORK_FILE}" "${text}")
execute_process(COMMAND "${MEDS_MAIN}" "${WORK_FILE}" WORKING_DIRECTORY "${SOURCE_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
if(rc EQUAL 0 OR NOT log MATCHES "fast\\.rtol > 1e-4 with time_integrator = \"rk45\"")
   message(FATAL_ERROR "meds_main did not refuse rk45 at rtol 1e-3 as expected\n--- output:\n${log}")
endif()
message(STATUS "meds_main refused rk45 at rtol 1e-3")
