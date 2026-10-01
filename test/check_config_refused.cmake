# SPDX-License-Identifier: Apache-2.0
# Runs meds_main on the reference config with three keys added that MEDS does not read, and fails
# unless the run stops before simulating and names each one (N-10): a misspelling with the key it
# meant, a key that moved section with its new home, and a retired key with what replaced it.
#
# Usage: cmake -DMEDS_MAIN=<meds_main> -DSOURCE_DIR=<repo> -DWORK_FILE=<scratch toml>
#              -P check_config_refused.cmake
cmake_minimum_required(VERSION 3.20)
file(READ "${SOURCE_DIR}/meds_config_main.toml" text)
string(APPEND text "\n[forcing]\ntimestpe = \"3600s\"\n[state]\ncohort_max = 2048\n[output]\nstrict_caps = true\n")
file(WRITE "${WORK_FILE}" "${text}")
execute_process(COMMAND "${MEDS_MAIN}" "${WORK_FILE}" WORKING_DIRECTORY "${SOURCE_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
set(problems "")
if(rc EQUAL 0)
   string(APPEND problems "\n  the run did not stop")
endif()
foreach(expect "forcing.timestpe .*did you mean forcing.timestep\\?"
               "state.cohort_max .*did you mean output.cohort_max\\?"
               "output.strict_caps .*is retired: .*raise the cap"
               "the configuration has missing, unknown or retired keys")
   if(NOT log MATCHES "${expect}")
      string(APPEND problems "\n  the report lacks: ${expect}")
   endif()
endforeach()
if(problems)
   message(FATAL_ERROR "meds_main did not refuse the unknown keys as expected:${problems}\n--- output:\n${log}")
endif()
message(STATUS "meds_main refused the three keys, naming each")
