# SPDX-License-Identifier: Apache-2.0
# Runs meds_main on the reference config set to start from a census (INIT_MODE 1) or a state file
# (INIT_MODE 2) that cannot start the run, and fails unless the run stops before simulating with a
# message naming the mode, the file and the REASON. A run that asked for a stand never starts from
# bare ground instead. With EMPTY_CENSUS the script first writes a census with a header and no row.
#
# Usage: cmake -DMEDS_MAIN=<meds_main> -DSOURCE_DIR=<repo> -DWORK_FILE=<scratch toml>
#              -DINIT_MODE=<1|2> -DINIT_FILE=<path> [-DEMPTY_CENSUS=ON] -DREASON=<text>
#              -P check_init_refused.cmake
cmake_minimum_required(VERSION 3.20)
if(INIT_MODE EQUAL 1)
   set(key "census_file")
   set(expect "census (mode 1): the census file ${INIT_FILE} ${REASON}")
else()
   set(key "restart_file")
   set(expect "restart (mode 2): the state file ${INIT_FILE} ${REASON}")
endif()
if(EMPTY_CENSUS)
   file(WRITE "${INIT_FILE}" "site_id,patch_id,cohort_id,dbh,height,pft,nplant\n")
else()
   file(REMOVE "${INIT_FILE}")
endif()
file(READ "${SOURCE_DIR}/meds_config_main.toml" text)
string(REGEX REPLACE "\ninit_mode *= *[0-9]" "\ninit_mode = ${INIT_MODE}\n${key} = \"${INIT_FILE}\"" text "${text}")
file(WRITE "${WORK_FILE}" "${text}")
execute_process(COMMAND "${MEDS_MAIN}" "${WORK_FILE}" WORKING_DIRECTORY "${SOURCE_DIR}"
                OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE rc)
set(log "${out}${err}")
set(problems "")
if(rc EQUAL 0)
   string(APPEND problems "\n  the run did not stop")
endif()
string(FIND "${log}" "bare ground" at)
if(NOT at EQUAL -1)
   string(APPEND problems "\n  the run started from bare ground")
endif()
string(FIND "${log}" "${expect}" at)
if(at EQUAL -1)
   string(APPEND problems "\n  the report lacks: ${expect}")
endif()
if(problems)
   message(FATAL_ERROR "meds_main did not refuse the initial-state file as expected:${problems}\n--- output:\n${log}")
endif()
message(STATUS "meds_main stopped, naming the file: ${expect}")
