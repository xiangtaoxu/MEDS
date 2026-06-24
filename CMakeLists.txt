#==========================================================================================#
# MEDS -- Modular Ecosystem Dynamics Simulator / Modern Ecosystem Demography Simulator        #
#         (standalone demographic core).                                                     #
#                                                                                          #
# CMake handles Fortran module-dependency ordering automatically (no hand-maintained object #
# lists, no repeated-build hack). Parallelism is expressed with `do concurrent`; the build  #
# selects the target via the NVHPC `-stdpar` flag (MEDS_STDPAR = none|multicore|gpu).        #
#==========================================================================================#
cmake_minimum_required(VERSION 3.20)
project(MEDS VERSION 0.1.0 LANGUAGES Fortran
        DESCRIPTION "Modular Ecosystem Dynamics Simulator (Modern Ecosystem Demography Simulator)")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(MEDS_BUILD_TESTS "Build the CTest test suite" ON)
set(MEDS_STDPAR "none" CACHE STRING "do concurrent target for NVHPC: none|multicore|gpu")
set_property(CACHE MEDS_STDPAR PROPERTY STRINGS none multicore gpu)

# All Fortran .mod files in one directory, visible to every target.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

#----- Per-compiler, per-config flags. ----------------------------------------------------#
function(meds_fortran_flags tgt)
   set(id ${CMAKE_Fortran_COMPILER_ID})
   if(id STREQUAL "IntelLLVM" OR id STREQUAL "Intel")
      target_compile_options(${tgt} PRIVATE -stand f18 -warn all -diag-disable=7712
         $<$<CONFIG:Debug>:-O0;-g;-check;all;-fpe0;-traceback>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "NVHPC" OR id STREQUAL "PGI")
      target_compile_options(${tgt} PRIVATE -Minform=warn
         $<$<CONFIG:Debug>:-O0;-g;-Mbounds;-Ktrap=fp>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "GNU")
      target_compile_options(${tgt} PRIVATE -std=f2018 -Wall -Wno-unused-dummy-argument
         $<$<CONFIG:Debug>:-O0;-g;-fcheck=all;-ffpe-trap=invalid,zero,overflow;-fbacktrace>
         $<$<CONFIG:Release>:-O2>)
   endif()
endfunction()

#----- Library. ---------------------------------------------------------------------------#
file(GLOB MEDS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90)
add_library(meds STATIC ${MEDS_SOURCES})
meds_fortran_flags(meds)

#----- do concurrent offload target (NVHPC only; PUBLIC so consumers inherit). ------------#
if(NOT MEDS_STDPAR STREQUAL "none")
   if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      target_compile_options(meds PUBLIC -stdpar=${MEDS_STDPAR})
      target_link_options(meds PUBLIC -stdpar=${MEDS_STDPAR})
      message(STATUS "MEDS: do concurrent -> -stdpar=${MEDS_STDPAR} (NVHPC)")
   else()
      message(WARNING "MEDS_STDPAR=${MEDS_STDPAR} ignored for ${CMAKE_Fortran_COMPILER_ID} (NVHPC only)")
   endif()
endif()

#----- Demo. ------------------------------------------------------------------------------#
add_executable(meds_demo app/meds_demo.f90)
target_link_libraries(meds_demo PRIVATE meds)
meds_fortran_flags(meds_demo)

#----- Tests. -----------------------------------------------------------------------------#
if(MEDS_BUILD_TESTS)
   enable_testing()
   add_library(meds_testsupport STATIC test/meds_test_support.f90)
   target_link_libraries(meds_testsupport PUBLIC meds)
   meds_fortran_flags(meds_testsupport)

   foreach(t containers fusion_cohort rates patch spinup)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()
endif()
