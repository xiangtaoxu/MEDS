#==========================================================================================#
# MEDS -- Modular Ecosystem Dynamics Simulator / Modern Ecosystem Demography Simulator        #
#         (standalone demographic core).                                                     #
#                                                                                          #
# CMake handles Fortran module-dependency ordering automatically (no hand-maintained object #
# lists, no repeated-build hack). The hot kernels carry explicit OpenMP `target` regions; the #
# build picks the device via the NVHPC `-mp` flag (MEDS_GPU = none|multicore|gpu). Explicit   #
# `map` + `-gpu=mem:separate` keep all state in normal host memory (no CUDA managed allocator).#
#==========================================================================================#
cmake_minimum_required(VERSION 3.20)
project(MEDS VERSION 0.1.0 LANGUAGES Fortran
        DESCRIPTION "Modular Ecosystem Dynamics Simulator (Modern Ecosystem Demography Simulator)")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(MEDS_BUILD_TESTS "Build the CTest test suite" ON)
option(MEDS_ENABLE_IO "Build the netCDF state-output layer (requires the netCDF C library)" OFF)
set(MEDS_GPU "none" CACHE STRING "OpenMP offload target for the kernels (NVHPC): none|multicore|gpu")
set_property(CACHE MEDS_GPU PROPERTY STRINGS none multicore gpu)

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

#----- Cross-cutting utilities (precision, constants, PFT params, config). ----------------#
file(GLOB SHARED_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/shared_util/*.f90)
add_library(meds_shared STATIC ${SHARED_SOURCES})
meds_fortran_flags(meds_shared)

#----- Self-contained demographic engine (building-block operations; links shared). -------#
file(GLOB DEMOG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/demography/*.f90)
add_library(meds_demography STATIC ${DEMOG_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_shared)
meds_fortran_flags(meds_demography)

#----- OpenMP-target offload of the kernels (NVHPC; the kernels live in meds_demography). --#
#       PUBLIC so the demo/tests/aux that link meds_demography inherit the offload flags.    #
#       `-gpu=mem:separate` forces explicit data movement (no CUDA managed allocator), which #
#       is the whole point: host state stays in normal memory.                               #
if(NOT MEDS_GPU STREQUAL "none")
   if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      if(MEDS_GPU STREQUAL "gpu")
         set(MEDS_MP_FLAGS -mp=gpu -gpu=mem:separate)
      else()
         set(MEDS_MP_FLAGS -mp)                         # multicore (host OpenMP)
      endif()
      target_compile_options(meds_demography PUBLIC ${MEDS_MP_FLAGS})
      target_link_options(meds_demography PUBLIC ${MEDS_MP_FLAGS})
      message(STATUS "MEDS: OpenMP target -> ${MEDS_MP_FLAGS} (NVHPC)")
   else()
      message(WARNING "MEDS_GPU=${MEDS_GPU} ignored for ${CMAKE_Fortran_COMPILER_ID} (NVHPC offload only)")
   endif()
endif()

#----- Always-on configuration INPUT: a minimal TOML reader + the meds_config.toml overlay.  #
#      No external dependency (unlike the netCDF OUTPUT layer below), so it is always built.  #
add_library(meds_config_io STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_toml.f90
                                  ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_config_io.f90)
target_link_libraries(meds_config_io PUBLIC meds_shared)
meds_fortran_flags(meds_config_io)

#----- Top-level driver + support utilities that wire the process modules together: the     #
#      master stepper (src/driver) plus the example rate provider, setup and diagnostics     #
#      (src/ root). Both globs are non-recursive.                                            #
file(GLOB AUX_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/*.f90)
add_library(meds_aux STATIC ${AUX_SOURCES})
target_link_libraries(meds_aux PUBLIC meds_demography meds_config_io)
meds_fortran_flags(meds_aux)

#----- Demo. ------------------------------------------------------------------------------#
add_executable(meds_demo app/meds_demo.f90)
target_link_libraries(meds_demo PRIVATE meds_aux)
meds_fortran_flags(meds_demo)

#----- Optional netCDF state output (src/io). Gated so the core build stays netCDF-free.     #
#      Uses the netCDF C library via ISO_C_BINDING -> pure C linkage, identical under ifx and #
#      nvfortran. Point CMake at the install with -DCMAKE_PREFIX_PATH=<netcdf-c prefix>.      #
if(MEDS_ENABLE_IO)
   find_package(netCDF CONFIG REQUIRED)
   add_library(meds_io STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_c.f90
                              ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_io.f90)
   target_link_libraries(meds_io PUBLIC meds_demography netCDF::netcdf)
   meds_fortran_flags(meds_io)

   add_executable(meds_io_demo app/meds_io_demo.f90)
   target_link_libraries(meds_io_demo PRIVATE meds_aux meds_io)
   meds_fortran_flags(meds_io_demo)
   #----- The C-binding layer passes array temporaries (NUL-terminated strings, contiguous   #
   #      slabs) to the netCDF C `(*)` dummies. That is correct for C interop, but Intel's    #
   #      Debug `-check all` emits an (harmless) arg_temp_created warning per call -- so build #
   #      the IO layer in Release (-DCMAKE_BUILD_TYPE=Release) for clean, production output.    #
   message(STATUS "MEDS: netCDF output enabled")
endif()

#----- Tests. -----------------------------------------------------------------------------#
if(MEDS_BUILD_TESTS)
   enable_testing()
   add_library(meds_testsupport STATIC test/meds_test_support.f90)
   target_link_libraries(meds_testsupport PUBLIC meds_demography)
   meds_fortran_flags(meds_testsupport)

   foreach(t containers fusion_cohort rates patch spinup allometry disturbance)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport meds_aux)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()
endif()
