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
# netCDF output is ON by default (the standard build needs the netCDF C library). For a
# netCDF-free build -- quick compiler checks, CI, or debugging the engine on a machine without
# netCDF -- configure with -DMEDS_ENABLE_IO=OFF, which links a no-op stub I/O layer instead.
option(MEDS_ENABLE_IO "netCDF state output (needs the netCDF C library); -DMEDS_ENABLE_IO=OFF builds a netCDF-free stub for testing/debugging" ON)
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

#----- Foundation: precision, constants, PFT trait table, config, time (no heavy deps). ---#
file(GLOB SHARED_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/shared/*.f90)
add_library(meds_shared STATIC ${SHARED_SOURCES})
meds_fortran_flags(meds_shared)

#----- Plant allometry: pure size relationships (its own library under plant/; links shared).#
#      hgt_max is now a per-PFT argument, so config no longer depends on allometry -- allometry #
#      sits just above shared and is depended on by state / demography / vitality.              #
add_library(meds_allometry STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/meds_allometry.f90)
target_link_libraries(meds_allometry PUBLIC meds_shared)
meds_fortran_flags(meds_allometry)

#----- The cohort/patch STATE hub (SoA + CSR + lockstep machinery). Lifted OUT of the        #
#      demography process into its own layer: every science library depends on it; it pulls    #
#      no heavy deps.                                                                           #
add_library(meds_state STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/state/meds_demography_types.f90)
target_link_libraries(meds_state PUBLIC meds_shared meds_allometry)
meds_fortran_flags(meds_state)

#----- Demographic PROCESS engine (dynamics/structure/interface/diagnostics; acts on state). #
file(GLOB DEMOG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/demography/*.f90)
add_library(meds_demography STATIC ${DEMOG_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_state meds_allometry)
meds_fortran_flags(meds_demography)

#----- Self-contained leaf-level photosynthesis + stomatal conductance (links shared only).--#
#       Sealed like the demography engine and orthogonal to it: NO dependency on the cohort/  #
#       patch state (no site_t), so it compiles standalone (cmake --build ... --target          #
#       meds_leaf_physiology). The driver/test/demo link it directly.                           #
file(GLOB LEAF_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/leaf/*.f90)
#----- The C-API shim is compiled ONLY into the optional Python shared library (below), not   #
#      into the core static engine.                                                            #
list(REMOVE_ITEM LEAF_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/leaf/meds_leaf_capi.f90)
add_library(meds_leaf_physiology STATIC ${LEAF_SOURCES})
target_link_libraries(meds_leaf_physiology PUBLIC meds_shared)
meds_fortran_flags(meds_leaf_physiology)

#----- Self-contained canopy radiative transfer (links shared only, like the leaf module).  --#
#       A faithful reimplementation of ED2's two-stream (icanrad=2), modernized: one unified    #
#       multi-band solver (VIS/NIR/LW), SCOPE/4SAIL leaf-angle scattering (Beta leaf-angle       #
#       distribution), and an O(N) adding-method solve. State-free (no site_t), so it compiles   #
#       standalone (cmake --build ... --target meds_biophys).                                    #
file(GLOB BIOPHYS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/biophys/*.f90)
add_library(meds_biophys STATIC ${BIOPHYS_SOURCES})
target_link_libraries(meds_biophys PUBLIC meds_shared)
meds_fortran_flags(meds_biophys)

#----- Optional: a SHARED C-API library for calling the leaf model from Python (ctypes).     #
#      OFF by default (keeps normal builds foreign-call-free). Build with -DMEDS_BUILD_PYLIB=ON#
#      -> libmeds_leaf_c; the static deps are rebuilt position-independent so they link in.     #
option(MEDS_BUILD_PYLIB "Build the leaf-physiology C-API shared library for Python (ctypes)" OFF)
if(MEDS_BUILD_PYLIB)
   set_target_properties(meds_shared meds_leaf_physiology PROPERTIES POSITION_INDEPENDENT_CODE ON)
   add_library(meds_leaf_c SHARED ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/leaf/meds_leaf_capi.f90)
   target_link_libraries(meds_leaf_c PRIVATE meds_leaf_physiology)
   meds_fortran_flags(meds_leaf_c)
   message(STATUS "MEDS: leaf-physiology Python C-API ENABLED -> libmeds_leaf_c")
endif()

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
target_link_libraries(meds_config_io PUBLIC meds_shared meds_allometry)
meds_fortran_flags(meds_config_io)

#----- Plant VITALITY: the phenomenological vital-rate provider -- growth / mortality /       #
#      recruitment (reproduction folded into recruitment for now), split into per-process       #
#      formula modules behind one orchestrator (the vital_rates seam). Links state + allometry. #
file(GLOB VITALITY_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/vitality/*.f90)
add_library(meds_vitality STATIC ${VITALITY_SOURCES})
target_link_libraries(meds_vitality PUBLIC meds_state meds_allometry)
meds_fortran_flags(meds_vitality)

#----- Top-level driver + support utilities that wire the process modules together: the     #
#      master stepper (src/driver) and the initial-community builders incl. the census reader #
#      (src/init). meds_main.f90 is the PROGRAM (built as an executable below), so it is        #
#      excluded from the library.                                                            #
file(GLOB AUX_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/*.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/src/init/*.f90)
list(REMOVE_ITEM AUX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/meds_main.f90)
add_library(meds_aux STATIC ${AUX_SOURCES})
target_link_libraries(meds_aux PUBLIC meds_demography meds_vitality meds_config_io)
meds_fortran_flags(meds_aux)

#----- Output layer (src/io). meds_main always links a module named `meds_io`: the real        #
#      netCDF writer by default (the netCDF C library via ISO_C_BINDING -> pure C linkage,      #
#      identical under ifx and nvfortran; point CMake at it with -DCMAKE_PREFIX_PATH=<prefix>), #
#      or a no-op STUB under -DMEDS_ENABLE_IO=OFF -- so the single entry point builds and runs   #
#      with or without netCDF, and the core/tests stay netCDF-free in the OFF (debug) build.     #
if(MEDS_ENABLE_IO)
   #----- netCDF-C's CMake config pulls in HDF5 -> FindThreads, which needs a C/CXX language.   #
   #      Enable C here (only for the netCDF build; the stub/debug build stays Fortran-only).    #
   enable_language(C)
   find_package(netCDF CONFIG REQUIRED)
   add_library(meds_io STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_c.f90
                              ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_io.f90)
   target_link_libraries(meds_io PUBLIC meds_demography netCDF::netcdf)
   meds_fortran_flags(meds_io)
   set(MEDS_IO_LIB meds_io)
   #----- The C-binding layer passes array temporaries (NUL-terminated strings, contiguous   #
   #      slabs) to the netCDF C `(*)` dummies. That is correct for C interop, but Intel's    #
   #      Debug `-check all` emits an (harmless) arg_temp_created warning per call -- so build #
   #      the IO layer in Release (-DCMAKE_BUILD_TYPE=Release) for clean, production output.    #
   message(STATUS "MEDS: netCDF output ENABLED (default)")
else()
   add_library(meds_io_stub STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_io_stub.f90)
   target_link_libraries(meds_io_stub PUBLIC meds_demography)
   meds_fortran_flags(meds_io_stub)
   set(MEDS_IO_LIB meds_io_stub)
   message(STATUS "MEDS: netCDF output DISABLED -- stub I/O (testing/debugging build, no netCDF needed)")
endif()

#----- The main entry point: read config, run the simulation, save output, exit. -----------#
add_executable(meds_main ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/meds_main.f90)
target_link_libraries(meds_main PRIVATE meds_aux ${MEDS_IO_LIB})
meds_fortran_flags(meds_main)

#----- Tests. -----------------------------------------------------------------------------#
if(MEDS_BUILD_TESTS)
   enable_testing()
   add_library(meds_testsupport STATIC test/meds_test_support.f90)
   target_link_libraries(meds_testsupport PUBLIC meds_demography)
   meds_fortran_flags(meds_testsupport)

   foreach(t containers fusion_cohort rates patch spinup allometry disturbance init_census time)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport meds_aux)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Leaf-physiology test links the leaf library (not meds_aux). ----------------------#
   add_executable(test_leaf_physiology test/test_leaf_physiology.f90)
   target_link_libraries(test_leaf_physiology PRIVATE meds_testsupport meds_leaf_physiology)
   meds_fortran_flags(test_leaf_physiology)
   add_test(NAME leaf_physiology COMMAND test_leaf_physiology)

   #----- Canopy radiative-transfer test links the biophys library (state-free). -----------#
   add_executable(test_canopy_radiation test/test_canopy_radiation.f90)
   target_link_libraries(test_canopy_radiation PRIVATE meds_biophys)
   meds_fortran_flags(test_canopy_radiation)
   add_test(NAME canopy_radiation COMMAND test_canopy_radiation)
endif()
