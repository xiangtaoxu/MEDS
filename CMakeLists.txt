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
# netCDF output + I/O are ALWAYS compiled: the standard (and only) build needs the netCDF C library.
# Point CMake at it with -DCMAKE_PREFIX_PATH=<prefix> (find_package(netCDF CONFIG REQUIRED), below).
# (The former -DMEDS_ENABLE_IO=OFF netCDF-free stub build + meds_io_stub.f90 were removed by design.)
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

#----- Foundation: precision, constants, PFT trait table, config, time, temperature response. #
#      No heavy deps. meds_temp_response (Arrhenius / peaked deactivation) was promoted here     #
#      from the leaf module so leaf, respiration, and any other tissue share one code path        #
#      without a plant-domain -> plant-domain library edge.                                        #
file(GLOB SHARED_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/shared/*.f90)
add_library(meds_shared STATIC ${SHARED_SOURCES})
meds_fortran_flags(meds_shared)

#----- Structural allometry: pure size relationships. Its OWN library (src/allometry/), a      #
#      shared foundation BELOW state: state (cohort geometry caching / fusion) AND the plant     #
#      ecophysiology library both depend on it, so it cannot live inside libmeds_plant -- that   #
#      would make the demographic core depend on all of ecophysiology (see issue #11).           #
add_library(meds_allometry STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/allometry/meds_allometry.f90)
target_link_libraries(meds_allometry PUBLIC meds_shared)
meds_fortran_flags(meds_allometry)

#----- The cohort/patch STATE hub (SoA + CSR + lockstep machinery). Lifted OUT of the        #
#      demography process into its own layer: every science library depends on it; it pulls    #
#      no heavy deps.                                                                           #
add_library(meds_state STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/state/meds_demography_types.f90)
target_link_libraries(meds_state PUBLIC meds_shared meds_allometry)
meds_fortran_flags(meds_state)

#----- SELF-CONTAINED demography model: the PROCESS engine (dynamics / structure / interface /  #
#      diagnostics) PLUS the phenomenological empirical vital-rate provider (growth / mortality /  #
#      recruitment + the vital_rates orchestrator), moved here from the former plant/vitality/.    #
#      Acts on state; has NO plant-ecophysiology dependency, so it runs the spin-up on its own.    #
file(GLOB DEMOG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/demography/*.f90)
add_library(meds_demography STATIC ${DEMOG_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_state meds_allometry)
meds_fortran_flags(meds_demography)

#----- Diagnostic-aggregation CORE (src/io): the netCDF-FREE half of the output subsystem --      #
#      variable registry + per-tier temporal integrators + the extract switchboard. Links         #
#      meds_demography (site_t + the site reductions) ONLY -- NO netCDF -- so the stepper edge      #
#      that ticks the integrators pulls no C dependency (the §2 DAG-hygiene wall). The netCDF        #
#      serializer half lives in meds_io (below). Design: archive/MEDS_IO_DESIGN.md.                   #
add_library(meds_output_core STATIC
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_types.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_integrate.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_registry.f90)
target_link_libraries(meds_output_core PUBLIC meds_demography)
meds_fortran_flags(meds_output_core)

#----- ONE flat plant-ECOPHYSIOLOGY library: leaf gas exchange (photosynthesis + stomata +       #
#      solver), hydraulics, phenology (and respiration, when it lands). Links meds_shared ONLY --  #
#      no cohort/patch state (no site_t), so it compiles + unit-tests standalone                    #
#      (cmake --build ... --target meds_plant). The *_capi.f90 shims are compiled ONLY into the      #
#      optional Python shared library (below), never into the static core.                           #
file(GLOB PLANT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/*.f90)
file(GLOB PLANT_CAPI    CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/*_capi.f90)
list(REMOVE_ITEM PLANT_SOURCES ${PLANT_CAPI})
add_library(meds_plant STATIC ${PLANT_SOURCES})
target_link_libraries(meds_plant PUBLIC meds_shared)
meds_fortran_flags(meds_plant)

#----- Self-contained canopy radiative transfer (links shared only, like the plant module).  --#
#       A faithful reimplementation of ED2's two-stream (icanrad=2), modernized: one unified    #
#       multi-band solver (VIS/NIR/LW), SCOPE/4SAIL leaf-angle scattering (Beta leaf-angle       #
#       distribution), and an O(N) adding-method solve. State-free (no site_t), so it compiles   #
#       standalone (cmake --build ... --target meds_biophysics).                                 #
file(GLOB BIOPHYSICS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/biophysics/*.f90)
add_library(meds_biophysics STATIC ${BIOPHYSICS_SOURCES})
target_link_libraries(meds_biophysics PUBLIC meds_shared)
meds_fortran_flags(meds_biophysics)

#----- Self-contained biogeochemistry kernels (the column CARBON cycle: fast canopy-air CO2      #
#      exchange + slow soil-carbon pools). Links meds_shared ONLY, like biophysics, so it stays  #
#      stateless/pure/GPU-eligible and compiles standalone (cmake --build ... --target            #
#      meds_biogeochemistry). Design: archive/MEDS_COLUMN_CO2_BALANCE_DESIGN.md.                   #
file(GLOB BIOGEOCHEM_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/biogeochemistry/*.f90)
add_library(meds_biogeochemistry STATIC ${BIOGEOCHEM_SOURCES})
target_link_libraries(meds_biogeochemistry PUBLIC meds_shared)
meds_fortran_flags(meds_biogeochemistry)

#----- Optional: a SHARED C-API library for calling the plant ecophysiology from Python.     #
#      OFF by default (keeps normal builds foreign-call-free). Build with -DMEDS_BUILD_PYLIB=ON#
#      -> libmeds_plant_c; the static deps are rebuilt position-independent so they link in.    #
#      The pylib links meds_plant only (shared-only), so the cohort SoA is NOT pulled in.       #
option(MEDS_BUILD_PYLIB "Build the plant-ecophysiology C-API shared library for Python (ctypes)" OFF)
if(MEDS_BUILD_PYLIB)
   set_target_properties(meds_shared meds_plant PROPERTIES POSITION_INDEPENDENT_CODE ON)
   add_library(meds_plant_c SHARED ${PLANT_CAPI})
   target_link_libraries(meds_plant_c PRIVATE meds_plant)
   meds_fortran_flags(meds_plant_c)
   message(STATUS "MEDS: plant-ecophysiology Python C-API ENABLED -> libmeds_plant_c")
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

#----- netCDF C bindings, PROMOTED to their own minimal target so BOTH the output layer         #
#      (meds_io, state serialization) AND the forcing reader (meds_forcing) link the netCDF C     #
#      API WITHOUT meds_io's -> meds_demography edge. netCDF is a hard dependency (always built);  #
#      its CMake config pulls in HDF5 -> FindThreads, which needs the C language enabled.          #
enable_language(C)
find_package(netCDF CONFIG REQUIRED)
add_library(meds_netcdf_c STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_c.f90)
target_link_libraries(meds_netcdf_c PUBLIC meds_shared netCDF::netcdf)
meds_fortran_flags(meds_netcdf_c)

#----- Meteorological (and future disturbance / land-use) FORCING (src/forcing). The reader +    #
#      disaggregation kernels + runtime types. Links meds_shared + the netCDF bindings ONLY (NOT   #
#      demography): a prescribed external driver, low in the DAG. Design MEDS_FORCING_DESIGN.md.    #
file(GLOB FORCING_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/forcing/*.f90)
add_library(meds_forcing STATIC ${FORCING_SOURCES})
target_link_libraries(meds_forcing PUBLIC meds_shared meds_netcdf_c)
meds_fortran_flags(meds_forcing)

#----- Top-level driver + support utilities that wire the process modules together: the     #
#      master stepper (src/driver) and the initial-community builders incl. the census reader #
#      (src/init). meds_main.f90 is the PROGRAM (built as an executable below), so it is        #
#      excluded from the library.                                                            #
file(GLOB AUX_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/*.f90
                                        ${CMAKE_CURRENT_SOURCE_DIR}/src/init/*.f90)
list(REMOVE_ITEM AUX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/driver/meds_main.f90)
add_library(meds_aux STATIC ${AUX_SOURCES})
#----- Links meds_plant: the vegetation-dynamics driver's CARBON path calls the plant carbon    #
#      seam (get_plant_flux_slow) to assemble NPP; demography stays plant-free (the driver weaves).#
#      Also links meds_biophysics + meds_biogeochemistry: the fast-loop column-dynamics driver     #
#      (meds_column_dynamics) weaves the biophysics + biogeochemistry kernels each dt_fast.         #
target_link_libraries(meds_aux PUBLIC meds_demography meds_config_io meds_plant
                                      meds_biophysics meds_biogeochemistry meds_forcing
                                      meds_output_core)
meds_fortran_flags(meds_aux)

#----- Output layer (src/io): the netCDF state/diagnostic writer. Links the promoted netCDF        #
#      bindings (meds_netcdf_c, above) + demography (state serialization). netCDF is always         #
#      compiled. (The C-binding layer passes array temporaries to the netCDF C `(*)` dummies --      #
#      correct for interop, but Intel Debug `-check all` emits a harmless arg_temp_created warning   #
#      per call; build Release (-DCMAKE_BUILD_TYPE=Release) for clean output.)                       #
add_library(meds_io STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_io.f90
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_stream.f90
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_manager.f90)
target_link_libraries(meds_io PUBLIC meds_demography meds_netcdf_c meds_output_core)
meds_fortran_flags(meds_io)
set(MEDS_IO_LIB meds_io)
message(STATUS "MEDS: netCDF output (always compiled)")

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

   #----- Shared numerical primitives + conservation checker (links meds_shared only). -----#
   add_executable(test_numerics test/test_numerics.f90)
   target_link_libraries(test_numerics PRIVATE meds_shared)
   meds_fortran_flags(test_numerics)
   add_test(NAME numerics COMMAND test_numerics)

   foreach(t containers fusion_cohort rates patch spinup allometry disturbance init_census time carbon_growth fast_loop phenology_driver)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport meds_aux)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Leaf-physiology test links the flat plant library + test support (assertions). ---#
   add_executable(test_leaf_physiology test/test_leaf_physiology.f90)
   target_link_libraries(test_leaf_physiology PRIVATE meds_testsupport meds_plant)
   meds_fortran_flags(test_leaf_physiology)
   add_test(NAME leaf_physiology COMMAND test_leaf_physiology)

   #----- Hydraulics / phenology / respiration / carbon tests link the flat plant library (state-free). #
   foreach(t plant_hydraulics plant_phenology plant_respiration plant_carbon_dynamics)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_plant)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Canopy radiative-transfer test links the biophysics library (state-free). --------#
   add_executable(test_canopy_radiation test/test_canopy_radiation.f90)
   target_link_libraries(test_canopy_radiation PRIVATE meds_biophysics)
   meds_fortran_flags(test_canopy_radiation)
   add_test(NAME canopy_radiation COMMAND test_canopy_radiation)

   #----- Canopy-aerodynamics test links the biophysics library (state-free). --------------#
   add_executable(test_aerodynamics test/test_aerodynamics.f90)
   target_link_libraries(test_aerodynamics PRIVATE meds_biophysics)
   meds_fortran_flags(test_aerodynamics)
   add_test(NAME aerodynamics COMMAND test_aerodynamics)

   #----- Vertical soil-hydrology test links the biophysics library (state-free). ----------#
   add_executable(test_column_hydrology test/test_column_hydrology.f90)
   target_link_libraries(test_column_hydrology PRIVATE meds_biophysics)
   meds_fortran_flags(test_column_hydrology)
   add_test(NAME column_hydrology COMMAND test_column_hydrology)

   #----- Energy-balance tests link the biophysics library (state-free). -------------------#
   foreach(t column_energy surface_energy)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_biophysics)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Column CO2-balance test links the biogeochemistry library (state-free). ----------#
   add_executable(test_column_co2 test/test_column_co2.f90)
   target_link_libraries(test_column_co2 PRIVATE meds_biogeochemistry)
   meds_fortran_flags(test_column_co2)
   add_test(NAME column_co2 COMMAND test_column_co2)

   #----- Slow soil-carbon matrix test links the biogeochemistry library (state-free). -----#
   add_executable(test_soil_biogeochem test/test_soil_biogeochem.f90)
   target_link_libraries(test_soil_biogeochem PRIVATE meds_biogeochemistry)
   meds_fortran_flags(test_soil_biogeochem)
   add_test(NAME soil_biogeochem COMMAND test_soil_biogeochem)

   #----- Meteorological-forcing reader test links the forcing library (writes+reads a NetCDF). #
   add_executable(test_met_driver test/test_met_driver.f90)
   target_link_libraries(test_met_driver PRIVATE meds_forcing)
   meds_fortran_flags(test_met_driver)
   add_test(NAME met_driver COMMAND test_met_driver)

   #----- Fast-loop column-dynamics integration test (links the aux driver + test support). #
   add_executable(test_column_dynamics test/test_column_dynamics.f90)
   target_link_libraries(test_column_dynamics PRIVATE meds_aux meds_testsupport)
   meds_fortran_flags(test_column_dynamics)
   add_test(NAME column_dynamics COMMAND test_column_dynamics)

   #----- Pure fast-loop RHS (IMEX-ARK P0); links meds_aux for the meds_column_derivs module. ----#
   add_executable(test_column_derivs test/test_column_derivs.f90)
   target_link_libraries(test_column_derivs PRIVATE meds_aux)
   meds_fortran_flags(test_column_derivs)
   add_test(NAME column_derivs COMMAND test_column_derivs)

   add_executable(test_picard_coupling test/test_picard_coupling.f90)
   target_link_libraries(test_picard_coupling PRIVATE meds_aux meds_testsupport)
   meds_fortran_flags(test_picard_coupling)
   add_test(NAME picard_coupling COMMAND test_picard_coupling)

   #----- Diagnostic-output netCDF-free core: integrator kernels + registry (state-free of C). #
   foreach(t output_integrate output_registry)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_output_core meds_testsupport)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Diagnostic-output serializer round-trip (writes + reopens netCDF; links the io layer). #
   add_executable(test_output_roundtrip test/test_output_roundtrip.f90)
   target_link_libraries(test_output_roundtrip PRIVATE meds_io meds_testsupport)
   meds_fortran_flags(test_output_roundtrip)
   add_test(NAME output_roundtrip COMMAND test_output_roundtrip)
endif()
