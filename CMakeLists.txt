# SPDX-License-Identifier: Apache-2.0
#==========================================================================================#
# MEDS -- Modular Ecosystem Dynamics Simulator / Modern Ecosystem Demography Simulator        #
#         (standalone demographic core).                                                     #
#                                                                                          #
# CMake handles Fortran module-dependency ordering automatically (no hand-maintained object #
# lists, no repeated-build hack).                                                            #
#                                                                                          #
# THE PARALLEL PATH IS THE CPU, not the GPU. Patch-axis OpenMP threading on the host is what  #
# ships and what is byte-identical at any thread count. GPU offload was evaluated and is NOT   #
# viable as scoped: measured, the offload build runs 1.4x SLOWER than the CPU (49.8 s against  #
# 36.2 s), one kernel sits at 0.4% occupancy, and the device treated as 20 slow cores is 26x   #
# slower than 4 CPU cores. MEDS_GPU=gpu is kept as a reproducible experiment, not a fast path. #
# See docs/dev_plans/archive/MEDS_GPU_EVALUATION.md.                                           #
#==========================================================================================#
cmake_minimum_required(VERSION 3.20)
project(MEDS VERSION 0.3.5 LANGUAGES Fortran
        DESCRIPTION "Modular Ecosystem Dynamics Simulator (Modern Ecosystem Demography Simulator)")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(MEDS_BUILD_TESTS "Build the CTest test suite" ON)
# netCDF output + I/O are ALWAYS compiled: the standard (and only) build needs the netCDF C library.
# Point CMake at it with -DCMAKE_PREFIX_PATH=<prefix> (find_package(netCDF CONFIG REQUIRED), below).
# `multicore` is host OpenMP under NVHPC. `gpu` is the REFUTED offload experiment (see the header).
set(MEDS_GPU "none" CACHE STRING "NVHPC OpenMP target: none|multicore|gpu (gpu is measured SLOWER than CPU)")
set_property(CACHE MEDS_GPU PROPERTY STRINGS none multicore gpu)

# Optional Python C-API shared libraries (ctypes). OFF by default -> normal builds are foreign-
# call-free. When ON, ALL static libs are built position-independent so they link into the .so.
option(MEDS_BUILD_PYLIB "Build the C-API shared libraries for Python (ctypes)" OFF)
if(MEDS_BUILD_PYLIB)
   set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# All Fortran .mod files in one directory, visible to every target.
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

#----- Per-compiler, per-config flags. ----------------------------------------------------#
#
#  THE AUTOMATIC-ARRAY STACK CEILING (opt-in flag below).
#
#  The fast loop carries ~740 B of automatic (stack) arrays per COHORT, live simultaneously across
#  four nested frames (fast_step -> build_column_frozen -> column_prepass ->
#  aero_bottom_to_top). On an 8 MB stack that is a hard SIGSEGV ceiling near 10-12k cohorts, and in
#  a threaded run the ceiling follows each thread's stack, OMP_STACKSIZE, instead.
#
#  PREFERRED FIX IS A RUNTIME ONE -- raise the stack rather than change codegen:
#      ulimit -s unlimited          (or e.g. `ulimit -s 65536` for ~90k cohorts)
#      export OMP_STACKSIZE=512M    (a threaded run)
#
#  -DMEDS_HEAP_ARRAYS=ON moves large automatic arrays to the heap instead. It is OFF by default
#  because ifx's -heap-arrays SEGFAULTS test_output_registry at the exit of apply_variable_override
#  in a Debug build -- a compiler codegen interaction, not a MEDS bug (the routine only assigns
#  plain integer components, and the same code is clean without the flag). Turn it on only if the
#  runtime fix is unavailable, and re-run the Debug suite when you do.
#-----------------------------------------------------------------------------------------#
#
#  AUTOMATIC (STACK) LOCALS ARE MANDATORY WHEN THE PATCH LOOP IS THREADED.
#
#  It is the hardest trap here to find, so it is written down. Intel Fortran defaults to
#  `-auto-scalar`: local SCALARS go on the stack, but local ARRAYS and DERIVED TYPES go in STATIC
#  storage. Static locals are shared by every thread. The fast loop calls kernels in meds_fast_kernels,
#  meds_slow_kernels, meds_share and meds_demography, and every one with a local array or a local
#  derived type (e.g. `rad_forcing_t rf` in apply_rt_forcing -> alloc_rad_forcing) then has all
#  threads writing the same object. It does NOT crash: it produces plausible, silently
#  thread-count-dependent answers -- precisely what the byte-identity-at-any-thread-count requirement
#  exists to catch. ThreadSanitizer finds it; a green test suite does not.
#
#  `-qopenmp`/`-fopenmp` implies automatic locals, but only for the target it is applied to, so
#  flagging just the target that owns meds_fast_dynamics is NOT enough -- the whole call tree needs
#  it. MEDS_OPENMP therefore adds the per-compiler "all locals on the stack" flag to EVERY target.
#  Corollary: more goes on the stack, so the OMP_STACKSIZE note above applies with more force.
#-----------------------------------------------------------------------------------------#
option(MEDS_HEAP_ARRAYS "Put large automatic arrays on the heap (see the note above)" OFF)
option(MEDS_OPENMP "Compile the patch (one polygon) and polygon (several) parallelism; serial if the compiler has none" ON)
#      ON BY DEFAULT. The thread count is a RUNTIME setting ([run].n_threads, default 1), so an OpenMP
#      build still runs serially until a config asks for threads, and its output is byte-identical at
#      any thread count. Only an OpenMP build lets test_fast_loop's 4-thread check test anything: in a
#      serial build the directives are comments and it passes trivially. The price is the automatic-
#      locals flag on every target, measured in docs/building.md "Parallel builds". -DMEDS_OPENMP=OFF
#      builds serial; a compiler with no Fortran OpenMP falls back to serial with a warning.
set(MEDS_OPENMP_ACTIVE OFF)
if(MEDS_OPENMP)
   if(NOT MEDS_GPU STREQUAL "none" AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      set(MEDS_OPENMP_ACTIVE ON)                  # -mp is already PUBLIC on meds_demography (below)
   else()
      find_package(OpenMP COMPONENTS Fortran)
      if(OpenMP_Fortran_FOUND)
         set(MEDS_OPENMP_ACTIVE ON)
      else()
         message(WARNING "MEDS: ${CMAKE_Fortran_COMPILER_ID} offers no Fortran OpenMP, so this build is "
                         "SERIAL and [run].n_threads has no effect. Configure with -DMEDS_OPENMP=OFF "
                         "to ask for that and silence this warning.")
      endif()
   endif()
endif()
#      Any build that can run host threads needs automatic locals -- which includes NVHPC with
#      MEDS_GPU != none, because that puts -mp PUBLIC on meds_demography and meds_fast inherits it, so
#      patch-loop directives are live there whether or not MEDS_OPENMP was asked for.
set(MEDS_AUTO_LOCALS OFF)
if(MEDS_OPENMP_ACTIVE OR NOT MEDS_GPU STREQUAL "none")
   set(MEDS_AUTO_LOCALS ON)
endif()
function(meds_fortran_flags tgt)
   set(id ${CMAKE_Fortran_COMPILER_ID})
   if(id STREQUAL "IntelLLVM" OR id STREQUAL "Intel")
      target_compile_options(${tgt} PRIVATE -stand f18 -warn all -diag-disable=7712
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-heap-arrays;256>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-auto>
         $<$<CONFIG:Debug>:-O0;-g;-check;all;-fpe0;-traceback>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "NVHPC" OR id STREQUAL "PGI")
      target_compile_options(${tgt} PRIVATE -Minform=warn
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-Mnostack_arrays>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-Mrecursive>
         $<$<CONFIG:Debug>:-O0;-g;-Mbounds;-Ktrap=fp>
         $<$<CONFIG:Release>:-O2>)
   elseif(id STREQUAL "GNU")
      target_compile_options(${tgt} PRIVATE -std=f2018 -Wall -Wno-unused-dummy-argument
         $<$<BOOL:${MEDS_HEAP_ARRAYS}>:-fmax-stack-var-size=262144>
         $<$<BOOL:${MEDS_AUTO_LOCALS}>:-frecursive>
         $<$<CONFIG:Debug>:-O0;-g;-fcheck=all;-ffpe-trap=invalid,zero,overflow;-fbacktrace>
         $<$<CONFIG:Release>:-O2>)
   endif()
endfunction()

#----- FOUNDATION (src/share/): the kinds, the constants, the numerics, the calendar and the sun's #
#      position, and the budget checks at the top, and the stateless constitutive laws in         #
#      functions/, one library per domain: allometry, plant hydraulics, soil hydraulics,           #
#      thermodynamics (moist air and its lapse, with the temperature response of enzyme rates),    #
#      canopy optics and the radiation above the canopy.                                           #
#                                                                                                  #
#      One folder, one library, and it links NOTHING. A module belongs in src/share/ only if it    #
#      `use`s nothing outside src/share/ -- no model state, no configuration, no external          #
#      library. That is what lets every layer above, and the tests that link one narrow library,   #
#      reach it freely, and it is the test that keeps this folder from becoming the place a file   #
#      goes when no other place fits. The temperature response sits here so leaf, respiration and  #
#      any other tissue share one code path without a plant-domain -> plant-domain library edge;   #
#      meds_allometry_lib because it is a geometric CONSTRAINT that state, demography, plant and   #
#      config all read.                                                                            #
file(GLOB SHARE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/share/*.f90
                                          ${CMAKE_CURRENT_SOURCE_DIR}/src/share/functions/*.f90)
add_library(meds_share STATIC ${SHARE_SOURCES})
meds_fortran_flags(meds_share)

#----- CONFIGURATION: the meds_config_t aggregator, the loader that fills it from two TOML files, #
#      the checks it passes and the minimal TOML reader; in params/, one leaf per block (the PFT   #
#      table, the soil, forcing, output, the fast loop, snow and aerodynamics, leaf,               #
#      biogeochemistry, maps, region, sites). `toml` and `config_io` live here,                    #
#      not in src/io/:                                                                             #
#      they are configuration INPUT, not I/O of model results, and they carry no external          #
#      dependency, unlike the netCDF output layer. Links meds_share only; nothing here imports     #
#      model state or a kernel.                                                                    #
file(GLOB CONFIG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/config/*.f90
                                           ${CMAKE_CURRENT_SOURCE_DIR}/src/config/params/*.f90)
add_library(meds_config STATIC ${CONFIG_SOURCES})
#      The keys a config may hold are the ones the two reference configs list (meds_config_check): #
#      written into an include at configure time, and re-written when either reference changes.   #
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/meds_config_keys.cmake)
meds_write_config_keys(${CMAKE_CURRENT_SOURCE_DIR}/meds_config_main.toml
                       ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml
                       ${CMAKE_CURRENT_BINARY_DIR}/generated/meds_config_keys.inc)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
             ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_main.toml ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml)
target_include_directories(meds_config PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(meds_config PUBLIC meds_share)
meds_fortran_flags(meds_config)

#----- STATE, one LAYER in five files by model hierarchy (src/README.md, rule 3). -------------#
#                                                                                            #
#  meds_cohort_state_vars  -- the flat cohort SoA and the cohort lockstep.                    #
#  meds_patch_state_vars   -- ONE patch's stores (canopy air, snow, soil water and heat,      #
#                             soil carbon), their blends, and the patch SoA with its CSR map. #
#  meds_polygon_state_vars -- polygon_t: its site, patch and cohort blocks, and its numbers.  #
#  meds_site_state_vars    -- the site block: each site's area and its slice of the patches.  #
#  meds_diag_vars          -- the diagnostic accumulators, per cohort, patch and polygon.     #
#                                                                                            #
#  The kernels take only the per-patch store types (from meds_patch_state_vars), never        #
#  polygon_t, which keeps them OpenMP-`target` device-eligible; ctest kernels_no_site_state    #
#  checks it, because the shared module directory cannot. The parameters that describe the    #
#  stores are config (meds_soil_params): a parameter is not state.                            #
file(GLOB STATE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/state/*.f90)
add_library(meds_state STATIC ${STATE_SOURCES})
target_link_libraries(meds_state PUBLIC meds_share meds_config)
meds_fortran_flags(meds_state)



#==========================================================================================#
#  PROCESSES, TIMESCALE-FIRST (src/README.md, rule 2).                                        #
#                                                                                            #
#  MEDS is an operator-split fast/slow model -- that is its dominant organizing axis -- and    #
#  the tree says so. Each half holds its domain kernels, and the fast half additionally         #
#  holds the integrator machinery (integrator/) and the loop that drives it (driver/).            #
#==========================================================================================#

#----- SLOW kernels (daily / annual): phenology, carbon allocation, trait dynamics, the CENTURY #
#      soil-carbon matrix and the litter-destination split. Stateless and `pure`-where-possible; #
#      links the patch store types it reads/writes plus the config leaves. Builds standalone.   #
file(GLOB SLOW_KERNEL_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/plant/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/soil/*.f90)
add_library(meds_slow_kernels STATIC ${SLOW_KERNEL_SOURCES})
target_link_libraries(meds_slow_kernels PUBLIC meds_state meds_config)
meds_fortran_flags(meds_slow_kernels)

#----- FAST kernels (sub-daily): canopy radiative transfer + aerodynamics + the canopy air space #
#      (canopy/ = the medium), leaf gas exchange + hydraulics + maintenance respiration + the     #
#      tissue energy balance (plant/ = the organisms), and the soil water / soil energy / ground  #
#      + snow stores (soil/). State-free apart from the patch store types; builds standalone.     #
#                                                                                                  #
#      It links meds_state and meds_config only: the fast kernels import nothing from                #
#      slow_dynamics/, and the one fast->slow edge (the heterotrophic-respiration seam,              #
#      docs/science/soil_carbon.md §7) sits on meds_fast.                                            #
file(GLOB FAST_KERNEL_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/canopy/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/plant/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/soil/*.f90)
add_library(meds_fast_kernels STATIC ${FAST_KERNEL_SOURCES})
target_link_libraries(meds_fast_kernels PUBLIC meds_state meds_config)
meds_fortran_flags(meds_fast_kernels)

#----- DEMOGRAPHY: the vital-rate laws (growth / mortality / recruitment) and the apply-operators #
#      that change the population -- state update, cohort and patch fuse-fission, termination,     #
#      recruitment, treefall disturbance -- plus the public facade. Links the site state it         #
#      mutates and the litter law its kills feed. src/README.md rule 4 keeps the two halves apart:  #
#      the operators take rate ARRAYS and never import meds_demography_rates.                        #
file(GLOB DEMOGRAPHY_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/demography/*.f90)
add_library(meds_demography STATIC ${DEMOGRAPHY_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_state meds_slow_kernels)
meds_fortran_flags(meds_demography)

#----- Diagnostic output (src/io/output), its netCDF-free half: every module of the folder but the  #
#      writer -- the types, the registry, the file set's setup, the extract, the reductions and the  #
#      per-tier integrators. Links meds_demography only, so the stepper edge that ticks the          #
#      integrators pulls no C dependency (src/README.md, "The library graph"). The writer is         #
#      meds_output_stream (below). The stages: docs/science/diagnostics.md §1.                      #
file(GLOB OUTPUT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/*.f90)
list(REMOVE_ITEM OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/meds_output_stream.f90)
add_library(meds_output STATIC ${OUTPUT_SOURCES})
target_link_libraries(meds_output PUBLIC meds_demography)
meds_fortran_flags(meds_output)





#----- NVHPC OpenMP target selection (the flags live on meds_demography). -------------------#
#       PUBLIC so the demo/tests/drivers that link meds_demography inherit them.                       #
#       `-gpu=mem:separate` forces explicit data movement (no CUDA managed allocator), so host #
#       state stays in normal memory.                                                          #
#       MEDS_GPU=gpu is a REPRODUCIBLE EXPERIMENT, not the fast path: it measured 1.4x slower  #
#       than the CPU build. Use MEDS_GPU=multicore, or plain MEDS_OPENMP, for real runs.       #
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


#----- netCDF C bindings, in their own minimal target so BOTH the output layer                  #
#      (meds_init_restart, state serialization) AND the forcing reader (meds_met) link the netCDF C     #
#      API WITHOUT meds_output_stream's -> meds_demography edge. netCDF is a hard dependency;            #
#      its CMake config pulls in HDF5 -> FindThreads, which needs the C language enabled.          #
enable_language(C)
find_package(netCDF CONFIG REQUIRED)
add_library(meds_netcdf_c STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_c.f90
                                 ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_netcdf_write.f90)
target_link_libraries(meds_netcdf_c PUBLIC meds_share netCDF::netcdf)
meds_fortran_flags(meds_netcdf_c)

#----- Meteorological FORCING (src/io/met): the reader, its samples at a site and at a patch's   #
#      canopy-air top, the grid and the runtime types; the laws are share/'s. Links meds_config   #
#      + the netCDF bindings ONLY (NOT demography): a prescribed external driver, low in the DAG. #
#      The science: docs/science/forcing.md.                                                     #
file(GLOB FORCING_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/met/*.f90)
add_library(meds_met STATIC ${FORCING_SOURCES})
target_link_libraries(meds_met PUBLIC meds_config meds_netcdf_c)
meds_fortran_flags(meds_met)
#----- Static per-cell MAPS (src/io/maps): the one module that opens a map -- the cell       #
#      mask, and each cell's soil and height-allometry inputs.                                     #
#      Links meds_met for the archive grid it checks maps against and the cell list it reads   #
#      them at; the driver links it.                                                                 #
file(GLOB MAPS_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/maps/*.f90)
add_library(meds_maps STATIC ${MAPS_SOURCES})
target_link_libraries(meds_maps PUBLIC meds_met meds_config meds_netcdf_c)
meds_fortran_flags(meds_maps)


#----- The FAST loop: the integrator machinery (integrator/ -- the state vector, the frozen work  #
#      records, the ARK and RK45 marches, the step controller and the pre-pass) and the         #
#      driver that walks a slow step in dt_fast sub-steps over the patch axis (driver/).          #
#                                                                                                  #
#      meds_fast_prepass sits in integrator/, not driver/, because the step calls it: both schemes   #
#      build their frozen record through build_column_frozen, which runs column_prepass once per   #
#      dt_fast, so filing it under driver/ would make a cycle. It links meds_slow_kernels for the  #
#      one fast->slow seam                                                                          #
#      (docs/science/soil_carbon.md §7): the sub-daily heterotrophic respiration must respire       #
#      the SAME CENTURY matrix the daily step debits.                                               #
file(GLOB FAST_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/integrator/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/driver/*.f90)
add_library(meds_fast STATIC ${FAST_SOURCES})
target_link_libraries(meds_fast PUBLIC meds_fast_kernels meds_slow_kernels meds_state
                                       meds_met meds_output)
meds_fortran_flags(meds_fast)

#----- The SLOW loop: the thin coordinator and the two peer domain drivers (vegetation dynamics, #
#      soil biogeochemistry). It is the ONE place a demographic rate meets its application.       #
file(GLOB SLOW_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/driver/*.f90)
add_library(meds_slow STATIC ${SLOW_SOURCES})
target_link_libraries(meds_slow PUBLIC meds_demography meds_slow_kernels meds_fast_kernels)
meds_fortran_flags(meds_slow)

#----- The initial state (src/io/init): bare ground, a census or a restart file, and           #
#      init_polygon, which starts a polygon from the one [init].init_mode names. It links       #
#      meds_slow for a census stand's restructuring and trait acclimation, and the netCDF        #
#      bindings for the restart file, which the checkpoints write.                               #
file(GLOB INIT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/io/init/*.f90)
add_library(meds_init STATIC ${INIT_SOURCES})
target_link_libraries(meds_init PUBLIC meds_demography meds_slow meds_config meds_netcdf_c)
meds_fortran_flags(meds_init)

#----- src/main: the master stepper (the cadence owner -- one call advances the coupled model by  #
#      one slow step, sub-stepping the fast loop inside it) and the program that drives it. The    #
#      stepper is a library because the integration tests drive it directly; meds_main.f90 is the   #
#      PROGRAM and is the only thing in the executable.                                             #
add_library(meds_stepper STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_stepper.f90)
target_link_libraries(meds_stepper PUBLIC meds_fast meds_slow)
meds_fortran_flags(meds_stepper)

#----- The whole coupled model, as one name for everything above it to link. -----------------#
add_library(meds_model INTERFACE)
target_link_libraries(meds_model INTERFACE meds_stepper meds_fast meds_slow meds_init)

#----- HOST OpenMP over the fast-loop PATCH axis (src/fast_dynamics/driver/meds_fast_dynamics).  #
#      ON by default (the option, beside meds_fortran_flags, says why). Applied to meds_fast -- the #
#      target that owns the fast loop -- and PUBLIC so meds_main and the tests inherit the link     #
#      flags. find_package, run with the option, resolves the per-compiler flag (-qopenmp         #
#      IntelLLVM, -fopenmp GNU, -mp NVHPC).                                                         #
#      NVHPC with MEDS_GPU != none ALREADY carries -mp PUBLIC on meds_demography, inherited here,   #
#      so adding it again would fight -mp=gpu; that case is left alone. Note the consequence: such a #
#      build compiles the patch-loop directives whether or not MEDS_OPENMP is set, which is exactly  #
#      why the thread count is a RUNTIME knob ([run].n_threads, default 1) and not a build-time one. #
#      The option itself is declared beside meds_fortran_flags(), which is where its OTHER (and more #
#      important) half lives: the automatic-locals flag on every target.                             #
if(MEDS_OPENMP_ACTIVE)
   if(NOT MEDS_GPU STREQUAL "none" AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      message(STATUS "MEDS: fast-loop OpenMP inherited from MEDS_GPU=${MEDS_GPU} (-mp is already PUBLIC on meds_demography)")
   else()
      target_link_libraries(meds_fast PUBLIC OpenMP::OpenMP_Fortran)
      message(STATUS "MEDS: fast-loop OpenMP -> ${OpenMP_Fortran_FLAGS} (patch axis; [run].n_threads)")
   endif()
else()
   message(STATUS "MEDS: serial build; [run].n_threads has no effect")
endif()

#----- Diagnostic output's writer (src/io/output/meds_output_stream.f90): the netCDF files. Links  #
#      the netCDF bindings (meds_netcdf_c, above). netCDF is always compiled.                        #
#      (The C-binding layer passes array temporaries to the netCDF C `(*)` dummies --                #
#      correct for interop, but Intel Debug `-check all` emits a harmless arg_temp_created warning   #
#      per call; build Release (-DCMAKE_BUILD_TYPE=Release) for clean output.)                       #
add_library(meds_output_stream STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/output/meds_output_stream.f90)
target_link_libraries(meds_output_stream PUBLIC meds_demography meds_netcdf_c meds_output)
meds_fortran_flags(meds_output_stream)
message(STATUS "MEDS: netCDF output (always compiled)")

#----- The coupled model as OPEN / STEP / FINALIZE objects: one polygon's state and step           #
#      (meds_polygon) and the run of one polygon or many (meds_region). A LIBRARY rather than part  #
#      of the program because it has two callers: the executable below, and the C-API shim         #
#      src/c_api/meds_c_api_run.f90 that drives the full coupled model from Python. It sits above   #
#      meds_model + the IO layer -- the one place that owns configuration, the calendar, the output #
#      files and the closing conservation reports together.                                         #
add_library(meds_driver STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_polygon.f90
                        ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_region.f90)
target_link_libraries(meds_driver PUBLIC meds_model meds_maps meds_output_stream)
meds_fortran_flags(meds_driver)

#----- The main entry point: parse the command line, drive the model, exit. ----------------#
add_executable(meds_main ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_main.f90)
target_link_libraries(meds_main PRIVATE meds_driver)

#----- THE Python C-API shared library (opt-in): ONE libmeds.so for the whole package.        #
#                                                                                              #
#      Per-subsystem independence lives at the CALL level (verbs, not libraries), not the LINK #
#      level: `meds.plant`, `meds.canopy`, `meds.demography` and `meds.model` are sub-packages #
#      over one backend. A kernels-only library beside a site-driver one would force every     #
#      domain kernel library to link meds_share ONLY, and the column state types into          #
#      src/share/. The accepted cost: the library links meds_model + IO, so even the leaf      #
#      shim carries the whole model, and the wheel declares netCDF4 a runtime dependency       #
#      rather than vendoring it.                                                               #
#                                                                                              #
#      Plain GLOB, no exclusions: every src/c_api/*.f90 belongs in the one library.            #
if(MEDS_BUILD_PYLIB)
   file(GLOB C_API_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/*.f90)
   add_library(meds_py SHARED ${C_API_SOURCES})
   set_target_properties(meds_py PROPERTIES OUTPUT_NAME meds)
   target_link_libraries(meds_py PRIVATE meds_driver)
   meds_fortran_flags(meds_py)
   #----- The same numbers as meds_main, in any process. An Intel build calls `exp`, `sin`, `pow`  #
   #      and the rest by their C names, and meds_main finds them in Intel's math library. A       #
   #      library loaded into Python would find glibc's first, because Python has already loaded  #
   #      it, so without this a run through `meds.model` differs from the same run of meds_main    #
   #      in the last digits (BCI: up to 7e-4 W m-2 in LE over ten days). Linking Intel's          #
   #      runtime into the library and keeping its names inside it makes every call land where     #
   #      meds_main's does.                                                                        #
   #      gfortran builds use glibc's math in both, so they need nothing here.                    #
   if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
      target_link_options(meds_py PRIVATE -static-intel LINKER:--exclude-libs,ALL)
   endif()
   #----- RPATH so the INSTALLED library resolves its own dependencies.                         #
   #                                                                                            #
   #      libmeds.so links netCDF, HDF5 and the Fortran runtime. Declaring netCDF4 rather than     #
   #      vendoring the ~67 MB closure is necessary, but on its own NOT enough: the pip netCDF4    #
   #      wheel bundles its libnetcdf under a mangled soname inside its own package, so it cannot  #
   #      satisfy a `libnetcdf.so.19` NEEDED by a different .so. Freshly installed into a bare     #
   #      venv, the library without this RPATH fails to load on ten unresolved sonames.            #
   #                                                                                            #
   #      USE_LINK_PATH writes the directories of the actually-linked externals into the install  #
   #      RPATH, so the wheel finds them with no LD_LIBRARY_PATH and no `source setvars.sh`.      #
   #      $ORIGIN comes first so a co-bundled copy would still win.                                #
   #                                                                                            #
   #      THIS MAKES THE WHEEL MACHINE-LOCAL: those paths are the build machine's conda prefix    #
   #      and oneAPI install. That is the right trade for `pip install python/` from a source      #
   #      checkout, which is how MEDS is used; it is NOT a redistributable manylinux wheel, and    #
   #      building one would mean auditwheel plus the vendoring declined above.                    #
   #      USE_LINK_PATH alone is NOT enough: it captures the directories of libraries CMake was  #
   #      told to link (netCDF, HDF5), but the compiler's own shared runtime (the OpenMP runtime  #
   #      in an Intel build, which stays shared; libgfortran in a GNU one) is added implicitly by #
   #      the compiler driver, so it is absent from the RPATH and the library still fails to       #
   #      dlopen outside a `source setvars.sh` shell. CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES is   #
   #      exactly that set. (`ldd` will NOT show this: run from a sourced shell it inherits        #
   #      LD_LIBRARY_PATH and reports everything resolved. Check with `env -u LD_LIBRARY_PATH      #
   #      ldd`, or by importing.)                                                                  #
   set_target_properties(meds_py PROPERTIES
                         INSTALL_RPATH "$ORIGIN;${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
                         INSTALL_RPATH_USE_LINK_PATH TRUE)
   #----- INSTALL it INTO the package directory, so a wheel carries meds/libmeds.so and         #
   #      meds._libmeds finds it as a sibling of __init__.py -- the first entry in its search    #
   #      order, needing no env var. The `pylib` component is what python/pyproject.toml         #
   #      installs; nothing else in this project is part of the wheel.                            #
   install(TARGETS meds_py LIBRARY DESTINATION meds COMPONENT pylib)
   message(STATUS "MEDS: Python C-API ENABLED -> libmeds.so (meds.plant, meds.canopy, meds.demography, meds.model)")
endif()
meds_fortran_flags(meds_main)

#----- Tests. -----------------------------------------------------------------------------#
if(MEDS_BUILD_TESTS)
   enable_testing()
   #----- Assertions only, depending on nothing but a kind. Kept separate from meds_testsupport  #
   #      because many tests deliberately link ONE narrow library each to keep that layer             #
   #      standalone-buildable, and they must be able to assert without linking the whole model.   #
   add_library(meds_testassert STATIC test/meds_test_assert.f90 test/meds_test_soil_grid.f90)
   target_link_libraries(meds_testassert PUBLIC meds_share)
   meds_fortran_flags(meds_testassert)

   #----- The in-code config fixture; re-exports the assertions so most tests need one `use`. ----#
   add_library(meds_testsupport STATIC test/meds_test_support.f90)
   target_link_libraries(meds_testsupport PUBLIC meds_demography meds_testassert)
   meds_fortran_flags(meds_testsupport)


   #----- Shared numerical primitives + conservation checker (links meds_share only). ------#
   add_executable(test_numerics test/test_numerics.f90)
   target_link_libraries(test_numerics PRIVATE meds_share meds_testassert)
   meds_fortran_flags(test_numerics)
   add_test(NAME numerics COMMAND test_numerics)

   #----- Soil pedotransfer function: Hodnett & Tomasella (2002)'s worked example, Ks, coarse  #
   #      fragments (links meds_share only). ---------------------------------------------------#
   add_executable(test_pedotransfer test/test_pedotransfer.f90)
   target_link_libraries(test_pedotransfer PRIVATE meds_share meds_testassert)
   meds_fortran_flags(test_pedotransfer)
   add_test(NAME pedotransfer COMMAND test_pedotransfer)

   #----- No procedure is passed as an argument anywhere in src/: the ifx cost of doing        #
   #      so serialises threads. The script explains the rule it checks.                       #
   add_test(NAME no_procedure_arguments
            COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_no_procedure_arguments.cmake)

   #----- The two reference configs list exactly the keys the loader reads: a run stops on      #
   #      any key its reference does not list. The script explains the rule it checks.          #
   add_test(NAME config_keys_listed
            COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src
                    -DMAIN_TOML=${CMAKE_CURRENT_SOURCE_DIR}/meds_config_main.toml
                    -DPFT_TOML=${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_config_keys.cmake)
   #----- No bare literal clamps in src/: a guard uses tiny_num, a threshold is a named             #
   #      setting, or the line says why. scripts/lint/check_clamps.py explains the rule.             #
   find_package(Python3 COMPONENTS Interpreter QUIET)
   if (Python3_Interpreter_FOUND)
      add_test(NAME clamps_named
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/lint/check_clamps.py
                       ${CMAKE_CURRENT_SOURCE_DIR}/src)
      #----- No issue or pull-request numbers outside CHANGELOG.md, docs/ROADMAP.md and docs/dev_plans/:  #
      #      a comment says why the code is as it is; the history is in the CHANGELOG.                  #
      add_test(NAME no_issue_refs
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/lint/check_issue_refs.py
                       ${CMAKE_CURRENT_SOURCE_DIR})
      #----- The kernels take only the per-patch store types from the state layer, never polygon_t: #
      #      the shared module directory cannot enforce it, so this does, on every run.             #
      add_test(NAME kernels_no_site_state
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/lint/check_kernel_state.py
                       ${CMAKE_CURRENT_SOURCE_DIR}/src)
      #----- A named constant is spelled by its kind: lower case in meds_constants and meds_kinds,     #
      #      UPPER everywhere else. scripts/lint/check_constant_case.py explains the rule.             #
      add_test(NAME constants_named
               COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/lint/check_constant_case.py
                       ${CMAKE_CURRENT_SOURCE_DIR})
   endif()
   add_executable(test_config_keys test/test_config_keys.f90)
   target_link_libraries(test_config_keys PRIVATE meds_config meds_testsupport)
   meds_fortran_flags(test_config_keys)
   add_test(NAME config_keys COMMAND test_config_keys)
   add_test(NAME config_unknown_keys_refused
            COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main>
                    -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                    -DWORK_FILE=${CMAKE_CURRENT_BINARY_DIR}/config_unknown_keys.toml
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_config_refused.cmake)
   #----- A census or state file that cannot start the run stops it, naming the mode and the     #
   #      file: one that is not there, and a census with no row. The script explains the rule.    #
   foreach(case "census_missing;1;no_such_census.csv;OFF;cannot be opened"
                "census_empty;1;empty_census.csv;ON;holds no cohort for the site"
                "restart_missing;2;no_such_state.nc;OFF;cannot be opened")
      list(GET case 0 name)
      list(GET case 1 mode)
      list(GET case 2 file)
      list(GET case 3 empty)
      list(GET case 4 reason)
      add_test(NAME init_${name}_refused
               COMMAND ${CMAKE_COMMAND} -DMEDS_MAIN=$<TARGET_FILE:meds_main>
                       -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                       -DWORK_FILE=${CMAKE_CURRENT_BINARY_DIR}/init_${name}.toml
                       -DINIT_MODE=${mode} -DINIT_FILE=${CMAKE_CURRENT_BINARY_DIR}/${file}
                       -DEMPTY_CENSUS=${empty} "-DREASON=${reason}"
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_init_refused.cmake)
   endforeach()

   #----- The [pft] canopy optics reach the two-stream table, and differentiate PFTs. --------#
   add_executable(test_pft_optics_config test/test_pft_optics_config.f90)
   target_link_libraries(test_pft_optics_config PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_pft_optics_config)
   add_test(NAME pft_optics_config COMMAND test_pft_optics_config)

   #----- The [soil] column keys reach the column they describe. ----------------------------#
   add_executable(test_soil_column_config test/test_soil_column_config.f90)
   target_link_libraries(test_soil_column_config PRIVATE meds_config meds_testsupport)
   meds_fortran_flags(test_soil_column_config)
   add_test(NAME soil_column_config COMMAND test_soil_column_config)

   #----- Field-coverage of the fast state vector's combinators. ---------------------------#
   add_executable(test_state_combinators test/test_state_combinators.f90)
   target_link_libraries(test_state_combinators PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_state_combinators)
   add_test(NAME state_combinators COMMAND test_state_combinators)

   #----- The slow-loop site conservation ledger: store coverage + declaration arithmetic. ------#
   add_executable(test_slow_ledger test/test_slow_ledger.f90)
   target_link_libraries(test_slow_ledger PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_slow_ledger)
   add_test(NAME slow_ledger COMMAND test_slow_ledger)

   foreach(t containers fusion_cohort patch allometry disturbance init_census init_soil_state time carbon_growth fast_loop phenology_driver biogeochem_dynamics slow_diag_units mortality_pathways tissue_energy)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport meds_model)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Census files the reader must refuse: each passes on the reason it returns. ------#
   add_executable(test_init_census_refuses test/test_init_census_refuses.f90)
   target_link_libraries(test_init_census_refuses PRIVATE meds_testsupport meds_model)
   meds_fortran_flags(test_init_census_refuses)
   add_test(NAME init_census_refuses_area COMMAND test_init_census_refuses area)
   set_tests_properties(init_census_refuses_area PROPERTIES
                        PASS_REGULAR_EXPRESSION "two different patch_area values"
                        FAIL_REGULAR_EXPRESSION "NOT REFUSED")
   add_test(NAME init_census_refuses_column COMMAND test_init_census_refuses column)
   set_tests_properties(init_census_refuses_column PROPERTIES
                        PASS_REGULAR_EXPRESSION "unknown column \"basal\""
                        FAIL_REGULAR_EXPRESSION "NOT REFUSED")

   #----- Leaf-physiology test links meds_fast + test support (assertions). ---------------#
   add_executable(test_leaf_physiology test/test_leaf_physiology.f90)
   target_link_libraries(test_leaf_physiology PRIVATE meds_testsupport meds_fast)
   meds_fortran_flags(test_leaf_physiology)
   add_test(NAME leaf_physiology COMMAND test_leaf_physiology)

   #----- C-API shims. Compiling each C-API source INTO a ctest target is the point: a shim        #
   #      built only by the optional -DMEDS_BUILD_PYLIB=ON shared library breaks unseen -- a       #
   #      component inserted mid-type in leaf_photo_args_t, or a signature changed under the       #
   #      demography shim, leaves the whole suite green on both back ends while the C API is       #
   #      broken. Registering it here makes that class of break a normal build failure.            #
   #      Needs no .so and no Python: the test calls the bind(c) procedures directly. EVERY shim   #
   #      in src/c_api/ gets one, not just the leaf one.                                           #
   add_executable(test_c_api_leaf test/test_c_api_leaf.f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/meds_c_api_leaf.f90)
   target_link_libraries(test_c_api_leaf PRIVATE meds_testsupport meds_fast_kernels)
   meds_fortran_flags(test_c_api_leaf)
   add_test(NAME c_api_leaf COMMAND test_c_api_leaf)

   add_executable(test_c_api_phenology test/test_c_api_phenology.f90
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/meds_c_api_phenology.f90)
   target_link_libraries(test_c_api_phenology PRIVATE meds_testsupport meds_slow_kernels)
   meds_fortran_flags(test_c_api_phenology)
   add_test(NAME c_api_phenology COMMAND test_c_api_phenology)

   #      The demography shim needs the whole model (config + init + both loops), same as the .so.  #
   add_executable(test_c_api_demography test/test_c_api_demography.f90
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/meds_c_api_demography.f90)
   target_link_libraries(test_c_api_demography PRIVATE meds_testsupport meds_model meds_output_stream)
   meds_fortran_flags(test_c_api_demography)
   #      Runs from the SOURCE dir: it loads examples/example03_demography's shipped config, which
   #      gives that config a build-time consumer. A config nothing loads drifts behind the schema
   #      exactly like a shim nothing compiles.
   add_test(NAME c_api_demography COMMAND test_c_api_demography
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   #      The FULL-MODEL shim: it holds a run (meds_region), so it links the driver library (which
   #      is meds_model + IO plus the open/step/finalize seam) rather than the model directly.
   add_executable(test_c_api_run test/test_c_api_run.f90
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/c_api/meds_c_api_run.f90)
   target_link_libraries(test_c_api_run PRIVATE meds_testsupport meds_driver)
   meds_fortran_flags(test_c_api_run)
   #      Runs from examples/example04_column_biophysics: it loads both of that example's shipped
   #      configs (the default and the calibrated run), so the build notices when they fall behind
   #      the config reader.
   add_test(NAME c_api_run COMMAND test_c_api_run
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/example04_column_biophysics)

   #      A run's sites from [[sites]]: derived from example 04's config, it runs from that example's
   #      directory and writes under the build tree; meds_main checks what a config may not hold.
   add_executable(test_site_definitions test/test_site_definitions.f90)
   target_link_libraries(test_site_definitions PRIVATE meds_testsupport meds_driver)
   meds_fortran_flags(test_site_definitions)
   add_test(NAME site_definitions
            COMMAND test_site_definitions ${CMAKE_CURRENT_BINARY_DIR}/site_definitions $<TARGET_FILE:meds_main>
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples/example04_column_biophysics)

   #      A region against its site runs, on the synthetic archive.
   #      Runs from test/region, whose site config it derives its configs from; it writes
   #      everything under the work directory in the build tree.
   add_executable(test_region test/test_region.f90)
   target_link_libraries(test_region PRIVATE meds_driver meds_testassert meds_testarchive)
   meds_fortran_flags(test_region)
   add_test(NAME region COMMAND test_region ${CMAKE_CURRENT_BINARY_DIR}/region_test_tmp $<TARGET_FILE:meds_main>
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/region)

   #----- Hydraulics / phenology / respiration / carbon tests link the kernel libraries (state-free). #
   #     The list spans both tiers -- hydraulics and respiration are sub-daily, phenology,
   #     allocation and trait dynamics are daily -- so it links both kernel libraries: there is no
   #     edge between them to ride on.
   foreach(t plant_hydraulics plant_phenology plant_respiration plant_carbon_allocation plant_trait_dynamics)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_fast_kernels meds_slow_kernels meds_testassert)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Canopy radiative-transfer test links meds_fast_kernels (state-free). ------------#
   add_executable(test_canopy_radiation test/test_canopy_radiation.f90)
   target_link_libraries(test_canopy_radiation PRIVATE meds_fast_kernels meds_testassert)
   meds_fortran_flags(test_canopy_radiation)
   add_test(NAME canopy_radiation COMMAND test_canopy_radiation)

   #----- Canopy-aerodynamics test links meds_fast_kernels (state-free). ------------------#
   add_executable(test_aerodynamics test/test_aerodynamics.f90)
   target_link_libraries(test_aerodynamics PRIVATE meds_fast_kernels meds_testassert)
   meds_fortran_flags(test_aerodynamics)
   add_test(NAME aerodynamics COMMAND test_aerodynamics)

   #----- Vertical soil-hydrology test links meds_fast (it runs the soil-plant water solve). #
   add_executable(test_column_hydrology test/test_column_hydrology.f90)
   target_link_libraries(test_column_hydrology PRIVATE meds_fast meds_testassert)
   meds_fortran_flags(test_column_hydrology)
   add_test(NAME column_hydrology COMMAND test_column_hydrology)

   #----- The soil-plant water network's water relations and flows: values and derivatives. --#
   add_executable(test_water_processes test/test_water_processes.f90)
   target_link_libraries(test_water_processes PRIVATE meds_fast_kernels meds_testassert)
   meds_fortran_flags(test_water_processes)
   add_test(NAME water_processes COMMAND test_water_processes)

   #----- The soil-plant water network on synthetic stands and soils. ------------------------#
   add_executable(test_water_network test/test_water_network.f90)
   target_link_libraries(test_water_network PRIVATE meds_fast meds_testassert)
   meds_fortran_flags(test_water_network)
   add_test(NAME water_network COMMAND test_water_network)

   #----- Energy-balance + snow tests link meds_fast_kernels (state-free). ----------------#
   foreach(t column_energy surface_energy snow)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_fast_kernels meds_testassert)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Column CO2-balance test: the fast CAS CO2 twin lives in meds_fast_kernels.      #
   add_executable(test_column_co2 test/test_column_co2.f90)
   target_link_libraries(test_column_co2 PRIVATE meds_fast_kernels meds_slow_kernels meds_testassert)
   meds_fortran_flags(test_column_co2)
   add_test(NAME column_co2 COMMAND test_column_co2)

   #----- Slow soil-carbon matrix test links biogeochemistry; the fast/slow-seam check also    #
   #      exercises the fast CAS-CO2 Rh kernel, which meds_fast_kernels carries.                #
   add_executable(test_soil_biogeochem test/test_soil_biogeochem.f90)
   target_link_libraries(test_soil_biogeochem PRIVATE meds_fast_kernels meds_slow_kernels meds_testassert)
   meds_fortran_flags(test_soil_biogeochem)
   add_test(NAME soil_biogeochem COMMAND test_soil_biogeochem)

   #----- A site's forcing within its polygon: the sun and the beam on a slope, the lapse to the site. #
   add_executable(test_site_terrain test/test_site_terrain.f90)
   target_link_libraries(test_site_terrain PRIVATE meds_met meds_testassert)
   meds_fortran_flags(test_site_terrain)
   add_test(NAME site_terrain COMMAND test_site_terrain)

   #----- Meteorological-forcing reader test links the forcing library (writes+reads a NetCDF). #
   add_executable(test_met_driver test/test_met_driver.f90)
   target_link_libraries(test_met_driver PRIVATE meds_met meds_testassert)
   meds_fortran_flags(test_met_driver)
   add_test(NAME met_driver COMMAND test_met_driver
            ${CMAKE_CURRENT_SOURCE_DIR}/data/co2/co2_cmip7_global_annual_1000-2022.txt)
   #      The synthetic ED_ERA5land archive, shared by the reader test and the region test.
   add_library(meds_testarchive STATIC test/meds_test_era5land_archive.f90)
   target_link_libraries(meds_testarchive PUBLIC meds_met)
   meds_fortran_flags(meds_testarchive)
   add_executable(test_met_era5land test/test_met_era5land.f90)
   target_link_libraries(test_met_era5land PRIVATE meds_met meds_testassert meds_testarchive)
   meds_fortran_flags(test_met_era5land)
   add_test(NAME met_era5land COMMAND test_met_era5land)
   #----- The static-map reader and the cell mask, on the synthetic archive's grid. ----------#
   add_executable(test_surface_maps test/test_surface_maps.f90)
   target_link_libraries(test_surface_maps PRIVATE meds_maps meds_testassert meds_testarchive)
   meds_fortran_flags(test_surface_maps)
   add_test(NAME surface_maps COMMAND test_surface_maps ${CMAKE_CURRENT_BINARY_DIR}/surface_maps_tmp)
   #      The flux-tower contract of a MEDS forcing file (docs/science/forcing.md §1).
   add_executable(test_met_tower test/test_met_tower.f90)
   target_link_libraries(test_met_tower PRIVATE meds_met meds_testassert)
   meds_fortran_flags(test_met_tower)
   add_test(NAME met_tower COMMAND test_met_tower)
   #      The flux-tower forcing tool's own tests (scripts/prepare_flux_tower/tests), when the Python
   #      CMake finds has what the tool needs; otherwise they are left out, and configure says so.
   find_package(Python3 COMPONENTS Interpreter QUIET)
   if (Python3_Interpreter_FOUND)
      execute_process(COMMAND ${Python3_EXECUTABLE} -c
                      "import sys, numpy, pandas, netCDF4, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                      RESULT_VARIABLE meds_tower_tool_deps OUTPUT_QUIET ERROR_QUIET)
      if (meds_tower_tool_deps EQUAL 0)
         add_test(NAME prepare_flux_tower
                  COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                          ${CMAKE_CURRENT_SOURCE_DIR}/scripts/prepare_flux_tower/tests)
         set_tests_properties(prepare_flux_tower PROPERTIES ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
         #----- A restart continues exactly, and re-acclimates its traits on request                #
         #      (test/python/test_restart_exact.py; meds_main on the constant reference climate).  #
         add_test(NAME restart_exact
                  COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                          ${CMAKE_CURRENT_SOURCE_DIR}/test/python/test_restart_exact.py)
         set_tests_properties(restart_exact PROPERTIES
                              ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;MEDS_MAIN=$<TARGET_FILE:meds_main>")
         #----- The fast-parameter calibration tool (scripts/calibrate_fast/tests): the unit tests,  #
         #      a smoke test that runs a chain, a trial and a Jacobian column through meds_main (and, #
         #      when libmeds.so is built, through the Python API), and a short end-to-end fit.        #
         set(meds_calib_env "PYTHONDONTWRITEBYTECODE=1" "MEDS_MAIN=$<TARGET_FILE:meds_main>")
         if (MEDS_BUILD_PYLIB)
            list(APPEND meds_calib_env "MEDS_LIB=$<TARGET_FILE:meds_py>")
         endif()
         add_test(NAME calibrate_fast
                  COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                          ${CMAKE_CURRENT_SOURCE_DIR}/scripts/calibrate_fast/tests)
         set_tests_properties(calibrate_fast PROPERTIES ENVIRONMENT "${meds_calib_env}")
         #----- The Python package's own tests (python/tests): the leaf and phenology kernels,     #
         #      meds.config, and meds.model against meds_main, bit for bit.                         #
         if (MEDS_BUILD_PYLIB)
            add_test(NAME python_api
                     COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                             ${CMAKE_CURRENT_SOURCE_DIR}/python/tests)
            set_tests_properties(python_api PROPERTIES ENVIRONMENT
                                 "PYTHONDONTWRITEBYTECODE=1;PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python;MEDS_LIB=$<TARGET_FILE:meds_py>;MEDS_MAIN=$<TARGET_FILE:meds_main>")
         endif()
      else()
         message(STATUS "MEDS: prepare_flux_tower tests skipped (${Python3_EXECUTABLE} lacks numpy, pandas, "
                        "netCDF4, pytest or tomli)")
      endif()
      #----- The census preparation tool's tests (scripts/prepare_census/tests), on the same terms. -#
      execute_process(COMMAND ${Python3_EXECUTABLE} -c
                      "import sys, numpy, pandas, pytest; sys.version_info >= (3, 11) or __import__('tomli')"
                      RESULT_VARIABLE meds_census_tool_deps OUTPUT_QUIET ERROR_QUIET)
      if (meds_census_tool_deps EQUAL 0)
         add_test(NAME prepare_census
                  COMMAND ${Python3_EXECUTABLE} -m pytest -q -p no:cacheprovider
                          ${CMAKE_CURRENT_SOURCE_DIR}/scripts/prepare_census/tests)
         set_tests_properties(prepare_census PROPERTIES ENVIRONMENT PYTHONDONTWRITEBYTECODE=1)
      else()
         message(STATUS "MEDS: prepare_census tests skipped (${Python3_EXECUTABLE} lacks numpy, pandas, "
                        "pytest or tomli)")
      endif()
   endif()

   #----- Fast-loop biophysics run-config TOML wiring ([soil]/[fast]/[snow]/[aerodynamics]).   #
   add_executable(test_biophysics_params_config test/test_biophysics_params_config.f90)
   target_link_libraries(test_biophysics_params_config PRIVATE meds_config meds_testassert)
   meds_fortran_flags(test_biophysics_params_config)
   add_test(NAME biophysics_params_config COMMAND test_biophysics_params_config)

   #----- Fast-loop column-dynamics integration test (links the whole model + test support). #
   add_executable(test_column_dynamics test/test_column_dynamics.f90)
   target_link_libraries(test_column_dynamics PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_dynamics)
   add_test(NAME column_dynamics COMMAND test_column_dynamics)

   #----- Pure fast-loop RHS; links meds_model for meds_fast_time_derivs.                    #
   add_executable(test_column_derivs test/test_column_derivs.f90)
   target_link_libraries(test_column_derivs PRIVATE meds_model meds_testassert)
   meds_fortran_flags(test_column_derivs)
   add_test(NAME column_derivs COMMAND test_column_derivs)

   #----- INTEG_ARK wiring test: the coupled ARK (ESDIRK2) path through fast_step. ----------#
   add_executable(test_column_ark test/test_column_ark.f90)
   target_link_libraries(test_column_ark PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_ark)
   add_test(NAME column_ark COMMAND test_column_ark)

   #----- INTEG_RK45 wiring test: the ED2-faithful adaptive Cash-Karp RK45 path. --------------#
   add_executable(test_column_rk45 test/test_column_rk45.f90)
   target_link_libraries(test_column_rk45 PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_rk45)
   add_test(NAME column_rk45 COMMAND test_column_rk45)

   #----- Diagnostic-output netCDF-free core: integrator kernels + registry (state-free of C). #
   foreach(t output_integrate output_registry phenology_config)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_output meds_testsupport)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- The shipped [allometry] block must equal the meds_allometry_lib initializers (the     #
   #      Python library's config-free values); biomass must be Chave 2014 in carbon, leaf      #
   #      ED2's over C2B.                                                                       #
   add_executable(test_allometry_defaults test/test_allometry_defaults.f90)
   target_link_libraries(test_allometry_defaults PRIVATE meds_testsupport meds_model)
   meds_fortran_flags(test_allometry_defaults)
   add_test(NAME allometry_defaults COMMAND test_allometry_defaults
            ${CMAKE_CURRENT_SOURCE_DIR}/meds_config_pft.toml)

   #----- The shipped per-variable override EXAMPLE is a generated file that is also tracked, so #
   #      it can drift silently from the registry. The source tree path is passed in because     #
   #      ctest runs from the build tree.                                                        #
   add_executable(test_io_config_example test/test_io_config_example.f90)
   target_link_libraries(test_io_config_example PRIVATE meds_output meds_testsupport)
   meds_fortran_flags(test_io_config_example)
   add_test(NAME io_config_example COMMAND test_io_config_example
            ${CMAKE_CURRENT_SOURCE_DIR}/meds_io_config.toml)

   #----- The Python wheel ships COPIES of LICENSE and NOTICE, because python/pyproject.toml     #
   #      cannot list license files above its own directory. A tracked copy drifts silently,     #
   #      for the same reason as the io_config_example above, so the build compares them.        #
   add_test(NAME license_copy COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
            ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/python/LICENSE)
   add_test(NAME notice_copy COMMAND ${CMAKE_COMMAND} -E compare_files --ignore-eol
            ${CMAKE_CURRENT_SOURCE_DIR}/NOTICE ${CMAKE_CURRENT_SOURCE_DIR}/python/NOTICE)

   #----- Diagnostic-output serializer round-trip (writes + reopens netCDF; links the io layer). #
   add_executable(test_output_roundtrip test/test_output_roundtrip.f90)
   target_link_libraries(test_output_roundtrip PRIVATE meds_output_stream meds_testsupport)
   meds_fortran_flags(test_output_roundtrip)
   add_test(NAME output_roundtrip COMMAND test_output_roundtrip)

   #----- The WEIGHTED scale-aggregation layer: extensive/intensive, PFT + DBH-class closure,  #
   #      empty-set handling, class edges, and the flattened 2-D (layer, patch) soil slab.      #
   #      netCDF-free, so it links meds_output rather than the serializer.                     #
   add_executable(test_diagnostic_reduce test/test_diagnostic_reduce.f90)
   target_link_libraries(test_diagnostic_reduce PRIVATE meds_output meds_testsupport)
   meds_fortran_flags(test_diagnostic_reduce)
   add_test(NAME diagnostic_reduce COMMAND test_diagnostic_reduce)

   #----- State restart round-trip (plastic-trait persistence) links the IO layer + model. #
   add_executable(test_state_roundtrip test/test_state_roundtrip.f90)
   target_link_libraries(test_state_roundtrip PRIVATE meds_model meds_output_stream meds_testsupport)
   meds_fortran_flags(test_state_roundtrip)
   add_test(NAME state_roundtrip COMMAND test_state_roundtrip)

   #----- A polygon of several sites: its stands, threads and state file (links the IO layer). #
   add_executable(test_several_sites test/test_several_sites.f90)
   target_link_libraries(test_several_sites PRIVATE meds_model meds_output_stream meds_testsupport)
   meds_fortran_flags(test_several_sites)
   add_test(NAME several_sites COMMAND test_several_sites)
endif()
