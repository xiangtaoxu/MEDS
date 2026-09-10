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
#  four nested frames (column_fast_step -> build_column_frozen -> column_prepass ->
#  aero_bottom_to_top). On an 8 MB stack that is a hard SIGSEGV ceiling near 10-12k cohorts, and it
#  scales DOWN with OMP_STACKSIZE once the patch loop is threaded.
#
#  PREFERRED FIX IS A RUNTIME ONE -- raise the stack rather than change codegen:
#      ulimit -s unlimited          (or e.g. `ulimit -s 65536` for ~90k cohorts)
#      export OMP_STACKSIZE=512M    (once the patch loop is threaded)
#
#  -DMEDS_HEAP_ARRAYS=ON moves large automatic arrays to the heap instead. It is OFF by default
#  because ifx's -heap-arrays SEGFAULTS test_output_registry at the exit of apply_variable_override
#  in a Debug build -- a compiler codegen interaction, not a MEDS bug (the routine only assigns
#  plain integer components, and the same code is clean without the flag). Turn it on only if the
#  runtime fix is unavailable, and re-run the Debug suite when you do.
#-----------------------------------------------------------------------------------------#
#
#  AUTOMATIC (STACK) LOCALS ARE MANDATORY ONCE THE PATCH LOOP IS THREADED -- section 7 C4.
#
#  This is the trap that cost the most to find, so it is written down. Intel Fortran defaults to
#  `-auto-scalar`: local SCALARS go on the stack, but local ARRAYS and DERIVED TYPES go in STATIC
#  storage. Static locals are shared by every thread. The fast loop calls kernels in meds_fast_kernels,
#  meds_slow_kernels, meds_shared and meds_demography, and every one with a local array or a local
#  derived type (e.g. `rad_forcing_t rf` in apply_rt_forcing -> alloc_rad_forcing) then has all
#  threads writing the same object. It does NOT crash: it produces plausible, silently
#  thread-count-dependent answers -- precisely what section 7's byte-identity requirement exists to
#  catch. ThreadSanitizer finds it; a green test suite does not.
#
#  `-qopenmp`/`-fopenmp` implies automatic locals, but only for the target it is applied to, so
#  flagging just the target that owns meds_fast_dynamics is NOT enough -- the whole call tree needs
#  it. MEDS_OPENMP therefore adds the per-compiler "all locals on the stack" flag to EVERY target.
#  Corollary: more goes on the stack, so the OMP_STACKSIZE note above applies with more force.
#-----------------------------------------------------------------------------------------#
option(MEDS_HEAP_ARRAYS "Put large automatic arrays on the heap (see the note above)" OFF)
option(MEDS_OPENMP "Compile the fast-loop patch parallelism (plan section 7)" OFF)
#      Any build that can run host threads needs automatic locals -- which includes NVHPC with
#      MEDS_GPU != none, because that puts -mp PUBLIC on meds_demography and meds_fast inherits it, so
#      patch-loop directives are live there whether or not MEDS_OPENMP was asked for.
set(MEDS_AUTO_LOCALS OFF)
if(MEDS_OPENMP OR NOT MEDS_GPU STREQUAL "none")
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

#----- Foundation: precision, constants, PFT trait table, config, time, temperature response. #
#      No heavy deps. meds_temp_response (Arrhenius / peaked deactivation) was promoted here     #
#      from the leaf module so leaf, respiration, and any other tissue share one code path        #
#      without a plant-domain -> plant-domain library edge.                                        #
#      `src/shared/` is GONE: decision #3 dissolved it into base/ + functions/ + config/ + state/, #
#      and the last three folders were still wearing it as a prefix. They are top-level now, so    #
#      "shared" is not a place anything can be filed for lack of a better one. The library target  #
#      keeps the name meds_shared: it IS the shared foundation, and that is a fact about the link  #
#      graph rather than a folder for leftovers.                                                    #
#      glob; module names + .mod resolution unchanged). meds_allometry (pure size relationships)  #
#      now lives in functions/ alongside meds_temp_response / meds_therm_lib -- a shared geometric  #
#      CONSTRAINT reachable by state, demography, plant, and config_io without a separate library. #
file(GLOB SHARED_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.f90
                                          ${CMAKE_CURRENT_SOURCE_DIR}/src/functions/*.f90
                                          ${CMAKE_CURRENT_SOURCE_DIR}/src/util/*.f90)
add_library(meds_shared STATIC ${SHARED_SOURCES})
meds_fortran_flags(meds_shared)

#----- CONFIGURATION: the PFT trait table, the per-domain option leaves (soil/energy/snow/aero/  #
#      decomp/output/forcing), the meds_config_t aggregator, the minimal TOML reader and the      #
#      loader that fills the aggregator from two TOML files. `toml` + `config_io` moved here out   #
#      of src/io/ (plan decision #9): they are configuration INPUT, not I/O of model results, and   #
#      they carry no external dependency -- unlike the netCDF output layer they used to sit beside. #
#      Depends on meds_shared only; nothing here may import state or a kernel (checked: it does not).#
file(GLOB CONFIG_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/config/*.f90)
add_library(meds_config STATIC ${CONFIG_SOURCES})
#      meds_state_column: meds_config takes `n_soil_layer_max` from meds_column_params for the
#      [soil_column] block's bounds check. The edge is acyclic (state/column links meds_shared
#      only) and was MISSING until a Ninja build exposed it -- the Makefile generator happened to
#      order the two .mod files correctly, so a parallel build with a different scheduler failed
#      where the everyday one passed.
target_link_libraries(meds_config PUBLIC meds_shared meds_state_column)
meds_fortran_flags(meds_config)

#----- STATE, as a LAYER in two halves (MEDS_CODE_STRUCTURE_DESIGN.md decision #2). ---------#
#                                                                                            #
#  state/column -- ONE patch's reservoirs, seen vertically: the soil water / soil energy /    #
#  snow / canopy-air-space / soil-carbon stores, their fusion blends, and the parameter        #
#  bundles that describe them. This is BOUNDARY state: the kernels mutate it and the           #
#  demographic hub owns it. It used to sit in `shared/state`, which by its own docstring is    #
#  "the foundation, NOT tied to any process" -- a workaround for the fact that `shared` was     #
#  the only node both sides could see. It is a layer now, so the kernels link it directly.      #
#                                                                                            #
#  state/site -- ALL the patches: the flat cohort SoA, the patch CSR, the lockstep machinery,   #
#  site_t and the diagnostic blocks. Drivers only. Kernels never link this, which is what keeps #
#  them site_t-free and OpenMP-`target` device-eligible, and what lets the 12 kernel-only        #
#  ctests keep linking one kernel library alone.                                                 #
file(GLOB STATE_COLUMN_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/state/column/*.f90)
add_library(meds_state_column STATIC ${STATE_COLUMN_SOURCES})
target_link_libraries(meds_state_column PUBLIC meds_shared)
meds_fortran_flags(meds_state_column)

file(GLOB STATE_SITE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/state/site/*.f90)
add_library(meds_state_site STATIC ${STATE_SITE_SOURCES})
target_link_libraries(meds_state_site PUBLIC meds_state_column meds_config)
meds_fortran_flags(meds_state_site)



#==========================================================================================#
#  PROCESSES, TIMESCALE-FIRST (MEDS_CODE_STRUCTURE_DESIGN.md decision #4).                    #
#                                                                                            #
#  MEDS is an operator-split fast/slow model -- that is its dominant organizing axis -- and    #
#  the tree now says so. Each half holds its domain kernels, and the fast half additionally     #
#  holds the integrator machinery (numerics/) and the loop that drives it (driver/).            #
#==========================================================================================#

#----- SLOW kernels (daily / annual): phenology, carbon allocation, trait dynamics, the CENTURY #
#      soil-carbon matrix and the litter-destination split. Stateless and `pure`-where-possible; #
#      links the column state it reads/writes plus the config leaves. Builds standalone.         #
file(GLOB SLOW_KERNEL_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/plant/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/soil/*.f90)
add_library(meds_slow_kernels STATIC ${SLOW_KERNEL_SOURCES})
target_link_libraries(meds_slow_kernels PUBLIC meds_state_column meds_config)
meds_fortran_flags(meds_slow_kernels)

#----- FAST kernels (sub-daily): canopy radiative transfer + aerodynamics + the canopy air space #
#      (canopy/ = the medium), leaf gas exchange + hydraulics + maintenance respiration + the     #
#      tissue energy balance (plant/ = the organisms), and the soil water / soil energy / ground  #
#      + snow stores (soil/). State-free apart from the column reservoirs; builds standalone.      #
#                                                                                                  #
#      It links meds_state_column and meds_config only. The meds_slow_kernels edge it used to carry  #
#      was an artefact of meds_plant_interface re-exporting phenology and allocation; deleting that  #
#      facade (step 9) removed it. The fast kernels now import nothing from slow_dynamics/, and the  #
#      one real fast->slow edge (the S1 Rh seam) sits where it belongs, on meds_fast.                #
file(GLOB FAST_KERNEL_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/canopy/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/plant/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/soil/*.f90)
add_library(meds_fast_kernels STATIC ${FAST_KERNEL_SOURCES})
target_link_libraries(meds_fast_kernels PUBLIC meds_state_column meds_config)
meds_fortran_flags(meds_fast_kernels)

#----- DEMOGRAPHY: the vital-rate laws (growth / mortality / recruitment) and the apply-operators #
#      that change the population -- state update, cohort and patch fuse-fission, termination,     #
#      recruitment, treefall disturbance -- plus the public facade. Links the site state it         #
#      mutates and the litter law its kills feed. Placement rule 8 keeps the two halves apart:      #
#      the operators take rate ARRAYS and never import meds_demography_rates.                        #
file(GLOB DEMOGRAPHY_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/demography/*.f90)
add_library(meds_demography STATIC ${DEMOGRAPHY_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_state_site meds_slow_kernels)
meds_fortran_flags(meds_demography)

#----- Output PREP (src/io -- netCDF + diagnostics only now that toml/config_io have left): the netCDF-FREE half of the output subsystem -- the site-level        #
#      reductions (meds_diagnostic_{kernels,reduce}) + variable registry + per-tier temporal integrators +   #
#      the extract switchboard. Links meds_demography (site_t + the reductions) ONLY -- NO netCDF --  #
#      the stepper edge that ticks the integrators pulls no C dependency (the §2 DAG-hygiene wall).  #
#      The netCDF serializer half lives in meds_io_stream (below). Design: docs/dev_plans/MEDS_IO_DESIGN.md. #
#      NOTE: this is an EXPLICIT file list, not a file(GLOB) -- adding a diagnostic module here is  #
#      a required manual edit (unlike src/plant / src/biophysics, which do glob).                    #
add_library(meds_io_prep STATIC
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_diagnostic_kernels.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_diagnostic_reduce.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_types.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_integrate.f90
   ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_registry.f90)
target_link_libraries(meds_io_prep PUBLIC meds_demography)
meds_fortran_flags(meds_io_prep)





#----- OpenMP-target offload of the hot kernel (NVHPC; it lives in meds_demography). -------#
#       PUBLIC so the demo/tests/drivers that link meds_demography inherit the offload flags.          #
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


#----- netCDF C bindings, PROMOTED to their own minimal target so BOTH the output layer         #
#      (meds_io, state serialization) AND the forcing reader (meds_forcing) link the netCDF C     #
#      API WITHOUT meds_io_stream's -> meds_demography edge. netCDF is a hard dependency;            #
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
target_link_libraries(meds_forcing PUBLIC meds_config meds_netcdf_c)
meds_fortran_flags(meds_forcing)


#----- The FAST loop: the integrator machinery (numerics/ -- the state vector, the frozen work  #
#      records, the ARK / RK45 / split marches, the step controller and the pre-pass) and the    #
#      driver that walks a slow step in dt_fast sub-steps over the patch axis (driver/).          #
#                                                                                                  #
#      meds_fast_prepass sits in numerics/, not driver/, because the marches call it: ark, rk45     #
#      and the oracle all invoke column_prepass mid-stage. The plan filed it under driver/ on the   #
#      strength of its name; the call graph says otherwise, and a cycle is not a naming question.    #
#      It links meds_slow_kernels for the documented fast->slow seam (S1): the sub-daily            #
#      heterotrophic respiration must respire the SAME CENTURY matrix the daily step debits.        #
file(GLOB FAST_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/numerics/*.f90
     ${CMAKE_CURRENT_SOURCE_DIR}/src/fast_dynamics/driver/*.f90)
add_library(meds_fast STATIC ${FAST_SOURCES})
target_link_libraries(meds_fast PUBLIC meds_fast_kernels meds_slow_kernels meds_state_site
                                       meds_forcing meds_io_prep)
meds_fortran_flags(meds_fast)

#----- The SLOW loop: the thin coordinator and the two peer domain drivers (vegetation dynamics, #
#      soil biogeochemistry). It is the ONE place a demographic rate meets its application.       #
file(GLOB SLOW_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/slow_dynamics/driver/*.f90)
add_library(meds_slow STATIC ${SLOW_SOURCES})
target_link_libraries(meds_slow PUBLIC meds_demography meds_slow_kernels meds_fast_kernels meds_io_prep)
meds_fortran_flags(meds_slow)

#----- Initial-community builders (bare ground, add_cohort, census). ------------------------#
file(GLOB INIT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/init/*.f90)
add_library(meds_init STATIC ${INIT_SOURCES})
target_link_libraries(meds_init PUBLIC meds_demography meds_config)
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

#----- Section 7 C4: HOST OpenMP over the fast-loop PATCH axis (src/driver/meds_fast_dynamics).  #
#      OFF by default, so no existing build changes. Applied to meds_fast -- the target that owns  #
#      the fast loop -- and PUBLIC so meds_main and the tests inherit the link flags.               #
#      find_package resolves the per-compiler flag (-qopenmp IntelLLVM, -fopenmp GNU, -mp NVHPC).   #
#      NVHPC with MEDS_GPU != none ALREADY carries -mp PUBLIC on meds_demography, inherited here,   #
#      so adding it again would fight -mp=gpu; that case is left alone. Note the consequence: such a #
#      build compiles the patch-loop directives whether or not MEDS_OPENMP is set, which is exactly  #
#      why the thread count is a RUNTIME knob ([run].n_threads, default 1) and not a build-time one. #
#      The option itself is declared beside meds_fortran_flags(), which is where its OTHER (and more #
#      important) half lives: the automatic-locals flag on every target.                             #
if(MEDS_OPENMP)
   if(NOT MEDS_GPU STREQUAL "none" AND CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      message(STATUS "MEDS: fast-loop OpenMP inherited from MEDS_GPU=${MEDS_GPU} (-mp is already PUBLIC on meds_demography)")
   else()
      find_package(OpenMP REQUIRED COMPONENTS Fortran)
      target_link_libraries(meds_fast PUBLIC OpenMP::OpenMP_Fortran)
      message(STATUS "MEDS: fast-loop OpenMP -> ${OpenMP_Fortran_FLAGS} (patch axis; [run].n_threads)")
   endif()
endif()

#----- Output layer (src/io): the netCDF state/diagnostic writer. Links the promoted netCDF        #
#      bindings (meds_netcdf_c, above) + demography (state serialization). netCDF is always         #
#      compiled. (The C-binding layer passes array temporaries to the netCDF C `(*)` dummies --      #
#      correct for interop, but Intel Debug `-check all` emits a harmless arg_temp_created warning   #
#      per call; build Release (-DCMAKE_BUILD_TYPE=Release) for clean output.)                       #
add_library(meds_io_stream STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_io.f90
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_stream.f90
                           ${CMAKE_CURRENT_SOURCE_DIR}/src/io/meds_output_manager.f90)
target_link_libraries(meds_io_stream PUBLIC meds_demography meds_netcdf_c meds_io_prep)
meds_fortran_flags(meds_io_stream)
set(MEDS_IO_LIB meds_io_stream)
message(STATUS "MEDS: netCDF output (always compiled)")

#----- The main entry point: read config, run the simulation, save output, exit. -----------#
add_executable(meds_main ${CMAKE_CURRENT_SOURCE_DIR}/src/main/meds_main.f90)
target_link_libraries(meds_main PRIVATE meds_model ${MEDS_IO_LIB})

#----- THE Python C-API shared library (opt-in): ONE libmeds.so for the whole package.        #
#                                                                                              #
#      Structure-plan decision #1. There used to be two -- libmeds_plant_c (leaf + phenology,   #
#      kernels only) and libmeds_c (the site driver) -- and that split was the ONLY reason for  #
#      the "each domain kernel library links meds_shared ONLY" rule, which in turn forced the   #
#      column state types into a shared/state folder. One .so dissolves it: per-subsystem       #
#      independence moves from the LINK level to the CALL level (verbs, not libraries), and     #
#      `meds.plant.leaf`, `meds.plant.pheno` and `meds.demography` are sub-modules over one     #
#      backend. It links meds_model + IO, so the leaf shim now carries the whole model; that is  #
#      the accepted cost of #1, and decision #11 is its consequence (netcdf4 is declared a       #
#      runtime dependency of the wheel rather than vendored).                                    #
#                                                                                              #
#      Plain GLOB, no exclusions: every src/capi/*.f90 belongs in the one library now.           #
if(MEDS_BUILD_PYLIB)
   file(GLOB CAPI_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/capi/*.f90)
   add_library(meds_py SHARED ${CAPI_SOURCES})
   set_target_properties(meds_py PROPERTIES OUTPUT_NAME meds)
   target_link_libraries(meds_py PRIVATE meds_model ${MEDS_IO_LIB})
   meds_fortran_flags(meds_py)
   #----- RPATH so the INSTALLED library resolves its own dependencies.                         #
   #                                                                                            #
   #      libmeds.so links netCDF, HDF5 and the Fortran runtime. Decision #11 says declare       #
   #      netcdf4 rather than vendor the ~67 MB closure -- necessary, but on its own NOT enough:  #
   #      the pip netCDF4 wheel bundles its libnetcdf under a mangled soname inside its own        #
   #      package, so it cannot satisfy a `libnetcdf.so.19` NEEDED by a different .so. Freshly     #
   #      installed into a bare venv, this library failed to load on ten unresolved sonames.       #
   #                                                                                            #
   #      USE_LINK_PATH writes the directories of the actually-linked externals into the install  #
   #      RPATH, so the wheel finds them with no LD_LIBRARY_PATH and no `source setvars.sh`.      #
   #      $ORIGIN comes first so a co-bundled copy would still win.                                #
   #                                                                                            #
   #      THIS MAKES THE WHEEL MACHINE-LOCAL: those paths are the build machine's conda prefix    #
   #      and oneAPI install. That is the right trade for `pip install python/` from a source      #
   #      checkout, which is how MEDS is used; it is NOT a redistributable manylinux wheel, and    #
   #      building one would mean auditwheel plus the vendoring §7.3 rejected.                     #
   #      USE_LINK_PATH alone is NOT enough: it captures the directories of libraries CMake was  #
   #      told to link (netCDF, HDF5), but the Fortran runtime -- libimf, libifport, libsvml --   #
   #      is added implicitly by the compiler driver, so it is absent from the RPATH and the      #
   #      library still fails to dlopen outside a `source setvars.sh` shell. CMAKE_Fortran_       #
   #      IMPLICIT_LINK_DIRECTORIES is exactly that set. (`ldd` will NOT show this: run from a    #
   #      sourced shell it inherits LD_LIBRARY_PATH and reports everything resolved. Check with   #
   #      `env -u LD_LIBRARY_PATH ldd`, or by importing.)                                          #
   set_target_properties(meds_py PROPERTIES
                         INSTALL_RPATH "$ORIGIN;${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}"
                         INSTALL_RPATH_USE_LINK_PATH TRUE)
   #----- INSTALL it INTO the package directory, so a wheel carries meds/libmeds.so and         #
   #      meds._libmeds finds it as a sibling of __init__.py -- the first entry in its search    #
   #      order, needing no env var. The `pylib` component is what python/pyproject.toml         #
   #      installs; nothing else in this project is part of the wheel.                            #
   install(TARGETS meds_py LIBRARY DESTINATION meds COMPONENT pylib)
   message(STATUS "MEDS: Python C-API ENABLED -> libmeds.so (leaf + phenology + demography)")
endif()
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

   #----- The [pft] canopy optics reach the two-stream table, and differentiate PFTs. --------#
   add_executable(test_pft_optics_config test/test_pft_optics_config.f90)
   target_link_libraries(test_pft_optics_config PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_pft_optics_config)
   add_test(NAME pft_optics_config COMMAND test_pft_optics_config)

   #----- The [soil_column] keys reach the column they describe. -----------------------------#
   add_executable(test_soil_column_config test/test_soil_column_config.f90)
   target_link_libraries(test_soil_column_config PRIVATE meds_config meds_state_column meds_testsupport)
   meds_fortran_flags(test_soil_column_config)
   add_test(NAME soil_column_config COMMAND test_soil_column_config)

   #----- Field-coverage of the fast state vector's combinators (structure plan §10.3). ------#
   add_executable(test_state_combinators test/test_state_combinators.f90)
   target_link_libraries(test_state_combinators PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_state_combinators)
   add_test(NAME state_combinators COMMAND test_state_combinators)

   #----- The slow-loop site conservation ledger: store coverage + declaration arithmetic. ------#
   add_executable(test_slow_ledger test/test_slow_ledger.f90)
   target_link_libraries(test_slow_ledger PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_slow_ledger)
   add_test(NAME slow_ledger COMMAND test_slow_ledger)

   foreach(t containers fusion_cohort patch allometry disturbance init_census time carbon_growth fast_loop phenology_driver biogeochem_dynamics)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_testsupport meds_model)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Leaf-physiology test links the flat plant library + test support (assertions). ---#
   add_executable(test_leaf_physiology test/test_leaf_physiology.f90)
   target_link_libraries(test_leaf_physiology PRIVATE meds_testsupport meds_fast)
   meds_fortran_flags(test_leaf_physiology)
   add_test(NAME leaf_physiology COMMAND test_leaf_physiology)

   #----- C-API shim (issue #100). Compiling the capi source INTO a ctest target is the point: the #
   #      shim used to be built only by the optional -DMEDS_BUILD_PYLIB=ON shared library, so a      #
   #      component inserted mid-type in leaf_photo_params_t could break the C API while the whole    #
   #      suite stayed green on both back ends -- which is exactly what happened (#95 -> #100).       #
   #      Registering it here makes that class of break a normal build failure. Needs no .so and no   #
   #      Python: the test calls the bind(c) procedures directly.                                     #
   #      EVERY shim in src/capi/ gets one (§7.6 #4), not just the leaf one. The demography shim   #
   #      had none, and it did not compile against main at all: PR #137 changed a signature under   #
   #      it and 42/42 stayed green because only the optional .so built it. That is the same hole,  #
   #      one subsystem over.                                                                        #
   add_executable(test_capi_leaf test/test_capi_leaf.f90 ${CMAKE_CURRENT_SOURCE_DIR}/src/capi/meds_capi_leaf.f90)
   target_link_libraries(test_capi_leaf PRIVATE meds_testsupport meds_fast_kernels)
   meds_fortran_flags(test_capi_leaf)
   add_test(NAME capi_leaf COMMAND test_capi_leaf)

   add_executable(test_capi_phenology test/test_capi_phenology.f90
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/capi/meds_capi_phenology.f90)
   target_link_libraries(test_capi_phenology PRIVATE meds_testsupport meds_slow_kernels)
   meds_fortran_flags(test_capi_phenology)
   add_test(NAME capi_phenology COMMAND test_capi_phenology)

   #      The demography shim needs the whole model (config + init + both loops), same as the .so.  #
   add_executable(test_capi_demography test/test_capi_demography.f90
                  ${CMAKE_CURRENT_SOURCE_DIR}/src/capi/meds_capi_demography.f90)
   target_link_libraries(test_capi_demography PRIVATE meds_testsupport meds_model ${MEDS_IO_LIB})
   meds_fortran_flags(test_capi_demography)
   #      Runs from the SOURCE dir: it loads examples/example_demography's shipped config, which
   #      gives that config a build-time consumer. It had none, and had drifted 37 required keys
   #      behind the schema -- a config nothing loads rots exactly like a shim nothing compiles.
   add_test(NAME capi_demography COMMAND test_capi_demography
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   #----- Hydraulics / phenology / respiration / carbon tests link the flat plant library (state-free). #
   #     The list spans both tiers -- hydraulics and respiration are sub-daily, phenology,
   #     allocation and trait dynamics are daily -- so it links both kernel libraries. Since
   #     step 9 there is no edge between them to ride on.
   foreach(t plant_hydraulics plant_phenology plant_respiration plant_carbon_allocation plant_trait_dynamics)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_fast_kernels meds_slow_kernels)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Canopy radiative-transfer test links the biophysics library (state-free). --------#
   add_executable(test_canopy_radiation test/test_canopy_radiation.f90)
   target_link_libraries(test_canopy_radiation PRIVATE meds_fast_kernels)
   meds_fortran_flags(test_canopy_radiation)
   add_test(NAME canopy_radiation COMMAND test_canopy_radiation)

   #----- Canopy-aerodynamics test links the biophysics library (state-free). --------------#
   add_executable(test_aerodynamics test/test_aerodynamics.f90)
   target_link_libraries(test_aerodynamics PRIVATE meds_fast_kernels)
   meds_fortran_flags(test_aerodynamics)
   add_test(NAME aerodynamics COMMAND test_aerodynamics)

   #----- Vertical soil-hydrology test links the biophysics library (state-free). ----------#
   add_executable(test_column_hydrology test/test_column_hydrology.f90)
   target_link_libraries(test_column_hydrology PRIVATE meds_fast_kernels)
   meds_fortran_flags(test_column_hydrology)
   add_test(NAME column_hydrology COMMAND test_column_hydrology)

   #----- Energy-balance + snow tests link the biophysics library (state-free). ------------#
   foreach(t column_energy surface_energy snow)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_fast_kernels)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Column CO2-balance test: the fast CAS CO2 twin now lives in the biophysics library. #
   add_executable(test_column_co2 test/test_column_co2.f90)
   target_link_libraries(test_column_co2 PRIVATE meds_fast_kernels meds_slow_kernels)
   meds_fortran_flags(test_column_co2)
   add_test(NAME column_co2 COMMAND test_column_co2)

   #----- Slow soil-carbon matrix test links biogeochemistry; the fast/slow-seam check also    #
   #      exercises the fast CAS-CO2 Rh kernel, which meds_fast_kernels carries.                #
   add_executable(test_soil_biogeochem test/test_soil_biogeochem.f90)
   target_link_libraries(test_soil_biogeochem PRIVATE meds_fast_kernels meds_slow_kernels)
   meds_fortran_flags(test_soil_biogeochem)
   add_test(NAME soil_biogeochem COMMAND test_soil_biogeochem)

   #----- Meteorological-forcing reader test links the forcing library (writes+reads a NetCDF). #
   add_executable(test_met_driver test/test_met_driver.f90)
   target_link_libraries(test_met_driver PRIVATE meds_forcing)
   meds_fortran_flags(test_met_driver)
   add_test(NAME met_driver COMMAND test_met_driver)

   #----- Fast-loop biophysics run-config TOML wiring ([soil]/[energy]/[snow]/[aerodynamics]). #
   add_executable(test_biophysics_opts_config test/test_biophysics_opts_config.f90)
   target_link_libraries(test_biophysics_opts_config PRIVATE meds_config)
   meds_fortran_flags(test_biophysics_opts_config)
   add_test(NAME biophysics_opts_config COMMAND test_biophysics_opts_config)

   #----- Fast-loop column-dynamics integration test (links the aux driver + test support). #
   add_executable(test_column_dynamics test/test_column_dynamics.f90)
   target_link_libraries(test_column_dynamics PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_dynamics)
   add_test(NAME column_dynamics COMMAND test_column_dynamics)

   #----- Pure fast-loop RHS (IMEX-ARK P0); links meds_model for meds_fast_time_derivs.         #
   add_executable(test_column_derivs test/test_column_derivs.f90)
   target_link_libraries(test_column_derivs PRIVATE meds_model)
   meds_fortran_flags(test_column_derivs)
   add_test(NAME column_derivs COMMAND test_column_derivs)

   #----- INTEG_ARK wiring test: the coupled IMEX-ARK path through column_fast_step. -----------#
   add_executable(test_column_ark test/test_column_ark.f90)
   target_link_libraries(test_column_ark PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_ark)
   add_test(NAME column_ark COMMAND test_column_ark)

   #----- INTEG_RK4 wiring test: the ED2-faithful adaptive Cash-Karp RK45 path. -----------------#
   add_executable(test_column_rk45 test/test_column_rk45.f90)
   target_link_libraries(test_column_rk45 PRIVATE meds_model meds_testsupport)
   meds_fortran_flags(test_column_rk45)
   add_test(NAME column_rk45 COMMAND test_column_rk45)

   #----- Diagnostic-output netCDF-free core: integrator kernels + registry (state-free of C). #
   foreach(t output_integrate output_registry)
      add_executable(test_${t} test/test_${t}.f90)
      target_link_libraries(test_${t} PRIVATE meds_io_prep meds_testsupport)
      meds_fortran_flags(test_${t})
      add_test(NAME ${t} COMMAND test_${t})
   endforeach()

   #----- Diagnostic-output serializer round-trip (writes + reopens netCDF; links the io layer). #
   add_executable(test_output_roundtrip test/test_output_roundtrip.f90)
   target_link_libraries(test_output_roundtrip PRIVATE meds_io_stream meds_testsupport)
   meds_fortran_flags(test_output_roundtrip)
   add_test(NAME output_roundtrip COMMAND test_output_roundtrip)

   #----- The WEIGHTED scale-aggregation layer: extensive/intensive, PFT + DBH-class closure,  #
   #      empty-set handling, class edges, and the flattened 2-D (layer, patch) soil slab.      #
   #      netCDF-free, so it links meds_io_prep rather than the serializer.                     #
   add_executable(test_diagnostic_reduce test/test_diagnostic_reduce.f90)
   target_link_libraries(test_diagnostic_reduce PRIVATE meds_io_prep meds_testsupport)
   meds_fortran_flags(test_diagnostic_reduce)
   add_test(NAME diagnostic_reduce COMMAND test_diagnostic_reduce)

   #----- State restart round-trip (plastic-trait persistence) links the full IO + aux libs. #
   add_executable(test_state_roundtrip test/test_state_roundtrip.f90)
   target_link_libraries(test_state_roundtrip PRIVATE meds_model meds_io_stream meds_testsupport)
   meds_fortran_flags(test_state_roundtrip)
   add_test(NAME state_roundtrip COMMAND test_state_roundtrip)
endif()
