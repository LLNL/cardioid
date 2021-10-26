
set(extra_definitions "" CACHE STRING "")
set(extra_includes "" CACHE STRING "")

if (ENABLE_CUDA)
   set(cuda_definitions -DUSE_CUDA -DUSE_NVTX)
   set(cuda cuda)
   set(cuda_runtime cuda_runtime)
   set(cuda_root_include ${CUDA_TOOLKIT_ROOT_DIR}/include)
   
   find_library(CUDA_NVTX_LIBRARY
      NAMES nvToolsExt nvTools nvtoolsext nvtools nvtx NVTX
      PATHS ${CUDA_TOOLKIT_ROOT_DIR}
      PATH_SUFFIXES "lib64" "common/lib64" "common/lib" "lib"
      DOC "Location of the CUDA Toolkit Extension (NVTX) library"
      NO_DEFAULT_PATH
   )
   mark_as_advanced(CUDA_NVTX_LIBRARY)
   find_library(CUDA_NVRTC_LIBRARY
      NAMES nvrtc NVRTC
      PATHS ${CUDA_TOOLKIT_ROOT_DIR}
      PATH_SUFFIXES "lib64" "common/lib64" "common/lib" "lib"
      DOC "Location of the CUDA NVidia Runtime Compilation library"
      NO_DEFAULT_PATH
   )
   mark_as_advanced(CUDA_NVRTC_LIBRARY)
else()
   set(cuda_definitions "")
   set(cuda "")
   set(cuda_runtime "")
   set(cuda_root_include "")
   set(CUDA_NVTX_LIBRARY "")
   set(CUDA_NVRTC_LIBRARY "")
   set(CUDA_CUDA_LIBRARY "")
endif()   

find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)

find_program(PERL perl)

if (BLAS_LIB AND LAPACK_LIB)
   set(lapack_definitions -DHAVE_LAPACK) 
   list(APPEND extra_libs ${LAPACK_LIB} ${BLAS_LIB})
else()
   set(lapack_definitions "")
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set (heart_transport_src
   SpaceAllocator.cc 
   lazy_array.hh
)

list(APPEND extra_libs ${lapack})

find_program(PERL perl)

blt_add_library(NAME heart_transport
                SOURCES ${heart_transport_src}
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${cuda_root_include}
                DEFINES ${cuda_definitions} ${lapack_definitions} ${extra_definitions}
                DEPENDS_ON ${cuda_runtime} openmp)

if (ENABLE_CUDA)
   blt_add_library(NAME heart_actual_cuda
                   SOURCES CUDADiffusion_ker.cu ECGSensor.cu
                   INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${cuda_root_include}
                   DEFINES ${cuda_definitions} ${lapack_definitions} ${extra_definitions}
                   DEPENDS_ON ${cuda} simUtil simdops openmp)
   set (heart_explicit_cuda_src
      CUDADiffusion.cc
      CUDADiffusion.hh
      ECGSensor.cc
      ECGSensor.hh
   )
   blt_add_library(NAME heart_explicit_cuda
                   SOURCES ${heart_explicit_cuda_src}
                   DEPENDS_ON heart_transport heart_actual_cuda ${cuda_runtime} openmp)
   set(heart_explicit_cuda heart_explicit_cuda)
else()
   set(heart_explicit_cuda "")
endif()

set(threading_src
    ThreadServer.cc
    ThreadUtils.cc
   )

blt_add_library(NAME threading
                SOURCES ${threading_src}
                DEPENDS_ON openmp
                )

set(profiling_src
    PerformanceTimers.cc
    unionOfStrings.cc
    )

blt_add_library(NAME profiling
                SOURCES ${profiling_src}
                DEPENDS_ON simUtil openmp 
                )

set (reaction_src
	#TT04_CellML_Reaction.cc
	#TT06_CellML_Reaction.cc 
	#TT06Dev_Reaction.cc
        Passive.cc
        #ReactionFHN.cc
        #sundnes_et_al_2016_FHN.cc
	#TT06Tau.cc
	#TT06Func.cc
        BetterTT06.cc
        #Grandi.cc
	#pade.cc
)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
        COMMAND ${PERL} ARGS registerBuiltinReactions.pl ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc ${reaction_src}
        MAIN_DEPENDENCY registerBuiltinReactions.pl
        DEPENDS ${reaction_src}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
 )


set(ode_gpu_aware_src
    ${reaction_src}
    Reaction.cc
    Reaction.hh
    ReactionManager.cc
    ReactionManager.hh
   ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
   Interpolation.cc
   reactionFactory.cc
)
blt_add_library(NAME ode_gpu_aware
                SOURCES ${ode_gpu_aware_src}
                DEPENDS_ON heart_transport simUtil simdops ode_cpu_only profiling threading ${CMAKE_DL_LIBS}
                           ${extra_libs}
                           ${CUDA_NVTX_LIBRARY} ${CUDA_NVRTC_LIBRARY} ${CUDA_CUDA_LIBRARY}
                           ${cuda} ${cuda_runtime} openmp mpi
                INCLUDES ${extra_includes} ${CMAKE_CURRENT_SOURCE_DIR}
                )


set(ode_cpu_only_src
    TT06GatesSimd2.c TT06GatesSimd.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c 
    svd.c
    )
blt_add_library(NAME ode_cpu_only
                SOURCES ${ode_cpu_only_src}
                DEPENDS_ON ${extra_libs} openmp
                INCLUDES ${extra_includes}
                )
    
set (heart_gpu_aware_src
   ActivationAndRecoverySensor.cc
   ActivationTimeSensor.cc
   BoxStimulus.cc
   BoxStimulus.hh
   DVThreshSensor.cc
   DataVoronoiCoarsening.cc
   Diffusion.hh
   GradientVoronoiCoarsening.cc
   GradientVoronoiCoarsening.hh
   HaloExchange.hh
   MaxDVSensor.cc
   MinMaxSensor.cc
   PointListSensor.cc
   PointStimulus.cc
   PointStimulus.hh
   simulationLoop.cc
   Simulate.cc
   Simulate.hh
   StateVariableSensor.cc
   Stimulus.hh
   TestStimulus.cc
   TestStimulus.hh
   checkpointIO.cc
   getRemoteCells.cc
	stringUtils.cc
	readCellList.cc
	PioHeaderData.cc
	Vector.cc SymmetricTensor.cc
	getUserInfo.cc
        GDLoadBalancer.cc BlockLoadBalancer.cc workBoundBalancer.cc
	pioBalancer.cc
	CommTable.cc
	BucketOfBits.cc
	stateLoader.cc
	readPioFile.cc
	AnatomyReader.cc
	Koradi.cc GridRouter.cc Grid3DStencil.cc writeCells.cc
	checkpointIO.cc
	DomainInfo.cc
	BoundingBox.cc
	initializeSimulate.cc
	initializeAnatomy.cc setConductivity.cc assignCellsToTasks.cc
	diffusionFactory.cc
	stimulusFactory.cc sensorFactory.cc
	FibreConductivity.cc JHUConductivity.cc
	DiffusionUtils.cc
	FGRDiffusion.cc
	FGRDiffusionOMP.cc
	FGRDiffusionThreads.cc
	FGRDiffusionOverlap.cc
	FGRDiffusionStrip.cc
	#OpenmpGpuRedblackDiffusion.cc
	#OpenmpGpuFlatDiffusion.cc
	FGRUtils.cc
	PeriodicPulse.cc RandomPulse.cc
	VoronoiCoarsening.cc
	CaAverageSensor.cc
        ProcBox.cc
   Version.cc
)
blt_add_library(NAME heart_gpu_aware
                SOURCES ${heart_gpu_aware_src}
                DEFINES ${cuda_definitions} ${lapack_definitions} ${extra_definitions}
                DEPENDS_ON ${heart_explicit_cuda} heart_transport ode_gpu_aware heart_cpu_only simUtil simdops 
                           ${CMAKE_DL_LIBS} ${extra_libs}
                           ${CUDA_NVTX_LIBRARY} ${CUDA_NVRTC_LIBRARY} ${CUDA_CUDA_LIBRARY}
                           ${cuda} ${cuda_runtime} openmp mpi
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${cuda_root_include} ${extra_includes}
                )

set (heart_cpu_only_src 
   spi_impl.c
   clooper.c workBound.c
)
blt_add_library(NAME heart_cpu_only
                SOURCES ${heart_cpu_only_src}
                DEFINES ${cuda_definitions} ${lapack_definitions} ${extra_definitions}
                DEPENDS_ON simUtil simdops ${extra_libs} openmp
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${cuda_root_include} ${extra_includes}
                )

blt_add_executable(NAME cardioid
                   SOURCES cardioid.cc
                   DEPENDS_ON heart_gpu_aware ${heart_explicit_cuda}
                              heart_gpu_aware ${heart_explicit_cuda}
                              ${cuda_runtime} openmp)

blt_add_executable(NAME singleCell
                   SOURCES singleCell.cc singleCellOptions.c
                   DEPENDS_ON ode_gpu_aware ${cuda} ${cuda_runtime} openmp)

if (LAPACK_LIB)
   blt_add_executable(NAME modifyAnatomyFile
                      SOURCES modifyAnatomyFile.cc
                      DEPENDS_ON heart_gpu_aware ${cuda_runtime} openmp)
endif ()

install(TARGETS cardioid singleCell
        RUNTIME DESTINATION bin
        )


### Getting compile flags for making plugins.

function(get_recursive_libs retvar target)
   set(${retvar} "" PARENT_SCOPE)
  if (NOT TARGET ${target})
    return()
  endif()

  get_target_property(libs ${target} INTERFACE_LINK_LIBRARIES)
  set(retval "${libs}")
  foreach(lib ${libs})
    get_recursive_libs(recursive_val ${lib} ${prop} ${searched})
    list(APPEND retval ${recursive_val})
  endforeach()

  list(REMOVE_DUPLICATES retval)
  set(${retvar} "${retval}" PARENT_SCOPE)
endfunction()

function(get_recursive_props retvar target prop)
  set(${retvar} "" PARENT_SCOPE)
  if (NOT TARGET ${target})
    return()
  endif()
  get_recursive_libs(libs target)

  set(retval "")
  foreach(lib ${libs})
    get_property(propval TARGET ${lib} PROPERTY ${prop} SET)
    if (propval)
      get_target_property(propval ${target} ${prop})
      list(APPEND retval ${propval})
    endif()
  endforeach()

  set(${retvar} ${retval} PARENT_SCOPE)
endfunction()

set(runtime_includes "")
get_recursive_props(val heart_gpu_aware INTERFACE_SYSTEM_INCLUDE_DIRECTORIES)
list(APPEND runtime_includes ${val})
get_recursive_props(val heart_gpu_aware INTERFACE_INCLUDE_DIRECTORIES)
list(APPEND runtime_includes ${val})

set(runtime_defs "")
get_recursive_props(val heart_gpu_aware INTERFACE_COMPILE_DEFINITIONS)
list(APPEND runtime_defs ${val})

set(runtime_flags "")
get_recursive_props(val heart_gpu_aware INTERFACE_COMPILE_OPTIONS)
list(APPEND runtime_flags ${val})

set(heart_compile_flags "")
foreach(dir ${runtime_includes})
  set(heart_compile_flags "${heart_compile_flags} -I${dir}")
endforeach()

foreach(flag ${runtime_flags})
  set(heart_compile_flags "${heart_compile_flags} ${flag}")
endforeach()

foreach(def ${runtime_defs})
  set(heart_compile_flags "${heart_compile_flags} -D${def}")
endforeach()

string(TOUPPER "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}" CMAKE_C_FLAGS_BUILDNAME)
set(CMAKE_C_FLAGS_BUILD "${${CMAKE_C_FLAGS_BUILDNAME}}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compileReaction.template
               ${CMAKE_BINARY_DIR}/bin/compileReaction
               @ONLY
               )
install(PROGRAMS ${CMAKE_BINARY_DIR}/bin/compileReaction
	DESTINATION bin
	)
