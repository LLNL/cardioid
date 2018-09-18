
project(cardioid_elec)
cmake_minimum_required(VERSION 3.1)

macro(use_c99)
  if (CMAKE_VERSION VERSION_LESS "3.1.3")
    if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_C_FLAGS "-std=gnu99 ${CMAKE_C_FLAGS}")
    endif ()
  else ()
    set (CMAKE_C_STANDARD 99)
  endif ()
endmacro(use_c99)

use_c99()

# Set a default build type if none was specified
set(default_build_type "Release")
#if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
#  set(default_build_type "Debug")
#endif()
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#Fixing bad default behavior for plugins on OSX
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  SET(CMAKE_SHARED_MODULE_CREATE_C_FLAGS "${CMAKE_SHARED_MODULE_CREATE_C_FLAGS} -undefined dynamic_lookup")
  SET(CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS "${CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()


add_definitions(-DDiff_Weight_Type_Single)
add_definitions(-DWITH_PIO -DWITH_MPI)



function(add_to_nvcc inString)
  set(_list "${inString}")
  separate_arguments(_list)
  foreach(flag IN LISTS _list)
    if (NOT "${flag}" STREQUAL "")
      set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler ${flag}" PARENT_SCOPE)
    endif()
  endforeach()
endfunction()

find_package(OpenMP)

if (ENABLE_CUDA)
   add_definitions(-DUSE_CUDA)
   add_definitions(-DUSE_NVTX)
   find_package(CUDA QUIET REQUIRED)
else()
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()   

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ddcMDutil/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ddcMDutil/include)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)

#find_package(CUDA QUIET REQUIRED)
#find_library(CUDA_NVTX_LIBRARY
#  NAMES nvToolsExt nvTools nvtoolsext nvtools nvtx NVTX
#  PATHS ${CUDA_TOOLKIT_ROOT_DIR}
#  PATH_SUFFIXES "lib64" "common/lib64" "common/lib" "lib"
#  DOC "Location of the CUDA Toolkit Extension (NVTX) library"
#  NO_DEFAULT_PATH
#  )
#mark_as_advanced(CUDA_NVTX_LIBRARY)
#set(CUDA_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_NVTX_LIBRARY})

find_program(PERL perl)

if (BLAS_LIB AND LAPACK_LIB)
   add_definitions(-DHAVE_LAPACK)
   list(APPEND extra_libs ${LAPACK_LIB} ${BLAS_LIB})
endif()

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
        COMMAND ${PERL} ARGS registerBuiltinReactions.pl ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc ${reaction_src}
        MAIN_DEPENDENCY registerBuiltinReactions.pl
        DEPENDS ${reaction_src}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
 )

set (heart_transport_src
   SpaceAllocator.cc ContextManager.cc 
   lazy_array.hh
)
add_library(heart_transport ${heart_transport_src})
if (ENABLE_CUDA)
   target_link_libraries(heart_transport ${cuda_libs})
endif()


if (ENABLE_CUDA)
   set (heart_explicit_cuda_src
      CUDADiffusion.cc
      CUDADiffusion.hh
      CUDADiffusion_ker.cu
      ECGSensor.cc
      ECGSensor.hh
      ECGSensor.cu
   )
   cuda_add_library(heart_explicit_cuda ${heart_explicit_cuda_src})
   target_link_libraries(heart_explicit_cuda heart_transport)
endif()

set (reaction_src
	#TT04_CellML_Reaction.cc
	#TT06_CellML_Reaction.cc 
	#TT06Dev_Reaction.cc
	#SimpleOHaraRudy.cc
	#SimpleGrandi.cc
        #Passive.cc
        #ReactionFHN.cc
        #sundnes_et_al_2016_FHN.cc
	TT06Tau.cc
	TT06Func.cc
        BetterTT06.cc
	pade.cc
        #TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc
        #TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc
	#TT06_RRG.cc
)   

set (heart_gpu_aware_src
   ${reaction_src}
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
   Reaction.cc
   Reaction.hh
   ReactionManager.cc
   ReactionManager.hh
   Simulate.cc
   Simulate.hh
   StateVariableSensor.cc
   Stimulus.hh
   TestStimulus.cc
   TestStimulus.hh
   checkpointIO.cc
   getRemoteCells.cc
   simulationLoop.cc
)
add_library(heart_gpu_aware ${heart_gpu_aware_src})
set(heart_gpu_aware_libs heart_transport ddcMDutil ${extra_libs})
if (ENABLE_CUDA)
  list(APPEND heart_gpu_aware_libs heart_explicit_cuda)
endif()
target_link_libraries(heart_gpu_aware ${heart_gpu_aware_libs})


set (heart_cpu_only_src 
   spi_impl.c
   clooper.c TT06GatesSimd2.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c svd.c
	ThreadServer.cc
	ThreadUtils.cc
	stringUtils.cc
	readCellList.cc
	PioHeaderData.cc
	Vector.cc SymmetricTensor.cc
	PerformanceTimers.cc unionOfStrings.cc getUserInfo.cc
        GDLoadBalancer.cc BlockLoadBalancer.cc workBoundBalancer.cc
	pioBalancer.cc
	CommTable.cc
	BucketOfBits.cc
	stateLoader.cc
	readPioFile.cc
	AnatomyReader.cc object_cc.cc
	Koradi.cc GridRouter.cc Grid3DStencil.cc writeCells.cc
	checkpointIO.cc
	DomainInfo.cc
	BoundingBox.cc
	initializeSimulate.cc
	initializeAnatomy.cc setConductivity.cc assignCellsToTasks.cc
	diffusionFactory.cc reactionFactory.cc
	stimulusFactory.cc sensorFactory.cc drugFactory.cc
   Drug.cc DrugChannel.cc
	FibreConductivity.cc JHUConductivity.cc
	DiffusionUtils.cc
	#FGRDiffusion.cc
	#FGRDiffusionOMP.cc
	#FGRDiffusionThreads.cc
	#FGRDiffusionOverlap.cc
	#FGRDiffusionStrip.cc
	#OpenmpGpuRedblackDiffusion.cc
	#OpenmpGpuFlatDiffusion.cc
	FGRUtils.cc
	PeriodicPulse.cc RandomPulse.cc
	VoronoiCoarsening.cc
	CaAverageSensor.cc
        ProcBox.cc
   Interpolation.cc
   Version.cc
   ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
)
if (NOT use_precompiled)
  list(APPEND heart_cpu_only_src TT06GatesSimd.c)
endif()
add_library(heart_cpu_only ${heart_cpu_only_src})
target_link_libraries(heart_cpu_only ddcMDutil ${extra_libs})


add_executable(cardioid cardioid.cc)
target_link_libraries(cardioid heart_gpu_aware heart_cpu_only heart_gpu_aware)

add_executable(singleCell singleCell.cc singleCellOptions.c)
target_link_libraries(singleCell heart_gpu_aware heart_cpu_only heart_gpu_aware)

if (LAPACK_LIB)
   add_executable(modifyAnatomyFile modifyAnatomyFile.cc)
   target_link_libraries(modifyAnatomyFile heart_gpu_aware heart_cpu_only heart_gpu_aware)
endif ()



### Getting compile flags for making plugins.

set(heart_compile_flags "")

get_property(p_dirs TARGET heart_gpu_aware PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${p_dirs})
  set(heart_compile_flags "${heart_compile_flags} -I${dir}")
endforeach()

get_property(p_flags TARGET heart_gpu_aware PROPERTY COMPILE_FLAGS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} ${flag}")
endforeach()

get_property(p_flags TARGET heart_gpu_aware PROPERTY COMPILE_DEFINITIONS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags TARGET heart_gpu_aware PROPERTY COMPILE_DEFINITIONS_RELEASE)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags DIRECTORY ${CMAKE_CURRENT_SOURCE_DIRECTORY} PROPERTY COMPILE_DEFINITIONS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags DIRECTORY ${CMAKE_CURRENT_SOURCE_DIRECTORY} PROPERTY COMPILE_FLAGS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} ${flag}")
endforeach()

string(TOUPPER "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}" CMAKE_C_FLAGS_BUILDNAME)
set(CMAKE_C_FLAGS_BUILD "${${CMAKE_C_FLAGS_BUILDNAME}}")


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compileReaction.template
               ${CMAKE_BINARY_DIR}/compileReaction
               @ONLY
               )

