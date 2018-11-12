

set(extra_definitions "" CACHE STRING "")
add_definitions(${extra_definitions})

set(extra_includes "" CACHE STRING "")

if (ENABLE_CUDA)
   include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include)
   add_definitions(-DUSE_CUDA)
   add_definitions(-DUSE_NVTX)
   set(cuda cuda)
   set(cuda_runtime cuda_runtime)

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
   set(cuda "")
   set(cuda_runtime "")
   set(CUDA_NVTX_LIBRARY "")
   set(CUDA_NVRTC_LIBRARY "")
   set(CUDA_CUDA_LIBRARY "")
endif()   

find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)

find_program(PERL perl)

if (BLAS_LIB AND LAPACK_LIB)
   add_definitions(-DHAVE_LAPACK)
   list(APPEND extra_libs ${LAPACK_LIB} ${BLAS_LIB})
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../simdops/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set (heart_transport_src
   SpaceAllocator.cc 
   lazy_array.hh
)

blt_add_library(NAME heart_transport
                SOURCES ${heart_transport_src}
                DEPENDS_ON ${cuda_runtime} openmp)

if (ENABLE_CUDA)
   blt_add_library(NAME heart_actual_cuda
                   SOURCES CUDADiffusion_ker.cu ECGSensor.cu
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

set (reaction_src
	#TT04_CellML_Reaction.cc
	#TT06_CellML_Reaction.cc 
	TT06Dev_Reaction.cc
        Passive.cc
        #ReactionFHN.cc
        #sundnes_et_al_2016_FHN.cc
	TT06Tau.cc
	TT06Func.cc
        BetterTT06.cc
        Grandi.cc
	pade.cc
)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
        COMMAND ${PERL} ARGS registerBuiltinReactions.pl ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc ${reaction_src}
        MAIN_DEPENDENCY registerBuiltinReactions.pl
        DEPENDS ${reaction_src}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
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
   simulationLoop.cc
   Simulate.cc
   Simulate.hh
   StateVariableSensor.cc
   Stimulus.hh
   TestStimulus.cc
   TestStimulus.hh
   checkpointIO.cc
   getRemoteCells.cc
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
	stimulusFactory.cc sensorFactory.cc
	FibreConductivity.cc JHUConductivity.cc
	DiffusionUtils.cc
	#FGRDiffusion.cc
	#FGRDiffusionOMP.cc
	#FGRDiffusionThreads.cc
	#FGRDiffusionOverlap.cc
	FGRDiffusionStrip.cc
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
blt_add_library(NAME heart_gpu_aware
                SOURCES ${heart_gpu_aware_src}
                DEPENDS_ON ${heart_explicit_cuda} heart_transport simUtil simdops ${extra_libs}
                           ${CUDA_NVTX_LIBRARY} ${CUDA_NVRTC_LIBRARY} ${CUDA_CUDA_LIBRARY}
                           ${cuda} ${cuda_runtime} openmp mpi
                INCLUDES ${extra_includes}
                )

set (heart_cpu_only_src 
   spi_impl.c
   clooper.c TT06GatesSimd2.c TT06GatesSimd.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c svd.c
)
blt_add_library(NAME heart_cpu_only
                SOURCES ${heart_cpu_only_src}
                DEPENDS_ON simUtil simdops ${extra_libs} openmp
                INCLUDES ${extra_includes}
                )

blt_add_executable(NAME cardioid
                   SOURCES cardioid.cc
                   DEPENDS_ON heart_gpu_aware heart_cpu_only ${heart_explicit_cuda}
                              heart_gpu_aware heart_cpu_only ${heart_explicit_cuda}
                              heart_gpu_aware heart_cpu_only ${heart_explicit_cuda}
                              ${cuda_runtime} openmp)

blt_add_executable(NAME singleCell
                   SOURCES singleCell.cc singleCellOptions.c
                   DEPENDS_ON heart_gpu_aware heart_cpu_only ${cuda} ${cuda_runtime} openmp)

if (LAPACK_LIB)
   blt_add_executable(NAME modifyAnatomyFile
                      SOURCES modifyAnatomyFile.cc
                      DEPENDS_ON heart_gpu_aware heart_cpu_only ${cuda_runtime} openmp)
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
