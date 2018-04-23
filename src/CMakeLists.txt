
project(cardioid_elec)
cmake_minimum_required(VERSION 3.0)

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

find_package(OpenMP)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_definitions(-DDiff_Weight_Type_Single)
add_definitions(-DWITH_PIO -DWITH_MPI)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ddcMDutil/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ddcMDutil/include)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_library(BLAS_LIB blas)
find_library(LAPACK_LIB lapack)

find_package(CUDA QUIET REQUIRED)

if (BLAS_LIB AND LAPACK_LIB)
   add_definitions(-DHAVE_LAPACK)
endif()

set (heart_src 
	ThreadServer.cc
	ThreadUtils.cc
	stringUtils.cc
	readCellList.cc
	PioHeaderData.cc
	Vector.cc SymmetricTensor.cc ConstantModel.cc ConstantReaction.cc
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
	Simulate.cc
	initializeAnatomy.cc setConductivity.cc assignCellsToTasks.cc
	getRemoteCells.cc
	simulationLoop.cc
	diffusionFactory.cc reactionFactory.cc
	stimulusFactory.cc sensorFactory.cc drugFactory.cc
   Drug.cc DrugChannel.cc
	FibreConductivity.cc JHUConductivity.cc
	DiffusionUtils.cc
	FGRDiffusion.cc
	FGRDiffusionOMP.cc
	FGRDiffusionThreads.cc
	FGRDiffusionOverlap.cc
	FGRDiffusionStrip.cc
	OpenmpGpuRedblackDiffusion.cc
	OpenmpGpuFlatDiffusion.cc
	CUDADiffusion.cc
	CUDADiffusion_ker.cu
	FGRUtils.cc
	Reaction.cc
	ReactionManager.cc
        #ReactionFHN.cc
	TT04_CellML_Reaction.cc
	TT06_CellML_Reaction.cc 
	TT06Dev_Reaction.cc
	#TT06_RRG_Reaction.cc
	#OHaraRudy_Reaction.cc
	SimpleOHaraRudy.cc
	SimpleGrandi.cc
        BetterTT06.cc
	SimpleTT06.cc
	SimpleTT06_ker.cu
        Passive.cc
	TT06Tau.cc
	TT06Func.cc
	pade.cc
   TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc
   TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc
	TT06_RRG.cc
	TestStimulus.cc PointStimulus.cc BoxStimulus.cc
	PeriodicPulse.cc RandomPulse.cc
	PointListSensor.cc ActivationTimeSensor.cc
	ActivationAndRecoverySensor.cc
	VoronoiCoarsening.cc DataVoronoiCoarsening.cc GradientVoronoiCoarsening.cc
	CaAverageSensor.cc
	ECGSensor.cc
	ECGSensor.cu
	MinMaxSensor.cc ProcBox.cc MaxDVSensor.cc DVThreshSensor.cc
   StateVariableSensor.cc Interpolation.cc
   Ledger.cc simulationLoopCuda.cu singleCellCuda.cu BoxStimulusCuda.cu HaloExchange.cu
   Version.cc
   spi_impl.c
   clooper.c TT06GatesSimd2.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c svd.c
)
if (NOT use_precompiled)
  list(APPEND heart_src TT06GatesSimd.c)
endif()
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_60 " )
cuda_add_library(heart
            ${heart_src})
target_link_libraries(heart ddcMDutil ${extra_libs})
if (BLAS_LIB AND LAPACK_LIB)
  target_link_libraries(heart ${LAPACK_LIB} ${BLAS_LIB})
endif()

cuda_add_executable(cardioid cardioid.cc)
target_link_libraries(cardioid heart)

cuda_add_executable(singleCell singleCell.cc singleCellOptions.c)
target_link_libraries(singleCell heart)

if (LAPACK_LIB)
   cuda_add_executable(modifyAnatomyFile modifyAnatomyFile.cc)
   target_link_libraries(modifyAnatomyFile heart)
endif ()
