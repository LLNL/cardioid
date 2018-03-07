
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

find_program(PERL perl)

if (BLAS_LIB AND LAPACK_LIB)
   add_definitions(-DHAVE_LAPACK)
endif()



set (reaction_src
        ReactionFHN.cc
	TT04_CellML_Reaction.cc
	TT06_CellML_Reaction.cc 
	TT06Dev_Reaction.cc
	TT06_RRG_Reaction.cc
	OHaraRudy_Reaction.cc
	SimpleOHaraRudy.cc
	Grandi_Reaction.cc
	SimpleGrandi.cc
        Passive.cc
	TT06Tau.cc
	TT06Func.cc
	pade.cc
   TT04_CellML_Endo.cc TT04_CellML_Mid.cc TT04_CellML_Epi.cc
   TT06_CellML_Endo.cc TT06_CellML_Mid.cc TT06_CellML_Epi.cc
	TT06_RRG.cc
   OHaraRudyScan.cc
   GrandiScan.cc 
   spi_impl.c
   clooper.c TT06GatesSimd2.c TT06Gates.c TT06NonGates.c TT06NonGatesSimd.c workBound.c svd.c
   OHaraRudy.c OHaraRudy_Voltage.c OHaraRudy_Concent.c OHaraRudy_INaCai.c OHaraRudy_INaCass.c
              OHaraRudy_INaK.c OHaraRudy_INab.c OHaraRudy_ICab.c OHaraRudy_IKb.c
              OHaraRudy_IpCa.c OHaraRudy_INaFast.c OHaraRudy_INaL.c OHaraRudy_Ito.c
              OHaraRudy_IKr.c OHaraRudy_IKs.c OHaraRudy_IK1.c OHaraRudy_ICa.c null_INull.c
              OHaraRudy_Fluxes.c OHaraRudy_CaMKtrap.c reversalPotentials.c
        OHaraRudyMod_INaFast.c RTYSC14A_IKr.c  MYBGBKC_INa.c
        Grandi.c         Grandi_IClCa.c   Grandi_IKs.c     Grandi_INaL.c    Grandi_Revpots.c
              Grandi_Concent.c Grandi_IClb.c    Grandi_IKur.c    Grandi_INab.c    Grandi_Voltage.c
              Grandi_Fluxes.c  Grandi_IK1.c     Grandi_INCX.c    Grandi_INull.c
              Grandi_ICa.c     Grandi_IKp.c     Grandi_INa.c     Grandi_IpCa.c
              Grandi_ICab.c    Grandi_IKr.c     Grandi_INaK.c    Grandi_Ito.c
)   

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
        COMMAND ${PERL} ARGS registerBuiltinReactions.pl ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc ${reaction_src}
        MAIN_DEPENDENCY registerBuiltinReactions.pl
        DEPENDS ${reaction_src}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
 )



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
	FGRUtils.cc
	Reaction.cc
	ReactionManager.cc
	TestStimulus.cc PointStimulus.cc BoxStimulus.cc
	PeriodicPulse.cc RandomPulse.cc
	PointListSensor.cc ActivationTimeSensor.cc
	ActivationAndRecoverySensor.cc
	VoronoiCoarsening.cc DataVoronoiCoarsening.cc GradientVoronoiCoarsening.cc
	CaAverageSensor.cc
	ECGSensor.cc
	MinMaxSensor.cc ProcBox.cc MaxDVSensor.cc DVThreshSensor.cc
   StateVariableSensor.cc Interpolation.cc
   Version.cc
   ${CMAKE_CURRENT_BINARY_DIR}/registerBuiltinReactions.cc
)
if (NOT use_precompiled)
  list(APPEND reaction_src TT06GatesSimd.c)
endif()
add_library(heart
            ${heart_src} ${reaction_src})
target_link_libraries(heart ddcMDutil ${extra_libs})
if (BLAS_LIB AND LAPACK_LIB)
  target_link_libraries(heart ${LAPACK_LIB} ${BLAS_LIB})
endif()

add_executable(cardioid cardioid.cc)
target_link_libraries(cardioid heart)

add_executable(singleCell singleCell.cc singleCellOptions.c)
target_link_libraries(singleCell heart)

if (LAPACK_LIB)
   add_executable(modifyAnatomyFile modifyAnatomyFile.cc)
   target_link_libraries(modifyAnatomyFile heart)
endif ()



### Getting compile flags for making plugins.

set(heart_compile_flags "")

get_property(p_dirs TARGET heart PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${p_dirs})
  set(heart_compile_flags "${heart_compile_flags} -I${dir}")
endforeach()

get_property(p_flags TARGET heart PROPERTY COMPILE_FLAGS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} ${flag}")
endforeach()

get_property(p_flags TARGET heart PROPERTY COMPILE_DEFINITIONS)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags TARGET heart PROPERTY COMPILE_DEFINITIONS_RELEASE)
foreach(flag ${p_flags})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags DIRECTORY ${CMAKE_CURRENT_SOURCE_DIRECTORY} PROPERTY COMPILE_DEFINITIONS)
foreach(flag ${p_flags})
  MESSAGE(ERROR ${flag})
  set(heart_compile_flags "${heart_compile_flags} -D${flag}")
endforeach()

get_property(p_flags DIRECTORY ${CMAKE_CURRENT_SOURCE_DIRECTORY} PROPERTY COMPILE_FLAGS)
foreach(flag ${p_flags})
  MESSAGE(ERROR ${flag})
  set(heart_compile_flags "${heart_compile_flags} ${flag}")
endforeach()

string(TOUPPER "CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}" CMAKE_C_FLAGS_BUILDNAME)
set(CMAKE_C_FLAGS_BUILD "${${CMAKE_C_FLAGS_BUILDNAME}}")


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/compileReaction.template
               ${CMAKE_BINARY_DIR}/compileReaction
               @ONLY
               )

