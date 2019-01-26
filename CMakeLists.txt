# Set a default build type if none was specified
set(default_build_type "Release")
 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build.")
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

SET(CMAKE_C_FLAGS_DEBUG "-g -ggdb -O0 " CACHE STRING "")
SET(CMAKE_C_FLAGS_MINSIZEREL "-g -pg -O3 -DPROFILE" CACHE STRING "")
SET(CMAKE_C_FLAGS_RELEASE "-g -O3" CACHE STRING "") 
SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-g -ggdb -O2" CACHE STRING "")

SET(CMAKE_CXX_FLAGS_DEBUG  ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "")
SET(CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_C_FLAGS_MINSIZEREL} CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" ) 
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CMAKE_C_FLAGS_RELWITHDEBINFO} CACHE STRING "")

SET(CMAKE_CUDA_FLAGS_DEBUG "-G -g -O0" CACHE STRING "")
SET(CMAKE_CUDA_FLAGS_RELEASE "-g -O3" CACHE STRING "")
SET(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -O2" CACHE STRING "")
SET(CMAKE_CUDA_FLAGS_MINSIZEREL "-g -O3" CACHE STRING "")

set (CUDA_ARCH "sm_70" CACHE PATH "")
set (CMAKE_CUDA_FLAGS "-restrict -arch ${CUDA_ARCH} -std=c++11 --expt-extended-lambda" CACHE STRING "" )
set (CMAKE_CUDA_HOST_COMPILER ${MPI_CXX_COMPILER} CACHE STRING "")

project(cardioid)
cmake_minimum_required(VERSION 3.2)

set (CMAKE_C_STANDARD 99)

set(ENABLE_GTEST OFF CACHE BOOL "")
set(BUILD_GTEST OFF CACHE BOOL "")
set(ENABLE_TESTS OFF CACHE BOOL "")
set(INSTALL_GTEST OFF CACHE BOOL "")
set(ENABLE_GTEST OFF CACHE BOOL "")
set(ENABLE_EXAMPLES OFF CACHE BOOL "")

set(ENABLE_ALL_WARNINGS OFF CACHE BOOL "")

if (NOT BLT_LOADED)
   include(${CMAKE_CURRENT_SOURCE_DIR}/blt/SetupBLT.cmake)
endif()


#Fixing bad default behavior for plugins on OSX
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  SET(CMAKE_SHARED_MODULE_CREATE_C_FLAGS "${CMAKE_SHARED_MODULE_CREATE_C_FLAGS} -undefined dynamic_lookup")
  SET(CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS "${CMAKE_SHARED_MODULE_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()


add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/simUtil/src)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/simdops)

blt_register_library(
  NAME kdtree
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libkdtree
  )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/elec)

macro(makefile_get_var makefile name var)
  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/echo.mak
"
echo-%:
	-@echo -n \${\$*}
")
  execute_process (COMMAND ${CMAKE_MAKE_PROGRAM} -f ${makefile} -f ${CMAKE_CURRENT_BINARY_DIR}/echo.mak --no-print-directory echo-${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE makefile_return)
  file (REMOVE ${CMAKE_CURRENT_BINARY_DIR}/echo.mak)
endmacro()

set(MFEM_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps" CACHE PATH "Path to the mfem directory")
if (EXISTS ${MFEM_DIR})
   makefile_get_var(${MFEM_DIR}/share/mfem/config.mk MFEM_LIB_FILE MFEM_LIBRARY)
   makefile_get_var(${MFEM_DIR}/share/mfem/config.mk MFEM_EXT_LIBS  MFEM_LINK_FLAGS)
   
   blt_register_library(
     NAME mfem
     INCLUDES ${MFEM_DIR}/include
     LIBRARIES ${MFEM_LIBRARY}
     LINK_FLAGS ${MFEM_LINK_FLAGS}
   )


   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ecg)
   add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fibers)
   option(BUILD_MECH "Should we build with mechanics?  Needs a non-public version of mfem." OFF)
   if (BUILD_MECH)
      add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mech)
   endif()
endif()
