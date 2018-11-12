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

#blt_register_library(
#  NAME kdtree
#  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/libkdtree
#  )

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/elec)

#macro(makefile_get_var makefile name var)
#  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/echo.mak
#"
#echo-%:
#	-@echo -n \${\$*}
#")
#  execute_process (COMMAND ${CMAKE_MAKE_PROGRAM} -f ${makefile} -f ${CMAKE_CURRENT_BINARY_DIR}/echo.mak --no-print-directory echo-${name}
#      OUTPUT_VARIABLE ${var}
#      RESULT_VARIABLE makefile_return)
#  file (REMOVE ${CMAKE_CURRENT_BINARY_DIR}/echo.mak)
#endmacro()

#makefile_get_var(${CMAKE_CURRENT_SOURCE_DIR}/deps/share/mfem/config.mk MFEM_LIB_FILE MFEM_LIBRARY)
#makefile_get_var(${CMAKE_CURRENT_SOURCE_DIR}/deps/share/mfem/config.mk MFEM_EXT_LIBS  MFEM_LINK_FLAGS)

#blt_register_library(
#  NAME mfem
#  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/deps/include
#  LIBRARIES ${MFEM_LIBRARY}
#  LINK_FLAGS ${MFEM_LINK_FLAGS}
#)


#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ecg)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/fibers)
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mech)
