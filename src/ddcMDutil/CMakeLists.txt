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


project(simUtil)
cmake_minimum_required(VERSION 3.2)

set (CMAKE_C_STANDARD 99)

set(ENABLE_GTEST OFF CACHE STRING "Cardioid" FORCE)
set(BUILD_GTEST OFF CACHE STRING "Cardioid" FORCE)
set(ENABLE_TESTS OFF CACHE STRING "Cardioid" FORCE)
set(INSTALL_GTEST OFF CACHE STRING "Cardioid" FORCE)
set(ENABLE_GTEST OFF CACHE STRING "Cardioid" FORCE)
set(ENABLE_EXAMPLES OFF CACHE STRING "Cardioid" FORCE)

if (NOT BLT_LOADED)
   include(${CMAKE_CURRENT_SOURCE_DIR}/blt/SetupBLT.cmake)
endif()

add_definitions(-DDiff_Weight_Type_Single)
add_definitions(-DWITH_PIO -DWITH_MPI)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
