cmake_minimum_required(VERSION 3.20)

project(METIS_INTERNAL)

if (NOT METIS_INTERNAL_LIBRARY_NAME)
  set(METIS_INTERNAL_LIBRARY_NAME metis_internal)
endif()

if (NOT PARMETIS_INTERNAL_LIBRARY_NAME)
  set(PARMETIS_INTERNAL_LIBRARY_NAME parmetis_internal)
endif()

#set default build type to relwithdebinfo, also sets flags O2 on linux
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  mark_as_advanced(CMAKE_BUILD_TYPE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
      "MinSizeRel" "RelWithDebInfo")
endif()

if(WIN32)
  add_definitions("/Zi")
  add_definitions("-DWINDOWS")
endif()

add_subdirectory(METISLib)
