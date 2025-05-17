cmake_minimum_required(VERSION 3.20)

project(GKLIB_INTERNAL)

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

include_directories(./)
include_directories(${MPI_C_INCLUDE_PATH})

set(CSRCS b64.c blas.c cache.c csr.c error.c evaluate.c
  fkvkselect.c fs.c getopt.c gk_util.c gkregex.c graph.c
  htable.c io.c itemsets.c mcore.c memory.c pqueue.c
  random.c rw.c seq.c sort.c string.c timers.c tokenizer.c
)

if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${GKLIB_INTERNAL_LIBRARY_NAME} SHARED ${CSRCS})
else()
  add_library(${GKLIB_INTERNAL_LIBRARY_NAME} STATIC ${CSRCS})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${GKLIB_INTERNAL_LIBRARY_NAME}
  RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
  LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
  ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries)
endif()

