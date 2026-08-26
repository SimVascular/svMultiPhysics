file(GLOB SVMP_CORE_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.inl"
)

add_library(svmp_core INTERFACE)

target_sources(svmp_core INTERFACE
  ${SVMP_CORE_HEADERS}
)

target_compile_features(svmp_core INTERFACE cxx_std_20)

target_include_directories(svmp_core INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(svmp_core INTERFACE
  MPI::MPI_CXX
)
