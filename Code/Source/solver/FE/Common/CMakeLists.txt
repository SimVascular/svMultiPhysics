file(GLOB SVMP_FE_COMMON_HEADERS CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.inl"
)

add_library(svmp_fe_common INTERFACE)

target_sources(svmp_fe_common INTERFACE
  ${SVMP_FE_COMMON_HEADERS}
)

target_compile_features(svmp_fe_common INTERFACE cxx_std_20)

target_include_directories(svmp_fe_common INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_link_libraries(svmp_fe_common INTERFACE
  svmp_core
)
