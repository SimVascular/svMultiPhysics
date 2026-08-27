file(GLOB SVMP_FE_BASIS_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.inl"
)

add_library(svmp_fe_basis OBJECT
  ${SVMP_FE_BASIS_FILES}
)

target_compile_features(svmp_fe_basis PUBLIC cxx_std_20)

target_include_directories(svmp_fe_basis PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
)

target_link_libraries(svmp_fe_basis PUBLIC
  svmp_fe_common
  svmp_fe_math
)
