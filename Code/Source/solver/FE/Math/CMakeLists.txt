file(GLOB SVMP_FE_MATH_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.inl"
)

add_library(svmp_fe_math OBJECT
  ${SVMP_FE_MATH_FILES}
)

target_compile_features(svmp_fe_math PUBLIC cxx_std_20)

target_include_directories(svmp_fe_math PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../..
  ${SV_SOURCE_DIR}/ThirdParty/eigen/include
  ${SV_SOURCE_DIR}/ThirdParty/eigen/include/eigen3
)

target_link_libraries(svmp_fe_math PUBLIC
  svmp_fe_common
)
