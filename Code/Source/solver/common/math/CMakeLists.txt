file(GLOB_RECURSE SVMP_MATH_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.inl"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_library(svmp_math OBJECT
  ${SVMP_MATH_FILES}
)

target_compile_features(svmp_math PUBLIC cxx_std_20)

target_include_directories(svmp_math PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${SV_SOURCE_DIR}/ThirdParty/eigen/include
  ${SV_SOURCE_DIR}/ThirdParty/eigen/include/eigen3
)

target_link_libraries(svmp_math PUBLIC
  svmp_core
)
