file(GLOB SVMP_FE_BASIS_FILES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.inl"
)

add_library(svmp_fe_basis ${SV_LIBRARY_TYPE}
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

if(SV_INSTALL_LIBS)
  install(TARGETS svmp_fe_basis
    RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT CoreExecutables
    LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT CoreLibraries
    ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT CoreLibraries
  )
endif()
