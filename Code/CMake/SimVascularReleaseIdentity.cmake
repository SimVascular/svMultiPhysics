# Resolve the identity from the source tree, never from the workflow checkout.
get_filename_component(SV_REPOSITORY_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
execute_process(
  COMMAND git -C "${SV_REPOSITORY_ROOT}" rev-parse --verify HEAD
  OUTPUT_VARIABLE SV_SOURCE_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
  RESULT_VARIABLE _sv_git_result)
if(NOT _sv_git_result EQUAL 0)
  set(SV_SOURCE_COMMIT "unknown")
endif()
execute_process(
  COMMAND git -C "${SV_REPOSITORY_ROOT}" status --porcelain --untracked-files=normal
  OUTPUT_VARIABLE _sv_status
  ERROR_QUIET
  RESULT_VARIABLE _sv_status_result)

set(SV_PACKAGE_VERSION "" CACHE STRING "Four-part version from release.json")
file(READ "${SV_REPOSITORY_ROOT}/release.json" _sv_release_json)
foreach(_part major minor patch build)
  string(JSON _value ERROR_VARIABLE _error GET "${_sv_release_json}" "${_part}")
  if(_error OR NOT _value MATCHES "^[0-9]+$")
    message(FATAL_ERROR "release.json must contain a numeric ${_part}")
  endif()
  list(APPEND _sv_version_parts "${_value}")
endforeach()
list(JOIN _sv_version_parts "." _sv_manifest_version)
if(SV_PACKAGE_VERSION AND NOT SV_PACKAGE_VERSION STREQUAL _sv_manifest_version)
  message(FATAL_ERROR "SV_PACKAGE_VERSION disagrees with release.json: ${SV_PACKAGE_VERSION} != ${_sv_manifest_version}")
endif()
set(SV_PACKAGE_VERSION "${_sv_manifest_version}" CACHE STRING "Four-part version from release.json" FORCE)

if(SV_SOURCE_COMMIT STREQUAL "unknown")
  set(SV_RELEASE_TAG "unknown")
elseif(_sv_status_result EQUAL 0 AND _sv_status STREQUAL "")
  execute_process(
    COMMAND git -C "${SV_REPOSITORY_ROOT}" describe --tags --exact-match HEAD
    OUTPUT_VARIABLE SV_RELEASE_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _sv_tag_result)
  if(NOT _sv_tag_result EQUAL 0)
    set(SV_RELEASE_TAG "untagged")
  endif()
else()
  set(SV_RELEASE_TAG "uncommitted")
endif()

if(SV_RELEASE_PACKAGE)
  if(NOT SV_ENABLE_DISTRIBUTION)
    message(FATAL_ERROR "SV_RELEASE_PACKAGE requires SV_ENABLE_DISTRIBUTION")
  endif()
  if(NOT _sv_status_result EQUAL 0 OR NOT _sv_status STREQUAL "")
    message(FATAL_ERROR "Release builds require a clean source tree")
  endif()
  execute_process(
    COMMAND git -C "${SV_REPOSITORY_ROOT}" rev-parse --verify "refs/tags/v${SV_PACKAGE_VERSION}^{commit}"
    OUTPUT_VARIABLE _sv_tag_commit
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
    RESULT_VARIABLE _sv_tag_result)
  if(NOT _sv_tag_result EQUAL 0 OR NOT _sv_tag_commit STREQUAL SV_SOURCE_COMMIT)
    message(FATAL_ERROR "Release tag v${SV_PACKAGE_VERSION} must exist and point to source commit ${SV_SOURCE_COMMIT}")
  endif()
  set(SV_RELEASE_TAG "v${SV_PACKAGE_VERSION}")
endif()

message(STATUS "svMultiPhysics source: ${SV_RELEASE_TAG} (${SV_SOURCE_COMMIT})")
