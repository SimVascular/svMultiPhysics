if(NOT CPACK_PACKAGE_FILES)
  message(FATAL_ERROR "CPack did not provide a macOS package to post-process")
endif()

if(CPACK_PRODUCTBUILD_IDENTITY_NAME)
  message(FATAL_ERROR "Sign the package after this script; repacking it discards the productbuild signature")
endif()

foreach(_package IN LISTS CPACK_PACKAGE_FILES)
  if(NOT _package MATCHES "[.]pkg$")
    continue()
  endif()

  set(_expanded "${_package}.expanded")
  set(_repacked "${_package}.arm64")
  file(REMOVE_RECURSE "${_expanded}")
  file(REMOVE "${_repacked}")

  execute_process(
    COMMAND /usr/sbin/pkgutil --expand "${_package}" "${_expanded}"
    RESULT_VARIABLE _result
    ERROR_VARIABLE _error)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Failed to expand ${_package}: ${_error}")
  endif()

  set(_distribution "${_expanded}/Distribution")
  if(NOT EXISTS "${_distribution}")
    message(FATAL_ERROR "Expanded package has no Distribution file: ${_package}")
  endif()

  file(READ "${_distribution}" _contents)
  if(_contents MATCHES "hostArchitectures=")
    message(FATAL_ERROR "Package Distribution already declares host architectures")
  endif()
  string(REGEX MATCHALL "<options[^>]*/>" _options "${_contents}")
  list(LENGTH _options _options_count)
  if(NOT _options_count EQUAL 1)
    message(FATAL_ERROR
      "Expected one options element in ${_distribution}, found ${_options_count}")
  endif()

  string(REPLACE "<options "
    "<options hostArchitectures=\"arm64\" " _contents "${_contents}")
  file(WRITE "${_distribution}" "${_contents}")

  execute_process(
    COMMAND /usr/sbin/pkgutil --flatten "${_expanded}" "${_repacked}"
    RESULT_VARIABLE _result
    ERROR_VARIABLE _error)
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR "Failed to rebuild ${_package}: ${_error}")
  endif()

  file(REMOVE "${_package}")
  file(RENAME "${_repacked}" "${_package}")
  file(REMOVE_RECURSE "${_expanded}")
endforeach()
