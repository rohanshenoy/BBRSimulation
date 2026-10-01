execute_process(COMMAND git -C "${REPO}" describe --always --dirty
  OUTPUT_VARIABLE version OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT version)
  set(version "${FALLBACK}")
endif()
# Reproducible source-content identity, including uncommitted edits. The
# library and each application call this separately so rebuilding one example
# never changes the library's source fingerprint.
if(NOT DEFINED SOURCE_DIR OR SOURCE_DIR STREQUAL "")
  set(SOURCE_DIR "${REPO}/library")
endif()
if(NOT DEFINED MACRO_PREFIX OR MACRO_PREFIX STREQUAL "")
  set(MACRO_PREFIX BBRSIM)
endif()
file(GLOB_RECURSE sources RELATIVE "${SOURCE_DIR}"
  "${SOURCE_DIR}/include/*.hh" "${SOURCE_DIR}/include/*.cc"
  "${SOURCE_DIR}/src/*.hh" "${SOURCE_DIR}/src/*.cc")
file(GLOB root_sources RELATIVE "${SOURCE_DIR}"
  "${SOURCE_DIR}/*.hh" "${SOURCE_DIR}/*.cc")
list(APPEND sources ${root_sources})
if(EXISTS "${SOURCE_DIR}/CMakeLists.txt")
  list(APPEND sources CMakeLists.txt)
endif()
list(REMOVE_DUPLICATES sources)
list(SORT sources)
set(manifest "")
foreach(source IN LISTS sources)
  # Generated build directories are not source provenance.
  if(NOT source MATCHES "(^|/)build[^/]*/")
    file(SHA256 "${SOURCE_DIR}/${source}" checksum)
    string(APPEND manifest "${source}:${checksum}\n")
  endif()
endforeach()
string(SHA256 fingerprint "${manifest}")
set(content "#pragma once\n#define ${MACRO_PREFIX}_BUILD_VERSION R\"bbr(${version})bbr\"\n#define ${MACRO_PREFIX}_SOURCE_FINGERPRINT \"sha256:${fingerprint}\"\n")
set(old "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" old)
endif()
if(NOT content STREQUAL old)
  file(WRITE "${OUTPUT}" "${content}")
endif()
