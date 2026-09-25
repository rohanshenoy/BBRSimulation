#----------------------------------------------------------------------------
# BBRsim — top-level project: builds and installs the BBRsim library.
# Layout mirrors G4CMP (library/, examples/, validation/, tools/); the examples
# and the validation smoke test are standalone projects built against the
# installed library (source bbrsim_env.sh or pass CMAKE_PREFIX_PATH).
cmake_minimum_required(VERSION 3.16...3.27)
# C is enabled too (as before the reorg) so the presets' CMAKE_C_COMPILER is used.
project(BBRsim VERSION 0.1.0 LANGUAGES C CXX)

option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
option(INSTALL_EXAMPLES "Copy examples/ to the installation area" OFF)
option(INSTALL_VALIDATION "Copy validation/ to the installation area" ON)
set(BBRSIM_VERSION "${PROJECT_VERSION}" CACHE STRING
    "Version string used when the source is not a git checkout")

if(WITH_GEANT4_UIVIS)
  set(BBRSIM_GEANT4_COMPONENTS ui_all vis_all)
  find_package(Geant4 REQUIRED ${BBRSIM_GEANT4_COMPONENTS})
else()
  set(BBRSIM_GEANT4_COMPONENTS "")
  find_package(Geant4 REQUIRED)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_subdirectory(library)

# Version stamp (G4CMP CMakeLists.txt:83-95). git -C also works in a worktree,
# where .git is a file; describe needs an annotated tag to report one.
find_program(SH sh HINTS /bin)
add_custom_target(version ALL
  COMMENT "Building version identification file"
  VERBATIM
  COMMAND ${CMAKE_COMMAND} -E env ${SH} -c
    "(git -C ${PROJECT_SOURCE_DIR} describe --always --dirty 2>/dev/null || echo ${BBRSIM_VERSION}) > ${PROJECT_BINARY_DIR}/.bbrsim-version")
install(FILES ${PROJECT_BINARY_DIR}/.bbrsim-version DESTINATION share/BBRsim)

# Package config for find_package(BBRsim).
include(CMakePackageConfigHelpers)
set(BBRSIM_DATA_INSTALL_DIR share/BBRsim/data)
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/BBRsimConfig.cmake.in
  ${PROJECT_BINARY_DIR}/BBRsimConfig.cmake
  INSTALL_DESTINATION lib/cmake/BBRsim
  PATH_VARS BBRSIM_DATA_INSTALL_DIR)
write_basic_package_version_file(${PROJECT_BINARY_DIR}/BBRsimConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/BBRsimConfig.cmake
              ${PROJECT_BINARY_DIR}/BBRsimConfigVersion.cmake
        DESTINATION lib/cmake/BBRsim)
if(INSTALL_EXAMPLES)
  install(DIRECTORY examples DESTINATION . PATTERN build EXCLUDE PATTERN .DS_Store EXCLUDE)
endif()
if(INSTALL_VALIDATION)
  # USE_SOURCE_PERMISSIONS keeps the runner and drift guards executable; the last
  # four patterns keep run artifacts (plots, ROOT files, build/output dirs) out.
  install(DIRECTORY validation DESTINATION . USE_SOURCE_PERMISSIONS
          PATTERN __pycache__ EXCLUDE PATTERN .DS_Store EXCLUDE
          PATTERN build EXCLUDE PATTERN output EXCLUDE
          PATTERN "*.root" EXCLUDE PATTERN "*.png" EXCLUDE)
endif()
