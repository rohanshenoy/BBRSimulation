#----------------------------------------------------------------------------
# bbrsimTestWorld — BBRsim test world (Cu slab, two HFSS crack volumes, Planck
# emitter or fixed gun). Builds against the parent source-tree target or, as a
# standalone project, an installed BBRsim package.
cmake_minimum_required(VERSION 3.16...3.27)
project(bbrsimTestWorld LANGUAGES C CXX)

if(TARGET BBRsim)
  set(testworld_BBRSIM_TARGET BBRsim)
  set(testworld_BUILD_INFO_SCRIPT ${CMAKE_SOURCE_DIR}/cmake/BBRBuildInfo.cmake)
  set(testworld_BUILD_VERSION ${BBRSIM_VERSION})
else()
  # BBRSIMINSTALL is <prefix>/share/BBRsim after sourcing an installed env
  # script. Prepended, so it wins over the preset's ../../install.
  if(DEFINED ENV{BBRSIMINSTALL} AND EXISTS "$ENV{BBRSIMINSTALL}/../../lib/cmake/BBRsim")
    list(PREPEND CMAKE_PREFIX_PATH "$ENV{BBRSIMINSTALL}/../..")
  endif()
  find_package(BBRsim REQUIRED)   # brings Geant4 through find_dependency
  set(testworld_BBRSIM_TARGET BBRsim::BBRsim)
  set(testworld_BUILD_INFO_SCRIPT ${BBRsim_BUILD_INFO_SCRIPT})
  set(testworld_BUILD_VERSION ${BBRsim_VERSION})
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# RPATH: installed binaries find their libraries without DYLD_LIBRARY_PATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" is_SYSTEM_DIR)
if("${is_SYSTEM_DIR}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

set(testworld_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestWorldActionInitialization.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestWorldDetectorConstruction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestWorldPrimaryGeneratorAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestWorldRunAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/TestWorldSteppingAction.cc
  )
set(testworld_MACROS vis.mac planck.mac planck_10K.mac planck_5M.mac planck_50M.mac
    reflectance.mac reflectance_OF_Cu.mac reflectance_HP_Cu.mac config_mt.mac)

add_library(testworldLib SHARED ${testworld_SOURCES})
set_target_properties(testworldLib PROPERTIES OUTPUT_NAME bbrsimTestWorld)
add_custom_target(testworldBuildInfo
  COMMAND ${CMAKE_COMMAND} "-DREPO=${CMAKE_CURRENT_SOURCE_DIR}"
    "-DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
    "-DOUTPUT=${CMAKE_CURRENT_BINARY_DIR}/BBRApplicationInfo.hh"
    "-DFALLBACK=${testworld_BUILD_VERSION}" "-DMACRO_PREFIX=BBRSIM_APP"
    -P ${testworld_BUILD_INFO_SCRIPT}
  BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/BBRApplicationInfo.hh
  VERBATIM)
add_dependencies(testworldLib testworldBuildInfo)
target_include_directories(testworldLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(testworldLib PUBLIC ${testworld_BBRSIM_TARGET})

add_executable(bbrsimTestWorld bbrsimTestWorld.cc)
target_link_libraries(bbrsimTestWorld PRIVATE testworldLib)

# Macros beside the binary, so it runs from its build directory (g4beginner pattern).
foreach(_m IN LISTS testworld_MACROS)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/G4Macros/${_m} ${PROJECT_BINARY_DIR}/${_m} COPYONLY)
endforeach()

# A top-level source-tree build installs the library only. Standalone example
# projects keep their existing install layout beside an installed BBRsim.
if(NOT TARGET BBRsim)
  install(TARGETS testworldLib DESTINATION lib)
  install(TARGETS bbrsimTestWorld DESTINATION bin)
  list(TRANSFORM testworld_MACROS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/G4Macros/ OUTPUT_VARIABLE _macro_paths)
  install(FILES ${_macro_paths} DESTINATION macros/bbrsimTestWorld)
endif()
