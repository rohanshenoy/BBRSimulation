#----------------------------------------------------------------------------
# bbrsimLightPipe — 4 K -> mixing-chamber light pipe (parametric or CAD/STL mode,
# /bbr/lightpipe/*). Standalone project built against an installed BBRsim:
# source <prefix>/share/BBRsim/bbrsim_env.sh or pass -DCMAKE_PREFIX_PATH=<prefix>.
cmake_minimum_required(VERSION 3.16...3.27)
project(bbrsimLightPipe LANGUAGES C CXX)

# BBRSIMINSTALL is <prefix>/share/BBRsim after sourcing an installed env script.
# Prepended, so the sourced prefix wins over the preset's ../../install.
if(DEFINED ENV{BBRSIMINSTALL} AND EXISTS "$ENV{BBRSIMINSTALL}/../../lib/cmake/BBRsim")
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{BBRSIMINSTALL}/../..")
endif()
find_package(BBRsim REQUIRED)   # brings Geant4 through find_dependency

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# RPATH: installed binaries find their libraries without DYLD_LIBRARY_PATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" is_SYSTEM_DIR)
if("${is_SYSTEM_DIR}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

set(lightpipe_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipeActionInitialization.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipeDetectorConstruction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipeMessenger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipePrimaryGeneratorAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipeRunAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/LightPipeSteppingAction.cc
  )
set(lightpipe_MACROS vis.mac lightpipe.mac)

add_library(lightpipeLib SHARED ${lightpipe_SOURCES})
set_target_properties(lightpipeLib PROPERTIES OUTPUT_NAME bbrsimLightPipe)
target_include_directories(lightpipeLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(lightpipeLib PUBLIC BBRsim::BBRsim)

add_executable(bbrsimLightPipe bbrsimLightPipe.cc)
target_link_libraries(bbrsimLightPipe PRIVATE lightpipeLib)

# Macros beside the binary, so it runs from its build directory (g4beginner pattern).
foreach(_m IN LISTS lightpipe_MACROS)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/G4Macros/${_m} ${PROJECT_BINARY_DIR}/${_m} COPYONLY)
endforeach()

install(TARGETS lightpipeLib DESTINATION lib)
install(TARGETS bbrsimLightPipe DESTINATION bin)
list(TRANSFORM lightpipe_MACROS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/G4Macros/ OUTPUT_VARIABLE _macro_paths)
install(FILES ${_macro_paths} DESTINATION macros/bbrsimLightPipe)
