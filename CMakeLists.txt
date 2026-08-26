#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(BBRSim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(BBRSim BBRSim.cc ${sources} ${headers})
target_link_libraries(BBRSim ${Geant4_LIBRARIES} )

add_executable(BBRLightPipe BBRLightPipe.cc ${sources} ${headers})
target_link_libraries(BBRLightPipe ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy macros to the build directory (the directory in which we build BBRSim)
# so the executable can be run directly from there: it reads these macros from
# the current working directory.
#
set(BBRSim_SCRIPTS
    vis.mac
    test.mac
    test_10K.mac
    test_50M.mac
    test_of_cu.mac
    test_hp_cu.mac
    reflectance.mac
    planck.mac
    verify_config_mt.mac
    lightpipe.mac
    )

foreach(_script ${BBRSim_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS BBRSim DESTINATION bin)

