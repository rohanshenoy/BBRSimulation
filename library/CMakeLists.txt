#----------------------------------------------------------------------------
# BBRsim library: target BBRsim -> libBBRsim, exported as BBRsim::BBRsim.
# Sources and headers are listed explicitly (no globs); the regression
# runner's drift guard checks that every file under include/ and src/ is here.

set(BBRsim_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBEvt.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBRConfigManager.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBRConfigMessenger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBRCrackLibrary.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBRHFSSData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBSimOpBoundaryProcess.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/BBSimPhysics.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GeometricSurface.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GetBBSpecCDF.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ThermalSurface.cc
  )
set(BBRsim_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBEvt.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBRConfigManager.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBRConfigMessenger.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBRCrackLibrary.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBRHFSSData.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBRMaterials.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBSimOpBoundaryProcess.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/BBSimPhysics.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/CADMesh.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/GeometricSurface.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/GetBBSpecCDF.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/ThermalSurface.hh
  )

# RPATH: installed binaries find libBBRsim and Geant4 without DYLD_LIBRARY_PATH.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" is_SYSTEM_DIR)
if("${is_SYSTEM_DIR}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

add_library(BBRsim SHARED ${BBRsim_SOURCES} ${BBRsim_HEADERS})
if(APPLE)
  set_property(TARGET BBRsim APPEND PROPERTY INSTALL_RPATH "@loader_path")
endif()
target_include_directories(BBRsim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/BBRsim>)
# Geant4 11.4 exports imported targets, so include dirs and compile
# definitions reach consumers through this PUBLIC link.
target_link_libraries(BBRsim PUBLIC ${Geant4_LIBRARIES})
target_compile_features(BBRsim PUBLIC cxx_std_${Geant4_CXX_STANDARD})
# Default data root, baked at configure time; $BBRSIMDATA and /bbr/dataDir
# override it at run time. Changing CMAKE_INSTALL_PREFIX requires a reconfigure.
target_compile_definitions(BBRsim PRIVATE
  "BBRSIM_DEFAULT_DATADIR=\"${CMAKE_INSTALL_PREFIX}/share/BBRsim/data\"")

install(TARGETS BBRsim EXPORT BBRsimTargets LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES ${BBRsim_HEADERS} DESTINATION include/BBRsim)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data DESTINATION share/BBRsim PATTERN .DS_Store EXCLUDE)
install(EXPORT BBRsimTargets NAMESPACE BBRsim:: DESTINATION lib/cmake/BBRsim)
install(FILES ${PROJECT_SOURCE_DIR}/bbrsim_env.sh ${PROJECT_SOURCE_DIR}/bbrsim_env.csh
        DESTINATION share/BBRsim)
