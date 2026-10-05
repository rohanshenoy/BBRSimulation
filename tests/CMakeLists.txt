#----------------------------------------------------------------------------
# BBRsim C++ tests, built with -DBUILD_BBRSIM_TESTS=ON (the regression runner
# forces it on). Plain executables linked to the in-tree BBRsim target, so the
# header-only library code is compiled with -I and its warnings count; one
# CTest test per case, because the library's singletons cannot be reset within
# a process. File lists are explicit, as elsewhere in the project.

# Every test runs without an inherited DYLD_LIBRARY_PATH / LD_LIBRARY_PATH
# (bbrsim_env.sh sets them, and they would load the INSTALLED libBBRsim instead
# of the build tree's) and without BBRSIMDATA. That takes the test property
# ENVIRONMENT_MODIFICATION, which needs CMake >= 3.22, so the tests need it too.
if(CMAKE_VERSION VERSION_LESS 3.22)
  message(FATAL_ERROR "BUILD_BBRSIM_TESTS needs CMake >= 3.22 (ENVIRONMENT_MODIFICATION keeps the tests hermetic)")
endif()
set(BBRSIM_TEST_HERMETIC "DYLD_LIBRARY_PATH=unset:;LD_LIBRARY_PATH=unset:;BBRSIMDATA=unset:")

# bbrsim_add_test_program(<name> LABEL <label> CASES <c1> <c2> ...
#                         [NEEDS_ENSDF] [ENV <VAR=value> ...]
#                         [ENVMOD <VAR=op:value> ...] [WILL_FAIL <case> ...])
# Builds <name>.cc and registers one test "<name>.<case>" per case, each in
# its own working directory. NEEDS_ENSDF: the case builds a G4OpticalPhoton or
# a run manager, which needs the G4ENSDFSTATE dataset. ENVMOD entries are
# applied after the hermetic ones, so ENVMOD "BBRSIMDATA=set:/x" wins.
# WILL_FAIL lists cases that must exit non-zero (the test inverts).
function(bbrsim_add_test_program name)
  cmake_parse_arguments(T "NEEDS_ENSDF" "LABEL" "CASES;ENV;ENVMOD;WILL_FAIL" ${ARGN})
  add_executable(${name} ${name}.cc BBRTestSupport.hh)
  target_link_libraries(${name} PRIVATE BBRsim)
  target_compile_options(${name} PRIVATE -Wall -Wextra)
  set(env ${T_ENV})
  if(T_NEEDS_ENSDF)
    list(APPEND env "G4ENSDFSTATEDATA=${Geant4_DATASET_G4ENSDFSTATE_PATH}")
  endif()
  foreach(c ${T_CASES})
    set(wd ${CMAKE_CURRENT_BINARY_DIR}/work/${name}.${c})
    file(MAKE_DIRECTORY ${wd})
    add_test(NAME ${name}.${c} COMMAND ${name} ${c})
    set_tests_properties(${name}.${c} PROPERTIES
      LABELS "${T_LABEL}" TIMEOUT 120 WORKING_DIRECTORY ${wd})
    if(env)
      set_tests_properties(${name}.${c} PROPERTIES ENVIRONMENT "${env}")
    endif()
    set(mod ${BBRSIM_TEST_HERMETIC} ${T_ENVMOD})
    set_tests_properties(${name}.${c} PROPERTIES ENVIRONMENT_MODIFICATION "${mod}")
    if(c IN_LIST T_WILL_FAIL)
      set_tests_properties(${name}.${c} PROPERTIES WILL_FAIL TRUE)
    endif()
  endforeach()
endfunction()

# bbrsim_add_abort_test(<test-name> <program> <case> <G4Exception code>)
# The case must die in a fatal G4Exception with the stock handler (abort());
# CTest cannot express that itself, so ExpectAbort.cmake checks it.
function(bbrsim_add_abort_test tname prog case code)
  add_test(NAME ${tname}
    COMMAND ${CMAKE_COMMAND} -DCMD=$<TARGET_FILE:${prog}> -DARGS=${case} -DCODE=${code}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/ExpectAbort.cmake)
  set_tests_properties(${tname} PROPERTIES LABELS "unit" TIMEOUT 60
    ENVIRONMENT_MODIFICATION "${BBRSIM_TEST_HERMETIC}")
endfunction()

# --- test programs ----------------------------------------------------------

bbrsim_add_test_program(testSupport LABEL unit
  CASES handler_throws handler_records_warning scripted_engine tempdir
        check_fails unexpected_exception message_mismatch_fails
  WILL_FAIL check_fails unexpected_exception message_mismatch_fails)
bbrsim_add_abort_test(testSupport.abort_stock testSupport abort_stock BBRT03)

bbrsim_add_test_program(testHFSSData LABEL unit
  CASES parse_frequency mini_loads transmittance_mixing load_cap_first_row nearest_key
        direction_cdf zero_weight_rows coherent_cross_term_in_direction exit_position
        err_missing_files err_no_keys err_key_mismatch err_freq_column err_ephi1_row_mismatch
        err_bad_number real_data_smoke transmittance_coherent real_data_diagonal_key
        sidecar_checks_pass sidecar_c1_header sidecar_c1_bom sidecar_c2_key_off_grid sidecar_c3_far_field
        sidecar_c4_exit_rows sidecar_c5_outside_section sidecar_c5_zero_outside sidecar_real_data
  ENV "BBRSIM_TEST_DATA=${PROJECT_SOURCE_DIR}/data")
bbrsim_add_test_program(testSidecar LABEL unit
  CASES parses_fixture forward_compatible f1_f3_identity f4_frames f5_f10_conventions f13_modes invariant_physics_only field_types f12_fits_solid)

bbrsim_add_test_program(testCrackLibrary LABEL unit
  CASES selection_grid5 token_verbatim discovery_errors unreadable_root concurrent_first_load)

bbrsim_add_test_program(testWrapper LABEL component NEEDS_ENSDF
  CASES fold_six_azimuths polarization_filter frequency_recorded rng_budget snap_and_rotation
        fatal_bbr004_bbr005 mirror_equivariance side_face_entry_documented
        reflectance_decision reflectance_specular reflectance_crack_wall reflectance_table_interp)

bbrsim_add_test_program(testPhysicsList LABEL regression NEEDS_ENSDF
  CASES wrapped_once_at_stock_index registered_twice_is_fatal wrong_order_is_fatal event_loop_fatal
  ENV "BBRSIM_TEST_DATA=${PROJECT_SOURCE_DIR}/data")

# testPassthrough: stock vs wrapped boundary optics, byte for byte. The program
# takes <mode> <nEvents> <outfile>, so it is registered by hand: four producer
# runs (fixture setups) and compare_files steps that require them.
add_executable(testPassthrough testPassthrough.cc BBRTestSupport.hh)
target_link_libraries(testPassthrough PRIVATE BBRsim)
target_compile_options(testPassthrough PRIVATE -Wall -Wextra)
set(_pt_n 2000)
set(_pt_wd ${CMAKE_CURRENT_BINARY_DIR}/work/passthrough)
set(_pt_ensdf "G4ENSDFSTATEDATA=${Geant4_DATASET_G4ENSDFSTATE_PATH}")
file(MAKE_DIRECTORY ${_pt_wd})
foreach(m stock wrapped stock2 mutant)
  set(mode ${m})
  if(m STREQUAL stock2)
    set(mode stock)
  endif()
  add_test(NAME passthrough.${m} COMMAND testPassthrough ${mode} ${_pt_n} pt_${m}.txt)
  set_tests_properties(passthrough.${m} PROPERTIES FIXTURES_SETUP pt_${m} ENVIRONMENT "${_pt_ensdf}")
endforeach()
add_test(NAME passthrough.identical COMMAND ${CMAKE_COMMAND} -E compare_files pt_stock.txt pt_wrapped.txt)
add_test(NAME passthrough.repeatable COMMAND ${CMAKE_COMMAND} -E compare_files pt_stock.txt pt_stock2.txt)
add_test(NAME passthrough.mutant_differs COMMAND ${CMAKE_COMMAND} -E compare_files pt_stock.txt pt_mutant.txt)
set_tests_properties(passthrough.identical PROPERTIES FIXTURES_REQUIRED "pt_stock;pt_wrapped")
set_tests_properties(passthrough.repeatable PROPERTIES FIXTURES_REQUIRED "pt_stock;pt_stock2")
set_tests_properties(passthrough.mutant_differs PROPERTIES FIXTURES_REQUIRED "pt_stock;pt_mutant" WILL_FAIL TRUE)
# With WLS active in box B the wrapper must still be bit-identical (D1: the
# wrapper keeps the stock process's place in the process list).
foreach(m stock wrapped)
  add_test(NAME passthrough.wls_${m} COMMAND testPassthrough ${m} ${_pt_n} wls_${m}.txt)
  set_tests_properties(passthrough.wls_${m} PROPERTIES FIXTURES_SETUP wls_${m} ENVIRONMENT "PT_WLS=1;${_pt_ensdf}")
endforeach()
add_test(NAME passthrough.wls_identical COMMAND ${CMAKE_COMMAND} -E compare_files wls_stock.txt wls_wrapped.txt)
set_tests_properties(passthrough.wls_identical PROPERTIES FIXTURES_REQUIRED "wls_stock;wls_wrapped")
# Common properties; the runs without WLS also drop an inherited PT_WLS.
foreach(t stock wrapped stock2 mutant identical repeatable mutant_differs wls_stock wls_wrapped wls_identical)
  set(mod ${BBRSIM_TEST_HERMETIC})
  if(NOT t MATCHES "^wls_")
    list(APPEND mod "PT_WLS=unset:")
  endif()
  set_tests_properties(passthrough.${t} PROPERTIES LABELS regression TIMEOUT 120 WORKING_DIRECTORY ${_pt_wd}
    ENVIRONMENT_MODIFICATION "${mod}")
endforeach()

bbrsim_add_test_program(testMaterials LABEL unit
  CASES cu_table_shape cu_drude_grid cu_limits cu_serov_of cu_offgrid_interp cu_cache cu_temperature_identity cu_invalid_input
        dielectrics dielectric_offgrid_interp dielectric_name_collision flag_materials
  ENV "G4LEDATA=/nonexistent;G4LEVELGAMMADATA=/nonexistent;G4NEUTRONHPDATA=/nonexistent;G4RADIOACTIVEDATA=/nonexistent;G4PARTICLEXSDATA=/nonexistent;G4PIIDATA=/nonexistent;G4REALSURFACEDATA=/nonexistent;G4SAIDXSDATA=/nonexistent;G4ABLADATA=/nonexistent;G4INCLDATA=/nonexistent;G4ENSDFSTATEDATA=/nonexistent;G4CHANNELINGDATA=/nonexistent")
bbrsim_add_test_program(testEmission LABEL unit
  CASES planck_cdf_shape planck_cdf_invalid planck_cdf_cold_band planck_narrow_band planck_failed_reinitialization planck_band_loss genevt_nonrandom_exact genevt_rotation genevt_statistics genevt_guards surface_input_guards zero_emissivity geometric_surface_area)
bbrsim_add_test_program(testConfig LABEL unit
  CASES defaults command_table datadir_quotes setT_units validation_paths datadir_default datadir_env datadir_env_late
        worker_clone_isolation worker_before_master print_flushes
  ENVMOD "BBRSIMDATA=unset:")
bbrsim_add_test_program(testBoundaryStatuses LABEL unit CASES stock_statuses)
target_compile_definitions(testConfig PRIVATE "BBRSIM_TEST_PREFIX=\"${CMAKE_INSTALL_PREFIX}\"")
set_tests_properties(testConfig.datadir_env testConfig.datadir_env_late PROPERTIES ENVIRONMENT_MODIFICATION "${BBRSIM_TEST_HERMETIC};BBRSIMDATA=set:/env/data")

# The event-loop fatal on a worker thread aborts the process; ExpectAbort.cmake
# checks the abort and the code. Same environment as the testPhysicsList cases.
bbrsim_add_abort_test(testPhysicsList.event_loop_fatal_mt testPhysicsList event_loop_fatal_mt BBR011)
set_tests_properties(testPhysicsList.event_loop_fatal_mt PROPERTIES LABELS regression
  ENVIRONMENT "BBRSIM_TEST_DATA=${PROJECT_SOURCE_DIR}/data;G4ENSDFSTATEDATA=${Geant4_DATASET_G4ENSDFSTATE_PATH}")
