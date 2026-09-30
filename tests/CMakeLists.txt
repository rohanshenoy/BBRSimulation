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
        check_fails unexpected_exception
  WILL_FAIL check_fails unexpected_exception)
bbrsim_add_abort_test(testSupport.abort_stock testSupport abort_stock BBRT03)

bbrsim_add_test_program(testHFSSData LABEL unit
  CASES parse_frequency mini_loads transmittance_mixing load_cap_first_row nearest_key
        direction_cdf zero_weight_rows coherent_cross_term_in_direction exit_position
        err_missing_files err_no_keys err_key_mismatch err_freq_column err_ephi1_row_mismatch
        err_bad_number real_data_smoke
  ENV "BBRSIM_TEST_DATA=${PROJECT_SOURCE_DIR}/data")

bbrsim_add_test_program(testCrackLibrary LABEL unit
  CASES selection_grid5 token_verbatim discovery_errors unreadable_root concurrent_first_load)
