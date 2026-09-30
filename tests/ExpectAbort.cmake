# ExpectAbort.cmake — passes when a command dies in a fatal G4Exception.
#   cmake -DCMD=<program> -DARGS=<case> -DCODE=<code> -P ExpectAbort.cmake
# The program must exit non-zero (a fatal G4Exception with the stock handler
# calls abort()) and its output must name the code. CTest's WILL_FAIL and
# PASS_REGULAR_EXPRESSION treat a signal as a failure, so this driver does it.
execute_process(COMMAND ${CMD} ${ARGS}
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "exit status: '${rc}'")
if(rc STREQUAL "0")
  message(FATAL_ERROR "expected an abort, got exit 0")
endif()
string(FIND "${out}${err}" "G4Exception : ${CODE}" pos)
if(pos EQUAL -1)
  message(FATAL_ERROR "no 'G4Exception : ${CODE}' in the output:\n${out}${err}")
endif()
message(STATUS "aborted with ${CODE} as expected")
