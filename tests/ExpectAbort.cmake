# ExpectAbort.cmake — passes when a command dies in a fatal G4Exception.
#   cmake -DCMD=<program> -DARGS=<case> -DCODE=<code> -P ExpectAbort.cmake
# The program must abort (a fatal G4Exception with the stock handler prints
# "*** G4Exception: Aborting execution ***" and calls abort()) and its output
# must name the code. Any other exit fails, a crash included: a fatal turned
# into a warning that then crashes would name the code as well. CTest's
# WILL_FAIL and PASS_REGULAR_EXPRESSION treat a signal as a failure, so this
# driver does it.
execute_process(COMMAND ${CMD} ${ARGS}
                RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
message(STATUS "exit status: '${rc}'")
string(FIND "${out}${err}" "*** G4Exception: Aborting execution ***" banner)
if(NOT rc STREQUAL "Subprocess aborted" AND banner EQUAL -1)
  message(FATAL_ERROR "expected an abort, got exit status '${rc}' and no abort banner:\n${out}${err}")
endif()
string(FIND "${out}${err}" "G4Exception : ${CODE}" pos)
if(pos EQUAL -1)
  message(FATAL_ERROR "no 'G4Exception : ${CODE}' in the output:\n${out}${err}")
endif()
message(STATUS "aborted with ${CODE} as expected")
