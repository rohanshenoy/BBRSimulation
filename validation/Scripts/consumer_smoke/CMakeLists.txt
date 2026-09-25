# Consumer smoke test: an external project that finds an installed BBRsim with
# nothing but CMAKE_PREFIX_PATH (plus the Apple Clang pin), links it and runs.
cmake_minimum_required(VERSION 3.16...3.27)
project(bbrsimConsumerSmoke LANGUAGES C CXX)
find_package(BBRsim REQUIRED)
add_executable(consumer_smoke main.cc)
target_link_libraries(consumer_smoke PRIVATE BBRsim::BBRsim)
