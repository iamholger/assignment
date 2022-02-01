cmake_minimum_required(VERSION 3.14)
project(my_project)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  assignment
  main.cxx
)

add_executable(
  test_matlib
  test_matlib.cxx
)

target_link_libraries(
  test_matlib
  gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_matlib)
