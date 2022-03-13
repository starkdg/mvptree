cmake_minimum_required(VERSION 3.13)

project(mvptree VERSION 0.1.0.0
  DESCRIPTION "distance-based multiple vantage point tree data structure")

set(CMAKE_BUILD_TYPE RelWithDebInfo)

add_library(mvptree INTERFACE)
target_include_directories(mvptree INTERFACE include/)

add_executable(testmvpnode tests/test_mvpnode.cpp)
target_compile_options(testmvpnode PUBLIC -g -O0 -Wall -Wno-unused-variable)
target_link_libraries(testmvpnode mvptree)

add_executable(testmvptree tests/test_mvptree.cpp)
target_compile_options(testmvptree PUBLIC -g -O0 -Wall -Wno-unused-variable)
target_link_libraries(testmvptree mvptree)

add_executable(runmvptree tests/run_mvptree.cpp)
target_compile_options(runmvptree PUBLIC -g -O0 -Wall -Wno-unused-variable)
target_link_libraries(runmvptree mvptree)

include(CTest)
add_test(NAME testnodes COMMAND  testmvpnode)
add_test(NAME testtree  COMMAND testmvptree)


install(TARGETS mvptree PUBLIC_HEADER DESTINATION include)



