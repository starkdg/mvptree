cmake_minimum_required(VERSION 3.13)

project(mvptree VERSION 0.1.0.0
  DESCRIPTION "distance-based multiple vantage point tree data structure")

set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(MVPTREE_SRCS mvpnode.cpp mvptree.cpp)


add_library(mvptree STATIC ${MVPTREE_SRCS})
set_target_properties(mvptree PROPERTIES PUBLIC_HEADER mvptree.hpp)

add_executable(testmvpnode test_mvpnode.cpp)
target_compile_options(testmvpnode PUBLIC -g -O0 -Wall)
target_link_libraries(testmvpnode mvptree)

add_executable(testmvptree test_mvptree.cpp)
target_compile_options(testmvptree PUBLIC -g -O0 -Wall)
target_link_libraries(testmvptree mvptree)


include(CTest)
add_test(NAME testnodes COMMAND  testmvpnode)
add_test(NAME testtree  COMMAND testmvptree)


install(TARGETS mvptree ARCHIVE DESTINATION lib PUBLIC_HEADER DESTINATION include)



