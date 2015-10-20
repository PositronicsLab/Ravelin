cmake_minimum_required(VERSION 2.8.3)

project(RavelinTest CXX C)

# you may have to add additional includes in here...
add_definitions(-DBUILD_DOUBLE)
#include_directories(/usr/local/include ./ ${Boost_INCLUDE_DIR} /usr/include/eigen3 ../include)
include_directories(/usr/local/include ./ /usr/include/eigen3 ../include)
link_directories(/usr/local/lib ../release/)

#add_executable(${PROJECT_NAME} Test.cpp Arithmetic.cpp TestUtils.cpp BlockOperations.cpp LinearAlgebra.cpp)
add_executable(RavelinTest LinearAlgebra.cpp BlockOperations.cpp Arithmetic.cpp TestUtils.cpp)
target_link_libraries(RavelinTest Ravelin gtest gtest_main)

