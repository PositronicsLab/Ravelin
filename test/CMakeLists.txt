cmake_minimum_required(VERSION 2.8.3)

project(Ravelin_Test)

find_package(Boost REQUIRED)

option(BUILD_SHARED_LIBS "Build Controller as a shared library?" ON)

# you may have to add additional includes in here...
add_definitions(-DBUILD_DOUBLE)
add_definitions(-DUSE_OSG)
add_definitions(-DSAFESTATIC=static)
add_definitions(-DUSE_GLPK)
include_directories(/usr/local/include ./ ${Boost_INCLUDE_DIR} /usr/include/eigen3 ../include /usr/local/Cellar/octave/3.6.2/include/octave-3.6.2)
link_directories(/usr/local/Cellar/octave/3.6.2/lib/octave/3.6.2 /usr/local/lib /usr/include/eigen3 ../debug/)

add_executable(${PROJECT_NAME} Arithmetic.cpp TestUtils.cpp BlockOperations.cpp LinearAlgebra.cpp)
target_link_libraries(${PROJECT_NAME} Ravelin)
