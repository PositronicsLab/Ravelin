# setup some things initially
cmake_minimum_required (VERSION 2.6)
project (Ravelin)
include (CheckIncludeFiles)
include (CheckLibraryExists)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

# setup include directories
include_directories ("include")

# setup library sources
set (SOURCES AAnglef.cpp AAngled.cpp ArticulatedBodyf.cpp ArticulatedBodyd.cpp cblas.cpp CRBAlgorithmd.cpp CRBAlgorithmf.cpp FixedJointd.cpp FixedJointf.cpp FSABAlgorithmd.cpp FSABAlgorithmf.cpp Jointd.cpp Jointf.cpp LinAlgf.cpp LinAlgd.cpp Log.cpp Matrix2d.cpp Matrix2f.cpp Matrix3d.cpp Matrix3f.cpp MatrixNf.cpp MatrixNd.cpp MovingTransform3f.cpp MovingTransform3d.cpp Origin2d.cpp Origin2f.cpp Origin3d.cpp Origin3f.cpp Pose2d.cpp Pose2f.cpp Pose3f.cpp Pose3d.cpp Quatf.cpp Quatd.cpp PrismaticJointf.cpp PrismaticJointd.cpp RCArticulatedBodyf.cpp RCArticulatedBodyd.cpp RevoluteJointf.cpp RevoluteJointd.cpp RNEAlgorithmf.cpp RNEAlgorithmd.cpp SpatialArithmeticd.cpp SpatialArithmeticf.cpp RigidBodyf.cpp RigidBodyd.cpp SForcef.cpp SForced.cpp SharedMatrixNf.cpp SharedMatrixNd.cpp SharedVectorNf.cpp SharedVectorNd.cpp SingleBodyf.cpp SingleBodyd.cpp SMomentumf.cpp SMomentumd.cpp SparseMatrixNf.cpp SparseMatrixNd.cpp SparseVectorNf.cpp SparseVectorNd.cpp SpatialABInertiad.cpp SpatialABInertiaf.cpp SpatialRBInertiaf.cpp SpatialRBInertiad.cpp SphericalJointd.cpp SphericalJointf.cpp SVector6f.cpp SVector6d.cpp SVelocityd.cpp SVelocityf.cpp Transform2d.cpp Transform2f.cpp Transform3d.cpp Transform3f.cpp UniversalJointd.cpp UniversalJointf.cpp URDFReaderd.cpp URDFReaderf.cpp Vector2f.cpp Vector2d.cpp Vector3f.cpp Vector3d.cpp VectorNf.cpp VectorNd.cpp XMLTree.cpp)

# build options 
option (BUILD_SHARED_LIBS "Build Ravelin as a shared library?" ON)
option (PROFILE "Build for profiling?" OFF)
option (REENTRANT "Build Ravelin to be reentrant ? (slower)" OFF)
option (DISABLE_EXCEPT "Disable user-level exceptions for extra speed (not recommended)?" OFF)
option (BUILD_EXAMPLES "Build example program binaries?" ON)
option (BUILD_TESTS "Build test program binaries?" ON)

# modify C++ flags
if (REENTRANT)
  add_definitions (-DREENTRANT)
endif (REENTRANT)
if (PROFILE)
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-pg -g")
  set (CMAKE_CXX_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} "-pg -g")
  set (CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} "-pg -g")
  set (CMAKE_CXX_FLAGS_RELDEB ${CMAKE_C_FLAGS_RELDEB} "pg -g")
  set (CMAKE_CXX_FLAGS_MINSIZEREL ${CMAKE_C_FLAGS_MINSIZEREL} "pg -g")
endif (PROFILE)

# find libraries
find_package (Boost REQUIRED)
find_package (LibXml2 REQUIRED)
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (APPLE)
  find_package (BLAS REQUIRED)
  find_package (LAPACK REQUIRED)
else (APPLE)
  find_package (CBLAS REQUIRED)
  find_package (LAPACKLite REQUIRED)
  set (BLAS_LIBRARIES ${CBLAS_LIBRARIES})
endif (APPLE)

# setup include directories, compiler flags, and libraries for required pkgs
include_directories (${Boost_INCLUDE_DIR})

# setup include directories, compiler flags, and libraries for optional pkgs
if (LIBXML2_FOUND)
  set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${LIBXML2_DEFINITIONS})
  include_directories (${LIBXML2_INCLUDE_DIR})
endif (LIBXML2_FOUND)

# prepend "src/" to each source file
foreach (i ${SOURCES})
  set (LIBSOURCES ${LIBSOURCES} "${CMAKE_SOURCE_DIR}/src/${i}")
endforeach (i)

# create the library
add_library(Ravelin "" "" ${LIBSOURCES})
target_link_libraries (Ravelin ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LIBXML2_LIBRARIES} ${EXTRA_LIBRARIES})

# build examples
if (BUILD_EXAMPLES)
  add_executable(Ravelin-block example/block.cpp)
  add_executable(Ravelin-pendulum example/pendulum.cpp)
  add_executable(Ravelin-double-pendulum example/doublependulum.cpp)
  add_executable(Ravelin-urdf example/urdf.cpp)
  target_link_libraries(Ravelin-block Ravelin)
  target_link_libraries(Ravelin-pendulum Ravelin)
  target_link_libraries(Ravelin-double-pendulum Ravelin)
  target_link_libraries(Ravelin-urdf Ravelin)
endif (BUILD_EXAMPLES)

# build examples
if (BUILD_TESTS)
include_directories(test /usr/include/eigen3 include)
link_directories(${PROJECT_BINARY_DIR})
add_executable(RavelinTest test/LinearAlgebra.cpp test/BlockOperations.cpp test/Arithmetic.cpp test/TestUtils.cpp)
target_link_libraries(RavelinTest Ravelin gtest gtest_main)
endif (BUILD_TESTS)

# setup install locations
install (TARGETS Ravelin DESTINATION lib)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/Ravelin DESTINATION include)

