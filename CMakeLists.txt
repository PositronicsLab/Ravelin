# setup some things initially
cmake_minimum_required (VERSION 2.6)
project (Ravelin)
include (CheckIncludeFiles)
include (CheckLibraryExists)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

# setup include directories
include_directories ("include")

# setup library sources
set (SOURCES AAnglef.cpp AAngled.cpp cblas.cpp Jointd.cpp Jointf.cpp UniversalJointd.cpp UniversalJointf.cpp SphericalJointd.cpp SphericalJointf.cpp FixedJointd.cpp FixedJointf.cpp RevoluteJointf.cpp RevoluteJointd.cpp PrismaticJointf.cpp PrismaticJointd.cpp LinAlgf.cpp LinAlgd.cpp Matrix2d.cpp Matrix2f.cpp Matrix3d.cpp Matrix3f.cpp MatrixNf.cpp MatrixNd.cpp MovingTransform3f.cpp MovingTransform3d.cpp Origin2d.cpp Origin2f.cpp Origin3d.cpp Origin3f.cpp Pose2d.cpp Pose2f.cpp Pose3f.cpp Pose3d.cpp Quatf.cpp Quatd.cpp SForcef.cpp SForced.cpp SharedMatrixNf.cpp SharedMatrixNd.cpp SharedVectorNf.cpp SharedVectorNd.cpp SMomentumf.cpp SMomentumd.cpp SparseMatrixNf.cpp SparseMatrixNd.cpp SparseVectorNf.cpp SparseVectorNd.cpp SpatialABInertiad.cpp SpatialABInertiaf.cpp SpatialRBInertiaf.cpp SpatialRBInertiad.cpp SVector6f.cpp SVector6d.cpp SVelocityd.cpp SVelocityf.cpp Transform2d.cpp Transform2f.cpp Transform3d.cpp Transform3f.cpp Vector2f.cpp Vector2d.cpp Vector3f.cpp Vector3d.cpp VectorNf.cpp VectorNd.cpp)

# build options 
option (BUILD_SHARED_LIBS "Build Ravelin as a shared library?" ON)
option (PROFILE "Build for profiling?" OFF)
option (REENTRANT "Build Ravelin to be reentrant ? (slower)" OFF)
option (DISABLE_EXCEPT "Disable user-level exceptions for extra speed (not recommended)?" OFF)

# superLU directories
set (SUPERLU_PATH /usr/local)

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

# prepend "src/" to each source file
foreach (i ${SOURCES})
  set (LIBSOURCES ${LIBSOURCES} "${CMAKE_SOURCE_DIR}/src/${i}")
endforeach (i)

# create the library
add_library(Ravelin "" "" ${LIBSOURCES})
target_link_libraries (Ravelin ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${EXTRA_LIBRARIES})

# setup install locations
install (TARGETS Ravelin DESTINATION lib)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/include/Ravelin DESTINATION include)

