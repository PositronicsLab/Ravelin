/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
#include <cstring>
#include <list>
#include <cmath>
#include <iostream>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SharedVectorNd.h>
#include <Ravelin/SharedMatrixNd.h>

using namespace Ravelin;
using boost::shared_array;
using std::vector;

#define REAL double
#define MATRIXN MatrixNd
#define MATRIX3 Matrix3d
#define SHAREDVECTORN SharedVectorNd
#define CONST_SHAREDVECTORN SharedConstVectorNd
#define SHAREDMATRIXN SharedMatrixNd
#define CONST_SHAREDMATRIXN SharedConstMatrixNd
#define VECTORN VectorNd
#define VECTOR3 Vector3d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "SharedMatrixN.cpp"

#undef REAL
#undef MATRIXN
#undef MATRIX3
#undef SHAREDVECTORN
#undef SHAREDMATRIXN
#undef CONST_SHAREDVECTORN
#undef CONST_SHAREDMATRIXN
#undef VECTORN
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR

