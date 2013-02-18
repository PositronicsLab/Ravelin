/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _DITERATOR_H
#define _DITERATOR_H

#include <stdexcept>

namespace Ravelin {

#define REAL double 
#define CONST_ITERATOR dIterator_const
#define ITERATOR dIterator
#define MATRIXN MatrixNd
#define VECTORN VectorNd
#define VECTOR3 Vector3d
#define VECTOR2 Vector2d
#define MATRIX3 Matrix3d
#define MATRIX2 Matrix2d
#define SHAREDMATRIXN SharedMatrixNd
#define SHAREDVECTORN SharedVectorNd

#include "rIterator.h"

#undef REAL
#undef ITERATOR 
#undef MATRIXN
#undef VECTORN
#undef VECTOR3
#undef VECTOR2
#undef MATRIX3
#undef MATRIX2
#undef SHAREDMATRIXN
#undef SHAREDVECTORN

} // end namespace

#endif

