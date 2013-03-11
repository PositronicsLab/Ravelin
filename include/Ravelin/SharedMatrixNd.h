/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SHARED_MATRIXND_H
#define _SHARED_MATRIXND_H

#include <Ravelin/cblas.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/DataMismatchException.h>
//#include <Ravelin/Posef.h>

namespace Ravelin {

#define REAL double
#define MATRIXN MatrixNd
#define MATRIX3 Matrix3d
#define MATRIX2 Matrix2d
#define SHAREDVECTORN SharedVectorNd
#define CONST_SHAREDVECTORN SharedConstVectorNd
#define SHAREDMATRIXN SharedMatrixNd
#define CONST_SHAREDMATRIXN SharedConstMatrixNd
#define VECTORN VectorNd
#define POSE Posed
#define CONST_ITERATOR dIterator_const
#define ITERATOR dIterator

#include "SharedMatrixN.h"

#undef REAL 
#undef MATRIXN 
#undef MATRIX3
#undef MATRIX2
#undef SHAREDVECTORN
#undef CONST_SHAREDVECTORN
#undef SHAREDMATRIXN
#undef CONST_SHAREDMATRIXN
#undef VECTORN 
#undef POSE 
#undef CONST_ITERATOR
#undef ITERATOR

} // end namespace

#endif

