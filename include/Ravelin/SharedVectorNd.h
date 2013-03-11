/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SHARED_VECTORND_H_ 
#define _SHARED_VECTORND_H_ 

#include <vector>
#include <boost/shared_array.hpp>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Operators.h>

namespace Ravelin {

#define SHAREDVECTORN SharedVectorNd
#define SHAREDMATRIXN SharedMatrixNd
#define CONST_SHAREDVECTORN SharedConstVectorNd
#define CONST_SHAREDMATRIXN SharedConstMatrixNd
#define EPS EPS_DOUBLE
#define VECTOR3 Vector3d
#define VECTORN VectorNd
#define MATRIXN MatrixNd
#define CONST_ITERATOR dIterator_const
#define ITERATOR dIterator
#define REAL double

#include "SharedVectorN.h"

#undef SHAREDVECTORN
#undef SHAREDMATRIXN
#undef CONST_SHAREDVECTORN
#undef CONST_SHAREDMATRIXN
#undef EPS
#undef VECTOR3
#undef VECTORN
#undef MATRIXN
#undef CONST_ITERATOR
#undef ITERATOR
#undef REAL

} // end namespace

#endif 

