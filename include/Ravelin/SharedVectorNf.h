/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SHARED_VECTORNF_H_ 
#define _SHARED_VECTORNF_H_ 

#include <vector>
#include <boost/shared_array.hpp>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Operators.h>

namespace Ravelin {

#define SHAREDVECTORN SharedVectorNf
#define SHAREDMATRIXN SharedMatrixNf
#define EPS EPS_FLOAT
#define VECTOR3 Vector3f
#define VECTORN VectorNf
#define MATRIXN MatrixNf
#define CONST_ITERATOR fIterator_const
#define ITERATOR fIterator
#define REAL float

#include "SharedVectorN.h"

#undef SHAREDVECTORN
#undef SHAREDMATRIXN
#undef EPS
#undef VECTOR3
#undef VECTORN
#undef MATRIXN
#undef CONST_ITERATOR
#undef ITERATOR
#undef REAL

} // end namespace

#endif 

