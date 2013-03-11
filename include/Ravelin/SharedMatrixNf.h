/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SHARED_MATRIXNF_H
#define _SHARED_MATRIXNF_H

#include <Ravelin/cblas.h>
#include <Ravelin/VectorNf.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/DataMismatchException.h>
//#include <Ravelin/Posef.h>

namespace Ravelin {

#define REAL float 
#define MATRIXN MatrixNf
#define MATRIX3 Matrix3f
#define MATRIX2 Matrix2f
#define SHAREDVECTORN SharedVectorNf
#define CONST_SHAREDVECTORN SharedConstVectorNf
#define SHAREDMATRIXN SharedMatrixNf
#define CONST_SHAREDMATRIXN SharedConstMatrixNf
#define VECTORN VectorNf
#define POSE Posef
#define CONST_ITERATOR fIterator_const
#define ITERATOR fIterator

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

