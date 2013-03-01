/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_MATRIX_NF_H_
#define _SPARSE_MATRIX_NF_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Ravelin/SparseVectorNf.h>
#include <Ravelin/MatrixNf.h>
#include <Ravelin/VectorNf.h>

namespace Ravelin {

#define SPARSEMATRIXN SparseMatrixNf
#define SPARSEVECTORN SparseVectorNf
#define MATRIXN MatrixNf
#define VECTORN VectorNf
#define REAL float

#include "SparseMatrixN.h"

#undef SPARSEMATRIXN
#undef SPARSEVECTORN
#undef MATRIXN
#undef VECTORN
#undef REAL

} // end namespace

#endif

