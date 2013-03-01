/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_MATRIX_ND_H_
#define _SPARSE_MATRIX_ND_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Ravelin/SparseVectorNd.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>

namespace Ravelin {

#define SPARSEMATRIXN SparseMatrixNd
#define SPARSEVECTORN SparseVectorNd
#define MATRIXN MatrixNd
#define VECTORN VectorNd
#define REAL double

#include "SparseMatrixN.h"

#undef SPARSEMATRIXN
#undef SPARSEVECTORN
#undef MATRIXN
#undef VECTORN
#undef REAL

} // end namespace

#endif

