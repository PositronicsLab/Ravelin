/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/FastThreadable.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/SparseMatrixNd.h>
#include <Ravelin/MatrixNd.h>

using std::pair;
using boost::shared_array;
using std::map;
using std::make_pair;
using std::vector;

using namespace Ravelin;

#define REAL double
#define SPARSEMATRIXN SparseMatrixNd
#define SPARSEVECTORN SparseVectorNd
#define EPS EPS_DOUBLE
#define MATRIXN MatrixNd
#define VECTORN VectorNd

#include "SparseMatrixN.cpp"

#undef REAL 
#undef SPARSEMATRIXN 
#undef SPARSEVECTORN 
#undef EPS
#undef MATRIXN 
#undef VECTORN 

