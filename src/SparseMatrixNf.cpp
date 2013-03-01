/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/FastThreadable.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/SparseMatrixNf.h>
#include <Ravelin/MatrixNf.h>

using std::pair;
using boost::shared_array;
using std::map;
using std::make_pair;
using std::vector;

using namespace Ravelin;

#define REAL float
#define SPARSEMATRIXN SparseMatrixNf
#define SPARSEVECTORN SparseVectorNf
#define EPS EPS_FLOAT
#define MATRIXN MatrixNf
#define VECTORN VectorNf

#include "SparseMatrixN.cpp"

#undef REAL 
#undef SPARSEMATRIXN 
#undef SPARSEVECTORN 
#undef EPS
#undef MATRIXN 
#undef VECTORN 

