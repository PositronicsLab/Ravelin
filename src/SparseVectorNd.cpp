/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/MissizeException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/SparseVectorNd.h>

using std::map;
using boost::shared_array;
using namespace Ravelin;

#define EPS EPS_DOUBLE
#define REAL double
#define SPARSEVECTORN SparseVectorNd
#define VECTORN VectorNd

#include "SparseVectorN.cpp"

#undef EPS
#undef REAL
#undef SPARSEVECTORN
#undef VECTORN

