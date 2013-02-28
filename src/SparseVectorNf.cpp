/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/MissizeException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/SparseVectorNf.h>

using std::map;
using boost::shared_array;
using namespace Ravelin;

#define EPS EPS_FLOAT
#define REAL float
#define SPARSEVECTORN SparseVectorNf
#define VECTORN VectorNf

#include "SparseVectorN.cpp"

#undef EPS
#undef REAL
#undef SPARSEVECTORN
#undef VECTORN

