/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <cmath>
#include <iomanip>
#include <Ravelin/FastThreadable.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/MatrixNd.h>

using std::vector;
using namespace Ravelin;
using boost::shared_array;

#define MATRIXN MatrixNd
#define REAL double
#define VECTORN VectorNd
#define MATRIX3 Matrix3d
#define ITERATOR dIterator 
#define CONST_ITERATOR dIterator_const 

#include "MatrixN.cpp"

#undef MATRIXN
#undef REAL
#undef VECTORN
#undef MATRIX3
#undef ITERATOR
#undef CONST_ITERATOR

