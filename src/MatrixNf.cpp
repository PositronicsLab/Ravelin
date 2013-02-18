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
#include <Ravelin/MatrixNf.h>

using std::vector;
using namespace Ravelin;
using boost::shared_array;

#define MATRIXN MatrixNf
#define REAL float
#define VECTORN VectorNf
#define MATRIX3 Matrix3f
#define ITERATOR fIterator 
#define CONST_ITERATOR fIterator_const 

#include "MatrixN.cpp"

#undef MATRIXN
#undef REAL
#undef VECTORN
#undef MATRIX3
#undef ITERATOR
#undef CONST_ITERATOR

