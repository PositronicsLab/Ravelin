/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
#include <cstring>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <Ravelin/Constants.h>
#include <Ravelin/Operators.h>
#include <Ravelin/Matrix2f.h>

using namespace Ravelin;

#define MATRIX2 Matrix2f
#define VECTOR2 Vector2f
#define REAL float
#define EPS EPS_FLOAT
#define CONST_ITERATOR fIterator_const
#define ITERATOR fIterator
#include "Matrix2.cpp"
#undef MATRIX2
#undef VECTOR2
#undef REAL
#undef EPS
#undef CONST_ITERATOR
#undef ITERATOR

