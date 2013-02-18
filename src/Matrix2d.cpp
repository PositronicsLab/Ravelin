/****************************************************************************
 * Copyright 2013 Evan Drumwright
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
#include <Ravelin/Matrix2d.h>

using namespace Ravelin;

#define MATRIX2 Matrix2d
#define VECTOR2 Vector2d
#define REAL double
#define EPS EPS_DOUBLE
#define CONST_ITERATOR dIterator_const
#define ITERATOR dIterator
#include "Matrix2.cpp"
#undef MATRIX2
#undef VECTOR2
#undef REAL
#undef EPS
#undef CONST_ITERATOR
#undef ITERATOR

