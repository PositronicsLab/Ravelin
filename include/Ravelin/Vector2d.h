/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR2D_H
#define _VECTOR2D_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <Ravelin/dIterator.h>

namespace Ravelin {

#define REAL double
#define VECTOR2 Vector2d
#define MATRIX2 Matrix2d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "Vector2.h"

#undef REAL
#undef VECTOR2
#undef MATRIX2
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif

