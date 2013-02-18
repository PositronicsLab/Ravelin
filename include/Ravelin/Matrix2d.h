/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX2D_H
#define _MATRIX2D_H

#include <vector>
#include <algorithm>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/Vector2d.h>

namespace Ravelin {

#define REAL double
#define MATRIX2 Matrix2d
#define VECTOR2 Vector2d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "Matrix2.h"

#undef REAL
#undef MATRIX2
#undef VECTOR2
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif
