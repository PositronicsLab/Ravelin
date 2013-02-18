/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR2F_H
#define _VECTOR2F_H

#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <Ravelin/fIterator.h>

namespace Ravelin {

#define REAL float 
#define VECTOR2 Vector2f
#define MATRIX2 Matrix2f
#define ITERATOR fIterator 
#define CONST_ITERATOR fIterator_const 

#include "Vector2.h"

#undef REAL
#undef VECTOR2
#undef MATRIX2
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif

