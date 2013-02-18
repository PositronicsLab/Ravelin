/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Matrix2f.h>
#include <Ravelin/Vector2f.h>

using namespace Ravelin;

#define VECTOR2 Vector2f
#define REAL float
#define ITERATOR fIterator
#define CONST_ITERATOR fIterator_const

#include "Vector2.cpp"

#undef VECTOR2
#undef REAL
#undef ITERATOR
#undef CONST_ITERATOR

