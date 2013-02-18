/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Vector2d.h>
#include <Ravelin/InvalidIndexException.h>

using namespace Ravelin;

#define VECTOR2 Vector2d
#define REAL double
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "Vector2.cpp"

#undef VECTOR2
#undef REAL
#undef ITERATOR
#undef CONST_ITERATOR

