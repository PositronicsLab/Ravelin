/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <Ravelin/Operators.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Vector3f.h>

using namespace Ravelin;

#define VECTOR3 Vector3f
#define REAL float
#define ITERATOR fIterator
#define CONST_ITERATOR fIterator_const

#include "Vector3.cpp"

#undef VECTOR3
#undef REAL
#undef ITERATOR
#undef CONST_ITERATOR

