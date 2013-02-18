/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <Ravelin/Operators.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Vector3d.h>

using namespace Ravelin;

#define VECTOR3 Vector3d
#define REAL double
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "Vector3.cpp"

#undef VECTOR3
#undef REAL
#undef ITERATOR
#undef CONST_ITERATOR

