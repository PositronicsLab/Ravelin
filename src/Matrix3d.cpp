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
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Operators.h>
#include <Ravelin/Matrix3d.h>

using namespace Ravelin;

#define MATRIX3 Matrix3d
#define VECTOR3 Vector3d 
#define REAL double
#define AANGLE AAngled
#define QUAT Quatd 
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "Matrix3.cpp"

#undef MATRIX3
#undef VECTOR3
#undef REAL
#undef AANGLE
#undef QUAT
#undef ITERATOR
#undef CONST_ITERATOR

