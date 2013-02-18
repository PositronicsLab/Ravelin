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
#include <Ravelin/Operators.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>

using namespace Ravelin;

#define MATRIX3 Matrix3f
#define REAL float
#define AANGLE AAnglef
#define QUAT Quatf 
#define VECTOR3 Vector3f 
#define ITERATOR fIterator 
#define CONST_ITERATOR fIterator_const 

#include "Matrix3.cpp"

#undef MATRIX3
#undef REAL
#undef AANGLE
#undef QUAT
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR

