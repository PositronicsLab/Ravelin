/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <assert.h>
#include <cmath>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Operators.h>
#include <Ravelin/MissizeException.h>

using namespace Ravelin;

#define REAL float
#define EPS EPS_FLOAT
#define AANGLE AAnglef
#define QUAT Quatf
#define VECTOR3 Vector3f
#define MATRIX3 Matrix3f

#include "AAngle.cpp"

#undef REAL
#undef EPS
#undef AANGLE
#undef QUAT
#undef VECTOR3
#undef MATRIX3

