/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <assert.h>
#include <cmath>
#include <Ravelin/AAngled.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Operators.h>
#include <Ravelin/MissizeException.h>

using namespace Ravelin;

#define REAL double 
#define EPS EPS_DOUBLE
#define AANGLE AAngled
#define QUAT Quatd
#define VECTOR3 Vector3d
#define MATRIX3 Matrix3d

#include "AAngle.cpp"

#undef REAL
#undef EPS
#undef AANGLE
#undef QUAT
#undef VECTOR3
#undef MATRIX3

