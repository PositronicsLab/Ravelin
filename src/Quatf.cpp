/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Ravelin/Constants.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Quatf.h>

using namespace Ravelin;

#define QUAT Quatf
#define REAL float
#define EPS EPS_FLOAT
#define MATRIX3 Matrix3f
#define AANGLE AAnglef
#define VECTOR3 Vector3f 

#include "Quat.cpp"

#undef QUAT
#undef REAL
#undef EPS
#undef MATRIX3
#undef AANGLE
#undef VECTOR3

