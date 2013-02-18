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
#include <Ravelin/AAngled.h>
#include <Ravelin/Quatd.h>

using namespace Ravelin;

#define QUAT Quatd
#define REAL double
#define EPS EPS_DOUBLE
#define MATRIX3 Matrix3d
#define AANGLE AAngled
#define VECTOR3 Vector3d 

#include "Quat.cpp"

#undef QUAT
#undef REAL
#undef EPS
#undef MATRIX3
#undef AANGLE
#undef VECTOR3

