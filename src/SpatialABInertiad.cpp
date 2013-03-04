/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Constants.h>
#include <Ravelin/SpatialABInertiad.h>

using namespace Ravelin;

#define SPATIAL_AB_INERTIA SpatialABInertiad
#define SPATIAL_RB_INERTIA SpatialRBInertiad
#define MATRIX3 Matrix3d
#define VECTOR3 Vector3d
#define TWIST Twistd
#define WRENCH Wrenchd
#define REAL double

#include "SpatialABInertia.cpp"

#undef SPATIAL_AB_INERTIA 
#undef SPATIAL_RB_INERTIA 
#undef MATRIX3 
#undef VECTOR3 
#undef TWIST 
#undef WRENCH
#undef REAL


