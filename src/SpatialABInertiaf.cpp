/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Constants.h>
#include <Ravelin/SpatialABInertiaf.h>

using namespace Ravelin;

#define SPATIAL_AB_INERTIA SpatialABInertiaf
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define MATRIX3 Matrix3f
#define VECTOR3 Vector3f
#define TWIST Twistf
#define WRENCH Wrenchf
#define REAL float

#include "SpatialABInertia.cpp"

#undef SPATIAL_AB_INERTIA 
#undef SPATIAL_RB_INERTIA 
#undef MATRIX3 
#undef VECTOR3 
#undef TWIST 
#undef WRENCH
#undef REAL


