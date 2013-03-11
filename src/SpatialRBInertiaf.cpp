/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Constants.h>
#include <Ravelin/FrameException.h>
#include <Ravelin/SpatialRBInertiaf.h>

using namespace Ravelin;

#define REAL float
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define TWIST Twistf
#define WRENCH Wrenchf
#define MATRIX3 Matrix3f
#define VECTOR3 Vector3f 

#include "SpatialRBInertia.cpp"

#undef REAL
#undef SPATIAL_RB_INERTIA
#undef TWIST
#undef WRENCH
#undef MATRIX3
#undef VECTOR3
