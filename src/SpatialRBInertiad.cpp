/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/Constants.h>
#include <Ravelin/FrameException.h>
#include <Ravelin/SpatialRBInertiad.h>

using namespace Ravelin;

#define REAL double
#define SPATIAL_RB_INERTIA SpatialRBInertiad
#define TWIST Twistd
#define WRENCH Wrenchd
#define MATRIX3 Matrix3d
#define VECTOR3 Vector3d 

#include "SpatialRBInertia.cpp"

#undef REAL
#undef SPATIALRBINERTIA
#undef TWIST
#undef WRENCH
#undef MATRIX3
#undef VECTOR3
