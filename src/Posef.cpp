/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/Posef.h>

using namespace Ravelin;

#define REAL float
#define POSE Posef
#define MATRIX3 Matrix3f
#define VECTOR3 Vector3f
#define QUAT Quatf
#define AANGLE AAnglef
#define WRENCH Wrenchf
#define TWIST Twistf
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define SPATIAL_AB_INERTIA SpatialABInertiaf 

#include "Pose.cpp"

#undef REAL
#undef ZEROS_3
#undef POSE
#undef MATRIX3
#undef VECTOR3
#undef QUAT
#undef AANGLE
#undef WRENCH
#undef TWIST
#undef SPATIAL_RB_INERTIA
#undef SPATIAL_AB_INERTIA

