/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/cblas.h>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <Ravelin/Constants.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Posed.h>

using namespace Ravelin;

#define REAL double
#define ZEROS_3 ZEROS_3D
#define POSE Posed
#define MATRIX3 Matrix3d
#define VECTOR3 Vector3d
#define QUAT Quatd
#define AANGLE AAngled

#include "Pose.cpp"

#undef REAL
#undef ZEROS_3
#undef POSE
#undef MATRIX3
#undef VECTOR3
#undef QUAT
#undef AANGLE

