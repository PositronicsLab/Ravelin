/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSEF_H
#define _POSEF_H

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/Wrenchf.h>
#include <Ravelin/Twistf.h>
#include <Ravelin/SpatialRBInertiaf.h>
#include <Ravelin/SpatialABInertiaf.h>

namespace Ravelin {

#define POSE Posef
#define REAL float
#define VECTOR3 Vector3f
#define QUAT Quatf
#define AANGLE AAnglef
#define MATRIX3 Matrix3f
#define WRENCH Wrenchf
#define TWIST Twistf
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define SPATIAL_AB_INERTIA SpatialABInertiaf

#include "Pose.h"

#undef POSE
#undef REAL 
#undef VECTOR3
#undef QUAT 
#undef AANGLE 
#undef MATRIX3 
#undef WRENCH
#undef TWIST
#undef SPATIAL_RB_INERTIA
#undef SPATIAL_AB_INERTIA

} // end namespace

#endif

