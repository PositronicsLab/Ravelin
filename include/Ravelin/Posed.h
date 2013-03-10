/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSED_H
#define _POSED_H

#include <boost/shared_ptr.hpp>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Wrenchd.h>
#include <Ravelin/Twistd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/SpatialABInertiad.h>

namespace Ravelin {

#define POSE Posed
#define REAL double
#define VECTOR3 Vector3d
#define QUAT Quatd
#define AANGLE AAngled
#define MATRIX3 Matrix3d
#define WRENCH Wrenchd
#define TWIST Twistd
#define SPATIAL_RB_INERTIA SpatialRBInertiad
#define SPATIAL_AB_INERTIA SpatialABInertiad

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

