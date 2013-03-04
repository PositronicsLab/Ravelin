/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_AB_INERTIAD_H
#define _SPATIAL_AB_INERTIAD_H

#include <boost/shared_ptr.hpp>
#include <Ravelin/Wrenchd.h>
#include <Ravelin/Twistd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Pose.h>

namespace Ravelin {

#define REAL double 
#define SPATIAL_AB_INERTIA SpatialABInertiad
#define SPATIAL_RB_INERTIA SpatialRBInertiad
#define WRENCH Wrenchd
#define TWIST Twistd
#define MATRIX3 Matrix3d
#define POSE Posed

#include "SpatialABInertia.h"

#undef REAL
#undef SPATIAL_AB_INERTIA
#undef SPATIAL_RB_INERTIA
#undef WRENCH
#undef TWIST
#undef MATRIX3
#undef POSE

} // end namespace

#endif

