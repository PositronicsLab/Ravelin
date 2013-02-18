/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_RB_INERTIAD_H
#define _SPATIAL_RB_INERTIAD_H

#include <boost/shared_ptr.hpp>
#include <Moby/Wrenchd.h>
#include <Moby/Twistd.h>
#include <Moby/Posed.h>

namespace Moby {

#define REAL double
#define SPATIAL_RB_INERTIA SpatialRBInertiad
#define WRENCH Wrenchd
#define TWIST Twistd
#define VECTOR3 Vector3d
#define MATRIX3 Matrix3d
#define POSE Posed

#include "SpatialRBInertia.h"

#undef REAL 
#undef SPATIAL_RB_INERTIA 
#undef WRENCH 
#undef TWIST 
#undef VECTOR3 
#undef MATRIX3 
#undef POSE

} // end namespace

#endif

