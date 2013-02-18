/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_RB_INERTIAF_H
#define _SPATIAL_RB_INERTIAF_H

#include <boost/shared_ptr.hpp>
#include <Moby/Wrenchf.h>
#include <Moby/Twistf.h>
#include <Moby/Posef.h>

namespace Moby {

#define REAL float
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define WRENCH Wrenchf
#define TWIST Twistf
#define VECTOR3 Vector3f
#define MATRIX3 Matrix3f
#define POSE Posef

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

