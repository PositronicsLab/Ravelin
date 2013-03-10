/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_AB_INERTIAF_H
#define _SPATIAL_AB_INERTIAF_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Wrenchf.h>
#include <Ravelin/Twistf.h>
#include <Ravelin/SpatialRBInertiaf.h>
#include <Ravelin/Matrix3f.h>

namespace Ravelin {

#define REAL float
#define SPATIAL_AB_INERTIA SpatialABInertiaf
#define SPATIAL_RB_INERTIA SpatialRBInertiaf
#define WRENCH Wrenchf
#define TWIST Twistf
#define MATRIX3 Matrix3f
#define POSE Posef

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

