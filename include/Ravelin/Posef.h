/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSEF_H
#define _POSEF_H

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>

namespace Ravelin {

#define REAL float
#define VECTOR3 Vector3f
#define QUAT Quatf
#define AANGLE AAnglef
#define MATRIX3 Matrix3f
#define POSE Posef

#include "Pose.h"

#undef REAL 
#undef VECTOR3
#undef QUAT 
#undef AANGLE 
#undef MATRIX3 
#undef POSE

} // end namespace

#endif

