/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSED_H
#define _POSED_H

#include <Ravelin/Vector3d.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Matrix3d.h>

namespace Ravelin {

#define REAL double
#define VECTOR3 Vector3d
#define QUAT Quatd
#define AANGLE AAngled
#define MATRIX3 Matrix3d
#define POSE Posed

#include "Pose.h"

#undef REAL
#undef VECTOR3
#undef QUAT 
#undef AANGLE 
#undef MATRIX3 
#undef POSE

} // end namespace

#endif

