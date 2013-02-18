/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _AXIS_ANGLED_H
#define _AXIS_ANGLED_H

#include <Ravelin/Constants.h>
#include <Ravelin/Matrix3d.h>

namespace Ravelin {

#define REAL double
#define EPS EPS_DOUBLE
#define AANGLE AAngled
#define VECTOR3 Vector3d
#define QUAT Quatd
#define MATRIX3 Matrix3d 

#include "AAngle.h"

#undef REAL 
#undef EPS 
#undef AANGLE 
#undef VECTOR3 
#undef QUAT
#undef MATRIX3

} // end namespace

#endif

