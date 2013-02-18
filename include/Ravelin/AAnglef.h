/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _AXIS_ANGLEF_H
#define _AXIS_ANGLEF_H

#include <Ravelin/Constants.h>
#include <Ravelin/Matrix3f.h>

namespace Ravelin {

#define REAL float 
#define EPS EPS_FLOAT
#define AANGLE AAnglef
#define VECTOR3 Vector3f
#define QUAT Quatf
#define MATRIX3 Matrix3f 

#include "AAngle.h"

#undef REAL 
#undef EPS 
#undef AANGLE 
#undef VECTOR3 
#undef QUAT 
#undef MATRIX3

} // end namespace

#endif

