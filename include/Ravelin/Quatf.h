/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RAVELIN_QUATF_H
#define _RAVELIN_QUATF_H

//#include <Ravelin/VectorNf.h>
#include <Ravelin/Vector3f.h>
//#include <Ravelin/MatrixNf.h>

namespace Ravelin {

class Matrix3f;
class Posef;
class AAnglef;

#define REAL float 
#define QUAT Quatf
#define POSE Posef
#define AANGLE AAnglef
#define VECTOR3 Vector3f
#define MATRIX3 Matrix3f
#define VECTORN VectorNf
#define MATRIXN MatrixNf

#include "Quat.h"

#undef REAL
#undef QUAT
#undef POSE
#undef AANGLE
#undef VECTOR3
#undef MATRIX3
#undef VECTORN
#undef MATRIXN

} // end namespace

#endif

