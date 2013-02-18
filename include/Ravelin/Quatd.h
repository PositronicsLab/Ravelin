/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RAVELIN_QUATD_H
#define _RAVELIN_QUATD_H

//#include <Ravelin/VectorNf.h>
#include <Ravelin/Vector3d.h>
//#include <Ravelin/MatrixNf.h>

namespace Ravelin {

#define REAL double 
#define QUAT Quatd 
#define VECTOR3 Vector3d
#define VECTORN VectorNd
#define MATRIXN MatrixNd
#define MATRIX3 Matrix3d
#define AANGLE AAngled 

#include "Quat.h"

#undef REAL
#undef QUAT
#undef VECTOR3
#undef VECTORN
#undef MATRIXN
#undef MATRIX3 
#undef AANGLE 

} // end namespace

#endif

