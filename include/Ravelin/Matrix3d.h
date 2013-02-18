/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX3D_H
#define _MATRIX3D_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <Ravelin/cblas.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Quatd.h>

namespace Ravelin {

#define REAL double 
#define QUAT Quatd
#define AANGLE AAngled
#define MATRIX3 Matrix3d
#define VECTOR3 Vector3d
#define ITERATOR dIterator 
#define CONST_ITERATOR dIterator_const

#include "Matrix3.h"

#undef REAL
#undef QUAT
#undef AANGLE
#undef MATRIX3
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif
