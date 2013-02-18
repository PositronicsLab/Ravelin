/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR3D_H
#define _VECTOR3D_H

#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <Ravelin/dIterator.h>

namespace Ravelin {

#define VECTOR3 Vector3d
#define MATRIX3 Matrix3d 
#define REAL double 
#define CONST_ITERATOR dIterator_const
#define ITERATOR dIterator

#include "Vector3.h"

#undef VECTOR3
#undef MATRIX3
#undef REAL
#undef CONST_ITERATOR
#undef ITERATOR

} // end namespace

#endif

