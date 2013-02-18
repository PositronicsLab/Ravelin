/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTOR3F_H
#define _VECTOR3F_H

#include <assert.h>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <iostream>
#include <Ravelin/fIterator.h>

namespace Ravelin {

#define VECTOR3 Vector3f
#define MATRIX3 Matrix3f 
#define REAL float
#define CONST_ITERATOR fIterator_const
#define ITERATOR fIterator

#include "Vector3.h"

#undef VECTOR3
#undef MATRIX3
#undef REAL
#undef CONST_ITERATOR
#undef ITERATOR

} // end namespace

#endif

