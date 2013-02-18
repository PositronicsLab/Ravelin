/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX3F_H
#define _MATRIX3F_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Quatf.h>

namespace Ravelin {

#define REAL float
#define QUAT Quatf
#define AANGLE AAnglef
#define MATRIX3 Matrix3f
#define VECTOR3 Vector3f 
#define ITERATOR fIterator 
#define CONST_ITERATOR fIterator_const

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
