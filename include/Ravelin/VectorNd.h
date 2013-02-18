/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTORND_H
#define _VECTORND_H

#include <cstdarg>
#include <boost/shared_array.hpp>
#include <Ravelin/cblas.h>
#include <Ravelin/Vector2d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SharedVectorNd.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>

namespace Ravelin {

#define REAL double 
#define EPS EPS_DOUBLE
#define VECTORN VectorNd
#define VECTOR2 Vector2d
#define VECTOR3 Vector3d
#define SHAREDVECTORN SharedVectorNd
#define CONST_ITERATOR dIterator_const 
#define ITERATOR dIterator

// include main vector functionality
#include "VectorN.h"

#undef REAL
#undef EPS
#undef VECTORN
#undef VECTOR2
#undef VECTOR3
#undef SHAREDVECTORN
#undef CONST_ITERATOR
#undef ITERATOR

} // end namespace

#endif

