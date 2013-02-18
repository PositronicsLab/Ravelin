/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _VECTORNF_H
#define _VECTORNF_H

#include <cstdarg>
#include <boost/shared_array.hpp>
#include <Ravelin/cblas.h>
#include <Ravelin/Vector2f.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/SharedVectorNf.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>

namespace Ravelin {

#define REAL float
#define EPS EPS_FLOAT
#define VECTORN VectorNf
#define VECTOR2 Vector2f
#define VECTOR3 Vector3f
#define SHAREDVECTORN SharedVectorNf
#define CONST_ITERATOR fIterator_const 
#define ITERATOR fIterator

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

