/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SVECTOR6F_H
#define _SVECTOR6F_H

#include <boost/shared_ptr.hpp>
#include <Ravelin/MissizeException.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/fIterator.h>

namespace Ravelin {

#define REAL float
#define SVECTOR6 SVector6f
#define VECTOR3 Vector3f
#define ITERATOR fIterator
#define CONST_ITERATOR fIterator_const

#include "SVector6.h"

#undef REAL
#undef SVECTOR6
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif

