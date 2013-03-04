/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SVECTOR6D_H
#define _SVECTOR6D_H

#include <boost/shared_ptr.hpp>
#include <Ravelin/MissizeException.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/dIterator.h>

namespace Ravelin {

#define REAL double
#define SVECTOR6 SVector6d
#define VECTOR3 Vector3d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "SVector6.h"

#undef REAL
#undef SVECTOR6
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif

