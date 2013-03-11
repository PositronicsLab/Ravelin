/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <functional>
#include <algorithm>
#include <Ravelin/cblas.h>
#include <Ravelin/FrameException.h>
#include <Ravelin/SVector6d.h>

using namespace Ravelin;

#define SVECTOR6 SVector6d
#define REAL double
#define VECTOR3 Vector3d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const

#include "SVector6.cpp"

#undef SVECTOR6
#undef REAL
#undef VECTOR3
#undef ITERATOR
#undef CONST_ITERATOR


