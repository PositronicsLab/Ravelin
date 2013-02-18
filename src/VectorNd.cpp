/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <list>
#include <Ravelin/Vector2d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>

using namespace Ravelin;

#define REAL double 
#define VECTORN VectorNd
#define VECTOR2 Vector2d
#define VECTOR3 Vector3d
#define SHAREDVECTORN SharedVectorNd
#include "VectorN.cpp"
#undef REAL
#undef VECTORN
#undef VECTOR2
#undef VECTOR3
#undef SHAREDVECTORN

