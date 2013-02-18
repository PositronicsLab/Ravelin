/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cstring>
#include <list>
#include <Ravelin/Vector2f.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/VectorNf.h>

using namespace Ravelin;

#define REAL float
#define VECTORN VectorNf
#define VECTOR2 Vector2f
#define VECTOR3 Vector3f
#define SHAREDVECTORN SharedVectorNf
#include "VectorN.cpp"
#undef REAL
#undef VECTORN
#undef VECTOR2
#undef VECTOR3
#undef SHAREDVECTORN

