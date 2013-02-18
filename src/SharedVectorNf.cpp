/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
#include <cstring>
#include <list>
#include <cmath>
#include <iostream>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/VectorNf.h>
#include <Ravelin/SharedVectorNf.h>

using namespace Ravelin;
using boost::shared_array;
using std::vector;

#define REAL float
#define SHAREDVECTORN SharedVectorNf
#define VECTORN VectorNf
#define VECTOR3 Vector3f
#include "SharedVectorN.cpp"
#undef REAL
#undef SHAREDVECTORN
#undef VECTORN
#undef VECTOR3

