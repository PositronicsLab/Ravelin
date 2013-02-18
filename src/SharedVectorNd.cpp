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
#include <Ravelin/VectorNd.h>
#include <Ravelin/SharedVectorNd.h>

using namespace Ravelin;
using boost::shared_array;
using std::vector;

#define REAL double 
#define SHAREDVECTORN SharedVectorNd
#define VECTORN VectorNd
#define VECTOR3 Vector3d
#include "SharedVectorN.cpp"
#undef REAL
#undef SHAREDVECTORN
#undef VECTORN
#undef VECTOR3

