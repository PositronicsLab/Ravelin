/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RAVELIN_CONSTANTS_H
#define _RAVELIN_CONSTANTS_H

#include <limits>
#include <cmath>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Vector3d.h>

namespace Ravelin {

// enums
enum Transposition { eNoTranspose, eTranspose };

// constants
const double EPS_DOUBLE = std::sqrt(std::numeric_limits<double>::epsilon());
const float EPS_FLOAT = std::sqrt(std::numeric_limits<float>::epsilon());

}

#endif
