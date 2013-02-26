/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
const Vector3f ZEROS_3F = Vector3f::zero();
//const Matrix3f ZEROS_3x3F = Matrix3f::zero();
//const Matrix3f IDENTITY_3x3F = Matrix3f::identity();
//const VectorNf EMPTY_VECF(0);
const Vector3d ZEROS_3D = Vector3d::zero();
//const Matrix3d ZEROS_3x3D = Matrix3d::zero();
//const Matrix3d IDENTITY_3x3D = Matrix3d::identity();
//const VectorNd EMPTY_VECD(0);

}

#endif
