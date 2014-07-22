/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _INVALID_TRANSFORM_EXCEPTION_H_
#define _INVALID_TRANSFORML_EXCEPTION_H_

#include <stdexcept>

namespace Ravelin {

/// Exception thrown when general numerical error occurs 
class InvalidTransformException : public std::runtime_error
{
  public:
    InvalidTransformException(const Matrix4& T) : std::runtime_error("Rotation component of matrix is invalid") {}
    InvalidTransformException(const MatrixN& T) : std::runtime_error("Rotation component of matrix is invalid, or bottom row of matrix is not 0 0 0 1") {}
}; // end class


} // end namespace

#endif

