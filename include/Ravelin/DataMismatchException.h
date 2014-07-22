/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _DATA_MISMATCH_EXCEPTION_H_
#define _DATA_MISMATCH_EXCEPTION_H_

#include <stdexcept>

namespace Ravelin {

/// Exception thrown when trying to initialize a fixed size vector/matrix with the wrong size
class DataMismatchException : public std::runtime_error
{
  public:
    DataMismatchException() : std::runtime_error("float / double incompatibility") {}
}; // end class


} // end namespace

#endif

