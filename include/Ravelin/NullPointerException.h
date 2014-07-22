/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _NULL_POINTER_EXCEPTION_H_
#define _NULL_POINTER_EXCEPTION_H_

#include <stdexcept>

namespace Ravelin {

/// Exception thrown when a pointer should be set but is not 
class NullPointerException : public std::runtime_error
{
  public:
    NullPointerException() : std::runtime_error("Pointer is not set") {}
    NullPointerException(const char* error) : std::runtime_error(error) {}
}; // end class


} // end namespace

#endif

