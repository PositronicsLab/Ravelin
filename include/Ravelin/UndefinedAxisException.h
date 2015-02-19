/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 License 
 ****************************************************************************/

#ifndef _UNDEFINED_AXIS_EXCEPTION_H_
#define _UNDEFINED_AXIS_EXCEPTION_H_

#include <stdexcept>

namespace Ravelin {

/// Exception thrown when axis has not been defined but user calls operation that requires its definition 
class UndefinedAxisException : public std::runtime_error
{
  public:
    UndefinedAxisException() : std::runtime_error("Axis was not defined prior to operation that requires a defined access") {}
}; // end class


} // end namespace

#endif

