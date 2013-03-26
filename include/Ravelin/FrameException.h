/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _FRAME_EXCEPTION_H_
#define _FRAME_EXCEPTION_H_

#include <stdexcept>

namespace Ravelin {

/// Exception thrown when trying to mix types from two different frames 
class FrameException : public std::runtime_error
{
  public:
    FrameException() : std::runtime_error("Frame mismatch") {}
}; // end class


} // end namespace

#endif

