/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/**
 * A bunch of typedefs to make things more readable.
 */

#ifndef _RAVELIN_TYPES_H
#define _RAVELIN_TYPES_H

#include <utility>
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>

namespace Ravelin {

/// reference frame type for reduced-coordinate dynamics computations
enum ReferenceFrameType { eGlobal, eLink, eLinkInertia, eLinkCOM, eJoint };

} // end namespace

#endif

