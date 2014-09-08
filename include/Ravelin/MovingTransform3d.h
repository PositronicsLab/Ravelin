/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOVING_TRANSFORM3D_H
#define _MOVING_TRANSFORM3D_H

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Origin3d.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/SVelocityd.h>
#include <Ravelin/SAcceld.h>

namespace Ravelin {

#include "ddefs.h"
#include "MovingTransform3.h"
#include "undefs.h"

} // end namespace

#endif

