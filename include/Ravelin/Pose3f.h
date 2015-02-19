/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _POSE3F_H
#define _POSE3F_H

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Origin3f.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/SForcef.h>
#include <Ravelin/SMomentumf.h>
#include <Ravelin/SAccelf.h> 
#include <Ravelin/SVelocityf.h>
#include <Ravelin/SpatialRBInertiaf.h>
#include <Ravelin/SpatialABInertiaf.h>
#include <Ravelin/Transform3f.h>

namespace Ravelin {

#include "fdefs.h"
#include "Pose3.h"
#include "undefs.h"

} // end namespace

#endif

