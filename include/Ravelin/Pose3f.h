/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSE3F_H
#define _POSE3F_H

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3f.h>
#include <Ravelin/Point3f.h>
#include <Ravelin/Origin3f.h>
#include <Ravelin/Quatf.h>
#include <Ravelin/AAnglef.h>
#include <Ravelin/Matrix3f.h>
#include <Ravelin/Wrenchf.h>
#include <Ravelin/Twistf.h>
#include <Ravelin/SpatialRBInertiaf.h>
#include <Ravelin/SpatialABInertiaf.h>

namespace Ravelin {

#include "fdefs.h"
#include "Pose3.h"
#include "undefs.h"

} // end namespace

#endif

