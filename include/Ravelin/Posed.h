/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _POSED_H
#define _POSED_H

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/FrameException.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Wrenchd.h>
#include <Ravelin/Twistd.h>
#include <Ravelin/SpatialRBInertiad.h>
#include <Ravelin/SpatialABInertiad.h>

namespace Ravelin {

#include "ddefs.h"
#include "Pose.h"
#include "undefs.h"

} // end namespace

#endif

