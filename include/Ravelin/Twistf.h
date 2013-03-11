/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _TWISTF_H
#define _TWISTF_H

#include <Ravelin/SVector6f.h>

namespace Ravelin {

#define REAL float
#define TWIST Twistf
#define WRENCH Wrenchf
#define SVECTOR6 SVector6f
#define VECTOR3 Vector3f
#define POSE Posef 

#include "Twist.h"

#undef REAL
#undef TWIST
#undef WRENCH
#undef SVECTOR6
#undef VECTOR3
#undef POSE

} // end namespace

#endif

