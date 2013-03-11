/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _WRENCHD_H
#define _WRENCHD_H

#include <Ravelin/SVector6d.h>

namespace Ravelin {

#define REAL double 
#define TWIST Twistd
#define WRENCH Wrenchd
#define SVECTOR6 SVector6d
#define VECTOR3 Vector3d
#define POSE Posed 

#include "Wrench.h"

#undef REAL
#undef TWIST
#undef WRENCH
#undef SVECTOR6
#undef VECTOR3
#undef POSE

} // end namespace

#endif

