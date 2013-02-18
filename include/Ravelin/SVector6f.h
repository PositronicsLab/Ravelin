/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SVECTOR6F_H
#define _SVECTOR6F_H

#include <Moby/Vector3f.h>

namespace Moby {

#define REAL float
#define SVECTOR6 SVector6f
#define POSE Posef
#define VECTOR3 Vector3f

#include "SVector6.h"

#undef REAL
#undef SVECTOR6
#undef POSE
#undef VECTOR3

} // end namespace

#endif

