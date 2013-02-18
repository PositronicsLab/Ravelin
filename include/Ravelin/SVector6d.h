/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SVECTOR6D_H
#define _SVECTOR6D_H

#include <Moby/Vector3d.h>

namespace Moby {

#define REAL double
#define SVECTOR6 SVector6d
#define VECTOR3 Vector3d

#include "SVector6.h"

#undef REAL
#undef SVECTOR6
#undef VECTOR3

} // end namespace

#endif

