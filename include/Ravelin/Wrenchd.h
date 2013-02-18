/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _WRENCHD_H
#define _WRENCHD_H

#include <Moby/SVector6d.h>

namespace Moby {

#define REAL double 
#define WRENCH Wrenchd
#define SVECTOR6 SVector6d
#define VECTOR3 Vector3d

#include "Wrench.h"

#undef REAL
#undef WRENCH
#undef SVECTOR6
#undef VECTOR3

} // end namespace

#endif

