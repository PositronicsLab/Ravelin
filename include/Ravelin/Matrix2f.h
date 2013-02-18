/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MATRIX2F_H
#define _MATRIX2F_H

#include <vector>
#include <algorithm>
#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/Vector2f.h>

namespace Ravelin {

#define REAL float 
#define MATRIX2 Matrix2f
#define VECTOR2 Vector2f
#define ITERATOR fIterator 
#define CONST_ITERATOR fIterator_const 

#include "Matrix2.h"

#undef REAL
#undef MATRIX2
#undef VECTOR2
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif

