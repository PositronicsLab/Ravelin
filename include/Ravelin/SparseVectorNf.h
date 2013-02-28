/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_VECTOR_NF_H_
#define _SPARSE_VECTOR_NF_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Ravelin/VectorNf.h>

namespace Ravelin {

#define REAL float
#define SPARSEVECTORN SparseVectorNf
#define VECTORN VectorNf

#include "SparseVectorN.h"

#undef REAL
#undef SPARSEVECTORN
#undef VECTORN

} // end namespace

#endif

