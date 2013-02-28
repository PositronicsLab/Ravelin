/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPARSE_VECTOR_ND_H_
#define _SPARSE_VECTOR_ND_H_

#include <map>
#include <boost/shared_ptr.hpp>
#include <Ravelin/VectorNd.h>

namespace Ravelin {

#define REAL double
#define SPARSEVECTORN SparseVectorNd
#define VECTORN VectorNd

#include "SparseVectorN.h"

#undef REAL
#undef SPARSEVECTORN
#undef VECTORN

} // end namespace

#endif

