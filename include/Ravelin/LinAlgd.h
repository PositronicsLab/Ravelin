/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RAVELIN_LINALGD_H
#define _RAVELIN_LINALGD_H

#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/minmax.hpp>
#include <Ravelin/SingularException.h>
#include <Ravelin/NumericalException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/Matrix2d.h>
#include <Ravelin/FastThreadable.h>
#include <Ravelin/SparseMatrixNd.h>
#include <Ravelin/VectorNd.h>

namespace Ravelin {

#define INTEGER int
#include "ddefs.h"
#include "LinAlg.h"
#undef INTEGER
#include "undefs.h"

} // end namespace

#endif
