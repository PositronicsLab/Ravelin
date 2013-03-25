/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
