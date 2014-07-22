/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RAVELIN_LINALGF_H
#define _RAVELIN_LINALGF_H

#include <boost/tuple/tuple.hpp>
#include <boost/algorithm/minmax.hpp>
#include <vector>
#include <Ravelin/SingularException.h>
#include <Ravelin/NumericalException.h>
#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MatrixNf.h>
#include <Ravelin/Matrix2f.h>
#include <Ravelin/FastThreadable.h>
#include <Ravelin/SparseMatrixNf.h>
#include <Ravelin/VectorNf.h>

namespace Ravelin {

#define INTEGER int
#include "fdefs.h"
#include "LinAlg.h"
#include "undefs.h"
#undef INTEGER

} // end namespace

#endif

