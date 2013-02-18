/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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
#include <Ravelin/VectorNf.h>

namespace Ravelin {

#define INTEGER int
#define LINALG LinAlgf
#define REAL float
#define MATRIXN MatrixNf
#define VECTORN VectorNf
#define MATRIX2 Matrix2f
#define CONST_ITERATOR fIterator_const
#define ITERATOR fIterator

#include "LinAlg.h"

#undef INTEGER
#undef LINALG
#undef REAL
#undef MATRIXN
#undef VECTORN
#undef MATRIX2
#undef CONST_ITERATOR
#undef ITERATOR

} // end namespace

#endif

