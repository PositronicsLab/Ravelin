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
#include <Ravelin/VectorNd.h>

namespace Ravelin {

#define INTEGER int
#define LINALG LinAlgd
#define REAL double
#define MATRIXN MatrixNd
#define VECTORN VectorNd
#define MATRIX2 Matrix2d
#define ITERATOR dIterator
#define CONST_ITERATOR dIterator_const
#include "LinAlg.h"
#undef INTEGER
#undef LINALG
#undef REAL
#undef MATRIXN
#undef VECTORN
#undef MATRIX2
#undef ITERATOR
#undef CONST_ITERATOR

} // end namespace

#endif
