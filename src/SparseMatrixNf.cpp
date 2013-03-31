/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU  General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <Ravelin/FastThreadable.h>
#include <Ravelin/Constants.h>
#include <Ravelin/MissizeException.h>
#include <Ravelin/InvalidIndexException.h>
#include <Ravelin/SparseMatrixNf.h>
#include <Ravelin/MatrixNf.h>

using std::pair;
using boost::shared_array;
using std::map;
using std::make_pair;
using std::vector;
using namespace boost::lambda;
using namespace Ravelin;


#include <Ravelin/fdefs.h>
#include "SparseMatrixN.cpp"
#include <Ravelin/undefs.h>


