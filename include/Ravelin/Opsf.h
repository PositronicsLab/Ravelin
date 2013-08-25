#ifndef _RAVELIN_OPERATORS_F_H_
#define _RAVELIN_OPERATORS_F_H_

#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/MissizeException.h>

namespace Ravelin {

#include "fdefs.h"
#include "Ops.h"
#include "undefs.h"

inline bool Opsf::rel_equal(float x, float y)
{
  return (std::fabs(x-y) <= EPS_FLOAT * std::max(std::fabs(x), std::max(std::fabs(y), 1.0f)));
}

/// Determines whether two floats are equal
inline bool Opsf::rel_equal(float x, float y, float tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), 
1.0f)));
}

} // end namespace

#endif

