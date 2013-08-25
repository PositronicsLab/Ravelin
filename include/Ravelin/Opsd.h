#ifndef _RAVELIN_OPERATORS_D_H_
#define _RAVELIN_OPERATORS_D_H_

#include <Ravelin/cblas.h>
#include <Ravelin/Constants.h>
#include <Ravelin/DataMismatchException.h>
#include <Ravelin/MissizeException.h>

namespace Ravelin {

#include "ddefs.h"
#include "Ops.h"
#include "undefs.h"

/// Determines whether two doubles are equal
inline bool Opsd::rel_equal(double x, double y)
{
  return (std::fabs(x-y) <= EPS_DOUBLE * std::max(std::fabs(x), std::max(std::fabs(y), 1.0)));
}

/// Determines whether two doubles are equal
inline bool Opsd::rel_equal(double x, double y, double tol)
{
  return (std::fabs(x-y) <= tol * std::max(std::fabs(x), std::max(std::fabs(y), 1.0)));
}


} // end namespace

#endif

