/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 *
 * This file contains code for VectorNf/VectorNd and 
 * SharedVectorNf/SharedVectorNd.
 ****************************************************************************/

#include <iomanip>

// Reads a VECTORN from the specified stream
std::istream& Ravelin::operator>>(std::istream& in, XVECTORN& v)
{
  unsigned n;
  in >> n;
  v.resize(n);

  for (unsigned i=0; i< n; i++)
    in >> v[i];

  return in;
}

/// Writes a XVECTORN to the specified stream
std::ostream& Ravelin::operator<<(std::ostream& out, const XVECTORN& v)
{
  const unsigned OUTPUT_PRECISION = 8;

  if (v.size() == 0)
  {
    out << "(empty) ";
    return out;
  }
  out << "[";
  for (unsigned i=0; i< v.size()-1; i++)
    out << std::setprecision(OUTPUT_PRECISION) << v[i] << ", ";
  out << std::setprecision(OUTPUT_PRECISION) << v[v.size()-1] << "] ";
  return out;
}


