/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 *
 * This file contains code specific to SharedVectorNf/SharedVectorNd.
 ****************************************************************************/

/// Default constructor - constructs an empty vector
SHAREDVECTORN::SHAREDVECTORN()
{
  _len = 0;
  _inc = 1;
  _start = 0;
}

/// Constructs a shared vector from another shared vector 
SHAREDVECTORN::SHAREDVECTORN(const SHAREDVECTORN& v)
{
  _len = v._len;
  _inc = v._inc;
  _data = v._data;
  _start = v._start;
}

/// Gets a shared subvector of this subvector
SHAREDVECTORN SHAREDVECTORN::get_sub_vec(unsigned start, unsigned end)
{
  #ifndef NEXCEPT
  if (start > end || end > _len)
    throw InvalidIndexException();
  #endif

  SHAREDVECTORN x;
  x._data = _data;
  x._start = _start + start;
  x._len = end - start;
  x._inc = 1;
  return x;
}

/// Does nothing
SHAREDVECTORN& SHAREDVECTORN::resize(unsigned N, bool preserve)
{
  // if the array already is the proper size, exit
  #ifndef NEXCEPT
  if (_len != N)
    throw std::runtime_error("Attempt to resize shared vector!");
  #endif

  return *this;
}

/// Copies another vector
SHAREDVECTORN& SHAREDVECTORN::operator=(const VECTOR3& source)
{
  // check size
  #ifndef NEXCEPT
  if (_len != source.size())
    throw MissizeException();
  #endif

  // get data
  REAL* x = data();

  // don't even worry about BLAS for copying
  *x = source.x(); x += inc();
  *x = source.y(); x += inc();
  *x = source.z();

  return *this;
}

/// Copies another vector 
SHAREDVECTORN& SHAREDVECTORN::operator=(const VECTORN& source)
{  
  // check size 
  if (_len != source.size())
    throw MissizeException();

  // use the BLAS routine for copying
  if (_len > 0)
    CBLAS::copy(source.size(),source.data(),source.inc(),data(),inc());

  return *this;
}

/// Copies another vector 
SHAREDVECTORN& SHAREDVECTORN::operator=(const CONST_SHAREDVECTORN& source)
{  
  // check size 
  if (_len != source.size())
    throw MissizeException();

  // use the BLAS routine for copying
  if (_len > 0)
    CBLAS::copy(source.size(),source.data(),source.inc(),data(),inc());

  return *this;
}

/// Copies another vector 
SHAREDVECTORN& SHAREDVECTORN::operator=(const SHAREDVECTORN& source)
{  
  // check size 
  if (_len != source.size())
    throw MissizeException();

  // use the BLAS routine for copying
  if (_len > 0)
    CBLAS::copy(source.size(),source.data(),source.inc(),data(),inc());

  return *this;
}

/// Default constructor - constructs an empty vector
CONST_SHAREDVECTORN::CONST_SHAREDVECTORN()
{
  _len = 0;
  _inc = 1;
  _start = 0;
}

/// Constructs a shared vector from another shared vector 
CONST_SHAREDVECTORN::CONST_SHAREDVECTORN(const CONST_SHAREDVECTORN& v)
{
  _len = v._len;
  _inc = v._inc;
  _data = v._data;
  _start = v._start;
}

/// Constructs a shared vector from another shared vector 
CONST_SHAREDVECTORN::CONST_SHAREDVECTORN(const SHAREDVECTORN& v)
{
  _len = v._len;
  _inc = v._inc;
  _data = v._data;
  _start = v._start;
}

/// Gets a shared subvector of this subvector
CONST_SHAREDVECTORN CONST_SHAREDVECTORN::get_sub_vec(unsigned start, unsigned end)
{
  #ifndef NEXCEPT
  if (start > end || end > _len)
    throw InvalidIndexException();
  #endif

  CONST_SHAREDVECTORN x;
  x._data = _data;
  x._start = _start + start;
  x._len = end - start;
  x._inc = 1;
  return x;
}

/// Does nothing
CONST_SHAREDVECTORN& CONST_SHAREDVECTORN::resize(unsigned N, bool preserve)
{
  // if the array already is the proper size, exit
  #ifndef NEXCEPT
  if (_len != N)
    throw std::runtime_error("Attempt to resize shared vector!");
  #endif

  return *this;
}


#define XVECTORN SHAREDVECTORN
#include "XVectorN.cpp"
#undef XVECTORN

/// Writes a XVECTORN to the specified stream
std::ostream& Ravelin::operator<<(std::ostream& out, const CONST_SHAREDVECTORN& v)
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

