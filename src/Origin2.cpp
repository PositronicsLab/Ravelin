/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
ORIGIN2::ORIGIN2(REAL x, REAL y)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
ORIGIN2::ORIGIN2(const REAL* array)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
}

/// Does nothing
ORIGIN2& ORIGIN2::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 2 || n != 1)
    throw MissizeException();
  #endif
  return *this;
}

/// Constructs this vector with the given values
inline ORIGIN2& ORIGIN2::operator=(const POINT2& p)
{
  const unsigned X = 0, Y = 1;
  _data[X] = p.x();
  _data[Y] = p.y();
  return *this;
}

/// Constructs this vector with the given values
inline ORIGIN2& ORIGIN2::operator=(const VECTOR2& v)
{
  const unsigned X = 0, Y = 1;
  _data[X] = v.x();
  _data[Y] = v.y();
  return *this;
}

/// Adds a point and an origin to yield a vector 
inline VECTOR2 ORIGIN2::operator+(const POINT2& p) const
{
  return p + *this;
}

/// Subtract a point from this origin to yield a vector 
inline VECTOR2 ORIGIN2::operator-(const POINT2& p) const
{
  return -p + *this;
}

/// Adds a vector and an origin to yield a vector 
inline VECTOR2 ORIGIN2::operator+(const VECTOR2& v) const
{
  return v + *this; 
}

/// Subtract a vector from this origin to yield a vector 
VECTOR2 ORIGIN2::operator-(const VECTOR2& v) const
{
  return -v + *this;
}

REAL& ORIGIN2::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

const REAL& ORIGIN2::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

REAL* ORIGIN2::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

const REAL* ORIGIN2::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}


