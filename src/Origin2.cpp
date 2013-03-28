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

/// Constructs this vector with the given values
ORIGIN2::ORIGIN2(const POINT2& p)
{
  const unsigned X = 0, Y = 1;
  _data[X] = p.x();
  _data[Y] = p.y();
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

/// Assigns the values from a point to an origin
ORIGIN2& ORIGIN2::operator=(const POINT2& p)
{
  x() = p.x(); 
  y() = p.y(); 
  return *this;
}

/// Adds a point and an origin to yield a point
POINT2 ORIGIN2::operator+(const POINT2& p) const
{
  POINT2 result;
  result.x() = x() + p.x(); 
  result.y() = y() + p.y();
  result.pose = p.pose; 
  return result;
}

/// Subtract a point from this origin to yield a point
POINT2 ORIGIN2::operator-(const POINT2& p) const
{
  POINT2 result;
  result.x() = x() - p.x(); 
  result.y() = y() - p.y();
  result.pose = p.pose; 
  return result;
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


