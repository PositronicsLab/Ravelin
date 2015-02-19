/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
ORIGIN2& ORIGIN2::resize(unsigned m, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 2)
    throw MissizeException();
  #endif
  return *this;
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
ORIGIN2& ORIGIN2::operator=(const VECTOR2& v)
{
  const unsigned X = 0, Y = 1;
  _data[X] = v.x();
  _data[Y] = v.y();
  return *this;
}

/// Adds a vector and an origin to yield a vector 
VECTOR2 ORIGIN2::operator+(const VECTOR2& v) const
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


