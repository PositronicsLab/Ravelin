/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
VECTOR2::VECTOR2(REAL x, REAL y)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
}

/// Constructs this vector with the given values
VECTOR2::VECTOR2(REAL x, REAL y, boost::shared_ptr<const POSE2> rpose)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
  pose = rpose;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
VECTOR2::VECTOR2(const REAL* array)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
}

/// Constructs a vector from a point
VECTOR2& VECTOR2::operator=(const POINT2& p)
{
  pose = p.pose;
  x() = p.x();
  y() = p.y();
  return *this;
}

REAL& VECTOR2::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

const REAL& VECTOR2::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

REAL* VECTOR2::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

const REAL* VECTOR2::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

VECTOR2& VECTOR2::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 2 || n != 1)
    throw std::runtime_error("Attempt to resize fixed-length vector!");
  #endif

  return *this;
}
 
ITERATOR VECTOR2::begin()
{
  ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_ITERATOR VECTOR2::begin() const
{
  CONST_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

ITERATOR VECTOR2::end()
{
  ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_ITERATOR VECTOR2::end() const
{
  CONST_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

/// Adds two vectors together
VECTOR2 VECTOR2::operator+(const VECTOR2& v) const
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif
  VECTOR2 result = *this;
  result += v;
  return result;
}

/// Adds a vector to this
VECTOR2& VECTOR2::operator+=(const VECTOR2& v) 
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  // do the addition 
  _data[0] += v._data[0];
  _data[1] += v._data[1];

  return *this;
}

/// Subtracts a vector from this
VECTOR2 VECTOR2::operator-(const VECTOR2& v) const
{
  VECTOR2 result = *this;
  result -= v;
  return result;
}

/// Subtracts a vector from this
VECTOR2& VECTOR2::operator-=(const VECTOR2& v) 
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  // do the subtraction
  _data[0] -= v._data[0];
  _data[1] -= v._data[1];

  return *this;
}

