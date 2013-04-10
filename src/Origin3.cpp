/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this origin with the given values
ORIGIN3::ORIGIN3(REAL x, REAL y, REAL z)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
}

/// Constructs this origin with the given values
ORIGIN3::ORIGIN3(const POINT3& p)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = p.x();
  _data[Y] = p.y();
  _data[Z] = p.z();
}

/// Constructs this origin from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
ORIGIN3::ORIGIN3(const REAL* array)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
}

/// Assigns the values from one origin to this origin
ORIGIN3& ORIGIN3::operator=(const ORIGIN3& o)
{
  x() = o.x(); 
  y() = o.y(); 
  z() = o.z(); 
  return *this;
}

/// Assigns the values from a point to an origin
ORIGIN3& ORIGIN3::operator=(const POINT3& p)
{
  x() = p.x(); 
  y() = p.y(); 
  z() = p.z(); 
  return *this;
}

/// Adds two origins
ORIGIN3 ORIGIN3::operator+(const ORIGIN3& o) const
{
  ORIGIN3 result;
  result.x() = x() + o.x(); 
  result.y() = y() + o.y();
  result.z() = z() + o.z();
  return result;
}

/// Adds two origins
ORIGIN3& ORIGIN3::operator+=(const ORIGIN3& o)
{
  x() += o.x(); 
  y() += o.y();
  z() += o.z();
  return *this;
}

/// Adds a point and an origin to yield a point
POINT3 ORIGIN3::operator+(const POINT3& p) const
{
  POINT3 result;
  result.x() = x() + p.x(); 
  result.y() = y() + p.y();
  result.z() = z() + p.z();
  result.pose = p.pose; 
  return result;
}

/// Subtracts one origin from another
ORIGIN3 ORIGIN3::operator-(const ORIGIN3& o) const
{
  ORIGIN3 result;
  result.x() = x() - o.x(); 
  result.y() = y() - o.y();
  result.z() = z() - o.z();
  return result;
}

/// Subtracts an origin from this 
ORIGIN3& ORIGIN3::operator-=(const ORIGIN3& o)
{
  x() -= o.x(); 
  y() -= o.y();
  z() -= o.z();
  return *this;
}

/// Subtract a point from this origin to yield a point
POINT3 ORIGIN3::operator-(const POINT3& p) const
{
  POINT3 result;
  result.x() = x() - p.x(); 
  result.y() = y() - p.y();
  result.z() = z() - p.z();
  result.pose = p.pose; 
  return result;
}

REAL& ORIGIN3::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

const REAL& ORIGIN3::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

REAL* ORIGIN3::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

const REAL* ORIGIN3::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

ITERATOR ORIGIN3::begin()
{
  ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_ITERATOR ORIGIN3::begin() const
{
  CONST_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

ITERATOR ORIGIN3::end()
{
  ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_ITERATOR ORIGIN3::end() const
{
  CONST_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

