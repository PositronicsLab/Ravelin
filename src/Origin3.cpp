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

/// Assigns this origin using the 3D vector 
ORIGIN3& ORIGIN3::operator=(const VECTOR3& v)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = v.x();
  _data[Y] = v.y();
  _data[Z] = v.z();
  return *this;
}

/// Assigns the values from one origin to this origin
ORIGIN3& ORIGIN3::operator=(const ORIGIN3& o)
{
  x() = o.x(); 
  y() = o.y(); 
  z() = o.z(); 
  return *this;
}

/// Does nothing
ORIGIN3& ORIGIN3::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 3 || n != 1)
    throw MissizeException();
  #endif
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

/// Adds a vector to an origin
VECTOR3 ORIGIN3::operator+(const VECTOR3& v) const
{
  VECTOR3 result(v.pose);
  result.x() = x() + v.x(); 
  result.y() = y() + v.y();
  result.z() = z() + v.z();
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

/// Subtracts a vector from an origin
VECTOR3 ORIGIN3::operator-(const VECTOR3& v) const
{
  VECTOR3 result(v.pose);
  result.x() = x() - v.x(); 
  result.y() = y() - v.y();
  result.z() = z() - v.z();
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

