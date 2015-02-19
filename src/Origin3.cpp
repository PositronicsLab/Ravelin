/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
ORIGIN3& ORIGIN3::resize(unsigned m, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 3)
    throw MissizeException();
  #endif
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

COLUMN_ITERATOR ORIGIN3::column_iterator_begin()
{
  COLUMN_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_COLUMN_ITERATOR ORIGIN3::column_iterator_begin() const
{
  CONST_COLUMN_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

COLUMN_ITERATOR ORIGIN3::column_iterator_end()
{
  COLUMN_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_COLUMN_ITERATOR ORIGIN3::column_iterator_end() const
{
  CONST_COLUMN_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

ROW_ITERATOR ORIGIN3::row_iterator_begin()
{
  ROW_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_ROW_ITERATOR ORIGIN3::row_iterator_begin() const
{
  CONST_ROW_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

ROW_ITERATOR ORIGIN3::row_iterator_end()
{
  ROW_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

CONST_ROW_ITERATOR ORIGIN3::row_iterator_end() const
{
  CONST_ROW_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+3;
  i._count = 3;
  i._sz = 3;
  i._rows = 3;
  i._columns = 1;
  i._ld = 3;
  return i;
}

/// Computes the cross-product of two vectors
ORIGIN3 ORIGIN3::cross(const ORIGIN3& v1, const ORIGIN3& v2)
{
  ORIGIN3 w;
  w[0] = v1[1] * v2[2] - v1[2] * v2[1];
  w[1] = v1[2] * v2[0] - v1[0] * v2[2];
  w[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return w;
}
