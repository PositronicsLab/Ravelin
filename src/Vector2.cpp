/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Constructs this vector with the given values
VECTOR2::VECTOR2(REAL x, REAL y, boost::shared_ptr<POSE2> pose)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
  this->pose = boost::const_pointer_cast<const POSE2>(pose);
}

/// Constructs this vector with the given values
VECTOR2::VECTOR2(REAL x, REAL y, boost::shared_ptr<const POSE2> pose)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
  this->pose = pose;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
VECTOR2::VECTOR2(const REAL* array, boost::shared_ptr<const POSE2> pose)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
  this->pose = pose;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
VECTOR2::VECTOR2(const REAL* array, boost::shared_ptr<POSE2> pose)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
  this->pose = boost::const_pointer_cast<const POSE2>(pose);
}

/// Constructs a vector from an origin object
VECTOR2& VECTOR2::operator=(const ORIGIN2& o)
{
  x() = o.x();
  y() = o.y();
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
  if (!((m == 2 && n == 1) || (m == 1 && n == 2)))
    throw std::runtime_error("Attempt to resize fixed-length vector!");
  #endif

  return *this;
}
 
COLUMN_ITERATOR VECTOR2::column_iterator_begin()
{
  COLUMN_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_COLUMN_ITERATOR VECTOR2::column_iterator_begin() const
{
  CONST_COLUMN_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

COLUMN_ITERATOR VECTOR2::column_iterator_end()
{
  COLUMN_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_COLUMN_ITERATOR VECTOR2::column_iterator_end() const
{
  CONST_COLUMN_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

ROW_ITERATOR VECTOR2::row_iterator_begin()
{
  ROW_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_ROW_ITERATOR VECTOR2::row_iterator_begin() const
{
  CONST_ROW_ITERATOR i;
  i._data_start = i._current_data = _data;
  i._count = 0;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

ROW_ITERATOR VECTOR2::row_iterator_end()
{
  ROW_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
}

CONST_ROW_ITERATOR VECTOR2::row_iterator_end() const
{
  CONST_ROW_ITERATOR i;
  i._data_start = _data;
  i._current_data = _data+2;
  i._count = 2;
  i._sz = 2;
  i._rows = 2;
  i._columns = 1;
  i._ld = 2;
  return i;
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

/// Adds an origin to this
VECTOR2& VECTOR2::operator+=(const ORIGIN2& o) 
{
  // do the addition 
  _data[0] += o.x();
  _data[1] += o.y();

  return *this;
}

/// Subtracts an origin from this
VECTOR2& VECTOR2::operator-=(const ORIGIN2& o) 
{
  // do the subtraction
  _data[0] -= o.x();
  _data[1] -= o.y();

  return *this;
}

